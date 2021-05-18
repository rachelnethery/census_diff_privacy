## pseudo-simulations using 2010 ACS data ##
## Generate data with "true" census pop and use ACS pop estimate in models ##
## this code uses raw (non-standardized) population denominators ##
## created by Rachel Nethery ##
## date: 3/1/20 ##

## read command line arguments ##
args<-commandArgs(TRUE)
for (i in 1:length(args)) { eval (parse (text = args[[i]] )) }


setwd(wd)

set.seed(2)

library(maptools)
library(spdep)
library(MASS)
library(msm)
library(tigris)
library(CARBayes)

#reps<-100
#nburn<-10000
#nsamp<-30000

####################################################
## 1. population size and covariate info from ACS ##
####################################################

## load the pre-prepped MA CT-level ACS/census/dp data ##
load('merged_denom_cov_data.RData')

## subset to only black and white ##
adat<-subset(adat,race %in% c('white','black'))

## aggregate across age/sex groups ##
acs_pop <- aggregate(cbind(acs_pop,ce_pop)~GEOID+race,data=adat, 
                     FUN = sum,na.rm=T)

acs_cov <- aggregate(cbind(acs_pov_pct,acs_ice_raceinc)~GEOID+race,data=adat,
                     FUN = mean)

## merge these aggregated datasets ##
adat<-merge(acs_pop,acs_cov,by=c('GEOID','race'))

## remove CTs with 0 pop ##
adat<-subset(adat,acs_pop>0 & ce_pop>0)

adat<-adat[order(adat$GEOID),]

## create the matrices simulate ##
adat$bw_ind<-0
adat$bw_ind[which(adat$race=='black')]<-1

X<-cbind(1,adat$bw_ind,adat$acs_pov_pct)
colnames(X)<-c('intercept','bw','pov')
P<-matrix(adat$acs_pop)
Ptrue<-matrix(adat$ce_pop)

N<-nrow(X)
Nct<-length(unique(adat$GEOID))

###########################################
## 2. adjacency matrix for census tracts ##
###########################################

## extract shapefile ##
ma_shp<-tracts(state = 'MA',year=2010)

## get adjacency matrix ##
foo = poly2nb(ma_shp, queen=TRUE, row.names=ma_shp@data$GEOID10)
W=nb2mat(foo,style='B')

##########################################################################
## 3. align ordering of adjacency matrix and covariates/population data ##
##########################################################################

## census tracts in W with matches in the data ##
tractrm<-which(rownames(W) %in% as.character(adat$GEOID))

W<-W[tractrm,tractrm]

W<-W[order(as.numeric(rownames(W))),order(as.numeric(rownames(W)))]

identical(unique(as.character(adat$GEOID)),rownames(W))

########################################
## 4. create simulated data structure ##
########################################

## spatial information ##
Dw<-diag(rowSums(W))
#rho_bds<-eigen(sqrt(solve(Dw))%*%W%*%sqrt(solve(Dw)))$values
#print(1/min(rho_bds))
#print(1/max(rho_bds))
rho<-.2

## lambda=scalar background rate ##
lambda<-.03

## beta coefficients for covariates ##
beta<-c(0,0.4,0.01)

###################################
## 5. set up storage for results ##
###################################

cars<-matrix(NA,nrow=reps,ncol=length(beta)*3+10)

########################
## 6. run simulations ##
########################

set.seed(simnum)

for (xx in 1:reps){
  
  ## U (spatial component) ##
  U<-matrix(mvrnorm(n=1,mu=rep(0,Nct),Sigma=solve(Dw-(rho*W))))
  U<-U[as.numeric(as.factor(adat$GEOID)),]
  
  ## V (unsturctured variance component) ##
  V<-matrix(rnorm(n=N,mean=0,sd=.5))
  
  ## generate outcome ##
  
  theta<-X%*%beta+U+V
  
  Y<-matrix(rpois(n=N,lambda=lambda*Ptrue*exp(theta)))
  
  ####################################
  ## FIT STANDARD SPATIAL CAR MODEL ##
  ####################################
  adat$Y<-Y
  adat$offs<-log(lambda*P)
  results_cars<-S.CARmultilevel(Y~bw_ind+acs_pov_pct+offset(offs),family="poisson",data=adat,
                                ind.area=as.numeric(as.factor(adat$GEOID)),
                                ind.re=as.factor(1:N),
                                W=W,burnin=nburn,n.sample=nsamp)
  
  cars[xx,]<-c(apply(results_cars$samples$beta,2,mean),
               apply(results_cars$samples$beta,2,quantile,.025),
               apply(results_cars$samples$beta,2,quantile,.975),
               apply(results_cars$samples$tau2,2,mean),
               apply(results_cars$samples$tau2,2,quantile,.025),
               apply(results_cars$samples$tau2,2,quantile,.975),
               apply(results_cars$samples$sigma2,2,mean),
               apply(results_cars$samples$sigma2,2,quantile,.025),
               apply(results_cars$samples$sigma2,2,quantile,.975),
               apply(results_cars$samples$rho,2,mean),
               apply(results_cars$samples$rho,2,quantile,.025),
               apply(results_cars$samples$rho,2,quantile,.975))
  
  
}

cars<-as.data.frame(cars)

names(cars)<-c(paste0('beta',0:(length(beta)-1)),paste0('beta',0:(length(beta)-1),'_lo'),
               paste0('beta',0:(length(beta)-1),'_hi'),'tau2','tau2_lo','tau2_hi',
               'sigma2','sigma2_lo','sigma2_hi','rho','rho_lo','rho_hi','mse')

save(cars,file=paste0('results/simnum',simnum,'.RData'))
