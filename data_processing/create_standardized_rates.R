#########################################################################
## Using the long-form MA census tract population and mortality data   ##
## create CT-level crude and age- and sex-standardized mortality rates ##
#########################################################################

library(maptools)
library(spdep)
library(MASS)
library(msm)
library(tigris)
library(CARBayes)
library(ggplot2)
library(sf)
library(tmap)
library(gtools)
library(reshape2)
library(rstudioapi)

## set working directory to wherever this file is located ##
current_path <- getActiveDocumentContext()$path 
setwd(dirname(current_path ))

################################################################
## 1. read in the long-form data with merged mortality counts ##
################################################################

load('merged_pm_denom_cov_data.RData')
adat<-subset(adat,race %in% c('total','white','black','hispanic'))

## get age- and sex-stratified mortality counts for standardization ##
mort_ref_1yr<-aggregate(ndeaths_1yr~agecat+sex,data=subset(adat,race=='total'),FUN=sum,na.rm=T)
names(mort_ref_1yr)[3]<-'deaths_agexsex'

mort_ref_5yr<-aggregate(ndeaths_5yr~agecat+sex,data=subset(adat,race=='total'),FUN=sum,na.rm=T)
names(mort_ref_5yr)[3]<-'deaths_agexsex'

#############################################################
## 2. function to create both crude and standardized rates ##
#############################################################

agestand<-function(dat,popnm,incnm,mort_ref,d,nyrs){
  datnew<-dat[,c('GEOID','race','agecat','sex',popnm,incnm)]
  datnew[[popnm]]<-datnew[[popnm]]*nyrs
  dattot<-datnew[which(datnew$race=='total'),]
  
  ## get pop counts in each age/sex group for standardization ##
  pop_ref<-aggregate(dattot[[popnm]],by=list(dattot$agecat,dattot$sex),FUN=sum,na.rm=T)
  names(pop_ref)<-c('agecat','sex','pop_agexsex')
  
  ## get reference mortality rates ##
  mort_ref<-merge(mort_ref,pop_ref,by=c('agecat','sex'))
  mort_ref$mr<-mort_ref$deaths_agexsex/mort_ref$pop_agexsex
  mort_ref$pop_prop<-mort_ref$pop_agexsex/sum(mort_ref$pop_agexsex,na.rm=T)
  
  ## make zero populations missing (otherwise we get infinite rates) ##
  datnew[[popnm]][which(datnew[[popnm]]==0)]<-NA

  ## data to compute crude rates by CT and race ##
  crudedat<-aggregate(list(datnew[[incnm]],datnew[[popnm]]),by=list(datnew$GEOID,datnew$race),FUN=sum,na.rm=T)
  names(crudedat)<-c('GEOID','race','ndeaths','pop')
  ## compute crude rates ##
  crudedat$pmrate_cr<-crudedat$ndeaths/crudedat$pop
  
  ## data to do direct age standardization ##
  dstdat<-merge(datnew,mort_ref,by=c('agecat','sex'))
  dstdat$pmrate_dst<-(dstdat[[incnm]]/dstdat[[popnm]])*dstdat[['pop_prop']]
  dstdat<-aggregate(pmrate_dst~GEOID+race,data=dstdat,FUN=sum,na.rm=T)
  
  ## data to do indirect age standardization ##
  ## also add both observed counts and expected (separately) to the data ##
  istdat<-merge(datnew,mort_ref,by=c('agecat','sex'))
  istdat$istexp<-(istdat[[popnm]]*istdat[['mr']])
  istdat<-aggregate(list(istdat$ndeaths,istdat$istexp),by=list(istdat$GEOID,istdat$race),FUN=sum,na.rm=T)
  names(istdat)<-c('GEOID','race','ndeaths','istexp')
  istdat$pmrate_ist<-istdat$ndeaths/istdat$istexp

  ## merge crude and standardized rates ##
  m1<-merge(crudedat,dstdat,by=c('GEOID','race'),all.x=T)
  m2<-merge(m1,istdat[,-3],by=c('GEOID','race'),all.x=T)
  
  ## make NA rates in places with 0 population ##
  m2[which(m2$pop==0),5:ncol(m2)]<-NA
  
  names(m2)[3:ncol(m2)]<-c(paste0('ndeaths_',nyrs,'yr'),paste0(d,'_',nyrs,'yr_',names(m2)[4:ncol(m2)]))

  return(m2)
}

##############################################################################################################
## 3. use the function to create 1-year and 5-year mortality rates with each set of population denominators ##
##############################################################################################################

pmrate_ce_1yr<-agestand(dat=adat,popnm='ce_pop',incnm='ndeaths_1yr',mort_ref=mort_ref_1yr,d='ce',nyrs=1)
pmrate_acs_1yr<-agestand(dat=adat,popnm='acs_pop',incnm='ndeaths_1yr',mort_ref=mort_ref_1yr,d='acs',nyrs=1)
pmrate_dp_1yr<-agestand(dat=adat,popnm='dp_pop',incnm='ndeaths_1yr',mort_ref=mort_ref_1yr,d='dp',nyrs=1)

pmrate_ce_5yr<-agestand(dat=adat,popnm='ce_pop',incnm='ndeaths_5yr',mort_ref=mort_ref_5yr,d='ce',nyrs=5)
pmrate_acs_5yr<-agestand(dat=adat,popnm='acs_pop',incnm='ndeaths_5yr',mort_ref=mort_ref_5yr,d='acs',nyrs=5)
pmrate_dp_5yr<-agestand(dat=adat,popnm='dp_pop',incnm='ndeaths_5yr',mort_ref=mort_ref_5yr,d='dp',nyrs=5)

pmrate_list<-list(pmrate_ce_1yr,pmrate_ce_5yr,pmrate_acs_1yr[,-3],pmrate_acs_5yr[,-3],pmrate_dp_1yr[,-3],pmrate_dp_5yr[,-3])

pmrate<-Reduce(function(...) merge(..., by=c('GEOID','race'), all=T), pmrate_list)

##################################################
## 4. merge back in the covariate data and save ##
##################################################

covdat<-adat[!duplicated(adat[,c('GEOID','race')]),c('GEOID','race','acs_pov_pct','acs_ice_racewb',
                                                     'acs_ice_inc','acs_ice_raceinc','ce_ice_racewb','dp_ice_racewb')]

adat<-merge(pmrate,covdat,by=c('GEOID','race'))

## sort by geoid ##
adat<-adat[order(adat$GEOID),]

save(adat,file='pm_rates.RData')
