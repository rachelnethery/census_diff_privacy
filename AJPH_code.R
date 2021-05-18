#################################################
## Data processing and analyses for AJPH paper ##
#################################################

library(tidyverse)
library(tidycensus)
library(tigris)
library(sp)
library(stringr)
library(haven)
library(reshape2)
library(rstudioapi)
library(dplyr)
library(RColorBrewer)
library(purrr)

## set working directory to wherever this file is located ##
current_path <- getActiveDocumentContext()$path 
setwd(dirname(current_path ))

#################
## 1. ACS DATA ##
#################

## variable names ##
#v12 <- load_variables(2012, "acs5")

##extract using tidycensus ##
ma_acs<- get_acs(geography = "tract",
                 variables = c(## total pop <65
                   paste0('B01001_',str_pad(as.character(c(3:19,27:43)),width=3,side='left',pad='0')),
                   ## non-hispanic white <65
                   paste0('B01001H_',str_pad(as.character(c(3:13,18:28)),width=3,side='left',pad='0')),
                   ## black <65
                   paste0('B01001B_',str_pad(as.character(c(3:13,18:28)),width=3,side='left',pad='0')),
                   ## american indian/alaska native <65
                   paste0('B01001C_',str_pad(as.character(c(3:13,18:28)),width=3,side='left',pad='0')),
                   ## asian <65
                   paste0('B01001D_',str_pad(as.character(c(3:13,18:28)),width=3,side='left',pad='0')),
                   ## native hawaiian/pacific islander
                   paste0('B01001E_',str_pad(as.character(c(3:13,18:28)),width=3,side='left',pad='0')),
                   ## hispanic
                   paste0('B01001I_',str_pad(as.character(c(3:13,18:28)),width=3,side='left',pad='0'))
                 ),
                 state = "MA",year=2012)

## organize ##
ma_acs<-as.data.frame(ma_acs)

## racevar is a variable that tells us what racegroup the row represents ##
ma_acs$racevar<-substr(ma_acs$variable, 1, 7)
## vnum is a variable that tells us the sex and age group the row represents ##
ma_acs$vnum<-as.numeric(substr(ma_acs$variable,nchar(ma_acs$variable)-2,nchar(ma_acs$variable)))

## create a dataset with proper race and sex variables and an age group variable with consistent groupings (called rxsxa) based on the racevar and vnum variables ##
## age ranges ##
agecat<-paste0('Age',c('0-4','5-9','10-14','15-17','18-19','20-24','25-29','30-34','35-44','45-54','55-64'))

rxsxa<-ma_acs[!duplicated(ma_acs[,c('racevar','vnum')]),c('racevar','vnum')]
rxsxa<-rxsxa[order(rxsxa$racevar,rxsxa$vnum),]

realign_fineage<-c(1:5,6,6,6,7,8,9,9,10,10,11,11,11)
rxsxa$age<-c(rep(realign_fineage,2),rep(1:length(agecat),2*(length(unique(rxsxa$racevar))-1)))
rxsxa$sex<-c(rep(c('m','f'),each=length(realign_fineage)),rep(rep(c('m','f'),each=length(agecat)),length(unique(rxsxa$racevar))-1))
rxsxa<-merge(rxsxa,data.frame('age'=1:length(agecat),'agecat'=agecat),by='age')
rxsxa<-merge(rxsxa,data.frame(
  'race'=c('total','white','black','am_indian','asian','pac_islander','hispanic'),
  'racevar'=paste0('B01001',c('_','H','B','C','D','E','I'))),by='racevar')

## merge this new dataset into ma_acs to add proper race, age, sex variables ##
ma_acs<-merge(ma_acs,rxsxa,by=c('racevar','vnum'))

## need to now aggregate population size and MOE where applicable due to the inconsistent age groupings ##
temp1<-aggregate(estimate~GEOID+race+agecat+sex,data=ma_acs,
                 FUN=sum,na.rm=T)
names(temp1)[5]<-'acs_pop'

temp2<-aggregate(moe~GEOID+race+agecat+sex,data=ma_acs,
                 FUN = function(x) 1.96*(sqrt(sum((x/1.645)^2,na.rm=T))))
names(temp2)[5]<-'acs_moe'

ma_acs<-merge(temp1,temp2,by=c('GEOID','race','agecat','sex'))

ma_acs<-ma_acs[order(ma_acs$GEOID),]

## add in the poverty and ICE measures computed from ACS ##
ice<-read_sas(data_file = 'icemeasures_acs_0812_12_12_18.sas7bdat')
ice<-as.data.frame(ice)
ice<-ice[,c('GEOid2','perc_belowpov','ICEracewb','ICEinc','ICEwbinc')]
names(ice)<-c('GEOID','acs_pov_pct','acs_ice_racewb','acs_ice_inc','acs_ice_raceinc')

ma_acs<-merge(ma_acs,ice,by='GEOID')

save(ma_acs,file='acs_data.RData')
rm(list=ls())

####################
## 2. CENSUS DATA ##
####################

## census 2010 data ##
## variable names ##
#v10 <- load_variables(2010, "sf1")
racemg<-data.frame('race'=c('total','white','black','am_indian','asian','pac_islander','hispanic'),'racevar'=paste0('P012',c('0','I','B','C','D','E','H')))
nums65<-str_pad(c(3:19,27:43),width=3,side='left',pad='0')
vnames<-c(paste0('P012',nums65),paste0(rep(racemg$racevar[2:nrow(racemg)],each=length(nums65)),nums65))

##extract using tidycensus ##
ma_ce<-get_decennial(geography = "tract",
                     variables=vnames,
                     state = "MA",year=2010,sumfile='sf1')

## organize ##
ma_ce<-as.data.frame(ma_ce)

## racevar is a variable that tells us what racegroup the row represents ##
ma_ce$racevar<-substr(ma_ce$variable,1,5)
## vnum is a variable that tells us the sex and age group the row represents ##
ma_ce$vnum<-as.numeric(substr(ma_ce$variable,nchar(ma_ce$variable)-2,nchar(ma_ce$variable)))

## create a dataset with proper race and sex variables and an age group variable with consistent groupings (called rxsxa) based on the racevar and vnum variables ##
## age ranges ##
agecat<-paste0('Age',c('0-4','5-9','10-14','15-17','18-19','20-24','25-29','30-34','35-44','45-54','55-64'))

rxsxa<-ma_ce[!duplicated(ma_ce[,c('racevar','vnum')]),c('racevar','vnum')]
rxsxa<-rxsxa[order(rxsxa$racevar,rxsxa$vnum),]

realign_fineage<-c(1:5,6,6,6,7,8,9,9,10,10,11,11,11)
rxsxa$age<-rep(realign_fineage,2*(length(unique(rxsxa$racevar))))
rxsxa$sex<-rep(rep(c('m','f'),each=length(realign_fineage)),length(unique(rxsxa$racevar)))
rxsxa<-merge(rxsxa,data.frame('age'=1:length(agecat),'agecat'=agecat),by='age')
rxsxa<-merge(rxsxa,racemg,by='racevar')

ma_ce<-merge(ma_ce,rxsxa,by=c('racevar','vnum'))

## need to now aggregate population size and MOE where applicable due to the inconsistent age groupings ##
ma_ce<-aggregate(value~GEOID+race+agecat+sex,data=ma_ce,
                 FUN=sum,na.rm=T)
names(ma_ce)[5]<-'ce_pop'

## merge in the ice for race from the census (P003002-P003003)/P001001 ##
ice<-get_decennial(geography = "tract",
                   variables=c('P001001','P003002','P003003'),
                   state = "MA",year=2010,sumfile='sf1',output = 'wide')
ice<-as.data.frame(ice)
ice<-data.frame('GEOID'=ice$GEOID,'ce_ice_racewb'=(ice$P003002-ice$P003003)/ice$P001001)

ma_ce<-merge(ma_ce,ice,by='GEOID')

ma_ce<-ma_ce[order(ma_ce$GEOID),]

save(ma_ce,file='ce_data.RData')

rm(list=ls())

################
## 3. DP DATA ##
################

racemg<-data.frame('race'=c('total','white','black','am_indian','asian','pac_islander','hispanic'),'racevar'=paste0('P012',c('0','I','B','C','D','E','H')))
nums65<-str_pad(c(3:19,27:43),width=3,side='left',pad='0')
vnames<-c(paste0('P0120',nums65),paste0(rep(racemg$racevar[2:nrow(racemg)],each=length(nums65)),nums65))

## download/unzip the MA DP data from https://ciser.cornell.edu/data/data-archive/census-2010-dhc-download-center/ ##

## read in the DP data, takes a few mins to read in ##
ma_dp<-read.csv('MA2010DHC.CSV',stringsAsFactors = F)

## only census tract level info ##
ma_dp<-subset(ma_dp,SUMLEV==140)

## create fips code ##
ma_dp$GEOID<-as.character(paste0(ma_dp$STATE,str_pad(ma_dp$COUNTY,width=3,side='left',pad='0'),str_pad(ma_dp$TRACT,width=6,side='left',pad='0')))

## pull out the variables needed to create ice_racewb ##
ice<-ma_dp[,c('GEOID','P0010001','P0030002','P0030003')]

## subset to only variables we will use for denominator ##
ma_dp<-ma_dp[,c('GEOID',vnames)]

## put in long form and organize ##
ma_dp<-melt(data = ma_dp, id.vars = "GEOID", measure.vars = vnames,factorsAsStrings = T)
ma_dp$variable<-as.character(ma_dp$variable)

## racevar is a variable that tells us what racegroup the row represents ##
ma_dp$racevar<-substr(ma_dp$variable,1,5)
## vnum is a variable that tells us the sex and age group the row represents ##
ma_dp$vnum<-as.numeric(substr(ma_dp$variable,nchar(ma_dp$variable)-2,nchar(ma_dp$variable)))

## create a dataset with proper race and sex variables and an age group variable with consistent groupings (called rxsxa) based on the racevar and vnum variables ##
## age ranges ##
agecat<-paste0('Age',c('0-4','5-9','10-14','15-17','18-19','20-24','25-29','30-34','35-44','45-54','55-64'))

rxsxa<-ma_dp[!duplicated(ma_dp[,c('racevar','vnum')]),c('racevar','vnum')]
rxsxa<-rxsxa[order(rxsxa$racevar,rxsxa$vnum),]

realign_fineage<-c(1:5,6,6,6,7,8,9,9,10,10,11,11,11)
rxsxa$age<-rep(realign_fineage,2*(length(unique(rxsxa$racevar))))
rxsxa$sex<-rep(rep(c('m','f'),each=length(realign_fineage)),length(unique(rxsxa$racevar)))
rxsxa<-merge(rxsxa,data.frame('age'=1:length(agecat),'agecat'=agecat),by='age')
rxsxa<-merge(rxsxa,racemg,by='racevar')

ma_dp<-merge(ma_dp,rxsxa,by=c('racevar','vnum'))

## need to now aggregate population size and MOE where applicable due to the inconsistent age groupings ##
ma_dp<-aggregate(value~GEOID+race+agecat+sex,data=ma_dp,
                 FUN=sum,na.rm=T)
names(ma_dp)[5]<-'dp_pop'

## merge in the ice for race from the dp data (P0030002-P0030003)/P0010001 ##
ice<-data.frame('GEOID'=ice$GEOID,'dp_ice_racewb'=(ice$P0030002-ice$P0030003)/ice$P0010001)

ma_dp<-merge(ma_dp,ice,by='GEOID')

ma_dp<-ma_dp[order(ma_dp$GEOID),]

save(ma_dp,file='dp_data.RData')

rm(list=ls())

#######################################
## 4. merge ACS, census, and dp data ##
#######################################

## read each processed dataset ##
load('acs_data.RData')

load('ce_data.RData')

load('dp_data.RData')

## merge them by fips code, race, age, and sex ##
mg1<-merge(ma_acs,ma_ce,by=c('GEOID','race','agecat','sex'))

adat<-merge(mg1,ma_dp,by=c('GEOID','race','agecat','sex'))

save(adat,file='merged_denom_cov_data.RData')
rm(list=ls())

###################################
## 5. merge in mortality counts  ##
###################################

## make a function to process mortality data ##
process_mort<-function(mort){
  ## remove deaths that are not assigned to a CT or gender or race or ethnicity ##
  mort<-subset(mort,areakey10>0 & !(racecat=='unknown') & !(gencat=='unknow') & !(hispanic==9))
  
  ## subset to premature mortalities (<65) ##
  pmort<-subset(mort,age<65)
  
  ## add age categories ##
  agecat<-paste0('Age',c('0-4','5-9','10-14','15-17','18-19','20-24','25-29','30-34','35-44','45-54','55-64'))
  
  pmort$agecat<-cut(pmort$age,breaks=c(0,5,10,15,18,20,25,30,35,45,55,65),right = F)
  levels(pmort$agecat)<-agecat
  
  ## count deaths within CT, age, sex groups (for total across races) ##
  cpmort<-as.data.frame(table(pmort$areakey10,pmort$agecat,pmort$gencat))
  names(cpmort)<-c('GEOID','agecat','sex','ndeaths')
  cpmort$race<-'total'
  
  ## count deaths within CT, age, sex groups for hispanic ethnicity ##
  pmort_eth<-subset(pmort,hispanic %in% 1:7 | race==15)
  
  cpmort_eth<-as.data.frame(table(pmort_eth$areakey10,pmort_eth$agecat,pmort_eth$gencat))
  names(cpmort_eth)<-c('GEOID','agecat','sex','ndeaths')
  cpmort_eth$race<-'hispanic'
  
  ## subset to only races non-hispanic white, black, asian/pacislander, american indian ##
  pmort$race2<-NA
  pmort$race2[which(pmort$hispanic==0 & pmort$race==1)]<-'white'
  pmort$race2[which(pmort$race==2)]<-'black'
  pmort$race2[which(pmort$race==3)]<-'am_indian'
  pmort$race2[which(pmort$race %in% c(4, 5, 7, 8, 9, 10, 12))]<-'asian'
  pmort$race2[which(pmort$race %in% c(6,11))]<-'pac_islander'
  
  pmort_sub<-subset(pmort,!is.na(race2))
  
  ## count deaths within CT, race, age, sex groups ##
  cpmort_sub<-as.data.frame(table(pmort_sub$areakey10,pmort_sub$race2,pmort_sub$agecat,pmort_sub$gencat))
  names(cpmort_sub)<-c('GEOID','race','agecat','sex','ndeaths')
  
  adat_pmort<-rbind(cpmort,cpmort_sub,cpmort_eth)
  levels(adat_pmort$sex)<-c('f','m')
  
  return(adat_pmort)
}

## read/process 2008-2012 MA mortality data ##
mort5yr<-NULL
for (i in paste0('mort_geo_',c('08','09',10:12),'.csv')){
  temp<-read.csv(i,stringsAsFactors = F,header=T)
  mort5yr<-rbind(mort5yr,temp[,c('areakey10','age','racecat','gencat','race','hispanic')])
}

adat_mort5yr<-process_mort(mort=mort5yr)
names(adat_mort5yr)[grep('ndeaths',names(adat_mort5yr))]<-'ndeaths_5yr'

## merge with the denominator and covariate data ##
load('merged_denom_cov_data.RData')

adat<-merge(adat,adat_mort5yr,by=c('GEOID','race','agecat','sex'),all.x=T)

## missing ndeaths ==> 0 deaths ##
adat$ndeaths_5yr[which(is.na(adat$ndeaths_5yr))]<-0

save(adat,file='merged_pm_denom_cov_data.RData')

rm(list=ls())

#########################
## 6. analyze the data ##
#########################

load("merged_pm_denom_cov_data.Rdata")

## download SEER standard population data from https://seer.cancer.gov/stdpopulations/stdpop.19ages.txt ##

# read in SEER standard population data
seer.age19 <- read_fwf("stdpop.19ages.txt",
                       fwf_widths(c(3,3,8),
                                  c("standard","age","std_raw"))) %>%
  filter(standard=="201") %>%
  mutate(agecat=recode(age,
                       '000'="Age0-4",                   
                       '001'="Age0-4",
                       '002'="Age5-9",
                       '003'="Age10-14",
                       '004'="Age15-19",
                       '005'="Age20-24",
                       '006'="Age25-29",
                       '007'="Age30-34",
                       '008'="Age35-44",
                       '009'="Age35-44",
                       '010'="Age45-54",
                       '011'="Age45-54",
                       '012'="Age55-64",
                       '013'="Age55-64",
                       '014'="Age65-74",
                       '015'="Age65-74",
                       '016'="Age75-84",
                       '017'="Age75-84",
                       '018'="Age85+"),
         std.pop=as.numeric(std_raw)) %>%
  group_by(agecat) %>%
  summarise(std=sum(std.pop))



# pull out ABSMs from adat
absm.data <- adat %>% select(GEOID,acs_pov_pct, acs_ice_racewb, acs_ice_inc, acs_ice_raceinc,
                             ce_pop, ce_ice_racewb, dp_pop, dp_ice_racewb) %>%
  group_by(GEOID) %>%
  summarise(indpov=first(acs_pov_pct),
            ICEinc=first(acs_ice_inc),
            ICEracewb=first(acs_ice_racewb),
            ICEwbinc=first(acs_ice_raceinc),
            ceICEracewb=first(ce_ice_racewb),
            dpICEracewb=first(dp_ice_racewb))

# Get decennial census ICE and also save P001001 for weighting
raw.ice <- get_decennial(geography="tract",variables=c("P001001","P003002","P003003"), 
                         year=2010, state="ma", output="wide")

raw.ice <- mutate(raw.ice, ceICEracewb=(P003002-P003003)/P001001) %>%
  select(GEOID, P001001)

absm.data2 <- inner_join(absm.data, raw.ice, by="GEOID") %>%
  mutate(apINDPOV=case_when(
    0<=indpov & indpov<5 ~ 1,
    5<=indpov & indpov<10 ~ 2,
    10<=indpov& indpov<20 ~ 3,
    20<=indpov& indpov<=100 ~ 4),
    qICEinc=cut(ICEinc, wtd.quantile(ICEinc, weight=P001001,
                                     probs=c(0,0.2,0.4,0.6,0.8,1))),
    qICEracewb=cut(ICEracewb, wtd.quantile(ICEracewb, weight=P001001,
                                           probs=c(0,0.2,0.4,0.6,0.8,1))),
    qICEwbinc=cut(ICEwbinc, wtd.quantile(ICEwbinc, weights=P001001,
                                         probs=c(0,0.2,0.4,0.6,0.8,1))),
    qceICEracewb=cut(ceICEracewb, wtd.quantile(ceICEracewb, weights=P001001,
                                               probs=c(0,0.2,0.4,0.6,0.8,1))),
    qdpICEracewb=cut(dpICEracewb, wtd.quantile(ceICEracewb, weights=P001001,
                                               probs=c(0,0.2,0.4,0.6,0.8,1))))


# output cutpoints for quintile variables
write.csv(rbind(c("qICEracewb (census)", names(table(absm.data2$qceICEracewb))),
                c("qICEracewb (dp)", names(table(absm.data2$qdpICEracewb))),
                c("qICEracewb (acs)", names(table(absm.data2$qICEracewb))),
                c("qICEinc (acs)", names(table(absm.data2$qICEinc))),
                c("qICEwbinc (acs)", names(table(absm.data2$qICEwbinc)))),
          "ICEquintiles_cutpoints.csv")


# adat seems to have an extra age category that doesn't
# exist in the SEER age standard
# so collapse into agecat="Age15-19" 
# then remerge in the absm data
adat2 <- adat %>% mutate(agecat=recode(agecat,
                                       'Age15-17'="Age15-19",
                                       'Age18-19'="Age15-19"),
                         race=recode(race,
                                     'asian' = "api",
                                     'pac_islander' = "api")) %>%
  group_by(GEOID, race, agecat) %>%
  summarise(numdeaths=sum(ndeaths_5yr),
            acs.pop=sum(acs_pop),
            ce.pop = sum(ce_pop),
            dp.pop = sum(dp_pop)) %>%
  left_join(absm.data,by="GEOID") 

# remerge ABSM data
merged.data <- left_join(adat2, absm.data2, by="GEOID")

# this function will do the race aggregation
# Step 1: Aggregate data over GEOID into age strata
#   within race x sex categories. nnn=numerator, ddd=denominator
#   right now this is coded for ce.pop
# Step 2: Create w.iri=weight times age-specific incidence rate
#   Create w.var.iri= weight squared times var(iri)
#   Keep track of weight squared for variance calculation
# Step 3: Group by race and sex so that we can sum over agecat
# Step 4: sum over age cat
# Create ir.std = age standardized incidence rate
# Create var.ir.std = variance of age standardized incidence rate

# note that we apply enquo() to the variable name we supply in popdata
# then we use !!pdata when we call summarise()
# This is some dark magic!
f.raceagg <- function(popdata){
  pdata <- enquo(popdata)
  
  race.agg<- group_by(merged.data, agecat, race) %>%
    summarise(nnn=sum(numdeaths), dd=sum(!!pdata)) %>%
    mutate(ddd=5*dd) %>%
    left_join(seer.age19, by="agecat") %>%
    mutate(w.iri=(std*nnn)/ddd,
           w.var.iri=std^2*(nnn)/ddd^2,
           w.sq=std^2) %>%
    group_by(race) %>%
    summarise(cases=sum(nnn), pop=sum(ddd), sum.w.iri=sum(w.iri),
              sum.w.var.iri=sum(w.var.iri),
              sum.wt = sum(std),
              sum.wt.sq = sum(w.sq)) %>%
    mutate(ir.std = sum.w.iri/sum.wt,
           var.ir.std = sum.w.var.iri/sum.wt.sq,
           dummyMerge=1) %>% ungroup()
  # extract the reference rate for race, which is white
  # and rename rate to ref.ir and variance to ref.var
  ref.race.agg <- race.agg %>%
    filter(race=="white") %>%
    mutate(ref.ir = ir.std, ref.var = var.ir.std, dummyMerge=1) %>%
    select(dummyMerge, ref.ir, ref.var) %>% ungroup()
  
  race.agg2 <- left_join(race.agg, ref.race.agg, by="dummyMerge") %>%
    mutate(ir.std.lo95 = ir.std - 1.96*sqrt(var.ir.std),
           ir.std.up95 = ir.std + 1.96*sqrt(var.ir.std),
           IRR = ir.std / ref.ir,
           varlogirr = (var.ir.std/ir.std^2) + (ref.var/ref.ir^2),
           IRRlo = exp(log(IRR) - 1.96*sqrt(varlogirr)),
           IRRup = exp(log(IRR) + 1.96*sqrt(varlogirr)),
           denSource=paste(enexpr(popdata))) %>%
    mutate(var.log.irr=ifelse(race=="white",NA,varlogirr),
           IRRlo95=ifelse(race=="white",NA,IRRlo),
           IRRup95=ifelse(race=="white",NA,IRRup)) %>%
    select(race, cases, pop, ir.std, var.ir.std, ir.std.lo95, ir.std.up95,
           IRR, var.log.irr,IRRlo95, IRRup95, denSource)
  
  return(race.agg2)
}

race.agg.ce.pop <- f.raceagg(ce.pop)
race.agg.acs.pop <- f.raceagg(acs.pop)                 
race.agg.dp.pop <- f.raceagg(dp.pop)                 

race.results <- bind_rows(race.agg.ce.pop, race.agg.acs.pop, race.agg.dp.pop)
write.csv(race.results, file="race_results.csv")



# this function will do the ABSM aggregation
# note that we apply enquo() to the variable name we supply in popdata
# then we use !!pdata when we call summarise()
# This is some dark magic!
f.absmagg <- function(popdata, absm, refval){
  pdata <- enquo(popdata)
  absm.var <- enquo(absm)
  
  race.agg<- group_by(merged.data, agecat, !!absm.var, race) %>%
    summarise(nnn=sum(numdeaths), dd=sum(!!pdata)) %>%
    mutate(ddd=5*dd) %>%
    left_join(seer.age19, by="agecat") %>%
    mutate(w.iri=(std*nnn)/ddd,
           w.var.iri=std^2*(nnn)/ddd^2,
           w.sq=std^2) %>%
    group_by(!!absm.var,race) %>%
    summarise(cases=sum(nnn), pop=sum(ddd), sum.w.iri=sum(w.iri),
              sum.w.var.iri=sum(w.var.iri),
              sum.wt = sum(std),
              sum.wt.sq = sum(w.sq)) %>%
    mutate(ir.std = sum.w.iri/sum.wt,
           var.ir.std = sum.w.var.iri/sum.wt.sq) %>% ungroup()
  # extract the reference rate for race, which is white
  # and rename rate to ref.ir and variance to ref.var
  
  ref.race.agg <- race.agg %>%
    filter(!!absm.var==refval) %>%
    mutate(ref.ir = ir.std, ref.var = var.ir.std) %>%
    select(race, ref.ir, ref.var) %>% ungroup()
  
  race.agg2 <- left_join(race.agg, ref.race.agg, by="race") %>%
    mutate(ir.std.lo95 = ir.std - 1.96*sqrt(var.ir.std),
           ir.std.up95 = ir.std + 1.96*sqrt(var.ir.std),
           IRR = ir.std / ref.ir,
           varlogirr = (var.ir.std/ir.std^2) + (ref.var/ref.ir^2),
           IRRlo = exp(log(IRR) - 1.96*sqrt(varlogirr)),
           IRRup = exp(log(IRR) + 1.96*sqrt(varlogirr)),
           absmVar = paste(enexpr(absm)),
           absmVal = as.numeric(!!absm.var),
           denSource=paste(enexpr(popdata))) %>%
    mutate(var.log.irr=ifelse(!!absm.var==refval,NA,varlogirr),
           IRRlo95=ifelse(!!absm.var==refval,NA,IRRlo),
           IRRup95=ifelse(!!absm.var==refval,NA,IRRup)) %>%
    select(race, cases, pop, ir.std, var.ir.std, ir.std.lo95, ir.std.up95,
           IRR, var.log.irr,IRRlo95, IRRup95, absmVar, absmVal, denSource)
  return(race.agg2)
}

apINDPOV.ce.pop <- f.absmagg(ce.pop, apINDPOV, 1)
apINDPOV.acs.pop <- f.absmagg(acs.pop, apINDPOV, 1)
apINDPOV.dp.pop <- f.absmagg(dp.pop, apINDPOV, 1)

qICEwbinc.ce.pop <- f.absmagg(ce.pop, qICEwbinc, levels(absm.data2$qICEwbinc)[5])
qICEwbinc.acs.pop <- f.absmagg(acs.pop, qICEwbinc, levels(absm.data2$qICEwbinc)[5])
qICEwbinc.dp.pop <- f.absmagg(dp.pop, qICEwbinc, levels(absm.data2$qICEwbinc)[5])

qICEracewb.ce.pop <- f.absmagg(ce.pop, qICEracewb, levels(absm.data2$qICEracewb)[5])
qICEracewb.acs.pop <- f.absmagg(acs.pop, qICEracewb, levels(absm.data2$qICEracewb)[5])
qICEracewb.dp.pop <- f.absmagg(dp.pop, qICEracewb, levels(absm.data2$qICEracewb)[5])

qICEinc.ce.pop <- f.absmagg(ce.pop, qICEinc, levels(absm.data2$qICEinc)[5])
qICEinc.acs.pop <- f.absmagg(acs.pop, qICEinc, levels(absm.data2$qICEinc)[5])
qICEinc.dp.pop <- f.absmagg(dp.pop, qICEinc, levels(absm.data2$qICEinc)[5])

qceICEracewb.ce.pop <- f.absmagg(ce.pop, qceICEracewb, levels(absm.data2$qceICEracewb)[5])
qceICEracewb.acs.pop <- f.absmagg(acs.pop, qceICEracewb, levels(absm.data2$qceICEracewb)[5])
qceICEracewb.dp.pop <- f.absmagg(dp.pop, qceICEracewb, levels(absm.data2$qceICEracewb)[5])

qdpICEracewb.ce.pop <- f.absmagg(ce.pop, qdpICEracewb, levels(absm.data2$qdpICEracewb)[5])
qdpICEracewb.acs.pop <- f.absmagg(acs.pop, qdpICEracewb, levels(absm.data2$qdpICEracewb)[5])
qdpICEracewb.dp.pop <- f.absmagg(dp.pop, qdpICEracewb, levels(absm.data2$qdpICEracewb)[5])



absm.results <- bind_rows(apINDPOV.ce.pop, apINDPOV.acs.pop, apINDPOV.dp.pop,
                          qICEinc.ce.pop, qICEinc.acs.pop, qICEinc.dp.pop,
                          qICEwbinc.ce.pop, qICEwbinc.acs.pop, qICEwbinc.dp.pop,
                          qICEracewb.ce.pop, qICEracewb.acs.pop, qICEracewb.dp.pop,
                          qceICEracewb.ce.pop, qceICEracewb.acs.pop, qceICEracewb.dp.pop,
                          qdpICEracewb.ce.pop, qdpICEracewb.acs.pop, qdpICEracewb.dp.pop)


write.csv(absm.results, file="absm_results.csv")
