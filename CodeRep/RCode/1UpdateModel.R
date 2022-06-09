rm(list=ls())
options(warn=-1) #0 to turn on
verbose <- T

setwd("/home/group/main/COVID/private/simulazioni/SIR_AgeDep/SIR_MAY/CodeRep/")

argss <- as.numeric(commandArgs(TRUE))

# Main global variables:
args <- argss[1]
end_arg <- argss[2]
framework <- argss[3]
mechanism <- argss[4]

suppressMessages(library(zoo))
suppressMessages(library(EpiEstim))
source('RCode/_LoadData.R')
source(paste0('RCode/_Mechanism_',mechanism,'.R'))

#############################
# Create/Update directories #
#############################

old_input <- paste0('input_',args,'/')
new_input <- paste0('input_',args+1)

if(!dir.exists(new_input)){
  dir.create(new_input)}

scenn <- paste0('Outs/Mechanism_',mechanism,'_Framework_',framework)
if(!dir.exists(scenn)){
  dir.create(scenn,recursive = T)}

# Read observations from previous week:
regions <- read.table('input/RegCodes',header=T)
regions <- data.frame(regions,Pop=read.table('input_0/popR')$V1)

HistCasi <- haven::read_dta('input/casi_effettivi_rilevati.dta')
HistCasi <- merge(HistCasi,regions,by.x = 'region',by.y='Code',all.x=T)
HistCasi$data <- as.Date(HistCasi$data)

zonas <- c('White','Yellow','Orange','Red') # Available tiers
nagecl <- 5 # 0-12, 13-19, 20-64, 64-79, 80+
jj <- 1 # For regional iterations
popR <- read.table(paste0(old_input,'/popR'))$V1 # Regional population (aggregate)

# Load last available data
out <- LoadData(args)
out$Giorno <- as.Date(out$Giorno)

days <- unique(out$Giorno)

# Define parameters for Rt computation
if(mechanism>2){
  # Estimate Rt
  library(EpiEstim)
  shape.est <- 2.46
  rate.est <- .265
  N <- 300
  serial.int <- dgamma(0:N, shape=shape.est, rate=rate.est) 
  SI <- (serial.int/sum(serial.int))
}

#######################
# Update Temperatures #
#######################

ttemp <- readxl::read_excel('input/temp_proj_jun_Dec_2021.xlsx')
ttemp <- merge(ttemp,regions,by.x='region',by.y='Code',all.x=T)
ttemp <- ttemp[order(ttemp$date,ttemp$codice_regione),]
ttemp$date <- as.Date(ttemp$date)

medbyreg <- aggregate(ttemp$diff_temp,by=list(ttemp$codice_regione),median)
names(medbyreg) <- c('codice_regione','median')
ttemp <- merge(ttemp,medbyreg,by='codice_regione',all.x=T)
ttemp$newdiff <- (ttemp$diff_temp-ttemp$median)*10
ttemp <- ttemp[order(ttemp$date,ttemp$codice_regione),]

thisweek <- seq(as.Date('2021-06-14'),length.out = as.numeric(end_arg)+3,by='week')
thisweek <- as.Date(thisweek[(args+1):(args+3)])

Temperaturesw <- Temperaturesen <- c()
for(i in unique(regions$codice_regione)){
  Temperaturesw <- rbind(Temperaturesw,-.003*ttemp$diff_temp[ttemp$codice_regione==i & 
                                                               (ttemp$date %in% thisweek)])
  Temperaturesen <- rbind(Temperaturesen,-.015*ttemp$diff_temp[ttemp$codice_regione==i & 
                                                                 (ttemp$date %in% thisweek)])
}

write.table(Temperaturesen,paste0(new_input,'/Temperatures1'),row.names = F,col.names = F)
write.table(Temperaturesw,paste0(new_input,'/Temperatures2'),row.names = F,col.names = F)

###############################
# Prepare next initialization #
###############################

#Current policy matrix
curr_pol <- cbind(seq(1,22)[-4],as.matrix(read.table(paste0(old_input,'/pol'),header=F)[1:21,]))

# Proportion first dose vaccinees
propfirstdose <- aggregate(out[out$Giorno==days[length(days-14)],c('Vacc1_1','Vacc2_1','Vacc1_2','Vacc2_2')],
                           by=list(out$ClasseEta[out$Giorno==days[6]]),sum)
propfirstdose <- apply(propfirstdose[,2:3],1,sum)/apply(propfirstdose[,-1],1,sum)
propfirstdose[is.nan(propfirstdose)] <- 1
write.table(propfirstdose,paste0(new_input,'/FirstDose'),row.names = F,col.names = F)

# Proportion breakthrough infections
dfg <- 7
propBTI <- apply(aggregate(out[out$Giorno==days[length(days)-dfg],c('Symp','H','ICU','BTIsymp','BTIH','BTICU')],
                           by=list(out$ClasseEta[out$Giorno==days[length(days)-dfg]]),sum)[,-1],1,sum)
propBTI <- apply(aggregate(out[out$Giorno==days[length(days)-dfg],c('BTIsymp','BTIH','BTICU')],
                           by=list(out$ClasseEta[out$Giorno==days[length(days)-dfg]]),sum)[,-1],1,sum)/propBTI
write.table(propBTI,paste0(new_input,'/propBTI'),row.names = F,col.names = F)

#########################
# Derive regional tiers #
#########################
infR <- RR <-Speo <-  prevvariant1 <- policy <- School <-vcomp1 <- vcomp2 <-RemB <- icu <- H <-
  dummy <- BTIinf <- flag <- c()
for(reg in 1:22){
  if(reg==4){next}
  
  # regional outcomes
  out_ <- out[out$Regione==reg,]

  ##############################################
  #INCIDENCE
  ##############################################
  # Daily regional incidence
  num0 <- aggregate(out_$IncidenzaRilevati, by=list(out_$Giorno), FUN=sum)
  nr <- rev(num0$Group.1)[1]
  num <- sum(num0$x[num0$Group.1>=(nr-21) & num0$Group.1<(nr-10)]) 
  
  # Weekly Incidence
  weeklyInc <- num*1e5/regions$Pop[regions$codice_regione==reg]
  

    if(mechanism %in% 3:4){
      IndicatorInc <- ifelse(weeklyInc<50,1,0)
    }

    inc0 <- aggregate(out_$IncidenzaRilevati,by=list(out_$Giorno),'sum')

    # Extend series if args<2
    if(args<2){
      inc <- HistCasi[HistCasi$codice_regione==reg & HistCasi$data<min(inc0$Group.1),]
      inc <- aggregate(inc$casi_effettivi/inc$ratio,by=list(inc$data),'sum')
      inc <- rbind(inc,inc0)
    }
    else{inc <- inc0}
    
    # Number of (official) cases over the last 5 days
    n_ <- nrow(inc)

    cases5days <- sum(inc$x[(n_-(10+5-1)):(n_-10)])
    # Occupancy rates:
    # https://www.agenas.gov.it/covid19/web/index.php?r=site%2Ftab2
    
    PL <- read.table('input/PL_hosp_ICU_upd',sep='\t',header=T)
    
    # MA -> "hosp_rate": occupation rate of MA (last available)
    totH <- aggregate(out_[c('H','BTIH')],by=list(out_$Giorno),sum)
    totH$x <- apply(totH[,c('H','BTIH')],1,sum)
    
    nr <- max(totH$Group.1)
    hosp_rate <- round(totH$x[totH$Group.1==(nr-3)]/PL$PLH[PL$Regioni==reg],2)
    
    # ICU -> "icu_rate": occupation rate of ICU (last available)
    totICU <- aggregate(out_[,c('ICU','BTICU')],by=list(out_$Giorno),sum)
    totICU$x <- apply(totICU[,c('ICU','BTICU')],1,sum)
    
    icu_rate <- round(totICU$x[totICU$Group.1==(nr-3)]/
                        (PL$PLICU[PL$Regioni==reg]+PL$PLICUpot[PL$Regioni==reg]),2)
    
    if(mechanism <3){
    HospInd <- ifelse(hosp_rate>.4 & icu_rate>.3,3,
                             ifelse(hosp_rate<=.15 | icu_rate<=.1,0,
                                    ifelse(hosp_rate<=.3 | icu_rate<=.2,1,2)))
    }

    if(mechanism %in% 3:4){
      ########################################################
      # SCENARIO: Mechanism 3 and 4 only
      # Scenario 1, Rt<1;
      # Scenario 2, 1<=Rt<1.25;
      # Scenario 3, 1,25<=Rt<1.5;
      # Scenario 4, Rt>=1.5
        
      if(mechanism==4){ # Replace series with new hospitalizations if mechanism is Rt-Hospital Admissions
          inc_ <- aggregate(out_$NewHosp,by=list(out_$Giorno),'sum')
          inc_ <- inc_[1:(nrow(inc_)-3),]
        }
      else{
        inc_ <- inc[(n_-34):(n_-14),]
        }
          
      t_start <- seq(2, nrow(inc_)-6)
      t_end <- t_start+6
      Rt <- estimate_R(incid = inc_$x,method = "non_parametric_si",
                       config = make_config(list(si_distr = SI,t_start=t_start,
                                                 t_end=t_end)))
      m <- nrow(Rt$R)
      min_rt_outlook <- Rt$R$`Mean(R)`[(m-6):m]
      
      Scenario <- ifelse(min(min_rt_outlook)<1,1,
                         ifelse(min(min_rt_outlook)<1.25,2,
                                ifelse(min(min_rt_outlook)<1.5,3,4)))
      

    # ################################################
    # # RISK
    # # Risk is derived from 'prob' and 'impact' variables
    # 
    # #############################
    # # PROB
    # # 4 levels: (1) very low, (2) low, (3) moderate, (4) high
    # # Levels are derived based on:

    # # a. "cases5days": absolute number of official cases over the last 5 days
    
    # # b. "delta_wc": variation in the number of weekly (official) cases:
    delta_wc <- sum(inc$x[(n_-(10+7-1)):(n_-10)])-sum(inc$x[(n_-(10+14-1)):(n_-(10+7))])
    
    # # c. "rt_outlook": average Rt value for the reference week;
    rt_outlook <- mean(min_rt_outlook)
    
    # d. "icu_rate": occupation rate of ICU (last available)
    
    prob <- ifelse(cases5days<1,1,
                    ifelse(rt_outlook<1 & delta_wc<0,2,
                           ifelse(icu_rate>=.4 & weeklyInc>=330 & rt_outlook>=1.4,4,3)))
    
    ############################################
    # # IMPACT
    # # 4 levels: (1) very low, (2) low, (3) moderate, (4) high
    # # Levels are derived based on:
    #
    # a. Official cases over the last 5 days among elderly individuals
    
    out__ <- out_[out_$ClasseEta>2 & out_$Giorno<=(nr-10) & out_$Giorno>(nr-(10+5-1)),]
    cases5days_over65 <- sum(aggregate(out__$IncidenzaRilevati, by=list(out__$Giorno), FUN=sum)$x)
    
    # b. "icu_rate": occupation rate of ICU (last available)
    
    # c. "hosp_rate": occupation rate of MA (last available)
    
    impact <- ifelse(cases5days_over65<1,1,
                      ifelse(hosp_rate<.4 & icu_rate<.3,2,
                           3))

    # RISK:
    (risk <- ifelse((impact==1 & prob>1 & prob<4) | (impact==2 & prob<3) | (impact==3 & prob==1),2,
                    ifelse((impact==1 & prob==4) | (impact==2 & prob>2) | (impact==3 & prob>1 & prob<4), 3,
                           ifelse((impact==3 & prob==4),4,1))))
    
    }
    #########################
    # Derive next tier:
    if(mechanism<3){
      Dummies50 <- read.table(paste0(old_input,'/flag'),header=T)
      outCalczona <- CalcZona(curr_pol=curr_pol,regione=reg,weeklyInc=weeklyInc,HospInd=HospInd,dummies=Dummies50)
    }
    else{
      outCalczona <- CalcZona(curr_pol=curr_pol,regione=reg,Scenario = Scenario, risk = risk, IndicatorInc = IndicatorInc,weeklyInc = weeklyInc)
    }

    region_class <- outCalczona[[1]]
    flag <- c(flag,outCalczona[[2]])
    if(verbose){
      nr <- max(out_$Giorno)-10
      if(reg==1){
      cat('\n',as.character(nr))}
      cat('\n -> ',regions$reg[regions$codice_regione==reg],': ',zonas[region_class+1],'\n')
    }

  policy <- c(policy,region_class)
  
  # Regional values to initialize the new iteration
  # Infectious
  infR <- rbind(infR,apply(out_[out_$Giorno==(nr-10),c('Iw','Ie','Symp')],1,sum)) #H, ICU
  # Breakthrough infectious
  BTIinf <- rbind(BTIinf,apply(out_[out_$Giorno==(nr-10),c('BTIw','BTIe','BTIsymp')],1,sum))
  # Recovered
  RR <- rbind(RR,out_$Rwe[out_$Giorno==(nr-10)])
  # Prevalence of variant 1
  prevvariant1 <- c(prevvariant1,sum(out_[out_$Giorno==(nr-10),c('Ie','BTIe')])/sum(out_[out_$Giorno==(nr-10),c('Iw','Ie','BTIw','BTIe')]))
  # Vaccinated compartments
  vcomp1 <- rbind(vcomp1,apply(out_[out_$Giorno==(nr-10),c('Vacc1_1','Vacc1_2')],1,sum))
  vcomp2 <- rbind(vcomp2,apply(out_[out_$Giorno==(nr-10),c('Vacc2_1','Vacc2_2')],1,sum))
  # Regional population
  popR[jj] <- popR[jj]-sum(out_$Deceased[out_$Giorno<=(nr-10)])
  # Hospitalized population
  icu <- rbind(icu,out_$ICU[out_$Giorno==(nr-10)])
  H <- rbind(H,out_$H[out_$Giorno==(nr-10)])
  # Vaccinated & Recovered
  RemB <- rbind(RemB,out_$Rbtiwe[out_$Giorno==(nr-10)])
  # To become infectious
  Speo <- rbind(Speo,out_$Symp[out_$Giorno==(nr-10)])
  
  # Schools re-opening
  iniz_scuole_reg <- read.csv('input/inizio_scuole.csv',sep='\t',header=T)
  
  SchoolClose <- ifelse(args<iniz_scuole_reg$Sett[iniz_scuole_reg$codice_regione==reg] | args>28,.6,
                        .75)
  
  School <- c(School,SchoolClose)
  jj <- jj+1
  }

flag_2 <- read.table(paste0(old_input,'/flag'),header=T)$flag_lag1
flag <- data.frame('Regione'=regions$codice_regione,'flag_lag2'=flag_2,'flag_lag1'=flag)
write.table(flag,paste0(new_input,'/flag'),row.names = F)

pol <- cbind(curr_pol[,-c(1:2)],policy)
# Exceptions European Football, Digital Certificate
imprmob <- ifelse(args>2 & args<6,1.07,
                  ifelse(args>17,1.05,1))

pol <- rbind(pol,rep(imprmob,ncol(pol)))
write.table(pol,paste0(new_input,'/pol'),row.names = F,col.names = F)

# PopUpd
write.table(popR,paste0(new_input,'/popR'),row.names = F,col.names = F)
# Sypeople
write.table(Speo,paste0(new_input,'/Speople'),row.names = F,col.names = F)
# ICUUpd
write.table(icu,paste0(new_input,'/ICUpeople'),row.names = F,col.names = F)
# HUpd
write.table(H,paste0(new_input,'/Hpeople'),row.names = F,col.names = F)
# Inf
write.table(infR,paste0(new_input,'/inf'),row.names = F,col.names = F)
# BTInf
write.table(BTIinf,paste0(new_input,'/BTInf'),row.names = F,col.names = F)
# RemB
write.table(RemB,paste0(new_input,'/RemB'),row.names = F,col.names = F)
# R
write.table(RR,paste0(new_input,'/serop'),row.names = F,col.names = F)
# preveng
write.table(prevvariant1,paste0(new_input,'/prevvariant1'),row.names = F,col.names = F)
# School


if(args=='3'){
    # Release of social restrictions
    ppp <- read.table('input/mobil')$V1
    ppp[1] <- .95
    ppp[2]<-.9
    write.table(cbind(ppp,ppp),'input/mobil',row.names = F,col.names = F)
  }
  if(args=='7'){
    # Digital COVID-19 certificate
    ppp <- read.table('input/mobil')$V1
    write.table(cbind(ppp,c(.9,.9,ppp[3:length(ppp)])),'input/mobil',row.names = F,col.names = F)
  }
  if(args=='11'){
    # Tightening of social restrictions
    ppp <- read.table('input/mobil')$V1
    ppp[1] <- .75
    ppp[2]<-.75
    pp0 <- read.table('input//mobil_orig')$V1
    write.table(cbind(ppp,ppp),'input/mobil',row.names = F,col.names = F)
  }
  if(args=='18'){
    # Regime change for workplaces
    ppp <- read.table('input/mobil')$V1
    ppp[1] <- 1
    ppp[2] <- .98
    
    write.table(cbind(ppp,c(ppp[2:length(ppp)],ppp[length(ppp)])),'input/mobil',row.names = F,col.names = F)
  }

Sc <- read.table(paste0(old_input,'/School'))
write.table(cbind(Sc[,2:3],School),paste0(new_input,'/School'),row.names = F,col.names = F)
# Vcomp1
write.table(vcomp1,paste0(new_input,'/vcomp1'),row.names = F,col.names = F)
# Vcomp2
write.table(vcomp2,paste0(new_input,'/vcomp2'),row.names = F,col.names = F)
# Vaccine:
suppressMessages(rollout <- readxl::read_excel('input/COVERAGES.xlsx',sheet = framework))

rollout$fornitore_code <- 0
rollout$fornitore_code[rollout$fornitore=="AstraZeneca/Janssen"] <- 1
rollout$fornitore_code[rollout$fornitore=="AstraZeneca/Janssen(BOOSTER)"] <- 2

names(rollout)[ncol(rollout)-3] <- 'efficacy'
cov1 <- cov2 <- cov3 <- matrix(0,nrow=5,ncol=21)
effcov1 <- effcov2 <- matrix(0,nrow=21,ncol=5)
curr_year <- lubridate::year(max(out$Giorno)-7)
curr_month <- lubridate::month(max(out$Giorno)-7)

for(a in 0:4){
  j <- 1
  for(r in sort(unique(rollout$codice_regione))){
    cand1 <- rollout[rollout$mese==curr_month& rollout$anno==curr_year & rollout$age_group==(a+1) & 
                       rollout$codice_regione==r & rollout$fornitore_code==0,]
    cand2 <- rollout[rollout$mese==curr_month & rollout$anno==curr_year & rollout$age_group==(a+1) & 
                       rollout$codice_regione==r & rollout$fornitore_code==1,]
    cand3 <- rollout[rollout$mese==curr_month & rollout$anno==curr_year & rollout$age_group==(a+1) & 
                       rollout$codice_regione==r & rollout$fornitore_code==2,]
    candeff1 <- rollout$efficacy[rollout$mese==curr_month & rollout$anno==curr_year & rollout$age_group==(a+1) & 
                                   rollout$codice_regione==r & rollout$fornitore_code==0]
    candeff2 <- rollout$efficacy[rollout$mese==curr_month & rollout$anno==curr_year & rollout$age_group==(a+1) & 
                                   rollout$codice_regione==r & rollout$fornitore_code==1]
    cov1[a+1,j] <- ifelse(nrow(cand1)==1,cand1$coverage_of_the_age_group,0)
    cov2[a+1,j] <- ifelse(nrow(cand2)==1,cand2$coverage_of_the_age_group,0)
    cov3[a+1,j] <- ifelse(nrow(cand3)==1,cand3$coverage_of_the_age_group,0)
    effcov1[j,a+1] <- candeff1
    effcov2[j,a+1] <- candeff2
    j <- j+1
  }
}

write.table(cov1,paste0(new_input,'/cov1'),row.names = F,col.names = F)
write.table(effcov1,paste0(new_input,'/efficacy1'),row.names = F,col.names = F)

write.table(cov2,paste0(new_input,'/cov2'),row.names = F,col.names = F)
write.table(effcov2,paste0(new_input,'/efficacy2'),row.names = F,col.names = F)

write.table(cov3,paste0(new_input,'/cov3'),row.names = F,col.names = F)

# Final Output
if(args==end_arg){
  nr <- max(out$Giorno)
  write.table(out[out$Giorno<=(nr-10),],paste0(scenn,"/OUT"),row.names = F)
}

######################### End