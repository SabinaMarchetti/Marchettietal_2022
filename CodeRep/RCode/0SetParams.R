###############################################################################
################## SET PARAMETERS FOR ILLUSTRATIVE SCENARIOS ##################
###############################################################################

argss <- as.numeric(commandArgs(TRUE))

framework <- argss[1]
mechanism <- argss[2]

cat('\n Framework: ',c('Baseline','Optimistic','Pessimistic')[framework],
    '\n Policy Mechanism: ',
    c('Occupancy Rates','Incidence','Rt - New Positives','Rt - Hospital Admissions')[mechanism])

setwd('/home/group/main/COVID/private/simulazioni/SIR_AgeDep/SIR_MAY/CodeRep/')
allparams <- read.table('input/framework')[,framework]

# Delta prevalence
preveng <- data.frame(readxl::read_excel('input/AllSimuls_.xlsx',sheet='Prev_delta'))
write.table(preveng$Preveng,'input_0/prevvariant1',row.names = F,col.names = F)

# Temperatures
regions <- read.table('input/RegCodes',header=T)
regions <- data.frame(regions,Pop=read.table('input/popR')$V1)
ttemp <- readxl::read_excel('input/temp_proj_jun_Dec_2021.xlsx')
ttemp <- merge(ttemp,regions,by.x='region',by.y='Code',all.x=T)
ttemp <- ttemp[order(ttemp$date,ttemp$codice_regione),]

ttemp$date <- as.Date(ttemp$date)

thisweek <- as.Date(seq(as.Date('2021-06-07'),length.out = 3,by='week'))

Temperaturesen <- Temperaturesw <- c()
for(i in unique(regions$codice_regione)){
  Temperaturesw <- rbind(Temperaturesw,-.003*
                           ttemp$diff_temp[ttemp$codice_regione==i & (ttemp$date %in% thisweek)])
  Temperaturesen <- rbind(Temperaturesen,-.015*ttemp$diff_temp[ttemp$codice_regione==i & (ttemp$date %in% thisweek)])
}

write.table(Temperaturesen,'input_0/Temperatures1',row.names = F,col.names = F)
write.table(Temperaturesw,'input_0/Temperatures2',row.names = F,col.names = F)

suppressMessages(rollout <- read.csv('input/coverages.csv'))

rollout$fornitore_code <- 0
rollout$fornitore_code[rollout$fornitore=="AstraZeneca/Janssen"] <- 1
rollout$fornitore_code[rollout$fornitore=="AstraZeneca/Janssen(BOOSTER)"] <- 2
names(rollout)[ncol(rollout)-3] <- 'efficacy'

cov1 <- cov2 <- cov3 <- matrix(0,nrow=5,ncol=21)
effcov1 <- effcov2 <- matrix(0,nrow=21,ncol=5)
curr_year <- 2021
curr_month <- 6
for(a in 0:4){
  j <- 1
  for(r in sort(unique(rollout$codice_regione))){
    
    cand1 <- rollout[rollout$mese==curr_month & rollout$anno==curr_year& rollout$age_group==(a+1) & 
                       rollout$codice_regione==r & rollout$fornitore_code==0,]
    cand2 <- rollout[rollout$mese==curr_month & rollout$anno==curr_year& rollout$age_group==(a+1) & 
                       rollout$codice_regione==r & rollout$fornitore_code==1,]
    cand3 <- rollout[rollout$mese==curr_month & rollout$anno==curr_year& rollout$age_group==(a+1) & 
                       rollout$codice_regione==r & rollout$fornitore_code==2,]
    candeff1 <- rollout$efficacy[rollout$mese==curr_month & rollout$anno==curr_year& rollout$age_group==(a+1) & 
                                   rollout$codice_regione==r & rollout$fornitore_code==0]
    candeff2 <- rollout$efficacy[rollout$mese==curr_month& rollout$anno==curr_year & rollout$age_group==(a+1) & 
                                   rollout$codice_regione==r & rollout$fornitore_code==1]
    cov1[a+1,j] <- ifelse(nrow(cand1)==1,cand1$coverage_of_the_age_group,0)
    cov2[a+1,j] <- ifelse(nrow(cand2)==1,cand2$coverage_of_the_age_group,0)
    cov3[a+1,j] <- ifelse(nrow(cand3)==1,cand3$coverage_of_the_age_group,0)
    effcov1[j,a+1] <- candeff1
    effcov2[j,a+1] <- candeff2
    j <- j+1
  }
}

write.table(cov1,paste0('input_0/cov1'),row.names = F,col.names = F)
write.table(effcov1,paste0('input_0/efficacy1'),row.names = F,col.names = F)

write.table(cov2,paste0('input_0/cov2'),row.names = F,col.names = F)
write.table(effcov2,paste0('input_0/efficacy2'),row.names = F,col.names = F)

write.table(cov3,paste0('input_0/cov3'),row.names = F,col.names = F)

mobil <- read.table('input/mobil_orig')$V1
write.table(cbind(mobil,mobil),'input/mobil',row.names = F,col.names = F)

cat('\n Parameters set. Simulation starts...\n')
