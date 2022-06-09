###########################################################################
# LOAD DATA FUNCTION
# An epidemic model for SARS-CoV-2 with self-adaptive containment measures #
# S. Marchetti, A. Borin, F.P. Conteduca et al., 2022

LoadData <- function(args){
  colNames <- c("Regione",'Giorno','ClasseEta','PopT','Deceased',
                'Incidenza','Incid_Variant1','NewHosp','S','Iw','Ie','Symp','H','ICU',
                'Rwe','Vacc1_1','Vacc2_1','Vacc1_2','Vacc2_2','BTIw','BTIe','BTIsymp',
                'BTIH','BTICU','Rbtiwe','00')
  out <- read.table(paste0("output_",args),sep="\t",col.names = colNames)[,-length(colNames)]
  out[out<0] <- 0
  if(args!='0'){
    out0 <- read.table(paste0("Outs/outRes_",args-1))[,-length(colNames)]
    nn <- nrow(out0)
    names(out0) <- names(out)
    out0$Giorno <- 0
    out <- rbind(out0[1:(nn-14*5*21),],out)
  }
  out$Giorno <- rep(seq(as.Date('2021-06-07'),by='day',length.out = nrow(out)/(21*5)),each=21*5)
  
  # From total to reported cases:
  ratio <- HistCasi[HistCasi$data==max(HistCasi$data),]
  ratio <- ratio[order(ratio$codice_regione,ratio$age_group,decreasing = F),'ratio']
  out$IncidenzaRilevati <- out$Incidenza/rep(ratio,length(unique(out$Giorno)))
  # Save
  write.table(out,paste0("Outs/outRes_",args),row.names = F,col.names=F)
  # remove out
  file.remove(paste0("output_",args))
  file.remove(paste0("Outs/outRes_",args-1))
  return(out)
}
