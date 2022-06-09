# MECHANISM 2: INCIDENCE

# Decree May 18, 2021
# See "S2 Appendix. The tier system in Italy" for details


CalcZona <- function(curr_pol, regione, weeklyInc, HospInd, dummies, Scenario=NULL, risk=NULL, IndicatorInc=NULL){
  
  current <- curr_pol[curr_pol[,1]==regione,3:5]
  dummbin <- ifelse(dummies[dummies[,1]==regione,2]==dummies[dummies[,1]==regione,3] &
                      dummies[dummies[,1]==regione,2]==1,1,0)
  
  # Color consistency: at least two weeks in orange and red
  if((current[3] %in% 2:3) & current[2]!=current[3]){return(list(current[3],0))}
  
  # WHITE
  if(weeklyInc<50){
    if(dummbin==0){
      return(list(max(current[3]-1,1),1))
    }
    else{
      return(list(max(current[3]-1,0),1))
    }
  }
  # Yellow
  if((weeklyInc>=150 & weeklyInc<250 & HospInd==0) | (weeklyInc<150 & weeklyInc>=50)){
    return(list(max(current[3]-1,1),0))
  }
  # Orange
  if((weeklyInc>=150 & weeklyInc<250 & HospInd==1)){
    return(list(2,0))
  }
  # Red
  if(weeklyInc>=250 | (weeklyInc>=150 & weeklyInc<250 & HospInd==2)){
    return(list(3,0))
  }
}