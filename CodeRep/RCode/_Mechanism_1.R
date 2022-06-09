# MECHANISM 1: OCCUPANCY RATES

# From July 23, 2021 to March 31, 2022
# See "S2 Appendix. The tier system in Italy" for details


CalcZona <- function(curr_pol, regione, weeklyInc, HospInd, dummies, Scenario=NULL, risk=NULL, IndicatorInc=NULL){
  current <- curr_pol[curr_pol[,1]==regione,3:5]
  dummbin <- ifelse(dummies[dummies[,1]==regione,2]==dummies[dummies[,1]==regione,3] &
                      dummies[dummies[,1]==regione,2]==1,1,0)
  
  # Color consistency: at least two weeks in orange and red
  if((current[3] %in% 2:3) & current[2]!=current[3]){return(list(current[3],0))}

  if((HospInd==0 | weeklyInc<50) & dummbin==1){
    return(list(0,(weeklyInc<50)*1))
    }
    else{
      if(weeklyInc<150){return(list(max(current[3]-1,1),(weeklyInc<50)*1))}
      else{
        return(list(max(current[3]-1,HospInd),0))
      }
  }
}
