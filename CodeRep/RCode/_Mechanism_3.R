# MECHANISM 3: RT-NEW POSITIVES

# Decree January, 14 2021 and February 23, 2021 (as revised on March, 2021)
# See "S2 Appendix. The tier system in Italy" for details

CalcZona <- function(curr_pol, regione, Scenario, risk, IndicatorInc, weeklyInc, HospInd=NULL, dummies=NULL){

  current <- curr_pol[curr_pol[,1]==regione,3:5]
  
  # Shortcut: RED if incidence > 250x100,000 (March update)
  if(weeklyInc>250){
    return(list(3,1))
  }
  # Color consistency: at least two weeks in orange and red
  if((current[3] %in% 2:3) & current[2]!=current[3]){return(list(current[3],0))}

  # Rt<1
  if(Scenario==1){
    if(risk<3 & IndicatorInc==1 & sum(current)<4){return(list(0,0))}
    else{
    if(risk==4){return(list(2,0))}
      else{return(list(1,0))}
  }
  }
  # 1 < Rt < 1.25
  if(Scenario==2){
    if(IndicatorInc==0 & risk>2){return(list(2,0))}
    else{return(list(1,0))}
  }
  # 1.25 < Rt < 1.5
  if(Scenario==3){
    if(IndicatorInc==1 & risk>3){return(list(2,0))}
    else{
      if(IndicatorInc==0 & risk>2){return(list(3,0))}
      else{return(list(1,0))}
    }
  }
  # Rt > 1.5
  if(Scenario==4){
    if((IndicatorInc==1 & risk==4) | (IndicatorInc==0 & risk>2)){return(list(3,0))}
    else{return(list(1,0))}
  }
}
