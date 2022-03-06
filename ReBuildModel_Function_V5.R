######################################################################
########rebuild the starting ind for phase II#########################
######################################################################

StartBuild <- function(Ind, n_factors, n_variables){
  
  Start <- as.vector(matrix(as.numeric(unlist(strsplit(Ind,split = ""))), ncol = n_factors*n_variables,byrow = TRUE))
  
  return(Start)
}




######################################################################
#########This function was used to rebuild model form the#############
#########AccList######################################################
ModelReBuild <- function(ModInd, n_factors, n_variables, Data = NULL){
  ind <- matrix(as.numeric(unlist(strsplit(ModInd,split = ""))), n_factors, n_variables,byrow = TRUE)
  model_ind <- as.numeric(unlist(strsplit(ModInd,split = "")))
  if (is.null(names(Data))){
    x <- rep("X",n_variables)
    x <- paste(x,c(1:n_variables),sep = "")
  }else{
    x <- names(Data)
  }
  
  f <- rep("F",n_factors)
  f <- paste(f,c(1:n_factors),sep = "")
  mod <- NULL
  for (i in 1:n_factors) {
    mod_sub <- paste(f[i],"=~",paste(x[which(ind[i,]==1)],collapse = "+"),sep = "")
    mod <- c(mod,mod_sub)
  }
  model_cfa <- paste(mod,collapse = "
                     ")
  ModInd <- as.numeric(unlist(strsplit(ModInd,split = "")))
  return(list(ModelInd = ModInd, CFAModel = model_cfa, ModelInd = model_ind))
}


#######################################################################
#######################################################################
################Esstimate a model from AccList#########################
ModelIndEst <- function(ModInd, Data, n_factors, n_variables, Fit = c(5,2,3)) {
  model <- ModelReBuild(ModInd = ModInd, n_factors = n_factors, n_variables = n_variables, Data = Data)
  resModelIndEst <- ModelEst(Model = model, Data = Data, Fit = Fit)
  return(list(Fit1 = resModelIndEst$Fit1, Fit2 = resModelIndEst$Fit2, Fit3 = resModelIndEst$Fit3, 
              Modelres = resModelIndEst$Modelres))
}



##################################################################################
###########Beyond average method to get model from pheromone level ###############
BeyondModel <- function(PheLevel, n_factors, n_variables){
  PheLevel <- round(PheLevel, 4)
  PheLevel[which(PheLevel < sum(PheLevel)/(n_factors*n_variables + 1))] <- 0
  PheLevel[which(PheLevel != 0)] <- 1
  Ind_BeyondModel <- matrix(PheLevel,n_factors,n_variables,byrow = TRUE)
  x <- rep("X",n_variables)
  x <- paste(x,c(1:n_variables),sep = "")
  f <- rep("F",n_factors)
  f <- paste(f,c(1:n_factors),sep = "")
  Beyondmod <- NULL
  for (i in 1:n_factors) {
    Beyondmod_sub <- paste(f[i],"=~",paste(x[which(Ind_BeyondModel[i,]==1)],collapse = "+"),sep = "")
    Beyondmod <- c(Beyondmod,Beyondmod_sub)
  }
  BeyondModelRun <- paste(Beyondmod,collapse = "
                     ")
  BeyondModelShow <- paste(Beyondmod,collapse = "   ")
  return(list(BeyondModelRun = BeyondModelRun, BeyondModelShow = BeyondModelShow))
}

BeyondModelEst <- function(PheLevel, Data, n_factors, n_variables, Fit = c(5,2,3)){
  beyondmodel <- BeyondModel(PheLevel = PheLevel, n_factors = n_factors, n_variables = n_variables)
  beyondmodelres <- OptModelEst(Model = beyondmodel$BeyondModelRun, Data = Data, Fit = Fit)
  return(list(Fit1 = beyondmodelres$Fit1, Fit2 = beyondmodelres$Fit2, Fit3 = beyondmodelres$Fit3, 
              Modelres = beyondmodelres$Modelres))
}

