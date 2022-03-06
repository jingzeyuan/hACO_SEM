#####################################################################
####################Build CFA model##################################
#####################################################################


###Incorporation orior knowledge###***
#build model probability 
ProbBuild <- function(n_factors, n_variables,loaded,unloaded){
  
  baseMatrix <- matrix(rep(0,n_variables*n_factors),n_factors,n_variables)
  loadMatrix <- baseMatrix
  for (i in 1:length(loaded)){
    loadMatrix[i,as.numeric(unlist(loaded[i]))] <- 1
  }
  unloadMatrix <- baseMatrix
  for (i in 1:length(unloaded)){
    unloadMatrix[i,as.numeric(unlist(unloaded[i]))] <- 1
  }
  if (is.element(n_factors,colSums(unloadMatrix))){
    stop(print("External misspecification"))
  }
  if (is.element(2,loadMatrix+unloadMatrix)){
    stop(print("Contradiction Setting"))
  }  
  
  pre_loads <- list(load = loadMatrix, unload = unloadMatrix)
  
  prob <- matrix(rep(0.55,n_variables*n_factors),n_factors,n_variables,byrow = TRUE)
  prob[which(pre_loads$load == 1)] <- 1.0
  prob[which(pre_loads$unload == 1)] <- 0.0
  return(prob)
} 

#build model indicators to build model 
IndBuild <- function(n_factors, n_variables, probability){
  
  ind <- rbinom(n_variables*n_factors,1,as.vector(t(probability)))
  ind_use <- matrix(ind,n_factors,n_variables,byrow = TRUE)
  repair_factor <- n_factors
  #if there is unused variable 
  #repair the model from bottom to top
  while (is.element(0,colSums(ind_use))){
    repair <- which(colSums(ind_use) == 0)
    repair_line <- which(probability[repair_factor,] != 0)
    repair_element <-  repair_line[which(repair_line %in% repair)]
    ind_use[repair_factor,][repair_element] <- 1
    repair_factor <- repair_factor - 1
  }
  ind <- c(t(ind_use))
  return(list(test = ind,ind_use = ind_use))
}

################################
#function to build model
#########################
ModelBuilder <- function(Data, n_factors, n_variables, probability, BanList, CutChange = 0, ChangeRun = ChangeRun){
  Data <- as.data.frame(Data)
  ind <- IndBuild(n_factors, n_variables, probability)
  model_preTest <- paste(ind$test,collapse="")
  Change = FALSE
  while ((is.element(model_preTest, BanList))) {
    ind <- IndBuild(n_factors, n_variables, probability)
    model_preTest <- paste(ind$test,collapse="")
    CutChange <- CutChange + 1
    if (CutChange > ChangeRun){
      Change <- TRUE
      BanList <- NULL
    }
  }
  x <- names(Data)
  f <- rep("F",n_factors)
  f <- paste(f,c(1:n_factors),sep = "")
  mod <- NULL
  for (i in 1:n_factors) {
    mod_sub <- paste(f[i],"=~",paste(x[which(ind$ind_use[i,]==1)],collapse = "+"),sep = "")
    mod <- c(mod,mod_sub)
  }
  model_cfa <- paste(mod,collapse = "
                     ")
  return(list(ModelInd = ind$test, CFAModel = model_cfa, Change = Change))
} 
##########################################################################################
##########################################################################################
