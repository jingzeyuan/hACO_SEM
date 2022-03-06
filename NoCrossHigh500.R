getwd()
#set new path for my desktop
setwd("D:/Dropbox (UFL)/ACO_paper/hACO_rerun_Modified/hACO Codes fit change/Simulation Runs/ACOResults500/NoCrossHigh500")

library(MASS)
library(lavaan)
library(ShortForm)
library(matrixcalc)
library(lubridate)
library(svMisc)
library(stringr)


#set factor loading
loading <-
  c(0.6,
    0,
    0,
    0.5,
    0.6,
    0,
    0.7,
    0,
    0.8,
    0,
    0,
    0.7,
    0.0,
    0.9,
    0.7,
    0,
    0,
    0.5,
    0,
    0.6)
loading_m <- matrix(data = loading, ncol = 2, byrow = TRUE)
loading_m

#set factor variance_covariance matrix
fac_co <- matrix(data = c(1, 0.15, 0.15, 1),
                 nrow = 2,
                 byrow = TRUE)

pre_load <- list()
pre_unload <- list(c(7),c(3))
recored <- NULL
Fit_In <- c("srmr","rmsea","ifi","tli","cfi","mfi")
Fit_C_I_R <- c(5,3,2)
Fit_M_T_S <- c(6,4,1)

for (inter in 1:100) {
  
  simu_data <- SimuData(Loadings = loading_m, FactorCov = fac_co, N = 500)
  
  initial.start.time <- Sys.time()
  test <- hACO(Data = simu_data,n_factors = 2, n_variables = 10, loaded = pre_load, 
               unloaded = pre_unload,AccList = NULL, BanList = NULL, Phe_level = NULL, 
               standards = c(0.93,1-0.07,0.93), maxRun = 1000, Punish_rate = 0.1, Fit = Fit_C_I_R)
  
  (final.time <- difftime(Sys.time(), initial.start.time, tz,
                          units = "secs"))
  
  #BeModel <- BeyondModel(PheLevel = test$PheLevel,n_factors = 2, n_variables = 10)
  #beyondmodel <- BeyondModelEst(PheLevel = test$PheLevel, Data = simu_data,n_factors = 2, n_variables = 10, Fit = Fit)
  
  #population model
  
  mod2 <- paste("F1 =~ X1+X3+X4+X5+X8", '
                ', "F2 =~ X2+X6+X7+X9+X10")
  mod_2 <-
    paste("F1 =~ X1+X3+X4+X5+X8", '/', "F2 =~ X2+X6+X7+X9+X10")
  
  fit_cfa_2 <- try(cfa(mod2, data = simu_data), silent = TRUE)
  
  
  ####################################################################
  ######################tabu##########################################
  initial.start.time <- Sys.time()
  pro <- ProbBuild(n_factors = 2, n_variables = 10, loaded = pre_load, unloaded = pre_unload)
  mod_tabu <- ModelBuilder(Data = simu_data, n_factors = 2,n_variables = 10,probability = pro,BanList = NULL, CutChange = 0, ChangeRun = 1000)
  
  fit.start <- try(lavaan(mod_tabu$CFAModel, data=simu_data,
                          auto.var=TRUE, auto.fix.first=FALSE,
                          std.lv=TRUE,auto.cov.lv.x=TRUE), silent = TRUE)
  
  while (class(fit.start)[1] == "try-error") {
    mod_tabu <- ModelBuilder(Data = simu_data, n_factors = 2,n_variables = 10,probability = pro,BanList = NULL, CutChange = 0, ChangeRun = 1000)
    
    fit.start <- try(lavaan(mod_tabu$CFAModel, data=simu_data,
                            auto.var=TRUE, auto.fix.first=FALSE,
                            std.lv=TRUE,auto.cov.lv.x=TRUE), silent = TRUE)
  }
  
  # prepare for specification search
  spec.table<-search.prep(fit.start,
                          loadings=TRUE,fcov=TRUE,errors=FALSE)
  
  search.result <- tabu.sem(
    fit.start,
    spec.table,
    obj = BIC,
    niter = 50,
    tabu.size = 5
  )
  
  
  (final.time_tabu <- difftime(Sys.time(), initial.start.time, tz,
                               units = "secs"))
  
  search.result$best.binvec <- search.result$best.binvec[-which(search.result$best.binvec$rhs == "F2"),]
  
  for (i in (1:(nrow(search.result$best.binvec)))) {
    if (search.result$best.binvec$lhs[i] == "F1") {
      if (search.result$best.binvec$free[i] == 1) {
        mod_tabu_1 <- search.result$best.binvec$rhs[i]
        break
      }
    }
  }
  
  for (j in ((i+1):(nrow(search.result$best.binvec)))){
    if (search.result$best.binvec$lhs[j] == "F1"){
      if (search.result$best.binvec$free[j] == 1){
        mod_tabu_1 <- paste(mod_tabu_1,search.result$best.binvec$rhs[j],sep = "+")
      }
    }
  }
  mod_tabu_1 <- paste('F1 =~', mod_tabu_1)
  
  
  for (i in (1:(nrow(search.result$best.binvec)))) {
    if (search.result$best.binvec$lhs[i] == "F2") {
      if (search.result$best.binvec$free[i] == 1) {
        mod_tabu_2 <- search.result$best.binvec$rhs[i]
        break
      }
    }
  }
  
  for (j in ((i+1):(nrow(search.result$best.binvec)))){
    if (search.result$best.binvec$lhs[j] == "F2"){
      if (search.result$best.binvec$free[j] == 1){
        mod_tabu_2 <- paste(mod_tabu_2,search.result$best.binvec$rhs[j],sep = "+")
      }
    }
  }
  mod_tabu_2 <- paste('F2 =~', mod_tabu_2)
  
  mod_tabu <- paste(mod_tabu_1, '/', mod_tabu_2)
  tabu_fit <- lavaan(
    paste(mod_tabu_1, '
                 ', mod_tabu_2),
    data = simu_data,
    auto.var = TRUE,
    auto.fix.first = FALSE,
    std.lv = TRUE,
    auto.cov.lv.x = TRUE
  )
  
  ###########################################################################################
  #Compare with different Fit indices stardard values
  
  initial.start.time <- Sys.time()
  test_Compare <- hACO(Data = simu_data,n_factors = 2, n_variables = 10, loaded = pre_load, 
                       unloaded = pre_unload,AccList = NULL, BanList = NULL, Phe_level = NULL, 
                       standards = c(0.90,1-0.10,0.90), maxRun = 1000, Punish_rate = 0.1, Fit = Fit_C_I_R)
  
  (final.time.Compare <- difftime(Sys.time(), initial.start.time, tz,
                                  units = "secs"))
  
  #BeModel_Compare <- BeyondModel(PheLevel = test_Compare$PheLevel,n_factors = 2, n_variables = 10)
  #beyondmodel_Compare <- BeyondModelEst(PheLevel = test_Compare$PheLevel, Data = simu_data,n_factors = 2, n_variables = 10, Fit = Fit)
  
  ###########################################################################################
  #Fit indices II
  initial.start.time <- Sys.time()
  test_Fit2 <- hACO(Data = simu_data,n_factors = 2, n_variables = 10, loaded = pre_load, 
                    unloaded = pre_unload,AccList = NULL, BanList = NULL, Phe_level = NULL, 
                    standards = c(0.93,1-0.07,0.93), maxRun = 1000, Punish_rate = 0.1, Fit = Fit_M_T_S)
  
  (final.time.Fit2 <- difftime(Sys.time(), initial.start.time, tz,
                               units = "secs"))
  
  
  ###########################################################################################
  #Fit indices II compare
  initial.start.time <- Sys.time()
  test_Fit2_Compare <- hACO(Data = simu_data,n_factors = 2, n_variables = 10, loaded = pre_load, 
                            unloaded = pre_unload,AccList = NULL, BanList = NULL, Phe_level = NULL, 
                            standards = c(0.90,1-0.10,0.90), maxRun = 1000, Punish_rate = 0.1, Fit = Fit_M_T_S)
  
  (final.time.Fit2.Compare <- difftime(Sys.time(), initial.start.time, tz,
                                       units = "secs"))
  
  
  #recored results
  options(warn = 2)
  #Fit1 <- try(fitMeasures(test$OptModelResults$Modelres, Fit_In[Fit_C_T_S[1]]), silent = TRUE)
  #Fit2 <- try(fitMeasures(test$OptModelResults$Modelres, Fit_In[Fit_C_T_S[2]]), silent = TRUE)
  #Fit3 <- try(fitMeasures(test$OptModelResults$Modelres, Fit_In[Fit_C_T_S[3]]), silent = TRUE)
  tabu_try <- try(fitMeasures(tabu_fit, Fit_In[1]),silent = TRUE)
  if (class(tabu_try)[1]== "try-error"){
    recored <-
      rbind2 (
        recored,
        c(
          test$OptModel,
          fitMeasures(test$OptModelResults$Modelres, c("chisq", "df", Fit_In[Fit_C_I_R])),
          BIC(test$OptModelResults$Modelres),
          AIC(test$OptModelResults$Modelres),
          #
          mod_2,
          fitMeasures(fit_cfa_2, c("chisq", "df", Fit_In)),
          BIC(fit_cfa_2),
          AIC(fit_cfa_2),
          final.time,
          inter,
          #
          mod_tabu,
          c(0,0,0,0,0,0,0,0),
          0,
          0,
          final.time_tabu,
          #
          test_Compare$OptModel,
          fitMeasures(test_Compare$OptModelResults$Modelres, c("chisq", "df", Fit_In[Fit_C_I_R])),
          BIC(test_Compare$OptModelResults$Modelres),
          AIC(test_Compare$OptModelResults$Modelres),
          final.time.Compare,
          #
          test_Fit2$OptModel,
          fitMeasures(test_Fit2$OptModelResults$Modelres, c("chisq", "df", Fit_In[Fit_M_T_S])),
          BIC(test_Fit2$OptModelResults$Modelres),
          AIC(test_Fit2$OptModelResults$Modelres),
          final.time.Fit2,
          #
          test_Fit2_Compare$OptModel,
          fitMeasures(test_Fit2_Compare$OptModelResults$Modelres, c("chisq", "df", Fit_In[Fit_M_T_S])),
          BIC(test_Fit2_Compare$OptModelResults$Modelres),
          AIC(test_Fit2_Compare$OptModelResults$Modelres),
          final.time.Fit2.Compare
        )
      )  
    #assign(paste("DataRedo", inter, sep = ""), simu_data)
  }else{
    recored <-
      rbind2 (
        recored,
        c(
          test$OptModel,
          fitMeasures(test$OptModelResults$Modelres, c("chisq", "df", Fit_In[Fit_C_I_R])),
          BIC(test$OptModelResults$Modelres),
          AIC(test$OptModelResults$Modelres),
          #
          mod_2,
          fitMeasures(fit_cfa_2, c("chisq", "df", Fit_In)),
          BIC(fit_cfa_2),
          AIC(fit_cfa_2),
          final.time,
          inter,
          #
          mod_tabu,
          fitMeasures(tabu_fit, c("chisq", "df", Fit_In)),
          BIC(tabu_fit),
          AIC(tabu_fit),
          final.time_tabu,
          #
          test_Compare$OptModel,
          fitMeasures(test_Compare$OptModelResults$Modelres, c("chisq", "df", Fit_In[Fit_C_I_R])),
          BIC(test_Compare$OptModelResults$Modelres),
          AIC(test_Compare$OptModelResults$Modelres),
          final.time.Compare,
          #
          test_Fit2$OptModel,
          fitMeasures(test_Fit2$OptModelResults$Modelres, c("chisq", "df", Fit_In[Fit_M_T_S])),
          BIC(test_Fit2$OptModelResults$Modelres),
          AIC(test_Fit2$OptModelResults$Modelres),
          final.time.Fit2,
          #
          test_Fit2_Compare$OptModel,
          fitMeasures(test_Fit2_Compare$OptModelResults$Modelres, c("chisq", "df", Fit_In[Fit_M_T_S])),
          BIC(test_Fit2_Compare$OptModelResults$Modelres),
          AIC(test_Fit2_Compare$OptModelResults$Modelres),
          final.time.Fit2.Compare
        )
      ) 
    if (signif(BIC(test$OptModelResults$Modelres)) != signif(BIC(fit_cfa_2))){
      assign(paste("DataRedo", inter, sep = ""), simu_data)
    }
  }
  options(warn = 0)  
  print(paste("Iteration =",inter))
}

Names <-
  c(
    "OptModel",
    "chisq",
    "df",
    Fit_In[Fit_C_I_R],
    "BIC",
    "AIC",
    "PopModel",
    "chisq",
    "df",
    Fit_In,
    "BIC",
    "AIC",
    "RunTime",
    "Iter",
    "TabuModel",
    "chisq",
    "df",
    Fit_In,
    "BIC",
    "AIC",
    "RunTimeTabu",
    "ComModel",
    "chisq",
    "df",
    Fit_In[Fit_C_I_R],
    "BIC",
    "AIC",
    "RunTimeCom",
    "OptFit2",
    "chisq",
    "df",
    Fit_In[Fit_M_T_S],
    "BIC",
    "AIC",
    "RunTimeFit2",
    "OptFit2Comp",
    "chisq",
    "df",
    Fit_In[Fit_M_T_S],
    "BIC",
    "AIC",
    "RunTimeFit2Comp"
  )

recored <- data.frame(recored)

names(recored) <- Names

write.csv(recored,file = "D:/Dropbox (UFL)/ACO_paper/hACO_rerun_Modified/hACO Codes fit change/Simulation Runs/ACOResults500/NoCrossHigh500/NoCrossHighTabu500.csv", row.names = FALSE)

correct_rate <- c(0,0,0,0,0)
correct_rate[1] <- ACC_Rate(standard = recored$PopModel[1], results = recored$OptModel)
correct_rate[2] <- ACC_Rate(standard = recored$PopModel[1], results = recored$TabuModel)
correct_rate[3] <- ACC_Rate(standard = recored$PopModel[1], results = recored$ComModel)
correct_rate[4] <- ACC_Rate(standard = recored$PopModel[1], results = recored$OptFit2)
correct_rate[5] <- ACC_Rate(standard = recored$PopModel[1], results = recored$OptFit2Comp)

write.csv(correct_rate,file = "D:/Dropbox (UFL)/ACO_paper/hACO_rerun_Modified/hACO Codes fit change/Simulation Runs/ACOResults500/NoCrossHigh500/NH500correct_rate.csv", row.names = FALSE)

AveSearchTime <- c(0,0,0,0,0)
AveSearchTime[1] <- mean(as.numeric(recored$RunTime))
AveSearchTime[2] <- mean(as.numeric(recored$RunTimeTabu))
AveSearchTime[3] <- mean(as.numeric(recored$RunTimeCom))
AveSearchTime[4] <- mean(as.numeric(recored$RunTimeFit2))
AveSearchTime[5] <- mean(as.numeric(recored$RunTimeFit2Comp))

write.csv(AveSearchTime,file = "D:/Dropbox (UFL)/ACO_paper/hACO_rerun_Modified/hACO Codes fit change/Simulation Runs/ACOResults500/NoCrossHigh500/NH500AveSearchTime.csv", row.names = FALSE)

