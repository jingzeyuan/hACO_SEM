library(MASS)
library(lavaan)
library(ShortForm)
library(matrixcalc)
library(lubridate)
library(svMisc)
library(stringr)
#################################################################
#################################################################
#####get the the optimum model###################################
OptModel <- function(Data, PheLevel, n_factors, n_variables, rate) {
  PheLevel <- signif(PheLevel, 3)
  PheLevel[which(PheLevel <= rate * max(PheLevel))] <- 0
  PheLevel[which(PheLevel != 0)] <- 1
  Ind_OptModel <-
    matrix(PheLevel, n_factors, n_variables, byrow = TRUE)
  x <- names(Data)
  f <- rep("F", n_factors)
  f <- paste(f, c(1:n_factors), sep = "")
  Optmod <- NULL
  for (i in 1:n_factors) {
    Optmod_sub <-
      paste(f[i], "=~", paste(x[which(Ind_OptModel[i, ] == 1)], collapse = "+"), sep = "")
    Optmod <- c(Optmod, Optmod_sub)
  }
  OptModelRun <- paste(Optmod, collapse = "
                     ")
  OptModelShow <- paste(Optmod, collapse = "   ")
  return(list(CFAModel = OptModelRun, OptModelShow = OptModelShow, ModelInd = PheLevel, ModelPreTest = Ind_OptModel))
}


##############################################################
##############################################################
###########optimal model estimation###########################
OptModelEst <- function(Model, Data, Fit, dicrimination = dicrimination, standards = standards, 
                        scopes= scopes, n_factors, n_variables, Punish_rate) {
  options(warn = 2)
  fit_cfa <-
    try(lavaan(
      Model$CFAModel,
      data = Data,
      auto.var = TRUE,
      auto.fix.first = FALSE,
      std.lv = TRUE,
      auto.cov.lv.x = TRUE
    ),
    silent = TRUE)
  Fit_In <- c("srmr","rmsea","ifi","tli","cfi","mfi")
  Fit1 <- try(fitMeasures(fit_cfa, Fit_In[Fit[1]]), silent = TRUE)
  Fit2 <- try(fitMeasures(fit_cfa, Fit_In[Fit[2]]), silent = TRUE)
  Fit3 <- try(fitMeasures(fit_cfa, Fit_In[Fit[3]]), silent = TRUE)
  bic <- try(BIC(fit_cfa), silent = TRUE)
  Fits <- cbind(Fit1, Fit2, Fit3)  
  
  if (class(Fit1)[1] == "try-error"|class(Fit2)[1] == "try-error"|class(Fit3)[1] == "try-error"|class(bic)[1] == "try-error"){
    
    if (length(which(Fit < 3)) > 0){
      Fits[which(Fit < 3)] <- 1
      Fits[which(Fit >= 3)] <- 0
      Fits <- as.numeric(Fits)
      bic <- Inf
      ph <- 0
    }else{
      Fits <- 0
      bic <- -Inf
      ph <- 0
    }
  }else{
    fit_temp <- Fits
    
    if (length(which(fit_temp > 1)) > 0) {
      fit_temp[which(fit_temp > 1)] <- 1
    }
    
    if (length(which(Fit < 3)) > 0) {
      fit_temp[which(Fit < 3)] <- 1 - fit_temp[which(Fit < 3)]
    }
    
    if (fit_temp[1] < 0.80 | fit_temp[2] < 0.80 | fit_temp[3] < 0.80) {
      ph <- 0
    } else{
      ph <- 1
      for (i in 1:(length(fit_temp))) {
        ph_sub <- (exp((
          dicrimination[i] *
            (as.numeric(fit_temp[i]) - standards[i])
        ) / scopes[i])) /
          (1 + exp((
            dicrimination[i] * (as.numeric(fit_temp[i]) - standards[i])
          ) / scopes[i]))
        ph <- ph * ph_sub
      }
      
      if (sum(Model$ModelInd) < n_variables) {
        punish <- 0
      } else{
        punish <-
          Punish_rate * (log(1 + (sum(Model$ModelInd) - n_variables) / ((n_factors - 1) *
                                                                          n_variables)))
      }
      
      ph <- ph - punish
      if (ph < 0) {
        ph <- 0
      }
    }
    
    if (length(which(Fits > 1)) > 0) {
      Fits[which(Fits > 1)] <- 1
    }
    
    if (length(which(Fit < 3)) > 0) {
      Fits[which(Fit < 3)] <- 1 - Fits[which(Fit < 3)]
    }
    
  }
  
  options(warn = 0)
  return(list(
    Fit1 = Fits[1], Fit2 = Fits[2], Fit3 = Fits[3],
    Modelres = fit_cfa,
    Ph = ph
  ))
}


############################################################################
############################################################################
######################Function for one iteration############################

hACO <-
  function(Data,
           n_factors,
           n_variables,
           BanList = NULL,
           AccList = NULL,
           Phe_level = NULL,
           loaded,
           unloaded,
           dicrimination = c(1, 1, 1),
           standards = c(0.93, 1-0.07, 0.93),
           scopes = c(0.015, 0.015, 0.015),
           P_dicrimination = 0.5,
           P_standards = 0.95,
           P_scopes = 0.005,
           maxRun = 15000,
           LowCut = 0.001,
           HighCut = 0.75,
           ChangeRun = 1000,
           Punish_rate = 0.1,
           Ban_neibor_length = 3,
           Fit = c(5,2,3)) {
    
    
    #Calculate the Probility Matrix
    #Which incorperate the prior knowledge 
    probability <-
      ProbBuild(
        n_factors = n_factors,
        n_variables = n_variables,
        loaded = loaded,
        unloaded = unloaded
      )
    
    entropy_pre <- 100
    diff <- 100
    numRuns <- 1
    numAcc <- 0
    numTotal <- 0
    numRuns_pre <- 0
    diffRuns <- 0
    MaxPheromone <- -Inf
    Ban_neibor <- NULL
    
    
    Data <- as.data.frame(Data)
    
    if (is.null(names(Data))) {
      x <- rep("X", n_variables)
      x <- paste(x, c(1:n_variables), sep = "")
      names(Data) <- x
    }
    
    #Check the input error
    if (HighCut < LowCut) {
      stop("HighCut can not be lower than LowCut", call. = FALSE)
    }
    EntropyStart <- FALSE
    
    
    model <-
      ModelBuilder(
        Data = Data,
        n_factors = n_factors,
        n_variables = n_variables,
        probability = probability,
        BanList = BanList,
        ChangeRun = ChangeRun
      )
    
    
    
    if (model$Change == TRUE) {
      BanList <- NULL
      if (HighCut > LowCut) {
        LowCut <- LowCut
        HighCut <- HighCut - 0.10
      } else{
        return(print("Check number of factors, loaded, or unloaded"))
      }
    }
    
    
    res <-
      ModelRun(
        Model = model,
        Data = Data,
        dicrimination = dicrimination,
        standards = standards,
        scopes = scopes,
        P_dicrimination = P_dicrimination,
        P_standards = P_standards,
        P_scopes = P_scopes,
        AccList = AccList,
        BanList = BanList,
        Phe_level = Phe_level,
        probability = probability,
        n_factors = n_factors,
        n_variables = n_variables,
        numRuns = numRuns,
        LowCut = LowCut,
        HighCut = HighCut,
        Punish_rate = Punish_rate,
        Fit = Fit
      )
    
    cat("Preparetion, Patience from Old Drives!\n")
    
    Start <- model$ModelInd
    BestSearch <- res$ModelInd
    
    patience_limit <- 0
    patience_Nolimit <- 10 * n_factors*n_variables
    
    while (res$Pheromone == 0) {
      model <-
        ModelBuilder(
          Data = Data,
          n_factors = n_factors,
          n_variables = n_variables,
          probability = probability,
          BanList = BanList,
          ChangeRun = ChangeRun
        )
      

      if (model$Change == TRUE) {
        BanList <- NULL
        if (HighCut > LowCut) {
          LowCut <- LowCut
          HighCut <- HighCut - 0.10
        } else{
          return(print("Check n_factor, loaded, or unloaded"))
        }
      }
      
      res <-
        ModelRun(
          Model = model,
          Data = Data,
          dicrimination = dicrimination,
          standards = standards,
          scopes = scopes,
          P_dicrimination = P_dicrimination,
          P_standards = P_standards,
          P_scopes = P_scopes,
          AccList = AccList,
          BanList = BanList,
          Phe_level = Phe_level,
          probability = probability,
          n_factors = n_factors,
          n_variables = n_variables,
          numRuns = numRuns,
          LowCut = LowCut,
          HighCut = HighCut,
          Punish_rate = Punish_rate,
          Fit = Fit
        )
      
      if (res$Pheromone > MaxPheromone) {
        BestSearch <- res$ModelInd
        MaxPheromone <- res$Pheromone
      }
      Phe_level <- res$PheLevel
      BanList <- res$BanList
      AccList <- res$AccList
      probability <- res$Prob
      Start <- model$ModelInd
      BestSearch <- res$ModelInd
      patience_limit <- patience_limit + 1
      if (patience_limit >= patience_Nolimit){
        break
      }
    }
    neiber_res <- res 
    Good_neibor <- NULL
    
    while ((numRuns < 3000) &
           (numAcc < 200) &
           ((abs(diff) > 0.0001) | (diffRuns < 1000))  & (numTotal < maxRun)) {
        # print(paste('
        # ', "numRuns:",numRuns, ", numAcc:", numAcc,"; diffRuns:", diffRuns, "; numTotal:", numTotal, '
        # '))
        Good_ph <- NULL
        Good_ind <- NULL
        Good_bic <- NULL
      for (i in 1:length(Start)) {
        iter <- Start
       
      
        
        if (is.element(i, which(t(probability) == 1 |
                                t(probability) == 0))) {
          next
        } else{
          
          #get the neighbor of the starting model
          iter[i] <- ifelse(Start[i] == 1, 0, 1)
          #iter[which(probability == 1)] <- 1
          #iter[which(probability == 0)] <- 0
          
          modelind <- paste(iter, collapse = "")
          
          ind <- matrix(iter, n_factors, n_variables, byrow = TRUE)
          
          
          if ((!is.element(0, colSums(ind))) &
              (!is.element(paste(iter, collapse = ""), BanList)) &
              (!is.element(modelind,Ban_neibor))) {
            
            #record the good neighbers
            
            model <-
              ModelReBuild(
                ModInd = paste(iter, collapse = ""),
                n_factors = n_factors,
                n_variables = n_variables,
                Data = Data
              )
            
            neiber_res <-
              ModelRun(
                Model = model,
                Data = Data,
                dicrimination = dicrimination,
                standards = standards,
                scopes = scopes,
                P_dicrimination = P_dicrimination,
                P_standards = P_standards,
                P_scopes = P_scopes,
                AccList = AccList,
                BanList = BanList,
                Phe_level = Phe_level,
                probability = probability,
                n_factors = n_factors,
                n_variables = n_variables,
                numRuns = numRuns,
                LowCut = LowCut,
                HighCut = HighCut,
                Punish_rate = Punish_rate,
                Fit = Fit
              )
            
            if (neiber_res$Pheromone >= MaxPheromone) {
              BestSearch <- neiber_res$ModelInd
              MaxPheromone <- neiber_res$Pheromone
            }
            
            Good_ph <- rbind(Good_ph,neiber_res$Pheromone)
            Good_ind <- rbind(Good_ind, neiber_res$ModelInd)
            Good_bic <- rbind(Good_bic, neiber_res$BIC)
            
            
            #Update parameters
            Phe_level <- neiber_res$PheLevel
            BanList <- neiber_res$BanList
            AccList <- neiber_res$AccList
            probability <- neiber_res$Prob
            entropy <- neiber_res$Entropy
            diff <- entropy_pre - neiber_res$Entropy
            
            #check stop criteria 
            if (((diff != 0) & (numTotal > 0))) {
              EntropyStart <- TRUE
            }
            if ((abs(diff) < 0.0001) & (EntropyStart == TRUE)) {
              diffRuns <- diffRuns + (neiber_res$NumRuns - numRuns_pre)
            }
            entropy_pre <- neiber_res$Entropy
            progress(numTotal, maxRun)
            if (numTotal == maxRun) {
              cat("Done!\n")
            }
          } else{
            next
          }
        }      
      }
        
        
        if (!is.null(Good_ph)) {
          Good_neibor <-
            data.frame(Good_ph, Good_ind, Good_bic, stringsAsFactors = FALSE)
          
          if (max(Good_neibor$Good_ph) == 0) {
            neiber_max <-
              which(Good_neibor$Good_bic == min(Good_neibor$Good_bic))
            if (length(neiber_max) > 1) {
              neiber_choose <- sample(neiber_max, 1)
            } else{
              neiber_choose <- neiber_max
            }
          } else{
            neiber_max <-
              which(Good_neibor$Good_ph == max(Good_neibor$Good_ph))
            if (length(neiber_max) > 1) {
              neiber_choose <- sample(neiber_max, 1)
            } else{
              neiber_choose <- neiber_max
            }
          }

        
        # if (length(neiber_max) > 1){
        #   
        #   # if (length(neiber_max) > 2){
        #   #   neiber_start <- sample(neiber_max[-neiber_choose],1)
        #   # }else{
        #   #   neiber_start <- neiber_max[-neiber_choose]
        #   # }
        # }else{
        #   # good_neiber_secend <- Good_neibor[-(neiber_max),]
        #   # neiber_second <- which(good_neiber_secend$Good_ph == max(good_neiber_secend$Good_ph))
        #   # if (length(neiber_second) > 1){
        #   #   neiber_start <- sample(neiber_second,1)
        #   # }else{
        #   #   neiber_start <- neiber_second
        #   # }
        # }
        
        neiber_start <- neiber_choose
        
        
        if (length(Ban_neibor) > Ban_neibor_length) {
          for (i in 1:(length(Ban_neibor) - 1)) {
            Ban_neibor[i] <- Ban_neibor[i + 1]
          }
          Ban_neibor[Ban_neibor_length] <- as.character(Good_neibor$Good_ind[neiber_choose])
        } else{
          Ban_neibor <- c(Ban_neibor,as.character(Good_neibor$Good_ind[neiber_choose]))
        }
        
        Start <- StartBuild(Good_neibor$Good_ind[neiber_start], n_factors = n_factors, n_variables = n_variables)

      }else{
        modelbest <- ModelReBuild(BestSearch,n_factors,n_variables, Data = Data)
        Start <- modelbest$ModelInd
      }
        
      

      numRuns <- neiber_res$NumRuns
      numRuns_pre <- numRuns
      numAcc <- length(neiber_res$AccList)
      numTotal <- numTotal + 1
      progress(numTotal, maxRun)
      
      
      if (numTotal == maxRun) {
        cat("Done!\n")
      }
    }
    
    # S1 <- 0
    # S2 <- 0
    # S3 <- 0
    # Pre_sum <- 0
    # rate_opt <- 0.8
    # for (i in seq(0.98, 0.8, by = -0.01)) {
    #   Optmodel <-
    #     OptModel(
    #       Data = Data,
    #       PheLevel = Phe_level,
    #       n_factors = n_factors,
    #       n_variables = n_variables,
    #       rate = i
    #     )
    #   OptEst <- OptModelEst(Model = Optmodel, Data = Data, Fit = Fit, dicrimination = dicrimination, standards = standards, 
    #                         scopes= scopes, n_factors, n_variables, Punish_rate)
    #   S1_temp <- OptEst$Fit1
    #   S2_temp <- OptEst$Fit2
    #   S3_temp <- OptEst$Fit3
    #   if (class(OptEst$Modelres)[1] == "try-error"){
    #     #print("1")
    #     S1 <- S1_temp
    #     S2 <- S2_temp
    #     S3 <- S3_temp
    #     Pre_sum <- sum(Optmodel$ModelPreTest)
    #     next
    #   }else{
    #     if (length(which(colSums(Optmodel$ModelPreTest) == 0)) == 0) {
    #       if (((S1_temp - S1) < 0.0) |
    #           ((S2_temp - S2) < 0.0) | ((S3_temp - S3) < 0.0)) {
    #         rate_opt <- i + 0.01
    #         #print("2")
    #         break
    #       } else{
    #         if (((S1_temp - S1) < 0.02) &
    #             ((S2_temp - S2) < 0.02) &
    #             ((S3_temp - S3) < 0.02) & (Pre_sum != sum(Optmodel$ModelPreTest))) {
    #           rate_opt <- i + 0.01
    #           #print("3")
    #           break
    #         }
    #       }
    #       #print("4")
    #       Pre_sum <- sum(Optmodel$ModelPreTest)
    #     }
    #   }
    #   #print("5")
    #   S1 <- S1_temp
    #   S2 <- S2_temp
    #   S3 <- S3_temp
    #   Pre_sum <- sum(Optmodel$ModelPreTest)
    # }
    # 
    
    rate_total <- NULL
    model_sum <- NULL
    opt_fit1 <- NULL
    opt_fit2 <- NULL
    opt_fit3 <- NULL
    
    for (i in seq(0.98,0.80, by = -0.01)){
      Optmodel <-
        OptModel(
          Data = Data,
          PheLevel = Phe_level,
          n_factors = n_factors,
          n_variables = n_variables,
          rate = i
        )
      OptEst <- OptModelEst(Model = Optmodel, Data = Data, Fit = Fit, dicrimination = dicrimination, standards = standards, 
                            scopes= scopes, n_factors, n_variables, Punish_rate)
      
      if (class(OptEst$Modelres)[1] == "try-error"){
        #print("1")
        next
      }else{
        if (length(which(colSums(Optmodel$ModelPreTest) == 0)) == 0) {
          rate_total <- c(rate_total,i)
          model_sum <- c(model_sum, sum(Optmodel$ModelPreTest))
          opt_fit1 <- c(opt_fit1, OptEst$Fit1)
          opt_fit2 <- c(opt_fit2, OptEst$Fit2)
          opt_fit3 <- c(opt_fit3, OptEst$Fit3)
        }
      }
    }
    
    if (is.null(rate_total)) {
      rate_opt <- 0.8
    } else{
      opt_sel <-
        as.data.frame(cbind(rate_total, model_sum, opt_fit1, opt_fit2, opt_fit3))
      opt_sel <- opt_sel[-which(duplicated(opt_sel[, 2])), ]
      rate_opt <- 0
      for (i in nrow(opt_sel):1) {
        if (i == 1) {
          rate_opt <- opt_sel[1, 1]
          break
        } else{
          if (((opt_sel[i, 3] - opt_sel[i - 1, 3]) < 0.02) |
              ((opt_sel[i, 4] - opt_sel[i - 1, 4]) < 0.02) |
              ((opt_sel[i, 5] - opt_sel[i - 1, 5]) < 0.02)) {
            next
          } else{
            rate_opt <- opt_sel[i, 1]
            break
          }
        }
      }
    }
    
    
 
    Optmodel <-
      OptModel(
        Data = Data,
        PheLevel = Phe_level,
        n_factors = n_factors,
        n_variables = n_variables,
        rate = rate_opt
      )
    OptEst <- OptModelEst(Model = Optmodel, Data = Data, Fit = Fit, dicrimination = dicrimination, standards = standards, 
                          scopes= scopes, n_factors, n_variables, Punish_rate)
    
    
    BestModel <-
      ModelReBuild(ModInd = BestSearch,
                   n_factors = n_factors,
                   n_variables = n_variables,
                   Data = Data)
    BestEst <- OptModelEst(Model = BestModel, Data = Data, Fit = Fit, dicrimination = dicrimination, standards = standards, 
                           scopes= scopes, n_factors, n_variables, Punish_rate)
    print("Iteration Done")
    return(
      list(
        OptModel = Optmodel$OptModelShow,
        OptModelResults = OptEst,
        OptModelRun = Optmodel$CFAModel,
        PheLevel = Phe_level,
        Prob = probability,
        BanList = BanList,
        AccList = unique(AccList),
        BestSearchModel = BestModel$CFAModel,
        BestEst = BestEst
      )
    )
  }
