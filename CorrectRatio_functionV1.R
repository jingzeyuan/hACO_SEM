#get the percentage of correct results
ACC_Rate <- function(standard, results){
  total <- length(results)
  #maker standard results to compare
  stan <- as.vector(as.data.frame(strsplit(standard,"F"))[,1])
  stan <- as.character(tidyr::extract_numeric(stan))
  number <- 0
  for (i in c(1:length(stan))) {
    stan_tem <- as.vector(as.data.frame(stringr::str_split(stan[i],"")[[1]]))[,1]
    if (length(which(!is.na(stan_tem))) > 1){
      number <- number + 1
      assign(paste("stan",number,sep = ""), stan_tem[2:length(stan_tem)])
    }
  }
  #compare each model with standards
  num_correct <- 0
  for (i in c(1:length(results))) {
    #get the information of the one results model
    res_tem <- results[i]
    res_tem <- as.vector(as.data.frame(strsplit(res_tem,"F"))[,1])
    res_tem <- as.character(tidyr::extract_numeric(res_tem))
    number <- 0
    for (j in c(1:length(res_tem))) {
      tem <- as.vector(as.data.frame(stringr::str_split(res_tem[j],"")[[1]]))[,1]
      if (length(which(!is.na(tem))) > 1){
        number <- number + 1
        assign(paste("tem",number,sep = ""), tem[2:length(tem)])
      }
    }
    #compare this model with standards
    correct_count <- 0
    for (j in c(1:number)){
      tem <- eval(parse(text =paste("tem",j,sep = "")))
      for (k in c(1:number)) {
        stan_tem <- eval(parse(text =paste("stan",k,sep = "")))
        if ((length(tem) == length(stan_tem))&(length(which(tem %in% stan_tem))==length(tem))){
          correct_count <- correct_count + 1
        }
      }
    }
    # if all factors are matched
    if (correct_count == number){
      num_correct <- num_correct + 1
    }
  }
  ratio <- num_correct/total
  return(ratio)
}




