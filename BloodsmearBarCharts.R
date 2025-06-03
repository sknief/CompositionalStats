# Barplots mini section

#Q1 stuff

mean_comp_bargraaphs <- function(Y, X1, X2) {
  mean_data <- data.frame(
    Day = rep(levels(X1), each = length(levels(X2))),
    Treatment = rep(levels(X2), times = length(levels(X1))))
  

  #for each loop
  for(i in 1:nrow(mean_data)) {
    mean_data$Granu[i] = mean(Y[X1 == mean_data[i,1] & X2 == mean_data[i,2]])[1]
    mean_data$Lymph[i] = mean(Y[X1 == mean_data[i,1] & X2 == mean_data[i,2]])[2]
    mean_data$Mono[i] = mean(Y[X1 == mean_data[i,1] & X2 == mean_data[i,2]])[3]
  }
  
  assign("Tester", mean_data, envir = .GlobalEnv) # assign the cleaned object to a new variable name
  return(mean_data)
}

mean_comp_bargraaphs_FBC <- function(Y, X1, X2) {
  mean_data <- data.frame(
    Day = rep(levels(X1), each = length(levels(X2))),
    Treatment = rep(levels(X2), times = length(levels(X1))))
  
  
  #for each loop
  for(i in 1:nrow(mean_data)) {
    mean_data$Granu[i] = mean(Y[X1 == mean_data[i,1] & X2 == mean_data[i,2]])[1]
    mean_data$Lymph[i] = mean(Y[X1 == mean_data[i,1] & X2 == mean_data[i,2]])[2]
    mean_data$Mono[i] = mean(Y[X1 == mean_data[i,1] & X2 == mean_data[i,2]])[3]
    mean_data$RBC[i] = mean(Y[X1 == mean_data[i,1] & X2 == mean_data[i,2]])[4]
  }
  
  assign("Tester", mean_data, envir = .GlobalEnv) # assign the cleaned object to a new variable name
  return(mean_data)
}

########

mean_comp_bargraaphs_FBC(Y, X1, X2)

Tester_Acomp <- acomp(Tester[,c("Granu", "Lymph", "Mono", 'RBC')]) #make a composition of the mean data
Tester_Acomp


#money maker
barplot(Tester_Acomp, beside = FALSE, col = c('#feedde','#fdbe85','#fd8d3c','#d94701'), 
        names.arg = c("D2 Con.", "D2 Chal.", "D7 Con.", "D7 Chal.", "D14 Con.", "D14 Chal.", "D21 Con.", "D21 Chal."), las = 2, ylim = c(0, 1), 
        main = "Mean Composition by Day and Treatment (OT FBC)", ylab = "Relative Contribution",
)



#Q2 stuff
mean_comp_bargraaphs_Q2 <- function(Y, X1, X2) {
  mean_data <- data.frame(
    Temp = rep(levels(X1), each = length(levels(X2))),
    Treatment = rep(levels(X2), times = length(levels(X1))))
  
  
  #for each loop
  for(i in 1:nrow(mean_data)) {
    mean_data$Granu[i] = mean(Y[X1 == mean_data[i,1] & X2 == mean_data[i,2]])[1]
    mean_data$Lymph[i] = mean(Y[X1 == mean_data[i,1] & X2 == mean_data[i,2]])[2]
    mean_data$Mono[i] = mean(Y[X1 == mean_data[i,1] & X2 == mean_data[i,2]])[3]
  }
  
  assign("Tester", mean_data, envir = .GlobalEnv) # assign the cleaned object to a new variable name
  return(mean_data)
}


mean_comp_bargraaphs_Q2_FBC <- function(Y, X1, X2) {
  mean_data <- data.frame(
    Temp = rep(levels(X1), each = length(levels(X2))),
    Treatment = rep(levels(X2), times = length(levels(X1))))
  
  
  #for each loop
  for(i in 1:nrow(mean_data)) {
    mean_data$Granu[i] = mean(Y[X1 == mean_data[i,1] & X2 == mean_data[i,2]])[1]
    mean_data$Lymph[i] = mean(Y[X1 == mean_data[i,1] & X2 == mean_data[i,2]])[2]
    mean_data$Mono[i] = mean(Y[X1 == mean_data[i,1] & X2 == mean_data[i,2]])[3]
    mean_data$RBC[i] = mean(Y[X1 == mean_data[i,1] & X2 == mean_data[i,2]])[4]
  }
  
  assign("Tester", mean_data, envir = .GlobalEnv) # assign the cleaned object to a new variable name
  return(mean_data)
}

########

mean_comp_bargraaphs_Q2_FBC(Y, X1, X2)

Tester_Acomp <- acomp(Tester[,c("Granu", "Lymph", "Mono", 'RBC')]) #make a composition of the mean data
Tester_Acomp


#money maker
barplot(Tester_Acomp, beside = FALSE, col = c('#feedde','#fdbe85','#fd8d3c','#d94701'), 
        names.arg = c("CC Con.", "CC Chal.", "FC Con.", "FC Chal.", "OT Con.", "OT Chal."), las = 2, ylim = c(0, 1), 
        main = "Mean Composition by Temp and Treatment (Day 21 FBC)", ylab = "Relative Contribution",
)

#other colour scheme: '#fee0d2','#fc9272','#de2d26'