###############################################
# Day by Temp (Q1) Analysis Code             ##
## SMSK 2025                                 ##
###############################################

#### First: Packages and Prereqs #####
library(compositions)
library(rgl)
library(tensorA)
library(colorspace)
library(robustbase)
library(energy)
library(cramer)
library(combinat)
library(MCMCpack)
library(tidyr)
library(MASS)
library(plyr)
library(dplyr)
library(car)
library(emmeans)

#### Define Functions ####

# this one explores the raw acomp and its three main transformations
comp_explore <- function(acomp_object, identifier) {
  comp_Ident <- paste0(identifier, "_comp_explore_graphs_") #this is to identify the output from your         specific run or hypothesis
  pdf(file = paste(comp_Ident, ".pdf", sep = "_"))
  #actual plot code
  plot(acomp_object)
  title("Untransformed Composition",outer=TRUE,line=-1)
  plot(clr(acomp_object)) #clr transformation (central log)
  title("CLR-transformed Composition",outer=TRUE,line=-1)
  plot(ilr(acomp_object)) #ilr transformation (isometric log) (check)
  title("ILR-transformed Composition",outer=TRUE,line=-1)
  dev.off()
  print("Graphs Made!")
}

# this one deals with any missing values
comp_clean <- function(acomp_object, identifier){
  BDL_MNAR_REPLACE <<- 0.001 #assign a value to replace BDL and MNARs
  acomp_object_replaced <<- zeroreplace(acomp_object, BDL_MNAR_REPLACE) #replace zeros
  acomp_object_replaced[is.na(acomp_object_replaced)] <<- BDL_MNAR_REPLACE #replace NAs
  
  original_name <<- deparse1(substitute(acomp_object)) # get the original variable name
  new_name <<- paste0(original_name, "_cleaned")
  assign(new_name, acomp_object_replaced, envir = .GlobalEnv) # assign the cleaned object to a new variable name
  
  clean_Ident <- paste0(identifier, "_cleaned_comp_explore_graphs_") 
  pdf(file = paste(clean_Ident, ".pdf", sep = "_"))
  #PDF code
  plot(0:10, type = "n", xaxt="n", yaxt="n", bty="n", xlab = "", ylab = "")
  text(5, 8, "Composition after treating missing values")
  text(5, 7, identifier)
  
  plot(acomp_object_replaced)
  title("Untransformed Cleaned Composition",outer=TRUE,line=-1)
  plot(clr(acomp_object_replaced)) #clr transformation (central log)
  title("CLR-transformed Cleaned Composition",outer=TRUE,line=-1)
  plot(ilr(acomp_object_replaced)) #ilr transformation (isometric log) (check)
  title("ILR-transformed Cleaned Composition",outer=TRUE,line=-1)
  
  dev.off()
  print("Graphs Made!")
}


#better panel making code from stackoverflow
panel.qq <- function(x, y, ...) {
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1), new = TRUE)
  qqplot(x, y, xlab = deparse(substitute(x)), ylab = deparse(substitute(y)))
  abline(c(0,1), ...)
}

#this one centers the acomp for you and makes some nice explorative graphs
comp_center <- function(acomp_object, identifier) {
  
  acomp_object_centered <- acomp_object - mean(acomp_object) #centering by subtracting the mean
  print(c("Centered mean is ", mean(acomp_object_centered))) #check if the centered mean is the neutral element
  
  original_name <<- deparse1(substitute(acomp_object)) # get the original variable name
  new_name <<- paste0(original_name, "_centered")
  assign(new_name, acomp_object_centered, envir = .GlobalEnv) # assign the cleaned object to a new variable name
  
  #PDF export code for this section: 
  center_Ident <- paste0(identifier, "_centered_comp_explore_graphs_") 
  pdf(file = paste(center_Ident, ".pdf", sep = "_"))
  plot(0:10, type = "n", xaxt="n", yaxt="n", bty="n", xlab = "", ylab = "")
  text(5, 8, "Centering the Data (Normality Fix)")
  text(5, 7, identifier)
  
  plot(acomp_object)
  title("Current Aitch. Composition",outer=TRUE,line=-1)
  pairs(acomp_object, lower.panel = panel.qq) #QQplots and histos
  title("Current Aitch. Composition",outer=TRUE,line=-1)
  plot(acomp_object_centered)
  title("Centered Aitch. Composition",outer=TRUE,line=-1)
  pairs(acomp_object_centered, lower.panel = panel.qq)
  title("Centered Aitch. Composition",outer=TRUE,line=-1)
  
  dev.off()
  
}

# bulk make anova plots if your focal group has two levels
comp_ANOVA_plots_2 <- function(model, Y, Group, identifier) {
  #plotting the variances of estimated parameters as ellipses
  vars <- vcovAcomp(model)
  dim(vars)
  alpha=0.05
  a = ilrInv(coef(model)[1,],orig=Y) #intercept
  b = ilrInv(rbind(0,coef(model)[-1,]),orig=Y) #slope
  mu = a + b #mean
  Sigma=ilrvar2clr(var(model)) #variance
  
  anova_Ident <- paste0(identifier, "_anova_graphs_") 
  pdf(file = paste(anova_Ident, ".pdf", sep = "_"))
  plot(0:10, type = "n", xaxt="n", yaxt="n", bty="n", xlab = "", ylab = "")
  text(5, 8, "ANOVA model plots")
  text(5, 7, identifier)
  
  #plotting the means and the variance ellipses based on the model data
  plot(mu,pch=20, col= c("blue","green"))
  plot(Y,pch=".",add=TRUE)
  ellipses(mu[1,],Sigma,2, col = "blue")
  ellipses(mu[2,],Sigma,2, col = "green")
  legend(x=0.75,y=0.65,abbreviate(levels(Group), minlength=1), pch=20,col=c("blue","green"),yjust=0)
  title("Means and Variance Ellipses",outer=TRUE,line=-1)
  
  coefs <- ilrInv(coef(model),orig=Y)
  
  #plotting the variance
  plot(coefs,pch=as.character(1:nrow(coefs)))
  plot(coefs-coefs,pch=20,add=TRUE)
  for(i in 1:nrow(coefs)){
    ellipses(acomp(coefs[i,,drop=FALSE]),
             ilrvar2clr(vars[,,i,i]),
             r=ConfRadius(model,1-alpha),
             col="gray")
  }
  par(xpd=FALSE)
  
  dev.off()  
  print("Graphs Made!")
}


# bulk make anova plots if your focal group has THREE levels
comp_ANOVA_plots_3 <- function(model, Y, Group, identifier) {
  #plotting the variances of estimated parameters as ellipses
  vars <- vcovAcomp(model)
  dim(vars)
  alpha=0.05
  a = ilrInv(coef(model)[1,],orig=Y) #intercept
  b = ilrInv(rbind(0,coef(model)[-1,]),orig=Y) #slope
  mu = a + b #mean
  Sigma=ilrvar2clr(var(model)) #variance
  
  anova_Ident <- paste0(identifier, "_anova_graphs_") 
  pdf(file = paste(anova_Ident, ".pdf", sep = "_"))
  plot(0:10, type = "n", xaxt="n", yaxt="n", bty="n", xlab = "", ylab = "")
  text(5, 8, "ANOVA model plots")
  text(5, 7, identifier)
  
  #plotting the means and the variance ellipses based on the model data
  plot(mu,pch=20, col= c("blue","green", "red"))
  plot(Y,pch=".",add=TRUE)
  ellipses(mu[1,],Sigma,2, col = "blue")
  ellipses(mu[2,],Sigma,2, col = "green")
  ellipses(mu[3,],Sigma,2, col = "red")
  legend(x=0.75,y=0.65,abbreviate(levels(Group), minlength=1), pch=20,col=c("blue","green", "red"),yjust=0)
  title("Means and Variance Ellipses",outer=TRUE,line=-1)
  
  coefs <- ilrInv(coef(model),orig=Y)
  
  #plotting the variance
  plot(coefs,pch=as.character(1:nrow(coefs)))
  plot(coefs-coefs,pch=20,add=TRUE)
  for(i in 1:nrow(coefs)){
    ellipses(acomp(coefs[i,,drop=FALSE]),
             ilrvar2clr(vars[,,i,i]),
             r=ConfRadius(model,1-alpha),
             col="gray")
  }
  par(xpd=FALSE)
  
  dev.off()  
}


# bulk make anova plots if your focal group has FOUR levels
comp_ANOVA_plots_4 <- function(model, Y, Group, identifier) {
  #plotting the variances of estimated parameters as ellipses
  vars <- vcovAcomp(model)
  dim(vars)
  alpha=0.05
  a = ilrInv(coef(model)[1,],orig=Y) #intercept
  b = ilrInv(rbind(0,coef(model)[-1,]),orig=Y) #slope
  mu = a + b #mean
  Sigma=ilrvar2clr(var(model)) #variance
  
  anova_Ident <- paste0(identifier, "_anova_graphs_") 
  pdf(file = paste(anova_Ident, ".pdf", sep = "_"))
  plot(0:10, type = "n", xaxt="n", yaxt="n", bty="n", xlab = "", ylab = "")
  text(5, 8, "ANOVA model plots")
  text(5, 7, identifier)
  
  #plotting the means and the variance ellipses based on the model data
  plot(mu,pch=20, col= c("blue","green", "red", "black"))
  plot(Y,pch=".",add=TRUE)
  ellipses(mu[1,],Sigma,2, col = "blue")
  ellipses(mu[2,],Sigma,2, col = "green")
  ellipses(mu[3,],Sigma,2, col = "red")
  ellipses(mu[4,],Sigma,2, col = "black")
  legend(x=0.75,y=0.65,abbreviate(levels(Group), minlength=1), pch=20,col=c("blue","green", "red", "black"),yjust=0)
  title("Means and Variance Ellipses",outer=TRUE,line=-1)
  
  coefs <- ilrInv(coef(model),orig=Y)
  
  #plotting the variance
  plot(coefs,pch=as.character(1:nrow(coefs)))
  plot(coefs-coefs,pch=20,add=TRUE)
  for(i in 1:nrow(coefs)){
    ellipses(acomp(coefs[i,,drop=FALSE]),
             ilrvar2clr(vars[,,i,i]),
             r=ConfRadius(model,1-alpha),
             col="gray")
  }
  par(xpd=FALSE)
  
  dev.off()  
}

comp_ANOVA_diag <- function(model, Y_variable, Group, identifier) {
  diag_Ident <- paste0(identifier, "_anova_diagnostics_graphs_") 
  pdf(file = paste(diag_Ident, ".pdf", sep = "_"))
  plot(0:10, type = "n", xaxt="n", yaxt="n", bty="n", xlab = "", ylab = "")
  text(5, 8, "ANOVA diagnostics plots")
  text(5, 7, identifier)
  
  Resid = ilrInv(resid(model),orig=Y_variable)
  plot(Resid, cex=0.5)
  title("Ternary Diagrams of Residuals",outer=TRUE,line=-1)
  
  boxplot(Resid) #also reveals outliers
  title("Boxplots of Residuals",outer=TRUE,line=-1)
  
  #QQ plot (also normally distributed residuals)
  qqnorm(Resid) 
  
  #homoschedascity (equal variance)
  Pred = ilrInv(predict(model),orig=Y_variable)
  opar <- par(oma=c(3,3,0,0),mar=c(4,4,1,0))
  pairwisePlot(clr(Pred),clr(Resid))
  mtext(text=c("predicted values (clr)","residuals (clr)"),
        side=c(1,2),at=0.5,line=2,outer=TRUE)
  par(opar)
  
  #color-coordinate the above diagram by group
  opar <- par(oma=c(0,0,2,0),mar=c(4,4,1,0))
  pairwisePlot(clr(Pred),clr(Resid),col=as.numeric(X1))
  legend(x=0.85,y=0.65,levels(Group),col=as.numeric(1:3),pch=1,xpd=NA) #starts a click to place mechanism for the legend
  par(opar)
  
  dev.off()
  print("Graphs Made!")
}




#### After packages and functions, add and structure your data ####
DATA <- as.data.frame(Bloodsmear_data_full)
#$ [Optional: Data Cleaning]
DATA$Timepoint <- replace(DATA$Timepoint, DATA$Timepoint == 17, 7) #okay that worked

#Subsetting
DATA$Timepoint <- as.factor(DATA$Timepoint)
levels(DATA$Timepoint)

#Splitting by temp
DATA$TempTreat <- as.factor(DATA$TempTreat)
levels(DATA$TempTreat)

DATA <- DATA[1:13] 

Cold_Data = subset(DATA, DATA$TempTreat == "Cold")
Cold_x100s = Cold_Data[c(1:6, 11:13)]
Cold_x1000s = Cold_Data[1:10]

Fluc_Data = subset(DATA, DATA$TempTreat == "Fluc")
Fluc_x100s = Fluc_Data[c(1:6, 11:13)]
Fluc_x1000s = Fluc_Data[1:10]

Opt_Data = subset(DATA, DATA$TempTreat == "Opt")
Opt_x100s = Opt_Data[c(1:6, 11:13)]
Opt_x1000s = Opt_Data[1:10]


# make a composition
Cold_x100s_acomp = acomp(Cold_x100s[7:9])
Cold_x100s_acomp
Cold_x1000s_acomp = acomp(Cold_x1000s[7:10])
Cold_x1000s_acomp

Fluc_x100s_acomp = acomp(Fluc_x100s[7:9])
Fluc_x100s_acomp
Fluc_x1000s_acomp = acomp(Fluc_x1000s[7:10])
Fluc_x1000s_acomp

Opt_x100s_acomp = acomp(Opt_x100s[7:9])
Opt_x100s_acomp
Opt_x1000s_acomp = acomp(Opt_x1000s[7:10])
Opt_x1000s_acomp


#### Then, use the custom functions to explore, clean and center your comps #####

#x100s Explore
comp_explore(Cold_x100s_acomp, "Cold_x100s")
comp_explore(Fluc_x100s_acomp, "Fluc_x100s")
comp_explore(Opt_x100s_acomp, "Opt_x100s")

#x1000s Explore
comp_explore(Cold_x1000s_acomp, "Cold_x1000s")
comp_explore(Fluc_x1000s_acomp, "Fluc_x1000s")
comp_explore(Opt_x1000s_acomp, "Opt_x1000s")


#x100s Clean
comp_clean(Cold_x100s_acomp, "Cold_x100s")
comp_clean(Fluc_x100s_acomp, "Fluc_x100s")
comp_clean(Opt_x100s_acomp, "Opt_x100s")

#x1000s clean
comp_clean(Cold_x1000s_acomp, "Cold_x1000s")
comp_clean(Fluc_x1000s_acomp, "Fluc_x1000s")
comp_clean(Opt_x1000s_acomp, "Opt_x1000s")


#x100s Center
comp_center(Cold_x100s_acomp_cleaned, "Cold_x100s")
comp_center(Fluc_x100s_acomp_cleaned, "Fluc_x100s")
comp_center(Opt_x100s_acomp_cleaned, "Opt_x100s")

#x1000s clean
comp_center(Cold_x1000s_acomp_cleaned, "Cold_x1000s")
comp_center(Fluc_x1000s_acomp_cleaned, "Fluc_x1000s")
comp_center(Opt_x1000s_acomp_cleaned, "Opt_x1000s")

#### Test for normality (look at graph output from comp_center as well) ####

#x100s
acompNormalGOF.test(Cold_x100s_acomp_cleaned_centered, R = 116) #sig, not normal
acompNormalGOF.test(Fluc_x100s_acomp_cleaned_centered, R = 103) #sig, not normal
acompNormalGOF.test(Opt_x100s_acomp_cleaned_centered, R = 92) #sig, not normal


#x1000s
acompNormalGOF.test(Cold_x1000s_acomp_cleaned_centered, R = 116) #sig, not normal
acompNormalGOF.test(Fluc_x1000s_acomp_cleaned_centered, R = 103) #sig, not normal
acompNormalGOF.test(Opt_x1000s_acomp_cleaned_centered, R = 92) #sig, not normal

#........ I'll just ignore those tests for now

#   MODEL 1: Cold x100s #########################################
##### Univariate Stats 101: Model Set-Up + Visualisation ###########
#Defining the independent and dependant variables
Y = Cold_x100s_acomp_cleaned_centered #centered and cleaned composition from above
X1 = Cold_x100s$Timepoint
X2 = Cold_x100s$Treatment #no need for temp-treat, as it is subsetted

#convert these to factors (and check their orders) if needed
X1 <- as.factor(X1)
X2 <- as.factor(X2)

#checking the levels orders
levels(X1)
levels(X2) 

dev.off()

#Visual exploration of day
plot(0:10, type = "n", xaxt="n", yaxt="n", bty="n", xlab = "", ylab = "")
text(5, 8, "Univariate Approaches")
text(5, 7, keyword)


#plotting the composition and color-coding/shape-coding by the factors
opar = par(xpd=NA,no.readonly=TRUE)
plot(Y,pch=c(1,3)[X2],col=c("#0041b9",
                             "#629cc4",
                             "#f4777f",
                             "#93003a")[X1]) #pch picks the symbol, col picks the color, matched it to the 2 variable X3 for pch and X1 for col with three levels
legend(x=0.85,y=0.65,abbreviate(levels(X2), minlength=1), pch=c(1,3)) #legend for X3
legend(x=0.75,y=0.65,abbreviate(levels(X1), minlength=1), pch=20,col=c("#0041b9",
                                                                       "#629cc4",
                                                                       "#f4777f",
                                                                       "#93003a"),yjust=0) #legend for X1
title("Cold Treatment (x100s)",outer=TRUE,line=-1)
#its an interesting graph, but subsetting it by day within cold might look more interesting...

 ------------

##### SUBSETTING BY DAY GRAPH SIDETRACK ##########

#Day 2: 
Day2_Y = Y[X1 == "2"] #copilot suggestion (thank u)
Day2_X1 = X1[X1 == "2"]
Day2_X2 = X2[X1 == "2"]

#Day 7: 
Day7_Y = Y[X1 == "7"] 
Day7_X1 = X1[X1 == "7"]
Day7_X2 = X2[X1 == "7"]

#Day 14: 
Day14_Y = Y[X1 == "14"] 
Day14_X1 = X1[X1 == "14"]
Day14_X2 = X2[X1 == "14"]

#Day 21
Day21_Y = Y[X1 == "21"] 
Day21_X1 = X1[X1 == "21"]
Day21_X2 = X2[X1 == "21"]

#panel the graph view to show me four graphs at a time
par(mfrow=c(2,2))

#plot no. 1
plot(Day2_Y, pch=c(1,3)[Day2_X2], col = c("blue", "red")[Day2_X2])
#plot no. 2
plot(Day7_Y, pch=c(1,3)[Day7_X2], col = c("blue", "red")[Day7_X2])
#plot no. 3
plot(Day14_Y, pch=c(1,3)[Day14_X2], col = c("blue", "red")[Day14_X2])
#plot no. 4
plot(Day21_Y, pch=c(1,3)[Day21_X2], col = c("blue", "red")[Day21_X2])
#Extras
legend(locator(1),abbreviate(levels(Sub_X2), minlength=1), pch=c(1,3), col = c("blue", "red")) #legend for X3
text(locator(1), "Day 2", col = "black")
text(locator(1), "Day 7", col = "black")
text(locator(1), "Day 14", col = "black")
text(locator(1), "Day 21", col = "black")
text(locator(1), "Change in BC composition over time (Cold Temp, x100s)")


boxplot(Y,X3,notch=TRUE) #boxplot of the data
title("Composition against ...",outer=TRUE,line=-1)


 --------
##### Univariate Stats 102: Model Making + ANOVAs #####
#e.g:
#setting the contrasts
contrasts(X2) <- "contr.treatment" #using the vaccine variable for this

fullModel = lm(ilr(Y) ~ X1 * X2) #Comp = day * treatment
coef(fullModel) #parameters
anova(fullModel) #Type 1 anova
Anova(fullModel) #Type 2 anova
#sig difference between days, not between treat and vax, no sign interaction

###NEW: Now with post-hoc testing
# i need a multivariate equivalent of Tukeys test
# Perform post-hoc tests
post_hoc <- emmeans(fullModel, pairwise ~ X1)
summary(post_hoc)


##### OPTIONAL CODE: to get the parameters of the model made in OG acomp scale #####
(a = ilrInv(coef(single_X_model)[1,],orig=Y)) #backtransforming the intercept
(b = ilrInv(rbind(0,coef(single_X_model)[-1,]),orig=Y)) #backtransforming the slope? slash other constant
#the above output shows the increase between levels of the treatment variable
(mu = a + b) #mean
(Sigma=ilrvar2clr(var(single_X_model))) #variance
vcov(model) #variance-covariance matrix
(coefs <- ilrInv(coef(fullModel),orig=Y)) #backtransforming the coefficients
barplot(coefs,las=2,xlim=c(0,11)) #barplot of the coefficients


## PLOTTING the models
#with one of the following
#Looking at the model visually
comp_ANOVA_plots_4(fullModel, Y, X1, "fullModel") #for 2 levels, for example
comp_ANOVA_plots_2(fullModel, Y, X2, "fullModel") #for 2 levels, for example
# A variable may be removed from the model if the confidence ellipse around its parameter always contains the barycenter of the composition.
#gonna have to update this code cause I am getting six ellipses


##### Univariate Stats 103: Model Diagnostics #####
#another round of plots
comp_ANOVA_diag(fullModel, Y, X1, "fullModel") #for X3 as a grouping variable, for example

#the next suggestion is to test for model redundancy by making a bunch of differently ordered models and seeing if the last variable is still significant
model1 = lm(ilr(Y) ~ X2 * X1)
anova(model1)
anova(fullModel, model1) #exactly the same





#BAr graphs attempts
mean(Y[X1 == "2" & X2 == "Vax"])[1] #mean of day 2, vaxxed groups

#make a dataset of the mean values per day and per treatment

# got it to work, had to remove the return mean data bit of the for loop
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

mean_comp_bargraaphs(Y, X1, X2)

Tester$Factorial <- paste0(Tester$Day, "_", Tester$Treatment) #make a new variable that combines the two factors

Tester_Acomp <- acomp(Tester[,c("Granu", "Lymph", "Mono")]) #make a composition of the mean data
Tester_Acomp

#money maker
barplot(Tester_Acomp, beside = FALSE, col = c( '#fee0d2','#fc9272','#de2d26'), 
        names.arg = c("D2 Con.", "D2 Chal.", "D7 Con.", "D7 Chal.", "D14 Con.", "D14 Chal.", "D21 Con.", "D21 Chal."), las = 2, ylim = c(0, 1), 
        main = "Mean Composition by Day and Treatment (CC WBC)", ylab = "Relative Contribution",
        )


#cleaner example provided by copilot
mean_comp_bargraphs <- function(Y, X1, X2) {
  mean_data <- data.frame(
    Day = rep(levels(X1), each = length(levels(X2))),
    Treatment = rep(levels(X2), times = length(levels(X1))),
    Granu = NA,
    Lymph = NA,
    Mono = NA
  )
  
  for(i in 1:nrow(mean_data)) {
    group_means <- mean(Y[X1 == mean_data$Day[i] & X2 == mean_data$Treatment[i], ])
    mean_data$Granu[i] <- group_means[1]
    mean_data$Lymph[i] <- group_means[2]
    mean_data$Mono[i]  <- group_means[3]
  }
  
  return(mean_data)
}

mean_comp_bargraphs(Y, X1, X2)
#both output the same nice, 




#   MODEL 2: Fluc x100s   ###############################################
##### Univariate Stats 101: Model Set-Up + Visualisation ###########
#Defining the independent and dependant variables
Y = Fluc_x100s_acomp_cleaned_centered #centered and cleaned composition from above
X1 = Fluc_x100s$Timepoint
X2 = Fluc_x100s$Treatment #no need for temp-treat, as it is subsetted

#convert these to factors (and check their orders) if needed
X1 <- as.factor(X1)
X2 <- as.factor(X2)

#checking the levels orders
levels(X1)
levels(X2) 

dev.off()

#Visual exploration of day
plot(0:10, type = "n", xaxt="n", yaxt="n", bty="n", xlab = "", ylab = "")
text(5, 8, "Univariate Approaches")
text(5, 7, keyword)


#plotting the composition and color-coding/shape-coding by the factors
opar = par(xpd=NA,no.readonly=TRUE)
plot(Y,pch=c(1,3)[X2],col=c("#0041b9",
                            "#629cc4",
                            "#f4777f",
                            "#93003a")[X1]) #pch picks the symbol, col picks the color, matched it to the 2 variable X3 for pch and X1 for col with three levels
legend(x=0.85,y=0.65,abbreviate(levels(X2), minlength=1), pch=c(1,3)) #legend for X3
legend(x=0.75,y=0.65,abbreviate(levels(X1), minlength=1), pch=20,col=c("#0041b9",
                                                                       "#629cc4",
                                                                       "#f4777f",
                                                                       "#93003a"),yjust=0) #legend for X1
title("Fluc Treatment (x100s)",outer=TRUE,line=-1)


##### SUBSETTING BY DAY GRAPH SIDETRACK ##########

#Day 2: 
Day2_Y = Y[X1 == "2"] #copilot suggestion (thank u)
Day2_X1 = X1[X1 == "2"]
Day2_X2 = X2[X1 == "2"]

#Day 7: 
Day7_Y = Y[X1 == "7"] 
Day7_X1 = X1[X1 == "7"]
Day7_X2 = X2[X1 == "7"]

#Day 14: 
Day14_Y = Y[X1 == "14"] 
Day14_X1 = X1[X1 == "14"]
Day14_X2 = X2[X1 == "14"]

#Day 21
Day21_Y = Y[X1 == "21"] 
Day21_X1 = X1[X1 == "21"]
Day21_X2 = X2[X1 == "21"]

#panel the graph view to show me four graphs at a time
par(mfrow=c(2,2))

#plot no. 1
plot(Day2_Y, pch=c(1,3)[Day2_X2], col = c("blue", "red")[Day2_X2])
#plot no. 2
plot(Day7_Y, pch=c(1,3)[Day7_X2], col = c("blue", "red")[Day7_X2])
#plot no. 3
plot(Day14_Y, pch=c(1,3)[Day14_X2], col = c("blue", "red")[Day14_X2])
#plot no. 4
plot(Day21_Y, pch=c(1,3)[Day21_X2], col = c("blue", "red")[Day21_X2])
#Extras
legend(locator(1),abbreviate(levels(Sub_X2), minlength=1), pch=c(1,3), col = c("blue", "red")) #legend for X3
text(locator(1), "Day 2", col = "black")
text(locator(1), "Day 7", col = "black")
text(locator(1), "Day 14", col = "black")
text(locator(1), "Day 21", col = "black")
text(locator(1), "Change in BC composition over time (Fluc Temp, x100s)")

#reset graph view
par(mfrow=c(1,1))

#### Univariate Stats 102: Model Making + ANOVAs #####
#e.g:
#setting the contrasts
contrasts(X2) <- "contr.treatment" #using the vaccine variable for this

fullModel = lm(ilr(Y) ~ X1 * X2) #Comp = day * treatment
coef(fullModel) #parameters
anova(fullModel)
Anova(fullModel) #Type 2 anova
#sig difference between days, not between treat and vax, no sign interaction

###NEW: Now with post-hoc testing
# i need a multivariate equivalent of Tukeys test
# Perform post-hoc tests
post_hoc <- emmeans(fullModel, pairwise ~ X1)
summary(post_hoc)


##### OPTIONAL CODE: to get the parameters of the model made in OG acomp scale #####
(a = ilrInv(coef(single_X_model)[1,],orig=Y)) #backtransforming the intercept
(b = ilrInv(rbind(0,coef(single_X_model)[-1,]),orig=Y)) #backtransforming the slope? slash other constant
#the above output shows the increase between levels of the treatment variable
(mu = a + b) #mean
(Sigma=ilrvar2clr(var(single_X_model))) #variance
vcov(model) #variance-covariance matrix
(coefs <- ilrInv(coef(fullModel),orig=Y)) #backtransforming the coefficients
barplot(coefs,las=2,xlim=c(0,11)) #barplot of the coefficients


## PLOTTING the models
#with one of the following
#Looking at the model visually
comp_ANOVA_plots_4(fullModel, Y, X1, "fullModel") #for 2 levels, for example
comp_ANOVA_plots_2(fullModel, Y, X2, "fullModel") #for 2 levels, for example
# A variable may be removed from the model if the confidence ellipse around its parameter always contains the barycenter of the composition.
#gonna have to update this code cause I am getting six ellipses


##### Univariate Stats 103: Model Diagnostics #####
#another round of plots
comp_ANOVA_diag(fullModel, Y, X1, "fullModel") #for X3 as a grouping variable, for example

#the next suggestion is to test for model redundancy by making a bunch of differently ordered models and seeing if the last variable is still significant
model1 = lm(ilr(Y) ~ X2 * X1)
anova(model1) #different significant results
anova(fullModel, model1) #exactly the same zapparently





#   MODEL 3: Opt x100s   ###############################################
##### Univariate Stats 101: Model Set-Up + Visualisation ###########
#Defining the independent and dependant variables
Y = Opt_x100s_acomp_cleaned_centered #centered and cleaned composition from above
X1 = Opt_x100s$Timepoint
X2 = Opt_x100s$Treatment #no need for temp-treat, as it is subsetted

#convert these to factors (and check their orders) if needed
X1 <- as.factor(X1)
X2 <- as.factor(X2)

#checking the levels orders
levels(X1)
levels(X2) 

dev.off()

#Visual exploration of day
plot(0:10, type = "n", xaxt="n", yaxt="n", bty="n", xlab = "", ylab = "")
text(5, 8, "Univariate Approaches")
text(5, 7, keyword)


#plotting the composition and color-coding/shape-coding by the factors
opar = par(xpd=NA,no.readonly=TRUE)
plot(Y,pch=c(1,3)[X2],col=c("#0041b9",
                            "#629cc4",
                            "#f4777f",
                            "#93003a")[X1]) #pch picks the symbol, col picks the color, matched it to the 2 variable X3 for pch and X1 for col with three levels
legend(x=0.85,y=0.65,abbreviate(levels(X2), minlength=1), pch=c(1,3)) #legend for X3
legend(x=0.75,y=0.65,abbreviate(levels(X1), minlength=1), pch=20,col=c("#0041b9",
                                                                       "#629cc4",
                                                                       "#f4777f",
                                                                       "#93003a"),yjust=0) #legend for X1
title("Opt Treatment (x100s)",outer=TRUE,line=-1)


#### SUBSETTING BY DAY GRAPH SIDETRACK ##########

#subsetting by day (tester section)
Sub_Y = Y[X1 == "7"] #copilot suggestion
Sub_Y #that should not have worked in my brain, but ok
Sub_X1 = X1[X1 == "7"]
Sub_X1
Sub_X2 = X2[X1 == "7"]
Sub_X2 #manually check this one later

#more tester section
plot(Sub_Y,col=c("#0041b9", "black","red", "green")[Sub_X1]) #now its only doing the color corresponding to 7, very good
plot(Sub_Y, pch=c(1,3)[Sub_X2])
legend(x=0.85,y=0.65,abbreviate(levels(Sub_X2), minlength=1), pch=c(1,3)) #legend for X3
title("Day 7, Opt Temp",outer=TRUE,line=-1)



#Day 2: 
Day2_Y = Y[X1 == "2"] #copilot suggestion (thank u)
Day2_X1 = X1[X1 == "2"]
Day2_X2 = X2[X1 == "2"]

#Day 7: 
Day7_Y = Y[X1 == "7"] 
Day7_X1 = X1[X1 == "7"]
Day7_X2 = X2[X1 == "7"]

#Day 14: 
Day14_Y = Y[X1 == "14"] 
Day14_X1 = X1[X1 == "14"]
Day14_X2 = X2[X1 == "14"]

#Day 21
Day21_Y = Y[X1 == "21"] 
Day21_X1 = X1[X1 == "21"]
Day21_X2 = X2[X1 == "21"]

#panel the graph view to show me four graphs at a time
par(mfrow=c(2,2))

#plot no. 1
plot(Day2_Y, pch=c(1,3)[Day2_X2], col = c("blue", "red")[Day2_X2])
#plot no. 2
plot(Day7_Y, pch=c(1,3)[Day7_X2], col = c("blue", "red")[Day7_X2])
#plot no. 3
plot(Day14_Y, pch=c(1,3)[Day14_X2], col = c("blue", "red")[Day14_X2])
#plot no. 4
plot(Day21_Y, pch=c(1,3)[Day21_X2], col = c("blue", "red")[Day21_X2])
#Extras
legend(locator(1),abbreviate(levels(Sub_X2), minlength=1), pch=c(1,3), col = c("blue", "red")) #legend for X3
text(locator(1), "Day 2", col = "black")
text(locator(1), "Day 7", col = "black")
text(locator(1), "Day 14", col = "black")
text(locator(1), "Day 21", col = "black")
text(locator(1), "Change in BC composition over time (Opt Temp, x100s)")


boxplot(Y,X3,notch=TRUE) #boxplot of the data
title("Composition against ...",outer=TRUE,line=-1)



#### Univariate Stats 102: Model Making + ANOVAs #####
#e.g:
#setting the contrasts
contrasts(X2) <- "contr.treatment" #using the vaccine variable for this

fullModel = lm(ilr(Y) ~ X1 * X2) #Comp = day * treatment
coef(fullModel) #parameters
anova(fullModel)
Anova(fullModel) #Type 2 anova
#sig difference between days, not between treat and vax, no sign interaction

###NEW: Now with post-hoc testing
# i need a multivariate equivalent of Tukeys test
# Perform post-hoc tests
post_hoc <- emmeans(fullModel, pairwise ~ X2 | X1)
summary(post_hoc)





##### OPTIONAL CODE: to get the parameters of the model made in OG acomp scale #####
(a = ilrInv(coef(single_X_model)[1,],orig=Y)) #backtransforming the intercept
(b = ilrInv(rbind(0,coef(single_X_model)[-1,]),orig=Y)) #backtransforming the slope? slash other constant
#the above output shows the increase between levels of the treatment variable
(mu = a + b) #mean
(Sigma=ilrvar2clr(var(single_X_model))) #variance
vcov(model) #variance-covariance matrix
(coefs <- ilrInv(coef(fullModel),orig=Y)) #backtransforming the coefficients
barplot(coefs,las=2,xlim=c(0,11)) #barplot of the coefficients


## PLOTTING the models
#with one of the following
#Looking at the model visually
comp_ANOVA_plots_4(fullModel, Y, X1, "fullModel") #for 2 levels, for example
comp_ANOVA_plots_2(fullModel, Y, X2, "fullModel") #for 2 levels, for example
# A variable may be removed from the model if the confidence ellipse around its parameter always contains the barycenter of the composition.
#gonna have to update this code cause I am getting six ellipses


##### Univariate Stats 103: Model Diagnostics #####
#another round of plots
comp_ANOVA_diag(fullModel, Y, X1, "fullModel") #for X3 as a grouping variable, for example

#the next suggestion is to test for model redundancy by making a bunch of differently ordered models and seeing if the last variable is still significant
model1 = lm(ilr(Y) ~ X2 * X1)
anova(model1) #same significant results
anova(fullModel, model1) #exactly the same apparently


#   MODEL 4: Cold x1000s ###############################################
##### Univariate Stats 101: Model Set-Up + Visualisation ###########
#Defining the independent and dependant variables
Y = Cold_x1000s_acomp_cleaned_centered #centered and cleaned composition from above
X1 = Cold_x1000s$Timepoint
X2 = Cold_x1000s$Treatment #no need for temp-treat, as it is subsetted

#convert these to factors (and check their orders) if needed
X1 <- as.factor(X1)
X2 <- as.factor(X2)

#checking the levels orders
levels(X1)
levels(X2) 

dev.off()

#Visual exploration of day
plot(0:10, type = "n", xaxt="n", yaxt="n", bty="n", xlab = "", ylab = "")
text(5, 8, "Univariate Approaches")
text(5, 7, keyword)


#plotting the composition and color-coding/shape-coding by the factors
opar = par(xpd=NA,no.readonly=TRUE)
plot(Y,pch=c(1,3)[X2],col=c("#0041b9",
                            "#629cc4",
                            "#f4777f",
                            "#93003a")[X1]) #pch picks the symbol, col picks the color, matched it to the 2 variable X3 for pch and X1 for col with three levels
legend((locator(1)),abbreviate(levels(X2), minlength=1), pch=c(1,3)) #legend for X3
legend((locator(1)),abbreviate(levels(X1), minlength=1), pch=20,col=c("#0041b9",
                                                                       "#629cc4",
                                                                       "#f4777f",
                                                                       "#93003a"),yjust=0) #legend for X1
title("Cold Treatment (x1000s)",outer=TRUE,line=-1)


##### SUBSETTING BY DAY GRAPH SIDETRACK ##########

#Day 2: 
Day2_Y = Y[X1 == "2"] #copilot suggestion (thank u)
Day2_X1 = X1[X1 == "2"]
Day2_X2 = X2[X1 == "2"]

#Day 7: 
Day7_Y = Y[X1 == "7"] 
Day7_X1 = X1[X1 == "7"]
Day7_X2 = X2[X1 == "7"]

#Day 14: 
Day14_Y = Y[X1 == "14"] 
Day14_X1 = X1[X1 == "14"]
Day14_X2 = X2[X1 == "14"]

#Day 21
Day21_Y = Y[X1 == "21"] 
Day21_X1 = X1[X1 == "21"]
Day21_X2 = X2[X1 == "21"]


#plot no. 1
plot(Day2_Y, pch=c(1,3)[Day2_X2], col = c("blue", "red")[Day2_X2])
text(locator(1), "Day 2", col = "black")

#plot no. 2
plot(Day7_Y, pch=c(1,3)[Day7_X2], col = c("blue", "red")[Day7_X2])
text(locator(1), "Day 7", col = "black")

#plot no. 3
plot(Day14_Y, pch=c(1,3)[Day14_X2], col = c("blue", "red")[Day14_X2])
text(locator(1), "Day 14", col = "black")

#plot no. 4
plot(Day21_Y, pch=c(1,3)[Day21_X2], col = c("blue", "red")[Day21_X2])
text(locator(1), "Day 21", col = "black")

#Extras
legend(locator(1),abbreviate(levels(Sub_X2), minlength=1), pch=c(1,3), col = c("blue", "red")) #legend for X3








#### Univariate Stats 102: Model Making + ANOVAs #####
#e.g:
#setting the contrasts
contrasts(X2) <- "contr.treatment" #using the vaccine variable for this

fullModel = lm(ilr(Y) ~ X1 * X2) #Comp = day * treatment
coef(fullModel) #parameters
anova(fullModel)
Anova(fullModel) #Type 2 anova
#sig difference between days, not between treat and vax, no sign interaction


###NEW: Now with post-hoc testing
# i need a multivariate equivalent of Tukeys test
# Perform post-hoc tests
post_hoc <- emmeans(fullModel, pairwise ~ X1)
summary(post_hoc)





##### OPTIONAL CODE: to get the parameters of the model made in OG acomp scale #####
(a = ilrInv(coef(single_X_model)[1,],orig=Y)) #backtransforming the intercept
(b = ilrInv(rbind(0,coef(single_X_model)[-1,]),orig=Y)) #backtransforming the slope? slash other constant
#the above output shows the increase between levels of the treatment variable
(mu = a + b) #mean
(Sigma=ilrvar2clr(var(single_X_model))) #variance
vcov(model) #variance-covariance matrix
(coefs <- ilrInv(coef(fullModel),orig=Y)) #backtransforming the coefficients
barplot(coefs,las=2,xlim=c(0,11)) #barplot of the coefficients


## PLOTTING the models
#with one of the following
#Looking at the model visually
comp_ANOVA_plots_4(fullModel, Y, X1, "fullModel") #for 2 levels, for example
comp_ANOVA_plots_2(fullModel, Y, X2, "fullModel") #for 2 levels, for example
# A variable may be removed from the model if the confidence ellipse around its parameter always contains the barycenter of the composition.
#gonna have to update this code cause I am getting six ellipses


##### Univariate Stats 103: Model Diagnostics #####
#another round of plots
comp_ANOVA_diag(fullModel, Y, X1, "fullModel") #for X3 as a grouping variable, for example



#   MODEL 5: Fluc x1000s ###############################################
##### Univariate Stats 101: Model Set-Up + Visualisation ###########
#Defining the independent and dependant variables
Y = Fluc_x1000s_acomp_cleaned_centered #centered and cleaned composition from above
X1 = Fluc_x1000s$Timepoint
X2 = Fluc_x1000s$Treatment #no need for temp-treat, as it is subsetted

#convert these to factors (and check their orders) if needed
X1 <- as.factor(X1)
X2 <- as.factor(X2)

#checking the levels orders
levels(X1)
levels(X2) 

dev.off()

#Visual exploration of day
plot(0:10, type = "n", xaxt="n", yaxt="n", bty="n", xlab = "", ylab = "")
text(5, 8, "Univariate Approaches")
text(5, 7, keyword)


#plotting the composition and color-coding/shape-coding by the factors
opar = par(xpd=NA,no.readonly=TRUE)
plot(Y,pch=c(1,3)[X2],col=c("#0041b9",
                            "#629cc4",
                            "#f4777f",
                            "#93003a")[X1]) #pch picks the symbol, col picks the color, matched it to the 2 variable X3 for pch and X1 for col with three levels
legend((locator(1)),abbreviate(levels(X2), minlength=1), pch=c(1,3)) #legend for X3
legend((locator(1)),abbreviate(levels(X1), minlength=1), pch=20,col=c("#0041b9",
                                                                      "#629cc4",
                                                                      "#f4777f",
                                                                      "#93003a"),yjust=0) #legend for X1
title("Fluc Treatment (x1000s)",outer=TRUE,line=-1)


##### SUBSETTING BY DAY GRAPH SIDETRACK ##########

#Day 2: 
Day2_Y = Y[X1 == "2"] #copilot suggestion (thank u)
Day2_X1 = X1[X1 == "2"]
Day2_X2 = X2[X1 == "2"]

#Day 7: 
Day7_Y = Y[X1 == "7"] 
Day7_X1 = X1[X1 == "7"]
Day7_X2 = X2[X1 == "7"]

#Day 14: 
Day14_Y = Y[X1 == "14"] 
Day14_X1 = X1[X1 == "14"]
Day14_X2 = X2[X1 == "14"]

#Day 21
Day21_Y = Y[X1 == "21"] 
Day21_X1 = X1[X1 == "21"]
Day21_X2 = X2[X1 == "21"]


#plot no. 1
plot(Day2_Y, pch=c(1,3)[Day2_X2], col = c("blue", "red")[Day2_X2])
text(locator(1), "Day 2", col = "black")

#plot no. 2
plot(Day7_Y, pch=c(1,3)[Day7_X2], col = c("blue", "red")[Day7_X2])
text(locator(1), "Day 7", col = "black")

#plot no. 3
plot(Day14_Y, pch=c(1,3)[Day14_X2], col = c("blue", "red")[Day14_X2])
text(locator(1), "Day 14", col = "black")

#plot no. 4
plot(Day21_Y, pch=c(1,3)[Day21_X2], col = c("blue", "red")[Day21_X2])
text(locator(1), "Day 21", col = "black")

#Extras
legend(locator(1),abbreviate(levels(Sub_X2), minlength=1), pch=c(1,3), col = c("blue", "red")) #legend for X3



#### Univariate Stats 102: Model Making + ANOVAs #####
#e.g:
#setting the contrasts
contrasts(X2) <- "contr.treatment" #using the vaccine variable for this

fullModel = lm(ilr(Y) ~ X1 * X2) #Comp = day * treatment
coef(fullModel) #parameters
anova(fullModel)
Anova(fullModel) #Type 2 anova
#sig difference between days, not between treat and vax, no sign interaction


###NEW: Now with post-hoc testing
# i need a multivariate equivalent of Tukeys test
# Perform post-hoc tests
post_hoc <- emmeans(fullModel, pairwise ~ X1)
summary(post_hoc)



##### OPTIONAL CODE: to get the parameters of the model made in OG acomp scale #####
(a = ilrInv(coef(single_X_model)[1,],orig=Y)) #backtransforming the intercept
(b = ilrInv(rbind(0,coef(single_X_model)[-1,]),orig=Y)) #backtransforming the slope? slash other constant
#the above output shows the increase between levels of the treatment variable
(mu = a + b) #mean
(Sigma=ilrvar2clr(var(single_X_model))) #variance
vcov(model) #variance-covariance matrix
(coefs <- ilrInv(coef(fullModel),orig=Y)) #backtransforming the coefficients
barplot(coefs,las=2,xlim=c(0,11)) #barplot of the coefficients


## PLOTTING the models
#with one of the following
#Looking at the model visually
comp_ANOVA_plots_4(fullModel, Y, X1, "fullModel") #for 2 levels, for example
comp_ANOVA_plots_2(fullModel, Y, X2, "fullModel") #for 2 levels, for example
# A variable may be removed from the model if the confidence ellipse around its parameter always contains the barycenter of the composition.
#gonna have to update this code cause I am getting six ellipses


##### Univariate Stats 103: Model Diagnostics #####
#another round of plots
comp_ANOVA_diag(fullModel, Y, X1, "fullModel") #for X3 as a grouping variable, for example




#   MODEL 6: Opt x1000s ###############################################
##### Univariate Stats 101: Model Set-Up + Visualisation ###########
#Defining the independent and dependant variables
Y = Opt_x1000s_acomp_cleaned_centered #centered and cleaned composition from above
X1 = Opt_x1000s$Timepoint
X2 = Opt_x1000s$Treatment #no need for temp-treat, as it is subsetted

#convert these to factors (and check their orders) if needed
X1 <- as.factor(X1)
X2 <- as.factor(X2)

#checking the levels orders
levels(X1)
levels(X2) 

dev.off()

#Visual exploration of day
plot(0:10, type = "n", xaxt="n", yaxt="n", bty="n", xlab = "", ylab = "")
text(5, 8, "Univariate Approaches")
text(5, 7, keyword)


#plotting the composition and color-coding/shape-coding by the factors
opar = par(xpd=NA,no.readonly=TRUE)
plot(Y,pch=c(1,3)[X2],col=c("#0041b9",
                            "#629cc4",
                            "#f4777f",
                            "#93003a")[X1]) #pch picks the symbol, col picks the color, matched it to the 2 variable X3 for pch and X1 for col with three levels
legend((locator(1)),abbreviate(levels(X2), minlength=1), pch=c(1,3)) #legend for X3
legend((locator(1)),abbreviate(levels(X1), minlength=1), pch=20,col=c("#0041b9",
                                                                      "#629cc4",
                                                                      "#f4777f",
                                                                      "#93003a"),yjust=0) #legend for X1
title("Opt Treatment (x1000s)",outer=TRUE,line=-1)


##### SUBSETTING BY DAY GRAPH SIDETRACK ##########

#Day 2: 
Day2_Y = Y[X1 == "2"] #copilot suggestion (thank u)
Day2_X1 = X1[X1 == "2"]
Day2_X2 = X2[X1 == "2"]

#Day 7: 
Day7_Y = Y[X1 == "7"] 
Day7_X1 = X1[X1 == "7"]
Day7_X2 = X2[X1 == "7"]

#Day 14: 
Day14_Y = Y[X1 == "14"] 
Day14_X1 = X1[X1 == "14"]
Day14_X2 = X2[X1 == "14"]

#Day 21
Day21_Y = Y[X1 == "21"] 
Day21_X1 = X1[X1 == "21"]
Day21_X2 = X2[X1 == "21"]


#plot no. 1
plot(Day2_Y, pch=c(1,3)[Day2_X2], col = c("blue", "red")[Day2_X2])
text(locator(1), "Day 2", col = "black")

#plot no. 2
plot(Day7_Y, pch=c(1,3)[Day7_X2], col = c("blue", "red")[Day7_X2])
text(locator(1), "Day 7", col = "black")

#plot no. 3
plot(Day14_Y, pch=c(1,3)[Day14_X2], col = c("blue", "red")[Day14_X2])
text(locator(1), "Day 14", col = "black")

#plot no. 4
plot(Day21_Y, pch=c(1,3)[Day21_X2], col = c("blue", "red")[Day21_X2])
text(locator(1), "Day 21", col = "black")

#Extras
legend(locator(1),abbreviate(levels(Sub_X2), minlength=1), pch=c(1,3), col = c("blue", "red")) #legend for X3


#### Univariate Stats 102: Model Making + ANOVAs #####
#e.g:
#setting the contrasts
contrasts(X2) <- "contr.treatment" #using the vaccine variable for this

fullModel = lm(ilr(Y) ~ X1 * X2) #Comp = day * treatment
coef(fullModel) #parameters
anova(fullModel)
Anova(fullModel) #Type 2 anova
#sig difference between days, not between treat and vax, no sign interaction


###NEW: Now with post-hoc testing
# i need a multivariate equivalent of Tukeys test
# Perform post-hoc tests
post_hoc <- emmeans(fullModel, pairwise ~ X1)
summary(post_hoc)





##### OPTIONAL CODE: to get the parameters of the model made in OG acomp scale #####
(a = ilrInv(coef(single_X_model)[1,],orig=Y)) #backtransforming the intercept
(b = ilrInv(rbind(0,coef(single_X_model)[-1,]),orig=Y)) #backtransforming the slope? slash other constant
#the above output shows the increase between levels of the treatment variable
(mu = a + b) #mean
(Sigma=ilrvar2clr(var(single_X_model))) #variance
vcov(model) #variance-covariance matrix
(coefs <- ilrInv(coef(fullModel),orig=Y)) #backtransforming the coefficients
barplot(coefs,las=2,xlim=c(0,11)) #barplot of the coefficients


## PLOTTING the models
#with one of the following
#Looking at the model visually
comp_ANOVA_plots_4(fullModel, Y, X1, "fullModel") #for 2 levels, for example
comp_ANOVA_plots_2(fullModel, Y, X2, "fullModel") #for 2 levels, for example
# A variable may be removed from the model if the confidence ellipse around its parameter always contains the barycenter of the composition.
#gonna have to update this code cause I am getting six ellipses


##### Univariate Stats 103: Model Diagnostics #####
#another round of plots
comp_ANOVA_diag(fullModel, Y, X1, "fullModel") #for X3 as a grouping variable, for example



##### Then: Mulivariate Stats, if wanted #####