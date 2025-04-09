###############################################
# Full Model Play Code                       ##
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


###### DATA #########
DATA <- as.data.frame(Bloodsmear_data_full)
#$ [Optional: Data Cleaning]
DATA$Timepoint <- replace(DATA$Timepoint, DATA$Timepoint == 17, 7) #okay that worked

#Turning thins into factors
DATA$Timepoint <- as.factor(DATA$Timepoint)
levels(DATA$Timepoint)
DATA$TempTreat <- as.factor(DATA$TempTreat)
levels(DATA$TempTreat)

Full_x100s = DATA[c(1:6, 11:13)]
Full_x1000s = DATA[1:10]


# make a composition
Full_x100s_comp <- acomp(Full_x100s[7:9])
Full_x1000s_comp <- acomp(Full_x1000s[7:10])

#explorative graphs
comp_explore(Full_x100s_comp, "Full_x100s")
comp_explore(Full_x1000s_comp, "Full_x1000s")

#cleaning
comp_clean(Full_x100s_comp, "Full_x100s")
comp_clean(Full_x1000s_comp, "Full_x1000s")

#center
comp_center(Full_x100s_comp_cleaned, "Full_x100s")
comp_center(Full_x1000s_comp_cleaned, "Full_x1000s")


#   MODEL 1: Full x100s #########################################
##### Univariate Stats 101: Model Set-Up + Visualisation ###########
#Defining the independent and dependant variables
Y = Full_x100s_comp_cleaned_centered #centered and cleaned composition from above
X1 = Full_x100s$Timepoint
X2 = Full_x100s$Treatment #no need for temp-treat, as it is subsetted
X3 = Full_x100s$TempTreat

#convert these to factors (and check their orders) if needed
X1 <- as.factor(X1)
X2 <- as.factor(X2)
X3 <- as.factor(X3)

#checking the levels orders
levels(X1)
levels(X2) 
levels(X3)

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
                                                                       "#93003a"),yjust=0) 
title("Full Data (x100s)",outer=TRUE,line=-1)

#ANOVA time
contrasts(X2) <- "contr.treatment" #using the vaccine variable for this

fullModel = lm(ilr(Y) ~ X1 * X2 * X3) #full model
coef(fullModel) #parameters
anova(fullModel) #Type 1 anova
Anova(fullModel) #Type 2 anova
#sig difference between days, not

comp_ANOVA_diag(fullModel, Y, X2, "fullModel")


#   MODEL 2: Full x1000s #########################################
##### Univariate Stats 101: Model Set-Up + Visualisation ###########
#Defining the independent and dependant variables
Y = Full_x1000s_comp_cleaned_centered #centered and cleaned composition from above
X1 = Full_x1000s$Timepoint
X2 = Full_x1000s$Treatment #no need for temp-treat, as it is subsetted
X3 = Full_x1000s$TempTreat

#convert these to factors (and check their orders) if needed
X1 <- as.factor(X1)
X2 <- as.factor(X2)
X3 <- as.factor(X3)

#checking the levels orders
levels(X1)
levels(X2) 
levels(X3)

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
                                                                       "#93003a"),yjust=0) 
title("Full Data (x1000s)",outer=TRUE,line=-1)

#ANOVA time
contrasts(X2) <- "contr.treatment" #using the vaccine variable for this

fullModel = lm(ilr(Y) ~ X1 * X2 * X3) #full model
coef(fullModel) #parameters
anova(fullModel) #Type 1 anova
Anova(fullModel) #Type 2 anova
#sig difference between days, not

comp_ANOVA_diag(fullModel, Y, X2, "fullModel")
