###############################################
# Function-based Compositional Analysis Code ##
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

comp_ANOVA_diag <- function(model, Y, Group, identifier) {
  
  diag_Ident <- paste0(identifier, "_anova_diagnostics_graphs_") 
  pdf(file = paste(diag_Ident, ".pdf", sep = "_"))
  plot(0:10, type = "n", xaxt="n", yaxt="n", bty="n", xlab = "", ylab = "")
  text(5, 8, "ANOVA diagnostics plots")
  text(5, 7, identifier)
  
  Resid = ilrInv(resid(model),orig=Y)
  plot(Resid, cex=0.5)
  title("Ternary Diagrams of Residuals",outer=TRUE,line=-1)
  
  boxplot(Resid) #also reveals outliers
  title("Boxplots of Residuals",outer=TRUE,line=-1)
  
  #QQ plot (also normally distributed residuals)
  qqnorm(Resid) 
  
  #homoschedascity (equal variance)
  Pred = ilrInv(predict(model),orig=Y)
  opar <- par(oma=c(3,3,0,0),mar=c(4,4,1,0))
  pairwisePlot(clr(Pred),clr(Resid))
  mtext(text=c("predicted values (clr)","residuals (clr)"),
        side=c(1,2),at=0.5,line=2,outer=TRUE)
  par(opar)
  
  #color-coordinate the above diagram by group
  opar <- par(oma=c(0,0,2,0),mar=c(4,4,1,0))
  pairwisePlot(clr(Pred),clr(Resid),col=as.numeric(X1))
  legend(locator(1),levels(Group),col=as.numeric(1:3),pch=1,xpd=NA) #starts a click to place mechanism for the legend
  par(opar)
  
  dev.off()
  print("Graphs Made!")
}




#### After packages and functions, add and structure your data ####






#### Then, use the custom functions to explore, clean and center your comps #####





#### Test for normality (look at graph output from comp_center as well) ####

acompNormalGOF.test(..., R = 116)


#### Univariate Stats 101: Model Set-Up + Visualisation #####
#e.g:
#Defining the independent and dependant variables
Y = Spooky_centered_acomp #centered and cleaned composition from above
X1 = DATA$TempTreat
X2 = DATA$Timepoint
X3 = DATA$Treatment

#convert these to factors (and check their orders) if needed
X1 <- as.factor(X1)
X2 <- as.factor(X2)
X3 <- as.factor(X3)

#checking the levels orders
levels(X1)
levels(X2) 
levels(X3)

plot(0:10, type = "n", xaxt="n", yaxt="n", bty="n", xlab = "", ylab = "")
text(5, 8, "Univariate Approaches")
text(5, 7, Keyword)


#Visual exploration of your groups
#plotting the composition and color-coding/shape-coding by the factors
opar = par(xpd=NA,no.readonly=TRUE)
plot(Y,pch=c(3,20)[X3],col=c("blue","green","red")[X1]) #pch picks the symbol, col picks the color, matched it to the 2 variable X3 for pch and X1 for col with three levels
legend(x=0.85,y=0.65,abbreviate(levels(X3), minlength=1), pch=c(3,20)) #legend for X3
legend(x=0.75,y=0.65,abbreviate(levels(X1), minlength=1), pch=20,col=c("blue","green","red"),yjust=0) #legend for X1
title("Composition by factor",outer=TRUE,line=-1)

boxplot(Y,X3,notch=TRUE) #boxplot of the data
title("Composition against ...",outer=TRUE,line=-1)

#### Univariate Stats 102: Model Making + ANOVAs #####
#e.g:
#setting the contrasts
contrasts(X3) <- "contr.treatment" #using the vaccine variable for this

fullModel = lm(ilr(Y) ~ X1 + X2 + X3) #model, the first level of each factor is assumed to be zero
coef(fullModel) #parameters
anova(fullModel)

## OPTIONAL CODE: to get the parameters of the model made in OG acomp scale
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
comp_ANOVA_plots_2(fullModel, Y, X3, "fullModel") #for 2 levels, for example
comp_ANOVA_plots_3()
comp_ANOVA_plots_4()
# A variable may be removed from the model if the confidence ellipse around its parameter always contains the barycenter of the composition.


##### Univariate Stats 103: Model Diagnostics #####
#another round of plots
comp_ANOVA_diag(fullModel, Y, X3, "fullModel") #for X3 as a grouping variable, for example

#the next suggestion is to test for model redundancy by making a bunch of differently ordered models and seeing if the last variable is still significant
model1 = lm(ilr(Y) ~ X3 + X2 + X1)
anova(model1)
model2 = lm(ilr(Y) ~ X3 + X1 + X2)
anova(model2)
model3 = lm(ilr(Y) ~ X2 + X3 + X1)
anova(model3)

##### Then: Mulivariate Stats, if wanted #####