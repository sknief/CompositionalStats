###############################################
## Clean, Basic, Compositional Analysis Code ##
## SMSK 2025                                 ##
###############################################

#### First: Packages and Prereqs
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

#### Step One: Setting up exports and the whole code
#All variables, output and models start with the term "Spooky" - you need to ctrl + f and replace this term with a keyword relevant to your current research question

#Let's save both of these pieces of information
researchQuestion <- "Enter Research Q" #keep brief
keyword <- "Spooky" #Spooky by default


#we will use the keyword in the export pdf titles as well
PDF_Ident <- paste0("BaseCompAnalysis_", keyword, "_") #this is to identify the output from your specific run or hypothesis
pdf(file = paste(PDF_Ident, "COLLAT.pdf", sep = "_"))

#Let's give the PDF some basic text, some info about scales for example
plot.new()
ScaleHeader <- "A note on scale"
ScaleBody1 <- "We will use the Aitchison compositional scale (acomp):"
ScaleBody2 <- "1. The stats textbook recommends using acomp for count data,"
ScaleBody3 <- "2. The underlying distribution is likely acomp"
ScaleBody4 <- "3. The total of the composition is large enough to avoid"
ScaleBody5 <- "random selection error in counts"

text(x=.5, y=1, ScaleHeader, font=4, cex=2, col="blue")
text(x=.5, y=.7, ScaleBody1)  # first 2 numbers are xy-coordinates within [0, 1]
text(x=.5, y=.6, ScaleBody2)
text(x=.5, y=.5, ScaleBody3)
text(x=.5, y=.4, ScaleBody4)
text(x=.5, y=.3, ScaleBody5)

# we can also note our research question and keyword
plot.new()
text(x=.5, y=1, keyword, font=3, cex=2, col="blue")
text(x=.5, y=.7, researchQuestion)  

#dev.off() #this closes the pdf file (shuts off the device)
#going forward, most plots will have the export to pdf text at the top and bottom of the relevant section, but commented out so that you can visualize plots in R first



#### Step Two: Data Import
#$ [Insert Data Here]
DATA <- as.dataframe(yourdatahere)
#$ [Optional: Data Cleaning]

# make a composition
Spooky_acomp = acomp(DATASET, c('Y1', 'Y2', 'Y3'))
Spooky_acomp

#explore the composition

#PDF export code
plot(0:10, type = "n", xaxt="n", yaxt="n", bty="n", xlab = "", ylab = "")
text(5, 8, "Graphical Exploration of Composition")
text(5, 7, Keyword)

plot(Spooky_acomp)
title("Untransformed Composition",outer=TRUE,line=-1)
plot(clr(Spooky_acomp)) #clr transformation (central log)
title("CLR-transformed Composition",outer=TRUE,line=-1)
plot(ilr(Spooky_acomp)) #ilr transformation (isometric log) (check)
title("ILR-transformed Composition",outer=TRUE,line=-1)


#### Step Three: Checking for and dealing with missing values
# notice: you need to work with an acomp object for the following code
missingSummary(Spooky_acomp) #table of missing values
# we will use the imputation method of the book in order to deal with missing values, which is to replace them with the smallest non-zero value in the composition (NOT OG DATASET) 
BDL_MNAR_REPLACE <- 0.001 #assign a value to replace BDL and MNARs
Spooky_acomp <- zeroreplace(Spooky_acomp, BDL_MNAR_REPLACE) #replace zeros
missingSummary(Spooky_acomp) #check for missing values again
Spooky_acomp[is.na(Spooky_acomp)] <- BDL_MNAR_REPLACE #replace NAs
missingSummary(Spooky_acomp) # final check

#explore the composition again
#PDF code
plot(0:10, type = "n", xaxt="n", yaxt="n", bty="n", xlab = "", ylab = "")
text(5, 8, "Composition after treating missing values")
text(5, 7, Keyword)


plot(Spooky_acomp)
title("Untransformed Composition",outer=TRUE,line=-1)
plot(clr(Spooky_acomp)) #clr transformation (central log)
title("CLR-transformed Composition",outer=TRUE,line=-1)
plot(ilr(Spooky_acomp)) #ilr transformation (isometric log) (check)
title("ILR-transformed Composition",outer=TRUE,line=-1)


### Step Five: Testing for normality + other assumptions (DATA diagnostics)
# First, testing for normality (normal distribution)
mvnorm.etest(ilr(Spooky_acomp)) #complete multivariate test
acompNormalGOF.test(Spooky_acomp, R = 311) #alternative multivariate test, works on non-transformed data

#exploring through graphs
#better panel making code from stackoverflow
panel.qq <- function(x, y, ...) {
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1), new = TRUE)
  qqplot(x, y, xlab = deparse(substitute(x)), ylab = deparse(substitute(y)))
  abline(c(0,1), ...)
}

#QQplots
pairs(Spooky_acomp, lower.panel = panel.qq) #use acomp object for these


#Alternative tests to check if a count comp distribution fits better 
Spooky_ccomp <- ccomp(DATASET, c('Y1', 'Y2', 'Y3'))
Spooky_ccomp
Spooky_ccomp #explore the composition

#Alternative ccomp distributions: multinomial and multi-poisson
ccompMultinomialGOF.test(Spooky_ccomp, R = 310) #not multinomially distributed, no use in using ccomp
ccompPoissonGOF.test(Spooky_ccomp) # not a multi-Poission distribution

#If composition failed normality test, we can try to transform or center
#centering the (acomp)
mean(Spooky_acomp)
Spooky_centered_acomp <- Spooky_acomp - mean(Spooky_acomp) #centering by subtracting the mean
mean(Spooky_centered_acomp) #check if the centered mean is the neutral element

#examining centered data
plot(Spooky_centered_acomp)
pairs(Spooky_centered_acomp, lower.panel = panel.qq)
mvnorm.etest(Spooky_centered_acomp, R = 310) #complete multivariate test, 
#according to the book, scaling is often used to compare subcompositions, but centering is used when date is often squashed into one corner of the simplex, which is what we had 

#PDF export code for this section: 
plot(0:10, type = "n", xaxt="n", yaxt="n", bty="n", xlab = "", ylab = "")
text(5, 8, "Testing for normality and other assumptions")
text(5, 7, Keyword)

pairs(Spooky_acomp, lower.panel = panel.qq) #QQplots and histos
title("Current Aitch. Composition",outer=TRUE,line=-1)
Spooky_ccomp 
title("Count Composition (ccomp)",outer=TRUE,line=-1)
plot(Spooky_centered_acomp)
title("Centered Aitch. Composition",outer=TRUE,line=-1)
pairs(Spooky_centered_acomp, lower.panel = panel.qq)
title("Centered Aitch. Composition",outer=TRUE,line=-1)

#### Step Six: Descriptive Statistics
#variance
mvar(Spooky_centered_acomp)

#metric standard deviation
msd(Spooky_centered_acomp)

#variation matrix
variation(Spooky_centered_acomp)

summary(Spooky_centered_acomp)
#mean.ratio is the geometric mean of the ratios and can be interpreted as a central value of each pairwise ratio

#PDF export code 
plot(0:10, type = "n", xaxt="n", yaxt="n", bty="n", xlab = "", ylab = "")
text(5, 8, "Descriptive Stats")
text(5, 7, Keyword)

#boxplot of pairwise ration
boxplot(Spooky_centered_acomp)
title("Boxplot of pairwise ration",outer=TRUE,line=-1)

#plotting ternary diagrans with the margin call
plot(Spooky_acomp, margin = "rcomp") #(CHECK) that these need to be non-centered
title("Boxplot of pairwise ration",outer=TRUE,line=-1)


#### Step Seven: Univariate Approaches (Treating the composition as one singular Y variable)
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

#plotting the composition and color-coding/shape-coding by the factors
opar = par(xpd=NA,no.readonly=TRUE)
plot(Y,pch=c(3,20)[X3],col=c("blue","green","red")[X1]) #pch picks the symbol, col picks the color, matched it to the 2 variable X3 for pch and X1 for col with three levels
legend(x=0.85,y=0.65,abbreviate(levels(X3), minlength=1), pch=c(3,20)) #legend for X3
legend(x=0.75,y=0.65,abbreviate(levels(X1), minlength=1), pch=20,col=c("blue","green","red"),yjust=0) #legend for X1
title("Composition by factor",outer=TRUE,line=-1)

boxplot(Y,X3,notch=TRUE) #boxplot of the data
title("Composition against ...",outer=TRUE,line=-1)

### ANOVAS and making the model
#For anovas, the book recommends using the ilr transformation, so we will do that
#the book also recommends setting the contrast to be the treatment level 
#(tbh i still dont know what contrasts are, but i trust the book)

#setting the contrasts
contrasts(X3) <- "contr.treatment" #using the vaccine variable for this

#making a model and getting its parameters (intercept, slope, mean, sigma (variance))
(single_X_model = lm(ilr(Y) ~ X3)) #model
(a = ilrInv(coef(single_X_model)[1,],orig=Y)) #backtransforming the intercept
(b = ilrInv(rbind(0,coef(single_X_model)[-1,]),orig=Y)) #backtransforming the slope? slash other constant
#the above output shows the increase between levels of the treatment variable
(mu = a + b) #mean
(Sigma=ilrvar2clr(var(single_X_model))) #variance

#plotting the means and the variance ellipses based on the model data
plot(mu,pch=20, col= c("blue","green"))
plot(Y,pch=".",add=TRUE)
ellipses(mu[1,],Sigma,2, col = "blue")
ellipses(mu[2,],Sigma,2, col = "green")
legend(x=0.75,y=0.65,abbreviate(levels(X3), minlength=1), pch=20,col=c("blue","green"),yjust=0)
title("Means and Variances from ilr(Y)~X",outer=TRUE,line=-1)

# !!! Remember to check your anova and linear model assumptions at this point 

anova(single_X_model) #actually running the anova

# Including mutliple predictors
(fullModel = lm(ilr(Y) ~ X1 + X2 + X3)) #model, the first level of each factor is assumed to be zero
coef(fullModel) #parameters
anova(fullModel)

# the above model output is all ilr transformed, so lets backtransform
(coefs <- ilrInv(coef(fullModel),orig=Y))
barplot(coefs,las=2,xlim=c(0,11)) #barplot of the coefficients (which I think is honestly not that useful)

#determining variances and plotting them as ellipses around the estimated parameters
vcov(fullModel) #variance-covariance matrix
vars <- vcovAcomp(fullModel)
dim(vars)
alpha=0.05


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
#  A variable may be removed from the model if the confidence ellipse around its parameter always contains the barycenter of the composition.

#the next suggestion is to test for model redundancy by making a bunch of differently ordered models and seeing if the last variable is still significant
model1 = lm(ilr(Y) ~ X3 + X2 + X1)
anova(model1)
model2 = lm(ilr(Y) ~ X3 + X1 + X2)
anova(model2)
model3 = lm(ilr(Y) ~ X2 + X3 + X1)
anova(model3)

#we can also look at multiplicative models
multiModel = lm(ilr(Y) ~ X1 * X2 * X3)
anova(multiModel)


## Checking the assumptions of the model (MODEL Diagnostics)
#This should occur before the model is run, and i'll move the code later when I review it again
#normally distributed residuals
Resid = ilrInv(resid(fullModel),orig=Y)
plot(Resid, cex=0.5)
title("Ternary Diagrams of Residuals",outer=TRUE,line=-1)

boxplot(Resid) #also reveals outliers
title("Boxplots of Residuals",outer=TRUE,line=-1)

#QQ plot (also normally distributed residuals)
qqnorm(Resid) 

#homoschedascity (equal variance)
Pred = ilrInv(predict(fullModel),orig=Y)
opar <- par(oma=c(3,3,0,0),mar=c(4,4,1,0))
pairwisePlot(clr(Pred),clr(Resid))
mtext(text=c("predicted values (clr)","residuals (clr)"),
      side=c(1,2),at=0.5,line=2,outer=TRUE)
par(opar)

#color-coordinate the above diagram by group
opar <- par(oma=c(0,0,2,0),mar=c(4,4,1,0))
pairwisePlot(clr(Pred),clr(Resid),col=as.numeric(X1))
legend(locator(1),levels(X1),col=as.numeric(1:3),pch=1,xpd=NA) #starts a click to place mechanism for the legend
par(opar)

#### Step Eight: Multivariate Approaches (treating the composition as a collection of multiple Ys)

#PDF code
plot(0:10, type = "n", xaxt="n", yaxt="n", bty="n", xlab = "", ylab = "")
text(5, 8, "Multivariate Approaches (PCA and LDA)")
text(5, 7, Keyword)

##PCA:
Spooky_PCA <- princomp(Spooky_centered_acomp) #This function returns a princomp object, containing the full result of a PCA on the covariance matrix of the clr-transformed datase, and requires an acomp object
Spooky_PCA #this is the PCA object

#plotting the scree plot
plot(Spooky_PCA) #this is the scree plot, showing the proportion of variance explained by each PC
#hideous, but informative

#plotting the biplot
sum(Spooky_PCA$sdev[1:2]^2)/mvar(Spooky_centered_acomp) #this is the proportion of variance explained by the biplot (first two PCs)
opar <- par(mar=c(1,1,1,1))
dots = rep(".",times=nrow(Spooky_centered_acomp))
biplot(Spooky_PCA, xlabs=dots)
par(opar)

##LDA (Linear Discriminant Analysis)
#defining the group of interest
Group = X1
mean(Spooky_acomp) #use non-centered data, because you dont care for normality here
table(Group) #okay if groups are not equal
res = lda( x=data.frame(ilr(Spooky_acomp)), grouping=Group )
res

#backtransforming from ilr coordinates to our data scales
ilrInv(res$means, orig=Spooky_acomp)

# Variance? Maybe?
V = ilrBase(Spooky_acomp) 
rownames(V)= colnames(Spooky_acomp)
t(ilr2clr(t(res$scaling), V=V)) #not 100% sure what this is showing me, but it is showing me something

#Showing the original groups and their alignment across the LDAs
grint = as.integer(Temp)
pairs(res, abbr=1, col=(1:4)[grint], cex=1.2)


#final PDF code
dev.off()