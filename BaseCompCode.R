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


#### Step 2: Data Import
#$ [Insert Data Here]
#$ [Optional: Data Cleaning]

# turn the composition into an acomp object
xxca <- acomp(xxc)
xxca

#explore the composition
plot(xxca)
plot(clr(xxca)) #clr transformation (central log)
plot(ilr(xxca)) #ilr transformation (isometric log) (check)


#### Step Three: Checking for and dealing with missing values
# notice: you need to work with an acomp object for the following code
missingSummary(xxca) #table of missing values
# we will use the imputation method of the book in order to deal with missing values, which is to replace them with the smallest non-zero value in the composition (NOT OG DATASET) 
BDL_MNAR_REPLACE <- 0.001 #assign a value to replace BDL and MNARs
xxca <- zeroreplace(xxca, BDL_MNAR_REPLACE) #replace zeros
missingSummary(xxca) #check for missing values again
xxca[is.na(xxca)] <- BDL_MNAR_REPLACE #replace NAs
missingSummary(xxca) # final check

#explore the composition again
plot(xxca)
plot(clr(xxca)) #clr transformation (central log)
plot(ilr(xxca)) #ilr transformation (isometric log) (check)


#### Step Four: Setting up exports and the whole code
#All variables, output and models start with the term "Spooky" - you need to ctrl + f and replace this term with a keyword relevant to your current research question

#Let's save both of these pieces of information
researchQuestion <- "Enter Research Q" #keep brief
keyword <- "Spooky" #Spooky by default


#we will use the keyword in the export pdf titles as well
PDF_Ident <- paste0("BaseCompAnalysis_", keyword, "_") #this is to identify the output from your specific run or hypothesis
pdf(file = paste(PDF_Ident, "1.pdf", sep = "_"))

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

#we can also note our research question and keyword
plot.new()
text(x=.5, y=1, keyword, font=3, cex=2, col="blue")
text(x=.5, y=.7, researchQuestion)  

dev.off() #this closes the pdf file (shuts off the device)
#going forward, most plots will have the export to pdf text at the top and bottom of the relevant section, but commented out so that you can visualize plots in R first


### Step Five: Testing for normality + other assumptions (DATA diagnostics)
# First, testing for normality (normal distribution)
mvnorm.etest(ilr(xc)) #complete multivariate test
acompNormalGOF.test(xc, R = 311) #alternative multivariate test, works on non-transformed data

#exploring through graphs
#better panel making code from stackoverflow
panel.qq <- function(x, y, ...) {
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1), new = TRUE)
  qqplot(x, y, xlab = deparse(substitute(x)), ylab = deparse(substitute(y)))
  abline(c(0,1), ...)
}

#QQplots
pairs(xxca, lower.panel = panel.qq) #use acomp object for these


#Alternative tests to check if a count comp distribution fits better 
xxcc <- ccomp(xxc)
xxcc
plot(xxcc) #explore the composition

#Alternative alternatives: multinomial and multi-poisson
ccompMultinomialGOF.test(xxcc, R = 310) #not multinomially distributed, no use in using ccomp
ccompPoissonGOF.test(xxcc) # not a multi-Poission distribution

#If composition failed normality test, we can try to transform or center
#centering the data
mean(xxca)
cxxca <- xxca - mean(xxca) #centering by subtracting the mean
mean(cxxca) #check if the centered mean is the neutral element

#examining centered data
plot(cxxca)
pairs(cxxca, lower.panel = panel.qq) #looks better
mvnorm.etest(cxxca, R = 310) #complete multivariate test, 

#according to the book, scaling is often used to compare subcompositions, but centering is used when date is often squashed into one corner of the simplex, which is what we had 

#### Step Six: Descriptive Statistics
#variance
mvar(cxxca)

#metric standard deviation
msd(cxxca)

#variation matrix
variation(cxxca)

summary(cxxca)
#mean.ratio is the geometric mean of the ratios and can be interpreted as a central value of each pairwise ratio

#boxplot of pairwise ration
boxplot(cxxca)

#plotting ternary diagrans with the margin call
plot(xxca, margin = "rcomp") 


#### Step Seven: Univariate Approaches (Treating the composition as one singular Y variable)
#Defining the independent and dependant variables
Y = cxxca #centered and cleaned composition from above
X1 = BD100s$TempTreat
X2 = BD100s$Timepoint
X3 = BD100s$Treatment

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
plot(Y,pch=c(3,20)[X3],col=c("blue","green","red")[X1]) #pch picks the symbol, col picks the color, matched it to the 2 variable X3 for pch and X1 for col with three levels
legend(x=0.85,y=0.65,abbreviate(levels(X3), minlength=1), pch=c(3,20)) #legend for X3
legend(x=0.75,y=0.65,abbreviate(levels(X1), minlength=1), pch=20,col=c("blue","green","red"),yjust=0) #legend for X1

boxplot(Y,X3,notch=TRUE) #boxplot of the data


### ANOVAS and making the model
#For anovas, the book recommends using the ilr transformation, so we will do that
#the book also recommends setting the contrast to be the treatment level 
#(tbh i still dont know what contrasts are, but i trust the book)

#setting the contrasts
contrasts(X3) <- "contr.treatment" #using the vaccine variable for this

#making a model and getting its parameters (intercept, slope, mean, sigma (variance))
(model = lm(ilr(Y) ~ X3)) #model
(a = ilrInv(coef(model)[1,],orig=Y)) #backtransforming the intercept
(b = ilrInv(rbind(0,coef(model)[-1,]),orig=Y)) #backtransforming the slope? slash other constant
#the above output shows the increase between levels of the treatment variable
(mu = a + b) #mean
(Sigma=ilrvar2clr(var(model))) #variance

#plotting the means and the variance ellipses based on the model data
plot(mu,pch=20, col= c("blue","green"))
plot(Y,pch=".",add=TRUE)
ellipses(mu[1,],Sigma,2, col = "blue")
ellipses(mu[2,],Sigma,2, col = "green")
legend(x=0.75,y=0.65,abbreviate(levels(X3), minlength=1), pch=20,col=c("blue","green"),yjust=0)

# !!! Remember to check your anova and linear model assumptions at this point 

anova(model) #actually running the anova

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
model4 = lm(ilr(Y) ~ X1 * X2 * X3)
anova(model4)


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

##PCA:
pccxxca <- princomp(cxxca) #This function returns a princomp object, containing the full result of a PCA on the covariance matrix of the clr-transformed datase, and requires an acomp object
pccxxca #this is the PCA object

#plotting the scree plot
plot(pccxxca) #this is the scree plot, showing the proportion of variance explained by each PC
#hideous, but informative

#plotting the biplot
sum(pccxxca$sdev[1:2]^2)/mvar(cxxca) #this is the proportion of variance explained by the biplot (first two PCs)
opar <- par(mar=c(1,1,1,1))
dots = rep(".",times=nrow(cxxca))
biplot(pccxxca, xlabs=dots)
par(opar)

##LDA (Linear Discriminant Analysis)
#defining the group of interest
Group = X1
mean(xxca) #use non-centered data
table(Temp) #okay if groups are not equal
res = lda( x=data.frame(ilr(xxca)), grouping=Group )
res

#backtransforming from ilr coordinates to our data scales
ilrInv(res$means, orig=xxca)

# Variance? Maybe?
V = ilrBase(xxca) 
rownames(V)= colnames(xxca)
t(ilr2clr(t(res$scaling), V=V)) #not 100% sure what this is showing me, but it is showing me something

#Showing the original groups and their alignment across the LDAs
grint = as.integer(Temp)
pairs(res, abbr=1, col=(1:4)[grint], cex=1.2)
