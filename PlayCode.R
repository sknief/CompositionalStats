### NOTES: This file is me going through the book and its exercises on the x100s data, irrespective of treatment and timepoints. I understand that I will likely have to repeat this and reorganize this for my actual research questions to be answered at a later timepoint. 


#### prereqs #####
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

### Basics ######
# acomp tells R to consider the its argument as a compositional dataset, forcing closure to 1 
# a composition which elements sum up to a constant is called a closed composition, this is what we are working with, we can close a vector or dataset using the command clo
# the set of possible closed compositions is called the D-part simplex
# pertubation is sum/subtraction or translation. command perturbe(x,y)
# powering is the product of components, multiplication, power(x,y)
# centered log-ratio transformation (CLR), clr() and clrInv()
#isometric log-ratio transformation (IRL), ilr(), and ilrInv() to undo,


#clean datafile
BloodData <- Bloodsmear_data_full[1:13]
BD100s <- BloodData[c(1:6, 11:13)]
colnames(BD100s) <- c('TempTreat', 'Timepoint', 'Subgroup', 'SlideID', 'FishID', 'Treatment', 'Granu', 'Lymph', 'Mono')  
BD1000s <- BloodData[1:10]
colnames(BD1000s) <- c('TempTreat', 'Timepoint', 'Subgroup', 'SlideID', 'FishID', 'Treatment', 'RBC', 'Granu', 'Lymph', 'Mono')  


# tester to see if compositions is working
xc = acomp(BD100s, c('Granu', 'Lymph', 'Mono'))
xc
plot(xc) #it is

plot(ilr(xc)) #isometric log-ratio transformation
xc

## Picking an appropriate scale
# We will use the Aitchison compositional scale (acomp), rather than the count compositional scale, due to the following arguments
# 1- the stats textbook I am reading recommends using acomp for count data, and 
# 2 - the underlying true (population) composition of the count data is likely an Aitchisons composition
# 3 - I believe the total of these compositions to be large enough to avoid random selection error in counts, thus avoiding the issues raised at 2.4.7 of the book

## Testing for compositional normality
# because the additive lognormal distribution is the typical reference for all linear models, time to work in coordinates (transforming datasets)
mvnorm.etest(ilr(xc)) #complete multivariate test
#alternative multivariate test below
acompNormalGOF.test(xc, R = 311) # not ALN-distributed

#trying to find the source of the nonnormality

#....

#troubleshooting (reset button)
xxc <- as.data.frame(BD100s[,c('Granu', 'Lymph', 'Mono')])
xxca <- acomp(xxc)
xxca
######

#$$$$ Intermission $$$$ 

##### dealing with missing values #####
missingSummary(xxca)
min(xxca, na.rm = TRUE) #great, thats just zero
BDL_MNAR_REPLACE <- 0.001
# a bunch of below detection limits (BDL) and one Missing Not at Random
# we will use the imputation method of the book, which is to replace the BDL with the smallest non-zero value in the dataset (1), 
xxca <- zeroreplace(xxca, BDL_MNAR_REPLACE)
missingSummary(xxca)
xxca #this left me with a bunch of MNARs which are NAs
# we will replace the MNARs
MNARvalue
xxca[is.na(xxca)] <- BDL_MNAR_REPLACE
xxca
missingSummary(xxca) # we good now
# we will also remove the rogue MNAR value
xxca <- xxca[-c(56),]
xxca

## i fucked up: i should have been using 0.1 for the replacements as i was already dealing with a closed composition, not the raw counts
# i will redo the above with 0.1 (fixed now)

#lets try the whole shebang again


### EXAMINING THE DATA AND TESTING FOR NORMALITY #####

plot(xxca)
plot(clr(xxca))
plot(ilr(xxca))
mvnorm.etest(xxca, R = 310) #complete multivariate test, still not normal (test should be acomp and not transformed)
qqnorm(xxca[, 'Granu'])

#better panel making code from stackoverflow
panel.qq <- function(x, y, ...) {
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1), new = TRUE)
  qqplot(x, y, xlab = deparse(substitute(x)), ylab = deparse(substitute(y)))
  abline(c(0,1), ...)
}

#QQplots
pairs(xxca, lower.panel = panel.qq) #missing values error, fixed by using the non-acomp version, but i think i was meant to be doing this on ile transformed a comp? 
#FIXED BY FIXING THE NON MISSING VALUES
pairs(xxc, lower.panel = panel.qq) #non-acomp version for comparision
pairs(ilr(xxca), lower.panel = panel.qq) #ilr version for comparision, looks uh.... weird.
#the lymphocytes are clearly the issue here

#testing for other distributions
acompDirichletGOF.test(acomp(xxca),R=310) #seems this function has been removed since

#testing for count comp distributions
xxcc <- as.data.frame(BD100s[,c('Granu', 'Lymph', 'Mono')])
xxcc <- ccomp(xxc)
xxcc
plot(xxcc)
xxcc <- xxcc[-c(56),]
xxcc

ccompMultinomialGOF.test(xxcc, R = 310) #not multinomially distributed, no use in using ccomp
 
ccompPoissonGOF.test(xxcc) # not a multi-Poission distribution

#next logical step is to center the data and try again? 
#centering the data
mean(xxca)
cxxca <- xxca - mean(xxca) #centering by subtracting the mean
mean(cxxca) #check if the centered mean is the neutral element (OMG IT IS!)

#examining centered data
plot(cxxca)
pairs(cxxca, lower.panel = panel.qq) #looks better
mvnorm.etest(cxxca, R = 310) #complete multivariate test, still not normal (test should be acomp and not transformed) 
#hopefully once i seperate out my date by treatment and timepoint, it will be normal
# from the qqnorm plots though it looks normal enough for what I want to do tbh

#according to the book, scaling is often used to compare subcompositions, but centering is used when date is often squashed into one corner of the simplex, which is what we had 

#### DESCRIPTIVE STATS #####
#from here on everything needs to be an acomp object, i think

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

#I'm gonna skip the sections on predictions and confidence regions as I dont think they will be relevant to me 

#plotting ternary diagrans with the margin call
plot(xxca, margin = "rcomp") #doesnt do anything cause i only have 3 variables but will be useful for the x1000s

###### LINEAR MODEL TIME #######

#I will only be following the sections on compositions as dependant variables, for obvious reasons

#Defining the independent and dependant variables
Y = cxxca #centered and cleaned composition from above
X1 = BD100s$TempTreat
X2 = BD100s$Timepoint
X3 = BD100s$Treatment

#for all of these, I need to remove the same line of data as I did above: 
X1 <- X1[-c(56)]
X2 <- X2[-c(56)]
X3 <- X3[-c(56)]

#fix an issue with X2
X2 <- replace(X2, X2 == 17, 7) #okay that worked

#I will also need to convert these to factors (and check their orders)
X1 <- as.factor(X1)
X2 <- as.factor(X2)
X3 <- as.factor(X3)

#checking the levels orders
levels(X1)
levels(X2) 
levels(X3)
#all good

#End for today, at 5.3.1.3 in the book

#Back to it on the 29th
#all of my predictors are factors, we we'll use colorcoding for the plots
opar = par(xpd=NA,no.readonly=TRUE)
plot(Y,pch=c(3,20)[X3],col=c("blue","green","red")[X1]) #pch picks the symbol, col picks the color, matched it to the 2 variable X3 for pch and X1 for col with three levels
legend(x=0.85,y=0.65,abbreviate(levels(X3), minlength=1), pch=c(3,20)) #legend for X3
legend(x=0.75,y=0.65,abbreviate(levels(X1), minlength=1), pch=20,col=c("blue","green","red"),yjust=0) #legend for X1

#plot is obviously not a set of a bunch of triangles, only one cause we only have 3 variables, but it works well 

#trying a different plot
boxplot(Y,X3,notch=TRUE) #boxplot of the data

#just out curiosity, lets see do the first plot again with time instead of vax
plot(Y,pch=c(3,20, 8)[X1],col=c("blue","green","red", "black")[X2]) #pch picks the symbol, col picks the color, matched it to the 2 variable X3 for pch and X1 for col with three levels
legend(x=0.85,y=0.65,abbreviate(levels(X1), minlength=1), pch=c(3,20, 8)) #legend for X3
legend(x=0.75,y=0.65,abbreviate(levels(X2), minlength=1), pch=20,col=c("blue","green","red", "black"),yjust=0)
#looks awful 

## Regression time!! OR rather, ANOVA TIME

#For anovas, the book recommends using the ilr transformation, so we will do that
#the book also recommends setting the contrast to be the treatment level 
#tbh i still dont know what contrasts are, but i trust the book

#setting the contrasts
contrasts(X3) <- "contr.treatment"#using the vaccine variable for this

#making a model and getting its parameters (intercept, slope, mean, sigma (variance))
(model = lm(ilr(Y) ~ X3)) #model
(a = ilrInv(coef(model)[1,],orig=Y)) #backtransforming the intercept
(b = ilrInv(rbind(0,coef(model)[-1,]),orig=Y)) #backtransforming the slope? slash other constant
#the above output shows the increase between levels of the treatment variable
(mu = a + b)
(Sigma=ilrvar2clr(var(model)))

#plotting the means and the ellipses based on the model data
plot(mu,pch=20, col= c("blue","green"))
plot(Y,pch=".",add=TRUE)
ellipses(mu[1,],Sigma,2, col = "blue")
ellipses(mu[2,],Sigma,2, col = "green")
legend(x=0.75,y=0.65,abbreviate(levels(X3), minlength=1), pch=20,col=c("blue","green"),yjust=0)
#this is a plot of the means of the different treatments, with the ellipses showing the 95% confidence intervals

#actually running a proper ANOVA
anova(model) #highly significant??? even though I would not have guessed that from the above graph

##Full Model
#time to put in all predictor variables - but I admit i'm not sure if this model is the most useful
(fullModel = lm(ilr(Y) ~ X1 + X2 + X3)) #model 
#the first level of each factor is assumed to be zero
#the same output can also be gotten with
coef(fullModel)
#out of curiosity, lets see the anova
anova(fullModel)
# ALL SIGNIFICANT??????

#the above model output is all ilr transformed, so lets backtransform
(coefs <- ilrInv(coef(fullModel),orig=Y))
barplot(coefs,las=2,xlim=c(0,11)) #barplot of the coefficients, which I think is honestly not that useful

#if we display the uncertainty of these parameters, we need to estimate their variance
vcov(fullModel) #variance-covariance matrix, its kinda big
#then we can plot the variance and uncertainty
vars <- vcovAcomp(fullModel)
dim(vars)
alpha=0.05
plot(coefs,pch=as.character(1:nrow(coefs)))
plot(coefs-coefs,pch=20,add=TRUE)
for(i in 1:nrow(coefs)){
  ellipses(acomp(coefs[i,,drop=FALSE]),
  ilrvar2clr(vars[,,i,i]),
  r=ConfRadius(model,1-alpha),
  col="gray")
}
par(xpd=FALSE)
#  A variable may be removed from the model if the confidence ellipse around its parameter always contains the barycenter of the composition. - so according to this graph, X1Fluc is not different from the general mean, but Opt and Cold must be different from each other at least (given the significant result earlier)

#another plot option: the biplot
#from the book: Each arrow represents a one-unit change in the direction of each explanatory variable. Projecting these arrows onto the links between the compositional variables, we obtain the change in the response composition predicted by the linear model

#I cannot get this code to work, so I will skip it for now
B = clr(coefs[-1,])
svdlm = svd(B)

# Verify dimensions
print(dim(svdlm$u))  # Should be 6x3
print(length(svdlm$d))  # Should be 3
print(dim(svdlm$v))  # Should be 3x3

opar <- par(xpd=NA, mar=c(2,2,0,0))

# Assuming coloredBiplot expects two matrices of the same number of columns
coloredBiplot(svdlm$u, y = (svdlm$v %*% diag(svdlm$d)), scale = 0, xarrows = TRUE, ypch = 4, ynames = colnames(B), xnames = rownames(B))

sum(svdlm$d[1:2]^2)/sum(svdlm$d^2)

# Reset the plot parameters to original
par(opar)

## end skip 


#the next suggestion is to test for model redundancy by making a bunch of differently ordered models and seeing if the last variable is still significant
model1 = lm(ilr(Y) ~ X3 + X2 + X1)
anova(model1)
model2 = lm(ilr(Y) ~ X3 + X1 + X2)
anova(model2)
model3 = lm(ilr(Y) ~ X2 + X3 + X1)
anova(model3)
#all looks good. didnt do any multiplicative models but thats a worry for the legit runs

##More diagnostics
#residuals
Resid = ilrInv(resid(fullModel),orig=Y)
plot(Resid, cex=0.5)
title("Ternary Diagrams of Residuals",outer=TRUE,line=-1)

boxplot(Resid) #decent amount of outliers
title("Boxplots of Residuals",outer=TRUE,line=-1)

#QQ plot
qqnorm(Resid) #WOAH looking weird

#we already dia whole bunch of diagnostics earlier though, so not quite sure why its popping up now again

#homoschedascity
Pred = ilrInv(predict(fullModel),orig=Y)

opar <- par(oma=c(3,3,0,0),mar=c(4,4,1,0))
pairwisePlot(clr(Pred),clr(Resid))
mtext(text=c("predicted values (clr)","residuals (clr)"),
        side=c(1,2),at=0.5,line=2,outer=TRUE)
par(opar)
#honestly pretty ok?

#color-coordinate by group
opar <- par(oma=c(0,0,2,0),mar=c(4,4,1,0))
pairwisePlot(clr(Pred),clr(Resid),col=as.numeric(X1))
legend(locator(1),levels(X1),col=as.numeric(1:3),pch=1,xpd=NA)
par(opar)
# ehhh still ok

#last step of model optimisation: testing an interactive model
#i will not go thru the whole shebang but ill make one base model
model4 = lm(ilr(Y) ~ X1 * X2 * X3)
anova(model4)
#my laptop has been thinking about this model for a looong time.... not good.

#We are at the end of Ch 5 now, so multivariate stuff will happen tmrw
#once i am done with this script, I'll make a cleaner base / instructional one about the process, and then the relevant ones for my research question 
#add an export to pdf function like Liss had been using
#should be able to finish this tmrw? 
