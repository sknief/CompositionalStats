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
