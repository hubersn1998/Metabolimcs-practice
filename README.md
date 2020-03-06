# Metabolimcs-practice
# This is the first script used to test the metabolomics data on the aggression

# Firts attempt to open the test metabolomics files

# Firts attempt to open the test metabolomics files

#libaries

library(foreign)
library(mixOmics)
library(dplyr)
library(base)

# Get the data

setwd("C:/Users/Gebruiker/Documents/Master GBH/Stage BehPsy/Data/Metabolomics")
dataset = read.spss("ACTION_NTR_atcollection.sav", to.data.frame = TRUE)
View(dataset)
datasetdemo = read.spss("ACTION_NTR_demogr.sav", to.data.frame = TRUE)
View(datasetdemo)

# Not possible to merge them, because to big, so first remove the missings.

dataset2 <- na.omit(dataset)
View(dataset2)

# In demo there are two columns without any variables, so lets remove them

datasetdemo[16:17] <- list(NULL)
head(datasetdemo)
datasetdemo2 <- na.omit(datasetdemo)
View(datasetdemo2)

# Let's check which names to use to merge the data on

colnames(dataset2)
colnames(datasetdemo2)

# Let's now merge it

dataT <- merge(x=dataset2, y=datasetdemo2, by.x = "SubjectFISNumber", by.y = "FISNumber", all.y = TRUE)
View(dataT)

# Looks better, but now there are agains missing data, which cannot be used in mixOmics
# So let's remove these

dataT2 <- na.omit(dataT)
View(dataT2)

# Ok, now we have 175 samples left, with 373 columns.
# lets inspect the column names

colnames(dataT2)


# let's first remove all the unimportant colums
# We want to be left with ID, sex, age_urine1 and aggression status
# these are columns 1,9, 353 and 367

dataT3 <- dataT2 %>%
  select(-(2:8), -(10:352), -(354:366), -(368:373))
View(dataT3)

# Ok, now lets get in the metabolomics data

metadata <- read.table(file = "C:/Users/Gebruiker/Documents/Master GBH/Stage BehPsy/Data/Metabolomics/20190808_QCedDNormalizedRINTTransformedValuesDF_NoCurium.tsv", sep = "\t", header = TRUE,)
View(metadata)
metadata2 <- na.omit(metadata)
View(metadata2)

# The subject names have an extra F, lets remove this

colnames(metadata2)
metadata2$SubjectFISNumber <- gsub("F", "", as.character(metadata2$SubjectFISNumber))
View(metadata2)

# Now we can merge the metadata2 and the dataT3

colnames(metadata2)
colnames(dataT3)
comdata <- merge(x = dataT3, y = metadata2, by.x = "SubjectFISNumber", by.y = "SubjectFISNumber", all.y = TRUE )
View(comdata)

# Now we have to remove the missings again

comdata2 <- na.omit(comdata)
View(comdata2)


# Lets now try the PCA, Probably wont work, but lets anyway try

aggrestatus <- as.numeric(comdata2$aggression_status.x)
tune.pca(aggrestatus, center = TRUE, scale = TRUE)

# Nopp, need to make it a seperate numeric matrix

comdatameta <- comdata2 %>%
  select(-(1:4), -(89:103))
View(comdatameta)

# Now lets make it into a data frame

# Not necesarry, you can just use the data with only the metabolomics data

tune.pca(comdatameta, center = TRUE, scale = TRUE)

# It works!
# So let's do a pca now based on the first three components

pca.metabol <- pca(comdatameta, ncomp = 5, center = TRUE)

# And let's make some plots
plotVar(pca.metabol)

# Individuals plot : same plots, different colors
plotIndiv(pca.metabol)
plotIndiv(pca.metabol, group = comdata2$MedicationUse, legend=TRUE)
plotIndiv(pca.metabol, group = comdata2$aggression_status.x, legend=TRUE)
plotIndiv(pca.metabol, ind.names = comdata2$aggression_status.x, group = comdata2$sex.x, legend=TRUE)
plotIndiv(pca.metabol, group = comdata2$aggression_status.x, ind.names = comdata2$sex.x, legend=TRUE,col=color.mixo(4:5))

View(comdata2)

# Seems like i have removed all of the males
# Let's find out why
