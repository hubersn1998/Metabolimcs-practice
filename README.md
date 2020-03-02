# Metabolimcs-practice
This is the first script used to test the metabolomics data on the aggression

# Firts attempt to open the test metabolomics files

#libaries

library(foreign)
library(mixOmics)
library(dplyr)

# Get the data

setwd("C:/Users/Gebruiker/Documents/Master GBH/Stage BehPsy/Data/Metabolomics")
dataset = read.spss("ACTION_NTR_atcollection.sav", to.data.frame = TRUE)
View(dataset)
datasetdemo = read.spss("ACTION_NTR_demogr.sav", to.data.frame = TRUE)
View(datasetdemo)

# Not possible to merge them, because to big, so first remove the missings.

?complete.cases()
complete.cases(dataset)
dataset2 <- na.omit(dataset)
View(dataset2)

# In demo there are two columns without any variables, so lets remove them

datasetdemo[16:17] <- list(NULL)
head(datasetdemo)
datasetdemo2 <- na.omit(datasetdemo)
View(datasetdemo2)

# Ok, now lets see if we can merge  them

merge(dataset2, datasetdemo2)
dataT <- merge(dataset2, datasetdemo2)
View(dataT)

# much more rows were added, so lets try something else

?merge

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
