# Metabolimcs-practice
This is the first script used to test the metabolomics data on the aggression

# Firts attempt to open the test metabolomics files

# libaries

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
# lets inspect the column names

colnames(dataT2)

# Lets try a first PCA analysis for the aggresion status

aggrestatus <- dataT2$aggression_status
tune.pca(aggrestatus, center = TRUE, scale = TRUE)

# Does not work because the object has to be a numeric matrix and still have to add the metabolites
# let's first remove all the unimportant colums
# We want to be left with ID, sex, age_urine1 and aggression status
# these are columns 1,9, 353 and 367

dataT3 <- dataT2 %>%
select(-(2:8), -(10:352), -(354:366), -(368:373))
View(dataT3)
# I think we have removed to many, so let's try again
dataT2 <- dataT %>%
  select(-(2:8), -(10:352), -(354:366), -(368:373))
View(dataT2)
dataT3 <- na.omit(dataT2)
View(dataT3)

# Nopp, same result
# Ok, now lets get in the metabolomics data

metadata <- read.table(file = "C:/Users/Gebruiker/Documents/Master GBH/Stage BehPsy/Data/Metabolomics/20190808_QCedDNormalizedRINTTransformedValuesDF_NoCurium.tsv", sep = "\t", header = TRUE,)
View(metadata)
metadata2 <- na.omit(metadata)
View(metadata2)

#The subject names have an extra F, lets remove this

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

# Now let's get the aggression status in "0" for low and "1" for high aggression
colnames(comdata2)
comdata2$aggression_status.x <- gsub("low", "0", as.character(comdata2$aggression_status.x))
View(comdata2)
comdata2$aggression_status.x <- gsub("high", "1", as.character(comdata2$aggression_status.x))
View(comdata2)

# Lets now try again the PCA, Probably wont work, but lets anyway try

aggrestatus <- comdata2$aggression_status.x
tune.pca(aggrestatus, center = TRUE, scale = TRUE)

# Nopp, still need to make it a numeric matrix

?data.matrix
comdata3 <- data.matrix(comdata2, rownames.force = NA)
View(comdata3)

# Dont know if it has worked so lets try

colnames(comdata3)
aggrestatus <- comdata3["aggression_status.x"]
tune.pca(aggrestatus, center = TRUE, scale = TRUE)

# Still needs to be in numeric matrix, so did not work
# Lets try somethings

comdata4 <- data.matrix(frame = comdata2, rownames.force = NA)
View(comdata4)
aggrestatus <- comdata4["aggression_status.x"]
tune.pca(aggrestatus, center = TRUE, scale = TRUE)

comdata5 <- data.matrix(frame = comdata2)
View(comdata5)
aggrestatus <- comdata5["aggression_status.x"]
tune.pca(aggrestatus, center = TRUE, scale = TRUE)

?matrix
comdata6 <- matrix(data = comdata2, nrow = 175, ncol = 103, byrow = FALSE, dimnames = NULL)
View(comdata6)
aggrestatus <- comdata6[[3]]
tune.pca(aggrestatus, center = TRUE, scale = TRUE)
