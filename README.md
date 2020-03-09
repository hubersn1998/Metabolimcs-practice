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
dataT$sex <- gsub("male", "0", as.character(dataT$sex))
dataT$sex <- gsub("female", "1", as.character(dataT$sex))
View(dataT)
# First do female, otherwise the male already gets replaced
dataT2 <- na.omit(dataT3)
View(dataT2)

# Ok, now we have 175 samples left, with 373 columns.
# This is the step when all the men get deleted. Let's see why
# lets inspect the column names
?na.omit
colnames(dataT2)


# let's first remove all the unimportant colums
# We want to be left with ID, sex, age_urine1 and aggression status
# these are columns 1,9, 353 and 367

dataT3 <- dataT %>%
  select(-(2:8), -(10:352), -(354:366), -(368:373))
View(dataT3)
# turns out we only have the aggrassion status for females,
# This is why all the man get deleted. 

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
comdata <- merge(x = dataT2, y = metadata2, by.x = "SubjectFISNumber", by.y = "SubjectFISNumber", all.y = TRUE )
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
?plotVar
# Individuals plot : same plots, different colors
plotIndiv(pca.metabol)
plotIndiv(pca.metabol, group = comdata2$aggression_status.x, legend=TRUE)
plotIndiv(pca.metabol, group = comdata2$aggression_status.x, legend=TRUE)
plotIndiv(pca.metabol, ind.names = comdata2$aggression_status.x, group = comdata2$sex.x, legend=TRUE)
plotIndiv(pca.metabol, group = comdata2$aggression_status.x, ind.names = comdata2$sex.x, legend=TRUE,col=color.mixo(4:5))

View(comdata2)

# Seems like i have removed all of the males
# Let's find out why --> it's because of the aggression status

# let's try a pls
?pls
pls.metabol <- pls(comdatameta, comdata2$aggression_status.y, ncomp = 7) 
pls.metabol
plotVar(pls.metabol)
plotIndiv(pls.metabol)
setwd("C:/Users/Gebruiker/Documents/Master GBH/Stage BehPsy/Data/Metabolomics")
bmp("cim_plot_meta")
cim(pls.metabol)
dev.off()
?dev.off

# The cim plot does nt yet work, rest does, still have to look up how this would work
# Now let's make seperate files for the different metabolic subsets
View(comdatameta)
# We have the amines, the lipo acids and the biomarkers
# Let's maka a seperate one for amines and lipoacids
# once we will get the steroïd, we can also split those

amines <- comdatameta %>%
  select(1:64)
View(amines)
lipids <- comdatameta %>%
  select(65:84)
View(lipids)

# Now we can perform pcs and pls on the differnt types of lipid seperatly
tune.pca(amines, center = TRUE, scale = TRUE)
pca.amines <- pca(amines, ncomp = 2, center = TRUE)
View(pca.amines)
plotVar(pca.amines, var.names = FALSE)
plotIndiv(pca.amines)
plotIndiv(pca.amines, group = comdata2$aggression_status.x, legend = TRUE)

tune.pca(lipids, center = TRUE, scale = TRUE)
pca.lipids <- pca(lipids, ncomp = 2, center = TRUE)
View(pca.lipids)
plotVar(pca.lipids)
plotIndiv(pca.lipids, group = comdata2$aggression_status.x, legend = TRUE)

# let's try to combine them in a pls
pls.Amines_lipids <- pls(amines, lipids, ncomp = 2)
plotVar(pls.Amines_lipids)
plotIndiv(pls.Amines_lipids, rep.space = "XY-variate", group = comdata2$aggression_status.x, legend = TRUE)
plotLoadings(pls.Amines_lipids, comp=1)
plotLoadings(pls.Amines_lipids, comp=2)

# All fun and all, but now we must try to get the residuals 
# And then perform these analysis

?glm
# Lets make a dataframe, with the amines, age, agression and ID
# Sex does not have to be included, because we only have females
amines.reg <- comdata2 %>%
  select((2:3), (5:68))
View(amines.reg)
colnames(amines.reg)
amines.reg$age_at_collection.x <- as.numeric(amines.reg$age_at_collection.x)

model.amines <- glm(aggression_status.x ~., family = binomial(link = 'logit'), data=amines.reg)
summary(model.amines)
str(amines.reg)
# Let's see how the residuals are shown and what we can do with them

tune.pca(amines.reg, center = TRUE, scale = TRUE)
df.amines <- data.frame(matrix(unlist(model.amines["coefficients"]), nrow = 66, byrow = T))
tune.pca(df.amines, center = TRUE, scale = TRUE)
?tune.pca
?pca
?pls

# Let's try something else
# Let's make the PCA's and then use the residuals from those.
# Is that even what we have to do?
# Like do we have to use the residuals of the pcs/pls analysis or of the logistic regression
# If it is the logistic regression, then we have to figure out how to get them back into the matrix
# Let's ask Fiona tomorrow, because today she is very busy
# For today we will figure out how the multi thing works
