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
# For this first remove the unnesecary columns, so we do not delete to many persons
dataset2 <- dataset %>%
  select(1, 4, 9, (352:353))
View(dataset2)
dataset3 <- na.omit(dataset2)
View(dataset3)

# In demo there are two columns without any variables, so lets remove them

datasetdemo2 <- datasetdemo %>%
  select(-3, -14, -(16:17), -(18:20), -(21:23))
View(datasetdemo2)
datasetdemo3 <- na.omit(datasetdemo2)
View(datasetdemo3)

# Let's check which names to use to merge the data on

colnames(dataset2)
colnames(datasetdemo2)

# Let's now merge it

dataT <- merge(x=dataset3, y=datasetdemo3, by.x = "SubjectFISNumber", by.y = "FISNumber", all.y = TRUE)
View(dataT)
dataT2 <- dataT %>%
  select(1, 5, 8, 17)
View(dataT2)
# Looks better, but now there are agains missing data, which cannot be used in mixOmics
# So let's remove these

# First do female, otherwise the male already gets replaced
dataT3 <- na.omit(dataT2)
View(dataT3)

# Ok, now we have 1245 samples left, with 4 columns.
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

# Now we have 1089 samples left with 103 columns

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
plotIndiv(pca.metabol, group = comdata2$sex.x, legend=TRUE)
plotIndiv(pca.metabol, group = comdata2$aggression_status.x, legend=TRUE)
plotIndiv(pca.metabol, ind.names = comdata2$aggression_status.x, group = comdata2$sex.x, legend=TRUE)
plotIndiv(pca.metabol, group = comdata2$aggression_status.x, ind.names = comdata2$sex.x, legend=TRUE,col=color.mixo(4:5))

View(comdata2)

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

?resid
datares <-glm(aggression_status.x ~ SES +sex.x, data = comdata2$ACTION_amines.1.Methylhistidine, family = "binomial")
View(data.res)
resid(datares)
str(comdata2)
comdata3 <- resid(datares)
View(comdata3)

# Look on Thursday furter into the residuals etc. 
