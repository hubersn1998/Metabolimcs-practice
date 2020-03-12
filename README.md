# Metabolimcs-practice
# This is the first script used to test the metabolomics data on the aggression

# Firts attempt to open the test metabolomics files

# libaries

library(foreign)
library(mixOmics)
library(dplyr)
library(base)

# Get the aggression data

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

colnames(dataset3)
colnames(datasetdemo3)

# Let's now merge it

dataT <- merge(x=dataset3, y=datasetdemo3, by.x = "SubjectFISNumber", by.y = "FISNumber", all.y = TRUE)
View(dataT)
dataT2 <- dataT %>%
  select(1, 3, 5, 8, 17)
View(dataT2)

# Looks better, but now there are agains missing data, which cannot be used in mixOmics
# So let's remove these

dataT3 <- na.omit(dataT2)
View(dataT3)

# Ok, now we have 1245 samples left, with 5 columns.
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

# Now we have 1089 samples left with 104 columns
# let's make different data frames with only the metabolimcs data

comdatameta <- comdata2 %>%
  select(-(1:5), -(90:104))
View(comdatameta)

# Now let's make seperate files for the different metabolic subsets
# We have the amines, the lipo acids and the biomarkers
# Let's maka a seperate one for amines and lipoacids
# once we will get the steroïd, we can also split those

amines <- comdatameta %>%
  select(1:64)
View(amines)
lipids <- comdatameta %>%
  select(65:84)
View(lipids)

# Now we must try to get the residuals 
# And then perform these analysis

for (i in 1:ncol(amines)) {
  amines[,i] <- resid(lm(amines[,i] ~ comdata2$sex.x + comdata2$age_at_collection.x))
}
View(amines)

for (i in 1:ncol(lipids)) {
  lipids[,i] <- resid(lm(lipids[,i] ~ comdata2$sex.x + comdata2$aggression_status.x))
}
View(lipids)

# Ok this works
# Now lets try to run the pcs, pls and multiblock analysis

# Amines
tune.pca(amines)
pca.amines <- pca(amines, ncomp=2)
setwd("C:/Users/Gebruiker/Documents/Master GBH/Stage BehPsy/Voorbereiding/Test Figures")
plotVar(pca.amines)
dev.off()
plotIndiv(pca.amines, group = comdata2$aggression_status.x, legend = TRUE)

# Lipids
tune.pca(lipids)
pca.lipids <- pca(lipids, ncomp=3)
plotVar(pca.lipids)
plotIndiv(pca.lipids, group = comdata2$aggression_status.x)

# pls together
pls.meta <- pls(amines, lipids, ncomp = 10)
plotVar(pls.meta, var.names = FALSE)
plotIndiv(pls.meta, rep.space = "XY-variate", group = comdata2$aggression_status.x, legend = TRUE)
jpeg("C:/Users/Gebruiker/Documents/Master GBH/Stage BehPsy/Voorbereiding/Test Figures/plsmeta2.jpg")
cim(pls.meta)
dev.off()

# The cim still does work, so we save this way.
# lets for now just go on with the multiblock
# first the list and the design

data.meta <- list(aminesdata = amines, lipidsdata = lipids)
design <- matrix(1, ncol = length(data.meta), nrow = length(data.meta),
                 dimnames = list(names(data.meta), names(data.meta)))
diag(design) <-  0
design

# start of block pls

block.pls.meta <- block.pls(X = data.meta, indY = 2, design = design)
plotVar(block.pls.meta, var.names = FALSE, legend = TRUE)
plotIndiv(block.pls.meta, group = comdata2$aggression_status.x, legend = TRUE)

# start of block spls

?block.spls
keepX.list <- list(aminesdata = rep(20,2), lipidsdata = rep(10,2))
block.spls.meta <- block.spls(X = data.meta, indY = 2, keepX = keepX.list, design = design)

# gives errors, and I do not know why
# Lets try pls-da

block.plsda.meta <- block.plsda(X = data.meta, Y = comdata2$aggression_status.x,
                                ncomp = 2, design = design)
block.plsda.meta

# Individual plots

plotIndiv(block.plsda.meta, ind.names = FALSE, legend = TRUE)
plotArrow(block.plsda.meta)
?plotArrow
plotDiablo(block.plsda.meta)

# Variable plots

plotVar(block.plsda.meta, var.names = FALSE, legend = TRUE)
plotLoadings(block.plsda.meta, ncomp = 1, contrib = 'max')

# Lets try the spls-da

block.splsda.meta <- block.splsda(X = data.meta, Y = comdata2$aggression_status.x,
                                  ncomp = 2, keepX = keepX.list, design = design) 
block.splsda.meta

# Individual plots

plotIndiv(block.splsda.meta)
jpeg("C:/Users/Gebruiker/Documents/Master GBH/Stage BehPsy/Voorbereiding/Test Figures/cimdiablo.jpg")
cimDiablo(block.splsda.meta)
dev.off()
plotDiablo(block.splsda.meta)

# Variable plots

plotVar(block.splsda.meta)
jpeg("C:/Users/Gebruiker/Documents/Master GBH/Stage BehPsy/Voorbereiding/Test Figures/circos.jpg")
circosPlot(block.splsda.meta, cutoff = 0.5)
dev.off()

pdf("C:/Users/Gebruiker/Documents/Master GBH/Stage BehPsy/Voorbereiding/Test Figures/network.pdf")
network(block.splsda.meta)
dev.off()

pdf("C:/Users/Gebruiker/Documents/Master GBH/Stage BehPsy/Voorbereiding/Test Figures/loadings.pdf")
plotLoadings(block.splsda.meta, ncomp = 1)
dev.off()

# So all seems to be working well
# For now



