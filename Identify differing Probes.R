library(VennDiagram)
library(ggplot2)

##############################
#############################
############################
#Find Differing Probes using the Delta Beta Data
############################
#############################
##############################


#load SWAN preprocessed Data (delta Betavalues)
data.swan <- read.csv("abs diff swan.csv")
extra.data.swan <- read.csv("Extra Data abs diff swan.csv")

#load probe names from UMAP Index 
index <- read.csv("index.csv", col.names = FALSE)

#calculate Median for every Row(so for every probe)
Median <- apply(data.swan[, 2:17], 1, median)
data.swan$Median <- Median

extra.Median <- apply(extra.data.swan[, 2:12], 1, median)
extra.data.swan$Median <- extra.Median

########calculate cutoff###################
#cutoff is defined as Median plus 2 standard deviations 

#cutoff for Original Data
median.data <- median(data.swan$Median)
sd.data <- sd(data.swan$Median)
cutoff.data <- median.data + (2*sd.data)


#cutoff for extra Data
median.extra <- median(extra.data.swan$Median)
sd.extra <- sd(extra.data.swan$Median)
cutoff.extra <- median.extra + (2*sd.extra)


#numbers of sample above cutoff
num.data <- sum(data.swan$Median >= cutoff.data)
num.extra <- sum(extra.data.swan$Median >= cutoff.extra)

#order data by median (highest median at the top lowestat the bottom of the DF)
data.swan <- data.swan[order(data.swan$Median,decreasing = TRUE),]
extra.data.swan <- extra.data.swan[order(extra.data.swan$Median,decreasing = TRUE),]

#get Top Data
Top.data <- data.swan[1:num.data,]
Top.extra <- extra.data.swan[1:num.extra,]

#get Overlap
Overlap <- subset(Top.data, Top.data$X %in% Top.extra$X)

#remove data not in index 
Overlap <- subset(Overlap, substr(Overlap$X,1,nchar(Overlap$X)-5) %in% index$FALSE.)

#read original UMAP Blacklist and compare to Overlap
Blacklist <- read.csv("infiniumProbeBlacklist.csv", col.names = FALSE)
sum(!(substr(Overlap$X,1,nchar(Overlap$X)-5) %in% Blacklist))
#change to EPIC v1 annotation  and save as csv
delta.overlap <- substr(Overlap$X,1,nchar(Overlap$X)-5)
delta.overlap <- as.data.frame(delta.overlap)
names(delta.overlap) <- names(Blacklist)
delta.blacklist <- rbind(Blacklist, delta.overlap)
write.csv(delta.blacklist, "Delta Beta Blacklist.csv", row.names = FALSE, col.names = FALSE)



##############################
#############################
############################
#Find Differing Probes using the Binarized Data
############################
#############################
##############################


#load Binarized Data with the number of times a Probe differs between EPIC v1 and EPIC v2 in the sum column
data <- read.csv("Binarized_Original_Probes_Ranked.csv")
extra <- read.csv("Binarized_EXTRA_Probes_Ranked.csv")

data <- subset(data, substr(data$X, 1, nchar(data$X)-5) %in% index$FALSE.)
extra <- subset(extra, substr(extra$X, 1, nchar(extra$X)-5) %in% index$FALSE.)



num = 2
data.sub <- subset(data, data$Sum >= num)
extra.sub <- subset(extra, extra$Sum >= num)
head(data.sub)
dim(data.sub)
dim(extra.sub)

#get Overlap and write out as CSV file 
binarization.overlap <- subset(data.sub, data.sub$X %in% extra.sub$X)
binarization.overlap <- substr(binarization.overlap$X,1,nchar(binarization.overlap$X)-5)
binarization.overlap <- as.data.frame(binarization.overlap)
names(binarization.overlap) <- names(Blacklist)
Binarization.Blacklist <- rbind(Blacklist, binarization.overlap)

write.csv(Binarization.Blacklist, "Binarization Blacklist.csv", row.names = FALSE)


##############################
#############################
############################
#Find Differing Probes using the Trimerized Data
############################
#############################
##############################


#load Trimerized Data with the number of times a Probe differs between EPIC v1 and EPIC v2 in the sum column
data.tri <- read.csv("Trimerized_Original_Probes_Ranked.csv")
extra.tri <- read.csv("Trimerized_EXTRA_Probes_Ranked.csv")

data.tri <- subset(data.tri, substr(data.tri$X, 1, nchar(data.tri$X)-5) %in% index$FALSE.)
extra.tri <- subset(extra.tri, substr(extra.tri$X, 1, nchar(extra.tri$X)-5) %in% index$FALSE.)



num = 2
data.tri.sub <- subset(data.tri, data.tri$Sum >= num)
extra.tri.sub <- subset(extra.tri, extra.tri$Sum >= num)
head(data.tri.sub)
dim(data.tri.sub)
dim(extra.tri.sub)

trimerization.overlap <- subset(data.tri.sub, data.tri.sub$X %in% extra.tri.sub$X)
trimerization.overlap <- substr(trimerization.overlap$X,1,nchar(trimerization.overlap$X)-5)
trimerization.overlap <- as.data.frame(trimerization.overlap)
names(trimerization.overlap) <- names(Blacklist)
Trimerization.Blacklist <- rbind(Blacklist, trimerization.overlap)

write.csv(Trimerization.Blacklist, "Trimerization Blacklist.csv", row.names = FALSE)


