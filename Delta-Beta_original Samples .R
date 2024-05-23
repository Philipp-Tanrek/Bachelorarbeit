#load library s 
library(minfi)
library(tidyverse)

# convert Red-Green Channel Set to Genomic Ratio Set and extract beta values without any preprocessing 
combinedpreprocessing.Raw <- function(RG_Channel_Set){
  Mset <-preprocessRaw(RG_Channel_Set)
  Rset <-ratioConvert(Mset)
  Gset <-mapToGenome(Rset)
  betas <- getBeta(Gset)
  return(betas)
}

# convert Red-Green Channel Set to Genomic Ratio Set and extract beta values with SWAN preprocessing 
combinedpreprocessing.SWAN <- function(RG_Channel_Set){
  Mset <-preprocessSWAN(RG_Channel_Set)
  Rset <-ratioConvert(Mset)
  Gset <-mapToGenome(Rset)
  betas <- getBeta(Gset)
  return(betas)
}

#remove all beta values in a vector to only get the overlap (overlap already determined and saved in a CSV)
filter.betas <- function(BetaValues, Names){
  filtered.data <- subset(BetaValues, rownames(BetaValues) %in% Names)
  return(filtered.data)
}

#calculate the Absolute difference of beta values between EPIC v1 and EPIC v2 for every  probe 
get.abs.diff.as.df <- function(Betas.Ev1, Betas.Ev2){
  ordered.ev1 <- Betas.Ev1[order(rownames(Betas.Ev1), decreasing = TRUE),]
  ordered.ev2 <- Betas.Ev2[order(rownames(Betas.Ev2), decreasing = TRUE),]
  difference <- abs(ordered.ev1-ordered.ev2)
  difference <- data.frame(difference)
  return(difference)
}

#get the mean of each row and orer the data frame accordig to the mean row from highest to lowest 
get.row.mean <- function(dataframe){
  Means = apply(dataframe, 1,mean)
  dataframe <- cbind(dataframe, Means)
  dataframe <- dataframe[order(dataframe$Means, decreasing = TRUE),]
  return(dataframe)
}

#save a csv
save.csv <- function(csv, file.Name){
  write.csv(csv, file.Name)
}




#read csv of overlapping Names to keep for EPIC v1 and EPIC v2
probe.names.epic.v1 <- read.csv("Names of Overlapping Probes EPIC v1.csv")
probe.names.epic.v1 <- probe.names.epic.v1$x

probe.names.epic.v2 <- read.csv("Names of Overlapping Probes EPIC v2.csv")
probe.names.epic.v2 <- probe.names.epic.v2$x

#initialize Dataframe for saving the data 
abs.diff.raw <- data.frame()
abs.diff.swan <- data.frame()

# List all directories of EPIC v1 Samples
ListEv1 <- list.dirs("/home/minknow/EPIC_V1/", full.names = TRUE)
ListEv2 <- list.dirs("/home/minknow/EPIC_V2/", full.names = TRUE)
#remove first entry because its the Parent folder and it would load it all twice 
ListEv1 <- ListEv1[-1]
ListEv2 <- ListEv2[-1]

#get the range for the loop 
range <- 1:length(ListEv1)



for (i in range) {
  # Get the name of the folder + names of the file 
  StringEv1 <- paste(ListEv1[i],"/" , sep = "")
  StringEv2 <- paste(ListEv2[i],"/" , sep = "")
  #just checking progress and that its loading the right sample
  print(StringEv1)
  print(StringEv2)
  #read sample using minfi into an RGset
  Ev1 <- read.metharray.exp(StringEv1)
  Ev2 <- read.metharray.exp(StringEv2)
  #manualy set the annotation package for EPIC v2
  Ev2@annotation <-c(array = "IlluminaHumanMethylationEPICv2", annotation = "20a1.hg38")
  
  #get the minfi Beta Values once raw and once with SWAN preprocessing 
  betas.ev1.raw <- combinedpreprocessing.Raw(Ev1)
  betas.ev2.raw <- combinedpreprocessing.Raw(Ev2)

  betas.ev1.swan <- combinedpreprocessing.SWAN(Ev1)
  betas.ev2.swan <- combinedpreprocessing.SWAN(Ev2)
  
  #filter out all probes that aren't in the overlap of HM450K EPIC v1 and EPIC v2 plus remove all cross reactive probes
  betas.ev1.raw <- filter.betas(betas.ev1.raw, probe.names.epic.v1)
  betas.ev2.raw <- filter.betas(betas.ev2.raw, probe.names.epic.v2)
  
  betas.ev1.swan <- filter.betas(betas.ev1.swan, probe.names.epic.v1)
  betas.ev2.swan <- filter.betas(betas.ev2.swan, probe.names.epic.v2)
  
  #get the Absolute difference for every probe
  diff.raw <-get.abs.diff.as.df(betas.ev1.raw, betas.ev2.raw)
  diff.swan <- get.abs.diff.as.df(betas.ev1.swan, betas.ev2.swan)

  if (i == 1) {
    #for the first iteration initialize the data frame with abs difference 
    abs.diff.raw <- diff.raw
    abs.diff.swan <- diff.swan
    rownames(abs.diff.raw) <- rownames(betas.ev2.raw)
    rownames(abs.diff.swan) <- rownames(betas.ev2.swan)
  } else {
    # add abs difference as a column 
    abs.diff.raw <- cbind(abs.diff.raw, diff.raw)
    abs.diff.swan <- cbind(abs.diff.swan, diff.swan)
  }
}

#create empty list for names 
result_list <- vector("list", length = length(ListEv1))

# Loop over each element in the List
for (i in seq_along(ListEv1)) {
  # Get the names of the Samples form the Folders and save them in a list
  result_list[[i]] <- substr(ListEv1[i], 30, nchar(ListEv1[i])-8)
}

#Names the columns with the according Sample Name 
names(abs.diff.raw) <- result_list
names(abs.diff.swan) <- result_list

#get rowmean plus sort according to row mean to get most differing probes 
abs.diff.raw <- get.row.mean(abs.diff.raw)
abs.diff.swan <- get.row.mean(abs.diff.swan)


#order according to means again just to be safe 
abs.diff.raw <- abs.diff.raw[order(abs.diff.raw$Means, decreasing = TRUE),]
abs.diff.swan <- abs.diff.swan[order(abs.diff.swan$Means, decreasing = TRUE),]

#save data as a csv
save.csv(abs.diff.raw, "abs diff raw.csv")
save.csv(abs.diff.swan, "abs diff swan.csv")


