

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

#save a csv
save.csv <- function(csv, file.Name){
  write.csv(csv, file.Name)
}



#create trimerization function to always get either methylated,mixed or unmethylated 
trimerize <- function(df) {
  df[df >= 0.666, 1] <-1
  df[df < 0.666 & df >= 0.333, 1] <- 0.5
  df[df < 0.333, 1] <-0
  return(df)
}



#read csv of overlapping Names to keep for EPIC v1 and EPIC v2
probe.names.epic.v1 <- read.csv("Names of Overlapping Probes EPIC v1.csv")
probe.names.epic.v1 <- probe.names.epic.v1$x

probe.names.epic.v2 <- read.csv("Names of Overlapping Probes EPIC v2.csv")
probe.names.epic.v2 <- probe.names.epic.v2$x

#initialize Dataframe for saving the data 
Ev1.trimerized.raw <- data.frame()
Ev2.trimerized.raw <- data.frame()
Ev1.trimerized.swan <- data.frame()
Ev2.trimerized.swan <- data.frame()

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
  
  #binarize the beta values 
  betas.ev1.raw <- trimerize(betas.ev1.raw)
  betas.ev2.raw <- trimerize(betas.ev2.raw)

  betas.ev1.swan <- trimerize(betas.ev1.swan)
  betas.ev2.swan <- trimerize(betas.ev2.swan)
  
  #order the data frame by probe name so that when saving one row corresponds to one column 
  betas.ev1.raw <- betas.ev1.raw[order(rownames(betas.ev1.raw), decreasing = TRUE),]
  betas.ev2.raw <- betas.ev2.raw[order(substr(rownames(betas.ev2.raw),1,nchar(rownames(betas.ev2.raw))-5), decreasing = TRUE),]
  betas.ev1.swan <- betas.ev1.swan[order(rownames(betas.ev1.swan), decreasing = TRUE),]
  betas.ev2.swan <- betas.ev2.swan[order(substr(rownames(betas.ev2.swan),1,nchar(rownames(betas.ev2.swan))-5), decreasing = TRUE),]
  
  
  if (i == 1) {
    #for the first iteration initialize the data frame with abs difference 
    Ev1.trimerized.raw <- as.data.frame(betas.ev1.raw)
    Ev2.trimerized.raw <- as.data.frame(betas.ev2.raw) 
    
    Ev1.trimerized.swan <- as.data.frame(betas.ev1.swan) 
    Ev2.trimerized.swan <- as.data.frame(betas.ev2.swan) 
    
    rownames(Ev1.trimerized.raw) <- rownames(betas.ev2.raw)
    rownames(Ev2.trimerized.raw) <- rownames(betas.ev2.swan)
    
    rownames(Ev1.trimerized.swan) <- rownames(betas.ev2.raw)
    rownames(Ev2.trimerized.swan) <- rownames(betas.ev2.swan)
  } else {
    # add abs difference as a column 
    Ev1.trimerized.raw <- cbind(Ev1.trimerized.raw, betas.ev1.raw)
    Ev2.trimerized.raw <- cbind(Ev2.trimerized.raw, betas.ev2.raw) 

    Ev1.trimerized.swan <- cbind(Ev1.trimerized.swan, betas.ev1.raw) 
    Ev2.trimerized.swan <- cbind(Ev2.trimerized.swan, betas.ev2.raw) 
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
names(Ev1.trimerized.raw) <- result_list
names(Ev2.trimerized.raw) <- result_list
names(Ev1.trimerized.swan) <- result_list
names(Ev2.trimerized.swan) <- result_list



#save data as a csv
save.csv(Ev1.trimerized.raw, "Ev1_Trimerized_raw.csv")
save.csv(Ev2.trimerized.raw, "Ev2_Trimerized_raw.csv")
save.csv(Ev1.trimerized.swan, "Ev1_Trimerized_swan.csv")
save.csv(Ev2.trimerized.swan, "Ev2_Trimerized_swan.csv")


#find differing probes 
#create a copy (only for the structure and order of the df to make it easier)
Copy <- Ev2.trimerized.swan
#loop overeach column in both dfs and if the are equal it is set to 0
#if they differ it is set to 1
for (column in 1:16) {
  temp <- (Ev1.trimerized.swan[,column] == Ev2.trimerized.swan[,column])
  print(head(temp))
  Copy[,column] <- ifelse(temp, 0, 1)
}
head(Copy)
#get the sum of each row and add it as a column 
temp <- rowSums(Copy,na.rm = TRUE)
Copy <- cbind(Copy, Sum = temp)
#order the df by how many times the trimerized beta value differs between EPIC v1 and EPIC v2
Copy <- Copy[order(Copy$Sum, decreasing = TRUE),]
#save as csv
save.csv(Copy, "Trimerized_Original_Probes_Ranked.csv")


