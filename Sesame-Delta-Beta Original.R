#load library s 
library(sesame)
library(tidyverse)

get_betas <- function(path_to_IDAT_pair){
  betas <- openSesame(path_to_IDAT_pair, prep = "CDB", func = getBetas, collapseToPfx = TRUE)
  return(betas)
}

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
  #BetaValues <- as.data.frame(BetaValues)
  #print(head(BetaValues))
  filtered.data <- subset(BetaValues, names(BetaValues) %in% Names)
  return(filtered.data)
}

#calculate the Absolute difference of beta values between EPIC v1 and EPIC v2 for every  probe 
get.abs.diff.as.df <- function(Betas.Ev1, Betas.Ev2){
  #ordered.ev1 <- Betas.Ev1[order(rownames(Betas.Ev1), decreasing = TRUE),]
  #ordered.ev2 <- Betas.Ev2[order(rownames(Betas.Ev2), decreasing = TRUE),]
  difference <- abs(Betas.Ev1-Betas.Ev2)
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

#initialize Dataframe for saving the data 
abs.diff.sesame <- data.frame()

# List all directories of EPIC v1 Samples
ListEv1 <- list.dirs("/home/minknow/EPIC_V1/", full.names = TRUE)
ListEv2 <- list.dirs("/home/minknow/EPIC_V2/", full.names = TRUE)
#remove first entry because its the Parent folder and it would load it all twice 
ListEv1 <- ListEv1[-1]
ListEv2 <- ListEv2[-1]

#get the range for the loop 
range <- 1:length(ListEv1)



for (i in range) {
  print(ListEv1[i])
  print(ListEv2[i])
  #get the minfi Beta Values once raw and once with SWAN preprocessing 
  
  betas.ev1 <- get_betas(ListEv1[i])
  betas.ev2 <- get_betas(ListEv2[i])
  print("Raw Beta values")
  print(head(betas.ev1, n= 10))
  print(head(betas.ev2, n= 10))
  
  #filter out all probes that aren't in the overlap of HM450K EPIC v1 and EPIC v2 plus remove all cross reactive probes
  
  betas.ev1 <- filter.betas(betas.ev1, probe.names.epic.v1)
  betas.ev2 <- filter.betas(betas.ev2, probe.names.epic.v1)
  print("filtered Beta values")
  print(head(betas.ev1, n= 10))
  print(head(betas.ev2, n= 10))
  
  #get the Absolute difference for every probe
  diff.sesame <- get.abs.diff.as.df(betas.ev1, betas.ev2)
  print("Abs diff values")
  print(head(diff.sesame))
  
  if (i == 1) {
    #for the first iteration initialize the data frame with abs difference 
    abs.diff.sesame <- diff.sesame
    rownames(abs.diff.sesame) <- rownames(betas.ev2)
    print("DF")
    print(head(abs.diff.sesame))
  } else {
    # add abs difference as a column 
    print("DF")
    abs.diff.sesame <- cbind(abs.diff.sesame, diff.sesame)
    print(head(abs.diff.sesame))
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
names(abs.diff.sesame) <- result_list

#get rowmean plus sort according to row mean to get most differing probes 
abs.diff.sesame <- get.row.mean(abs.diff.sesame)


#order according to means again just to be safe 
abs.diff.sesame <- abs.diff.sesame[order(abs.diff.sesame$Means, decreasing = TRUE),]

#save data as a csv
save.csv(abs.diff.sesame, "abs diff sesame.csv")


