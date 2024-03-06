#load librarys
library(sesame)
library(tidyverse)

# Initialize two dfs for saving the beta values of the Type I and Type II probes
df <- data.frame()
df2 <- data.frame()

#Load csv of indexes into a vector 
EPICv2.Indexes<- read.csv("Index EPIC v2 changed probes.csv")
EPICv2.Indexes <- EPICv2.Indexes$x

# List all directories of EPIC v2 Samples
List <- list.dirs("/home/minknow/EPIC_V2/", full.names = TRUE)
#remove first entry because ist the Parent folder and it would load it all twice 
List <- List[-1]

#create a range  over all 16 samples
range <- 1:16

# Loop through each directory
for (i in range) {
    # Read data from each directory
    String <- paste(List[i],list.files(List[i])[1], sep = "/")
    #remove the _Grn.idat /_Red.idat to be able to use readIDATpair
    smaller_string <- substr(String, 1, nchar(String)-9)
    #just checkin progress and that its loading the right sample
    print(smaller_string)
    #read sample using sesame into an sdf
    sdf <- readIDATpair(smaller_string)
    #get Probes that changed by Index
    sdf <- sdf[c(EPICv2.Indexes),]
    #subset by type to save them seperatly 
    type1<- subset(sdf, col !=2)
    type2<- subset(sdf, col ==2)
    #orderby Probe_ID so the probes are always in the same order
    #makes it easier to make scatteplots to see correlations
    type1 <- type1[order(type1$Probe_ID),]
    type2 <- type2[order(type2$Probe_ID),]
    #get the Beta Values
    betas1 <- getBetas(type1)
    betas2 <- getBetas(type2)
    # Create a data frame from the read data
    frame1 <- data.frame(betas1)
    frame2 <- data.frame(betas2)
    #if its the first itteration of the loop initialize the df
    #else just add a coulumn with the beta values 
    if (i == 1) {
        df <- frame1
        df2 <- frame2
    } else {
        # add beta values as a column
        df <- cbind(df, frame1)
        df2 <- cbind(df2, frame2)
    }
}

#create empty list 
result_list <- vector("list", length = length(List))

# Loop over each element in the List
for (i in seq_along(List)) {
  # Get the names of the Samples form the Folders and save them in a list
  result_list[[i]] <- substr(List[i], 30, nchar(List[i])-8)
}

#Names the columns with the according Sample Name 
names(df) <- result_list
names(df2) <- result_list

# quick check of the dfs 
head(df)
head(df2)

#save the dfs as a csv to later plot 
write.csv(df, "bbbbbbbbbbbbbb.csv")#Changed Type I Probes EPIC v2
write.csv(df2, "bbbbbbbbbbbbb.csv")#Changed Type II Probes EPIC v2
