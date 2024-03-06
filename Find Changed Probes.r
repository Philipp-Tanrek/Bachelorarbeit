#load librarys
library(sesame)
library(dplyr)

#load the same sample with using sesame as an EPIC v1 and an EPIC v2 array
Ev1 <- readIDATpair("/home/minknow/EPIC_V1/B2022.2747_V1_IDAT/205756360010_R01C01")
Ev2 <- readIDATpair("/home/minknow/EPIC_V2/B2022.2747_V2_IDAT/207097420059_R02C01")

#make the Index into a column so after finding all Probes that changed Type
#its easyer to just select the Probes by Index
Ev1$Index <- rownames(Ev1)
Ev2$Index <- rownames(Ev2)

#remove last 5 characters from Probe_ID because othewise its more difficult
#to compare EPIC v1 and v2 (Ev2 has duplicates which differ in these
#last 5 characters )
Ev2$Probe_ID <-substr(Ev2$Probe_ID, 1, nchar(Ev2$Probe_ID)-5)

#join dfs by Probe_ID and filter out all entries which have the same Probe Type
#the column col holds the Probe Type  (either R , G or 2 ) R and G = Type I
joined_df <- left_join(Ev1, Ev2, by = "Probe_ID")
changed_rows <- filter(joined_df, col.x != col.y)

#Remove duplicates and safe just the Probe_IDs in ChangedProbes 
changed_rows <- changed_rows[!duplicated(changed_rows$Probe_ID),]
ChangedProbes <-changed_rows$Probe_ID

#subset the original Array with the Probe Names 
#Test1 already has the right size but Test 2 has the bevor mentioned duplicates
#which also can have different types.
Test1<- subset(Ev1, Probe_ID %in% ChangedProbes)
Test2<- subset(Ev2, Probe_ID %in% ChangedProbes)

#Initialize DFs
Ev1.changed <- data.frame()
Ev2.changed <- data.frame()

#loop is to filter out duplicates with unchanged Probe Types and only keep 
#the Probes with changed Probe Types 

#Each of the ProbesIn Changed Probes is looped Over
for (i in ChangedProbes) {
  #checks if there areduplicates of this Probe in the EPIC v2
  if (sum(Test2$Probe_ID == i) > 1) {  
    #if there are duplicates loops over these duplicates and checks if they
    #are a different Probe Type to the same Probe in EPIC v1
    # if the Probe is it gets appended to the Dataframe
    for (j in 1:dim(Test2[Test2$Probe_ID == i,])[1]) { 
      #checks if types are different
      if (Test1[Test1$Probe_ID == i,]$col != Test2[Test2$Probe_ID == i,][j,]$col) {
        #adds the row to the df
        Ev1.changed <- rbind(Ev1.changed, Test1[Test1$Probe_ID == i,])
        Ev2.changed <- rbind(Ev2.changed, Test2[Test2$Probe_ID == i,][j,])
      }  
    }
  #for probes with no duplicates the Probe type gets checked
  #and if they are different they get added to the Dataframe
  } else {
    if (Test1[Test1$Probe_ID == i,]$col != Test2[Test2$Probe_ID == i,]$col) {
      Ev1.changed <- rbind(Ev1.changed, Test1[Test1$Probe_ID == i,])
      Ev2.changed <- rbind(Ev2.changed, Test2[Test2$Probe_ID == i,])
    }  
  }
}

# Remove duplicates again to make comparison easier
# probe cg12981137 (Index 450164 & 450167) and
# cg14194875 (Index 495677 & 495678)

# from probe cg12981137 the one with index 450164 is kept 
# from probe cg14194875 the one with index 495677 is kept because the probe 
# with index 495678 has only changed color and not Probe Type
Ev1.changed <- Ev1.changed[!duplicated(Ev1.changed$Probe_ID),]
Ev2.changed <- Ev2.changed[!duplicated(Ev2.changed$Probe_ID),]

# remove last entry which only changed color but not Probe Type 
Ev1.changed <- slice(Ev1.changed, 1:(n()-1))
Ev2.changed <- slice(Ev2.changed, 1:(n()-1))

# I saved the Indexes as a csv (one for EPIC v1 and another for EPIC v2)



