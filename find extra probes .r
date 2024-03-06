#load librarys
library(sesame)
library(minfi)

#load the same sample with minfi and sesame 
sesame <- readIDATpair("/home/minknow/EPIC_V1/B2022.2747_V1_IDAT/205756360010_R01C01")
minfi <- read.metharray.exp("/home/minknow/EPIC_V1/B2022.2747_V1_IDAT/")

#get the Type I and Type II probes 
minfi.typeI <-getProbeInfo(minfi, type = ("I"))
minfi.typeII <-getProbeInfo(minfi, type = ("II"))
#if you compare the size now the sesame DataFrame appears bigger
#but that is because the ctrl and the SNP Probes are also stored in the
#sesame DataFrame. If subtracted the minfi DF is bigger by 379 Probes.


#Get the Probe Info of the xtra Probes that are in minfi but not in sesame
#If you subset the other way ( subset(sesame, !(Probe_ID %in% minfi.typeI$Name))) 
# you get the the SNP and ctrl Probes
sub1 <- subset(minfi.typeI, !(Name %in% sesame$Probe_ID))
sub2 <- subset(minfi.typeII, !(Name %in% sesame$Probe_ID))

#get the minfi Beta Values 
Mset <- preprocessRaw(minfi)
RSet <- ratioConvert(Mset)
minfi.betas <- getBeta(RSet)

#get the sesame Beta Values
sesame.betas <- getBetas(sesame)

#to check you can use the same subset method and you get the same probes
#with beta values 
sesame.names <- names(sesame.betas)
minfi.names<- rownames(minfi.betas)
names.sub <-subset(minfi.names, !(minfi.names %in% sesame.names))
#get the beta values of the extra Probes 
betas.extra <- minfi.betas[c(names.sub),]

##NOTE!!
# I still need to check all the older Manifest files for the Probes


