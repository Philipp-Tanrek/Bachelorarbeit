Code for my Bachelor thesis

find extra probes.r :
Finds Probes that are in the minfi Manifest file but not the sesame manifest File 

R: 4.3.2
Bioconductor: 3.18
sesame: 1.20.0
sesameData: 1.20.0
ExperimentHub: 2.10.0
minfi_1.48.0
IlluminaHumanMethylationEPICmanifest_0.3.0

Find Changed Probes.r :
Finds and saves the Indexes of Probes that changed type in the newer EPICv2
R: 4.3.2
Bioconductor: 3.18
sesame: 1.20.0
sesameData: 1.20.0
ExperimentHub: 2.10.0
tidyverse 2.0.0

Create CSV of Beta Values from changed Probes (EPICv1).r and 
Create CSV of Beta Values from changed Probes (EPICv2).r :
Read all Idat Files in a directory with either EPIC v1 or EPIC v2 IDAT files.
Extract the changed Probes by Index and save in a Df.
Save the DF as a CSV for plotting and analysis.
R: 4.3.2
Bioconductor: 3.18
sesame: 1.20.0
sesameData: 1.20.0
ExperimentHub: 2.10.0
tidyverse 2.0.0


Delta-Beta_original Sample.r:
Compares EPIC v1 to EPIC v2 files, calculates the delta beta value, and stores them in a DF where each column is a sample. Lastly, calculates the Mean of each row/ CpG/ probe.


Binarize_original Sample.r:
Binarizes EPIC v1 and EPIC v2 Beta-Values compares them and stores the result in a DF, where each column is a sample. Lastly, sums each row/ CpG/ probe.

Trimerize_original Sample.r:
Trimerizes EPIC v1 and EPIC v2 Beta-Values, compares them, and stores the result in a DF, where each column is a sample. Lastly, sums each row/ CpG/ probe.

Identify differing Probes.R: 
Loads the DFs produced in   Delta-Beta_original Sample.r, Binarize_original Sample.r and Trimerize_original Sample.r (for both the Kaur and USB data set) and Identifies differing probes and saves them in a CSV

UMAP Analysis Code.ipynb:
Code used to load and analyze Excel files generated by the EpiDip UMAP algorithm.

