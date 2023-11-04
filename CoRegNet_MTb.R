# Reading and saving .csv files into their respective format for CoRegNet

Mtb_Exp_data <- as.matrix(read.csv(“Exp_matrix_GSE59086.csv”, sep= “,”, header= TRUE, row.names=1))
Mtb_TFs <- read.csv(“Myco_TFs.csv”, header= F)$V1
Mtb_TFs <- as.character(Mtb_TFs)
Mtb_ChIP_OMICS <- as.matrix(read.csv(“Mtb_Network_ChIP.csv”, sep= “,”, header=F))

Save (Mtb_Exp_data, file= “Mtb_Exp_data.RData”)
Save (Mtb_TFs, file= “Mtb_TFs.RData”
Save (Mtb_ChIP_OMICS, file= “Mtb_ChIP_OMICS.RData”)
      
#Executing the hLICORN algorithm in the CoRegNet 
      
GRN_Mtb= hLICORN( Mtb_Exp_data, TFlist= Mtb_TFs)
      
# Refining the inferred regulatory network
      
En_GRN_Mtb = addEvidences(GRN_Mtb,Mtb_ChIP_OMICS)
      
# Supervised refinement of Gene Regulatory Network with ChIP data from MTb Network Portal 
      
refinedGRN = refine(En_GRN_Mtb, integration="supervised", referenceEvidence="Mtb_ChIP_OMICS")
      
#Identification of additional transcriptional programs
      
Mtb_influence= regulatorInfluence(refinedGRN, Mtb_Exp_data)
      
#Display and Visualization of the network
display(refineGRN, expressionData = Mtb_Exp_data, TFA = Mtb_influence
              
#Identifying of master regulator for module detected target class of genes
Target_MEblack = as.character(read.csv("Target_MEblack.csv", header = F)$V1)
save(Target_MEblack, file = "Target_MEblack.RData")
Master_MEblack = masterRegulator(refinedGRN, Target_MEblack)
print(Master_MEblack)
