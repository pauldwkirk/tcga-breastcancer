# tcga-breastcancer
Datasets from TCGA (see: https://tcga-data.nci.nih.gov/docs/publications/brca_2012/ )

Raw datasets:
=============
BRCA.348.precursor.txt               -- miRNA dataset

BRCA.exp.348.med.csv                 -- gene expression dataset

BRCA.Methylation.574probes.802.txt   -- methylation dataset

rppaData-403Samp-171Ab-Trimmed.txt   -- protein dataset 

Processing script:
==================
TCGA_BrCa_dataProcessing.r

Processed datasets:
===================
processedExpression.RDa    -- gene expression  (data matrix dimensions: 645 x 348)

processedMethylation.RDa   -- methylation      (data matrix dimensions: 574 x 348)

processedmiRNA.RDa         -- miRNA            (data matrix dimensions: 423 x 348)  

processedProtein.RDa       -- protein          (data matrix dimensions: 171 x 348)
