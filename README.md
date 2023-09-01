# Consistency-Index
This repository contains the required script and associated metadata to run a consistency index analysis between different sets of genes using the same phylogeny.
The script has been provided as an rscript and can be run locally in Rstudio or Visual Studio Core.

## Publication
The methodology underlying the use of the script has been detailed in the article "Gene pool transmission of multidrug resistance among _Campylobacter_ from livestock, sewage and human disease" published in Environmental Microbiology,
and in the article "Genome evolution and the emergence of pathogenicity in avian _Escherichia coli_". Please refer to citation information at the bottom of this document.

## Files summary
Two different datasets have been used to quantify the consistency index of two different sets of genes on the same phylogeny. The provided datasets used in the manuscript have been provided as input multi-fasta alignments for each gene of interest.
Each directory has the following:
* **_Core_genes_**: directory with core gene fasta files
* **_AMR_genes_**: directory with AMR gene multi-fasta files
* **_Core_genes.txt_**: path to Core_genes directory
* **_AMR_genes.txt_**: path to AMR_genes directory
* **_ML_tree.nwk_**: ML tree rooted at midpoint root

# Consistency Index script
## Libraries
* ape
* phangorn
* binr

## Set working directory
```setwd("/Path/to/working/directory")```

## Import and read files
```
core_genes_file <- "Core_genes.txt"
f1 <- read.table(core_genes_file, stringsAsFactors = F)

AMR_fasta_file <- "AMR_genes.txt"
f2 <- read.table(AMR_fasta_file, stringsAsFactors = F)

ML_tree <- "ML_tree.nwk"
phyloTree = midpoint(read.tree(ML_tree))
```
### Core genes section
```
outputs_core <- NULL

for(i in 1:length(f1$V1)){
  
  core_gene_file <- f1$V1[i]
  
  seqdata = read.FASTA(core_gene_file)
  seqdataPhy = phyDat(seqdata, type='DNA')
  
  phyloTree = midpoint(read.tree(ML_tree))
  
  p1 <- phyloTree$tip.label
  p2 <- names(seqdata)
  p3 <- setdiff(p1,p2)
  
  t2 <- drop.tip(phyloTree,as.character(p3))
  plot(phyloTree)
  plot(t2)
  
  ci <- CI(t2, seqdataPhy, cost = NULL, sitewise = FALSE)
  outputs_core[i] <- ci
  
}

write.table(outputs_core, file = "Core_CI_Outputs.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
```
### AMR genes section
```
outputs <- NULL

for(i in 1:length(f2$V1)){
  
  AMR_fasta_file <- sprintf("%s.fas", as.character(f2$V1[i]))
  
  seqdata = read.FASTA(AMR_fasta_file)
  seqdataPhy = phyDat(seqdata, type='DNA')
  
  phyloTree = midpoint(read.tree(ML_tree))
  
  p1 <- phyloTree$tip.label
  p2 <- names(seqdata)
  p3 <- setdiff(p1,p2)
  
  t2 <- drop.tip(phyloTree,as.character(p3))
  plot(phyloTree)
  plot(t2)
  
  ci <- CI(t2, seqdataPhy, cost = NULL, sitewise = FALSE)
  outputs[i] <- ci
  
}

write.table(outputs, file = "AMR_CI_Outputs.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
```
## Example 1 - Consistency index AMR vs core genes in _Campylobacter_ 
![CI_Mourkas_et_al_2019](https://github.com/evangelosmourkas/Consistency-Index/assets/73548463/a0a9e030-bb01-404b-b4de-e14f4262d446)

## Example 2 - Consistency index virulence vs core genes in _Escherichia coli_
![Mageiros_et_al](https://github.com/evangelosmourkas/Consistency-Index/assets/73548463/6b3351ed-b9ee-4b5d-89f9-c74b9040bf41)

# How to cite
Mourkas, E., Florez-Cuadrado, D., Pascoe, B., Calland, J.K., Bayliss, S.C., Mageiros, L., Méric, G., Hitchings, M.D., Quesada, A., Porrero, C., Ugarte-Ruiz, M., Gutiérrez-Fernández, J., Domínguez, L. and Sheppard, S.K. (2019), **Gene pool transmission of multidrug resistance among _Campylobacter_ from livestock, sewage and human disease**. Environ Microbiol, 21: 4597-4613. https://doi.org/10.1111/1462-2920.14760

Mageiros, L., Méric, G., Bayliss, S.C., Pensar, J., Pascoe, B., Mourkas, E., Calland, J.K., Yahara, K., Murray, S., Wilkinson, T.S., Williams, L.K., Hitchings, M.D., Porter, J., Kemmett, K., Feil, E.J., Jolley, K.A., Williams, N.J., Corander, J. and Sheppard, S.K. (2021) **Genome evolution and the emergence of pathogenicity in avian _Escherichia coli_**. Nat Commun 12, 765. https://doi.org/10.1038/s41467-021-20988-w
