#Consistency index

library(ape)
library(phangorn)
library(binr)

setwd("C:/Users/biol0023/OneDrive - Nexus365/Documents/Oxford_Files/1.Projects/1.PostDoc/8.Peru/1.cmeB/1.CI_analysis")

AMR_fasta_file <- "AMR_fasta.txt"

f1 <- read.table(AMR_fasta_file, stringsAsFactors = F)

core_genes_file <- "core_genes.txt"

f2 <- read.table(core_genes_file, stringsAsFactors = F)

Campy_tree <- "Renamed_Rooted_RAxML_bestTree.cmeB_new_core_T14.nwk"
phyloTree = midpoint(read.tree(Campy_tree))


#AMR gene section
outputs <- NULL

for(i in 1:length(f1$V1)){
  
  AMR_fasta_file <- sprintf("%s.fas", as.character(f1$V1[i]))
  
  seqdata = read.FASTA(AMR_fasta_file)
  seqdataPhy = phyDat(seqdata, type='DNA')
  
  phyloTree = midpoint(read.tree(Campy_tree))
  
  p1 <- phyloTree$tip.label
  p2 <- names(seqdata)
  p3 <- setdiff(p1,p2)
  
  t2 <- drop.tip(phyloTree,as.character(p3))
  plot(phyloTree)
  plot(t2)
  
  ci <- CI(t2, seqdataPhy, cost = NULL, sitewise = FALSE)
  outputs[i] <- ci
  
}

write.table(outputs, file = "Peru_AMR_CI_Outputs", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

#core gene section
outputs_core <- NULL

for(i in 1:length(f2$V1)){
  
  core_gene_file <- f2$V1[i]
  
  seqdata = read.FASTA(core_gene_file)
  seqdataPhy = phyDat(seqdata, type='DNA')
  
  phyloTree = midpoint(read.tree(Campy_tree))
  
  p1 <- phyloTree$tip.label
  p2 <- names(seqdata)
  p3 <- setdiff(p1,p2)
  
  t2 <- drop.tip(phyloTree,as.character(p3))
  plot(phyloTree)
  plot(t2)
  
  ci <- CI(t2, seqdataPhy, cost = NULL, sitewise = FALSE)
  outputs_core[i] <- ci
  
}

write.table(outputs_core, file = "Peru_core_genes_CI_Outputs", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
