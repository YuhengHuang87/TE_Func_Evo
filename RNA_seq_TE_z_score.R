library(tidyverse)
rna <- read_tsv("/Users/yuhenghuang/Documents/Postdoc_UCI/RNA-seq/A4_ref_RSEM_rpkm.combined.minRPKM2.txt")
as_tibble(rna)

############ PCA analysis
pca_matrix <- rna %>% 
  column_to_rownames("gene_name") %>% 
  as.matrix() %>% 
  t()
sample_pca <- prcomp(pca_matrix)
pc_scores <- sample_pca$x
pc_scores <- pc_scores %>% 
  # convert to a tibble retaining the sample names as a new column
  as_tibble(rownames = "sample")

# print the result
pc_scores
pc_scores %>% 
  # create the plot
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point()
###########

rna[paste0("A4_", 1:4, "rank")] <- apply(-rna[,2:5], 2, rank);rna[paste0("A7_", 1:4, "rank")] <- apply(-rna[,6:9], 2, rank)
for (i in 1:nrow(rna)){
v_A4=c(rna$A4_1[i],rna$A4_2[i],rna$A4_3[i],rna$A4_4[i]); v_A7=c(rna$A7_1[i],rna$A7_2[i],rna$A7_3[i],rna$A7_4[i]);  r_A4=c(rna$A4_1rank[i],rna$A4_2rank[i],rna$A4_3rank[i],rna$A4_4rank[i]); r_A7=c(rna$A7_1rank[i],rna$A7_2rank[i],rna$A7_3rank[i],rna$A7_4rank[i])
rna$z_rpkm[i] = (mean(v_A4)-mean(v_A7))/((var(v_A4)/2)+(var(v_A7)/2))^0.5;rna$z_rank[i] = (mean(r_A4)-mean(r_A7))/((var(r_A4)/2)+(var(r_A7)/2))^0.5
}
write.table(rna, file="/Users/yuhenghuang/Documents/Postdoc_UCI/RNA-seq/A4_ref_RSEM_rpkm_rank_z.txt",sep="\t",eol="\n",row.names = F,col.names = F)

rna <- read_tsv("/Users/yuhenghuang/Documents/Postdoc_UCI/RNA-seq/A7_ref_RSEM_rpkm.combined.minRPKM2.txt")
as_tibble(rna)
rna[paste0("A7_", 1:4, "rank")] <- apply(-rna[,2:5], 2, rank);rna[paste0("A4_", 1:4, "rank")] <- apply(-rna[,6:9], 2, rank)
for (i in 1:nrow(rna)){
  v_A4=c(rna$A4_1[i],rna$A4_2[i],rna$A4_3[i],rna$A4_4[i]); v_A7=c(rna$A7_1[i],rna$A7_2[i],rna$A7_3[i],rna$A7_4[i]);  r_A4=c(rna$A4_1rank[i],rna$A4_2rank[i],rna$A4_3rank[i],rna$A4_4rank[i]); r_A7=c(rna$A7_1rank[i],rna$A7_2rank[i],rna$A7_3rank[i],rna$A7_4rank[i])
  rna$z_rpkm[i] = (mean(v_A7)-mean(v_A4))/((var(v_A4)/2)+(var(v_A7)/2))^0.5;rna$z_rank[i] = (mean(r_A7)-mean(r_A4))/((var(r_A4)/2)+(var(r_A7)/2))^0.5
}
write.table(rna, file="/Users/yuhenghuang/Documents/Postdoc_UCI/RNA-seq/A7_ref_RSEM_rpkm_rank_z.txt",sep="\t",eol="\n",row.names = F,col.names = F)


