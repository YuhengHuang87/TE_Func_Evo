## scripts for the functional and evolutionary consequence of TEs' 3D interaction with PCH
#### TE-medicated H3K9me3 enrichment analysis
1. cut_tag_mapping_bowtie2_bed_hifi_genome.sub
2. fragment_CUT_Tag_bedgraph_Hifi.pl
3. adjacent_H3K9me3_enrichment_estimate_Hifi.pl
4. avg_K9_per_TE_hifi_genome.pl


#### expression rank for TE adjacent genes
1. RNA_seq_TE_z_score.R
2. assign_TE_exp_z_Hifi_genome.pl


#### pannagram to estimate TE ages
1. pannagram_shared.sub
2. pangenome_hubplot.R #generates the Connection graph. Code was written with assistence from Gemini 2.5 Pro


#### Figure generation and statistical test
3D_TE_functional_similarity.qmd #generate Fig. 3C & E, Fig. 4A, C & D, Supplementary Fig 11-15
