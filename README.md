## scripts for the functional and evolutionary consequence of TEs' 3D interaction with PCH
#### TE-medicated H3K9me3 enrichment analysis
1. cut_tag_mapping_bowtie2_bed_hifi_genome.sub #map H3K9me3 CUT&Tag data to their respective Hifi genomes.
2. fragment_CUT_Tag_bedgraph_Hifi.pl #calculate K9 fragments per 25 bps region
3. adjacent_H3K9me3_enrichment_estimate_Hifi.pl #estimate H3K9me3 enrichment for TE flanking regions.
4. avg_K9_per_TE_hifi_genome.pl #average the enrichment across replicates


#### expression rank for TE adjacent genes
1. RNA_seq_TE_z_score.R # calculate the expression rank based on RPKM
2. assign_TE_exp_z_Hifi_genome.pl # assign the expression measures to the nearby TEs


#### pannagram to estimate TE ages
1. pannagram_shared.sub
2. pangenome_hubplot.R #generates the Connection graph. Code was written with assistence from Gemini 2.5 Pro


#### Figure generation and statistical test
3D_TE_functional_similarity.qmd #generate Fig. 3C & E, Fig. 4A, C & D, Supplementary Fig 11-15
