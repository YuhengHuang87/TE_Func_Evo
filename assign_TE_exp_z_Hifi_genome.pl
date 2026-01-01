#!/usr/bin/perl
use strict;
use warnings;

my $k=0;
my @name=("A4","A7");
my $infile1="/dfs7/grylee/hgshukla/HiC_Project/A4/sample1/UU_minmapq_0/code_process_UU_MU_parallel/LNorm/UU_MU/A4_A7.UU_MU.Zscore.SigFlag.txt";
#my $infile1="/dfs7/grylee/hgshukla/HiC_Project/A7/sample1/code_process_MU_UU_Parallel/LNorm/UU_MU/A7_A4.UU_MU.Zscore.SigFlag.txt";

my $path = "/dfs7/grylee/yuhenh3/Six_species_2019";
my %TE_enrich;
my $win_leng=999999999;
my $TE_num=0;

my $outfile=$name[$k]."_interacting_TE_nearby_gene_expression.txt";
#my $outfile=$name[$k]."_Zita_K9_interacting_TE_nearby_gene_expression.txt";
open(OUT, ">$outfile");

open(FILE1,"<", "$infile1")||die"$!";
	while(my $count = <FILE1>){
		chomp($count);
		my @b = split("\t", $count);
		if ($b[0] eq "CHR"){}else{
			#my $left=$b[1]; my $right=$b[2];
my $locat=$b[0]."\t".$b[1]."\t".$b[2];#for raw TE position
$TE_enrich{$locat}=$count;
$TE_num++;
}}
#print "$name[$k]\t","$TE_num\n";
my %TE_K9_height;
my $infile3="/dfs7/grylee/yuhenh3/HC_HiFi/".$name[$k]."_mean_K9_5kb_TE_flanking_left_right.txt";
#my $infile3="/dfs7/grylee/yuhenh3/HC_HiFi/".$name[$k]."_Zita_mean_K9_5kb_TE_flanking_left_right.txt";

open(FILE1,"<", "$infile3")||die"$!";
	while(my $count = <FILE1>){
		chomp($count);
		my @a = split("\t", $count);
        my $TE_posi=$a[0]."\t".$a[1]."\t".$a[2];
		#my $average_K9 = ($a[-2]+$a[-1])/2;
		if (exists $TE_enrich{$TE_posi}){
        $TE_K9_height{$TE_posi}=$TE_enrich{$TE_posi}."\t".$a[-2]."\t".$a[-1];
        }
}

my %anno;
$infile3="/dfs7/grylee/hgshukla/HiC_Project/".$name[$k]."/RNA_seq/annotation/".$name[$k].".gff";
open(FILE1,"<", "$infile3")||die"$!";
	while(my $count = <FILE1>){
		chomp($count);
		my @a = split("\t", $count);
        if ($a[2] eq "gene"){
        my @b = split('"', $a[8]);
        my $gene_id = $b[3];
        my $gene_posi=$a[0]."\t".$a[3]."\t".$a[4];
        $anno{$gene_id}=$gene_posi;
        }
}


#my $infile2="/pub/yuhenh3/ChIP-seq_2019/repeat_modeller_update/".$name[$k]."_".$p2."_gene_body_height_HMD_25_1kb_unique_mapped_10measures_update.txt";
my $infile2=$name[$k]."_ref_RSEM_rpkm_rank_z.txt";
open(FILE1,"<", "$infile2")||die"$!";
	while(my $count = <FILE1>){
		chomp($count);
		my @b = split("\t", $count);
        my @a = split('"', $b[0]);
		my $id=$a[1]; 

		if (exists $anno{$id}){
            #print "$id\n";
			my $TE_name="NA"; my $gene_dis=$win_leng; my $gene_num=0; 			
            my @g_posi=split("\t", $anno{$id});
            my $chr = $g_posi[0]; my $left = $g_posi[1]; my $right = $g_posi[2];
			foreach my $key (keys %TE_K9_height) {
			my @posi = split("\t", $key);
			my $dis;
        if ($chr eq $posi[0]){
		if ($posi[2] < $left){
			$dis=$left-$posi[2];
		}elsif($posi[1]>$right){
			$dis=$posi[1]-$right;
		}else{
			$dis=0;
		}
		if ($dis<$gene_dis){
			$gene_dis=$dis;
			$TE_name=$key;
			$gene_num++;
		}
}
}

#	if(($gene_dis>=$c)&&($gene_dis<$d)&&($gene_num>0)){
if ($TE_name ne "NA"){
	print OUT "$count\t","$anno{$id}\t","$TE_K9_height{$TE_name}\t","$gene_dis\t","$gene_num\n";
	#print OUT "$TE_enrich{$TE_name}\t","$id\t","$gene_posi\t","$hmd_1{$id}\t","$b[7]\t","$gene_dis\t","$gene_num\n";
}
}
}
