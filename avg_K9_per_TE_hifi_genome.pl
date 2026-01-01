#!/usr/bin/perl
use strict;
use warnings;
use List::Util qw( min max );
my $measure=50; my $win_leng=5000;
#my $measure=10; my $win_leng=1000;
my $win_size=5000; my $control_win_size=$win_size*2; my $dis1=30000;my $dis2=40000;

my $name="A4";
my @r=("1","2");
#my @r=("1","2","3");

my $outfile1="/dfs7/grylee/yuhenh3/HC_HiFi/".$name."_mean_K9_5kb_TE_flanking_left_right.txt"; 
#my $outfile1="/dfs7/grylee/yuhenh3/HC_HiFi/".$name."_Zita_mean_K9_5kb_TE_flanking_left_right.txt"; 

open(OUT, ">$outfile1");

my %hmd_left;my %hmd_right;
for (my $i=0; $i<@r; $i+=1){
    my $p= $r[$i];

my $infile1="/dfs7/grylee/yuhenh3/HC_HiFi/focal_".$name."_HiFi_".$p."_height_HMD_25_".$win_leng."_".$measure."_measures_mean.txt"; 
#my $infile1="/dfs7/grylee/yuhenh3/CUT_tag_experiment/Zita_A4_A7_CUT_tag/".$name."_MA_parental_".$p."_height_HMD_25_5000_50_measures_mean.txt"; 

open(FILE1,"<", "$infile1")||die"$!";
	while(my $count = <FILE1>){
		chomp($count);
		my @b = split("\t", $count);
		if ($b[0] ne "Chr"){
my $locat=$b[0]."\t".$b[1]."\t".$b[2];
my $left=$b[7]; my $right=$b[8]; #for mass

if (exists $hmd_left{$locat}){
$hmd_left{$locat}=$hmd_left{$locat}."\t".$left;
}else{
	$hmd_left{$locat}=$left;
}
if (exists $hmd_right{$locat}){
$hmd_right{$locat}=$hmd_right{$locat}."\t".$right;
}else{
	$hmd_right{$locat}=$right;
}
}
}
}

my $infile2="/dfs7/grylee/hgshukla/HiC_Project/A4/sample1/UU_minmapq_0/code_process_UU_MU_parallel/LNorm/UU_MU/A4_A7.UU_MU.Zscore.SigFlag.txt";

#my $infile2= "/dfs7/grylee/hgshukla/HiC_Project/A7/sample1/code_process_MU_UU_Parallel/LNorm/A7_A4.combined.LNorm.Zscore.txt";
open(FILE1,"<", "$infile2")||die"$!";
	while(my $count = <FILE1>){
		chomp($count);
		my @b = split("\t", $count);
my $chr=$b[0]; my $star=$b[1]; my $end=$b[2];
my $TE= $chr."\t".$star."\t".$end;

if (exists $hmd_left{$TE}){
    my @left=split("\t", $hmd_left{$TE});my @right=split("\t", $hmd_right{$TE});
    my $left_mean=mean(@left); my $right_mean=mean(@right);
print OUT "$count\t","$left_mean\t","$right_mean\n"; #for betw strain comparisons
#print OUT "$count\t","alternative_K9\t","$hmd{$TE}\n"; #for betw strain comparisons
}
}




sub mean {
my @array = @_; # save the array passed to this function
my $sum; # create a variable to hold the sum of the array's values
foreach (@array) { $sum += $_; } # add each element of the array
# to the sum
return $sum/@array; # divide sum by the number of elements in the
# array to find the mean
}


sub median
{
    my @vals = sort {$a <=> $b} @_;
    my $len = @vals;
    if($len%2) #odd?
    {
        return $vals[int($len/2)];
    }
    else #even
    {
        return ($vals[int($len/2)-1] + $vals[int($len/2)])/2;
    }
}
