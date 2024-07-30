#!/usr/bin/perl


use warnings;


my $INTRONS = shift; # Intron Bed Files
my $EXONS = shift; # Exon Bed Files
my $bam = shift;  # Bam Input  # Bam Input
#my $b1 = shift; # Genes bottom index for // computing
#my $b2 = shift; # Genes up index for // computing 

open IN1,$EXONS or die;

my %gene_e;
my %s_e;
#### Load Exons information from bed files 
while(<IN1>){

 my @a = split;
 my $exon_s = $a[1];
 my $exon_e = $a[2];
 my $chr = $a[0];
 my $gene_n = $a[3];
 my $comp = "$chr\:$exon_s\-$exon_e";
 my $strand = $a[5];
 $gene_e{$gene_n}{$comp} = $strand;
if(exists($s_e{$gene_n})){
	if($exon_s < $s_e{$gene_n}{'s'}){
	$s_e{$gene_n}{'s'} = $exon_s;
	}
	if($exon_e > $s_e{$gene_n}{'e'}){
        $s_e{$gene_n}{'e'} = $exon_s;
        }

}else{
$s_e{$gene_n}{'s'}=$exon_s;
$s_e{$gene_n}{'e'}=$exon_e;
}
}
close IN1;

#### Load Introns information from bed files 
open IN2,$INTRONS or die;

my %gene_i;


while(<IN2>){

 my @a = split;

 my $intron_s = $a[1] + 1;
 my $intron_e = $a[2] - 1;
 my $chr = $a[0];
 my $gene_n = $a[3];
 my $comp = "$chr\:$intron_s\-$intron_e";
 my $strand = $a[5];
 $gene_i{$gene_n}{$comp} = $strand;
}
                
close IN2;               
                                
#######################################

my $strand_1 = "";
my $strand_2 = "";
my $index_g=1;
foreach my $key1 (reverse sort keys %gene_e)
{
#if(($index_g >= $b1) & ($index_g <= $b2)){
	my %read_e = ();
	my $length_i = 0;
	my $length_e = 0;
	foreach my $key2 (sort keys %{$gene_e{$key1}})
    {
    	my @IN2="";

 		if($gene_e{$key1}{$key2} eq '-') # A changer en + si zero count  
		{
			$strand_1 = 99;  # read paired, read mapped in proper pair, mate reverse strand, first  in pair
			$strand_2 = 147; # read paired, read mapped in proper pair, read reverse strand, second in pair
		}
		else
		{
			$strand_1 = 83;  # read paired, read mapped in proper pair, read reverse strand, first  in pair
			$strand_2 = 163; # read paired, read mapped in proper pair, mate reverse strand, second in pair
		}
		@IN2 =`samtools view  $bam $key2`;
		$key2 =~ m/\w+\:(\d+)\-(\d+)/;
		$length_e = $length_e + $2 -$1;
		foreach(@IN2){
			my @c = split;
			if(($c[11] eq "NH:i:1") & ($c[3] > $s_e{$key1}{'s'}) & ($c[3] < $s_e{$key1}{'e'}) & (($c[1] == $strand_1) | ($c[1] == $strand_2)) ){
				$read_e{$c[0]}{$c[1]} = 1;
			}
        }
	}

	my %read_i = ();
	foreach my $key2 (sort keys %{$gene_i{$key1}})
    {
	
    	my @IN2="";


		@IN2 =`samtools view  $bam $key2`;
		$key2 =~ m/\w+\:(\d+)\-(\d+)/;
		$length_i = $length_i + $2 -$1;
        foreach(@IN2){
            my @c = split;
			my $l = 0;
			$c[5] =~ m/(\d+)M/;
		 	my $le = $1;
			if(($c[5] !~ /N/) & ($le > 95) & ($c[11] eq 'NH:i:1') & (($c[1] == $strand_1) | ($c[1] == $strand_2)) ){
                $read_i{$c[0]}{$c[1]}  = 1;
			}
        }
    }

	foreach my $kk (sort keys %read_e)
	{
		unless((exists($read_e{$kk}{99}) & exists($read_e{$kk}{147})) or (exists($read_e{$kk}{83}) & exists($read_e{$kk}{163}))){
			delete $read_e{$kk};
		}
	
		if(exists($read_i{$kk})){
			#delete $read_e{$kk}{$strand_1};
			#delete $read_e{$kk}{$strand_2};
			delete $read_e{$kk};
		}
	}
	my $nb_keys_e = keys %read_e;
	my $nb_keys_i = keys %read_i;
	print("$key1\t$nb_keys_e\t$nb_keys_i\t$length_e\t$length_i\n");

#}
#$index_g=$index_g + 1;
}















