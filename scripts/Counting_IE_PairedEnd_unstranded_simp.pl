#!/usr/bin/perl


use warnings;

my $INTRONS = shift; # Intron Bed Files
my $EXONS = shift; # Exon Bed Files
my $bam = shift;  # Bam Input
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
    $s_e{$gene_n}{'e'} = $exon_e;
    }
}else{
$s_e{$gene_n}{'s'}=$exon_s;
$s_e{$gene_n}{'e'}=$exon_e;
}
my $ss=$s_e{$gene_n}{'s'};
my $ee=$s_e{$gene_n}{'e'};
$s_e{$gene_n}{'comp'}= "$chr\:$ss\-$ee";
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
foreach my $key1 (reverse sort keys %gene_e){
#if(($index_g >= $b1) & ($index_g <= $b2)){
my %read_e = ();

my @IN2="";
my $length_e = 0; 
@IN2 =`samtools view $bam $s_e{$key1}{'comp'}`;
$key2 =~ m/\w+\:(\d+)\-(\d+)/;
$length_e =  $2 - $1;
foreach(@IN2){
	my @c = split;
		if(($c[11] eq "NH:i:1") & ($c[3] > $s_e{$key1}{'s'}) & ($c[3] < $s_e{$key1}{'e'})){
		$read_e{$c[0]}{$c[1]}  = 1;
		}
	}


my %read_i = ();
my $length_i = 0;

	foreach my $key2 (sort keys %{$gene_i{$key1}}){
    my @IN2="";
	@IN2 =`samtools view $bam $key2`;
	$key2 =~ m/\w+\:(\d+)\-(\d+)/;
    $length_i = $length_i + $2 - $1;
    	foreach(@IN2){
        my @c = split;
        my $l = 0;
        $c[5] =~ m/(\d+)M/;
		my $le = $1;
			if(($c[5] !~ /N/) & ($le > 70) & ($c[11] eq 'NH:i:1')){
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
my $l_e=$length_e - $length_i;
print("$key1\t$nb_keys_e\t$nb_keys_i\t$l_e\t$length_i\n");

#}
#$index_g=$index_g + 1;
}















