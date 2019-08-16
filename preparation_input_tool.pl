#!/usr/bin/perl

#open the files
$num_args = $#ARGV + 1;
if ($num_args != 4) {
    print "\nUsage: preparation_input_tool.pl file_orthologies file_int_species file_ref_species output\n";
    exit;
}
$orthofile=$ARGV[0];
$intfile=$ARGV[1];
$reffile=$ARGV[2];
$outfile=$ARGV[3];
open(IN1, $orthofile) or die 'Fail in Orthology file';
open(IN2, $intfile) or die 'Fail in Int Species file';
open(IN3, $reffile) or die 'Fail in Ref Species file';
open(OUT1, '>',$outfile) or die 'Could not create output file';
$i=0;
$familiy=0;

#Set the orthologic relations of every ref gene to each int gene
while (<IN1>){
    @colum= split /\t/, $_;
    if ($family == $colum[0]){
            if ($colum[1] eq "Hsa"){
            $hsatablort[$i]=$colum[2];
            $i++;
        }
            else{
                $hashort{$colum[2]}="@hsatablort";
            }
    }
    else{
    $i=0;
    undef @hsatablort;
    if ($colum[1] eq "Hsa"){
            $hsatablort[$i]=$colum[2];
            $i++;
        }
    $family=$colum[0];
    }

}
$i=1;
close IN1;
undef @colum;
#Prepare a hash with the position, chromosome and strand for every Gene_ID 
while (<IN2>){
    @colum= split /\t/, $_;
    @subcolum= split /"/, $colum[8];
    if ($colum[0]!= $currentcrom){
        $i=1;
        $currentcrom=$colum[0];
    }
    $sign= "-" if $colum[6] eq "-";
    $sign= "+" if $colum[6] eq "+";
	$chromo=$colum[0];
	$chromo= 1000 if $colum[0] eq "Y";
	$chromo= 0 if $colum[0] eq "X";
    $hashpos{$subcolum[1]}="$i|$chromo|$sign";
    $i++;
}
$i=1;
close IN2;
undef @colum;
undef @hsatablort;
#Add To every Ref gene, its correspondent Int Orthologs
while(<IN3>){
    $_ =~ s/\r|\n//g;
    @colum= split /,/, $_;
    if ($colum[0] ne $currentcrom){
        $i=1;
        $currentcrom=$colum[0];
    }
    @hsatablort= split /\s/, $hashort{$colum[3]};
    $sign= -1 if $colum[6] eq "-";
    $sign= 1 if $colum[6] eq "+";
    print OUT1 "$colum[3]";
    $i++;
    print OUT1 "\tHypnc_$colum[3]" if $colum[6] eq "linc";
    for ($j=0; $j <= @hsatablort; $j++){
        print OUT1 "\t$hsatablort[$j]|$hashpos{$hsatablort[$j]}" if exists $hashpos{$hsatablort[$j]};
    }
    print OUT1 "\n";
    
}
close OUT1;
