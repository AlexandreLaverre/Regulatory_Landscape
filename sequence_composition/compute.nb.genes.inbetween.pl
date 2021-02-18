#!/usr/bin/perl

use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use strict;

#######################################################################

sub readGeneCoords{
    my $pathin=$_[0];
    my $genecoords=$_[1];

    open(my $input, $pathin);
    
    my $line=<$input>; ## header
    chomp $line;
    
    my @s=split("\t", $line);
    my %header;
    
    for(my $i=0; $i<@s; $i++){
	$header{$s[$i]}=$i;
    }

    $line=<$input>;

    while($line){
	chomp $line;
	my @s=split("\t", $line);

	my $gene=$s[$header{"GeneID"}];
	my $chr="chr".$s[$header{"Chr"}];
	my $start=$s[$header{"Start"}]+0;
	my $end=$s[$header{"End"}]+0;
	my $strand=$s[$header{"Strand"}];

	my $tss=$start;

	if($strand eq "-1" || $strand eq "-"){
	    $tss=$end;
	}

	$genecoords->{$gene}={"chr"=>$chr, "start"=>$start, "end"=>$end, "strand"=>$strand, "tss"=>$tss};
	
	$line=<$input>;
    }

    close($input);
}

#######################################################################

sub readCoordinates{
    my $pathin=$_[0];
    my $type=$_[1];
    my $genes=$_[2];
    my $coords=$_[3];
   
    open(my $input, $pathin);
    
    my $line=<$input>; ## header
    chomp $line;
    
    my @s=split("\t", $line);
    my %header;
        
    for(my $i=0; $i<@s; $i++){
	$header{$s[$i]}=$i;
    }

    $line=<$input>;

    my %unordered;

    while($line){
	chomp $line;

	my @s=split("\t", $line);
	
	my $chr="NA"; 
	my $start="NA"; 
	my $end="NA"; 

	if($type eq "genes"){
	    my $biotype=$s[$header{"GeneBiotype"}];

	    if($biotype ne "protein_coding"){
		$line=<$input>;
		next;
	    }
	    
	    $chr="chr".$s[$header{"Chr"}];
	    $start=$s[$header{"Start"}]+0;
	    $end=$s[$header{"End"}]+0;
	} else{
	    ## interactions
	    my $origin_gene=$s[$header{"origin_gene"}];

	    if(!exists $genes->{$origin_gene}){
		print "Weird! cannot find ".$origin_gene." in gene coordinates.\n";
		exit(1);
	    }

	    my $chr_gene=$genes->{$origin_gene}{"chr"};
	    my $tss_gene=$genes->{$origin_gene}{"tss"};
	   
	    my $origin_enh=$s[$header{"origin_enh"}];
	    my @u=split(":", $origin_enh);
	    
	    my $chr_enh=$u[0];
	    my $start_enh=$u[1]+0;
	    my $end_enh=$u[2]+0;

	    if($chr_gene ne $chr_enh){
		print "weird! different chromosomes found in ".$line."\n";
		exit(1);
	    }

	    $chr=$chr_enh;
	    
	    $start=min($tss_gene, $start_enh);
	    $end=max($tss_gene, $end_enh);
	}
	
	if($chr ne "NA"){
	    if(exists $unordered{$chr}){
		if(exists $unordered{$chr}{$start}){
		    push(@{$unordered{$chr}{$start}}, $end);
		} else{
		    $unordered{$chr}{$start}=[$end];
		}
	    } else{
		$unordered{$chr}={$start=>[$end]};
	    }
	}
	
	$line=<$input>;
    }

    foreach my $chr (keys %unordered){
	$coords->{$chr}={"start"=>[], "end"=>[]};

	my @uniquestart=keys %{$unordered{$chr}};
	my @sortedstart=sort {$a<=>$b} @uniquestart;

	foreach my $start (@sortedstart){
	    foreach my $end (sort {$a<=>$b} @{$unordered{$chr}{$start}}){
		push(@{$coords->{$chr}{"start"}}, $start);
		push(@{$coords->{$chr}{"end"}}, $end);
	    }
	}
    }
    
    close($input);
}

#######################################################################

sub computeIntersection{
    my $coords1=$_[0];
    my $coords2=$_[1];
    my $margin=$_[2];
    my $intersect=$_[3];
    my $intersectcoords=$_[4];

    foreach my $chr (keys %{$coords1}){
	my $nb1=@{$coords1->{$chr}{"start"}};
	
	if(exists $coords2->{$chr}){
	    my $nb2=@{$coords2->{$chr}{"start"}};

	    my $firstj=0;

	    for(my $i=0; $i<$nb1; $i++){
		my $start1=${$coords1->{$chr}{"start"}}[$i];
		my $end1=${$coords1->{$chr}{"end"}}[$i];

		my $id1=$chr.":".$start1."-".$end1;
		
		$intersect->{$id1}=[];
		$intersectcoords->{$id1}=[];

		my $j=$firstj;
		
		while($j<$nb2 && ${$coords2->{$chr}{"end"}}[$j]<($start1-$margin)){
		    $j++;
		}
		
		$firstj=$j;
		
		while($j<$nb2 && ${$coords2->{$chr}{"start"}}[$j]<($end1+$margin)){
		    my $start2=${$coords2->{$chr}{"start"}}[$j];
		    my $end2=${$coords2->{$chr}{"end"}}[$j];
		    
		    my $M=max($start1, ($start2-$margin));
		    my $m=min($end1, ($end2+$margin));

		    if($M<=$m){
			my $id2=$chr.":".$start2."-".$end2;
			my $idoverlap=$chr.":".$M."-".$m;

			push(@{$intersect->{$id1}}, $id2);
			push(@{$intersectcoords->{$id1}}, $idoverlap);
		    }
		    
		    $j++;
		}
	    }
	    
	} else{
	    for(my $i=0; $i<$nb1; $i++){
		my $start1=${$coords1->{$chr}{"start"}}[$i];
		my $end1=${$coords1->{$chr}{"end"}}[$i];
		my $id1=$chr.":".$start1."-".$end1;

		$intersect->{$id1}=[];
		$intersectcoords->{$id1}=[];
	    }
	}
    }
}

#######################################################################

sub computeIntersectionBlocks{
    my $intersect=$_[0];
    my $blocks=$_[1];
    
    foreach my $id (keys %{$intersect}){
	my %hashcoords;

	foreach my $id2 (@{$intersect->{$id}}){
	    my @s=split(":", $id2);
	    my @t=split("-", $s[1]);
	    
	    my $start=$t[0]+0;
	    my $end=$t[1]+0;

	    if(exists $hashcoords{$start}){
		if($end>$hashcoords{$start}){
		    $hashcoords{$start}=$end;
		}
	    } else{
		$hashcoords{$start}=$end;
	    }
	}

	$blocks->{$id}={"start"=>[], "end"=>[], "length"=>0};

	my @uniquestart=keys %hashcoords;
	my @sortedstart=sort {$a<=>$b} @uniquestart;

	my $currentstart=$sortedstart[0];
	my $currentend=$hashcoords{$currentstart};

	for(my $i=1; $i<@sortedstart; $i++){
	    my $thisstart=$sortedstart[$i];
	    my $thisend=$hashcoords{$thisstart};

	    if($thisstart>=$currentstart && $thisstart<=($currentend+1)){
		if($thisend>$currentend){
		    $currentend=$thisend;
		}
	    } else{
		push(@{$blocks->{$id}{"start"}}, $currentstart);
		push(@{$blocks->{$id}{"end"}}, $currentend);
		$blocks->{$id}{"length"}+=($currentend-$currentstart+1);

		$currentstart=$thisstart;
		$currentend=$thisend;
	    }
	}

	push(@{$blocks->{$id}{"start"}}, $currentstart);
	push(@{$blocks->{$id}{"end"}}, $currentend);
	$blocks->{$id}{"length"}+=($currentend-$currentstart+1);
    }
}

#######################################################################

sub printHelp{

    my $parnames=$_[0];
    my $parvalues=$_[1];
    
    print "\n";
    print "This script computes numbers of baits or genes between two fragments in contact. \n";
    print "\n";
    print "Options:\n";
    
    foreach my $par (@{$parnames}){
	print "--".$par."  [  default value: ".$parvalues->{$par}."  ]\n";
    }
    print "\n";
}

##########################################################################################
##########################################################################################

my %parameters;

$parameters{"pathContacts"}="NA";
$parameters{"pathReferenceGenes"}="NA";
$parameters{"pathTargetGenes"}="NA";
$parameters{"pathOutput"}="NA";

my %defaultvalues;
my @defaultpars=("pathContacts", "pathReferenceGenes", "pathTargetGenes", "pathOutput");

my %numericpars;

foreach my $par (keys %parameters){
    $defaultvalues{$par}=$parameters{$par};
}

## check if help was asked 

foreach my $arg (@ARGV){
    if($arg eq "--help"){
	printHelp(\@defaultpars, \%defaultvalues);
	exit(0);
    }
}

## check new parameters

my $nbargs=@ARGV;

for(my $i=0; $i<$nbargs; $i++){
    my $arg=$ARGV[$i];
    $arg=substr $arg,2;
    my @s=split("=",$arg);
    my $parname=$s[0];
    my $parval=$s[1];
    
    if(exists $parameters{$parname}){
	$parameters{$parname}=$parval;
	
	if(exists $numericpars{$parname}){
	    $parameters{$parname}=$parval+0.0;
	}
    }
    else{
	print "Error: parameter ".$parname." was not recognized!!!\n";
	printHelp(\@defaultpars, \%defaultvalues);
	exit(1);
    }
}

## show parameters

print "\n";

print "Running program with the following parameters:\n";

foreach my $par (@defaultpars){
    print "--".$par."=".$parameters{$par}."\n";
}

print "\n";

######################################################################
######################################################################

print "Reading gene coordinates...\n";

my %refgenes;
readGeneCoords($parameters{"pathReferenceGenes"}, \%refgenes);

my %tggenes;
readGeneCoords($parameters{"pathTargetGenes"}, \%tggenes);

print "Done.\n";

######################################################################

print "Reading interacting element coordinates...\n";

my %coords1;
readCoordinates($parameters{"pathContacts"}, "contacts", \%refgenes, \%coords1);

my %coords2;
readCoordinates($parameters{"pathReferenceGenes"}, "genes", \%refgenes, \%coords2);

print "Done.\n";

######################################################################

print "Computing intersections...\n";

my %intersections;
my %intersectcoords;

my $margin=0;

print "margin ".$margin."\n";

$intersections{$margin}={};
$intersectcoords{$margin}={};

computeIntersection(\%coords1, \%coords2, $margin, $intersections{$margin},$intersectcoords{$margin});

print "Done.\n";

######################################################################

print "Constructing intersection blocks...\n";

my %blocks;

foreach my $margin (keys %intersectcoords){
    $blocks{$margin}={};
    computeIntersectionBlocks($intersectcoords{$margin}, $blocks{$margin});
}

print "Done.\n";

######################################################################

print "Writing output...\n";

open(my $input, $parameters{"pathContacts"});
open(my $output, ">".$parameters{"pathOutput"});

my $line=<$input>;
chomp $line;

my @s=split("\t", $line);

my %header;

for(my $i=0; $i<@s; $i++){
    $header{$s[$i]}=$i;
}

print $output "origin_gene\torigin_gene_coord\torigin_enh\torigin_dist\ttarget_gene\ttarget_gene_coord\ttarget_enh\ttarget_dist\tnb_genes_inbetween\tfr_length_genes_inbetween\n";

$line=<$input>;

while($line){
    chomp $line;
    my @s=split("\t", $line);
    
    my $origin_gene=$s[$header{"origin_gene"}];
    
    if(!exists $refgenes{$origin_gene}){
	print "Weird! cannot find ".$origin_gene." in gene coordinates.\n";
	exit(1);
    }
    
    my $chr_gene=$refgenes{$origin_gene}{"chr"};
    my $start_gene=$refgenes{$origin_gene}{"start"};
    my $end_gene=$refgenes{$origin_gene}{"end"};
    my $strand_gene=$refgenes{$origin_gene}{"strand"};
    my $tss_gene=$refgenes{$origin_gene}{"tss"};
    
    my $origin_enh=$s[$header{"origin_enh"}];
    my @u=split(":", $origin_enh);
	    
    my $chr_enh=$u[0];
    my $start_enh=$u[1]+0;
    my $end_enh=$u[2]+0;
    
    if($chr_gene ne $chr_enh){
	print "weird! different chromosomes found in ".$line."\n";
	exit(1);
    }
    
    my $chr=$chr_enh;
    
    my $start=min($tss_gene, $start_enh);
    my $end=max($tss_gene, $end_enh);
            
    my $margin=0;
    
    my $id=$chr.":".$start."-".$end;
    my $nb=@{$intersections{$margin}{$id}}; ## nb elements in this window
    my $length=$blocks{$margin}{$id}{"length"};

    my $frlength=$length/($end-$start+1);

    ## modify output slightly

    my $origin_gene_coord=$chr_gene.":".$start_gene.":".$end_gene.":".$strand_gene;
    my $origin_dist=$s[$header{"origin_dist"}];
    my $target_gene=$s[$header{"target_gene"}];
    my $target_enh=$s[$header{"target_enh"}];
    my $target_dist=$s[$header{"target_dist"}];

    if(!exists $tggenes{$target_gene}){
	print "Weird! cannot find ".$target_gene." in gene coordinates.\n";
	exit(1);
    }

    my $chr_gene_tg=$tggenes{$target_gene}{"chr"};
    my $start_gene_tg=$tggenes{$target_gene}{"start"};
    my $end_gene_tg=$tggenes{$target_gene}{"end"};
    my $strand_gene_tg=$tggenes{$target_gene}{"strand"};
    my $target_gene_coord=$chr_gene_tg.":".$start_gene_tg.":".$end_gene_tg.":".$strand_gene_tg;
    
    my $lineout=$origin_gene."\t".$origin_gene_coord."\t".$origin_enh."\t".$origin_dist."\t".$target_gene."\t".$target_gene_coord."\t".$target_enh."\t".$target_dist."\t".$nb."\t".$frlength;
    
    print $output $lineout."\n";
    
    $line=<$input>;
}

close($output);

print "Done.\n";

######################################################################
