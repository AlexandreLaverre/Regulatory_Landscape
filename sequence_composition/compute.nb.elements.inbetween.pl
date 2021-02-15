#!/usr/bin/perl

use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use strict;

#######################################################################

sub readCoordinates{
    my $pathin=$_[0];
    my $type=$_[1];
    my $coords=$_[2];
   
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

	    my $origin_gene_coord=$s[$header{"origin_gene_coord"}];
	    my @t=split(":", $origin_gene_coord);
	    
	    my $chr_gene=$t[0];
	    my $start_gene=$t[1]+0;
	    my $end_gene=$t[2]+0;

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
	    
	    $start=min($start_gene, $start_enh);
	    $end=max($end_gene, $end_enh);
	
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

			push(@{$intersect->{$id1}}, $id2);
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
	    }
	}
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

$parameters{"pathCoordinates1"}="NA";
$parameters{"pathCoordinates2"}="NA";
$parameters{"type"}="NA";
$parameters{"pathOutput"}="NA";

my %defaultvalues;
my @defaultpars=("pathCoordinates1", "pathCoordinates2","type", "pathOutput");

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

print "Reading interacting element coordinates...\n";

my $type=$parameters{"type"};

print "type ".$type."\n";

my %fragments;
readCoordinates($parameters{"pathCoordinates1"}, "fragments", \%fragments);

print "Done.\n";

print "Reading ".$type." coordinates...\n";

my %baits;
readCoordinates($parameters{"pathCoordinates2"}, $type, \%baits);

print "Done.\n";

######################################################################

print "Computing intersections...\n";

my %intersections;

my $margin=0;

print "margin ".$margin."\n";

$intersections{$margin}={};

computeIntersection(\%fragments, \%baits, $margin, $intersections{$margin});

print "Done.\n";

######################################################################

print "Writing output...\n";

open(my $input, $parameters{"pathCoordinates1"});
open(my $output, ">".$parameters{"pathOutput"});

my $line=<$input>;
chomp $line;

my @s=split("\t", $line);

my %header;

for(my $i=0; $i<@s; $i++){
    $header{$s[$i]}=$i;
}

print $output $line."\tnb_".$type."_inbetween\n";

$line=<$input>;

while($line){
    chomp $line;
    my @s=split("\t", $line);

    my $origin_gene_coord=$s[$header{"origin_gene_coord"}];
    my @t=split(":", $origin_gene_coord);
    
    my $chr_gene=$t[0];
    my $start_gene=$t[1]+0;
    my $end_gene=$t[2]+0;
    
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
    
    my $start=min($start_gene, $start_enh);
    my $end=max($end_gene, $end_enh);
    
    my $lineout=$line;
    my $margin=0;
    
    my $id=$chr.":".$start."-".$end;
    my $nb=@{$intersections{$margin}{$id}}; ## nb elements in this window
    
    $lineout.="\t".$nb;
    
    print $output $lineout."\n";
    
    $line=<$input>;
}

close($output);

print "Done.\n";

######################################################################
