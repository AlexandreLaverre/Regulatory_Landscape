#!/usr/bin/perl

use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use strict;

#######################################################################

sub readTSSCoords{
    my $pathin=$_[0];
    my $baited=$_[1];
    my $tsscoords=$_[2];

    open(my $input, $pathin);

    my $line=<$input>; ## header
    chomp $line;
    
    my @s=split("\t", $line);
    my %header;
    
    for(my $i=0; $i<@s; $i++){
	$header{$s[$i]}=$i;
    }
    
    while($line){
	chomp $line;
	my @s=split("\t", $line);

	my $txid=$s[$header{"TranscriptID"}];
	
	if($baited eq "any" || exists $baited->{$txid}){
	    my $gene=$s[$header{"GeneID"}];
	    my $chr="chr".$s[$header{"Chr"}];
	    my $start=$s[$header{"Start"}]+0;
	    my $end=$s[$header{"End"}]+0;
	    my $strand=$s[$header{"Strand"}];
	    
	    my $tss=$start;
	    
	    if($strand eq "-1" || $strand eq "-"){
		$tss=$end;
	    }
	    
	    if(exists $tsscoords->{$gene}){
		$tsscoords->{$gene}{$tss}=1;
	    } else {
		$tsscoords->{$gene}={$tss=>1};
	    }
	}
	
	$line=<$input>;
    }
    
    close($input);
}

#######################################################################

sub readBaitedTranscripts{
    my $pathin=$_[0];
    my $baited=$_[1];
 
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

	my $tx=$s[$header{"transcript_ID"}];
	my @t=split(",", $tx);

	foreach my $txid (@t){
	    $baited->{$txid}=1;
	}
	
	$line=<$input>;
    }
    
    close($input);
}

#######################################################################

sub printHelp{

    my $parnames=$_[0];
    my $parvalues=$_[1];
    
    print "\n";
    print "This script computes the minimum distance between TSS and enhancers and replaces these columns in contact datasets. \n";
    print "\n";
    print "Options:\n";
    
    foreach my $par (@{$parnames}){
	print "--".$par."  [  default value: ".$parvalues->{$par}."  ]\n";
    }
    print "\n";
}

#########################################################################
#########################################################################

my %parameters;

$parameters{"pathContacts"}="NA";
$parameters{"pathReferenceBaitAnnotation"}="NA";
$parameters{"pathReferenceTranscriptCoordinates"}="NA";
$parameters{"pathTargetTranscriptCoordinates"}="NA";
$parameters{"pathOutput"}="NA";

my %defaultvalues;
my @defaultpars=("pathContacts", "pathReferenceBaitAnnotation", "pathReferenceTranscriptCoordinates", "pathTargetTranscriptCoordinates", "pathOutput");

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

print "Reading reference baited transcripts...\n";

my %refbaitedtx;
readBaitedTranscripts($parameters{"pathReferenceBaitAnnotation"}, \%refbaitedtx);

my $nbbaited=keys %refbaitedtx;

print "Found ".$nbbaited." baited transcripts.\n";
print "Done.\n";

######################################################################

print "Reading TSS coordinates...\n";

my %reftss; 
readTSSCoords($parameters{"pathReferenceTranscriptCoordinates"}, "any", \%reftss);

my %refbaitedtss;
readTSSCoords($parameters{"pathReferenceTranscriptCoordinates"}, \%refbaitedtx, \%refbaitedtss);

my %tgtss; 
readTSSCoords($parameters{"pathTargetTranscriptCoordinates"}, "any", \%tgtss);

print "Done.\n";

######################################################################

print "Reading contact data and writing output...\n";

open(my $input, $parameters{"pathContacts"});
open(my $output, ">".$parameters{"pathOutput"});

my $line=<$input>; ## header
chomp $line;

my @s=split("\t", $line);
my %header;

for(my $i=0; $i<@s; $i++){
    $header{$s[$i]}=$i;
}

my $lineout="";

my $indexorigindist="NA";
my $indextargetdist="NA";
    
for(my $i=0; $i<@s; $i++){
    my $field=$s[$i];
    
    if($field eq "origin_dist"){
	$indexorigindist=$i;
	$lineout.="origin_min_dist_allTSS\torigin_min_dist_baitedTSS\t";
    } else{
	if($field eq "target_dist"){
	    $indextargetdist=$i;
	    $lineout.="target_min_dist_allTSS\t";
	} else{
	    $lineout.=$field."\t";
	}
    }
}

chop $lineout;

print $output $lineout."\n";

if($indextargetdist ne "NA" && $indextargetdist<$indexorigindist){
    print "Weird! target distance appears first.\n";
    exit(1);
}

$line=<$input>;

while($line){
    chomp $line;
    my @s=split("\t", $line);

    ## distances for reference species

    my $origingene=$s[$header{"origin_gene"}];

    if(!exists $reftss{$origingene}){
	print "Weird! we don't have transcript coordinates for ".$origingene."\n";
	exit(1);
    }

     if(!exists $refbaitedtss{$origingene}){
	print "Weird! we don't have baited transcript coordinates for ".$origingene."\n";
	exit(1);
    }
    
    my $originenh=$s[$header{"origin_enh"}];
    my @t=split(":", $originenh);

    my $startenh=$t[1]+0;
    my $endenh=$t[2]+0;

    my @distances;
    
    foreach my $pos (keys %{$reftss{$origingene}}){
	my $dist=abs($pos-($startenh+$endenh)/2);
	push(@distances, $dist);
    }

    my $mindistallref=min @distances;

    my @baiteddistances;

    foreach my $pos (keys %{$refbaitedtss{$origingene}}){
	my $dist=abs($pos-($startenh+$endenh)/2);
	push(@baiteddistances, $dist);
    }

    my $mindistbaitedref=min @baiteddistances;

    ## target species, only if needed

    if($indextargetdist ne "NA"){
	my $targetgene=$s[$header{"target_gene"}];
	my $targetgenecoord=$s[$header{"target_gene_coord"}];
	my $oldtargetdist=$s[$header{"target_dist"}];

	my @u=split(":", $targetgenecoord);
	my $chrgene=$u[0];

	 if(!exists $tgtss{$targetgene}){
	     print "Weird! we don't have transcript coordinates for ".$targetgene."\n";
	     exit(1);
	 }
	 
	 my $tgenh=$s[$header{"target_enh"}];
	 my @t=split(":", $tgenh);
	 
	 my $chrenh=$t[0];
	 my $startenh=$t[1]+0;
	 my $endenh=$t[2]+0;
	 
	 my @distances;
	 
	 foreach my $pos (keys %{$tgtss{$targetgene}}){
	     my $dist=abs($pos-($startenh+$endenh)/2);
	     push(@distances, $dist);
	 }
	 
	my $mindistalltg=min @distances;

	if($oldtargetdist eq "trans"){ 
	    $mindistalltg="trans"; 
	}
	
	## indexorigindist should always be smaller than indextargetdist
	
	print $output join("\t", @s[0..($indexorigindist-1)])."\t".$mindistallref."\t".$mindistbaitedref."\t".join("\t", @s[($indexorigindist+1)..($indextargetdist-1)])."\t".$mindistalltg."\t".join("\t", @s[($indextargetdist+1)..$#s])."\n";
    	
    } else{
	## there is no target dist in this file
        print $output join("\t", @s[0..($indexorigindist-1)])."\t".$mindistallref."\t".$mindistbaitedref."\t".join("\t", @s[($indexorigindist+1)..$#s])."\n";
    }
    
    $line=<$input>;
}

close($input);
close($output);

print "Done.\n";

######################################################################
