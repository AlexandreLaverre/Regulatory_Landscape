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
    print "This script computes the minimum distance between TSS and enhancers. \n";
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
$parameters{"pathBaitAnnotation"}="NA";
$parameters{"pathTranscriptCoordinates"}="NA";
$parameters{"pathOutput"}="NA";

my %defaultvalues;
my @defaultpars=("pathContacts", "pathBaitAnnotation", "pathTranscriptCoordinates", "pathOutput");

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

print "Reading baited transcripts...\n";

my %baitedtx;
readBaitedTranscripts($parameters{"pathBaitAnnotation"}, \%baitedtx);

my $nbbaited=keys %baitedtx;

print "Found ".$nbbaited." baited transcripts.\n";
print "Done.\n";

######################################################################

print "Reading TSS coordinates...\n";

my %tss; 
readTSSCoords($parameters{"pathTranscriptCoordinates"}, "any", \%tss);

my %baitedtss;
readTSSCoords($parameters{"pathTranscriptCoordinates"}, \%baitedtx, \%baitedtss);

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

my $indexdist="NA";
    
for(my $i=0; $i<@s; $i++){
    my $field=$s[$i];
    
    if($field eq "dist"){
	$indexdist=$i;
	$lineout.="mean_dist_allTSS\tmin_dist_allTSS\tmin_dist_baitedTSS\t";
    } else{
	$lineout.=$field."\t";
    }
}

chop $lineout;

print $output $lineout."\n";

$line=<$input>;

while($line){
    chomp $line;
    my @s=split("\t", $line);

    my $gene=$s[$header{"gene"}];

    if(!exists $tss{$gene}){
	print "Weird! we don't have transcript coordinates for ".$gene."\n";
	exit(1);
    }

     if(!exists $baitedtss{$gene}){
	print "Weird! we don't have baited transcript coordinates for ".$gene."\n";
	exit(1);
    }
    
    my $enh=$s[$header{"enhancer"}];
    my @t=split(":", $enh);

    my $startenh=$t[1]+0;
    my $endenh=$t[2]+0;

    my @distances;
    
    foreach my $pos (keys %{$tss{$gene}}){
	my $dist=abs($pos-($startenh+$endenh)/2);
	push(@distances, $dist);
    }

    my $mindistall=min @distances;

    my @baiteddistances;

    foreach my $pos (keys %{$baitedtss{$gene}}){
	my $dist=abs($pos-($startenh+$endenh)/2);
	push(@baiteddistances, $dist);
    }

    my $mindistbaited=min @baiteddistances;
    
    print $output join("\t", @s[0..$indexdist])."\t".$mindistall."\t".$mindistbaited."\t".join("\t", @s[($indexdist+1)..$#s])."\n";
    
    $line=<$input>;
}

close($input);
close($output);

print "Done.\n";

######################################################################
