#!/usr/bin/perl

use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use strict;

#######################################################################

sub readFasta{
    my $path=$_[0];
    my $reffasta=$_[1];
   
    my @s=split("\\.",$path);
    my $ext=$s[-1];

    my $input;

    if($ext eq "gz"){
	open($input,"zcat $path |");
    }
    else{
	open($input, $path);
    }
    
    my $line=<$input>;

    while($line){
	my $b=substr $line,0,1;
	
	if($b eq ">"){
	    chomp $line;
	    my $id=substr $line,1;

	    my @s=split(" ",$id);
	    $id=$s[0];

	    # print "saw chromosome ".$id."\n";
	    
	    $reffasta->{$id}="";

	    $line=<$input>;
	    $b=substr $line,0,1;
	    
	    while($line && !($b eq ">")){
		chomp $line;
		$reffasta->{$id}.=$line;
		$line=<$input>;
		$b=substr $line,0,1;
	    }
	}
    }

    close($input);
}

#########################################################################################

sub computeGCContent{
    my $seq=$_[0];

    my $nbA = ($seq =~ tr/A//);
    my $nbC = ($seq =~ tr/C//);
    my $nbG = ($seq =~ tr/G//);
    my $nbT = ($seq =~ tr/T//);

    my $GC="NA";

    my $nbtot=$nbA+$nbC+$nbG+$nbT;
    
    if($nbtot>0){
	$GC=($nbG+$nbC)/$nbtot;
    }

    return $GC;
}

#########################################################################################

sub printHelp{

    my $parnames=$_[0];
    my $parvalues=$_[1];
    
    print "\n";
    print "This script computes GC content for genomic sequences. \n";
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

$parameters{"pathGenomeSequence"}="NA";
$parameters{"pathCoordinates"}="NA";
$parameters{"pathOutput"}="NA";

my %defaultvalues;
my @defaultpars=("pathGenomeSequence", "pathCoordinates", "pathOutput");

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

#####################################################################################
#####################################################################################

# my $gc1=computeGCContent("ACGT");
# print "0.5 = ".$gc1."\n";

# my $gc2=computeGCContent("ACGGGT");
# print "0.66 = ".$gc2."\n";

# my $gc3=computeGCContent("ATTTT");
# print "0 = ".$gc3."\n";

# my $gc4=computeGCContent("GGCCGG");
# print "1 = ".$gc4."\n";

#####################################################################################

print "Reading genome sequence...\n";

my %genome;
readFasta($parameters{"pathGenomeSequence"}, \%genome);

print "Done.\n";

#####################################################################################

print "Reading input coordinates and writing output for GC content...\n";

open(my $input, $parameters{"pathCoordinates"});
open(my $output, ">".$parameters{"pathOutput"});

my $line=<$input>; ## header
chomp $line;

my @s=split("\t", $line);
my %header;

my $lineout="";

for(my $i=0; $i<@s; $i++){
    if($s[$i] eq "GC_bp"){
	$lineout.="GC_bp_norepeat\t";
    } else{
	$lineout.=$s[$i]."\t";
    }
    
    $header{$s[$i]}=$i;
}

$lineout.="GC_content"; ## we have a trailing \t
print $output $lineout."\n";

$line=<$input>;

while($line){
    chomp $line;
    my @s=split("\t", $line);
    
    my $chr=$s[$header{"chr"}];
    my $start=$s[$header{"start"}]+0;
    my $end=$s[$header{"end"}]+0;

    my $gc="NA";

    if(exists $genome{$chr}){
	my $sequence=substr $genome{$chr}, $start-1, ($end-$start+1);
	$sequence=uc $sequence;
	$gc=computeGCContent($sequence);
    } else{
	my $prefix=substr $chr, 0, 3;

	if($prefix eq "chr"){
	    $chr=substr $chr, 3;
	    
	    if(exists $genome{$chr}){
		my $sequence=substr $genome{$chr}, $start-1, ($end-$start+1);
		$sequence=uc $sequence;
		$gc=computeGCContent($sequence);
	    } else{
		print "cannot find ".$chr." in genome sequence.\n";
		exit(1);
	    }
	} else{
	    print "cannot find ".$chr." in genome sequence.\n";
	    exit(1);
	}
    }

    print $output $line."\t".$gc."\n";
    
    $line=<$input>;
}

close($output);
close($input);

print "Done.\n";

#####################################################################################
