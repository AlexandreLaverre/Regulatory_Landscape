use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use strict;

####################################################################################
####################################################################################

sub readCoverage{
    my $pathin=$_[0];
    my $sample=$_[1];
    my $combinedcoverage=$_[2];

    open(my $input, $pathin);
    my %header;
    my $line=<$input>;
    chomp $line;
    my @s=split("\t", $line);

    for(my $i=0; $i<@s; $i++){
	$header{$s[$i]}=$i;
    }

    $line=<$input>;

    while($line){
	chomp $line;
	my @s=split("\t", $line);

	my $chr=$s[$header{"Chr"}];
	my $start=$s[$header{"Start"}];
	my $end=$s[$header{"End"}];
	my $coverage=$s[$header{"Coverage"}]+0.0;

	my $id=$chr.",".$start.",".$end;

	if(exists $combinedcoverage->{$id}){
	    $combinedcoverage->{$id}{$sample}=$coverage; 
	} else{
	    $combinedcoverage->{$id}={$sample=>$coverage}; 
	}
	
	$line=<$input>;
    }

    close($input);
}

####################################################################################

sub printHelp{
    my $parnames=$_[0];
    my $parvalues=$_[1];
    
    print "\n";
    print "This script combines read coverage for multiple samples.\n";
    print "\n";
    
    print "Options:\n";
    
    foreach my $par (@{$parnames}){
	print "--".$par."  [  default value: ".$parvalues->{$par}."  ]\n";
    }

    print "\n";
}

####################################################################################
####################################################################################

my %parameters;

$parameters{"pathsCoverage"}="NA";
$parameters{"samples"}="NA";
$parameters{"pathOutput"}="NA";

my @defaultpars=("pathsCoverage", "samples", "pathOutput");

my %defaultvalues;

foreach my $par (keys %parameters){
    $defaultvalues{$par}=$parameters{$par};
}

## update arguments

my $nbargs=@ARGV;

for(my $i=0;$i<$nbargs; $i++){
    my $arg=$ARGV[$i];
    $arg=substr $arg,2;
    
    my @s=split("=",$arg);
    my $parname=$s[0];
    my $parval=$s[1];
    
    if(exists $parameters{$parname}){
	$parameters{$parname}=$parval;
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

####################################################################################
####################################################################################

print "Reading coverage for each sample...\n";

my %coverage;

my @paths=split(",", $parameters{"pathsCoverage"});
my @samples=split(",", $parameters{"samples"});

my $nbp=@paths;
my $nbs=@samples;

if($nbp!=$nbs){
    print "Saw ".$nbp." paths and ".$nbs." samples!\n";
    exit(1);
}

for(my $i=0; $i<$nbp; $i++){
    my $thispath=$paths[$i];
    my $thissample=$samples[$i];

    print "Reading coverage from ".$thispath." for ".$thissample."\n";
    
    readCoverage($thispath, $thissample, \%coverage);
}

print "Done.\n";

####################################################################################
####################################################################################

print "Writing output...\n";

open(my $output, ">".$parameters{"pathOutput"});

my $line="ID";

foreach my $sample (@samples){
    $line.="\t".$sample;
}
print $output $line."\n";

my @ids=keys %coverage;
my @sortedids=sort @ids;

foreach my $id (@sortedids){
    my $line=$id;

    foreach my $sample (@samples){
	$line.="\t".$coverage{$id}{$sample};
    }
    print $output $line."\n";
}

close($output);

print "Done.\n";

####################################################################################
####################################################################################
