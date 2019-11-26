use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use strict;

####################################################################################
####################################################################################

sub printHelp{
    my $parnames=$_[0];
    my $parvalues=$_[1];
    
    print "\n";
    print "This script selects samples based on ChIPSeq signal strength.\n";
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

$parameters{"pathCombinedCoverage"}="NA";
$parameters{"minCoverage"}="NA";
$parameters{"minSamples"}="NA";
$parameters{"pathOutput"}="NA";

my @defaultpars=("pathCombinedCoverage", "minCoverage", "minSamples", "pathOutput");

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

print "Reading input and filtering regions...\n";

my $mincoverage=$parameters{"minCoverage"}+0;
my $minsamples=$parameters{"minSamples"}+0;

print "Selecting regions that have at least ".$mincoverage." coverage for at least ".$minsamples." samples.\n";

open(my $input, $parameters{"pathCombinedCoverage"});
open(my $output, ">".$parameters{"pathOutput"});

my $line=<$input>; ## header

$line=<$input>;

my $nbtotal=0;
my $nbkept=0;

while($line){
    chomp $line;
    my @s=split("\t", $line);

    my $nbok=0;
    
    for(my $i=1; $i<@s; $i++){ ## first column is the ID
	if($s[$i]>=$mincoverage){
	    $nbok++;
	}

	if($nbok>=$minsamples){
	    last;
	}
    }

    if($nbok>=$minsamples){
	my $id=$s[0];
	my @t=split(",", $id);
	print $output join("\t", @t)."\t".$id."\n"; ##bed format
	
	$nbkept++;
    }

    $nbtotal++;
    
    $line=<$input>;
}

close($input);
close($output);


print "Kept ".$nbkept." regions out ouf ".$nbtotal."\n";

print "Done.\n";


####################################################################################
####################################################################################
