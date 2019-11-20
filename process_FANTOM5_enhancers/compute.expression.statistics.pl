use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use strict;

################################################################################
################################################################################

sub computeExpressionStatistics{
    my $pathin=$_[0];
    my $minvalue=$_[1];
    my $pathoutnbexppeaks=$_[2];
    my $pathoutnbexpsamples=$_[3];

    my @s=split("\\.", $pathin);
    my $ext=$s[-1];
    my $input;

    if($ext eq "gz"){
	open($input, "zcat $pathin |")
    } else{
	open($input, $pathin);
    }

    my $line=<$input>;
    my $firstchar=substr $line, 0, 1;

    while($firstchar eq "#"){
	$line=<$input>;
	$firstchar=substr $line, 0, 1;
    }

    ## sample header
    chomp $line;
    my @s=split("\t", $line);
    my %header;
    $header{"Annotation"}=0;
    
    for(my $i=1; $i<@s; $i++){
	my @t=split("\\.", $s[$i]);
	my @ss = grep(/CNhs/, @t);
	my $nbss=@ss;
	
	if($nbss!=1){
	    print "found ".$nbss." samples.\n";
	    exit(1);
	}
	
	my $sample=$ss[0];
	$header{$sample}=$i;
    }

    ## compute number of expressed genes by sample

    my %nbexpressed;
    foreach my $sample (keys %header){
	if($sample ne "Annotation"){
	    $nbexpressed{$sample}=0;
	}
    }
    
    $line=<$input>;

    open(my $outputnbsamples, ">".$pathoutnbexpsamples);
    print $outputnbsamples "PeakID\tNbSamples\n";

    while($line){
	chomp $line;
	my @s=split("\t", $line);

	my $annot=$s[$header{"Annotation"}];
	my $firstchar=substr $annot, 0, 1;

	if($firstchar ne "0"){
	    my $nbexp=0;
	    
	    foreach my $sample (keys %header){
		if($sample ne "Annotation"){
		    my $val=$s[$header{$sample}]+0;
		    
		    if($val>=$minvalue){
			$nbexpressed{$sample}++;
			$nbexp++;
		    }
		}
	    }
	    print $outputnbsamples $annot."\t".$nbexp."\n"; 
	}
	
	$line=<$input>;
    }
   
    close($input);
    close($outputnbsamples);

    # write output for number of expressed peaks

    open(my $outputnbpeaks, ">". $pathoutnbexppeaks);
    print $outputnbpeaks "SampleID\tNbPeaks\n";
    
    foreach my $sample (keys %nbexpressed){
	print $outputnbpeaks $sample."\t".$nbexpressed{$sample}."\n";
    }

    close($outputnbpeaks);
}

################################################################################

sub printHelp{

    my $parnames=$_[0];
    my $parvalues=$_[1];
    
    print "\n";
    print "This script computes numbers of peaks expressed above a certain threshold in each sample. \n";
    print "\n";
    print "Options:\n";
    
    foreach my $par (@{$parnames}){
	print "--".$par."  [  default value: ".$parvalues->{$par}."  ]\n";
    }
    print "\n";
}

################################################################################
################################################################################

my %parameters;

$parameters{"pathTPM"}="NA";
$parameters{"minTPM"}="NA";
$parameters{"pathOutputStatisticsSamples"}="NA";
$parameters{"pathOutputStatisticsPeaks"}="NA";

my %defaultvalues;
my @defaultpars=("pathTPM", "minTPM", "pathOutputStatisticsSamples", "pathOutputStatisticsPeaks");

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

#####################################################################
#####################################################################

print "Computing expression statistics and writing output...\n";

my $minval=$parameters{"minTPM"}+0;

print "minimum value: ".$minval."\n";

computeExpressionStatistics($parameters{"pathTPM"}, $minval, $parameters{"pathOutputStatisticsSamples"}, $parameters{"pathOutputStatisticsPeaks"});

print "Done.\n";

#####################################################################
