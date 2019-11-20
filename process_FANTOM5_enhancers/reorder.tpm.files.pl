use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use strict;

################################################################################
################################################################################

sub readSamples{
    my $pathin=$_[0];
    my $samples=$_[1];

    open(my $input, $pathin);
    my $line=<$input>;
    
    while($line){
	chomp $line;
	push(@{$samples}, $line);
	
	$line=<$input>;
    }
    close($input);
}

################################################################################

sub readPeaks{
    my $pathin=$_[0];
    my $peaks=$_[1];

    open(my $input, $pathin);
    my $line=<$input>;
    
    while($line){
	chomp $line;
	$peaks->{$line}=1;
	
	$line=<$input>;
    }
    close($input);
}

################################################################################

sub reorderTPM{
    my $pathin=$_[0];
    my $sampleorder=$_[1];
    my $peaklist=$_[2];
    my $pathout=$_[3];

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

    ## check that all samples are there

    foreach my $sample (@{$sampleorder}){
	if(!exists $header{$sample}){
	    print "cannot find ".$sample."\n";
	    exit(1);
	}
    }
    
    
    open(my $output, ">".$pathout);
    my $lineout="ID\t".join("\t", @{$sampleorder})."\n";
    print $output $lineout;

    $line=<$input>;

    while($line){
	chomp $line;
	my @s=split("\t", $line);
	my $id=$s[0];

	if(exists $peaklist->{$id}){
	    my $lineout=$s[0];
	    
	    foreach my $sample (@{$sampleorder}){
		$lineout.="\t".$s[$header{$sample}];
	    }
	    
	    print $output $lineout."\n";
	}
	
	$line=<$input>;
    }
    
    close($output);
    close($input);

}

################################################################################

sub printHelp{

    my $parnames=$_[0];
    my $parvalues=$_[1];
    
    print "\n";
    print "This script reorders TPM files. \n";
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

$parameters{"pathSelectedSamples"}="NA";
$parameters{"pathSelectedPeaks"}="NA";
$parameters{"pathOriginalTPM"}="NA";
$parameters{"pathReorderedTPM"}="NA";

my %defaultvalues;
my @defaultpars=("pathSelectedSamples", "pathSelectedPeaks", "pathOriginalTPM", "pathReorderedTPM");

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

print "Reading samples...\n";

my @samples;
readSamples($parameters{"pathSelectedSamples"}, \@samples);

my $nbs=@samples;

print "Found ". $nbs." samples.\n";

print "Done.\n";

#####################################################################

print "Reading peaks...\n";

my %peaks;
readPeaks($parameters{"pathSelectedPeaks"}, \%peaks);

my $nbp=keys %peaks;

print "Found ". $nbp." peaks.\n";

print "Done.\n";

#####################################################################

print "Reordering TPM file...\n";

reorderTPM($parameters{"pathOriginalTPM"}, \@samples, \%peaks, $parameters{"pathReorderedTPM"});
 
print "Done.\n";

#####################################################################
