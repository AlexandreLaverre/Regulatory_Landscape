use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use strict;

################################################################################

sub printHelp{
    my $parnames=$_[0];
    my $parvalues=$_[1];
    
    print "\n";
    print "This script removes strand info from sequence ids. \n";
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
$parameters{"pathInput"}="NA";
$parameters{"pathOutput"}="NA";

my %defaultvalues;
my @defaultpars=("pathInput", "pathOutput");

my @numericpars=();


my %numericpars;

foreach my $par (@numericpars){
    $numericpars{$par}=1;
}

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

print "Running program with the following parameters:\n";

foreach my $par (@defaultpars){
    print "--".$par."=".$parameters{$par}."\n";
}

print "\n";

#####################################################################
#####################################################################

if(-e $parameters{"pathInput"}) {

    print "Reading input and writing cleaned up output...\n";

    my $pathin=$parameters{"pathInput"};

    my @s=split("\t", $pathin);
    my $ext=$s[-1];

    my $input;

    if($ext eq "gz"){
	open($input, "zcat $pathin |");
    } else{
	open($input, $pathin);
    }

    open(my $output, ">".$parameters{"pathOutput"});

    my $line=<$input>; ## header
    print $output $line;

    $line=<$input>;

    while($line){
	chomp $line;
	my @s=split("\t", $line);

	my $id1=shift @s;
	my @t1=split(":", $id1);
	$id1=join(":", @t1[0..2]);

	my $id2=shift @s;
	my @t2=split(":", $id2);
	$id2=join(":", @t2[0..2]);

	print $output $id1."\t".$id2."\t".join("\t", @s)."\n";
	
	$line=<$input>;
    }

    close($input);
    close($output);

    print "Done.\n";
} else{
    print "Cannot find input file, not doing anything.\n";
}

#####################################################################
#####################################################################
