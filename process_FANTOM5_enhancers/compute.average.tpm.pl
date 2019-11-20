use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use strict;

################################################################################
################################################################################

sub readSamplesCategories{
    my $pathin=$_[0];
    my $samplecategories=$_[1];
   
    open(my $input, $pathin);
    my $line=<$input>;
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
	my $sample=$s[$header{"SampleID"}];
	my $type=$s[$header{"TypeID"}];

	if(exists $samplecategories->{$type}){
	    push(@{$samplecategories->{$type}}, $sample);
	} else{
	    $samplecategories->{$type}=[$sample];
	}
	
	$line=<$input>;
    }
    close($input);
}

################################################################################

sub computeAverageTPM{
    my $pathin=$_[0];
    my $samplecategories=$_[1];
    my $pathout=$_[2];

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

    my @types=keys %{$samplecategories};
    my @orderedtypes=sort @types;
  
    open(my $output, ">".$pathout);
    my $lineout="ID\t".join("\t", @orderedtypes)."\n";
    print $output $lineout;

    $line=<$input>;

    while($line){
	chomp $line;
	my @s=split("\t", $line);
	my $id=$s[0];

	my $lineout=$id;

	foreach my $type (@orderedtypes){
	    my @tpmvalues=();
	    
	    foreach my $sample (@{$samplecategories->{$type}}){
		if(!exists $header{$sample}){
		    print "Cannot find ".$sample."!!\n";
		    exit(1);
		}
		
		my $index=$header{$sample};
		my $thistpm=$s[$index]+0;
		push(@tpmvalues, $thistpm);
	    }
	    
	    my $nbval=@tpmvalues;
	    
	    if($nbval>1){
		my $sumtpm=sum (@tpmvalues);
		my $meantpm=($sumtpm+0.0)/($nbval+0.0);
		
		$lineout.="\t".$meantpm;
	    } else{
		$lineout.="\t".$tpmvalues[0];
	    }
	}
	
	print $output $lineout."\n";
	
	
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
    print "This script computes average TPM values per sample category. \n";
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

$parameters{"pathSampleCategories"}="NA";
$parameters{"pathOriginalTPM"}="NA";
$parameters{"pathAverageTPM"}="NA";

my %defaultvalues;
my @defaultpars=("pathSampleCategories",  "pathOriginalTPM", "pathAverageTPM");

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

print "Reading sample categories...\n";

my %samplecategories;

readSamplesCategories($parameters{"pathSampleCategories"}, \%samplecategories);
  
my $nbt=keys %samplecategories;

print "Found ". $nbt." sample categories.\n";

print "Done.\n";

#####################################################################

print "Reading TPM file and computing average values...\n";

computeAverageTPM($parameters{"pathOriginalTPM"}, \%samplecategories, $parameters{"pathAverageTPM"});
 
print "Done.\n";

#####################################################################
