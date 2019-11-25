use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use strict;

####################################################################################
####################################################################################

sub readContacts{
    my $pathin=$_[0];
    my $promoters=$_[1];
    my $enhancers=$_[2];
    my $contacts=$_[3];

    open(my $input, $pathin);
    my $line=<$input>;
    chomp $line;
    my @s=split("\t", $line);
    my %header;

    for(my $i=0; $i<@s; $i++){
	$header{$s[$i]}=$i;
    }
    
    $line=<$input>;

    my $nbread=0;
    
    while($line){
	chomp $line;
	my @s=split("\t", $line);
	my $prom=$s[$header{"IDPromoter"}];
	my $en=$s[$header{"IDEnhancer"}];

	$promoters->{$prom}=1;
	$enhancers->{$en}=1;

	if(exists $contacts->{$prom}){
	    $contacts->{$prom}{$en}=1;
	} else{
	    $contacts->{$prom}={$en=>1};
	}

	$line=<$input>;

	$nbread++;

	if($nbread%100==0){
	    print "Read ".$nbread." contacts.\n";
	}
    }

    close($input);
}

####################################################################################

sub readExpressionData{
    my $pathin=$_[0];
    my $selectedregions=$_[1];
    my $expdata=$_[2]; ## we assume samples are in the same order! 

    open(my $input, $pathin);
    my $line=<$input>; ## header
   
    $line=<$input>;

    my $nbread=0;

    while($line){
	chomp $line;
	my @s=split("\t", $line);

	my $id=$s[0];

	if(exists $selectedregions->{$id}){
	    $expdata->{$id}=[];

	    for(my $i=1; $i<@s; $i++){
		push(@{$expdata->{$id}}, ($s[$i]+0.0));
	    }
	}
	
	$line=<$input>;

	$nbread++;

	if($nbread%10000==0){
	    print "Read ".$nbread." lines.\n";
	}
    }

    close($input);
}

####################################################################################

sub computeSD{
    my $x=$_[0];
    
    my $n=@{$x};
    my $meanx=(sum @{$x})/$n;

    my $scem=0;

    for(my $i=0; $i<$n; $i++){
	$scem+=($x->[$i]-$meanx)^2;
    }

    my $sd=sqrt($scem/$n);

    return($sd);
}

####################################################################################

sub computeCovariance{
    my $x=$_[0];
    my $y=$_[1];

    my $nx=@{$x};
    my $ny=@{$y};

    if($nx!=$ny){
	print "x and y must have the same length.\n";
	exit(1);
    }

    my $meanx=(sum @{$x})/$nx;
    my $meany=(sum @{$y})/$ny;

    my $cov=0;

    for(my $i=0; $i<$nx; $i++){
	$cov+=($x->[$i]-$meanx)*($y->[$i]-$meany);
    }

    $cov=$cov/$nx;

    return($cov);
}

####################################################################################

sub computeCorrelation{
    my $x=$_[0];
    my $y=$_[1];

    my $cov=computeCovariance($x,$y);
    my $sdx=computeSD($x);
    my $sdy=computeSD($y);
    
    my $corr=$cov/($sdx*$sdy);

    return($corr);
}

####################################################################################

sub printHelp{
    my $parnames=$_[0];
    my $parvalues=$_[1];
    
    print "\n";
    print "This script computes expression correlations for promoter-enhancer pairs.\n";
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

$parameters{"pathContacts"}="NA";
$parameters{"pathPromoterExpression"}="NA";
$parameters{"pathEnhancerExpression"}="NA";
$parameters{"pathOutput"}="NA";

my @defaultpars=("pathContacts", "pathPromoterExpression", "pathEnhancerExpression", "pathOutput");

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

print "Reading contacts...\n";

my %promoters;
my %enhancers;
my %contacts;

readContacts($parameters{"pathContacts"}, \%promoters, \%enhancers, \%contacts);

my $nbp=keys %promoters;
my $nbe=keys %enhancers;

print "Found ".$nbp." promoters and ".$nbe." enhancers.\n";
  
print "Done.\n";

####################################################################################

print "Reading expression data for promoters...\n";

my %promoterexp;

readExpressionData($parameters{"pathPromoterExpression"}, \%promoters, \%promoterexp);

my $nbp=keys %promoterexp;

print "Found expression data for ".$nbp." promoters.\n";
   
print "Done.\n";

####################################################################################

print "Reading expression data for enhancers...\n";

my %enhancerexp;

readExpressionData($parameters{"pathEnhancerExpression"}, \%enhancers, \%enhancerexp);

my $nbe=keys %enhancerexp;

print "Found expression data for ".$nbe." enhancers.\n";
   
print "Done.\n";

####################################################################################

print "Computing correlations and writing output...\n";

open(my $output, ">".$parameters{"pathOutput"});

print $output "IDPromoter\tIDEnhancer\tPearsonCorrelation\n";

my $nbdone=0;

foreach my $prom (keys %contacts){
    
    my @pexp=@{$promoterexp{$prom}};
    
    foreach my $en (keys %{$contacts{$prom}}){
	my @eexp=@{$enhancerexp{$en}};
	
	my $c=computeCorrelation(\@pexp, \@eexp);

	print $output $prom."\t".$en."\t".$c."\n";
    }
    
    $nbdone++;

    if($nbdone%100==0){
	print $nbdone." enhancers done.\n";
    }
}

close($output);

print "Done.\n";

####################################################################################
