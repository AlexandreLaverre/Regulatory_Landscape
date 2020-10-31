use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use strict;

################################################################################

sub readCoordinates{
    my $pathin=$_[0];
    my $coords=$_[1];

    my @s=split("\\.", $pathin);
    my $ext=$s[-1];

    my $input;

    if($ext eq "gz"){
	open($input, "zcat $pathin");
    } else{
	open($input, $pathin);
    }
   
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

	my $chr=$s[$header{"chr"}];
	my $start=$s[$header{"start"}];
	my $end=$s[$header{"end"}];

	my $id=$chr.":".$start.":".$end;

	$coords->{$id}=1;

	$line=<$input>;
    }
    
    close $input;
}

################################################################################

sub readAlignmentStatistics{
    my $pathin=$_[0];
    my $ref=$_[1];
    my $tg=$_[2];
    my $aln=$_[3];

    my @s=split("\\.", $pathin);
    my $ext=$s[-1];

    my $input;

    if($ext eq "gz"){
	open($input, "zcat $pathin");
    } else{
	open($input, $pathin);
    }
   
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

	my $idref=$s[$header{"ID.".$ref}];
	my $idtg=$s[$header{"ID.".$tg}];

	my @t=split(":", $idref);
	my $id=join(":", @t[0..2]);

	my $identicallength=$s[$header{"FilteredIdenticalLength"}]+0;
	my $ungappedlength=$s[$header{"FilteredUngappedLength"}]+0;
	my $alnlength=$s[$header{"FilteredAlignmentLength"}]+0;

	my $frid="NA";
	my $frungapped="NA";

	if($alnlength!=0){
	    $frid=$identicallength/$alnlength;
	}
	if($alnlength!=0){
	    $frungapped=$ungappedlength/$alnlength;
	}

	$aln->{$id}={"ungapped"=>$frungapped, "identity"=>$frid};
	
	$line=<$input>;
    }
    
    close $input;
}

################################################################################

sub printHelp{
    my $parnames=$_[0];
    my $parvalues=$_[1];
    
    print "\n";
    print "This script combine alignment statistics for multiple species. \n";
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
$parameters{"pathCoordinates"}="NA";
$parameters{"refSpecies"}="NA";
$parameters{"targetSpecies"}="NA";
$parameters{"pathsAlignmentStatistics"}="NA";
$parameters{"pathOutputUngappedAlignment"}="NA";
$parameters{"pathOutputSequenceIdentity"}="NA";

my %defaultvalues;
my @defaultpars=("pathCoordinates", "refSpecies", "targetSpecies", "pathsAlignmentStatistics", "pathOutputUngappedAlignment", "pathOutputSequenceIdentity");

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

print "Reading element coordinates...\n";

my %coords;

readCoordinates($parameters{"pathCoordinates"}, \%coords);

print "Done.\n";

#####################################################################

print "Reading alignment statistics...\n";

my $ref=$parameters{"refSpecies"};

my @tg=split(",", $parameters{"targetSpecies"});
my @paths=split(",", $parameters{"pathsAlignmentStatistics"});

my $nbt=@tg;
my $nbp=@paths;

if($nbt!=$nbp){
    print "Weird, there are ".$nbt." target species and ".$nbp." paths.\n";
    exit(1);
}

my %alnstats;

for(my $i=0; $i<$nbt; $i++){
    my $sp=$tg[$i];
    
    $alnstats{$sp}={};
    
    readAlignmentStatistics($paths[$i], $ref, $sp, $alnstats{$sp});
}

print "Done.\n";

#####################################################################

print "Writing output for percentage identity...\n";

open(my $output, ">".$parameters{"pathOutputSequenceIdentity"});

my $line="ID\t".join("\t", @tg)."\n";

print $output $line;

foreach my $id (keys %coords){
    my $line=$id;

    foreach my $sp (@tg){
	if(exists $alnstats{$sp}{$id}){
	    $line.="\t".$alnstats{$sp}{$id}{"identity"};
	} else{
	    $line.="\tNA";
	}
    }

    print $output $line."\n";
}

close($output);

print "Done.\n";

#####################################################################

print "Writing output for ungapped alignment..\n";

open(my $output, ">".$parameters{"pathOutputUngappedAlignment"});

my $line="ID\t".join("\t", @tg)."\n";

print $output $line;

foreach my $id (keys %coords){
    my $line=$id;

    foreach my $sp (@tg){
	if(exists $alnstats{$sp}{$id}){
	    $line.="\t".$alnstats{$sp}{$id}{"ungapped"};
	} else{
	    $line.="\tNA";
	}
    }

    print $output $line."\n";
}

close($output);

print "Done.\n";


#####################################################################
#####################################################################
