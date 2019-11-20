use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use strict;

##############################################################
##############################################################

sub readCoordinates{
    my $pathin=$_[0];
    my $coords=$_[1];
    
    open(my $input, $pathin);
    
    my $line=<$input>;

    while($line){
	chomp $line;
	my @s=split("\t", $line);
	my $chr=$s[0];
	my $start=$s[1]+0;
	my $end=$s[2]+0;

	## coordinates are already ordered
	
	if(exists $coords->{$chr}){
	    my $lastend=${$coords->{$chr}{"end"}}[-1];
	    
	    if($start<$lastend){
		print "coordinates are not sorted!\n";
		exit(1);
	    }

	    push(@{$coords->{$chr}{"start"}}, $start);
	    push(@{$coords->{$chr}{"end"}}, $end);
	    
	} else{
	    $coords->{$chr}={"start"=>[$start], "end"=>[$end]};
	}
	
	       	
	$line=<$input>;
    }
    
    close($input);
}

##############################################################

sub readCoverage{
    my $pathin=$_[0];
    my $coverage=$_[1];

    my $input;
    my @s=split("\\.",$pathin);
    my $ext=$s[-1];
    
    if($ext eq "gz"){
	open($input, "zcat $pathin |");
    }
    else{
	open($input, $pathin);
    }
    
    my $line=<$input>;

    while($line){
	chomp $line;
	my @s=split("\t" , $line);
	
	my $chr=$s[0];
	my $start=$s[1]+1; ## 1-based, included
	my $end=$s[2]+1;
	my $score=$s[3]+0.0;
	
	if(exists $coverage->{$chr}){
	    my $laststart=${$coverage->{$chr}{"start"}}[-1];
	    
	    if($start<$laststart){
		print "Data are not ordered! ".$start." ".$laststart."\n";
		exit(1);
	    }
	    
	    push(@{$coverage->{$chr}{"start"}}, $start);
	    push(@{$coverage->{$chr}{"end"}}, $end);
	    push(@{$coverage->{$chr}{"score"}}, $score);
	}
	else{
	    $coverage->{$chr}={"start"=>[$start], "end"=>[$end], "score"=>[$score]};
	}
	
	
	$line=<$input>;
    }
    
    close($input);
}

##############################################################

sub computeCoverage{
    my $coords=$_[0];
    my $coverage=$_[1];
    my $covexons=$_[2];

    foreach my $chr (keys %{$coords}){
	if(exists $coverage->{$chr}){
	    my $nbex=@{$coords->{$chr}{"start"}};
	    my $nbreg=@{$coverage->{$chr}{"start"}};
	    
	    my $firstreg=0;
	    
	    for(my $i=0; $i<$nbex; $i++){
		my $startex=${$coords->{$chr}{"start"}}[$i];
		my $endex=${$coords->{$chr}{"end"}}[$i];
		my $idex=$chr.",".$startex.",".$endex;
		
		my $lenex=($endex-$startex+1.0);
		my $sumcov=0;
		
		my $j=$firstreg;
		
		while($j<$nbreg && ${$coverage->{$chr}{"end"}}[$j] < $startex){
		    $j++;
		}
		    
		$firstreg=$j;
		
		while($j<$nbreg && ${$coverage->{$chr}{"start"}}[$j] <= $endex){
		    my $startreg=${$coverage->{$chr}{"start"}}[$j];
		    my $endreg=${$coverage->{$chr}{"end"}}[$j];
			
		    my $M=max($startex, $startreg);
		    my $m=min($endex, $endreg);
		    
		    my $lenov=($m-$M+1);
		    
		    if($lenov>=1){
			my $score=${$coverage->{$chr}{"score"}}[$j];
			
			$sumcov+=($lenov+0.0)*($score+0.0);
		    }
		    
		    $j++;
		}
		
		my $meancov=($sumcov+0.0)/($lenex);
		
		$covexons->{$idex}=$meancov;
	    }
	}
    }
}

##############################################################

sub printHelp{
    my $parnames=$_[0];
    my $parvalues=$_[1];
    
    print "\n";
    print "This script computes read coverage for input regions.\n";
    print "\n";
    
    print "Options:\n";
    
    foreach my $par (@{$parnames}){
	print "--".$par."  [  default value: ".$parvalues->{$par}."  ]\n";
    }

    print "\n";
}

##############################################################
##############################################################

my %parameters;
$parameters{"pathCoordinates"}="NA";
$parameters{"pathCoverage"}="NA";
$parameters{"pathOutput"}="NA";

my @defaultpars=("pathCoordinates", "pathCoverage", "pathOutput");

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

##############################################################
##############################################################

print "Reading coordinates...\n";

my %coords;
readCoordinates($parameters{"pathCoordinates"}, \%coords);

my $nbchr=keys %coords;

print "Found regions on ".$nbchr." chromosomes .\n";

print "Done.\n";

##############################################################

print "Reading coverage...\n";

my %coverage;

readCoverage($parameters{"pathCoverage"}, \%coverage);

print "Done.\n";

##############################################################

print "Computing coverage for input regions...\n";
 
my %coverageregions;

computeCoverage(\%coords, \%coverage, \%coverageregions);

print "Done.\n";

##############################################################

print "Writing output...\n";

open(my $output, ">".$parameters{"pathOutput"});

print $output "Chr\tStart\tEnd\tCoverage\n";

foreach my $chr (keys %coords){
    my $nb=@{$coords{$chr}{"start"}};

    for(my $i=0; $i<$nb; $i++){
	my $start=${$coords{$chr}{"start"}}[$i];
	my $end=${$coords{$chr}{"end"}}[$i];
	my $id=$chr.",".$start.",".$end;
	
	if(exists $coverageregions{$id}){
	    print $output $chr."\t".$start."\t".$end."\t".$coverageregions{$id}."\n";
	} else{
	    print $output $chr."\t".$start."\t".$end."\t0\n";
	}
    }
}

close($output);

print "Done.\n";

##############################################################
