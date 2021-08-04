use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use strict;

####################################################################################
####################################################################################

sub readCoordinatesBED{
    my $pathin=$_[0];
    my $coords=$_[1];
    my $hashcoords=$_[2];

    open(my $input, $pathin);
    my $line=<$input>;

    my %unordered;
    
    while($line){
	chomp $line;
	my @s=split("\t", $line);

	my $chr=$s[0];
	my $start=$s[1]+0;
	my $end=$s[2]+0;
	my $id=$s[3];

	$hashcoords->{$id}={"chr"=>$chr, "start"=>$start, "end"=>$end};

	if(exists $unordered{$chr}){
	    push(@{$unordered{$chr}{"start"}}, $start);
	    push(@{$unordered{$chr}{"end"}}, $end);
	    push(@{$unordered{$chr}{"id"}}, $id);
	} else{
	    $unordered{$chr}={"start"=>[$start], "end"=>[$end], "id"=>[$id]};
	}
	
	$line=<$input>;
    }
    
    close($input);

    orderCoordinates(\%unordered, $coords);
}

####################################################################################

sub orderCoordinates{
    my $unordered=$_[0];
    my $ordered=$_[1];

    foreach my $chr (keys %{$unordered}){
	$ordered->{$chr}={"start"=>[], "end"=>[], "id"=>[]};
	
	my %hashcoords;
	my $nbpos=@{$unordered->{$chr}{"start"}};

	for(my $i=0; $i<$nbpos; $i++){
	    my $start=${$unordered->{$chr}{"start"}}[$i];
	    my $end=${$unordered->{$chr}{"end"}}[$i];
	    my $id=${$unordered->{$chr}{"id"}}[$i];

	    if(exists $hashcoords{$start}){
		push(@{$hashcoords{$start}{"end"}}, $end);
		push(@{$hashcoords{$start}{"id"}}, $id);
	    } else{
		$hashcoords{$start}={"end"=>[$end], "id"=>[$id]};
	    }
	}

	my @startpos=keys %hashcoords;
	my @orderedstart=sort {$a<=>$b} @startpos;

	foreach my $start (@orderedstart){
	    my $nbpos=@{$hashcoords{$start}{"end"}};

	    for(my $i=0; $i<$nbpos; $i++){
		my $end=${$hashcoords{$start}{"end"}}[$i];
		my $id=${$hashcoords{$start}{"id"}}[$i];

		push(@{$ordered->{$chr}{"start"}}, $start);
		push(@{$ordered->{$chr}{"end"}}, $end);
		push(@{$ordered->{$chr}{"id"}}, $id);
	    }
	}
    }
}

####################################################################################

sub extractRegulatoryRegions{
    my $promoters=$_[0];
    my $maxdist=$_[1];
    my $hashregions=$_[2];
    my $regions=$_[3];

    my %unordered;

    foreach my $chr (keys %{$promoters}){
	my $nbp=@{$promoters->{$chr}{"start"}};

	for(my $i=0; $i<$nbp; $i++){
	    my $id=${$promoters->{$chr}{"id"}}[$i];
	    my $thisstart=${$promoters->{$chr}{"start"}}[$i];
	    my $thisend=${$promoters->{$chr}{"end"}}[$i];

	    my $startregion="NA";
	    my $endregion="NA";
	    
	    if($i<($nbp-1)){
		my $nextstart=${$promoters->{$chr}{"start"}}[$i+1];
		$endregion=min($nextstart-1, $thisend+$maxdist);
	    } else{
		$endregion=$thisend+$maxdist;
	    }

	    if($i>0){
		my $prevend=${$promoters->{$chr}{"end"}}[$i-1];
		$startregion=max($prevend+1, $thisstart-$maxdist);
	    } else{
		$startregion=max(1, $thisstart-$maxdist);
	    }

	    if(exists $unordered{$chr}){
		push(@{$unordered{$chr}{"start"}}, $startregion);
		push(@{$unordered{$chr}{"end"}}, $endregion);
		push(@{$unordered{$chr}{"id"}}, $id);
	    } else{
		$unordered{$chr}={"start"=>[$startregion], "end"=>[$endregion], "id"=>[$id]};
	    }

	    $hashregions->{$id}={"chr"=>$chr, "start"=>$startregion, "end"=>$endregion};
	}
    }

    orderCoordinates(\%unordered, $regions);
}

####################################################################################

sub overlapCoordinates{
    my $coords1=$_[0];
    my $coords2=$_[1];
    my $overlap=$_[2];

    foreach my $chr (keys %{$coords1}){
	if(exists $coords2->{$chr}){
	    my $nb1=@{$coords1->{$chr}{"start"}};
	    my $nb2=@{$coords2->{$chr}{"start"}};

	    my $firstj=0;
	    
	    for(my $i=0; $i<$nb1; $i++){
		my $start1=${$coords1->{$chr}{"start"}}[$i];
		my $end1=${$coords1->{$chr}{"end"}}[$i];
		my $id1=${$coords1->{$chr}{"id"}}[$i];

		my $j=$firstj;

		while($j<$nb2 && ${$coords2->{$chr}{"end"}}[$j]<$start1){
		    $j++;
		}

		$firstj=$j;

		while($j<$nb2 && ${$coords2->{$chr}{"start"}}[$j]<=$end1){
		    my $start2=${$coords2->{$chr}{"start"}}[$j];
		    my $end2=${$coords2->{$chr}{"end"}}[$j];
		    my $id2=${$coords2->{$chr}{"id"}}[$j];

		    my $M=max($start1, $start2);
		    my $m=min($end1, $end2);

		    if($M<=$m){
			if(exists $overlap->{$id1}){
			    $overlap->{$id1}{$id2}=1;
			} else{
			    $overlap->{$id1}={$id2=>1};
			}
		    }
		    
		    $j++;
		}
	    }
	}	
    }
}

####################################################################################
####################################################################################

sub printHelp{
    my $parnames=$_[0];
    my $parvalues=$_[1];
    
    print "\n";
    print "This script extracts neighbor enhancers.\n";
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

$parameters{"pathPromoters"}="NA";
$parameters{"pathEnhancers"}="NA";
$parameters{"maxDistance"}="NA";
$parameters{"pathOutput"}="NA";

my @defaultpars=("pathPromoters", "pathEnhancers", "maxDistance", "pathOutput");

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

print "Reading promoter regions...\n";

my %promoters;
my %hashpromoters;

readCoordinatesBED($parameters{"pathPromoters"}, \%promoters, \%hashpromoters);

print "Done.\n";
    
####################################################################################

print "Reading enhancer regions...\n";

my %enhancers;
my %hashenhancers;

readCoordinatesBED($parameters{"pathEnhancers"}, \%enhancers, \%hashenhancers);

print "Done.\n";
    
####################################################################################

print "Extracting regulatory regions...\n";

my %regions;
my %hashregions;
my $maxdist=$parameters{"maxDistance"}+0;

print "maximum distance: ".$maxdist."\n";

extractRegulatoryRegions(\%promoters, $maxdist, \%hashregions, \%regions);

print "Done.\n";

####################################################################################

print "Extracting overlap between regulatory regions and enhancers...\n";

my %overlap;

overlapCoordinates(\%regions, \%enhancers, \%overlap);

my $nbob=keys %overlap;

print "Found ".$nbob." regulatory regions that overlap with enhancers.\n";

print "Done.\n";

####################################################################################

if($nbob>0){

    print "Writing output...\n";
    
    open(my $output, ">".$parameters{"pathOutput"});
    
    print $output "IDPromoter\tChrPromoter\tStartPromoter\tEndPromoter\tStartRegion\tEndRegion\tIDEnhancer\tStartEnhancer\tEndEnhancer\tDistance\n";
    
    foreach my $idprom (keys %overlap){
	my $chr=$hashpromoters{$idprom}{"chr"};
	
	my $startprom=$hashpromoters{$idprom}{"start"};
	my $endprom=$hashpromoters{$idprom}{"end"};
	my $midposprom=($startprom+$endprom)/2;
	
	my $startregion=$hashregions{$idprom}{"start"};
	my $endregion=$hashregions{$idprom}{"end"};
	
	if(exists $overlap{$idprom}){
	    foreach my $idenh (keys %{$overlap{$idprom}}){
		my $startenh=$hashenhancers{$idenh}{"start"};
		my $endenh=$hashenhancers{$idenh}{"end"};
		my $midposenh=($startenh+$endenh)/2;
		
		my $dist=abs($midposprom-$midposenh);
		
		print $output $idprom."\t".$chr."\t".$startprom."\t".$endprom."\t".$startregion."\t".$endregion."\t".$idenh."\t".$startenh."\t".$endenh."\t".$dist."\n";
	    }
	}
    }
    
    print "Done.\n";
    
    close($output);
    
}

####################################################################################
####################################################################################
