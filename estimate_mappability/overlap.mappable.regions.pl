use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use strict;

####################################################################################
####################################################################################

sub readCoordinatesBED{
    my $pathin=$_[0];
    my $coords=$_[1];

    open(my $input, $pathin);
    my $line=<$input>;

    my %unordered;
    
    while($line){
	chomp $line;
	my @s=split("\t", $line);
	my $nbf=@s;
	
	my $chr=$s[0];
	my $prefix=substr $chr, 0, 3;
	
	if($prefix eq "chr"){
	    $chr=substr $chr, 3;
	}
	
	my $start=$s[1]+0;
	my $end=$s[2]+0;
	
	my $id="NA";

	if($nbf>=4){
	    $id=$s[3];
	} else{
	    $id=$chr.":".$start.":".$end;
	}

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

		    ## we count the total number of bases
		    ## intervals cannot overlap
		    if($M<=$m){
			my $ovlen=$m-$M+1;
			
			if(exists $overlap->{$id1}){
			    $overlap->{$id1}{"length"}+=$ovlen;
			    my $ml=$overlap->{$id1}{"maxstretch"};

			    if($ovlen>$ml){
				$overlap->{$id1}{"maxstretch"}=$ovlen;
			    }
			} else{
			    $overlap->{$id1}={"length"=>$ovlen, "maxstretch"=>$ovlen};
			}
		    }
		    
		    $j++;
		}
	    }
	}	
    }
}

####################################################################################

sub printHelp{
    my $parnames=$_[0];
    my $parvalues=$_[1];
    
    print "\n";
    print "This script extracts mappability statistics for restriction fragments.\n";
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

$parameters{"pathRestrictionFragments"}="NA";
$parameters{"pathMappedRegions"}="NA";
$parameters{"margin"}="NA";
$parameters{"step"}="NA";
$parameters{"pathOutput"}="NA";

my @defaultpars=("pathRestrictionFragments", "pathMappedRegions", "margin","step", "pathOutput");

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

print "Reading restriction fragments...\n";

my %fragments;
readCoordinatesBED($parameters{"pathRestrictionFragments"}, \%fragments);

print "Done.\n";

####################################################################################

print "Reading mapped regions...\n";

my %mappedregions;
readCoordinatesBED($parameters{"pathMappedRegions"}, \%mappedregions);

print "Done.\n";

####################################################################################

print "Computing total overlap...\n";

my %totoverlap;
overlapCoordinates(\%fragments, \%mappedregions, \%totoverlap);

print "Done.\n";

####################################################################################

print "Constructing additional intervals and computing overlaps...\n";

my $margin=$parameters{"margin"}+0;
my $step=$parameters{"step"}+0;

if($step==0){
    print "stopping, step cannot be 0.\n";
    exit(1);
}

print "margin ".$margin." step ".$step."\n";

my %overlapleft;
my %overlapright;

my $maxoffset=0;

for(my $offset=0; $offset<$margin; $offset+=$step){

    print "offset ".$offset."\n";
 
    my %coordsleft;
    my %coordsright;

    foreach my $chr (keys %mappedregions){
	$coordsleft{$chr}={"start"=>[], "end"=>[], "id"=>[]};
	$coordsright{$chr}={"start"=>[], "end"=>[], "id"=>[]};
	
	my $nb=@{$fragments{$chr}{"start"}};

	for(my $i=0; $i<$nb; $i++){
	    my $start=${$fragments{$chr}{"start"}}[$i];
	    my $end=${$fragments{$chr}{"end"}}[$i];
	    my $id=${$fragments{$chr}{"id"}}[$i];
	    my $midpos=($start+$end)/2;
	    
	    my $newstart=$start+$offset;
	    my $newend=$newstart+$step-1;

	    if($newend<$midpos){
		push(@{$coordsleft{$chr}{"start"}}, $newstart);
		push(@{$coordsleft{$chr}{"end"}}, $newend);
		push(@{$coordsleft{$chr}{"id"}}, $id);
	    }

	    $newend=$end-$offset;
	    $newstart=$newend-$step+1;

	    if($newstart>$midpos){
		push(@{$coordsright{$chr}{"start"}}, $newstart);
		push(@{$coordsright{$chr}{"end"}}, $newend);
		push(@{$coordsright{$chr}{"id"}}, $id);
	    }
	}
    }

    $overlapleft{$offset}={};
    $overlapright{$offset}={};

    overlapCoordinates(\%coordsleft, \%mappedregions, $overlapleft{$offset});
    overlapCoordinates(\%coordsright, \%mappedregions, $overlapright{$offset});

    $maxoffset=$offset;
}

print "Done.\n";

####################################################################################

print "Writing output...\n";

open(my $output, ">".$parameters{"pathOutput"});

my $line="IDFragment\tChr\tStart\tEnd\tTotalLength\tTotalMappableLength\tMaxMappableStretch";

for(my $offset=0; $offset<$margin; $offset+=$step){
    $line.="\tL".$offset;
}

for(my $offset=$maxoffset; $offset>=0; $offset-=$step){
    $line.="\tR".$offset;
}

print $output $line."\n";

 foreach my $chr (keys %mappedregions){
     
     my $nb=@{$fragments{$chr}{"start"}};
     
     for(my $i=0; $i<$nb; $i++){
	 my $start=${$fragments{$chr}{"start"}}[$i];
	 my $end=${$fragments{$chr}{"end"}}[$i];
	 my $id=${$fragments{$chr}{"id"}}[$i];
	 my $midpos=($start+$end)/2;
	 my $len=$end-$start+1;

	 my $lenov=0;
	 my $maxstretch=0;
	 
	 if(exists $totoverlap{$id}){
	     $lenov=$totoverlap{$id}{"length"};
	     $maxstretch=$totoverlap{$id}{"maxstretch"};
	 }

	 my $line=$id."\t".$chr."\t".$start."\t".$end."\t".$len."\t".$lenov."\t".$maxstretch;

	 for(my $offset=0; $offset<$margin; $offset+=$step){
 	    my $newstart=$start+$offset;
	    my $newend=$newstart+$step-1;

	    my $lenovleft="NA";
	    
	    if($newend<$midpos){
		$lenovleft=0;

		if(exists $overlapleft{$offset}{$id}){
		    $lenovleft=$overlapleft{$offset}{$id}{"length"};
		}
	    }
	    
	    $line.="\t".$lenovleft;
	 }

	 for(my $offset=$maxoffset; $offset>=0; $offset-=$step){
	     my $newend=$end-$offset;
	     my $newstart=$newend-$step+1;

	     my $lenovright="NA";
	     
	     if($newstart>$midpos){
		 $lenovright=0;
		 
		 if(exists $overlapright{$offset}{$id}){
		    $lenovright=$overlapright{$offset}{$id}{"length"};
		}
	    }
	    
	    $line.="\t".$lenovright;
	 }

	 print $output $line."\n";
     }
}


close($output);

print "Done.\n";

####################################################################################

