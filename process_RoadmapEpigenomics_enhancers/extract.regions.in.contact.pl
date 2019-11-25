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

	my $chr=$s[0];
	my $start=$s[1]+0;
	my $end=$s[2]+0;
	my $id=$s[3];

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

sub readContacts{
    my $pathin=$_[0];
    my $baitcoords=$_[1];
    my $fragmentcoords=$_[2];
    my $contacts=$_[3];

    open(my $input, $pathin);

    my %header;
    my $line=<$input>;
    chomp $line;
    my @s=split("\t", $line);

    for(my $i=0; $i<@s; $i++){
	$header{$s[$i]}=$i;
    }

    $line=<$input>;

    my %unorderedbaits;
    my %unorderedfragments;
    
    while($line){
	chomp $line;
	my @s=split("\t", $line);
	
	my $baitchr=$s[$header{"bait_chr"}];
	my $otherchr=$s[$header{"chr"}];

	if($baitchr eq $otherchr){ ## only cis
	    my $baitedfrag=$s[$header{"baited_frag"}];
	
	    if($baitedfrag eq "unbaited"){ ## only contacts between baits and unbaited fragments
		my $baitstart=$s[$header{"bait_start"}]+0;
		my $baitend=$s[$header{"bait_end"}]+0;
		my $otherstart=$s[$header{"start"}]+0;
		my $otherend=$s[$header{"end"}]+0;

		my $idbait=$baitchr.",".$baitstart.",".$baitend;
		my $idother=$otherchr.",".$otherstart.",".$otherend;
		
		if(exists $unorderedbaits{$baitchr}){
		    push(@{$unorderedbaits{$baitchr}{"start"}}, $baitstart);
		    push(@{$unorderedbaits{$baitchr}{"end"}}, $baitend);
		    push(@{$unorderedbaits{$baitchr}{"id"}}, $idbait);
		} else{
		    $unorderedbaits{$baitchr}={"start"=>[$baitstart], "end"=>[$baitend], "id"=>[$idbait]};
		}
		
		if(exists $unorderedfragments{$otherchr}){
		    push(@{$unorderedfragments{$otherchr}{"start"}}, $otherstart);
		    push(@{$unorderedfragments{$otherchr}{"end"}}, $otherend);
		    push(@{$unorderedfragments{$otherchr}{"id"}}, $idother);
		} else{
		    $unorderedfragments{$otherchr}={"start"=>[$otherstart], "end"=>[$otherend], "id"=>[$idother]};
		}

		if(exists $contacts->{$idbait}){
		    push(@{$contacts->{$idbait}}, $idother); 
		} else{
		    $contacts->{$idbait}=[$idother];
		}
	    }
	}
    }
    
    close($input);

    ## order coordinates

    orderCoordinates(\%unorderedbaits, $baitcoords);
    orderCoordinates(\%unorderedfragments, $fragmentcoords);

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
			    push(@{$overlap->{$id1}}, $id2);
			} else{
			    $overlap->{$id1}=[$id2];
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
    print "This script extracts promoters and enhancers in contact.\n";
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
$parameters{"pathPromoters"}="NA";
$parameters{"pathEnhancers"}="NA";
$parameters{"pathOutput"}="NA";

my @defaultpars=("pathContacts", "pathPromoters", "pathEnhancers", "pathOutput");

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

print "Reading regions in contact...\n";

my %baitcoords;
my %contactedcoords;
my %contacts;

readContacts($parameters{"pathContacts"}, \%baitcoords, \%contactedcoords, \%contacts);

my $nbbaits=keys %contacts;

print "Found contacts for ".$nbbaits."\n";

print "Done.\n";

####################################################################################

print "Reading promoter regions...\n";

my %promoters;
readCoordinatesBED($parameters{"pathPromoters"}, \%promoters);

print "Done.\n";
    
####################################################################################

print "Reading enhancer regions...\n";

my %enhancers;
readCoordinatesBED($parameters{"pathEnhancers"}, \%enhancers);

print "Done.\n";
    
####################################################################################

print "Extracting overlap between baits and promoter regions...\n";

my %overlapbaitpromoters;

overlapCoordinates(\%baitcoords, \%promoters, \%overlapbaitpromoters);

my $nbob=keys %overlapbaitpromoters;

print "Found ".$nbob." baits that overlap with promoters.\n";

print "Done.\n";

####################################################################################

print "Extracting overlap between contacted fragments and enhancer regions...\n";

my %overlapfragmentenhancers;

overlapCoordinates(\%contactedcoords, \%enhancers, \%overlapfragmentenhancers);

my $nboe=keys %overlapfragmentenhancers;

print "Found ".$nboe." fragments that overlap with enhancers.\n";

print "Done.\n";

####################################################################################

print "Writing output...\n";

open(my $output, ">".$parameters{"pathOutput"});

print $output "IDBait\tIDContactedFragment\tIDPromoter\tIDEnhancer\n";

foreach my $idbait (keys %overlapbaitpromoters){
    foreach my $idfrag (@{$contacts{$idbait}}){
	if(exists $overlapfragmentenhancers{$idfrag}){
	    foreach my $idprom (@{$overlapbaitpromoters{$idbait}}){
		foreach my $idenh (@{$overlapfragmentenhancers{$idfrag}}){
		    print $output $idbait."\t".$idfrag."\t".$idprom."\t".$idenh."\n";
		}
	    }
	}
    }
    
}

print "Done.\n";

close($output);

####################################################################################
####################################################################################
