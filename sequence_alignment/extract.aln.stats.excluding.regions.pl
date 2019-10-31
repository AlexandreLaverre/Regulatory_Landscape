use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use strict;

#########################################################################

sub readAlignmentsMAF{
    my $pathin=$_[0];
    my $minlen=$_[1]; ## minimum alignment length (excluding gaps, for each species)
    my $aln=$_[2];
    
    my $input;

    my @ext=split("\\.",$pathin);
    my $nbext=@ext;
    my $ext=$ext[$nbext-1];

    if($ext eq "gz"){
	open($input,"zcat $pathin | ");
    }
    else{
	open($input,$pathin);
    }

    my $line=<$input>;
    my $firstchar=substr $line,0,1;

    while($firstchar eq "#"){
	$line=<$input>;
	$firstchar=substr $line,0,1;
    }
 
    my $currentscore="NA";
    
    my $indexaln=0;

    while($line){

	chomp $line;
	$firstchar=substr $line,0,1;
	
	if($firstchar eq "a"){
	    
	    $indexaln++;

	    my @s=split(" ",$line);
	    my $score=$s[1];
	    my @t=split("=",$score);
	    $score=$t[1]+0;
	    $currentscore=$score;

	    $aln->{$indexaln}={};
	}

	if($firstchar eq "s"){
	    my @s=split(" ",$line);
	    
	    my $spchr=$s[1];
	    my @t=split("\\.",$spchr); 
	    my $sp=$t[0];
	    my $chr=$t[1];

	    my $alnstart=$s[2]+0; ## starts at 0 
	    
	    my $ungappedlen=$s[3]+0;
	    my $strand=$s[4];
	    my $chrlen=$s[5]+0; ## length of the entire sequence that was aligned
	    my $sequence=$s[6];

	    my $start="NA";
	    my $increment="NA";
	    
	    if($strand eq "+"){
		$start=$alnstart;##included, starts at 0
		$increment=1;
	    } else{
		if($strand eq "-"){
		    $start=$chrlen-$alnstart-1; ## included, starts at 0
		    $increment=-1;
		} else{
		    print "Weird alignment strand ".$strand."\n";
		    exit(1);
		}
	    }
	    
	    $aln->{$indexaln}{$sp}={"start"=>$start, "increment"=>$increment, "strand"=>$strand,"sequence"=>$sequence,"chrlen"=>$chrlen,"ungappedlen"=>$ungappedlen};
	    
	}
		
	$line=<$input>;
    }

    close($input);
 

    ## remove alignments that are too short

    my @indexes=keys %{$aln};
    
    foreach my $index (@indexes){
	my $nbsp=keys %{$aln->{$index}};

	my $alllong=1;

	foreach my $sp (keys %{$aln->{$index}}){
	    my $len=$aln->{$index}{$sp}{"ungappedlen"};
	    
	    if($len<$minlen){
		$alllong=0;
		last;
	    }
	}

	if($alllong==0){
	    delete $aln->{$index};
	}
    }
   
    my $nbkept=keys %{$aln};

    # print "Kept ".$nbkept." alignments.\n";
}

################################################################################

sub readAlignmentsFasta{
    my $pathin=$_[0];
    my $minlen=$_[1]; ## minimum alignment length (excluding gaps, for each species)
    my $aln=$_[2];

    my @s=split("\\.",$pathin);
    my $ext=$s[-1];

    my $input;

    if($ext eq "gz"){
	open($input,"zcat $pathin |");
    }
    else{
	open($input, $pathin);
    }
    
    my $line=<$input>;

    my %fasta;
    
    while($line){
	my $b=substr $line,0,1;
	
	if($b eq ">"){
	    chomp $line;
	    my $id=substr $line,1;

	    my @s=split(" ",$id);
	    $id=$s[0];
	    
	    $fasta{$id}="";

	    $line=<$input>;
	    $b=substr $line,0,1;
	    
	    while($line && !($b eq ">")){
		chomp $line;
		$fasta{$id}.=$line;
		$line=<$input>;
		$b=substr $line,0,1;
	    }
	}
    }

    close($input);

    ## now format alignments like for MAF

    my $alllong=1;

    my $index=0;

    $aln->{$index}={};
    
    foreach my $sp (keys %fasta){
	my $sequence=$fasta{$sp};
	my $nbgaps = $sequence =~ tr/-//;
	my $totlen = length $fasta{$sp};
	my $ungappedlen=$totlen-$nbgaps;
	
	## these are stitched multiple alignments! always on the + strand, always starting at 0, covering the entire length of the sequence

	$aln->{$index}{$sp}={"start"=>0, "increment"=>1, "strand"=>"+", "sequence"=>$sequence, "chrlen"=>$ungappedlen, "ungappedlen"=>$ungappedlen};
	
	if($ungappedlen<$minlen){
	    $alllong=0;
	    last;
	}
    }
    	
    if($alllong==0){
	delete $aln->{$index};
    }
   
}

################################################################################

sub extractAlignmentStats{
    my $aln=$_[0];
    my $requiredspecies=$_[1];
    my $excludedhash=$_[2]; ## relative coordinates, for each species in the alignemnt
    my $alnstats=$_[3];
  
    my @allindexes=keys %{$aln};
    my $firstindex=$allindexes[0];
    my @splist=keys %{$aln->{$firstindex}};
    my $firstsp=$splist[0];

    my $nbungappedtotal=0;
    my $nbidenticaltotal=0;

    my $nbungappedfiltered=0;
    my $nbidenticalfiltered=0;

    my $alnlentotal=0;
    my $alnlenfiltered=0;
    
    foreach my $index (keys %{$aln}){
	## go over the alignment base by base

	my $alnlength=length $aln->{$index}{$firstsp}{"sequence"};

	$alnlentotal+=$alnlength;

	my %currentbase;

	foreach my $sp (keys %{$aln->{$index}}){
	    my $thisstart=$aln->{$index}{$sp}{"start"};
	    my $thisincrement=$aln->{$index}{$sp}{"increment"};
	    
	    $currentbase{$sp}=$thisstart-$thisincrement;
	}
	
	for(my $i=0; $i<$alnlength; $i++){
	    my %bases;

	    my $allungap=1;
	    my $isexcluded=0;
	    my $hasallsp=1;
	    
	    foreach my $sp (@{$requiredspecies}){
		if(exists $aln->{$index}{$sp}){
		    my $base=uc (substr $aln->{$index}{$sp}{"sequence"}, $i, 1);
		    
		    if($base eq "-"){
			$allungap=0;
		    }
		    else{
			$bases{$base}=1;
			$currentbase{$sp}+=$aln->{$index}{$sp}{"increment"}; ## not gapped, so we increment the position
			
			if(exists $excludedhash->{$sp}{$currentbase{$sp}}){
			    $isexcluded=1;
			}
		    }
		} else{
		    $hasallsp=0;
		}
	    }

	    if($isexcluded==0){
		$alnlenfiltered++;
	    }
	    
	    if($allungap==1 && $hasallsp==1){
		$nbungappedtotal++;

		if($isexcluded==0){
		    $nbungappedfiltered++;
		}

		my $nbdiffbases=keys %bases;

		if($nbdiffbases==1){
		    $nbidenticaltotal++;

		    if($isexcluded==0){
			$nbidenticalfiltered++;
		    }
		}
	    }
	}
    }

    $alnstats->{"ungappedtotal"}=$nbungappedtotal;
    $alnstats->{"identicaltotal"}=$nbidenticaltotal;

    $alnstats->{"ungappedfiltered"}=$nbungappedfiltered;
    $alnstats->{"identicalfiltered"}=$nbidenticalfiltered;

    $alnstats->{"alnlentotal"}=$alnlentotal;
    $alnstats->{"alnlenfiltered"}=$alnlenfiltered;
}

################################################################################

sub readPairs{
    my $pathin=$_[0];
    my $sp1=$_[1];
    my $sp2=$_[2];
    my $clusters=$_[3];
    
    open(my $input, $pathin);
    my $line=<$input>;
    chomp $line;
    my @s=split("\t", $line);
    my %header;
    for(my $i=0; $i<@s; $i++){
	$header{$s[$i]}=$i;
    }

    $line=<$input>;
    
    my $index=0;

    while($line){
	chomp $line;
	my @s=split("\t", $line);
	
	my $g1=$s[$header{"ID.".$sp1}];
	my $g2=$s[$header{"ID.".$sp2}];
	
	if($g1 ne "NA" && $g2 ne "NA"){
	    my @genes1=split(";", $g1);
	    my @genes2=split(";", $g2);
	    
	    $clusters->{$index}={$sp1=>[], $sp2=>[]};

	    push(@{$clusters->{$index}{$sp1}}, @genes1);
	    push(@{$clusters->{$index}{$sp2}}, @genes2);
		
	    $index++;
	}
	
	$line=<$input>;
    }
        
    close($input);
}

################################################################################

sub readExcludedRegions{
    my $pathin=$_[0];
    my $overlaps=$_[1];

    open(my $input, $pathin);
    
    my $line=<$input>; ## header

    chomp $line;
    my @s=split("\t", $line);
    my %header;

    for(my $i=0; $i<@s; $i++){
	$header{$s[$i]}=$i;
    }

    if((!exists $header{"ID"}) || (!exists $header{"overlap_ID"}) || (!exists $header{"start"}) || (!exists $header{"end"})){
	print "Weird header for excluded regions file: ".$line."\n";
	print "Expected header: ID      chr     start   end     overlap_ID\n";
	exit(1);
    }
    
    $line=<$input>;

    my %unordered;
        
    while($line){
	chomp $line;
	my @s=split("\t", $line);

	my $regionid=$s[$header{"ID"}];

	my @u=split(":", $regionid);
	my $regionstrand=$u[3];

	if($regionstrand ne "+" && $regionstrand ne "-"){
	    print "Weird strand for ".$regionid."\n";
	    exit(1);
	}
	
	my $regionstart=$s[$header{"start"}]+0;
	my $regionend=$s[$header{"end"}]+0;
	my $overlap=$s[$header{"overlap_ID"}];

	if($overlap ne "NA"){
	    if(exists $unordered{$regionid}){
		print "Weird! already saw ".$regionid." in excluded regions file!\n";
		exit(1);
	    }
	    
	    my @allov=split(",", $overlap);

	    foreach my $ov (@allov){
		my @t=split(":", $ov);
		my $absstart=$t[1]+0;
		my $absend=$t[2]+0;

		if($absend<$absstart){
		    print "Weird overlap coordinates, start>end, in ".$line."\n";
		    exit(1);
		}

		if($absend<=0 || $absstart<=0){
		    print "Weird overlap coordinates, they should be strictly positive: ".$absstart." ".$absend."\n";
		    exit(1);
		}

		if($absstart<$regionstart){
		    $absstart=$regionstart;	
		}
		
		if($absend<$regionstart){
		    print "Weird overlap coordinates, they should be within the initial region coordinates.\n";
		    print "region coordinates: ".$regionstart." ".$regionend."\n";
		    print "overlap coordinates: ".$absstart." ".$absend."\n";
		    exit(1);
		}
		
		
		if($absend>$regionend){
		    $absend=$regionend;	
		}
		
		if($absstart>$regionend){
		    print "Weird overlap coordinates, they should be within the initial region coordinates.\n";
		    print "region coordinates: ".$regionstart." ".$regionend."\n";
		    print "overlap coordinates: ".$absstart." ".$absend."\n";
		    exit(1);
		}
		
		
		## convert to relative coordinates, with respect to the sequence that is given to TBA or PECAN (reverse complemented if on strand -)

		my $relstart=$absstart-$regionstart; ## starts now at 0
		my $relend=$absend-$regionstart; ## starts now at 0, is included

		if($regionstrand eq "-"){
		    $relstart=$regionend-$absend; ## starts now at 0
		    $relend=$regionend-$absstart; ## starts now at 0, is included
		}

		my $maxpos=$regionend-$regionstart;

		if($relstart<0 || $relend<0 || $relstart>$maxpos || $relend>$maxpos || $relstart>$relend){
		    print "Weird relative coordinates: ".$relstart." ".$relend."\n";
		    exit(1);
		}

		#print "Saw excluded coordinates: ".$relstart." ".$relend." for ".$regionid."\n";
		
		if(exists $unordered{$regionid}){
		    if(exists $unordered{$regionid}{$relstart}){
			if($relend > $unordered{$regionid}{$relstart}){
			    $unordered{$regionid}{$relstart}=$relend;
			}
		    } else{
			$unordered{$regionid}{$relstart}=$relend;
		    }
		} else{
		    $unordered{$regionid}={$relstart=>$relend};
		}
	    }
	}
	
	$line=<$input>;
    }
    
    close($input);

    ## now order the overlaps

    foreach my $regionid (keys %unordered){
	$overlaps->{$regionid}={"start"=>[], "end"=>[]};
	
	my @uniquestart=keys %{$unordered{$regionid}};
	my @sortedstart=sort @uniquestart;

	my $nbov=@sortedstart;

	my $currentstart=$sortedstart[0];
	my $currentend=$unordered{$regionid}{$currentstart};

	for(my $i=1; $i<$nbov; $i++){
	    my $thisstart=$sortedstart[$i];
	    my $thisend=$unordered{$regionid}{$thisstart};

	    if($thisstart>=$currentstart && $thisstart<=($currentend+1)){
		if($thisend>$currentend){
		    $currentend=$thisend;
		}
	    } else{
		push(@{$overlaps->{$regionid}{"start"}}, $currentstart);
		push(@{$overlaps->{$regionid}{"end"}}, $currentend);

		$currentstart=$thisstart;
		$currentend=$thisend;
	    }
	}

	## don't forget the last block

	push(@{$overlaps->{$regionid}{"start"}}, $currentstart);
	push(@{$overlaps->{$regionid}{"end"}}, $currentend);
    }
}

################################################################################

sub makeCoordinatesHash{
    my $startpos=$_[0];
    my $endpos=$_[1];
    my $hashcoords=$_[2];
    
    my $nbstart=@{$startpos};
    my $nbend=@{$endpos};

    if($nbstart!=$nbend){
	print "Weird! different array sizes for start and end.\n";
	exit(1);
    }

    for(my $i=0; $i<$nbstart; $i++){
	my $start=$startpos->[$i];
	my $end=$endpos->[$i];

	for(my $pos=$start; $pos<=$end; $pos++){
	    
	    if(exists $hashcoords->{$pos}){
		"Weird, already seen position when making coordinates dictionary. The regions are supposed to be non-overlapping.\n";
		exit(1);
	    }
	    
	    $hashcoords->{$pos}=1;
	}
    }
}

################################################################################

sub printHelp{
    my $parnames=$_[0];
    my $parvalues=$_[1];
    
    print "\n";
    print "This script extracts exon alignment stats from TBA or PECAN alignments. \n";
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
$parameters{"species1"}="NA";
$parameters{"species2"}="NA";
$parameters{"pathPairs"}="NA";
$parameters{"pathExcludedRegionsSpecies1"}="NA";
$parameters{"pathExcludedRegionsSpecies2"}="NA";
$parameters{"dirAlignments"}="NA";
$parameters{"suffixAlignments"}="NA";
$parameters{"format"}="NA";
$parameters{"minAlignmentLength"}="NA";
$parameters{"pathOutput"}="NA";

my %defaultvalues;
my @defaultpars=("species1", "species2", "pathPairs", "pathExcludedRegionsSpecies1", "pathExcludedRegionsSpecies2",  "dirAlignments", "suffixAlignments", "format", "minAlignmentLength", "pathOutput");

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

my $sp1=$parameters{"species1"};
my $sp2=$parameters{"species2"};
my @species=($sp1, $sp2);

#####################################################################

print "Reading pairs...\n";

my %clusters;
readPairs($parameters{"pathPairs"}, $sp1, $sp2, \%clusters);

my $nbclust=keys %clusters;

print "Found ".$nbclust." non-empty pairs.\n";
print "Done.\n";

#####################################################################

print "Reading excluded regions...\n";

my %excludedregions1;

if(-e $parameters{"pathExcludedRegionsSpecies1"}){
    readExcludedRegions($parameters{"pathExcludedRegionsSpecies1"}, \%excludedregions1);
    
    my $nbreg1=keys %excludedregions1;
    
    print "Found ".$nbreg1." fragments with regions to be excluded for ".$sp1."\n";
} else{
    print "We don't exclude anything for ".$sp1."\n";
}


my %excludedregions2;

if(-e $parameters{"pathExcludedRegionsSpecies2"}){
    readExcludedRegions($parameters{"pathExcludedRegionsSpecies2"}, \%excludedregions2);
    
    my $nbreg2=keys %excludedregions2;
    
    print "Found ".$nbreg2." fragments with regions to be excluded for ".$sp2."\n";
} else{
    print "We don't exclude anything for ".$sp2."\n";
}
    
print "Done.\n";

#####################################################################

print "Reading alignments and writing output...\n";

my $format=$parameters{"format"};

if($format ne "MAF" && $format ne "MFA" && $format ne "fasta"){
    print "Unsupported alignment format: ".$format."\n";
    print "Accepted formats: MAF, MFA (multi-fasta), fasta.\n";
    exit(1);
}

my $minlength=$parameters{"minAlignmentLength"}+0;

print "minimum length: ".$minlength."\n";

open(my $output, ">".$parameters{"pathOutput"});
print $output "ID.".$sp1."\tID.".$sp2."\tNbExcludedBases1\tNbExcludedBases2\tTotalUngappedLength\tTotalIdenticalLength\tFilteredUngappedLength\tFilteredIdenticalLength\tTotalAlignmentLength\tFilteredAlignmentLength\n";

foreach my $idclust (keys %clusters){
    foreach my $gene1 (@{$clusters{$idclust}{$sp1}}){
	my @s1=split(":", $gene1);
	my $start1=$s1[1]+0;
	my $end1=$s1[2]+0;
	my $seqlen1=$end1-$start1;
	
	foreach my $gene2 (@{$clusters{$idclust}{$sp2}}){
	    my @s2=split(":", $gene2);
	    my $start2=$s2[1]+0;
	    my $end2=$s2[2]+0;
	    my $seqlen2=$end2-$start2;
	    
	    my $pathAln=$parameters{"dirAlignments"}."/".$gene1.".".$gene2.$parameters{"suffixAlignments"};
	    
	    if(-e $pathAln){
				
		my %aln;

		if($format eq "MAF"){
		    readAlignmentsMAF($pathAln, $minlength, \%aln);
		}
		
		if($format eq "MFA" || $format eq "fasta"){
		    readAlignmentsFasta($pathAln, $minlength, \%aln);

		    my $nbaln=keys %aln;

		    if($nbaln==1){
			my $ungappedlen1=$aln{0}{$sp1}{"chrlen"};
			my $ungappedlen2=$aln{0}{$sp2}{"chrlen"};
			
			if($seqlen1!=$ungappedlen1){
			    print "Weird! different sequence length for ".$sp1." sequence ".$gene1." length ".$seqlen1." ungapped alignment length ".$ungappedlen1."\n";
			    exit(1);
			}
			
			if($seqlen2!=$ungappedlen2){
			    print "Weird! different sequence length for ".$sp2." sequence ".$gene2." length ".$seqlen2." ungapped alignment length ".$ungappedlen2."\n";
			    exit(1);
			}
		    } else{
			print $nbaln." alignments found for ".$gene1." and ".$gene2."\n";
		    }
		}
		
		my %excludedregionshash;
		$excludedregionshash{$sp1}={};
		$excludedregionshash{$sp2}={};

		if(exists $excludedregions1{$gene1}){
		    makeCoordinatesHash($excludedregions1{$gene1}{"start"}, $excludedregions1{$gene1}{"end"}, $excludedregionshash{$sp1});
		}

		if(exists $excludedregions2{$gene2}){
		    makeCoordinatesHash($excludedregions2{$gene2}{"start"}, $excludedregions2{$gene2}{"end"}, $excludedregionshash{$sp2});
		}
		
		my $nbexcluded1=keys %{$excludedregionshash{$sp1}};
		my $nbexcluded2=keys %{$excludedregionshash{$sp2}};
		    
		my %alnstats;
		extractAlignmentStats(\%aln, \@species, \%excludedregionshash, \%alnstats);
		
		print $output $gene1."\t".$gene2."\t".$nbexcluded1."\t".$nbexcluded2."\t".$alnstats{"ungappedtotal"}."\t".$alnstats{"identicaltotal"}."\t".$alnstats{"ungappedfiltered"}."\t".$alnstats{"identicalfiltered"}."\t".$alnstats{"alnlentotal"}."\t".$alnstats{"alnlenfiltered"}."\n";

		
	    } else{
		print "Weird! cannot find ".$pathAln." for ".$gene1." and ".$gene2."\n";
	    }
	}
    }
}

close($output);

print "Done.\n";

#####################################################################
#####################################################################
