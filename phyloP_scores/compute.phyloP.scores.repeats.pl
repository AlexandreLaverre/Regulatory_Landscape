use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use strict;

##############################################################

sub readScores{
    my $pathin=$_[0];
    my $refphast=$_[1];
    
    my @s=split("\\.",$pathin);
    my $nbs=@s;
    my $ext=$s[$nbs-1];

    my $input;
    
    if($ext eq "gz"){
	open($input,"zcat $pathin |");
    }
    else{
	open($input,$pathin);
    }
    
    my $line=<$input>;
    my $currentpos="NA";
    
    while($line){
	my $first=substr $line,0,1;
	
	if($first eq "f"){
	    my @s=split(" ",$line);
	    my $start=$s[2];
	    my @t=split("=",$start);
	    $currentpos=$t[1]+0;
	    $line=<$input>;
	    next;
	}
	else{
	    chomp $line;
	    my $val=$line+0.0;

	    $refphast->{$currentpos}=$val;
	    
	    $currentpos++;
	    $line=<$input>;
	}
    }
    
    close($input);
}

##############################################################

sub computeScore{
    my $coords=$_[0];
    my $phast=$_[1];
    my $maskedcoords=$_[2];
    my $repcoords=$_[3];
    my $chr=$_[4];

    my $nbmotifs=@{$coords->{$chr}{"start"}};

    for(my $i=0; $i<$nbmotifs; $i++){
	my $start=${$coords->{$chr}{"start"}}[$i];
	my $end=${$coords->{$chr}{"end"}}[$i];

	my $sumscoreRepeatOnly=0;
	my $nbanalyzedRepeatOnly=0;
	my $nbcoveredRepeatOnly=0;

	my $sumscoreRepeatFree=0;
	my $nbanalyzedRepeatFree=0;
	my $nbcoveredRepeatFree=0;

	for(my $j=$start; $j<=$end; $j++){
	    if(!(exists $maskedcoords->{$j})){

		if(exists $repcoords->{$j}){
		    $nbanalyzedRepeatOnly++;
		} else{
		    $nbanalyzedRepeatFree++;
		}
		
		if(exists $phast->{$j}){
		    if(exists $repcoords->{$j}){
			$nbcoveredRepeatOnly++;
			$sumscoreRepeatOnly+=$phast->{$j};
		    } else{
			$nbcoveredRepeatFree++;
			$sumscoreRepeatFree+=$phast->{$j};
		    }
		}
	    }
	}
	
	if($nbcoveredRepeatOnly>0){
	    ${$coords->{$chr}{"scoreRepeatOnly"}}[$i]=($sumscoreRepeatOnly+0.0)/($nbcoveredRepeatOnly+0.0);
	}

	if($nbcoveredRepeatFree>0){
	    ${$coords->{$chr}{"scoreRepeatFree"}}[$i]=($sumscoreRepeatFree+0.0)/($nbcoveredRepeatFree+0.0);
	} 
	
	${$coords->{$chr}{"coveredbasesRepeatOnly"}}[$i]=$nbcoveredRepeatOnly+0.0;
	${$coords->{$chr}{"analyzedbasesRepeatOnly"}}[$i]=$nbanalyzedRepeatOnly+0.0;

	${$coords->{$chr}{"coveredbasesRepeatFree"}}[$i]=$nbcoveredRepeatFree+0.0;
	${$coords->{$chr}{"analyzedbasesRepeatFree"}}[$i]=$nbanalyzedRepeatFree+0.0;
    }
    
}

##############################################################

sub readExonBlocks{
    my $pathin=$_[0];
    my $chr=$_[1];
    my $hashcoords=$_[2];

    open(my $input, $pathin);
    my $line=<$input>; ## header
    
    my $ensembl=$chr;
    my $prefix=substr $chr, 0, 3;

    if($prefix eq "chr"){
	$ensembl=substr $chr, 3;
    }
    
    my $ucsc="chr".$ensembl;
    
    while($line){
	chomp $line;
	my @s=split("\t",$line);
	
	my $thischr=$s[2];

	if($thischr eq $ensembl || $thischr eq $ucsc){
	    my $start=$s[3]+0; ## 1-based
	    my $end=$s[4]+0; ## 1-based, included

	    for(my $pos=$start; $pos<=$end; $pos++){
		$hashcoords->{$pos}=1;
	    }
	}
	    
	$line=<$input>;
    }
    close($input);
}

##########################################################################

sub readRepeatMasker{
    my $pathin=$_[0];
    my $chr=$_[1];
    my $hashcoords=$_[2];

    my @s=split("\\.", $pathin);
    my $ext=$s[-1];

    my $input;

    if($ext eq "gz"){
	open($input, "zcat $pathin |");	
    } else{
	open($input, $pathin);	
    }
    
    my $line=<$input>;

    my $ensembl=$chr;
    my $prefix=substr $chr, 0, 3;

    if($prefix eq "chr"){
	$ensembl=substr $chr, 3;
    }
    
    my $ucsc="chr".$ensembl;
    
    while($line){
	chomp $line;
	my @s=split("\t",$line);
	
	my $thischr=$s[5];

	if($thischr eq $ensembl || $thischr eq $ucsc){
	    my $start=$s[6]+1; ## now 1-based
	    my $end=$s[7]+0; ## now 1-based, included

	    for(my $pos=$start; $pos<=$end; $pos++){
		$hashcoords->{$pos}=1;
	    }
	}
	    
	$line=<$input>;
    }
    close($input);
}

##########################################################################

sub readCoords{
    my $pathin=$_[0];
    my $coordtype=$_[1];
    my $refcoords=$_[2];

    my $offsetstart=0;
    my $offsetend=0;

    if($coordtype eq "0_open_end"){
	$offsetstart=1;
	$offsetend=0;
    } else{
	if($coordtype eq "1_closed_end"){
	    $offsetstart=0;
	    $offsetend=0;
	} else{
	    print "Weird! unknown coordinate convention: ".$coordtype."\n";
	    exit(1);
	}
    }

    print "Adding offset ".$offsetstart." for start coordinates.\n";
    print "Adding offset ".$offsetend." for end coordinates.\n";
    
    open(my $input, $pathin);
    my $line=<$input>;

    my $nbel=0;
    
    while($line){
	chomp $line;
	my @s=split("\t",$line);
	
	my $chr=$s[0];
	my $start=$s[1]+$offsetstart;
	my $end=$s[2]+$offsetend; 
	my $id=$s[3];
	my $strand="NA";

	$refcoords->{$id}={"chr"=>$chr, "strand"=>$strand, "start"=>[$start], "end"=>[$end]};

	$nbel++;
	 	
	$line=<$input>;
    }
    
    close($input);

    print "Found ".$nbel." elements.\n";
}

##########################################################################

sub orderCoords{
    my $refblocks=$_[0];
    my $refordered=$_[1];
    
    my %hashstart;
    
    foreach my $gene (keys %{$refblocks}){
	my $nbblocks=@{$refblocks->{$gene}{"start"}};
	my $chr=$refblocks->{$gene}{"chr"};
	my $strand=$refblocks->{$gene}{"strand"};

	for(my $i=0;$i<$nbblocks;$i++){
	    my $tb=${$refblocks->{$gene}{"start"}}[$i];
	    my $te=${$refblocks->{$gene}{"end"}}[$i];
	
	    if(exists $hashstart{$chr}){
		if(exists $hashstart{$chr}{$tb}){
		    if(exists $hashstart{$chr}{$tb}{$te}){
			if(exists $hashstart{$chr}{$tb}{$te}{$strand}){
			    push(@{$hashstart{$chr}{$tb}{$te}{$strand}},$gene);
			}
			else{
			    $hashstart{$chr}{$tb}{$te}{$strand}=[$gene];
			}
		    }
		    else{
			$hashstart{$chr}{$tb}{$te}={$strand=>[$gene]};
		    }
		}
		else{
		    $hashstart{$chr}{$tb}={$te=>{$strand=>[$gene]}};
		}
	    }
	    else{
		$hashstart{$chr}={$tb=>{$te=>{$strand=>[$gene]}}};
	    }
	}
	
    }

    foreach my $chr (keys %hashstart){
	$refordered->{$chr}={"start"=>[],"end"=>[],"strand"=>[],"gene"=>[],"scoreRepeatOnly"=>[], "scoreRepeatFree"=>[], "coveredbasesRepeatOnly"=>[], "coveredbasesRepeatFree"=>[], "analyzedbasesRepeatOnly"=>[], "analyzedbasesRepeatFree"=>[]};

	my @uniquestart=keys %{$hashstart{$chr}};
	
	my @sortedstart = sort {$a <=> $b} @uniquestart;

	foreach my $start (@sortedstart){
	    my @uniqueend=keys %{$hashstart{$chr}{$start}};
	    my @sortedend = sort {$a <=> $b} @uniqueend;

	    foreach my $end (@sortedend){
		foreach my $strand (keys %{$hashstart{$chr}{$start}{$end}}){
		    foreach my $gene (@{$hashstart{$chr}{$start}{$end}{$strand}}){
			push(@{$refordered->{$chr}{"start"}},$start);
			push(@{$refordered->{$chr}{"end"}},$end);
			push(@{$refordered->{$chr}{"strand"}},$strand);
			push(@{$refordered->{$chr}{"gene"}},$gene);
			
			push(@{$refordered->{$chr}{"scoreRepeatOnly"}},"NA");
			push(@{$refordered->{$chr}{"coveredbasesRepeatOnly"}},0);
			push(@{$refordered->{$chr}{"analyzedbasesRepeatOnly"}},0);

			push(@{$refordered->{$chr}{"scoreRepeatFree"}},"NA");
			push(@{$refordered->{$chr}{"coveredbasesRepeatFree"}},0);
			push(@{$refordered->{$chr}{"analyzedbasesRepeatFree"}},0);
		
		    }
		} 
	    }
	} 
    }

}

##########################################################################

sub printHelp{

    my $parnames=$_[0];
    my $parvalues=$_[1];
    
    print "\n";
    print "This script computes phastCons or phyloP scores for genomic coordinates. \n";
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

$parameters{"pathCoords"}="NA";
$parameters{"coordConvention"}="NA";
$parameters{"pathMaskExonBlocks"}="NA";
$parameters{"pathScores"}="NA";
$parameters{"pathRepeatMasker"}="NA";
$parameters{"chr"}="NA";
$parameters{"pathOutput"}="NA";


my %defaultvalues;
my @defaultpars=("pathCoords", "coordConvention", "pathMaskExonBlocks","pathRepeatMasker", "pathScores", "chr", "pathOutput");

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


##############################################################
##############################################################

my $chr=$parameters{"chr"};

print "Analyzing chromosome ".$chr."\n";

my %maskedcoords;

##############################################################

my $mask="no";

if(-e $parameters{"pathMaskExonBlocks"}){
    $mask="yes";

    print "Reading masked coordinates...\n";
    
    readExonBlocks($parameters{"pathMaskExonBlocks"}, $chr, \%maskedcoords);
        
    my $nbmask=keys %maskedcoords;
    
    if($nbmask==0){
	print "Weird! there is nothing to mask on this chromosome...\n";
	exit(1);
    }
    
    print "Masking ".$nbmask. " positions.\n";
    
    print "Done.\n";
}

##############################################################

print "Reading RepeatMasker coordinates...\n";

my %repeatcoords;
readRepeatMasker($parameters{"pathRepeatMasker"}, $chr, \%repeatcoords);

my $nbrep=keys %repeatcoords;

if($nbrep==0){
    print "Weird! there are no repeats on this chromosome...\n";
    exit(1);
}

print "Masking ".$nbrep. " positions.\n";

print "Done.\n";

##############################################################

print "Reading genomic coordinates...\n";
my %elements;

my $convention=$parameters{"coordConvention"};

print "coordinate convention: ".$convention."\n";

readCoords($parameters{"pathCoords"}, $convention, \%elements);
print "Done.\n";

print "Ordering regions...\n";
my %orderedelements;
orderCoords(\%elements, \%orderedelements);
print "Done\n";

print "Computing score and writing output \n";

open(my $output,">".$parameters{"pathOutput"});
print $output "ID\tChr\tStart\tEnd\tScoreRepeatOnly\tCoveredLengthRepeatOnly\tAnalyzedLengthRepeatOnly\tScoreRepeatFree\tCoveredLengthRepeatFree\tAnalyzedLengthRepeatFree\n";

my $path=$parameters{"pathScores"};


if(-e $path){
    if(exists $orderedelements{$chr}){
	
	my %phastCons;
	readScores($path,\%phastCons);
	
	computeScore(\%orderedelements,\%phastCons,\%maskedcoords, \%repeatcoords, $chr);
        
	my $nborderedelements=@{$orderedelements{$chr}{"start"}};
	
	my %scoregenes;
	
	
	for(my $i=0;$i<$nborderedelements;$i++){
	    
	    my $gene=${$orderedelements{$chr}{"gene"}}[$i];
	    
	    print $output ${$orderedelements{$chr}{"gene"}}[$i]."\t".$chr."\t".${$orderedelements{$chr}{"start"}}[$i]."\t".${$orderedelements{$chr}{"end"}}[$i]."\t".${$orderedelements{$chr}{"scoreRepeatOnly"}}[$i]."\t".${$orderedelements{$chr}{"coveredbasesRepeatOnly"}}[$i]."\t".${$orderedelements{$chr}{"analyzedbasesRepeatOnly"}}[$i]."\t".${$orderedelements{$chr}{"scoreRepeatFree"}}[$i]."\t".${$orderedelements{$chr}{"coveredbasesRepeatFree"}}[$i]."\t".${$orderedelements{$chr}{"analyzedbasesRepeatFree"}}[$i]."\n";
	    
	}
    }
}

close($output);

print "Done.\n";


##############################################################
