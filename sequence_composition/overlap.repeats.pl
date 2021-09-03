use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use strict;

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

	if(exists $refcoords->{$chr}){
	    push(@{$refcoords->{$chr}{"id"}}, $id);
	    push(@{$refcoords->{$chr}{"start"}}, $start);
	    push(@{$refcoords->{$chr}{"end"}}, $end);
	    push(@{$refcoords->{$chr}{"exon_total_bp"}}, "NA");
	    push(@{$refcoords->{$chr}{"exon_repeat_bp"}}, "NA");
	    push(@{$refcoords->{$chr}{"repeat_total_bp"}}, "NA");
	} else{
	    $refcoords->{$chr}={"id"=>[$id], "start"=>[$start], "end"=>[$end], "exon_total_bp"=>["NA"], "exon_repeat_bp"=>["NA"], "repeat_total_bp"=>["NA"]};
	}

	$nbel++;
	 	
	$line=<$input>;
    }
    
    close($input);

    print "Found ".$nbel." elements.\n";
}

##########################################################################

sub readRepeatMasker{
    my $pathin=$_[0];
    my $coords=$_[1];

    my @s=split("\\.", $pathin);
    my $ext=$s[-1];

    my $input;

    if($ext eq "gz"){
	open($input, "zcat $pathin |");	
    } else{
	open($input, $pathin);	
    }
    
    my $line=<$input>;

    while($line){
	chomp $line;
	my @s=split("\t",$line);
	
	my $chr=$s[5];
	my $start=$s[6]+1; ## now 1-based
	my $end=$s[7]+0; ## now 1-based, included

	if(exists $coords->{$chr}){
	    push(@{$coords->{$chr}{"start"}}, $start);
	    push(@{$coords->{$chr}{"end"}}, $end);
	} else{
	    $coords->{$chr}={"start"=>[$start], "end"=>[$end]};
	}
	
	$line=<$input>;
    }
    close($input);
}

##########################################################################

sub readExonBlocks{
    my $pathin=$_[0];
    my $refblocks=$_[1];

    open(my $input, $pathin);
    my $line=<$input>;
    
    while($line){
	chomp $line;
	my @s=split("\t",$line);

	## gene id, exon id, strand ignored
	my $chr=$s[2];
	my $start=$s[3];
	my $end=$s[4]+0; 

	if(exists $refblocks->{$chr}){
	    push(@{$refblocks->{$chr}{"start"}},$start);
	    push(@{$refblocks->{$chr}{"end"}},$end);
	}
	else{
	    $refblocks->{$chr}={"start"=>[$start],"end"=>[$end]};
	}
	    
	$line=<$input>;
    }
    close($input);
}

##########################################################################

sub makeHashCoords{
    my $coords=$_[0];
    my $unfold=$_[1];
    
    ## single chromosome, start end
    my $nbpos=@{$coords->{"start"}};

    for(my $i=0; $i<$nbpos; $i++){
	my $start=${$coords->{"start"}}[$i];
	my $end=${$coords->{"end"}}[$i];

	for(my $k=$start; $k<=$end; $k++){
	    $unfold->{$k}=1;
	}
    }
}

##########################################################################

sub computeOverlap{
    my $coords=$_[0]; ## single chromosome
    my $exons=$_[1];
    my $repeats=$_[2];

    my $nbpos=@{$coords->{"start"}};

    for(my $i=0; $i<$nbpos; $i++){
	my $start=${$coords->{"start"}}[$i];
	my $end=${$coords->{"end"}}[$i];

	my $nbrep=0;
	my $nbex=0;
	my $nbexrep=0;

	for(my $k=$start; $k<=$end; $k++){
	    
	    if(exists $repeats->{$k}){
		$nbrep++;
	    }
	    
	    if(exists $exons->{$k}){
		$nbex++;
		
		if(exists $repeats->{$k}){
		    $nbexrep++;
		}
	    } 
	}

	${$coords->{"exon_total_bp"}}[$i]=$nbex;
	${$coords->{"repeat_total_bp"}}[$i]=$nbrep;
	${$coords->{"exon_repeat_bp"}}[$i]=$nbexrep;
    }
}
    
##########################################################################

sub printHelp{

    my $parnames=$_[0];
    my $parvalues=$_[1];
    
    print "\n";
    print "This script computes overlap with repeats for restriction fragments and enhancers. \n";
    print "\n";
    print "Options:\n";
    
    foreach my $par (@{$parnames}){
	print "--".$par."  [  default value: ".$parvalues->{$par}."  ]\n";
    }
    print "\n";
}

######################################################################
######################################################################

my %parameters;

$parameters{"pathCoords"}="NA";
$parameters{"coordConvention"}="NA";
$parameters{"pathExonBlocks"}="NA";
$parameters{"pathRepeatMasker"}="NA";
$parameters{"pathOutput"}="NA";


my %defaultvalues;
my @defaultpars=("pathCoords", "coordConvention", "pathExonBlocks", "pathRepeatMasker", "pathOutput");

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

print "Reading exon coordinates...\n";

my %exoncoords;
readExonBlocks($parameters{"pathExonBlocks"}, \%exoncoords);

print "Done.\n";

##############################################################

print "Reading RepeatMasker coordinates...\n";

my %repeatcoords;
readRepeatMasker($parameters{"pathRepeatMasker"}, \%repeatcoords);

print "Done.\n";

##############################################################

print "Reading genomic coordinates...\n";
my %elements;

my $convention=$parameters{"coordConvention"};

print "coordinate convention: ".$convention."\n";

readCoords($parameters{"pathCoords"}, $convention, \%elements);
print "Done.\n";

##############################################################

print "Computing overlap and writing output \n";

open(my $output,">".$parameters{"pathOutput"});

print $output "ID\tChr\tStart\tEnd\TotalLength\tExonicLength\tRepeatLength\tExonRepeatLength\n";

foreach my $chr (keys %elements){

    print "chr ".$chr."\n";

    print "formatting coordinates for exons...\n";
    
    my %hashexons;

    if(exists $exoncoords{$chr}){
	makeHashCoords($exoncoords{$chr}, \%hashexons);
    }

    print "formatting coordinates for repeats...\n";
    
    my %hashrepeats;

    if(exists $repeatcoords{$chr}){
	makeHashCoords($repeatcoords{$chr}, \%hashrepeats);
    }

    print "computing overlap...\n";

    computeOverlap($elements{$chr}, \%hashexons, \%hashrepeats);

    my $nbpos=@{$elements{$chr}{"start"}};

    for(my $i=0; $i<$nbpos; $i++){
	my $id=${$elements{$chr}{"id"}}[$i];
	my $start=${$elements{$chr}{"start"}}[$i];
	my $end=${$elements{$chr}{"end"}}[$i];
	my $exonbp=${$elements{$chr}{"exon_total_bp"}}[$i];
	my $repeatbp=${$elements{$chr}{"repeat_total_bp"}}[$i];
	my $repexonbp=${$elements{$chr}{"exon_repeat_bp"}}[$i];

	my $totalbp=$end-$start+1;
	
	print $output $id."\t".$chr."\t".$start."\t".$end."\t".$totalbp."\t".$exonbp."\t".$repeatbp."\t".$repexonbp."\n";
    }
}
    
close($output);

print "Done.\n";


##############################################################
