use List::Util qw(first max maxstr min minstr reduce shuffle sum);

use strict;



##############################################################

##############################################################



sub readTSS{

    my $pathin=$_[0];

    my $orderedtss=$_[1];



    open(my $input, $pathin);

    my $line=<$input>; # header

    $line=<$input>; # first actual line

    

    my %unordered;



    while($line){

	chomp $line;

	my @s=split("\t", $line);

	my $geneid=$s[0];

	my $chr=$s[3];

	my $start=$s[4]+0;

	my $end=$s[5]+0;

	my $strand=$s[6];



	my $tss;



	if($strand eq "1"){

	    $tss=$start;

	} else{

	    if($strand eq "-1"){

		$tss=$end;

	    }

	    else{

		print "Weird strand!! ".$line."\n";

		exit(1);

	    }

	}



	if(exists $unordered{$chr}){

	    if(exists $unordered{$chr}{$tss}){

		push(@{$unordered{$chr}{$tss}{"gene"}}, $geneid);

		push(@{$unordered{$chr}{$tss}{"strand"}}, $strand);

	    }

	    else{

		$unordered{$chr}{$tss}={"gene"=>[$geneid], "strand"=>[$strand]};

	    }

	}

	else{

	    $unordered{$chr}={$tss=>{"gene"=>[$geneid], "strand"=>[$strand]}};

	}



	$line=<$input>;	

    }



    foreach my $chr (keys %unordered){

	$orderedtss->{$chr}={"position"=>[], "strand"=>[], "gene"=>[]};



	my @positions=keys %{$unordered{$chr}};

	my @sortedpositions=sort {$a<=>$b} @positions;



	foreach my $pos (@sortedpositions){

	    my $nbg=@{$unordered{$chr}{$pos}{"gene"}};



	    for(my $i=0; $i<$nbg; $i++){

		push(@{$orderedtss->{$chr}{"position"}}, $pos);

		push(@{$orderedtss->{$chr}{"strand"}}, ${$unordered{$chr}{$pos}{"strand"}}[$i]);

		push(@{$orderedtss->{$chr}{"gene"}}, ${$unordered{$chr}{$pos}{"gene"}}[$i]);

	    }

	}

	

    }



    close($input);

}



##############################################################



sub readBaits{

    my $pathin=$_[0];

    my $acceptablechromo=$_[1];

    my $orderedbaits=$_[2];



    my $input;

    my @s=split("\\.",$pathin);

    my $ext=$s[-1];

    

    if($ext eq "gz"){

	open($input, "zcat $pathin |");

    }

    else{

	open($input, $pathin);

    }



    my $line=<$input>;  # header

    $line=<$input>; # first actual line

    

    my %unordered;



    while($line){

	chomp $line;

	my @s=split("\t", $line);



	my $chr=$s[0];

	

	if(exists $acceptablechromo->{$chr}){

	    $chr=substr $chr, 3;

	

	    my $start=$s[1]+0;

	    my $end=$s[2]+0;

	    

	    if(exists $unordered{$chr}){

		if(exists $unordered{$chr}{$start}){

		    push(@{$unordered{$chr}{$start}{"end"}}, $end);

		}

		else{

		    $unordered{$chr}{$start}={"end"=>[$end]};

		}

	    }

	    else{

		$unordered{$chr}={$start=>{"end"=>[$end]}};

	    }

	}

	else{

	    print "Discarding ".$line.", unknown chromosome\n";

	}

		

	$line=<$input>;	

    }



    foreach my $chr (keys %unordered){

	$orderedbaits->{$chr}={"start"=>[], "end"=>[]};

	

	my @uniquestart=keys %{$unordered{$chr}};

	my @sortedstart=sort {$a<=>$b} @uniquestart;

	

	foreach my $pos (@sortedstart){

	    my $nbe=@{$unordered{$chr}{$pos}{"end"}};

	    

	    for(my $i=0; $i<$nbe; $i++){

		push(@{$orderedbaits->{$chr}{"start"}}, $pos);

		push(@{$orderedbaits->{$chr}{"end"}}, ${$unordered{$chr}{$pos}{"end"}}[$i]);

	    }

	}

    }

    

    close($input);

}



##############################################################



sub readGeneNames{

    my $pathin=$_[0];

    my $genenames=$_[1];



    open(my $input, $pathin);

    my $line=<$input>;

    $line=<$input>;



    while($line){

	chomp $line;

	my @s=split("\t", $line);

	my $id=$s[0];

	my $name=$s[1];



	$genenames->{$id}=$name;



	$line=<$input>;

    }



    close($input);

}



##############################################################



sub overlapCoords{

    my $baits=$_[0];

    my $tss=$_[1];

    my $maxdist=$_[2];

    my $overlap=$_[3];



    foreach my $chr (keys %{$baits}){

	if(exists $tss->{$chr}){

	    my $nbbaits=@{$baits->{$chr}{"start"}};

	    my $nbtss=@{$tss->{$chr}{"position"}};



	    my $firstj=0;



	    for(my $i=0; $i<$nbbaits; $i++){

		my $startbait=${$baits->{$chr}{"start"}}[$i];

		my $endbait=${$baits->{$chr}{"end"}}[$i];

		my $keybait=$chr.",".$startbait.",".$endbait;



		my $j=$firstj;



		while($j<$nbtss && ${$tss->{$chr}{"position"}}[$j]<($startbait-$maxdist)){

		    $j++;

		}

		

		$firstj=$j;

		

		while($j<$nbtss && ${$tss->{$chr}{"position"}}[$j]<=($endbait+$maxdist)){

		   my $pos=${$tss->{$chr}{"position"}}[$j];

		    

		   if($pos>=($startbait-$maxdist) && $pos<=($endbait+$maxdist)){ ## double check, probably not necessary

		       my $gene=${$tss->{$chr}{"gene"}}[$j];

		       my $strand=${$tss->{$chr}{"strand"}}[$j];

		    

		       if(exists $overlap->{$keybait}){

			   if(exists $overlap->{$keybait}{$gene}){

			       $overlap->{$keybait}{$gene}{"tss"}{$pos}=1;

			       $overlap->{$keybait}{$gene}{"strand"}{$strand}=1;

			   }

			   else{

			       $overlap->{$keybait}{$gene}={"tss"=>{$pos=>1}, "strand"=>{$strand=>1}};

			   }

		       }

		       else{

			   $overlap->{$keybait}={$gene=>{"tss"=>{$pos=>1}, "strand"=>{$strand=>1}}};

		       }

		   }

		    

		    $j++;

		}

	    }

	}

    }

}





##############################################################



sub printHelp{



    my $parnames=$_[0];

    my $parvalues=$_[1];

    

    print "\n";

    print "This script annotates promoter capture HiC baits with respect to genes. \n";

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

$parameters{"pathPromoterCaptureHiC"}="NA";

$parameters{"pathTranscriptInfo"}="NA";

$parameters{"pathGeneNames"}="NA";

$parameters{"maxDistance"}="NA";

$parameters{"pathOutput"}="NA";



my %defaultvalues;

my @defaultpars=("pathPromoterCaptureHiC", "pathTranscriptInfo", "pathGeneNames", "maxDistance", "pathOutput");



my %numericpars;

my @numericpars=();



foreach my $par (@numericpars){

    $numericpars{$par}=1;

}



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



print "Reading gene names...\n";



my %genenames;

readGeneNames($parameters{"pathGeneNames"}, \%genenames);

my $nbg=keys %genenames;



print "Found ".$nbg." genes in the annotations.\n";



print "Done.\n";



#####################################################################



print "Reading TSS positions...\n";



my %tss;

readTSS($parameters{"pathTranscriptInfo"}, \%tss);



print "Done.\n";



#####################################################################



print "Reading promoter capture HiC baits...\n";



my %acceptablechromo;



for(my $i=1; $i<=50; $i++){

    $acceptablechromo{"chr".$i}=1;

}



$acceptablechromo{"chrX"}=1;

$acceptablechromo{"chrY"}=1;

$acceptablechromo{"chrZ"}=1;

$acceptablechromo{"chrW"}=1;



my %baits;

readBaits($parameters{"pathPromoterCaptureHiC"}, \%acceptablechromo, \%baits);

 

print "Done.\n";



#####################################################################



print "Computing overlap...\n";



my $maxdist=$parameters{"maxDistance"};

print "maximum distance: ".$maxdist."\n";



my %overlap;

overlapCoords(\%baits, \%tss, $maxdist, \%overlap);

 

print "Done.\n";



#####################################################################



print "Writing output...\n";



my $pathin=$parameters{"pathPromoterCaptureHiC"};

my $input;

my @s=split("\\.",$pathin);

my $ext=$s[-1];



if($ext eq "gz"){

    open($input, "zcat $pathin |");

}

else{

    open($input, $pathin);

}



open(my $output, ">".$parameters{"pathOutput"});

my $line=<$input>;

chomp $line;

my @s=split("\t", $line);

my $len=@s;

my $lineout=join("\t", @s[0..2])."\tSymbol\tEnsembl Gene ID\t".join("\t",@s[3..($len-1)])."\n";

print $output $lineout;



$line=<$input>;



while($line){

    chomp $line;

    my @s=split("\t", $line);

    my $chr=$s[0];

    

    if(exists $acceptablechromo{$chr}){

	$chr=substr $chr, 3;

	my $start=$s[1]+0;

	my $end=$s[2]+0;

	my $key=$chr.",".$start.",".$end;



	if(exists $overlap{$key}){

	    my %names;

	    my $name="NA";

	    

	    foreach my $gene (keys %{$overlap{$key}}){

		if(exists $genenames{$gene}){

		    $names{$genenames{$gene}}=1;

		}

	    }

	    

	    my $idgene=join("|", keys  %{$overlap{$key}});

	    

	    my $lenname=keys %names;



	    if($lenname>0){

		$name=join("|", keys %names);

	    }

	    

	    my $lineout=join("\t", @s[0..2])."\t".$name."\t".$idgene."\t".join("\t",@s[3..($len-1)])."\n";

	    

	    print $output $lineout;

	}

	else{

	    print $output join("\t", @s[0..2])."\tNA\tnot_promoter\t".join("\t",@s[3..($len-1)])."\n";

	}

    }

    else{

	print $output join("\t", @s[0..2])."\tNA\tnot_promoter\t".join("\t",@s[3..($len-1)])."\n";

    }



    $line=<$input>;

}





close($output);

close($input);



print "Done.\n";



#####################################################################

