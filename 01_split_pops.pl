#!/usr/bin/perl
#
# splits a genotype likelihood file by population and year
# (or whatever you decide to put in the parentheses for matching)
#
# usage splitPops.pl in.gl
#
# in order to use this script, you should modify the matching line to match 1) the full name of the individual from the 2nd line of the .gl file (everything between the forward slashes of m/stuff/) and 2) modify the location of the parentheses in the matching line to match only the bits that correspond to your treatment group. 
#
# for example, my beetles are named something like CA.3.F9.C.5.G8, following the format population.replicate.generation.host.plate.well. But the treatment group I want to split by is each pop.rep.gen.host combo (i.e. I do not care about plate or well). So for my script, the parenthese for $1 should capture pop.rep.gen.host, but exclude plate and well.  

#takes a file passed to it in @ARGV and names it $in: 
my $in = shift (@ARGV);


open(IN, $in) or die "failed to open the infile\n";

## get number of loci and number of individuals from first line of .gl file: 
$line = <IN>;
chomp($line);
$line =~ m/(\d+)\s+(\d+)/;
$nind = $1;
$nloc = $2;

## get ind. and pop. ids (second line of the .gl file): 
$line = <IN>;
chomp($line);
@line = split (" ",$line);
foreach $ind (@line){

	#The matching line below should be modified
	#such that it matches your .gl file ID line names: 
	# example ID: CA.3.F9.C.5.G8	
	# I am matching a series of A-Z, 0-9, and dots
	# followed by a dot, some digits (\d+), a dot...
	# and finally some combination of A-Z and 0-9

	$ind =~ m/([A-Z0-9\.\_]+)\.\d+\.[A-Z0-9]+/;
	$id = $1;
	
	#print the id of the treatment group to the out file:
	print "$id\n";
	#save the treatment group name to a matrix/scaffold thing:
	push (@id,$id);
	push (@{$popids{$id}},$ind);
	$ids{$id} = 1;
	if(defined $popn{$id}){
		$popn{$id}++;
	}
	else {
		$popn{$id} = 1;
	}
}

## open one file per population
foreach $id (sort keys %ids){
	$fh = "F"."$id";
	$out = "$id"."_$in";
	print "$out\n";
	open ($fh, "> $out") or die "Could not write $id\n";
	$files{$id} = $fh;
	print {$files{$id}} "$popn{$id} $nloc\n";
	$pids = join (" ",@{$popids{$id}});
	print {$files{$id}} "$pids\n";
	@ones = ();
	for($i=0;$i<$popn{$id}; $i++){
		push (@ones,1);
	}
	$ones = join (" ", @ones);
	print {$files{$id}} "$ones\n";
}

## read and write
while (<IN>){
	chomp;
	@line = split (" ",$_);
	$a = shift(@line); ## locus info
	foreach $id (sort keys %ids){
		print {$files{$id}} "$a";
	}
	for ($i=0; $i<$nind; $i++){
		$id = $id[$i];
		for ($j=0; $j<3; $j++){
			$a = shift(@line); 
			print {$files{$id}} " $a";
		}
	}	
	foreach $id (sort keys %ids){
		print {$files{$id}} "\n";
	}
			
}
close (IN);
