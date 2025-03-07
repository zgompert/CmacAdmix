#!/usr/bin/perl
#
# then subset of files from a given set


foreach $in_t (@ARGV){
	
	#change number of loci below to match your own files:
	$nloc = 79079; 

	#example name to match: 
	# input_allele_freq_files/AS.BF.1.F1.C_varne_infile.txt
	
	#matching a part of the infile that contains the pop-rep-gen-host information: 
	#This is called $in_t because it's the input file for time t
	$in_t =~ m/input_allele_freq_files\/AS\.([A-Z_]+)\.([0-9]+)\.F([0-9]+)\.([A-Z]+)_varne_infile\.txt/;
	#pulling out the pop, rep, gen, and host from the input file name above:
	$pop = $1;
	$rep = $2;
	$gen = $3;
	$host = $4;

	#making a new object called $in_0 (input file for time point zero) that is equal to $in_t
	$in_0 = "input_allele_freq_files/AS.$pop.1.F1.C_varne_infile.txt";
	
	#making input file for the startpoint-midpoint and startpoint-endpoint Ne calculation (F1 to F7 or 20):
	#$in_0 = "input_allele_freq_files/AS.$pop.1.F1.C_varne_infile.txt";

	#making an input file for the midpoint-endpoint Ne calculation (F7 to F20):
	#$in_0 = "input_allele_freq_files/AS.$pop.$rep.F[567].${host}_varne_infile.txt";


	#substitute the end timepoint (i.e. F20 or F7) with the starting timepoint (i.e. F1)
	#This will create the appropriate initial allele freq file that belongs with the ending allele freqs 
	#$in_0 =~ s/20/1/ or die "failed sub $in_t\n";
	

	#running varne for the 19-gen timepoint (F1 to F20) with census size of 2000: 
	$lag = $gen +1;
	#$lag = $gen;
	#$lag = $gen -1;
	system "/uufs/chpc.utah.edu/common/home/u6000989/bin/varne -a $in_0 -b $in_t -l $nloc -t $lag -n 2000 -x 1000 -m 0.01 > output_Ne_files/Ne_AS.${pop}.${rep}.F1_to_F${gen}.${host}_census_2000.txt\n";

	#running varne for the 6-gen timepoint (F1 to F7): 
	#system "/uufs/chpc.utah.edu/common/home/u6000989/bin/varne -a $in_0 -b $in_t -l $nloc -t 6 -n 2000 -x 1000 > output_Ne_files/Ne_AS.${pop}.${rep}.F1_to_F${gen}.${host}_census_2000.txt\n";}

	##running varne for the mid-to-end timepoint (F7 to F20):
	#system "/uufs/chpc.utah.edu/common/home/u6000989/bin/varne -a $in_0 -b $in_t -l $nloc -t 8 -n 2000 -x 1000 > output_Ne_files/Ne_AS.${pop}.${rep}.F7ish_to_F${gen}.${host}_census_2000.txt\n";}
}
