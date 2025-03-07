#!/usr/bin/perl
#
# convert sam to bam, then sort and index 
#


use Parallel::ForkManager;
my $max = 36;
my $pm = Parallel::ForkManager->new($max);

FILES:
foreach $sam (@ARGV){ #calling the sam files I passed from shell script
	$pm->start and next FILES; ## fork
		
	#example name to match is aln_BZ.5.F20.L.16.C3.sam, can also be underscores 		
	#matching the beetle id out of the sam file name for use in naming future files:
	$sam =~ m/^output_sam_files\/([A-Za-z0-9_\.]+)\.sam/ or die "failed to match $sam\n";
	$base = $1;

	#converting sam to bam: 
	system "samtools view -b -O bam $sam > output_bam_files/$base.bam\n"; 
	#-b option specifies that output should be in bam format... redundant??
	#-O option specifies the format you want the output in
	#That said, BCF tools will automatically select what it thinks the appropriate output format is based on your filename extension
        
	#sorting file to match the order of sequences in my reference genome:
	system "samtools sort -O bam output_bam_files/$base.bam > output_sorted_bam_files/sorted_$base.bam\n"; 
	
	#creates an index for fast look up of particular sequences within the bam files:
        system "samtools index -b output_sorted_bam_files/sorted_$base.bam\n"; 
	#-b option creates a BAI index. This is the default, so not using this option would yield the same result
        
	$pm->finish;
}

$pm->wait_all_children;



