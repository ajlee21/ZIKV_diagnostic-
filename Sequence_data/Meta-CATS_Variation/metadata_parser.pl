############################
# Written by: Brett Pickett at JCVI (May 2016)
# Input 1: multiple sequence alignment in fasta format and 
# Input 2 TSV file with ID (matching fasta header) in first column and each metadata category in additional columns
# Output 1: meta-CATS Goodness of Fit test results for each metadata category
# Output 2: meta-CATS Test of Independence results for each metadata category
# Output 3: ID-metadata-sequence input file required for meta-CATS (for each metadata category)
# Output 4: Summary table showing number of significant sequence positions for each metadata category and for the combination of all categories
############################
#!/usr/bin/perl
use strict;
use warnings;
use Statistics::R;

my $metaDataFile = 'RPRC_Human_RSV_MasterMetadata_v2.txt';
my $seqFile = 'RPRC_completeGenome_Batch1-5_nogap-Concatenated.pep.fasta.afa'; #aligned sequences

my %categories;
my $lineNum = 0;
my $numCategories = 0;
my @headers;
my @data;
my @fileNames;

#load metadata into hash
open (METADATA, "<$metaDataFile") || die "$metaDataFile: $!\n";  #check to see if file exists
print "Loading Metadata...\n";

while (my $line = <METADATA>) {     #read all lines of the file
	chomp $line;
	if($lineNum ==1){ # read metadata lines and automatically populate 2D-hash`
		@data=split /\t/,$line;
		#print "$_\n" foreach(@data);
		#my $dataLength = @data;
		#print $dataLength;
		for (my $i = 1; $i < $numCategories; $i++){
			$categories{$data[0]}{$headers[$i]} = $data[$i]; #1D = seqID, 2D = metadata category, value = attribute
		}
	}
	elsif ($lineNum == 0){ #read header line and count number of metadata categories
		chomp $line;
		$lineNum = 1;
		@headers=split /\t/,$line;
		$numCategories = @headers;
		#print "$_\n" foreach(@headers);
#		for (my $h = 0; $h < $numCategories; $h++){
#			print "$h = $headers[$h]\n";
#		}
#		print "$numCategories\n";
	}
}
close (METADATA);
# #verify contents of 2D-hash for metadata
# foreach my $outer (keys %categories) {
# 	for (my $j = 1; $j < $numCategories; $j++){
# 		print "$outer, $headers[$j] = ", $categories{$outer}{$headers[$j]}, "\n";
# 	}
# }

# read in sequence data and store in same hash as metadata
open (SEQS, "<$seqFile") || die "$seqFile: $!\n"; 	#check to see if file exists
my $tempName;
print "Loading Sequence Data...\n";

while (my $line1 = <SEQS>) {		#read all lines of the file
	chomp $line1;
	if ($line1 =~m/^>(.+)/){
		$tempName = $1;
		#print "$tempName\n";
	}
	else {
		$categories{$tempName}{Seq} .= $line1;
	}
}
close (SEQS);
# #validate sequence key in hash
# foreach my $outer (keys %categories) {
# 	for (my $i = 1; $i < $numCategories; $i++){
# 		print "$outer, $headers[$i] = ", $categories{$outer}{$headers[$i]}, "\n";
# 	}
# 	print "$outer, $categories{$outer}{Seq}\n";
# }

#add 1 space between each sequence letter to format for downstream analysis
foreach my $outer (keys %categories) {
	$categories{$outer}{Seq} = join(" ",split(//,$categories{$outer}{Seq}));
	my $space = ' ';
	$categories{$outer}{Seq} = $space.$categories{$outer}{Seq};
}

# calculate number of different fields for each category
my %categoryCounts;
for (my $k = 1; $k < $numCategories; $k++){
	#print "$k\n";
	foreach my $outer (keys %categories) {
		$categoryCounts{$headers[$k]}{$categories{$outer}{$headers[$k]}}++;
		#print "$headers[$k],$categories{$outer}{$headers[$k]}:$categoryCounts{$headers[$k]}{$categories{$outer}{$headers[$k]}}\n";
	}
}
#validate proper loading of %categoryCounts
# foreach my $outer1 (keys %categoryCounts){ #outer = header name, inner = attribute field, value = counts
# 	foreach my $inner1 (keys $categoryCounts{$outer1}) {
# 		print "$outer1,$inner1: $categoryCounts{$outer1}{$inner1}\n";
# 	}
# }

#assign group numbers in new hash
my %assignments;
my @sortedOuter = sort keys %categoryCounts;
my $outerLength = @sortedOuter;
my @sortedInner;
my $innerLength;
my $number = 1;
#Assign values in %assignments using alphabetically-organized header and attribute names
foreach my $outer2 (@sortedOuter){
	@sortedInner = sort keys $categoryCounts{$outer2};
	#print @sortedInner, "\n";
	$innerLength = @sortedInner;
	$number = 1;
	foreach my $inner2 (@sortedInner) {
		$assignments{$outer2}{$inner2}= $number;
		#***#print "$outer2,$inner2:$assignments{$outer2}{$inner2}\n";
		$number++;
	}
}
#print "SortedOuter: @sortedOuter\n";

# construct separate output files for each category and sequence
foreach my $category1 (@sortedOuter){#metadata category
	#print "\n$category1\n";
	push(@fileNames, "$category1-rMsaInput.txt");
	open (OUTFILE, ">$category1-rMsaInput.txt") || die "$category1-rMsaInput.txt: $!\n";
	my @seqLabels = sort keys %categories;
	foreach my $seq1 (@seqLabels){
		my $attribute = $categories{$seq1}{$category1}; #points to numerical assignment representing metadata attribute
		#***#print "$seq1, $category1, $attribute, $assignments{$category1}{$attribute}\n";
		print OUTFILE "$seq1\t$assignments{$category1}{$attribute}\t$categories{$seq1}{Seq}\n";
	}
	close (OUTFILE);
}

#Create a communication bridge with R and start R to run meta-CATS code (in R)
my @results;
my $R = Statistics::R->new();
print "Running Statistics...\n";
#send @sortedOuter & @fileNames
foreach my $file1 (@fileNames){#metadata category
	if ($file1 =~m/(.+)-rMsaInput.txt/){
		my $name1 = $1;
		print "\t$name1...\n";
	}
	$R->set('inFilename', $file1);
	my $Rrun1 = $R->run_from_file("chisq_test.R");
	my $returnedFilename = $R->get('outfilename1');
	push(@results, "$returnedFilename");
}
#$R->stop();
#print "Results: @results\n";

#Parse all results files and save summaries of those that have highest number of significant positions
print "Parsing Results...\n";
open (SUMMARY, ">AA_Summary_Results-$metaDataFile.txt") || die "AA_Summary_Results-$metaDataFile.txt: $!\n";
print SUMMARY "Category Counts\nMetadata Category\t# Significant Positions\t# Total Positions\t% for each Category (# Sig / # Total )\n";
my %resultsStats;
my @resultsArray;
my $resultsLength = @results;
my $numSigResults =0;
for (my $p = 0; $p < $resultsLength; $p++){#metadata category
	#print "$file2\n";
	open (RESULTSDATA, "<$results[$p]") || die "$results[$p]: $!\n";
	#**#print "\t$results[$p]...\n";
	while (my $line4 = <RESULTSDATA>) {     #read all lines of the file
		#print "$line4\n";
		if ($line4=~m/^ .+/){
			@resultsArray=split /\t/,$line4;
			$resultsArray[0] =~ s/^\s+|\s+$//g;
			$resultsArray[2] =~ s/^\s+|\s+$//g;
			if($resultsArray[2] <= 0.05){
				#print "$resultsArray[0]\t$resultsArray[2]\n";
				$resultsStats{$results[$p]}{sig}++;
			}
		}
		else{
			next;
		}
		$resultsStats{$results[$p]}{total}++;
	}
	my $quotient = 100*($resultsStats{$results[$p]}{sig}/$resultsStats{$results[$p]}{total});
	$numSigResults += $resultsStats{$results[$p]}{sig};
	#my $quotient1 = printf("%.2f", $quotient);
	print SUMMARY "$sortedOuter[$p]\t$resultsStats{$results[$p]}{sig}\t$resultsStats{$results[$p]}{total}\t";
	printf SUMMARY ("%.1f",$quotient);
	print SUMMARY "%\n"; 
	close (RESULTSDATA);
}
my @percValues;
#add all significant sites for each category and calculate percentages for each
print SUMMARY "\nComprehensive Counts:\nMetadata Category\t# Significant Positions in Single Category\t# Significant Positions from ALL Categories\tPercentage for ALL Categories\n";
for (my $q = 0; $q < $resultsLength; $q++){#metadata category
	print SUMMARY "$sortedOuter[$q]\t$resultsStats{$results[$q]}{sig}\t$numSigResults\t";
	$resultsStats{$results[$q]}{totalPerc} = 100*($resultsStats{$results[$q]}{sig}/$numSigResults);
	#print "$resultsStats{$results[$q]}{totalPerc}\n";
	$resultsStats{$results[$q]}{totalPerc} = sprintf "%.1f", $resultsStats{$results[$q]}{totalPerc};
	print SUMMARY "$resultsStats{$results[$q]}{totalPerc} %\n";
	push(@percValues, "$resultsStats{$results[$q]}{totalPerc}");
}
close (SUMMARY);
$R->set('lbls', \@sortedOuter);
$R->set('slices', \@percValues);
$R->set('inFilename3', $metaDataFile);
my $Rrun2 = $R->run_from_file("pieChart.R");
$R->stop();

print "Complete\n";
print "Number of significant results for each metadata category can be found in:\n\"AA_Summary_Results-$metaDataFile\"\n\n";
print "Pie chart of results can be found in:\n\"AA-$metaDataFile-PieChart.png\"\n";