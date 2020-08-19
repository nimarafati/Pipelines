#!/usr/bin/perl -W
# 150224 Changed to print first/second parent instead of ref/alt allele. Added extra field in the end with reference allele.
# 150226 Changed to check the correct column for ACGT (= skip rows with only gaps, but ACGT in quality column)
#
#use strict;
use Math::CDF qw(:all);

my ($vcf_file,$bam_file) = @ARGV;

#print "VCF: $vcf_file\n";
#print "BAM: $bam_file\n";

open(IN,$vcf_file) or die "Cannot open infile: $vcf_file\n";

while (my $row = <IN>) {
	next if $row =~ /^#/;

	my @line = split "\t",$row;		
	my ($chr,$coord) = ($line[0],$line[1]);
	my ($first,$second) = ($line[3],$line[4]);

	# Determine which parent has what genotype
#	my @par1 = split ":",$line[9];
#	my $ref = $first;

	# Check if first parent has alternative allele and swap
#	if ($par1[0] eq "1/1") {
#		($first,$second) = ($second,$first);
#	}

	my $results = `samtools mpileup -Q 20 $bam_file -r $chr:$coord-$coord 2> /dev/null`;
	my @res = split "\t",$results;
	if (defined($res[4]) && $res[4] =~ /[ACGTacgt]/) {	# Only consider rows with at least one base sequenced (=not only gaps)
		my $N = $res[3];
		my $reads = $res[4];	
		my @chars = split "",$reads;

		my %ACGT;
		my $i = 0;
		while ($i < @chars) {	# Read char by char
			if ($chars[$i] eq "\$") { $i++; next; }	# End marker, skip
			if ($chars[$i] eq "\^") { $i+=2; next; }	# Start marker, skip
			if ($chars[$i] eq "+") {	# Insertion
				my $size;
				my $pos = 1;
				while ($chars[$i+$pos] =~ /\d/) {
					$size .= $chars[$i+$pos];
					$pos++;
				}
				my $insert = join "",@chars[$i+$pos .. $i+$pos+$size-1];
				$ACGT{uc $chars[$i-1].uc $insert}++;
				$ACGT{uc $chars[$i-1]}--;
				$i += $pos+$size;
				next;
			}
			if ($chars[$i] eq "-") {   # Deletion
			my $size;
			my $pos = 1;
			while ($chars[$i+$pos] =~ /\d/) {
				#print "$i\t$pos\t$chars[$i+$pos]";<STDIN>;
				$size .= $chars[$i+$pos];
				$pos++;
			}
		        my $del = join "",@chars[$i+$pos .. $i+$pos+$size-1];#concatenate cahracters
#			print "$chars[$i-1]\t$ACGT{$chars[$i-1]}";<STDIN>;
		        $ACGT{uc $chars[$i-1]."_".uc $del}++;
			$ACGT{uc $chars[$i-1]}--;
			$i += $pos+$size;
			next;
			}
#			print "$chars[$i-1]\t$ACGT{$chars[$i-1]}";<STDIN>;
			$ACGT{uc $chars[$i]}++;	# Count all other characters
			$i++;	# Increase position counter
		}

		# Compute p-value (two-sided binomial)
#		my $trials = ($ACGT{$first}//0)+($ACGT{$second}//0);
#		my $success = min($ACGT{$first}//0,$ACGT{$second}//0);
#		my $p = 1;
#		if ($trials) {
#			$p = 2*pbinom($success,$trials,0.5);
#		}

#		my $fail = $N-$trials;	# Non-reference and non-alternative allele
		# Do not count aligned gaps
#		my $gaps = ($ACGT{"<"}//0) + ($ACGT{">"}//0);
#		$N -= $gaps;
#		$fail -= $gaps;

#		my $p_fail = pbinom($N-$fail,$N,0.99);

		$start=$coord-1;
		print "$chr\t$start\t$coord\t$first\t$second\t".($ACGT{$first}//0)."\t".($ACGT{$second}//0)."\t";#$p\t";
		foreach my $key (sort {$ACGT{$b} <=> $ACGT{$a}} keys %ACGT) {
			print "$key=$ACGT{$key},";
		}
		print "\n";
#		print "\t".($p_fail//1)."\t$ref\n";
	}
}

close IN;

sub min {
	return $_[0] if $_[1] > $_[0];
	return $_[1];
}
