#!/usr/bin/perl -w

# General purpose pattern searcher, that produces output that can be used in mathematica directly.

use strict;
use 5.010;

my $sUsage = "Usage: perl $0 <input fasta filename> <patterns>";
if (@ARGV <= 1) {die $sUsage;}
my ($input_contig_file) = $ARGV[0];
my @readseq_arr;	# Array containing read sequences
my $next_seqID;
my $fullseq;
my $motnum	= length(@ARGV);
my @patterns  = @ARGV[1, $motnum+1];
my @patt_copy = @patterns;
my %match_num;

my $n = 0;

# modifying patterns to expand ambiguity codes and taking reverse complementarity into account

while ($patterns[$n]) {
	my $exp_pattern = replace_ambiguous($patterns[$n]);
	my $rev_pattern = revdnacomp($exp_pattern);
    $patterns[$n] = $exp_pattern."|".$rev_pattern;
	$n++;
}

open(FASTA, $input_contig_file) || die "Cannot open input file.\n";

# Processing align file
open(FASTA, $ARGV[0]);

my $seqID = <FASTA>;
die "Input file does not appear to be in Fasta format" if ($seqID !~ /^>/);
$seqID =~ s/^>//;
chomp($seqID);

print "patternout={";
while (my $seq = <FASTA>) {
	chomp($seq);
	until ($seq =~ /^>/ || eof(FASTA)) {
		$fullseq .= $seq;
		$seq = <FASTA>;
		chomp($seq);
	}
	if ($seq =~ /^>/) {
		$seq =~ s/^>//;	# Removing Fasta character
		chomp($seq);
		$next_seqID = $seq;
	}
	if (eof(FASTA)) {
		chomp($seq);
		$fullseq .= $seq;
	}

	# extract range
	my @id = split(/\s/, "$seqID");
	my $found_any = 0;
	$n = 0;
	my $fullout = "";

	while ($patterns[$n]) {
		my $tot	      = 0;
		my $repcount	 = 0;
		my $curr_pattern = $patterns[$n];
		my $outline	  = "";
		while ($fullseq =~ m/$curr_pattern/ig) {
			$outline .= "$-[0], ";
			$repcount++;
			$match_num{$patt_copy[$n]}++;
		}
		$outline =~ s/, $//;
		$found_any += $repcount;
		if ($repcount) {$fullout .= "{{\"$patt_copy[$n]\"}, {$outline}}, ";}
		$n++;
	}
	$fullout =~ s/, $//;
	print "{\"$id[0]\",$fullout},\n" if ($found_any);

	#	   print "$fullseq\n";
	$seqID   = $next_seqID;
	$fullseq = "";
}
close(FASTA);
foreach my $pattern (keys %match_num) {
	print "\"$pattern\" $match_num{$pattern}\n";
}

sub revdnacomp {
	my $dna	 = shift;
	$dna =~ tr/ACGTacgtN/TGCAtgca./;
	my @elem = $dna =~ /((?:\[.*?\]|.)(?:\{.*?})?)/g;
	my $rev = join '', reverse @elem;
	return $rev;
}

sub replace_ambiguous {
	my $patt = shift;
	$patt =~ s/R/\[AG\]/ig;
	$patt =~ s/Y/\[CT\]/ig;
	$patt =~ s/S/\[GC\]/ig;
	$patt =~ s/W/\[AT\]/ig;
	$patt =~ s/K/\[GT\]/ig;
	$patt =~ s/M/\[AC\]/ig;
	$patt =~ s/B/\[CGT\]/ig;
	$patt =~ s/D/\[AGT\]/ig;
	$patt =~ s/H/\[ACT\]/ig;
	$patt =~ s/V/\[ACG\]/ig;
	return $patt;
}
