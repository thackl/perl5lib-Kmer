#!/usr/bin/env perl
use warnings;
no warnings 'qw';
use strict;

use Carp;
use Getopt::Long;
use Pod::Usage;
use FindBin;
use File::Basename;
use Log::Log4perl qw(:easy :no_extra_logdie_message);
use Data::Dumper;
$Data::Dumper::Sortkeys = 1;

use Fastq::Parser;
use Fastq::Seq;
use Jellyfish;
use Kmer;
use Verbose::ProgressBar;


our $VERSION = '0.02';


##------------------------------------------------------------------------##

=head1 NAME 

read_kmer_filter.pl

=cut

=head1 DESCRIPTION

Filter reads based on kmer count information.

=cut


=head1 CHANGELOG

=cut

=head2 0.02

=over

=item [Change] Only retrieve _2 counts if _1 did not fail.

=item [Feature] --kmer-shift - only take every n-th kmer into account

=back

=head2 0.01

=over

=item [Initial]

=back

=cut

=head1 TODO

--help

-q 0 <= q <= 100

#=item --both [OFF]
#
#Only in paired mode. By default, reads pair is discarded if either one read 
# fails the filter. With --both pairs, where both reads fail are discarded.
#
#=item --dump-discarded
#
#Dump discarded reads in extra files.
#
#=item -o|--out <FILE>
#
#Output file name. There will be thread wise temporary output files created
# at --out location and finally merged to --out.
#
#=item [-t|--threads <INT>] [1]
#
#Number of threads to use

=cut

=head1 SYNOPSIS

  read_kmer_filter.pl [<OPTIONS>] -h <HASH> -k <KMERSIZE> -1 <FASTQ> [ <FASTQ> ... ] [-2 <FASTQ> <FASTQ> ... ]
  
=cut

=head1 OPTIONS

=over

=item -h|--kmer-hash <FILE>

Jellyfish kmer count hash to use for count assignment.

=item -1|--reads <FASTQ> [<FASTQ> ...]

Read files in FASTQ format, single or first of pair.

=item -2|-mates <FASTQ> [<FASTQ> ...]

Read files in FASTQ format, second of pair. Use the same order as for --reads.

=item -k|--kmer-size

Size of kmers used in Jellyfish hash.

=item -u|--upper <INT> [0]

Upper kmer count cutoff. 0 deactivates cutoff.

=item -l|--lower <INT> [0]

Lower kmer count cutoff. 0 deactivates cutoff.

=item -q|--quantile/--quantile-lower/--quantile-upper <INT> [50 (median)]

Quantile of kmers in reads above/below cutoff to trigger filter. Use -q to
 set for both, -u and -l, --quantile-lower/--quantile-upper to set options
 individually.

=item [--debug]

Turn on debug messages.

=item [--quiet]

=item [--help]

Show this help screen.

Supress verbose information messages.

=back

=cut

Log::Log4perl->init(\<<'CFG');
	log4perl.logger.main				= DEBUG, Screen
	log4perl.appender.Screen			= Log::Log4perl::Appender::Screen
	log4perl.appender.Screen.stderr		= 0
	log4perl.appender.Screen.layout		= PatternLayout
	log4perl.appender.Screen.layout.ConversionPattern = [%d{MM-dd HH:mm:ss}] [rkf] %m%n

CFG

my $L = Log::Log4perl::get_logger();
$L->level($INFO);

my @fq_suffixes=(qw/.fq .fastq .FQ .FASTQ/);
my @opt_reads;
my @opt_mates;
my %opt = (
	reads => \@opt_reads,
	mates => \@opt_mates,
	lower => \(my $opt_l = 0),
	upper => \(my $opt_u = 0),
	quantile => \(my $opt_q = 50),
	'quantile-lower' => \(my $opt_qlf = undef),
	'quantile-upper' => \(my $opt_quf = undef),
	'kmer-shift' => 1,
);

Getopt::Long::Configure("no_ignore_case");
GetOptions(\%opt, qw(
	reads|1=s@{,}
	mates|2=s@{,}
	kmer-hash|h=s
	kmer-size|k=i
	kmer-shift|x=i
	lower|l=i
	upper|u=i
	quantile|q=i
	quantile-lower=i
	quantile-upper=i
	quiet
	debug
	help
)) or $L->logcroak($!);

pod2usage(1) if $opt{help};

$opt{quiet} && $L->level($WARN);
$opt{debug} && $L->level($DEBUG);

$L->debug("GetOptions:\n", Dumper(\%opt));

##------------------------------------------------------------------------##	
# required	
for(qw(reads kmer-hash kmer-size)){
	pod2usage("required: --$_") unless defined ($opt{$_}) 
};

if(@opt_mates){
	$L->logcroak("Number of -1 and -2 files differs!") if @opt_mates != @opt_reads;
}

# check db files -e -s
foreach my $file(@opt_reads, @opt_mates){
	$L->logcroak("Cannot find file: $file ") unless -e $file && -f $file;
}

# cutoff and quantiles
$opt_qlf //= $opt_q; #//
$opt_quf //= $opt_q; #//
$opt_qlf /= 100;
$opt_quf /= 100;

$L->debug("factor quantile lower: $opt_qlf");
$L->debug("factor quantile upper: $opt_quf");


$L->warn("--lower($opt_l) > --upper($opt_u), are you sure that's what you want?") if $opt_u && $opt_l > $opt_u;

##------------------------------------------------------------------------##	

my $jf = Jellyfish->new();
my $km = Kmer->new(
	kmer_size => $opt{'kmer-size'},
	shift_by => $opt{'kmer-shift'}
);

my $FC;
unless(@opt_mates){
	$L->info("Mode: single end");
	
	for($FC=0; $FC < @opt_reads;$FC++){
		
		$L->info("File: $opt_reads[$FC]");
		my $fp1 = Fastq::Parser->new(file => $opt_reads[$FC]);
		open (FQ1, '>', basename($opt_reads[$FC], @fq_suffixes).".fil.fq") or $L->logcroak("$!");

		my $pgc=0;
		my $pg = Verbose::ProgressBar->new(
			size => $fp1->fh,
			level => 2,
			report_level => $opt{quiet} ? 0 : 2,
		);

		while(my $fq1 = $fp1->next_seq){
			$pg->update unless $pgc++%10000;
			
			my @counts = sort{$a<=>$b} $jf->query(
	    		['--both-strands', $opt{'kmer-hash'}], 
				kmers => [$km->cmerize($fq1->seq)],
				table=>0 
			);
			
			if(
				$opt_u && $counts[int($#counts * $opt_quf)] > $opt_u
				||
				$counts[int($#counts * $opt_qlf)] < $opt_l
			){ 
				# TODO: dumping of discarded ...
			}else{
				print FQ1 "$fq1";
			}
		}
		$pg->finish;
		close FQ1;
	}	
}else{
	$L->info("Mode: paired end");
	
	for($FC=0; $FC < @opt_reads;$FC++){
		$L->info("Files: $opt_reads[$FC] $opt_mates[$FC]");
		my $fp1 = Fastq::Parser->new(file => $opt_reads[$FC]);
		open (FQ1, '>', basename($opt_reads[$FC], @fq_suffixes).".fil.fq") or $L->logcroak("$!");
		my $fp2 = Fastq::Parser->new(file => $opt_mates[$FC]);
		open (FQ2, '>', basename($opt_mates[$FC], @fq_suffixes).".fil.fq") or $L->logcroak("$!");

		my $pgc=0;
		my $pg = Verbose::ProgressBar->new(
			size => $fp1->fh,
			level => 2,
			report_level => $opt{quiet} ? 0 : 2
		);
		
		while(
			(my $fq1 = $fp1->next_seq) &&
			(my $fq2 = $fp2->next_seq)
		){
			$pg->update unless $pgc++%10000;
			my @counts1 = sort{$a<=>$b} $jf->query(
	    		['--both-strands', $opt{'kmer-hash'}], 
				kmers => [$km->cmerize($fq1->seq)],
				table=>0 
			);
			
			if(
				$opt_u && $counts1[int($#counts1 * $opt_quf)] > $opt_u
				||
				$counts1[int($#counts1 * $opt_qlf)] < $opt_l
			){ 
				$L->debug("discarded by _1: ".$fq1->id);
			}else{
				# only retrieve _2 counts if _1 did not fail
				my @counts2 = sort{$a<=>$b} $jf->query(
		    		['--both-strands', $opt{'kmer-hash'}], 
					kmers => [$km->cmerize($fq2->seq)],
					table=>0 
				);
				if(
					$opt_u && $counts2[int($#counts2 * $opt_quf)] > $opt_u
					||
					$counts2[int($#counts2 * $opt_qlf)] < $opt_l
				){
					$L->debug("discarded by _2: ".$fq1->id);
				}else{
					print FQ1 "$fq1";
					print FQ2 "$fq2";
				}
			}
		}
		$pg->finish;
		close FQ1;
		close FQ2;
	}
}
















