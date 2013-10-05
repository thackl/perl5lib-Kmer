#!/usr/bin/env perl
use warnings;
use strict;

use Benchmark qw(:all);
use Data::Dumper;

use Kmer;

my $seq_n = 1000;
my $seq_l = 100;
my @c = qw(A T G C);

my @const_len = map{rand_seq($seq_l)}1..$seq_n;
my @var_len = map{rand_seq( int(rand($seq_l)+$seq_l/2))}1..$seq_n;

my $kc = Kmer->new(kmer_size => 19);
my $k5 = Kmer->new(kmer_size => 19, shift_by => 5);


cmpthese(-1, {
	con_c => sub{
		$kc->kmerize($_) for @const_len;
	},
	con_c5 => sub{
		$k5->kmerize($_) for @const_len;
	},
	var_c => sub{
		$kc->kmerize($_) for @var_len;
	}
});



















sub rand_seq{
	my ($i, $seq) = (shift, '');	
	$seq .= $c[int(rand(4))] while $i--;
	return $seq;
}




