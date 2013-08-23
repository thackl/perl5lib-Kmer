#!/usr/bin/env perl

# $Id$

use strict;
use warnings;

use Test::More;
use Data::Dumper;

use FindBin qw($RealBin);
use lib "$RealBin/../lib/";


#--------------------------------------------------------------------------#
=head2 sample data

=cut

my $seq = "ATTATATATCGACTAGCC";
my @kmers = qw(
	ATTA
	TTAT
	TATA
	ATAT
	TATA
	ATAT
	TATC
	ATCG
	TCGA
	CGAC
	GACT
	ACTA
	CTAG
	TAGC
	AGCC
);

my @kmers_nr = (
	'ATTA',	# TAAT
	'ATAA',	# TTAT rc
	'TATA',	# TATA
	'ATAT',	# ATAT
	'TATA',	# TATA
	'ATAT',	# ATAT
	'GATA',	# TATC rc	
	'ATCG',	# CGAT
	'TCGA',	# TCGA
	'CGAC',	# GTCG
	'AGTC',	# GACT rc
	'ACTA',	# TAGT
	'CTAG',	# CTAG
	'GCTA',	# TAGC rc
	'AGCC',	# GGCT
);


#--------------------------------------------------------------------------#
=head2 load module

=cut

BEGIN { use_ok('Kmer'); }

my $Class = 'Kmer';

#--------------------------------------------------------------------------#

my $obj;
subtest 'new object' => sub{
	$obj = new_ok($Class, [kmer_size => 19]);
	is($obj->{kmer_size}, 19, "attribute kmer_size")
};

subtest '$obj->kmer_size' => sub{
	can_ok($obj, 'kmer_size');
	is($obj->kmer_size, 19, 'kmer_size get');
	is($obj->kmer_size(4), 4, 'kmer_size set');
};

subtest '$obj->kmerize' => sub{
	can_ok($obj, 'kmerize');
	is_deeply([$obj->kmerize($seq)], \@kmers, 'kmerize');
};

subtest '$obj->kmerize_nr' => sub{
	can_ok($obj, 'kmerize_nr');
	is_deeply([$obj->kmerize_nr($seq)], \@kmers_nr, 'kmerize_nr');
};

done_testing();




















