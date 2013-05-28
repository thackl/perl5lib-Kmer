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

map{my $krc = reverse $_; $krc =~ tr/ATGC/TACG/; print $krc, "\n"}@kmers;


#--------------------------------------------------------------------------#
=head2 load module

=cut

BEGIN { use_ok('Kmer'); }

my $class = 'Kmer';

#--------------------------------------------------------------------------#


subtest "$class->KmerSize" => sub{
	can_ok($class, 'KmerSize');
	is(Kmer->KmerSize, 19, 'KmerSize default');
	is($Kmer::KmerSize, 19, '$KmerSize default');
	is(Kmer->KmerSize(4), 4, 'KmerSize set');
	is($Kmer::KmerSize, 4, '$KmerSize set');

};

subtest "$class->Kmerize" => sub{
	can_ok($class, 'Kmerize');
	is_deeply([$class->Kmerize($seq)], \@kmers, 'Kmerize');
};


subtest "$class->Kmerize_nr" => sub{
	can_ok($class, 'Kmerize_nr');
	is_deeply([$class->Kmerize_nr($seq)], \@kmers_nr, 'Kmerize nr');
};

done_testing();

__END__


