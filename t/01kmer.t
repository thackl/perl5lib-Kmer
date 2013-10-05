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
	my @attr = (
		kmer_size => 19,
		shift_by => 1,
		_u_shift => -18,
		_u_patt => '(A19X18)%dA19',
		_u_size => 1,
	);
	while(@attr){
		my ($attr, $re) = (shift @attr, shift @attr);
		is($obj->{$attr}, $re, "attribute $attr");
	}
	is_deeply($obj->{_u_tpl},{},"attribute _u_tpl");
};

subtest '$obj->kmer_size' => sub{
	can_ok($obj, 'kmer_size');
	is($obj->kmer_size, 19, 'kmer_size get');
	is($obj->kmer_size(4), 4, 'kmer_size set');
};

subtest '$obj->shift_by' => sub{
	can_ok($obj, 'shift_by');
	is($obj->shift_by, 1, 'shift_by get');
	is($obj->shift_by(4), 4, 'shift_by set');
	is($obj->shift_by(1), 1, 'shift_by reset');
};

subtest '$obj->kmerize' => sub{
	can_ok($obj, 'kmerize');
	is_deeply([$obj->kmerize($seq)], \@kmers, 'kmerize');
	foreach my $incr (2..20){
		$obj->shift_by($incr);
		my $i=0;
		is_deeply([$obj->kmerize($seq)], [grep{! ($i++%$incr)}@kmers], 'kmerize shift_by => '.$incr);
	}
	$obj->shift_by(1);
};

subtest '$obj->cmerize' => sub{
	can_ok($obj, 'cmerize');
	is_deeply([$obj->cmerize($seq)], \@kmers_nr, 'cmerize');
	foreach my $incr (2..20){
		$obj->shift_by($incr);
		my $i=0;
		is_deeply([$obj->cmerize($seq)], [grep{! ($i++%$incr)}@kmers_nr], 'cmerize shift_by => '.$incr);
	}
	$obj->shift_by(1);
};

done_testing();




















