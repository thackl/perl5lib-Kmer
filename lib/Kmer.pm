package Kmer;

use warnings;
use strict;

use Verbose;

our $VERSION = '0.01';


##------------------------------------------------------------------------##

=head1 NAME 

Kmer.pm

=cut

=head1 DESCRIPTION

Extract kmers from sequence strings.

=cut

=head1 SYNOPSIS

=cut

=head1 CHANGELOG

=cut

=head2 0.01

=over

=item [Refacture] OO handle and OO METHODS instead of Class METHODS

=item [PODfix]

=item [Initial]

=back

=cut

=head1 TODO

=over

=back

=cut


##------------------------------------------------------------------------##

=head1 Class Attributes

=cut

##------------------------------------------------------------------------##

=head1 Class METHODS

=cut

##------------------------------------------------------------------------##

=head1 Constructor METHOD

=head2 new

  $kh = Kmer->new(kmer_size => 19);

=cut

sub new{
	my $proto = shift;
	my $self;
	my $class;
	
	# object method -> clone + overwrite
	if($class = ref $proto){ 
		return bless ({%$proto, @_}, $class);
	}

	# class method -> construct + overwrite
	# init empty obj
	$self = {
		kmer_size => undef,
		@_
	};
	
	die "kmer_size required" unless $self->{kmer_size};
	
	$self->{_unpack} = "(A".$self->{kmer_size}."X".($self->{kmer_size}-1).")"; 
	
	return bless $self, $proto;
}



##------------------------------------------------------------------------##

=head1 Object METHODS

=cut

=head2

  $kh->kmer_size()
    # 19
  $kh->kmer_size(4)
    # 4

=cut

sub kmer_size{
	my ($self, $ks) = @_;
	if($ks){
		$self->{kmer_size} = $ks;
		$self->{_unpack} = "(A".$ks."X".($ks-1).")";
	};
	return $self->{kmer_size};
}

=head2 kmerize

Factor a STRING into a LIST of overlapping kmers. Kmers are returned in 
 their literal version.

  $kh->kmerize("ATAGG");
    # ATAG,TAGG

=cut

sub kmerize{
	my ($self,$seq) = @_;
	return unpack($self->{_unpack}.((length $seq) - $self->{kmer_size}+1 ), $seq);
}

=head2 cmerize

Factor a STRING into a LIST of overlapping kmers. Kmers are returned in 
 their canonical representation. The canonical representation is the
 lexically first of the literal and the reverse complement version of any
 kmer.

  $kh->cmerize("ATAGG");
    # ATAG,CCTA

=cut

sub cmerize{
	my ($self,$seq) = @_;
	map{my $krc = reverse $_; $krc =~ tr/ATGC/TACG/; $_ gt $krc ? $krc : $_}unpack($self->{_unpack}.((length $seq) - $self->{kmer_size}+1 ), $seq);
}


##------------------------------------------------------------------------##

=head1 Accessor METHODS

=cut

=head1 AUTHOR

Thomas Hackl S<thomas.hackl@uni-wuerzburg.de>

=cut



1;



