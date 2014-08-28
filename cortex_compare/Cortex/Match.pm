package Cortex::Match;

use strict;
use warnings FATAL => 'all';
use Data::Dumper;


our $VERSION = '0.01';




=head2 new

Create a new object representing the Cortex::Match output. usually called from Cortex::Parse->next

	my $b = Cortex::Match->new(
			-seq_obj => $seq_obj_hash #Hashref for hash of path numbers and Bio::Seq objects e.g  (1 => Bio::Seq object)
			-seq_headers => $headers_hash #array ref of headers info in the fasta file
			-seq_header_info => #hashref for hash of path numbers and header info from fasta gule
			-paths => \@paths #array ref of path names through the bubble

	);

=cut

sub new {
    my $class_name = shift;
    my $self = {};
	my %args = @_;
    bless ($self, $class_name);

	$$self{_seq_obj} = $args{-seq_obj};
	$$self{_seq_headers} = $args{-seq_headers};
	$$self{_seq_header_info} = $args{-seq_header_info};
	$$self{_paths} = $args{-paths};
	$$self{_coverages} = $args{-coverages};
	$$self{_name} = $args{-name};
	return $self;
}


sub paths{
	my $self = shift;
	$$self{_paths};
}

sub name{
	my $self = shift;
	$$self{_name};
}

sub variant_identity {
	my ($self, $path) = @_;
	return substr($self->seq($path)->seq, $self->get('pre_length',$path), 1)
}

sub get {
	my ($self,$thing,$path) = @_;
	if(grep /$thing/i, @{$$self{_seq_headers} } ){
		return $$self{_seq_header_info}{$path}{$thing};
	}
	else{
		return undef;
	}
}

=head2 seq

returns Bio::Seq object of path sequence

	$b->seq('1') #returns Bio::Seq object of path 1
	$b->seq('1')->seq returns sequence string of path 1

=cut

sub seq{
	my ($self,$path) = @_;	
	return $$self{_seq_obj}{$path};
}


1; # End of Cortex::Match