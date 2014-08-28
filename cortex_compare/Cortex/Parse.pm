package Cortex::Parse;

use strict;
use warnings FATAL => 'all';
use FileHandle;
use Data::Dumper;
use Cortex::Match;
use Bio::Seq;
use Carp;

=head1 NAME

Cortex::Parse 

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';


=head1 SYNOPSIS

Parser for Cortex Match output

    use Cortex::Parse;

    my $foo = Cortex::Parse->new();
    ...

=head1 SUBROUTINES/METHODS

=head2 new

Create a new file object representing the Bubbleparse output

	my $bp = Cortex::Parse->new(

			-fasta => "somefile.fasta",
			
	);

=cut

sub new {
    my $class_name = shift;
    my $self = {};
    bless ($self, $class_name);
	my %arg = @_;
	
	#now do the fasta file
	$$self{_fastafile} = $arg{-fasta};
	##a hash of positions of the sequence header in the fasta files ..
	($$self{_match_file_index}, $$self{_fasta_fh}) = $self->_get_file_positions($$self{_fastafile});
	my @t =  keys %{$$self{_match_file_index}};
	#print Dumper @t;
	$$self{_matches} = \@t;
	#print Dumper $$self{_matches};
	return $self;
}

#get the position in bytes of each fasta entry (for coverage and match files) and make an index
sub _get_file_positions{
	
	my ($self,$file) = @_;
	my $info = {};
	my $fileh = FileHandle->new( $file, "r") || croak "could not open file $file";
	while (my $line = $fileh->getline){
		if ($line =~ m/^>/){
			my $string_length = _length_in_bytes($line);
			$line =~ m/^>match_(\d+)_path_(\d+)/;
			my ($match_num,$match_path) = ($1,$2);
			my $pos = tell($fileh) - $string_length;
			$$info{$match_num}{$match_path} = $pos; 
		}
	}
	return ($info, $fileh);
	#$$self{_match_file_index} = $info;
	#warn Dumper $info;
	#$$self{_fasta_fh} = $file;
}

#returns a hash of attributes for the match sought including sequence
sub _seek_sequence{
	my ($self,$match_num) = @_;
	my $seqs =  {};
	my $headers = {};
		
	foreach my $path (keys %{$$self{_match_file_index}{$match_num}}){
		my $pos = $$self{_match_file_index}{$match_num}{$path};
		seek($$self{_fasta_fh},$pos,0);
		my $header = $$self{_fasta_fh}->getline;
		chomp $header;
		$header =~ s/^>//;
		
		my @info = split(/\s+/,$header);
		my $name = shift @info;
		foreach my $pair (@info){
			my @pr = split(/:/,$pair);
			$$headers{$path}{$pr[0]} = $pr[1];
		}
		
		my $seq = "";
		while(my $line = $$self{_fasta_fh}->getline){
			last if $line =~ m/^>/;
			chomp $line;
			$seq .= $line;
		}
		#warn $seq;
		#warn $name;
		$$seqs{$path} = Bio::Seq->new(-seq => $seq, -id => $name);
	}
	
	
	return ($seqs, $headers);
}


sub _length_in_bytes{

	use bytes;
	return length shift;

}




=head2 next

Returns the next Bubble::Bubble object representing a single bubble

	while (my $bub = $bp->next){
	 ##do stuff
	}

=cut

sub next{
	my $self = shift;
	my $match = shift @$self{_matches};
	return 0 unless $match;
	my ($seq_obj_hash, $headers_hash) = $self->_seek_sequence($match);
	
	#work out the names of the attrs in the fasta header
	#and how many paths
	my %tmp;
	my %paths;
	foreach my $path (keys %{$headers_hash}){
		$paths{$path} = 1;
		foreach my $att (keys %{$$headers_hash{$path} }){
			$tmp{$att} = 1; 
		}
	} 
	
	my @seq_headers = keys %tmp;
	my @paths = keys %paths;
	return Cortex::Match->new(
		-seq_obj => $seq_obj_hash, #will be hashref of Bio::Seq objects (one for each path)
		-seq_header_info => $headers_hash, #will be hashref of info from fasta headers, (one set for each path)
		-seq_headers => \@seq_headers, #arrayref of header info in fasta file
		-paths => \@paths,
		-name => $match
	);
}



1; # End of Cortex::Parse
