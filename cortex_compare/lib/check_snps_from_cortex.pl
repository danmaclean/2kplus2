#!/usr/bin/perl
# encoding: utf-8

use strict;
use warnings;
use Carp;
use Cortex::Parse;
use Cortex::Match;
use Bio::SearchIO;
#
#  check_snps_from_cortex.pl
#
#  Created by Dan MacLean (TSL) on 2014-08-27.
#  Copyright (c). All rights reserved.
#

##check database files exist
my $BLAST_DATABASE = "db/reference_genome.fa";
croak "Need the Arabidopsis TAIR10 sequence as a blastn database in $BLAST_DATABASE\n" unless -e $BLAST_DATABASE;


my $c = Cortex::Parse->new( -fasta => $ARGV[0] );

my %lines;
while (my $bubble = $c->next){
	foreach my $path (@{$bubble->paths()} ){
		my $pre_length =  $bubble->get('pre_length',$path);
		my $var = $bubble->variant_identity($path);
		my $contig = $bubble->name($path);
		my %blast_info = %{blast_path($bubble->seq($path)->seq,$pre_length,$var)};
		foreach my $hit (keys %blast_info){
			next if $hit eq 'hits' or $hit eq 'no_hit';
			if ($blast_info{$hit}{'hsps'} eq 1  and $blast_info{$hit}{'is_snp_path'} eq "TRUE" and $blast_info{$hit}{'snp_position_in_reference'} ne 'NA'){
				my $result = join("\t", $hit, $blast_info{$hit}{'snp_position_in_reference'}, $contig) . "\n";
				$lines{$result} = 0;
			}
		}
	}
	
}
foreach my $line (keys %lines){
	print $line;
}

##############################################################################

sub blast_path{
	my $seq = uc(shift);
	my $snp_pos_in_path = uc(shift);
	my $snp_id = uc(shift);
	my $command = "echo  $seq | blastn -db $BLAST_DATABASE -query -";
	my $fh;
	open $fh, "$command |" || die("cannot run blast cmd of $command: $!\n");
	my $in = new Bio::SearchIO(-format => 'blast',  -fh  => $fh);
	my %summary;
	while (my $result = $in->next_result){
		$summary{'hits'} = $result->num_hits;
		if ($result->num_hits == 0){
			$summary{'no_hit'}{'snp_position_in_reference'} = 'NA';
			$summary{'no_hit'}{'is_snp_path'} = 'NA';
			$summary{'no_hit'}{'hsps'} = 0;
			$summary{'no_hit'}{'strand_query'} = "NA";
			$summary{'no_hit'}{'strand_hit'} = "NA";
		}
		else{
			while (my $hit = $result->next_hit){
				while (my $hsp = $hit->next_hsp){
					$summary{$hit->name}{'hsps'} = $hit->num_hsps;
					my ($is_snp_path, $snp_position) = get_snp_position($hit, $hsp, $snp_pos_in_path, $snp_id);
					$summary{$hit->name}{'snp_position_in_reference'} = $snp_position;
					$summary{$hit->name}{'is_snp_path'} = $is_snp_path;
					$summary{$hit->name}{'strand_query'} = $hsp->strand('query');
					$summary{$hit->name}{'strand_hit'} = $hsp->strand('hit');
				}
			}
		}	
	}
	return \%summary;
}

sub get_snp_position{
	my ($hit,$hsp,$pos_in_branch, $snp_id) = @_;
	my $pos = "NA";
	$pos = $hsp->start('hit') + $pos_in_branch - 1;
	if ($hsp->strand('hit') == -1 ){
		$pos = $hsp->end('hit') - $pos_in_branch -1;		
	}
	my $is_snp_path = 'FALSE';
	if (substr($hsp->hit_string, $pos_in_branch, 1) eq $snp_id){
		$is_snp_path = 'TRUE';
	}
	return ( $is_snp_path, $pos);
}

