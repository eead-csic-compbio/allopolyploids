#!/usr/bin/perl -w
use strict;

# Script that reads trimmed, aligned clusters of sequences,
# shortens taxon names and puts outgroup sequences (1/per outgroup) ahead
# 
# B Contreras-Moreira, R Sancho, EEAD-CSIC & EPS-UNIZAR 2017

die "# usage: $0 <trimmed filed>" if(!$ARGV[0]);

my ($taxon,$short_taxon);
my %long2short = qw(
[B422_80_75] _B422
[Bdistachyon_314_v3.1] _Bdis
[Osativa_323_v7.0.transcript] _Osat
[Sbicolor_313_v3.1.transcript] _Sbic
[Hvulgare_IBSC2016.HQ.cds] _Hvul
[arb_80_75] _Barb
[boi_80_75] _Bboi
[hyb_80_75] _Bhyb
[mex_80_75] _Bmex
[pho_80_75] _Bpho
[pin_80_75] _Bpin
[ret_80_75] _Bret
[rup_80_75] _Brup
[sta_80_75] _Bsta
[syl_Cor_80_75] _BsyC
[syl_Esp_80_75] _BsyE
[syl_Gre_80_75] _BsyG
);

my %outgroup = qw( _Osat 0 _Hvul 1 ); # _Sbic 0

# read trimmed alignment
my $trimmed_file = $ARGV[0];
my $outfile = $trimmed_file;
my $printOK = 0;
my %printed_outgroups;
my ($ingroups,@outgroups);
$outfile =~ s/fna.aln.fna/fixname.fna/;

open(RAW,$trimmed_file);
while(<RAW>)
{
	if(/^>.*?(\[\S+?\])$/)
	{
		$taxon = $1;
		$short_taxon = $long2short{$taxon};
		s/\Q$taxon\E/$short_taxon/;

		# init this sequence
		$printOK = 1;

		if(defined($outgroup{$short_taxon}))
		{
			if(!defined($printed_outgroups{$short_taxon})){
				$printed_outgroups{$short_taxon}++;
      }
      else{ $printOK = 0; next; }

      # add this header to vector in position $outgroup{$short_taxon}
      $outgroups[$outgroup{$short_taxon}] .= $_;
    }  
    else{
      $ingroups .= $_;
    }
	}
	else
  {
    if($printOK)
    {
      if(defined($outgroup{$short_taxon})){ $outgroups[$outgroup{$short_taxon}] .= $_ }
      else{ $ingroups .= $_ } 
    }
  }	
}
close(RAW);

# print cluster starting with sorted outgroups
open(GRAMPA,">",$outfile) || die "# cannot create $outfile\n";

print GRAMPA join('',@outgroups);
print GRAMPA $ingroups;

close(GRAMPA);
