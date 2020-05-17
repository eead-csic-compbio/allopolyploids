#!/usr/bin/perl 

# This script takes a FASTA multiple alignment of [nucleotide] sequences of
# diploid and polyploid species and produces a trimmed MSA suitable for 
# phylogenetic tree inference.
#
# The goal is to define a solid diploid backbone, which should be covered by
# outgroup sequences as well, and then use it to filter out polyploid 
# sequences with diploid block overlap < $min_block_overlap
#
# Uses diploid, outgroup and polyploid species defined in polyconfig.pm
# Also takes default block values, such as $min_block_overlap, from polyconfig.pm
#
# Only the longest sequence is taken for diploid species. 
#
# Note that outgroups are not used to compute diploid block.
#
# B Contreras-Moreira, R Sancho, EEAD-CSIC & EPS-UNIZAR 2018-20

use warnings;
use Getopt::Std;
use FindBin '$Bin';
use lib "$Bin";
use polyconfig;

my $min_block_length = $polyconfig::MINBLOCKLENGTH;
my $max_gaps_block   = $polyconfig::MAXGAPSPERBLOCK;
my $min_block_overlap= $polyconfig::MINBLOCKOVERLAP; 
my ($input_MSA_file,$block_MSA_file,$verbose,%opts) = ('','',0);
my ($remove_short_polyploids,$rename_file,$regex) = ('','','');

getopts('hi:o:O:m:M:s:R:t:v', \%opts);

if(($opts{'h'})||(scalar(keys(%opts))==0))
{
	print "\nusage: $0 [options]\n\n";
	print "-h this message\n";
	print "-i input FASTA file with aligned DNA sequences\n";
	print "-o output FASTA file with aligned block of sequences\n";   
	print "-v verbose output, useful for debugging              (optional)\n";
	print "-m min diploid block length                          ".
				"(optional, default:$min_block_length)\n";      
	print "-M max total gaps in block                           ".
				"(optional, default:$max_gaps_block)\n";
	print "-O fraction of block covered by outgroup/polyploids  ".
				"(optional, default:$min_block_overlap)\n";
	print "-R remove short polyploid isoforms using regex       ".
				"(optional, example: -R 'c\\d+_g\\d+_i')\n";
	print "-t TSV file with pairs long species name, abbrev     ".
				"(optional, example line: Oryza sativa\tOsat)\n";
	exit(0);
} 

if(defined($opts{'i'}) && -s $opts{'i'}){ 
	$input_MSA_file = $opts{'i'}; 
}
else{ die "# EXIT : need a valid input file\n" }

if(defined($opts{'o'})){
   $block_MSA_file = $opts{'o'};
}
else{ die "# EXIT : need a valid output file\n" }

if(defined($opts{'m'})){ $min_block_length = $opts{'m'} }

if(defined($opts{'M'})){ $max_gaps_block = $opts{'M'} }

if(defined($opts{'R'})){ 
	$remove_short_polyploids = 1;
	$regex = $opts{'R'};
}

if(defined($opts{'t'}) && -s $opts{'t'}){ $rename_file = $opts{'t'} }

if(defined($opts{'v'})){ $verbose = 1 }

# print parsed arguments
print "# $0 -i $input_MSA_file -o $block_MSA_file -m $min_block_length ".
			"-M $max_gaps_block -O $min_block_overlap -R $remove_short_polyploids ".
			"-t $rename_file -v $verbose\n\n";

print "# diploids: ".join(',',@polyconfig::diploids)."\n";
print "# outgroups: ".join(',',(keys(%polyconfig::outgroups)))."\n";
print "# polyploids: ".join(',',@polyconfig::polyploids)."\n\n";

################################################################################

my $MSAwidth = 0;
my $total_seqs_in_block = 0;
my $blockMSA = '';
my (%FASTA,%length,%inblock,%gaps,%header2taxon,%isoform,%long2short);
my ($ngaps,$bases,$header,$seq,$previous_seq,$taxon,$short_taxon);
my ($coord,$coord3,$coord5,$overlap,$is_redundant);
my %to_be_removed = (); #species names to be removed, used while debugging only

# read TSV file if required
if($rename_file){
	open(TSV,"<",$rename_file) || die "# ERROR: cannot read $rename_file\n";
	while(<TSV>){
		next if(/^#/);
		chomp;
		my @data = split(/\t/,$_);
		if($data[1]){
			$long2short{ $data[0] } = $data[1];
			print "# $data[0] -> $data[1]\n" if($verbose);
		}
	}
	close(TSV); 
	printf("# read %d pairs from file %s\n\n",scalar(keys(%long2short)),$rename_file);
}

# count how many diploids there are, excluding outgroups
# these should be included in block
my $diploids_minus_outgroups = 0;
foreach $taxon (@polyconfig::diploids){
	next if(defined($polyconfig::outgroups{$taxon}));
   $diploids_minus_outgroups++;
}


# Simplify headers and exclude sequences to be removed 
open(MSA,"<",$input_MSA_file) || 
	die "# ERROR: cannot open input file $input_MSA_file, exit\n";
while(my $line = <MSA>)
{
	# 1st regex=get_homologes 2nd regex=abbreviated sp names in place
	# note this regex might need to be edited with your own data
   if($line =~ /^>(.*?)(\[\S+\])$/ || $line =~ /^>(.*?)_([^_\n]+)$/)
   {
		$header = $1;
		$taxon = $2;
			
		if(!$to_be_removed{ $taxon })
		{ 
			$MSAwidth = 0;

			# take short taxon from TSV file if available
			$short_taxon = $long2short{$taxon} || $taxon;

			$header .= "_$short_taxon"; 
			$header2taxon{ $header } = $short_taxon; 
		}
   }
   else
	{	
		if($short_taxon ne '')
		{	
			chomp($line);
			$FASTA{ $short_taxon }{ $header } .= $line;
			
			# record sequence length
         $bases = ($line =~ tr/[ACGTN]//);
			$length{ $header } += $bases;
			$MSAwidth += length($line);
		}
	}
}
close(MSA);  

# define diploid block by checking the overlap of the longest sequences of taxons 
my $diploid_block_ok = 0;
my $diploid5prime = 0;
my $diploid3prime = $MSAwidth-1;
my $diploid_midcoord = 0;
my $diploid_block_length = 0;
foreach $taxon (@polyconfig::diploids)
{ 
	#skip outgroups
	next if(defined($polyconfig::outgroups{$taxon}));

	foreach $header (sort {$length{$b}<=>$length{$a}} keys(%{ $FASTA{$taxon} }))
	{	
		print "$taxon $length{$header} $header\n" if($verbose);
		
		# to know which diploids have been considered for this calculation
		$diploid_block_ok++;

		# record 5' side of block of aligned diploids
      $coord = $diploid5prime;
      while(substr($FASTA{$taxon}{$header},$coord,1) eq '-'){ $coord++ }
      if($coord > $diploid5prime){ $diploid5prime = $coord }

      # record 3' side of block of aligned diploids
      $coord = $diploid3prime;
      while(substr($FASTA{$taxon}{$header},$coord,1) eq '-'){ $coord-- }
      if($coord < $diploid3prime){ $diploid3prime = $coord }
      print "# after $taxon : $diploid5prime to $diploid3prime\n" if($verbose);

      $diploid_block_length = $diploid3prime - $diploid5prime + 1;

		$inblock{$header} = $total_seqs_in_block++;
		
		last; #consider only longest sequence of $short_taxon
	} 				
}

# compute block length
$diploid_block_length = $diploid3prime - $diploid5prime + 1;

if($diploid_block_ok < $diploids_minus_outgroups)
{
   print "# MSA skipped as not all required diploids are aligned: $diploid_block_ok\n";
   exit;
}
elsif($diploid_block_length < $min_block_length)
{
	print "# MSA skipped due to short diploid block: $diploid_block_length\n";
   exit;
}

# compute midpoint of diploid block
$diploid_midcoord = $diploid3prime - (($diploid3prime - $diploid5prime)/2);

print "# aligned diploid block: $diploid5prime < $diploid_midcoord > $diploid3prime\n\n";

# check gap fraction of sequences in diploid block
foreach $header (sort {$inblock{$a}<=>$inblock{$b}} keys(%inblock))
{
	$seq = substr( $FASTA{ $header2taxon{$header} } { $header },
							$diploid5prime,
                     $diploid3prime - $diploid5prime + 1 );

	$ngaps = ($seq =~ tr/\-//);

	if($ngaps > $max_gaps_block)
	{
		print "# MSA skipped due to gappy diploid block: $ngaps ($header)\n";
		exit;
	}
}

# check overlap of outgroup sequences
foreach $taxon (keys(%polyconfig::outgroups))
{
   foreach $header (sort {$length{$b}<=>$length{$a}} keys(%{ $FASTA{$taxon} }))
   {
      print "outg $taxon $length{$header} $header\n" if($verbose);

		# make sure this sequence overlaps with diploid block
      # 5'
      $coord = $diploid5prime;
      while(substr($FASTA{$taxon}{$header},$coord,1) eq '-'){ $coord++ }
      $coord5 = $coord;
      # 3'
      $coord = $diploid3prime;
      while(substr($FASTA{$taxon}{$header},$coord,1) eq '-'){ $coord-- }
      $coord3 = $coord;

		$overlap = ($coord3 - $coord5 + 1) / $diploid_block_length;
	
		if( $overlap >= $min_block_overlap ){
			$inblock{ $header } = $total_seqs_in_block++;			
		}
		else
		{
			print "# MSA skipped as outgroup $taxon has poor block coverage: $overlap\n";
			exit;
		}

		last;
	}
}

# identify polyploid isoforms if requested
if($remove_short_polyploids)
{
	my ($isof);
	foreach $taxon (@polyconfig::polyploids)
	{
		my %tx_isoforms;
		foreach $header (sort {$length{$b}<=>$length{$a}} keys(%{ $FASTA{$taxon} }))
		{
			if($header =~ m/($regex)/)
			{
				$isoform = $1;
				$tx_isoforms{ $isoform }++;
				if($tx_isoforms{ $isoform } > 1)
				{
					$isoform{$taxon}{$header} = 1;
				}
			}
		}
	}
}

## check overlap of polyploid sequences
foreach $taxon (@polyconfig::polyploids)
{
	my (@seqs_in_block); # only for this poly species

   foreach $header (sort {$length{$b}<=>$length{$a}} keys(%{ $FASTA{$taxon} }))
   {
      print "poly $taxon $length{$header} $header\n" if($verbose);
		
		if($isoform{$taxon}{$header})
		{
			print "# polyploid sequence $header skipped for being a redundant Trinity isoform\n";
			next;
		}

		# make sure this sequence overlaps with diploid block
      # 5'
      $coord = $diploid5prime;
      while(substr($FASTA{$taxon}{$header},$coord,1) eq '-'){ $coord++ }
      $coord5 = $coord;
      # 3'
      $coord = $diploid3prime;;
      while(substr($FASTA{$taxon}{$header},$coord,1) eq '-'){ $coord-- }
      $coord3 = $coord;

		$overlap = ($coord3 - $coord5 + 1) / $diploid_block_length;

		if( $overlap >= $min_block_overlap )
		{
			print "$header $overlap\n" if($verbose);

			# 1st seq of this poly species to be considered
			if(!@seqs_in_block)
			{
				$inblock{ $header } = $total_seqs_in_block++;
				push(@seqs_in_block, $header); # header of this sequence
			}
			else # make sure 2nd, 3rd seqs are not redundant
			{
				$is_redundant = 0;	
				$seq = substr( $FASTA{ $header2taxon{$header} } { $header },
                        $coord5,
                        $coord3 - $coord5 + 1 );				
				
				foreach my $previous_header (@seqs_in_block)
				{
					$previous_seq = substr( $FASTA{ $header2taxon{$header} } { $previous_header },
													$coord5,
													$coord3 - $coord5 + 1 );

					print "$seq $header\n$previous_seq $previous_header\n\n" if($verbose);

					if($previous_seq eq $seq)
					{
						print "# polyploid redundant sequence skipped: $header\n";			
						$is_redundant = 1;
						last;
					}
				}	

				if($is_redundant == 0)
				{
					$inblock{ $header } = $total_seqs_in_block++;
					push(@seqs_in_block, $header);
				}
			}
      }
      else
      {
         print "# polyploid sequence $header skipped due to poor overlap: $overlap\n";
      }
	}
}


# create MSA of sequences trimmed to cover diploid block
open(BLOCKMSA,">",$block_MSA_file) || die "# ERROR: cannot create $block_MSA_file\n";

foreach $header (sort {$inblock{$a}<=>$inblock{$b}} keys(%inblock))
{
	print BLOCKMSA ">$header\n";
	print BLOCKMSA substr( $FASTA{ $header2taxon{$header} } { $header }, 
							$diploid5prime,
							$diploid3prime - $diploid5prime + 1 ) ."\n";		
}

close(BLOCKMSA);

print "\n# outfile: $block_MSA_file\n";
