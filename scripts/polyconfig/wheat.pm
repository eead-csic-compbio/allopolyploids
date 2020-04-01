package wheat;

# B Contreras-Moreira, R Sancho, EEAD-CSIC & EPS-UNIZAR 2018-20

use strict;
require Exporter;

our @ISA = qw( Exporter );

our @EXPORT = qw(
	get_label_from_rules
	@diploids @polyploids %outgroups $ROOT %sister_clades @CODES $NODEORDER
	$MINBLOCKLENGTH $MAXGAPSPERBLOCK $MINBLOCKOVERLAP
);

# Abbreviated names of diploid species as found in FASTA and tree files.
# See %sister_clades below for how to define sister species.
our @diploids = qw( Osat Hvul Bdis Tura Tmon Asha Atau Aspe ); 

# note these are diploids as well; hash instead of list
our %outgroups = ( 
'Osat',1, 
'Hvul',1,
'Bdis',1 
);

# diploid used to root trees 
our $ROOT = 'Osat';

# Abbreviated names of polyploid species as found in FASTA and tree files.
our @polyploids = ('Taes', 'Ttur');

# Abbreviated names of labelled polyploid sequences as found in FASTA and tree files.
# Note that the capital letters correpond to @CODES below
our @polyploids_labelled = (
   'Taes_A','Taes_B','Taes_C','Taes_D','Taes_E','Taes_F','Taes_G','Taes_H','Taes_I',
   'Ttur_A','Ttur_B','Ttur_C','Ttur_D','Ttur_E','Ttur_F','Ttur_G','Ttur_H','Ttur_I'
);

# Optional custom definition of clades that contain >1 diploids, if any.
# Leave empty or comment out the examples otherwise.
# MRCA nodes are internal node names, can be used in rules below; there are two keys:
# # i) MRCA for all species in all sister clades, also called a bifurcation in the code
# # ii) MRCA for each explicitely defined clade (there should be two)
# # finally the arrays contain lists of species in each clade, should be diploid
our %sister_clades = ();
$sister_clades{'MRCAAeTr'}{'MRCAAe'} = ['Asha','Atau','Aspe']; # clade 1
$sister_clades{'MRCAAeTr'}{'MRCATr'} = ['Tura','Tmon'];        # clade 2

# Default values for paremeters controlling how aligned blocks of sequences are
# produced in script _trim_MSA_block.pl
our $MINBLOCKLENGTH = 200;
our $MAXGAPSPERBLOCK = 200; # tolerated gaps for diploids in block
our $MINBLOCKOVERLAP = 0.50; # fraction of diploid block covered by outgroups & polyploids

# Optional user-defined contribution of diploid species to block width calculations
# in script _trim_MSA_block.pl (outgroup species are not used)
# Can be used to indicate that only a member of a clade is required, see Bsyl example 
#my %diploids4width = (
# Example wheat values
#   'Tura'=>1,
#	'Tmon'=>1,
#	'Asha'=>1,
#	'Atau'=>1,
#	'Aspe'=>1,
# Brachypodium values
#'Bsta'=>1,
#'Bdis'=>1,
#'Barb'=>1,
#'Bpin'=>1,
#'Bsyl'=>0.2,
#);


# see rules defined below
our @CODES = qw( A B C D E F G H I all ); 

# use to ladderize trees
our $NODEORDER = 0; # 1:increasing, 0:decreasing

# Takes 4 parameters: 
# 1) string with abbreviated name of diploid/MRCA taxon (ancestor)
# 2) string with abbreviated name of diploid/MRCA taxon (descendant)
# 3) boolean scalar to indicate ancestor is sister
# 4) boolean scalar to indicate descendat is sister
# Returns a label, which is either a code from @CODES or '-' otherwise
sub get_label_from_rules {

	my ($anc_dip_taxon, $desc_dip_taxon, $anc_is_sister, $desc_is_sister) = @_;
	
	my  $lineage_code = '-';

	# check input params
	my $ancOK = 0;
	if(grep(/^$anc_dip_taxon/,@diploids)){ $ancOK = 1 }
	if(defined($sister_clades{$anc_dip_taxon})){ $ancOK = 1 }
	else {
		foreach my $MRCA (keys(%sister_clades)){
			if(defined($sister_clades{$MRCA}{$anc_dip_taxon})){
				$ancOK = 1; 
				last;
			}
		}
	}	
	if($ancOK == 0){
		print "# ERROR (get_label_from_rules): unrecognized ancestor $anc_dip_taxon\n";
		return $lineage_code;
	}	

	my $descOK = 0;
   if(grep(/^$desc_dip_taxon/,@diploids)){ $descOK = 1 }
   if(defined($sister_clades{$desc_dip_taxon})){ $descOK = 1 }
   else {
      foreach my $MRCA (keys(%sister_clades)){
         if(defined($sister_clades{$MRCA}{$desc_dip_taxon})){
            $descOK = 1;
            last;
         }
		}
   }
   if($descOK == 0){
      print "# ERROR (get_label_from_rules): unrecognized descendant $desc_dip_taxon\n";
      return $lineage_code;
   }

	# start applying rules
	
	## ancestor is sister or descendant is empty, 
	## only ancestor diploid is looked up
	if($anc_is_sister == 1 || $desc_dip_taxon eq ''){
		if($anc_dip_taxon eq 'Aspe'){ $lineage_code = 'C' }
		elsif($anc_dip_taxon eq 'Asha'){ $lineage_code = 'E' }
		elsif($anc_dip_taxon eq 'Atau'){ $lineage_code = 'F' }
		elsif($anc_dip_taxon eq 'Tmon'){ $lineage_code = 'H' }
		elsif($anc_dip_taxon eq 'Tura'){ $lineage_code = 'I' }
	}
	else { ## both ancestor/descendant diploids/clades are considered

		if($anc_dip_taxon eq 'Hvul' && $desc_dip_taxon eq 'MRCAAeTr'){ 
			$lineage_code = 'A' 
		}
		elsif($anc_dip_taxon eq 'MRCAAeTr' && $desc_dip_taxon eq 'MRCAAe'){ 
			$lineage_code = 'B' 
		}
		elsif($anc_dip_taxon eq 'MRCAAeTr' && $desc_dip_taxon eq 'MRCATr'){ 
			$lineage_code = 'G' 
		}
		elsif($anc_dip_taxon eq 'MRCAAe' || $anc_dip_taxon eq 'Aspe'){
			if($desc_dip_taxon eq 'Asha' || $desc_dip_taxon eq 'Atau'){ 
				$lineage_code = 'D' 
			}
		}
      elsif($anc_dip_taxon eq 'MRCATr' && $desc_dip_taxon eq 'Tura'){
            $lineage_code = 'I'
		}
      elsif($anc_dip_taxon eq 'MRCATr' && $desc_dip_taxon eq 'Tmon'){
            $lineage_code = 'H'
      }
      elsif($anc_dip_taxon eq 'Atau' && $desc_dip_taxon eq 'Asha'){
         $lineage_code = 'E'
      }
		elsif($anc_dip_taxon eq 'Asha' && $desc_dip_taxon eq 'Atau'){ 
			$lineage_code = 'F' 
		}
      elsif($anc_dip_taxon eq 'Tura' && $desc_dip_taxon eq 'Tmon'){
         $lineage_code = 'H'
      }
      elsif($anc_dip_taxon eq 'Tmon' && $desc_dip_taxon eq 'Tura'){
         $lineage_code = 'I'
      }

	}

	return $lineage_code;
}

1;
