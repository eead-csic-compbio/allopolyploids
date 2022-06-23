package polyconfig;

# B Contreras-Moreira, R Sancho, EEAD-CSIC & EPS-UNIZAR 2018-20

use strict;
require Exporter;

our @ISA = qw( Exporter );

our @EXPORT = qw( 
  get_label_from_rules
  @diploids @polyploids @subgenomes 
  %outgroups $ROOT %sister_clades @CODES $NODEORDER
  $MINBLOCKLENGTH $MAXGAPSPERBLOCK $MINBLOCKOVERLAP
);

# Abbreviated names of diploid species as found in FASTA and tree files.
# See %sister_clades below for how to define sister species.
our @diploids = qw( Osat Hvul Bsta Bdis Barb Bpin Bsyl ); 

# note these are diploids as well; hash instead of list
our %outgroups = ( 
  'Osat',1, 
  'Hvul',1 
);

# diploid used to root trees 
our $ROOT = 'Osat';

# Abbreviated names of polyploid species as found in FASTA and tree files.
our @polyploids = ('Bmex', 'Bboi', 'Bret', 'Bhyb', 'Brup', 'Bpho', 'B422');

# Abbreviated names of polyploid subgenomes, defined by user.
# These are used to annotate alleles from the subgenomes as ancestral.
# Leave empty if not required
our @subgenomes = ();


# Abbreviated names of labelled polyploid sequences as found in FASTA and tree files.
# Note that the capital letters correpond to @CODES below
our @polyploids_labelled = (
   'Bmex_A','Bmex_B','Bmex_C','Bmex_D','Bmex_E','Bmex_F','Bmex_G','Bmex_H','Bmex_I',
   'Bboi_A','Bboi_B','Bboi_C','Bboi_D','Bboi_E','Bboi_F','Bboi_G','Bboi_H','Bboi_I',
   'Bret_A','Bret_B','Bret_C','Bret_D','Bret_E','Bret_F','Bret_G','Bret_H','Bret_I',
   'Bhyb_A','Bhyb_B','Bhyb_C','Bhyb_D','Bhyb_E','Bhyb_F','Bhyb_G','Bhyb_H','Bhyb_I',
   'Brup_A','Brup_B','Brup_C','Brup_D','Brup_E','Brup_F','Brup_G','Brup_H','Brup_I',
   'Bpho_A','Bpho_B','Bpho_C','Bpho_D','Bpho_E','Bpho_F','Bpho_G','Bpho_H','Bpho_I',
   'B422_A','B422_B','B422_C','B422_D','B422_E','B422_F','B422_G','B422_H','B422_I'
);

# Optional custom definition of clades that contain >1 diploids, if any.
# Leave empty or comment out the examples otherwise.
# MRCA nodes are internal node names, can be used in rules below; there are two keys:
# # i) MRCA for all species in all sister clades, also called a bifurcation in the code
# # ii) MRCA for each explicitely defined clade (there should be two)
# # finally the arrays contain lists of species in each clade, should be diploid
our %sister_clades = ();
$sister_clades{'MRCABsylBpin'}{'MRCA1'} = ['Bsyl','Bpin']; # clade 1
$sister_clades{'MRCABsylBpin'}{'MRCA2'} = ['Bsyl','Bpin']; # clade1 again with a different name, hack

# Default values for paremeters controlling how aligned blocks of sequences are
# produced in script _trim_MSA_block.pl
our $MINBLOCKLENGTH = 100;
our $MAXGAPSPERBLOCK = 100; # tolerated gaps for diploids in block
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
    if($anc_dip_taxon eq 'Bsta'){ $lineage_code = 'B' }
    elsif($anc_dip_taxon eq 'Bdis'){ $lineage_code = 'D' }
    elsif($anc_dip_taxon eq 'Barb'){ $lineage_code = 'F' }
    elsif($anc_dip_taxon eq 'Bsyl'){ $lineage_code = 'H' }
    elsif($anc_dip_taxon eq 'Bpin'){ $lineage_code = 'I' }
  }
  else { ## both ancestor/descendant diploids/clades are considered

    if($anc_dip_taxon eq 'Hvul' && $desc_dip_taxon eq 'Bsta'){ 
      $lineage_code = 'A' 
    }
    elsif($anc_dip_taxon eq 'Bsta' && $desc_dip_taxon eq 'Bdis'){ 
      $lineage_code = 'C' 
    }
    elsif($anc_dip_taxon eq 'Bdis' && $desc_dip_taxon eq 'Barb'){ 
      $lineage_code = 'E' 
    }
    elsif($anc_dip_taxon eq 'Barb' && $desc_dip_taxon eq 'MRCABsylBpin'){
      $lineage_code = 'G'
    }
    elsif($anc_dip_taxon eq 'MRCABsylBpin' && $desc_dip_taxon eq 'Bsyl'){
      $lineage_code = 'H'
    }
    elsif($anc_dip_taxon eq 'MRCABsylBpin' && $desc_dip_taxon eq 'Bpin'){
      $lineage_code = 'I'
    }
    elsif($anc_dip_taxon eq 'Bsyl' && $desc_dip_taxon eq 'Bpin'){
      $lineage_code = 'I'
    }
    elsif($anc_dip_taxon eq 'Bpin' && $desc_dip_taxon eq 'Bsyl'){
      $lineage_code = 'H'
    }
    # swapped Barb and Bsyl/Bpin
    elsif($anc_dip_taxon eq 'MRCABsylBpin' || $anc_dip_taxon eq 'Bpin'
          || $anc_dip_taxon eq 'Bsyl'){
      if($desc_dip_taxon eq 'Barb'){
        $lineage_code = 'F'
      }
    }
    elsif($anc_dip_taxon eq 'Bdis'){
      if($desc_dip_taxon eq 'MRCABsylBpin' || $desc_dip_taxon eq 'Bpin' 
         || $desc_dip_taxon eq 'Bsyl'){
        $lineage_code = 'E'
      }
    }
  }

  return $lineage_code;
}

1;
