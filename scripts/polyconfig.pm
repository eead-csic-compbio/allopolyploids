package polyconfig;

# B Contreras-Moreira, R Sancho, EEAD-CSIC & EPS-UNIZAR 2018-20

use strict;
require Exporter;

our @ISA = qw( Exporter );

our @EXPORT = qw(
	get_label_from_rules
	@diploids @polyploids %sister_clades @CODES $NODEORDER
);

# Abbreviated names of diploid species as found in FASTA and tree files.
# See %sister_clades below for how to define sister species.
our @diploids = qw( Osat Hvul Bdis Tura Tmon Asha Atau Aspe ); 

# Abbreviated names of polyploid species as found in FASTA and tree files.
our @polyploids = ('Taes', 'Ttur');

# Abbreviated names of labelled polyploid sequences as found in FASTA and tree files.
# See @CODES below for the labels.
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

# see rules defined below
our @CODES = qw( A B C D E F G H I all ); 

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
	
	## sister rules
   if($anc_dip_taxon eq 'Aspe' && $anc_is_sister == 1){$lineage_code = 'C'}
   if($anc_dip_taxon eq 'Asha' && $anc_is_sister == 1){$lineage_code = 'E'}
   if($anc_dip_taxon eq 'Atau' && $anc_is_sister == 1){$lineage_code = 'F'}
   if($anc_dip_taxon eq 'Tmon' && $anc_is_sister == 1){$lineage_code = 'H'}
   if($anc_dip_taxon eq 'Tura' && $anc_is_sister == 1){$lineage_code = 'I'}


	## ancestor/descendant diploids/clades rules
   if($anc_dip_taxon eq 'Hvul' && $desc_dip_taxon eq 'MRCAAeTr'){ $lineage_code = 'A' }

   if($anc_dip_taxon eq 'MRCAAeTr' && $desc_dip_taxon eq 'MRCAAe'){ $lineage_code = 'B' }

   if($anc_dip_taxon eq 'MRCAAeTr' && $desc_dip_taxon eq 'MRCATr'){ $lineage_code = 'G' }

   if($anc_dip_taxon eq 'Aspe'){ # or $anc_dip_taxon eq 'MRCAAe'
		if($desc_dip_taxon eq 'Asha' || $desc_dip_taxon eq 'Atau'){ $lineage_code = 'D' }
		}

   if($anc_dip_taxon eq 'Asha'){
		if($desc_dip_taxon eq 'Atau'){ $lineage_code = 'E' }
		#elsif
      }












	if($anc_dip_taxon eq 'Hvul'){
		if($desc_dip_taxon eq 'Aspe'){ $lineage_code = 'B' }
			elsif($desc_dip_taxon eq 'Tura' || $desc_dip_taxon eq 'Tmon'){ $lineage_code = 'G' }
         else{ $lineage_code = 'A' }
   }


      elsif($anc_dip_taxon eq 'Aspe'){
         if($desc_dip_taxon eq 'Asha' || $desc_dip_taxon eq 'Atau'){ $lineage_code = 'D' }
         elsif($desc_dip_taxon eq 'Tura' || $desc_dip_taxon eq 'Tmon'){ $lineage_code = 'G' }
         elsif($desc_dip_taxon eq 'Aspe' && $desc_is_sister == 1){ $lineage_code = 'C' }
			else{ $lineage_code = 'C' }

      }

		elsif($anc_dip_taxon eq 'Asha'){
         if($desc_dip_taxon eq 'Tura' || $desc_dip_taxon eq 'Tmon'){ $lineage_code = 'G' }
         elsif($desc_dip_taxon eq 'Asha' && $desc_is_sister == 1){ $lineage_code = 'E' }
         elsif($desc_dip_taxon eq 'Atau'){ $lineage_code = 'F' }
         else{ $lineage_code = 'E' }

      }

		elsif($anc_dip_taxon eq 'Atau'){
         if($desc_dip_taxon eq 'Tura' || $desc_dip_taxon eq 'Tmon'){ $lineage_code = 'G' }
         elsif($desc_dip_taxon eq 'Atau' && $desc_is_sister == 1){ $lineage_code = 'F' }
         elsif($desc_dip_taxon eq 'Asha'){ $lineage_code = 'E' }
         else{ $lineage_code = 'F' }
      }

      elsif($anc_dip_taxon eq 'Tmon'){
         if($desc_dip_taxon eq 'Tura'){ $lineage_code = 'I' }
			elsif($desc_dip_taxon eq 'Tmon' && $desc_is_sister == 1){ $lineage_code = 'H' }
			else{ $lineage_code = 'H' }
      }

		elsif($anc_dip_taxon eq 'Tura'){
         if($desc_dip_taxon eq 'Tmon'){ $lineage_code = 'H' }
			elsif($desc_dip_taxon eq 'Tura' && $desc_is_sister == 1){ $lineage_code = 'I' }
         elsif($desc_dip_taxon eq 'Aspe'){ $lineage_code = 'B' }
			else{ $lineage_code = 'I' }
		}

	return $lineage_code;
}

1;
