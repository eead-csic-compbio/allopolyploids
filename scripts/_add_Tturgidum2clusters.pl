#!/usr/bin/perl
use strict;
use warnings;

# Take precomputed wheat FASTA files (Marcussen), 
# add the best Triticum turgidum (Ttur) sequences A & B from Ensembl Plants cDNA and 
# build multiple sequence alignments.

# B Contreras-Moreira, R Sancho, EEAD-CSIC & EPS-UNIZAR 2020

my ($seqid,$sp,$name,$shortsp,%seqs);

# https://science.sciencemag.org/content/345/6194/1250092.full
my $fastadir  = 'Marcussen_et_al_2014/275genes_final_alignments/';
my $blastdir  = 'blast/Marcussen_et_al_2014/275genes_final_alignments/';

# ftp://ftp.ensemblgenomes.org/pub/release-46/plants/fasta/triticum_turgidum/cdna/
my $FASTATtu  = 'Triticum_turgidum.Svevo.v1.cdna.all.fa';

my $aligndir  = '01_align/';
my $clusterdir= 'clusters/';

my $MSAexe= 'mafft'; # add full path if necessary

# add 4-letter abbreviations
# see https://science.sciencemag.org/content/345/6194/1250092.full
my %name2sp = (
	'Ta'  => 'Taes',
	'Tm'  => 'Tmon',
	'Tu'  => 'Tura',
	'AE'  => 'Atau',
	'Ash' => 'Asha',
	'Asp' => 'Aspe',
	'Hv'  => 'Hvul',
	'Bra' => 'Bdis',
	'LOC' => 'Osat'
); 

# change to 1 if you want to add missing Ttur sequences padded as Ns
my $ADDMISSINGSEQS = 0; 

# read Ttu sequences
open(FASTA,"<",$FASTATtu) ||die "#cannot read $FASTATtu\n";
while(<FASTA>){
   next if(/^$/);
   chomp;
	if(/^>(\S+)/){ $seqid = $1 }
	else{ $seqs{$seqid} .= lc($_) }
}
close(FASTA);

opendir(BLASTDIR,$blastdir) || die "# cannot list $blastdir\n";
my @files = grep{/\.fa/} readdir(BLASTDIR);
closedir(BLASTDIR);

my ($qseqid,$sseqid,$pident,$len,$mism,$gap,$qstart,$qend,$sstart,$send,$eval,$bit);

foreach my $file (@files){

	#next if($file ne '1L_OrthoMCL_group5261_0.7_NoGaps.fa'); # debug

	# find out best Tturgidum sequences for this cluster by parsing BLASTN results
	my %stats;
	open(BLAST,"<",$blastdir.$file) || die "# cannot read $blastdir$file\n";
	while(<BLAST>){
		#Hv_MLOC_44488.1 TRITD7Av1G247000.2      95.522  67      3       0       1497    1563    1404  1470     1.08e-21        108
		chomp:
		($qseqid,$sseqid,$pident,$len,$mism,$gap,$qstart,$qend,$sstart,$send,$eval,$bit) = split(/\t/,$_);
		$stats{$sseqid}{'tot'}++;
		$stats{$sseqid}{'bit'}+=$bit;
	}
	close(BLAST);  

   # take only 1 A gene and 1 B
   my ($subg,$A,$B,%subgenomes) = ('','','');
   foreach $sseqid (sort {$stats{$b}{'tot'}<=>$stats{$a}{'tot'} ||
									$stats{$b}{'bit'}<=>$stats{$a}{'bit'} } keys(%stats)){
		
		$subg = '';
		if($sseqid =~ m/\d+([ABU])v/){ $subg = $1 }
		$subgenomes{$subg}++;
		next if( $subgenomes{$subg} > 1 );
		if($A eq '' && ($subg eq 'A' || $subg eq 'U')){ $A = $sseqid }
		if($B eq '' && ($subg eq 'B' || $subg eq 'U')){ $B = $sseqid }
		#print "$file $sseqid $subg $stats{$sseqid}{'tot'} $stats{$sseqid}{'bit'}\n";
	}

	# create cluster with Marcussen plus Ttu sequence(s)
	open(CLUSTER,">",$clusterdir.$file) || die "# cannot write to $clusterdir$file\n";

	open(FASTA,"<",$fastadir.$file) || die "# cannot read $fastadir$file\n"; 
	while(<FASTA>){
		if(/^>(\S+)/){
			$name = $1;
			$shortsp = '';
			foreach $sp (keys(%name2sp)){
				if($name =~ /^$sp/){
					$shortsp = $name2sp{$sp};
					last;
				}
			}
			print CLUSTER ">$1_$shortsp\n";
		} 
		else { print CLUSTER }
	}
	close(FASTA);

	if($A || $B){
		print CLUSTER ">$A\_Ttur\n$seqs{$A}\n" if($A ne '' && $seqs{$A});
		print CLUSTER ">$B\_Ttur\n$seqs{$B}\n" if($B ne '' && $seqs{$B});	
	}
	close(CLUSTER);
	
	# align seqs
	system("$MSAexe -i $clusterdir$file > $aligndir$file");

	# gauge alignment width
	my $width=0;
	open(ALIGN,"<",$aligndir.$file) || die "# cannot read $aligndir$file\n";
	while(<ALIGN>){
		next if(/^>/);
		chomp;
		$width += length($_);
	}
	close(ALIGN);
	
	# add missing seqs so that all alignments have same num of seqs 
	if(!$A || !$B){
		if(!$A){
			warn "# $aligndir$file no A gene\n";
			if($ADDMISSINGSEQS){
				$A ='TRITD_Ttur';
				$seqs{$A} = 'n' x $width;
				open(ALIGN,">>",$aligndir.$file) || die "# cannot write to $aligndir$file\n";
				print ALIGN ">$A\n$seqs{$A}\n";
				close(ALIGN);
			}
		}

		if(!$B){
         warn "# $aligndir$file no B gene\n";
			if($ADDMISSINGSEQS){
				$B ='TRITD_Ttur';
	         $seqs{$B} = 'n' x $width;
	         open(ALIGN,">>",$aligndir.$file) || die "# cannot write to $aligndir$file\n";
		      print ALIGN ">$B\n$seqs{$B}\n";
			   close(ALIGN);
			}
      }		
	}

	#last;
}

