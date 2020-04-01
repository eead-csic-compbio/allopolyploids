package polyutils;

# B Contreras-Moreira, R Sancho, EEAD-CSIC & EPS-UNIZAR 2018-20

use strict;
require Exporter;

our @ISA = qw( Exporter );

our @EXPORT = qw(
	get_all_ancestors get_MRCA get_shared_ancestors get_next_MRCA
);

# Takes a tree node object and recursively gets all its ancestors.
# Returns a list of internal ids sorted from immediate to toplevel/root ancestor.
sub get_all_ancestors
{
	my $node = $_[0];
	my @ancestors;

	while(defined($node->ancestor()))
   {
		$node = $node->ancestor();
      push(@ancestors,$node->internal_id());
	}

   return @ancestors;
}

# Receives two referencs to sorted lists of ancestors
# Returns the Most Recent Common Ancestor (MRCA)
sub get_MRCA
{
	my ($ref_ancestorsA, $ref_ancestorsB) = @_;
   my ($nodeA,$nodeB,%hashB);
	my $MRCA = undef;

	# hash elements in list B
	foreach $nodeB (@$ref_ancestorsB)
	{
		$hashB{$nodeB} = 1;	
	}

	# iterate along A
	foreach $nodeA (@$ref_ancestorsA)
   {
		if($hashB{$nodeA})
		{
			$MRCA = $nodeA;
			last;
		}
	}

   return $MRCA;
}

# Takes two sorted lists of ancestors (diploid, polyploid) and
# a reference (ancestor) node.
# Returns a list of shared contiguous descendent nodes that
# terminates with the reference node.
sub get_shared_ancestors
{
   my ($ref_ancestorsA, $ref_ancestorsB, $reference_node) = @_;
   my ($nodeA,$nodeB,%hashB);
   my @shared = ();

	# hash elements in list B (polyploid)
	# Taes 26,27,28,29
	foreach $nodeB (@$ref_ancestorsB)
	{
		$hashB{$nodeB} = 1;
   }

	# iterate along A (diploid) until reference node is found
	# Aspe 25,26,27,28,29
	foreach $nodeA (@$ref_ancestorsA)
   {
		last if($nodeA > $reference_node);

      if($hashB{$nodeA}) # is this a shared node?
      {
			push(@shared, $nodeA);
      }
   }

	return @shared;
}

# Takes a list of sorted ancestors, a reference node and
# a hash reference with (MRCA) nodes as keys.
# Returns the next descendant MRCA node, if any
sub get_next_MRCA
{
	my ($ref_ancestorsA, $reference_node, $ref_MRCA) = @_;
   my $next_MRCA = '';

   # loop backwards along ancestral nodes starting from $reference_node
   # until a MRCA node is found
   foreach my $nodeA (reverse(@$ref_ancestorsA))
   {	
		next if($nodeA >= $reference_node);
      if($ref_MRCA->{$nodeA}){
			$next_MRCA = $nodeA;
			last;
		}
	}

	return $next_MRCA;
}

1;
