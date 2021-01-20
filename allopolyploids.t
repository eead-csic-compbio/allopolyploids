use strict;
use warnings;
use Test::More tests => 3;
use lib "scripts/";

BEGIN { use_ok('polyutils') };

BEGIN { use_ok('polyconfig') };

ok( eval{ `perl scripts/_check_lineages_polyploids.pl -t scripts/controls/brachy.ph` } =~ 
	/[1\s]{9}/ , '_check_lineages_polyploids.pl standalone' );


