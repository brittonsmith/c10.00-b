#!perl -w

# this is the path to the executable version of the code
# change the string between the double quotes to the location
# of the cloudy executable on your system
#
# this path is valid on my account on sdx
# $exe = "/u/home1/gary/cloudy/trunk/main.exe";
#
# this is the path on my XP box
$exe = "c:/projects/cloudy/trunk/release/trunk.exe";
$exe = "c:/projects/cloudy/trunk/debug/trunk.exe";

$nMod = 0;
$nskip = 0;

# this loops over all the *.in files in the current directory
# and runs the code to produce *.out files
while ( defined( $input = glob("coll*.in") ) )
{
	if ( $nMod>=$nskip )
	{
                # $output =~ s/.out//gi;
                # $output =~ s/.out//gi;
		print( "$input going to " );
		$output = $input;
		$output =~ s/\.in//gi;
		if( -e "$output".".out" )
		{
			rename( "$output".".out" , "$output".".bak" );
		}
                $out = "$output".".out";
		print("$output\n");
		# actually execute the code
		system "nice -n 5 $exe < $input  > $out";
#		system "$exe < $input  > $out";
	}
	++$nMod;
}

print("\n=========================\n");
printf( "\n %i models were computed, and %i were skipped.\n ",$nMod,$nskip);
print("Now use the checkall.pl script to check results.\n");
print("=========================\n");
