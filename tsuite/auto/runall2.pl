#!/usr/bin/perl
#

#$exe = "c:/Users/gary/cloudy/trunk/icl/release/cloudy_icl.exe";
#$exe = "c:/Users/gary/cloudy/trunk/source/sys_icl/cloudy.exe";

#$exe = "c:/Users/gary/Cloudy/trunk/debug/trunk.exe";
#$exe = "c:/Users/gary/cloudy/trunk/release/trunk.exe";
#$exe = "c:/Users/gary/Cloudy/trunk/icl64/release/cloudy_icl64.exe";

# Gary's Mac
$exe = "/Users/gary/cloudy/trunk/source_hot/sys_gcc/cloudy.exe";

# count total number of sims we compute
$nMod = 0;

# sets number to skip - will start tests after skipping the first
# nSkip in alphabetical order.  set to 0 to so all sims
$nSkip = 0;

# total number to do, usually much larger than the total
# to do them all
$nLimit = 700;

# loop over the *.in files
# and runs the test, producinig the *.out files
while ( defined( $input = glob("*.in") ) )
{
	if ( $nMod>=$nSkip  && $nMod < $nLimit )
	{
		print( "$input going to " );
		$output = $input;
		$output =~ s/\.in//gi;
		if( -e "$output".".out" )
		{
			rename( "$output".".out" , "$output".".bak" );
		}
		$out = "$output".".out";
		print("$output\n");
		# actually execute the code, first uses nice
		system "nice -n 5 $exe < $input  > $out";
#		system "$exe < $input  > $out";
	}
	++$nMod;
}

print("\n=========================\n");
printf( "\n %i simulations were computed, and %i were skipped.\n ",$nMod,$nSkip);
print("Now use the checkall.pl script to check results.\n");
print("=========================\n");
