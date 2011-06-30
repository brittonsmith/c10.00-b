#!perl -w
$exe = "c:/Projects/Cloudy/trunk/debug/trunk.exe";

# count number of sims we computed
$nMod = 0;

# set number of sims to skip
$nSkip = 0;

# this loops over all the *.in files in the current directory
# and runs the code to produce *.out files
while ( defined( $input = glob("*.in") ) )
{
	if ( $nMod>$nSkip )
	{
		print( "$input going to " );
		$output = $input;
		$output =~ s/\.in/.out/gi;
		print("$output\n");
		# actually execute the code using purify
		system "nice -n 5 coverage $exe < $input  > $output";
		# actually execute the code using insure
		# system "nice -n 5 inject $exe < $input  > $output";
	}
	++$nMod;
}

print("\n=========================\n");
printf( "\n %i models were computed, and %i were skipped.\n ",$nMod,$nSkip);
print("Now use the checkall.pl script to check results.\n");
print("=========================\n");
