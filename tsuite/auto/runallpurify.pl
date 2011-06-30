#!perl -w
#   #!/usr/bin/perl/
# first of the above is generic form of #! for perl,
# second is needed to submit perl script in batch mode on ncx
#
# This is a perl script to run all the input files in a directory.
# this IS NOT the main auto run script that happens every night - that is autorun.pl
#
# to submit this script in batch mode on the ncx make the #!perl the first line,
# submit to serial queue with bsub -q serial runall.pl

# It is designed to be executed from the directory where the input files live.
# The string "$exe" must be modified to point to your executable version
# of cloudy.

# under most cumstances this script could be run as simply
# perl runall.pl, or even 
# runall.pl

# after all the models are run the "checkall.pl" script can be run
# to check that everything went ok

# in perl a comment starts with the "#" character, so this is a comment
# variable names start with a "$", so $exe below is a variable, the
# path to the executable

# this is the path to the executable version of the code
# change the string between the double quotes to the location
# of the cloudy executable on your system
#
# this path is valid on my account on ncx
# $exe = "/home/gary/cloudy/trunk/main.exe";
#
# this is the path on my XP box
$exe = "c:/Projects/Cloudy/trunk/debug/trunk.exe";

$nMod = 0;
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
		system "nice -n 5 purify $exe < $input  > $output";
		# actually execute the code using insure
		# system "nice -n 5 inject $exe < $input  > $output";
	}
	++$nMod;
}

print("\n=========================\n");
printf( "\n %i models were computed, and %i were skipped.\n ",$nMod,$nSkip);
print("Now use the checkall.pl script to check results.\n");
print("=========================\n");
