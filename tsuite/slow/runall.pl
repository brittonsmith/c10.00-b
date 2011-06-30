#!/usr/bin/perl
#
# This runs all the *.in input files in a directory.
# and should be executed from the directory containing the input files.

# The string "$exe" must be modified to point to your executable version
# of Cloudy.

# This can be done using an argument to runall.pl.  To test the build in source
#
# runall.pl
#
# To test the build in the sys_XXX subdirectory of source
#
# runall.pl sys_XXX

# after all the models are run the "checkall.pl" script can be run
# to check that everything went ok

# This gives the path to the executable version of the code
# This will work if the directory structure was not changed from the
# original download and if the executable was built in source
# using g++ and the makefiles provided
$exe = '../../source/'.$ARGV[0].'/cloudy.exe';

# these are the various compiler builds that occur in sys_XXX
# directories below source.  
# uncomment the one that you built if you did not use the default
#
# sys_gcc
# $exe = "../../source/sys_gcc/cloudy.exe";
# $exe = "../../source/sys_gprof/cloudy.exe";
# $exe = "../../source/sys_IBMxSeries/cloudy.exe";
# $exe = "../../source/sys_icc/cloudy.exe";
# $exe = "../../source/sys_pgcc/cloudy.exe";
# $exe = "../../source/sys_pgccBounds/cloudy.exe";

# if this file does not exist we have problems
if( ! -e "$exe" )
{
	die "The executable file $exe does not exist\n";
}

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
		# actually execute the code including nice
		system "nice -n 5 \"$exe\" < $input  > $out";
	}
	++$nMod;
}

print("\n=========================\n");
printf( "\n %i simulations were computed, and %i were skipped.\n ",$nMod,$nSkip);
print("Now use the checkall.pl script to check results.\n");
print("=========================\n");

system("perl mail.pl\" ");

