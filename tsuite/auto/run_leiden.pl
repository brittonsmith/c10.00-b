#!perl -w
# the above is generic form of #! for perl, with warnings
#
# uncommnent following line to run on unix machine, where *.pl is not
# associated with perl 
# this is needed to submit perl script in batch mode on sdx
# #!/usr/bin/perl/
#
# This is a perl script to run all the input files in a directory.
# this IS NOT the main auto run script that happens every night - that is autorun.pl
#
# to submit this script in batch mode on the sdx make the #!perl the first line,
# submit to serial queue with bsub -q serial runall.pl

# It is designed to be executed from the directory where the input files live.
# The string "$exe" must be modified to point to your executable version
# of cloudy.

# this script could be run as simply
# perl runall.pl, or even 
# runall.pl

# after all the models are run the "checkall.pl" script can be run
# to check that everything went ok

# in perl a comment starts with the "#" character, so this is a comment
# variable names start with a "$", so $exe below is a variable, the
# path to the executable


# if there is already an exe file in this dir, back it up before
# bringing in a new one, in case we need to run the old models

if(-e "trunk.exe" )
{
   system("copy \"trunk.exe\" \"cloudy_bak.exe\" ");
}

# this is the path to the executable version of the code
# change the string between the double quotes to the location
# of the cloudy executable on your system
#
# this path is valid on my account on sdx
# $exe = "/u/home1/gary/cloudy/trunk/main.exe";
#
# this is the path on my XP box
#$exe = "c:/Projects/Cloudy/trunk/debug/trunk.exe";
#system("cp c:/projects/cloudy/trunk/debug/trunk.exe . ");
system("cp c:/projects/cloudy/trunk/release/trunk.exe . ");
#$exe = "c:/projects/Cloudy/trunk/release/trunk.exe";
$exe = "trunk.exe";

$nMod = 0;
$nskip = 0;

# this loops over all the *.in files in the current directory
# and runs the code to produce *.out files
while ( defined( $input = glob("pdr_leiden*.in") ) )
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
