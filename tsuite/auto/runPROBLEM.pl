#!perl -w
#   #!/usr/bin/perl/
# first of the above is generic form of #! for perl,
# second is needed to submit perl script in batch mode on ncx
#
# all statement starting with a # are comments
#
# look for models that had botched asserts and recompute them
#

# this is the path to the executable version of the code
# change the string between the double quotes to the location
# of the cloudy executable on your system
#
# this path is valid on my account on ncx
# $exe = "/home/gary/cloudy/trunk/main.exe";
#
# the debug version of the code
#$exe = "c:/Projects/Cloudy/trunk/debug/trunk.exe";
# the fast version of the code
$exe = "c:/projects/Cloudy/trunk/release/trunk.exe";


# check all output files to see whether they contain "PROBLEM"
$nMod = 0;
while(defined($output= glob("*.out")) )
{
        system "grep PROBLEM  $output >checkend.txt";

#	if non-zero length, model had botched asserts
#       -s option returns number of bytes if file exists and has finite size
	if(-s "checkend.txt")
	{
           ++$nMod;

          $input = $output;
          $input =~ s/\.out/.in/gi;
#          print("$output\n");
	  printf( STDERR "%s going to %s\n" , $input,$output ); 
#         actually execute the code
          system "nice -n 5 $exe < $input  > $output";
	}
}
printf( "\n %i models had PROBLEM and were recomputed.\n ",$nMod);

