#!/usr/bin/perl

# look for models that had botched asserts and recompute them

# path to the executable version of the code
$exe = "c:/users/gary/Cloudy/trunk/release/trunk.exe";

# check all output files to see whether they contain "Cloudy ends:"
# string and write corresponding line into checkend.txt file.
# if size of file $o_dir/checkend.txt is zero, then program crashed,
# and we want to rerun it.  

$nMod = 0;
while(defined($output= glob("*.out")) )
{
	#system "grep 'Cloudy ends'  $output >checkend.txt";
        system "grep Botched  $output >checkend.txt";
        system "grep 'W-'  $output >>checkend.txt";

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
#	 system "nice -n 5 $exe < $input  > $output";
          system "$exe < $input  > $output";
	}
}
printf( "\n %i models had botched monitors, and were recomputed.\n ",$nMod);

