#!usr/bin/perl -w
#This progarm reads in the names of all the input "*.in" files 
#in runsome.dat copies its models and saves them as "*.out"

#Assigning a variable to 'runsome.dat'
$rsome='runsome_purify.dat';
# the existance of this file prevents the automatic test suite from running */
$blockfile = "blockfile.txt";

# pure coverage options
# coverage /SaveMergeData /SaveMergeTextData exename.exe

#opening 'runsome.dat'
open(RSOM,"$rsome");

# open blocking file, will delete it at end 
open(BLOCK,">$blockfile");
printf( BLOCK "blocking file opened\n");
close(BLOCK);

#The command to be executed
$exe='cp';
$exe = "c:/Projects/Cloudy/trunk/debug/trunk.exe";

print "\n";

#Reading in the input file names that are to be executed
while($r=<RSOM>)
{
   # check is beginning of line (^) is not a word
   if( $r !~/^\W/ )
   {
     $r=~s/\n//; 	  #removing any newline characters
     $out=$r;	  #giving name of output file same as the input file.
     $out=~s/\.in/.out/gi;
     print "$r goes to $out\n";
     system "purify $exe < $r > $out";   #copy the contents of the .in files to its counterpart .out files
   }
}
close(RSOM);
# remove the old unix source files
# unlink <"$blockfile"> ; 
system( "rm $blockfile");

 

