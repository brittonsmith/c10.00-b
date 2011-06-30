#!/usr/bin/perl -w
#This progarm reads in the names of all the input "*.in" files 
#in runsome.dat copies its models and saves them as "*.out"

#Assigning a variable to 'runsome.dat'
$rsome='runsome.dat';

#opening 'runsome.dat'
open(RSOM,"$rsome");

#The command to be executed
#
# this path is valid on my account on sdx
# $exe = "/u/home1/gary/cloudy/trunk/main.exe";
#
# path on the alpha
# this is the path on my XP box
#$exe = "c:/Projects/Cloudy/trunk/icl/release/cloudy_icl.exe";
#$exe = "c:/Projects/Cloudy/trunk/debug/trunk.exe";
#$exe = "c:/projects/Cloudy/trunk/release/trunk.exe";
$exe = "../../source/sys_gcc/cloudy.exe";

print "\n";

#Reading in the input file names that are to be executed
while($r=<RSOM>)
{ 
  $r=~s/\n//; 	  #removing any newline characters
  $out=$r;	  #giving name of output file same as the input file.
  $out=~s/\.in/.out/gi;
  print "$r goes to $out\n";
  system "nice -n 5 $exe < $r > $out";   #copy the contents of the .in files to its counterpart .out files
}
close(RSOM);
 

