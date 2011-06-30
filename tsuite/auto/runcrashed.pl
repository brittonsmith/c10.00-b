#!/usr/bin/perl
#
# recompute models that did not end
#
# various executables
$exe = "/Users/gary/cloudy/trunk/source/sys_gcc/cloudy.exe";
if( not -e $exe )
{
    print "exec not found at  $exe\n";
    exit;
}

# check all input files for a corresponding output file, and 
# check whether output file contains "Cloudy ends:" string
# write corresponding line into checkend.txt file.
# if size of file $o_dir/checkend.txt is zero, then program crashed,
# and we want to rerun it.  

$nMod = 0;
while ( defined( $input = glob("*.in") ) )
{
	$output = $input;
	$output =~ s/\.in/.out/gi;
	system "grep 'ChkMonitor'  $output > checkend.txt";

#	if zero length, code did  not end
	if ( -z "checkend.txt" )
	{
		++$nMod;

#		print("$output\n");
		printf( STDERR "%s going to %s\n" , $input,$output ); 
#		actually execute the code
#		system "nice -n 5 $exe < $input  > $output";
		system "$exe < $input  > $output";
	}
}
printf( "\n %i models had crashed, and were recomputed.\n ",$nMod);

