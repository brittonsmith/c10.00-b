#!/usr/bin/perl

# This is a little perl script that runs both the auto and slow test suite in
# parallel on suitable machines or clusters. See the file auto/run_parallel.pl
# for a full description of how this script should be used.
#
# Peter van Hoof

do "./run_service.pl";

# this sets the path to various executables, sets up the list of CPU slots,
# and determines the number of CPUs to be used
$nproc = RunService::initialize( "." );

# clean up the test suite
system "./clean_tsuite.pl";

# create list of input scripts to be executed
chomp( $jobListAuto = `cd auto ; ./run_sequence.pl` );
chomp( $jobListSlow = `cd slow ; ./run_sequence.pl` );

$jobList = "";
foreach $input ( split( / +/, $jobListSlow ) ) {
    $jobList .= "slow/$input ";
}
foreach $input ( split( / +/, $jobListAuto ) ) {
    $jobList .= "auto/$input ";
}

print "\nNow I will run the test suite using $nproc processors\n\n";

# this does the heavy lifting...
RunService::run_jobs( $jobList );

if( $nproc > 1 ) {
#     just in case func_trans_read was run before func_trans_punch was complete
    $jobList = "auto/func_trans_read.in";
    RunService::run_jobs( $jobList );
}

# now do the checking
print "\n\n\n";
print "========================================\n";
print "      Checking auto test suite          \n";
print "========================================\n";
system "cd auto ; ./checkall.pl";
print "\n\n\n";
print "========================================\n";
print "      Checking slow test suite          \n";
print "========================================\n";
system "cd slow ; ./checkall.pl";
