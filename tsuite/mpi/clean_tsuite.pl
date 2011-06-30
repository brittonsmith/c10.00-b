#!/usr/bin/perl

# This perl script cleans out the test suite after a test run, it will only
# leave files that are part of the official distribution and the cloudy
# executable (with a name ending in .exe) if that was present!
#
# Peter van Hoof

system "../auto/clean_tsuite.pl";

# when a grid run crashes, these will be left behind...
while( defined( $input = glob("grid?????????_*.in") ) ) {
    print "deleting $input...\n";
    unlink "$input";
}
