#!/usr/bin/perl -w
#
# this script looks for problems in the Cloudy test suite it looks for botched
# monitors, warnings and for simulations that did not end. all serious errors
# will be reported in serious.txt all lines with PROBLEM will be in minor.txt
# sims that are close to thowing a monitor (but OK) are in close.txt those
# with DEBUG print statements or TestCode calls are in debug.txt
#
# execute this script from within the directory where the *.out files live.
# The command is: ./checkall.pl
#
# Peter van Hoof

# this was produced by the optimizers and is not really part of the test suite
# its presence will fool this script....
if( -r "optimal.in" )
{
    unlink "optimal.in";
}

print( "\nThis script checks for problems in the output from the Cloudy test suite.\n");
print( "  It checks for three things in the simulations:\n");
print( "  1) Botched monitors or warnings:  These are serious problems.\n");
print( "  2) Whether the PROBLEM or DEBUG string were printed.\n");
print( "  3) Crashes, in which the code did not end - these are also serious problems.\n");

$nMod = 0;
$nOK = 0;
$nSkip = 0;
$nFail = 0;
$nCrash = 0;
$nBotch = 0;
$nWarn = 0;

open( CRASHED, ">crashed.txt" );
open( SERIOUS, ">serious.txt" );
open( BOTCH, ">serious.asr" );
open( PERFORMANCE, ">performance.asr" );
open( MINOR, ">minor.txt" );
open( CLOSE, ">close.txt" );
open( DEBUG, ">debug.txt" );
open( SKIP, ">skip.txt" );

# loop over all input scripts
while( defined( $input = glob("*.in") ) )
{
    # count the number of simulations
    ++$nMod;

    # form name of output file
    $output = $input;
    $output =~ s/\.in/\.out/gi;

    if( -r $output )
    {
	$checkEnd = 0;
	$nWarnMessage = 0;
	$lgCrash = 0;
	$lgFailure = 0;
	$lgBotch = 0;
	$lgAtmosNotFound = 0;
	open( OUTP, $output );
	while( <OUTP> )
	{
# find "DISASTER", "Botched" and "W-" strings in all the output files and write
# corresponding line into serious.txt file.
	    if( /DISASTER/ || /Botched/ || /W-/ )
	    {
		print SERIOUS "$output: $_";
		if( /W-/ )
		{
		    ++$nWarnMessage;
		}
	    }
	    if( /botch/ )
	    {
		if( /itrz/ || /zone/ || /iter/ )
		{
		    print PERFORMANCE "$output: $_";
		}
		else
		{
		    print BOTCH "$output: $_";
		}
		$lgBotch = 1;
	    }
# now look for string PROBLEM, which indicates an internal problem
# during the calculation - this is not a serious problem
	    if( /PROBLEM/ && ! /DISASTER/ )
	    {
		print MINOR "$output: $_";
	    }
	    if( /ChkMonitor ----/ )
	    {
		print CLOSE "$output: $_";
	    }
	    if( /DEBUG/ || /Test code/ )
	    {
		print DEBUG "$output: $_";
	    }
# check output files to see whether they contain any monitors output
# also determine if the program crashed.  
	    if( /ChkMonitor/ || /Cloudy exited OK\]/ || /something went wrong\]/ )
	    {
		++$checkEnd;
	    }
	    if( /something bad has happened./ )
	    {
		$lgCrash = 1;
	    }
# this is the last line, so lgBotch and nWarnMessage are already determined.
	    if( /something went wrong\]/ && !$lgBotch && $nWarnMessage == 0 )
	    {
		$lgFailure = 1;
	    }
	    if( /Error: stellar atmosphere file not found./ )
	    {
		$lgAtmosNotFound = 1;
	    }
	}
	close( OUTP );
	if( $lgAtmosNotFound )
	{
	    ++$nSkip;
	    print SKIP "$input: stellar atmosphere file is not installed\n";
	}
	elsif( $checkEnd == 0 || $lgCrash )
	{
	    ++$nCrash;
	    print SERIOUS "$output: crashed\n";
	    print CRASHED "$output: crashed\n";
	    $lgSeriousProblems = 1;
	}
	elsif( $lgFailure )
	{
	    ++$nFail;
	    print SERIOUS "$output: returned with FAILURE status\n";
	    $lgSeriousProblems = 1;
	}
	elsif( $nWarnMessage > 0 )
	{
	    ++$nWarn;
	}
	elsif( $lgBotch )
	{
	    ++$nBotch;
	}
	else
	{
	    ++$nOK;
	}
    }
    else
    {
	++$nSkip;
	print SKIP "$input: no corresponding output was found\n";
    }
}

close( CRASHED );
close( SERIOUS );
close( BOTCH );
close( PERFORMANCE );
close( MINOR );
close( CLOSE );
close( DEBUG );
close( SKIP );

if( $nMod != $nOK+$nSkip+$nCrash+$nFail+$nWarn+$nBotch )
{
    print STDERR "\nINTERNAL ERROR: the number of sims don't add up!!!!!\n\n";
}

print( "\nLooking for crashes, botched results and warnings:");
# if this file has non-zero length, we detected a problem
if( -s "serious.txt" )
{
    print "\nWARNING! crashes, failures, botches or warnings were found. This is a serious problem.\n";
    if( $nCrash > 0 )
    {
	print "Cloudy crashed in $nCrash sims.\n";
    }
    print "Check the file serious.txt for names, serious.asr for botched monitors.\n";
    system "grep Cloudy < serious.txt";
    $lgSeriousProblems = 1;
}
else
{
    print " good, none found.\n";
    $lgSeriousProblems = 0;
}

print( "\nLooking for PROBLEM string in the output files:");
if( -s "minor.txt" )
{
    print "\nPROBLEM string was found.  Check minor.txt for a list of names.\n";      
    print "This can occur a few times and is not a serious concern.\n";      
    $lgMinorProblems = 1;
}
else
{
    print " good, none found.\n";
    $lgMinorProblems = 0;
}

if( -s "close.txt" )
{
    print "\nSome sims are close to a monitor.  Check close.txt for a list of names.\n";      
    print "This is for information only, not a problem.\n";      
}

print( "\nLooking for DEBUG string or TestCode call:");
if( -s "debug.txt" )
{
    print "\nDEBUG string or TestCode calls were found.  Check debug.txt for names.\n";      
    $lgDebugPrint = 1;
}
else
{
    print " good, none found.\n";
    $lgDebugPrint = 0;
}
    
# check whether any sims were skipped
print "\nLooking for skipped sims:";
if( $nSkip > 0 )
{
    print " $nSkip sims were skipped - a list is in skip.txt\n" ;
}
else
{
    print " good: no sims were skipped.\n";
}

printf( "\n%i sims are present", $nMod );
if( $nSkip > 0 )
{
    printf( ", %i were skipped", $nSkip );
}
if( $nCrash > 0 )
{
    printf( ", %i crashed", $nCrash );
}
if( $nFail > 0 )
{
    printf( ", %i failed", $nFail );
}
if( $nWarn > 0 )
{
    printf( ", %i have warnings", $nWarn );
}
if( $nBotch > 0 )
{
    printf( ", %i botched", $nBotch );
}
if( $nOK > 0 )
{
    printf( ", %i are OK", $nOK );
}
print ".\n";

if( $lgMinorProblems )
{
    print "Minor problems were found, this is normal, list is in minor.txt\n";
}

if( $lgDebugPrint )
{
    print "DEBUG prints or TestCode calls were found, list is in debug.txt\n";
}

if( $lgSeriousProblems )
{
    print "Serious problems were found. See serious.txt for a summary\n";
}
else
{
    if( $nOK > 0 )
    {
	print "The test suite was successful; you have a valid executable.\n";
    }
}
