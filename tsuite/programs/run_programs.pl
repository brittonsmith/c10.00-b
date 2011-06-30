#!/usr/bin/perl
# This progarm reads in the names of all the program files 
# in the file run_programs.dat, compiles and links them, then runs the tests
#
# Some examples of the calling sequence for this script are:
#
# run_programs.pl
# run_programs.pl sys_icc
#
# in the first form it will use the binary in the source directory itself, in
# the second form you can supply the name of one of the sys_xxx subdirectories.

# location of source directory
$DirSource  = "../../source";

if( ! -r "$DirSource/Makefile" )
{
    die "could not find file $DirSource/Makefile\n";
}

# this is where the object files should be
$DirObject  = "$DirSource/$ARGV[0]";

if( ! -d "$DirObject" )
{
    die "could not find directory $DirObject\n";
}

# now determine the compiler and flags to use by reading Makefile and Makefile.conf
$cxx = "";
$opt = "";
$cxxflags = "";

open( FOO, "<$DirSource/Makefile" );
while( <FOO> )
{
    chomp;
    if( /^CXX *\.?=/ )
    {
	s/CXX *\.?= *//;
	$cxx = $_;
    }
    elsif( /^OPT *\.?=/ )
    {
	s/OPT *\.?= *//;
	$opt = $_;
    }
    elsif( /^CXXFLAGS *\.?=/ )
    {
	s/CXXFLAGS *\.?= *//;
	$cxxflags = $_;
    }
}
close( FOO );

if( -r "$DirObject/Makefile.conf" )
{
    open( FOO, "<$DirObject/Makefile.conf" );
    while( <FOO> )
    {
	chomp;
	if( /^CXX *\.?=/ )
	{
	    s/CXX *\.?= *//;
	    $cxx = $_;
	}
	elsif( /^OPT *[\.\+]?=/ )
	{
	    if( /\+=/ )
	    {
		s/OPT *\+= *//;
		$opt = "$opt $_";
	    }
	    else
	    {
		s/OPT *\.?= *//;
		$opt = $_;
	    }
	}
	elsif( /^CXXFLAGS *[\.\+]?=/ )
	{
	    if( /\+=/ )
	    {
		s/CXXFLAGS *\+= *//;
		$cxxflags = "$cxxflags $_";
	    }
	    else
	    {
		s/CXXFLAGS *\.?= *//;
		$cxxflags = $_;
	    }
	}
    }
    close( FOO );
}

# substitute optimization flags in compiler flags
$cxxflags =~ s/\$[{\[]OPT[}\]]/$opt/;

# this file contains the names of the directoriesc containing the programs
$rsome = 'run_programs.dat';
open( RSOM, "<$rsome" ) or die "could not open $rsome\n";

# read the program names that are to be executed, each is in a directory of the same name
while($r=<RSOM>)
{
    # check is beginning of line (^) is not a word
    if( $r !~ /^\W/ )
    {
	# remove trailing newline
	chomp( $r );
	chdir( $r ) or die "could not change to $r\n";
	$command = "$cxx $cxxflags $r.cpp -o $r.exe -I../$DirSource -L../$DirObject -lcloudy";
	print "$command\n";
	system( "$command" );
	print "compilation done.\nnow execute $r.exe\n";
	system( "./$r.exe");
	chdir( ".." ) or die "could not return home\n";
    }
}

close(RSOM);


 

