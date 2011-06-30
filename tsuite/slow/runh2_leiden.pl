#!/usr/bin/perl
# the above is generic form of #! for perl, with warnings
#
# of the cloudy executable on your system
#
# this path is valid on my account on sdx
# $exe = "/u/home1/gary/cloudy/trunk/main.exe";
#

if(-e "trunk.exe" )
{
   system("copy \"trunk.exe\" \"cloudy_bak.exe\" ");
}

#system("copy l:/projects/cloudy/trunk/release/trunk.exe c:/projects/cloudy/trunk/release/trunk.exe ");
#system("copy l:\\projects\\cloudy\\trunk\\release\\trunk.exe c:\\projects\\cloudy\\trunk\\release\\trunk.exe ");
#system("copy l:\\projects\\cloudy\\trunk\\release\\trunk.exe . ");
system("copy c:\\projects\\cloudy\\trunk\\release\\trunk.exe . ");

# copy fast exe from lenticular to here
#system("copy l:/projects/cloudy/trunk/release/trunk.exe c:/projects/cloudy/trunk/release/trunk.exe ");
#system("copy l:\\projects\\cloudy\\trunk\\release\\trunk.exe c:\\projects\\cloudy\\trunk\\release\\trunk.exe ");
#$exe = "c:/projects/cloudy/trunk/vs6/release/vs6.exe";
# this is the path on my XP box
#$exe = "c:/Projects/Cloudy/trunk/debug/trunk.exe";
#$exe = "c:/projects/Cloudy/trunk/release/trunk.exe";
$exe = "trunk.exe";

$nMod = 0;
$nSkip = 0;

# this loops over all the *.in files in the current directory
# and runs the code to produce *.out files
while ( defined( $input = glob("h2_leiden*.in") ) )
{
	if ( $nMod>=$nSkip )
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
printf( "\n %i models were computed, and %i were skipped.\n ",$nMod,$nSkip);
print("Now use the checkall.pl script to check results.\n");
print("=========================\n");

system("perl mail.pl\" ");


