#!/usr/bin/perl
# this script checks whether the installed g++ supports precompiling header files
$res = "";
if( $ARGV[0] =~ /^g\+\+/ ) {
    $version = `$ARGV[0] --version`;
#     remove comments in parentheses, they may or may not contain spaces
    $version =~ s/\(.*?\)//g;
#     at this point $version should look like "g++  x.y.z ...."
    @v1 = split( /\s+/, $version );
    @v2 = split( /\./, $v1[1] );
#     $v2[0] is major revision, $v2[1] is minor revision, $v2[2] is microversion
#     precompiling headers is supported from g++ 3.4.0 onwards
    if( ( $v2[0] == 3 && $v2[1] >= 4 ) || $v2[0] >= 4 ) {
	$res .= "precompile ";
    }
#     vectorization is supported from g++ 4.0.0 onwards
    if( $v2[0] >= 4 ) {
	$res .= "vectorize ";
    }
}
# remove trailing spaces
$res =~ s/ +$//;
if( $res ne "" ) {
    print "$res\n";
}
