#!/usr/bin/perl -w
#
# this script looks for problems in the Cloudy test suite it looks for botched
# asserts, warnings and for simulations that did not end. all serious errors
# will be reported in serious.txt all lines with PROBLEM will be in minor.txt
# sims that are close to an assert (but OK) are in close.txt those with DEBUG
# print statements or TestCode calls are in debug.txt
#
# execute this script from within the directory where the *.out files live.
# The command is: ./checkall.pl
#
# Peter van Hoof

system "../auto/checkall.pl";
