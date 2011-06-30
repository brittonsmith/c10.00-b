#
# this is called by several scripts in the slow directory to
# announce completion of a set of models
# 
open( ioEMAIL , ">message.txt" );
printf( ioEMAIL " slow run has finished\n");  
close( ioEMAIL );
system("c:\\u\\blat\\blat.exe message.txt -t gary\@pa.uky.edu -s \"slow test suite has completed\" " );

