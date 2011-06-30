/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland 
and
 * others.  For conditions of distribution and use see copyright notice in licen
se.txt */
/* out-of-line constructor for assert -- put breakpoint in this
   routine to trap assert throws for IDEs without built-in facility. */
#include "cddefines.h"

bad_assert::bad_assert(const char* file, long line, const char* comment):
	p_file(file), p_line(line), p_comment(comment)
{
}
