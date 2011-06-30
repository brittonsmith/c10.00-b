/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

/* Configuration file specifically for Visual Studio -- every other system
   can generate this automatically so it's sure to be correct... */

#undef HAVE_POWI
#define HAVE_POW_DOUBLE_INT 1
#define HAVE_POW_DOUBLE_LONG 1
#define HAVE_POW_FLOAT_INT 1
#define HAVE_POW_FLOAT_LONG 1
#undef HAVE_POW_DOUBLE_FLOAT
#undef HAVE_POW_FLOAT_DOUBLE
#undef HAVE_FUNC
