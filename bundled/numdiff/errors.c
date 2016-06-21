/*
    Numdiff - compare putatively similar files, 
    ignoring small numeric differences
    Copyright (C) 2005, 2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013  Ivano Primi  <ivprimi@libero.it>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include "number.h"

#ifndef STDC_HEADERS
#define VARARGS 1
#else
#undef VARARGS
#endif

#ifndef VARARGS
#include <stdarg.h>
#else
#include <varargs.h>
#endif

/* Other things for number.c. */
/* int std_only; */

void
out_of_memory()
{
  fprintf (stderr, 
	   _("%s: insufficient memory for new allocation,\nthe execution of the program ends now\n"),
	   PACKAGE);
  exit (EXIT_TROUBLE);
}

/* Runtime error will  print a message and stop the machine. */

#ifndef VARARGS
#ifdef __STDC__
void
rt_error (char *mesg, ...)
#else
void
rt_error (mesg)
     char *mesg;
#endif
#else
void
rt_error (mesg, va_alist)
     char *mesg;
#endif
{
  va_list args;
  char error_mesg [255];

#ifndef VARARGS   
  va_start (args, mesg);
#else
  va_start (args);
#endif
  vsprintf (error_mesg, mesg, args);
  va_end (args);
  
  fprintf (stderr, _("Runtime error: %s\n"), error_mesg);
  exit (EXIT_TROUBLE);
}

/* A runtime warning tells of some action taken by the processor that
   may change the program execution but was not enough of a problem
   to stop the execution. */

#ifndef VARARGS
#ifdef __STDC__
void
rt_warn (char *mesg, ...)
#else
void
rt_warn (mesg)
     char *mesg;
#endif
#else
void
rt_warn (mesg, va_alist)
     char *mesg;
#endif
{
  va_list args;
  char error_mesg [255];

#ifndef VARARGS   
  va_start (args, mesg);
#else
  va_start (args);
#endif
  vsprintf (error_mesg, mesg, args);
  va_end (args);

  fprintf (stderr, _("Runtime warning: %s\n"), error_mesg);
}

