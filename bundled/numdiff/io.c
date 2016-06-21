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

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include"numdiff.h"

#ifdef _DMALLOC_
#include <dmalloc.h> /* Useful only for the debugging */
#endif

#include"read_line.c"

static
void print_head (const char* str, size_t length)
{
  const char* t;

  for (t = str; t - str < length; putchar(*t), t++);
}

static
void print_ws (unsigned n)
{
  while ((n--))
    putchar (' ');
}

static
void writeln (const char* line, int addemptyline)
{
  const char *ptr;

  for (ptr = line; *ptr != '\0'; putchar(*ptr), ptr++);
  if (ptr == line || *(ptr-1) != '\n')
    puts (EOF_INDICATOR);
  if ((addemptyline))
    putchar('\n');
}

/*
  This function assumes that at least one between line1 and line2 is
  not empty (!NULL).
*/
void print_lines (const char* line1, const char* line2, 
		  unsigned long lineno1, unsigned long lineno2, 
		  int delimiter_only)
{
  puts (LINE_SEP); 
  if (!delimiter_only)
    {
      if (!line1)
	{
	  printf ("          <==\n##%-7lu ==> ", lineno2);
	  writeln (line2, 1);
	}
      else if (!line2)
	{
	  printf ("##%-7lu <== ", lineno1);
	  writeln (line1, 0);
	  printf ("          ==>\n\n");
	}
      else
	{
	  printf ("##%-7lu <== ", lineno1);
	  writeln (line1, 0);
	  printf ("##%-7lu ==> ", lineno2);
	  writeln (line2, 1);
	}
    }
}

/*
  This function assumes that at least one between field1 and field2 is
  not empty.
*/
void print_fields (const char* field1, const char* field2,
		   size_t l1, size_t l2, 
		   unsigned long lineno1, unsigned long lineno2, 
		   unsigned long fieldno1, unsigned long fieldno2)
{
  fieldno1++; /* The field number must start from 1, not from 0 */
  fieldno2++; /* The field number must start from 1, not from 0 */
  if (*field1 == '\0')
    {
      printf ("##%-7lu       <==\n##%-7lu #>%-3lu ==> ", lineno1, lineno2, fieldno2);
      /* length(" #> ") + 3 = 7 white spaces before the arrow */
      writeln (field2, 0);
    }
  else if (*field2 == '\0')
    {
      printf ("##%-7lu #>%-3lu <== ", lineno1, fieldno1);
      writeln (field1, 0);
      /* length(" #> ") + 3 = 7 white spaces before the arrow */
      printf ("##%-7lu       ==>\n", lineno2);
    }
  else
    {
      printf ("##%-7lu #:%-3lu <== ", lineno1, fieldno1);
      print_head (field1, l1);
      putchar ('\n');
      printf ("##%-7lu #:%-3lu ==> ", lineno2, fieldno2);
      print_head (field2, l2);
      putchar ('\n');
    }
}

void print_errors (Real abserr, Real relerr)
{
  fputs (_("@ Absolute error = "), stdout);
  printno (abserr, DEF_LIM);
  fputs (_(", Relative error = "), stdout);
  printno (relerr, DEF_LIM); 
  putchar ('\n');
}

void print_separator (void)
{
  putchar ('@');
  print_ws (53);
  puts ("@@");
}
