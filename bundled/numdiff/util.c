/* Support routines for GNU DIFF.

   Copyright (C) 1988, 1989, 1992, 1993, 1994, 1995, 1998, 2001, 2002
   Free Software Foundation, Inc.

   This file is part of GNU DIFF.

   GNU DIFF is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2, or (at your option)
   any later version.

   GNU DIFF is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; see the file COPYING.
   If not, write to the Free Software Foundation,
   59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.  */

/*
    This file, coming from the source code of GNU DIFF,
    has been modified by Ivano Primi  <ivprimi@libero.it>
    so that it could be merged into the source code of Numdiff.

    Numdiff is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.

    Numdiff is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
*/

#include "numdiff.h"
#include "linesplit.h"
#include <error.h>
#include <xalloc.h>

/* Use when a system call returns non-zero status.
   NAME should normally be the file name.  */

void
perror_with_name (char const *name)
{
  error (0, errno, "%s", name);
}

/* Use when a system call returns non-zero status and that is fatal.  */

void
pfatal_with_name (char const *name)
{
  int e = errno;
  error (EXIT_TROUBLE, e, "%s", name);
  abort ();
}

/* Compare two lines (typically one from each input file)
   according to the command line options.
   For efficiency, this is invoked only when the lines do not match exactly
   but an option like -i might cause us to ignore the difference.
   Return nonzero if the lines differ.  */

/* 
 * lines_differ() is used in "inout.c" to compare lines coming possibly
 * from the same file.
 */

extern char** def_ifs;

bool
lines_differ (char const *s1, char const *s2, int index1, int index2, argslist* argl)
{
  const unsigned long fieldno_upper_limit = 8*FIELDMASK_SIZE;
  register char const *f1 = s1;
  register char const *f2 = s2;
  char *e1, *e2, ch1, ch2; 
  unsigned long fieldno1, fieldno2;
  int f1_is_num, f2_is_num, f1_is_blurred, f2_is_blurred;

  /* Cache often-used quantities in local variables to help the compiler.  */
  int ignore_case = argl->optmask & _SI_MASK;
  const char** ifs = (const char**)((index1) ? argl->ifs2 : argl->ifs1); 
  const struct numfmt* pnf = (index1) ? &argl->nf2 : &argl->nf1; 
  unsigned char* ghostmask = (index1) ? argl->ghostmask2 : argl->ghostmask1;
  unsigned char* pblurmask = (index1) ? argl->pblurmask2 : argl->pblurmask1;
  unsigned char* tblurmask = (index1) ? argl->tblurmask2 : argl->tblurmask1; /* s1 */

  const char** Ifs = (const char**)((index2) ? argl->ifs2 : argl->ifs1); 
  const struct numfmt* Pnf = (index2) ? &argl->nf2 : &argl->nf1; 
  unsigned char* Ghostmask = (index2) ? argl->ghostmask2 : argl->ghostmask1;
  unsigned char* Pblurmask = (index2) ? argl->pblurmask2 : argl->pblurmask1;
  unsigned char* Tblurmask = (index2) ? argl->tblurmask2 : argl->tblurmask1; /* s2 */

  ifs = (!ifs) ? (const char**)def_ifs : ifs;
  Ifs = (!Ifs) ? (const char**)def_ifs : Ifs;
  f1 = string_spn (f1, ifs, '\n');
  f2 = string_spn (f2, Ifs, '\n');
  fieldno1 = fieldno2 = 0;

  while (*f1 != '\n' && *f2 != '\n')
    {
      /*
	Ignore the fields selected through the option -X
      */
      while ( *f1 != '\n' && fieldno1 < fieldno_upper_limit && (ghostmask[fieldno1 >> 3] & 0x80 >> (fieldno1 & 0x7)) )
	{
	  /* First move `f1' to the begin of the next field */
          e1 = string_cspn (f1, ifs, '\n');
          f1 = string_spn (e1, ifs, '\n');
	  /* and then increment the field index */
	  fieldno1++; 
	}
      if ( fieldno1 >= fieldno_upper_limit )
	{
	  fprintf (stderr, _("***  Fatal error occurred in function %s:\n%s"),
		   __FUNCTION__,
		   _("***  a very long line has been encountered which contains\n***  too many fields to be correctly handled\n"));
	  exit (EXIT_TROUBLE);
	}
      while ( *f2 != '\n' && fieldno2 < fieldno_upper_limit && (Ghostmask[fieldno2 >> 3] & 0x80 >> (fieldno2 & 0x7)) )
	{
	  /* First move `f2' to the begin of the next field */
          e2 = string_cspn (f2, Ifs, '\n');
          f2 = string_spn (e2, Ifs, '\n');
	  /* and then increment the field index */
	  fieldno2++; 
	}
      if ( fieldno2 >= fieldno_upper_limit )
	{
	  fprintf (stderr, _("***  Fatal error occurred in function %s:\n%s"),
		   __FUNCTION__,
		   _("***  a very long line has been encountered which contains\n***  too many fields to be correctly handled\n"));
	  exit (EXIT_TROUBLE);
	}
      if (*f1 != '\n' && *f2 != '\n')
	{
	  /* Find the ends of the fields */
          e1 = string_cspn (f1, ifs, '\n');
          e2 = string_cspn (f2, Ifs, '\n');
	  /*
	    Mark the ends of the fields before calling acxnum().
	    But before doing this, save the original characters placed
	    just after the fields.
	  */
	  ch1 = *e1; ch2 = *e2;
	  *e1 = *e2 = '\0';
	  /* Determine the types of the fields */
	  f1_is_num = acxnum (f1, pnf) >= e1;
	  f2_is_num = acxnum (f2, Pnf) >= e2;
	  /* Determine whether the fields are blurred or not */
	  f1_is_blurred = (tblurmask[fieldno1 >> 3] & 0x80 >> (fieldno1 & 0x7)) ||
	    ((pblurmask[fieldno1 >> 3] & 0x80 >> (fieldno1 & 0x7)) && (f1_is_num));
	  f2_is_blurred = (Tblurmask[fieldno2 >> 3] & 0x80 >> (fieldno2 & 0x7)) ||
	    ((Pblurmask[fieldno2 >> 3] & 0x80 >> (fieldno2 & 0x7)) && (f2_is_num));
	  if (f1_is_blurred != f2_is_blurred)
	    return 1;
	  else if ((f1_is_blurred))
	    {
	      /*
		`f1' and `f2' are both blurred 
	       */
	      *e1 = ch1;
	      *e2 = ch2;
	    }
	  else
	    {
	      /*
		Neither `f1' nor `f2' is blurred
	       */
	      if (f1_is_num != f2_is_num)
		/*
		  If one field is numeric but the other one not,
		  then the two lines are surely different :)
		 */
		return 1;
	      else if ( (f1_is_num) )
		{
		  /*
		    If both fields are numeric,
		    then perform a reduction
		    to a standard numerical format before
		    the byte-by-byte comparison.
		  */
		  if ( (compare_numeric_strings (f1, pnf, f2, Pnf)) )
		    return 1;
		  else
		    {
		      *e1 = ch1;
		      *e2 = ch2;
		    }
		}
	      else
		{
		  /*
		    If the fields are not numeric,
		    then go on with byte-by-byte comparison,
		  */
		  *e1 = ch1;
		  *e2 = ch2;
		  if ((ignore_case))
		    for (; f1 < e1 && f2 < e2 && TOLOWER(*f1) == TOLOWER(*f2); f1++, f2++);
		  else
		    for (; f1 < e1 && f2 < e2 && *f1 == *f2; f1++, f2++);
		  if (f1 < e1 || f2 < e2)
		    return 1;
		  /*
		    else: We have automatically
		    f1 == e1 && f2 == e2
		  */
		}
	    } 
	  /*
	    Move to the next field
	    and increase the field number
	   */
          f1 = string_spn (e1, ifs, '\n');
          f2 = string_spn (e2, Ifs, '\n');
	  fieldno1++; 
	  fieldno2++;
	} /* End else --> if (*f1 != '\n' && *f2!= '\n') */
    } /* end  while (*f1 != '\n' && *f2 != '\n') */

  /*
    Ignore the fields selected through the option -X
  */
  while ( *f1 != '\n' && fieldno1 < fieldno_upper_limit && (ghostmask[fieldno1 >> 3] & 0x80 >> (fieldno1 & 0x7)) )
    {
      /* First move `f1' to the begin of the next field */
      e1 = string_cspn (f1, ifs, '\n');
      f1 = string_spn (e1, ifs, '\n');
      /* and then increment the field index */
      fieldno1++; 
    }
  if ( fieldno1 >= fieldno_upper_limit )
    {
      fprintf (stderr, _("***  Fatal error occurred in function %s:\n%s"),
	       __FUNCTION__,
	       _("***  a very long line has been encountered which contains\n***  too many fields to be correctly handled\n"));
      exit (EXIT_TROUBLE);
    }  

  while ( *f2 != '\n' && fieldno2 < fieldno_upper_limit && (Ghostmask[fieldno2 >> 3] & 0x80 >> (fieldno2 & 0x7)) )
    {
      /* First move `f2' to the begin of the next field */
      e2 = string_cspn (f2, Ifs, '\n');
      f2 = string_spn (e2, Ifs, '\n');
      /* and then increment the field index */
      fieldno2++; 
    }
  if ( fieldno2 >= fieldno_upper_limit )
    {
      fprintf (stderr, _("***  Fatal error occurred in function %s:\n%s"),
	       __FUNCTION__,
	       _("***  a very long line has been encountered which contains\n***  too many fields to be correctly handled\n"));
      exit (EXIT_TROUBLE);
    }

  return (*f1 != *f2);
}

/* Divide SCRIPT into pieces by calling HUNKFUN and
   print each piece with PRINTFUN.
   Both functions take one arg, an edit script.

   HUNKFUN is called with the tail of the script
   and returns the last link that belongs together with the start
   of the tail.

   PRINTFUN takes a subscript which belongs together (with a null
   link at the end) and prints it.  */

void
print_script (struct change *script,
	      void (*printfun) (struct change *))
{
  struct change *next = script;
  struct change *this;

  while ((next))
    {
      this = next;

      /* Disconnect them from the rest of the changes,
	 making them a hunk, and remember the rest for next iteration.  */
      next = this->link;
      this->link = 0;
#ifdef _DEBUG_SCRIPT_
      debug_script (this);
#endif

      /* Print this hunk.  */
      (*printfun) (this);

      /* Reconnect the script so it will all be freed properly.  */
      this->link = next;
    }
}

/* Look at a hunk of edit script and report the range of lines in each file
   that it applies to.  HUNK is the start of the hunk, which is a chain
   of `struct change'.  The first and last line numbers of file 0 are stored in
   *FIRST0 and *LAST0, and likewise for file 1 in *FIRST1 and *LAST1.
   Note that these are internal line numbers that count from 0.

   If no lines from file 0 are deleted, then FIRST0 is LAST0+1.

   Return UNCHANGED if only ignorable lines are inserted or deleted,
   OLD if lines of file 0 are deleted,
   NEW if lines of file 1 are inserted,
   and CHANGED if both kinds of changes are found. */

enum changes
analyze_hunk (struct change *hunk,
	      lin *first0, lin *last0,
	      lin *first1, lin *last1)
{
  struct change *next;
  lin l0, l1;
  lin show_from, show_to;

  show_from = show_to = 0;

  *first0 = hunk->line0;
  *first1 = hunk->line1;

  next = hunk;
  do
    {
      l0 = next->line0 + next->deleted - 1;
      l1 = next->line1 + next->inserted - 1;
      show_from += next->deleted;
      show_to += next->inserted;
    }
  while ((next = next->link) != 0);

  *last0 = l0;
  *last1 = l1;

  return (show_from ? OLD : UNCHANGED) | (show_to ? NEW : UNCHANGED);
}

/* Yield a new block of SIZE bytes, initialized to zero.  */

void *
zalloc (size_t size)
{
  void *p = xmalloc (size);
  memset (p, 0, size);
  return p;
}

#ifdef _DEBUG_SCRIPT_

void
debug_script (struct change *sp)
{
  fflush (stdout);

  for (; sp; sp = sp->link)
    {
      long line0 = sp->line0;
      long line1 = sp->line1;
      long deleted = sp->deleted;
      long inserted = sp->inserted;
      fprintf (stderr, "%3ld %3ld delete %ld insert %ld\n",
	       line0, line1, deleted, inserted);
    }

  fflush (stderr);
}

#endif
