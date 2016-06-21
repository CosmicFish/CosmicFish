/* sdiff-format output routines for GNU DIFF.

   Copyright (C) 1991, 1992, 1993, 1998, 2001, 2002 Free Software
   Foundation, Inc.

   This file is part of GNU DIFF.

   GNU DIFF is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY.  No author or distributor
   accepts responsibility to anyone for the consequences of using it
   or for whether it serves any particular purpose or works at all,
   unless he says so in writing.  Refer to the GNU DIFF General Public
   License for full details.

   Everyone is granted permission to copy, modify and redistribute
   GNU DIFF, but only under the conditions described in the
   GNU DIFF General Public License.   A copy of this license is
   supposed to have been given to you along with GNU DIFF so you
   can know your rights and responsibilities.  It should be in a
   file named COPYING.  Among other things, the copyright notice
   and this notice must be preserved on all copies.  */

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

static void print_sdiff_common_lines (lin, lin);
static void print_sdiff_hunk (struct change *);

/* Next line number to be printed in the two input files.  */
static lin next0, next1;

/* Print the edit-script SCRIPT as a sdiff style output.  */

void
print_sdiff_script (struct change *script)
{
  next0 = next1 = - files[0].prefix_lines;
  print_script (script, print_sdiff_hunk);

  print_sdiff_common_lines (files[0].valid_lines, files[1].valid_lines);
}

/* Tab from column FROM to column TO, where FROM <= TO.  Yield TO.  */

static unsigned int
tab_from_to (unsigned int from, unsigned int to)
{
  FILE *out = outfile;
  unsigned int tab;

  if (!expand_tabs)
    for (tab = from + TAB_WIDTH - from % TAB_WIDTH;  tab <= to;  tab += TAB_WIDTH)
      {
	putc ('\t', out);
	from = tab;
      }
  while (from++ < to)
    putc (' ', out);
  return to;
}

/*
 * Print the text for half an sdiff line.  This means truncate to width
 * observing tabs, and trim a trailing newline.  Returns the last column
 * written (not the number of chars).
 */
static unsigned int
print_half_line (char const *const *line, unsigned int indent,
		 unsigned int out_bound)
{
  FILE *out = outfile;
  register unsigned int in_position = 0;
  register unsigned int out_position = 0;
  register char const *text_pointer = line[0];
  register char const *text_limit = line[1];

  while (text_pointer < text_limit)
    {
      register unsigned char c = *text_pointer++;

      switch (c)
	{
	case '\t':
	  {
	    unsigned int spaces = TAB_WIDTH - in_position % TAB_WIDTH;
	    if (in_position == out_position)
	      {
		unsigned int tabstop = out_position + spaces;
		if ((expand_tabs))
		  {
		    if (out_bound < tabstop)
		      tabstop = out_bound;
		    for (;  out_position < tabstop;  out_position++)
		      putc (' ', out);
		  }
		else
		  if (tabstop < out_bound)
		    {
		      out_position = tabstop;
		      putc (c, out);
		    }
	      }
	    in_position += spaces;
	  }
	  break;

	case '\r':
	  {
	    putc (c, out);
	    tab_from_to (0, indent);
	    in_position = out_position = 0;
	  }
	  break;

	case '\b':
	  if (in_position != 0 && --in_position < out_bound)
	    {
	      if (out_position <= in_position)
		/* Add spaces to make up for suppressed tab past out_bound.  */
		for (;  out_position < in_position;  out_position++)
		  putc (' ', out);
	      else
		{
		  out_position = in_position;
		  putc (c, out);
		}
	    }
	  break;

	case '\f':
	case '\v':
	control_char:
	  if (in_position < out_bound)
	    putc (c, out);
	  break;

	default:
	  if (! ISPRINT (c))
	    goto control_char;
	  /* falls through */
	case ' ':
	  if (in_position++ < out_bound)
	    {
	      out_position = in_position;
	      putc (c, out);
	    }
	  break;

	case '\n':
	  return out_position;
	}
    }

  return out_position;
}

/*
 * Print side by side lines with a separator in the middle.
 * 0 parameters are taken to indicate white space text.
 * Blank lines that can easily be caught are reduced to a single newline.
 */

static void
print_1sdiff_line (char const *const *left, char sep,
		   char const *const *right)
{
  FILE *out = outfile;
  unsigned int hw = sdiff_half_width, c2o = sdiff_column2_offset;
  unsigned int col = 0;
  bool put_newline = 0;

  if (left)
    {
      put_newline |= left[1][-1] == '\n';
      col = print_half_line (left, 0, hw);
    }

  if (sep != ' ')
    {
      col = tab_from_to (col, (hw + c2o - 1) / 2) + 1;
      if (sep == '|' && put_newline != (right[1][-1] == '\n'))
	sep = put_newline ? '/' : '\\';
      putc (sep, out);
    }

  if (right)
    {
      put_newline |= right[1][-1] == '\n';
      if (**right != '\n')
	{
	  col = tab_from_to (col, c2o);
	  print_half_line (right, col, hw);
	}
    }

  if (put_newline)
    putc ('\n', out);
}

/* Print lines common to both files in side-by-side format.  */
static void
print_sdiff_common_lines (lin limit0, lin limit1)
{
  lin i0 = next0, i1 = next1;
/*   long len0, len1; */

  if (i0 != limit0 || i1 != limit1)
    {
      if (!suppress_common_lines)
	{
	  while (i0 != limit0 && i1 != limit1)
	    print_1sdiff_line (&files[0].linbuf[i0++], ' ',
			       &files[1].linbuf[i1++]);
	  while (i1 != limit1)
	    print_1sdiff_line (0, ')', &files[1].linbuf[i1++]);
	  
	  while (i0 != limit0)
	    print_1sdiff_line (&files[0].linbuf[i0++], '(', 0);
	}
/*       else */
/* 	{ */
/* 	  len0 = limit0 - i0; */
/* 	  len1 = limit1 - i1; */
/* 	  fprintf (outfile, "i%ld,%ld\n", len0, len1); */
/* 	} */
    }

  next0 = limit0;
  next1 = limit1;
}

/* Print a hunk of an sdiff diff.
   This is a contiguous portion of a complete edit script,
   describing changes in consecutive lines.  */

static void
print_sdiff_hunk (struct change *hunk)
{
  lin first0, last0, first1, last1;
  register lin i, j;

  /* Determine range of line numbers involved in each file.  */
  enum changes changes =
    analyze_hunk (hunk, &first0, &last0, &first1, &last1);
  if (!changes)
    return;

  /* Print out lines up to this change.  */
  print_sdiff_common_lines (first0, first1);

/*   if ((suppress_common_lines)) */
/*     { */
/*       long len0 = last0 - first0 + 1; */
/*       long len1 = last1 - first1 + 1; */
/*       fprintf (outfile, "c%ld,%ld\n", len0, len1); */
/*     } */

  /* Print ``xxx  |  xxx '' lines */
  if (changes == CHANGED)
    {
      for (i = first0, j = first1;  i <= last0 && j <= last1;  i++, j++)
	print_1sdiff_line (&files[0].linbuf[i], '|', &files[1].linbuf[j]);
      changes = (i <= last0 ? OLD : 0) + (j <= last1 ? NEW : 0);
      next0 = first0 = i;
      next1 = first1 = j;
    }

  /* Print ``     >  xxx '' lines */
  if (changes & NEW)
    {
      for (j = first1; j <= last1; ++j)
	print_1sdiff_line (0, '>', &files[1].linbuf[j]);
      next1 = j;
    }

  /* Print ``xxx  <     '' lines */
  if (changes & OLD)
    {
      for (i = first0; i <= last0; ++i)
	print_1sdiff_line (&files[0].linbuf[i], '<', 0);
      next0 = i;
    }
}

/*
 * This function does the same actions as print_half_line(),
 * but takes as first argument a simple pointer to
 * const char
 */
static unsigned int
display_half_line (const char* line, unsigned int indent,
		   unsigned int out_bound)
{
  FILE *out = outfile;
  register unsigned int in_position = 0;
  register unsigned int out_position = 0;
  register const char *text_pointer = line;
  register const char *text_limit = line;

  for (; *text_limit != '\0'; text_limit++);
  while (text_pointer < text_limit)
    {
      register unsigned char c = *text_pointer++;

      switch (c)
	{
	case '\t':
	  {
	    unsigned int spaces = TAB_WIDTH - in_position % TAB_WIDTH;
	    if (in_position == out_position)
	      {
		unsigned int tabstop = out_position + spaces;
		if ((expand_tabs))
		  {
		    if (out_bound < tabstop)
		      tabstop = out_bound;
		    for (;  out_position < tabstop;  out_position++)
		      putc (' ', out);
		  }
		else
		  if (tabstop < out_bound)
		    {
		      out_position = tabstop;
		      putc (c, out);
		    }
	      }
	    in_position += spaces;
	  }
	  break;

	case '\r':
	  {
	    putc (c, out);
	    tab_from_to (0, indent);
	    in_position = out_position = 0;
	  }
	  break;

	case '\b':
	  if (in_position != 0 && --in_position < out_bound)
	    {
	      if (out_position <= in_position)
		/* Add spaces to make up for suppressed tab past out_bound.  */
		for (;  out_position < in_position;  out_position++)
		  putc (' ', out);
	      else
		{
		  out_position = in_position;
		  putc (c, out);
		}
	    }
	  break;

	case '\f':
	case '\v':
	control_character:
	  if (in_position < out_bound)
	    putc (c, out);
	  break;

	default:
	  if (! ISPRINT (c))
	    goto control_character;
	  /* falls through */
	case ' ':
	  if (in_position++ < out_bound)
	    {
	      out_position = in_position;
	      putc (c, out);
	    }
	  break;

	case '\n':
	  return out_position;
	}
    }

  return out_position;
}

/*
 * Print side by side lines with a separator in the middle.
 * 0 parameters are taken to indicate white space text.
 * Blank lines that can easily be caught are reduced to a single newline.
 * This function works almost like  print_1sdiff_line() (see above),
 * but the parameter list is slightly different.
 * In addition, different separators are used
 */

void
print_1overview_line (const char *left, int are_different,
		      const char *right)
{
  FILE *out = outfile;
  unsigned int hw = sdiff_half_width, c2o = sdiff_column2_offset;
  unsigned int col = 0;
  bool put_newline = 0; 
  register const char *end_of_left = left;
  register const char *end_of_right = right;
  const char *sep = ":!:";

  if (left)
    {
      for (; *end_of_left != '\0'; end_of_left++);
      if (end_of_left != left)
	end_of_left--;
      put_newline |= *end_of_left == '\n';
      col = display_half_line (left, 0, hw);
    }
  else
    sep = ":>:";

  if (!right)
    sep = ":<:";

  if ((are_different))
    {
      col = tab_from_to (col, (hw + c2o - 3) / 2) + 3;
      fputs (sep, out);
    }

  if (right)
    {
      for (; *end_of_right != '\0'; end_of_right++);
      if (end_of_right != right)
	end_of_right--;
      put_newline |= (*end_of_right == '\n');
      if (*right != '\n')
	{
	  col = tab_from_to (col, c2o);
	  display_half_line (right, col, hw);
	}
    }

  if (put_newline)
    putc ('\n', out);
}
