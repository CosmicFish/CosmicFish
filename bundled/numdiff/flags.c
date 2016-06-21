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

#include"numdiff.h"
#include<xalloc.h>

/* Constants and default values */
#define CHUNK_ALLOC  20480

static flg_array internal_table;

/*
  `table' must be the address of a `flg_array' structure.
*/
static int init_flg_array (flg_array* table)
{
  if ((table))
    {
      table->ptr = NULL;
      table->size = table->len = 0;
      return 0;
    }
  else
    return -1;
}

int init_flags (void)
{
  return init_flg_array (&internal_table);
}

/*
  `table' must be the address of a `flg_array' structure.
*/
static int print_flg_array (FILE* fp, const flg_array* table)
{
  size_t n;
  unsigned short i;
  unsigned char byte, elem;

  if (!table || !table->ptr)
    return (fputs (_("<The array is empty>\n"), fp) == EOF);
  else
    {
      int outerror = 0;

      for (n = 0; !outerror && n < table->len; n += 4U)
	{
	  byte = table->ptr[n/4U];
	  for (i = 0; i < 4U; i++)
	    {
	      if ( (elem = (byte & 0x03 << 2*i) >> 2 * i) )
		{
		  if ( fprintf (fp, "%u", elem) < 0 )
		    {
		      outerror = 1;
		      break;
		    }
		}
	      else
		break;
	    }
	}
      if (!outerror)
	outerror = (fputc ('\n', fp) == EOF);
      return outerror;
    }
}

int print_flags (FILE* fp)
{
  return print_flg_array (fp, &internal_table);
}

/*
  array must be the address of a `flg_array' structure.
*/
static void addnewelem (flg_array* array, unsigned char elem)
{
  if (!array)
    return;
  if ( !array->ptr )
    {
      array->ptr = xcalloc (CHUNK_ALLOC, sizeof (unsigned char));
      array->size = CHUNK_ALLOC;
      array->ptr[0] = elem;
      array->len = 1;
#ifdef __MEMDEBUG__
      fprintf (stderr, "size = %6zu, len = %6zu, string:",
	       array->size, array->len);
      print_flg_array (stderr, array);
#endif
    }
  else
    {
      if (array->len == 4 * array->size - 1)
	{
	  unsigned char* p;

	  array->ptr = xrealloc (array->ptr, array->size + CHUNK_ALLOC);
	  array->size += CHUNK_ALLOC;
	  array->ptr[array->len / 4U] |= elem << 2 * (array->len % 4U);
	  array->len++; /* Now array->len == 4 * array->"old"size */
	  for (p = array->ptr + array->len / 4; p < array->ptr + array->size; *p = 0, p++);
#ifdef __MEMDEBUG__
	  fprintf (stderr, "size = %6zu, len = %6zu, string:",
		   array->size, array->len);
	  print_flg_array (stderr, array);
#endif
	}
      else
	{
	  array->ptr[array->len / 4U] |= elem << 2 * (array->len % 4U);
	  array->len++;
#ifdef __MEMDEBUG__
	  fprintf (stderr, "size = %6zu, len = %6zu, string:",
		   array->size, array->len);
	  print_flg_array (stderr, array);
#endif
	}
    }
}

flg_array copy_of_intflagtab (void)
{
  return internal_table;
}

/*
  array must be the address of a `flg_array' structure.
*/
static void destroy_flg_array (flg_array* array)
{
  if ((array) && (array->ptr))
    {
      free ((void*)array->ptr);
      array->ptr = NULL;
      array->size = array->len = 0;
    }
} 

void erase_flags (void)
{
  destroy_flg_array (&internal_table);
}

static void notedown_sdiff_common_lines (lin, lin);
static void notedown_sdiff_hunk (struct change *);

/* Next line number to be printed in the two input files.  */
static lin next0, next1;

/* 
 * Note down the edit-script SCRIPT in the INTERNAL_TABLE.
 */

void
notedown_sdiff_script (struct change *script)
{
  next0 = next1 = - files[0].prefix_lines;
  print_script (script, notedown_sdiff_hunk);

  notedown_sdiff_common_lines (files[0].valid_lines, files[1].valid_lines);
#ifdef _DEBUG_FLAGSTABLE_
  print_flg_array (stderr, &internal_table);
#endif
}

/* Print lines common to both files in side-by-side format.  */
static void
notedown_sdiff_common_lines (lin limit0, lin limit1)
{
  lin i0 = next0, i1 = next1;

  if (i0 != limit0 || i1 != limit1)
    {
      for (; i0 != limit0 && i1 != limit1; i0++, i1++)
	addnewelem (&internal_table, 3);

      for (; i1 != limit1; i1++)
	addnewelem (&internal_table, 2);

      for (; i0 != limit0; i0++)
	addnewelem (&internal_table, 1);
    }

  next0 = limit0;
  next1 = limit1;
}

/* Note down a hunk of an sdiff diff.
   This is a contiguous portion of a complete edit script,
   describing changes in consecutive lines.  */

static void
notedown_sdiff_hunk (struct change *hunk)
{
  lin first0, last0, first1, last1;
  register lin i, j;

  /* Determine range of line numbers involved in each file.  */
  enum changes changes =
    analyze_hunk (hunk, &first0, &last0, &first1, &last1);
  if (!changes)
    return;

  /* Note down lines up to this change.  */
  notedown_sdiff_common_lines (first0, first1);

  /* Note down ``xxx  |  xxx '' lines */
  if (changes == CHANGED)
    {
      for (i = first0, j = first1;  i <= last0 && j <= last1;  i++, j++)
	addnewelem (&internal_table, 3);
      changes = (i <= last0 ? OLD : 0) + (j <= last1 ? NEW : 0);
      next0 = first0 = i;
      next1 = first1 = j;
    }

  /* Note down ``     >  xxx '' lines */
  if (changes & NEW)
    {
      for (j = first1; j <= last1; ++j)
	addnewelem (&internal_table, 2);
      next1 = j;
    }

  /* Note down ``xxx  <     '' lines */
  if (changes & OLD)
    {
      for (i = first0; i <= last0; ++i)
	addnewelem (&internal_table, 1);
      next0 = i;
    }
}
