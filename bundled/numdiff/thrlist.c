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
#include "numdiff.h"
#include "xalloc.h"

extern Real Zero;

/* Care that internally the field numbers start from zero, not from one */
const long FIELDNO_MAX = (8 * FIELDMASK_SIZE);
const char separator = ':';

static const struct numfmt defaults = {CURRENCY,
				       DP, 
				       THSEP, 
				       GROUPING, 
				       POS_SIGN, 
				       NEG_SIGN, 
				       ECH, 
				       IU};

/*
  Create a new stack of threshold values, then
  push onto the stack the element 

  {
   threshold = 0.0, beg1 = 0, end1 = FIELDNO_MAX-1,
                    beg2 = 0, end2 = FIELDNO_MAX-1,
   double_range_spec = 0
  }

  and return a pointer to the top of the stack.
 */

thrlist thrlist_new (void)
{
  thrlist list;

  list = (thrlist) xmalloc (sizeof (thrlist_node));
  initR (&list->threshold);
  /* copyR (&list->threshold, Zero); */
  list->beg1 = list->beg2 = 0;
  list->end1 = list->end2 = FIELDNO_MAX - 1;
  list->double_range_spec = 0;
  list->next = NULL;
  return list;
}


/*
  If the string STR begins with a valid range of positive integer values, then
  set *B to (start value - 1), *E to (end value - 1), and, if TAIL != NULL, 
  set *TAIL so that it points to the first character of STR after the range
  specification. Finally, return THRLIST_OK (0).
  For example, if STR == "1-10#abc", then *B will be set to 0, *E to 9, and
  *TAIL will point to the character '#'.
  If the string STR does not begin with a valid range of positive integer values,
  then leave *B, *E and *TAIL unchanged and return the error code THRLIST_INVALID_FORMAT.
  
  Remarks:

  1. The following abbreviated range specifications are admitted, where N is a positive
     integer value and M == FIELDNO_MAX :

    N    stays for N-N
   -N    stays for 1-N
    N-   stays for N-M .

  2. The extremal values of a valid range are integer numbers between
     1 and FIELDNO_MAX. In addition, the start value may not be greater than
     the end value. In case these conditions are not both true
     the function returns THRLIST_INVALID_FORMAT.
*/

static int set_interval (const char *str, unsigned long *b, unsigned long *e, char **tail)
{
  long beg, end;
  char *ptr, *endptr;

  beg = end = -1;
  if (!str || !*str)
    return THRLIST_INVALID_FORMAT; /* illegal input */

  /* If we arrive here we are sure that *str != '\0' ! */
  if ((beg = strtol (str, &endptr, 10)) == 0
      || beg > FIELDNO_MAX || beg < -FIELDNO_MAX)
    return THRLIST_INVALID_FORMAT; /* illegal input */
  else if (beg < 0)
    {
      if (*endptr == '\0' || *endptr == separator)
	{
	  end = -beg;
	  beg = 1;
	}
      else
	return THRLIST_INVALID_FORMAT;
    }
  else if (*endptr == '\0' || *endptr == separator)
    end = beg;
  else if (*endptr == '-')
    {
      ptr = endptr + 1;
      if (*ptr == '\0' || *ptr == separator)
	end = FIELDNO_MAX; 
      else
	{
	  if ((end = strtol (ptr, &endptr, 10)) <= 0
	      || (*endptr != '\0' && *endptr != separator) || end > FIELDNO_MAX)
	    return THRLIST_INVALID_FORMAT; /* illegal input */
	}
    }
  if (beg > end)
    return THRLIST_INVALID_FORMAT;
  else
    {
      *b = (unsigned long)(beg-1);
      *e = (unsigned long)(end-1);
      if (tail != NULL)
	*tail = endptr;
      return THRLIST_OK;
    }
}


/*
  If the string DEF contains a specification of the form

  RANGE1:RANGE2,

  where the ranges RANGE1 and RANGE2 are both in the form accepted by
  the function set_interval() and have the same length, then
  set *BEG1, *END1, *BEG2 and *END2 to the extremal
  values of RANGE1 and RANGE2, respectively, after
  decrementing them of 1.
  Then set *DBL_RNG_SPEC to 1 and return THRLIST_OK 
  to the calling function.

  If the string DEF contains a specification of the form

  RANGE,

  where RANGE is in the form accepted by the function 
  set_interval(), then set *BEG1 and *BEG2 to 
  (start_value_of_RANGE - 1), *END1 and *END2 
  (end_value_of_RANGE - 1), and *DBL_RNG_SPEC to zero.
  Then return THRLIST_OK to the calling function.
  
  If DEF does not contain a valid specification of the form

  RANGE  or  RANGE1:RANGE2

  then return THRLIST_INVALID_FORMAT.
  If the given specification is valid, but
  RANGE1 and RANGE2 do not have the same length,
  then return THRLIST_INVALID_RANGES.

  Remark: In the last two cases after returning 
  from the function the values of

  *BEG1, *BEG2, *END1, *END2 and *DBL_RNG_SPEC

  can not be trusted.
 */

static int set_intervals (const char *def, 
			  unsigned long *beg1, unsigned long *beg2,
			  unsigned long *end1, unsigned long *end2,
			  int *dbl_rng_spec)
{
  char *endptr;
  int rv;
  
  rv = set_interval (def, beg1, end1, &endptr);
  if (rv != 0)
    return rv;
  else if (*endptr == separator)
    {
      rv = set_interval (endptr + 1, beg2, end2, &endptr);
      if ((rv != 0) || (*endptr != '\0'))
	return THRLIST_INVALID_FORMAT;
      else
	{
	  *dbl_rng_spec = 1;
	  return (*end2 + *beg1 == *end1 + *beg2) ? THRLIST_OK : THRLIST_INVALID_RANGES;
	}
    }
  else 
    {
      /* *endptr == '\0' */
      *beg2 = *beg1;
      *end2 = *end1;
      *dbl_rng_spec = 0;
      return THRLIST_OK;
    }
}


/*
  Push onto the stack pointed to by PLIST a new element
  according to the specification contained in the string
  DEF.
  
  This specification must have the form

  1.          THRESHOLD             or

  2.          THRESHOLD:RANGE       or

  3.          THRESHOLD:RANGE1:RANGE2

  where THRESHOLD is a non-negative real (decimal) number,
  and RANGE or RANGE1:RANGE2 must be a valid specification for 
  the function set_intervals().

  In case 1. and 2. set the field DOUBLE_RANGE_SPEC of
  the just pushed element to zero, otherwise set it to 1. 

  If DEF does not contain a valid specification, then
  do nothing but return a suitable error code
  (either THRLIST_INVALID_FORMAT or THRLIST_INVALID_RANGES).
  Otherwise, return THRLIST_OK.
 */

int thrlist_add (thrlist *plist, const char* def)
{
  int rv, dbl_rng_spec;
  Real r;
  unsigned long beg1, beg2, end1, end2;
  thrlist_node *pnode;
  char *endptr;
  
  if (!def || !*def)
    return THRLIST_INVALID_FORMAT; /* no input */

  /* If we arrive here we are sure that *def != '\0' ! */
  initR (&r);
  str2R (def, &endptr, ISCALE, &defaults, &r);
  if (endptr == def || cmp(r, Zero) < 0)
    {
      delR (&r);
      return THRLIST_INVALID_FORMAT; /* no valid positive number */
    }
  else
    {
      if (*endptr == '\0')
	{
	  beg1 = beg2 = 0;
	  end1 = end2 = FIELDNO_MAX - 1;
	  dbl_rng_spec = 0;
	}
      else if (*endptr == separator)
	{
	  if ( (rv = set_intervals (endptr + 1, &beg1, &beg2, &end1, &end2, &dbl_rng_spec)) != 0)
	    {
	      delR (&r);
	      return rv;
	    }
	}
      else
	{
	  delR (&r);
	  return THRLIST_INVALID_FORMAT;
	}
      pnode = (thrlist_node*) xmalloc (sizeof (thrlist_node));
      initR (&pnode->threshold);
      copyR (&pnode->threshold, r);
      delR (&r);
      pnode->beg1 = beg1;
      pnode->beg2 = beg2;
      pnode->end1 = end1;
      pnode->end2 = end2;
      pnode->double_range_spec = dbl_rng_spec;
      pnode->next = *plist;
      *plist = pnode;
      return THRLIST_OK;
    }
}


/*
  Look in the stack LIST for the first element E from the top such that
  either E.DOUBLE_RANGE_SPEC == 0 and
  E.BEG1 <= FIELDNO1 <= E.END1, E.BEG2 <= FIELDNO2 <= E.END2, or
  E.DOUBLE_RANGE_SPEC == 1 and
  E.BEG1 <= FIELDNO1 <= E.END1, FIELDNO2 == E.BEG2 + FIELDNO1 - E.BEG1 .
  Then return the value of the comparison test between R and E.THRESHOLD.

  If there is no element E in LIST satisfying one of the previous conditions
  (this should never happen if LIST has been created by thrlist_new()), 
  then interrupt the main program via exit(EXIT_TROUBLE) after printing
  a suitable error message on stdout.
 */

int thrlist_cmp (Real r, thrlist list, unsigned long fieldno1, unsigned long fieldno2)
{
  thrlist_node *pnode, *pnext;

  pnode = list;
  while (pnode != NULL)
    {
      pnext = pnode->next; 
      if ( (pnode->double_range_spec) )
	{
	  if (fieldno1 >= pnode->beg1 && fieldno1 <= pnode->end1 &&
	      fieldno2 == pnode->beg2 + fieldno1 - pnode->beg1)
	    {
	      return cmp (r, pnode->threshold); /* found */
	    }
	}
      else
	{
	  if (fieldno1 >= pnode->beg1 && fieldno1 <= pnode->end1 &&
	      fieldno2 >= pnode->beg2 && fieldno2 <= pnode->end2)
	    {
	      return cmp (r, pnode->threshold); /* found */
	    }
	}
      pnode = pnext;
    } 
  /*
    The code execution should NEVER reach this point.
   */
  printf (_("Fatal error occurred during comparison of two numerical fields\n"));
  exit (EXIT_TROUBLE);
}


/*
  Dispose the stack pointed to by PLIST (free and clean the memory).
 */

void thrlist_dispose (thrlist *plist)
{
  thrlist_node *pnode, *pnext;

  if (plist != NULL)
    {
      pnode = *plist;
      while (pnode != NULL)
	{
	  pnext = pnode->next; 
	  delR (&pnode->threshold);
	  free ((void*)pnode);
	  pnode = pnext;
	}
      *plist = NULL;
    }
}
