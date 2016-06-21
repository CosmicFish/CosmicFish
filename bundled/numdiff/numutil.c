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

#include "numdiff.h"
#include<stdio.h>
#include<stdlib.h> /* for exit() */
#include<string.h>
#if defined(HAVE_LOCALECONV) && defined(USE_GMP)
#include<locale.h>
#endif /* HAVE_LOCALECONV and USE_GMP */

/*
  The following structure is used to store 
  the canonical form of a real number. 
  The canonical form corresponding to zero is
  .sgn = 0
  .ipart = NULL
  .dpart = NULL
  .expn  = 0     .
  
  If X is a non-null number with sign S (either + or -),
  integer part II...II (each I is a digit but the digits have
  not to be necessarily all equal), decimal
  part DD...DD (each D is a digit but the digits have
  not to be necessarily all equal) and exponent N,
  i.e.  X = SII...II(decimal point)DD...DD * 10^N,
  then the canonical form of X is given by:

  .sgn = S
  .ipart = "II...II"
  .dpart = "DD...DD"

  and, denoted by NI the number of digits
  of the integer part and by NN the number
  of null digits after the decimal point
  (NN can be zero, for instance if the first
   digit of the decimal part is not null),

           +----   N + NI - 1, if NI > 0  
           |
  .expn  = <
           |
           +----   N - NN - 1, if NI == 0

  The exponent stored in the canonical form
  is then equal to the one that you would find
  in the scientific representation of the number.
  If the integer (decimal) part of X is zero,
  then .ipart (.dpart respectively) is the empty
  string.

  The string that one obtains by appending the decimal 
  part to the integer part is what I call mantissa.

  Given two non-null numbers X and Y both
  in canonical form, they are equal if and only if
  they have the same sign, the same exponents 
  and the same mantissa.
*/
struct canform {
  unsigned char sgn;
  char *ipart, *dpart;
  long expn;
};

struct canform null = { 0, NULL, NULL, 0};

static char* 
anum (const char *str, const struct numfmt* pnf)
{
  char *ptr;
  long digits, strscale;
  size_t length;

  ptr = (char*)str;
  digits = 0;
  strscale = 0;
  length = strlen(pnf->currency); 

  if ( (*ptr == pnf->pos_sign) || (*ptr == pnf->neg_sign) )
    move_ahead(ptr);
  if ( length > 0 && strncmp (ptr, pnf->currency, length) == 0 )
    ptr += length; /* Skip the currency name, if specified */
  if (pnf->grouping > 0)
    {
      while ( (is_digit((int)*ptr)) || *ptr == pnf->thsep )
	{
	  int first_sep = 1;

	  if (*ptr == pnf->thsep)
	    {
	      unsigned long i;
	      char* ptr2;

	      if ((first_sep))
		{
		  for (ptr2 = ptr;
		       ptr2 > str && (is_digit ((int)*(ptr2-1)));
		       ptr2--);
		  if ((i=ptr-ptr2) == 0 || i > pnf->grouping)
		    {
		      return (char*) str;
		    }
		  first_sep = 0;
		}
	      for (ptr2 = ptr + 1; (is_digit((int)*ptr2)); ptr2++);
	      if ((i=ptr2-ptr-1) != pnf->grouping)
		{
		  return (char*) str;
		}
	      ptr++;
	    }
	  else
	    {
	      move_ahead(ptr), digits++;
	    }
	}
    }
  else
    {
      while ( (is_digit((int)*ptr)) )
	{
	  move_ahead(ptr), digits++;
	}
    }
  if (*ptr == pnf->dp)
    move_ahead(ptr);

  while ( (is_digit((int)*ptr)) )
    move_ahead(ptr), strscale++;
  if (digits+strscale == 0)
    return (char*)str;
  if (TOLOWER(*ptr) == TOLOWER(pnf->ech) && !is_space (*(ptr+1)))
    {
      char *tail;
      long expn;

      expn = strtol (ptr + 1, &tail, 10);
      if (expn < MIN_EXPN && tail != ptr + 1)
	{
	  fprintf (stderr, _("%s: a number with a too small exponent has been found,\nnamely \"%s\".\n"), PACKAGE, str);
	  fprintf (stderr, _("Exponents smaller than %ld are not accepted,\n"), MIN_EXPN);
	  fprintf (stderr, _("the execution of the program ends now\n"));
	  exit (EXIT_TROUBLE);
	}
      if (expn > MAX_EXPN && tail != ptr + 1)
	{
	  fprintf (stderr, _("%s: a number with a too large exponent has been found,\nnamely \"%s\".\n"), PACKAGE, str);
	  fprintf (stderr, _("Exponents larger than %ld are not accepted,\n"), MAX_EXPN);
	  fprintf (stderr, _("the execution of the program ends now\n"));
	  exit (EXIT_TROUBLE);
	}
      return (tail != ptr + 1 ? tail : ptr);
    }
  else
    return ptr;
}

/*
  This function is going to be used in inout.c and util.c
*/
char* 
acxnum (const char *str, const struct numfmt* pnf)
{
  char *ptr, *endrp, *ptr2;
  endrp = ptr = anum (str, pnf);

  if (ptr != str)
    {
      if (*ptr == pnf->iu)
	{
	  return ptr + 1;
	}
      else
	{
	  while (is_space (*ptr))
	    ptr++;
	  if (*ptr == POS_SIGN || *ptr == NEG_SIGN)
	    {
	      ptr2 = anum (ptr, pnf);
	      return (*ptr2 != pnf->iu ? endrp : ptr2 + 1);
	    }
	  else
	    return endrp;
	}
    }
  else
    return ptr;
}

static char* 
snum (const char *str, const struct numfmt* pnf,
      struct canform* pcform)
{
  char *ptr = (char*)str;
  char *it = pcform->ipart;
  int nonnull_digit_found = 0;
  long digits, strscale, nzeros, expn;
  char* tail;
  size_t length;

  digits = strscale = 0;
  length = strlen(pnf->currency); 
  if ( (*ptr == pnf->pos_sign) )
    {
      pcform->sgn = POS_SIGN;
      move_ahead(ptr);
    }
  else if ( (*ptr == pnf->neg_sign) )
    {
      pcform->sgn = NEG_SIGN;
      move_ahead(ptr);
    }
  if ( length > 0 && strncmp (ptr, pnf->currency, length) == 0 )
    ptr += length; /* Skip the currency name, if specified */
  if (pnf->grouping > 0)
    {
      while ( (is_digit((int)*ptr)) || *ptr == pnf->thsep )
	{
	  int first_sep = 1;

	  if (*ptr == pnf->thsep)
	    {
	      unsigned long i;
	      char* ptr2;

	      if ((first_sep))
		{
		  for (ptr2 = ptr;
		       ptr2 > str && (is_digit ((int)*(ptr2-1)));
		       ptr2--);
		  if ((i=ptr-ptr2) == 0 || i > pnf->grouping)
		    {
		      return (char*) str;
		    }
		  first_sep = 0;
		}
	      for (ptr2 = ptr + 1; (is_digit((int)*ptr2)); ptr2++);
	      if ((i=ptr2-ptr-1) != pnf->grouping)
		{
		  return (char*) str;
		}
	      ptr++;
	    }
	  else
	    {
	      if ((nonnull_digit_found |= *ptr > CHAR_ZERO))
		*it++ = *ptr;
	      move_ahead(ptr), digits++;
	    }
	}
    }
  else
    {
      while ( (is_digit((int)*ptr)) )
	{
	  if ((nonnull_digit_found |= *ptr > CHAR_ZERO))
	    *it++ = *ptr;
	  move_ahead(ptr), digits++;
	}
    }
  if (*ptr == pnf->dp)
    {
      it = pcform->dpart;
      move_ahead(ptr);
    }
  while ( (is_digit((int)*ptr)) )
    {
      for (nzeros = 0; *ptr == CHAR_ZERO; nzeros++, move_ahead(ptr));
      strscale += nzeros;
      if ( (is_digit((int)*ptr)) )
	{
	  nonnull_digit_found = 1;
	  for (; nzeros > 0; nzeros--)
	    *it++ = CHAR_ZERO;
	  *it++ = *ptr;
	  move_ahead(ptr), strscale++;
	}
    }
  if (digits+strscale == 0)
    return (char*)str;
  if (!nonnull_digit_found)
    {
      free((void*)pcform->ipart);
      free((void*)pcform->dpart);
      pcform->ipart = pcform->dpart = NULL;
      pcform->sgn = 0;
    }
  if (TOLOWER(*ptr) == TOLOWER(pnf->ech) && !is_space (*(ptr+1)))
    {
      expn = strtol (ptr + 1, &tail, 10);
      if (tail != ptr + 1)
	{
	  if (expn < MIN_EXPN)
	    {
	      fprintf (stderr, _("%s: a number with a too small exponent has been found,\nnamely \"%s\".\n"), PACKAGE, str);
	      fprintf (stderr, _("Exponents smaller than %ld are not accepted,\n"), MIN_EXPN);
	      fprintf (stderr, _("the execution of the program ends now\n"));
	      exit (EXIT_TROUBLE);
	    }
	  if (expn > MAX_EXPN)
	    {
	      fprintf (stderr, _("%s: a number with a too large exponent has been found,\nnamely \"%s\".\n"), PACKAGE, str);
	      fprintf (stderr, _("Exponents larger than %ld are not accepted,\n"), MAX_EXPN);
	      fprintf (stderr, _("the execution of the program ends now\n"));
	      exit (EXIT_TROUBLE);
	    }
	  pcform->expn = (nonnull_digit_found) ? expn : 0;
	  return tail;
	}
      else
	{
	  pcform->expn = 0;
	  return ptr;
	}
    }
  else
    return ptr;
}

static int 
inc_with_check (long* pl, unsigned long r)
{
  long l = *pl;
  unsigned long L;
  
  if (r == 0)
    return 0;
  else if (r > LONG_MAX)
    return -1;
  /*
    After this we are sure that 1 <= r <= LONG_MAX 
   */
  else if ( (long)r >= -l ) /* This is surely the case if l >= 0 */
    {
      L = l + r;
      if (L > LONG_MAX)
	return -1;
      else
	{
	  *pl = L;
	  return 0;
	}
    }
  else
    {
      /* r < -l  ==>  l < r+l < 0 */
      *pl = l + (long)r;
      return 0;
    }
}

static int 
dec_with_check (long* pl, unsigned long r)
{
  long l = *pl;
  unsigned long L;
  
  if (r == 0)
    return 0;
  else if (r > LONG_MAX)
    return -1;
  /*
    After this we are sure that 1 <= r <= LONG_MAX 
   */
  else if ( (long)r >= l ) /* This is surely the case if l <= 0 */
    {
      L = r - l;
      if (L > LONG_MAX)
	return -1;
      else
	{
	  /* 0 <= r - l <= LONG_MAX   ==>  LONG_MIN <= -LONG_MAX <= l - r <= 0 */
	  *pl = -(long)L;
	  return 0;
	}
    }
  else
    {
      /* r < l  ==>  0 < l-r < l */
      *pl = l - (long)r;
      return 0;
    }
}

static void 
msg_and_abort (const char* progname)
{
  fprintf (stderr, 
	   "%s: number with too big mantissa or\nwith exponent out of the allowed range\n",
	   progname);
  exit (EXIT_TROUBLE);
}

static void 
normalize (struct canform* pnum)
{
  size_t lip, lnullprefix;
  char *ptr;

  if (!pnum || !pnum->ipart || !pnum->dpart)
    return;
  if (pnum->sgn == 0)
    pnum->sgn = POS_SIGN;
  if ( (lip = strlen (pnum->ipart)) > ULONG_MAX )
    msg_and_abort (PACKAGE);

  if (lip != 0)
    {
      if ( inc_with_check (&pnum->expn, lip-1) != 0 )
	msg_and_abort (PACKAGE);
      if ( *pnum->dpart == '\0' )
	for (ptr = pnum->ipart + lip - 1; *ptr == CHAR_ZERO; *ptr-- = '\0');
    }
  else
    {
      /* The integer part has zero length ==> decimal part with positive length */
      for (lnullprefix = 0, ptr = pnum->dpart;
	   *ptr == CHAR_ZERO; lnullprefix++, ptr++);
      if ( lnullprefix > ULONG_MAX )
	msg_and_abort (PACKAGE);
      if ( dec_with_check (&pnum->expn, lnullprefix+1) != 0 )
	msg_and_abort (PACKAGE);
      if ((lnullprefix))
	{
	  for (ptr = pnum->dpart + lnullprefix; *ptr != '\0'; ptr++)
	    *(ptr-lnullprefix) = *ptr;
	  *(ptr-lnullprefix) = '\0';
	}
    }
}

static char* 
scxnum (const char *str, const struct numfmt* pnf, 
	struct canform* pre, struct canform* pim)
{
  char *ptr, *endrp, *ptr2;
  struct canform cform1, cform2;
  size_t slen = strlen (str);

  cform1.sgn   = 0;
  cform1.ipart = stralloc (slen);
  cform1.dpart = stralloc (slen);
  cform1.expn  = 0;  
  endrp = ptr = snum (str, pnf, &cform1);
#ifdef _SCXNUM_DEBUG_
  fprintf (stderr, "%s -->\n sign = \'%c\'\n ip = \"%s\"\n dp = \"%s\"\n expn = %ld, tail = \"%s\n",
	   str, cform1.sgn, cform1.ipart, cform1.dpart, cform1.expn, endrp);
#endif
  if (ptr != str)
    {
      normalize (&cform1);
#ifdef _SCXNUM_DEBUG_
      fprintf (stderr, "After normalization:\n");    
      fprintf (stderr, "%s -->\n sign = \'%c\'\n ip = \"%s\"\n dp = \"%s\"\n expn = %ld, tail = \"%s\n",
	       str, cform1.sgn, cform1.ipart, cform1.dpart, cform1.expn, endrp);
#endif
      if (*ptr == pnf->iu)
	{
	  /*
	    We have read a pure imaginary number.
	   */
	  *pre = null;
	  *pim = cform1;
	  return ptr + 1;
	}
      else
	{
	  /*
	    We have to check if we have read 
	    a pure real number or if we can read a full
	    complex number.
	  */
	  *pre = cform1;
	  while (is_space (*ptr))
	    ptr++;
	  if (*ptr == POS_SIGN || *ptr == NEG_SIGN)
	    {
	      cform2.sgn   = 0;
	      cform2.ipart = stralloc (slen);
	      cform2.dpart = stralloc (slen);
	      cform2.expn  = 0;
	      ptr2 = snum (ptr, pnf, &cform2);
#ifdef _SCXNUM_DEBUG_
	      fprintf (stderr, "%s -->\n sign = \'%c\'\n ip = \"%s\"\n dp = \"%s\"\n expn = %ld, tail = \"%s\n",
		       ptr, cform2.sgn, cform2.ipart, cform2.dpart, cform2.expn, ptr2);
#endif
	      if (*ptr2 == pnf->iu)
		{
		  /*
		    We have read a full complex number.
		   */
		  normalize (&cform2);
#ifdef _SCXNUM_DEBUG_
		  fprintf (stderr, "After normalization:\n");    
		  fprintf (stderr, "%s -->\n sign = \'%c\'\n ip = \"%s\"\n dp = \"%s\"\n expn = %ld, tail = \"%s\n",
			   ptr, cform2.sgn, cform2.ipart, cform2.dpart, cform2.expn, ptr2);
#endif
		  *pim = cform2;
		  return ptr2 + 1;
		}
	      else
		{
		  /*
		    We have read a pure real number.
		   */
		  free ((void*) cform2.ipart);
		  free ((void*) cform2.dpart);
		  *pim = null;
		  return endrp;
		}
	    }
	  else
	    {
	      /*
		We have read a pure real number.
	      */
	      *pim = null;
	      return endrp;
	    }
	}
    }
  else
    {
      /*
	We have read no valid number,
	clean and return.
       */
      free ((void*) cform1.ipart);
      free ((void*) cform1.dpart);
      *pre = *pim = null;
      return ptr;
    }
}

/*
  The following function compares two canonic forms of numerical
  values. It returns 0 if they are equal, 1 if they differ.
*/
static int 
compare_canonic_forms (struct canform x, struct canform y)
{
  char *px, *py;
  int  ppx, ppy;
  
  if (!x.ipart && !x.dpart)
    {
      /* If `x' is zero, then */
      return ((y.ipart) || (y.dpart) ? 1 : 0);
      /* return 1 if `y' is not zero, */
      /* 0 if also `y' is zero        */
    }
  else
    {
      /* If `x' is not zero, then */
      if (!y.ipart && !y.dpart)
	return 1;
      /* return 1 if `y' is zero */

      /* If both `x' and `y' are not zero */
      /* we need a detailed comparison    */
      if (x.sgn != y.sgn)
	return 1; 
      if (x.expn != y.expn)
	return 1;
      /* After comparing signs and exponents */
      /* we compare the mantissas            */
      ppx = ppy = 0;
      px = x.ipart;
      py = y.ipart;
      /*
	Explanation about the meaning of `ppx' and `ppy':
	`ppx' (`ppy') is 0 as long as `px' (`py') points to
	a character in the string `x.ipart' (`y.ipart'),
	`ppx' (`ppy') is 1 when `px' (`py') points to
	a character in the string `x.dpart' (`y.dpart'),
	and becomes 2 after that both `x.ipart' and `x.dpart'
	(`y.ipart' and `y.dpart') have been scanned
      */
      do
	{
	  for ( ;
		*px != '\0' && *py != '\0' && *px == *py;
		px++, py++);
	  if (*px == '\0')
	    {
	      ppx++;
	      if (ppx == 1)
		px = x.dpart;
	    }
	  if (*py == '\0')
	    {
	      ppy++;
	      if (ppy == 1)
		py = y.dpart;
	    }
	} while (*px == *py && ppx < 2 && ppy < 2);
      return (*px != *py);
    }
}

/*
  The following function compares `str1' and `str2' assuming that they
  both contain a valid (possibly complex) number and neglecting differences
  in the numeric format. To be sure that the contents of `str1' and `str2'
  correspond to valid numbers, call before acxnum() on `str1' and `str2'.
  This function returns !0 to mean that the numeric contents of `str1' and
  `str2' differ, else 0.

  This function is going to be used in util.c
*/
int 
compare_numeric_strings (const char *str1, const struct numfmt* pnf1,
			 const char *str2, const struct numfmt* pnf2)
{
  struct canform re1, im1, re2, im2;
  int re_differ, im_differ;

  scxnum (str1, pnf1, &re1, &im1);
  scxnum (str2, pnf2, &re2, &im2);
  re_differ = compare_canonic_forms (re1, re2);
  im_differ = compare_canonic_forms (im1, im2);
  /*
    Clean the memory before returning the
    result of the comparison
  */
  if ((re1.ipart))
    free ((void*)re1.ipart);
  if ((re1.dpart))
    free ((void*)re1.dpart);
  if ((im1.ipart))
    free ((void*)im1.ipart);
  if ((im1.dpart))
    free ((void*)im1.dpart);
  if ((re2.ipart))
    free ((void*)re2.ipart);
  if ((re2.dpart))
    free ((void*)re2.dpart);
  if ((im2.ipart))
    free ((void*)im2.ipart);
  if ((im2.dpart))
    free ((void*)im2.dpart);
  return (re_differ || im_differ);
}

#define HNDDULONG 80

static void 
hash_canonic_form (struct canform num, int pos_flag, hash_value* ph)
{
  char *ptr, exp_ch[HNDDULONG+1];
  size_t li, ld, le;
  int s_exp = (num.expn >= 0) ? 1 : -1;
  unsigned long abs_exp = num.expn * s_exp;

  li = (num.ipart) ? strlen(num.ipart) : 0;
  ld = (num.dpart) ? strlen(num.dpart) : 0;
  exp_ch[HNDDULONG] = '\0';
  for (le = 0; abs_exp != 0; le++, abs_exp /= 10)
    exp_ch[HNDDULONG-1-le] = (char)(abs_exp % 10) + CHAR_ZERO;
  if (le > 0)
    exp_ch[HNDDULONG-1-le] = (s_exp < 0) ? NEG_SIGN : POS_SIGN;
  else
    exp_ch[HNDDULONG-1] = CHAR_ZERO;
  s_exp = HNDDULONG-2-le;
  exp_ch[s_exp] = ECH;

  if (li == 0 && ld == 0)
    {
      *ph = HASH (*ph, (unsigned char) CHAR_ZERO);
      *ph = HASH (*ph, (unsigned char) ECH);
      *ph = HASH (*ph, (unsigned char) CHAR_ZERO);
    }
  else
    {
      if (num.sgn == NEG_SIGN)
	*ph = HASH (*ph, (unsigned char) NEG_SIGN);
      else
	{
	  if ((pos_flag))
	    *ph = HASH (*ph, (unsigned char) POS_SIGN);
	}
      if (li > 0)
	{
	  *ph = HASH (*ph, (unsigned char) *num.ipart);
	  if (li > 1 || ld > 0)
	    *ph = HASH (*ph, (unsigned char) DP);
	  for (ptr = num.ipart+1; *ptr != '\0'; ptr++)
	    *ph = HASH (*ph, (unsigned char) *ptr);
	  if (ld > 0)
	    {
	      for (ptr = num.dpart; *ptr != '\0'; ptr++)
		*ph = HASH (*ph, (unsigned char) *ptr);
	    }
	}
      else /* li == 0 but ld > 0 */
	{
	  *ph = HASH (*ph, (unsigned char) *num.dpart);
	  if (ld > 1)
	    {
	      *ph = HASH (*ph, (unsigned char) DP);
	      for (ptr = num.dpart+1; *ptr != '\0'; ptr++)
		*ph = HASH (*ph, (unsigned char) *ptr);
	    }
	}
      for (ptr = exp_ch + s_exp; *ptr != '\0'; ptr++)
	*ph = HASH (*ph, (unsigned char) *ptr);
    }
}

/*
  This function is going to be used in inout.c
*/
char* 
hcxnum (const char *str, const struct numfmt* pnf, hash_value *ph)
{
  struct canform re, im;
  char *endptr;

  endptr = scxnum (str, pnf, &re, &im);
  if (endptr != str)
    {
      if ((im.ipart == NULL && im.dpart == NULL))
	/*
	  If the imaginary part is zero, hash only
	  the real part
	*/
	hash_canonic_form (re, 0, ph);
      else
	{
	  /*
	    If the imaginary part is not zero but
	    the real part is null, then hash
            only the imaginary part
	  */
	  if ((re.ipart == NULL && re.dpart == NULL))
	    hash_canonic_form (im, 0, ph);
	  else
	    {
	      /*
		If both real and imaginary part are not
		zero, hash both of them. The imaginary
		part must always be hashed with its sign
	      */
	      hash_canonic_form (re, 0, ph);
	      hash_canonic_form (im, 1, ph);
	    }
	  *ph = HASH (*ph, IU);
	}
    }
  /*
    Clean the memory and return
  */
  if ((re.ipart))
    free ((void*)re.ipart);
  if ((re.dpart))
    free ((void*)re.dpart);
  if ((im.ipart))
    free ((void*)im.ipart);
  if ((im.dpart))
    free ((void*)im.dpart);
  return endptr;
}

#ifdef USE_GMP

int mpf_a2num (Real* pr, const char *q, char** endptr, const struct numfmt* pnf)
{
  struct canform cform;
  size_t slen = strlen (q);

  cform.sgn   = 0;
  cform.ipart = stralloc (slen);
  cform.dpart = stralloc (slen);
  cform.expn  = 0;  
  *endptr = snum (q, pnf, &cform);
  if (*endptr != q)
    {
      char *ptr, *str; 
      long e;
      size_t len, explen;
#ifdef HAVE_LOCALECONV
      struct lconv *plconv = localeconv(); 
      char* dec_point = plconv->decimal_point;
#else /* not HAVE_LOCALECONV */ 
      char* dec_point = ".";
#endif /* not HAVE_LOCALECONV */ 

      if (cform.ipart == NULL || cform.dpart == NULL)
	{
	  /*
	    The number contained in the string Q is zero
	  */
	  mpf_set_ui (*pr, 0);
	  return 0;
	}
      normalize (&cform);
      e = cform.expn >= 0 ? cform.expn : -cform.expn;
      for (len = 0; e != 0; len++, e /= 10);
      explen = len;
      if (cform.expn < 0)
	/* We need space for the minus sign in front of the exponent */
	len++; 
      if (cform.sgn == NEG_SIGN)
	/* We need space for the minus sign in front of the whole number */
	len++;
      /* We need space also for the decimal point and the exponent letter */
      len += strlen (cform.ipart) + strlen (cform.dpart) + strlen(dec_point) + 1;
      ptr = str = stralloc (len);

      if (cform.sgn == NEG_SIGN)
	{
	  *ptr = NEG_SIGN;
	  ptr++;
	}
      if (*cform.ipart != '\0')
	{
	  *ptr = *cform.ipart;
	  strcat (str, dec_point);
	  strcat (str, cform.ipart+1);
	  strcat (str, cform.dpart);
	}
      else
	{
	  *ptr = *cform.dpart;
	  strcat (str, dec_point);
	  strcat (str, cform.dpart+1);
	}
      if (explen > 0)
	{
	  for (ptr++; *ptr != '\0'; ptr++);
	  *ptr = ECH;
	  if (cform.expn < 0)
	    {
	      ptr++;
	      *ptr = NEG_SIGN;
	    }
	  ptr += explen;
	  for (e = cform.expn >= 0 ? cform.expn : -cform.expn; e != 0; ptr--, e /= 10)
	    *ptr = CHAR_ZERO + e % 10;
	}
      if ( mpf_set_str (*pr, str, 10) == -1 )
	{
	  /* This should never happen.  If mpf_set_str() returns -1, then */
	  /* there is something wrong with the code above.                */
	  fprintf (stderr,
		   _("The string \"%s\"\nis not a valid number, the execution of the program ends now\n"),
		   str);
	  free ((void*) cform.ipart);
	  free ((void*) cform.dpart);
	  free ((void*) str);
	  exit (EXIT_TROUBLE);
	}
      else
	{
	  free ((void*) cform.ipart);
	  free ((void*) cform.dpart);
	  free ((void*) str);
	  return 0;
	}
    } /* *endptr != q */
  else
    {
      /*
	We have read no valid number, then
	free the previously allocated memory,
	set the value pointed to by PR
	to zero and return -1
       */
      free ((void*) cform.ipart);
      free ((void*) cform.dpart);
      mpf_set_ui (*pr, 0);
      return -1;
    } /* *endptr == q */
}

#endif /* USE_GMP */
