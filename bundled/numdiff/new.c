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

#include<stdlib.h>
#include<ctype.h> /* for isspace() */
#include<assert.h>
#include"number.h"
#include"numdiff.h"

#ifdef _DMALLOC_
#include <dmalloc.h> /* Useful only for debugging */
#endif

#ifdef BC_A2NUM_DEBUG
#define move_ahead(ptr) out_char_stderr(*ptr), ptr++
#define set(endptr,str) *endptr = (char*)str, fprintf (stderr, "tail: \"%s\"\n", *endptr)
#else
#define move_ahead(ptr) ptr++
#define set(endptr,str) *endptr = (char*)str
#endif

static
bc_num bc_10_raised_to (long expn, int add_d_dig)
{
  bc_num pw;

  if (expn >= 0)
    {
      pw = bc_new_num (1 + expn, 0);
      pw->n_value[0] = 1;
    }
  else
    {
      pw = bc_new_num (1, -expn + add_d_dig);
      pw->n_value[-expn] = 1;
    }
#ifdef BC_10_DEBUG
  fprintf (stderr, "10^%ld = ", expn);
  pn_stderr (pw);
#endif
  return pw;
}

static
int bc_a2num (bc_num *num, const char *str, char **endptr, int scale,
	      const struct numfmt* pnf)
{
  char *ptr, *nptr;
  int digits, strscale;
  long expn;
  char zero_int;
  size_t length;

#ifdef BC_A2NUM_DEBUG
  fprintf (stderr, "\n\"%s\"\n", str);
#endif
  /* Check for valid number and count digits. */
  ptr = (char*)str;
  digits = 0;
  strscale = 0;
  zero_int = FALSE;
  length = strlen(pnf->currency);
  if ( (*ptr == pnf->pos_sign) || (*ptr == pnf->neg_sign) )
    move_ahead(ptr); /* Sign */
  if ( length > 0 && strncmp (ptr, pnf->currency, length) == 0 )
    ptr += length;   /* Skip the currency name, if specified */
  if (pnf->grouping > 0)
    {
      while ( (is_digit((int)*ptr)) || *ptr == pnf->thsep ) 
	{
	  int first_sep = 1;

	  if (*ptr == pnf->thsep)
	    {
	      unsigned i;
	      char* ptr2;

	      if ((first_sep))
		{
		  /* We have to check that before the separator there is */
		  /* at least one digit but no more than 'pnf->grouping' */
		  for (i=0, ptr2 = ptr; 
		       ptr2 > str && (is_digit ((int)*(ptr2-1)));
		       i++, ptr2--);
		  if (i == 0 || i > pnf->grouping)
		    {
		      if ((endptr))
			set(endptr,str);
		      return -1;
		    }
		  first_sep = 0;
		}
	      /* We have to check that after the separator   */
	      /* there are exactly 'pnf->grouping' digits    */
	      for (i = 0, ptr2 = ptr + 1; (is_digit((int)*ptr2)); 
		   i++, ptr2++);
	      if (i != pnf->grouping)
		{
		  if ((endptr))
		    set(endptr,str);
		  return -1;
		}
	      ptr++;
	    }
	  else /* We have just found a new digit */
	    move_ahead(ptr), digits++;
	}
    }
  else /* pnf->grouping == 0 */
    {
      while ( (is_digit((int)*ptr)) ) 
	move_ahead(ptr), digits++;	/* digits */
    }
  /*   OLD CODE                         */
  /*   while ( (is_digit((int)*ptr)) )  */
  /*     ptr++, digits++;	        */
  if (*ptr == pnf->dp) 
    move_ahead(ptr);		/* decimal point */
  while (is_digit((int)*ptr)) 
    move_ahead(ptr), strscale++;	/* digits */
  if (digits+strscale == 0)
    {
      if ((endptr))
	set(endptr,str);
      *num = bc_copy_num (_zero_);
      return -1;
    }
  if (TOLOWER(*ptr) == TOLOWER(pnf->ech) && !is_space (*(ptr+1)))
    {
      char *tail;

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
      if ((endptr))
	*endptr = tail != ptr + 1 ? tail : ptr; 
    }
  else
    {
      expn = 0;
      if ((endptr))
	*endptr = ptr;
    }
#ifdef BC_A2NUM_DEBUG
  fprintf (stderr, "%c%ld\n", pnf->ech, expn);
#endif
  /* Adjust numbers and allocate storage and initialize fields. */
  strscale = MIN(strscale, scale);
  if (digits == 0)
    {
      zero_int = TRUE;
      digits = 1;
    }
  *num = bc_new_num (digits, strscale);

  /* Build the whole number. */
  ptr = (char*)str;
  if (*ptr == pnf->neg_sign)
    {
      (*num)->n_sign = MINUS;
      ptr++;
    }
  else
    {
      (*num)->n_sign = PLUS;
      if (*ptr == pnf->pos_sign) ptr++;
    }
  if ( length > 0 && strncmp (ptr, pnf->currency, length) == 0 )
    ptr += length;   /* Skip the currency name, if specified */
  nptr = (*num)->n_value;
  if ((zero_int))
    {
      *nptr++ = 0;
      digits = 0;
    }
  for (; digits > 0; ptr++)
    {
      /*
	We have to take into account the
	presence of the thousands separator !
      */
      if (*ptr != pnf->thsep)
	{
	  /* *ptr is a digit! */
	  *nptr++ = CH_VAL(*ptr);
	  digits--;
	}
    }
  /*
    OLD CODE (to build the integer part)
  
    for (;digits > 0; digits--)
    *nptr++ = CH_VAL(*ptr++);
  */

  /* Build the fractional part. */
  if (strscale > 0)
    {
      ptr++;  /* skip the decimal point! */
      for (;strscale > 0; strscale--)
	*nptr++ = CH_VAL(*ptr++);
    }
  
  /* Take into account the exponent */
  if((expn))
    {
      bc_num Factor, Prod;

      bc_init_num (&Prod);
      Factor = expn > 0 ? 
	bc_10_raised_to (expn, 0) : bc_10_raised_to (expn, scale);

      bc_multiply (*num, Factor, &Prod, strscale);
      *num = bc_copy_num (Prod);

      bc_free_num (&Factor);
      bc_free_num (&Prod);
    }
#ifdef BC_A2NUM_DEBUG
  pn_stderr (*num);
  if((endptr))
    fprintf (stderr, "tail: \"%s\"\n", *endptr);
#endif
  return 0;
} 

#define DIM 50

static
void bc_print_num (bc_num num, void (* out_char)(int), int prec)
{
  /* The negative sign if needed. */
  if (num->n_sign == MINUS) (*out_char) ('-');
  /* Output the number. */
  if (num->n_value[0] == 10)
    {
      (*out_char) ('I');
      (*out_char) ('n');
      (*out_char) ('f');
    }
  else if (bc_is_zero (num))
    {
      (*out_char) ('0');
      (*out_char) ('.');
      for (; prec > 0; prec--)
	(*out_char) ('0');
    }
  else
    {
      char *nptr;
      int index, ndig;
      signed char dig_list[DIM], *buffer, *p;
      long expn;

      assert ((buffer = calloc (prec + 1, sizeof (signed char))));
      nptr = num->n_value;
#ifdef __DEBUG__
      pv ("nptr", nptr, num->n_len + num->n_scale);
      (*out_char)('\n');
#endif
      for (expn = num->n_len; expn > 0 && !*nptr; nptr++, expn--);
      /* Now expn == number of non-null digits */
      /* in the integer part of the number     */
      if (expn == 0)
	{
	  /* Now nptr points to the first digit */
	  /* of the fractional part.            */
	  for (index = num->n_scale, expn = -1;
	       !*nptr; 
	       index--, expn--, nptr++);
	  /* Now expn is the exponent of the number */
	  /* in base 10.                            */
	}  
      else
	{
	  expn--;
	  /* Now expn is the exponent of the number */
	  /* in base 10.                            */
	  
	  
	  for (index = num->n_len + num->n_scale, nptr = num->n_value; 
	       !*nptr; index--, nptr++);
	}  
      
      /* First we print the mantissa */ 
      for (p = buffer, ndig = prec + 1;
	   ndig > 0 && index > 0; 
	   ndig--, index--, p++, nptr++)
	*p = *nptr;
      if (ndig == 0 && index > 0)
	{
	  /* Rounding */
	  if (*nptr >= 5)
	    {
	      for (--p; p > buffer && *p == 9; *p = 0, p--);
	      if (p == buffer)
		{
		  if (*p == 9)
		    {
		      *p = 1;
		      expn++;
		    }
		  else
		    (*p)++;
		}
	      else
		(*p)++;
	    }
	}
      (*out_char)(BCD_CHAR(*buffer));
      (*out_char)('.');
      for (index = 1; index <= prec; index++)
	(*out_char)(BCD_CHAR(buffer[index]));
      free((void*)buffer);
      /* Now it is the moment to print the exponent... */
      (*out_char)('e');
      if (expn < 0)
	{
	  (*out_char) ('-');
	  expn = -expn;
	}
      else
	(*out_char) ('+');
      for (index = 0; index < DIM; dig_list[index] = -1, index++);
      index = 0;
      do
	{
	  dig_list[index] = expn % 10;
	  expn /= 10;
	  index++;
	} while ((expn));
      for (index = 0; dig_list[index] >= 0; index++);
      for (--index; index >= 0; index--)
	(*out_char)(BCD_CHAR(dig_list[index]));
      /* .. done ! */
    } /* End of num != 0 */
}

