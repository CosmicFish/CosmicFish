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
#include"numdiff.h"

#ifdef USE_GMP

#include<limits.h>

Real Zero, Ten, Inf;

void init_mpa(int iscale)
{
  char tmpbuff[1024];

  /*
    Hopefully the highest possible
    number of digits for a LONG
    will never become larger than 1000!! :)
   */
  sprintf (tmpbuff, "1e%ld", LONG_MAX);
  mpf_set_default_prec ((iscale * 32) / 10);
  mpf_init (Zero);
  mpf_init_set_ui (Ten, 10);
  mpf_init_set_str (Inf, tmpbuff, 10);
}

void initR (Real* px)
{
  mpf_init (*px);
}

void initC (Complex* pz)
{
  mpf_init (pz->re);
  mpf_init (pz->im);
}

void copyR (Real* dst, Real src)
{
  mpf_set (*dst, src);
}

void copyC (Complex* dst, Complex src)
{
  mpf_set (dst->re, src.re);
  mpf_set (dst->im, src.im);
}

void str2R (const char *q, char **endptr, int iscale,
	    const struct numfmt* pnf, Real* pr)
{
  mpf_a2num (pr, q, endptr, pnf);
}

void str2C (const char *q, char **endptr, int iscale,
	    const struct numfmt* pnf, Complex* pc)
{
  char *ptr, *ptr2;

  mpf_init (pc->re);
  mpf_init (pc->im);
  mpf_a2num (&pc->re, q, &ptr, pnf);
  if ((endptr))
    *endptr = ptr;
  if (ptr != q)
    {
      if (*ptr == pnf->iu)
	{
	  mpf_set (pc->im, pc->re);
	  mpf_set (pc->re, Zero);
	  if ((endptr))
	    *endptr = ptr + 1;
	}
      else
	{
	  while (is_space (*ptr))
	    ptr++;
	  if (*ptr == POS_SIGN || *ptr == NEG_SIGN)
	    {
	      mpf_a2num (&pc->im, ptr, &ptr2, pnf);
	      if (*ptr2 != pnf->iu)
		mpf_set (pc->im, Zero);
	      else
		{
		  if ((endptr))
		    *endptr = ptr2 + 1;
		}
	    }
	  /*
	     else
	     : we have successfully read a real number
	     but there is no another number after it.
	     : So, we leave pc->im set to zero.
	  */
	}
    }
  else
    /* 
       : q does not contain any valid number
       ==> pc->re is 0. Then we set pc->im to 0.
       We remark that, if endptr is
       not NULL, then *endptr == q.
    */
    mpf_set (pc->im, Zero);
}

void add (Real s, Real t, Real* q, int iscale)
{
  mpf_add (*q, s, t);
}

void square (Real s, Real* q, int iscale)
{
  mpf_pow_ui (*q, s, 2);
}

void divide (Real s, Real t, Real* q, int iscale)
{
  mpf_div (*q, s, t);
}

void divide_by_int (Real* q, int d, int iscale)
{
  Real div;

  if (d == 0)
    mpf_set (*q, Zero);
  else
    {
      mpf_init_set_si (div, d);
      mpf_div (*q, *q, div);
      mpf_clear (div);
    }
}

void square_root (Real* q, int iscale)
{
  mpf_sqrt (*q, *q);
}

void Cabs (Complex z, Real* pm, int iscale)
{
  mpf_t a, b, q;

  if ( mpf_sgn (z.re) >= 0 )
    mpf_init_set (a, z.re);
  else
    {
      mpf_init (a);
      mpf_sub (a, Zero, z.re); 
    }
  if ( (mpf_sgn(z.im) == 0) )
    {
      mpf_set (*pm, a);
      mpf_clear (a);
      return;
    }
  else if ( (mpf_sgn(z.im) > 0) )
    mpf_init_set (b, z.im);
  else
    {
      mpf_init (b);
      mpf_sub (b, Zero, z.im);
    }
  mpf_init (q);
  if ( mpf_cmp (b, a) > 0 )
    {
      mpf_div (q, a, b);
      mpf_mul (a, q, q);
      mpf_add_ui (q, a, 1);
      mpf_sqrt (q, q);
      mpf_mul (*pm, b, q);
    }
  else
    {
      /* a >= b ===> a > 0 */
      mpf_div (q, b, a);
      mpf_mul (b, q, q);
      mpf_add_ui (q, b, 1);
      mpf_sqrt (q, q);
      mpf_mul (*pm, a, q);
    }
  mpf_clear (q);
  mpf_clear (b);
  mpf_clear (a);
}

void Csub (Complex z1, Complex z2, Complex* pw, int iscale)
{
  mpf_sub (pw->re, z1.re, z2.re); 
  mpf_sub (pw->im, z1.im, z2.im); 
}

int cmp (Real p, Real q)
{
  return mpf_cmp (p, q);
}

int is0 (Real u)
{
  return (mpf_cmp (u, Zero) == 0 ? 1 : 0);
}

int smart_cmp (const Complex* pz1, const Complex* pz2, int flag)
{
  if (flag == 0)
    return 1;
  else if (flag > 0)
    return (mpf_cmp (pz1->re, pz2->re) >= 0 &&
	    mpf_cmp (pz1->im, pz2->im) >= 0);
  else /* flag < 0 */
    return (mpf_cmp (pz1->re, pz2->re) <= 0 &&
	    mpf_cmp (pz1->im, pz2->im) <= 0);
}

static int round_far_from_zero (char* mantissa, int prec)
{
  size_t length;
  char *ptr, *abs_mantissa;

  /*
    abs_mantissa is a pointer to the absolute value
    of MANTISSA. If MANTISSA starts with a minus sign,
    then abs_mantissa has to point to the
    digit immediately after this sign.
  */
  if (*mantissa == NEG_SIGN)
    abs_mantissa = mantissa + 1;
  else
    abs_mantissa = mantissa;
  length = strlen (abs_mantissa);
  /*
    If the length of the string ABS_MANTISSA
    is less or equal than PREC+1, then no
    rounding is required and zero has to
    be returned.
  */
  if (length <= prec + 1)
    return 0;
  else
    {
      /*
	If length > prec + 1, then the character
	of ABS_MANTISSA in the position PREC+1 is
	non-null. If this character is a number between
	0 and 4, then the rounding procedure ends hier.
      */
      ptr = abs_mantissa + (prec + 1);
      if ((*ptr - CHAR_ZERO) >= 5)
	{
	  /*
	    Otherwise we have to "increment" the previous
	    character by 1, taking into account
	    each time the amount to be carried.
	  */
	  for (ptr--; ptr > abs_mantissa && *ptr == CHAR_NINE; *ptr = CHAR_ZERO, ptr--);
	  if (ptr == abs_mantissa)
	    {
	      /*
		If the first PREC+1 digits of ABS_MANTISSA
		were all equal to 9, then the rounding
		procedure has to set all of them to zero
		and return 1, i.e. the last carried amount.
	      */
	      if (*ptr == CHAR_NINE)
		{
		  *ptr = CHAR_ZERO;
		  return 1;
		}
	      /*
		If the first digit of ABS_MANTISSA
		is less than 9, then the rounding procedure ends
		by incrementing this digit by 1 and returning 0.
	      */
	      else
		*ptr = *ptr + 1;
	    }
	  /*
	    If at least one of the first PREC+1 digits of ABS_MANTISSA
	    is less than 9, then the rounding procedure ends
	    by incrementing this digit by 1 and returning 0.
	  */
	  else
	    *ptr = *ptr + 1;
	}
      return 0;
    } /* LENGTH > PREC + 1 */
}

static void fprintno (FILE *fp, Real u, int prec)
{
  if (mpf_cmp (u, Inf) == 0)
    printf ("Inf");
  else
    {
      char *ptr, *mantissa;
      mp_exp_t expn;

      ptr = mantissa = mpf_get_str (NULL, &expn, 10, 0, u);
      if (*mantissa == '\0')
	{
	  /*
	    If MANTISSA is the empty string, then
	    print 0 with the required number of decimal
	    digits.
	  */
	  fputc (CHAR_ZERO, fp);
	  fputc (DP, fp);
	  for (; prec > 0; prec--)
	    fputc (CHAR_ZERO, fp);
	  fputc (ECH, fp);
	  fputc (POS_SIGN, fp);
	  fputc (CHAR_ZERO, fp);
	}
      else
	{
	  int amount_to_carry = round_far_from_zero (mantissa, prec);

	  if (*ptr == NEG_SIGN)
	    {
	      /*
		Print a minus sign if the number is negative
	      */
	      fputc (NEG_SIGN, fp);
	      ptr++;
	    }
	  /*
	    The mantissa has been already rounded to the
	    given precision PREC. If MANTISSA before
	    the rounding was "99...9", then now it
	    is given by "00...0" but there is a leading
	    1 which has to be carried.
	  */
	  if (amount_to_carry)
	    {
	      fputc (CHAR_ONE, fp);
	      expn++;
	    }
	  else
	    {
	      fputc (*ptr, fp);
	      ptr++;
	    }
	  /*
	    We print now the decimal point
	  */
	  fputc (DP, fp);
	  /*
	    and the remaining digits of the mantissa.
	    After the decimal point has been printed,
	    one must print no more than PREC digits.
	  */
	  for (; *ptr != '\0' && prec > 0; ptr++, prec--)
	    fputc (*ptr, fp);

	  /*
	    If the remaining digits of the mantissa 
	    are less than PREC, then fill up with zeros.
	  */
	  if (*ptr == '\0')
	    {
	      for (; prec > 0; prec--)
		fputc (CHAR_ZERO, fp);
	    }
	  /*
	    Finally we print the exponent with its sign
	  */
	  fprintf (fp, "%c%+ld", ECH, expn - 1);  
	}
      free ((void*) mantissa);
    } /* u != Inf */ 
}

void printno (Real u, int m)
{
  fprintno (stdout, u, m);
}

#ifdef _MPA_DEBUG
void debug_printno (Real u, int m)
{
  fprintno (stderr, u, m);
}
#endif

void delR (Real* px)
{
  mpf_clear (*px);
}

void delC (Complex* pz)
{
  mpf_clear (pz->re); 
  mpf_clear (pz->im);
}

void end_mpa(void)
{
  mpf_clear (Inf);
  mpf_clear (Ten);
  mpf_clear (Zero);
}

#else /* not USE_GMP */

#include"number.c"
#include"errors.c"
#include"new.c"

Real Zero, Inf;

void init_mpa(int iscale)
{
  bc_init_numbers();
  bc_init_num(&Zero);
  Inf = bc_new_num (1,0);
  Inf->n_value[0] = 10;
}

void initR (Real* px)
{
  bc_init_num (px);
}

void initC (Complex* pz)
{
  bc_init_num (&pz->re);
  bc_init_num (&pz->im);
}

void copyR (Real* dst, Real src)
{
  *dst = bc_copy_num (src);
}

void copyC (Complex* dst, Complex src)
{
  dst->re = bc_copy_num (src.re);
  dst->im = bc_copy_num (src.im);
}

void str2R (const char *q, char **endptr, int iscale,
	    const struct numfmt* pnf, Real* pr)
{
  bc_a2num (pr, q, endptr, iscale, pnf);
}

void str2C (const char *q, char **endptr, int iscale,
	    const struct numfmt* pnf, Complex* pc)
{
  char *ptr, *ptr2;

  bc_init_num (&pc->re);
  bc_init_num (&pc->im);
  bc_a2num (&pc->re, q, &ptr, iscale, pnf);
  if ((endptr))
    *endptr = ptr;
  if (ptr != q)
    {
      if (*ptr == pnf->iu)
	{
	  pc->im = bc_copy_num (pc->re);
	  pc->re = bc_copy_num (_zero_);
	  if ((endptr))
	    *endptr = ptr + 1;
	}
      else
	{
	  while (is_space (*ptr))
	    ptr++;
	  if (*ptr == POS_SIGN || *ptr == NEG_SIGN)
	    {
	      bc_a2num (&pc->im, ptr, &ptr2, iscale, pnf);
	      if (*ptr2 != pnf->iu)
		pc->im = bc_copy_num (_zero_);
	      else
		{
		  if ((endptr))
		    *endptr = ptr2 + 1;
		}
	    }
	  /*
	     else
	     : we have successfully read a real number
	     but there is no another number after it.
	     : So, we leave pc->im set to zero.
	  */
	}
    }
  else
    /* 
       : q does not contain any valid number
       ==> pc->re is 0. Then we set pc->im to 0.
       We remark that, if endptr is
       not NULL, then *endptr == q.
    */
    pc->im = bc_copy_num (_zero_);
}

void add (Real s, Real t, Real* q, int iscale)
{
  bc_add (s, t, q, iscale);
}

void square (Real s, Real* q, int iscale)
{
  bc_multiply (s, s, q, iscale);
}

void divide (Real s, Real t, Real* q, int iscale)
{
  bc_divide (s, t, q, iscale);
}

void divide_by_int (Real* q, int d, int iscale)
{
  Real num, div;

  if (d == 0)
    {
      bc_free_num (q);
      bc_init_num (q);
    }
  else
    {
      bc_init_num (&div);
      bc_int2num (&div, d);
      num = bc_copy_num (*q);
      bc_divide (num, div, q, iscale);
      bc_free_num (&num);
      bc_free_num (&div);
    }
}

void square_root (Real* q, int iscale)
{
  bc_sqrt (q, iscale);
}

void Cabs (Complex z, Real* pm, int iscale)
{
  bc_num a, b, q;

  if ( z.re->n_sign == PLUS )
    a = bc_copy_num (z.re);
  else
    {
      bc_init_num (&a);
      bc_sub (Zero, z.re, &a, iscale); 
    }
  if ( (bc_is_zero(z.im)) )
    {
      *pm = a;
      return;
    }
  else if ( z.im->n_sign == PLUS )
    b = bc_copy_num (z.im);
  else
    {
      bc_init_num (&b);
      bc_sub (Zero, z.im, &b, iscale); 
    }
  bc_init_num (&q);
  if ( bc_compare (b, a) > 0 )
    {
      bc_divide (a, b, &q, iscale);
      bc_multiply (q, q, &a, iscale);
      bc_add (a, _one_, &q, iscale);
      bc_sqrt (&q, iscale);
      bc_multiply (b, q, pm, iscale);
    }
  else
    {
      /* a >= b ===> a > 0 */
      bc_divide (b, a, &q, iscale);
      bc_multiply (q, q, &b, iscale);
      bc_add (b, _one_, &q, iscale);
      bc_sqrt (&q, iscale);
      bc_multiply (a, q, pm, iscale);
    }
  bc_free_num (&q);
  bc_free_num (&b);
  bc_free_num (&a);
}

void Csub (Complex z1, Complex z2, Complex* pw, int iscale)
{
  bc_sub (z1.re, z2.re, &pw->re, iscale);
  bc_sub (z1.im, z2.im, &pw->im, iscale);
}

int cmp (Real p, Real q)
{
  return bc_compare (p, q);
}

int is0 (Real u)
{
  return bc_is_zero (u);
}

int smart_cmp (const Complex* pz1, const Complex* pz2, int flag)
{
  if (flag == 0)
    return 1;
  else if (flag > 0)
    return (bc_compare (pz1->re, pz2->re) >= 0 &&
	    bc_compare (pz1->im, pz2->im) >= 0);
  else /* flag < 0 */
    return (bc_compare (pz1->re, pz2->re) <= 0 &&
	    bc_compare (pz1->im, pz2->im) <= 0);
}

void printno (Real u, int m)
{
  bc_print_num (u, out_char, m);
}

#ifdef _MPA_DEBUG
void debug_printno (Real u, int m)
{
  bc_print_num (u, out_char_stderr, m);
}
#endif

void delR (Real* px)
{
  bc_free_num (px);
}

void delC (Complex* pz)
{
  bc_free_num (&pz->re);
  bc_free_num (&pz->im);
}

void end_mpa(void)
{
  bc_free_num (&Inf);
  bc_free_num (&Zero);
  bc_end();
}

#endif /* not USE_GMP */
