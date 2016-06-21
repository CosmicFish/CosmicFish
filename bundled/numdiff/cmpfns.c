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
#include<ctype.h>
#include"numdiff.h"
#include"linesplit.h"

/* See io.c */
extern char* read_line (FILE* pf, int* errcode); 
extern void print_lines (const char* line1, const char* line2,
			 unsigned long lineno1, unsigned long lineno2, 
			 int delimiter_only);
extern void print_fields (const char* field1, const char* field2,
			  size_t l1, size_t l2, 
			  unsigned long lineno1, unsigned long lineno2, 
			  unsigned long fieldno1, unsigned long fieldno2);
extern void print_errors (Real abserr, Real relerr);
extern void print_separator (void);

static
void field2cx (const char* field, size_t length, char** tail, 
	       int iscale, const struct numfmt* pnf, Complex* pz)
{
  char* ptr = (char*) field;
  char ch;

  ch = ptr[length];
  ptr[length] = '\0';
  str2C (ptr, tail, iscale, pnf, pz);
  ptr[length] = ch;
}

/*
  Be careful ! strNcmp() and strNcasecmp() are only used
  when n <= strlen (s) == strlen (t).
*/

static int strNcmp (const char* s, const char* t, size_t n)
{
  const char *p, *q;

  for (p = s, q = t; p < s + n && *p == *q; p++, q++);
  return p < s + n;
}

static int strNcasecmp (const char* s, const char* t, size_t n)
{
  const char *p, *q;

  for (p = s, q = t; p < s + n && TOLOWER(*p) == TOLOWER(*q); p++, q++);
  return p < s + n;
}

static
int cmp_fields (const char* field1, const char* field2,
		unsigned long fieldno1, unsigned long fieldno2, 
		size_t l1, size_t l2, argslist* argl,
		Real* abserr, Real* relerr)
{
  char *tail1, *tail2;
  Complex z1, z2;
  
  field2cx (field1, l1, &tail1, argl->iscale, &argl->nf1, &z1);
  field2cx (field2, l2, &tail2, argl->iscale, &argl->nf2, &z2);

#ifdef __DEBUG__
  fprintf (stderr, "l1 = %zu, tail1 - field1 = %zu\n", l1, tail1 - field1);
  fprintf (stderr, "l2 = %zu, tail2 - field2 = %zu\n", l2, tail2 - field2);
#endif  

  if (tail1 - field1 == l1 && tail2 - field2 == l2)
    {
      /*
	This second test manages the options -P
	and -N . If neither of them has been set,
	then the condition is always TRUE.
      */
      if ( (smart_cmp (&z2, &z1, argl->flag)) )
	{
	  int exit_code;

	  /* Numeric comparison */
	  Complex w;
	  Real x1, x2;
	  int iscale = argl->iscale;

	  initC (&w);
	  initR (&x1);
	  initR (&x2);
	  Csub (z1, z2, &w, iscale);
	  Cabs (w, abserr, iscale);
#ifdef _MPA_DEBUG
	  fputs ("*** MPA Debug output\n", stderr);
	  fputs ("1st number= ( ", stderr);
	  debug_printno (z1.re, 20);
	  fputs (", ", stderr);
	  debug_printno (z1.im, 20);
	  fputs (" )\n", stderr);
	  fputs ("2nd number= ( ", stderr);
	  debug_printno (z2.re, 20);
	  fputs (", ", stderr);
	  debug_printno (z2.im, 20);
	  fputs (" )\n", stderr);
	  fputs ("abs. err= ", stderr);
	  debug_printno (*abserr, 20);
	  fputs ("\n***     ***     ***\n", stderr);
#endif
	  if (argl->relerr_formula == CLASSIC_FORMULA)
	    {
	      Cabs (z1, &x1, iscale);
	      Cabs (z2, &x2, iscale);
	      if ( cmp (x1, x2) > 0 )
		copyR (&x1, x2);
	      if ( (is0(x1)) )
		(is0(*abserr)) ? copyR(relerr, Zero) : copyR(relerr, Inf);
	      else
		divide (*abserr, x1, relerr, iscale);
	    }
	  else if (argl->relerr_formula == WR_TO_FIRST_FILE)
	    {
	      Cabs (z1, &x1, iscale);
	      if ( (is0(x1)) )
		(is0(*abserr)) ? copyR(relerr, Zero) : copyR(relerr, Inf);
	      else
		divide (*abserr, x1, relerr, iscale);
	    }
	  else if (argl->relerr_formula == WR_TO_SECOND_FILE)
	    {
	      Cabs (z2, &x2, iscale);
	      if ( (is0(x2)) )
		(is0(*abserr)) ? copyR(relerr, Zero) : copyR(relerr, Inf);
	      else
		divide (*abserr, x2, relerr, iscale);
	    }
	  delR (&x2);
	  delR (&x1);
	  delC (&w);
	  delC (&z2);
	  delC (&z1);

	  if (!(argl->optmask & _2_MASK))
	    exit_code = 
	      thrlist_cmp (*relerr, argl->maxrelerr, fieldno1, fieldno2) > 0 && 
	      thrlist_cmp (*abserr, argl->maxabserr, fieldno1, fieldno2) > 0 ? 2:0;
	  else
	    exit_code =
	      thrlist_cmp (*relerr, argl->maxrelerr, fieldno1, fieldno2) > 0 || 
	      thrlist_cmp (*abserr, argl->maxabserr, fieldno1, fieldno2) > 0 ? 2:0;
	  if ((argl->optmask & _SS_MASK))
	    {
	      argl->Nentries++;
	      initR (&x1);
	      initR (&x2);
	      /* To compute the 1-norm of all errors */
	      add (*abserr, argl->N1abserr, &x1, iscale);
	      copyR (&argl->N1abserr, x1);
	      /* To compute the 2-norm of all errors */
	      square (*abserr, &x1, iscale);
	      add (x1, argl->N2abserr, &x2, iscale);
	      copyR (&argl->N2abserr, x2);
	      if ((exit_code))
		{
		  int test;

		  argl->Ndisperr++;
		  /* To compute the 1-norm of the displayed errors */
		  add (*abserr, argl->N1disperr, &x1, iscale);
		  copyR (&argl->N1disperr, x1);
		  /* To compute the 2-norm of the displayed errors */
		  square (*abserr, &x1, iscale);
		  add (x1, argl->N2disperr, &x2, iscale);
		  copyR (&argl->N2disperr, x2);
		  if ((test = cmp (*abserr, argl->Labserr)) > 0)
		    {
		      copyR (&argl->Labserr, *abserr);
		      copyR (&argl->Crelerr, *relerr);
		    }
		  else if (test == 0 && cmp (*relerr, argl->Crelerr) > 0)
		    copyR (&argl->Crelerr, *relerr);
		  if ((test = cmp (*relerr, argl->Lrelerr)) > 0)
		    {
		      copyR (&argl->Lrelerr, *relerr);
		      copyR (&argl->Cabserr, *abserr);
		    }
		  else if (test == 0 && cmp (*abserr, argl->Cabserr) > 0)
		    copyR (&argl->Cabserr, *abserr);
		}
	      delR (&x2);
	      delR (&x1);
	    }
	  return exit_code;
	}
      else
	return 0;
    }
  else
    {
      delC (&z2);
      delC (&z1);
      /* Byte by byte comparison */
      if ((argl->optmask & _SI_MASK))
	return (l1 != l2 || strNcasecmp (field1, field2, l1) != 0 ? 1 : 0);
      else
	return (l1 != l2 || strNcmp (field1, field2, l1) != 0 ? 1 : 0);
    }
}

static
int cmp_lines (const char* line1, const char* line2, 
	       unsigned long lineno1, unsigned long lineno2,  
	       const char** ifs1, const char** ifs2, int output_mode, argslist* argl)
{
  const unsigned long fieldno_upper_limit = 8*FIELDMASK_SIZE;

  if (!line1 && !line2)
    return 0;
  else if (!line1)
    {
      if (output_mode >= OUTMODE_NORMAL)
	print_lines (line1, line2, lineno1, lineno2, 0);
      return 1;
    }
  else if (!line2)
    {
      if (output_mode >= OUTMODE_NORMAL)
	print_lines (line1, line2, lineno1, lineno2, 0);
      return 1;
    }
  else
    {
      const char *field1, *field2;
      char *end_field1, *end_field2;
      size_t l1, l2;
      unsigned long fieldno1, fieldno2;
      Real abserr, relerr;
      int rv, lines_differ = 0, _1sttime = 1;

      initR (&abserr);
      initR (&relerr);
      field1 = string_spn (line1, ifs1, '\0');
      field2 = string_spn (line2, ifs2, '\0');
      fieldno1 = fieldno2 = 0;

      while (*field1 != '\0' && *field2 != '\0')
	{
	  /*
	    Ignore the fields selected through the option -X 1:
	   */
	  while ( *field1 != '\0' && fieldno1 < fieldno_upper_limit && 
		  (argl->ghostmask1[fieldno1 >> 3] & 0x80 >> (fieldno1 & 0x7)) )
	    {
              end_field1 = string_cspn (field1, ifs1, '\0');
	      field1 = string_spn (end_field1, ifs1, '\0');
	      fieldno1++; 
	    }
	  if ( fieldno1 >= fieldno_upper_limit )
	    {
	      printf (_("@ Line %lu in file \"%s\"\n  contains too many fields to be properly processed!\n"),
		      lineno1, argl->file1);
	      exit (EXIT_TROUBLE);
	    }
	  /*
	    Ignore the fields selected through the option -X 2:
	   */
	  while ( *field2 != '\0' && fieldno2 < fieldno_upper_limit && 
		  (argl->ghostmask2[fieldno2 >> 3] & 0x80 >> (fieldno2 & 0x7)) )
	    {
	      end_field2 = string_cspn (field2, ifs2, '\0');
	      field2 = string_spn (end_field2, ifs2, '\0');
	      fieldno2++; 
	    }
	  if ( fieldno2 >= fieldno_upper_limit )
	    {
	      printf (_("@ Line %lu in file \"%s\"\n  contains too many fields to be properly processed!\n"), 
		      lineno2, argl->file2);
	      exit (EXIT_TROUBLE);
	    }
	  if (*field1 != '\0' && *field2 != '\0')
	    {
              end_field1 = string_cspn (field1, ifs1, '\0');
	      end_field2 = string_cspn (field2, ifs2, '\0');
	      l1 = end_field1 - field1;
	      l2 = end_field2 - field2;
	      rv = cmp_fields (field1, field2, fieldno1, fieldno2, l1, l2, argl, &abserr, &relerr);
	      if ( rv >= 1 )
		{
		  if ( (output_mode > OUTMODE_QUIET) && 
		       !(rv == 1 && argl->optmask & _SE_MASK) &&
		       !(rv == 2 && argl->optmask & _SU_MASK))
		    {
		      if ((_1sttime))
			{
			  print_lines (line1, line2, lineno1, lineno2,  
				       output_mode != OUTMODE_VERBOSE &&
				       output_mode != OUTMODE_COINCISE);
			  _1sttime = 0;
			}
		      print_fields (field1, field2, l1, l2, lineno1, lineno2, fieldno1, fieldno2);
		      if ( rv == 2 )
			print_errors (abserr, relerr);
		      else
			print_separator ();
		    }
		  lines_differ = 1;
		}
	      field1 = string_spn (end_field1, ifs1, '\0');
	      fieldno1++; 
	      field2 = string_spn (end_field2, ifs2, '\0');
	      fieldno2++; 
	    }
	} /* end  while (*field1 != '\0' && *field2 != '\0') */

      delR (&abserr);
      delR (&relerr);

      /*
	Ignore the fields selected through the option -X 1:
      */
      while ( *field1 != '\0' && fieldno1 < fieldno_upper_limit && 
	      (argl->ghostmask1[fieldno1 >> 3] & 0x80 >> (fieldno1 & 0x7)) )
	{
	  end_field1 = string_cspn (field1, ifs1, '\0');
	  field1 = string_spn (end_field1, ifs1, '\0');
	  fieldno1++; 
	}
      if ( fieldno1 >= fieldno_upper_limit )
	{
	  printf (_("@ Line %lu in file \"%s\"\n  contains too many fields to be properly processed!\n"),
		  lineno1, argl->file1);
	  exit (EXIT_TROUBLE);
	}
      /*
	Ignore the fields selected through the option -X 2:
      */
      while ( *field2 != '\0' && fieldno2 < fieldno_upper_limit && 
	      (argl->ghostmask2[fieldno2 >> 3] & 0x80 >> (fieldno2 & 0x7)) )
	{
	  end_field2 = string_cspn (field2, ifs2, '\0');
	  field2 = string_spn (end_field2, ifs2, '\0');
	  fieldno2++; 
	}
      if ( fieldno2 >= fieldno_upper_limit )
	{
	  printf (_("@ Line %lu in file \"%s\"\n  contains too many fields to be properly processed!\n"),
		  lineno2, argl->file2);
	  exit (EXIT_TROUBLE);
	}
      
      if (*field1 != '\0')
	{
	  if (output_mode >= OUTMODE_NORMAL)
	    {
	      if ((_1sttime))
		{
		  print_lines (line1, line2, lineno1, lineno2, 
			       output_mode < OUTMODE_VERBOSE);
		  _1sttime = 0;
		}
	      print_fields (field1, field2, 0, 0, lineno1, lineno2, fieldno1, fieldno2);
	      printf (_("@ Line %lu in file \"%s\" is shorter than expected!\n"), 
		      lineno2, argl->file2);
	    }
	  return 1;
	}
      else if (*field2 != '\0')
	{
	  if (output_mode >= OUTMODE_NORMAL)
	    {
	      if ((_1sttime))
		{
		  print_lines (line1, line2, lineno1, lineno2,
			       output_mode < OUTMODE_VERBOSE);
		  _1sttime = 0;
		}
	      print_fields (field1, field2, 0, 0, lineno1, lineno2, fieldno1, fieldno2);
	      printf (_("@ Line %lu in file \"%s\" is shorter than expected!\n"), 
		      lineno1, argl->file1);
	    }
	  return 1;
	}
      else
	return lines_differ;
    }
}

extern char** def_ifs;

int cmp_files (FILE* pf1, FILE* pf2, argslist* argl)
{
  char *line1, *line2, **ifs1, **ifs2;
  int err1, err2, files_differ = 0;
  unsigned long lineno1 = 1, lineno2 = 1;
  size_t n;
  unsigned short i;
  unsigned char byte, elem;
  flg_array table = copy_of_intflagtab();

  ifs1 = (!argl->ifs1) ? def_ifs : argl->ifs1;
  ifs2 = (!argl->ifs2) ? def_ifs : argl->ifs2;
  if ( (table.ptr) )
    {
      /*
	Filter on
       */
      for (n = 0, err1 = err2 = OK; n < table.len; n += 4U)
	{
	  byte = table.ptr[n/4U];
	  for (i = 0; i < 4U; i++)
	    {
	      elem = (byte & 0x03 << 2*i) >> 2 * i;
	      switch (elem)
		{
		case 1:
		  line1 = read_line (pf1, &err1);
		  line2 = NULL;
		  break;
		case 2:
		  line1 = NULL;
		  line2 = read_line (pf2, &err2);
		  break;
		case 3:
		  line1 = read_line (pf1, &err1);
		  line2 = read_line (pf2, &err2);
		  break;
		case 0:
		  goto catch_error;
		}
	      if ( (cmp_lines (line1, line2, 
                               lineno1, lineno2, 
                               (const char**) ifs1, (const char**) ifs2, 
                               argl->output_mode, argl)) )
		{
		  files_differ = 1;
		  if (argl->output_mode == OUTMODE_OVERVIEW)
		    print_1overview_line (line1, 1, line2);
		}
	      else
		{
		  if (argl->output_mode == OUTMODE_OVERVIEW && !suppress_common_lines)
		    print_1overview_line (line1, 0, line2);
		}
	      if ((line1))
		{
		  lineno1++;
		  free ((void*)line1);
		}
	      if ((line2))
		{
		  lineno2++;
		  free ((void*)line2);
		}
	      if (err1 == OUT_OF_MEMORY || 
                  err2 == OUT_OF_MEMORY || 
                  err1 == FILE_IS_BINARY || 
                  err2 == FILE_IS_BINARY || 
                  err1 == READING_ERROR || 
                  err2 == READING_ERROR)
		goto catch_error;
	    }
	}
    } /* end of `if ( (table.ptr) )'*/
  else
    {
      /*
	Filter off
       */
      do
	{
	  line1 = read_line (pf1, &err1);
	  line2 = read_line (pf2, &err2);
	  if ( (cmp_lines (line1, line2, 
                           lineno1, lineno2, 
                           (const char**) ifs1, (const char**) ifs2, 
                           argl->output_mode, argl)) )
	    {
	      files_differ = 1;
	      if (argl->output_mode == OUTMODE_OVERVIEW)
		print_1overview_line (line1, 1, line2);
	    }
	  else
	    {
	      if (argl->output_mode == OUTMODE_OVERVIEW && !suppress_common_lines)
		print_1overview_line (line1, 0, line2);
	    }
	  if ((line1))
	    free ((void*)line1);
	  if ((line2))
	    free ((void*)line2);
	  lineno1++, lineno2++;
	} while (err1 == OK && err2 == OK);
    }

 catch_error:
  /* 
   * If we arrive here, then
   *
   * either ist n == table.len,
   * or elem == 0,
   * or either err1 or err2 is different from OK.
   */
  fflush (stdout);
  if ( (table.ptr) )
    {
      if (n >= table.len || !elem)
	return files_differ;
    }
  if (err1 == OK)
    {
      switch (err2)
	{
        case FILE_IS_BINARY:
          fprintf (stderr, _("\n***  File \"%s\" is binary,\n***  cannot read from it\n"), argl->file2);
	  return FILE_IS_BINARY;          
	case READING_ERROR:
	  fprintf (stderr, _("\n***  Error while reading from file \"%s\"\n"), argl->file2);
	  return READING_ERROR;
	case OUT_OF_MEMORY:
	  fprintf (stderr, _("\n***  Out of memory while reading from file \"%s\"\n"), argl->file2);
	  return OUT_OF_MEMORY;
	case LINE_INTERR:
	  if (getc (pf1) == EOF)
	    break;
	  /* default == EOF_REACHED */
	default: 
	  fprintf (stderr, _("\n***  End of file \"%s\" reached while trying to read line %lu.\n"), argl->file2, lineno2-1);
	  fprintf (stderr, _("***  File \"%s\" has more lines than file \"%s\",\n"),
		   argl->file1, argl->file2);
	  fprintf (stderr, _("***  line %lu is the last one read from file \"%s\"\n\n"), lineno1-1, argl->file1);
	}
      return files_differ;
    }
  else if (err1 == LINE_INTERR)
    {
      switch (err2)
	{
        case FILE_IS_BINARY:
          fprintf (stderr, _("\n***  File \"%s\" is binary,\n***  cannot read from it\n"), argl->file2);
	  return FILE_IS_BINARY;          
	case READING_ERROR:
	  fprintf (stderr, _("\n***  Error while reading from file \"%s\"\n"), argl->file2);
	  return READING_ERROR;
	case OUT_OF_MEMORY:
	  fprintf (stderr, _("\n***  Out of memory while reading from file \"%s\"\n"), argl->file2);
	  return OUT_OF_MEMORY;
	case OK:
	  if (getc (pf2) != EOF)
	    {
	      fprintf (stderr, _("\n***  End of file \"%s\" reached while trying to read line %lu.\n"), argl->file1, lineno1-1);
	      fprintf (stderr, _("***  File \"%s\" has more lines than file \"%s\",\n"),
		       argl->file2, argl->file1);
	      fprintf (stderr, _("***  line %lu is the last one read from file \"%s\"\n\n"), lineno2-1, argl->file2);
	    }
	  break;
	case EOF_REACHED:
	  fprintf (stderr, _("\n***  End of file \"%s\" reached while trying to read line %lu.\n"), argl->file2, lineno2-1);
	  fprintf (stderr, _("***  File \"%s\" has more lines than file \"%s\",\n"),
		   argl->file1, argl->file2);
	  fprintf (stderr, _("***  line %lu is the last one read from file \"%s\"\n\n"), lineno1-1, argl->file1);
	  /*
	    No particular action to do if err2 == LINE_INTERR
	  */
	}
      return files_differ;
    }
  else if (err1 == EOF_REACHED)
    {
      switch (err2)
	{
        case FILE_IS_BINARY:
          fprintf (stderr, _("\n***  File \"%s\" is binary,\n***  cannot read from it\n"), argl->file2);
	  return FILE_IS_BINARY;          
	case READING_ERROR:
	  fprintf (stderr, _("\n***  Error while reading from file \"%s\"\n"), argl->file2);
	  return READING_ERROR;
	case OUT_OF_MEMORY:
	  fprintf (stderr, _("\n***  Out of memory while reading from file \"%s\"\n"), argl->file2);
	  return OUT_OF_MEMORY;
	case OK:
	case LINE_INTERR:
	  fprintf (stderr, _("\n***  End of file \"%s\" reached while trying to read line %lu.\n"), argl->file1, lineno1-1);
	  fprintf (stderr, _("***  File \"%s\" has more lines than file \"%s\",\n"),
		   argl->file2, argl->file1);
	  fprintf (stderr, _("***  line %lu is the last one read from file \"%s\"\n\n"), lineno2-1, argl->file2);
	}
      return files_differ;
    }
  else if (err1 == FILE_IS_BINARY)
    {
      fprintf (stderr, _("\n***  File \"%s\" is binary,\n***  cannot read from it\n"), argl->file1);
      return FILE_IS_BINARY;          
    }
  else if (err1 == READING_ERROR)
    {
      fprintf (stderr, _("\n***  Error while reading from file \"%s\"\n"), argl->file1);
      return READING_ERROR;
    }
  else /* err1 == OUT_OF_MEMORY */
    {
      fprintf (stderr, _("\n***  Out of memory while reading from file \"%s\"\n"), argl->file1);
      return OUT_OF_MEMORY;
    }
}
