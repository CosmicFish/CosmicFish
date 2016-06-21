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
#include"linesplit.h"

#define O_DEF_SEP '\n'
#define ESC_CHAR  '\\'

const unsigned char InvDigit = (unsigned char) -1;

static
unsigned char is_hex_digit (char ch)
{
  if (ch >= '0' && ch <= '9')
    return ch - '0';
  else
    {
      switch (ch)
	{
	case 'a':
	case 'A':
	  return 10;
	case 'b':
	case 'B':
	  return 11;	  
	case 'c':
	case 'C':
	  return 12;
	case 'd':
	case 'D':
	  return 13;
	case 'e':
	case 'E':
	  return 14;	  
	case 'f':
	case 'F':
	  return 15;
	default:
	  return InvDigit;
	}
    }
}

static
unsigned char is_oct_digit (char ch)
{
  if (ch >= '0' && ch <= '7')
    return ch - '0';
  else
    return InvDigit;
}

static
int is_hex (const char* pafter_escape, char* byte)
{
  unsigned char hb, lb;

  /*
    The return value is the number of hexadecimal digits successfully read
  */  
  if ( (hb = is_hex_digit (*pafter_escape)) != InvDigit )
    {
      int rv = 1;

      lb = is_hex_digit (*(pafter_escape+1));
      if (lb == InvDigit)
	{
	  lb = hb;
	  hb = 0;
	}
      else
	rv = 2;
      
      if ((byte))
	*byte = (char)(16 * hb + lb);
      return rv;
    }
  else
    {
      if ((byte))
	*byte = -1;
      return 0;
    }
}

static
int is_oct (const char* pafter_escape, char* byte)
{
  unsigned char ho, mo, lo;

  /*
    The return value is the number of octal digits successfully read
  */
  if ( (ho = is_oct_digit (*pafter_escape)) != InvDigit )
    {
      int rv = 1;

      mo = is_oct_digit (*(pafter_escape+1));
      lo = is_oct_digit (*(pafter_escape+2));
      if (mo == InvDigit)
	{
	  lo = ho;
	  ho = 0;
	  mo = 0;
	}
      else
	{
	  rv++;
	  if (lo == InvDigit)
	    {
	      lo = mo;
	      mo = ho;
	      ho = 0;
	    }
	  else
	    rv++;
	}
      if ((byte))
	*byte = (char)(64 * ho + 8 * mo + lo);
      return (64 * (int)ho + 8 * (int)mo + (int)lo < 256 ? rv : 0);
    }
  else
    {
      if ((byte))
	*byte = -1;
      return 0;
    }
}

static
char process_character (const char* chp, const char** new_chp)
{
  const char* nchp;
  char byte;
  int r;

  if (*chp == ESC_CHAR)
    {
      switch (*(chp+1))
	{
	case 'a':
	  byte = '\a';
	  nchp = chp+2;
	  break;
	case 'b':
	  byte = '\b';
	  nchp = chp+2;
	  break;
	case 'f':
	  byte = '\f';
	  nchp = chp+2;
	  break;
	case 'n':
	  byte = '\n';
	  nchp = chp+2;
	  break;
	case 'r':
	  byte = '\r';
	  nchp = chp+2;
	  break;
	case 't':
	  byte = '\t';
	  nchp = chp+2;
	  break;
	case 'v':
	  byte = '\v';
	  nchp = chp+2;
	  break;
	case 's':
	  byte = ' ';
	  nchp = chp+2;
	  break;	    
	case 'x':
	  if ( (r = is_hex (chp+2, &byte)) && byte != '\0' )
	    nchp = chp + (r + 2);
	  else
	    {
	      byte = *(chp+1);
	      nchp = chp+2;
	    }
	  break;
	case '0':
	case '1':
	case '2':
	case '3':
	case '4':
	case '5':
	case '6':
	case '7':
	  if ( (r = is_oct (chp+1, &byte)) && byte != '\0' )
	    nchp = chp + (r + 1);
	  else
	    {
	      byte = *(chp+1);
	      nchp = chp+2;
	    }
	  break;
	default:
	  byte = *(chp+1);
	  nchp = chp+2;
	  break;
	}
    }
  else
    {
      byte = *chp;
      nchp = chp+1;
    }
  if ((new_chp))
    *new_chp = nchp;
  return byte;
}

/*
  Remark: process the substring [BPTR, EPTR) and return
          the result (NULL in case of error while allocating memory for the result). 
          Precondition is EPTR >= BPTR.
*/
static
char* process_substring (const char* bptr, const char* eptr)
{
  size_t subssize = eptr - bptr + 1;
  const char* ptr;
  char *pstr, *pstrp;
    
  if ( !(pstr = (char*) calloc(subssize, sizeof(char))) )
    return NULL;
  else
    {
      for (pstrp = pstr, ptr = bptr; ptr < eptr; pstrp++)
	{
	  *pstrp = process_character (ptr, &ptr);
	}
      return pstr;
    }
}

/*
  Create and return a vector of strings using the description
  contained in the string pointed to by STR.
  The items in STR (to each of which should correspond
  a string in the returned vector) are separated by the
  character SEPARATOR.
  In case of error while allocating memory for the vector
  and its elements return NULL.
  Return NULL also if STR == NULL.

  Remark: SEPARATOR cannot be the nul character.
          Return NULL if SEPARATOR is the nul character.  
*/
char** ssplit (const char* str, char separator)
{
  size_t i, n;
  const char *beg, *ptr, *ptr2sep;
  char** sv;

  if (!str || separator == '\0')
    return NULL;
  for (beg = str; *beg == separator; beg++);
  /* 
     Now BEG points to the first charatacrer of the buffer 
     pointed to by STR which is not equal to SEPARATOR.
  */

  /* 
     First count the substrings contained 
     in the buffer pointed to by STR.
  */
  for (n = 1, ptr = beg; (ptr2sep = strchr (ptr, separator)) != NULL;
       n++)
    {
      for (ptr = ptr2sep+1; *ptr == separator; ptr++);
    }
  /*
    Now allocate memory for a vector of N+1 char*.
    If the allocation fails, return NULL.
  */
  if ( !(sv = (char**) malloc ((n+1)*sizeof(char*))) )
    return NULL;
  sv[n] = NULL;
  
  for (i = 0, ptr = beg; (ptr2sep = strchr (ptr, separator)) != NULL;
       i++)
    {
      sv[i] = process_substring (ptr, ptr2sep);
      if (!sv[i])
        {
          delete_string_vector (sv);
          return NULL;
        }
      for (ptr = ptr2sep+1; *ptr == separator; ptr++);
    }
  if (*ptr != '\0')
    {
      ptr2sep = strchr (ptr, '\0');
      sv[i] = process_substring (ptr, ptr2sep);
      if (!sv[i])
        {
          delete_string_vector (sv);
          return NULL;
        }
    }
  return sv;
}

/*
  Create and return a vector of strings using the characters
  contained in the string pointed to by STR. To each
  (eventually escaped) character in this string will 
  correspond exactly one string in the returned vector.
  Return NULL if STR == NULL or in case of out of memory.
*/
char** ssplit_former_way (const char* str)
{
  if ((str))
    {
      size_t n, ls = strlen(str);
      char **sv;
      const char *ptr, *nptr;

      sv = (char**) calloc (ls + 1, sizeof(char*));
      if (!sv)
        return NULL;
      for (n = 0, ptr = str; *ptr != '\0'; ptr = nptr, n++)
        {
          sv[n] = (char*) malloc (2 *sizeof(char));
          if ((sv[n]))
            {
              sv[n][0] = process_character (ptr, &nptr);
              sv[n][1] = '\0';
            }
          else
            {
              delete_string_vector (sv);
              return NULL;
            }
        }
      return sv;
    }
  else
    return NULL;
}

/*
  Process the string pointed to by ISTR and return the result
  (NULL in case of error while allocating memory for the result).
*/
char* get_separating_string (const char* istr)
{
  return process_substring (istr, istr+strlen(istr));
}

/*
  Write to the file pointed to by FP the strings contained in
  the vector SV. Use SEPARATOR to separate each string from
  the following one.
*/
void print_string_vector (FILE* fp, const char** sv, char separator)
{
  size_t n;

  if (!sv)
    {
      fputs ("<Empty>", fp);
      fputc (separator, fp);
    }
  else
    {
      for (n = 0; sv[n] != NULL; n++)
	{
	  fprintf (fp, "\"%s\"%c", sv[n], separator);
	}
    }
}

/*
  Rearrange the strings of the vector SV in descending order
  with respect to their length.

  Rem.: Pre-condition is that SV is NULL-terminated.
        This function is suitable only for small vectors,
	since it uses a bubble-sort algorithm.
*/
void sort_string_vector (char** sv)
{
  if ((sv))
    {
      size_t n, m, l, lmax, poslmax;
      char *tmp;
      
      for (n = 0; sv[n] != NULL; n++)
	{
	  lmax = strlen(sv[n]);
	  poslmax = n;
	  for (m = n+1; sv[m] != NULL; m++)
	    {
	      if ( (l = strlen(sv[m])) > lmax )
		{
		  lmax = l;
		  poslmax = m;
		}
	    }
	  tmp = sv[n];
	  sv[n] = sv[poslmax];
	  sv[poslmax] = tmp;
	}
    }
}

/*
  Remove duplicates from the vector SV.
*/
void remove_duplicates_from_string_vector (char** sv)
{
  if ((sv))
    {
      size_t k, m, n;

      for (n = 0; sv[n] != NULL; n++)
	{
          m = n+1;
	  while (sv[m] != NULL)
	    {
              if (strcmp (sv[m], sv[n]) == 0)
                {
                  free((void*)sv[m]);
                  for (k = m+1; sv[k] != NULL; k++)
                    sv[k-1] = sv[k];
                  sv[k-1] = NULL;
                }
              else
                m++;
	    }
        }
    }
}

/*
  Return 1 if the string pointed to by STR is found in the vector SV,
  otherwise 0. 0 should be also returned if STR or SV is NULL.
*/
int is_string_in_vector (const char* str, const char** sv)
{
  if ((sv) && (str))
    {
      size_t n;

      for (n = 0; sv[n] != NULL && strcmp(str, sv[n]) != 0; n++);
      return (sv[n] == NULL ? 0 : 1);
    }
  else
    return 0;
}

/*
  Return 0 if there is no string in the vector SV which contains the
  character CH, otherwise return the length of the longest string
  between those ones which contain the character CH.
  0 should also be returned if SV is null.
*/
size_t is_char_in_vector (int ch, const char** sv)
{
  if ((sv))
    {
      size_t l, lm, n;

      for (lm = n = 0; sv[n] != NULL; n++)
        {
          if ( (strchr (sv[n], ch)) && (l = strlen (sv[n])) > lm)
            lm = l;
        }
      return lm;
    }
  else
    return 0;
}

/*
  Remove the memory allocated for the strings of the vector SV
  and then free the memory allocated for the vector itself.
*/
void delete_string_vector (char** sv)
{
  size_t n;

  if ((sv))
    {
      for (n = 0; sv[n] != NULL; n++)
	{
	  free ((void*)sv[n]);
	}
      free((void*)sv);
    }
}

/*
  Return a pointer to the position following the initial
  segment of STR that does not contain any string
  from the vector SV. If such an initial segment does not
  exist, return a pointer to STR.
  Consider the string STR as ending at the first occurrence of EOS.

  SV must be NULL terminated, it cannot contain the empty ("") string
  nor a string of length > 1 with EOS being one of its non-null characters
  (but SV may well contain the string of length 1 having EOS as its 
  only non-null character). 

  Rem.: EOS can be the null character.
        If the string pointed to by STR does not contain any EOS
        character, a buffer overrun will occur.
*/
char* string_cspn (const char* str, const char** sv, int eos)
{
  register const char *sviptr;
  register const char *endptr;
  register const char *nendptr;
  register const char *ptr;
  register size_t n;

  if (!str || !sv)
    {
      /* security check */
      return NULL;
    }
  else
    {
      for (endptr = str; *endptr != eos; endptr++);
      for (nendptr = str; nendptr < endptr; nendptr++)
        {
          for (n = 0; sv[n] != NULL; n++)
            {
              for (ptr = nendptr, sviptr = sv[n];
                   *sviptr != '\0' && *sviptr == *ptr; 
                   sviptr++, ptr++);                   
              if (*sviptr == '\0')
                return (char*)nendptr;
            }
        }
      return (char*)nendptr;
    }
}

/*
  Return a pointer to the position following the initial
  segment of STR that consists entirely of strings
  from the vector SV. If such an initial segment does not
  exist, return a pointer to STR.
  Consider the string STR as ending at the first occurrence of EOS.

  SV must be NULL terminated, it cannot contain the empty ("") string
  nor a string of length > 1 with EOS being one of its non-null characters
  (but SV may well contain the string of length 1 having EOS as its 
  only non-null character). 

  Rem.: this function works under the assumption that the strings in 
        the vector SV are ordered according to their lengths, where SV[0]
        is the string with the greatest length.

        EOS can be the null character.
*/
char* string_spn (const char* str, const char** sv, int eos)
{
  register const char *ptr;
  register const char *nptr;
  register const char *sviptr;
  register size_t n;

  if (!str || !sv)
    {
      /* security check */
      return NULL;
    }
  else
    {
      ptr = str;
      while (*ptr != eos)
	{
          /*
            Rem.: if strlen(sv[n])== 1 and sv[n][0] == EOS, then
                  strstr(ptr, sv[n]) != ptr. Thus, whenever the following
                  for cycle terminates, sv[n] can not be equal to the
                  string "<EOS>" (i.e. the string having EOS as its only
                  null character).
          */
	  for (n = 0; sv[n] != NULL; n++)
            {
              for (nptr = ptr, sviptr = sv[n]; 
                   *sviptr!='\0' && *sviptr == *nptr; 
                   sviptr++, nptr++);
              if (*sviptr == '\0')
                {
                  /*
                    Rem.: if sv[n] does not contain any EOS, then
                    by setting PTR to NPTR we do not
                    skip any EOS.
                  */
                  ptr = nptr;
                  break;
                }
            }
          if (!sv[n])
	    break;
	} 
      return (char*)ptr;
    }
}

#ifdef _TEST_LINE_SPLIT_

#define I_DEF_SEP ' '

static
void print_help (const char* progname)
{
  printf ("Usage: %s STRING\n\n", progname);
}

static
void print_substring (FILE* fp, const char* bptr, const char* eptr, int nl)
{
  const char *ptr;

  if (eptr > bptr)
    {
      for (ptr = bptr; ptr != eptr; ptr++)
	{
	  putc (*ptr, fp);
	}
      if ((nl))
	putc ('\n', fp);	
    }
}

#define BUFFSIZE 1024

int main (int argc, char* argv[])
{
  if (argc != 2)
    {
      print_help(argv[0]);
      return 1;
    }
  else
    {
      char** string_vector = NULL;
      char** sv = NULL;
      char linebuff[BUFFSIZE] = "";
      char *rv, *ptr, *endptr;
      size_t l;

      string_vector = ssplit (argv[1], I_DEF_SEP);
      sv = ssplit_former_way (argv[1]);
      remove_duplicates_from_string_vector (string_vector);
      remove_duplicates_from_string_vector (sv);
      sort_string_vector (string_vector);
      print_string_vector (stdout, (const char**)string_vector, O_DEF_SEP);
      l = is_char_in_vector (':', string_vector);
      printf ("Length of the longest string containing \':\' = %zu\n", l);
      puts ("\n\nSplitting the string in the former way produces the following result:");
      print_string_vector (stdout, (const char**)sv, O_DEF_SEP);
      do
	{
#ifdef _MINOR_TEST_
	  puts ("\nEnter a line of text (Ctrl+D to terminate)");
#endif
	  rv = fgets (linebuff, BUFFSIZE, stdin);
	  if ((rv))
	    {
#ifdef _MINOR_TEST_
	      ptr = string_cspn (linebuff, (const char**)string_vector, '\0');
	      fputs ("Cspn =", stdout);
	      print_substring (stdout, linebuff, ptr, 0);
	      fputs ("|EoS|\n", stdout);

	      ptr = string_spn (linebuff, (const char**)string_vector, '\0');
	      fputs ("Spn  =", stdout);
	      print_substring (stdout, linebuff, ptr, 0);	      
	      fputs ("|EoS|\n", stdout);
#else
	      unsigned long fieldno;


	      for (fieldno = 1, ptr = linebuff; *ptr != '\0'; fieldno++)
		{
		  ptr = string_spn (ptr, (const char**)string_vector, '\0');
		  endptr = string_cspn (ptr, (const char**)string_vector, '\0');
		  if ((*ptr))
		    {
		      printf ("%3lu.>", fieldno);
		      print_substring (stdout, ptr, endptr, 0);
		      puts ("<");
		    }
		  ptr = endptr;
		}
	      putchar ('\n');
#endif /* _MINOR_TEST_ */
	    }
	} while ((rv));
      delete_string_vector (string_vector);
      return 0;
    }
}

#endif /* _TEST_LINE_SPLIT_ */
