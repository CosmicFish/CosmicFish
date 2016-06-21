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

#ifdef _TEST_READ_LINE_

#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#define BUFF_SIZE 32

/* Error codes */
#define OK             0
#define LINE_INTERR    1
#define EOF_REACHED    2
#define READING_ERROR  3
#define OUT_OF_MEMORY  4
/* *** */
#define OPEN_FAILED    5
#define WRONG_USAGE    6
#define FILE_IS_BINARY 7

#endif /* _TEST_READ_LINE_ */

char* read_line (FILE* pf, int* errcode)
{
  char buffer[BUFF_SIZE];
  char *ptr, *line = NULL;
  size_t lline = 0;

  register size_t n;
  register int ch;
  int exception_occurred = 0;

  do
    {
      for (n = 0; n < BUFF_SIZE-1 && (ch = getc(pf)) != EOF && ch != '\0' && ch != '\n'; n++)
	buffer[n] = ch;
      if ( !(exception_occurred = ((ch == EOF && n == 0) || ch == '\0')) )
	{
	  /* 
	     If we enter this block, then CH is not the nul character. 
	     Thus either is n == BUFF_SIZE-1, or is CH == EOF or is CH == NEWLINE.
	   */
	  if (n == BUFF_SIZE-1 || ch == EOF)
	    {
	      buffer[n] = '\0';
	      lline += n;
	    }
	  else
	    {
	      buffer[n] = '\n';
	      buffer[n+1] = '\0';
	      lline += (n+1);	      
	    }
	  if (!line)
	    ptr = (char*) calloc (lline + 1, sizeof(char));
	  else
	    ptr = (char*) realloc ((void*)line, (lline + 1) * sizeof(char));
	  if (!ptr)
	    {
	      if ((line))
		free ((void*)line);
	      *errcode = OUT_OF_MEMORY;
	      return NULL;
	    }
	  else
	    {
	      line = ptr;
	      strcat (line, buffer);
	    }
	  if (lline > 0 && line[lline-1] == '\n')
	    break;
	}
    } while (!exception_occurred);
  if ((exception_occurred))
    {
      if ( (ferror(pf)) )
	*errcode = READING_ERROR;
      else if ( ch == '\0' )
	*errcode = FILE_IS_BINARY;
      else if (lline > 0)
	*errcode = LINE_INTERR;
      else
	*errcode = EOF_REACHED; 
    }
  else
    *errcode = OK;
  return line;
} 

#ifdef _TEST_READ_LINE_

int main (int argc, char* argv[])
{
  FILE* fp;

  if (argc != 2)
    {
      fprintf (stderr, "\n***  Usage: %s FILEPATH\n", argv[0]);
      return 1;
    }
  else if ( !(fp = fopen (argv[1], "r")) )
    {
      fprintf (stderr, "\n*** Cannot open \"%s\" for reading:\n", argv[1]);
      perror(NULL);
      return 1;
    }
  else
    {
      char* line = NULL;
      int errcode = OK;

      do
	{
	  line = read_line (fp, &errcode);
	  switch (errcode)
	    {
	    case FILE_IS_BINARY:
	      fputs ("\n***  Cannot read from a binary file\n", stderr);
	      break;	      
	    case LINE_INTERR:
	      fputs (line, stdout);
	      puts ("<LNINTERR>");
	      break;
	    case READING_ERROR:
	      fputs ("\n***  Error while reading from file\n", stderr);
	      break;
	    case OUT_OF_MEMORY:
	      fputs ("\n***  Memory exhausted while reading from file, abort now...\n", stderr);
	      break;
	    case EOF_REACHED:
	      puts ("<EOF>");
	      break;
	    default: /* OK */
	      fputs (line, stdout);
	      break;
	    }
	  if ((line))
	    free ((void*)line);
	} while (errcode == OK);
      if ( (fclose(fp)) )
	{
	  fprintf (stderr, "\n*** Could not close \"%s\":\n", argv[1]);
	  perror(NULL);
	  return 1;
	}
      else
	{
	  return 0;
	}
    }
}

#endif /* _TEST_READ_LINE_ */
