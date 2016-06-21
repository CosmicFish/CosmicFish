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
#include<ctype.h>
#include<errno.h>
#ifdef ENABLE_NLS
#include<locale.h>
#endif
#include"getopt.h"
#include"linesplit.h"
#include"numdiff.h" /* For NEWLINE and LINE_INTERR definitions */
#include"ndselect.h"

#include"read_line.c"

static 
int we_can_go_on (const char *field, 
                  unsigned long fieldno, unsigned long last_field)
{
  if (!last_field)
    return (*field != '\0');
  else
    return (fieldno <= last_field && *field != '\0');
}

static
void print_line (const char *line, const char **ifs, const char* osep,
		 unsigned long first_field, unsigned long last_field,
		 unsigned long increment, int omit_if_empty)
{
  unsigned long fieldno, lf;
  int print_nl, at_least_one_field_printed;
  const char *field, *endfield, *nextfield, *stop, *ptr;
  
  if (!last_field)
    lf = 0;
  else
    {
      lf = (last_field - first_field) / increment;
      lf = lf * increment + first_field;
    }
  at_least_one_field_printed = 0;
  print_nl = strchr (line, NEWLINE) != NULL ? 1 : 0;

  field = string_spn (line, ifs, '\0');
  for (fieldno = 1; (we_can_go_on (field, fieldno, last_field)) ; fieldno++)
    {
      endfield = string_cspn (field, ifs, '\0');
      nextfield = string_spn (endfield, ifs, '\0');
      if (fieldno >= first_field && (fieldno - first_field) % increment == 0)
	{
          if ((osep))
            stop = endfield;
          else
            stop = (fieldno == lf || *nextfield == '\0') ? endfield : nextfield;
	  for (ptr = field; ptr < stop; ptr++)
	    fputc(*ptr, stdout);
          if ((osep))
            {
              fputs (osep, stdout);
            }
	  at_least_one_field_printed = 1;
	}
      field = nextfield;
    }
  if ((omit_if_empty))
    {
      if ((print_nl) && (at_least_one_field_printed))
	fputc(NEWLINE, stdout);
    }
  else
    {
      if ((print_nl))
	fputc(NEWLINE, stdout);
      else
	fputc(' ', stdout);
    }
}

static
int scan_file (char** def_ifs, const Argslist* data)
{
  FILE *fp;
  char *line_buffer;
  char **ifs;
  unsigned long lineno;
  int errcode, omit_empty_lines;

  ifs = (!data->ifs) ? def_ifs : data->ifs;
  omit_empty_lines = (data->optmask & ___X_MASK) ? 1 : 0;
  if (!data->file || !*data->file)
    fp = stdin;
  else
    {
      if ( !(fp = fopen (data->file, "r")) )
	return OPEN_ERROR;
    }

  if (!data->end_line)
    {
      for (lineno = 1; 
	   (line_buffer = read_line(fp, &errcode), errcode <= LINE_INTERR); 
	   lineno++)
	{
	  if (lineno >= data->begin_line && (lineno - data->begin_line) % data->step == 0)
	    print_line (line_buffer, (const char**) ifs, data->osep, 
                        data->first_field, data->last_field, data->increment, omit_empty_lines);
	  free ((void*)line_buffer);
	}
    }
  else
    {
      for (lineno = 1; 
	   lineno <= data->end_line && 
	     (line_buffer = read_line(fp, &errcode), errcode <= LINE_INTERR); 
	   lineno++)
	{
	  if (lineno >= data->begin_line && (lineno - data->begin_line) % data->step == 0)
	    print_line (line_buffer, (const char**) ifs, data->osep,
                        data->first_field, data->last_field, data->increment, omit_empty_lines);
	  free ((void*)line_buffer);
	}
    }

  if (errcode <= EOF_REACHED)
    {
      if (fp == stdin)
	return OK;
      else
	return fclose (fp) == EOF ? CLOSE_ERROR : OK;
    }
  else
    {
      if((line_buffer))
	free ((void*)line_buffer);
      if (fp != stdin)
	fclose (fp);
      return READ_ERROR;
    }
}

static
void print_selversion (const char* progname)
{
  printf ("%s %s\n", progname, VERSION);
  printf ("Copyright (C) 2005, 2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013  %s <ivprimi@libero.it>\n", 
	  /* TRANSLATORS: This is a proper name.  See the gettext
	     manual, section Names.
	     Pronounciation is like "evaa-no pree-me".  */
	  _("Ivano Primi"));
  printf (_("\
License GPLv3+: GNU GPL version 3 or later,\n\
see <http://gnu.org/licenses/gpl.html>.\n\
This is free software: you are free to change and redistribute it.\n\
There is NO WARRANTY, to the extent permitted by law.\n"));
}

static
void print_selhelp (const char* progname)
{
  puts (_("Usage:"));
  printf ("%s -h|--help|-v|--version   %s\n\n", progname, _("or"));
  printf ("%s %s\n", progname, "[-b N][-e N][-s N][-F N][-L N][-I N][-S IFS][-D DELIMS][-O OSEP][-x][-l PATH][-o PATH] [FILE]");
  /* %%% */
  printf(_("\nPrint to standard output a subset of lines and fields from a given file.\n"));
  printf ("\n%s\n%s\n%s\n\n%s\n\n",
	  _("The argument after the options is the name of the file to read from."),
	  _("The complete path of the file should be given,\na directory name is not accepted."),
	  _("If no input file is specified, the program reads from the standard input."),
	  _("Exit status: 0 in case of normal termination, -1 (255) in case of error"));
  printf ("-b, --beginning, --start=N\n    %s\n    %s\n", 
	  _("Set to N the number of the first line to print"),
	  _("(The default behavior is to start with line number 1)"));
  printf ("-e, --end=N\n    %s\n    %s\n", 
	  _("Set to N the number of the last line that can be printed"),
	  _("(The default behavior is to arrive till to the end of the file)"));
  printf ("-s, --step=N\n    %s\n    %s\n", 
	  _("Set to N the increment to use when selecting the lines to print"),
	  _("(The default value for the increment is 1)"));
  printf ("-F, --first-field=N\n    %s\n    %s\n", 
	  _("Set to N the number of the first field to print"),
	  _("(The default behavior is to start with field number 1)"));
  printf ("-L, --last-field=N\n    %s\n    %s\n", 
	  _("Set to N the number of the last field that can be printed"),
	  _("(The default behavior is to arrive till to the end of every line)"));
  printf ("-I, --increment=N\n    %s\n    %s\n", 
	  _("Set to N the increment to use when selecting the fields to print"),
	  _("(The default value for the increment is 1)"));
  printf ("-S, --separators=IFS\n    %s\n    %s\n",
          _("Specify the set of characters to use as delimiters\n    while splitting the input lines into fields"),
          _("(The default set of delimiters is space, tab and newline)"));
  printf ("-D, --delimiters=DELIMS\n    %s\n    %s\n",
          _("Specify the set of strings to use as delimiters\n    while splitting the input lines into fields"),
          _("(The default set of delimiters is space, tab and newline)"));
  printf ("-O, --output-separator=OSEP\n    %s\n    %s\n",
          _("Specify the string to use as separator\n    while writing the selected fields to the standard output"),
          _("(The default behavior consists in reusing\n     the delimiters found in the input lines)"));
  printf ("-x, --omit-empty-lines\n    %s\n",
	  _("Do not print empty lines"));
  printf ("-l, --warnings-to=PATH\n    %s\n",
	  _("Redirect warning and error messages from stderr to the indicated file"));
  printf ("-o, --output=PATH\n    %s\n",
	  _("Redirect output from stdout to the indicated file"));
  printf ("-h, --help\n    %s\n", _("Show this help message"));
  printf ("-v, --version\n    %s\n\n", _("Show version number, Copyright, Distribution Terms and NO-Warranty"));
}

extern int errno;
extern char *optarg; 
extern int optind;

static
int set_args (int argc, char* argv[], Argslist *list)
{
  const char *optstring = "hb:e:s:F:L:I:S:D:O:xl:o:v";
  struct option long_options[] = {
    {"help",             0, NULL, 'h'},
    {"beginning",        1, NULL, 'b'},
    {"start",            1, NULL, 'b'},
    {"end",              1, NULL, 'e'},
    {"step",             1, NULL, 's'},
    {"first-field",      1, NULL, 'F'},
    {"last-field",       1, NULL, 'L'},
    {"increment",        1, NULL, 'I'},
    {"separators",       1, NULL, 'S'},
    {"delimiters",       1, NULL, 'D'},
    {"output-separator", 1, NULL, 'O'},
    {"omit-empty-lines", 0, NULL, 'x'},
    {"warnings-to",      1, NULL, 'l'},
    {"output",           1, NULL, 'o'},
    {"version",          0, NULL, 'v'},
    {0, 0, 0, 0}
  };
  int option_index=0;
  char *endptr;
  int optch; 

  /*
    We start by loading the default values
    for the user settable options.
  */
  list->optmask = 0x0;
  list->begin_line=1;
  list->end_line=0;
  list->step=1;
  list->first_field=1;
  list->last_field=0;
  list->increment=1;
  list->ifs = NULL;
  list->osep = NULL;
  list->file = NULL;

  while ( (optch = getopt_long (argc, argv, optstring, long_options, &option_index)) != -1 )
    {
      switch (optch)
	{
	case 'h':
	  list->optmask |= ___H_MASK;
	  break;
	case 'b':
	  list->optmask |= ___B_MASK;
	  errno = 0;
	  for (endptr = optarg; is_space (*endptr) != 0; endptr++);
	  if (*endptr == '-' ||
	      (list->begin_line = strtoul (optarg, &endptr, 10), 
	       errno == ERANGE) || list->begin_line == 0 ||
	      *endptr != '\0')
	    {
	      fprintf (stderr, _("%s: invalid argument after `-%c\' option\n"),
		       PACKAGE2, optch);
	      return -1;
	    }
	  break;
	case 'e':
	  list->optmask |= ___E_MASK;
	  errno = 0;
	  for (endptr = optarg; is_space (*endptr) != 0; endptr++);
	  if (*endptr == '-' ||
	      (list->end_line = strtoul (optarg, &endptr, 10), 
	       errno == ERANGE) ||
	      *endptr != '\0')
	    {
	      fprintf (stderr, _("%s: invalid argument after `-%c\' option\n"),
		       PACKAGE2, optch);
	      return -1;
	    }
	  break;
	case 's':
	  list->optmask |= ___S_MASK;
	  errno = 0;
	  for (endptr = optarg; is_space (*endptr) != 0; endptr++);
	  if (*endptr == '-' ||
	      (list->step = strtoul (optarg, &endptr, 10), 
	       errno == ERANGE) || list->step == 0 ||
	      *endptr != '\0')
	    {
	      fprintf (stderr, _("%s: invalid argument after `-%c\' option\n"),
		       PACKAGE2, optch);
	      return -1;
	    }
	  break;
	case 'F':
	  list->optmask |= ___SF_MASK;
	  errno = 0;
	  for (endptr = optarg; is_space (*endptr) != 0; endptr++);
	  if (*endptr == '-' ||
	      (list->first_field = strtoul (optarg, &endptr, 10), 
	       errno == ERANGE) || list->first_field == 0 ||
	      *endptr != '\0')
	    {
	      fprintf (stderr, _("%s: invalid argument after `-%c\' option\n"),
		       PACKAGE2, optch);
	      return -1;
	    }
	  break;
	case 'L':
	  list->optmask |= ___SL_MASK;
	  errno = 0;
	  for (endptr = optarg; is_space (*endptr) != 0; endptr++);
	  if (*endptr == '-' ||
	      (list->last_field = strtoul (optarg, &endptr, 10), 
	       errno == ERANGE) ||
	      *endptr != '\0')
	    {
	      fprintf (stderr, _("%s: invalid argument after `-%c\' option\n"),
		       PACKAGE2, optch);
	      return -1;
	    }
	  break;
	case 'I':
	  list->optmask |= ___SI_MASK;
	  errno = 0;
	  for (endptr = optarg; is_space (*endptr) != 0; endptr++);
	  if (*endptr == '-' ||
	      (list->increment = strtoul (optarg, &endptr, 10), 
	       errno == ERANGE) || list->increment == 0 ||
	      *endptr != '\0')
	    {
	      fprintf (stderr, _("%s: invalid argument after `-%c\' option\n"),
		       PACKAGE2, optch);
	      return -1;
	    }
	  break;
        case 'S':
          if ((list->ifs))
	    {
              delete_string_vector (list->ifs);
	      list->ifs = NULL;
	    }
          list->ifs = ssplit_former_way (optarg);
          if (!list->ifs)
            {
              fprintf (stderr, _("%s: memory exhausted\n"), PACKAGE2);
              return -1;
            }
	  else if ( !is_string_in_vector(NEWLINE_STR, (const char**) list->ifs) )
            {
              fprintf (stderr, _("%s: invalid argument after `-%c\' option:\n"),
                       PACKAGE2, optch);
              fprintf (stderr, _("  The list of field delimiters cannot be empty and\n  must always include the newline character (\'\\n\')\n"));
              return -1;
            }
          else
            {
              remove_duplicates_from_string_vector (list->ifs);
              sort_string_vector (list->ifs); /* This is not strictly necessary */
              list->optmask |= ___SS_MASK;
            }
          break;
        case 'D':
          if ((list->ifs))
	    {
              delete_string_vector (list->ifs);
	      list->ifs = NULL;
	    }
          list->ifs = ssplit (optarg, I_DEF_SEP);
          if (!list->ifs)
            {
              fprintf (stderr, _("%s: memory exhausted\n"), PACKAGE2);
              return -1;
            }
	  else if ( !is_string_in_vector(NEWLINE_STR, (const char**) list->ifs) )
            {
              fprintf (stderr, _("%s: invalid argument after `-%c\' option:\n"),
                       PACKAGE2, optch);
              fprintf (stderr, _("  The list of field delimiters cannot be empty and\n  must always include the newline string (\"\\n\")\n"));
              return -1;
            }
          else
            {
              remove_duplicates_from_string_vector (list->ifs);
              sort_string_vector (list->ifs);
              list->optmask |= ___SD_MASK;
            }
          break;
        case 'O':
          if ((list->osep))
	    {
              free((void*)list->osep);
	      list->osep = NULL;
	    }
          list->osep = get_separating_string (optarg);
          if (!list->osep)
            {
              fprintf (stderr, _("%s: memory exhausted\n"), PACKAGE2);
              return -1;
            }
          else
            list->optmask |= ___SO_MASK;          
          break;
	case 'x':
	  list->optmask |= ___X_MASK;
	  break;
	case 'l':
	  if (!freopen (optarg, "w", stderr))
	    {
	      fprintf (stderr, _("%s: cannot open file \"%s\":\n"),
		       PACKAGE2, optarg);
	      perror(0);
	      return -1;
	    }
	  break;
	case 'o':
	  if (!freopen (optarg, "w", stdout))
	    {
	      fprintf (stderr, _("%s: cannot open file \"%s\":\n"),
		       PACKAGE2, optarg);
	      perror(0);
	      return -1;
	    }
	  break;
	case 'v':
	  list->optmask |= ___V_MASK;
	  break;
	default:
/*  	  fprintf (stderr,  */
/*  		   _("%s: unrecognized option `-%c\' \n"), PACKAGE2, optch);  */
	  return -1;
	}
    }
  if (!(list->optmask & (___H_MASK | ___V_MASK)) && argc - optind > 1)
    {
      print_selhelp (PACKAGE2);
      return -1;
    }
  else
    {
      if( !(list->optmask & (___H_MASK | ___V_MASK)) )
	list->file = (const char*) argv[optind];
      return 0;
    }
}

static
char** def_ifs = NULL;

static 
void clean_up_memory (void)
{
  delete_string_vector (def_ifs);
}

int main (int argc, char** argv)
{
  Argslist arg_list;

  def_ifs = ssplit (DEF_IFS, I_DEF_SEP);
  if ((atexit(clean_up_memory)) || !def_ifs)
    {
      fprintf (stderr, _("%s: memory exhausted\n"), PACKAGE2);
      return -1;
    }
#ifdef ENABLE_NLS
  setlocale (LC_CTYPE, "");
  setlocale (LC_MESSAGES, "");
#endif
  bindtextdomain (PACKAGE2, LOCALEDIR);
  textdomain (PACKAGE2);
  if ( set_args (argc, argv, &arg_list) != 0 )
    return -1;
  else if ( (arg_list.optmask & (___H_MASK | ___V_MASK)) )
    {
      if ((arg_list.optmask & ___V_MASK))
      	print_selversion(PACKAGE2);
      if ((arg_list.optmask & ___H_MASK))
      	print_selhelp(PACKAGE2);
      if (argc > 2)
	return -1;
      else
	return 0;
    }
  else
    {
      switch (scan_file (def_ifs, &arg_list))
	{
	case OPEN_ERROR:
	  fprintf (stderr, _("%s: cannot open file \"%s\":\n"), PACKAGE2, arg_list.file);
	  perror(0);
	  return -1;
	case CLOSE_ERROR:
	  fprintf (stderr, _("%s: cannot close file \"%s\":\n"), PACKAGE2, arg_list.file);
	  perror(0);
	  return -1;
	case READ_ERROR:
	  fprintf (stderr, _("%s: Error occurred while reading from file \"%s\"\n\n"), PACKAGE2, arg_list.file);
	  return -1;
	default:
	  return 0;
	}
    }
}
