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

#define GDIFF_OPTIONS 1

/* Leave this inclusion at the begin, otherwise problems */
/* with the symbol __USE_FILE_OFFSET64                   */
#include"numdiff.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<ctype.h>
#include<float.h>
#include"getopt.h"
#include"error.h"
#include"linesplit.h"
#include"xalloc.h"

#ifdef _DMALLOC_
#include <dmalloc.h> /* Useful only for the debugging */
#endif

void print_version (const char* progname)
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
#ifdef USE_GMP
  printf ("\n%s %s.\n", _("The software has been linked against\n\
the GNU Multiple Precision Arithmetic Library,\n\
version number"), gmp_version);
#else  /* not USE_GMP */
  printf ("\n%s.\n", _("The software has been built with\n\
its own internal support for multiple precision arithmetic"));
#endif /* not USE_GMP */
}

void print_help (const char* progname)
{
  puts (_("Usage:"));
  printf ("%s -h|--help|-v|--version   %s\n\n", progname, _("or"));
  printf ("%s %s\n", progname, "[-s IFS][-D DELIMS][-a THRVAL[:RANGE|:RANGE1:RANGE2]][-r THRVAL[:RANGE|:RANGE1:RANGE2]][-2][-F NUM][-# NUM][-P][-N][-I][-c CURRNAME][-d C1C2][-t C1C2][-g N1N2][-p C1C2][-n C1C2][-e C1C2][-i C1C2][-X 1:RANGE][-X 2:RANGE][-E][-U][-b][-V][-O[NUM]][-q][-S][-z 1:RANGE][-z 2:RANGE][-Z 1:RANGE][-Z 2:RANGE][-m][-H][-f[NUM]][-T][-B][-l PATH][-o PATH] FILE1 FILE2");
  printf (_("\nCompare putatively similar files line by line and field by field,\nignoring small numeric differences or/and different numeric formats.\n\n"));
  printf (_("RANGE, RANGE1 and RANGE2 stay for a positive integer value or\nfor a range of integer values, like 1-, 3-5 or -7.\n"));
  printf ("%s\n%s\n%s\n\n%s\n\n",
	  _("The two arguments after the options are the names of the files to compare."),
	  _("The complete paths of the files should be given,\na directory name is not accepted."),
	  _("The given paths cannot refer to the same file\nbut one of them can be \"-\", which refers to stdin."),
	  _("Exit status: 1 if files differ, 0 if they are equal, -1 (255) in case of error"));
  /* %%% */
  printf ("-s, --separators=IFS\n    %s\n    %s\n    %s\n",
	  _("Specify the set of characters to use as delimiters\n    while splitting the input lines into fields"),
	  _("(The default set of delimiters is space, tab and newline)."),
	  _("If IFS is prefixed with 1: or 2: then use the given delimiter set\n    only for the lines from the first or the second file respectively"));
  printf ("-D, --delimiters=DELIMS\n    %s\n    %s\n    %s\n",
	  _("Specify the set of strings to use as delimiters\n    while splitting the input lines into fields"),
	  _("(The default set of delimiters is space, tab and newline)."),
	  _("If DELIMS is prefixed with 1: or 2: then use the given delimiter set\n    only for the lines from the first or the second file respectively"));
  printf ("-a, --absolute-tolerance=THRVAL[:RANGE|:RANGE1:RANGE2]\n    %s\n    %s\n    %s\n", 
	  _("Set to THRVAL the maximum absolute difference permitted\n    before that two numeric fields are regarded as different\n    (The default value is zero)."),
	  _("If a RANGE is given, use the specified\n    threshold only when comparing fields whose positions lie in RANGE."),
	  _("If both RANGE1 and RANGE2 are given and have the same length,\n    then use the specified threshold when comparing a field of FILE1\n    lying in RANGE1 with the corresponding field of FILE2 in RANGE2"));
  printf ("-r, --relative-tolerance=THRVAL[:RANGE|:RANGE1:RANGE2]\n    %s\n    %s\n    %s\n", 
	  _("Set to THRVAL the maximum relative difference permitted\n    before that two numeric fields are regarded as different\n    (The default value is zero)."),
	  _("If a RANGE is given, use the specified\n    threshold only when comparing fields whose positions lie in RANGE."),
	  _("If both RANGE1 and RANGE2 are given and have the same length,\n    then use the specified threshold when comparing a field of FILE1\n    lying in RANGE1 with the corresponding field of FILE2 in RANGE2"));
  printf ("-2, --strict\n    %s\n",
	  _("Consider two numerical values as equal only if\n    both absolute and relative difference do not exceed\n    the corresponding tolerance threshold"));
  printf ("-F, --formula=NUM\n    %s\n    %s\n    %s\n    %s\n",
	  _("Use the formula indicated by NUM to compute the relative errors."),
	  _("If \'NUM\' is 0 use the classic formula."),
	  _("If \'NUM\' is 1 compute the relative errors by considering\n    the values in FILE1 as sample values."),
	  _("If \'NUM\' is 2 compute the relative errors by considering\n    the values in FILE2 as sample values."));
  printf ("-#, --digits=NUM\n    %s\n",
	  _("Set to NUM the number of digits in the significands\n    used in multiple precision arithmetic"));
  printf ("-P, --positive-differences\n    %s\n",
	  _("Ignore all differences due to numeric fields of the second file that\n    are less than the corresponding numeric fields in the first file"));
  printf ("-N, --negative-differences\n    %s\n",
	  _("Ignore all differences due to numeric fields of the second file that\n    are greater than the corresponding numeric fields in the first file"));
  printf ("-I, --ignore-case\n    %s\n",
	  _("Ignore changes in case while doing literal comparisons"));
  printf ("-c, --currency=CURRNAME\n    %s\n    %s\n",
	  _("Set to CURRNAME the currency name for the two files to compare."),
	  _("CURRNAME must be prefixed with 1: or 2: to specify the\n    currency name only for the first or the second file"));
  printf ("-d, --decimal-point=C1C2\n    %s\n",
	  _("Specify the characters representing the decimal point\n    in the two files to compare"));
  printf ("-t, --thousands-separator=C1C2\n    %s\n",
	  _("Specify the characters representing the thousands separator\n    in the two files to compare"));
  printf ("-g, --group-length=N1N2\n    %s\n",
	  _("Specify the number of digits forming each group of thousands\n    in the two files to compare"));
  printf ("-p, --plus-prefix=C1C2\n    %s\n",
	  _("Specify the (optional) prefixes for positive values\n    used in the two files to compare"));
  printf ("-n, --minus-prefix=C1C2\n    %s\n",
	  _("Specify the prefixes for negative values\n    used in the two files to compare"));
  printf ("-e, --exponent-letter=C1C2\n    %s\n",
	  _("Specify the exponent letters\n    used in the two files to compare"));
  printf ("-i, --imaginary-unit=C1C2\n    %s\n",
	  _("Specify the characters representing the imaginary unit\n    in the two files to compare"));
  printf ("-X, --exclude=1:RANGE\n    %s\n",
	  _("Select the fields of the first file that have to be ignored"));
  printf ("-X, --exclude=2:RANGE\n    %s\n",
	  _("Select the fields of the second file that have to be ignored"));
  printf ("-E, --essential\n    %s\n",
	  _("While printing the differences between the two compared files\n    show only the numerical ones"));
  printf ("-U, --dummy\n    %s\n",
	  _("While printing the differences between the two compared files\n    neglect all the numerical ones (dummy mode)"));
  printf ("-b, --brief\n    %s\n",
	  _("Suppress all messages concerning the differences discovered\n    in the structures of the two files"));
  printf ("-V, --verbose\n    %s\n",
	  _("For every couple of lines which differ in at least one field print\n    an header to show how these lines appear in the two compared files"));
  printf ("-O, --overview[=NUM]\n    %s\n    %s\n    %s\n    %s\n",
	  _("Display a side by side difference listing of the two files\n    showing which lines are present only in one file, which\n    lines are present in both files but with one or more differing fields,\n    and which lines are identical."),
	  _("If \'NUM\' is zero or is not specified, output at most 130 columns per line."),
	  _("If \'NUM\' is a positive number, output at most \'NUM\' columns per line."),
	  _("If \'NUM\' is a negative number, do not output common lines\n    and display at most -\'NUM\' columns per line."));
  printf ("-q, --quiet, --silent\n    %s\n",
	  _("Suppress all the standard output"));
  printf ("-S, --statistics\n    %s\n",
	  _("Add some statistics to the standard output"));
  printf ("-z, --blur-if-numerical=1:RANGE\n    %s\n",
	  _("Select the fields of the first file that have to be\n    blurred during the synchronization procedure\n    only if they turn out to be numeric"));
  printf ("-z, --blur-if-numerical=2:RANGE\n    %s\n",
	  _("Select the fields of the second file that have to be\n    blurred during the synchronization procedure\n    only if they turn out to be numeric"));
  printf ("-Z, --blur-unconditionally=1:RANGE\n    %s\n",
	  _("Select the fields of the first file that have to be\n    unconditionally blurred during the synchronization procedure"));
  printf ("-Z, --blur-unconditionally=2:RANGE\n    %s\n",
	  _("Select the fields of the second file that have to be\n    unconditionally blurred during the synchronization procedure"));
  printf ("-m, --minimal\n    %s\n",
	  _("During synchronization try hard to find a smaller set of changes"));
  printf ("-H, --speed-large-files\n    %s\n",
	  _("During synchronization assume large files and\n    many scattered small changes"));
  printf ("-f, --test-filter[=NUM]\n    %s\n    %s\n    %s\n    %s\n",
	  _("Run only the filter and then show the results of its\n    attempt to synchronize the two files."),
	  _("If \'NUM\' is zero or is not specified, output at most 130 columns per line."),
	  _("If \'NUM\' is a positive number, output at most \'NUM\' columns per line."),
	  _("If \'NUM\' is a negative number, do not output common lines\n    and display at most -\'NUM\' columns per line."));
  printf ("-T, --expand-tabs\n    %s\n",
	  _("Expand tabs to spaces in output while displaying the results of the\n    synchronization procedure (meaningful only together with option -O or -f)"));
  printf ("-B, --binary\n    %s\n",
	  _("Treat both files as binary files (only meaningful under Doz/Windoz)"));
  printf ("-l, --warnings-to=PATH\n    %s\n",
	  _("Redirect warning and error messages from stderr to the indicated file"));
  printf ("-o, --output=PATH\n    %s\n",
	  _("Redirect output from stdout to the indicated file"));
  printf ("-h, --help\n    %s\n", _("Show help message and predefined settings"));
  printf ("-v, --version\n    %s\n", _("Show version number, Copyright, Distribution Terms and NO-Warranty"));
  /* %%% */
  puts (_("\n  Default numeric format (for both files to compare):\n"));
  printf (_("Currency name = \"%s\"\n"), CURRENCY);
  printf (_("Decimal point = `%c\'\n"), DP);

  printf (_("Thousands separator = `%c\'\n"), THSEP);
  printf (_("Number of digits in each thousands group = %u\n"), GROUPING);

  printf (_("Leading positive sign = `%c\'\n"), POS_SIGN);
  printf (_("Leading negative sign = `%c\'\n"), NEG_SIGN);
  printf (_("Prefix for decimal exponent = `%c\'\n"), ECH);
  printf (_("Symbol used to denote the imaginary unit = `%c\'\n\n"), IU);
}

static
int nfset (int opt_ch, const char* opt_arg, argslist* arg_list)
{
  if (strlen(opt_arg) <= 2)
    {
      char _1st = *opt_arg, _2nd = *(opt_arg+1);

      switch (opt_ch)
	{
	case 'd':
	  if ( (is_punct(_1st)) && (_2nd == '\0' || is_punct(_2nd)) )
	    {
	      arg_list->optmask |= _D_MASK;
	      arg_list->nf1.dp = _1st;
	      arg_list->nf2.dp = (_2nd) ? _2nd : _1st; 
	      return 0;
	    }
	  break;
	case 't':
	  if ( (is_punct(_1st)) && (_2nd == '\0' || is_punct(_2nd)) )
	    {
	      arg_list->optmask |= _T_MASK;
	      arg_list->nf1.thsep = _1st;
	      arg_list->nf2.thsep = (_2nd) ? _2nd : _1st;
	      return 0;
	    }
	  break;
	case 'e':
	  if ( is_print(_1st) && (_2nd == '\0' || is_print(_2nd)) )
	    {
	      arg_list->optmask |= _E_MASK;
	      arg_list->nf1.ech = _1st;
	      arg_list->nf2.ech = (_2nd) ? _2nd : _1st;
	      return 0;
	    }
	  break;
	case 'n':
	  if ( is_print(_1st) && (_2nd == '\0' || is_print(_2nd)) )
	    {
	      arg_list->optmask |= _N_MASK;
	      arg_list->nf1.neg_sign = _1st;
	      arg_list->nf2.neg_sign = (_2nd) ? _2nd : _1st;
	      return 0;
	    }
	  break;
	case 'i':
	  if ( is_print(_1st) && (_2nd == '\0' || is_print(_2nd)) )
	    {
	      arg_list->optmask |= _I_MASK;
	      arg_list->nf1.iu = _1st;
	      arg_list->nf2.iu = (_2nd) ? _2nd : _1st;
	      return 0;
	    }
	  break;
	case 'p':
	  if ( is_print(_1st) && (_2nd == '\0' || is_print(_2nd)) )
	    {
	      arg_list->optmask |= _P_MASK;
	      arg_list->nf1.pos_sign = _1st;
	      arg_list->nf2.pos_sign = (_2nd) ? _2nd : _1st;
	      return 0;
	    }
	  break;
	case 'g':
	  if ( (is_digit(_1st)) && (_2nd == '\0' || is_digit(_2nd)) )
	    {
	      arg_list->optmask |= _G_MASK;
	      arg_list->nf1.grouping = _1st - '0';
	      arg_list->nf2.grouping = (_2nd) ? _2nd - '0': _1st - '0';
	      return 0;
	    }
	  break;
	}
    }
  return -1;
}

static
int fselect (const char* str, unsigned char* mask, int mask_size)
{
  long beg, end;
  unsigned long n;
  char *ptr, *endptr;

  beg = end = -1;
  if (!str || !*str)
    return 0; /* no field selected */
  /* If we arrive here we are sure that *str != '\0' ! */
  if ( strcmp (str, "@") == 0 )
    {
      /* select all fields */
      for (mask_size /= 8; mask_size > 0; mask_size--, mask[mask_size] = 0xFF);       
      return 1;
    }
  if ((beg = strtol (str, &endptr, 10)) == 0
      || beg > mask_size || beg < -mask_size)
    return -1; /* illegal input */
  else if (beg < 0)
    {
      if (*endptr == '\0')
	{
	  end = -beg;
	  beg = 1;
	}
      else
	return -1;
    }
  else if (*endptr == '\0')
    end = beg;
  else if (*endptr == '-')
    {
      if (*(ptr = endptr + 1) == '\0')
	end = mask_size; 
      else
	{
	  if ((end = strtol (ptr, &endptr, 10)) <= 0
	      || *endptr != '\0' || end > mask_size)
	    return -1; /* illegal input */
	}
    }
  if (beg > end)
    return -1;
  else
    {
      /* Remark: internally the field numbers
	 start from zero, not from one */
      for (n = beg - 1; n <= end - 1; n++)
	mask[n >> 3] |= 0x80 >> (n & 0x7);
      return 1;
    }
}

#define VALID_NUMFMT      0
#define INVALID_NUMFMT   -1
#define INVALID_CURRENCY -2

static
int is_numfmt_valid (const struct numfmt* pnf)
{
  char store[NUMFMT_CHARS];
  char *ptr;
  int i, j;

  for (ptr = pnf->currency; 
       *ptr != '\0' && 
	 (!is_digit(*ptr)) &&
	 *ptr != pnf->dp &&
	 *ptr != pnf->thsep &&
	 *ptr != pnf->pos_sign &&
	 *ptr != pnf->neg_sign; ptr++);
  if (*ptr != '\0')
    return INVALID_CURRENCY;
  store[0] = pnf->dp;
  store[1] = pnf->thsep;
  store[2] = pnf->pos_sign;
  store[3] = pnf->neg_sign;
  store[4] = pnf->ech;
  store[5] = pnf->iu;
  for (i=0; i < NUMFMT_CHARS; i++)
    {
      for (j = i+1; j < NUMFMT_CHARS; j++)
	if (store[i] == store[j])
	  return INVALID_NUMFMT;
    }
  return VALID_NUMFMT;
}

extern int optind;

int setargs (int argc, char* argv[], argslist *list)
{
  const int mask_size = FIELDMASK_SIZE*8;
  const char *optstring = "h2F:bBVO::qUESIPNz:Z:mHT#:s:D:a:r:c:d:t:g:p:n:e:i:f::X:l:o:v";
  struct option long_options[] = {
    {"help",                 0, NULL, 'h'},
    {"strict",               0, NULL, '2'},
    {"formula",              1, NULL, 'F'},
    {"brief",                0, NULL, 'b'},
    {"binary",               0, NULL, 'B'},
    {"verbose",              0, NULL, 'V'},
    {"overview",             2, NULL, 'O'},
    {"quiet",                0, NULL, 'q'},
    {"silent",               0, NULL, 'q'},
    {"dummy",                0, NULL, 'U'},
    {"essential",            0, NULL, 'E'},
    {"statistics",           0, NULL, 'S'},
    {"ignore-case",          0, NULL, 'I'},
    {"positive-differences", 0, NULL, 'P'},
    {"negative-differences", 0, NULL, 'N'},
    {"blur-if-numerical",    1, NULL, 'z'},
    {"blur-unconditionally", 1, NULL, 'Z'},
    {"minimal",              0, NULL, 'm'},
    {"speed-large-files",    0, NULL, 'H'},
    {"expand-tabs",          0, NULL, 'T'},
    {"digits",               1, NULL, '#'},
    {"separators",           1, NULL, 's'},
    {"delimiters",           1, NULL, 'D'},
    {"absolute-tolerance",   1, NULL, 'a'},
    {"relative-tolerance",   1, NULL, 'r'},
    {"currency",             1, NULL, 'c'},
    {"decimal-point",        1, NULL, 'd'},
    {"thousands-separator",  1, NULL, 't'},
    {"group-length",         1, NULL, 'g'},
    {"plus-prefix",          1, NULL, 'p'},
    {"minus-prefix",         1, NULL, 'n'},
    {"exponent-letter",      1, NULL, 'e'},
    {"imaginary-unit",       1, NULL, 'i'},
    {"test-filter",          2, NULL, 'f'},
    {"exclude",              1, NULL, 'X'},
    {"warnings-to",          1, NULL, 'l'},
    {"output",               1, NULL, 'o'},
    {"version",              0, NULL, 'v'},
    {0, 0, 0, 0}
  };
  int option_index=0;
  char *tail;
  int i, optch, off, rv; 
  unsigned int t, file_id;
  long w;
  unsigned char *bitmask;  

  /*
    We start by loading the default values
    for the user settable options.

    The initialization of the variables

    list->maxrelerr, list->maxabserr,
    list->Labserr, list->Crelerr, list->Lrelerr, list->Cabserr,
    list->N1abserr, list->N1disperr, list->N2abserr, list->N2disperr 

    is done within main() through init_mpa_support().
  */

  binary = 0;
  suppress_common_lines = 0;
  ignore_white_space = IGNORE_NO_WHITE_SPACE;
  expand_tabs = 0;
  w = DEF_ATMOST_NCOLS;
  speed_large_files = 0;
  program_name = PACKAGE;

  list->optmask = 0x0;
  list->output_mode = OUTMODE_NORMAL;
  for (i=0; i < FIELDMASK_SIZE; 
       list->ghostmask1[i] = list->ghostmask2[i] = list->tblurmask1[i] = list->tblurmask2[i] = list->pblurmask1[i] = list->pblurmask2[i] = 0x0, i++); 
  list->relerr_formula = CLASSIC_FORMULA;
  list->Nentries = list->Ndisperr = 0;
  list->flag = 0;
  list->ifs1 = list->ifs2 = NULL;
  list->iscale = ISCALE;
  list->nf1.dp = DP;
  list->nf1.thsep = THSEP;
  list->nf1.grouping = GROUPING;
  list->nf1.pos_sign = POS_SIGN;
  list->nf1.neg_sign = NEG_SIGN;
  list->nf1.ech = ECH;
  list->nf1.iu = IU;
  list->file1 = list->file2 = NULL;
  list->nf2 = list->nf1;
  list->nf1.currency = get_separating_string (CURRENCY);
  list->nf2.currency = get_separating_string (CURRENCY);

  while ( (optch = getopt_long (argc, argv, optstring, long_options, &option_index)) != -1 )
    {
      switch (optch)
	{
	case 'h':
	  list->optmask |= _H_MASK;
	  break;
	case '2':
	  list->optmask |= _2_MASK;
	  break;	  
	case 'F':
	  if (strncmp ("0", optarg, 1) == 0)
	    list->relerr_formula = CLASSIC_FORMULA;	      
	  else if (strncmp ("1", optarg, 1) == 0)
	    list->relerr_formula = WR_TO_FIRST_FILE;	      
	  else if (strncmp ("2", optarg, 1) == 0)
	    list->relerr_formula = WR_TO_SECOND_FILE;
	  else
	    {
	      fprintf (stderr, _("%s: invalid argument after `-%c\' option\n"),
		       PACKAGE, optch);
	      return -1;
	    }
	  break;
	case 'b':
	  list->optmask |= _B_MASK;
	  break;
	case 'B':
	  binary = 1;
	  break;
	case 'V':
	  list->optmask |= _SV_MASK;
	  break;
	case 'q':
	  list->optmask |= _Q_MASK;
	  break;
	case 'U':
	  list->optmask |= _SU_MASK;
	  break;
	case 'E':
	  list->optmask |= _SE_MASK;
	  break;
	case 'S':
	  list->optmask |= _SS_MASK;
	  break;
	case 'I':
	  list->optmask |= _SI_MASK;
	  break;
	case 'P':
	  list->optmask |= _SP_MASK;
	  list->flag = 1;
	  break;
	case 'N':
	  list->optmask |= _SN_MASK;
	  list->flag = -1;
	  break;
	case 'z':
	  if ( (i = strncmp (optarg, "1:", 2)) && (strncmp (optarg, "2:", 2)) )
	    {
	      /*
		None of the prefixes 1: and 2: has been used,
		then we have to select fields for both files
	      */
	      if (fselect (optarg, list->pblurmask1, mask_size) <= 0)
		{
		  fprintf (stderr, _("%s: invalid argument after `-%c\' option\n"),
			   PACKAGE, optch);
		  return -1;
		}
	      else
		{
		  fselect (optarg, list->pblurmask2, mask_size);
		  list->optmask |= _Z_MASK;
		}
	    }
	  else 
	    {
	      bitmask = i == 0 ? list->pblurmask1 : list->pblurmask2;
	      if (fselect (optarg+2, bitmask, mask_size) <= 0)
		{
		  fprintf (stderr, _("%s: invalid argument after `-%c\' option\n"),
			   PACKAGE, optch);
		  return -1;
		}
	      else
		list->optmask |= _Z_MASK;
	    }
	  break;
	case 'Z':
	  if ( (i = strncmp (optarg, "1:", 2)) && (strncmp (optarg, "2:", 2)) )
	    {
	      /*
		None of the prefixes 1: and 2: has been used,
		then we have to select fields for both files
	      */
	      if (fselect (optarg, list->tblurmask1, mask_size) <= 0)
		{
		  fprintf (stderr, _("%s: invalid argument after `-%c\' option\n"),
			   PACKAGE, optch);
		  return -1;
		}
	      else
		{
		  fselect (optarg, list->tblurmask2, mask_size);
		  list->optmask |= _SZ_MASK;
		}
	    }
	  else 
	    {
	      bitmask = i == 0 ? list->tblurmask1 : list->tblurmask2;
	      if (fselect (optarg+2, bitmask, mask_size) <= 0)
		{
		  fprintf (stderr, _("%s: invalid argument after `-%c\' option\n"),
			   PACKAGE, optch);
		  return -1;
		}
	      else
		list->optmask |= _SZ_MASK;
	    }
	  break;
	case 'm':
	  list->optmask |= _M_MASK;
	  break;
	case 'H':
	  list->optmask |= _SH_MASK;
	  speed_large_files = 1;
	  break;
	case 'T':
	  expand_tabs = 1;
	  break;
	case '#':
	  list->iscale = strtol (optarg, &tail, 10);
	  if (*tail != '\0' || list->iscale < 0 || list->iscale > MAX_ISCALE)
	    {
	      fprintf (stderr, _("%s: invalid argument after `-%c\' option\n"),
		       PACKAGE, optch);
	      return -1;
	    }
	  else
	    list->optmask |= _X_MASK;	  
	  break;
	case 's':
	  if (*optarg == '1' && *(optarg+1) == ':')
	    {
              if ((list->ifs1))
                {
                  delete_string_vector (list->ifs1);
                  list->ifs1 = NULL;
                }
              list->ifs1 = ssplit_former_way (optarg+2);
              if (!list->ifs1)
                {
                  fprintf (stderr, _("%s: memory exhausted\n"), PACKAGE);
                  return -1;
                }
	    }
	  else if (*optarg == '2' && *(optarg+1) == ':')
	    {
              if ((list->ifs2))
                {
                  delete_string_vector (list->ifs2);
                  list->ifs2 = NULL;
                }
              list->ifs2 = ssplit_former_way (optarg+2);
              if (!list->ifs2)
                {
                  fprintf (stderr, _("%s: memory exhausted\n"), PACKAGE);
                  return -1;
                }
	    }
	  else
	    {
              if ((list->ifs1))
                {
                  delete_string_vector (list->ifs1);
                  list->ifs1 = NULL;
                }
              if ((list->ifs2))
                {
                  delete_string_vector (list->ifs2);
                  list->ifs2 = NULL;
                }
              list->ifs1 = ssplit_former_way (optarg);
              list->ifs2 = ssplit_former_way (optarg);
              if (!list->ifs1 || !list->ifs2)
                {
                  fprintf (stderr, _("%s: memory exhausted\n"), PACKAGE);
                  return -1;
                }
	    }

          if ( ((list->ifs1) && !is_string_in_vector(NEWLINE_STR, (const char**) list->ifs1)) ||
               ((list->ifs2) && !is_string_in_vector(NEWLINE_STR, (const char**) list->ifs2)) )
	    {
	      fprintf (stderr, _("%s: invalid argument after `-%c\' option:\n"),
		       PACKAGE, optch);
	      fprintf (stderr, _("  The list of field delimiters can not be empty and\n  must always include the newline character (\'\\n\')\n"));
	      return -1;
	    }
	  else
            {
              remove_duplicates_from_string_vector (list->ifs1);
              remove_duplicates_from_string_vector (list->ifs2);
              sort_string_vector (list->ifs1); /* This is not strictly necessary */
              sort_string_vector (list->ifs2); /* This is not strictly necessary */              
              list->optmask |= _S_MASK;
            }
	  break;
        case 'D':
	  if (*optarg == '1' && *(optarg+1) == ':')
	    {
              if ((list->ifs1))
                {
                  delete_string_vector (list->ifs1);
                  list->ifs1 = NULL;
                }
              list->ifs1 = ssplit (optarg+2, I_DEF_SEP);
              if (!list->ifs1)
                {
                  fprintf (stderr, _("%s: memory exhausted\n"), PACKAGE);
                  return -1;
                }
	    }
	  else if (*optarg == '2' && *(optarg+1) == ':')
	    {
              if ((list->ifs2))
                {
                  delete_string_vector (list->ifs2);
                  list->ifs2 = NULL;
                }
              list->ifs2 = ssplit (optarg+2, I_DEF_SEP);
              if (!list->ifs2)
                {
                  fprintf (stderr, _("%s: memory exhausted\n"), PACKAGE);
                  return -1;
                }
	    }
	  else
	    {
              if ((list->ifs1))
                {
                  delete_string_vector (list->ifs1);
                  list->ifs1 = NULL;
                }
              if ((list->ifs2))
                {
                  delete_string_vector (list->ifs2);
                  list->ifs2 = NULL;
                }
              list->ifs1 = ssplit (optarg, I_DEF_SEP);
              list->ifs2 = ssplit (optarg, I_DEF_SEP);
              if (!list->ifs1 || !list->ifs2)
                {
                  fprintf (stderr, _("%s: memory exhausted\n"), PACKAGE);
                  return -1;
                }
	    }
          if ( ((list->ifs1) && is_char_in_vector(NEWLINE, (const char**) list->ifs1) != 1) ||
               ((list->ifs2) && is_char_in_vector(NEWLINE, (const char**) list->ifs2) != 1) )
	    {
	      fprintf (stderr, _("%s: invalid argument after `-%c\' option:\n"),
		       PACKAGE, optch);
	      fprintf (stderr, _("  The list of field delimiters cannot be empty and\n  must always include the newline string (\"\\n\").\n  Care that the newline character cannot appear\n  in any other delimiter than the newline string\n"));
	      return -1;
	    }
	  else
            {
              remove_duplicates_from_string_vector (list->ifs1);
              remove_duplicates_from_string_vector (list->ifs2);
              sort_string_vector (list->ifs1); 
              sort_string_vector (list->ifs2); 
              list->optmask |= _S_MASK;
            }
          break;
	case 'a':
	  rv = thrlist_add (&list->maxabserr, optarg);
	  if (rv == THRLIST_INVALID_FORMAT)
	    {
	      fprintf (stderr, _("%s: invalid argument after `-%c\' option:\n"),
		       PACKAGE, optch);
	      fprintf (stderr, _("  The format specification has not been respected\n"));
	      return -1;
	    }
	  else if (rv == THRLIST_INVALID_RANGES)
	    {
	      fprintf (stderr, _("%s: invalid argument after `-%c\' option:\n"),
		       PACKAGE, optch);
	      fprintf (stderr, _("  The specified ranges do not have the same length\n"));
	      return -1;
	    }
	  else
	    list->optmask |= _A_MASK;
	  break;
	case 'r':
	  rv = thrlist_add (&list->maxrelerr, optarg);
	  if (rv == THRLIST_INVALID_FORMAT)
	    {
	      fprintf (stderr, _("%s: invalid argument after `-%c\' option:\n"),
		       PACKAGE, optch);
	      fprintf (stderr, _("  The format specification has not been respected\n"));
	      return -1;
	    }
	  else if (rv == THRLIST_INVALID_RANGES)
	    {
	      fprintf (stderr, _("%s: invalid argument after `-%c\' option:\n"),
		       PACKAGE, optch);
	      fprintf (stderr, _("  The specified ranges do not have the same length\n"));
	      return -1;
	    }
	  else
	    list->optmask |= _R_MASK;
	  break;
	case 'c':
	  file_id = 0x0;
	  if (*optarg == '1' && *(optarg+1) == ':')
	    {
	      if ((list->nf1.currency))
		free((void*)list->nf1.currency);
	      list->nf1.currency = get_separating_string (optarg+2);
              file_id = 0x1;
	    }
	  else if (*optarg == '2' && *(optarg+1) == ':')
	    {
	      if ((list->nf2.currency))
		free((void*)list->nf2.currency);
	      list->nf2.currency = get_separating_string (optarg+2);
	      file_id = 0x2;
	    }
	  else
	    {
	      if ((list->nf1.currency))
		free((void*)list->nf1.currency);
	      list->nf1.currency = get_separating_string (optarg);
	      if ((list->nf2.currency))
		free((void*)list->nf2.currency);
	      list->nf2.currency = get_separating_string (optarg);
	      file_id = 0x3;
	    }
	  if ( (strlen(list->nf1.currency) == 0 && ((file_id & 0x3) == 0x1))
	    || (strlen(list->nf2.currency) == 0 && ((file_id & 0x3) == 0x2)) )
	    {
	      fprintf (stderr, _("%s: invalid argument after `-%c\' option:\n"),
		       PACKAGE, optch);
	      fprintf (stderr, _("  you have missed to specify the currency name\n"));
	      return -1;
	    }
	  break;
	case 'd':
	case 't':
	case 'g':
	case 'p':
	case 'n':
	case 'e':
	case 'i':
	  if (nfset (optch, optarg, list) < 0)
	    {
	      fprintf (stderr, _("%s: invalid argument after `-%c\' option\n"),
		       PACKAGE, optch);
	      return -1;
	    }
	  break;
	case 'f':
	case 'O':
	  if (!(list->optmask & (_F_MASK | _SO_MASK)))
	    {
	      if(!optarg)
		{
		  /* There is no optional argument, then set */
		  /* 'w' to 'DEF_ATMOST_NCOLS'.              */
		  if (optch == 'f')
		    list->optmask |= _F_MASK;
		  else
		    list->optmask |= _SO_MASK;
		  w = DEF_ATMOST_NCOLS;
		}
	      else
		{
		  /* An argument follows */
		  w = strtol (optarg, &tail, 10);
		  /* If the argument of the option is not a valid number, */
		  /* then exit after printing a suitable error message.   */
		  if (*tail != '\0')
		    {
		      fprintf (stderr, _("%s: invalid argument after `-%c\' option\n"),
			       PACKAGE, optch);
		      return -1;
		    }
		  else
		    {
		      if (optch == 'f')
			list->optmask |= _F_MASK;
		      else
			list->optmask |= _SO_MASK;
		    }
		  /* Otherwise you have to set 'w' appropriately. */
		  /* If the given argument is less than -MAX_ATMOST_NCOLS */
		  /* then set 'w' to 'DEF_ATMOST_NCOLS' and 'suppress_common_lines' to 'TRUE'. */
		  if (w < -MAX_ATMOST_NCOLS)
		    {
		      w = DEF_ATMOST_NCOLS;
		      suppress_common_lines = 1;
		    }
		  /* If the argument were negative, then remove the sign */
		  /* and set 'suppress_common_lines' to 'TRUE'.          */
		  if (w < 0)
		    {
		      w *= -1;
		      suppress_common_lines = 1;
		    }
		  /* If the given argument is too small or too big in absolute value, */
		  /* then set 'w' to 'DEF_ATMOST_NCOLS'.                              */
		  if (w < MIN_ATMOST_NCOLS || w > MAX_ATMOST_NCOLS)
		    w = DEF_ATMOST_NCOLS;
		  /* Otherwise leave 'w' set to the value of the argument. */
		} /* end optarg != 0 */
	    }
	  break;
	case 'X':
	  if ( (i = strncmp (optarg, "1:", 2)) && (strncmp (optarg, "2:", 2)) )
	    {
	      /*
		None of the prefixes 1: and 2: has been used,
		then we have to select fields for both files
	      */
	      if (fselect (optarg, list->ghostmask1, mask_size) <= 0)
		{
		  fprintf (stderr, _("%s: invalid argument after `-%c\' option\n"),
			   PACKAGE, optch);
		  return -1;
		}
	      else
		{
		  fselect (optarg, list->ghostmask2, mask_size);
		  list->optmask |= _SX_MASK;
		}
	    }
	  else 
	    {
	      bitmask = i == 0 ? list->ghostmask1 : list->ghostmask2;
	      if (fselect (optarg+2, bitmask, mask_size) <= 0)
		{
		  fprintf (stderr, _("%s: invalid argument after `-%c\' option\n"),
			   PACKAGE, optch);
		  return -1;
		}
	      else
		list->optmask |= _SX_MASK; 
	    }
	  break;
	case 'l':
	  if (!freopen (optarg, "w", stderr))
	    {
	      fprintf (stderr, _("%s: cannot open file \"%s\":\n"),
		       PACKAGE, optarg);
	      perror(0);
	      return -1;
	    }
	  break;
	case 'o':
	  if (!freopen (optarg, "w", stdout))
	    {
	      fprintf (stderr, _("%s: cannot open file \"%s\":\n"),
		       PACKAGE, optarg);
	      perror(0);
	      return -1;
	    }
	  break;
	case 'v':
	  list->optmask |= _V_MASK;
	  break;
	default:
	  /* 	  
		  fprintf (stderr, 
		  _("%s: unrecognized option `-%c\' \n"), PACKAGE, optch); 
	  */
	  return -1;
	}
    }

  t = expand_tabs ? 1 : TAB_WIDTH;
  off = (w + t + 3) / (2 * t)  *  t;
  sdiff_half_width = MAX (0, MIN (off - 3, w - off)),
    sdiff_column2_offset = sdiff_half_width ? off : w;

  if ( list->optmask & _SV_MASK )
    list->output_mode = OUTMODE_VERBOSE;
  if ( list->optmask & _B_MASK )
    list->output_mode = OUTMODE_BRIEF;
  if (list->optmask & _B_MASK && list->optmask & _SV_MASK)
    list->output_mode = OUTMODE_COINCISE;
  if ( list->optmask & _SO_MASK )
    list->output_mode = OUTMODE_OVERVIEW;
  if ( list->optmask & _Q_MASK )
    list->output_mode = OUTMODE_QUIET;

  if (!(list->optmask & (_H_MASK | _V_MASK)) && argc - optind != 2)
    {
      print_help (PACKAGE);
      return -1;
    }
  else if ( (rv = is_numfmt_valid(&list->nf1)) == INVALID_NUMFMT )
    {
      fprintf (stderr, 
	       _("The numeric format specified for the first file is illegal,\n"));
      fprintf (stderr,
	       _("the following symbols should be all different\nwhile two or more of them are actually equal:\n"));
      fprintf (stderr, _("\nDecimal point = `%c\'\n"), list->nf1.dp);
      fprintf (stderr, _("Thousands separator = `%c\'\n"), list->nf1.thsep);
      fprintf (stderr, _("Leading positive sign = `%c\'\n"), list->nf1.pos_sign);
      fprintf (stderr, _("Leading negative sign = `%c\'\n"), list->nf1.neg_sign);
      fprintf (stderr, _("Prefix for decimal exponent = `%c\'\n"), 
	       list->nf1.ech);
      fprintf (stderr, 
	       _("Symbol used to denote the imaginary unit = `%c\'\n\n"), 
	       list->nf1.iu);
      return -1;
    }
  else if ( rv == INVALID_CURRENCY )
    {
      fprintf (stderr, 
	       _("The numeric format specified for the first file is illegal:\n"));
      fprintf (stderr,
	       _("the name of the currency may not contain digits,\n"));
      fprintf (stderr,
	       _("the symbol for the leading positive sign (`%c\'),\n"), list->nf1.pos_sign);
      fprintf (stderr,
	       _("the symbol for the leading negative sign (`%c\'),\n"), list->nf1.neg_sign);
      fprintf (stderr,
	       _("the decimal point (`%c\'), or the thousands separator (`%c\')\n"), list->nf1.dp, list->nf1.thsep);
      return -1;
    }
  else if ( (rv = is_numfmt_valid(&list->nf2)) == INVALID_NUMFMT )
    {
      fprintf (stderr, 
	       _("The numeric format specified for the second file is illegal,\n"));
      fprintf (stderr,
	       _("the following symbols should be all different\nwhile two or more of them are actually equal:\n"));
      fprintf (stderr, _("\nDecimal point = `%c\'\n"), list->nf2.dp);
      fprintf (stderr, _("Thousands separator = `%c\'\n"), list->nf2.thsep);
      fprintf (stderr, _("Leading positive sign = `%c\'\n"), list->nf2.pos_sign);
      fprintf (stderr, _("Leading negative sign = `%c\'\n"), list->nf2.neg_sign);
      fprintf (stderr, _("Prefix for decimal exponent = `%c\'\n"), 
	       list->nf2.ech);
      fprintf (stderr, 
	       _("Symbol used to denote the imaginary unit = `%c\'\n\n"), 
	       list->nf2.iu);
      return -1;
    }
  else if ( rv == INVALID_CURRENCY )
    {
      fprintf (stderr, 
	       _("The numeric format specified for the second file is illegal:\n"));
      fprintf (stderr,
	       _("the name of the currency may not contain digits,\n"));
      fprintf (stderr,
	       _("the symbol for the leading positive sign (`%c\'),\n"), list->nf2.pos_sign);
      fprintf (stderr,
	       _("the symbol for the leading negative sign (`%c\'),\n"), list->nf2.neg_sign);
      fprintf (stderr,
	       _("the decimal point (`%c\'), or the thousands separator (`%c\')\n"), list->nf2.dp, list->nf2.thsep);
      return -1;
    }
  else
    {
      if( !(list->optmask & (_H_MASK | _V_MASK)) )
	{
	  list->file1 = (const char*) argv[optind];
	  list->file2 = (const char*) argv[optind+1];
	}
      return 0;
    }
}  
