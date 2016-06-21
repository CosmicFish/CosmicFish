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

/* Leave this inclusion at the begin, otherwise problems */
/* with the symbol __USE_FILE_OFFSET64                   */
#include"numdiff.h" 
#include"linesplit.h"
#include<stdio.h>
#include<stdlib.h> /* for free() */
#include<string.h>
#if HAVE_GETTIMEOFDAY
#include<sys/time.h>
#endif
#ifdef ENABLE_NLS
#include<locale.h>
#endif

/* See cmpfns.c */
extern int cmp_files (FILE* pf1, FILE* pf2, argslist* argl);

/* See options.c */
extern void print_version (const char* progname);
extern void print_help (const char* progname);
extern int setargs (int argc, char* argv[], argslist *list);

static
void init_mpa_support (argslist* list)
{
  init_mpa(list->iscale);
  initR (&list->Labserr);
  initR (&list->Crelerr);
  initR (&list->Lrelerr);
  initR (&list->Cabserr);
  initR (&list->N1abserr);
  initR (&list->N1disperr);
  initR (&list->N2abserr);
  initR (&list->N2disperr);
  list->maxabserr = thrlist_new ();
  list->maxrelerr = thrlist_new ();
}

static
void dismiss_mpa_support (argslist* list)
{
  delR (&list->Labserr);
  delR (&list->Crelerr);
  delR (&list->Lrelerr);
  delR (&list->Cabserr);
  delR (&list->N1abserr);
  delR (&list->N1disperr);
  delR (&list->N2abserr);
  delR (&list->N2disperr);
  thrlist_dispose (&list->maxabserr);
  thrlist_dispose (&list->maxrelerr);
  end_mpa();
}

static void
set_mtime_to_now (struct stat *st)
{
#ifdef ST_MTIM_NSEC

# if HAVE_CLOCK_GETTIME && defined CLOCK_REALTIME
  if (clock_gettime (CLOCK_REALTIME, &st->st_mtim) == 0)
    return;
# endif

# if HAVE_GETTIMEOFDAY
  {
    struct timeval timeval;
    if (gettimeofday (&timeval, NULL) == 0)
      {
	st->st_mtime = timeval.tv_sec;
	st->st_mtim.ST_MTIM_NSEC = timeval.tv_usec * 1000;
	return;
      }
  }
# endif

#endif /* ST_MTIM_NSEC */

  time (&st->st_mtime);
}

/* cmp.file[f].desc markers */
#define NONEXISTENT (-1)   /* nonexistent file */
#define UNOPENED (-2)      /* unopened file (e.g. directory) */
#define ERRNO_ENCODE(errno) (-3 - (errno)) /* encoded errno value */

#define ERRNO_DECODE(desc) (-3 - (desc))   /* inverse of ERRNO_ENCODE */
#define DIR_P(f) (S_ISDIR (files[f].stat.st_mode) != 0)

static 
int open_files (const char* name0, const char* name1)
{
  register int f;
  int status = EXIT_SUCCESS;
  bool same_files;

  if (!name0 || !name1)
    return EXIT_TROUBLE;

  memset (files, 0, sizeof files);
  files[0].desc = UNOPENED;
  files[1].desc = UNOPENED;
  files[0].name = name0;
  files[1].name = name1;

  /* Stat the files.  */

  for (f = 0; f < 2; f++)
    {
      if ((f) && file_name_cmp (files[f].name, files[0].name) == 0)
	{
	  files[f].desc = files[0].desc;
	  files[f].stat = files[0].stat;
	}
      else if (strcmp (files[f].name, "-") == 0)
	{
	  files[f].desc = STDIN_FILENO;
	  if (fstat (STDIN_FILENO, &files[f].stat) != 0)
	    files[f].desc = ERRNO_ENCODE (errno);
	  else
	    {
	      if (S_ISREG (files[f].stat.st_mode))
		{
		  off_t pos = lseek (STDIN_FILENO, (off_t) 0, SEEK_CUR);
		  if (pos < 0)
		    files[f].desc = ERRNO_ENCODE (errno);
		  else
		    files[f].stat.st_size =
		      MAX (0, files[f].stat.st_size - pos);
		}

	      /* POSIX 1003.1-2001 requires current time for
		 stdin.  */
	      set_mtime_to_now (&files[f].stat);
	    }
	}
      else if (stat (files[f].name, &files[f].stat) != 0)
	files[f].desc = ERRNO_ENCODE (errno);
    }

  for (f = 0; f < 2; f++)
    {
      int e = ERRNO_DECODE (files[f].desc);
      if (0 <= e)
	{
	  errno = e;
	  perror_with_name (files[f].name);
	  status = EXIT_TROUBLE;
	}
    }

  if (status != EXIT_SUCCESS)
    /* One of the files should exist but does not.  */
    return status;
  else if (DIR_P (0) | DIR_P (1))
    return EXIT_TROUBLE;
  else
    {
      /* Both exist and neither is a directory.  */
      /* Are they the same file ?                */
      same_files
	= (files[0].desc != NONEXISTENT
	   && files[1].desc != NONEXISTENT
	   && 0 < same_file (&files[0].stat, &files[1].stat)
	   && same_file_attributes (&files[0].stat,
				    &files[1].stat));

      /* Open the files and record their descriptors.  */

      if (files[0].desc == UNOPENED)
	if ((files[0].desc = open (files[0].name, O_RDONLY, 0)) < 0)
	  {
	    perror_with_name (files[0].name);
	    status = EXIT_TROUBLE;
	  }
      if (files[1].desc == UNOPENED)
	{
	  if ((same_files))
	    files[1].desc = files[0].desc;
	  else if ((files[1].desc = open (files[1].name, O_RDONLY, 0))
		   < 0)
	    {
	      perror_with_name (files[1].name);
	      status = EXIT_TROUBLE;
	    }
	}

#if HAVE_SETMODE_DOS
      if (binary)
	for (f = 0; f < 2; f++)
	  if (0 <= files[f].desc)
	    set_binary_mode (files[f].desc, 1);
#endif
      return status;
    }
}

static 
int compare_files (argslist* list, int* is_same_physical_file)
{
  if ((files[0].desc != NONEXISTENT
       && files[1].desc != NONEXISTENT
       && 0 < same_file (&files[0].stat, &files[1].stat)
       && same_file_attributes (&files[0].stat,
				   &files[1].stat)))
    {
      /* The two named files are actually the same physical file.
	 We know they are identical without actually reading them.  */
      *is_same_physical_file = 1;
      return 0;
    }
  else
    {
      int status = diff_2_files (files, list);
      /* 
         STATUS is 0 if no changes have been found,
         1 in case of detected changes, -1 if either 
         file is binary. 
      */

      *is_same_physical_file = 0;
      return status; 
    }
}

static
int rewind_files (void)
{
  off_t pos0, pos1;
  int status = EXIT_SUCCESS;

  if ((pos0 = lseek(files[0].desc, (off_t) 0, SEEK_SET)) < 0)
    {
      perror_with_name (files[0].name);
      status = EXIT_TROUBLE;
    }
  if ((pos1 = lseek(files[1].desc, (off_t) 0, SEEK_SET)) < 0)
    {
      perror_with_name (files[1].name);
      status = EXIT_TROUBLE;
    }
  return status;
}

static
int set_file_pointers (FILE** fpp1, FILE** fpp2)
{
  int status = EXIT_SUCCESS;

  if ( !(*fpp1 = fdopen (files[0].desc, "r")) )
    {
      perror_with_name (files[0].name);
      status = EXIT_TROUBLE;
    }
  if ( !(*fpp2 = fdopen (files[1].desc, "r")) )
    {
      perror_with_name (files[1].name);
      status = EXIT_TROUBLE;
    }
  return status;
}

static
int close_files (void)
{
  /* Close the file descriptors.  */
  
  if (0 <= files[0].desc && close (files[0].desc) != 0)
    {
      perror_with_name (files[0].name);
      return EXIT_TROUBLE;
    }
  if (0 <= files[1].desc && files[0].desc != files[1].desc
      && close (files[1].desc) != 0)
    {
      perror_with_name (files[1].name);
      return EXIT_TROUBLE;
    }
  return EXIT_SUCCESS;
}

static
void print_statistics (argslist* list)
{
  Real qm_abserr, qm_relerr;

#ifdef USE_GMP
  initR (&qm_abserr);
  initR (&qm_relerr);
#endif /* USE_GMP */
  if (list->flag > 0)
    {
      fputs (_("\n  In the computation of the following quantities\n  only the errors with positive sign are considered:\n"),
	     stdout);
      fputs (_("  differences due to numeric fields of the second file that are\n  less than the corresponding fields in the first file are neglected\n\n"),
	     stdout);
    }
  if (list->flag < 0)
    {
      fputs (_("\n  In the computation of the following quantities\n  only the errors with negative sign are considered:\n"),
	     stdout);
      fputs (_("  differences due to numeric fields of the second file that are\n  greater than the corresponding fields in the first file are neglected\n\n"),
	     stdout);
    }
  if ( list->Ndisperr == 0 )
    {
      if ( list->Nentries == 0 )
	fputs (_("\nNo numeric comparison has been done\n"),
	       stdout);
      else
	printf(ngettext (
			 "\nOne numeric comparison has been done and\nthe resulting numeric difference is negligible\n", 
			 "\n%d numeric comparisons have been done and\nthe resulting numeric differences are all negligible\n",
			 list->Nentries), list->Nentries);
    }
  else if ( list->Ndisperr == list->Nentries )
    {
      printf(ngettext (
		       "\nOne numeric comparison has been done and\nhas produced an outcome beyond the tolerance threshold\n",
		       "\n%d numeric comparisons have been done, all of them\nhave produced an outcome beyond the tolerance threshold\n",
		       list->Nentries), list->Nentries);
    }
  else
    {
      /* Case  0 < LIST->NDISPERR < LIST->NENTRIES */
      printf (ngettext (
			"\nOne numeric comparison has been done,\n",
			"\n%d numeric comparisons have been done,\n",
			list->Nentries), list->Nentries);
      
      printf (ngettext (
			"only one numeric comparison has produced an outcome\nbeyond the tolerance threshold\n",
			"%d numeric comparisons have produced an outcome\nbeyond the tolerance threshold\n",
			list->Ndisperr), list->Ndisperr);
    }

  fputs (_("\nLargest absolute error in the set of relevant numerical differences:\n"),
	 stdout);
  printno (list->Labserr, DEF_LIM);
  fputs (_("\nCorresponding relative error:\n"), stdout);
  printno (list->Crelerr, DEF_LIM);

  fputs (_("\nLargest relative error in the set of relevant numerical differences:\n"),
	 stdout);
  printno (list->Lrelerr, DEF_LIM);
  fputs (_("\nCorresponding absolute error:\n"), stdout);
  printno (list->Cabserr, DEF_LIM);
  putchar ('\n');

  fputs (_("\nSum of all absolute errors:\n"),
	 stdout);
  printno (list->N1abserr, DEF_LIM);
  fputs (_("\nSum of the relevant absolute errors:\n"),
	 stdout);
  printno (list->N1disperr, DEF_LIM);
  /* Arithmetic means */
  divide_by_int (&list->N1abserr, list->Nentries, list->iscale);
  divide_by_int (&list->N1disperr, list->Ndisperr, list->iscale);
  fputs (_("\nArithmetic mean of all absolute errors:\n"),
	 stdout);
  printno (list->N1abserr, DEF_LIM);
  fputs (_("\nArithmetic mean of the relevant absolute errors:\n"),
	 stdout);
  printno (list->N1disperr, DEF_LIM);

  /* 2-norms and quadratic means of the errors */
  copyR (&qm_abserr, list->N2abserr);
  divide_by_int (&qm_abserr, list->Nentries, list->iscale);
  square_root (&qm_abserr, list->iscale);
  square_root (&list->N2abserr, list->iscale);
  fputs (_("\nSquare root of the sum of the squares of all absolute errors:\n"),
	 stdout);
  printno (list->N2abserr, DEF_LIM);
  fputs (_("\nQuadratic mean of all absolute errors:\n"),
	 stdout);
  printno (qm_abserr, DEF_LIM);

  copyR (&qm_relerr, list->N2disperr);
  divide_by_int (&qm_relerr, list->Ndisperr, list->iscale);
  square_root (&qm_relerr, list->iscale);
  square_root (&list->N2disperr, list->iscale);
  fputs (_("\nSquare root of the sum of the squares\nof the relevant absolute errors:\n"),
	 stdout);
  printno (list->N2disperr, DEF_LIM);
  fputs (_("\nQuadratic mean of the relevant absolute errors:\n"),
	 stdout);
  printno (qm_relerr, DEF_LIM);
  putchar ('\n');
  delR (&qm_relerr);
  delR (&qm_abserr);
}

char **def_ifs = NULL;

int main (int argc, char* argv[])
{
  argslist list;

#ifdef ENABLE_NLS
  setlocale (LC_CTYPE, "");
  setlocale (LC_MESSAGES, "");
#endif
  bindtextdomain (PACKAGE, LOCALEDIR);
  textdomain (PACKAGE);

  def_ifs = ssplit (DEF_IFS, I_DEF_SEP);
  if (!def_ifs)
    {
      fprintf (stderr, _("***  %s: memory exhausted\n"), PACKAGE);
      return -1;
    }

  /* This code was used to discover the reason of a bug */  
  /*
    #ifdef __USE_FILE_OFFSET64
    printf ("\n %s: FILE OFFSET 64 in use, sizeof (struct stat) = %u\n", __FILE__,
    sizeof (struct stat));
    #else
    printf ("\n %s: FILE OFFSET 64 NOT in use, sizeof (struct stat) = %u\n", __FILE__,
    sizeof (struct stat));
    #endif
  */

  init_mpa_support (&list);
  init_flags ();

  if ( setargs (argc, argv, &list) != 0 )
    {
      delete_string_vector (def_ifs);
      dismiss_mpa_support (&list);
      if ((list.ifs1))
	delete_string_vector (list.ifs1);
      if ((list.ifs2))
	delete_string_vector (list.ifs2);
      return -1;
    }
  else if ( (list.optmask & (_H_MASK | _V_MASK)) )
    {
      if ((list.optmask & _V_MASK))
	print_version(PACKAGE);
      if ((list.optmask & _H_MASK))
	print_help(PACKAGE);
      delete_string_vector (def_ifs);
      dismiss_mpa_support (&list);
      if ((list.ifs1))
	delete_string_vector (list.ifs1);
      if ((list.ifs2))
	delete_string_vector (list.ifs2);
      if (argc > 2)
	return -1;
      else
	return 0;
    }
  else
    {
      int test = 0, ident_files = 0;
      FILE *fp1, *fp2;
      int qm = list.optmask & _Q_MASK;

      if ( open_files (list.file1, list.file2) != EXIT_SUCCESS )
        {
          delete_string_vector (def_ifs);
          dismiss_mpa_support (&list);
          if ((list.ifs1))
            delete_string_vector (list.ifs1);
          if ((list.ifs2))
            delete_string_vector (list.ifs2);
          return EXIT_TROUBLE;
        }
      if ( list.optmask & (_F_MASK | _Z_MASK | _SZ_MASK) )
	test = compare_files (&list, &ident_files);

      if (test < 0)
        {
          fputs (_("\n***  The requested comparison cannot be performed:\n"), stdout);
	  printf (_("***  At least one between \"%s\" and \"%s\" is a binary file\n"),
                  list.file1, list.file2);
          return EXIT_TROUBLE;
        }

      if ( !(list.optmask & _F_MASK) && !ident_files )
	{
	  if ( rewind_files () != EXIT_SUCCESS ||
	       set_file_pointers (&fp1, &fp2) != EXIT_SUCCESS )
            {
              delete_string_vector (def_ifs);
              dismiss_mpa_support (&list);
              if ((list.ifs1))
                delete_string_vector (list.ifs1);
              if ((list.ifs2))
                delete_string_vector (list.ifs2);
              close_files ();
              return EXIT_TROUBLE;
            }
	  test = cmp_files (fp1, fp2, &list);
	  if ((list.optmask & _SS_MASK) && test <= 1)
	    print_statistics (&list);
	}
      if (test == 0 && !qm)
	{
	  if ((list.optmask & _F_MASK))
	    printf (_("\n+++  Files \"%s\" and \"%s\" have the same structure\n"),
		    list.file1, list.file2);
	  else
	    printf (_("\n+++  Files \"%s\" and \"%s\" are equal\n"),
		    list.file1, list.file2);
	}
      if (test == 1 && !qm)
	printf (_("\n+++  File \"%s\" differs from file \"%s\"\n"),
		list.file1, list.file2);
      erase_flags ();
      delete_string_vector (def_ifs);
      dismiss_mpa_support (&list);
      if ((list.ifs1))
	delete_string_vector (list.ifs1);
      if ((list.ifs2))
	delete_string_vector (list.ifs2);
      close_files ();
      return test;
    }
}
