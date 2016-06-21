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

#ifndef _NUMDIFF_H_
#define _NUMDIFF_H_

#include "system.h"

#if defined(HAVE_LIBGMP) && !defined(DISABLE_GMP)
#define USE_GMP 1
#endif /* defined(HAVE_LIBGMP) && !defined(DISABLE_GMP) */

/* The type of a hash value.  */
typedef size_t hash_value;
verify (hash_value_is_unsigned, ! TYPE_SIGNED (hash_value));

/* Rotate an unsigned value to the left.  */
#define ROL(v, n) ((v) << (n) | (v) >> (sizeof (v) * CHAR_BIT - (n)))

/* Given a hash value and a new character, return a new hash value.  */
#ifdef _DEBUG_HASHING_
#define HASH(h, c) (putc((c), stderr), (c) + ROL (h, 7))
#else
#define HASH(h, c) ((c) + ROL (h, 7))
#endif

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

/*
  Begin section -- Math types
*/

#include"number.h"

#ifdef USE_GMP

#include<gmp.h>

typedef mpf_t Real;
typedef struct {
  mpf_t re, im;
} Complex;

#else /* not USE_GMP */

typedef bc_num Real;
typedef struct {
  bc_num re, im;
} Complex;

#endif /* USE_GMP */

/*
  End section -- Math types
*/

#define THRLIST_OK              0
#define THRLIST_INVALID_FORMAT -1
#define THRLIST_INVALID_RANGES -2

/*
  A structure of type 'thrlist_node' is used
  to construct lists of threshold values.
*/

struct __thrlist_node {
  Real threshold;               /* the threshold value */

  /* BEG1-END1 and BEG2-END2 are the ranges of fields    */
  /* to which the specified THRESHOLD value applies.     */
  /* The two ranges must have the same length, but refer */
  /* to different files.                                 */

  unsigned long beg1, end1;     
  unsigned long beg2, end2;

  /* DOUBLE_RANGE_SPEC is a boolean value, which is 1 if */
  /* BEG1-END1 and BEG2-END2  have been both explicitely */
  /* set via the command line through a specification of */
  /* the form:                                           */
  /*                 BEG1-END1:BEG2-END2                 */
  /* and otherwise is equal to 0.                        */

  int double_range_spec;

  /* pointer to the next node in the list */
  struct __thrlist_node *next;  
};

typedef struct __thrlist_node thrlist_node;

/* A pointer to a list of threshold values */
typedef thrlist_node *thrlist;

struct numfmt {
  char *currency;    /* currency            */
  char dp;           /* decimal point       */
  char thsep;        /* thousands separator */
  unsigned grouping; /* Number of digits in each group */
  char pos_sign;     /* positive sign */
  char neg_sign;     /* negative sign */
  char ech;          /* prefix for decimal exponent  */
  char iu;           /* symbol of the imaginary unit */ 
}; 
/* A structure of this type is used to store  */
/* information about the legal format for the */
/* numbers in input.                          */

/*
  This is the number of the fields in the 'numftm' structure
  having "char" type. None of them can be a digit and 
  they must have all different values
  (see the code of the function valid_numfmt() in the file options.c).
*/
#define NUMFMT_CHARS 6

typedef struct {
  unsigned char* ptr;
  size_t len, size;
} flg_array; /* A structure of this type is used to store the information
		retrieved from the execution of a diff command */

typedef struct {
  /* Mask of the options */
  unsigned long optmask;

  /* Output mode. This field can take any of the values */
  /* OUTMODE_* (see below)                              */
  int output_mode;

  /* This is a mask of bits specifying the fields of the first file 
     which must be ignored. */
  unsigned char ghostmask1[FIELDMASK_SIZE];

  /* This is a mask of bits specifying the fields of the second file 
     which must be ignored. */
  unsigned char ghostmask2[FIELDMASK_SIZE];

  /* This is a mask of bits specifying the fields in the 
     first file for which partial blurring must be enabled during 
     the filtering procedure. */
  unsigned char pblurmask1[FIELDMASK_SIZE];

  /* This is a mask of bits specifying the fields in the
     second file for which partial blurring must be enabled during 
     the filtering procedure. */
  unsigned char pblurmask2[FIELDMASK_SIZE];

  /* This is a mask of bits specifying the fields in the
     first file for which total blurring must be enabled during 
     the filtering procedure. */
  unsigned char tblurmask1[FIELDMASK_SIZE];

  /* This is a mask of bits specifying the fields in the
     second file for which total blurring must be enabled during 
     the filtering procedure. */
  unsigned char tblurmask2[FIELDMASK_SIZE];

  /* This parameter specifies how the relative errors */
  /* have to be computed. It may have one of the      */
  /* following values:                                */
  /* CLASSIC_FORMULA            0  (the default one)  */
  /* WR_TO_FIRST_FILE           1                     */
  /* WR_TO_SECOND_FILE          2                     */
  int relerr_formula;

  /* Tolerance thresholds for absolute and relative errors */
  thrlist maxabserr, maxrelerr;

  /* These variables are used to print statistics */
  Real Labserr,  Crelerr,  Lrelerr,  Cabserr, N1abserr, N1disperr, N2abserr, N2disperr;
  int Nentries, Ndisperr;

  /* Flag > 0 --> If a numeric field in the first file is greater than the  */
  /*              corresponding numeric field in the second file, then the  */
  /*              related difference is ignored, i.e. it is never output    */
  /* Flag < 0 --> If a numeric field in the first file is less than the     */
  /*              corresponding numeric field in the second file, then the  */
  /*              related difference is ignored, i.e. it is never output    */
  /* Flag = 0 --> Standard behavior: the difference between 2 corresponding */
  /*              numeric fields (one in the first file, the other one in   */
  /*              the second file) is always considered and it is output    */
  /*              whenever its absolute value is greater than the given     */
  /*              tolerance thresholds for the absolute and relative errors.*/
  signed char flag;

  /* Internal scale (number of digits of accuracy) */
  int iscale;

  /* Files to be compared */
  const char *file1, *file2;

  /* Internal fields separators (IFS) for file1 and file2 */
  char **ifs1, **ifs2;

  /* Numeric conventions for file1 (.nf1) and file2 (.nf2) */
  struct numfmt nf1, nf2;

} argslist ; /* A structure of this type is used to store the options */
/* set by the user                                                    */ 

#define _H_MASK   0x00000001 /* -h option, used to recall help */
#define _A_MASK   0x00000002 /* -a option, used to set tolerance for abs. error  */
#define _R_MASK   0x00000004 /* -r option, used to set tolerance for rel. error  */
#define _2_MASK   0x00000008 /* -2 option, used to enable the "strict" control */
#define _S_MASK   0x00000010 /* -s option, used to explicitly set IFS  */
#define _B_MASK   0x00000020 /* -b option, used to enable the "brief" mode */
#define _F_MASK   0x00000040 /* -f option, used to enable the "filter-only" mode */
#define _Q_MASK   0x00000080 /* -q option, used to enable "quiet" mode */
#define _X_MASK   0x00000100 /* -# option, used to set the precision   */
#define _D_MASK   0x00000200 /* -d option, used to set the decimal point */
#define _T_MASK   0x00000400 /* -t option, used to set the thousands separator */
#define _G_MASK   0x00000800 /* -g option, used to set the 'grouping' */
#define _P_MASK   0x00001000 /* -p option, used to set the character 
		  	    	'positive sign' */
#define _N_MASK   0x00002000 /* -n option, used to set the character 
		  	    	'negative sign' */
#define _E_MASK   0x00004000 /* -e option, used to set prefix for 
		  	    	decimal exponent */
#define _I_MASK   0x00008000 /* -i option, used to set the symbol
		  	    	of the imaginary unit */
#define _L_MASK   0x00010000 /* -l option, used to redirect the standard
		  	    	error on a file */
#define _O_MASK   0x00020000 /* -o option, used to redirect the standard
		  	    	output on a file */
#define _Z_MASK   0x00040000 /* -z option, used to activate the filter 
			        (normal mode)    */
#define _SZ_MASK  0x00080000 /* -Z option, used to activate the filter 
                                (alternative mode) */
#define _SX_MASK  0x00100000 /* -X option, used to select which fields in the
		  	    	lines of the files must be ignored */
#define _SP_MASK  0x00200000 /* -P option, used to ignore negative errors */
#define _SN_MASK  0x00400000 /* -N option, used to ignore positive errors */
#define _SU_MASK  0x00800000 /* -U option, used to enable the "dummy" mode */
#define _SE_MASK  0x01000000 /* -E option, used to enable the "essential" 
				mode */
#define _SV_MASK  0x02000000 /* -V option, used to enable the "verbose" mode */
#define _SO_MASK  0x04000000 /* -O option, used to enable the "overview" mode */
#define _SS_MASK  0x08000000 /* -S option, used to print statistics */
#define _SI_MASK  0x10000000 /* -I option, used to ignore case while comparing
			        non numerical fields */
#define _SH_MASK  0x20000000 /* -H option, by filtering assume large files and
				many scattered small changes */
#define _M_MASK   0x40000000 /* -m option, by filtering try hard to find a smaller 
				set of changes */
#define _V_MASK   0x80000000 /* -v option, used to show version number,
			        Copyright and No-Warrany */
/* Remark: Intentionally there are no masks for the options -c and -F */

/* Output modes: verbose, normal, brief, and quiet.       */
/* Do not change the relative order of the values of      */
/* these macros, the code in cmp_lines()  (see file       */
/* cmpfns.c) relies on the fact that:                     */
/* OUTMODE_VERBOSE > OUTMODE_NORMAL > OUTMODE_COINCISE    */
/* > OUTMODE_BRIEF > OUTMODE_QUIET  > OUTMODE_OVERVIEW    */

#define OUTMODE_VERBOSE     4
#define OUTMODE_NORMAL      3
#define OUTMODE_COINCISE    2
#define OUTMODE_BRIEF       1
#define OUTMODE_QUIET       0
#define OUTMODE_OVERVIEW   -1

/* Methods to compute the relative differences */

#define CLASSIC_FORMULA         0
#define WR_TO_FIRST_FILE        1
#define WR_TO_SECOND_FILE       2

#ifndef PACKAGE
#define PACKAGE "numdiff"
#endif

#ifndef LOCALEDIR
#define LOCALEDIR "/usr/local/share/locale/"
#endif

/* The character representing the number zero */

#define CHAR_ZERO '0'

/* The character representing the number one */

#define CHAR_ONE '1'

/* The character representing the number nine */

#define CHAR_NINE '9'

/* newline character */

#define NEWLINE '\n'

/*
  Predefined values for
  .nf*.currency (currency name)
  .nf*.dp       (decimal point)
  .nf*.thsep    (thousands separator)
  .nf*.grouping (number of digits in each thousands group)
  .nf*.pos_sign (positive sign)
  .nf*.neg_sign (negative sign)
  .nf*.ech      (prefix for decimal exponent)
  .nf*.iu       (symbol of the imaginary unit)
  .iscale      (decimal digits of accuracy)
*/
#define CURRENCY  ""
#define DP        '.'
#define THSEP     ','
#define GROUPING   3
#define POS_SIGN  '+'
#define NEG_SIGN  '-'
#define ECH       'e'
#define IU        'i'
#define ISCALE    35


/*
  Largest possible value for .iscale
*/
#define MAX_ISCALE  180

/*
  Largest possible exponent accepted by Numdiff
  when a number is written in scientific notation
*/
#define MAX_EXPN    +1073741824L

/*
  Lowest possible exponent accepted by Numdiff
  when a number is written in scientific notation
*/
#define MIN_EXPN    -1073741824L

/*
  Macro to move ahead a pointer
*/
#define move_ahead(ptr) ptr++

/*
  Character classification macros.
  The macro CTYPE_DOMAIN is defined in "system.h".
*/
#define is_digit(c) ((unsigned int) (c) - '0' <= 9 ? 1 : 0)
#define is_punct(c) (CTYPE_DOMAIN((unsigned char)(c)) && ispunct((unsigned char)(c)))
#define is_print(c) (CTYPE_DOMAIN((unsigned char)(c)) && isgraph((unsigned char)(c)) && ((unsigned int) (c) - '0' > 9))
#define is_space(c) (CTYPE_DOMAIN((unsigned char)(c)) && isspace((unsigned char)(c)))

/*
  Mathematical functions
*/

int      cmp (Real p, Real q);
int      is0 (Real u);
int      smart_cmp (const Complex* pz1, const Complex* pz2, int flag);
void     printno (Real u, int m);

extern Real Zero, Inf;

void     init_mpa(int iscale);

void     initR (Real* px);
void     initC (Complex* pz);

void copyR (Real* dst, Real src);
void copyC (Complex* dst, Complex src);

#ifdef _MPA_DEBUG
void     debug_printno (Real u, int m);
#endif /* _MPA_DEBUG */

void     str2R (const char *q, char **endptr, int iscale,
		const struct numfmt* pnf, Real* pr);
void     str2C (const char *q, char **endptr, int iscale,
		const struct numfmt* pnf, Complex* pc);

void     add (Real s, Real t, Real* q, int iscale);
void     square (Real s, Real* q, int iscale);
void     divide (Real s, Real t, Real* q, int iscale);
void     divide_by_int (Real* q, int d, int iscale);
void     square_root (Real* q, int iscale);

void     Cabs (Complex z, Real* pm, int iscale);
void     Csub (Complex z1, Complex z2, Complex* pw, int iscale);

void     delR (Real* px);
void     delC (Complex* pz);

void     end_mpa(void);

/*
  Functions used to manipulate lists of threshold specifications (see thrlist.c)
*/

thrlist thrlist_new (void);
int thrlist_add (thrlist *plist, const char* def);
int thrlist_cmp (Real r, thrlist list, unsigned long fieldno1, unsigned long fieldno2);
void thrlist_dispose (thrlist *plist);


/* Shared definitions coming from GNU DIFF

   Copyright (C) 1988, 1989, 1991, 1992, 1993, 1994, 1995, 1998, 2001,
   2002 Free Software Foundation, Inc.

   This file is part of GNU DIFF.

   GNU DIFF is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2, or (at your option)
   any later version.

   GNU DIFF is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; see the file COPYING.
   If not, write to the Free Software Foundation,
   59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.  */

#include <stdio.h>

#define TAB_WIDTH 8

/*
  Added by Ivano Primi, July 31 2008
*/
#define MIN_ATMOST_NCOLS 16
#define DEF_ATMOST_NCOLS 130
#define MAX_ATMOST_NCOLS 512

/* What kind of changes a hunk contains.  */
enum changes
{
  /* No changes: lines common to both files.  */
  UNCHANGED,

  /* Deletes only: lines taken from just the first file.  */
  OLD,

  /* Inserts only: lines taken from just the second file.  */
  NEW,

  /* Both deletes and inserts: a hunk containing both old and new lines.  */
  CHANGED
};

/* Variables for command line options */

#ifndef GDIFF_OPTIONS
# define XTERN extern
#else
# define XTERN
#endif

/* The significance of white space during comparisons.  */
XTERN enum
{
  /* All white space is significant (the default).  */
  IGNORE_NO_WHITE_SPACE,

  /* Ignore changes due to tab expansion.  */
  IGNORE_TAB_EXPANSION,

  /* Ignore changes in horizontal white space.  */
  IGNORE_SPACE_CHANGE,

  /* Ignore all horizontal white space.  */
  IGNORE_ALL_SPACE
} ignore_white_space;

/* Treat both files as binary files (meaningful only under Doz/Windoz) */
XTERN bool binary;

/* Nonzero means to not show common lines.  */
XTERN bool suppress_common_lines;

/* Expand tabs in the output so the text lines up properly
   despite the characters added to the front of each line (-T).  */
XTERN bool expand_tabs;

/* The half line width and column 2 offset for OUTPUT_SDIFF.  */
XTERN unsigned int sdiff_half_width;
XTERN unsigned int sdiff_column2_offset;

/* Use heuristics for better speed with large files with a small
   density of changes.  */
XTERN bool speed_large_files;

/* Name of program the user invoked (for error messages).  */
XTERN char *program_name;

/* The result of comparison is an "edit script": a chain of `struct change'.
   Each `struct change' represents one place where some lines are deleted
   and some are inserted.

   LINE0 and LINE1 are the first affected lines in the two files (origin 0).
   DELETED is the number of lines deleted here from file 0.
   INSERTED is the number of lines inserted here in file 1.

   If DELETED is 0 then LINE0 is the number of the line before
   which the insertion was done; vice versa for INSERTED and LINE1.  */

struct change
{
  struct change *link;		/* Previous or next edit command  */
  lin inserted;			/* # lines of file 1 changed here.  */
  lin deleted;			/* # lines of file 0 changed here.  */
  lin line0;			/* Line number of 1st deleted line.  */
  lin line1;			/* Line number of 1st inserted line.  */
  bool ignore;			/* Flag used in context.c.  */
};

/* Structures that describe the input files.  */

/* Data on one input file being compared.  */

struct file_data {
    int             desc;	/* File descriptor  */
    char const      *name;	/* File name  */
    struct stat     stat;	/* File status */

    /* Buffer in which text of file is read.  */
    word *buffer;

    /* Allocated size of buffer, in bytes.  Always a multiple of
       sizeof *buffer.  */
    size_t bufsize;

    /* Number of valid bytes now in the buffer.  */
    size_t buffered;

    /* Array of pointers to lines in the file.  */
    char const **linbuf;

    /* linbuf_base <= buffered_lines <= valid_lines <= alloc_lines.
       linebuf[linbuf_base ... buffered_lines - 1] are possibly differing.
       linebuf[linbuf_base ... valid_lines - 1] contain valid data.
       linebuf[linbuf_base ... alloc_lines - 1] are allocated.  */
    lin linbuf_base, buffered_lines, valid_lines, alloc_lines;

    /* Pointer to end of prefix of this file to ignore when hashing.  */
    char const *prefix_end;

    /* Count of lines in the prefix.
       There are this many lines in the file before linbuf[0].  */
    lin prefix_lines;

    /* Pointer to start of suffix of this file to ignore when hashing.  */
    char const *suffix_begin;

    /* Vector, indexed by line number, containing an equivalence code for
       each line.  It is this vector that is actually compared with that
       of another file to generate differences.  */
    lin *equivs;

    /* Vector, like the previous one except that
       the elements for discarded lines have been squeezed out.  */
    lin *undiscarded;

    /* Vector mapping virtual line numbers (not counting discarded lines)
       to real ones (counting those lines).  Both are origin-0.  */
    lin *realindexes;

    /* Total number of nondiscarded lines.  */
    lin nondiscarded_lines;

    /* Vector, indexed by real origin-0 line number,
       containing TRUE for a line that is an insertion or a deletion.
       The results of comparison are stored here.  */
    bool *changed;

    /* 1 if file ends in a line with no final newline.  */
    bool missing_newline;

    /* 1 if at end of file.  */
    bool eof;

    /* 1 more than the maximum equivalence value used for this or its
       sibling file.  */
    lin equiv_max;
};

/* The file buffer, considered as an array of bytes rather than
   as an array of words.  */
#define FILE_BUFFER(f) ((char *) (f)->buffer)

/* Describe the two files currently being compared.  */

XTERN struct file_data files[2];

/* Stdio stream to output diffs to.  */

#define outfile stdout

/* Declare various functions.  */

/* analyze.c */
int diff_2_files (struct file_data[], argslist*);

/* inout.c */
bool read_files (struct file_data[], argslist*);

/* numutil.c */
char* acxnum (const char *str, const struct numfmt* pnf);
int compare_numeric_strings (const char *str1, const struct numfmt* pnf1,
			     const char *str2, const struct numfmt* pnf2);
char* hcxnum (const char *str, const struct numfmt* pnf, hash_value *ph);
#ifdef USE_GMP
int mpf_a2num (Real* pr, const char *q, char** endptr, const struct numfmt* pnf);
#endif /* USE_GMP */

/* side.c */
void print_sdiff_script (struct change *);
void print_1overview_line (const char *left, int are_different, const char *right);

/* util.c */
bool lines_differ (char const *, char const *, int, int, argslist*);
void *zalloc (size_t);

#define stralloc(length) zalloc ((length)+1)

enum changes analyze_hunk (struct change *, lin *, lin *, lin *, lin *);
#ifdef _DEBUG_SCRIPT_
void debug_script (struct change *);
#endif
void perror_with_name (char const *);
void pfatal_with_name (char const *) __attribute__((noreturn));
void print_script (struct change *, void (*) (struct change *));

/* flags.c */
/* This functions were added by Ivano Primi, 14-02-08 */

int init_flags (void);
int print_flags (FILE* fp);
flg_array copy_of_intflagtab (void);
void erase_flags (void);
void notedown_sdiff_script (struct change *script);

/* End Section "Shared definitions coming from GNU DIFF" */

#endif /* _NUMDIFF_H_ */
