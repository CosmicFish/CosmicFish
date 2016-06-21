/* Failure exit status

   Copyright (C) 2002, 2003, 2005, 2006, 2007 Free Software Foundation, Inc.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; see the file COPYING.
   If not, write to the Free Software Foundation,
   51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.  */

#include "system.h"

#include "exitfail.h"

#include <stdlib.h>

/*
  Modified by Ivano Primi on April 18, 2009.
  The value of the constant EXIT_FAILURE (1)
  is used by Numdiff to mean that the given files
  are different.
  To signal troubles Numdiff uses EXIT_TROUBLE (-1).  
*/

int volatile exit_failure = EXIT_TROUBLE;
