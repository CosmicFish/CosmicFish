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

#define I_DEF_SEP ' '
#define DEF_IFS   "\\s \\t \\n"
#define NEWLINE_STR "\n"

char** ssplit (const char* str, char separator);
char** ssplit_former_way (const char* str);
void print_string_vector (FILE* fp, const char** sv, char separator);
void sort_string_vector (char** sv);
void remove_duplicates_from_string_vector (char** sv);
int is_string_in_vector (const char* str, const char** sv);
size_t is_char_in_vector (int ch, const char** sv);
void delete_string_vector (char** sv);

char* get_separating_string (const char* istr);
char* string_cspn (const char* str, const char** svi, int eos);
char* string_spn (const char* str, const char** svi, int eos);
