#include "kwdset.h"
#include "iotools.h"
#include <string.h>

kwdset::kwdset(const long num, const enumstr *setptr)
{
  n = num;
  set = setptr;
}



/**
  The function returns string with the enum alias representation of the given enum integer value i.

  @param i - integer value of enum that would be converted to string representation

  Returns:
   @returns The function returns a constant string representation of the given enum integer value.

  Created by Tomas Koudelka, 4.7.2014
*/
const char* kwdset::get_str(const int i) const
{
  long j;
  for(j=0; j < n; j++)
  {
    if (set[j].id == i)
      return (set[j].alias);
  }
  print_err("Cannot convert enum integer value %d to alias string", __FILE__, __LINE__, __func__, i);
  return NULL;
}



/**
  The function returns enum integer value for the given constant string with the enum alias.

  @param s - constant string with the enum alias name

  Returns:
   @returns The function returns enum integer value for the given string with enum alias.

  Created by Tomas Koudelka, 4.7.2014
*/
int kwdset::get_id(const char *s) const
{
  long j;
  for(j=0; j < n; j++)
  {
    if (strcmp(s, set[j].alias) == 0)
      return set[j].id;
  }
  print_err("Cannot convert string '%s' to enum integer value", __FILE__, __LINE__, __func__, s);
  return -1;
}



/**
  The function returns enum integer value for the given string with the enum alias.

  @param s - string with the enum alias name

  Returns:
   @returns The function returns enum integer value for the given string with enum alias.

  Created by Tomas Koudelka, 4.7.2014
*/
int kwdset::get_id(char *s) const
{
  long j;
  for(j=0; j < n; j++)
  {
    if (strcmp(s, set[j].alias) == 0)
      return set[j].id;
  }
  print_err("Cannot convert string '%s' to enum integer value", __FILE__, __LINE__, __func__, s);
  return -1;
}



/**
  The function checks whether the val is defined in the kwdset

  @param val - tested enum value

  @retval 0 - on success - the val was found in the given kwdset
  @retval 1 - on error - the val is not defined in the given enum

  Created by Tomas Koudelka, 11.8.2014
*/
long kwdset::check_int(int val) const
{
  long j;

  for(j=0; j < n; j++)
  {
    if (set[j].id == val)
      return 0;
  }
  return 1;
}



/**
  The function returns order of the given enum in the set.

  @param val - tested enum value

  @return The function returns index of the given val in the enum set array.

  Created by Tomas Koudelka, 11.8.2014
*/
long kwdset::get_ordid(int val) const
{
  long j;

  for(j=0; j < n; j++)
  {
    if (set[j].id == val)
      return j;
  }
  return -1;
}
