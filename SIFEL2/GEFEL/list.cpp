#include "list.h"
#include "stdlib.h"



/**
  The constructor initializes all data members

  Created by TKo, 09.2010 
*/
list::list()
{
  num = 0;
  limit = 20;
  data = new void*[limit];
  for(long i = 0; i < limit; i++)
    data[i] = NULL;
}



/**
  The destructor releases memory used by data array.
  Particular elements of the array are NOT deleted.

  Created by TKo, 09.2010 
*/
list::~list()
{
  delete [] data;
}



/**
  The function adds new item to the list. If it is necessary,
  the array data is reallocated. The new item is given by the 
  pointer p.

  @param p - pointer to the added object

  @returns The function does not return anything.

  Created by TKo, 09.2010
*/
void list::append(void* p)
{
  if (num == limit)
  {
    limit += 20;
    void **n_data = new void*[limit];
    for (long i = 0; i < num; i++)
      n_data[i] = data[i];
    delete [] data;
    data = n_data;
  }
  data[num] = p;
  num++;
}



/**
  The function returns i-th item from the list.

  @param i - index of the returned item

  @return The function returns required i-th item.

  Created by TKo, 09.2010
*/
void* list::at(long i)
{
  if (i < 0 || i >= num)
    return NULL;
  return data[i];
}



/**
  The function returns actual number of inserted items.

  @return The function returns number of items

  Created by TKo, 09.2010
*/
long list::count()
{
  return num;
}
