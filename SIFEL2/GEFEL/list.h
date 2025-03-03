#ifndef LIST_H
#define LIST_H

/**
  The class creates list of arbitrary items. Pointers to
  particular items are stored in the array of pointers data

  Created by TKo, 09.2010 
*/
class list
{
  /// number of items
  long num;
  /// limit of array data
  long limit;
  /// array with pointers to particular items of the list
  void **data;

  public:
   list();
   ~list();
   void append(void* p);
   void* at(long i);
   long  count();
//   void  set_act(long i);
//   long  get_act();
};

#endif
