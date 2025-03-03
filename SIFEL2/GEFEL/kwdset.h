#ifndef KWDSET_H
#define KWDSET_H

struct enumstr
{
  const char *alias; /// string with enum alias
  int id;            /// enum alias id
};

class kwdset
{
  public :
    kwdset(const long num, const enumstr *setptr);

    /// returns string with alias for the given enum integer value
    const char* get_str(const int i) const;

    /// returns integer value for the given string alias
    int get_id(const char *) const;

    /// returns integer value for the given string alias
    int get_id(char *) const;

    /// checks the val whether is defined in the kwdset
    long check_int(int val) const;

    /// returns index of val in the array of enum set
    long get_ordid(int val) const;

    long n;             /// number of keywords in set
    const enumstr *set; /// array of enumstr
};

#endif
