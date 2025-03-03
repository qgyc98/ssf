template <class type> class QList
{
  int num, limit, ptr, autoDelete;
  type **data;

  public:
  QList()
  {
    autoDelete = 0;
    num = 0;
    limit = 20;
    ptr = -1;
    data = new type*[20];
    for(int i = 0; i < limit; i++)
      data[i] = NULL;
  }

  void append(type* p)
  {
    if (num == limit)
    {
      limit += 20;
      type **n_data = new type*[limit];
      for (int i = 0; i < num; i++)
        n_data[i] = data[i];
      delete data;
      data = n_data;
    }
    data[num++] = p;
  }

  type* at(int i)
  {
    if (i < 0 || i >= num)
      return NULL;
    return data[i];
  }

  int at()
  {
    return ptr;
  }

  type* first()
  {
    return data[ptr = 0];
  }

  type* next()
  {
    return data[++ptr];
  }

  void setAutoDelete(int d)
  {
    autoDelete = d;
  }

  void clear()
  {
    if (autoDelete)
    {
      for (int i = 0; i < num; i++)
        delete data[i];
    }
    delete [] data;
  }

  int count()
  {
    return num;
  }
};
