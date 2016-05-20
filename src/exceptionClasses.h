#ifndef MPScode_EXCEPTION_CLASSES
#define MPScode_EXCEPTION_CLASSES

#include <exception>

class empty_table{
 public:
 explicit empty_table(int site_arg):position(site_arg){}
  int site() const {return position;}
 private:
  int position;
};

#endif
