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

class svd_failure{
 public:
 svd_failure(int site_arg, int dir_arg):position(site_arg),dir(dir_arg){}
  int site() const {return position;}
  int direction() const {return dir;}
 private:
  int position;
  int dir;
};

class critical_error{
 public: 
  critical_error(){}
};

#endif
