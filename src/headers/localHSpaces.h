#ifndef LOCAL_HILBERT_SPACE_DIMENSION_LIST
#define LOCAL_HILBERT_SPACE_DIMENSION_LIST

#include <vector>

class localHSpaces{
 public:
  localHSpaces();
  localHSpaces(int d);
  localHSpaces(std::vector<int> const &d);
  int locd(int const i) const;
  int maxd() const;
 private:
  int constantDimension;
  std::vector<int> localHDims;
};

#endif
