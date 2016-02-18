#ifndef LOCAL_HILBERT_SPACE_DIMENSION_LIST
#define LOCAL_HILBERT_SPACE_DIMENSION_LIST

#include <vector>

class localHSpaces{
 public:
  localHSpaces();
  localHSpaces(int const d);
  localHSpaces(std::vector<int> d);
  int locd(int const i) const;
  int maxd() const;
 private:
  int constantDimension;
  std::vector<int> localHDims;
};

#endif
