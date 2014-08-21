#ifndef SPINOR_H
#define SPINOR_H

#include "cpp_dslash_scalar.h"
#include "shift_table_scalar.h"


// A type for a Blocked Spinor
template<typename T, int V>
class SpinorCB {
public:

  // Traits
  typedef T  BaseType;
  typedef T  FourSpinor[4][3][2][V];
  typedef T  TwoSpinor[2][3][2][V];

  Spinor(const ShiftTable& s) {


  }; // constructor 

  ~Spinor() { 
    // ALIGNED_FREE(data);

  } // Destructor
 
  FourSpinor operator[](const int i) { return data[i]; }
  const FourSpinor operator[](const int i) const { return data[i]; }

private: 
  FourSpinor* data;

};




#endif
