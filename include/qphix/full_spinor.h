/*
 * full_spinor.h
 *
 *  Created on: Oct 10, 2017
 *      Author: bjoo
 */

#pragma once
#ifndef INCLUDE_QPHIX_FULL_SPINOR_H_
#define INCLUDE_QPHIX_FULL_SPINOR_H_


#include <qphix/geometry.h>
namespace QPhiX {

  template<typename T, int V, int S, bool compress12>
  class FullSpinor {
  public:
    using CBSpinor = typename FourSpinorHandle<T,V,S,compress12>::ValueType;

    FullSpinor(Geometry<T,V,S, compress12>& geom) : _geom(geom),
        _data_cb0(FourSpinorHandle<T,V,S,compress12>(geom)),
        _data_cb1(FourSpinorHandle<T,V,S,compress12>(geom))
      {
          _data[0]=&_data_cb0;
          _data[1]=&_data_cb1;
      }


    // Destructor is automatic.

    FourSpinorHandle<T,V,S,compress12>& getCB(int cb) {
      return *(_data[cb]);
    }

    const FourSpinorHandle<T,V,S,compress12>& getCB(int cb) const {
      return *(_data[cb]);
    }

    // Dunno if we will need this:
    CBSpinor* getCBData(int cb) const {
      return _data[cb]->get();
    }

    Geometry<T,V,S,compress12>& getGeometry(void) {
      return _geom;
    }

  private:
    Geometry<T,V,S,compress12>& _geom;
    FourSpinorHandle<T,V,S,compress12> _data_cb0;
    FourSpinorHandle<T,V,S,compress12> _data_cb1;
    FourSpinorHandle<T,V,S,compress12>* _data[2];
  };

}




#endif /* INCLUDE_QPHIX_FULL_SPINOR_H_ */
