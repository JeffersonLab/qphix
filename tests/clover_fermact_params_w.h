// -*- C++ -*-

#ifndef CLOVER_FERMACT_PARAMS_W_H
#define CLOVER_FERMACT_PARAMS_W_H

#include <qdp.h>

// These are parameters for setting up a clover term
// I nicked them from chroma, and cut out all the XML
// IO stuff

//! Parameters for anisotropy
struct AnisoParam_t {
  AnisoParam_t()
  {
    anisoP = false; // Isotropic by default
    t_dir = 3; // direction 3 is time (or anisotropic direction anyhow)
    xi_0 = QDP::Real(1); // No anisotropy
    nu = QDP::Real(1); // No anisotropy
  }
  ~AnisoParam_t() {}

  bool anisoP;
  int t_dir;
  QDP::Real xi_0;
  QDP::Real nu;
};

//! Params for clover ferm acts
/*! \ingroup fermacts */
struct CloverFermActParams {
  CloverFermActParams()
  {
    Mass = QDP::Real(0); // Default zero mass
    clovCoeffR = QDP::Real(0);
    clovCoeffT = QDP::Real(0);
    u0 = QDP::Real(1);
  }

  QDP::Real Mass;
  QDP::Real clovCoeffR;
  QDP::Real clovCoeffT;
  QDP::Real u0;

  // Optional Anisotropy
  AnisoParam_t anisoParam;
};

#endif
