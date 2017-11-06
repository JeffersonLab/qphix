#pragma once

#include "twisted_mass_enum.h"

#include <cstdio>
#include <cstdlib>

using namespace std;

#include "instructions.h"
#include "data_types.h"

typedef struct {
  const char *name;
  int s[2][2];
  void (*CVecFunc[2])(InstVector &ivector,
                      FVec *r,
                      FVec *s1,
                      FVec *s2,
                      string &mask);
  //	void (*CVecFunc2)(InstVector& ivector, FVec *r, FVec *s1,
  //FVec *s2, string &mask);
} proj_ops;

typedef struct {
  const char *name;
  int s2;
  int s3;
  void (*CVecFuncTop2)(InstVector &ivector,
                       FVec *r,
                       FVec *s1,
                       FVec *s2,
                       string &mask);
  void (*CVecFunc1)(InstVector &ivector,
                    FVec *r,
                    FVec *s1,
                    FVec *s2,
                    string &mask);
  void (*CVecFunc2)(InstVector &ivector,
                    FVec *r,
                    FVec *s1,
                    FVec *s2,
                    string &mask);
} recons_ops;

// Uses psi[][] as temp and b_spinor[][] as implicit output
void project(InstVector &ivector,
             string base,
             string offset,
             proj_ops &ops,
             bool isFace,
             string mask,
             int dir);
// Serial Spin version
void project(InstVector &ivector,
             string base,
             string offset,
             proj_ops &ops,
             bool isFace,
             string mask,
             int dir,
             int s);
void recons_add(InstVector &ivector,
                recons_ops &ops,
                FVec outspinor[4][3][2],
                string &mask);
// Serial Spin version
void recons_add(InstVector &ivector,
                recons_ops &ops,
                FVec outspinor[4][3][2],
                string &mask,
                int s);
void zeroResult(InstVector &ivector, FVec *outspinor);
void two_flav_tm_inverse_clover_term(InstVector &ivector,
                                     string _mask = "");
void clover_term(InstVector &ivector,
                 FVec in_spinor[4][3][2],
                 bool face,
                 string _mask = "");
void full_clover_term(InstVector &ivector,
                      FVec in_spinor[4][3][2],
                      bool face,
                      string _mask = "");
void twisted_term(InstVector &ivector, bool isPlus);
void inverse_twisted_term(InstVector &ivector,
                          FVec in_spinor[4][3][2],
                          bool face,
                          bool isPlus,
                          string _mask = "");
void applyTwistedBoundaryConditions(InstVector &ivector,
                                    bool const adjMul,
                                    bool const has_tbc);
// void achiResult(InstVector& ivector, bool clover);

void achiResult(InstVector &ivector,
                bool const clover,
                TwistedMassVariant const twisted_mass,
                bool const isPlus);

void loadGaugeDir(InstVector &ivector, int dir, bool compress12);
void matMultVec(InstVector &ivector, bool adjMul, int s);
void matMultVec(InstVector &ivector, bool adjMul);

void dslash_plain_body(InstVector &ivector,
                       bool const compress12,
                       bool const clover,
                       TwistedMassVariant const twisted_mass,
                       bool const isPlus,
                       bool const *const tbc);

// ***** ------- a chi - b D psi versions

void dslash_achimbdpsi_body(InstVector &ivector,
                            bool const compress12,
                            bool const clover,
                            TwistedMassVariant const twisted_mass,
                            bool const isPlus,
                            bool const *const tbc);

void pack_face_to_dir_dim_vec(InstVector &ivector,
                              bool isPlus,
                              int dir,
                              int dim);
void recons_add_face_from_dir_dim_vec(InstVector &ivector,
                                      bool compress12,
                                      bool isPlus,
                                      int dir,
                                      int dim,
                                      bool clover,
                                      bool twisted_mass,
                                      bool const use_tbc);

void dslash_body(InstVector &ivector,
                 bool compress12,
                 proj_ops *ops,
                 recons_ops *rec_ops_bw,
                 recons_ops *rec_ops_fw,
                 FVec outspinor[4][3][2],
                 bool const *const tbc);

void pack_face_vec(InstVector &ivector,
                   FVec spinor[2][3][2],
                   proj_ops proj[],
                   int dir);
void recons_add_face_vec(InstVector &ivector,
                         bool compress12,
                         bool adjMul,
                         recons_ops rops[],
                         int dir,
                         int dim,
                         bool clover,
                         bool twisted_mass,
                         bool isPlus,
                         bool const use_tbc);
