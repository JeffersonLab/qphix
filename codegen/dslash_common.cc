#include <cstdio>
#include <cstdlib>
using namespace std;

#include "dslash.h"

#include "unsupported_values.h"

extern string beta_names[8];

extern string alpha_name;
extern string outBase;
extern string outOffs;
extern string gBase;
extern string gOffs;
extern string chiBase;
extern string chiOffs;
extern string clBase;
extern string clOffs;

extern string mu_name;
extern string mu_inv_name;

// kernel parameter names for the two-flavour kernel(s)
extern string chi2Base;
extern string out2Base;
extern string fclBase;
extern string fcl2Base;

extern string prec_mass_rho_name;

FVec b_S0_C0_RE("b_S0_C0_RE");
FVec b_S0_C0_IM("b_S0_C0_IM");
FVec b_S0_C1_RE("b_S0_C1_RE");
FVec b_S0_C1_IM("b_S0_C1_IM");
FVec b_S0_C2_RE("b_S0_C2_RE");
FVec b_S0_C2_IM("b_S0_C2_IM");
FVec b_S1_C0_RE("b_S1_C0_RE");
FVec b_S1_C0_IM("b_S1_C0_IM");
FVec b_S1_C1_RE("b_S1_C1_RE");
FVec b_S1_C1_IM("b_S1_C1_IM");
FVec b_S1_C2_RE("b_S1_C2_RE");
FVec b_S1_C2_IM("b_S1_C2_IM");

FVec b_spinor[2][3][2] = {{{b_S0_C0_RE, b_S0_C0_IM},
                           {b_S0_C1_RE, b_S0_C1_IM},
                           {b_S0_C2_RE, b_S0_C2_IM}},
                          {{b_S1_C0_RE, b_S1_C0_IM},
                           {b_S1_C1_RE, b_S1_C1_IM},
                           {b_S1_C2_RE, b_S1_C2_IM}}};

FVec ub_S0_C0_RE("ub_S0_C0_RE");
FVec ub_S0_C0_IM("ub_S0_C0_IM");
FVec ub_S0_C1_RE("ub_S0_C1_RE");
FVec ub_S0_C1_IM("ub_S0_C1_IM");
FVec ub_S0_C2_RE("ub_S0_C2_RE");
FVec ub_S0_C2_IM("ub_S0_C2_IM");
FVec ub_S1_C0_RE("ub_S1_C0_RE");
FVec ub_S1_C0_IM("ub_S1_C0_IM");
FVec ub_S1_C1_RE("ub_S1_C1_RE");
FVec ub_S1_C1_IM("ub_S1_C1_IM");
FVec ub_S1_C2_RE("ub_S1_C2_RE");
FVec ub_S1_C2_IM("ub_S1_C2_IM");

FVec ub_spinor[2][3][2] = {{{ub_S0_C0_RE, ub_S0_C0_IM},
                            {ub_S0_C1_RE, ub_S0_C1_IM},
                            {ub_S0_C2_RE, ub_S0_C2_IM}},
                           {{ub_S1_C0_RE, ub_S1_C0_IM},
                            {ub_S1_C1_RE, ub_S1_C1_IM},
                            {ub_S1_C2_RE, ub_S1_C2_IM}}};

FVec out_S0_C0_RE("out_S0_C0_RE");
FVec out_S0_C0_IM("out_S0_C0_IM");
FVec out_S0_C1_RE("out_S0_C1_RE");
FVec out_S0_C1_IM("out_S0_C1_IM");
FVec out_S0_C2_RE("out_S0_C2_RE");
FVec out_S0_C2_IM("out_S0_C2_IM");
FVec out_S1_C0_RE("out_S1_C0_RE");
FVec out_S1_C0_IM("out_S1_C0_IM");
FVec out_S1_C1_RE("out_S1_C1_RE");
FVec out_S1_C1_IM("out_S1_C1_IM");
FVec out_S1_C2_RE("out_S1_C2_RE");
FVec out_S1_C2_IM("out_S1_C2_IM");
FVec out_S2_C0_RE("out_S2_C0_RE");
FVec out_S2_C0_IM("out_S2_C0_IM");
FVec out_S2_C1_RE("out_S2_C1_RE");
FVec out_S2_C1_IM("out_S2_C1_IM");
FVec out_S2_C2_RE("out_S2_C2_RE");
FVec out_S2_C2_IM("out_S2_C2_IM");
FVec out_S3_C0_RE("out_S3_C0_RE");
FVec out_S3_C0_IM("out_S3_C0_IM");
FVec out_S3_C1_RE("out_S3_C1_RE");
FVec out_S3_C1_IM("out_S3_C1_IM");
FVec out_S3_C2_RE("out_S3_C2_RE");
FVec out_S3_C2_IM("out_S3_C2_IM");

FVec out_spinor[4][3][2] = {{{out_S0_C0_RE, out_S0_C0_IM},
                             {out_S0_C1_RE, out_S0_C1_IM},
                             {out_S0_C2_RE, out_S0_C2_IM}},
                            {{out_S1_C0_RE, out_S1_C0_IM},
                             {out_S1_C1_RE, out_S1_C1_IM},
                             {out_S1_C2_RE, out_S1_C2_IM}},
                            {{out_S2_C0_RE, out_S2_C0_IM},
                             {out_S2_C1_RE, out_S2_C1_IM},
                             {out_S2_C2_RE, out_S2_C2_IM}},
                            {{out_S3_C0_RE, out_S3_C0_IM},
                             {out_S3_C1_RE, out_S3_C1_IM},
                             {out_S3_C2_RE, out_S3_C2_IM}}};

// for two-flavour kernels we have a second output spinor
FVec out2_S0_C0_RE("out2_S0_C0_RE");
FVec out2_S0_C0_IM("out2_S0_C0_IM");
FVec out2_S0_C1_RE("out2_S0_C1_RE");
FVec out2_S0_C1_IM("out2_S0_C1_IM");
FVec out2_S0_C2_RE("out2_S0_C2_RE");
FVec out2_S0_C2_IM("out2_S0_C2_IM");
FVec out2_S1_C0_RE("out2_S1_C0_RE");
FVec out2_S1_C0_IM("out2_S1_C0_IM");
FVec out2_S1_C1_RE("out2_S1_C1_RE");
FVec out2_S1_C1_IM("out2_S1_C1_IM");
FVec out2_S1_C2_RE("out2_S1_C2_RE");
FVec out2_S1_C2_IM("out2_S1_C2_IM");
FVec out2_S2_C0_RE("out2_S2_C0_RE");
FVec out2_S2_C0_IM("out2_S2_C0_IM");
FVec out2_S2_C1_RE("out2_S2_C1_RE");
FVec out2_S2_C1_IM("out2_S2_C1_IM");
FVec out2_S2_C2_RE("out2_S2_C2_RE");
FVec out2_S2_C2_IM("out2_S2_C2_IM");
FVec out2_S3_C0_RE("out2_S3_C0_RE");
FVec out2_S3_C0_IM("out2_S3_C0_IM");
FVec out2_S3_C1_RE("out2_S3_C1_RE");
FVec out2_S3_C1_IM("out2_S3_C1_IM");
FVec out2_S3_C2_RE("out2_S3_C2_RE");
FVec out2_S3_C2_IM("out2_S3_C2_IM");

FVec out2_spinor[4][3][2] = {{{out2_S0_C0_RE, out2_S0_C0_IM},
                              {out2_S0_C1_RE, out2_S0_C1_IM},
                              {out2_S0_C2_RE, out2_S0_C2_IM}},
                             {{out2_S1_C0_RE, out2_S1_C0_IM},
                              {out2_S1_C1_RE, out2_S1_C1_IM},
                              {out2_S1_C2_RE, out2_S1_C2_IM}},
                             {{out2_S2_C0_RE, out2_S2_C0_IM},
                              {out2_S2_C1_RE, out2_S2_C1_IM},
                              {out2_S2_C2_RE, out2_S2_C2_IM}},
                             {{out2_S3_C0_RE, out2_S3_C0_IM},
                              {out2_S3_C1_RE, out2_S3_C1_IM},
                              {out2_S3_C2_RE, out2_S3_C2_IM}}};

FVec clout_spinor[2][6][2] = {{{out_S0_C0_RE, out_S0_C0_IM},
                               {out_S0_C1_RE, out_S0_C1_IM},
                               {out_S0_C2_RE, out_S0_C2_IM},
                               {out_S1_C0_RE, out_S1_C0_IM},
                               {out_S1_C1_RE, out_S1_C1_IM},
                               {out_S1_C2_RE, out_S1_C2_IM}},
                              {{out_S2_C0_RE, out_S2_C0_IM},
                               {out_S2_C1_RE, out_S2_C1_IM},
                               {out_S2_C2_RE, out_S2_C2_IM},
                               {out_S3_C0_RE, out_S3_C0_IM},
                               {out_S3_C1_RE, out_S3_C1_IM},
                               {out_S3_C2_RE, out_S3_C2_IM}}};

FVec chi_S0_C0_RE("chi_S0_C0_RE");
FVec chi_S0_C0_IM("chi_S0_C0_IM");
FVec chi_S0_C1_RE("chi_S0_C1_RE");
FVec chi_S0_C1_IM("chi_S0_C1_IM");
FVec chi_S0_C2_RE("chi_S0_C2_RE");
FVec chi_S0_C2_IM("chi_S0_C2_IM");
FVec chi_S1_C0_RE("chi_S1_C0_RE");
FVec chi_S1_C0_IM("chi_S1_C0_IM");
FVec chi_S1_C1_RE("chi_S1_C1_RE");
FVec chi_S1_C1_IM("chi_S1_C1_IM");
FVec chi_S1_C2_RE("chi_S1_C2_RE");
FVec chi_S1_C2_IM("chi_S1_C2_IM");
FVec chi_S2_C0_RE("chi_S2_C0_RE");
FVec chi_S2_C0_IM("chi_S2_C0_IM");
FVec chi_S2_C1_RE("chi_S2_C1_RE");
FVec chi_S2_C1_IM("chi_S2_C1_IM");
FVec chi_S2_C2_RE("chi_S2_C2_RE");
FVec chi_S2_C2_IM("chi_S2_C2_IM");
FVec chi_S3_C0_RE("chi_S3_C0_RE");
FVec chi_S3_C0_IM("chi_S3_C0_IM");
FVec chi_S3_C1_RE("chi_S3_C1_RE");
FVec chi_S3_C1_IM("chi_S3_C1_IM");
FVec chi_S3_C2_RE("chi_S3_C2_RE");
FVec chi_S3_C2_IM("chi_S3_C2_IM");

FVec chi_spinor[4][3][2] = {{{chi_S0_C0_RE, chi_S0_C0_IM},
                             {chi_S0_C1_RE, chi_S0_C1_IM},
                             {chi_S0_C2_RE, chi_S0_C2_IM}},
                            {{chi_S1_C0_RE, chi_S1_C0_IM},
                             {chi_S1_C1_RE, chi_S1_C1_IM},
                             {chi_S1_C2_RE, chi_S1_C2_IM}},
                            {{chi_S2_C0_RE, chi_S2_C0_IM},
                             {chi_S2_C1_RE, chi_S2_C1_IM},
                             {chi_S2_C2_RE, chi_S2_C2_IM}},
                            {{chi_S3_C0_RE, chi_S3_C0_IM},
                             {chi_S3_C1_RE, chi_S3_C1_IM},
                             {chi_S3_C2_RE, chi_S3_C2_IM}}};

// for two-flavour kernels, we have a second input spinor
FVec chi2_S0_C0_RE("chi2_S0_C0_RE");
FVec chi2_S0_C0_IM("chi2_S0_C0_IM");
FVec chi2_S0_C1_RE("chi2_S0_C1_RE");
FVec chi2_S0_C1_IM("chi2_S0_C1_IM");
FVec chi2_S0_C2_RE("chi2_S0_C2_RE");
FVec chi2_S0_C2_IM("chi2_S0_C2_IM");
FVec chi2_S1_C0_RE("chi2_S1_C0_RE");
FVec chi2_S1_C0_IM("chi2_S1_C0_IM");
FVec chi2_S1_C1_RE("chi2_S1_C1_RE");
FVec chi2_S1_C1_IM("chi2_S1_C1_IM");
FVec chi2_S1_C2_RE("chi2_S1_C2_RE");
FVec chi2_S1_C2_IM("chi2_S1_C2_IM");
FVec chi2_S2_C0_RE("chi2_S2_C0_RE");
FVec chi2_S2_C0_IM("chi2_S2_C0_IM");
FVec chi2_S2_C1_RE("chi2_S2_C1_RE");
FVec chi2_S2_C1_IM("chi2_S2_C1_IM");
FVec chi2_S2_C2_RE("chi2_S2_C2_RE");
FVec chi2_S2_C2_IM("chi2_S2_C2_IM");
FVec chi2_S3_C0_RE("chi2_S3_C0_RE");
FVec chi2_S3_C0_IM("chi2_S3_C0_IM");
FVec chi2_S3_C1_RE("chi2_S3_C1_RE");
FVec chi2_S3_C1_IM("chi2_S3_C1_IM");
FVec chi2_S3_C2_RE("chi2_S3_C2_RE");
FVec chi2_S3_C2_IM("chi2_S3_C2_IM");

FVec chi2_spinor[4][3][2] = {{{chi2_S0_C0_RE, chi2_S0_C0_IM},
                              {chi2_S0_C1_RE, chi2_S0_C1_IM},
                              {chi2_S0_C2_RE, chi2_S0_C2_IM}},
                             {{chi2_S1_C0_RE, chi2_S1_C0_IM},
                              {chi2_S1_C1_RE, chi2_S1_C1_IM},
                              {chi2_S1_C2_RE, chi2_S1_C2_IM}},
                             {{chi2_S2_C0_RE, chi2_S2_C0_IM},
                              {chi2_S2_C1_RE, chi2_S2_C1_IM},
                              {chi2_S2_C2_RE, chi2_S2_C2_IM}},
                             {{chi2_S3_C0_RE, chi2_S3_C0_IM},
                              {chi2_S3_C1_RE, chi2_S3_C1_IM},
                              {chi2_S3_C2_RE, chi2_S3_C2_IM}}};

FVec dout_S0_C0_RE("dout_S0_C0_RE");
FVec dout_S0_C0_IM("dout_S0_C0_IM");
FVec dout_S0_C1_RE("dout_S0_C1_RE");
FVec dout_S0_C1_IM("dout_S0_C1_IM");
FVec dout_S0_C2_RE("dout_S0_C2_RE");
FVec dout_S0_C2_IM("dout_S0_C2_IM");
FVec dout_S1_C0_RE("dout_S1_C0_RE");
FVec dout_S1_C0_IM("dout_S1_C0_IM");
FVec dout_S1_C1_RE("dout_S1_C1_RE");
FVec dout_S1_C1_IM("dout_S1_C1_IM");
FVec dout_S1_C2_RE("dout_S1_C2_RE");
FVec dout_S1_C2_IM("dout_S1_C2_IM");
FVec dout_S2_C0_RE("dout_S2_C0_RE");
FVec dout_S2_C0_IM("dout_S2_C0_IM");
FVec dout_S2_C1_RE("dout_S2_C1_RE");
FVec dout_S2_C1_IM("dout_S2_C1_IM");
FVec dout_S2_C2_RE("dout_S2_C2_RE");
FVec dout_S2_C2_IM("dout_S2_C2_IM");
FVec dout_S3_C0_RE("dout_S3_C0_RE");
FVec dout_S3_C0_IM("dout_S3_C0_IM");
FVec dout_S3_C1_RE("dout_S3_C1_RE");
FVec dout_S3_C1_IM("dout_S3_C1_IM");
FVec dout_S3_C2_RE("dout_S3_C2_RE");
FVec dout_S3_C2_IM("dout_S3_C2_IM");

FVec dout_spinor[4][3][2] = {{{dout_S0_C0_RE, dout_S0_C0_IM},
                              {dout_S0_C1_RE, dout_S0_C1_IM},
                              {dout_S0_C2_RE, dout_S0_C2_IM}},
                             {{dout_S1_C0_RE, dout_S1_C0_IM},
                              {dout_S1_C1_RE, dout_S1_C1_IM},
                              {dout_S1_C2_RE, dout_S1_C2_IM}},
                             {{dout_S2_C0_RE, dout_S2_C0_IM},
                              {dout_S2_C1_RE, dout_S2_C1_IM},
                              {dout_S2_C2_RE, dout_S2_C2_IM}},
                             {{dout_S3_C0_RE, dout_S3_C0_IM},
                              {dout_S3_C1_RE, dout_S3_C1_IM},
                              {dout_S3_C2_RE, dout_S3_C2_IM}}};

FVec cl_diag_0("cl_diag_0");
FVec cl_diag_1("cl_diag_1");
FVec cl_diag_2("cl_diag_2");
FVec cl_diag_3("cl_diag_3");
FVec cl_diag_4("cl_diag_4");
FVec cl_diag_5("cl_diag_5");

FVec clov_diag[6] = {cl_diag_0,
                     cl_diag_1,
                     cl_diag_2,
                     cl_diag_3,
                     cl_diag_4,
                     cl_diag_5};

FVec cl_offdiag_0_RE("cl_offdiag_0_RE");
FVec cl_offdiag_0_IM("cl_offdiag_0_IM");
FVec cl_offdiag_1_RE("cl_offdiag_1_RE");
FVec cl_offdiag_1_IM("cl_offdiag_1_IM");
FVec cl_offdiag_2_RE("cl_offdiag_2_RE");
FVec cl_offdiag_2_IM("cl_offdiag_2_IM");
FVec cl_offdiag_3_RE("cl_offdiag_3_RE");
FVec cl_offdiag_3_IM("cl_offdiag_3_IM");
FVec cl_offdiag_4_RE("cl_offdiag_4_RE");
FVec cl_offdiag_4_IM("cl_offdiag_4_IM");
FVec cl_offdiag_5_RE("cl_offdiag_5_RE");
FVec cl_offdiag_5_IM("cl_offdiag_5_IM");
FVec cl_offdiag_6_RE("cl_offdiag_6_RE");
FVec cl_offdiag_6_IM("cl_offdiag_6_IM");
FVec cl_offdiag_7_RE("cl_offdiag_7_RE");
FVec cl_offdiag_7_IM("cl_offdiag_7_IM");
FVec cl_offdiag_8_RE("cl_offdiag_8_RE");
FVec cl_offdiag_8_IM("cl_offdiag_8_IM");
FVec cl_offdiag_9_RE("cl_offdiag_9_RE");
FVec cl_offdiag_9_IM("cl_offdiag_9_IM");
FVec cl_offdiag_10_RE("cl_offdiag_10_RE");
FVec cl_offdiag_10_IM("cl_offdiag_10_IM");
FVec cl_offdiag_11_RE("cl_offdiag_11_RE");
FVec cl_offdiag_11_IM("cl_offdiag_11_IM");
FVec cl_offdiag_12_RE("cl_offdiag_12_RE");
FVec cl_offdiag_12_IM("cl_offdiag_12_IM");
FVec cl_offdiag_13_RE("cl_offdiag_13_RE");
FVec cl_offdiag_13_IM("cl_offdiag_13_IM");
FVec cl_offdiag_14_RE("cl_offdiag_14_RE");
FVec cl_offdiag_14_IM("cl_offdiag_14_IM");

FVec clov_offdiag[15][2] = {{cl_offdiag_0_RE, cl_offdiag_0_IM},
                            {cl_offdiag_1_RE, cl_offdiag_1_IM},
                            {cl_offdiag_2_RE, cl_offdiag_2_IM},
                            {cl_offdiag_3_RE, cl_offdiag_3_IM},
                            {cl_offdiag_4_RE, cl_offdiag_4_IM},
                            {cl_offdiag_5_RE, cl_offdiag_5_IM},
                            {cl_offdiag_6_RE, cl_offdiag_6_IM},
                            {cl_offdiag_7_RE, cl_offdiag_7_IM},
                            {cl_offdiag_8_RE, cl_offdiag_8_IM},
                            {cl_offdiag_9_RE, cl_offdiag_9_IM},
                            {cl_offdiag_10_RE, cl_offdiag_10_IM},
                            {cl_offdiag_11_RE, cl_offdiag_11_IM},
                            {cl_offdiag_12_RE, cl_offdiag_12_IM},
                            {cl_offdiag_13_RE, cl_offdiag_13_IM},
                            {cl_offdiag_14_RE, cl_offdiag_14_IM}};

FVec cl_full_00_RE("cl_full_00_RE");
FVec cl_full_00_IM("cl_full_00_IM");
FVec cl_full_01_RE("cl_full_01_RE");
FVec cl_full_01_IM("cl_full_01_IM");
FVec cl_full_02_RE("cl_full_02_RE");
FVec cl_full_02_IM("cl_full_02_IM");
FVec cl_full_03_RE("cl_full_03_RE");
FVec cl_full_03_IM("cl_full_03_IM");
FVec cl_full_04_RE("cl_full_04_RE");
FVec cl_full_04_IM("cl_full_04_IM");
FVec cl_full_05_RE("cl_full_05_RE");
FVec cl_full_05_IM("cl_full_05_IM");
FVec cl_full_10_RE("cl_full_10_RE");
FVec cl_full_10_IM("cl_full_10_IM");
FVec cl_full_11_RE("cl_full_11_RE");
FVec cl_full_11_IM("cl_full_11_IM");
FVec cl_full_12_RE("cl_full_12_RE");
FVec cl_full_12_IM("cl_full_12_IM");
FVec cl_full_13_RE("cl_full_13_RE");
FVec cl_full_13_IM("cl_full_13_IM");
FVec cl_full_14_RE("cl_full_14_RE");
FVec cl_full_14_IM("cl_full_14_IM");
FVec cl_full_15_RE("cl_full_15_RE");
FVec cl_full_15_IM("cl_full_15_IM");
FVec cl_full_20_RE("cl_full_20_RE");
FVec cl_full_20_IM("cl_full_20_IM");
FVec cl_full_21_RE("cl_full_21_RE");
FVec cl_full_21_IM("cl_full_21_IM");
FVec cl_full_22_RE("cl_full_22_RE");
FVec cl_full_22_IM("cl_full_22_IM");
FVec cl_full_23_RE("cl_full_23_RE");
FVec cl_full_23_IM("cl_full_23_IM");
FVec cl_full_24_RE("cl_full_24_RE");
FVec cl_full_24_IM("cl_full_24_IM");
FVec cl_full_25_RE("cl_full_25_RE");
FVec cl_full_25_IM("cl_full_25_IM");
FVec cl_full_30_RE("cl_full_30_RE");
FVec cl_full_30_IM("cl_full_30_IM");
FVec cl_full_31_RE("cl_full_31_RE");
FVec cl_full_31_IM("cl_full_31_IM");
FVec cl_full_32_RE("cl_full_32_RE");
FVec cl_full_32_IM("cl_full_32_IM");
FVec cl_full_33_RE("cl_full_33_RE");
FVec cl_full_33_IM("cl_full_33_IM");
FVec cl_full_34_RE("cl_full_34_RE");
FVec cl_full_34_IM("cl_full_34_IM");
FVec cl_full_35_RE("cl_full_35_RE");
FVec cl_full_35_IM("cl_full_35_IM");
FVec cl_full_40_RE("cl_full_40_RE");
FVec cl_full_40_IM("cl_full_40_IM");
FVec cl_full_41_RE("cl_full_41_RE");
FVec cl_full_41_IM("cl_full_41_IM");
FVec cl_full_42_RE("cl_full_42_RE");
FVec cl_full_42_IM("cl_full_42_IM");
FVec cl_full_43_RE("cl_full_43_RE");
FVec cl_full_43_IM("cl_full_43_IM");
FVec cl_full_44_RE("cl_full_44_RE");
FVec cl_full_44_IM("cl_full_44_IM");
FVec cl_full_45_RE("cl_full_45_RE");
FVec cl_full_45_IM("cl_full_45_IM");
FVec cl_full_50_RE("cl_full_50_RE");
FVec cl_full_50_IM("cl_full_50_IM");
FVec cl_full_51_RE("cl_full_51_RE");
FVec cl_full_51_IM("cl_full_51_IM");
FVec cl_full_52_RE("cl_full_52_RE");
FVec cl_full_52_IM("cl_full_52_IM");
FVec cl_full_53_RE("cl_full_53_RE");
FVec cl_full_53_IM("cl_full_53_IM");
FVec cl_full_54_RE("cl_full_54_RE");
FVec cl_full_54_IM("cl_full_54_IM");
FVec cl_full_55_RE("cl_full_55_RE");
FVec cl_full_55_IM("cl_full_55_IM");

FVec clov_full[6][6][2] = {{
                               {cl_full_00_RE, cl_full_00_IM},
                               {cl_full_01_RE, cl_full_01_IM},
                               {cl_full_02_RE, cl_full_02_IM},
                               {cl_full_03_RE, cl_full_03_IM},
                               {cl_full_04_RE, cl_full_04_IM},
                               {cl_full_05_RE, cl_full_05_IM},
                           },
                           {
                               {cl_full_10_RE, cl_full_10_IM},
                               {cl_full_11_RE, cl_full_11_IM},
                               {cl_full_12_RE, cl_full_12_IM},
                               {cl_full_13_RE, cl_full_13_IM},
                               {cl_full_14_RE, cl_full_14_IM},
                               {cl_full_15_RE, cl_full_15_IM},
                           },
                           {
                               {cl_full_20_RE, cl_full_20_IM},
                               {cl_full_21_RE, cl_full_21_IM},
                               {cl_full_22_RE, cl_full_22_IM},
                               {cl_full_23_RE, cl_full_23_IM},
                               {cl_full_24_RE, cl_full_24_IM},
                               {cl_full_25_RE, cl_full_25_IM},
                           },
                           {
                               {cl_full_30_RE, cl_full_30_IM},
                               {cl_full_31_RE, cl_full_31_IM},
                               {cl_full_32_RE, cl_full_32_IM},
                               {cl_full_33_RE, cl_full_33_IM},
                               {cl_full_34_RE, cl_full_34_IM},
                               {cl_full_35_RE, cl_full_35_IM},
                           },
                           {
                               {cl_full_40_RE, cl_full_40_IM},
                               {cl_full_41_RE, cl_full_41_IM},
                               {cl_full_42_RE, cl_full_42_IM},
                               {cl_full_43_RE, cl_full_43_IM},
                               {cl_full_44_RE, cl_full_44_IM},
                               {cl_full_45_RE, cl_full_45_IM},
                           },
                           {
                               {cl_full_50_RE, cl_full_50_IM},
                               {cl_full_51_RE, cl_full_51_IM},
                               {cl_full_52_RE, cl_full_52_IM},
                               {cl_full_53_RE, cl_full_53_IM},
                               {cl_full_54_RE, cl_full_54_IM},
                               {cl_full_55_RE, cl_full_55_IM},
                           }};

FVec zero("zero");
FVec alpha_vec("alpha_vec");
FVec beta_vec("beta_vec");

FVec mu_vec("mu_vec");
FVec mu_inv_vec("mu_inv_vec");

FVec prec_mass_rho_vec("prec_mass_rho_vec");

FVec psi_S0_RE("psi_S0_RE");
FVec psi_S0_IM("psi_S0_IM");
FVec psi_S1_RE("psi_S1_RE");
FVec psi_S1_IM("psi_S1_IM");

FVec psi[2][2] = {{psi_S0_RE, psi_S0_IM}, {psi_S1_RE, psi_S1_IM}};

FVec tmp_1_re("tmp_1_re");
FVec tmp_1_im("tmp_1_im");

FVec tmp_2_re("tmp_2_re");
FVec tmp_2_im("tmp_2_im");

FVec tmp_3_re("tmp_3_re");
FVec tmp_3_im("tmp_3_im");

FVec tmp_4_re("tmp_4_re");
FVec tmp_4_im("tmp_4_im");

FVec tmp[4] = {tmp_1_re, tmp_2_re, tmp_3_re, tmp_4_re};

// for storing twisted boundary conditions phases
FVec tbc_phase_re("tbc_phase_re");
FVec tbc_phase_im("tbc_phase_im");
FVec tbc_phase[2] = {tbc_phase_re, tbc_phase_im};

FVec u_00_re("u_00_re");
FVec u_00_im("u_00_im");
FVec u_10_re("u_10_re");
FVec u_10_im("u_10_im");
FVec u_20_re("u_20_re");
FVec u_20_im("u_20_im");

FVec u_01_re("u_01_re");
FVec u_01_im("u_01_im");
FVec u_11_re("u_11_re");
FVec u_11_im("u_11_im");
FVec u_21_re("u_21_re");
FVec u_21_im("u_21_im");

FVec u_02_re("u_02_re");
FVec u_02_im("u_02_im");
FVec u_12_re("u_12_re");
FVec u_12_im("u_12_im");
FVec u_22_re("u_22_re");
FVec u_22_im("u_22_im");

FVec u_gauge[3][3][2] = {
    {{u_00_re, u_00_im}, {u_01_re, u_01_im}, {u_02_re, u_02_im}},
    {{u_10_re, u_10_im}, {u_11_re, u_11_im}, {u_12_re, u_12_im}},
    {{u_20_re, u_20_im}, {u_21_re, u_21_im}, {u_22_re, u_22_im}}};

void declare_b_Spins(InstVector &ivector)
{
  for (int s = 0; s < 2; s++) {
    for (int c = 0; c < 3; c++) {
      declareFVecFromFVec(ivector, b_spinor[s][c][RE]);
      declareFVecFromFVec(ivector, b_spinor[s][c][IM]);
    }
  }
}

void declare_ub_Spins(InstVector &ivector)
{
  for (int s = 0; s < 2; s++) {
    for (int c = 0; c < 3; c++) {
      declareFVecFromFVec(ivector, ub_spinor[s][c][RE]);
      declareFVecFromFVec(ivector, ub_spinor[s][c][IM]);
    }
  }
}

void declare_u_gaus(InstVector &ivector)
{
  for (int c1 = 0; c1 < 3; c1++) {
    for (int c2 = 0; c2 < 3; c2++) {
      declareFVecFromFVec(ivector, u_gauge[c1][c2][RE]);
      declareFVecFromFVec(ivector, u_gauge[c1][c2][IM]);
    }
  }
}

void declare_outs(InstVector &ivector, bool two_flav = false)
{
  for (int s = 0; s < 4; s++) {
    for (int c = 0; c < 3; c++) {
      declareFVecFromFVec(ivector, out_spinor[s][c][RE]);
      declareFVecFromFVec(ivector, out_spinor[s][c][IM]);
    }
  }
  if(two_flav){
    for (int s = 0; s < 4; s++) {
      for (int c = 0; c < 3; c++) {
        declareFVecFromFVec(ivector, out2_spinor[s][c][RE]);
        declareFVecFromFVec(ivector, out2_spinor[s][c][IM]);
      }
    }
  }
}

void declare_douts(InstVector &ivector)
{
  for (int s = 0; s < 4; s++) {
    for (int c = 0; c < 3; c++) {
      declareFVecFromFVec(ivector, dout_spinor[s][c][RE]);
      declareFVecFromFVec(ivector, dout_spinor[s][c][IM]);
    }
  }
}

void declare_chi(InstVector &ivector, bool two_flav = false)
{
  for (int s = 0; s < 4; s++) {
    for (int c = 0; c < 3; c++) {
      declareFVecFromFVec(ivector, chi_spinor[s][c][RE]);
      declareFVecFromFVec(ivector, chi_spinor[s][c][IM]);
    }
  }
  if(two_flav){
    for (int s = 0; s < 4; s++) {
      for (int c = 0; c < 3; c++) {
        declareFVecFromFVec(ivector, chi2_spinor[s][c][RE]);
        declareFVecFromFVec(ivector, chi2_spinor[s][c][IM]);
      }
    }
  }
}

void declare_clover(InstVector &ivector)
{
  for (int s = 0; s < 6; s++) {
    declareFVecFromFVec(ivector, clov_diag[s]);
  }

  for (int s = 0; s < 15; s++) {
    declareFVecFromFVec(ivector, clov_offdiag[s][RE]);
    declareFVecFromFVec(ivector, clov_offdiag[s][IM]);
  }
}

void declare_full_clover(InstVector &ivector)
{
  for (int cs1 = 0; cs1 < 6; cs1++) {
    for (int cs2 = 0; cs2 < 6; cs2++) {
      declareFVecFromFVec(ivector, clov_full[cs1][cs2][RE]);
      declareFVecFromFVec(ivector, clov_full[cs1][cs2][IM]);
    }
  }
}

void declare_misc(InstVector &ivector)
{
  declareFVecFromFVec(ivector, psi_S0_RE);
  declareFVecFromFVec(ivector, psi_S0_IM);
  declareFVecFromFVec(ivector, psi_S1_RE);
  declareFVecFromFVec(ivector, psi_S1_IM);

  declareFVecFromFVec(ivector, tmp_1_re);
  declareFVecFromFVec(ivector, tmp_1_im);
  declareFVecFromFVec(ivector, tmp_2_re);
  declareFVecFromFVec(ivector, tmp_2_im);
  declareFVecFromFVec(ivector, tmp_3_re);
  declareFVecFromFVec(ivector, tmp_3_im);
  declareFVecFromFVec(ivector, tmp_4_re);
  declareFVecFromFVec(ivector, tmp_4_im);

  declareFVecFromFVec(ivector, zero);
  setZero(ivector, zero);
}

void movCVec(InstVector &ivector, FVec *r, FVec *s1, string &mask)
{
  movFVec(ivector, r[RE], s1[RE], mask);
  movFVec(ivector, r[IM], s1[IM], mask);
}

void addCVec(
    InstVector &ivector, FVec *r, FVec *s1, FVec *s2, string &mask)
{
  addFVec(ivector, r[RE], s1[RE], s2[RE], mask);
  addFVec(ivector, r[IM], s1[IM], s2[IM], mask);
}

void subCVec(
    InstVector &ivector, FVec *r, FVec *s1, FVec *s2, string &mask)
{
  subFVec(ivector, r[RE], s1[RE], s2[RE], mask);
  subFVec(ivector, r[IM], s1[IM], s2[IM], mask);
}

void addiCVec(
    InstVector &ivector, FVec *r, FVec *s1, FVec *s2, string &mask)
{
  subFVec(ivector, r[RE], s1[RE], s2[IM], mask);
  addFVec(ivector, r[IM], s1[IM], s2[RE], mask);
}

void subiCVec(
    InstVector &ivector, FVec *r, FVec *s1, FVec *s2, string &mask)
{
  addFVec(ivector, r[RE], s1[RE], s2[IM], mask);
  subFVec(ivector, r[IM], s1[IM], s2[RE], mask);
}

// r[RE] = s1[RE]-beta_vec*s2[RE] = fnmadd(beta_vec,s2[RE],s1[RE])
// r[IM] = s1[IM]-beta_vec*s2[IM] = fnamdd(beta_vec,s2[IM],s1[IM])
void addCVec_mbeta(
    InstVector &ivector, FVec *r, FVec *s1, FVec *s2, string &mask)
{
  fnmaddFVec(ivector, r[RE], beta_vec, s2[RE], s1[RE], mask);
  fnmaddFVec(ivector, r[IM], beta_vec, s2[IM], s1[IM], mask);
}
// r[RE] = s1[RE] + beta_vec*s2[RE] = fmadd(beta_vec, s2[RE],
// s1[RE]);
// r[IM] = s1[IM] + beta_vec*s2[IM] = fmadd(beta_vec, s2[IM],
// s1[IM]);
void subCVec_mbeta(
    InstVector &ivector, FVec *r, FVec *s1, FVec *s2, string &mask)
{
  fmaddFVec(ivector, r[RE], beta_vec, s2[RE], s1[RE], mask);
  fmaddFVec(ivector, r[IM], beta_vec, s2[IM], s1[IM], mask);
}

// r[RE] = s1[RE] + beta_vec * s2[IM] = fmadd(beta_vec,s2[IM],
// s1[RE])
// r[IM] = s1[IM] - beta_vec * s2[RE] = fnmadd(beta_vec, s2[RE],
// s1[IM])

void addiCVec_mbeta(
    InstVector &ivector, FVec *r, FVec *s1, FVec *s2, string &mask)
{
  fmaddFVec(ivector, r[RE], beta_vec, s2[IM], s1[RE], mask);
  fnmaddFVec(ivector, r[IM], beta_vec, s2[RE], s1[IM], mask);
}

// r[RE] = s1[RE] - beta_vec*s2[IM] = fnmadd( beta_vec, s2[IM],
// s1[RE]);
// r[IM] = s1[IM] + beta_vec*s2[RE] = fmadd ( beta_vec, s2[RE],
// s1[IM]);
void subiCVec_mbeta(
    InstVector &ivector, FVec *r, FVec *s1, FVec *s2, string &mask)
{
  fnmaddFVec(ivector, r[RE], beta_vec, s2[IM], s1[RE], mask);
  fmaddFVec(ivector, r[IM], beta_vec, s2[RE], s1[IM], mask);
}

// r[RE] = s1[RE]+beta_vec*s2[RE] = fmadd(beta_vec,s2[RE],s1[RE])
// r[IM] = s1[IM]+beta_vec*s2[IM] = fmadd(beta_vec,s2[IM],s1[IM])
void addCVec_pbeta(
    InstVector &ivector, FVec *r, FVec *s1, FVec *s2, string &mask)
{
  fmaddFVec(ivector, r[RE], beta_vec, s2[RE], s1[RE], mask);
  fmaddFVec(ivector, r[IM], beta_vec, s2[IM], s1[IM], mask);
}
// r[RE] = s1[RE] - beta_vec*s2[RE] = fnmadd(beta_vec, s2[RE],
// s1[RE]);
// r[IM] = s1[IM] - beta_vec*s2[IM] = fnmadd(beta_vec, s2[IM],
// s1[IM]);
void subCVec_pbeta(
    InstVector &ivector, FVec *r, FVec *s1, FVec *s2, string &mask)
{
  fnmaddFVec(ivector, r[RE], beta_vec, s2[RE], s1[RE], mask);
  fnmaddFVec(ivector, r[IM], beta_vec, s2[IM], s1[IM], mask);
}

// r[RE] = s1[RE] - beta_vec * s2[IM] = fnmadd(beta_vec,s2[IM],
// s1[RE])
// r[IM] = s1[IM] + beta_vec * s2[RE] = fmadd(beta_vec, s2[RE],
// s1[IM])
void addiCVec_pbeta(
    InstVector &ivector, FVec *r, FVec *s1, FVec *s2, string &mask)
{
  fnmaddFVec(ivector, r[RE], beta_vec, s2[IM], s1[RE], mask);
  fmaddFVec(ivector, r[IM], beta_vec, s2[RE], s1[IM], mask);
}

// r[RE] = s1[RE] + beta_vec*s2[IM] = fmadd( beta_vec, s2[IM],
// s1[RE]);
// r[IM] = s1[IM] - beta_vec*s2[RE] = fnmadd ( beta_vec, s2[RE],
// s1[IM]);
void subiCVec_pbeta(
    InstVector &ivector, FVec *r, FVec *s1, FVec *s2, string &mask)
{
  fmaddFVec(ivector, r[RE], beta_vec, s2[IM], s1[RE], mask);
  fnmaddFVec(ivector, r[IM], beta_vec, s2[RE], s1[IM], mask);
}

// r = s1*s2
void mulCVec(
    InstVector &ivector, FVec *r, FVec *s1, FVec *s2, string &mask)
{
  mulFVec(ivector, r[RE], s1[RE], s2[RE], mask);
  fnmaddFVec(ivector, r[RE], s1[IM], s2[IM], r[RE], mask);
  mulFVec(ivector, r[IM], s1[RE], s2[IM], mask);
  fmaddFVec(ivector, r[IM], s1[IM], s2[RE], r[IM], mask);
}

// r = s1*s2+s3
void fmaddCVec(InstVector &ivector,
               FVec *r,
               FVec *s1,
               FVec *s2,
               FVec *s3,
               string &mask)
{
  fmaddFVec(ivector, r[RE], s1[RE], s2[RE], s3[RE], mask);
  fnmaddFVec(ivector, r[RE], s1[IM], s2[IM], r[RE], mask);
  fmaddFVec(ivector, r[IM], s1[RE], s2[IM], s3[IM], mask);
  fmaddFVec(ivector, r[IM], s1[IM], s2[RE], r[IM], mask);
}

// r = s3-s1*s2
// r[RE] = (s3[RE]-s1[RE]*s2[RE])+(s1[IM]*s2[IM])
// r[IM] = (s3[IM]-s1[RE]*s2[IM])-(s1[IM]*s2[RE])
void fnmaddCVec(InstVector &ivector,
                FVec *r,
                FVec *s1,
                FVec *s2,
                FVec *s3,
                string &mask)
{
  fnmaddFVec(ivector, r[RE], s1[RE], s2[RE], s3[RE], mask);
  fmaddFVec(ivector, r[RE], s1[IM], s2[IM], r[RE], mask);
  fnmaddFVec(ivector, r[IM], s1[RE], s2[IM], s3[IM], mask);
  fnmaddFVec(ivector, r[IM], s1[IM], s2[RE], r[IM], mask);
}

// r = (s1*s2-s3*s4)'
// r[RE] =
// (s1[RE]*s2[RE])-(s1[IM]*s2[IM])-(s3[RE]*s4[RE])+(s3[IM]*s4[IM])
// r[IM] =
// (s3[RE]*s4[IM])+(s3[IM]*s4[RE])-(s1[RE]*s2[IM])-(s1[IM]*s2[RE])
void Conj_CrossProd(InstVector &ivector,
                    FVec *r,
                    FVec *s1,
                    FVec *s2,
                    FVec *s3,
                    FVec *s4,
                    string &mask)
{
  mulFVec(ivector, r[RE], s1[RE], s2[RE], mask);
  fnmaddFVec(ivector, r[RE], s1[IM], s2[IM], r[RE], mask);
  fnmaddFVec(ivector, r[RE], s3[RE], s4[RE], r[RE], mask);
  fmaddFVec(ivector, r[RE], s3[IM], s4[IM], r[RE], mask);

  mulFVec(ivector, r[IM], s3[RE], s4[IM], mask);
  fmaddFVec(ivector, r[IM], s3[IM], s4[RE], r[IM], mask);
  fnmaddFVec(ivector, r[IM], s1[RE], s2[IM], r[IM], mask);
  fnmaddFVec(ivector, r[IM], s1[IM], s2[RE], r[IM], mask);
}

// r = s1'*s2
// r[RE] = (s1[RE]*s2[RE])+(s1[IM]*s2[IM])
// r[IM] = (s1[RE]*s2[IM])-(s1[IM]*s2[RE])
void mulConjCVec(
    InstVector &ivector, FVec *r, FVec *s1, FVec *s2, string &mask)
{
  mulFVec(ivector, r[RE], s1[RE], s2[RE], mask);
  fmaddFVec(ivector, r[RE], s1[IM], s2[IM], r[RE], mask);
  mulFVec(ivector, r[IM], s1[RE], s2[IM], mask);
  fnmaddFVec(ivector, r[IM], s1[IM], s2[RE], r[IM], mask);
}

// r = s1'*s2+s3
// r[RE] = (s3[RE]+s1[RE]*s2[RE])+(s1[IM]*s2[IM])
// r[IM] = (s3[IM]+s1[RE]*s2[IM])-(s1[IM]*s2[RE])
void fmaddConjCVec(InstVector &ivector,
                   FVec *r,
                   FVec *s1,
                   FVec *s2,
                   FVec *s3,
                   string &mask)
{
  fmaddFVec(ivector, r[RE], s1[RE], s2[RE], s3[RE], mask);
  fmaddFVec(ivector, r[RE], s1[IM], s2[IM], r[RE], mask);
  fmaddFVec(ivector, r[IM], s1[RE], s2[IM], s3[IM], mask);
  fnmaddFVec(ivector, r[IM], s1[IM], s2[RE], r[IM], mask);
}

proj_ops proj_ops_plus[] = {
    {"plus_X_back", {{0, 3}, {1, 2}}, {addiCVec, addiCVec}},
    {"minus_X_forw", {{0, 3}, {1, 2}}, {subiCVec, subiCVec}},
    {"plus_Y", {{0, 3}, {1, 2}}, {subCVec, addCVec}},
    {"minus_Y", {{0, 3}, {1, 2}}, {addCVec, subCVec}},
    {"plus_Z", {{0, 2}, {1, 3}}, {addiCVec, subiCVec}},
    {"minus_Z", {{0, 2}, {1, 3}}, {subiCVec, addiCVec}},
    {"plus_T", {{0, 2}, {1, 3}}, {addCVec, addCVec}},
    {"minus_T", {{0, 2}, {1, 3}}, {subCVec, subCVec}}};

proj_ops proj_ops_minus[] = {
    {"minus_X_back", {{0, 3}, {1, 2}}, {subiCVec, subiCVec}},
    {"plus_X_forw", {{0, 3}, {1, 2}}, {addiCVec, addiCVec}},
    {"minus_Y", {{0, 3}, {1, 2}}, {addCVec, subCVec}},
    {"plus_Y", {{0, 3}, {1, 2}}, {subCVec, addCVec}},
    {"minus_Z", {{0, 2}, {1, 3}}, {subiCVec, addiCVec}},
    {"plus_Z", {{0, 2}, {1, 3}}, {addiCVec, subiCVec}},
    {"minus_T", {{0, 2}, {1, 3}}, {subCVec, subCVec}},
    {"plus_T", {{0, 2}, {1, 3}}, {addCVec, addCVec}}};

/*
recons_ops rec_plus_ops[] = {
        {"recons_plus_X", 1,0, addCVec, subiCVec, subiCVec},
        {"recons_plus_Y", 1,0, addCVec, addCVec,  subCVec},
        {"recons_plus_Z", 0,1, addCVec, subiCVec, addiCVec},
        {"recons_plus_T", 0,1, addCVec, addCVec,  addCVec},
};

recons_ops rec_minus_ops[] = {
        {"recons_minus_X", 1,0, addCVec, addiCVec, addiCVec},
        {"recons_minus_Y", 1,0, addCVec, subCVec,  addCVec},
        {"recons_minus_Z", 0,1, addCVec, addiCVec, subiCVec},
        {"recons_minus_T", 0,1, addCVec, subCVec,  subCVec},
};
*/

recons_ops rec_plus_pbeta_ops[] = {
    {"recons_plus_X_pbeta",
     1,
     0,
     addCVec_pbeta,
     subiCVec_pbeta,
     subiCVec_pbeta},
    {"recons_plus_Y_pbeta",
     1,
     0,
     addCVec_pbeta,
     addCVec_pbeta,
     subCVec_pbeta},
    {"recons_plus_Z_pbeta",
     0,
     1,
     addCVec_pbeta,
     subiCVec_pbeta,
     addiCVec_pbeta},
    {"recons_plus_T_pbeta",
     0,
     1,
     addCVec_pbeta,
     addCVec_pbeta,
     addCVec_pbeta},
};

recons_ops rec_minus_pbeta_ops[] = {
    {"recons_minus_X",
     1,
     0,
     addCVec_pbeta,
     addiCVec_pbeta,
     addiCVec_pbeta},
    {"recons_minus_Y",
     1,
     0,
     addCVec_pbeta,
     subCVec_pbeta,
     addCVec_pbeta},
    {"recons_minus_Z",
     0,
     1,
     addCVec_pbeta,
     addiCVec_pbeta,
     subiCVec_pbeta},
    {"recons_minus_T",
     0,
     1,
     addCVec_pbeta,
     subCVec_pbeta,
     subCVec_pbeta},
};

recons_ops rec_plus_mbeta_ops[] = {
    {"recons_plus_X_mbeta",
     1,
     0,
     addCVec_mbeta,
     subiCVec_mbeta,
     subiCVec_mbeta},
    {"recons_plus_Y_mbeta",
     1,
     0,
     addCVec_mbeta,
     addCVec_mbeta,
     subCVec_mbeta},
    {"recons_plus_Z_mbeta",
     0,
     1,
     addCVec_mbeta,
     subiCVec_mbeta,
     addiCVec_mbeta},
    {"recons_plus_T_mbeta",
     0,
     1,
     addCVec_mbeta,
     addCVec_mbeta,
     addCVec_mbeta},
};

recons_ops rec_minus_mbeta_ops[] = {
    {"recons_minus_X",
     1,
     0,
     addCVec_mbeta,
     addiCVec_mbeta,
     addiCVec_mbeta},
    {"recons_minus_Y",
     1,
     0,
     addCVec_mbeta,
     subCVec_mbeta,
     addCVec_mbeta},
    {"recons_minus_Z",
     0,
     1,
     addCVec_mbeta,
     addiCVec_mbeta,
     subiCVec_mbeta},
    {"recons_minus_T",
     0,
     1,
     addCVec_mbeta,
     subCVec_mbeta,
     subCVec_mbeta},
};

// Uses ub_spinor as implicit input
// Uses psi[][] as temp and b_spinor[][] as implicit output
void project(InstVector &ivector,
             string base,
             string offset,
             proj_ops &ops,
             bool isFace,
             string mask,
             int dir)
{
  string tmask("");
  PrefetchL1FullSpinorDirIn(ivector, base, offset, dir);

  for (int s = 0; s < 2; s++) {
    for (int c = 0; c < 3; c++) {
      LoadSpinorElement(ivector,
                        psi[0][RE],
                        base,
                        offset,
                        ops.s[s][0],
                        c,
                        RE,
                        isFace,
                        mask,
                        dir);
      LoadSpinorElement(ivector,
                        psi[0][IM],
                        base,
                        offset,
                        ops.s[s][0],
                        c,
                        IM,
                        isFace,
                        mask,
                        dir);
      LoadSpinorElement(ivector,
                        psi[1][RE],
                        base,
                        offset,
                        ops.s[s][1],
                        c,
                        RE,
                        isFace,
                        mask,
                        dir);
      LoadSpinorElement(ivector,
                        psi[1][IM],
                        base,
                        offset,
                        ops.s[s][1],
                        c,
                        IM,
                        isFace,
                        mask,
                        dir);

      ops.CVecFunc[s](ivector,
                      b_spinor[s][c],
                      psi[0],
                      psi[1],
                      /*mask*/ tmask); // Not using mask here
    }
  }
}

// Serial Spin version
void project(InstVector &ivector,
             string base,
             string offset,
             proj_ops &ops,
             bool isFace,
             string mask,
             int dir,
             int s)
{
  string tmask("");

  if (s == 0) {
    PrefetchL1FullSpinorDirIn(ivector, base, offset, dir);
  }

  for (int c = 0; c < 3; c++) {
    LoadSpinorElement(ivector,
                      psi[0][RE],
                      base,
                      offset,
                      ops.s[s][0],
                      c,
                      RE,
                      isFace,
                      mask,
                      dir);
    LoadSpinorElement(ivector,
                      psi[0][IM],
                      base,
                      offset,
                      ops.s[s][0],
                      c,
                      IM,
                      isFace,
                      mask,
                      dir);
    LoadSpinorElement(ivector,
                      psi[1][RE],
                      base,
                      offset,
                      ops.s[s][1],
                      c,
                      RE,
                      isFace,
                      mask,
                      dir);
    LoadSpinorElement(ivector,
                      psi[1][IM],
                      base,
                      offset,
                      ops.s[s][1],
                      c,
                      IM,
                      isFace,
                      mask,
                      dir);

    ops.CVecFunc[s](ivector,
                    b_spinor[s][c],
                    psi[0],
                    psi[1],
                    /*mask*/ tmask); // Not using mask here
  }
}

void recons_add(InstVector &ivector,
                recons_ops &ops,
                FVec outspinor[4][3][2],
                string &mask)
{
  for (int s = 0; s < 2; s++) {
    for (int c = 0; c < 3; c++) {
      ops.CVecFuncTop2(ivector,
                       outspinor[s][c],
                       outspinor[s][c],
                       ub_spinor[s][c],
                       mask);
    }

    if (ops.s2 == s) {
      for (int c = 0; c < 3; c++) {
        ops.CVecFunc1(ivector,
                      outspinor[2][c],
                      outspinor[2][c],
                      ub_spinor[s][c],
                      mask);
      }
    } else {
      for (int c = 0; c < 3; c++) {
        ops.CVecFunc2(ivector,
                      outspinor[3][c],
                      outspinor[3][c],
                      ub_spinor[s][c],
                      mask);
      }
    }
  }
}

// Serial Spin version
void recons_add(InstVector &ivector,
                recons_ops &ops,
                FVec outspinor[4][3][2],
                string &mask,
                int s)
{
  for (int c = 0; c < 3; c++) {
    ops.CVecFuncTop2(ivector,
                     outspinor[s][c],
                     outspinor[s][c],
                     ub_spinor[s][c],
                     mask);
  }

  if (ops.s2 == s) {
    for (int c = 0; c < 3; c++) {
      ops.CVecFunc1(ivector,
                    outspinor[2][c],
                    outspinor[2][c],
                    ub_spinor[s][c],
                    mask);
    }
  } else {
    for (int c = 0; c < 3; c++) {
      ops.CVecFunc2(ivector,
                    outspinor[3][c],
                    outspinor[3][c],
                    ub_spinor[s][c],
                    mask);
    }
  }
}

void zeroResult(InstVector &ivector, FVec *outspinor)
{
  for (int i = 0; i < 24; i++) {
    setZero(ivector, outspinor[i]);
  }
}

void two_flav_clover_term(InstVector &ivector,
                          FVec in_spinor[4][3][2],
                          FVec in2_spinor[4][3][2],
                          FVec out_spinor[4][3][2],
                          FVec out2_spinor[4][3][2],
                          string _mask,
                          bool acc){

  FVec clout_tmp[2] = {tmp_1_re, tmp_1_im};
  for (int block = 0; block < 2; block++) {
    // BaKo, August 2017: we leave out all prefetches because they don's seem to
    // do anything here
    //PrefetchL1FullCloverBlockIn(ivector, clBase, clOffs, block);
    LoadFullCloverBlock(
        ivector, clov_diag, clov_offdiag, clBase, clOffs, block);

    for (int c1 = 0; c1 < 6; c1++) {
      int spin = 2 * block + c1 / 3;
      int col = c1 % 3;
      string mask = _mask;

      // tau1 in flavour imlpemented on the in spinors
      FVec *clout = out_spinor[spin][col];
      FVec *clout2 = out2_spinor[spin][col];
      FVec *clin = in2_spinor[spin][col];
      FVec *clin2 = in_spinor[spin][col];
#ifdef NO_HW_MASKING

      if (_mask != "") {
        acc = false;
        clout = clout_tmp;
        mask = "";
      }

#endif

      if (acc) {
        fmaddFVec(ivector, clout[RE], clov_diag[c1], clin[RE], clout[RE], mask);
        fmaddFVec(ivector, clout2[RE], clov_diag[c1], clin2[RE], clout2[RE], mask);
        fmaddFVec(ivector, clout[IM], clov_diag[c1], clin[IM], clout[IM], mask);
        fmaddFVec(ivector, clout2[IM], clov_diag[c1], clin2[IM], clout2[IM], mask);
      } else {
        mulFVec(ivector, clout[RE], clov_diag[c1], clin[RE], mask);
        mulFVec(ivector, clout2[RE], clov_diag[c1], clin2[RE], mask);
        mulFVec(ivector, clout[IM], clov_diag[c1], clin[IM], mask);
        mulFVec(ivector, clout2[IM], clov_diag[c1], clin2[IM], mask);
      }

      for (int c2 = 0; c2 < 6; c2++) {
        if (c1 == c2) {
          continue; // diagonal case
        }

        if (c1 < c2) {
          int od = c2 * (c2 - 1) / 2 + c1;
          // note the tau1 in flavour
          fmaddConjCVec(ivector, clout, clov_offdiag[od], in2_spinor[2 * block + c2 / 3][c2 % 3], clout, mask);
          fmaddConjCVec(ivector, clout2, clov_offdiag[od], in_spinor[2 * block + c2 / 3][c2 % 3], clout2, mask);
        } else {
          int od = c1 * (c1 - 1) / 2 + c2;
          // note tha tau1 in flavour
          fmaddCVec(ivector, clout, clov_offdiag[od], in2_spinor[2 * block + c2 / 3][c2 % 3], clout, mask);
          fmaddCVec(ivector, clout2, clov_offdiag[od], in_spinor[2 * block + c2 / 3][c2 % 3], clout2, mask);
        }
      } // c2

#ifdef NO_HW_MASKING

      if (_mask != "") {
        if(acc){
          addCVec(ivector, out_spinor[spin][col], clout, out_spinor[block][c1], _mask);
          addCVec(ivector, out2_spinor[spin][col], clout2, out2_spinor[block][c1], _mask);
        }else{
          movCVec(ivector, out_spinor[spin][col], clout, _mask);
          movCVec(ivector, out2_spinor[spin][col], clout2, _mask);
        }
      }

#endif
    } // c1
  } // block

}

void two_flav_full_clover_term(InstVector &ivector,
                               FVec in_spinor[4][3][2],
                               FVec out_spinor[4][3][2],
                               string cloverBase,
                               string _mask,
                               bool acc){
  for (int block = 0; block < 2; block++) {

    // BaKo, August 2017: we leave out all prefetches because they don's seem to
    // do anything here
    //PrefetchL1FullCloverFullBlockIn(ivector, cloverBase, clOffs, block);
    LoadFullCloverFullBlock(
        ivector, clov_full, cloverBase, clOffs, block);

    for (int sc1 = 0; sc1 < 6; sc1++) { // half-spin-colour row

      int spin_out = 2 * block + sc1 / 3;
      int col_out = sc1 % 3;
      FVec *clout = out_spinor[spin_out][col_out];

      for (int sc2 = 0; sc2 < 6; sc2++) { // half-spin-colour column

        int spin_in = 2 * block + sc2 / 3;
        int col_in = sc2 % 3;
        FVec *clin = in_spinor[spin_in][col_in];

        if (sc2 == 0 && !acc) {
          mulCVec(ivector, clout, clov_full[sc1][sc2], clin, _mask);
        } else {
          fmaddCVec(ivector,
                    clout,
                    clov_full[sc1][sc2],
                    clin,
                    clout,
                    _mask);
        }
      } // half-spin-colour column
    } // half-spin-colour row
  } // block
}

void two_flav_tm_inverse_clover_term(InstVector &ivector,
                                     string _mask)
{
  declare_outs(ivector, true);
  declare_chi(ivector, true);
  declare_clover(ivector);
  declare_full_clover(ivector);

  LoadFullSpinor(ivector, chi_spinor, chiBase, chiOffs, _mask);
  LoadFullSpinor(ivector, chi2_spinor, chi2Base, chiOffs, _mask);
  // first we apply the inverse full clover term on the diagonal
  // we do so flavour by flavour because otherwise we would jump around in memory a lot
  two_flav_full_clover_term(ivector, chi_spinor, out_spinor, fclBase, _mask, false);
  two_flav_full_clover_term(ivector, chi2_spinor, out2_spinor, fcl2Base, _mask, false);
  // and then we add the flavour-off-diagonal contribution which comes with
  // just a Wilson-type clover term with a real spin-colour diagonal
  // it makes sense to do this in one go because the same term is applied to both flavours
  two_flav_clover_term(ivector, chi_spinor, chi2_spinor,
                       out_spinor, out2_spinor, _mask, true);

}


void clover_term(InstVector &ivector,
                 FVec in_spinor[4][3][2],
                 bool face,
                 string _mask)
{
  FVec clout_tmp[2] = {tmp_1_re, tmp_1_im};

  for (int block = 0; block < 2; block++) {
    PrefetchL1FullCloverBlockIn(ivector, clBase, clOffs, block);
    LoadFullCloverBlock(
        ivector, clov_diag, clov_offdiag, clBase, clOffs, block);

    for (int c1 = 0; c1 < 6; c1++) {
      int spin = 2 * block + c1 / 3;
      int col = c1 % 3;
      bool acc = face;
      string mask = _mask;
      FVec *clout = out_spinor[spin][col];
      FVec *clin = in_spinor[spin][col];
#ifdef NO_HW_MASKING

      if (_mask != "") {
        acc = false;
        clout = clout_tmp;
        mask = "";
      }

#endif

      if (acc) {
        fmaddFVec(ivector,
                  clout[RE],
                  clov_diag[c1],
                  clin[RE],
                  clout[RE],
                  mask);
        fmaddFVec(ivector,
                  clout[IM],
                  clov_diag[c1],
                  clin[IM],
                  clout[IM],
                  mask);
      } else {
        mulFVec(ivector, clout[RE], clov_diag[c1], clin[RE], mask);
        mulFVec(ivector, clout[IM], clov_diag[c1], clin[IM], mask);
      }

      for (int c2 = 0; c2 < 6; c2++) {
        if (c1 == c2) {
          continue; // diagonal case
        }

        if (c1 < c2) {
          int od = c2 * (c2 - 1) / 2 + c1;
          fmaddConjCVec(ivector,
                        clout,
                        clov_offdiag[od],
                        in_spinor[2 * block + c2 / 3][c2 % 3],
                        clout,
                        mask);
        } else {
          int od = c1 * (c1 - 1) / 2 + c2;
          fmaddCVec(ivector,
                    clout,
                    clov_offdiag[od],
                    in_spinor[2 * block + c2 / 3][c2 % 3],
                    clout,
                    mask);
        }
      }

#ifdef NO_HW_MASKING

      if (_mask != "") {
        if (face) {
          addCVec(ivector,
                  out_spinor[spin][col],
                  clout,
                  clout_spinor[block][c1],
                  _mask);
        } else {
          movCVec(ivector, out_spinor[spin][col], clout, _mask);
        }
      }

#endif
    }
  }
}

void full_clover_term(InstVector &ivector,
                      FVec in_spinor[4][3][2],
                      bool face,
                      string mask)
{
  // FIXME: Function is never called with non-zero mask! Do we need
  // this argument here?
  if (mask != "") {
    printf("full_clover_term:: ERROR: Non-empty mask is not "
           "implemented!\n");
  }

  for (int block = 0; block < 2; block++) {

    PrefetchL1FullCloverFullBlockIn(ivector, clBase, clOffs, block);
    LoadFullCloverFullBlock(
        ivector, clov_full, clBase, clOffs, block);

    for (int sc1 = 0; sc1 < 6; sc1++) { // half-spin-colour row

      int spin_out = 2 * block + sc1 / 3;
      int col_out = sc1 % 3;
      FVec *clout = out_spinor[spin_out][col_out];

      for (int sc2 = 0; sc2 < 6; sc2++) { // half-spin-colour column

        int spin_in = 2 * block + sc2 / 3;
        int col_in = sc2 % 3;
        FVec *clin = in_spinor[spin_in][col_in];

        if (sc2 == 0 && !face) {
          mulCVec(ivector, clout, clov_full[sc1][sc2], clin, mask);
        } else {
          fmaddCVec(ivector,
                    clout,
                    clov_full[sc1][sc2],
                    clin,
                    clout,
                    mask);
        }

      } // half-spin-colour column
    } // half-spin-colour row
  } // block
}

void inverse_twisted_term(InstVector &ivector,
                          FVec in_spinor[4][3][2],
                          bool face,
                          bool isPlus,
                          string _mask)
{
  /**

    This routine generates the result spinor of the matrix
    multiplication

      \chi = A^{-1} \psi

    for a tile of sites (which is of size of a SIMD vector). A
    denotes the
    sum of the Wilson mass and the twisted mass term. Its inverse is
    given by

               \alpha 1I - i \mu \gamma_5
      A^{-1} = --------------------------- 1I^{color} .
                    \alpha^2 + \mu^2

    Here, \mu is the twisted mass and \alpha = 4 + M_0, the Wilson
    mass term.
    Note that the TMDlash member function will pass the following
    parameters:

      mu     <---| \mu / \alpha
      mu_inv <---| \alpha / (\alpha^2 + \mu^2)

    This allows to write the matrix multiplication as a fused
    multiply add, and
    a rescaling with mu_inv:

     spin 0, 1:
     ----------
      Re \chi = mu_inv * ( Re \psi + mu * Im \psi )
      Im \chi = mu_inv * ( - mu * Re \psi + Im \psi )

     spin 2, 3:
     ----------
      Re \chi = mu_inv * ( Re \psi - mu * Im \psi )
      Im \chi = mu_inv * ( mu * Re \psi + Im \psi )

    For the hermitian-conjugate multiplication (isPlus==false) the
    case for
    spins 0, 1 is interchanged with the case for spins 2, 3.


    Input:
      \psi = in_spinor[spin][color][RE/IM]

    Output:
      \chi = out_spinor[spin][color][RE/IM]

   */

  bool acc = face;

  // Declare vector variables and set them to the scalar value
  // passed to the kernel
  declareFVecFromFVec(ivector, mu_vec);
  declareFVecFromFVec(ivector, mu_inv_vec);
  loadBroadcastScalar(ivector, mu_vec, mu_name, SpinorType);
  loadBroadcastScalar(ivector, mu_inv_vec, mu_inv_name, SpinorType);

  for (int col = 0; col < 3; col++) {

    // Upper half spinor
    for (int spin = 0; spin < 2; spin++) {
      FVec *in = in_spinor[spin][col];
      FVec *out = out_spinor[spin][col];

      // fused multiply add using mu
      if (isPlus) { // normal case
        fmaddFVec(ivector, tmp_1_re, mu_vec, in[IM], in[RE], _mask);
        fnmaddFVec(
            ivector, tmp_1_im, mu_vec, in[RE], in[IM], _mask);
      } else { // hermitian-conjugate case
        fnmaddFVec(
            ivector, tmp_1_re, mu_vec, in[IM], in[RE], _mask);
        fmaddFVec(ivector, tmp_1_im, mu_vec, in[RE], in[IM], _mask);
      }

      // rescaling with mu_inv
      if (acc) { // face processing
        fmaddFVec(
            ivector, out[RE], mu_inv_vec, tmp_1_re, out[RE], _mask);
        fmaddFVec(
            ivector, out[IM], mu_inv_vec, tmp_1_im, out[IM], _mask);
      } else { // body processing
        mulFVec(ivector, in[RE], mu_inv_vec, tmp_1_re, _mask);
        mulFVec(ivector, in[IM], mu_inv_vec, tmp_1_im, _mask);
      }
    }

    // Lower half spinor
    for (int spin = 2; spin < 4; spin++) {
      FVec *in = in_spinor[spin][col];
      FVec *out = out_spinor[spin][col];

      // fused multiply add using mu
      if (isPlus) { // normal case
        fnmaddFVec(
            ivector, tmp_1_re, mu_vec, in[IM], in[RE], _mask);
        fmaddFVec(ivector, tmp_1_im, mu_vec, in[RE], in[IM], _mask);
      } else { // hermitian-conjugate case
        fmaddFVec(ivector, tmp_1_re, mu_vec, in[IM], in[RE], _mask);
        fnmaddFVec(
            ivector, tmp_1_im, mu_vec, in[RE], in[IM], _mask);
      }

      // rescaling with mu_inv
      if (acc) { // face processing
        fmaddFVec(
            ivector, out[RE], mu_inv_vec, tmp_1_re, out[RE], _mask);
        fmaddFVec(
            ivector, out[IM], mu_inv_vec, tmp_1_im, out[IM], _mask);
      } else { // body processing
        mulFVec(ivector, in[RE], mu_inv_vec, tmp_1_re, _mask);
        mulFVec(ivector, in[IM], mu_inv_vec, tmp_1_im, _mask);
      }
    }

  } // color
}

// TODO: REFACTOR: This function does the same as the hermitian
// conjugate of inverse_twisted_term, where mu_inv is replaced by
// alpha.
void twisted_term(InstVector &ivector, bool isPlus)
{
  /**

    This routine generates the result spinor of the matrix
    multiplication

      out_spinor = A \chi

    for a tile of sites (which is of size of a SIMD vector). A
    denotes the
    sum of the Wilson mass and the twisted mass term:

      A = ( \alpha 1I + i \mu \gamma_5 ) 1I^{color} .

    Here, \mu is the twisted mass and \alpha = 4 + M_0, the Wilson
    mass term.
    Note that the TMDlash member function will pass the following
    parameters:

      alpha <---| \alpha
      beta  <---| \beta  [not needed here]
      mu    <---| \mu / \alpha

    This allows to write the matrix multiplication as a fused
    multiply add, and
    a rescaling with mu_inv:

     spin 0, 1:
     ----------
      Re out_spinor = alpha * ( Re \chi - mu * Im \chi )
      Im out_spinor = alpha * ( mu * Re \chi + Im \chi )

     spin 2, 3:
     ----------
      Re out_spinor = alpha * ( Re \chi + mu * Im \chi )
      Im out_spinor = alpha * ( - mu * Re \chi + Im \chi )

    For the hermitian-conjugate multiplication (isPlus==false) the
    case for
    spins 0, 1 is interchanged with the case for spins 2, 3.


    Input:
      \chi read from the \chi spinor field passed to the kernel
    routine

    Output:
      out_spinor[spin][color][RE/IM]

   */

  // Declare vector variable mu and set it to the scalar value
  // passed to the kernel.
  // The variable alpha_vec has already been declared in
  // dslash_achimbdpsi_body
  declareFVecFromFVec(ivector, mu_vec);
  loadBroadcastScalar(ivector, mu_vec, mu_name, SpinorType);

  // Load all relevant elements of the chi input spinor
  for (int col = 0; col < 3; col++) {
    for (int spin = 0; spin < 4; spin++) {
      LoadSpinorElement(ivector,
                        chi_spinor[spin][col][RE],
                        chiBase,
                        chiOffs,
                        spin,
                        col,
                        RE,
                        false,
                        "");
      LoadSpinorElement(ivector,
                        chi_spinor[spin][col][IM],
                        chiBase,
                        chiOffs,
                        spin,
                        col,
                        IM,
                        false,
                        "");
    }
  }

  // Carry out the matrix multiplication in two steps using
  // the buffers tmp_1_re and tmp_1_im, and storing the final
  // result in out_spinor
  for (int col = 0; col < 3; col++) {
    // Upper half spinor
    for (int spin = 0; spin < 2; spin++) {
      FVec *in = chi_spinor[spin][col];
      FVec *out = out_spinor[spin][col];

      // fused multiply add using mu
      if (isPlus) { // normal case
        fnmaddFVec(ivector, tmp_1_re, mu_vec, in[IM], in[RE]);
        fmaddFVec(ivector, tmp_1_im, mu_vec, in[RE], in[IM]);
      } else { // hermitian-conjugate case
        fmaddFVec(ivector, tmp_1_re, mu_vec, in[IM], in[RE]);
        fnmaddFVec(ivector, tmp_1_im, mu_vec, in[RE], in[IM]);
      }

      // rescaling with alpha
      mulFVec(ivector, out[RE], alpha_vec, tmp_1_re);
      mulFVec(ivector, out[IM], alpha_vec, tmp_1_im);
    }

    // Lower half spinor
    for (int spin = 2; spin < 4; spin++) {
      FVec *in = chi_spinor[spin][col];
      FVec *out = out_spinor[spin][col];

      // fused multiply add using mu
      if (isPlus) { // normal case
        fmaddFVec(ivector, tmp_1_re, mu_vec, in[IM], in[RE]);
        fnmaddFVec(ivector, tmp_1_im, mu_vec, in[RE], in[IM]);
      } else { // hermitian-conjugate case
        fnmaddFVec(ivector, tmp_1_re, mu_vec, in[IM], in[RE]);
        fmaddFVec(ivector, tmp_1_im, mu_vec, in[RE], in[IM]);
      }

      // rescaling with alpha
      mulFVec(ivector, out[RE], alpha_vec, tmp_1_re);
      mulFVec(ivector, out[IM], alpha_vec, tmp_1_im);
    }
  } // color
}

void applyTwistedBoundaryConditions(InstVector &ivector,
                                    bool const adjMul,
                                    bool const has_tbc)
{
  string mask;
  FVec tbc_tmp[2] = {tmp_1_re, tmp_1_im};
  if (has_tbc) {
    for (int s = 0; s < 2; s++) {
      for (int c1 = 0; c1 < 3; c1++) {
        if (!adjMul) {
          mulCVec(ivector,
                  tbc_tmp,
                  tbc_phase,
                  ub_spinor[s][c1],
                  mask);
          movCVec(ivector,
                  ub_spinor[s][c1],
                  tbc_tmp,
                  mask);
        } else {
          mulConjCVec(ivector,
                      tbc_tmp,
                      tbc_phase,
                      ub_spinor[s][c1],
                      mask);
          movCVec(ivector,
                  ub_spinor[s][c1],
                  tbc_tmp,
                  mask);
        }
      }
    }
  }
}

#if 0
void achiResult(InstVector& ivector, bool clover)
{
    PrefetchL1FullSpinorDirIn(ivector, chiBase, chiOffs, -1, 1 /*NTA*/);

    if(!clover) {
        for(int col=0; col < 3; col++) {
            for(int spin=0; spin < 4; spin++) {
                LoadSpinorElement(ivector, tmp_1_re, chiBase, chiOffs, spin, col, RE, false, "");
                LoadSpinorElement(ivector, tmp_1_im, chiBase, chiOffs, spin, col, IM, false, "");
                mulFVec(ivector, out_spinor[spin][col][RE], alpha_vec, tmp_1_re);
                mulFVec(ivector, out_spinor[spin][col][IM], alpha_vec, tmp_1_im);
            }
        }
    }
    else {
        for(int col=0; col < 3; col++) {
            for(int spin=0; spin < 4; spin++) {
                LoadSpinorElement(ivector, chi_spinor[spin][col][RE], chiBase, chiOffs, spin, col, RE, false, "");
                LoadSpinorElement(ivector, chi_spinor[spin][col][IM], chiBase, chiOffs, spin, col, IM, false, "");
            }
        }

        // Apply clover term, and store result in out spinor.
        // This is only on the AChi - bDPsi op (achimbdpsi = true)
        // This is only in body kernel (face = false)
        clover_term(ivector, chi_spinor, false);
    }
}
#endif

// TODO: Eliminate concrete implementations and use this function as
// a case switcher only
void achiResult(InstVector &ivector,
                bool const clover,
                TwistedMassVariant const twisted_mass,
                bool const isPlus)
{
  PrefetchL1FullSpinorDirIn(
      ivector, chiBase, chiOffs, -1, 1 /*NTA*/);

  if (clover) { // CLOVER [with or without twisted-mass]
    // Load all relevant elements of the chi input spinor
    for (int col = 0; col < 3; col++) {
      for (int spin = 0; spin < 4; spin++) {
        LoadSpinorElement(ivector,
                          chi_spinor[spin][col][RE],
                          chiBase,
                          chiOffs,
                          spin,
                          col,
                          RE,
                          false,
                          "");
        LoadSpinorElement(ivector,
                          chi_spinor[spin][col][IM],
                          chiBase,
                          chiOffs,
                          spin,
                          col,
                          IM,
                          false,
                          "");
      }
    }
    
    // BaKo, August 2017: for AChiMBDPsi we can use the standard clover
    // term also for twisted mass operators. We pass in the twisted mass
    // via the pre-conditioning mass "rho"
    clover_term(ivector, chi_spinor, false);
  } else { // NO CLOVER
    if (twisted_mass == TwistedMassVariant::none) {
      for (int col = 0; col < 3; col++) {
        for (int spin = 0; spin < 4; spin++) {
          LoadSpinorElement(ivector,
                            tmp_1_re,
                            chiBase,
                            chiOffs,
                            spin,
                            col,
                            RE,
                            false,
                            "");
          LoadSpinorElement(ivector,
                            tmp_1_im,
                            chiBase,
                            chiOffs,
                            spin,
                            col,
                            IM,
                            false,
                            "");
          mulFVec(ivector,
                  out_spinor[spin][col][RE],
                  alpha_vec,
                  tmp_1_re);
          mulFVec(ivector,
                  out_spinor[spin][col][IM],
                  alpha_vec,
                  tmp_1_im);
        }
      }
    } else if (twisted_mass == TwistedMassVariant::degenerate) {
      // Loads spinor elements from chiBase, multiplies with A
      // for pure twisted mass and stores result in out_spinor
      twisted_term(ivector, isPlus);
    } else if (twisted_mass == TwistedMassVariant::non_degenerate) {
      // BaKo, August 2017:
      // the full ND twisted mass term with mass splitting is implemented
      // as a linalg functor in QPhiX, for now this branch is never
      // entered

      // Loads spinor elements from chiBase, multiplies with A
      // for pure twisted mass and stores result in out_spinor
      twisted_term(ivector, isPlus);
    } else {
      unsupported_twisted_mass_variant();
    }
  }
}

void loadGaugeDir(InstVector &ivector, int dir, bool compress12)
{
  string mask;

  PrefetchL1FullGaugeDirIn(ivector, gBase, gOffs, dir, compress12);
  LoadFullGaugeDir(ivector, u_gauge, gBase, gOffs, dir, compress12);

  if (compress12) {
    // printf("Using Compressed Gauges\n");
    for (int c = 0; c < 3; c++) {
      Conj_CrossProd(ivector,
                     u_gauge[2][c],
                     u_gauge[0][(c + 1) % 3],
                     u_gauge[1][(c + 2) % 3],
                     u_gauge[0][(c + 2) % 3],
                     u_gauge[1][(c + 1) % 3],
                     mask);
    }
  }
}

void matMultVec(InstVector &ivector, bool adjMul, int s)
{
  string mask;

  for (int c1 = 0; c1 < 3; c1++) {
    if (!adjMul) {
      mulCVec(ivector,
              ub_spinor[s][c1],
              u_gauge[0][c1],
              b_spinor[s][0],
              mask);
      fmaddCVec(ivector,
                ub_spinor[s][c1],
                u_gauge[1][c1],
                b_spinor[s][1],
                ub_spinor[s][c1],
                mask);
      fmaddCVec(ivector,
                ub_spinor[s][c1],
                u_gauge[2][c1],
                b_spinor[s][2],
                ub_spinor[s][c1],
                mask);
    } else {
      mulConjCVec(ivector,
                  ub_spinor[s][c1],
                  u_gauge[c1][0],
                  b_spinor[s][0],
                  mask);
      fmaddConjCVec(ivector,
                    ub_spinor[s][c1],
                    u_gauge[c1][1],
                    b_spinor[s][1],
                    ub_spinor[s][c1],
                    mask);
      fmaddConjCVec(ivector,
                    ub_spinor[s][c1],
                    u_gauge[c1][2],
                    b_spinor[s][2],
                    ub_spinor[s][c1],
                    mask);
    }
  }
}

void matMultVec(InstVector &ivector, bool adjMul)
{
  matMultVec(ivector, adjMul, 0);
  matMultVec(ivector, adjMul, 1);
}

void dslash_plain_body(InstVector &ivector,
                       bool const compress12,
                       bool const clover,
                       TwistedMassVariant const twisted_mass,
                       bool const isPlus,
                       bool const *const tbc)
{
  declare_b_Spins(ivector);
  declare_ub_Spins(ivector);
  declare_u_gaus(ivector);
  declare_misc(ivector);

  declare_outs(ivector);

  if (clover) {
    // declare temporary spinor, s.t.
    // dout = Dslash * psi
    // out_spinor = clover * dout
    declare_douts(ivector);
  }

  if (clover) {
    if (twisted_mass == TwistedMassVariant::none) {
      declare_clover(ivector);
    } else if (twisted_mass == TwistedMassVariant::degenerate) {
      declare_full_clover(ivector);
    } else if (twisted_mass == TwistedMassVariant::non_degenerate) {
      // BaKo, August 2017:
      // the full ND twisted clover term with mass splitting is implemented
      // as a separate operator which applies the flavour-diagonal
      // and off-diagonal contributions of the inverse clover term
      // after the hopping matrix has been called on both
      // flavours seperately
      //
      // this branch of the code generator is presently never entered
      
      // TODO Here something new for the ND case has to be
      // implemented.
      // Currently this is just copied from the degenerate case.
      declare_full_clover(ivector);
    } else {
      unsupported_twisted_mass_variant();
    }
  }

  FVec(*outspinor)[4][3][2];

  if (clover) {
    outspinor = &dout_spinor;
  } else {
    outspinor = &out_spinor;
  }

  zeroResult(ivector, (*outspinor)[0][0]);

  proj_ops *p_ops;
  recons_ops *rec_ops_bw;
  recons_ops *rec_ops_fw;

  if (isPlus) {
    p_ops = proj_ops_plus;
    rec_ops_bw = rec_plus_pbeta_ops;
    rec_ops_fw = rec_minus_pbeta_ops;
  } else {
    p_ops = proj_ops_minus;
    rec_ops_bw = rec_minus_pbeta_ops;
    rec_ops_fw = rec_plus_pbeta_ops;
  }

  dslash_body(ivector,
              compress12,
              p_ops,
              rec_ops_bw,
              rec_ops_fw,
              *outspinor, tbc);

  if (clover) {
    if (twisted_mass == TwistedMassVariant::none) {
      clover_term(ivector, *outspinor, false);
    } else if (twisted_mass == TwistedMassVariant::degenerate) {
      full_clover_term(ivector, *outspinor, false);
    } else if (twisted_mass == TwistedMassVariant::non_degenerate) {
      // BaKo, August 2017:
      // the full ND twisted clover term with mass splitting is implemented
      // as a separate operator which applies the flavour-diagonal
      // and off-diagonal contributions of the inverse clover term
      // after the hopping matrix has been called on both
      // flavours seperately
      //
      // this branch of the code generator is presently never entered
      full_clover_term(ivector, *outspinor, false);
      // TODO Here something new for the ND case has to be
      // implemented.
      // Currently this is just copied from the degenerate case.
    } else {
      unsupported_twisted_mass_variant();
    }
  } else {
    if (twisted_mass == TwistedMassVariant::none) {
      // Nothing to do.
    } else if (twisted_mass == TwistedMassVariant::degenerate) {
      inverse_twisted_term(ivector, *outspinor, false, isPlus);
    } else if (twisted_mass == TwistedMassVariant::non_degenerate) {
      // BaKo, August 2017:
      // the full ND twisted mass term with mass splitting is implemented
      // as a linalg functor in QPhiX, for now this branch is never
      // entered
      // TODO Here something new for the ND case has to be
      // implemented.
      // Currently this is just copied from the degenerate case.
      inverse_twisted_term(ivector, *outspinor, false, isPlus);
    } else {
      unsupported_twisted_mass_variant();
    }
  }

  // Store
  StreamFullSpinor(ivector, out_spinor, outBase, outOffs);
}

// ***** ------- a chi - b D psi versions

void dslash_achimbdpsi_body(InstVector &ivector,
                            bool const compress12,
                            bool const clover,
                            TwistedMassVariant const twisted_mass,
                            bool const isPlus,
                            bool const *const tbc)
{
  declare_b_Spins(ivector);
  declare_ub_Spins(ivector);
  declare_u_gaus(ivector);
  declare_misc(ivector);

  declare_outs(ivector);
  declare_chi(ivector);

  if (clover) {
    // declare temporary spinor, s.t.
    // dout = Dslash * psi
    // out_spinor = clover * dout
    declare_douts(ivector);
  } else {
    declareFVecFromFVec(ivector, alpha_vec);
    loadBroadcastScalar(ivector, alpha_vec, alpha_name, SpinorType);
  }

  if (clover) {
    // whether we are using twisted mass or not, the clover term in AChiMinusBDPsi
    // is always of the Wilson form. The twisted mass is added via the preconditioning
    // mass below.
    declare_clover(ivector);
  }

  // Fill result with a*chi
  achiResult(ivector, clover, twisted_mass, isPlus);

  // Add twisted preconditioning mass. The following block will
  // generate instructions that compute
  //
  //     out := i \rho \gamma_5 \xhi + out.
  if (clover) {
    declareFVecFromFVec(ivector, prec_mass_rho_vec);
    loadBroadcastScalar(
        ivector, prec_mass_rho_vec, prec_mass_rho_name, SpinorType);
    for (int col = 0; col < 3; col++) {
      for (int spin = 0; spin < 4; spin++) {
        bool const isLower = (spin >= 2);
        FVec *in = chi_spinor[spin][col];
        FVec *out = out_spinor[spin][col];

        // The `\gamma_5` will give a minus sign to the lower spin components.
        // Also if we do the adjoint version, we need a minus for the upper
        // spin components. Therefore we need a plus sign if we do the normal
        // version (`isPlus == true`) and we are on the upper component. The
        // following exclusive or (XOR) operation should take care of both
        // cases.
        if (isPlus ^ isLower) {
          fnmaddFVec(
              ivector, out[RE], prec_mass_rho_vec, in[IM], out[RE]);
          fmaddFVec(
              ivector, out[IM], prec_mass_rho_vec, in[RE], out[IM]);
        } else {
          fmaddFVec(
              ivector, out[RE], prec_mass_rho_vec, in[IM], out[RE]);
          fnmaddFVec(
              ivector, out[IM], prec_mass_rho_vec, in[RE], out[IM]);
        }
      }
    }
  }

  proj_ops *p_ops;
  recons_ops *rec_ops_bw;
  recons_ops *rec_ops_fw;

  if (isPlus) {
    p_ops = proj_ops_plus;
    rec_ops_bw = rec_plus_mbeta_ops;
    rec_ops_fw = rec_minus_mbeta_ops;
  } else {
    p_ops = proj_ops_minus;
    rec_ops_bw = rec_minus_mbeta_ops;
    rec_ops_fw = rec_plus_mbeta_ops;
  }

  dslash_body(ivector,
              compress12,
              p_ops,
              rec_ops_bw,
              rec_ops_fw,
              out_spinor,
              tbc);

  // Store
  StreamFullSpinor(ivector, out_spinor, outBase, outOffs);
}

void pack_face_to_dir_dim_vec(InstVector &ivector,
                              bool isPlus,
                              int dir,
                              int dim)
{
  declare_b_Spins(ivector);
  declare_misc(ivector);

  proj_ops *p_ops =
      (isPlus == true ? proj_ops_plus : proj_ops_minus);
  pack_face_vec(ivector, b_spinor, p_ops, 2 * dim + dir);
}

void recons_add_face_from_dir_dim_vec(InstVector &ivector,
                                      bool compress12,
                                      bool isPlus,
                                      int dir,
                                      int dim,
                                      bool clover,
                                      bool twisted_mass,
                                      bool const use_tbc)
{
  declare_b_Spins(ivector);
  declare_ub_Spins(ivector);
  declare_outs(ivector);
  declare_u_gaus(ivector);
  declare_misc(ivector);

  if (clover || twisted_mass) {
    declare_douts(ivector);
  }

  if (clover && !twisted_mass) {
    declare_clover(ivector);
  } else if (clover && twisted_mass) {
    declare_full_clover(ivector);
  }

  bool isBack = (dir == 0 ? true : false);
  recons_ops *rec_ops;

  if (clover) {
    rec_ops = (isPlus == isBack ? rec_plus_pbeta_ops
                                : rec_minus_pbeta_ops);
  } else {
    rec_ops = (isPlus == isBack ? rec_plus_mbeta_ops
                                : rec_minus_mbeta_ops);
  }

  recons_add_face_vec(ivector,
                      compress12,
                      isBack,
                      rec_ops,
                      dir,
                      dim,
                      clover,
                      twisted_mass,
                      isPlus,
                      use_tbc);
}
