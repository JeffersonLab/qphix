#ifndef QPHIX_QDP_PACKER_QDPJIT_H
#define QPHIX_QDP_PACKER_QDPJIT_H

#warning "building with QDPJIT packers"


#include "qdp.h"
#include "qphix/geometry.h"

#include "qphix/dslash_def.h"
#include "qphix/qphix_config.h"

#if defined(QPHIX_MIC_SOURCE) || defined(QPHIX_AVX512_SOURCE)
#include <immintrin.h>
#endif

#ifdef QPHIX_BUILD_CLOVER
#include "qphix/clover_dslash_def.h"
#endif

using namespace QDP;

namespace QPhiX { 


  template<typename FT, int veclen, int soalen, bool compress, typename QDPGauge>
    void qdp_pack_gauge(const QDPGauge& u,  
      typename Geometry<FT,veclen,soalen,compress>::SU3MatrixBlock *u_cb0, 
      typename Geometry<FT,veclen,soalen,compress>::SU3MatrixBlock *u_cb1, 
      Geometry<FT,veclen,soalen,compress>& s)
  {
    // Get the subgrid latt size.
    int Nt = s.Nt();
    int Nz = s.Nz();
    int Ny = s.Ny();
    int nvecs = s.nVecs();
    int nyg = s.nGY();
    int Pxy = s.getPxy();
    int Pxyz = s.getPxyz();
	int qdp_inner_dim = (int)getDataLayoutInnerSize();
	const int n_complex_dim = 2;

    // Shift the lattice to get U(x-mu)
    QDPGauge u_minus(4);
    for(int mu=0; mu < 4; mu++) {
      u_minus[mu] = shift(u[mu], BACKWARD, mu);
    }

    // This ought to be the underlying precision type (float/double) in use
    // by QDP-JIT and may well be different from FT
    // So casting will be needed
    using WT = typename WordType<QDPGauge>::Type_t;

#pragma omp parallel for collapse(4)
		for(int t = 0; t < Nt; t++) {
			for(int z = 0; z < Nz; z++) {
				for(int y = 0; y < Ny; y++) {
					for(int s = 0; s < nvecs; s++) {
						for(int mu = 0; mu < 4; mu++) {

						    const WT* u_ptr = (const WT*)((u[mu]).getFjit());
							const WT* u_minus_ptr = (const WT*)((u_minus[mu]).getFjit());

							int outer_c = 3;
							if ( compress ) {
								outer_c = 2;
							}
							for(int c = 0; c < outer_c; c++) {
								for(int c2 = 0; c2 < 3; c2++) {
									for(int x = 0; x < soalen; x++) {


										int block = (t*Pxyz+z*Pxy)/nyg+(y/nyg)*nvecs+s;
										int xx = (y%nyg)*soalen+x;

										// QDP-JIT in OCSRI layout:
										//  Outer x Color x Spin x Real x Inner
										//
										// However, it hardly matters maybe because
										// There is no spin on the gauge fields

										int qdpsite = x + soalen*(s + nvecs*(y + Ny*(z + Nz*t)));

										int qdpsite_rb0 = qdpsite+rb[0].start();
										int qdp_outer_rb0 = qdpsite_rb0 / qdp_inner_dim;
										int qdp_inner_rb0 = qdpsite_rb0 % qdp_inner_dim;
										int rb0_offset = qdp_inner_rb0 + qdp_inner_dim*n_complex_dim*Nc*Nc*qdp_outer_rb0;

										u_cb0[block][2*mu][c][c2][RE][xx] = (FT)u_minus_ptr[rb0_offset + qdp_inner_dim*(RE+n_complex_dim*(c + Nc*c2)) ];
										u_cb0[block][2*mu][c][c2][IM][xx] = (FT)u_minus_ptr[rb0_offset + qdp_inner_dim*(IM+n_complex_dim*(c + Nc*c2)) ];
										u_cb0[block][2*mu+1][c][c2][RE][xx] = (FT)u_ptr[rb0_offset + qdp_inner_dim*(RE+n_complex_dim*(c + Nc*c2)) ];
										u_cb0[block][2*mu+1][c][c2][IM][xx] = (FT)u_ptr[rb0_offset + qdp_inner_dim*(IM+n_complex_dim*(c + Nc*c2)) ];

										int qdpsite_rb1 = qdpsite+rb[1].start();
										int qdp_outer_rb1 = qdpsite_rb1 / qdp_inner_dim;
										int qdp_inner_rb1 = qdpsite_rb1 % qdp_inner_dim;
										int rb1_offset = qdp_inner_rb1 + qdp_inner_dim*n_complex_dim*Nc*Nc*qdp_outer_rb1;

										u_cb1[block][2*mu][c][c2][RE][xx] = (FT)u_minus_ptr[rb1_offset + qdp_inner_dim*(RE+n_complex_dim*(c + Nc*c2)) ];
										u_cb1[block][2*mu][c][c2][IM][xx] = (FT)u_minus_ptr[rb1_offset + qdp_inner_dim*(IM+n_complex_dim*(c + Nc*c2)) ];
										u_cb1[block][2*mu+1][c][c2][RE][xx] = (FT)u_ptr[rb1_offset + qdp_inner_dim*(RE+n_complex_dim*(c + Nc*c2)) ];
										u_cb1[block][2*mu+1][c][c2][IM][xx] = (FT)u_ptr[rb1_offset + qdp_inner_dim*(IM+n_complex_dim*(c + Nc*c2)) ];

									}
								}
							}
						}
					}
				}
			}
		}


	}

#ifdef QPHIX_BUILD_CLOVER

	// This accesses the Internals of the LLVMCloverTerm
	template<typename FT, int veclen, int soalen, bool compress, typename ClovTerm>
	void qdp_pack_clover(const ClovTerm& qdp_clov_in,
			typename ClovDslash<FT,veclen,soalen,compress>::CloverBlock* cl_out,Geometry<FT,veclen,soalen,compress>& s, int cb)
	{
		// Get the subgrid latt size.
		int Nt = s.Nt();
		int Nz = s.Nz();
		int Ny = s.Ny();
		int nvecs = s.nVecs();
		int nyg = s.nGY();
		int Pxy = s.getPxy();
		int Pxyz = s.getPxyz();

		// Sanity Check
		// QDP Type is
		// Outer x 2 chiral blocks x 6 floats x inner sites
		// JIT-LLVM CLOVER will need to expose DiagType and OffDiagType
		// It does in the testcase here but make sure chroma version does it too.

		using DiagType = typename ClovTerm::DiagType;
		using OffDiagType = typename ClovTerm::OffDiagType;

		// This is the base type used here. Either double or float
		using WT = typename WordType<DiagType>::Type_t;

		const DiagType& diag_term = qdp_clov_in.getDiagBuffer();
		const OffDiagType& off_diag_term = qdp_clov_in.getOffDiagBuffer();

		const WT* diag_buf = (const WT*)diag_term.getFjit();
		const WT* off_diag_buf = (const WT*)off_diag_term.getFjit();

		int qdp_inner_dim = (int)getDataLayoutInnerSize();
		const int n_complex_dim = 2;
		const int diag_dim = 6;
		const int offdiag_dim = 15;
		const int chiral_dim = 2;
		const int diag_offset = qdp_inner_dim*diag_dim*chiral_dim;
		const int offdiag_offset = qdp_inner_dim*offdiag_dim*n_complex_dim*chiral_dim;

		// No elem calls in parallel region
#pragma omp parallel for collapse(4)
		for(int t = 0; t < Nt; t++) {
			for(int z = 0; z < Nz; z++) {
				for(int y = 0; y < Ny; y++) {
					for(int s = 0; s < nvecs; s++) {
						for(int x = 0; x < soalen; x++) {

							int block = (t*Pxyz+z*Pxy)/nyg+(y/nyg)*nvecs+s;
							int xx = (y%nyg)*soalen+x;

							int qdpsite = x + soalen*(s + nvecs*(y + Ny*(z + Nz*t)))+rb[cb].start();
							int qdp_outer_idx = qdpsite / qdp_inner_dim;
							int qdp_inner_idx = qdpsite % qdp_inner_dim;

							const WT* diag_base = &diag_buf[diag_offset*qdp_outer_idx];
							const WT* offdiag_base = &off_diag_buf[offdiag_offset*qdp_outer_idx];

							// WARNING: THIS WORKS ONLY IN OCSRI Layout in QDP_JIT
							//
							// For OSCRI Layout: chiral component and diagonal must be transposed.
							//                   chiral component and off diagonal must also be transposed

							// Logical Outer x Spin x Color x Inner => Outer x Comp x Diag x Inner
							// But we are using OCSRI => Outer Diag Comp Inner

							for(int d=0; d < 6; d++) {
								cl_out[block].diag1[d][xx] = (FT)(diag_base[ qdp_inner_idx + qdp_inner_dim*(0 + chiral_dim*d ) ]);
							}

							for(int od=0; od < 15; od++) {

								cl_out[block].off_diag1[od][RE][xx] = (FT)(offdiag_base[ qdp_inner_idx +
										qdp_inner_dim*(RE + n_complex_dim*(0 + chiral_dim*od))]);

								cl_out[block].off_diag1[od][IM][xx] = (FT)(offdiag_base[ qdp_inner_idx +
										qdp_inner_dim*(IM + n_complex_dim*(0 + chiral_dim*od))]);

							}

							for(int d=0; d < 6; d++) {
								cl_out[block].diag2[d][xx] = (FT)(diag_base[ qdp_inner_idx + qdp_inner_dim*(1 + chiral_dim*d ) ]);
							}

							for(int od=0; od < 15; od++) {
								cl_out[block].off_diag2[od][RE][xx] = (FT)(offdiag_base[ qdp_inner_idx +
										qdp_inner_dim*(RE + n_complex_dim*(1 + chiral_dim*od))]);

								cl_out[block].off_diag2[od][IM][xx] = (FT)(offdiag_base[ qdp_inner_idx +
										qdp_inner_dim*(IM + n_complex_dim*(1 + chiral_dim*od))]);

							}
						}
					}
				}
			}
		}
	}



#endif  // IFDEF BUILD CLOVER


	template<typename FT, int veclen, int soalen, bool compress, typename QDPSpinor>
	void qdp_pack_cb_spinor(const QDPSpinor& psi_in,
			typename Geometry<FT,veclen,soalen, compress>::FourSpinorBlock* psi,
			Geometry<FT,veclen,soalen,compress>& s,
			int cb)
	{
		// Get the subgrid latt size.
		int Nt = s.Nt();
		int Nz = s.Nz();
		int Ny = s.Ny();
		int Nxh = s.Nxh();
		int nvecs = s.nVecs();
		int Pxy = s.getPxy();
		int Pxyz = s.getPxyz();

		int qdp_inner_dim = (int)getDataLayoutInnerSize();
		const int n_complex_dim = 2;
		using WT =typename WordType<QDPSpinor>::Type_t;

		const int outer_block_size = qdp_inner_dim*n_complex_dim*Nc*Ns;

#pragma omp parallel for collapse(4)
		for(int t=0; t < Nt; t++) {
			for(int z=0; z < Nz; z++) {
				for(int y=0; y < Ny; y++) {
					for(int s=0; s < nvecs; s++) {
						for(int col=0; col < 3; col++) {
							for(int spin=0; spin < 4; spin++) {
								for(int x=0; x < soalen; x++) {

									int ind = t*Pxyz+z*Pxy+y*nvecs+s; //((t*Nz+z)*Ny+y)*nvecs+s;
									int x_coord = s*soalen + x;

									int qdp_index = ((t*Nz + z)*Ny + y)*Nxh + x_coord + rb[cb].start();
									int qdp_inner_idx = qdp_index % qdp_inner_dim;
									int qdp_outer_idx = qdp_index / qdp_inner_dim;

									WT* psiptr=(WT *)(psi_in.getFjit());

									// ASSUME OCSRI order. So the fixed points our
									// qdp_outer*block-size + qdp_inner
									// and the iteration is over the spinor as Inner x Complex x Spin x Color

									int offset =  qdp_inner_idx + outer_block_size*qdp_outer_idx;
									psi[ind][col][spin][0][x] = (FT)psiptr[offset+qdp_inner_dim*(RE + n_complex_dim*(spin + Ns*col)) ];

									psi[ind][col][spin][1][x] = (FT)psiptr[offset+qdp_inner_dim*(IM + n_complex_dim*(spin + Ns*col)) ];
								}
							}
						}
					}
				}
			}
		}

	}


  template<typename FT, int veclen, int soalen, bool compress, typename QDPSpinor>
    void qdp_unpack_cb_spinor(typename Geometry<FT,veclen,soalen,compress>::FourSpinorBlock* chi_packed, 
            QDPSpinor& chi,
            Geometry<FT,veclen,soalen,compress>& s,
            int cb) 
  { 
    int Nt = s.Nt();
    int Nz = s.Nz();
    int Ny = s.Ny();
    int Nxh = s.Nxh();
    int nvecs = s.nVecs();
    int Pxy = s.getPxy();
    int Pxyz = s.getPxyz();

    int qdp_inner_dim = (int)getDataLayoutInnerSize();
    const int n_complex_dim = 2;
    using WT =typename WordType<QDPSpinor>::Type_t;

    const int outer_block_size = qdp_inner_dim*n_complex_dim*Nc*Ns;


#pragma omp parallel for collapse(4)
		for(int t=0; t < Nt; t++) {
			for(int z=0; z < Nz; z++) {
				for(int y=0; y < Ny; y++) {
					for(int s=0; s < nvecs; s++) {
						for(int spin=0; spin < 4; spin++) {
							for(int col=0; col < 3; col++) {
								for(int x=0; x < soalen; x++) {

									int ind = t*Pxyz+z*Pxy+y*nvecs+s; //((t*Nz+z)*Ny+y)*nvecs+s;
									int x_coord = s*soalen + x;

									int qdp_index = ((t*Nz + z)*Ny + y)*Nxh + x_coord + rb[cb].start();
									int qdp_inner_idx = qdp_index % qdp_inner_dim;
									int qdp_outer_idx = qdp_index / qdp_inner_dim;

									WT* chiptr=(WT *)(chi.getFjit());
									// ASSUME OCSRI order. So the fixed points our
									// qdp_outer*block-size + qdp_inner
									// and the iteration is over the spinor as Inner x Complex x Spin x Color

									int offset = qdp_inner_idx + outer_block_size*qdp_outer_idx;

									chiptr[offset+qdp_inner_dim*(RE + n_complex_dim*(spin + Ns*col)) ] = (WT)chi_packed[ind][col][spin][0][x];
									chiptr[offset+qdp_inner_dim*(IM + n_complex_dim*(spin + Ns*col)) ] = (WT)chi_packed[ind][col][spin][1][x];

								}
							}
						}
					}
				}
			}
		}
	}

};


#endif
