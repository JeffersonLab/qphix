
#include "testDslashNoQDP.h"
#include <omp.h>
#include "cpp_dslash.h"


#include <cstdlib>

using namespace std;
using namespace CPlusPlusWilsonDslash;


#if defined(MIC_SOURCE)
#warning MIC_VECLEN
#define VECLEN 16 
#elif defined(SSE_SOURCE) 
#warning SSE VECLEN
#define VECLEN 4
#include <xmmintrin.h>
#include <mmintrin.h>
#include <smmintrin.h>

#elif defined(AVX_SOURCE) 
#warning AVX VECLEN
#define VECLEN 8
#else
#warning SCALAR_SOURCE
#define VECLEN 1
#endif

#define VECLEN_SCALAR 1

#include <qmp.h>
void
testDslashNoQDP::run(const int lattSize[], const int qmp_geom[]) 
{
  // Work out local lattice size
  int subLattSize[4];
  for(int mu=0; mu < 4; mu++){ 
    subLattSize[mu]=lattSize[mu]/qmp_geom[mu];
  }

  // Work out the size of checkerboarded X-dimension
  int X1h = subLattSize[0]/2;

  // Work out the number of vectors in a 
  int N_blocks_scalar = X1h/VECLEN_SCALAR; if (X1h % VECLEN_SCALAR != 0 ) N_blocks_scalar++;
  N_blocks_scalar *= subLattSize[1]*subLattSize[2]*subLattSize[3];
  if( QMP_is_primary_node() ) { 
    printf("Local volume is %d sites\n", N_blocks_scalar);
  }
  

  // Diagnostic information:
  if( QMP_is_primary_node() ) { 
    printf("Global Lattice Size = ");
    for(int mu=0; mu < 4; mu++){ 
      printf(" %d", lattSize[mu]);
    }
    printf("\n");

    printf("Local Lattice Size = ");
    for(int mu=0; mu < 4; mu++){ 
      printf(" %d", subLattSize[mu]);
    }
    printf("\n");
  

    printf("Block Sizes: By= %d Bz=%d\n", By, Bz);
    printf("Core Grid:  Cy= %d Cz=%d Ct=%d\n", Cy, Cz, Ct);
    int N_cores = Cy*Cz*Ct;
    printf("Cores = %d\n", N_cores);
    printf("Threads_per_core = %d\n", N_simt);

    
  }

  if( QMP_is_primary_node()) { 
    printf("Setting up Scalar Dslash to use as reference\n");
  }

  Dslash<float,VECLEN_SCALAR> D32_scalar(subLattSize, By, Bz, Cy, Cz, Ct, N_simt, compress12);

  // Allocate data for the gauges
  Dslash<float,VECLEN_SCALAR>::SU3MatrixBlockF* packed_gauge_cb0_scalar = (Dslash<float,VECLEN_SCALAR>::SU3MatrixBlockF*)D32_scalar.allocCB(sizeof(Dslash<float,VECLEN_SCALAR>::SU3MatrixBlockF));
  Dslash<float,VECLEN_SCALAR>::SU3MatrixBlockF* packed_gauge_cb1_scalar = (Dslash<float,VECLEN_SCALAR>::SU3MatrixBlockF*)D32_scalar.allocCB(sizeof(Dslash<float,VECLEN_SCALAR>::SU3MatrixBlockF));


  Dslash<float,VECLEN_SCALAR>::SU3MatrixBlockF* u_packed_scalar[2];
  u_packed_scalar[0] = packed_gauge_cb0_scalar;
  u_packed_scalar[1] = packed_gauge_cb1_scalar;

  // Allocate data for the spinors
  Dslash<float,VECLEN_SCALAR>::FourSpinorBlockF* p_even_scalar=(Dslash<float,VECLEN_SCALAR>::FourSpinorBlockF*)D32_scalar.allocCB(sizeof(Dslash<float,VECLEN_SCALAR>::FourSpinorBlockF));
  Dslash<float,VECLEN_SCALAR>::FourSpinorBlockF* p_odd_scalar=(Dslash<float,VECLEN_SCALAR>::FourSpinorBlockF*)D32_scalar.allocCB(sizeof(Dslash<float,VECLEN_SCALAR>::FourSpinorBlockF));
  Dslash<float,VECLEN_SCALAR>::FourSpinorBlockF* c_even_scalar=(Dslash<float,VECLEN_SCALAR>::FourSpinorBlockF*)D32_scalar.allocCB(sizeof(Dslash<float,VECLEN_SCALAR>::FourSpinorBlockF));
  Dslash<float,VECLEN_SCALAR>::FourSpinorBlockF* c_odd_scalar=(Dslash<float,VECLEN_SCALAR>::FourSpinorBlockF*)D32_scalar.allocCB(sizeof(Dslash<float,VECLEN_SCALAR>::FourSpinorBlockF));

 
  // Point to the second block of the array. Now there is padding on both ends.
  Dslash<float,VECLEN_SCALAR>::FourSpinorBlockF* psi_even_scalar=p_even_scalar+1;
  Dslash<float,VECLEN_SCALAR>::FourSpinorBlockF* psi_odd_scalar=p_odd_scalar+1;
  Dslash<float,VECLEN_SCALAR>::FourSpinorBlockF* chi_even_scalar=c_even_scalar+1;
  Dslash<float,VECLEN_SCALAR>::FourSpinorBlockF* chi_odd_scalar=c_odd_scalar+1;
 

  Dslash<float,VECLEN_SCALAR>::FourSpinorBlockF *psi_s_scalar[2] = { psi_even_scalar, psi_odd_scalar };
  Dslash<float,VECLEN_SCALAR>::FourSpinorBlockF *chi_s_scalar[2] = { chi_even_scalar, chi_odd_scalar };

  if( QMP_is_primary_node()) { 
    printf("Filling gauge field: ");
  }

  double start = omp_get_wtime();
#pragma omp parallel for 
  for(int site=0; site < N_blocks_scalar; site++) {

    // 8 directions
    for(int mu=0; mu < 8; mu++) {

      // For a given mu. Fill 2 rows with drand48 noise
      for(int row=0; row < 2; row++) { 
	for(int col=0; col < 3; col++) { 
	  u_packed_scalar[0][site][mu][row][col][RE][0]=drand48();
	  u_packed_scalar[0][site][mu][row][col][IM][0]=drand48();

	  u_packed_scalar[1][site][mu][row][col][RE][0]=drand48();
	  u_packed_scalar[1][site][mu][row][col][IM][0]=drand48();

	}
      } // row
      
      // Normalize the rows
      for(int row=0; row < 2; row++) { 

	double norm_row_cb0 = 0;
	double norm_row_cb1 = 0;

	// Accumulate the norms
	for(int col=0; col < 3; col++) { 
	  norm_row_cb0 += (u_packed_scalar[0][site][mu][row][col][RE][0]
			   *u_packed_scalar[0][site][mu][row][col][RE][0])
	    +(u_packed_scalar[0][site][mu][row][col][IM][0]
	      *u_packed_scalar[0][site][mu][row][col][IM][0]);

	  norm_row_cb1 += (u_packed_scalar[1][site][mu][row][col][RE][0]
			   *u_packed_scalar[1][site][mu][row][col][RE][0])
	    +(u_packed_scalar[1][site][mu][row][col][IM][0]
	      *u_packed_scalar[1][site][mu][row][col][IM][0]);
	}

	norm_row_cb0 = sqrt(norm_row_cb0);
	norm_row_cb1 = sqrt(norm_row_cb1);

	// Normalize each component.
	for(int col=0; col < 3; col++) { 
	  u_packed_scalar[0][site][mu][row][col][RE][0]/=norm_row_cb0;
	  u_packed_scalar[0][site][mu][row][col][IM][0]/=norm_row_cb0;

	  u_packed_scalar[1][site][mu][row][col][RE][0]/=norm_row_cb1;
	  u_packed_scalar[1][site][mu][row][col][IM][0]/=norm_row_cb1;
	}
      }

      {
      // 3rd row reconstruction.
	float ar=u_packed_scalar[0][site][mu][0][0][RE][0];
	float ai=u_packed_scalar[0][site][mu][0][0][IM][0];
	
	float br=u_packed_scalar[0][site][mu][0][1][RE][0];
	float bi=u_packed_scalar[0][site][mu][0][1][IM][0];
	
	float cr=u_packed_scalar[0][site][mu][0][2][RE][0];
	float ci=u_packed_scalar[0][site][mu][0][2][IM][0];
	
	float dr=u_packed_scalar[0][site][mu][1][0][RE][0];
	float di=u_packed_scalar[0][site][mu][1][0][IM][0];
	
	float er=u_packed_scalar[0][site][mu][1][1][RE][0];
	float ei=u_packed_scalar[0][site][mu][1][1][IM][0];
	
	float fr=u_packed_scalar[0][site][mu][1][2][RE][0];
	float fi=u_packed_scalar[0][site][mu][1][2][IM][0];
	
	u_packed_scalar[0][site][mu][2][0][RE][0]=br*fr-bi*fi-er*cr+ei*ci;
	u_packed_scalar[0][site][mu][2][0][IM][0]=er*ci+ei*cr-br*fi-bi*fr;
	u_packed_scalar[0][site][mu][2][1][RE][0]=dr*cr-di*ci-ar*fr+ai*fi;
	u_packed_scalar[0][site][mu][2][1][IM][0]=ar*fi+ai*fr-dr*ci-di*cr;
	u_packed_scalar[0][site][mu][2][2][RE][0]=ar*er-ai*ei-dr*br+di*bi;
	u_packed_scalar[0][site][mu][2][2][IM][0]=dr*bi+di*br-ar*ei-ai*er;
      }

      {
	// 3rd row reconstruction.
	float ar=u_packed_scalar[1][site][mu][0][0][RE][0];
	float ai=u_packed_scalar[1][site][mu][0][0][IM][0];
	
	float br=u_packed_scalar[1][site][mu][0][1][RE][0];
	float bi=u_packed_scalar[1][site][mu][0][1][IM][0];
	
	float cr=u_packed_scalar[1][site][mu][0][2][RE][0];
	float ci=u_packed_scalar[1][site][mu][0][2][IM][0];
	
	float dr=u_packed_scalar[1][site][mu][1][0][RE][0];
	float di=u_packed_scalar[1][site][mu][1][0][IM][0];
	
	float er=u_packed_scalar[1][site][mu][1][1][RE][0];
	float ei=u_packed_scalar[1][site][mu][1][1][IM][0];
	
	float fr=u_packed_scalar[1][site][mu][1][2][RE][0];
	float fi=u_packed_scalar[1][site][mu][1][2][IM][0];
	
	u_packed_scalar[1][site][mu][2][0][RE][0]=br*fr-bi*fi-er*cr+ei*ci;
	u_packed_scalar[1][site][mu][2][0][IM][0]=er*ci+ei*cr-br*fi-bi*fr;
	u_packed_scalar[1][site][mu][2][1][RE][0]=dr*cr-di*ci-ar*fr+ai*fi;
	u_packed_scalar[1][site][mu][2][1][IM][0]=ar*fi+ai*fr-dr*ci-di*cr;
	u_packed_scalar[1][site][mu][2][2][RE][0]=ar*er-ai*ei-dr*br+di*bi;
	u_packed_scalar[1][site][mu][2][2][IM][0]=dr*bi+di*br-ar*ei-ai*er;
      }
      
    }// Mu
    
  } // site
  double end = omp_get_wtime();
  if (QMP_is_primary_node()) printf(" %g sec\n", end - start);


  if (QMP_is_primary_node()) printf("Filling input spinor: ");
  start = omp_get_wtime();

  // Now we need to fill the arrays with drand48 numbers
#pragma omp parallel for
  for(int site=0; site < N_blocks_scalar; site++) {
    for(int spin=0; spin < 4; spin++) { 
      for(int col=0; col < 3; col++) { 
	psi_s_scalar[0][site][spin][col][RE][0]=drand48();
	psi_s_scalar[0][site][spin][col][IM][0]=drand48();
      }
    }
  }

#pragma omp parallel for
  for(int site=0; site < N_blocks_scalar; site++) {
    for(int spin=0; spin < 4; spin++) { 
      for(int col=0; col < 3; col++) { 
	psi_s_scalar[1][site][spin][col][RE][0]=drand48();
	psi_s_scalar[1][site][spin][col][IM][0]=drand48();
      }
    }
  }

  end = omp_get_wtime();
  if( QMP_is_primary_node() ) printf(" %g sec \n", end -start);


  int Nblocks = X1h/VECLEN; if (X1h % VECLEN != 0 ) Nblocks++;
  Nblocks *= subLattSize[1]*subLattSize[2]*subLattSize[3];
  if ( QMP_is_primary_node() ) { 
    printf("Nblocks= %d\n", Nblocks);
  }


  if( QMP_is_primary_node()) printf("Creating Vectorized Dslash\n");
	  // Create Scalar Dslash Class
  Dslash<float,VECLEN> D32(subLattSize, By, Bz, Cy, Cz, Ct, N_simt, compress12);

  // Allocate data for the gauges
  Dslash<float,VECLEN>::SU3MatrixBlockF* packed_gauge_cb0 = (Dslash<float,VECLEN>::SU3MatrixBlockF*)D32.allocCB(sizeof(Dslash<float,VECLEN>::SU3MatrixBlockF));
  Dslash<float,VECLEN>::SU3MatrixBlockF* packed_gauge_cb1 = (Dslash<float,VECLEN>::SU3MatrixBlockF*)D32.allocCB(sizeof(Dslash<float,VECLEN>::SU3MatrixBlockF));


  Dslash<float,VECLEN>::SU3MatrixBlockF* u_packed[2];
  u_packed[0] = packed_gauge_cb0;
  u_packed[1] = packed_gauge_cb1;

  // Allocate data for the spinors
  Dslash<float,VECLEN>::FourSpinorBlockF* p_even=(Dslash<float,VECLEN>::FourSpinorBlockF*)D32.allocCB(sizeof(Dslash<float,VECLEN>::FourSpinorBlockF));
  Dslash<float,VECLEN>::FourSpinorBlockF* p_odd=(Dslash<float,VECLEN>::FourSpinorBlockF*)D32.allocCB(sizeof(Dslash<float,VECLEN>::FourSpinorBlockF));
  Dslash<float,VECLEN>::FourSpinorBlockF* c_even=(Dslash<float,VECLEN>::FourSpinorBlockF*)D32.allocCB(sizeof(Dslash<float,VECLEN>::FourSpinorBlockF));
  Dslash<float,VECLEN>::FourSpinorBlockF* c_odd=(Dslash<float,VECLEN>::FourSpinorBlockF*)D32.allocCB(sizeof(Dslash<float,VECLEN>::FourSpinorBlockF));

 
  // Point to the second block of the array. Now there is padding on both ends.
  Dslash<float,VECLEN>::FourSpinorBlockF* psi_even=p_even+1;
  Dslash<float,VECLEN>::FourSpinorBlockF* psi_odd=p_odd+1;
  Dslash<float,VECLEN>::FourSpinorBlockF* chi_even=c_even+1;
  Dslash<float,VECLEN>::FourSpinorBlockF* chi_odd=c_odd+1;
 

  Dslash<float,VECLEN>::FourSpinorBlockF *psi_s[2] = { psi_even, psi_odd };
  Dslash<float,VECLEN>::FourSpinorBlockF *chi_s[2] = { chi_even, chi_odd };


  // Now we have to pack.
  int Ny = subLattSize[1];
  int Nz = subLattSize[2];
  int Nt = subLattSize[3];
  int Nlines = Ny*Nz*Nt;
  int Nvecs = X1h/VECLEN;
  if (X1h % VECLEN != 0 ) Nvecs++;

  if ( QMP_is_primary_node()) printf("Repacking gauge field in vector format: ");
  start = omp_get_wtime();
  for(int line=0; line < Nlines; line++) {
    int startBlock = Nvecs*line;
    int startBlockScalar = X1h*line;

    for(int vec=0; vec < Nvecs; vec++) {
      int block = startBlock+vec;
    
      
      for(int mu=0; mu < 8; mu++) { 
	for(int row=0; row < 3; row++) { 
	  for(int col=0; col < 3; col++) { 
	    
	    for(int site =0 ; site < VECLEN; site++) { 
	    
	      int blockScalar = site + vec*VECLEN + startBlockScalar;

	      int site_in_line = site + vec*VECLEN;
	      if( site_in_line < X1h ) { 
		packed_gauge_cb0[block][mu][row][col][RE][site] = 
		  u_packed_scalar[0][blockScalar][mu][row][col][RE][0];
		
		packed_gauge_cb0[block][mu][row][col][IM][site] = 
		  u_packed_scalar[0][blockScalar][mu][row][col][IM][0];
		
		packed_gauge_cb1[block][mu][row][col][RE][site] = 
		  u_packed_scalar[1][blockScalar][mu][row][col][RE][0];
		
		packed_gauge_cb1[block][mu][row][col][IM][site] = 
		  u_packed_scalar[1][blockScalar][mu][row][col][IM][0];
	      }
	      else { 
		packed_gauge_cb0[block][mu][row][col][RE][site] = 0;
		packed_gauge_cb0[block][mu][row][col][IM][site] = 0;
		packed_gauge_cb1[block][mu][row][col][RE][site] = 0;
		packed_gauge_cb1[block][mu][row][col][IM][site] = 0;
	      }
	       
	    }
	  }
	}
      }
    }
  }
  end=omp_get_wtime();
  if(QMP_is_primary_node()) printf(" %g sec\n", end-start);


  if(QMP_is_primary_node()) printf("Repacking Input Spinor in vector format ");
  start = omp_get_wtime();

  for(int line=0; line < Nlines; line++) {
    int startBlock = Nvecs*line;
    int startBlockScalar = X1h*line;

    for(int vec=0; vec < Nvecs; vec++) {
      int block = startBlock+vec;

      for(int spin=0; spin < 4; spin++) { 
	for(int col=0; col < 3; col++) { 
	  
	  for(int site =0; site < VECLEN; site++) { 
	    int blockScalar = site + vec*VECLEN + startBlockScalar;
	    int site_in_line = site + vec*VECLEN;
	    if( site_in_line < X1h ) { 
	      psi_s[0][block][spin][col][RE][site] = psi_s_scalar[0][blockScalar][spin][col][RE][0];
	      psi_s[0][block][spin][col][IM][site] = psi_s_scalar[0][blockScalar][spin][col][IM][0];
	      
	      psi_s[1][block][spin][col][RE][site] = psi_s_scalar[1][blockScalar][spin][col][RE][0];
	      psi_s[1][block][spin][col][IM][site] = psi_s_scalar[1][blockScalar][spin][col][IM][0];
	    }
	    else { 
	      psi_s[0][block][spin][col][RE][site] = 0;
	      psi_s[0][block][spin][col][IM][site] = 0;
	      
	      psi_s[1][block][spin][col][RE][site] = 0;
	      psi_s[1][block][spin][col][IM][site] = 0;
	    }
	  }
	}
      }
    }
  }
  end = omp_get_wtime();
  if( QMP_is_primary_node()) printf(" %g sec\n", end - start);

  // Go through the test cases -- apply SSE dslash versus, QDP Dslash 
  for(int isign=1; isign >= -1; isign -=2) {
    for(int cb=0; cb < 2; cb++) { 
      int source_cb = 1 - cb;
      int target_cb = cb;

      if( QMP_is_primary_node()) printf("Zeroing vector output spinor: ");
      start = omp_get_wtime();

#pragma omp parallel for
      for(int block=0; block < Nblocks; block++) { 
	for(int spin=0; spin < 4; spin++) { 
	  for(int color=0; color < 3; color++) { 
	    for(int site=0; site < VECLEN; site++) { 
	      chi_s[0][block][spin][color][0][site] = 0; 
	      chi_s[0][block][spin][color][1][site] = 0;
	      
	      chi_s[1][block][spin][color][0][site] = 0;
	      chi_s[1][block][spin][color][1][site] = 0;
	    }
	  }
	}
      }
      end = omp_get_wtime();
      if (QMP_is_primary_node()) printf( " %g sec\n", end-start);


      if (QMP_is_primary_node()) printf( "Zeroing scalar output spinor: ");
      start = omp_get_wtime();
#pragma omp parallel for
      for(int block=0; block < N_blocks_scalar; block++) { 
	for(int spin=0; spin < 4; spin++) { 
	  for(int color=0; color < 3; color++) { 
	      chi_s_scalar[0][block][spin][color][0][0] = 0; 
	      chi_s_scalar[0][block][spin][color][1][0] = 0;
	      
	      chi_s_scalar[1][block][spin][color][0][0] = 0;
	      chi_s_scalar[1][block][spin][color][1][0] = 0;
	  }
	}
      }
      end = omp_get_wtime();
      if ( QMP_is_primary_node()) printf(" %g sec\n", end - start);
     

      if( QMP_is_primary_node()) printf("Applying Vector Dslash: ");
      start = omp_get_wtime();
      // Apply Optimized Dslash
      D32.dslash((float *)(chi_s[target_cb]),	
		 (float *)(psi_s[source_cb]),
		 (void *)(u_packed[target_cb]),
		 isign, 
		 target_cb);
      end=omp_get_wtime();
      if( QMP_is_primary_node()) printf(" %g sec\n", end - start);

      if( QMP_is_primary_node()) printf("Applying scalar Dslash for reference: ");
      start = omp_get_wtime();

      D32_scalar.dslash((float *)(chi_s_scalar[target_cb]),
	  (float *)(psi_s_scalar[source_cb]),
	  (void *)(u_packed_scalar[target_cb]),
	  isign,
	  target_cb);

      end = omp_get_wtime();
      
      if( QMP_is_primary_node()) printf(" %g sec\n", end - start);

      if( QMP_is_primary_node()) printf("Diffing...");
      start = omp_get_wtime();
      double norm_diff = 0;

#pragma omp parallel for reduction(+:norm_diff)
      for(int line=0; line < Nlines; line++) {
	int startBlock = Nvecs*line;
	int startBlockScalar = X1h*line;

	for(int vec=0; vec < Nvecs; vec++) {
	  int block = startBlock+vec;

	  for(int spin=0; spin < 4; spin++) { 
	    for(int col=0; col < 3; col++) { 
	      for(int cmpx = 0; cmpx < 2; cmpx++) { 

		for(int site =0; site < VECLEN; site++) { 
		  int blockScalar = site + vec*VECLEN + startBlockScalar;
		  int site_in_line= site+vec*VECLEN;
		  double diff = (chi_s[target_cb][block][spin][col][cmpx][site]
				  - chi_s_scalar[target_cb][blockScalar][spin][col][cmpx][0]);
		  norm_diff += diff*diff;
		}
	      }
	    }
	  }
	}
      }
      end = omp_get_wtime();
      if( QMP_is_primary_node()) printf(" %g sec \n", end - start);

      QMP_sum_double(&norm_diff);
      norm_diff=sqrt(norm_diff);
      norm_diff /= (double)(4*3*X1h*Ny*Nz*Nt);
      if(QMP_is_primary_node()) { 
	printf( "\t isign=%d cb=%d : Norm Diff = %e\n", isign, cb, norm_diff);
      }


      if( norm_diff > 1.0e-7 ) { 
#pragma omp parallel for
	for(int line=0; line < Nlines; line++) {
	  int startBlock = Nvecs*line;
	  int startBlockScalar = X1h*line;
	  
	  for(int vec=0; vec < Nvecs; vec++) {
	    int block = startBlock+vec;
	    
	    for(int spin=0; spin < 4; spin++) { 
	      for(int col=0; col < 3; col++) { 
		for(int cmpx = 0; cmpx < 2; cmpx++) { 
		  
		  for(int site =0; site < VECLEN; site++) { 
		    int blockScalar = site + vec*VECLEN + startBlockScalar;
		    int site_in_line= site+vec*VECLEN;
		    if( site_in_line < X1h ) { 
		      double diff = fabs(chi_s[target_cb][block][spin][col][cmpx][site]
					 - chi_s_scalar[target_cb][blockScalar][spin][col][cmpx][0]);
		      
		      if (diff > 5e-6) { 
			printf("Norm diff = %e, block=%d site=%d line=%d vec=%d isign=%d cb=%d, \n", 
			       diff, block, site, line,vec, isign, cb);
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
	
    }
	
	
    }
  }

  ALIGNED_FREE(packed_gauge_cb0_scalar);
  ALIGNED_FREE(packed_gauge_cb1_scalar);
  ALIGNED_FREE(p_even_scalar);
  ALIGNED_FREE(p_odd_scalar);
  ALIGNED_FREE(c_even_scalar);
  ALIGNED_FREE(c_odd_scalar);

  ALIGNED_FREE(packed_gauge_cb0);
  ALIGNED_FREE(packed_gauge_cb1);
  ALIGNED_FREE(p_even);
  ALIGNED_FREE(p_odd);
  ALIGNED_FREE(c_even);
  ALIGNED_FREE(c_odd);
						   

}
