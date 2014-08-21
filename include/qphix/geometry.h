#ifndef QPHIX_GEOMETRY_H
#define QPHIX_GEOMETRY_H

#include "qphix/dslash_utils.h"
#include "qphix/print_utils.h"

#ifdef QPHIX_PARSCALAR
#warning Enabling QMP in Dslash
#undef SEEK_SET
#undef SEEK_END
#undef SEEK_CUR
#include "mpi.h"
#include "qmp.h"

// If this is defined, QMP will be used for geometry information, but explicit
// MPI Isends and Irecvs will be used here with tag 12
//
#define QPHIX_MPI_COMMS_CALLS  

// If this is defined send/recv start/wait in the T-direction will print source, destination and length info.
// #define QMP_DIAGNOSTICS
#endif // QPHIX_PARSCALAR


#include <cstdlib>
#include <iostream>

#ifdef QPHIX_MPI_COMMS_CALLS

// This is a completely arbitrary tag
#define QPHIX_DSLASH_MPI_TAG (12)
#endif


using namespace std;

namespace QPhiX {

  struct CorePhase {
    int Ct;
    int Cyz;
    int startBlock;
    //    char cache_pad[ 64-3*sizeof(int) ]; // Pad to cacheline
  }; 

  typedef unsigned short half;
  
#if defined(QPHIX_MIC_SOURCE)
  float cvtHalf2Float(half val) {
    float ret;
    _mm512_mask_packstorelo_ps(&ret, 0x1, _mm512_mask_extloadunpacklo_ps(_mm512_undefined_ps(), 0x1, &val, _MM_UPCONV_PS_FLOAT16, _MM_HINT_NONE));
    return ret;
  }
  
  half cvtFloat2Half(float val) {
    half ret;
    _mm512_mask_extpackstorelo_ps(&ret, 0x1, _mm512_mask_loadunpacklo_ps(_mm512_undefined_ps(), 0x1, &val), _MM_DOWNCONV_PS_FLOAT16, _MM_HINT_NONE);
    return ret;
  }
  

  // rep: cast 'in' of type T2 into a 'T1' and return it. 
  template <typename T1, typename T2>
  T1 rep(const T2& in) 
  {
    if(sizeof(T1) != sizeof(T2)) {
      
      if(sizeof(T1) == 2) // we are converting float/double to half
	return cvtFloat2Half((float)in);
      else if(sizeof(T2) == 2) // we are converting half to float/double
	return (T1)cvtHalf2Float(in);
      else
	return (T1)in; // we are converting between float and double so just cast is enough
    }
    else {
      return static_cast<T1>(in);  // both T1 and T2 are same
    }
  }
  
#else 

  // rep: cast 'in' of type T2 into a 'T1' and return it. 
  template <typename T1, typename T2>
  T1 rep(const T2& in) 
  {
    return (T1)(in); 
  }
 

#endif


  template<typename T, int V, int S, const bool compressP>
  class Geometry {
  public:
    // Later change this to depend on compressP
    typedef T  FourSpinorBlock[3][4][2][S];
    typedef T  TwoSpinorBlock[3][2][2][V];
    typedef T  SU3MatrixBlock[8][ ( compressP ? 2 : 3 ) ][3][2][V];


    struct CloverBlock {
      T diag1[6][V];         // Real Diagonal part of block 1
      T off_diag1[15][2][V];  // Complex, off diagonal part of block 1
      T diag2[6][V];          // Real Diagonal part of block 2
      T off_diag2[15][2][V];  // Complex, off diagonal part of block 2
    };

    Geometry(const int latt_size[],
	     int By_,
	     int Bz_,
	     int NCores_,
	     int Sy_,
	     int Sz_,
	     int PadXY_,
	     int PadXYZ_,
	     int MinCt_,
	     bool doAlloc_=true)
      : Nd(4),  By(By_), Bz(Bz_), num_cores(NCores_), Sy(Sy_), Sz(Sz_), PadXY(PadXY_), PadXYZ(PadXYZ_), MinCt(MinCt_), nsimt(Sy_*Sz_),  num_threads(NCores_*Sy_*Sz_), doAlloc(doAlloc_)
    {   
      Nx_ = latt_size[0];
      Ny_ = latt_size[1];
      Nz_ = latt_size[2];
      Nt_ = latt_size[3];
      Nxh_ = Nx_/2;
      
      nvecs_ = Nxh()/ S;
      if (Nxh()% S != 0) nvecs_++;
      
      if ( V % S != 0 ) { 
	cerr << "Error: Geometry constructor: SOALEN="<< S <<" does not divide V=" << V << endl;
	abort();
      }
      ngy_ = V/S;
      
      // Padding constants
      Pxy = (nvecs_*Ny_+ PadXY);
      Pxyz = (Pxy*Nz_+ PadXYZ);
      
      
      // Deal with the faces -- in terms of 2 spinors
      NFaceDir_[0] = (Ny_ * Nz_ * Nt_)/2;
      NFaceDir_[1] = (Nx_ * Nz_ * Nt_)/2;
      NFaceDir_[2] = (Nx_ * Ny_ * Nt_)/2;
      NFaceDir_[3] = (Nx_ * Ny_ * Nz_)/2;
            
      // Allos sizes 
      spinor_bytes = (Pxyz * Nt_ + 1)*sizeof(FourSpinorBlock);
      gauge_bytes = ((Pxyz*Nt_*S)/V)*sizeof(SU3MatrixBlock);
      clover_bytes = ((Pxyz*Nt_*S)/V)*sizeof(CloverBlock);

      // Later one will work this out from QMP. While we are testing the buffers, I can set by hand.
      
	  for(int d = 0; d < 4; d++) localDir_[d] = true; 
      
      // This works out the phase breakdown
      int ly = Ny_ / By;
      int lz = Nz_ / Bz;
      int rem = ly * lz;
      int stblk = 0;
      n_phases = 0;
      int n_cores_per_minct = num_cores / MinCt;
      while(rem > 0) {
	int ctd = n_cores_per_minct / rem;
	int ctu = (n_cores_per_minct + rem - 1) / rem;
	CorePhase& p = getCorePhase(n_phases);
	p.Ct = (ctu <= 4 ? ctu : ctd)*MinCt;
	p.Cyz = num_cores / p.Ct;
	if(p.Cyz > rem) p.Cyz = rem;
	p.startBlock = stblk;
	stblk += p.Cyz;
	rem -= p.Cyz;
	//	masterPrintf("Phase %d: Cyz = %d Ct = %d, start = %d\n", n_phases, p.Cyz, p.Ct, p.startBlock);
	n_phases++;
      }
      
#ifdef QPHIX_QMP_COMMS 
      if ( doAlloc) { 
	// We have QMP
	// Decide which directions are local by appealing to 
	// QMP geometry
	const int* machine_size = QMP_get_logical_dimensions();
	if( QMP_get_logical_number_of_dimensions() != 4 ) {
	  QMP_error("Number of QMP logical dimensions must be 4");
	  QMP_abort(1);
	}
	
	for(int d = 0; d < 4; d++)
	  if( machine_size[d] > 1 ) localDir_[d] = false;
	
	myRank = QMP_get_node_number();
	const int *qmp_dims = QMP_get_logical_dimensions();
	const int  *qmp_coords = QMP_get_logical_coordinates();
	int fw_neigh_coords[4];
	int bw_neigh_coords[4];
	
	totalBufSize = 0;
	for(int d = 0; d < 4; d++) {
	  if ( !localDir(d) ) {
	    for(int dim=0; dim < 4; dim++) { 
	      fw_neigh_coords[dim]=bw_neigh_coords[dim]=qmp_coords[dim];
	    }
	    bw_neigh_coords[d]--; if (bw_neigh_coords[d] < 0  ) bw_neigh_coords[d] = qmp_dims[d]-1;
	    fw_neigh_coords[d]++; if (fw_neigh_coords[d] == qmp_dims[d] ) fw_neigh_coords[d] = 0;
	    myNeighboursInDir[2*d+0] = QMP_get_node_number_from( bw_neigh_coords);
	    myNeighboursInDir[2*d+1] = QMP_get_node_number_from( fw_neigh_coords);
	    faceInBytes[d] = NFaceDir(d)*12*sizeof(T); //12 T elements of the half spinor
	    totalBufSize += faceInBytes[d];
	  }
	  else {
	    myNeighboursInDir[2*d+0]= myRank;
	    myNeighboursInDir[2*d+1]= myRank;
	    faceInBytes[d] = 0;
	  }
	}
	totalBufSize *= 4; // 2 bufs for sends & 2 for recvs
	
	/*
	  if(totalBufSize > 0) {
	  commsBuf = (char*)BUFFER_MALLOC(totalBufSize, 4096);
	  if ( commsBuf == 0x0 ) { 
	  cerr << "Couldnt allocate comms buffer" << endl;
	  abort();
	  }
	  }
	  else {
	  commsBuf = NULL;
	  }
	  int tmpptr = 0;
	*/
	for(int d = 0; d < 4; d++) {
	  if(!localDir(d)) {
	    sendToDir[2*d+0]   = (T*)ALIGNED_MALLOC(faceInBytes[d], 4096);
	    sendToDir[2*d+1]   = (T*)ALIGNED_MALLOC(faceInBytes[d], 4096);
	    recvFromDir[2*d+0] = (T*)ALIGNED_MALLOC(faceInBytes[d], 4096);
	    recvFromDir[2*d+1] = (T*)ALIGNED_MALLOC(faceInBytes[d], 4096);
	    
#ifndef QPHIX_MPI_COMMS_CALLS
	    msgmem_sendToDir[2*d+0] = QMP_declare_msgmem(sendToDir[2*d+0], faceInBytes[d]);
	    msgmem_sendToDir[2*d+1] = QMP_declare_msgmem(sendToDir[2*d+1], faceInBytes[d]);
	    msgmem_recvFromDir[2*d+0] = QMP_declare_msgmem(recvFromDir[2*d+0], faceInBytes[d]);
	    msgmem_recvFromDir[2*d+1] = QMP_declare_msgmem(recvFromDir[2*d+1], faceInBytes[d]);
	    
	    mh_sendToDir[2*d+1] = QMP_declare_send_to(msgmem_sendToDir[2*d+1], myNeighboursInDir[2*d+1], 0);
	    mh_recvFromDir[2*d+0] = QMP_declare_receive_from(msgmem_recvFromDir[2*d+0], myNeighboursInDir[2*d+0], 0);		
	    mh_sendToDir[2*d+0] = QMP_declare_send_to(msgmem_sendToDir[2*d+0], myNeighboursInDir[2*d+0], 0);
	    mh_recvFromDir[2*d+1] = QMP_declare_receive_from(msgmem_recvFromDir[2*d+1], myNeighboursInDir[2*d+1], 0);
#endif
	  }
	  else {
	    sendToDir[2*d+0]   = NULL;
	    sendToDir[2*d+1]   = NULL;
	    recvFromDir[2*d+0] = NULL;
	    recvFromDir[2*d+1] = NULL;
#ifndef QPHIX_MPI_COMMS_CALLS
	    msgmem_sendToDir[2*d+0] = NULL;
	    msgmem_sendToDir[2*d+1] = NULL;
	    msgmem_recvFromDir[2*d+0] = NULL;
	    msgmem_recvFromDir[2*d+1] = NULL;
	    
	    mh_sendToDir[2*d+1] = NULL;
	    mh_recvFromDir[2*d+0] = NULL;
	    mh_sendToDir[2*d+0] = NULL;
	    mh_recvFromDir[2*d+1] = NULL;
#endif
	  }
#ifdef QPHIX_MPI_COMMS_CALLS
	  reqSendToDir[2*d+0] = reqRecvFromDir[2*d+0] = MPI_REQUEST_NULL;
	  reqSendToDir[2*d+1] = reqRecvFromDir[2*d+1] = MPI_REQUEST_NULL;
#endif
	} // End loop over dir

	// Determine if I am minimum/maximum in the time direction in the processor grid:
	const int* logical_dimensions = QMP_get_logical_dimensions();
	const int* logical_coordinates = QMP_get_logical_coordinates();
	
	amIPtMin_ = (logical_coordinates[3] == 0);
	amIPtMax_ = (logical_coordinates[3] == (logical_dimensions[3]-1) );
      }
#else
      // On a single node I am always both min and max
      amIPtMin_ = true;
      amIPtMax_ = true;
#endif
    }
    
    
    ~Geometry() {
      
#ifdef QPHIX_QMP_COMMS
	for(int d = 0; d < 4; d++) {
		if(!localDir(d)) {
#ifndef QPHIX_MPI_COMMS_CALLS
			QMP_free_msghandle(mh_sendToDir[2*d+1]);
			QMP_free_msghandle(mh_sendToDir[2*d+0]);
			QMP_free_msghandle(mh_recvFromDir[2*d+1]);
			QMP_free_msghandle(mh_recvFromDir[2*d+0]);
			
			QMP_free_msgmem(msgmem_sendToDir[2*d+1]);
			QMP_free_msgmem(msgmem_sendToDir[2*d+0]);
			QMP_free_msgmem(msgmem_recvFromDir[2*d+1]);
			QMP_free_msgmem(msgmem_recvFromDir[2*d+0]);
#endif
		}
	}

	// Dont free
	// if(commsBuf) BUFFER_FREE( commsBuf, totalBufSize  );
#endif
    }

    inline   int Nxh() const   { return Nxh_; } // Keep
    inline   int Nx()  const   { return Nx_; } // Keep
    inline   int Ny()  const   { return Ny_; } // Keep
    inline   int Nz()  const   { return Nz_; } // Keep
    inline   int Nt()  const   { return Nt_; }  //Keep

    //inline   int NFaceZ() const { return NFaceZ_; }
    //inline   int NFaceT() const { return NFaceT_; }
    inline   int NFaceDir(int d) const { return NFaceDir_[d]; }


    // Part of refactoring -- emulate old behaviour 
    inline   bool localX() const { return localDir_[0]; }
    inline   bool localY() const { return localDir_[1]; }

    inline   bool localZ() const { return localDir_[2]; }
    inline   bool localT() const { return localDir_[3]; }
    inline   bool localDir(int d) const { return localDir_[d]; }

    inline int nVecs() const { return nvecs_; }
    inline int nGY() const { return ngy_; }

    /* Am I the processor with smallest (t=0) in time */
    inline bool amIPtMin() const 
    {
      return amIPtMin_;
    }

    /* Am I the process with the largest (t=P_t - 1)  coordinate in time */ 
    inline bool amIPtMax() const 
    {
      return amIPtMax_;
    }


    /*! \brief Checkerboarded FourSpinor Allocator
     *
     * Allocates a single checkerboard of a Four Spinor.
     * An extra spinor is allocated beyond what is required.
     */
    FourSpinorBlock* allocCBFourSpinor()
    {
            
      FourSpinorBlock *ret_val = (FourSpinorBlock *)BUFFER_MALLOC(spinor_bytes, 64);
      if ( ret_val == (FourSpinorBlock *)0x0 ) { 
	masterPrintf("Failed to allocate FourSpinorBlock\n");
	abort();
      }

      // Zero the field.
      // Cast the pointer.
      T *ret_val_ft = (T *)ret_val;
      
      // change from number of bytes to number of T type elements
      size_t num_ft = spinor_bytes / sizeof(T);
      
      // Zero it all (including) (especially) the pad regions.
      // FIXME: this is not NUMA friendly necessarily
      
#pragma simd 
#pragma vector nontemporal(ret_val_ft)
#pragma omp parallel for
      for(int i=0; i < num_ft; i++) {
	ret_val_ft[i] =rep<T,double>(0.0);
      }
      
      return ret_val+1;
    }


    
    void free(FourSpinorBlock* p) 
    {
      FourSpinorBlock* freeme=p-1;
      BUFFER_FREE(freeme,spinor_bytes);
    }
  
  /*! \brief Checkerboard Gauge Field Allocation
   *
   * This function allocates memory for a single checkerboard of 
   * a gauge field
   */
    SU3MatrixBlock* allocCBGauge()
    {
      SU3MatrixBlock *ret_val = (SU3MatrixBlock *)BUFFER_MALLOC(gauge_bytes, 64);
      if ( ret_val == (SU3MatrixBlock *)0x0 ) { 
	masterPrintf("Failed to allocate SU3MatrixBlock\n");
	abort();
      }

      // For AVX we should loop and zero it here....
      // later on.
      
      // Zero the field.
      // Cast the pointer.
      T *ret_val_ft = (T *)ret_val;
      
      // change from number of bytes to number of T type elements
      size_t num_ft = gauge_bytes / sizeof(T);
      
      // Zero it all (including) (especially) the pad regions.
      // FIXME: this is not NUMA friendly necessarily
#pragma simd
#pragma vector nontemporal(ret_val_ft)
#pragma omp parallel for
      for(int i=0; i < num_ft; i++) {
	ret_val_ft[i] = rep<T,double>(0.0);
      }
      
      return ret_val;
    }

    void free(SU3MatrixBlock *p)
    {
      BUFFER_FREE(p, gauge_bytes);
    }

    CloverBlock* allocCBClov()
    {
      CloverBlock *ret_val = (CloverBlock *)BUFFER_MALLOC(clover_bytes, 64);
      if ( ret_val == (CloverBlock *)0x0 ) { 
	masterPrintf("Failed to allocate CloverBlock\n");
	abort();
      }

      // For AVX we should loop and zero it here....
      // later on.

      // Zero the field.
      // Cast the pointer.
      T *ret_val_ft = (T *)ret_val;

      // change from number of bytes to number of T type elements
      size_t num_ft = clover_bytes / sizeof(T);

      // Zero it all (including) (especially) the pad regions.
      // FIXME: this is not NUMA friendly necessarily
#pragma simd
#pragma vector nontemporal(ret_val_ft)
#pragma omp parallel for
      for(int i=0; i < num_ft; i++) {
	ret_val_ft[i] = rep<T,double>(0.0);
      }

      return ret_val;
    }

    void free(CloverBlock* p) 
    {
      BUFFER_FREE(p,clover_bytes);
    }

    void startSendDir(int d) {
#ifdef QPHIX_QMP_COMMS
#ifndef QPHIX_MPI_COMMS_CALLS
      if( QMP_start(mh_sendToDir[d]) != QMP_SUCCESS ) { 
	QMP_error("Failed to start send\n");
	QMP_abort(1);
      }
#else
      /* **** MPI HERE ******* */
      if (  MPI_Isend( (void *)sendToDir[d], faceInBytes[d/2], MPI_BYTE, myNeighboursInDir[d],  QPHIX_DSLASH_MPI_TAG, MPI_COMM_WORLD, &reqSendToDir[d] ) != MPI_SUCCESS ) { 
	QMP_error("Failed to start send in forward T direction\n");
	QMP_abort(1);
      }
#endif

#ifdef QMP_DIAGNOSTICS
      printf("My Rank: %d, start send dir %d,  My Records: srce=%d, dest=%d len=%d\n", myRank, d, myRank, myNeighboursInDir[d], faceInBytes[d/2]);
#endif
     
#endif // QPHIX_QMP_COMMS
    }
    
    void finishSendDir(int d) {
#ifdef QPHIX_QMP_COMMS

#ifdef QMP_DIAGNOSTICS
      printf("My Rank: %d, finish send dir %d,  My Records: srce=%d, dest=%d len=%d\n", myRank, d, myRank, myNeighboursInDir[d], faceInBytes[d/2]);
#endif

#ifndef QPHIX_MPI_COMMS_CALLS
      if( QMP_wait(mh_sendToDir[d]) != QMP_SUCCESS ) { 
	QMP_error("Failed to finish send\n");
	QMP_abort(1);
      }
#else

      /* **** MPI HERE ******* */
      if (  MPI_Wait(&reqSendToDir[d], MPI_STATUS_IGNORE) != MPI_SUCCESS ) { 
	QMP_error("Wait on send failed \n");
	QMP_abort(1);
      }
#endif
#endif // QMP COMMS
    }

    void startRecvFromDir(int d) { 
#ifdef QPHIX_QMP_COMMS
#ifndef QPHIX_MPI_COMMS_CALLS
      if( QMP_start(mh_recvFromDir[d]) != QMP_SUCCESS ) { 
	QMP_error("Failed to start recv\n");
	QMP_abort(1);
      }
#else
      if ( MPI_Irecv((void *)recvFromDir[d], faceInBytes[d/2], MPI_BYTE, myNeighboursInDir[d], QPHIX_DSLASH_MPI_TAG, MPI_COMM_WORLD, &reqRecvFromDir[d]) != MPI_SUCCESS ) { 
	QMP_error("Recv from dir failed\n");
	QMP_abort(1);
      }
#endif

#ifdef QMP_DIAGNOSTICS
      printf("My Rank: %d, start recv from dir %d,  My Records: srce=%d, dest=%d len=%d\n", myRank, d, myNeighboursInDir[d], myRank,  faceInBytes[d/2]);
#endif
#endif // QPHIX_QMP_COMMS
    }

    void finishRecvFromDir(int d) {
#ifdef QPHIX_QMP_COMMS

#ifdef QMP_DIAGNOSTICS
      printf("My Rank: %d, finish recv from dir %d,  My Records: srce=%d, dest=%d len=%d\n", myRank, d, myNeighboursInDir[d], myRank,  faceInBytes[d/2]);
#endif

#ifndef QPHIX_MPI_COMMS_CALLS
      if( QMP_wait(mh_recvFromDir[d]) != QMP_SUCCESS ) { 
	QMP_error("Failed to finish recv dir\n");
	QMP_abort(1);
      }
#else
      if ( MPI_Wait(&reqRecvFromDir[d], MPI_STATUS_IGNORE) != QMP_SUCCESS ) { 
	QMP_error("Wait on recv from dir failed\n");
	QMP_abort(1);
      }
#endif
#endif // QPHIX_QMP_COMMS
    }

    void progressComms() {
#if 1
#ifdef QPHIX_QMP_COMMS
      int flag = 0;
      MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &flag, MPI_STATUS_IGNORE);
#endif
#endif
    }

#ifdef QPHIX_QMP_COMMS 
    // Send buffers (the actual buffers)
	char *commsBuf;
	
    T* sendToDir[8];

    // Receive buffers (the actual buffers)
    T*  recvFromDir[8];

    // Ranks of the neighbours in the Y, Z and T directions
    int myRank;
    int myNeighboursInDir[8];

    unsigned int faceInBytes[4];
    size_t totalBufSize;

    //  Hack for Karthik here: (Handles for the requests) 
#ifdef QPHIX_MPI_COMMS_CALLS
    MPI_Request reqSendToDir[8];
    MPI_Request reqRecvFromDir[8];
#else
    QMP_msgmem_t msgmem_sendToDir[8];
    QMP_msgmem_t msgmem_recvFromDir[8];
    QMP_msghandle_t mh_sendToDir[8];
    QMP_msghandle_t mh_recvFromDir[8];
#endif // else MPI COMMS_CALLS
#endif // QMP COMMS

    int getBy() const { return By; }
    int getBz() const { return Bz; }
    int getSy() const { return Sy; }
    int getSz() const { return Sz; }
    int getPadXY() const { return PadXY; }
    int getPadXYZ() const { return PadXYZ; }
    int getPxy() const { return Pxy; }
    int getPxyz() const { return Pxyz; }
    int getNSIMT() const { return nsimt; }
    int getNumThreads() const { return num_threads; }
    int getNumCores() const { return num_cores; }
    int getMinCt() const { return MinCt; }
    int getVolCB() const { return Nxh_*Ny_*Nz_*Nt_ ; }
    CorePhase& getCorePhase(int i) { return phase[i]; }
    const CorePhase& getCorePhase(int i) const { return phase[i]; }
    int getNumPhases() const { return n_phases; }
  private:
    

    int Nxh_;
    int Nx_;
    int Ny_;
    int Nz_;
    int Nt_;

    const int Nd;
    const int By;
    const int Bz;
    const int Sy;
    const int Sz;
    const int PadXY;
    const int PadXYZ;
    int Pxy;
    int Pxyz;
    int MinCt; // Minimum no of cores in T dir
    //  MinCt = 1 for single socket/Xeon Phi
    //  MinCt = 2 for dual socket
    //  MinCt = 4 for quad socket

    const int nsimt;
    const int num_threads;
    const int num_cores;

    int NFaceDir_[4];

    bool localDir_[4];

    bool amIPtMin_;
    bool amIPtMax_;

    bool doAlloc;

    int nvecs_;
    int ngy_;

    // Dhiraj's new core mapping
    static const int MAX_PHASES=128;
    CorePhase phase[MAX_PHASES];
    int n_phases;
    int minCt; // Minimum no of cores in T dir

    size_t gauge_bytes;
    size_t spinor_bytes;
    size_t clover_bytes;

    //  minCt = 1 for single socket/Xeon Phi
    //  minCt = 2 for dual socket
    //  minCt = 4 for quad socket
  };
} // Namespace

#endif
