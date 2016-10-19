#include <sstream>
#include <fstream>
#include <iostream>
#include <vector>
#include <typeinfo>
#include <string>

using namespace std;

#include "dslash.h"

#if PRECISION == 1
#if VECLEN == 16
#ifdef AVX512
std::string ARCH_NAME="avx512";
#else
std::string ARCH_NAME="mic";
#endif
#elif VECLEN == 8
#ifdef AVX2
std::string ARCH_NAME="avx2";
#else
std::string ARCH_NAME="avx";
#endif
#elif VECLEN == 4
std::string ARCH_NAME="sse";
#elif VECLEN == 1
std::string ARCH_NAME="scalar";
#endif
#elif PRECISION == 2
#if VECLEN == 8
#ifdef AVX512
std::string ARCH_NAME="avx512";
#else
std::string ARCH_NAME="mic";
#endif
#elif VECLEN == 4
#ifdef AVX2
std::string ARCH_NAME="avx2";
#else
std::string ARCH_NAME="avx";
#endif
#elif VECLEN == 2
std::string ARCH_NAME="sse";
#elif VECLEN == 1
std::string ARCH_NAME="scalar";
#endif
#endif //PRECISION

void mergeIvectorWithL2Prefetches(InstVector& ivector, InstVector& l2prefs);
void dumpIVector(InstVector& ivector, string filename);

string dirname[2] = {"back", "forw"};
string dimchar[4] = {"X", "Y", "Z", "T"};



string basenames[8] = {"xyBase", "xyBase", "xyBase", "xyBase", "zbBase", "zfBase", "tbBase", "tfBase"};
string offsnames[8] = {"xbOffs", "xfOffs", "ybOffs", "yfOffs", "offs",   "offs",   "offs",   "offs"  };
string beta_names[8] = {"coeff_s", "coeff_s", "coeff_s", "coeff_s", "coeff_s", "coeff_s", "coeff_t_b", "coeff_t_f"};

extern FVec beta_vec;
string beta_name("beta");
string alpha_name("alpha");
string outBase("oBase");
string outOffs("offs");
string gBase("gBase");
string gOffs("gOffs");
string chiBase("chiBase");
string chiOffs("offs");
string clBase("clBase");
string clOffs("gOffs");

string mu_name("mu");
string mu_inv_name("muinv");

// Defines which dimensions are involved in SIMD blocking
// Currently just X and Y
bool requireAllOneCheck[4] = {true, true, false, false};
//bool requireAllOneCheck[4] = {false, false, false, false};

void generateFacePackL2Prefetches(InstVector& ivector, int dir)
{
    PrefetchL2HalfSpinorDir(ivector, "outbuf", "hsprefdist", dir, true, 2 /* Ex*/);
    PrefetchL2FullSpinorDirIn(ivector, "xyBase", "offs", "si_prefdist");
}

void generateFaceUnpackL2Prefetches(InstVector& ivector, int dir, bool compress12, bool clover, bool twisted_mass)
{
    PrefetchL2HalfSpinorDir(ivector, "inbuf", "hsprefdist", dir, false, 0 /* None*/);
    PrefetchL2FullGaugeDirIn(ivector, "gBase", "gOffs", dir, "gprefdist", compress12);

    if(clover && !twisted_mass) {
        PrefetchL2FullCloverIn(ivector, "clBase", "gOffs", "clprefdist");
    }
		else if(clover && twisted_mass) {
				// TODO: call the routines for full, non-hermitian, complex clover term
		}

    PrefetchL2FullSpinorDirIn(ivector, outBase, "offs", "soprefdist");
}


// Generate all L2 prefetches
void generateL2Prefetches(InstVector& ivector, bool compress12, bool chi, bool clover, bool twisted_mass)
{
    PrefetchL2FullSpinorDirIn(ivector, "xyBase", "pfyOffs", "siprefdist1");
    //PrefetchL2FullSpinorDirIn(ivector, "pfBase1", "offs", "siprefdist1");
    PrefetchL2FullSpinorDirIn(ivector, "pfBase2", "offs", "siprefdist2");
    PrefetchL2FullSpinorDirIn(ivector, "pfBase3", "offs", "siprefdist3");
    PrefetchL2FullSpinorDirIn(ivector, "pfBase4", "offs", "siprefdist4");

    if(clover && !twisted_mass) {
        PrefetchL2FullCloverIn(ivector, "clBase", "gOffs", "clprefdist");
    }
		else if(clover && twisted_mass) {
				// TODO: call the routines for full, non-hermitian, complex clover term
		}

    if(chi) {
        PrefetchL2FullSpinorDirIn(ivector, "pfBaseChi", "offs", "chiprefdist");
    }

    PrefetchL2FullGaugeIn(ivector, "gBase", "gOffs", "gprefdist", compress12);
    PrefetchL2FullSpinorOut(ivector, outBase, "offs", "siprefdist4");
}



#ifdef SERIAL_SPIN
void dslash_body(InstVector& ivector, bool compress12, proj_ops *ops, recons_ops *rec_ops_bw, recons_ops *rec_ops_fw, FVec outspinor[4][3][2])
{
    for(int dim = 0; dim < 4; dim++) {
        for(int dir = 0; dir < 2; dir++) {
            int d = dim * 2 + dir;
            stringstream d_str;
            d_str << d;
            string mask;
            bool adjMul;
            recons_ops rec_op;

            if(dir == 0) {
                adjMul = true;
                rec_op = rec_ops_bw[dim];
            }
            else {
                adjMul = false;
                rec_op = rec_ops_fw[dim];
            }

            ifStatement(ivector,"accumulate[" + d_str.str() + "]");
            {
                declareFVecFromFVec(ivector, beta_vec);
                loadBroadcastScalar(ivector, beta_vec, beta_names[d], SpinorType);

#ifdef NO_HW_MASKING

                if(requireAllOneCheck[dim]) {
                    ifAllOneStatement(ivector,"accumulate[" + d_str.str() + "]");
                    {
                        for(int s = 0; s < 2; s++) {
                            project(ivector, basenames[d], offsnames[d], ops[d], false, mask, d, s);

                            if(s==0) {
                                loadGaugeDir(ivector, d, compress12);
                            }

                            matMultVec(ivector, adjMul, s);
                            recons_add(ivector, rec_op, outspinor, mask, s);
                        }
                    }
                    elseStatement(ivector);
                }

#endif

                if(requireAllOneCheck[dim]) {
                    mask = "accMask";
                    declareMask(ivector, mask);
                    intToMask(ivector, mask, "accumulate[" + d_str.str() + "]");
                }

                for(int s = 0; s < 2; s++) {
                    project(ivector, basenames[d], offsnames[d], ops[d], false, mask, d, s);

                    if(s==0) {
                        loadGaugeDir(ivector, d, compress12);
                    }

                    matMultVec(ivector, adjMul, s);
                    recons_add(ivector, rec_op, outspinor, mask, s);
                }

#ifdef NO_HW_MASKING

                if(requireAllOneCheck[dim]) {
                    endScope(ivector);
                }

#endif
            }
            endScope(ivector);
        }
    }
}
#else // NO SERIAL_SPIN
void dslash_body(InstVector& ivector, bool compress12, proj_ops *ops, recons_ops *rec_ops_bw, recons_ops *rec_ops_fw, FVec outspinor[4][3][2])
{
    for(int dim = 0; dim < 4; dim++) {
        for(int dir = 0; dir < 2; dir++) {
            int d = dim * 2 + dir;
            stringstream d_str;
            d_str << d;
            string mask;
            bool adjMul;
            recons_ops rec_op;

            adjMul = (dir == 0 ? true : false);
            rec_op = (dir == 0 ? rec_ops_bw[dim] : rec_ops_fw[dim]);

            ifStatement(ivector,"accumulate[" + d_str.str() + "]");
            {
                declareFVecFromFVec(ivector, beta_vec);
                loadBroadcastScalar(ivector, beta_vec, beta_names[d], SpinorType);

#ifdef NO_HW_MASKING

                if(requireAllOneCheck[dim]) {
                    ifAllOneStatement(ivector,"accumulate[" + d_str.str() + "]");
                    {
                        project(ivector, basenames[d], offsnames[d], ops[d], false, mask, d);
                        loadGaugeDir(ivector, d, compress12);
                        matMultVec(ivector, adjMul);
                        recons_add(ivector, rec_op, outspinor, mask);
                    }
                    elseStatement(ivector);
                }

#endif

                if(requireAllOneCheck[dim]) {
                    mask = "accMask";
                    declareMask(ivector, mask);
                    intToMask(ivector, mask, "accumulate[" + d_str.str() + "]");
                }

                project(ivector, basenames[d], offsnames[d], ops[d], false, mask, d);
                loadGaugeDir(ivector, d, compress12);
                matMultVec(ivector, adjMul);
                recons_add(ivector, rec_op, outspinor, mask);
#ifdef NO_HW_MASKING

                if(requireAllOneCheck[dim]) {
                    endScope(ivector);
                }

#endif
            }
            endScope(ivector);
        }
    }
}
#endif // SERIAL_SPIN



// need xyBase, and offs to specify input spinor
// need outbuf for output half spinor
void pack_face_vec(InstVector& ivector, FVec spinor[2][3][2], proj_ops proj[], int dir)
{
    std::string intMask, mask;

    // Check if this dir has mask argument
    if(requireAllOneCheck[dir/2]) {
        intMask = "mask";
        mask = "accMask";
        declareMask(ivector, mask);
        intToMask(ivector, mask, "mask");
    }

    std::string out("outbuf");
    PrefetchL1HalfSpinorDir(ivector, out, dir, true, 2 /*Exclusive*/);

    // We need to reverse direction of projection for our neighbor
    int fb = (dir % 2 == 0 ? 1 : -1);
    project(ivector, "xyBase","offs", proj[dir+fb], true, mask, dir);

    // This will write it to outbuf
    PackHalfSpinor(ivector, spinor, out, dir, intMask);
}

// need inbuf pointer to half spinor
// need gBase and goffs to point to gauge
// need obase and offs to point to spinor to scatter.
void recons_add_face_vec(InstVector& ivector, bool compress12, bool adjMul, recons_ops rops[], int dir, int dim, bool clover, bool twisted_mass, bool isPlus)
{

    std::string in("inbuf");
    std::string mask, intMask;

    extern FVec out_spinor[4][3][2];
    extern FVec dout_spinor[4][3][2];
    extern FVec b_spinor[2][3][2];

    int gauge_index = dim * 2 + dir;

    // Check if this dir has mask argument
    if(requireAllOneCheck[dim]) {
        intMask = "mask";
        mask = "accMask";
        declareMask(ivector, mask);
        intToMask(ivector, mask, "mask");
    }

    declareFVecFromFVec(ivector, beta_vec);
    loadBroadcastScalar(ivector, beta_vec, beta_name, SpinorType);

    FVec (*outspinor)[4][3][2];

    if(clover || twisted_mass) {
        outspinor = &dout_spinor;
        zeroResult(ivector, (*outspinor)[0][0]);
    }
    else {
        outspinor = &out_spinor;
    }

    PrefetchL1HalfSpinorDir(ivector, in, dir, false, 0 /*None*/);
    // Gather in the partial result
    PrefetchL1FullSpinorDirIn(ivector, outBase, outOffs, -1);
    LoadFullSpinor(ivector, out_spinor, outBase, outOffs, "");

    // load b-from inbuf
    UnpackHalfSpinor(ivector, b_spinor, in, gauge_index, intMask);

    loadGaugeDir(ivector, gauge_index, compress12);
    matMultVec(ivector, adjMul);
    recons_add(ivector, rops[dim], *outspinor, mask);

    if(clover && !twisted_mass) clover_term(ivector, *outspinor, true);
		else if(clover && twisted_mass) full_clover_term(ivector, *outspinor, true);
    else if(!clover && twisted_mass) twisted_term(ivector, *outspinor, true, isPlus);
		else if(!clover && !twisted_mass) {};

    // scatter it out
    StoreFullSpinor(ivector, out_spinor, outBase, outOffs);
}


string getTypeName(size_t s)
{
    if(s == 2) {
        return "half";
    }
    else if(s == 4) {
        return "float";
    }
    else if(s == 8) {
        return "double";
    }
    else {
        return "Unknown";
    }
}

void generate_code(void)
{
    if(SOALEN == VECLEN) {
        requireAllOneCheck[1] = false;
    }

#ifdef NO_MASKS
    for(int i = 0; i < 4; i++) {
        requireAllOneCheck[i] = false;
    }
#endif

    const std::string SpinorTypeName = getTypeName(sizeof(SpinorBaseType));
    const std::string GaugeTypeName  = getTypeName(sizeof(GaugeBaseType));
    const std::string CloverTypeName = getTypeName(sizeof(CloverBaseType));


		// DSLASH and DSLASH_ACHIMBDPSI ROUTINES
		// =====================================
		for(auto twisted_mass : {true, false}) {
			for(auto clover : {true, false}) {
				for(auto kernel : {"dslash", "dslash_achimbdpsi"}) {
					for(auto isPlus : {true, false}) {
						for(auto compress12 : {true, false}) {

							if(twisted_mass && clover) continue; // TODO: still to be implemented

							InstVector ivector;
							InstVector l2prefs;
							std::ostringstream filename;

							std::string tmf_prefix  = twisted_mass ? "tmf_" : "";
							std::string clov_prefix = clover ? "clov_"+CloverTypeName+"_" : "";
							std::string plusminus   = isPlus ? "plus" : "minus";
							int num_components = (compress12 ? 12 : 18);
							bool chi_prefetches = false;
							if(kernel == "dslash_achimbdpsi") chi_prefetches = true;

							filename << "./" << ARCH_NAME << "/" << tmf_prefix << clov_prefix << kernel << "_" << plusminus << "_"
								<< "body" << "_" << SpinorTypeName << "_" << GaugeTypeName << "_v" << VECLEN << "_s" << SOALEN << "_" << num_components;

							cout << "GENERATING " << kernel << "_" << plusminus << "_" << "vec body" << endl;

							// Generate instructions
							generateL2Prefetches(l2prefs, compress12, chi_prefetches, clover, twisted_mass);
							if(kernel == "dslash")
								dslash_plain_body(ivector, compress12, clover, twisted_mass, isPlus);
							else if(kernel == "dslash_achimbdpsi")
								dslash_achimbdpsi_body(ivector, compress12, clover, twisted_mass, isPlus);
							mergeIvectorWithL2Prefetches(ivector, l2prefs);
							dumpIVector(ivector, filename.str());

						} // gauge compression
					} // plus/minus
				} // kernel
			} // clover
		} // twisted_mass


		// FACE UNPACK ROUTINES
		// ====================
		for(auto twisted_mass : {true, false}) {
			for(auto clover : {true, false}) {
				for(int dir = 0; dir < 2; dir++) {
					for(int dim = 0; dim < 4; dim++) {
						for(auto isPlus : {true, false}) {
							for(auto compress12 : {true, false}) {

								if(twisted_mass && clover) continue; // TODO: still to be implemented

								InstVector ivector;
								InstVector l2prefs;
								std::ostringstream filename;

								std::string tmf_prefix  = twisted_mass ? "tmf_" : "";
								std::string clov_prefix = clover ? "clov_"+CloverTypeName+"_" : "";
								std::string plusminus   = isPlus ? "plus" : "minus";
								int num_components = (compress12 ? 12 : 18);

								filename << "./" << ARCH_NAME << "/" << tmf_prefix << clov_prefix << "dslash_face_unpack_from_"
									<< dirname[dir] << "_" << dimchar[dim] << "_" << plusminus << "_"
									<< SpinorTypeName << "_" << GaugeTypeName << "_v" << VECLEN << "_s" << SOALEN << "_" << num_components;

								cout << "GENERATING face unpack file " << filename.str() << endl;

								// Generate instructions
								generateFaceUnpackL2Prefetches(l2prefs, 2*dim+dir, compress12, clover, twisted_mass);
								recons_add_face_from_dir_dim_vec(ivector, compress12, isPlus, dir, dim, clover, twisted_mass);
								mergeIvectorWithL2Prefetches(ivector, l2prefs);
								dumpIVector(ivector, filename.str());

							} // gauge compression
						} // plus/minus
					} // dimension
				} // direction
			} // clover
		} // twisted_mass


		// FACE PACK ROUTINES
		// ==================
		for(auto isPlus : {true, false}) {
			for(int dir = 0; dir < 2; dir++) {
				for(int dim = 0; dim < 4; dim++) {
					InstVector ivector;
					InstVector l2prefs;
					std::ostringstream filename;

					string plusminus = isPlus ? "plus" : "minus";

					filename << "./" << ARCH_NAME << "/dslash_face_pack_to_" << dirname[dir] << "_" << dimchar[dim] << "_"
						<< plusminus << "_" << SpinorTypeName << "_" << GaugeTypeName << "_v" << VECLEN << "_s" << SOALEN;

					cout << "GENERATING face pack file " << filename.str() << endl;

					// Generate instructions
					generateFacePackL2Prefetches(l2prefs, 2*dim+dir);
					pack_face_to_dir_dim_vec(ivector, isPlus, dir, dim);
					mergeIvectorWithL2Prefetches(ivector, l2prefs);
					dumpIVector(ivector, filename.str());
				}
			}
		}

    data_types<float,VECLEN,SOALEN,true>::Gauge cmped;
    data_types<float,VECLEN,SOALEN,false>::Gauge uncmped;

    cout << "Compressed Gauge size is " << sizeof(cmped) << endl;
    cout << "Uncompressed Gauge size is " << sizeof(uncmped) << endl;
}
