/*
 * qphix_cli_args.h
 *
 *  Created on: Aug 1, 2016
 *      Author: bjoo
 */

#include <string>
#include <cstdlib>

#include "qphix/qphix_cli_args.h"
#include "qphix/print_utils.h"

#include <iostream>

using namespace std;

namespace QPhiX {

void QPhiXCLIArgs::init(int argc, char *argv[])
{
	bool FoundBy=false;
	bool FoundBz=false;
	bool FoundPxy=false;
	bool FoundPxyz=false;
	bool FoundNCores = false;
	bool FoundSy = false;
	bool FoundSz = false;
	bool FoundMinCt = false;

	int i=0;

	while( i < argc)  {

		if ( string(argv[i]).compare("-by") == 0 ) {
			By=atoi(argv[i+1]);
			i+=2;
			FoundBy=true;
	    }
		else if (string(argv[i]).compare("-bz") == 0 ) {
			Bz=atoi(argv[i+1]);
			i+=2;
			FoundBz=true;
	    }
		else if ( string(argv[i]).compare("-c") == 0 ) {
			NCores=atoi(argv[i+1]);
	        i+=2;
	        FoundNCores=true;
	    }
	    else if ( string(argv[i]).compare("-sy") == 0 ) {
	    	Sy=atoi(argv[i+1]);
	    	i+=2;
	    	FoundSy=true;
	    }
	    else if (string(argv[i]).compare("-sz") == 0 ) {
	    	Sz=atoi(argv[i+1]);
	    	i+=2;
	    	FoundSz=true;
	    }
		else if ( string(argv[i]).compare("-pxy") == 0 ) {
			Pxy=atoi(argv[i+1]);
			i+=2;
			FoundPxy=true;
	    }
	    else if (string(argv[i]).compare("-pxyz") == 0 ) {
	    	Pxyz=atoi(argv[i+1]);
	        i+=2;
	        FoundPxyz=true;
	    }
		else if (string(argv[i]).compare("-minct") == 0 ) {
			MinCt=atoi(argv[i+1]);
			i+=2;
			FoundMinCt=true;
	    }
		else {
			i++;
		}

	} // While loop.
	if( !FoundPxy ) { Pxy = 0; FoundPxy=true; }
	if( !FoundPxyz) { Pxyz = 0; FoundPxyz=true;}
	if( !FoundMinCt) { MinCt = 1; FoundMinCt=true;}

	bool FoundAll = FoundBy && FoundBz && FoundPxy && FoundPxyz
			&& FoundNCores && FoundSy && FoundSz & FoundMinCt;

	if (!FoundAll ) { printHelp();abort();}
	initedP = true;
}

void QPhiXCLIArgs::printHelp() const
{
	cout <<"QPhiX Geometry Args: -by BY -bz BZ -c NCores  -sy SY -sz SZ  -pxy Pxy -pxyz Pxyz -minct MinCt\n";
	cout << "   BY is the block size in Y \n" ;
	cout << "   BZ is the block size in Z \n" ;
	cout << "   NCores is the number of cores \n";
	cout << "   Sy is the no of simt threads in y\n";
	cout << "   Sz is the no of simt threads in z\n";
	cout << "   Pxy is the extra pad in the XY plane\n";
	cout << "   Pxyz is the extra pad in the XYZ plane\n";
	cout << "   MinCt is the MinCt in the blocking scheme\n";
}


}
