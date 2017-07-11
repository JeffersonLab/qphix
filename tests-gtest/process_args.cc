/*
 * process_argc.cpp
 *
 *  Created on: Jul 10, 2017
 *      Author: traveluser
 */

#include "process_args.h"
#include <cstdlib>
#include <string>

namespace QPhiXTesting {

void processArgs(int argc, char *argv[], int* nrow_in)
{
	  int i = 1;
	  while (i < argc) {
	    if (std::string(argv[i]).compare("-x") == 0) {
	      nrow_in[0] = std::atoi(argv[i + 1]);
	      i += 2;
	    } else if (std::string(argv[i]).compare("-y") == 0) {
	      nrow_in[1] = std::atoi(argv[i + 1]);
	      i += 2;
	    } else if (std::string(argv[i]).compare("-z") == 0) {
	      nrow_in[2] = std::atoi(argv[i + 1]);
	      i += 2;
	    } else if (std::string(argv[i]).compare("-t") == 0) {
	      nrow_in[3] = std::atoi(argv[i + 1]);
	      i += 2;
	    }
	    else {
	      i++;
	    }
	  }
}

} // Namespace QPhiX Testing



