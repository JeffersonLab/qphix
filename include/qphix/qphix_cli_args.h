#pragma once

#include <cstdlib>
#include <iostream>

namespace QPhiX
{

class QPhiXCLIArgs
{
 public:
  int getBy() const { return retIfInited(By); }
  int getBz() const { return retIfInited(Bz); }
  int getPxy() const { return retIfInited(Pxy); }
  int getPxyz() const { return retIfInited(Pxyz); }
  int getNCores() const { return retIfInited(NCores); }
  int getSy() const { return retIfInited(Sy); }
  int getSz() const { return retIfInited(Sz); }
  int getMinCt() const { return retIfInited(MinCt); }

  QPhiXCLIArgs() : initedP(false) {}

  int retIfInited(int value) const
  {
    if (!initedP) {
      std::cout << "ERROR: QPhiX CLI Arg Not Inited\n";
      printHelp();
      std::abort();
    } else {
      return value;
    }
  }

  void init(int &argc, char **&argv);
  void printHelp() const;
  void printArgHelp() const;

 private:
  bool initedP;
  int By;
  int Bz;
  int Pxy;
  int Pxyz;
  int NCores;
  int Sy;
  int Sz;
  int MinCt;
};
}
