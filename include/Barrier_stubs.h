#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
using namespace std;


class Barrier {
public:

  Barrier(int num_threads, int threads_per_core) {}
  void init(int tid) { }
  ~Barrier() {}
  void wait(int tid) {}
 private:
  int dummy;
};
