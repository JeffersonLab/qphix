#include "unittest.h"
#include "testTWMCloverFull.h"

using namespace QDP;

int main(int argc, char **argv)
{
  TestRunner tests(&argc, &argv);

  tests.addTest(new testTWMCloverFull(tests.args()), "testTWMCloverFull");

  tests.run();
  tests.summary();
}
