#include "unittest.h"
#include "testMInvCG.h"

using namespace QDP;

int main(int argc, char **argv)
{
  TestRunner tests(&argc, &argv);

  tests.addTest(new TestMultishift(tests.args()), "MinvCGTester\n");

  tests.run();
  tests.summary();
}
