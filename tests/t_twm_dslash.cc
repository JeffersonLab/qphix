#include "unittest.h"
#include "testTWMDslashFull.h"

using namespace QDP;

int main(int argc, char **argv)
{
  TestRunner tests(&argc, &argv);

  tests.addTest(new TestTMDslash(tests.args()), "testDslashFull");

  tests.run();
  tests.summary();
}
