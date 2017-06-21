#include "unittest.h"
#include "testClovDslashFull.h"

using namespace QDP;

int main(int argc, char **argv)
{
  TestRunner tests(&argc, &argv);

  tests.addTest(new TestClover(tests.args()), "TestClover\n");

  tests.run();
  tests.summary();
}
