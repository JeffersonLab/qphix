#include "unittest.h"

#include "testDslashFull.h"

using namespace QDP;

int main(int argc, char **argv)
{
  TestRunner tests(&argc, &argv);

  tests.addTest(new TestDslash(tests.args()), "TestDslash\n");

  tests.run();
  tests.summary();
}
