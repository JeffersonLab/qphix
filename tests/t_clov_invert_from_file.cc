#include "testClovInvertFromFile.h"

using namespace QDP;

int main(int argc, char **argv)
{
  TestRunner tests(&argc, &argv);

  tests.addTest(new TestClovFile(tests.args()), "TestClovFile");

  tests.run();
  tests.summary();
}
