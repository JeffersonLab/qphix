#include "unittest.h"
#include "testClovInvertFromFile.h"

using namespace QDP;

int main(int argc, char **argv)
{
  TestRunner tests(&argc, &argv);

  tests.addTest(new testClovInvertFromFile(tests.args()),
                "testClovInvertFromFile");

  tests.run();
  tests.summary();
}
