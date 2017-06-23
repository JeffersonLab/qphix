#include "qmp_context.h"
#include "testBlas.h"

int main(int argc, char **argv)
{
  QmpContext qmp_context(argc, argv);
  TestBlas test(qmp_context.args());
  test.run();
}
