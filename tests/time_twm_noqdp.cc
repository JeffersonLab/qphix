#include "qphix/qphix_config.h"

#include "timeTWMNoQDP.h"

#include "qmp_context.h"

int main(int argc, char **argv)
{
  QmpContext qmp_context(argc, argv);
  auto args = processArgs(argc, argv, true);

  TimeTMDslash test(qmp_context.args());
  test.run();
}
