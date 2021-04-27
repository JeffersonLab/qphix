#include "qphix/qphix_config.h"

#include "timeDslashNoQDP.h"

#include "qmp_context.h"

int main(int argc, char **argv)
{
  QmpContext qmp_context(argc, argv);
  TimeDslash test(qmp_context.args());
  test.run();
}
