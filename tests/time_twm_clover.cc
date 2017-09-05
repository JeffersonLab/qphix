#include "qphix/qphix_config.h"

#include "timeTWMClover.h"

#include "qmp_context.h"

int main(int argc, char **argv)
{
  QmpContext qmp_context(argc, argv);
  TimeTMClover test(qmp_context.args());
  test.run();
}
