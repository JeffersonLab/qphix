#include "timeClovNoQDP.h"
#include "qmp_context.h"

#include "qphix/qphix_config.h"

int main(int argc, char **argv)
{
  QmpContext qmp_context(argc, argv);
  TimeClover test(qmp_context.args());
  test.run();
}
