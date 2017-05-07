#pragma once

/**
  \file

  Conditionally print diagnostic messages during compilation.

  In the various USQCD projects, the `#warning` directive has been used to
  signal the various compile time options to the user. This can clutter the
  output and hide real warnings. Therefore this macro will enable toggling those
  warnings.

  With help from http://stackoverflow.com/a/43796429/653152.
  */

#include "qphix/qphix_config.h"

#ifdef QPHIX_EMIT_MESSAGES
#define QPHIX_MESSAGE_I(s) _Pragma(#s)
#define QPHIX_MESSAGE(s) QPHIX_MESSAGE_I(message(s))
#else
#define QPHIX_MESSAGE(s)
#endif
