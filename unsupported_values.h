// Copyright © 2017 Martin Ueding <dev@martin-ueding.de>

#pragma once

#include <sstream>
#include <stdexcept>

#define UNSUPPORTED_SOALEN(soalen) (unsupported_soalen((soalen), __FILE__, __LINE__))

/**
  Aborts the program and raises an error due to a non-supported SOALEN value.

  Consider using the macro UNSUPPORTED_SOALEN such that the filename and line
  number are automatically passed to this function.

  @param[in] soalen The soalen in question
  @param[in] file Filename where this function was called, usually __FILE__
  @param[in] line Line number where this function was called, usually __LINE__
  */
void unsupported_soalen(int const soalen,
                        std::string const &file,
                        int const line) {
    std::ostringstream oss;
    oss << "SOALEN = " << soalen << " is not supported at " << file << ":"
        << line << ".";
    throw std::domain_error(oss.str());
}
