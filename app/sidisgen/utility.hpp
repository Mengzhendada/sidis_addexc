#ifndef SIDISGEN_UTILITY_HPP
#define SIDISGEN_UTILITY_HPP

#include <ostream>
#include <type_traits>

#include <iostream>

// Draws a progress bar in the terminal.
bool write_progress_bar(std::ostream& os, unsigned percent, unsigned width=70);

#endif

