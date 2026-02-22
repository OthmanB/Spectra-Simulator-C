#pragma once

#include <cstdint>
#include <random>

void set_global_seed(uint64_t seed);
bool has_global_seed();
std::mt19937& global_rng();
