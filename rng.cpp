#include "rng.h"

#include <random>

namespace {
bool g_seeded = false;
std::mt19937& rng_instance() {
	static std::mt19937 gen{std::random_device{}()};
	return gen;
}
}

void set_global_seed(uint64_t seed) {
	const uint32_t seed_lo = static_cast<uint32_t>(seed & 0xFFFFFFFFu);
	const uint32_t seed_hi = static_cast<uint32_t>((seed >> 32) & 0xFFFFFFFFu);
	std::seed_seq seq{seed_lo, seed_hi};
	rng_instance().seed(seq);
	g_seeded = true;
}

bool has_global_seed() {
	return g_seeded;
}

std::mt19937& global_rng() {
	return rng_instance();
}
