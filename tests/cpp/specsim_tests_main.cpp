#include "tests/cpp/tinytest.h"

#include <iostream>

int main() {
	int passed = 0;
	int failed = 0;

	for (const auto& t : tinytest::registry()) {
		try {
			t.fn();
			passed++;
		} catch (const tinytest::Failure& e) {
			failed++;
			std::cerr << "[FAIL] " << t.name << "\n" << e.what() << "\n";
		} catch (const std::exception& e) {
			failed++;
			std::cerr << "[ERROR] " << t.name << "\n" << e.what() << "\n";
		} catch (...) {
			failed++;
			std::cerr << "[ERROR] " << t.name << "\n" << "unknown exception\n";
		}
	}

	std::cerr << "Ran " << (passed + failed) << " tests: " << passed << " passed, " << failed << " failed\n";
	return failed == 0 ? 0 : 1;
}
