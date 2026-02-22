#include "tinytest.h"

#include "function_rot.h"

#include <cmath>

TEST_CASE("amplitude_ratio basics") {
	for (int l = 0; l <= 3; l++) {
		for (double beta : {0.0, 30.0, 60.0, 90.0}) {
			const Eigen::VectorXd v = amplitude_ratio(l, beta);
			REQUIRE(v.size() == 2 * l + 1);
			double sum = 0.0;
			for (int i = 0; i < v.size(); i++) {
				CHECK_FINITE(v[i]);
				CHECK(v[i] >= 0.0);
				sum += v[i];
			}
			// Function_rot is used as a set of m-component visibility weights; it should be normalized.
			CHECK_NEAR(sum, 1.0, 1e-12, 1e-12);
		}
	}
}
