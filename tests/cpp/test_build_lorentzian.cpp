#include "tinytest.h"

#include "build_lorentzian.h"
#include "function_rot.h"
#include "rng.h"

#include <random>

static double urand(double a, double b) {
	std::mt19937& gen = global_rng();
	std::uniform_real_distribution<double> d(a, b);
	return d(gen);
}

TEST_CASE("build_l_mode_aj non-negative finite") {
	set_global_seed(123);

	for (int trial = 0; trial < 50; trial++) {
		const int l = static_cast<int>(urand(0.0, 4.0));
		const double fc = urand(10.0, 4000.0);
		const double gamma = urand(0.01, 5.0);
		const double H = urand(1e-6, 1e3);
		const double asym = urand(-0.2, 0.2);
		const double a1 = urand(-1.0, 1.0);
		const double a2 = urand(-0.2, 0.2);
		const double a3 = urand(-0.2, 0.2);
		const double a4 = urand(-0.2, 0.2);
		const double a5 = urand(-0.2, 0.2);
		const double a6 = urand(-0.2, 0.2);
		const double eta0 = 0.0;
		const Eigen::VectorXd V = amplitude_ratio(l, urand(0.0, 90.0));

		Eigen::VectorXd x(4096);
		const double span = 50.0;
		x.setLinSpaced(x.size(), fc - span, fc + span);
		const Eigen::VectorXd y = build_l_mode_aj(x, H, fc, a1, a2, a3, a4, a5, a6, eta0, asym, gamma, l, V);

		REQUIRE(y.size() == x.size());
		for (int i = 0; i < y.size(); i++) {
			CHECK_FINITE(y[i]);
			CHECK(y[i] >= 0.0);
		}
	}
}

TEST_CASE("build_l_mode_aj symmetry when no splitting and no asymmetry") {
	const int l = 1;
	const double fc = 1000.0;
	const double gamma = 1.5;
	const double H = 10.0;
	const double a1 = 0.0, a2 = 0.0, a3 = 0.0, a4 = 0.0, a5 = 0.0, a6 = 0.0;
	const double eta0 = 0.0;
	const double asym = 0.0;
	const Eigen::VectorXd V = amplitude_ratio(l, 45.0);

	Eigen::VectorXd x(2001);
	x.setLinSpaced(x.size(), fc - 40.0, fc + 40.0);
	const Eigen::VectorXd y = build_l_mode_aj(x, H, fc, a1, a2, a3, a4, a5, a6, eta0, asym, gamma, l, V);

	for (int i = 0; i < y.size(); i++) {
		CHECK_FINITE(y[i]);
		CHECK(y[i] >= 0.0);
		const int j = y.size() - 1 - i;
		CHECK_NEAR(y[i], y[j], 1e-12, 1e-10);
	}
}
