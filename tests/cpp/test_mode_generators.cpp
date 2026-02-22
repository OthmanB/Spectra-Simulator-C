#include "tinytest.h"

#include "io_star_params.h"
#include "models_database.h"
#include "rng.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <map>
#include <random>
#include <string>

namespace fs = std::filesystem;

static fs::path repo_root() {
	// This file lives at <repo>/tests/cpp/test_mode_generators.cpp
	fs::path p = fs::path(__FILE__).lexically_normal();
	return p.parent_path().parent_path().parent_path();
}

static fs::path make_temp_dir(const std::string& prefix) {
	fs::path base = fs::temp_directory_path();
	fs::path dir = base / fs::path(prefix + std::to_string(static_cast<unsigned long long>(std::random_device{}())));
	fs::create_directories(dir);
	return dir;
}

static void check_mode_table_basic(const Eigen::MatrixXd& m) {
	REQUIRE(m.rows() > 0);
	REQUIRE(m.cols() >= 4);
	for (int i = 0; i < m.rows(); i++) {
		const double l = m(i, 0);
		const double nu = m(i, 1);
		const double H = m(i, 2);
		const double W = m(i, 3);
		CHECK_FINITE(l);
		CHECK_FINITE(nu);
		CHECK_FINITE(H);
		CHECK_FINITE(W);
		const int li = static_cast<int>(std::llround(l));
		CHECK(std::fabs(l - li) < 1e-6);
		CHECK(li >= 0 && li <= 3);
		CHECK(nu > 0.0);
		CHECK(H >= 0.0);
		CHECK(W > 0.0);
	}

	// Monotonicity within each l (tolerant to solver rounding).
	for (int el = 0; el <= 3; el++) {
		std::vector<double> nus;
		nus.reserve(static_cast<size_t>(m.rows()));
		for (int i = 0; i < m.rows(); i++) {
			const int li = static_cast<int>(std::llround(m(i, 0)));
			if (li == el) {
				nus.push_back(m(i, 1));
			}
		}
		if (nus.size() < 2) {
			continue;
		}
		std::sort(nus.begin(), nus.end());
		for (size_t i = 0; i + 1 < nus.size(); i++) {
			CHECK(nus[i + 1] - nus[i] > 1e-8);
		}
	}
}

static Eigen::VectorXd kallinger_noise_means_from_cfg(const fs::path& noise_cfg_path) {
	const Config_Noise cfg = readNoiseConfigFile(noise_cfg_path.string());
	std::map<std::string, double> mean_by_name;
	for (size_t i = 0; i < cfg.name_random.size() && i < cfg.x1_random.size(); i++) {
		mean_by_name[cfg.name_random[i]] = cfg.x1_random[i];
	}
	const std::vector<std::string> names = {
		"k_Agran",
		"s_Agran",
		"k_taugran",
		"s_taugran",
		"c0",
		"ka",
		"ks",
		"k1",
		"s1",
		"c1",
		"k2",
		"s2",
		"c2",
		"N0",
	};
	Eigen::VectorXd v(static_cast<int>(names.size()));
	for (int i = 0; i < v.size(); i++) {
		auto it = mean_by_name.find(names[static_cast<size_t>(i)]);
		REQUIRE(it != mean_by_name.end());
		v[i] = it->second;
	}
	return v;
}

static double urand(double a, double b) {
	std::mt19937& gen = global_rng();
	std::uniform_real_distribution<double> d(a, b);
	return d(gen);
}

TEST_CASE("generate_cfg_from_synthese_file_Wscaled_aj outputs sane") {
	const fs::path root = repo_root();
	const fs::path infile = root / "Configurations" / "infiles" / "8379927.in";
	REQUIRE(fs::exists(infile));

	const fs::path tmp = make_temp_dir("specsim-aj-");
	fs::current_path(tmp);

	const fs::path modes_out = tmp / "modes_tmp.cfg";
	const fs::path noise_out = tmp / "noise_tmp.cfg";

	Eigen::VectorXd p(17);
	// Dnu, epsilon, delta0l_percent, HNR, a1ovGamma, Gamma_at_numax, a2..a6, beta_asym, i, H_spread, nu_spread, Gamma_spread, do_flat_noise
	p << 70.0, 0.5, 1.0, 10.0, 0.6, 1.0,
		0.1, -0.1, 0.15, 0.2, 0.05,
		10.0, 60.0,
		0.0, 0.0, 0.0, 0.0;

	set_global_seed(123);
	generate_cfg_from_synthese_file_Wscaled_aj(p, modes_out.string(), noise_out.string(), infile.string());
	REQUIRE(fs::exists(modes_out));
	REQUIRE(fs::exists(noise_out));

	const Data_Nd modes = read_data_ascii_Ncols(modes_out.string(), " ", false);
	check_mode_table_basic(modes.data);
	CHECK(modes.data.cols() == 16);
}

TEST_CASE("asymptotic_mm_freeDp_numaxspread_curvepmodes_v3_GRANscaled_Kallinger2014 always succeeds in safe range") {
	const fs::path root = repo_root();
	const fs::path template_file = root / "Configurations" / "templates" / "12508433.template";
	REQUIRE(fs::exists(template_file));
	const fs::path noise_cfg = root / "Configurations" / "noise_Kallinger2014.cfg";
	REQUIRE(fs::exists(noise_cfg));
	const Eigen::VectorXd noise_means = kallinger_noise_means_from_cfg(noise_cfg);

	const fs::path tmp = make_temp_dir("specsim-mm-");
	fs::current_path(tmp);

	const fs::path modes_out = tmp / "modes_tmp.cfg";
	const fs::path noise_out = tmp / "noise_tmp.cfg";
	const fs::path cfg_mm = tmp / "star_params.theoretical"; // deliberately absent
	const std::string external_path = (root / "external").string();

	set_global_seed(1000);
	for (int trial = 0; trial < 50; trial++) {
		set_global_seed(static_cast<uint64_t>(1000 + trial));

		Eigen::VectorXd p(45);
		p.setZero();
		// Use the example-cfg parameter envelope (Configurations/examples_cfg/main.cfg.freeDP_curvepmodes.v3_GRANscaled)
		// to keep the property sweep within a domain known to be stable in practice.
		p[0] = urand(0.1286, 0.3858); // nurot_env
		p[1] = urand(0.1286, 3.858);  // nurot_core
		p[2] = 0.0;               // a2_l1_core (unused)
		p[3] = 0.0;               // a2_l1_env (unused)
		p[4] = urand(-0.1, 0.1);    // a2_l2_env
		p[5] = urand(-0.1, 0.1);    // a2_l3_env
		p[6] = urand(-0.05, 0.05);  // a3_l2_env
		p[7] = urand(-0.05, 0.05);  // a3_l3_env
		p[8] = urand(-0.1, 0.1);    // a4_l2_env
		p[9] = urand(-0.1, 0.1);    // a4_l3_env
		p[10] = urand(-0.05, 0.05); // a5_l3_env
		p[11] = urand(-0.05, 0.05); // a6_l3_env
		p[12] = urand(20.0, 30.0);  // Dnu
		p[13] = urand(0.0, 1.0);    // epsilon
		p[14] = urand(0.0, 0.5);    // delta0l_percent (keep non-negative)
		p[15] = urand(0.0, 0.10);   // beta_p_star
		p[16] = 5.0;                // nmax_spread
		p[17] = urand(75.0, 400.0); // DP1
		p[18] = urand(0.0, 1.0);    // alpha
		p[19] = urand(0.1, 0.5);    // q
		p[20] = urand(1000.0, 5000.0); // SNR
		p[21] = urand(0.1, 1.5);    // maxGamma
		p[22] = 15.0;               // numax_spread (percent)
		p[23] = urand(1.0, 1.75);   // Vl1
		p[24] = urand(0.25, 0.80);  // Vl2
		p[25] = urand(0.0, 0.10);   // Vl3
		p[26] = 10.0;               // H0_spread
		// Noise parameters from Kallinger2014 cfg (use means for stability).
		p.segment(27, 14) = noise_means;
		p[41] = urand(0.0, 1.0);  // Hfactor
		p[42] = urand(0.0, 0.9);  // Wfactor must stay <1 to avoid negative mixed-mode widths (gamma_l_fct2)
		p[43] = 100.0;           // Tobs (days)
		p[44] = 120.0;           // Cadence (s)

		const bool failed = asymptotic_mm_freeDp_numaxspread_curvepmodes_v3_GRANscaled_Kallinger2014(
			p,
			modes_out.string(),
			noise_out.string(),
			cfg_mm.string(),
			external_path,
			template_file.string());
		REQUIRE(!failed);
		REQUIRE(fs::exists(modes_out));
		REQUIRE(fs::exists(noise_out));

		const Data_Nd modes = read_data_ascii_Ncols(modes_out.string(), " ", false);
		check_mode_table_basic(modes.data);
		CHECK(modes.data.cols() >= 12);
	}
}
