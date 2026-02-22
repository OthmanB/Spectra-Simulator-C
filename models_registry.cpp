#include "models_registry.h"

#include <filesystem>
#include <sstream>

#include "artificial_spectrum.h"
#include "logging.h"
#include "models_database.h"

namespace {
void append_tobs_cadence(Eigen::VectorXd& input_params, const Config_Data& cfg) {
	const size_t old_size = static_cast<size_t>(input_params.size());
	input_params.conservativeResize(static_cast<int>(old_size + 2));
	input_params[static_cast<int>(old_size)] = cfg.Tobs;
	input_params[static_cast<int>(old_size + 1)] = cfg.Cadence;
}

void copy_mixed_mode_cfg(const ModelRuntime& ctx, const std::string& id_str) {
	if (ctx.file_cfg_mm.empty()) {
		return;
	}
	std::filesystem::path src(ctx.file_cfg_mm);
	std::filesystem::path dst = std::filesystem::path(ctx.data_path) / "Spectra_info" / (strtrim(id_str) + ".global");
	std::error_code ec;
	std::filesystem::copy_file(src, dst, std::filesystem::copy_options::overwrite_existing, ec);
	if (ec) {
		LOG_ERROR("Warning: failed to copy mixed-mode cfg file: " << src.string());
		LOG_ERROR("         to: " << dst.string());
		LOG_ERROR("         reason: " << ec.message());
	}
}

bool run_alm_random(Eigen::VectorXd input_params, const ModelRuntime& ctx, const std::string& id_str) {
	const Config_Data& cfg = *ctx.cfg;
	generate_cfg_from_synthese_file_Wscaled_Alm(input_params, ctx.file_out_modes, ctx.file_out_noise, cfg.extra_params);
	artificial_spectrum_a1Alma3(cfg.Tobs, cfg.Cadence, cfg.Nspectra, cfg.Nrealisation, ctx.dir_core, id_str, cfg.doplots,
			cfg.write_inmodel, cfg.do_modelfiles, cfg.limit_data_range, cfg.modefile_modelname, ctx.data_path);
	return true;
}

bool run_alm_grid(Eigen::VectorXd input_params, const ModelRuntime& ctx, const std::string& id_str, const Model_data&) {
	return run_alm_random(input_params, ctx, id_str);
}

bool run_aj_random(Eigen::VectorXd input_params, const ModelRuntime& ctx, const std::string& id_str) {
	const Config_Data& cfg = *ctx.cfg;
	generate_cfg_from_synthese_file_Wscaled_aj(input_params, ctx.file_out_modes, ctx.file_out_noise, cfg.extra_params);
	artificial_spectrum_aj(cfg.Tobs, cfg.Cadence, cfg.Nspectra, cfg.Nrealisation, ctx.dir_core, id_str, cfg.doplots,
			cfg.write_inmodel, cfg.do_modelfiles, cfg.limit_data_range, cfg.modefile_modelname, ctx.data_path, "harvey_like");
	return true;
}

bool run_aj_grid(Eigen::VectorXd input_params, const ModelRuntime& ctx, const std::string& id_str, const Model_data&) {
	const Config_Data& cfg = *ctx.cfg;
	generate_cfg_from_synthese_file_Wscaled_aj(input_params, ctx.file_out_modes, ctx.file_out_noise, cfg.extra_params);
	artificial_spectrum_aj(cfg.Tobs, cfg.Cadence, cfg.Nspectra, cfg.Nrealisation, ctx.dir_core, id_str, cfg.doplots,
			cfg.write_inmodel, cfg.do_modelfiles, cfg.limit_data_range, cfg.modefile_modelname, ctx.data_path);
	return true;
}

bool run_aj_kallinger_random(Eigen::VectorXd input_params, const ModelRuntime& ctx, const std::string& id_str) {
	const Config_Data& cfg = *ctx.cfg;
	append_tobs_cadence(input_params, cfg);
	generate_cfg_from_synthese_file_Wscaled_aj_GRANscaled_Kallinger2014(
			input_params,
			ctx.file_out_modes,
			ctx.file_out_noise,
			cfg.extra_params);
	artificial_spectrum_aj(cfg.Tobs, cfg.Cadence, cfg.Nspectra, cfg.Nrealisation, ctx.dir_core, id_str, cfg.doplots,
			cfg.write_inmodel, cfg.do_modelfiles, cfg.limit_data_range, cfg.modefile_modelname, ctx.data_path, "harvey_like");
	return true;
}

bool run_mm_v3_kallinger_random(Eigen::VectorXd input_params, const ModelRuntime& ctx, const std::string& id_str) {
	const Config_Data& cfg = *ctx.cfg;
	append_tobs_cadence(input_params, cfg);
	const bool failed = asymptotic_mm_freeDp_numaxspread_curvepmodes_v3_GRANscaled_Kallinger2014(
			input_params,
			ctx.file_out_modes,
			ctx.file_out_noise,
			ctx.file_cfg_mm,
			ctx.external_path,
			ctx.template_file);
	if (!failed) {
		artificial_spectrum_aj(cfg.Tobs, cfg.Cadence, cfg.Nspectra, cfg.Nrealisation, ctx.dir_core, id_str, cfg.doplots,
				cfg.write_inmodel, cfg.do_modelfiles, cfg.limit_data_range, cfg.modefile_modelname, ctx.data_path, "harvey_like");
	}
	copy_mixed_mode_cfg(ctx, id_str);
	return true;
}

const std::vector<ModelSpec>& model_specs() {
	static const std::vector<ModelSpec> specs = {
		{
			"generate_cfg_from_synthese_file_Wscaled_Alm",
			ForestSupport::random_and_grid,
			{
				"HNR",
				"a1ovGamma",
				"Gamma_at_numax",
				"epsilon_nl",
				"theta0",
				"delta",
				"a3",
				"beta_asym",
				"i",
			},
			true,
			false,
			false,
			false,
			run_alm_random,
			run_alm_grid,
		},
		{
			"generate_cfg_from_synthese_file_Wscaled_aj",
			ForestSupport::random_and_grid,
			{
				"Dnu",
				"epsilon",
				"delta0l_percent",
				"HNR",
				"a1ovGamma",
				"Gamma_at_numax",
				"a2",
				"a3",
				"a4",
				"a5",
				"a6",
				"beta_asym",
				"i",
				"H_spread",
				"nu_spread",
				"Gamma_spread",
				"do_flat_noise",
			},
			true,
			false,
			false,
			false,
			run_aj_random,
			run_aj_grid,
		},
		{
			"generate_cfg_from_synthese_file_Wscaled_aj_GRANscaled_Kallinger2014",
			ForestSupport::random_only,
			{
				"Dnu",
				"epsilon",
				"delta0l_percent",
				"HNR",
				"a1ovGamma",
				"Gamma_at_numax",
				"a2",
				"a3",
				"a4",
				"a5",
				"a6",
				"beta_asym",
				"i",
				"numax_spread",
				"H_spread",
				"nu_spread",
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
			},
			true,
			false,
			true,
			true,
			run_aj_kallinger_random,
			nullptr,
		},
		{
			"asymptotic_mm_freeDp_numaxspread_curvepmodes_v3_GRANscaled_Kallinger2014",
			ForestSupport::random_only,
			{
				"nurot_env",
				"nurot_core",
				"a2_l1_core",
				"a2_l1_env",
				"a2_l2_env",
				"a2_l3_env",
				"a3_l2_env",
				"a3_l3_env",
				"a4_l2_env",
				"a4_l3_env",
				"a5_l3_env",
				"a6_l3_env",
				"Dnu",
				"epsilon",
				"delta0l_percent",
				"beta_p_star",
				"nmax_spread",
				"DP1",
				"alpha",
				"q",
				"SNR",
				"maxGamma",
				"numax_spread",
				"Vl1",
				"Vl2",
				"Vl3",
				"H0_spread",
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
				"Hfactor",
				"Wfactor",
			},
			false,
			true,
			true,
			true,
			run_mm_v3_kallinger_random,
			nullptr,
		},
	};
	return specs;
}
}

const ModelSpec* find_model_spec(const std::string& name) {
	const auto& specs = model_specs();
	for (const auto& spec : specs) {
		if (spec.name == name) {
			return &spec;
		}
	}
	return nullptr;
}

std::vector<const ModelSpec*> list_model_specs() {
	const auto& specs = model_specs();
	std::vector<const ModelSpec*> out;
	out.reserve(specs.size());
	for (const auto& spec : specs) {
		out.push_back(&spec);
	}
	return out;
}

std::string forest_support_to_string(ForestSupport support) {
	switch (support) {
		case ForestSupport::random_only:
			return "random";
		case ForestSupport::grid_only:
			return "grid";
		case ForestSupport::random_and_grid:
			return "random,grid";
	}
	return "unknown";
}

bool model_supports_random(const ModelSpec& spec) {
	return spec.forest_support == ForestSupport::random_only || spec.forest_support == ForestSupport::random_and_grid;
}

bool model_supports_grid(const ModelSpec& spec) {
	return spec.forest_support == ForestSupport::grid_only || spec.forest_support == ForestSupport::random_and_grid;
}

std::string describe_model(const ModelSpec& spec) {
	std::ostringstream oss;
	oss << "Model: " << spec.name << "\n";
	oss << "Forest support: " << forest_support_to_string(spec.forest_support) << "\n";
	oss << "Needs extra_params: " << (spec.needs_extra_params ? "yes" : "no") << "\n";
	oss << "Needs templates: " << (spec.needs_templates ? "yes" : "no") << "\n";
	oss << "Needs noise cfg: " << (spec.needs_noise_cfg ? "yes" : "no") << "\n";
	oss << "Appends Tobs/Cadence: " << (spec.append_tobs_cadence ? "yes" : "no") << "\n";
	oss << "Parameter order (" << spec.param_names.size() << "):";
	for (const auto& name : spec.param_names) {
		oss << "\n  " << name;
	}
	return oss.str();
}
