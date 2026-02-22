#pragma once

#include <string>
#include <vector>
#include <Eigen/Dense>

struct Config_Data;
struct Model_data;

enum class ForestSupport {
	random_only,
	grid_only,
	random_and_grid,
};

struct ModelRuntime {
	std::string file_out_modes;
	std::string file_out_noise;
	std::string file_cfg_mm;
	std::string dir_core;
	std::string data_path;
	std::string external_path;
	std::string template_file;
	const Config_Data* cfg;
};

using ModelRunRandom = bool (*)(Eigen::VectorXd input_params, const ModelRuntime& ctx, const std::string& id_str);
using ModelRunGrid = bool (*)(Eigen::VectorXd input_params, const ModelRuntime& ctx, const std::string& id_str, const Model_data& model);

struct ModelSpec {
	std::string name;
	ForestSupport forest_support;
	std::vector<std::string> param_names;
	bool needs_extra_params;
	bool needs_templates;
	bool needs_noise_cfg;
	bool append_tobs_cadence;
	ModelRunRandom run_random;
	ModelRunGrid run_grid;
};

const ModelSpec* find_model_spec(const std::string& name);
std::vector<const ModelSpec*> list_model_specs();
std::string describe_model(const ModelSpec& spec);
std::string forest_support_to_string(ForestSupport support);
bool model_supports_random(const ModelSpec& spec);
bool model_supports_grid(const ModelSpec& spec);
