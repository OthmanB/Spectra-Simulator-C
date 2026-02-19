#pragma once

#include <sstream>
#include <string>

enum class LogLevel {
	debug = 0,
	info = 1,
	warn = 2,
	error = 3,
};

void init_logging(LogLevel level);
void set_log_level(LogLevel level);
LogLevel get_log_level();
LogLevel log_level_from_string(const std::string& value);
bool try_parse_log_level(const std::string& value, LogLevel* level);
bool is_log_level_enabled(LogLevel level);
void log_message(LogLevel level, const char* file, int line, const char* func, const std::string& message);
void enable_std_stream_logging();

#define LOG_DEBUG(msg) do { if (is_log_level_enabled(LogLevel::debug)) { std::ostringstream _oss; _oss << msg; log_message(LogLevel::debug, __FILE__, __LINE__, __func__, _oss.str()); } } while(0)
#define LOG_INFO(msg) do { if (is_log_level_enabled(LogLevel::info)) { std::ostringstream _oss; _oss << msg; log_message(LogLevel::info, __FILE__, __LINE__, __func__, _oss.str()); } } while(0)
#define LOG_WARN(msg) do { if (is_log_level_enabled(LogLevel::warn)) { std::ostringstream _oss; _oss << msg; log_message(LogLevel::warn, __FILE__, __LINE__, __func__, _oss.str()); } } while(0)
#define LOG_ERROR(msg) do { if (is_log_level_enabled(LogLevel::error)) { std::ostringstream _oss; _oss << msg; log_message(LogLevel::error, __FILE__, __LINE__, __func__, _oss.str()); } } while(0)
