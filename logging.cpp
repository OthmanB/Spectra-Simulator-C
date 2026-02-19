#include "logging.h"

#include <chrono>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <sstream>
#include <cctype>

#if defined(__unix__) || defined(__APPLE__)
#include <unistd.h>
#endif

namespace {
std::mutex g_log_mutex;
LogLevel g_level = LogLevel::info;
std::streambuf* g_stdout_buf = nullptr;
std::streambuf* g_stderr_buf = nullptr;
bool g_streams_redirected = false;

bool contains_word(const std::string& haystack, const std::string& needle) {
	if (needle.empty()) {
		return false;
	}
	const size_t n = needle.size();
	for (size_t pos = haystack.find(needle); pos != std::string::npos; pos = haystack.find(needle, pos + 1)) {
		const bool left_ok = (pos == 0) || !std::isalnum(static_cast<unsigned char>(haystack[pos - 1]));
		const size_t after = pos + n;
		const bool right_ok = (after >= haystack.size()) || !std::isalnum(static_cast<unsigned char>(haystack[after]));
		if (left_ok && right_ok) {
			return true;
		}
	}
	return false;
}

class LogStreamBuf : public std::stringbuf {
public:
	explicit LogStreamBuf(LogLevel level) : level_(level) {}

	int sync() override {
		std::string s = str();
		if (!s.empty()) {
			if (!s.empty() && s.back() == '\n') {
				s.pop_back();
			}
			LogLevel out_level = level_;
			std::string lower;
			lower.reserve(s.size());
			for (size_t i = 0; i < s.size(); i++) {
				lower.push_back(static_cast<char>(std::tolower(static_cast<unsigned char>(s[i]))));
			}
			if (contains_word(lower, "error") || contains_word(lower, "fatal")) {
				out_level = LogLevel::error;
			} else if (contains_word(lower, "warning") || contains_word(lower, "warn")) {
				out_level = LogLevel::warn;
			}
			log_message(out_level, "stream", 0, "", s);
			str("");
		}
		return 0;
	}

private:
	LogLevel level_;
};

LogStreamBuf* g_cout_buf = nullptr;
LogStreamBuf* g_cerr_buf = nullptr;

bool use_color() {
#if defined(__unix__) || defined(__APPLE__)
	return isatty(fileno(stderr)) != 0;
#else
	return false;
#endif
}

const char* level_name(LogLevel level) {
	switch (level) {
		case LogLevel::debug:
			return "DEBUG";
		case LogLevel::info:
			return "INFO";
		case LogLevel::warn:
			return "WARN";
		case LogLevel::error:
			return "ERROR";
	}
	return "INFO";
}

const char* level_color(LogLevel level) {
	switch (level) {
		case LogLevel::debug:
			return "\033[36m";
		case LogLevel::info:
			return "\033[32m";
		case LogLevel::warn:
			return "\033[33m";
		case LogLevel::error:
			return "\033[31m";
	}
	return "\033[0m";
}

std::string timestamp_now() {
	using namespace std::chrono;
	auto now = system_clock::now();
	auto now_time = system_clock::to_time_t(now);
	auto ms = duration_cast<milliseconds>(now.time_since_epoch()) % 1000;
	std::tm local_tm{};
#if defined(_WIN32)
	localtime_s(&local_tm, &now_time);
#else
	localtime_r(&now_time, &local_tm);
#endif
	std::ostringstream oss;
	oss << std::put_time(&local_tm, "%Y-%m-%d %H:%M:%S")
		<< '.' << std::setw(3) << std::setfill('0') << ms.count();
	return oss.str();
}

std::string basename_only(const char* path) {
	if (path == nullptr) {
		return "";
	}
	std::string s(path);
	const size_t pos = s.find_last_of("/\\");
	if (pos == std::string::npos) {
		return s;
	}
	return s.substr(pos + 1);
}
}

void init_logging(LogLevel level) {
	set_log_level(level);
	enable_std_stream_logging();
}

void set_log_level(LogLevel level) {
	g_level = level;
}

LogLevel get_log_level() {
	return g_level;
}

LogLevel log_level_from_string(const std::string& value) {
	std::string v;
	v.reserve(value.size());
	for (size_t i = 0; i < value.size(); i++) {
		v.push_back(static_cast<char>(std::tolower(static_cast<unsigned char>(value[i]))));
	}
	if (v == "debug") {
		return LogLevel::debug;
	}
	if (v == "warn" || v == "warning") {
		return LogLevel::warn;
	}
	if (v == "error") {
		return LogLevel::error;
	}
	return LogLevel::info;
}

bool try_parse_log_level(const std::string& value, LogLevel* level) {
	if (level == nullptr) {
		return false;
	}
	std::string v;
	v.reserve(value.size());
	for (size_t i = 0; i < value.size(); i++) {
		v.push_back(static_cast<char>(std::tolower(static_cast<unsigned char>(value[i]))));
	}
	if (v == "debug") {
		*level = LogLevel::debug;
		return true;
	}
	if (v == "info") {
		*level = LogLevel::info;
		return true;
	}
	if (v == "warn" || v == "warning") {
		*level = LogLevel::warn;
		return true;
	}
	if (v == "error") {
		*level = LogLevel::error;
		return true;
	}
	return false;
}

bool is_log_level_enabled(LogLevel level) {
	return static_cast<int>(level) >= static_cast<int>(g_level);
}

void log_message(LogLevel level, const char* file, int line, const char* func, const std::string& message) {
	if (!is_log_level_enabled(level)) {
		return;
	}
	std::lock_guard<std::mutex> guard(g_log_mutex);
	const bool color = use_color();
	std::ostringstream oss;
	oss << timestamp_now() << " [";
	if (color) {
		oss << level_color(level) << level_name(level) << "\033[0m";
	} else {
		oss << level_name(level);
	}
	oss << "] " << basename_only(file) << ":" << line;
	if (func != nullptr && std::string(func).size() > 0) {
		oss << " " << func;
	}
	oss << " | " << message;

	std::streambuf* target_buf = nullptr;
	if (level == LogLevel::warn || level == LogLevel::error) {
		target_buf = g_stderr_buf ? g_stderr_buf : std::cerr.rdbuf();
	} else {
		target_buf = g_stdout_buf ? g_stdout_buf : std::cout.rdbuf();
	}
	std::ostream out(target_buf);
	out << oss.str() << std::endl;
}

void enable_std_stream_logging() {
	if (g_streams_redirected) {
		return;
	}
	if (!g_stdout_buf) {
		g_stdout_buf = std::cout.rdbuf();
	}
	if (!g_stderr_buf) {
		g_stderr_buf = std::cerr.rdbuf();
	}
	if (!g_cout_buf) {
		g_cout_buf = new LogStreamBuf(LogLevel::info);
	}
	if (!g_cerr_buf) {
		g_cerr_buf = new LogStreamBuf(LogLevel::error);
	}
	std::cout.rdbuf(g_cout_buf);
	std::cerr.rdbuf(g_cerr_buf);
	g_streams_redirected = true;
}
