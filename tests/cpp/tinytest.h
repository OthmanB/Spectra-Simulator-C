#pragma once

#include <cmath>
#include <cstddef>
#include <exception>
#include <functional>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace tinytest {

struct Test {
	const char* name;
	std::function<void()> fn;
};

inline std::vector<Test>& registry() {
	static std::vector<Test> tests;
	return tests;
}

inline void register_test(const char* name, std::function<void()> fn) {
	registry().push_back(Test{name, std::move(fn)});
}

struct Failure : public std::exception {
	std::string msg;
	explicit Failure(std::string m) : msg(std::move(m)) {}
	const char* what() const noexcept override { return msg.c_str(); }
};

inline std::string loc_str(const char* file, int line) {
	std::ostringstream oss;
	oss << file << ":" << line;
	return oss.str();
}

inline void fail(const char* file, int line, const std::string& msg) {
	throw Failure(loc_str(file, line) + ": " + msg);
}

inline bool is_finite(double v) {
	return std::isfinite(v);
}

inline bool near(double a, double b, double atol, double rtol) {
	const double diff = std::fabs(a - b);
	const double scale = std::max(std::fabs(a), std::fabs(b));
	return diff <= (atol + rtol * scale);
}

} // namespace tinytest

// Macro helpers
#define TINYTEST_CONCAT_INNER(a, b) a##b
#define TINYTEST_CONCAT(a, b) TINYTEST_CONCAT_INNER(a, b)

#define TEST_CASE(name_literal) \
	static void TINYTEST_CONCAT(tinytest_case_, __LINE__)(); \
	namespace { \
	struct TINYTEST_CONCAT(TinyTestReg_, __LINE__) { \
		TINYTEST_CONCAT(TinyTestReg_, __LINE__)() { \
			tinytest::register_test(name_literal, &TINYTEST_CONCAT(tinytest_case_, __LINE__)); \
		} \
	}; \
	static TINYTEST_CONCAT(TinyTestReg_, __LINE__) TINYTEST_CONCAT(tinytest_reg_instance_, __LINE__); \
	} \
	static void TINYTEST_CONCAT(tinytest_case_, __LINE__)()

#define REQUIRE(expr) \
	do { \
		if (!(expr)) { \
			tinytest::fail(__FILE__, __LINE__, std::string("REQUIRE failed: ") + #expr); \
		} \
	} while (0)

#define CHECK(expr) \
	do { \
		if (!(expr)) { \
			tinytest::fail(__FILE__, __LINE__, std::string("CHECK failed: ") + #expr); \
		} \
	} while (0)

#define CHECK_FINITE(x) \
	do { \
		double _v = static_cast<double>(x); \
		if (!tinytest::is_finite(_v)) { \
			tinytest::fail(__FILE__, __LINE__, std::string("not finite: ") + #x); \
		} \
	} while (0)

#define CHECK_NEAR(a, b, atol, rtol) \
	do { \
		double _a = static_cast<double>(a); \
		double _b = static_cast<double>(b); \
		double _at = static_cast<double>(atol); \
		double _rt = static_cast<double>(rtol); \
		if (!tinytest::near(_a, _b, _at, _rt)) { \
			std::ostringstream _oss; \
			_oss << "CHECK_NEAR failed: " << #a << "=" << _a << ", " << #b << "=" << _b \
				 << " (atol=" << _at << ", rtol=" << _rt << ")"; \
			tinytest::fail(__FILE__, __LINE__, _oss.str()); \
		} \
	} while (0)
