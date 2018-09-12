// logic.cpp : tests for logic operators between posits
//
// Copyright (C) 2017-2018 Stillwater Supercomputing, Inc.
//
// This file is part of the universal numbers project, which is released under an MIT Open Source license.

#include "common.hpp"

// minimum set of include files to reflect source code dependencies
#define POSIT_ENABLE_LITERALS 1
#include "../../posit/posit.hpp"
#include "../tests/test_helpers.hpp"

// Posit equal diverges from IEEE float in dealing with INFINITY/NAN
// Posit NaR can be checked for equality/inequality
template<size_t nbits, size_t es>
int ValidatePositLogicEqual() {
	const size_t NR_TEST_CASES = (unsigned(1) << nbits);
	int nrOfFailedTestCases = 0;
	sw::unum::posit<nbits, es> a, b;
	bool ref, presult;

	for (unsigned i = 0; i < NR_TEST_CASES; i++) {
		a.set_raw_bits(i);
		for (unsigned j = 0; j < NR_TEST_CASES; j++) {
			b.set_raw_bits(j);
			// set the golden reference
			if (a.isNaR() && b.isNaR()) {
				// special case of posit equality
				ref = true;
			}
			else {
				// initially, we thought this would be the same behavior as IEEE floats
				// ref = double(a) == double(b);
				// but we have found that some compilers (MSVC) take liberty with NaN
				// \fp:fast		floating point model set to fast
				//	NaN == NaN  : IEEE = true    Posit = true
				//	NaN == real : IEEE = true    Posit = false
				// \fp:strict	floating point model set to strict
				//	NaN == NaN  : IEEE = false    Posit = true
				//	NaN == real : IEEE = false    Posit = false
				// and thus we can't relay on IEEE float as reference

				// instead, use the bit pattern as reference
				ref = (i == j ? true : false);
			}

			presult = a == b;
			if (ref != presult) {
				nrOfFailedTestCases++;
				std::cout << a << " == " << b << " fails: reference is " << ref << " actual is " << presult << std::endl;
			}
		}
	}
	return nrOfFailedTestCases;
}

// Posit not-equal diverges from IEEE float in dealing with INFINITY/NAN
// Posit NaR can be checked for equality/inequality
template<size_t nbits, size_t es>
int ValidatePositLogicNotEqual() {
	const size_t NR_TEST_CASES = (unsigned(1) << nbits);
	int nrOfFailedTestCases = 0;
	sw::unum::posit<nbits, es> a, b;
	bool ref, presult;

	for (unsigned i = 0; i < NR_TEST_CASES; i++) {
		a.set_raw_bits(i);
		for (unsigned j = 0; j < NR_TEST_CASES; j++) {
			b.set_raw_bits(j);

			// set the golden reference
			if (a.isNaR() && b.isNaR()) {
				// special case of posit equality
				ref = false;
			}
			else {
				// initially, we thought this would be the same behavior as IEEE floats
				// ref = double(a) == double(b);
				// but we have found that some compilers (MSVC) take liberty with NaN
				// \fp:fast		floating point model set to fast
				//	NaN == NaN  : IEEE = true    Posit = true
				//	NaN == real : IEEE = true    Posit = false
				// \fp:strict	floating point model set to strict
				//	NaN == NaN  : IEEE = false    Posit = true
				//	NaN == real : IEEE = false    Posit = false
				// and thus we can't relay on IEEE float as reference

				// instead, use the bit pattern as reference
				ref = (i != j ? true : false);
			}

			presult = a != b;

			if (ref != presult) {
				nrOfFailedTestCases++;
				std::cout << a << " != " << b << " fails: reference is " << ref << " actual is " << presult << std::endl;
			}
		}
	}
	return nrOfFailedTestCases;
}

// Posit less-than diverges from IEEE float in dealing with INFINITY/NAN
// Posit NaR is smaller than any other value
template<size_t nbits, size_t es>
int ValidatePositLogicLessThan() {
	const size_t NR_TEST_CASES = (unsigned(1) << nbits);
	int nrOfFailedTestCases = 0;
	sw::unum::posit<nbits, es> a, b;
	bool ref, presult;

	for (unsigned i = 0; i < NR_TEST_CASES; i++) {
		a.set_raw_bits(i);
		for (unsigned j = 0; j < NR_TEST_CASES; j++) {
			b.set_raw_bits(j);

			// generate the golden reference
			if (a.isNaR() && !b.isNaR()) {
				// special case of posit NaR
				ref = true;
			}
			else {
				// same behavior as IEEE floats
				ref = double(a) < double(b);
			}

			presult = a < b;
			if (ref != presult) {
				nrOfFailedTestCases++;
				std::cout << a << " < " << b << " fails: reference is " << ref << " actual is " << presult << std::endl;
			}
		}
	}
	return nrOfFailedTestCases;
}

// Posit greater-than diverges from IEEE float in dealing with INFINITY/NAN
// Any number is greater-than posit NaR
template<size_t nbits, size_t es>
int ValidatePositLogicGreaterThan() {
	const size_t NR_TEST_CASES = (unsigned(1) << nbits);
	int nrOfFailedTestCases = 0;
	sw::unum::posit<nbits, es> a, b;
	bool ref, presult;

	for (unsigned i = 0; i < NR_TEST_CASES; i++) {
		a.set_raw_bits(i);
		for (unsigned j = 0; j < NR_TEST_CASES; j++) {
			b.set_raw_bits(j);

			// generate the golden reference
			if (!a.isNaR() && b.isNaR()) {
				// special case of posit NaR
				ref = true;
			}
			else {
				// same behavior as IEEE floats
				ref = double(a) > double(b);
			}

			presult = a > b;
			if (ref != presult) {
				nrOfFailedTestCases++;
				std::cout << a << " > " << b << " fails: reference is " << ref << " actual is " << presult << std::endl;
			}
		}
	}
	return nrOfFailedTestCases;
}

// Posit less-or-equal-than diverges from IEEE float in dealing with INFINITY/NAN
// Posit NaR is smaller or equal than any other value
template<size_t nbits, size_t es>
int ValidatePositLogicLessOrEqualThan() {
	const size_t NR_TEST_CASES = (unsigned(1) << nbits);
	int nrOfFailedTestCases = 0;
	sw::unum::posit<nbits, es> a, b;
	bool ref, presult;

	for (unsigned i = 0; i < NR_TEST_CASES; i++) {
		a.set_raw_bits(i);
		for (unsigned j = 0; j < NR_TEST_CASES; j++) {
			b.set_raw_bits(j);

			// set the golden reference
			if (a.isNaR()) {
				// special case of posit <= for NaR
				ref = true;
			}
			else {
				// same behavior as IEEE floats
				ref = double(a) <= double(b);
			}

			presult = a <= b;

			if (ref != presult) {
				nrOfFailedTestCases++;
				std::cout << a << " <= " << b << " fails: reference is " << ref << " actual is " << presult << std::endl;
			}
		}
	}
	return nrOfFailedTestCases;
}

// Posit greater-or-equal-than diverges from IEEE float in dealing with INFINITY/NAN
// Any number is greater-or-equal-than posit NaR
template<size_t nbits, size_t es>
int ValidatePositLogicGreaterOrEqualThan() {
	const size_t NR_TEST_CASES = (unsigned(1) << nbits);
	int nrOfFailedTestCases = 0;
	sw::unum::posit<nbits, es> a, b;
	bool ref, presult;

	for (unsigned i = 0; i < NR_TEST_CASES; i++) {
		a.set_raw_bits(i);
		for (unsigned j = 0; j < NR_TEST_CASES; j++) {
			b.set_raw_bits(j);

			// set the golden reference
			if (b.isNaR()) {
				// special case of posit >= for NaR
				ref = true;
			}
			else {
				// same behavior as IEEE floats
				ref = double(a) >= double(b);
			}

			presult = a >= b;

			if (ref != presult) {
				nrOfFailedTestCases++;
				std::cout << a << " >= " << b << " fails: reference is " << ref << " actual is " << presult << std::endl;
			}
		}
	}
	return nrOfFailedTestCases;
}

#define MANUAL_TESTING 0
#define STRESS_TESTING 0

int main()
try {
	using namespace std;
	using namespace sw::unum;

	int nrOfFailedTestCases = 0;

#if MANUAL_TESTING
	constexpr size_t nbits = 8;
	constexpr size_t es    = 0;

	double nan    = NAN;
	double inf    = INFINITY;
	double normal = 0;
	posit<nbits, es> pa(nan), pb(inf), pc(normal);
	std::cout << pa << " " << pb << " " << pc << std::endl;
	
	// showcasing the differences between posit and IEEE float
	std::cout << "NaN ==  NaN: IEEE=" << (nan == nan ? "true" : "false")    << "    Posit=" << (pa == pa ? "true" : "false") << std::endl;
	std::cout << "NaN == real: IEEE=" << (nan == normal ? "true" : "false") << "    Posit=" << (pa == pc ? "true" : "false") << std::endl;
	std::cout << "INF ==  INF: IEEE=" << (inf == inf ? "true" : "false")    << "    Posit=" << (pb == pb ? "true" : "false") << std::endl;
	std::cout << "NaN !=  NaN: IEEE=" << (nan != nan ? "true" : "false")    << "   Posit=" << (pa != pb ? "true" : "false") << std::endl;
	std::cout << "INF !=  INF: IEEE=" << (inf != inf ? "true" : "false")    << "   Posit=" << (pb != pb ? "true" : "false") << std::endl;
	std::cout << "NaN <= real: IEEE=" << (nan <= normal ? "true" : "false") << "    Posit=" << (pa <= pc ? "true" : "false") << std::endl;
	std::cout << "NaN >= real: IEEE=" << (nan >= normal ? "true" : "false") << "    Posit=" << (pa >= pc ? "true" : "false") << std::endl;
	std::cout << "INF <  real: IEEE=" << (inf <  normal ? "true" : "false") << "   Posit=" << (pa <  pc ? "true" : "false") << std::endl;
	std::cout << "INF >  real: IEEE=" << (inf  > normal ? "true" : "false") << "    Posit=" << (pa  > pc ? "true" : "false") << std::endl;

	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicEqual<3, 0>(), "posit<3,0>", "==");
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicNotEqual<3, 0>(), "posit<3,0>", "!=");
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicLessThan<3, 0>(), "posit<3,0>", "<");
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicGreaterThan<3, 0>(), "posit<3,0>", ">");
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicLessOrEqualThan<3, 0>(), "posit<3,0>", "<=");
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicGreaterOrEqualThan<3, 0>(), "posit<3,0>", ">=");

#else
	posit<16, 1> p;

	cout << "Logic: operator==()" << endl;
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicEqual<3, 0>(), "posit<3,0>", "==");
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicEqual<4, 0>(), "posit<4,0>", "==");
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicEqual<4, 1>(), "posit<4,1>", "==");
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicEqual<5, 0>(), "posit<5,0>", "==");
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicEqual<5, 1>(), "posit<5,1>", "==");
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicEqual<5, 2>(), "posit<5,2>", "==");
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicEqual<6, 0>(), "posit<6,0>", "==");
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicEqual<6, 1>(), "posit<6,1>", "==");
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicEqual<6, 2>(), "posit<6,2>", "==");
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicEqual<6, 3>(), "posit<6,3>", "==");
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicEqual<7, 0>(), "posit<7,0>", "==");
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicEqual<7, 1>(), "posit<7,1>", "==");
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicEqual<7, 2>(), "posit<7,2>", "==");
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicEqual<7, 3>(), "posit<7,3>", "==");
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicEqual<7, 4>(), "posit<7,4>", "==");
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicEqual<8, 0>(), "posit<8,0>", "==");
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicEqual<8, 1>(), "posit<8,1>", "==");
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicEqual<8, 2>(), "posit<8,2>", "==");
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicEqual<8, 3>(), "posit<8,3>", "==");
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicEqual<8, 4>(), "posit<8,4>", "==");
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicEqual<8, 5>(), "posit<8,5>", "==");
	if (!(p == 0)) { 
		nrOfFailedTestCases += ReportTestResult(1, "posit<16,1> == 0", "== int literal");
	}
	else {
		ReportTestResult(0, "posit<16,1> == 0", "== int literal");
	}
	if (!(p == 0.0f)) {
		nrOfFailedTestCases += ReportTestResult(1, "posit<16,1> == 0.0", "== float literal");
	}
	else {
		ReportTestResult(0, "posit<16,1> == 0.0", "== float literal");
	}
	if (!(p == 0.0)) { 
		nrOfFailedTestCases += ReportTestResult(1, "posit<16,1> == 0.0", "== double literal");
	}
	else {
		ReportTestResult(0, "posit<16,1> == 0.0", "== double literal");
	}
	if (!(p == 0.0l)) {
		nrOfFailedTestCases += ReportTestResult(1, "posit<16,1> == 0.0", "== long double literal");
	}
	else {
		ReportTestResult(0, "posit<16,1> == 0.0", "== long double literal");
	}
	
	cout << "Logic: operator!=()" << endl;
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicNotEqual<3, 0>(), "posit<3,0>", "!=");
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicNotEqual<4, 0>(), "posit<4,0>", "!=");
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicNotEqual<4, 1>(), "posit<4,1>", "!=");
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicNotEqual<5, 0>(), "posit<5,0>", "!=");
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicNotEqual<5, 1>(), "posit<5,1>", "!=");
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicNotEqual<5, 2>(), "posit<5,2>", "!=");
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicNotEqual<6, 0>(), "posit<6,0>", "!=");
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicNotEqual<6, 1>(), "posit<6,1>", "!=");
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicNotEqual<6, 2>(), "posit<6,2>", "!=");
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicNotEqual<6, 3>(), "posit<6,3>", "!=");
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicNotEqual<7, 0>(), "posit<7,0>", "!=");
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicNotEqual<7, 1>(), "posit<7,1>", "!=");
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicNotEqual<7, 2>(), "posit<7,2>", "!=");
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicNotEqual<7, 3>(), "posit<7,3>", "!=");
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicNotEqual<8, 0>(), "posit<8,0>", "!=");
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicNotEqual<8, 1>(), "posit<8,1>", "!=");
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicNotEqual<8, 2>(), "posit<8,2>", "!=");
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicNotEqual<8, 3>(), "posit<8,3>", "!=");
	if (p != 0) { 
		nrOfFailedTestCases += ReportTestResult(1, "posit<16,1> != 0", "!= int literal");
	}
	else {
		ReportTestResult(0, "posit<16,1> != 0", "!= int literal");
	}
	if (p != 0.0f) {
		nrOfFailedTestCases += ReportTestResult(1, "posit<16,1> != 0.0", "!= float literal");
	}
	else {
		ReportTestResult(0, "posit<16,1> != 0.0", "!= float literal");
	}
	if (p != 0.0) {
		nrOfFailedTestCases += ReportTestResult(1, "posit<16,1> != 0.0", "!= double literal");
	}
	else {
		ReportTestResult(0, "posit<16,1> != 0.0", "!= double literal");
	}
	if (p != 0.0l) {
		nrOfFailedTestCases += ReportTestResult(1, "posit<16,1> != 0.0", "!= long double literal");
	}
	else {
		ReportTestResult(0, "posit<16,1> != 0.0", "!= long double literal");
	}

	std::cout << "Logic: operator<()" << endl;
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicLessThan<3, 0>(), "posit<3,0>", "<");
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicLessThan<4, 0>(), "posit<4,0>", "<");
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicLessThan<4, 1>(), "posit<4,1>", "<");
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicLessThan<5, 0>(), "posit<5,0>", "<");
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicLessThan<5, 1>(), "posit<5,1>", "<");
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicLessThan<5, 2>(), "posit<5,2>", "<");
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicLessThan<6, 0>(), "posit<6,0>", "<");
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicLessThan<6, 1>(), "posit<6,1>", "<");
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicLessThan<6, 2>(), "posit<6,2>", "<");
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicLessThan<6, 3>(), "posit<6,3>", "<");
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicLessThan<7, 0>(), "posit<7,0>", "<");
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicLessThan<7, 1>(), "posit<7,1>", "<");
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicLessThan<7, 2>(), "posit<7,2>", "<");
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicLessThan<7, 3>(), "posit<7,3>", "<");
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicLessThan<8, 0>(), "posit<8,0>", "<");
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicLessThan<8, 1>(), "posit<8,1>", "<");
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicLessThan<8, 2>(), "posit<8,2>", "<");
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicLessThan<8, 3>(), "posit<8,3>", "<");
	if (p < 0) { 
		nrOfFailedTestCases += ReportTestResult(1, "posit<16,1> < 0", "< int literal");
	}
	else {
		ReportTestResult(0, "posit<16,1> < 0", "< int literal");
	}
	if (p < 0.0f) {
		nrOfFailedTestCases += ReportTestResult(1, "posit<16,1> < 0.0", "< float literal");
	}
	else {
		ReportTestResult(0, "posit<16,1> < 0.0", "< float literal");
	}
	if (p < 0.0) { 
		nrOfFailedTestCases += ReportTestResult(1, "posit<16,1> < 0.0", "< double literal");
	}
	else {
		ReportTestResult(0, "posit<16,1> < 0.0", "< double literal");
	}
	if (p < 0.0l) {
		nrOfFailedTestCases += ReportTestResult(1, "posit<16,1> < 0.0", "< long double literal");
	}
	else {
		ReportTestResult(0, "posit<16,1> < 0.0", "< long double literal");
	}

	std::cout << "Logic: operator<=()" << endl;
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicLessOrEqualThan<3, 0>(), "posit<3,0>", "<=");
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicLessOrEqualThan<4, 0>(), "posit<4,0>", "<=");
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicLessOrEqualThan<4, 1>(), "posit<4,1>", "<=");
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicLessOrEqualThan<5, 0>(), "posit<5,0>", "<=");
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicLessOrEqualThan<5, 1>(), "posit<5,1>", "<=");
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicLessOrEqualThan<5, 2>(), "posit<5,2>", "<=");
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicLessOrEqualThan<6, 0>(), "posit<6,0>", "<=");
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicLessOrEqualThan<6, 1>(), "posit<6,1>", "<=");
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicLessOrEqualThan<6, 2>(), "posit<6,2>", "<=");
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicLessOrEqualThan<6, 3>(), "posit<6,3>", "<=");
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicLessOrEqualThan<7, 0>(), "posit<7,0>", "<=");
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicLessOrEqualThan<7, 1>(), "posit<7,1>", "<=");
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicLessOrEqualThan<7, 2>(), "posit<7,2>", "<=");
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicLessOrEqualThan<7, 3>(), "posit<7,3>", "<=");
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicLessOrEqualThan<8, 0>(), "posit<8,0>", "<=");
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicLessOrEqualThan<8, 1>(), "posit<8,1>", "<=");
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicLessOrEqualThan<8, 2>(), "posit<8,2>", "<=");
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicLessOrEqualThan<8, 3>(), "posit<8,3>", "<=");
	if (!(p <= 0)) {
		nrOfFailedTestCases += ReportTestResult(1, "posit<16,1> <= 0", "<= int literal");
	}
	else {
		ReportTestResult(0, "posit<16,1> <= 0", "<= int literal");
	}
	if (!(p <= 0.0f)) {
		nrOfFailedTestCases += ReportTestResult(1, "posit<16,1> <= 0.0", "<= float literal");
	}
	else {
		ReportTestResult(0, "posit<16,1> <= 0.0", "<= float literal");
	}
	if (!(p <= 0.0)) {
		nrOfFailedTestCases += ReportTestResult(1, "posit<16,1> <= 0.0", "<= double literal");
	}
	else {
		ReportTestResult(0, "posit<16,1> <= 0.0", "<= double literal");
	}
	if (!(p <= 0.0l)) {
		nrOfFailedTestCases += ReportTestResult(1, "posit<16,1> <= 0.0", "<= long double literal");
	}
	else {
		ReportTestResult(0, "posit<16,1> <= 0.0", "<= long double literal");
	}

	std::cout << "Logic: operator>()" << endl;
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicGreaterThan<3, 0>(), "posit<3,0>", ">");
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicGreaterThan<4, 0>(), "posit<4,0>", ">");
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicGreaterThan<4, 1>(), "posit<4,1>", ">");
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicGreaterThan<5, 0>(), "posit<5,0>", ">");
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicGreaterThan<5, 1>(), "posit<5,1>", ">");
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicGreaterThan<5, 2>(), "posit<5,2>", ">");
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicGreaterThan<6, 0>(), "posit<6,0>", ">");
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicGreaterThan<6, 1>(), "posit<6,1>", ">");
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicGreaterThan<6, 2>(), "posit<6,2>", ">");
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicGreaterThan<6, 3>(), "posit<6,3>", ">");
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicGreaterThan<7, 0>(), "posit<7,0>", ">");
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicGreaterThan<7, 1>(), "posit<7,1>", ">");
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicGreaterThan<7, 2>(), "posit<7,2>", ">");
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicGreaterThan<7, 3>(), "posit<7,3>", ">");
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicGreaterThan<8, 0>(), "posit<8,0>", ">");
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicGreaterThan<8, 1>(), "posit<8,1>", ">");
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicGreaterThan<8, 2>(), "posit<8,2>", ">");
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicGreaterThan<8, 3>(), "posit<8,3>", ">");
	if (p > 0) { 
		nrOfFailedTestCases += ReportTestResult(1, "posit<16,1> > 0", "> int literal");
	}
	else {
		ReportTestResult(0, "posit<16,1> > 0", "> int literal");
	}
	if (p > 0.0f) {
		nrOfFailedTestCases += ReportTestResult(1, "posit<16,1> > 0.0", "> float literal");
	}
	else {
		ReportTestResult(0, "posit<16,1> > 0.0", "> float literal");
	}
	if (p > 0.0) {
		nrOfFailedTestCases += ReportTestResult(1, "posit<16,1> > 0.0", "> double literal");
	}
	else {
		ReportTestResult(0, "posit<16,1> > 0.0", "> double literal");
	}
	if (p > 0.0l) {
		nrOfFailedTestCases += ReportTestResult(1, "posit<16,1> > 0.0", "> long double literal");
	}
	else {
		ReportTestResult(0, "posit<16,1> > 0.0", "> long double literal");
	}

	std::cout << "Logic: operator>=()" << endl;
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicGreaterOrEqualThan<3, 0>(), "posit<3,0>", ">=");
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicGreaterOrEqualThan<4, 0>(), "posit<4,0>", ">=");
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicGreaterOrEqualThan<4, 1>(), "posit<4,1>", ">=");
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicGreaterOrEqualThan<5, 0>(), "posit<5,0>", ">=");
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicGreaterOrEqualThan<5, 1>(), "posit<5,1>", ">=");
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicGreaterOrEqualThan<5, 2>(), "posit<5,2>", ">=");
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicGreaterOrEqualThan<6, 0>(), "posit<6,0>", ">=");
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicGreaterOrEqualThan<6, 1>(), "posit<6,1>", ">=");
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicGreaterOrEqualThan<6, 2>(), "posit<6,2>", ">=");
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicGreaterOrEqualThan<6, 3>(), "posit<6,3>", ">=");
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicGreaterOrEqualThan<7, 0>(), "posit<7,0>", ">=");
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicGreaterOrEqualThan<7, 1>(), "posit<7,1>", ">=");
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicGreaterOrEqualThan<7, 2>(), "posit<7,2>", ">=");
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicGreaterOrEqualThan<7, 3>(), "posit<7,3>", ">=");
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicGreaterOrEqualThan<8, 0>(), "posit<8,0>", ">=");
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicGreaterOrEqualThan<8, 1>(), "posit<8,1>", ">=");
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicGreaterOrEqualThan<8, 2>(), "posit<8,2>", ">=");
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicGreaterOrEqualThan<8, 3>(), "posit<8,3>", ">=");
	if (!(p >= 0)) { 
		nrOfFailedTestCases += ReportTestResult(1, "posit<16,1> >= 0", ">= int literal");
	}
	else {
		ReportTestResult(0, "posit<16,1> >= 0", ">= int literal");
	}
	if (!(p >= 0.0f)) {
		nrOfFailedTestCases += ReportTestResult(1, "posit<16,1> >= 0.0", ">= float literal");
	}
	else {
		ReportTestResult(0, "posit<16,1> >= 0.0", ">= float literal");
	}
	if (!(p >= 0.0)) { 
		nrOfFailedTestCases += ReportTestResult(1, "posit<16,1> >= 0.0", ">= double literal");
	}
	else {
		ReportTestResult(0, "posit<16,1> >= 0.0", ">= double literal");
	}
	if (!(p >= 0.0l)) {
		nrOfFailedTestCases += ReportTestResult(1, "posit<16,1> >= 0.0", ">= long double literal");
	}
	else {
		ReportTestResult(0, "posit<16,1> >= 0.0", ">= long double literal");
	}

#if STRESS_TESTING

	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicEqual<16, 1>(), "posit<16,1>", "==");
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicNotEqual<16, 1>(), "posit<16,1>", "!=");
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicLessThan<16, 1>(), "posit<16,1>", "<");
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicLessOrEqualThan<16, 1>(), "posit<16,1>", "<=");
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicGreaterThan<16, 1>(), "posit<16,1>", ">");
	nrOfFailedTestCases += ReportTestResult(ValidatePositLogicGreaterOrEqualThan<16, 1>(), "posit<16,1>", ">=");

#endif // STRESS_TESTING

#endif // MANUAL_TESTING

	return (nrOfFailedTestCases > 0 ? EXIT_FAILURE : EXIT_SUCCESS);
}
catch (char const* msg) {
	std::cerr << msg << std::endl;
	return EXIT_FAILURE;
}
catch (...) {
	std::cerr << "Caught unknown exception" << std::endl;
	return EXIT_FAILURE;
}



