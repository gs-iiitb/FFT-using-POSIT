// 48bit_posit.cpp: Functionality tests for extended standard 48-bit posits
//
// Copyright (C) 2017-2018 Stillwater Supercomputing, Inc.
//
// This file is part of the universal numbers project, which is released under an MIT Open Source license.

#include "common.hpp"

#include <vector>
#include <posit>

#include "../tests/test_helpers.hpp"
#include "../tests/posit_test_helpers.hpp"

/*
Extended Standard posit with nbits = 48 have es = 2 exponent bits.
*/

int main(int argc, char** argv)
try {
	using namespace std;
	using namespace sw::unum;

	const size_t RND_TEST_CASES = 0; // 150000;

	const size_t nbits = 48;
	const size_t es = 2;

	int nrOfFailedTestCases = 0;
	bool bReportIndividualTestCases = false;
	std::string tag = " posit<48,2>";

	cout << "Extended Standard posit<48,2> configuration tests" << endl;

	posit<nbits, es> p;
	cout << spec_to_string(p) << endl << endl;

	cout << "Arithmetic tests " << RND_TEST_CASES << " randoms each" << endl;
	nrOfFailedTestCases += ReportTestResult(ValidateThroughRandoms<nbits, es>(tag, bReportIndividualTestCases, OPCODE_ADD, RND_TEST_CASES), tag, "addition      ");
	nrOfFailedTestCases += ReportTestResult(ValidateThroughRandoms<nbits, es>(tag, bReportIndividualTestCases, OPCODE_SUB, RND_TEST_CASES), tag, "subtraction   ");
	nrOfFailedTestCases += ReportTestResult(ValidateThroughRandoms<nbits, es>(tag, bReportIndividualTestCases, OPCODE_MUL, RND_TEST_CASES), tag, "multiplication");
	nrOfFailedTestCases += ReportTestResult(ValidateThroughRandoms<nbits, es>(tag, bReportIndividualTestCases, OPCODE_DIV, RND_TEST_CASES), tag, "division      ");

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

