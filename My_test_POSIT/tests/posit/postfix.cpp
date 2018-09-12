// postfix.cpp functional tests for postfix operators
//
// Copyright (C) 2017-2018 Stillwater Supercomputing, Inc.
//
// This file is part of the universal numbers project, which is released under an MIT Open Source license.

#include "common.hpp"
#include <vector>
#include <algorithm>

#include "../../posit/posit.hpp"
#include "../../posit/posit_manipulators.hpp"
#include "../tests/test_helpers.hpp"
#include "../tests/posit_test_helpers.hpp"

int main(int argc, char** argv)
try {
	using namespace std;
	using namespace sw::unum;

	bool bReportIndividualTestCases = false;
	int nrOfFailedTestCases = 0;

	nrOfFailedTestCases += ReportTestResult(ValidatePostfix<3, 0>("Postfix failed", bReportIndividualTestCases), "posit<3,0>", "posit++");

	nrOfFailedTestCases += ReportTestResult(ValidatePostfix<4, 0>("Postfix failed", bReportIndividualTestCases), "posit<4,0>", "posit++");
	nrOfFailedTestCases += ReportTestResult(ValidatePostfix<4, 1>("Postfix failed", bReportIndividualTestCases), "posit<4,1>", "posit++");

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

#ifdef TEACHING_MOMENT
// just because you can, doesn't mean you should
void DoNotDoStuffLikeThis() {
	using namespace std;
	using namespace sw::unum;

	// order of function evaluation is not defined, so there is no
	// mechanism for these methods to evaluate left to right
	// integer example
	// DON'T
	int i = 0;
	cout << i << " " << --i << " " << i << " " << i++ << " " << i << endl;
	cout << i << " " << i++ << " " << i << endl;
	i = 0;
	cout << --i << " " << --i << " " << --i << endl;
	i = 0;
	cout << --(--(--i)) << endl;
	i = 0;
	cout << ------i << " " << endl;

	// equivalent posit example
	const size_t nbits = 4;
	const size_t es = 0;
	posit<nbits, es> result, p = 0.0f;
	cout << p << " " << --p << " " << p << " " << p++ << " " << p << endl;
	cout << p << " " << p++ << " " << p << endl;
	p = 0.0f;
	cout << --p << " " << --p << " " << --p << endl;
	p = 0.0f;
	cout << --(--(--p)) << endl;


	p = 0.0f;
	result = --p++;
	cout << "result " << result << endl;

	int nrOfFailedTestCases = 0;
	p = 0.0f;
	if (!p.isZero()) {
		cout << "FAIL 1 " << p << endl; nrOfFailedTestCases++;
	}
	p = 0.0f; --(--(--(p++)++)++);
	if (!p.isZero()) {
		cout << "FAIL 2 " << p << endl; nrOfFailedTestCases++;
	}
	p = 0.0f; ++(++(++(p--)--)--);
	if (!p.isZero()) {
		cout << "FAIL 3 " << p << endl; nrOfFailedTestCases++;
	}
	p = 0.0f; ----------p++++++++++;
	if (!p.isZero()) {
		cout << "FAIL 4 " << p << endl; nrOfFailedTestCases++;
	}
	p = 0.0f; p++++++++++;
	if (p != posit<nbits, es>(1.0f)) {
		cout << "FAIL 5 " << p << endl; nrOfFailedTestCases++;
	}
}
#endif
