// arithmetic_reciprocate.cpp: functional tests for arithmetic reciprocation
//
// Copyright (C) 2017-2018 Stillwater Supercomputing, Inc.
//
// This file is part of the universal numbers project, which is released under an MIT Open Source license.

#include "common.hpp"

// when you define POSIT_VERBOSE_OUTPUT executing an reciprocate the code will print intermediate results
//#define POSIT_VERBOSE_OUTPUT
#define POSIT_TRACE_RECIPROCATE 0
#define POSIT_TRACE_CONVERT

// minimum set of include files to reflect source code dependencies
#include "../../posit/posit.hpp"
#include "../../posit/posit_manipulators.hpp"
#include "/home/gps/Downloads/My_test_POSIT/tests/test_helpers.hpp"
#include "/home/gps/Downloads/My_test_POSIT/tests/posit_test_helpers.hpp"
//#include "/home/gps/Downloads/My_test_POSIT/posit/trace_constants.hpp"
// generate specific test case that you can trace with the trace conditions in posit.h
// for most bugs they are traceable with _trace_conversion and _trace_add
template<size_t nbits, size_t es, typename Ty>
void GenerateTestCase(Ty a) {
	Ty reference;
	sw::unum::posit<nbits, es> pa, pref, preciprocal;
	pa = a;
	reference = (Ty)1.0 / a;
	pref = reference;
	preciprocal = pa.reciprocate();
	std::cout << "input " << a << " reference 1/fa " << reference << " pref " << pref << " result " << preciprocal << std::endl << std::endl;
}

#define MANUAL_TESTING 0
#define STRESS_TESTING 0

int main(int argc, char** argv)
try {
	using namespace std;
	using namespace sw::unum;

	bool bReportIndividualTestCases = false;
	int nrOfFailedTestCases = 0;

	cout << "Posit reciprocate validation" << endl;

	std::string tag = "Reciprocation failed: ";

#if MANUAL_TESTING

	// generate individual testcases to hand trace/debug

	GenerateTestCase<4, 0, double>(0.75);
	GenerateTestCase<5, 0, double>(0.75);
	GenerateTestCase<6, 0, double>(0.75);
	GenerateTestCase<16, 0, double>(0.75);
	posit<16, 0> p(1 / 0.75);
	cout << p.get() << " " << pretty_print(p, 17) << endl;

	nrOfFailedTestCases += ReportTestResult(ValidateReciprocation<3, 0>("Manual testing", true), "posit<3,0>", "reciprocation");
	nrOfFailedTestCases += ReportTestResult(ValidateReciprocation<4, 0>("Manual testing", true), "posit<4,0>", "reciprocation");
	nrOfFailedTestCases += ReportTestResult(ValidateReciprocation<5, 0>("Manual testing", true), "posit<5,0>", "reciprocation");
	nrOfFailedTestCases += ReportTestResult(ValidateReciprocation<6, 0>("Manual testing", true), "posit<6,0>", "reciprocation");

	nrOfFailedTestCases += ReportTestResult(ValidateReciprocation<5, 1>("Manual testing", true), "posit<5,1>", "reciprocation");
	nrOfFailedTestCases += ReportTestResult(ValidateReciprocation<6, 1>("Manual testing", true), "posit<6,1>", "reciprocation");
	nrOfFailedTestCases += ReportTestResult(ValidateReciprocation<7, 1>("Manual testing", true), "posit<7,1>", "reciprocation");

	nrOfFailedTestCases += ReportTestResult(ValidateReciprocation<8, 2>(tag, true), "posit<8,2>", "reciprocation");

#else

	nrOfFailedTestCases += ReportTestResult(ValidateReciprocation<3, 0>(tag, bReportIndividualTestCases), "posit<3,0>", "reciprocation");

	nrOfFailedTestCases += ReportTestResult(ValidateReciprocation<4, 0>(tag, bReportIndividualTestCases), "posit<4,0>", "reciprocation");
	nrOfFailedTestCases += ReportTestResult(ValidateReciprocation<4, 1>(tag, bReportIndividualTestCases), "posit<4,1>", "reciprocation");

	nrOfFailedTestCases += ReportTestResult(ValidateReciprocation<5, 0>(tag, bReportIndividualTestCases), "posit<5,0>", "reciprocation");
	nrOfFailedTestCases += ReportTestResult(ValidateReciprocation<5, 1>(tag, bReportIndividualTestCases), "posit<5,1>", "reciprocation");
	nrOfFailedTestCases += ReportTestResult(ValidateReciprocation<5, 2>(tag, bReportIndividualTestCases), "posit<5,2>", "reciprocation");

	nrOfFailedTestCases += ReportTestResult(ValidateReciprocation<6, 0>(tag, bReportIndividualTestCases), "posit<6,0>", "reciprocation");
	nrOfFailedTestCases += ReportTestResult(ValidateReciprocation<6, 1>(tag, bReportIndividualTestCases), "posit<6,1>", "reciprocation");
	nrOfFailedTestCases += ReportTestResult(ValidateReciprocation<6, 2>(tag, bReportIndividualTestCases), "posit<6,2>", "reciprocation");
	nrOfFailedTestCases += ReportTestResult(ValidateReciprocation<6, 3>(tag, bReportIndividualTestCases), "posit<6,3>", "reciprocation");

	nrOfFailedTestCases += ReportTestResult(ValidateReciprocation<7, 0>(tag, bReportIndividualTestCases), "posit<7,0>", "reciprocation");
	nrOfFailedTestCases += ReportTestResult(ValidateReciprocation<7, 1>(tag, bReportIndividualTestCases), "posit<7,1>", "reciprocation");
	nrOfFailedTestCases += ReportTestResult(ValidateReciprocation<7, 2>(tag, bReportIndividualTestCases), "posit<7,2>", "reciprocation");
	nrOfFailedTestCases += ReportTestResult(ValidateReciprocation<7, 3>(tag, bReportIndividualTestCases), "posit<7,3>", "reciprocation");
	nrOfFailedTestCases += ReportTestResult(ValidateReciprocation<7, 4>(tag, bReportIndividualTestCases), "posit<7,4>", "reciprocation");

	nrOfFailedTestCases += ReportTestResult(ValidateReciprocation<8, 0>(tag, bReportIndividualTestCases), "posit<8,0>", "reciprocation");
	nrOfFailedTestCases += ReportTestResult(ValidateReciprocation<8, 1>(tag, bReportIndividualTestCases), "posit<8,1>", "reciprocation");
	nrOfFailedTestCases += ReportTestResult(ValidateReciprocation<8, 2>(tag, bReportIndividualTestCases), "posit<8,2>", "reciprocation");
	nrOfFailedTestCases += ReportTestResult(ValidateReciprocation<8, 3>(tag, bReportIndividualTestCases), "posit<8,3>", "reciprocation");
	nrOfFailedTestCases += ReportTestResult(ValidateReciprocation<8, 4>(tag, bReportIndividualTestCases), "posit<8,4>", "reciprocation");
	nrOfFailedTestCases += ReportTestResult(ValidateReciprocation<8, 5>(tag, bReportIndividualTestCases), "posit<8,5>", "reciprocation");

	nrOfFailedTestCases += ReportTestResult(ValidateReciprocation<10, 1>(tag, bReportIndividualTestCases), "posit<10,1>", "reciprocation");
	nrOfFailedTestCases += ReportTestResult(ValidateReciprocation<12, 1>(tag, bReportIndividualTestCases), "posit<12,1>", "reciprocation");
	nrOfFailedTestCases += ReportTestResult(ValidateReciprocation<16, 1>(tag, bReportIndividualTestCases), "posit<16,1>", "reciprocation");

#if STRESS_TESTING

	nrOfFailedTestCases += ReportTestResult(ValidateReciprocation<20, 1>(tag, bReportIndividualTestCases), "posit<20,1>", "reciprocation");

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

