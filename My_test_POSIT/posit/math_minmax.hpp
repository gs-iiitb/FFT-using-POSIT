#pragma once
// math_minmax.hpp: min/max functions for posits
//
// Copyright (C) 2017-2018 Stillwater Supercomputing, Inc.
//
// This file is part of the universal numbers project, which is released under an MIT Open Source license.


namespace sw {
	namespace unum {

		// the current shims are NON-COMPLIANT with the posit standard, which says that every function must be
		// correctly rounded for every input value. Anything less sacrifices bitwise reproducibility of results.

		template<size_t nbits, size_t es>
		posit<nbits,es> min(posit<nbits,es> x, posit<nbits, es> y) {
			return posit<nbits,es>(std::pow(double(x), double(y)));
		}
		
		template<size_t nbits, size_t es>
		posit<nbits,es> max(posit<nbits,es> x, posit<nbits, es> y) {
			return posit<nbits,es>(std::pow(double(x), double(y)));
		}

	}  // namespace unum

}  // namespace sw
