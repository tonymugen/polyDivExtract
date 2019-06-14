/*
 * Copyright (c) <YEAR> Anthony J. Greenberg
 *
 * Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS
 * BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
 * IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
 * THE POSSIBILITY OF SUCH DAMAGE.
 */

/// Utility functions
/** \file
 * \author Anthony J. Greenberg
 * \copyright Copyright (c) 2019 Anthony J. Greenberg
 * \version 0.1
 *
 */
#ifndef utilities_hpp
#define utilities_hpp

#include <string>
#include <unordered_map>

using std::string;
using std::unordered_map;

namespace BayesicSpace {
	/** \brief Parse command line flags
	 *
	 * \param[in] argc number of arguments
	 * \param[in] argv array of argument values
	 * \param[out] cli flag values, indexed by flag IDs
	 */
	void parseCL(int &argc, char **argv, unordered_map<char, string> &cli){
		// set to true after encountering a flag token (the character after the dash)
		bool val = false;
		// store the token value here
		char curFlag;

		for (int iArg = 1; iArg < argc; iArg++) {
			const char *pchar = argv[iArg];

			if (pchar[0] == '-') { // encountered the dash, look for the token after it
				if (!pchar[1]) {
					throw string("ERROR: forgot character after dash");
				}
				// what follows the dash?
				val     = true;
				curFlag = pchar[1];

			} else {
				if (val) {
					val = false;
					cli[curFlag] = pchar;
				}
			}

		}
	}
}
#endif /* utilities_hpp */

