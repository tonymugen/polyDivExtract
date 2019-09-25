/*
 * Copyright (c) 2019 Anthony J. Greenberg
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

/// Extract four-fold synonymous sites from a CDS FASTA file
/** \file
 * \author Anthony J. Greenberg
 * \copyright Copyright (c) 2019 Anthony J. Greenberg
 * \version 0.1
 *
 * Takes a FASTA file with coding sequences (CDS), sorted by chromosome and position by `fastaSort`, and outputs a list of four-fold synonymous sites. Regions that are covered by overlapping CDS are discarded.
 *
 * The flags are:
 *
 * -i input file name
 * -l log file name
 * -o output file name
 *
 */

#include <string>
#include <vector>
#include <unordered_map>
#include <iostream>
#include <fstream>

#include "utilities.hpp"
#include "ffExtract.hpp"

using std::vector;
using std::unordered_map;
using std::cerr;
using std::endl;
using std::fstream;
using std::ios;

using namespace BayesicSpace;

int main(int argc, char *argv[]){
	unordered_map<char, string> clInfo;
	parseCL(argc, argv, clInfo);
	if ( clInfo['i'].empty() ) {
		cerr << "Must specify a FASTA input file with flag -i" << endl;
		exit(1);
	} else if ( clInfo['o'].empty() ) {
		cerr << "Must specify output file name with flag -o" << endl;
		exit(2);
	} else if ( clInfo['l'].empty() ) {
		cerr << "Must specify the log file name with flag -l" << endl;
	}
	try {
		FFextract fasta(clInfo['i'], clInfo['l']);
		vector<string> out;
		fasta.extractFFsites(out);
		fstream oFS;
		oFS.open(clInfo['o'].c_str(), ios::out|ios::trunc);
		oFS << "chr\tFBgn\tpos" << endl;
		for (auto &r : out) {
			oFS << r << endl;
		}
		oFS.close();

	} catch(string error) {
		cerr << error << endl;
		exit(1);
	}
	exit(0);
}

