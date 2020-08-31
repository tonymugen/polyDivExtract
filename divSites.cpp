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

/// Extract diverged sites
/** \file
 * \author Anthony J. Greenberg
 * \copyright Copyright (c) 2019 Anthony J. Greenberg
 * \version 0.1
 *
 * Extracts divergent sites from MSL complex peak ranges and the four-fold silent site file list. The divergence is either to _D. simulans_ or _D. yakuba_.
 * The flags are:
 *
 * -q query file name (binding locations or four-fold sites)
 * -a .axt file name
 * -o output file name
 *
 */

#include <string>
#include <vector>
#include <unordered_map>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cctype>

#include "parseAXT.hpp"
#include "utilities.hpp"

using namespace BayesicSpace;
using std::vector;
using std::unordered_map;
using std::cout;
using std::cerr;
using std::endl;
using std::flush;
using std::fstream;
using std::ofstream;
using std::stringstream;
using std::ios;

using namespace BayesicSpace;

int main(int argc, char *argv[]){
	try {
		unordered_map<char, string> clInfo;
		parseCL(argc, argv, clInfo);
		if ( clInfo['a'].empty() ) {
			throw string("Must specify .axt file with flag -a");
		} else if ( clInfo['q'].empty() ) {
			throw string("Must specify input file with flag -q");
		} else if ( clInfo['o'].empty() ) {
			throw string("Must specify output file name with flag -o");
		}

		ParseAXT axt(clInfo['a']);

		if (clInfo['q'].find("_ff") != string::npos) {  // found the signature of the four-fold silent site file
			string qLine;
			vector<string> chrNams;
			vector<uint64_t> positions;

			fstream queryFile;
			queryFile.open(clInfo['q'].c_str(), ios::in);

			while( getline(queryFile, qLine) ){
				stringstream lnSS(qLine);
				vector<string> fields;
				string value;
				while(lnSS >> value){
					fields.push_back(value);
				}
				if (fields.size() != 2) {
					throw string("Four-fold file should have two fields");
				}
				chrNams.push_back(fields[0]);
				positions.push_back( strtoul(fields[1].c_str(), NULL, 0) );
			}
			queryFile.close();

			vector<string> divergedSites;
			unordered_map<string, uint64_t> lengths;
			axt.getDivergedSites(chrNams, positions, divergedSites, lengths);

			fstream outFile;
			outFile.open(clInfo['o'].c_str(), ios::out | ios::trunc);

			// first put meta-data (total number of good sites) in commented lines at the beginning of the file
			for (auto &c : lengths) {
				outFile << "#\t" << c.first << "\t" << c.second << endl;
			}

			// now output the results
			outFile << "chr\tposition\tprNuc\talNuc\tsameCHR\tgoodQual" << endl;
			for (auto &ds : divergedSites) {
				outFile << ds << endl;
			}
			outFile.close();
		} else {
			string qLine;

			fstream queryFile;
			queryFile.open(clInfo['q'].c_str(), ios::in);

			fstream outFile;
			outFile.open(clInfo['o'].c_str(), ios::out | ios::trunc);
			outFile << "peakID\trealLen\tchr\tposition\tprNuc\talNuc\tsameCHR\tgoodQual" << endl;

			uint32_t peakID = 1;
			while( getline(queryFile, qLine) ){
				//getline(queryFile, qLine);
				stringstream lnSS(qLine);
				vector<string> fields;
				string value;
				while(lnSS >> value){
					fields.push_back(value);
				}
				if (fields.size() < 3) {
					throw string("Peak files should have at least three fields");
				}
				uint64_t length = 0;
				vector<string> divergedSites;
				axt.getDivergedSites(fields[0], strtoul(fields[1].c_str(), NULL, 0), strtoul(fields[2].c_str(), NULL, 0), divergedSites, length);
				for (auto &ds : divergedSites) {
					outFile << "P" << peakID << "\t" << length << "\t" << ds << endl;
				}
				peakID++;
			}
			queryFile.close();
			outFile.close();
		}

	} catch(string error) {
		cerr << error << endl;
	}
}


