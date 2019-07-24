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

/// Extract polymorphic sites
/** \file
 * \author Anthony J. Greenberg
 * \copyright Copyright (c) 2019 Anthony J. Greenberg
 * \version 0.1
 *
 * Extracts polymorphic sites from MSL complex peak ranges and the four-fold silent site file list. Outgroups for ancestral state determination is either _D. simulans_ or _D. yakuba_.
 * The flags are:
 *
 * -q query file name (binding locations or four-fold sites)
 * -a .axt file name (for the outgroup)
 * -v VCF file name
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

#include "parseVCF.hpp"
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
		}  else if ( clInfo['v'].empty() ) {
			throw string("Must specify input file with flag -q");
		} else if ( clInfo['o'].empty() ) {
			throw string("Must specify output file name with flag -o");
		}

		ParseVCF vcf(clInfo['v'], clInfo['a']);

		if (clInfo['q'].find("fourFold.tsv") != string::npos) {  // found the signature of the four-fold silent site file
			string qLine;
			vector<string> chrNams;
			vector<uint64_t> positions;

			fstream queryFile;
			queryFile.open(clInfo['q'].c_str(), ios::in);

			while( getline(queryFile, qLine) ){
				//getline(queryFile, qLine);
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

			vector<string> polySites;
			vcf.getPolySites(chrNams, positions, polySites);

			fstream outFile;
			outFile.open(clInfo['o'].c_str(), ios::out | ios::trunc);

			// output the results
			outFile << "CHR\tPOS\tREF\tALT\tANC\tAC\tMLAC\tAF\tMLAF\tNMISS\tSAME_CHR\tOUTQUAL\tSITEQUAL" << endl;
			for (auto &ps : polySites) {
				outFile << ps << endl;
			}
			outFile.close();
		} else {
			string qLine;

			fstream queryFile;
			queryFile.open(clInfo['q'].c_str(), ios::in);

			fstream outFile;
			outFile.open(clInfo['o'].c_str(), ios::out | ios::trunc);
			outFile << "PEAK_ID\tCHR\tPOS\tREF\tALT\tANC\tAC\tMLAC\tAF\tMLAF\tNMISS\tSAME_CHR\tOUTQUAL\tSITEQUAL" << endl;

			uint32_t peakID = 1;
			while( getline(queryFile, qLine) ){
				//getline(queryFile, qLine);
				stringstream lnSS(qLine);
				vector<string> fields;
				string value;
				while(lnSS >> value){
					fields.push_back(value);
				}
				if (fields.size() != 5) {
					throw string("Peak files should have five fields");
				}
				vector<string> polySites;
				vcf.getPolySites(fields[0], strtoul(fields[1].c_str(), NULL, 0), strtoul(fields[2].c_str(), NULL, 0), polySites);
				for (auto &ps : polySites) {
					outFile << "P" << peakID << "\t" << ps << endl;
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


