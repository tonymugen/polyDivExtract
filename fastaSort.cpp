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

/// Sort a FASTA file by chromosome and position
/** \file
 * \author Anthony J. Greenberg
 * \copyright Copyright (c) 2019 Anthony J. Greenberg
 * \version 0.1
 *
 * Sorts a FASTA file that has a _loc=_ field in the header of each sequence by position (of the start nucleotide) within each chromosome. If records with the same position are found, the longest one is kept.
 * The flags are:
 *
 * -i input file name
 * -o output file name
 *
 */
#include <string>
#include <vector>
#include <map>
#include <unordered_map>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cctype>

#include "utilities.hpp"

using std::vector;
using std::map;
using std::unordered_map;
using std::pair;
using std::cout;
using std::cerr;
using std::endl;
using std::flush;
using std::fstream;
using std::ofstream;
using std::stringstream;
using std::ios;
using BayesicSpace::parseCL;

int main(int argc, char *argv[]){
	unordered_map<char, string> clInfo;
	parseCL(argc, argv, clInfo);
	if ( clInfo['i'].empty() ) {
		cerr << "Must specify a FASTA input file with flag -i" << endl;
		exit(1);
	} else if ( clInfo['o'].empty() ) {
		cerr << "Must specify output file name with flag -o" << endl;
		exit(2);
	}
	// map that relates chromosome arms to data
	// the data are a map of position and a vector of two strings that has the FASTA record (the header at 0 and sequence at 1, sequence all on one line)
	map<string, map<uint64_t, vector<string>>> outData;
	fstream fastaIn;
	fastaIn.open(clInfo['i'].c_str(), ios::in);
	string curLine;               // current line of the input file
	vector<string> curRecord(2);  // the current record (header and sequence; sequence all on one line) to be added to the map; first element is the header, second is the sequence
	uint64_t curPos;              // position of the start of the current sequence
	string curChr;                // current chromosome name
	while( getline(fastaIn, curLine) ){
		if (curLine[0] == '>') { // we are at the FASTA header
			if ( !curRecord[0].empty() ) {  // deal with the previous record (if any)
				pair<map<uint64_t, vector<string>>::iterator, bool> successTest; // success of insertion test
				if (curChr.compare(0, 4, "Scr_") == 0) { // soecial case of chromosome naming in the Dsim CDS FASTA file
					curChr.erase(0,4);
				}
				if ( (curChr == "X") ||
						(curChr == "2L") || (curChr == "2R") ||
						(curChr == "3L") || (curChr == "3R") ||
						(curChr == "4") ) {
					successTest = outData[curChr].insert( pair<uint64_t, vector<string>>(curPos, curRecord) );
					// if a record at this position already exists, insert a new record if the its sequence is longer
					if (!successTest.second & ( curRecord[1].size() > successTest.first->second[1].size() ) ) {
						successTest.first->second = curRecord;
					}
				}
			}
			curChr.clear();
			// now parse the header and populate the current variables
			curRecord[0] = curLine;
			curRecord[1].clear();  // reset the sequence so that it's ready for appending
			stringstream hSS(curLine);
			string field;
			while(getline(hSS, field, ' ')){
				if ( field.compare(0, 4, "loc=") == 0 ) {
					field.erase(0, 4);
					for (auto &f : field) {
						if (f == ':') {
							break;
						}
						curChr += f;
					}
					cout << curChr << "; ";
					string startPosition;
					for (auto &c : field) {
						if ( isdigit(c) ) {
							startPosition += c;
						} else if (c == '.'){
							break;
						}
					}
					curPos = strtoul(startPosition.c_str(), NULL, 0);
					cout << curPos << endl;
				}
			}
		} else {
			curRecord[1] += curLine;
		}
	}
	fastaIn.close();
	// now save the results
	fstream fastaOut;
	fastaOut.open(clInfo['o'].c_str(), ios::out|ios::trunc);
	for (auto &chr : outData) {
		for (auto &r : chr.second) {
			fastaOut << r.second[0] << endl;
			fastaOut << r.second[1] << endl;
		}
	}
	fastaOut.close();
}


