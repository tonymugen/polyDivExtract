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

/// Parse .axt aligniment files
/** \file
 * \author Anthony J. Greenberg
 * \copyright Copyright (c) 2019 Anthony J. Greenberg
 * \version 0.1
 *
 * Class implementation for parsing of .axt alignement files to extract ancestral nucleotide and divergence information
 *
 */

#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <system_error>

#include <iostream>

#include "parseAXT.hpp"

using std::fstream;
using std::ofstream;
using std::stringstream;
using std::string;
using std::vector;
using std::endl;
using std::flush;
using std::system_error;
using std::ios;
using std::bad_alloc;

using namespace BayesicSpace;

ParseAXT::ParseAXT(const string &fileName) : chrID_{""}, sameChr_{0}, primaryStart_{0}, primaryEnd_{0}, alignedStart_{0}, alignedEnd_{0}, primarySeq_{""}, alignSeq_{""} {
	if( axtFile_.is_open() ){
		axtFile_.close();
	}

	try {
		axtFile_.open(fileName.c_str(), ios::in);
	} catch(system_error &error) {
		string message = "ERROR: cannot open file " + fileName + " to read: " + error.code().message();
		throw message;
	}

	getNextRecord_();
}

string ParseAXT::getMetaData(){
	stringstream outLine;

	outLine << chrID_ + " ";
	outLine << sameChr_;
	outLine << " ";
	outLine << primaryStart_;
	outLine << " ";
	outLine << primaryEnd_;
	outLine << " ";
	outLine << alignedStart_;
	outLine << " ";
	outLine << alignedEnd_;

	return outLine.str();
}

void ParseAXT::getSiteStates(const string &chromosome, const uint64_t &position, char &primaryState, char &alignedState, uint16_t &sameChromosome){
	bool noneFound = true;
	while(axtFile_){
		if (chrID_ != chromosome) {
			getNextRecord_();
			continue;
		}
		if (primaryEnd_ >= position) {
			uint64_t truePos = primaryStart_;                   // this is the genomic position (with gaps eliminated)
			for (size_t i = 0; i < primarySeq_.size() ; i++) {  // string length equality already checked in getNextRecord_()
				if (primarySeq_[i] == '-') {
					continue;
				}
				if (truePos == position) {
					primaryState   = primarySeq_[i];
					alignedState   = alignSeq_[i];   // may be a gap, that can be checked in post-processing
					sameChromosome = sameChr_;
					noneFound      = false;
					break;
				}
			truePos++;
			}
			break;
		} else {
			getNextRecord_();
		}

	}
	if (noneFound) {
		stringstream wrongThing;
		wrongThing << "Reached the end of file before finding a record for postition ";
		wrongThing << position;
		wrongThing << " on chromosome ";
		wrongThing << chromosome;
		throw wrongThing.str();
	}
}

void ParseAXT::getNextRecord_(){
	string curLine("");
	while(axtFile_){
		getline(axtFile_, curLine);
		if (curLine[0] == '#') {
			continue;
		} else if (curLine == "") {
			continue;
		} else {
			break;
		}
	}
	if (curLine == ""){
		throw string("End of file");
	}

	// we have a non-empty line, presumably the meta-data header for .axt
	stringstream metaSS(curLine);
	vector<string> fields;
	string value;
	while(metaSS >> value){ // makes sure we don't have anything extra at the end
		fields.push_back(value);
	}
	if (fields.size() != 9) {
		throw string("Wrong number of fields in .axt metada");
	}
	if ( (fields[1][0] != 'c') || (fields[1][1] != 'h') || (fields[1][2] != 'r') ) { // do not have "chr" at the beginning of the chromosome field
		string wrongThing = "Wrong chromosome field: " + fields[1];
		throw wrongThing;
	}
	string tmpChrID = chrID_;
	chrID_ = fields[1];

	uint64_t tmpStart = primaryStart_;
	primaryStart_ = strtoul(fields[2].c_str(), NULL, 0);
	if (primaryStart_ == 0) {
		string wrongThing = "Wrong primary sequence start: " + fields[2];
		throw wrongThing;
	} else if ( (primaryStart_ <= tmpStart) && (tmpChrID == chrID_) ) { // the records should be in order of increasing primary sequence position, unless there is a chromosome switch
		string wrongThing = "Primary start of the current record (" + fields[2] + ") not greater than the perivous record";
		throw wrongThing;
	}
	primaryEnd_ = strtoul(fields[3].c_str(), NULL, 0);
	if (primaryEnd_ == 0) {
		string wrongThing = "Wrong primary sequence end: " + fields[3];
		throw wrongThing;
	} else if (primaryEnd_ <= primaryStart_) {
		stringstream wrongThing;
		wrongThing << "Position of the end of primary sequence (";
		wrongThing << fields[3];
		wrongThing << ") not greater than the position of the start: ";
		wrongThing << primaryStart_;
		wrongThing << " (if the last number is huge, the original number was likely negative)";
		throw wrongThing.str();
	}

	if ( (fields[4][0] != 'c') || (fields[4][1] != 'h') || (fields[4][2] != 'r') ) { // do not have "chr" at the beginning of the chromosome field
		string wrongThing = "Wrong aligned chromosome field: " + fields[4];
		throw wrongThing;
	}
	sameChr_ = (fields[4] == fields[1] ? 1 : 0);

	alignedStart_ = strtoul(fields[5].c_str(), NULL, 0);
	if (alignedStart_ == 0) {
		string wrongThing = "Wrong primary sequence start: " + fields[5];
		throw wrongThing;
	}
	alignedEnd_ = strtoul(fields[6].c_str(), NULL, 0);
	if (alignedEnd_ == 0) {
		string wrongThing = "Wrong primary sequence end: " + fields[6];
		throw wrongThing;
	}

	// now just read the sequences
	if(axtFile_.eof()){
		throw string("End of file reached before primary sequence read");
	}
	getline(axtFile_, primarySeq_);

	if(axtFile_.eof()){
		throw string("End of file reached before aligned sequence read");
	}
	getline(axtFile_, alignSeq_);
	if ( primarySeq_.size() != alignSeq_.size() ) {
		string wrongThing = "The sequence strings for record #" + fields[0] + " are not equal length";
		throw wrongThing;
	}
}

