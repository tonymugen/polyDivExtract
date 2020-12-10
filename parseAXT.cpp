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
#include <unordered_map>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cctype>
#include <system_error>

#include "parseAXT.hpp"

using std::fstream;
using std::ofstream;
using std::stringstream;
using std::string;
using std::vector;
using std::unordered_map;
using std::system_error;
using std::ios;

using namespace BayesicSpace;

ParseAXT::ParseAXT(const string &fileName) : sameChr_{0}, primaryStart_{0}, primaryEnd_{0}, alignedStart_{0}, alignedEnd_{0}, chrID_{""}, primarySeq_{""}, alignSeq_{""}, foundChr_{""} {
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
ParseAXT &ParseAXT::operator=(ParseAXT &&in){
	if(&in != this){
		axtFile_      = move(in.axtFile_);
		sameChr_      = in.sameChr_;
		primaryStart_ = in.primaryStart_;
		primaryEnd_   = in.primaryEnd_;
		alignedStart_ = in.alignedStart_;
		alignedEnd_   = in.alignedEnd_;
		chrID_        = move(in.chrID_);
		primarySeq_   = move(in.primarySeq_);
		alignSeq_     = move(in.alignSeq_);
		foundChr_     = move(in.foundChr_);

	}

	return(*this);
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

void ParseAXT::getDivergedSites(const string &chromName, const uint64_t &start, const uint64_t &end, vector<string> &sites, uint64_t &length){
	if (start >= end) {
		stringstream wrongThing;
		wrongThing << "ERROR: start position (";
		wrongThing << start;
		wrongThing << ") must come before the end postion (";
		wrongThing << end;
		wrongThing << ") in getDivergedSites()";
		throw wrongThing.str();
	}
	length = 0;
	if ( sites.size() ){
		sites.clear();
	}
	for (uint64_t iSite = start; iSite <= end; iSite++) {
		// if the current chromosome has already been explored to the end, no need to bother looking
		if (chromName == foundChr_) {
			return;
		}
		char primary;
		char aligned;
		uint16_t same;
		getSiteStates_(chromName, iSite, primary, aligned, same);  // will search .axt records
		if ( (primary == '-') || (aligned == '-') ) {  // gaps present; ignore
			continue;
		}
		if ( (primary == 'n') || (aligned == 'n') ) {  // unkown nucleotide present; ignore
			continue;
		}
		if ( (primary == 'N') || (aligned == 'N') ) {  // unkown nucleotide present; ignore
			continue;
		}
		if (primary == aligned) {
			length++;
		} else if ( toupper(primary) == toupper(aligned) ) {  // sometimes there are lower-case bases (low-quality I think)
			length++;
		} else {  // the sites are divergent
			stringstream siteInfo;
			siteInfo << chromName << "\t";
			siteInfo << iSite << "\t";
			siteInfo << primary << "\t" << aligned << "\t";
			siteInfo << same << "\t";
			if ( isupper(primary) && isupper(aligned) ) {
				siteInfo << "1";
			} else {
				siteInfo << "0";
			}
			sites.push_back( siteInfo.str() );
			length++;
		}
	}
}

void ParseAXT::getDivergedSites(const vector<string> &chromNames, const vector<uint64_t> &positions, vector<string> &sites, unordered_map<string, uint64_t> &lengths){
	if (positions.size() != chromNames.size()) {
		stringstream wrongThing;
		wrongThing << "ERROR: the vector of chromosome names (size = ";
		wrongThing << chromNames.size();
		wrongThing << ") not the same size as the vector of positions (size = ";
		wrongThing << positions.size();
		wrongThing << ") in getDivergedSites()";
		throw wrongThing.str();
	}
	for (uint64_t iPos = 0; iPos < positions.size(); iPos++) {
		// if the current chromosome has already been explored to the end, skip (don't return b/c other chromosomes still may be on the list)
		if (chromNames[iPos] == foundChr_) {
			continue;
		}
		char primary;
		char aligned;
		uint16_t same;
		getSiteStates_(chromNames[iPos], positions[iPos], primary, aligned, same); // will search the .axt records
		if ( (primary == '-') || (aligned == '-') ) {  // gaps present; ignore
			continue;
		}
		if ( (primary == 'n') || (aligned == 'n') ) {  // unkown nucleotide present; ignore
			continue;
		}
		if ( (primary == 'N') || (aligned == 'N') ) {  // unkown nucleotide present; ignore
			continue;
		}
		if (primary == aligned) {
			uint64_t tmpSz = lengths.size(); // save the pre-insertion size
			lengths[ chromNames[iPos] ]++;
			if (lengths.size() > tmpSz) {    // if a new element is added, make sure the count starts at 1
				lengths[ chromNames[iPos] ] = 1;
			}
		} else if ( toupper(primary) == toupper(aligned) ) {  // sometimes there are lower-case bases (low-quality I think)
			uint64_t tmpSz = lengths.size(); // save the pre-insertion size
			lengths[ chromNames[iPos] ]++;
			if (lengths.size() > tmpSz) {    // if a new element is added, make sure the count starts at 1
				lengths[ chromNames[iPos] ] = 1;
			}
		} else {  // the sites are divergent
			stringstream siteInfo;
			siteInfo << chromNames[iPos] << "\t";
			siteInfo << positions[iPos] << "\t";
			siteInfo << primary << "\t" << aligned << "\t";
			siteInfo << same << "\t";
			if ( isupper(primary) && isupper(aligned) ) {
				siteInfo << "1";
			} else {
				siteInfo << "0";
			}
			sites.push_back( siteInfo.str() );

			uint64_t tmpSz = lengths.size(); // save the pre-insertion size
			lengths[ chromNames[iPos] ]++;
			if (lengths.size() > tmpSz) {    // if a new element is added, make sure the count starts at 1
				lengths[ chromNames[iPos] ] = 1;
			}
		}
	}
}

void ParseAXT::getOutgroupState(const string &chromName, const uint64_t &position, string &site){
	if (chromName == foundChr_) { // this chromosome already completed; site unavailable
		site = "N00";
		return;
	}
	char primary;
	char aligned;
	uint16_t same;
	getSiteStates_(chromName, position, primary, aligned, same); // will search the .axt records
	if ( (aligned == '-') || (aligned == 'n') || (aligned == 'N') ) {
		site = "N0";
		if (same) {
			site += "1";
		} else {
			site += "0";
		}
	} else {
		site  = aligned;
		if (isupper(aligned)) {
			site += "1";
		} else {
			site += "0";
		}

		if (same) {
			site += "1";
		} else {
			site += "0";
		}
	}
}

void ParseAXT::getNextRecord_(){
	string curLine("");
	while(getline(axtFile_, curLine)){
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
		if ( (fields[1][0] != 'X') || (fields[1][0] != '2') || (fields[1][0] != '3') || (fields[1][0] != '4') ){
			string wrongThing = "Wrong chromosome field: " + fields[1];
			throw wrongThing;
		} else {
			fields[1] = "chr" + fields[1];
		}
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
	} else if (primaryEnd_ < primaryStart_) {  // strictly less because there are cases where only one nucleotide is in a record
		stringstream wrongThing;
		wrongThing << "Position of the end of primary sequence (";
		wrongThing << fields[3];
		wrongThing << ") not greater than the position of the start: ";
		wrongThing << primaryStart_;
		wrongThing << " (if the last number is huge, the original number was likely negative)";
		throw wrongThing.str();
	}

	if ( (fields[4][0] != 'c') || (fields[4][1] != 'h') || (fields[4][2] != 'r') ) { // do not have "chr" at the beginning of the chromosome field
		if ( (fields[4][0] != 'X') || (fields[4][0] != '2') || (fields[4][0] != '3') || (fields[4][0] != '4') ){
			string wrongThing = "Wrong aligned chromosome field: " + fields[4];
			throw wrongThing;
		} else {
			fields[4] = "chr" + fields[4];
		}
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

void ParseAXT::getSiteStates_(const string &chromosome, const uint64_t &position, char &primaryState, char &alignedState, uint16_t &sameChromosome){
	bool noneFound = true;
	bool correctChrFound = false;
	while(axtFile_){
		if (chrID_ != chromosome) {
			if (correctChrFound) { // blew past the correct chromosome without finding the position (maybe not covered)
				primaryState   = '-';
				alignedState   = '-';
				sameChromosome = 0;
				noneFound      = false;
				if (foundChr_ != chromosome) {
					foundChr_ = chromosome;
				}
				return;
			}
			getNextRecord_();
			continue;
		}
		correctChrFound = true;
		if (primaryEnd_ >= position) {
			if (position < primaryStart_) {  // our position may fall into a gap between alignment chunks; will return values that will be filtered downstream
				primaryState   = '-';
				alignedState   = '-';
				sameChromosome = 0;
				noneFound      = false;
				return;
			}
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
					return;
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
		wrongThing << "Reached the end of file before finding a record for positition ";
		wrongThing << position;
		wrongThing << " on chromosome ";
		wrongThing << chromosome;
		throw wrongThing.str();
	}
}

