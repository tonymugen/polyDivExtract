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

/// Parse .vcf files
/** \file
 * \author Anthony J. Greenberg
 * \copyright Copyright (c) 2019 Anthony J. Greenberg
 * \version 0.1
 *
 * Class implementation for parsing of .vcf files to extract ancestral nucleotide and polymorphism information. This is specific to the dosage compensation evolution project.
 *
 */
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <cstdlib>
#include <cctype>
#include <system_error>

#include "parseVCF.hpp"
#include "parseAXT.hpp"

using std::fstream;
using std::ofstream;
using std::stringstream;
using std::string;
using std::vector;
using std::endl;
using std::system_error;
using std::ios;
using std::bad_alloc;

using namespace BayesicSpace;


ParseVCF::ParseVCF(const string &vcfFileName, const string &axtFileName) : varPos_{0}, refID_{'\0'}, altID_{'\0'}, ancState_{'u'}, outQual_{0}, sameChr_{0}, numMissing_{0}, numCalled_{0}, refAC_{0}, refMLAC_{0}, refAF_{0.0}, refMLAF_{0.0}, quality_{0.0}, chrID_{""}, foundChr_{""}, fullRecord_{""} {
	if( vcfFile_.is_open() ){
		vcfFile_.close();
	}

	try {
		vcfFile_.open(vcfFileName.c_str(), ios::in);
	} catch(system_error &error) {
		string message = "ERROR: cannot open file " + vcfFileName + " to read: " + error.code().message();
		throw message;
	}

	axtObj_ = ParseAXT(axtFileName);

	while(getline(vcfFile_, fullRecord_)){
		if (fullRecord_[0] == '#') {
			continue;
		} else if (fullRecord_ == "") {
			continue;
		} else {
			break;
		}
	}
	if (fullRecord_ == ""){
		throw string("No non-empty non-comment lines in file ") + vcfFileName;
	}

}

void ParseVCF::getPolySites(const string &chromName, const uint64_t &start, const uint64_t &end, vector<string> &sites){
	if (start >= end) {
		stringstream wrongThing;
		wrongThing << "ERROR: start position (";
		wrongThing << start;
		wrongThing << ") must come before the end postion (";
		wrongThing << end;
		wrongThing << ") in getPolySites()";
		throw wrongThing.str();
	}
	if (chromName == foundChr_) {
		return;
	}
	bool foundChrom = false; // keep track if the target chromosome was found in the search; needed to test if we looked though the whole thing without finding our site(s)

	// process first record (already loaded at construction and guaranteed non-empty)
	if ( (chromName == "chrX") || (chromName == "chr4") ) {
		if (chromName[3] == fullRecord_[0]) {
			foundChrom = true;
			stringstream recSS(fullRecord_);
			string field;
			recSS >> field;
			recSS >> field;
			uint64_t curPos = strtoul(field.c_str(), NULL, 0);
			if ( (curPos >= start) && (curPos <= end) ) {
				parseCurrentRecord_();
				sites.push_back( exportCurRecord_() );
			}
		} else if (foundChrom) {
			return;
		}
	} else {
		if (chromName.compare(3, 2, fullRecord_, 0, 2) == 0) {
			foundChrom = true;
			stringstream recSS(fullRecord_);
			string field;
			recSS >> field;
			recSS >> field;
			uint64_t curPos = strtoul(field.c_str(), NULL, 0);
			if ( (curPos >= start) && (curPos <= end) ) {
				parseCurrentRecord_();
				sites.push_back( exportCurRecord_() );
			}
		} else if (foundChrom) {
			return;
		}
	}
	while(getline(vcfFile_, fullRecord_)){
		if (fullRecord_.size() == 0) {
			continue;
		}
		if ( (chromName == "chrX") || (chromName == "chr4") ) {
			if (chromName[3] == fullRecord_[0]) {
				foundChrom = true;
				stringstream recSS(fullRecord_);
				string field;
				recSS >> field;
				recSS >> field;
				uint64_t curPos = strtoul(field.c_str(), NULL, 0);
				if ( (curPos >= start) && (curPos <= end) ) {
					parseCurrentRecord_();
					sites.push_back( exportCurRecord_() );
				}
			} else if (foundChrom) {
				return;
			}
		} else {
			if (chromName.compare(3, 2, fullRecord_, 0, 2) == 0) {
				foundChrom = true;
				stringstream recSS(fullRecord_);
				string field;
				recSS >> field;
				recSS >> field;
				uint64_t curPos = strtoul(field.c_str(), NULL, 0);
				if ( (curPos >= start) && (curPos <= end) ) {
					parseCurrentRecord_();
					sites.push_back( exportCurRecord_() );
				}
			} else if (foundChrom) {
				return;
			}
		}
	}
}

//void ParseVCF::getPolySites(const vector<string> &chromNames, const vector<uint64_t> &positions, vector<string> &sites){
//
//}

void ParseVCF::parseCurrentRecord_(){
	// we have a non-empty line, presumably a VCF record
	stringstream metaSS(fullRecord_);
	vector<string> fields;
	string value;
	while(metaSS >> value){ // makes sure we don't have anything extra at the end
		fields.push_back(value);
	}

	varPos_  = strtoul(fields[1].c_str(), NULL, 0);
	refID_   = fields[3][0];
	altID_   = fields[4][0];
	chrID_   = "chr" + fields[0];
	quality_ = strtod(fields[5].c_str(), NULL);

	// count missing data; would be a bit faster to start from fields[9], but the overhead should be small at the expense of the better for loop
	numMissing_ = 0;
	for (auto &s : fields) {
		if (s == "./.") {
			numMissing_++;
		}
	}
	// parse the INFO field
	stringstream fSS(fields[7]);
	string info;
	while( getline(fSS, info, ';') ){
		if (info.compare(0, 2, "AC") == 0) {
			refAC_ = strtoul(info.c_str()+3, NULL, 0);
		} else if (info.compare(0, 2, "AF") == 0) {
			refAF_ = strtod(info.c_str()+3, NULL);
		} else if (info.compare(0, 2, "AN") == 0) {
			numCalled_ = strtoul(info.c_str()+3, NULL, 0);
		} else if (info.compare(0, 5, "MLEAC") == 0) {
			refMLAC_ = strtoul(info.c_str()+6, NULL, 0);
		} else if (info.compare(0, 5, "MLEAF") == 0) {
			refMLAF_ = strtod(info.c_str()+6, NULL);
		}
	}
	// Now find the ancestral state if we can
	string outInfo;
	axtObj_.getOutgroupState(chrID_, varPos_, outInfo);
	if (outInfo[0] == 'N') {
		ancState_ = 'u';
		sameChr_  = 0;
		outQual_  = 0;
	} else {
		ancState_ = (outInfo[0] == refID_ ? 'r' : 'a');
		outQual_  = (outInfo[1] == '1' ? 1 : 0);
		sameChr_  = (outInfo[2] == '1' ? 1 : 0);
	}
}

string ParseVCF::exportCurRecord_(){
	stringstream siteInfo;
	siteInfo << chrID_ << "\t";
	siteInfo << varPos_ << "\t";
	siteInfo << refID_ << "\t";
	siteInfo << altID_ << "\t";
	siteInfo << ancState_ << "\t";
	if (ancState_ == 'a') {
		siteInfo << numCalled_ - refAC_ << "\t";
		siteInfo << numCalled_ - refMLAC_ << "\t";
		siteInfo << 1.0 - refAF_ << "\t";
		siteInfo << 1.0 - refMLAF_ << "\t";
	} else {
		siteInfo << refAC_ << "\t";
		siteInfo << refMLAC_ << "\t";
		siteInfo << refAF_ << "\t";
		siteInfo << refMLAF_ << "\t";
	}
	siteInfo << numMissing_ << "\t";
	siteInfo << sameChr_ << "\t";
	siteInfo << outQual_ << "\t";
	siteInfo << quality_;

	return siteInfo.str();
}

