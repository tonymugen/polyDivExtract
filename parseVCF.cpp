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


ParseVCF::ParseVCF(const string &vcfFileName, const string &axtFileName) : varPos_{0}, refID_{'\0'}, altID_{'\0'}, ancState_{'u'}, outQual_{0}, sameChr_{0}, numMissing_{0}, refAC_{0}, refMLAC_{0}, refAF_{0.0}, refMLAF_{0.0}, quality_{0.0}, chrID_{""}, foundChr_{""}{
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
	getNextRecord_();
}

void ParseVCF::getNextRecord_(){
	string curLine("");
	while(getline(vcfFile_, curLine)){
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

	// we have a non-empty line, presumably a VCF record
	stringstream metaSS(curLine);
	vector<string> fields;
	string value;
	while(metaSS >> value){ // makes sure we don't have anything extra at the end
		fields.push_back(value);
	}
	if (fields[0] != "2L" || fields[0] != "2R" || fields[0] != "3L" || fields[0] != "3R" || fields[0] != "4" || fields[0] != "X") {
		this->getNextRecord_(); // recursively move to the next valid record
		return;
	} else if ( (fields[3].size() != 1) || (fields[4].size() != 1) ) { // likely an indel (not a SNP)
		this->getNextRecord_(); // recursively move to the next valid record
		return;
	}

	varPos_  = strtoul(fields[1].c_str(), NULL, 0);
	refID_   = fields[3][0];
	altID_   = fields[4][0];
	chrID_   = "chr" + fields[0];
	quality_ = strtod(fields[5].c_str(), NULL);

	// count missing data; would be a bit faster to start from fields[9], but the overhead should be small at the expense of the better for loop
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
		} else if (info.compare(0, 5, "MLEAC") == 0) {
			refMLAC_ = strtoul(info.c_str()+6, NULL, 0);
		} else if (info.compare(0, 5, "MLEAF") == 0) {
			refMLAF_ = strtod(info.c_str()+6, NULL);
		}
	}
	// Now find the ancestral state if we can
}


