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

/// Extract four-fold sites
/** \file
 * \author Anthony J. Greenberg
 * \copyright Copyright (c) 2019 Anthony J. Greenberg
 * \version 0.1
 *
 * Class implementation for extracting four-fold synonymous site positions from FASTA files.
 *
 */

#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cctype>
#include <system_error>

#include <iostream>

#include "ffExtract.hpp"

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

FFextract::FFextract(const string &fastaName) : nextHeader_{""}, sequence_{""}, end_{0}, chr_{""}, delStart_{0}, delLength_{0} {
	if (fastaFile_.is_open()) {
		fastaFile_.close();
	}
	try {
		fastaFile_.open(fastaName.c_str(), ios::in);
	} catch (system_error &error){
		string message = "ERROR: cannot open file " + fastaName + ": " + error.code().message();
		throw message;
	}

	getline(fastaFile_, nextHeader_);
	parseHeader_(positions_, chr_);
	end_ = ( (positions_[0] < positions_.back()) ? positions_.back() : positions_[0] );
	//getline(fastaFile_, sequence_);
	if (sequence_.size() != positions_.size()) {
		std::cerr << "sequence (" << sequence_.size() << ") and position (" << positions_.size() << ") lengths not equal" << endl;
		exit(1);
	}
	fstream tFS;
	tFS.open("tstOut.tsv", ios::out|ios::trunc);
	tFS << chr_ << " " << end_ << endl;
	for (size_t i = 0; i < sequence_.size(); ++i) {
		tFS << sequence_[i] << " " << positions_[i] << endl;
	}
	tFS.close();
	//getNextRecord_();
}

FFextract::~FFextract(){
	if (fastaFile_.is_open()) {
		fastaFile_.close();
	}
}

void FFextract::parseHeader_(vector<uint64_t> &positions, string &chr){
	positions.clear();
	stringstream hSS(nextHeader_);
	string field;
	while( getline(hSS, field, ' ') ){
		if (field.compare(0, 3, "loc") == 0) {
			bool complemented = false; // is the current record complemented?
			field.erase(0, 4);
			if (field.compare(0, 4, "Scf_") == 0) {
				field.erase(0, 4);
			}
			if ( (field[0] == 'X') || (field[0] == '4') ) {
				chr = field[0];
				field.erase(0, 2);
			} else if ( (field[0] == '3') || (field[0] == '2') ) {
				chr.assign(field, 0, 2);
				field.erase(0,3);
			} else {
				string error("ERROR: unkown chromosome ");
				size_t col = field.find_first_of(':');
				if (col > 0) {
					col--;
				}
				error += field.substr(0, col);
				throw error;
			}
			field.erase(field.end()-1);
			if (field[0] == 'c') { // the only way this occurs is when complement() is specified
				complemented = true;
				field.erase(0, 11);
				field.erase(field.end()-1);
			}
			if (isdigit(field[0])) { // only one exon, no complement
				if (field.size() <= 2) {
					string error("Cannot parse postion range in header\n");
					error += nextHeader_;
					error += "\n";
					throw error;
				}
				uint64_t start;
				uint64_t end;
				getNumbers_(field, start, end);
				if (start >= end) {
					string error("Start position is not before the end position in header ");
					error += nextHeader_;
					error += "\n";
					throw error;
				}
				if (complemented) {
					while (end != start){
						positions.push_back(end);
						end--;
					}
					positions.push_back(start); // necessary to top off
				} else {
					while (start != end){
						positions.push_back(start);
						start++;
					}
					positions.push_back(end); // necessary to top off
				}
			} else if (field[0] == 'j') { //  there is a join
				field.erase(0, 5);
				field.erase(field.end()-1);
				stringstream fSS(field);
				string curRange;
				vector<string> ranges;
				while ( getline(fSS, curRange, ',') ){
					ranges.push_back(curRange);
				}
				if (complemented) {
					for (auto rIt = ranges.rbegin(); rIt != ranges.rend(); ++rIt) {
						uint64_t start;
						uint64_t end;
						getNumbers_(*rIt, start, end);
						if (start >= end) {
							string error("Start position is not before the end position in header ");
							error += nextHeader_;
							error += "\n";
							throw error;
						}
						while (end != start){
							positions.push_back(end);
							end--;
						}
						positions.push_back(start); // necessary to top off
					}
				} else {
					for (auto it = ranges.begin(); it != ranges.end(); ++it) {
						uint64_t start;
						uint64_t end;
						getNumbers_(*it, start, end);
						if (start >= end) {
							string error("Start position is not before the end position in header ");
							error += nextHeader_;
							error += "\n";
							throw error;
						}
						while (start != end){
							positions.push_back(start);
							start++;
						}
						positions.push_back(end); // necessary to top off
					}
				}
			} else {
				string error("Unknown value in position list of field ");
				error += field;
				throw error;
			}
			break;
		}
	}
}

void FFextract::getNumbers_(const string &range, uint64_t &start, uint64_t &end){
	size_t pos = range.find_first_of('.'); // assuming that the range was checked for non-empty number fields
	start      = strtoul(range.substr(0, pos).c_str(), NULL, 0);
	pos        = range.find_last_of('.') + 1;
	end        = strtoul(range.substr(pos).c_str(), NULL, 0);
}

void FFextract::getNextRecord_(){
	// first header already read on construction
	string curLine;
	sequence_.clear();
	while(getline(fastaFile_, curLine)){
		if (curLine[0] == '>') {
			vector<uint64_t> tmpPos;
			string tmpChr;
			nextHeader_ = move(curLine);
			parseHeader_(tmpPos, tmpChr);
			if (chr_ != tmpChr) {
				break;
			}
			const uint64_t sNew = ( ( tmpPos[0] < tmpPos.back() ) ? tmpPos[0] : tmpPos.back() );
			if (sNew < end_){ // one nucleotide overlap is OK
				// test if the CDS with the new header is completely within the previous CDS
				const uint64_t eNew = ( ( tmpPos[0] >= tmpPos.back() ) ? tmpPos[0] : tmpPos.back() );
				if ( eNew <= (end_ + 3) ) { // no need to bother with a single codon, it's a stop
					if ( positions_[0] < positions_.back() ) {
						size_t iEnd = positions_.size() - 1;
						while ( (positions_[iEnd] > eNew) && (iEnd != 0) ){
							iEnd--;
						}
						size_t iBeg = iEnd;
						while ( (positions_[iBeg] > sNew) && (iBeg != 0) ){
							iBeg--;
						}
						if (iBeg != iEnd) {
							sequence_.erase(sequence_.begin() + iBeg - (iBeg%3), sequence_.begin() + iEnd + 3 - (iEnd%3));
							positions_.erase(positions_.begin() + iBeg - (iBeg%3), positions_.begin() + iEnd + 3 - (iEnd%3));
						}
						getline(fastaFile_, curLine); // read and discard the sequence for this record; it is completely within the previous locus. Discard even if within an intron.
						return;
					} else {
						size_t iBeg = 0;
						while ( (positions_[iBeg] > eNew) && (iBeg < positions_.size()) ){
							iBeg++;
						}
						size_t iEnd = iBeg;
						while ( (positions_[iEnd] > sNew) && (iEnd < positions_.size()) ){
							iEnd++;
						}
						if (iBeg != iEnd) {
							sequence_.erase(sequence_.begin() + iBeg - (iBeg%3), sequence_.begin() + iEnd + 3 - (iEnd%3));
							positions_.erase(positions_.begin() + iBeg - (iBeg%3), positions_.begin() + iEnd + 3 - (iEnd%3));
						}
						getline(fastaFile_, curLine); // read and discard the sequence for this record; it is completely within the previous locus. Discard even if within an intron.
						return;
					}
				} else if ( positions_[0] < positions_.back() ) { // there is overlap, but the next CDS is not within the new one. Have to truncate both
					uint64_t i = positions_.size() - 1;
					while (positions_[i] >= sNew){
						if (i == 0) {
							break;
						}
						i--;
					}

				}
			}
			break;
		} else {
			sequence_ = move(curLine);
		}
	}
}


