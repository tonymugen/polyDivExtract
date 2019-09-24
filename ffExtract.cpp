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

FFextract::FFextract(const string &fastaName, const string &logName) : header_{""}, sequence_{""}, end_{0}, chr_{""}, fbgn_{""}, delStart_{0}, delLength_{0} {
	if (fastaFile_.is_open()) {
		fastaFile_.close();
	}
	if (logFile_.is_open()) {
		logFile_.close();
	}
	try {
		fastaFile_.open(fastaName.c_str(), ios::in);
	} catch (system_error &error){
		string message = "ERROR: cannot open file " + fastaName + ": " + error.code().message();
		throw message;
	}
	try {
		logFile_.open(logName.c_str(), ios::out|ios::trunc);
	} catch(system_error &error) {
		string message = "ERROR: cannot open file " + logName + ": " + error.code().message();
		throw message;
	}
}

FFextract::~FFextract(){
	if (fastaFile_.is_open()) {
		fastaFile_.close();
	}
	if (logFile_.is_open()) {
		logFile_.close();
	}
}

void FFextract::parseHeader_(vector<uint64_t> &positions, string &chr, string &fbgn){
	positions.clear();
	stringstream hSS(header_);
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
					error += header_;
					error += "\n";
					throw error;
				}
				uint64_t start;
				uint64_t end;
				getNumbers_(field, start, end);
				if (start >= end) {
					string error("Start position is not before the end position in header ");
					error += header_;
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
							error += header_;
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
							error += header_;
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
		} else if (field.compare(0, 7, "parent=") == 0) {
			fbgn = field.substr(11, 7);
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
	string curLine;
	sequence_.clear();
	while(getline(fastaFile_, curLine)){
		if (curLine[0] == '>') {
			header_ = move(curLine);
			vector<uint64_t> curPos;
			string curChr;
			string curFBgn;
			parseHeader_(curPos, curChr, curFBgn);
			if (curChr != chr_) { // new chromosome; no need to check for overlap
				logFile_ << "Switched from chromosome " << chr_ << " to " << curChr << " at FBgn" << curFBgn << endl;
				getFFsites_(); // extracting from the previous record
				positions_ = move(curPos);
				end_       = ( ( positions_[0] < positions_.back() ) ? positions_.back() : positions_[0] );
				chr_       = move(curChr);
				fbgn_      = move(curFBgn);
				delLength_ = 0;
				continue;
			} else if ( positions_.empty() ) {
				logFile_ << "Previous record empty at FBgn" << curFBgn << endl;
				positions_ = move(curPos);
				end_       = ( ( positions_[0] < positions_.back() ) ? positions_.back() : positions_[0] );
				chr_       = move(curChr);
				fbgn_      = move(curFBgn);
				delLength_ = 0;
				continue;
			}
			if ( curPos[0] < curPos.back() ) { // current CDS not complemented
				if (curPos[0] < end_) { // there is overlap
					logFile_ << "Detected overlap between " << fbgn_ << " and " << curFBgn << endl;
					if ( positions_[0] < positions_.back() ) { // previous CDS not complemented, either
						delLength_ = 0;
						for (auto pRit = positions_.rbegin(); pRit != positions_.rend(); ++pRit) {
							if ( (*pRit) >= curPos[0] ) {
								delLength_++;
							} else {
								break;
							}
						}
						delLength_ = delLength_ + (delLength_%3);
						if ( (delLength_ < positions_.size()) && (delLength_ < curPos.size()) ) { // delLength_ can be larger than these because of the rounding to codons
							delStart_  = positions_.size() - delLength_;
							positions_.resize(positions_.size() - delLength_);
							sequence_.resize(sequence_.size() - delLength_);
							getFFsites_();
							curPos.erase(curPos.begin(), curPos.begin()+delLength_);
							positions_ = move(curPos);
							end_       = positions_.back();
							chr_       = move(curChr);
							fbgn_      = move(curFBgn);
							continue;
						} else if ( delLength_ >= positions_.size() ) {
							if ( delLength_ >= curPos.size() ) {
								logFile_ << fbgn_ << " deleted by overlapping " << curFBgn << ", which is also deleted" << endl;
								// do not extract FF sites, do not move the current fields
								delStart_ = 0;
								positions_.clear();
								sequence_.clear();
								continue;
							} else {
								logFile_ << fbgn_ << " deleted by overlapping " << curFBgn << endl;
								// do not extract FF sites
								delStart_ = curPos.size() - delLength_;
								curPos.erase(curPos.begin(), curPos.begin()+delLength_);
								positions_ = move(curPos);
								end_       = positions_.back();
								chr_       = move(curChr);
								fbgn_      = move(curFBgn);
								continue;
							}
						} else {
							logFile_ << fbgn_ << " deletes the overlapping " << curFBgn << endl;
							// extract sites, but do not move the new position info (chromosome and FBgn info is moved)
							positions_.resize(positions_.size() - delLength_);
							sequence_.resize(sequence_.size() - delLength_);
							getFFsites_();
							positions_.clear();
							delStart_ = 0;
							end_      = 0;
							chr_      = move(curChr);
							fbgn_     = move(curFBgn);
							continue;
						}
					} else { // previous CDS complemented
						delLength_ = 0;
						for (auto pIt = positions_.begin(); pIt != positions_.end(); ++pIt) {
							if ( (*pIt) >= curPos[0] ) {
								delLength_++;
							} else {
								break;
							}
						}
						delLength_ = delLength_ + (delLength_%3);
						delStart_  = 0;
						if ( (delLength_ < positions_.size()) && (delLength_ < curPos.size()) ) { // delLength_ can be larger than these because of the rounding to codons
							positions_.erase(positions_.begin(), positions_.begin() + delLength_);
							sequence_.erase(0, delLength_);
							getFFsites_();
							curPos.erase(curPos.begin(), curPos.begin()+delLength_);
							positions_ = move(curPos);
							end_       = positions_.back();
							chr_       = move(curChr);
							fbgn_      = move(curFBgn);
							continue;
						} else if ( delLength_ >= positions_.size() ) {
							if ( delLength_ >= curPos.size() ) {
								logFile_ << fbgn_ << " deleted by overlapping " << curFBgn << ", which is also deleted" << endl;
								// do not extract FF sites, do not move the current fields
								positions_.clear();
								sequence_.clear();
								continue;
							} else {
								logFile_ << fbgn_ << " deleted by overlapping " << curFBgn << endl;
								// do not extract FF sites
								curPos.erase(curPos.begin(), curPos.begin()+delLength_);
								positions_ = move(curPos);
								end_       = positions_.back();
								chr_       = move(curChr);
								fbgn_      = move(curFBgn);
								continue;
							}
						} else {
							logFile_ << fbgn_ << " deletes the overlapping " << curFBgn << endl;
							// extract sites, but do not move the new position info (chromosome and FBgn info is moved)
							positions_.erase(positions_.begin(), positions_.begin() + delLength_);
							sequence_.erase(0, delLength_);
							getFFsites_();
							positions_.clear();
							end_  = 0;
							chr_  = move(curChr);
							fbgn_ = move(curFBgn);
							continue;
						}
					}
				} else { // no overlap
					getFFsites_();
					positions_ = move(curPos);
					end_       = positions_.back();
					chr_       = move(curChr);
					fbgn_      = move(curFBgn);
					delLength_ = 0;
					continue;
				}
			} else { // current CDS complemented
				if (curPos.back() < end_) { // there is overlap
					logFile_ << "Detected overlap between " << fbgn_ << " and " << curFBgn << endl;
					if ( positions_[0] < positions_.back() ) { // previous CDS not complemented, either
						delLength_ = 0;
						for (auto pRit = positions_.rbegin(); pRit != positions_.rend(); ++pRit) {
							if ( (*pRit) >= curPos.back() ) {
								delLength_++;
							} else {
								break;
							}
						}
						delLength_ = delLength_ + (delLength_%3);
						if ( (delLength_ < positions_.size()) && (delLength_ < curPos.size()) ) { // delLength_ can be larger than these because of the rounding to codons
							delStart_  = positions_.size() - delLength_;
							positions_.resize(positions_.size() - delLength_);
							sequence_.resize(sequence_.size() - delLength_);
							getFFsites_();
							curPos.resize(curPos.size() - delLength_);
							positions_ = move(curPos);
							end_       = positions_[0];
							chr_       = move(curChr);
							fbgn_      = move(curFBgn);
							continue;
						} else if ( delLength_ >= positions_.size() ) {
							if ( delLength_ >= curPos.size() ) {
								logFile_ << fbgn_ << " deleted by overlapping " << curFBgn << ", which is also deleted" << endl;
								// do not extract FF sites, do not move the current fields
								delStart_ = 0;
								positions_.clear();
								sequence_.clear();
								continue;
							} else {
								logFile_ << fbgn_ << " deleted by overlapping " << curFBgn << endl;
								// do not extract FF sites
								curPos.resize(curPos.size() - delLength_);
								delStart_  = curPos.size() - delLength_;
								positions_ = move(curPos);
								end_       = positions_[0];
								chr_       = move(curChr);
								fbgn_      = move(curFBgn);
								continue;
							}
						} else {
							logFile_ << fbgn_ << " deletes the overlapping " << curFBgn << endl;
							// extract sites, but do not move the new position info (chromosome and FBgn info is moved)
							positions_.resize(positions_.size() - delLength_);
							sequence_.resize(sequence_.size() - delLength_);
							getFFsites_();
							positions_.clear();
							delStart_ = 0;
							end_      = 0;
							chr_      = move(curChr);
							fbgn_     = move(curFBgn);
							continue;
						}
					} else { // previous CDS complemented
						delLength_ = 0;
						for (auto pIt = positions_.begin(); pIt != positions_.end(); ++pIt) {
							if ( (*pIt) >= curPos.back() ) {
								delLength_++;
							} else {
								break;
							}
						}
						delLength_ = delLength_ + (delLength_%3);
						delStart_  = 0;
						if ( (delLength_ < positions_.size()) && (delLength_ < curPos.size()) ) { // delLength_ can be larger than these because of the rounding to codons
							positions_.erase(positions_.begin(), positions_.begin() + delLength_);
							sequence_.erase(0, delLength_);
							getFFsites_();
							curPos.resize(curPos.size() - delLength_);
							positions_ = move(curPos);
							end_       = positions_[0];
							chr_       = move(curChr);
							fbgn_      = move(curFBgn);
							continue;
						} else if ( delLength_ >= positions_.size() ) {
							if ( delLength_ >= curPos.size() ) {
								logFile_ << fbgn_ << " deleted by overlapping " << curFBgn << ", which is also deleted" << endl;
								// do not extract FF sites, do not move the current fields
								positions_.clear();
								sequence_.clear();
								continue;
							} else {
								logFile_ << fbgn_ << " deleted by overlapping " << curFBgn << endl;
								// do not extract FF sites
								curPos.resize(curPos.size() - delLength_);
								positions_ = move(curPos);
								end_       = positions_[0];
								chr_       = move(curChr);
								fbgn_      = move(curFBgn);
								continue;
							}
						} else {
							logFile_ << fbgn_ << " deletes the overlapping " << curFBgn << endl;
							// extract sites, but do not move the new position info (chromosome and FBgn info is moved)
							positions_.erase(positions_.begin(), positions_.begin() + delLength_);
							sequence_.erase(0, delLength_);
							getFFsites_();
							positions_.clear();
							end_  = 0;
							chr_  = move(curChr);
							fbgn_ = move(curFBgn);
							continue;
						}
					}
				} else { // no overlap
					getFFsites_();
					positions_ = move(curPos);
					end_       = positions_[0];
					chr_       = move(curChr);
					fbgn_      = move(curFBgn);
					delLength_ = 0;
					continue;
				}
			}
		} else {
			if ( positions_.empty() ) {
				break; // sequence_ would already have been emptied
			} else if (delLength_) {
				curLine.erase(delStart_, delLength_);
			}
			sequence_ = move(curLine);
			break;
		}
	}
}

void FFextract::getFFsites_(){
	for (size_t i = 0; i < sequence_.size(); i += 3) {
		string codon = sequence_.substr(i, 3);
		if (codon[1] == 'A') { // no codons with A in second position have four-fold sites
			continue;
		} else if (codon[1] == 'T') {
			if ( (codon[0] == 'C') || (codon[0] == 'G') ) {
				stringstream rSS(ios::out);
				rSS << chr_ << "\t" << fbgn_ << "\t" << positions_[i+3] << "\n";
				ffSites_.push_back( rSS.str() );
			}
		} else if (codon[1] == 'C') { // all codons with C in second position are four-fold
			stringstream rSS(ios::out);
			rSS << chr_ << "\t" << fbgn_ << "\t" << positions_[i+3] << "\n";
			ffSites_.push_back( rSS.str() );
		} else { // G
			if ( (codon[0] == 'C') || (codon[0] == 'G') ) {
				stringstream rSS(ios::out);
				rSS << chr_ << "\t" << fbgn_ << "\t" << positions_[i+3] << "\n";
				ffSites_.push_back( rSS.str() );
			}
		}
	}
}


