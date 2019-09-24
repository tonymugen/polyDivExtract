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
 * Class definitions and interface documentation for extracting four-fold synonymous site positions from FASTA files.
 *
 */

#ifndef ffExtract_hpp
#define ffExtract_hpp
#include <fstream>
#include <string>
#include <vector>

using std::fstream;
using std::string;
using std::vector;

namespace BayesicSpace {
	/** \brief Four-fold synonymous site extraction
	 *
	 * The class reads a FASTA file with coding sequences and extracts four-fold synonymous sites.
	 * The algorithm assumes that the records are sorted by chromosome position of the start site. This can be achieved by running the enclosed `fastaSort` program.
	 * It is also assumed that there are no CDS that completely within other genes. The `fastaSort` program eliminates such cases.
	 * We also assume that the sequence portions of FASTA records are all on one line. This is how `fastaSort` outputs them.
	 * The chromosome names are for Drosophila (2L, 2R, 3L, 3R, 4, X). The names may be preceded by the Scf_ prefix, which is used in the D. simulans genome.
	 * Output chromosome names are the Drosophila set preceded by `chr`.
	 *
	 */
	class FFextract {
	public:
		/** \brief Default constructor */
		FFextract() : header_{""}, sequence_{""}, end_{0}, chr_{""}, fbgn_{""}, delStart_{0}, delLength_{0} { fastaFile_.exceptions(fstream::badbit); logFile_.exceptions(fstream::badbit); };
		/** \brief Constructor
		 *
		 * \param[in] fastaName name of the FASTA file
		 * \param[in] logName name of the log file
		 */
		FFextract(const string &fastaName, const string &logName);

		/** \brief Destructor */
		~FFextract();

		/** \brief Copy constructor (deleted) */
		FFextract(const FFextract &in) = delete;
		/** \brief Move constructor
		 *
		 * \param[in] in the object to be moved
		 */
		FFextract(FFextract &&in) : fastaFile_{move(in.fastaFile_)}, logFile_{move(in.logFile_)}, header_{move(in.header_)}, sequence_{move(in.sequence_)}, positions_{move(in.positions_)}, end_{in.end_}, chr_{move(in.chr_)}, fbgn_{move(in.fbgn_)}, delStart_{in.delStart_}, delLength_{in.delLength_} {};
		/** \brief Extract four-fold sites from the current record
		 *
		 * The vector of positions is appended.
		 *
		 * \param[out] positionList vector of four-fold site positions
		 */
		void extractFFsites(vector<uint64_t> &positionList);
	private:
		/** \brief FASTA file to be parsed */
		fstream fastaFile_;
		/** \brief Log file */
		fstream logFile_;
		/** \brief Next FASTA header
		 *
		 * Because we need to scan the sequence until we hit the next header, this variable contains the header for the next FASTA record (if any).
		 */
		string header_;
		/** \brief Current FASTA sequence */
		string sequence_;
		/** \brief Position vector
		 *
		 * The same length as `sequence_`, each element is the genome position of the nucleotide in the sequence.
		 */
		vector <uint64_t> positions_;
		/** \brief Last CDS position
		 *
		 * Not necessarily contained in `positions_` because of possible truncation.
		 */
		uint64_t end_;
		/** \brief Chromosome */
		string chr_;
		/** \brief FBgn number */
		string fbgn_;
		/** \brief Start position of sequence truncation
		 *
		 * Where to start truncation of the newly-read sequence portion of the FASTA file. Saved from parsing `nextHeader_`. Is needed only if there is overlap with the previous sequence.
		 */
		size_t delStart_;
		/** \brief Sequence truncation length
		 *
		 * Length of the new sequence region that must be deleted. Is not 0 only if there is overlap with the previous sequence.
		 */
		size_t delLength_;
		/** \brief Vector of four-fold site records */
		vector<string> ffSites_;

		/** \brief Parse the FASTA header
		 *
		 *
		 * \param[out] positions vector of chromosome positions of the sites in the sequence; any contents are replaced
		 * \param[out] chr chromosome name
		 * \param[out] fbgn FBgn number
		 */
		void parseHeader_(vector<uint64_t> &positions, string &chr, string &fbgn);
		/** \brief Parse a string to range of numbers
		 *
		 * String of a STARTPOS..ENDPOS type is parsed and the start and end position returned. The string must be validated before calling.
		 *
		 * \param[in] range the string with the range
		 */
		void getNumbers_(const string &range, uint64_t &start, uint64_t &end);
		/** \brief Get next FASTA record */
		void getNextRecord_();
		/** \brief Identifies four-fold sites in the current record
		 *
		 * The results are appended to `ffSites_`.
		 */
		void getFFsites_();
	};
}
#endif /* ffExtract.hpp */

