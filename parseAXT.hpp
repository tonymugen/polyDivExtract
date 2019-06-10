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
 * Class definitions and interface documentation for parsing of .axt alignement files to extract ancestral nucleotide and divergence information
 *
 */

#ifndef parseAXT_hpp
#define parseAXT_hpp

#include <fstream>
#include <string>
#include <vector>

using std::fstream;
using std::string;
using std::vector;

namespace BayesicSpace {

	/** \brief .axt alignment parsing class
	 *
	 * Exatracts features from an .axt alignement file.
	 *
	 */
	class ParseAXT {
		public:
			/** \brief Default constructor */
			ParseAXT() : chrID_{""}, sameChr_{0}, primaryStart_{0}, primaryEnd_{0}, alignedStart_{0}, alignedEnd_{0}, primarySeq_{""}, alignSeq_{""} { axtFile_.exceptions(fstream::badbit); };
			/** \brief File name constructor
			 *
			 * Initializes the file stream and loads first AXT record.
			 *
			 * \param[in] fileName file name
			 */
			ParseAXT(const string &fileName);

			/** \brief Destructor */
			~ParseAXT(){ if(axtFile_.is_open()) axtFile_.close(); };

			/// Copy constructor
			ParseAXT(const ParseAXT &in) = delete;
			/// Move constructor
			ParseAXT(ParseAXT &&in) = delete;
			/// Copy assignment
			ParseAXT &operator=(const ParseAXT &in) = delete;
			/// Move assignment
			ParseAXT &operator=(ParseAXT &&in) = delete;

			/** \brief Record meta data
			 *
			 * Returns a string with space-delimited metadata for the curent record:
			 *
			 * - primary chromosome
			 * - is the aligned chromosome the same (0/1)?
			 * - primary start
			 * - primary end
			 * - aligned start
			 * - aligned end
			 *
			 *   \return string with metadata
			 */
			string getMetaData();
			/** \brief Get primary sequence
			 *
			 * \return string with the primary sequence
			 */
			string getPrimarySeq() {return primarySeq_; };
			/** \brief Get aligned sequence
			 *
			 * \return string with the aligned sequence
			 */
			string getAlignedSeq() {return alignSeq_; };
			/** \brief Get list of divergent sites from a range
			 *
			 * Get a list of divergent sites from a range of positions on a chromosome. Sites that are not covered or align to gaps are not counted in computing the overall length.
			 * The vector of sites is appended by the function, so any exsiting information will be preserved.
			 * The site description is in a tab-delimited string with the following fields:
			 *
			 * - chromosome name
			 * - position
			 * - primary nucleotide
			 * - aligned nucleotide
			 * - whether the aligned nucleotide is on the same chromosome
			 * - whether both nucleotides are in upper case (indicating high quality base calls)
			 *
			 * _NOTE_: The range must be confined to a single chromosome.
			 *
			 * \param[in] chromName chromosome name
			 * \param[in] start start position of the target range
			 * \param[in] end end position of the target range
			 * \param[out] sites vector of divergent site information (appended after execution)
			 * \param[out] length length not counting sites that are missing or align to gaps
			 *
			 */
			void getDivergedSites(const string &chromName, const uint64_t &start, const uint64_t &end, vector<string> &sites, uint64_t &length);
		private:
			/// The file stream
			fstream axtFile_;

			// variables for the current record
			/// Primary chromosome
			string chrID_;
			/// Is the aligned chromosome the same (1 for yes, 0 for no)?
			uint16_t sameChr_;
			/// Primary start position
			uint64_t primaryStart_;
			/// Primary end position
			uint64_t primaryEnd_;
			/// Aligned start posistion
			uint64_t alignedStart_;
			/// Aligned end position
			uint64_t alignedEnd_;
			/// Current record's primary sequence
			string primarySeq_;
			/// Current record's aligning sequence
			string alignSeq_;

			/** \brief Get next record */
			void getNextRecord_();
			/** \brief Extracts the nucleotides at a given position
			 *
			 * The query position references the primary sequence
			 *
			 * \param[in] chromosome primary chromosome
			 * \param[in] position site position in the primary sequence
			 * \param[out] primaryState the primary nucleotide at the query position
			 * \param[out] alignedState the aligned nucleotide at the query position
			 * \param[out] sameChromosome is the aligned sequence on the same chromosome as primary? (0: no, 1: yes)
			 *
			 */
			void getSiteStates_(const string &chromosome, const uint64_t &position, char &primaryState, char &alignedState, uint16_t &sameChromosome);
	};
}
#endif /* parseAXT_hpp */

