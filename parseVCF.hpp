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
 * Class definitions and interface documentation for parsing of .vcf files to extract ancestral nucleotide and polymorphism information. This is specific to the dosage compensation evolution project.
 *
 */

#ifndef parseVCF_hpp
#define parseVCF_hpp

#include <fstream>
#include <string>
#include <vector>

#include "parseAXT.hpp"

using std::fstream;
using std::string;
using std::vector;

namespace BayesicSpace {
	/** \brief VCF file parsing class
	 *
	 * Extracts information from a VCF file by position. Only SNPs are considered. The parsing is for the specific VCF files with fields defined in the dosage compensation project, may not be generally applicable.
	 *
	 */
	class ParseVCF {
		public:
			/** \brief Default constructor */
			ParseVCF() : varPos_{0}, refID_{'\0'}, altID_{'\0'}, ancState_{'u'}, outQual_{0}, sameChr_{0}, numMissing_{0}, numCalled_{0}, refAC_{0}, refMLAC_{0}, refAF_{0.0}, refMLAF_{0.0}, quality_{0.0}, chrID_{""}, foundChr_{""}, fullRecord_{""} { vcfFile_.exceptions(fstream::badbit); };
			/** \brief Constructor with file names
			 *
			 * Opens the VCF file and the corresponding .axt alignment file for ancestral state tracking.
			 *
			 * \param[in] vcfFileName name of the VCF file
			 * \param[in] axtFileName name of the .axt file
			 *
			 */
			ParseVCF(const string &vcfFileName, const string &axtFileName);
			
			/** \brief Destructor */
			~ParseVCF() { if(vcfFile_.is_open()) vcfFile_.close(); };
			/// Copy constructor
			ParseVCF(const ParseVCF &in) = delete;
			/// Move constructor
			ParseVCF(ParseVCF &&in) = delete;
			/// Copy assignment
			ParseVCF &operator=(const ParseVCF &in) = delete;
			/// Move assignment
			ParseAXT &operator=(ParseVCF &&in) = delete;

			/** \brief Get list of polymorphic sites from a range
			 *
			 * Get a list of polymorphic sites from a range of positions on a chromosome.
			 * The vector of sites is appended by the function, so any existing information will be preserved.
			 * The site description is in a tab-delimited string with the following fields:
			 *
			 * - chromosome name
			 * - position
			 * - reference nucleotide
			 * - alternative nucleotide
			 * - which nucleotide is ancestral ('r' for reference, 'a' alternative, 'u' unknown)
			 * - derived allele count
			 * - derived allele count (maximum likelihood)
			 * - derived allele frequency
			 * - derived allele frequencey (maximum likelihood)
			 * - number of missing genotypes
			 * - whether the outgroup nucleotide is on the same chromosome as the polymorphic site
			 * - whether the outgroup nucleotide is good quality
			 * - site quality score
			 *
			 *
			 * _NOTE_: The range must be confined to a single chromosome.
			 *
			 * \param[in] chromName chromosome name
			 * \param[in] start start position of the target range
			 * \param[in] end end position of the target range
			 * \param[out] sites vector of divergent site information (appended after execution)
			 *
			 */
			void getPolySites(const string &chromName, const uint64_t &start, const uint64_t &end, vector<string> &sites);
			/** \brief Get list of polymorphic sites from a vector of positions
			 *
			 * Get a list of polymorphic sites from a vector of positions. The provided vector of cromosome names must be the same length as the vector of genome positions.
			 * The chromosome names must be arranged in contiguous blocks, with the same order as in the target VCF file. This is to speed up file traversal.
			 * The vector of sites is appended by the function, so any existing information will be preserved.
			 * The site description is in a tab-delimited string with the following fields:
			 *
			 * - chromosome name
			 * - position
			 * - reference nucleotide
			 * - alternative nucleotide
			 * - which nucleotide is ancestral ('r' for reference, 'a' alternative, 'u' unknown)
			 * - derived allele count
			 * - derived allele count (maximum likelihood)
			 * - derived allele frequency
			 * - derived allele frequencey (maximum likelihood)
			 * - number of missing genotypes
			 * - whether the outgroup nucleotide is on the same chromosome as the polymorphic site
			 * - whether the outgroup nucleotide is good quality
			 * - site quality score
			 *
			 * \param[in] chromNames vector of chromosome names
			 * \param[in] positions vector of query site genome positions
			 * \param[out] sites vector of divergent site information (appended after execution)
			 *
			 */
			void getPolySites(const vector<string> &chromNames, const vector<uint64_t> &positions, vector<string> &sites);

		private:
			// Variables for the current record
			/// Current SNP position
			uint64_t varPos_;
			/// Reference nucleotide
			char refID_;
			/// Alternative nucleotide
			char altID_;
			/** \brief Which is ancestral
			 *
			 * 'r' if reference, 'a' if alternative, 'u' if unknown. Unknown could be because it's not biallelic or dvierged nucleotide is missing
			 */
			char ancState_;
			/** \brief Outgroup quality
			 *
			 * 1 if the outgroup nucleotide is good quality (uppercase), 0 otherwise.
			 */
			uint16_t outQual_;
			/// Is the divergent site on the same chromosome (1 yes, 0 no)?
			uint16_t sameChr_;
			/// Number missing
			uint32_t numMissing_;
			/// Number called
			uint32_t numCalled_;
			/// Reference allele count
			uint32_t refAC_;
			/// Reference allele count (maximum likelihood)
			uint32_t refMLAC_;
			/// Reference allele frequency
			double refAF_;
			/// Reference allele frequency (maximum likelihood)
			double refMLAF_;
			/// Site quality score
			double quality_;
			/// Chromosome ID
			string chrID_;
			/// Last completely searched chromosome
			string foundChr_;
			/// The full VCF line (record)
			string fullRecord_;

			/// The file stream
			fstream vcfFile_;
			/// The corresponding .axt object
			ParseAXT axtObj_;

			/// Parse current record
			void parseCurrentRecord_();
			/** Export current record
			 *
			 * \return string with the requisite site information
			 */
			string exportCurRecord_();
	};
}

#endif /* parseVCF_hpp */


