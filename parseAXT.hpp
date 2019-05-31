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

using std::fstream;
using std::string;

namespace BayesicSpace {

	/** \brief .axt alignment parsing class
	 *
	 * Exatracts features from an .axt alignement file.
	 *
	 */
	class ParseAXT {
		public:
			/** \brief Default constructor */
			ParseAXT(){ axtFile_.exceptions(fstream::badbit); };
			/** \brief File name constructor
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
		private:
			/// The file stream
			fstream axtFile_;
			
			/// Current record's chromosome
			string chrID_;
			/// Is the conuterpart chromosome the same for the current record?
			bool sameChr_;
			/// Current record's start position
			uint64_t startPos_;
			/// Current record's primary sequence
			string primarySeq_;
			/// Current record's aligning sequence
			string alignSeq_;
			
			/** \brief Get next record */
			void getNextRecord_();
	};
}
#endif /* parseAXT_hpp */
