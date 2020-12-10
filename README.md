This set of software extracts polymorphism/divergence data from various files. These data will be used to do MacDonald-Kreitman type tests. This is for a specific project and data set. Generalization is possible, but may require modifying the code.

# Installation

Clone this repository, change directory to `polyDivExtract` and type

```sh
make
make install clean
```

This assumes that you are using `g++` as the C++ compiler. To use a different compiler, specify it on the `make` command line, e.g.

```sh
make CXX=c++
```

## Dependencies

No dependencies other than a C++ compiler that understands the C++11 standard.

# Usage

The `divSites` program takes nucleotide positions or ranges and returns a file with sites that have diverged between two species, inferred from the provided AXT between-species alignment file. Run it with

```sh
divSites -q file_with_positions -a AXT_alignment_file -o output_file
```

The query files should have at least two fields (chromosome ID and position). Chromosome IDs must be chrX, chr2L, chr2R, chr3L, or chr3R (this is for a _Drosophila_ data set). It is OK to omit "chr" from the chromosome arm names. Chromosomes must be listed in the same order in all files. The AXT file should have the same chromosome names as the query file. If there are exactly two fields, it is assumed that the file provides individual site positions ("positions file"). If there are more than two fields, it is assumed that the query file contains ranges of positions, with the first field indicating the chromosome arm, the second the start of the range, and the third the end. If there are more that three fields, the rest are ignored. Commented (starting with "#") and empty lines are ignored. The number of fields is checked on the first uncommented non-empty line. This line can be a header (defined as having non-numeric values in the position or start and/or end fields). It must have two fields for a positions file or no fewer than three fields for a ranges file. If the query file contains positions, the output file has the chromosome ID, position, focal species nucleotide, alternative (diverged) nucleotide, whether the alternative is on the same chromosome (1 if yes), and whether both nucleotides are good quality (1 if yes). The total number of good quality nucleotides per chromosome is listed as meta-data (commented out with `#`) at the start if the file. Of the query file has ranges, the output is similar but lists the "peak ID" (corresponding to each range) and number of good quality nucleotides in the range before the fields listed above, and no meta-data.

The `polySites` program extracts polymorphic sites. Run it with

```sh
polySites -q query_file -a AXT_alignment_file -v VCF_file -o output_file
```

Query files are the same as for `divSites`, as are the AXT files (the latter are used to call outgroup states). The VCF file contains polymorphism information. Chromosomes must be labeled the same as in the AXT and query files (with or without "chr" in front). The output file for position queries has the chromosome ID, position, reference nucleotide, alternative nucleotide, ancestral state (`r` if reference, `a` if alternative), derived allele count, maximum likelihood derived allele count, derived allele frequency, maximum likelihood derived allele frequency, number of missing genotypes, whether the outgroup site is on the same chromosome (1 if yes), whether the outgroup nucleotide is good quality (1 if yes), and the site quality score. The output is similar for a range query file, but includes "peak ID" (i.e., range ID).

The `fastaSort` program sorts FASTA files that have _loc=_ fields in their headers by start nucleotide position. If there are records with the same start position, only the longest one is kept. Run with

```sh
fastaSort -i input_FASTA -o output_FASTA
```

Having sorted a FASTA file with coding sequences, run

```sh
getFFsites -i input_sorted_FASTA -l log_file_name -o output_file
```

This extracts four-fold silent sites from each CDS, discarding regions of CDS overlap. The output lists the chromosome, FBgn number of the CDS, and chromosome position of the site. The log file contains debugging information, flags overlapping CDS, and highlights potentially problematic records.

The software assumes _Drosophila_ chromosomes, labeled as chrX, chr2L, etc.
