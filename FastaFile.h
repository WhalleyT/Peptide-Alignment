//
// Created by tom on 29/08/18.
//

#include <string>
#include <unordered_set>
#include <vector>
#include <string>

#ifndef ALIGNMENT_HPC_FASTAFILE_H
#define ALIGNMENT_HPC_FASTAFILE_H


class FastaFile {
public:
    //variables
    int peplen;
    char *filename;
    int num_kmers;
    int total_indexes;
    char *peptide_char;
    unsigned long long total_kmers = 0ULL;
    std::vector<std::pair<std::string, std::string>> sequences;

    //functions
    void read_fasta();
    void allocate_char();
    FastaFile(char *name_in, int k_in);
    ~FastaFile();

protected:
    //variables
    int num_sequences;
    std::unordered_set<std::string> kmers;
};

#endif //ALIGNMENT_HPC_FASTAFILE_H
