//
// Created by tom on 29/08/18.
//

#include <iostream>
#include <fstream>
#include <unordered_set>
#include <string>

#include <stdio.h>
#include <cuda.h>
#include <signal.h>

#include "FastaFile.h"


void FastaFile::read_fasta() {

    std::ifstream input(filename);
    if(!input.good()){
        throw std::runtime_error("Error opening file");
    }

    std::string line, iter_name, content;
    while( std::getline( input, line ).good() ){
        if( line.empty() || line[0] == '>' ){
            if( !iter_name.empty() ) {

                int end = content.size() - peplen;
                for (int i = 0; i < end ; ++i) {
                    std::string peptide = content.substr(i, peplen);
                    kmers.insert(peptide);
                    std::pair <std::string, std::string> in = std::make_pair(iter_name, peptide);
                    sequences.push_back(in);
                    ++total_kmers;
                }

                iter_name.clear();
            }
            if(!line.empty()) {
                iter_name = line.substr(1);
            }
            content.clear();
        } else if( !iter_name.empty()) {
            if( line.find(' ') != std::string::npos ) {
                iter_name.clear();
                content.clear();
            } else {
                content += line;
            }
        }
    }

    if( !iter_name.empty()) {

        int end = content.size() - peplen;
        for (int i = 0; i < end ; ++i) {
            std::string peptide = content.substr(i, peplen);
            kmers.insert(peptide);
            std::pair <std::string, std::string> in = std::make_pair(iter_name, peptide);
            sequences.push_back(in);
            ++total_kmers;
        }
    }

    num_sequences = sequences.size();
    std::cout << "There are " << num_sequences << " sequences for " << filename << std::endl;
    std::cout << "There are " << kmers.size() << " unique  and "
              << total_kmers << " total " << peplen << "-mers" << std::endl;
    num_kmers = kmers.size();
}

FastaFile::FastaFile(char *name_in, int k_in) {
    filename = name_in;
    peplen = k_in;
    std::cout << "Initialising a FASTA object for " << filename << std::endl;
}

FastaFile::~FastaFile() {
    std::cout << "Cleaning up reference" << std::endl;
    delete [] peptide_char;
}

void FastaFile::allocate_char() {
    std::cout << "Allocating memory for " << filename << "'s array" << std::endl;
    total_indexes = total_kmers * peplen;
    peptide_char = new char[total_indexes];

    int index = 0;

    for (int i = 0; i < sequences.size(); ++i) {
        std::pair<std::string, std::string> pairing = sequences[i];
        std::string seq = pairing.second;
        for (int i = 0; i < peplen; ++i) {
            char AA = seq[i];
            peptide_char[index] = AA;
            ++index;
        }
    }
}
