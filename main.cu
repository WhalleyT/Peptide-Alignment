#include <iostream>
#include <fstream>

#include "FastaFile.h"
#include "AlignmentMatrix.h"
#include "device_functions.cuh"

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, char *file, int line, bool abort=true)
{
   if (code != cudaSuccess)
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}

int main(int argc, char **argv) {

    if(argc != 5){
        std::cerr << "4 command line arguments must be supplied." << std::endl;
        std::cerr << std::endl;
        std::cerr << "They are as follows:" << std::endl;
        std::cerr << "-Query FASTA" << std::endl;
        std::cerr << "-Reference FASTA" << std::endl;
        std::cerr << "-Peptide length" << std::endl;
        std::cerr << "-GPU Device" << std::endl;
        std::cerr << std::endl;
        return 1;
    }

    //cudaSetDevice(std::stoi(argv[4]));
    //std::cout << "CUDA set to device " << argv[4] << std::endl;

    FastaFile query(argv[1], std::stoi(argv[3]));
    FastaFile reference(argv[2], std::stoi(argv[3]));

    query.read_fasta();
    reference.read_fasta();

    //allocate peplen * numpeps array
    query.allocate_char();
    reference.allocate_char();

    int matdim = 25;

    char letter_order[] = {'A', 'R', 'N', 'D', 'C', 'Q',
                           'E', 'G', 'H', 'I', 'L', 'K',
                           'M', 'F', 'P', 'S', 'T', 'W',
                           'Y', 'V', 'B', 'J', 'Z', 'X',
                           '*', '\0'};

    std::vector<int> scores = {6, -7, -4, -3, -6, -4, -2, -2, -7, -5, -6, -7, -5, -8, -2, 0, -1, -13, -8, -2, -3, -6,
                               -3, -1,
                               -17, -7, 8, -6, -10, -8, -2, -9, -9, -2, -5, -8, 0, -4, -9, -4, -3, -6, -2, -10, -8, -7,
                               -7, -4, -1,
                               -17, -4, -6, 8, 2, -11, -3, -2, -3, 0, -5, -7, -1, -9, -9, -6, 0, -2, -8, -4, -8, 6, -6,
                               -3, -1,
                               -17, -3, -10, 2, 8, -14, -2, 2, -3, -4, -7, -12, -4, -11, -15, -8, -4, -5, -15, -11, -8,
                               6, -10, 1, -1,
                               -17, -6, -8, -11, -14, 10, -14, -14, -9, -7, -6, -15, -14, -13, -13, -8, -3, -8, -15, -4,
                               -6, -12,
                               -9, -14, -1, -17, -4, -2, -3, -2, -14, 8, 1, -7, 1, -8, -5, -3, -4, -13, -3, -5, -5, -13,
                               -12, -7,
                               -3, -5, 6, -1, -17, -2, -9, -2, 2, -14, 1, 8, -4, -5, -5, -9, -4, -7, -14, -5, -4, -6,
                               -17, -8, -6,
                               1, -7, 6, -1, -17, -2, -9, -3, -3, -9, -7, -4, 6, -9, -11, -10, -7, -8, -9, -6, -2, -6,
                               -15, -14,
                               -5, -3, -10, -5, -1, -17, -7, -2, 0, -4, -7, 1, -5, -9, 9, -9, -6, -6, -10, -6, -4, -6,
                               -7, -7, -3,
                               -6, -1, -7, -1, -1, -17, -5, -5, -5, -7, -6, -8, -5, -11, -9, 8, -1, -6, -1, -2, -8, -7,
                               -2, -14,
                               -6, 2, -6, 5, -6, -1, -17, -6, -8, -7, -12, -15, -5, -9, -10, -6, -1, 7, -8, 1, -3, -7,
                               -8, -7, -6,
                               -7, -2, -9, 6, -7, -1, -17, -7, 0, -1, -4, -14, -3, -4, -7, -6, -6, -8, 7, -2, -14, -6,
                               -4, -3, -12,
                               -9, -9, -2, -7, -4, -1, -17, -5, -4, -9, -11, -13, -4, -7, -8, -10, -1, 1, -2, 11, -4,
                               -8, -5, -4,
                               -13, -11, -1, -10, 0, -5, -1, -17, -8, -9, -9, -15, -13, -13, -14, -9, -6, -2, -3, -14,
                               -4, 9, -10,
                               -6, -9, -4, 2, -8, -10, -2, -13, -1, -17, -2, -4, -6, -8, -8, -3, -5, -6, -4, -8, -7, -6,
                               -8, -10,
                               8, -2, -4, -14, -13, -6, -7, -7, -4, -1, -17, 0, -3, 0, -4, -3, -5, -4, -2, -6, -7, -8,
                               -4, -5, -6,
                               -2, 6, 0, -5, -7, -6, -1, -8, -5, -1, -17, -1, -6, -2, -5, -8, -5, -6, -6, -7, -2, -7,
                               -3, -4, -9,
                               -4, 0, 7, -13, -6, -3, -3, -5, -6, -1, -17, -13, -2, -8, -15, -15, -13, -17, -15, -7,
                               -14, -6, -12,
                               -13, -4, -14, -5, -13, 13, -5, -15, -10, -7, -14, -1, -17, -8, -10, -4, -11, -4, -12, -8,
                               -14, -3,
                               -6, -7, -9, -11, 2, -13, -7, -6, -5, 10, -7, -6, -7, -9, -1, -17, -2, -8, -8, -8, -6, -7,
                               -6, -5,
                               -6, 2, -2, -9, -1, -8, -6, -6, -3, -15, -7, 7, -8, 0, -6, -1, -17, -3, -7, 6, 6, -12, -3,
                               1, -3,
                               -1, -6, -9, -2, -10, -10, -7, -1, -3, -10, -6, -8, 6, -8, 0, -1, -17, -6, -7, -6, -10,
                               -9, -5, -7, -10,
                               -7, 5, 6, -7, 0, -2, -7, -8, -5, -7, -7, 0, -8, 6, -6, -1, -17, -3, -4, -3, 1, -14, 6, 6,
                               -5, -1,
                               -6, -7, -4, -5, -13, -4, -5, -6, -14, -9, -6, 0, -6, 6, -1, -17, -1, -1, -1, -1, -1, -1,
                               -1, -1,
                               -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -17, -17, -17, -17, -17,
                               -17, -17,
                               -17, -17, -17, -17, -17, -17, -17, -17, -17, -17, -17, -17, -17, -17, -17, -17, -17, -17,
                               1};


    AlignmentMatrix PAM30(letter_order, letter_order, matdim, scores);


    //now we're more or less setup, let's start the CUDA stuff
    std::cout << "Allocating memory for CUDA kernel" << std::endl;

    int *average_scores = new int [query.total_kmers];

    for (int i = 0; i < query.total_kmers; i++) {
      average_scores[i] = 0;
    }

    int peplen = std::stoi(argv[3]);

    const int threads = 32;
    const int blocks = 65535;
    const int total = threads * blocks;

    int num_x_calls  = query.total_kmers / total + (query.total_kmers % total == 0 ? 0:1);
    int num_y_calls  = reference.total_kmers / total + (reference.total_kmers % total == 0 ? 0:1);

    printf("%i x %i CUDA calls will be made\n", num_x_calls, num_y_calls);

    for (int i = 0; i < num_x_calls; i++) {
      for (int j = 0; j < num_y_calls; j++) {
        int x_modifier = i * total;
        int y_modifier = j * total;

        std::cout << "Operating on chunk " << i << " of the query FASTA and " << j << " of the reference FASTA" << std::endl;

        char *d_query;
        char *d_ref;
        int  *d_mat;
        char *d_key;
        int *d_scores;
        unsigned long long int *d_query_numpeps;
        unsigned long long int *d_ref_numpeps;
        int *d_peplen;
        int *d_x_modifier;
        int *d_y_modifier;

        gpuErrchk(cudaMalloc(&d_query,    query.total_indexes * sizeof(char)));
        gpuErrchk(cudaMalloc(&d_ref,      reference.total_indexes * sizeof(char)));
        gpuErrchk(cudaMalloc(&d_mat,      PAM30.mat_size * sizeof(int)));
        gpuErrchk(cudaMalloc(&d_key,      matdim * sizeof(char)));
        gpuErrchk(cudaMalloc(&d_scores,   query.total_indexes * sizeof(int)));
        gpuErrchk(cudaMalloc(&d_query_numpeps,  sizeof(unsigned long long int)));
        gpuErrchk(cudaMalloc(&d_ref_numpeps,  sizeof(unsigned long long int)));
        gpuErrchk(cudaMalloc(&d_peplen, sizeof(int)));
        gpuErrchk(cudaMalloc(&d_x_modifier, sizeof(int)));
        gpuErrchk(cudaMalloc(&d_y_modifier, sizeof(int)));
        //copy to GPU
        std::cout << "Copying to GPU memory" << std::endl;

        gpuErrchk(cudaMemcpy(d_query,             query.peptide_char,     query.total_indexes     * sizeof(char),   cudaMemcpyHostToDevice));
        gpuErrchk(cudaMemcpy(d_ref,               reference.peptide_char, reference.total_indexes * sizeof(char),   cudaMemcpyHostToDevice));
        gpuErrchk(cudaMemcpy(d_mat,               PAM30.scoring_matrix,   PAM30.mat_size          * sizeof(int),    cudaMemcpyHostToDevice));
        gpuErrchk(cudaMemcpy(d_key,               PAM30.x_order,          matdim                  * sizeof(char),   cudaMemcpyHostToDevice));
        gpuErrchk(cudaMemcpy(d_scores,            average_scores,         query.total_kmers     * sizeof(int), cudaMemcpyHostToDevice));
        gpuErrchk(cudaMemcpy(d_query_numpeps,     &query.total_kmers,              sizeof(unsigned long long int),    cudaMemcpyHostToDevice));
        gpuErrchk(cudaMemcpy(d_ref_numpeps,       &reference.total_kmers,          sizeof(unsigned long long int),    cudaMemcpyHostToDevice));
        gpuErrchk(cudaMemcpy(d_peplen,            &peplen,                                          sizeof(int),    cudaMemcpyHostToDevice));
        gpuErrchk(cudaMemcpy(d_x_modifier, &x_modifier, sizeof(int), cudaMemcpyHostToDevice));
        gpuErrchk(cudaMemcpy(d_y_modifier, &y_modifier, sizeof(int), cudaMemcpyHostToDevice));

        dim3 thread_pair(threads, threads, 1);

        int xblock  = blocks;
        int yblock  = blocks;

        std::cout << "Kernel will be called across " << xblock << " blocks of " << threads <<
        " threads in the x dimension and " << yblock << " blocks of " << threads << " in the y dimension" << std::endl;

        dim3 block_pair(xblock, yblock, 1);

        std::cout << "Calling CUDA kernel..." << std::endl;

        score_alignment_2D<<<block_pair, thread_pair>>>(d_query, d_ref, d_mat, d_key, d_scores, d_query_numpeps,
                                                        d_ref_numpeps, d_peplen, d_x_modifier, d_y_modifier);
        cudaThreadSynchronize();

        cudaError_t err = cudaGetLastError();

        if(err == cudaSuccess){
          std::cout << "CUDA kernel successfully finished" << std::endl;
        }else{
          std::cerr << "CUDA error: " << err << std::endl;
          return 1;
        }

        gpuErrchk(cudaMemcpy(average_scores, d_scores, query.total_kmers * sizeof(int), cudaMemcpyDeviceToHost));
        gpuErrchk(cudaMemcpy(query.peptide_char, d_query,     query.total_indexes     * sizeof(char),   cudaMemcpyDeviceToHost));

        std::cout << "cudaMemcpy completed" << std::endl;
        cudaThreadSynchronize();
        cudaDeviceReset();
        std::cout << "Device reset" << std::endl;
      }
    }

    double *out_average_scores = new double [query.total_kmers];

    for (int i = 0; i < query.total_kmers; i++) {
      out_average_scores[i] =  double(average_scores[i]) / double(reference.total_kmers);
    }

    std::cout << "Writing..." << std::endl;

    std::ofstream outfile("alignment_scores.txt");
    outfile << "peptide\tprotein\tscore" << std::endl;

    long long unsigned int peptide_idx = 0ULL;

    for (int i = 0; i < query.total_kmers; i++) {
      for (int j = 0; j < peplen; j++) {
        outfile << query.peptide_char[peptide_idx];
        ++peptide_idx;
      }
      outfile << "\t" << query.sequences[i].first << "\t" << out_average_scores[i] << std::endl;
    }

    delete [] out_average_scores;
    delete [] average_scores;
}
