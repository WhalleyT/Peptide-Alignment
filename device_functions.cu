#include <stdio.h>
#include <cuda.h>
#include <signal.h>

__global__ void score_alignment_2D(char *query, char *reference, int *matrix, char *key, int *scores,
  unsigned long long *query_kmers, unsigned long long *ref_kmers, int *peplen, int *x_modifier, int *y_modifier){

  int idx = ((blockIdx.x * blockDim.x) + threadIdx.x) + *x_modifier;
  int idy = ((blockIdx.y * blockDim.y) + threadIdx.y) + *y_modifier;

  int ref_start = idx * *peplen;
  int que_start = idy * *peplen;

  if(idx < *ref_kmers && idy < *query_kmers){
    for (int j = 0; j < *peplen ; j++) {

      char ref_AA = reference[ref_start + j];
      char que_AA = query[que_start + j];

      bool ref_found = false;
      bool que_found = false;

      int ref_index;
      int que_index;

      for (int k = 0; k < 25; k++) {
          if(ref_AA == key[k]){
            ref_index = k;
            ref_found = true;
          }
          if(que_AA == key[k]){
            que_index = k;
            que_found = true;
          }
        }

        //now score by making index out of the two
        if(ref_found && que_found){
          int lookup_index = (ref_index * 25) + que_index;
          atomicAdd(&scores[idx], matrix[lookup_index]);
        }
      }
    }
  }
