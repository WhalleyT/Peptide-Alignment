#ifndef DEVICE_FUNCTIONS_CUH
#define DEVICE_FUNCTIONS_CUH

__global__ void score_alignment_2D(char *query, char *reference, int *matrix, char *key, int *scores,
  unsigned long long *query_kmers, unsigned long long *ref_kmers, int *peplen, int *x_modifier, int *y_modifier);

#endif
