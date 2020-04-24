/*
 * t8dg_tensor.c
 *
 *  Created on: Apr 24, 2020
 *      Author: lukas
 */

#include "t8dg.h"
#include <sc_containers.h>

void
t8dg_tensor_array_extract_vector (sc_array_t * tensor_array, const int ivector, const int stride, sc_array_t * vector)
{
  T8DG_ASSERT (stride > 0);
  T8DG_ASSERT (ivector >= 0);

  int                 ivec_element, tensor_array_index, vector_length;
  vector_length = vector->elem_count;
  tensor_array_index = (ivector / stride) * vector_length * stride + (ivector % stride);

  for (ivec_element = 0; ivec_element < vector_length; ivec_element++) {
    *(double *) sc_array_index_int (vector, ivec_element) = *(double *) sc_array_index_int (tensor_array, tensor_array_index);
    tensor_array_index += stride;
  }
}

void
t8dg_tensor_transform_tensoridx (int idx, const int tensordims[DIM3], int tensoridx[DIM3])
{
  /*TODO: Check */
  int                 numtensor = 0, itensor;
  while (numtensor < 3 && tensordims[numtensor] > 0) {
    numtensor++;
  }
  for (itensor = 0; itensor < numtensor; itensor++) {
    tensoridx[itensor] = idx % tensordims[itensor];
    idx /= tensordims[itensor];
  }
}

int
t8dg_tensor_mult_other_lengths (const int num_tensor, const int tensor_length[DIM3], const int itensor)
{
  int                 result = 1;
  int                 iothertensor;
  for (iothertensor = 0; iothertensor < num_tensor; iothertensor++) {
    if (iothertensor != itensor) {
      result *= tensor_length[iothertensor];
    }
  }
  return result;
}

void
t8dg_tensor_array_inject_vector (sc_array_t * vector, const int ivector, const int stride, sc_array_t * tensor_array)
{
  T8DG_ASSERT (ivector >= 0);
  T8DG_ASSERT (stride > 0);

  int                 ivec_element, tensor_array_index, vector_length;
  vector_length = vector->elem_count;
  tensor_array_index = (ivector / stride) * vector_length * stride + (ivector % stride);

  for (ivec_element = 0; ivec_element < vector_length; ivec_element++) {
    *(double *) sc_array_index_int (tensor_array, tensor_array_index) = *(double *) sc_array_index_int (vector, ivec_element);
    tensor_array_index += stride;
  }
}
