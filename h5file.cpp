/*
 *
 *    Copyright 2012 Tom Mullins
 *
 *
 *    This file is part of hdf5oct.
 *
 *    hdf5oct is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU Lesser General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    hdf5oct is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU Lesser General Public License for more details.
 *
 *    You should have received a copy of the GNU Lesser General Public License
 *    along with hdf5oct.  If not, see <http://www.gnu.org/licenses/>.
 *
 *
 */

#include <octave/oct.h>
#include <octave/lo-ieee.h>
#include <hdf5.h>
#include <cerrno>
#include <iostream>
#include <algorithm>
#include "h5file.h"
using namespace std;

H5File::H5File(const char *filename, const char *setname)
{
  dataset = -1;
  dataspace = -1;
  memspace = -1;
  h5_dims = NULL;
  rank = -1;

  file = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT); 
  if (file < 0)
  {
    error("%s: %s", filename, strerror(errno));
    return;
  }

  dataset = H5Dopen(file, setname, H5P_DEFAULT);
  if (dataset < 0)
  {
    error("Error opening dataset %s", setname);
    return;
  }

  dataspace = H5Dget_space(dataset);
  if (dataspace < 0)
  {
    return;
  }

  rank = H5Sget_simple_extent_ndims(dataspace);
  if (rank < 0)
  {
    return;
  }

  h5_dims = (hsize_t*)malloc(rank * sizeof(hsize_t));
  if (!h5_dims || H5Sget_simple_extent_dims(dataspace, h5_dims, NULL) < 0)
  {
    rank = -1;
    return;
  }

  // .resize(1) still leaves mat_dims with a length of 2 for some reason, so
  // we need at least 2 filled
  mat_dims.resize(max(rank, 2));
  mat_dims(0) = mat_dims(1) = 1;
}

H5File::~H5File()
{

  if (h5_dims)
  {
    free(h5_dims);
    h5_dims = NULL;
  }

  if (memspace >= 0)
  {
    H5Sclose(memspace);
    memspace = -1;
  }

  if (dataspace >= 0)
  {
    H5Sclose(dataspace);
    dataspace = -1;
  }

  if (dataset >= 0)
  {
    H5Dclose(dataset);
    dataset = -1;
  }

  if (file >= 0)
  {
    H5Fclose(file);
    file = -1;
  }
}

// T will be Matrix or dim_vector
template <typename T>
hsize_t* alloc_hsize(const T& dim)
{
  int rank = dim.length();
  hsize_t *hsize = (hsize_t*)malloc(rank * sizeof(hsize_t));
  for (int i = 0; i < rank; i++)
  {
    hsize[i] = dim(i);
  }
  return hsize;
}

int H5File::select_all()
{
  for (int i = 0; i < rank; i++)
  {
    mat_dims(i) = h5_dims[i];
  }
  return dataspace < 0 || H5Sselect_all(dataspace) < 0;
}

int H5File::select_hyperslab(const Matrix& start, const Matrix& count,
    const Matrix& stride, const Matrix& block)
{
  if (rank == 0)
  {
    return 1;
  }

  Matrix _count = count;
  for (int i = 0; i < rank; i++)
  {
    if (stride(i) < block(i))
    {
      error("In dimension %d, requested stride %d smaller than block size %d",
          i+1, (int)stride(i), (int)block(i));
      return 1;
    }
    if (_count(i) == 0)
    {
      _count(i) = (h5_dims[i] - start(i) - block(i)) / stride(i) + 1;
    }
    mat_dims(i) = _count(i)*block(i);
    int end = start(i) + stride(i)*(_count(i)-1) + block(i); // non-inclusive
    if (h5_dims[i] < end)
    {
      error("In dimension %d, dataset only has %d elements, but at least %d"
            " required for requested hyperslab", i+1, (int)h5_dims[i], end);
      return 1;
    }
  }

  hsize_t *hstart = alloc_hsize(start);
  hsize_t *hstride = alloc_hsize(stride);
  hsize_t *hcount = alloc_hsize(_count);
  hsize_t *hblock = alloc_hsize(block);
  // TODO check these (and hmem) for NULLs

  herr_t sel_result = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, hstart,
      hstride, hcount, hblock);

  free(hstart);
  free(hstride);
  free(hcount);
  free(hblock);

  if (sel_result < 0)
  {
    return 1;
  }

  return 0;
}

octave_value H5File::read()
{
    hid_t type_id = -1, type_class_id = -1;
    bool is_cmplx=false;
    type_id = H5Dget_type (dataset);
    type_class_id = H5Tget_class (type_id);
    if (type_class_id == H5T_COMPOUND) {
	hid_t complex_type = hdf5_make_complex_type (H5T_NATIVE_DOUBLE);
        if (hdf5_types_compatible (type_id, complex_type))
		is_cmplx=true;
    }
    if (is_cmplx)
	return octave_value(read_complex());
    else
	return octave_value(read_double());
}

ComplexNDArray H5File::read_complex()
{
  if (dataspace < 0)
  {
    return NDArray();
  }

  Matrix perm_vec(1, rank);
  if (rank >= 2)
  {
    dim_vector new_mat_dims;
    new_mat_dims.resize(rank);
    for (int i = 0; i < rank; i++)
    {
      int j = rank-i-1;
      new_mat_dims(i) = mat_dims(j);
      perm_vec(i) = j;
    }
    mat_dims = new_mat_dims;
  }

  hsize_t *hmem = alloc_hsize(mat_dims);
  hid_t memspace = H5Screate_simple(rank, hmem, hmem);
  free(hmem);
  if (memspace < 0)
  {
    return NDArray();
  }

  ComplexNDArray mat(mat_dims);

  herr_t read_result = H5Dread(dataset, hdf5_make_complex_type (H5T_NATIVE_DOUBLE), memspace, dataspace,
      H5P_DEFAULT, mat.fortran_vec());
  if (read_result < 0)
  {
    mat = NDArray();
  }

  else if (rank >= 2)
  {
    mat = mat.permute(perm_vec);
  }

  return mat;
}


NDArray H5File::read_double()
{
  if (dataspace < 0)
  {
    return NDArray();
  }

  Matrix perm_vec(1, rank);
  if (rank >= 2)
  {
    dim_vector new_mat_dims;
    new_mat_dims.resize(rank);
    for (int i = 0; i < rank; i++)
    {
      int j = rank-i-1;
      new_mat_dims(i) = mat_dims(j);
      perm_vec(i) = j;
    }
    mat_dims = new_mat_dims;
  }

  hsize_t *hmem = alloc_hsize(mat_dims);
  hid_t memspace = H5Screate_simple(rank, hmem, hmem);
  free(hmem);
  if (memspace < 0)
  {
    return NDArray();
  }

  NDArray mat(mat_dims);

  herr_t read_result = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace,
      H5P_DEFAULT, mat.fortran_vec());
  if (read_result < 0)
  {
    mat = NDArray();
  }

  else if (rank >= 2)
  {
    mat = mat.permute(perm_vec);
  }

  return mat;
}


bool
hdf5_types_compatible (hid_t t1, hid_t t2)
{
  int n;
  if ((n = H5Tget_nmembers (t1)) != H5Tget_nmembers (t2))
    return false;

  for (int i = 0; i < n; ++i)
    {
      hid_t mt1 = H5Tget_member_type (t1, i);
      hid_t mt2 = H5Tget_member_type (t2, i);

      if (H5Tget_class (mt1) != H5Tget_class (mt2))
        return false;

      H5Tclose (mt2);
      H5Tclose (mt1);
    }

  return true;
}

hid_t
hdf5_make_complex_type (hid_t num_type)
{
  hid_t type_id = H5Tcreate (H5T_COMPOUND, sizeof (double) * 2);

  H5Tinsert (type_id, "real", 0 * sizeof (double), num_type);
  H5Tinsert (type_id, "imag", 1 * sizeof (double), num_type);

  return type_id;
}

