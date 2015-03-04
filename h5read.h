/*
 *
 *    Copyright (C) 2012 Tom Mullins
 *    Copyright (C) 2015 Tom Mullins, Thorsten Liebig, Stefan Großhauser
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

#ifdef HAVE_CONFIG_H
#include <config.h>
// integrated into the GNU Octave build
#include "oct.h"
#else
// as a package
#include <octave/oct.h>
#endif


#if defined (HAVE_HDF5) && defined (HAVE_HDF5_18)
#include <hdf5.h>

class H5File
{
  
 public:
  
  H5File (const char *filename, const bool create_if_nonexisting);
  
  ~H5File ();
  
  octave_value read_dset_complete (const char *dsetname);
  octave_value read_dset_hyperslab (const char *dsetname,
                                    const Matrix& start, const Matrix& count,
                                    const Matrix& stride, const Matrix& block,
                                    int nargin);

  void write_dset (const char *location,
                   const octave_value ov_data);
  void write_dset_hyperslab (const char *location,
                             const octave_value ov_data,
                             const Matrix& start, const Matrix& count,
                             const Matrix& stride, const Matrix& block,
                             int nargin);
  octave_value read_att (const char *location, const char *attname);
  void write_att (const char *location, const char *attname,
                  const octave_value& attvalue);
  void create_dset (const char *location, const Matrix& size,
                    const char *datatype, Matrix& chunksize);
 private:
  const static int ALLOC_HSIZE_INFZERO_TO_UNLIMITED = 1;
  const static int ALLOC_HSIZE_INF_TO_ZERO = 2;
  const static int ALLOC_HSIZE_DEFAULT = 3;
  
  //rank of the hdf5 dataset
  int rank;
  //dimensions of the hdf5 dataset
  hsize_t *h5_dims = NULL;
  hsize_t *h5_maxdims = NULL;

  hid_t file;
  hid_t dset_id;
  hid_t dspace_id;
  hid_t memspace_id;
  hid_t obj_id;
  hid_t att_id;
  hid_t type_id;
  hid_t mem_type_id;

  //dimensions of the returned octave matrix
  dim_vector mat_dims;
  
  int open_dset (const char *dsetname);
  octave_value read_dset ();
  Matrix get_auto_chunksize (const Matrix& size, int typesize);

  template <typename T> hsize_t* alloc_hsize (const T& dim, const int mode, const bool reverse);

};



#endif
