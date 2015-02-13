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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <octave/oct.h>

#ifdef HAVE_HDF5
#include <hdf5.h>

class H5File
{
public:
  H5File(const char *filename, const char *setname);
  ~H5File();
  int get_rank() {return rank;}
  int select_all();
  int select_hyperslab(const Matrix& start, const Matrix& count,
      const Matrix& stride, const Matrix& block);
  octave_value read();
      
private:
  int rank;
  hsize_t *h5_dims;
  dim_vector mat_dims;
  hid_t file;
  hid_t dataset;
  hid_t dataspace;
  hid_t memspace;
  NDArray read_double();
  ComplexNDArray read_complex();
};

bool hdf5_types_compatible (hid_t t1, hid_t t2);
hid_t hdf5_make_complex_type (hid_t num_type);

#endif
