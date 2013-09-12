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
#include <cstdlib>
#include <string>
#include "h5file.h"
#include "h5read.doc.h"
using namespace std;

bool any_int_leq_zero(const Matrix& mat)
{
  for (int i = 0; i < mat.length(); i++)
  {
    if (mat(i) < 0.5)
    {
      return true;
    }
  }
  return false;
}

// if ALLOW_ZEROS, then Infs will be converted to 0s (used for COUNT)
// else, 0s will produce an error message (used for others)
int check_vec(const octave_value& val, Matrix& mat/*out*/,
    const char *name, int rank, bool allow_zeros)
{
  mat = val.matrix_value();
  if (error_state)
  {
    return 1;
  }

  if (rank == 0 && mat.is_empty())
  {
    return 0;
  }

  if (!mat.is_vector() || mat.nelem() != rank)
  {
    error("%s must be a vector of length %d, the dataset rank", name, rank);
    return 1;
  }

  double mind, maxd;
  if (allow_zeros)
  {
    for (int i = 0; i < rank; i++)
    {
      if (mat(i) == octave_Inf)
      {
        mat(i) = 0;
      }
    }
    if (!mat.all_integers(mind, maxd) || mat.any_element_is_negative())
    {
      error("%s can only contain non-negative integers", name);
      return 1;
    }
  }
  else if (!mat.all_integers(mind, maxd) || any_int_leq_zero(mat))
  {
    error("%s can only contain positive integers", name);
    return 1;
  } 

  return 0;
}

DEFUN_DLD(h5read, args, nargout, string((char*) h5read_doc))
{
  int nargin = args.length();

  if (nargin < 2 || nargin == 3 || nargin > 6 || nargout > 1)
  {
    print_usage();
    return octave_value_list();
  }

  string filename = args(0).string_value();
  string dataset = args(1).string_value();
  if (error_state)
    return octave_value_list();
    
  void *olderr;
  H5E_auto_t oef;
  
  H5Eget_auto(H5E_DEFAULT,&oef,&olderr);
  H5Eset_auto(H5E_DEFAULT,0,0);
  H5File file(filename.c_str(), dataset.c_str());
  H5Eset_auto(H5E_DEFAULT,oef,olderr);

  int rank = file.get_rank();
  if (rank < 0)
    return octave_value_list();

  if (nargin < 4)
  {
    if (file.select_all())
      return octave_value_list();
  }
  else if (rank == 0)
  {
    error("Cannot specify hyperslab for scalar datasets (rank 0)");
    return octave_value_list();
  }
  else
  {
    Matrix start, count, stride, block;
    int err = 0;

    err = err || check_vec(args(2), start, "START", rank, false);
    err = err || check_vec(args(3), count, "COUNT", rank, true);

    if (nargin < 5) stride = Matrix(dim_vector(1, rank), 1);
    else err = err || check_vec(args(4), stride, "STRIDE", rank, false);

    if (nargin < 6) block = Matrix(dim_vector(1, rank), 1);
    else err = err || check_vec(args(5), block, "BLOCK", rank, false);

    if (err)
      return octave_value_list();

    start -= 1;

    if (file.select_hyperslab(start, count, stride, block))
      return octave_value_list();
  }
  
//  return octave_value(file.read());
  return file.read();
}
