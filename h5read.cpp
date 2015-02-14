/*
 *
 *    Copyright 2012 Tom Mullins
 *    Copyright 2015 Tom Mullins, Stefan Großhauser
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
#include <octave/lo-ieee.h>
#include <cstdlib>
#include <string>
#include "gripes.h"

#include "h5read.doc.h"

using namespace std;

#if ((H5_VERS_MAJOR > 1) || (H5_VERS_MINOR >= 8))
#define HAVE_HDF5_18 1
#endif

#if defined(HAVE_HDF5) && defined(HAVE_HDF5_18)
#include "h5file.h"

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
    const char *name, bool allow_zeros)
{
  mat = val.matrix_value();
  if (error_state)
  {
    return 1;
  }

  if (!mat.is_vector())
  {
    error("%s must be a vector", name);
    return 1;
  }

  double mind, maxd;
  if (allow_zeros)
  {
    for (int i = 0; i < mat.nelem(); i++)
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

#endif

DEFUN_DLD(h5read, args, nargout, string((char*) h5read_doc))
{
#if !(defined(HAVE_HDF5) && defined(HAVE_HDF5_18))
  gripe_disabled_feature("h5read", "HDF5 IO");
  return octave_value_list();
#else
  int nargin = args.length();

  if (nargin < 2 || nargin == 3 || nargin > 6 || nargout > 1)
  {
    print_usage();
    return octave_value_list();
  }
  if ((args(0).is_string()==false) || (args(1).is_string()==false))
  {
    print_usage();
    return octave_value_list();
  }

  string filename = args(0).string_value();
  string dsetname = args(1).string_value();
  if (error_state)
    return octave_value_list();
  H5E_auto_t oef;
  void *olderr;
  H5Eget_auto(H5E_DEFAULT,&oef,&olderr);
  // suppress hdf5 error output
  //H5Eset_auto(H5E_DEFAULT,0,0);
  //open the hdf5 file
  H5File file(filename.c_str());
  // restore old setting
  //H5Eset_auto(H5E_DEFAULT,oef,olderr);
  if (nargin < 4)
  {
    octave_value retval = file.read_dset_complete(dsetname.c_str());
    return retval;
  }
  else
  {
    Matrix start, count, stride, block;
    int err = 0;
    
    err = err || check_vec(args(2), start, "START", false);
    start -= 1;

    err = err || check_vec(args(3), count, "COUNT", true);

    if (nargin < 5)
      stride = Matrix();
    else
      err = err || check_vec(args(4), stride, "STRIDE", false);

    if (nargin < 6)
      block = Matrix();
    else
      err = err || check_vec(args(5), block, "BLOCK", false);

    if (err)
      return octave_value_list();

    return file.read_dset_hyperslab(dsetname.c_str(),
				    start, count, stride, block, nargin-2);
  }
#endif
}

