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
#include <cerrno>
#include <iostream>
#include <algorithm>
#include <string>
#include "gripes.h"

using namespace std;

#if ((H5_VERS_MAJOR > 1) || (H5_VERS_MINOR >= 8))
// define this in case there is no configure script at work. This
// should not be necessary any more when integrated into core.
#define HAVE_HDF5_18 1
#endif

#if defined(HAVE_HDF5) && defined(HAVE_HDF5_18)
#include <hdf5.h>
#include "h5read.h"

bool
any_int_leq_zero(const Matrix& mat)
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
int
check_vec(const octave_value& val, Matrix& mat/*out*/,
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

DEFUN_DLD(h5read, args, nargout,
	  "-*- texinfo -*-\n\
@deftypefn {Loadable Function} {@var{data} =} h5read (@var{filename}, @var{dsetname})\n\
@deftypefnx {Loadable Function} {@var{data} =} h5read (@var{filename}, @var{dsetname}, @var{start}, @var{count})\n\
@deftypefnx {Loadable Function} {@var{data} =} h5read (@var{filename}, @var{dsetname}, @var{start}, @var{count}, @var{stride})\n\
@deftypefnx {Loadable Function} {@var{data} =} h5read (@var{filename}, @var{dsetname}, @var{start}, @var{count}, @var{stride}, @var{block})\n\
Read a hyperslab of data from an HDF5 file specified by its @var{filename}. \n\
The datatype will be coerced to double.\n\
For example:\n\
\n\
@example\n\
@group\n\
data = h5read (\"mydata.h5\", \"/grid/time\");\n\
@end group\n\
@end example\n\
\n\
The variable @var{dsetname} is the name of the dataset in the HDF5\n\
file to read. It has to be specified by its absolute path (or relative\n\
to the root group / ).\n\
\n\
The other four arguments are 1xn or nx1 matrices, where n is the number of dimensions\n\
in the dataset. If all are omitted, the entire dataset will be read.\n\
\n\
@var{start} is a 1-based starting offset\n\
@var{count} is the number of blocks to read. If 0 or Inf\n\
is specified in any dimension, as many blocks as possible \n\
are read in that dimension.\n\
@var{stride} is the offset between the start of each block. \n\
Defaults to a vector of ones.\n\
@var{block} is the size of each block to read. Defaults to a vector of ones.\n\
@seealso{h5write}\n\
@end deftypefn")
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

  //open the hdf5 file
  H5File file(filename.c_str());
  if (error_state)
    return octave_value_list();

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

DEFUN_DLD (h5readatt, args, nargout,
	   "-*- texinfo -*-\n\
@deftypefn {Loadable Function} {@var{data} =} h5readatt (@var{filename}, @var{objectname}, @var{attname})\n\
\n\
Reads one attribute of an object from an HDF5 file, specified by the\n\
@var{filename} and the @var{objectname}.\n\
The third argument @var{attname} is the name of the attribute which \n\
is to read.\n\
\n\
@seealso{h5writeatt}\n\
@end deftypefn")
{
  octave_value retval;
#if !(defined(HAVE_HDF5) && defined(HAVE_HDF5_18))
  gripe_disabled_feature("h5readatt", "HDF5 IO");
  return octave_value_list();
#else
  int nargin = args.length();
  if (nargin != 3)
    {
      print_usage();
      return retval;
    }
  if ((args(0).is_string()==false) || (args(1).is_string()==false) || (args(2).is_string()==false))
    {
      print_usage();
      return retval;
    }

  string filename = args(0).string_value();
  string objname = args(1).string_value();
  string attname = args(2).string_value();
  if (error_state)
    return octave_value_list();
  
  //open the hdf5 file
  H5File file(filename.c_str());
  if (error_state)
    return octave_value_list();
        
  retval = file.read_att(objname.c_str(), attname.c_str());
  return retval;

#endif
}


DEFUN_DLD(h5write, args, nargout,
	  "-*- texinfo -*-\n\
@deftypefn {Loadable Function} h5write (@var{filename}, @var{dsetname}, @var{data})\n\
\n\
Write a matrix @var{data} to the specified location @var{dsetname} in \n\
a HDF5 file specified by @var{filename}.\n\
\n\
@seealso{h5read}\n\
@end deftypefn")
{
#if !(defined(HAVE_HDF5) && defined(HAVE_HDF5_18))
  gripe_disabled_feature("h5write", "HDF5 IO");
  return octave_value_list();
#else
  int nargin = args.length();

  if (nargin != 3)
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
  string location = args(1).string_value();
  NDArray data = args(2).array_value();
  
  if (error_state)
    return octave_value_list();
    
  //open the hdf5 file
  H5File file(filename.c_str());
  if (error_state)
    return octave_value_list();

  file.write_dset(location.c_str(),
		  data);
  return octave_value_list();
#endif
}


DEFUN_DLD(h5writeatt, args, nargout,
	  "-*- texinfo -*-\n\
@deftypefn {Loadable Function} h5writeatt (@var{filename}, @var{objectname}, @var{attname}, @var{attvalue})\n\
\n\
Write an attribute with name @var{attname} and value @var{attvalue} to\n\
the object named @var{objectname} in the HDF5 file specified by @var{filename}.\n\
\n\
@seealso{h5readatt}\n\
@end deftypefn")
{
#if !(defined(HAVE_HDF5) && defined(HAVE_HDF5_18))
  gripe_disabled_feature("h5writeatt", "HDF5 IO");
  return octave_value_list();
#else
  int nargin = args.length();

  if (nargin != 4)
    {
      print_usage();
      return octave_value_list();
    }
  if ((args(0).is_string()==false) || (args(1).is_string()==false)
      || (args(2).is_string()==false))
    {
      print_usage();
      return octave_value_list();
    }

  string filename = args(0).string_value();
  string location = args(1).string_value();
  string attname = args(2).string_value();
  
  if (error_state)
    return octave_value_list();
    
  //open the hdf5 file
  H5File file(filename.c_str());
  if (error_state)
    return octave_value_list();

  file.write_att(location.c_str(), attname.c_str(),
  		 args(3));

  return octave_value_list();
#endif
}


#if defined(HAVE_HDF5) && defined(HAVE_HDF5_18)

H5File::H5File(const char *filename)
{
  H5E_auto_t oef;
  void *olderr;
  H5Eget_auto(H5E_DEFAULT,&oef,&olderr);
  //suppress hdf5 error output
  H5Eset_auto(H5E_DEFAULT,0,0);
  file = H5Fopen(filename, H5F_ACC_RDWR, H5P_DEFAULT); 
  // restore old setting
  H5Eset_auto(H5E_DEFAULT,oef,olderr);
  
  if (file < 0)
    {
      error("Opening the file failed, %s: %s", filename, strerror(errno));
      return;
    }
}

H5File::~H5File()
{
  if (H5Iis_valid(memspace_id))
    H5Sclose(memspace_id);

  if (H5Iis_valid(dspace_id))
    H5Sclose(dspace_id);

  if (H5Iis_valid(dset_id))
    H5Dclose(dset_id);
  
  if(H5Iis_valid(att_id))
    H5Aclose(att_id);

  if(H5Iis_valid(obj_id))
    H5Oclose(obj_id);

  if(H5Iis_valid(file))
    H5Fclose(file);

  if (h5_dims != NULL)
    {
      free(h5_dims);
      h5_dims = NULL;
    }
}

// T will be Matrix or dim_vector
template <typename T>
hsize_t*
alloc_hsize(const T& dim)
{
  int rank = dim.length();
  hsize_t *hsize = (hsize_t*)malloc(rank * sizeof(hsize_t));
  for (int i = 0; i < rank; i++)
    {
      hsize[i] = dim(i);
    }
  return hsize;
}


int
H5File::open_dset(const char *dsetname)
{
  dset_id = H5Dopen(file, dsetname, H5P_DEFAULT);
  if (dset_id < 0)
    {
      error("Error opening dataset %s", dsetname);
      return -1;
    }

  dspace_id = H5Dget_space(dset_id);
  if (dspace_id < 0)
    {
      error("Error opening dataspace %s", dsetname);
      return -1;
    }

  rank = H5Sget_simple_extent_ndims(dspace_id);
  if (rank < 0)
    {
      error("Error reading extent of %s", dsetname);
      return -1;
    }

  h5_dims = (hsize_t*)malloc(rank * sizeof(hsize_t));
  if (!h5_dims)
    {
      error("Error allocating memory for %s", dsetname);
      return -1;
    }
  if(H5Sget_simple_extent_dims(dspace_id, h5_dims, NULL) < 0)
    {
      error("Error determining dataset extent dimensions %s", dsetname);
      return -1;
    }
  
  return 0;
}

octave_value
H5File::read_dset_complete(const char *dsetname)
{
  if(open_dset(dsetname) < 0)
    return octave_value_list();

  mat_dims.resize(max(rank, 2));
  // .resize(1) still leaves mat_dims with a length of 2 for some reason, so
  // we need at least 2 filled
  mat_dims(0) = mat_dims(1) = 1;
  for (int i = 0; i < rank; i++)
    {
      mat_dims(i) = h5_dims[i];
    }

  if(H5Sselect_all(dspace_id) < 0)
    {
      error("Error selecting complete dataset %s", dsetname);
      return octave_value_list();
    }

  octave_value retval = read_dset();
  return retval;
}

octave_value
H5File::read_dset_hyperslab(const char *dsetname,
			    const Matrix& start, const Matrix& count,
			    const Matrix& stride, const Matrix& block,
			    int nargin)
{
  if(open_dset(dsetname) < 0)
    return octave_value_list();

  if (rank == 0 && !(start.is_empty() && count.is_empty()
		     && stride.is_empty() && block.is_empty()))
    {
      error("Cannot specify hyperslab for scalar datasets (rank 0)");
      return octave_value_list();
    }

  if(start.nelem() != rank)
    {
      error("start must be a vector of length %d, the dataset rank", rank);
      return octave_value_list();
    }
  if(count.nelem() != rank)
    {
      error("count must be a vector of length %d, the dataset rank", rank);
      return octave_value_list();
    }
  
  Matrix _stride = stride;
  if(nargin < 3)
    {
      _stride = Matrix(dim_vector(1, rank), 1);
    }
  if(_stride.nelem() != rank)
    {
      error("stride must be a vector of length %d, the dataset rank", rank);
      return octave_value_list();
    }
  Matrix _block = block;
  if(nargin < 4)
    {
      _block = Matrix(dim_vector(1, rank), 1);
    }
  if(_block.nelem() != rank)
    {
      error("block must be a vector of length %d, the dataset rank", rank);
      return octave_value_list();
    }
  
  // .resize(1) still leaves mat_dims with a length of 2 for some reason, so
  // we need at least 2 filled
  mat_dims.resize(max(rank, 2));
  mat_dims(0) = mat_dims(1) = 1;

  Matrix _count = count;
  for (int i = 0; i < rank; i++)
    {
      if (_stride(i) < _block(i))
	{
	  error("In dimension %d, requested stride %d smaller than block size %d",
		i+1, (int)_stride(i), (int)_block(i));
	  return octave_value_list();
	}
      if (_count(i) == 0)
	{
	  // a value of 0 (or Inf) means that as many blocks as possible
	  // shall be read in this dimension
	  _count(i) = (h5_dims[i] - start(i) - _block(i)) / _stride(i) + 1;
	}
      mat_dims(i) = _count(i)*_block(i);
      int end = start(i) + _stride(i)*(_count(i)-1) + _block(i); // exclusive
      if (h5_dims[i] < end)
	{
	  error("In dimension %d, dataset only has %d elements, but at least %d"
		" are required for requested hyperslab", i+1, (int)h5_dims[i],
		end);
	  return octave_value_list();
	}
    }

  hsize_t *hstart = alloc_hsize(start);
  hsize_t *hstride = alloc_hsize(_stride);
  hsize_t *hcount = alloc_hsize(_count);
  hsize_t *hblock = alloc_hsize(_block);
  // TODO check these (and hmem) for NULLs

  herr_t sel_result = H5Sselect_hyperslab(dspace_id, H5S_SELECT_SET, hstart,
					  hstride, hcount, hblock);

  free(hstart);
  free(hstride);
  free(hcount);
  free(hblock);

  if (sel_result < 0)
    {
      return octave_value_list();
    }

  octave_value retval = read_dset();
  return retval;
}

octave_value
H5File::read_dset()
{
  bool is_cmplx = false;
  hid_t type_id = H5Dget_type(dset_id);
  hid_t type_class_id = H5Tget_class(type_id);
  
  if (type_class_id == H5T_COMPOUND)
    {
      hid_t complex_type = hdf5_make_complex_type(H5T_NATIVE_DOUBLE);
      if (hdf5_types_compatible(type_id, complex_type))
	is_cmplx = true;
    }

  Matrix perm_vec(1, rank);
  if (rank >= 2)
    {
      // reverse the elements in mat_dims
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
  hid_t memspace_id = H5Screate_simple(rank, hmem, hmem);
  free(hmem);

  herr_t read_result;
  if (is_cmplx)
    {
      ComplexNDArray retval(mat_dims); 
      read_result = H5Dread(dset_id,
			    hdf5_make_complex_type(H5T_NATIVE_DOUBLE),
			    memspace_id, dspace_id,
			    H5P_DEFAULT, retval.fortran_vec());
      if (read_result < 0)
	return octave_value_list();
      if (rank >= 2)
	retval = retval.permute(perm_vec, false);
      return octave_value(retval);
    }
  else
    {
      NDArray retval(mat_dims);
      read_result = H5Dread(dset_id,
			    H5T_NATIVE_DOUBLE,
			    memspace_id, dspace_id,
			    H5P_DEFAULT, retval.fortran_vec());
      if (read_result < 0)
	return octave_value_list();
      else if (rank >= 2)
	retval = retval.permute(perm_vec, false);

      return octave_value(retval);
    }
}

void
H5File::write_dset(const char *dsetname,
			const NDArray& data)
{
  int rank = data.dims().length();
  // create a dims vector where the elements are reversed with respect
  // to the octave matrix.
  dim_vector new_mat_dims;
  Matrix perm_vec(1, rank);
  new_mat_dims.resize(rank);
  for (int i = 0; i < rank; i++)
    {
      int j = rank-i-1;
      if(rank >= 2)
	new_mat_dims(i) = data.dims()(j);
      else
	// this is not yet correct, a onedimensional matrix is transposed
	new_mat_dims(i) = data.dims()(i);
      
      perm_vec(i) = j;
    }
  hsize_t *dims = alloc_hsize(new_mat_dims);
  dspace_id = H5Screate_simple(rank, dims, NULL);

  //check if all groups in the path dsetname exist. if not, create them
  string path(dsetname);
  for(int i=1; i < path.length(); i++)
    {
      if(path[i] == '/')
	{
	  if(!H5Lexists(file, path.substr(0,i).c_str(), H5P_DEFAULT))
	    {
	      hid_t group_id = H5Gcreate(file, path.substr(0,i).c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	      H5Gclose(group_id);
	    }
	}
    }
  
  //check if the data set already exists. if it does, open it,
  //otherwise, create it.
  if(H5Lexists(file,dsetname,H5P_DEFAULT))
    {
      if(open_dset(dsetname) < 0)
	{
	  error("Could not open existing dataset in order to write to");
	}
    }
  else
    dset_id = H5Dcreate(file, dsetname, H5T_NATIVE_DOUBLE, dspace_id,
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if(rank >= 2)
    H5Dwrite(dset_id, H5T_NATIVE_DOUBLE,
	     H5S_ALL, H5S_ALL, H5P_DEFAULT,
	     data.permute(perm_vec,true).fortran_vec());
  else
    // this is not yet correct, a onedimensional matrix is transposed
    H5Dwrite(dset_id, H5T_NATIVE_DOUBLE,
	     H5S_ALL, H5S_ALL, H5P_DEFAULT,
	     data.fortran_vec());


  free(dims);
}

octave_value
H5File::read_att(const char *objname, const char *attname)
{
  octave_value retval;
  obj_id = H5Oopen(file, objname, H5P_DEFAULT);

  if (obj_id < 0)
    {
      error("h5readatt: opening the given object failed");
      return retval;
    }

  att_id = H5Aopen_name(obj_id, attname);
  if (att_id < 0)
    {
      error("h5readatt: opening the given attribute failed");
      return retval;
    }

  hid_t type = H5Aget_type(att_id);
  if (type < 0)
    {
      error("h5readatt: dataset type error");
      return retval;
    }

  size_t numVal = H5Aget_storage_size(att_id)/H5Tget_size(type);
  if (H5Tget_class(type) == H5T_FLOAT)
    {
      double value[numVal];
      if (H5Tget_size(type) == sizeof(float))
	{
	  float f_value[numVal];
	  if (H5Aread(att_id, H5T_NATIVE_FLOAT, f_value)<0)
	    {
	      error("h5readatt: reading the given float Attribute failed");
	      return retval;
	    }
	  for (size_t n = 0; n < numVal; ++n)
	    value[n] = f_value[n];
	}
      else if (H5Tget_size(type) == sizeof(double))
	{
	  if (H5Aread(att_id, H5T_NATIVE_DOUBLE, value)<0)
	    {
	      error("h5readatt: reading the given double Attribute failed");
	      return retval;
	    }
	}
      else
	{
	  error("h5readatt: reading the given float Attribute failed: \
cannot handle size of type");
	  return retval;
	}

      Matrix mat(numVal,1);
      for (size_t n=0;n<numVal;++n)
	mat(n)=value[n];
      retval = octave_value(mat);

    }
  else if (H5Tget_class(type)==H5T_INTEGER)
    {
      // Integer attributes are casted to floating point octave values
    
      double value[numVal];
      if (H5Tget_size(type)==sizeof(int))
	{
	  int f_value[numVal];
	  if (H5Aread(att_id, H5T_NATIVE_INT, f_value)<0)
	    {
	      error("h5readatt: reading the given integer Attribute failed");
	      return retval;
	    }
	  for (size_t n=0;n<numVal;++n)
	    value[n] = f_value[n]*1.0;
	}
      else
	{
	  error("h5readatt: reading the given integer Attribute failed: cannot handle size of type");
	  return retval;
	}

      Matrix mat(numVal,1);
      for (size_t n=0;n<numVal;++n)
	mat(n)=value[n];
      retval = octave_value(mat);
    
    }
  else if (H5Tget_class(type)==H5T_STRING)
    {
      // Size of each string:
      size_t size = H5Tget_size(type);
      // to read an array of strings (for future work): 
      //totsize = size*sdim[0]*sdim[1];
      // to read a single string:
      size_t totsize = size;
      // Set up read buffer for attribute
      char* buf = (char*)calloc(totsize, sizeof(char));
      if (H5Aread(att_id, type, buf)<0)
	{
	  error("h5readatt: reading the given string Attribute failed");
	  return retval;
	}
      retval = octave_value(buf);
    }
  else //none of the supported data types
    {
      error("h5readatt: attribute type not supported");
      return retval;
    }
  
  return retval;
}

void
H5File::write_att(const char *location, const char *attname,
		       const octave_value& attvalue)
{
  hsize_t *dims;
  if(attvalue.is_scalar_type())
    {
      dspace_id = H5Screate(H5S_SCALAR);
    }
  else if(attvalue.is_matrix_type())
    {
      error("matrix types are not yet supported.");
      return;
      // dims = alloc_hsize(attvalue.dims());
      // dspace_id = H5Screate_simple(attvalue.dims().length(), dims, NULL);
    }
  else
    {
      error("Only scalar or matrix types are supported.");
      return;
    }
  obj_id = H5Oopen(file, location, H5P_DEFAULT);

  hid_t attr_id;
  hid_t type_id;
  hid_t mem_type_id;
  void* buf;
  double attval_double;
  int attval_int;
    
  if(attvalue.is_real_type())
    {
      type_id = H5Tcopy(H5T_IEEE_F64LE);
      mem_type_id = H5Tcopy(H5T_NATIVE_DOUBLE);
      attval_double = attvalue.double_value();
      buf = (void *) &attval_double;
    }
  else if(attvalue.is_complex_type())
    {
      error("writing complex attributes is not yet supported");
    }
  else if(attvalue.is_integer_type())
    {
      type_id = H5Tcopy(H5T_STD_I64LE);
      mem_type_id = H5Tcopy(H5T_NATIVE_INT);
      attval_int = attvalue.int_value();
      buf = (void *) &attval_int;
    }
  else if(attvalue.is_string())
    {
      error("writing string attributes is not yet supported");
    }

  if(H5Aexists(obj_id,attname))
    {
      attr_id = H5Aopen(obj_id, attname, H5P_DEFAULT);
    }
  else
    {
      attr_id = H5Acreate(obj_id, attname, type_id,
			  dspace_id, H5P_DEFAULT, H5P_DEFAULT);
    }
  herr_t status = H5Awrite(attr_id, mem_type_id, buf);
  H5Aclose(attr_id);
  H5Tclose(type_id);
  H5Tclose(mem_type_id);
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

#endif
