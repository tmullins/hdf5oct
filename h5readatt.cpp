/*
 *
 *    Copyright 2015 Thorsten Liebig
 *
 *    Previous versions of this file were part of openEMS, and openEMS
 *    is free software: you can redistribute it and/or modify it under
 *    the terms of the GNU General Public License as published by the
 *    Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    With permission of the original author of this file, it is now
 *    is part of hdf5oct and the licence of this file is now the same
 *    as that of hdf5oct, the GNU Lesser General Public License.
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
#include <octave/ov-struct.h>
#include <cstdlib>
#include <string>
#include "gripes.h"

#include "h5readatt.doc.h"

using namespace std;

#if ((H5_VERS_MAJOR > 1) || (H5_VERS_MINOR >= 8))
#define HAVE_HDF5_18 1
#endif

#if defined(HAVE_HDF5) && defined(HAVE_HDF5_18)
#include "h5file.h"
#endif

DEFUN_DLD (h5readatt, args, nargout, string((char*) h5readatt_doc))
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
  if (error_state)
    return octave_value_list();
  
  H5E_auto_t oef;
  void *olderr;
  H5Eget_auto(H5E_DEFAULT,&oef,&olderr);
  //suppress hdf5 error output
  H5Eset_auto(H5E_DEFAULT,0,0);
  //open the hdf5 file
  hid_t file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  // restore old setting
  H5Eset_auto(H5E_DEFAULT,oef,olderr);
        
  if (file_id < 0)
  {
    error("h5readatt: opening the given file failed");
    return retval;
  }

  hid_t obj = H5Oopen(file_id, objname.c_str(), H5P_DEFAULT);

  if (obj < 0)
  {
    H5Oclose(obj);
    H5Fclose(file_id);
    error("h5readatt: opening the given object failed");
    return retval;
  }

  hid_t attr = H5Aopen_name(obj, args(2).string_value().c_str());
  if (attr < 0)
  {
    H5Oclose(obj);
    H5Fclose(file_id);
    error("h5readatt: opening the given attribute failed");
    return retval;
  }

  hid_t type = H5Aget_type(attr);
  if (type < 0)
  {
    H5Aclose(attr);
    H5Oclose(obj);
    H5Fclose(file_id);
    error("h5readatt: dataset type error");
    return retval;
  }

  size_t numVal = H5Aget_storage_size(attr)/H5Tget_size(type);
  if (H5Tget_class(type) == H5T_FLOAT)
  {
    double value[numVal];
    if (H5Tget_size(type) == sizeof(float))
    {
      float f_value[numVal];
      if (H5Aread(attr, H5T_NATIVE_FLOAT, f_value)<0)
      {
	H5Aclose(attr);
	H5Oclose(obj);
	H5Fclose(file_id);
	error("h5readatt: reading the given float Attribute failed");
	return retval;
      }
      for (size_t n = 0; n < numVal; ++n)
	value[n] = f_value[n];
    }
    else if (H5Tget_size(type) == sizeof(double))
    {
      if (H5Aread(attr, H5T_NATIVE_DOUBLE, value)<0)
      {
	H5Aclose(attr);
	H5Oclose(obj);
	H5Fclose(file_id);
	error("h5readatt: reading the given double Attribute failed");
	return retval;
      }
    }
    else
    {
      H5Aclose(attr);
      H5Oclose(obj);
      H5Fclose(file_id);
      error("h5readatt: reading the given float Attribute failed: cannot handle size of type");
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
      if (H5Aread(attr, H5T_NATIVE_INT, f_value)<0)
      {
	H5Aclose(attr);
	H5Oclose(obj);
	H5Fclose(file_id);
	error("h5readatt: reading the given integer Attribute failed");
	return retval;
      }
      for (size_t n=0;n<numVal;++n)
	value[n] = f_value[n]*1.0;
    }
    else
    {
      H5Aclose(attr);
      H5Oclose(obj);
      H5Fclose(file_id);
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
    if (H5Aread(attr, type, buf)<0)
    {
      H5Aclose(attr);
      H5Oclose(obj);
      H5Fclose(file_id);
      error("h5readatt: reading the given string Attribute failed");
      return retval;
    }
    retval = octave_value(buf);
  }
  else //none of the supported data types
  {
    H5Aclose(attr);
    H5Oclose(obj);
    H5Fclose(file_id);
    error("h5readatt: attribute type not supported");
    return retval;
  }
  

  H5Aclose(attr);
  H5Oclose(obj);
  H5Fclose(file_id);

  return retval;
#endif
}

