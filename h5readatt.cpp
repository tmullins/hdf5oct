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

#include <octave/oct.h>
#include <octave/ov-struct.h>
#include <cstdlib>
#include <string>
#include "hdf5.h"
#include "h5readatt.doc.h"
using namespace std;

// this special treatment is necessary because Win32-Octave ships with a very old hdf5 version (1.6.10)
void CloseH5Object(hid_t obj)
{
#if ((H5_VERS_MAJOR == 1) && (H5_VERS_MINOR == 6))
	// try group close, then Dataset close
	if (H5Gclose(obj)<0)
		H5Dclose(obj);
#else
	H5Oclose(obj);
#endif
}

DEFUN_DLD (h5readatt, args, nargout, string((char*) h5readatt_doc))
{
	octave_value retval;
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

	//suppress hdf5 error output
	H5Eset_auto(H5E_DEFAULT,0,0);

	//open the hdf5 file
	hid_t file = H5Fopen( args(0).string_value().c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );
	
	if (file==-1)
	{
	  error("h5readatt: opening the given File failed");
	  return retval;
	}

#if ((H5_VERS_MAJOR == 1) && (H5_VERS_MINOR == 6))
	// this special treatment is necessary because Win32-Octave ships with a very old hdf5 version (1.6.10)
	hid_t obj = -1;
	//try opening the group
	obj = H5Gopen(file, args(1).string_value().c_str());
	//try opening the dataset if group failed
	if (obj==-1)
		obj = H5Dopen(file, args(1).string_value().c_str());
#else
	hid_t obj = H5Oopen(file, args(1).string_value().c_str(), H5P_DEFAULT);
#endif

	if (obj==-1)
	{
	  CloseH5Object(obj);
	  H5Fclose(file);
	  error("h5readatt: opening the given Object failed");
	  return retval;
	}

	hid_t attr = H5Aopen_name(obj, args(2).string_value().c_str());
	if (attr==-1)
	{
	  CloseH5Object(obj);
	  H5Fclose(file);
	  error("h5readatt: opening the given Attribute failed");
	  return retval;
	}

	hid_t type = H5Aget_type(attr);
	if (type<0)
	{
	  H5Aclose(attr);
	  CloseH5Object(obj);
	  H5Fclose(file);
	  error("h5readatt: dataset type error");
	  return retval;
	}

	size_t numVal = H5Aget_storage_size(attr)/H5Tget_size(type);
	if(H5Tget_class(type)==H5T_FLOAT)
	  {
	    double value[numVal];
	    if (H5Tget_size(type)==sizeof(float))
	      {
		float f_value[numVal];
		if (H5Aread(attr, H5T_NATIVE_FLOAT, f_value)<0)
		  {
		    H5Aclose(attr);
		    CloseH5Object(obj);
		    H5Fclose(file);
		    error("h5readatt: reading the given float Attribute failed");
		    return retval;
		  }
		for (size_t n=0;n<numVal;++n)
		  value[n] = f_value[n];
	      }
	    else if (H5Tget_size(type)==sizeof(double))
	      {
		if (H5Aread(attr, H5T_NATIVE_DOUBLE, value)<0)
		  {
		    H5Aclose(attr);
		    CloseH5Object(obj);
		    H5Fclose(file);
		    error("h5readatt: reading the given double Attribute failed");
		    return retval;
		  }
	      }
	    else
	      {
		H5Aclose(attr);
		CloseH5Object(obj);
		H5Fclose(file);
		error("h5readatt: reading the given float Attribute failed: cannot handle size of type");
		return retval;
	      }

	    Matrix mat(numVal,1);
	    for (size_t n=0;n<numVal;++n)
	      mat(n)=value[n];
	    retval = octave_value(mat);

	  }
	else if(H5Tget_class(type)==H5T_INTEGER)
	  {
	    // Integer attributes are casted to floating point octave values
	    
	    double value[numVal];
	    if (H5Tget_size(type)==sizeof(int))
	      {
		int f_value[numVal];
		if(H5Aread(attr, H5T_NATIVE_INT, f_value)<0)
		  {
		    H5Aclose(attr);
		    CloseH5Object(obj);
		    H5Fclose(file);
		    error("h5readatt: reading the given integer Attribute failed");
		    return retval;
		  }
		for (size_t n=0;n<numVal;++n)
		  value[n] = f_value[n]*1.0;
	      }
	    else
	      {
		H5Aclose(attr);
		CloseH5Object(obj);
		H5Fclose(file);
		error("h5readatt: reading the given integer Attribute failed: cannot handle size of type");
		return retval;
	      }

	    Matrix mat(numVal,1);
	    for (size_t n=0;n<numVal;++n)
	      mat(n)=value[n];
	    retval = octave_value(mat);

	  }
	else if(H5Tget_class(type)==H5T_STRING)
	  {
	    // Size of each string:
	    size_t size = H5Tget_size(type);
	    // to read an array of strings (for future work): 
	    //totsize = size*sdim[0]*sdim[1];
	    // to read a single string:
	    size_t totsize = size;
	    // Set up read buffer for attribute
	    char* buf = (char*)calloc(totsize, sizeof(char));
	    if(H5Aread(attr, type, buf)<0)
	    {
	      H5Aclose(attr);
	      CloseH5Object(obj);
	      H5Fclose(file);
	      error("h5readatt: reading the given string Attribute failed");
	      return retval;
	    }
	    retval = octave_value(buf);
	  }
	else //none of the supported data types
	{
	  H5Aclose(attr);
	  CloseH5Object(obj);
	  H5Fclose(file);
	  error("h5readatt: attribute type not supported");
	  return retval;
	}


	H5Aclose(attr);
	CloseH5Object(obj);
	H5Fclose(file);

	return retval;
}

