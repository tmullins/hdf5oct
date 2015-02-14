/*
 *
 *    Copyright 2015 Stefan Großhauser
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

#include "h5writeatt.doc.h"

using namespace std;

#if ((H5_VERS_MAJOR > 1) || (H5_VERS_MINOR >= 8))
#define HAVE_HDF5_18 1
#endif

#if defined(HAVE_HDF5) && defined(HAVE_HDF5_18)
#include "h5file.h"
#endif

DEFUN_DLD(h5writeatt, args, nargout, string((char*) h5writeatt_doc))
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
    
  H5E_auto_t oef;
  void *olderr;
  H5Eget_auto(H5E_DEFAULT,&oef,&olderr);
  // suppress hdf5 error output
  H5Eset_auto(H5E_DEFAULT,0,0);
  //open the hdf5 file
  H5File file(filename.c_str());
  // restore old setting
  H5Eset_auto(H5E_DEFAULT,oef,olderr);

  file.write_att(location.c_str(), attname.c_str(),
  		 args(3));

  return octave_value_list();
#endif
}
