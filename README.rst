=======================================
hdf5oct - a HDF5 wrapper for GNU Octave
=======================================

Copyright 2012 Tom Mullins
Copyright 2015 Tom Mullins, Anton Starikov, Thorsten Liebig, Stefan Großhauser
Copyright 2008-2013 Andrew Collette

.. contents::

This is a library for GNU Octave for reading hdf5 files. At the moment it
provides the following functions:

h5read: modeled after Matlab's h5read, which can read subsets of a
 	 dataset.  Octave's load function will attempt to read an
 	 entire dataset, which for very large datasets is
 	 undesired. h5util's h5read function does slightly better,
 	 reading only 2d slices of 3d datasets, but that's still
 	 fairly limiting. This exposes libhdf5's H5Sselect_hyperslab
 	 in a way which tries to be compatible with Matlab.

h5readatt: Most of this function was written by thliebig. It allows
           to read scalar HDF5 attributes of some types.

h5write: Write a matrix to a dataset. This
         will either overwrite an already existing dataset, or allow to
         append hyperslabs to existing datasets.

h5writeatt: Attach an attribute to an object.

h5create: Create a dataset and specify its extent dimensions,
          datatype and chunk size.

h5delete: Delete a group, dataset, or attribute.

Note that only few of the HDF5 datatypes are supported by each of the
functions hdf5oct at the moment, typically one or several of double,
integer and string.

INSTALLATION
============

To install, just use::

    make
    
This will produce a package file named "hdf5oct-*.tar.gz" .  Then
you may either install the package with::

    make install

or you may start GNU Octave and install the package manually (using
the correct file name) with the command::

    pkg install hdf5oct-0.2.0.tar.gz

This will put the *.oct files somewhere where Octave will find them.
You can try running::

    make test


DEINSTALLATION
==============

To uninstall the package you may want to use::

   make uninstall

TROUBLESHOOTING
===============
* If you get an error like ``mkoctfile: not found``, you need ``liboctave-dev``


TODO
====

- write h5info, h5disp

- support compression flags for h5create

- read string typed datasets

- read string-array typed attributes

- write more comprehensive tests instead of a few random choices. Also
  test for error conditions.
