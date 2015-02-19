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

#include <hdf5.h>
#include <cstdlib>
#include <cstring>

#define FILENAME "test.h5"

void write_test(hid_t file, int rank)
{
  if (rank < 0 || rank > 4)
  {
    printf("Error: unimplemented rank %d\n", rank);
    return;
  }

  hsize_t dims[] = {4, 8, 3, 2};
  hsize_t idx[] = {1, 1, 1, 1};

  int size = 1;
  for (int i = 0; i < rank; i++)
  {
    size *= dims[i];
  }
  double* data = (double*)malloc(size * sizeof(double));
  for (int i = 0; i < size; i++)
  {
    data[i] = 0;
    int powten = 1;
    for (int d = 0; d < rank; d++)
    {
      data[i] += powten*idx[d];
      powten *= 10;
    }
    for (int d = rank-1; d >= 0; d--)
    {
      idx[d]++;
      if (idx[d] <= dims[d]) break;
      else idx[d] = 1;
    }
  }

  char setname[32];
  sprintf(setname, "t%d", rank);

  hid_t dataspace = H5Screate_simple(rank, dims, NULL);
  hid_t dataset = H5Dcreate(file, setname, H5T_NATIVE_DOUBLE, dataspace,
      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
  H5Dclose(dataset);
  H5Sclose(dataspace);

  free(data);
}

int main(int argc, char **argv)
{
  hid_t file = H5Fcreate(FILENAME, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  for (int i = 0; i < 5; i++)
  {
    write_test(file, i);
  }
  H5Fclose(file);
  return 0;
}
