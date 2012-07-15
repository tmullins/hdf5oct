%{

    Copyright 2012 Tom Mullins


    This file is part of hdf5oct.

    hdf5oct is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    hdf5oct is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with hdf5oct.  If not, see <http://www.gnu.org/licenses/>.


%}

dims = [4, 8, 3, 2];

global d;
d = {(1:dims(1))'};
for i = 2:4
  i_dims = dims(1:i);
  d{i} = zeros(i_dims);
  for j = 1:i
    perm = i_dims;
    perm(j) = 1;
    M = 10^(j-1) * (1:i_dims(j))';
    M = permute(M, [2:j, 1, j+1:i]);
    d{i} += repmat(M, perm);
  end
end

function mat = hyperslab(mat, rank, start, count, stride, block)
  start = [start, ones(1, 4-rank)];
  count = [count, ones(1, 4-rank)];
  stride = [stride, ones(1, 4-rank)];
  block = [block, ones(1, 4-rank)];
  idx = {};
  for i = 1:4
    idx{i} = [];
    for b = 0:(block(i)-1)
      idx{i} = [idx{i}, (start(i)+b) : stride(i) : (count(i)*stride(i)-1+start(i)+b)];
    end
    idx{i} = sort(idx{i});
  end
  mat = mat(idx{1}, idx{2}, idx{3}, idx{4});
end

function check_slice(filename, rank, start=[], count=[], stride=[], block=[])
  global d;
  dataset = strcat('t', num2str(rank));
  if start
    if block
      data = h5read(filename, dataset, start, count, stride, block);
    else
      block = ones(1, rank);
      if stride
        data = h5read(filename, dataset, start, count, stride);
      else
        stride = ones(1, rank);
        data = h5read(filename, dataset, start, count);
      end
    end
    comp = hyperslab(d{rank}, rank, start, count, stride, block);
  else
    data = h5read(filename, dataset);
    if rank == 0
      comp = 0;
    else
      comp = d{rank};
    end
  end
  res = data == comp;
  while ~isscalar(res)
    res = all(res);
  end
  if res
    printf('Success\n');
  else
    printf('Failed\n');
  end
end

check_slice("test.h5", 0);
check_slice("test.h5", 1);
check_slice("test.h5", 2);
check_slice("test.h5", 3);
check_slice("test.h5", 4);
check_slice("test.h5", 1, [2], [2], [2], [1]);
check_slice("test.h5", 2, [1, 2], [2, 1]);
check_slice("test.h5", 3, [1, 2, 1], [2, 1, 2], [2, 1, 1]);
check_slice("test.h5", 3, [2, 1, 1], [2, 2, 2], [2, 3, 1], [1, 2, 1]);
check_slice("test.h5", 4, [2, 1, 1, 2], [2, 2, 2, 1], [2, 3, 1, 3], [1, 2, 1, 1]);
