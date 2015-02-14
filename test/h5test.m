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

pkg load hdf5oct

%system (sprintf ("gnome-terminal --command 'gdb -p %d'", getpid ()), "async");

sleep(6);

% generate testdata %%%%%%%%%%%%%%%%%%%%%%%%%mm
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

disp("Test h5read...")
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


disp("Test h5write...")
disp("1D..")
h5write("test.h5", '/foo1', 0.1:0.1:1)
## disp("foo1b")
## h5write("test.h5", '/foo1b', (0.1:0.1:1)*2)
disp("2D..")
h5write("test.h5", '/foo2', reshape(0.1:0.1:10,[10 10]))
disp("3D..")
h5write("test.h5", '/foo3', reshape(0.1:0.1:100,[10 10 10]))
disp("4D..")
h5write("test.h5", '/foo4', reshape(0.1:0.1:1000,[10 10 10 10]))
s=5
h5write("test.h5", '/foo1_int', 1:s**1)
h5write("test.h5", '/foo2_int', reshape(1:s**2, [s s]))
h5write("test.h5", '/foo3_int', reshape(1:s**3, [s s s ]))
h5write("test.h5", '/foo4_int', reshape(1:s**4, [s s s s]))
				       

disp("Test h5writeatt...")
h5writeatt("test.h5", '/', 'testatt_double',12.34)
%h5writeatt("test.h5", '/', 'testatt2_double',0.1:0.1:0.5)
h5writeatt("test.h5", '/', 'testatt_int',7)
%h5writeatt("test.h5", '/', 'testatt_string','hallo')

h5writeatt("test.h5", '/foo1', 'testatt_double',12.34)
h5writeatt("test.h5", '/foo1', 'testatt_int',7)
%h5writeatt("test.h5", '/foo1', 'testatt_string','hallo')

disp("Test h5readatt...")
h5readatt("test.h5", '/', 'testatt_double')
h5readatt("test.h5", '/', 'testatt_int')
%h5readatt("test.h5", '/', 'testatt_string')

h5readatt("test.h5", '/foo1', 'testatt_double')
h5readatt("test.h5", '/foo1', 'testatt_int')
%h5readatt("test.h5", '/foo1', 'testatt_string')
