octs=h5read.oct
objs=$(octs:.oct=.o) h5file.o
docs=$(octs:.oct=.doc)
hdrs=$(docs:.doc=.doc.h)

all: $(octs) test.h5

h5file.o:
h5read.o: h5read.doc.h
h5write.o: h5write.doc.h

%.oct: %.o h5file.o
	mkoctfile -o $@ $< h5file.o -lhdf5

%.o: %.cpp h5file.h
	mkoctfile -c $<

%.doc.h: %.doc
	xxd -i $< | sed 's/\(0x..\)$$/\1, 0x00/' > $@

clean:
	rm -f $(objs) $(octs) $(hdrs) write_test test.h5

test.h5: write_test
	./write_test

write_test: write_test.cpp
	g++ -o $@ $< -lhdf5
