CXX=h5c++
#CXX=g++
CPPFLAGS=
octs=h5read.oct h5readatt.oct
objs=$(octs:.oct=.o) h5file.o
docs=$(octs:.oct=.doc)
hdrs=$(docs:.doc=.doc.h)
MKOCTFILE=CXX=$(CXX) mkoctfile

INSTALLDIR=~/octave/hdf5oct/

all: $(octs) test.h5

h5file.o:
h5read.o: h5read.doc.h
h5readatt.o: h5readatt.doc.h
h5write.o: h5write.doc.h

%.oct: %.o h5file.o
	$(MKOCTFILE) -o $@ $< h5file.o -lhdf5 -DH5_USE_16_API

%.o: %.cpp h5file.h
	$(MKOCTFILE) -c $<

%.doc.h: %.doc
	xxd -i $< | sed 's/\(0x..\)$$/\1, 0x00/' > $@

clean:
	rm -f $(objs) $(octs) $(hdrs) write_test test.h5

test.h5: write_test
	./write_test

write_test: write_test.cpp
	$(CXX) -o $@ $< -lhdf5

install: $(octs)
	mkdir -p $(INSTALLDIR)
	cp $(octs) $(INSTALLDIR)

uninstall:
	rm $(addprefix $(INSTALLDIR),$(octs))
	rmdir $(INSTALLDIR)
