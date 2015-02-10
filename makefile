CXX=h5c++
#CXX=g++
CPPFLAGS=
octs=h5read.oct h5readatt.oct
objs=$(octs:.oct=.o) h5file.o
docs=$(octs:.oct=.doc)
hdrs=$(docs:.doc=.doc.h)
MKOCTFILE=CXX=$(CXX) mkoctfile

VERSION=0.2.0
PACKAGEFILE=hdf5oct-$(VERSION).tar.gz


all: $(octs) test.h5 package

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
	rm -f $(objs) $(octs) $(hdrs) write_test test.h5 $(PACKAGEFILE)


install: $(PACKAGEFILE)
	@echo "-- Install Octave Package ------------"
	octave --silent --no-gui --eval "pkg install $(PACKAGEFILE)"

uninstall:
	@echo "-- Uninstall Octave Package ----------"
	octave --silent --no-gui --eval "pkg uninstall hdf5oct"

package: $(PACKAGEFILE)

$(PACKAGEFILE): $(octs)
	mkdir -p package/inst
	cp $(octs) package/inst
	tar -czf $(PACKAGEFILE) package/

# TESTING ###########

# a minimal program to generate some testdata
write_test: write_test.cpp
	$(CXX) -o $@ $< -lhdf5

test.h5: write_test
	./write_test

# a target to test the octave functions
test: install test.h5
	@echo "-- Perform Tests --------------"
	octave --silent --no-gui h5test.m
