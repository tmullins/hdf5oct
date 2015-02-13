CXX=h5c++
#CXX=g++
#CPPFLAGS=-DH5_USE_16_API
src=h5read.cpp h5readatt.cpp h5file.cpp
headers=h5file.h
octfile=hdf5oct.oct
objs=$(src:.cpp=.o)
testobjs=test/write_test.o
docs=$(src:.cpp=.doc)
hdrs=$(docs:.doc=.doc.h)
MKOCTFILE=CXX=$(CXX) mkoctfile

VERSION=0.2.0
PACKAGEFILE=hdf5oct-$(VERSION).tar.gz

.PHONY: test clean install uninstall package

all: $(octfile) package

h5file.o:
h5read.o: h5read.doc.h
h5readatt.o: h5readatt.doc.h
h5write.o: h5write.doc.h

$(octfile): $(objs)
	$(MKOCTFILE) -o $@ $(objs)

%.o: %.cpp $(headers)
	$(MKOCTFILE) -c $<

%.doc.h: %.doc
	xxd -i $< | sed 's/\(0x..\)$$/\1, 0x00/' > $@

clean:
	rm -f $(objs) $(octfile) $(hdrs) $(testobjs) write_test test.h5 $(PACKAGEFILE)

install: $(PACKAGEFILE)
	@echo "-- Install Octave Package ------------"
	octave --silent --no-gui --eval "pkg install $(PACKAGEFILE)"

uninstall:
	@echo "-- Uninstall Octave Package ----------"
	octave --silent --no-gui --eval "pkg uninstall hdf5oct"

package: $(PACKAGEFILE)

$(PACKAGEFILE): $(octfile)
	mkdir -p package/inst
	cp $(octfile) package/inst
	tar -czf $(PACKAGEFILE) package/

# TESTING ###########

# a minimal program to generate some testdata
write_test: test/write_test.cpp
	$(CXX) -o $@ $< -lhdf5

test/test.h5: test/write_test
	cd test && ./write_test

# a target to test the octave functions
test: test/test.h5
	@echo "-- Perform Tests --------------"
	cd test && octave --silent --no-gui h5test.m
