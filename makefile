CXX=h5c++
src= h5write.cpp h5writeatt.cpp h5read.cpp h5readatt.cpp
headers=h5file.h
octs=$(src:.cpp=.oct)
objs=$(src:.cpp=.o) h5file.o
testobjs=test/write_test.o
docs=$(src:.cpp=.doc)
hdrs=$(docs:.doc=.doc.h)
MKOCTFILE=CXX=$(CXX) mkoctfile -g

VERSION=0.2.0
PACKAGEFILE=hdf5oct-$(VERSION).tar.gz

.PHONY: test clean install uninstall package

all: $(octs) package

h5file.o:
h5read.o: h5read.doc.h
h5readatt.o: h5readatt.doc.h
h5write.o: h5write.doc.h
h5writeatt.o: h5writeatt.doc.h

%.oct: $(objs)
	$(MKOCTFILE) -o $@ $(objs)

%.o: %.cpp $(headers)
	$(MKOCTFILE) -c $<

%.doc.h: %.doc
	xxd -i $< | sed 's/\(0x..\)$$/\1, 0x00/' > $@

clean:
	rm -f *.o *.oct *.doc.h package/inst/* $(testobjs) write_test test.h5 $(PACKAGEFILE)

install: $(PACKAGEFILE)
	@echo "-- Install Octave Package ------------"
	octave --silent --no-gui --eval "pkg install $(PACKAGEFILE)"

uninstall:
	@echo "-- Uninstall Octave Package ----------"
	octave --silent --no-gui --eval "pkg uninstall hdf5oct"

package: $(PACKAGEFILE)

$(PACKAGEFILE): $(octs)
	@echo "-- Create Octave Package Archive ------------"
	mkdir -p package/inst
	cp *.oct package/inst
	tar -czf $(PACKAGEFILE) package/

# TESTING ###########

# a minimal program to generate some testdata
test/write_test: test/write_test.cpp
	$(CXX) -o $@ $< -lhdf5

# a target to test the octave functions
test: test/write_test
	@echo "-- Perform Tests --------------"
	rm test/test.h5
	cd test && ./write_test
	cd test && octave --silent --no-gui h5test.m
