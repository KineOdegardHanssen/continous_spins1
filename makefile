include Makefile.local


SOURCES =  gaussiandeviate.cpp bond.cpp site.cpp lattice.cpp gaussiandeviate.h site.h lattice.h makefile Makefile.local

MYPROGSOURCES = main.cpp makefile Makefile.local
MONTECARLOSOURCES = montecarlo.cpp montecarlo.h makefile Makefile.local
BONDSOURCES = bond.cpp bond.h makefile Makefile.local
SITESOURCES = site.cpp site.h makefile Makefile.local
LATTICESOURCES = lattice.cpp lattice.h makefile Makefile.local
GAUSSIANDEVIATESOURCES = gaussiandeviate.cpp gaussiandeviate.h makefile Makefile.local

MYOBJECTS = myprog.o montecarlo.o bond.o site.o lattice.o gaussiandeviate.o


#%.exec : %.o 
#	$(CCC) $(LINKOPTS) -o $@ $^ $(LIBDIR) $(LIBS)

myprog.exec : $(MYOBJECTS)
	$(CCC) $(LINKOPTS) -o $@ $^ $(LIBDIR) $(LIBS)


myprog.o	: $(MYPROGSOURCES)
	$(CCC) $(CCFLAGS) -c -o $@ $<

montecarlo.o	: $(MONTECARLOSOURCES)
	$(CCC) $(CCFLAGS) -c -o $@ $<
	
bond.o	: $(BONDSOURCES)
	$(CCC) $(CCFLAGS) -c -o $@ $<

site.o	: $(SITESOURCES)
	$(CCC) $(CCFLAGS) -c -o $@ $<

lattice.o	: $(LATTICESOURCES)
	$(CCC) $(CCFLAGS) -c -o $@ $<

gaussiandeviate.o	: $(GAUSSIANDEVIATESOURCES)
	$(CCC) $(CCFLAGS) -c -o $@ $<





spotless:	
	make clean



clean	:
	rm -f *.dvi
	rm -f *.aux
	rm -f *.log
	rm -f *~
	rm -f core
	rm -f *.o
	rm -f *.exec
	rm -f *.d








