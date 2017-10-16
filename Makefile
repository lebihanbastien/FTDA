#------------------------------------------------------------------------------#
# This makefile was generated by 'cbp2make' tool rev.147                       #
#------------------------------------------------------------------------------#


WORKDIR = `pwd`

CC = gcc
CXX = g++
AR = ar
LD = g++
WINDRES = windres

INC = -Isrc -Iinclude -Itest -I../OOFTDA -I/usr/local/include/gsl -Ifortran
CFLAGS = -Wall
RESINC = 
LIBDIR = -Lfortran -L../OOFTDA
LIB = 
LDFLAGS = fortran/routines0.o -lgfortran fortran/foun.o -lgfortran -fopenmp

INC_DEBUG = $(INC)
CFLAGS_DEBUG = $(CFLAGS) -pg -g -fopenmp
RESINC_DEBUG = $(RESINC)
RCFLAGS_DEBUG = $(RCFLAGS)
LIBDIR_DEBUG = $(LIBDIR)
LIB_DEBUG = $(LIB)
LDFLAGS_DEBUG = $(LDFLAGS) -pg -fopenmp
OBJDIR_DEBUG = obj/Debug
DEP_DEBUG = 
OUT_DEBUG = bin/Debug/OOFTDA

INC_RELEASE = $(INC)
CFLAGS_RELEASE = $(CFLAGS) -O2 -fopenmp
RESINC_RELEASE = $(RESINC)
RCFLAGS_RELEASE = $(RCFLAGS)
LIBDIR_RELEASE = $(LIBDIR)
LIB_RELEASE = $(LIB)
LDFLAGS_RELEASE = $(LDFLAGS) -s -fopenmp
OBJDIR_RELEASE = obj/Release
DEP_RELEASE = 
OUT_RELEASE = bin/Release/OOFTDA

OBJ_DEBUG = $(OBJDIR_DEBUG)/src/lpdyneq.o $(OBJDIR_DEBUG)/src/odezero.o $(OBJDIR_DEBUG)/src/ode.o $(OBJDIR_DEBUG)/src/nrutil.o $(OBJDIR_DEBUG)/src/nfo2.o $(OBJDIR_DEBUG)/src/multimin_test.o $(OBJDIR_DEBUG)/src/multimin.o $(OBJDIR_DEBUG)/src/init.o $(OBJDIR_DEBUG)/src/gslc.o $(OBJDIR_DEBUG)/src/gnuplot_i.o $(OBJDIR_DEBUG)/src/ftda.o $(OBJDIR_DEBUG)/src/ertbp.o $(OBJDIR_DEBUG)/src/eminsem.o $(OBJDIR_DEBUG)/src/diffcorr.o $(OBJDIR_DEBUG)/src/define_env.o $(OBJDIR_DEBUG)/src/qbcp.o $(OBJDIR_DEBUG)/test/oftsh_test.o $(OBJDIR_DEBUG)/test/ofts_test.o $(OBJDIR_DEBUG)/test/ofs_test.o $(OBJDIR_DEBUG)/src/timec.o $(OBJDIR_DEBUG)/src/qbtbp.o $(OBJDIR_DEBUG)/src/coc.o $(OBJDIR_DEBUG)/src/poincare.o $(OBJDIR_DEBUG)/src/pmt.o $(OBJDIR_DEBUG)/src/pmode.o $(OBJDIR_DEBUG)/src/pmcoc.o $(OBJDIR_DEBUG)/src/pm.o $(OBJDIR_DEBUG)/src/parameters.o $(OBJDIR_DEBUG)/src/Constants.o $(OBJDIR_DEBUG)/src/Config.o $(OBJDIR_DEBUG)/main.o

OBJ_RELEASE = $(OBJDIR_RELEASE)/src/lpdyneq.o $(OBJDIR_RELEASE)/src/odezero.o $(OBJDIR_RELEASE)/src/ode.o $(OBJDIR_RELEASE)/src/nrutil.o $(OBJDIR_RELEASE)/src/nfo2.o $(OBJDIR_RELEASE)/src/multimin_test.o $(OBJDIR_RELEASE)/src/multimin.o $(OBJDIR_RELEASE)/src/init.o $(OBJDIR_RELEASE)/src/gslc.o $(OBJDIR_RELEASE)/src/gnuplot_i.o $(OBJDIR_RELEASE)/src/ftda.o $(OBJDIR_RELEASE)/src/ertbp.o $(OBJDIR_RELEASE)/src/eminsem.o $(OBJDIR_RELEASE)/src/diffcorr.o $(OBJDIR_RELEASE)/src/define_env.o $(OBJDIR_RELEASE)/src/qbcp.o $(OBJDIR_RELEASE)/test/oftsh_test.o $(OBJDIR_RELEASE)/test/ofts_test.o $(OBJDIR_RELEASE)/test/ofs_test.o $(OBJDIR_RELEASE)/src/timec.o $(OBJDIR_RELEASE)/src/qbtbp.o $(OBJDIR_RELEASE)/src/coc.o $(OBJDIR_RELEASE)/src/poincare.o $(OBJDIR_RELEASE)/src/pmt.o $(OBJDIR_RELEASE)/src/pmode.o $(OBJDIR_RELEASE)/src/pmcoc.o $(OBJDIR_RELEASE)/src/pm.o $(OBJDIR_RELEASE)/src/parameters.o $(OBJDIR_RELEASE)/src/Constants.o $(OBJDIR_RELEASE)/src/Config.o $(OBJDIR_RELEASE)/main.o

all: debug release

clean: clean_debug clean_release

before_debug: 
	test -d bin/Debug || mkdir -p bin/Debug
	test -d $(OBJDIR_DEBUG)/src || mkdir -p $(OBJDIR_DEBUG)/src
	test -d $(OBJDIR_DEBUG)/test || mkdir -p $(OBJDIR_DEBUG)/test
	test -d $(OBJDIR_DEBUG) || mkdir -p $(OBJDIR_DEBUG)

after_debug: 

debug: before_debug out_debug after_debug

out_debug: before_debug $(OBJ_DEBUG) $(DEP_DEBUG)
	$(LD) $(LIBDIR_DEBUG) -o $(OUT_DEBUG) $(OBJ_DEBUG)  $(LDFLAGS_DEBUG) $(LIB_DEBUG)

$(OBJDIR_DEBUG)/src/lpdyneq.o: src/lpdyneq.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c src/lpdyneq.cpp -o $(OBJDIR_DEBUG)/src/lpdyneq.o

$(OBJDIR_DEBUG)/src/odezero.o: src/odezero.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c src/odezero.cpp -o $(OBJDIR_DEBUG)/src/odezero.o

$(OBJDIR_DEBUG)/src/ode.o: src/ode.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c src/ode.cpp -o $(OBJDIR_DEBUG)/src/ode.o

$(OBJDIR_DEBUG)/src/nrutil.o: src/nrutil.c
	$(CC) $(CFLAGS_DEBUG) $(INC_DEBUG) -c src/nrutil.c -o $(OBJDIR_DEBUG)/src/nrutil.o

$(OBJDIR_DEBUG)/src/nfo2.o: src/nfo2.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c src/nfo2.cpp -o $(OBJDIR_DEBUG)/src/nfo2.o

$(OBJDIR_DEBUG)/src/multimin_test.o: src/multimin_test.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c src/multimin_test.cpp -o $(OBJDIR_DEBUG)/src/multimin_test.o

$(OBJDIR_DEBUG)/src/multimin.o: src/multimin.c
	$(CC) $(CFLAGS_DEBUG) $(INC_DEBUG) -c src/multimin.c -o $(OBJDIR_DEBUG)/src/multimin.o

$(OBJDIR_DEBUG)/src/init.o: src/init.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c src/init.cpp -o $(OBJDIR_DEBUG)/src/init.o

$(OBJDIR_DEBUG)/src/gslc.o: src/gslc.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c src/gslc.cpp -o $(OBJDIR_DEBUG)/src/gslc.o

$(OBJDIR_DEBUG)/src/gnuplot_i.o: src/gnuplot_i.c
	$(CC) $(CFLAGS_DEBUG) $(INC_DEBUG) -c src/gnuplot_i.c -o $(OBJDIR_DEBUG)/src/gnuplot_i.o

$(OBJDIR_DEBUG)/src/ftda.o: src/ftda.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c src/ftda.cpp -o $(OBJDIR_DEBUG)/src/ftda.o

$(OBJDIR_DEBUG)/src/ertbp.o: src/ertbp.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c src/ertbp.cpp -o $(OBJDIR_DEBUG)/src/ertbp.o

$(OBJDIR_DEBUG)/src/eminsem.o: src/eminsem.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c src/eminsem.cpp -o $(OBJDIR_DEBUG)/src/eminsem.o

$(OBJDIR_DEBUG)/src/diffcorr.o: src/diffcorr.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c src/diffcorr.cpp -o $(OBJDIR_DEBUG)/src/diffcorr.o

$(OBJDIR_DEBUG)/src/define_env.o: src/define_env.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c src/define_env.cpp -o $(OBJDIR_DEBUG)/src/define_env.o

$(OBJDIR_DEBUG)/src/qbcp.o: src/qbcp.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c src/qbcp.cpp -o $(OBJDIR_DEBUG)/src/qbcp.o

$(OBJDIR_DEBUG)/test/oftsh_test.o: test/oftsh_test.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c test/oftsh_test.cpp -o $(OBJDIR_DEBUG)/test/oftsh_test.o

$(OBJDIR_DEBUG)/test/ofts_test.o: test/ofts_test.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c test/ofts_test.cpp -o $(OBJDIR_DEBUG)/test/ofts_test.o

$(OBJDIR_DEBUG)/test/ofs_test.o: test/ofs_test.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c test/ofs_test.cpp -o $(OBJDIR_DEBUG)/test/ofs_test.o

$(OBJDIR_DEBUG)/src/timec.o: src/timec.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c src/timec.cpp -o $(OBJDIR_DEBUG)/src/timec.o

$(OBJDIR_DEBUG)/src/qbtbp.o: src/qbtbp.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c src/qbtbp.cpp -o $(OBJDIR_DEBUG)/src/qbtbp.o

$(OBJDIR_DEBUG)/src/coc.o: src/coc.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c src/coc.cpp -o $(OBJDIR_DEBUG)/src/coc.o

$(OBJDIR_DEBUG)/src/poincare.o: src/poincare.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c src/poincare.cpp -o $(OBJDIR_DEBUG)/src/poincare.o

$(OBJDIR_DEBUG)/src/pmt.o: src/pmt.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c src/pmt.cpp -o $(OBJDIR_DEBUG)/src/pmt.o

$(OBJDIR_DEBUG)/src/pmode.o: src/pmode.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c src/pmode.cpp -o $(OBJDIR_DEBUG)/src/pmode.o

$(OBJDIR_DEBUG)/src/pmcoc.o: src/pmcoc.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c src/pmcoc.cpp -o $(OBJDIR_DEBUG)/src/pmcoc.o

$(OBJDIR_DEBUG)/src/pm.o: src/pm.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c src/pm.cpp -o $(OBJDIR_DEBUG)/src/pm.o

$(OBJDIR_DEBUG)/src/parameters.o: src/parameters.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c src/parameters.cpp -o $(OBJDIR_DEBUG)/src/parameters.o

$(OBJDIR_DEBUG)/src/Constants.o: src/Constants.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c src/Constants.cpp -o $(OBJDIR_DEBUG)/src/Constants.o

$(OBJDIR_DEBUG)/src/Config.o: src/Config.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c src/Config.cpp -o $(OBJDIR_DEBUG)/src/Config.o

$(OBJDIR_DEBUG)/main.o: main.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c main.cpp -o $(OBJDIR_DEBUG)/main.o

clean_debug: 
	rm -f $(OBJ_DEBUG) $(OUT_DEBUG)
	rm -rf bin/Debug
	rm -rf $(OBJDIR_DEBUG)/src
	rm -rf $(OBJDIR_DEBUG)/test
	rm -rf $(OBJDIR_DEBUG)

before_release: 
	test -d bin/Release || mkdir -p bin/Release
	test -d $(OBJDIR_RELEASE)/src || mkdir -p $(OBJDIR_RELEASE)/src
	test -d $(OBJDIR_RELEASE)/test || mkdir -p $(OBJDIR_RELEASE)/test
	test -d $(OBJDIR_RELEASE) || mkdir -p $(OBJDIR_RELEASE)

after_release: 

release: before_release out_release after_release

out_release: before_release $(OBJ_RELEASE) $(DEP_RELEASE)
	$(LD) $(LIBDIR_RELEASE) -o $(OUT_RELEASE) $(OBJ_RELEASE)  $(LDFLAGS_RELEASE) $(LIB_RELEASE)

$(OBJDIR_RELEASE)/src/lpdyneq.o: src/lpdyneq.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c src/lpdyneq.cpp -o $(OBJDIR_RELEASE)/src/lpdyneq.o

$(OBJDIR_RELEASE)/src/odezero.o: src/odezero.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c src/odezero.cpp -o $(OBJDIR_RELEASE)/src/odezero.o

$(OBJDIR_RELEASE)/src/ode.o: src/ode.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c src/ode.cpp -o $(OBJDIR_RELEASE)/src/ode.o

$(OBJDIR_RELEASE)/src/nrutil.o: src/nrutil.c
	$(CC) $(CFLAGS_RELEASE) $(INC_RELEASE) -c src/nrutil.c -o $(OBJDIR_RELEASE)/src/nrutil.o

$(OBJDIR_RELEASE)/src/nfo2.o: src/nfo2.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c src/nfo2.cpp -o $(OBJDIR_RELEASE)/src/nfo2.o

$(OBJDIR_RELEASE)/src/multimin_test.o: src/multimin_test.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c src/multimin_test.cpp -o $(OBJDIR_RELEASE)/src/multimin_test.o

$(OBJDIR_RELEASE)/src/multimin.o: src/multimin.c
	$(CC) $(CFLAGS_RELEASE) $(INC_RELEASE) -c src/multimin.c -o $(OBJDIR_RELEASE)/src/multimin.o

$(OBJDIR_RELEASE)/src/init.o: src/init.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c src/init.cpp -o $(OBJDIR_RELEASE)/src/init.o

$(OBJDIR_RELEASE)/src/gslc.o: src/gslc.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c src/gslc.cpp -o $(OBJDIR_RELEASE)/src/gslc.o

$(OBJDIR_RELEASE)/src/gnuplot_i.o: src/gnuplot_i.c
	$(CC) $(CFLAGS_RELEASE) $(INC_RELEASE) -c src/gnuplot_i.c -o $(OBJDIR_RELEASE)/src/gnuplot_i.o

$(OBJDIR_RELEASE)/src/ftda.o: src/ftda.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c src/ftda.cpp -o $(OBJDIR_RELEASE)/src/ftda.o

$(OBJDIR_RELEASE)/src/ertbp.o: src/ertbp.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c src/ertbp.cpp -o $(OBJDIR_RELEASE)/src/ertbp.o

$(OBJDIR_RELEASE)/src/eminsem.o: src/eminsem.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c src/eminsem.cpp -o $(OBJDIR_RELEASE)/src/eminsem.o

$(OBJDIR_RELEASE)/src/diffcorr.o: src/diffcorr.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c src/diffcorr.cpp -o $(OBJDIR_RELEASE)/src/diffcorr.o

$(OBJDIR_RELEASE)/src/define_env.o: src/define_env.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c src/define_env.cpp -o $(OBJDIR_RELEASE)/src/define_env.o

$(OBJDIR_RELEASE)/src/qbcp.o: src/qbcp.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c src/qbcp.cpp -o $(OBJDIR_RELEASE)/src/qbcp.o

$(OBJDIR_RELEASE)/test/oftsh_test.o: test/oftsh_test.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c test/oftsh_test.cpp -o $(OBJDIR_RELEASE)/test/oftsh_test.o

$(OBJDIR_RELEASE)/test/ofts_test.o: test/ofts_test.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c test/ofts_test.cpp -o $(OBJDIR_RELEASE)/test/ofts_test.o

$(OBJDIR_RELEASE)/test/ofs_test.o: test/ofs_test.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c test/ofs_test.cpp -o $(OBJDIR_RELEASE)/test/ofs_test.o

$(OBJDIR_RELEASE)/src/timec.o: src/timec.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c src/timec.cpp -o $(OBJDIR_RELEASE)/src/timec.o

$(OBJDIR_RELEASE)/src/qbtbp.o: src/qbtbp.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c src/qbtbp.cpp -o $(OBJDIR_RELEASE)/src/qbtbp.o

$(OBJDIR_RELEASE)/src/coc.o: src/coc.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c src/coc.cpp -o $(OBJDIR_RELEASE)/src/coc.o

$(OBJDIR_RELEASE)/src/poincare.o: src/poincare.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c src/poincare.cpp -o $(OBJDIR_RELEASE)/src/poincare.o

$(OBJDIR_RELEASE)/src/pmt.o: src/pmt.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c src/pmt.cpp -o $(OBJDIR_RELEASE)/src/pmt.o

$(OBJDIR_RELEASE)/src/pmode.o: src/pmode.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c src/pmode.cpp -o $(OBJDIR_RELEASE)/src/pmode.o

$(OBJDIR_RELEASE)/src/pmcoc.o: src/pmcoc.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c src/pmcoc.cpp -o $(OBJDIR_RELEASE)/src/pmcoc.o

$(OBJDIR_RELEASE)/src/pm.o: src/pm.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c src/pm.cpp -o $(OBJDIR_RELEASE)/src/pm.o

$(OBJDIR_RELEASE)/src/parameters.o: src/parameters.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c src/parameters.cpp -o $(OBJDIR_RELEASE)/src/parameters.o

$(OBJDIR_RELEASE)/src/Constants.o: src/Constants.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c src/Constants.cpp -o $(OBJDIR_RELEASE)/src/Constants.o

$(OBJDIR_RELEASE)/src/Config.o: src/Config.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c src/Config.cpp -o $(OBJDIR_RELEASE)/src/Config.o

$(OBJDIR_RELEASE)/main.o: main.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c main.cpp -o $(OBJDIR_RELEASE)/main.o

clean_release: 
	rm -f $(OBJ_RELEASE) $(OUT_RELEASE)
	rm -rf bin/Release
	rm -rf $(OBJDIR_RELEASE)/src
	rm -rf $(OBJDIR_RELEASE)/test
	rm -rf $(OBJDIR_RELEASE)

.PHONY: before_debug after_debug clean_debug before_release after_release clean_release

