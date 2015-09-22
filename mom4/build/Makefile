# Makefile created by mkmf $Id: mkmf,v 18.0 2010/03/02 23:26:08 fms Exp $ 



include template


.DEFAULT:
	-echo $@ does not exist.
all: letkf.MAKE.056
SFMT.o: ../../common/SFMT.f90
	ftn  -O3 -parallel -what   -c	../../common/SFMT.f90
common.o: ../../common/common.f90
	ftn  -O3 -parallel -what   -c	../../common/common.f90
common_letkf.o: ../../common/common_letkf.f90 common.o common_mtx.o params_letkf.o
	ftn  -O3 -parallel -what  -c	../../common/common_letkf.f90
common_mom4.o: ../common/common_mom4.f90 common.o params_model.o vars_model.o params_letkf.o
	ftn  -O3 -parallel -what   -I/opt/cray/netcdf/4.2.0/intel/120/include -c	../common/common_mom4.f90
common_mpi.o: ../../common/common_mpi.f90
	ftn  -O3 -parallel -what  -c	../../common/common_mpi.f90
common_mpi_mom4.o: ../common/common_mpi_mom4.f90 common.o common_mpi.o common_mom4.o params_letkf.o
	ftn  -O3 -parallel -what  -c -I/opt/cray/netcdf/4.2.0/intel/120/include	../common/common_mpi_mom4.f90
common_mtx.o: ../../common/common_mtx.f90 common.o
	ftn  -O3 -parallel -what   -c	../../common/common_mtx.f90
common_obs_mom4.o: ../common/common_obs_mom4.f90 common.o common_mom4.o params_obs.o
	ftn  -O3 -parallel -what  -c	../common/common_obs_mom4.f90
letkf.o: ../letkf/letkf.f90 common.o common_mpi.o common_mom4.o common_mpi_mom4.o common_letkf.o letkf_obs.o letkf_tools.o params_letkf.o params_model.o params_obs.o
	ftn  -O3 -parallel -what  -c	../letkf/letkf.f90
letkf_local.o: ../letkf/letkf_local.f90 common.o common_mpi.o common_mom4.o common_mpi_mom4.o common_letkf.o letkf_obs.o params_letkf.o vars_letkf.o
	ftn  -O3 -parallel -what  -c	../letkf/letkf_local.f90
letkf_obs.o: ../letkf/letkf_obs.f90 common.o common_mpi.o common_mom4.o common_obs_mom4.o common_mpi_mom4.o common_letkf.o params_letkf.o params_obs.o vars_obs.o
	ftn  -O3 -parallel -what  -c	../letkf/letkf_obs.f90
letkf_tools.o: ../letkf/letkf_tools.f90 common.o common_mpi.o common_mom4.o common_mpi_mom4.o common_letkf.o letkf_obs.o letkf_local.o params_letkf.o vars_letkf.o
	ftn  -O3 -parallel -what  -c	../letkf/letkf_tools.f90
netlib2.o: ./netlib2.f
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./netlib2.f
params_letkf.o: ../letkf/params_letkf.f90 common.o
	ftn  -O3 -parallel -what   -c	../letkf/params_letkf.f90
params_model.o: ../letkf/params_model.f90
	ftn  -O3 -parallel -what  -c	../letkf/params_model.f90
params_obs.o: ../obs/params_obs.f90 common.o
	ftn  -O3 -parallel -what  -c	../obs/params_obs.f90
vars_letkf.o: ../letkf/vars_letkf.f90 common.o params_obs.o params_model.o
	ftn  -O3 -parallel -what  -c	../letkf/vars_letkf.f90
vars_model.o: ../letkf/vars_model.f90 common.o params_model.o
	ftn  -O3 -parallel -what  -c	../letkf/vars_model.f90
vars_obs.o: ../obs/vars_obs.f90 common.o
	ftn  -O3 -parallel -what  -c	../obs/vars_obs.f90
./common.f90: ../../common/common.f90
	cp ../../common/common.f90 .
./letkf_tools.f90: ../letkf/letkf_tools.f90
	cp ../letkf/letkf_tools.f90 .
./letkf_local.f90: ../letkf/letkf_local.f90
	cp ../letkf/letkf_local.f90 .
./common_mtx.f90: ../../common/common_mtx.f90
	cp ../../common/common_mtx.f90 .
./SFMT.f90: ../../common/SFMT.f90
	cp ../../common/SFMT.f90 .
./params_letkf.f90: ../letkf/params_letkf.f90
	cp ../letkf/params_letkf.f90 .
./letkf_obs.f90: ../letkf/letkf_obs.f90
	cp ../letkf/letkf_obs.f90 .
./vars_letkf.f90: ../letkf/vars_letkf.f90
	cp ../letkf/vars_letkf.f90 .
./common_mpi_mom4.f90: ../common/common_mpi_mom4.f90
	cp ../common/common_mpi_mom4.f90 .
./common_mpi.f90: ../../common/common_mpi.f90
	cp ../../common/common_mpi.f90 .
./common_obs_mom4.f90: ../common/common_obs_mom4.f90
	cp ../common/common_obs_mom4.f90 .
./vars_model.f90: ../letkf/vars_model.f90
	cp ../letkf/vars_model.f90 .
./params_obs.f90: ../obs/params_obs.f90
	cp ../obs/params_obs.f90 .
./vars_obs.f90: ../obs/vars_obs.f90
	cp ../obs/vars_obs.f90 .
./params_model.f90: ../letkf/params_model.f90
	cp ../letkf/params_model.f90 .
./common_mom4.f90: ../common/common_mom4.f90
	cp ../common/common_mom4.f90 .
./common_letkf.f90: ../../common/common_letkf.f90
	cp ../../common/common_letkf.f90 .
./letkf.f90: ../letkf/letkf.f90
	cp ../letkf/letkf.f90 .
SRC = ../letkf/letkf_obs.f90 ./netlib2.f ../letkf/letkf_local.f90 ../common/common_mpi_mom4.f90 ../letkf/vars_letkf.f90 ../common/common_mom4.f90 ../letkf/params_letkf.f90 ../common/common_obs_mom4.f90 ../obs/vars_obs.f90 ../obs/params_obs.f90 ../../common/common.f90 ../../common/common_mpi.f90 ../../common/common_mtx.f90 ../../common/common_letkf.f90 ../letkf/params_model.f90 ../letkf/vars_model.f90 ../letkf/letkf_tools.f90 ../letkf/letkf.f90 ../../common/SFMT.f90
OBJ = letkf_obs.o netlib2.o letkf_local.o common_mpi_mom4.o vars_letkf.o common_mom4.o params_letkf.o common_obs_mom4.o vars_obs.o params_obs.o common.o common_mpi.o common_mtx.o common_letkf.o params_model.o vars_model.o letkf_tools.o letkf.o SFMT.o
OFF = ../../common/common.f90 ../letkf/letkf_tools.f90 ../letkf/letkf_local.f90 ../../common/common_mtx.f90 ../../common/SFMT.f90 ../letkf/params_letkf.f90 ../letkf/letkf_obs.f90 ../letkf/vars_letkf.f90 ../common/common_mpi_mom4.f90 ../../common/common_mpi.f90 ../common/common_obs_mom4.f90 ../letkf/vars_model.f90 ../obs/params_obs.f90 ../obs/vars_obs.f90 ../letkf/params_model.f90 ../common/common_mom4.f90 ../../common/common_letkf.f90 ../letkf/letkf.f90
clean: neat
	-rm -f .letkf.MAKE.056.cppdefs $(OBJ) letkf.MAKE.056
neat:
	-rm -f $(TMPFILES)
localize: $(OFF)
	cp $(OFF) .
TAGS: $(SRC)
	etags $(SRC)
tags: $(SRC)
	ctags $(SRC)
letkf.MAKE.056: $(OBJ) 
	$(LD) $(OBJ) -o letkf.MAKE.056  $(LDFLAGS)
