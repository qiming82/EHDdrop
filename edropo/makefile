FC=gfortran
#FFLAGS=-g
# edrop
#---------
#
OBJ = edrop.o arc_evl.o bess_I01K01.o df.o dgesv.o ell_int.o ell_int2.o gauss_leg.o gauss_log.o \
gsint.o lgf_ax_fs.o polyn.o regridf.o regridfadp.o\
seval.o sgf_ax_fs.o sgf_ax_g.o sgf_ax_str.o splc_clm.o strint.o strint2.o thomas.o \
geten2.o getef.o getef2.o
#
FFLAGS = -O3 -g
#
e_vdrop: $(OBJ)
	$(FC) $(FFLAGS) -o e_vdrop $(OBJ)	
edrop.o: edrop.for
	$(FC) $(FFLAGS) -c edrop.for
arc_evl.o: arc_evl.for
	$(FC) $(FFLAGS) -c arc_evl.for
bess_I01K01.o: bess_I01K01.f
	$(FC) $(FFLAGS) -c bess_I01K01.f
df.o: df.for
	$(FC) $(FFLAGS) -c df.for
dgesv.o: dgesv.f
	$(FC) $(FFLAGS) -c dgesv.f
ell_int.o: ell_int.f
	$(FC) $(FFLAGS) -c ell_int.f
ell_int2.o: ell_int2.for
	$(FC) $(FFLAGS) -c ell_int2.for
gauss_leg.o: gauss_leg.f
	$(FC) $(FFLAGS) -c gauss_leg.f
gauss_log.o: gauss_log.for
	$(FC) $(FFLAGS) -c gauss_log.for
geten2.o: geten2.for
	$(FC) $(FFLAGS) -c geten2.for
getef.o: getef.for
	$(FC) $(FFLAGS) -c getef.for
getef2.o: getef2.for
	$(FC) $(FFLAGS) -c getef2.for
gsint.o: gsint.for
	$(FC) $(FFLAGS) -c gsint.for
lgf_ax_fs.o: lgf_ax_fs.f
	$(FC) $(FFLAGS) -c lgf_ax_fs.f
polyn.o: polyn.for
	$(FC) $(FFLAGS) -c polyn.for
regridf.o: regridf.for
	$(FC) $(FFLAGS) -c regridf.for
regridfadp.o: regridfadp.for
	$(FC) $(FFLAGS) -c regridfadp.for
#regridf3.o: regridf3.for
#	$(FC) $(FFLAGS) -c regridf3.for
seval.o: seval.f
	$(FC) $(FFLAGS) -c seval.f
sgf_ax_fs.o: sgf_ax_fs.f
	$(FC) $(FFLAGS) -c sgf_ax_fs.f
sgf_ax_g.o: sgf_ax_g.for
	$(FC) $(FFLAGS) -c sgf_ax_g.for
sgf_ax_str.o: sgf_ax_str.for
	$(FC) $(FFLAGS) -c sgf_ax_str.for
splc_clm.o: splc_clm.f
	$(FC) $(FFLAGS) -c splc_clm.f
strint.o: strint.for
	$(FC) $(FFLAGS) -c strint.for
strint2.o: strint2.for
	$(FC) $(FFLAGS) -c strint2.for
thomas.o: thomas.f
	$(FC) $(FFLAGS) -c thomas.f
#
#clean:
#\rm -rf *.o pinch
#cleaner: clean
#\rm -rf output/*
#
