# Customization; see Options.mk.example
CONFIG ?= ../Options.mk

include $(CONFIG)

INCL=../libgadget/config.h

include ../Makefile.rules

OBJS = main.o params.o

OBJS2 = main_cicass.o params.o

OBJS_gadget = ../libgadget/.objs/gdbtools.o ../libgadget/.objs/hci.o ../libgadget/.objs/fof.o ../libgadget/.objs/fofpetaio.o ../libgadget/.objs/petaio.o ../libgadget/.objs/domain.o ../libgadget/.objs/exchange.o ../libgadget/.objs/slotsmanager.o ../libgadget/.objs/partmanager.o ../libgadget/.objs/blackhole.o ../libgadget/.objs/timebinmgr.o ../libgadget/.objs/run.o ../libgadget/.objs/drift.o ../libgadget/.objs/stats.o ../libgadget/.objs/timestep.o ../libgadget/.objs/init.o ../libgadget/.objs/checkpoint.o ../libgadget/.objs/sfr_eff.o ../libgadget/.objs/cooling.o ../libgadget/.objs/cooling_rates.o ../libgadget/.objs/cooling_uvfluc.o ../libgadget/.objs/cooling_qso_lightup.o ../libgadget/.objs/winds.o ../libgadget/.objs/density.o ../libgadget/.objs/metal_return.o ../libgadget/.objs/treewalk.o ../libgadget/.objs/cosmology.o ../libgadget/.objs/gravshort-tree.o ../libgadget/.objs/gravshort-pair.o ../libgadget/.objs/hydra.o ../libgadget/.objs/timefac.o ../libgadget/.objs/gravpm.o ../libgadget/.objs/powerspectrum.o ../libgadget/.objs/forcetree.o ../libgadget/.objs/petapm.o ../libgadget/.objs/gravity.o ../libgadget/.objs/densitykernel.o ../libgadget/.objs/lightcone.o ../libgadget/.objs/walltime.o ../libgadget/.objs/runtests.o ../libgadget/.objs/neutrinos_lra.o ../libgadget/.objs/omega_nu_single.o ../libgadget/.objs/config.o

all: MP-GenIC MP-GenIC-CICASS

MP-GenIC: $(OBJS) ../libgenic/libgenic.a $(OBJS_gadget) ../libgadget/libgadget-utils.a
	$(MPICC) $(OPTIMIZE) $^ $(LIBS) -o $@

MP-GenIC-CICASS: $(OBJS2) ../libgenic/libgenic.a ../libgadget/libgadget.a ../libgadget/libgadget-utils.a
	$(MPICC) $(OPTIMIZE) $^ $(LIBS) -o $@

clean:
	rm -rf $(OBJS) MP-GenIC
