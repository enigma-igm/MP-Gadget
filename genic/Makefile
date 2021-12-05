# Customization; see Options.mk.example
CONFIG ?= ../Options.mk

include $(CONFIG)

INCL=../libgadget/config.h

include ../Makefile.rules

OBJS = main.o params.o

OBJS2 = main_cicass.o params.o

OBJS := $(OBJS:%.o=.objs/%.o)

all: MP-GenIC MP-GenIC-CICASS

MP-GenIC: $(OBJS) ../libgenic/libgenic.a ../libgadget/libgadget.a ../libgadget/libgadget-utils.a
	$(MPICC) $(OPTIMIZE) $^ $(LIBS) -o $@

MP-GenIC-CICASS: $(OBJS2) ../libgenic/libgenic.a ../libgadget/libgadget.a ../libgadget/libgadget-utils.a
	$(MPICC) $(OPTIMIZE) $^ $(LIBS) -o $@

clean:
	rm -rf $(OBJS) MP-GenIC
