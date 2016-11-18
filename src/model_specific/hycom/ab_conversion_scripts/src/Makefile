#
.SUFFIXES:
.SUFFIXES: .c .F .f .o

.F:
	@echo "Must have an explicit rule for" $*
	@exit 1
.f:
	@echo "Must have an explicit rule for" $*
	@exit 1

include intelIFC_setup

MODS = mod_xc.o mod_za.o wtime.o
read_archive:  $(MODS) read_archive.o
	$(LD) $(LDFLAGS) read_archive.o  $(MODS) -o read_archive
write_archive:  $(MODS) write_archive.o
	$(LD) $(LDFLAGS) write_archive.o  $(MODS) -o write_archive
mod_xc.o:    mod_xc.F
mod_za.o:    mod_za.F mod_xc.o
