# Makefile: compile every SRC/*.F90 into BIN/*.so via elmerf90
#
# Usage:
#   make            # build all shared objects into BIN/
#   make clean      # remove BIN/ and its contents
#   make list       # show what would be built
#
# Requires a loaded Elmer environment (elmerf90 / ElmerSolver_mpi on PATH).

SRCDIR  := SRC
BINDIR  := BIN

# elmerf90 wrapper already supplies -shared -fPIC and the Elmer include/link flags
ELMERF90 ?= elmerf90

# All Fortran sources -> matching shared objects in BINDIR
SOURCES := $(wildcard $(SRCDIR)/*.F90)
OBJECTS := $(patsubst $(SRCDIR)/%.F90,$(BINDIR)/%.so,$(SOURCES))

.PHONY: all clean list check

all: check $(OBJECTS)

# Fail early with a clear message if the Elmer environment is not loaded
check:
	@command -v ElmerSolver_mpi >/dev/null 2>&1 || \
		{ echo "ERROR: Elmer environment not loaded (ElmerSolver_mpi not found on PATH)."; exit 1; }
	@command -v $(ELMERF90) >/dev/null 2>&1 || \
		{ echo "ERROR: $(ELMERF90) not found on PATH."; exit 1; }

# filename.F90 -> filename.so, placed in BINDIR
$(BINDIR)/%.so: $(SRCDIR)/%.F90 | $(BINDIR)
	$(ELMERF90) -o $@ $<

# Create the output directory on demand (order-only prerequisite)
$(BINDIR):
	mkdir -p $(BINDIR)

list:
	@echo "Sources: $(SOURCES)"
	@echo "Targets: $(OBJECTS)"

clean:
	rm -rf $(BINDIR)
