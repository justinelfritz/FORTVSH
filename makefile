NAME   = vsh_test.exe
STALIB = libvsh.a
SOLIB  = libvsh.so

BASE   = $(HOME)/Desktop/magfric/fortran_VSH
EXEC   = $(BASE)/$(NAME)
LIBDIR = $(BASE)/lib
MODDIR = $(BASE)/mod
OBJDIR = $(BASE)/obj
SRCDIR = $(BASE)/src

FC      = gfortran
FCFLAGS = -O3 -fPIC

LIB_OBJS = $(OBJDIR)/kinds.o $(OBJDIR)/globals.o $(OBJDIR)/vsh.o
TST_OBJS = $(OBJDIR)/tests.o $(OBJDIR)/main.o

.PHONY: all static shared clean

all: static $(EXEC)

# ── Static library ────────────────────────────────────────────────────────────
static: $(LIBDIR)/$(STALIB)

$(LIBDIR)/$(STALIB): $(LIB_OBJS) | $(LIBDIR)
	ar rcs $@ $^

# ── Shared library (optional) ─────────────────────────────────────────────────
shared: $(LIBDIR)/$(SOLIB)

$(LIBDIR)/$(SOLIB): $(LIB_OBJS) | $(LIBDIR)
	$(FC) $(FCFLAGS) -shared -o $@ $^

# ── Library objects (module output → mod/) ────────────────────────────────────
$(OBJDIR)/kinds.o: $(SRCDIR)/kinds.f90 | $(OBJDIR) $(MODDIR)
	@echo Building kinds.o
	$(FC) $(FCFLAGS) -J $(MODDIR) -c $< -o $@

$(OBJDIR)/globals.o: $(SRCDIR)/globals.f90 $(OBJDIR)/kinds.o | $(OBJDIR) $(MODDIR)
	@echo Building globals.o
	$(FC) $(FCFLAGS) -I $(MODDIR) -J $(MODDIR) -c $< -o $@

$(OBJDIR)/vsh.o: $(SRCDIR)/vsh.f90 $(OBJDIR)/globals.o | $(OBJDIR) $(MODDIR)
	@echo Building vsh.o
	$(FC) $(FCFLAGS) -I $(MODDIR) -J $(MODDIR) -c $< -o $@

# ── Test objects (read library mods from mod/, write tests.mod into obj/) ─────
$(OBJDIR)/tests.o: $(SRCDIR)/tests.f90 $(LIBDIR)/$(STALIB) | $(OBJDIR)
	@echo Building tests.o
	$(FC) $(FCFLAGS) -I $(MODDIR) -J $(OBJDIR) -c $< -o $@

$(OBJDIR)/main.o: $(SRCDIR)/main.f90 $(OBJDIR)/tests.o | $(OBJDIR)
	@echo Building main.o
	$(FC) $(FCFLAGS) -I $(MODDIR) -I $(OBJDIR) -c $< -o $@

# ── Test executable ───────────────────────────────────────────────────────────
$(EXEC): $(TST_OBJS) $(LIBDIR)/$(STALIB)
	@echo Building $(NAME)
	$(FC) $(FCFLAGS) -o $@ $(TST_OBJS) -L $(LIBDIR) -lvsh

# ── Directory bootstrap ───────────────────────────────────────────────────────
$(LIBDIR) $(MODDIR) $(OBJDIR):
	mkdir -p $@

# ── Clean ─────────────────────────────────────────────────────────────────────
clean:
	rm -f $(OBJDIR)/*.o $(OBJDIR)/*.mod \
	      $(MODDIR)/*.mod \
	      $(LIBDIR)/*.a  $(LIBDIR)/*.so \
	      validation/*.dat $(EXEC)
