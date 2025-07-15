FC = gfortran
FCFLAGS = -J./mod -I./mod 
LDFLAGS = -llapack -lblas

SRCDIR = ./src
MODDIR = ./mod
OBJDIR = ./obj
BINDIR = ./bin
OUTDIR = ./output

MODULE_SRCS = \
    $(SRCDIR)/modules/utils.f90 \
	$(SRCDIR)/modules/HistogramHandler.f90 \
	$(SRCDIR)/modules/FitterModule.f90 \
    $(SRCDIR)/modules/Minimizer.f90 \
    $(SRCDIR)/modules/FunctionModule.f90 \
    $(SRCDIR)/modules/GaussianGenerator.f90 \
	$(SRCDIR)/modules/InputOutput.f90 \

OBJS = $(patsubst $(SRCDIR)/modules/%.f90, $(OBJDIR)/%.o, $(MODULE_SRCS))


all: directories main

# The main target compiles the main program and links it with the object files

main: $(OBJS) $(SRCDIR)/main.f90
	$(FC) $(FCFLAGS) $(LDFLAGS) -o $(BINDIR)/$@ $^	

$(OBJDIR)/%.o: $(SRCDIR)/modules/%.f90
	$(FC) $(FCFLAGS) -c $< -o $@


# The directories target ensures that the necessary directories exist before building

directories: $(BINDIR) $(MODDIR) $(OBJDIR) $(OUTDIR)

$(BINDIR):
	mkdir -p $(BINDIR)

$(MODDIR):
	mkdir -p $(MODDIR)

$(OBJDIR):
	mkdir -p $(OBJDIR)

$(OUTDIR):
	mkdir -p $(OUTDIR)

# The clean target removes all compiled files and directories

.PHONY: clean
clean: 
	rm -rf $(BINDIR)/* $(MODDIR)/*.mod $(OBJDIR)/*.o 
	rmdir $(BINDIR) $(MODDIR) $(OBJDIR) 

# The clean_output removes all output files and plots

.PHONY: clean_output
clean_output:
	rm -f $(OUTDIR)/*
	rmdir $(OUTDIR)