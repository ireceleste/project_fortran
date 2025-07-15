FC = gfortran
FCFLAGS = -J./mod -I./mod 
SRCDIR = ./src
MODDIR = ./mod
OBJDIR = ./obj
BINDIR = ./bin
OUTDIR = ./output

MODULE_SRCS = \
    $(SRCDIR)/modules/utils.f08 \
	$(SRCDIR)/modules/HistogramHandler.f08 \
	$(SRCDIR)/modules/FitterModule.f08 \
    $(SRCDIR)/modules/Minimizer.f08 \
    $(SRCDIR)/modules/FunctionModule.f08 \
    $(SRCDIR)/modules/GaussianGenerator.f08 \
	$(SRCDIR)/modules/InputOutput.f08 \

OBJS = $(patsubst $(SRCDIR)/modules/%.f08, $(OBJDIR)/%.o, $(MODULE_SRCS))


all: directories main

# The main target compiles the main program and links it with the object files

main: $(OBJS) $(SRCDIR)/main.f08
	$(FC) $(FCFLAGS) -llapack -lblas -o $(BINDIR)/$@ $^	

$(OBJDIR)/%.o: $(SRCDIR)/modules/%.f08
	$(FC) $(FCFLAGS) -c $< -o $@


# The directories target ensures that the necessary directories exist before building

directories: $(BINDIR) $(MODDIR) $(OBJDIR) $(OUTDIR)

$(BINDIR):
	mkdir -p $(BINDIR)

$(MODDIR):
	mkdir -p $(MODDIR)

$(OBJDIR):
	mkdir -p $(OBJDIR)

${OUTDIR}:
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