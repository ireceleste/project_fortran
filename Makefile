FC = gfortran-14
FCFLAGS = -J./mod -I./mod -llapack -lblas
SRCDIR = ./src
MODDIR = ./mod
OBJDIR = ./obj
BINDIR = ./bin

# VPATH specifies the directories to be searched for modules
VPATH = $(SRCDIR)/modules
# The wildcard function is used to generate a list of .f08 files in the modules directory
# Collect all modules
# List source files in correct compilation order
MODULE_SRCS = \
    $(SRCDIR)/modules/utils.f08 \
	$(SRCDIR)/modules/HistogramHandler.f08 \
	$(SRCDIR)/modules/FitterModule.f08 \
    $(SRCDIR)/modules/Minimizer.f08 \
    $(SRCDIR)/modules/FunctionModule.f08 \
    $(SRCDIR)/modules/GaussianGenerator.f08 \
	$(SRCDIR)/modules/InputOutput.f08 \


# Object files in the desired order: utils.o, Minimizer.o, others...
OBJS = $(patsubst $(SRCDIR)/modules/%.f08, $(OBJDIR)/%.o, $(MODULE_SRCS))


all: directories main

# The main target compiles the main program and links it with the object files

main: $(OBJS) $(SRCDIR)/main.f08
	$(FC) $(FCFLAGS) -o $(BINDIR)/$@ $^	

$(OBJDIR)/%.o: $(SRCDIR)/modules/%.f08
	$(FC) $(FCFLAGS) -c $< -o $@

# The directories target ensures that the necessary directories exist before building

directories: $(BINDIR) $(MODDIR) $(OBJDIR)

$(BINDIR):
	mkdir -p $(BINDIR)

$(MODDIR):
	mkdir -p $(MODDIR)

$(OBJDIR):
	mkdir -p $(OBJDIR)

# The clean target removes all compiled files and directories
.PHONY: clean
clean: 
	rm -f $(BINDIR)/* $(MODDIR)/*.mod $(OBJDIR)/*.o