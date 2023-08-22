#----------------------------------------------------------------------------------------------
# makefile:  Controls the build of the primary ParaDiS executable
#            and the parallel FMM correction table generator (ctablegenp)
#
# Usage:
#    make              : builds the paradis executable and FMM correction table generator
#    make clean        : removes executables and object files
#    make depend       : creates the makefile.dep dependency file
#    make purify       : 
#    make prof         :
#    make predefines   : simple tool for listing all the defined preprocessor macros
#                        that are enabled/set under the current build settings.
#----------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------
# MAKEFLAGS = -j<num>   where <num> is the number of simultanueous targets to process
#----------------------------------------------------------------------------------------------

MAKEFLAGS  = -j1         # (use for single-threaded builds)
MAKEFLAGS  = -j16

MAKEFILES = makefile ../makefile.sys ../makefile.setup

#-----------------------------------------------------------------------------------------------
# makefile.sys   : contains system-specific macro definitions as well as user-settings 
# makefile.setup : contains user-specific settings and option flags 
#-----------------------------------------------------------------------------------------------

INCDIR = ../include
BINDIR = ../bin

include ../makefile.setup

#-----------------------------------------------------------------------------------------------
# The destination object directory depends on the compilation mode set in makefile.setup.
#    ../obj/p = PARALLEL objects
#    ../obj/s = SERIAL   objects
#-----------------------------------------------------------------------------------------------

ifeq ($(MODE),PARALLEL)
OBJDIR = ../obj/p
else
OBJDIR = ../obj/s
endif

LIB_PARADIS = $(OBJDIR)/libparadis.a 

include ../makefile.srcs

#-----------------------------------------------------------------------------------------------
# compilation suffix rules....
#-----------------------------------------------------------------------------------------------

.SUFFIXES : .cc .cpp 

$(OBJDIR)/%.o: ../src/%.cc   $(MAKEFILES)
	@echo $(CPP) $<
	@$(CPP) $(OPT) $(CPPFLAG) $(INCS) -c $< -o $@

$(OBJDIR)/%.o: ../src/%.cpp  $(MAKEFILES)
	@echo $(CPP) $<
	@$(CPP) $(OPT) $(CPPFLAG) $(INCS) -c $< -o $@

#----------------------------------------------------------------------------------------------
# Define all targets and dependencies below
#----------------------------------------------------------------------------------------------

.PHONY: all
all: info $(OBJDIR) $(BINDIR) paradis paradis_net ctablegenp 

.PHONY: clean
clean:
	rm -rf $(OBJDIR) .DS_Store *.core *.rslts *_results *.log
	rm -f  $(BINDIR)/ctablegenp
	rm -f  $(BINDIR)/paradis
	rm -f  $(BINDIR)/paradis_net

.PHONY: info
info:
	@echo " --------------------------------------------------------------------------------------------"
	@echo " operating system = " $(OPSYS)
	@echo " node             = " $(NODE)
	@echo " compile target   = " $(SYS)
	@echo " compile mode     = " $(MODE)
	@echo " cpp compiler     = " $(CPP)
	@echo " cpp options      = " $(OPT) $(CPPFLAG)
	@echo " cpp includes     = " $(INCS)
ifeq ($(GPU_ENABLED),ON)
	@echo " nvcc compiler    = " $(NVCC)
	@echo " nvcc options     = " $(NVCC_FLAGS) $(DEFS)
	@echo " nvcc includes    = " $(INCS) $(NVCC_INCLUDES)
endif
	@echo " lib options      = " $(LIB_PARADIS) $(LIB)
	@echo " --------------------------------------------------------------------------------------------"

#----------------------------------------------------------------------------------------------

$(OBJDIR) :
	mkdir -p $(OBJDIR)

$(BINDIR) :
	mkdir -p $(BINDIR) 

#----------------------------------------------------------------------------------------------
# The depend target will execute makedepend with the files in the current directory.
# Note that we need use sed to update the .o and .h locations to reflect the current
# INCDIR and OBJDIR locations.
#----------------------------------------------------------------------------------------------

.PHONY: depend
depend: *.cc $(INCDIR)/*.h $(MAKEFILES)
	makedepend -Y$(INCDIR) $(DEFS) *.cc -fmakefile.dep -w 200
	sed -i -e 's/\.\.\/include/\$$\(INCDIR\)/g' makefile.dep;
	sed -i -e 's/^/\$$(OBJDIR)\//g'             makefile.dep;
	sed -i -e '1,2d'                            makefile.dep;
	printf "# DO NOT DELETE\n\n" | cat -        makefile.dep > makefile.dep.tmp; mv makefile.dep.tmp makefile.dep;
	rm -f makefile.dep.bak
	rm -f makefile.dep-e

#----------------------------------------------------------------------------------------------
# main targets...
#----------------------------------------------------------------------------------------------

paradis               : $(BINDIR) $(OBJDIR) $(LIB_PARADIS) $(BINDIR)/paradis 
$(BINDIR)/paradis     : $(OBJDIR)/ParaDiS.o     $(LIB_PARADIS) 
	@echo creating application $@
	@$(CPP) $(OPT) $(OPENMP_FLAG) -o $@ $(OBJDIR)/ParaDiS.o    $(LIB_PARADIS) $(LIB)

paradis_net           : $(BINDIR) $(OBJDIR) $(LIB_PARADIS) $(BINDIR)/paradis_net
$(BINDIR)/paradis_net : $(OBJDIR)/ParaDiS_Net.o $(LIB_PARADIS) 
	@echo creating application $@
	@$(CPP) $(OPT) $(OPENMP_FLAG) -o $@ $(OBJDIR)/ParaDiS_Net.o $(LIB_PARADIS) $(LIB)

ctablegenp            : $(BINDIR) $(OBJDIR) $(LIB_PARADIS) $(BINDIR)/ctablegenp
$(BINDIR)/ctablegenp  : $(OBJDIR)/CTableGen.o   $(LIB_PARADIS) 
	@echo creating application $@
	@$(CPP) $(OPT) $(OPENMP_FLAG) -o $@ $(OBJDIR)/CTableGen.o   $(LIB_PARADIS) $(LIB)

#----------------------------------------------------------------------------------------------
# There is currently a bug in the compiler on LLNL/BGL which causes errors when the ghost comm 
# module is compiled with optimizations. This temporarily adds a target for CommSendGhosts.o 
# which forces CommSendGhosts.cc do be compiled with the optimizations disabled.
#----------------------------------------------------------------------------------------------

$(OBJDIR)/CommSendGhosts.o: ../src/CommSendGhosts.cc $(MAKEFILES)
	@echo $(CPP) $< "(n.b. optimization disabled)"
	@$(CPP) -g  $(CPPFLAG) $(INCS) -c $< -o $@

#----------------------------------------------------------------------------------------------
# If we are compiling with GPU-enabled, all the modules with GPU device code needs to be
# compiled using NVCC...
#----------------------------------------------------------------------------------------------

ifeq ($(GPU_ENABLED),ON)

$(OBJDIR)/GPU_Error.o                        : ../src/GPU_Error.cc                         $(MAKEFILES)
	@echo $(NVCC) $<  
	@$(NVCC) -x cu $(NVCC_FLAGS) $(DEFS) $(INCS) $(NVCC_INCLUDES) -c $< -o $@
$(OBJDIR)/GPU_Node.o                         : ../src/GPU_Node.cc                          $(MAKEFILES)
	@echo $(NVCC) $<  
	@$(NVCC) -x cu $(NVCC_FLAGS) $(DEFS) $(INCS) $(NVCC_INCLUDES) -c $< -o $@
$(OBJDIR)/GPU_Segment.o                      : ../src/GPU_Segment.cc                       $(MAKEFILES)
	@echo $(NVCC) $<  
	@$(NVCC) -x cu $(NVCC_FLAGS) $(DEFS) $(INCS) $(NVCC_INCLUDES) -c $< -o $@
$(OBJDIR)/PBC.o                              : ../src/PBC.cc                               $(MAKEFILES)
	@echo $(NVCC) $<  
	@$(NVCC) -x cu $(NVCC_FLAGS) $(DEFS) $(INCS) $(NVCC_INCLUDES) -c $< -o $@
$(OBJDIR)/SSF_Aniso.o                        : ../src/SSF_Aniso.cc                         $(MAKEFILES)
	@echo $(NVCC) $<  
	@$(NVCC) -x cu $(NVCC_FLAGS) $(DEFS) $(INCS) $(NVCC_INCLUDES) -c $< -o $@
$(OBJDIR)/SSF_Iso.o                          : ../src/SSF_Iso.cc                           $(MAKEFILES)
	@echo $(NVCC) $<  
	@$(NVCC) -x cu $(NVCC_FLAGS) $(DEFS) $(INCS) $(NVCC_INCLUDES) -c $< -o $@
$(OBJDIR)/SSF_Iso_Rational.o                 : ../src/SSF_Iso_Rational.cc                  $(MAKEFILES)
	@echo $(NVCC) $<  
	@$(NVCC) -x cu $(NVCC_FLAGS) $(DEFS) $(INCS) $(NVCC_INCLUDES) -c $< -o $@
$(OBJDIR)/SSF_Driver.o                       : ../src/SSF_Driver.cc                        $(MAKEFILES)
	@echo $(NVCC) $<  
	@$(NVCC) -x cu $(NVCC_FLAGS) $(DEFS) $(INCS) $(NVCC_INCLUDES) -c $< -o $@
$(OBJDIR)/SSF_Aniso_NonParallel.o            : ../src/SSF_Aniso_NonParallel.cc             $(MAKEFILES)
	@echo $(NVCC) $<  
	@$(NVCC) -x cu $(NVCC_FLAGS) $(DEFS) $(INCS) $(NVCC_INCLUDES) -c $< -o $@
$(OBJDIR)/SSF_Aniso_Parallel.o               : ../src/SSF_Aniso_Parallel.cc                $(MAKEFILES)
	@echo $(NVCC) $<  
	@$(NVCC) -x cu $(NVCC_FLAGS) $(DEFS) $(INCS) $(NVCC_INCLUDES) -c $< -o $@
$(OBJDIR)/SSF_Node_t.o                       : ../src/SSF_Node_t.cc                        $(MAKEFILES)
	@echo $(NVCC) $<  
	@$(NVCC) -x cu $(NVCC_FLAGS) $(DEFS) $(INCS) $(NVCC_INCLUDES) -c $< -o $@
$(OBJDIR)/SSF_FV_t.o                         : ../src/SSF_FV_t.cc                          $(MAKEFILES)
	@echo $(NVCC) $<  
	@$(NVCC) -x cu $(NVCC_FLAGS) $(DEFS) $(INCS) $(NVCC_INCLUDES) -c $< -o $@
$(OBJDIR)/SSF_PV_t.o                         : ../src/SSF_PV_t.cc                          $(MAKEFILES)
	@echo $(NVCC) $<  
	@$(NVCC) -x cu $(NVCC_FLAGS) $(DEFS) $(INCS) $(NVCC_INCLUDES) -c $< -o $@

endif

#----------------------------------------------------------------------------------------------

ifeq ($(OPSYS),Darwin)
   RANLIB_ARGS = -no_warning_for_no_symbols 
endif

$(LIB_PARADIS) : $(OBJDIR) $(PARADIS_OBJS)
	@echo building library $(LIB_PARADIS)
	@ar -rcS $(LIB_PARADIS) $(PARADIS_OBJS)
	@ranlib $(RANLIB_ARGS) $(LIB_PARADIS)

#----------------------------------------------------------------------------------------------
# deal with cuda relocatable device code...
#----------------------------------------------------------------------------------------------

$(OBJDIR)/nvcc_rdc_dlink.o : $(PARADIS_RDC_OBJS)
	@echo linking relocatable gpu device code
	@$(NVCC) -Wno-deprecated-gpu-targets -dlink -o $(OBJDIR)/nvcc_rdc_dlink.o $(PARADIS_RDC_OBJS)

#----------------------------------------------------------------------------------------------
# diagnostic targets...
#----------------------------------------------------------------------------------------------

.PHONY: purify
purify: $(OBJDIR)/ParaDiS.o $(PARADIS_OBJS) 
	purify ./gcc $(OBJDIR)/ParaDiS.o $(OPT) $(PARADIS_OBJS) -o paradis $(LIB) 

.PHONY: prof
prof: $(OBJDIR)/ParaDiS.o $(PARADIS_OBJS) 
	$(CPP) $(OBJDIR)/ParaDiS.o $(OPT) $(PARADIS_OBJS) -p  -o paradis $(LIB) 

#----------------------------------------------------------------------------------------------
# The predefines target is useful for determining the state of the currently defined 
# preprocessor macros within the compilation environment.  The target essentially invokes the 
# compiler with the current options and prints a human readable list of the resulting C/C++ 
# preprocessor macros that are defined. No code is actually compiled.
#----------------------------------------------------------------------------------------------

.PHONY: predefines
predefines:
	$(CPP) $(OPT) $(CPPFLAG) $(INCS) -dM -E -c ParaDiS.cc | cut -c9- | sort
	rm -f ParaDiS.o

#----------------------------------------------------------------------------------------------
# The includes target was added to obtain a list of the paths of the various include
# files that are pulled in for building the ParaDiS application.
#----------------------------------------------------------------------------------------------

.PHONY: includes
includes:
	$(CPP) $(OPT) $(CPPFLAG) $(INCS) -H -c ParaDiS.cc | cut -c9- | sort --unique --ignore-nonprinting
	rm -f ParaDiS.o

# The warnings target is useful for displaying the status of the current GCC compiler warnings.
#----------------------------------------------------------------------------------------------

.PHONY: warnings
warnings:
	gcc $(OPT) $(CPPFLAG) -Q --help=warning - < /dev/null 

#----------------------------------------------------------------------------------------------
# finally - include the actual dependencies...
#----------------------------------------------------------------------------------------------

include ./makefile.dep

