#-----------------------------------------------------------------------------------------------
#    makefile: controls the build of various ParaDiS support utilities
#
#    Builds the following utilities:
#
#        paradisgen     --  problem generator
#        paradisrepart  --  creates a new problem decomposition with a new
#                           domain geometry from an previous nodal data file
#        paradisconvert --  Converts older format restart files to the current format
#        ctablegen      --  Serial version of utility used to create the FMM correction tables.
#        calcdensity    --  Calculates a dislocation density 'grid' from
#                           a specified restart file and writes it to a
#                           file for visualization via an external tool
#        inclusiongen   --  Generates eshelby inclusion data in a format
#                           suitable for use in ParaDiS simulations.
#
#
#	NOTE: The utilities use various source modules from the parallel
#             code located in the parent directory.  These modules are
#             compiled in the parent directory with a different set of
#             preprocessor definitions than are needed here, so we create
#             links in this directory back to the source modlues and create
#             separate object modules for these sources.
#-----------------------------------------------------------------------------------------------

MAKEFLAGS  = -j1         # (use for single-threaded builds)
MAKEFLAGS  = -j16

MAKEFILES  = makefile ../makefile.sys ../makefile.setup

INCDIR = ../include
OBJDIR = ../obj/s
BINDIR = ../bin

LIB_PARADIS = $(OBJDIR)/libparadis.a 

include ../makefile.setup
include ../makefile.srcs

#-----------------------------------------------------------------------------------------------
# compilation suffix rules....
#-----------------------------------------------------------------------------------------------

.SUFFIXES : .cc .cpp 

$(OBJDIR)/%.o: %.cc          $(MAKEFILES)
	@echo $(CPP_SERIAL) $<
	@$(CPP_SERIAL) $(OPT) $(CPPFLAG_SERIAL) $(INCS_SERIAL) -c $< -o $@

$(OBJDIR)/%.o: %.cpp         $(MAKEFILES)
	@echo $(CPP_SERIAL) $<
	@$(CPP_SERIAL) $(OPT) $(CPPFLAG_SERIAL) $(INCS_SERIAL) -c $< -o $@

$(OBJDIR)/%.o: ../src/%.cc   $(MAKEFILES)
	@echo $(CPP_SERIAL) $<
	@$(CPP_SERIAL) $(OPT) $(CPPFLAG_SERIAL) $(INCS_SERIAL) -c $< -o $@

$(OBJDIR)/%.o: ../src/%.cpp  $(MAKEFILES)
	@echo $(CPP_SERIAL) $<
	@$(CPP_SERIAL) $(OPT) $(CPPFLAG_SERIAL) $(INCS_SERIAL) -c $< -o $@

#-----------------------------------------------------------------------------------------------

UTIL_APPS  =
UTIL_APPS += $(BINDIR)/paradisgen
UTIL_APPS += $(BINDIR)/paradisrepart
UTIL_APPS += $(BINDIR)/paradisconvert
UTIL_APPS += $(BINDIR)/calcdensity
UTIL_APPS += $(BINDIR)/ctablegen
UTIL_APPS += $(BINDIR)/inclusiongen
UTIL_APPS += $(BINDIR)/precipitates
UTIL_APPS += $(BINDIR)/paradis_extract
UTIL_APPS += $(BINDIR)/paradis_gnu
UTIL_APPS += $(BINDIR)/paradis_compress
UTIL_APPS += $(BINDIR)/paradis_refine
UTIL_APPS += $(BINDIR)/paradis_stat
UTIL_APPS += $(BINDIR)/paradis_bvdump
UTIL_APPS += $(BINDIR)/stresstablegen
UTIL_APPS += $(BINDIR)/gen_fmm_derivatives
UTIL_APPS += $(BINDIR)/StanfordFMMsTableGen

ifeq ($(GPU_ENABLED),ON)
UTIL_APPS += $(BINDIR)/gpu_info
endif

.PHONY: all
all:  info $(OBJDIR) $(BINDIR) $(UTIL_APPS)

.PHONY: clean
clean:
	rm -rf $(UTIL_APPS) $(OBJDIR)/*.o $(LIB_PARADIS) .DS_Store *.core

$(OBJDIR) :
	mkdir -p $(OBJDIR)

$(BINDIR) :
	mkdir -p $(BINDIR)

#-----------------------------------------------------------------------------------------------

$(BINDIR)/paradisgen           : $(OBJDIR)/ParadisGen.o $(OBJDIR)/CreateConfig.o $(OBJDIR)/InitRemesh.o $(LIB_PARADIS) $(MAKEFILES)
	@echo creating application $@
	@$(CPP_SERIAL) $(OPT) $(CPPFLAG_SERIAL) $(INCS_SERIAL) -o $@ $(OBJDIR)/ParadisGen.o $(OBJDIR)/CreateConfig.o $(OBJDIR)/InitRemesh.o $(LIB_PARADIS) $(LIB_SERIAL)

$(BINDIR)/paradisrepart        : $(OBJDIR)/ParadisRepart.o $(LIB_PARADIS) $(MAKEFILES)
	@echo creating application $@
	@$(CPP_SERIAL) $(OPT) $(CPPFLAG_SERIAL) $(INCS_SERIAL) -o $@ $(OBJDIR)/ParadisRepart.o $(LIB_PARADIS) $(LIB_SERIAL)

$(BINDIR)/paradisconvert       : $(OBJDIR)/ParadisConvert.o $(LIB_PARADIS) $(MAKEFILES)
	@echo creating application $@
	@$(CPP_SERIAL) $(OPT) $(CPPFLAG_SERIAL) $(INCS_SERIAL) -o $@ $(OBJDIR)/ParadisConvert.o $(LIB_PARADIS) $(LIB_SERIAL)

$(BINDIR)/calcdensity          : $(OBJDIR)/CalcDensity.o $(LIB_PARADIS) $(MAKEFILES)
	@echo creating application $@
	@$(CPP_SERIAL) $(OPT) $(CPPFLAG_SERIAL) $(INCS_SERIAL) -o $@ $(OBJDIR)/CalcDensity.o $(LIB_PARADIS) $(LIB_SERIAL)

$(BINDIR)/ctablegen            : $(OBJDIR)/CTableGen.o $(LIB_PARADIS) $(MAKEFILES)
	@echo creating application $@
	@$(CPP_SERIAL) $(OPT) $(CPPFLAG_SERIAL) $(INCS_SERIAL) -o $@ $(OBJDIR)/CTableGen.o $(LIB_PARADIS) $(LIB_SERIAL)

$(BINDIR)/inclusiongen         : $(OBJDIR)/InclusionGen.o $(LIB_PARADIS) $(MAKEFILES)
	@echo creating application $@
	@$(CPP_SERIAL) $(OPT) $(CPPFLAG_SERIAL) $(INCS_SERIAL) -o $@ $(OBJDIR)/InclusionGen.o $(LIB_PARADIS) $(LIB_SERIAL)

$(BINDIR)/stresstablegen       : $(OBJDIR)/StressTableGen.o $(LIB_PARADIS) $(MAKEFILES)
	@echo creating application $@
	@$(CPP_SERIAL) $(OPT) $(CPPFLAG_SERIAL) $(INCS_SERIAL) -o $@ $(OBJDIR)/StressTableGen.o $(LIB_PARADIS) $(LIB_SERIAL)

$(BINDIR)/gen_fmm_derivatives  : $(OBJDIR)/GenAnisoFMMDerivTbl.o $(LIB_PARADIS) $(MAKEFILES)
	@echo creating application $@
	@$(CPP_SERIAL) $(OPT) $(CPPFLAG_SERIAL) $(INCS_SERIAL) -o $@ $(OBJDIR)/GenAnisoFMMDerivTbl.o $(LIB_PARADIS) $(LIB_SERIAL)

$(BINDIR)/precipitates         : $(OBJDIR)/Precipitates.o $(LIB_PARADIS) $(MAKEFILES)
	@echo creating application $@
	@$(CPP_SERIAL) $(OPT) $(CPPFLAG_SERIAL) $(INCS_SERIAL) -o $@ $(OBJDIR)/Precipitates.o $(LIB_PARADIS) $(LIB_SERIAL)

$(BINDIR)/paradis_extract      : $(OBJDIR)/paradis_extract.o  $(LIB_PARADIS) $(MAKEFILES)
	@echo creating application $@
	@$(CPP_SERIAL) $(OPT) $(CPPFLAG_SERIAL) $(INCS_SERIAL) -o $@ $(OBJDIR)/paradis_extract.o $(LIB_PARADIS) $(LIB_SERIAL)

$(BINDIR)/paradis_gnu          : $(OBJDIR)/paradis_gnu.o  $(LIB_PARADIS) $(MAKEFILES)
	@echo creating application $@
	@$(CPP_SERIAL) $(OPT) $(CPPFLAG_SERIAL) $(INCS_SERIAL) -o $@ $(OBJDIR)/paradis_gnu.o $(LIB_PARADIS) $(LIB_SERIAL)

$(BINDIR)/paradis_compress     : $(OBJDIR)/paradis_compress.o  $(LIB_PARADIS) $(MAKEFILES)
	@echo creating application $@
	@$(CPP_SERIAL) $(OPT) $(CPPFLAG_SERIAL) $(INCS_SERIAL) -o $@ $(OBJDIR)/paradis_compress.o $(LIB_PARADIS) $(LIB_SERIAL)

$(BINDIR)/paradis_refine       : $(OBJDIR)/paradis_refine.o  $(LIB_PARADIS) $(MAKEFILES)
	@echo creating application $@
	@$(CPP_SERIAL) $(OPT) $(CPPFLAG_SERIAL) $(INCS_SERIAL) -o $@ $(OBJDIR)/paradis_refine.o $(LIB_PARADIS) $(LIB_SERIAL)

$(BINDIR)/paradis_stat         : $(OBJDIR)/paradis_stat.o  $(LIB_PARADIS) $(MAKEFILES)
	@echo creating application $@
	@$(CPP_SERIAL) $(OPT) $(CPPFLAG_SERIAL) $(INCS_SERIAL) -o $@ $(OBJDIR)/paradis_stat.o $(LIB_PARADIS) $(LIB_SERIAL)

$(BINDIR)/paradis_bvdump       : $(OBJDIR)/paradis_bvdump.o  $(LIB_PARADIS) $(MAKEFILES)
	@echo creating application $@
	@$(CPP_SERIAL) $(OPT) $(CPPFLAG_SERIAL) $(INCS_SERIAL) -o $@ $(OBJDIR)/paradis_bvdump.o $(LIB_PARADIS) $(LIB_SERIAL)

$(BINDIR)/gpu_info             : $(OBJDIR)/gpu_info.o $(LIB_PARADIS) $(MAKEFILES)
	@echo creating application $@
	@$(CPP_SERIAL) $(OPT) $(CPPFLAG_SERIAL) $(INCS_SERIAL) -o $@ $(OBJDIR)/gpu_info.o $(LIB_PARADIS) $(LIB_SERIAL)

$(BINDIR)/StanfordFMMsTableGen : $(OBJDIR)/StanfordFMMsTableGen.o $(LIB_PARADIS) $(MAKEFILES)
	@echo creating application $@
	@$(CPP_SERIAL) $(OPT) $(CPPFLAG_SERIAL) $(INCS_SERIAL) -o $@ $(OBJDIR)/StanfordFMMsTableGen.o $(LIB_PARADIS) $(LIB_SERIAL)

#----------------------------------------------------------------------------------------------
# There is currently a bug in the compiler on LLNL/BGL which causes errors when the ghost comm 
# module is compiled with optimizations. This temporarily adds a target for CommSendGhosts.o 
# which forces CommSendGhosts.cc do be compiled with the optimizations disabled.
#----------------------------------------------------------------------------------------------

$(OBJDIR)/CommSendGhosts.o: ../src/CommSendGhosts.cc $(MAKEFILES)
	@echo $(CPP_SERIAL) $< "(n.b. optimization disabled)"
	@$(CPP_SERIAL) -g  $(CPPFLAG_SERIAL) $(INCS_SERIAL) -c $< -o $@

#----------------------------------------------------------------------------------------------
# If we are compiling active cuda code, all the modules with GPU device code needs to be
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
# Object files containing relocatable device code need a separate RDC link using nvcc...
#----------------------------------------------------------------------------------------------

$(OBJDIR)/nvcc_rdc_dlink.o : $(PARADIS_RDC_OBJS)
	@echo linking relocatable gpu device code
	@$(NVCC) -Wno-deprecated-gpu-targets -dlink -o $(OBJDIR)/nvcc_rdc_dlink.o $(PARADIS_RDC_OBJS)

#-----------------------------------------------------------------------------------------------

ifeq ($(OPSYS),Darwin)
   RANLIB_ARGS = -no_warning_for_no_symbols 
endif

$(LIB_PARADIS) : $(OBJDIR) $(PARADIS_OBJS)
	@echo building library $(LIB_PARADIS)
	@ar -rcS $(LIB_PARADIS) $(PARADIS_OBJS)
	@ranlib $(RANLIB_ARGS) $(LIB_PARADIS)

#-----------------------------------------------------------------------------------------------

.PHONY: info
info:
	@echo " --------------------------------------------------------------------------------------------"
	@echo " operating system = " $(OPSYS)
	@echo " node             = " $(NODE)
	@echo " compile target   = " $(SYS)
	@echo " cpp compiler     = " $(CPP_SERIAL)
	@echo " cpp options      = " $(OPT) $(CPPFLAG_SERIAL)
	@echo " cpp includes     = " $(INCS_SERIAL)
ifeq ($(GPU_ENABLED),ON)
	@echo " nvcc compiler    = " $(NVCC)
	@echo " nvcc options     = " $(NVCC_FLAGS) $(DEFS)
	@echo " nvcc includes    = " $(INCS_SERIAL) $(NVCC_INCLUDES)
endif
	@echo " lib options      = " $(LIB_PARADIS) $(LIB_SERIAL)
	@echo " --------------------------------------------------------------------------------------------"

#-----------------------------------------------------------------------------------------------

include ../src/makefile.dep
