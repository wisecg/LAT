include $(MGDODIR)/buildTools/config.mk

# Give the list of applications, which must be the stems of cc files with 'main'.
APPS = skim_mjd_data wave_skim ds_livetime auto-thresh validate_skim

# Stuff needed by BasicMakefile
SHLIB =
ARCHIVE =
# TAMDIR ?= $(ROOTSYS)
SOURCESSCRATCH = $(wildcard *.cc)

# Add RooFit stuff
# INCLUDEFLAGS = $(shell root-config --libs) -lRooFit -lRooFitCore -lMinuit

# Add MJSW stuff

INCLUDEFLAGS += $(CLHEP_INCLUDE_FLAGS) -I$(MGDODIR)/Base -I$(MGDODIR)/Root -I$(MGDODIR)/Transforms

INCLUDEFLAGS += -I$(MGDODIR)/Majorana -I$(MGDODIR)/MJDB $(ROOT_INCLUDE_FLAGS)  -I$(MGDODIR)/Tabree -I$(TAMDIR)/inc -I$(TAMDIR)/include

INCLUDEFLAGS += -I$(GATDIR)/BaseClasses -I$(GATDIR)/MGTEventProcessing -I$(GATDIR)/MGOutputMCRunProcessing -I$(GATDIR)/Analysis -I$(GATDIR)/MJDAnalysis -I$(GATDIR)/DCProcs

LIBFLAGS = -L$(MGDODIR)/lib -lMGDORoot -lMGDOBase -lMGDOTransforms -lMGDOMajorana -lMGDOGerdaTransforms -lMGDOMJDB -lMGDOTabree

LIBFLAGS += -L$(GATDIR)/lib -lGATBaseClasses -lGATMGTEventProcessing -lGATMGOutputMCRunProcessing -lGATAnalysis -lGATMJDAnalysis -lGATDCProcs $(ROOT_LIB_FLAGS) -lSpectrum -lTreePlayer -L$(TAMDIR)/lib -lTAM

include $(MGDODIR)/buildTools/BasicMakefile
