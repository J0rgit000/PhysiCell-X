VERSION := $(shell grep . ../../../VERSION.txt 2>/dev/null | cut -f1 -d:)

PROGRAM_NAME_REF := mechanics_ref
PROGRAM_NAME_VAL := mechanics_val

CC := mpic++
ifdef PHYSICELL_CPP
	CC := $(PHYSICELL_CPP)
endif

ARCH := native
CFLAGS := -march=$(ARCH) -O0 -fomit-frame-pointer -mfpmath=both -m64 -std=c++17 -g -fopenmp
COMPILE_COMMAND := $(CC) $(CFLAGS)

DistPhy_OBJECTS := DistPhy_Environment.o DistPhy_Cartesian.o DistPhy_Utils.o DistPhy_Collective.o

BioFVM_OBJECTS := BioFVM_vector.o BioFVM_mesh.o BioFVM_microenvironment.o BioFVM_solvers.o BioFVM_matlab.o \
BioFVM_utilities.o BioFVM_basic_agent.o BioFVM_MultiCellDS.o BioFVM_agent_container.o

PhysiCell_core_OBJECTS := PhysiCell_phenotype.o PhysiCell_cell_container.o PhysiCell_standard_models.o PhysiCell_cell.o \
PhysiCell_custom.o PhysiCell_utilities.o PhysiCell_constants.o PhysiCell_cell_mpi.o PhysiCell_signal_behavior.o \
PhysiCell_rules.o PhysiCell_basic_signaling.o PhysiCell_cell_container_mechanics.o PhysiCell_standard_models_mpi.o

PhysiCell_module_OBJECTS := PhysiCell_SVG.o PhysiCell_pathology.o PhysiCell_MultiCellDS.o PhysiCell_various_outputs.o \
PhysiCell_pugixml.o PhysiCell_settings.o

PhysiCell_custom_module_OBJECTS := custom.o snapshot_writer.o

PhysiCell_OBJECTS := $(BioFVM_OBJECTS) pugixml.o $(PhysiCell_core_OBJECTS) $(PhysiCell_module_OBJECTS) $(DistPhy_OBJECTS)
ALL_OBJECTS := $(PhysiCell_OBJECTS) $(PhysiCell_custom_module_OBJECTS)

all: reference validate

reference: main_serial.cpp $(ALL_OBJECTS)
	$(COMPILE_COMMAND) -o $(PROGRAM_NAME_REF) $(ALL_OBJECTS) main_serial.cpp

validate: main.cpp $(ALL_OBJECTS)
	$(COMPILE_COMMAND) -o $(PROGRAM_NAME_VAL) $(ALL_OBJECTS) main.cpp

DistPhy_Environment.o: ../../../DistPhy/DistPhy_Environment.cpp
	$(COMPILE_COMMAND) -c ../../../DistPhy/DistPhy_Environment.cpp

DistPhy_Cartesian.o: ../../../DistPhy/DistPhy_Cartesian.cpp
	$(COMPILE_COMMAND) -c ../../../DistPhy/DistPhy_Cartesian.cpp

DistPhy_Utils.o: ../../../DistPhy/DistPhy_Utils.cpp
	$(COMPILE_COMMAND) -c ../../../DistPhy/DistPhy_Utils.cpp

DistPhy_Collective.o: ../../../DistPhy/DistPhy_Collective.cpp
	$(COMPILE_COMMAND) -c ../../../DistPhy/DistPhy_Collective.cpp

PhysiCell_phenotype.o: ../../../core/PhysiCell_phenotype.cpp
	$(COMPILE_COMMAND) -c ../../../core/PhysiCell_phenotype.cpp

PhysiCell_digital_cell_line.o: ../../../core/PhysiCell_digital_cell_line.cpp
	$(COMPILE_COMMAND) -c ../../../core/PhysiCell_digital_cell_line.cpp

PhysiCell_cell.o: ../../../core/PhysiCell_cell.cpp
	$(COMPILE_COMMAND) -c ../../../core/PhysiCell_cell.cpp

PhysiCell_cell_mpi.o: ../../../core/PhysiCell_cell_mpi.cpp
	$(COMPILE_COMMAND) -c ../../../core/PhysiCell_cell_mpi.cpp

PhysiCell_cell_container.o: ../../../core/PhysiCell_cell_container.cpp
	$(COMPILE_COMMAND) -c ../../../core/PhysiCell_cell_container.cpp

PhysiCell_standard_models.o: ../../../core/PhysiCell_standard_models.cpp
	$(COMPILE_COMMAND) -c ../../../core/PhysiCell_standard_models.cpp

PhysiCell_utilities.o: ../../../core/PhysiCell_utilities.cpp
	$(COMPILE_COMMAND) -c ../../../core/PhysiCell_utilities.cpp

PhysiCell_custom.o: ../../../core/PhysiCell_custom.cpp
	$(COMPILE_COMMAND) -c ../../../core/PhysiCell_custom.cpp

PhysiCell_constants.o: ../../../core/PhysiCell_constants.cpp
	$(COMPILE_COMMAND) -c ../../../core/PhysiCell_constants.cpp

PhysiCell_signal_behavior.o: ../../../core/PhysiCell_signal_behavior.cpp
	$(COMPILE_COMMAND) -c ../../../core/PhysiCell_signal_behavior.cpp

PhysiCell_rules.o: ../../../core/PhysiCell_rules.cpp
	$(COMPILE_COMMAND) -c ../../../core/PhysiCell_rules.cpp

PhysiCell_basic_signaling.o: ../../../core/PhysiCell_basic_signaling.cpp
	$(COMPILE_COMMAND) -c ../../../core/PhysiCell_basic_signaling.cpp

PhysiCell_cell_container_mechanics.o: ../../../core/PhysiCell_cell_container_mechanics.cpp
	$(COMPILE_COMMAND) -c ../../../core/PhysiCell_cell_container_mechanics.cpp

PhysiCell_standard_models_mpi.o: ../../../core/PhysiCell_standard_models_mpi.cpp
	$(COMPILE_COMMAND) -c ../../../core/PhysiCell_standard_models_mpi.cpp

BioFVM_vector.o: ../../../BioFVM/BioFVM_vector.cpp
	$(COMPILE_COMMAND) -c ../../../BioFVM/BioFVM_vector.cpp

BioFVM_agent_container.o: ../../../BioFVM/BioFVM_agent_container.cpp
	$(COMPILE_COMMAND) -c ../../../BioFVM/BioFVM_agent_container.cpp

BioFVM_mesh.o: ../../../BioFVM/BioFVM_mesh.cpp
	$(COMPILE_COMMAND) -c ../../../BioFVM/BioFVM_mesh.cpp

BioFVM_microenvironment.o: ../../../BioFVM/BioFVM_microenvironment.cpp
	$(COMPILE_COMMAND) -c ../../../BioFVM/BioFVM_microenvironment.cpp

BioFVM_solvers.o: ../../../BioFVM/BioFVM_solvers.cpp
	$(COMPILE_COMMAND) -c ../../../BioFVM/BioFVM_solvers.cpp

BioFVM_utilities.o: ../../../BioFVM/BioFVM_utilities.cpp
	$(COMPILE_COMMAND) -c ../../../BioFVM/BioFVM_utilities.cpp

BioFVM_basic_agent.o: ../../../BioFVM/BioFVM_basic_agent.cpp
	$(COMPILE_COMMAND) -c ../../../BioFVM/BioFVM_basic_agent.cpp

BioFVM_matlab.o: ../../../BioFVM/BioFVM_matlab.cpp
	$(COMPILE_COMMAND) -c ../../../BioFVM/BioFVM_matlab.cpp

BioFVM_MultiCellDS.o: ../../../BioFVM/BioFVM_MultiCellDS.cpp
	$(COMPILE_COMMAND) -c ../../../BioFVM/BioFVM_MultiCellDS.cpp

pugixml.o: ../../../BioFVM/pugixml.cpp
	$(COMPILE_COMMAND) -c ../../../BioFVM/pugixml.cpp

PhysiCell_SVG.o: ../../../modules/PhysiCell_SVG.cpp
	$(COMPILE_COMMAND) -c ../../../modules/PhysiCell_SVG.cpp

PhysiCell_pathology.o: ../../../modules/PhysiCell_pathology.cpp
	$(COMPILE_COMMAND) -c ../../../modules/PhysiCell_pathology.cpp

PhysiCell_MultiCellDS.o: ../../../modules/PhysiCell_MultiCellDS.cpp
	$(COMPILE_COMMAND) -c ../../../modules/PhysiCell_MultiCellDS.cpp

PhysiCell_various_outputs.o: ../../../modules/PhysiCell_various_outputs.cpp
	$(COMPILE_COMMAND) -c ../../../modules/PhysiCell_various_outputs.cpp

PhysiCell_pugixml.o: ../../../modules/PhysiCell_pugixml.cpp
	$(COMPILE_COMMAND) -c ../../../modules/PhysiCell_pugixml.cpp

PhysiCell_settings.o: ../../../modules/PhysiCell_settings.cpp
	$(COMPILE_COMMAND) -c ../../../modules/PhysiCell_settings.cpp

custom.o: ./custom_modules/custom.cpp
	$(COMPILE_COMMAND) -c ./custom_modules/custom.cpp

snapshot_writer.o: ./custom_modules/snapshot_writer.cpp
	$(COMPILE_COMMAND) -c ./custom_modules/snapshot_writer.cpp

clean:
	rm -f *.o
	rm -f $(PROGRAM_NAME_REF) $(PROGRAM_NAME_VAL)

data-cleanup:
	rm -f *.mat
	rm -f *.xml
	rm -f *.svg
	rm -f ./output/*
	touch ./output/empty.txt

checkpoint:
	zip -r $$(date +%b_%d_%Y_%H%M).zip Makefile *.cpp *.h config/*.xml custom_modules/*
