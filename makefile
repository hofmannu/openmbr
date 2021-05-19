CC=cuda-g++

NVCC=nvcc
CUDAFLAGS=-ccbin /usr/bin/cuda-g++ -g -Xptxas -O3 -rdc=true

CFLAGS=-O2
LINKLIBS=-lstdc++ -lhdf5 -lhdf5_cpp -lm -lcurses

ODIR = obj
_OBJS = transducerProperties.o modelMatrix.o cudaErrorCheck.o basicMathOp.o baseClass.o mbrecon.o volume.o transModel.o reconSettings.o vtkWriter.o fieldProperties.o
OBJS = $(patsubst %,$(ODIR)/%,$(_OBJS))

obj/baseClass.o :
	$(NVCC) $(CUDAFLAGS) -c src/baseClass.cpp -o obj/baseClass.o

obj/basicMathOp.o :
	$(NVCC) $(CUDAFLAGS) -c src/basicMathOp.cpp -o obj/basicMathOp.o

obj/modelMatrix.o : obj/basicMathOp.o obj/griddedData.o obj/vtkWriter.o obj/baseClass.o
	$(NVCC) $(CUDAFLAGS) -c src/modelMatrix.cpp -o obj/modelMatrix.o

obj/transducerProperties.o : obj/baseClass.o
	$(NVCC) $(CUDAFLAGS) -c src/transducerProperties.cpp -o obj/transducerProperties.o

obj/cudaErrorCheck.o :
	$(NVCC) $(CUDAFLAGS) -c src/cudaErrorCheck.cu -o obj/cudaErrorCheck.o

obj/volume.o: obj/baseClass.o
	$(NVCC) $(CUDAFLAGS) -c src/volume.cpp -o obj/volume.o

obj/griddedData.o:
	$(NVCC) $(CUDAFLAGS) -c src/griddedData.cu -o obj/griddedData.o

obj/transModel.o: obj/vtkWriter.o obj/baseClass.o obj/modelMatrix.o obj/transducerProperties.o obj/fieldProperties.o
	$(NVCC) $(CUDAFLAGS) -c src/transModel.cu -o obj/transModel.o

obj/reconSettings.o:
	$(NVCC) $(CUDAFLAGS) -c src/reconSettings.cpp -o obj/reconSettings.o

obj/vtkWriter.o: obj/griddedData.o
	$(NVCC) $(CUDAFLAGS) -c src/vtkWriter.cpp -o obj/vtkWriter.o

obj/fieldProperties.o:
	$(NVCC) $(CUDAFLAGS) -c src/fieldProperties.cpp -o obj/fieldProperties.o

obj/mbrecon.o: obj/modelMatrix.o obj/fieldProperties.o obj/transducerProperties.o obj/griddedData.o obj/cudaErrorCheck.o obj/transModel.o obj/volume.o obj/reconSettings.o obj/baseClass.o obj/basicMathOp.o
	$(NVCC) $(CUDAFLAGS) -c src/mbrecon.cu -o obj/mbrecon.o

# unit test functions
test_volume: obj/volume.o
	$(NVCC) $(CUDAFLAGS) obj/volume.OBJS unit_tests/test_volume.cpp -o bin/testVolume $(LINKLIBS)

test_transmodel: $(OBJS) 
	$(NVCC) $(CUDAFLAGS) $(OBJS) unit_tests/test_transmodel.cu -o bin/testTransmodel $(LINKLIBS)

test_mbrecon: obj/mbrecon.o
	$(NVCC) $(CUDAFLAGS) $(OBJS) unit_tests/test_mbrecon.cu -o bin/testMbrecon $(LINKLIBS)

test_vtkexport: obj/vtkWriter.o
	$(NVCC) $(CUDAFLAGS) $(OBJS) unit_tests/test_vtkexport.cpp -o bin/test_vtkexport $(LINKLIBS)

# subprocedures / helper scripts
defineNewTransducer: obj/transducerProperties.o
	$(NVCC) $(CUDAFLAGS) $(OBJS) src/defineNewTransducer.cpp -o bin/defineNewTransducer $(LINKLIBS)

# subprocedure to build a model and save it to a file
buildTransModel : obj/transModel.o obj/transducerProperties.o obj/modelMatrix.o obj/fieldProperties.o
	$(NVCC) $(CUDAFLAGS) $(OBJS) src/buildTransModel.cu -o bin/buildTransModel $(LINKLIBS)

# main procedure
all: $(OBJS) 
	$(NVCC) $(CUDAFLAGS) $(OBJS) src/main.cu -o bin/main $(LINKLIBS)

.PHONY: obj/modelMatrix.o obj/mbrecon.o obj/reconSettings.o obj/transducerProperties.o obj/basicMathOp.o obj/transModel.o

