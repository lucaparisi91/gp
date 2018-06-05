FFTW_DIR=~/software/fftw

all: tools.h GPMain.cpp tools.cpp GP.h
	g++  fourier.cpp tools.cpp GPMain.cpp  -I $(FFTW_DIR)/include -L $(FFTW_DIR)/lib -g -lfftw3 -lm -o gp

test: testFourier.cpp
	g++ tools.cpp testFourier.cpp -I $(FFTW_DIR)/include -L $(FFTW_DIR)/lib -g -lfftw3 -lm -o test
