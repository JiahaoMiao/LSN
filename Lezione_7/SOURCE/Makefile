CC = g++
CFLAGS = -O3 --std=c++11 -g -Wall -Wextra -pedantic 
AFLAGS = -larmadillo
all: simulator.exe
simulator.exe : system.o NSL_SIMULATOR.o particle.o random.o
	$(CC) system.o particle.o random.o NSL_SIMULATOR.o -o simulator.exe $(AFLAGS)
NSL_SIMULATOR.o : NSL_SIMULATOR.cpp system.h particle.h
	$(CC) -c NSL_SIMULATOR.cpp -o NSL_SIMULATOR.o $(CFLAGS)
system.o : system.cpp system.h random.h particle.h
	$(CC) -c system.cpp -o system.o $(CFLAGS)
particle.o : particle.cpp particle.h
	$(CC) -c particle.cpp -o particle.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean:
	rm *.o simulator.exe
remove:
	rm ../OUTPUT/*.*
	rm ../OUTPUT/CONFIG/*.*

gasMC:
	mv ../OUTPUT/*.dat ../DATA/gofr/MC/Gas/
	
liquidMC:
	mv ../OUTPUT/*.dat ../DATA/gofr/MC/Liquid/

solidMC:
	mv ../OUTPUT/*.dat ../DATA/gofr/MC/Solid/

gasMD:
	mv ../OUTPUT/*.dat ../DATA/gofr/MD/Gas/

liquidMD:
	mv ../OUTPUT/*.dat ../DATA/gofr/MD/Liquid/

solidMD:
	mv ../OUTPUT/*.dat ../DATA/gofr/MD/Solid/