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
	mv ../OUTPUT/temperature.dat ../DATA/gofr/MC/Gas/temperature.dat
	mv ../OUTPUT/pressure.dat ../DATA/gofr/MC/Gas/pressure.dat
	mv ../OUTPUT/kinetic_energy.dat ../DATA/gofr/MC/Gas/kinetic_energy.dat
	mv ../OUTPUT/potential_energy.dat ../DATA/gofr/MC/Gas/potential_energy.dat
	mv ../OUTPUT/total_energy.dat ../DATA/gofr/MC/Gas/total_energy.dat
	mv ../OUTPUT/gofr.dat ../DATA/gofr/MC/Gas/gofr.dat
	mv ../OUTPUT/output.dat ../DATA/gofr/MC/Gas/output.dat
	mv ../OUTPUT/acceptance.dat ../DATA/gofr/MC/Gas/acceptance.dat
	
liquidMC:
	mv ../OUTPUT/temperature.dat ../DATA/gofr/MC/Liquid/temperature.dat
	mv ../OUTPUT/pressure.dat ../DATA/gofr/MC/Liquid/pressure.dat
	mv ../OUTPUT/kinetic_energy.dat ../DATA/gofr/MC/Liquid/kinetic_energy.dat
	mv ../OUTPUT/potential_energy.dat ../DATA/gofr/MC/Liquid/potential_energy.dat
	mv ../OUTPUT/total_energy.dat ../DATA/gofr/MC/Liquid/total_energy.dat
	mv ../OUTPUT/gofr.dat ../DATA/gofr/MC/Liquid/gofr.dat
	mv ../OUTPUT/output.dat ../DATA/gofr/MC/Liquid/output.dat
	mv ../OUTPUT/acceptance.dat ../DATA/gofr/MC/Liquid/acceptance.dat

solidMC:
	mv ../OUTPUT/temperature.dat ../DATA/gofr/MC/Solid/temperature.dat
	mv ../OUTPUT/pressure.dat ../DATA/gofr/MC/Solid/pressure.dat
	mv ../OUTPUT/kinetic_energy.dat ../DATA/gofr/MC/Solid/kinetic_energy.dat
	mv ../OUTPUT/potential_energy.dat ../DATA/gofr/MC/Solid/potential_energy.dat
	mv ../OUTPUT/total_energy.dat ../DATA/gofr/MC/Solid/total_energy.dat
	mv ../OUTPUT/gofr.dat ../DATA/gofr/MC/Solid/gofr.dat
	mv ../OUTPUT/output.dat ../DATA/gofr/MC/Solid/output.dat
	mv ../OUTPUT/acceptance.dat ../DATA/gofr/MC/Solid/acceptance.dat

gasMD:
	mv ../OUTPUT/temperature.dat ../DATA/gofr/MD/Gas/temperature.dat
	mv ../OUTPUT/pressure.dat ../DATA/gofr/MD/Gas/pressure.dat
	mv ../OUTPUT/kinetic_energy.dat ../DATA/gofr/MD/Gas/kinetic_energy.dat
	mv ../OUTPUT/potential_energy.dat ../DATA/gofr/MD/Gas/potential_energy.dat
	mv ../OUTPUT/total_energy.dat ../DATA/gofr/MD/Gas/total_energy.dat
	mv ../OUTPUT/gofr.dat ../DATA/gofr/MD/Gas/gofr.dat
	mv ../OUTPUT/acceptance.dat ../DATA/gofr/MD/Gas/acceptance.dat
	mv ../OUTPUT/output.dat ../DATA/gofr/MD/Gas/output.dat

liquidMD:
	mv ../OUTPUT/temperature.dat ../DATA/gofr/MD/Liquid/temperature.dat
	mv ../OUTPUT/pressure.dat ../DATA/gofr/MD/Liquid/pressure.dat
	mv ../OUTPUT/kinetic_energy.dat ../DATA/gofr/MD/Liquid/kinetic_energy.dat
	mv ../OUTPUT/potential_energy.dat ../DATA/gofr/MD/Liquid/potential_energy.dat
	mv ../OUTPUT/total_energy.dat ../DATA/gofr/MD/Liquid/total_energy.dat
	mv ../OUTPUT/gofr.dat ../DATA/gofr/MD/Liquid/gofr.dat
	mv ../OUTPUT/acceptance.dat ../DATA/gofr/MD/Liquid/acceptance.dat
	mv ../OUTPUT/output.dat ../DATA/gofr/MD/Liquid/output.dat

solidMD:
	mv ../OUTPUT/temperature.dat ../DATA/gofr/MD/Solid/temperature.dat
	mv ../OUTPUT/pressure.dat ../DATA/gofr/MD/Solid/pressure.dat
	mv ../OUTPUT/kinetic_energy.dat ../DATA/gofr/MD/Solid/kinetic_energy.dat
	mv ../OUTPUT/potential_energy.dat ../DATA/gofr/MD/Solid/potential_energy.dat
	mv ../OUTPUT/total_energy.dat ../DATA/gofr/MD/Solid/total_energy.dat
	mv ../OUTPUT/gofr.dat ../DATA/gofr/MD/Solid/gofr.dat
	mv ../OUTPUT/acceptance.dat ../DATA/gofr/MD/Solid/acceptance.dat
	mv ../OUTPUT/output.dat ../DATA/gofr/MD/Solid/output.dat