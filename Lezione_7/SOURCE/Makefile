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

restart_gas:
	mv ../OUTPUT/CONFIG/config.xyz ../INPUT/CONFIG/config.xyz
	mv ../OUTPUT/CONFIG/config_old.xyz ../INPUT/CONFIG/config_old.xyz
	mv ../OUTPUT/CONFIG/velocities.out ../INPUT/CONFIG/velocities.in
	mv ../OUTPUT/seed.out ../INPUT/seed.in
	tail -n 20 ../OUTPUT/temperature.dat >>../DATA/Gas/Restart/gas_temperature.dat
	tail -n 20 ../OUTPUT/pressure.dat >>../DATA/Gas/Restart/gas_pressure.dat
	tail -n 20 ../OUTPUT/kinetic_energy.dat >>../DATA/Gas/Restart/gas_kinetic_energy.dat
	tail -n 20 ../OUTPUT/potential_energy.dat >>../DATA/Gas/Restart/gas_potential_energy.dat
	tail -n 20 ../OUTPUT/total_energy.dat >>../DATA/Gas/Restart/gas_total_energy.dat

restart_liquid:
	mv ../OUTPUT/CONFIG/config.xyz ../INPUT/CONFIG/config.xyz
	mv ../OUTPUT/CONFIG/config_old.xyz ../INPUT/CONFIG/config_old.xyz
	mv ../OUTPUT/CONFIG/velocities.out ../INPUT/CONFIG/velocities.in
	mv ../OUTPUT/seed.out ../INPUT/seed.in
	tail -n 20 ../OUTPUT/temperature.dat >>../DATA/Liquid/Restart/liquid_temperature.dat
	tail -n 20 ../OUTPUT/pressure.dat >>../DATA/Liquid/Restart/liquid_pressure.dat
	tail -n 20 ../OUTPUT/kinetic_energy.dat >>../DATA/Liquid/Restart/liquid_kinetic_energy.dat
	tail -n 20 ../OUTPUT/potential_energy.dat >>../DATA/Liquid/Restart/liquid_potential_energy.dat
	tail -n 20 ../OUTPUT/total_energy.dat >>../DATA/Liquid/Restart/liquid_total_energy.dat
restart_solid:
	mv ../OUTPUT/CONFIG/config.xyz ../INPUT/CONFIG/config.xyz
	mv ../OUTPUT/CONFIG/config_old.xyz ../INPUT/CONFIG/config_old.xyz
	mv ../OUTPUT/CONFIG/velocities.out ../INPUT/CONFIG/velocities.in
	mv ../OUTPUT/seed.out ../INPUT/seed.in
	tail -n 20 ../OUTPUT/temperature.dat >>../DATA/Solid/Restart/solid_temperature.dat
	tail -n 20 ../OUTPUT/pressure.dat >>../DATA/Solid/Restart/solid_pressure.dat
	tail -n 20 ../OUTPUT/kinetic_energy.dat >>../DATA/Solid/Restart/solid_kinetic_energy.dat
	tail -n 20 ../OUTPUT/potential_energy.dat >>../DATA/Solid/Restart/solid_potential_energy.dat
	tail -n 20 ../OUTPUT/total_energy.dat >>../DATA/Solid/Restart/solid_total_energy.dat


termal_gas:
	mv ../OUTPUT/temperature.dat ../DATA/Gas/Termal/temperature.dat
	mv ../OUTPUT/total_energy.dat ../DATA/Gas/Termal/total_energy.dat
termal_liquid:
	mv ../OUTPUT/temperature.dat ../DATA/Liquid/Termal/temperature.dat
	mv ../OUTPUT/total_energy.dat ../DATA/Liquid/Termal/total_energy.dat
termal_solid:
	mv ../OUTPUT/temperature.dat ../DATA/Solid/Termal/temperature.dat
	mv ../OUTPUT/total_energy.dat ../DATA/Solid/Termal/total_energy.dat

gas:
	mv ../OUTPUT/temperature.dat ../DATA/Gas/temperature.dat
	mv ../OUTPUT/pressure.dat ../DATA/Gas/pressure.dat
	mv ../OUTPUT/kinetic_energy.dat ../DATA/Gas/kinetic_energy.dat
	mv ../OUTPUT/potential_energy.dat ../DATA/Gas/potential_energy.dat
	mv ../OUTPUT/total_energy.dat ../DATA/Gas/total_energy.dat
liquid:
	mv ../OUTPUT/temperature.dat ../DATA/Liquid/temperature.dat
	mv ../OUTPUT/pressure.dat ../DATA/Liquid/pressure.dat
	mv ../OUTPUT/kinetic_energy.dat ../DATA/Liquid/kinetic_energy.dat
	mv ../OUTPUT/potential_energy.dat ../DATA/Liquid/potential_energy.dat
	mv ../OUTPUT/total_energy.dat ../DATA/Liquid/total_energy.dat
solid:

	mv ../OUTPUT/temperature.dat ../DATA/Solid/temperature.dat
	mv ../OUTPUT/pressure.dat ../DATA/Solid/pressure.dat
	mv ../OUTPUT/kinetic_energy.dat ../DATA/Solid/kinetic_energy.dat
	mv ../OUTPUT/potential_energy.dat ../DATA/Solid/potential_energy.dat
	mv ../OUTPUT/total_energy.dat ../DATA/Solid/total_energy.dat

eq_ising_MRT_h0:
	mv ../OUTPUT/total_energy.dat ../DATA/Equilibration/Metropolis/total_energy.dat
eq_ising_MRT_h002:
	mv ../OUTPUT/magnetization.dat ../DATA/Equilibration/Metropolis/magnetization.dat
eq_ising_Gibbs_h0:
	mv ../OUTPUT/total_energy.dat ../DATA/Equilibration/GIBBS/total_energy.dat

eq_ising_Gibbs_h002:
	mv ../OUTPUT/magnetization.dat ../DATA/Equilibration/GIBBS/magnetization.dat
ising_metro:
	mv ../OUTPUT/magnetization.dat ../DATA/Ising/Metropolis/magnetization.dat
	mv ../OUTPUT/total_energy.dat ../DATA/Ising/Metropolis/total_energy.dat
	mv ../OUTPUT/specific_heat.dat ../DATA/Ising/Metropolis/specific_heat.dat
	mv ../OUTPUT/susceptibility.dat ../DATA/Ising/Metropolis/susceptibility.dat
	mv ../OUTPUT/acceptance.dat ../DATA/Ising/Metropolis/acceptance.dat
ising_gibbs:
	mv ../OUTPUT/magnetization.dat ../DATA/Ising/GIBBS/magnetization.dat
	mv ../OUTPUT/total_energy.dat ../DATA/Ising/GIBBS/total_energy.dat
	mv ../OUTPUT/specific_heat.dat ../DATA/Ising/GIBBS/specific_heat.dat
	mv ../OUTPUT/susceptibility.dat ../DATA/Ising/GIBBS/susceptibility.dat
	mv ../OUTPUT/acceptance.dat ../DATA/Ising/GIBBS/acceptance.dat
restart_ising:
	mv ../OUTPUT/seed.out ../INPUT/seed.in
	mv ../OUTPUT/CONFIG/config.spin ../INPUT/CONFIG/config.spin
clean_ising:
	rm ../DATA/Ising/Metropolis/*
	rm ../DATA/Ising/GIBBS/*