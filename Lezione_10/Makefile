CC = mpicxx
CFLAGS = -Wall -O3 -g --std=c++17

main.exe: main.o TSP.o RndGen/random.o
	$(CC) -o $@ $^ $(CFLAGS)
test.exe: test.o TSP.o RndGen/random.o
	$(CC) -o $@ $^ $(CFLAGS)
%.exe: %.o simulation.o RndGen/random.o
	$(CC) -o $@ $^ $(CFLAGS)
%.o : %.cpp
	$(CC) -o $@ -c $<  $(CFLAGS)
clean:
	rm *.exe *.o
format:
	clang-format -i -style=google *.cpp