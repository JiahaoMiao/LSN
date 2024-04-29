CC = g++
CFLAGS = -Wall -O3 -g --std=c++11

%.exe: %.o wave_func.o RndGen/random.o
	$(CC) -o $@ $^ $(CFLAGS)
%.o : %.cpp
	$(CC) -o $@ -c $<  $(CFLAGS)
clean:
	rm *.exe
format:
	clang-format -i -style=google *.cpp
clean_output:
	rm ES/*.dat
	rm GS/*.dat