CC = g++
CFLAGS = -Wall -O3 --std=c++11

%.exe: %.o RndGen/random.o
	$(CC) -o $@ $^ $(CFLAGS) && ./$@
%.o : %.cpp
	$(CC) -o $@ -c $<  $(CFLAGS)
clean:
	rm *.exe
format:
	clang-format -i -style=google *.cpp