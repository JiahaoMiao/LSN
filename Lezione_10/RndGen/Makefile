CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.exe : main.o random.o
	$(CC) random.o main.o -o main.exe
main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
random_old.o : random_old.cpp random_old.h
	$(CC) -c random_old.cpp -o random_old.o $(CFLAGS)
clean :
	rm *.o main.exe seed.out
test.o: test.cpp
	$(CC) -c test.cpp -o test.o $(CFLAGS)
test.exe : test.o random.o random_old.o
	$(CC) random.o random_old.o test.o -o test.exe