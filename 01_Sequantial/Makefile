CFLAGS = -O3 -Wall -Wextra -std=c++11
CC = g++

compile: 
	$(CC) $(CFLAGS) main.cpp -o main

compile_check: 
	$(CC) $(CFLAGS) main_check.cpp -o main_check

run: compile
	./main

run_check: compile_check
	./main_check

clean:
	rm -f main main_check solution.txt