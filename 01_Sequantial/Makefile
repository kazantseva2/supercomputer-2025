CFLAGS = -o1 -Wall -Wextra -std=c++11
CC = g++

compile: 
	$(CC) $(CFLAGS) main.cpp -o main

compile_check: 
	$(CC) $(CFLAGS) main_check.cpp -o main_check

run: compile
	./main

run_check: compile_check
	./main_check

run_polus: compile
	for N in 1; do \
		echo "Submitting sequential job"; \
		bsub -J "seq_job" \
		     -n 1 \
		     -W 00:15 \
		     -q short \
		     -o "seq_job.out" \
		     -e "seq_job.err" \
		     "./main"; \
	done

clean:
	rm -f main main_check *.png *.txt *.out *.err 