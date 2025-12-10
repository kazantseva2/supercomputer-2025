CFLAGS = -o1 -Wall -Wextra -std=c++11 -fopenmp
CC = g++

compile: 
	$(CC) $(CFLAGS) main.cpp -o main

run_1: compile
	./main 1

run_4: compile
	./main 4

run_8: compile
	./main 8

run_16: compile
	./main 16

run_check: compile_check
	./main_check

# === Запуск на Polus ===
run_polus: compile
	for N in 2 4 8 16 32; do \
		M=$$(( ($$N + 7) / 8 )); \
		echo "Submitting job with $$N threads and $$M cores"; \
		bsub -J "OpenMP_job_$$N" \
		     -n $$M \
		     -W 00:15 \
		     -q short \
		     -o "OpenMP_job_$$N.out" \
		     -e "OpenMP_job_$$N.err" \
		     -R "span[hosts=1]" \
		     "OMP_NUM_THREADS=$$N ./main $$N"; \
	done

# === Показать результаты (out файлы) ===
show_results:
	@for N in 1 2 4 8 16 32; do \
		if [ -f "OpenMP_job_$$N.out" ]; then \
			echo "=== Results for $$N threads ==="; \
			cat "OpenMP_job_$$N.out"; \
			echo ""; \
		else \
			echo "File OpenMP_job_$$N.out not found"; \
		fi; \
	done

clean:
	rm -f main *.out *.err *.lsf *.txt