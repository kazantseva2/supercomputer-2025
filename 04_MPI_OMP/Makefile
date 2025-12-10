CFLAGS = -o1 -Wall -Wextra -std=c++11 -fopenmp
CC = mpic++ 
CC_polus = mpicxx


compile: 
	$(CC) $(CFLAGS) main.cpp -o main

compile_polus:
	module load SpectrumMPI && $(CC_polus) $(CFLAGS) main.cpp -o main

run_1_4: compile
	mpirun -np 1 ./main 4

run_2_4: compile
	mpirun -np 2 ./main 4

# === Запуск на Polus с mpisubmit.pl ===

run_polus_2_all: compile_polus
	for THREADS in 1 2 4 8; do \
		echo "Submitting 2 MPI processes * $$THREADS threads"; \
		bsub -n 2 \
		     -W 00:15 \
		     -o "hybrid_2_$${THREADS}.out" \
		     -e "hybrid_2_$${THREADS}.err" \
		     -R "affinity[core($$THREADS)]" \
		     "module load SpectrumMPI && OMP_NUM_THREADS=$$THREADS mpiexec ./main $$THREADS"; \
		sleep 1; \
	done


run_polus_4_all: compile_polus
	for THREADS in 1 2 4 8; do \
		echo "Submitting 2 MPI processes * $$THREADS threads"; \
		bsub -n 4 \
		     -W 00:15 \
		     -o "hybrid_4_$${THREADS}.out" \
		     -e "hybrid_4_$${THREADS}.err" \
		     -R "affinity[core($$THREADS)]" \
		     "module load SpectrumMPI && OMP_NUM_THREADS=$$THREADS mpiexec ./main $$THREADS"; \
		sleep 1; \
	done


# === Показать результаты гибридных заданий ===
show_hybrid:
	@for PROCS in 2 4; do \
		for THREADS in 1 2 4 8; do \
			FILE="hybrid_$${PROCS}_$${THREADS}.out"; \
			if [ -f "$$FILE" ]; then \
				echo "=== Results for $$PROCS MPI processes * $$THREADS OpenMP threads ==="; \
				cat "$$FILE"; \
				echo ""; \
			else \
				echo "File $$FILE not found"; \
			fi; \
		done; \
	done

show_hybrid_errors:
	@for PROCS in 2 4; do \
		for THREADS in 1 2 4 8; do \
			FILE="hybrid_$${PROCS}_$${THREADS}.err"; \
			if [ -f "$$FILE" ]; then \
				echo "=== Errors for $$PROCS MPI processes * $$THREADS OpenMP threads ==="; \
				cat "$$FILE"; \
				echo ""; \
			fi; \
		done; \
	done

clean:
	rm -f main *.png *.txt *.btr *.out *.err 