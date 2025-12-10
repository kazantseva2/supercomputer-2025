CFLAGS = -o1 -Wall -Wextra -std=c++11
CC = mpic++ 
CC_polus = mpicxx

# Локальная компиляция
compile: 
	$(CC) $(CFLAGS) main.cpp -o main

compile_polus:
	module load SpectrumMPI && $(CC_polus) $(CFLAGS) main.cpp -o main

run_1: compile
	mpirun -np 1 ./main

run_2: compile
	mpirun -np 2 ./main

run_4: compile
	mpirun -np 4 ./main

run_8: compile
	mpirun -np 8 -map-by :OVERSUBSCRIBE ./main


# === Запуск на Polus с mpisubmit.pl ===
run_polus: compile_polus
	mpisubmit.pl -p 1 -W 00:15 --stdout mpi_1.out --stderr mpi_1.err ./main

run_polus_all: compile_polus
	for N in 2 4 8 16 32; do \
		echo "Submitting MPI job with $$N processes"; \
		mpisubmit.pl -p $$N -W 00:15 --stdout mpi_$$N.out --stderr mpi_$$N.err ./main; \
	done

# === Показать результаты MPI заданий ===
show:
	@for N in 2 4 8 16 32; do \
		if [ -f "mpi_$$N.out" ]; then \
			echo "=== Results for $$N MPI processes ==="; \
			cat "mpi_$$N.out"; \
			echo ""; \
		else \
			echo "File mpi_$$N.out not found"; \
		fi; \
	done

# === Показать ошибки MPI заданий ===
show_errors:
	@for N in 2 4 8 16 32; do \
		if [ -f "mpi_$$N.err" ]; then \
			echo "=== Errors for $$N MPI processes ==="; \
			cat "mpi_$$N.err"; \
			echo ""; \
		fi; \
	done

clean:
	rm -f main *.png *.txt *.btr *.out *.err