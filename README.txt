Команды для сборки и запуска на IBM Polus
(Таже описаны в соответствующих Makefile'ах)

===============================
1. ПОСЛЕДОВАТЕЛЬНАЯ ВЕРСИЯ
===============================

Компиляция:
g++ -O1 -Wall -Wextra -std=c++11 main.cpp -o main

Запуск:
bsub -n 1 -W 00:15 -o seq_job.out -e seq_job.err ./main

===============================
2. OpenMP ВЕРСИЯ
===============================

Компиляция:
g++ -O1 -Wall -Wextra -std=c++11 -fopenmp main.cpp -o main

Запуск для разного числа потоков:
bsub -n 2 -W 00:15 -o OpenMP_job_2.out -e OpenMP_job_2.err -R "span[hosts=1]" "OMP_NUM_THREADS=2 ./main 2"
bsub -n 4 -W 00:15 -o OpenMP_job_4.out -e OpenMP_job_4.err -R "span[hosts=1]" "OMP_NUM_THREADS=4 ./main 4"
bsub -n 8 -W 00:15 -o OpenMP_job_8.out -e OpenMP_job_8.err -R "span[hosts=1]" "OMP_NUM_THREADS=8 ./main 8"
bsub -n 16 -W 00:15 -o OpenMP_job_16.out -e OpenMP_job_16.err -R "span[hosts=1]" "OMP_NUM_THREADS=16 ./main 16"
bsub -n 32 -W 00:15 -o OpenMP_job_32.out -e OpenMP_job_32.err -R "span[hosts=1]" "OMP_NUM_THREADS=32 ./main 32"

===============================
3. MPI ВЕРСИЯ
===============================

Компиляция:
module load SpectrumMPI
mpicxx -O1 -Wall -Wextra -std=c++11 main.cpp -o main

Запуск для разного числа процессов:
mpisubmit.pl -p 2 -W 00:15 --stdout mpi_2.out --stderr mpi_2.err ./main
mpisubmit.pl -p 4 -W 00:15 --stdout mpi_4.out --stderr mpi_4.err ./main
mpisubmit.pl -p 8 -W 00:15 --stdout mpi_8.out --stderr mpi_8.err ./main
mpisubmit.pl -p 16 -W 00:15 --stdout mpi_16.out --stderr mpi_16.err ./main
mpisubmit.pl -p 32 -W 00:15 --stdout mpi_32.out --stderr mpi_32.err ./main

===============================
4. ГИБРИДНАЯ MPI+OpenMP ВЕРСИЯ
===============================

Компиляция:
module load SpectrumMPI
mpicxx -O1 -Wall -Wextra -std=c++11 -fopenmp main.cpp -o main

Запуск для 2 MPI процессов с разным числом OpenMP потоков:
bsub -n 2 -W 00:15 -o hybrid_2_1.out -e hybrid_2_1.err -R "affinity[core(1)]" "module load SpectrumMPI && OMP_NUM_THREADS=1 mpiexec ./main 1"
bsub -n 2 -W 00:15 -o hybrid_2_2.out -e hybrid_2_2.err -R "affinity[core(2)]" "module load SpectrumMPI && OMP_NUM_THREADS=2 mpiexec ./main 2"
bsub -n 2 -W 00:15 -o hybrid_2_4.out -e hybrid_2_4.err -R "affinity[core(4)]" "module load SpectrumMPI && OMP_NUM_THREADS=4 mpiexec ./main 4"
bsub -n 2 -W 00:15 -o hybrid_2_8.out -e hybrid_2_8.err -R "affinity[core(8)]" "module load SpectrumMPI && OMP_NUM_THREADS=8 mpiexec ./main 8"

Запуск для 4 MPI процессов с разным числом OpenMP потоков:
bsub -n 4 -W 00:15 -o hybrid_4_1.out -e hybrid_4_1.err -R "affinity[core(1)]" "module load SpectrumMPI && OMP_NUM_THREADS=1 mpiexec ./main 1"
bsub -n 4 -W 00:15 -o hybrid_4_2.out -e hybrid_4_2.err -R "affinity[core(2)]" "module load SpectrumMPI && OMP_NUM_THREADS=2 mpiexec ./main 2"
bsub -n 4 -W 00:15 -o hybrid_4_4.out -e hybrid_4_4.err -R "affinity[core(4)]" "module load SpectrumMPI && OMP_NUM_THREADS=4 mpiexec ./main 4"
bsub -n 4 -W 00:15 -o hybrid_4_8.out -e hybrid_4_8.err -R "affinity[core(8)]" "module load SpectrumMPI && OMP_NUM_THREADS=8 mpiexec ./main 8"
