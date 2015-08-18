do: run

# test:
# 	g++ test.cpp -o test.out
# 	clear
# 	./test.out

run: compile
	clear
	./main.out

run_mpi: compile_mpi
	qsub -pe orte 1 cls1.sh
	qsub -pe orte 1 cls2.sh
	qsub -pe orte 1 cls3.sh
	qsub -pe orte 1 cls4.sh
	qsub -pe orte 16 cls5.sh
	qsub -pe orte 16 cls6.sh
	qsub -pe orte 16 cls7.sh
	qsub -pe orte 16 cls8.sh

run_elbme: compile_elbme
	clear
	./elbme.out

compile:
	@echo "==================================================================================="
	g++ gkls.c rnd_gen.c main.cpp -o main.out

compile_mpi:
	/opt/openmpi/bin/mpiCC -o mpi_main.out gkls.c rnd_gen.c mpi_main.cpp

compile_elbme:
	@echo "==================================================================================="
	g++ gkls.c elbme_main.cpp -o elbme.out

clean:
	rm -f *.sh.po* *.sh.o* core.*

cat:
	cat *.sh.po*
	cat *.sh.o*

queue:
	qstat -f
	# qdel

num:
	ls results/Disimpl-v/ | wc -l

mem_check: compile
	valgrind --tool=memcheck --leak-check=full -v ./main.out
