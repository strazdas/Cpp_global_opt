run:
	@echo "==================================================================================="
	clear
	g++ gkls.c rnd_gen.c main.cpp -o main.out
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

compile_mpi:
	/opt/openmpi/bin/mpiCC -o mpi_main.out gkls.c rnd_gen.c mpi_main.cpp

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
