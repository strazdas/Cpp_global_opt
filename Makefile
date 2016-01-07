test_worker: compile
	clear
	./worker.py -exp=2 -exe=asimpl.out &

do: run_asimpl

run_asimpl: compile
	clear
	./asimpl.out

compile: 
	@echo "==================================================================================="
	g++ gkls.c rnd_gen.c Asimpl_main.cpp -o asimpl.out
	# g++ -std=c++11 gkls.c rnd_gen.c Asimpl_main.cpp -o asimpl.out

test:
	g++ -std=c++11 test.cpp -o test.out
	clear
	./test.out

# Deprecated: add algorithm name next to command
# run: compile
# 	clear
# 	./main.out

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

# Deprecated: add algorithm name next to command
# compile: 
# 	@echo "==================================================================================="
# 	g++ gkls.c rnd_gen.c main.cpp -o main.out

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

profiler: compile
	valgrind --tool=callgrind ./asimpl.out
	# git clone https://github.com/jrfonseca/gprof2dot
	./gprof2dot/gprof2dot.py -f callgrind callgrind.out.X | dot -Tsvg -o profile.svg
