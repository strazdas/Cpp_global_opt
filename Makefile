do: run_asimpl

compile_asimpl:
	@echo "==================================================================================="
	g++ gkls.c rnd_gen.c Asimpl_main.cpp -o asimpl.out
	# g++ -std=c++11 gkls.c rnd_gen.c Asimpl_main.cpp -o asimpl.out

run_asimpl: compile_asimpl
	clear
	./asimpl.out --gkls_cls=5 --gkls_fid=84
	

test:
	g++ test.cpp -o test.out
	clear
	./test.out

run_elbme: compile_elbme
	clear
	./elbme.out


compile_elbme:
	@echo "==================================================================================="
	g++ gkls.c elbme_main.cpp -o elbme.out

cat:
	cat *.sh.po*
	cat *.sh.o*

queue:
	qstat -f
	# qdel

num:
	ls results/Disimpl-v/ | wc -l

mem_check:  compile_asimpl
	valgrind --tool=memcheck --leak-check=full -v ./main.out

profiler:  compile_asimpl
	valgrind --tool=callgrind ./asimpl.out --gkls_cls=1 --gkls_fid=1
	# git clone https://github.com/jrfonseca/gprof2dot bin/gprof2dot
	./bin/gprof2dot/gprof2dot.py -f callgrind callgrind.out.* | dot -Tsvg -o profile.svg
