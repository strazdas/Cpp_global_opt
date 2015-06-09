build:
	@echo "==================================================================================="
	clear
	g++ gkls.c rnd_gen.c avl.cpp -o avl.out
	# g++ gkls.c rnd_gen.c main.cpp -o main.out
	# ./main.out
	./avl.out
