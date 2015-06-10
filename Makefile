build:
	@echo "==================================================================================="
	clear
	g++ gkls.c rnd_gen.c main.cpp -o main.out
	./main.out
