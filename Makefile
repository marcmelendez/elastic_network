all:
	gcc -Wall -Ofast -march=native -o elastic_network elastic_network.c -lm
