a.out: gilbert.c
	gcc -o a.out gilbert.c -lgsl -lgslcblas -lm -O3

