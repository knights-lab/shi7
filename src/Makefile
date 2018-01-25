CC=gcc
.PHONY : all

all: shi7_trimmer gotta_split

clean:
	rm shi7_trimmer
	rm gotta_split

shi7_trimmer: shi7_trimmer.c ; $(CC) -m64 -std=gnu99 -O3 -fwhole-program -o shi7_trimmer shi7_trimmer.c

gotta_split: gotta_split_swapped.c ; $(CC) -m64 -std=gnu99 -O3 -o gotta_split gotta_split_swapped.c
