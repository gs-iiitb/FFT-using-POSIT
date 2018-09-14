
default:
	gcc  fft.c -o fft.o -lm -c
	gcc fft_test.c -o fft_test.o -lm -c
	gcc  -o test fft.o fft_test.o -lm
	./test > test1.txt
