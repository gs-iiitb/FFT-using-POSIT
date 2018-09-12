
default:
	g++ -std=c++11  fft2.cpp -o fft2.o -lm -c
	g++ -std=c++11 fft_test2.cpp -o fft_test2.o -lm -c
	g++ -std=c++11  -o test2 fft2.o fft_test2.o -lm
	./test2
