CC=g++
#CFLAGS=-Wall -std=c++11
CFLAGS= -std=c++17

run: PA0-b
	./PA0-b 

main: PA0-b.cpp
	$(CC) $(CFLAGS) PA0-b.cpp -o PA0-b

clean:
	rm -f PA0-b

test: PA0-b
	./PA0-b NLDM_lib_max2Inp cleaned_iscas89_99_circuits/c1908_.isc 1020 8 478 8001

	./PA0-b NLDM_lib_max2Inp cleaned_iscas89_99_circuits/b22.isc 55 8731 29797

.PHONY: test_PA0-b.cpp
