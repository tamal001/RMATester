all: RMATester

RMATester: RMATester.cpp
	g++ RMATester.cpp -o RMATester -std=c++14

clean: rm -f RMATester