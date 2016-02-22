#include <iostream>
#include <fstream>
#include <stdlib.h>

namespace std {

void read_input_file() {


	string line;

	ifstream file ("../INPUT/input.pfs"); 
	if(!file.is_open()) {
		cout << "Error opening input file\n";
		exit(1);
	}

	if(!getline(file,line))
		cout << "Nothing in input file\n";

	else {
		while(getline(file,line)) {
	
			cout << "Hi\n";
	
			if (line.compare("block limits") != 0)
	    			cout << line << " is not " << "block limits" << '\n';
	
	
		}
	}
	file.close();
}

}
