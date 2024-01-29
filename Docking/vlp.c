//to compile: gcc vlp.c -o vlp

#include <string.h>
#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <stdlib.h>

//Evaluates to true if string 'in' starts with string literal 'pref' 
#define matchesPrefix(in, pref) 	(!strncmp(in, pref, sizeof(pref)/sizeof(pref[0])-1))

//Evaluates to true if char c is a decimal digit in ASCII
#define isDigit(c)	( c >= '0' && c <= '9' )

//This function takes a string, and returns a pointer to the first character of
//the next word, as separated by one or more whitespace
char* nextWord(const char * s) {

	//nullptrchk
	if ( !s ) return (char*)s;

	//check for end of string
	if ( !(*s) ) return (char*)s;

	//advance to whitespace, but not past null
	while ( *s != ' ' && *s != '\0') s++;
	if ( !(*s) ) return (char*)s;

	//advance past whitespace
	while ( *s == ' ') s++;

	//cast to explicitly drop const
	return (char*)s;
}

int main(int argc, char* argv[]) {

	//Open file
	FILE* infile = fopen(argv[1], "r");
	if(!infile) {
		perror("Unable to open file!");
		return 1;
	}

	//Buffer for current line
	char line[1024];
	//Iterator over line
	const char* i;

	//State variable; are we currently within the energy/rmsd table?
	bool inTable = false;

	//Accumulators
	uint16_t confs = 0;
	double avg_energy = 0.0;
	double avg_rsmd = 0.0;
	double min_energy =  100000000.0;
	double max_energy = -100000000.0;
	//Temporary helper variable to save on atof() calls
	double tempnum;

	//Loop over lines of file
	while(fgets(line, sizeof(line), infile) != NULL) {

		if ( !inTable ) {
			inTable = matchesPrefix(line, "-----+------------+----------+----------");
			continue;
		}

		if ( line[0] != ' ' && !isDigit(line[0]) ) break;

		//i is now at start of line
		i=line;

		//i is now at configuration id number
		i = nextWord(i);

		//i is now at energy
		i = nextWord(i);
		tempnum = atof(i);
		if ( tempnum < min_energy ) min_energy = tempnum;
		if ( tempnum > max_energy ) max_energy = tempnum;
		avg_energy += tempnum;

		//i is now at rsmd l.b.
		i = nextWord(i);
		tempnum = atof(i);

		//i is now at rsmd u.b.
		i = nextWord(i);
		tempnum += atof(i);
		avg_rsmd += tempnum / 2;

		confs++;

	}

	if ( confs == 0 ) {
		puts("No configurations were found");
		return 1;
	}

	avg_energy /= confs;
	avg_rsmd /= confs;

	printf("%f\t%f\t%f\t%f\n", avg_energy, max_energy, min_energy, avg_rsmd);

	return 0;

}