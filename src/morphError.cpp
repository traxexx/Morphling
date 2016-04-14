#include "morphError.h"
#include "Globals.h"
#include <iostream>
#include <time.h> // log time

using std::cout;
using std::cerr;
using std::endl;
using std::string;

void morphMessageNoTime(string & message)
{
	cout << message << endl;
	LOG << message << endl;
}

void morphMessageNoTime(const char* message)
{
	cout << message << endl;
	LOG << message << endl;
}

void morphMessage( string & message )
{
	time_t raw_time;
	time(&raw_time);
	cout << message << ".    " << ctime(&raw_time) << endl;
	LOG << message << ".    " << ctime(&raw_time) << endl;
}

void morphMessage( const char* message )
{
	time_t raw_time;
	time(&raw_time);
	cout << message << ".    " << ctime(&raw_time) << endl;
	LOG << message << ".   " << ctime(&raw_time) << endl;
}

void morphWarning( string & message )
{
	cerr << "Warning: " << message << "!" << endl;
	LOG << "Warning: " << message << "!" << endl;
}

void morphWarning( const char* message )
{
	cerr << "Warning: " << message << "!" << endl;
	LOG << "Warning: " << message << "!" << endl;
}

void morphError( string & message, int code )
{
	cerr << "ERROR: " << message << "!\n" << endl;
	LOG << "ERROR: " << message << "!\n" << endl;
	exit(code);
}

void morphError( const char * message ) {
	cerr << "ERROR: " << message << "!\n" << endl;
	LOG << "ERROR: " << message << "!\n" << endl;
	exit(1);
}

void morphErrorFile( string & filename )
{
	cerr << "ERROR: cannot open " << filename << "!" << endl;
	LOG << "ERROR: cannot open " << filename << "!" << endl;
	exit(1);
}