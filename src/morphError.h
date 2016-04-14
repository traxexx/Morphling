#include <string>

void morphMessageNoTime(std::string & message);

void morphMessageNoTime(const char* message);

void morphMessage( std::string & message );

void morphMessage( const char* message );

void morphWarning( std::string & message );

void morphWarning( const char* message );

void morphError( std::string & message, int code );

void morphError( const char * message );

void morphErrorFile( std::string & filename );