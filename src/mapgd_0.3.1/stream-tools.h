#ifndef STREAMTOOLS_H_
#define STREAMTOOLS_H_

#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>

/// Takes a stream/string and returns a vector with one element for each field.
/* 
 * For example, if the string "A,B,C , D" is given, then {"A", "B", "C "," D"} is returned.
 * If the string "A,B,C,\nD,E,F" is given, then {"A", "B", "C"} is returned.
 */
inline std::vector<std::string> split(std::istream &in, const char &delim)
{
	std::vector<std::string> elements;
	std::string line;
	size_t last_delim=0, this_delim=0;
	getline(in, line);
	do { 
		this_delim=line.find(delim, last_delim);
		elements.push_back( line.substr(last_delim, this_delim-last_delim) );
		last_delim=this_delim+1;
	} while ( this_delim!=std::string::npos );
/*
	std::cerr << line << std::endl;
	for (size_t x=0; x<elements.size(); ++x){
		std::cerr << "[" << x << "]:" << elements[x] << std::endl;
	}
	std::cerr <<std::endl;
*/
	return elements;
}

inline std::vector<std::string> split(const std::string &s, const char &delim)
{
	std::stringstream ss;
	ss << s;
	return split(ss, delim);
}

///Default behavior is to split on white space and remove it
/*inline std::vector<std::string> split(std::istream &in)
{
	std::string line_in;
	std::string line_out;
	getline(in, line_in);
	line_out.reserve(line_in.size() );
	//TODO REMOVE WHITESPACE
	return split(line_out, ' ');
}*/

/// Takes a string and returns a vector with two elements, split on first. 
/* For example, if the string "A,B,C , D" is given, then {"A", "B,C , D"} is returned.
 */

/// An overloaded decleartion of split_first.
inline std::vector<std::string> split_first(std::istream &in, const char &delim)
{
	std::string first, second;
	std::getline(in, first, delim);
	std::getline(in, second); 
	return std::vector <std::string> {first, second};
}

inline std::vector<std::string> split_first(const std::string &s, const char &delim)
{
	std::stringstream ss;
	ss << s;
	return split_first(ss, delim);
}

/// Writes a vector out in a format . . . 
/* For example . . .
 */ 
/*
inline int write(iterator <const std:string> &key, <const std::string> &str, const char &delim)
{
	//SANITIZE STR!!
	out << key << delim << str;
}*/
#endif
