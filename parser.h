/*****************************************
*	Class Parser
*	Author : Bruno Blais
*
*	Role : This class takes care of parsing the number of mesh element
* 			in the x and y direction using the tinyxml library
*			Overly complicated for the goal, but this is, in a way
*			a demonstration
*
******************************************/


#include "tinyxml2.h"
#include <iostream> 
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <regex>
#include <sstream>

#ifndef PARSER_H
#define PARSER_H

class  Parser
{
 public:
    Parser();                       // Public constructor
    ~Parser(); 	                    // Public destructor
    int parse();		    // Parse the xml document
    double getDx() {return dx_;};
    double getTimeEnd() {return timeEnd_;};
	double getError() {return error_;};
    int    getOutputFrequency() {return outputFrequency_;};
    double getXmin() {return xmin_;};
    double getXmax() {return xmax_;};
    double getYmin() {return ymin_;};
    double getYmax() {return ymax_;};
    double getTau() {return tau_;};
    double getRho0() {return rho0_;};
    double getMu() {return mu_;};

 private:
        void readDomain();
        void readTime();
        void readMesh();
        void readPhysicalProperties();
        double dx_;    // number of x nodes
        double timeEnd_;
		double error_;
        double mu_;
        double rho0_;
        double tau_;
        int outputFrequency_;
        double xmin_;
        double xmax_;
        double ymin_;
        double ymax_;
        tinyxml2::XMLDocument doc_; // document being read
};

// Templated classes to detect if a number has what it takes to be an int or a double
// Serves as a sanity check for the xml parsing 
template<typename T>
bool isInt(T x)
{
    std::string s;
    std::regex e ("^-?\\d+");
    std::stringstream ss; 
    ss << x;
    ss >>s;
    if (std::regex_match (s,e)) return true;
    else return false;
}

template<typename T>
bool isDouble(T x)
{
    std::string s;
    std::regex e ("^-?\\d*\\.?\\d+");
    std::stringstream ss; 
    ss << x;
    ss >>s;
    if (std::regex_match (s,e)) return true;
    else return false;
}

#endif // PARSER DEFINITION
