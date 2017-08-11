/**********************************************************************************************************
*
*   LBM BURST CODE
*   
*   
*   Author  : Bruno Blais
*   File    : readwrite.h
*   Date    : November 9th 2012    
*   Began   : November 8th 2012
*
*   Description : This is the header for the reading and writing functions function
*
************************************************************************************************************/
#include <string>
#include "domain.h"
#include "lattice.h"

#ifndef SOLUTIONEXPORTEREULERIAN
#define SOLUTIONEXPORTEREULERIAN // Prevents multiple definition of the same class
  void writeScalar(Domain&, std::string, int, double*); // write a scalar variable
  void writeVector(Domain&, std::string, int, double**); // write a vector variable
  void writeLattice(Domain&, std::string, int, Lattice); // write a vector variable
  void writeTimeScalar(std::string,int, double**);
#endif
