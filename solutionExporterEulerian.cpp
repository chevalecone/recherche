/**********************************************************************************************************
*
*   solutionExporterEulerian
*
*
*   Author  : Bruno Blais
*   File    : solutionExporterEulerian.cpp
*
*   Description : Used to house functions for exportation of solution DATA to paraview LEGACY format
*                 This class is made around structured mesh data
*
*
************************************************************************************************************/

//*******************
//  GENERAL INCLUDES
//*******************

#include <iostream>
#include <cmath>
#include <string>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>


//*******************
// HEADER INCLUDES
//*******************
#include "domain.h"
#include "solutionExporterEulerian.h"
#include "lattice.h"

//*******************
// STD USAGE
//*******************

using namespace std;

// Write scalar variable to a file
void writeScalar(Domain& domain, string label, int iter, double* scalar)
{
    ostringstream oss;

    // Write iteration number to flux
    oss << iter;

    // Establish the file name
    string filename = "output5/"+label+"_"+oss.str()+".vtk";

    // Declare the output flux and open the file
    ofstream file(filename.c_str(), ios::out | ios::trunc);
    cout << "Saving file: " << filename << endl;

    if(file)  // File has been opened correctly
    {
        file << "# vtk DataFile Version 2.0" << endl;
        file << "Test" << endl;
        file << "ASCII" << endl;
        file << "DATASET STRUCTURED_POINTS" << endl;
        file << "DIMENSIONS" << " " << domain.getNx()+1 << " " << domain.getNy()+1 <<
                " " << 1 << endl;
        file << "ORIGIN" << " " << domain.getXMin() << " " << domain.getYMin() <<
                " " << 0. << endl;
        file << "SPACING" << " " << domain.getDX() << " " << domain.getDX() <<
                " " << 1. << endl;
        file << "CELL_DATA" << " " << domain.getNTot() << endl;
        file << "SCALARS" << " " << "ScalarValue" << " " << "double" << " " << 1 <<endl;
        file << "LOOKUP_TABLE default" << endl;
        for (int i=0 ; i<domain.getNTot() ; ++i)
        {
            file <<  scalar[i] <<endl ;

        }
        file.close();  // close the file
    }
}

// Write scalar variable to a file
void writeVector(Domain& domain, string label, int iter, double** vector)
{
    ostringstream oss;

    // Write iteration number to flux
    oss << iter;

    // Establish the file name
    string filename = "output/"+label+"_"+oss.str()+".vtk";

    // Declare the output flux and open the file
    ofstream file(filename.c_str(), ios::out | ios::trunc);
    cout << "Saving file: " << filename << endl;

    if(file)  // File has been opened correctly
    {
        file << "# vtk DataFile Version 2.0" << endl;
        file << "Test" << endl;
        file << "ASCII" << endl;
        file << "DATASET STRUCTURED_POINTS" << endl;
        file << "DIMENSIONS" << " " << domain.getNx() << " " << domain.getNy() <<
                " " << 1 << endl;
        file << "ORIGIN" << " " << domain.getXMin() << " " << domain.getYMin() <<
                " " << 0. << endl;
        file << "SPACING" << " " << domain.getDX() << " " << domain.getDX() <<
                " " << 1. << endl;
        file << "POINT_DATA" << " " << domain.getNTot() << endl;
        file << "VECTORS" << " " << "Velocity" << " " << "double" <<endl;
        for (int i=0 ; i<domain.getNTot() ; ++i)
        {
            file << vector[i][0] << " " << vector[i][1] << " " << 0. << endl ;
        }


        file.close();  // close the file
    }
}

void writeLattice(Domain& domain, string label, double Kn,double mu,char name, int iter, Lattice lat)
{
    ostringstream oss;

    // Write iteration number to flux
	//oss << name;
	oss << "test_ellipse";
	//oss << "DBB_Kn_0.1";
	//oss << "DBB_Kn_0.25";
	//oss << "DBB_Kn_0.5";
	//oss << "DBB_Kn_0.75";
	//oss << "DBB_Kn_1";
	//oss << "DBB_Kn_1.5";
	oss << "_";
	//oss << "DBB";
    oss << iter;
	//oss << "_poro";
	//oss << Kn;
	//oss << "_";
	//oss << mu;
	
	//oss << name;

    // Establish the file name
    string filename = "output/"+label+"_"+oss.str()+".vtk";

    // Declare the output flux and open the file
    ofstream file(filename.c_str(), ios::out | ios::trunc);
    cout << "Saving file: " << filename << endl;

    if(file)  // File has been opened correctly
    {
        file << "# vtk DataFile Version 2.0" << endl;
        file << "Test" << endl;
        file << "ASCII" << endl;
        file << "DATASET STRUCTURED_POINTS" << endl;
        file << "DIMENSIONS" << " " << domain.getNx()+1 << " " << domain.getNy()+1 <<
                " " << 1 << endl;
        file << "ORIGIN" << " " << domain.getXMin() << " " << domain.getYMin() <<
                " " << 0. << endl;
        file << "SPACING" << " " << domain.getDX() << " " << domain.getDX() <<
                " " << 1. << endl;
        file << "CELL_DATA" << " " << domain.getNTot() << endl;
        file << "SCALARS" << " " << "Density" << " " << "double" << " " << 1 <<endl;
        file << "LOOKUP_TABLE default" << endl;
        for (int i=0 ; i<domain.getNTot() ; ++i)
        {
            file <<  lat.rho_[i] <<endl ;
        }
        /*file << "VECTORS" << " " << "Vorticity" << " " << "double" <<endl;
        for (int i=0 ; i<domain.getNTot() ; ++i)
        {
            file << lat.vorticity_[i][0] << " " << lat.vorticity_[i][1] << " " << 0. << endl ;
        }*/
        file << "VECTORS" << " " << "Velocity" << " " << "double" <<endl;
        for (int i=0 ; i<domain.getNTot() ; ++i)
        {
            file << lat.u_[i][0] << " " << lat.u_[i][1] << " " << 0. << endl ;
        }

        /*file << "SCALARS" << " " << "f0" << " " << "double" << " " << 1 <<endl;
        file << "LOOKUP_TABLE default" << endl;
        for (int i=0 ; i<domain.getNTot() ; ++i)
        {
            file <<  lat.f_[i][0] <<endl ;
        }

        file << "SCALARS" << " " << "f1" << " " << "double" << " " << 1 <<endl;
        file << "LOOKUP_TABLE default" << endl;
        for (int i=0 ; i<domain.getNTot() ; ++i)
        {
            file <<  lat.f_[i][1] <<endl ;
        }

        file << "SCALARS" << " " << "f2" << " " << "double" << " " << 1 <<endl;
        file << "LOOKUP_TABLE default" << endl;
        for (int i=0 ; i<domain.getNTot() ; ++i)
        {
            file <<  lat.f_[i][2] <<endl ;
        }

        file << "SCALARS" << " " << "f3" << " " << "double" << " " << 1 <<endl;
        file << "LOOKUP_TABLE default" << endl;
        for (int i=0 ; i<domain.getNTot() ; ++i)
        {
            file <<  lat.f_[i][3] <<endl ;
        }

        file << "SCALARS" << " " << "f4" << " " << "double" << " " << 1 <<endl;
        file << "LOOKUP_TABLE default" << endl;
        for (int i=0 ; i<domain.getNTot() ; ++i)
        {
            file <<  lat.f_[i][4] <<endl ;
        }

        file << "SCALARS" << " " << "f5" << " " << "double" << " " << 1 <<endl;
        file << "LOOKUP_TABLE default" << endl;
        for (int i=0 ; i<domain.getNTot() ; ++i)
        {
            file <<  lat.f_[i][5] <<endl ;
        }

        file << "SCALARS" << " " << "f6" << " " << "double" << " " << 1 <<endl;
        file << "LOOKUP_TABLE default" << endl;
        for (int i=0 ; i<domain.getNTot() ; ++i)
        {
            file <<  lat.f_[i][6] <<endl ;
        }

        file << "SCALARS" << " " << "f7" << " " << "double" << " " << 1 <<endl;
        file << "LOOKUP_TABLE default" << endl;
        for (int i=0 ; i<domain.getNTot() ; ++i)
        {
            file <<  lat.f_[i][7] <<endl ;
        }

        file << "SCALARS" << " " << "f8" << " " << "double" << " " << 1 <<endl;
        file << "LOOKUP_TABLE default" << endl;
        for (int i=0 ; i<domain.getNTot() ; ++i)
        {
            file <<  lat.f_[i][8] <<endl ;
        }
		
		file << "SCALARS" << " " << "m0" << " " << "double" << " " << 1 <<endl;
        file << "LOOKUP_TABLE default" << endl;
        for (int i=0 ; i<domain.getNTot() ; ++i)
        {
            file <<  lat.m_[i][0] <<endl ;
        }

        file << "SCALARS" << " " << "m1" << " " << "double" << " " << 1 <<endl;
        file << "LOOKUP_TABLE default" << endl;
        for (int i=0 ; i<domain.getNTot() ; ++i)
        {
            file <<  lat.m_[i][1] <<endl ;
        }

        file << "SCALARS" << " " << "m2" << " " << "double" << " " << 1 <<endl;
        file << "LOOKUP_TABLE default" << endl;
        for (int i=0 ; i<domain.getNTot() ; ++i)
        {
            file <<  lat.m_[i][2] <<endl ;
        }

        file << "SCALARS" << " " << "m3" << " " << "double" << " " << 1 <<endl;
        file << "LOOKUP_TABLE default" << endl;
        for (int i=0 ; i<domain.getNTot() ; ++i)
        {
            file <<  lat.m_[i][3] <<endl ;
        }

        file << "SCALARS" << " " << "m4" << " " << "double" << " " << 1 <<endl;
        file << "LOOKUP_TABLE default" << endl;
        for (int i=0 ; i<domain.getNTot() ; ++i)
        {
            file <<  lat.m_[i][4] <<endl ;
        }

        file << "SCALARS" << " " << "m5" << " " << "double" << " " << 1 <<endl;
        file << "LOOKUP_TABLE default" << endl;
        for (int i=0 ; i<domain.getNTot() ; ++i)
        {
            file <<  lat.m_[i][5] <<endl ;
        }

        file << "SCALARS" << " " << "m6" << " " << "double" << " " << 1 <<endl;
        file << "LOOKUP_TABLE default" << endl;
        for (int i=0 ; i<domain.getNTot() ; ++i)
        {
            file <<  lat.m_[i][6] <<endl ;
        }

        file << "SCALARS" << " " << "m7" << " " << "double" << " " << 1 <<endl;
        file << "LOOKUP_TABLE default" << endl;
        for (int i=0 ; i<domain.getNTot() ; ++i)
        {
            file <<  lat.m_[i][7] <<endl ;
        }

        file << "SCALARS" << " " << "m8" << " " << "double" << " " << 1 <<endl;
        file << "LOOKUP_TABLE default" << endl;
        for (int i=0 ; i<domain.getNTot() ; ++i)
        {
            file <<  lat.m_[i][8] <<endl ;
        }*/
        file.close();  // close the file
    }
}


void writeTimeScalar(string label, int taille, double** vector)
{
    ostringstream oss;

    // Write array size to flux
    oss << taille;

    // Establish the file name
    string filename = "output/"+label+"_.vtk";

    // Declare the output flux and open the file
    ofstream file(filename.c_str(), ios::out | ios::trunc);
    cout << "Saving file: " << filename << endl;

    if(file)  // File has been opened correctly
    {
        file << "# vtk DataFile Version 2.0" << endl;
        file << "Test" << endl;
        file << "ASCII" << endl;
        file << "VECTORS" << " " << label << " " << "double" <<endl;
        for (int i=0 ; i<taille ; ++i)
        {
            file << vector[i][0] << " " << vector[i][1] << " " << 0. << endl ;
        }
        file.close();  // close the file
    }
}
