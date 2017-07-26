#include "tinyxml2.h"
#include "parser.h"
#include <stdlib.h>

// Constructor
Parser::Parser()
{

}

// Destructor
Parser::~Parser()
{
}

int Parser::parse()
{
    // Load the file
    doc_.LoadFile( "resources/parameter.xml" );

    readMesh();
    readDomain();
    readTime();
    readPhysicalProperties();

    return doc_.ErrorID();
}

void Parser::readDomain()
{
    tinyxml2::XMLElement* root = doc_.FirstChildElement("MAIN" );
    if (root == NULL)
    {
        std::cout << "Root element of XML file MAIN is not present \n";
        std::cout << "XML Reader will crash \n";
        abort();
    }
    // Pointer to mesh sub-category of the XML file
    tinyxml2::XMLElement* domainInfo = root->FirstChildElement("domain");
    if (domainInfo == NULL)
    {
        std::cout << "A domain information category is necessary \n" << std::endl;
        abort();
    }

    xmin_ = atof(domainInfo->FirstChildElement("xmin")->GetText());
    std::cout << "Value of xMin : " << xmin_ << std::endl;
    
    xmax_ = atof(domainInfo->FirstChildElement("xmax")->GetText());
    std::cout << "Value of xMax : " << xmax_ << std::endl;
    
    ymin_ = atof(domainInfo->FirstChildElement("ymin")->GetText());
    std::cout << "Value of yMin : " << ymin_ << std::endl;

    ymax_ = atof(domainInfo->FirstChildElement("ymax")->GetText());
    std::cout << "Value of yMax : " << ymax_ << std::endl;

}

void Parser::readMesh()
{
    tinyxml2::XMLElement* root = doc_.FirstChildElement("MAIN" );
    if (root == NULL)
    {
        std::cout << "Root element of XML file MAIN is not present \n";
        std::cout << "XML Reader will crash \n";
        abort();
    }
    // Pointer to mesh sub-category of the XML file
    tinyxml2::XMLElement* mesh = root->FirstChildElement("mesh");
    if (mesh == NULL)
    {
        std::cout << "A mesh information category is necessary \n" << std::endl;
        abort();
    }

    dx_ = atof(mesh->FirstChildElement("dx")->GetText());
    std::cout << "Dx value of  : " << dx_ << std::endl;
}

void Parser::readTime()
{
    tinyxml2::XMLElement* root = doc_.FirstChildElement("MAIN" );
    if (root == NULL)
    {
        std::cout << "Root element of XML file MAIN is not present \n";
        std::cout << "XML Reader will crash \n";
        abort();
    }
    // Pointer to mesh sub-category of the XML file
    tinyxml2::XMLElement* timeInfo = root->FirstChildElement("time");
    if (timeInfo== NULL)
    {
        std::cout << "A time information category is necessary \n" << std::endl;
        abort();
    }

    timeEnd_ = atof(timeInfo->FirstChildElement("timeEnd")->GetText());
    std::cout << "Endtime value of : " << timeEnd_ << std::endl;
	
	error_ = atof(timeInfo->FirstChildElement("error")->GetText());
	std::cout << "Error value of : " << error_ << std::endl;

    outputFrequency_ = atoi(timeInfo->FirstChildElement("outputFrequency")->GetText());
    std::cout << "Output frequency value of : " << outputFrequency_ << std::endl;
}

void Parser::readPhysicalProperties()
{

    tinyxml2::XMLElement* root = doc_.FirstChildElement("MAIN" );
    if (root == NULL)
    {
        std::cout << "Root element of XML file MAIN is not present \n";
        std::cout << "XML Reader will crash \n";
        abort();
    }
    // Pointer to mesh sub-category of the XML file
    tinyxml2::XMLElement* propInfo = root->FirstChildElement("properties");
    if (propInfo== NULL)
    {
        std::cout << "A property information category is necessary \n" << std::endl;
        abort();
    }

    mu_ = atof(propInfo->FirstChildElement("mu")->GetText());
    std::cout << "mu value of : " << mu_ << std::endl;

    rho0_ = atof(propInfo->FirstChildElement("rho0")->GetText());
    std::cout << "Rho0 value of : " << rho0_ << std::endl;

    tau_ = atof(propInfo->FirstChildElement("tau")->GetText());
    std::cout << "Tau value of : " << tau_ << std::endl;
}

