#include <iostream>
#include "starfile_error.h"

// Object Constructor
SimpleError::SimpleError(const std::string &what, const std::string &fileArg, const long lineArg)
{
    msg = "ERROR: \n" + what;
    file= fileArg;
    line=lineArg;

    std::cerr << "in: " << file << ", line " << line << "\n";
}

// Show message
std::ostream& operator << (std::ostream& o, SimpleError& XE)
{

    o << XE.msg << std::endl;

    return o;
}
