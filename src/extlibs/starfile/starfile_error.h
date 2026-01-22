#ifndef __STARFILE_ERROR_H_
#define __STARFILE_ERROR_H_

class SimpleError
{
public:
    /** Error code */
    int __errno;

    /** Message shown */
    std::string msg;

    /** File produstd::cing the error */
    std::string file;

    /** Line number */
    long line;

    SimpleError(const std::string& what, const std::string &fileArg, const long lineArg);
    friend std::ostream& operator<<(std::ostream& o, SimpleError& XE);
};

#endif
