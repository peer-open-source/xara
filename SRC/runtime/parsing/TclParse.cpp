#include <string>

struct Tcl_Interp {};
#define TCL_OK 0
#define TCL_ERROR 1

int 
Tcl_GetDouble(Tcl_Interp *interp, const char *arg, double *value)
{
    std::string str(arg);
    try {
        *value = std::stod(str);
        return TCL_OK;
    } catch (...) {
        return TCL_ERROR;
    }
}

int 
Tcl_GetInt(Tcl_Interp *interp, const char *arg, int *value)
{
    std::string str(arg);
    try {
        *value = std::stoi(str);
        return TCL_OK;
    } catch (...) {
        return TCL_ERROR;
    }
}

int 
Tcl_SplitList(Tcl_Interp *interp, const char *list, int *argcPtr, char ***argvPtr)
{
    std::string str(list);
    std::vector<std::string> elements;
    size_t pos = 0;
    while ((pos = str.find(' ')) != std::string::npos) {
        elements.push_back(str.substr(0, pos));
        str.erase(0, pos + 1);
    }
    if (!str.empty()) {
        elements.push_back(str);
    }

    *argcPtr = static_cast<int>(elements.size());
    *argvPtr = new char*[*argcPtr];
    for (int i = 0; i < *argcPtr; ++i) {
        (*argvPtr)[i] = new char[elements[i].size() + 1];
        std::strcpy((*argvPtr)[i], elements[i].c_str());
    }
    return TCL_OK;
}