//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation
//
//===----------------------------------------------------------------------===//
//
// written: cmp
//
#include <stdio.h>
// #include <unordered_map>
// #include <string>
// #include <vector>

#include <tcl.h>
#include <runtimeAPI.h>

// typedef std::unordered_map<std::string, std::vector<std::string>> G3_Config;

class Domain;
class BasicModelBuilder;

class AnalysisModel;
class ConstraintHandler;
class LinearSOE;
class EigenSOE;
class DOF_Numberer;


class G3_Runtime {
public:

  Tcl_Interp     *m_interp = nullptr;

// MODEL BUILDING
  BasicModelBuilder *m_builder = nullptr;
  Domain            *m_domain  = nullptr;
  bool            model_is_built=false;

// ANALYSIS
  AnalysisModel  *m_analysis_model     = nullptr;
  AnalysisModel **m_analysis_model_ptr = &m_analysis_model;

// IO
  FILE* streams[3] = {stdin,stdout,stderr};
};



