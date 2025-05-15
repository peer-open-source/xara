//===----------------------------------------------------------------------===//
//
//                                   xara
//
//===----------------------------------------------------------------------===//
//                              https://xara.so
//===----------------------------------------------------------------------===//
// written: cmp
//
#if 0
#include <string>
#include <vector>
#include "G3_Runtime.h"

#include <LinearSOE.h>
#include <EigenSOE.h>
#include <DOF_Numberer.h>
#include <ConstraintHandler.h>
#include <StaticIntegrator.h>
#include <TransientIntegrator.h>
#include <EquiSolnAlgo.h>

// DEFAULTS
#include <AnalysisModel.h>
#include <ProfileSPDLinSolver.h>
#include <ProfileSPDLinDirectSolver.h>
#include <NewtonRaphson.h>
#include <TransformationConstraintHandler.h>
#include <CTestNormUnbalance.h>
#include <ProfileSPDLinSOE.h>
#include <PlainHandler.h>
#include <Newmark.h>
#include <RCM.h>
#include <LoadControl.h>

class StaticIntegrator;

#define G3Config_keyExists(conf, key) ((conf).find((key)) != (conf).end())

template <typename T>
using G3_Parse = T* (*)(G3_Runtime*, int, const char **const);


// Wrap a function with signature
//     T* G3Parse_newT(G3_Runtime*, int argc, G3_Char** argv)
// so that it works with a std::vector<std::string>
template<typename T, T* (*fn)(G3_Runtime*, int, G3_Char **)>
T* G3Object_newParsed(G3_Runtime *rt, G3_Char* command, std::vector<std::string> args) {
    std::vector<G3_Char *> cstrs;
    cstrs.reserve(args.size()+1);
    cstrs.push_back(command);
    for (auto &s : args)
      cstrs.push_back(const_cast<char *>(s.c_str()));
    return (*fn)(rt, cstrs.size(), cstrs.data());
}

#endif