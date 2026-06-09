//===----------------------------------------------------------------------===//
//
//                                   xara
//                              https://xara.so
//
//===----------------------------------------------------------------------===//
//
// Copyright (c) 2025, OpenSees/Xara Developers
// All rights reserved.  No warranty, explicit or implicit, is provided.
//
// This source code is licensed under the BSD 2-Clause License.
// See LICENSE file or https://opensource.org/licenses/BSD-2-Clause
//
//===----------------------------------------------------------------------===//
//
// BasicAnalysisBuilder is an aggregate class which manages the analysis 
// objects:
//
// - LinearSOE
// - Domain                    *theDomain;
// - ConstraintHandler 	       *theHandler;
// - DOF_Numberer 	           *theNumberer;
// - AnalysisModel 	           *theAnalysisModel;
// - EquiSolnAlgo 	           *theAlgorithm;
// - EigenSOE 		             *theEigenSOE;
// - StaticIntegrator          *theStaticIntegrator;
// - TransientIntegrator       *theTransientIntegrator;
// - ConvergenceTest           *theTest;
//
// The BasicAnalysisBuilder assumes responsibility for
// deleting these objects, but ownership of the SOE may be
// given up.
//
// Written: cmp
//
#ifndef BasicAnalysisBulider_h
#define BasicAnalysisBulider_h

class Domain;
class ModelRegistry;
class G3_Table;
class ConstraintHandler;
class DOF_Numberer;
class AnalysisModel;
class EquiSolnAlgo;
class LinearSOE;
class EigenSOE;
class StaticIntegrator;
class TransientIntegrator;
class ConvergenceTest;
namespace OpenSees {
class LoadCase;
}

#include <string>

class BasicAnalysisBuilder
{
public:
    BasicAnalysisBuilder(ModelRegistry&);
    ~BasicAnalysisBuilder();

    enum CurrentAnalysis {
      EMPTY_ANALYSIS,
      EIGEN_ANALYSIS,
      STATIC_ANALYSIS, 
      TRANSIENT_ANALYSIS
    };

    enum Perform : int {
      Increment = 1<<0,
      Iterate   = 1<<1,
      Commit    = 1<<2,
      Update    = 1<<3
    };

    void set(ConstraintHandler*);
    void set(DOF_Numberer*);
    void set(EquiSolnAlgo*);
    void set(LinearSOE*, bool free=true);
    void unset(LinearSOE&) {
      theSOE = nullptr;
    }
    void set(StaticIntegrator&);
    void set(TransientIntegrator&, bool free=true);
    void unset(TransientIntegrator&) {
      theTransientIntegrator = nullptr;
    }
    void set(ConvergenceTest*);
    void set(EigenSOE&);
    LinearSOE* getLinearSOE();
    const EquiSolnAlgo*  getAlgorithm() const;
    StaticIntegrator*    getStaticIntegrator();
    TransientIntegrator* getTransientIntegrator();
    // for getCTestIter command
    ConvergenceTest*     getConvergenceTest();


    Domain* getDomain();
    AnalysisModel& getAnalysisModel() {
      return *theAnalysisModel;
    }

    int  initialize();
    int  domainChanged();
    int  setStaticAnalysis();
    int  setTransientAnalysis();

    //   Eigen
    void newEigenAnalysis(int typeSolver, double shift);
    int  eigen(int numMode, bool generalized, bool findSmallest);
    int  getNumEigen() {return numEigen;}

    int formUnbalance();

    // Performing analysis
    int analyze(int num_steps, double size_steps, int flag=Increment|Iterate|Commit);
    int analyzeStatic(int num_steps, int flag);
    
    int analyzeTransient(int numSteps, double dT);
    int analyzeVariable(int numSteps, double dT, double dtMin, double dtMax, int Jd);
    void Print(OPS_Stream &s, int flag);
private:
    int analyzeStep(double dT);
    int analyzeSubLevel(int level, double dT);

public:
    int analyzeGradient();
    int setGradientType(int flag);
    void wipe();

    int setLoadCase(std::string& name);
    int newLoadCase(std::string& name);

    enum CurrentAnalysis  CurrentAnalysisFlag = EMPTY_ANALYSIS;

private:
    void setLinks(CurrentAnalysis flag = EMPTY_ANALYSIS);
    void fillDefaults(enum CurrentAnalysis flag);
    int  number();

    AnalysisModel             *theAnalysisModel;
    //
    ModelRegistry&             context;
    Domain                    *theDomain;
    //
    EquiSolnAlgo              *theAlgorithm;
    StaticIntegrator          *theStaticIntegrator;
    TransientIntegrator       *theTransientIntegrator;
    //
    ConstraintHandler         *theHandler;
    DOF_Numberer              *theNumberer;
    LinearSOE                 *theSOE;
    EigenSOE                  *theEigenSOE;
    ConvergenceTest           *theTest;

    int domainStamp;
    int numEigen = 0;

    int numSubLevels = 0;
    int numSubSteps  = 0;

    bool freeSOE = true;
    bool freeTI  = true;
};

#endif
