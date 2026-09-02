//===----------------------------------------------------------------------===//
//
//                                   xara
//                              https://xara.so
//
//===----------------------------------------------------------------------===//
//
// Copyright (c) 2025, Claudio M. Perez
// All rights reserved.  No warranty, explicit or implicit, is provided.
//
// This source code is licensed under the BSD 2-Clause License.
// See LICENSE file or https://opensource.org/licenses/BSD-2-Clause
//
//===----------------------------------------------------------------------===//
//
// Purpose: This file contains the class definition for XaraClassBroker.
// XaraClassBroker is is an object broker class that is meant to become
// a threadsafe replacement for the BrokerAllClasses class.
// All methods are virtual to allow for subclasses; which can be
// used by programmers when introducing new subclasses of the main objects.
//
//===----------------------------------------------------------------------===//
//

#include <Hash.h>
using namespace OpenSees::Hash;
using namespace OpenSees::Hash::literals;

#define DISPATCH(symbol) case  hasher<std::string>()(#symbol): return new symbol();

#include "packages.h"
#include <XaraClassBroker.h>

// Convergence tests
#include "analysis/criteria/CTestNormUnbalance.h"
#include "analysis/criteria/CTestRelativeNormUnbalance.h"
#include "analysis/criteria/CTestNormDispIncr.h"
#include "analysis/criteria/CTestRelativeNormDispIncr.h"
#include "analysis/criteria/CTestRelativeTotalNormDispIncr.h"
#include "analysis/criteria/CTestEnergyIncr.h"
#include "analysis/criteria/CTestRelativeEnergyIncr.h"
#include "analysis/criteria/CTestFixedNumIter.h"

// graph numbering schemes
#include "graph/numberer/RCM.h"
#include "graph/numberer/SimpleNumberer.h"

// uniaxial material model header files

// PY springs: RWBoulanger and BJeremic
#include "PY/PySimple1.h"


// Sections
#include "ElasticSection2d.h"
#include "ElasticSection3d.h"
#include "GenericSection1d.h"
//#include "GenericSectionNd.h"
#include "SectionAggregator.h"
//#include "FiberSection.h"
#include "FiberSection2d.h"
#include "FiberSection3d.h"
#include "FiberSectionAsym3d.h" //Xinlong Du
#include "ElasticMembranePlateSection.h"
#include "MembranePlateFiberSection.h"
#include "Bidirectional.h"
#include "LayeredShellFiberSection.h" // Yuli Huang & Xinzheng Lu

// NDMaterials
// start Yuli Huang & Xinzheng L
#include "PlateRebarMaterial.h"
#include "PlateFromPlaneStressMaterial.h"
#include "PlaneStressUserMaterial.h"
// end Yuli Huang & Xinzheng Lu
#include "feap/FeapMaterial03.h"
#include "CycLiqCP3D.h"
#include "CycLiqCPPlaneStrain.h"
#include "CycLiqCPSP3D.h"
#include "CycLiqCPSPPlaneStrain.h"

#include "soil/FluidSolidPorousMaterial.h"
#include "soil/PressureDependMultiYield.h"
#include "soil/PressureDependMultiYield02.h"
#include "soil/PressureIndependMultiYield.h"
//
// element header files
//
#include "Element.h"
#include "Truss/Truss.h"
#include "Truss/Truss2.h"
#include "Truss/TrussSection.h"
#include "Truss/CorotTruss.h"
#include "Truss/CorotTrussSection.h"
#include "Truss/InertiaTruss.h"
#include "Point/ZeroLength.h"
#include "Point/ZeroLengthSection.h"
#include "Point/ZeroLengthContact2D.h"
#include "Point/ZeroLengthContact3D.h"
#include "Point/ZeroLengthContactNTS2D.h"
#include "Point/ZeroLengthInterface2D.h"
#include "Point/ZeroLengthContactASDimplex.h"
//#include "ZeroLengthND.h"

#include "Plane/FourNodeQuad.h"
#include "Plane/EnhancedQuad.h"
#include "Plane/NineNodeMixedQuad.h"
#include "Plane/NineNodeQuad.h"
#include "Plane/EightNodeQuad.h"
#include "Plane/ConstantPressureVolumeQuad.h"
#include "Plane/Nine_Four_Node_QuadUP.h"
#include "Plane/BBarFourNodeQuadUP.h"
#include "Plane/FourNodeQuadUP.h"
#include "Plane/Tri31.h"

#include "Frame/Elastic/ElasticBeam2d.h"
#include "Frame/Elastic/ElasticBeam3d.h"
#include "Frame/Elastic/ModElasticBeam2d.h" //SAJalali
#include "Frame/Elastic/ModElasticBeam3d.h"
#include "Frame/Elastic/ElasticTimoshenkoBeam2d.h"
#include "Frame/Elastic/ElasticTimoshenkoBeam3d.h"
#include "Frame/Other/Force/ForceBeamColumn2d.h"
#include "Frame/Other/Force/ForceBeamColumn3d.h"
#include "Frame/Other/Displ/DispBeamColumn2d.h"
#include "Frame/Other/Displ/DispBeamColumn3d.h"
#include "Frame/Other/Displ/DispBeamColumnAsym3d.h"   // Xinlong Du
#include "Frame/Other/Mixed/MixedBeamColumnAsym3d.h"  // Xinlong Du


#include "UWelements/SSPquad.h"
#include "UWelements/SSPquadUP.h"
#include "UWelements/SSPbrick.h"
#include "UWelements/SSPbrickUP.h"
#include "UWelements/BeamContact2D.h"
#include "UWelements/BeamContact2Dp.h"
#include "UWelements/BeamContact3D.h"
#include "UWelements/BeamContact3Dp.h"
#include "UWelements/BeamEndContact3D.h"
#include "UWelements/BeamEndContact3Dp.h"
#include "UWelements/QuadBeamEmbedContact.h"

#include "Other/PML/PML2D.h"
#include "Other/PML/PML3D.h"

#include "Shell/ShellMITC4.h"
#include "Shell/ShellMITC9.h"
#include "Shell/ShellDKGQ.h"
#include "Shell/ShellNLDKGQ.h"
#include "Shell/ASDShellQ4.h"
#include "Brick/Brick.h"
#include "Brick/BbarBrick.h"
#include "Brick/BrickUP.h"
#include "Brick/BBarBrickUP.h"
#include "Brick/Twenty_Eight_Node_BrickUP.h"
#include "Joint/Joint2D.h"
#include "Link/TwoNodeLink.h"
#include "Link/LinearElasticSpring.h"
#include "Link/Inerter.h"

// Bearings
#include "Bearing/elastomeric/ElastomericBearingBoucWen2d.h"
#include "Bearing/elastomeric/ElastomericBearingBoucWen3d.h"
#include "Bearing/elastomeric/ElastomericBearingPlasticity2d.h"
#include "Bearing/elastomeric/ElastomericBearingPlasticity3d.h"
#include "Bearing/elastomeric/ElastomericBearingUFRP2d.h"
#include "Bearing/elastomeric/ElastomericX.h"
#include "Bearing/elastomeric/HDR.h"
#include "Bearing/elastomeric/LeadRubberX.h"
#include "Bearing/friction/FlatSliderSimple2d.h"
#include "Bearing/friction/FlatSliderSimple3d.h"
#include "Bearing/friction/FPBearingPTV.h"
#include "Bearing/friction/RJWatsonEQS2d.h"
#include "Bearing/friction/RJWatsonEQS3d.h"
#include "Bearing/friction/SingleFPSimple2d.h"
#include "Bearing/friction/SingleFPSimple3d.h"
#include "Bearing/friction/TripleFrictionPendulum.h"


#include "mvlem/MVLEM.h"       
#include "mvlem/SFI_MVLEM.h"   
#include "mvlem/MVLEM_3D.h"    
#include "mvlem/SFI_MVLEM_3D.h"

#include "Boundary/RockingBC.h"

#include "ASDEA/CEqElement/ASDEmbeddedNodeElement.h"
#include "ASDEA/absorbentBoundaries/ASDAbsorbingBoundary2D.h"
#include "ASDEA/absorbentBoundaries/ASDAbsorbingBoundary3D.h"

#include "LinearCrdTransf2d.h"
#include "LinearCrdTransf3d.h"
#include "PDeltaCrdTransf2d.h"
#include "PDeltaCrdTransf3d.h"
#include "CorotCrdTransf2d.h"
#include "CorotCrdTransf3d.h"

#include "quadrature/Frame/HingeMidpointBeamIntegration.h"
#include "quadrature/Frame/HingeEndpointBeamIntegration.h"
#include "quadrature/Frame/HingeRadauBeamIntegration.h"
#include "quadrature/Frame/HingeRadauTwoBeamIntegration.h"
#include "quadrature/Frame/UserDefinedHingeIntegration.h"
#include "quadrature/Frame/DistHingeIntegration.h"
#include "quadrature/Frame/RegularizedHingeIntegration.h"

#include "quadrature/Frame/LobattoBeamIntegration.h"
#include "quadrature/Frame/LegendreBeamIntegration.h"
#include "quadrature/Frame/RadauBeamIntegration.h"
#include "quadrature/Frame/NewtonCotesBeamIntegration.h"
#include "quadrature/Frame/TrapezoidalBeamIntegration.h"
#include "quadrature/Frame/UserDefinedBeamIntegration.h"
#include "quadrature/Frame/FixedLocationBeamIntegration.h"
#include "quadrature/Frame/LowOrderBeamIntegration.h"
#include "quadrature/Frame/MidDistanceBeamIntegration.h"
#include "quadrature/Frame/CompositeSimpsonBeamIntegration.h"

// Node header files
#include "Node.h"
#ifdef HEAP_NODE
#include "HeapNode.h"
#endif

#include "FileStream.h"
#include "StandardStream.h"
#include "XmlFileStream.h"
#include "DataFileStream.h"
#include "DataFileStreamAdd.h"
#include "BinaryFileStream.h"
#include "DatabaseStream.h"
#include "DummyStream.h"

#include "NodeRecorder.h"
#include "ElementRecorder.h"
#include "EnvelopeNodeRecorder.h"
#include "EnvelopeElementRecorder.h"
#include "DriftRecorder.h"
//#include "MPCORecorder.h"
#include "VTK_Recorder.h"
#include "GmshRecorder.h"

// mp_constraint header files
#include "MP_Constraint.h"
#include "Joint/MP_Joint2D.h"

// sp_constraint header files
#include "SP_Constraint.h"
#include "SP_Constraint.h"
#include "ImposedMotionSP.h"
#include "ImposedMotionSP1.h"

// Pressure_Constraint header file
#include "Pressure_Constraint.h"

// nodal load header files
#include "NodalLoad.h"

// elemental load header files
#include "ElementalLoad.h"
#include "Beam2dUniformLoad.h"
#include "Beam2dPointLoad.h"
#include "Beam3dUniformLoad.h"
#include "Beam3dPointLoad.h"
#include "BrickSelfWeight.h"
#include "SelfWeight.h"
#include "SurfaceLoader.h"

// matrix, vector & id header files
#include "Matrix.h"
#include "Vector.h"
#include "ID.h"

// subdomain header files
#include "Subdomain.h"

// constraint handler header files
#include "ConstraintHandler.h"
#include "PlainHandler.h"
#include "PenaltyConstraintHandler.h"
#include "LagrangeConstraintHandler.h"
#include "TransformationConstraintHandler.h"


// equi soln algo
#include "EquiSolnAlgo.h"
#include "Linear.h"
#include "NewtonRaphson.h"
#include "Broyden.h"
#include "NewtonLineSearch.h"
#include "ModifiedNewton.h"


#include "BisectionLineSearch.h"
#include "InitialInterpolatedLineSearch.h"
#include "RegulaFalsiLineSearch.h"
#include "SecantLineSearch.h"

// domain decomp soln algo header files
#include "DomainDecompAlgo.h"

// integrator header files
#include "DisplacementControl.h"
#ifdef _PARALLEL_PROCESSING
#include "DistributedDisplacementControl.h"
#endif
#include "LoadControl.h"


// system of eqn header files
#include "LinearSOE.h"
#include "DomainSolver.h"
#include "fullGEN/FullGenLinSOE.h"
#include "bandGEN/BandGenLinSOE.h"
#include "bandSPD/BandSPDLinSOE.h"
#include "profileSPD/ProfileSPDLinSOE.h"
#include "profileSPD/ProfileSPDLinSubstrSolver.h"
#include "sparseGEN/SparseGenColLinSOE.h"
#include "DomainDecompositionAnalysis.h"

// load patterns
#include "StaticPattern.h"
#include "UniformExcitation.h"
#include "MultiSupportPattern.h"
#include "GroundMotion.h"
#include "InterpolatedGroundMotion.h"
#ifdef OPSDEF_DRM
#  include "drm/DRMLoadPatternWrapper.h"
#endif // OPSDEF_DRM

#ifdef _H5DRM
#  include "drm/H5DRM.h"
#endif

#include "Parameter.h"
#include "ElementParameter.h"
#include "MaterialStageParameter.h"
#include "MatParameter.h"
#include "InitialStateParameter.h"
#include "ElementStateParameter.h"

// time series
#include "LinearSeries.h"
#include "PathSeries.h"
#include "PathTimeSeries.h"
#include "RectangularSeries.h"
#include "ConstantSeries.h"
#include "TrigSeries.h"
#include "TriangleSeries.h"

// time series integrators
#include "TrapezoidalTimeSeriesIntegrator.h"

#include "eigenSOE/ArpackSOE.h"

#ifdef _PETSC
#  include "PetscSOE.h"
#  include "PetscSolver.h"
#  include "SparseGenColLinSOE.h"
#endif

#ifdef _MUMPS
#  include "MumpsSOE.h"
#  ifdef _PARALLEL_PROCESSING
#  include "MumpsParallelSOE.h"
#  endif
#endif

#ifdef _PARALLEL_PROCESSING
#  include "DistributedBandSPDLinSOE.h"
#  include "DistributedProfileSPDLinSOE.h"
#  include "DistributedSparseGenColLinSOE.h"
#  include "DistributedSparseGenRowLinSOE.h"
#  include "DistributedBandGenLinSOE.h"
#  include "DistributedSuperLU.h"
#  include "ParallelNumberer.h"
#  include "StaticDomainDecompositionAnalysis.h"
#  include "TransientDomainDecompositionAnalysis.h"
#  include "DistributedDiagonalSOE.h"
#endif

using namespace OpenSees;

typedef struct uniaxialPackage {
  int classTag;
  char *libName;
  char *funcName;
  UniaxialMaterial *(*funcPtr)(void);
  struct uniaxialPackage *next;
} UniaxialPackage;

static UniaxialPackage *theUniaxialPackage = NULL;

XaraClassBroker::XaraClassBroker() : lastDomainSolver(0) {}

XaraClassBroker::~XaraClassBroker() {}

Actor *
XaraClassBroker::getNewActor(int classTag, Channel *theChannel)
{
  return nullptr;
}


GraphNumberer *
XaraClassBroker::getPtrNewGraphNumberer(int classTag)
{
  switch (classTag) {
  case GraphNUMBERER_TAG_RCM:
    return new RCM();

  case GraphNUMBERER_TAG_SimpleNumberer:
    return new SimpleNumberer();

  default:
    opserr << "XaraClassBroker::getPtrNewGraphNumberer - ";
    opserr << " - no GraphNumberer type exists for class tag ";
    opserr << classTag << "\n";
    return 0;
  }
}

/*****************************************
 *
 * METHODS TO GET NEW MODELLING CLASSES
 *
 *****************************************/

Element *
XaraClassBroker::getNewElement(int classTag)
{
  return nullptr;
}

Node *
XaraClassBroker::getNewNode(int classTag)
{
  return nullptr;
}

MP_Constraint *
XaraClassBroker::getNewMP(int classTag)
{
  switch (classTag) {
  case CNSTRNT_TAG_MP_Constraint:
    return new MP_Constraint(classTag);

  case CNSTRNT_TAG_MP_Joint2D:
    return new MP_Joint2D();

  default:
    opserr << "XaraClassBroker::getNewMP - ";
    opserr << " - no MP_Constraint type exists for class tag ";
    opserr << classTag << "\n";
    return 0;
  }
}

SP_Constraint *
XaraClassBroker::getNewSP(int classTag)
{
  return nullptr;
}

Pressure_Constraint *
XaraClassBroker::getNewPC(int classTag)
{
  switch (classTag) {
  case CNSTRNT_TAG_Pressure_Constraint:
    return new Pressure_Constraint(classTag);

  default:
    opserr << "XaraClassBroker::getNewPC - ";
    opserr << " - no Pressure_Constraint type exists for class tag ";
    opserr << classTag << "\n";
    return 0;
  }
}


ElementalLoad *
XaraClassBroker::getNewElementalLoad(int classTag)
{
  return nullptr;
}

CrdTransf *
XaraClassBroker::getNewCrdTransf(int classTag)
{
  switch (classTag) {
  case CRDTR_TAG_LinearCrdTransf2d:
    return new LinearCrdTransf2d();
  case CRDTR_TAG_PDeltaCrdTransf2d:
    return new PDeltaCrdTransf2d();
  case CRDTR_TAG_CorotCrdTransf2d:
    return new CorotCrdTransf2d();
  case CRDTR_TAG_LinearCrdTransf3d:
    return new LinearCrdTransf3d();
  case CRDTR_TAG_PDeltaCrdTransf3d:
    return new PDeltaCrdTransf3d();
  default:
    opserr << "XaraClassBroker::getCrdTransf - ";
    opserr << " - no CrdTransf type exists for class tag ";
    opserr << classTag << "\n";
    return 0;
  }
}

BeamIntegration *
XaraClassBroker::getNewBeamIntegration(int classTag)
{
  switch (classTag) {
  case BEAM_INTEGRATION_TAG_Lobatto:
    return new LobattoBeamIntegration();

  case BEAM_INTEGRATION_TAG_Legendre:
    return new LegendreBeamIntegration();

  case BEAM_INTEGRATION_TAG_Radau:
    return new RadauBeamIntegration();

  case BEAM_INTEGRATION_TAG_NewtonCotes:
    return new NewtonCotesBeamIntegration();

  case BEAM_INTEGRATION_TAG_Trapezoidal:
    return new TrapezoidalBeamIntegration();

  case BEAM_INTEGRATION_TAG_UserDefined:
    return new UserDefinedBeamIntegration();

  case BEAM_INTEGRATION_TAG_FixedLocation:
    return new FixedLocationBeamIntegration();

  case BEAM_INTEGRATION_TAG_LowOrder:
    return new LowOrderBeamIntegration();

  case BEAM_INTEGRATION_TAG_MidDistance:
    return new MidDistanceBeamIntegration();

  case BEAM_INTEGRATION_TAG_CompositeSimpson:
    return new CompositeSimpsonBeamIntegration();

  case BEAM_INTEGRATION_TAG_HingeMidpoint:
    return new HingeMidpointBeamIntegration();

  case BEAM_INTEGRATION_TAG_HingeRadau:
    return new HingeRadauBeamIntegration();

  case BEAM_INTEGRATION_TAG_HingeRadauTwo:
    return new HingeRadauTwoBeamIntegration();

  case BEAM_INTEGRATION_TAG_HingeEndpoint:
    return new HingeEndpointBeamIntegration();

  case BEAM_INTEGRATION_TAG_UserHinge:
    return new UserDefinedHingeIntegration();

  case BEAM_INTEGRATION_TAG_DistHinge:
    return new DistHingeIntegration();

  case BEAM_INTEGRATION_TAG_RegularizedHinge:
    return new RegularizedHingeIntegration();

  default:
    opserr << "XaraClassBroker::getBeamIntegration - ";
    opserr << " - no BeamIntegration type exists for class tag ";
    opserr << classTag << "\n";
    return 0;
  }
}

UniaxialMaterial *
XaraClassBroker::getNewUniaxialMaterial(int classTag)
{
  return nullptr;
}

SectionForceDeformation *
XaraClassBroker::getNewSection(int classTag)
{
  return nullptr;
}

NDMaterial *
XaraClassBroker::getNewNDMaterial(int classTag)
{
  return nullptr;
}

Fiber *
XaraClassBroker::getNewFiber(int classTag)
{
  return nullptr;
}


FrictionModel *
XaraClassBroker::getNewFrictionModel(int classTag)
{
  return nullptr;
}

ConvergenceTest *
XaraClassBroker::getNewConvergenceTest(int classTag)
{
  return nullptr;
}

LoadPattern *
XaraClassBroker::getNewLoadPattern(int classTag)
{
  switch (classTag) {

  case PATTERN_TAG_MultiSupportPattern:
    return new MultiSupportPattern();

#ifdef OPSDEF_DRM
  case PATTERN_TAG_DRMLoadPattern:
    return new DRMLoadPatternWrapper();
#endif // OPSDEF_DRM

#ifdef _H5DRM
  case PATTERN_TAG_H5DRM:
    return new H5DRM();
#endif
  default:
    opserr << "XaraClassBroker::getPtrLoadPattern - ";
    opserr << " - no Load type exists for class tag ";
    opserr << classTag << "\n";
    return 0;
  }
}

GroundMotion *
XaraClassBroker::getNewGroundMotion(int classTag)
{
  switch (classTag) {

  case GROUND_MOTION_TAG_GroundMotion:
    return new GroundMotion(GROUND_MOTION_TAG_GroundMotion);

  case GROUND_MOTION_TAG_InterpolatedGroundMotion:
    return new GroundMotion(GROUND_MOTION_TAG_InterpolatedGroundMotion);

  default:
    opserr << "XaraClassBroker::getPtrGroundMotion - ";
    opserr << " - no Load type exists for class tag ";
    opserr << classTag << "\n";
    return 0;
  }
}

TimeSeries *
XaraClassBroker::getNewTimeSeries(int classTag)
{
  switch (classTag) {
  case TSERIES_TAG_LinearSeries:
    return new LinearSeries;

  case TSERIES_TAG_RectangularSeries:
    return new RectangularSeries;

  case TSERIES_TAG_PathTimeSeries:
    return new PathTimeSeries;

  case TSERIES_TAG_PathSeries:
    return new PathSeries;

  case TSERIES_TAG_ConstantSeries:
    return new ConstantSeries;

  case TSERIES_TAG_TriangleSeries:
    return new TriangleSeries;

  case TSERIES_TAG_TrigSeries:
    return new TrigSeries;

  default:
    opserr << "XaraClassBroker::getPtrTimeSeries - ";
    opserr << " - no Load type exists for class tag ";
    opserr << classTag << "\n";
    return 0;
  }
}

TimeSeriesIntegrator *
XaraClassBroker::getNewTimeSeriesIntegrator(int classTag)
{
  switch (classTag) {
  case TIMESERIES_INTEGRATOR_TAG_Trapezoidal:
    return new TrapezoidalTimeSeriesIntegrator();

  default:
    opserr << "XaraClassBroker::getPtrTimeSeriesIntegrator - ";
    opserr << " - no Load type exists for class tag ";
    opserr << classTag << "\n";
    return 0;
  }
}

Matrix *
XaraClassBroker::getPtrNewMatrix(int classTag, int noRows, int noCols)
{
  switch (classTag) {
  case MATRIX_TAG_Matrix:
    return new Matrix(noRows, noCols);

  default:
    opserr << "XaraClassBroker::getPtrNewMatrix - ";
    opserr << " - no NodalLoad type exists for class tag ";
    opserr << classTag << "\n";
    return 0;
  }
}

Vector *
XaraClassBroker::getPtrNewVector(int classTag, int size)
{
  switch (classTag) {
  case VECTOR_TAG_Vector:
    return new Vector(size);

  default:
    opserr << "XaraClassBroker::getPtrNewVector - ";
    opserr << " - no Vector type exists for class tag ";
    opserr << classTag << "\n";
    return 0;
  }
}

ID *
XaraClassBroker::getPtrNewID(int classTag, int size)
{
  switch (classTag) {
  case ID_TAG_ID:
    return new ID(size);

  default:
    opserr << "XaraClassBroker::getPtrNewID - ";
    opserr << " - no ID type exists for class tag ";
    opserr << classTag << "\n";
    return 0;
  }
}

/*****************************************
 *
 * METHODS TO GET NEW OUTPUT CLASS OBJECTS
 *
 *****************************************/

OPS_Stream *
XaraClassBroker::getPtrNewStream(int classTag)
{
  switch (classTag) {
  case OPS_STREAM_TAGS_StandardStream:
    return new StandardStream();

  case OPS_STREAM_TAGS_FileStream:
    return new FileStream();

  case OPS_STREAM_TAGS_XmlFileStream:
    return new XmlFileStream();

  case OPS_STREAM_TAGS_DataFileStream:
    return new DataFileStream();

  case OPS_STREAM_TAGS_DataFileStreamAdd:
    return new DataFileStreamAdd();

  case OPS_STREAM_TAGS_BinaryFileStream:
    return new BinaryFileStream();

  case OPS_STREAM_TAGS_DatabaseStream:
    return new DatabaseStream();

  case OPS_STREAM_TAGS_DummyStream:
    return new DummyStream();

  default:
    opserr << "XaraClassBroker::getPtrNewStream - ";
    opserr << " - no DataOutputHandler type exists for class tag ";
    opserr << classTag << "\n";
    return 0;
  }
}

Recorder *
XaraClassBroker::getPtrNewRecorder(int classTag)
{
  switch (classTag) {
  case RECORDER_TAGS_ElementRecorder:
    return new ElementRecorder();

  case RECORDER_TAGS_NodeRecorder:
    return new NodeRecorder();

  case RECORDER_TAGS_EnvelopeNodeRecorder:
    return new EnvelopeNodeRecorder();

  case RECORDER_TAGS_EnvelopeElementRecorder:
    return new EnvelopeElementRecorder();

  case RECORDER_TAGS_VTK_Recorder:
    return new VTK_Recorder();

  case RECORDER_TAGS_DriftRecorder:
    return new DriftRecorder();

  case RECORDER_TAGS_TclFeViewer:
    return 0;

  case RECORDER_TAGS_GmshRecorder:
    return new GmshRecorder();

    //        case RECORDER_TAGS_MPCORecorder:
    //          return new MPCORecorder();

  default:
    opserr << "XaraClassBroker::getNewRecordr - ";
    opserr << " - no Recorder type exists for class tag ";
    opserr << classTag << "\n";
    return 0;
  }
}

/*****************************************
 *
 * METHODS TO GET NEW ANALYSIS CLASSES
 *
 *****************************************/

ConstraintHandler *
XaraClassBroker::getNewConstraintHandler(int classTag)
{
  return nullptr;
}

DOF_Numberer *
XaraClassBroker::getNewNumberer(int classTag)
{
  return nullptr;
}

AnalysisModel *
XaraClassBroker::getNewAnalysisModel(int classTag)
{
  return nullptr;
}

EquiSolnAlgo *
XaraClassBroker::getNewEquiSolnAlgo(int classTag)
{
  return nullptr;
}

Accelerator *
XaraClassBroker::getAccelerator(int classTag)
{
  return nullptr;
}

LineSearch *
XaraClassBroker::getLineSearch(int classTag)
{
  return nullptr;
}

DomainDecompAlgo *
XaraClassBroker::getNewDomainDecompAlgo(int classTag)
{
  return nullptr;
}

StaticIntegrator *
XaraClassBroker::getNewStaticIntegrator(int classTag)
{
  switch (classTag) {
  case INTEGRATOR_TAGS_LoadControl:
    return new LoadControl(1.0, 1, 1.0, .10);

  default:
    opserr << "XaraClassBroker::getNewStaticIntegrator - ";
    opserr << " - no StaticIntegrator type exists for class tag ";
    opserr << classTag << "\n";
    return nullptr;
  }
}


TransientIntegrator *
XaraClassBroker::getNewTransientIntegrator(int classTag)
{
  return nullptr;
}

IncrementalIntegrator *
XaraClassBroker::getNewIncrementalIntegrator(int classTag)
{
  return nullptr;
}


LinearSOE *
XaraClassBroker::getNewLinearSOE(int classTagSOE)
{
  return nullptr;
}

EigenSOE *
XaraClassBroker::getNewEigenSOE(int classTagSOE)
{
  return nullptr;
}

DomainSolver *
XaraClassBroker::getNewDomainSolver()
{
  return lastDomainSolver;
}


LinearSOE *
XaraClassBroker::getPtrNewDDLinearSOE(int classTagSOE,
                                            int classTagDDSolver)
{
  ProfileSPDLinSubstrSolver *theProfileSPDSolver = 0;

  switch (classTagSOE) {
  case LinSOE_TAGS_ProfileSPDLinSOE:

    if (classTagDDSolver == SOLVER_TAGS_ProfileSPDLinSubstrSolver) {
      theProfileSPDSolver = new ProfileSPDLinSubstrSolver();
      LinearSOE *theSOE = new ProfileSPDLinSOE(*theProfileSPDSolver);
      lastDomainSolver = theProfileSPDSolver;
      return theSOE;
    } else {
      opserr << "XaraClassBroker::getNewLinearSOE - ";
      opserr << " - no ProfileSPD Domain Solver type exists for class tag ";
      opserr << classTagDDSolver << "\n";
      return 0;
    }

  default:
    opserr << "XaraClassBroker::getNewLinearSOE - ";
    opserr << " - no LinearSOE type exists for class tag ";
    opserr << classTagSOE << "\n";
    return 0;
  }
}


DomainDecompositionAnalysis *
XaraClassBroker::getNewDomainDecompAnalysis(int classTag,
                                                  [[maybe_unused]] Subdomain &theSubdomain)
{
  switch (classTag) {
//case DomDecompANALYSIS_TAGS_DomainDecompositionAnalysis:
//  return new DomainDecompositionAnalysis(theSubdomain);

#ifdef _PARALLEL_PROCESSING
  case ANALYSIS_TAGS_StaticDomainDecompositionAnalysis:
    return new StaticDomainDecompositionAnalysis(theSubdomain);

  case ANALYSIS_TAGS_TransientDomainDecompositionAnalysis:
    return new TransientDomainDecompositionAnalysis(theSubdomain);
#endif

  default:
    opserr << "XaraClassBroker::getNewDomainDecompAnalysis ";
    opserr << " - no DomainDecompAnalysis type exists for class tag ";
    opserr << classTag << "\n";
    return 0;
  }
}

Subdomain *
XaraClassBroker::getSubdomainPtr(int classTag)
{
  opserr << "XaraClassBroker: NOT IMPLEMENTED YET";
  return 0;
}

int
XaraClassBroker::addUniaxialMaterial(int classTag, const char *lib,
                                           const char *funcName,
                                           UniaxialMaterial *(*funcPtr)(void))
{
  // check to see if it's already added

  UniaxialPackage *matCommands = theUniaxialPackage;
  bool found = false;
  while (matCommands != NULL && found == false) {
    if ((strcmp(lib, matCommands->libName) == 0) &&
        (strcmp(funcName, matCommands->funcName) == 0)) {
      return 0;
    }
  }

  //
  // if funPtr == 0; go get the handle
  //

  void *libHandle;
  if (funcPtr == 0) {
    if (getLibraryFunction(lib, funcName, &libHandle, (void **)&funcPtr) != 0) {
      opserr << "XaraClassBroker::addUniaxialMaterial - could not find "
                "function\n";
      return -1;
    }
  }

  //
  // add the new funcPtr
  //

  char *libNameCopy = new char[strlen(lib) + 1];
  char *funcNameCopy = new char[strlen(funcName) + 1];
  UniaxialPackage *theMat = new UniaxialPackage;
  if (libNameCopy == 0 || funcNameCopy == 0 || theMat == 0) {
    opserr << "XaraClassBroker::addUniaxialMaterial - could not add lib, "
              "out of memory\n";
    return -1;
  }
  strcpy(libNameCopy, lib);
  strcpy(funcNameCopy, funcName);

  theMat->classTag = classTag;
  theMat->funcName = funcNameCopy;
  theMat->libName = libNameCopy;
  theMat->funcPtr = funcPtr;
  theMat->next = theUniaxialPackage;
  theUniaxialPackage = theMat;

  return 0;
}

Parameter *
XaraClassBroker::getParameter(int classTag)
{
  Parameter *theRes = 0;

  switch (classTag) {
  case PARAMETER_TAG_Parameter:
    theRes = new Parameter;
    break;

  case PARAMETER_TAG_ElementParameter:
    theRes = new ElementParameter;
    break;

  case PARAMETER_TAG_MaterialStageParameter:
    theRes = new MaterialStageParameter();
    break;

  case PARAMETER_TAG_MatParameter:
    theRes = new MatParameter();
    break;

  case PARAMETER_TAG_InitialStateParameter:
    theRes = new InitialStateParameter();
    break;

  case PARAMETER_TAG_ElementStateParameter:
    theRes = new ElementStateParameter();
    break;

  default:;
  }

  return theRes;
}
