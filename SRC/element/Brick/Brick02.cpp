
//
// Eight node Brick element
//
#include <stdlib.h>
#include <cmath>

#include <ID.h>
#include <Vector.h>
#include <Matrix.h>
#include <Element.h>
#include <Node.h>
#include <Domain.h>
#include <Brick02.h>
#include <shp3d.h>
#include <ElementResponse.h>
#include <Parameter.h>
#include <ElementalLoad.h>
#include <isoparametric.tpp>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <analysis/fe_ele/TemplateElementFE.h>
#include <analysis/fe_ele/ElementFE.h>


using namespace OpenSees;

// static data
Matrix Brick02::mass(24,24);



Brick02::Brick02()
: Element(0, ELE_TAG_Brick),
  conn(NEN), applyLoad(0), load(0), Ki(0),
  K_wrap(0,0), p_wrap(0), response_wrap(0), inertia_wrap(0)
{
  stiff.zero();
  resid.zero();
  response.zero();
  inertia.zero();

  K_wrap.setData(stiff);
  p_wrap.setData(resid);
  response_wrap.setData(response);
  inertia_wrap.setData(inertia);

  for (int i=0; i<NEN; i++) {
    materialPointers[i] = nullptr;
    theNodes[i] = nullptr;
  }

  b[0] = 0.0;
  b[1] = 0.0;
  b[2] = 0.0;

  appliedB[0] = 0.0;
  appliedB[1] = 0.0;
  appliedB[2] = 0.0;
}


Brick02::Brick02(int tag,
                 const std::array<int, 8>& nodes,
                 NDMaterial &theMaterial,
                 double b1, double b2, double b3)
  : Element(tag, ELE_TAG_Brick),
    conn(NEN), applyLoad(0), load(0), Ki(0),
    K_wrap(0,0), p_wrap(0), response_wrap(0), inertia_wrap(0)
{
  stiff.zero();
  resid.zero();
  response.zero();
  inertia.zero();

  K_wrap.setData(stiff);
  p_wrap.setData(resid);
  response_wrap.setData(response);
  inertia_wrap.setData(inertia);

  for (int i=0; i<NEN; i++) {
    conn(i) = nodes[i];
    theNodes[i] = nullptr;
  }

  for (int i=0; i<NIP; i++)
    materialPointers[i] = theMaterial.getCopy("ThreeDimensional");

  // Body forces
  b[0] = b1;
  b[1] = b2;
  b[2] = b3;

  appliedB[0] = 0.0;
  appliedB[1] = 0.0;
  appliedB[2] = 0.0;
}


Brick02::~Brick02()
{
  for (int i=0; i<NIP; i++)
    delete materialPointers[i];

  if (load != nullptr)
    delete load;

  if (Ki != nullptr)
    delete Ki;
}


FE_Element*
Brick02::createFE_Element(int tag)
{
  // if (getenv("OLD_FE") != nullptr)
  //   return new ElementFE(tag, this);
  // else
    return new TemplateElementFE<NDOF>(tag, *this);
}


void
Brick02::setDomain(Domain *theDomain)
{
  for (int i=0; i<NEN; i++)
    theNodes[i] = theDomain->getNode(conn(i));

  if (theDomain != nullptr)
    this->Element::link(*theDomain);

  computeBasis();
}


int
Brick02::getNumExternalNodes() const
{
  return NEN;
}


const ID&
Brick02::getExternalNodes()
{
  return conn;
}


Node **
Brick02::getNodePtrs()
{
  return &theNodes[0];
}


int
Brick02::getNumDOF()
{
  return NDOF;
}


int
Brick02::commitState()
{
  int success = 0;

  if ((success = this->Element::commitState()) != 0)
    opserr << "Brick02::commitState () - failed in base class";

  for (int i=0; i<NIP; i++)
    success += materialPointers[i]->commitState();

  return success;
}


int
Brick02::revertToLastCommit()
{
  int success = 0;

  for (int i=0; i<NIP; i++)
    success += materialPointers[i]->revertToLastCommit();

  return success;
}


int
Brick02::revertToStart()
{
  int success = 0;

  for (int i=0; i<NIP; i++)
    success += materialPointers[i]->revertToStart();

  return success;
}


const Matrix&
Brick02::getTangentStiff()
{
  return K_wrap;
}


const Matrix&
Brick02::getInitialStiff()
{
  if (Ki != nullptr)
    return *Ki;

  MatrixND<NDOF,NDOF> initial{};

  // compute basis vectors and local nodal coordinates
  computeBasis();

  int count = 0;
  for (int i=0; i<2; i++) {
    for (int j=0; j<2; j++) {
      for (int k=0; k<2; k++) {
        const double gaussPoint[NDM] = {sg[i], sg[j], sg[k]};
        double xsj;
        double shp[4][NEN];
        shp3d(gaussPoint, xsj, shp, xl);

        const int ip = count++;
        const double dvol = wg[ip] * xsj;
        const Matrix& D = materialPointers[ip]->getInitialTangent();

        MatrixND<NST,NDF> B[NEN];
        MatrixND<NDF,NST> BTD[NEN];
        for (int node=0; node<NEN; node++) {
          computeB(node, shp, B[node]);
          BTD[node].addMatrixTransposeProduct(0.0, B[node], D, dvol);
        }

        for (int jnode=0; jnode<NEN; jnode++) {
          for (int knode=0; knode<NEN; knode++) {
            MatrixND<NDF,NDF> stiffJK;
            stiffJK.addMatrixProduct(0.0, BTD[jnode], B[knode], 1.0);
            initial.assemble(stiffJK, NDF*jnode, NDF*knode, 1.0);
          }
        }
      }
    }
  }

  Ki = new Matrix(NDOF,NDOF);
  for (int j=0; j<NDOF; j++)
    for (int i=0; i<NDOF; i++)
      (*Ki)(i,j) = initial(i,j);

  return *Ki;
}


const Matrix&
Brick02::getMass()
{
  inertia.zero();
  formInertiaTerms(1, inertia_wrap);

  return mass;
}


// form residual and tangent
int
Brick02::update()
{
  // strains ordered : eps11, eps22, eps33, 2*eps12, 2*eps23, 2*eps31

  stiff.zero();
  resid.zero();

  int count = 0;
  for (int i=0; i<2; i++) {
    for (int j=0; j<2; j++) {
      for (int k=0; k<2; k++) {
        const double gaussPoint[NDM] = {sg[i], sg[j], sg[k]};
        double xsj;
        double shp[4][NEN];
        shp3d(gaussPoint, xsj, shp, xl);

        const double dvol = wg[count] * xsj;

        VectorND<NST> strain{};
        MatrixND<NST,NDF> B[NEN];

        // node loop to compute strain
        for (int node=0; node<NEN; node++) {
          const double b00 = shp[0][node];
          const double b11 = shp[1][node];
          const double b22 = shp[2][node];
          const double b30 = shp[1][node];
          const double b31 = shp[0][node];
          const double b41 = shp[2][node];
          const double b42 = shp[1][node];
          const double b50 = shp[2][node];
          const double b52 = shp[0][node];

          const Vector &ul = theNodes[node]->getTrialDisp();

          const double ul0 = ul(0);
          const double ul1 = ul(1);
          const double ul2 = ul(2);

          strain(0) += b00 * ul0;
          strain(1) += b11 * ul1;
          strain(2) += b22 * ul2;
          strain(3) += b30 * ul0 + b31 * ul1;
          strain(4) += b41 * ul1 + b42 * ul2;
          strain(5) += b50 * ul0 + b52 * ul2;

          computeB(node, shp, B[node]);
        }

        materialPointers[count]->setTrialStrain(strain);

        const Vector& stress = materialPointers[count]->getStress();
        const Matrix& D = materialPointers[count]->getTangent();

        const double stress0 = stress(0)*dvol;
        const double stress1 = stress(1)*dvol;
        const double stress2 = stress(2)*dvol;
        const double stress3 = stress(3)*dvol;
        const double stress4 = stress(4)*dvol;
        const double stress5 = stress(5)*dvol;

        MatrixND<NDF,NST> BTD[NEN];
        for (int node=0; node<NEN; node++)
          BTD[node].addMatrixTransposeProduct(0.0, B[node], D, dvol);

        for (int jnode=0; jnode<NEN; jnode++) {
          const double b00 = shp[0][jnode];
          const double b11 = shp[1][jnode];
          const double b22 = shp[2][jnode];
          const double b30 = shp[1][jnode];
          const double b31 = shp[0][jnode];
          const double b41 = shp[2][jnode];
          const double b42 = shp[1][jnode];
          const double b50 = shp[2][jnode];
          const double b52 = shp[0][jnode];

          const VectorND<NDF> residJ {
            b00*stress0 + b30*stress3 + b50*stress5,
            b11*stress1 + b31*stress3 + b41*stress4,
            b22*stress2 + b42*stress4 + b52*stress5
          };

          const int jj = NDF*jnode;
          resid.assemble(jj, residJ, 1.0);

          for (int p=0; p<NDF; p++) {
            if (applyLoad == 0)
              resid(jj+p) -= dvol*b[p]*shp[3][jnode];
            else
              resid(jj+p) -= dvol*appliedB[p]*shp[3][jnode];
          }

          for (int knode=0; knode<NEN; knode++) {
            MatrixND<NDF,NDF> stiffJK;
            stiffJK.addMatrixProduct(0.0, BTD[jnode], B[knode], 1.0);
            stiff.assemble(stiffJK, jj, NDF*knode, 1.0);
          }
        }

        count++;
      }
    }
  }

  return 0;
}


void
Brick02::formInertiaTerms(int tangFlag, Vector &resid)
{
  mass.Zero();

  int count = 0;
  for (int i=0; i<2; i++) {
    for (int j=0; j<2; j++) {
      for (int k=0; k<2; k++) {
        const double gaussPoint[NDM] = {sg[i], sg[j], sg[k]};
        double xsj;
        double shp[4][NEN];
        shp3d(gaussPoint, xsj, shp, xl);

        const int ip = count++;
        const double dvol = wg[ip] * xsj;

        VectorND<NDF> momentum{};
        for (int node=0; node<NEN; node++)
          momentum.addVector(1.0, theNodes[node]->getTrialAccel(), shp[3][node]);

        const double rho = materialPointers[ip]->getRho();
        momentum *= rho;

        for (int jnode=0; jnode<NEN; jnode++) {
          const int jj = NDF*jnode;
          double temp = shp[3][jnode] * dvol;

          for (int p=0; p<NDF; p++)
            resid(jj+p) += temp * momentum(p);

          if (tangFlag == 1) {
            temp *= rho;

            for (int knode=0; knode<NEN; knode++) {
              const double massJK = temp * shp[3][knode];
              const int kk = NDF*knode;

              for (int p=0; p<NDF; p++)
                mass(jj+p, kk+p) += massJK;
            }
          }
        }
      }
    }
  }
}


void
Brick02::zeroLoad()
{
  if (load != nullptr)
    load->Zero();

  applyLoad = 0;

  appliedB[0] = 0.0;
  appliedB[1] = 0.0;
  appliedB[2] = 0.0;

  return;
}


int
Brick02::addLoad(ElementalLoad *theLoad, double loadFactor)
{
  int type;
  const Vector &data = theLoad->getData(type, loadFactor);

  if (type == LOAD_TAG_BrickSelfWeight) {
    applyLoad = 1;
    appliedB[0] += loadFactor * b[0];
    appliedB[1] += loadFactor * b[1];
    appliedB[2] += loadFactor * b[2];
    return 0;
  }
  else if (type == LOAD_TAG_SelfWeight) {
    // added compatibility with selfWeight class implemented for all continuum elements, C.McGann, U.W.
    applyLoad = 1;
    appliedB[0] += loadFactor*data(0)*b[0];
    appliedB[1] += loadFactor*data(1)*b[1];
    appliedB[2] += loadFactor*data(2)*b[2];
    return 0;
  }
  else {
    opserr << "Brick02::addLoad() - ele with tag: " << this->getTag()
           << " does not deal with load type: " << type << "\n";
    return -1;
  }
}


int
Brick02::addInertiaLoadToUnbalance(const Vector &accel)
{
  int haveRho = 0;
  for (int i=0; i<NIP; i++) {
    if (materialPointers[i]->getRho() != 0.0)
      haveRho = 1;
  }

  if (haveRho == 0)
    return 0;

  inertia.zero();
  formInertiaTerms(1, inertia_wrap);

  // store computed RV for nodes in inertia vector
  int count = 0;
  for (int i=0; i<NEN; i++) {
    const Vector &Raccel = theNodes[i]->getRV(accel);
    for (int j=0; j<NDF; j++)
      inertia(count++) = Raccel(j);
  }

  // create the load vector if one does not exist
  if (load == nullptr)
    load = new Vector(NDOF);

  // add -M * RV(accel) to the load vector
  load->addMatrixVector(1.0, mass, inertia_wrap, -1.0);

  return 0;
}


const Vector&
Brick02::getResistingForce()
{
  if (load == nullptr)
    return p_wrap;

  response = resid;
  response_wrap -= *load;
  return response_wrap;
}


const Vector&
Brick02::getResistingForceIncInertia()
{
  response = resid;
  formInertiaTerms(0, response_wrap);

  // add the damping forces if rayleigh damping
  if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
    response_wrap += this->getRayleighDampingForces();

  if (load != nullptr)
    response_wrap -= *load;

  return response_wrap;
}


// compute local coordinates and basis

void
Brick02::computeBasis()
{
  // nodal coordinates
  for (int i=0; i<NEN; i++) {
    const Vector &coorI = theNodes[i]->getCrds();

    xl[0][i] = coorI(0);
    xl[1][i] = coorI(1);
    xl[2][i] = coorI(2);
  }
}


const MatrixND<6,3>&
Brick02::computeB(int node,
                  const double shp[4][8],
                  MatrixND<6,3> &B
                ) const noexcept
{

//---B Matrix in standard {1,2,3} mechanics notation---------
//
//                -                   -
//               | N,1      0     0    |
//   B       =   |   0     N,2    0    |
//               |   0      0     N,3  |   (6x3)
//               | N,2     N,1     0   |
//               |   0     N,3    N,2  |
//               | N,3      0     N,1  |
//                -                   -
//
//-------------------------------------------------------------------
  B.zero();
  B(0,0) = shp[0][node];
  B(1,1) = shp[1][node];
  B(2,2) = shp[2][node];

  B(3,0) = shp[1][node];
  B(3,1) = shp[0][node];

  B(4,1) = shp[2][node];
  B(4,2) = shp[1][node];

  B(5,0) = shp[2][node];
  B(5,2) = shp[0][node];

  return B;
}



Response*
Brick02::setResponse(const char **argv, int argc, OPS_Stream &output)
{
  Response *theResponse = nullptr;

  char outputData[32];

  output.tag("ElementOutput");
  output.attr("eleType","Brick02");
  output.attr("eleTag",this->getTag());
  for (int i=1; i<=8; i++) {
    sprintf(outputData,"node%d",i);
    output.attr(outputData, theNodes[i-1]->getTag());
  }


  if (strcmp(argv[0],"force") == 0 || strcmp(argv[0],"forces") == 0) {

    for (int i=1; i<=8; i++) {
      sprintf(outputData,"P1_%d",i);
      output.tag("ResponseType",outputData);
      sprintf(outputData,"P2_%d",i);
      output.tag("ResponseType",outputData);
      sprintf(outputData,"P3_%d",i);
      output.tag("ResponseType",outputData);
    }

    theResponse = new ElementResponse(this, 1, Vector(NDOF));
  }
  else if (strcmp(argv[0],"material") == 0 || strcmp(argv[0],"integrPoint") == 0) {

    int pointNum = atoi(argv[1]);
    if (pointNum > 0 && pointNum <= NIP) {

      output.tag("GaussPoint");
      output.attr("number",pointNum);

      theResponse = materialPointers[pointNum-1]->setResponse(&argv[2], argc-2, output);

      output.endTag(); // GaussPoint
    }
  }

  else if ((strcmp(argv[0],"stresses") ==0) || (strcmp(argv[0],"stress") ==0)) {

    for (int i=0; i<NIP; i++) {
      output.tag("GaussPoint");
      output.attr("number",i+1);
      output.tag("NdMaterialOutput");
      output.attr("classType", materialPointers[i]->getClassTag());
      output.attr("tag",       materialPointers[i]->getTag());

      output.tag("ResponseType","sigma11");
      output.tag("ResponseType","sigma22");
      output.tag("ResponseType","sigma33");
      output.tag("ResponseType","sigma12");
      output.tag("ResponseType","sigma23");
      output.tag("ResponseType","sigma13");

      output.endTag(); // NdMaterialOutput
      output.endTag(); // GaussPoint
    }
    theResponse = new ElementResponse(this, 3, Vector(NST*NIP));

  }

  else if ((strcmp(argv[0],"strains") ==0) || (strcmp(argv[0],"strain") ==0)) {

    for (int i=0; i<NIP; i++) {
      output.tag("GaussPoint");
      output.attr("number",i+1);
      output.tag("NdMaterialOutput");
      output.attr("classType", materialPointers[i]->getClassTag());
      output.attr("tag", materialPointers[i]->getTag());

      output.tag("ResponseType","eps11");
      output.tag("ResponseType","eps22");
      output.tag("ResponseType","eps33");
      output.tag("ResponseType","eps12");
      output.tag("ResponseType","eps23");
      output.tag("ResponseType","eps13");

      output.endTag(); // NdMaterialOutput
      output.endTag(); // GaussPoint
    }
    theResponse = new ElementResponse(this, 4, Vector(NST*NIP));
  }

  else if (strcmp(argv[0], "stressAtNodes") == 0) {
    theResponse = new ElementResponse(this, 11, Vector(NST*NEN));
  }

  else if (strcmp(argv[0], "shape") == 0) {
    output.tag("Shape");
    output.attr("number", 1);
    output.attr("type", "Brick02");
    output.attr("tag", this->getTag());
    theResponse = new ElementResponse(this, 500 + atoi(argv[1])*10 + atoi(argv[2]), Vector(NEN));
  }

  output.endTag(); // ElementOutput
  return theResponse;
}


int
Brick02::getResponse(int responseID, Information &eleInfo)
{
  static Vector stresses(NST*NIP);

  if (responseID == 1)
    return eleInfo.setVector(this->getResistingForce());

  else if (responseID == 2)
    return eleInfo.setMatrix(this->getTangentStiff());

  else if (responseID == 3) {

    // Loop over the integration points
    int cnt = 0;
    for (int i=0; i<NIP; i++) {

      // Get material stress response
      const Vector &sigma = materialPointers[i]->getStress();
      stresses(cnt++) = sigma(0);
      stresses(cnt++) = sigma(1);
      stresses(cnt++) = sigma(2);
      stresses(cnt++) = sigma(3);
      stresses(cnt++) = sigma(4);
      stresses(cnt++) = sigma(5);
    }
    return eleInfo.setVector(stresses);

  }
  else if (responseID == 4) {

    // Loop over the integration points
    int cnt = 0;
    for (int i=0; i<NIP; i++) {

      // Get material strain response
      const Vector &sigma = materialPointers[i]->getStrain();
      stresses(cnt++) = sigma(0);
      stresses(cnt++) = sigma(1);
      stresses(cnt++) = sigma(2);
      stresses(cnt++) = sigma(3);
      stresses(cnt++) = sigma(4);
      stresses(cnt++) = sigma(5);
    }
    return eleInfo.setVector(stresses);
  }

  else if (responseID == 11) {
    constexpr OpenSees::MatrixND<8,8> We {{
        2.549038105676660, -0.683012701892219,  0.183012701892219, -0.683012701892219, -0.683012701892219, 0.183012701892219, -0.049038105676658, 0.183012701892219,
       -0.683012701892219,  2.549038105676660, -0.683012701892219,  0.183012701892219, 0.183012701892219, -0.683012701892219, 0.183012701892219, -0.049038105676658,
        0.183012701892219, -0.683012701892219,  2.549038105676660, -0.683012701892219, -0.049038105676658, 0.183012701892219, -0.683012701892219, 0.183012701892219,
       -0.683012701892219,  0.183012701892219, -0.683012701892219,  2.549038105676660, 0.183012701892219, -0.049038105676658, 0.183012701892219, -0.683012701892219,
       -0.683012701892219,  0.183012701892219, -0.049038105676658,  0.183012701892219, 2.549038105676660, -0.683012701892219, 0.183012701892219, -0.683012701892219,
        0.183012701892219, -0.683012701892219,  0.183012701892219, -0.049038105676658, -0.683012701892219, 2.549038105676660, -0.683012701892219, 0.183012701892219,
       -0.049038105676658,  0.183012701892219, -0.683012701892219,  0.183012701892219, 0.183012701892219, -0.683012701892219, 2.549038105676660, -0.683012701892219,
        0.183012701892219, -0.049038105676658,  0.183012701892219, -0.683012701892219, -0.683012701892219, 0.183012701892219, -0.683012701892219, 2.549038105676660
    }};

    static VectorND<NST*NEN> stressAtNodes;
    static Vector output(stressAtNodes);

    stressAtNodes.zero();
    Xara::StressExtrapolation<NEN,NIP,NST>(materialPointers, We, stressAtNodes);
    return eleInfo.setVector(output);
  }

  else
    return -1;
}


int
Brick02::setParameter(const char **argv, int argc, Parameter &param)
{
  if (argc < 1)
    return -1;

  int res = -1;

  if ((strstr(argv[0],"material") != 0) && (strcmp(argv[0],"materialState") != 0)) {

    if (argc < 3)
      return -1;

    int pointNum = atoi(argv[1]);
    if (pointNum > 0 && pointNum <= NIP)
      return materialPointers[pointNum-1]->setParameter(&argv[2], argc-2, param);
    else
      return -1;
  }

  // otherwise it could be just a forall material parameter
  else {
    int matRes = res;
    for (int i=0; i<NIP; i++) {
      matRes = materialPointers[i]->setParameter(argv, argc, param);
      if (matRes != -1)
        res = matRes;
    }
  }

  return res;
}


int
Brick02::updateParameter(int parameterID, Information &info)
{
  int res = -1;
  int matRes = res;

  if (parameterID == res) {
    return -1;
  }
  else {
    for (int i=0; i<NIP; i++)
      matRes = materialPointers[i]->updateParameter(parameterID, info);

    if (matRes != -1)
      res = matRes;

    return res;
  }
}


int
Brick02::sendSelf(int commitTag, Channel &theChannel)
{
  int res = 0;

  // note: we don't check for dataTag == 0 for Element
  // objects as that is taken care of in a commit by the Domain
  // object - don't want to have to do the check if sending data
  int dataTag = this->getDbTag();

  // Brick packs its data into an ID and a Vector and sends them to theChannel
  // along with its dbTag and the commitTag passed in the arguments

  // Now Brick sends the ids of its materials
  int matDbTag;

  static ID idData(26);

  idData(24) = this->getTag();
  if (alphaM != 0 || betaK != 0 || betaK0 != 0 || betaKc != 0)
    idData(25) = 1;
  else
    idData(25) = 0;


  for (int i=0; i<NIP; i++) {
    idData(i) = materialPointers[i]->getClassTag();
    matDbTag = materialPointers[i]->getDbTag();
    // NOTE: we do have to ensure that the material has a database
    // tag if we are sending to a database channel.
    if (matDbTag == 0) {
      matDbTag = theChannel.getDbTag();
      if (matDbTag != 0)
        materialPointers[i]->setDbTag(matDbTag);
    }
    idData(i+8) = matDbTag;
  }

  idData(16) = conn(0);
  idData(17) = conn(1);
  idData(18) = conn(2);
  idData(19) = conn(3);
  idData(20) = conn(4);
  idData(21) = conn(5);
  idData(22) = conn(6);
  idData(23) = conn(7);

  res += theChannel.sendID(dataTag, commitTag, idData);
  if (res < 0) {
    opserr << "WARNING Brick02::sendSelf() - " << this->getTag() << " failed to send ID\n";
    return res;
  }

  static Vector dData(7);
  dData(0) = alphaM;
  dData(1) = betaK;
  dData(2) = betaK0;
  dData(3) = betaKc;
  dData(4) = b[0];
  dData(5) = b[1];
  dData(6) = b[2];

  if (theChannel.sendVector(dataTag, commitTag, dData) < 0) {
    opserr << "Brick02::sendSelf() - failed to send double data\n";
    return -1;
  }

  // Finally, Brick asks its material objects to send themselves
  for (int i=0; i<NIP; i++) {
    res += materialPointers[i]->sendSelf(commitTag, theChannel);
    if (res < 0) {
      opserr << "WARNING Brick02::sendSelf() - " << this->getTag() << " failed to send its Material\n";
      return res;
    }
  }

  return res;
}


int
Brick02::recvSelf(int commitTag,
                  Channel &theChannel,
                  FEM_ObjectBroker &theBroker)
{
  int res = 0;

  int dataTag = this->getDbTag();

  static ID idData(26);
  res += theChannel.recvID(dataTag, commitTag, idData);
  if (res < 0) {
    opserr << "WARNING Brick02::recvSelf() - " << this->getTag() << " failed to receive ID\n";
    return res;
  }

  this->setTag(idData(24));

  static Vector dData(7);
  if (theChannel.recvVector(dataTag, commitTag, dData) < 0) {
    opserr << "Brick02::recvSelf() - failed to recv double data\n";
    return -1;
  }
  alphaM = dData(0);
  betaK = dData(1);
  betaK0 = dData(2);
  betaKc = dData(3);
  b[0] = dData(4);
  b[1] = dData(5);
  b[2] = dData(6);


  conn(0) = idData(16);
  conn(1) = idData(17);
  conn(2) = idData(18);
  conn(3) = idData(19);
  conn(4) = idData(20);
  conn(5) = idData(21);
  conn(6) = idData(22);
  conn(7) = idData(23);


  if (materialPointers[0] == nullptr) {
    for (int i=0; i<NIP; i++) {
      int matClassTag = idData(i);
      int matDbTag = idData(i+8);
      // Allocate new material with the sent class tag
      materialPointers[i] = theBroker.getNewNDMaterial(matClassTag);
      if (materialPointers[i] == nullptr) {
        opserr << "Brick02::recvSelf() - Broker could not create NDMaterial of class type " << matClassTag << "\n";
        return -1;
      }
      // Now receive materials into the newly allocated space
      materialPointers[i]->setDbTag(matDbTag);
      res += materialPointers[i]->recvSelf(commitTag, theChannel, theBroker);
      if (res < 0) {
        opserr << "Brick02::recvSelf() - material " << i << " failed to recv itself\n";
        return res;
      }
    }
  }
  // materials exist, ensure materials of correct type and recvSelf on them
  else {
    for (int i=0; i<NIP; i++) {
      int matClassTag = idData(i);
      int matDbTag = idData(i+8);
      // Check that material is of the right type; if not,
      // delete it and create a new one of the right type
      if (materialPointers[i]->getClassTag() != matClassTag) {
        delete materialPointers[i];
        materialPointers[i] = theBroker.getNewNDMaterial(matClassTag);
        if (materialPointers[i] == nullptr) {
          opserr << "Broker could not create NDMaterial of class type "
                 << matClassTag << "\n";
          return -1;
        }
        materialPointers[i]->setDbTag(matDbTag);
      }
      // Receive the material

      res += materialPointers[i]->recvSelf(commitTag, theChannel, theBroker);
      if (res < 0) {
        opserr << "Brick02::recvSelf() - material " << i << " failed to recv itself\n";
        return res;
      }
    }
  }

  return res;
}


void
Brick02::Print(OPS_Stream &s, int flag)
{
  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << OPS_PRINT_JSON_ELEM_INDENT << "{";
    s << "\"name\": " << this->getTag() << ", ";
    s << "\"type\": \"Brick02\", ";
    s << "\"nodes\": ["
      << conn(0) << ", ";
    for (int i=1; i<7; i++)
      s << conn(i) << ", ";
    s << conn(7) << "], ";
    s << "\"bodyForces\": [" << b[0] << ", " << b[1] << ", " << b[2] << "], ";
    s << "\"material\": [" << materialPointers[0]->getTag() << "]}";

    return;
  }

  if (flag == 2) {

    s << "#Brick02\n";

    int i;
    const int numNodes = 8;
    const int nstress = 6;

    for (i=0; i<numNodes; i++) {
      const Vector &nodeCrd = theNodes[i]->getCrds();
      const Vector &nodeDisp = theNodes[i]->getDisp();
      s << "#NODE " << nodeCrd(0) << " " << nodeCrd(1) << " " << nodeCrd(2)
        << " " << nodeDisp(0) << " " << nodeDisp(1) << " " << nodeDisp(2) << "\n";
    }

    // spit out the section location & invoke print on the section
    const int numMaterials = 8;

    static Vector avgStress(nstress);
    static Vector avgStrain(nstress);
    avgStress.Zero();
    avgStrain.Zero();
    for (i=0; i<numMaterials; i++) {
      avgStress += materialPointers[i]->getStress();
      avgStrain += materialPointers[i]->getStrain();
    }
    avgStress /= numMaterials;
    avgStrain /= numMaterials;

    s << "#AVERAGE_STRESS ";
    for (i=0; i<nstress; i++)
      s << avgStress(i) << " ";
    s << "\n";

    s << "#AVERAGE_STRAIN ";
    for (i=0; i<nstress; i++)
      s << avgStrain(i) << " ";
    s << "\n";
  }

  if (flag == OPS_PRINT_CURRENTSTATE) {
    s << "Standard Eight Node Brick02 \n";
    s << "Element Number: " << this->getTag() << "\n";
    s << "Nodes: " << conn;

    s << "Material Information : \n ";
    materialPointers[0]->Print(s, flag);

    s << "\n";
    s << this->getTag() << " " << conn(0)
      << " " << conn(1)
      << " " << conn(2)
      << " " << conn(3)
      << " " << conn(4)
      << " " << conn(5)
      << " " << conn(6)
      << " " << conn(7)
      << "\n";

    s << "Body Forces: " << b[0] << " " << b[1] << " " << b[2] << "\n";
    s << "Resisting Force (no inertia): " << this->getResistingForce();
  }
}
