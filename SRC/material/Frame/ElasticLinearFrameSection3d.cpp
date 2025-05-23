//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation    
//
//===----------------------------------------------------------------------===//
//
#include <ElasticLinearFrameSection3d.h>
#include <Matrix.h>
#include <Vector.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Parameter.h>
#include <classTags.h>

using namespace OpenSees;

constexpr int SEC_TAG_ElasticLinearFrame3d = 0;

static int layout_array[] = {
    FrameStress::N,
    FrameStress::Vy,
    FrameStress::Vz,
    FrameStress::T,
    FrameStress::My,
    FrameStress::Mz,
    FrameStress::Bimoment,
    FrameStress::By,
    FrameStress::Bz,
    FrameStress::Bishear,
    FrameStress::Qy,
    FrameStress::Qz
};

ID ElasticLinearFrameSection3d::layout(layout_array, nr);


ElasticLinearFrameSection3d::ElasticLinearFrameSection3d()
: FrameSection(0, SEC_TAG_ElasticLinearFrame3d, 0, false),
  E(0.0),
  G(0.0),
  Ks(new MatrixND<nr,nr> {}),
  Ksen(nullptr)
{

}

ElasticLinearFrameSection3d::ElasticLinearFrameSection3d(int tag,
    double E_in,
    double G_in,
    const FrameSectionConstants& cons,
    //
    double mass_,
    bool use_mass
)
: FrameSection(tag, SEC_TAG_ElasticLinearFrame3d, mass_, use_mass),
  E(E_in),
  G(G_in),
  Ksen(nullptr),
  Ks(new MatrixND<nr,nr> {})
{

  // n-n
  double A  = cons.A;
  double Ay = cons.Ay;
  double Az = cons.Az;
  // m-m
  double Iy = cons.Iy;     //   \int z^2
  double Iz = cons.Iz;     //   \int y^2
  double Iyz = cons.Iyz;   //
  // w-w 
  double Cw = cons.Cw;     // 
  double Ca = cons.Ca;     // 
  // n-m
  double Qy = cons.Qy;     //  int y = A*zs
  double Qz = cons.Qz;     //  int z

  double Rw = cons.Rw;     // n-w
  double Ry = cons.Ry;     // n-v
  double Rz = cons.Rz;

  double Sa = cons.Sa;     // m-v
  double Sy = cons.Sy;     // m-w, int  z warp(x,y)
  double Sz = cons.Sz;     // m-w, int -y warp(x,y) 

  centroid = {Qz/A, Qy/A};

  // Polar moment of inertia
  const double I0   = Iy + Iz;

  double GSnyy = G*(Ay - A);
  double GSnzz = G*(Az - A);
  double GSy   = -GSnyy; 
  double GSz   = -GSnzz;

  // Centroidal moments of inertia
  // const double Ic33 =  Iy  - A*zc*zc;
  // const double Ic22 =  Iz  - A*yc*yc;
  // const double Icyz =  Iyz - A*yc*zc;

  // const double Qw3  = -Ic33*yc + Icyz*zc;
  // const double Qw2  =  Ic22*zc - Icyz*yc;

  *Ks = {
    //            N                       M                      W                 Q    
    //                        |                        |                 |                |
    {{  E*A,      0.,     0.,      0.,   E*Qy,  -E*Qz,    E*Rw ,  0.,  0.,    0.,     0.,     0.},
     {    0.,    G*A,     0.,   -G*Qy,     0.,     0.,     0.  ,  0.,  0.,  G*Ry,  GSnyy,     0.},
     {    0.,     0.,    G*A,    G*Qz,     0.,     0.,     0.  ,  0.,  0.,  G*Rz,     0.,  GSnzz},
    //                        |                        |          0.,  0.,      ,     0.,     0.|
     {    0.,  -G*Qy,   G*Qz,    G*I0,     0.,     0.,     0.  ,  0.,  0.,  G*Sa,     0.,     0.},
     {  E*Qy,     0.,     0.,      0.,  E*Iy , -E*Iyz,     E*Sy,  0.,  0.,    0.,     0.,     0.},
     { -E*Qz,     0.,     0.,      0., -E*Iyz,  E*Iz ,     E*Sz,  0.,  0.,    0.,     0.,     0.},
    //                        |                        |          0.,  0.,      ,     0.,     0.|
     {  E*Rw,     0.,     0.,      0.,   E*Sy,   E*Sz,    E*Cw ,  0.,  0.,    0.,     0.,     0.}, // Bimoment
     {    0.,     0.,     0.,      0.,     0.,     0.,      0. ,  0.,  0.,    0.,     0.,     0.},
     {    0.,     0.,     0.,      0.,     0.,     0.,      0. ,  0.,  0.,    0.,     0.,     0.},
     {    0.,   G*Ry,   G*Rz,    G*Sa,     0.,     0.,      0. ,  0.,  0.,  G*Ca,     0.,     0.}, // Bishear
     {    0.,  GSnyy,     0.,      0.,     0.,     0.,      0. ,  0.,  0.,    0.,    GSy,     0.},
     {    0.,     0.,  GSnzz,      0.,     0.,     0.,      0. ,  0.,  0.,    0.,     0.,    GSz}}
  };

    //                        |                        |                |
    //    0       1       2   |    3       4       5   |      6     7   |
}

void
ElasticLinearFrameSection3d::getConstants(FrameSectionConstants& consts) const
{
  const MatrixND<nr,nr>& K = *Ks;
  // TODO
  consts.A   =  K(0,0)/E;
  consts.Iy  =  K(4,4)/E;
  consts.Iyz = -K(4,5)/E;
  consts.Iz  =  K(5,5)/E;

  if (G != 0) {
    consts.Ay  =  (K(1,10) + K(1,1))/G;
    consts.Az  =  (K(2,11) + K(2,2))/G;
    consts.Ca  =  K(9,9)/G;
    consts.Cw  =  K(6,6)/E;
  } else {
    consts.Ay  =  0;
    consts.Az  =  0;
    consts.Ca  =  0;
  }
}

int
ElasticLinearFrameSection3d::getIntegral(Field field, State state, double& value) const
{
  FrameSectionConstants consts{};
  getConstants(consts);

  switch (field) {
    case Field::Unit:
      value = consts.A;
      return 0;

    case Field::UnitY:
      value = consts.Qy;
      return 0;

    case Field::UnitZ:
      value = consts.Qz;
      return 0;

    case Field::Density:
      if (this->FrameSection::getIntegral(field, state, value) != 0) 
        return -1;
      else
        return  0;
      // Density may be specified for the section

    case Field::UnitYY:
    case Field::UnitCentroidYY:
      value = consts.Iz;
      if (field == Field::UnitCentroidYY) {
        double yc = centroid[0];
        value -= consts.A*yc*yc;
      }
      return 0;

    case Field::UnitZZ:
    case Field::UnitCentroidZZ:
      value = consts.Iy;
      if (field == Field::UnitCentroidZZ) {
        double zc = centroid[1];
        value -= consts.A*zc*zc;
      }
      return 0;

    default:
      return -1;
  }
}

ElasticLinearFrameSection3d::~ElasticLinearFrameSection3d()
{
  if (Fs != nullptr)
    delete Fs;

  if (Ksen != nullptr)
    delete Ksen;

  return;
}

FrameSection*
ElasticLinearFrameSection3d::getFrameCopy(const FrameStressLayout& layout)
{
  // TODO: 
  // - take layout as argument
  // - overload 
  //   template<int n> ID::operator==(std::array<int, n>)
  //
  // OR
  // - add FrameSection::setLayout()
  //
  ElasticLinearFrameSection3d *theCopy = new ElasticLinearFrameSection3d();

  // Copy over all data
  *theCopy = *this;

  // Revoke any pointers that are owned by this instance
  theCopy->Ksen = nullptr;
  // return theCopy;

  int ni=0;
  bool ind[nr]{};
  double data[nr][nr]{};


  bool uniform_twist = true;
  Matrix Kc(*(theCopy->Ks));


  // TODO: check layout here
  if (uniform_twist) {
    for (int i=0; i<nr; i++)
      Kc(i,3) += Kc(i, 7);
    for (int i=0; i<nr; i++)
      Kc(3,i) += Kc(7, i);
  }

  // Count number of independent variables
  for (int i=0; i<nr; i++)
    if (Kc(i,i) != 0.0) {
      // // dont include twist term if its being condensed
      if (i == 7 && uniform_twist)
        continue;
      ind[i] = true;
      ni ++;
    }

  // Form Ki with only independent variables
  Matrix Ki(&data[0][0], ni, ni);
  int ii=0;
  for (int i=0; i<nr; i++)
    if (ind[i]) {
      int jj=0;
      for (int j=0; j<nr; j++)
        if (ind[j])
          Ki(ii,jj++) = Kc(i,j);

      ii++;
    }


  Ki.Invert();


  theCopy->Fs = new Matrix(nr,nr);
  Matrix& Fc = *(theCopy->Fs);
  Fc.Zero();
  ii = 0;
  for (int i=0; i<nr; i++)
    if (ind[i]) {
      int jj=0;
      for (int j=0; j<nr; j++)
        if (ind[j])
          Fc(i,j) = Ki(ii,jj++);

      ii++;
    }

  return theCopy;
}

FrameSection*
ElasticLinearFrameSection3d::getFrameCopy()
{
  // TODO: 
  // - take layout as argument
  // - overload 
  //   template<int n> ID::operator==(std::array<int, n>)
  //
  // OR
  // - add FrameSection::setLayout()
  //
  ElasticLinearFrameSection3d *theCopy = new ElasticLinearFrameSection3d();

  // Copy over all data
  *theCopy = *this;

  // Revoke any pointers that are owned by this instance
  theCopy->Ksen = nullptr;

  return theCopy;
}


int
ElasticLinearFrameSection3d::commitState()
{
  return 0;
}


int
ElasticLinearFrameSection3d::revertToLastCommit()
{
  return 0;
}


int
ElasticLinearFrameSection3d::revertToStart()
{
  return 0;
}


int
ElasticLinearFrameSection3d::setTrialSectionDeformation(const Vector &def)
{
  assert(def.Size() == nr);
  e = def;
  return 0;
}


const Vector &
ElasticLinearFrameSection3d::getSectionDeformation() // TODO: needed ?
{
  v.setData(e);
  return v;
}


const Vector &
ElasticLinearFrameSection3d::getStressResultant()
{
  #ifdef THREAD_SAFE 
  #else 
  static VectorND<nr> s;
  static Vector s_wrap(s);
  #endif
  s = (*Ks)*e;
  return s_wrap;
}


const Matrix &
ElasticLinearFrameSection3d::getSectionTangent()
{
  return getInitialTangent();
}


const Matrix &
ElasticLinearFrameSection3d::getInitialTangent()
{
  M.setData(*Ks);
  return M;
}


const Matrix &
ElasticLinearFrameSection3d::getSectionFlexibility()
{
    return getInitialFlexibility();
}


const Matrix &
ElasticLinearFrameSection3d::getInitialFlexibility()
{
  if (Fs == nullptr) {
    // This only happens when getCopy is called without
    // the layout  argument
    Fs = new Matrix(nr,nr);
    Matrix Kwrap(*Ks);
    Kwrap.Invert(*Fs);
  }

  return *Fs;
}


const ID&
ElasticLinearFrameSection3d::getType()
{
  return layout;
}


int
ElasticLinearFrameSection3d::getOrder() const
{
  return nr;
}

int
ElasticLinearFrameSection3d::sendSelf(int commitTag, Channel &theChannel)
{
    int res = 0;

    static Vector data(7);

    int dataTag = this->getDbTag();

    data(0) = this->getTag();
    data(1) = E;
//  data(2) = A;
//  data(3) = Iz;
//  data(4) = Iy;
    data(5) = G;
//  data(6) = J;

    res += theChannel.sendVector(dataTag, commitTag, data);
    if(res < 0) {
      opserr << "ElasticLinearFrameSection3d::sendSelf -- failed to send data\n";
      return res;
    }

    return res;
}


int
ElasticLinearFrameSection3d::recvSelf(int commitTag, Channel &theChannel,
                                          FEM_ObjectBroker &theBroker)
{
    int res = 0;

      static Vector data(7);

    int dataTag = this->getDbTag();

    res += theChannel.recvVector(dataTag, commitTag, data);
    if(res < 0) {
      opserr << "ElasticLinearFrameSection3d::recvSelf -- failed to receive data\n";
      return res;
    }

    this->setTag((int)data(0));
    E  = data(1);
//  A  = data(2);
//  Iz = data(3);
//  Iy = data(4);
    G  = data(5);
//  J  = data(6);

    return res;
}

void
ElasticLinearFrameSection3d::Print(OPS_Stream &s, int flag)
{

  FrameSectionConstants consts{};
  getConstants(consts);

  double J = consts.Iy + consts.Iz - consts.Ca;

  if (flag == OPS_PRINT_PRINTMODEL_SECTION) {
    s << "ElasticLinearFrameSection3d, tag: " << this->getTag() << "\n";
    s << "\t E: " << E << "\n";
    s << "\t G: " << G         << "\n";
    s << "\t A: " << consts.A << "\n";
    s << "\tIz: " << consts.Iz << "\n";
    s << "\tIy: " << consts.Iy << "\n";
    s << "\t J: " <<        J << "\n";
  }

  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << OPS_PRINT_JSON_MATE_INDENT << "{";
    s << "\"name\": " << this->getTag() << ", ";
    s << "\"type\": \"ElasticFrameSection3d\", ";
    s << "\"E\": "   << E  << ", ";
    s << "\"G\": "   << G  << ", ";
    s << "\"A\": "   << consts.A   << ", ";
    s << "\"Ay\": "  << consts.Ay  << ", ";
    s << "\"Az\": "  << consts.Az  << ", ";
    s << "\"Iy\": " << consts.Iy << ", ";
    s << "\"Iz\": " << consts.Iz << ", ";

    s << "\"Jx\": " <<        J << ", ";
    s << "\"Ca\": " << consts.Ca << ", ";
    s << "\"Cw\": " << consts.Cw ;

    double mass;
    if (this->FrameSection::getIntegral(Field::Density, State::Init, mass) == 0) {
      s << ", ";
      s << "\"mass\": " << mass;
    }
    s << "}";
  }
}

int
ElasticLinearFrameSection3d::setParameter(const char **argv, int argc, Parameter &param)
{
  if (argc < 1)
    return -1;

  FrameSectionConstants consts;
  getConstants(consts);

  if (strcmp(argv[0],"E") == 0) {
    param.setValue(E);
    return param.addObject(1, this);
  }
  if (strcmp(argv[0],"G") == 0) {
    param.setValue(G);
    return param.addObject(5, this);
  }
  if (strcmp(argv[0],"A") == 0) {
    param.setValue(consts.A);
    return param.addObject(2, this);
  }
  if (strcmp(argv[0],"Iz") == 0) {
    param.setValue(consts.Iz);
    return param.addObject(3, this);
  }
  if (strcmp(argv[0],"Iy") == 0) {
    param.setValue(consts.Iy);
    return param.addObject(4, this);
  }
  if (strcmp(argv[0],"J") == 0) {
    double J = consts.Iy + consts.Iz - consts.Ca;
    param.setValue(J);
    return param.addObject(6, this);
  }
  return -1;
}

int
ElasticLinearFrameSection3d::updateParameter(int paramID, Information &info)
{


  FrameSectionConstants consts;
  getConstants(consts);

  // n-n
  double &A=consts.A;
  double &Ay=consts.Ay, &Az=consts.Az;
  // m-m
  double &Iy=consts.Iy, &Iz=consts.Iz, &Iyz=consts.Iyz;
  // w-w
  double &Cw=consts.Cw, &Ca=consts.Ca;
  // n-m
  double &Qy=consts.Qy, &Qz=consts.Qz;
  // n-w
  double &Rw=consts.Rw, &Ry=consts.Ry, &Rz=consts.Rz;
  // m-w
  double &Sa=consts.Sa, &Sy=consts.Sy, &Sz=consts.Sz;


  if (paramID == 1)
    E = info.theDouble;
  if (paramID == 2)
    A = info.theDouble;
  if (paramID == 3)
    Iz = info.theDouble;
  if (paramID == 4)
    Iy = info.theDouble;
  if (paramID == 5)
    G = info.theDouble;
  if (paramID == 6) {
    double J = info.theDouble;
    Ca = Iy + Iz - J;
  }

  double I0 = Iy + Iz;
  *Ks = {
    //            N                       M                      W                 Q    
    //                        |                        |                 |                |
    {{  E*A,      0.,     0.,      0.,   E*Qy,  -E*Qz,   E*Rw  ,  0.,  0.,    0.,  0.,  0.},
     {    0.,    G*A,     0.,   -G*Qy,     0.,     0.,     0.  ,  0.,  0.,  G*Ry,  0.,  0.},
     {    0.,     0.,    G*A,    G*Qz,     0.,     0.,     0.  ,  0.,  0.,  G*Rz,  0.,  0.},
    //                        |                        |          0.,  0.,      ,  0.,  0.|
     {    0.,  -G*Qy,   G*Qz,    G*I0,     0.,     0.,     0.  ,  0.,  0.,  G*Sa,  0.,  0.},
     {  E*Qy,     0.,     0.,      0.,  E*Iy , -E*Iyz,     E*Sy,  0.,  0.,    0.,  0.,  0.},
     { -E*Qz,     0.,     0.,      0., -E*Iyz,  E*Iz ,     E*Sz,  0.,  0.,    0.,  0.,  0.},
    //                        |                        |          0.,  0.,      ,  0.,  0.|
     {  E*Rw,     0.,     0.,      0.,   E*Sy,   E*Sz,    E*Cw ,  0.,  0.,    0.,  0.,  0.},
     {    0.,     0.,     0.,      0.,     0.,     0.,      0. ,  0.,  0.,    0.,  0.,  0.},
     {    0.,     0.,     0.,      0.,     0.,     0.,      0. ,  0.,  0.,    0.,  0.,  0.},
     {    0.,   G*Ry,   G*Rz,    G*Sa,     0.,     0.,      0. ,  0.,  0.,  G*Ca,  0.,  0.},
     {    0.,     0.,     0.,      0.,     0.,     0.,      0. ,  0.,  0.,    0.,  0.,  0.},
     {    0.,     0.,     0.,      0.,     0.,     0.,      0. ,  0.,  0.,    0.,  0.,  0.},}
  };

  return 0;
}

int
ElasticLinearFrameSection3d::activateParameter(int paramID)
{
  parameterID = paramID;

  return 0;
}

const Vector&
ElasticLinearFrameSection3d::getStressResultantSensitivity(int gradIndex,
                                                           bool conditional)
{
  s.zero();

  FrameSectionConstants consts;
  getConstants(consts);

  if (parameterID == 1) { // E
    s[0] = consts.A*e[0];
    s[1] = consts.Iz*e[1];
    s[2] = consts.Iy*e[2];
  }
  if (parameterID == 2) // A
    s[0] = E*e[0];
  if (parameterID == 3) // Iz
    s[1] = E*e[1];
  if (parameterID == 4) // Iy
    s[2] = E*e[2];
  if (parameterID == 5) // G
    s[3] = consts.Ca*e[3];
  if (parameterID == 6) // J
    s[3] = G*e[3];

  v.setData(s);
  return v;
}

const Matrix&
ElasticLinearFrameSection3d::getInitialTangentSensitivity(int gradIndex)
{
  if (Ksen == nullptr)
    Ksen = new Matrix(nr,nr);

  return *Ksen;
}

