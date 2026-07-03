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


namespace {

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

enum class Parameters : int {
  E  = 1,
  G  = 5,
  //
  A  = 2,
  Ay = 17,
  Az = 18,
  //
  Iy = 4,
  Iz = 3,
  Iyz = 19,
  //
  J = 6,
  // n-m
  Qy = 7,
  Qz = 8,
  Rw = 9,
  Ry = 10,
  Rz = 11,
  Sa = 12,
  Sy = 13,
  Sz = 14,
  Cw = 15,
  None = 0
};

void
SetupTangent(MatrixND<12,12>& Ks, const Frame::Shape& cons)
{
  double E = *cons.E;
  double G = *cons.G;
  // n-n
  double A   =  cons.A;
  // m-m
  double Iy  = *cons.Iy;    //   \int z^2
  double Iz  = *cons.Iz;    //   \int y^2
  double Iyz =  cons.Iyz;   //
  // w-w 
  double Cw = cons.Cw;     // 

  // v-v
  double Ca = cons.Ca;     //
  double Sy = cons.Sy;     // m-w, int  z warp(x,y)
  double Sz = cons.Sz;     // m-w, int -y warp(x,y) 
  // n-m
  double Qy = cons.Qy;     //  int y = A*zs
  double Qz = cons.Qz;     //  int z

  double Rw = cons.Rw;     // n-w
  double Ry = cons.Ry;     // n-v
  double Rz = cons.Rz;

  double Sa = cons.Sa;     // m-v


  // Polar moment of inertia
  const double I0   = Iy + Iz;

  double Ay  = cons.Ay? *cons.Ay : A;
  double Az  = cons.Az? *cons.Az : A;
  double GSnyy = G*(Ay - A);
  double GSnzz = G*(Az - A);
  double GSy   = -GSnyy;
  double GSz   = -GSnzz;

  double GAy,GAz,GJ;
  if (!getenv("XARA_OLD_WARP")) {
    switch (cons.mixed_form) {
      case Frame::Shape::MixedType::Equilibrium:
      case Frame::Shape::MixedType::None:
        GAy = G*A;
        GAz = G*A;
        GJ  = G*I0;
        break;
      case Frame::Shape::MixedType::UT:
        GAy = G*A;
        GAz = G*A;
        GJ = G*(*cons.J);
        break;
      case Frame::Shape::MixedType::Energetic:
      case Frame::Shape::MixedType::Constant:
        GAy = G*Ay;
        GAz = G*Az;
        GJ  = G*(*cons.J);
        break;
    }
  }
  else {
    GAy = G*A;
    GAz = G*A;
    GJ  = G*I0;
  }


  Ks = {
    //            N                       M                       W                   Q    
    //                        |                        |                  |                     |
    {    E*A,     0.,     0.,      0.,   E*Qy,  -E*Qz,    E*Rw ,  0.,  0.,    0.,     0.,     0.,
          0.,    GAy,     0.,   -G*Qy,     0.,     0.,     0.  ,  0.,  0.,  G*Ry,  GSnyy,     0.,
          0.,     0.,    GAz,    G*Qz,     0.,     0.,     0.  ,  0.,  0.,  G*Rz,     0.,  GSnzz,
    //                        |                        |
          0.,  -G*Qy,   G*Qz,      GJ,     0.,     0.,     0.  ,  0.,  0.,  G*Sa,     0.,     0.,
        E*Qy,     0.,     0.,      0.,  E*Iy , -E*Iyz,     E*Sy,  0.,  0.,    0.,     0.,     0.,
       -E*Qz,     0.,     0.,      0., -E*Iyz,  E*Iz ,     E*Sz,  0.,  0.,    0.,     0.,     0.,
    //                        |                        | 
        E*Rw,     0.,     0.,      0.,   E*Sy,   E*Sz,    E*Cw ,  0.,  0.,    0.,     0.,     0., // Bimoment
          0.,     0.,     0.,      0.,     0.,     0.,      0. ,  0.,  0.,    0.,     0.,     0.,
          0.,     0.,     0.,      0.,     0.,     0.,      0. ,  0.,  0.,    0.,     0.,     0.,
    //                        |                        |          
          0.,   G*Ry,   G*Rz,    G*Sa,     0.,     0.,      0. ,  0.,  0.,  G*Ca,     0.,     0., // Bishear
          0.,  GSnyy,     0.,      0.,     0.,     0.,      0. ,  0.,  0.,    0.,    GSy,     0.,
          0.,     0.,  GSnzz,      0.,     0.,     0.,      0. ,  0.,  0.,    0.,     0.,    GSz}
  };
    //                        |                        |                  |
    //    0       1       2   |    3       4       5   |     6     7   8  |  9      10      11
}
}

ID ElasticLinearFrameSection3d::layout(layout_array, nr, false);



ElasticLinearFrameSection3d::ElasticLinearFrameSection3d(int tag,
    const Frame::Shape& shape,
    double mass_,
    bool use_mass
)
: FrameSection(tag, SEC_TAG_ElasticLinearFrame3d, mass_, use_mass),
  Ksen(nullptr),
  shape_data(std::make_shared<Frame::Shape>(shape)),
  Ks(new MatrixND<nr,nr> {}),
  e{},
  s{},
  parameterID(0)
{
  centroid = {shape.Qz/shape.A, 
              shape.Qy/shape.A};
  SetupTangent(*Ks, shape);
}


ElasticLinearFrameSection3d::ElasticLinearFrameSection3d(ElasticLinearFrameSection3d& other)
: FrameSection(other.getTag(), SEC_TAG_ElasticLinearFrame3d),
  shape_data(other.shape_data),
  Ks(other.Ks),
  Ksen(nullptr),
  e{},
  s{},
  parameterID(0)
{

}


void
ElasticLinearFrameSection3d::getConstants(Frame::Shape& consts) const
{
  consts = *shape_data;
  return;
}

int
ElasticLinearFrameSection3d::getIntegral(Field field, State state, double& value) const
{
  const Frame::Shape& consts = *shape_data;

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
      value = *consts.Iz;
      if (field == Field::UnitCentroidYY) {
        double yc = centroid[0];
        value -= consts.A*yc*yc;
      }
      return 0;

    case Field::UnitZZ:
    case Field::UnitCentroidZZ:
      value = *consts.Iy;
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
  // - overload 
  //   template<int n> ID::operator==(std::array<int, n>)
  //
  // OR
  // - add FrameSection::setLayout()
  
  ElasticLinearFrameSection3d *theCopy = new ElasticLinearFrameSection3d(*this);


  int ni=0;
  bool ind[nr]{};
  double K_data[nr][nr]{}, F_data[nr][nr]{};


  bool uniform_twist = true; // NOTE
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
  Matrix Ki(&K_data[0][0], ni, ni),
         Fi(&F_data[0][0], ni, ni);
  int ii=0;
  for (int i=0; i<nr; i++)
    if (ind[i]) {
      int jj=0;
      for (int j=0; j<nr; j++)
        if (ind[j])
          Ki(ii,jj++) = Kc(i,j);

      ii++;
    }

  Ki.Invert(Fi);

  theCopy->Fs = new Matrix(nr,nr);
  Matrix& Fc = *(theCopy->Fs);
  Fc.Zero();
  ii = 0;
  for (int i=0; i<nr; i++)
    if (ind[i]) {
      int jj=0;
      for (int j=0; j<nr; j++)
        if (ind[j])
          Fc(i,j) = Fi(ii,jj++);

      ii++;
    }

  return theCopy;
}

FrameSection*
ElasticLinearFrameSection3d::getFrameCopy()
{
  ElasticLinearFrameSection3d *theCopy = new ElasticLinearFrameSection3d(*this);

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
ElasticLinearFrameSection3d::setParameter(const char **argv, int argc, Parameter &param)
{
  if (argc < 1)
    return -1;

  const Frame::Shape& consts = *shape_data;

  if (strcmp(argv[0],"E") == 0) {
    param.setValue(*consts.E);
    return param.addObject(static_cast<int>(Parameters::E), this);
  }
  if (strcmp(argv[0],"G") == 0) {
    param.setValue(*consts.G);
    return param.addObject(static_cast<int>(Parameters::G), this);
  }
  if (strcmp(argv[0],"A") == 0) {
    param.setValue(consts.A);
    return param.addObject(static_cast<int>(Parameters::A), this);
  }
  if (strcmp(argv[0],"Iz") == 0) {
    param.setValue(*consts.Iz);
    return param.addObject(static_cast<int>(Parameters::Iz), this);
  }
  if (strcmp(argv[0],"Iy") == 0) {
    param.setValue(*consts.Iy);
    return param.addObject(static_cast<int>(Parameters::Iy), this);
  }
  if (strcmp(argv[0],"J") == 0) {
    param.setValue(*consts.J);
    return param.addObject(static_cast<int>(Parameters::J), this);
  }
  return -1;
}

int
ElasticLinearFrameSection3d::updateParameter(int paramID, Information &info)
{
  Ks = std::make_shared<MatrixND<nr,nr>>(*Ks); // copy on write
  shape_data = std::make_shared<Frame::Shape>(*shape_data); // copy on write
  Frame::Shape& consts = *shape_data;

  Parameters parameter = static_cast<Parameters>(paramID);

  if (parameter == Parameters::E)
    consts.E = info.theDouble;
  if (parameter == Parameters::G)
    consts.G = info.theDouble;

  // n-n
  if (parameter == Parameters::A)
    consts.A = info.theDouble;
  if (parameter == Parameters::Ay)
    consts.Ay = info.theDouble;
  if (parameter == Parameters::Az)
    consts.Az = info.theDouble;
  // m-m
  if (parameter == Parameters::Iyz)
    consts.Iyz = info.theDouble;
  if (parameter == Parameters::Iy)
    consts.Iy = info.theDouble;
  if (parameter == Parameters::Iz)
    consts.Iz = info.theDouble;
  
  //
  if (parameter == Parameters::J) {
    consts.J = info.theDouble;
  }
  // n-m
  if (parameter == Parameters::Qy)
    consts.Qy = info.theDouble;
  if (parameter == Parameters::Qz)
    consts.Qz = info.theDouble;

  SetupTangent(*Ks, consts);

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
  static VectorND<nr> ds;
  static Vector wrapper(ds);
  ds.zero();

  if (parameterID == 0)
    return wrapper; // no sensitivity

  
  Frame::Shape& C = *shape_data;
  Frame::Shape dC(shape_data->ndm, shape_data->ndf);
  if (C.Ay)  dC.Ay = 0.0;
  if (C.Az)  dC.Az = 0.0;
  if (C.J)   dC.J  = 0.0;
  if (C.Iy)  dC.Iy = 0.0;
  if (C.Iz)  dC.Iz = 0.0;
  Parameters parameter = static_cast<Parameters>(parameterID);
  if (parameter == Parameters::E) {
    dC = C;
    dC.E = 1.0;
    dC.G = 0.0;
  }
  else if (parameter == Parameters::G) {
    dC = C;
    dC.E = 0.0;
    dC.G = 1.0;
  }
  else {
    dC.E = C.E;
    dC.G = C.G;
    if (parameter == Parameters::A) {
      dC.A  = 1.0;
    }
    if (parameter == Parameters::Iz)
      dC.Iz = 1.0;
    if (parameter == Parameters::Iy)
      dC.Iy = 1.0;
    if (parameter == Parameters::J) {
      dC.J  = 1.0;
    }
    if (parameter == Parameters::Qy)
      dC.Qy = 1.0;
    if (parameter == Parameters::Qz)
      dC.Qz = 1.0;
    if (parameter == Parameters::Cw)
      dC.Cw = 1.0;
  }

  MatrixND<nr,nr> dK{};
  SetupTangent(dK, dC);
  ds = dK*e;

  return wrapper;
}

const Matrix&
ElasticLinearFrameSection3d::getInitialTangentSensitivity(int gradIndex)
{
  if (Ksen == nullptr)
    Ksen = new Matrix(nr,nr);

  return *Ksen;
}


int
ElasticLinearFrameSection3d::sendSelf(int commitTag, Channel &theChannel)
{
  return -1;
}


int
ElasticLinearFrameSection3d::recvSelf(int commitTag, Channel &theChannel,
                                      FEM_ObjectBroker &theBroker)
{
  return -1;
}

void
ElasticLinearFrameSection3d::Print(OPS_Stream &s, int flag)
{

  const Frame::Shape& consts = *shape_data;

  if (flag == OPS_PRINT_PRINTMODEL_SECTION) {
    s << "ElasticLinearFrameSection3d, tag: " << this->getTag() << "\n";
    s << "\t E: " << *consts.E << "\n";
    s << "\t G: " << *consts.G         << "\n";
    s << "\t A: " << consts.A << "\n";
    s << "\tIz: " << *consts.Iz << "\n";
    s << "\tIy: " << *consts.Iy << "\n";
    s << "\t J: " << *consts.J << "\n";
  }

  else if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << OPS_PRINT_JSON_MATE_INDENT << "{";
    s << "\"name\": " << this->getTag() << ", ";
    s << "\"type\": \"" << this->getClassType() <<"\", ";
    s << "\"E\": "   << *consts.E  << ", ";
    s << "\"G\": "   << *consts.G  << ", ";
    s << "\"A\": "   << consts.A   << ", ";
    s << "\"Ay\": "  << *consts.Ay  << ", ";
    s << "\"Az\": "  << *consts.Az  << ", ";
    s << "\"Iy\": "  << *consts.Iy << ", ";
    s << "\"Iz\": "  << *consts.Iz << ", ";

    s << "\"Jx\": " << *consts.J << ", ";
    s << "\"Ca\": " << consts.Ca << ", ";
    s << "\"Cw\": " << consts.Cw ;

    double mass;
    if (this->FrameSection::getIntegral(Field::Density, State::Init, mass) == 0) {
      s << ", ";
      s << "\"mass\": " << mass;
    }
    s << "}";
  }

  else if (flag == OPS_PRINT_CURRENTSTATE) {
    s << "  e: " << Vector(e);
    s << "  s: " << Vector((*Ks)*e);
  }
}