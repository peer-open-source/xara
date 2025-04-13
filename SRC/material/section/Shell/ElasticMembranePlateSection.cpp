/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
** ****************************************************************** */
//
// Elastic Plate Section with membrane
//
// Ed "C++" Love
//
// Out-of-Plane stiffness modifier added by Pearl Ranchal
// Supported by Degenkolb Engineers
//

#include <ElasticMembranePlateSection.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>


const double ElasticMembranePlateSection::five6 = 5.0/6.0 ; // shear correction

Vector  ElasticMembranePlateSection::stress(8) ;
Matrix  ElasticMembranePlateSection::tangent(8,8) ;
ID      ElasticMembranePlateSection::array(8) ;


ElasticMembranePlateSection::ElasticMembranePlateSection( ) : 
SectionForceDeformation( 0, SEC_TAG_ElasticMembranePlateSection ), 
strain(8) 
{ 

}


ElasticMembranePlateSection::ElasticMembranePlateSection(int    tag,
                                                         double young_membrane,
                                                         double poisson,
                                                         double thickness,
                                                         double r, 
                                                         double young_plate_mod)
  :
  SectionForceDeformation(tag, SEC_TAG_ElasticMembranePlateSection),
  strain(8)
{
    this->Em    = young_membrane;
    this->Ep    = young_membrane * young_plate_mod;
    this->nu    = poisson;
    this->h     = thickness;
    this->rhoH  = r * thickness;
}


ElasticMembranePlateSection::~ElasticMembranePlateSection() 
{ 

} 


SectionForceDeformation*  
ElasticMembranePlateSection::getCopy() 
{
  ElasticMembranePlateSection *clone ;   

  clone = new ElasticMembranePlateSection(this->getTag(), Em, nu, h, rhoH, Ep/Em) ; //new instance of this class

  //    *clone = *this ; //assignment to make copy
  clone->rhoH = this->rhoH ;
  clone->strain = this->strain;

  return clone ;
}

// density per unit area
double
ElasticMembranePlateSection::getRho()
{
  return rhoH;
}


int
ElasticMembranePlateSection::getOrder() const
{
  return 8 ;
}


//send back order of strain in vector form
const ID& ElasticMembranePlateSection::getType()
{
    static bool initialized = false;
    if (!initialized) {
        array(0) = SECTION_RESPONSE_FXX;
        array(1) = SECTION_RESPONSE_FYY;
        array(2) = SECTION_RESPONSE_FXY;
        array(3) = SECTION_RESPONSE_MXX;
        array(4) = SECTION_RESPONSE_MYY;
        array(5) = SECTION_RESPONSE_MXY;
        array(6) = SECTION_RESPONSE_VXZ;
        array(7) = SECTION_RESPONSE_VYZ;
        initialized = true;
    }
    return array;
}



int
ElasticMembranePlateSection::commitState() 
{
  return 0 ;
}


int
ElasticMembranePlateSection::revertToLastCommit()
{
  return 0 ;
}


int
ElasticMembranePlateSection::revertToStart()
{
  return 0 ;
}



int
ElasticMembranePlateSection::setTrialSectionDeformation(const Vector &strain_from_element)
{
  this->strain = strain_from_element ;

  return 0 ;
}


const Vector&
ElasticMembranePlateSection::getSectionDeformation()
{
  return this->strain ;
}


const Vector&
ElasticMembranePlateSection::getStressResultant()
{

  double M  = Em / ( 1.0 - nu*nu ) ; // membrane modulus

  double G  =  0.5 * Em / ( 1.0 + nu ) ; // shear modulus
 
  G *= h ;
  M *= h ;

  //membrane resultants

  stress(0) =  M*strain(0) + (nu*M)*strain(1)  ;
 
  stress(1) =  (nu*M)*strain(0) +  M*strain(1)  ;

  stress(2) =  G*strain(2) ;

 

  G *= (five6 * (Ep / Em));  //multiply by product of shear correction factor and ratio of bending to membrane moduli

  double D  =  Ep * (h*h*h) / 12.0 / ( 1.0 - nu*nu ) ;  // bending modulus

  // bending resultants

  stress(3) = -( D*strain(3) + nu*D*strain(4) ) ;
 
  stress(4) = -( nu*D*strain(3) + D*strain(4) ) ;

  stress(5) = -0.5*D*( 1.0 - nu )*strain(5) ;

  stress(6) = G*strain(6) ;

  stress(7) = G*strain(7) ;

 
  return this->stress ;
}


const Matrix&
ElasticMembranePlateSection::getSectionTangent()
{

  double M  = Em / ( 1.0 - nu*nu ) ; //membrane modulus

  double G  =  0.5 * Em / ( 1.0 + nu ) ; //shear modulus

  G *= h ;  //multiply by thickness
  M *= h ;

  tangent.Zero() ;

  //membrane tangent terms

  tangent(0,0) = M ;
  tangent(1,1) = M ;

  tangent(0,1) = nu*M ;
  tangent(1,0) = tangent(0,1) ;

  tangent(2,2) = G ;



  G *= (five6 * (Ep / Em));  //multiply by product of shear correction factor and ratio of bending to membrane moduli

  double D  =  Ep * (h*h*h) / 12.0 / ( 1.0 - nu*nu ) ;  //bending modulus

  // bending tangent terms

  tangent(3,3) = -D ;
  tangent(4,4) = -D ;

  tangent(3,4) = -nu*D ;
  tangent(4,3) = tangent(3,4) ;

  tangent(5,5) = -0.5 * D * ( 1.0 - nu ) ;

  tangent(6,6) = G ;

  tangent(7,7) = G ;

  return this->tangent ;
}


//send back the initial tangent 
const Matrix&
ElasticMembranePlateSection::getInitialTangent()
{

  double M  = Em / ( 1.0 - nu*nu ) ; // membrane modulus

  double G  =  0.5 * Em / ( 1.0 + nu ) ; // shear modulus

  G *= h ;  //multiply by thickness
  M *= h ;

  tangent.Zero() ;

  //membrane tangent terms

  tangent(0,0) = M ;
  tangent(1,1) = M ;

  tangent(0,1) = nu*M ;
  tangent(1,0) = tangent(0,1) ;

  tangent(2,2) = G ;



  G *= (five6 * (Ep / Em));  //multiply by product of shear correction factor and ratio of bending to membrane moduli

  double D  =  Ep * (h*h*h) / 12.0 / ( 1.0 - nu*nu ) ;  //bending modulus

  //bending tangent terms

  tangent(3,3) = -D ;
  tangent(4,4) = -D ;

  tangent(3,4) = -nu*D ;
  tangent(4,3) = tangent(3,4) ;

  tangent(5,5) = -0.5 * D * ( 1.0 - nu ) ;

  tangent(6,6) = G ;

  tangent(7,7) = G ;

  return this->tangent ;
}


//print out data
void  ElasticMembranePlateSection::Print( OPS_Stream &s, int flag )
{
    if (flag == OPS_PRINT_PRINTMODEL_SECTION) {
        s << "ElasticMembranePlateSection: \n ";
        s << "  Young's Modulus for Membrane (in-plane) Action, Em = " << Em << endln;
        s << "  Young's Modulus for Plate (out-of-plane) Action, Ep = " << Ep << endln;
        s << "  Poisson's Ratio nu = " << nu << "\n";
        s << "  Thickness h = " << h << "\n";
        s << "  Density rho = " << (rhoH/h) << "\n";
    }
    
    if (flag == OPS_PRINT_PRINTMODEL_JSON) {
        s << OPS_PRINT_JSON_ELEM_INDENT << "{";
        s << "\"name\": \"" << this->getTag() << "\", ";
        s << "\"type\": \"ElasticMembranePlateSection\", ";
        s << "\"Em\": " << Em << ", ";
        s << "\"Ep\": " << Ep << ", ";
        s << "\"nu\": " << nu << ", ";
        s << "\"thickness\": " << h << ", ";
        s << "\"masspervolume\": " << (rhoH/h) << "}";
    }
}


int
ElasticMembranePlateSection::sendSelf(int cTag, Channel &theChannel) 
{
  int res = 0;
  static Vector data(6);
  data(0) = this->getTag();
  data(1) = Em;
  data(2) = nu;
  data(3) = h;
  data(4) = rhoH;
  data(5) = Ep/Em;

  res = theChannel.sendVector(this->getDbTag(), cTag, data);
  if (res < 0) 
    opserr << "ElasticMembranePlateSection::sendSelf() - failed to send data\n";

  return res;
}


int
ElasticMembranePlateSection::recvSelf(int cTag, Channel &theChannel, 
				      FEM_ObjectBroker &theBroker)
{
  int res = 0;
  static Vector data(6);
  res = theChannel.recvVector(this->getDbTag(), cTag, data);
  if (res < 0) 
    opserr << "ElasticMembranePlateSection::recvSelf() - failed to recv data\n";
  else {
    this->setTag(data(0));
    Em   = data(1);
    Ep   = data(1) * data(5);
    nu   = data(2);
    h    = data(3);
    rhoH = data(4);
  }

  return res;
}
