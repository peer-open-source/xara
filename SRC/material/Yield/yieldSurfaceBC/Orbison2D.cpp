// Orbison2D.cpp: implementation of the YieldSurfaceBC class.
//
//////////////////////////////////////////////////////////////////////

#include "Orbison2D.h"
#include <math.h>
#define ORBISON_CLASS_TAG -1
//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

Orbison2D::Orbison2D(int tag, double xcap, double ycap,
                     YS_Evolution &model)
:YieldSurface_BC2D(tag, ORBISON_CLASS_TAG, xcap, ycap, model)
{

}

Orbison2D::~Orbison2D()
{

}

//////////////////////////////////////////////////////////////////////
// YS specific methods
//////////////////////////////////////////////////////////////////////

void Orbison2D::setExtent()
{
	// Extent along the axis
	xPos =  1;
	xNeg = -1;
	yPos =  1;
	yNeg = -1;
}


void Orbison2D::getGradient(double &gx, double &gy, double x, double y)
{
// check if the point is on the surface
    double drift =  getDrift(x, y);
    double capx = capXdim;
    double capy = capYdim;
    
    if(forceLocation(drift)!=0)
    {
     	opserr << "ERROR - Orbison2D::getGradient(double &gx, double &gy, double x, double y)\n";
        opserr << "Force point not on the yield surface\n";
		opserr << " fx = " << x << ", fy = " << y  << " drift = " << drift << "\n";
        opserr << "\a";
    }
    else
    {
    	gx = 2*x/(capx) + 7.34*pow(y, 2)*(x/(capx));
    	gy = 2.3*y/(capy) - 0.9*pow(y, 5)/(capy) + 7.34*pow(x, 2)*(y/(capy));
//  	p1 = 2.3*p - 0.9*pow(p, 5) + 7.34*pow(m, 2)*(p);
//  	m1 = 2*m + 7.34*pow(p, 2)*(m);

//      gx = 2*x + 7.34*pow(y, 2)*x;
//      gy = 2.3*y - 0.9*pow(y, 5) + 7.34*pow(x, 2)*y;

    }

}

double Orbison2D::getSurfaceDrift(double x, double y)
{
double phi = 1.15*y*y - 0.15*pow(y, 6) + x*x + 3.67*y*y*x*x;
double drift = phi - 1;
	return drift;
}

YieldSurface_BC *Orbison2D::getCopy(void)
{
    Orbison2D *theCopy = new Orbison2D(this->getTag(), capX, capY, *hModel);
    if(theCopy==0)
    {
     	opserr << "Orbison2D::getCopy(void) - unable to make copy\n";
     	opserr << "\a";
    }
    //later  copy all the state variables
    return theCopy;
}


void Orbison2D::Print(OPS_Stream &s, int flag)
{
    s << "\nYield Surface No: " << this->getTag() << " type: Attalla2D\n";
}
