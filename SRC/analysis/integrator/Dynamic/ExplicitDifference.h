
#ifndef ExplicitDifference_h
#define ExplicitDifference_h


#include<TransientIntegrator.h>

class DOF_Group;
class FE_Element;
class Vector;

class ExplicitDifference : public TransientIntegrator
{
public:
	ExplicitDifference();
	ExplicitDifference(double alphaM, double betaK, double betaKi, double betaKc);
	~ExplicitDifference();

	                                                 

	int formEleTangent(FE_Element *theEle);

	int formNodTangent(DOF_Group *theDof);
	
	const Vector & getVel();    //added for Modal damping

	int domainChanged();
	int newStep(double deltaT);
	int update(const Vector &U);

	int commit();

	void Print(OPS_Stream &s, int flag);


private:
	double deltaT;
	static double deltaT1;
	double alphaM;
	double betaK;
	double betaKi;
	double betaKc;

	int updateCount;
	Vector *U, *Ut;
	Vector  *Utdotdot, *Utdotdot1;
	Vector *Udot, *Utdot, *Utdot1;

};

#endif
