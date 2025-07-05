// MultiLinear.h: interface for the MultiLinear class.
//
//////////////////////////////////////////////////////////////////////


#ifndef MultiLinearKp_h
#define MultiLinearKp_h

#include <PlasticHardeningMaterial.h>
#include <math.h>

class MultiLinearKp : public PlasticHardeningMaterial
{
public:
	MultiLinearKp(int tag, Vector &sum_plas_defo, Vector &kp);
	virtual ~MultiLinearKp();
	
	double getTrialPlasticStiffness();
    PlasticHardeningMaterial *getCopy();
    void Print(OPS_Stream &, int flag) final;

  private:
    Vector sumPlasDefo;
    Vector Kp;
    int numPoints;
};

#endif
