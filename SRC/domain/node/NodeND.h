#include "Node.h"
#include <VectorND.h>
#include <assert.h>

template <int ndm, int ndf>
class NodeND: public Node {
  public:
  NodeND(int tag, double crd)                           : Node(tag, ndf) {static_assert(ndm == 1); createDisp();}
  NodeND(int tag, double crd1, double crd2)             : Node(tag, ndf) {static_assert(ndm == 2); createDisp();}
  NodeND(int tag, double crd1, double crd2, double crd3): Node(tag, ndf) {static_assert(ndm == 3); createDisp();}


  virtual int setTrialDisp  (double value, int dof) final {
    assert(dof >= 0 && dof < ndf);
    double tDisp = value;
    displ[dof+2*ndf] = tDisp - displ[dof+ndf];
    displ[dof+3*ndf] = tDisp - displ[dof];	
    displ[dof]       = tDisp;
    return 0;
  }

  virtual int setTrialDisp  (const Vector & newTrialDisp) final {
    assert(newTrialDisp.Size() == ndf);

    for (int i=0; i<ndf; i++) {
        double tDisp = newTrialDisp(i);
        displ[i+2*ndf] = tDisp - displ[i+ndf];
        displ[i+3*ndf] = tDisp - displ[i];	
        displ[i] = tDisp;
    }
    return 0;
  }

  virtual int incrTrialDisp (const Vector &incrDispl) override final {
    assert(incrDispl.Size() == ndf);
    // set trial = incr + trial
    for (int i = 0; i<ndf; i++) {
        double incrDispI = incrDispl(i);
        displ[i]       += incrDispI;
        displ[i+2*ndf] += incrDispI;
        displ[i+3*ndf]  = incrDispI;
    }
    return 0;
  }

  virtual int commitState() override final {
      // set commit = trial, incr = 0.0
      for (int i=0; i<ndf; i++) {
        displ[i+ndf] = displ[i];  
        displ[i+2*ndf] = 0.0;
        displ[i+3*ndf] = 0.0;
      }

      // check vel exists, if does set commit = trial    
      if (trialVel != 0) {
        for (int i=0; i<ndf; i++)
          vel[i+ndf] = vel[i];
      }

      // check accel exists, if does set commit = trial        
      if (trialAccel != 0) {
        for (int i=0; i<ndf; i++)
          accel[i+ndf] = accel[i];
      }

      return 0;
  }

  virtual int revertToLastCommit() final {
      // check disp exists, if does set trial = last commit, incr = 0
      for (int i=0 ; i<ndf; i++) {
        displ[i] = displ[i+ndf];
        displ[i+2*ndf] = 0.0;
        displ[i+3*ndf] = 0.0;
      }
    
      // check vel exists, if does set trial = last commit
      if (vel != nullptr) {
        for (int i=0 ; i<ndf; i++)
          vel[i] = vel[ndf+i];
      }

      // check accel exists, if does set trial = last commit
      if (accel != nullptr) {
        for (int i=0 ; i<ndf; i++)
          accel[i] = accel[ndf+i];
      }

      return 0;
  }

  virtual int revertToStart() override final {
      // set all to zero
      for (int i=0 ; i<4*ndf; i++)
        displ[i] = 0.0;

      // check vel exists, if does set all to zero
      if (vel != 0) {
        for (int i=0 ; i<2*ndf; i++)
          vel[i] = 0.0;
      }

      // check accel exists, if does set all to zero
      if (accel != 0) {    
        for (int i=0 ; i<2*ndf; i++)
          accel[i] = 0.0;
      }
      
      if (unbalLoad != nullptr)
          (*unbalLoad) *= 0;

      return 0;
  }

  // MovableObject
  virtual int sendSelf(int commitTag, Channel &theChannel) {
    return -1;
  }

  virtual int recvSelf(int commitTag, Channel &theChannel, 
                       FEM_ObjectBroker &theBroker) {
    return -1;
  }

  private:
    double *disp;
    double *vel, *accel;              // double arrays holding the vel and accel values
    Vector *trialDisp, *commitDisp, *incrDisp, *incrDeltaDisp;
    Vector *commitVel, *commitAccel;  // committed quantities
    Vector *trialVel,  *trialAccel;   // trial quantities
    Vector *unbalLoad;                // unbalanced load

    int createDisp() {
      trialDisp     = new Vector(displ, ndf);
      commitDisp    = new Vector(&displ[ndf], ndf); 
      incrDisp      = new Vector(&displ[2*ndf], ndf);
      incrDeltaDisp = new Vector(&displ[3*ndf], ndf);
      return 0;
    }

//  OpenSees::VectorND<4*ndf> displ = {{0.0}};
    double displ[4*ndf] = {{0.0}};
};

