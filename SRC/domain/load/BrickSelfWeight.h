                                                                        
#ifndef BrickSelfWeight_h
#define BrickSelfWeight_h

// Written: ZHYang, UCDavis

// Purpose: This file contains the class definition for 8 node brick self weight load.

#include <ElementalLoad.h>

class BrickSelfWeight : public ElementalLoad
{
  public:
    BrickSelfWeight(int tag, int eleTag);
    BrickSelfWeight();    
    ~BrickSelfWeight();

    const Vector &getData(int &type, double loadFactor);

    int sendSelf(int commitTag, Channel &);  
    int recvSelf(int commitTag, Channel &,  FEM_ObjectBroker &);
    void Print(OPS_Stream &s, int flag);
	
  private:
    static Vector data;
};

#endif

