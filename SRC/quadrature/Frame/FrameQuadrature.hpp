#include <BeamIntegration.h>
#include <element/Frame/for_int.tpp>

template<typename T>
class FrameQuadrature : public BeamIntegration
{
  public:
    FrameQuadrature() : BeamIntegration(0) {}

    BeamIntegration *getCopy(void) override {
        return new FrameQuadrature<T>();
    }

    virtual void 
    getSectionLocations(int numSections, double L, double *xi) const override {
        
        static_loop<0, T::nip>([&](auto i) {
            if ((int)i.value >= numSections)
                return;
            xi[i.value] = 0.5*(1.0 + T::pts[i.value]);
        });
    }

    virtual void 
    getSectionWeights(int numSections, double L, double *wt) const override {
        static_loop<0, T::nip>([&](auto i) {
            if ((int)i.value >= numSections) 
                return;
            wt[i.value] = T::wts[i.value]*0.5;
        });
    }
    int sendSelf(int cTag, Channel &) {return 0;}
    int recvSelf(int cTag, Channel &, FEM_ObjectBroker &) {return 0;}

    virtual void Print(OPS_Stream &s, int flag) override {
        if (flag == OPS_PRINT_PRINTMODEL_JSON) {
            s << "{\"type\": \"FrameQuadrature\"}";
        }
        else {
            s << "FrameQuadrature" << endln;
        }
    }
};
