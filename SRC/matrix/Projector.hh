//===----------------------------------------------------------------------===//
//
//                                   xara
//
//===----------------------------------------------------------------------===//
//                              https://xara.so
//===----------------------------------------------------------------------===//
//
#include <MatrixND.h>

namespace OpenSees {

    // 2nd order Identity Tensor
    static constexpr MatrixND<3,3> I1 {{
        {1.0, 0.0, 0.0},
        {0.0, 1.0, 0.0},
        {0.0, 0.0, 1.0}}};

    // 4th order mixed variant identity tensor (51b)
    static constexpr MatrixND<6,6> IImix {{
        {1, 0, 0,  0, 0, 0},
        {0, 1, 0,  0, 0, 0},
        {0, 0, 1,  0, 0, 0},
        {0, 0, 0,  1, 0, 0},
        {0, 0, 0,  0, 1, 0},
        {0, 0, 0,  0, 0, 1},
    }};

    // 4th order covariant identity tensor (51a)
    static constexpr MatrixND<6,6> IIco {{
        {1, 0, 0,  0, 0, 0},
        {0, 1, 0,  0, 0, 0},
        {0, 0, 1,  0, 0, 0},
        {0, 0, 0,  2, 0, 0},
        {0, 0, 0,  0, 2, 0},
        {0, 0, 0,  0, 0, 2},
    }};

    // 4th order contravariant identity tensor (51a)
    static constexpr MatrixND<6,6> IIcon {{
        {1, 0, 0,   0 ,  0 ,  0 },
        {0, 1, 0,   0 ,  0 ,  0 },
        {0, 0, 1,   0 ,  0 ,  0 },
        {0, 0, 0,  0.5,  0 ,  0 },
        {0, 0, 0,   0 , 0.5,  0 },
        {0, 0, 0,   0 ,  0 , 0.5},
    }};

    // 4th order Volumetric Tensor (57)
    // IIvol = I1 tensor I1
    static constexpr MatrixND<6,6> IoI {{
        {1, 1, 1,  0, 0, 0},
        {1, 1, 1,  0, 0, 0},
        {1, 1, 1,  0, 0, 0},
        {0, 0, 0,  0, 0, 0},
        {0, 0, 0,  0, 0, 0},
        {0, 0, 0,  0, 0, 0},
    }};
    static constexpr MatrixND<6,6> IIvol {{
        {1, 1, 1,  0, 0, 0},
        {1, 1, 1,  0, 0, 0},
        {1, 1, 1,  0, 0, 0},
        {0, 0, 0,  0, 0, 0},
        {0, 0, 0,  0, 0, 0},
        {0, 0, 0,  0, 0, 0},
    }};

    static constexpr VectorND<6> ivol {{1, 1, 1, 0, 0, 0}};

    // 4th order Deviatoric Tensor
    //
    // Note:  this is the contravariant form!
    //        usable for s^a = 2G * IIdev^ab * epsilon_b
    // (Need a different form for s^a = IIdev ^a_b * sigma^a)
    static constexpr MatrixND<6,6> IIdev {{
        {  2./3.,  -1./3.,  -1./3.,   0 ,  0 ,  0 },
        { -1./3.,   2./3.,  -1./3.,   0 ,  0 ,  0 },
        { -1./3.,  -1./3.,   2./3.,   0 ,  0 ,  0 },
        {     0.,      0.,      0.,  0.5,  0 ,  0 },
        {     0.,      0.,      0.,   0 , 0.5,  0 },
        {     0.,      0.,      0.,   0 ,  0 , 0.5}
    }};

    // 4th order contravariant deviatoric tensor (Id)
    static constexpr MatrixND<6,6> IIdevCon = IIcon - 1./3.*IIvol;
    // 4th order covariant deviatoric tensor
    static constexpr MatrixND<6,6> IIdevCo  = IIco  - 1./3.*IIvol;
    // 4th order mixed variant deviatoric tensor (Idp)
    static constexpr MatrixND<6,6> IIdevMix = IImix - 1./3.*IIvol;
}
