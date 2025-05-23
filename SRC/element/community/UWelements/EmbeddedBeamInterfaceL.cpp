/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */

// Written: Alborz Ghofrani, Diego Turello, Pedro Arduino, U.Washington 
// Created: May 2017
// Description: This file contains the class definition for EmbeddedBeamInterfaceL.

#include <EmbeddedBeamInterfaceL.h>
#include <Node.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
#include <Renderer.h>
#include <Domain.h>
#include <string.h>
#include <Information.h>
#include <Parameter.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <ElementResponse.h>
#include <CrdTransf.h>
#include <elementAPI.h>
#include <cmath>
#include <NodeIter.h>
#include <FileStream.h>

#ifdef _PARALLEL_PROCESSING
#include <PartitionedDomain.h>
#endif

static int num_EmbeddedBeamInterfaceL = 0;
static const double m_Pi = 3.141592653589793;

void * OPS_ADD_RUNTIME_VPV(OPS_EmbeddedBeamInterfaceL)
{
    if (num_EmbeddedBeamInterfaceL == 0) {
        num_EmbeddedBeamInterfaceL++;
        opslog << "EmbeddedBeamInterfaceL element - Written: A.Ghofrani, D.Turello, P.Arduino, U.Washington\n";
    }

    Element *theElement = 0;

    int numArgs = OPS_GetNumRemainingInputArgs();
    if (numArgs < 1) {
        opserr << "Want: EmbeddedBeamInterfaceL tag? \n";
        return 0;
    }

    int iData[1];
    int numData = 1;
    if (OPS_GetIntInput(&numData, iData) != 0) {
        opserr << "WARNING invalid integer data: element EmbeddedBeamInterfaceL" << endln;
        return 0;
    }

    int eleTag = iData[0];

    theElement = new EmbeddedBeamInterfaceL(eleTag);

    if (theElement == 0) {
        opserr << "WARNING could not create element of type EmbeddedBeamInterfaceL\n";
        return 0;
    }

    return theElement;
}


EmbeddedBeamInterfaceL::EmbeddedBeamInterfaceL(int tag) :
  Element(tag, ELE_TAG_EmbeddedBeamInterfaceL),
  theSolidTags(0), solidNodeTags(0), theBeamTags(0), beamNodeTags(0), theNodes(0),
  crdTransf(0)
{

}

EmbeddedBeamInterfaceL::EmbeddedBeamInterfaceL(int tag, std::vector <int> beamTag, std::vector <int> solidTag, int crdTransfTag,
    std::vector <double>  beamRho,
    std::vector <double>  beamTheta,
    std::vector <double>  solidXi, 
    std::vector <double>  solidEta,
    std::vector <double>  solidZeta, 
    double radius, 
    std::vector <double> area, 
    std::vector <double> length, 
    bool writeConnectivity,
    const char * connectivityFN,
    Domain& theDomain
) : 
    Element(tag, ELE_TAG_EmbeddedBeamInterfaceL),
    m_beam_radius(radius), mQa(3, 3), mQb(3, 3), mQc(3, 3),
    mBphi(3, 12), mBu(3, 12), mHf(3, 12), m_Ns(8)
{
  
    // get domain to access element tags and their nodes
#if 0
#ifdef _PARALLEL_PROCESSING
    extern PartitionedDomain theDomain;
#else
    extern Domain theDomain;
#endif
#endif

    m_numEmbeddedPoints = solidTag.size();
    // theSolidTags  = new int[m_numEmbeddedPoints];
    solidNodeTags = new int[8 * m_numEmbeddedPoints];
    theBeamTags   = new int[m_numEmbeddedPoints];
    beamNodeTags  = new int[2 * m_numEmbeddedPoints];
    memoryallocated  = true;
    m_beam_rho = m_beam_theta = m_solid_xi = m_solid_eta = m_solid_zeta = m_area = m_beamLength = Vector(m_numEmbeddedPoints);

    std::set <int> uniqueSolidNodeTags;
    std::set <int> uniqueBeamNodeTags;
    std::set <int> uniqueBeamTags;
    Element *theElement;
    for (int ii = 0; ii < m_numEmbeddedPoints; ii++)
    {
        // theSolidTags[ii]    = solidTag[ii];
        theBeamTags[ii]     = beamTag[ii];
        m_solid_xi(ii)      = solidXi[ii];
        m_solid_eta(ii)     = solidEta[ii];
        m_solid_zeta(ii)    = solidZeta[ii];
        m_beam_rho(ii)      = beamRho[ii];
        m_beam_theta(ii)    = beamTheta[ii];
        m_area(ii)          = area[ii];
        m_beamLength(ii)    = length[ii];

        theElement = theDomain.getElement(solidTag[ii]);
        if (ii == 0)
            m_numSolidDOF = theElement->getNodePtrs()[0]->getNumberDOF();
        // opserr << "Point " << ii +1 << " : element " << solidTag[ii] << " at (" << solidXi[ii] << "," << solidEta[ii] << "," << solidZeta[ii] << ") , beam: " << beamTag << " at (" << beamRho[ii] << "," << beamTheta[ii] << ")" << endln;
        for (int jj = 0; jj < 8; jj++)
        {
            uniqueSolidNodeTags.insert(theElement->getNodePtrs()[jj]->getTag());
            solidNodeTags[ii * 8 + jj] = theElement->getNodePtrs()[jj]->getTag();
        }
        uniqueBeamTags.insert(beamTag[ii]);
        theElement = theDomain.getElement(beamTag[ii]);
        // opserr << "Point " << ii +1 << " : element " << solidTag[ii] << " at (" << solidXi[ii] << "," << solidEta[ii] << "," << solidZeta[ii] << ") , beam: " << beamTag << " at (" << beamRho[ii] << "," << beamTheta[ii] << ")" << endln;
        for (int jj = 0; jj < 2; jj++)
        {
            uniqueBeamNodeTags.insert(theElement->getNodePtrs()[jj]->getTag());
            beamNodeTags[ii * 2 + jj] = theElement->getNodePtrs()[jj]->getTag();
        }
    }


    m_numSolidNodes = (int)uniqueSolidNodeTags.size();
    m_numBeamNodes  = (int)uniqueBeamNodeTags.size();
    EBIL_numNodes   = m_numSolidNodes + m_numBeamNodes * 2;
    EBIL_numDOF     = m_numSolidNodes * 3 + m_numBeamNodes * 12;

    m_Lambda        = Vector(6 * m_numBeamNodes);
    m_solidInitDisp = Vector(m_numSolidNodes * 3);
    m_beamInitDisp  = Vector(m_numBeamNodes * 6);

    // create the Lagrange multiplier nodes
    NodeIter& theNodeIter = theDomain.getNodes();
    Node * theNode = theNodeIter();
    int maxTag = theNode->getTag();
    while ((theNode = theNodeIter()) != 0)
        if (maxTag < theNode->getTag())
            maxTag = theNode->getTag();

    for (int ii = 0; ii < m_numBeamNodes; ii++)
        theDomain.addNode(new Node(maxTag + ii + 1, 6, 0.0, 0.0, 0.0));



    externalNodes = ID(EBIL_numNodes);
    theNodes = new Node*[EBIL_numNodes];
    theNodesStatus = true;

    int count = 0;
    for (std::set<int>::iterator it = uniqueSolidNodeTags.begin(); it != uniqueSolidNodeTags.end(); ++it)
    {
        m_solidNodeMap[*it] = count;
        externalNodes(count) = *it;

        // theNodes[count] = theDomain.getNode(*it);

        // Vector tempDisp = theNodes[count]->getDisp();
        // m_solidInitDisp(count * 3 + 0) = tempDisp(0);
        // m_solidInitDisp(count * 3 + 1) = tempDisp(1);
        // m_solidInitDisp(count * 3 + 2) = tempDisp(2);

        count++;
    }

    int curCount = count;
    for (std::set<int>::iterator it = uniqueBeamNodeTags.begin(); it != uniqueBeamNodeTags.end(); ++it)
    {
        m_beamNodeMap[*it] = count - curCount;
        externalNodes(count) = *it;

        // theNodes[count] = theDomain.getNode(*it);

        // Vector tempDisp = theNodes[count]->getDisp();
        // m_beamInitDisp((count - curCount) * 6 + 0) = tempDisp(0);
        // m_beamInitDisp((count - curCount) * 6 + 1) = tempDisp(1);
        // m_beamInitDisp((count - curCount) * 6 + 2) = tempDisp(2);
        // m_beamInitDisp((count - curCount) * 6 + 3) = tempDisp(3);
        // m_beamInitDisp((count - curCount) * 6 + 4) = tempDisp(4);
        // m_beamInitDisp((count - curCount) * 6 + 5) = tempDisp(5);

        count++;
    }

    for (int ii = 0; ii < m_numBeamNodes; ii++)
    {
        int lagNodeTag = maxTag + ii + 1;
        externalNodes(count) = lagNodeTag;

        // theNodes[count] = theDomain.getNode(lagNodeTag);
        count++;
    }

    m_InterfaceForces = Vector(EBIL_numDOF);
    m_InterfaceStiffness = Matrix(EBIL_numDOF, EBIL_numDOF);
    mA = Matrix(3 * m_numSolidNodes, 6 * m_numBeamNodes);
    mB = Matrix(6 * m_numBeamNodes,  6 * m_numBeamNodes);


    // get the coordinate transformation object
    // TODO: Wont this segfault on a bad transform tag? do it before constructor - cmp
    crdTransf = G3_getSafeBuilder(rt)->getTypedObject<CrdTransf>(crdTransfTag)->getCopy3d();

    if (writeConnectivity)
    {
        FileStream connFile(connectivityFN, APPEND);
        connFile << m_numSolidNodes << " " << (int)uniqueBeamTags.size();
        for (int ii = 0; ii < m_numSolidNodes; ii++)
            connFile << " " << externalNodes(ii);
        int curBeamTag = theBeamTags[0];
        connFile << " " << beamNodeTags[0] << " " << beamNodeTags[1];
        for (int ii = 1; ii < m_numEmbeddedPoints; ii++)
            if (theBeamTags[ii] == curBeamTag)
                continue;
            else
            {
                curBeamTag = theBeamTags[ii];
                connFile << " " << beamNodeTags[2 * ii] << " " << beamNodeTags[2 * ii + 1];
            }
        connFile << endln;
        connFile.close();
    }

}

EmbeddedBeamInterfaceL::EmbeddedBeamInterfaceL()
    :
  Element(0, ELE_TAG_EmbeddedBeamInterfaceL),
  theSolidTags(0), solidNodeTags(0), theBeamTags(0), beamNodeTags(0), theNodes(0),
  m_beam_radius(0), mQa(3, 3), mQb(3, 3), mQc(3, 3),
  mBphi(3, 12), mBu(3, 12), mHf(3, 12), m_Ns(8), crdTransf(0)
{
    // mInitialize = false;
    theNodes = new Node*[2];
}

EmbeddedBeamInterfaceL::~EmbeddedBeamInterfaceL()
{
  if (theSolidTags != 0)
    delete [] theSolidTags;
  if (solidNodeTags != 0)
    delete [] solidNodeTags;
  if (theBeamTags != 0)
    delete [] theBeamTags;
  if (beamNodeTags != 0)
    delete [] beamNodeTags;

  if (crdTransf != 0)
    delete crdTransf;

  if (theNodes != 0)
    delete [] theNodes;
}

int
EmbeddedBeamInterfaceL::getNumExternalNodes(void) const
{
    return EBIL_numNodes;
}

const ID&
EmbeddedBeamInterfaceL::getExternalNodes(void)
{
    return externalNodes;
}

Node **
EmbeddedBeamInterfaceL::getNodePtrs(void)
{
    return theNodes;
}

int
EmbeddedBeamInterfaceL::getNumDOF(void)
{
    return EBIL_numDOF;
}

int
EmbeddedBeamInterfaceL::revertToLastCommit(void)
{
    return 0;
}

int
EmbeddedBeamInterfaceL::revertToStart(void)
{
    return 0;
}


const Matrix&
EmbeddedBeamInterfaceL::getTangentStiff(void)
{

    return m_InterfaceStiffness;

}

const Matrix&
EmbeddedBeamInterfaceL::getInitialStiff(void)
{
    return this->getTangentStiff();
}

const Vector&
EmbeddedBeamInterfaceL::getResistingForce(void)
{
    m_InterfaceForces.Zero();
    Vector temp2(3 * m_numSolidNodes), temp3(6 * m_numBeamNodes);

    temp2 = mA * m_Lambda;
    temp3 = -1.0 * (mB * m_Lambda);

    int II;
    for (int ii = 0; ii < 3 * m_numSolidNodes; ii++)
    {
        II = m_numSolidDOF * (ii / 3) + (ii % 3);
        m_InterfaceForces(II) = temp2(ii);
    }

    for (int ii = 0; ii < 6 * m_numBeamNodes; ii++)
        m_InterfaceForces(m_numSolidDOF * m_numSolidNodes + ii) = temp3(ii);

    return m_InterfaceForces;
}

int
EmbeddedBeamInterfaceL::sendSelf(int commitTag, Channel &theChannel)
{
    // opserr << "==================== sendlsef =======================\n";
    int res = 0;
    static ID idData(9);
    idData.Zero();
    idData(0) = this->getTag(); 
    idData(1) = m_numEmbeddedPoints;
    idData(2) = m_numSolidDOF;
    idData(3) = m_numSolidNodes;
    idData(4) = m_numBeamNodes;  
    idData(5) = EBIL_numNodes;  
    idData(6) = EBIL_numDOF; 
    idData(7) = crdTransf->getClassTag();
    int crdTransfDbTag  = crdTransf->getDbTag();
    if (crdTransfDbTag  == 0) {
        crdTransfDbTag = theChannel.getDbTag();
    if (crdTransfDbTag  != 0) 
        crdTransf->setDbTag(crdTransfDbTag);
    }
    idData(8) = crdTransfDbTag;
    
    // opserr <<"eleTag             :"<<this->getTag()      <<" "<<idData(0)<<"\n"; 
    // opserr <<"m_numEmbeddedPoints:"<<m_numEmbeddedPoints <<" "<<idData(1)<<"\n"; 
    // opserr <<"m_numSolidDOF      :"<<m_numSolidDOF       <<" "<<idData(2)<<"\n"; 
    // opserr <<"m_numSolidNodes    :"<<m_numSolidNodes     <<" "<<idData(3)<<"\n"; 
    // opserr <<"m_numBeamNodes     :"<<m_numBeamNodes      <<" "<<idData(4)<<"\n"; 
    // opserr <<"EBIL_numNodes      :"<<EBIL_numNodes       <<" "<<idData(5)<<"\n"; 
    // opserr <<"EBIL_numDOF        :"<<EBIL_numDOF         <<" "<<idData(6)<<"\n"; 
    // opserr <<"crdTransfClassTag  :"<<crdTransf->getClassTag()  <<"\n" ;
    // opserr <<"crdTransfDbTag     :"<<crdTransf->getDbTag()     <<"\n" ;
    

    res = theChannel.sendID(this->getDbTag(), commitTag, idData);
    if (res < 0) {
        opserr << "EmbeddedBeamInterfaceL::sendSelf -- could not send ID\n";
        return res;
    }

    

    
    int sizeofvector = (18*m_numEmbeddedPoints)+(6*m_numBeamNodes) + EBIL_numNodes + 1;
    static Vector data(sizeofvector);


    int index = 0;
    for (int i =0; i<8*m_numEmbeddedPoints;i++) {  
        data(index++) = solidNodeTags[i];
        // opserr << solidNodeTags[i] << " ";
    }

    for (int i =0; i<m_numEmbeddedPoints;i++) {  
        data(index++) = theBeamTags[i];
        // opserr << theBeamTags[i] << " ";
    }

    for (int i =0; i<2*m_numEmbeddedPoints;i++) {  
        data(index++) = beamNodeTags[i];
        // opserr << beamNodeTags[i] << " ";
    }

    for (int i =0; i<m_numEmbeddedPoints;i++) {  
        data(index++) = m_beam_rho[i];
        // opserr << m_beam_rho[i] << " ";
    }

    for (int i =0; i<m_numEmbeddedPoints;i++) {  
        data(index++) = m_beam_theta[i];
        // opserr << m_beam_theta[i] << " ";
    }

    for (int i =0; i<m_numEmbeddedPoints;i++) {  
        data(index++) = m_solid_xi[i];
        // opserr << m_solid_xi[i] << " ";
    }

    for (int i =0; i<m_numEmbeddedPoints;i++) {  
        data(index++) = m_solid_eta[i];
        // opserr << m_solid_eta[i] << " ";
    }

    for (int i =0; i<m_numEmbeddedPoints;i++) {  
        data(index++) = m_solid_zeta[i];
        // opserr << m_solid_zeta[i] << " ";
    }

    for (int i =0; i<m_numEmbeddedPoints;i++) {  
        data(index++) = m_area[i];
        // opserr << m_area[i] << " ";
    }

    for (int i =0; i<m_numEmbeddedPoints;i++) {  
        data(index++) = m_beamLength[i];
        // opserr << m_area[i] << " ";
    }

    for (int i =0; i<6*m_numBeamNodes;i++) {  
        data(index++) = m_Lambda[i];
        // opserr << m_Lambda[i] << " ";
    }
    
    for (int i =0; i<EBIL_numNodes;i++) {
        data(index++) = externalNodes[i];
        // opserr << externalNodes[i] << " ";
    }
    data(index++) = m_beam_radius;

    // opserr << "==================solidtags=====================\n";
    // for (int i =0; i<m_numEmbeddedPoints;i++) {
    //     data(index++) = theSolidTags[i];
    //     opserr << theSolidTags[i] << " ";
    // }


    // for (int i =0; i<EBIL_numNodes;i++) {  
    //     if (i<=m_numSolidNodes) {

    //         opserr << "solid: " << externalNodes[i]<< " "<< m_solidNodeMap[externalNodes[i]]<< " \n";
    //     } else {
    //         opserr << "beam:  " << externalNodes[i]<< " "<< m_beamNodeMap[externalNodes[i]] << " \n";
    //     }
    // }
    

    // int arrSize = sizeof(theNodes)/sizeof(theNodes[0]);
    // opserr << "theNodes size is equal to: " << arrSize<< "\n"; 







    res = theChannel.sendVector(this->getDbTag(), commitTag, data);
    if (res < 0) {
        opserr << "WARNING EmbeddedBeamInterfaceL::sendSelf() - " << this->getTag() << " failed to send Vector\n";
        return res;
    }

    // send the coordinate transformation
    if (crdTransf->sendSelf(commitTag, theChannel) < 0) {
        opserr << "EmbeddedBeamInterfaceL::sendSelf() - failed to send crdTranf\n";
        return -1;
    }


    // for (int i =0; i<8;i++) {  
    //     idData(index) = solidNodeTags[i];
    //     index++;
    //     opserr << solidNodeTags[i] << " ";
    // }

    // for (int i =0; i<m_numEmbeddedPoints;i++) {  
    //     data(index+i) = theBeamTags[i];
    //     index++;
    // }

    // for (int i =0; i<2*m_numEmbeddedPoints;i++) {  
    //     data(index+i) = beamNodeTags[i];
    //     index++;
    // }

    // for (int i =0; i<m_numEmbeddedPoints;i++) {  
    //     data(index+i) = m_beam_rho[i];
    //     index++;
    // }

    // for (int i =0; i<m_numEmbeddedPoints;i++) {  
    //     data(index+i) = m_beam_theta[i];
    //     index++;
    // }

    // for (int i =0; i<m_numEmbeddedPoints;i++) {  
    //     data(index+i) = m_solid_xi[i];
    //     index++;
    // }

    // for (int i =0; i<m_numEmbeddedPoints;i++) {  
    //     data(index+i) = m_solid_eta[i];
    //     index++;
    // }

    // for (int i =0; i<m_numEmbeddedPoints;i++) {  
    //     data(index+i) = m_solid_zeta[i];
    //     index++;
    // }

    // for (int i =0; i<m_numEmbeddedPoints;i++) {  
    //     data(index+i) = m_area[i];
    //     index++;
    // }

    // for (int i =0; i<m_numEmbeddedPoints;i++) {  
    //     data(index+i) = m_beamLength[i];
    //     index++;
    // }

    // for (int i =0; i<6*m_numBeamNodes;i++) {  
        // data(index+i) = m_Lambda[i];
    //     index++;
    // }

    
    // opserr << "yes\n";


    // // EmbeddedBeamInterfaceL then sends the tags of its nodes
    // res = theChannel.sendID(dataTag, commitTag, externalNodes);
    // if (res < 0) {
    // opserr << "WARNING EmbeddedBeamInterfaceL::sendSelf() - " << this->getTag() << " failed to send ID\n";
    //     return res;
    // }

    // int dbTag = crdTransf->getDbTag();
    // if (dbTag == 0) {
    //   dbTag = theChannel.getDbTag();
    //   if (dbTag != 0)
	// crdTransf->setDbTag(dbTag);
    // }

    // // Ask the CoordTransf to send itself
    // res = crdTransf->sendSelf(commitTag, theChannel);
    // if (res < 0) {
    //   opserr << "EmbeddedBeamInterfaceL::sendSelf -- could not send CoordTransf\n";
    //   return res;
    // }

    return 0;
}

int
EmbeddedBeamInterfaceL::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker
    &theBroker)
{
    int res = 0;
    static ID idData(9);
    res = theChannel.recvID(this->getDbTag(), commitTag, idData);
    if (res < 0) {
        opserr << "EmbeddedBeamInterfaceL::recvSelf -- could not recv ID\n";
        return res;
    }

    this->setTag((int)idData(0)); 
    m_numEmbeddedPoints = idData(1); 
    m_numSolidDOF       = idData(2); 
    m_numSolidNodes     = idData(3);
    m_numBeamNodes      = idData(4);  
    EBIL_numNodes       = idData(5);  
    EBIL_numDOF         = idData(6);  
    int crdTransfClassTag = idData(7);
    int crdTransfDbTag    = idData(8);    
    // opserr <<"eleTag             :"<<idData(0)          <<"\n";  
    // opserr <<"m_numEmbeddedPoints:"<<m_numEmbeddedPoints<<"\n" ; 
    // opserr <<"m_numSolidDOF      :"<<m_numSolidDOF      <<"\n" ; 
    // opserr <<"m_numSolidNodes    :"<<m_numSolidNodes    <<"\n" ; 
    // opserr <<"m_numBeamNodes     :"<<m_numBeamNodes     <<"\n" ; 
    // opserr <<"EBIL_numNodes      :"<<EBIL_numNodes      <<"\n" ; 
    // opserr <<"EBIL_numDOF        :"<<EBIL_numDOF        <<"\n" ; 
    // opserr <<"crdTransfClassTag  :"<<crdTransfClassTag  <<"\n" ;
    // opserr <<"crdTransfDbTag     :"<<crdTransfDbTag     <<"\n" ;



    int sizeofvector = (18*m_numEmbeddedPoints)+(6*m_numBeamNodes) + EBIL_numNodes+1;
    static Vector data(sizeofvector);

    res = theChannel.recvVector(this->getDbTag(), commitTag, data);
    if (res < 0) {
      opserr << "WARNING EmbeddedBeamInterfaceL::recvSelf() - failed to receive Vector\n";
      return res;
    }
    
    // for (int i =0; i<8*m_numEmbeddedPoints;i++) {  
    //     opserr << data(i) << " ";  
    // }
    // opserr << "yes it is working\n";



    if (!memoryallocated) {
        // theSolidTags  = new int[m_numEmbeddedPoints];
        solidNodeTags = new int[8 * m_numEmbeddedPoints];
        theBeamTags   = new int[m_numEmbeddedPoints];
        beamNodeTags  = new int[2 * m_numEmbeddedPoints];
        memoryallocated = true;
        m_beam_rho = m_beam_theta = m_solid_xi = m_solid_eta = m_solid_zeta = m_area = m_beamLength = Vector(m_numEmbeddedPoints);
        m_Lambda        = Vector(6 * m_numBeamNodes);
        m_solidInitDisp = Vector(m_numSolidNodes * 3);
        m_beamInitDisp  = Vector(m_numBeamNodes * 6);
        m_InterfaceForces = Vector(EBIL_numDOF);
        m_InterfaceStiffness = Matrix(EBIL_numDOF, EBIL_numDOF);
        mA = Matrix(3 * m_numSolidNodes, 6 * m_numBeamNodes);
        mB = Matrix(6 * m_numBeamNodes,  6 * m_numBeamNodes);
        
        
        // theNodes = new Node*[EBIL_numNodes];
    }


    int index = 0;
    for (int i =0; i<8*m_numEmbeddedPoints;i++) {  
        solidNodeTags[i] = (int)data(index++);
        // opserr << solidNodeTags[i] << " ";  
    }
    
    for (int i =0; i<m_numEmbeddedPoints;i++) {  
        theBeamTags[i] = (int)data(index++);
        // opserr << theBeamTags[i] << " ";
    }


    for (int i =0; i<2*m_numEmbeddedPoints;i++) {  
        beamNodeTags[i] = (int)data(index++);
        // opserr << beamNodeTags[i] << " ";
    }


    for (int i =0; i<m_numEmbeddedPoints;i++) {  
        m_beam_rho[i] = data(index++);
        // opserr << m_beam_rho[i] << " ";
    }

        for (int i =0; i<m_numEmbeddedPoints;i++) {  
        m_beam_theta[i]= data(index++);
        // opserr << m_beam_theta[i] << " ";
    }

    for (int i =0; i<m_numEmbeddedPoints;i++) {  
        m_solid_xi[i] = data(index++);
        // opserr << m_solid_xi[i] << " ";
    }

    for (int i =0; i<m_numEmbeddedPoints;i++) {  
        m_solid_eta[i] = data(index++);
        // opserr << m_solid_eta[i] << " ";
    }

    for (int i =0; i<m_numEmbeddedPoints;i++) {  
        m_solid_zeta[i] = data(index++);
        // opserr << m_solid_zeta[i] << " ";
    }

    for (int i =0; i<m_numEmbeddedPoints;i++) {  
        m_area[i] = data(index++);
        // opserr << m_area[i] << " ";
    }

    for (int i =0; i<m_numEmbeddedPoints;i++) {  
        m_beamLength[i] = data(index++);
        // opserr << m_beamLength[i] << " ";
    }

    
    for (int i =0; i<6*m_numBeamNodes;i++) {  
        m_Lambda[i] = data(index++);
        // opserr << m_Lambda[i] << " ";
    }


    // allocating externalNodes, m_solidNodeMap and m_beamNodeMap
    if (externalNodes.Size() != EBIL_numNodes) {
        // opserr << "size of externalNodes " << externalNodes.Size() <<"\n";
        externalNodes.resize(EBIL_numNodes);
        // opserr << "size of externalNodes " << externalNodes.Size() <<"\n";
    }
    for (int i =0; i<EBIL_numNodes;i++) {  
        externalNodes[i] = data(index++);
        // opserr << externalNodes[i] << " ";
        if (i<m_numSolidNodes) 
            m_solidNodeMap[externalNodes[i]] = i;
        if ((i>=m_numSolidNodes) && (i<(m_numSolidNodes+m_numBeamNodes)))
            m_beamNodeMap[externalNodes[i]] = i-m_numSolidNodes;  
    }
    
    m_beam_radius = data(index++);

    // opserr << "==================solidtags=====================\n";
    // for (int i =0; i<m_numEmbeddedPoints;i++) {
    //     theSolidTags[i] = (int)data(index++);
    //     opserr << theSolidTags[i] << " ";
    // }



    // create a new crdTransf object if one needed
    if (crdTransf == 0 || (crdTransf->getClassTag() != crdTransfClassTag)) {
        if (crdTransf != 0) {
	        delete crdTransf;
        }
        crdTransf = theBroker.getNewCrdTransf(crdTransfClassTag);
        if (crdTransf == 0) {
	        opserr << "EmbeddedBeamInterfaceL::recvSelf() - " <<
	        "failed to obtain a CrdTrans object with classTag" <<
	        crdTransfClassTag << endln;
	        return -2;	  
        }
    }
    crdTransf->setDbTag(crdTransfDbTag);
    res = crdTransf->recvSelf(commitTag, theChannel, theBroker);
    // invoke recvSelf on the crdTransf object
    if ( res < 0) {
        opserr << "EmbeddedBeamInterfaceL::sendSelf() - failed to recv CrdTranf\n";
        return -3;
    } 



    // for (int i =0; i<EBIL_numNodes;i++) {  
    //     if (i<=m_numSolidNodes) {
    //         opserr << "solid: " << externalNodes[i]<< " "<< m_solidNodeMap[externalNodes[i]]<< " \n";
    //     } else {
    //         opserr << "beam:  " << externalNodes[i]<< " "<< m_beamNodeMap[externalNodes[i]] << " \n";
    //     }
    // }


    // int arrSize = sizeof(theNodes)/sizeof(theNodes[0]);
    // opserr << "theNodes size is equal to: " << arrSize<< "\n"; 


    return 0;
}

int
EmbeddedBeamInterfaceL::displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes, int numMode)
{
    return 0;
}

void
EmbeddedBeamInterfaceL::Print(OPS_Stream &s, int flag)
{
    return;
}

Response*
EmbeddedBeamInterfaceL::setResponse(const char **argv, int argc,
    OPS_Stream &s)
{
    if (strcmp(argv[0], "force") == 0 || strcmp(argv[0], "globalForce") == 0)
    {
        return new ElementResponse(this, 1, Vector(3 * m_numEmbeddedPoints));
    }
    else if (strcmp(argv[0], "displacement") == 0 || strcmp(argv[0], "disp") == 0)
    {
        return new ElementResponse(this, 2, Vector(3 * m_numEmbeddedPoints));
    }
    else if (strcmp(argv[0], "beamForce") == 0 || strcmp(argv[0], "beamInteractionForce") == 0)
    {
        return new ElementResponse(this, 3, Vector(12 * (m_numBeamNodes - 1)));
    }
    else if (strcmp(argv[0], "solidForce") == 0 || strcmp(argv[0], "solidInteractionForce") == 0)
    {
        return new ElementResponse(this, 4, Vector(3 * m_numSolidNodes));
    }
    else if (strcmp(argv[0], "localForce") == 0 || strcmp(argv[0], "locForce") == 0)
    {
        return new ElementResponse(this, 5, Vector(3 * m_numEmbeddedPoints));
    }
    else if (strcmp(argv[0], "axialForce") == 0 || strcmp(argv[0], "locForceAxial") == 0)
    {
        return new ElementResponse(this, 6, Vector(3 * m_numEmbeddedPoints));
    }
    else if (strcmp(argv[0], "radialForce") == 0 || strcmp(argv[0], "locForceNormal") == 0)
    {
        return new ElementResponse(this, 7, Vector(3 * m_numEmbeddedPoints));
    }
    else if (strcmp(argv[0], "tangentialForce") == 0 || strcmp(argv[0], "locForceTangent") == 0)
    {
        return new ElementResponse(this, 8, Vector(3 * m_numEmbeddedPoints));
    }
    else if (strcmp(argv[0], "localDisplacement") == 0 || strcmp(argv[0], "locDisp") == 0)
    {
        return new ElementResponse(this, 9, Vector(3 * m_numEmbeddedPoints));
    }
    else if (strcmp(argv[0], "axialDisp") == 0 || strcmp(argv[0], "locDispAxial") == 0)
    {
        return new ElementResponse(this, 10, Vector(3 * m_numEmbeddedPoints));
    }
    else if (strcmp(argv[0], "radialDisp") == 0 || strcmp(argv[0], "locDispNormal") == 0)
    {
        return new ElementResponse(this, 11, Vector(3 * m_numEmbeddedPoints));
    }
    else if (strcmp(argv[0], "tangentialDisp") == 0 || strcmp(argv[0], "locDispTangent") == 0)
    {
        return new ElementResponse(this, 12, Vector(3 * m_numEmbeddedPoints));
    }
    else if (strcmp(argv[0], "beamLocalForce") == 0 || strcmp(argv[0], "beamInteractionLocalForce") == 0)
    {
        return new ElementResponse(this, 13, Vector(12 * (m_numBeamNodes - 1)));
    }
    else
    {
        opserr << "EmbeddedBeamInterfaceL Recorder, " << argv[0] << "is an unknown recorder request"
            << "  Element tag : " << this->getTag() << endln;
        return 0;
    }
}

int
EmbeddedBeamInterfaceL::getResponse(int responseID, Information &eleInformation)
{
    if (responseID == 1) // force
        return eleInformation.setVector(GetInteractionPtForce(0));
    else if (responseID == 2) // displacement
        return eleInformation.setVector(GetInteractionPtDisp(0));
    else if (responseID == 3) // beamForce
        return eleInformation.setVector(GetElementalForce(1));
    else if (responseID == 4) // solidForce
        return eleInformation.setVector(GetElementalForce(2));
    else if (responseID == 5) // localForce
        return eleInformation.setVector(GetInteractionPtForce(1));
    else if (responseID == 6) // axialForce
        return eleInformation.setVector(GetInteractionPtForce(2));
    else if (responseID == 7) // radialForce
        return eleInformation.setVector(GetInteractionPtForce(3));
    else if (responseID == 8) // tangentialForce
        return eleInformation.setVector(GetInteractionPtForce(4));
    else if (responseID == 9) // localDisplacement
        return eleInformation.setVector(GetInteractionPtDisp(1));
    else if (responseID == 10) // axialDisp
        return eleInformation.setVector(GetInteractionPtDisp(2));
    else if (responseID == 11) // radialDisp
        return eleInformation.setVector(GetInteractionPtDisp(3));
    else if (responseID == 12) // tangentialDisp
        return eleInformation.setVector(GetInteractionPtDisp(4));
    else if (responseID == 13) // beamLocalForce
        return eleInformation.setVector(GetElementalForce(3));
    else
    {
        opserr << "EmbeddedBeamInterfaceL, tag = " << this->getTag()
            << " -- unknown request" << endln;
        return -1;
    }
}

int
EmbeddedBeamInterfaceL::setParameter(const char **argv, int argc, Parameter &param)
{
    return 0;
}

int
EmbeddedBeamInterfaceL::updateParameter(int parameterID, Information &info)
{
    return 0;
}

void
EmbeddedBeamInterfaceL::setDomain(Domain *theDomain)
{
    for (int i = 0; i<EBIL_numNodes;i++) {
        if(!theNodesStatus) {
            delete [] theNodes;
            theNodes = new Node*[EBIL_numNodes];
            theNodesStatus = true;
            // opserr << "NEW MEMORY ALLOCATED to theNodes pointer";
        }
        theNodes[i] = theDomain->getNode(externalNodes(i));
        // opserr << theNodes[i]->getNumberDOF() << "\n";
    }


    for (int ii = 0; ii < m_numSolidNodes + 2 * m_numBeamNodes; ii++)
    {
        if (theNodes[ii] == 0)
        {
            opserr << "Could not find node " << externalNodes(ii) << "." << endln;
            return;
        }
        if (!((theNodes[ii]->getNumberDOF() == 3) || (theNodes[ii]->getNumberDOF() == 4)) && (ii < m_numSolidNodes))
        {   
            opserr << "Solid node " << externalNodes(ii) << " has to have 3 or 4 degrees of freedom." << endln;
            return;
        }
        if ((theNodes[ii]->getNumberDOF() != 6) && (ii > m_numSolidNodes - 1))
        {
            opserr << "Beam (or Lagrange) node " << externalNodes(ii) << " has to have 6 degrees of freedom." << endln;
            return;
        }
    }
    
    // initialize the transformation
    if (crdTransf->initialize(theDomain->getNode(beamNodeTags[0]), theDomain->getNode(beamNodeTags[1])))
    {
        opserr << "EmbeddedBeamInterfaceL::setDomain(): Error initializing coordinate transformation";
        return;
    }
    m_beam_length = crdTransf->getInitialLength();
    if (m_beam_length < 1.0e-12) {
        opserr << "FATAL ERROR EmbeddedBeamInterfaceL (tag: " << this->getTag() << ") : "
            << "Beam element has zero length." << endln;
        return;
    }
    Vector initXAxis(3);
    Vector initYAxis(3);
    Vector initZAxis(3);
    crdTransf->getLocalAxes(initXAxis, initYAxis, initZAxis);
    // fill mQa
    for (int i = 0; i < 3; i++) 
    {
        mQa(i, 0) = initXAxis(i);
        mQa(i, 1) = initYAxis(i);
        mQa(i, 2) = initZAxis(i);
    }
    // set mQb = mQa : beam column element requires zero initial twist
    mQc = mQb = mQa;    

    // calculate A and B
    mA.Zero();
    mB.Zero();
    for (int ii = 0; ii < m_numEmbeddedPoints; ii++)
    {
        updateShapeFuncs(m_solid_xi(ii), m_solid_eta(ii), m_solid_zeta(ii), m_beam_rho(ii), m_beamLength(ii));
        ComputeHf(mHf, m_beam_theta(ii));
        ComputeBphiAndBu(mBphi, mBu, m_beamLength(ii));

        Vector c2(3), c3(3);
        for (int jj = 0; jj < 3; jj++)
        {
            c2(jj) = mQc(jj, 1);
            c3(jj) = mQc(jj, 2);
        }
        Matrix Hb(3, 12);
        Hb  = mBu - (m_beam_radius*(cos(m_beam_theta(ii))*ComputeSkew(c2) + sin(m_beam_theta(ii))*ComputeSkew(c3))) * mBphi;


        //Matrix mM(3, 3);
        //mM(0, 1) = -m_beam_radius * sin(m_beam_theta(ii));
        //mM(0, 2) = m_beam_radius * cos(m_beam_theta(ii));
        //mM(1, 0) = m_beam_radius * sin(m_beam_theta(ii));
        //mM(2, 0) = -m_beam_radius * cos(m_beam_theta(ii));
        //Matrix Hb(3, 12);
        ////Hb  = mBu - (m_beam_radius*(cos(m_beam_theta(ii))*ComputeSkew(c2) + sin(m_beam_theta(ii))*ComputeSkew(c3))) * mBphi;
        //Hb = mBu - mQc * mM * mBphi;

        int beamNodeInA1 = m_beamNodeMap[beamNodeTags[2 * ii + 0]];
        int beamNodeInA2 = m_beamNodeMap[beamNodeTags[2 * ii + 1]];

        for (int jj = 0; jj < 8; jj++)
        {
            int solidNodeInA = m_solidNodeMap[solidNodeTags[8 * ii + jj]];

            double curNs = m_Ns(jj);
            
            for (int kk = 0; kk < 6; kk++)
            {
                mA(3 * solidNodeInA + 0, 6 * beamNodeInA1 + kk) += curNs * mHf(0, kk) * m_area(ii);
                mA(3 * solidNodeInA + 1, 6 * beamNodeInA1 + kk) += curNs * mHf(1, kk) * m_area(ii);
                mA(3 * solidNodeInA + 2, 6 * beamNodeInA1 + kk) += curNs * mHf(2, kk) * m_area(ii);

                mA(3 * solidNodeInA + 0, 6 * beamNodeInA2 + kk) += curNs * mHf(0, kk + 6) * m_area(ii);
                mA(3 * solidNodeInA + 1, 6 * beamNodeInA2 + kk) += curNs * mHf(1, kk + 6) * m_area(ii);
                mA(3 * solidNodeInA + 2, 6 * beamNodeInA2 + kk) += curNs * mHf(2, kk + 6) * m_area(ii);
            }
        }

        Matrix HbT  = Transpose(3, 12, Hb);
        Matrix temp = HbT * mHf * m_area(ii);

        for (int jj = 0; jj < 6; jj++)
        {
            for (int kk = 0; kk < 6; kk++)
            {
                mB(6 * beamNodeInA1 + jj, 6 * beamNodeInA1 + kk) += temp(jj    , kk    );
                mB(6 * beamNodeInA1 + jj, 6 * beamNodeInA2 + kk) += temp(jj    , kk + 6);
                mB(6 * beamNodeInA2 + jj, 6 * beamNodeInA1 + kk) += temp(jj + 6, kk    );
                mB(6 * beamNodeInA2 + jj, 6 * beamNodeInA2 + kk) += temp(jj + 6, kk + 6);
            }
        }

        /*FileStream tempFile("tempA.dat",APPEND);
        tempFile.setPrecision(15);
        tempFile << mA;
        FileStream tempFile2("tempB.dat",APPEND);
        tempFile2.setPrecision(15);
        tempFile2  << mB;*/

    }
    m_InterfaceStiffness.Zero();

    int II;
    for (int ii = 0; ii < 3 * m_numSolidNodes; ii++)
        for (int jj = 0; jj < 6 * m_numBeamNodes; jj++)
        {
            II = m_numSolidDOF * (ii / 3) + (ii % 3);
            m_InterfaceStiffness(II, m_numSolidDOF * m_numSolidNodes + 6 * m_numBeamNodes + jj) = mA(ii, jj);
            m_InterfaceStiffness(m_numSolidDOF * m_numSolidNodes + 6 * m_numBeamNodes + jj, II) = mA(ii, jj);
        }

    for (int ii = 0; ii < 6 * m_numBeamNodes; ii++)
        for (int jj = 0; jj < 6 * m_numBeamNodes; jj++)
        {
            m_InterfaceStiffness(m_numSolidDOF * m_numSolidNodes + ii, m_numSolidDOF * m_numSolidNodes + 6 * m_numBeamNodes + jj) = -mB(ii, jj);
            m_InterfaceStiffness(m_numSolidDOF * m_numSolidNodes + 6 * m_numBeamNodes + jj, m_numSolidDOF * m_numSolidNodes + ii) = -mB(ii, jj);
        }

    this->DomainComponent::setDomain(theDomain);
    return;
}

int
EmbeddedBeamInterfaceL::update(void)
{
    for (int ii = 0; ii < m_numBeamNodes; ii++)
    {
        Vector lambDisp = theNodes[m_numSolidNodes + m_numBeamNodes + ii]->getTrialDisp();
        for (int jj = 0; jj < 6; jj++)
            m_Lambda(6 * ii + jj) = lambDisp(jj);
    }

    return 0;
}

int
EmbeddedBeamInterfaceL::commitState(void)
{
    int retVal = 0;

    // call element commitState to do any base class stuff
    if ((retVal = this->Element::commitState()) != 0) {
        opserr << "EmbeddedBeamInterfaceL::commitState() - failed in base class";
    }

    return retVal;
}

int EmbeddedBeamInterfaceL::updateShapeFuncs(double xi, double eta, double zeta, double rho, double L)
{
    if ((xi < -1.0) || (xi > 1.0) || (eta < -1.0) || (eta > 1.0) || (zeta < -1.0) || (zeta > 1.0))
    {
        opserr << "Error in shape function." << endln;
        return -1;
    }

    if ((rho < -1.0) || (rho > 1.0))
    {
        opserr << "Error in shape function." << endln;
        return -1;
    }

    double rho2 = rho * rho;
    double rho3 = rho * rho2;

    m_Ns(0) = 0.125 * (1 - xi) * (1 - eta) * (1 - zeta);
    m_Ns(1) = 0.125 * (1 + xi) * (1 - eta) * (1 - zeta);
    m_Ns(2) = 0.125 * (1 + xi) * (1 + eta) * (1 - zeta);
    m_Ns(3) = 0.125 * (1 - xi) * (1 + eta) * (1 - zeta);
    m_Ns(4) = 0.125 * (1 - xi) * (1 - eta) * (1 + zeta);
    m_Ns(5) = 0.125 * (1 + xi) * (1 - eta) * (1 + zeta);
    m_Ns(6) = 0.125 * (1 + xi) * (1 + eta) * (1 + zeta);
    m_Ns(7) = 0.125 * (1 - xi) * (1 + eta) * (1 + zeta);

    m_Hb1 = 0.125 * (4.0 - 6.0 * rho + 2.0 * rho3);
    m_Hb3 = 0.125 * (4.0 + 6.0 * rho - 2.0 * rho3);
    m_Hb2 = 0.125 * L * (1.0 - rho - rho2 + rho3);
    m_Hb4 = 0.125 * L * (-1.0 - rho + rho2 + rho3);

    m_Nb1 = 0.5 * (1 - rho);
    m_Nb2 = 0.5 * (1 + rho); 

    m_dH1 = 1.5 * (-1.0 + rho2);
    m_dH3 = 1.5 * (1.0 - rho2);
    m_dH2 = 0.25 * L * (-1.0 - 2.0 * rho + 3.0 * rho2);
    m_dH4 = 0.25 * L * (-1.0 + 2.0 * rho + 3.0 * rho2);
    
    return 0;
}

Vector
EmbeddedBeamInterfaceL::CrossProduct(const Vector &V1, const Vector &V2)
{
    Vector V3(3);

    V3(0) = V1(1)*V2(2) - V1(2)*V2(1);
    V3(1) = V1(2)*V2(0) - V1(0)*V2(2);
    V3(2) = V1(0)*V2(1) - V1(1)*V2(0);

    return V3;
}

Matrix
EmbeddedBeamInterfaceL::Transpose(int dim1, int dim2, const Matrix &M)
{
    // copied from transpose function in Brick.cpp

    Matrix Mtran(dim2, dim1);

    for (int i = 0; i < dim1; i++)
        for (int j = 0; j < dim2; j++)
            Mtran(j, i) = M(i, j);

    return Mtran;
}


Matrix
EmbeddedBeamInterfaceL::ComputeSkew(Vector th)
{
    Matrix skew_th(3, 3);

    skew_th(0, 0) = 0.0;
    skew_th(0, 1) = -th(2);
    skew_th(0, 2) = th(1);
    skew_th(1, 0) = th(2);
    skew_th(1, 1) = 0.0;
    skew_th(1, 2) = -th(0);
    skew_th(2, 0) = -th(1);
    skew_th(2, 1) = th(0);
    skew_th(2, 2) = 0.0;

    return skew_th;
}

void
EmbeddedBeamInterfaceL::ComputeBphiAndBu(Matrix &Bphi, Matrix &Bu, double L)
{
    Matrix dummy1(3, 3);
    Matrix dummy2(3, 3);
    Matrix dummy3(3, 3);
    Matrix dummy4(3, 3);

    Bphi.Zero();
    Bu.Zero();
    
    // Compute Bphi(0:2, 3:5)
    dummy1.Zero();
    dummy2.Zero();
    dummy3.Zero();
    dummy4.Zero();

    // dummy1 = N1 * Qc*(E1 dyadic E1)
    dummy1(0, 0) = m_Nb1*mQc(0, 0);
    dummy1(1, 0) = m_Nb1*mQc(1, 0);
    dummy1(2, 0) = m_Nb1*mQc(2, 0);
    // dummy1 += dH2 * Qc*P1
    dummy1(0, 1) = m_dH2*mQc(0, 1)/L;  // dH2 * mQc(0:2,1:2)
    dummy1(1, 1) = m_dH2*mQc(1, 1)/L;
    dummy1(2, 1) = m_dH2*mQc(2, 1)/L;
    dummy1(0, 2) = m_dH2*mQc(0, 2)/L;
    dummy1(1, 2) = m_dH2*mQc(1, 2)/L;
    dummy1(2, 2) = m_dH2*mQc(2, 2)/L;
    // dummy2 = Qa^T
    dummy2 = Transpose(3, 3, mQa);
    // dummy3 = Qc * (N1*(E1 dyadic E1)+ dH2 * P1) * Qa^T
    dummy3 = dummy1*dummy2;
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            Bphi(i, 3 + j) = dummy3(i, j);

    // Reuse parts of dummy1 and dummy2 to calculate Bu(0:2,0:2)
    // dummy1 += H1 * Qc*P1
    dummy1(0, 1) = m_Hb1*mQc(0, 1);  // H1 * mQc(0:2,1:2)
    dummy1(1, 1) = m_Hb1*mQc(1, 1);
    dummy1(2, 1) = m_Hb1*mQc(2, 1);
    dummy1(0, 2) = m_Hb1*mQc(0, 2);
    dummy1(1, 2) = m_Hb1*mQc(1, 2);
    dummy1(2, 2) = m_Hb1*mQc(2, 2);
    // dummy3 = Qc * (N1*(E1 dyadic E1)+ H1 * P1) * Qa^T
    dummy3 = dummy1*dummy2;
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            Bu(i, j) = dummy3(i, j);

    // Reuse dummy2 and Compute Bphi(0:2, 0:2) and Bu(0:2, 3:5)
    dummy1.Zero();
    dummy3.Zero();
    // dummy1 = Qc*E^R*P1 (E^R is the skew symmetric meatrix for E1 cross product => E1 x a = [E^R].a)
    dummy1(0, 1) =  mQc(0, 2);
    dummy1(0, 2) = -mQc(0, 1);
    dummy1(1, 1) =  mQc(1, 2);
    dummy1(1, 2) = -mQc(1, 1);
    dummy1(2, 1) =  mQc(2, 2);
    dummy1(2, 2) = -mQc(2, 1);
    // dummy3 = Qc*E^R*P1*Qa^T
    dummy3 = dummy1*dummy2;
    // Compute Bphi(0:2, 0:2)
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            Bphi(i, j) = m_dH1 / L * dummy3(i, j);
    // Compute Bu(0:2, 3:5)
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            Bu(i, 3 + j) = -m_Hb2 * dummy3(i, j);

    // Reuse dummy1 and Compute Bphi(0:2, 6:8) and Bu(0:2, 9:11)
    dummy2.Zero();
    dummy3.Zero();

    // dummy2 = Qb^T
    dummy2 = Transpose(3, 3, mQb);
    // dummy3 = Qc*E^R*P1*Qb^T
    dummy3 = dummy1*dummy2;
    // Compute Bphi(0:2, 6 : 8)
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            Bphi(i, 6 + j) = m_dH3 / L * dummy3(i, j);
    // Compute Bu(0:2, 9:11)
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            Bu(i, 9 + j) = -m_Hb4 * dummy3(i, j);

    // Reuse dummy2 and Compute Bphi(0:2, 9:11)
    dummy1.Zero();
    dummy3.Zero();
    // dummy1 = N2 * Qc*(E1 dyadic E1)
    dummy1(0, 0) = m_Nb2*mQc(0, 0);  // N2 * mQc(0:2,0)
    dummy1(1, 0) = m_Nb2*mQc(1, 0);
    dummy1(2, 0) = m_Nb2*mQc(2, 0);
    // dummy1 += dH4 * Qc*P1
    dummy1(0, 1) = m_dH4*mQc(0, 1)/L;     // dH4 * mQc(0:2,1:2)
    dummy1(1, 1) = m_dH4*mQc(1, 1)/L;
    dummy1(2, 1) = m_dH4*mQc(2, 1)/L;
    dummy1(0, 2) = m_dH4*mQc(0, 2)/L;
    dummy1(1, 2) = m_dH4*mQc(1, 2)/L;
    dummy1(2, 2) = m_dH4*mQc(2, 2)/L;
    // dummy3 = Qc * (N2*(E1 dyadic E1)+ dH4 * P1) * Qb^T
    dummy3 = dummy1*dummy2;
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) 
            Bphi(i, 9 + j) = dummy3(i, j);

    // Reuse parts of dummy1 and dummy2 to calculate Bu(0:2,6:8)
    // dummy1 += dH4 * Qc*P1
    dummy1(0, 1) = m_Hb3*mQc(0, 1);     // H3 * mQc(0:2,1:2)
    dummy1(1, 1) = m_Hb3*mQc(1, 1);
    dummy1(2, 1) = m_Hb3*mQc(2, 1);
    dummy1(0, 2) = m_Hb3*mQc(0, 2);
    dummy1(1, 2) = m_Hb3*mQc(1, 2);
    dummy1(2, 2) = m_Hb3*mQc(2, 2);
    // dummy3 = Qc * (N2*(E1 dyadic E1)+ H3 * P1) * Qb^T
    dummy3 = dummy1*dummy2;
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            Bu(i, 6 + j) = dummy3(i, j);

    //// new way to calculate Bu
    //Bu.Zero();
    //dummy1.Zero();
    //dummy2.Zero();
    //dummy3.Zero();
    //dummy4.Zero();

    //dummy2 = Transpose(3, 3, mQa);
    //dummy3 = Transpose(3, 3, mQb);

    //dummy1(0, 0) = m_Nb1;
    //dummy1(1, 1) = m_Hb1;
    //dummy1(2, 2) = m_Hb1;

    //dummy4 = mQc * dummy1 * dummy2;

    //for (int i = 0; i < 3; i++)
    //    for (int j = 0; j < 3; j++)
    //        Bu(i, j) = dummy4(i, j);

    //dummy1.Zero();
    //dummy1(1, 2) = m_Hb2;
    //dummy1(2, 1) = -m_Hb2;

    //dummy4 = mQc * dummy1 * dummy2;

    //for (int i = 0; i < 3; i++)
    //    for (int j = 0; j < 3; j++)
    //        Bu(i, j + 3) = dummy4(i, j);


    //dummy1.Zero();
    //dummy1(0, 0) = m_Nb2;
    //dummy1(1, 1) = m_Hb3;
    //dummy1(2, 2) = m_Hb3;

    //dummy4 = mQc * dummy1 * dummy3;

    //for (int i = 0; i < 3; i++)
    //    for (int j = 0; j < 3; j++)
    //        Bu(i, j + 6) = dummy4(i, j);

    //dummy1.Zero();
    //dummy1(1, 2) = m_Hb4;
    //dummy1(2, 1) = -m_Hb4;

    //dummy4 = mQc * dummy1 * dummy3;

    //for (int i = 0; i < 3; i++)
    //    for (int j = 0; j < 3; j++)
    //        Bu(i, j + 9) = dummy4(i, j);

    //// new way to calculate Bphi
    //Bphi.Zero();
    //dummy1.Zero();
    //dummy2.Zero();
    //dummy3.Zero();
    //dummy4.Zero();

    //dummy2 = Transpose(3, 3, mQa);
    //dummy3 = Transpose(3, 3, mQb);

    //dummy1(1, 2) = -m_dH1 / L;
    //dummy1(2, 1) =  m_dH1 / L;

    //dummy4 = dummy1 * dummy2;

    //for (int i = 0; i < 3; i++)
    //    for (int j = 0; j < 3; j++)
    //        Bphi(i, j) = dummy4(i, j);

    //dummy1.Zero();
    //dummy1(0, 0) = m_Nb1;
    //dummy1(1, 1) = m_dH2 / L;
    //dummy1(2, 2) = m_dH2 / L;

    //dummy4 = dummy1 * dummy2;

    //for (int i = 0; i < 3; i++)
    //    for (int j = 0; j < 3; j++)
    //        Bphi(i, j + 3) = dummy4(i, j);


    //dummy1.Zero();
    //dummy1(1, 2) = -m_dH3 / L;
    //dummy1(2, 1) =  m_dH3 / L;

    //dummy4 = dummy1 * dummy3;

    //for (int i = 0; i < 3; i++)
    //    for (int j = 0; j < 3; j++)
    //        Bphi(i, j + 6) = dummy4(i, j);

    //dummy1.Zero();
    //dummy1(0, 0) = m_Nb2;
    //dummy1(1, 1) = m_dH4 / L;
    //dummy1(2, 2) = m_dH4 / L;

    //dummy4 = dummy1 * dummy3;

    //for (int i = 0; i < 3; i++)
    //    for (int j = 0; j < 3; j++)
    //        Bphi(i, j + 9) = dummy4(i, j);

    return;
}

void EmbeddedBeamInterfaceL::ComputeHf(Matrix & Hf, double theta)
{
    Hf.Zero();

    double oneOver2PiR = 0.5 / m_Pi / m_beam_radius;
    double oneOver2PiR2 = oneOver2PiR / m_beam_radius;
    for (int ii = 0; ii < 3; ii++)
    {
        for (int jj = 0; jj < 3; jj++)
        {
            Hf(ii, jj)     = oneOver2PiR * m_Nb1 * mQa(jj, ii);
            Hf(ii, jj + 6) = oneOver2PiR * m_Nb2 * mQb(jj, ii);
        }
        Hf(0, ii + 3)  = 2.0 * oneOver2PiR2 * m_Nb1 * (mQa(ii, 1) * sin(theta) - mQa(ii, 2) * cos(theta));
        Hf(1, ii + 3)  = -oneOver2PiR2 * mQa(ii, 0) * m_Nb1 * sin(theta);
        Hf(2, ii + 3)  =  oneOver2PiR2 * mQa(ii, 0) * m_Nb1 * cos(theta);
        Hf(0, ii + 9)  = 2.0 * oneOver2PiR2 * m_Nb2 * (mQb(ii, 1) * sin(theta) - mQb(ii, 2) * cos(theta));
        Hf(1, ii + 9)  = -oneOver2PiR2 * mQb(ii, 0) * m_Nb2 * sin(theta);
        Hf(2, ii + 9)  =  oneOver2PiR2 * mQb(ii, 0) * m_Nb2 * cos(theta);
    }
    
    Hf = mQc * Hf;
    return;
}

Vector EmbeddedBeamInterfaceL::GetInteractionPtDisp(int flag)
{
    Vector res(3 * m_numEmbeddedPoints);
    Vector c1(3), c2(3), c3(3);
    Vector bDisp(6 * m_numBeamNodes);
    Matrix Hb(3, 12);
    Vector ptDisp(3);

    // update local coordinate system
    for (int ii = 0; ii < 3; ii++)
    {
        c1(ii) = mQc(ii, 0);
        c2(ii) = mQc(ii, 1);
        c3(ii) = mQc(ii, 2);
    }

    for (int ii = 0; ii < m_numBeamNodes; ii++)
    {
        Vector bDispCur = theNodes[m_numSolidNodes + ii]->getTrialDisp();
        for (int jj = 0; jj < 6; jj++)
            bDisp(6 * ii + jj) = bDispCur(jj);
    }


    for (int ii = 0; ii < m_numEmbeddedPoints; ii++)
    {
        // update the interpolation functions for displacements and interaction forces
        updateShapeFuncs(m_solid_xi(ii), m_solid_eta(ii), m_solid_zeta(ii), m_beam_rho(ii), m_beamLength(ii));

        // update the matrices that define kinematics of the points on the beam surface
        ComputeBphiAndBu(mBphi, mBu, m_beamLength(ii));

        Hb = mBu - (m_beam_radius*(cos(m_beam_theta(ii))*ComputeSkew(c2) + sin(m_beam_theta(ii))*ComputeSkew(c3))) * mBphi;

        int beamNodeInA1 = m_beamNodeMap[beamNodeTags[2 * ii + 0]];
        int beamNodeInA2 = m_beamNodeMap[beamNodeTags[2 * ii + 1]];

        Vector thisBeamDisp(12);
        for (int jj = 0; jj < 6; jj++)
        {
            thisBeamDisp(jj    ) = bDisp(6 * beamNodeInA1 + jj);
            thisBeamDisp(jj + 6) = bDisp(6 * beamNodeInA2 + jj);
        }

        ptDisp = Hb * thisBeamDisp;

        Vector tempVec(3); double temp;
        switch (flag)
        {
        case 0:
            for (int jj = 0; jj < 3; jj++)
                res(3 * ii + jj) = ptDisp(jj);
            break;
        case 1:
            ptDisp = Transpose(3, 3, mQc) * ptDisp;
            for (int jj = 0; jj < 3; jj++)
                res(3 * ii + jj) = ptDisp(jj);
            break;
        case 2:
            temp = ptDisp(0)*c1(0) + ptDisp(1)*c1(1) + ptDisp(2)*c1(2);
            for (int jj = 0; jj < 3; jj++)
                res(3 * ii + jj) = temp*c1(jj);
            break;
        case 3:
            tempVec = cos(m_beam_theta(ii)) * c2 + sin(m_beam_theta(ii)) * c3;
            temp = ptDisp(0)*tempVec(0) + ptDisp(1)*tempVec(1) + ptDisp(2)*tempVec(2);
            for (int jj = 0; jj < 3; jj++)
                res(3 * ii + jj) = temp*tempVec(jj);
            break;
        case 4:
            tempVec = cos(m_beam_theta(ii)) * c3 - sin(m_beam_theta(ii)) * c2;
            temp = ptDisp(0)*tempVec(0) + ptDisp(1)*tempVec(1) + ptDisp(2)*tempVec(2);
            for (int jj = 0; jj < 3; jj++)
                res(3 * ii + jj) = temp*tempVec(jj);
            break;
        }
    }

    return res;
}

Vector EmbeddedBeamInterfaceL::GetInteractionPtForce(int flag)
{
    Vector res(3 * m_numEmbeddedPoints);
    Vector ptForces(3);
    Matrix Hf(3, 12);
    Vector c1(3), c2(3), c3(3);

    // update local coordinate system
    for (int ii = 0; ii < 3; ii++)
    {
        c1(ii) = mQc(ii, 0);
        c2(ii) = mQc(ii, 1);
        c3(ii) = mQc(ii, 2);
    }

    for (int ii = 0; ii < m_numEmbeddedPoints; ii++)
    {
        // update the interpolation functions for displacements and interaction forces
        updateShapeFuncs(m_solid_xi(ii), m_solid_eta(ii), m_solid_zeta(ii), m_beam_rho(ii), m_beamLength(ii));

        // update the matrices that define kinematics of the points on the beam surface
        ComputeHf(Hf, m_beam_theta(ii));
        
        int beamNodeInA1 = m_beamNodeMap[beamNodeTags[2 * ii + 0]];
        int beamNodeInA2 = m_beamNodeMap[beamNodeTags[2 * ii + 1]];

        Vector thisBeamForce(12);
        for (int jj = 0; jj < 6; jj++)
        {
            thisBeamForce(jj    ) = m_Lambda(6 * beamNodeInA1 + jj);
            thisBeamForce(jj + 6) = m_Lambda(6 * beamNodeInA2 + jj);
        }

        ptForces = Hf * thisBeamForce * m_area(ii);

        Vector tempVec(3); double temp;
        switch (flag)
        {
        case 0:
            for (int jj = 0; jj < 3; jj++)
                res(3 * ii + jj) = ptForces(jj);
            break;
        case 1:
            ptForces = Transpose(3, 3, mQc) * ptForces;
            for (int jj = 0; jj < 3; jj++)
                res(3 * ii + jj) = ptForces(jj);
            break;
        case 2:
            temp = ptForces(0)*c1(0) + ptForces(1)*c1(1) + ptForces(2)*c1(2);
            for (int jj = 0; jj < 3; jj++)
                res(3 * ii + jj) = temp*c1(jj);
            break;
        case 3:
            tempVec = cos(m_beam_theta(ii)) * c2 + sin(m_beam_theta(ii)) * c3;
            temp = ptForces(0)*tempVec(0) + ptForces(1)*tempVec(1) + ptForces(2)*tempVec(2);
            for (int jj = 0; jj < 3; jj++)
                res(3 * ii + jj) = temp*tempVec(jj);
            break;
        case 4:
            tempVec = cos(m_beam_theta(ii)) * c3 - sin(m_beam_theta(ii)) * c2;
            temp = ptForces(0)*tempVec(0) + ptForces(1)*tempVec(1) + ptForces(2)*tempVec(2);
            for (int jj = 0; jj < 3; jj++)
                res(3 * ii + jj) = temp*tempVec(jj);
            break;
        }
    }

    return res;
}

Vector EmbeddedBeamInterfaceL::GetElementalForce(int flag)
{
    Vector res;

    if (flag == 1)
    {
        res.resize(12 * (m_numBeamNodes - 1));
        Matrix eleB(12, 12);
        Matrix Hf(3, 12), Hu(3,12);
        Vector c2(3), c3(3);

        // update local coordinate system
        for (int ii = 0; ii < 3; ii++)
        {
            c2(ii) = mQc(ii, 1);
            c3(ii) = mQc(ii, 2);
        }

        int curBeamTag   = theBeamTags[0];
        int curBeamCount = 0;
        for (int ii = 0; ii < m_numEmbeddedPoints; ii++)
        {
            if (theBeamTags[ii] != curBeamTag)
            {
                curBeamTag = theBeamTags[ii];
                Vector thisBeamLambda(12);

                int beamNodeInA1 = m_beamNodeMap[beamNodeTags[2 * (ii - 1) + 0]];
                int beamNodeInA2 = m_beamNodeMap[beamNodeTags[2 * (ii - 1) + 1]];

                for (int jj = 0; jj < 6; jj++)
                {
                    thisBeamLambda(jj    ) = m_Lambda(6 * beamNodeInA1 + jj);
                    thisBeamLambda(jj + 6) = m_Lambda(6 * beamNodeInA2 + jj);
                }
                
                Vector thisBeamForce = -1.0  * eleB * thisBeamLambda;
                
                for (int jj = 0; jj < 6; jj++)
                {
                    res(12 * curBeamCount + jj + 0) = thisBeamForce(jj    );
                    res(12 * curBeamCount + jj + 6) = thisBeamForce(jj + 6);
                }

                eleB.Zero();
                curBeamCount++;
            }

            updateShapeFuncs(m_solid_xi(ii), m_solid_eta(ii), m_solid_zeta(ii), m_beam_rho(ii), m_beamLength(ii));
            ComputeHf(Hf, m_beam_theta(ii));
            ComputeBphiAndBu(mBphi, mBu, m_beamLength(ii));

            Hu = mBu - (m_beam_radius*(cos(m_beam_theta(ii))*ComputeSkew(c2) + sin(m_beam_theta(ii))*ComputeSkew(c3))) * mBphi;
            eleB += Transpose(3, 12, Hu) * Hf * m_area(ii);
        }
        Vector thisBeamLambda(12);

        int beamNodeInA1 = m_beamNodeMap[beamNodeTags[2 * (m_numEmbeddedPoints - 1) + 0]];
        int beamNodeInA2 = m_beamNodeMap[beamNodeTags[2 * (m_numEmbeddedPoints - 1) + 1]];

        for (int jj = 0; jj < 6; jj++)
        {
            thisBeamLambda(jj    ) = m_Lambda(6 * beamNodeInA1 + jj);
            thisBeamLambda(jj + 6) = m_Lambda(6 * beamNodeInA2 + jj);
        }

        Vector thisBeamForce = -1.0  * eleB * thisBeamLambda;

        for (int jj = 0; jj < 6; jj++)
        {
            res(12 * curBeamCount + jj + 0) = thisBeamForce(jj);
            res(12 * curBeamCount + jj + 6) = thisBeamForce(jj + 6);
        }
    }
    else if (flag == 2)
    {
        res = mA * m_Lambda;
    }
    else if (flag == 3)
    {
        res.resize(12 * (m_numBeamNodes - 1));
        Matrix eleB(12, 12);
        Matrix Hf(3, 12), Hu(3, 12);
        Vector c2(3), c3(3);

        // update local coordinate system
        for (int ii = 0; ii < 3; ii++)
        {
            c2(ii) = mQc(ii, 1);
            c3(ii) = mQc(ii, 2);
        }

        int curBeamTag = theBeamTags[0];
        int curBeamCount = 0;
        for (int ii = 0; ii < m_numEmbeddedPoints; ii++)
        {
            if (theBeamTags[ii] != curBeamTag)
            {
                curBeamTag = theBeamTags[ii];
                Vector thisBeamLambda(12);

                int beamNodeInA1 = m_beamNodeMap[beamNodeTags[2 * (ii - 1) + 0]];
                int beamNodeInA2 = m_beamNodeMap[beamNodeTags[2 * (ii - 1) + 1]];

                for (int jj = 0; jj < 6; jj++)
                {
                    thisBeamLambda(jj) = m_Lambda(6 * beamNodeInA1 + jj);
                    thisBeamLambda(jj + 6) = m_Lambda(6 * beamNodeInA2 + jj);
                }

                Vector thisBeamForce = -1.0  * eleB * thisBeamLambda;

                for (int jj = 0; jj < 4; jj++)
                {
                    Vector temp1(3), temp2(3); Matrix Q(3, 3);
                    for (int kk = 0; kk < 3; kk++)
                        temp1(kk) = thisBeamForce(3 * jj + kk);
                    if (jj < 2)
                        Q = Transpose(3, 3, mQa);
                    else
                        Q = Transpose(3, 3, mQb);
                    temp2 = Q * temp1;
                    for (int kk = 0; kk < 3; kk++)
                        thisBeamForce(3 * jj + kk) = temp2(kk);
                }

                for (int jj = 0; jj < 6; jj++)
                {
                    res(12 * curBeamCount + jj + 0) = thisBeamForce(jj);
                    res(12 * curBeamCount + jj + 6) = thisBeamForce(jj + 6);
                }

                eleB.Zero();
                curBeamCount++;
            }

            updateShapeFuncs(m_solid_xi(ii), m_solid_eta(ii), m_solid_zeta(ii), m_beam_rho(ii), m_beamLength(ii));
            ComputeHf(Hf, m_beam_theta(ii));
            ComputeBphiAndBu(mBphi, mBu, m_beamLength(ii));

            Vector c2(3), c3(3);
            for (int jj = 0; jj < 3; jj++)
            {
                c2(jj) = mQc(jj, 1);
                c3(jj) = mQc(jj, 2);
            }

            Hu = mBu - (m_beam_radius*(cos(m_beam_theta(ii))*ComputeSkew(c2) + sin(m_beam_theta(ii))*ComputeSkew(c3))) * mBphi;
            eleB += Transpose(3, 12, Hu) * Hf * m_area(ii);
        }
        Vector thisBeamLambda(12);

        int beamNodeInA1 = m_beamNodeMap[beamNodeTags[2 * (m_numEmbeddedPoints - 1) + 0]];
        int beamNodeInA2 = m_beamNodeMap[beamNodeTags[2 * (m_numEmbeddedPoints - 1) + 1]];

        for (int jj = 0; jj < 6; jj++)
        {
            thisBeamLambda(jj) = m_Lambda(6 * beamNodeInA1 + jj);
            thisBeamLambda(jj + 6) = m_Lambda(6 * beamNodeInA2 + jj);
        }

        Vector thisBeamForce = -1.0  * eleB * thisBeamLambda;
        
        for (int jj = 0; jj < 4; jj++)
        {
            Vector temp1(3), temp2(3); Matrix Q(3, 3);
            for (int kk = 0; kk < 3; kk++)
                temp1(kk) = thisBeamForce(3 * jj + kk);
            if (jj < 2)
                Q = Transpose(3, 3, mQa);
            else
                Q = Transpose(3, 3, mQb);
            temp2 = Q * temp1;
            for (int kk = 0; kk < 3; kk++)
                thisBeamForce(3 * jj + kk) = temp2(kk);
        }

        for (int jj = 0; jj < 6; jj++)
        {
            res(12 * curBeamCount + jj + 0) = thisBeamForce(jj);
            res(12 * curBeamCount + jj + 6) = thisBeamForce(jj + 6);
        }
    }
    else
    {
        opserr << "Unknown request ..." << endln;
        return Vector();
    }
    return res;
}

