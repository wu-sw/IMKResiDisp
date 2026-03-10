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
** redistribution, and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */

//**********************************************************************
// Shaowei Wu -- 2025 -- IMKResiDisp
//**********************************************************************

#ifndef IMKResiDisp_h
#define IMKResiDisp_h

#include <UniaxialMaterial.h>

class IMKResiDisp : public UniaxialMaterial
{
public:
    IMKResiDisp(int tag, double Ke,
        double posUy_0, double posUcap_0, double posUu_0, double posFy_0, double posFcapFy_0, double posFresFy_0,
        double negUy_0, double negUcap_0, double negUu_0, double negFy_0, double negFcapFy_0, double negFresFy_0,
        double LAMBDA_S, double LAMBDA_C, double LAMBDA_A, double c_S, double c_C, double c_A, 
        double alpha, double Offset, double theta_F, double theta_U, double beta_F, double beta_U);
    IMKResiDisp();
    ~IMKResiDisp();
    const char *getClassType(void) const { return "IMKResiDisp"; };
    int setTrialStrain(double strain, double strainRate = 0.0);
    double  getStrain(void);
    double  getStress(void);
    double  getTangent(void);
    double  getInitialTangent(void);
    int commitState(void);
    int revertToLastCommit(void);
    int revertToStart(void);
    UniaxialMaterial *getCopy(void);
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
    void Print(OPS_Stream &s, int flag = 0);


protected:

private:
// 25 Fixed input material parameters
    double  Ke;
    double  posUp_0;
    double  posUpc_0;
    double  posUu_0;
    double  posFy_0;
    double  posFcapFy_0;
    double  posFresFy_0;
    double  negUp_0;
    double  negUpc_0;
    double  negUu_0;
    double  negFy_0;
    double  negFcapFy_0;
    double  negFresFy_0;
    double  LAMBDA_S;
    double  LAMBDA_C;
    double  LAMBDA_A;
    double  c_S;
    double  c_C;
    double  c_A;
    double  alpha;
    double  Offset;
    double  theta_F;
    double  theta_U;
    double  beta_F;
    double  beta_U;
// 14 Initial Variables
    double  posUy_0,        negUy_0;
    double  posUcap_0,      negUcap_0;
    double  posFcap_0,      negFcap_0;
    double  posKp_0,        negKp_0;
    double  posKpc_0,       negKpc_0;
    double  engRefS;
    double  engRefC;
    double  engRefA;
// History Variables
// 12 Positive U and F
    double  posUy,          cPosUy;
    double  posFy,          cPosFy;
    double  posUcap,        cPosUcap;
    double  posFcap,        cPosFcap;
    double  posUlocal,      cPosUlocal;
    double  posFlocal,      cPosFlocal;
    double  posUglobal,     cPosUglobal;
    double  posFglobal,     cPosFglobal;
    double  posUres,        cPosUres;
    double  posFres,        cPosFres;
    double  posKp,          cPosKp;
    double  posKpc,         cPosKpc;
    double  pos_Uinter_local,   cPos_Uinter_local; 
    double  pos_Uinter_global,  cPos_Uinter_global;
    double  beta_F_pos,  cbeta_F_pos;
    double  beta_U_pos,  cbeta_U_pos;
    double  theta_F_cyc_pos,    ctheta_F_cyc_pos;
    double  theta_U_cyc_pos,    ctheta_U_cyc_pos;
// 12 Negative U and F
    double  negUy,          cNegUy;
    double  negFy,          cNegFy;
    double  negUcap,        cNegUcap;
    double  negFcap,        cNegFcap;
    double  negUlocal,      cNegUlocal;
    double  negFlocal,      cNegFlocal;
    double  negUglobal,     cNegUglobal;
    double  negFglobal,     cNegFglobal;
    double  negUres,        cNegUres;
    double  negFres,        cNegFres;
    double  negKp,          cNegKp;
    double  negKpc,         cNegKpc;
    double  neg_Uinter_local,   cNeg_Uinter_local;  
    double  neg_Uinter_global,  cNeg_Uinter_global; 
    double  beta_F_neg,  cbeta_F_neg;
    double  beta_U_neg,  cbeta_U_neg;
    double  theta_F_cyc_neg,    ctheta_F_cyc_neg;
    double  theta_U_cyc_neg,    ctheta_U_cyc_neg;
// 2 Pinching
// 4 Offset
    double  AsymPos,        cAsymPos;
    double  AsymNeg,        cAsymNeg;
    double  pos_Offset,     cPos_Offset;
    double  neg_Offset,     cNeg_Offset;
// 3 State Variables
    double  Ui,             cUi;
    double  Fi,             cFi;
// 3 Stiffness
    double  Kreload,        cKreload, KgetTangent;
    double  posKunload,     cPosKunload;
    double  negKunload,     cNegKunload;
// 2 Energy
    double  engAcml,        cEngAcml;
    double  engDspt,        cEngDspt;
// 2 Flag
    bool    Failure_Flag,   cFailure_Flag;
    int     Branch,         cBranch;
    bool    pos_ExtreSmallDisp, cPos_ExtreSmallDisp;
    bool    neg_ExtreSmallDisp, cNeg_ExtreSmallDisp;
};

#endif