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

#include <math.h>
#include <IMKResiDisp.h>
#include <elementAPI.h>
#include <Vector.h>
#include <Channel.h>
#include <OPS_Globals.h>
#include <algorithm>

#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
using namespace std;

static int numIMKResiDispMaterials = 0;

void *
OPS_IMKResiDisp()
{
    if (numIMKResiDispMaterials == 0) {
        numIMKResiDispMaterials++;
        opserr << "Improved IMK model for Residual Displacement calculation - by Shaowei Wu\n";
    }

    // Pointer to a uniaxial material that will be returned
    UniaxialMaterial *theMaterial = 0;

    int    iData[1];
    double dData[25];
    int numInt = 1;

    if (OPS_GetIntInput(&numInt, iData) != 0) {
        opserr << "WARNING invalid uniaxialMaterial IMKResiDisp tag" << endln;
        return 0;
    }

    int numDouble = 25;


    if (OPS_GetDoubleInput(&numDouble, dData) != 0) {
        opserr << "Invalid Args want: uniaxialMaterial IMKResiDisp tag? Ke? ";
        opserr << "posUp_0? posUpc_0? posUu_0? posFy_0? posFcapFy_0? posFresFy_0? ";
        opserr << "negUp_0? negUpc_0? negUu_0? negFy_0? negFcapFy_0? negFresFy_0? ";
        opserr << "LamdaS? LamdaC? LamdaA? Cs? Cc? Ca? ";
        opserr << "alpha? Offset? thetaF? thetaU? betaF? betaU?";
        return 0;
    }


    // Parsing was successful, allocate the material
    theMaterial = new IMKResiDisp(iData[0], dData[0],
        dData[1], dData[2], dData[3], dData[4], dData[5], dData[6],
        dData[7], dData[8], dData[9], dData[10], dData[11], dData[12],
        dData[13], dData[14], dData[15], dData[16], dData[17], dData[18], 
        dData[19], dData[20],dData[21], dData[22], dData[23], dData[24]);

    if (theMaterial == 0) {
        opserr << "WARNING could not create uniaxialMaterial of type IMKResiDisp Material\n";
        return 0;
    }

    return theMaterial;
}

IMKResiDisp::IMKResiDisp(int tag, double p_Ke,
    double p_posUp_0, double p_posUpc_0, double p_posUu_0, double p_posFy_0, double p_posFcapFy_0, double p_posFresFy_0,
    double p_negUp_0, double p_negUpc_0, double p_negUu_0, double p_negFy_0, double p_negFcapFy_0, double p_negFresFy_0,
    double p_LAMBDA_S, double p_LAMBDA_C, double p_LAMBDA_A, double p_c_S, double p_c_C, double p_c_A, 
    double p_alpha, double p_Offset, double p_theta_F, double p_theta_U, double p_beta_F, double p_beta_U)
    : UniaxialMaterial(tag, 0), Ke(p_Ke),
    posUp_0(p_posUp_0), posUpc_0(p_posUpc_0), posUu_0(p_posUu_0), posFy_0(p_posFy_0), posFcapFy_0(p_posFcapFy_0), posFresFy_0(p_posFresFy_0),
    negUp_0(p_negUp_0), negUpc_0(p_negUpc_0), negUu_0(p_negUu_0), negFy_0(p_negFy_0), negFcapFy_0(p_negFcapFy_0), negFresFy_0(p_negFresFy_0),
    LAMBDA_S(p_LAMBDA_S), LAMBDA_C(p_LAMBDA_C), LAMBDA_A(p_LAMBDA_A), c_S(p_c_S), c_C(p_c_C), c_A(p_c_A), 
    alpha(p_alpha), Offset(p_Offset), theta_F(p_theta_F), theta_U(p_theta_U), beta_F(p_beta_F), beta_U(p_beta_U)
{
    // Make sure these are all positive 
    if (negUp_0 < 0)
       negUp_0 = -negUp_0;
    if (negUpc_0 < 0)
       negUpc_0 = -negUpc_0;
    if (negUu_0 < 0)
       negUu_0 = -negUu_0;
    if (negFy_0 < 0)
       negFy_0 = -negFy_0;
    
    this->revertToStart();
}

IMKResiDisp::IMKResiDisp()
    :UniaxialMaterial(0, 0), Ke(0),
    posUp_0(0), posUpc_0(0), posUu_0(0), posFy_0(0), posFcapFy_0(0), posFresFy_0(0),
    negUp_0(0), negUpc_0(0), negUu_0(0), negFy_0(0), negFcapFy_0(0), negFresFy_0(0),
    LAMBDA_S(0), LAMBDA_C(0), LAMBDA_A(0),  c_S(0), c_C(0), c_A(0), 
    alpha(0), Offset(0), theta_F(0), theta_U(0), beta_F(0), beta_U(0)
{
    this->revertToStart();
}

IMKResiDisp::~IMKResiDisp()
{
    // does nothing
}

int IMKResiDisp::setTrialStrain(double strain, double strainRate)   //main function
{
   //all variables to the last commit
    this->revertToLastCommit();

    //state determination algorithm: defines the current force and tangent stiffness
    const double Ui_1 = Ui;
    const double Fi_1 = Fi;
    Ui = strain; //set trial displacement
    const double dU = Ui - Ui_1;    // Incremental deformation at current step
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////  MAIN CODE //////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    if (Failure_Flag) {     // When failure has already occured
        Fi = 0;
    } else if (dU == 0) {   // When deformation doesn't change from the last
        Fi = Fi_1;
    } else {
        double betaS = 0, betaC = 0, betaK = 0, betaA = 0;
        bool FailS = false, FailC = false, FailK = false, FailA = false;
        const bool onBackbone = (Branch > 1 && Branch != 11 && Branch != 2 && Branch != 12);
    ///////////////////////////////////////////////////////////////////////////////////////////
    /////////////////// WHEN REVERSAL /////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////
        if ( (onBackbone && Fi_1 * dU < 0) || (onBackbone && Fi_1 == 0 && Ui_1 * dU <= 0) ) {
            if (dU < 0) {
                Branch = 1;
            } else {
                Branch = 11;
            }
    /////////////////////////// UPDATE PEAK POINTS ////////////////////////////////////////////
            if ( Fi_1 > 0 ) {
                theta_F_cyc_pos = theta_F; theta_U_cyc_pos = theta_U;

                posUlocal = Ui_1;           // UPDATE LOCAL
                posFlocal = Fi_1;
                pos_Uinter_local = posUlocal - posFlocal / posKunload;
                
                const double pos_Uinter_local_1 = posUlocal - posFlocal / Ke; //using Ke
                neg_Offset = - (posFy_0*Offset) * (((pos_Uinter_local_1 + posUy_0) / posUy_0 - 1) / 4); 
                neg_Offset = ((pos_Uinter_local_1 + posUy_0) / posUy_0 > 1 ) ? neg_Offset : 0.0;
                if (neg_Offset < beta_F * negFglobal) {
                    neg_Offset = beta_F * negFglobal;
                }
                if ( Ui_1 > posUglobal) {    // UPDATE GLOBAL
                    posUglobal = Ui_1;
                    posFglobal = Fi_1;
                    posKunload = Ke / (pow((posUglobal / posUy_0), alpha));  
                    pos_Uinter_global = posUglobal - posFglobal / posKunload;  

                    beta_U_pos = beta_F_pos = 1.0;
                    if (posUglobal / posUy_0 > 2.0) { 

                        beta_U_pos = beta_U;
                        beta_F_pos = beta_F;
                    }
                }
                if (pos_ExtreSmallDisp == 1) { 
                    pos_ExtreSmallDisp = 0;
                    posKunload = Ke / (pow((posUglobal / posUy_0), alpha));
                    pos_Uinter_local = posUlocal - posFlocal / posKunload;
                }
                if ( pos_Uinter_local < neg_Uinter_local ) { 
                    pos_ExtreSmallDisp = 1;
                    posKunload = 1.1* posFlocal / (posUlocal - neg_Uinter_local);
                    pos_Uinter_local = posUlocal - posFlocal / posKunload;

                    theta_F_cyc_pos = 0.5; theta_U_cyc_pos = 0.5;  
                }
            } else {
                theta_F_cyc_neg = theta_F; theta_U_cyc_neg = theta_U;

                negUlocal = Ui_1;           // UPDATE LOCAL
                negFlocal = Fi_1;
                neg_Uinter_local = negUlocal - negFlocal / negKunload;
                
                const double neg_Uinter_local_1 = negUlocal - negFlocal / Ke; //using Ke
                pos_Offset = (negFy_0*Offset) * (((neg_Uinter_local_1+(-negUy_0)) / (-negUy_0) - 1) / 4); 
                pos_Offset = (((neg_Uinter_local_1+(-negUy_0)) / (-negUy_0)) > 1 ) ? pos_Offset : 0.0;
                if (pos_Offset > beta_F * posFglobal) {
                    pos_Offset = beta_F * posFglobal;
                }
                if ( Ui_1 < negUglobal) {    // UPDATE GLOBAL
                    negUglobal = Ui_1;
                    negFglobal = Fi_1;
                    negKunload = Ke / (pow((negUglobal / (-negUy_0)), alpha)); 
                    neg_Uinter_global = negUglobal - negFglobal / negKunload;  

                    beta_U_neg = beta_F_neg = 1.0;
                    if (negUglobal/(-negUy_0) > 2.0) {

                        beta_U_neg = beta_U;
                        beta_F_neg = beta_F;
                    }
                }
                if (neg_ExtreSmallDisp == 1) { //
                    neg_ExtreSmallDisp = 0;
                    negKunload = Ke / (pow((negUglobal / (-negUy_0)), alpha));
                    neg_Uinter_local = negUlocal - negFlocal / negKunload;
                }
                if ( pos_Uinter_local < neg_Uinter_local ) { //
                    neg_ExtreSmallDisp = 1;
                    negKunload = 1.1* negFlocal / (negUlocal - pos_Uinter_local);
                    neg_Uinter_local = negUlocal - negFlocal / negKunload;

                    theta_F_cyc_neg = 0.5; theta_U_cyc_neg = 0.5;
                }
            }
        }
        ///////////////////////////////////////////////////////////////////////////////////////////
        /////////////////// Calculate the degree of asymmetry  /////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////
        if (pos_Uinter_global > 0) { //initial value of pos_Uinter_global is 0, this 'If statement' is to assure that if no displacement happened in pos direction, the pos_Offset value takes its initial value (Offset)
            AsymPos = (pos_Uinter_global + neg_Uinter_local) / pos_Uinter_global;  //degree of asymmetry
        }
        if (neg_Uinter_global < 0) { //initial value of neg_Uinter_global is 0, this 'If statement' is to assure that if no displacement happened in neg direction, the neg_Offset value takes its initial value (-Offset)
            AsymNeg = (neg_Uinter_global + pos_Uinter_local) / neg_Uinter_global;
        }
    ///////////////////  preparation for judging if passing zero force line or offset line  ///////////////
        if (Branch == 1)  {                 //quadrant 1 : posKunload
            Fi = Fi_1 + posKunload * dU;
        } else if (Branch == 11) {          //quadrant 3 : negKunload
            Fi = Fi_1 + negKunload * dU;
        } else if (Branch == 12) {          //quadrant 4 : posKunload
            Fi = Fi_1 + posKunload * dU;
        } else if (Branch == 2) {           //quadrant 2 : negKunload
            Fi = Fi_1 + negKunload * dU;
        }
    ///////////////////////////////////////////////////////////////////////////////////////////
    /////////////////// WHEN PASS ZERO FORCE LINE/////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////
        if (Branch == 1 && (Fi_1  * Fi <= 0.0) && dU < 0 && Fi_1 > 0) {
            Branch = 12;
            Kreload = posKunload;
        }
        if (Branch == 11 && (Fi_1  * Fi <= 0.0) && dU > 0 && Fi_1 < 0) {
            Branch = 2;
            Kreload = negKunload;
        }
        if (Branch == 12 && (Fi_1  * Fi <= 0.0) && dU > 0) {  //also define the path return to Branch 1
            Branch = 1;
        }
        if (Branch == 2 && (Fi_1  * Fi <= 0.0) && dU < 0) {  //also define the path return to Branch 1
            Branch = 11;
        }
    ///////////////////////////////////////////////////////////////////////////////////////////
    /////////////////// WHEN PASS OFFSET LINE/////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////
        if (Branch == 2 && ((Fi_1 - pos_Offset) * (Fi - pos_Offset) <= 0.0) && dU > 0) {   //modified by Wu : dU > 0
    /////////////////// UPDATE BACKBONE CURVE /////////////////////////////////////////////////
            const double Ei = max(0.0, engAcml - engDspt);
            betaS = pow((Ei / (engRefS - engAcml)), c_S);
            betaC = pow((Ei / (engRefC - engAcml)), c_C);
            betaA = pow((Ei / (engRefA - engAcml)), c_A);
            FailS = (betaS > 1);
            FailC = (betaC > 1);
            FailA = (betaA > 1);
            betaS = betaS < 0 ? 0 : (betaS > 1 ? 1 : betaS);
            betaC = betaC < 0 ? 0 : (betaC > 1 ? 1 : betaC);
            betaA = betaA < 0 ? 0 : (betaA > 1 ? 1 : betaA);
            engDspt = engAcml;
        // Positive
            double FcapProj = posFcap - posKpc * posUcap;
        // Yield Point
            posFy *= (1 - betaS);
            posKp *= (1 - betaS); // Post-Yield Stiffness
            FcapProj *= (1 - betaC);
            posUglobal *= (1 + betaA); // Accelerated Reloading Stiffness
            posUy = posFy / Ke;
        // Capping Point
            const double FyProj = posFy - posKp*posUy;
            posUcap = posKp <= posKpc ? 0 : (FcapProj - FyProj) / (posKp - posKpc);
            posFcap = FyProj + posKp*posUcap;
        // When a part of backbone is beneath the residual strength
        // Global Peak on the Updated Backbone
            if (posUglobal < posUy) {           // Elastic Branch
                posFglobal = Ke * posUglobal;
            }
            else if (posUglobal < posUcap) {    // Post-Yield Branch
                posFglobal = posFy + posKp * (posUglobal - posUy);
            }
            else {                              // Post-Capping Branch
                posFglobal = posFcap + posKpc * (posUglobal - posUcap);
            }
            if (posFglobal < posFres) {     // Residual Branch
                posFglobal = posFres;
            }
            posUres = (posFres - posFcap + posKpc * posUcap) / posKpc;
            ////////////////////////// RELOADING TARGET DETERMINATION /////////////////////////////////
            const double u0 = Ui_1 - ((Fi_1 - pos_Offset) / negKunload);   // the displacement of the OFFSET line
            if (AsymPos < 1.0) {  //small asym degree (opposite side)
                const double Kglobal = (beta_F_pos * posFglobal - pos_Offset) / (beta_U_pos * posUglobal - u0);  //Add beta_F & beta_U by Wu
                const double Klocal  = (posFlocal - pos_Offset) / (posUlocal - u0);
                if ( (u0 < posUlocal) && (posFlocal < beta_F_pos * posFglobal) && (Klocal > Kglobal) ) {
                    Branch = 5;
                    Kreload = Klocal;
                } else {
                    Branch = 4;
                    Kreload = Kglobal;
                }
            } else {            //large asym degree (same side)，after shifting into branch 2, it will immediately shift into branch 5 or 4
                const double Kglobal = (theta_F * posFglobal) / (theta_U * (posUglobal - u0));
                const double Klocal  = (theta_F_cyc_pos * posFlocal) / (theta_U_cyc_pos * (posUlocal - u0));
                if ( u0 < posUlocal && posFlocal < posFglobal && Klocal > Kglobal) {
                    Branch = 51;    //Branch 5-1
                    Kreload = Klocal;
                } else {
                    Branch = 41;
                    Kreload = Kglobal;
                }
            }
        }
        if (Branch == 12 && ((Fi_1 - neg_Offset) * (Fi - neg_Offset) <= 0.0) && dU < 0) {   //modified by Wu : dU < 0
            const double Ei = max(0.0, engAcml - engDspt);
            betaS = pow((Ei / (engRefS - engAcml)), c_S);
            betaC = pow((Ei / (engRefC - engAcml)), c_C);
            betaA = pow((Ei / (engRefA - engAcml)), c_A);
            FailS = (betaS > 1);
            FailC = (betaC > 1);
            FailA = (betaA > 1);
            betaS = betaS < 0 ? 0 : (betaS > 1 ? 1 : betaS);
            betaC = betaC < 0 ? 0 : (betaC > 1 ? 1 : betaC);
            betaA = betaA < 0 ? 0 : (betaA > 1 ? 1 : betaA);
            engDspt = engAcml;
    
            double FcapProj = negFcap - negKpc * negUcap;
        // Yield Point
            negFy *= (1 - betaS);
            negKp *= (1 - betaS); // Post-Yield Stiffness
            FcapProj *= (1 - betaC);
            negUglobal *= (1 + betaA); // Accelerated Reloading Stiffness
            negUy = negFy / Ke;
        // Capping Point
            const double FyProj = negFy - negKp * negUy;
            negUcap = negKp <= negKpc ? 0 : (FcapProj - FyProj) / (negKp - negKpc);
            negFcap = FyProj + negKp * negUcap;
        // When a part of backbone is beneath the residual strength
        // Global Peak on the Updated Backbone
            if (negUy < negUglobal) {           // Elastic Branch
                negFglobal = Ke * negUglobal;
            }
            else if (negUcap < negUglobal) {    // Post-Yield Branch
                negFglobal = negFy + negKp * (negUglobal - negUy);
            }
            else {                              // Post-Capping Branch
                negFglobal = negFcap + negKpc * (negUglobal - negUcap);
            }
            if (negFres < negFglobal) {     // Residual Branch
                negFglobal = negFres;
            }
            negUres = (negFres - negFcap + negKpc * negUcap) / negKpc;
            ////////////////////////// RELOADING TARGET DETERMINATION /////////////////////////////////
            const double u0 = Ui_1 - ((Fi_1 - neg_Offset) / posKunload);   //Add an 'Offset' by Wu
            if (AsymNeg < 1.0) {  //small asym degree (opposite side)
                const double Kglobal = (beta_F_neg * negFglobal - neg_Offset) / (beta_U_neg * negUglobal - u0);  //Add beta_F & beta_U by Wu
                const double Klocal  = (negFlocal - neg_Offset) / (negUlocal - u0);
                if ( (u0 > negUlocal) && (negFlocal > beta_F_neg * negFglobal) && (Klocal > Kglobal) ) {
                    Branch = 15;
                    Kreload = Klocal;
                } else {
                    Branch = 14;
                    Kreload = Kglobal;
                }
            } else {            //large asym degree (same side)，after shifting into branch 2, it will immediately shift into branch 15 or 14
                const double Kglobal = (theta_F * negFglobal) / (theta_U * (negUglobal - u0));
                const double Klocal  = (theta_F_cyc_neg * negFlocal) / (theta_U_cyc_neg * (negUlocal - u0));
                if ( u0 > negUlocal && negFlocal > negFglobal && Klocal > Kglobal) {
                    Branch = 151;    //Branch 15-1
                    Kreload = Klocal;
                } else {
                    Branch = 141;
                    Kreload = Kglobal;
                }
            }
        }
    ///////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////
    ////////////////// BRANCH SHIFT CHECK AND TANGENT STIFFNESS UPDATE ////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////
        // int exBranch = Branch;
        if (Branch == 0 && Ui > posUy) {            // Yield in Positive
            Branch = 7;
        } else if (Branch == 0 && Ui < negUy) {     // Yield in Negative
            Branch = 17;
        } else if (Branch == 1 && Ui > posUlocal && Fi_1 > pos_Offset) {    // Back to Reloading (Positive)
            if ( AsymPos < 1 ){   //small asym degree (opposite side)
                if ( Ui < beta_U_pos * posUglobal ){
                    Kreload = (beta_F_pos * posFglobal - posFlocal) / (beta_U_pos * posUglobal - posUlocal);
                    Branch = 4;
                } else {
                    Kreload = (posFglobal - posFlocal) / (posUglobal - posUlocal);
                    Branch = 6;
                }
            } else {              //large asym degree (same side)
                if ( Fi < theta_F * posFglobal ){
                    Kreload = (theta_F * posFglobal - posFlocal) / ((neg_Uinter_local+theta_U*(posUglobal-neg_Uinter_local)) - posUlocal);
                    Branch = 41;
                } else {
                    Kreload = (posFglobal - posFlocal) / (posUglobal - posUlocal);
                    Branch = 61;
                }
            }
        } else if (Branch == 11 && Ui < negUlocal && Fi_1 < neg_Offset) {    // Back to Reloading (Negative)
            if ( AsymNeg < 1 ){   //small asym degree (opposite side)
                if ( Ui > beta_U_neg * negUglobal ){
                    Kreload = (beta_F_neg * negFglobal - negFlocal) / (beta_U_neg * negUglobal - negUlocal);
                    Branch = 14;
                } else {
                    Kreload = (negFglobal - negFlocal) / (negUglobal - negUlocal);
                    Branch = 16;
                }
            } else {              //large asym degree (same side)
                if ( Fi > theta_F * negFglobal ){
                    Kreload = (theta_F * negFglobal - negFlocal) / ((pos_Uinter_local+theta_U*(negUglobal-pos_Uinter_local)) - negUlocal);
                    Branch = 141;
                } else {
                    Kreload = (negFglobal - negFlocal) / (negUglobal - negUlocal);
                    Branch = 161;
                }
            }
        }
    // Positive
        if (Branch == 5 && Ui > posUlocal) {   //small asym degree (opposite side)
            Kreload = (beta_F_pos * posFglobal - posFlocal) / (beta_U_pos * posUglobal - posUlocal);
            Branch = 4;
        }
        if (Branch == 4 && Ui > beta_U_pos * posUglobal) {  //small asym degree (opposite side)
            Kreload = (beta_F_pos * posFglobal - posFglobal) / (beta_U_pos * posUglobal - posUglobal);
            Branch = 6;
        }
        if (Branch == 41 && Ui > neg_Uinter_local+theta_U*(posUglobal-neg_Uinter_local) ) {
            Kreload = (theta_F * posFglobal - posFglobal) / (neg_Uinter_local+theta_U*(posUglobal-neg_Uinter_local) - posUglobal);
            Branch = 61;
        }
        if (Branch == 51 && Ui > neg_Uinter_local+theta_U_cyc_pos*(posUlocal-neg_Uinter_local) ) {
            Kreload = (theta_F_cyc_pos * posFlocal - posFlocal) / (neg_Uinter_local+theta_U_cyc_pos*(posUlocal-neg_Uinter_local) - posUlocal);
            Branch = 52;
        }
        if (Branch == 52 && Ui > posUlocal ) {
            Kreload = (posFlocal+theta_F*(posFglobal-posFlocal) - posFlocal) / (posUlocal+theta_U*(posUglobal-posUlocal) - posUlocal);
            Branch = 53;
        }
        if (Branch == 53 && Ui > posUlocal+theta_U*(posUglobal-posUlocal) ) {
            Kreload = (posFlocal+theta_F*(posFglobal-posFlocal) - posFglobal) / (posUlocal+theta_U*(posUglobal-posUlocal) - posUglobal);
            Branch = 54;
        }
        if ((Branch == 6 || Branch == 61 || Branch == 54) && Ui > posUglobal) {
            Branch = 7;
        }
        if (Branch == 7 && Ui > posUcap) {
            Branch = 8;
        }
        if (Branch == 8 && Ui > posUres) {
            Branch = 9;
        }
    // Negative
        if (Branch == 15 && Ui < negUlocal) {   //small asym degree (opposite side)
            Kreload = (beta_F_neg * negFglobal - negFlocal) / (beta_U_neg * negUglobal - negUlocal);
            Branch = 14;
        }
        if (Branch == 14 && Ui < beta_U_neg * negUglobal) {  //small asym degree (opposite side)
            Kreload = (beta_F_neg * negFglobal - negFglobal) / (beta_U_neg * negUglobal - negUglobal);
            Branch = 16;
        }
        if (Branch == 141 && Ui < pos_Uinter_local+theta_U*(negUglobal-pos_Uinter_local) ) {
            Kreload = (theta_F * negFglobal - negFglobal) / (pos_Uinter_local+theta_U*(negUglobal-pos_Uinter_local) - negUglobal);
            Branch = 161;
        }
        if (Branch == 151 && Ui < pos_Uinter_local+theta_U_cyc_neg*(negUlocal-pos_Uinter_local) ) {
            Kreload = (theta_F_cyc_neg * negFlocal - negFlocal) / (pos_Uinter_local+theta_U_cyc_neg*(negUlocal-pos_Uinter_local) - negUlocal);
            Branch = 152;
        }
        if (Branch == 152 && Ui < negUlocal ) {
            Kreload = (negFlocal+theta_F*(negFglobal-negFlocal) - negFlocal) / (negUlocal+theta_U*(negUglobal-negUlocal) - negUlocal);
            Branch = 153;
        }
        if (Branch == 153 && Ui < negUlocal+theta_U*(negUglobal-negUlocal) ) {
            Kreload = (negFlocal+theta_F*(negFglobal-negFlocal) - negFglobal) / (negUlocal+theta_U*(negUglobal-negUlocal) - negUglobal);
            Branch = 154;
        }
        if ((Branch == 16 || Branch == 161 || Branch == 154) && Ui < negUglobal) {
            Branch = 17;
        }
        if (Branch == 17 && Ui < negUcap) {
            Branch = 18;
        }
        if (Branch == 18 && Ui < negUres) {
            Branch = 19;
        }
    // Branch Change check
        // if (Branch!=exBranch) {
        //     std::cout << exBranch << " -> " << Branch << "\n";
        // }
    ///////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////// COMPUTE FORCE BASED ON BRANCH /////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////
        if (Branch == 0) {
            Fi = Ke * Ui;
        } else if (Branch == 1) {
            Fi = Fi_1 + posKunload * dU;
        } else if (Branch == 11) {
            Fi = Fi_1 + negKunload * dU;
    // Positive
        } else if (Branch == 2) {
            Fi = Fi_1 + Kreload * dU;
        } else if (Branch == 4) {
            Fi = beta_F_pos * posFglobal + Kreload * (Ui - beta_U_pos * posUglobal);
        } else if ((Branch == 5) || (Branch == 52)) {
            Fi = posFlocal + Kreload * (Ui - posUlocal);
        } else if ((Branch == 6) || (Branch == 61) || (Branch == 54)) {
            Fi = posFglobal + Kreload * (Ui - posUglobal);
        } else if (Branch == 41) {
            Fi = (theta_F * posFglobal) + Kreload * (Ui - (neg_Uinter_local+theta_U*(posUglobal-neg_Uinter_local)));
        } else if (Branch == 51) {
            Fi = (theta_F_cyc_pos * posFlocal) + Kreload * (Ui - (neg_Uinter_local+theta_U_cyc_pos*(posUlocal-neg_Uinter_local)));
        } else if (Branch == 53) {
            Fi = (posFlocal+theta_F*(posFglobal-posFlocal)) + Kreload * (Ui - (posUlocal+theta_U*(posUglobal-posUlocal)));
        } else if (Branch == 7) {
            Fi = posFcap + posKp * (Ui - posUcap);
        } else if (Branch == 8) {
            Fi = posFcap + posKpc * (Ui - posUcap);
        } else if (Branch == 9) {
            Fi = posFres;
    // Negative
        } else if (Branch == 12) {
            Fi = Fi_1 + Kreload * dU;
        } else if (Branch == 14) {
            Fi = beta_F_neg * negFglobal + Kreload * (Ui - beta_U_neg * negUglobal);
        } else if ((Branch == 15) || (Branch == 152)) {
            Fi = negFlocal + Kreload * (Ui - negUlocal);
        } else if ((Branch == 16) || (Branch == 161) || (Branch == 154)) {
            Fi = negFglobal + Kreload * (Ui - negUglobal);
        } else if (Branch == 141) {
            Fi = (theta_F * negFglobal) + Kreload * (Ui - (pos_Uinter_local+theta_U*(negUglobal-pos_Uinter_local)));
        } else if (Branch == 151) {
            Fi = (theta_F_cyc_neg * negFlocal) + Kreload * (Ui - (pos_Uinter_local+theta_U_cyc_neg*(negUlocal-pos_Uinter_local)));
        } else if (Branch == 153) {
            Fi = (negFlocal+theta_F*(negFglobal-negFlocal)) + Kreload * (Ui - (negUlocal+theta_U*(negUglobal-negUlocal)));
        } else if (Branch == 17) {
            Fi = negFcap + negKp * (Ui - negUcap);
        } else if (Branch == 18) {
            Fi = negFcap + negKpc * (Ui - negUcap);
        } else if (Branch == 19) {
            Fi = negFres;
        }
    ///////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////
    // CHECK FOR FAILURE
    ///////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////
        const bool FailPp = ( posFglobal == 0 );
        const bool FailPn = ( negFglobal == 0 );
        const bool FailDp = ( dU > 0 && Ui >=  posUu_0 );
        const bool FailDn = ( dU < 0 && Ui <= -negUu_0 );
        const bool FailRp = ( Branch == 9 && Fi <= 0);
        const bool FailRn = ( Branch == 19 && Fi >= 0);
        if (FailS || FailC || FailA || FailK || FailPp || FailPn || FailRp || FailRn || FailDp || FailDn) {
            Fi = 0;
            Failure_Flag = true;
        }

        engAcml += 0.5 * (Fi + Fi_1) * dU;   // Internal energy increment

        KgetTangent = (Fi - Fi_1) / dU;
    }
    if (KgetTangent == 0) {
        KgetTangent = 1e-6;
    }
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////// END OF MAIN CODE ///////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    return 0;
}

double IMKResiDisp::getStress(void)
{
    //cout << " getStress" << endln;
    return (Fi);
}

double IMKResiDisp::getTangent(void)
{
    //cout << " getTangent" << endln;
    return (KgetTangent);
}

double IMKResiDisp::getInitialTangent(void)
{
    //cout << " getInitialTangent" << endln;
    return (Ke);
}

double IMKResiDisp::getStrain(void)
{
    //cout << " getStrain" << endln;
    return (Ui);
}

int IMKResiDisp::commitState(void)
{
    //cout << " commitState" << endln;

    //commit trial  variables
// 12 Pos U and F
    cPosUy = posUy;
    cPosFy = posFy;
    cPosUcap = posUcap;
    cPosFcap = posFcap;
    cPosUlocal = posUlocal;
    cPosFlocal = posFlocal;
    cPosUglobal = posUglobal;
    cPosFglobal = posFglobal;
    cPosUres = posUres;
    cPosFres = posFres;
    cPosKp = posKp;
    cPosKpc = posKpc;
    cPos_Uinter_local  = pos_Uinter_local;
    cPos_Uinter_global = pos_Uinter_global;
    cbeta_F_pos = beta_F_pos;
    cbeta_U_pos = beta_U_pos;
    ctheta_F_cyc_pos = theta_F_cyc_pos;
    ctheta_U_cyc_pos = theta_U_cyc_pos;
// 12 Neg U and F
    cNegUy = negUy;
    cNegFy = negFy;
    cNegUcap = negUcap;
    cNegFcap = negFcap;
    cNegUlocal = negUlocal;
    cNegFlocal = negFlocal;
    cNegUglobal = negUglobal;
    cNegFglobal = negFglobal;
    cNegUres = negUres;
    cNegFres = negFres;
    cNegKp = negKp;
    cNegKpc = negKpc;
    cNeg_Uinter_local  = neg_Uinter_local;
    cNeg_Uinter_global = neg_Uinter_global;
    cbeta_F_neg = beta_F_neg;
    cbeta_U_neg = beta_U_neg;
    ctheta_F_cyc_neg = theta_F_cyc_neg;
    ctheta_U_cyc_neg = theta_U_cyc_neg;
// 2 Pinching
// 4 Offset
    cAsymPos = AsymPos;
    cAsymNeg = AsymNeg;
    cPos_Offset = pos_Offset;
    cNeg_Offset = neg_Offset;
// 3 State
    cUi = Ui;
    cFi = Fi;
// 2 Stiffness
    cKreload = Kreload;
    cPosKunload = posKunload;
    cNegKunload = negKunload;
// 2 Energy
    cEngAcml = engAcml;
    cEngDspt = engDspt;
// 2 Flag
    cFailure_Flag = Failure_Flag;
    cBranch = Branch;
    cPos_ExtreSmallDisp = pos_ExtreSmallDisp;
    cNeg_ExtreSmallDisp = neg_ExtreSmallDisp;
    return 0;
}

int IMKResiDisp::revertToLastCommit(void)
{
    //cout << " revertToLastCommit" << endln;
    //the opposite of commit trial history variables
// 12 Positive U and F
    posUy = cPosUy;
    posFy = cPosFy;
    posUcap = cPosUcap;
    posFcap = cPosFcap;
    posUlocal = cPosUlocal;
    posFlocal = cPosFlocal;
    posUglobal = cPosUglobal;
    posFglobal = cPosFglobal;
    posUres = cPosUres;
    posFres = cPosFres;
    posKp = cPosKp;
    posKpc = cPosKpc;
    pos_Uinter_local  = cPos_Uinter_local;
    pos_Uinter_global = cPos_Uinter_global;
    beta_F_pos = cbeta_F_pos;
    beta_U_pos = cbeta_U_pos;
    theta_F_cyc_pos = ctheta_F_cyc_pos;
    theta_U_cyc_pos = ctheta_U_cyc_pos;
// 12 Negative U and F
    negUy = cNegUy;
    negFy = cNegFy;
    negUcap = cNegUcap;
    negFcap = cNegFcap;
    negUlocal = cNegUlocal;
    negFlocal = cNegFlocal;
    negUglobal = cNegUglobal;
    negFglobal = cNegFglobal;
    negUres = cNegUres;
    negFres = cNegFres;
    negKp = cNegKp;
    negKpc = cNegKpc;
    neg_Uinter_local  = cNeg_Uinter_local;
    neg_Uinter_global = cNeg_Uinter_global;
    beta_F_neg = cbeta_F_neg;
    beta_U_neg = cbeta_U_neg;
    theta_F_cyc_neg = ctheta_F_cyc_neg;
    theta_U_cyc_neg = ctheta_U_cyc_neg;
// 2 Pinching
// 2 Offset
    AsymPos = cAsymPos;
    AsymNeg = cAsymNeg;
    pos_Offset = cPos_Offset;
    neg_Offset = cNeg_Offset;
// 3 State Variables
    Ui = cUi;
    Fi = cFi;
// 2 Stiffness
    Kreload = cKreload;
    posKunload = cPosKunload;
    negKunload = cNegKunload;
// 2 Energy
    engAcml = cEngAcml;
    engDspt = cEngDspt;
// 2 Flag
    Failure_Flag = cFailure_Flag;
    Branch = cBranch;
    pos_ExtreSmallDisp = cPos_ExtreSmallDisp;
    neg_ExtreSmallDisp = cNeg_ExtreSmallDisp;
    return 0;
}

int IMKResiDisp::revertToStart(void)
{
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////\\
    //////////////////////////////////////////////////////////////////// ONE TIME CALCULATIONS ////////////////////////////////////////////////////////////////////\\
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
// 14 Initial Values
    posUy_0 = posFy_0 / Ke;
    posUcap_0 = posUy_0 + posUp_0;
    posFcap_0 = posFcapFy_0*posFy_0;
    posKp_0 = (posFcap_0 - posFy_0) / posUp_0;
    posKpc_0 = posFcap_0 / posUpc_0;
    negUy_0 = negFy_0 / Ke;
    negUcap_0 = negUy_0 + negUp_0;
    negFcap_0 = negFcapFy_0*negFy_0;
    negKp_0 = (negFcap_0 - negFy_0) / negUp_0;
    negKpc_0 = negFcap_0 / negUpc_0;
    engRefS = LAMBDA_S * posFy_0;
    engRefC = LAMBDA_C * posFy_0;
    engRefA = LAMBDA_A * posFy_0;
// 12 Positive U and F
    posUy = cPosUy = posUy_0;
    posFy = cPosFy = posFy_0;
    posUcap = cPosUcap = posUcap_0;
    posFcap = cPosFcap = posFcap_0;
    posUlocal = cPosUlocal = posUy_0;
    posFlocal = cPosFlocal = posFy_0;
    posUglobal = cPosUglobal = posUy_0;
    posFglobal = cPosFglobal = posFy_0;
    posFres = cPosFres = posFy_0*posFresFy_0;
    posKp = cPosKp =  posKp_0;
    posKpc = cPosKpc = -posKpc_0;
    posUres = cPosUres = (posFres - posFcap) / posKpc + posUcap;
    pos_Uinter_local  = cPos_Uinter_local  = 0;
    pos_Uinter_global = cPos_Uinter_global = 0;
    beta_F_pos = cbeta_F_pos = 1.0; beta_U_pos = cbeta_U_pos = 1.0;
    theta_F_cyc_pos = ctheta_F_cyc_pos = theta_F; theta_U_cyc_pos = ctheta_U_cyc_pos = theta_U;
// 12 Negative U and F
    negUy = cNegUy = -negUy_0;
    negFy = cNegFy = -negFy_0;
    negUcap = cNegUcap = -negUcap_0;
    negFcap = cNegFcap = -negFcap_0;
    negUlocal = cNegUlocal = -negUy_0;
    negFlocal = cNegFlocal = -negFy_0;
    negUglobal = cNegUglobal = -negUy_0;
    negFglobal = cNegFglobal = -negFy_0;
    negFres = cNegFres = -negFy_0*negFresFy_0;
    negKp = cNegKp =  negKp_0;
    negKpc = cNegKpc = -negKpc_0;
    negUres = cNegUres = (negFres - negFcap) / negKpc + negUcap;
    neg_Uinter_local  = cNeg_Uinter_local  = 0;
    neg_Uinter_global = cNeg_Uinter_global = 0;
    beta_F_neg = cbeta_F_neg = 1.0; beta_U_neg = cbeta_U_neg = 1.0;
    theta_F_cyc_neg = ctheta_F_cyc_neg = theta_F; theta_U_cyc_neg = ctheta_U_cyc_neg = theta_U;
// 3 State Values
    Ui = cUi = 0;
    Fi = cFi = 0;
// 2 Stiffness
    Kreload = cKreload = Ke;
    posKunload = cPosKunload = Ke;
    negKunload = cNegKunload = Ke;
    KgetTangent = Ke;
// 2 Energy
    engAcml = cEngAcml = 0.0;
    engDspt = cEngDspt = 0.0;
// 2 Flag
    Failure_Flag = cFailure_Flag = false;
    Branch = cBranch = 0;
    pos_ExtreSmallDisp = cPos_ExtreSmallDisp = 0;
    neg_ExtreSmallDisp = cNeg_ExtreSmallDisp = 0;
// 2 Pinching
// 2 Offset
    AsymPos = cAsymPos = 0.0;
    AsymNeg = cAsymNeg = 0.0;
    pos_Offset = cPos_Offset = posFy_0*Offset;
    neg_Offset = cNeg_Offset = -posFy_0*Offset;
    return 0;
}

UniaxialMaterial *
IMKResiDisp::getCopy(void)
{
    IMKResiDisp *theCopy = new IMKResiDisp(this->getTag(), Ke,
        posUy_0, posUcap_0, posUu_0, posFy_0, posFcapFy_0, posFresFy_0,
        negUy_0, negUcap_0, negUu_0, negFy_0, negFcapFy_0, negFresFy_0,
        LAMBDA_S, LAMBDA_C, LAMBDA_A, c_S, c_C, c_A,  alpha, Offset,
        beta_F, beta_U, theta_F, theta_U);

    //cout << " getCopy" << endln;
// 12 Positive U and F
    theCopy->posUy = posUy;
    theCopy->posFy = posFy;
    theCopy->posUcap = posUcap;
    theCopy->posFcap = posFcap;
    theCopy->posUlocal = posUlocal;
    theCopy->posFlocal = posFlocal;
    theCopy->posUglobal = posUglobal;
    theCopy->posFglobal = posFglobal;
    theCopy->posUres = posUres;
    theCopy->posFres = posFres;
    theCopy->posKp = posKp;
    theCopy->posKpc = posKpc;
    theCopy->pos_Uinter_local  = pos_Uinter_local ;
    theCopy->pos_Uinter_global = pos_Uinter_global;
    theCopy->beta_F_pos = beta_F_pos;
    theCopy->beta_U_pos = beta_U_pos;
    theCopy->theta_F_cyc_pos   = theta_F_cyc_pos  ;
    theCopy->theta_U_cyc_pos   = theta_U_cyc_pos  ;
// 12 Negative U and F
    theCopy->negUy = negUy;
    theCopy->negFy = negFy;
    theCopy->negUcap = negUcap;
    theCopy->negFcap = negFcap;
    theCopy->negUlocal = negUlocal;
    theCopy->negFlocal = negFlocal;
    theCopy->negUglobal = negUglobal;
    theCopy->negFglobal = negFglobal;
    theCopy->negUres = negUres;
    theCopy->negFres = negFres;
    theCopy->negKp = negKp;
    theCopy->negKpc = negKpc;
    theCopy->neg_Uinter_local  = neg_Uinter_local ;
    theCopy->neg_Uinter_global = neg_Uinter_global;
    theCopy->beta_F_neg = beta_F_neg;
    theCopy->beta_U_neg = beta_U_neg;
    theCopy->theta_F_cyc_neg   = theta_F_cyc_neg  ;
    theCopy->theta_U_cyc_neg   = theta_U_cyc_neg  ;
// 2 Pinching
// 2 Offset
    theCopy->AsymPos = AsymPos;
    theCopy->AsymNeg = AsymNeg;
    theCopy->pos_Offset = pos_Offset;
    theCopy->neg_Offset = neg_Offset;
// 3 State Values
    theCopy->Ui = Ui;
    theCopy->Fi = Fi;
// 2 Stiffness
    theCopy->Kreload = Kreload;
    theCopy->posKunload = posKunload;
    theCopy->negKunload = negKunload;
// 2 Energy
    theCopy->engAcml = engAcml;
    theCopy->engDspt = engDspt;
// 2 Flag
    theCopy->Failure_Flag = Failure_Flag;
    theCopy->Branch = Branch;
    theCopy->pos_ExtreSmallDisp = pos_ExtreSmallDisp;
    theCopy->neg_ExtreSmallDisp = neg_ExtreSmallDisp;
// 12 Positive U and F
    theCopy->cPosUy = cPosUy;
    theCopy->cPosFy = cPosFy;
    theCopy->cPosUcap = cPosUcap;
    theCopy->cPosFcap = cPosFcap;
    theCopy->cPosUlocal = cPosUlocal;
    theCopy->cPosFlocal = cPosFlocal;
    theCopy->cPosUglobal= cPosUglobal;
    theCopy->cPosFglobal= cPosFglobal;
    theCopy->cPosUres = cPosUres;
    theCopy->cPosFres = cPosFres;
    theCopy->cPosKp = cPosKp;
    theCopy->cPosKpc = cPosKpc;
    theCopy->cPos_Uinter_local  = cPos_Uinter_local ;
    theCopy->cPos_Uinter_global = cPos_Uinter_global;
    theCopy->cbeta_F_pos = cbeta_F_pos;
    theCopy->cbeta_U_pos = cbeta_U_pos;
    theCopy->ctheta_F_cyc_pos   = ctheta_F_cyc_pos  ;
    theCopy->ctheta_U_cyc_pos   = ctheta_U_cyc_pos  ;
// 12 Negative U and F
    theCopy->cNegUy = cNegUy;
    theCopy->cNegFy = cNegFy;
    theCopy->cNegUcap = cNegUcap;
    theCopy->cNegFcap = cNegFcap;
    theCopy->cNegUglobal= cNegUglobal;
    theCopy->cNegFglobal= cNegFglobal;
    theCopy->cNegUlocal = cNegUlocal;
    theCopy->cNegFlocal = cNegFlocal;
    theCopy->cNegUres = cNegUres;
    theCopy->cNegFres = cNegFres;
    theCopy->cNegKp = cNegKp;
    theCopy->cNegKpc = cNegKpc;
    theCopy->cNeg_Uinter_local  = cNeg_Uinter_local ;
    theCopy->cNeg_Uinter_global = cNeg_Uinter_global;
    theCopy->cbeta_F_neg = cbeta_F_neg;
    theCopy->cbeta_U_neg = cbeta_U_neg;
    theCopy->ctheta_F_cyc_neg   = ctheta_F_cyc_neg  ;
    theCopy->ctheta_U_cyc_neg   = ctheta_U_cyc_neg  ;
// 2 Pinching
// 2 Offset
    theCopy->cAsymPos    = cAsymPos;
    theCopy->cAsymNeg    = cAsymNeg;
    theCopy->cPos_Offset = cPos_Offset;
    theCopy->cNeg_Offset = cNeg_Offset;
// 3 State
    theCopy->cUi = cUi;
    theCopy->cFi = cFi;
// 2 Stiffness
    theCopy->cKreload = cKreload;
    theCopy->cPosKunload = cPosKunload;
    theCopy->cNegKunload = cNegKunload;
// 2 Energy
    theCopy->cEngAcml = cEngAcml;
    theCopy->cEngDspt = cEngDspt;
// 2 Flag
    theCopy->cFailure_Flag = cFailure_Flag;
    theCopy->cBranch = cBranch;
    theCopy->cPos_ExtreSmallDisp = cPos_ExtreSmallDisp;
    theCopy->cNeg_ExtreSmallDisp = cNeg_ExtreSmallDisp;
    return theCopy;
}

int IMKResiDisp::sendSelf(int cTag, Channel &theChannel)
{
    int res = 0;
    cout << " sendSelf" << endln;

    static Vector data(149);
    data(0) = this->getTag();
    data(1  ) = Ke;
    data(2  ) = posUp_0;
    data(3  ) = posUpc_0;
    data(4  ) = posUu_0;
    data(5  ) = posFy_0;
    data(6  ) = posFcapFy_0;
    data(7  ) = posFresFy_0;
    data(8  ) = negUp_0;
    data(9  ) = negUpc_0;
    data(10 ) = negUu_0;
    data(11 ) = negFy_0;
    data(12 ) = negFcapFy_0;
    data(13 ) = negFresFy_0;
    data(14 ) = LAMBDA_S;
    data(15 ) = LAMBDA_C;
    data(16 ) = LAMBDA_A;
    data(18 ) = c_S;
    data(19 ) = c_C;
    data(20 ) = c_A;
    data(24 ) = alpha;
    data(25 ) = Offset;
    data(26 ) = beta_F;
    data(27 ) = beta_U;
    data(28 ) = theta_F;
    data(29 ) = theta_U;
    data(32 ) = posUy_0;
    data(33 ) = negUy_0;
    data(34 ) = posUcap_0;
    data(35 ) = negUcap_0;
    data(36 ) = posFcap_0;
    data(37 ) = negFcap_0;
    data(38 ) = posKp_0;
    data(39 ) = negKp_0;
    data(40 ) = posKpc_0;
    data(41 ) = negKpc_0;
    data(42 ) = engRefS;
    data(43 ) = engRefC;
    data(44 ) = engRefA;
    data(46 ) = posUy;
    data(47 ) = cPosUy;
    data(48 ) = posFy;
    data(49 ) = cPosFy;
    data(50 ) = posUcap;
    data(51 ) = cPosUcap;
    data(52 ) = posFcap;
    data(53 ) = cPosFcap;
    data(54 ) = posUlocal;
    data(55 ) = cPosUlocal;
    data(56 ) = posFlocal;
    data(57 ) = cPosFlocal;
    data(58 ) = posUglobal;
    data(59 ) = cPosUglobal;
    data(60 ) = posFglobal;
    data(61 ) = cPosFglobal;
    data(62 ) = posUres;
    data(63 ) = cPosUres;
    data(64 ) = posFres;
    data(65 ) = cPosFres;
    data(66 ) = posKp;
    data(67 ) = cPosKp;
    data(68 ) = posKpc;
    data(69 ) = cPosKpc;
    data(70 ) = pos_Uinter_local;
    data(71 ) = cPos_Uinter_local;
    data(72 ) = pos_Uinter_global;
    data(73 ) = cPos_Uinter_global;
    data(74 ) = beta_F_pos;
    data(75 ) = cbeta_F_pos;
    data(76 ) = beta_U_pos;
    data(77 ) = cbeta_U_pos;
    data(78 ) = theta_F_cyc_pos;
    data(79 ) = ctheta_F_cyc_pos;
    data(80 ) = theta_U_cyc_pos;
    data(81 ) = ctheta_U_cyc_pos;
    data(82 ) = negUy;
    data(83 ) = cNegUy;
    data(84 ) = negFy;
    data(85 ) = cNegFy;
    data(86 ) = negUcap;
    data(87 ) = cNegUcap;
    data(88 ) = negFcap;
    data(89 ) = cNegFcap;
    data(90 ) = negUlocal;
    data(91 ) = cNegUlocal;
    data(92 ) = negFlocal;
    data(93 ) = cNegFlocal;
    data(94 ) = negUglobal;
    data(95 ) = cNegUglobal;
    data(96 ) = negFglobal;
    data(97 ) = cNegFglobal;
    data(98 ) = negUres;
    data(99 ) = cNegUres;
    data(100) = negFres;
    data(101) = cNegFres;
    data(102) = negKp;
    data(103) = cNegKp;
    data(104) = negKpc;
    data(105) = cNegKpc;
    data(106) = neg_Uinter_local;
    data(107) = cNeg_Uinter_local;
    data(108) = neg_Uinter_global;
    data(109) = cNeg_Uinter_global;
    data(110) = beta_F_neg;
    data(111) = cbeta_F_neg;
    data(112) = beta_U_neg;
    data(113) = cbeta_U_neg;
    data(114) = theta_F_cyc_neg;
    data(115) = ctheta_F_cyc_neg;
    data(116) = theta_U_cyc_neg;
    data(117) = ctheta_U_cyc_neg;
    data(118) = AsymPos;
    data(119) = cAsymPos;
    data(120) = AsymNeg;
    data(121) = cAsymNeg;
    data(122) = pos_Offset;
    data(123) = cPos_Offset;
    data(124) = neg_Offset;
    data(125) = cNeg_Offset;
    data(126) = Ui;
    data(127) = cUi;
    data(128) = Fi;
    data(129) = cFi;
    data(130) = Kreload;
    data(131) = cKreload;
    data(132) = KgetTangent;
    data(133) = posKunload;
    data(134) = cPosKunload;
    data(135) = negKunload;
    data(136) = cNegKunload;
    data(137) = engAcml;
    data(138) = cEngAcml;
    data(139) = engDspt;
    data(140) = cEngDspt;
    data(141) = Failure_Flag;
    data(142) = cFailure_Flag;
    data(143) = Branch;
    data(144) = cBranch;
    data(145) = pos_ExtreSmallDisp;
    data(146) = cPos_ExtreSmallDisp;
    data(147) = neg_ExtreSmallDisp;
    data(148) = cNeg_ExtreSmallDisp;
    res = theChannel.sendVector(this->getDbTag(), cTag, data);
    if (res < 0)
        opserr << "IMKResiDisp::sendSelf() - failed to send data\n";
    return res;
}

int IMKResiDisp::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    int res = 0;
    static Vector data(149);
    res = theChannel.recvVector(this->getDbTag(), cTag, data);

    if (res < 0) {
        opserr << "IMKResiDisp::recvSelf() - failed to receive data\n";
        this->setTag(0);
    }
    else {
        cout << " recvSelf" << endln;
        this->setTag((int)data(0));
        Ke = data(1  );
        posUp_0 = data(2  );
        posUpc_0 = data(3  );
        posUu_0 = data(4  );
        posFy_0 = data(5  );
        posFcapFy_0 = data(6  );
        posFresFy_0 = data(7  );
        negUp_0 = data(8  );
        negUpc_0 = data(9  );
        negUu_0 = data(10 );
        negFy_0 = data(11 );
        negFcapFy_0 = data(12 );
        negFresFy_0 = data(13 );
        LAMBDA_S = data(14 );
        LAMBDA_C = data(15 );
        LAMBDA_A = data(16 );
        c_S = data(18 );
        c_C = data(19 );
        c_A = data(20 );
        alpha = data(24 );
        Offset = data(25 );
        beta_F = data(26 );
        beta_U = data(27 );
        theta_F = data(28 );
        theta_U = data(29 );
        posUy_0 = data(32 );
        negUy_0 = data(33 );
        posUcap_0 = data(34 );
        negUcap_0 = data(35 );
        posFcap_0 = data(36 );
        negFcap_0 = data(37 );
        posKp_0 = data(38 );
        negKp_0 = data(39 );
        posKpc_0 = data(40 );
        negKpc_0 = data(41 );
        engRefS = data(42 );
        engRefC = data(43 );
        engRefA = data(44 );
        posUy = data(46 );
        cPosUy = data(47 );
        posFy = data(48 );
        cPosFy = data(49 );
        posUcap = data(50 );
        cPosUcap = data(51 );
        posFcap = data(52 );
        cPosFcap = data(53 );
        posUlocal = data(54 );
        cPosUlocal = data(55 );
        posFlocal = data(56 );
        cPosFlocal = data(57 );
        posUglobal = data(58 );
        cPosUglobal = data(59 );
        posFglobal = data(60 );
        cPosFglobal = data(61 );
        posUres = data(62 );
        cPosUres = data(63 );
        posFres = data(64 );
        cPosFres = data(65 );
        posKp = data(66 );
        cPosKp = data(67 );
        posKpc = data(68 );
        cPosKpc = data(69 );
        pos_Uinter_local = data(70 );
        cPos_Uinter_local = data(71 );
        pos_Uinter_global = data(72 );
        cPos_Uinter_global = data(73 );
        beta_F_pos = data(74 );
        cbeta_F_pos = data(75 );
        beta_U_pos = data(76 );
        cbeta_U_pos = data(77 );
        theta_F_cyc_pos = data(78 );
        ctheta_F_cyc_pos = data(79 );
        theta_U_cyc_pos = data(80 );
        ctheta_U_cyc_pos = data(81 );
        negUy = data(82 );
        cNegUy = data(83 );
        negFy = data(84 );
        cNegFy = data(85 );
        negUcap = data(86 );
        cNegUcap = data(87 );
        negFcap = data(88 );
        cNegFcap = data(89 );
        negUlocal = data(90 );
        cNegUlocal = data(91 );
        negFlocal = data(92 );
        cNegFlocal = data(93 );
        negUglobal = data(94 );
        cNegUglobal = data(95 );
        negFglobal = data(96 );
        cNegFglobal = data(97 );
        negUres = data(98 );
        cNegUres = data(99 );
        negFres = data(100);
        cNegFres = data(101);
        negKp = data(102);
        cNegKp = data(103);
        negKpc = data(104);
        cNegKpc = data(105);
        neg_Uinter_local = data(106);
        cNeg_Uinter_local = data(107);
        neg_Uinter_global = data(108);
        cNeg_Uinter_global = data(109);
        beta_F_neg = data(110);
        cbeta_F_neg = data(111);
        beta_U_neg = data(112);
        cbeta_U_neg = data(113);
        theta_F_cyc_neg = data(114);
        ctheta_F_cyc_neg = data(115);
        theta_U_cyc_neg = data(116);
        ctheta_U_cyc_neg = data(117);
        AsymPos = data(118);
        cAsymPos = data(119);
        AsymNeg = data(120);
        cAsymNeg = data(121);
        pos_Offset = data(122);
        cPos_Offset = data(123);
        neg_Offset = data(124);
        cNeg_Offset = data(125);
        Ui = data(126);
        cUi = data(127);
        Fi = data(128);
        cFi = data(129);
        Kreload = data(130);
        cKreload = data(131);
        KgetTangent = data(132);
        posKunload = data(133);
        cPosKunload = data(134);
        negKunload = data(135);
        cNegKunload = data(136);
        engAcml = data(137);
        cEngAcml = data(138);
        engDspt = data(139);
        cEngDspt = data(140);
        Failure_Flag = data(141);
        cFailure_Flag = data(142);
        Branch = data(143);
        cBranch = data(144);
        pos_ExtreSmallDisp = data(145);
        cPos_ExtreSmallDisp = data(146);
        neg_ExtreSmallDisp = data(147);
        cNeg_ExtreSmallDisp = data(148);
    }
    return res;
}

void IMKResiDisp::Print(OPS_Stream &s, int flag)
{
    cout << "IMKResiDisp tag: " << this->getTag() << endln;
}
