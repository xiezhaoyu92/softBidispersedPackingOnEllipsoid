//
//  particle.c
//  softBidispersityOnEllipsoid
//
//  Created by Zhaoyu Xie on 3/3/20.
//  Copyright © 2020 Zhaoyu Xie. All rights reserved.
//

#include "particle.h"
#include "operation.h"
#include "surface.h"
#include "hessian.h"
#include <math.h>
#include <stdlib.h>

int overlapQ(particle *p1, particle *p2) {
    if(abs(p1->coord[0]-p2->coord[0])>1||abs(p1->coord[1]-p2->coord[1])>1||abs(p1->coord[2]-p2->coord[2])>1)
        return FALSE;
    
    double r[3];
    vectorSubtract(p1->position, p2->position, r);
    double R = vectorNorm(r);
    if(R<p1->rad+p2->rad)
        return TRUE;
    else
        return FALSE;
}

double repulsiveSpringPotential(particle *p1, particle *p2) {
    if(abs(p1->coord[0]-p2->coord[0])>1||abs(p1->coord[1]-p2->coord[1])>1||abs(p1->coord[2]-p2->coord[2])>1)
        return 0;
    
    double epsilon = 1;
    double sigma = p1->rad+p2->rad;
    double r[3];
    vectorSubtract(p1->position, p2->position, r);
    double R = vectorNorm(r);
    if(R<sigma)
        return epsilon/2*pow(1-R/sigma, 2);
    else
        return 0;
}

void repulsiveSpringForce(particle *p1, particle *p2) {
    if(abs(p1->coord[0]-p2->coord[0])>1||abs(p1->coord[1]-p2->coord[1])>1||abs(p1->coord[2]-p2->coord[2])>1)
        return;
    
    double epsilon = 1;
    double sigma = p1->rad+p2->rad;
    double r[3];
    vectorSubtract(p1->position, p2->position, r);
    double R = vectorNorm(r);
    if(R<sigma) {
        double magnitude = epsilon/sigma*(1-R/sigma);
        for(int i=0; i<3; i++) {
            p1->force[i] += magnitude/R*(p1->position[i]-p2->position[i]);
            p2->force[i] += magnitude/R*(p2->position[i]-p1->position[i]);
        }
    }
    else
        return;
}

void addRepulsiveSpringForce(particle p[], int np) {
    for(int i=0; i<np; i++)
        for(int j=i+1; j<np; j++)
            repulsiveSpringForce(&(p[i]), &(p[j]));
}

double totalRepulsiveSpringEnergy(particle p[], int np) {
    double energy = 0;
    for(int i=0; i<np; i++)
        for(int j=i+1; j<np; j++)
            energy += repulsiveSpringPotential(&(p[i]), &(p[j]));
    return energy;
}

//calculate conjugate gradient force
void calculateConjugateForce(particle p[], int np, double a, double b) {
    addRepulsiveSpringForce(p, np);
    double gammaNum = 0;
    double gammaDenom = 0;
    for(int i=0; i<np; i++)
        for(int j=0; j<3; j++) {
            gammaNum += p[i].force[j]*p[i].force[j];
            gammaDenom += p[i].oldForce[j]*p[i].oldForce[j];
        }
    for(int i=0; i<np; i++)
        for(int j=0; j<3; j++) {
            p[i].conjGrad[j] = p[i].force[j]+gammaNum/gammaDenom*p[i].conjGrad[j];
            p[i].oldForce[j] = p[i].force[j];
        }
}
//find the component of conjugate gradient force on the tangent plane
void projectConjugateForce(particle p[], int np, double a, double b) {
    for(int i=0; i<np; i++) {
        //double n[3];
        //normal(a, b, p[i].position, n);
        //double projection = vectorDotProduct(p[i].conjGrad, n);
        for(int j=0; j<3; j++)
            p[i].conjGradProj[j] = p[i].conjGrad[j];
            //p[i].conjGradProj[j] = p[i].conjGrad[j]-projection*n[j];
    }
}

void partitionOneParticle(particle *p, int nPart[], double partBoundX[], double partBoundY[], double partBoundZ[]) {
    int i=0, j=0, k=0;
    while(i<nPart[0]&&p->position[0]>=partBoundX[i])
        i++;
    p->coord[0] = i-1;
    while(j<nPart[1]&&p->position[1]>=partBoundY[j])
        j++;
    p->coord[1] = j-1;
    while(k<nPart[2]&&p->position[2]>=partBoundZ[k])
        k++;
    p->coord[2] = k-1;
}

/* Use the force on a particle to update its position, subject to the surface constraint
 * (enforced by Lagrange multiplier)
 * Input:  - a              - ellipsoid semi major axis length
 *         - b              - ellipsoid semi minor axis length
 *         - p	            - particle to be integrated
 * Output:  p.position	    - particle's new position
 * Returns:	Nothing
 */
void conjGradDescentConstrainedIntegrate(double a, double b, particle *p, double dt) {
    double newPosition[3];
    double sigmaNew; // constraint function value
    double lambda; // Lagrange muliplier
    double gradOld[3]; // gradient of constraint function at old position
    double gradNew[3]; // gradient of constrain function at new position
    double updateForce[3];
    
    for (int i=0; i<3; i++) {
        newPosition[i] = p->position[i] + (p->conjGradProj[i])*dt;
    }
    
    constraintGradient(a,b,p->position,gradOld);
    sigmaNew = constraintFunction(a, b, newPosition);
    do {
        constraintGradient(a,b,newPosition,gradNew);
        lambda = sigmaNew/vectorDotProduct(gradNew, gradNew);
        vectorScale(gradNew, -lambda, updateForce);
        vectorAdd(newPosition, updateForce, newPosition);
        sigmaNew = constraintFunction(a, b, newPosition);
    } while (fabs(sigmaNew) > MACHINE_EPSILON);
    vectorCopy(newPosition,p->position);
}

//Two options can be used to move particles onto the surface: scale when the shape shrink size, project along normal directions when shape deform.
void scaleOntoNewSurface(particle *p, double aOld, double bOld, double a, double b) {
    for(int j=0; j<3; j++)
        p->position[j] = p->position[j]/aOld*a;
}
void projectOntoNewSurface(particle *p, double aOld, double bOld, double a, double b) {
    double newPosition[3];
    double sigmaNew;
    double lambda;
    double gradOld[3];
    double gradNew[3];
    double updateForce[3];
    
    for (int i=0; i<3; i++) {
        newPosition[i] = p->position[i];
    }
    
    normal(aOld, bOld, p->position, gradOld);
    //constraintGradient(aOld, bOld, p->position, gradOld);
    sigmaNew = constraintFunction(a, b, p->position);
    do {
        constraintGradient(a, b, newPosition, gradNew);
        lambda = sigmaNew/vectorDotProduct(gradNew, gradNew);
        vectorScale(gradNew, -lambda, updateForce);
        vectorAdd(newPosition, updateForce, newPosition);
        sigmaNew = constraintFunction(a, b, newPosition);
    } while (fabs(sigmaNew) > MACHINE_EPSILON);
    vectorCopy(newPosition, p->position);
}

double energyForTimestep(particle p[], int np, double dt, double a, double b, int nPart[], double partBoundX[], double partBoundY[], double partBoundZ[]) {
    for (int i=0; i<np; i++)
        for(int j=0; j<3; j++) {
            p[i].oldPosition[j] = p[i].position[j];
            p[i].oldCoord[j] = p[i].coord[j];
        }
    for(int i=0; i<np; i++) {
        conjGradDescentConstrainedIntegrate(a, b, &(p[i]), dt);
        partitionOneParticle(&(p[i]), nPart, partBoundX, partBoundY, partBoundZ);
    }
    double energy = totalRepulsiveSpringEnergy(p, np);
    for (int i=0; i<np; i++) {
        for(int j=0; j<3; j++)
            p[i].position[j] = p[i].oldPosition[j];
        for(int j=0; j<3; j++) {
            p[i].coord[j] = p[i].oldCoord[j];
        }
    }
    return energy;
}

double goldenSearch(particle p[], int np, double rad, double a, double b, int nPart[], double partBoundX[], double partBoundY[], double partBoundZ[]) {
    double phi = 1.618033988749895;//goldenRatio+1
    double w = 0.3819660112501051;//1-goldenRatio
    double tol = MACHINE_EPSILON;
    double currentEnergy = totalRepulsiveSpringEnergy(p, np);
    
    double maxForce = 0;
    for(int i=0; i<np; i++)
        if(vectorNorm(p[i].conjGradProj)>maxForce)
            maxForce = vectorNorm(p[i].conjGradProj);
    double bracketLimit;
    double initialBracketGuess;
    if(maxForce>0) {
        bracketLimit = rad/maxForce;
        initialBracketGuess = 0.1*rad/maxForce;
    }
    else
        return 1;// return any finite real number
    
    double dtA, dtB, enA, enB;
    dtA = initialBracketGuess;
    enA = energyForTimestep(p, np, dtA, a, b, nPart, partBoundX, partBoundY, partBoundZ);
    dtB = dtA*phi;
    enB = energyForTimestep(p, np, dtB, a, b, nPart, partBoundX, partBoundY, partBoundZ);
    if(enA == 0)
        return dtA;
    do {
        if(enA >= currentEnergy) {
            if(dtA >= MACHINE_EPSILON){
                dtB = dtA;
                enB = enA;
                dtA = dtA/phi;
                enA = energyForTimestep(p, np, dtA, a, b, nPart, partBoundX, partBoundY, partBoundZ);
                if(enA == 0)
                    return dtA;
            } else
                return 0;
        } else if(enB <= enA) {
            dtA = dtB;
            enA = enB;
            dtB = dtB*phi;
            enB = energyForTimestep(p, np, dtB, a, b, nPart, partBoundX, partBoundY, partBoundZ);
            if(enA == 0)
                return dtA;
            if(dtA > bracketLimit)
                return bracketLimit;
        }
    } while(enA>=currentEnergy || enB<=enA);
    //Now the initial bracket is (0, dtA, dtB). Find the minimum next.
    double dtx, dty;
    double enx, eny;
    dtx = dtA;
    enx = enA;
    dtA = 0;
    enA = currentEnergy;
    //(dtA, dtx, dty, dtB)
    do {
        if(dtB-dtx > dtx-dtA) {
            dty = dtx + w*(dtB-dtx);
            eny = energyForTimestep(p, np, dty, a, b, nPart, partBoundX, partBoundY, partBoundZ);
        }
        else {
            dty = dtx;
            eny = enx;
            dtx = dtx - w*(dtx-dtA);
            enx = energyForTimestep(p, np, dtx, a, b, nPart, partBoundX, partBoundY, partBoundZ);
        }
        if(enx<eny) {
            dtB = dty;
            enB = eny;
        } else {
            dtA = dtx;
            enA = enx;
            dtx = dty;
            enx = eny;
        }
    } while(dtB-dtA > tol*(dtx+dty));
    //printf("dt: %lf\n", dtx);
    return dtx;
}

void oneConjGradDescentConstrainedIntegrateMove(particle p[], int np, double dt, double a, double b, int nPart[], double partBoundX[], double partBoundY[], double partBoundZ[]) {
    for(int i=0; i<np; i++) {
        conjGradDescentConstrainedIntegrate(a, b, &(p[i]), dt);
        partitionOneParticle(&(p[i]), nPart, partBoundX, partBoundY, partBoundZ);
    }
}

void addStochasticForce(particle *p, double D, double a, double b, double rad) {
    double t1[3];
    double t2[3];
    tangent1(a, b, p->position, t1);
    tangent2(a, b, p->position, t2);
    double ETA1, ETA2;
    ETA1 = randNormal();
    ETA2 = randNormal();
    
    for (int i=0; i<3; i++){
        p->conjGradProj[i] = sqrt(2*D)/rad*(p->rad)*randNormal();
    }
}
double oneParticleEnergy(int i, particle p[], int np) {
    double energy = 0;
    for(int j=0; j<np; j++) {
        if(j==i)
            continue;
        energy += repulsiveSpringPotential(&(p[i]), &(p[j]));
    }
    return energy;
}
double MonteCarloMoves(particle p[], int np, double dt, double a, double b, int nPart[], double partBoundX[], double partBoundY[], double partBoundZ[]) {
    // going to loop through particles in a random order - set up shuffled array of indices
    int randIndices[np];
    for (int k=0; k<np; k++) {
        randIndices[k]=k;
    }
    for (int k=np-1; k>0; k--) {
        int swapWith = randInteger(k-1);
        int hold = randIndices[k];
        randIndices[k] = randIndices[swapWith];
        randIndices[swapWith] = hold;
    }
    
    int acceptance = 0;
    for (int k0 = 0; k0 < np; k0++) {
        int k = randIndices[k0];
        for(int j=0; j<3; j++) {
            p[k].oldPosition[j] = p[k].position[j];
            p[k].oldCoord[j] = p[k].coord[j];
        }
        double energyOneOld = oneParticleEnergy(k, p, np);
        
        conjGradDescentConstrainedIntegrate(a, b, &(p[k]), sqrt(dt));
        partitionOneParticle(&(p[k]), nPart, partBoundX, partBoundY, partBoundZ);
        double energyOne = oneParticleEnergy(k, p, np);
        if(energyOne>energyOneOld) {
            for(int j=0; j<3; j++) {
                p[k].position[j] = p[k].oldPosition[j];
                p[k].coord[j] = p[k].oldCoord[j];
            }
        }
        else
            acceptance += 1;
    }
    double acceptanceRatio = acceptance/np;
    return acceptanceRatio;
}

void conjGradDescentConstrainedIntegrateMinimize(particle p[], int np, double rad, double a, double b, int nPart[], double partBoundX[], double partBoundY[], double partBoundZ[]) {
    double energy = totalRepulsiveSpringEnergy(p, np);
    double energyOld = energy;
    double dt;
    double diffusionCoeff = pow(2*rad/1000,2)/2;
    
    for(int i=0; i<np; i++)
        for(int j=0; j<3; j++) {
            p[i].force[j] = 0;
            p[i].conjGrad[j] = 0;
            p[i].conjGradProj[j] = 0;
        }
    addRepulsiveSpringForce(p, np);
/*
    for(int i=0; i<np; i++) {
        printf("%d: %.20lf, %.20lf, %.20lf\n", i, p[i].force[0], p[i].force[1], p[i].force[2]);
    }
*/
    for(int i=0; i<np; i++)
        for(int j=0; j<3; j++) {
            p[i].conjGrad[j] = p[i].force[j];
            p[i].oldForce[j] = p[i].force[j];
        }
    projectConjugateForce(p, np, a, b);
/*
    for(int i=0; i<np; i++) {
        printf("%d: %.20lf, %.20lf, %.20lf\n", i, p[i].conjGradProj[0], p[i].conjGradProj[1], p[i].conjGradProj[2]);
    }
*/    
    int conjGradCount = 0;
    
    while(1) {
        dt = goldenSearch(p, np, rad, a, b, nPart, partBoundX, partBoundY, partBoundZ);
        if(dt==0) {
            if(conjGradCount==0)
                break;
            else {
                conjGradCount = 0;
                for(int i=0; i<np; i++)
                    for(int j=0; j<3; j++) {
                        p[i].conjGrad[j] = p[i].force[j];
                        p[i].oldForce[j] = p[i].force[j];
                    }
                projectConjugateForce(p, np, a, b);
            }
        }
        else {
            oneConjGradDescentConstrainedIntegrateMove(p, np, dt, a, b, nPart, partBoundX, partBoundY, partBoundZ);
            conjGradCount++;
            if(conjGradCount>=1000){
/*
                do{
                int info;
                solveHessian(p, np, nPart, partBoundX, partBoundY, partBoundZ, &info);
                if(info==0) {
                    printf("Do one newton move!\n");
                    for(int j=0; j<np; j++)
                        for(int k=0; k<3; k++)
                            p[j].conjGrad[k] = p[j].force[k];
                    for(int i=0; i<np; i++) {
                        double n[3];
                        normal(a, b, p[i].position, n);
                        double projection = vectorDotProduct(p[i].conjGrad, n);
                        for(int j=0; j<3; j++)
                            //p[i].conjGradProj[j] = p[i].conjGrad[j];
                            p[i].conjGradProj[j] = p[i].conjGrad[j]-projection*n[j];
                    }                    
                    dt = goldenSearch(p, np, rad, a, b, nPart, partBoundX, partBoundY, partBoundZ);
                    oneConjGradDescentConstrainedIntegrateMove(p, np, dt, a, b, nPart, partBoundX, partBoundY, partBoundZ);
                }
                else {
                    printf("newton method failed!\n");
                    break;
                }
                energy = totalRepulsiveSpringEnergy(p, np);
                printf("%.25lf\n", energy/np);
                } while(dt!=0);
*/
                
/*
                for(int i=0; i<np; i++)
                    for(int j=0; j<3; j++) {
                        p[i].force[j] = 0;
                    }
                calculateConjugateForce(p, np, a, b);
                projectConjugateForce(p, np, a, b);
                FILE *configRecord=fopen("configRecord.dat","a");
                for(int i=0; i<np; i++)
                    fprintf(configRecord, "%.15lf %.15lf %.15lf %.15lf %.15lf %.15lf\n", p[i].position[0], p[i].position[1], p[i].position[2], p[i].conjGrad[0], p[i].conjGrad[1], p[i].conjGrad[2]);
                fclose(configRecord);
 */
                conjGradCount = 0;
                energyOld = totalRepulsiveSpringEnergy(p, np);
                for(int i=0; i<np; i++)
                    for(int j=0; j<3; j++) {
                        p[i].force[j] = 0;
                        p[i].conjGrad[j] = 0;
                        p[i].conjGradProj[j] = 0;
                    }
                addRepulsiveSpringForce(p, np);
                for(int i=0; i<np; i++)
                    for(int j=0; j<3; j++) {
                        p[i].conjGrad[j] = p[i].force[j];
                        p[i].oldForce[j] = p[i].force[j];
                    }
                projectConjugateForce(p, np, a, b);
                continue;
            }
            energy = totalRepulsiveSpringEnergy(p, np);
            //printf("%d: %.25lf\n", conjGradCount, energy/np);
            if(energy/np<1e-16 || fabs(energy-energyOld)/energy<1e-16)
                break;
            energyOld = energy;
            for(int i=0; i<np; i++)
                for(int j=0; j<3; j++) {
                    p[i].force[j] = 0;
                }
            calculateConjugateForce(p, np, a, b);
            projectConjugateForce(p, np, a, b);
/*
            FILE *configurationTest=NULL;
            configurationTest = fopen("configurationTest.asc","w");
            for (int j = 0; j < np; j++)
                fprintf(configurationTest, "%.15lf %.15lf %.15lf\n", p[j].position[0], p[j].position[1], p[j].position[2]);
            fclose(configurationTest);
            FILE *conjForceTest=NULL;
            conjForceTest = fopen("conjForceTest.dat","w");
            for(int j=0; j<np; j++)
                fprintf(conjForceTest,"%.15lf %.15lf %.15lf\n", p[j].conjGradProj[0], p[j].conjGradProj[1], p[j].conjGradProj[2]);
            fclose(conjForceTest);
*/            
        }
    }
}
void conjGradDescentConstrainedIntegrateMinimizeSteps(particle p[], int np, double rad, double a, double b, int nPart[], double partBoundX[], double partBoundY[], double partBoundZ[], int steps) {
    double energy = totalRepulsiveSpringEnergy(p, np);
    double energyOld = energy;
    double dt;
    
    for(int i=0; i<np; i++)
        for(int j=0; j<3; j++) {
            p[i].force[j] = 0;
            p[i].conjGrad[j] = 0;
            p[i].conjGradProj[j] = 0;
        }
    addRepulsiveSpringForce(p, np);
    for(int i=0; i<np; i++)
        for(int j=0; j<3; j++) {
            p[i].conjGrad[j] = p[i].force[j];
            p[i].oldForce[j] = p[i].force[j];
        }
    projectConjugateForce(p, np, a, b);

    int conjGradCount = 0;
    int totalCount = 0;
    
    while(1) {
        dt = goldenSearch(p, np, rad, a, b, nPart, partBoundX, partBoundY, partBoundZ);
        if(dt==0) {
            if(conjGradCount==0)
                break;
            else {
                conjGradCount = 0;
                for(int i=0; i<np; i++)
                    for(int j=0; j<3; j++) {
                        p[i].conjGrad[j] = p[i].force[j];
                        p[i].oldForce[j] = p[i].force[j];
                    }
                projectConjugateForce(p, np, a, b);
            }
        }
        else {
            oneConjGradDescentConstrainedIntegrateMove(p, np, dt, a, b, nPart, partBoundX, partBoundY, partBoundZ);
            conjGradCount++;
            totalCount++;
            if(conjGradCount>=1000){
                conjGradCount = 0;
                energyOld = totalRepulsiveSpringEnergy(p, np);
                for(int i=0; i<np; i++)
                    for(int j=0; j<3; j++) {
                        p[i].force[j] = 0;
                        p[i].conjGrad[j] = 0;
                        p[i].conjGradProj[j] = 0;
                    }
                addRepulsiveSpringForce(p, np);
                for(int i=0; i<np; i++)
                    for(int j=0; j<3; j++) {
                        p[i].conjGrad[j] = p[i].force[j];
                        p[i].oldForce[j] = p[i].force[j];
                    }
                projectConjugateForce(p, np, a, b);
                continue;
            }
            energy = totalRepulsiveSpringEnergy(p, np);
            //printf("%d: %.25lf\n", conjGradCount, energy/np);
            if(energy/np<1e-16 || fabs(energy-energyOld)/energyOld<1e-8 || totalCount>steps)
                break;
            energyOld = energy;
            for(int i=0; i<np; i++)
                for(int j=0; j<3; j++) {
                    p[i].force[j] = 0;
                }
            calculateConjugateForce(p, np, a, b);
            projectConjugateForce(p, np, a, b);
        }
    }
}

//functions to remove overlaps in the initial state
void overlapForce(particle *p1, particle *p2) {
    double r[3];
    double R;
    vectorSubtract(p1->position, p2->position, r);
    R = vectorNorm(r);
    for(int i=0; i<3; i++){
        p1->conjGradProj[i] += (p1->rad+p2->rad)*r[i]/R;
        p2->conjGradProj[i] += -(p1->rad+p2->rad)*r[i]/R;
    }
}

void undoOverlaps(particle p[],  int np, double a, double b, int nPart[], double partBoundX[], double partBoundY[], double partBoundZ[], int *nearJammed) {
    int maxUndoSteps = 10000;
    double dtOverlap = 0.001;
    
    int undoSteps = 0;
    *nearJammed = FALSE;
    int totalOverlapQ;
    do {
        totalOverlapQ = FALSE;
        undoSteps++;
        //printf("%d\n", undoSteps);
        if(undoSteps>maxUndoSteps) {
            *nearJammed = TRUE;
            break;
        }
        for(int i=0; i<np; i++)
            for(int j=0; j<3; j++)
                p[i].conjGradProj[j] = 0;
        for(int i=0; i<np; i++)
            for(int j=i+1; j<np; j++)
                if (overlapQ(&(p[i]), &(p[j]))) {
                    overlapForce(&(p[i]), &(p[j]));
                    totalOverlapQ = TRUE;
                }
        if(!totalOverlapQ)
            break;
        
        oneConjGradDescentConstrainedIntegrateMove(p, np, dtOverlap, a, b, nPart, partBoundX, partBoundY, partBoundZ);        
    } while(totalOverlapQ);
    //printf("undoOverlapSteps are %d\n", undoSteps);
}
