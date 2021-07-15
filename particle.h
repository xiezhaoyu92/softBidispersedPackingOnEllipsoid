//
//  particle.h
//  softBidispersityOnEllipsoid
//
//  Created by Zhaoyu Xie on 3/3/20.
//  Copyright © 2020 Zhaoyu Xie. All rights reserved.
//

#ifndef particle_h
#define particle_h

#ifndef MACHINE_EPSILON
#define MACHINE_EPSILON 1e-15
#endif

#ifndef TRUE
#define TRUE -1
#endif

#ifndef FALSE
#define FALSE 0
#endif

#ifndef PI
#define PI 3.141592653589793
#endif

#include <stdio.h>

typedef struct {
    double rad;
    double rr;
    
    double position[3];
    double oldPosition[3];
    double postPreRelaxPosition[3];
    double positionTemp[3];
    
    double force[3];
    double forceStoc[3];
    double oldForce[3];
    double conjGrad[3];
    double conjGradProj[3];
    
    int oldCoord[3];
    int coord[3];
    int postPreRelaxCoord[3];
    int coordTemp[3];
} particle;

int overlapQ(particle *p1, particle *p2);
void addRepulsiveSpringForce(particle p[], int np);
double totalRepulsiveSpringEnergy(particle p[], int np);
void partitionOneParticle(particle *p, int nPart[], double partBoundX[], double partBoundY[], double partBoundZ[]);
void conjGradDescentConstrainedIntegrate(double a, double b, particle *p, double dt);
void scaleOntoNewSurface(particle *p, double aOld, double bOld, double a, double b);
void projectOntoNewSurface(particle *p, double aOld, double bOld, double a, double b);
void conjGradDescentConstrainedIntegrateMinimize(particle p[], int np, double rad, double a, double b, int nPart[], double partBoundX[], double partBoundY[], double partBoundZ[]);
void conjGradDescentConstrainedIntegrateMinimizeSteps(particle p[], int np, double rad, double a, double b, int nPart[], double partBoundX[], double partBoundY[], double partBoundZ[], int steps);
void undoOverlaps(particle p[],  int np, double a, double b, int nPart[], double partBoundX[], double partBoundY[], double partBoundZ[], int *nearJammed);
#endif /* particle_h */
