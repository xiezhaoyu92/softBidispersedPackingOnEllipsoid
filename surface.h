//
//  surface.h
//  softBidispersityOnEllipsoid
//
//  Created by Zhaoyu Xie on 2/24/20.
//  Copyright © 2020 Zhaoyu Xie. All rights reserved.
//

#ifndef surface_h
#define surface_h

#ifndef MACHINE_EPSILON
#define MACHINE_EPSILON 1e-15
#endif

#ifndef PI
#define PI 3.141592653589793
#endif

#include "particle.h"
#include <stdio.h>

void tangent1(double a, double b, double pos[], double t1[]);
void tangent2(double a, double b, double pos[], double t2[]);
void normal(double a, double b, double pos[], double n[]);

void constraintGradient(double a, double b, double *x, double *grad);
double constraintFunction(double a, double b, double *x);

double areaAtTime(double a0, double b0, double f, double t);
double evoParameterAtSurfaceArea(double a0, double b0, double area);
void ellipsoidUpdate(double s, double a0, double b0, double *a, double *b);

double conformalInvert(double a, double b, double u2);

double packingFraction(particle p[], int np, double a, double b);
#endif /* surface_h */
