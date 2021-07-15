//
//  surface.c
//  softBidispersityOnEllipsoid
//
//  Created by Zhaoyu Xie on 2/24/20.
//  Copyright © 2020 Zhaoyu Xie. All rights reserved.
//

#include "surface.h"
#include "operation.h"
#include <math.h>

//a is the long axis, b is the short axis
void tangent1(double a, double b, double pos[], double t1[]) {
    double x = pos[0];
    double y = pos[1];
    double z = pos[2];
    double rho = sqrt(x*x+y*y);
    
    if(rho < MACHINE_EPSILON) {
        t1[0] = 1;
        t1[1] = 0;
        t1[2] = 0;
    }
    else {
        t1[0] = -y/rho;
        t1[1] = x/rho;
        t1[2] = 0;
    }
}

void tangent2(double a, double b, double pos[], double t2[]) {
    double x = pos[0];
    double y = pos[1];
    double z = pos[2];
    double rho = sqrt(x*x+y*y);
    
    if(a-z < MACHINE_EPSILON) {
        t2[0] = 0;
        t2[1] = 1;
        t2[2] = 0;
    }
    else if(a+z < MACHINE_EPSILON) {
        t2[0] = 0;
        t2[1] = -1;
        t2[2] = 0;
    }
    else {
        t2[0] = -b*x*z/a/rho;
        t2[1] = -b*y*z/a/rho;
        t2[2] = sqrt(a*a-z*z);
        vectorNormalize(t2, t2);
    }
}

void normal(double a, double b, double pos[], double n[]) {
    double x = pos[0];
    double y = pos[1];
    double z = pos[2];
    double rho = sqrt(x*x+y*y);
    
    if(a-z < MACHINE_EPSILON) {
        n[0] = 0;
        n[1] = 0;
        n[2] = 1;
    }
    else if(a+z < MACHINE_EPSILON) {
        n[0] = 0;
        n[1] = 0;
        n[2] = -1;
    }
    else {
        n[0] = a*x/rho;
        n[1] = a*y/rho;
        n[2] = b*z/sqrt(a*a-z*z);
        vectorNormalize(n, n);
    }
}

/* return the gradient of the constraint function at a given position
 * In      - a      - ellipsoid semi major axis length
 *         - b      - ellipsoid semi minor axis length
 *         - x      - position (length 3 array)
 * Out     - grad   - constraint function gradient
 * Returns - none
 */
void constraintGradient(double a, double b, double *x, double *grad) {
    grad[0]=2*x[0]/(b*b);
    grad[1]=2*x[1]/(b*b);
    grad[2]=2*x[2]/(a*a);
}

double constraintFunction(double a, double b, double *x) {
    return (x[0]/b)*(x[0]/b)+(x[1]/b)*(x[1]/b)+(x[2]/a)*(x[2]/a)-1;
}

/* Give the ellipsoid axis parameters for a given shape evolution parameter
 * Input:	s		- evolution parameter, runs from 0 to 1
 *          a0      - inital semi-major axis length
 *          b0      - ''
 * Output:  a		- semi-major axis length at evolution parameter s
 *          b       - ''
 * Returns:	Nothing
 */
void ellipsoidUpdate(double s, double a0, double b0, double *a, double *b) {
    *b = b0+(pow(b0*b0*a0,1.0/3.0)-b0)*s;
    *a = a0*(b0*b0)/((*b)*(*b));
}

/* Calculates the surface area of an ellipsoid given the semi-major/minor axis lengths
 * In      - a      - semi-major axis of the ellipsoid
 *           b      - semi-minor
 * Out     - none
 * Returns - surface area
 */
double ellipsoidSurfaceArea(double a, double b) {
    double e;
    if (fabs(a-b) < MACHINE_EPSILON) {
        // sphere
        return 4*PI*a*a;
    } else if (a>b) {
        // prolate
        e = sqrt(a*a-b*b)/a;
        return 2*PI*b*(b+a/e*asin(e));
    } else {
        // oblate
        e = sqrt(b*b-a*a)/b;
        return PI*(2*b*b+a*a/e*log((1+e)/(1-e)));
    }
}

/* Give the surface area for a shape evolution parameter s given an initial shape
 * Input:	a0      - inital semi-major axis length
 *          b0      - ''
 *          s       - shape evolution parameter, between 0 and 1
 * Output:  Nothing
 * Returns:	surface area at the given evolution parameter
 */
double surfaceAreaAtEvoParameter(double a0, double b0, double s) {
    double a;
    double b;
    ellipsoidUpdate(s, a0, b0, &a, &b);
    double area = ellipsoidSurfaceArea(a,b);
    return area;
}

/* Given a desired area and an initial shape, return the shape evolution parameter
 * at which that area occurs. Uses bisection method to invert surfaceAreaAtEvoParameter.
 * Input:	a0      - inital semi-major axis length
 *          b0      - ''
 *          area    - the surface area corresponding to the desired evolution parameter
 * Output:  Nothing
 * Returns:	evolution parameter corresponding to the given area
 */
double evoParameterAtSurfaceArea(double a0, double b0, double area) {
    double eps = MACHINE_EPSILON; // precision goal
    double x1 = 0.0; //brackets --- we know the function is monotonically decreasing between these values
    double x2 = 1.0;
    int maxIter = 50; // maximum iterations
    double dx;
    double xMid;
    double root;
    
    double f = surfaceAreaAtEvoParameter(a0,b0,x1)-area;
    double fMid = surfaceAreaAtEvoParameter(a0,b0,x2)-area;
    
    if (f*fMid>0.0) {
        printf("not bracketed correctly: %f %f\n",f,fMid);
    }
    
    //orient search from lower f value
    if (f<=0.0) {
        root = x1;
        dx = x2-x1;
    } else {
        root = x2;
        dx = x1-x2;
    }
    
    for (int i=0; i<maxIter; i++) {
        dx*=0.5;
        xMid = root+dx;
        fMid = surfaceAreaAtEvoParameter(a0,b0,xMid)-area;
        if (fMid<=0.0) {
            root = xMid;
        }
        if (fabs(dx)<eps || fMid==0.0) {
            break;
        }
        if (i==maxIter-1) {
            printf("Too many bisections. dx: %f\n",dx);
        }
    }
    return root;
}

/* Give the surface area at a given time such that the surface area
 * decays exponentially. At t=1, the area is the surface area of the final sphere, which will have decreased by 
 * a fraction f of the total ammount it would have decreased by at t=infinity
 * Input:	a0      - inital semi-major axis length
 *          b0      - ''
 *          f       - fraction of the t=inf decay that will occur by t=1
 *          t       - time, between 0 and 1
 * Output:  Nothing
 * Returns:	surface area at the given time
 */
double areaAtTime(double a0, double b0, double f, double t) {
    double area0 = ellipsoidSurfaceArea(a0,b0);
    double area1 = surfaceAreaAtEvoParameter(a0,b0,1.0);
    double areaInf = area0-(area0-area1)/f;
    double tau = 1/log(1/(1-f));
    return areaInf+(area0-areaInf)*exp(-t/tau);    
}

/* Calculate the conformal factor at a point on the surface
 * In      - a      - semi-major axis of the ellipsoid
 *           b      -    "                            "
 *           u      - scaled z-axis coordinate
 * Out     - none
 * Returns - conformal factor
 */
double conformalFactor(double a, double b, double u) {
    return sqrt(b*b*(a*a+(b*b-a*a)*u*u));
}

/* Given a conformal coordinate, give the corresponding scaled z-axis coordinate.
 * Used for generating an even surface distribution.
 * This is done by inverting an integral equation numerically - this is good enough
 * to generate an evenly distributed starting configuration, which will then be allowed to diffuse,
 * but does not have enough precision for general usage.
 * Includes exact solution for the special case of a sphere.
 * In      - a      - semi-major/minor axis of the ellipsoid
 *           b      -    "                            "
 *           u2     - conformal coordinate - must be in range [0,1]
 * Out     - none
 * Returns - scaled z-axis coordinate
 */
double conformalInvert(double a, double b, double u2) {
    double du=0.00001;
    double sum=0;
    double intMax=0;
    double result=1.0; // will return 1 in case the integral never reaches its max for numerical reasons
    
    if (u2<0-MACHINE_EPSILON || u2>1+MACHINE_EPSILON) {
        printf("conformal coordinate out of range - must be between 0 and 1\n");
        return(0);
    }
    
    if (fabs(a-b) < MACHINE_EPSILON) { //special case: sphere
        result = 2*u2-1;
    } else { // general ellipsoid case
        for (double i=-1.0+du/2; i<1.0; i+=du) {
            intMax += conformalFactor(a,b,i)*du;
        }
        for (double i=-1.0+du/2; i<1.0; i+=du) {
            sum += conformalFactor(a,b,i)*du;
            if (sum > u2*intMax) {
                result = i;
                break;
            }
        }
    }
    
    return(result);
    
}


double packingFraction(particle p[], int np, double a, double b) {
    double sum = 0;
    for(int i=0; i<np; i++)
        sum += PI*p[i].rad*p[i].rad;
    return sum/ellipsoidSurfaceArea(a, b);
}













