//
//  hessian.c
//  softBidispersityOnEllipsoid
//
//  Created by Zhaoyu Xie on 6/6/20.
//  Copyright © 2020 Zhaoyu Xie. All rights reserved.
//

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "hessian.h"
#include "solveHessian.h"
#include "particle.h"
#include "operation.h"


double Vii(int i, particle p1, particle p2) {
    double epsilon = 1;
    double r[3];
    vectorSubtract(p1.position, p2.position, r);
    double R = vectorNorm(r);
    double sigma = p1.rad+p2.rad;
    return epsilon/(sigma*sigma)*pow(p1.position[i]-p2.position[i],2)/(R*R)-epsilon/sigma*(1-R/sigma)*(1/R-pow(p1.position[i]-p2.position[i],2)/(R*R*R));
}
double Vij(int i, int j, particle p1, particle p2) {
    double epsilon = 1;
    double r[3];
    vectorSubtract(p1.position, p2.position, r);
    double R = vectorNorm(r);
    double sigma = p1.rad+p2.rad;
    return epsilon/(sigma*sigma)*(p1.position[i]-p2.position[i])*(p1.position[j]-p2.position[j])/(R*R)+epsilon/sigma*(1-R/sigma)*(p1.position[i]-p2.position[i])*(p1.position[j]-p2.position[j])/(R*R*R);
}
double ViVi(int i, particle p1, particle p2) {
    double epsilon = 1;
    double r[3];
    vectorSubtract(p1.position, p2.position, r);
    double R = vectorNorm(r);
    double sigma = p1.rad+p2.rad;
    return -epsilon/(sigma*sigma)*pow(p1.position[i]-p2.position[i],2)/(R*R)-epsilon/sigma*(1-R/sigma)*(-1/R+pow(p1.position[i]-p2.position[i],2)/(R*R*R));
}
double ViVj(int i, int j, particle p1, particle p2) {
    double epsilon = 1;
    double r[3];
    vectorSubtract(p1.position, p2.position, r);
    double R = vectorNorm(r);
    double sigma = p1.rad+p2.rad;
    return -epsilon/(sigma*sigma)*(p1.position[i]-p2.position[i])*(p1.position[j]-p2.position[j])/(R*R)-epsilon/sigma*(1-R/sigma)*(p1.position[i]-p2.position[i])*(p1.position[j]-p2.position[j])/(R*R*R);
}
double secondTotalD(int m, int i, particle p[], int np, int nPart[], double partBoundX[], double partBoundY[], double partBoundZ[]) {
    double result = 0;
    for(int n=0; n<np; n++) {
        if(n==m||!overlapQ(&(p[m]), &(p[n])))
            continue;
        result += Vii(i, p[m], p[n]);
    }
    return result;
}
double secondPartialDSingle(int m, int i, int j, particle p[], int np, int nPart[], double partBoundX[], double partBoundY[], double partBoundZ[]){
    double result = 0;
    for(int n=0; n<np; n++) {
        if(n==m||!overlapQ(&(p[m]), &(p[n])))
            continue;
        result += Vij(i, j, p[m], p[n]);
    }
    return result;
}
double secondPartialDPair(int m, int n, int i, int j, particle p[], int np, int nPart[], double partBoundX[], double partBoundY[], double partBoundZ[]) {
    if(i==j)
        return ViVi(i, p[m], p[n]);
    else
        return ViVj(i, j, p[m], p[n]);
}

int indexColumn(int np, int row, int column) {
    return column*3*np+row;
}
void hessianMatrix(particle p[], int np, int nPart[], double partBoundX[], double partBoundY[], double partBoundZ[], double *hessian) {
    for(int m=0; m<np; m++) {
        for(int i=0; i<3; i++)
            hessian[indexColumn(np,3*m+i,3*m+i)] = secondTotalD(m, i, p, np, nPart, partBoundX, partBoundY, partBoundZ);
        for(int i=0; i<3; i++)
            for(int j=i+1; j<3; j++) {
                hessian[indexColumn(np,3*m+i,3*m+j)] = secondPartialDSingle(m, i, j, p, np, nPart, partBoundX, partBoundY, partBoundZ);
                hessian[indexColumn(np,3*m+j,3*m+i)] = hessian[indexColumn(np,3*m+i,3*m+j)];
            }
       
        for(int n=m+1; n<np; n++) {
            if(!overlapQ(&(p[m]), &(p[n])))
                continue;
            for(int i=0; i<3; i++)
                for(int j=0; j<3; j++) {
                    hessian[indexColumn(np,3*m+i,3*n+j)] = secondPartialDPair(m, n, i, j, p, np, nPart, partBoundX, partBoundY, partBoundZ);
                    hessian[indexColumn(np,3*n+j,3*m+i)] = hessian[indexColumn(np,3*m+i,3*n+j)];
                }
        }
    }
}

void solveHessian(particle p[], int np, int nPart[], double partBoundX[], double partBoundY[], double partBoundZ[], int *info) {
    int nn=0;
    int pi[np];
    for(int i=0; i<np; i++)
        for(int j=i+1; j<np; j++)
            if(overlapQ(&(p[i]), &(p[j]))){
                if(nn==0){
                    pi[0]=i;
                    pi[1]=j;
                    nn=2;
                    continue;
                }
                int k;
                for(k=0; k<nn; k++)
                    if(pi[k]==i)
                        break;
                if(k==nn){
                    pi[k]=i;
                    nn++;
                }
                for(k=0; k<nn; k++)
                    if(pi[k]==j)
                        break;
                if(k==nn){
                    pi[k]=j;
                    nn++;
                }
            }
    particle pp[nn];
    for(int i=0; i<nn; i++)
        pp[i] = p[pi[i]];
    
    double *hessian = malloc(3*nn*3*nn*sizeof(double));
    for(int i=0; i<3*nn*3*nn; i++)
        hessian[i] = 0;
    hessianMatrix(pp, nn, nPart, partBoundX, partBoundY, partBoundZ, hessian);
    /*FILE *hes=fopen("hessian.txt","w");
     for(int i=0; i<2*nn*2*nn; i++)
     fprintf(hes,"%lf\n",hessian[i]);
     fclose(hes);*/
    
    for(int j=0; j<np; j++)
        for(int k=0; k<3; k++)
            p[j].force[k] = 0;
    addRepulsiveSpringForce(p, np);
    double forces[3*nn];
    for(int i=0; i<nn; i++)
        for(int j=0; j<3; j++)
            forces[3*i+j] = p[pi[i]].force[j];
    
    /*FILE *force=fopen("force.txt","w");
     for(int i=0; i<2*nn; i++)
     fprintf(force,"%.15lf\n",forces[i]);
     fclose(force);*/
    
    /*double oldForces[3*nn];
    for(int i=0; i<3*nn; i++)
        oldForces[i]=forces[i];
    double *oldHessian = malloc(3*nn*3*nn*sizeof(double));
    for(int i=0; i<3*nn*3*nn; i++)
        oldHessian[i] = 0;
    for(int i=0; i<3*nn*3*nn; i++)
        oldHessian[i]=hessian[i];*/
    
    int n = 3*nn;
    int nrhs = 1;
    int lda = n;
    int ldb = n;
    int ipiv[n];
    dgesv_(&n,&nrhs,hessian,&lda,ipiv,forces,&ldb,info);
    /*for(int i=0; i<2*nn; i++)
     printf("%d: %.15lf\n", i, forces[i]);*/
    
    /*for(int i=0; i<3*nn; i++){
        double sum=0;
        for(int j=0; j<3*nn; j++)
            sum += oldHessian[i*3*nn+j]*forces[j];
        printf("%.15lf\n", sum-oldForces[i]);
     }*/
    
    for(int j=0; j<np; j++)
        for(int k=0; k<3; k++)
            p[j].force[k] = 0;
    for(int i=0; i<nn; i++)
        for(int j=0; j<3; j++)
            p[pi[i]].force[j]=forces[3*i+j];
    
    //free(oldHessian);
    free(hessian);
}

/*
int main() {
    int np;
    double a0, b0;
    double rad;
    FILE *npFile=NULL;
    npFile = fopen("npts.dat","r");
    if (npFile) {
        fscanf(npFile, "%i", &np);
        fclose(npFile);
    } else {
        printf("npFile pointer is null\n");
        exit(1);
    }
    FILE *abInFile=NULL;
    abInFile = fopen("ellipsoidParamsInit.dat","r");
    if (abInFile) {
        fscanf(abInFile, "%lf %lf", &a0, &b0);
        fclose(abInFile);
    } else {
        printf("abInFile pointer is null\n");
        exit(1);
    }
    FILE *radFile=NULL;
    radFile = fopen("rad.dat","r");
    if (radFile) {
        fscanf(radFile, "%lf", &rad);
        fclose(radFile);
    } else {
        printf("radFile pointer is null\n");
        exit(1);
    }
    
    particle p[np];
    
    int nPart[3];
    nPart[0] = (int)(a0/(1.1*rad));
    nPart[2] = (int)(a0/(1.1*rad));
    nPart[1] = (int)(a0/(1.1*rad));
    double partBoundX[nPart[0]];
    double partBoundY[nPart[1]];
    double partBoundZ[nPart[2]];
    for(int i=0; i<nPart[0]; i++)
        partBoundX[i] = -a0+i*2*1.1*rad;
    for(int i=0; i<nPart[1]; i++)
        partBoundY[i] = -a0+i*2*1.1*rad;
    for(int i=0; i<nPart[2]; i++)
        partBoundZ[i] = -a0+i*2*1.1*rad;
    
    FILE *configurationTemp=NULL;
    FILE *radiusFile=NULL;
    configurationTemp = fopen("configurationTemp.asc","r");
    radiusFile = fopen("radii.dat", "r");
    if (configurationTemp&&radiusFile) {
        for (int j = 0; j < np; j++) {
            fscanf(configurationTemp, "%lf %lf %lf", &(p[j].position[0]), &(p[j].position[1]), &(p[j].position[2]));
            fscanf(radiusFile, "%lf", &(p[j].rad));
            partitionOneParticle(&(p[j]), nPart,partBoundX, partBoundY, partBoundZ);
            for(int k=0; k<3; k++) {
                p[j].postPreRelaxPosition[k] = p[j].position[k];
                p[j].postPreRelaxCoord[k] = p[j].coord[k];
            }
        }
        fclose(configurationTemp);
        fclose(radiusFile);
    }
    else {
        printf("configuration or radii pointer is null\n");
        exit(1);
    }
    
    double *hessian = malloc(3*np*3*np*sizeof(double));
    for(int i=0; i<3*np*3*np; i++)
        hessian[i] = 0;
    hessianMatrix(p, np, nPart, partBoundX, partBoundY, partBoundZ, hessian);

    for(int i=0; i<np; i++)
        for(int j=0; j<3; j++)
            p[i].force[j] = 0;
    addRepulsiveSpringForce(p, np);
    double oldForces[3*np];
    for(int i=0; i<np; i++)
        for(int j=0; j<3; j++)
            oldForces[3*i+j] = p[i].force[j];
    
    //for(int i=0; i<np; i++)
    //    printf("%d: %.12lf %.12lf %.12lf\n", i, p[i].force[0], p[i].force[1], p[i].force[2]);

    int info;
    solveHessian(p, np, nPart, partBoundX, partBoundY, partBoundZ, &info);
    
    double forces[3*np];
    for(int i=0; i<np; i++)
        for(int j=0; j<3; j++)
            forces[3*i+j] = p[i].force[j];
    for(int i=0; i<3*np; i++){
        double sum=0;
        for(int j=0; j<3*np; j++)
            sum += hessian[i*3*np+j]*forces[j];
        printf("%.20lf\n", sum-oldForces[i]);
    }
    
    for(int i=0; i<np; i++)
        printf("%d: %.12lf %.12lf %.12lf\n", i, p[i].force[0], p[i].force[1], p[i].force[2]);
        
}
*/
