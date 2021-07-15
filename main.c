//
//  main.c
//  softBidispersityOnEllipsoid
//
//  Created by Zhaoyu Xie on 2/24/20.
//  Copyright © 2020 Zhaoyu Xie. All rights reserved.
//

#include <stdio.h>
#include <time.h>
#include <getopt.h>
#include <stdlib.h>
#include <math.h>
#include "particle.h"
#include "surface.h"
#include "operation.h"
#include "mt19937ar.h"

int main(int argc, char * argv[]) {
    
    int np = 100;
    double rad = 0.1;
    double radRatio = 1; //large to small
    double mixRatio = 0.5;
    double a0=1;//long axis
    double b0=1;//short axis
    double diffusionTimescale = 0.25;
    double expDecayFrac = 0.9;
    int outputSteps = 1;
    
    int animationQ = FALSE;
    int animationRate = 100;
    int readUnfinishedQ = FALSE;
    int energyRecordQ = TRUE;
    int dynamicsQ = FALSE;
    
    static struct option longOptions[]={
        {"radius",              required_argument,  NULL, 'r'},
        {"radRatio",            required_argument,  NULL, 'R'},
        {"particleNumber",      required_argument,  NULL, 'n'},
        {"mixRatio",            required_argument,  NULL, 'i'},
        {"major axis",          required_argument,  NULL, 'a'},
        {"minor axis",          required_argument,  NULL, 'b'},
        {"animate",             no_argument,        NULL, 'm'},
        {"readUnfinished",      no_argument,        NULL, 'u'},
        {"dynamics",            no_argument,        NULL, 'd'},
        {0,                     0,                  0,     0 }
    };
    int optIndex = 0;
    int opt;
    while ((opt = getopt_long(argc, argv,"r:R:n:i:a:b:mud",
                              longOptions, &optIndex )) != -1) {
        switch (opt) {
            case 0:
                break;
            case 'r' :
                rad = atof(optarg);
                break;
            case 'R' :
                radRatio = atof(optarg);
                break;
            case 'n' :
                np = atoi(optarg);
                break;
            case 'i' :
                mixRatio = atof(optarg);
                break;
            case 'a' :
                a0 = atof(optarg);
                break;
            case 'b' :
                b0 = atof(optarg);
                break;
            case 'm' :
                animationQ = TRUE;
                break;
            case 'u' :
                readUnfinishedQ = TRUE;
                break;
            case 'd' :
                dynamicsQ = TRUE;
                break;
            default:
                exit(1);
        }
    }
    
    double a = a0;
    double b = b0;
    
    if(readUnfinishedQ) {
        FILE *radFile=NULL;
        radFile = fopen("rad.dat","r");
        if (radFile) {
            fscanf(radFile, "%lf", &rad);
            fclose(radFile);
        } else {
            printf("radFile pointer is null\n");
            exit(1);
        }
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
        FILE *abTempFile=NULL;
        abTempFile = fopen("ellipsoidParamsTemp.dat","r");
        if (abTempFile) {
            fscanf(abTempFile, "%lf %lf", &a, &b);
            fclose(abTempFile);
        } else {
            printf("abTempFile pointer is null\n");
            exit(1);
        }
    }
    else {
        FILE *radFile=NULL;
        radFile = fopen("rad.dat","w");
        if (radFile) {
            fprintf(radFile, "%.15lf\n", rad);
            fclose(radFile);
        } else {
            printf("abFile pointer is null\n");
            exit(1);
        }
        FILE *npFile=NULL;
        npFile = fopen("npts.dat","w");
        if (npFile) {
            fprintf(npFile, "%i\n", np);
            fclose(npFile);
        } else {
            printf("npFile pointer is null\n");
            exit(1);
        }
        FILE *abInFile=NULL;
        abInFile = fopen("ellipsoidParamsInit.dat","w");
        if (abInFile) {
            fprintf(abInFile, "%.15lf %.15lf\n", a0, b0);
            fclose(abInFile);
        } else {
            printf("abInFile pointer is null\n");
            exit(1);
        }
        FILE *abTempFile=NULL;
        abTempFile = fopen("ellipsoidParamsTemp.dat","w");
        if (abTempFile) {
            fprintf(abTempFile, "%.15lf %.15lf\n", a, b);
            fclose(abTempFile);
        } else {
            printf("abTempFile pointer is null\n");
            exit(1);
        }
    }
    
    double D = pow(2*rad/1000, 2)/2; //diffusion coefficient along long axis
    double delta = sqrt(2*D);
    double simTimeFinal = 2*rad*rad/(D*diffusionTimescale);
    double relaxationSteps = 1e5;
    double dtRelaxationMax = simTimeFinal/relaxationSteps;
    double dtTol = 1e-5*dtRelaxationMax;
    double dtRelaxation;
    
    double simTimeStart = 0;
    double nextRelaxationTime = 0;
    int relaxationStep = 0;
    
    if(readUnfinishedQ) {
        FILE *simTimeTempFile=NULL;
        simTimeTempFile = fopen("simTimeTemp.dat","r");
        if (simTimeTempFile) {
            fscanf(simTimeTempFile, "%lf", &simTimeStart);
            fclose(simTimeTempFile);
        } else {
            printf("simTimeTempFile pointer is null\n");
            exit(1);
        }
        FILE *nextRelaxationFile=NULL;
        nextRelaxationFile = fopen("nextRelaxationTime.dat","r");
        if (nextRelaxationFile) {
            fscanf(nextRelaxationFile, "%lf", &nextRelaxationTime);
            fclose(nextRelaxationFile);
        } else {
            printf("nextRelaxationFile pointer is null\n");
            exit(1);
        }
        FILE *dtRelaxationFile=NULL;
        dtRelaxationFile = fopen("dtRelaxationTime.dat", "r");
        if(dtRelaxationFile) {
            fscanf(dtRelaxationFile, "%lf", &dtRelaxation);
            fclose(dtRelaxationFile);
        } else {
            printf("dtRelaxationFile pointer is null\n");
            exit(1);
        }
        FILE *steps = NULL;
        steps = fopen("relaxationStepsTemp.dat","r");
        if(steps) {
            fscanf(steps, "%d", &relaxationStep);
            fclose(steps);
        }
        else {
            printf("steps pointer is null\n");
            exit(1);
        }
    }
    else {
        dtRelaxation = dtRelaxationMax;
        
        FILE *simTimeTempFile=NULL;
        simTimeTempFile = fopen("simTimeTemp.dat","w");
        if (simTimeTempFile) {
            fprintf(simTimeTempFile, "%.15lf\n", simTimeStart);
            fclose(simTimeTempFile);
        } else {
            printf("simTimeTempFile pointer is null\n");
            exit(1);
        }
        FILE *nextRelaxationFile=NULL;
        nextRelaxationFile = fopen("nextRelaxationTime.dat","w");
        if (nextRelaxationFile) {
            fprintf(nextRelaxationFile, "%.15lf\n", nextRelaxationTime);
            fclose(nextRelaxationFile);
        } else {
            printf("nextRelaxationFile pointer is null\n");
            exit(1);
        }
        FILE *dtRelaxationFile=NULL;
        dtRelaxationFile = fopen("dtRelaxationTime.dat", "w");
        if(dtRelaxationFile) {
            fprintf(dtRelaxationFile, "%.15lf\n", dtRelaxation);
            fclose(dtRelaxationFile);
        } else {
            printf("dtRelaxationFile pointer is null\n");
            exit(1);
        }
        FILE *steps = NULL;
        steps = fopen("relaxationStepsTemp.dat","w");
        if(steps) {
            fprintf(steps, "%d\n", relaxationStep);
            fclose(steps);
        }
        else {
            printf("steps pointer is null\n");
            exit(1);
        }
    }
    
    double simTimeLastRelax = simTimeStart;
    double nextRelaxationTimeLastRelax = nextRelaxationTime;
    nextRelaxationTime += dtRelaxation;
    
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
    
    particle p[np];
    
    randinitialize();
    
    if(readUnfinishedQ){
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
    }
    else {
        for (int i=0; i< (int)(np*mixRatio); i++) {
            p[i].rad=rad;
        }
        for (int i=(int)(np*mixRatio); i<np; i++) {
            p[i].rad=rad/radRatio;
        }
        FILE *radiusFile=NULL;
        radiusFile = fopen("radii.dat","w");
        if (radiusFile) {
            for(int i=0; i<np; i++)
                fprintf(radiusFile, "%lf\n", p[i].rad);
            fclose(radiusFile);
        } else {
            printf("radiusFile pointer is null\n");
            exit(1);
        }
        
        for(int i=0; i<np; i++) {
            double u = conformalInvert(a,b,genrand_real1());
            double v = 2*PI*genrand_real2();
            p[i].position[0]=b*sqrt(1-u*u)*cos(v);
            p[i].position[1]=b*sqrt(1-u*u)*sin(v);
            p[i].position[2]=a*u;
            partitionOneParticle(&(p[i]), nPart, partBoundX, partBoundY, partBoundZ);            
        }
        int nearJammed = FALSE;
        undoOverlaps(p, np, a, b, nPart, partBoundX, partBoundY, partBoundZ, &nearJammed);
        if (nearJammed) {
            printf("overcrowded starting point - use fewer particles or smaller radius\n");
            exit(1);
        }
        FILE *init=NULL; // after projection
        init = fopen("init.asc","w");
        if (!init) {
            printf("file init.asc failed to open");
            exit(1);
        }
        for(int i=0; i<np; i++)
            for(int k=0; k<3; k++) {
                p[i].postPreRelaxPosition[k] = p[i].position[k];
                p[i].postPreRelaxCoord[k] = p[i].coord[k];
                fprintf(init, "%.15lf %.15lf %.15lf\n", p[i].position[0], p[i].position[1], p[i].position[2]);
            }
        if(init) fclose(init);
        
        FILE *configurationTemp=NULL;
        configurationTemp = fopen("configurationTemp.asc","w");
        if (configurationTemp) {
            for (int j = 0; j < np; j++) {
                fprintf(configurationTemp, "%.15lf %.15lf %.15lf\n", p[j].position[0], p[j].position[1], p[j].position[2]);
            }
            fclose(configurationTemp);
        }
        else {
            printf("configuration pointer is null\n");
            exit(1);
        }
    }
    
    FILE *animationFile=NULL;
    if (animationQ && relaxationStep%animationRate==0) {
        animationFile = fopen("animation.dat", "a");
        if (animationFile) {
            fprintf(animationFile, "%.15lf %.15lf\n", a, b);
            for (int j=0; j<np; j++)
                fprintf(animationFile, "%.15lf %.15lf %.15lf\n", p[j].position[0], p[j].position[1], p[j].position[2]);
            fclose(animationFile);
        }
        else {
            printf("animationFile pointer is null\n");
            exit(1);
        }
    }
    
    double Vmin = 1e-16;
    
    FILE *energyFile = NULL;
    
    double configEnergy = totalRepulsiveSpringEnergy(p, np);
    double phi = packingFraction(p, np, a, b);
    int rollbackQ = FALSE;
    double aOld;
    double bOld;
    double dtRelaxationOld = 0;
    double simTime;
    for(simTime=simTimeStart; simTime<=simTimeFinal; simTime += dtRelaxation) {
        if(simTime >= nextRelaxationTime) {
            aOld = a;
            bOld = b;
            if(dynamicsQ) {
                double area = areaAtTime(a0, b0, expDecayFrac, simTime/simTimeFinal);
                double evo = evoParameterAtSurfaceArea(a0, b0, area);
                ellipsoidUpdate(evo, a0, b0, &a, &b);
            }
            else {
                a = a0*(1-simTime/simTimeFinal);
                b = b0*(1-simTime/simTimeFinal);
            }
            for(int i=0; i<np; i++) {
                if(dynamicsQ)
                    projectOntoNewSurface(&(p[i]), aOld, bOld, a, b);
                else
                    scaleOntoNewSurface(&(p[i]), aOld, bOld, a, b);
                partitionOneParticle(&(p[i]), nPart, partBoundX, partBoundY, partBoundZ);
            }
/*
            FILE *abFileTest=NULL;
            abFileTest = fopen("ellipsoidParamsTest.dat","w");
            fprintf(abFileTest, "%.15lf %.15lf\n", a, b);
            fclose(abFileTest);
            FILE *configFileBefore=NULL;
            configFileBefore = fopen("configurationBefore.asc","w");
            for (int j = 0; j < np; j++) {
                fprintf(configFileBefore, "%.15lf %.15lf %.15lf\n", p[j].position[0], p[j].position[1], p[j].position[2]);
            }
            fclose(configFileBefore);
*/
/*            for(int i=0; i<np; i++){
                p[i].rr = p[i].rad;
                p[i].rad = 1.01*p[i].rad;
            }
            conjGradDescentConstrainedIntegrateMinimizeSteps(p, np, rad, a, b, nPart, partBoundX, partBoundY, partBoundZ,1000);
            for(int i=0; i<np; i++)
                p[i].rad = p[i].rr;
*/
            conjGradDescentConstrainedIntegrateMinimizeSteps(p, np, rad, a, b, nPart, partBoundX, partBoundY, partBoundZ, 1e6);
            
            configEnergy = totalRepulsiveSpringEnergy(p, np);
/*
            FILE *configFileAfter=NULL;
            configFileAfter = fopen("configurationAfter.asc","w");
            for (int j = 0; j < np; j++) {
                fprintf(configFileAfter, "%.15lf %.15lf %.15lf\n", p[j].position[0], p[j].position[1], p[j].position[2]);
            }
            fclose(configFileAfter);
*/
/*
            for(int i=0; i<np; i++)
                for(int j=0; j<3; j++) {
                    p[i].positionTemp[j] = p[i].position[j];
                    p[i].coordTemp[j] = p[i].coord[j];
                }
            
            if(configEnergy/np>=Vmin) {
                for(int i=0; i<np; i++){
                    p[i].rr = p[i].rad;
                    p[i].rad = 1.01*p[i].rad;
                }
                conjGradDescentConstrainedIntegrateMinimizeSteps(p, np, rad, a, b, nPart, partBoundX, partBoundY, partBoundZ, 1e7);
                for(int i=0; i<np; i++)
                    p[i].rad = p[i].rr;
                conjGradDescentConstrainedIntegrateMinimizeSteps(p, np, rad, a, b, nPart, partBoundX, partBoundY, partBoundZ, 1e7);
                double configEnergyTemp = totalRepulsiveSpringEnergy(p, np);
                if(configEnergyTemp>configEnergy)
                    for(int i=0; i<np; i++)
                        for(int j=0; j<3; j++) {
                            p[i].position[j] = p[i].positionTemp[j];
                            p[i].coord[j] = p[i].coordTemp[j];
                        }
                else
                    configEnergy = configEnergyTemp;
            }
*/
            if(energyRecordQ){
                energyFile = fopen("energy.dat","a");
                if(energyFile) {
                    fprintf(energyFile, "%.15lf %.20lf\n", simTime, configEnergy);
                    fclose(energyFile);
                }
                else {
                    printf("energyFile pointer is null\n");
                    exit(1);
                }
            }
            
            if(configEnergy/np<Vmin&&aOld!=a) {
                for (int j=0; j<np; j++)
                    for(int k=0; k<3; k++) {
                        p[j].postPreRelaxPosition[k] = p[j].position[k];
                        p[j].postPreRelaxCoord[k] = p[j].coord[k];
                    }
                simTimeLastRelax = simTime;
                nextRelaxationTimeLastRelax = nextRelaxationTime;
                relaxationStep++;
                //dtRelaxation = dtRelaxation*2;
                if(rollbackQ) {
                    dtRelaxation = dtRelaxation/2;
                    rollbackQ = FALSE;
                    dtRelaxationOld = 0;
                }
                else {
                    if(dtRelaxation == dtRelaxationOld)
                        dtRelaxation = dtRelaxation*2;
                    else
                        dtRelaxationOld = dtRelaxation;
                }
                if(dtRelaxation>dtRelaxationMax) {
                    dtRelaxation = dtRelaxationMax;
                    dtRelaxationOld = dtRelaxation;
                }
                
                if (animationQ && relaxationStep%animationRate==0) {
                    animationFile = fopen("animation.dat", "a");
                    if (animationFile) {
                        fprintf(animationFile, "%.15lf %.15lf\n", a, b);
                        for (int j=0; j<np; j++)
                            fprintf(animationFile, "%.15lf %.15lf %.15lf\n", p[j].position[0], p[j].position[1], p[j].position[2]);
                        fclose(animationFile);
                    }
                    else {
                        printf("animationFile pointer is null\n");
                        exit(1);
                    }
                }
                
                if (relaxationStep%outputSteps==0) {
                    FILE *configurationTemp=NULL;
                    configurationTemp = fopen("configurationTemp.asc","w");
                    if (configurationTemp) {
                        for (int j = 0; j < np; j++)
                            fprintf(configurationTemp, "%.15lf %.15lf %.15lf\n", p[j].position[0], p[j].position[1], p[j].position[2]);
                        fclose(configurationTemp);
                    }
                    else {
                        printf("configuration pointer is null\n");
                        exit(1);
                    }
                    FILE *abFileTemp=NULL;
                    abFileTemp = fopen("ellipsoidParamsTemp.dat","w");
                    if (abFileTemp) {
                        fprintf(abFileTemp, "%.15lf %.15lf\n", a, b);
                        fclose(abFileTemp);
                    } else {
                        printf("abFile pointer is null\n");
                        exit(1);
                    }
                    FILE *simTimeTempFile=NULL;
                    simTimeTempFile = fopen("simTimeTemp.dat","w");
                    if (simTimeTempFile) {
                        fprintf(simTimeTempFile, "%.15lf\n", simTime);
                        fclose(simTimeTempFile);
                    } else {
                        printf("simTimeTempFile pointer is null\n");
                        exit(1);
                    }
                    FILE *nextRelaxationFile=NULL;
                    nextRelaxationFile = fopen("nextRelaxationTime.dat","w");
                    if (nextRelaxationFile) {
                        fprintf(nextRelaxationFile, "%.15lf\n", nextRelaxationTime);
                        fclose(nextRelaxationFile);
                    } else {
                        printf("nextRelaxationFile pointer is null\n");
                        exit(1);
                    }
                    FILE *stepsTemp = NULL;
                    stepsTemp = fopen("relaxationStepsTemp.dat","w");
                    if(stepsTemp) {
                        fprintf(stepsTemp, "%d\n", relaxationStep);
                        fclose(stepsTemp);
                    }
                    else {
                        printf("stepsTemp pointer is null\n");
                        exit(1);
                    }
                }
            }
            else {
                if(configEnergy/np<2*Vmin) {
                    relaxationStep++;
                    break;
                }
                rollbackQ = TRUE;
                dtRelaxation = dtRelaxation/2;
                a = aOld;
                b = bOld;
                for (int j=0; j<np; j++)
                    for(int k=0; k<3; k++) {
                        p[j].position[k] = p[j].postPreRelaxPosition[k];
                        p[j].coord[k] = p[j].postPreRelaxCoord[k];
                    }
                simTime = simTimeLastRelax;
                nextRelaxationTime = nextRelaxationTimeLastRelax;
            }
            if (relaxationStep%outputSteps==0) {
                FILE *dtTempFile=NULL;
                dtTempFile = fopen("dtRelaxationTime.dat","w");
                if (dtTempFile) {
                    fprintf(dtTempFile, "%.15lf\n", dtRelaxation);
                    fclose(dtTempFile);
                } else {
                    printf("dtTempFile pointer is null\n");
                    exit(1);
                }
            }
            nextRelaxationTime += dtRelaxation;
        }
        printf("simTime/simTimeFinal: %.16lf\n", simTime/simTimeFinal);
    }
    
    if (animationQ && relaxationStep%animationRate==0) {
        animationFile = fopen("animation.dat", "a");
        if (animationFile) {
            fprintf(animationFile, "%.15lf %.15lf\n", a, b);
            for (int j=0; j<np; j++)
                fprintf(animationFile, "%.15lf %.15lf %.15lf\n", p[j].position[0], p[j].position[1], p[j].position[2]);
            fclose(animationFile);
        }
        else {
            printf("animationFile pointer is null\n");
            exit(1);
        }
    }
    
    FILE *abFile=NULL;
    abFile = fopen("ellipsoidParams.dat","w");
    if (abFile) {
        fprintf(abFile, "%.15lf %.15lf\n", a, b);
        fclose(abFile);
    } else {
        printf("abFile pointer is null\n");
        exit(1);
    }
    FILE *configuration=NULL;
    configuration = fopen("configuration.asc","w");
    if (configuration) {
        for (int j = 0; j < np; j++) {
            fprintf(configuration, "%.15lf %.15lf %.15lf\n", p[j].position[0], p[j].position[1], p[j].position[2]);
        }
        fclose(configuration);
    }
    else {
        printf("configuration pointer is null\n");
        exit(1);
    }
    FILE *simTimeFile=NULL;
    simTimeFile = fopen("simTimeArrested.dat","w");
    if (simTimeFile) {
        fprintf(simTimeFile, "%.15lf\n", simTime);
        fclose(simTimeFile);
    } else {
        printf("simTimeTempFile pointer is null\n");
        exit(1);
    }
    FILE *nextRelaxationFile=NULL;
    nextRelaxationFile = fopen("nextRelaxationTimeArrested.dat","w");
    if (nextRelaxationFile) {
        fprintf(nextRelaxationFile, "%.15lf\n", nextRelaxationTime);
        fclose(nextRelaxationFile);
    } else {
        printf("nextRelaxationFile pointer is null\n");
        exit(1);
    }
    FILE *steps = NULL;
    steps = fopen("relaxationSteps.dat","w");
    if(steps) {
        fprintf(steps, "%d\n", relaxationStep);
        fclose(steps);
    }
    else {
        printf("steps pointer is null\n");
        exit(1);
    }
    FILE *dtFile=NULL;
    dtFile = fopen("dtRelaxationTimeArrested.dat","w");
    if (dtFile) {
        fprintf(dtFile, "%.15lf\n", dtRelaxation);
        fclose(dtFile);
    } else {
        printf("dtTempFile pointer is null\n");
        exit(1);
    }
    
    printf("Hello, World!\n");
    return 0;
}
