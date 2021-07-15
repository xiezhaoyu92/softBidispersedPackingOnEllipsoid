//
//  hessian.h
//  softBidispersityOnEllipsoid
//
//  Created by Zhaoyu Xie on 6/6/20.
//  Copyright © 2020 Zhaoyu Xie. All rights reserved.
//

#ifndef hessian_h
#define hessian_h

#include <stdio.h>
#include "particle.h"

void solveHessian(particle p[], int np, int nPart[], double partBoundX[], double partBoundY[], double partBoundZ[], int *info);

#endif /* hessian_h */
