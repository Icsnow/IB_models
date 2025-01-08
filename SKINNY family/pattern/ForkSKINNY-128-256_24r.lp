\ Model ./ForkSKINNY-128-256_24r
\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize multi-objectives
  Min_Tc: Priority=4 Weight=1 AbsTol=1e-06 RelTol=0
   Tc
  Min_T32: Priority=3 Weight=1 AbsTol=1e-06 RelTol=0
   T32
  Min_Dc: Priority=2 Weight=1 AbsTol=1e-06 RelTol=0
   Dc
  Min_Mc: Priority=1 Weight=1 AbsTol=1e-06 RelTol=0
   Mc
Subject To
 uBasic: uDX[1,0] - udx[1,0] >= 0
 uBasic: uDX[2,0] - udx[2,0] >= 0
 uBasic: uDX[3,0] - udx[3,0] >= 0
 uBasic: uDX[4,0] - udx[4,0] >= 0
 uBasic: uDX[5,0] - udx[5,0] >= 0
 uBasic: uDX[6,0] - udx[6,0] >= 0
 uBasic: uDX[7,0] - udx[7,0] >= 0
 uBasic: uDX[8,0] - udx[8,0] >= 0
 uBasic: uDX[9,0] - udx[9,0] >= 0
 uBasic: uDX[10,0] - udx[10,0] >= 0
 uBasic: uDX[11,0] - udx[11,0] >= 0
 uBasic: uDX[12,0] - udx[12,0] >= 0
 uBasic: uDY[1,0] - udy[1,0] >= 0
 uBasic: uDZ[1,0] - udz[1,0] >= 0
 uBasic: uDY[2,0] - udy[2,0] >= 0
 uBasic: uDZ[2,0] - udz[2,0] >= 0
 uBasic: uDY[3,0] - udy[3,0] >= 0
 uBasic: uDZ[3,0] - udz[3,0] >= 0
 uBasic: uDY[4,0] - udy[4,0] >= 0
 uBasic: uDZ[4,0] - udz[4,0] >= 0
 uBasic: uDY[5,0] - udy[5,0] >= 0
 uBasic: uDZ[5,0] - udz[5,0] >= 0
 uBasic: uDY[6,0] - udy[6,0] >= 0
 uBasic: uDZ[6,0] - udz[6,0] >= 0
 uBasic: uDY[7,0] - udy[7,0] >= 0
 uBasic: uDZ[7,0] - udz[7,0] >= 0
 uBasic: uDY[8,0] - udy[8,0] >= 0
 uBasic: uDZ[8,0] - udz[8,0] >= 0
 uBasic: uDY[9,0] - udy[9,0] >= 0
 uBasic: uDZ[9,0] - udz[9,0] >= 0
 uBasic: uDY[10,0] - udy[10,0] >= 0
 uBasic: uDZ[10,0] - udz[10,0] >= 0
 uBasic: uDY[11,0] - udy[11,0] >= 0
 uBasic: uDZ[11,0] - udz[11,0] >= 0
 uBasic: uDW[0,0] - udw[0,0] >= 0
 uBasic: uDW[1,0] - udw[1,0] >= 0
 uBasic: uDW[2,0] - udw[2,0] >= 0
 uBasic: uDW[3,0] - udw[3,0] >= 0
 uBasic: uDW[4,0] - udw[4,0] >= 0
 uBasic: uDW[5,0] - udw[5,0] >= 0
 uBasic: uDW[6,0] - udw[6,0] >= 0
 uBasic: uDW[7,0] - udw[7,0] >= 0
 uBasic: uDW[8,0] - udw[8,0] >= 0
 uBasic: uDW[9,0] - udw[9,0] >= 0
 uBasic: uDW[10,0] - udw[10,0] >= 0
 uBasic: uDW[11,0] - udw[11,0] >= 0
 uBasic: uDX[1,1] - udx[1,1] >= 0
 uBasic: uDX[2,1] - udx[2,1] >= 0
 uBasic: uDX[3,1] - udx[3,1] >= 0
 uBasic: uDX[4,1] - udx[4,1] >= 0
 uBasic: uDX[5,1] - udx[5,1] >= 0
 uBasic: uDX[6,1] - udx[6,1] >= 0
 uBasic: uDX[7,1] - udx[7,1] >= 0
 uBasic: uDX[8,1] - udx[8,1] >= 0
 uBasic: uDX[9,1] - udx[9,1] >= 0
 uBasic: uDX[10,1] - udx[10,1] >= 0
 uBasic: uDX[11,1] - udx[11,1] >= 0
 uBasic: uDX[12,1] - udx[12,1] >= 0
 uBasic: uDY[1,1] - udy[1,1] >= 0
 uBasic: uDZ[1,1] - udz[1,1] >= 0
 uBasic: uDY[2,1] - udy[2,1] >= 0
 uBasic: uDZ[2,1] - udz[2,1] >= 0
 uBasic: uDY[3,1] - udy[3,1] >= 0
 uBasic: uDZ[3,1] - udz[3,1] >= 0
 uBasic: uDY[4,1] - udy[4,1] >= 0
 uBasic: uDZ[4,1] - udz[4,1] >= 0
 uBasic: uDY[5,1] - udy[5,1] >= 0
 uBasic: uDZ[5,1] - udz[5,1] >= 0
 uBasic: uDY[6,1] - udy[6,1] >= 0
 uBasic: uDZ[6,1] - udz[6,1] >= 0
 uBasic: uDY[7,1] - udy[7,1] >= 0
 uBasic: uDZ[7,1] - udz[7,1] >= 0
 uBasic: uDY[8,1] - udy[8,1] >= 0
 uBasic: uDZ[8,1] - udz[8,1] >= 0
 uBasic: uDY[9,1] - udy[9,1] >= 0
 uBasic: uDZ[9,1] - udz[9,1] >= 0
 uBasic: uDY[10,1] - udy[10,1] >= 0
 uBasic: uDZ[10,1] - udz[10,1] >= 0
 uBasic: uDY[11,1] - udy[11,1] >= 0
 uBasic: uDZ[11,1] - udz[11,1] >= 0
 uBasic: uDW[0,1] - udw[0,1] >= 0
 uBasic: uDW[1,1] - udw[1,1] >= 0
 uBasic: uDW[2,1] - udw[2,1] >= 0
 uBasic: uDW[3,1] - udw[3,1] >= 0
 uBasic: uDW[4,1] - udw[4,1] >= 0
 uBasic: uDW[5,1] - udw[5,1] >= 0
 uBasic: uDW[6,1] - udw[6,1] >= 0
 uBasic: uDW[7,1] - udw[7,1] >= 0
 uBasic: uDW[8,1] - udw[8,1] >= 0
 uBasic: uDW[9,1] - udw[9,1] >= 0
 uBasic: uDW[10,1] - udw[10,1] >= 0
 uBasic: uDW[11,1] - udw[11,1] >= 0
 uBasic: uDX[1,2] - udx[1,2] >= 0
 uBasic: uDX[2,2] - udx[2,2] >= 0
 uBasic: uDX[3,2] - udx[3,2] >= 0
 uBasic: uDX[4,2] - udx[4,2] >= 0
 uBasic: uDX[5,2] - udx[5,2] >= 0
 uBasic: uDX[6,2] - udx[6,2] >= 0
 uBasic: uDX[7,2] - udx[7,2] >= 0
 uBasic: uDX[8,2] - udx[8,2] >= 0
 uBasic: uDX[9,2] - udx[9,2] >= 0
 uBasic: uDX[10,2] - udx[10,2] >= 0
 uBasic: uDX[11,2] - udx[11,2] >= 0
 uBasic: uDX[12,2] - udx[12,2] >= 0
 uBasic: uDY[1,2] - udy[1,2] >= 0
 uBasic: uDZ[1,2] - udz[1,2] >= 0
 uBasic: uDY[2,2] - udy[2,2] >= 0
 uBasic: uDZ[2,2] - udz[2,2] >= 0
 uBasic: uDY[3,2] - udy[3,2] >= 0
 uBasic: uDZ[3,2] - udz[3,2] >= 0
 uBasic: uDY[4,2] - udy[4,2] >= 0
 uBasic: uDZ[4,2] - udz[4,2] >= 0
 uBasic: uDY[5,2] - udy[5,2] >= 0
 uBasic: uDZ[5,2] - udz[5,2] >= 0
 uBasic: uDY[6,2] - udy[6,2] >= 0
 uBasic: uDZ[6,2] - udz[6,2] >= 0
 uBasic: uDY[7,2] - udy[7,2] >= 0
 uBasic: uDZ[7,2] - udz[7,2] >= 0
 uBasic: uDY[8,2] - udy[8,2] >= 0
 uBasic: uDZ[8,2] - udz[8,2] >= 0
 uBasic: uDY[9,2] - udy[9,2] >= 0
 uBasic: uDZ[9,2] - udz[9,2] >= 0
 uBasic: uDY[10,2] - udy[10,2] >= 0
 uBasic: uDZ[10,2] - udz[10,2] >= 0
 uBasic: uDY[11,2] - udy[11,2] >= 0
 uBasic: uDZ[11,2] - udz[11,2] >= 0
 uBasic: uDW[0,2] - udw[0,2] >= 0
 uBasic: uDW[1,2] - udw[1,2] >= 0
 uBasic: uDW[2,2] - udw[2,2] >= 0
 uBasic: uDW[3,2] - udw[3,2] >= 0
 uBasic: uDW[4,2] - udw[4,2] >= 0
 uBasic: uDW[5,2] - udw[5,2] >= 0
 uBasic: uDW[6,2] - udw[6,2] >= 0
 uBasic: uDW[7,2] - udw[7,2] >= 0
 uBasic: uDW[8,2] - udw[8,2] >= 0
 uBasic: uDW[9,2] - udw[9,2] >= 0
 uBasic: uDW[10,2] - udw[10,2] >= 0
 uBasic: uDW[11,2] - udw[11,2] >= 0
 uBasic: uDX[1,3] - udx[1,3] >= 0
 uBasic: uDX[2,3] - udx[2,3] >= 0
 uBasic: uDX[3,3] - udx[3,3] >= 0
 uBasic: uDX[4,3] - udx[4,3] >= 0
 uBasic: uDX[5,3] - udx[5,3] >= 0
 uBasic: uDX[6,3] - udx[6,3] >= 0
 uBasic: uDX[7,3] - udx[7,3] >= 0
 uBasic: uDX[8,3] - udx[8,3] >= 0
 uBasic: uDX[9,3] - udx[9,3] >= 0
 uBasic: uDX[10,3] - udx[10,3] >= 0
 uBasic: uDX[11,3] - udx[11,3] >= 0
 uBasic: uDX[12,3] - udx[12,3] >= 0
 uBasic: uDY[1,3] - udy[1,3] >= 0
 uBasic: uDZ[1,3] - udz[1,3] >= 0
 uBasic: uDY[2,3] - udy[2,3] >= 0
 uBasic: uDZ[2,3] - udz[2,3] >= 0
 uBasic: uDY[3,3] - udy[3,3] >= 0
 uBasic: uDZ[3,3] - udz[3,3] >= 0
 uBasic: uDY[4,3] - udy[4,3] >= 0
 uBasic: uDZ[4,3] - udz[4,3] >= 0
 uBasic: uDY[5,3] - udy[5,3] >= 0
 uBasic: uDZ[5,3] - udz[5,3] >= 0
 uBasic: uDY[6,3] - udy[6,3] >= 0
 uBasic: uDZ[6,3] - udz[6,3] >= 0
 uBasic: uDY[7,3] - udy[7,3] >= 0
 uBasic: uDZ[7,3] - udz[7,3] >= 0
 uBasic: uDY[8,3] - udy[8,3] >= 0
 uBasic: uDZ[8,3] - udz[8,3] >= 0
 uBasic: uDY[9,3] - udy[9,3] >= 0
 uBasic: uDZ[9,3] - udz[9,3] >= 0
 uBasic: uDY[10,3] - udy[10,3] >= 0
 uBasic: uDZ[10,3] - udz[10,3] >= 0
 uBasic: uDY[11,3] - udy[11,3] >= 0
 uBasic: uDZ[11,3] - udz[11,3] >= 0
 uBasic: uDW[0,3] - udw[0,3] >= 0
 uBasic: uDW[1,3] - udw[1,3] >= 0
 uBasic: uDW[2,3] - udw[2,3] >= 0
 uBasic: uDW[3,3] - udw[3,3] >= 0
 uBasic: uDW[4,3] - udw[4,3] >= 0
 uBasic: uDW[5,3] - udw[5,3] >= 0
 uBasic: uDW[6,3] - udw[6,3] >= 0
 uBasic: uDW[7,3] - udw[7,3] >= 0
 uBasic: uDW[8,3] - udw[8,3] >= 0
 uBasic: uDW[9,3] - udw[9,3] >= 0
 uBasic: uDW[10,3] - udw[10,3] >= 0
 uBasic: uDW[11,3] - udw[11,3] >= 0
 uBasic: uDX[1,4] - udx[1,4] >= 0
 uBasic: uDX[2,4] - udx[2,4] >= 0
 uBasic: uDX[3,4] - udx[3,4] >= 0
 uBasic: uDX[4,4] - udx[4,4] >= 0
 uBasic: uDX[5,4] - udx[5,4] >= 0
 uBasic: uDX[6,4] - udx[6,4] >= 0
 uBasic: uDX[7,4] - udx[7,4] >= 0
 uBasic: uDX[8,4] - udx[8,4] >= 0
 uBasic: uDX[9,4] - udx[9,4] >= 0
 uBasic: uDX[10,4] - udx[10,4] >= 0
 uBasic: uDX[11,4] - udx[11,4] >= 0
 uBasic: uDX[12,4] - udx[12,4] >= 0
 uBasic: uDY[1,4] - udy[1,4] >= 0
 uBasic: uDZ[1,4] - udz[1,4] >= 0
 uBasic: uDY[2,4] - udy[2,4] >= 0
 uBasic: uDZ[2,4] - udz[2,4] >= 0
 uBasic: uDY[3,4] - udy[3,4] >= 0
 uBasic: uDZ[3,4] - udz[3,4] >= 0
 uBasic: uDY[4,4] - udy[4,4] >= 0
 uBasic: uDZ[4,4] - udz[4,4] >= 0
 uBasic: uDY[5,4] - udy[5,4] >= 0
 uBasic: uDZ[5,4] - udz[5,4] >= 0
 uBasic: uDY[6,4] - udy[6,4] >= 0
 uBasic: uDZ[6,4] - udz[6,4] >= 0
 uBasic: uDY[7,4] - udy[7,4] >= 0
 uBasic: uDZ[7,4] - udz[7,4] >= 0
 uBasic: uDY[8,4] - udy[8,4] >= 0
 uBasic: uDZ[8,4] - udz[8,4] >= 0
 uBasic: uDY[9,4] - udy[9,4] >= 0
 uBasic: uDZ[9,4] - udz[9,4] >= 0
 uBasic: uDY[10,4] - udy[10,4] >= 0
 uBasic: uDZ[10,4] - udz[10,4] >= 0
 uBasic: uDY[11,4] - udy[11,4] >= 0
 uBasic: uDZ[11,4] - udz[11,4] >= 0
 uBasic: uDW[0,4] - udw[0,4] >= 0
 uBasic: uDW[1,4] - udw[1,4] >= 0
 uBasic: uDW[2,4] - udw[2,4] >= 0
 uBasic: uDW[3,4] - udw[3,4] >= 0
 uBasic: uDW[4,4] - udw[4,4] >= 0
 uBasic: uDW[5,4] - udw[5,4] >= 0
 uBasic: uDW[6,4] - udw[6,4] >= 0
 uBasic: uDW[7,4] - udw[7,4] >= 0
 uBasic: uDW[8,4] - udw[8,4] >= 0
 uBasic: uDW[9,4] - udw[9,4] >= 0
 uBasic: uDW[10,4] - udw[10,4] >= 0
 uBasic: uDW[11,4] - udw[11,4] >= 0
 uBasic: uDX[1,5] - udx[1,5] >= 0
 uBasic: uDX[2,5] - udx[2,5] >= 0
 uBasic: uDX[3,5] - udx[3,5] >= 0
 uBasic: uDX[4,5] - udx[4,5] >= 0
 uBasic: uDX[5,5] - udx[5,5] >= 0
 uBasic: uDX[6,5] - udx[6,5] >= 0
 uBasic: uDX[7,5] - udx[7,5] >= 0
 uBasic: uDX[8,5] - udx[8,5] >= 0
 uBasic: uDX[9,5] - udx[9,5] >= 0
 uBasic: uDX[10,5] - udx[10,5] >= 0
 uBasic: uDX[11,5] - udx[11,5] >= 0
 uBasic: uDX[12,5] - udx[12,5] >= 0
 uBasic: uDY[1,5] - udy[1,5] >= 0
 uBasic: uDZ[1,5] - udz[1,5] >= 0
 uBasic: uDY[2,5] - udy[2,5] >= 0
 uBasic: uDZ[2,5] - udz[2,5] >= 0
 uBasic: uDY[3,5] - udy[3,5] >= 0
 uBasic: uDZ[3,5] - udz[3,5] >= 0
 uBasic: uDY[4,5] - udy[4,5] >= 0
 uBasic: uDZ[4,5] - udz[4,5] >= 0
 uBasic: uDY[5,5] - udy[5,5] >= 0
 uBasic: uDZ[5,5] - udz[5,5] >= 0
 uBasic: uDY[6,5] - udy[6,5] >= 0
 uBasic: uDZ[6,5] - udz[6,5] >= 0
 uBasic: uDY[7,5] - udy[7,5] >= 0
 uBasic: uDZ[7,5] - udz[7,5] >= 0
 uBasic: uDY[8,5] - udy[8,5] >= 0
 uBasic: uDZ[8,5] - udz[8,5] >= 0
 uBasic: uDY[9,5] - udy[9,5] >= 0
 uBasic: uDZ[9,5] - udz[9,5] >= 0
 uBasic: uDY[10,5] - udy[10,5] >= 0
 uBasic: uDZ[10,5] - udz[10,5] >= 0
 uBasic: uDY[11,5] - udy[11,5] >= 0
 uBasic: uDZ[11,5] - udz[11,5] >= 0
 uBasic: uDW[0,5] - udw[0,5] >= 0
 uBasic: uDW[1,5] - udw[1,5] >= 0
 uBasic: uDW[2,5] - udw[2,5] >= 0
 uBasic: uDW[3,5] - udw[3,5] >= 0
 uBasic: uDW[4,5] - udw[4,5] >= 0
 uBasic: uDW[5,5] - udw[5,5] >= 0
 uBasic: uDW[6,5] - udw[6,5] >= 0
 uBasic: uDW[7,5] - udw[7,5] >= 0
 uBasic: uDW[8,5] - udw[8,5] >= 0
 uBasic: uDW[9,5] - udw[9,5] >= 0
 uBasic: uDW[10,5] - udw[10,5] >= 0
 uBasic: uDW[11,5] - udw[11,5] >= 0
 uBasic: uDX[1,6] - udx[1,6] >= 0
 uBasic: uDX[2,6] - udx[2,6] >= 0
 uBasic: uDX[3,6] - udx[3,6] >= 0
 uBasic: uDX[4,6] - udx[4,6] >= 0
 uBasic: uDX[5,6] - udx[5,6] >= 0
 uBasic: uDX[6,6] - udx[6,6] >= 0
 uBasic: uDX[7,6] - udx[7,6] >= 0
 uBasic: uDX[8,6] - udx[8,6] >= 0
 uBasic: uDX[9,6] - udx[9,6] >= 0
 uBasic: uDX[10,6] - udx[10,6] >= 0
 uBasic: uDX[11,6] - udx[11,6] >= 0
 uBasic: uDX[12,6] - udx[12,6] >= 0
 uBasic: uDY[1,6] - udy[1,6] >= 0
 uBasic: uDZ[1,6] - udz[1,6] >= 0
 uBasic: uDY[2,6] - udy[2,6] >= 0
 uBasic: uDZ[2,6] - udz[2,6] >= 0
 uBasic: uDY[3,6] - udy[3,6] >= 0
 uBasic: uDZ[3,6] - udz[3,6] >= 0
 uBasic: uDY[4,6] - udy[4,6] >= 0
 uBasic: uDZ[4,6] - udz[4,6] >= 0
 uBasic: uDY[5,6] - udy[5,6] >= 0
 uBasic: uDZ[5,6] - udz[5,6] >= 0
 uBasic: uDY[6,6] - udy[6,6] >= 0
 uBasic: uDZ[6,6] - udz[6,6] >= 0
 uBasic: uDY[7,6] - udy[7,6] >= 0
 uBasic: uDZ[7,6] - udz[7,6] >= 0
 uBasic: uDY[8,6] - udy[8,6] >= 0
 uBasic: uDZ[8,6] - udz[8,6] >= 0
 uBasic: uDY[9,6] - udy[9,6] >= 0
 uBasic: uDZ[9,6] - udz[9,6] >= 0
 uBasic: uDY[10,6] - udy[10,6] >= 0
 uBasic: uDZ[10,6] - udz[10,6] >= 0
 uBasic: uDY[11,6] - udy[11,6] >= 0
 uBasic: uDZ[11,6] - udz[11,6] >= 0
 uBasic: uDW[0,6] - udw[0,6] >= 0
 uBasic: uDW[1,6] - udw[1,6] >= 0
 uBasic: uDW[2,6] - udw[2,6] >= 0
 uBasic: uDW[3,6] - udw[3,6] >= 0
 uBasic: uDW[4,6] - udw[4,6] >= 0
 uBasic: uDW[5,6] - udw[5,6] >= 0
 uBasic: uDW[6,6] - udw[6,6] >= 0
 uBasic: uDW[7,6] - udw[7,6] >= 0
 uBasic: uDW[8,6] - udw[8,6] >= 0
 uBasic: uDW[9,6] - udw[9,6] >= 0
 uBasic: uDW[10,6] - udw[10,6] >= 0
 uBasic: uDW[11,6] - udw[11,6] >= 0
 uBasic: uDX[1,7] - udx[1,7] >= 0
 uBasic: uDX[2,7] - udx[2,7] >= 0
 uBasic: uDX[3,7] - udx[3,7] >= 0
 uBasic: uDX[4,7] - udx[4,7] >= 0
 uBasic: uDX[5,7] - udx[5,7] >= 0
 uBasic: uDX[6,7] - udx[6,7] >= 0
 uBasic: uDX[7,7] - udx[7,7] >= 0
 uBasic: uDX[8,7] - udx[8,7] >= 0
 uBasic: uDX[9,7] - udx[9,7] >= 0
 uBasic: uDX[10,7] - udx[10,7] >= 0
 uBasic: uDX[11,7] - udx[11,7] >= 0
 uBasic: uDX[12,7] - udx[12,7] >= 0
 uBasic: uDY[1,7] - udy[1,7] >= 0
 uBasic: uDZ[1,7] - udz[1,7] >= 0
 uBasic: uDY[2,7] - udy[2,7] >= 0
 uBasic: uDZ[2,7] - udz[2,7] >= 0
 uBasic: uDY[3,7] - udy[3,7] >= 0
 uBasic: uDZ[3,7] - udz[3,7] >= 0
 uBasic: uDY[4,7] - udy[4,7] >= 0
 uBasic: uDZ[4,7] - udz[4,7] >= 0
 uBasic: uDY[5,7] - udy[5,7] >= 0
 uBasic: uDZ[5,7] - udz[5,7] >= 0
 uBasic: uDY[6,7] - udy[6,7] >= 0
 uBasic: uDZ[6,7] - udz[6,7] >= 0
 uBasic: uDY[7,7] - udy[7,7] >= 0
 uBasic: uDZ[7,7] - udz[7,7] >= 0
 uBasic: uDY[8,7] - udy[8,7] >= 0
 uBasic: uDZ[8,7] - udz[8,7] >= 0
 uBasic: uDY[9,7] - udy[9,7] >= 0
 uBasic: uDZ[9,7] - udz[9,7] >= 0
 uBasic: uDY[10,7] - udy[10,7] >= 0
 uBasic: uDZ[10,7] - udz[10,7] >= 0
 uBasic: uDY[11,7] - udy[11,7] >= 0
 uBasic: uDZ[11,7] - udz[11,7] >= 0
 uBasic: uDW[0,7] - udw[0,7] >= 0
 uBasic: uDW[1,7] - udw[1,7] >= 0
 uBasic: uDW[2,7] - udw[2,7] >= 0
 uBasic: uDW[3,7] - udw[3,7] >= 0
 uBasic: uDW[4,7] - udw[4,7] >= 0
 uBasic: uDW[5,7] - udw[5,7] >= 0
 uBasic: uDW[6,7] - udw[6,7] >= 0
 uBasic: uDW[7,7] - udw[7,7] >= 0
 uBasic: uDW[8,7] - udw[8,7] >= 0
 uBasic: uDW[9,7] - udw[9,7] >= 0
 uBasic: uDW[10,7] - udw[10,7] >= 0
 uBasic: uDW[11,7] - udw[11,7] >= 0
 uBasic: uDX[1,8] - udx[1,8] >= 0
 uBasic: uDX[2,8] - udx[2,8] >= 0
 uBasic: uDX[3,8] - udx[3,8] >= 0
 uBasic: uDX[4,8] - udx[4,8] >= 0
 uBasic: uDX[5,8] - udx[5,8] >= 0
 uBasic: uDX[6,8] - udx[6,8] >= 0
 uBasic: uDX[7,8] - udx[7,8] >= 0
 uBasic: uDX[8,8] - udx[8,8] >= 0
 uBasic: uDX[9,8] - udx[9,8] >= 0
 uBasic: uDX[10,8] - udx[10,8] >= 0
 uBasic: uDX[11,8] - udx[11,8] >= 0
 uBasic: uDX[12,8] - udx[12,8] >= 0
 uBasic: uDY[1,8] - udy[1,8] >= 0
 uBasic: uDZ[1,8] - udz[1,8] >= 0
 uBasic: uDY[2,8] - udy[2,8] >= 0
 uBasic: uDZ[2,8] - udz[2,8] >= 0
 uBasic: uDY[3,8] - udy[3,8] >= 0
 uBasic: uDZ[3,8] - udz[3,8] >= 0
 uBasic: uDY[4,8] - udy[4,8] >= 0
 uBasic: uDZ[4,8] - udz[4,8] >= 0
 uBasic: uDY[5,8] - udy[5,8] >= 0
 uBasic: uDZ[5,8] - udz[5,8] >= 0
 uBasic: uDY[6,8] - udy[6,8] >= 0
 uBasic: uDZ[6,8] - udz[6,8] >= 0
 uBasic: uDY[7,8] - udy[7,8] >= 0
 uBasic: uDZ[7,8] - udz[7,8] >= 0
 uBasic: uDY[8,8] - udy[8,8] >= 0
 uBasic: uDZ[8,8] - udz[8,8] >= 0
 uBasic: uDY[9,8] - udy[9,8] >= 0
 uBasic: uDZ[9,8] - udz[9,8] >= 0
 uBasic: uDY[10,8] - udy[10,8] >= 0
 uBasic: uDZ[10,8] - udz[10,8] >= 0
 uBasic: uDY[11,8] - udy[11,8] >= 0
 uBasic: uDZ[11,8] - udz[11,8] >= 0
 uBasic: uDW[0,8] - udw[0,8] >= 0
 uBasic: uDW[1,8] - udw[1,8] >= 0
 uBasic: uDW[2,8] - udw[2,8] >= 0
 uBasic: uDW[3,8] - udw[3,8] >= 0
 uBasic: uDW[4,8] - udw[4,8] >= 0
 uBasic: uDW[5,8] - udw[5,8] >= 0
 uBasic: uDW[6,8] - udw[6,8] >= 0
 uBasic: uDW[7,8] - udw[7,8] >= 0
 uBasic: uDW[8,8] - udw[8,8] >= 0
 uBasic: uDW[9,8] - udw[9,8] >= 0
 uBasic: uDW[10,8] - udw[10,8] >= 0
 uBasic: uDW[11,8] - udw[11,8] >= 0
 uBasic: uDX[1,9] - udx[1,9] >= 0
 uBasic: uDX[2,9] - udx[2,9] >= 0
 uBasic: uDX[3,9] - udx[3,9] >= 0
 uBasic: uDX[4,9] - udx[4,9] >= 0
 uBasic: uDX[5,9] - udx[5,9] >= 0
 uBasic: uDX[6,9] - udx[6,9] >= 0
 uBasic: uDX[7,9] - udx[7,9] >= 0
 uBasic: uDX[8,9] - udx[8,9] >= 0
 uBasic: uDX[9,9] - udx[9,9] >= 0
 uBasic: uDX[10,9] - udx[10,9] >= 0
 uBasic: uDX[11,9] - udx[11,9] >= 0
 uBasic: uDX[12,9] - udx[12,9] >= 0
 uBasic: uDY[1,9] - udy[1,9] >= 0
 uBasic: uDZ[1,9] - udz[1,9] >= 0
 uBasic: uDY[2,9] - udy[2,9] >= 0
 uBasic: uDZ[2,9] - udz[2,9] >= 0
 uBasic: uDY[3,9] - udy[3,9] >= 0
 uBasic: uDZ[3,9] - udz[3,9] >= 0
 uBasic: uDY[4,9] - udy[4,9] >= 0
 uBasic: uDZ[4,9] - udz[4,9] >= 0
 uBasic: uDY[5,9] - udy[5,9] >= 0
 uBasic: uDZ[5,9] - udz[5,9] >= 0
 uBasic: uDY[6,9] - udy[6,9] >= 0
 uBasic: uDZ[6,9] - udz[6,9] >= 0
 uBasic: uDY[7,9] - udy[7,9] >= 0
 uBasic: uDZ[7,9] - udz[7,9] >= 0
 uBasic: uDY[8,9] - udy[8,9] >= 0
 uBasic: uDZ[8,9] - udz[8,9] >= 0
 uBasic: uDY[9,9] - udy[9,9] >= 0
 uBasic: uDZ[9,9] - udz[9,9] >= 0
 uBasic: uDY[10,9] - udy[10,9] >= 0
 uBasic: uDZ[10,9] - udz[10,9] >= 0
 uBasic: uDY[11,9] - udy[11,9] >= 0
 uBasic: uDZ[11,9] - udz[11,9] >= 0
 uBasic: uDW[0,9] - udw[0,9] >= 0
 uBasic: uDW[1,9] - udw[1,9] >= 0
 uBasic: uDW[2,9] - udw[2,9] >= 0
 uBasic: uDW[3,9] - udw[3,9] >= 0
 uBasic: uDW[4,9] - udw[4,9] >= 0
 uBasic: uDW[5,9] - udw[5,9] >= 0
 uBasic: uDW[6,9] - udw[6,9] >= 0
 uBasic: uDW[7,9] - udw[7,9] >= 0
 uBasic: uDW[8,9] - udw[8,9] >= 0
 uBasic: uDW[9,9] - udw[9,9] >= 0
 uBasic: uDW[10,9] - udw[10,9] >= 0
 uBasic: uDW[11,9] - udw[11,9] >= 0
 uBasic: uDX[1,10] - udx[1,10] >= 0
 uBasic: uDX[2,10] - udx[2,10] >= 0
 uBasic: uDX[3,10] - udx[3,10] >= 0
 uBasic: uDX[4,10] - udx[4,10] >= 0
 uBasic: uDX[5,10] - udx[5,10] >= 0
 uBasic: uDX[6,10] - udx[6,10] >= 0
 uBasic: uDX[7,10] - udx[7,10] >= 0
 uBasic: uDX[8,10] - udx[8,10] >= 0
 uBasic: uDX[9,10] - udx[9,10] >= 0
 uBasic: uDX[10,10] - udx[10,10] >= 0
 uBasic: uDX[11,10] - udx[11,10] >= 0
 uBasic: uDX[12,10] - udx[12,10] >= 0
 uBasic: uDY[1,10] - udy[1,10] >= 0
 uBasic: uDZ[1,10] - udz[1,10] >= 0
 uBasic: uDY[2,10] - udy[2,10] >= 0
 uBasic: uDZ[2,10] - udz[2,10] >= 0
 uBasic: uDY[3,10] - udy[3,10] >= 0
 uBasic: uDZ[3,10] - udz[3,10] >= 0
 uBasic: uDY[4,10] - udy[4,10] >= 0
 uBasic: uDZ[4,10] - udz[4,10] >= 0
 uBasic: uDY[5,10] - udy[5,10] >= 0
 uBasic: uDZ[5,10] - udz[5,10] >= 0
 uBasic: uDY[6,10] - udy[6,10] >= 0
 uBasic: uDZ[6,10] - udz[6,10] >= 0
 uBasic: uDY[7,10] - udy[7,10] >= 0
 uBasic: uDZ[7,10] - udz[7,10] >= 0
 uBasic: uDY[8,10] - udy[8,10] >= 0
 uBasic: uDZ[8,10] - udz[8,10] >= 0
 uBasic: uDY[9,10] - udy[9,10] >= 0
 uBasic: uDZ[9,10] - udz[9,10] >= 0
 uBasic: uDY[10,10] - udy[10,10] >= 0
 uBasic: uDZ[10,10] - udz[10,10] >= 0
 uBasic: uDY[11,10] - udy[11,10] >= 0
 uBasic: uDZ[11,10] - udz[11,10] >= 0
 uBasic: uDW[0,10] - udw[0,10] >= 0
 uBasic: uDW[1,10] - udw[1,10] >= 0
 uBasic: uDW[2,10] - udw[2,10] >= 0
 uBasic: uDW[3,10] - udw[3,10] >= 0
 uBasic: uDW[4,10] - udw[4,10] >= 0
 uBasic: uDW[5,10] - udw[5,10] >= 0
 uBasic: uDW[6,10] - udw[6,10] >= 0
 uBasic: uDW[7,10] - udw[7,10] >= 0
 uBasic: uDW[8,10] - udw[8,10] >= 0
 uBasic: uDW[9,10] - udw[9,10] >= 0
 uBasic: uDW[10,10] - udw[10,10] >= 0
 uBasic: uDW[11,10] - udw[11,10] >= 0
 uBasic: uDX[1,11] - udx[1,11] >= 0
 uBasic: uDX[2,11] - udx[2,11] >= 0
 uBasic: uDX[3,11] - udx[3,11] >= 0
 uBasic: uDX[4,11] - udx[4,11] >= 0
 uBasic: uDX[5,11] - udx[5,11] >= 0
 uBasic: uDX[6,11] - udx[6,11] >= 0
 uBasic: uDX[7,11] - udx[7,11] >= 0
 uBasic: uDX[8,11] - udx[8,11] >= 0
 uBasic: uDX[9,11] - udx[9,11] >= 0
 uBasic: uDX[10,11] - udx[10,11] >= 0
 uBasic: uDX[11,11] - udx[11,11] >= 0
 uBasic: uDX[12,11] - udx[12,11] >= 0
 uBasic: uDY[1,11] - udy[1,11] >= 0
 uBasic: uDZ[1,11] - udz[1,11] >= 0
 uBasic: uDY[2,11] - udy[2,11] >= 0
 uBasic: uDZ[2,11] - udz[2,11] >= 0
 uBasic: uDY[3,11] - udy[3,11] >= 0
 uBasic: uDZ[3,11] - udz[3,11] >= 0
 uBasic: uDY[4,11] - udy[4,11] >= 0
 uBasic: uDZ[4,11] - udz[4,11] >= 0
 uBasic: uDY[5,11] - udy[5,11] >= 0
 uBasic: uDZ[5,11] - udz[5,11] >= 0
 uBasic: uDY[6,11] - udy[6,11] >= 0
 uBasic: uDZ[6,11] - udz[6,11] >= 0
 uBasic: uDY[7,11] - udy[7,11] >= 0
 uBasic: uDZ[7,11] - udz[7,11] >= 0
 uBasic: uDY[8,11] - udy[8,11] >= 0
 uBasic: uDZ[8,11] - udz[8,11] >= 0
 uBasic: uDY[9,11] - udy[9,11] >= 0
 uBasic: uDZ[9,11] - udz[9,11] >= 0
 uBasic: uDY[10,11] - udy[10,11] >= 0
 uBasic: uDZ[10,11] - udz[10,11] >= 0
 uBasic: uDY[11,11] - udy[11,11] >= 0
 uBasic: uDZ[11,11] - udz[11,11] >= 0
 uBasic: uDW[0,11] - udw[0,11] >= 0
 uBasic: uDW[1,11] - udw[1,11] >= 0
 uBasic: uDW[2,11] - udw[2,11] >= 0
 uBasic: uDW[3,11] - udw[3,11] >= 0
 uBasic: uDW[4,11] - udw[4,11] >= 0
 uBasic: uDW[5,11] - udw[5,11] >= 0
 uBasic: uDW[6,11] - udw[6,11] >= 0
 uBasic: uDW[7,11] - udw[7,11] >= 0
 uBasic: uDW[8,11] - udw[8,11] >= 0
 uBasic: uDW[9,11] - udw[9,11] >= 0
 uBasic: uDW[10,11] - udw[10,11] >= 0
 uBasic: uDW[11,11] - udw[11,11] >= 0
 uBasic: uDX[1,12] - udx[1,12] >= 0
 uBasic: uDX[2,12] - udx[2,12] >= 0
 uBasic: uDX[3,12] - udx[3,12] >= 0
 uBasic: uDX[4,12] - udx[4,12] >= 0
 uBasic: uDX[5,12] - udx[5,12] >= 0
 uBasic: uDX[6,12] - udx[6,12] >= 0
 uBasic: uDX[7,12] - udx[7,12] >= 0
 uBasic: uDX[8,12] - udx[8,12] >= 0
 uBasic: uDX[9,12] - udx[9,12] >= 0
 uBasic: uDX[10,12] - udx[10,12] >= 0
 uBasic: uDX[11,12] - udx[11,12] >= 0
 uBasic: uDX[12,12] - udx[12,12] >= 0
 uBasic: uDY[1,12] - udy[1,12] >= 0
 uBasic: uDZ[1,12] - udz[1,12] >= 0
 uBasic: uDY[2,12] - udy[2,12] >= 0
 uBasic: uDZ[2,12] - udz[2,12] >= 0
 uBasic: uDY[3,12] - udy[3,12] >= 0
 uBasic: uDZ[3,12] - udz[3,12] >= 0
 uBasic: uDY[4,12] - udy[4,12] >= 0
 uBasic: uDZ[4,12] - udz[4,12] >= 0
 uBasic: uDY[5,12] - udy[5,12] >= 0
 uBasic: uDZ[5,12] - udz[5,12] >= 0
 uBasic: uDY[6,12] - udy[6,12] >= 0
 uBasic: uDZ[6,12] - udz[6,12] >= 0
 uBasic: uDY[7,12] - udy[7,12] >= 0
 uBasic: uDZ[7,12] - udz[7,12] >= 0
 uBasic: uDY[8,12] - udy[8,12] >= 0
 uBasic: uDZ[8,12] - udz[8,12] >= 0
 uBasic: uDY[9,12] - udy[9,12] >= 0
 uBasic: uDZ[9,12] - udz[9,12] >= 0
 uBasic: uDY[10,12] - udy[10,12] >= 0
 uBasic: uDZ[10,12] - udz[10,12] >= 0
 uBasic: uDY[11,12] - udy[11,12] >= 0
 uBasic: uDZ[11,12] - udz[11,12] >= 0
 uBasic: uDW[0,12] - udw[0,12] >= 0
 uBasic: uDW[1,12] - udw[1,12] >= 0
 uBasic: uDW[2,12] - udw[2,12] >= 0
 uBasic: uDW[3,12] - udw[3,12] >= 0
 uBasic: uDW[4,12] - udw[4,12] >= 0
 uBasic: uDW[5,12] - udw[5,12] >= 0
 uBasic: uDW[6,12] - udw[6,12] >= 0
 uBasic: uDW[7,12] - udw[7,12] >= 0
 uBasic: uDW[8,12] - udw[8,12] >= 0
 uBasic: uDW[9,12] - udw[9,12] >= 0
 uBasic: uDW[10,12] - udw[10,12] >= 0
 uBasic: uDW[11,12] - udw[11,12] >= 0
 uBasic: uDX[1,13] - udx[1,13] >= 0
 uBasic: uDX[2,13] - udx[2,13] >= 0
 uBasic: uDX[3,13] - udx[3,13] >= 0
 uBasic: uDX[4,13] - udx[4,13] >= 0
 uBasic: uDX[5,13] - udx[5,13] >= 0
 uBasic: uDX[6,13] - udx[6,13] >= 0
 uBasic: uDX[7,13] - udx[7,13] >= 0
 uBasic: uDX[8,13] - udx[8,13] >= 0
 uBasic: uDX[9,13] - udx[9,13] >= 0
 uBasic: uDX[10,13] - udx[10,13] >= 0
 uBasic: uDX[11,13] - udx[11,13] >= 0
 uBasic: uDX[12,13] - udx[12,13] >= 0
 uBasic: uDY[1,13] - udy[1,13] >= 0
 uBasic: uDZ[1,13] - udz[1,13] >= 0
 uBasic: uDY[2,13] - udy[2,13] >= 0
 uBasic: uDZ[2,13] - udz[2,13] >= 0
 uBasic: uDY[3,13] - udy[3,13] >= 0
 uBasic: uDZ[3,13] - udz[3,13] >= 0
 uBasic: uDY[4,13] - udy[4,13] >= 0
 uBasic: uDZ[4,13] - udz[4,13] >= 0
 uBasic: uDY[5,13] - udy[5,13] >= 0
 uBasic: uDZ[5,13] - udz[5,13] >= 0
 uBasic: uDY[6,13] - udy[6,13] >= 0
 uBasic: uDZ[6,13] - udz[6,13] >= 0
 uBasic: uDY[7,13] - udy[7,13] >= 0
 uBasic: uDZ[7,13] - udz[7,13] >= 0
 uBasic: uDY[8,13] - udy[8,13] >= 0
 uBasic: uDZ[8,13] - udz[8,13] >= 0
 uBasic: uDY[9,13] - udy[9,13] >= 0
 uBasic: uDZ[9,13] - udz[9,13] >= 0
 uBasic: uDY[10,13] - udy[10,13] >= 0
 uBasic: uDZ[10,13] - udz[10,13] >= 0
 uBasic: uDY[11,13] - udy[11,13] >= 0
 uBasic: uDZ[11,13] - udz[11,13] >= 0
 uBasic: uDW[0,13] - udw[0,13] >= 0
 uBasic: uDW[1,13] - udw[1,13] >= 0
 uBasic: uDW[2,13] - udw[2,13] >= 0
 uBasic: uDW[3,13] - udw[3,13] >= 0
 uBasic: uDW[4,13] - udw[4,13] >= 0
 uBasic: uDW[5,13] - udw[5,13] >= 0
 uBasic: uDW[6,13] - udw[6,13] >= 0
 uBasic: uDW[7,13] - udw[7,13] >= 0
 uBasic: uDW[8,13] - udw[8,13] >= 0
 uBasic: uDW[9,13] - udw[9,13] >= 0
 uBasic: uDW[10,13] - udw[10,13] >= 0
 uBasic: uDW[11,13] - udw[11,13] >= 0
 uBasic: uDX[1,14] - udx[1,14] >= 0
 uBasic: uDX[2,14] - udx[2,14] >= 0
 uBasic: uDX[3,14] - udx[3,14] >= 0
 uBasic: uDX[4,14] - udx[4,14] >= 0
 uBasic: uDX[5,14] - udx[5,14] >= 0
 uBasic: uDX[6,14] - udx[6,14] >= 0
 uBasic: uDX[7,14] - udx[7,14] >= 0
 uBasic: uDX[8,14] - udx[8,14] >= 0
 uBasic: uDX[9,14] - udx[9,14] >= 0
 uBasic: uDX[10,14] - udx[10,14] >= 0
 uBasic: uDX[11,14] - udx[11,14] >= 0
 uBasic: uDX[12,14] - udx[12,14] >= 0
 uBasic: uDY[1,14] - udy[1,14] >= 0
 uBasic: uDZ[1,14] - udz[1,14] >= 0
 uBasic: uDY[2,14] - udy[2,14] >= 0
 uBasic: uDZ[2,14] - udz[2,14] >= 0
 uBasic: uDY[3,14] - udy[3,14] >= 0
 uBasic: uDZ[3,14] - udz[3,14] >= 0
 uBasic: uDY[4,14] - udy[4,14] >= 0
 uBasic: uDZ[4,14] - udz[4,14] >= 0
 uBasic: uDY[5,14] - udy[5,14] >= 0
 uBasic: uDZ[5,14] - udz[5,14] >= 0
 uBasic: uDY[6,14] - udy[6,14] >= 0
 uBasic: uDZ[6,14] - udz[6,14] >= 0
 uBasic: uDY[7,14] - udy[7,14] >= 0
 uBasic: uDZ[7,14] - udz[7,14] >= 0
 uBasic: uDY[8,14] - udy[8,14] >= 0
 uBasic: uDZ[8,14] - udz[8,14] >= 0
 uBasic: uDY[9,14] - udy[9,14] >= 0
 uBasic: uDZ[9,14] - udz[9,14] >= 0
 uBasic: uDY[10,14] - udy[10,14] >= 0
 uBasic: uDZ[10,14] - udz[10,14] >= 0
 uBasic: uDY[11,14] - udy[11,14] >= 0
 uBasic: uDZ[11,14] - udz[11,14] >= 0
 uBasic: uDW[0,14] - udw[0,14] >= 0
 uBasic: uDW[1,14] - udw[1,14] >= 0
 uBasic: uDW[2,14] - udw[2,14] >= 0
 uBasic: uDW[3,14] - udw[3,14] >= 0
 uBasic: uDW[4,14] - udw[4,14] >= 0
 uBasic: uDW[5,14] - udw[5,14] >= 0
 uBasic: uDW[6,14] - udw[6,14] >= 0
 uBasic: uDW[7,14] - udw[7,14] >= 0
 uBasic: uDW[8,14] - udw[8,14] >= 0
 uBasic: uDW[9,14] - udw[9,14] >= 0
 uBasic: uDW[10,14] - udw[10,14] >= 0
 uBasic: uDW[11,14] - udw[11,14] >= 0
 uBasic: uDX[1,15] - udx[1,15] >= 0
 uBasic: uDX[2,15] - udx[2,15] >= 0
 uBasic: uDX[3,15] - udx[3,15] >= 0
 uBasic: uDX[4,15] - udx[4,15] >= 0
 uBasic: uDX[5,15] - udx[5,15] >= 0
 uBasic: uDX[6,15] - udx[6,15] >= 0
 uBasic: uDX[7,15] - udx[7,15] >= 0
 uBasic: uDX[8,15] - udx[8,15] >= 0
 uBasic: uDX[9,15] - udx[9,15] >= 0
 uBasic: uDX[10,15] - udx[10,15] >= 0
 uBasic: uDX[11,15] - udx[11,15] >= 0
 uBasic: uDX[12,15] - udx[12,15] >= 0
 uBasic: uDY[1,15] - udy[1,15] >= 0
 uBasic: uDZ[1,15] - udz[1,15] >= 0
 uBasic: uDY[2,15] - udy[2,15] >= 0
 uBasic: uDZ[2,15] - udz[2,15] >= 0
 uBasic: uDY[3,15] - udy[3,15] >= 0
 uBasic: uDZ[3,15] - udz[3,15] >= 0
 uBasic: uDY[4,15] - udy[4,15] >= 0
 uBasic: uDZ[4,15] - udz[4,15] >= 0
 uBasic: uDY[5,15] - udy[5,15] >= 0
 uBasic: uDZ[5,15] - udz[5,15] >= 0
 uBasic: uDY[6,15] - udy[6,15] >= 0
 uBasic: uDZ[6,15] - udz[6,15] >= 0
 uBasic: uDY[7,15] - udy[7,15] >= 0
 uBasic: uDZ[7,15] - udz[7,15] >= 0
 uBasic: uDY[8,15] - udy[8,15] >= 0
 uBasic: uDZ[8,15] - udz[8,15] >= 0
 uBasic: uDY[9,15] - udy[9,15] >= 0
 uBasic: uDZ[9,15] - udz[9,15] >= 0
 uBasic: uDY[10,15] - udy[10,15] >= 0
 uBasic: uDZ[10,15] - udz[10,15] >= 0
 uBasic: uDY[11,15] - udy[11,15] >= 0
 uBasic: uDZ[11,15] - udz[11,15] >= 0
 uBasic: uDW[0,15] - udw[0,15] >= 0
 uBasic: uDW[1,15] - udw[1,15] >= 0
 uBasic: uDW[2,15] - udw[2,15] >= 0
 uBasic: uDW[3,15] - udw[3,15] >= 0
 uBasic: uDW[4,15] - udw[4,15] >= 0
 uBasic: uDW[5,15] - udw[5,15] >= 0
 uBasic: uDW[6,15] - udw[6,15] >= 0
 uBasic: uDW[7,15] - udw[7,15] >= 0
 uBasic: uDW[8,15] - udw[8,15] >= 0
 uBasic: uDW[9,15] - udw[9,15] >= 0
 uBasic: uDW[10,15] - udw[10,15] >= 0
 uBasic: uDW[11,15] - udw[11,15] >= 0
 ueSC: uDX[1,0] - uDY[1,0] = 0
 ueSC: udx[1,0] - uDY[1,0] = 0
 ueSC: uDX[1,1] - uDY[1,1] = 0
 ueSC: udx[1,1] - uDY[1,1] = 0
 ueSC: uDX[1,2] - uDY[1,2] = 0
 ueSC: udx[1,2] - uDY[1,2] = 0
 ueSC: uDX[1,3] - uDY[1,3] = 0
 ueSC: udx[1,3] - uDY[1,3] = 0
 ueSC: uDX[1,4] - uDY[1,4] = 0
 ueSC: udx[1,4] - uDY[1,4] = 0
 ueSC: uDX[1,5] - uDY[1,5] = 0
 ueSC: udx[1,5] - uDY[1,5] = 0
 ueSC: uDX[1,6] - uDY[1,6] = 0
 ueSC: udx[1,6] - uDY[1,6] = 0
 ueSC: uDX[1,7] - uDY[1,7] = 0
 ueSC: udx[1,7] - uDY[1,7] = 0
 ueSC: uDX[1,8] - uDY[1,8] = 0
 ueSC: udx[1,8] - uDY[1,8] = 0
 ueSC: uDX[1,9] - uDY[1,9] = 0
 ueSC: udx[1,9] - uDY[1,9] = 0
 ueSC: uDX[1,10] - uDY[1,10] = 0
 ueSC: udx[1,10] - uDY[1,10] = 0
 ueSC: uDX[1,11] - uDY[1,11] = 0
 ueSC: udx[1,11] - uDY[1,11] = 0
 ueSC: uDX[1,12] - uDY[1,12] = 0
 ueSC: udx[1,12] - uDY[1,12] = 0
 ueSC: uDX[1,13] - uDY[1,13] = 0
 ueSC: udx[1,13] - uDY[1,13] = 0
 ueSC: uDX[1,14] - uDY[1,14] = 0
 ueSC: udx[1,14] - uDY[1,14] = 0
 ueSC: uDX[1,15] - uDY[1,15] = 0
 ueSC: udx[1,15] - uDY[1,15] = 0
 ueSC: uDX[2,0] - uDY[2,0] = 0
 ueSC: udx[2,0] - uDY[2,0] = 0
 ueSC: uDX[2,1] - uDY[2,1] = 0
 ueSC: udx[2,1] - uDY[2,1] = 0
 ueSC: uDX[2,2] - uDY[2,2] = 0
 ueSC: udx[2,2] - uDY[2,2] = 0
 ueSC: uDX[2,3] - uDY[2,3] = 0
 ueSC: udx[2,3] - uDY[2,3] = 0
 ueSC: uDX[2,4] - uDY[2,4] = 0
 ueSC: udx[2,4] - uDY[2,4] = 0
 ueSC: uDX[2,5] - uDY[2,5] = 0
 ueSC: udx[2,5] - uDY[2,5] = 0
 ueSC: uDX[2,6] - uDY[2,6] = 0
 ueSC: udx[2,6] - uDY[2,6] = 0
 ueSC: uDX[2,7] - uDY[2,7] = 0
 ueSC: udx[2,7] - uDY[2,7] = 0
 ueSC: uDX[2,8] - uDY[2,8] = 0
 ueSC: udx[2,8] - uDY[2,8] = 0
 ueSC: uDX[2,9] - uDY[2,9] = 0
 ueSC: udx[2,9] - uDY[2,9] = 0
 ueSC: uDX[2,10] - uDY[2,10] = 0
 ueSC: udx[2,10] - uDY[2,10] = 0
 ueSC: uDX[2,11] - uDY[2,11] = 0
 ueSC: udx[2,11] - uDY[2,11] = 0
 ueSC: uDX[2,12] - uDY[2,12] = 0
 ueSC: udx[2,12] - uDY[2,12] = 0
 ueSC: uDX[2,13] - uDY[2,13] = 0
 ueSC: udx[2,13] - uDY[2,13] = 0
 ueSC: uDX[2,14] - uDY[2,14] = 0
 ueSC: udx[2,14] - uDY[2,14] = 0
 ueSC: uDX[2,15] - uDY[2,15] = 0
 ueSC: udx[2,15] - uDY[2,15] = 0
 ueSC: uDX[3,0] - uDY[3,0] = 0
 ueSC: udx[3,0] - uDY[3,0] = 0
 ueSC: uDX[3,1] - uDY[3,1] = 0
 ueSC: udx[3,1] - uDY[3,1] = 0
 ueSC: uDX[3,2] - uDY[3,2] = 0
 ueSC: udx[3,2] - uDY[3,2] = 0
 ueSC: uDX[3,3] - uDY[3,3] = 0
 ueSC: udx[3,3] - uDY[3,3] = 0
 ueSC: uDX[3,4] - uDY[3,4] = 0
 ueSC: udx[3,4] - uDY[3,4] = 0
 ueSC: uDX[3,5] - uDY[3,5] = 0
 ueSC: udx[3,5] - uDY[3,5] = 0
 ueSC: uDX[3,6] - uDY[3,6] = 0
 ueSC: udx[3,6] - uDY[3,6] = 0
 ueSC: uDX[3,7] - uDY[3,7] = 0
 ueSC: udx[3,7] - uDY[3,7] = 0
 ueSC: uDX[3,8] - uDY[3,8] = 0
 ueSC: udx[3,8] - uDY[3,8] = 0
 ueSC: uDX[3,9] - uDY[3,9] = 0
 ueSC: udx[3,9] - uDY[3,9] = 0
 ueSC: uDX[3,10] - uDY[3,10] = 0
 ueSC: udx[3,10] - uDY[3,10] = 0
 ueSC: uDX[3,11] - uDY[3,11] = 0
 ueSC: udx[3,11] - uDY[3,11] = 0
 ueSC: uDX[3,12] - uDY[3,12] = 0
 ueSC: udx[3,12] - uDY[3,12] = 0
 ueSC: uDX[3,13] - uDY[3,13] = 0
 ueSC: udx[3,13] - uDY[3,13] = 0
 ueSC: uDX[3,14] - uDY[3,14] = 0
 ueSC: udx[3,14] - uDY[3,14] = 0
 ueSC: uDX[3,15] - uDY[3,15] = 0
 ueSC: udx[3,15] - uDY[3,15] = 0
 ueSC: uDX[4,0] - uDY[4,0] = 0
 ueSC: udx[4,0] - uDY[4,0] = 0
 ueSC: uDX[4,1] - uDY[4,1] = 0
 ueSC: udx[4,1] - uDY[4,1] = 0
 ueSC: uDX[4,2] - uDY[4,2] = 0
 ueSC: udx[4,2] - uDY[4,2] = 0
 ueSC: uDX[4,3] - uDY[4,3] = 0
 ueSC: udx[4,3] - uDY[4,3] = 0
 ueSC: uDX[4,4] - uDY[4,4] = 0
 ueSC: udx[4,4] - uDY[4,4] = 0
 ueSC: uDX[4,5] - uDY[4,5] = 0
 ueSC: udx[4,5] - uDY[4,5] = 0
 ueSC: uDX[4,6] - uDY[4,6] = 0
 ueSC: udx[4,6] - uDY[4,6] = 0
 ueSC: uDX[4,7] - uDY[4,7] = 0
 ueSC: udx[4,7] - uDY[4,7] = 0
 ueSC: uDX[4,8] - uDY[4,8] = 0
 ueSC: udx[4,8] - uDY[4,8] = 0
 ueSC: uDX[4,9] - uDY[4,9] = 0
 ueSC: udx[4,9] - uDY[4,9] = 0
 ueSC: uDX[4,10] - uDY[4,10] = 0
 ueSC: udx[4,10] - uDY[4,10] = 0
 ueSC: uDX[4,11] - uDY[4,11] = 0
 ueSC: udx[4,11] - uDY[4,11] = 0
 ueSC: uDX[4,12] - uDY[4,12] = 0
 ueSC: udx[4,12] - uDY[4,12] = 0
 ueSC: uDX[4,13] - uDY[4,13] = 0
 ueSC: udx[4,13] - uDY[4,13] = 0
 ueSC: uDX[4,14] - uDY[4,14] = 0
 ueSC: udx[4,14] - uDY[4,14] = 0
 ueSC: uDX[4,15] - uDY[4,15] = 0
 ueSC: udx[4,15] - uDY[4,15] = 0
 udSC: - uDX[5,0] + uDY[5,0] = 0
 udSC: - uDX[5,0] + udy[5,0] = 0
 udSC: - uDX[5,1] + uDY[5,1] = 0
 udSC: - uDX[5,1] + udy[5,1] = 0
 udSC: - uDX[5,2] + uDY[5,2] = 0
 udSC: - uDX[5,2] + udy[5,2] = 0
 udSC: - uDX[5,3] + uDY[5,3] = 0
 udSC: - uDX[5,3] + udy[5,3] = 0
 udSC: - uDX[5,4] + uDY[5,4] = 0
 udSC: - uDX[5,4] + udy[5,4] = 0
 udSC: - uDX[5,5] + uDY[5,5] = 0
 udSC: - uDX[5,5] + udy[5,5] = 0
 udSC: - uDX[5,6] + uDY[5,6] = 0
 udSC: - uDX[5,6] + udy[5,6] = 0
 udSC: - uDX[5,7] + uDY[5,7] = 0
 udSC: - uDX[5,7] + udy[5,7] = 0
 udSC: - uDX[5,8] + uDY[5,8] = 0
 udSC: - uDX[5,8] + udy[5,8] = 0
 udSC: - uDX[5,9] + uDY[5,9] = 0
 udSC: - uDX[5,9] + udy[5,9] = 0
 udSC: - uDX[5,10] + uDY[5,10] = 0
 udSC: - uDX[5,10] + udy[5,10] = 0
 udSC: - uDX[5,11] + uDY[5,11] = 0
 udSC: - uDX[5,11] + udy[5,11] = 0
 udSC: - uDX[5,12] + uDY[5,12] = 0
 udSC: - uDX[5,12] + udy[5,12] = 0
 udSC: - uDX[5,13] + uDY[5,13] = 0
 udSC: - uDX[5,13] + udy[5,13] = 0
 udSC: - uDX[5,14] + uDY[5,14] = 0
 udSC: - uDX[5,14] + udy[5,14] = 0
 udSC: - uDX[5,15] + uDY[5,15] = 0
 udSC: - uDX[5,15] + udy[5,15] = 0
 udSC: - uDX[6,0] + uDY[6,0] = 0
 udSC: - uDX[6,0] + udy[6,0] = 0
 udSC: - uDX[6,1] + uDY[6,1] = 0
 udSC: - uDX[6,1] + udy[6,1] = 0
 udSC: - uDX[6,2] + uDY[6,2] = 0
 udSC: - uDX[6,2] + udy[6,2] = 0
 udSC: - uDX[6,3] + uDY[6,3] = 0
 udSC: - uDX[6,3] + udy[6,3] = 0
 udSC: - uDX[6,4] + uDY[6,4] = 0
 udSC: - uDX[6,4] + udy[6,4] = 0
 udSC: - uDX[6,5] + uDY[6,5] = 0
 udSC: - uDX[6,5] + udy[6,5] = 0
 udSC: - uDX[6,6] + uDY[6,6] = 0
 udSC: - uDX[6,6] + udy[6,6] = 0
 udSC: - uDX[6,7] + uDY[6,7] = 0
 udSC: - uDX[6,7] + udy[6,7] = 0
 udSC: - uDX[6,8] + uDY[6,8] = 0
 udSC: - uDX[6,8] + udy[6,8] = 0
 udSC: - uDX[6,9] + uDY[6,9] = 0
 udSC: - uDX[6,9] + udy[6,9] = 0
 udSC: - uDX[6,10] + uDY[6,10] = 0
 udSC: - uDX[6,10] + udy[6,10] = 0
 udSC: - uDX[6,11] + uDY[6,11] = 0
 udSC: - uDX[6,11] + udy[6,11] = 0
 udSC: - uDX[6,12] + uDY[6,12] = 0
 udSC: - uDX[6,12] + udy[6,12] = 0
 udSC: - uDX[6,13] + uDY[6,13] = 0
 udSC: - uDX[6,13] + udy[6,13] = 0
 udSC: - uDX[6,14] + uDY[6,14] = 0
 udSC: - uDX[6,14] + udy[6,14] = 0
 udSC: - uDX[6,15] + uDY[6,15] = 0
 udSC: - uDX[6,15] + udy[6,15] = 0
 udSC: - uDX[7,0] + uDY[7,0] = 0
 udSC: - uDX[7,0] + udy[7,0] = 0
 udSC: - uDX[7,1] + uDY[7,1] = 0
 udSC: - uDX[7,1] + udy[7,1] = 0
 udSC: - uDX[7,2] + uDY[7,2] = 0
 udSC: - uDX[7,2] + udy[7,2] = 0
 udSC: - uDX[7,3] + uDY[7,3] = 0
 udSC: - uDX[7,3] + udy[7,3] = 0
 udSC: - uDX[7,4] + uDY[7,4] = 0
 udSC: - uDX[7,4] + udy[7,4] = 0
 udSC: - uDX[7,5] + uDY[7,5] = 0
 udSC: - uDX[7,5] + udy[7,5] = 0
 udSC: - uDX[7,6] + uDY[7,6] = 0
 udSC: - uDX[7,6] + udy[7,6] = 0
 udSC: - uDX[7,7] + uDY[7,7] = 0
 udSC: - uDX[7,7] + udy[7,7] = 0
 udSC: - uDX[7,8] + uDY[7,8] = 0
 udSC: - uDX[7,8] + udy[7,8] = 0
 udSC: - uDX[7,9] + uDY[7,9] = 0
 udSC: - uDX[7,9] + udy[7,9] = 0
 udSC: - uDX[7,10] + uDY[7,10] = 0
 udSC: - uDX[7,10] + udy[7,10] = 0
 udSC: - uDX[7,11] + uDY[7,11] = 0
 udSC: - uDX[7,11] + udy[7,11] = 0
 udSC: - uDX[7,12] + uDY[7,12] = 0
 udSC: - uDX[7,12] + udy[7,12] = 0
 udSC: - uDX[7,13] + uDY[7,13] = 0
 udSC: - uDX[7,13] + udy[7,13] = 0
 udSC: - uDX[7,14] + uDY[7,14] = 0
 udSC: - uDX[7,14] + udy[7,14] = 0
 udSC: - uDX[7,15] + uDY[7,15] = 0
 udSC: - uDX[7,15] + udy[7,15] = 0
 udSC: - uDX[8,0] + uDY[8,0] = 0
 udSC: - uDX[8,0] + udy[8,0] = 0
 udSC: - uDX[8,1] + uDY[8,1] = 0
 udSC: - uDX[8,1] + udy[8,1] = 0
 udSC: - uDX[8,2] + uDY[8,2] = 0
 udSC: - uDX[8,2] + udy[8,2] = 0
 udSC: - uDX[8,3] + uDY[8,3] = 0
 udSC: - uDX[8,3] + udy[8,3] = 0
 udSC: - uDX[8,4] + uDY[8,4] = 0
 udSC: - uDX[8,4] + udy[8,4] = 0
 udSC: - uDX[8,5] + uDY[8,5] = 0
 udSC: - uDX[8,5] + udy[8,5] = 0
 udSC: - uDX[8,6] + uDY[8,6] = 0
 udSC: - uDX[8,6] + udy[8,6] = 0
 udSC: - uDX[8,7] + uDY[8,7] = 0
 udSC: - uDX[8,7] + udy[8,7] = 0
 udSC: - uDX[8,8] + uDY[8,8] = 0
 udSC: - uDX[8,8] + udy[8,8] = 0
 udSC: - uDX[8,9] + uDY[8,9] = 0
 udSC: - uDX[8,9] + udy[8,9] = 0
 udSC: - uDX[8,10] + uDY[8,10] = 0
 udSC: - uDX[8,10] + udy[8,10] = 0
 udSC: - uDX[8,11] + uDY[8,11] = 0
 udSC: - uDX[8,11] + udy[8,11] = 0
 udSC: - uDX[8,12] + uDY[8,12] = 0
 udSC: - uDX[8,12] + udy[8,12] = 0
 udSC: - uDX[8,13] + uDY[8,13] = 0
 udSC: - uDX[8,13] + udy[8,13] = 0
 udSC: - uDX[8,14] + uDY[8,14] = 0
 udSC: - uDX[8,14] + udy[8,14] = 0
 udSC: - uDX[8,15] + uDY[8,15] = 0
 udSC: - uDX[8,15] + udy[8,15] = 0
 udSC: - uDX[9,0] + uDY[9,0] = 0
 udSC: - uDX[9,0] + udy[9,0] = 0
 udSC: - uDX[9,1] + uDY[9,1] = 0
 udSC: - uDX[9,1] + udy[9,1] = 0
 udSC: - uDX[9,2] + uDY[9,2] = 0
 udSC: - uDX[9,2] + udy[9,2] = 0
 udSC: - uDX[9,3] + uDY[9,3] = 0
 udSC: - uDX[9,3] + udy[9,3] = 0
 udSC: - uDX[9,4] + uDY[9,4] = 0
 udSC: - uDX[9,4] + udy[9,4] = 0
 udSC: - uDX[9,5] + uDY[9,5] = 0
 udSC: - uDX[9,5] + udy[9,5] = 0
 udSC: - uDX[9,6] + uDY[9,6] = 0
 udSC: - uDX[9,6] + udy[9,6] = 0
 udSC: - uDX[9,7] + uDY[9,7] = 0
 udSC: - uDX[9,7] + udy[9,7] = 0
 udSC: - uDX[9,8] + uDY[9,8] = 0
 udSC: - uDX[9,8] + udy[9,8] = 0
 udSC: - uDX[9,9] + uDY[9,9] = 0
 udSC: - uDX[9,9] + udy[9,9] = 0
 udSC: - uDX[9,10] + uDY[9,10] = 0
 udSC: - uDX[9,10] + udy[9,10] = 0
 udSC: - uDX[9,11] + uDY[9,11] = 0
 udSC: - uDX[9,11] + udy[9,11] = 0
 udSC: - uDX[9,12] + uDY[9,12] = 0
 udSC: - uDX[9,12] + udy[9,12] = 0
 udSC: - uDX[9,13] + uDY[9,13] = 0
 udSC: - uDX[9,13] + udy[9,13] = 0
 udSC: - uDX[9,14] + uDY[9,14] = 0
 udSC: - uDX[9,14] + udy[9,14] = 0
 udSC: - uDX[9,15] + uDY[9,15] = 0
 udSC: - uDX[9,15] + udy[9,15] = 0
 udSC: - uDX[10,0] + uDY[10,0] = 0
 udSC: - uDX[10,0] + udy[10,0] = 0
 udSC: - uDX[10,1] + uDY[10,1] = 0
 udSC: - uDX[10,1] + udy[10,1] = 0
 udSC: - uDX[10,2] + uDY[10,2] = 0
 udSC: - uDX[10,2] + udy[10,2] = 0
 udSC: - uDX[10,3] + uDY[10,3] = 0
 udSC: - uDX[10,3] + udy[10,3] = 0
 udSC: - uDX[10,4] + uDY[10,4] = 0
 udSC: - uDX[10,4] + udy[10,4] = 0
 udSC: - uDX[10,5] + uDY[10,5] = 0
 udSC: - uDX[10,5] + udy[10,5] = 0
 udSC: - uDX[10,6] + uDY[10,6] = 0
 udSC: - uDX[10,6] + udy[10,6] = 0
 udSC: - uDX[10,7] + uDY[10,7] = 0
 udSC: - uDX[10,7] + udy[10,7] = 0
 udSC: - uDX[10,8] + uDY[10,8] = 0
 udSC: - uDX[10,8] + udy[10,8] = 0
 udSC: - uDX[10,9] + uDY[10,9] = 0
 udSC: - uDX[10,9] + udy[10,9] = 0
 udSC: - uDX[10,10] + uDY[10,10] = 0
 udSC: - uDX[10,10] + udy[10,10] = 0
 udSC: - uDX[10,11] + uDY[10,11] = 0
 udSC: - uDX[10,11] + udy[10,11] = 0
 udSC: - uDX[10,12] + uDY[10,12] = 0
 udSC: - uDX[10,12] + udy[10,12] = 0
 udSC: - uDX[10,13] + uDY[10,13] = 0
 udSC: - uDX[10,13] + udy[10,13] = 0
 udSC: - uDX[10,14] + uDY[10,14] = 0
 udSC: - uDX[10,14] + udy[10,14] = 0
 udSC: - uDX[10,15] + uDY[10,15] = 0
 udSC: - uDX[10,15] + udy[10,15] = 0
 udSC: - uDX[11,0] + uDY[11,0] = 0
 udSC: - uDX[11,0] + udy[11,0] = 0
 udSC: - uDX[11,1] + uDY[11,1] = 0
 udSC: - uDX[11,1] + udy[11,1] = 0
 udSC: - uDX[11,2] + uDY[11,2] = 0
 udSC: - uDX[11,2] + udy[11,2] = 0
 udSC: - uDX[11,3] + uDY[11,3] = 0
 udSC: - uDX[11,3] + udy[11,3] = 0
 udSC: - uDX[11,4] + uDY[11,4] = 0
 udSC: - uDX[11,4] + udy[11,4] = 0
 udSC: - uDX[11,5] + uDY[11,5] = 0
 udSC: - uDX[11,5] + udy[11,5] = 0
 udSC: - uDX[11,6] + uDY[11,6] = 0
 udSC: - uDX[11,6] + udy[11,6] = 0
 udSC: - uDX[11,7] + uDY[11,7] = 0
 udSC: - uDX[11,7] + udy[11,7] = 0
 udSC: - uDX[11,8] + uDY[11,8] = 0
 udSC: - uDX[11,8] + udy[11,8] = 0
 udSC: - uDX[11,9] + uDY[11,9] = 0
 udSC: - uDX[11,9] + udy[11,9] = 0
 udSC: - uDX[11,10] + uDY[11,10] = 0
 udSC: - uDX[11,10] + udy[11,10] = 0
 udSC: - uDX[11,11] + uDY[11,11] = 0
 udSC: - uDX[11,11] + udy[11,11] = 0
 udSC: - uDX[11,12] + uDY[11,12] = 0
 udSC: - uDX[11,12] + udy[11,12] = 0
 udSC: - uDX[11,13] + uDY[11,13] = 0
 udSC: - uDX[11,13] + udy[11,13] = 0
 udSC: - uDX[11,14] + uDY[11,14] = 0
 udSC: - uDX[11,14] + udy[11,14] = 0
 udSC: - uDX[11,15] + uDY[11,15] = 0
 udSC: - uDX[11,15] + udy[11,15] = 0
 eqk0: - ustk[0,13] + eqk[0] >= 0
 eqk0: eqk[0] - eqk[12] >= 0
 eqk0: - ustk[0,13] + eqk[0] - eqk[12] <= 0
 eqk1: - ustk[0,0] + eqk[4] = 0
 eqk2: - ustk[0,7] + eqk[8] >= 0
 eqk2: - ustk[0,10] + eqk[8] >= 0
 eqk2: - ustk[0,7] - ustk[0,10] + eqk[8] <= 0
 eqk2: - ustk[0,0] + eqk[12] >= 0
 eqk2: - ustk[0,10] + eqk[12] >= 0
 eqk2: - ustk[0,0] - ustk[0,10] + eqk[12] <= 0
 eqk0: - ustk[0,14] + eqk[1] >= 0
 eqk0: eqk[1] - eqk[13] >= 0
 eqk0: - ustk[0,14] + eqk[1] - eqk[13] <= 0
 eqk1: - ustk[0,1] + eqk[5] = 0
 eqk2: - ustk[0,4] + eqk[9] >= 0
 eqk2: - ustk[0,11] + eqk[9] >= 0
 eqk2: - ustk[0,4] - ustk[0,11] + eqk[9] <= 0
 eqk2: - ustk[0,1] + eqk[13] >= 0
 eqk2: - ustk[0,11] + eqk[13] >= 0
 eqk2: - ustk[0,1] - ustk[0,11] + eqk[13] <= 0
 eqk0: - ustk[0,15] + eqk[2] >= 0
 eqk0: eqk[2] - eqk[14] >= 0
 eqk0: - ustk[0,15] + eqk[2] - eqk[14] <= 0
 eqk1: - ustk[0,2] + eqk[6] = 0
 eqk2: - ustk[0,5] + eqk[10] >= 0
 eqk2: - ustk[0,8] + eqk[10] >= 0
 eqk2: - ustk[0,5] - ustk[0,8] + eqk[10] <= 0
 eqk2: - ustk[0,2] + eqk[14] >= 0
 eqk2: - ustk[0,8] + eqk[14] >= 0
 eqk2: - ustk[0,2] - ustk[0,8] + eqk[14] <= 0
 eqk0: - ustk[0,12] + eqk[3] >= 0
 eqk0: eqk[3] - eqk[15] >= 0
 eqk0: - ustk[0,12] + eqk[3] - eqk[15] <= 0
 eqk1: - ustk[0,3] + eqk[7] = 0
 eqk2: - ustk[0,6] + eqk[11] >= 0
 eqk2: - ustk[0,9] + eqk[11] >= 0
 eqk2: - ustk[0,6] - ustk[0,9] + eqk[11] <= 0
 eqk2: - ustk[0,3] + eqk[15] >= 0
 eqk2: - ustk[0,9] + eqk[15] >= 0
 eqk2: - ustk[0,3] - ustk[0,9] + eqk[15] <= 0
 uAEqk: - udx[1,0] + udw[0,0] = 0
 uAEqk: - uDX[1,0] + uDW[0,0] + eqk[0] >= 0
 uAEqk: uDX[1,0] - uDW[0,0] + eqk[0] >= 0
 uAEqk: uDX[1,0] - udx[1,0] + uDW[0,0] - eqk[0] >= 0
 uAEqk: - udx[1,0] + uDW[0,0] >= 0
 uAEqk: - udx[1,1] + udw[0,1] = 0
 uAEqk: - uDX[1,1] + uDW[0,1] + eqk[1] >= 0
 uAEqk: uDX[1,1] - uDW[0,1] + eqk[1] >= 0
 uAEqk: uDX[1,1] - udx[1,1] + uDW[0,1] - eqk[1] >= 0
 uAEqk: - udx[1,1] + uDW[0,1] >= 0
 uAEqk: - udx[1,2] + udw[0,2] = 0
 uAEqk: - uDX[1,2] + uDW[0,2] + eqk[2] >= 0
 uAEqk: uDX[1,2] - uDW[0,2] + eqk[2] >= 0
 uAEqk: uDX[1,2] - udx[1,2] + uDW[0,2] - eqk[2] >= 0
 uAEqk: - udx[1,2] + uDW[0,2] >= 0
 uAEqk: - udx[1,3] + udw[0,3] = 0
 uAEqk: - uDX[1,3] + uDW[0,3] + eqk[3] >= 0
 uAEqk: uDX[1,3] - uDW[0,3] + eqk[3] >= 0
 uAEqk: uDX[1,3] - udx[1,3] + uDW[0,3] - eqk[3] >= 0
 uAEqk: - udx[1,3] + uDW[0,3] >= 0
 uAEqk: - udx[1,4] + udw[0,4] = 0
 uAEqk: - uDX[1,4] + uDW[0,4] + eqk[4] >= 0
 uAEqk: uDX[1,4] - uDW[0,4] + eqk[4] >= 0
 uAEqk: uDX[1,4] - udx[1,4] + uDW[0,4] - eqk[4] >= 0
 uAEqk: - udx[1,4] + uDW[0,4] >= 0
 uAEqk: - udx[1,5] + udw[0,5] = 0
 uAEqk: - uDX[1,5] + uDW[0,5] + eqk[5] >= 0
 uAEqk: uDX[1,5] - uDW[0,5] + eqk[5] >= 0
 uAEqk: uDX[1,5] - udx[1,5] + uDW[0,5] - eqk[5] >= 0
 uAEqk: - udx[1,5] + uDW[0,5] >= 0
 uAEqk: - udx[1,6] + udw[0,6] = 0
 uAEqk: - uDX[1,6] + uDW[0,6] + eqk[6] >= 0
 uAEqk: uDX[1,6] - uDW[0,6] + eqk[6] >= 0
 uAEqk: uDX[1,6] - udx[1,6] + uDW[0,6] - eqk[6] >= 0
 uAEqk: - udx[1,6] + uDW[0,6] >= 0
 uAEqk: - udx[1,7] + udw[0,7] = 0
 uAEqk: - uDX[1,7] + uDW[0,7] + eqk[7] >= 0
 uAEqk: uDX[1,7] - uDW[0,7] + eqk[7] >= 0
 uAEqk: uDX[1,7] - udx[1,7] + uDW[0,7] - eqk[7] >= 0
 uAEqk: - udx[1,7] + uDW[0,7] >= 0
 uAEqk: - udx[1,8] + udw[0,8] = 0
 uAEqk: - uDX[1,8] + uDW[0,8] + eqk[8] >= 0
 uAEqk: uDX[1,8] - uDW[0,8] + eqk[8] >= 0
 uAEqk: uDX[1,8] - udx[1,8] + uDW[0,8] - eqk[8] >= 0
 uAEqk: - udx[1,8] + uDW[0,8] >= 0
 uAEqk: - udx[1,9] + udw[0,9] = 0
 uAEqk: - uDX[1,9] + uDW[0,9] + eqk[9] >= 0
 uAEqk: uDX[1,9] - uDW[0,9] + eqk[9] >= 0
 uAEqk: uDX[1,9] - udx[1,9] + uDW[0,9] - eqk[9] >= 0
 uAEqk: - udx[1,9] + uDW[0,9] >= 0
 uAEqk: - udx[1,10] + udw[0,10] = 0
 uAEqk: - uDX[1,10] + uDW[0,10] + eqk[10] >= 0
 uAEqk: uDX[1,10] - uDW[0,10] + eqk[10] >= 0
 uAEqk: uDX[1,10] - udx[1,10] + uDW[0,10] - eqk[10] >= 0
 uAEqk: - udx[1,10] + uDW[0,10] >= 0
 uAEqk: - udx[1,11] + udw[0,11] = 0
 uAEqk: - uDX[1,11] + uDW[0,11] + eqk[11] >= 0
 uAEqk: uDX[1,11] - uDW[0,11] + eqk[11] >= 0
 uAEqk: uDX[1,11] - udx[1,11] + uDW[0,11] - eqk[11] >= 0
 uAEqk: - udx[1,11] + uDW[0,11] >= 0
 uAEqk: - udx[1,12] + udw[0,12] = 0
 uAEqk: - uDX[1,12] + uDW[0,12] + eqk[12] >= 0
 uAEqk: uDX[1,12] - uDW[0,12] + eqk[12] >= 0
 uAEqk: uDX[1,12] - udx[1,12] + uDW[0,12] - eqk[12] >= 0
 uAEqk: - udx[1,12] + uDW[0,12] >= 0
 uAEqk: - udx[1,13] + udw[0,13] = 0
 uAEqk: - uDX[1,13] + uDW[0,13] + eqk[13] >= 0
 uAEqk: uDX[1,13] - uDW[0,13] + eqk[13] >= 0
 uAEqk: uDX[1,13] - udx[1,13] + uDW[0,13] - eqk[13] >= 0
 uAEqk: - udx[1,13] + uDW[0,13] >= 0
 uAEqk: - udx[1,14] + udw[0,14] = 0
 uAEqk: - uDX[1,14] + uDW[0,14] + eqk[14] >= 0
 uAEqk: uDX[1,14] - uDW[0,14] + eqk[14] >= 0
 uAEqk: uDX[1,14] - udx[1,14] + uDW[0,14] - eqk[14] >= 0
 uAEqk: - udx[1,14] + uDW[0,14] >= 0
 uAEqk: - udx[1,15] + udw[0,15] = 0
 uAEqk: - uDX[1,15] + uDW[0,15] + eqk[15] >= 0
 uAEqk: uDX[1,15] - uDW[0,15] + eqk[15] >= 0
 uAEqk: uDX[1,15] - udx[1,15] + uDW[0,15] - eqk[15] >= 0
 uAEqk: - udx[1,15] + uDW[0,15] >= 0
 uATKri: udy[1,0] - udz[1,0] = 0
 uATKri: uDY[1,0] - uDZ[1,0] + ustk[1,0] >= 0
 uATKri: - uDY[1,0] + uDZ[1,0] + ustk[1,0] >= 0
 uATKri: uDY[1,0] + uDZ[1,0] - udz[1,0] - ustk[1,0] >= 0
 uATKri: uDY[1,0] - udz[1,0] >= 0
 uATKri: udy[1,1] - udz[1,1] = 0
 uATKri: uDY[1,1] - uDZ[1,1] + ustk[1,1] >= 0
 uATKri: - uDY[1,1] + uDZ[1,1] + ustk[1,1] >= 0
 uATKri: uDY[1,1] + uDZ[1,1] - udz[1,1] - ustk[1,1] >= 0
 uATKri: uDY[1,1] - udz[1,1] >= 0
 uATKri: udy[1,2] - udz[1,2] = 0
 uATKri: uDY[1,2] - uDZ[1,2] + ustk[1,2] >= 0
 uATKri: - uDY[1,2] + uDZ[1,2] + ustk[1,2] >= 0
 uATKri: uDY[1,2] + uDZ[1,2] - udz[1,2] - ustk[1,2] >= 0
 uATKri: uDY[1,2] - udz[1,2] >= 0
 uATKri: udy[1,3] - udz[1,3] = 0
 uATKri: uDY[1,3] - uDZ[1,3] + ustk[1,3] >= 0
 uATKri: - uDY[1,3] + uDZ[1,3] + ustk[1,3] >= 0
 uATKri: uDY[1,3] + uDZ[1,3] - udz[1,3] - ustk[1,3] >= 0
 uATKri: uDY[1,3] - udz[1,3] >= 0
 uATKri: udy[1,4] - udz[1,4] = 0
 uATKri: uDY[1,4] - uDZ[1,4] + ustk[1,4] >= 0
 uATKri: - uDY[1,4] + uDZ[1,4] + ustk[1,4] >= 0
 uATKri: uDY[1,4] + uDZ[1,4] - udz[1,4] - ustk[1,4] >= 0
 uATKri: uDY[1,4] - udz[1,4] >= 0
 uATKri: udy[1,5] - udz[1,5] = 0
 uATKri: uDY[1,5] - uDZ[1,5] + ustk[1,5] >= 0
 uATKri: - uDY[1,5] + uDZ[1,5] + ustk[1,5] >= 0
 uATKri: uDY[1,5] + uDZ[1,5] - udz[1,5] - ustk[1,5] >= 0
 uATKri: uDY[1,5] - udz[1,5] >= 0
 uATKri: udy[1,6] - udz[1,6] = 0
 uATKri: uDY[1,6] - uDZ[1,6] + ustk[1,6] >= 0
 uATKri: - uDY[1,6] + uDZ[1,6] + ustk[1,6] >= 0
 uATKri: uDY[1,6] + uDZ[1,6] - udz[1,6] - ustk[1,6] >= 0
 uATKri: uDY[1,6] - udz[1,6] >= 0
 uATKri: udy[1,7] - udz[1,7] = 0
 uATKri: uDY[1,7] - uDZ[1,7] + ustk[1,7] >= 0
 uATKri: - uDY[1,7] + uDZ[1,7] + ustk[1,7] >= 0
 uATKri: uDY[1,7] + uDZ[1,7] - udz[1,7] - ustk[1,7] >= 0
 uATKri: uDY[1,7] - udz[1,7] >= 0
 uATKri: uDY[1,8] - uDZ[1,8] = 0
 uATKri: udy[1,8] - udz[1,8] = 0
 uATKri: uDY[1,9] - uDZ[1,9] = 0
 uATKri: udy[1,9] - udz[1,9] = 0
 uATKri: uDY[1,10] - uDZ[1,10] = 0
 uATKri: udy[1,10] - udz[1,10] = 0
 uATKri: uDY[1,11] - uDZ[1,11] = 0
 uATKri: udy[1,11] - udz[1,11] = 0
 uATKri: uDY[1,12] - uDZ[1,12] = 0
 uATKri: udy[1,12] - udz[1,12] = 0
 uATKri: uDY[1,13] - uDZ[1,13] = 0
 uATKri: udy[1,13] - udz[1,13] = 0
 uATKri: uDY[1,14] - uDZ[1,14] = 0
 uATKri: udy[1,14] - udz[1,14] = 0
 uATKri: uDY[1,15] - uDZ[1,15] = 0
 uATKri: udy[1,15] - udz[1,15] = 0
 uATKri: udy[2,0] - udz[2,0] = 0
 uATKri: uDY[2,0] - uDZ[2,0] + ustk[2,0] >= 0
 uATKri: - uDY[2,0] + uDZ[2,0] + ustk[2,0] >= 0
 uATKri: uDY[2,0] + uDZ[2,0] - udz[2,0] - ustk[2,0] >= 0
 uATKri: uDY[2,0] - udz[2,0] >= 0
 uATKri: udy[2,1] - udz[2,1] = 0
 uATKri: uDY[2,1] - uDZ[2,1] + ustk[2,1] >= 0
 uATKri: - uDY[2,1] + uDZ[2,1] + ustk[2,1] >= 0
 uATKri: uDY[2,1] + uDZ[2,1] - udz[2,1] - ustk[2,1] >= 0
 uATKri: uDY[2,1] - udz[2,1] >= 0
 uATKri: udy[2,2] - udz[2,2] = 0
 uATKri: uDY[2,2] - uDZ[2,2] + ustk[2,2] >= 0
 uATKri: - uDY[2,2] + uDZ[2,2] + ustk[2,2] >= 0
 uATKri: uDY[2,2] + uDZ[2,2] - udz[2,2] - ustk[2,2] >= 0
 uATKri: uDY[2,2] - udz[2,2] >= 0
 uATKri: udy[2,3] - udz[2,3] = 0
 uATKri: uDY[2,3] - uDZ[2,3] + ustk[2,3] >= 0
 uATKri: - uDY[2,3] + uDZ[2,3] + ustk[2,3] >= 0
 uATKri: uDY[2,3] + uDZ[2,3] - udz[2,3] - ustk[2,3] >= 0
 uATKri: uDY[2,3] - udz[2,3] >= 0
 uATKri: udy[2,4] - udz[2,4] = 0
 uATKri: uDY[2,4] - uDZ[2,4] + ustk[2,4] >= 0
 uATKri: - uDY[2,4] + uDZ[2,4] + ustk[2,4] >= 0
 uATKri: uDY[2,4] + uDZ[2,4] - udz[2,4] - ustk[2,4] >= 0
 uATKri: uDY[2,4] - udz[2,4] >= 0
 uATKri: udy[2,5] - udz[2,5] = 0
 uATKri: uDY[2,5] - uDZ[2,5] + ustk[2,5] >= 0
 uATKri: - uDY[2,5] + uDZ[2,5] + ustk[2,5] >= 0
 uATKri: uDY[2,5] + uDZ[2,5] - udz[2,5] - ustk[2,5] >= 0
 uATKri: uDY[2,5] - udz[2,5] >= 0
 uATKri: udy[2,6] - udz[2,6] = 0
 uATKri: uDY[2,6] - uDZ[2,6] + ustk[2,6] >= 0
 uATKri: - uDY[2,6] + uDZ[2,6] + ustk[2,6] >= 0
 uATKri: uDY[2,6] + uDZ[2,6] - udz[2,6] - ustk[2,6] >= 0
 uATKri: uDY[2,6] - udz[2,6] >= 0
 uATKri: udy[2,7] - udz[2,7] = 0
 uATKri: uDY[2,7] - uDZ[2,7] + ustk[2,7] >= 0
 uATKri: - uDY[2,7] + uDZ[2,7] + ustk[2,7] >= 0
 uATKri: uDY[2,7] + uDZ[2,7] - udz[2,7] - ustk[2,7] >= 0
 uATKri: uDY[2,7] - udz[2,7] >= 0
 uATKri: uDY[2,8] - uDZ[2,8] = 0
 uATKri: udy[2,8] - udz[2,8] = 0
 uATKri: uDY[2,9] - uDZ[2,9] = 0
 uATKri: udy[2,9] - udz[2,9] = 0
 uATKri: uDY[2,10] - uDZ[2,10] = 0
 uATKri: udy[2,10] - udz[2,10] = 0
 uATKri: uDY[2,11] - uDZ[2,11] = 0
 uATKri: udy[2,11] - udz[2,11] = 0
 uATKri: uDY[2,12] - uDZ[2,12] = 0
 uATKri: udy[2,12] - udz[2,12] = 0
 uATKri: uDY[2,13] - uDZ[2,13] = 0
 uATKri: udy[2,13] - udz[2,13] = 0
 uATKri: uDY[2,14] - uDZ[2,14] = 0
 uATKri: udy[2,14] - udz[2,14] = 0
 uATKri: uDY[2,15] - uDZ[2,15] = 0
 uATKri: udy[2,15] - udz[2,15] = 0
 uATKri: udy[3,0] - udz[3,0] = 0
 uATKri: uDY[3,0] - uDZ[3,0] + ustk[3,0] >= 0
 uATKri: - uDY[3,0] + uDZ[3,0] + ustk[3,0] >= 0
 uATKri: uDY[3,0] + uDZ[3,0] - udz[3,0] - ustk[3,0] >= 0
 uATKri: uDY[3,0] - udz[3,0] >= 0
 uATKri: udy[3,1] - udz[3,1] = 0
 uATKri: uDY[3,1] - uDZ[3,1] + ustk[3,1] >= 0
 uATKri: - uDY[3,1] + uDZ[3,1] + ustk[3,1] >= 0
 uATKri: uDY[3,1] + uDZ[3,1] - udz[3,1] - ustk[3,1] >= 0
 uATKri: uDY[3,1] - udz[3,1] >= 0
 uATKri: udy[3,2] - udz[3,2] = 0
 uATKri: uDY[3,2] - uDZ[3,2] + ustk[3,2] >= 0
 uATKri: - uDY[3,2] + uDZ[3,2] + ustk[3,2] >= 0
 uATKri: uDY[3,2] + uDZ[3,2] - udz[3,2] - ustk[3,2] >= 0
 uATKri: uDY[3,2] - udz[3,2] >= 0
 uATKri: udy[3,3] - udz[3,3] = 0
 uATKri: uDY[3,3] - uDZ[3,3] + ustk[3,3] >= 0
 uATKri: - uDY[3,3] + uDZ[3,3] + ustk[3,3] >= 0
 uATKri: uDY[3,3] + uDZ[3,3] - udz[3,3] - ustk[3,3] >= 0
 uATKri: uDY[3,3] - udz[3,3] >= 0
 uATKri: udy[3,4] - udz[3,4] = 0
 uATKri: uDY[3,4] - uDZ[3,4] + ustk[3,4] >= 0
 uATKri: - uDY[3,4] + uDZ[3,4] + ustk[3,4] >= 0
 uATKri: uDY[3,4] + uDZ[3,4] - udz[3,4] - ustk[3,4] >= 0
 uATKri: uDY[3,4] - udz[3,4] >= 0
 uATKri: udy[3,5] - udz[3,5] = 0
 uATKri: uDY[3,5] - uDZ[3,5] + ustk[3,5] >= 0
 uATKri: - uDY[3,5] + uDZ[3,5] + ustk[3,5] >= 0
 uATKri: uDY[3,5] + uDZ[3,5] - udz[3,5] - ustk[3,5] >= 0
 uATKri: uDY[3,5] - udz[3,5] >= 0
 uATKri: udy[3,6] - udz[3,6] = 0
 uATKri: uDY[3,6] - uDZ[3,6] + ustk[3,6] >= 0
 uATKri: - uDY[3,6] + uDZ[3,6] + ustk[3,6] >= 0
 uATKri: uDY[3,6] + uDZ[3,6] - udz[3,6] - ustk[3,6] >= 0
 uATKri: uDY[3,6] - udz[3,6] >= 0
 uATKri: udy[3,7] - udz[3,7] = 0
 uATKri: uDY[3,7] - uDZ[3,7] + ustk[3,7] >= 0
 uATKri: - uDY[3,7] + uDZ[3,7] + ustk[3,7] >= 0
 uATKri: uDY[3,7] + uDZ[3,7] - udz[3,7] - ustk[3,7] >= 0
 uATKri: uDY[3,7] - udz[3,7] >= 0
 uATKri: uDY[3,8] - uDZ[3,8] = 0
 uATKri: udy[3,8] - udz[3,8] = 0
 uATKri: uDY[3,9] - uDZ[3,9] = 0
 uATKri: udy[3,9] - udz[3,9] = 0
 uATKri: uDY[3,10] - uDZ[3,10] = 0
 uATKri: udy[3,10] - udz[3,10] = 0
 uATKri: uDY[3,11] - uDZ[3,11] = 0
 uATKri: udy[3,11] - udz[3,11] = 0
 uATKri: uDY[3,12] - uDZ[3,12] = 0
 uATKri: udy[3,12] - udz[3,12] = 0
 uATKri: uDY[3,13] - uDZ[3,13] = 0
 uATKri: udy[3,13] - udz[3,13] = 0
 uATKri: uDY[3,14] - uDZ[3,14] = 0
 uATKri: udy[3,14] - udz[3,14] = 0
 uATKri: uDY[3,15] - uDZ[3,15] = 0
 uATKri: udy[3,15] - udz[3,15] = 0
 uATKri: udy[4,0] - udz[4,0] = 0
 uATKri: uDY[4,0] - uDZ[4,0] + ustk[4,0] >= 0
 uATKri: - uDY[4,0] + uDZ[4,0] + ustk[4,0] >= 0
 uATKri: uDY[4,0] + uDZ[4,0] - udz[4,0] - ustk[4,0] >= 0
 uATKri: uDY[4,0] - udz[4,0] >= 0
 uATKri: udy[4,1] - udz[4,1] = 0
 uATKri: uDY[4,1] - uDZ[4,1] + ustk[4,1] >= 0
 uATKri: - uDY[4,1] + uDZ[4,1] + ustk[4,1] >= 0
 uATKri: uDY[4,1] + uDZ[4,1] - udz[4,1] - ustk[4,1] >= 0
 uATKri: uDY[4,1] - udz[4,1] >= 0
 uATKri: udy[4,2] - udz[4,2] = 0
 uATKri: uDY[4,2] - uDZ[4,2] + ustk[4,2] >= 0
 uATKri: - uDY[4,2] + uDZ[4,2] + ustk[4,2] >= 0
 uATKri: uDY[4,2] + uDZ[4,2] - udz[4,2] - ustk[4,2] >= 0
 uATKri: uDY[4,2] - udz[4,2] >= 0
 uATKri: udy[4,3] - udz[4,3] = 0
 uATKri: uDY[4,3] - uDZ[4,3] + ustk[4,3] >= 0
 uATKri: - uDY[4,3] + uDZ[4,3] + ustk[4,3] >= 0
 uATKri: uDY[4,3] + uDZ[4,3] - udz[4,3] - ustk[4,3] >= 0
 uATKri: uDY[4,3] - udz[4,3] >= 0
 uATKri: udy[4,4] - udz[4,4] = 0
 uATKri: uDY[4,4] - uDZ[4,4] + ustk[4,4] >= 0
 uATKri: - uDY[4,4] + uDZ[4,4] + ustk[4,4] >= 0
 uATKri: uDY[4,4] + uDZ[4,4] - udz[4,4] - ustk[4,4] >= 0
 uATKri: uDY[4,4] - udz[4,4] >= 0
 uATKri: udy[4,5] - udz[4,5] = 0
 uATKri: uDY[4,5] - uDZ[4,5] + ustk[4,5] >= 0
 uATKri: - uDY[4,5] + uDZ[4,5] + ustk[4,5] >= 0
 uATKri: uDY[4,5] + uDZ[4,5] - udz[4,5] - ustk[4,5] >= 0
 uATKri: uDY[4,5] - udz[4,5] >= 0
 uATKri: udy[4,6] - udz[4,6] = 0
 uATKri: uDY[4,6] - uDZ[4,6] + ustk[4,6] >= 0
 uATKri: - uDY[4,6] + uDZ[4,6] + ustk[4,6] >= 0
 uATKri: uDY[4,6] + uDZ[4,6] - udz[4,6] - ustk[4,6] >= 0
 uATKri: uDY[4,6] - udz[4,6] >= 0
 uATKri: udy[4,7] - udz[4,7] = 0
 uATKri: uDY[4,7] - uDZ[4,7] + ustk[4,7] >= 0
 uATKri: - uDY[4,7] + uDZ[4,7] + ustk[4,7] >= 0
 uATKri: uDY[4,7] + uDZ[4,7] - udz[4,7] - ustk[4,7] >= 0
 uATKri: uDY[4,7] - udz[4,7] >= 0
 uATKri: uDY[4,8] - uDZ[4,8] = 0
 uATKri: udy[4,8] - udz[4,8] = 0
 uATKri: uDY[4,9] - uDZ[4,9] = 0
 uATKri: udy[4,9] - udz[4,9] = 0
 uATKri: uDY[4,10] - uDZ[4,10] = 0
 uATKri: udy[4,10] - udz[4,10] = 0
 uATKri: uDY[4,11] - uDZ[4,11] = 0
 uATKri: udy[4,11] - udz[4,11] = 0
 uATKri: uDY[4,12] - uDZ[4,12] = 0
 uATKri: udy[4,12] - udz[4,12] = 0
 uATKri: uDY[4,13] - uDZ[4,13] = 0
 uATKri: udy[4,13] - udz[4,13] = 0
 uATKri: uDY[4,14] - uDZ[4,14] = 0
 uATKri: udy[4,14] - udz[4,14] = 0
 uATKri: uDY[4,15] - uDZ[4,15] = 0
 uATKri: udy[4,15] - udz[4,15] = 0
 uATKri: udy[5,0] - udz[5,0] = 0
 uATKri: uDY[5,0] - uDZ[5,0] + ustk[5,0] >= 0
 uATKri: - uDY[5,0] + uDZ[5,0] + ustk[5,0] >= 0
 uATKri: uDY[5,0] + uDZ[5,0] - udz[5,0] - ustk[5,0] >= 0
 uATKri: uDY[5,0] - udz[5,0] >= 0
 uATKri: udy[5,1] - udz[5,1] = 0
 uATKri: uDY[5,1] - uDZ[5,1] + ustk[5,1] >= 0
 uATKri: - uDY[5,1] + uDZ[5,1] + ustk[5,1] >= 0
 uATKri: uDY[5,1] + uDZ[5,1] - udz[5,1] - ustk[5,1] >= 0
 uATKri: uDY[5,1] - udz[5,1] >= 0
 uATKri: udy[5,2] - udz[5,2] = 0
 uATKri: uDY[5,2] - uDZ[5,2] + ustk[5,2] >= 0
 uATKri: - uDY[5,2] + uDZ[5,2] + ustk[5,2] >= 0
 uATKri: uDY[5,2] + uDZ[5,2] - udz[5,2] - ustk[5,2] >= 0
 uATKri: uDY[5,2] - udz[5,2] >= 0
 uATKri: udy[5,3] - udz[5,3] = 0
 uATKri: uDY[5,3] - uDZ[5,3] + ustk[5,3] >= 0
 uATKri: - uDY[5,3] + uDZ[5,3] + ustk[5,3] >= 0
 uATKri: uDY[5,3] + uDZ[5,3] - udz[5,3] - ustk[5,3] >= 0
 uATKri: uDY[5,3] - udz[5,3] >= 0
 uATKri: udy[5,4] - udz[5,4] = 0
 uATKri: uDY[5,4] - uDZ[5,4] + ustk[5,4] >= 0
 uATKri: - uDY[5,4] + uDZ[5,4] + ustk[5,4] >= 0
 uATKri: uDY[5,4] + uDZ[5,4] - udz[5,4] - ustk[5,4] >= 0
 uATKri: uDY[5,4] - udz[5,4] >= 0
 uATKri: udy[5,5] - udz[5,5] = 0
 uATKri: uDY[5,5] - uDZ[5,5] + ustk[5,5] >= 0
 uATKri: - uDY[5,5] + uDZ[5,5] + ustk[5,5] >= 0
 uATKri: uDY[5,5] + uDZ[5,5] - udz[5,5] - ustk[5,5] >= 0
 uATKri: uDY[5,5] - udz[5,5] >= 0
 uATKri: udy[5,6] - udz[5,6] = 0
 uATKri: uDY[5,6] - uDZ[5,6] + ustk[5,6] >= 0
 uATKri: - uDY[5,6] + uDZ[5,6] + ustk[5,6] >= 0
 uATKri: uDY[5,6] + uDZ[5,6] - udz[5,6] - ustk[5,6] >= 0
 uATKri: uDY[5,6] - udz[5,6] >= 0
 uATKri: udy[5,7] - udz[5,7] = 0
 uATKri: uDY[5,7] - uDZ[5,7] + ustk[5,7] >= 0
 uATKri: - uDY[5,7] + uDZ[5,7] + ustk[5,7] >= 0
 uATKri: uDY[5,7] + uDZ[5,7] - udz[5,7] - ustk[5,7] >= 0
 uATKri: uDY[5,7] - udz[5,7] >= 0
 uATKri: uDY[5,8] - uDZ[5,8] = 0
 uATKri: udy[5,8] - udz[5,8] = 0
 uATKri: uDY[5,9] - uDZ[5,9] = 0
 uATKri: udy[5,9] - udz[5,9] = 0
 uATKri: uDY[5,10] - uDZ[5,10] = 0
 uATKri: udy[5,10] - udz[5,10] = 0
 uATKri: uDY[5,11] - uDZ[5,11] = 0
 uATKri: udy[5,11] - udz[5,11] = 0
 uATKri: uDY[5,12] - uDZ[5,12] = 0
 uATKri: udy[5,12] - udz[5,12] = 0
 uATKri: uDY[5,13] - uDZ[5,13] = 0
 uATKri: udy[5,13] - udz[5,13] = 0
 uATKri: uDY[5,14] - uDZ[5,14] = 0
 uATKri: udy[5,14] - udz[5,14] = 0
 uATKri: uDY[5,15] - uDZ[5,15] = 0
 uATKri: udy[5,15] - udz[5,15] = 0
 uATKr1: udy[6,0] - udz[6,0] = 0
 uATKr1: uDY[6,0] - uDZ[6,0] + ustk[33,0] >= 0
 uATKr1: - uDY[6,0] + uDZ[6,0] + ustk[33,0] >= 0
 uATKr1: uDY[6,0] + uDZ[6,0] - udz[6,0] - ustk[33,0] >= 0
 uATKr1: uDY[6,0] - udz[6,0] >= 0
 uATKr1: udy[6,1] - udz[6,1] = 0
 uATKr1: uDY[6,1] - uDZ[6,1] + ustk[33,1] >= 0
 uATKr1: - uDY[6,1] + uDZ[6,1] + ustk[33,1] >= 0
 uATKr1: uDY[6,1] + uDZ[6,1] - udz[6,1] - ustk[33,1] >= 0
 uATKr1: uDY[6,1] - udz[6,1] >= 0
 uATKr1: udy[6,2] - udz[6,2] = 0
 uATKr1: uDY[6,2] - uDZ[6,2] + ustk[33,2] >= 0
 uATKr1: - uDY[6,2] + uDZ[6,2] + ustk[33,2] >= 0
 uATKr1: uDY[6,2] + uDZ[6,2] - udz[6,2] - ustk[33,2] >= 0
 uATKr1: uDY[6,2] - udz[6,2] >= 0
 uATKr1: udy[6,3] - udz[6,3] = 0
 uATKr1: uDY[6,3] - uDZ[6,3] + ustk[33,3] >= 0
 uATKr1: - uDY[6,3] + uDZ[6,3] + ustk[33,3] >= 0
 uATKr1: uDY[6,3] + uDZ[6,3] - udz[6,3] - ustk[33,3] >= 0
 uATKr1: uDY[6,3] - udz[6,3] >= 0
 uATKr1: udy[6,4] - udz[6,4] = 0
 uATKr1: uDY[6,4] - uDZ[6,4] + ustk[33,4] >= 0
 uATKr1: - uDY[6,4] + uDZ[6,4] + ustk[33,4] >= 0
 uATKr1: uDY[6,4] + uDZ[6,4] - udz[6,4] - ustk[33,4] >= 0
 uATKr1: uDY[6,4] - udz[6,4] >= 0
 uATKr1: udy[6,5] - udz[6,5] = 0
 uATKr1: uDY[6,5] - uDZ[6,5] + ustk[33,5] >= 0
 uATKr1: - uDY[6,5] + uDZ[6,5] + ustk[33,5] >= 0
 uATKr1: uDY[6,5] + uDZ[6,5] - udz[6,5] - ustk[33,5] >= 0
 uATKr1: uDY[6,5] - udz[6,5] >= 0
 uATKr1: udy[6,6] - udz[6,6] = 0
 uATKr1: uDY[6,6] - uDZ[6,6] + ustk[33,6] >= 0
 uATKr1: - uDY[6,6] + uDZ[6,6] + ustk[33,6] >= 0
 uATKr1: uDY[6,6] + uDZ[6,6] - udz[6,6] - ustk[33,6] >= 0
 uATKr1: uDY[6,6] - udz[6,6] >= 0
 uATKr1: udy[6,7] - udz[6,7] = 0
 uATKr1: uDY[6,7] - uDZ[6,7] + ustk[33,7] >= 0
 uATKr1: - uDY[6,7] + uDZ[6,7] + ustk[33,7] >= 0
 uATKr1: uDY[6,7] + uDZ[6,7] - udz[6,7] - ustk[33,7] >= 0
 uATKr1: uDY[6,7] - udz[6,7] >= 0
 uATKr1: uDY[6,8] - uDZ[6,8] = 0
 uATKr1: udy[6,8] - udz[6,8] = 0
 uATKr1: uDY[6,9] - uDZ[6,9] = 0
 uATKr1: udy[6,9] - udz[6,9] = 0
 uATKr1: uDY[6,10] - uDZ[6,10] = 0
 uATKr1: udy[6,10] - udz[6,10] = 0
 uATKr1: uDY[6,11] - uDZ[6,11] = 0
 uATKr1: udy[6,11] - udz[6,11] = 0
 uATKr1: uDY[6,12] - uDZ[6,12] = 0
 uATKr1: udy[6,12] - udz[6,12] = 0
 uATKr1: uDY[6,13] - uDZ[6,13] = 0
 uATKr1: udy[6,13] - udz[6,13] = 0
 uATKr1: uDY[6,14] - uDZ[6,14] = 0
 uATKr1: udy[6,14] - udz[6,14] = 0
 uATKr1: uDY[6,15] - uDZ[6,15] = 0
 uATKr1: udy[6,15] - udz[6,15] = 0
 uATKr1: udy[7,0] - udz[7,0] = 0
 uATKr1: uDY[7,0] - uDZ[7,0] + ustk[34,0] >= 0
 uATKr1: - uDY[7,0] + uDZ[7,0] + ustk[34,0] >= 0
 uATKr1: uDY[7,0] + uDZ[7,0] - udz[7,0] - ustk[34,0] >= 0
 uATKr1: uDY[7,0] - udz[7,0] >= 0
 uATKr1: udy[7,1] - udz[7,1] = 0
 uATKr1: uDY[7,1] - uDZ[7,1] + ustk[34,1] >= 0
 uATKr1: - uDY[7,1] + uDZ[7,1] + ustk[34,1] >= 0
 uATKr1: uDY[7,1] + uDZ[7,1] - udz[7,1] - ustk[34,1] >= 0
 uATKr1: uDY[7,1] - udz[7,1] >= 0
 uATKr1: udy[7,2] - udz[7,2] = 0
 uATKr1: uDY[7,2] - uDZ[7,2] + ustk[34,2] >= 0
 uATKr1: - uDY[7,2] + uDZ[7,2] + ustk[34,2] >= 0
 uATKr1: uDY[7,2] + uDZ[7,2] - udz[7,2] - ustk[34,2] >= 0
 uATKr1: uDY[7,2] - udz[7,2] >= 0
 uATKr1: udy[7,3] - udz[7,3] = 0
 uATKr1: uDY[7,3] - uDZ[7,3] + ustk[34,3] >= 0
 uATKr1: - uDY[7,3] + uDZ[7,3] + ustk[34,3] >= 0
 uATKr1: uDY[7,3] + uDZ[7,3] - udz[7,3] - ustk[34,3] >= 0
 uATKr1: uDY[7,3] - udz[7,3] >= 0
 uATKr1: udy[7,4] - udz[7,4] = 0
 uATKr1: uDY[7,4] - uDZ[7,4] + ustk[34,4] >= 0
 uATKr1: - uDY[7,4] + uDZ[7,4] + ustk[34,4] >= 0
 uATKr1: uDY[7,4] + uDZ[7,4] - udz[7,4] - ustk[34,4] >= 0
 uATKr1: uDY[7,4] - udz[7,4] >= 0
 uATKr1: udy[7,5] - udz[7,5] = 0
 uATKr1: uDY[7,5] - uDZ[7,5] + ustk[34,5] >= 0
 uATKr1: - uDY[7,5] + uDZ[7,5] + ustk[34,5] >= 0
 uATKr1: uDY[7,5] + uDZ[7,5] - udz[7,5] - ustk[34,5] >= 0
 uATKr1: uDY[7,5] - udz[7,5] >= 0
 uATKr1: udy[7,6] - udz[7,6] = 0
 uATKr1: uDY[7,6] - uDZ[7,6] + ustk[34,6] >= 0
 uATKr1: - uDY[7,6] + uDZ[7,6] + ustk[34,6] >= 0
 uATKr1: uDY[7,6] + uDZ[7,6] - udz[7,6] - ustk[34,6] >= 0
 uATKr1: uDY[7,6] - udz[7,6] >= 0
 uATKr1: udy[7,7] - udz[7,7] = 0
 uATKr1: uDY[7,7] - uDZ[7,7] + ustk[34,7] >= 0
 uATKr1: - uDY[7,7] + uDZ[7,7] + ustk[34,7] >= 0
 uATKr1: uDY[7,7] + uDZ[7,7] - udz[7,7] - ustk[34,7] >= 0
 uATKr1: uDY[7,7] - udz[7,7] >= 0
 uATKr1: uDY[7,8] - uDZ[7,8] = 0
 uATKr1: udy[7,8] - udz[7,8] = 0
 uATKr1: uDY[7,9] - uDZ[7,9] = 0
 uATKr1: udy[7,9] - udz[7,9] = 0
 uATKr1: uDY[7,10] - uDZ[7,10] = 0
 uATKr1: udy[7,10] - udz[7,10] = 0
 uATKr1: uDY[7,11] - uDZ[7,11] = 0
 uATKr1: udy[7,11] - udz[7,11] = 0
 uATKr1: uDY[7,12] - uDZ[7,12] = 0
 uATKr1: udy[7,12] - udz[7,12] = 0
 uATKr1: uDY[7,13] - uDZ[7,13] = 0
 uATKr1: udy[7,13] - udz[7,13] = 0
 uATKr1: uDY[7,14] - uDZ[7,14] = 0
 uATKr1: udy[7,14] - udz[7,14] = 0
 uATKr1: uDY[7,15] - uDZ[7,15] = 0
 uATKr1: udy[7,15] - udz[7,15] = 0
 uATKr1: udy[8,0] - udz[8,0] = 0
 uATKr1: uDY[8,0] - uDZ[8,0] + ustk[35,0] >= 0
 uATKr1: - uDY[8,0] + uDZ[8,0] + ustk[35,0] >= 0
 uATKr1: uDY[8,0] + uDZ[8,0] - udz[8,0] - ustk[35,0] >= 0
 uATKr1: uDY[8,0] - udz[8,0] >= 0
 uATKr1: udy[8,1] - udz[8,1] = 0
 uATKr1: uDY[8,1] - uDZ[8,1] + ustk[35,1] >= 0
 uATKr1: - uDY[8,1] + uDZ[8,1] + ustk[35,1] >= 0
 uATKr1: uDY[8,1] + uDZ[8,1] - udz[8,1] - ustk[35,1] >= 0
 uATKr1: uDY[8,1] - udz[8,1] >= 0
 uATKr1: udy[8,2] - udz[8,2] = 0
 uATKr1: uDY[8,2] - uDZ[8,2] + ustk[35,2] >= 0
 uATKr1: - uDY[8,2] + uDZ[8,2] + ustk[35,2] >= 0
 uATKr1: uDY[8,2] + uDZ[8,2] - udz[8,2] - ustk[35,2] >= 0
 uATKr1: uDY[8,2] - udz[8,2] >= 0
 uATKr1: udy[8,3] - udz[8,3] = 0
 uATKr1: uDY[8,3] - uDZ[8,3] + ustk[35,3] >= 0
 uATKr1: - uDY[8,3] + uDZ[8,3] + ustk[35,3] >= 0
 uATKr1: uDY[8,3] + uDZ[8,3] - udz[8,3] - ustk[35,3] >= 0
 uATKr1: uDY[8,3] - udz[8,3] >= 0
 uATKr1: udy[8,4] - udz[8,4] = 0
 uATKr1: uDY[8,4] - uDZ[8,4] + ustk[35,4] >= 0
 uATKr1: - uDY[8,4] + uDZ[8,4] + ustk[35,4] >= 0
 uATKr1: uDY[8,4] + uDZ[8,4] - udz[8,4] - ustk[35,4] >= 0
 uATKr1: uDY[8,4] - udz[8,4] >= 0
 uATKr1: udy[8,5] - udz[8,5] = 0
 uATKr1: uDY[8,5] - uDZ[8,5] + ustk[35,5] >= 0
 uATKr1: - uDY[8,5] + uDZ[8,5] + ustk[35,5] >= 0
 uATKr1: uDY[8,5] + uDZ[8,5] - udz[8,5] - ustk[35,5] >= 0
 uATKr1: uDY[8,5] - udz[8,5] >= 0
 uATKr1: udy[8,6] - udz[8,6] = 0
 uATKr1: uDY[8,6] - uDZ[8,6] + ustk[35,6] >= 0
 uATKr1: - uDY[8,6] + uDZ[8,6] + ustk[35,6] >= 0
 uATKr1: uDY[8,6] + uDZ[8,6] - udz[8,6] - ustk[35,6] >= 0
 uATKr1: uDY[8,6] - udz[8,6] >= 0
 uATKr1: udy[8,7] - udz[8,7] = 0
 uATKr1: uDY[8,7] - uDZ[8,7] + ustk[35,7] >= 0
 uATKr1: - uDY[8,7] + uDZ[8,7] + ustk[35,7] >= 0
 uATKr1: uDY[8,7] + uDZ[8,7] - udz[8,7] - ustk[35,7] >= 0
 uATKr1: uDY[8,7] - udz[8,7] >= 0
 uATKr1: uDY[8,8] - uDZ[8,8] = 0
 uATKr1: udy[8,8] - udz[8,8] = 0
 uATKr1: uDY[8,9] - uDZ[8,9] = 0
 uATKr1: udy[8,9] - udz[8,9] = 0
 uATKr1: uDY[8,10] - uDZ[8,10] = 0
 uATKr1: udy[8,10] - udz[8,10] = 0
 uATKr1: uDY[8,11] - uDZ[8,11] = 0
 uATKr1: udy[8,11] - udz[8,11] = 0
 uATKr1: uDY[8,12] - uDZ[8,12] = 0
 uATKr1: udy[8,12] - udz[8,12] = 0
 uATKr1: uDY[8,13] - uDZ[8,13] = 0
 uATKr1: udy[8,13] - udz[8,13] = 0
 uATKr1: uDY[8,14] - uDZ[8,14] = 0
 uATKr1: udy[8,14] - udz[8,14] = 0
 uATKr1: uDY[8,15] - uDZ[8,15] = 0
 uATKr1: udy[8,15] - udz[8,15] = 0
 uATKr1: udy[9,0] - udz[9,0] = 0
 uATKr1: uDY[9,0] - uDZ[9,0] + ustk[36,0] >= 0
 uATKr1: - uDY[9,0] + uDZ[9,0] + ustk[36,0] >= 0
 uATKr1: uDY[9,0] + uDZ[9,0] - udz[9,0] - ustk[36,0] >= 0
 uATKr1: uDY[9,0] - udz[9,0] >= 0
 uATKr1: udy[9,1] - udz[9,1] = 0
 uATKr1: uDY[9,1] - uDZ[9,1] + ustk[36,1] >= 0
 uATKr1: - uDY[9,1] + uDZ[9,1] + ustk[36,1] >= 0
 uATKr1: uDY[9,1] + uDZ[9,1] - udz[9,1] - ustk[36,1] >= 0
 uATKr1: uDY[9,1] - udz[9,1] >= 0
 uATKr1: udy[9,2] - udz[9,2] = 0
 uATKr1: uDY[9,2] - uDZ[9,2] + ustk[36,2] >= 0
 uATKr1: - uDY[9,2] + uDZ[9,2] + ustk[36,2] >= 0
 uATKr1: uDY[9,2] + uDZ[9,2] - udz[9,2] - ustk[36,2] >= 0
 uATKr1: uDY[9,2] - udz[9,2] >= 0
 uATKr1: udy[9,3] - udz[9,3] = 0
 uATKr1: uDY[9,3] - uDZ[9,3] + ustk[36,3] >= 0
 uATKr1: - uDY[9,3] + uDZ[9,3] + ustk[36,3] >= 0
 uATKr1: uDY[9,3] + uDZ[9,3] - udz[9,3] - ustk[36,3] >= 0
 uATKr1: uDY[9,3] - udz[9,3] >= 0
 uATKr1: udy[9,4] - udz[9,4] = 0
 uATKr1: uDY[9,4] - uDZ[9,4] + ustk[36,4] >= 0
 uATKr1: - uDY[9,4] + uDZ[9,4] + ustk[36,4] >= 0
 uATKr1: uDY[9,4] + uDZ[9,4] - udz[9,4] - ustk[36,4] >= 0
 uATKr1: uDY[9,4] - udz[9,4] >= 0
 uATKr1: udy[9,5] - udz[9,5] = 0
 uATKr1: uDY[9,5] - uDZ[9,5] + ustk[36,5] >= 0
 uATKr1: - uDY[9,5] + uDZ[9,5] + ustk[36,5] >= 0
 uATKr1: uDY[9,5] + uDZ[9,5] - udz[9,5] - ustk[36,5] >= 0
 uATKr1: uDY[9,5] - udz[9,5] >= 0
 uATKr1: udy[9,6] - udz[9,6] = 0
 uATKr1: uDY[9,6] - uDZ[9,6] + ustk[36,6] >= 0
 uATKr1: - uDY[9,6] + uDZ[9,6] + ustk[36,6] >= 0
 uATKr1: uDY[9,6] + uDZ[9,6] - udz[9,6] - ustk[36,6] >= 0
 uATKr1: uDY[9,6] - udz[9,6] >= 0
 uATKr1: udy[9,7] - udz[9,7] = 0
 uATKr1: uDY[9,7] - uDZ[9,7] + ustk[36,7] >= 0
 uATKr1: - uDY[9,7] + uDZ[9,7] + ustk[36,7] >= 0
 uATKr1: uDY[9,7] + uDZ[9,7] - udz[9,7] - ustk[36,7] >= 0
 uATKr1: uDY[9,7] - udz[9,7] >= 0
 uATKr1: uDY[9,8] - uDZ[9,8] = 0
 uATKr1: udy[9,8] - udz[9,8] = 0
 uATKr1: uDY[9,9] - uDZ[9,9] = 0
 uATKr1: udy[9,9] - udz[9,9] = 0
 uATKr1: uDY[9,10] - uDZ[9,10] = 0
 uATKr1: udy[9,10] - udz[9,10] = 0
 uATKr1: uDY[9,11] - uDZ[9,11] = 0
 uATKr1: udy[9,11] - udz[9,11] = 0
 uATKr1: uDY[9,12] - uDZ[9,12] = 0
 uATKr1: udy[9,12] - udz[9,12] = 0
 uATKr1: uDY[9,13] - uDZ[9,13] = 0
 uATKr1: udy[9,13] - udz[9,13] = 0
 uATKr1: uDY[9,14] - uDZ[9,14] = 0
 uATKr1: udy[9,14] - udz[9,14] = 0
 uATKr1: uDY[9,15] - uDZ[9,15] = 0
 uATKr1: udy[9,15] - udz[9,15] = 0
 uATKr1: udy[10,0] - udz[10,0] = 0
 uATKr1: uDY[10,0] - uDZ[10,0] + ustk[37,0] >= 0
 uATKr1: - uDY[10,0] + uDZ[10,0] + ustk[37,0] >= 0
 uATKr1: uDY[10,0] + uDZ[10,0] - udz[10,0] - ustk[37,0] >= 0
 uATKr1: uDY[10,0] - udz[10,0] >= 0
 uATKr1: udy[10,1] - udz[10,1] = 0
 uATKr1: uDY[10,1] - uDZ[10,1] + ustk[37,1] >= 0
 uATKr1: - uDY[10,1] + uDZ[10,1] + ustk[37,1] >= 0
 uATKr1: uDY[10,1] + uDZ[10,1] - udz[10,1] - ustk[37,1] >= 0
 uATKr1: uDY[10,1] - udz[10,1] >= 0
 uATKr1: udy[10,2] - udz[10,2] = 0
 uATKr1: uDY[10,2] - uDZ[10,2] + ustk[37,2] >= 0
 uATKr1: - uDY[10,2] + uDZ[10,2] + ustk[37,2] >= 0
 uATKr1: uDY[10,2] + uDZ[10,2] - udz[10,2] - ustk[37,2] >= 0
 uATKr1: uDY[10,2] - udz[10,2] >= 0
 uATKr1: udy[10,3] - udz[10,3] = 0
 uATKr1: uDY[10,3] - uDZ[10,3] + ustk[37,3] >= 0
 uATKr1: - uDY[10,3] + uDZ[10,3] + ustk[37,3] >= 0
 uATKr1: uDY[10,3] + uDZ[10,3] - udz[10,3] - ustk[37,3] >= 0
 uATKr1: uDY[10,3] - udz[10,3] >= 0
 uATKr1: udy[10,4] - udz[10,4] = 0
 uATKr1: uDY[10,4] - uDZ[10,4] + ustk[37,4] >= 0
 uATKr1: - uDY[10,4] + uDZ[10,4] + ustk[37,4] >= 0
 uATKr1: uDY[10,4] + uDZ[10,4] - udz[10,4] - ustk[37,4] >= 0
 uATKr1: uDY[10,4] - udz[10,4] >= 0
 uATKr1: udy[10,5] - udz[10,5] = 0
 uATKr1: uDY[10,5] - uDZ[10,5] + ustk[37,5] >= 0
 uATKr1: - uDY[10,5] + uDZ[10,5] + ustk[37,5] >= 0
 uATKr1: uDY[10,5] + uDZ[10,5] - udz[10,5] - ustk[37,5] >= 0
 uATKr1: uDY[10,5] - udz[10,5] >= 0
 uATKr1: udy[10,6] - udz[10,6] = 0
 uATKr1: uDY[10,6] - uDZ[10,6] + ustk[37,6] >= 0
 uATKr1: - uDY[10,6] + uDZ[10,6] + ustk[37,6] >= 0
 uATKr1: uDY[10,6] + uDZ[10,6] - udz[10,6] - ustk[37,6] >= 0
 uATKr1: uDY[10,6] - udz[10,6] >= 0
 uATKr1: udy[10,7] - udz[10,7] = 0
 uATKr1: uDY[10,7] - uDZ[10,7] + ustk[37,7] >= 0
 uATKr1: - uDY[10,7] + uDZ[10,7] + ustk[37,7] >= 0
 uATKr1: uDY[10,7] + uDZ[10,7] - udz[10,7] - ustk[37,7] >= 0
 uATKr1: uDY[10,7] - udz[10,7] >= 0
 uATKr1: uDY[10,8] - uDZ[10,8] = 0
 uATKr1: udy[10,8] - udz[10,8] = 0
 uATKr1: uDY[10,9] - uDZ[10,9] = 0
 uATKr1: udy[10,9] - udz[10,9] = 0
 uATKr1: uDY[10,10] - uDZ[10,10] = 0
 uATKr1: udy[10,10] - udz[10,10] = 0
 uATKr1: uDY[10,11] - uDZ[10,11] = 0
 uATKr1: udy[10,11] - udz[10,11] = 0
 uATKr1: uDY[10,12] - uDZ[10,12] = 0
 uATKr1: udy[10,12] - udz[10,12] = 0
 uATKr1: uDY[10,13] - uDZ[10,13] = 0
 uATKr1: udy[10,13] - udz[10,13] = 0
 uATKr1: uDY[10,14] - uDZ[10,14] = 0
 uATKr1: udy[10,14] - udz[10,14] = 0
 uATKr1: uDY[10,15] - uDZ[10,15] = 0
 uATKr1: udy[10,15] - udz[10,15] = 0
 uATKr1: udy[11,0] - udz[11,0] = 0
 uATKr1: uDY[11,0] - uDZ[11,0] + ustk[38,0] >= 0
 uATKr1: - uDY[11,0] + uDZ[11,0] + ustk[38,0] >= 0
 uATKr1: uDY[11,0] + uDZ[11,0] - udz[11,0] - ustk[38,0] >= 0
 uATKr1: uDY[11,0] - udz[11,0] >= 0
 uATKr1: udy[11,1] - udz[11,1] = 0
 uATKr1: uDY[11,1] - uDZ[11,1] + ustk[38,1] >= 0
 uATKr1: - uDY[11,1] + uDZ[11,1] + ustk[38,1] >= 0
 uATKr1: uDY[11,1] + uDZ[11,1] - udz[11,1] - ustk[38,1] >= 0
 uATKr1: uDY[11,1] - udz[11,1] >= 0
 uATKr1: udy[11,2] - udz[11,2] = 0
 uATKr1: uDY[11,2] - uDZ[11,2] + ustk[38,2] >= 0
 uATKr1: - uDY[11,2] + uDZ[11,2] + ustk[38,2] >= 0
 uATKr1: uDY[11,2] + uDZ[11,2] - udz[11,2] - ustk[38,2] >= 0
 uATKr1: uDY[11,2] - udz[11,2] >= 0
 uATKr1: udy[11,3] - udz[11,3] = 0
 uATKr1: uDY[11,3] - uDZ[11,3] + ustk[38,3] >= 0
 uATKr1: - uDY[11,3] + uDZ[11,3] + ustk[38,3] >= 0
 uATKr1: uDY[11,3] + uDZ[11,3] - udz[11,3] - ustk[38,3] >= 0
 uATKr1: uDY[11,3] - udz[11,3] >= 0
 uATKr1: udy[11,4] - udz[11,4] = 0
 uATKr1: uDY[11,4] - uDZ[11,4] + ustk[38,4] >= 0
 uATKr1: - uDY[11,4] + uDZ[11,4] + ustk[38,4] >= 0
 uATKr1: uDY[11,4] + uDZ[11,4] - udz[11,4] - ustk[38,4] >= 0
 uATKr1: uDY[11,4] - udz[11,4] >= 0
 uATKr1: udy[11,5] - udz[11,5] = 0
 uATKr1: uDY[11,5] - uDZ[11,5] + ustk[38,5] >= 0
 uATKr1: - uDY[11,5] + uDZ[11,5] + ustk[38,5] >= 0
 uATKr1: uDY[11,5] + uDZ[11,5] - udz[11,5] - ustk[38,5] >= 0
 uATKr1: uDY[11,5] - udz[11,5] >= 0
 uATKr1: udy[11,6] - udz[11,6] = 0
 uATKr1: uDY[11,6] - uDZ[11,6] + ustk[38,6] >= 0
 uATKr1: - uDY[11,6] + uDZ[11,6] + ustk[38,6] >= 0
 uATKr1: uDY[11,6] + uDZ[11,6] - udz[11,6] - ustk[38,6] >= 0
 uATKr1: uDY[11,6] - udz[11,6] >= 0
 uATKr1: udy[11,7] - udz[11,7] = 0
 uATKr1: uDY[11,7] - uDZ[11,7] + ustk[38,7] >= 0
 uATKr1: - uDY[11,7] + uDZ[11,7] + ustk[38,7] >= 0
 uATKr1: uDY[11,7] + uDZ[11,7] - udz[11,7] - ustk[38,7] >= 0
 uATKr1: uDY[11,7] - udz[11,7] >= 0
 uATKr1: uDY[11,8] - uDZ[11,8] = 0
 uATKr1: udy[11,8] - udz[11,8] = 0
 uATKr1: uDY[11,9] - uDZ[11,9] = 0
 uATKr1: udy[11,9] - udz[11,9] = 0
 uATKr1: uDY[11,10] - uDZ[11,10] = 0
 uATKr1: udy[11,10] - udz[11,10] = 0
 uATKr1: uDY[11,11] - uDZ[11,11] = 0
 uATKr1: udy[11,11] - udz[11,11] = 0
 uATKr1: uDY[11,12] - uDZ[11,12] = 0
 uATKr1: udy[11,12] - udz[11,12] = 0
 uATKr1: uDY[11,13] - uDZ[11,13] = 0
 uATKr1: udy[11,13] - udz[11,13] = 0
 uATKr1: uDY[11,14] - uDZ[11,14] = 0
 uATKr1: udy[11,14] - udz[11,14] = 0
 uATKr1: uDY[11,15] - uDZ[11,15] = 0
 uATKr1: udy[11,15] - udz[11,15] = 0
 uSR: uDZ[1,0] - uDW[1,0] = 0
 uSR: udz[1,0] - udw[1,0] = 0
 uSR: uDZ[1,1] - uDW[1,1] = 0
 uSR: udz[1,1] - udw[1,1] = 0
 uSR: uDZ[1,2] - uDW[1,2] = 0
 uSR: udz[1,2] - udw[1,2] = 0
 uSR: uDZ[1,3] - uDW[1,3] = 0
 uSR: udz[1,3] - udw[1,3] = 0
 uSR: uDZ[1,4] - uDW[1,5] = 0
 uSR: udz[1,4] - udw[1,5] = 0
 uSR: uDZ[1,5] - uDW[1,6] = 0
 uSR: udz[1,5] - udw[1,6] = 0
 uSR: uDZ[1,6] - uDW[1,7] = 0
 uSR: udz[1,6] - udw[1,7] = 0
 uSR: uDZ[1,7] - uDW[1,4] = 0
 uSR: udz[1,7] - udw[1,4] = 0
 uSR: uDZ[1,8] - uDW[1,10] = 0
 uSR: udz[1,8] - udw[1,10] = 0
 uSR: uDZ[1,9] - uDW[1,11] = 0
 uSR: udz[1,9] - udw[1,11] = 0
 uSR: uDZ[1,10] - uDW[1,8] = 0
 uSR: udz[1,10] - udw[1,8] = 0
 uSR: uDZ[1,11] - uDW[1,9] = 0
 uSR: udz[1,11] - udw[1,9] = 0
 uSR: uDZ[1,12] - uDW[1,15] = 0
 uSR: udz[1,12] - udw[1,15] = 0
 uSR: uDZ[1,13] - uDW[1,12] = 0
 uSR: udz[1,13] - udw[1,12] = 0
 uSR: uDZ[1,14] - uDW[1,13] = 0
 uSR: udz[1,14] - udw[1,13] = 0
 uSR: uDZ[1,15] - uDW[1,14] = 0
 uSR: udz[1,15] - udw[1,14] = 0
 uSR: uDZ[2,0] - uDW[2,0] = 0
 uSR: udz[2,0] - udw[2,0] = 0
 uSR: uDZ[2,1] - uDW[2,1] = 0
 uSR: udz[2,1] - udw[2,1] = 0
 uSR: uDZ[2,2] - uDW[2,2] = 0
 uSR: udz[2,2] - udw[2,2] = 0
 uSR: uDZ[2,3] - uDW[2,3] = 0
 uSR: udz[2,3] - udw[2,3] = 0
 uSR: uDZ[2,4] - uDW[2,5] = 0
 uSR: udz[2,4] - udw[2,5] = 0
 uSR: uDZ[2,5] - uDW[2,6] = 0
 uSR: udz[2,5] - udw[2,6] = 0
 uSR: uDZ[2,6] - uDW[2,7] = 0
 uSR: udz[2,6] - udw[2,7] = 0
 uSR: uDZ[2,7] - uDW[2,4] = 0
 uSR: udz[2,7] - udw[2,4] = 0
 uSR: uDZ[2,8] - uDW[2,10] = 0
 uSR: udz[2,8] - udw[2,10] = 0
 uSR: uDZ[2,9] - uDW[2,11] = 0
 uSR: udz[2,9] - udw[2,11] = 0
 uSR: uDZ[2,10] - uDW[2,8] = 0
 uSR: udz[2,10] - udw[2,8] = 0
 uSR: uDZ[2,11] - uDW[2,9] = 0
 uSR: udz[2,11] - udw[2,9] = 0
 uSR: uDZ[2,12] - uDW[2,15] = 0
 uSR: udz[2,12] - udw[2,15] = 0
 uSR: uDZ[2,13] - uDW[2,12] = 0
 uSR: udz[2,13] - udw[2,12] = 0
 uSR: uDZ[2,14] - uDW[2,13] = 0
 uSR: udz[2,14] - udw[2,13] = 0
 uSR: uDZ[2,15] - uDW[2,14] = 0
 uSR: udz[2,15] - udw[2,14] = 0
 uSR: uDZ[3,0] - uDW[3,0] = 0
 uSR: udz[3,0] - udw[3,0] = 0
 uSR: uDZ[3,1] - uDW[3,1] = 0
 uSR: udz[3,1] - udw[3,1] = 0
 uSR: uDZ[3,2] - uDW[3,2] = 0
 uSR: udz[3,2] - udw[3,2] = 0
 uSR: uDZ[3,3] - uDW[3,3] = 0
 uSR: udz[3,3] - udw[3,3] = 0
 uSR: uDZ[3,4] - uDW[3,5] = 0
 uSR: udz[3,4] - udw[3,5] = 0
 uSR: uDZ[3,5] - uDW[3,6] = 0
 uSR: udz[3,5] - udw[3,6] = 0
 uSR: uDZ[3,6] - uDW[3,7] = 0
 uSR: udz[3,6] - udw[3,7] = 0
 uSR: uDZ[3,7] - uDW[3,4] = 0
 uSR: udz[3,7] - udw[3,4] = 0
 uSR: uDZ[3,8] - uDW[3,10] = 0
 uSR: udz[3,8] - udw[3,10] = 0
 uSR: uDZ[3,9] - uDW[3,11] = 0
 uSR: udz[3,9] - udw[3,11] = 0
 uSR: uDZ[3,10] - uDW[3,8] = 0
 uSR: udz[3,10] - udw[3,8] = 0
 uSR: uDZ[3,11] - uDW[3,9] = 0
 uSR: udz[3,11] - udw[3,9] = 0
 uSR: uDZ[3,12] - uDW[3,15] = 0
 uSR: udz[3,12] - udw[3,15] = 0
 uSR: uDZ[3,13] - uDW[3,12] = 0
 uSR: udz[3,13] - udw[3,12] = 0
 uSR: uDZ[3,14] - uDW[3,13] = 0
 uSR: udz[3,14] - udw[3,13] = 0
 uSR: uDZ[3,15] - uDW[3,14] = 0
 uSR: udz[3,15] - udw[3,14] = 0
 uSR: uDZ[4,0] - uDW[4,0] = 0
 uSR: udz[4,0] - udw[4,0] = 0
 uSR: uDZ[4,1] - uDW[4,1] = 0
 uSR: udz[4,1] - udw[4,1] = 0
 uSR: uDZ[4,2] - uDW[4,2] = 0
 uSR: udz[4,2] - udw[4,2] = 0
 uSR: uDZ[4,3] - uDW[4,3] = 0
 uSR: udz[4,3] - udw[4,3] = 0
 uSR: uDZ[4,4] - uDW[4,5] = 0
 uSR: udz[4,4] - udw[4,5] = 0
 uSR: uDZ[4,5] - uDW[4,6] = 0
 uSR: udz[4,5] - udw[4,6] = 0
 uSR: uDZ[4,6] - uDW[4,7] = 0
 uSR: udz[4,6] - udw[4,7] = 0
 uSR: uDZ[4,7] - uDW[4,4] = 0
 uSR: udz[4,7] - udw[4,4] = 0
 uSR: uDZ[4,8] - uDW[4,10] = 0
 uSR: udz[4,8] - udw[4,10] = 0
 uSR: uDZ[4,9] - uDW[4,11] = 0
 uSR: udz[4,9] - udw[4,11] = 0
 uSR: uDZ[4,10] - uDW[4,8] = 0
 uSR: udz[4,10] - udw[4,8] = 0
 uSR: uDZ[4,11] - uDW[4,9] = 0
 uSR: udz[4,11] - udw[4,9] = 0
 uSR: uDZ[4,12] - uDW[4,15] = 0
 uSR: udz[4,12] - udw[4,15] = 0
 uSR: uDZ[4,13] - uDW[4,12] = 0
 uSR: udz[4,13] - udw[4,12] = 0
 uSR: uDZ[4,14] - uDW[4,13] = 0
 uSR: udz[4,14] - udw[4,13] = 0
 uSR: uDZ[4,15] - uDW[4,14] = 0
 uSR: udz[4,15] - udw[4,14] = 0
 uSR: uDZ[5,0] - uDW[5,0] = 0
 uSR: udz[5,0] - udw[5,0] = 0
 uSR: uDZ[5,1] - uDW[5,1] = 0
 uSR: udz[5,1] - udw[5,1] = 0
 uSR: uDZ[5,2] - uDW[5,2] = 0
 uSR: udz[5,2] - udw[5,2] = 0
 uSR: uDZ[5,3] - uDW[5,3] = 0
 uSR: udz[5,3] - udw[5,3] = 0
 uSR: uDZ[5,4] - uDW[5,5] = 0
 uSR: udz[5,4] - udw[5,5] = 0
 uSR: uDZ[5,5] - uDW[5,6] = 0
 uSR: udz[5,5] - udw[5,6] = 0
 uSR: uDZ[5,6] - uDW[5,7] = 0
 uSR: udz[5,6] - udw[5,7] = 0
 uSR: uDZ[5,7] - uDW[5,4] = 0
 uSR: udz[5,7] - udw[5,4] = 0
 uSR: uDZ[5,8] - uDW[5,10] = 0
 uSR: udz[5,8] - udw[5,10] = 0
 uSR: uDZ[5,9] - uDW[5,11] = 0
 uSR: udz[5,9] - udw[5,11] = 0
 uSR: uDZ[5,10] - uDW[5,8] = 0
 uSR: udz[5,10] - udw[5,8] = 0
 uSR: uDZ[5,11] - uDW[5,9] = 0
 uSR: udz[5,11] - udw[5,9] = 0
 uSR: uDZ[5,12] - uDW[5,15] = 0
 uSR: udz[5,12] - udw[5,15] = 0
 uSR: uDZ[5,13] - uDW[5,12] = 0
 uSR: udz[5,13] - udw[5,12] = 0
 uSR: uDZ[5,14] - uDW[5,13] = 0
 uSR: udz[5,14] - udw[5,13] = 0
 uSR: uDZ[5,15] - uDW[5,14] = 0
 uSR: udz[5,15] - udw[5,14] = 0
 uSR: uDZ[6,0] - uDW[6,0] = 0
 uSR: udz[6,0] - udw[6,0] = 0
 uSR: uDZ[6,1] - uDW[6,1] = 0
 uSR: udz[6,1] - udw[6,1] = 0
 uSR: uDZ[6,2] - uDW[6,2] = 0
 uSR: udz[6,2] - udw[6,2] = 0
 uSR: uDZ[6,3] - uDW[6,3] = 0
 uSR: udz[6,3] - udw[6,3] = 0
 uSR: uDZ[6,4] - uDW[6,5] = 0
 uSR: udz[6,4] - udw[6,5] = 0
 uSR: uDZ[6,5] - uDW[6,6] = 0
 uSR: udz[6,5] - udw[6,6] = 0
 uSR: uDZ[6,6] - uDW[6,7] = 0
 uSR: udz[6,6] - udw[6,7] = 0
 uSR: uDZ[6,7] - uDW[6,4] = 0
 uSR: udz[6,7] - udw[6,4] = 0
 uSR: uDZ[6,8] - uDW[6,10] = 0
 uSR: udz[6,8] - udw[6,10] = 0
 uSR: uDZ[6,9] - uDW[6,11] = 0
 uSR: udz[6,9] - udw[6,11] = 0
 uSR: uDZ[6,10] - uDW[6,8] = 0
 uSR: udz[6,10] - udw[6,8] = 0
 uSR: uDZ[6,11] - uDW[6,9] = 0
 uSR: udz[6,11] - udw[6,9] = 0
 uSR: uDZ[6,12] - uDW[6,15] = 0
 uSR: udz[6,12] - udw[6,15] = 0
 uSR: uDZ[6,13] - uDW[6,12] = 0
 uSR: udz[6,13] - udw[6,12] = 0
 uSR: uDZ[6,14] - uDW[6,13] = 0
 uSR: udz[6,14] - udw[6,13] = 0
 uSR: uDZ[6,15] - uDW[6,14] = 0
 uSR: udz[6,15] - udw[6,14] = 0
 uSR: uDZ[7,0] - uDW[7,0] = 0
 uSR: udz[7,0] - udw[7,0] = 0
 uSR: uDZ[7,1] - uDW[7,1] = 0
 uSR: udz[7,1] - udw[7,1] = 0
 uSR: uDZ[7,2] - uDW[7,2] = 0
 uSR: udz[7,2] - udw[7,2] = 0
 uSR: uDZ[7,3] - uDW[7,3] = 0
 uSR: udz[7,3] - udw[7,3] = 0
 uSR: uDZ[7,4] - uDW[7,5] = 0
 uSR: udz[7,4] - udw[7,5] = 0
 uSR: uDZ[7,5] - uDW[7,6] = 0
 uSR: udz[7,5] - udw[7,6] = 0
 uSR: uDZ[7,6] - uDW[7,7] = 0
 uSR: udz[7,6] - udw[7,7] = 0
 uSR: uDZ[7,7] - uDW[7,4] = 0
 uSR: udz[7,7] - udw[7,4] = 0
 uSR: uDZ[7,8] - uDW[7,10] = 0
 uSR: udz[7,8] - udw[7,10] = 0
 uSR: uDZ[7,9] - uDW[7,11] = 0
 uSR: udz[7,9] - udw[7,11] = 0
 uSR: uDZ[7,10] - uDW[7,8] = 0
 uSR: udz[7,10] - udw[7,8] = 0
 uSR: uDZ[7,11] - uDW[7,9] = 0
 uSR: udz[7,11] - udw[7,9] = 0
 uSR: uDZ[7,12] - uDW[7,15] = 0
 uSR: udz[7,12] - udw[7,15] = 0
 uSR: uDZ[7,13] - uDW[7,12] = 0
 uSR: udz[7,13] - udw[7,12] = 0
 uSR: uDZ[7,14] - uDW[7,13] = 0
 uSR: udz[7,14] - udw[7,13] = 0
 uSR: uDZ[7,15] - uDW[7,14] = 0
 uSR: udz[7,15] - udw[7,14] = 0
 uSR: uDZ[8,0] - uDW[8,0] = 0
 uSR: udz[8,0] - udw[8,0] = 0
 uSR: uDZ[8,1] - uDW[8,1] = 0
 uSR: udz[8,1] - udw[8,1] = 0
 uSR: uDZ[8,2] - uDW[8,2] = 0
 uSR: udz[8,2] - udw[8,2] = 0
 uSR: uDZ[8,3] - uDW[8,3] = 0
 uSR: udz[8,3] - udw[8,3] = 0
 uSR: uDZ[8,4] - uDW[8,5] = 0
 uSR: udz[8,4] - udw[8,5] = 0
 uSR: uDZ[8,5] - uDW[8,6] = 0
 uSR: udz[8,5] - udw[8,6] = 0
 uSR: uDZ[8,6] - uDW[8,7] = 0
 uSR: udz[8,6] - udw[8,7] = 0
 uSR: uDZ[8,7] - uDW[8,4] = 0
 uSR: udz[8,7] - udw[8,4] = 0
 uSR: uDZ[8,8] - uDW[8,10] = 0
 uSR: udz[8,8] - udw[8,10] = 0
 uSR: uDZ[8,9] - uDW[8,11] = 0
 uSR: udz[8,9] - udw[8,11] = 0
 uSR: uDZ[8,10] - uDW[8,8] = 0
 uSR: udz[8,10] - udw[8,8] = 0
 uSR: uDZ[8,11] - uDW[8,9] = 0
 uSR: udz[8,11] - udw[8,9] = 0
 uSR: uDZ[8,12] - uDW[8,15] = 0
 uSR: udz[8,12] - udw[8,15] = 0
 uSR: uDZ[8,13] - uDW[8,12] = 0
 uSR: udz[8,13] - udw[8,12] = 0
 uSR: uDZ[8,14] - uDW[8,13] = 0
 uSR: udz[8,14] - udw[8,13] = 0
 uSR: uDZ[8,15] - uDW[8,14] = 0
 uSR: udz[8,15] - udw[8,14] = 0
 uSR: uDZ[9,0] - uDW[9,0] = 0
 uSR: udz[9,0] - udw[9,0] = 0
 uSR: uDZ[9,1] - uDW[9,1] = 0
 uSR: udz[9,1] - udw[9,1] = 0
 uSR: uDZ[9,2] - uDW[9,2] = 0
 uSR: udz[9,2] - udw[9,2] = 0
 uSR: uDZ[9,3] - uDW[9,3] = 0
 uSR: udz[9,3] - udw[9,3] = 0
 uSR: uDZ[9,4] - uDW[9,5] = 0
 uSR: udz[9,4] - udw[9,5] = 0
 uSR: uDZ[9,5] - uDW[9,6] = 0
 uSR: udz[9,5] - udw[9,6] = 0
 uSR: uDZ[9,6] - uDW[9,7] = 0
 uSR: udz[9,6] - udw[9,7] = 0
 uSR: uDZ[9,7] - uDW[9,4] = 0
 uSR: udz[9,7] - udw[9,4] = 0
 uSR: uDZ[9,8] - uDW[9,10] = 0
 uSR: udz[9,8] - udw[9,10] = 0
 uSR: uDZ[9,9] - uDW[9,11] = 0
 uSR: udz[9,9] - udw[9,11] = 0
 uSR: uDZ[9,10] - uDW[9,8] = 0
 uSR: udz[9,10] - udw[9,8] = 0
 uSR: uDZ[9,11] - uDW[9,9] = 0
 uSR: udz[9,11] - udw[9,9] = 0
 uSR: uDZ[9,12] - uDW[9,15] = 0
 uSR: udz[9,12] - udw[9,15] = 0
 uSR: uDZ[9,13] - uDW[9,12] = 0
 uSR: udz[9,13] - udw[9,12] = 0
 uSR: uDZ[9,14] - uDW[9,13] = 0
 uSR: udz[9,14] - udw[9,13] = 0
 uSR: uDZ[9,15] - uDW[9,14] = 0
 uSR: udz[9,15] - udw[9,14] = 0
 uSR: uDZ[10,0] - uDW[10,0] = 0
 uSR: udz[10,0] - udw[10,0] = 0
 uSR: uDZ[10,1] - uDW[10,1] = 0
 uSR: udz[10,1] - udw[10,1] = 0
 uSR: uDZ[10,2] - uDW[10,2] = 0
 uSR: udz[10,2] - udw[10,2] = 0
 uSR: uDZ[10,3] - uDW[10,3] = 0
 uSR: udz[10,3] - udw[10,3] = 0
 uSR: uDZ[10,4] - uDW[10,5] = 0
 uSR: udz[10,4] - udw[10,5] = 0
 uSR: uDZ[10,5] - uDW[10,6] = 0
 uSR: udz[10,5] - udw[10,6] = 0
 uSR: uDZ[10,6] - uDW[10,7] = 0
 uSR: udz[10,6] - udw[10,7] = 0
 uSR: uDZ[10,7] - uDW[10,4] = 0
 uSR: udz[10,7] - udw[10,4] = 0
 uSR: uDZ[10,8] - uDW[10,10] = 0
 uSR: udz[10,8] - udw[10,10] = 0
 uSR: uDZ[10,9] - uDW[10,11] = 0
 uSR: udz[10,9] - udw[10,11] = 0
 uSR: uDZ[10,10] - uDW[10,8] = 0
 uSR: udz[10,10] - udw[10,8] = 0
 uSR: uDZ[10,11] - uDW[10,9] = 0
 uSR: udz[10,11] - udw[10,9] = 0
 uSR: uDZ[10,12] - uDW[10,15] = 0
 uSR: udz[10,12] - udw[10,15] = 0
 uSR: uDZ[10,13] - uDW[10,12] = 0
 uSR: udz[10,13] - udw[10,12] = 0
 uSR: uDZ[10,14] - uDW[10,13] = 0
 uSR: udz[10,14] - udw[10,13] = 0
 uSR: uDZ[10,15] - uDW[10,14] = 0
 uSR: udz[10,15] - udw[10,14] = 0
 uSR: uDZ[11,0] - uDW[11,0] = 0
 uSR: udz[11,0] - udw[11,0] = 0
 uSR: uDZ[11,1] - uDW[11,1] = 0
 uSR: udz[11,1] - udw[11,1] = 0
 uSR: uDZ[11,2] - uDW[11,2] = 0
 uSR: udz[11,2] - udw[11,2] = 0
 uSR: uDZ[11,3] - uDW[11,3] = 0
 uSR: udz[11,3] - udw[11,3] = 0
 uSR: uDZ[11,4] - uDW[11,5] = 0
 uSR: udz[11,4] - udw[11,5] = 0
 uSR: uDZ[11,5] - uDW[11,6] = 0
 uSR: udz[11,5] - udw[11,6] = 0
 uSR: uDZ[11,6] - uDW[11,7] = 0
 uSR: udz[11,6] - udw[11,7] = 0
 uSR: uDZ[11,7] - uDW[11,4] = 0
 uSR: udz[11,7] - udw[11,4] = 0
 uSR: uDZ[11,8] - uDW[11,10] = 0
 uSR: udz[11,8] - udw[11,10] = 0
 uSR: uDZ[11,9] - uDW[11,11] = 0
 uSR: udz[11,9] - udw[11,11] = 0
 uSR: uDZ[11,10] - uDW[11,8] = 0
 uSR: udz[11,10] - udw[11,8] = 0
 uSR: uDZ[11,11] - uDW[11,9] = 0
 uSR: udz[11,11] - udw[11,9] = 0
 uSR: uDZ[11,12] - uDW[11,15] = 0
 uSR: udz[11,12] - udw[11,15] = 0
 uSR: uDZ[11,13] - uDW[11,12] = 0
 uSR: udz[11,13] - udw[11,12] = 0
 uSR: uDZ[11,14] - uDW[11,13] = 0
 uSR: udz[11,14] - udw[11,13] = 0
 uSR: uDZ[11,15] - uDW[11,14] = 0
 uSR: udz[11,15] - udw[11,14] = 0
 ueMC0: - udx[2,4] + udw[1,0] = 0
 ueMC0: - uDX[2,4] + uDW[1,0] = 0
 ueMC1: udw[1,4] - udw[1,8] >= 0
 ueMC1: - udx[2,8] + udw[1,4] >= 0
 ueMC1: - udx[2,8] + udw[1,4] - udw[1,8] <= 0
 ueMC1: uDX[2,8] - uDW[1,4] + uDW[1,8] >= 0
 ueMC1: uDX[2,8] + uDW[1,4] - uDW[1,8] >= 0
 ueMC1: - uDX[2,8] + uDW[1,4] + uDW[1,8] >= 0
 ueMC1: - uDX[2,8] + uDW[1,4] - uDW[1,8] <= 0
 ueMC2: - udx[2,4] + udw[1,8] >= 0
 ueMC2: - udx[2,12] + udw[1,8] >= 0
 ueMC2: - udx[2,4] - udx[2,12] + udw[1,8] <= 0
 ueMC2: uDX[2,4] + uDX[2,12] - uDW[1,8] >= 0
 ueMC2: - uDX[2,4] + uDX[2,12] + uDW[1,8] >= 0
 ueMC2: uDX[2,4] - uDX[2,12] + uDW[1,8] >= 0
 ueMC2: - uDX[2,4] - uDX[2,12] + uDW[1,8] <= 0
 ueMC3: - udx[2,0] + udw[1,12] >= 0
 ueMC3: - udx[2,12] + udw[1,12] >= 0
 ueMC3: - udx[2,0] - udx[2,12] + udw[1,12] <= 0
 ueMC3: uDX[2,0] + uDX[2,12] - uDW[1,12] >= 0
 ueMC3: - uDX[2,0] + uDX[2,12] + uDW[1,12] >= 0
 ueMC3: uDX[2,0] - uDX[2,12] + uDW[1,12] >= 0
 ueMC3: - uDX[2,0] - uDX[2,12] + uDW[1,12] <= 0
 ueMC0: - udx[2,5] + udw[1,1] = 0
 ueMC0: - uDX[2,5] + uDW[1,1] = 0
 ueMC1: udw[1,5] - udw[1,9] >= 0
 ueMC1: - udx[2,9] + udw[1,5] >= 0
 ueMC1: - udx[2,9] + udw[1,5] - udw[1,9] <= 0
 ueMC1: uDX[2,9] - uDW[1,5] + uDW[1,9] >= 0
 ueMC1: uDX[2,9] + uDW[1,5] - uDW[1,9] >= 0
 ueMC1: - uDX[2,9] + uDW[1,5] + uDW[1,9] >= 0
 ueMC1: - uDX[2,9] + uDW[1,5] - uDW[1,9] <= 0
 ueMC2: - udx[2,5] + udw[1,9] >= 0
 ueMC2: - udx[2,13] + udw[1,9] >= 0
 ueMC2: - udx[2,5] - udx[2,13] + udw[1,9] <= 0
 ueMC2: uDX[2,5] + uDX[2,13] - uDW[1,9] >= 0
 ueMC2: - uDX[2,5] + uDX[2,13] + uDW[1,9] >= 0
 ueMC2: uDX[2,5] - uDX[2,13] + uDW[1,9] >= 0
 ueMC2: - uDX[2,5] - uDX[2,13] + uDW[1,9] <= 0
 ueMC3: - udx[2,1] + udw[1,13] >= 0
 ueMC3: - udx[2,13] + udw[1,13] >= 0
 ueMC3: - udx[2,1] - udx[2,13] + udw[1,13] <= 0
 ueMC3: uDX[2,1] + uDX[2,13] - uDW[1,13] >= 0
 ueMC3: - uDX[2,1] + uDX[2,13] + uDW[1,13] >= 0
 ueMC3: uDX[2,1] - uDX[2,13] + uDW[1,13] >= 0
 ueMC3: - uDX[2,1] - uDX[2,13] + uDW[1,13] <= 0
 ueMC0: - udx[2,6] + udw[1,2] = 0
 ueMC0: - uDX[2,6] + uDW[1,2] = 0
 ueMC1: udw[1,6] - udw[1,10] >= 0
 ueMC1: - udx[2,10] + udw[1,6] >= 0
 ueMC1: - udx[2,10] + udw[1,6] - udw[1,10] <= 0
 ueMC1: uDX[2,10] - uDW[1,6] + uDW[1,10] >= 0
 ueMC1: uDX[2,10] + uDW[1,6] - uDW[1,10] >= 0
 ueMC1: - uDX[2,10] + uDW[1,6] + uDW[1,10] >= 0
 ueMC1: - uDX[2,10] + uDW[1,6] - uDW[1,10] <= 0
 ueMC2: - udx[2,6] + udw[1,10] >= 0
 ueMC2: - udx[2,14] + udw[1,10] >= 0
 ueMC2: - udx[2,6] - udx[2,14] + udw[1,10] <= 0
 ueMC2: uDX[2,6] + uDX[2,14] - uDW[1,10] >= 0
 ueMC2: - uDX[2,6] + uDX[2,14] + uDW[1,10] >= 0
 ueMC2: uDX[2,6] - uDX[2,14] + uDW[1,10] >= 0
 ueMC2: - uDX[2,6] - uDX[2,14] + uDW[1,10] <= 0
 ueMC3: - udx[2,2] + udw[1,14] >= 0
 ueMC3: - udx[2,14] + udw[1,14] >= 0
 ueMC3: - udx[2,2] - udx[2,14] + udw[1,14] <= 0
 ueMC3: uDX[2,2] + uDX[2,14] - uDW[1,14] >= 0
 ueMC3: - uDX[2,2] + uDX[2,14] + uDW[1,14] >= 0
 ueMC3: uDX[2,2] - uDX[2,14] + uDW[1,14] >= 0
 ueMC3: - uDX[2,2] - uDX[2,14] + uDW[1,14] <= 0
 ueMC0: - udx[2,7] + udw[1,3] = 0
 ueMC0: - uDX[2,7] + uDW[1,3] = 0
 ueMC1: udw[1,7] - udw[1,11] >= 0
 ueMC1: - udx[2,11] + udw[1,7] >= 0
 ueMC1: - udx[2,11] + udw[1,7] - udw[1,11] <= 0
 ueMC1: uDX[2,11] - uDW[1,7] + uDW[1,11] >= 0
 ueMC1: uDX[2,11] + uDW[1,7] - uDW[1,11] >= 0
 ueMC1: - uDX[2,11] + uDW[1,7] + uDW[1,11] >= 0
 ueMC1: - uDX[2,11] + uDW[1,7] - uDW[1,11] <= 0
 ueMC2: - udx[2,7] + udw[1,11] >= 0
 ueMC2: - udx[2,15] + udw[1,11] >= 0
 ueMC2: - udx[2,7] - udx[2,15] + udw[1,11] <= 0
 ueMC2: uDX[2,7] + uDX[2,15] - uDW[1,11] >= 0
 ueMC2: - uDX[2,7] + uDX[2,15] + uDW[1,11] >= 0
 ueMC2: uDX[2,7] - uDX[2,15] + uDW[1,11] >= 0
 ueMC2: - uDX[2,7] - uDX[2,15] + uDW[1,11] <= 0
 ueMC3: - udx[2,3] + udw[1,15] >= 0
 ueMC3: - udx[2,15] + udw[1,15] >= 0
 ueMC3: - udx[2,3] - udx[2,15] + udw[1,15] <= 0
 ueMC3: uDX[2,3] + uDX[2,15] - uDW[1,15] >= 0
 ueMC3: - uDX[2,3] + uDX[2,15] + uDW[1,15] >= 0
 ueMC3: uDX[2,3] - uDX[2,15] + uDW[1,15] >= 0
 ueMC3: - uDX[2,3] - uDX[2,15] + uDW[1,15] <= 0
 ueMC0: - udx[3,4] + udw[2,0] = 0
 ueMC0: - uDX[3,4] + uDW[2,0] = 0
 ueMC1: udw[2,4] - udw[2,8] >= 0
 ueMC1: - udx[3,8] + udw[2,4] >= 0
 ueMC1: - udx[3,8] + udw[2,4] - udw[2,8] <= 0
 ueMC1: uDX[3,8] - uDW[2,4] + uDW[2,8] >= 0
 ueMC1: uDX[3,8] + uDW[2,4] - uDW[2,8] >= 0
 ueMC1: - uDX[3,8] + uDW[2,4] + uDW[2,8] >= 0
 ueMC1: - uDX[3,8] + uDW[2,4] - uDW[2,8] <= 0
 ueMC2: - udx[3,4] + udw[2,8] >= 0
 ueMC2: - udx[3,12] + udw[2,8] >= 0
 ueMC2: - udx[3,4] - udx[3,12] + udw[2,8] <= 0
 ueMC2: uDX[3,4] + uDX[3,12] - uDW[2,8] >= 0
 ueMC2: - uDX[3,4] + uDX[3,12] + uDW[2,8] >= 0
 ueMC2: uDX[3,4] - uDX[3,12] + uDW[2,8] >= 0
 ueMC2: - uDX[3,4] - uDX[3,12] + uDW[2,8] <= 0
 ueMC3: - udx[3,0] + udw[2,12] >= 0
 ueMC3: - udx[3,12] + udw[2,12] >= 0
 ueMC3: - udx[3,0] - udx[3,12] + udw[2,12] <= 0
 ueMC3: uDX[3,0] + uDX[3,12] - uDW[2,12] >= 0
 ueMC3: - uDX[3,0] + uDX[3,12] + uDW[2,12] >= 0
 ueMC3: uDX[3,0] - uDX[3,12] + uDW[2,12] >= 0
 ueMC3: - uDX[3,0] - uDX[3,12] + uDW[2,12] <= 0
 ueMC0: - udx[3,5] + udw[2,1] = 0
 ueMC0: - uDX[3,5] + uDW[2,1] = 0
 ueMC1: udw[2,5] - udw[2,9] >= 0
 ueMC1: - udx[3,9] + udw[2,5] >= 0
 ueMC1: - udx[3,9] + udw[2,5] - udw[2,9] <= 0
 ueMC1: uDX[3,9] - uDW[2,5] + uDW[2,9] >= 0
 ueMC1: uDX[3,9] + uDW[2,5] - uDW[2,9] >= 0
 ueMC1: - uDX[3,9] + uDW[2,5] + uDW[2,9] >= 0
 ueMC1: - uDX[3,9] + uDW[2,5] - uDW[2,9] <= 0
 ueMC2: - udx[3,5] + udw[2,9] >= 0
 ueMC2: - udx[3,13] + udw[2,9] >= 0
 ueMC2: - udx[3,5] - udx[3,13] + udw[2,9] <= 0
 ueMC2: uDX[3,5] + uDX[3,13] - uDW[2,9] >= 0
 ueMC2: - uDX[3,5] + uDX[3,13] + uDW[2,9] >= 0
 ueMC2: uDX[3,5] - uDX[3,13] + uDW[2,9] >= 0
 ueMC2: - uDX[3,5] - uDX[3,13] + uDW[2,9] <= 0
 ueMC3: - udx[3,1] + udw[2,13] >= 0
 ueMC3: - udx[3,13] + udw[2,13] >= 0
 ueMC3: - udx[3,1] - udx[3,13] + udw[2,13] <= 0
 ueMC3: uDX[3,1] + uDX[3,13] - uDW[2,13] >= 0
 ueMC3: - uDX[3,1] + uDX[3,13] + uDW[2,13] >= 0
 ueMC3: uDX[3,1] - uDX[3,13] + uDW[2,13] >= 0
 ueMC3: - uDX[3,1] - uDX[3,13] + uDW[2,13] <= 0
 ueMC0: - udx[3,6] + udw[2,2] = 0
 ueMC0: - uDX[3,6] + uDW[2,2] = 0
 ueMC1: udw[2,6] - udw[2,10] >= 0
 ueMC1: - udx[3,10] + udw[2,6] >= 0
 ueMC1: - udx[3,10] + udw[2,6] - udw[2,10] <= 0
 ueMC1: uDX[3,10] - uDW[2,6] + uDW[2,10] >= 0
 ueMC1: uDX[3,10] + uDW[2,6] - uDW[2,10] >= 0
 ueMC1: - uDX[3,10] + uDW[2,6] + uDW[2,10] >= 0
 ueMC1: - uDX[3,10] + uDW[2,6] - uDW[2,10] <= 0
 ueMC2: - udx[3,6] + udw[2,10] >= 0
 ueMC2: - udx[3,14] + udw[2,10] >= 0
 ueMC2: - udx[3,6] - udx[3,14] + udw[2,10] <= 0
 ueMC2: uDX[3,6] + uDX[3,14] - uDW[2,10] >= 0
 ueMC2: - uDX[3,6] + uDX[3,14] + uDW[2,10] >= 0
 ueMC2: uDX[3,6] - uDX[3,14] + uDW[2,10] >= 0
 ueMC2: - uDX[3,6] - uDX[3,14] + uDW[2,10] <= 0
 ueMC3: - udx[3,2] + udw[2,14] >= 0
 ueMC3: - udx[3,14] + udw[2,14] >= 0
 ueMC3: - udx[3,2] - udx[3,14] + udw[2,14] <= 0
 ueMC3: uDX[3,2] + uDX[3,14] - uDW[2,14] >= 0
 ueMC3: - uDX[3,2] + uDX[3,14] + uDW[2,14] >= 0
 ueMC3: uDX[3,2] - uDX[3,14] + uDW[2,14] >= 0
 ueMC3: - uDX[3,2] - uDX[3,14] + uDW[2,14] <= 0
 ueMC0: - udx[3,7] + udw[2,3] = 0
 ueMC0: - uDX[3,7] + uDW[2,3] = 0
 ueMC1: udw[2,7] - udw[2,11] >= 0
 ueMC1: - udx[3,11] + udw[2,7] >= 0
 ueMC1: - udx[3,11] + udw[2,7] - udw[2,11] <= 0
 ueMC1: uDX[3,11] - uDW[2,7] + uDW[2,11] >= 0
 ueMC1: uDX[3,11] + uDW[2,7] - uDW[2,11] >= 0
 ueMC1: - uDX[3,11] + uDW[2,7] + uDW[2,11] >= 0
 ueMC1: - uDX[3,11] + uDW[2,7] - uDW[2,11] <= 0
 ueMC2: - udx[3,7] + udw[2,11] >= 0
 ueMC2: - udx[3,15] + udw[2,11] >= 0
 ueMC2: - udx[3,7] - udx[3,15] + udw[2,11] <= 0
 ueMC2: uDX[3,7] + uDX[3,15] - uDW[2,11] >= 0
 ueMC2: - uDX[3,7] + uDX[3,15] + uDW[2,11] >= 0
 ueMC2: uDX[3,7] - uDX[3,15] + uDW[2,11] >= 0
 ueMC2: - uDX[3,7] - uDX[3,15] + uDW[2,11] <= 0
 ueMC3: - udx[3,3] + udw[2,15] >= 0
 ueMC3: - udx[3,15] + udw[2,15] >= 0
 ueMC3: - udx[3,3] - udx[3,15] + udw[2,15] <= 0
 ueMC3: uDX[3,3] + uDX[3,15] - uDW[2,15] >= 0
 ueMC3: - uDX[3,3] + uDX[3,15] + uDW[2,15] >= 0
 ueMC3: uDX[3,3] - uDX[3,15] + uDW[2,15] >= 0
 ueMC3: - uDX[3,3] - uDX[3,15] + uDW[2,15] <= 0
 ueMC0: - udx[4,4] + udw[3,0] = 0
 ueMC0: - uDX[4,4] + uDW[3,0] = 0
 ueMC1: udw[3,4] - udw[3,8] >= 0
 ueMC1: - udx[4,8] + udw[3,4] >= 0
 ueMC1: - udx[4,8] + udw[3,4] - udw[3,8] <= 0
 ueMC1: uDX[4,8] - uDW[3,4] + uDW[3,8] >= 0
 ueMC1: uDX[4,8] + uDW[3,4] - uDW[3,8] >= 0
 ueMC1: - uDX[4,8] + uDW[3,4] + uDW[3,8] >= 0
 ueMC1: - uDX[4,8] + uDW[3,4] - uDW[3,8] <= 0
 ueMC2: - udx[4,4] + udw[3,8] >= 0
 ueMC2: - udx[4,12] + udw[3,8] >= 0
 ueMC2: - udx[4,4] - udx[4,12] + udw[3,8] <= 0
 ueMC2: uDX[4,4] + uDX[4,12] - uDW[3,8] >= 0
 ueMC2: - uDX[4,4] + uDX[4,12] + uDW[3,8] >= 0
 ueMC2: uDX[4,4] - uDX[4,12] + uDW[3,8] >= 0
 ueMC2: - uDX[4,4] - uDX[4,12] + uDW[3,8] <= 0
 ueMC3: - udx[4,0] + udw[3,12] >= 0
 ueMC3: - udx[4,12] + udw[3,12] >= 0
 ueMC3: - udx[4,0] - udx[4,12] + udw[3,12] <= 0
 ueMC3: uDX[4,0] + uDX[4,12] - uDW[3,12] >= 0
 ueMC3: - uDX[4,0] + uDX[4,12] + uDW[3,12] >= 0
 ueMC3: uDX[4,0] - uDX[4,12] + uDW[3,12] >= 0
 ueMC3: - uDX[4,0] - uDX[4,12] + uDW[3,12] <= 0
 ueMC0: - udx[4,5] + udw[3,1] = 0
 ueMC0: - uDX[4,5] + uDW[3,1] = 0
 ueMC1: udw[3,5] - udw[3,9] >= 0
 ueMC1: - udx[4,9] + udw[3,5] >= 0
 ueMC1: - udx[4,9] + udw[3,5] - udw[3,9] <= 0
 ueMC1: uDX[4,9] - uDW[3,5] + uDW[3,9] >= 0
 ueMC1: uDX[4,9] + uDW[3,5] - uDW[3,9] >= 0
 ueMC1: - uDX[4,9] + uDW[3,5] + uDW[3,9] >= 0
 ueMC1: - uDX[4,9] + uDW[3,5] - uDW[3,9] <= 0
 ueMC2: - udx[4,5] + udw[3,9] >= 0
 ueMC2: - udx[4,13] + udw[3,9] >= 0
 ueMC2: - udx[4,5] - udx[4,13] + udw[3,9] <= 0
 ueMC2: uDX[4,5] + uDX[4,13] - uDW[3,9] >= 0
 ueMC2: - uDX[4,5] + uDX[4,13] + uDW[3,9] >= 0
 ueMC2: uDX[4,5] - uDX[4,13] + uDW[3,9] >= 0
 ueMC2: - uDX[4,5] - uDX[4,13] + uDW[3,9] <= 0
 ueMC3: - udx[4,1] + udw[3,13] >= 0
 ueMC3: - udx[4,13] + udw[3,13] >= 0
 ueMC3: - udx[4,1] - udx[4,13] + udw[3,13] <= 0
 ueMC3: uDX[4,1] + uDX[4,13] - uDW[3,13] >= 0
 ueMC3: - uDX[4,1] + uDX[4,13] + uDW[3,13] >= 0
 ueMC3: uDX[4,1] - uDX[4,13] + uDW[3,13] >= 0
 ueMC3: - uDX[4,1] - uDX[4,13] + uDW[3,13] <= 0
 ueMC0: - udx[4,6] + udw[3,2] = 0
 ueMC0: - uDX[4,6] + uDW[3,2] = 0
 ueMC1: udw[3,6] - udw[3,10] >= 0
 ueMC1: - udx[4,10] + udw[3,6] >= 0
 ueMC1: - udx[4,10] + udw[3,6] - udw[3,10] <= 0
 ueMC1: uDX[4,10] - uDW[3,6] + uDW[3,10] >= 0
 ueMC1: uDX[4,10] + uDW[3,6] - uDW[3,10] >= 0
 ueMC1: - uDX[4,10] + uDW[3,6] + uDW[3,10] >= 0
 ueMC1: - uDX[4,10] + uDW[3,6] - uDW[3,10] <= 0
 ueMC2: - udx[4,6] + udw[3,10] >= 0
 ueMC2: - udx[4,14] + udw[3,10] >= 0
 ueMC2: - udx[4,6] - udx[4,14] + udw[3,10] <= 0
 ueMC2: uDX[4,6] + uDX[4,14] - uDW[3,10] >= 0
 ueMC2: - uDX[4,6] + uDX[4,14] + uDW[3,10] >= 0
 ueMC2: uDX[4,6] - uDX[4,14] + uDW[3,10] >= 0
 ueMC2: - uDX[4,6] - uDX[4,14] + uDW[3,10] <= 0
 ueMC3: - udx[4,2] + udw[3,14] >= 0
 ueMC3: - udx[4,14] + udw[3,14] >= 0
 ueMC3: - udx[4,2] - udx[4,14] + udw[3,14] <= 0
 ueMC3: uDX[4,2] + uDX[4,14] - uDW[3,14] >= 0
 ueMC3: - uDX[4,2] + uDX[4,14] + uDW[3,14] >= 0
 ueMC3: uDX[4,2] - uDX[4,14] + uDW[3,14] >= 0
 ueMC3: - uDX[4,2] - uDX[4,14] + uDW[3,14] <= 0
 ueMC0: - udx[4,7] + udw[3,3] = 0
 ueMC0: - uDX[4,7] + uDW[3,3] = 0
 ueMC1: udw[3,7] - udw[3,11] >= 0
 ueMC1: - udx[4,11] + udw[3,7] >= 0
 ueMC1: - udx[4,11] + udw[3,7] - udw[3,11] <= 0
 ueMC1: uDX[4,11] - uDW[3,7] + uDW[3,11] >= 0
 ueMC1: uDX[4,11] + uDW[3,7] - uDW[3,11] >= 0
 ueMC1: - uDX[4,11] + uDW[3,7] + uDW[3,11] >= 0
 ueMC1: - uDX[4,11] + uDW[3,7] - uDW[3,11] <= 0
 ueMC2: - udx[4,7] + udw[3,11] >= 0
 ueMC2: - udx[4,15] + udw[3,11] >= 0
 ueMC2: - udx[4,7] - udx[4,15] + udw[3,11] <= 0
 ueMC2: uDX[4,7] + uDX[4,15] - uDW[3,11] >= 0
 ueMC2: - uDX[4,7] + uDX[4,15] + uDW[3,11] >= 0
 ueMC2: uDX[4,7] - uDX[4,15] + uDW[3,11] >= 0
 ueMC2: - uDX[4,7] - uDX[4,15] + uDW[3,11] <= 0
 ueMC3: - udx[4,3] + udw[3,15] >= 0
 ueMC3: - udx[4,15] + udw[3,15] >= 0
 ueMC3: - udx[4,3] - udx[4,15] + udw[3,15] <= 0
 ueMC3: uDX[4,3] + uDX[4,15] - uDW[3,15] >= 0
 ueMC3: - uDX[4,3] + uDX[4,15] + uDW[3,15] >= 0
 ueMC3: uDX[4,3] - uDX[4,15] + uDW[3,15] >= 0
 ueMC3: - uDX[4,3] - uDX[4,15] + uDW[3,15] <= 0
 udMC0: udx[5,0] - udw[4,12] >= 0
 udMC0: udx[5,0] - udx[5,12] >= 0
 udMC0: udx[5,0] - udx[5,12] - udw[4,12] <= 0
 udMC0: - uDX[5,0] + uDX[5,12] + uDW[4,12] >= 0
 udMC0: uDX[5,0] + uDX[5,12] - uDW[4,12] >= 0
 udMC0: uDX[5,0] - uDX[5,12] + uDW[4,12] >= 0
 udMC0: uDX[5,0] - uDX[5,12] - uDW[4,12] <= 0
 udMC1: udx[5,4] - udw[4,0] = 0
 udMC1: uDX[5,4] - uDW[4,0] = 0
 udMC2: udx[5,8] - udw[4,4] >= 0
 udMC2: udx[5,8] - udw[4,8] >= 0
 udMC2: udx[5,8] - udw[4,4] - udw[4,8] <= 0
 udMC2: - uDX[5,8] + uDW[4,4] + uDW[4,8] >= 0
 udMC2: uDX[5,8] - uDW[4,4] + uDW[4,8] >= 0
 udMC2: uDX[5,8] + uDW[4,4] - uDW[4,8] >= 0
 udMC2: uDX[5,8] - uDW[4,4] - uDW[4,8] <= 0
 udMC3: udx[5,12] - udw[4,0] >= 0
 udMC3: udx[5,12] - udw[4,8] >= 0
 udMC3: udx[5,12] - udw[4,0] - udw[4,8] <= 0
 udMC3: - uDX[5,12] + uDW[4,0] + uDW[4,8] >= 0
 udMC3: uDX[5,12] - uDW[4,0] + uDW[4,8] >= 0
 udMC3: uDX[5,12] + uDW[4,0] - uDW[4,8] >= 0
 udMC3: uDX[5,12] - uDW[4,0] - uDW[4,8] <= 0
 udMC0: udx[5,1] - udw[4,13] >= 0
 udMC0: udx[5,1] - udx[5,13] >= 0
 udMC0: udx[5,1] - udx[5,13] - udw[4,13] <= 0
 udMC0: - uDX[5,1] + uDX[5,13] + uDW[4,13] >= 0
 udMC0: uDX[5,1] + uDX[5,13] - uDW[4,13] >= 0
 udMC0: uDX[5,1] - uDX[5,13] + uDW[4,13] >= 0
 udMC0: uDX[5,1] - uDX[5,13] - uDW[4,13] <= 0
 udMC1: udx[5,5] - udw[4,1] = 0
 udMC1: uDX[5,5] - uDW[4,1] = 0
 udMC2: udx[5,9] - udw[4,5] >= 0
 udMC2: udx[5,9] - udw[4,9] >= 0
 udMC2: udx[5,9] - udw[4,5] - udw[4,9] <= 0
 udMC2: - uDX[5,9] + uDW[4,5] + uDW[4,9] >= 0
 udMC2: uDX[5,9] - uDW[4,5] + uDW[4,9] >= 0
 udMC2: uDX[5,9] + uDW[4,5] - uDW[4,9] >= 0
 udMC2: uDX[5,9] - uDW[4,5] - uDW[4,9] <= 0
 udMC3: udx[5,13] - udw[4,1] >= 0
 udMC3: udx[5,13] - udw[4,9] >= 0
 udMC3: udx[5,13] - udw[4,1] - udw[4,9] <= 0
 udMC3: - uDX[5,13] + uDW[4,1] + uDW[4,9] >= 0
 udMC3: uDX[5,13] - uDW[4,1] + uDW[4,9] >= 0
 udMC3: uDX[5,13] + uDW[4,1] - uDW[4,9] >= 0
 udMC3: uDX[5,13] - uDW[4,1] - uDW[4,9] <= 0
 udMC0: udx[5,2] - udw[4,14] >= 0
 udMC0: udx[5,2] - udx[5,14] >= 0
 udMC0: udx[5,2] - udx[5,14] - udw[4,14] <= 0
 udMC0: - uDX[5,2] + uDX[5,14] + uDW[4,14] >= 0
 udMC0: uDX[5,2] + uDX[5,14] - uDW[4,14] >= 0
 udMC0: uDX[5,2] - uDX[5,14] + uDW[4,14] >= 0
 udMC0: uDX[5,2] - uDX[5,14] - uDW[4,14] <= 0
 udMC1: udx[5,6] - udw[4,2] = 0
 udMC1: uDX[5,6] - uDW[4,2] = 0
 udMC2: udx[5,10] - udw[4,6] >= 0
 udMC2: udx[5,10] - udw[4,10] >= 0
 udMC2: udx[5,10] - udw[4,6] - udw[4,10] <= 0
 udMC2: - uDX[5,10] + uDW[4,6] + uDW[4,10] >= 0
 udMC2: uDX[5,10] - uDW[4,6] + uDW[4,10] >= 0
 udMC2: uDX[5,10] + uDW[4,6] - uDW[4,10] >= 0
 udMC2: uDX[5,10] - uDW[4,6] - uDW[4,10] <= 0
 udMC3: udx[5,14] - udw[4,2] >= 0
 udMC3: udx[5,14] - udw[4,10] >= 0
 udMC3: udx[5,14] - udw[4,2] - udw[4,10] <= 0
 udMC3: - uDX[5,14] + uDW[4,2] + uDW[4,10] >= 0
 udMC3: uDX[5,14] - uDW[4,2] + uDW[4,10] >= 0
 udMC3: uDX[5,14] + uDW[4,2] - uDW[4,10] >= 0
 udMC3: uDX[5,14] - uDW[4,2] - uDW[4,10] <= 0
 udMC0: udx[5,3] - udw[4,15] >= 0
 udMC0: udx[5,3] - udx[5,15] >= 0
 udMC0: udx[5,3] - udx[5,15] - udw[4,15] <= 0
 udMC0: - uDX[5,3] + uDX[5,15] + uDW[4,15] >= 0
 udMC0: uDX[5,3] + uDX[5,15] - uDW[4,15] >= 0
 udMC0: uDX[5,3] - uDX[5,15] + uDW[4,15] >= 0
 udMC0: uDX[5,3] - uDX[5,15] - uDW[4,15] <= 0
 udMC1: udx[5,7] - udw[4,3] = 0
 udMC1: uDX[5,7] - uDW[4,3] = 0
 udMC2: udx[5,11] - udw[4,7] >= 0
 udMC2: udx[5,11] - udw[4,11] >= 0
 udMC2: udx[5,11] - udw[4,7] - udw[4,11] <= 0
 udMC2: - uDX[5,11] + uDW[4,7] + uDW[4,11] >= 0
 udMC2: uDX[5,11] - uDW[4,7] + uDW[4,11] >= 0
 udMC2: uDX[5,11] + uDW[4,7] - uDW[4,11] >= 0
 udMC2: uDX[5,11] - uDW[4,7] - uDW[4,11] <= 0
 udMC3: udx[5,15] - udw[4,3] >= 0
 udMC3: udx[5,15] - udw[4,11] >= 0
 udMC3: udx[5,15] - udw[4,3] - udw[4,11] <= 0
 udMC3: - uDX[5,15] + uDW[4,3] + uDW[4,11] >= 0
 udMC3: uDX[5,15] - uDW[4,3] + uDW[4,11] >= 0
 udMC3: uDX[5,15] + uDW[4,3] - uDW[4,11] >= 0
 udMC3: uDX[5,15] - uDW[4,3] - uDW[4,11] <= 0
 udMC0: udx[6,0] - udw[5,12] >= 0
 udMC0: udx[6,0] - udx[6,12] >= 0
 udMC0: udx[6,0] - udx[6,12] - udw[5,12] <= 0
 udMC0: - uDX[6,0] + uDX[6,12] + uDW[5,12] >= 0
 udMC0: uDX[6,0] + uDX[6,12] - uDW[5,12] >= 0
 udMC0: uDX[6,0] - uDX[6,12] + uDW[5,12] >= 0
 udMC0: uDX[6,0] - uDX[6,12] - uDW[5,12] <= 0
 udMC1: udx[6,4] - udw[5,0] = 0
 udMC1: uDX[6,4] - uDW[5,0] = 0
 udMC2: udx[6,8] - udw[5,4] >= 0
 udMC2: udx[6,8] - udw[5,8] >= 0
 udMC2: udx[6,8] - udw[5,4] - udw[5,8] <= 0
 udMC2: - uDX[6,8] + uDW[5,4] + uDW[5,8] >= 0
 udMC2: uDX[6,8] - uDW[5,4] + uDW[5,8] >= 0
 udMC2: uDX[6,8] + uDW[5,4] - uDW[5,8] >= 0
 udMC2: uDX[6,8] - uDW[5,4] - uDW[5,8] <= 0
 udMC3: udx[6,12] - udw[5,0] >= 0
 udMC3: udx[6,12] - udw[5,8] >= 0
 udMC3: udx[6,12] - udw[5,0] - udw[5,8] <= 0
 udMC3: - uDX[6,12] + uDW[5,0] + uDW[5,8] >= 0
 udMC3: uDX[6,12] - uDW[5,0] + uDW[5,8] >= 0
 udMC3: uDX[6,12] + uDW[5,0] - uDW[5,8] >= 0
 udMC3: uDX[6,12] - uDW[5,0] - uDW[5,8] <= 0
 udMC0: udx[6,1] - udw[5,13] >= 0
 udMC0: udx[6,1] - udx[6,13] >= 0
 udMC0: udx[6,1] - udx[6,13] - udw[5,13] <= 0
 udMC0: - uDX[6,1] + uDX[6,13] + uDW[5,13] >= 0
 udMC0: uDX[6,1] + uDX[6,13] - uDW[5,13] >= 0
 udMC0: uDX[6,1] - uDX[6,13] + uDW[5,13] >= 0
 udMC0: uDX[6,1] - uDX[6,13] - uDW[5,13] <= 0
 udMC1: udx[6,5] - udw[5,1] = 0
 udMC1: uDX[6,5] - uDW[5,1] = 0
 udMC2: udx[6,9] - udw[5,5] >= 0
 udMC2: udx[6,9] - udw[5,9] >= 0
 udMC2: udx[6,9] - udw[5,5] - udw[5,9] <= 0
 udMC2: - uDX[6,9] + uDW[5,5] + uDW[5,9] >= 0
 udMC2: uDX[6,9] - uDW[5,5] + uDW[5,9] >= 0
 udMC2: uDX[6,9] + uDW[5,5] - uDW[5,9] >= 0
 udMC2: uDX[6,9] - uDW[5,5] - uDW[5,9] <= 0
 udMC3: udx[6,13] - udw[5,1] >= 0
 udMC3: udx[6,13] - udw[5,9] >= 0
 udMC3: udx[6,13] - udw[5,1] - udw[5,9] <= 0
 udMC3: - uDX[6,13] + uDW[5,1] + uDW[5,9] >= 0
 udMC3: uDX[6,13] - uDW[5,1] + uDW[5,9] >= 0
 udMC3: uDX[6,13] + uDW[5,1] - uDW[5,9] >= 0
 udMC3: uDX[6,13] - uDW[5,1] - uDW[5,9] <= 0
 udMC0: udx[6,2] - udw[5,14] >= 0
 udMC0: udx[6,2] - udx[6,14] >= 0
 udMC0: udx[6,2] - udx[6,14] - udw[5,14] <= 0
 udMC0: - uDX[6,2] + uDX[6,14] + uDW[5,14] >= 0
 udMC0: uDX[6,2] + uDX[6,14] - uDW[5,14] >= 0
 udMC0: uDX[6,2] - uDX[6,14] + uDW[5,14] >= 0
 udMC0: uDX[6,2] - uDX[6,14] - uDW[5,14] <= 0
 udMC1: udx[6,6] - udw[5,2] = 0
 udMC1: uDX[6,6] - uDW[5,2] = 0
 udMC2: udx[6,10] - udw[5,6] >= 0
 udMC2: udx[6,10] - udw[5,10] >= 0
 udMC2: udx[6,10] - udw[5,6] - udw[5,10] <= 0
 udMC2: - uDX[6,10] + uDW[5,6] + uDW[5,10] >= 0
 udMC2: uDX[6,10] - uDW[5,6] + uDW[5,10] >= 0
 udMC2: uDX[6,10] + uDW[5,6] - uDW[5,10] >= 0
 udMC2: uDX[6,10] - uDW[5,6] - uDW[5,10] <= 0
 udMC3: udx[6,14] - udw[5,2] >= 0
 udMC3: udx[6,14] - udw[5,10] >= 0
 udMC3: udx[6,14] - udw[5,2] - udw[5,10] <= 0
 udMC3: - uDX[6,14] + uDW[5,2] + uDW[5,10] >= 0
 udMC3: uDX[6,14] - uDW[5,2] + uDW[5,10] >= 0
 udMC3: uDX[6,14] + uDW[5,2] - uDW[5,10] >= 0
 udMC3: uDX[6,14] - uDW[5,2] - uDW[5,10] <= 0
 udMC0: udx[6,3] - udw[5,15] >= 0
 udMC0: udx[6,3] - udx[6,15] >= 0
 udMC0: udx[6,3] - udx[6,15] - udw[5,15] <= 0
 udMC0: - uDX[6,3] + uDX[6,15] + uDW[5,15] >= 0
 udMC0: uDX[6,3] + uDX[6,15] - uDW[5,15] >= 0
 udMC0: uDX[6,3] - uDX[6,15] + uDW[5,15] >= 0
 udMC0: uDX[6,3] - uDX[6,15] - uDW[5,15] <= 0
 udMC1: udx[6,7] - udw[5,3] = 0
 udMC1: uDX[6,7] - uDW[5,3] = 0
 udMC2: udx[6,11] - udw[5,7] >= 0
 udMC2: udx[6,11] - udw[5,11] >= 0
 udMC2: udx[6,11] - udw[5,7] - udw[5,11] <= 0
 udMC2: - uDX[6,11] + uDW[5,7] + uDW[5,11] >= 0
 udMC2: uDX[6,11] - uDW[5,7] + uDW[5,11] >= 0
 udMC2: uDX[6,11] + uDW[5,7] - uDW[5,11] >= 0
 udMC2: uDX[6,11] - uDW[5,7] - uDW[5,11] <= 0
 udMC3: udx[6,15] - udw[5,3] >= 0
 udMC3: udx[6,15] - udw[5,11] >= 0
 udMC3: udx[6,15] - udw[5,3] - udw[5,11] <= 0
 udMC3: - uDX[6,15] + uDW[5,3] + uDW[5,11] >= 0
 udMC3: uDX[6,15] - uDW[5,3] + uDW[5,11] >= 0
 udMC3: uDX[6,15] + uDW[5,3] - uDW[5,11] >= 0
 udMC3: uDX[6,15] - uDW[5,3] - uDW[5,11] <= 0
 udMC0: udx[7,0] - udw[6,12] >= 0
 udMC0: udx[7,0] - udx[7,12] >= 0
 udMC0: udx[7,0] - udx[7,12] - udw[6,12] <= 0
 udMC0: - uDX[7,0] + uDX[7,12] + uDW[6,12] >= 0
 udMC0: uDX[7,0] + uDX[7,12] - uDW[6,12] >= 0
 udMC0: uDX[7,0] - uDX[7,12] + uDW[6,12] >= 0
 udMC0: uDX[7,0] - uDX[7,12] - uDW[6,12] <= 0
 udMC1: udx[7,4] - udw[6,0] = 0
 udMC1: uDX[7,4] - uDW[6,0] = 0
 udMC2: udx[7,8] - udw[6,4] >= 0
 udMC2: udx[7,8] - udw[6,8] >= 0
 udMC2: udx[7,8] - udw[6,4] - udw[6,8] <= 0
 udMC2: - uDX[7,8] + uDW[6,4] + uDW[6,8] >= 0
 udMC2: uDX[7,8] - uDW[6,4] + uDW[6,8] >= 0
 udMC2: uDX[7,8] + uDW[6,4] - uDW[6,8] >= 0
 udMC2: uDX[7,8] - uDW[6,4] - uDW[6,8] <= 0
 udMC3: udx[7,12] - udw[6,0] >= 0
 udMC3: udx[7,12] - udw[6,8] >= 0
 udMC3: udx[7,12] - udw[6,0] - udw[6,8] <= 0
 udMC3: - uDX[7,12] + uDW[6,0] + uDW[6,8] >= 0
 udMC3: uDX[7,12] - uDW[6,0] + uDW[6,8] >= 0
 udMC3: uDX[7,12] + uDW[6,0] - uDW[6,8] >= 0
 udMC3: uDX[7,12] - uDW[6,0] - uDW[6,8] <= 0
 udMC0: udx[7,1] - udw[6,13] >= 0
 udMC0: udx[7,1] - udx[7,13] >= 0
 udMC0: udx[7,1] - udx[7,13] - udw[6,13] <= 0
 udMC0: - uDX[7,1] + uDX[7,13] + uDW[6,13] >= 0
 udMC0: uDX[7,1] + uDX[7,13] - uDW[6,13] >= 0
 udMC0: uDX[7,1] - uDX[7,13] + uDW[6,13] >= 0
 udMC0: uDX[7,1] - uDX[7,13] - uDW[6,13] <= 0
 udMC1: udx[7,5] - udw[6,1] = 0
 udMC1: uDX[7,5] - uDW[6,1] = 0
 udMC2: udx[7,9] - udw[6,5] >= 0
 udMC2: udx[7,9] - udw[6,9] >= 0
 udMC2: udx[7,9] - udw[6,5] - udw[6,9] <= 0
 udMC2: - uDX[7,9] + uDW[6,5] + uDW[6,9] >= 0
 udMC2: uDX[7,9] - uDW[6,5] + uDW[6,9] >= 0
 udMC2: uDX[7,9] + uDW[6,5] - uDW[6,9] >= 0
 udMC2: uDX[7,9] - uDW[6,5] - uDW[6,9] <= 0
 udMC3: udx[7,13] - udw[6,1] >= 0
 udMC3: udx[7,13] - udw[6,9] >= 0
 udMC3: udx[7,13] - udw[6,1] - udw[6,9] <= 0
 udMC3: - uDX[7,13] + uDW[6,1] + uDW[6,9] >= 0
 udMC3: uDX[7,13] - uDW[6,1] + uDW[6,9] >= 0
 udMC3: uDX[7,13] + uDW[6,1] - uDW[6,9] >= 0
 udMC3: uDX[7,13] - uDW[6,1] - uDW[6,9] <= 0
 udMC0: udx[7,2] - udw[6,14] >= 0
 udMC0: udx[7,2] - udx[7,14] >= 0
 udMC0: udx[7,2] - udx[7,14] - udw[6,14] <= 0
 udMC0: - uDX[7,2] + uDX[7,14] + uDW[6,14] >= 0
 udMC0: uDX[7,2] + uDX[7,14] - uDW[6,14] >= 0
 udMC0: uDX[7,2] - uDX[7,14] + uDW[6,14] >= 0
 udMC0: uDX[7,2] - uDX[7,14] - uDW[6,14] <= 0
 udMC1: udx[7,6] - udw[6,2] = 0
 udMC1: uDX[7,6] - uDW[6,2] = 0
 udMC2: udx[7,10] - udw[6,6] >= 0
 udMC2: udx[7,10] - udw[6,10] >= 0
 udMC2: udx[7,10] - udw[6,6] - udw[6,10] <= 0
 udMC2: - uDX[7,10] + uDW[6,6] + uDW[6,10] >= 0
 udMC2: uDX[7,10] - uDW[6,6] + uDW[6,10] >= 0
 udMC2: uDX[7,10] + uDW[6,6] - uDW[6,10] >= 0
 udMC2: uDX[7,10] - uDW[6,6] - uDW[6,10] <= 0
 udMC3: udx[7,14] - udw[6,2] >= 0
 udMC3: udx[7,14] - udw[6,10] >= 0
 udMC3: udx[7,14] - udw[6,2] - udw[6,10] <= 0
 udMC3: - uDX[7,14] + uDW[6,2] + uDW[6,10] >= 0
 udMC3: uDX[7,14] - uDW[6,2] + uDW[6,10] >= 0
 udMC3: uDX[7,14] + uDW[6,2] - uDW[6,10] >= 0
 udMC3: uDX[7,14] - uDW[6,2] - uDW[6,10] <= 0
 udMC0: udx[7,3] - udw[6,15] >= 0
 udMC0: udx[7,3] - udx[7,15] >= 0
 udMC0: udx[7,3] - udx[7,15] - udw[6,15] <= 0
 udMC0: - uDX[7,3] + uDX[7,15] + uDW[6,15] >= 0
 udMC0: uDX[7,3] + uDX[7,15] - uDW[6,15] >= 0
 udMC0: uDX[7,3] - uDX[7,15] + uDW[6,15] >= 0
 udMC0: uDX[7,3] - uDX[7,15] - uDW[6,15] <= 0
 udMC1: udx[7,7] - udw[6,3] = 0
 udMC1: uDX[7,7] - uDW[6,3] = 0
 udMC2: udx[7,11] - udw[6,7] >= 0
 udMC2: udx[7,11] - udw[6,11] >= 0
 udMC2: udx[7,11] - udw[6,7] - udw[6,11] <= 0
 udMC2: - uDX[7,11] + uDW[6,7] + uDW[6,11] >= 0
 udMC2: uDX[7,11] - uDW[6,7] + uDW[6,11] >= 0
 udMC2: uDX[7,11] + uDW[6,7] - uDW[6,11] >= 0
 udMC2: uDX[7,11] - uDW[6,7] - uDW[6,11] <= 0
 udMC3: udx[7,15] - udw[6,3] >= 0
 udMC3: udx[7,15] - udw[6,11] >= 0
 udMC3: udx[7,15] - udw[6,3] - udw[6,11] <= 0
 udMC3: - uDX[7,15] + uDW[6,3] + uDW[6,11] >= 0
 udMC3: uDX[7,15] - uDW[6,3] + uDW[6,11] >= 0
 udMC3: uDX[7,15] + uDW[6,3] - uDW[6,11] >= 0
 udMC3: uDX[7,15] - uDW[6,3] - uDW[6,11] <= 0
 udMC0: udx[8,0] - udw[7,12] >= 0
 udMC0: udx[8,0] - udx[8,12] >= 0
 udMC0: udx[8,0] - udx[8,12] - udw[7,12] <= 0
 udMC0: - uDX[8,0] + uDX[8,12] + uDW[7,12] >= 0
 udMC0: uDX[8,0] + uDX[8,12] - uDW[7,12] >= 0
 udMC0: uDX[8,0] - uDX[8,12] + uDW[7,12] >= 0
 udMC0: uDX[8,0] - uDX[8,12] - uDW[7,12] <= 0
 udMC1: udx[8,4] - udw[7,0] = 0
 udMC1: uDX[8,4] - uDW[7,0] = 0
 udMC2: udx[8,8] - udw[7,4] >= 0
 udMC2: udx[8,8] - udw[7,8] >= 0
 udMC2: udx[8,8] - udw[7,4] - udw[7,8] <= 0
 udMC2: - uDX[8,8] + uDW[7,4] + uDW[7,8] >= 0
 udMC2: uDX[8,8] - uDW[7,4] + uDW[7,8] >= 0
 udMC2: uDX[8,8] + uDW[7,4] - uDW[7,8] >= 0
 udMC2: uDX[8,8] - uDW[7,4] - uDW[7,8] <= 0
 udMC3: udx[8,12] - udw[7,0] >= 0
 udMC3: udx[8,12] - udw[7,8] >= 0
 udMC3: udx[8,12] - udw[7,0] - udw[7,8] <= 0
 udMC3: - uDX[8,12] + uDW[7,0] + uDW[7,8] >= 0
 udMC3: uDX[8,12] - uDW[7,0] + uDW[7,8] >= 0
 udMC3: uDX[8,12] + uDW[7,0] - uDW[7,8] >= 0
 udMC3: uDX[8,12] - uDW[7,0] - uDW[7,8] <= 0
 udMC0: udx[8,1] - udw[7,13] >= 0
 udMC0: udx[8,1] - udx[8,13] >= 0
 udMC0: udx[8,1] - udx[8,13] - udw[7,13] <= 0
 udMC0: - uDX[8,1] + uDX[8,13] + uDW[7,13] >= 0
 udMC0: uDX[8,1] + uDX[8,13] - uDW[7,13] >= 0
 udMC0: uDX[8,1] - uDX[8,13] + uDW[7,13] >= 0
 udMC0: uDX[8,1] - uDX[8,13] - uDW[7,13] <= 0
 udMC1: udx[8,5] - udw[7,1] = 0
 udMC1: uDX[8,5] - uDW[7,1] = 0
 udMC2: udx[8,9] - udw[7,5] >= 0
 udMC2: udx[8,9] - udw[7,9] >= 0
 udMC2: udx[8,9] - udw[7,5] - udw[7,9] <= 0
 udMC2: - uDX[8,9] + uDW[7,5] + uDW[7,9] >= 0
 udMC2: uDX[8,9] - uDW[7,5] + uDW[7,9] >= 0
 udMC2: uDX[8,9] + uDW[7,5] - uDW[7,9] >= 0
 udMC2: uDX[8,9] - uDW[7,5] - uDW[7,9] <= 0
 udMC3: udx[8,13] - udw[7,1] >= 0
 udMC3: udx[8,13] - udw[7,9] >= 0
 udMC3: udx[8,13] - udw[7,1] - udw[7,9] <= 0
 udMC3: - uDX[8,13] + uDW[7,1] + uDW[7,9] >= 0
 udMC3: uDX[8,13] - uDW[7,1] + uDW[7,9] >= 0
 udMC3: uDX[8,13] + uDW[7,1] - uDW[7,9] >= 0
 udMC3: uDX[8,13] - uDW[7,1] - uDW[7,9] <= 0
 udMC0: udx[8,2] - udw[7,14] >= 0
 udMC0: udx[8,2] - udx[8,14] >= 0
 udMC0: udx[8,2] - udx[8,14] - udw[7,14] <= 0
 udMC0: - uDX[8,2] + uDX[8,14] + uDW[7,14] >= 0
 udMC0: uDX[8,2] + uDX[8,14] - uDW[7,14] >= 0
 udMC0: uDX[8,2] - uDX[8,14] + uDW[7,14] >= 0
 udMC0: uDX[8,2] - uDX[8,14] - uDW[7,14] <= 0
 udMC1: udx[8,6] - udw[7,2] = 0
 udMC1: uDX[8,6] - uDW[7,2] = 0
 udMC2: udx[8,10] - udw[7,6] >= 0
 udMC2: udx[8,10] - udw[7,10] >= 0
 udMC2: udx[8,10] - udw[7,6] - udw[7,10] <= 0
 udMC2: - uDX[8,10] + uDW[7,6] + uDW[7,10] >= 0
 udMC2: uDX[8,10] - uDW[7,6] + uDW[7,10] >= 0
 udMC2: uDX[8,10] + uDW[7,6] - uDW[7,10] >= 0
 udMC2: uDX[8,10] - uDW[7,6] - uDW[7,10] <= 0
 udMC3: udx[8,14] - udw[7,2] >= 0
 udMC3: udx[8,14] - udw[7,10] >= 0
 udMC3: udx[8,14] - udw[7,2] - udw[7,10] <= 0
 udMC3: - uDX[8,14] + uDW[7,2] + uDW[7,10] >= 0
 udMC3: uDX[8,14] - uDW[7,2] + uDW[7,10] >= 0
 udMC3: uDX[8,14] + uDW[7,2] - uDW[7,10] >= 0
 udMC3: uDX[8,14] - uDW[7,2] - uDW[7,10] <= 0
 udMC0: udx[8,3] - udw[7,15] >= 0
 udMC0: udx[8,3] - udx[8,15] >= 0
 udMC0: udx[8,3] - udx[8,15] - udw[7,15] <= 0
 udMC0: - uDX[8,3] + uDX[8,15] + uDW[7,15] >= 0
 udMC0: uDX[8,3] + uDX[8,15] - uDW[7,15] >= 0
 udMC0: uDX[8,3] - uDX[8,15] + uDW[7,15] >= 0
 udMC0: uDX[8,3] - uDX[8,15] - uDW[7,15] <= 0
 udMC1: udx[8,7] - udw[7,3] = 0
 udMC1: uDX[8,7] - uDW[7,3] = 0
 udMC2: udx[8,11] - udw[7,7] >= 0
 udMC2: udx[8,11] - udw[7,11] >= 0
 udMC2: udx[8,11] - udw[7,7] - udw[7,11] <= 0
 udMC2: - uDX[8,11] + uDW[7,7] + uDW[7,11] >= 0
 udMC2: uDX[8,11] - uDW[7,7] + uDW[7,11] >= 0
 udMC2: uDX[8,11] + uDW[7,7] - uDW[7,11] >= 0
 udMC2: uDX[8,11] - uDW[7,7] - uDW[7,11] <= 0
 udMC3: udx[8,15] - udw[7,3] >= 0
 udMC3: udx[8,15] - udw[7,11] >= 0
 udMC3: udx[8,15] - udw[7,3] - udw[7,11] <= 0
 udMC3: - uDX[8,15] + uDW[7,3] + uDW[7,11] >= 0
 udMC3: uDX[8,15] - uDW[7,3] + uDW[7,11] >= 0
 udMC3: uDX[8,15] + uDW[7,3] - uDW[7,11] >= 0
 udMC3: uDX[8,15] - uDW[7,3] - uDW[7,11] <= 0
 udMC0: udx[9,0] - udw[8,12] >= 0
 udMC0: udx[9,0] - udx[9,12] >= 0
 udMC0: udx[9,0] - udx[9,12] - udw[8,12] <= 0
 udMC0: - uDX[9,0] + uDX[9,12] + uDW[8,12] >= 0
 udMC0: uDX[9,0] + uDX[9,12] - uDW[8,12] >= 0
 udMC0: uDX[9,0] - uDX[9,12] + uDW[8,12] >= 0
 udMC0: uDX[9,0] - uDX[9,12] - uDW[8,12] <= 0
 udMC1: udx[9,4] - udw[8,0] = 0
 udMC1: uDX[9,4] - uDW[8,0] = 0
 udMC2: udx[9,8] - udw[8,4] >= 0
 udMC2: udx[9,8] - udw[8,8] >= 0
 udMC2: udx[9,8] - udw[8,4] - udw[8,8] <= 0
 udMC2: - uDX[9,8] + uDW[8,4] + uDW[8,8] >= 0
 udMC2: uDX[9,8] - uDW[8,4] + uDW[8,8] >= 0
 udMC2: uDX[9,8] + uDW[8,4] - uDW[8,8] >= 0
 udMC2: uDX[9,8] - uDW[8,4] - uDW[8,8] <= 0
 udMC3: udx[9,12] - udw[8,0] >= 0
 udMC3: udx[9,12] - udw[8,8] >= 0
 udMC3: udx[9,12] - udw[8,0] - udw[8,8] <= 0
 udMC3: - uDX[9,12] + uDW[8,0] + uDW[8,8] >= 0
 udMC3: uDX[9,12] - uDW[8,0] + uDW[8,8] >= 0
 udMC3: uDX[9,12] + uDW[8,0] - uDW[8,8] >= 0
 udMC3: uDX[9,12] - uDW[8,0] - uDW[8,8] <= 0
 udMC0: udx[9,1] - udw[8,13] >= 0
 udMC0: udx[9,1] - udx[9,13] >= 0
 udMC0: udx[9,1] - udx[9,13] - udw[8,13] <= 0
 udMC0: - uDX[9,1] + uDX[9,13] + uDW[8,13] >= 0
 udMC0: uDX[9,1] + uDX[9,13] - uDW[8,13] >= 0
 udMC0: uDX[9,1] - uDX[9,13] + uDW[8,13] >= 0
 udMC0: uDX[9,1] - uDX[9,13] - uDW[8,13] <= 0
 udMC1: udx[9,5] - udw[8,1] = 0
 udMC1: uDX[9,5] - uDW[8,1] = 0
 udMC2: udx[9,9] - udw[8,5] >= 0
 udMC2: udx[9,9] - udw[8,9] >= 0
 udMC2: udx[9,9] - udw[8,5] - udw[8,9] <= 0
 udMC2: - uDX[9,9] + uDW[8,5] + uDW[8,9] >= 0
 udMC2: uDX[9,9] - uDW[8,5] + uDW[8,9] >= 0
 udMC2: uDX[9,9] + uDW[8,5] - uDW[8,9] >= 0
 udMC2: uDX[9,9] - uDW[8,5] - uDW[8,9] <= 0
 udMC3: udx[9,13] - udw[8,1] >= 0
 udMC3: udx[9,13] - udw[8,9] >= 0
 udMC3: udx[9,13] - udw[8,1] - udw[8,9] <= 0
 udMC3: - uDX[9,13] + uDW[8,1] + uDW[8,9] >= 0
 udMC3: uDX[9,13] - uDW[8,1] + uDW[8,9] >= 0
 udMC3: uDX[9,13] + uDW[8,1] - uDW[8,9] >= 0
 udMC3: uDX[9,13] - uDW[8,1] - uDW[8,9] <= 0
 udMC0: udx[9,2] - udw[8,14] >= 0
 udMC0: udx[9,2] - udx[9,14] >= 0
 udMC0: udx[9,2] - udx[9,14] - udw[8,14] <= 0
 udMC0: - uDX[9,2] + uDX[9,14] + uDW[8,14] >= 0
 udMC0: uDX[9,2] + uDX[9,14] - uDW[8,14] >= 0
 udMC0: uDX[9,2] - uDX[9,14] + uDW[8,14] >= 0
 udMC0: uDX[9,2] - uDX[9,14] - uDW[8,14] <= 0
 udMC1: udx[9,6] - udw[8,2] = 0
 udMC1: uDX[9,6] - uDW[8,2] = 0
 udMC2: udx[9,10] - udw[8,6] >= 0
 udMC2: udx[9,10] - udw[8,10] >= 0
 udMC2: udx[9,10] - udw[8,6] - udw[8,10] <= 0
 udMC2: - uDX[9,10] + uDW[8,6] + uDW[8,10] >= 0
 udMC2: uDX[9,10] - uDW[8,6] + uDW[8,10] >= 0
 udMC2: uDX[9,10] + uDW[8,6] - uDW[8,10] >= 0
 udMC2: uDX[9,10] - uDW[8,6] - uDW[8,10] <= 0
 udMC3: udx[9,14] - udw[8,2] >= 0
 udMC3: udx[9,14] - udw[8,10] >= 0
 udMC3: udx[9,14] - udw[8,2] - udw[8,10] <= 0
 udMC3: - uDX[9,14] + uDW[8,2] + uDW[8,10] >= 0
 udMC3: uDX[9,14] - uDW[8,2] + uDW[8,10] >= 0
 udMC3: uDX[9,14] + uDW[8,2] - uDW[8,10] >= 0
 udMC3: uDX[9,14] - uDW[8,2] - uDW[8,10] <= 0
 udMC0: udx[9,3] - udw[8,15] >= 0
 udMC0: udx[9,3] - udx[9,15] >= 0
 udMC0: udx[9,3] - udx[9,15] - udw[8,15] <= 0
 udMC0: - uDX[9,3] + uDX[9,15] + uDW[8,15] >= 0
 udMC0: uDX[9,3] + uDX[9,15] - uDW[8,15] >= 0
 udMC0: uDX[9,3] - uDX[9,15] + uDW[8,15] >= 0
 udMC0: uDX[9,3] - uDX[9,15] - uDW[8,15] <= 0
 udMC1: udx[9,7] - udw[8,3] = 0
 udMC1: uDX[9,7] - uDW[8,3] = 0
 udMC2: udx[9,11] - udw[8,7] >= 0
 udMC2: udx[9,11] - udw[8,11] >= 0
 udMC2: udx[9,11] - udw[8,7] - udw[8,11] <= 0
 udMC2: - uDX[9,11] + uDW[8,7] + uDW[8,11] >= 0
 udMC2: uDX[9,11] - uDW[8,7] + uDW[8,11] >= 0
 udMC2: uDX[9,11] + uDW[8,7] - uDW[8,11] >= 0
 udMC2: uDX[9,11] - uDW[8,7] - uDW[8,11] <= 0
 udMC3: udx[9,15] - udw[8,3] >= 0
 udMC3: udx[9,15] - udw[8,11] >= 0
 udMC3: udx[9,15] - udw[8,3] - udw[8,11] <= 0
 udMC3: - uDX[9,15] + uDW[8,3] + uDW[8,11] >= 0
 udMC3: uDX[9,15] - uDW[8,3] + uDW[8,11] >= 0
 udMC3: uDX[9,15] + uDW[8,3] - uDW[8,11] >= 0
 udMC3: uDX[9,15] - uDW[8,3] - uDW[8,11] <= 0
 udMC0: udx[10,0] - udw[9,12] >= 0
 udMC0: udx[10,0] - udx[10,12] >= 0
 udMC0: udx[10,0] - udx[10,12] - udw[9,12] <= 0
 udMC0: - uDX[10,0] + uDX[10,12] + uDW[9,12] >= 0
 udMC0: uDX[10,0] + uDX[10,12] - uDW[9,12] >= 0
 udMC0: uDX[10,0] - uDX[10,12] + uDW[9,12] >= 0
 udMC0: uDX[10,0] - uDX[10,12] - uDW[9,12] <= 0
 udMC1: udx[10,4] - udw[9,0] = 0
 udMC1: uDX[10,4] - uDW[9,0] = 0
 udMC2: udx[10,8] - udw[9,4] >= 0
 udMC2: udx[10,8] - udw[9,8] >= 0
 udMC2: udx[10,8] - udw[9,4] - udw[9,8] <= 0
 udMC2: - uDX[10,8] + uDW[9,4] + uDW[9,8] >= 0
 udMC2: uDX[10,8] - uDW[9,4] + uDW[9,8] >= 0
 udMC2: uDX[10,8] + uDW[9,4] - uDW[9,8] >= 0
 udMC2: uDX[10,8] - uDW[9,4] - uDW[9,8] <= 0
 udMC3: udx[10,12] - udw[9,0] >= 0
 udMC3: udx[10,12] - udw[9,8] >= 0
 udMC3: udx[10,12] - udw[9,0] - udw[9,8] <= 0
 udMC3: - uDX[10,12] + uDW[9,0] + uDW[9,8] >= 0
 udMC3: uDX[10,12] - uDW[9,0] + uDW[9,8] >= 0
 udMC3: uDX[10,12] + uDW[9,0] - uDW[9,8] >= 0
 udMC3: uDX[10,12] - uDW[9,0] - uDW[9,8] <= 0
 udMC0: udx[10,1] - udw[9,13] >= 0
 udMC0: udx[10,1] - udx[10,13] >= 0
 udMC0: udx[10,1] - udx[10,13] - udw[9,13] <= 0
 udMC0: - uDX[10,1] + uDX[10,13] + uDW[9,13] >= 0
 udMC0: uDX[10,1] + uDX[10,13] - uDW[9,13] >= 0
 udMC0: uDX[10,1] - uDX[10,13] + uDW[9,13] >= 0
 udMC0: uDX[10,1] - uDX[10,13] - uDW[9,13] <= 0
 udMC1: udx[10,5] - udw[9,1] = 0
 udMC1: uDX[10,5] - uDW[9,1] = 0
 udMC2: udx[10,9] - udw[9,5] >= 0
 udMC2: udx[10,9] - udw[9,9] >= 0
 udMC2: udx[10,9] - udw[9,5] - udw[9,9] <= 0
 udMC2: - uDX[10,9] + uDW[9,5] + uDW[9,9] >= 0
 udMC2: uDX[10,9] - uDW[9,5] + uDW[9,9] >= 0
 udMC2: uDX[10,9] + uDW[9,5] - uDW[9,9] >= 0
 udMC2: uDX[10,9] - uDW[9,5] - uDW[9,9] <= 0
 udMC3: udx[10,13] - udw[9,1] >= 0
 udMC3: udx[10,13] - udw[9,9] >= 0
 udMC3: udx[10,13] - udw[9,1] - udw[9,9] <= 0
 udMC3: - uDX[10,13] + uDW[9,1] + uDW[9,9] >= 0
 udMC3: uDX[10,13] - uDW[9,1] + uDW[9,9] >= 0
 udMC3: uDX[10,13] + uDW[9,1] - uDW[9,9] >= 0
 udMC3: uDX[10,13] - uDW[9,1] - uDW[9,9] <= 0
 udMC0: udx[10,2] - udw[9,14] >= 0
 udMC0: udx[10,2] - udx[10,14] >= 0
 udMC0: udx[10,2] - udx[10,14] - udw[9,14] <= 0
 udMC0: - uDX[10,2] + uDX[10,14] + uDW[9,14] >= 0
 udMC0: uDX[10,2] + uDX[10,14] - uDW[9,14] >= 0
 udMC0: uDX[10,2] - uDX[10,14] + uDW[9,14] >= 0
 udMC0: uDX[10,2] - uDX[10,14] - uDW[9,14] <= 0
 udMC1: udx[10,6] - udw[9,2] = 0
 udMC1: uDX[10,6] - uDW[9,2] = 0
 udMC2: udx[10,10] - udw[9,6] >= 0
 udMC2: udx[10,10] - udw[9,10] >= 0
 udMC2: udx[10,10] - udw[9,6] - udw[9,10] <= 0
 udMC2: - uDX[10,10] + uDW[9,6] + uDW[9,10] >= 0
 udMC2: uDX[10,10] - uDW[9,6] + uDW[9,10] >= 0
 udMC2: uDX[10,10] + uDW[9,6] - uDW[9,10] >= 0
 udMC2: uDX[10,10] - uDW[9,6] - uDW[9,10] <= 0
 udMC3: udx[10,14] - udw[9,2] >= 0
 udMC3: udx[10,14] - udw[9,10] >= 0
 udMC3: udx[10,14] - udw[9,2] - udw[9,10] <= 0
 udMC3: - uDX[10,14] + uDW[9,2] + uDW[9,10] >= 0
 udMC3: uDX[10,14] - uDW[9,2] + uDW[9,10] >= 0
 udMC3: uDX[10,14] + uDW[9,2] - uDW[9,10] >= 0
 udMC3: uDX[10,14] - uDW[9,2] - uDW[9,10] <= 0
 udMC0: udx[10,3] - udw[9,15] >= 0
 udMC0: udx[10,3] - udx[10,15] >= 0
 udMC0: udx[10,3] - udx[10,15] - udw[9,15] <= 0
 udMC0: - uDX[10,3] + uDX[10,15] + uDW[9,15] >= 0
 udMC0: uDX[10,3] + uDX[10,15] - uDW[9,15] >= 0
 udMC0: uDX[10,3] - uDX[10,15] + uDW[9,15] >= 0
 udMC0: uDX[10,3] - uDX[10,15] - uDW[9,15] <= 0
 udMC1: udx[10,7] - udw[9,3] = 0
 udMC1: uDX[10,7] - uDW[9,3] = 0
 udMC2: udx[10,11] - udw[9,7] >= 0
 udMC2: udx[10,11] - udw[9,11] >= 0
 udMC2: udx[10,11] - udw[9,7] - udw[9,11] <= 0
 udMC2: - uDX[10,11] + uDW[9,7] + uDW[9,11] >= 0
 udMC2: uDX[10,11] - uDW[9,7] + uDW[9,11] >= 0
 udMC2: uDX[10,11] + uDW[9,7] - uDW[9,11] >= 0
 udMC2: uDX[10,11] - uDW[9,7] - uDW[9,11] <= 0
 udMC3: udx[10,15] - udw[9,3] >= 0
 udMC3: udx[10,15] - udw[9,11] >= 0
 udMC3: udx[10,15] - udw[9,3] - udw[9,11] <= 0
 udMC3: - uDX[10,15] + uDW[9,3] + uDW[9,11] >= 0
 udMC3: uDX[10,15] - uDW[9,3] + uDW[9,11] >= 0
 udMC3: uDX[10,15] + uDW[9,3] - uDW[9,11] >= 0
 udMC3: uDX[10,15] - uDW[9,3] - uDW[9,11] <= 0
 udMC0: udx[11,0] - udw[10,12] >= 0
 udMC0: udx[11,0] - udx[11,12] >= 0
 udMC0: udx[11,0] - udx[11,12] - udw[10,12] <= 0
 udMC0: - uDX[11,0] + uDX[11,12] + uDW[10,12] >= 0
 udMC0: uDX[11,0] + uDX[11,12] - uDW[10,12] >= 0
 udMC0: uDX[11,0] - uDX[11,12] + uDW[10,12] >= 0
 udMC0: uDX[11,0] - uDX[11,12] - uDW[10,12] <= 0
 udMC1: udx[11,4] - udw[10,0] = 0
 udMC1: uDX[11,4] - uDW[10,0] = 0
 udMC2: udx[11,8] - udw[10,4] >= 0
 udMC2: udx[11,8] - udw[10,8] >= 0
 udMC2: udx[11,8] - udw[10,4] - udw[10,8] <= 0
 udMC2: - uDX[11,8] + uDW[10,4] + uDW[10,8] >= 0
 udMC2: uDX[11,8] - uDW[10,4] + uDW[10,8] >= 0
 udMC2: uDX[11,8] + uDW[10,4] - uDW[10,8] >= 0
 udMC2: uDX[11,8] - uDW[10,4] - uDW[10,8] <= 0
 udMC3: udx[11,12] - udw[10,0] >= 0
 udMC3: udx[11,12] - udw[10,8] >= 0
 udMC3: udx[11,12] - udw[10,0] - udw[10,8] <= 0
 udMC3: - uDX[11,12] + uDW[10,0] + uDW[10,8] >= 0
 udMC3: uDX[11,12] - uDW[10,0] + uDW[10,8] >= 0
 udMC3: uDX[11,12] + uDW[10,0] - uDW[10,8] >= 0
 udMC3: uDX[11,12] - uDW[10,0] - uDW[10,8] <= 0
 udMC0: udx[11,1] - udw[10,13] >= 0
 udMC0: udx[11,1] - udx[11,13] >= 0
 udMC0: udx[11,1] - udx[11,13] - udw[10,13] <= 0
 udMC0: - uDX[11,1] + uDX[11,13] + uDW[10,13] >= 0
 udMC0: uDX[11,1] + uDX[11,13] - uDW[10,13] >= 0
 udMC0: uDX[11,1] - uDX[11,13] + uDW[10,13] >= 0
 udMC0: uDX[11,1] - uDX[11,13] - uDW[10,13] <= 0
 udMC1: udx[11,5] - udw[10,1] = 0
 udMC1: uDX[11,5] - uDW[10,1] = 0
 udMC2: udx[11,9] - udw[10,5] >= 0
 udMC2: udx[11,9] - udw[10,9] >= 0
 udMC2: udx[11,9] - udw[10,5] - udw[10,9] <= 0
 udMC2: - uDX[11,9] + uDW[10,5] + uDW[10,9] >= 0
 udMC2: uDX[11,9] - uDW[10,5] + uDW[10,9] >= 0
 udMC2: uDX[11,9] + uDW[10,5] - uDW[10,9] >= 0
 udMC2: uDX[11,9] - uDW[10,5] - uDW[10,9] <= 0
 udMC3: udx[11,13] - udw[10,1] >= 0
 udMC3: udx[11,13] - udw[10,9] >= 0
 udMC3: udx[11,13] - udw[10,1] - udw[10,9] <= 0
 udMC3: - uDX[11,13] + uDW[10,1] + uDW[10,9] >= 0
 udMC3: uDX[11,13] - uDW[10,1] + uDW[10,9] >= 0
 udMC3: uDX[11,13] + uDW[10,1] - uDW[10,9] >= 0
 udMC3: uDX[11,13] - uDW[10,1] - uDW[10,9] <= 0
 udMC0: udx[11,2] - udw[10,14] >= 0
 udMC0: udx[11,2] - udx[11,14] >= 0
 udMC0: udx[11,2] - udx[11,14] - udw[10,14] <= 0
 udMC0: - uDX[11,2] + uDX[11,14] + uDW[10,14] >= 0
 udMC0: uDX[11,2] + uDX[11,14] - uDW[10,14] >= 0
 udMC0: uDX[11,2] - uDX[11,14] + uDW[10,14] >= 0
 udMC0: uDX[11,2] - uDX[11,14] - uDW[10,14] <= 0
 udMC1: udx[11,6] - udw[10,2] = 0
 udMC1: uDX[11,6] - uDW[10,2] = 0
 udMC2: udx[11,10] - udw[10,6] >= 0
 udMC2: udx[11,10] - udw[10,10] >= 0
 udMC2: udx[11,10] - udw[10,6] - udw[10,10] <= 0
 udMC2: - uDX[11,10] + uDW[10,6] + uDW[10,10] >= 0
 udMC2: uDX[11,10] - uDW[10,6] + uDW[10,10] >= 0
 udMC2: uDX[11,10] + uDW[10,6] - uDW[10,10] >= 0
 udMC2: uDX[11,10] - uDW[10,6] - uDW[10,10] <= 0
 udMC3: udx[11,14] - udw[10,2] >= 0
 udMC3: udx[11,14] - udw[10,10] >= 0
 udMC3: udx[11,14] - udw[10,2] - udw[10,10] <= 0
 udMC3: - uDX[11,14] + uDW[10,2] + uDW[10,10] >= 0
 udMC3: uDX[11,14] - uDW[10,2] + uDW[10,10] >= 0
 udMC3: uDX[11,14] + uDW[10,2] - uDW[10,10] >= 0
 udMC3: uDX[11,14] - uDW[10,2] - uDW[10,10] <= 0
 udMC0: udx[11,3] - udw[10,15] >= 0
 udMC0: udx[11,3] - udx[11,15] >= 0
 udMC0: udx[11,3] - udx[11,15] - udw[10,15] <= 0
 udMC0: - uDX[11,3] + uDX[11,15] + uDW[10,15] >= 0
 udMC0: uDX[11,3] + uDX[11,15] - uDW[10,15] >= 0
 udMC0: uDX[11,3] - uDX[11,15] + uDW[10,15] >= 0
 udMC0: uDX[11,3] - uDX[11,15] - uDW[10,15] <= 0
 udMC1: udx[11,7] - udw[10,3] = 0
 udMC1: uDX[11,7] - uDW[10,3] = 0
 udMC2: udx[11,11] - udw[10,7] >= 0
 udMC2: udx[11,11] - udw[10,11] >= 0
 udMC2: udx[11,11] - udw[10,7] - udw[10,11] <= 0
 udMC2: - uDX[11,11] + uDW[10,7] + uDW[10,11] >= 0
 udMC2: uDX[11,11] - uDW[10,7] + uDW[10,11] >= 0
 udMC2: uDX[11,11] + uDW[10,7] - uDW[10,11] >= 0
 udMC2: uDX[11,11] - uDW[10,7] - uDW[10,11] <= 0
 udMC3: udx[11,15] - udw[10,3] >= 0
 udMC3: udx[11,15] - udw[10,11] >= 0
 udMC3: udx[11,15] - udw[10,3] - udw[10,11] <= 0
 udMC3: - uDX[11,15] + uDW[10,3] + uDW[10,11] >= 0
 udMC3: uDX[11,15] - uDW[10,3] + uDW[10,11] >= 0
 udMC3: uDX[11,15] + uDW[10,3] - uDW[10,11] >= 0
 udMC3: uDX[11,15] - uDW[10,3] - uDW[10,11] <= 0
 udMC0: udx[12,0] - udw[11,12] >= 0
 udMC0: udx[12,0] - udx[12,12] >= 0
 udMC0: udx[12,0] - udx[12,12] - udw[11,12] <= 0
 udMC0: - uDX[12,0] + uDX[12,12] + uDW[11,12] >= 0
 udMC0: uDX[12,0] + uDX[12,12] - uDW[11,12] >= 0
 udMC0: uDX[12,0] - uDX[12,12] + uDW[11,12] >= 0
 udMC0: uDX[12,0] - uDX[12,12] - uDW[11,12] <= 0
 udMC1: udx[12,4] - udw[11,0] = 0
 udMC1: uDX[12,4] - uDW[11,0] = 0
 udMC2: udx[12,8] - udw[11,4] >= 0
 udMC2: udx[12,8] - udw[11,8] >= 0
 udMC2: udx[12,8] - udw[11,4] - udw[11,8] <= 0
 udMC2: - uDX[12,8] + uDW[11,4] + uDW[11,8] >= 0
 udMC2: uDX[12,8] - uDW[11,4] + uDW[11,8] >= 0
 udMC2: uDX[12,8] + uDW[11,4] - uDW[11,8] >= 0
 udMC2: uDX[12,8] - uDW[11,4] - uDW[11,8] <= 0
 udMC3: udx[12,12] - udw[11,0] >= 0
 udMC3: udx[12,12] - udw[11,8] >= 0
 udMC3: udx[12,12] - udw[11,0] - udw[11,8] <= 0
 udMC3: - uDX[12,12] + uDW[11,0] + uDW[11,8] >= 0
 udMC3: uDX[12,12] - uDW[11,0] + uDW[11,8] >= 0
 udMC3: uDX[12,12] + uDW[11,0] - uDW[11,8] >= 0
 udMC3: uDX[12,12] - uDW[11,0] - uDW[11,8] <= 0
 udMC0: udx[12,1] - udw[11,13] >= 0
 udMC0: udx[12,1] - udx[12,13] >= 0
 udMC0: udx[12,1] - udx[12,13] - udw[11,13] <= 0
 udMC0: - uDX[12,1] + uDX[12,13] + uDW[11,13] >= 0
 udMC0: uDX[12,1] + uDX[12,13] - uDW[11,13] >= 0
 udMC0: uDX[12,1] - uDX[12,13] + uDW[11,13] >= 0
 udMC0: uDX[12,1] - uDX[12,13] - uDW[11,13] <= 0
 udMC1: udx[12,5] - udw[11,1] = 0
 udMC1: uDX[12,5] - uDW[11,1] = 0
 udMC2: udx[12,9] - udw[11,5] >= 0
 udMC2: udx[12,9] - udw[11,9] >= 0
 udMC2: udx[12,9] - udw[11,5] - udw[11,9] <= 0
 udMC2: - uDX[12,9] + uDW[11,5] + uDW[11,9] >= 0
 udMC2: uDX[12,9] - uDW[11,5] + uDW[11,9] >= 0
 udMC2: uDX[12,9] + uDW[11,5] - uDW[11,9] >= 0
 udMC2: uDX[12,9] - uDW[11,5] - uDW[11,9] <= 0
 udMC3: udx[12,13] - udw[11,1] >= 0
 udMC3: udx[12,13] - udw[11,9] >= 0
 udMC3: udx[12,13] - udw[11,1] - udw[11,9] <= 0
 udMC3: - uDX[12,13] + uDW[11,1] + uDW[11,9] >= 0
 udMC3: uDX[12,13] - uDW[11,1] + uDW[11,9] >= 0
 udMC3: uDX[12,13] + uDW[11,1] - uDW[11,9] >= 0
 udMC3: uDX[12,13] - uDW[11,1] - uDW[11,9] <= 0
 udMC0: udx[12,2] - udw[11,14] >= 0
 udMC0: udx[12,2] - udx[12,14] >= 0
 udMC0: udx[12,2] - udx[12,14] - udw[11,14] <= 0
 udMC0: - uDX[12,2] + uDX[12,14] + uDW[11,14] >= 0
 udMC0: uDX[12,2] + uDX[12,14] - uDW[11,14] >= 0
 udMC0: uDX[12,2] - uDX[12,14] + uDW[11,14] >= 0
 udMC0: uDX[12,2] - uDX[12,14] - uDW[11,14] <= 0
 udMC1: udx[12,6] - udw[11,2] = 0
 udMC1: uDX[12,6] - uDW[11,2] = 0
 udMC2: udx[12,10] - udw[11,6] >= 0
 udMC2: udx[12,10] - udw[11,10] >= 0
 udMC2: udx[12,10] - udw[11,6] - udw[11,10] <= 0
 udMC2: - uDX[12,10] + uDW[11,6] + uDW[11,10] >= 0
 udMC2: uDX[12,10] - uDW[11,6] + uDW[11,10] >= 0
 udMC2: uDX[12,10] + uDW[11,6] - uDW[11,10] >= 0
 udMC2: uDX[12,10] - uDW[11,6] - uDW[11,10] <= 0
 udMC3: udx[12,14] - udw[11,2] >= 0
 udMC3: udx[12,14] - udw[11,10] >= 0
 udMC3: udx[12,14] - udw[11,2] - udw[11,10] <= 0
 udMC3: - uDX[12,14] + uDW[11,2] + uDW[11,10] >= 0
 udMC3: uDX[12,14] - uDW[11,2] + uDW[11,10] >= 0
 udMC3: uDX[12,14] + uDW[11,2] - uDW[11,10] >= 0
 udMC3: uDX[12,14] - uDW[11,2] - uDW[11,10] <= 0
 udMC0: udx[12,3] - udw[11,15] >= 0
 udMC0: udx[12,3] - udx[12,15] >= 0
 udMC0: udx[12,3] - udx[12,15] - udw[11,15] <= 0
 udMC0: - uDX[12,3] + uDX[12,15] + uDW[11,15] >= 0
 udMC0: uDX[12,3] + uDX[12,15] - uDW[11,15] >= 0
 udMC0: uDX[12,3] - uDX[12,15] + uDW[11,15] >= 0
 udMC0: uDX[12,3] - uDX[12,15] - uDW[11,15] <= 0
 udMC1: udx[12,7] - udw[11,3] = 0
 udMC1: uDX[12,7] - uDW[11,3] = 0
 udMC2: udx[12,11] - udw[11,7] >= 0
 udMC2: udx[12,11] - udw[11,11] >= 0
 udMC2: udx[12,11] - udw[11,7] - udw[11,11] <= 0
 udMC2: - uDX[12,11] + uDW[11,7] + uDW[11,11] >= 0
 udMC2: uDX[12,11] - uDW[11,7] + uDW[11,11] >= 0
 udMC2: uDX[12,11] + uDW[11,7] - uDW[11,11] >= 0
 udMC2: uDX[12,11] - uDW[11,7] - uDW[11,11] <= 0
 udMC3: udx[12,15] - udw[11,3] >= 0
 udMC3: udx[12,15] - udw[11,11] >= 0
 udMC3: udx[12,15] - udw[11,3] - udw[11,11] <= 0
 udMC3: - uDX[12,15] + uDW[11,3] + uDW[11,11] >= 0
 udMC3: uDX[12,15] - uDW[11,3] + uDW[11,11] >= 0
 udMC3: uDX[12,15] + uDW[11,3] - uDW[11,11] >= 0
 udMC3: uDX[12,15] - uDW[11,3] - uDW[11,11] <= 0
 uKS: - ustk[0,0] + uLANE[0] >= 0
 uKS: - ustk[1,8] + uLANE[0] >= 0
 uKS: - ustk[2,2] + uLANE[0] >= 0
 uKS: - ustk[3,10] + uLANE[0] >= 0
 uKS: - ustk[4,4] + uLANE[0] >= 0
 uKS: - ustk[5,12] + uLANE[0] >= 0
 uKS: - ustk[6,6] + uLANE[0] >= 0
 uKS: - ustk[7,14] + uLANE[0] >= 0
 uKS: - ustk[8,5] + uLANE[0] >= 0
 uKS: - ustk[9,13] + uLANE[0] >= 0
 uKS: - ustk[10,3] + uLANE[0] >= 0
 uKS: - ustk[11,11] + uLANE[0] >= 0
 uKS: - ustk[12,7] + uLANE[0] >= 0
 uKS: - ustk[13,15] + uLANE[0] >= 0
 uKS: - ustk[14,1] + uLANE[0] >= 0
 uKS: - ustk[15,9] + uLANE[0] >= 0
 uKS: - ustk[16,0] + uLANE[0] >= 0
 uKS: - ustk[17,8] + uLANE[0] >= 0
 uKS: - ustk[18,2] + uLANE[0] >= 0
 uKS: - ustk[19,10] + uLANE[0] >= 0
 uKS: - ustk[20,4] + uLANE[0] >= 0
 uKS: - ustk[21,12] + uLANE[0] >= 0
 uKS: - ustk[22,6] + uLANE[0] >= 0
 uKS: - ustk[23,14] + uLANE[0] >= 0
 uKS: - ustk[24,5] + uLANE[0] >= 0
 uKS: - ustk[25,13] + uLANE[0] >= 0
 uKS: - ustk[26,3] + uLANE[0] >= 0
 uKS: - ustk[27,11] + uLANE[0] >= 0
 uKS: - ustk[28,7] + uLANE[0] >= 0
 uKS: - ustk[29,15] + uLANE[0] >= 0
 uKS: - ustk[30,1] + uLANE[0] >= 0
 uKS: - ustk[31,9] + uLANE[0] >= 0
 uKS: - ustk[32,0] + uLANE[0] >= 0
 uKS: - ustk[33,8] + uLANE[0] >= 0
 uKS: - ustk[34,2] + uLANE[0] >= 0
 uKS: - ustk[35,10] + uLANE[0] >= 0
 uKS: - ustk[36,4] + uLANE[0] >= 0
 uKS: - ustk[37,12] + uLANE[0] >= 0
 uKS: - ustk[38,6] + uLANE[0] >= 0
 R3227: - ustk[0,0] - ustk[1,8] - ustk[2,2] - ustk[3,10] - ustk[4,4]
   - ustk[5,12] - ustk[6,6] - ustk[7,14] - ustk[8,5] - ustk[9,13]
   - ustk[10,3] - ustk[11,11] - ustk[12,7] - ustk[13,15] - ustk[14,1]
   - ustk[15,9] - ustk[16,0] - ustk[17,8] - ustk[18,2] - ustk[19,10]
   - ustk[20,4] - ustk[21,12] - ustk[22,6] - ustk[23,14] - ustk[24,5]
   - ustk[25,13] - ustk[26,3] - ustk[27,11] - ustk[28,7] - ustk[29,15]
   + 30 uLANE[0] <= 1
 R3228: - ustk[30,1] - ustk[31,9] - ustk[32,0] - ustk[33,8] - ustk[34,2]
   - ustk[35,10] - ustk[36,4] - ustk[37,12] - ustk[38,6] + 9 uLANE[0] <= 1
 uKS: - ustk[0,1] + uLANE[1] >= 0
 uKS: - ustk[1,9] + uLANE[1] >= 0
 uKS: - ustk[2,0] + uLANE[1] >= 0
 uKS: - ustk[3,8] + uLANE[1] >= 0
 uKS: - ustk[4,2] + uLANE[1] >= 0
 uKS: - ustk[5,10] + uLANE[1] >= 0
 uKS: - ustk[6,4] + uLANE[1] >= 0
 uKS: - ustk[7,12] + uLANE[1] >= 0
 uKS: - ustk[8,6] + uLANE[1] >= 0
 uKS: - ustk[9,14] + uLANE[1] >= 0
 uKS: - ustk[10,5] + uLANE[1] >= 0
 uKS: - ustk[11,13] + uLANE[1] >= 0
 uKS: - ustk[12,3] + uLANE[1] >= 0
 uKS: - ustk[13,11] + uLANE[1] >= 0
 uKS: - ustk[14,7] + uLANE[1] >= 0
 uKS: - ustk[15,15] + uLANE[1] >= 0
 uKS: - ustk[16,1] + uLANE[1] >= 0
 uKS: - ustk[17,9] + uLANE[1] >= 0
 uKS: - ustk[18,0] + uLANE[1] >= 0
 uKS: - ustk[19,8] + uLANE[1] >= 0
 uKS: - ustk[20,2] + uLANE[1] >= 0
 uKS: - ustk[21,10] + uLANE[1] >= 0
 uKS: - ustk[22,4] + uLANE[1] >= 0
 uKS: - ustk[23,12] + uLANE[1] >= 0
 uKS: - ustk[24,6] + uLANE[1] >= 0
 uKS: - ustk[25,14] + uLANE[1] >= 0
 uKS: - ustk[26,5] + uLANE[1] >= 0
 uKS: - ustk[27,13] + uLANE[1] >= 0
 uKS: - ustk[28,3] + uLANE[1] >= 0
 uKS: - ustk[29,11] + uLANE[1] >= 0
 uKS: - ustk[30,7] + uLANE[1] >= 0
 uKS: - ustk[31,15] + uLANE[1] >= 0
 uKS: - ustk[32,1] + uLANE[1] >= 0
 uKS: - ustk[33,9] + uLANE[1] >= 0
 uKS: - ustk[34,0] + uLANE[1] >= 0
 uKS: - ustk[35,8] + uLANE[1] >= 0
 uKS: - ustk[36,2] + uLANE[1] >= 0
 uKS: - ustk[37,10] + uLANE[1] >= 0
 uKS: - ustk[38,4] + uLANE[1] >= 0
 R3268: - ustk[0,1] - ustk[1,9] - ustk[2,0] - ustk[3,8] - ustk[4,2]
   - ustk[5,10] - ustk[6,4] - ustk[7,12] - ustk[8,6] - ustk[9,14]
   - ustk[10,5] - ustk[11,13] - ustk[12,3] - ustk[13,11] - ustk[14,7]
   - ustk[15,15] - ustk[16,1] - ustk[17,9] - ustk[18,0] - ustk[19,8]
   - ustk[20,2] - ustk[21,10] - ustk[22,4] - ustk[23,12] - ustk[24,6]
   - ustk[25,14] - ustk[26,5] - ustk[27,13] - ustk[28,3] - ustk[29,11]
   + 30 uLANE[1] <= 1
 R3269: - ustk[30,7] - ustk[31,15] - ustk[32,1] - ustk[33,9] - ustk[34,0]
   - ustk[35,8] - ustk[36,2] - ustk[37,10] - ustk[38,4] + 9 uLANE[1] <= 1
 uKS: - ustk[0,2] + uLANE[2] >= 0
 uKS: - ustk[1,10] + uLANE[2] >= 0
 uKS: - ustk[2,4] + uLANE[2] >= 0
 uKS: - ustk[3,12] + uLANE[2] >= 0
 uKS: - ustk[4,6] + uLANE[2] >= 0
 uKS: - ustk[5,14] + uLANE[2] >= 0
 uKS: - ustk[6,5] + uLANE[2] >= 0
 uKS: - ustk[7,13] + uLANE[2] >= 0
 uKS: - ustk[8,3] + uLANE[2] >= 0
 uKS: - ustk[9,11] + uLANE[2] >= 0
 uKS: - ustk[10,7] + uLANE[2] >= 0
 uKS: - ustk[11,15] + uLANE[2] >= 0
 uKS: - ustk[12,1] + uLANE[2] >= 0
 uKS: - ustk[13,9] + uLANE[2] >= 0
 uKS: - ustk[14,0] + uLANE[2] >= 0
 uKS: - ustk[15,8] + uLANE[2] >= 0
 uKS: - ustk[16,2] + uLANE[2] >= 0
 uKS: - ustk[17,10] + uLANE[2] >= 0
 uKS: - ustk[18,4] + uLANE[2] >= 0
 uKS: - ustk[19,12] + uLANE[2] >= 0
 uKS: - ustk[20,6] + uLANE[2] >= 0
 uKS: - ustk[21,14] + uLANE[2] >= 0
 uKS: - ustk[22,5] + uLANE[2] >= 0
 uKS: - ustk[23,13] + uLANE[2] >= 0
 uKS: - ustk[24,3] + uLANE[2] >= 0
 uKS: - ustk[25,11] + uLANE[2] >= 0
 uKS: - ustk[26,7] + uLANE[2] >= 0
 uKS: - ustk[27,15] + uLANE[2] >= 0
 uKS: - ustk[28,1] + uLANE[2] >= 0
 uKS: - ustk[29,9] + uLANE[2] >= 0
 uKS: - ustk[30,0] + uLANE[2] >= 0
 uKS: - ustk[31,8] + uLANE[2] >= 0
 uKS: - ustk[32,2] + uLANE[2] >= 0
 uKS: - ustk[33,10] + uLANE[2] >= 0
 uKS: - ustk[34,4] + uLANE[2] >= 0
 uKS: - ustk[35,12] + uLANE[2] >= 0
 uKS: - ustk[36,6] + uLANE[2] >= 0
 uKS: - ustk[37,14] + uLANE[2] >= 0
 uKS: - ustk[38,5] + uLANE[2] >= 0
 R3309: - ustk[0,2] - ustk[1,10] - ustk[2,4] - ustk[3,12] - ustk[4,6]
   - ustk[5,14] - ustk[6,5] - ustk[7,13] - ustk[8,3] - ustk[9,11]
   - ustk[10,7] - ustk[11,15] - ustk[12,1] - ustk[13,9] - ustk[14,0]
   - ustk[15,8] - ustk[16,2] - ustk[17,10] - ustk[18,4] - ustk[19,12]
   - ustk[20,6] - ustk[21,14] - ustk[22,5] - ustk[23,13] - ustk[24,3]
   - ustk[25,11] - ustk[26,7] - ustk[27,15] - ustk[28,1] - ustk[29,9]
   + 30 uLANE[2] <= 1
 R3310: - ustk[30,0] - ustk[31,8] - ustk[32,2] - ustk[33,10] - ustk[34,4]
   - ustk[35,12] - ustk[36,6] - ustk[37,14] - ustk[38,5] + 9 uLANE[2] <= 1
 uKS: - ustk[0,3] + uLANE[3] >= 0
 uKS: - ustk[1,11] + uLANE[3] >= 0
 uKS: - ustk[2,7] + uLANE[3] >= 0
 uKS: - ustk[3,15] + uLANE[3] >= 0
 uKS: - ustk[4,1] + uLANE[3] >= 0
 uKS: - ustk[5,9] + uLANE[3] >= 0
 uKS: - ustk[6,0] + uLANE[3] >= 0
 uKS: - ustk[7,8] + uLANE[3] >= 0
 uKS: - ustk[8,2] + uLANE[3] >= 0
 uKS: - ustk[9,10] + uLANE[3] >= 0
 uKS: - ustk[10,4] + uLANE[3] >= 0
 uKS: - ustk[11,12] + uLANE[3] >= 0
 uKS: - ustk[12,6] + uLANE[3] >= 0
 uKS: - ustk[13,14] + uLANE[3] >= 0
 uKS: - ustk[14,5] + uLANE[3] >= 0
 uKS: - ustk[15,13] + uLANE[3] >= 0
 uKS: - ustk[16,3] + uLANE[3] >= 0
 uKS: - ustk[17,11] + uLANE[3] >= 0
 uKS: - ustk[18,7] + uLANE[3] >= 0
 uKS: - ustk[19,15] + uLANE[3] >= 0
 uKS: - ustk[20,1] + uLANE[3] >= 0
 uKS: - ustk[21,9] + uLANE[3] >= 0
 uKS: - ustk[22,0] + uLANE[3] >= 0
 uKS: - ustk[23,8] + uLANE[3] >= 0
 uKS: - ustk[24,2] + uLANE[3] >= 0
 uKS: - ustk[25,10] + uLANE[3] >= 0
 uKS: - ustk[26,4] + uLANE[3] >= 0
 uKS: - ustk[27,12] + uLANE[3] >= 0
 uKS: - ustk[28,6] + uLANE[3] >= 0
 uKS: - ustk[29,14] + uLANE[3] >= 0
 uKS: - ustk[30,5] + uLANE[3] >= 0
 uKS: - ustk[31,13] + uLANE[3] >= 0
 uKS: - ustk[32,3] + uLANE[3] >= 0
 uKS: - ustk[33,11] + uLANE[3] >= 0
 uKS: - ustk[34,7] + uLANE[3] >= 0
 uKS: - ustk[35,15] + uLANE[3] >= 0
 uKS: - ustk[36,1] + uLANE[3] >= 0
 uKS: - ustk[37,9] + uLANE[3] >= 0
 uKS: - ustk[38,0] + uLANE[3] >= 0
 R3350: - ustk[0,3] - ustk[1,11] - ustk[2,7] - ustk[3,15] - ustk[4,1]
   - ustk[5,9] - ustk[6,0] - ustk[7,8] - ustk[8,2] - ustk[9,10]
   - ustk[10,4] - ustk[11,12] - ustk[12,6] - ustk[13,14] - ustk[14,5]
   - ustk[15,13] - ustk[16,3] - ustk[17,11] - ustk[18,7] - ustk[19,15]
   - ustk[20,1] - ustk[21,9] - ustk[22,0] - ustk[23,8] - ustk[24,2]
   - ustk[25,10] - ustk[26,4] - ustk[27,12] - ustk[28,6] - ustk[29,14]
   + 30 uLANE[3] <= 1
 R3351: - ustk[30,5] - ustk[31,13] - ustk[32,3] - ustk[33,11] - ustk[34,7]
   - ustk[35,15] - ustk[36,1] - ustk[37,9] - ustk[38,0] + 9 uLANE[3] <= 1
 uKS: - ustk[0,4] + uLANE[4] >= 0
 uKS: - ustk[1,12] + uLANE[4] >= 0
 uKS: - ustk[2,6] + uLANE[4] >= 0
 uKS: - ustk[3,14] + uLANE[4] >= 0
 uKS: - ustk[4,5] + uLANE[4] >= 0
 uKS: - ustk[5,13] + uLANE[4] >= 0
 uKS: - ustk[6,3] + uLANE[4] >= 0
 uKS: - ustk[7,11] + uLANE[4] >= 0
 uKS: - ustk[8,7] + uLANE[4] >= 0
 uKS: - ustk[9,15] + uLANE[4] >= 0
 uKS: - ustk[10,1] + uLANE[4] >= 0
 uKS: - ustk[11,9] + uLANE[4] >= 0
 uKS: - ustk[12,0] + uLANE[4] >= 0
 uKS: - ustk[13,8] + uLANE[4] >= 0
 uKS: - ustk[14,2] + uLANE[4] >= 0
 uKS: - ustk[15,10] + uLANE[4] >= 0
 uKS: - ustk[16,4] + uLANE[4] >= 0
 uKS: - ustk[17,12] + uLANE[4] >= 0
 uKS: - ustk[18,6] + uLANE[4] >= 0
 uKS: - ustk[19,14] + uLANE[4] >= 0
 uKS: - ustk[20,5] + uLANE[4] >= 0
 uKS: - ustk[21,13] + uLANE[4] >= 0
 uKS: - ustk[22,3] + uLANE[4] >= 0
 uKS: - ustk[23,11] + uLANE[4] >= 0
 uKS: - ustk[24,7] + uLANE[4] >= 0
 uKS: - ustk[25,15] + uLANE[4] >= 0
 uKS: - ustk[26,1] + uLANE[4] >= 0
 uKS: - ustk[27,9] + uLANE[4] >= 0
 uKS: - ustk[28,0] + uLANE[4] >= 0
 uKS: - ustk[29,8] + uLANE[4] >= 0
 uKS: - ustk[30,2] + uLANE[4] >= 0
 uKS: - ustk[31,10] + uLANE[4] >= 0
 uKS: - ustk[32,4] + uLANE[4] >= 0
 uKS: - ustk[33,12] + uLANE[4] >= 0
 uKS: - ustk[34,6] + uLANE[4] >= 0
 uKS: - ustk[35,14] + uLANE[4] >= 0
 uKS: - ustk[36,5] + uLANE[4] >= 0
 uKS: - ustk[37,13] + uLANE[4] >= 0
 uKS: - ustk[38,3] + uLANE[4] >= 0
 R3391: - ustk[0,4] - ustk[1,12] - ustk[2,6] - ustk[3,14] - ustk[4,5]
   - ustk[5,13] - ustk[6,3] - ustk[7,11] - ustk[8,7] - ustk[9,15]
   - ustk[10,1] - ustk[11,9] - ustk[12,0] - ustk[13,8] - ustk[14,2]
   - ustk[15,10] - ustk[16,4] - ustk[17,12] - ustk[18,6] - ustk[19,14]
   - ustk[20,5] - ustk[21,13] - ustk[22,3] - ustk[23,11] - ustk[24,7]
   - ustk[25,15] - ustk[26,1] - ustk[27,9] - ustk[28,0] - ustk[29,8]
   + 30 uLANE[4] <= 1
 R3392: - ustk[30,2] - ustk[31,10] - ustk[32,4] - ustk[33,12] - ustk[34,6]
   - ustk[35,14] - ustk[36,5] - ustk[37,13] - ustk[38,3] + 9 uLANE[4] <= 1
 uKS: - ustk[0,5] + uLANE[5] >= 0
 uKS: - ustk[1,13] + uLANE[5] >= 0
 uKS: - ustk[2,3] + uLANE[5] >= 0
 uKS: - ustk[3,11] + uLANE[5] >= 0
 uKS: - ustk[4,7] + uLANE[5] >= 0
 uKS: - ustk[5,15] + uLANE[5] >= 0
 uKS: - ustk[6,1] + uLANE[5] >= 0
 uKS: - ustk[7,9] + uLANE[5] >= 0
 uKS: - ustk[8,0] + uLANE[5] >= 0
 uKS: - ustk[9,8] + uLANE[5] >= 0
 uKS: - ustk[10,2] + uLANE[5] >= 0
 uKS: - ustk[11,10] + uLANE[5] >= 0
 uKS: - ustk[12,4] + uLANE[5] >= 0
 uKS: - ustk[13,12] + uLANE[5] >= 0
 uKS: - ustk[14,6] + uLANE[5] >= 0
 uKS: - ustk[15,14] + uLANE[5] >= 0
 uKS: - ustk[16,5] + uLANE[5] >= 0
 uKS: - ustk[17,13] + uLANE[5] >= 0
 uKS: - ustk[18,3] + uLANE[5] >= 0
 uKS: - ustk[19,11] + uLANE[5] >= 0
 uKS: - ustk[20,7] + uLANE[5] >= 0
 uKS: - ustk[21,15] + uLANE[5] >= 0
 uKS: - ustk[22,1] + uLANE[5] >= 0
 uKS: - ustk[23,9] + uLANE[5] >= 0
 uKS: - ustk[24,0] + uLANE[5] >= 0
 uKS: - ustk[25,8] + uLANE[5] >= 0
 uKS: - ustk[26,2] + uLANE[5] >= 0
 uKS: - ustk[27,10] + uLANE[5] >= 0
 uKS: - ustk[28,4] + uLANE[5] >= 0
 uKS: - ustk[29,12] + uLANE[5] >= 0
 uKS: - ustk[30,6] + uLANE[5] >= 0
 uKS: - ustk[31,14] + uLANE[5] >= 0
 uKS: - ustk[32,5] + uLANE[5] >= 0
 uKS: - ustk[33,13] + uLANE[5] >= 0
 uKS: - ustk[34,3] + uLANE[5] >= 0
 uKS: - ustk[35,11] + uLANE[5] >= 0
 uKS: - ustk[36,7] + uLANE[5] >= 0
 uKS: - ustk[37,15] + uLANE[5] >= 0
 uKS: - ustk[38,1] + uLANE[5] >= 0
 R3432: - ustk[0,5] - ustk[1,13] - ustk[2,3] - ustk[3,11] - ustk[4,7]
   - ustk[5,15] - ustk[6,1] - ustk[7,9] - ustk[8,0] - ustk[9,8]
   - ustk[10,2] - ustk[11,10] - ustk[12,4] - ustk[13,12] - ustk[14,6]
   - ustk[15,14] - ustk[16,5] - ustk[17,13] - ustk[18,3] - ustk[19,11]
   - ustk[20,7] - ustk[21,15] - ustk[22,1] - ustk[23,9] - ustk[24,0]
   - ustk[25,8] - ustk[26,2] - ustk[27,10] - ustk[28,4] - ustk[29,12]
   + 30 uLANE[5] <= 1
 R3433: - ustk[30,6] - ustk[31,14] - ustk[32,5] - ustk[33,13] - ustk[34,3]
   - ustk[35,11] - ustk[36,7] - ustk[37,15] - ustk[38,1] + 9 uLANE[5] <= 1
 uKS: - ustk[0,6] + uLANE[6] >= 0
 uKS: - ustk[1,14] + uLANE[6] >= 0
 uKS: - ustk[2,5] + uLANE[6] >= 0
 uKS: - ustk[3,13] + uLANE[6] >= 0
 uKS: - ustk[4,3] + uLANE[6] >= 0
 uKS: - ustk[5,11] + uLANE[6] >= 0
 uKS: - ustk[6,7] + uLANE[6] >= 0
 uKS: - ustk[7,15] + uLANE[6] >= 0
 uKS: - ustk[8,1] + uLANE[6] >= 0
 uKS: - ustk[9,9] + uLANE[6] >= 0
 uKS: - ustk[10,0] + uLANE[6] >= 0
 uKS: - ustk[11,8] + uLANE[6] >= 0
 uKS: - ustk[12,2] + uLANE[6] >= 0
 uKS: - ustk[13,10] + uLANE[6] >= 0
 uKS: - ustk[14,4] + uLANE[6] >= 0
 uKS: - ustk[15,12] + uLANE[6] >= 0
 uKS: - ustk[16,6] + uLANE[6] >= 0
 uKS: - ustk[17,14] + uLANE[6] >= 0
 uKS: - ustk[18,5] + uLANE[6] >= 0
 uKS: - ustk[19,13] + uLANE[6] >= 0
 uKS: - ustk[20,3] + uLANE[6] >= 0
 uKS: - ustk[21,11] + uLANE[6] >= 0
 uKS: - ustk[22,7] + uLANE[6] >= 0
 uKS: - ustk[23,15] + uLANE[6] >= 0
 uKS: - ustk[24,1] + uLANE[6] >= 0
 uKS: - ustk[25,9] + uLANE[6] >= 0
 uKS: - ustk[26,0] + uLANE[6] >= 0
 uKS: - ustk[27,8] + uLANE[6] >= 0
 uKS: - ustk[28,2] + uLANE[6] >= 0
 uKS: - ustk[29,10] + uLANE[6] >= 0
 uKS: - ustk[30,4] + uLANE[6] >= 0
 uKS: - ustk[31,12] + uLANE[6] >= 0
 uKS: - ustk[32,6] + uLANE[6] >= 0
 uKS: - ustk[33,14] + uLANE[6] >= 0
 uKS: - ustk[34,5] + uLANE[6] >= 0
 uKS: - ustk[35,13] + uLANE[6] >= 0
 uKS: - ustk[36,3] + uLANE[6] >= 0
 uKS: - ustk[37,11] + uLANE[6] >= 0
 uKS: - ustk[38,7] + uLANE[6] >= 0
 R3473: - ustk[0,6] - ustk[1,14] - ustk[2,5] - ustk[3,13] - ustk[4,3]
   - ustk[5,11] - ustk[6,7] - ustk[7,15] - ustk[8,1] - ustk[9,9]
   - ustk[10,0] - ustk[11,8] - ustk[12,2] - ustk[13,10] - ustk[14,4]
   - ustk[15,12] - ustk[16,6] - ustk[17,14] - ustk[18,5] - ustk[19,13]
   - ustk[20,3] - ustk[21,11] - ustk[22,7] - ustk[23,15] - ustk[24,1]
   - ustk[25,9] - ustk[26,0] - ustk[27,8] - ustk[28,2] - ustk[29,10]
   + 30 uLANE[6] <= 1
 R3474: - ustk[30,4] - ustk[31,12] - ustk[32,6] - ustk[33,14] - ustk[34,5]
   - ustk[35,13] - ustk[36,3] - ustk[37,11] - ustk[38,7] + 9 uLANE[6] <= 1
 uKS: - ustk[0,7] + uLANE[7] >= 0
 uKS: - ustk[1,15] + uLANE[7] >= 0
 uKS: - ustk[2,1] + uLANE[7] >= 0
 uKS: - ustk[3,9] + uLANE[7] >= 0
 uKS: - ustk[4,0] + uLANE[7] >= 0
 uKS: - ustk[5,8] + uLANE[7] >= 0
 uKS: - ustk[6,2] + uLANE[7] >= 0
 uKS: - ustk[7,10] + uLANE[7] >= 0
 uKS: - ustk[8,4] + uLANE[7] >= 0
 uKS: - ustk[9,12] + uLANE[7] >= 0
 uKS: - ustk[10,6] + uLANE[7] >= 0
 uKS: - ustk[11,14] + uLANE[7] >= 0
 uKS: - ustk[12,5] + uLANE[7] >= 0
 uKS: - ustk[13,13] + uLANE[7] >= 0
 uKS: - ustk[14,3] + uLANE[7] >= 0
 uKS: - ustk[15,11] + uLANE[7] >= 0
 uKS: - ustk[16,7] + uLANE[7] >= 0
 uKS: - ustk[17,15] + uLANE[7] >= 0
 uKS: - ustk[18,1] + uLANE[7] >= 0
 uKS: - ustk[19,9] + uLANE[7] >= 0
 uKS: - ustk[20,0] + uLANE[7] >= 0
 uKS: - ustk[21,8] + uLANE[7] >= 0
 uKS: - ustk[22,2] + uLANE[7] >= 0
 uKS: - ustk[23,10] + uLANE[7] >= 0
 uKS: - ustk[24,4] + uLANE[7] >= 0
 uKS: - ustk[25,12] + uLANE[7] >= 0
 uKS: - ustk[26,6] + uLANE[7] >= 0
 uKS: - ustk[27,14] + uLANE[7] >= 0
 uKS: - ustk[28,5] + uLANE[7] >= 0
 uKS: - ustk[29,13] + uLANE[7] >= 0
 uKS: - ustk[30,3] + uLANE[7] >= 0
 uKS: - ustk[31,11] + uLANE[7] >= 0
 uKS: - ustk[32,7] + uLANE[7] >= 0
 uKS: - ustk[33,15] + uLANE[7] >= 0
 uKS: - ustk[34,1] + uLANE[7] >= 0
 uKS: - ustk[35,9] + uLANE[7] >= 0
 uKS: - ustk[36,0] + uLANE[7] >= 0
 uKS: - ustk[37,8] + uLANE[7] >= 0
 uKS: - ustk[38,2] + uLANE[7] >= 0
 R3514: - ustk[0,7] - ustk[1,15] - ustk[2,1] - ustk[3,9] - ustk[4,0]
   - ustk[5,8] - ustk[6,2] - ustk[7,10] - ustk[8,4] - ustk[9,12]
   - ustk[10,6] - ustk[11,14] - ustk[12,5] - ustk[13,13] - ustk[14,3]
   - ustk[15,11] - ustk[16,7] - ustk[17,15] - ustk[18,1] - ustk[19,9]
   - ustk[20,0] - ustk[21,8] - ustk[22,2] - ustk[23,10] - ustk[24,4]
   - ustk[25,12] - ustk[26,6] - ustk[27,14] - ustk[28,5] - ustk[29,13]
   + 30 uLANE[7] <= 1
 R3515: - ustk[30,3] - ustk[31,11] - ustk[32,7] - ustk[33,15] - ustk[34,1]
   - ustk[35,9] - ustk[36,0] - ustk[37,8] - ustk[38,2] + 9 uLANE[7] <= 1
 uKS: - ustk[0,8] + uLANE[8] >= 0
 uKS: - ustk[1,2] + uLANE[8] >= 0
 uKS: - ustk[2,10] + uLANE[8] >= 0
 uKS: - ustk[3,4] + uLANE[8] >= 0
 uKS: - ustk[4,12] + uLANE[8] >= 0
 uKS: - ustk[5,6] + uLANE[8] >= 0
 uKS: - ustk[6,14] + uLANE[8] >= 0
 uKS: - ustk[7,5] + uLANE[8] >= 0
 uKS: - ustk[8,13] + uLANE[8] >= 0
 uKS: - ustk[9,3] + uLANE[8] >= 0
 uKS: - ustk[10,11] + uLANE[8] >= 0
 uKS: - ustk[11,7] + uLANE[8] >= 0
 uKS: - ustk[12,15] + uLANE[8] >= 0
 uKS: - ustk[13,1] + uLANE[8] >= 0
 uKS: - ustk[14,9] + uLANE[8] >= 0
 uKS: - ustk[15,0] + uLANE[8] >= 0
 uKS: - ustk[16,8] + uLANE[8] >= 0
 uKS: - ustk[17,2] + uLANE[8] >= 0
 uKS: - ustk[18,10] + uLANE[8] >= 0
 uKS: - ustk[19,4] + uLANE[8] >= 0
 uKS: - ustk[20,12] + uLANE[8] >= 0
 uKS: - ustk[21,6] + uLANE[8] >= 0
 uKS: - ustk[22,14] + uLANE[8] >= 0
 uKS: - ustk[23,5] + uLANE[8] >= 0
 uKS: - ustk[24,13] + uLANE[8] >= 0
 uKS: - ustk[25,3] + uLANE[8] >= 0
 uKS: - ustk[26,11] + uLANE[8] >= 0
 uKS: - ustk[27,7] + uLANE[8] >= 0
 uKS: - ustk[28,15] + uLANE[8] >= 0
 uKS: - ustk[29,1] + uLANE[8] >= 0
 uKS: - ustk[30,9] + uLANE[8] >= 0
 uKS: - ustk[31,0] + uLANE[8] >= 0
 uKS: - ustk[32,8] + uLANE[8] >= 0
 uKS: - ustk[33,2] + uLANE[8] >= 0
 uKS: - ustk[34,10] + uLANE[8] >= 0
 uKS: - ustk[35,4] + uLANE[8] >= 0
 uKS: - ustk[36,12] + uLANE[8] >= 0
 uKS: - ustk[37,6] + uLANE[8] >= 0
 uKS: - ustk[38,14] + uLANE[8] >= 0
 R3555: - ustk[0,8] - ustk[1,2] - ustk[2,10] - ustk[3,4] - ustk[4,12]
   - ustk[5,6] - ustk[6,14] - ustk[7,5] - ustk[8,13] - ustk[9,3]
   - ustk[10,11] - ustk[11,7] - ustk[12,15] - ustk[13,1] - ustk[14,9]
   - ustk[15,0] - ustk[16,8] - ustk[17,2] - ustk[18,10] - ustk[19,4]
   - ustk[20,12] - ustk[21,6] - ustk[22,14] - ustk[23,5] - ustk[24,13]
   - ustk[25,3] - ustk[26,11] - ustk[27,7] - ustk[28,15] - ustk[29,1]
   + 30 uLANE[8] <= 1
 R3556: - ustk[30,9] - ustk[31,0] - ustk[32,8] - ustk[33,2] - ustk[34,10]
   - ustk[35,4] - ustk[36,12] - ustk[37,6] - ustk[38,14] + 9 uLANE[8] <= 1
 uKS: - ustk[0,9] + uLANE[9] >= 0
 uKS: - ustk[1,0] + uLANE[9] >= 0
 uKS: - ustk[2,8] + uLANE[9] >= 0
 uKS: - ustk[3,2] + uLANE[9] >= 0
 uKS: - ustk[4,10] + uLANE[9] >= 0
 uKS: - ustk[5,4] + uLANE[9] >= 0
 uKS: - ustk[6,12] + uLANE[9] >= 0
 uKS: - ustk[7,6] + uLANE[9] >= 0
 uKS: - ustk[8,14] + uLANE[9] >= 0
 uKS: - ustk[9,5] + uLANE[9] >= 0
 uKS: - ustk[10,13] + uLANE[9] >= 0
 uKS: - ustk[11,3] + uLANE[9] >= 0
 uKS: - ustk[12,11] + uLANE[9] >= 0
 uKS: - ustk[13,7] + uLANE[9] >= 0
 uKS: - ustk[14,15] + uLANE[9] >= 0
 uKS: - ustk[15,1] + uLANE[9] >= 0
 uKS: - ustk[16,9] + uLANE[9] >= 0
 uKS: - ustk[17,0] + uLANE[9] >= 0
 uKS: - ustk[18,8] + uLANE[9] >= 0
 uKS: - ustk[19,2] + uLANE[9] >= 0
 uKS: - ustk[20,10] + uLANE[9] >= 0
 uKS: - ustk[21,4] + uLANE[9] >= 0
 uKS: - ustk[22,12] + uLANE[9] >= 0
 uKS: - ustk[23,6] + uLANE[9] >= 0
 uKS: - ustk[24,14] + uLANE[9] >= 0
 uKS: - ustk[25,5] + uLANE[9] >= 0
 uKS: - ustk[26,13] + uLANE[9] >= 0
 uKS: - ustk[27,3] + uLANE[9] >= 0
 uKS: - ustk[28,11] + uLANE[9] >= 0
 uKS: - ustk[29,7] + uLANE[9] >= 0
 uKS: - ustk[30,15] + uLANE[9] >= 0
 uKS: - ustk[31,1] + uLANE[9] >= 0
 uKS: - ustk[32,9] + uLANE[9] >= 0
 uKS: - ustk[33,0] + uLANE[9] >= 0
 uKS: - ustk[34,8] + uLANE[9] >= 0
 uKS: - ustk[35,2] + uLANE[9] >= 0
 uKS: - ustk[36,10] + uLANE[9] >= 0
 uKS: - ustk[37,4] + uLANE[9] >= 0
 uKS: - ustk[38,12] + uLANE[9] >= 0
 R3596: - ustk[0,9] - ustk[1,0] - ustk[2,8] - ustk[3,2] - ustk[4,10]
   - ustk[5,4] - ustk[6,12] - ustk[7,6] - ustk[8,14] - ustk[9,5]
   - ustk[10,13] - ustk[11,3] - ustk[12,11] - ustk[13,7] - ustk[14,15]
   - ustk[15,1] - ustk[16,9] - ustk[17,0] - ustk[18,8] - ustk[19,2]
   - ustk[20,10] - ustk[21,4] - ustk[22,12] - ustk[23,6] - ustk[24,14]
   - ustk[25,5] - ustk[26,13] - ustk[27,3] - ustk[28,11] - ustk[29,7]
   + 30 uLANE[9] <= 1
 R3597: - ustk[30,15] - ustk[31,1] - ustk[32,9] - ustk[33,0] - ustk[34,8]
   - ustk[35,2] - ustk[36,10] - ustk[37,4] - ustk[38,12] + 9 uLANE[9] <= 1
 uKS: - ustk[0,10] + uLANE[10] >= 0
 uKS: - ustk[1,4] + uLANE[10] >= 0
 uKS: - ustk[2,12] + uLANE[10] >= 0
 uKS: - ustk[3,6] + uLANE[10] >= 0
 uKS: - ustk[4,14] + uLANE[10] >= 0
 uKS: - ustk[5,5] + uLANE[10] >= 0
 uKS: - ustk[6,13] + uLANE[10] >= 0
 uKS: - ustk[7,3] + uLANE[10] >= 0
 uKS: - ustk[8,11] + uLANE[10] >= 0
 uKS: - ustk[9,7] + uLANE[10] >= 0
 uKS: - ustk[10,15] + uLANE[10] >= 0
 uKS: - ustk[11,1] + uLANE[10] >= 0
 uKS: - ustk[12,9] + uLANE[10] >= 0
 uKS: - ustk[13,0] + uLANE[10] >= 0
 uKS: - ustk[14,8] + uLANE[10] >= 0
 uKS: - ustk[15,2] + uLANE[10] >= 0
 uKS: - ustk[16,10] + uLANE[10] >= 0
 uKS: - ustk[17,4] + uLANE[10] >= 0
 uKS: - ustk[18,12] + uLANE[10] >= 0
 uKS: - ustk[19,6] + uLANE[10] >= 0
 uKS: - ustk[20,14] + uLANE[10] >= 0
 uKS: - ustk[21,5] + uLANE[10] >= 0
 uKS: - ustk[22,13] + uLANE[10] >= 0
 uKS: - ustk[23,3] + uLANE[10] >= 0
 uKS: - ustk[24,11] + uLANE[10] >= 0
 uKS: - ustk[25,7] + uLANE[10] >= 0
 uKS: - ustk[26,15] + uLANE[10] >= 0
 uKS: - ustk[27,1] + uLANE[10] >= 0
 uKS: - ustk[28,9] + uLANE[10] >= 0
 uKS: - ustk[29,0] + uLANE[10] >= 0
 uKS: - ustk[30,8] + uLANE[10] >= 0
 uKS: - ustk[31,2] + uLANE[10] >= 0
 uKS: - ustk[32,10] + uLANE[10] >= 0
 uKS: - ustk[33,4] + uLANE[10] >= 0
 uKS: - ustk[34,12] + uLANE[10] >= 0
 uKS: - ustk[35,6] + uLANE[10] >= 0
 uKS: - ustk[36,14] + uLANE[10] >= 0
 uKS: - ustk[37,5] + uLANE[10] >= 0
 uKS: - ustk[38,13] + uLANE[10] >= 0
 R3637: - ustk[0,10] - ustk[1,4] - ustk[2,12] - ustk[3,6] - ustk[4,14]
   - ustk[5,5] - ustk[6,13] - ustk[7,3] - ustk[8,11] - ustk[9,7]
   - ustk[10,15] - ustk[11,1] - ustk[12,9] - ustk[13,0] - ustk[14,8]
   - ustk[15,2] - ustk[16,10] - ustk[17,4] - ustk[18,12] - ustk[19,6]
   - ustk[20,14] - ustk[21,5] - ustk[22,13] - ustk[23,3] - ustk[24,11]
   - ustk[25,7] - ustk[26,15] - ustk[27,1] - ustk[28,9] - ustk[29,0]
   + 30 uLANE[10] <= 1
 R3638: - ustk[30,8] - ustk[31,2] - ustk[32,10] - ustk[33,4] - ustk[34,12]
   - ustk[35,6] - ustk[36,14] - ustk[37,5] - ustk[38,13] + 9 uLANE[10]
   <= 1
 uKS: - ustk[0,11] + uLANE[11] >= 0
 uKS: - ustk[1,7] + uLANE[11] >= 0
 uKS: - ustk[2,15] + uLANE[11] >= 0
 uKS: - ustk[3,1] + uLANE[11] >= 0
 uKS: - ustk[4,9] + uLANE[11] >= 0
 uKS: - ustk[5,0] + uLANE[11] >= 0
 uKS: - ustk[6,8] + uLANE[11] >= 0
 uKS: - ustk[7,2] + uLANE[11] >= 0
 uKS: - ustk[8,10] + uLANE[11] >= 0
 uKS: - ustk[9,4] + uLANE[11] >= 0
 uKS: - ustk[10,12] + uLANE[11] >= 0
 uKS: - ustk[11,6] + uLANE[11] >= 0
 uKS: - ustk[12,14] + uLANE[11] >= 0
 uKS: - ustk[13,5] + uLANE[11] >= 0
 uKS: - ustk[14,13] + uLANE[11] >= 0
 uKS: - ustk[15,3] + uLANE[11] >= 0
 uKS: - ustk[16,11] + uLANE[11] >= 0
 uKS: - ustk[17,7] + uLANE[11] >= 0
 uKS: - ustk[18,15] + uLANE[11] >= 0
 uKS: - ustk[19,1] + uLANE[11] >= 0
 uKS: - ustk[20,9] + uLANE[11] >= 0
 uKS: - ustk[21,0] + uLANE[11] >= 0
 uKS: - ustk[22,8] + uLANE[11] >= 0
 uKS: - ustk[23,2] + uLANE[11] >= 0
 uKS: - ustk[24,10] + uLANE[11] >= 0
 uKS: - ustk[25,4] + uLANE[11] >= 0
 uKS: - ustk[26,12] + uLANE[11] >= 0
 uKS: - ustk[27,6] + uLANE[11] >= 0
 uKS: - ustk[28,14] + uLANE[11] >= 0
 uKS: - ustk[29,5] + uLANE[11] >= 0
 uKS: - ustk[30,13] + uLANE[11] >= 0
 uKS: - ustk[31,3] + uLANE[11] >= 0
 uKS: - ustk[32,11] + uLANE[11] >= 0
 uKS: - ustk[33,7] + uLANE[11] >= 0
 uKS: - ustk[34,15] + uLANE[11] >= 0
 uKS: - ustk[35,1] + uLANE[11] >= 0
 uKS: - ustk[36,9] + uLANE[11] >= 0
 uKS: - ustk[37,0] + uLANE[11] >= 0
 uKS: - ustk[38,8] + uLANE[11] >= 0
 R3678: - ustk[0,11] - ustk[1,7] - ustk[2,15] - ustk[3,1] - ustk[4,9]
   - ustk[5,0] - ustk[6,8] - ustk[7,2] - ustk[8,10] - ustk[9,4]
   - ustk[10,12] - ustk[11,6] - ustk[12,14] - ustk[13,5] - ustk[14,13]
   - ustk[15,3] - ustk[16,11] - ustk[17,7] - ustk[18,15] - ustk[19,1]
   - ustk[20,9] - ustk[21,0] - ustk[22,8] - ustk[23,2] - ustk[24,10]
   - ustk[25,4] - ustk[26,12] - ustk[27,6] - ustk[28,14] - ustk[29,5]
   + 30 uLANE[11] <= 1
 R3679: - ustk[30,13] - ustk[31,3] - ustk[32,11] - ustk[33,7] - ustk[34,15]
   - ustk[35,1] - ustk[36,9] - ustk[37,0] - ustk[38,8] + 9 uLANE[11] <= 1
 uKS: - ustk[0,12] + uLANE[12] >= 0
 uKS: - ustk[1,6] + uLANE[12] >= 0
 uKS: - ustk[2,14] + uLANE[12] >= 0
 uKS: - ustk[3,5] + uLANE[12] >= 0
 uKS: - ustk[4,13] + uLANE[12] >= 0
 uKS: - ustk[5,3] + uLANE[12] >= 0
 uKS: - ustk[6,11] + uLANE[12] >= 0
 uKS: - ustk[7,7] + uLANE[12] >= 0
 uKS: - ustk[8,15] + uLANE[12] >= 0
 uKS: - ustk[9,1] + uLANE[12] >= 0
 uKS: - ustk[10,9] + uLANE[12] >= 0
 uKS: - ustk[11,0] + uLANE[12] >= 0
 uKS: - ustk[12,8] + uLANE[12] >= 0
 uKS: - ustk[13,2] + uLANE[12] >= 0
 uKS: - ustk[14,10] + uLANE[12] >= 0
 uKS: - ustk[15,4] + uLANE[12] >= 0
 uKS: - ustk[16,12] + uLANE[12] >= 0
 uKS: - ustk[17,6] + uLANE[12] >= 0
 uKS: - ustk[18,14] + uLANE[12] >= 0
 uKS: - ustk[19,5] + uLANE[12] >= 0
 uKS: - ustk[20,13] + uLANE[12] >= 0
 uKS: - ustk[21,3] + uLANE[12] >= 0
 uKS: - ustk[22,11] + uLANE[12] >= 0
 uKS: - ustk[23,7] + uLANE[12] >= 0
 uKS: - ustk[24,15] + uLANE[12] >= 0
 uKS: - ustk[25,1] + uLANE[12] >= 0
 uKS: - ustk[26,9] + uLANE[12] >= 0
 uKS: - ustk[27,0] + uLANE[12] >= 0
 uKS: - ustk[28,8] + uLANE[12] >= 0
 uKS: - ustk[29,2] + uLANE[12] >= 0
 uKS: - ustk[30,10] + uLANE[12] >= 0
 uKS: - ustk[31,4] + uLANE[12] >= 0
 uKS: - ustk[32,12] + uLANE[12] >= 0
 uKS: - ustk[33,6] + uLANE[12] >= 0
 uKS: - ustk[34,14] + uLANE[12] >= 0
 uKS: - ustk[35,5] + uLANE[12] >= 0
 uKS: - ustk[36,13] + uLANE[12] >= 0
 uKS: - ustk[37,3] + uLANE[12] >= 0
 uKS: - ustk[38,11] + uLANE[12] >= 0
 R3719: - ustk[0,12] - ustk[1,6] - ustk[2,14] - ustk[3,5] - ustk[4,13]
   - ustk[5,3] - ustk[6,11] - ustk[7,7] - ustk[8,15] - ustk[9,1]
   - ustk[10,9] - ustk[11,0] - ustk[12,8] - ustk[13,2] - ustk[14,10]
   - ustk[15,4] - ustk[16,12] - ustk[17,6] - ustk[18,14] - ustk[19,5]
   - ustk[20,13] - ustk[21,3] - ustk[22,11] - ustk[23,7] - ustk[24,15]
   - ustk[25,1] - ustk[26,9] - ustk[27,0] - ustk[28,8] - ustk[29,2]
   + 30 uLANE[12] <= 1
 R3720: - ustk[30,10] - ustk[31,4] - ustk[32,12] - ustk[33,6] - ustk[34,14]
   - ustk[35,5] - ustk[36,13] - ustk[37,3] - ustk[38,11] + 9 uLANE[12]
   <= 1
 uKS: - ustk[0,13] + uLANE[13] >= 0
 uKS: - ustk[1,3] + uLANE[13] >= 0
 uKS: - ustk[2,11] + uLANE[13] >= 0
 uKS: - ustk[3,7] + uLANE[13] >= 0
 uKS: - ustk[4,15] + uLANE[13] >= 0
 uKS: - ustk[5,1] + uLANE[13] >= 0
 uKS: - ustk[6,9] + uLANE[13] >= 0
 uKS: - ustk[7,0] + uLANE[13] >= 0
 uKS: - ustk[8,8] + uLANE[13] >= 0
 uKS: - ustk[9,2] + uLANE[13] >= 0
 uKS: - ustk[10,10] + uLANE[13] >= 0
 uKS: - ustk[11,4] + uLANE[13] >= 0
 uKS: - ustk[12,12] + uLANE[13] >= 0
 uKS: - ustk[13,6] + uLANE[13] >= 0
 uKS: - ustk[14,14] + uLANE[13] >= 0
 uKS: - ustk[15,5] + uLANE[13] >= 0
 uKS: - ustk[16,13] + uLANE[13] >= 0
 uKS: - ustk[17,3] + uLANE[13] >= 0
 uKS: - ustk[18,11] + uLANE[13] >= 0
 uKS: - ustk[19,7] + uLANE[13] >= 0
 uKS: - ustk[20,15] + uLANE[13] >= 0
 uKS: - ustk[21,1] + uLANE[13] >= 0
 uKS: - ustk[22,9] + uLANE[13] >= 0
 uKS: - ustk[23,0] + uLANE[13] >= 0
 uKS: - ustk[24,8] + uLANE[13] >= 0
 uKS: - ustk[25,2] + uLANE[13] >= 0
 uKS: - ustk[26,10] + uLANE[13] >= 0
 uKS: - ustk[27,4] + uLANE[13] >= 0
 uKS: - ustk[28,12] + uLANE[13] >= 0
 uKS: - ustk[29,6] + uLANE[13] >= 0
 uKS: - ustk[30,14] + uLANE[13] >= 0
 uKS: - ustk[31,5] + uLANE[13] >= 0
 uKS: - ustk[32,13] + uLANE[13] >= 0
 uKS: - ustk[33,3] + uLANE[13] >= 0
 uKS: - ustk[34,11] + uLANE[13] >= 0
 uKS: - ustk[35,7] + uLANE[13] >= 0
 uKS: - ustk[36,15] + uLANE[13] >= 0
 uKS: - ustk[37,1] + uLANE[13] >= 0
 uKS: - ustk[38,9] + uLANE[13] >= 0
 R3760: - ustk[0,13] - ustk[1,3] - ustk[2,11] - ustk[3,7] - ustk[4,15]
   - ustk[5,1] - ustk[6,9] - ustk[7,0] - ustk[8,8] - ustk[9,2]
   - ustk[10,10] - ustk[11,4] - ustk[12,12] - ustk[13,6] - ustk[14,14]
   - ustk[15,5] - ustk[16,13] - ustk[17,3] - ustk[18,11] - ustk[19,7]
   - ustk[20,15] - ustk[21,1] - ustk[22,9] - ustk[23,0] - ustk[24,8]
   - ustk[25,2] - ustk[26,10] - ustk[27,4] - ustk[28,12] - ustk[29,6]
   + 30 uLANE[13] <= 1
 R3761: - ustk[30,14] - ustk[31,5] - ustk[32,13] - ustk[33,3] - ustk[34,11]
   - ustk[35,7] - ustk[36,15] - ustk[37,1] - ustk[38,9] + 9 uLANE[13] <= 1
 uKS: - ustk[0,14] + uLANE[14] >= 0
 uKS: - ustk[1,5] + uLANE[14] >= 0
 uKS: - ustk[2,13] + uLANE[14] >= 0
 uKS: - ustk[3,3] + uLANE[14] >= 0
 uKS: - ustk[4,11] + uLANE[14] >= 0
 uKS: - ustk[5,7] + uLANE[14] >= 0
 uKS: - ustk[6,15] + uLANE[14] >= 0
 uKS: - ustk[7,1] + uLANE[14] >= 0
 uKS: - ustk[8,9] + uLANE[14] >= 0
 uKS: - ustk[9,0] + uLANE[14] >= 0
 uKS: - ustk[10,8] + uLANE[14] >= 0
 uKS: - ustk[11,2] + uLANE[14] >= 0
 uKS: - ustk[12,10] + uLANE[14] >= 0
 uKS: - ustk[13,4] + uLANE[14] >= 0
 uKS: - ustk[14,12] + uLANE[14] >= 0
 uKS: - ustk[15,6] + uLANE[14] >= 0
 uKS: - ustk[16,14] + uLANE[14] >= 0
 uKS: - ustk[17,5] + uLANE[14] >= 0
 uKS: - ustk[18,13] + uLANE[14] >= 0
 uKS: - ustk[19,3] + uLANE[14] >= 0
 uKS: - ustk[20,11] + uLANE[14] >= 0
 uKS: - ustk[21,7] + uLANE[14] >= 0
 uKS: - ustk[22,15] + uLANE[14] >= 0
 uKS: - ustk[23,1] + uLANE[14] >= 0
 uKS: - ustk[24,9] + uLANE[14] >= 0
 uKS: - ustk[25,0] + uLANE[14] >= 0
 uKS: - ustk[26,8] + uLANE[14] >= 0
 uKS: - ustk[27,2] + uLANE[14] >= 0
 uKS: - ustk[28,10] + uLANE[14] >= 0
 uKS: - ustk[29,4] + uLANE[14] >= 0
 uKS: - ustk[30,12] + uLANE[14] >= 0
 uKS: - ustk[31,6] + uLANE[14] >= 0
 uKS: - ustk[32,14] + uLANE[14] >= 0
 uKS: - ustk[33,5] + uLANE[14] >= 0
 uKS: - ustk[34,13] + uLANE[14] >= 0
 uKS: - ustk[35,3] + uLANE[14] >= 0
 uKS: - ustk[36,11] + uLANE[14] >= 0
 uKS: - ustk[37,7] + uLANE[14] >= 0
 uKS: - ustk[38,15] + uLANE[14] >= 0
 R3801: - ustk[0,14] - ustk[1,5] - ustk[2,13] - ustk[3,3] - ustk[4,11]
   - ustk[5,7] - ustk[6,15] - ustk[7,1] - ustk[8,9] - ustk[9,0]
   - ustk[10,8] - ustk[11,2] - ustk[12,10] - ustk[13,4] - ustk[14,12]
   - ustk[15,6] - ustk[16,14] - ustk[17,5] - ustk[18,13] - ustk[19,3]
   - ustk[20,11] - ustk[21,7] - ustk[22,15] - ustk[23,1] - ustk[24,9]
   - ustk[25,0] - ustk[26,8] - ustk[27,2] - ustk[28,10] - ustk[29,4]
   + 30 uLANE[14] <= 1
 R3802: - ustk[30,12] - ustk[31,6] - ustk[32,14] - ustk[33,5] - ustk[34,13]
   - ustk[35,3] - ustk[36,11] - ustk[37,7] - ustk[38,15] + 9 uLANE[14]
   <= 1
 uKS: - ustk[0,15] + uLANE[15] >= 0
 uKS: - ustk[1,1] + uLANE[15] >= 0
 uKS: - ustk[2,9] + uLANE[15] >= 0
 uKS: - ustk[3,0] + uLANE[15] >= 0
 uKS: - ustk[4,8] + uLANE[15] >= 0
 uKS: - ustk[5,2] + uLANE[15] >= 0
 uKS: - ustk[6,10] + uLANE[15] >= 0
 uKS: - ustk[7,4] + uLANE[15] >= 0
 uKS: - ustk[8,12] + uLANE[15] >= 0
 uKS: - ustk[9,6] + uLANE[15] >= 0
 uKS: - ustk[10,14] + uLANE[15] >= 0
 uKS: - ustk[11,5] + uLANE[15] >= 0
 uKS: - ustk[12,13] + uLANE[15] >= 0
 uKS: - ustk[13,3] + uLANE[15] >= 0
 uKS: - ustk[14,11] + uLANE[15] >= 0
 uKS: - ustk[15,7] + uLANE[15] >= 0
 uKS: - ustk[16,15] + uLANE[15] >= 0
 uKS: - ustk[17,1] + uLANE[15] >= 0
 uKS: - ustk[18,9] + uLANE[15] >= 0
 uKS: - ustk[19,0] + uLANE[15] >= 0
 uKS: - ustk[20,8] + uLANE[15] >= 0
 uKS: - ustk[21,2] + uLANE[15] >= 0
 uKS: - ustk[22,10] + uLANE[15] >= 0
 uKS: - ustk[23,4] + uLANE[15] >= 0
 uKS: - ustk[24,12] + uLANE[15] >= 0
 uKS: - ustk[25,6] + uLANE[15] >= 0
 uKS: - ustk[26,14] + uLANE[15] >= 0
 uKS: - ustk[27,5] + uLANE[15] >= 0
 uKS: - ustk[28,13] + uLANE[15] >= 0
 uKS: - ustk[29,3] + uLANE[15] >= 0
 uKS: - ustk[30,11] + uLANE[15] >= 0
 uKS: - ustk[31,7] + uLANE[15] >= 0
 uKS: - ustk[32,15] + uLANE[15] >= 0
 uKS: - ustk[33,1] + uLANE[15] >= 0
 uKS: - ustk[34,9] + uLANE[15] >= 0
 uKS: - ustk[35,0] + uLANE[15] >= 0
 uKS: - ustk[36,8] + uLANE[15] >= 0
 uKS: - ustk[37,2] + uLANE[15] >= 0
 uKS: - ustk[38,10] + uLANE[15] >= 0
 R3842: - ustk[0,15] - ustk[1,1] - ustk[2,9] - ustk[3,0] - ustk[4,8]
   - ustk[5,2] - ustk[6,10] - ustk[7,4] - ustk[8,12] - ustk[9,6]
   - ustk[10,14] - ustk[11,5] - ustk[12,13] - ustk[13,3] - ustk[14,11]
   - ustk[15,7] - ustk[16,15] - ustk[17,1] - ustk[18,9] - ustk[19,0]
   - ustk[20,8] - ustk[21,2] - ustk[22,10] - ustk[23,4] - ustk[24,12]
   - ustk[25,6] - ustk[26,14] - ustk[27,5] - ustk[28,13] - ustk[29,3]
   + 30 uLANE[15] <= 1
 R3843: - ustk[30,11] - ustk[31,7] - ustk[32,15] - ustk[33,1] - ustk[34,9]
   - ustk[35,0] - ustk[36,8] - ustk[37,2] - ustk[38,10] + 9 uLANE[15] <= 1
 uKS: ustk[0,0] - ustk[30,1] = 0
 uKS: ustk[0,1] - ustk[30,7] = 0
 uKS: ustk[0,2] - ustk[30,0] = 0
 uKS: ustk[0,3] - ustk[30,5] = 0
 uKS: ustk[0,4] - ustk[30,2] = 0
 uKS: ustk[0,5] - ustk[30,6] = 0
 uKS: ustk[0,6] - ustk[30,4] = 0
 uKS: ustk[0,7] - ustk[30,3] = 0
 uKS: ustk[0,8] - ustk[30,9] = 0
 uKS: ustk[0,9] - ustk[30,15] = 0
 uKS: ustk[0,10] - ustk[30,8] = 0
 uKS: ustk[0,11] - ustk[30,13] = 0
 uKS: ustk[0,12] - ustk[30,10] = 0
 uKS: ustk[0,13] - ustk[30,14] = 0
 uKS: ustk[0,14] - ustk[30,12] = 0
 uKS: ustk[0,15] - ustk[30,11] = 0
 uKS: ustk[1,0] - ustk[31,1] = 0
 uKS: ustk[1,1] - ustk[31,7] = 0
 uKS: ustk[1,2] - ustk[31,0] = 0
 uKS: ustk[1,3] - ustk[31,5] = 0
 uKS: ustk[1,4] - ustk[31,2] = 0
 uKS: ustk[1,5] - ustk[31,6] = 0
 uKS: ustk[1,6] - ustk[31,4] = 0
 uKS: ustk[1,7] - ustk[31,3] = 0
 uKS: ustk[1,8] - ustk[31,9] = 0
 uKS: ustk[1,9] - ustk[31,15] = 0
 uKS: ustk[1,10] - ustk[31,8] = 0
 uKS: ustk[1,11] - ustk[31,13] = 0
 uKS: ustk[1,12] - ustk[31,10] = 0
 uKS: ustk[1,13] - ustk[31,14] = 0
 uKS: ustk[1,14] - ustk[31,12] = 0
 uKS: ustk[1,15] - ustk[31,11] = 0
 uKS: ustk[2,0] - ustk[32,1] = 0
 uKS: ustk[2,1] - ustk[32,7] = 0
 uKS: ustk[2,2] - ustk[32,0] = 0
 uKS: ustk[2,3] - ustk[32,5] = 0
 uKS: ustk[2,4] - ustk[32,2] = 0
 uKS: ustk[2,5] - ustk[32,6] = 0
 uKS: ustk[2,6] - ustk[32,4] = 0
 uKS: ustk[2,7] - ustk[32,3] = 0
 uKS: ustk[2,8] - ustk[32,9] = 0
 uKS: ustk[2,9] - ustk[32,15] = 0
 uKS: ustk[2,10] - ustk[32,8] = 0
 uKS: ustk[2,11] - ustk[32,13] = 0
 uKS: ustk[2,12] - ustk[32,10] = 0
 uKS: ustk[2,13] - ustk[32,14] = 0
 uKS: ustk[2,14] - ustk[32,12] = 0
 uKS: ustk[2,15] - ustk[32,11] = 0
 uKS: ustk[3,0] - ustk[33,1] = 0
 uKS: ustk[3,1] - ustk[33,7] = 0
 uKS: ustk[3,2] - ustk[33,0] = 0
 uKS: ustk[3,3] - ustk[33,5] = 0
 uKS: ustk[3,4] - ustk[33,2] = 0
 uKS: ustk[3,5] - ustk[33,6] = 0
 uKS: ustk[3,6] - ustk[33,4] = 0
 uKS: ustk[3,7] - ustk[33,3] = 0
 uKS: ustk[3,8] - ustk[33,9] = 0
 uKS: ustk[3,9] - ustk[33,15] = 0
 uKS: ustk[3,10] - ustk[33,8] = 0
 uKS: ustk[3,11] - ustk[33,13] = 0
 uKS: ustk[3,12] - ustk[33,10] = 0
 uKS: ustk[3,13] - ustk[33,14] = 0
 uKS: ustk[3,14] - ustk[33,12] = 0
 uKS: ustk[3,15] - ustk[33,11] = 0
 uKS: ustk[4,0] - ustk[34,1] = 0
 uKS: ustk[4,1] - ustk[34,7] = 0
 uKS: ustk[4,2] - ustk[34,0] = 0
 uKS: ustk[4,3] - ustk[34,5] = 0
 uKS: ustk[4,4] - ustk[34,2] = 0
 uKS: ustk[4,5] - ustk[34,6] = 0
 uKS: ustk[4,6] - ustk[34,4] = 0
 uKS: ustk[4,7] - ustk[34,3] = 0
 uKS: ustk[4,8] - ustk[34,9] = 0
 uKS: ustk[4,9] - ustk[34,15] = 0
 uKS: ustk[4,10] - ustk[34,8] = 0
 uKS: ustk[4,11] - ustk[34,13] = 0
 uKS: ustk[4,12] - ustk[34,10] = 0
 uKS: ustk[4,13] - ustk[34,14] = 0
 uKS: ustk[4,14] - ustk[34,12] = 0
 uKS: ustk[4,15] - ustk[34,11] = 0
 uKS: ustk[5,0] - ustk[35,1] = 0
 uKS: ustk[5,1] - ustk[35,7] = 0
 uKS: ustk[5,2] - ustk[35,0] = 0
 uKS: ustk[5,3] - ustk[35,5] = 0
 uKS: ustk[5,4] - ustk[35,2] = 0
 uKS: ustk[5,5] - ustk[35,6] = 0
 uKS: ustk[5,6] - ustk[35,4] = 0
 uKS: ustk[5,7] - ustk[35,3] = 0
 uKS: ustk[5,8] - ustk[35,9] = 0
 uKS: ustk[5,9] - ustk[35,15] = 0
 uKS: ustk[5,10] - ustk[35,8] = 0
 uKS: ustk[5,11] - ustk[35,13] = 0
 uKS: ustk[5,12] - ustk[35,10] = 0
 uKS: ustk[5,13] - ustk[35,14] = 0
 uKS: ustk[5,14] - ustk[35,12] = 0
 uKS: ustk[5,15] - ustk[35,11] = 0
 uKS: ustk[6,0] - ustk[36,1] = 0
 uKS: ustk[6,1] - ustk[36,7] = 0
 uKS: ustk[6,2] - ustk[36,0] = 0
 uKS: ustk[6,3] - ustk[36,5] = 0
 uKS: ustk[6,4] - ustk[36,2] = 0
 uKS: ustk[6,5] - ustk[36,6] = 0
 uKS: ustk[6,6] - ustk[36,4] = 0
 uKS: ustk[6,7] - ustk[36,3] = 0
 uKS: ustk[6,8] - ustk[36,9] = 0
 uKS: ustk[6,9] - ustk[36,15] = 0
 uKS: ustk[6,10] - ustk[36,8] = 0
 uKS: ustk[6,11] - ustk[36,13] = 0
 uKS: ustk[6,12] - ustk[36,10] = 0
 uKS: ustk[6,13] - ustk[36,14] = 0
 uKS: ustk[6,14] - ustk[36,12] = 0
 uKS: ustk[6,15] - ustk[36,11] = 0
 uKS: ustk[7,0] - ustk[37,1] = 0
 uKS: ustk[7,1] - ustk[37,7] = 0
 uKS: ustk[7,2] - ustk[37,0] = 0
 uKS: ustk[7,3] - ustk[37,5] = 0
 uKS: ustk[7,4] - ustk[37,2] = 0
 uKS: ustk[7,5] - ustk[37,6] = 0
 uKS: ustk[7,6] - ustk[37,4] = 0
 uKS: ustk[7,7] - ustk[37,3] = 0
 uKS: ustk[7,8] - ustk[37,9] = 0
 uKS: ustk[7,9] - ustk[37,15] = 0
 uKS: ustk[7,10] - ustk[37,8] = 0
 uKS: ustk[7,11] - ustk[37,13] = 0
 uKS: ustk[7,12] - ustk[37,10] = 0
 uKS: ustk[7,13] - ustk[37,14] = 0
 uKS: ustk[7,14] - ustk[37,12] = 0
 uKS: ustk[7,15] - ustk[37,11] = 0
 uKS: ustk[8,0] - ustk[38,1] = 0
 uKS: ustk[8,1] - ustk[38,7] = 0
 uKS: ustk[8,2] - ustk[38,0] = 0
 uKS: ustk[8,3] - ustk[38,5] = 0
 uKS: ustk[8,4] - ustk[38,2] = 0
 uKS: ustk[8,5] - ustk[38,6] = 0
 uKS: ustk[8,6] - ustk[38,4] = 0
 uKS: ustk[8,7] - ustk[38,3] = 0
 uKS: ustk[8,8] - ustk[38,9] = 0
 uKS: ustk[8,9] - ustk[38,15] = 0
 uKS: ustk[8,10] - ustk[38,8] = 0
 uKS: ustk[8,11] - ustk[38,13] = 0
 uKS: ustk[8,12] - ustk[38,10] = 0
 uKS: ustk[8,13] - ustk[38,14] = 0
 uKS: ustk[8,14] - ustk[38,12] = 0
 uKS: ustk[8,15] - ustk[38,11] = 0
 lBasic: lDX[13,0] - ldx[13,0] >= 0
 lBasic: lDY[13,0] - ldy[13,0] >= 0
 lBasic: lDZ[13,0] - ldz[13,0] >= 0
 lBasic: lDW[13,0] - ldw[13,0] >= 0
 lBasic: lDX[13,1] - ldx[13,1] >= 0
 lBasic: lDY[13,1] - ldy[13,1] >= 0
 lBasic: lDZ[13,1] - ldz[13,1] >= 0
 lBasic: lDW[13,1] - ldw[13,1] >= 0
 lBasic: lDX[13,2] - ldx[13,2] >= 0
 lBasic: lDY[13,2] - ldy[13,2] >= 0
 lBasic: lDZ[13,2] - ldz[13,2] >= 0
 lBasic: lDW[13,2] - ldw[13,2] >= 0
 lBasic: lDX[13,3] - ldx[13,3] >= 0
 lBasic: lDY[13,3] - ldy[13,3] >= 0
 lBasic: lDZ[13,3] - ldz[13,3] >= 0
 lBasic: lDW[13,3] - ldw[13,3] >= 0
 lBasic: lDX[13,4] - ldx[13,4] >= 0
 lBasic: lDY[13,4] - ldy[13,4] >= 0
 lBasic: lDZ[13,4] - ldz[13,4] >= 0
 lBasic: lDW[13,4] - ldw[13,4] >= 0
 lBasic: lDX[13,5] - ldx[13,5] >= 0
 lBasic: lDY[13,5] - ldy[13,5] >= 0
 lBasic: lDZ[13,5] - ldz[13,5] >= 0
 lBasic: lDW[13,5] - ldw[13,5] >= 0
 lBasic: lDX[13,6] - ldx[13,6] >= 0
 lBasic: lDY[13,6] - ldy[13,6] >= 0
 lBasic: lDZ[13,6] - ldz[13,6] >= 0
 lBasic: lDW[13,6] - ldw[13,6] >= 0
 lBasic: lDX[13,7] - ldx[13,7] >= 0
 lBasic: lDY[13,7] - ldy[13,7] >= 0
 lBasic: lDZ[13,7] - ldz[13,7] >= 0
 lBasic: lDW[13,7] - ldw[13,7] >= 0
 lBasic: lDX[13,8] - ldx[13,8] >= 0
 lBasic: lDY[13,8] - ldy[13,8] >= 0
 lBasic: lDZ[13,8] - ldz[13,8] >= 0
 lBasic: lDW[13,8] - ldw[13,8] >= 0
 lBasic: lDX[13,9] - ldx[13,9] >= 0
 lBasic: lDY[13,9] - ldy[13,9] >= 0
 lBasic: lDZ[13,9] - ldz[13,9] >= 0
 lBasic: lDW[13,9] - ldw[13,9] >= 0
 lBasic: lDX[13,10] - ldx[13,10] >= 0
 lBasic: lDY[13,10] - ldy[13,10] >= 0
 lBasic: lDZ[13,10] - ldz[13,10] >= 0
 lBasic: lDW[13,10] - ldw[13,10] >= 0
 lBasic: lDX[13,11] - ldx[13,11] >= 0
 lBasic: lDY[13,11] - ldy[13,11] >= 0
 lBasic: lDZ[13,11] - ldz[13,11] >= 0
 lBasic: lDW[13,11] - ldw[13,11] >= 0
 lBasic: lDX[13,12] - ldx[13,12] >= 0
 lBasic: lDY[13,12] - ldy[13,12] >= 0
 lBasic: lDZ[13,12] - ldz[13,12] >= 0
 lBasic: lDW[13,12] - ldw[13,12] >= 0
 lBasic: lDX[13,13] - ldx[13,13] >= 0
 lBasic: lDY[13,13] - ldy[13,13] >= 0
 lBasic: lDZ[13,13] - ldz[13,13] >= 0
 lBasic: lDW[13,13] - ldw[13,13] >= 0
 lBasic: lDX[13,14] - ldx[13,14] >= 0
 lBasic: lDY[13,14] - ldy[13,14] >= 0
 lBasic: lDZ[13,14] - ldz[13,14] >= 0
 lBasic: lDW[13,14] - ldw[13,14] >= 0
 lBasic: lDX[13,15] - ldx[13,15] >= 0
 lBasic: lDY[13,15] - ldy[13,15] >= 0
 lBasic: lDZ[13,15] - ldz[13,15] >= 0
 lBasic: lDW[13,15] - ldw[13,15] >= 0
 lBasic: lDX[14,0] - ldx[14,0] >= 0
 lBasic: lDY[14,0] - ldy[14,0] >= 0
 lBasic: lDZ[14,0] - ldz[14,0] >= 0
 lBasic: lDW[14,0] - ldw[14,0] >= 0
 lBasic: lDX[14,1] - ldx[14,1] >= 0
 lBasic: lDY[14,1] - ldy[14,1] >= 0
 lBasic: lDZ[14,1] - ldz[14,1] >= 0
 lBasic: lDW[14,1] - ldw[14,1] >= 0
 lBasic: lDX[14,2] - ldx[14,2] >= 0
 lBasic: lDY[14,2] - ldy[14,2] >= 0
 lBasic: lDZ[14,2] - ldz[14,2] >= 0
 lBasic: lDW[14,2] - ldw[14,2] >= 0
 lBasic: lDX[14,3] - ldx[14,3] >= 0
 lBasic: lDY[14,3] - ldy[14,3] >= 0
 lBasic: lDZ[14,3] - ldz[14,3] >= 0
 lBasic: lDW[14,3] - ldw[14,3] >= 0
 lBasic: lDX[14,4] - ldx[14,4] >= 0
 lBasic: lDY[14,4] - ldy[14,4] >= 0
 lBasic: lDZ[14,4] - ldz[14,4] >= 0
 lBasic: lDW[14,4] - ldw[14,4] >= 0
 lBasic: lDX[14,5] - ldx[14,5] >= 0
 lBasic: lDY[14,5] - ldy[14,5] >= 0
 lBasic: lDZ[14,5] - ldz[14,5] >= 0
 lBasic: lDW[14,5] - ldw[14,5] >= 0
 lBasic: lDX[14,6] - ldx[14,6] >= 0
 lBasic: lDY[14,6] - ldy[14,6] >= 0
 lBasic: lDZ[14,6] - ldz[14,6] >= 0
 lBasic: lDW[14,6] - ldw[14,6] >= 0
 lBasic: lDX[14,7] - ldx[14,7] >= 0
 lBasic: lDY[14,7] - ldy[14,7] >= 0
 lBasic: lDZ[14,7] - ldz[14,7] >= 0
 lBasic: lDW[14,7] - ldw[14,7] >= 0
 lBasic: lDX[14,8] - ldx[14,8] >= 0
 lBasic: lDY[14,8] - ldy[14,8] >= 0
 lBasic: lDZ[14,8] - ldz[14,8] >= 0
 lBasic: lDW[14,8] - ldw[14,8] >= 0
 lBasic: lDX[14,9] - ldx[14,9] >= 0
 lBasic: lDY[14,9] - ldy[14,9] >= 0
 lBasic: lDZ[14,9] - ldz[14,9] >= 0
 lBasic: lDW[14,9] - ldw[14,9] >= 0
 lBasic: lDX[14,10] - ldx[14,10] >= 0
 lBasic: lDY[14,10] - ldy[14,10] >= 0
 lBasic: lDZ[14,10] - ldz[14,10] >= 0
 lBasic: lDW[14,10] - ldw[14,10] >= 0
 lBasic: lDX[14,11] - ldx[14,11] >= 0
 lBasic: lDY[14,11] - ldy[14,11] >= 0
 lBasic: lDZ[14,11] - ldz[14,11] >= 0
 lBasic: lDW[14,11] - ldw[14,11] >= 0
 lBasic: lDX[14,12] - ldx[14,12] >= 0
 lBasic: lDY[14,12] - ldy[14,12] >= 0
 lBasic: lDZ[14,12] - ldz[14,12] >= 0
 lBasic: lDW[14,12] - ldw[14,12] >= 0
 lBasic: lDX[14,13] - ldx[14,13] >= 0
 lBasic: lDY[14,13] - ldy[14,13] >= 0
 lBasic: lDZ[14,13] - ldz[14,13] >= 0
 lBasic: lDW[14,13] - ldw[14,13] >= 0
 lBasic: lDX[14,14] - ldx[14,14] >= 0
 lBasic: lDY[14,14] - ldy[14,14] >= 0
 lBasic: lDZ[14,14] - ldz[14,14] >= 0
 lBasic: lDW[14,14] - ldw[14,14] >= 0
 lBasic: lDX[14,15] - ldx[14,15] >= 0
 lBasic: lDY[14,15] - ldy[14,15] >= 0
 lBasic: lDZ[14,15] - ldz[14,15] >= 0
 lBasic: lDW[14,15] - ldw[14,15] >= 0
 lBasic: lDX[15,0] - ldx[15,0] >= 0
 lBasic: lDY[15,0] - ldy[15,0] >= 0
 lBasic: lDZ[15,0] - ldz[15,0] >= 0
 lBasic: lDW[15,0] - ldw[15,0] >= 0
 lBasic: lDX[15,1] - ldx[15,1] >= 0
 lBasic: lDY[15,1] - ldy[15,1] >= 0
 lBasic: lDZ[15,1] - ldz[15,1] >= 0
 lBasic: lDW[15,1] - ldw[15,1] >= 0
 lBasic: lDX[15,2] - ldx[15,2] >= 0
 lBasic: lDY[15,2] - ldy[15,2] >= 0
 lBasic: lDZ[15,2] - ldz[15,2] >= 0
 lBasic: lDW[15,2] - ldw[15,2] >= 0
 lBasic: lDX[15,3] - ldx[15,3] >= 0
 lBasic: lDY[15,3] - ldy[15,3] >= 0
 lBasic: lDZ[15,3] - ldz[15,3] >= 0
 lBasic: lDW[15,3] - ldw[15,3] >= 0
 lBasic: lDX[15,4] - ldx[15,4] >= 0
 lBasic: lDY[15,4] - ldy[15,4] >= 0
 lBasic: lDZ[15,4] - ldz[15,4] >= 0
 lBasic: lDW[15,4] - ldw[15,4] >= 0
 lBasic: lDX[15,5] - ldx[15,5] >= 0
 lBasic: lDY[15,5] - ldy[15,5] >= 0
 lBasic: lDZ[15,5] - ldz[15,5] >= 0
 lBasic: lDW[15,5] - ldw[15,5] >= 0
 lBasic: lDX[15,6] - ldx[15,6] >= 0
 lBasic: lDY[15,6] - ldy[15,6] >= 0
 lBasic: lDZ[15,6] - ldz[15,6] >= 0
 lBasic: lDW[15,6] - ldw[15,6] >= 0
 lBasic: lDX[15,7] - ldx[15,7] >= 0
 lBasic: lDY[15,7] - ldy[15,7] >= 0
 lBasic: lDZ[15,7] - ldz[15,7] >= 0
 lBasic: lDW[15,7] - ldw[15,7] >= 0
 lBasic: lDX[15,8] - ldx[15,8] >= 0
 lBasic: lDY[15,8] - ldy[15,8] >= 0
 lBasic: lDZ[15,8] - ldz[15,8] >= 0
 lBasic: lDW[15,8] - ldw[15,8] >= 0
 lBasic: lDX[15,9] - ldx[15,9] >= 0
 lBasic: lDY[15,9] - ldy[15,9] >= 0
 lBasic: lDZ[15,9] - ldz[15,9] >= 0
 lBasic: lDW[15,9] - ldw[15,9] >= 0
 lBasic: lDX[15,10] - ldx[15,10] >= 0
 lBasic: lDY[15,10] - ldy[15,10] >= 0
 lBasic: lDZ[15,10] - ldz[15,10] >= 0
 lBasic: lDW[15,10] - ldw[15,10] >= 0
 lBasic: lDX[15,11] - ldx[15,11] >= 0
 lBasic: lDY[15,11] - ldy[15,11] >= 0
 lBasic: lDZ[15,11] - ldz[15,11] >= 0
 lBasic: lDW[15,11] - ldw[15,11] >= 0
 lBasic: lDX[15,12] - ldx[15,12] >= 0
 lBasic: lDY[15,12] - ldy[15,12] >= 0
 lBasic: lDZ[15,12] - ldz[15,12] >= 0
 lBasic: lDW[15,12] - ldw[15,12] >= 0
 lBasic: lDX[15,13] - ldx[15,13] >= 0
 lBasic: lDY[15,13] - ldy[15,13] >= 0
 lBasic: lDZ[15,13] - ldz[15,13] >= 0
 lBasic: lDW[15,13] - ldw[15,13] >= 0
 lBasic: lDX[15,14] - ldx[15,14] >= 0
 lBasic: lDY[15,14] - ldy[15,14] >= 0
 lBasic: lDZ[15,14] - ldz[15,14] >= 0
 lBasic: lDW[15,14] - ldw[15,14] >= 0
 lBasic: lDX[15,15] - ldx[15,15] >= 0
 lBasic: lDY[15,15] - ldy[15,15] >= 0
 lBasic: lDZ[15,15] - ldz[15,15] >= 0
 lBasic: lDW[15,15] - ldw[15,15] >= 0
 lBasic: lDX[16,0] - ldx[16,0] >= 0
 lBasic: lDY[16,0] - ldy[16,0] >= 0
 lBasic: lDZ[16,0] - ldz[16,0] >= 0
 lBasic: lDW[16,0] - ldw[16,0] >= 0
 lBasic: lDX[16,1] - ldx[16,1] >= 0
 lBasic: lDY[16,1] - ldy[16,1] >= 0
 lBasic: lDZ[16,1] - ldz[16,1] >= 0
 lBasic: lDW[16,1] - ldw[16,1] >= 0
 lBasic: lDX[16,2] - ldx[16,2] >= 0
 lBasic: lDY[16,2] - ldy[16,2] >= 0
 lBasic: lDZ[16,2] - ldz[16,2] >= 0
 lBasic: lDW[16,2] - ldw[16,2] >= 0
 lBasic: lDX[16,3] - ldx[16,3] >= 0
 lBasic: lDY[16,3] - ldy[16,3] >= 0
 lBasic: lDZ[16,3] - ldz[16,3] >= 0
 lBasic: lDW[16,3] - ldw[16,3] >= 0
 lBasic: lDX[16,4] - ldx[16,4] >= 0
 lBasic: lDY[16,4] - ldy[16,4] >= 0
 lBasic: lDZ[16,4] - ldz[16,4] >= 0
 lBasic: lDW[16,4] - ldw[16,4] >= 0
 lBasic: lDX[16,5] - ldx[16,5] >= 0
 lBasic: lDY[16,5] - ldy[16,5] >= 0
 lBasic: lDZ[16,5] - ldz[16,5] >= 0
 lBasic: lDW[16,5] - ldw[16,5] >= 0
 lBasic: lDX[16,6] - ldx[16,6] >= 0
 lBasic: lDY[16,6] - ldy[16,6] >= 0
 lBasic: lDZ[16,6] - ldz[16,6] >= 0
 lBasic: lDW[16,6] - ldw[16,6] >= 0
 lBasic: lDX[16,7] - ldx[16,7] >= 0
 lBasic: lDY[16,7] - ldy[16,7] >= 0
 lBasic: lDZ[16,7] - ldz[16,7] >= 0
 lBasic: lDW[16,7] - ldw[16,7] >= 0
 lBasic: lDX[16,8] - ldx[16,8] >= 0
 lBasic: lDY[16,8] - ldy[16,8] >= 0
 lBasic: lDZ[16,8] - ldz[16,8] >= 0
 lBasic: lDW[16,8] - ldw[16,8] >= 0
 lBasic: lDX[16,9] - ldx[16,9] >= 0
 lBasic: lDY[16,9] - ldy[16,9] >= 0
 lBasic: lDZ[16,9] - ldz[16,9] >= 0
 lBasic: lDW[16,9] - ldw[16,9] >= 0
 lBasic: lDX[16,10] - ldx[16,10] >= 0
 lBasic: lDY[16,10] - ldy[16,10] >= 0
 lBasic: lDZ[16,10] - ldz[16,10] >= 0
 lBasic: lDW[16,10] - ldw[16,10] >= 0
 lBasic: lDX[16,11] - ldx[16,11] >= 0
 lBasic: lDY[16,11] - ldy[16,11] >= 0
 lBasic: lDZ[16,11] - ldz[16,11] >= 0
 lBasic: lDW[16,11] - ldw[16,11] >= 0
 lBasic: lDX[16,12] - ldx[16,12] >= 0
 lBasic: lDY[16,12] - ldy[16,12] >= 0
 lBasic: lDZ[16,12] - ldz[16,12] >= 0
 lBasic: lDW[16,12] - ldw[16,12] >= 0
 lBasic: lDX[16,13] - ldx[16,13] >= 0
 lBasic: lDY[16,13] - ldy[16,13] >= 0
 lBasic: lDZ[16,13] - ldz[16,13] >= 0
 lBasic: lDW[16,13] - ldw[16,13] >= 0
 lBasic: lDX[16,14] - ldx[16,14] >= 0
 lBasic: lDY[16,14] - ldy[16,14] >= 0
 lBasic: lDZ[16,14] - ldz[16,14] >= 0
 lBasic: lDW[16,14] - ldw[16,14] >= 0
 lBasic: lDX[16,15] - ldx[16,15] >= 0
 lBasic: lDY[16,15] - ldy[16,15] >= 0
 lBasic: lDZ[16,15] - ldz[16,15] >= 0
 lBasic: lDW[16,15] - ldw[16,15] >= 0
 lBasic: lDX[17,0] - ldx[17,0] >= 0
 lBasic: lDY[17,0] - ldy[17,0] >= 0
 lBasic: lDZ[17,0] - ldz[17,0] >= 0
 lBasic: lDW[17,0] - ldw[17,0] >= 0
 lBasic: lDX[17,1] - ldx[17,1] >= 0
 lBasic: lDY[17,1] - ldy[17,1] >= 0
 lBasic: lDZ[17,1] - ldz[17,1] >= 0
 lBasic: lDW[17,1] - ldw[17,1] >= 0
 lBasic: lDX[17,2] - ldx[17,2] >= 0
 lBasic: lDY[17,2] - ldy[17,2] >= 0
 lBasic: lDZ[17,2] - ldz[17,2] >= 0
 lBasic: lDW[17,2] - ldw[17,2] >= 0
 lBasic: lDX[17,3] - ldx[17,3] >= 0
 lBasic: lDY[17,3] - ldy[17,3] >= 0
 lBasic: lDZ[17,3] - ldz[17,3] >= 0
 lBasic: lDW[17,3] - ldw[17,3] >= 0
 lBasic: lDX[17,4] - ldx[17,4] >= 0
 lBasic: lDY[17,4] - ldy[17,4] >= 0
 lBasic: lDZ[17,4] - ldz[17,4] >= 0
 lBasic: lDW[17,4] - ldw[17,4] >= 0
 lBasic: lDX[17,5] - ldx[17,5] >= 0
 lBasic: lDY[17,5] - ldy[17,5] >= 0
 lBasic: lDZ[17,5] - ldz[17,5] >= 0
 lBasic: lDW[17,5] - ldw[17,5] >= 0
 lBasic: lDX[17,6] - ldx[17,6] >= 0
 lBasic: lDY[17,6] - ldy[17,6] >= 0
 lBasic: lDZ[17,6] - ldz[17,6] >= 0
 lBasic: lDW[17,6] - ldw[17,6] >= 0
 lBasic: lDX[17,7] - ldx[17,7] >= 0
 lBasic: lDY[17,7] - ldy[17,7] >= 0
 lBasic: lDZ[17,7] - ldz[17,7] >= 0
 lBasic: lDW[17,7] - ldw[17,7] >= 0
 lBasic: lDX[17,8] - ldx[17,8] >= 0
 lBasic: lDY[17,8] - ldy[17,8] >= 0
 lBasic: lDZ[17,8] - ldz[17,8] >= 0
 lBasic: lDW[17,8] - ldw[17,8] >= 0
 lBasic: lDX[17,9] - ldx[17,9] >= 0
 lBasic: lDY[17,9] - ldy[17,9] >= 0
 lBasic: lDZ[17,9] - ldz[17,9] >= 0
 lBasic: lDW[17,9] - ldw[17,9] >= 0
 lBasic: lDX[17,10] - ldx[17,10] >= 0
 lBasic: lDY[17,10] - ldy[17,10] >= 0
 lBasic: lDZ[17,10] - ldz[17,10] >= 0
 lBasic: lDW[17,10] - ldw[17,10] >= 0
 lBasic: lDX[17,11] - ldx[17,11] >= 0
 lBasic: lDY[17,11] - ldy[17,11] >= 0
 lBasic: lDZ[17,11] - ldz[17,11] >= 0
 lBasic: lDW[17,11] - ldw[17,11] >= 0
 lBasic: lDX[17,12] - ldx[17,12] >= 0
 lBasic: lDY[17,12] - ldy[17,12] >= 0
 lBasic: lDZ[17,12] - ldz[17,12] >= 0
 lBasic: lDW[17,12] - ldw[17,12] >= 0
 lBasic: lDX[17,13] - ldx[17,13] >= 0
 lBasic: lDY[17,13] - ldy[17,13] >= 0
 lBasic: lDZ[17,13] - ldz[17,13] >= 0
 lBasic: lDW[17,13] - ldw[17,13] >= 0
 lBasic: lDX[17,14] - ldx[17,14] >= 0
 lBasic: lDY[17,14] - ldy[17,14] >= 0
 lBasic: lDZ[17,14] - ldz[17,14] >= 0
 lBasic: lDW[17,14] - ldw[17,14] >= 0
 lBasic: lDX[17,15] - ldx[17,15] >= 0
 lBasic: lDY[17,15] - ldy[17,15] >= 0
 lBasic: lDZ[17,15] - ldz[17,15] >= 0
 lBasic: lDW[17,15] - ldw[17,15] >= 0
 lBasic: lDX[18,0] - ldx[18,0] >= 0
 lBasic: lDY[18,0] - ldy[18,0] >= 0
 lBasic: lDZ[18,0] - ldz[18,0] >= 0
 lBasic: lDW[18,0] - ldw[18,0] >= 0
 lBasic: lDX[18,1] - ldx[18,1] >= 0
 lBasic: lDY[18,1] - ldy[18,1] >= 0
 lBasic: lDZ[18,1] - ldz[18,1] >= 0
 lBasic: lDW[18,1] - ldw[18,1] >= 0
 lBasic: lDX[18,2] - ldx[18,2] >= 0
 lBasic: lDY[18,2] - ldy[18,2] >= 0
 lBasic: lDZ[18,2] - ldz[18,2] >= 0
 lBasic: lDW[18,2] - ldw[18,2] >= 0
 lBasic: lDX[18,3] - ldx[18,3] >= 0
 lBasic: lDY[18,3] - ldy[18,3] >= 0
 lBasic: lDZ[18,3] - ldz[18,3] >= 0
 lBasic: lDW[18,3] - ldw[18,3] >= 0
 lBasic: lDX[18,4] - ldx[18,4] >= 0
 lBasic: lDY[18,4] - ldy[18,4] >= 0
 lBasic: lDZ[18,4] - ldz[18,4] >= 0
 lBasic: lDW[18,4] - ldw[18,4] >= 0
 lBasic: lDX[18,5] - ldx[18,5] >= 0
 lBasic: lDY[18,5] - ldy[18,5] >= 0
 lBasic: lDZ[18,5] - ldz[18,5] >= 0
 lBasic: lDW[18,5] - ldw[18,5] >= 0
 lBasic: lDX[18,6] - ldx[18,6] >= 0
 lBasic: lDY[18,6] - ldy[18,6] >= 0
 lBasic: lDZ[18,6] - ldz[18,6] >= 0
 lBasic: lDW[18,6] - ldw[18,6] >= 0
 lBasic: lDX[18,7] - ldx[18,7] >= 0
 lBasic: lDY[18,7] - ldy[18,7] >= 0
 lBasic: lDZ[18,7] - ldz[18,7] >= 0
 lBasic: lDW[18,7] - ldw[18,7] >= 0
 lBasic: lDX[18,8] - ldx[18,8] >= 0
 lBasic: lDY[18,8] - ldy[18,8] >= 0
 lBasic: lDZ[18,8] - ldz[18,8] >= 0
 lBasic: lDW[18,8] - ldw[18,8] >= 0
 lBasic: lDX[18,9] - ldx[18,9] >= 0
 lBasic: lDY[18,9] - ldy[18,9] >= 0
 lBasic: lDZ[18,9] - ldz[18,9] >= 0
 lBasic: lDW[18,9] - ldw[18,9] >= 0
 lBasic: lDX[18,10] - ldx[18,10] >= 0
 lBasic: lDY[18,10] - ldy[18,10] >= 0
 lBasic: lDZ[18,10] - ldz[18,10] >= 0
 lBasic: lDW[18,10] - ldw[18,10] >= 0
 lBasic: lDX[18,11] - ldx[18,11] >= 0
 lBasic: lDY[18,11] - ldy[18,11] >= 0
 lBasic: lDZ[18,11] - ldz[18,11] >= 0
 lBasic: lDW[18,11] - ldw[18,11] >= 0
 lBasic: lDX[18,12] - ldx[18,12] >= 0
 lBasic: lDY[18,12] - ldy[18,12] >= 0
 lBasic: lDZ[18,12] - ldz[18,12] >= 0
 lBasic: lDW[18,12] - ldw[18,12] >= 0
 lBasic: lDX[18,13] - ldx[18,13] >= 0
 lBasic: lDY[18,13] - ldy[18,13] >= 0
 lBasic: lDZ[18,13] - ldz[18,13] >= 0
 lBasic: lDW[18,13] - ldw[18,13] >= 0
 lBasic: lDX[18,14] - ldx[18,14] >= 0
 lBasic: lDY[18,14] - ldy[18,14] >= 0
 lBasic: lDZ[18,14] - ldz[18,14] >= 0
 lBasic: lDW[18,14] - ldw[18,14] >= 0
 lBasic: lDX[18,15] - ldx[18,15] >= 0
 lBasic: lDY[18,15] - ldy[18,15] >= 0
 lBasic: lDZ[18,15] - ldz[18,15] >= 0
 lBasic: lDW[18,15] - ldw[18,15] >= 0
 lBasic: lDX[19,0] - ldx[19,0] >= 0
 lBasic: lDY[19,0] - ldy[19,0] >= 0
 lBasic: lDZ[19,0] - ldz[19,0] >= 0
 lBasic: lDW[19,0] - ldw[19,0] >= 0
 lBasic: lDX[19,1] - ldx[19,1] >= 0
 lBasic: lDY[19,1] - ldy[19,1] >= 0
 lBasic: lDZ[19,1] - ldz[19,1] >= 0
 lBasic: lDW[19,1] - ldw[19,1] >= 0
 lBasic: lDX[19,2] - ldx[19,2] >= 0
 lBasic: lDY[19,2] - ldy[19,2] >= 0
 lBasic: lDZ[19,2] - ldz[19,2] >= 0
 lBasic: lDW[19,2] - ldw[19,2] >= 0
 lBasic: lDX[19,3] - ldx[19,3] >= 0
 lBasic: lDY[19,3] - ldy[19,3] >= 0
 lBasic: lDZ[19,3] - ldz[19,3] >= 0
 lBasic: lDW[19,3] - ldw[19,3] >= 0
 lBasic: lDX[19,4] - ldx[19,4] >= 0
 lBasic: lDY[19,4] - ldy[19,4] >= 0
 lBasic: lDZ[19,4] - ldz[19,4] >= 0
 lBasic: lDW[19,4] - ldw[19,4] >= 0
 lBasic: lDX[19,5] - ldx[19,5] >= 0
 lBasic: lDY[19,5] - ldy[19,5] >= 0
 lBasic: lDZ[19,5] - ldz[19,5] >= 0
 lBasic: lDW[19,5] - ldw[19,5] >= 0
 lBasic: lDX[19,6] - ldx[19,6] >= 0
 lBasic: lDY[19,6] - ldy[19,6] >= 0
 lBasic: lDZ[19,6] - ldz[19,6] >= 0
 lBasic: lDW[19,6] - ldw[19,6] >= 0
 lBasic: lDX[19,7] - ldx[19,7] >= 0
 lBasic: lDY[19,7] - ldy[19,7] >= 0
 lBasic: lDZ[19,7] - ldz[19,7] >= 0
 lBasic: lDW[19,7] - ldw[19,7] >= 0
 lBasic: lDX[19,8] - ldx[19,8] >= 0
 lBasic: lDY[19,8] - ldy[19,8] >= 0
 lBasic: lDZ[19,8] - ldz[19,8] >= 0
 lBasic: lDW[19,8] - ldw[19,8] >= 0
 lBasic: lDX[19,9] - ldx[19,9] >= 0
 lBasic: lDY[19,9] - ldy[19,9] >= 0
 lBasic: lDZ[19,9] - ldz[19,9] >= 0
 lBasic: lDW[19,9] - ldw[19,9] >= 0
 lBasic: lDX[19,10] - ldx[19,10] >= 0
 lBasic: lDY[19,10] - ldy[19,10] >= 0
 lBasic: lDZ[19,10] - ldz[19,10] >= 0
 lBasic: lDW[19,10] - ldw[19,10] >= 0
 lBasic: lDX[19,11] - ldx[19,11] >= 0
 lBasic: lDY[19,11] - ldy[19,11] >= 0
 lBasic: lDZ[19,11] - ldz[19,11] >= 0
 lBasic: lDW[19,11] - ldw[19,11] >= 0
 lBasic: lDX[19,12] - ldx[19,12] >= 0
 lBasic: lDY[19,12] - ldy[19,12] >= 0
 lBasic: lDZ[19,12] - ldz[19,12] >= 0
 lBasic: lDW[19,12] - ldw[19,12] >= 0
 lBasic: lDX[19,13] - ldx[19,13] >= 0
 lBasic: lDY[19,13] - ldy[19,13] >= 0
 lBasic: lDZ[19,13] - ldz[19,13] >= 0
 lBasic: lDW[19,13] - ldw[19,13] >= 0
 lBasic: lDX[19,14] - ldx[19,14] >= 0
 lBasic: lDY[19,14] - ldy[19,14] >= 0
 lBasic: lDZ[19,14] - ldz[19,14] >= 0
 lBasic: lDW[19,14] - ldw[19,14] >= 0
 lBasic: lDX[19,15] - ldx[19,15] >= 0
 lBasic: lDY[19,15] - ldy[19,15] >= 0
 lBasic: lDZ[19,15] - ldz[19,15] >= 0
 lBasic: lDW[19,15] - ldw[19,15] >= 0
 lBasic: lDX[20,0] - ldx[20,0] >= 0
 lBasic: lDY[20,0] - ldy[20,0] >= 0
 lBasic: lDZ[20,0] - ldz[20,0] >= 0
 lBasic: lDW[20,0] - ldw[20,0] >= 0
 lBasic: lDX[20,1] - ldx[20,1] >= 0
 lBasic: lDY[20,1] - ldy[20,1] >= 0
 lBasic: lDZ[20,1] - ldz[20,1] >= 0
 lBasic: lDW[20,1] - ldw[20,1] >= 0
 lBasic: lDX[20,2] - ldx[20,2] >= 0
 lBasic: lDY[20,2] - ldy[20,2] >= 0
 lBasic: lDZ[20,2] - ldz[20,2] >= 0
 lBasic: lDW[20,2] - ldw[20,2] >= 0
 lBasic: lDX[20,3] - ldx[20,3] >= 0
 lBasic: lDY[20,3] - ldy[20,3] >= 0
 lBasic: lDZ[20,3] - ldz[20,3] >= 0
 lBasic: lDW[20,3] - ldw[20,3] >= 0
 lBasic: lDX[20,4] - ldx[20,4] >= 0
 lBasic: lDY[20,4] - ldy[20,4] >= 0
 lBasic: lDZ[20,4] - ldz[20,4] >= 0
 lBasic: lDW[20,4] - ldw[20,4] >= 0
 lBasic: lDX[20,5] - ldx[20,5] >= 0
 lBasic: lDY[20,5] - ldy[20,5] >= 0
 lBasic: lDZ[20,5] - ldz[20,5] >= 0
 lBasic: lDW[20,5] - ldw[20,5] >= 0
 lBasic: lDX[20,6] - ldx[20,6] >= 0
 lBasic: lDY[20,6] - ldy[20,6] >= 0
 lBasic: lDZ[20,6] - ldz[20,6] >= 0
 lBasic: lDW[20,6] - ldw[20,6] >= 0
 lBasic: lDX[20,7] - ldx[20,7] >= 0
 lBasic: lDY[20,7] - ldy[20,7] >= 0
 lBasic: lDZ[20,7] - ldz[20,7] >= 0
 lBasic: lDW[20,7] - ldw[20,7] >= 0
 lBasic: lDX[20,8] - ldx[20,8] >= 0
 lBasic: lDY[20,8] - ldy[20,8] >= 0
 lBasic: lDZ[20,8] - ldz[20,8] >= 0
 lBasic: lDW[20,8] - ldw[20,8] >= 0
 lBasic: lDX[20,9] - ldx[20,9] >= 0
 lBasic: lDY[20,9] - ldy[20,9] >= 0
 lBasic: lDZ[20,9] - ldz[20,9] >= 0
 lBasic: lDW[20,9] - ldw[20,9] >= 0
 lBasic: lDX[20,10] - ldx[20,10] >= 0
 lBasic: lDY[20,10] - ldy[20,10] >= 0
 lBasic: lDZ[20,10] - ldz[20,10] >= 0
 lBasic: lDW[20,10] - ldw[20,10] >= 0
 lBasic: lDX[20,11] - ldx[20,11] >= 0
 lBasic: lDY[20,11] - ldy[20,11] >= 0
 lBasic: lDZ[20,11] - ldz[20,11] >= 0
 lBasic: lDW[20,11] - ldw[20,11] >= 0
 lBasic: lDX[20,12] - ldx[20,12] >= 0
 lBasic: lDY[20,12] - ldy[20,12] >= 0
 lBasic: lDZ[20,12] - ldz[20,12] >= 0
 lBasic: lDW[20,12] - ldw[20,12] >= 0
 lBasic: lDX[20,13] - ldx[20,13] >= 0
 lBasic: lDY[20,13] - ldy[20,13] >= 0
 lBasic: lDZ[20,13] - ldz[20,13] >= 0
 lBasic: lDW[20,13] - ldw[20,13] >= 0
 lBasic: lDX[20,14] - ldx[20,14] >= 0
 lBasic: lDY[20,14] - ldy[20,14] >= 0
 lBasic: lDZ[20,14] - ldz[20,14] >= 0
 lBasic: lDW[20,14] - ldw[20,14] >= 0
 lBasic: lDX[20,15] - ldx[20,15] >= 0
 lBasic: lDY[20,15] - ldy[20,15] >= 0
 lBasic: lDZ[20,15] - ldz[20,15] >= 0
 lBasic: lDW[20,15] - ldw[20,15] >= 0
 lBasic: lDX[21,0] - ldx[21,0] >= 0
 lBasic: lDY[21,0] - ldy[21,0] >= 0
 lBasic: lDZ[21,0] - ldz[21,0] >= 0
 lBasic: lDW[21,0] - ldw[21,0] >= 0
 lBasic: lDX[21,1] - ldx[21,1] >= 0
 lBasic: lDY[21,1] - ldy[21,1] >= 0
 lBasic: lDZ[21,1] - ldz[21,1] >= 0
 lBasic: lDW[21,1] - ldw[21,1] >= 0
 lBasic: lDX[21,2] - ldx[21,2] >= 0
 lBasic: lDY[21,2] - ldy[21,2] >= 0
 lBasic: lDZ[21,2] - ldz[21,2] >= 0
 lBasic: lDW[21,2] - ldw[21,2] >= 0
 lBasic: lDX[21,3] - ldx[21,3] >= 0
 lBasic: lDY[21,3] - ldy[21,3] >= 0
 lBasic: lDZ[21,3] - ldz[21,3] >= 0
 lBasic: lDW[21,3] - ldw[21,3] >= 0
 lBasic: lDX[21,4] - ldx[21,4] >= 0
 lBasic: lDY[21,4] - ldy[21,4] >= 0
 lBasic: lDZ[21,4] - ldz[21,4] >= 0
 lBasic: lDW[21,4] - ldw[21,4] >= 0
 lBasic: lDX[21,5] - ldx[21,5] >= 0
 lBasic: lDY[21,5] - ldy[21,5] >= 0
 lBasic: lDZ[21,5] - ldz[21,5] >= 0
 lBasic: lDW[21,5] - ldw[21,5] >= 0
 lBasic: lDX[21,6] - ldx[21,6] >= 0
 lBasic: lDY[21,6] - ldy[21,6] >= 0
 lBasic: lDZ[21,6] - ldz[21,6] >= 0
 lBasic: lDW[21,6] - ldw[21,6] >= 0
 lBasic: lDX[21,7] - ldx[21,7] >= 0
 lBasic: lDY[21,7] - ldy[21,7] >= 0
 lBasic: lDZ[21,7] - ldz[21,7] >= 0
 lBasic: lDW[21,7] - ldw[21,7] >= 0
 lBasic: lDX[21,8] - ldx[21,8] >= 0
 lBasic: lDY[21,8] - ldy[21,8] >= 0
 lBasic: lDZ[21,8] - ldz[21,8] >= 0
 lBasic: lDW[21,8] - ldw[21,8] >= 0
 lBasic: lDX[21,9] - ldx[21,9] >= 0
 lBasic: lDY[21,9] - ldy[21,9] >= 0
 lBasic: lDZ[21,9] - ldz[21,9] >= 0
 lBasic: lDW[21,9] - ldw[21,9] >= 0
 lBasic: lDX[21,10] - ldx[21,10] >= 0
 lBasic: lDY[21,10] - ldy[21,10] >= 0
 lBasic: lDZ[21,10] - ldz[21,10] >= 0
 lBasic: lDW[21,10] - ldw[21,10] >= 0
 lBasic: lDX[21,11] - ldx[21,11] >= 0
 lBasic: lDY[21,11] - ldy[21,11] >= 0
 lBasic: lDZ[21,11] - ldz[21,11] >= 0
 lBasic: lDW[21,11] - ldw[21,11] >= 0
 lBasic: lDX[21,12] - ldx[21,12] >= 0
 lBasic: lDY[21,12] - ldy[21,12] >= 0
 lBasic: lDZ[21,12] - ldz[21,12] >= 0
 lBasic: lDW[21,12] - ldw[21,12] >= 0
 lBasic: lDX[21,13] - ldx[21,13] >= 0
 lBasic: lDY[21,13] - ldy[21,13] >= 0
 lBasic: lDZ[21,13] - ldz[21,13] >= 0
 lBasic: lDW[21,13] - ldw[21,13] >= 0
 lBasic: lDX[21,14] - ldx[21,14] >= 0
 lBasic: lDY[21,14] - ldy[21,14] >= 0
 lBasic: lDZ[21,14] - ldz[21,14] >= 0
 lBasic: lDW[21,14] - ldw[21,14] >= 0
 lBasic: lDX[21,15] - ldx[21,15] >= 0
 lBasic: lDY[21,15] - ldy[21,15] >= 0
 lBasic: lDZ[21,15] - ldz[21,15] >= 0
 lBasic: lDW[21,15] - ldw[21,15] >= 0
 lBasic: lDX[22,0] - ldx[22,0] >= 0
 lBasic: lDY[22,0] - ldy[22,0] >= 0
 lBasic: lDZ[22,0] - ldz[22,0] >= 0
 lBasic: lDW[22,0] - ldw[22,0] >= 0
 lBasic: lDX[22,1] - ldx[22,1] >= 0
 lBasic: lDY[22,1] - ldy[22,1] >= 0
 lBasic: lDZ[22,1] - ldz[22,1] >= 0
 lBasic: lDW[22,1] - ldw[22,1] >= 0
 lBasic: lDX[22,2] - ldx[22,2] >= 0
 lBasic: lDY[22,2] - ldy[22,2] >= 0
 lBasic: lDZ[22,2] - ldz[22,2] >= 0
 lBasic: lDW[22,2] - ldw[22,2] >= 0
 lBasic: lDX[22,3] - ldx[22,3] >= 0
 lBasic: lDY[22,3] - ldy[22,3] >= 0
 lBasic: lDZ[22,3] - ldz[22,3] >= 0
 lBasic: lDW[22,3] - ldw[22,3] >= 0
 lBasic: lDX[22,4] - ldx[22,4] >= 0
 lBasic: lDY[22,4] - ldy[22,4] >= 0
 lBasic: lDZ[22,4] - ldz[22,4] >= 0
 lBasic: lDW[22,4] - ldw[22,4] >= 0
 lBasic: lDX[22,5] - ldx[22,5] >= 0
 lBasic: lDY[22,5] - ldy[22,5] >= 0
 lBasic: lDZ[22,5] - ldz[22,5] >= 0
 lBasic: lDW[22,5] - ldw[22,5] >= 0
 lBasic: lDX[22,6] - ldx[22,6] >= 0
 lBasic: lDY[22,6] - ldy[22,6] >= 0
 lBasic: lDZ[22,6] - ldz[22,6] >= 0
 lBasic: lDW[22,6] - ldw[22,6] >= 0
 lBasic: lDX[22,7] - ldx[22,7] >= 0
 lBasic: lDY[22,7] - ldy[22,7] >= 0
 lBasic: lDZ[22,7] - ldz[22,7] >= 0
 lBasic: lDW[22,7] - ldw[22,7] >= 0
 lBasic: lDX[22,8] - ldx[22,8] >= 0
 lBasic: lDY[22,8] - ldy[22,8] >= 0
 lBasic: lDZ[22,8] - ldz[22,8] >= 0
 lBasic: lDW[22,8] - ldw[22,8] >= 0
 lBasic: lDX[22,9] - ldx[22,9] >= 0
 lBasic: lDY[22,9] - ldy[22,9] >= 0
 lBasic: lDZ[22,9] - ldz[22,9] >= 0
 lBasic: lDW[22,9] - ldw[22,9] >= 0
 lBasic: lDX[22,10] - ldx[22,10] >= 0
 lBasic: lDY[22,10] - ldy[22,10] >= 0
 lBasic: lDZ[22,10] - ldz[22,10] >= 0
 lBasic: lDW[22,10] - ldw[22,10] >= 0
 lBasic: lDX[22,11] - ldx[22,11] >= 0
 lBasic: lDY[22,11] - ldy[22,11] >= 0
 lBasic: lDZ[22,11] - ldz[22,11] >= 0
 lBasic: lDW[22,11] - ldw[22,11] >= 0
 lBasic: lDX[22,12] - ldx[22,12] >= 0
 lBasic: lDY[22,12] - ldy[22,12] >= 0
 lBasic: lDZ[22,12] - ldz[22,12] >= 0
 lBasic: lDW[22,12] - ldw[22,12] >= 0
 lBasic: lDX[22,13] - ldx[22,13] >= 0
 lBasic: lDY[22,13] - ldy[22,13] >= 0
 lBasic: lDZ[22,13] - ldz[22,13] >= 0
 lBasic: lDW[22,13] - ldw[22,13] >= 0
 lBasic: lDX[22,14] - ldx[22,14] >= 0
 lBasic: lDY[22,14] - ldy[22,14] >= 0
 lBasic: lDZ[22,14] - ldz[22,14] >= 0
 lBasic: lDW[22,14] - ldw[22,14] >= 0
 lBasic: lDX[22,15] - ldx[22,15] >= 0
 lBasic: lDY[22,15] - ldy[22,15] >= 0
 lBasic: lDZ[22,15] - ldz[22,15] >= 0
 lBasic: lDW[22,15] - ldw[22,15] >= 0
 lBasic: lDX[23,0] - ldx[23,0] >= 0
 lBasic: lDY[23,0] - ldy[23,0] >= 0
 lBasic: lDZ[23,0] - ldz[23,0] >= 0
 lBasic: lDW[23,0] - ldw[23,0] >= 0
 lBasic: lDX[23,1] - ldx[23,1] >= 0
 lBasic: lDY[23,1] - ldy[23,1] >= 0
 lBasic: lDZ[23,1] - ldz[23,1] >= 0
 lBasic: lDW[23,1] - ldw[23,1] >= 0
 lBasic: lDX[23,2] - ldx[23,2] >= 0
 lBasic: lDY[23,2] - ldy[23,2] >= 0
 lBasic: lDZ[23,2] - ldz[23,2] >= 0
 lBasic: lDW[23,2] - ldw[23,2] >= 0
 lBasic: lDX[23,3] - ldx[23,3] >= 0
 lBasic: lDY[23,3] - ldy[23,3] >= 0
 lBasic: lDZ[23,3] - ldz[23,3] >= 0
 lBasic: lDW[23,3] - ldw[23,3] >= 0
 lBasic: lDX[23,4] - ldx[23,4] >= 0
 lBasic: lDY[23,4] - ldy[23,4] >= 0
 lBasic: lDZ[23,4] - ldz[23,4] >= 0
 lBasic: lDW[23,4] - ldw[23,4] >= 0
 lBasic: lDX[23,5] - ldx[23,5] >= 0
 lBasic: lDY[23,5] - ldy[23,5] >= 0
 lBasic: lDZ[23,5] - ldz[23,5] >= 0
 lBasic: lDW[23,5] - ldw[23,5] >= 0
 lBasic: lDX[23,6] - ldx[23,6] >= 0
 lBasic: lDY[23,6] - ldy[23,6] >= 0
 lBasic: lDZ[23,6] - ldz[23,6] >= 0
 lBasic: lDW[23,6] - ldw[23,6] >= 0
 lBasic: lDX[23,7] - ldx[23,7] >= 0
 lBasic: lDY[23,7] - ldy[23,7] >= 0
 lBasic: lDZ[23,7] - ldz[23,7] >= 0
 lBasic: lDW[23,7] - ldw[23,7] >= 0
 lBasic: lDX[23,8] - ldx[23,8] >= 0
 lBasic: lDY[23,8] - ldy[23,8] >= 0
 lBasic: lDZ[23,8] - ldz[23,8] >= 0
 lBasic: lDW[23,8] - ldw[23,8] >= 0
 lBasic: lDX[23,9] - ldx[23,9] >= 0
 lBasic: lDY[23,9] - ldy[23,9] >= 0
 lBasic: lDZ[23,9] - ldz[23,9] >= 0
 lBasic: lDW[23,9] - ldw[23,9] >= 0
 lBasic: lDX[23,10] - ldx[23,10] >= 0
 lBasic: lDY[23,10] - ldy[23,10] >= 0
 lBasic: lDZ[23,10] - ldz[23,10] >= 0
 lBasic: lDW[23,10] - ldw[23,10] >= 0
 lBasic: lDX[23,11] - ldx[23,11] >= 0
 lBasic: lDY[23,11] - ldy[23,11] >= 0
 lBasic: lDZ[23,11] - ldz[23,11] >= 0
 lBasic: lDW[23,11] - ldw[23,11] >= 0
 lBasic: lDX[23,12] - ldx[23,12] >= 0
 lBasic: lDY[23,12] - ldy[23,12] >= 0
 lBasic: lDZ[23,12] - ldz[23,12] >= 0
 lBasic: lDW[23,12] - ldw[23,12] >= 0
 lBasic: lDX[23,13] - ldx[23,13] >= 0
 lBasic: lDY[23,13] - ldy[23,13] >= 0
 lBasic: lDZ[23,13] - ldz[23,13] >= 0
 lBasic: lDW[23,13] - ldw[23,13] >= 0
 lBasic: lDX[23,14] - ldx[23,14] >= 0
 lBasic: lDY[23,14] - ldy[23,14] >= 0
 lBasic: lDZ[23,14] - ldz[23,14] >= 0
 lBasic: lDW[23,14] - ldw[23,14] >= 0
 lBasic: lDX[23,15] - ldx[23,15] >= 0
 lBasic: lDY[23,15] - ldy[23,15] >= 0
 lBasic: lDZ[23,15] - ldz[23,15] >= 0
 lBasic: lDW[23,15] - ldw[23,15] >= 0
 lBasic: lDY[12,0] - ldy[12,0] >= 0
 lBasic: lDZ[12,0] - ldz[12,0] >= 0
 lBasic: lDW[12,0] - ldw[12,0] >= 0
 lBasic: lDY[12,1] - ldy[12,1] >= 0
 lBasic: lDZ[12,1] - ldz[12,1] >= 0
 lBasic: lDW[12,1] - ldw[12,1] >= 0
 lBasic: lDY[12,2] - ldy[12,2] >= 0
 lBasic: lDZ[12,2] - ldz[12,2] >= 0
 lBasic: lDW[12,2] - ldw[12,2] >= 0
 lBasic: lDY[12,3] - ldy[12,3] >= 0
 lBasic: lDZ[12,3] - ldz[12,3] >= 0
 lBasic: lDW[12,3] - ldw[12,3] >= 0
 lBasic: lDY[12,4] - ldy[12,4] >= 0
 lBasic: lDZ[12,4] - ldz[12,4] >= 0
 lBasic: lDW[12,4] - ldw[12,4] >= 0
 lBasic: lDY[12,5] - ldy[12,5] >= 0
 lBasic: lDZ[12,5] - ldz[12,5] >= 0
 lBasic: lDW[12,5] - ldw[12,5] >= 0
 lBasic: lDY[12,6] - ldy[12,6] >= 0
 lBasic: lDZ[12,6] - ldz[12,6] >= 0
 lBasic: lDW[12,6] - ldw[12,6] >= 0
 lBasic: lDY[12,7] - ldy[12,7] >= 0
 lBasic: lDZ[12,7] - ldz[12,7] >= 0
 lBasic: lDW[12,7] - ldw[12,7] >= 0
 lBasic: lDY[12,8] - ldy[12,8] >= 0
 lBasic: lDZ[12,8] - ldz[12,8] >= 0
 lBasic: lDW[12,8] - ldw[12,8] >= 0
 lBasic: lDY[12,9] - ldy[12,9] >= 0
 lBasic: lDZ[12,9] - ldz[12,9] >= 0
 lBasic: lDW[12,9] - ldw[12,9] >= 0
 lBasic: lDY[12,10] - ldy[12,10] >= 0
 lBasic: lDZ[12,10] - ldz[12,10] >= 0
 lBasic: lDW[12,10] - ldw[12,10] >= 0
 lBasic: lDY[12,11] - ldy[12,11] >= 0
 lBasic: lDZ[12,11] - ldz[12,11] >= 0
 lBasic: lDW[12,11] - ldw[12,11] >= 0
 lBasic: lDY[12,12] - ldy[12,12] >= 0
 lBasic: lDZ[12,12] - ldz[12,12] >= 0
 lBasic: lDW[12,12] - ldw[12,12] >= 0
 lBasic: lDY[12,13] - ldy[12,13] >= 0
 lBasic: lDZ[12,13] - ldz[12,13] >= 0
 lBasic: lDW[12,13] - ldw[12,13] >= 0
 lBasic: lDY[12,14] - ldy[12,14] >= 0
 lBasic: lDZ[12,14] - ldz[12,14] >= 0
 lBasic: lDW[12,14] - ldw[12,14] >= 0
 lBasic: lDY[12,15] - ldy[12,15] >= 0
 lBasic: lDZ[12,15] - ldz[12,15] >= 0
 lBasic: lDW[12,15] - ldw[12,15] >= 0
 lBasic: lDX[24,0] - ldx[24,0] >= 0
 lBasic: lDX[24,1] - ldx[24,1] >= 0
 lBasic: lDX[24,2] - ldx[24,2] >= 0
 lBasic: lDX[24,3] - ldx[24,3] >= 0
 lBasic: lDX[24,4] - ldx[24,4] >= 0
 lBasic: lDX[24,5] - ldx[24,5] >= 0
 lBasic: lDX[24,6] - ldx[24,6] >= 0
 lBasic: lDX[24,7] - ldx[24,7] >= 0
 lBasic: lDX[24,8] - ldx[24,8] >= 0
 lBasic: lDX[24,9] - ldx[24,9] >= 0
 lBasic: lDX[24,10] - ldx[24,10] >= 0
 lBasic: lDX[24,11] - ldx[24,11] >= 0
 lBasic: lDX[24,12] - ldx[24,12] >= 0
 lBasic: lDX[24,13] - ldx[24,13] >= 0
 lBasic: lDX[24,14] - ldx[24,14] >= 0
 lBasic: lDX[24,15] - ldx[24,15] >= 0
 ldSC: lDX[13,0] - lDY[13,0] = 0
 ldSC: ldx[13,0] - lDY[13,0] = 0
 ldSC: lDX[13,1] - lDY[13,1] = 0
 ldSC: ldx[13,1] - lDY[13,1] = 0
 ldSC: lDX[13,2] - lDY[13,2] = 0
 ldSC: ldx[13,2] - lDY[13,2] = 0
 ldSC: lDX[13,3] - lDY[13,3] = 0
 ldSC: ldx[13,3] - lDY[13,3] = 0
 ldSC: lDX[13,4] - lDY[13,4] = 0
 ldSC: ldx[13,4] - lDY[13,4] = 0
 ldSC: lDX[13,5] - lDY[13,5] = 0
 ldSC: ldx[13,5] - lDY[13,5] = 0
 ldSC: lDX[13,6] - lDY[13,6] = 0
 ldSC: ldx[13,6] - lDY[13,6] = 0
 ldSC: lDX[13,7] - lDY[13,7] = 0
 ldSC: ldx[13,7] - lDY[13,7] = 0
 ldSC: lDX[13,8] - lDY[13,8] = 0
 ldSC: ldx[13,8] - lDY[13,8] = 0
 ldSC: lDX[13,9] - lDY[13,9] = 0
 ldSC: ldx[13,9] - lDY[13,9] = 0
 ldSC: lDX[13,10] - lDY[13,10] = 0
 ldSC: ldx[13,10] - lDY[13,10] = 0
 ldSC: lDX[13,11] - lDY[13,11] = 0
 ldSC: ldx[13,11] - lDY[13,11] = 0
 ldSC: lDX[13,12] - lDY[13,12] = 0
 ldSC: ldx[13,12] - lDY[13,12] = 0
 ldSC: lDX[13,13] - lDY[13,13] = 0
 ldSC: ldx[13,13] - lDY[13,13] = 0
 ldSC: lDX[13,14] - lDY[13,14] = 0
 ldSC: ldx[13,14] - lDY[13,14] = 0
 ldSC: lDX[13,15] - lDY[13,15] = 0
 ldSC: ldx[13,15] - lDY[13,15] = 0
 ldSC: lDX[14,0] - lDY[14,0] = 0
 ldSC: ldx[14,0] - lDY[14,0] = 0
 ldSC: lDX[14,1] - lDY[14,1] = 0
 ldSC: ldx[14,1] - lDY[14,1] = 0
 ldSC: lDX[14,2] - lDY[14,2] = 0
 ldSC: ldx[14,2] - lDY[14,2] = 0
 ldSC: lDX[14,3] - lDY[14,3] = 0
 ldSC: ldx[14,3] - lDY[14,3] = 0
 ldSC: lDX[14,4] - lDY[14,4] = 0
 ldSC: ldx[14,4] - lDY[14,4] = 0
 ldSC: lDX[14,5] - lDY[14,5] = 0
 ldSC: ldx[14,5] - lDY[14,5] = 0
 ldSC: lDX[14,6] - lDY[14,6] = 0
 ldSC: ldx[14,6] - lDY[14,6] = 0
 ldSC: lDX[14,7] - lDY[14,7] = 0
 ldSC: ldx[14,7] - lDY[14,7] = 0
 ldSC: lDX[14,8] - lDY[14,8] = 0
 ldSC: ldx[14,8] - lDY[14,8] = 0
 ldSC: lDX[14,9] - lDY[14,9] = 0
 ldSC: ldx[14,9] - lDY[14,9] = 0
 ldSC: lDX[14,10] - lDY[14,10] = 0
 ldSC: ldx[14,10] - lDY[14,10] = 0
 ldSC: lDX[14,11] - lDY[14,11] = 0
 ldSC: ldx[14,11] - lDY[14,11] = 0
 ldSC: lDX[14,12] - lDY[14,12] = 0
 ldSC: ldx[14,12] - lDY[14,12] = 0
 ldSC: lDX[14,13] - lDY[14,13] = 0
 ldSC: ldx[14,13] - lDY[14,13] = 0
 ldSC: lDX[14,14] - lDY[14,14] = 0
 ldSC: ldx[14,14] - lDY[14,14] = 0
 ldSC: lDX[14,15] - lDY[14,15] = 0
 ldSC: ldx[14,15] - lDY[14,15] = 0
 ldSC: lDX[15,0] - lDY[15,0] = 0
 ldSC: ldx[15,0] - lDY[15,0] = 0
 ldSC: lDX[15,1] - lDY[15,1] = 0
 ldSC: ldx[15,1] - lDY[15,1] = 0
 ldSC: lDX[15,2] - lDY[15,2] = 0
 ldSC: ldx[15,2] - lDY[15,2] = 0
 ldSC: lDX[15,3] - lDY[15,3] = 0
 ldSC: ldx[15,3] - lDY[15,3] = 0
 ldSC: lDX[15,4] - lDY[15,4] = 0
 ldSC: ldx[15,4] - lDY[15,4] = 0
 ldSC: lDX[15,5] - lDY[15,5] = 0
 ldSC: ldx[15,5] - lDY[15,5] = 0
 ldSC: lDX[15,6] - lDY[15,6] = 0
 ldSC: ldx[15,6] - lDY[15,6] = 0
 ldSC: lDX[15,7] - lDY[15,7] = 0
 ldSC: ldx[15,7] - lDY[15,7] = 0
 ldSC: lDX[15,8] - lDY[15,8] = 0
 ldSC: ldx[15,8] - lDY[15,8] = 0
 ldSC: lDX[15,9] - lDY[15,9] = 0
 ldSC: ldx[15,9] - lDY[15,9] = 0
 ldSC: lDX[15,10] - lDY[15,10] = 0
 ldSC: ldx[15,10] - lDY[15,10] = 0
 ldSC: lDX[15,11] - lDY[15,11] = 0
 ldSC: ldx[15,11] - lDY[15,11] = 0
 ldSC: lDX[15,12] - lDY[15,12] = 0
 ldSC: ldx[15,12] - lDY[15,12] = 0
 ldSC: lDX[15,13] - lDY[15,13] = 0
 ldSC: ldx[15,13] - lDY[15,13] = 0
 ldSC: lDX[15,14] - lDY[15,14] = 0
 ldSC: ldx[15,14] - lDY[15,14] = 0
 ldSC: lDX[15,15] - lDY[15,15] = 0
 ldSC: ldx[15,15] - lDY[15,15] = 0
 ldSC: lDX[16,0] - lDY[16,0] = 0
 ldSC: ldx[16,0] - lDY[16,0] = 0
 ldSC: lDX[16,1] - lDY[16,1] = 0
 ldSC: ldx[16,1] - lDY[16,1] = 0
 ldSC: lDX[16,2] - lDY[16,2] = 0
 ldSC: ldx[16,2] - lDY[16,2] = 0
 ldSC: lDX[16,3] - lDY[16,3] = 0
 ldSC: ldx[16,3] - lDY[16,3] = 0
 ldSC: lDX[16,4] - lDY[16,4] = 0
 ldSC: ldx[16,4] - lDY[16,4] = 0
 ldSC: lDX[16,5] - lDY[16,5] = 0
 ldSC: ldx[16,5] - lDY[16,5] = 0
 ldSC: lDX[16,6] - lDY[16,6] = 0
 ldSC: ldx[16,6] - lDY[16,6] = 0
 ldSC: lDX[16,7] - lDY[16,7] = 0
 ldSC: ldx[16,7] - lDY[16,7] = 0
 ldSC: lDX[16,8] - lDY[16,8] = 0
 ldSC: ldx[16,8] - lDY[16,8] = 0
 ldSC: lDX[16,9] - lDY[16,9] = 0
 ldSC: ldx[16,9] - lDY[16,9] = 0
 ldSC: lDX[16,10] - lDY[16,10] = 0
 ldSC: ldx[16,10] - lDY[16,10] = 0
 ldSC: lDX[16,11] - lDY[16,11] = 0
 ldSC: ldx[16,11] - lDY[16,11] = 0
 ldSC: lDX[16,12] - lDY[16,12] = 0
 ldSC: ldx[16,12] - lDY[16,12] = 0
 ldSC: lDX[16,13] - lDY[16,13] = 0
 ldSC: ldx[16,13] - lDY[16,13] = 0
 ldSC: lDX[16,14] - lDY[16,14] = 0
 ldSC: ldx[16,14] - lDY[16,14] = 0
 ldSC: lDX[16,15] - lDY[16,15] = 0
 ldSC: ldx[16,15] - lDY[16,15] = 0
 ldSC: lDX[17,0] - lDY[17,0] = 0
 ldSC: ldx[17,0] - lDY[17,0] = 0
 ldSC: lDX[17,1] - lDY[17,1] = 0
 ldSC: ldx[17,1] - lDY[17,1] = 0
 ldSC: lDX[17,2] - lDY[17,2] = 0
 ldSC: ldx[17,2] - lDY[17,2] = 0
 ldSC: lDX[17,3] - lDY[17,3] = 0
 ldSC: ldx[17,3] - lDY[17,3] = 0
 ldSC: lDX[17,4] - lDY[17,4] = 0
 ldSC: ldx[17,4] - lDY[17,4] = 0
 ldSC: lDX[17,5] - lDY[17,5] = 0
 ldSC: ldx[17,5] - lDY[17,5] = 0
 ldSC: lDX[17,6] - lDY[17,6] = 0
 ldSC: ldx[17,6] - lDY[17,6] = 0
 ldSC: lDX[17,7] - lDY[17,7] = 0
 ldSC: ldx[17,7] - lDY[17,7] = 0
 ldSC: lDX[17,8] - lDY[17,8] = 0
 ldSC: ldx[17,8] - lDY[17,8] = 0
 ldSC: lDX[17,9] - lDY[17,9] = 0
 ldSC: ldx[17,9] - lDY[17,9] = 0
 ldSC: lDX[17,10] - lDY[17,10] = 0
 ldSC: ldx[17,10] - lDY[17,10] = 0
 ldSC: lDX[17,11] - lDY[17,11] = 0
 ldSC: ldx[17,11] - lDY[17,11] = 0
 ldSC: lDX[17,12] - lDY[17,12] = 0
 ldSC: ldx[17,12] - lDY[17,12] = 0
 ldSC: lDX[17,13] - lDY[17,13] = 0
 ldSC: ldx[17,13] - lDY[17,13] = 0
 ldSC: lDX[17,14] - lDY[17,14] = 0
 ldSC: ldx[17,14] - lDY[17,14] = 0
 ldSC: lDX[17,15] - lDY[17,15] = 0
 ldSC: ldx[17,15] - lDY[17,15] = 0
 ldSC: lDX[18,0] - lDY[18,0] = 0
 ldSC: ldx[18,0] - lDY[18,0] = 0
 ldSC: lDX[18,1] - lDY[18,1] = 0
 ldSC: ldx[18,1] - lDY[18,1] = 0
 ldSC: lDX[18,2] - lDY[18,2] = 0
 ldSC: ldx[18,2] - lDY[18,2] = 0
 ldSC: lDX[18,3] - lDY[18,3] = 0
 ldSC: ldx[18,3] - lDY[18,3] = 0
 ldSC: lDX[18,4] - lDY[18,4] = 0
 ldSC: ldx[18,4] - lDY[18,4] = 0
 ldSC: lDX[18,5] - lDY[18,5] = 0
 ldSC: ldx[18,5] - lDY[18,5] = 0
 ldSC: lDX[18,6] - lDY[18,6] = 0
 ldSC: ldx[18,6] - lDY[18,6] = 0
 ldSC: lDX[18,7] - lDY[18,7] = 0
 ldSC: ldx[18,7] - lDY[18,7] = 0
 ldSC: lDX[18,8] - lDY[18,8] = 0
 ldSC: ldx[18,8] - lDY[18,8] = 0
 ldSC: lDX[18,9] - lDY[18,9] = 0
 ldSC: ldx[18,9] - lDY[18,9] = 0
 ldSC: lDX[18,10] - lDY[18,10] = 0
 ldSC: ldx[18,10] - lDY[18,10] = 0
 ldSC: lDX[18,11] - lDY[18,11] = 0
 ldSC: ldx[18,11] - lDY[18,11] = 0
 ldSC: lDX[18,12] - lDY[18,12] = 0
 ldSC: ldx[18,12] - lDY[18,12] = 0
 ldSC: lDX[18,13] - lDY[18,13] = 0
 ldSC: ldx[18,13] - lDY[18,13] = 0
 ldSC: lDX[18,14] - lDY[18,14] = 0
 ldSC: ldx[18,14] - lDY[18,14] = 0
 ldSC: lDX[18,15] - lDY[18,15] = 0
 ldSC: ldx[18,15] - lDY[18,15] = 0
 leSC: - lDX[19,0] + lDY[19,0] = 0
 leSC: - lDX[19,0] + ldy[19,0] = 0
 leSC: - lDX[19,1] + lDY[19,1] = 0
 leSC: - lDX[19,1] + ldy[19,1] = 0
 leSC: - lDX[19,2] + lDY[19,2] = 0
 leSC: - lDX[19,2] + ldy[19,2] = 0
 leSC: - lDX[19,3] + lDY[19,3] = 0
 leSC: - lDX[19,3] + ldy[19,3] = 0
 leSC: - lDX[19,4] + lDY[19,4] = 0
 leSC: - lDX[19,4] + ldy[19,4] = 0
 leSC: - lDX[19,5] + lDY[19,5] = 0
 leSC: - lDX[19,5] + ldy[19,5] = 0
 leSC: - lDX[19,6] + lDY[19,6] = 0
 leSC: - lDX[19,6] + ldy[19,6] = 0
 leSC: - lDX[19,7] + lDY[19,7] = 0
 leSC: - lDX[19,7] + ldy[19,7] = 0
 leSC: - lDX[19,8] + lDY[19,8] = 0
 leSC: - lDX[19,8] + ldy[19,8] = 0
 leSC: - lDX[19,9] + lDY[19,9] = 0
 leSC: - lDX[19,9] + ldy[19,9] = 0
 leSC: - lDX[19,10] + lDY[19,10] = 0
 leSC: - lDX[19,10] + ldy[19,10] = 0
 leSC: - lDX[19,11] + lDY[19,11] = 0
 leSC: - lDX[19,11] + ldy[19,11] = 0
 leSC: - lDX[19,12] + lDY[19,12] = 0
 leSC: - lDX[19,12] + ldy[19,12] = 0
 leSC: - lDX[19,13] + lDY[19,13] = 0
 leSC: - lDX[19,13] + ldy[19,13] = 0
 leSC: - lDX[19,14] + lDY[19,14] = 0
 leSC: - lDX[19,14] + ldy[19,14] = 0
 leSC: - lDX[19,15] + lDY[19,15] = 0
 leSC: - lDX[19,15] + ldy[19,15] = 0
 leSC: - lDX[20,0] + lDY[20,0] = 0
 leSC: - lDX[20,0] + ldy[20,0] = 0
 leSC: - lDX[20,1] + lDY[20,1] = 0
 leSC: - lDX[20,1] + ldy[20,1] = 0
 leSC: - lDX[20,2] + lDY[20,2] = 0
 leSC: - lDX[20,2] + ldy[20,2] = 0
 leSC: - lDX[20,3] + lDY[20,3] = 0
 leSC: - lDX[20,3] + ldy[20,3] = 0
 leSC: - lDX[20,4] + lDY[20,4] = 0
 leSC: - lDX[20,4] + ldy[20,4] = 0
 leSC: - lDX[20,5] + lDY[20,5] = 0
 leSC: - lDX[20,5] + ldy[20,5] = 0
 leSC: - lDX[20,6] + lDY[20,6] = 0
 leSC: - lDX[20,6] + ldy[20,6] = 0
 leSC: - lDX[20,7] + lDY[20,7] = 0
 leSC: - lDX[20,7] + ldy[20,7] = 0
 leSC: - lDX[20,8] + lDY[20,8] = 0
 leSC: - lDX[20,8] + ldy[20,8] = 0
 leSC: - lDX[20,9] + lDY[20,9] = 0
 leSC: - lDX[20,9] + ldy[20,9] = 0
 leSC: - lDX[20,10] + lDY[20,10] = 0
 leSC: - lDX[20,10] + ldy[20,10] = 0
 leSC: - lDX[20,11] + lDY[20,11] = 0
 leSC: - lDX[20,11] + ldy[20,11] = 0
 leSC: - lDX[20,12] + lDY[20,12] = 0
 leSC: - lDX[20,12] + ldy[20,12] = 0
 leSC: - lDX[20,13] + lDY[20,13] = 0
 leSC: - lDX[20,13] + ldy[20,13] = 0
 leSC: - lDX[20,14] + lDY[20,14] = 0
 leSC: - lDX[20,14] + ldy[20,14] = 0
 leSC: - lDX[20,15] + lDY[20,15] = 0
 leSC: - lDX[20,15] + ldy[20,15] = 0
 leSC: - lDX[21,0] + lDY[21,0] = 0
 leSC: - lDX[21,0] + ldy[21,0] = 0
 leSC: - lDX[21,1] + lDY[21,1] = 0
 leSC: - lDX[21,1] + ldy[21,1] = 0
 leSC: - lDX[21,2] + lDY[21,2] = 0
 leSC: - lDX[21,2] + ldy[21,2] = 0
 leSC: - lDX[21,3] + lDY[21,3] = 0
 leSC: - lDX[21,3] + ldy[21,3] = 0
 leSC: - lDX[21,4] + lDY[21,4] = 0
 leSC: - lDX[21,4] + ldy[21,4] = 0
 leSC: - lDX[21,5] + lDY[21,5] = 0
 leSC: - lDX[21,5] + ldy[21,5] = 0
 leSC: - lDX[21,6] + lDY[21,6] = 0
 leSC: - lDX[21,6] + ldy[21,6] = 0
 leSC: - lDX[21,7] + lDY[21,7] = 0
 leSC: - lDX[21,7] + ldy[21,7] = 0
 leSC: - lDX[21,8] + lDY[21,8] = 0
 leSC: - lDX[21,8] + ldy[21,8] = 0
 leSC: - lDX[21,9] + lDY[21,9] = 0
 leSC: - lDX[21,9] + ldy[21,9] = 0
 leSC: - lDX[21,10] + lDY[21,10] = 0
 leSC: - lDX[21,10] + ldy[21,10] = 0
 leSC: - lDX[21,11] + lDY[21,11] = 0
 leSC: - lDX[21,11] + ldy[21,11] = 0
 leSC: - lDX[21,12] + lDY[21,12] = 0
 leSC: - lDX[21,12] + ldy[21,12] = 0
 leSC: - lDX[21,13] + lDY[21,13] = 0
 leSC: - lDX[21,13] + ldy[21,13] = 0
 leSC: - lDX[21,14] + lDY[21,14] = 0
 leSC: - lDX[21,14] + ldy[21,14] = 0
 leSC: - lDX[21,15] + lDY[21,15] = 0
 leSC: - lDX[21,15] + ldy[21,15] = 0
 leSC: - lDX[22,0] + lDY[22,0] = 0
 leSC: - lDX[22,0] + ldy[22,0] = 0
 leSC: - lDX[22,1] + lDY[22,1] = 0
 leSC: - lDX[22,1] + ldy[22,1] = 0
 leSC: - lDX[22,2] + lDY[22,2] = 0
 leSC: - lDX[22,2] + ldy[22,2] = 0
 leSC: - lDX[22,3] + lDY[22,3] = 0
 leSC: - lDX[22,3] + ldy[22,3] = 0
 leSC: - lDX[22,4] + lDY[22,4] = 0
 leSC: - lDX[22,4] + ldy[22,4] = 0
 leSC: - lDX[22,5] + lDY[22,5] = 0
 leSC: - lDX[22,5] + ldy[22,5] = 0
 leSC: - lDX[22,6] + lDY[22,6] = 0
 leSC: - lDX[22,6] + ldy[22,6] = 0
 leSC: - lDX[22,7] + lDY[22,7] = 0
 leSC: - lDX[22,7] + ldy[22,7] = 0
 leSC: - lDX[22,8] + lDY[22,8] = 0
 leSC: - lDX[22,8] + ldy[22,8] = 0
 leSC: - lDX[22,9] + lDY[22,9] = 0
 leSC: - lDX[22,9] + ldy[22,9] = 0
 leSC: - lDX[22,10] + lDY[22,10] = 0
 leSC: - lDX[22,10] + ldy[22,10] = 0
 leSC: - lDX[22,11] + lDY[22,11] = 0
 leSC: - lDX[22,11] + ldy[22,11] = 0
 leSC: - lDX[22,12] + lDY[22,12] = 0
 leSC: - lDX[22,12] + ldy[22,12] = 0
 leSC: - lDX[22,13] + lDY[22,13] = 0
 leSC: - lDX[22,13] + ldy[22,13] = 0
 leSC: - lDX[22,14] + lDY[22,14] = 0
 leSC: - lDX[22,14] + ldy[22,14] = 0
 leSC: - lDX[22,15] + lDY[22,15] = 0
 leSC: - lDX[22,15] + ldy[22,15] = 0
 leSC: - lDX[23,0] + lDY[23,0] = 0
 leSC: - lDX[23,0] + ldy[23,0] = 0
 leSC: - lDX[23,1] + lDY[23,1] = 0
 leSC: - lDX[23,1] + ldy[23,1] = 0
 leSC: - lDX[23,2] + lDY[23,2] = 0
 leSC: - lDX[23,2] + ldy[23,2] = 0
 leSC: - lDX[23,3] + lDY[23,3] = 0
 leSC: - lDX[23,3] + ldy[23,3] = 0
 leSC: - lDX[23,4] + lDY[23,4] = 0
 leSC: - lDX[23,4] + ldy[23,4] = 0
 leSC: - lDX[23,5] + lDY[23,5] = 0
 leSC: - lDX[23,5] + ldy[23,5] = 0
 leSC: - lDX[23,6] + lDY[23,6] = 0
 leSC: - lDX[23,6] + ldy[23,6] = 0
 leSC: - lDX[23,7] + lDY[23,7] = 0
 leSC: - lDX[23,7] + ldy[23,7] = 0
 leSC: - lDX[23,8] + lDY[23,8] = 0
 leSC: - lDX[23,8] + ldy[23,8] = 0
 leSC: - lDX[23,9] + lDY[23,9] = 0
 leSC: - lDX[23,9] + ldy[23,9] = 0
 leSC: - lDX[23,10] + lDY[23,10] = 0
 leSC: - lDX[23,10] + ldy[23,10] = 0
 leSC: - lDX[23,11] + lDY[23,11] = 0
 leSC: - lDX[23,11] + ldy[23,11] = 0
 leSC: - lDX[23,12] + lDY[23,12] = 0
 leSC: - lDX[23,12] + ldy[23,12] = 0
 leSC: - lDX[23,13] + lDY[23,13] = 0
 leSC: - lDX[23,13] + ldy[23,13] = 0
 leSC: - lDX[23,14] + lDY[23,14] = 0
 leSC: - lDX[23,14] + ldy[23,14] = 0
 leSC: - lDX[23,15] + lDY[23,15] = 0
 leSC: - lDX[23,15] + ldy[23,15] = 0
 lATK: ldy[12,0] - ldz[12,0] = 0
 lATK: lDY[12,0] - lDZ[12,0] + lstk[39,0] >= 0
 lATK: - lDY[12,0] + lDZ[12,0] + lstk[39,0] >= 0
 lATK: lDY[12,0] + lDZ[12,0] - ldz[12,0] - lstk[39,0] >= 0
 lATK: lDY[12,0] - ldz[12,0] >= 0
 lATK: ldy[12,1] - ldz[12,1] = 0
 lATK: lDY[12,1] - lDZ[12,1] + lstk[39,1] >= 0
 lATK: - lDY[12,1] + lDZ[12,1] + lstk[39,1] >= 0
 lATK: lDY[12,1] + lDZ[12,1] - ldz[12,1] - lstk[39,1] >= 0
 lATK: lDY[12,1] - ldz[12,1] >= 0
 lATK: ldy[12,2] - ldz[12,2] = 0
 lATK: lDY[12,2] - lDZ[12,2] + lstk[39,2] >= 0
 lATK: - lDY[12,2] + lDZ[12,2] + lstk[39,2] >= 0
 lATK: lDY[12,2] + lDZ[12,2] - ldz[12,2] - lstk[39,2] >= 0
 lATK: lDY[12,2] - ldz[12,2] >= 0
 lATK: ldy[12,3] - ldz[12,3] = 0
 lATK: lDY[12,3] - lDZ[12,3] + lstk[39,3] >= 0
 lATK: - lDY[12,3] + lDZ[12,3] + lstk[39,3] >= 0
 lATK: lDY[12,3] + lDZ[12,3] - ldz[12,3] - lstk[39,3] >= 0
 lATK: lDY[12,3] - ldz[12,3] >= 0
 lATK: ldy[12,4] - ldz[12,4] = 0
 lATK: lDY[12,4] - lDZ[12,4] + lstk[39,4] >= 0
 lATK: - lDY[12,4] + lDZ[12,4] + lstk[39,4] >= 0
 lATK: lDY[12,4] + lDZ[12,4] - ldz[12,4] - lstk[39,4] >= 0
 lATK: lDY[12,4] - ldz[12,4] >= 0
 lATK: ldy[12,5] - ldz[12,5] = 0
 lATK: lDY[12,5] - lDZ[12,5] + lstk[39,5] >= 0
 lATK: - lDY[12,5] + lDZ[12,5] + lstk[39,5] >= 0
 lATK: lDY[12,5] + lDZ[12,5] - ldz[12,5] - lstk[39,5] >= 0
 lATK: lDY[12,5] - ldz[12,5] >= 0
 lATK: ldy[12,6] - ldz[12,6] = 0
 lATK: lDY[12,6] - lDZ[12,6] + lstk[39,6] >= 0
 lATK: - lDY[12,6] + lDZ[12,6] + lstk[39,6] >= 0
 lATK: lDY[12,6] + lDZ[12,6] - ldz[12,6] - lstk[39,6] >= 0
 lATK: lDY[12,6] - ldz[12,6] >= 0
 lATK: ldy[12,7] - ldz[12,7] = 0
 lATK: lDY[12,7] - lDZ[12,7] + lstk[39,7] >= 0
 lATK: - lDY[12,7] + lDZ[12,7] + lstk[39,7] >= 0
 lATK: lDY[12,7] + lDZ[12,7] - ldz[12,7] - lstk[39,7] >= 0
 lATK: lDY[12,7] - ldz[12,7] >= 0
 lATK: lDY[12,8] - lDZ[12,8] = 0
 lATK: ldy[12,8] - ldz[12,8] = 0
 lATK: lDY[12,9] - lDZ[12,9] = 0
 lATK: ldy[12,9] - ldz[12,9] = 0
 lATK: lDY[12,10] - lDZ[12,10] = 0
 lATK: ldy[12,10] - ldz[12,10] = 0
 lATK: lDY[12,11] - lDZ[12,11] = 0
 lATK: ldy[12,11] - ldz[12,11] = 0
 lATK: lDY[12,12] - lDZ[12,12] = 0
 lATK: ldy[12,12] - ldz[12,12] = 0
 lATK: lDY[12,13] - lDZ[12,13] = 0
 lATK: ldy[12,13] - ldz[12,13] = 0
 lATK: lDY[12,14] - lDZ[12,14] = 0
 lATK: ldy[12,14] - ldz[12,14] = 0
 lATK: lDY[12,15] - lDZ[12,15] = 0
 lATK: ldy[12,15] - ldz[12,15] = 0
 lATK: ldy[13,0] - ldz[13,0] = 0
 lATK: lDY[13,0] - lDZ[13,0] + lstk[40,0] >= 0
 lATK: - lDY[13,0] + lDZ[13,0] + lstk[40,0] >= 0
 lATK: lDY[13,0] + lDZ[13,0] - ldz[13,0] - lstk[40,0] >= 0
 lATK: lDY[13,0] - ldz[13,0] >= 0
 lATK: ldy[13,1] - ldz[13,1] = 0
 lATK: lDY[13,1] - lDZ[13,1] + lstk[40,1] >= 0
 lATK: - lDY[13,1] + lDZ[13,1] + lstk[40,1] >= 0
 lATK: lDY[13,1] + lDZ[13,1] - ldz[13,1] - lstk[40,1] >= 0
 lATK: lDY[13,1] - ldz[13,1] >= 0
 lATK: ldy[13,2] - ldz[13,2] = 0
 lATK: lDY[13,2] - lDZ[13,2] + lstk[40,2] >= 0
 lATK: - lDY[13,2] + lDZ[13,2] + lstk[40,2] >= 0
 lATK: lDY[13,2] + lDZ[13,2] - ldz[13,2] - lstk[40,2] >= 0
 lATK: lDY[13,2] - ldz[13,2] >= 0
 lATK: ldy[13,3] - ldz[13,3] = 0
 lATK: lDY[13,3] - lDZ[13,3] + lstk[40,3] >= 0
 lATK: - lDY[13,3] + lDZ[13,3] + lstk[40,3] >= 0
 lATK: lDY[13,3] + lDZ[13,3] - ldz[13,3] - lstk[40,3] >= 0
 lATK: lDY[13,3] - ldz[13,3] >= 0
 lATK: ldy[13,4] - ldz[13,4] = 0
 lATK: lDY[13,4] - lDZ[13,4] + lstk[40,4] >= 0
 lATK: - lDY[13,4] + lDZ[13,4] + lstk[40,4] >= 0
 lATK: lDY[13,4] + lDZ[13,4] - ldz[13,4] - lstk[40,4] >= 0
 lATK: lDY[13,4] - ldz[13,4] >= 0
 lATK: ldy[13,5] - ldz[13,5] = 0
 lATK: lDY[13,5] - lDZ[13,5] + lstk[40,5] >= 0
 lATK: - lDY[13,5] + lDZ[13,5] + lstk[40,5] >= 0
 lATK: lDY[13,5] + lDZ[13,5] - ldz[13,5] - lstk[40,5] >= 0
 lATK: lDY[13,5] - ldz[13,5] >= 0
 lATK: ldy[13,6] - ldz[13,6] = 0
 lATK: lDY[13,6] - lDZ[13,6] + lstk[40,6] >= 0
 lATK: - lDY[13,6] + lDZ[13,6] + lstk[40,6] >= 0
 lATK: lDY[13,6] + lDZ[13,6] - ldz[13,6] - lstk[40,6] >= 0
 lATK: lDY[13,6] - ldz[13,6] >= 0
 lATK: ldy[13,7] - ldz[13,7] = 0
 lATK: lDY[13,7] - lDZ[13,7] + lstk[40,7] >= 0
 lATK: - lDY[13,7] + lDZ[13,7] + lstk[40,7] >= 0
 lATK: lDY[13,7] + lDZ[13,7] - ldz[13,7] - lstk[40,7] >= 0
 lATK: lDY[13,7] - ldz[13,7] >= 0
 lATK: lDY[13,8] - lDZ[13,8] = 0
 lATK: ldy[13,8] - ldz[13,8] = 0
 lATK: lDY[13,9] - lDZ[13,9] = 0
 lATK: ldy[13,9] - ldz[13,9] = 0
 lATK: lDY[13,10] - lDZ[13,10] = 0
 lATK: ldy[13,10] - ldz[13,10] = 0
 lATK: lDY[13,11] - lDZ[13,11] = 0
 lATK: ldy[13,11] - ldz[13,11] = 0
 lATK: lDY[13,12] - lDZ[13,12] = 0
 lATK: ldy[13,12] - ldz[13,12] = 0
 lATK: lDY[13,13] - lDZ[13,13] = 0
 lATK: ldy[13,13] - ldz[13,13] = 0
 lATK: lDY[13,14] - lDZ[13,14] = 0
 lATK: ldy[13,14] - ldz[13,14] = 0
 lATK: lDY[13,15] - lDZ[13,15] = 0
 lATK: ldy[13,15] - ldz[13,15] = 0
 lATK: ldy[14,0] - ldz[14,0] = 0
 lATK: lDY[14,0] - lDZ[14,0] + lstk[41,0] >= 0
 lATK: - lDY[14,0] + lDZ[14,0] + lstk[41,0] >= 0
 lATK: lDY[14,0] + lDZ[14,0] - ldz[14,0] - lstk[41,0] >= 0
 lATK: lDY[14,0] - ldz[14,0] >= 0
 lATK: ldy[14,1] - ldz[14,1] = 0
 lATK: lDY[14,1] - lDZ[14,1] + lstk[41,1] >= 0
 lATK: - lDY[14,1] + lDZ[14,1] + lstk[41,1] >= 0
 lATK: lDY[14,1] + lDZ[14,1] - ldz[14,1] - lstk[41,1] >= 0
 lATK: lDY[14,1] - ldz[14,1] >= 0
 lATK: ldy[14,2] - ldz[14,2] = 0
 lATK: lDY[14,2] - lDZ[14,2] + lstk[41,2] >= 0
 lATK: - lDY[14,2] + lDZ[14,2] + lstk[41,2] >= 0
 lATK: lDY[14,2] + lDZ[14,2] - ldz[14,2] - lstk[41,2] >= 0
 lATK: lDY[14,2] - ldz[14,2] >= 0
 lATK: ldy[14,3] - ldz[14,3] = 0
 lATK: lDY[14,3] - lDZ[14,3] + lstk[41,3] >= 0
 lATK: - lDY[14,3] + lDZ[14,3] + lstk[41,3] >= 0
 lATK: lDY[14,3] + lDZ[14,3] - ldz[14,3] - lstk[41,3] >= 0
 lATK: lDY[14,3] - ldz[14,3] >= 0
 lATK: ldy[14,4] - ldz[14,4] = 0
 lATK: lDY[14,4] - lDZ[14,4] + lstk[41,4] >= 0
 lATK: - lDY[14,4] + lDZ[14,4] + lstk[41,4] >= 0
 lATK: lDY[14,4] + lDZ[14,4] - ldz[14,4] - lstk[41,4] >= 0
 lATK: lDY[14,4] - ldz[14,4] >= 0
 lATK: ldy[14,5] - ldz[14,5] = 0
 lATK: lDY[14,5] - lDZ[14,5] + lstk[41,5] >= 0
 lATK: - lDY[14,5] + lDZ[14,5] + lstk[41,5] >= 0
 lATK: lDY[14,5] + lDZ[14,5] - ldz[14,5] - lstk[41,5] >= 0
 lATK: lDY[14,5] - ldz[14,5] >= 0
 lATK: ldy[14,6] - ldz[14,6] = 0
 lATK: lDY[14,6] - lDZ[14,6] + lstk[41,6] >= 0
 lATK: - lDY[14,6] + lDZ[14,6] + lstk[41,6] >= 0
 lATK: lDY[14,6] + lDZ[14,6] - ldz[14,6] - lstk[41,6] >= 0
 lATK: lDY[14,6] - ldz[14,6] >= 0
 lATK: ldy[14,7] - ldz[14,7] = 0
 lATK: lDY[14,7] - lDZ[14,7] + lstk[41,7] >= 0
 lATK: - lDY[14,7] + lDZ[14,7] + lstk[41,7] >= 0
 lATK: lDY[14,7] + lDZ[14,7] - ldz[14,7] - lstk[41,7] >= 0
 lATK: lDY[14,7] - ldz[14,7] >= 0
 lATK: lDY[14,8] - lDZ[14,8] = 0
 lATK: ldy[14,8] - ldz[14,8] = 0
 lATK: lDY[14,9] - lDZ[14,9] = 0
 lATK: ldy[14,9] - ldz[14,9] = 0
 lATK: lDY[14,10] - lDZ[14,10] = 0
 lATK: ldy[14,10] - ldz[14,10] = 0
 lATK: lDY[14,11] - lDZ[14,11] = 0
 lATK: ldy[14,11] - ldz[14,11] = 0
 lATK: lDY[14,12] - lDZ[14,12] = 0
 lATK: ldy[14,12] - ldz[14,12] = 0
 lATK: lDY[14,13] - lDZ[14,13] = 0
 lATK: ldy[14,13] - ldz[14,13] = 0
 lATK: lDY[14,14] - lDZ[14,14] = 0
 lATK: ldy[14,14] - ldz[14,14] = 0
 lATK: lDY[14,15] - lDZ[14,15] = 0
 lATK: ldy[14,15] - ldz[14,15] = 0
 lATK: ldy[15,0] - ldz[15,0] = 0
 lATK: lDY[15,0] - lDZ[15,0] + lstk[42,0] >= 0
 lATK: - lDY[15,0] + lDZ[15,0] + lstk[42,0] >= 0
 lATK: lDY[15,0] + lDZ[15,0] - ldz[15,0] - lstk[42,0] >= 0
 lATK: lDY[15,0] - ldz[15,0] >= 0
 lATK: ldy[15,1] - ldz[15,1] = 0
 lATK: lDY[15,1] - lDZ[15,1] + lstk[42,1] >= 0
 lATK: - lDY[15,1] + lDZ[15,1] + lstk[42,1] >= 0
 lATK: lDY[15,1] + lDZ[15,1] - ldz[15,1] - lstk[42,1] >= 0
 lATK: lDY[15,1] - ldz[15,1] >= 0
 lATK: ldy[15,2] - ldz[15,2] = 0
 lATK: lDY[15,2] - lDZ[15,2] + lstk[42,2] >= 0
 lATK: - lDY[15,2] + lDZ[15,2] + lstk[42,2] >= 0
 lATK: lDY[15,2] + lDZ[15,2] - ldz[15,2] - lstk[42,2] >= 0
 lATK: lDY[15,2] - ldz[15,2] >= 0
 lATK: ldy[15,3] - ldz[15,3] = 0
 lATK: lDY[15,3] - lDZ[15,3] + lstk[42,3] >= 0
 lATK: - lDY[15,3] + lDZ[15,3] + lstk[42,3] >= 0
 lATK: lDY[15,3] + lDZ[15,3] - ldz[15,3] - lstk[42,3] >= 0
 lATK: lDY[15,3] - ldz[15,3] >= 0
 lATK: ldy[15,4] - ldz[15,4] = 0
 lATK: lDY[15,4] - lDZ[15,4] + lstk[42,4] >= 0
 lATK: - lDY[15,4] + lDZ[15,4] + lstk[42,4] >= 0
 lATK: lDY[15,4] + lDZ[15,4] - ldz[15,4] - lstk[42,4] >= 0
 lATK: lDY[15,4] - ldz[15,4] >= 0
 lATK: ldy[15,5] - ldz[15,5] = 0
 lATK: lDY[15,5] - lDZ[15,5] + lstk[42,5] >= 0
 lATK: - lDY[15,5] + lDZ[15,5] + lstk[42,5] >= 0
 lATK: lDY[15,5] + lDZ[15,5] - ldz[15,5] - lstk[42,5] >= 0
 lATK: lDY[15,5] - ldz[15,5] >= 0
 lATK: ldy[15,6] - ldz[15,6] = 0
 lATK: lDY[15,6] - lDZ[15,6] + lstk[42,6] >= 0
 lATK: - lDY[15,6] + lDZ[15,6] + lstk[42,6] >= 0
 lATK: lDY[15,6] + lDZ[15,6] - ldz[15,6] - lstk[42,6] >= 0
 lATK: lDY[15,6] - ldz[15,6] >= 0
 lATK: ldy[15,7] - ldz[15,7] = 0
 lATK: lDY[15,7] - lDZ[15,7] + lstk[42,7] >= 0
 lATK: - lDY[15,7] + lDZ[15,7] + lstk[42,7] >= 0
 lATK: lDY[15,7] + lDZ[15,7] - ldz[15,7] - lstk[42,7] >= 0
 lATK: lDY[15,7] - ldz[15,7] >= 0
 lATK: lDY[15,8] - lDZ[15,8] = 0
 lATK: ldy[15,8] - ldz[15,8] = 0
 lATK: lDY[15,9] - lDZ[15,9] = 0
 lATK: ldy[15,9] - ldz[15,9] = 0
 lATK: lDY[15,10] - lDZ[15,10] = 0
 lATK: ldy[15,10] - ldz[15,10] = 0
 lATK: lDY[15,11] - lDZ[15,11] = 0
 lATK: ldy[15,11] - ldz[15,11] = 0
 lATK: lDY[15,12] - lDZ[15,12] = 0
 lATK: ldy[15,12] - ldz[15,12] = 0
 lATK: lDY[15,13] - lDZ[15,13] = 0
 lATK: ldy[15,13] - ldz[15,13] = 0
 lATK: lDY[15,14] - lDZ[15,14] = 0
 lATK: ldy[15,14] - ldz[15,14] = 0
 lATK: lDY[15,15] - lDZ[15,15] = 0
 lATK: ldy[15,15] - ldz[15,15] = 0
 lATK: ldy[16,0] - ldz[16,0] = 0
 lATK: lDY[16,0] - lDZ[16,0] + lstk[43,0] >= 0
 lATK: - lDY[16,0] + lDZ[16,0] + lstk[43,0] >= 0
 lATK: lDY[16,0] + lDZ[16,0] - ldz[16,0] - lstk[43,0] >= 0
 lATK: lDY[16,0] - ldz[16,0] >= 0
 lATK: ldy[16,1] - ldz[16,1] = 0
 lATK: lDY[16,1] - lDZ[16,1] + lstk[43,1] >= 0
 lATK: - lDY[16,1] + lDZ[16,1] + lstk[43,1] >= 0
 lATK: lDY[16,1] + lDZ[16,1] - ldz[16,1] - lstk[43,1] >= 0
 lATK: lDY[16,1] - ldz[16,1] >= 0
 lATK: ldy[16,2] - ldz[16,2] = 0
 lATK: lDY[16,2] - lDZ[16,2] + lstk[43,2] >= 0
 lATK: - lDY[16,2] + lDZ[16,2] + lstk[43,2] >= 0
 lATK: lDY[16,2] + lDZ[16,2] - ldz[16,2] - lstk[43,2] >= 0
 lATK: lDY[16,2] - ldz[16,2] >= 0
 lATK: ldy[16,3] - ldz[16,3] = 0
 lATK: lDY[16,3] - lDZ[16,3] + lstk[43,3] >= 0
 lATK: - lDY[16,3] + lDZ[16,3] + lstk[43,3] >= 0
 lATK: lDY[16,3] + lDZ[16,3] - ldz[16,3] - lstk[43,3] >= 0
 lATK: lDY[16,3] - ldz[16,3] >= 0
 lATK: ldy[16,4] - ldz[16,4] = 0
 lATK: lDY[16,4] - lDZ[16,4] + lstk[43,4] >= 0
 lATK: - lDY[16,4] + lDZ[16,4] + lstk[43,4] >= 0
 lATK: lDY[16,4] + lDZ[16,4] - ldz[16,4] - lstk[43,4] >= 0
 lATK: lDY[16,4] - ldz[16,4] >= 0
 lATK: ldy[16,5] - ldz[16,5] = 0
 lATK: lDY[16,5] - lDZ[16,5] + lstk[43,5] >= 0
 lATK: - lDY[16,5] + lDZ[16,5] + lstk[43,5] >= 0
 lATK: lDY[16,5] + lDZ[16,5] - ldz[16,5] - lstk[43,5] >= 0
 lATK: lDY[16,5] - ldz[16,5] >= 0
 lATK: ldy[16,6] - ldz[16,6] = 0
 lATK: lDY[16,6] - lDZ[16,6] + lstk[43,6] >= 0
 lATK: - lDY[16,6] + lDZ[16,6] + lstk[43,6] >= 0
 lATK: lDY[16,6] + lDZ[16,6] - ldz[16,6] - lstk[43,6] >= 0
 lATK: lDY[16,6] - ldz[16,6] >= 0
 lATK: ldy[16,7] - ldz[16,7] = 0
 lATK: lDY[16,7] - lDZ[16,7] + lstk[43,7] >= 0
 lATK: - lDY[16,7] + lDZ[16,7] + lstk[43,7] >= 0
 lATK: lDY[16,7] + lDZ[16,7] - ldz[16,7] - lstk[43,7] >= 0
 lATK: lDY[16,7] - ldz[16,7] >= 0
 lATK: lDY[16,8] - lDZ[16,8] = 0
 lATK: ldy[16,8] - ldz[16,8] = 0
 lATK: lDY[16,9] - lDZ[16,9] = 0
 lATK: ldy[16,9] - ldz[16,9] = 0
 lATK: lDY[16,10] - lDZ[16,10] = 0
 lATK: ldy[16,10] - ldz[16,10] = 0
 lATK: lDY[16,11] - lDZ[16,11] = 0
 lATK: ldy[16,11] - ldz[16,11] = 0
 lATK: lDY[16,12] - lDZ[16,12] = 0
 lATK: ldy[16,12] - ldz[16,12] = 0
 lATK: lDY[16,13] - lDZ[16,13] = 0
 lATK: ldy[16,13] - ldz[16,13] = 0
 lATK: lDY[16,14] - lDZ[16,14] = 0
 lATK: ldy[16,14] - ldz[16,14] = 0
 lATK: lDY[16,15] - lDZ[16,15] = 0
 lATK: ldy[16,15] - ldz[16,15] = 0
 lATK: ldy[17,0] - ldz[17,0] = 0
 lATK: lDY[17,0] - lDZ[17,0] + lstk[44,0] >= 0
 lATK: - lDY[17,0] + lDZ[17,0] + lstk[44,0] >= 0
 lATK: lDY[17,0] + lDZ[17,0] - ldz[17,0] - lstk[44,0] >= 0
 lATK: lDY[17,0] - ldz[17,0] >= 0
 lATK: ldy[17,1] - ldz[17,1] = 0
 lATK: lDY[17,1] - lDZ[17,1] + lstk[44,1] >= 0
 lATK: - lDY[17,1] + lDZ[17,1] + lstk[44,1] >= 0
 lATK: lDY[17,1] + lDZ[17,1] - ldz[17,1] - lstk[44,1] >= 0
 lATK: lDY[17,1] - ldz[17,1] >= 0
 lATK: ldy[17,2] - ldz[17,2] = 0
 lATK: lDY[17,2] - lDZ[17,2] + lstk[44,2] >= 0
 lATK: - lDY[17,2] + lDZ[17,2] + lstk[44,2] >= 0
 lATK: lDY[17,2] + lDZ[17,2] - ldz[17,2] - lstk[44,2] >= 0
 lATK: lDY[17,2] - ldz[17,2] >= 0
 lATK: ldy[17,3] - ldz[17,3] = 0
 lATK: lDY[17,3] - lDZ[17,3] + lstk[44,3] >= 0
 lATK: - lDY[17,3] + lDZ[17,3] + lstk[44,3] >= 0
 lATK: lDY[17,3] + lDZ[17,3] - ldz[17,3] - lstk[44,3] >= 0
 lATK: lDY[17,3] - ldz[17,3] >= 0
 lATK: ldy[17,4] - ldz[17,4] = 0
 lATK: lDY[17,4] - lDZ[17,4] + lstk[44,4] >= 0
 lATK: - lDY[17,4] + lDZ[17,4] + lstk[44,4] >= 0
 lATK: lDY[17,4] + lDZ[17,4] - ldz[17,4] - lstk[44,4] >= 0
 lATK: lDY[17,4] - ldz[17,4] >= 0
 lATK: ldy[17,5] - ldz[17,5] = 0
 lATK: lDY[17,5] - lDZ[17,5] + lstk[44,5] >= 0
 lATK: - lDY[17,5] + lDZ[17,5] + lstk[44,5] >= 0
 lATK: lDY[17,5] + lDZ[17,5] - ldz[17,5] - lstk[44,5] >= 0
 lATK: lDY[17,5] - ldz[17,5] >= 0
 lATK: ldy[17,6] - ldz[17,6] = 0
 lATK: lDY[17,6] - lDZ[17,6] + lstk[44,6] >= 0
 lATK: - lDY[17,6] + lDZ[17,6] + lstk[44,6] >= 0
 lATK: lDY[17,6] + lDZ[17,6] - ldz[17,6] - lstk[44,6] >= 0
 lATK: lDY[17,6] - ldz[17,6] >= 0
 lATK: ldy[17,7] - ldz[17,7] = 0
 lATK: lDY[17,7] - lDZ[17,7] + lstk[44,7] >= 0
 lATK: - lDY[17,7] + lDZ[17,7] + lstk[44,7] >= 0
 lATK: lDY[17,7] + lDZ[17,7] - ldz[17,7] - lstk[44,7] >= 0
 lATK: lDY[17,7] - ldz[17,7] >= 0
 lATK: lDY[17,8] - lDZ[17,8] = 0
 lATK: ldy[17,8] - ldz[17,8] = 0
 lATK: lDY[17,9] - lDZ[17,9] = 0
 lATK: ldy[17,9] - ldz[17,9] = 0
 lATK: lDY[17,10] - lDZ[17,10] = 0
 lATK: ldy[17,10] - ldz[17,10] = 0
 lATK: lDY[17,11] - lDZ[17,11] = 0
 lATK: ldy[17,11] - ldz[17,11] = 0
 lATK: lDY[17,12] - lDZ[17,12] = 0
 lATK: ldy[17,12] - ldz[17,12] = 0
 lATK: lDY[17,13] - lDZ[17,13] = 0
 lATK: ldy[17,13] - ldz[17,13] = 0
 lATK: lDY[17,14] - lDZ[17,14] = 0
 lATK: ldy[17,14] - ldz[17,14] = 0
 lATK: lDY[17,15] - lDZ[17,15] = 0
 lATK: ldy[17,15] - ldz[17,15] = 0
 lATK: ldy[18,0] - ldz[18,0] = 0
 lATK: lDY[18,0] - lDZ[18,0] + lstk[45,0] >= 0
 lATK: - lDY[18,0] + lDZ[18,0] + lstk[45,0] >= 0
 lATK: lDY[18,0] + lDZ[18,0] - ldz[18,0] - lstk[45,0] >= 0
 lATK: lDY[18,0] - ldz[18,0] >= 0
 lATK: ldy[18,1] - ldz[18,1] = 0
 lATK: lDY[18,1] - lDZ[18,1] + lstk[45,1] >= 0
 lATK: - lDY[18,1] + lDZ[18,1] + lstk[45,1] >= 0
 lATK: lDY[18,1] + lDZ[18,1] - ldz[18,1] - lstk[45,1] >= 0
 lATK: lDY[18,1] - ldz[18,1] >= 0
 lATK: ldy[18,2] - ldz[18,2] = 0
 lATK: lDY[18,2] - lDZ[18,2] + lstk[45,2] >= 0
 lATK: - lDY[18,2] + lDZ[18,2] + lstk[45,2] >= 0
 lATK: lDY[18,2] + lDZ[18,2] - ldz[18,2] - lstk[45,2] >= 0
 lATK: lDY[18,2] - ldz[18,2] >= 0
 lATK: ldy[18,3] - ldz[18,3] = 0
 lATK: lDY[18,3] - lDZ[18,3] + lstk[45,3] >= 0
 lATK: - lDY[18,3] + lDZ[18,3] + lstk[45,3] >= 0
 lATK: lDY[18,3] + lDZ[18,3] - ldz[18,3] - lstk[45,3] >= 0
 lATK: lDY[18,3] - ldz[18,3] >= 0
 lATK: ldy[18,4] - ldz[18,4] = 0
 lATK: lDY[18,4] - lDZ[18,4] + lstk[45,4] >= 0
 lATK: - lDY[18,4] + lDZ[18,4] + lstk[45,4] >= 0
 lATK: lDY[18,4] + lDZ[18,4] - ldz[18,4] - lstk[45,4] >= 0
 lATK: lDY[18,4] - ldz[18,4] >= 0
 lATK: ldy[18,5] - ldz[18,5] = 0
 lATK: lDY[18,5] - lDZ[18,5] + lstk[45,5] >= 0
 lATK: - lDY[18,5] + lDZ[18,5] + lstk[45,5] >= 0
 lATK: lDY[18,5] + lDZ[18,5] - ldz[18,5] - lstk[45,5] >= 0
 lATK: lDY[18,5] - ldz[18,5] >= 0
 lATK: ldy[18,6] - ldz[18,6] = 0
 lATK: lDY[18,6] - lDZ[18,6] + lstk[45,6] >= 0
 lATK: - lDY[18,6] + lDZ[18,6] + lstk[45,6] >= 0
 lATK: lDY[18,6] + lDZ[18,6] - ldz[18,6] - lstk[45,6] >= 0
 lATK: lDY[18,6] - ldz[18,6] >= 0
 lATK: ldy[18,7] - ldz[18,7] = 0
 lATK: lDY[18,7] - lDZ[18,7] + lstk[45,7] >= 0
 lATK: - lDY[18,7] + lDZ[18,7] + lstk[45,7] >= 0
 lATK: lDY[18,7] + lDZ[18,7] - ldz[18,7] - lstk[45,7] >= 0
 lATK: lDY[18,7] - ldz[18,7] >= 0
 lATK: lDY[18,8] - lDZ[18,8] = 0
 lATK: ldy[18,8] - ldz[18,8] = 0
 lATK: lDY[18,9] - lDZ[18,9] = 0
 lATK: ldy[18,9] - ldz[18,9] = 0
 lATK: lDY[18,10] - lDZ[18,10] = 0
 lATK: ldy[18,10] - ldz[18,10] = 0
 lATK: lDY[18,11] - lDZ[18,11] = 0
 lATK: ldy[18,11] - ldz[18,11] = 0
 lATK: lDY[18,12] - lDZ[18,12] = 0
 lATK: ldy[18,12] - ldz[18,12] = 0
 lATK: lDY[18,13] - lDZ[18,13] = 0
 lATK: ldy[18,13] - ldz[18,13] = 0
 lATK: lDY[18,14] - lDZ[18,14] = 0
 lATK: ldy[18,14] - ldz[18,14] = 0
 lATK: lDY[18,15] - lDZ[18,15] = 0
 lATK: ldy[18,15] - ldz[18,15] = 0
 lATK: ldy[19,0] - ldz[19,0] = 0
 lATK: lDY[19,0] - lDZ[19,0] + lstk[46,0] >= 0
 lATK: - lDY[19,0] + lDZ[19,0] + lstk[46,0] >= 0
 lATK: lDY[19,0] + lDZ[19,0] - ldz[19,0] - lstk[46,0] >= 0
 lATK: lDY[19,0] - ldz[19,0] >= 0
 lATK: ldy[19,1] - ldz[19,1] = 0
 lATK: lDY[19,1] - lDZ[19,1] + lstk[46,1] >= 0
 lATK: - lDY[19,1] + lDZ[19,1] + lstk[46,1] >= 0
 lATK: lDY[19,1] + lDZ[19,1] - ldz[19,1] - lstk[46,1] >= 0
 lATK: lDY[19,1] - ldz[19,1] >= 0
 lATK: ldy[19,2] - ldz[19,2] = 0
 lATK: lDY[19,2] - lDZ[19,2] + lstk[46,2] >= 0
 lATK: - lDY[19,2] + lDZ[19,2] + lstk[46,2] >= 0
 lATK: lDY[19,2] + lDZ[19,2] - ldz[19,2] - lstk[46,2] >= 0
 lATK: lDY[19,2] - ldz[19,2] >= 0
 lATK: ldy[19,3] - ldz[19,3] = 0
 lATK: lDY[19,3] - lDZ[19,3] + lstk[46,3] >= 0
 lATK: - lDY[19,3] + lDZ[19,3] + lstk[46,3] >= 0
 lATK: lDY[19,3] + lDZ[19,3] - ldz[19,3] - lstk[46,3] >= 0
 lATK: lDY[19,3] - ldz[19,3] >= 0
 lATK: ldy[19,4] - ldz[19,4] = 0
 lATK: lDY[19,4] - lDZ[19,4] + lstk[46,4] >= 0
 lATK: - lDY[19,4] + lDZ[19,4] + lstk[46,4] >= 0
 lATK: lDY[19,4] + lDZ[19,4] - ldz[19,4] - lstk[46,4] >= 0
 lATK: lDY[19,4] - ldz[19,4] >= 0
 lATK: ldy[19,5] - ldz[19,5] = 0
 lATK: lDY[19,5] - lDZ[19,5] + lstk[46,5] >= 0
 lATK: - lDY[19,5] + lDZ[19,5] + lstk[46,5] >= 0
 lATK: lDY[19,5] + lDZ[19,5] - ldz[19,5] - lstk[46,5] >= 0
 lATK: lDY[19,5] - ldz[19,5] >= 0
 lATK: ldy[19,6] - ldz[19,6] = 0
 lATK: lDY[19,6] - lDZ[19,6] + lstk[46,6] >= 0
 lATK: - lDY[19,6] + lDZ[19,6] + lstk[46,6] >= 0
 lATK: lDY[19,6] + lDZ[19,6] - ldz[19,6] - lstk[46,6] >= 0
 lATK: lDY[19,6] - ldz[19,6] >= 0
 lATK: ldy[19,7] - ldz[19,7] = 0
 lATK: lDY[19,7] - lDZ[19,7] + lstk[46,7] >= 0
 lATK: - lDY[19,7] + lDZ[19,7] + lstk[46,7] >= 0
 lATK: lDY[19,7] + lDZ[19,7] - ldz[19,7] - lstk[46,7] >= 0
 lATK: lDY[19,7] - ldz[19,7] >= 0
 lATK: lDY[19,8] - lDZ[19,8] = 0
 lATK: ldy[19,8] - ldz[19,8] = 0
 lATK: lDY[19,9] - lDZ[19,9] = 0
 lATK: ldy[19,9] - ldz[19,9] = 0
 lATK: lDY[19,10] - lDZ[19,10] = 0
 lATK: ldy[19,10] - ldz[19,10] = 0
 lATK: lDY[19,11] - lDZ[19,11] = 0
 lATK: ldy[19,11] - ldz[19,11] = 0
 lATK: lDY[19,12] - lDZ[19,12] = 0
 lATK: ldy[19,12] - ldz[19,12] = 0
 lATK: lDY[19,13] - lDZ[19,13] = 0
 lATK: ldy[19,13] - ldz[19,13] = 0
 lATK: lDY[19,14] - lDZ[19,14] = 0
 lATK: ldy[19,14] - ldz[19,14] = 0
 lATK: lDY[19,15] - lDZ[19,15] = 0
 lATK: ldy[19,15] - ldz[19,15] = 0
 lATK: ldy[20,0] - ldz[20,0] = 0
 lATK: lDY[20,0] - lDZ[20,0] + lstk[47,0] >= 0
 lATK: - lDY[20,0] + lDZ[20,0] + lstk[47,0] >= 0
 lATK: lDY[20,0] + lDZ[20,0] - ldz[20,0] - lstk[47,0] >= 0
 lATK: lDY[20,0] - ldz[20,0] >= 0
 lATK: ldy[20,1] - ldz[20,1] = 0
 lATK: lDY[20,1] - lDZ[20,1] + lstk[47,1] >= 0
 lATK: - lDY[20,1] + lDZ[20,1] + lstk[47,1] >= 0
 lATK: lDY[20,1] + lDZ[20,1] - ldz[20,1] - lstk[47,1] >= 0
 lATK: lDY[20,1] - ldz[20,1] >= 0
 lATK: ldy[20,2] - ldz[20,2] = 0
 lATK: lDY[20,2] - lDZ[20,2] + lstk[47,2] >= 0
 lATK: - lDY[20,2] + lDZ[20,2] + lstk[47,2] >= 0
 lATK: lDY[20,2] + lDZ[20,2] - ldz[20,2] - lstk[47,2] >= 0
 lATK: lDY[20,2] - ldz[20,2] >= 0
 lATK: ldy[20,3] - ldz[20,3] = 0
 lATK: lDY[20,3] - lDZ[20,3] + lstk[47,3] >= 0
 lATK: - lDY[20,3] + lDZ[20,3] + lstk[47,3] >= 0
 lATK: lDY[20,3] + lDZ[20,3] - ldz[20,3] - lstk[47,3] >= 0
 lATK: lDY[20,3] - ldz[20,3] >= 0
 lATK: ldy[20,4] - ldz[20,4] = 0
 lATK: lDY[20,4] - lDZ[20,4] + lstk[47,4] >= 0
 lATK: - lDY[20,4] + lDZ[20,4] + lstk[47,4] >= 0
 lATK: lDY[20,4] + lDZ[20,4] - ldz[20,4] - lstk[47,4] >= 0
 lATK: lDY[20,4] - ldz[20,4] >= 0
 lATK: ldy[20,5] - ldz[20,5] = 0
 lATK: lDY[20,5] - lDZ[20,5] + lstk[47,5] >= 0
 lATK: - lDY[20,5] + lDZ[20,5] + lstk[47,5] >= 0
 lATK: lDY[20,5] + lDZ[20,5] - ldz[20,5] - lstk[47,5] >= 0
 lATK: lDY[20,5] - ldz[20,5] >= 0
 lATK: ldy[20,6] - ldz[20,6] = 0
 lATK: lDY[20,6] - lDZ[20,6] + lstk[47,6] >= 0
 lATK: - lDY[20,6] + lDZ[20,6] + lstk[47,6] >= 0
 lATK: lDY[20,6] + lDZ[20,6] - ldz[20,6] - lstk[47,6] >= 0
 lATK: lDY[20,6] - ldz[20,6] >= 0
 lATK: ldy[20,7] - ldz[20,7] = 0
 lATK: lDY[20,7] - lDZ[20,7] + lstk[47,7] >= 0
 lATK: - lDY[20,7] + lDZ[20,7] + lstk[47,7] >= 0
 lATK: lDY[20,7] + lDZ[20,7] - ldz[20,7] - lstk[47,7] >= 0
 lATK: lDY[20,7] - ldz[20,7] >= 0
 lATK: lDY[20,8] - lDZ[20,8] = 0
 lATK: ldy[20,8] - ldz[20,8] = 0
 lATK: lDY[20,9] - lDZ[20,9] = 0
 lATK: ldy[20,9] - ldz[20,9] = 0
 lATK: lDY[20,10] - lDZ[20,10] = 0
 lATK: ldy[20,10] - ldz[20,10] = 0
 lATK: lDY[20,11] - lDZ[20,11] = 0
 lATK: ldy[20,11] - ldz[20,11] = 0
 lATK: lDY[20,12] - lDZ[20,12] = 0
 lATK: ldy[20,12] - ldz[20,12] = 0
 lATK: lDY[20,13] - lDZ[20,13] = 0
 lATK: ldy[20,13] - ldz[20,13] = 0
 lATK: lDY[20,14] - lDZ[20,14] = 0
 lATK: ldy[20,14] - ldz[20,14] = 0
 lATK: lDY[20,15] - lDZ[20,15] = 0
 lATK: ldy[20,15] - ldz[20,15] = 0
 lATK: ldy[21,0] - ldz[21,0] = 0
 lATK: lDY[21,0] - lDZ[21,0] + lstk[48,0] >= 0
 lATK: - lDY[21,0] + lDZ[21,0] + lstk[48,0] >= 0
 lATK: lDY[21,0] + lDZ[21,0] - ldz[21,0] - lstk[48,0] >= 0
 lATK: lDY[21,0] - ldz[21,0] >= 0
 lATK: ldy[21,1] - ldz[21,1] = 0
 lATK: lDY[21,1] - lDZ[21,1] + lstk[48,1] >= 0
 lATK: - lDY[21,1] + lDZ[21,1] + lstk[48,1] >= 0
 lATK: lDY[21,1] + lDZ[21,1] - ldz[21,1] - lstk[48,1] >= 0
 lATK: lDY[21,1] - ldz[21,1] >= 0
 lATK: ldy[21,2] - ldz[21,2] = 0
 lATK: lDY[21,2] - lDZ[21,2] + lstk[48,2] >= 0
 lATK: - lDY[21,2] + lDZ[21,2] + lstk[48,2] >= 0
 lATK: lDY[21,2] + lDZ[21,2] - ldz[21,2] - lstk[48,2] >= 0
 lATK: lDY[21,2] - ldz[21,2] >= 0
 lATK: ldy[21,3] - ldz[21,3] = 0
 lATK: lDY[21,3] - lDZ[21,3] + lstk[48,3] >= 0
 lATK: - lDY[21,3] + lDZ[21,3] + lstk[48,3] >= 0
 lATK: lDY[21,3] + lDZ[21,3] - ldz[21,3] - lstk[48,3] >= 0
 lATK: lDY[21,3] - ldz[21,3] >= 0
 lATK: ldy[21,4] - ldz[21,4] = 0
 lATK: lDY[21,4] - lDZ[21,4] + lstk[48,4] >= 0
 lATK: - lDY[21,4] + lDZ[21,4] + lstk[48,4] >= 0
 lATK: lDY[21,4] + lDZ[21,4] - ldz[21,4] - lstk[48,4] >= 0
 lATK: lDY[21,4] - ldz[21,4] >= 0
 lATK: ldy[21,5] - ldz[21,5] = 0
 lATK: lDY[21,5] - lDZ[21,5] + lstk[48,5] >= 0
 lATK: - lDY[21,5] + lDZ[21,5] + lstk[48,5] >= 0
 lATK: lDY[21,5] + lDZ[21,5] - ldz[21,5] - lstk[48,5] >= 0
 lATK: lDY[21,5] - ldz[21,5] >= 0
 lATK: ldy[21,6] - ldz[21,6] = 0
 lATK: lDY[21,6] - lDZ[21,6] + lstk[48,6] >= 0
 lATK: - lDY[21,6] + lDZ[21,6] + lstk[48,6] >= 0
 lATK: lDY[21,6] + lDZ[21,6] - ldz[21,6] - lstk[48,6] >= 0
 lATK: lDY[21,6] - ldz[21,6] >= 0
 lATK: ldy[21,7] - ldz[21,7] = 0
 lATK: lDY[21,7] - lDZ[21,7] + lstk[48,7] >= 0
 lATK: - lDY[21,7] + lDZ[21,7] + lstk[48,7] >= 0
 lATK: lDY[21,7] + lDZ[21,7] - ldz[21,7] - lstk[48,7] >= 0
 lATK: lDY[21,7] - ldz[21,7] >= 0
 lATK: lDY[21,8] - lDZ[21,8] = 0
 lATK: ldy[21,8] - ldz[21,8] = 0
 lATK: lDY[21,9] - lDZ[21,9] = 0
 lATK: ldy[21,9] - ldz[21,9] = 0
 lATK: lDY[21,10] - lDZ[21,10] = 0
 lATK: ldy[21,10] - ldz[21,10] = 0
 lATK: lDY[21,11] - lDZ[21,11] = 0
 lATK: ldy[21,11] - ldz[21,11] = 0
 lATK: lDY[21,12] - lDZ[21,12] = 0
 lATK: ldy[21,12] - ldz[21,12] = 0
 lATK: lDY[21,13] - lDZ[21,13] = 0
 lATK: ldy[21,13] - ldz[21,13] = 0
 lATK: lDY[21,14] - lDZ[21,14] = 0
 lATK: ldy[21,14] - ldz[21,14] = 0
 lATK: lDY[21,15] - lDZ[21,15] = 0
 lATK: ldy[21,15] - ldz[21,15] = 0
 lATK: ldy[22,0] - ldz[22,0] = 0
 lATK: lDY[22,0] - lDZ[22,0] + lstk[49,0] >= 0
 lATK: - lDY[22,0] + lDZ[22,0] + lstk[49,0] >= 0
 lATK: lDY[22,0] + lDZ[22,0] - ldz[22,0] - lstk[49,0] >= 0
 lATK: lDY[22,0] - ldz[22,0] >= 0
 lATK: ldy[22,1] - ldz[22,1] = 0
 lATK: lDY[22,1] - lDZ[22,1] + lstk[49,1] >= 0
 lATK: - lDY[22,1] + lDZ[22,1] + lstk[49,1] >= 0
 lATK: lDY[22,1] + lDZ[22,1] - ldz[22,1] - lstk[49,1] >= 0
 lATK: lDY[22,1] - ldz[22,1] >= 0
 lATK: ldy[22,2] - ldz[22,2] = 0
 lATK: lDY[22,2] - lDZ[22,2] + lstk[49,2] >= 0
 lATK: - lDY[22,2] + lDZ[22,2] + lstk[49,2] >= 0
 lATK: lDY[22,2] + lDZ[22,2] - ldz[22,2] - lstk[49,2] >= 0
 lATK: lDY[22,2] - ldz[22,2] >= 0
 lATK: ldy[22,3] - ldz[22,3] = 0
 lATK: lDY[22,3] - lDZ[22,3] + lstk[49,3] >= 0
 lATK: - lDY[22,3] + lDZ[22,3] + lstk[49,3] >= 0
 lATK: lDY[22,3] + lDZ[22,3] - ldz[22,3] - lstk[49,3] >= 0
 lATK: lDY[22,3] - ldz[22,3] >= 0
 lATK: ldy[22,4] - ldz[22,4] = 0
 lATK: lDY[22,4] - lDZ[22,4] + lstk[49,4] >= 0
 lATK: - lDY[22,4] + lDZ[22,4] + lstk[49,4] >= 0
 lATK: lDY[22,4] + lDZ[22,4] - ldz[22,4] - lstk[49,4] >= 0
 lATK: lDY[22,4] - ldz[22,4] >= 0
 lATK: ldy[22,5] - ldz[22,5] = 0
 lATK: lDY[22,5] - lDZ[22,5] + lstk[49,5] >= 0
 lATK: - lDY[22,5] + lDZ[22,5] + lstk[49,5] >= 0
 lATK: lDY[22,5] + lDZ[22,5] - ldz[22,5] - lstk[49,5] >= 0
 lATK: lDY[22,5] - ldz[22,5] >= 0
 lATK: ldy[22,6] - ldz[22,6] = 0
 lATK: lDY[22,6] - lDZ[22,6] + lstk[49,6] >= 0
 lATK: - lDY[22,6] + lDZ[22,6] + lstk[49,6] >= 0
 lATK: lDY[22,6] + lDZ[22,6] - ldz[22,6] - lstk[49,6] >= 0
 lATK: lDY[22,6] - ldz[22,6] >= 0
 lATK: ldy[22,7] - ldz[22,7] = 0
 lATK: lDY[22,7] - lDZ[22,7] + lstk[49,7] >= 0
 lATK: - lDY[22,7] + lDZ[22,7] + lstk[49,7] >= 0
 lATK: lDY[22,7] + lDZ[22,7] - ldz[22,7] - lstk[49,7] >= 0
 lATK: lDY[22,7] - ldz[22,7] >= 0
 lATK: lDY[22,8] - lDZ[22,8] = 0
 lATK: ldy[22,8] - ldz[22,8] = 0
 lATK: lDY[22,9] - lDZ[22,9] = 0
 lATK: ldy[22,9] - ldz[22,9] = 0
 lATK: lDY[22,10] - lDZ[22,10] = 0
 lATK: ldy[22,10] - ldz[22,10] = 0
 lATK: lDY[22,11] - lDZ[22,11] = 0
 lATK: ldy[22,11] - ldz[22,11] = 0
 lATK: lDY[22,12] - lDZ[22,12] = 0
 lATK: ldy[22,12] - ldz[22,12] = 0
 lATK: lDY[22,13] - lDZ[22,13] = 0
 lATK: ldy[22,13] - ldz[22,13] = 0
 lATK: lDY[22,14] - lDZ[22,14] = 0
 lATK: ldy[22,14] - ldz[22,14] = 0
 lATK: lDY[22,15] - lDZ[22,15] = 0
 lATK: ldy[22,15] - ldz[22,15] = 0
 lATK: ldy[23,0] - ldz[23,0] = 0
 lATK: lDY[23,0] - lDZ[23,0] + lstk[50,0] >= 0
 lATK: - lDY[23,0] + lDZ[23,0] + lstk[50,0] >= 0
 lATK: lDY[23,0] + lDZ[23,0] - ldz[23,0] - lstk[50,0] >= 0
 lATK: lDY[23,0] - ldz[23,0] >= 0
 lATK: ldy[23,1] - ldz[23,1] = 0
 lATK: lDY[23,1] - lDZ[23,1] + lstk[50,1] >= 0
 lATK: - lDY[23,1] + lDZ[23,1] + lstk[50,1] >= 0
 lATK: lDY[23,1] + lDZ[23,1] - ldz[23,1] - lstk[50,1] >= 0
 lATK: lDY[23,1] - ldz[23,1] >= 0
 lATK: ldy[23,2] - ldz[23,2] = 0
 lATK: lDY[23,2] - lDZ[23,2] + lstk[50,2] >= 0
 lATK: - lDY[23,2] + lDZ[23,2] + lstk[50,2] >= 0
 lATK: lDY[23,2] + lDZ[23,2] - ldz[23,2] - lstk[50,2] >= 0
 lATK: lDY[23,2] - ldz[23,2] >= 0
 lATK: ldy[23,3] - ldz[23,3] = 0
 lATK: lDY[23,3] - lDZ[23,3] + lstk[50,3] >= 0
 lATK: - lDY[23,3] + lDZ[23,3] + lstk[50,3] >= 0
 lATK: lDY[23,3] + lDZ[23,3] - ldz[23,3] - lstk[50,3] >= 0
 lATK: lDY[23,3] - ldz[23,3] >= 0
 lATK: ldy[23,4] - ldz[23,4] = 0
 lATK: lDY[23,4] - lDZ[23,4] + lstk[50,4] >= 0
 lATK: - lDY[23,4] + lDZ[23,4] + lstk[50,4] >= 0
 lATK: lDY[23,4] + lDZ[23,4] - ldz[23,4] - lstk[50,4] >= 0
 lATK: lDY[23,4] - ldz[23,4] >= 0
 lATK: ldy[23,5] - ldz[23,5] = 0
 lATK: lDY[23,5] - lDZ[23,5] + lstk[50,5] >= 0
 lATK: - lDY[23,5] + lDZ[23,5] + lstk[50,5] >= 0
 lATK: lDY[23,5] + lDZ[23,5] - ldz[23,5] - lstk[50,5] >= 0
 lATK: lDY[23,5] - ldz[23,5] >= 0
 lATK: ldy[23,6] - ldz[23,6] = 0
 lATK: lDY[23,6] - lDZ[23,6] + lstk[50,6] >= 0
 lATK: - lDY[23,6] + lDZ[23,6] + lstk[50,6] >= 0
 lATK: lDY[23,6] + lDZ[23,6] - ldz[23,6] - lstk[50,6] >= 0
 lATK: lDY[23,6] - ldz[23,6] >= 0
 lATK: ldy[23,7] - ldz[23,7] = 0
 lATK: lDY[23,7] - lDZ[23,7] + lstk[50,7] >= 0
 lATK: - lDY[23,7] + lDZ[23,7] + lstk[50,7] >= 0
 lATK: lDY[23,7] + lDZ[23,7] - ldz[23,7] - lstk[50,7] >= 0
 lATK: lDY[23,7] - ldz[23,7] >= 0
 lATK: lDY[23,8] - lDZ[23,8] = 0
 lATK: ldy[23,8] - ldz[23,8] = 0
 lATK: lDY[23,9] - lDZ[23,9] = 0
 lATK: ldy[23,9] - ldz[23,9] = 0
 lATK: lDY[23,10] - lDZ[23,10] = 0
 lATK: ldy[23,10] - ldz[23,10] = 0
 lATK: lDY[23,11] - lDZ[23,11] = 0
 lATK: ldy[23,11] - ldz[23,11] = 0
 lATK: lDY[23,12] - lDZ[23,12] = 0
 lATK: ldy[23,12] - ldz[23,12] = 0
 lATK: lDY[23,13] - lDZ[23,13] = 0
 lATK: ldy[23,13] - ldz[23,13] = 0
 lATK: lDY[23,14] - lDZ[23,14] = 0
 lATK: ldy[23,14] - ldz[23,14] = 0
 lATK: lDY[23,15] - lDZ[23,15] = 0
 lATK: ldy[23,15] - ldz[23,15] = 0
 lSR: lDZ[12,0] - lDW[12,0] = 0
 lSR: ldz[12,0] - ldw[12,0] = 0
 lSR: lDZ[12,1] - lDW[12,1] = 0
 lSR: ldz[12,1] - ldw[12,1] = 0
 lSR: lDZ[12,2] - lDW[12,2] = 0
 lSR: ldz[12,2] - ldw[12,2] = 0
 lSR: lDZ[12,3] - lDW[12,3] = 0
 lSR: ldz[12,3] - ldw[12,3] = 0
 lSR: lDZ[12,4] - lDW[12,5] = 0
 lSR: ldz[12,4] - ldw[12,5] = 0
 lSR: lDZ[12,5] - lDW[12,6] = 0
 lSR: ldz[12,5] - ldw[12,6] = 0
 lSR: lDZ[12,6] - lDW[12,7] = 0
 lSR: ldz[12,6] - ldw[12,7] = 0
 lSR: lDZ[12,7] - lDW[12,4] = 0
 lSR: ldz[12,7] - ldw[12,4] = 0
 lSR: lDZ[12,8] - lDW[12,10] = 0
 lSR: ldz[12,8] - ldw[12,10] = 0
 lSR: lDZ[12,9] - lDW[12,11] = 0
 lSR: ldz[12,9] - ldw[12,11] = 0
 lSR: lDZ[12,10] - lDW[12,8] = 0
 lSR: ldz[12,10] - ldw[12,8] = 0
 lSR: lDZ[12,11] - lDW[12,9] = 0
 lSR: ldz[12,11] - ldw[12,9] = 0
 lSR: lDZ[12,12] - lDW[12,15] = 0
 lSR: ldz[12,12] - ldw[12,15] = 0
 lSR: lDZ[12,13] - lDW[12,12] = 0
 lSR: ldz[12,13] - ldw[12,12] = 0
 lSR: lDZ[12,14] - lDW[12,13] = 0
 lSR: ldz[12,14] - ldw[12,13] = 0
 lSR: lDZ[12,15] - lDW[12,14] = 0
 lSR: ldz[12,15] - ldw[12,14] = 0
 lSR: lDZ[13,0] - lDW[13,0] = 0
 lSR: ldz[13,0] - ldw[13,0] = 0
 lSR: lDZ[13,1] - lDW[13,1] = 0
 lSR: ldz[13,1] - ldw[13,1] = 0
 lSR: lDZ[13,2] - lDW[13,2] = 0
 lSR: ldz[13,2] - ldw[13,2] = 0
 lSR: lDZ[13,3] - lDW[13,3] = 0
 lSR: ldz[13,3] - ldw[13,3] = 0
 lSR: lDZ[13,4] - lDW[13,5] = 0
 lSR: ldz[13,4] - ldw[13,5] = 0
 lSR: lDZ[13,5] - lDW[13,6] = 0
 lSR: ldz[13,5] - ldw[13,6] = 0
 lSR: lDZ[13,6] - lDW[13,7] = 0
 lSR: ldz[13,6] - ldw[13,7] = 0
 lSR: lDZ[13,7] - lDW[13,4] = 0
 lSR: ldz[13,7] - ldw[13,4] = 0
 lSR: lDZ[13,8] - lDW[13,10] = 0
 lSR: ldz[13,8] - ldw[13,10] = 0
 lSR: lDZ[13,9] - lDW[13,11] = 0
 lSR: ldz[13,9] - ldw[13,11] = 0
 lSR: lDZ[13,10] - lDW[13,8] = 0
 lSR: ldz[13,10] - ldw[13,8] = 0
 lSR: lDZ[13,11] - lDW[13,9] = 0
 lSR: ldz[13,11] - ldw[13,9] = 0
 lSR: lDZ[13,12] - lDW[13,15] = 0
 lSR: ldz[13,12] - ldw[13,15] = 0
 lSR: lDZ[13,13] - lDW[13,12] = 0
 lSR: ldz[13,13] - ldw[13,12] = 0
 lSR: lDZ[13,14] - lDW[13,13] = 0
 lSR: ldz[13,14] - ldw[13,13] = 0
 lSR: lDZ[13,15] - lDW[13,14] = 0
 lSR: ldz[13,15] - ldw[13,14] = 0
 lSR: lDZ[14,0] - lDW[14,0] = 0
 lSR: ldz[14,0] - ldw[14,0] = 0
 lSR: lDZ[14,1] - lDW[14,1] = 0
 lSR: ldz[14,1] - ldw[14,1] = 0
 lSR: lDZ[14,2] - lDW[14,2] = 0
 lSR: ldz[14,2] - ldw[14,2] = 0
 lSR: lDZ[14,3] - lDW[14,3] = 0
 lSR: ldz[14,3] - ldw[14,3] = 0
 lSR: lDZ[14,4] - lDW[14,5] = 0
 lSR: ldz[14,4] - ldw[14,5] = 0
 lSR: lDZ[14,5] - lDW[14,6] = 0
 lSR: ldz[14,5] - ldw[14,6] = 0
 lSR: lDZ[14,6] - lDW[14,7] = 0
 lSR: ldz[14,6] - ldw[14,7] = 0
 lSR: lDZ[14,7] - lDW[14,4] = 0
 lSR: ldz[14,7] - ldw[14,4] = 0
 lSR: lDZ[14,8] - lDW[14,10] = 0
 lSR: ldz[14,8] - ldw[14,10] = 0
 lSR: lDZ[14,9] - lDW[14,11] = 0
 lSR: ldz[14,9] - ldw[14,11] = 0
 lSR: lDZ[14,10] - lDW[14,8] = 0
 lSR: ldz[14,10] - ldw[14,8] = 0
 lSR: lDZ[14,11] - lDW[14,9] = 0
 lSR: ldz[14,11] - ldw[14,9] = 0
 lSR: lDZ[14,12] - lDW[14,15] = 0
 lSR: ldz[14,12] - ldw[14,15] = 0
 lSR: lDZ[14,13] - lDW[14,12] = 0
 lSR: ldz[14,13] - ldw[14,12] = 0
 lSR: lDZ[14,14] - lDW[14,13] = 0
 lSR: ldz[14,14] - ldw[14,13] = 0
 lSR: lDZ[14,15] - lDW[14,14] = 0
 lSR: ldz[14,15] - ldw[14,14] = 0
 lSR: lDZ[15,0] - lDW[15,0] = 0
 lSR: ldz[15,0] - ldw[15,0] = 0
 lSR: lDZ[15,1] - lDW[15,1] = 0
 lSR: ldz[15,1] - ldw[15,1] = 0
 lSR: lDZ[15,2] - lDW[15,2] = 0
 lSR: ldz[15,2] - ldw[15,2] = 0
 lSR: lDZ[15,3] - lDW[15,3] = 0
 lSR: ldz[15,3] - ldw[15,3] = 0
 lSR: lDZ[15,4] - lDW[15,5] = 0
 lSR: ldz[15,4] - ldw[15,5] = 0
 lSR: lDZ[15,5] - lDW[15,6] = 0
 lSR: ldz[15,5] - ldw[15,6] = 0
 lSR: lDZ[15,6] - lDW[15,7] = 0
 lSR: ldz[15,6] - ldw[15,7] = 0
 lSR: lDZ[15,7] - lDW[15,4] = 0
 lSR: ldz[15,7] - ldw[15,4] = 0
 lSR: lDZ[15,8] - lDW[15,10] = 0
 lSR: ldz[15,8] - ldw[15,10] = 0
 lSR: lDZ[15,9] - lDW[15,11] = 0
 lSR: ldz[15,9] - ldw[15,11] = 0
 lSR: lDZ[15,10] - lDW[15,8] = 0
 lSR: ldz[15,10] - ldw[15,8] = 0
 lSR: lDZ[15,11] - lDW[15,9] = 0
 lSR: ldz[15,11] - ldw[15,9] = 0
 lSR: lDZ[15,12] - lDW[15,15] = 0
 lSR: ldz[15,12] - ldw[15,15] = 0
 lSR: lDZ[15,13] - lDW[15,12] = 0
 lSR: ldz[15,13] - ldw[15,12] = 0
 lSR: lDZ[15,14] - lDW[15,13] = 0
 lSR: ldz[15,14] - ldw[15,13] = 0
 lSR: lDZ[15,15] - lDW[15,14] = 0
 lSR: ldz[15,15] - ldw[15,14] = 0
 lSR: lDZ[16,0] - lDW[16,0] = 0
 lSR: ldz[16,0] - ldw[16,0] = 0
 lSR: lDZ[16,1] - lDW[16,1] = 0
 lSR: ldz[16,1] - ldw[16,1] = 0
 lSR: lDZ[16,2] - lDW[16,2] = 0
 lSR: ldz[16,2] - ldw[16,2] = 0
 lSR: lDZ[16,3] - lDW[16,3] = 0
 lSR: ldz[16,3] - ldw[16,3] = 0
 lSR: lDZ[16,4] - lDW[16,5] = 0
 lSR: ldz[16,4] - ldw[16,5] = 0
 lSR: lDZ[16,5] - lDW[16,6] = 0
 lSR: ldz[16,5] - ldw[16,6] = 0
 lSR: lDZ[16,6] - lDW[16,7] = 0
 lSR: ldz[16,6] - ldw[16,7] = 0
 lSR: lDZ[16,7] - lDW[16,4] = 0
 lSR: ldz[16,7] - ldw[16,4] = 0
 lSR: lDZ[16,8] - lDW[16,10] = 0
 lSR: ldz[16,8] - ldw[16,10] = 0
 lSR: lDZ[16,9] - lDW[16,11] = 0
 lSR: ldz[16,9] - ldw[16,11] = 0
 lSR: lDZ[16,10] - lDW[16,8] = 0
 lSR: ldz[16,10] - ldw[16,8] = 0
 lSR: lDZ[16,11] - lDW[16,9] = 0
 lSR: ldz[16,11] - ldw[16,9] = 0
 lSR: lDZ[16,12] - lDW[16,15] = 0
 lSR: ldz[16,12] - ldw[16,15] = 0
 lSR: lDZ[16,13] - lDW[16,12] = 0
 lSR: ldz[16,13] - ldw[16,12] = 0
 lSR: lDZ[16,14] - lDW[16,13] = 0
 lSR: ldz[16,14] - ldw[16,13] = 0
 lSR: lDZ[16,15] - lDW[16,14] = 0
 lSR: ldz[16,15] - ldw[16,14] = 0
 lSR: lDZ[17,0] - lDW[17,0] = 0
 lSR: ldz[17,0] - ldw[17,0] = 0
 lSR: lDZ[17,1] - lDW[17,1] = 0
 lSR: ldz[17,1] - ldw[17,1] = 0
 lSR: lDZ[17,2] - lDW[17,2] = 0
 lSR: ldz[17,2] - ldw[17,2] = 0
 lSR: lDZ[17,3] - lDW[17,3] = 0
 lSR: ldz[17,3] - ldw[17,3] = 0
 lSR: lDZ[17,4] - lDW[17,5] = 0
 lSR: ldz[17,4] - ldw[17,5] = 0
 lSR: lDZ[17,5] - lDW[17,6] = 0
 lSR: ldz[17,5] - ldw[17,6] = 0
 lSR: lDZ[17,6] - lDW[17,7] = 0
 lSR: ldz[17,6] - ldw[17,7] = 0
 lSR: lDZ[17,7] - lDW[17,4] = 0
 lSR: ldz[17,7] - ldw[17,4] = 0
 lSR: lDZ[17,8] - lDW[17,10] = 0
 lSR: ldz[17,8] - ldw[17,10] = 0
 lSR: lDZ[17,9] - lDW[17,11] = 0
 lSR: ldz[17,9] - ldw[17,11] = 0
 lSR: lDZ[17,10] - lDW[17,8] = 0
 lSR: ldz[17,10] - ldw[17,8] = 0
 lSR: lDZ[17,11] - lDW[17,9] = 0
 lSR: ldz[17,11] - ldw[17,9] = 0
 lSR: lDZ[17,12] - lDW[17,15] = 0
 lSR: ldz[17,12] - ldw[17,15] = 0
 lSR: lDZ[17,13] - lDW[17,12] = 0
 lSR: ldz[17,13] - ldw[17,12] = 0
 lSR: lDZ[17,14] - lDW[17,13] = 0
 lSR: ldz[17,14] - ldw[17,13] = 0
 lSR: lDZ[17,15] - lDW[17,14] = 0
 lSR: ldz[17,15] - ldw[17,14] = 0
 lSR: lDZ[18,0] - lDW[18,0] = 0
 lSR: ldz[18,0] - ldw[18,0] = 0
 lSR: lDZ[18,1] - lDW[18,1] = 0
 lSR: ldz[18,1] - ldw[18,1] = 0
 lSR: lDZ[18,2] - lDW[18,2] = 0
 lSR: ldz[18,2] - ldw[18,2] = 0
 lSR: lDZ[18,3] - lDW[18,3] = 0
 lSR: ldz[18,3] - ldw[18,3] = 0
 lSR: lDZ[18,4] - lDW[18,5] = 0
 lSR: ldz[18,4] - ldw[18,5] = 0
 lSR: lDZ[18,5] - lDW[18,6] = 0
 lSR: ldz[18,5] - ldw[18,6] = 0
 lSR: lDZ[18,6] - lDW[18,7] = 0
 lSR: ldz[18,6] - ldw[18,7] = 0
 lSR: lDZ[18,7] - lDW[18,4] = 0
 lSR: ldz[18,7] - ldw[18,4] = 0
 lSR: lDZ[18,8] - lDW[18,10] = 0
 lSR: ldz[18,8] - ldw[18,10] = 0
 lSR: lDZ[18,9] - lDW[18,11] = 0
 lSR: ldz[18,9] - ldw[18,11] = 0
 lSR: lDZ[18,10] - lDW[18,8] = 0
 lSR: ldz[18,10] - ldw[18,8] = 0
 lSR: lDZ[18,11] - lDW[18,9] = 0
 lSR: ldz[18,11] - ldw[18,9] = 0
 lSR: lDZ[18,12] - lDW[18,15] = 0
 lSR: ldz[18,12] - ldw[18,15] = 0
 lSR: lDZ[18,13] - lDW[18,12] = 0
 lSR: ldz[18,13] - ldw[18,12] = 0
 lSR: lDZ[18,14] - lDW[18,13] = 0
 lSR: ldz[18,14] - ldw[18,13] = 0
 lSR: lDZ[18,15] - lDW[18,14] = 0
 lSR: ldz[18,15] - ldw[18,14] = 0
 lSR: lDZ[19,0] - lDW[19,0] = 0
 lSR: ldz[19,0] - ldw[19,0] = 0
 lSR: lDZ[19,1] - lDW[19,1] = 0
 lSR: ldz[19,1] - ldw[19,1] = 0
 lSR: lDZ[19,2] - lDW[19,2] = 0
 lSR: ldz[19,2] - ldw[19,2] = 0
 lSR: lDZ[19,3] - lDW[19,3] = 0
 lSR: ldz[19,3] - ldw[19,3] = 0
 lSR: lDZ[19,4] - lDW[19,5] = 0
 lSR: ldz[19,4] - ldw[19,5] = 0
 lSR: lDZ[19,5] - lDW[19,6] = 0
 lSR: ldz[19,5] - ldw[19,6] = 0
 lSR: lDZ[19,6] - lDW[19,7] = 0
 lSR: ldz[19,6] - ldw[19,7] = 0
 lSR: lDZ[19,7] - lDW[19,4] = 0
 lSR: ldz[19,7] - ldw[19,4] = 0
 lSR: lDZ[19,8] - lDW[19,10] = 0
 lSR: ldz[19,8] - ldw[19,10] = 0
 lSR: lDZ[19,9] - lDW[19,11] = 0
 lSR: ldz[19,9] - ldw[19,11] = 0
 lSR: lDZ[19,10] - lDW[19,8] = 0
 lSR: ldz[19,10] - ldw[19,8] = 0
 lSR: lDZ[19,11] - lDW[19,9] = 0
 lSR: ldz[19,11] - ldw[19,9] = 0
 lSR: lDZ[19,12] - lDW[19,15] = 0
 lSR: ldz[19,12] - ldw[19,15] = 0
 lSR: lDZ[19,13] - lDW[19,12] = 0
 lSR: ldz[19,13] - ldw[19,12] = 0
 lSR: lDZ[19,14] - lDW[19,13] = 0
 lSR: ldz[19,14] - ldw[19,13] = 0
 lSR: lDZ[19,15] - lDW[19,14] = 0
 lSR: ldz[19,15] - ldw[19,14] = 0
 lSR: lDZ[20,0] - lDW[20,0] = 0
 lSR: ldz[20,0] - ldw[20,0] = 0
 lSR: lDZ[20,1] - lDW[20,1] = 0
 lSR: ldz[20,1] - ldw[20,1] = 0
 lSR: lDZ[20,2] - lDW[20,2] = 0
 lSR: ldz[20,2] - ldw[20,2] = 0
 lSR: lDZ[20,3] - lDW[20,3] = 0
 lSR: ldz[20,3] - ldw[20,3] = 0
 lSR: lDZ[20,4] - lDW[20,5] = 0
 lSR: ldz[20,4] - ldw[20,5] = 0
 lSR: lDZ[20,5] - lDW[20,6] = 0
 lSR: ldz[20,5] - ldw[20,6] = 0
 lSR: lDZ[20,6] - lDW[20,7] = 0
 lSR: ldz[20,6] - ldw[20,7] = 0
 lSR: lDZ[20,7] - lDW[20,4] = 0
 lSR: ldz[20,7] - ldw[20,4] = 0
 lSR: lDZ[20,8] - lDW[20,10] = 0
 lSR: ldz[20,8] - ldw[20,10] = 0
 lSR: lDZ[20,9] - lDW[20,11] = 0
 lSR: ldz[20,9] - ldw[20,11] = 0
 lSR: lDZ[20,10] - lDW[20,8] = 0
 lSR: ldz[20,10] - ldw[20,8] = 0
 lSR: lDZ[20,11] - lDW[20,9] = 0
 lSR: ldz[20,11] - ldw[20,9] = 0
 lSR: lDZ[20,12] - lDW[20,15] = 0
 lSR: ldz[20,12] - ldw[20,15] = 0
 lSR: lDZ[20,13] - lDW[20,12] = 0
 lSR: ldz[20,13] - ldw[20,12] = 0
 lSR: lDZ[20,14] - lDW[20,13] = 0
 lSR: ldz[20,14] - ldw[20,13] = 0
 lSR: lDZ[20,15] - lDW[20,14] = 0
 lSR: ldz[20,15] - ldw[20,14] = 0
 lSR: lDZ[21,0] - lDW[21,0] = 0
 lSR: ldz[21,0] - ldw[21,0] = 0
 lSR: lDZ[21,1] - lDW[21,1] = 0
 lSR: ldz[21,1] - ldw[21,1] = 0
 lSR: lDZ[21,2] - lDW[21,2] = 0
 lSR: ldz[21,2] - ldw[21,2] = 0
 lSR: lDZ[21,3] - lDW[21,3] = 0
 lSR: ldz[21,3] - ldw[21,3] = 0
 lSR: lDZ[21,4] - lDW[21,5] = 0
 lSR: ldz[21,4] - ldw[21,5] = 0
 lSR: lDZ[21,5] - lDW[21,6] = 0
 lSR: ldz[21,5] - ldw[21,6] = 0
 lSR: lDZ[21,6] - lDW[21,7] = 0
 lSR: ldz[21,6] - ldw[21,7] = 0
 lSR: lDZ[21,7] - lDW[21,4] = 0
 lSR: ldz[21,7] - ldw[21,4] = 0
 lSR: lDZ[21,8] - lDW[21,10] = 0
 lSR: ldz[21,8] - ldw[21,10] = 0
 lSR: lDZ[21,9] - lDW[21,11] = 0
 lSR: ldz[21,9] - ldw[21,11] = 0
 lSR: lDZ[21,10] - lDW[21,8] = 0
 lSR: ldz[21,10] - ldw[21,8] = 0
 lSR: lDZ[21,11] - lDW[21,9] = 0
 lSR: ldz[21,11] - ldw[21,9] = 0
 lSR: lDZ[21,12] - lDW[21,15] = 0
 lSR: ldz[21,12] - ldw[21,15] = 0
 lSR: lDZ[21,13] - lDW[21,12] = 0
 lSR: ldz[21,13] - ldw[21,12] = 0
 lSR: lDZ[21,14] - lDW[21,13] = 0
 lSR: ldz[21,14] - ldw[21,13] = 0
 lSR: lDZ[21,15] - lDW[21,14] = 0
 lSR: ldz[21,15] - ldw[21,14] = 0
 lSR: lDZ[22,0] - lDW[22,0] = 0
 lSR: ldz[22,0] - ldw[22,0] = 0
 lSR: lDZ[22,1] - lDW[22,1] = 0
 lSR: ldz[22,1] - ldw[22,1] = 0
 lSR: lDZ[22,2] - lDW[22,2] = 0
 lSR: ldz[22,2] - ldw[22,2] = 0
 lSR: lDZ[22,3] - lDW[22,3] = 0
 lSR: ldz[22,3] - ldw[22,3] = 0
 lSR: lDZ[22,4] - lDW[22,5] = 0
 lSR: ldz[22,4] - ldw[22,5] = 0
 lSR: lDZ[22,5] - lDW[22,6] = 0
 lSR: ldz[22,5] - ldw[22,6] = 0
 lSR: lDZ[22,6] - lDW[22,7] = 0
 lSR: ldz[22,6] - ldw[22,7] = 0
 lSR: lDZ[22,7] - lDW[22,4] = 0
 lSR: ldz[22,7] - ldw[22,4] = 0
 lSR: lDZ[22,8] - lDW[22,10] = 0
 lSR: ldz[22,8] - ldw[22,10] = 0
 lSR: lDZ[22,9] - lDW[22,11] = 0
 lSR: ldz[22,9] - ldw[22,11] = 0
 lSR: lDZ[22,10] - lDW[22,8] = 0
 lSR: ldz[22,10] - ldw[22,8] = 0
 lSR: lDZ[22,11] - lDW[22,9] = 0
 lSR: ldz[22,11] - ldw[22,9] = 0
 lSR: lDZ[22,12] - lDW[22,15] = 0
 lSR: ldz[22,12] - ldw[22,15] = 0
 lSR: lDZ[22,13] - lDW[22,12] = 0
 lSR: ldz[22,13] - ldw[22,12] = 0
 lSR: lDZ[22,14] - lDW[22,13] = 0
 lSR: ldz[22,14] - ldw[22,13] = 0
 lSR: lDZ[22,15] - lDW[22,14] = 0
 lSR: ldz[22,15] - ldw[22,14] = 0
 lSR: lDZ[23,0] - lDW[23,0] = 0
 lSR: ldz[23,0] - ldw[23,0] = 0
 lSR: lDZ[23,1] - lDW[23,1] = 0
 lSR: ldz[23,1] - ldw[23,1] = 0
 lSR: lDZ[23,2] - lDW[23,2] = 0
 lSR: ldz[23,2] - ldw[23,2] = 0
 lSR: lDZ[23,3] - lDW[23,3] = 0
 lSR: ldz[23,3] - ldw[23,3] = 0
 lSR: lDZ[23,4] - lDW[23,5] = 0
 lSR: ldz[23,4] - ldw[23,5] = 0
 lSR: lDZ[23,5] - lDW[23,6] = 0
 lSR: ldz[23,5] - ldw[23,6] = 0
 lSR: lDZ[23,6] - lDW[23,7] = 0
 lSR: ldz[23,6] - ldw[23,7] = 0
 lSR: lDZ[23,7] - lDW[23,4] = 0
 lSR: ldz[23,7] - ldw[23,4] = 0
 lSR: lDZ[23,8] - lDW[23,10] = 0
 lSR: ldz[23,8] - ldw[23,10] = 0
 lSR: lDZ[23,9] - lDW[23,11] = 0
 lSR: ldz[23,9] - ldw[23,11] = 0
 lSR: lDZ[23,10] - lDW[23,8] = 0
 lSR: ldz[23,10] - ldw[23,8] = 0
 lSR: lDZ[23,11] - lDW[23,9] = 0
 lSR: ldz[23,11] - ldw[23,9] = 0
 lSR: lDZ[23,12] - lDW[23,15] = 0
 lSR: ldz[23,12] - ldw[23,15] = 0
 lSR: lDZ[23,13] - lDW[23,12] = 0
 lSR: ldz[23,13] - ldw[23,12] = 0
 lSR: lDZ[23,14] - lDW[23,13] = 0
 lSR: ldz[23,14] - ldw[23,13] = 0
 lSR: lDZ[23,15] - lDW[23,14] = 0
 lSR: ldz[23,15] - ldw[23,14] = 0
 ldMC0: - ldx[13,4] + ldw[12,0] = 0
 ldMC0: - lDX[13,4] + lDW[12,0] = 0
 ldMC1: ldw[12,4] - ldw[12,8] >= 0
 ldMC1: - ldx[13,8] + ldw[12,4] >= 0
 ldMC1: - ldx[13,8] + ldw[12,4] - ldw[12,8] <= 0
 ldMC1: lDX[13,8] - lDW[12,4] + lDW[12,8] >= 0
 ldMC1: lDX[13,8] + lDW[12,4] - lDW[12,8] >= 0
 ldMC1: - lDX[13,8] + lDW[12,4] + lDW[12,8] >= 0
 ldMC1: - lDX[13,8] + lDW[12,4] - lDW[12,8] <= 0
 ldMC2: - ldx[13,4] + ldw[12,8] >= 0
 ldMC2: - ldx[13,12] + ldw[12,8] >= 0
 ldMC2: - ldx[13,4] - ldx[13,12] + ldw[12,8] <= 0
 ldMC2: lDX[13,4] + lDX[13,12] - lDW[12,8] >= 0
 ldMC2: - lDX[13,4] + lDX[13,12] + lDW[12,8] >= 0
 ldMC2: lDX[13,4] - lDX[13,12] + lDW[12,8] >= 0
 ldMC2: - lDX[13,4] - lDX[13,12] + lDW[12,8] <= 0
 ldMC3: - ldx[13,0] + ldw[12,12] >= 0
 ldMC3: - ldx[13,12] + ldw[12,12] >= 0
 ldMC3: - ldx[13,0] - ldx[13,12] + ldw[12,12] <= 0
 ldMC3: lDX[13,0] + lDX[13,12] - lDW[12,12] >= 0
 ldMC3: - lDX[13,0] + lDX[13,12] + lDW[12,12] >= 0
 ldMC3: lDX[13,0] - lDX[13,12] + lDW[12,12] >= 0
 ldMC3: - lDX[13,0] - lDX[13,12] + lDW[12,12] <= 0
 ldMC0: - ldx[13,5] + ldw[12,1] = 0
 ldMC0: - lDX[13,5] + lDW[12,1] = 0
 ldMC1: ldw[12,5] - ldw[12,9] >= 0
 ldMC1: - ldx[13,9] + ldw[12,5] >= 0
 ldMC1: - ldx[13,9] + ldw[12,5] - ldw[12,9] <= 0
 ldMC1: lDX[13,9] - lDW[12,5] + lDW[12,9] >= 0
 ldMC1: lDX[13,9] + lDW[12,5] - lDW[12,9] >= 0
 ldMC1: - lDX[13,9] + lDW[12,5] + lDW[12,9] >= 0
 ldMC1: - lDX[13,9] + lDW[12,5] - lDW[12,9] <= 0
 ldMC2: - ldx[13,5] + ldw[12,9] >= 0
 ldMC2: - ldx[13,13] + ldw[12,9] >= 0
 ldMC2: - ldx[13,5] - ldx[13,13] + ldw[12,9] <= 0
 ldMC2: lDX[13,5] + lDX[13,13] - lDW[12,9] >= 0
 ldMC2: - lDX[13,5] + lDX[13,13] + lDW[12,9] >= 0
 ldMC2: lDX[13,5] - lDX[13,13] + lDW[12,9] >= 0
 ldMC2: - lDX[13,5] - lDX[13,13] + lDW[12,9] <= 0
 ldMC3: - ldx[13,1] + ldw[12,13] >= 0
 ldMC3: - ldx[13,13] + ldw[12,13] >= 0
 ldMC3: - ldx[13,1] - ldx[13,13] + ldw[12,13] <= 0
 ldMC3: lDX[13,1] + lDX[13,13] - lDW[12,13] >= 0
 ldMC3: - lDX[13,1] + lDX[13,13] + lDW[12,13] >= 0
 ldMC3: lDX[13,1] - lDX[13,13] + lDW[12,13] >= 0
 ldMC3: - lDX[13,1] - lDX[13,13] + lDW[12,13] <= 0
 ldMC0: - ldx[13,6] + ldw[12,2] = 0
 ldMC0: - lDX[13,6] + lDW[12,2] = 0
 ldMC1: ldw[12,6] - ldw[12,10] >= 0
 ldMC1: - ldx[13,10] + ldw[12,6] >= 0
 ldMC1: - ldx[13,10] + ldw[12,6] - ldw[12,10] <= 0
 ldMC1: lDX[13,10] - lDW[12,6] + lDW[12,10] >= 0
 ldMC1: lDX[13,10] + lDW[12,6] - lDW[12,10] >= 0
 ldMC1: - lDX[13,10] + lDW[12,6] + lDW[12,10] >= 0
 ldMC1: - lDX[13,10] + lDW[12,6] - lDW[12,10] <= 0
 ldMC2: - ldx[13,6] + ldw[12,10] >= 0
 ldMC2: - ldx[13,14] + ldw[12,10] >= 0
 ldMC2: - ldx[13,6] - ldx[13,14] + ldw[12,10] <= 0
 ldMC2: lDX[13,6] + lDX[13,14] - lDW[12,10] >= 0
 ldMC2: - lDX[13,6] + lDX[13,14] + lDW[12,10] >= 0
 ldMC2: lDX[13,6] - lDX[13,14] + lDW[12,10] >= 0
 ldMC2: - lDX[13,6] - lDX[13,14] + lDW[12,10] <= 0
 ldMC3: - ldx[13,2] + ldw[12,14] >= 0
 ldMC3: - ldx[13,14] + ldw[12,14] >= 0
 ldMC3: - ldx[13,2] - ldx[13,14] + ldw[12,14] <= 0
 ldMC3: lDX[13,2] + lDX[13,14] - lDW[12,14] >= 0
 ldMC3: - lDX[13,2] + lDX[13,14] + lDW[12,14] >= 0
 ldMC3: lDX[13,2] - lDX[13,14] + lDW[12,14] >= 0
 ldMC3: - lDX[13,2] - lDX[13,14] + lDW[12,14] <= 0
 ldMC0: - ldx[13,7] + ldw[12,3] = 0
 ldMC0: - lDX[13,7] + lDW[12,3] = 0
 ldMC1: ldw[12,7] - ldw[12,11] >= 0
 ldMC1: - ldx[13,11] + ldw[12,7] >= 0
 ldMC1: - ldx[13,11] + ldw[12,7] - ldw[12,11] <= 0
 ldMC1: lDX[13,11] - lDW[12,7] + lDW[12,11] >= 0
 ldMC1: lDX[13,11] + lDW[12,7] - lDW[12,11] >= 0
 ldMC1: - lDX[13,11] + lDW[12,7] + lDW[12,11] >= 0
 ldMC1: - lDX[13,11] + lDW[12,7] - lDW[12,11] <= 0
 ldMC2: - ldx[13,7] + ldw[12,11] >= 0
 ldMC2: - ldx[13,15] + ldw[12,11] >= 0
 ldMC2: - ldx[13,7] - ldx[13,15] + ldw[12,11] <= 0
 ldMC2: lDX[13,7] + lDX[13,15] - lDW[12,11] >= 0
 ldMC2: - lDX[13,7] + lDX[13,15] + lDW[12,11] >= 0
 ldMC2: lDX[13,7] - lDX[13,15] + lDW[12,11] >= 0
 ldMC2: - lDX[13,7] - lDX[13,15] + lDW[12,11] <= 0
 ldMC3: - ldx[13,3] + ldw[12,15] >= 0
 ldMC3: - ldx[13,15] + ldw[12,15] >= 0
 ldMC3: - ldx[13,3] - ldx[13,15] + ldw[12,15] <= 0
 ldMC3: lDX[13,3] + lDX[13,15] - lDW[12,15] >= 0
 ldMC3: - lDX[13,3] + lDX[13,15] + lDW[12,15] >= 0
 ldMC3: lDX[13,3] - lDX[13,15] + lDW[12,15] >= 0
 ldMC3: - lDX[13,3] - lDX[13,15] + lDW[12,15] <= 0
 ldMC0: - ldx[14,4] + ldw[13,0] = 0
 ldMC0: - lDX[14,4] + lDW[13,0] = 0
 ldMC1: ldw[13,4] - ldw[13,8] >= 0
 ldMC1: - ldx[14,8] + ldw[13,4] >= 0
 ldMC1: - ldx[14,8] + ldw[13,4] - ldw[13,8] <= 0
 ldMC1: lDX[14,8] - lDW[13,4] + lDW[13,8] >= 0
 ldMC1: lDX[14,8] + lDW[13,4] - lDW[13,8] >= 0
 ldMC1: - lDX[14,8] + lDW[13,4] + lDW[13,8] >= 0
 ldMC1: - lDX[14,8] + lDW[13,4] - lDW[13,8] <= 0
 ldMC2: - ldx[14,4] + ldw[13,8] >= 0
 ldMC2: - ldx[14,12] + ldw[13,8] >= 0
 ldMC2: - ldx[14,4] - ldx[14,12] + ldw[13,8] <= 0
 ldMC2: lDX[14,4] + lDX[14,12] - lDW[13,8] >= 0
 ldMC2: - lDX[14,4] + lDX[14,12] + lDW[13,8] >= 0
 ldMC2: lDX[14,4] - lDX[14,12] + lDW[13,8] >= 0
 ldMC2: - lDX[14,4] - lDX[14,12] + lDW[13,8] <= 0
 ldMC3: - ldx[14,0] + ldw[13,12] >= 0
 ldMC3: - ldx[14,12] + ldw[13,12] >= 0
 ldMC3: - ldx[14,0] - ldx[14,12] + ldw[13,12] <= 0
 ldMC3: lDX[14,0] + lDX[14,12] - lDW[13,12] >= 0
 ldMC3: - lDX[14,0] + lDX[14,12] + lDW[13,12] >= 0
 ldMC3: lDX[14,0] - lDX[14,12] + lDW[13,12] >= 0
 ldMC3: - lDX[14,0] - lDX[14,12] + lDW[13,12] <= 0
 ldMC0: - ldx[14,5] + ldw[13,1] = 0
 ldMC0: - lDX[14,5] + lDW[13,1] = 0
 ldMC1: ldw[13,5] - ldw[13,9] >= 0
 ldMC1: - ldx[14,9] + ldw[13,5] >= 0
 ldMC1: - ldx[14,9] + ldw[13,5] - ldw[13,9] <= 0
 ldMC1: lDX[14,9] - lDW[13,5] + lDW[13,9] >= 0
 ldMC1: lDX[14,9] + lDW[13,5] - lDW[13,9] >= 0
 ldMC1: - lDX[14,9] + lDW[13,5] + lDW[13,9] >= 0
 ldMC1: - lDX[14,9] + lDW[13,5] - lDW[13,9] <= 0
 ldMC2: - ldx[14,5] + ldw[13,9] >= 0
 ldMC2: - ldx[14,13] + ldw[13,9] >= 0
 ldMC2: - ldx[14,5] - ldx[14,13] + ldw[13,9] <= 0
 ldMC2: lDX[14,5] + lDX[14,13] - lDW[13,9] >= 0
 ldMC2: - lDX[14,5] + lDX[14,13] + lDW[13,9] >= 0
 ldMC2: lDX[14,5] - lDX[14,13] + lDW[13,9] >= 0
 ldMC2: - lDX[14,5] - lDX[14,13] + lDW[13,9] <= 0
 ldMC3: - ldx[14,1] + ldw[13,13] >= 0
 ldMC3: - ldx[14,13] + ldw[13,13] >= 0
 ldMC3: - ldx[14,1] - ldx[14,13] + ldw[13,13] <= 0
 ldMC3: lDX[14,1] + lDX[14,13] - lDW[13,13] >= 0
 ldMC3: - lDX[14,1] + lDX[14,13] + lDW[13,13] >= 0
 ldMC3: lDX[14,1] - lDX[14,13] + lDW[13,13] >= 0
 ldMC3: - lDX[14,1] - lDX[14,13] + lDW[13,13] <= 0
 ldMC0: - ldx[14,6] + ldw[13,2] = 0
 ldMC0: - lDX[14,6] + lDW[13,2] = 0
 ldMC1: ldw[13,6] - ldw[13,10] >= 0
 ldMC1: - ldx[14,10] + ldw[13,6] >= 0
 ldMC1: - ldx[14,10] + ldw[13,6] - ldw[13,10] <= 0
 ldMC1: lDX[14,10] - lDW[13,6] + lDW[13,10] >= 0
 ldMC1: lDX[14,10] + lDW[13,6] - lDW[13,10] >= 0
 ldMC1: - lDX[14,10] + lDW[13,6] + lDW[13,10] >= 0
 ldMC1: - lDX[14,10] + lDW[13,6] - lDW[13,10] <= 0
 ldMC2: - ldx[14,6] + ldw[13,10] >= 0
 ldMC2: - ldx[14,14] + ldw[13,10] >= 0
 ldMC2: - ldx[14,6] - ldx[14,14] + ldw[13,10] <= 0
 ldMC2: lDX[14,6] + lDX[14,14] - lDW[13,10] >= 0
 ldMC2: - lDX[14,6] + lDX[14,14] + lDW[13,10] >= 0
 ldMC2: lDX[14,6] - lDX[14,14] + lDW[13,10] >= 0
 ldMC2: - lDX[14,6] - lDX[14,14] + lDW[13,10] <= 0
 ldMC3: - ldx[14,2] + ldw[13,14] >= 0
 ldMC3: - ldx[14,14] + ldw[13,14] >= 0
 ldMC3: - ldx[14,2] - ldx[14,14] + ldw[13,14] <= 0
 ldMC3: lDX[14,2] + lDX[14,14] - lDW[13,14] >= 0
 ldMC3: - lDX[14,2] + lDX[14,14] + lDW[13,14] >= 0
 ldMC3: lDX[14,2] - lDX[14,14] + lDW[13,14] >= 0
 ldMC3: - lDX[14,2] - lDX[14,14] + lDW[13,14] <= 0
 ldMC0: - ldx[14,7] + ldw[13,3] = 0
 ldMC0: - lDX[14,7] + lDW[13,3] = 0
 ldMC1: ldw[13,7] - ldw[13,11] >= 0
 ldMC1: - ldx[14,11] + ldw[13,7] >= 0
 ldMC1: - ldx[14,11] + ldw[13,7] - ldw[13,11] <= 0
 ldMC1: lDX[14,11] - lDW[13,7] + lDW[13,11] >= 0
 ldMC1: lDX[14,11] + lDW[13,7] - lDW[13,11] >= 0
 ldMC1: - lDX[14,11] + lDW[13,7] + lDW[13,11] >= 0
 ldMC1: - lDX[14,11] + lDW[13,7] - lDW[13,11] <= 0
 ldMC2: - ldx[14,7] + ldw[13,11] >= 0
 ldMC2: - ldx[14,15] + ldw[13,11] >= 0
 ldMC2: - ldx[14,7] - ldx[14,15] + ldw[13,11] <= 0
 ldMC2: lDX[14,7] + lDX[14,15] - lDW[13,11] >= 0
 ldMC2: - lDX[14,7] + lDX[14,15] + lDW[13,11] >= 0
 ldMC2: lDX[14,7] - lDX[14,15] + lDW[13,11] >= 0
 ldMC2: - lDX[14,7] - lDX[14,15] + lDW[13,11] <= 0
 ldMC3: - ldx[14,3] + ldw[13,15] >= 0
 ldMC3: - ldx[14,15] + ldw[13,15] >= 0
 ldMC3: - ldx[14,3] - ldx[14,15] + ldw[13,15] <= 0
 ldMC3: lDX[14,3] + lDX[14,15] - lDW[13,15] >= 0
 ldMC3: - lDX[14,3] + lDX[14,15] + lDW[13,15] >= 0
 ldMC3: lDX[14,3] - lDX[14,15] + lDW[13,15] >= 0
 ldMC3: - lDX[14,3] - lDX[14,15] + lDW[13,15] <= 0
 ldMC0: - ldx[15,4] + ldw[14,0] = 0
 ldMC0: - lDX[15,4] + lDW[14,0] = 0
 ldMC1: ldw[14,4] - ldw[14,8] >= 0
 ldMC1: - ldx[15,8] + ldw[14,4] >= 0
 ldMC1: - ldx[15,8] + ldw[14,4] - ldw[14,8] <= 0
 ldMC1: lDX[15,8] - lDW[14,4] + lDW[14,8] >= 0
 ldMC1: lDX[15,8] + lDW[14,4] - lDW[14,8] >= 0
 ldMC1: - lDX[15,8] + lDW[14,4] + lDW[14,8] >= 0
 ldMC1: - lDX[15,8] + lDW[14,4] - lDW[14,8] <= 0
 ldMC2: - ldx[15,4] + ldw[14,8] >= 0
 ldMC2: - ldx[15,12] + ldw[14,8] >= 0
 ldMC2: - ldx[15,4] - ldx[15,12] + ldw[14,8] <= 0
 ldMC2: lDX[15,4] + lDX[15,12] - lDW[14,8] >= 0
 ldMC2: - lDX[15,4] + lDX[15,12] + lDW[14,8] >= 0
 ldMC2: lDX[15,4] - lDX[15,12] + lDW[14,8] >= 0
 ldMC2: - lDX[15,4] - lDX[15,12] + lDW[14,8] <= 0
 ldMC3: - ldx[15,0] + ldw[14,12] >= 0
 ldMC3: - ldx[15,12] + ldw[14,12] >= 0
 ldMC3: - ldx[15,0] - ldx[15,12] + ldw[14,12] <= 0
 ldMC3: lDX[15,0] + lDX[15,12] - lDW[14,12] >= 0
 ldMC3: - lDX[15,0] + lDX[15,12] + lDW[14,12] >= 0
 ldMC3: lDX[15,0] - lDX[15,12] + lDW[14,12] >= 0
 ldMC3: - lDX[15,0] - lDX[15,12] + lDW[14,12] <= 0
 ldMC0: - ldx[15,5] + ldw[14,1] = 0
 ldMC0: - lDX[15,5] + lDW[14,1] = 0
 ldMC1: ldw[14,5] - ldw[14,9] >= 0
 ldMC1: - ldx[15,9] + ldw[14,5] >= 0
 ldMC1: - ldx[15,9] + ldw[14,5] - ldw[14,9] <= 0
 ldMC1: lDX[15,9] - lDW[14,5] + lDW[14,9] >= 0
 ldMC1: lDX[15,9] + lDW[14,5] - lDW[14,9] >= 0
 ldMC1: - lDX[15,9] + lDW[14,5] + lDW[14,9] >= 0
 ldMC1: - lDX[15,9] + lDW[14,5] - lDW[14,9] <= 0
 ldMC2: - ldx[15,5] + ldw[14,9] >= 0
 ldMC2: - ldx[15,13] + ldw[14,9] >= 0
 ldMC2: - ldx[15,5] - ldx[15,13] + ldw[14,9] <= 0
 ldMC2: lDX[15,5] + lDX[15,13] - lDW[14,9] >= 0
 ldMC2: - lDX[15,5] + lDX[15,13] + lDW[14,9] >= 0
 ldMC2: lDX[15,5] - lDX[15,13] + lDW[14,9] >= 0
 ldMC2: - lDX[15,5] - lDX[15,13] + lDW[14,9] <= 0
 ldMC3: - ldx[15,1] + ldw[14,13] >= 0
 ldMC3: - ldx[15,13] + ldw[14,13] >= 0
 ldMC3: - ldx[15,1] - ldx[15,13] + ldw[14,13] <= 0
 ldMC3: lDX[15,1] + lDX[15,13] - lDW[14,13] >= 0
 ldMC3: - lDX[15,1] + lDX[15,13] + lDW[14,13] >= 0
 ldMC3: lDX[15,1] - lDX[15,13] + lDW[14,13] >= 0
 ldMC3: - lDX[15,1] - lDX[15,13] + lDW[14,13] <= 0
 ldMC0: - ldx[15,6] + ldw[14,2] = 0
 ldMC0: - lDX[15,6] + lDW[14,2] = 0
 ldMC1: ldw[14,6] - ldw[14,10] >= 0
 ldMC1: - ldx[15,10] + ldw[14,6] >= 0
 ldMC1: - ldx[15,10] + ldw[14,6] - ldw[14,10] <= 0
 ldMC1: lDX[15,10] - lDW[14,6] + lDW[14,10] >= 0
 ldMC1: lDX[15,10] + lDW[14,6] - lDW[14,10] >= 0
 ldMC1: - lDX[15,10] + lDW[14,6] + lDW[14,10] >= 0
 ldMC1: - lDX[15,10] + lDW[14,6] - lDW[14,10] <= 0
 ldMC2: - ldx[15,6] + ldw[14,10] >= 0
 ldMC2: - ldx[15,14] + ldw[14,10] >= 0
 ldMC2: - ldx[15,6] - ldx[15,14] + ldw[14,10] <= 0
 ldMC2: lDX[15,6] + lDX[15,14] - lDW[14,10] >= 0
 ldMC2: - lDX[15,6] + lDX[15,14] + lDW[14,10] >= 0
 ldMC2: lDX[15,6] - lDX[15,14] + lDW[14,10] >= 0
 ldMC2: - lDX[15,6] - lDX[15,14] + lDW[14,10] <= 0
 ldMC3: - ldx[15,2] + ldw[14,14] >= 0
 ldMC3: - ldx[15,14] + ldw[14,14] >= 0
 ldMC3: - ldx[15,2] - ldx[15,14] + ldw[14,14] <= 0
 ldMC3: lDX[15,2] + lDX[15,14] - lDW[14,14] >= 0
 ldMC3: - lDX[15,2] + lDX[15,14] + lDW[14,14] >= 0
 ldMC3: lDX[15,2] - lDX[15,14] + lDW[14,14] >= 0
 ldMC3: - lDX[15,2] - lDX[15,14] + lDW[14,14] <= 0
 ldMC0: - ldx[15,7] + ldw[14,3] = 0
 ldMC0: - lDX[15,7] + lDW[14,3] = 0
 ldMC1: ldw[14,7] - ldw[14,11] >= 0
 ldMC1: - ldx[15,11] + ldw[14,7] >= 0
 ldMC1: - ldx[15,11] + ldw[14,7] - ldw[14,11] <= 0
 ldMC1: lDX[15,11] - lDW[14,7] + lDW[14,11] >= 0
 ldMC1: lDX[15,11] + lDW[14,7] - lDW[14,11] >= 0
 ldMC1: - lDX[15,11] + lDW[14,7] + lDW[14,11] >= 0
 ldMC1: - lDX[15,11] + lDW[14,7] - lDW[14,11] <= 0
 ldMC2: - ldx[15,7] + ldw[14,11] >= 0
 ldMC2: - ldx[15,15] + ldw[14,11] >= 0
 ldMC2: - ldx[15,7] - ldx[15,15] + ldw[14,11] <= 0
 ldMC2: lDX[15,7] + lDX[15,15] - lDW[14,11] >= 0
 ldMC2: - lDX[15,7] + lDX[15,15] + lDW[14,11] >= 0
 ldMC2: lDX[15,7] - lDX[15,15] + lDW[14,11] >= 0
 ldMC2: - lDX[15,7] - lDX[15,15] + lDW[14,11] <= 0
 ldMC3: - ldx[15,3] + ldw[14,15] >= 0
 ldMC3: - ldx[15,15] + ldw[14,15] >= 0
 ldMC3: - ldx[15,3] - ldx[15,15] + ldw[14,15] <= 0
 ldMC3: lDX[15,3] + lDX[15,15] - lDW[14,15] >= 0
 ldMC3: - lDX[15,3] + lDX[15,15] + lDW[14,15] >= 0
 ldMC3: lDX[15,3] - lDX[15,15] + lDW[14,15] >= 0
 ldMC3: - lDX[15,3] - lDX[15,15] + lDW[14,15] <= 0
 ldMC0: - ldx[16,4] + ldw[15,0] = 0
 ldMC0: - lDX[16,4] + lDW[15,0] = 0
 ldMC1: ldw[15,4] - ldw[15,8] >= 0
 ldMC1: - ldx[16,8] + ldw[15,4] >= 0
 ldMC1: - ldx[16,8] + ldw[15,4] - ldw[15,8] <= 0
 ldMC1: lDX[16,8] - lDW[15,4] + lDW[15,8] >= 0
 ldMC1: lDX[16,8] + lDW[15,4] - lDW[15,8] >= 0
 ldMC1: - lDX[16,8] + lDW[15,4] + lDW[15,8] >= 0
 ldMC1: - lDX[16,8] + lDW[15,4] - lDW[15,8] <= 0
 ldMC2: - ldx[16,4] + ldw[15,8] >= 0
 ldMC2: - ldx[16,12] + ldw[15,8] >= 0
 ldMC2: - ldx[16,4] - ldx[16,12] + ldw[15,8] <= 0
 ldMC2: lDX[16,4] + lDX[16,12] - lDW[15,8] >= 0
 ldMC2: - lDX[16,4] + lDX[16,12] + lDW[15,8] >= 0
 ldMC2: lDX[16,4] - lDX[16,12] + lDW[15,8] >= 0
 ldMC2: - lDX[16,4] - lDX[16,12] + lDW[15,8] <= 0
 ldMC3: - ldx[16,0] + ldw[15,12] >= 0
 ldMC3: - ldx[16,12] + ldw[15,12] >= 0
 ldMC3: - ldx[16,0] - ldx[16,12] + ldw[15,12] <= 0
 ldMC3: lDX[16,0] + lDX[16,12] - lDW[15,12] >= 0
 ldMC3: - lDX[16,0] + lDX[16,12] + lDW[15,12] >= 0
 ldMC3: lDX[16,0] - lDX[16,12] + lDW[15,12] >= 0
 ldMC3: - lDX[16,0] - lDX[16,12] + lDW[15,12] <= 0
 ldMC0: - ldx[16,5] + ldw[15,1] = 0
 ldMC0: - lDX[16,5] + lDW[15,1] = 0
 ldMC1: ldw[15,5] - ldw[15,9] >= 0
 ldMC1: - ldx[16,9] + ldw[15,5] >= 0
 ldMC1: - ldx[16,9] + ldw[15,5] - ldw[15,9] <= 0
 ldMC1: lDX[16,9] - lDW[15,5] + lDW[15,9] >= 0
 ldMC1: lDX[16,9] + lDW[15,5] - lDW[15,9] >= 0
 ldMC1: - lDX[16,9] + lDW[15,5] + lDW[15,9] >= 0
 ldMC1: - lDX[16,9] + lDW[15,5] - lDW[15,9] <= 0
 ldMC2: - ldx[16,5] + ldw[15,9] >= 0
 ldMC2: - ldx[16,13] + ldw[15,9] >= 0
 ldMC2: - ldx[16,5] - ldx[16,13] + ldw[15,9] <= 0
 ldMC2: lDX[16,5] + lDX[16,13] - lDW[15,9] >= 0
 ldMC2: - lDX[16,5] + lDX[16,13] + lDW[15,9] >= 0
 ldMC2: lDX[16,5] - lDX[16,13] + lDW[15,9] >= 0
 ldMC2: - lDX[16,5] - lDX[16,13] + lDW[15,9] <= 0
 ldMC3: - ldx[16,1] + ldw[15,13] >= 0
 ldMC3: - ldx[16,13] + ldw[15,13] >= 0
 ldMC3: - ldx[16,1] - ldx[16,13] + ldw[15,13] <= 0
 ldMC3: lDX[16,1] + lDX[16,13] - lDW[15,13] >= 0
 ldMC3: - lDX[16,1] + lDX[16,13] + lDW[15,13] >= 0
 ldMC3: lDX[16,1] - lDX[16,13] + lDW[15,13] >= 0
 ldMC3: - lDX[16,1] - lDX[16,13] + lDW[15,13] <= 0
 ldMC0: - ldx[16,6] + ldw[15,2] = 0
 ldMC0: - lDX[16,6] + lDW[15,2] = 0
 ldMC1: ldw[15,6] - ldw[15,10] >= 0
 ldMC1: - ldx[16,10] + ldw[15,6] >= 0
 ldMC1: - ldx[16,10] + ldw[15,6] - ldw[15,10] <= 0
 ldMC1: lDX[16,10] - lDW[15,6] + lDW[15,10] >= 0
 ldMC1: lDX[16,10] + lDW[15,6] - lDW[15,10] >= 0
 ldMC1: - lDX[16,10] + lDW[15,6] + lDW[15,10] >= 0
 ldMC1: - lDX[16,10] + lDW[15,6] - lDW[15,10] <= 0
 ldMC2: - ldx[16,6] + ldw[15,10] >= 0
 ldMC2: - ldx[16,14] + ldw[15,10] >= 0
 ldMC2: - ldx[16,6] - ldx[16,14] + ldw[15,10] <= 0
 ldMC2: lDX[16,6] + lDX[16,14] - lDW[15,10] >= 0
 ldMC2: - lDX[16,6] + lDX[16,14] + lDW[15,10] >= 0
 ldMC2: lDX[16,6] - lDX[16,14] + lDW[15,10] >= 0
 ldMC2: - lDX[16,6] - lDX[16,14] + lDW[15,10] <= 0
 ldMC3: - ldx[16,2] + ldw[15,14] >= 0
 ldMC3: - ldx[16,14] + ldw[15,14] >= 0
 ldMC3: - ldx[16,2] - ldx[16,14] + ldw[15,14] <= 0
 ldMC3: lDX[16,2] + lDX[16,14] - lDW[15,14] >= 0
 ldMC3: - lDX[16,2] + lDX[16,14] + lDW[15,14] >= 0
 ldMC3: lDX[16,2] - lDX[16,14] + lDW[15,14] >= 0
 ldMC3: - lDX[16,2] - lDX[16,14] + lDW[15,14] <= 0
 ldMC0: - ldx[16,7] + ldw[15,3] = 0
 ldMC0: - lDX[16,7] + lDW[15,3] = 0
 ldMC1: ldw[15,7] - ldw[15,11] >= 0
 ldMC1: - ldx[16,11] + ldw[15,7] >= 0
 ldMC1: - ldx[16,11] + ldw[15,7] - ldw[15,11] <= 0
 ldMC1: lDX[16,11] - lDW[15,7] + lDW[15,11] >= 0
 ldMC1: lDX[16,11] + lDW[15,7] - lDW[15,11] >= 0
 ldMC1: - lDX[16,11] + lDW[15,7] + lDW[15,11] >= 0
 ldMC1: - lDX[16,11] + lDW[15,7] - lDW[15,11] <= 0
 ldMC2: - ldx[16,7] + ldw[15,11] >= 0
 ldMC2: - ldx[16,15] + ldw[15,11] >= 0
 ldMC2: - ldx[16,7] - ldx[16,15] + ldw[15,11] <= 0
 ldMC2: lDX[16,7] + lDX[16,15] - lDW[15,11] >= 0
 ldMC2: - lDX[16,7] + lDX[16,15] + lDW[15,11] >= 0
 ldMC2: lDX[16,7] - lDX[16,15] + lDW[15,11] >= 0
 ldMC2: - lDX[16,7] - lDX[16,15] + lDW[15,11] <= 0
 ldMC3: - ldx[16,3] + ldw[15,15] >= 0
 ldMC3: - ldx[16,15] + ldw[15,15] >= 0
 ldMC3: - ldx[16,3] - ldx[16,15] + ldw[15,15] <= 0
 ldMC3: lDX[16,3] + lDX[16,15] - lDW[15,15] >= 0
 ldMC3: - lDX[16,3] + lDX[16,15] + lDW[15,15] >= 0
 ldMC3: lDX[16,3] - lDX[16,15] + lDW[15,15] >= 0
 ldMC3: - lDX[16,3] - lDX[16,15] + lDW[15,15] <= 0
 ldMC0: - ldx[17,4] + ldw[16,0] = 0
 ldMC0: - lDX[17,4] + lDW[16,0] = 0
 ldMC1: ldw[16,4] - ldw[16,8] >= 0
 ldMC1: - ldx[17,8] + ldw[16,4] >= 0
 ldMC1: - ldx[17,8] + ldw[16,4] - ldw[16,8] <= 0
 ldMC1: lDX[17,8] - lDW[16,4] + lDW[16,8] >= 0
 ldMC1: lDX[17,8] + lDW[16,4] - lDW[16,8] >= 0
 ldMC1: - lDX[17,8] + lDW[16,4] + lDW[16,8] >= 0
 ldMC1: - lDX[17,8] + lDW[16,4] - lDW[16,8] <= 0
 ldMC2: - ldx[17,4] + ldw[16,8] >= 0
 ldMC2: - ldx[17,12] + ldw[16,8] >= 0
 ldMC2: - ldx[17,4] - ldx[17,12] + ldw[16,8] <= 0
 ldMC2: lDX[17,4] + lDX[17,12] - lDW[16,8] >= 0
 ldMC2: - lDX[17,4] + lDX[17,12] + lDW[16,8] >= 0
 ldMC2: lDX[17,4] - lDX[17,12] + lDW[16,8] >= 0
 ldMC2: - lDX[17,4] - lDX[17,12] + lDW[16,8] <= 0
 ldMC3: - ldx[17,0] + ldw[16,12] >= 0
 ldMC3: - ldx[17,12] + ldw[16,12] >= 0
 ldMC3: - ldx[17,0] - ldx[17,12] + ldw[16,12] <= 0
 ldMC3: lDX[17,0] + lDX[17,12] - lDW[16,12] >= 0
 ldMC3: - lDX[17,0] + lDX[17,12] + lDW[16,12] >= 0
 ldMC3: lDX[17,0] - lDX[17,12] + lDW[16,12] >= 0
 ldMC3: - lDX[17,0] - lDX[17,12] + lDW[16,12] <= 0
 ldMC0: - ldx[17,5] + ldw[16,1] = 0
 ldMC0: - lDX[17,5] + lDW[16,1] = 0
 ldMC1: ldw[16,5] - ldw[16,9] >= 0
 ldMC1: - ldx[17,9] + ldw[16,5] >= 0
 ldMC1: - ldx[17,9] + ldw[16,5] - ldw[16,9] <= 0
 ldMC1: lDX[17,9] - lDW[16,5] + lDW[16,9] >= 0
 ldMC1: lDX[17,9] + lDW[16,5] - lDW[16,9] >= 0
 ldMC1: - lDX[17,9] + lDW[16,5] + lDW[16,9] >= 0
 ldMC1: - lDX[17,9] + lDW[16,5] - lDW[16,9] <= 0
 ldMC2: - ldx[17,5] + ldw[16,9] >= 0
 ldMC2: - ldx[17,13] + ldw[16,9] >= 0
 ldMC2: - ldx[17,5] - ldx[17,13] + ldw[16,9] <= 0
 ldMC2: lDX[17,5] + lDX[17,13] - lDW[16,9] >= 0
 ldMC2: - lDX[17,5] + lDX[17,13] + lDW[16,9] >= 0
 ldMC2: lDX[17,5] - lDX[17,13] + lDW[16,9] >= 0
 ldMC2: - lDX[17,5] - lDX[17,13] + lDW[16,9] <= 0
 ldMC3: - ldx[17,1] + ldw[16,13] >= 0
 ldMC3: - ldx[17,13] + ldw[16,13] >= 0
 ldMC3: - ldx[17,1] - ldx[17,13] + ldw[16,13] <= 0
 ldMC3: lDX[17,1] + lDX[17,13] - lDW[16,13] >= 0
 ldMC3: - lDX[17,1] + lDX[17,13] + lDW[16,13] >= 0
 ldMC3: lDX[17,1] - lDX[17,13] + lDW[16,13] >= 0
 ldMC3: - lDX[17,1] - lDX[17,13] + lDW[16,13] <= 0
 ldMC0: - ldx[17,6] + ldw[16,2] = 0
 ldMC0: - lDX[17,6] + lDW[16,2] = 0
 ldMC1: ldw[16,6] - ldw[16,10] >= 0
 ldMC1: - ldx[17,10] + ldw[16,6] >= 0
 ldMC1: - ldx[17,10] + ldw[16,6] - ldw[16,10] <= 0
 ldMC1: lDX[17,10] - lDW[16,6] + lDW[16,10] >= 0
 ldMC1: lDX[17,10] + lDW[16,6] - lDW[16,10] >= 0
 ldMC1: - lDX[17,10] + lDW[16,6] + lDW[16,10] >= 0
 ldMC1: - lDX[17,10] + lDW[16,6] - lDW[16,10] <= 0
 ldMC2: - ldx[17,6] + ldw[16,10] >= 0
 ldMC2: - ldx[17,14] + ldw[16,10] >= 0
 ldMC2: - ldx[17,6] - ldx[17,14] + ldw[16,10] <= 0
 ldMC2: lDX[17,6] + lDX[17,14] - lDW[16,10] >= 0
 ldMC2: - lDX[17,6] + lDX[17,14] + lDW[16,10] >= 0
 ldMC2: lDX[17,6] - lDX[17,14] + lDW[16,10] >= 0
 ldMC2: - lDX[17,6] - lDX[17,14] + lDW[16,10] <= 0
 ldMC3: - ldx[17,2] + ldw[16,14] >= 0
 ldMC3: - ldx[17,14] + ldw[16,14] >= 0
 ldMC3: - ldx[17,2] - ldx[17,14] + ldw[16,14] <= 0
 ldMC3: lDX[17,2] + lDX[17,14] - lDW[16,14] >= 0
 ldMC3: - lDX[17,2] + lDX[17,14] + lDW[16,14] >= 0
 ldMC3: lDX[17,2] - lDX[17,14] + lDW[16,14] >= 0
 ldMC3: - lDX[17,2] - lDX[17,14] + lDW[16,14] <= 0
 ldMC0: - ldx[17,7] + ldw[16,3] = 0
 ldMC0: - lDX[17,7] + lDW[16,3] = 0
 ldMC1: ldw[16,7] - ldw[16,11] >= 0
 ldMC1: - ldx[17,11] + ldw[16,7] >= 0
 ldMC1: - ldx[17,11] + ldw[16,7] - ldw[16,11] <= 0
 ldMC1: lDX[17,11] - lDW[16,7] + lDW[16,11] >= 0
 ldMC1: lDX[17,11] + lDW[16,7] - lDW[16,11] >= 0
 ldMC1: - lDX[17,11] + lDW[16,7] + lDW[16,11] >= 0
 ldMC1: - lDX[17,11] + lDW[16,7] - lDW[16,11] <= 0
 ldMC2: - ldx[17,7] + ldw[16,11] >= 0
 ldMC2: - ldx[17,15] + ldw[16,11] >= 0
 ldMC2: - ldx[17,7] - ldx[17,15] + ldw[16,11] <= 0
 ldMC2: lDX[17,7] + lDX[17,15] - lDW[16,11] >= 0
 ldMC2: - lDX[17,7] + lDX[17,15] + lDW[16,11] >= 0
 ldMC2: lDX[17,7] - lDX[17,15] + lDW[16,11] >= 0
 ldMC2: - lDX[17,7] - lDX[17,15] + lDW[16,11] <= 0
 ldMC3: - ldx[17,3] + ldw[16,15] >= 0
 ldMC3: - ldx[17,15] + ldw[16,15] >= 0
 ldMC3: - ldx[17,3] - ldx[17,15] + ldw[16,15] <= 0
 ldMC3: lDX[17,3] + lDX[17,15] - lDW[16,15] >= 0
 ldMC3: - lDX[17,3] + lDX[17,15] + lDW[16,15] >= 0
 ldMC3: lDX[17,3] - lDX[17,15] + lDW[16,15] >= 0
 ldMC3: - lDX[17,3] - lDX[17,15] + lDW[16,15] <= 0
 ldMC0: - ldx[18,4] + ldw[17,0] = 0
 ldMC0: - lDX[18,4] + lDW[17,0] = 0
 ldMC1: ldw[17,4] - ldw[17,8] >= 0
 ldMC1: - ldx[18,8] + ldw[17,4] >= 0
 ldMC1: - ldx[18,8] + ldw[17,4] - ldw[17,8] <= 0
 ldMC1: lDX[18,8] - lDW[17,4] + lDW[17,8] >= 0
 ldMC1: lDX[18,8] + lDW[17,4] - lDW[17,8] >= 0
 ldMC1: - lDX[18,8] + lDW[17,4] + lDW[17,8] >= 0
 ldMC1: - lDX[18,8] + lDW[17,4] - lDW[17,8] <= 0
 ldMC2: - ldx[18,4] + ldw[17,8] >= 0
 ldMC2: - ldx[18,12] + ldw[17,8] >= 0
 ldMC2: - ldx[18,4] - ldx[18,12] + ldw[17,8] <= 0
 ldMC2: lDX[18,4] + lDX[18,12] - lDW[17,8] >= 0
 ldMC2: - lDX[18,4] + lDX[18,12] + lDW[17,8] >= 0
 ldMC2: lDX[18,4] - lDX[18,12] + lDW[17,8] >= 0
 ldMC2: - lDX[18,4] - lDX[18,12] + lDW[17,8] <= 0
 ldMC3: - ldx[18,0] + ldw[17,12] >= 0
 ldMC3: - ldx[18,12] + ldw[17,12] >= 0
 ldMC3: - ldx[18,0] - ldx[18,12] + ldw[17,12] <= 0
 ldMC3: lDX[18,0] + lDX[18,12] - lDW[17,12] >= 0
 ldMC3: - lDX[18,0] + lDX[18,12] + lDW[17,12] >= 0
 ldMC3: lDX[18,0] - lDX[18,12] + lDW[17,12] >= 0
 ldMC3: - lDX[18,0] - lDX[18,12] + lDW[17,12] <= 0
 ldMC0: - ldx[18,5] + ldw[17,1] = 0
 ldMC0: - lDX[18,5] + lDW[17,1] = 0
 ldMC1: ldw[17,5] - ldw[17,9] >= 0
 ldMC1: - ldx[18,9] + ldw[17,5] >= 0
 ldMC1: - ldx[18,9] + ldw[17,5] - ldw[17,9] <= 0
 ldMC1: lDX[18,9] - lDW[17,5] + lDW[17,9] >= 0
 ldMC1: lDX[18,9] + lDW[17,5] - lDW[17,9] >= 0
 ldMC1: - lDX[18,9] + lDW[17,5] + lDW[17,9] >= 0
 ldMC1: - lDX[18,9] + lDW[17,5] - lDW[17,9] <= 0
 ldMC2: - ldx[18,5] + ldw[17,9] >= 0
 ldMC2: - ldx[18,13] + ldw[17,9] >= 0
 ldMC2: - ldx[18,5] - ldx[18,13] + ldw[17,9] <= 0
 ldMC2: lDX[18,5] + lDX[18,13] - lDW[17,9] >= 0
 ldMC2: - lDX[18,5] + lDX[18,13] + lDW[17,9] >= 0
 ldMC2: lDX[18,5] - lDX[18,13] + lDW[17,9] >= 0
 ldMC2: - lDX[18,5] - lDX[18,13] + lDW[17,9] <= 0
 ldMC3: - ldx[18,1] + ldw[17,13] >= 0
 ldMC3: - ldx[18,13] + ldw[17,13] >= 0
 ldMC3: - ldx[18,1] - ldx[18,13] + ldw[17,13] <= 0
 ldMC3: lDX[18,1] + lDX[18,13] - lDW[17,13] >= 0
 ldMC3: - lDX[18,1] + lDX[18,13] + lDW[17,13] >= 0
 ldMC3: lDX[18,1] - lDX[18,13] + lDW[17,13] >= 0
 ldMC3: - lDX[18,1] - lDX[18,13] + lDW[17,13] <= 0
 ldMC0: - ldx[18,6] + ldw[17,2] = 0
 ldMC0: - lDX[18,6] + lDW[17,2] = 0
 ldMC1: ldw[17,6] - ldw[17,10] >= 0
 ldMC1: - ldx[18,10] + ldw[17,6] >= 0
 ldMC1: - ldx[18,10] + ldw[17,6] - ldw[17,10] <= 0
 ldMC1: lDX[18,10] - lDW[17,6] + lDW[17,10] >= 0
 ldMC1: lDX[18,10] + lDW[17,6] - lDW[17,10] >= 0
 ldMC1: - lDX[18,10] + lDW[17,6] + lDW[17,10] >= 0
 ldMC1: - lDX[18,10] + lDW[17,6] - lDW[17,10] <= 0
 ldMC2: - ldx[18,6] + ldw[17,10] >= 0
 ldMC2: - ldx[18,14] + ldw[17,10] >= 0
 ldMC2: - ldx[18,6] - ldx[18,14] + ldw[17,10] <= 0
 ldMC2: lDX[18,6] + lDX[18,14] - lDW[17,10] >= 0
 ldMC2: - lDX[18,6] + lDX[18,14] + lDW[17,10] >= 0
 ldMC2: lDX[18,6] - lDX[18,14] + lDW[17,10] >= 0
 ldMC2: - lDX[18,6] - lDX[18,14] + lDW[17,10] <= 0
 ldMC3: - ldx[18,2] + ldw[17,14] >= 0
 ldMC3: - ldx[18,14] + ldw[17,14] >= 0
 ldMC3: - ldx[18,2] - ldx[18,14] + ldw[17,14] <= 0
 ldMC3: lDX[18,2] + lDX[18,14] - lDW[17,14] >= 0
 ldMC3: - lDX[18,2] + lDX[18,14] + lDW[17,14] >= 0
 ldMC3: lDX[18,2] - lDX[18,14] + lDW[17,14] >= 0
 ldMC3: - lDX[18,2] - lDX[18,14] + lDW[17,14] <= 0
 ldMC0: - ldx[18,7] + ldw[17,3] = 0
 ldMC0: - lDX[18,7] + lDW[17,3] = 0
 ldMC1: ldw[17,7] - ldw[17,11] >= 0
 ldMC1: - ldx[18,11] + ldw[17,7] >= 0
 ldMC1: - ldx[18,11] + ldw[17,7] - ldw[17,11] <= 0
 ldMC1: lDX[18,11] - lDW[17,7] + lDW[17,11] >= 0
 ldMC1: lDX[18,11] + lDW[17,7] - lDW[17,11] >= 0
 ldMC1: - lDX[18,11] + lDW[17,7] + lDW[17,11] >= 0
 ldMC1: - lDX[18,11] + lDW[17,7] - lDW[17,11] <= 0
 ldMC2: - ldx[18,7] + ldw[17,11] >= 0
 ldMC2: - ldx[18,15] + ldw[17,11] >= 0
 ldMC2: - ldx[18,7] - ldx[18,15] + ldw[17,11] <= 0
 ldMC2: lDX[18,7] + lDX[18,15] - lDW[17,11] >= 0
 ldMC2: - lDX[18,7] + lDX[18,15] + lDW[17,11] >= 0
 ldMC2: lDX[18,7] - lDX[18,15] + lDW[17,11] >= 0
 ldMC2: - lDX[18,7] - lDX[18,15] + lDW[17,11] <= 0
 ldMC3: - ldx[18,3] + ldw[17,15] >= 0
 ldMC3: - ldx[18,15] + ldw[17,15] >= 0
 ldMC3: - ldx[18,3] - ldx[18,15] + ldw[17,15] <= 0
 ldMC3: lDX[18,3] + lDX[18,15] - lDW[17,15] >= 0
 ldMC3: - lDX[18,3] + lDX[18,15] + lDW[17,15] >= 0
 ldMC3: lDX[18,3] - lDX[18,15] + lDW[17,15] >= 0
 ldMC3: - lDX[18,3] - lDX[18,15] + lDW[17,15] <= 0
 ldMC0: - ldx[19,4] + ldw[18,0] = 0
 ldMC0: - lDX[19,4] + lDW[18,0] = 0
 ldMC1: ldw[18,4] - ldw[18,8] >= 0
 ldMC1: - ldx[19,8] + ldw[18,4] >= 0
 ldMC1: - ldx[19,8] + ldw[18,4] - ldw[18,8] <= 0
 ldMC1: lDX[19,8] - lDW[18,4] + lDW[18,8] >= 0
 ldMC1: lDX[19,8] + lDW[18,4] - lDW[18,8] >= 0
 ldMC1: - lDX[19,8] + lDW[18,4] + lDW[18,8] >= 0
 ldMC1: - lDX[19,8] + lDW[18,4] - lDW[18,8] <= 0
 ldMC2: - ldx[19,4] + ldw[18,8] >= 0
 ldMC2: - ldx[19,12] + ldw[18,8] >= 0
 ldMC2: - ldx[19,4] - ldx[19,12] + ldw[18,8] <= 0
 ldMC2: lDX[19,4] + lDX[19,12] - lDW[18,8] >= 0
 ldMC2: - lDX[19,4] + lDX[19,12] + lDW[18,8] >= 0
 ldMC2: lDX[19,4] - lDX[19,12] + lDW[18,8] >= 0
 ldMC2: - lDX[19,4] - lDX[19,12] + lDW[18,8] <= 0
 ldMC3: - ldx[19,0] + ldw[18,12] >= 0
 ldMC3: - ldx[19,12] + ldw[18,12] >= 0
 ldMC3: - ldx[19,0] - ldx[19,12] + ldw[18,12] <= 0
 ldMC3: lDX[19,0] + lDX[19,12] - lDW[18,12] >= 0
 ldMC3: - lDX[19,0] + lDX[19,12] + lDW[18,12] >= 0
 ldMC3: lDX[19,0] - lDX[19,12] + lDW[18,12] >= 0
 ldMC3: - lDX[19,0] - lDX[19,12] + lDW[18,12] <= 0
 ldMC0: - ldx[19,5] + ldw[18,1] = 0
 ldMC0: - lDX[19,5] + lDW[18,1] = 0
 ldMC1: ldw[18,5] - ldw[18,9] >= 0
 ldMC1: - ldx[19,9] + ldw[18,5] >= 0
 ldMC1: - ldx[19,9] + ldw[18,5] - ldw[18,9] <= 0
 ldMC1: lDX[19,9] - lDW[18,5] + lDW[18,9] >= 0
 ldMC1: lDX[19,9] + lDW[18,5] - lDW[18,9] >= 0
 ldMC1: - lDX[19,9] + lDW[18,5] + lDW[18,9] >= 0
 ldMC1: - lDX[19,9] + lDW[18,5] - lDW[18,9] <= 0
 ldMC2: - ldx[19,5] + ldw[18,9] >= 0
 ldMC2: - ldx[19,13] + ldw[18,9] >= 0
 ldMC2: - ldx[19,5] - ldx[19,13] + ldw[18,9] <= 0
 ldMC2: lDX[19,5] + lDX[19,13] - lDW[18,9] >= 0
 ldMC2: - lDX[19,5] + lDX[19,13] + lDW[18,9] >= 0
 ldMC2: lDX[19,5] - lDX[19,13] + lDW[18,9] >= 0
 ldMC2: - lDX[19,5] - lDX[19,13] + lDW[18,9] <= 0
 ldMC3: - ldx[19,1] + ldw[18,13] >= 0
 ldMC3: - ldx[19,13] + ldw[18,13] >= 0
 ldMC3: - ldx[19,1] - ldx[19,13] + ldw[18,13] <= 0
 ldMC3: lDX[19,1] + lDX[19,13] - lDW[18,13] >= 0
 ldMC3: - lDX[19,1] + lDX[19,13] + lDW[18,13] >= 0
 ldMC3: lDX[19,1] - lDX[19,13] + lDW[18,13] >= 0
 ldMC3: - lDX[19,1] - lDX[19,13] + lDW[18,13] <= 0
 ldMC0: - ldx[19,6] + ldw[18,2] = 0
 ldMC0: - lDX[19,6] + lDW[18,2] = 0
 ldMC1: ldw[18,6] - ldw[18,10] >= 0
 ldMC1: - ldx[19,10] + ldw[18,6] >= 0
 ldMC1: - ldx[19,10] + ldw[18,6] - ldw[18,10] <= 0
 ldMC1: lDX[19,10] - lDW[18,6] + lDW[18,10] >= 0
 ldMC1: lDX[19,10] + lDW[18,6] - lDW[18,10] >= 0
 ldMC1: - lDX[19,10] + lDW[18,6] + lDW[18,10] >= 0
 ldMC1: - lDX[19,10] + lDW[18,6] - lDW[18,10] <= 0
 ldMC2: - ldx[19,6] + ldw[18,10] >= 0
 ldMC2: - ldx[19,14] + ldw[18,10] >= 0
 ldMC2: - ldx[19,6] - ldx[19,14] + ldw[18,10] <= 0
 ldMC2: lDX[19,6] + lDX[19,14] - lDW[18,10] >= 0
 ldMC2: - lDX[19,6] + lDX[19,14] + lDW[18,10] >= 0
 ldMC2: lDX[19,6] - lDX[19,14] + lDW[18,10] >= 0
 ldMC2: - lDX[19,6] - lDX[19,14] + lDW[18,10] <= 0
 ldMC3: - ldx[19,2] + ldw[18,14] >= 0
 ldMC3: - ldx[19,14] + ldw[18,14] >= 0
 ldMC3: - ldx[19,2] - ldx[19,14] + ldw[18,14] <= 0
 ldMC3: lDX[19,2] + lDX[19,14] - lDW[18,14] >= 0
 ldMC3: - lDX[19,2] + lDX[19,14] + lDW[18,14] >= 0
 ldMC3: lDX[19,2] - lDX[19,14] + lDW[18,14] >= 0
 ldMC3: - lDX[19,2] - lDX[19,14] + lDW[18,14] <= 0
 ldMC0: - ldx[19,7] + ldw[18,3] = 0
 ldMC0: - lDX[19,7] + lDW[18,3] = 0
 ldMC1: ldw[18,7] - ldw[18,11] >= 0
 ldMC1: - ldx[19,11] + ldw[18,7] >= 0
 ldMC1: - ldx[19,11] + ldw[18,7] - ldw[18,11] <= 0
 ldMC1: lDX[19,11] - lDW[18,7] + lDW[18,11] >= 0
 ldMC1: lDX[19,11] + lDW[18,7] - lDW[18,11] >= 0
 ldMC1: - lDX[19,11] + lDW[18,7] + lDW[18,11] >= 0
 ldMC1: - lDX[19,11] + lDW[18,7] - lDW[18,11] <= 0
 ldMC2: - ldx[19,7] + ldw[18,11] >= 0
 ldMC2: - ldx[19,15] + ldw[18,11] >= 0
 ldMC2: - ldx[19,7] - ldx[19,15] + ldw[18,11] <= 0
 ldMC2: lDX[19,7] + lDX[19,15] - lDW[18,11] >= 0
 ldMC2: - lDX[19,7] + lDX[19,15] + lDW[18,11] >= 0
 ldMC2: lDX[19,7] - lDX[19,15] + lDW[18,11] >= 0
 ldMC2: - lDX[19,7] - lDX[19,15] + lDW[18,11] <= 0
 ldMC3: - ldx[19,3] + ldw[18,15] >= 0
 ldMC3: - ldx[19,15] + ldw[18,15] >= 0
 ldMC3: - ldx[19,3] - ldx[19,15] + ldw[18,15] <= 0
 ldMC3: lDX[19,3] + lDX[19,15] - lDW[18,15] >= 0
 ldMC3: - lDX[19,3] + lDX[19,15] + lDW[18,15] >= 0
 ldMC3: lDX[19,3] - lDX[19,15] + lDW[18,15] >= 0
 ldMC3: - lDX[19,3] - lDX[19,15] + lDW[18,15] <= 0
 leMC0: ldx[20,0] - ldw[19,12] >= 0
 leMC0: ldx[20,0] - ldx[20,12] >= 0
 leMC0: ldx[20,0] - ldx[20,12] - ldw[19,12] <= 0
 leMC0: - lDX[20,0] + lDX[20,12] + lDW[19,12] >= 0
 leMC0: lDX[20,0] + lDX[20,12] - lDW[19,12] >= 0
 leMC0: lDX[20,0] - lDX[20,12] + lDW[19,12] >= 0
 leMC0: lDX[20,0] - lDX[20,12] - lDW[19,12] <= 0
 leMC1: ldx[20,4] - ldw[19,0] = 0
 leMC1: lDX[20,4] - lDW[19,0] = 0
 leMC2: ldx[20,8] - ldw[19,4] >= 0
 leMC2: ldx[20,8] - ldw[19,8] >= 0
 leMC2: ldx[20,8] - ldw[19,4] - ldw[19,8] <= 0
 leMC2: - lDX[20,8] + lDW[19,4] + lDW[19,8] >= 0
 leMC2: lDX[20,8] - lDW[19,4] + lDW[19,8] >= 0
 leMC2: lDX[20,8] + lDW[19,4] - lDW[19,8] >= 0
 leMC2: lDX[20,8] - lDW[19,4] - lDW[19,8] <= 0
 leMC3: ldx[20,12] - ldw[19,0] >= 0
 leMC3: ldx[20,12] - ldw[19,8] >= 0
 leMC3: ldx[20,12] - ldw[19,0] - ldw[19,8] <= 0
 leMC3: - lDX[20,12] + lDW[19,0] + lDW[19,8] >= 0
 leMC3: lDX[20,12] - lDW[19,0] + lDW[19,8] >= 0
 leMC3: lDX[20,12] + lDW[19,0] - lDW[19,8] >= 0
 leMC3: lDX[20,12] - lDW[19,0] - lDW[19,8] <= 0
 leMC0: ldx[20,1] - ldw[19,13] >= 0
 leMC0: ldx[20,1] - ldx[20,13] >= 0
 leMC0: ldx[20,1] - ldx[20,13] - ldw[19,13] <= 0
 leMC0: - lDX[20,1] + lDX[20,13] + lDW[19,13] >= 0
 leMC0: lDX[20,1] + lDX[20,13] - lDW[19,13] >= 0
 leMC0: lDX[20,1] - lDX[20,13] + lDW[19,13] >= 0
 leMC0: lDX[20,1] - lDX[20,13] - lDW[19,13] <= 0
 leMC1: ldx[20,5] - ldw[19,1] = 0
 leMC1: lDX[20,5] - lDW[19,1] = 0
 leMC2: ldx[20,9] - ldw[19,5] >= 0
 leMC2: ldx[20,9] - ldw[19,9] >= 0
 leMC2: ldx[20,9] - ldw[19,5] - ldw[19,9] <= 0
 leMC2: - lDX[20,9] + lDW[19,5] + lDW[19,9] >= 0
 leMC2: lDX[20,9] - lDW[19,5] + lDW[19,9] >= 0
 leMC2: lDX[20,9] + lDW[19,5] - lDW[19,9] >= 0
 leMC2: lDX[20,9] - lDW[19,5] - lDW[19,9] <= 0
 leMC3: ldx[20,13] - ldw[19,1] >= 0
 leMC3: ldx[20,13] - ldw[19,9] >= 0
 leMC3: ldx[20,13] - ldw[19,1] - ldw[19,9] <= 0
 leMC3: - lDX[20,13] + lDW[19,1] + lDW[19,9] >= 0
 leMC3: lDX[20,13] - lDW[19,1] + lDW[19,9] >= 0
 leMC3: lDX[20,13] + lDW[19,1] - lDW[19,9] >= 0
 leMC3: lDX[20,13] - lDW[19,1] - lDW[19,9] <= 0
 leMC0: ldx[20,2] - ldw[19,14] >= 0
 leMC0: ldx[20,2] - ldx[20,14] >= 0
 leMC0: ldx[20,2] - ldx[20,14] - ldw[19,14] <= 0
 leMC0: - lDX[20,2] + lDX[20,14] + lDW[19,14] >= 0
 leMC0: lDX[20,2] + lDX[20,14] - lDW[19,14] >= 0
 leMC0: lDX[20,2] - lDX[20,14] + lDW[19,14] >= 0
 leMC0: lDX[20,2] - lDX[20,14] - lDW[19,14] <= 0
 leMC1: ldx[20,6] - ldw[19,2] = 0
 leMC1: lDX[20,6] - lDW[19,2] = 0
 leMC2: ldx[20,10] - ldw[19,6] >= 0
 leMC2: ldx[20,10] - ldw[19,10] >= 0
 leMC2: ldx[20,10] - ldw[19,6] - ldw[19,10] <= 0
 leMC2: - lDX[20,10] + lDW[19,6] + lDW[19,10] >= 0
 leMC2: lDX[20,10] - lDW[19,6] + lDW[19,10] >= 0
 leMC2: lDX[20,10] + lDW[19,6] - lDW[19,10] >= 0
 leMC2: lDX[20,10] - lDW[19,6] - lDW[19,10] <= 0
 leMC3: ldx[20,14] - ldw[19,2] >= 0
 leMC3: ldx[20,14] - ldw[19,10] >= 0
 leMC3: ldx[20,14] - ldw[19,2] - ldw[19,10] <= 0
 leMC3: - lDX[20,14] + lDW[19,2] + lDW[19,10] >= 0
 leMC3: lDX[20,14] - lDW[19,2] + lDW[19,10] >= 0
 leMC3: lDX[20,14] + lDW[19,2] - lDW[19,10] >= 0
 leMC3: lDX[20,14] - lDW[19,2] - lDW[19,10] <= 0
 leMC0: ldx[20,3] - ldw[19,15] >= 0
 leMC0: ldx[20,3] - ldx[20,15] >= 0
 leMC0: ldx[20,3] - ldx[20,15] - ldw[19,15] <= 0
 leMC0: - lDX[20,3] + lDX[20,15] + lDW[19,15] >= 0
 leMC0: lDX[20,3] + lDX[20,15] - lDW[19,15] >= 0
 leMC0: lDX[20,3] - lDX[20,15] + lDW[19,15] >= 0
 leMC0: lDX[20,3] - lDX[20,15] - lDW[19,15] <= 0
 leMC1: ldx[20,7] - ldw[19,3] = 0
 leMC1: lDX[20,7] - lDW[19,3] = 0
 leMC2: ldx[20,11] - ldw[19,7] >= 0
 leMC2: ldx[20,11] - ldw[19,11] >= 0
 leMC2: ldx[20,11] - ldw[19,7] - ldw[19,11] <= 0
 leMC2: - lDX[20,11] + lDW[19,7] + lDW[19,11] >= 0
 leMC2: lDX[20,11] - lDW[19,7] + lDW[19,11] >= 0
 leMC2: lDX[20,11] + lDW[19,7] - lDW[19,11] >= 0
 leMC2: lDX[20,11] - lDW[19,7] - lDW[19,11] <= 0
 leMC3: ldx[20,15] - ldw[19,3] >= 0
 leMC3: ldx[20,15] - ldw[19,11] >= 0
 leMC3: ldx[20,15] - ldw[19,3] - ldw[19,11] <= 0
 leMC3: - lDX[20,15] + lDW[19,3] + lDW[19,11] >= 0
 leMC3: lDX[20,15] - lDW[19,3] + lDW[19,11] >= 0
 leMC3: lDX[20,15] + lDW[19,3] - lDW[19,11] >= 0
 leMC3: lDX[20,15] - lDW[19,3] - lDW[19,11] <= 0
 leMC0: ldx[21,0] - ldw[20,12] >= 0
 leMC0: ldx[21,0] - ldx[21,12] >= 0
 leMC0: ldx[21,0] - ldx[21,12] - ldw[20,12] <= 0
 leMC0: - lDX[21,0] + lDX[21,12] + lDW[20,12] >= 0
 leMC0: lDX[21,0] + lDX[21,12] - lDW[20,12] >= 0
 leMC0: lDX[21,0] - lDX[21,12] + lDW[20,12] >= 0
 leMC0: lDX[21,0] - lDX[21,12] - lDW[20,12] <= 0
 leMC1: ldx[21,4] - ldw[20,0] = 0
 leMC1: lDX[21,4] - lDW[20,0] = 0
 leMC2: ldx[21,8] - ldw[20,4] >= 0
 leMC2: ldx[21,8] - ldw[20,8] >= 0
 leMC2: ldx[21,8] - ldw[20,4] - ldw[20,8] <= 0
 leMC2: - lDX[21,8] + lDW[20,4] + lDW[20,8] >= 0
 leMC2: lDX[21,8] - lDW[20,4] + lDW[20,8] >= 0
 leMC2: lDX[21,8] + lDW[20,4] - lDW[20,8] >= 0
 leMC2: lDX[21,8] - lDW[20,4] - lDW[20,8] <= 0
 leMC3: ldx[21,12] - ldw[20,0] >= 0
 leMC3: ldx[21,12] - ldw[20,8] >= 0
 leMC3: ldx[21,12] - ldw[20,0] - ldw[20,8] <= 0
 leMC3: - lDX[21,12] + lDW[20,0] + lDW[20,8] >= 0
 leMC3: lDX[21,12] - lDW[20,0] + lDW[20,8] >= 0
 leMC3: lDX[21,12] + lDW[20,0] - lDW[20,8] >= 0
 leMC3: lDX[21,12] - lDW[20,0] - lDW[20,8] <= 0
 leMC0: ldx[21,1] - ldw[20,13] >= 0
 leMC0: ldx[21,1] - ldx[21,13] >= 0
 leMC0: ldx[21,1] - ldx[21,13] - ldw[20,13] <= 0
 leMC0: - lDX[21,1] + lDX[21,13] + lDW[20,13] >= 0
 leMC0: lDX[21,1] + lDX[21,13] - lDW[20,13] >= 0
 leMC0: lDX[21,1] - lDX[21,13] + lDW[20,13] >= 0
 leMC0: lDX[21,1] - lDX[21,13] - lDW[20,13] <= 0
 leMC1: ldx[21,5] - ldw[20,1] = 0
 leMC1: lDX[21,5] - lDW[20,1] = 0
 leMC2: ldx[21,9] - ldw[20,5] >= 0
 leMC2: ldx[21,9] - ldw[20,9] >= 0
 leMC2: ldx[21,9] - ldw[20,5] - ldw[20,9] <= 0
 leMC2: - lDX[21,9] + lDW[20,5] + lDW[20,9] >= 0
 leMC2: lDX[21,9] - lDW[20,5] + lDW[20,9] >= 0
 leMC2: lDX[21,9] + lDW[20,5] - lDW[20,9] >= 0
 leMC2: lDX[21,9] - lDW[20,5] - lDW[20,9] <= 0
 leMC3: ldx[21,13] - ldw[20,1] >= 0
 leMC3: ldx[21,13] - ldw[20,9] >= 0
 leMC3: ldx[21,13] - ldw[20,1] - ldw[20,9] <= 0
 leMC3: - lDX[21,13] + lDW[20,1] + lDW[20,9] >= 0
 leMC3: lDX[21,13] - lDW[20,1] + lDW[20,9] >= 0
 leMC3: lDX[21,13] + lDW[20,1] - lDW[20,9] >= 0
 leMC3: lDX[21,13] - lDW[20,1] - lDW[20,9] <= 0
 leMC0: ldx[21,2] - ldw[20,14] >= 0
 leMC0: ldx[21,2] - ldx[21,14] >= 0
 leMC0: ldx[21,2] - ldx[21,14] - ldw[20,14] <= 0
 leMC0: - lDX[21,2] + lDX[21,14] + lDW[20,14] >= 0
 leMC0: lDX[21,2] + lDX[21,14] - lDW[20,14] >= 0
 leMC0: lDX[21,2] - lDX[21,14] + lDW[20,14] >= 0
 leMC0: lDX[21,2] - lDX[21,14] - lDW[20,14] <= 0
 leMC1: ldx[21,6] - ldw[20,2] = 0
 leMC1: lDX[21,6] - lDW[20,2] = 0
 leMC2: ldx[21,10] - ldw[20,6] >= 0
 leMC2: ldx[21,10] - ldw[20,10] >= 0
 leMC2: ldx[21,10] - ldw[20,6] - ldw[20,10] <= 0
 leMC2: - lDX[21,10] + lDW[20,6] + lDW[20,10] >= 0
 leMC2: lDX[21,10] - lDW[20,6] + lDW[20,10] >= 0
 leMC2: lDX[21,10] + lDW[20,6] - lDW[20,10] >= 0
 leMC2: lDX[21,10] - lDW[20,6] - lDW[20,10] <= 0
 leMC3: ldx[21,14] - ldw[20,2] >= 0
 leMC3: ldx[21,14] - ldw[20,10] >= 0
 leMC3: ldx[21,14] - ldw[20,2] - ldw[20,10] <= 0
 leMC3: - lDX[21,14] + lDW[20,2] + lDW[20,10] >= 0
 leMC3: lDX[21,14] - lDW[20,2] + lDW[20,10] >= 0
 leMC3: lDX[21,14] + lDW[20,2] - lDW[20,10] >= 0
 leMC3: lDX[21,14] - lDW[20,2] - lDW[20,10] <= 0
 leMC0: ldx[21,3] - ldw[20,15] >= 0
 leMC0: ldx[21,3] - ldx[21,15] >= 0
 leMC0: ldx[21,3] - ldx[21,15] - ldw[20,15] <= 0
 leMC0: - lDX[21,3] + lDX[21,15] + lDW[20,15] >= 0
 leMC0: lDX[21,3] + lDX[21,15] - lDW[20,15] >= 0
 leMC0: lDX[21,3] - lDX[21,15] + lDW[20,15] >= 0
 leMC0: lDX[21,3] - lDX[21,15] - lDW[20,15] <= 0
 leMC1: ldx[21,7] - ldw[20,3] = 0
 leMC1: lDX[21,7] - lDW[20,3] = 0
 leMC2: ldx[21,11] - ldw[20,7] >= 0
 leMC2: ldx[21,11] - ldw[20,11] >= 0
 leMC2: ldx[21,11] - ldw[20,7] - ldw[20,11] <= 0
 leMC2: - lDX[21,11] + lDW[20,7] + lDW[20,11] >= 0
 leMC2: lDX[21,11] - lDW[20,7] + lDW[20,11] >= 0
 leMC2: lDX[21,11] + lDW[20,7] - lDW[20,11] >= 0
 leMC2: lDX[21,11] - lDW[20,7] - lDW[20,11] <= 0
 leMC3: ldx[21,15] - ldw[20,3] >= 0
 leMC3: ldx[21,15] - ldw[20,11] >= 0
 leMC3: ldx[21,15] - ldw[20,3] - ldw[20,11] <= 0
 leMC3: - lDX[21,15] + lDW[20,3] + lDW[20,11] >= 0
 leMC3: lDX[21,15] - lDW[20,3] + lDW[20,11] >= 0
 leMC3: lDX[21,15] + lDW[20,3] - lDW[20,11] >= 0
 leMC3: lDX[21,15] - lDW[20,3] - lDW[20,11] <= 0
 leMC0: ldx[22,0] - ldw[21,12] >= 0
 leMC0: ldx[22,0] - ldx[22,12] >= 0
 leMC0: ldx[22,0] - ldx[22,12] - ldw[21,12] <= 0
 leMC0: - lDX[22,0] + lDX[22,12] + lDW[21,12] >= 0
 leMC0: lDX[22,0] + lDX[22,12] - lDW[21,12] >= 0
 leMC0: lDX[22,0] - lDX[22,12] + lDW[21,12] >= 0
 leMC0: lDX[22,0] - lDX[22,12] - lDW[21,12] <= 0
 leMC1: ldx[22,4] - ldw[21,0] = 0
 leMC1: lDX[22,4] - lDW[21,0] = 0
 leMC2: ldx[22,8] - ldw[21,4] >= 0
 leMC2: ldx[22,8] - ldw[21,8] >= 0
 leMC2: ldx[22,8] - ldw[21,4] - ldw[21,8] <= 0
 leMC2: - lDX[22,8] + lDW[21,4] + lDW[21,8] >= 0
 leMC2: lDX[22,8] - lDW[21,4] + lDW[21,8] >= 0
 leMC2: lDX[22,8] + lDW[21,4] - lDW[21,8] >= 0
 leMC2: lDX[22,8] - lDW[21,4] - lDW[21,8] <= 0
 leMC3: ldx[22,12] - ldw[21,0] >= 0
 leMC3: ldx[22,12] - ldw[21,8] >= 0
 leMC3: ldx[22,12] - ldw[21,0] - ldw[21,8] <= 0
 leMC3: - lDX[22,12] + lDW[21,0] + lDW[21,8] >= 0
 leMC3: lDX[22,12] - lDW[21,0] + lDW[21,8] >= 0
 leMC3: lDX[22,12] + lDW[21,0] - lDW[21,8] >= 0
 leMC3: lDX[22,12] - lDW[21,0] - lDW[21,8] <= 0
 leMC0: ldx[22,1] - ldw[21,13] >= 0
 leMC0: ldx[22,1] - ldx[22,13] >= 0
 leMC0: ldx[22,1] - ldx[22,13] - ldw[21,13] <= 0
 leMC0: - lDX[22,1] + lDX[22,13] + lDW[21,13] >= 0
 leMC0: lDX[22,1] + lDX[22,13] - lDW[21,13] >= 0
 leMC0: lDX[22,1] - lDX[22,13] + lDW[21,13] >= 0
 leMC0: lDX[22,1] - lDX[22,13] - lDW[21,13] <= 0
 leMC1: ldx[22,5] - ldw[21,1] = 0
 leMC1: lDX[22,5] - lDW[21,1] = 0
 leMC2: ldx[22,9] - ldw[21,5] >= 0
 leMC2: ldx[22,9] - ldw[21,9] >= 0
 leMC2: ldx[22,9] - ldw[21,5] - ldw[21,9] <= 0
 leMC2: - lDX[22,9] + lDW[21,5] + lDW[21,9] >= 0
 leMC2: lDX[22,9] - lDW[21,5] + lDW[21,9] >= 0
 leMC2: lDX[22,9] + lDW[21,5] - lDW[21,9] >= 0
 leMC2: lDX[22,9] - lDW[21,5] - lDW[21,9] <= 0
 leMC3: ldx[22,13] - ldw[21,1] >= 0
 leMC3: ldx[22,13] - ldw[21,9] >= 0
 leMC3: ldx[22,13] - ldw[21,1] - ldw[21,9] <= 0
 leMC3: - lDX[22,13] + lDW[21,1] + lDW[21,9] >= 0
 leMC3: lDX[22,13] - lDW[21,1] + lDW[21,9] >= 0
 leMC3: lDX[22,13] + lDW[21,1] - lDW[21,9] >= 0
 leMC3: lDX[22,13] - lDW[21,1] - lDW[21,9] <= 0
 leMC0: ldx[22,2] - ldw[21,14] >= 0
 leMC0: ldx[22,2] - ldx[22,14] >= 0
 leMC0: ldx[22,2] - ldx[22,14] - ldw[21,14] <= 0
 leMC0: - lDX[22,2] + lDX[22,14] + lDW[21,14] >= 0
 leMC0: lDX[22,2] + lDX[22,14] - lDW[21,14] >= 0
 leMC0: lDX[22,2] - lDX[22,14] + lDW[21,14] >= 0
 leMC0: lDX[22,2] - lDX[22,14] - lDW[21,14] <= 0
 leMC1: ldx[22,6] - ldw[21,2] = 0
 leMC1: lDX[22,6] - lDW[21,2] = 0
 leMC2: ldx[22,10] - ldw[21,6] >= 0
 leMC2: ldx[22,10] - ldw[21,10] >= 0
 leMC2: ldx[22,10] - ldw[21,6] - ldw[21,10] <= 0
 leMC2: - lDX[22,10] + lDW[21,6] + lDW[21,10] >= 0
 leMC2: lDX[22,10] - lDW[21,6] + lDW[21,10] >= 0
 leMC2: lDX[22,10] + lDW[21,6] - lDW[21,10] >= 0
 leMC2: lDX[22,10] - lDW[21,6] - lDW[21,10] <= 0
 leMC3: ldx[22,14] - ldw[21,2] >= 0
 leMC3: ldx[22,14] - ldw[21,10] >= 0
 leMC3: ldx[22,14] - ldw[21,2] - ldw[21,10] <= 0
 leMC3: - lDX[22,14] + lDW[21,2] + lDW[21,10] >= 0
 leMC3: lDX[22,14] - lDW[21,2] + lDW[21,10] >= 0
 leMC3: lDX[22,14] + lDW[21,2] - lDW[21,10] >= 0
 leMC3: lDX[22,14] - lDW[21,2] - lDW[21,10] <= 0
 leMC0: ldx[22,3] - ldw[21,15] >= 0
 leMC0: ldx[22,3] - ldx[22,15] >= 0
 leMC0: ldx[22,3] - ldx[22,15] - ldw[21,15] <= 0
 leMC0: - lDX[22,3] + lDX[22,15] + lDW[21,15] >= 0
 leMC0: lDX[22,3] + lDX[22,15] - lDW[21,15] >= 0
 leMC0: lDX[22,3] - lDX[22,15] + lDW[21,15] >= 0
 leMC0: lDX[22,3] - lDX[22,15] - lDW[21,15] <= 0
 leMC1: ldx[22,7] - ldw[21,3] = 0
 leMC1: lDX[22,7] - lDW[21,3] = 0
 leMC2: ldx[22,11] - ldw[21,7] >= 0
 leMC2: ldx[22,11] - ldw[21,11] >= 0
 leMC2: ldx[22,11] - ldw[21,7] - ldw[21,11] <= 0
 leMC2: - lDX[22,11] + lDW[21,7] + lDW[21,11] >= 0
 leMC2: lDX[22,11] - lDW[21,7] + lDW[21,11] >= 0
 leMC2: lDX[22,11] + lDW[21,7] - lDW[21,11] >= 0
 leMC2: lDX[22,11] - lDW[21,7] - lDW[21,11] <= 0
 leMC3: ldx[22,15] - ldw[21,3] >= 0
 leMC3: ldx[22,15] - ldw[21,11] >= 0
 leMC3: ldx[22,15] - ldw[21,3] - ldw[21,11] <= 0
 leMC3: - lDX[22,15] + lDW[21,3] + lDW[21,11] >= 0
 leMC3: lDX[22,15] - lDW[21,3] + lDW[21,11] >= 0
 leMC3: lDX[22,15] + lDW[21,3] - lDW[21,11] >= 0
 leMC3: lDX[22,15] - lDW[21,3] - lDW[21,11] <= 0
 leMC0: ldx[23,0] - ldw[22,12] >= 0
 leMC0: ldx[23,0] - ldx[23,12] >= 0
 leMC0: ldx[23,0] - ldx[23,12] - ldw[22,12] <= 0
 leMC0: - lDX[23,0] + lDX[23,12] + lDW[22,12] >= 0
 leMC0: lDX[23,0] + lDX[23,12] - lDW[22,12] >= 0
 leMC0: lDX[23,0] - lDX[23,12] + lDW[22,12] >= 0
 leMC0: lDX[23,0] - lDX[23,12] - lDW[22,12] <= 0
 leMC1: ldx[23,4] - ldw[22,0] = 0
 leMC1: lDX[23,4] - lDW[22,0] = 0
 leMC2: ldx[23,8] - ldw[22,4] >= 0
 leMC2: ldx[23,8] - ldw[22,8] >= 0
 leMC2: ldx[23,8] - ldw[22,4] - ldw[22,8] <= 0
 leMC2: - lDX[23,8] + lDW[22,4] + lDW[22,8] >= 0
 leMC2: lDX[23,8] - lDW[22,4] + lDW[22,8] >= 0
 leMC2: lDX[23,8] + lDW[22,4] - lDW[22,8] >= 0
 leMC2: lDX[23,8] - lDW[22,4] - lDW[22,8] <= 0
 leMC3: ldx[23,12] - ldw[22,0] >= 0
 leMC3: ldx[23,12] - ldw[22,8] >= 0
 leMC3: ldx[23,12] - ldw[22,0] - ldw[22,8] <= 0
 leMC3: - lDX[23,12] + lDW[22,0] + lDW[22,8] >= 0
 leMC3: lDX[23,12] - lDW[22,0] + lDW[22,8] >= 0
 leMC3: lDX[23,12] + lDW[22,0] - lDW[22,8] >= 0
 leMC3: lDX[23,12] - lDW[22,0] - lDW[22,8] <= 0
 leMC0: ldx[23,1] - ldw[22,13] >= 0
 leMC0: ldx[23,1] - ldx[23,13] >= 0
 leMC0: ldx[23,1] - ldx[23,13] - ldw[22,13] <= 0
 leMC0: - lDX[23,1] + lDX[23,13] + lDW[22,13] >= 0
 leMC0: lDX[23,1] + lDX[23,13] - lDW[22,13] >= 0
 leMC0: lDX[23,1] - lDX[23,13] + lDW[22,13] >= 0
 leMC0: lDX[23,1] - lDX[23,13] - lDW[22,13] <= 0
 leMC1: ldx[23,5] - ldw[22,1] = 0
 leMC1: lDX[23,5] - lDW[22,1] = 0
 leMC2: ldx[23,9] - ldw[22,5] >= 0
 leMC2: ldx[23,9] - ldw[22,9] >= 0
 leMC2: ldx[23,9] - ldw[22,5] - ldw[22,9] <= 0
 leMC2: - lDX[23,9] + lDW[22,5] + lDW[22,9] >= 0
 leMC2: lDX[23,9] - lDW[22,5] + lDW[22,9] >= 0
 leMC2: lDX[23,9] + lDW[22,5] - lDW[22,9] >= 0
 leMC2: lDX[23,9] - lDW[22,5] - lDW[22,9] <= 0
 leMC3: ldx[23,13] - ldw[22,1] >= 0
 leMC3: ldx[23,13] - ldw[22,9] >= 0
 leMC3: ldx[23,13] - ldw[22,1] - ldw[22,9] <= 0
 leMC3: - lDX[23,13] + lDW[22,1] + lDW[22,9] >= 0
 leMC3: lDX[23,13] - lDW[22,1] + lDW[22,9] >= 0
 leMC3: lDX[23,13] + lDW[22,1] - lDW[22,9] >= 0
 leMC3: lDX[23,13] - lDW[22,1] - lDW[22,9] <= 0
 leMC0: ldx[23,2] - ldw[22,14] >= 0
 leMC0: ldx[23,2] - ldx[23,14] >= 0
 leMC0: ldx[23,2] - ldx[23,14] - ldw[22,14] <= 0
 leMC0: - lDX[23,2] + lDX[23,14] + lDW[22,14] >= 0
 leMC0: lDX[23,2] + lDX[23,14] - lDW[22,14] >= 0
 leMC0: lDX[23,2] - lDX[23,14] + lDW[22,14] >= 0
 leMC0: lDX[23,2] - lDX[23,14] - lDW[22,14] <= 0
 leMC1: ldx[23,6] - ldw[22,2] = 0
 leMC1: lDX[23,6] - lDW[22,2] = 0
 leMC2: ldx[23,10] - ldw[22,6] >= 0
 leMC2: ldx[23,10] - ldw[22,10] >= 0
 leMC2: ldx[23,10] - ldw[22,6] - ldw[22,10] <= 0
 leMC2: - lDX[23,10] + lDW[22,6] + lDW[22,10] >= 0
 leMC2: lDX[23,10] - lDW[22,6] + lDW[22,10] >= 0
 leMC2: lDX[23,10] + lDW[22,6] - lDW[22,10] >= 0
 leMC2: lDX[23,10] - lDW[22,6] - lDW[22,10] <= 0
 leMC3: ldx[23,14] - ldw[22,2] >= 0
 leMC3: ldx[23,14] - ldw[22,10] >= 0
 leMC3: ldx[23,14] - ldw[22,2] - ldw[22,10] <= 0
 leMC3: - lDX[23,14] + lDW[22,2] + lDW[22,10] >= 0
 leMC3: lDX[23,14] - lDW[22,2] + lDW[22,10] >= 0
 leMC3: lDX[23,14] + lDW[22,2] - lDW[22,10] >= 0
 leMC3: lDX[23,14] - lDW[22,2] - lDW[22,10] <= 0
 leMC0: ldx[23,3] - ldw[22,15] >= 0
 leMC0: ldx[23,3] - ldx[23,15] >= 0
 leMC0: ldx[23,3] - ldx[23,15] - ldw[22,15] <= 0
 leMC0: - lDX[23,3] + lDX[23,15] + lDW[22,15] >= 0
 leMC0: lDX[23,3] + lDX[23,15] - lDW[22,15] >= 0
 leMC0: lDX[23,3] - lDX[23,15] + lDW[22,15] >= 0
 leMC0: lDX[23,3] - lDX[23,15] - lDW[22,15] <= 0
 leMC1: ldx[23,7] - ldw[22,3] = 0
 leMC1: lDX[23,7] - lDW[22,3] = 0
 leMC2: ldx[23,11] - ldw[22,7] >= 0
 leMC2: ldx[23,11] - ldw[22,11] >= 0
 leMC2: ldx[23,11] - ldw[22,7] - ldw[22,11] <= 0
 leMC2: - lDX[23,11] + lDW[22,7] + lDW[22,11] >= 0
 leMC2: lDX[23,11] - lDW[22,7] + lDW[22,11] >= 0
 leMC2: lDX[23,11] + lDW[22,7] - lDW[22,11] >= 0
 leMC2: lDX[23,11] - lDW[22,7] - lDW[22,11] <= 0
 leMC3: ldx[23,15] - ldw[22,3] >= 0
 leMC3: ldx[23,15] - ldw[22,11] >= 0
 leMC3: ldx[23,15] - ldw[22,3] - ldw[22,11] <= 0
 leMC3: - lDX[23,15] + lDW[22,3] + lDW[22,11] >= 0
 leMC3: lDX[23,15] - lDW[22,3] + lDW[22,11] >= 0
 leMC3: lDX[23,15] + lDW[22,3] - lDW[22,11] >= 0
 leMC3: lDX[23,15] - lDW[22,3] - lDW[22,11] <= 0
 leMC0: ldx[24,0] - ldw[23,12] >= 0
 leMC0: ldx[24,0] - ldx[24,12] >= 0
 leMC0: ldx[24,0] - ldx[24,12] - ldw[23,12] <= 0
 leMC0: - lDX[24,0] + lDX[24,12] + lDW[23,12] >= 0
 leMC0: lDX[24,0] + lDX[24,12] - lDW[23,12] >= 0
 leMC0: lDX[24,0] - lDX[24,12] + lDW[23,12] >= 0
 leMC0: lDX[24,0] - lDX[24,12] - lDW[23,12] <= 0
 leMC1: ldx[24,4] - ldw[23,0] = 0
 leMC1: lDX[24,4] - lDW[23,0] = 0
 leMC2: ldx[24,8] - ldw[23,4] >= 0
 leMC2: ldx[24,8] - ldw[23,8] >= 0
 leMC2: ldx[24,8] - ldw[23,4] - ldw[23,8] <= 0
 leMC2: - lDX[24,8] + lDW[23,4] + lDW[23,8] >= 0
 leMC2: lDX[24,8] - lDW[23,4] + lDW[23,8] >= 0
 leMC2: lDX[24,8] + lDW[23,4] - lDW[23,8] >= 0
 leMC2: lDX[24,8] - lDW[23,4] - lDW[23,8] <= 0
 leMC3: ldx[24,12] - ldw[23,0] >= 0
 leMC3: ldx[24,12] - ldw[23,8] >= 0
 leMC3: ldx[24,12] - ldw[23,0] - ldw[23,8] <= 0
 leMC3: - lDX[24,12] + lDW[23,0] + lDW[23,8] >= 0
 leMC3: lDX[24,12] - lDW[23,0] + lDW[23,8] >= 0
 leMC3: lDX[24,12] + lDW[23,0] - lDW[23,8] >= 0
 leMC3: lDX[24,12] - lDW[23,0] - lDW[23,8] <= 0
 leMC0: ldx[24,1] - ldw[23,13] >= 0
 leMC0: ldx[24,1] - ldx[24,13] >= 0
 leMC0: ldx[24,1] - ldx[24,13] - ldw[23,13] <= 0
 leMC0: - lDX[24,1] + lDX[24,13] + lDW[23,13] >= 0
 leMC0: lDX[24,1] + lDX[24,13] - lDW[23,13] >= 0
 leMC0: lDX[24,1] - lDX[24,13] + lDW[23,13] >= 0
 leMC0: lDX[24,1] - lDX[24,13] - lDW[23,13] <= 0
 leMC1: ldx[24,5] - ldw[23,1] = 0
 leMC1: lDX[24,5] - lDW[23,1] = 0
 leMC2: ldx[24,9] - ldw[23,5] >= 0
 leMC2: ldx[24,9] - ldw[23,9] >= 0
 leMC2: ldx[24,9] - ldw[23,5] - ldw[23,9] <= 0
 leMC2: - lDX[24,9] + lDW[23,5] + lDW[23,9] >= 0
 leMC2: lDX[24,9] - lDW[23,5] + lDW[23,9] >= 0
 leMC2: lDX[24,9] + lDW[23,5] - lDW[23,9] >= 0
 leMC2: lDX[24,9] - lDW[23,5] - lDW[23,9] <= 0
 leMC3: ldx[24,13] - ldw[23,1] >= 0
 leMC3: ldx[24,13] - ldw[23,9] >= 0
 leMC3: ldx[24,13] - ldw[23,1] - ldw[23,9] <= 0
 leMC3: - lDX[24,13] + lDW[23,1] + lDW[23,9] >= 0
 leMC3: lDX[24,13] - lDW[23,1] + lDW[23,9] >= 0
 leMC3: lDX[24,13] + lDW[23,1] - lDW[23,9] >= 0
 leMC3: lDX[24,13] - lDW[23,1] - lDW[23,9] <= 0
 leMC0: ldx[24,2] - ldw[23,14] >= 0
 leMC0: ldx[24,2] - ldx[24,14] >= 0
 leMC0: ldx[24,2] - ldx[24,14] - ldw[23,14] <= 0
 leMC0: - lDX[24,2] + lDX[24,14] + lDW[23,14] >= 0
 leMC0: lDX[24,2] + lDX[24,14] - lDW[23,14] >= 0
 leMC0: lDX[24,2] - lDX[24,14] + lDW[23,14] >= 0
 leMC0: lDX[24,2] - lDX[24,14] - lDW[23,14] <= 0
 leMC1: ldx[24,6] - ldw[23,2] = 0
 leMC1: lDX[24,6] - lDW[23,2] = 0
 leMC2: ldx[24,10] - ldw[23,6] >= 0
 leMC2: ldx[24,10] - ldw[23,10] >= 0
 leMC2: ldx[24,10] - ldw[23,6] - ldw[23,10] <= 0
 leMC2: - lDX[24,10] + lDW[23,6] + lDW[23,10] >= 0
 leMC2: lDX[24,10] - lDW[23,6] + lDW[23,10] >= 0
 leMC2: lDX[24,10] + lDW[23,6] - lDW[23,10] >= 0
 leMC2: lDX[24,10] - lDW[23,6] - lDW[23,10] <= 0
 leMC3: ldx[24,14] - ldw[23,2] >= 0
 leMC3: ldx[24,14] - ldw[23,10] >= 0
 leMC3: ldx[24,14] - ldw[23,2] - ldw[23,10] <= 0
 leMC3: - lDX[24,14] + lDW[23,2] + lDW[23,10] >= 0
 leMC3: lDX[24,14] - lDW[23,2] + lDW[23,10] >= 0
 leMC3: lDX[24,14] + lDW[23,2] - lDW[23,10] >= 0
 leMC3: lDX[24,14] - lDW[23,2] - lDW[23,10] <= 0
 leMC0: ldx[24,3] - ldw[23,15] >= 0
 leMC0: ldx[24,3] - ldx[24,15] >= 0
 leMC0: ldx[24,3] - ldx[24,15] - ldw[23,15] <= 0
 leMC0: - lDX[24,3] + lDX[24,15] + lDW[23,15] >= 0
 leMC0: lDX[24,3] + lDX[24,15] - lDW[23,15] >= 0
 leMC0: lDX[24,3] - lDX[24,15] + lDW[23,15] >= 0
 leMC0: lDX[24,3] - lDX[24,15] - lDW[23,15] <= 0
 leMC1: ldx[24,7] - ldw[23,3] = 0
 leMC1: lDX[24,7] - lDW[23,3] = 0
 leMC2: ldx[24,11] - ldw[23,7] >= 0
 leMC2: ldx[24,11] - ldw[23,11] >= 0
 leMC2: ldx[24,11] - ldw[23,7] - ldw[23,11] <= 0
 leMC2: - lDX[24,11] + lDW[23,7] + lDW[23,11] >= 0
 leMC2: lDX[24,11] - lDW[23,7] + lDW[23,11] >= 0
 leMC2: lDX[24,11] + lDW[23,7] - lDW[23,11] >= 0
 leMC2: lDX[24,11] - lDW[23,7] - lDW[23,11] <= 0
 leMC3: ldx[24,15] - ldw[23,3] >= 0
 leMC3: ldx[24,15] - ldw[23,11] >= 0
 leMC3: ldx[24,15] - ldw[23,3] - ldw[23,11] <= 0
 leMC3: - lDX[24,15] + lDW[23,3] + lDW[23,11] >= 0
 leMC3: lDX[24,15] - lDW[23,3] + lDW[23,11] >= 0
 leMC3: lDX[24,15] + lDW[23,3] - lDW[23,11] >= 0
 leMC3: lDX[24,15] - lDW[23,3] - lDW[23,11] <= 0
 lKS: - lstk[39,0] + lLANE[0] >= 0
 lKS: - lstk[40,8] + lLANE[0] >= 0
 lKS: - lstk[41,2] + lLANE[0] >= 0
 lKS: - lstk[42,10] + lLANE[0] >= 0
 lKS: - lstk[43,4] + lLANE[0] >= 0
 lKS: - lstk[44,12] + lLANE[0] >= 0
 lKS: - lstk[45,6] + lLANE[0] >= 0
 lKS: - lstk[46,14] + lLANE[0] >= 0
 lKS: - lstk[47,5] + lLANE[0] >= 0
 lKS: - lstk[48,13] + lLANE[0] >= 0
 lKS: - lstk[49,3] + lLANE[0] >= 0
 lKS: - lstk[50,11] + lLANE[0] >= 0
 lKS: lstk[39,0] + lstk[40,8] + lstk[41,2] + lstk[42,10] + lstk[43,4]
   + lstk[44,12] + lstk[45,6] + lstk[46,14] + lstk[47,5] + lstk[48,13]
   + lstk[49,3] + lstk[50,11] - 12 lLANE[0] >= -1
 lKS: - lstk[39,1] + lLANE[1] >= 0
 lKS: - lstk[40,9] + lLANE[1] >= 0
 lKS: - lstk[41,0] + lLANE[1] >= 0
 lKS: - lstk[42,8] + lLANE[1] >= 0
 lKS: - lstk[43,2] + lLANE[1] >= 0
 lKS: - lstk[44,10] + lLANE[1] >= 0
 lKS: - lstk[45,4] + lLANE[1] >= 0
 lKS: - lstk[46,12] + lLANE[1] >= 0
 lKS: - lstk[47,6] + lLANE[1] >= 0
 lKS: - lstk[48,14] + lLANE[1] >= 0
 lKS: - lstk[49,5] + lLANE[1] >= 0
 lKS: - lstk[50,13] + lLANE[1] >= 0
 lKS: lstk[39,1] + lstk[40,9] + lstk[41,0] + lstk[42,8] + lstk[43,2]
   + lstk[44,10] + lstk[45,4] + lstk[46,12] + lstk[47,6] + lstk[48,14]
   + lstk[49,5] + lstk[50,13] - 12 lLANE[1] >= -1
 lKS: - lstk[39,2] + lLANE[2] >= 0
 lKS: - lstk[40,10] + lLANE[2] >= 0
 lKS: - lstk[41,4] + lLANE[2] >= 0
 lKS: - lstk[42,12] + lLANE[2] >= 0
 lKS: - lstk[43,6] + lLANE[2] >= 0
 lKS: - lstk[44,14] + lLANE[2] >= 0
 lKS: - lstk[45,5] + lLANE[2] >= 0
 lKS: - lstk[46,13] + lLANE[2] >= 0
 lKS: - lstk[47,3] + lLANE[2] >= 0
 lKS: - lstk[48,11] + lLANE[2] >= 0
 lKS: - lstk[49,7] + lLANE[2] >= 0
 lKS: - lstk[50,15] + lLANE[2] >= 0
 lKS: lstk[39,2] + lstk[40,10] + lstk[41,4] + lstk[42,12] + lstk[43,6]
   + lstk[44,14] + lstk[45,5] + lstk[46,13] + lstk[47,3] + lstk[48,11]
   + lstk[49,7] + lstk[50,15] - 12 lLANE[2] >= -1
 lKS: - lstk[39,3] + lLANE[3] >= 0
 lKS: - lstk[40,11] + lLANE[3] >= 0
 lKS: - lstk[41,7] + lLANE[3] >= 0
 lKS: - lstk[42,15] + lLANE[3] >= 0
 lKS: - lstk[43,1] + lLANE[3] >= 0
 lKS: - lstk[44,9] + lLANE[3] >= 0
 lKS: - lstk[45,0] + lLANE[3] >= 0
 lKS: - lstk[46,8] + lLANE[3] >= 0
 lKS: - lstk[47,2] + lLANE[3] >= 0
 lKS: - lstk[48,10] + lLANE[3] >= 0
 lKS: - lstk[49,4] + lLANE[3] >= 0
 lKS: - lstk[50,12] + lLANE[3] >= 0
 lKS: lstk[39,3] + lstk[40,11] + lstk[41,7] + lstk[42,15] + lstk[43,1]
   + lstk[44,9] + lstk[45,0] + lstk[46,8] + lstk[47,2] + lstk[48,10]
   + lstk[49,4] + lstk[50,12] - 12 lLANE[3] >= -1
 lKS: - lstk[39,4] + lLANE[4] >= 0
 lKS: - lstk[40,12] + lLANE[4] >= 0
 lKS: - lstk[41,6] + lLANE[4] >= 0
 lKS: - lstk[42,14] + lLANE[4] >= 0
 lKS: - lstk[43,5] + lLANE[4] >= 0
 lKS: - lstk[44,13] + lLANE[4] >= 0
 lKS: - lstk[45,3] + lLANE[4] >= 0
 lKS: - lstk[46,11] + lLANE[4] >= 0
 lKS: - lstk[47,7] + lLANE[4] >= 0
 lKS: - lstk[48,15] + lLANE[4] >= 0
 lKS: - lstk[49,1] + lLANE[4] >= 0
 lKS: - lstk[50,9] + lLANE[4] >= 0
 lKS: lstk[39,4] + lstk[40,12] + lstk[41,6] + lstk[42,14] + lstk[43,5]
   + lstk[44,13] + lstk[45,3] + lstk[46,11] + lstk[47,7] + lstk[48,15]
   + lstk[49,1] + lstk[50,9] - 12 lLANE[4] >= -1
 lKS: - lstk[39,5] + lLANE[5] >= 0
 lKS: - lstk[40,13] + lLANE[5] >= 0
 lKS: - lstk[41,3] + lLANE[5] >= 0
 lKS: - lstk[42,11] + lLANE[5] >= 0
 lKS: - lstk[43,7] + lLANE[5] >= 0
 lKS: - lstk[44,15] + lLANE[5] >= 0
 lKS: - lstk[45,1] + lLANE[5] >= 0
 lKS: - lstk[46,9] + lLANE[5] >= 0
 lKS: - lstk[47,0] + lLANE[5] >= 0
 lKS: - lstk[48,8] + lLANE[5] >= 0
 lKS: - lstk[49,2] + lLANE[5] >= 0
 lKS: - lstk[50,10] + lLANE[5] >= 0
 lKS: lstk[39,5] + lstk[40,13] + lstk[41,3] + lstk[42,11] + lstk[43,7]
   + lstk[44,15] + lstk[45,1] + lstk[46,9] + lstk[47,0] + lstk[48,8]
   + lstk[49,2] + lstk[50,10] - 12 lLANE[5] >= -1
 lKS: - lstk[39,6] + lLANE[6] >= 0
 lKS: - lstk[40,14] + lLANE[6] >= 0
 lKS: - lstk[41,5] + lLANE[6] >= 0
 lKS: - lstk[42,13] + lLANE[6] >= 0
 lKS: - lstk[43,3] + lLANE[6] >= 0
 lKS: - lstk[44,11] + lLANE[6] >= 0
 lKS: - lstk[45,7] + lLANE[6] >= 0
 lKS: - lstk[46,15] + lLANE[6] >= 0
 lKS: - lstk[47,1] + lLANE[6] >= 0
 lKS: - lstk[48,9] + lLANE[6] >= 0
 lKS: - lstk[49,0] + lLANE[6] >= 0
 lKS: - lstk[50,8] + lLANE[6] >= 0
 lKS: lstk[39,6] + lstk[40,14] + lstk[41,5] + lstk[42,13] + lstk[43,3]
   + lstk[44,11] + lstk[45,7] + lstk[46,15] + lstk[47,1] + lstk[48,9]
   + lstk[49,0] + lstk[50,8] - 12 lLANE[6] >= -1
 lKS: - lstk[39,7] + lLANE[7] >= 0
 lKS: - lstk[40,15] + lLANE[7] >= 0
 lKS: - lstk[41,1] + lLANE[7] >= 0
 lKS: - lstk[42,9] + lLANE[7] >= 0
 lKS: - lstk[43,0] + lLANE[7] >= 0
 lKS: - lstk[44,8] + lLANE[7] >= 0
 lKS: - lstk[45,2] + lLANE[7] >= 0
 lKS: - lstk[46,10] + lLANE[7] >= 0
 lKS: - lstk[47,4] + lLANE[7] >= 0
 lKS: - lstk[48,12] + lLANE[7] >= 0
 lKS: - lstk[49,6] + lLANE[7] >= 0
 lKS: - lstk[50,14] + lLANE[7] >= 0
 lKS: lstk[39,7] + lstk[40,15] + lstk[41,1] + lstk[42,9] + lstk[43,0]
   + lstk[44,8] + lstk[45,2] + lstk[46,10] + lstk[47,4] + lstk[48,12]
   + lstk[49,6] + lstk[50,14] - 12 lLANE[7] >= -1
 lKS: - lstk[39,8] + lLANE[8] >= 0
 lKS: - lstk[40,2] + lLANE[8] >= 0
 lKS: - lstk[41,10] + lLANE[8] >= 0
 lKS: - lstk[42,4] + lLANE[8] >= 0
 lKS: - lstk[43,12] + lLANE[8] >= 0
 lKS: - lstk[44,6] + lLANE[8] >= 0
 lKS: - lstk[45,14] + lLANE[8] >= 0
 lKS: - lstk[46,5] + lLANE[8] >= 0
 lKS: - lstk[47,13] + lLANE[8] >= 0
 lKS: - lstk[48,3] + lLANE[8] >= 0
 lKS: - lstk[49,11] + lLANE[8] >= 0
 lKS: - lstk[50,7] + lLANE[8] >= 0
 lKS: lstk[39,8] + lstk[40,2] + lstk[41,10] + lstk[42,4] + lstk[43,12]
   + lstk[44,6] + lstk[45,14] + lstk[46,5] + lstk[47,13] + lstk[48,3]
   + lstk[49,11] + lstk[50,7] - 12 lLANE[8] >= -1
 lKS: - lstk[39,9] + lLANE[9] >= 0
 lKS: - lstk[40,0] + lLANE[9] >= 0
 lKS: - lstk[41,8] + lLANE[9] >= 0
 lKS: - lstk[42,2] + lLANE[9] >= 0
 lKS: - lstk[43,10] + lLANE[9] >= 0
 lKS: - lstk[44,4] + lLANE[9] >= 0
 lKS: - lstk[45,12] + lLANE[9] >= 0
 lKS: - lstk[46,6] + lLANE[9] >= 0
 lKS: - lstk[47,14] + lLANE[9] >= 0
 lKS: - lstk[48,5] + lLANE[9] >= 0
 lKS: - lstk[49,13] + lLANE[9] >= 0
 lKS: - lstk[50,3] + lLANE[9] >= 0
 lKS: lstk[39,9] + lstk[40,0] + lstk[41,8] + lstk[42,2] + lstk[43,10]
   + lstk[44,4] + lstk[45,12] + lstk[46,6] + lstk[47,14] + lstk[48,5]
   + lstk[49,13] + lstk[50,3] - 12 lLANE[9] >= -1
 lKS: - lstk[39,10] + lLANE[10] >= 0
 lKS: - lstk[40,4] + lLANE[10] >= 0
 lKS: - lstk[41,12] + lLANE[10] >= 0
 lKS: - lstk[42,6] + lLANE[10] >= 0
 lKS: - lstk[43,14] + lLANE[10] >= 0
 lKS: - lstk[44,5] + lLANE[10] >= 0
 lKS: - lstk[45,13] + lLANE[10] >= 0
 lKS: - lstk[46,3] + lLANE[10] >= 0
 lKS: - lstk[47,11] + lLANE[10] >= 0
 lKS: - lstk[48,7] + lLANE[10] >= 0
 lKS: - lstk[49,15] + lLANE[10] >= 0
 lKS: - lstk[50,1] + lLANE[10] >= 0
 lKS: lstk[39,10] + lstk[40,4] + lstk[41,12] + lstk[42,6] + lstk[43,14]
   + lstk[44,5] + lstk[45,13] + lstk[46,3] + lstk[47,11] + lstk[48,7]
   + lstk[49,15] + lstk[50,1] - 12 lLANE[10] >= -1
 lKS: - lstk[39,11] + lLANE[11] >= 0
 lKS: - lstk[40,7] + lLANE[11] >= 0
 lKS: - lstk[41,15] + lLANE[11] >= 0
 lKS: - lstk[42,1] + lLANE[11] >= 0
 lKS: - lstk[43,9] + lLANE[11] >= 0
 lKS: - lstk[44,0] + lLANE[11] >= 0
 lKS: - lstk[45,8] + lLANE[11] >= 0
 lKS: - lstk[46,2] + lLANE[11] >= 0
 lKS: - lstk[47,10] + lLANE[11] >= 0
 lKS: - lstk[48,4] + lLANE[11] >= 0
 lKS: - lstk[49,12] + lLANE[11] >= 0
 lKS: - lstk[50,6] + lLANE[11] >= 0
 lKS: lstk[39,11] + lstk[40,7] + lstk[41,15] + lstk[42,1] + lstk[43,9]
   + lstk[44,0] + lstk[45,8] + lstk[46,2] + lstk[47,10] + lstk[48,4]
   + lstk[49,12] + lstk[50,6] - 12 lLANE[11] >= -1
 lKS: - lstk[39,12] + lLANE[12] >= 0
 lKS: - lstk[40,6] + lLANE[12] >= 0
 lKS: - lstk[41,14] + lLANE[12] >= 0
 lKS: - lstk[42,5] + lLANE[12] >= 0
 lKS: - lstk[43,13] + lLANE[12] >= 0
 lKS: - lstk[44,3] + lLANE[12] >= 0
 lKS: - lstk[45,11] + lLANE[12] >= 0
 lKS: - lstk[46,7] + lLANE[12] >= 0
 lKS: - lstk[47,15] + lLANE[12] >= 0
 lKS: - lstk[48,1] + lLANE[12] >= 0
 lKS: - lstk[49,9] + lLANE[12] >= 0
 lKS: - lstk[50,0] + lLANE[12] >= 0
 lKS: lstk[39,12] + lstk[40,6] + lstk[41,14] + lstk[42,5] + lstk[43,13]
   + lstk[44,3] + lstk[45,11] + lstk[46,7] + lstk[47,15] + lstk[48,1]
   + lstk[49,9] + lstk[50,0] - 12 lLANE[12] >= -1
 lKS: - lstk[39,13] + lLANE[13] >= 0
 lKS: - lstk[40,3] + lLANE[13] >= 0
 lKS: - lstk[41,11] + lLANE[13] >= 0
 lKS: - lstk[42,7] + lLANE[13] >= 0
 lKS: - lstk[43,15] + lLANE[13] >= 0
 lKS: - lstk[44,1] + lLANE[13] >= 0
 lKS: - lstk[45,9] + lLANE[13] >= 0
 lKS: - lstk[46,0] + lLANE[13] >= 0
 lKS: - lstk[47,8] + lLANE[13] >= 0
 lKS: - lstk[48,2] + lLANE[13] >= 0
 lKS: - lstk[49,10] + lLANE[13] >= 0
 lKS: - lstk[50,4] + lLANE[13] >= 0
 lKS: lstk[39,13] + lstk[40,3] + lstk[41,11] + lstk[42,7] + lstk[43,15]
   + lstk[44,1] + lstk[45,9] + lstk[46,0] + lstk[47,8] + lstk[48,2]
   + lstk[49,10] + lstk[50,4] - 12 lLANE[13] >= -1
 lKS: - lstk[39,14] + lLANE[14] >= 0
 lKS: - lstk[40,5] + lLANE[14] >= 0
 lKS: - lstk[41,13] + lLANE[14] >= 0
 lKS: - lstk[42,3] + lLANE[14] >= 0
 lKS: - lstk[43,11] + lLANE[14] >= 0
 lKS: - lstk[44,7] + lLANE[14] >= 0
 lKS: - lstk[45,15] + lLANE[14] >= 0
 lKS: - lstk[46,1] + lLANE[14] >= 0
 lKS: - lstk[47,9] + lLANE[14] >= 0
 lKS: - lstk[48,0] + lLANE[14] >= 0
 lKS: - lstk[49,8] + lLANE[14] >= 0
 lKS: - lstk[50,2] + lLANE[14] >= 0
 lKS: lstk[39,14] + lstk[40,5] + lstk[41,13] + lstk[42,3] + lstk[43,11]
   + lstk[44,7] + lstk[45,15] + lstk[46,1] + lstk[47,9] + lstk[48,0]
   + lstk[49,8] + lstk[50,2] - 12 lLANE[14] >= -1
 lKS: - lstk[39,15] + lLANE[15] >= 0
 lKS: - lstk[40,1] + lLANE[15] >= 0
 lKS: - lstk[41,9] + lLANE[15] >= 0
 lKS: - lstk[42,0] + lLANE[15] >= 0
 lKS: - lstk[43,8] + lLANE[15] >= 0
 lKS: - lstk[44,2] + lLANE[15] >= 0
 lKS: - lstk[45,10] + lLANE[15] >= 0
 lKS: - lstk[46,4] + lLANE[15] >= 0
 lKS: - lstk[47,12] + lLANE[15] >= 0
 lKS: - lstk[48,6] + lLANE[15] >= 0
 lKS: - lstk[49,14] + lLANE[15] >= 0
 lKS: - lstk[50,5] + lLANE[15] >= 0
 lKS: lstk[39,15] + lstk[40,1] + lstk[41,9] + lstk[42,0] + lstk[43,8]
   + lstk[44,2] + lstk[45,10] + lstk[46,4] + lstk[47,12] + lstk[48,6]
   + lstk[49,14] + lstk[50,5] - 12 lLANE[15] >= -1
 uGstkL: uGstk[0,8] = 0
 uGstkL: uGstk[0,9] = 0
 uGstkL: uGstk[0,10] = 0
 uGstkL: uGstk[0,11] = 0
 uGstkL: uGstk[0,12] = 0
 uGstkL: uGstk[0,13] = 0
 uGstkL: uGstk[0,14] = 0
 uGstkL: uGstk[0,15] = 0
 uGstkL: uGstk[1,8] = 0
 uGstkL: uGstk[1,9] = 0
 uGstkL: uGstk[1,10] = 0
 uGstkL: uGstk[1,11] = 0
 uGstkL: uGstk[1,12] = 0
 uGstkL: uGstk[1,13] = 0
 uGstkL: uGstk[1,14] = 0
 uGstkL: uGstk[1,15] = 0
 uGstkL: uGstk[2,8] = 0
 uGstkL: uGstk[2,9] = 0
 uGstkL: uGstk[2,10] = 0
 uGstkL: uGstk[2,11] = 0
 uGstkL: uGstk[2,12] = 0
 uGstkL: uGstk[2,13] = 0
 uGstkL: uGstk[2,14] = 0
 uGstkL: uGstk[2,15] = 0
 uGstkL: uGstk[3,8] = 0
 uGstkL: uGstk[3,9] = 0
 uGstkL: uGstk[3,10] = 0
 uGstkL: uGstk[3,11] = 0
 uGstkL: uGstk[3,12] = 0
 uGstkL: uGstk[3,13] = 0
 uGstkL: uGstk[3,14] = 0
 uGstkL: uGstk[3,15] = 0
 ugstkLa: ugstka[0,8] = 0
 ugstkLa: ugstka[0,9] = 0
 ugstkLa: ugstka[0,10] = 0
 ugstkLa: ugstka[0,11] = 0
 ugstkLa: ugstka[0,12] = 0
 ugstkLa: ugstka[0,13] = 0
 ugstkLa: ugstka[0,14] = 0
 ugstkLa: ugstka[0,15] = 0
 ugstkLa: ugstka[1,8] = 0
 ugstkLa: ugstka[1,9] = 0
 ugstkLa: ugstka[1,10] = 0
 ugstkLa: ugstka[1,11] = 0
 ugstkLa: ugstka[1,12] = 0
 ugstkLa: ugstka[1,13] = 0
 ugstkLa: ugstka[1,14] = 0
 ugstkLa: ugstka[1,15] = 0
 ugstkLa: ugstka[2,8] = 0
 ugstkLa: ugstka[2,9] = 0
 ugstkLa: ugstka[2,10] = 0
 ugstkLa: ugstka[2,11] = 0
 ugstkLa: ugstka[2,12] = 0
 ugstkLa: ugstka[2,13] = 0
 ugstkLa: ugstka[2,14] = 0
 ugstkLa: ugstka[2,15] = 0
 ugstkLa: ugstka[3,8] = 0
 ugstkLa: ugstka[3,9] = 0
 ugstkLa: ugstka[3,10] = 0
 ugstkLa: ugstka[3,11] = 0
 ugstkLa: ugstka[3,12] = 0
 ugstkLa: ugstka[3,13] = 0
 ugstkLa: ugstka[3,14] = 0
 ugstkLa: ugstka[3,15] = 0
 VP: uDetW[0,0] + uDetW[0,1] + uDetW[0,2] + uDetW[0,3] + uDetW[0,4]
   + uDetW[0,5] + uDetW[0,6] + uDetW[0,7] + uDetW[0,8] + uDetW[0,9]
   + uDetW[0,10] + uDetW[0,11] + uDetW[0,12] + uDetW[0,13] + uDetW[0,14]
   + uDetW[0,15] = 16
 W0X0a: uDetX[1,0] - uDetW[0,0] = 0
 W0X0a: uDetX[1,1] - uDetW[0,1] = 0
 W0X0a: uDetX[1,2] - uDetW[0,2] = 0
 W0X0a: uDetX[1,3] - uDetW[0,3] = 0
 W0X0a: uDetX[1,4] - uDetW[0,4] = 0
 W0X0a: uDetX[1,5] - uDetW[0,5] = 0
 W0X0a: uDetX[1,6] - uDetW[0,6] = 0
 W0X0a: uDetX[1,7] - uDetW[0,7] = 0
 W0X0a: uDetX[1,8] - uDetW[0,8] = 0
 W0X0a: uDetX[1,9] - uDetW[0,9] = 0
 W0X0a: uDetX[1,10] - uDetW[0,10] = 0
 W0X0a: uDetX[1,11] - uDetW[0,11] = 0
 W0X0a: uDetX[1,12] - uDetW[0,12] = 0
 W0X0a: uDetX[1,13] - uDetW[0,13] = 0
 W0X0a: uDetX[1,14] - uDetW[0,14] = 0
 W0X0a: uDetX[1,15] - uDetW[0,15] = 0
 VPa: udetWa[0,0] + udetWa[0,1] + udetWa[0,2] + udetWa[0,3] + udetWa[0,4]
   + udetWa[0,5] + udetWa[0,6] + udetWa[0,7] + udetWa[0,8] + udetWa[0,9]
   + udetWa[0,10] + udetWa[0,11] + udetWa[0,12] + udetWa[0,13]
   + udetWa[0,14] + udetWa[0,15] = 16
 W0X0a: - uDetW[0,0] + udetXa[1,0] = 0
 W0X0a: - uDetW[0,1] + udetXa[1,1] = 0
 W0X0a: - uDetW[0,2] + udetXa[1,2] = 0
 W0X0a: - uDetW[0,3] + udetXa[1,3] = 0
 W0X0a: - uDetW[0,4] + udetXa[1,4] = 0
 W0X0a: - uDetW[0,5] + udetXa[1,5] = 0
 W0X0a: - uDetW[0,6] + udetXa[1,6] = 0
 W0X0a: - uDetW[0,7] + udetXa[1,7] = 0
 W0X0a: - uDetW[0,8] + udetXa[1,8] = 0
 W0X0a: - uDetW[0,9] + udetXa[1,9] = 0
 W0X0a: - uDetW[0,10] + udetXa[1,10] = 0
 W0X0a: - uDetW[0,11] + udetXa[1,11] = 0
 W0X0a: - uDetW[0,12] + udetXa[1,12] = 0
 W0X0a: - uDetW[0,13] + udetXa[1,13] = 0
 W0X0a: - uDetW[0,14] + udetXa[1,14] = 0
 W0X0a: - uDetW[0,15] + udetXa[1,15] = 0
 geqk0: ugEQK[0,0] - uGstk[0,0] = 0
 geqk0A: ugeqka[0,0] - ugstka[0,0] = 0
 geqk1: ugEQK[0,4] - uGstk[0,0] = 0
 geqk1A: ugeqka[0,4] - ugstka[0,0] = 0
 geqk2: ugEQK[0,8] - uGstk[0,7] = 0
 geqk2A: ugeqka[0,8] - ugstka[0,7] = 0
 geqk3: ugEQK[0,12] - uGstk[0,0] = 0
 geqk3A: ugeqka[0,12] - ugstka[0,0] = 0
 geqk0: ugEQK[0,1] - uGstk[0,1] = 0
 geqk0A: ugeqka[0,1] - ugstka[0,1] = 0
 geqk1: ugEQK[0,5] - uGstk[0,1] = 0
 geqk1A: ugeqka[0,5] - ugstka[0,1] = 0
 geqk2: ugEQK[0,9] - uGstk[0,4] = 0
 geqk2A: ugeqka[0,9] - ugstka[0,4] = 0
 geqk3: ugEQK[0,13] - uGstk[0,1] = 0
 geqk3A: ugeqka[0,13] - ugstka[0,1] = 0
 geqk0: ugEQK[0,2] - uGstk[0,2] = 0
 geqk0A: ugeqka[0,2] - ugstka[0,2] = 0
 geqk1: ugEQK[0,6] - uGstk[0,2] = 0
 geqk1A: ugeqka[0,6] - ugstka[0,2] = 0
 geqk2: ugEQK[0,10] - uGstk[0,5] = 0
 geqk2A: ugeqka[0,10] - ugstka[0,5] = 0
 geqk3: ugEQK[0,14] - uGstk[0,2] = 0
 geqk3A: ugeqka[0,14] - ugstka[0,2] = 0
 geqk0: ugEQK[0,3] - uGstk[0,3] = 0
 geqk0A: ugeqka[0,3] - ugstka[0,3] = 0
 geqk1: ugEQK[0,7] - uGstk[0,3] = 0
 geqk1A: ugeqka[0,7] - ugstka[0,3] = 0
 geqk2: ugEQK[0,11] - uGstk[0,6] = 0
 geqk2A: ugeqka[0,11] - ugstka[0,6] = 0
 geqk3: ugEQK[0,15] - uGstk[0,3] = 0
 geqk3A: ugeqka[0,15] - ugstka[0,3] = 0
 geqk0: ugEQK[1,0] - uGstk[1,0] = 0
 geqk0A: ugeqka[1,0] - ugstka[1,0] = 0
 geqk1: ugEQK[1,4] - uGstk[1,0] = 0
 geqk1A: ugeqka[1,4] - ugstka[1,0] = 0
 geqk2: ugEQK[1,8] - uGstk[1,7] = 0
 geqk2A: ugeqka[1,8] - ugstka[1,7] = 0
 geqk3: ugEQK[1,12] - uGstk[1,0] = 0
 geqk3A: ugeqka[1,12] - ugstka[1,0] = 0
 geqk0: ugEQK[1,1] - uGstk[1,1] = 0
 geqk0A: ugeqka[1,1] - ugstka[1,1] = 0
 geqk1: ugEQK[1,5] - uGstk[1,1] = 0
 geqk1A: ugeqka[1,5] - ugstka[1,1] = 0
 geqk2: ugEQK[1,9] - uGstk[1,4] = 0
 geqk2A: ugeqka[1,9] - ugstka[1,4] = 0
 geqk3: ugEQK[1,13] - uGstk[1,1] = 0
 geqk3A: ugeqka[1,13] - ugstka[1,1] = 0
 geqk0: ugEQK[1,2] - uGstk[1,2] = 0
 geqk0A: ugeqka[1,2] - ugstka[1,2] = 0
 geqk1: ugEQK[1,6] - uGstk[1,2] = 0
 geqk1A: ugeqka[1,6] - ugstka[1,2] = 0
 geqk2: ugEQK[1,10] - uGstk[1,5] = 0
 geqk2A: ugeqka[1,10] - ugstka[1,5] = 0
 geqk3: ugEQK[1,14] - uGstk[1,2] = 0
 geqk3A: ugeqka[1,14] - ugstka[1,2] = 0
 geqk0: ugEQK[1,3] - uGstk[1,3] = 0
 geqk0A: ugeqka[1,3] - ugstka[1,3] = 0
 geqk1: ugEQK[1,7] - uGstk[1,3] = 0
 geqk1A: ugeqka[1,7] - ugstka[1,3] = 0
 geqk2: ugEQK[1,11] - uGstk[1,6] = 0
 geqk2A: ugeqka[1,11] - ugstka[1,6] = 0
 geqk3: ugEQK[1,15] - uGstk[1,3] = 0
 geqk3A: ugeqka[1,15] - ugstka[1,3] = 0
 geqk0: ugEQK[2,0] - uGstk[2,0] = 0
 geqk0A: ugeqka[2,0] - ugstka[2,0] = 0
 geqk1: ugEQK[2,4] - uGstk[2,0] = 0
 geqk1A: ugeqka[2,4] - ugstka[2,0] = 0
 geqk2: ugEQK[2,8] - uGstk[2,7] = 0
 geqk2A: ugeqka[2,8] - ugstka[2,7] = 0
 geqk3: ugEQK[2,12] - uGstk[2,0] = 0
 geqk3A: ugeqka[2,12] - ugstka[2,0] = 0
 geqk0: ugEQK[2,1] - uGstk[2,1] = 0
 geqk0A: ugeqka[2,1] - ugstka[2,1] = 0
 geqk1: ugEQK[2,5] - uGstk[2,1] = 0
 geqk1A: ugeqka[2,5] - ugstka[2,1] = 0
 geqk2: ugEQK[2,9] - uGstk[2,4] = 0
 geqk2A: ugeqka[2,9] - ugstka[2,4] = 0
 geqk3: ugEQK[2,13] - uGstk[2,1] = 0
 geqk3A: ugeqka[2,13] - ugstka[2,1] = 0
 geqk0: ugEQK[2,2] - uGstk[2,2] = 0
 geqk0A: ugeqka[2,2] - ugstka[2,2] = 0
 geqk1: ugEQK[2,6] - uGstk[2,2] = 0
 geqk1A: ugeqka[2,6] - ugstka[2,2] = 0
 geqk2: ugEQK[2,10] - uGstk[2,5] = 0
 geqk2A: ugeqka[2,10] - ugstka[2,5] = 0
 geqk3: ugEQK[2,14] - uGstk[2,2] = 0
 geqk3A: ugeqka[2,14] - ugstka[2,2] = 0
 geqk0: ugEQK[2,3] - uGstk[2,3] = 0
 geqk0A: ugeqka[2,3] - ugstka[2,3] = 0
 geqk1: ugEQK[2,7] - uGstk[2,3] = 0
 geqk1A: ugeqka[2,7] - ugstka[2,3] = 0
 geqk2: ugEQK[2,11] - uGstk[2,6] = 0
 geqk2A: ugeqka[2,11] - ugstka[2,6] = 0
 geqk3: ugEQK[2,15] - uGstk[2,3] = 0
 geqk3A: ugeqka[2,15] - ugstka[2,3] = 0
 geqk0: ugEQK[3,0] - uGstk[3,0] = 0
 geqk0A: ugeqka[3,0] - ugstka[3,0] = 0
 geqk1: ugEQK[3,4] - uGstk[3,0] = 0
 geqk1A: ugeqka[3,4] - ugstka[3,0] = 0
 geqk2: ugEQK[3,8] - uGstk[3,7] = 0
 geqk2A: ugeqka[3,8] - ugstka[3,7] = 0
 geqk3: ugEQK[3,12] - uGstk[3,0] = 0
 geqk3A: ugeqka[3,12] - ugstka[3,0] = 0
 geqk0: ugEQK[3,1] - uGstk[3,1] = 0
 geqk0A: ugeqka[3,1] - ugstka[3,1] = 0
 geqk1: ugEQK[3,5] - uGstk[3,1] = 0
 geqk1A: ugeqka[3,5] - ugstka[3,1] = 0
 geqk2: ugEQK[3,9] - uGstk[3,4] = 0
 geqk2A: ugeqka[3,9] - ugstka[3,4] = 0
 geqk3: ugEQK[3,13] - uGstk[3,1] = 0
 geqk3A: ugeqka[3,13] - ugstka[3,1] = 0
 geqk0: ugEQK[3,2] - uGstk[3,2] = 0
 geqk0A: ugeqka[3,2] - ugstka[3,2] = 0
 geqk1: ugEQK[3,6] - uGstk[3,2] = 0
 geqk1A: ugeqka[3,6] - ugstka[3,2] = 0
 geqk2: ugEQK[3,10] - uGstk[3,5] = 0
 geqk2A: ugeqka[3,10] - ugstka[3,5] = 0
 geqk3: ugEQK[3,14] - uGstk[3,2] = 0
 geqk3A: ugeqka[3,14] - ugstka[3,2] = 0
 geqk0: ugEQK[3,3] - uGstk[3,3] = 0
 geqk0A: ugeqka[3,3] - ugstka[3,3] = 0
 geqk1: ugEQK[3,7] - uGstk[3,3] = 0
 geqk1A: ugeqka[3,7] - ugstka[3,3] = 0
 geqk2: ugEQK[3,11] - uGstk[3,6] = 0
 geqk2A: ugeqka[3,11] - ugstka[3,6] = 0
 geqk3: ugEQK[3,15] - uGstk[3,3] = 0
 geqk3A: ugeqka[3,15] - ugstka[3,3] = 0
 udetSC+ek: - uDetX[1,0] + uDetY[1,0] <= 0
 udetSC+ek: - ugEQK[0,0] + uDetY[1,0] <= 0
 udetSC+ek: - ugEQK[0,0] - uDetX[1,0] + uDetY[1,0] >= -1
 udetSC+ekA: - udetXa[1,0] + udetYa[1,0] <= 0
 udetSC+ekA: - ugeqka[0,0] + udetYa[1,0] <= 0
 udetSC+ekA: - ugeqka[0,0] - udetXa[1,0] + udetYa[1,0] >= -1
 ufrSC: uDY[1,0] - udy[1,0] - uFrSC[1,0] >= 0
 ufrSC: uDetY[1,0] - uFrSC[1,0] >= 0
 ufrSC: - uDY[1,0] + udy[1,0] - uDetY[1,0] + uFrSC[1,0] >= -1
 ufrSCa: uDY[1,0] - udy[1,0] - ufrSCa[1,0] >= 0
 ufrSCa: udetYa[1,0] - ufrSCa[1,0] >= 0
 ufrSCa: - uDY[1,0] + udy[1,0] - udetYa[1,0] + ufrSCa[1,0] >= -1
 udetSC+ek: - uDetX[1,1] + uDetY[1,1] <= 0
 udetSC+ek: - ugEQK[0,1] + uDetY[1,1] <= 0
 udetSC+ek: - ugEQK[0,1] - uDetX[1,1] + uDetY[1,1] >= -1
 udetSC+ekA: - udetXa[1,1] + udetYa[1,1] <= 0
 udetSC+ekA: - ugeqka[0,1] + udetYa[1,1] <= 0
 udetSC+ekA: - ugeqka[0,1] - udetXa[1,1] + udetYa[1,1] >= -1
 ufrSC: uDY[1,1] - udy[1,1] - uFrSC[1,1] >= 0
 ufrSC: uDetY[1,1] - uFrSC[1,1] >= 0
 ufrSC: - uDY[1,1] + udy[1,1] - uDetY[1,1] + uFrSC[1,1] >= -1
 ufrSCa: uDY[1,1] - udy[1,1] - ufrSCa[1,1] >= 0
 ufrSCa: udetYa[1,1] - ufrSCa[1,1] >= 0
 ufrSCa: - uDY[1,1] + udy[1,1] - udetYa[1,1] + ufrSCa[1,1] >= -1
 udetSC+ek: - uDetX[1,2] + uDetY[1,2] <= 0
 udetSC+ek: - ugEQK[0,2] + uDetY[1,2] <= 0
 udetSC+ek: - ugEQK[0,2] - uDetX[1,2] + uDetY[1,2] >= -1
 udetSC+ekA: - udetXa[1,2] + udetYa[1,2] <= 0
 udetSC+ekA: - ugeqka[0,2] + udetYa[1,2] <= 0
 udetSC+ekA: - ugeqka[0,2] - udetXa[1,2] + udetYa[1,2] >= -1
 ufrSC: uDY[1,2] - udy[1,2] - uFrSC[1,2] >= 0
 ufrSC: uDetY[1,2] - uFrSC[1,2] >= 0
 ufrSC: - uDY[1,2] + udy[1,2] - uDetY[1,2] + uFrSC[1,2] >= -1
 ufrSCa: uDY[1,2] - udy[1,2] - ufrSCa[1,2] >= 0
 ufrSCa: udetYa[1,2] - ufrSCa[1,2] >= 0
 ufrSCa: - uDY[1,2] + udy[1,2] - udetYa[1,2] + ufrSCa[1,2] >= -1
 udetSC+ek: - uDetX[1,3] + uDetY[1,3] <= 0
 udetSC+ek: - ugEQK[0,3] + uDetY[1,3] <= 0
 udetSC+ek: - ugEQK[0,3] - uDetX[1,3] + uDetY[1,3] >= -1
 udetSC+ekA: - udetXa[1,3] + udetYa[1,3] <= 0
 udetSC+ekA: - ugeqka[0,3] + udetYa[1,3] <= 0
 udetSC+ekA: - ugeqka[0,3] - udetXa[1,3] + udetYa[1,3] >= -1
 ufrSC: uDY[1,3] - udy[1,3] - uFrSC[1,3] >= 0
 ufrSC: uDetY[1,3] - uFrSC[1,3] >= 0
 ufrSC: - uDY[1,3] + udy[1,3] - uDetY[1,3] + uFrSC[1,3] >= -1
 ufrSCa: uDY[1,3] - udy[1,3] - ufrSCa[1,3] >= 0
 ufrSCa: udetYa[1,3] - ufrSCa[1,3] >= 0
 ufrSCa: - uDY[1,3] + udy[1,3] - udetYa[1,3] + ufrSCa[1,3] >= -1
 udetSC+ek: - uDetX[1,4] + uDetY[1,4] <= 0
 udetSC+ek: - ugEQK[0,4] + uDetY[1,4] <= 0
 udetSC+ek: - ugEQK[0,4] - uDetX[1,4] + uDetY[1,4] >= -1
 udetSC+ekA: - udetXa[1,4] + udetYa[1,4] <= 0
 udetSC+ekA: - ugeqka[0,4] + udetYa[1,4] <= 0
 udetSC+ekA: - ugeqka[0,4] - udetXa[1,4] + udetYa[1,4] >= -1
 ufrSC: uDY[1,4] - udy[1,4] - uFrSC[1,4] >= 0
 ufrSC: uDetY[1,4] - uFrSC[1,4] >= 0
 ufrSC: - uDY[1,4] + udy[1,4] - uDetY[1,4] + uFrSC[1,4] >= -1
 ufrSCa: uDY[1,4] - udy[1,4] - ufrSCa[1,4] >= 0
 ufrSCa: udetYa[1,4] - ufrSCa[1,4] >= 0
 ufrSCa: - uDY[1,4] + udy[1,4] - udetYa[1,4] + ufrSCa[1,4] >= -1
 udetSC+ek: - uDetX[1,5] + uDetY[1,5] <= 0
 udetSC+ek: - ugEQK[0,5] + uDetY[1,5] <= 0
 udetSC+ek: - ugEQK[0,5] - uDetX[1,5] + uDetY[1,5] >= -1
 udetSC+ekA: - udetXa[1,5] + udetYa[1,5] <= 0
 udetSC+ekA: - ugeqka[0,5] + udetYa[1,5] <= 0
 udetSC+ekA: - ugeqka[0,5] - udetXa[1,5] + udetYa[1,5] >= -1
 ufrSC: uDY[1,5] - udy[1,5] - uFrSC[1,5] >= 0
 ufrSC: uDetY[1,5] - uFrSC[1,5] >= 0
 ufrSC: - uDY[1,5] + udy[1,5] - uDetY[1,5] + uFrSC[1,5] >= -1
 ufrSCa: uDY[1,5] - udy[1,5] - ufrSCa[1,5] >= 0
 ufrSCa: udetYa[1,5] - ufrSCa[1,5] >= 0
 ufrSCa: - uDY[1,5] + udy[1,5] - udetYa[1,5] + ufrSCa[1,5] >= -1
 udetSC+ek: - uDetX[1,6] + uDetY[1,6] <= 0
 udetSC+ek: - ugEQK[0,6] + uDetY[1,6] <= 0
 udetSC+ek: - ugEQK[0,6] - uDetX[1,6] + uDetY[1,6] >= -1
 udetSC+ekA: - udetXa[1,6] + udetYa[1,6] <= 0
 udetSC+ekA: - ugeqka[0,6] + udetYa[1,6] <= 0
 udetSC+ekA: - ugeqka[0,6] - udetXa[1,6] + udetYa[1,6] >= -1
 ufrSC: uDY[1,6] - udy[1,6] - uFrSC[1,6] >= 0
 ufrSC: uDetY[1,6] - uFrSC[1,6] >= 0
 ufrSC: - uDY[1,6] + udy[1,6] - uDetY[1,6] + uFrSC[1,6] >= -1
 ufrSCa: uDY[1,6] - udy[1,6] - ufrSCa[1,6] >= 0
 ufrSCa: udetYa[1,6] - ufrSCa[1,6] >= 0
 ufrSCa: - uDY[1,6] + udy[1,6] - udetYa[1,6] + ufrSCa[1,6] >= -1
 udetSC+ek: - uDetX[1,7] + uDetY[1,7] <= 0
 udetSC+ek: - ugEQK[0,7] + uDetY[1,7] <= 0
 udetSC+ek: - ugEQK[0,7] - uDetX[1,7] + uDetY[1,7] >= -1
 udetSC+ekA: - udetXa[1,7] + udetYa[1,7] <= 0
 udetSC+ekA: - ugeqka[0,7] + udetYa[1,7] <= 0
 udetSC+ekA: - ugeqka[0,7] - udetXa[1,7] + udetYa[1,7] >= -1
 ufrSC: uDY[1,7] - udy[1,7] - uFrSC[1,7] >= 0
 ufrSC: uDetY[1,7] - uFrSC[1,7] >= 0
 ufrSC: - uDY[1,7] + udy[1,7] - uDetY[1,7] + uFrSC[1,7] >= -1
 ufrSCa: uDY[1,7] - udy[1,7] - ufrSCa[1,7] >= 0
 ufrSCa: udetYa[1,7] - ufrSCa[1,7] >= 0
 ufrSCa: - uDY[1,7] + udy[1,7] - udetYa[1,7] + ufrSCa[1,7] >= -1
 udetSC+ek: - uDetX[1,8] + uDetY[1,8] <= 0
 udetSC+ek: - ugEQK[0,8] + uDetY[1,8] <= 0
 udetSC+ek: - ugEQK[0,8] - uDetX[1,8] + uDetY[1,8] >= -1
 udetSC+ekA: - udetXa[1,8] + udetYa[1,8] <= 0
 udetSC+ekA: - ugeqka[0,8] + udetYa[1,8] <= 0
 udetSC+ekA: - ugeqka[0,8] - udetXa[1,8] + udetYa[1,8] >= -1
 ufrSC: uDY[1,8] - udy[1,8] - uFrSC[1,8] >= 0
 ufrSC: uDetY[1,8] - uFrSC[1,8] >= 0
 ufrSC: - uDY[1,8] + udy[1,8] - uDetY[1,8] + uFrSC[1,8] >= -1
 ufrSCa: uDY[1,8] - udy[1,8] - ufrSCa[1,8] >= 0
 ufrSCa: udetYa[1,8] - ufrSCa[1,8] >= 0
 ufrSCa: - uDY[1,8] + udy[1,8] - udetYa[1,8] + ufrSCa[1,8] >= -1
 udetSC+ek: - uDetX[1,9] + uDetY[1,9] <= 0
 udetSC+ek: - ugEQK[0,9] + uDetY[1,9] <= 0
 udetSC+ek: - ugEQK[0,9] - uDetX[1,9] + uDetY[1,9] >= -1
 udetSC+ekA: - udetXa[1,9] + udetYa[1,9] <= 0
 udetSC+ekA: - ugeqka[0,9] + udetYa[1,9] <= 0
 udetSC+ekA: - ugeqka[0,9] - udetXa[1,9] + udetYa[1,9] >= -1
 ufrSC: uDY[1,9] - udy[1,9] - uFrSC[1,9] >= 0
 ufrSC: uDetY[1,9] - uFrSC[1,9] >= 0
 ufrSC: - uDY[1,9] + udy[1,9] - uDetY[1,9] + uFrSC[1,9] >= -1
 ufrSCa: uDY[1,9] - udy[1,9] - ufrSCa[1,9] >= 0
 ufrSCa: udetYa[1,9] - ufrSCa[1,9] >= 0
 ufrSCa: - uDY[1,9] + udy[1,9] - udetYa[1,9] + ufrSCa[1,9] >= -1
 udetSC+ek: - uDetX[1,10] + uDetY[1,10] <= 0
 udetSC+ek: - ugEQK[0,10] + uDetY[1,10] <= 0
 udetSC+ek: - ugEQK[0,10] - uDetX[1,10] + uDetY[1,10] >= -1
 udetSC+ekA: - udetXa[1,10] + udetYa[1,10] <= 0
 udetSC+ekA: - ugeqka[0,10] + udetYa[1,10] <= 0
 udetSC+ekA: - ugeqka[0,10] - udetXa[1,10] + udetYa[1,10] >= -1
 ufrSC: uDY[1,10] - udy[1,10] - uFrSC[1,10] >= 0
 ufrSC: uDetY[1,10] - uFrSC[1,10] >= 0
 ufrSC: - uDY[1,10] + udy[1,10] - uDetY[1,10] + uFrSC[1,10] >= -1
 ufrSCa: uDY[1,10] - udy[1,10] - ufrSCa[1,10] >= 0
 ufrSCa: udetYa[1,10] - ufrSCa[1,10] >= 0
 ufrSCa: - uDY[1,10] + udy[1,10] - udetYa[1,10] + ufrSCa[1,10] >= -1
 udetSC+ek: - uDetX[1,11] + uDetY[1,11] <= 0
 udetSC+ek: - ugEQK[0,11] + uDetY[1,11] <= 0
 udetSC+ek: - ugEQK[0,11] - uDetX[1,11] + uDetY[1,11] >= -1
 udetSC+ekA: - udetXa[1,11] + udetYa[1,11] <= 0
 udetSC+ekA: - ugeqka[0,11] + udetYa[1,11] <= 0
 udetSC+ekA: - ugeqka[0,11] - udetXa[1,11] + udetYa[1,11] >= -1
 ufrSC: uDY[1,11] - udy[1,11] - uFrSC[1,11] >= 0
 ufrSC: uDetY[1,11] - uFrSC[1,11] >= 0
 ufrSC: - uDY[1,11] + udy[1,11] - uDetY[1,11] + uFrSC[1,11] >= -1
 ufrSCa: uDY[1,11] - udy[1,11] - ufrSCa[1,11] >= 0
 ufrSCa: udetYa[1,11] - ufrSCa[1,11] >= 0
 ufrSCa: - uDY[1,11] + udy[1,11] - udetYa[1,11] + ufrSCa[1,11] >= -1
 udetSC+ek: - uDetX[1,12] + uDetY[1,12] <= 0
 udetSC+ek: - ugEQK[0,12] + uDetY[1,12] <= 0
 udetSC+ek: - ugEQK[0,12] - uDetX[1,12] + uDetY[1,12] >= -1
 udetSC+ekA: - udetXa[1,12] + udetYa[1,12] <= 0
 udetSC+ekA: - ugeqka[0,12] + udetYa[1,12] <= 0
 udetSC+ekA: - ugeqka[0,12] - udetXa[1,12] + udetYa[1,12] >= -1
 ufrSC: uDY[1,12] - udy[1,12] - uFrSC[1,12] >= 0
 ufrSC: uDetY[1,12] - uFrSC[1,12] >= 0
 ufrSC: - uDY[1,12] + udy[1,12] - uDetY[1,12] + uFrSC[1,12] >= -1
 ufrSCa: uDY[1,12] - udy[1,12] - ufrSCa[1,12] >= 0
 ufrSCa: udetYa[1,12] - ufrSCa[1,12] >= 0
 ufrSCa: - uDY[1,12] + udy[1,12] - udetYa[1,12] + ufrSCa[1,12] >= -1
 udetSC+ek: - uDetX[1,13] + uDetY[1,13] <= 0
 udetSC+ek: - ugEQK[0,13] + uDetY[1,13] <= 0
 udetSC+ek: - ugEQK[0,13] - uDetX[1,13] + uDetY[1,13] >= -1
 udetSC+ekA: - udetXa[1,13] + udetYa[1,13] <= 0
 udetSC+ekA: - ugeqka[0,13] + udetYa[1,13] <= 0
 udetSC+ekA: - ugeqka[0,13] - udetXa[1,13] + udetYa[1,13] >= -1
 ufrSC: uDY[1,13] - udy[1,13] - uFrSC[1,13] >= 0
 ufrSC: uDetY[1,13] - uFrSC[1,13] >= 0
 ufrSC: - uDY[1,13] + udy[1,13] - uDetY[1,13] + uFrSC[1,13] >= -1
 ufrSCa: uDY[1,13] - udy[1,13] - ufrSCa[1,13] >= 0
 ufrSCa: udetYa[1,13] - ufrSCa[1,13] >= 0
 ufrSCa: - uDY[1,13] + udy[1,13] - udetYa[1,13] + ufrSCa[1,13] >= -1
 udetSC+ek: - uDetX[1,14] + uDetY[1,14] <= 0
 udetSC+ek: - ugEQK[0,14] + uDetY[1,14] <= 0
 udetSC+ek: - ugEQK[0,14] - uDetX[1,14] + uDetY[1,14] >= -1
 udetSC+ekA: - udetXa[1,14] + udetYa[1,14] <= 0
 udetSC+ekA: - ugeqka[0,14] + udetYa[1,14] <= 0
 udetSC+ekA: - ugeqka[0,14] - udetXa[1,14] + udetYa[1,14] >= -1
 ufrSC: uDY[1,14] - udy[1,14] - uFrSC[1,14] >= 0
 ufrSC: uDetY[1,14] - uFrSC[1,14] >= 0
 ufrSC: - uDY[1,14] + udy[1,14] - uDetY[1,14] + uFrSC[1,14] >= -1
 ufrSCa: uDY[1,14] - udy[1,14] - ufrSCa[1,14] >= 0
 ufrSCa: udetYa[1,14] - ufrSCa[1,14] >= 0
 ufrSCa: - uDY[1,14] + udy[1,14] - udetYa[1,14] + ufrSCa[1,14] >= -1
 udetSC+ek: - uDetX[1,15] + uDetY[1,15] <= 0
 udetSC+ek: - ugEQK[0,15] + uDetY[1,15] <= 0
 udetSC+ek: - ugEQK[0,15] - uDetX[1,15] + uDetY[1,15] >= -1
 udetSC+ekA: - udetXa[1,15] + udetYa[1,15] <= 0
 udetSC+ekA: - ugeqka[0,15] + udetYa[1,15] <= 0
 udetSC+ekA: - ugeqka[0,15] - udetXa[1,15] + udetYa[1,15] >= -1
 ufrSC: uDY[1,15] - udy[1,15] - uFrSC[1,15] >= 0
 ufrSC: uDetY[1,15] - uFrSC[1,15] >= 0
 ufrSC: - uDY[1,15] + udy[1,15] - uDetY[1,15] + uFrSC[1,15] >= -1
 ufrSCa: uDY[1,15] - udy[1,15] - ufrSCa[1,15] >= 0
 ufrSCa: udetYa[1,15] - ufrSCa[1,15] >= 0
 ufrSCa: - uDY[1,15] + udy[1,15] - udetYa[1,15] + ufrSCa[1,15] >= -1
 udetSC+ek: - uDetX[2,0] + uDetY[2,0] <= 0
 udetSC+ek: - ugEQK[1,0] + uDetY[2,0] <= 0
 udetSC+ek: - ugEQK[1,0] - uDetX[2,0] + uDetY[2,0] >= -1
 udetSC+ekA: - udetXa[2,0] + udetYa[2,0] <= 0
 udetSC+ekA: - ugeqka[1,0] + udetYa[2,0] <= 0
 udetSC+ekA: - ugeqka[1,0] - udetXa[2,0] + udetYa[2,0] >= -1
 ufrSC: uDY[2,0] - udy[2,0] - uFrSC[2,0] >= 0
 ufrSC: uDetY[2,0] - uFrSC[2,0] >= 0
 ufrSC: - uDY[2,0] + udy[2,0] - uDetY[2,0] + uFrSC[2,0] >= -1
 ufrSCa: uDY[2,0] - udy[2,0] - ufrSCa[2,0] >= 0
 ufrSCa: udetYa[2,0] - ufrSCa[2,0] >= 0
 ufrSCa: - uDY[2,0] + udy[2,0] - udetYa[2,0] + ufrSCa[2,0] >= -1
 udetSC+ek: - uDetX[2,1] + uDetY[2,1] <= 0
 udetSC+ek: - ugEQK[1,1] + uDetY[2,1] <= 0
 udetSC+ek: - ugEQK[1,1] - uDetX[2,1] + uDetY[2,1] >= -1
 udetSC+ekA: - udetXa[2,1] + udetYa[2,1] <= 0
 udetSC+ekA: - ugeqka[1,1] + udetYa[2,1] <= 0
 udetSC+ekA: - ugeqka[1,1] - udetXa[2,1] + udetYa[2,1] >= -1
 ufrSC: uDY[2,1] - udy[2,1] - uFrSC[2,1] >= 0
 ufrSC: uDetY[2,1] - uFrSC[2,1] >= 0
 ufrSC: - uDY[2,1] + udy[2,1] - uDetY[2,1] + uFrSC[2,1] >= -1
 ufrSCa: uDY[2,1] - udy[2,1] - ufrSCa[2,1] >= 0
 ufrSCa: udetYa[2,1] - ufrSCa[2,1] >= 0
 ufrSCa: - uDY[2,1] + udy[2,1] - udetYa[2,1] + ufrSCa[2,1] >= -1
 udetSC+ek: - uDetX[2,2] + uDetY[2,2] <= 0
 udetSC+ek: - ugEQK[1,2] + uDetY[2,2] <= 0
 udetSC+ek: - ugEQK[1,2] - uDetX[2,2] + uDetY[2,2] >= -1
 udetSC+ekA: - udetXa[2,2] + udetYa[2,2] <= 0
 udetSC+ekA: - ugeqka[1,2] + udetYa[2,2] <= 0
 udetSC+ekA: - ugeqka[1,2] - udetXa[2,2] + udetYa[2,2] >= -1
 ufrSC: uDY[2,2] - udy[2,2] - uFrSC[2,2] >= 0
 ufrSC: uDetY[2,2] - uFrSC[2,2] >= 0
 ufrSC: - uDY[2,2] + udy[2,2] - uDetY[2,2] + uFrSC[2,2] >= -1
 ufrSCa: uDY[2,2] - udy[2,2] - ufrSCa[2,2] >= 0
 ufrSCa: udetYa[2,2] - ufrSCa[2,2] >= 0
 ufrSCa: - uDY[2,2] + udy[2,2] - udetYa[2,2] + ufrSCa[2,2] >= -1
 udetSC+ek: - uDetX[2,3] + uDetY[2,3] <= 0
 udetSC+ek: - ugEQK[1,3] + uDetY[2,3] <= 0
 udetSC+ek: - ugEQK[1,3] - uDetX[2,3] + uDetY[2,3] >= -1
 udetSC+ekA: - udetXa[2,3] + udetYa[2,3] <= 0
 udetSC+ekA: - ugeqka[1,3] + udetYa[2,3] <= 0
 udetSC+ekA: - ugeqka[1,3] - udetXa[2,3] + udetYa[2,3] >= -1
 ufrSC: uDY[2,3] - udy[2,3] - uFrSC[2,3] >= 0
 ufrSC: uDetY[2,3] - uFrSC[2,3] >= 0
 ufrSC: - uDY[2,3] + udy[2,3] - uDetY[2,3] + uFrSC[2,3] >= -1
 ufrSCa: uDY[2,3] - udy[2,3] - ufrSCa[2,3] >= 0
 ufrSCa: udetYa[2,3] - ufrSCa[2,3] >= 0
 ufrSCa: - uDY[2,3] + udy[2,3] - udetYa[2,3] + ufrSCa[2,3] >= -1
 udetSC+ek: - uDetX[2,4] + uDetY[2,4] <= 0
 udetSC+ek: - ugEQK[1,4] + uDetY[2,4] <= 0
 udetSC+ek: - ugEQK[1,4] - uDetX[2,4] + uDetY[2,4] >= -1
 udetSC+ekA: - udetXa[2,4] + udetYa[2,4] <= 0
 udetSC+ekA: - ugeqka[1,4] + udetYa[2,4] <= 0
 udetSC+ekA: - ugeqka[1,4] - udetXa[2,4] + udetYa[2,4] >= -1
 ufrSC: uDY[2,4] - udy[2,4] - uFrSC[2,4] >= 0
 ufrSC: uDetY[2,4] - uFrSC[2,4] >= 0
 ufrSC: - uDY[2,4] + udy[2,4] - uDetY[2,4] + uFrSC[2,4] >= -1
 ufrSCa: uDY[2,4] - udy[2,4] - ufrSCa[2,4] >= 0
 ufrSCa: udetYa[2,4] - ufrSCa[2,4] >= 0
 ufrSCa: - uDY[2,4] + udy[2,4] - udetYa[2,4] + ufrSCa[2,4] >= -1
 udetSC+ek: - uDetX[2,5] + uDetY[2,5] <= 0
 udetSC+ek: - ugEQK[1,5] + uDetY[2,5] <= 0
 udetSC+ek: - ugEQK[1,5] - uDetX[2,5] + uDetY[2,5] >= -1
 udetSC+ekA: - udetXa[2,5] + udetYa[2,5] <= 0
 udetSC+ekA: - ugeqka[1,5] + udetYa[2,5] <= 0
 udetSC+ekA: - ugeqka[1,5] - udetXa[2,5] + udetYa[2,5] >= -1
 ufrSC: uDY[2,5] - udy[2,5] - uFrSC[2,5] >= 0
 ufrSC: uDetY[2,5] - uFrSC[2,5] >= 0
 ufrSC: - uDY[2,5] + udy[2,5] - uDetY[2,5] + uFrSC[2,5] >= -1
 ufrSCa: uDY[2,5] - udy[2,5] - ufrSCa[2,5] >= 0
 ufrSCa: udetYa[2,5] - ufrSCa[2,5] >= 0
 ufrSCa: - uDY[2,5] + udy[2,5] - udetYa[2,5] + ufrSCa[2,5] >= -1
 udetSC+ek: - uDetX[2,6] + uDetY[2,6] <= 0
 udetSC+ek: - ugEQK[1,6] + uDetY[2,6] <= 0
 udetSC+ek: - ugEQK[1,6] - uDetX[2,6] + uDetY[2,6] >= -1
 udetSC+ekA: - udetXa[2,6] + udetYa[2,6] <= 0
 udetSC+ekA: - ugeqka[1,6] + udetYa[2,6] <= 0
 udetSC+ekA: - ugeqka[1,6] - udetXa[2,6] + udetYa[2,6] >= -1
 ufrSC: uDY[2,6] - udy[2,6] - uFrSC[2,6] >= 0
 ufrSC: uDetY[2,6] - uFrSC[2,6] >= 0
 ufrSC: - uDY[2,6] + udy[2,6] - uDetY[2,6] + uFrSC[2,6] >= -1
 ufrSCa: uDY[2,6] - udy[2,6] - ufrSCa[2,6] >= 0
 ufrSCa: udetYa[2,6] - ufrSCa[2,6] >= 0
 ufrSCa: - uDY[2,6] + udy[2,6] - udetYa[2,6] + ufrSCa[2,6] >= -1
 udetSC+ek: - uDetX[2,7] + uDetY[2,7] <= 0
 udetSC+ek: - ugEQK[1,7] + uDetY[2,7] <= 0
 udetSC+ek: - ugEQK[1,7] - uDetX[2,7] + uDetY[2,7] >= -1
 udetSC+ekA: - udetXa[2,7] + udetYa[2,7] <= 0
 udetSC+ekA: - ugeqka[1,7] + udetYa[2,7] <= 0
 udetSC+ekA: - ugeqka[1,7] - udetXa[2,7] + udetYa[2,7] >= -1
 ufrSC: uDY[2,7] - udy[2,7] - uFrSC[2,7] >= 0
 ufrSC: uDetY[2,7] - uFrSC[2,7] >= 0
 ufrSC: - uDY[2,7] + udy[2,7] - uDetY[2,7] + uFrSC[2,7] >= -1
 ufrSCa: uDY[2,7] - udy[2,7] - ufrSCa[2,7] >= 0
 ufrSCa: udetYa[2,7] - ufrSCa[2,7] >= 0
 ufrSCa: - uDY[2,7] + udy[2,7] - udetYa[2,7] + ufrSCa[2,7] >= -1
 udetSC+ek: - uDetX[2,8] + uDetY[2,8] <= 0
 udetSC+ek: - ugEQK[1,8] + uDetY[2,8] <= 0
 udetSC+ek: - ugEQK[1,8] - uDetX[2,8] + uDetY[2,8] >= -1
 udetSC+ekA: - udetXa[2,8] + udetYa[2,8] <= 0
 udetSC+ekA: - ugeqka[1,8] + udetYa[2,8] <= 0
 udetSC+ekA: - ugeqka[1,8] - udetXa[2,8] + udetYa[2,8] >= -1
 ufrSC: uDY[2,8] - udy[2,8] - uFrSC[2,8] >= 0
 ufrSC: uDetY[2,8] - uFrSC[2,8] >= 0
 ufrSC: - uDY[2,8] + udy[2,8] - uDetY[2,8] + uFrSC[2,8] >= -1
 ufrSCa: uDY[2,8] - udy[2,8] - ufrSCa[2,8] >= 0
 ufrSCa: udetYa[2,8] - ufrSCa[2,8] >= 0
 ufrSCa: - uDY[2,8] + udy[2,8] - udetYa[2,8] + ufrSCa[2,8] >= -1
 udetSC+ek: - uDetX[2,9] + uDetY[2,9] <= 0
 udetSC+ek: - ugEQK[1,9] + uDetY[2,9] <= 0
 udetSC+ek: - ugEQK[1,9] - uDetX[2,9] + uDetY[2,9] >= -1
 udetSC+ekA: - udetXa[2,9] + udetYa[2,9] <= 0
 udetSC+ekA: - ugeqka[1,9] + udetYa[2,9] <= 0
 udetSC+ekA: - ugeqka[1,9] - udetXa[2,9] + udetYa[2,9] >= -1
 ufrSC: uDY[2,9] - udy[2,9] - uFrSC[2,9] >= 0
 ufrSC: uDetY[2,9] - uFrSC[2,9] >= 0
 ufrSC: - uDY[2,9] + udy[2,9] - uDetY[2,9] + uFrSC[2,9] >= -1
 ufrSCa: uDY[2,9] - udy[2,9] - ufrSCa[2,9] >= 0
 ufrSCa: udetYa[2,9] - ufrSCa[2,9] >= 0
 ufrSCa: - uDY[2,9] + udy[2,9] - udetYa[2,9] + ufrSCa[2,9] >= -1
 udetSC+ek: - uDetX[2,10] + uDetY[2,10] <= 0
 udetSC+ek: - ugEQK[1,10] + uDetY[2,10] <= 0
 udetSC+ek: - ugEQK[1,10] - uDetX[2,10] + uDetY[2,10] >= -1
 udetSC+ekA: - udetXa[2,10] + udetYa[2,10] <= 0
 udetSC+ekA: - ugeqka[1,10] + udetYa[2,10] <= 0
 udetSC+ekA: - ugeqka[1,10] - udetXa[2,10] + udetYa[2,10] >= -1
 ufrSC: uDY[2,10] - udy[2,10] - uFrSC[2,10] >= 0
 ufrSC: uDetY[2,10] - uFrSC[2,10] >= 0
 ufrSC: - uDY[2,10] + udy[2,10] - uDetY[2,10] + uFrSC[2,10] >= -1
 ufrSCa: uDY[2,10] - udy[2,10] - ufrSCa[2,10] >= 0
 ufrSCa: udetYa[2,10] - ufrSCa[2,10] >= 0
 ufrSCa: - uDY[2,10] + udy[2,10] - udetYa[2,10] + ufrSCa[2,10] >= -1
 udetSC+ek: - uDetX[2,11] + uDetY[2,11] <= 0
 udetSC+ek: - ugEQK[1,11] + uDetY[2,11] <= 0
 udetSC+ek: - ugEQK[1,11] - uDetX[2,11] + uDetY[2,11] >= -1
 udetSC+ekA: - udetXa[2,11] + udetYa[2,11] <= 0
 udetSC+ekA: - ugeqka[1,11] + udetYa[2,11] <= 0
 udetSC+ekA: - ugeqka[1,11] - udetXa[2,11] + udetYa[2,11] >= -1
 ufrSC: uDY[2,11] - udy[2,11] - uFrSC[2,11] >= 0
 ufrSC: uDetY[2,11] - uFrSC[2,11] >= 0
 ufrSC: - uDY[2,11] + udy[2,11] - uDetY[2,11] + uFrSC[2,11] >= -1
 ufrSCa: uDY[2,11] - udy[2,11] - ufrSCa[2,11] >= 0
 ufrSCa: udetYa[2,11] - ufrSCa[2,11] >= 0
 ufrSCa: - uDY[2,11] + udy[2,11] - udetYa[2,11] + ufrSCa[2,11] >= -1
 udetSC+ek: - uDetX[2,12] + uDetY[2,12] <= 0
 udetSC+ek: - ugEQK[1,12] + uDetY[2,12] <= 0
 udetSC+ek: - ugEQK[1,12] - uDetX[2,12] + uDetY[2,12] >= -1
 udetSC+ekA: - udetXa[2,12] + udetYa[2,12] <= 0
 udetSC+ekA: - ugeqka[1,12] + udetYa[2,12] <= 0
 udetSC+ekA: - ugeqka[1,12] - udetXa[2,12] + udetYa[2,12] >= -1
 ufrSC: uDY[2,12] - udy[2,12] - uFrSC[2,12] >= 0
 ufrSC: uDetY[2,12] - uFrSC[2,12] >= 0
 ufrSC: - uDY[2,12] + udy[2,12] - uDetY[2,12] + uFrSC[2,12] >= -1
 ufrSCa: uDY[2,12] - udy[2,12] - ufrSCa[2,12] >= 0
 ufrSCa: udetYa[2,12] - ufrSCa[2,12] >= 0
 ufrSCa: - uDY[2,12] + udy[2,12] - udetYa[2,12] + ufrSCa[2,12] >= -1
 udetSC+ek: - uDetX[2,13] + uDetY[2,13] <= 0
 udetSC+ek: - ugEQK[1,13] + uDetY[2,13] <= 0
 udetSC+ek: - ugEQK[1,13] - uDetX[2,13] + uDetY[2,13] >= -1
 udetSC+ekA: - udetXa[2,13] + udetYa[2,13] <= 0
 udetSC+ekA: - ugeqka[1,13] + udetYa[2,13] <= 0
 udetSC+ekA: - ugeqka[1,13] - udetXa[2,13] + udetYa[2,13] >= -1
 ufrSC: uDY[2,13] - udy[2,13] - uFrSC[2,13] >= 0
 ufrSC: uDetY[2,13] - uFrSC[2,13] >= 0
 ufrSC: - uDY[2,13] + udy[2,13] - uDetY[2,13] + uFrSC[2,13] >= -1
 ufrSCa: uDY[2,13] - udy[2,13] - ufrSCa[2,13] >= 0
 ufrSCa: udetYa[2,13] - ufrSCa[2,13] >= 0
 ufrSCa: - uDY[2,13] + udy[2,13] - udetYa[2,13] + ufrSCa[2,13] >= -1
 udetSC+ek: - uDetX[2,14] + uDetY[2,14] <= 0
 udetSC+ek: - ugEQK[1,14] + uDetY[2,14] <= 0
 udetSC+ek: - ugEQK[1,14] - uDetX[2,14] + uDetY[2,14] >= -1
 udetSC+ekA: - udetXa[2,14] + udetYa[2,14] <= 0
 udetSC+ekA: - ugeqka[1,14] + udetYa[2,14] <= 0
 udetSC+ekA: - ugeqka[1,14] - udetXa[2,14] + udetYa[2,14] >= -1
 ufrSC: uDY[2,14] - udy[2,14] - uFrSC[2,14] >= 0
 ufrSC: uDetY[2,14] - uFrSC[2,14] >= 0
 ufrSC: - uDY[2,14] + udy[2,14] - uDetY[2,14] + uFrSC[2,14] >= -1
 ufrSCa: uDY[2,14] - udy[2,14] - ufrSCa[2,14] >= 0
 ufrSCa: udetYa[2,14] - ufrSCa[2,14] >= 0
 ufrSCa: - uDY[2,14] + udy[2,14] - udetYa[2,14] + ufrSCa[2,14] >= -1
 udetSC+ek: - uDetX[2,15] + uDetY[2,15] <= 0
 udetSC+ek: - ugEQK[1,15] + uDetY[2,15] <= 0
 udetSC+ek: - ugEQK[1,15] - uDetX[2,15] + uDetY[2,15] >= -1
 udetSC+ekA: - udetXa[2,15] + udetYa[2,15] <= 0
 udetSC+ekA: - ugeqka[1,15] + udetYa[2,15] <= 0
 udetSC+ekA: - ugeqka[1,15] - udetXa[2,15] + udetYa[2,15] >= -1
 ufrSC: uDY[2,15] - udy[2,15] - uFrSC[2,15] >= 0
 ufrSC: uDetY[2,15] - uFrSC[2,15] >= 0
 ufrSC: - uDY[2,15] + udy[2,15] - uDetY[2,15] + uFrSC[2,15] >= -1
 ufrSCa: uDY[2,15] - udy[2,15] - ufrSCa[2,15] >= 0
 ufrSCa: udetYa[2,15] - ufrSCa[2,15] >= 0
 ufrSCa: - uDY[2,15] + udy[2,15] - udetYa[2,15] + ufrSCa[2,15] >= -1
 udetSC+ek: - uDetX[3,0] + uDetY[3,0] <= 0
 udetSC+ek: - ugEQK[2,0] + uDetY[3,0] <= 0
 udetSC+ek: - ugEQK[2,0] - uDetX[3,0] + uDetY[3,0] >= -1
 udetSC+ekA: - udetXa[3,0] + udetYa[3,0] <= 0
 udetSC+ekA: - ugeqka[2,0] + udetYa[3,0] <= 0
 udetSC+ekA: - ugeqka[2,0] - udetXa[3,0] + udetYa[3,0] >= -1
 ufrSC: uDY[3,0] - udy[3,0] - uFrSC[3,0] >= 0
 ufrSC: uDetY[3,0] - uFrSC[3,0] >= 0
 ufrSC: - uDY[3,0] + udy[3,0] - uDetY[3,0] + uFrSC[3,0] >= -1
 ufrSCa: uDY[3,0] - udy[3,0] - ufrSCa[3,0] >= 0
 ufrSCa: udetYa[3,0] - ufrSCa[3,0] >= 0
 ufrSCa: - uDY[3,0] + udy[3,0] - udetYa[3,0] + ufrSCa[3,0] >= -1
 udetSC+ek: - uDetX[3,1] + uDetY[3,1] <= 0
 udetSC+ek: - ugEQK[2,1] + uDetY[3,1] <= 0
 udetSC+ek: - ugEQK[2,1] - uDetX[3,1] + uDetY[3,1] >= -1
 udetSC+ekA: - udetXa[3,1] + udetYa[3,1] <= 0
 udetSC+ekA: - ugeqka[2,1] + udetYa[3,1] <= 0
 udetSC+ekA: - ugeqka[2,1] - udetXa[3,1] + udetYa[3,1] >= -1
 ufrSC: uDY[3,1] - udy[3,1] - uFrSC[3,1] >= 0
 ufrSC: uDetY[3,1] - uFrSC[3,1] >= 0
 ufrSC: - uDY[3,1] + udy[3,1] - uDetY[3,1] + uFrSC[3,1] >= -1
 ufrSCa: uDY[3,1] - udy[3,1] - ufrSCa[3,1] >= 0
 ufrSCa: udetYa[3,1] - ufrSCa[3,1] >= 0
 ufrSCa: - uDY[3,1] + udy[3,1] - udetYa[3,1] + ufrSCa[3,1] >= -1
 udetSC+ek: - uDetX[3,2] + uDetY[3,2] <= 0
 udetSC+ek: - ugEQK[2,2] + uDetY[3,2] <= 0
 udetSC+ek: - ugEQK[2,2] - uDetX[3,2] + uDetY[3,2] >= -1
 udetSC+ekA: - udetXa[3,2] + udetYa[3,2] <= 0
 udetSC+ekA: - ugeqka[2,2] + udetYa[3,2] <= 0
 udetSC+ekA: - ugeqka[2,2] - udetXa[3,2] + udetYa[3,2] >= -1
 ufrSC: uDY[3,2] - udy[3,2] - uFrSC[3,2] >= 0
 ufrSC: uDetY[3,2] - uFrSC[3,2] >= 0
 ufrSC: - uDY[3,2] + udy[3,2] - uDetY[3,2] + uFrSC[3,2] >= -1
 ufrSCa: uDY[3,2] - udy[3,2] - ufrSCa[3,2] >= 0
 ufrSCa: udetYa[3,2] - ufrSCa[3,2] >= 0
 ufrSCa: - uDY[3,2] + udy[3,2] - udetYa[3,2] + ufrSCa[3,2] >= -1
 udetSC+ek: - uDetX[3,3] + uDetY[3,3] <= 0
 udetSC+ek: - ugEQK[2,3] + uDetY[3,3] <= 0
 udetSC+ek: - ugEQK[2,3] - uDetX[3,3] + uDetY[3,3] >= -1
 udetSC+ekA: - udetXa[3,3] + udetYa[3,3] <= 0
 udetSC+ekA: - ugeqka[2,3] + udetYa[3,3] <= 0
 udetSC+ekA: - ugeqka[2,3] - udetXa[3,3] + udetYa[3,3] >= -1
 ufrSC: uDY[3,3] - udy[3,3] - uFrSC[3,3] >= 0
 ufrSC: uDetY[3,3] - uFrSC[3,3] >= 0
 ufrSC: - uDY[3,3] + udy[3,3] - uDetY[3,3] + uFrSC[3,3] >= -1
 ufrSCa: uDY[3,3] - udy[3,3] - ufrSCa[3,3] >= 0
 ufrSCa: udetYa[3,3] - ufrSCa[3,3] >= 0
 ufrSCa: - uDY[3,3] + udy[3,3] - udetYa[3,3] + ufrSCa[3,3] >= -1
 udetSC+ek: - uDetX[3,4] + uDetY[3,4] <= 0
 udetSC+ek: - ugEQK[2,4] + uDetY[3,4] <= 0
 udetSC+ek: - ugEQK[2,4] - uDetX[3,4] + uDetY[3,4] >= -1
 udetSC+ekA: - udetXa[3,4] + udetYa[3,4] <= 0
 udetSC+ekA: - ugeqka[2,4] + udetYa[3,4] <= 0
 udetSC+ekA: - ugeqka[2,4] - udetXa[3,4] + udetYa[3,4] >= -1
 ufrSC: uDY[3,4] - udy[3,4] - uFrSC[3,4] >= 0
 ufrSC: uDetY[3,4] - uFrSC[3,4] >= 0
 ufrSC: - uDY[3,4] + udy[3,4] - uDetY[3,4] + uFrSC[3,4] >= -1
 ufrSCa: uDY[3,4] - udy[3,4] - ufrSCa[3,4] >= 0
 ufrSCa: udetYa[3,4] - ufrSCa[3,4] >= 0
 ufrSCa: - uDY[3,4] + udy[3,4] - udetYa[3,4] + ufrSCa[3,4] >= -1
 udetSC+ek: - uDetX[3,5] + uDetY[3,5] <= 0
 udetSC+ek: - ugEQK[2,5] + uDetY[3,5] <= 0
 udetSC+ek: - ugEQK[2,5] - uDetX[3,5] + uDetY[3,5] >= -1
 udetSC+ekA: - udetXa[3,5] + udetYa[3,5] <= 0
 udetSC+ekA: - ugeqka[2,5] + udetYa[3,5] <= 0
 udetSC+ekA: - ugeqka[2,5] - udetXa[3,5] + udetYa[3,5] >= -1
 ufrSC: uDY[3,5] - udy[3,5] - uFrSC[3,5] >= 0
 ufrSC: uDetY[3,5] - uFrSC[3,5] >= 0
 ufrSC: - uDY[3,5] + udy[3,5] - uDetY[3,5] + uFrSC[3,5] >= -1
 ufrSCa: uDY[3,5] - udy[3,5] - ufrSCa[3,5] >= 0
 ufrSCa: udetYa[3,5] - ufrSCa[3,5] >= 0
 ufrSCa: - uDY[3,5] + udy[3,5] - udetYa[3,5] + ufrSCa[3,5] >= -1
 udetSC+ek: - uDetX[3,6] + uDetY[3,6] <= 0
 udetSC+ek: - ugEQK[2,6] + uDetY[3,6] <= 0
 udetSC+ek: - ugEQK[2,6] - uDetX[3,6] + uDetY[3,6] >= -1
 udetSC+ekA: - udetXa[3,6] + udetYa[3,6] <= 0
 udetSC+ekA: - ugeqka[2,6] + udetYa[3,6] <= 0
 udetSC+ekA: - ugeqka[2,6] - udetXa[3,6] + udetYa[3,6] >= -1
 ufrSC: uDY[3,6] - udy[3,6] - uFrSC[3,6] >= 0
 ufrSC: uDetY[3,6] - uFrSC[3,6] >= 0
 ufrSC: - uDY[3,6] + udy[3,6] - uDetY[3,6] + uFrSC[3,6] >= -1
 ufrSCa: uDY[3,6] - udy[3,6] - ufrSCa[3,6] >= 0
 ufrSCa: udetYa[3,6] - ufrSCa[3,6] >= 0
 ufrSCa: - uDY[3,6] + udy[3,6] - udetYa[3,6] + ufrSCa[3,6] >= -1
 udetSC+ek: - uDetX[3,7] + uDetY[3,7] <= 0
 udetSC+ek: - ugEQK[2,7] + uDetY[3,7] <= 0
 udetSC+ek: - ugEQK[2,7] - uDetX[3,7] + uDetY[3,7] >= -1
 udetSC+ekA: - udetXa[3,7] + udetYa[3,7] <= 0
 udetSC+ekA: - ugeqka[2,7] + udetYa[3,7] <= 0
 udetSC+ekA: - ugeqka[2,7] - udetXa[3,7] + udetYa[3,7] >= -1
 ufrSC: uDY[3,7] - udy[3,7] - uFrSC[3,7] >= 0
 ufrSC: uDetY[3,7] - uFrSC[3,7] >= 0
 ufrSC: - uDY[3,7] + udy[3,7] - uDetY[3,7] + uFrSC[3,7] >= -1
 ufrSCa: uDY[3,7] - udy[3,7] - ufrSCa[3,7] >= 0
 ufrSCa: udetYa[3,7] - ufrSCa[3,7] >= 0
 ufrSCa: - uDY[3,7] + udy[3,7] - udetYa[3,7] + ufrSCa[3,7] >= -1
 udetSC+ek: - uDetX[3,8] + uDetY[3,8] <= 0
 udetSC+ek: - ugEQK[2,8] + uDetY[3,8] <= 0
 udetSC+ek: - ugEQK[2,8] - uDetX[3,8] + uDetY[3,8] >= -1
 udetSC+ekA: - udetXa[3,8] + udetYa[3,8] <= 0
 udetSC+ekA: - ugeqka[2,8] + udetYa[3,8] <= 0
 udetSC+ekA: - ugeqka[2,8] - udetXa[3,8] + udetYa[3,8] >= -1
 ufrSC: uDY[3,8] - udy[3,8] - uFrSC[3,8] >= 0
 ufrSC: uDetY[3,8] - uFrSC[3,8] >= 0
 ufrSC: - uDY[3,8] + udy[3,8] - uDetY[3,8] + uFrSC[3,8] >= -1
 ufrSCa: uDY[3,8] - udy[3,8] - ufrSCa[3,8] >= 0
 ufrSCa: udetYa[3,8] - ufrSCa[3,8] >= 0
 ufrSCa: - uDY[3,8] + udy[3,8] - udetYa[3,8] + ufrSCa[3,8] >= -1
 udetSC+ek: - uDetX[3,9] + uDetY[3,9] <= 0
 udetSC+ek: - ugEQK[2,9] + uDetY[3,9] <= 0
 udetSC+ek: - ugEQK[2,9] - uDetX[3,9] + uDetY[3,9] >= -1
 udetSC+ekA: - udetXa[3,9] + udetYa[3,9] <= 0
 udetSC+ekA: - ugeqka[2,9] + udetYa[3,9] <= 0
 udetSC+ekA: - ugeqka[2,9] - udetXa[3,9] + udetYa[3,9] >= -1
 ufrSC: uDY[3,9] - udy[3,9] - uFrSC[3,9] >= 0
 ufrSC: uDetY[3,9] - uFrSC[3,9] >= 0
 ufrSC: - uDY[3,9] + udy[3,9] - uDetY[3,9] + uFrSC[3,9] >= -1
 ufrSCa: uDY[3,9] - udy[3,9] - ufrSCa[3,9] >= 0
 ufrSCa: udetYa[3,9] - ufrSCa[3,9] >= 0
 ufrSCa: - uDY[3,9] + udy[3,9] - udetYa[3,9] + ufrSCa[3,9] >= -1
 udetSC+ek: - uDetX[3,10] + uDetY[3,10] <= 0
 udetSC+ek: - ugEQK[2,10] + uDetY[3,10] <= 0
 udetSC+ek: - ugEQK[2,10] - uDetX[3,10] + uDetY[3,10] >= -1
 udetSC+ekA: - udetXa[3,10] + udetYa[3,10] <= 0
 udetSC+ekA: - ugeqka[2,10] + udetYa[3,10] <= 0
 udetSC+ekA: - ugeqka[2,10] - udetXa[3,10] + udetYa[3,10] >= -1
 ufrSC: uDY[3,10] - udy[3,10] - uFrSC[3,10] >= 0
 ufrSC: uDetY[3,10] - uFrSC[3,10] >= 0
 ufrSC: - uDY[3,10] + udy[3,10] - uDetY[3,10] + uFrSC[3,10] >= -1
 ufrSCa: uDY[3,10] - udy[3,10] - ufrSCa[3,10] >= 0
 ufrSCa: udetYa[3,10] - ufrSCa[3,10] >= 0
 ufrSCa: - uDY[3,10] + udy[3,10] - udetYa[3,10] + ufrSCa[3,10] >= -1
 udetSC+ek: - uDetX[3,11] + uDetY[3,11] <= 0
 udetSC+ek: - ugEQK[2,11] + uDetY[3,11] <= 0
 udetSC+ek: - ugEQK[2,11] - uDetX[3,11] + uDetY[3,11] >= -1
 udetSC+ekA: - udetXa[3,11] + udetYa[3,11] <= 0
 udetSC+ekA: - ugeqka[2,11] + udetYa[3,11] <= 0
 udetSC+ekA: - ugeqka[2,11] - udetXa[3,11] + udetYa[3,11] >= -1
 ufrSC: uDY[3,11] - udy[3,11] - uFrSC[3,11] >= 0
 ufrSC: uDetY[3,11] - uFrSC[3,11] >= 0
 ufrSC: - uDY[3,11] + udy[3,11] - uDetY[3,11] + uFrSC[3,11] >= -1
 ufrSCa: uDY[3,11] - udy[3,11] - ufrSCa[3,11] >= 0
 ufrSCa: udetYa[3,11] - ufrSCa[3,11] >= 0
 ufrSCa: - uDY[3,11] + udy[3,11] - udetYa[3,11] + ufrSCa[3,11] >= -1
 udetSC+ek: - uDetX[3,12] + uDetY[3,12] <= 0
 udetSC+ek: - ugEQK[2,12] + uDetY[3,12] <= 0
 udetSC+ek: - ugEQK[2,12] - uDetX[3,12] + uDetY[3,12] >= -1
 udetSC+ekA: - udetXa[3,12] + udetYa[3,12] <= 0
 udetSC+ekA: - ugeqka[2,12] + udetYa[3,12] <= 0
 udetSC+ekA: - ugeqka[2,12] - udetXa[3,12] + udetYa[3,12] >= -1
 ufrSC: uDY[3,12] - udy[3,12] - uFrSC[3,12] >= 0
 ufrSC: uDetY[3,12] - uFrSC[3,12] >= 0
 ufrSC: - uDY[3,12] + udy[3,12] - uDetY[3,12] + uFrSC[3,12] >= -1
 ufrSCa: uDY[3,12] - udy[3,12] - ufrSCa[3,12] >= 0
 ufrSCa: udetYa[3,12] - ufrSCa[3,12] >= 0
 ufrSCa: - uDY[3,12] + udy[3,12] - udetYa[3,12] + ufrSCa[3,12] >= -1
 udetSC+ek: - uDetX[3,13] + uDetY[3,13] <= 0
 udetSC+ek: - ugEQK[2,13] + uDetY[3,13] <= 0
 udetSC+ek: - ugEQK[2,13] - uDetX[3,13] + uDetY[3,13] >= -1
 udetSC+ekA: - udetXa[3,13] + udetYa[3,13] <= 0
 udetSC+ekA: - ugeqka[2,13] + udetYa[3,13] <= 0
 udetSC+ekA: - ugeqka[2,13] - udetXa[3,13] + udetYa[3,13] >= -1
 ufrSC: uDY[3,13] - udy[3,13] - uFrSC[3,13] >= 0
 ufrSC: uDetY[3,13] - uFrSC[3,13] >= 0
 ufrSC: - uDY[3,13] + udy[3,13] - uDetY[3,13] + uFrSC[3,13] >= -1
 ufrSCa: uDY[3,13] - udy[3,13] - ufrSCa[3,13] >= 0
 ufrSCa: udetYa[3,13] - ufrSCa[3,13] >= 0
 ufrSCa: - uDY[3,13] + udy[3,13] - udetYa[3,13] + ufrSCa[3,13] >= -1
 udetSC+ek: - uDetX[3,14] + uDetY[3,14] <= 0
 udetSC+ek: - ugEQK[2,14] + uDetY[3,14] <= 0
 udetSC+ek: - ugEQK[2,14] - uDetX[3,14] + uDetY[3,14] >= -1
 udetSC+ekA: - udetXa[3,14] + udetYa[3,14] <= 0
 udetSC+ekA: - ugeqka[2,14] + udetYa[3,14] <= 0
 udetSC+ekA: - ugeqka[2,14] - udetXa[3,14] + udetYa[3,14] >= -1
 ufrSC: uDY[3,14] - udy[3,14] - uFrSC[3,14] >= 0
 ufrSC: uDetY[3,14] - uFrSC[3,14] >= 0
 ufrSC: - uDY[3,14] + udy[3,14] - uDetY[3,14] + uFrSC[3,14] >= -1
 ufrSCa: uDY[3,14] - udy[3,14] - ufrSCa[3,14] >= 0
 ufrSCa: udetYa[3,14] - ufrSCa[3,14] >= 0
 ufrSCa: - uDY[3,14] + udy[3,14] - udetYa[3,14] + ufrSCa[3,14] >= -1
 udetSC+ek: - uDetX[3,15] + uDetY[3,15] <= 0
 udetSC+ek: - ugEQK[2,15] + uDetY[3,15] <= 0
 udetSC+ek: - ugEQK[2,15] - uDetX[3,15] + uDetY[3,15] >= -1
 udetSC+ekA: - udetXa[3,15] + udetYa[3,15] <= 0
 udetSC+ekA: - ugeqka[2,15] + udetYa[3,15] <= 0
 udetSC+ekA: - ugeqka[2,15] - udetXa[3,15] + udetYa[3,15] >= -1
 ufrSC: uDY[3,15] - udy[3,15] - uFrSC[3,15] >= 0
 ufrSC: uDetY[3,15] - uFrSC[3,15] >= 0
 ufrSC: - uDY[3,15] + udy[3,15] - uDetY[3,15] + uFrSC[3,15] >= -1
 ufrSCa: uDY[3,15] - udy[3,15] - ufrSCa[3,15] >= 0
 ufrSCa: udetYa[3,15] - ufrSCa[3,15] >= 0
 ufrSCa: - uDY[3,15] + udy[3,15] - udetYa[3,15] + ufrSCa[3,15] >= -1
 udetSC+ek: - uDetX[4,0] + uDetY[4,0] <= 0
 udetSC+ek: - ugEQK[3,0] + uDetY[4,0] <= 0
 udetSC+ek: - ugEQK[3,0] - uDetX[4,0] + uDetY[4,0] >= -1
 udetSC+ekA: - udetXa[4,0] + udetYa[4,0] <= 0
 udetSC+ekA: - ugeqka[3,0] + udetYa[4,0] <= 0
 udetSC+ekA: - ugeqka[3,0] - udetXa[4,0] + udetYa[4,0] >= -1
 ufrSC: uDY[4,0] - udy[4,0] - uFrSC[4,0] >= 0
 ufrSC: uDetY[4,0] - uFrSC[4,0] >= 0
 ufrSC: - uDY[4,0] + udy[4,0] - uDetY[4,0] + uFrSC[4,0] >= -1
 ufrSCa: uDY[4,0] - udy[4,0] - ufrSCa[4,0] >= 0
 ufrSCa: udetYa[4,0] - ufrSCa[4,0] >= 0
 ufrSCa: - uDY[4,0] + udy[4,0] - udetYa[4,0] + ufrSCa[4,0] >= -1
 udetSC+ek: - uDetX[4,1] + uDetY[4,1] <= 0
 udetSC+ek: - ugEQK[3,1] + uDetY[4,1] <= 0
 udetSC+ek: - ugEQK[3,1] - uDetX[4,1] + uDetY[4,1] >= -1
 udetSC+ekA: - udetXa[4,1] + udetYa[4,1] <= 0
 udetSC+ekA: - ugeqka[3,1] + udetYa[4,1] <= 0
 udetSC+ekA: - ugeqka[3,1] - udetXa[4,1] + udetYa[4,1] >= -1
 ufrSC: uDY[4,1] - udy[4,1] - uFrSC[4,1] >= 0
 ufrSC: uDetY[4,1] - uFrSC[4,1] >= 0
 ufrSC: - uDY[4,1] + udy[4,1] - uDetY[4,1] + uFrSC[4,1] >= -1
 ufrSCa: uDY[4,1] - udy[4,1] - ufrSCa[4,1] >= 0
 ufrSCa: udetYa[4,1] - ufrSCa[4,1] >= 0
 ufrSCa: - uDY[4,1] + udy[4,1] - udetYa[4,1] + ufrSCa[4,1] >= -1
 udetSC+ek: - uDetX[4,2] + uDetY[4,2] <= 0
 udetSC+ek: - ugEQK[3,2] + uDetY[4,2] <= 0
 udetSC+ek: - ugEQK[3,2] - uDetX[4,2] + uDetY[4,2] >= -1
 udetSC+ekA: - udetXa[4,2] + udetYa[4,2] <= 0
 udetSC+ekA: - ugeqka[3,2] + udetYa[4,2] <= 0
 udetSC+ekA: - ugeqka[3,2] - udetXa[4,2] + udetYa[4,2] >= -1
 ufrSC: uDY[4,2] - udy[4,2] - uFrSC[4,2] >= 0
 ufrSC: uDetY[4,2] - uFrSC[4,2] >= 0
 ufrSC: - uDY[4,2] + udy[4,2] - uDetY[4,2] + uFrSC[4,2] >= -1
 ufrSCa: uDY[4,2] - udy[4,2] - ufrSCa[4,2] >= 0
 ufrSCa: udetYa[4,2] - ufrSCa[4,2] >= 0
 ufrSCa: - uDY[4,2] + udy[4,2] - udetYa[4,2] + ufrSCa[4,2] >= -1
 udetSC+ek: - uDetX[4,3] + uDetY[4,3] <= 0
 udetSC+ek: - ugEQK[3,3] + uDetY[4,3] <= 0
 udetSC+ek: - ugEQK[3,3] - uDetX[4,3] + uDetY[4,3] >= -1
 udetSC+ekA: - udetXa[4,3] + udetYa[4,3] <= 0
 udetSC+ekA: - ugeqka[3,3] + udetYa[4,3] <= 0
 udetSC+ekA: - ugeqka[3,3] - udetXa[4,3] + udetYa[4,3] >= -1
 ufrSC: uDY[4,3] - udy[4,3] - uFrSC[4,3] >= 0
 ufrSC: uDetY[4,3] - uFrSC[4,3] >= 0
 ufrSC: - uDY[4,3] + udy[4,3] - uDetY[4,3] + uFrSC[4,3] >= -1
 ufrSCa: uDY[4,3] - udy[4,3] - ufrSCa[4,3] >= 0
 ufrSCa: udetYa[4,3] - ufrSCa[4,3] >= 0
 ufrSCa: - uDY[4,3] + udy[4,3] - udetYa[4,3] + ufrSCa[4,3] >= -1
 udetSC+ek: - uDetX[4,4] + uDetY[4,4] <= 0
 udetSC+ek: - ugEQK[3,4] + uDetY[4,4] <= 0
 udetSC+ek: - ugEQK[3,4] - uDetX[4,4] + uDetY[4,4] >= -1
 udetSC+ekA: - udetXa[4,4] + udetYa[4,4] <= 0
 udetSC+ekA: - ugeqka[3,4] + udetYa[4,4] <= 0
 udetSC+ekA: - ugeqka[3,4] - udetXa[4,4] + udetYa[4,4] >= -1
 ufrSC: uDY[4,4] - udy[4,4] - uFrSC[4,4] >= 0
 ufrSC: uDetY[4,4] - uFrSC[4,4] >= 0
 ufrSC: - uDY[4,4] + udy[4,4] - uDetY[4,4] + uFrSC[4,4] >= -1
 ufrSCa: uDY[4,4] - udy[4,4] - ufrSCa[4,4] >= 0
 ufrSCa: udetYa[4,4] - ufrSCa[4,4] >= 0
 ufrSCa: - uDY[4,4] + udy[4,4] - udetYa[4,4] + ufrSCa[4,4] >= -1
 udetSC+ek: - uDetX[4,5] + uDetY[4,5] <= 0
 udetSC+ek: - ugEQK[3,5] + uDetY[4,5] <= 0
 udetSC+ek: - ugEQK[3,5] - uDetX[4,5] + uDetY[4,5] >= -1
 udetSC+ekA: - udetXa[4,5] + udetYa[4,5] <= 0
 udetSC+ekA: - ugeqka[3,5] + udetYa[4,5] <= 0
 udetSC+ekA: - ugeqka[3,5] - udetXa[4,5] + udetYa[4,5] >= -1
 ufrSC: uDY[4,5] - udy[4,5] - uFrSC[4,5] >= 0
 ufrSC: uDetY[4,5] - uFrSC[4,5] >= 0
 ufrSC: - uDY[4,5] + udy[4,5] - uDetY[4,5] + uFrSC[4,5] >= -1
 ufrSCa: uDY[4,5] - udy[4,5] - ufrSCa[4,5] >= 0
 ufrSCa: udetYa[4,5] - ufrSCa[4,5] >= 0
 ufrSCa: - uDY[4,5] + udy[4,5] - udetYa[4,5] + ufrSCa[4,5] >= -1
 udetSC+ek: - uDetX[4,6] + uDetY[4,6] <= 0
 udetSC+ek: - ugEQK[3,6] + uDetY[4,6] <= 0
 udetSC+ek: - ugEQK[3,6] - uDetX[4,6] + uDetY[4,6] >= -1
 udetSC+ekA: - udetXa[4,6] + udetYa[4,6] <= 0
 udetSC+ekA: - ugeqka[3,6] + udetYa[4,6] <= 0
 udetSC+ekA: - ugeqka[3,6] - udetXa[4,6] + udetYa[4,6] >= -1
 ufrSC: uDY[4,6] - udy[4,6] - uFrSC[4,6] >= 0
 ufrSC: uDetY[4,6] - uFrSC[4,6] >= 0
 ufrSC: - uDY[4,6] + udy[4,6] - uDetY[4,6] + uFrSC[4,6] >= -1
 ufrSCa: uDY[4,6] - udy[4,6] - ufrSCa[4,6] >= 0
 ufrSCa: udetYa[4,6] - ufrSCa[4,6] >= 0
 ufrSCa: - uDY[4,6] + udy[4,6] - udetYa[4,6] + ufrSCa[4,6] >= -1
 udetSC+ek: - uDetX[4,7] + uDetY[4,7] <= 0
 udetSC+ek: - ugEQK[3,7] + uDetY[4,7] <= 0
 udetSC+ek: - ugEQK[3,7] - uDetX[4,7] + uDetY[4,7] >= -1
 udetSC+ekA: - udetXa[4,7] + udetYa[4,7] <= 0
 udetSC+ekA: - ugeqka[3,7] + udetYa[4,7] <= 0
 udetSC+ekA: - ugeqka[3,7] - udetXa[4,7] + udetYa[4,7] >= -1
 ufrSC: uDY[4,7] - udy[4,7] - uFrSC[4,7] >= 0
 ufrSC: uDetY[4,7] - uFrSC[4,7] >= 0
 ufrSC: - uDY[4,7] + udy[4,7] - uDetY[4,7] + uFrSC[4,7] >= -1
 ufrSCa: uDY[4,7] - udy[4,7] - ufrSCa[4,7] >= 0
 ufrSCa: udetYa[4,7] - ufrSCa[4,7] >= 0
 ufrSCa: - uDY[4,7] + udy[4,7] - udetYa[4,7] + ufrSCa[4,7] >= -1
 udetSC+ek: - uDetX[4,8] + uDetY[4,8] <= 0
 udetSC+ek: - ugEQK[3,8] + uDetY[4,8] <= 0
 udetSC+ek: - ugEQK[3,8] - uDetX[4,8] + uDetY[4,8] >= -1
 udetSC+ekA: - udetXa[4,8] + udetYa[4,8] <= 0
 udetSC+ekA: - ugeqka[3,8] + udetYa[4,8] <= 0
 udetSC+ekA: - ugeqka[3,8] - udetXa[4,8] + udetYa[4,8] >= -1
 ufrSC: uDY[4,8] - udy[4,8] - uFrSC[4,8] >= 0
 ufrSC: uDetY[4,8] - uFrSC[4,8] >= 0
 ufrSC: - uDY[4,8] + udy[4,8] - uDetY[4,8] + uFrSC[4,8] >= -1
 ufrSCa: uDY[4,8] - udy[4,8] - ufrSCa[4,8] >= 0
 ufrSCa: udetYa[4,8] - ufrSCa[4,8] >= 0
 ufrSCa: - uDY[4,8] + udy[4,8] - udetYa[4,8] + ufrSCa[4,8] >= -1
 udetSC+ek: - uDetX[4,9] + uDetY[4,9] <= 0
 udetSC+ek: - ugEQK[3,9] + uDetY[4,9] <= 0
 udetSC+ek: - ugEQK[3,9] - uDetX[4,9] + uDetY[4,9] >= -1
 udetSC+ekA: - udetXa[4,9] + udetYa[4,9] <= 0
 udetSC+ekA: - ugeqka[3,9] + udetYa[4,9] <= 0
 udetSC+ekA: - ugeqka[3,9] - udetXa[4,9] + udetYa[4,9] >= -1
 ufrSC: uDY[4,9] - udy[4,9] - uFrSC[4,9] >= 0
 ufrSC: uDetY[4,9] - uFrSC[4,9] >= 0
 ufrSC: - uDY[4,9] + udy[4,9] - uDetY[4,9] + uFrSC[4,9] >= -1
 ufrSCa: uDY[4,9] - udy[4,9] - ufrSCa[4,9] >= 0
 ufrSCa: udetYa[4,9] - ufrSCa[4,9] >= 0
 ufrSCa: - uDY[4,9] + udy[4,9] - udetYa[4,9] + ufrSCa[4,9] >= -1
 udetSC+ek: - uDetX[4,10] + uDetY[4,10] <= 0
 udetSC+ek: - ugEQK[3,10] + uDetY[4,10] <= 0
 udetSC+ek: - ugEQK[3,10] - uDetX[4,10] + uDetY[4,10] >= -1
 udetSC+ekA: - udetXa[4,10] + udetYa[4,10] <= 0
 udetSC+ekA: - ugeqka[3,10] + udetYa[4,10] <= 0
 udetSC+ekA: - ugeqka[3,10] - udetXa[4,10] + udetYa[4,10] >= -1
 ufrSC: uDY[4,10] - udy[4,10] - uFrSC[4,10] >= 0
 ufrSC: uDetY[4,10] - uFrSC[4,10] >= 0
 ufrSC: - uDY[4,10] + udy[4,10] - uDetY[4,10] + uFrSC[4,10] >= -1
 ufrSCa: uDY[4,10] - udy[4,10] - ufrSCa[4,10] >= 0
 ufrSCa: udetYa[4,10] - ufrSCa[4,10] >= 0
 ufrSCa: - uDY[4,10] + udy[4,10] - udetYa[4,10] + ufrSCa[4,10] >= -1
 udetSC+ek: - uDetX[4,11] + uDetY[4,11] <= 0
 udetSC+ek: - ugEQK[3,11] + uDetY[4,11] <= 0
 udetSC+ek: - ugEQK[3,11] - uDetX[4,11] + uDetY[4,11] >= -1
 udetSC+ekA: - udetXa[4,11] + udetYa[4,11] <= 0
 udetSC+ekA: - ugeqka[3,11] + udetYa[4,11] <= 0
 udetSC+ekA: - ugeqka[3,11] - udetXa[4,11] + udetYa[4,11] >= -1
 ufrSC: uDY[4,11] - udy[4,11] - uFrSC[4,11] >= 0
 ufrSC: uDetY[4,11] - uFrSC[4,11] >= 0
 ufrSC: - uDY[4,11] + udy[4,11] - uDetY[4,11] + uFrSC[4,11] >= -1
 ufrSCa: uDY[4,11] - udy[4,11] - ufrSCa[4,11] >= 0
 ufrSCa: udetYa[4,11] - ufrSCa[4,11] >= 0
 ufrSCa: - uDY[4,11] + udy[4,11] - udetYa[4,11] + ufrSCa[4,11] >= -1
 udetSC+ek: - uDetX[4,12] + uDetY[4,12] <= 0
 udetSC+ek: - ugEQK[3,12] + uDetY[4,12] <= 0
 udetSC+ek: - ugEQK[3,12] - uDetX[4,12] + uDetY[4,12] >= -1
 udetSC+ekA: - udetXa[4,12] + udetYa[4,12] <= 0
 udetSC+ekA: - ugeqka[3,12] + udetYa[4,12] <= 0
 udetSC+ekA: - ugeqka[3,12] - udetXa[4,12] + udetYa[4,12] >= -1
 ufrSC: uDY[4,12] - udy[4,12] - uFrSC[4,12] >= 0
 ufrSC: uDetY[4,12] - uFrSC[4,12] >= 0
 ufrSC: - uDY[4,12] + udy[4,12] - uDetY[4,12] + uFrSC[4,12] >= -1
 ufrSCa: uDY[4,12] - udy[4,12] - ufrSCa[4,12] >= 0
 ufrSCa: udetYa[4,12] - ufrSCa[4,12] >= 0
 ufrSCa: - uDY[4,12] + udy[4,12] - udetYa[4,12] + ufrSCa[4,12] >= -1
 udetSC+ek: - uDetX[4,13] + uDetY[4,13] <= 0
 udetSC+ek: - ugEQK[3,13] + uDetY[4,13] <= 0
 udetSC+ek: - ugEQK[3,13] - uDetX[4,13] + uDetY[4,13] >= -1
 udetSC+ekA: - udetXa[4,13] + udetYa[4,13] <= 0
 udetSC+ekA: - ugeqka[3,13] + udetYa[4,13] <= 0
 udetSC+ekA: - ugeqka[3,13] - udetXa[4,13] + udetYa[4,13] >= -1
 ufrSC: uDY[4,13] - udy[4,13] - uFrSC[4,13] >= 0
 ufrSC: uDetY[4,13] - uFrSC[4,13] >= 0
 ufrSC: - uDY[4,13] + udy[4,13] - uDetY[4,13] + uFrSC[4,13] >= -1
 ufrSCa: uDY[4,13] - udy[4,13] - ufrSCa[4,13] >= 0
 ufrSCa: udetYa[4,13] - ufrSCa[4,13] >= 0
 ufrSCa: - uDY[4,13] + udy[4,13] - udetYa[4,13] + ufrSCa[4,13] >= -1
 udetSC+ek: - uDetX[4,14] + uDetY[4,14] <= 0
 udetSC+ek: - ugEQK[3,14] + uDetY[4,14] <= 0
 udetSC+ek: - ugEQK[3,14] - uDetX[4,14] + uDetY[4,14] >= -1
 udetSC+ekA: - udetXa[4,14] + udetYa[4,14] <= 0
 udetSC+ekA: - ugeqka[3,14] + udetYa[4,14] <= 0
 udetSC+ekA: - ugeqka[3,14] - udetXa[4,14] + udetYa[4,14] >= -1
 ufrSC: uDY[4,14] - udy[4,14] - uFrSC[4,14] >= 0
 ufrSC: uDetY[4,14] - uFrSC[4,14] >= 0
 ufrSC: - uDY[4,14] + udy[4,14] - uDetY[4,14] + uFrSC[4,14] >= -1
 ufrSCa: uDY[4,14] - udy[4,14] - ufrSCa[4,14] >= 0
 ufrSCa: udetYa[4,14] - ufrSCa[4,14] >= 0
 ufrSCa: - uDY[4,14] + udy[4,14] - udetYa[4,14] + ufrSCa[4,14] >= -1
 udetSC+ek: - uDetX[4,15] + uDetY[4,15] <= 0
 udetSC+ek: - ugEQK[3,15] + uDetY[4,15] <= 0
 udetSC+ek: - ugEQK[3,15] - uDetX[4,15] + uDetY[4,15] >= -1
 udetSC+ekA: - udetXa[4,15] + udetYa[4,15] <= 0
 udetSC+ekA: - ugeqka[3,15] + udetYa[4,15] <= 0
 udetSC+ekA: - ugeqka[3,15] - udetXa[4,15] + udetYa[4,15] >= -1
 ufrSC: uDY[4,15] - udy[4,15] - uFrSC[4,15] >= 0
 ufrSC: uDetY[4,15] - uFrSC[4,15] >= 0
 ufrSC: - uDY[4,15] + udy[4,15] - uDetY[4,15] + uFrSC[4,15] >= -1
 ufrSCa: uDY[4,15] - udy[4,15] - ufrSCa[4,15] >= 0
 ufrSCa: udetYa[4,15] - ufrSCa[4,15] >= 0
 ufrSCa: - uDY[4,15] + udy[4,15] - udetYa[4,15] + ufrSCa[4,15] >= -1
 udetY2Z: uDetY[1,0] - uDetZ[1,0] = 0
 udety2zA: udetYa[1,0] - udetZa[1,0] = 0
 udetY2Z: uDetY[1,1] - uDetZ[1,1] = 0
 udety2zA: udetYa[1,1] - udetZa[1,1] = 0
 udetY2Z: uDetY[1,2] - uDetZ[1,2] = 0
 udety2zA: udetYa[1,2] - udetZa[1,2] = 0
 udetY2Z: uDetY[1,3] - uDetZ[1,3] = 0
 udety2zA: udetYa[1,3] - udetZa[1,3] = 0
 udetY2Z: uDetY[1,4] - uDetZ[1,4] = 0
 udety2zA: udetYa[1,4] - udetZa[1,4] = 0
 udetY2Z: uDetY[1,5] - uDetZ[1,5] = 0
 udety2zA: udetYa[1,5] - udetZa[1,5] = 0
 udetY2Z: uDetY[1,6] - uDetZ[1,6] = 0
 udety2zA: udetYa[1,6] - udetZa[1,6] = 0
 udetY2Z: uDetY[1,7] - uDetZ[1,7] = 0
 udety2zA: udetYa[1,7] - udetZa[1,7] = 0
 udetY2Z: uDetY[1,8] - uDetZ[1,8] = 0
 udety2zA: udetYa[1,8] - udetZa[1,8] = 0
 udetY2Z: uDetY[1,9] - uDetZ[1,9] = 0
 udety2zA: udetYa[1,9] - udetZa[1,9] = 0
 udetY2Z: uDetY[1,10] - uDetZ[1,10] = 0
 udety2zA: udetYa[1,10] - udetZa[1,10] = 0
 udetY2Z: uDetY[1,11] - uDetZ[1,11] = 0
 udety2zA: udetYa[1,11] - udetZa[1,11] = 0
 udetY2Z: uDetY[1,12] - uDetZ[1,12] = 0
 udety2zA: udetYa[1,12] - udetZa[1,12] = 0
 udetY2Z: uDetY[1,13] - uDetZ[1,13] = 0
 udety2zA: udetYa[1,13] - udetZa[1,13] = 0
 udetY2Z: uDetY[1,14] - uDetZ[1,14] = 0
 udety2zA: udetYa[1,14] - udetZa[1,14] = 0
 udetY2Z: uDetY[1,15] - uDetZ[1,15] = 0
 udety2zA: udetYa[1,15] - udetZa[1,15] = 0
 udetY2Z: uDetY[2,0] - uDetZ[2,0] = 0
 udety2zA: udetYa[2,0] - udetZa[2,0] = 0
 udetY2Z: uDetY[2,1] - uDetZ[2,1] = 0
 udety2zA: udetYa[2,1] - udetZa[2,1] = 0
 udetY2Z: uDetY[2,2] - uDetZ[2,2] = 0
 udety2zA: udetYa[2,2] - udetZa[2,2] = 0
 udetY2Z: uDetY[2,3] - uDetZ[2,3] = 0
 udety2zA: udetYa[2,3] - udetZa[2,3] = 0
 udetY2Z: uDetY[2,4] - uDetZ[2,4] = 0
 udety2zA: udetYa[2,4] - udetZa[2,4] = 0
 udetY2Z: uDetY[2,5] - uDetZ[2,5] = 0
 udety2zA: udetYa[2,5] - udetZa[2,5] = 0
 udetY2Z: uDetY[2,6] - uDetZ[2,6] = 0
 udety2zA: udetYa[2,6] - udetZa[2,6] = 0
 udetY2Z: uDetY[2,7] - uDetZ[2,7] = 0
 udety2zA: udetYa[2,7] - udetZa[2,7] = 0
 udetY2Z: uDetY[2,8] - uDetZ[2,8] = 0
 udety2zA: udetYa[2,8] - udetZa[2,8] = 0
 udetY2Z: uDetY[2,9] - uDetZ[2,9] = 0
 udety2zA: udetYa[2,9] - udetZa[2,9] = 0
 udetY2Z: uDetY[2,10] - uDetZ[2,10] = 0
 udety2zA: udetYa[2,10] - udetZa[2,10] = 0
 udetY2Z: uDetY[2,11] - uDetZ[2,11] = 0
 udety2zA: udetYa[2,11] - udetZa[2,11] = 0
 udetY2Z: uDetY[2,12] - uDetZ[2,12] = 0
 udety2zA: udetYa[2,12] - udetZa[2,12] = 0
 udetY2Z: uDetY[2,13] - uDetZ[2,13] = 0
 udety2zA: udetYa[2,13] - udetZa[2,13] = 0
 udetY2Z: uDetY[2,14] - uDetZ[2,14] = 0
 udety2zA: udetYa[2,14] - udetZa[2,14] = 0
 udetY2Z: uDetY[2,15] - uDetZ[2,15] = 0
 udety2zA: udetYa[2,15] - udetZa[2,15] = 0
 udetY2Z: uDetY[3,0] - uDetZ[3,0] = 0
 udety2zA: udetYa[3,0] - udetZa[3,0] = 0
 udetY2Z: uDetY[3,1] - uDetZ[3,1] = 0
 udety2zA: udetYa[3,1] - udetZa[3,1] = 0
 udetY2Z: uDetY[3,2] - uDetZ[3,2] = 0
 udety2zA: udetYa[3,2] - udetZa[3,2] = 0
 udetY2Z: uDetY[3,3] - uDetZ[3,3] = 0
 udety2zA: udetYa[3,3] - udetZa[3,3] = 0
 udetY2Z: uDetY[3,4] - uDetZ[3,4] = 0
 udety2zA: udetYa[3,4] - udetZa[3,4] = 0
 udetY2Z: uDetY[3,5] - uDetZ[3,5] = 0
 udety2zA: udetYa[3,5] - udetZa[3,5] = 0
 udetY2Z: uDetY[3,6] - uDetZ[3,6] = 0
 udety2zA: udetYa[3,6] - udetZa[3,6] = 0
 udetY2Z: uDetY[3,7] - uDetZ[3,7] = 0
 udety2zA: udetYa[3,7] - udetZa[3,7] = 0
 udetY2Z: uDetY[3,8] - uDetZ[3,8] = 0
 udety2zA: udetYa[3,8] - udetZa[3,8] = 0
 udetY2Z: uDetY[3,9] - uDetZ[3,9] = 0
 udety2zA: udetYa[3,9] - udetZa[3,9] = 0
 udetY2Z: uDetY[3,10] - uDetZ[3,10] = 0
 udety2zA: udetYa[3,10] - udetZa[3,10] = 0
 udetY2Z: uDetY[3,11] - uDetZ[3,11] = 0
 udety2zA: udetYa[3,11] - udetZa[3,11] = 0
 udetY2Z: uDetY[3,12] - uDetZ[3,12] = 0
 udety2zA: udetYa[3,12] - udetZa[3,12] = 0
 udetY2Z: uDetY[3,13] - uDetZ[3,13] = 0
 udety2zA: udetYa[3,13] - udetZa[3,13] = 0
 udetY2Z: uDetY[3,14] - uDetZ[3,14] = 0
 udety2zA: udetYa[3,14] - udetZa[3,14] = 0
 udetY2Z: uDetY[3,15] - uDetZ[3,15] = 0
 udety2zA: udetYa[3,15] - udetZa[3,15] = 0
 udetSR: - uDetZ[1,0] + uDetW[1,0] = 0
 udetsrA: - udetZa[1,0] + udetWa[1,0] = 0
 udetSR: - uDetZ[1,1] + uDetW[1,1] = 0
 udetsrA: - udetZa[1,1] + udetWa[1,1] = 0
 udetSR: - uDetZ[1,2] + uDetW[1,2] = 0
 udetsrA: - udetZa[1,2] + udetWa[1,2] = 0
 udetSR: - uDetZ[1,3] + uDetW[1,3] = 0
 udetsrA: - udetZa[1,3] + udetWa[1,3] = 0
 udetSR: - uDetZ[1,7] + uDetW[1,4] = 0
 udetsrA: - udetZa[1,7] + udetWa[1,4] = 0
 udetSR: - uDetZ[1,4] + uDetW[1,5] = 0
 udetsrA: - udetZa[1,4] + udetWa[1,5] = 0
 udetSR: - uDetZ[1,5] + uDetW[1,6] = 0
 udetsrA: - udetZa[1,5] + udetWa[1,6] = 0
 udetSR: - uDetZ[1,6] + uDetW[1,7] = 0
 udetsrA: - udetZa[1,6] + udetWa[1,7] = 0
 udetSR: - uDetZ[1,10] + uDetW[1,8] = 0
 udetsrA: - udetZa[1,10] + udetWa[1,8] = 0
 udetSR: - uDetZ[1,11] + uDetW[1,9] = 0
 udetsrA: - udetZa[1,11] + udetWa[1,9] = 0
 udetSR: - uDetZ[1,8] + uDetW[1,10] = 0
 udetsrA: - udetZa[1,8] + udetWa[1,10] = 0
 udetSR: - uDetZ[1,9] + uDetW[1,11] = 0
 udetsrA: - udetZa[1,9] + udetWa[1,11] = 0
 udetSR: - uDetZ[1,13] + uDetW[1,12] = 0
 udetsrA: - udetZa[1,13] + udetWa[1,12] = 0
 udetSR: - uDetZ[1,14] + uDetW[1,13] = 0
 udetsrA: - udetZa[1,14] + udetWa[1,13] = 0
 udetSR: - uDetZ[1,15] + uDetW[1,14] = 0
 udetsrA: - udetZa[1,15] + udetWa[1,14] = 0
 udetSR: - uDetZ[1,12] + uDetW[1,15] = 0
 udetsrA: - udetZa[1,12] + udetWa[1,15] = 0
 udetSR: - uDetZ[2,0] + uDetW[2,0] = 0
 udetsrA: - udetZa[2,0] + udetWa[2,0] = 0
 udetSR: - uDetZ[2,1] + uDetW[2,1] = 0
 udetsrA: - udetZa[2,1] + udetWa[2,1] = 0
 udetSR: - uDetZ[2,2] + uDetW[2,2] = 0
 udetsrA: - udetZa[2,2] + udetWa[2,2] = 0
 udetSR: - uDetZ[2,3] + uDetW[2,3] = 0
 udetsrA: - udetZa[2,3] + udetWa[2,3] = 0
 udetSR: - uDetZ[2,7] + uDetW[2,4] = 0
 udetsrA: - udetZa[2,7] + udetWa[2,4] = 0
 udetSR: - uDetZ[2,4] + uDetW[2,5] = 0
 udetsrA: - udetZa[2,4] + udetWa[2,5] = 0
 udetSR: - uDetZ[2,5] + uDetW[2,6] = 0
 udetsrA: - udetZa[2,5] + udetWa[2,6] = 0
 udetSR: - uDetZ[2,6] + uDetW[2,7] = 0
 udetsrA: - udetZa[2,6] + udetWa[2,7] = 0
 udetSR: - uDetZ[2,10] + uDetW[2,8] = 0
 udetsrA: - udetZa[2,10] + udetWa[2,8] = 0
 udetSR: - uDetZ[2,11] + uDetW[2,9] = 0
 udetsrA: - udetZa[2,11] + udetWa[2,9] = 0
 udetSR: - uDetZ[2,8] + uDetW[2,10] = 0
 udetsrA: - udetZa[2,8] + udetWa[2,10] = 0
 udetSR: - uDetZ[2,9] + uDetW[2,11] = 0
 udetsrA: - udetZa[2,9] + udetWa[2,11] = 0
 udetSR: - uDetZ[2,13] + uDetW[2,12] = 0
 udetsrA: - udetZa[2,13] + udetWa[2,12] = 0
 udetSR: - uDetZ[2,14] + uDetW[2,13] = 0
 udetsrA: - udetZa[2,14] + udetWa[2,13] = 0
 udetSR: - uDetZ[2,15] + uDetW[2,14] = 0
 udetsrA: - udetZa[2,15] + udetWa[2,14] = 0
 udetSR: - uDetZ[2,12] + uDetW[2,15] = 0
 udetsrA: - udetZa[2,12] + udetWa[2,15] = 0
 udetSR: - uDetZ[3,0] + uDetW[3,0] = 0
 udetsrA: - udetZa[3,0] + udetWa[3,0] = 0
 udetSR: - uDetZ[3,1] + uDetW[3,1] = 0
 udetsrA: - udetZa[3,1] + udetWa[3,1] = 0
 udetSR: - uDetZ[3,2] + uDetW[3,2] = 0
 udetsrA: - udetZa[3,2] + udetWa[3,2] = 0
 udetSR: - uDetZ[3,3] + uDetW[3,3] = 0
 udetsrA: - udetZa[3,3] + udetWa[3,3] = 0
 udetSR: - uDetZ[3,7] + uDetW[3,4] = 0
 udetsrA: - udetZa[3,7] + udetWa[3,4] = 0
 udetSR: - uDetZ[3,4] + uDetW[3,5] = 0
 udetsrA: - udetZa[3,4] + udetWa[3,5] = 0
 udetSR: - uDetZ[3,5] + uDetW[3,6] = 0
 udetsrA: - udetZa[3,5] + udetWa[3,6] = 0
 udetSR: - uDetZ[3,6] + uDetW[3,7] = 0
 udetsrA: - udetZa[3,6] + udetWa[3,7] = 0
 udetSR: - uDetZ[3,10] + uDetW[3,8] = 0
 udetsrA: - udetZa[3,10] + udetWa[3,8] = 0
 udetSR: - uDetZ[3,11] + uDetW[3,9] = 0
 udetsrA: - udetZa[3,11] + udetWa[3,9] = 0
 udetSR: - uDetZ[3,8] + uDetW[3,10] = 0
 udetsrA: - udetZa[3,8] + udetWa[3,10] = 0
 udetSR: - uDetZ[3,9] + uDetW[3,11] = 0
 udetsrA: - udetZa[3,9] + udetWa[3,11] = 0
 udetSR: - uDetZ[3,13] + uDetW[3,12] = 0
 udetsrA: - udetZa[3,13] + udetWa[3,12] = 0
 udetSR: - uDetZ[3,14] + uDetW[3,13] = 0
 udetsrA: - udetZa[3,14] + udetWa[3,13] = 0
 udetSR: - uDetZ[3,15] + uDetW[3,14] = 0
 udetsrA: - udetZa[3,15] + udetWa[3,14] = 0
 udetSR: - uDetZ[3,12] + uDetW[3,15] = 0
 udetsrA: - udetZa[3,12] + udetWa[3,15] = 0
 udetMC0: uDetX[2,0] - uDetW[1,12] <= 0
 udetMC0: uDetX[2,0] - uDetX[2,12] <= 0
 udetMC0: uDetX[2,0] - uDetX[2,12] - uDetW[1,12] >= -1
 udetmc0A: udetXa[2,0] - udetWa[1,12] <= 0
 udetmc0A: udetXa[2,0] - udetXa[2,12] <= 0
 udetmc0A: udetXa[2,0] - udetXa[2,12] - udetWa[1,12] >= -1
 udetMC1: uDetX[2,4] - uDetW[1,0] = 0
 udetmc1A: udetXa[2,4] - udetWa[1,0] = 0
 udetMC2: uDetX[2,8] - uDetW[1,4] <= 0
 udetMC2: uDetX[2,8] - uDetW[1,8] <= 0
 udetMC2: uDetX[2,8] - uDetW[1,4] - uDetW[1,8] >= -1
 udetmc2A: udetXa[2,8] - udetWa[1,4] <= 0
 udetmc2A: udetXa[2,8] - udetWa[1,8] <= 0
 udetmc2A: udetXa[2,8] - udetWa[1,4] - udetWa[1,8] >= -1
 udetMC3: uDetX[2,12] - uDetW[1,0] <= 0
 udetMC3: uDetX[2,12] - uDetW[1,8] <= 0
 udetMC3: uDetX[2,12] - uDetW[1,0] - uDetW[1,8] >= -1
 udetmc3A: udetXa[2,12] - udetWa[1,0] <= 0
 udetmc3A: udetXa[2,12] - udetWa[1,8] <= 0
 udetmc3A: udetXa[2,12] - udetWa[1,0] - udetWa[1,8] >= -1
 udetMC0: uDetX[2,1] - uDetW[1,13] <= 0
 udetMC0: uDetX[2,1] - uDetX[2,13] <= 0
 udetMC0: uDetX[2,1] - uDetX[2,13] - uDetW[1,13] >= -1
 udetmc0A: udetXa[2,1] - udetWa[1,13] <= 0
 udetmc0A: udetXa[2,1] - udetXa[2,13] <= 0
 udetmc0A: udetXa[2,1] - udetXa[2,13] - udetWa[1,13] >= -1
 udetMC1: uDetX[2,5] - uDetW[1,1] = 0
 udetmc1A: udetXa[2,5] - udetWa[1,1] = 0
 udetMC2: uDetX[2,9] - uDetW[1,5] <= 0
 udetMC2: uDetX[2,9] - uDetW[1,9] <= 0
 udetMC2: uDetX[2,9] - uDetW[1,5] - uDetW[1,9] >= -1
 udetmc2A: udetXa[2,9] - udetWa[1,5] <= 0
 udetmc2A: udetXa[2,9] - udetWa[1,9] <= 0
 udetmc2A: udetXa[2,9] - udetWa[1,5] - udetWa[1,9] >= -1
 udetMC3: uDetX[2,13] - uDetW[1,1] <= 0
 udetMC3: uDetX[2,13] - uDetW[1,9] <= 0
 udetMC3: uDetX[2,13] - uDetW[1,1] - uDetW[1,9] >= -1
 udetmc3A: udetXa[2,13] - udetWa[1,1] <= 0
 udetmc3A: udetXa[2,13] - udetWa[1,9] <= 0
 udetmc3A: udetXa[2,13] - udetWa[1,1] - udetWa[1,9] >= -1
 udetMC0: uDetX[2,2] - uDetW[1,14] <= 0
 udetMC0: uDetX[2,2] - uDetX[2,14] <= 0
 udetMC0: uDetX[2,2] - uDetX[2,14] - uDetW[1,14] >= -1
 udetmc0A: udetXa[2,2] - udetWa[1,14] <= 0
 udetmc0A: udetXa[2,2] - udetXa[2,14] <= 0
 udetmc0A: udetXa[2,2] - udetXa[2,14] - udetWa[1,14] >= -1
 udetMC1: uDetX[2,6] - uDetW[1,2] = 0
 udetmc1A: udetXa[2,6] - udetWa[1,2] = 0
 udetMC2: uDetX[2,10] - uDetW[1,6] <= 0
 udetMC2: uDetX[2,10] - uDetW[1,10] <= 0
 udetMC2: uDetX[2,10] - uDetW[1,6] - uDetW[1,10] >= -1
 udetmc2A: udetXa[2,10] - udetWa[1,6] <= 0
 udetmc2A: udetXa[2,10] - udetWa[1,10] <= 0
 udetmc2A: udetXa[2,10] - udetWa[1,6] - udetWa[1,10] >= -1
 udetMC3: uDetX[2,14] - uDetW[1,2] <= 0
 udetMC3: uDetX[2,14] - uDetW[1,10] <= 0
 udetMC3: uDetX[2,14] - uDetW[1,2] - uDetW[1,10] >= -1
 udetmc3A: udetXa[2,14] - udetWa[1,2] <= 0
 udetmc3A: udetXa[2,14] - udetWa[1,10] <= 0
 udetmc3A: udetXa[2,14] - udetWa[1,2] - udetWa[1,10] >= -1
 udetMC0: uDetX[2,3] - uDetW[1,15] <= 0
 udetMC0: uDetX[2,3] - uDetX[2,15] <= 0
 udetMC0: uDetX[2,3] - uDetX[2,15] - uDetW[1,15] >= -1
 udetmc0A: udetXa[2,3] - udetWa[1,15] <= 0
 udetmc0A: udetXa[2,3] - udetXa[2,15] <= 0
 udetmc0A: udetXa[2,3] - udetXa[2,15] - udetWa[1,15] >= -1
 udetMC1: uDetX[2,7] - uDetW[1,3] = 0
 udetmc1A: udetXa[2,7] - udetWa[1,3] = 0
 udetMC2: uDetX[2,11] - uDetW[1,7] <= 0
 udetMC2: uDetX[2,11] - uDetW[1,11] <= 0
 udetMC2: uDetX[2,11] - uDetW[1,7] - uDetW[1,11] >= -1
 udetmc2A: udetXa[2,11] - udetWa[1,7] <= 0
 udetmc2A: udetXa[2,11] - udetWa[1,11] <= 0
 udetmc2A: udetXa[2,11] - udetWa[1,7] - udetWa[1,11] >= -1
 udetMC3: uDetX[2,15] - uDetW[1,3] <= 0
 udetMC3: uDetX[2,15] - uDetW[1,11] <= 0
 udetMC3: uDetX[2,15] - uDetW[1,3] - uDetW[1,11] >= -1
 udetmc3A: udetXa[2,15] - udetWa[1,3] <= 0
 udetmc3A: udetXa[2,15] - udetWa[1,11] <= 0
 udetmc3A: udetXa[2,15] - udetWa[1,3] - udetWa[1,11] >= -1
 udetMC0: uDetX[3,0] - uDetW[2,12] <= 0
 udetMC0: uDetX[3,0] - uDetX[3,12] <= 0
 udetMC0: uDetX[3,0] - uDetX[3,12] - uDetW[2,12] >= -1
 udetmc0A: udetXa[3,0] - udetWa[2,12] <= 0
 udetmc0A: udetXa[3,0] - udetXa[3,12] <= 0
 udetmc0A: udetXa[3,0] - udetXa[3,12] - udetWa[2,12] >= -1
 udetMC1: uDetX[3,4] - uDetW[2,0] = 0
 udetmc1A: udetXa[3,4] - udetWa[2,0] = 0
 udetMC2: uDetX[3,8] - uDetW[2,4] <= 0
 udetMC2: uDetX[3,8] - uDetW[2,8] <= 0
 udetMC2: uDetX[3,8] - uDetW[2,4] - uDetW[2,8] >= -1
 udetmc2A: udetXa[3,8] - udetWa[2,4] <= 0
 udetmc2A: udetXa[3,8] - udetWa[2,8] <= 0
 udetmc2A: udetXa[3,8] - udetWa[2,4] - udetWa[2,8] >= -1
 udetMC3: uDetX[3,12] - uDetW[2,0] <= 0
 udetMC3: uDetX[3,12] - uDetW[2,8] <= 0
 udetMC3: uDetX[3,12] - uDetW[2,0] - uDetW[2,8] >= -1
 udetmc3A: udetXa[3,12] - udetWa[2,0] <= 0
 udetmc3A: udetXa[3,12] - udetWa[2,8] <= 0
 udetmc3A: udetXa[3,12] - udetWa[2,0] - udetWa[2,8] >= -1
 udetMC0: uDetX[3,1] - uDetW[2,13] <= 0
 udetMC0: uDetX[3,1] - uDetX[3,13] <= 0
 udetMC0: uDetX[3,1] - uDetX[3,13] - uDetW[2,13] >= -1
 udetmc0A: udetXa[3,1] - udetWa[2,13] <= 0
 udetmc0A: udetXa[3,1] - udetXa[3,13] <= 0
 udetmc0A: udetXa[3,1] - udetXa[3,13] - udetWa[2,13] >= -1
 udetMC1: uDetX[3,5] - uDetW[2,1] = 0
 udetmc1A: udetXa[3,5] - udetWa[2,1] = 0
 udetMC2: uDetX[3,9] - uDetW[2,5] <= 0
 udetMC2: uDetX[3,9] - uDetW[2,9] <= 0
 udetMC2: uDetX[3,9] - uDetW[2,5] - uDetW[2,9] >= -1
 udetmc2A: udetXa[3,9] - udetWa[2,5] <= 0
 udetmc2A: udetXa[3,9] - udetWa[2,9] <= 0
 udetmc2A: udetXa[3,9] - udetWa[2,5] - udetWa[2,9] >= -1
 udetMC3: uDetX[3,13] - uDetW[2,1] <= 0
 udetMC3: uDetX[3,13] - uDetW[2,9] <= 0
 udetMC3: uDetX[3,13] - uDetW[2,1] - uDetW[2,9] >= -1
 udetmc3A: udetXa[3,13] - udetWa[2,1] <= 0
 udetmc3A: udetXa[3,13] - udetWa[2,9] <= 0
 udetmc3A: udetXa[3,13] - udetWa[2,1] - udetWa[2,9] >= -1
 udetMC0: uDetX[3,2] - uDetW[2,14] <= 0
 udetMC0: uDetX[3,2] - uDetX[3,14] <= 0
 udetMC0: uDetX[3,2] - uDetX[3,14] - uDetW[2,14] >= -1
 udetmc0A: udetXa[3,2] - udetWa[2,14] <= 0
 udetmc0A: udetXa[3,2] - udetXa[3,14] <= 0
 udetmc0A: udetXa[3,2] - udetXa[3,14] - udetWa[2,14] >= -1
 udetMC1: uDetX[3,6] - uDetW[2,2] = 0
 udetmc1A: udetXa[3,6] - udetWa[2,2] = 0
 udetMC2: uDetX[3,10] - uDetW[2,6] <= 0
 udetMC2: uDetX[3,10] - uDetW[2,10] <= 0
 udetMC2: uDetX[3,10] - uDetW[2,6] - uDetW[2,10] >= -1
 udetmc2A: udetXa[3,10] - udetWa[2,6] <= 0
 udetmc2A: udetXa[3,10] - udetWa[2,10] <= 0
 udetmc2A: udetXa[3,10] - udetWa[2,6] - udetWa[2,10] >= -1
 udetMC3: uDetX[3,14] - uDetW[2,2] <= 0
 udetMC3: uDetX[3,14] - uDetW[2,10] <= 0
 udetMC3: uDetX[3,14] - uDetW[2,2] - uDetW[2,10] >= -1
 udetmc3A: udetXa[3,14] - udetWa[2,2] <= 0
 udetmc3A: udetXa[3,14] - udetWa[2,10] <= 0
 udetmc3A: udetXa[3,14] - udetWa[2,2] - udetWa[2,10] >= -1
 udetMC0: uDetX[3,3] - uDetW[2,15] <= 0
 udetMC0: uDetX[3,3] - uDetX[3,15] <= 0
 udetMC0: uDetX[3,3] - uDetX[3,15] - uDetW[2,15] >= -1
 udetmc0A: udetXa[3,3] - udetWa[2,15] <= 0
 udetmc0A: udetXa[3,3] - udetXa[3,15] <= 0
 udetmc0A: udetXa[3,3] - udetXa[3,15] - udetWa[2,15] >= -1
 udetMC1: uDetX[3,7] - uDetW[2,3] = 0
 udetmc1A: udetXa[3,7] - udetWa[2,3] = 0
 udetMC2: uDetX[3,11] - uDetW[2,7] <= 0
 udetMC2: uDetX[3,11] - uDetW[2,11] <= 0
 udetMC2: uDetX[3,11] - uDetW[2,7] - uDetW[2,11] >= -1
 udetmc2A: udetXa[3,11] - udetWa[2,7] <= 0
 udetmc2A: udetXa[3,11] - udetWa[2,11] <= 0
 udetmc2A: udetXa[3,11] - udetWa[2,7] - udetWa[2,11] >= -1
 udetMC3: uDetX[3,15] - uDetW[2,3] <= 0
 udetMC3: uDetX[3,15] - uDetW[2,11] <= 0
 udetMC3: uDetX[3,15] - uDetW[2,3] - uDetW[2,11] >= -1
 udetmc3A: udetXa[3,15] - udetWa[2,3] <= 0
 udetmc3A: udetXa[3,15] - udetWa[2,11] <= 0
 udetmc3A: udetXa[3,15] - udetWa[2,3] - udetWa[2,11] >= -1
 udetMC0: uDetX[4,0] - uDetW[3,12] <= 0
 udetMC0: uDetX[4,0] - uDetX[4,12] <= 0
 udetMC0: uDetX[4,0] - uDetX[4,12] - uDetW[3,12] >= -1
 udetmc0A: udetXa[4,0] - udetWa[3,12] <= 0
 udetmc0A: udetXa[4,0] - udetXa[4,12] <= 0
 udetmc0A: udetXa[4,0] - udetXa[4,12] - udetWa[3,12] >= -1
 udetMC1: uDetX[4,4] - uDetW[3,0] = 0
 udetmc1A: udetXa[4,4] - udetWa[3,0] = 0
 udetMC2: uDetX[4,8] - uDetW[3,4] <= 0
 udetMC2: uDetX[4,8] - uDetW[3,8] <= 0
 udetMC2: uDetX[4,8] - uDetW[3,4] - uDetW[3,8] >= -1
 udetmc2A: udetXa[4,8] - udetWa[3,4] <= 0
 udetmc2A: udetXa[4,8] - udetWa[3,8] <= 0
 udetmc2A: udetXa[4,8] - udetWa[3,4] - udetWa[3,8] >= -1
 udetMC3: uDetX[4,12] - uDetW[3,0] <= 0
 udetMC3: uDetX[4,12] - uDetW[3,8] <= 0
 udetMC3: uDetX[4,12] - uDetW[3,0] - uDetW[3,8] >= -1
 udetmc3A: udetXa[4,12] - udetWa[3,0] <= 0
 udetmc3A: udetXa[4,12] - udetWa[3,8] <= 0
 udetmc3A: udetXa[4,12] - udetWa[3,0] - udetWa[3,8] >= -1
 udetMC0: uDetX[4,1] - uDetW[3,13] <= 0
 udetMC0: uDetX[4,1] - uDetX[4,13] <= 0
 udetMC0: uDetX[4,1] - uDetX[4,13] - uDetW[3,13] >= -1
 udetmc0A: udetXa[4,1] - udetWa[3,13] <= 0
 udetmc0A: udetXa[4,1] - udetXa[4,13] <= 0
 udetmc0A: udetXa[4,1] - udetXa[4,13] - udetWa[3,13] >= -1
 udetMC1: uDetX[4,5] - uDetW[3,1] = 0
 udetmc1A: udetXa[4,5] - udetWa[3,1] = 0
 udetMC2: uDetX[4,9] - uDetW[3,5] <= 0
 udetMC2: uDetX[4,9] - uDetW[3,9] <= 0
 udetMC2: uDetX[4,9] - uDetW[3,5] - uDetW[3,9] >= -1
 udetmc2A: udetXa[4,9] - udetWa[3,5] <= 0
 udetmc2A: udetXa[4,9] - udetWa[3,9] <= 0
 udetmc2A: udetXa[4,9] - udetWa[3,5] - udetWa[3,9] >= -1
 udetMC3: uDetX[4,13] - uDetW[3,1] <= 0
 udetMC3: uDetX[4,13] - uDetW[3,9] <= 0
 udetMC3: uDetX[4,13] - uDetW[3,1] - uDetW[3,9] >= -1
 udetmc3A: udetXa[4,13] - udetWa[3,1] <= 0
 udetmc3A: udetXa[4,13] - udetWa[3,9] <= 0
 udetmc3A: udetXa[4,13] - udetWa[3,1] - udetWa[3,9] >= -1
 udetMC0: uDetX[4,2] - uDetW[3,14] <= 0
 udetMC0: uDetX[4,2] - uDetX[4,14] <= 0
 udetMC0: uDetX[4,2] - uDetX[4,14] - uDetW[3,14] >= -1
 udetmc0A: udetXa[4,2] - udetWa[3,14] <= 0
 udetmc0A: udetXa[4,2] - udetXa[4,14] <= 0
 udetmc0A: udetXa[4,2] - udetXa[4,14] - udetWa[3,14] >= -1
 udetMC1: uDetX[4,6] - uDetW[3,2] = 0
 udetmc1A: udetXa[4,6] - udetWa[3,2] = 0
 udetMC2: uDetX[4,10] - uDetW[3,6] <= 0
 udetMC2: uDetX[4,10] - uDetW[3,10] <= 0
 udetMC2: uDetX[4,10] - uDetW[3,6] - uDetW[3,10] >= -1
 udetmc2A: udetXa[4,10] - udetWa[3,6] <= 0
 udetmc2A: udetXa[4,10] - udetWa[3,10] <= 0
 udetmc2A: udetXa[4,10] - udetWa[3,6] - udetWa[3,10] >= -1
 udetMC3: uDetX[4,14] - uDetW[3,2] <= 0
 udetMC3: uDetX[4,14] - uDetW[3,10] <= 0
 udetMC3: uDetX[4,14] - uDetW[3,2] - uDetW[3,10] >= -1
 udetmc3A: udetXa[4,14] - udetWa[3,2] <= 0
 udetmc3A: udetXa[4,14] - udetWa[3,10] <= 0
 udetmc3A: udetXa[4,14] - udetWa[3,2] - udetWa[3,10] >= -1
 udetMC0: uDetX[4,3] - uDetW[3,15] <= 0
 udetMC0: uDetX[4,3] - uDetX[4,15] <= 0
 udetMC0: uDetX[4,3] - uDetX[4,15] - uDetW[3,15] >= -1
 udetmc0A: udetXa[4,3] - udetWa[3,15] <= 0
 udetmc0A: udetXa[4,3] - udetXa[4,15] <= 0
 udetmc0A: udetXa[4,3] - udetXa[4,15] - udetWa[3,15] >= -1
 udetMC1: uDetX[4,7] - uDetW[3,3] = 0
 udetmc1A: udetXa[4,7] - udetWa[3,3] = 0
 udetMC2: uDetX[4,11] - uDetW[3,7] <= 0
 udetMC2: uDetX[4,11] - uDetW[3,11] <= 0
 udetMC2: uDetX[4,11] - uDetW[3,7] - uDetW[3,11] >= -1
 udetmc2A: udetXa[4,11] - udetWa[3,7] <= 0
 udetmc2A: udetXa[4,11] - udetWa[3,11] <= 0
 udetmc2A: udetXa[4,11] - udetWa[3,7] - udetWa[3,11] >= -1
 udetMC3: uDetX[4,15] - uDetW[3,3] <= 0
 udetMC3: uDetX[4,15] - uDetW[3,11] <= 0
 udetMC3: uDetX[4,15] - uDetW[3,3] - uDetW[3,11] >= -1
 udetmc3A: udetXa[4,15] - udetWa[3,3] <= 0
 udetmc3A: udetXa[4,15] - udetWa[3,11] <= 0
 udetmc3A: udetXa[4,15] - udetWa[3,3] - udetWa[3,11] >= -1
 ufrMC0: - uDetX[2,0] + uFrMC[2,0] <= 0
 ufrMC0: - udx[2,12] + uFrMC[2,0] <= 0
 ufrMC0: - udw[1,12] + uFrMC[2,0] <= 0
 ufrMC0: udx[2,0] + uFrMC[2,0] <= 1
 ufrMC0: udx[2,0] - udx[2,12] - udw[1,12] - uDetX[2,15] + uFrMC[2,0] >= -2
 ufrmc0A: - udetXa[2,0] + ufrMCa[2,0] <= 0
 ufrmc0A: - udx[2,12] + ufrMCa[2,0] <= 0
 ufrmc0A: - udw[1,12] + ufrMCa[2,0] <= 0
 ufrmc0A: udx[2,0] + ufrMCa[2,0] <= 1
 ufrmc0A: udx[2,0] - udx[2,12] - udw[1,12] - udetXa[2,15] + ufrMCa[2,0]
   >= -2
 ufrMC1: uFrMC[2,4] = 0
 ufrmc1A: ufrMCa[2,4] = 0
 ufrMC2: - uDetX[2,8] + uFrMC[2,8] <= 0
 ufrMC2: - udw[1,4] + uFrMC[2,8] <= 0
 ufrMC2: - udw[1,8] + uFrMC[2,8] <= 0
 ufrMC2: udx[2,8] + uFrMC[2,8] <= 1
 ufrMC2: udx[2,8] - udw[1,4] - udw[1,8] - uDetX[2,8] + uFrMC[2,8] >= -2
 ufrmc2A: - udetXa[2,8] + ufrMCa[2,8] <= 0
 ufrmc2A: - udw[1,4] + ufrMCa[2,8] <= 0
 ufrmc2A: - udw[1,8] + ufrMCa[2,8] <= 0
 ufrmc2A: udx[2,8] + ufrMCa[2,8] <= 1
 ufrmc2A: udx[2,8] - udw[1,4] - udw[1,8] - udetXa[2,8] + ufrMCa[2,8] >= -2
 ufrMC3: - uDetX[2,12] + uFrMC[2,12] <= 0
 ufrMC3: - udw[1,0] + uFrMC[2,12] <= 0
 ufrMC3: - udw[1,8] + uFrMC[2,12] <= 0
 ufrMC3: udx[2,12] + uFrMC[2,12] <= 1
 ufrMC3: udx[2,8] - udw[1,0] - udw[1,8] - uDetX[2,12] + uFrMC[2,12] >= -2
 ufrmc3A: - udetXa[2,12] + ufrMCa[2,12] <= 0
 ufrmc3A: - udw[1,0] + ufrMCa[2,12] <= 0
 ufrmc3A: - udw[1,8] + ufrMCa[2,12] <= 0
 ufrmc3A: udx[2,12] + ufrMCa[2,12] <= 1
 ufrmc3A: udx[2,8] - udw[1,0] - udw[1,8] - udetXa[2,12] + ufrMCa[2,12]
   >= -2
 ufrMC0: - uDetX[2,1] + uFrMC[2,1] <= 0
 ufrMC0: - udx[2,13] + uFrMC[2,1] <= 0
 ufrMC0: - udw[1,13] + uFrMC[2,1] <= 0
 ufrMC0: udx[2,1] + uFrMC[2,1] <= 1
 ufrMC0: udx[2,1] - udx[2,13] - udw[1,13] - uDetX[2,15] + uFrMC[2,1] >= -2
 ufrmc0A: - udetXa[2,1] + ufrMCa[2,1] <= 0
 ufrmc0A: - udx[2,13] + ufrMCa[2,1] <= 0
 ufrmc0A: - udw[1,13] + ufrMCa[2,1] <= 0
 ufrmc0A: udx[2,1] + ufrMCa[2,1] <= 1
 ufrmc0A: udx[2,1] - udx[2,13] - udw[1,13] - udetXa[2,15] + ufrMCa[2,1]
   >= -2
 ufrMC1: uFrMC[2,5] = 0
 ufrmc1A: ufrMCa[2,5] = 0
 ufrMC2: - uDetX[2,9] + uFrMC[2,9] <= 0
 ufrMC2: - udw[1,5] + uFrMC[2,9] <= 0
 ufrMC2: - udw[1,9] + uFrMC[2,9] <= 0
 ufrMC2: udx[2,9] + uFrMC[2,9] <= 1
 ufrMC2: udx[2,9] - udw[1,5] - udw[1,9] - uDetX[2,9] + uFrMC[2,9] >= -2
 ufrmc2A: - udetXa[2,9] + ufrMCa[2,9] <= 0
 ufrmc2A: - udw[1,5] + ufrMCa[2,9] <= 0
 ufrmc2A: - udw[1,9] + ufrMCa[2,9] <= 0
 ufrmc2A: udx[2,9] + ufrMCa[2,9] <= 1
 ufrmc2A: udx[2,9] - udw[1,5] - udw[1,9] - udetXa[2,9] + ufrMCa[2,9] >= -2
 ufrMC3: - uDetX[2,13] + uFrMC[2,13] <= 0
 ufrMC3: - udw[1,1] + uFrMC[2,13] <= 0
 ufrMC3: - udw[1,9] + uFrMC[2,13] <= 0
 ufrMC3: udx[2,13] + uFrMC[2,13] <= 1
 ufrMC3: udx[2,9] - udw[1,1] - udw[1,9] - uDetX[2,13] + uFrMC[2,13] >= -2
 ufrmc3A: - udetXa[2,13] + ufrMCa[2,13] <= 0
 ufrmc3A: - udw[1,1] + ufrMCa[2,13] <= 0
 ufrmc3A: - udw[1,9] + ufrMCa[2,13] <= 0
 ufrmc3A: udx[2,13] + ufrMCa[2,13] <= 1
 ufrmc3A: udx[2,9] - udw[1,1] - udw[1,9] - udetXa[2,13] + ufrMCa[2,13]
   >= -2
 ufrMC0: - uDetX[2,2] + uFrMC[2,2] <= 0
 ufrMC0: - udx[2,14] + uFrMC[2,2] <= 0
 ufrMC0: - udw[1,14] + uFrMC[2,2] <= 0
 ufrMC0: udx[2,2] + uFrMC[2,2] <= 1
 ufrMC0: udx[2,2] - udx[2,14] - udw[1,14] - uDetX[2,15] + uFrMC[2,2] >= -2
 ufrmc0A: - udetXa[2,2] + ufrMCa[2,2] <= 0
 ufrmc0A: - udx[2,14] + ufrMCa[2,2] <= 0
 ufrmc0A: - udw[1,14] + ufrMCa[2,2] <= 0
 ufrmc0A: udx[2,2] + ufrMCa[2,2] <= 1
 ufrmc0A: udx[2,2] - udx[2,14] - udw[1,14] - udetXa[2,15] + ufrMCa[2,2]
   >= -2
 ufrMC1: uFrMC[2,6] = 0
 ufrmc1A: ufrMCa[2,6] = 0
 ufrMC2: - uDetX[2,10] + uFrMC[2,10] <= 0
 ufrMC2: - udw[1,6] + uFrMC[2,10] <= 0
 ufrMC2: - udw[1,10] + uFrMC[2,10] <= 0
 ufrMC2: udx[2,10] + uFrMC[2,10] <= 1
 ufrMC2: udx[2,10] - udw[1,6] - udw[1,10] - uDetX[2,10] + uFrMC[2,10]
   >= -2
 ufrmc2A: - udetXa[2,10] + ufrMCa[2,10] <= 0
 ufrmc2A: - udw[1,6] + ufrMCa[2,10] <= 0
 ufrmc2A: - udw[1,10] + ufrMCa[2,10] <= 0
 ufrmc2A: udx[2,10] + ufrMCa[2,10] <= 1
 ufrmc2A: udx[2,10] - udw[1,6] - udw[1,10] - udetXa[2,10] + ufrMCa[2,10]
   >= -2
 ufrMC3: - uDetX[2,14] + uFrMC[2,14] <= 0
 ufrMC3: - udw[1,2] + uFrMC[2,14] <= 0
 ufrMC3: - udw[1,10] + uFrMC[2,14] <= 0
 ufrMC3: udx[2,14] + uFrMC[2,14] <= 1
 ufrMC3: udx[2,10] - udw[1,2] - udw[1,10] - uDetX[2,14] + uFrMC[2,14]
   >= -2
 ufrmc3A: - udetXa[2,14] + ufrMCa[2,14] <= 0
 ufrmc3A: - udw[1,2] + ufrMCa[2,14] <= 0
 ufrmc3A: - udw[1,10] + ufrMCa[2,14] <= 0
 ufrmc3A: udx[2,14] + ufrMCa[2,14] <= 1
 ufrmc3A: udx[2,10] - udw[1,2] - udw[1,10] - udetXa[2,14] + ufrMCa[2,14]
   >= -2
 ufrMC0: - uDetX[2,3] + uFrMC[2,3] <= 0
 ufrMC0: - udx[2,15] + uFrMC[2,3] <= 0
 ufrMC0: - udw[1,15] + uFrMC[2,3] <= 0
 ufrMC0: udx[2,3] + uFrMC[2,3] <= 1
 ufrMC0: udx[2,3] - udx[2,15] - udw[1,15] - uDetX[2,15] + uFrMC[2,3] >= -2
 ufrmc0A: - udetXa[2,3] + ufrMCa[2,3] <= 0
 ufrmc0A: - udx[2,15] + ufrMCa[2,3] <= 0
 ufrmc0A: - udw[1,15] + ufrMCa[2,3] <= 0
 ufrmc0A: udx[2,3] + ufrMCa[2,3] <= 1
 ufrmc0A: udx[2,3] - udx[2,15] - udw[1,15] - udetXa[2,15] + ufrMCa[2,3]
   >= -2
 ufrMC1: uFrMC[2,7] = 0
 ufrmc1A: ufrMCa[2,7] = 0
 ufrMC2: - uDetX[2,11] + uFrMC[2,11] <= 0
 ufrMC2: - udw[1,7] + uFrMC[2,11] <= 0
 ufrMC2: - udw[1,11] + uFrMC[2,11] <= 0
 ufrMC2: udx[2,11] + uFrMC[2,11] <= 1
 ufrMC2: udx[2,11] - udw[1,7] - udw[1,11] - uDetX[2,11] + uFrMC[2,11]
   >= -2
 ufrmc2A: - udetXa[2,11] + ufrMCa[2,11] <= 0
 ufrmc2A: - udw[1,7] + ufrMCa[2,11] <= 0
 ufrmc2A: - udw[1,11] + ufrMCa[2,11] <= 0
 ufrmc2A: udx[2,11] + ufrMCa[2,11] <= 1
 ufrmc2A: udx[2,11] - udw[1,7] - udw[1,11] - udetXa[2,11] + ufrMCa[2,11]
   >= -2
 ufrMC3: - uDetX[2,15] + uFrMC[2,15] <= 0
 ufrMC3: - udw[1,3] + uFrMC[2,15] <= 0
 ufrMC3: - udw[1,11] + uFrMC[2,15] <= 0
 ufrMC3: udx[2,15] + uFrMC[2,15] <= 1
 ufrMC3: udx[2,11] - udw[1,3] - udw[1,11] - uDetX[2,15] + uFrMC[2,15]
   >= -2
 ufrmc3A: - udetXa[2,15] + ufrMCa[2,15] <= 0
 ufrmc3A: - udw[1,3] + ufrMCa[2,15] <= 0
 ufrmc3A: - udw[1,11] + ufrMCa[2,15] <= 0
 ufrmc3A: udx[2,15] + ufrMCa[2,15] <= 1
 ufrmc3A: udx[2,11] - udw[1,3] - udw[1,11] - udetXa[2,15] + ufrMCa[2,15]
   >= -2
 ufrMC0: - uDetX[3,0] + uFrMC[3,0] <= 0
 ufrMC0: - udx[3,12] + uFrMC[3,0] <= 0
 ufrMC0: - udw[2,12] + uFrMC[3,0] <= 0
 ufrMC0: udx[3,0] + uFrMC[3,0] <= 1
 ufrMC0: udx[3,0] - udx[3,12] - udw[2,12] - uDetX[3,15] + uFrMC[3,0] >= -2
 ufrmc0A: - udetXa[3,0] + ufrMCa[3,0] <= 0
 ufrmc0A: - udx[3,12] + ufrMCa[3,0] <= 0
 ufrmc0A: - udw[2,12] + ufrMCa[3,0] <= 0
 ufrmc0A: udx[3,0] + ufrMCa[3,0] <= 1
 ufrmc0A: udx[3,0] - udx[3,12] - udw[2,12] - udetXa[3,15] + ufrMCa[3,0]
   >= -2
 ufrMC1: uFrMC[3,4] = 0
 ufrmc1A: ufrMCa[3,4] = 0
 ufrMC2: - uDetX[3,8] + uFrMC[3,8] <= 0
 ufrMC2: - udw[2,4] + uFrMC[3,8] <= 0
 ufrMC2: - udw[2,8] + uFrMC[3,8] <= 0
 ufrMC2: udx[3,8] + uFrMC[3,8] <= 1
 ufrMC2: udx[3,8] - udw[2,4] - udw[2,8] - uDetX[3,8] + uFrMC[3,8] >= -2
 ufrmc2A: - udetXa[3,8] + ufrMCa[3,8] <= 0
 ufrmc2A: - udw[2,4] + ufrMCa[3,8] <= 0
 ufrmc2A: - udw[2,8] + ufrMCa[3,8] <= 0
 ufrmc2A: udx[3,8] + ufrMCa[3,8] <= 1
 ufrmc2A: udx[3,8] - udw[2,4] - udw[2,8] - udetXa[3,8] + ufrMCa[3,8] >= -2
 ufrMC3: - uDetX[3,12] + uFrMC[3,12] <= 0
 ufrMC3: - udw[2,0] + uFrMC[3,12] <= 0
 ufrMC3: - udw[2,8] + uFrMC[3,12] <= 0
 ufrMC3: udx[3,12] + uFrMC[3,12] <= 1
 ufrMC3: udx[3,8] - udw[2,0] - udw[2,8] - uDetX[3,12] + uFrMC[3,12] >= -2
 ufrmc3A: - udetXa[3,12] + ufrMCa[3,12] <= 0
 ufrmc3A: - udw[2,0] + ufrMCa[3,12] <= 0
 ufrmc3A: - udw[2,8] + ufrMCa[3,12] <= 0
 ufrmc3A: udx[3,12] + ufrMCa[3,12] <= 1
 ufrmc3A: udx[3,8] - udw[2,0] - udw[2,8] - udetXa[3,12] + ufrMCa[3,12]
   >= -2
 ufrMC0: - uDetX[3,1] + uFrMC[3,1] <= 0
 ufrMC0: - udx[3,13] + uFrMC[3,1] <= 0
 ufrMC0: - udw[2,13] + uFrMC[3,1] <= 0
 ufrMC0: udx[3,1] + uFrMC[3,1] <= 1
 ufrMC0: udx[3,1] - udx[3,13] - udw[2,13] - uDetX[3,15] + uFrMC[3,1] >= -2
 ufrmc0A: - udetXa[3,1] + ufrMCa[3,1] <= 0
 ufrmc0A: - udx[3,13] + ufrMCa[3,1] <= 0
 ufrmc0A: - udw[2,13] + ufrMCa[3,1] <= 0
 ufrmc0A: udx[3,1] + ufrMCa[3,1] <= 1
 ufrmc0A: udx[3,1] - udx[3,13] - udw[2,13] - udetXa[3,15] + ufrMCa[3,1]
   >= -2
 ufrMC1: uFrMC[3,5] = 0
 ufrmc1A: ufrMCa[3,5] = 0
 ufrMC2: - uDetX[3,9] + uFrMC[3,9] <= 0
 ufrMC2: - udw[2,5] + uFrMC[3,9] <= 0
 ufrMC2: - udw[2,9] + uFrMC[3,9] <= 0
 ufrMC2: udx[3,9] + uFrMC[3,9] <= 1
 ufrMC2: udx[3,9] - udw[2,5] - udw[2,9] - uDetX[3,9] + uFrMC[3,9] >= -2
 ufrmc2A: - udetXa[3,9] + ufrMCa[3,9] <= 0
 ufrmc2A: - udw[2,5] + ufrMCa[3,9] <= 0
 ufrmc2A: - udw[2,9] + ufrMCa[3,9] <= 0
 ufrmc2A: udx[3,9] + ufrMCa[3,9] <= 1
 ufrmc2A: udx[3,9] - udw[2,5] - udw[2,9] - udetXa[3,9] + ufrMCa[3,9] >= -2
 ufrMC3: - uDetX[3,13] + uFrMC[3,13] <= 0
 ufrMC3: - udw[2,1] + uFrMC[3,13] <= 0
 ufrMC3: - udw[2,9] + uFrMC[3,13] <= 0
 ufrMC3: udx[3,13] + uFrMC[3,13] <= 1
 ufrMC3: udx[3,9] - udw[2,1] - udw[2,9] - uDetX[3,13] + uFrMC[3,13] >= -2
 ufrmc3A: - udetXa[3,13] + ufrMCa[3,13] <= 0
 ufrmc3A: - udw[2,1] + ufrMCa[3,13] <= 0
 ufrmc3A: - udw[2,9] + ufrMCa[3,13] <= 0
 ufrmc3A: udx[3,13] + ufrMCa[3,13] <= 1
 ufrmc3A: udx[3,9] - udw[2,1] - udw[2,9] - udetXa[3,13] + ufrMCa[3,13]
   >= -2
 ufrMC0: - uDetX[3,2] + uFrMC[3,2] <= 0
 ufrMC0: - udx[3,14] + uFrMC[3,2] <= 0
 ufrMC0: - udw[2,14] + uFrMC[3,2] <= 0
 ufrMC0: udx[3,2] + uFrMC[3,2] <= 1
 ufrMC0: udx[3,2] - udx[3,14] - udw[2,14] - uDetX[3,15] + uFrMC[3,2] >= -2
 ufrmc0A: - udetXa[3,2] + ufrMCa[3,2] <= 0
 ufrmc0A: - udx[3,14] + ufrMCa[3,2] <= 0
 ufrmc0A: - udw[2,14] + ufrMCa[3,2] <= 0
 ufrmc0A: udx[3,2] + ufrMCa[3,2] <= 1
 ufrmc0A: udx[3,2] - udx[3,14] - udw[2,14] - udetXa[3,15] + ufrMCa[3,2]
   >= -2
 ufrMC1: uFrMC[3,6] = 0
 ufrmc1A: ufrMCa[3,6] = 0
 ufrMC2: - uDetX[3,10] + uFrMC[3,10] <= 0
 ufrMC2: - udw[2,6] + uFrMC[3,10] <= 0
 ufrMC2: - udw[2,10] + uFrMC[3,10] <= 0
 ufrMC2: udx[3,10] + uFrMC[3,10] <= 1
 ufrMC2: udx[3,10] - udw[2,6] - udw[2,10] - uDetX[3,10] + uFrMC[3,10]
   >= -2
 ufrmc2A: - udetXa[3,10] + ufrMCa[3,10] <= 0
 ufrmc2A: - udw[2,6] + ufrMCa[3,10] <= 0
 ufrmc2A: - udw[2,10] + ufrMCa[3,10] <= 0
 ufrmc2A: udx[3,10] + ufrMCa[3,10] <= 1
 ufrmc2A: udx[3,10] - udw[2,6] - udw[2,10] - udetXa[3,10] + ufrMCa[3,10]
   >= -2
 ufrMC3: - uDetX[3,14] + uFrMC[3,14] <= 0
 ufrMC3: - udw[2,2] + uFrMC[3,14] <= 0
 ufrMC3: - udw[2,10] + uFrMC[3,14] <= 0
 ufrMC3: udx[3,14] + uFrMC[3,14] <= 1
 ufrMC3: udx[3,10] - udw[2,2] - udw[2,10] - uDetX[3,14] + uFrMC[3,14]
   >= -2
 ufrmc3A: - udetXa[3,14] + ufrMCa[3,14] <= 0
 ufrmc3A: - udw[2,2] + ufrMCa[3,14] <= 0
 ufrmc3A: - udw[2,10] + ufrMCa[3,14] <= 0
 ufrmc3A: udx[3,14] + ufrMCa[3,14] <= 1
 ufrmc3A: udx[3,10] - udw[2,2] - udw[2,10] - udetXa[3,14] + ufrMCa[3,14]
   >= -2
 ufrMC0: - uDetX[3,3] + uFrMC[3,3] <= 0
 ufrMC0: - udx[3,15] + uFrMC[3,3] <= 0
 ufrMC0: - udw[2,15] + uFrMC[3,3] <= 0
 ufrMC0: udx[3,3] + uFrMC[3,3] <= 1
 ufrMC0: udx[3,3] - udx[3,15] - udw[2,15] - uDetX[3,15] + uFrMC[3,3] >= -2
 ufrmc0A: - udetXa[3,3] + ufrMCa[3,3] <= 0
 ufrmc0A: - udx[3,15] + ufrMCa[3,3] <= 0
 ufrmc0A: - udw[2,15] + ufrMCa[3,3] <= 0
 ufrmc0A: udx[3,3] + ufrMCa[3,3] <= 1
 ufrmc0A: udx[3,3] - udx[3,15] - udw[2,15] - udetXa[3,15] + ufrMCa[3,3]
   >= -2
 ufrMC1: uFrMC[3,7] = 0
 ufrmc1A: ufrMCa[3,7] = 0
 ufrMC2: - uDetX[3,11] + uFrMC[3,11] <= 0
 ufrMC2: - udw[2,7] + uFrMC[3,11] <= 0
 ufrMC2: - udw[2,11] + uFrMC[3,11] <= 0
 ufrMC2: udx[3,11] + uFrMC[3,11] <= 1
 ufrMC2: udx[3,11] - udw[2,7] - udw[2,11] - uDetX[3,11] + uFrMC[3,11]
   >= -2
 ufrmc2A: - udetXa[3,11] + ufrMCa[3,11] <= 0
 ufrmc2A: - udw[2,7] + ufrMCa[3,11] <= 0
 ufrmc2A: - udw[2,11] + ufrMCa[3,11] <= 0
 ufrmc2A: udx[3,11] + ufrMCa[3,11] <= 1
 ufrmc2A: udx[3,11] - udw[2,7] - udw[2,11] - udetXa[3,11] + ufrMCa[3,11]
   >= -2
 ufrMC3: - uDetX[3,15] + uFrMC[3,15] <= 0
 ufrMC3: - udw[2,3] + uFrMC[3,15] <= 0
 ufrMC3: - udw[2,11] + uFrMC[3,15] <= 0
 ufrMC3: udx[3,15] + uFrMC[3,15] <= 1
 ufrMC3: udx[3,11] - udw[2,3] - udw[2,11] - uDetX[3,15] + uFrMC[3,15]
   >= -2
 ufrmc3A: - udetXa[3,15] + ufrMCa[3,15] <= 0
 ufrmc3A: - udw[2,3] + ufrMCa[3,15] <= 0
 ufrmc3A: - udw[2,11] + ufrMCa[3,15] <= 0
 ufrmc3A: udx[3,15] + ufrMCa[3,15] <= 1
 ufrmc3A: udx[3,11] - udw[2,3] - udw[2,11] - udetXa[3,15] + ufrMCa[3,15]
   >= -2
 ufrMC0: - uDetX[4,0] + uFrMC[4,0] <= 0
 ufrMC0: - udx[4,12] + uFrMC[4,0] <= 0
 ufrMC0: - udw[3,12] + uFrMC[4,0] <= 0
 ufrMC0: udx[4,0] + uFrMC[4,0] <= 1
 ufrMC0: udx[4,0] - udx[4,12] - udw[3,12] - uDetX[4,15] + uFrMC[4,0] >= -2
 ufrmc0A: - udetXa[4,0] + ufrMCa[4,0] <= 0
 ufrmc0A: - udx[4,12] + ufrMCa[4,0] <= 0
 ufrmc0A: - udw[3,12] + ufrMCa[4,0] <= 0
 ufrmc0A: udx[4,0] + ufrMCa[4,0] <= 1
 ufrmc0A: udx[4,0] - udx[4,12] - udw[3,12] - udetXa[4,15] + ufrMCa[4,0]
   >= -2
 ufrMC1: uFrMC[4,4] = 0
 ufrmc1A: ufrMCa[4,4] = 0
 ufrMC2: - uDetX[4,8] + uFrMC[4,8] <= 0
 ufrMC2: - udw[3,4] + uFrMC[4,8] <= 0
 ufrMC2: - udw[3,8] + uFrMC[4,8] <= 0
 ufrMC2: udx[4,8] + uFrMC[4,8] <= 1
 ufrMC2: udx[4,8] - udw[3,4] - udw[3,8] - uDetX[4,8] + uFrMC[4,8] >= -2
 ufrmc2A: - udetXa[4,8] + ufrMCa[4,8] <= 0
 ufrmc2A: - udw[3,4] + ufrMCa[4,8] <= 0
 ufrmc2A: - udw[3,8] + ufrMCa[4,8] <= 0
 ufrmc2A: udx[4,8] + ufrMCa[4,8] <= 1
 ufrmc2A: udx[4,8] - udw[3,4] - udw[3,8] - udetXa[4,8] + ufrMCa[4,8] >= -2
 ufrMC3: - uDetX[4,12] + uFrMC[4,12] <= 0
 ufrMC3: - udw[3,0] + uFrMC[4,12] <= 0
 ufrMC3: - udw[3,8] + uFrMC[4,12] <= 0
 ufrMC3: udx[4,12] + uFrMC[4,12] <= 1
 ufrMC3: udx[4,8] - udw[3,0] - udw[3,8] - uDetX[4,12] + uFrMC[4,12] >= -2
 ufrmc3A: - udetXa[4,12] + ufrMCa[4,12] <= 0
 ufrmc3A: - udw[3,0] + ufrMCa[4,12] <= 0
 ufrmc3A: - udw[3,8] + ufrMCa[4,12] <= 0
 ufrmc3A: udx[4,12] + ufrMCa[4,12] <= 1
 ufrmc3A: udx[4,8] - udw[3,0] - udw[3,8] - udetXa[4,12] + ufrMCa[4,12]
   >= -2
 ufrMC0: - uDetX[4,1] + uFrMC[4,1] <= 0
 ufrMC0: - udx[4,13] + uFrMC[4,1] <= 0
 ufrMC0: - udw[3,13] + uFrMC[4,1] <= 0
 ufrMC0: udx[4,1] + uFrMC[4,1] <= 1
 ufrMC0: udx[4,1] - udx[4,13] - udw[3,13] - uDetX[4,15] + uFrMC[4,1] >= -2
 ufrmc0A: - udetXa[4,1] + ufrMCa[4,1] <= 0
 ufrmc0A: - udx[4,13] + ufrMCa[4,1] <= 0
 ufrmc0A: - udw[3,13] + ufrMCa[4,1] <= 0
 ufrmc0A: udx[4,1] + ufrMCa[4,1] <= 1
 ufrmc0A: udx[4,1] - udx[4,13] - udw[3,13] - udetXa[4,15] + ufrMCa[4,1]
   >= -2
 ufrMC1: uFrMC[4,5] = 0
 ufrmc1A: ufrMCa[4,5] = 0
 ufrMC2: - uDetX[4,9] + uFrMC[4,9] <= 0
 ufrMC2: - udw[3,5] + uFrMC[4,9] <= 0
 ufrMC2: - udw[3,9] + uFrMC[4,9] <= 0
 ufrMC2: udx[4,9] + uFrMC[4,9] <= 1
 ufrMC2: udx[4,9] - udw[3,5] - udw[3,9] - uDetX[4,9] + uFrMC[4,9] >= -2
 ufrmc2A: - udetXa[4,9] + ufrMCa[4,9] <= 0
 ufrmc2A: - udw[3,5] + ufrMCa[4,9] <= 0
 ufrmc2A: - udw[3,9] + ufrMCa[4,9] <= 0
 ufrmc2A: udx[4,9] + ufrMCa[4,9] <= 1
 ufrmc2A: udx[4,9] - udw[3,5] - udw[3,9] - udetXa[4,9] + ufrMCa[4,9] >= -2
 ufrMC3: - uDetX[4,13] + uFrMC[4,13] <= 0
 ufrMC3: - udw[3,1] + uFrMC[4,13] <= 0
 ufrMC3: - udw[3,9] + uFrMC[4,13] <= 0
 ufrMC3: udx[4,13] + uFrMC[4,13] <= 1
 ufrMC3: udx[4,9] - udw[3,1] - udw[3,9] - uDetX[4,13] + uFrMC[4,13] >= -2
 ufrmc3A: - udetXa[4,13] + ufrMCa[4,13] <= 0
 ufrmc3A: - udw[3,1] + ufrMCa[4,13] <= 0
 ufrmc3A: - udw[3,9] + ufrMCa[4,13] <= 0
 ufrmc3A: udx[4,13] + ufrMCa[4,13] <= 1
 ufrmc3A: udx[4,9] - udw[3,1] - udw[3,9] - udetXa[4,13] + ufrMCa[4,13]
   >= -2
 ufrMC0: - uDetX[4,2] + uFrMC[4,2] <= 0
 ufrMC0: - udx[4,14] + uFrMC[4,2] <= 0
 ufrMC0: - udw[3,14] + uFrMC[4,2] <= 0
 ufrMC0: udx[4,2] + uFrMC[4,2] <= 1
 ufrMC0: udx[4,2] - udx[4,14] - udw[3,14] - uDetX[4,15] + uFrMC[4,2] >= -2
 ufrmc0A: - udetXa[4,2] + ufrMCa[4,2] <= 0
 ufrmc0A: - udx[4,14] + ufrMCa[4,2] <= 0
 ufrmc0A: - udw[3,14] + ufrMCa[4,2] <= 0
 ufrmc0A: udx[4,2] + ufrMCa[4,2] <= 1
 ufrmc0A: udx[4,2] - udx[4,14] - udw[3,14] - udetXa[4,15] + ufrMCa[4,2]
   >= -2
 ufrMC1: uFrMC[4,6] = 0
 ufrmc1A: ufrMCa[4,6] = 0
 ufrMC2: - uDetX[4,10] + uFrMC[4,10] <= 0
 ufrMC2: - udw[3,6] + uFrMC[4,10] <= 0
 ufrMC2: - udw[3,10] + uFrMC[4,10] <= 0
 ufrMC2: udx[4,10] + uFrMC[4,10] <= 1
 ufrMC2: udx[4,10] - udw[3,6] - udw[3,10] - uDetX[4,10] + uFrMC[4,10]
   >= -2
 ufrmc2A: - udetXa[4,10] + ufrMCa[4,10] <= 0
 ufrmc2A: - udw[3,6] + ufrMCa[4,10] <= 0
 ufrmc2A: - udw[3,10] + ufrMCa[4,10] <= 0
 ufrmc2A: udx[4,10] + ufrMCa[4,10] <= 1
 ufrmc2A: udx[4,10] - udw[3,6] - udw[3,10] - udetXa[4,10] + ufrMCa[4,10]
   >= -2
 ufrMC3: - uDetX[4,14] + uFrMC[4,14] <= 0
 ufrMC3: - udw[3,2] + uFrMC[4,14] <= 0
 ufrMC3: - udw[3,10] + uFrMC[4,14] <= 0
 ufrMC3: udx[4,14] + uFrMC[4,14] <= 1
 ufrMC3: udx[4,10] - udw[3,2] - udw[3,10] - uDetX[4,14] + uFrMC[4,14]
   >= -2
 ufrmc3A: - udetXa[4,14] + ufrMCa[4,14] <= 0
 ufrmc3A: - udw[3,2] + ufrMCa[4,14] <= 0
 ufrmc3A: - udw[3,10] + ufrMCa[4,14] <= 0
 ufrmc3A: udx[4,14] + ufrMCa[4,14] <= 1
 ufrmc3A: udx[4,10] - udw[3,2] - udw[3,10] - udetXa[4,14] + ufrMCa[4,14]
   >= -2
 ufrMC0: - uDetX[4,3] + uFrMC[4,3] <= 0
 ufrMC0: - udx[4,15] + uFrMC[4,3] <= 0
 ufrMC0: - udw[3,15] + uFrMC[4,3] <= 0
 ufrMC0: udx[4,3] + uFrMC[4,3] <= 1
 ufrMC0: udx[4,3] - udx[4,15] - udw[3,15] - uDetX[4,15] + uFrMC[4,3] >= -2
 ufrmc0A: - udetXa[4,3] + ufrMCa[4,3] <= 0
 ufrmc0A: - udx[4,15] + ufrMCa[4,3] <= 0
 ufrmc0A: - udw[3,15] + ufrMCa[4,3] <= 0
 ufrmc0A: udx[4,3] + ufrMCa[4,3] <= 1
 ufrmc0A: udx[4,3] - udx[4,15] - udw[3,15] - udetXa[4,15] + ufrMCa[4,3]
   >= -2
 ufrMC1: uFrMC[4,7] = 0
 ufrmc1A: ufrMCa[4,7] = 0
 ufrMC2: - uDetX[4,11] + uFrMC[4,11] <= 0
 ufrMC2: - udw[3,7] + uFrMC[4,11] <= 0
 ufrMC2: - udw[3,11] + uFrMC[4,11] <= 0
 ufrMC2: udx[4,11] + uFrMC[4,11] <= 1
 ufrMC2: udx[4,11] - udw[3,7] - udw[3,11] - uDetX[4,11] + uFrMC[4,11]
   >= -2
 ufrmc2A: - udetXa[4,11] + ufrMCa[4,11] <= 0
 ufrmc2A: - udw[3,7] + ufrMCa[4,11] <= 0
 ufrmc2A: - udw[3,11] + ufrMCa[4,11] <= 0
 ufrmc2A: udx[4,11] + ufrMCa[4,11] <= 1
 ufrmc2A: udx[4,11] - udw[3,7] - udw[3,11] - udetXa[4,11] + ufrMCa[4,11]
   >= -2
 ufrMC3: - uDetX[4,15] + uFrMC[4,15] <= 0
 ufrMC3: - udw[3,3] + uFrMC[4,15] <= 0
 ufrMC3: - udw[3,11] + uFrMC[4,15] <= 0
 ufrMC3: udx[4,15] + uFrMC[4,15] <= 1
 ufrMC3: udx[4,11] - udw[3,3] - udw[3,11] - uDetX[4,15] + uFrMC[4,15]
   >= -2
 ufrmc3A: - udetXa[4,15] + ufrMCa[4,15] <= 0
 ufrmc3A: - udw[3,3] + ufrMCa[4,15] <= 0
 ufrmc3A: - udw[3,11] + ufrMCa[4,15] <= 0
 ufrmc3A: udx[4,15] + ufrMCa[4,15] <= 1
 ufrmc3A: udx[4,11] - udw[3,3] - udw[3,11] - udetXa[4,15] + ufrMCa[4,15]
   >= -2
 rfUfr: udw[0,0] + udw[0,1] + udw[0,2] + udw[0,3] + udw[0,4] + udw[0,5]
   + udw[0,6] + udw[0,7] + udw[0,8] + udw[0,9] + udw[0,10] + udw[0,11]
   + udw[0,12] + udw[0,13] + udw[0,14] + udw[0,15] - ufrSCa[1,0]
   - ufrSCa[1,1] - ufrSCa[1,2] - ufrSCa[1,3] - ufrSCa[1,4] - ufrSCa[1,5]
   - ufrSCa[1,6] - ufrSCa[1,7] - ufrSCa[1,8] - ufrSCa[1,9] - ufrSCa[1,10]
   - ufrSCa[1,11] - ufrSCa[1,12] - ufrSCa[1,13] - ufrSCa[1,14]
   - ufrSCa[1,15] - ufrSCa[2,0] - ufrSCa[2,1] - ufrSCa[2,2] - ufrSCa[2,3]
   - ufrSCa[2,4] - ufrSCa[2,5] - ufrSCa[2,6] - ufrSCa[2,7] - ufrSCa[2,8]
   - ufrSCa[2,9] - ufrSCa[2,10] - ufrSCa[2,11] - ufrSCa[2,12]
   - ufrSCa[2,13] - ufrSCa[2,14] - ufrSCa[2,15] - ufrSCa[3,0] - ufrSCa[3,1]
   - ufrSCa[3,2] - ufrSCa[3,3] - ufrSCa[3,4] - ufrSCa[3,5] - ufrSCa[3,6]
   - ufrSCa[3,7] - ufrSCa[3,8] - ufrSCa[3,9] - ufrSCa[3,10] - ufrSCa[3,11]
   - ufrSCa[3,12] - ufrSCa[3,13] - ufrSCa[3,14] - ufrSCa[3,15]
   - ufrSCa[4,0] - ufrSCa[4,1] - ufrSCa[4,2] - ufrSCa[4,3] - ufrSCa[4,4]
   - ufrSCa[4,5] - ufrSCa[4,6] - ufrSCa[4,7] - ufrSCa[4,8] - ufrSCa[4,9]
   - ufrSCa[4,10] - ufrSCa[4,11] - ufrSCa[4,12] - ufrSCa[4,13]
   - ufrSCa[4,14] - ufrSCa[4,15] - ufrMCa[2,0] - ufrMCa[2,1] - ufrMCa[2,2]
   - ufrMCa[2,3] - ufrMCa[2,4] - ufrMCa[2,5] - ufrMCa[2,6] - ufrMCa[2,7]
   - ufrMCa[2,8] - ufrMCa[2,9] - ufrMCa[2,10] - ufrMCa[2,11] - ufrMCa[2,12]
   - ufrMCa[2,13] - ufrMCa[2,14] - ufrMCa[2,15] - ufrMCa[3,0] - ufrMCa[3,1]
   - ufrMCa[3,2] - ufrMCa[3,3] - ufrMCa[3,4] - ufrMCa[3,5] - ufrMCa[3,6]
   - ufrMCa[3,7] - ufrMCa[3,8] - ufrMCa[3,9] - ufrMCa[3,10] - ufrMCa[3,11]
   - ufrMCa[3,12] - ufrMCa[3,13] - ufrMCa[3,14] - ufrMCa[3,15]
   - ufrMCa[4,0] - ufrMCa[4,1] - ufrMCa[4,2] - ufrMCa[4,3] - ufrMCa[4,4]
   - ufrMCa[4,5] - ufrMCa[4,6] - ufrMCa[4,7] - ufrMCa[4,8] - ufrMCa[4,9]
   - ufrMCa[4,10] - ufrMCa[4,11] - ufrMCa[4,12] - ufrMCa[4,13]
   - ufrMCa[4,14] - ufrMCa[4,15] = 0
 R9287: lGstk[46,8] = 0
 R9288: lGstk[46,9] = 0
 R9289: lGstk[46,10] = 0
 R9290: lGstk[46,11] = 0
 R9291: lGstk[46,12] = 0
 R9292: lGstk[46,13] = 0
 R9293: lGstk[46,14] = 0
 R9294: lGstk[46,15] = 0
 R9295: lGstk[47,8] = 0
 R9296: lGstk[47,9] = 0
 R9297: lGstk[47,10] = 0
 R9298: lGstk[47,11] = 0
 R9299: lGstk[47,12] = 0
 R9300: lGstk[47,13] = 0
 R9301: lGstk[47,14] = 0
 R9302: lGstk[47,15] = 0
 R9303: lGstk[48,8] = 0
 R9304: lGstk[48,9] = 0
 R9305: lGstk[48,10] = 0
 R9306: lGstk[48,11] = 0
 R9307: lGstk[48,12] = 0
 R9308: lGstk[48,13] = 0
 R9309: lGstk[48,14] = 0
 R9310: lGstk[48,15] = 0
 R9311: lGstk[49,8] = 0
 R9312: lGstk[49,9] = 0
 R9313: lGstk[49,10] = 0
 R9314: lGstk[49,11] = 0
 R9315: lGstk[49,12] = 0
 R9316: lGstk[49,13] = 0
 R9317: lGstk[49,14] = 0
 R9318: lGstk[49,15] = 0
 R9319: lGstk[50,8] = 0
 R9320: lGstk[50,9] = 0
 R9321: lGstk[50,10] = 0
 R9322: lGstk[50,11] = 0
 R9323: lGstk[50,12] = 0
 R9324: lGstk[50,13] = 0
 R9325: lGstk[50,14] = 0
 R9326: lGstk[50,15] = 0
 R9327: lgstka[46,8] = 0
 R9328: lgstka[46,9] = 0
 R9329: lgstka[46,10] = 0
 R9330: lgstka[46,11] = 0
 R9331: lgstka[46,12] = 0
 R9332: lgstka[46,13] = 0
 R9333: lgstka[46,14] = 0
 R9334: lgstka[46,15] = 0
 R9335: lgstka[47,8] = 0
 R9336: lgstka[47,9] = 0
 R9337: lgstka[47,10] = 0
 R9338: lgstka[47,11] = 0
 R9339: lgstka[47,12] = 0
 R9340: lgstka[47,13] = 0
 R9341: lgstka[47,14] = 0
 R9342: lgstka[47,15] = 0
 R9343: lgstka[48,8] = 0
 R9344: lgstka[48,9] = 0
 R9345: lgstka[48,10] = 0
 R9346: lgstka[48,11] = 0
 R9347: lgstka[48,12] = 0
 R9348: lgstka[48,13] = 0
 R9349: lgstka[48,14] = 0
 R9350: lgstka[48,15] = 0
 R9351: lgstka[49,8] = 0
 R9352: lgstka[49,9] = 0
 R9353: lgstka[49,10] = 0
 R9354: lgstka[49,11] = 0
 R9355: lgstka[49,12] = 0
 R9356: lgstka[49,13] = 0
 R9357: lgstka[49,14] = 0
 R9358: lgstka[49,15] = 0
 R9359: lgstka[50,8] = 0
 R9360: lgstka[50,9] = 0
 R9361: lgstka[50,10] = 0
 R9362: lgstka[50,11] = 0
 R9363: lgstka[50,12] = 0
 R9364: lgstka[50,13] = 0
 R9365: lgstka[50,14] = 0
 R9366: lgstka[50,15] = 0
 VC: lDetZ[23,0] + lDetZ[23,1] + lDetZ[23,2] + lDetZ[23,3] + lDetZ[23,4]
   + lDetZ[23,5] + lDetZ[23,6] + lDetZ[23,7] + lDetZ[23,8] + lDetZ[23,9]
   + lDetZ[23,10] + lDetZ[23,11] + lDetZ[23,12] + lDetZ[23,13]
   + lDetZ[23,14] + lDetZ[23,15] = 16
 ldetSC: lDetX[19,0] - lDetY[19,0] = 0
 ldetSCa: ldetXa[19,0] - ldetYa[19,0] = 0
 lfrSC: lDX[19,0] - ldx[19,0] - lFrSC[19,0] >= 0
 lfrSC: lDetX[19,0] - lFrSC[19,0] >= 0
 lfrSC: - lDX[19,0] + ldx[19,0] - lDetX[19,0] + lFrSC[19,0] >= -1
 lfrSCa: lDX[19,0] - ldx[19,0] - lfrSCa[19,0] >= 0
 lfrSCa: ldetXa[19,0] - lfrSCa[19,0] >= 0
 lfrSCa: - lDX[19,0] + ldx[19,0] - ldetXa[19,0] + lfrSCa[19,0] >= -1
 ldetSC: lDetX[19,1] - lDetY[19,1] = 0
 ldetSCa: ldetXa[19,1] - ldetYa[19,1] = 0
 lfrSC: lDX[19,1] - ldx[19,1] - lFrSC[19,1] >= 0
 lfrSC: lDetX[19,1] - lFrSC[19,1] >= 0
 lfrSC: - lDX[19,1] + ldx[19,1] - lDetX[19,1] + lFrSC[19,1] >= -1
 lfrSCa: lDX[19,1] - ldx[19,1] - lfrSCa[19,1] >= 0
 lfrSCa: ldetXa[19,1] - lfrSCa[19,1] >= 0
 lfrSCa: - lDX[19,1] + ldx[19,1] - ldetXa[19,1] + lfrSCa[19,1] >= -1
 ldetSC: lDetX[19,2] - lDetY[19,2] = 0
 ldetSCa: ldetXa[19,2] - ldetYa[19,2] = 0
 lfrSC: lDX[19,2] - ldx[19,2] - lFrSC[19,2] >= 0
 lfrSC: lDetX[19,2] - lFrSC[19,2] >= 0
 lfrSC: - lDX[19,2] + ldx[19,2] - lDetX[19,2] + lFrSC[19,2] >= -1
 lfrSCa: lDX[19,2] - ldx[19,2] - lfrSCa[19,2] >= 0
 lfrSCa: ldetXa[19,2] - lfrSCa[19,2] >= 0
 lfrSCa: - lDX[19,2] + ldx[19,2] - ldetXa[19,2] + lfrSCa[19,2] >= -1
 ldetSC: lDetX[19,3] - lDetY[19,3] = 0
 ldetSCa: ldetXa[19,3] - ldetYa[19,3] = 0
 lfrSC: lDX[19,3] - ldx[19,3] - lFrSC[19,3] >= 0
 lfrSC: lDetX[19,3] - lFrSC[19,3] >= 0
 lfrSC: - lDX[19,3] + ldx[19,3] - lDetX[19,3] + lFrSC[19,3] >= -1
 lfrSCa: lDX[19,3] - ldx[19,3] - lfrSCa[19,3] >= 0
 lfrSCa: ldetXa[19,3] - lfrSCa[19,3] >= 0
 lfrSCa: - lDX[19,3] + ldx[19,3] - ldetXa[19,3] + lfrSCa[19,3] >= -1
 ldetSC: lDetX[19,4] - lDetY[19,4] = 0
 ldetSCa: ldetXa[19,4] - ldetYa[19,4] = 0
 lfrSC: lDX[19,4] - ldx[19,4] - lFrSC[19,4] >= 0
 lfrSC: lDetX[19,4] - lFrSC[19,4] >= 0
 lfrSC: - lDX[19,4] + ldx[19,4] - lDetX[19,4] + lFrSC[19,4] >= -1
 lfrSCa: lDX[19,4] - ldx[19,4] - lfrSCa[19,4] >= 0
 lfrSCa: ldetXa[19,4] - lfrSCa[19,4] >= 0
 lfrSCa: - lDX[19,4] + ldx[19,4] - ldetXa[19,4] + lfrSCa[19,4] >= -1
 ldetSC: lDetX[19,5] - lDetY[19,5] = 0
 ldetSCa: ldetXa[19,5] - ldetYa[19,5] = 0
 lfrSC: lDX[19,5] - ldx[19,5] - lFrSC[19,5] >= 0
 lfrSC: lDetX[19,5] - lFrSC[19,5] >= 0
 lfrSC: - lDX[19,5] + ldx[19,5] - lDetX[19,5] + lFrSC[19,5] >= -1
 lfrSCa: lDX[19,5] - ldx[19,5] - lfrSCa[19,5] >= 0
 lfrSCa: ldetXa[19,5] - lfrSCa[19,5] >= 0
 lfrSCa: - lDX[19,5] + ldx[19,5] - ldetXa[19,5] + lfrSCa[19,5] >= -1
 ldetSC: lDetX[19,6] - lDetY[19,6] = 0
 ldetSCa: ldetXa[19,6] - ldetYa[19,6] = 0
 lfrSC: lDX[19,6] - ldx[19,6] - lFrSC[19,6] >= 0
 lfrSC: lDetX[19,6] - lFrSC[19,6] >= 0
 lfrSC: - lDX[19,6] + ldx[19,6] - lDetX[19,6] + lFrSC[19,6] >= -1
 lfrSCa: lDX[19,6] - ldx[19,6] - lfrSCa[19,6] >= 0
 lfrSCa: ldetXa[19,6] - lfrSCa[19,6] >= 0
 lfrSCa: - lDX[19,6] + ldx[19,6] - ldetXa[19,6] + lfrSCa[19,6] >= -1
 ldetSC: lDetX[19,7] - lDetY[19,7] = 0
 ldetSCa: ldetXa[19,7] - ldetYa[19,7] = 0
 lfrSC: lDX[19,7] - ldx[19,7] - lFrSC[19,7] >= 0
 lfrSC: lDetX[19,7] - lFrSC[19,7] >= 0
 lfrSC: - lDX[19,7] + ldx[19,7] - lDetX[19,7] + lFrSC[19,7] >= -1
 lfrSCa: lDX[19,7] - ldx[19,7] - lfrSCa[19,7] >= 0
 lfrSCa: ldetXa[19,7] - lfrSCa[19,7] >= 0
 lfrSCa: - lDX[19,7] + ldx[19,7] - ldetXa[19,7] + lfrSCa[19,7] >= -1
 ldetSC: lDetX[19,8] - lDetY[19,8] = 0
 ldetSCa: ldetXa[19,8] - ldetYa[19,8] = 0
 lfrSC: lDX[19,8] - ldx[19,8] - lFrSC[19,8] >= 0
 lfrSC: lDetX[19,8] - lFrSC[19,8] >= 0
 lfrSC: - lDX[19,8] + ldx[19,8] - lDetX[19,8] + lFrSC[19,8] >= -1
 lfrSCa: lDX[19,8] - ldx[19,8] - lfrSCa[19,8] >= 0
 lfrSCa: ldetXa[19,8] - lfrSCa[19,8] >= 0
 lfrSCa: - lDX[19,8] + ldx[19,8] - ldetXa[19,8] + lfrSCa[19,8] >= -1
 ldetSC: lDetX[19,9] - lDetY[19,9] = 0
 ldetSCa: ldetXa[19,9] - ldetYa[19,9] = 0
 lfrSC: lDX[19,9] - ldx[19,9] - lFrSC[19,9] >= 0
 lfrSC: lDetX[19,9] - lFrSC[19,9] >= 0
 lfrSC: - lDX[19,9] + ldx[19,9] - lDetX[19,9] + lFrSC[19,9] >= -1
 lfrSCa: lDX[19,9] - ldx[19,9] - lfrSCa[19,9] >= 0
 lfrSCa: ldetXa[19,9] - lfrSCa[19,9] >= 0
 lfrSCa: - lDX[19,9] + ldx[19,9] - ldetXa[19,9] + lfrSCa[19,9] >= -1
 ldetSC: lDetX[19,10] - lDetY[19,10] = 0
 ldetSCa: ldetXa[19,10] - ldetYa[19,10] = 0
 lfrSC: lDX[19,10] - ldx[19,10] - lFrSC[19,10] >= 0
 lfrSC: lDetX[19,10] - lFrSC[19,10] >= 0
 lfrSC: - lDX[19,10] + ldx[19,10] - lDetX[19,10] + lFrSC[19,10] >= -1
 lfrSCa: lDX[19,10] - ldx[19,10] - lfrSCa[19,10] >= 0
 lfrSCa: ldetXa[19,10] - lfrSCa[19,10] >= 0
 lfrSCa: - lDX[19,10] + ldx[19,10] - ldetXa[19,10] + lfrSCa[19,10] >= -1
 ldetSC: lDetX[19,11] - lDetY[19,11] = 0
 ldetSCa: ldetXa[19,11] - ldetYa[19,11] = 0
 lfrSC: lDX[19,11] - ldx[19,11] - lFrSC[19,11] >= 0
 lfrSC: lDetX[19,11] - lFrSC[19,11] >= 0
 lfrSC: - lDX[19,11] + ldx[19,11] - lDetX[19,11] + lFrSC[19,11] >= -1
 lfrSCa: lDX[19,11] - ldx[19,11] - lfrSCa[19,11] >= 0
 lfrSCa: ldetXa[19,11] - lfrSCa[19,11] >= 0
 lfrSCa: - lDX[19,11] + ldx[19,11] - ldetXa[19,11] + lfrSCa[19,11] >= -1
 ldetSC: lDetX[19,12] - lDetY[19,12] = 0
 ldetSCa: ldetXa[19,12] - ldetYa[19,12] = 0
 lfrSC: lDX[19,12] - ldx[19,12] - lFrSC[19,12] >= 0
 lfrSC: lDetX[19,12] - lFrSC[19,12] >= 0
 lfrSC: - lDX[19,12] + ldx[19,12] - lDetX[19,12] + lFrSC[19,12] >= -1
 lfrSCa: lDX[19,12] - ldx[19,12] - lfrSCa[19,12] >= 0
 lfrSCa: ldetXa[19,12] - lfrSCa[19,12] >= 0
 lfrSCa: - lDX[19,12] + ldx[19,12] - ldetXa[19,12] + lfrSCa[19,12] >= -1
 ldetSC: lDetX[19,13] - lDetY[19,13] = 0
 ldetSCa: ldetXa[19,13] - ldetYa[19,13] = 0
 lfrSC: lDX[19,13] - ldx[19,13] - lFrSC[19,13] >= 0
 lfrSC: lDetX[19,13] - lFrSC[19,13] >= 0
 lfrSC: - lDX[19,13] + ldx[19,13] - lDetX[19,13] + lFrSC[19,13] >= -1
 lfrSCa: lDX[19,13] - ldx[19,13] - lfrSCa[19,13] >= 0
 lfrSCa: ldetXa[19,13] - lfrSCa[19,13] >= 0
 lfrSCa: - lDX[19,13] + ldx[19,13] - ldetXa[19,13] + lfrSCa[19,13] >= -1
 ldetSC: lDetX[19,14] - lDetY[19,14] = 0
 ldetSCa: ldetXa[19,14] - ldetYa[19,14] = 0
 lfrSC: lDX[19,14] - ldx[19,14] - lFrSC[19,14] >= 0
 lfrSC: lDetX[19,14] - lFrSC[19,14] >= 0
 lfrSC: - lDX[19,14] + ldx[19,14] - lDetX[19,14] + lFrSC[19,14] >= -1
 lfrSCa: lDX[19,14] - ldx[19,14] - lfrSCa[19,14] >= 0
 lfrSCa: ldetXa[19,14] - lfrSCa[19,14] >= 0
 lfrSCa: - lDX[19,14] + ldx[19,14] - ldetXa[19,14] + lfrSCa[19,14] >= -1
 ldetSC: lDetX[19,15] - lDetY[19,15] = 0
 ldetSCa: ldetXa[19,15] - ldetYa[19,15] = 0
 lfrSC: lDX[19,15] - ldx[19,15] - lFrSC[19,15] >= 0
 lfrSC: lDetX[19,15] - lFrSC[19,15] >= 0
 lfrSC: - lDX[19,15] + ldx[19,15] - lDetX[19,15] + lFrSC[19,15] >= -1
 lfrSCa: lDX[19,15] - ldx[19,15] - lfrSCa[19,15] >= 0
 lfrSCa: ldetXa[19,15] - lfrSCa[19,15] >= 0
 lfrSCa: - lDX[19,15] + ldx[19,15] - ldetXa[19,15] + lfrSCa[19,15] >= -1
 ldetSC: lDetX[20,0] - lDetY[20,0] = 0
 ldetSCa: ldetXa[20,0] - ldetYa[20,0] = 0
 lfrSC: lDX[20,0] - ldx[20,0] - lFrSC[20,0] >= 0
 lfrSC: lDetX[20,0] - lFrSC[20,0] >= 0
 lfrSC: - lDX[20,0] + ldx[20,0] - lDetX[20,0] + lFrSC[20,0] >= -1
 lfrSCa: lDX[20,0] - ldx[20,0] - lfrSCa[20,0] >= 0
 lfrSCa: ldetXa[20,0] - lfrSCa[20,0] >= 0
 lfrSCa: - lDX[20,0] + ldx[20,0] - ldetXa[20,0] + lfrSCa[20,0] >= -1
 ldetSC: lDetX[20,1] - lDetY[20,1] = 0
 ldetSCa: ldetXa[20,1] - ldetYa[20,1] = 0
 lfrSC: lDX[20,1] - ldx[20,1] - lFrSC[20,1] >= 0
 lfrSC: lDetX[20,1] - lFrSC[20,1] >= 0
 lfrSC: - lDX[20,1] + ldx[20,1] - lDetX[20,1] + lFrSC[20,1] >= -1
 lfrSCa: lDX[20,1] - ldx[20,1] - lfrSCa[20,1] >= 0
 lfrSCa: ldetXa[20,1] - lfrSCa[20,1] >= 0
 lfrSCa: - lDX[20,1] + ldx[20,1] - ldetXa[20,1] + lfrSCa[20,1] >= -1
 ldetSC: lDetX[20,2] - lDetY[20,2] = 0
 ldetSCa: ldetXa[20,2] - ldetYa[20,2] = 0
 lfrSC: lDX[20,2] - ldx[20,2] - lFrSC[20,2] >= 0
 lfrSC: lDetX[20,2] - lFrSC[20,2] >= 0
 lfrSC: - lDX[20,2] + ldx[20,2] - lDetX[20,2] + lFrSC[20,2] >= -1
 lfrSCa: lDX[20,2] - ldx[20,2] - lfrSCa[20,2] >= 0
 lfrSCa: ldetXa[20,2] - lfrSCa[20,2] >= 0
 lfrSCa: - lDX[20,2] + ldx[20,2] - ldetXa[20,2] + lfrSCa[20,2] >= -1
 ldetSC: lDetX[20,3] - lDetY[20,3] = 0
 ldetSCa: ldetXa[20,3] - ldetYa[20,3] = 0
 lfrSC: lDX[20,3] - ldx[20,3] - lFrSC[20,3] >= 0
 lfrSC: lDetX[20,3] - lFrSC[20,3] >= 0
 lfrSC: - lDX[20,3] + ldx[20,3] - lDetX[20,3] + lFrSC[20,3] >= -1
 lfrSCa: lDX[20,3] - ldx[20,3] - lfrSCa[20,3] >= 0
 lfrSCa: ldetXa[20,3] - lfrSCa[20,3] >= 0
 lfrSCa: - lDX[20,3] + ldx[20,3] - ldetXa[20,3] + lfrSCa[20,3] >= -1
 ldetSC: lDetX[20,4] - lDetY[20,4] = 0
 ldetSCa: ldetXa[20,4] - ldetYa[20,4] = 0
 lfrSC: lDX[20,4] - ldx[20,4] - lFrSC[20,4] >= 0
 lfrSC: lDetX[20,4] - lFrSC[20,4] >= 0
 lfrSC: - lDX[20,4] + ldx[20,4] - lDetX[20,4] + lFrSC[20,4] >= -1
 lfrSCa: lDX[20,4] - ldx[20,4] - lfrSCa[20,4] >= 0
 lfrSCa: ldetXa[20,4] - lfrSCa[20,4] >= 0
 lfrSCa: - lDX[20,4] + ldx[20,4] - ldetXa[20,4] + lfrSCa[20,4] >= -1
 ldetSC: lDetX[20,5] - lDetY[20,5] = 0
 ldetSCa: ldetXa[20,5] - ldetYa[20,5] = 0
 lfrSC: lDX[20,5] - ldx[20,5] - lFrSC[20,5] >= 0
 lfrSC: lDetX[20,5] - lFrSC[20,5] >= 0
 lfrSC: - lDX[20,5] + ldx[20,5] - lDetX[20,5] + lFrSC[20,5] >= -1
 lfrSCa: lDX[20,5] - ldx[20,5] - lfrSCa[20,5] >= 0
 lfrSCa: ldetXa[20,5] - lfrSCa[20,5] >= 0
 lfrSCa: - lDX[20,5] + ldx[20,5] - ldetXa[20,5] + lfrSCa[20,5] >= -1
 ldetSC: lDetX[20,6] - lDetY[20,6] = 0
 ldetSCa: ldetXa[20,6] - ldetYa[20,6] = 0
 lfrSC: lDX[20,6] - ldx[20,6] - lFrSC[20,6] >= 0
 lfrSC: lDetX[20,6] - lFrSC[20,6] >= 0
 lfrSC: - lDX[20,6] + ldx[20,6] - lDetX[20,6] + lFrSC[20,6] >= -1
 lfrSCa: lDX[20,6] - ldx[20,6] - lfrSCa[20,6] >= 0
 lfrSCa: ldetXa[20,6] - lfrSCa[20,6] >= 0
 lfrSCa: - lDX[20,6] + ldx[20,6] - ldetXa[20,6] + lfrSCa[20,6] >= -1
 ldetSC: lDetX[20,7] - lDetY[20,7] = 0
 ldetSCa: ldetXa[20,7] - ldetYa[20,7] = 0
 lfrSC: lDX[20,7] - ldx[20,7] - lFrSC[20,7] >= 0
 lfrSC: lDetX[20,7] - lFrSC[20,7] >= 0
 lfrSC: - lDX[20,7] + ldx[20,7] - lDetX[20,7] + lFrSC[20,7] >= -1
 lfrSCa: lDX[20,7] - ldx[20,7] - lfrSCa[20,7] >= 0
 lfrSCa: ldetXa[20,7] - lfrSCa[20,7] >= 0
 lfrSCa: - lDX[20,7] + ldx[20,7] - ldetXa[20,7] + lfrSCa[20,7] >= -1
 ldetSC: lDetX[20,8] - lDetY[20,8] = 0
 ldetSCa: ldetXa[20,8] - ldetYa[20,8] = 0
 lfrSC: lDX[20,8] - ldx[20,8] - lFrSC[20,8] >= 0
 lfrSC: lDetX[20,8] - lFrSC[20,8] >= 0
 lfrSC: - lDX[20,8] + ldx[20,8] - lDetX[20,8] + lFrSC[20,8] >= -1
 lfrSCa: lDX[20,8] - ldx[20,8] - lfrSCa[20,8] >= 0
 lfrSCa: ldetXa[20,8] - lfrSCa[20,8] >= 0
 lfrSCa: - lDX[20,8] + ldx[20,8] - ldetXa[20,8] + lfrSCa[20,8] >= -1
 ldetSC: lDetX[20,9] - lDetY[20,9] = 0
 ldetSCa: ldetXa[20,9] - ldetYa[20,9] = 0
 lfrSC: lDX[20,9] - ldx[20,9] - lFrSC[20,9] >= 0
 lfrSC: lDetX[20,9] - lFrSC[20,9] >= 0
 lfrSC: - lDX[20,9] + ldx[20,9] - lDetX[20,9] + lFrSC[20,9] >= -1
 lfrSCa: lDX[20,9] - ldx[20,9] - lfrSCa[20,9] >= 0
 lfrSCa: ldetXa[20,9] - lfrSCa[20,9] >= 0
 lfrSCa: - lDX[20,9] + ldx[20,9] - ldetXa[20,9] + lfrSCa[20,9] >= -1
 ldetSC: lDetX[20,10] - lDetY[20,10] = 0
 ldetSCa: ldetXa[20,10] - ldetYa[20,10] = 0
 lfrSC: lDX[20,10] - ldx[20,10] - lFrSC[20,10] >= 0
 lfrSC: lDetX[20,10] - lFrSC[20,10] >= 0
 lfrSC: - lDX[20,10] + ldx[20,10] - lDetX[20,10] + lFrSC[20,10] >= -1
 lfrSCa: lDX[20,10] - ldx[20,10] - lfrSCa[20,10] >= 0
 lfrSCa: ldetXa[20,10] - lfrSCa[20,10] >= 0
 lfrSCa: - lDX[20,10] + ldx[20,10] - ldetXa[20,10] + lfrSCa[20,10] >= -1
 ldetSC: lDetX[20,11] - lDetY[20,11] = 0
 ldetSCa: ldetXa[20,11] - ldetYa[20,11] = 0
 lfrSC: lDX[20,11] - ldx[20,11] - lFrSC[20,11] >= 0
 lfrSC: lDetX[20,11] - lFrSC[20,11] >= 0
 lfrSC: - lDX[20,11] + ldx[20,11] - lDetX[20,11] + lFrSC[20,11] >= -1
 lfrSCa: lDX[20,11] - ldx[20,11] - lfrSCa[20,11] >= 0
 lfrSCa: ldetXa[20,11] - lfrSCa[20,11] >= 0
 lfrSCa: - lDX[20,11] + ldx[20,11] - ldetXa[20,11] + lfrSCa[20,11] >= -1
 ldetSC: lDetX[20,12] - lDetY[20,12] = 0
 ldetSCa: ldetXa[20,12] - ldetYa[20,12] = 0
 lfrSC: lDX[20,12] - ldx[20,12] - lFrSC[20,12] >= 0
 lfrSC: lDetX[20,12] - lFrSC[20,12] >= 0
 lfrSC: - lDX[20,12] + ldx[20,12] - lDetX[20,12] + lFrSC[20,12] >= -1
 lfrSCa: lDX[20,12] - ldx[20,12] - lfrSCa[20,12] >= 0
 lfrSCa: ldetXa[20,12] - lfrSCa[20,12] >= 0
 lfrSCa: - lDX[20,12] + ldx[20,12] - ldetXa[20,12] + lfrSCa[20,12] >= -1
 ldetSC: lDetX[20,13] - lDetY[20,13] = 0
 ldetSCa: ldetXa[20,13] - ldetYa[20,13] = 0
 lfrSC: lDX[20,13] - ldx[20,13] - lFrSC[20,13] >= 0
 lfrSC: lDetX[20,13] - lFrSC[20,13] >= 0
 lfrSC: - lDX[20,13] + ldx[20,13] - lDetX[20,13] + lFrSC[20,13] >= -1
 lfrSCa: lDX[20,13] - ldx[20,13] - lfrSCa[20,13] >= 0
 lfrSCa: ldetXa[20,13] - lfrSCa[20,13] >= 0
 lfrSCa: - lDX[20,13] + ldx[20,13] - ldetXa[20,13] + lfrSCa[20,13] >= -1
 ldetSC: lDetX[20,14] - lDetY[20,14] = 0
 ldetSCa: ldetXa[20,14] - ldetYa[20,14] = 0
 lfrSC: lDX[20,14] - ldx[20,14] - lFrSC[20,14] >= 0
 lfrSC: lDetX[20,14] - lFrSC[20,14] >= 0
 lfrSC: - lDX[20,14] + ldx[20,14] - lDetX[20,14] + lFrSC[20,14] >= -1
 lfrSCa: lDX[20,14] - ldx[20,14] - lfrSCa[20,14] >= 0
 lfrSCa: ldetXa[20,14] - lfrSCa[20,14] >= 0
 lfrSCa: - lDX[20,14] + ldx[20,14] - ldetXa[20,14] + lfrSCa[20,14] >= -1
 ldetSC: lDetX[20,15] - lDetY[20,15] = 0
 ldetSCa: ldetXa[20,15] - ldetYa[20,15] = 0
 lfrSC: lDX[20,15] - ldx[20,15] - lFrSC[20,15] >= 0
 lfrSC: lDetX[20,15] - lFrSC[20,15] >= 0
 lfrSC: - lDX[20,15] + ldx[20,15] - lDetX[20,15] + lFrSC[20,15] >= -1
 lfrSCa: lDX[20,15] - ldx[20,15] - lfrSCa[20,15] >= 0
 lfrSCa: ldetXa[20,15] - lfrSCa[20,15] >= 0
 lfrSCa: - lDX[20,15] + ldx[20,15] - ldetXa[20,15] + lfrSCa[20,15] >= -1
 ldetSC: lDetX[21,0] - lDetY[21,0] = 0
 ldetSCa: ldetXa[21,0] - ldetYa[21,0] = 0
 lfrSC: lDX[21,0] - ldx[21,0] - lFrSC[21,0] >= 0
 lfrSC: lDetX[21,0] - lFrSC[21,0] >= 0
 lfrSC: - lDX[21,0] + ldx[21,0] - lDetX[21,0] + lFrSC[21,0] >= -1
 lfrSCa: lDX[21,0] - ldx[21,0] - lfrSCa[21,0] >= 0
 lfrSCa: ldetXa[21,0] - lfrSCa[21,0] >= 0
 lfrSCa: - lDX[21,0] + ldx[21,0] - ldetXa[21,0] + lfrSCa[21,0] >= -1
 ldetSC: lDetX[21,1] - lDetY[21,1] = 0
 ldetSCa: ldetXa[21,1] - ldetYa[21,1] = 0
 lfrSC: lDX[21,1] - ldx[21,1] - lFrSC[21,1] >= 0
 lfrSC: lDetX[21,1] - lFrSC[21,1] >= 0
 lfrSC: - lDX[21,1] + ldx[21,1] - lDetX[21,1] + lFrSC[21,1] >= -1
 lfrSCa: lDX[21,1] - ldx[21,1] - lfrSCa[21,1] >= 0
 lfrSCa: ldetXa[21,1] - lfrSCa[21,1] >= 0
 lfrSCa: - lDX[21,1] + ldx[21,1] - ldetXa[21,1] + lfrSCa[21,1] >= -1
 ldetSC: lDetX[21,2] - lDetY[21,2] = 0
 ldetSCa: ldetXa[21,2] - ldetYa[21,2] = 0
 lfrSC: lDX[21,2] - ldx[21,2] - lFrSC[21,2] >= 0
 lfrSC: lDetX[21,2] - lFrSC[21,2] >= 0
 lfrSC: - lDX[21,2] + ldx[21,2] - lDetX[21,2] + lFrSC[21,2] >= -1
 lfrSCa: lDX[21,2] - ldx[21,2] - lfrSCa[21,2] >= 0
 lfrSCa: ldetXa[21,2] - lfrSCa[21,2] >= 0
 lfrSCa: - lDX[21,2] + ldx[21,2] - ldetXa[21,2] + lfrSCa[21,2] >= -1
 ldetSC: lDetX[21,3] - lDetY[21,3] = 0
 ldetSCa: ldetXa[21,3] - ldetYa[21,3] = 0
 lfrSC: lDX[21,3] - ldx[21,3] - lFrSC[21,3] >= 0
 lfrSC: lDetX[21,3] - lFrSC[21,3] >= 0
 lfrSC: - lDX[21,3] + ldx[21,3] - lDetX[21,3] + lFrSC[21,3] >= -1
 lfrSCa: lDX[21,3] - ldx[21,3] - lfrSCa[21,3] >= 0
 lfrSCa: ldetXa[21,3] - lfrSCa[21,3] >= 0
 lfrSCa: - lDX[21,3] + ldx[21,3] - ldetXa[21,3] + lfrSCa[21,3] >= -1
 ldetSC: lDetX[21,4] - lDetY[21,4] = 0
 ldetSCa: ldetXa[21,4] - ldetYa[21,4] = 0
 lfrSC: lDX[21,4] - ldx[21,4] - lFrSC[21,4] >= 0
 lfrSC: lDetX[21,4] - lFrSC[21,4] >= 0
 lfrSC: - lDX[21,4] + ldx[21,4] - lDetX[21,4] + lFrSC[21,4] >= -1
 lfrSCa: lDX[21,4] - ldx[21,4] - lfrSCa[21,4] >= 0
 lfrSCa: ldetXa[21,4] - lfrSCa[21,4] >= 0
 lfrSCa: - lDX[21,4] + ldx[21,4] - ldetXa[21,4] + lfrSCa[21,4] >= -1
 ldetSC: lDetX[21,5] - lDetY[21,5] = 0
 ldetSCa: ldetXa[21,5] - ldetYa[21,5] = 0
 lfrSC: lDX[21,5] - ldx[21,5] - lFrSC[21,5] >= 0
 lfrSC: lDetX[21,5] - lFrSC[21,5] >= 0
 lfrSC: - lDX[21,5] + ldx[21,5] - lDetX[21,5] + lFrSC[21,5] >= -1
 lfrSCa: lDX[21,5] - ldx[21,5] - lfrSCa[21,5] >= 0
 lfrSCa: ldetXa[21,5] - lfrSCa[21,5] >= 0
 lfrSCa: - lDX[21,5] + ldx[21,5] - ldetXa[21,5] + lfrSCa[21,5] >= -1
 ldetSC: lDetX[21,6] - lDetY[21,6] = 0
 ldetSCa: ldetXa[21,6] - ldetYa[21,6] = 0
 lfrSC: lDX[21,6] - ldx[21,6] - lFrSC[21,6] >= 0
 lfrSC: lDetX[21,6] - lFrSC[21,6] >= 0
 lfrSC: - lDX[21,6] + ldx[21,6] - lDetX[21,6] + lFrSC[21,6] >= -1
 lfrSCa: lDX[21,6] - ldx[21,6] - lfrSCa[21,6] >= 0
 lfrSCa: ldetXa[21,6] - lfrSCa[21,6] >= 0
 lfrSCa: - lDX[21,6] + ldx[21,6] - ldetXa[21,6] + lfrSCa[21,6] >= -1
 ldetSC: lDetX[21,7] - lDetY[21,7] = 0
 ldetSCa: ldetXa[21,7] - ldetYa[21,7] = 0
 lfrSC: lDX[21,7] - ldx[21,7] - lFrSC[21,7] >= 0
 lfrSC: lDetX[21,7] - lFrSC[21,7] >= 0
 lfrSC: - lDX[21,7] + ldx[21,7] - lDetX[21,7] + lFrSC[21,7] >= -1
 lfrSCa: lDX[21,7] - ldx[21,7] - lfrSCa[21,7] >= 0
 lfrSCa: ldetXa[21,7] - lfrSCa[21,7] >= 0
 lfrSCa: - lDX[21,7] + ldx[21,7] - ldetXa[21,7] + lfrSCa[21,7] >= -1
 ldetSC: lDetX[21,8] - lDetY[21,8] = 0
 ldetSCa: ldetXa[21,8] - ldetYa[21,8] = 0
 lfrSC: lDX[21,8] - ldx[21,8] - lFrSC[21,8] >= 0
 lfrSC: lDetX[21,8] - lFrSC[21,8] >= 0
 lfrSC: - lDX[21,8] + ldx[21,8] - lDetX[21,8] + lFrSC[21,8] >= -1
 lfrSCa: lDX[21,8] - ldx[21,8] - lfrSCa[21,8] >= 0
 lfrSCa: ldetXa[21,8] - lfrSCa[21,8] >= 0
 lfrSCa: - lDX[21,8] + ldx[21,8] - ldetXa[21,8] + lfrSCa[21,8] >= -1
 ldetSC: lDetX[21,9] - lDetY[21,9] = 0
 ldetSCa: ldetXa[21,9] - ldetYa[21,9] = 0
 lfrSC: lDX[21,9] - ldx[21,9] - lFrSC[21,9] >= 0
 lfrSC: lDetX[21,9] - lFrSC[21,9] >= 0
 lfrSC: - lDX[21,9] + ldx[21,9] - lDetX[21,9] + lFrSC[21,9] >= -1
 lfrSCa: lDX[21,9] - ldx[21,9] - lfrSCa[21,9] >= 0
 lfrSCa: ldetXa[21,9] - lfrSCa[21,9] >= 0
 lfrSCa: - lDX[21,9] + ldx[21,9] - ldetXa[21,9] + lfrSCa[21,9] >= -1
 ldetSC: lDetX[21,10] - lDetY[21,10] = 0
 ldetSCa: ldetXa[21,10] - ldetYa[21,10] = 0
 lfrSC: lDX[21,10] - ldx[21,10] - lFrSC[21,10] >= 0
 lfrSC: lDetX[21,10] - lFrSC[21,10] >= 0
 lfrSC: - lDX[21,10] + ldx[21,10] - lDetX[21,10] + lFrSC[21,10] >= -1
 lfrSCa: lDX[21,10] - ldx[21,10] - lfrSCa[21,10] >= 0
 lfrSCa: ldetXa[21,10] - lfrSCa[21,10] >= 0
 lfrSCa: - lDX[21,10] + ldx[21,10] - ldetXa[21,10] + lfrSCa[21,10] >= -1
 ldetSC: lDetX[21,11] - lDetY[21,11] = 0
 ldetSCa: ldetXa[21,11] - ldetYa[21,11] = 0
 lfrSC: lDX[21,11] - ldx[21,11] - lFrSC[21,11] >= 0
 lfrSC: lDetX[21,11] - lFrSC[21,11] >= 0
 lfrSC: - lDX[21,11] + ldx[21,11] - lDetX[21,11] + lFrSC[21,11] >= -1
 lfrSCa: lDX[21,11] - ldx[21,11] - lfrSCa[21,11] >= 0
 lfrSCa: ldetXa[21,11] - lfrSCa[21,11] >= 0
 lfrSCa: - lDX[21,11] + ldx[21,11] - ldetXa[21,11] + lfrSCa[21,11] >= -1
 ldetSC: lDetX[21,12] - lDetY[21,12] = 0
 ldetSCa: ldetXa[21,12] - ldetYa[21,12] = 0
 lfrSC: lDX[21,12] - ldx[21,12] - lFrSC[21,12] >= 0
 lfrSC: lDetX[21,12] - lFrSC[21,12] >= 0
 lfrSC: - lDX[21,12] + ldx[21,12] - lDetX[21,12] + lFrSC[21,12] >= -1
 lfrSCa: lDX[21,12] - ldx[21,12] - lfrSCa[21,12] >= 0
 lfrSCa: ldetXa[21,12] - lfrSCa[21,12] >= 0
 lfrSCa: - lDX[21,12] + ldx[21,12] - ldetXa[21,12] + lfrSCa[21,12] >= -1
 ldetSC: lDetX[21,13] - lDetY[21,13] = 0
 ldetSCa: ldetXa[21,13] - ldetYa[21,13] = 0
 lfrSC: lDX[21,13] - ldx[21,13] - lFrSC[21,13] >= 0
 lfrSC: lDetX[21,13] - lFrSC[21,13] >= 0
 lfrSC: - lDX[21,13] + ldx[21,13] - lDetX[21,13] + lFrSC[21,13] >= -1
 lfrSCa: lDX[21,13] - ldx[21,13] - lfrSCa[21,13] >= 0
 lfrSCa: ldetXa[21,13] - lfrSCa[21,13] >= 0
 lfrSCa: - lDX[21,13] + ldx[21,13] - ldetXa[21,13] + lfrSCa[21,13] >= -1
 ldetSC: lDetX[21,14] - lDetY[21,14] = 0
 ldetSCa: ldetXa[21,14] - ldetYa[21,14] = 0
 lfrSC: lDX[21,14] - ldx[21,14] - lFrSC[21,14] >= 0
 lfrSC: lDetX[21,14] - lFrSC[21,14] >= 0
 lfrSC: - lDX[21,14] + ldx[21,14] - lDetX[21,14] + lFrSC[21,14] >= -1
 lfrSCa: lDX[21,14] - ldx[21,14] - lfrSCa[21,14] >= 0
 lfrSCa: ldetXa[21,14] - lfrSCa[21,14] >= 0
 lfrSCa: - lDX[21,14] + ldx[21,14] - ldetXa[21,14] + lfrSCa[21,14] >= -1
 ldetSC: lDetX[21,15] - lDetY[21,15] = 0
 ldetSCa: ldetXa[21,15] - ldetYa[21,15] = 0
 lfrSC: lDX[21,15] - ldx[21,15] - lFrSC[21,15] >= 0
 lfrSC: lDetX[21,15] - lFrSC[21,15] >= 0
 lfrSC: - lDX[21,15] + ldx[21,15] - lDetX[21,15] + lFrSC[21,15] >= -1
 lfrSCa: lDX[21,15] - ldx[21,15] - lfrSCa[21,15] >= 0
 lfrSCa: ldetXa[21,15] - lfrSCa[21,15] >= 0
 lfrSCa: - lDX[21,15] + ldx[21,15] - ldetXa[21,15] + lfrSCa[21,15] >= -1
 ldetSC: lDetX[22,0] - lDetY[22,0] = 0
 ldetSCa: ldetXa[22,0] - ldetYa[22,0] = 0
 lfrSC: lDX[22,0] - ldx[22,0] - lFrSC[22,0] >= 0
 lfrSC: lDetX[22,0] - lFrSC[22,0] >= 0
 lfrSC: - lDX[22,0] + ldx[22,0] - lDetX[22,0] + lFrSC[22,0] >= -1
 lfrSCa: lDX[22,0] - ldx[22,0] - lfrSCa[22,0] >= 0
 lfrSCa: ldetXa[22,0] - lfrSCa[22,0] >= 0
 lfrSCa: - lDX[22,0] + ldx[22,0] - ldetXa[22,0] + lfrSCa[22,0] >= -1
 ldetSC: lDetX[22,1] - lDetY[22,1] = 0
 ldetSCa: ldetXa[22,1] - ldetYa[22,1] = 0
 lfrSC: lDX[22,1] - ldx[22,1] - lFrSC[22,1] >= 0
 lfrSC: lDetX[22,1] - lFrSC[22,1] >= 0
 lfrSC: - lDX[22,1] + ldx[22,1] - lDetX[22,1] + lFrSC[22,1] >= -1
 lfrSCa: lDX[22,1] - ldx[22,1] - lfrSCa[22,1] >= 0
 lfrSCa: ldetXa[22,1] - lfrSCa[22,1] >= 0
 lfrSCa: - lDX[22,1] + ldx[22,1] - ldetXa[22,1] + lfrSCa[22,1] >= -1
 ldetSC: lDetX[22,2] - lDetY[22,2] = 0
 ldetSCa: ldetXa[22,2] - ldetYa[22,2] = 0
 lfrSC: lDX[22,2] - ldx[22,2] - lFrSC[22,2] >= 0
 lfrSC: lDetX[22,2] - lFrSC[22,2] >= 0
 lfrSC: - lDX[22,2] + ldx[22,2] - lDetX[22,2] + lFrSC[22,2] >= -1
 lfrSCa: lDX[22,2] - ldx[22,2] - lfrSCa[22,2] >= 0
 lfrSCa: ldetXa[22,2] - lfrSCa[22,2] >= 0
 lfrSCa: - lDX[22,2] + ldx[22,2] - ldetXa[22,2] + lfrSCa[22,2] >= -1
 ldetSC: lDetX[22,3] - lDetY[22,3] = 0
 ldetSCa: ldetXa[22,3] - ldetYa[22,3] = 0
 lfrSC: lDX[22,3] - ldx[22,3] - lFrSC[22,3] >= 0
 lfrSC: lDetX[22,3] - lFrSC[22,3] >= 0
 lfrSC: - lDX[22,3] + ldx[22,3] - lDetX[22,3] + lFrSC[22,3] >= -1
 lfrSCa: lDX[22,3] - ldx[22,3] - lfrSCa[22,3] >= 0
 lfrSCa: ldetXa[22,3] - lfrSCa[22,3] >= 0
 lfrSCa: - lDX[22,3] + ldx[22,3] - ldetXa[22,3] + lfrSCa[22,3] >= -1
 ldetSC: lDetX[22,4] - lDetY[22,4] = 0
 ldetSCa: ldetXa[22,4] - ldetYa[22,4] = 0
 lfrSC: lDX[22,4] - ldx[22,4] - lFrSC[22,4] >= 0
 lfrSC: lDetX[22,4] - lFrSC[22,4] >= 0
 lfrSC: - lDX[22,4] + ldx[22,4] - lDetX[22,4] + lFrSC[22,4] >= -1
 lfrSCa: lDX[22,4] - ldx[22,4] - lfrSCa[22,4] >= 0
 lfrSCa: ldetXa[22,4] - lfrSCa[22,4] >= 0
 lfrSCa: - lDX[22,4] + ldx[22,4] - ldetXa[22,4] + lfrSCa[22,4] >= -1
 ldetSC: lDetX[22,5] - lDetY[22,5] = 0
 ldetSCa: ldetXa[22,5] - ldetYa[22,5] = 0
 lfrSC: lDX[22,5] - ldx[22,5] - lFrSC[22,5] >= 0
 lfrSC: lDetX[22,5] - lFrSC[22,5] >= 0
 lfrSC: - lDX[22,5] + ldx[22,5] - lDetX[22,5] + lFrSC[22,5] >= -1
 lfrSCa: lDX[22,5] - ldx[22,5] - lfrSCa[22,5] >= 0
 lfrSCa: ldetXa[22,5] - lfrSCa[22,5] >= 0
 lfrSCa: - lDX[22,5] + ldx[22,5] - ldetXa[22,5] + lfrSCa[22,5] >= -1
 ldetSC: lDetX[22,6] - lDetY[22,6] = 0
 ldetSCa: ldetXa[22,6] - ldetYa[22,6] = 0
 lfrSC: lDX[22,6] - ldx[22,6] - lFrSC[22,6] >= 0
 lfrSC: lDetX[22,6] - lFrSC[22,6] >= 0
 lfrSC: - lDX[22,6] + ldx[22,6] - lDetX[22,6] + lFrSC[22,6] >= -1
 lfrSCa: lDX[22,6] - ldx[22,6] - lfrSCa[22,6] >= 0
 lfrSCa: ldetXa[22,6] - lfrSCa[22,6] >= 0
 lfrSCa: - lDX[22,6] + ldx[22,6] - ldetXa[22,6] + lfrSCa[22,6] >= -1
 ldetSC: lDetX[22,7] - lDetY[22,7] = 0
 ldetSCa: ldetXa[22,7] - ldetYa[22,7] = 0
 lfrSC: lDX[22,7] - ldx[22,7] - lFrSC[22,7] >= 0
 lfrSC: lDetX[22,7] - lFrSC[22,7] >= 0
 lfrSC: - lDX[22,7] + ldx[22,7] - lDetX[22,7] + lFrSC[22,7] >= -1
 lfrSCa: lDX[22,7] - ldx[22,7] - lfrSCa[22,7] >= 0
 lfrSCa: ldetXa[22,7] - lfrSCa[22,7] >= 0
 lfrSCa: - lDX[22,7] + ldx[22,7] - ldetXa[22,7] + lfrSCa[22,7] >= -1
 ldetSC: lDetX[22,8] - lDetY[22,8] = 0
 ldetSCa: ldetXa[22,8] - ldetYa[22,8] = 0
 lfrSC: lDX[22,8] - ldx[22,8] - lFrSC[22,8] >= 0
 lfrSC: lDetX[22,8] - lFrSC[22,8] >= 0
 lfrSC: - lDX[22,8] + ldx[22,8] - lDetX[22,8] + lFrSC[22,8] >= -1
 lfrSCa: lDX[22,8] - ldx[22,8] - lfrSCa[22,8] >= 0
 lfrSCa: ldetXa[22,8] - lfrSCa[22,8] >= 0
 lfrSCa: - lDX[22,8] + ldx[22,8] - ldetXa[22,8] + lfrSCa[22,8] >= -1
 ldetSC: lDetX[22,9] - lDetY[22,9] = 0
 ldetSCa: ldetXa[22,9] - ldetYa[22,9] = 0
 lfrSC: lDX[22,9] - ldx[22,9] - lFrSC[22,9] >= 0
 lfrSC: lDetX[22,9] - lFrSC[22,9] >= 0
 lfrSC: - lDX[22,9] + ldx[22,9] - lDetX[22,9] + lFrSC[22,9] >= -1
 lfrSCa: lDX[22,9] - ldx[22,9] - lfrSCa[22,9] >= 0
 lfrSCa: ldetXa[22,9] - lfrSCa[22,9] >= 0
 lfrSCa: - lDX[22,9] + ldx[22,9] - ldetXa[22,9] + lfrSCa[22,9] >= -1
 ldetSC: lDetX[22,10] - lDetY[22,10] = 0
 ldetSCa: ldetXa[22,10] - ldetYa[22,10] = 0
 lfrSC: lDX[22,10] - ldx[22,10] - lFrSC[22,10] >= 0
 lfrSC: lDetX[22,10] - lFrSC[22,10] >= 0
 lfrSC: - lDX[22,10] + ldx[22,10] - lDetX[22,10] + lFrSC[22,10] >= -1
 lfrSCa: lDX[22,10] - ldx[22,10] - lfrSCa[22,10] >= 0
 lfrSCa: ldetXa[22,10] - lfrSCa[22,10] >= 0
 lfrSCa: - lDX[22,10] + ldx[22,10] - ldetXa[22,10] + lfrSCa[22,10] >= -1
 ldetSC: lDetX[22,11] - lDetY[22,11] = 0
 ldetSCa: ldetXa[22,11] - ldetYa[22,11] = 0
 lfrSC: lDX[22,11] - ldx[22,11] - lFrSC[22,11] >= 0
 lfrSC: lDetX[22,11] - lFrSC[22,11] >= 0
 lfrSC: - lDX[22,11] + ldx[22,11] - lDetX[22,11] + lFrSC[22,11] >= -1
 lfrSCa: lDX[22,11] - ldx[22,11] - lfrSCa[22,11] >= 0
 lfrSCa: ldetXa[22,11] - lfrSCa[22,11] >= 0
 lfrSCa: - lDX[22,11] + ldx[22,11] - ldetXa[22,11] + lfrSCa[22,11] >= -1
 ldetSC: lDetX[22,12] - lDetY[22,12] = 0
 ldetSCa: ldetXa[22,12] - ldetYa[22,12] = 0
 lfrSC: lDX[22,12] - ldx[22,12] - lFrSC[22,12] >= 0
 lfrSC: lDetX[22,12] - lFrSC[22,12] >= 0
 lfrSC: - lDX[22,12] + ldx[22,12] - lDetX[22,12] + lFrSC[22,12] >= -1
 lfrSCa: lDX[22,12] - ldx[22,12] - lfrSCa[22,12] >= 0
 lfrSCa: ldetXa[22,12] - lfrSCa[22,12] >= 0
 lfrSCa: - lDX[22,12] + ldx[22,12] - ldetXa[22,12] + lfrSCa[22,12] >= -1
 ldetSC: lDetX[22,13] - lDetY[22,13] = 0
 ldetSCa: ldetXa[22,13] - ldetYa[22,13] = 0
 lfrSC: lDX[22,13] - ldx[22,13] - lFrSC[22,13] >= 0
 lfrSC: lDetX[22,13] - lFrSC[22,13] >= 0
 lfrSC: - lDX[22,13] + ldx[22,13] - lDetX[22,13] + lFrSC[22,13] >= -1
 lfrSCa: lDX[22,13] - ldx[22,13] - lfrSCa[22,13] >= 0
 lfrSCa: ldetXa[22,13] - lfrSCa[22,13] >= 0
 lfrSCa: - lDX[22,13] + ldx[22,13] - ldetXa[22,13] + lfrSCa[22,13] >= -1
 ldetSC: lDetX[22,14] - lDetY[22,14] = 0
 ldetSCa: ldetXa[22,14] - ldetYa[22,14] = 0
 lfrSC: lDX[22,14] - ldx[22,14] - lFrSC[22,14] >= 0
 lfrSC: lDetX[22,14] - lFrSC[22,14] >= 0
 lfrSC: - lDX[22,14] + ldx[22,14] - lDetX[22,14] + lFrSC[22,14] >= -1
 lfrSCa: lDX[22,14] - ldx[22,14] - lfrSCa[22,14] >= 0
 lfrSCa: ldetXa[22,14] - lfrSCa[22,14] >= 0
 lfrSCa: - lDX[22,14] + ldx[22,14] - ldetXa[22,14] + lfrSCa[22,14] >= -1
 ldetSC: lDetX[22,15] - lDetY[22,15] = 0
 ldetSCa: ldetXa[22,15] - ldetYa[22,15] = 0
 lfrSC: lDX[22,15] - ldx[22,15] - lFrSC[22,15] >= 0
 lfrSC: lDetX[22,15] - lFrSC[22,15] >= 0
 lfrSC: - lDX[22,15] + ldx[22,15] - lDetX[22,15] + lFrSC[22,15] >= -1
 lfrSCa: lDX[22,15] - ldx[22,15] - lfrSCa[22,15] >= 0
 lfrSCa: ldetXa[22,15] - lfrSCa[22,15] >= 0
 lfrSCa: - lDX[22,15] + ldx[22,15] - ldetXa[22,15] + lfrSCa[22,15] >= -1
 ldetSC: lDetX[23,0] - lDetY[23,0] = 0
 ldetSCa: ldetXa[23,0] - ldetYa[23,0] = 0
 lfrSC: lDX[23,0] - ldx[23,0] - lFrSC[23,0] >= 0
 lfrSC: lDetX[23,0] - lFrSC[23,0] >= 0
 lfrSC: - lDX[23,0] + ldx[23,0] - lDetX[23,0] + lFrSC[23,0] >= -1
 lfrSCa: lDX[23,0] - ldx[23,0] - lfrSCa[23,0] >= 0
 lfrSCa: ldetXa[23,0] - lfrSCa[23,0] >= 0
 lfrSCa: - lDX[23,0] + ldx[23,0] - ldetXa[23,0] + lfrSCa[23,0] >= -1
 ldetSC: lDetX[23,1] - lDetY[23,1] = 0
 ldetSCa: ldetXa[23,1] - ldetYa[23,1] = 0
 lfrSC: lDX[23,1] - ldx[23,1] - lFrSC[23,1] >= 0
 lfrSC: lDetX[23,1] - lFrSC[23,1] >= 0
 lfrSC: - lDX[23,1] + ldx[23,1] - lDetX[23,1] + lFrSC[23,1] >= -1
 lfrSCa: lDX[23,1] - ldx[23,1] - lfrSCa[23,1] >= 0
 lfrSCa: ldetXa[23,1] - lfrSCa[23,1] >= 0
 lfrSCa: - lDX[23,1] + ldx[23,1] - ldetXa[23,1] + lfrSCa[23,1] >= -1
 ldetSC: lDetX[23,2] - lDetY[23,2] = 0
 ldetSCa: ldetXa[23,2] - ldetYa[23,2] = 0
 lfrSC: lDX[23,2] - ldx[23,2] - lFrSC[23,2] >= 0
 lfrSC: lDetX[23,2] - lFrSC[23,2] >= 0
 lfrSC: - lDX[23,2] + ldx[23,2] - lDetX[23,2] + lFrSC[23,2] >= -1
 lfrSCa: lDX[23,2] - ldx[23,2] - lfrSCa[23,2] >= 0
 lfrSCa: ldetXa[23,2] - lfrSCa[23,2] >= 0
 lfrSCa: - lDX[23,2] + ldx[23,2] - ldetXa[23,2] + lfrSCa[23,2] >= -1
 ldetSC: lDetX[23,3] - lDetY[23,3] = 0
 ldetSCa: ldetXa[23,3] - ldetYa[23,3] = 0
 lfrSC: lDX[23,3] - ldx[23,3] - lFrSC[23,3] >= 0
 lfrSC: lDetX[23,3] - lFrSC[23,3] >= 0
 lfrSC: - lDX[23,3] + ldx[23,3] - lDetX[23,3] + lFrSC[23,3] >= -1
 lfrSCa: lDX[23,3] - ldx[23,3] - lfrSCa[23,3] >= 0
 lfrSCa: ldetXa[23,3] - lfrSCa[23,3] >= 0
 lfrSCa: - lDX[23,3] + ldx[23,3] - ldetXa[23,3] + lfrSCa[23,3] >= -1
 ldetSC: lDetX[23,4] - lDetY[23,4] = 0
 ldetSCa: ldetXa[23,4] - ldetYa[23,4] = 0
 lfrSC: lDX[23,4] - ldx[23,4] - lFrSC[23,4] >= 0
 lfrSC: lDetX[23,4] - lFrSC[23,4] >= 0
 lfrSC: - lDX[23,4] + ldx[23,4] - lDetX[23,4] + lFrSC[23,4] >= -1
 lfrSCa: lDX[23,4] - ldx[23,4] - lfrSCa[23,4] >= 0
 lfrSCa: ldetXa[23,4] - lfrSCa[23,4] >= 0
 lfrSCa: - lDX[23,4] + ldx[23,4] - ldetXa[23,4] + lfrSCa[23,4] >= -1
 ldetSC: lDetX[23,5] - lDetY[23,5] = 0
 ldetSCa: ldetXa[23,5] - ldetYa[23,5] = 0
 lfrSC: lDX[23,5] - ldx[23,5] - lFrSC[23,5] >= 0
 lfrSC: lDetX[23,5] - lFrSC[23,5] >= 0
 lfrSC: - lDX[23,5] + ldx[23,5] - lDetX[23,5] + lFrSC[23,5] >= -1
 lfrSCa: lDX[23,5] - ldx[23,5] - lfrSCa[23,5] >= 0
 lfrSCa: ldetXa[23,5] - lfrSCa[23,5] >= 0
 lfrSCa: - lDX[23,5] + ldx[23,5] - ldetXa[23,5] + lfrSCa[23,5] >= -1
 ldetSC: lDetX[23,6] - lDetY[23,6] = 0
 ldetSCa: ldetXa[23,6] - ldetYa[23,6] = 0
 lfrSC: lDX[23,6] - ldx[23,6] - lFrSC[23,6] >= 0
 lfrSC: lDetX[23,6] - lFrSC[23,6] >= 0
 lfrSC: - lDX[23,6] + ldx[23,6] - lDetX[23,6] + lFrSC[23,6] >= -1
 lfrSCa: lDX[23,6] - ldx[23,6] - lfrSCa[23,6] >= 0
 lfrSCa: ldetXa[23,6] - lfrSCa[23,6] >= 0
 lfrSCa: - lDX[23,6] + ldx[23,6] - ldetXa[23,6] + lfrSCa[23,6] >= -1
 ldetSC: lDetX[23,7] - lDetY[23,7] = 0
 ldetSCa: ldetXa[23,7] - ldetYa[23,7] = 0
 lfrSC: lDX[23,7] - ldx[23,7] - lFrSC[23,7] >= 0
 lfrSC: lDetX[23,7] - lFrSC[23,7] >= 0
 lfrSC: - lDX[23,7] + ldx[23,7] - lDetX[23,7] + lFrSC[23,7] >= -1
 lfrSCa: lDX[23,7] - ldx[23,7] - lfrSCa[23,7] >= 0
 lfrSCa: ldetXa[23,7] - lfrSCa[23,7] >= 0
 lfrSCa: - lDX[23,7] + ldx[23,7] - ldetXa[23,7] + lfrSCa[23,7] >= -1
 ldetSC: lDetX[23,8] - lDetY[23,8] = 0
 ldetSCa: ldetXa[23,8] - ldetYa[23,8] = 0
 lfrSC: lDX[23,8] - ldx[23,8] - lFrSC[23,8] >= 0
 lfrSC: lDetX[23,8] - lFrSC[23,8] >= 0
 lfrSC: - lDX[23,8] + ldx[23,8] - lDetX[23,8] + lFrSC[23,8] >= -1
 lfrSCa: lDX[23,8] - ldx[23,8] - lfrSCa[23,8] >= 0
 lfrSCa: ldetXa[23,8] - lfrSCa[23,8] >= 0
 lfrSCa: - lDX[23,8] + ldx[23,8] - ldetXa[23,8] + lfrSCa[23,8] >= -1
 ldetSC: lDetX[23,9] - lDetY[23,9] = 0
 ldetSCa: ldetXa[23,9] - ldetYa[23,9] = 0
 lfrSC: lDX[23,9] - ldx[23,9] - lFrSC[23,9] >= 0
 lfrSC: lDetX[23,9] - lFrSC[23,9] >= 0
 lfrSC: - lDX[23,9] + ldx[23,9] - lDetX[23,9] + lFrSC[23,9] >= -1
 lfrSCa: lDX[23,9] - ldx[23,9] - lfrSCa[23,9] >= 0
 lfrSCa: ldetXa[23,9] - lfrSCa[23,9] >= 0
 lfrSCa: - lDX[23,9] + ldx[23,9] - ldetXa[23,9] + lfrSCa[23,9] >= -1
 ldetSC: lDetX[23,10] - lDetY[23,10] = 0
 ldetSCa: ldetXa[23,10] - ldetYa[23,10] = 0
 lfrSC: lDX[23,10] - ldx[23,10] - lFrSC[23,10] >= 0
 lfrSC: lDetX[23,10] - lFrSC[23,10] >= 0
 lfrSC: - lDX[23,10] + ldx[23,10] - lDetX[23,10] + lFrSC[23,10] >= -1
 lfrSCa: lDX[23,10] - ldx[23,10] - lfrSCa[23,10] >= 0
 lfrSCa: ldetXa[23,10] - lfrSCa[23,10] >= 0
 lfrSCa: - lDX[23,10] + ldx[23,10] - ldetXa[23,10] + lfrSCa[23,10] >= -1
 ldetSC: lDetX[23,11] - lDetY[23,11] = 0
 ldetSCa: ldetXa[23,11] - ldetYa[23,11] = 0
 lfrSC: lDX[23,11] - ldx[23,11] - lFrSC[23,11] >= 0
 lfrSC: lDetX[23,11] - lFrSC[23,11] >= 0
 lfrSC: - lDX[23,11] + ldx[23,11] - lDetX[23,11] + lFrSC[23,11] >= -1
 lfrSCa: lDX[23,11] - ldx[23,11] - lfrSCa[23,11] >= 0
 lfrSCa: ldetXa[23,11] - lfrSCa[23,11] >= 0
 lfrSCa: - lDX[23,11] + ldx[23,11] - ldetXa[23,11] + lfrSCa[23,11] >= -1
 ldetSC: lDetX[23,12] - lDetY[23,12] = 0
 ldetSCa: ldetXa[23,12] - ldetYa[23,12] = 0
 lfrSC: lDX[23,12] - ldx[23,12] - lFrSC[23,12] >= 0
 lfrSC: lDetX[23,12] - lFrSC[23,12] >= 0
 lfrSC: - lDX[23,12] + ldx[23,12] - lDetX[23,12] + lFrSC[23,12] >= -1
 lfrSCa: lDX[23,12] - ldx[23,12] - lfrSCa[23,12] >= 0
 lfrSCa: ldetXa[23,12] - lfrSCa[23,12] >= 0
 lfrSCa: - lDX[23,12] + ldx[23,12] - ldetXa[23,12] + lfrSCa[23,12] >= -1
 ldetSC: lDetX[23,13] - lDetY[23,13] = 0
 ldetSCa: ldetXa[23,13] - ldetYa[23,13] = 0
 lfrSC: lDX[23,13] - ldx[23,13] - lFrSC[23,13] >= 0
 lfrSC: lDetX[23,13] - lFrSC[23,13] >= 0
 lfrSC: - lDX[23,13] + ldx[23,13] - lDetX[23,13] + lFrSC[23,13] >= -1
 lfrSCa: lDX[23,13] - ldx[23,13] - lfrSCa[23,13] >= 0
 lfrSCa: ldetXa[23,13] - lfrSCa[23,13] >= 0
 lfrSCa: - lDX[23,13] + ldx[23,13] - ldetXa[23,13] + lfrSCa[23,13] >= -1
 ldetSC: lDetX[23,14] - lDetY[23,14] = 0
 ldetSCa: ldetXa[23,14] - ldetYa[23,14] = 0
 lfrSC: lDX[23,14] - ldx[23,14] - lFrSC[23,14] >= 0
 lfrSC: lDetX[23,14] - lFrSC[23,14] >= 0
 lfrSC: - lDX[23,14] + ldx[23,14] - lDetX[23,14] + lFrSC[23,14] >= -1
 lfrSCa: lDX[23,14] - ldx[23,14] - lfrSCa[23,14] >= 0
 lfrSCa: ldetXa[23,14] - lfrSCa[23,14] >= 0
 lfrSCa: - lDX[23,14] + ldx[23,14] - ldetXa[23,14] + lfrSCa[23,14] >= -1
 ldetSC: lDetX[23,15] - lDetY[23,15] = 0
 ldetSCa: ldetXa[23,15] - ldetYa[23,15] = 0
 lfrSC: lDX[23,15] - ldx[23,15] - lFrSC[23,15] >= 0
 lfrSC: lDetX[23,15] - lFrSC[23,15] >= 0
 lfrSC: - lDX[23,15] + ldx[23,15] - lDetX[23,15] + lFrSC[23,15] >= -1
 lfrSCa: lDX[23,15] - ldx[23,15] - lfrSCa[23,15] >= 0
 lfrSCa: ldetXa[23,15] - lfrSCa[23,15] >= 0
 lfrSCa: - lDX[23,15] + ldx[23,15] - ldetXa[23,15] + lfrSCa[23,15] >= -1
 lgkiU: - lGstk[46,0] + lDetY[19,0] <= 0
 lgkiU: lDetY[19,0] - lDetZ[19,0] <= 0
 lgkiU: - lGstk[46,0] + lDetY[19,0] - lDetZ[19,0] >= -1
 lgkiUa: - lgstka[46,0] + ldetYa[19,0] <= 0
 lgkiUa: ldetYa[19,0] - ldetZa[19,0] <= 0
 lgkiUa: - lgstka[46,0] + ldetYa[19,0] - ldetZa[19,0] >= -1
 lgkiU: - lGstk[46,1] + lDetY[19,1] <= 0
 lgkiU: lDetY[19,1] - lDetZ[19,1] <= 0
 lgkiU: - lGstk[46,1] + lDetY[19,1] - lDetZ[19,1] >= -1
 lgkiUa: - lgstka[46,1] + ldetYa[19,1] <= 0
 lgkiUa: ldetYa[19,1] - ldetZa[19,1] <= 0
 lgkiUa: - lgstka[46,1] + ldetYa[19,1] - ldetZa[19,1] >= -1
 lgkiU: - lGstk[46,2] + lDetY[19,2] <= 0
 lgkiU: lDetY[19,2] - lDetZ[19,2] <= 0
 lgkiU: - lGstk[46,2] + lDetY[19,2] - lDetZ[19,2] >= -1
 lgkiUa: - lgstka[46,2] + ldetYa[19,2] <= 0
 lgkiUa: ldetYa[19,2] - ldetZa[19,2] <= 0
 lgkiUa: - lgstka[46,2] + ldetYa[19,2] - ldetZa[19,2] >= -1
 lgkiU: - lGstk[46,3] + lDetY[19,3] <= 0
 lgkiU: lDetY[19,3] - lDetZ[19,3] <= 0
 lgkiU: - lGstk[46,3] + lDetY[19,3] - lDetZ[19,3] >= -1
 lgkiUa: - lgstka[46,3] + ldetYa[19,3] <= 0
 lgkiUa: ldetYa[19,3] - ldetZa[19,3] <= 0
 lgkiUa: - lgstka[46,3] + ldetYa[19,3] - ldetZa[19,3] >= -1
 lgkiU: - lGstk[46,4] + lDetY[19,4] <= 0
 lgkiU: lDetY[19,4] - lDetZ[19,4] <= 0
 lgkiU: - lGstk[46,4] + lDetY[19,4] - lDetZ[19,4] >= -1
 lgkiUa: - lgstka[46,4] + ldetYa[19,4] <= 0
 lgkiUa: ldetYa[19,4] - ldetZa[19,4] <= 0
 lgkiUa: - lgstka[46,4] + ldetYa[19,4] - ldetZa[19,4] >= -1
 lgkiU: - lGstk[46,5] + lDetY[19,5] <= 0
 lgkiU: lDetY[19,5] - lDetZ[19,5] <= 0
 lgkiU: - lGstk[46,5] + lDetY[19,5] - lDetZ[19,5] >= -1
 lgkiUa: - lgstka[46,5] + ldetYa[19,5] <= 0
 lgkiUa: ldetYa[19,5] - ldetZa[19,5] <= 0
 lgkiUa: - lgstka[46,5] + ldetYa[19,5] - ldetZa[19,5] >= -1
 lgkiU: - lGstk[46,6] + lDetY[19,6] <= 0
 lgkiU: lDetY[19,6] - lDetZ[19,6] <= 0
 lgkiU: - lGstk[46,6] + lDetY[19,6] - lDetZ[19,6] >= -1
 lgkiUa: - lgstka[46,6] + ldetYa[19,6] <= 0
 lgkiUa: ldetYa[19,6] - ldetZa[19,6] <= 0
 lgkiUa: - lgstka[46,6] + ldetYa[19,6] - ldetZa[19,6] >= -1
 lgkiU: - lGstk[46,7] + lDetY[19,7] <= 0
 lgkiU: lDetY[19,7] - lDetZ[19,7] <= 0
 lgkiU: - lGstk[46,7] + lDetY[19,7] - lDetZ[19,7] >= -1
 lgkiUa: - lgstka[46,7] + ldetYa[19,7] <= 0
 lgkiUa: ldetYa[19,7] - ldetZa[19,7] <= 0
 lgkiUa: - lgstka[46,7] + ldetYa[19,7] - ldetZa[19,7] >= -1
 lgkiL: lDetY[19,8] - lDetZ[19,8] = 0
 lgkiLa: ldetYa[19,8] - ldetZa[19,8] = 0
 lgkiL: lDetY[19,9] - lDetZ[19,9] = 0
 lgkiLa: ldetYa[19,9] - ldetZa[19,9] = 0
 lgkiL: lDetY[19,10] - lDetZ[19,10] = 0
 lgkiLa: ldetYa[19,10] - ldetZa[19,10] = 0
 lgkiL: lDetY[19,11] - lDetZ[19,11] = 0
 lgkiLa: ldetYa[19,11] - ldetZa[19,11] = 0
 lgkiL: lDetY[19,12] - lDetZ[19,12] = 0
 lgkiLa: ldetYa[19,12] - ldetZa[19,12] = 0
 lgkiL: lDetY[19,13] - lDetZ[19,13] = 0
 lgkiLa: ldetYa[19,13] - ldetZa[19,13] = 0
 lgkiL: lDetY[19,14] - lDetZ[19,14] = 0
 lgkiLa: ldetYa[19,14] - ldetZa[19,14] = 0
 lgkiL: lDetY[19,15] - lDetZ[19,15] = 0
 lgkiLa: ldetYa[19,15] - ldetZa[19,15] = 0
 lgkiU: - lGstk[47,0] + lDetY[20,0] <= 0
 lgkiU: lDetY[20,0] - lDetZ[20,0] <= 0
 lgkiU: - lGstk[47,0] + lDetY[20,0] - lDetZ[20,0] >= -1
 lgkiUa: - lgstka[47,0] + ldetYa[20,0] <= 0
 lgkiUa: ldetYa[20,0] - ldetZa[20,0] <= 0
 lgkiUa: - lgstka[47,0] + ldetYa[20,0] - ldetZa[20,0] >= -1
 lgkiU: - lGstk[47,1] + lDetY[20,1] <= 0
 lgkiU: lDetY[20,1] - lDetZ[20,1] <= 0
 lgkiU: - lGstk[47,1] + lDetY[20,1] - lDetZ[20,1] >= -1
 lgkiUa: - lgstka[47,1] + ldetYa[20,1] <= 0
 lgkiUa: ldetYa[20,1] - ldetZa[20,1] <= 0
 lgkiUa: - lgstka[47,1] + ldetYa[20,1] - ldetZa[20,1] >= -1
 lgkiU: - lGstk[47,2] + lDetY[20,2] <= 0
 lgkiU: lDetY[20,2] - lDetZ[20,2] <= 0
 lgkiU: - lGstk[47,2] + lDetY[20,2] - lDetZ[20,2] >= -1
 lgkiUa: - lgstka[47,2] + ldetYa[20,2] <= 0
 lgkiUa: ldetYa[20,2] - ldetZa[20,2] <= 0
 lgkiUa: - lgstka[47,2] + ldetYa[20,2] - ldetZa[20,2] >= -1
 lgkiU: - lGstk[47,3] + lDetY[20,3] <= 0
 lgkiU: lDetY[20,3] - lDetZ[20,3] <= 0
 lgkiU: - lGstk[47,3] + lDetY[20,3] - lDetZ[20,3] >= -1
 lgkiUa: - lgstka[47,3] + ldetYa[20,3] <= 0
 lgkiUa: ldetYa[20,3] - ldetZa[20,3] <= 0
 lgkiUa: - lgstka[47,3] + ldetYa[20,3] - ldetZa[20,3] >= -1
 lgkiU: - lGstk[47,4] + lDetY[20,4] <= 0
 lgkiU: lDetY[20,4] - lDetZ[20,4] <= 0
 lgkiU: - lGstk[47,4] + lDetY[20,4] - lDetZ[20,4] >= -1
 lgkiUa: - lgstka[47,4] + ldetYa[20,4] <= 0
 lgkiUa: ldetYa[20,4] - ldetZa[20,4] <= 0
 lgkiUa: - lgstka[47,4] + ldetYa[20,4] - ldetZa[20,4] >= -1
 lgkiU: - lGstk[47,5] + lDetY[20,5] <= 0
 lgkiU: lDetY[20,5] - lDetZ[20,5] <= 0
 lgkiU: - lGstk[47,5] + lDetY[20,5] - lDetZ[20,5] >= -1
 lgkiUa: - lgstka[47,5] + ldetYa[20,5] <= 0
 lgkiUa: ldetYa[20,5] - ldetZa[20,5] <= 0
 lgkiUa: - lgstka[47,5] + ldetYa[20,5] - ldetZa[20,5] >= -1
 lgkiU: - lGstk[47,6] + lDetY[20,6] <= 0
 lgkiU: lDetY[20,6] - lDetZ[20,6] <= 0
 lgkiU: - lGstk[47,6] + lDetY[20,6] - lDetZ[20,6] >= -1
 lgkiUa: - lgstka[47,6] + ldetYa[20,6] <= 0
 lgkiUa: ldetYa[20,6] - ldetZa[20,6] <= 0
 lgkiUa: - lgstka[47,6] + ldetYa[20,6] - ldetZa[20,6] >= -1
 lgkiU: - lGstk[47,7] + lDetY[20,7] <= 0
 lgkiU: lDetY[20,7] - lDetZ[20,7] <= 0
 lgkiU: - lGstk[47,7] + lDetY[20,7] - lDetZ[20,7] >= -1
 lgkiUa: - lgstka[47,7] + ldetYa[20,7] <= 0
 lgkiUa: ldetYa[20,7] - ldetZa[20,7] <= 0
 lgkiUa: - lgstka[47,7] + ldetYa[20,7] - ldetZa[20,7] >= -1
 lgkiL: lDetY[20,8] - lDetZ[20,8] = 0
 lgkiLa: ldetYa[20,8] - ldetZa[20,8] = 0
 lgkiL: lDetY[20,9] - lDetZ[20,9] = 0
 lgkiLa: ldetYa[20,9] - ldetZa[20,9] = 0
 lgkiL: lDetY[20,10] - lDetZ[20,10] = 0
 lgkiLa: ldetYa[20,10] - ldetZa[20,10] = 0
 lgkiL: lDetY[20,11] - lDetZ[20,11] = 0
 lgkiLa: ldetYa[20,11] - ldetZa[20,11] = 0
 lgkiL: lDetY[20,12] - lDetZ[20,12] = 0
 lgkiLa: ldetYa[20,12] - ldetZa[20,12] = 0
 lgkiL: lDetY[20,13] - lDetZ[20,13] = 0
 lgkiLa: ldetYa[20,13] - ldetZa[20,13] = 0
 lgkiL: lDetY[20,14] - lDetZ[20,14] = 0
 lgkiLa: ldetYa[20,14] - ldetZa[20,14] = 0
 lgkiL: lDetY[20,15] - lDetZ[20,15] = 0
 lgkiLa: ldetYa[20,15] - ldetZa[20,15] = 0
 lgkiU: - lGstk[48,0] + lDetY[21,0] <= 0
 lgkiU: lDetY[21,0] - lDetZ[21,0] <= 0
 lgkiU: - lGstk[48,0] + lDetY[21,0] - lDetZ[21,0] >= -1
 lgkiUa: - lgstka[48,0] + ldetYa[21,0] <= 0
 lgkiUa: ldetYa[21,0] - ldetZa[21,0] <= 0
 lgkiUa: - lgstka[48,0] + ldetYa[21,0] - ldetZa[21,0] >= -1
 lgkiU: - lGstk[48,1] + lDetY[21,1] <= 0
 lgkiU: lDetY[21,1] - lDetZ[21,1] <= 0
 lgkiU: - lGstk[48,1] + lDetY[21,1] - lDetZ[21,1] >= -1
 lgkiUa: - lgstka[48,1] + ldetYa[21,1] <= 0
 lgkiUa: ldetYa[21,1] - ldetZa[21,1] <= 0
 lgkiUa: - lgstka[48,1] + ldetYa[21,1] - ldetZa[21,1] >= -1
 lgkiU: - lGstk[48,2] + lDetY[21,2] <= 0
 lgkiU: lDetY[21,2] - lDetZ[21,2] <= 0
 lgkiU: - lGstk[48,2] + lDetY[21,2] - lDetZ[21,2] >= -1
 lgkiUa: - lgstka[48,2] + ldetYa[21,2] <= 0
 lgkiUa: ldetYa[21,2] - ldetZa[21,2] <= 0
 lgkiUa: - lgstka[48,2] + ldetYa[21,2] - ldetZa[21,2] >= -1
 lgkiU: - lGstk[48,3] + lDetY[21,3] <= 0
 lgkiU: lDetY[21,3] - lDetZ[21,3] <= 0
 lgkiU: - lGstk[48,3] + lDetY[21,3] - lDetZ[21,3] >= -1
 lgkiUa: - lgstka[48,3] + ldetYa[21,3] <= 0
 lgkiUa: ldetYa[21,3] - ldetZa[21,3] <= 0
 lgkiUa: - lgstka[48,3] + ldetYa[21,3] - ldetZa[21,3] >= -1
 lgkiU: - lGstk[48,4] + lDetY[21,4] <= 0
 lgkiU: lDetY[21,4] - lDetZ[21,4] <= 0
 lgkiU: - lGstk[48,4] + lDetY[21,4] - lDetZ[21,4] >= -1
 lgkiUa: - lgstka[48,4] + ldetYa[21,4] <= 0
 lgkiUa: ldetYa[21,4] - ldetZa[21,4] <= 0
 lgkiUa: - lgstka[48,4] + ldetYa[21,4] - ldetZa[21,4] >= -1
 lgkiU: - lGstk[48,5] + lDetY[21,5] <= 0
 lgkiU: lDetY[21,5] - lDetZ[21,5] <= 0
 lgkiU: - lGstk[48,5] + lDetY[21,5] - lDetZ[21,5] >= -1
 lgkiUa: - lgstka[48,5] + ldetYa[21,5] <= 0
 lgkiUa: ldetYa[21,5] - ldetZa[21,5] <= 0
 lgkiUa: - lgstka[48,5] + ldetYa[21,5] - ldetZa[21,5] >= -1
 lgkiU: - lGstk[48,6] + lDetY[21,6] <= 0
 lgkiU: lDetY[21,6] - lDetZ[21,6] <= 0
 lgkiU: - lGstk[48,6] + lDetY[21,6] - lDetZ[21,6] >= -1
 lgkiUa: - lgstka[48,6] + ldetYa[21,6] <= 0
 lgkiUa: ldetYa[21,6] - ldetZa[21,6] <= 0
 lgkiUa: - lgstka[48,6] + ldetYa[21,6] - ldetZa[21,6] >= -1
 lgkiU: - lGstk[48,7] + lDetY[21,7] <= 0
 lgkiU: lDetY[21,7] - lDetZ[21,7] <= 0
 lgkiU: - lGstk[48,7] + lDetY[21,7] - lDetZ[21,7] >= -1
 lgkiUa: - lgstka[48,7] + ldetYa[21,7] <= 0
 lgkiUa: ldetYa[21,7] - ldetZa[21,7] <= 0
 lgkiUa: - lgstka[48,7] + ldetYa[21,7] - ldetZa[21,7] >= -1
 lgkiL: lDetY[21,8] - lDetZ[21,8] = 0
 lgkiLa: ldetYa[21,8] - ldetZa[21,8] = 0
 lgkiL: lDetY[21,9] - lDetZ[21,9] = 0
 lgkiLa: ldetYa[21,9] - ldetZa[21,9] = 0
 lgkiL: lDetY[21,10] - lDetZ[21,10] = 0
 lgkiLa: ldetYa[21,10] - ldetZa[21,10] = 0
 lgkiL: lDetY[21,11] - lDetZ[21,11] = 0
 lgkiLa: ldetYa[21,11] - ldetZa[21,11] = 0
 lgkiL: lDetY[21,12] - lDetZ[21,12] = 0
 lgkiLa: ldetYa[21,12] - ldetZa[21,12] = 0
 lgkiL: lDetY[21,13] - lDetZ[21,13] = 0
 lgkiLa: ldetYa[21,13] - ldetZa[21,13] = 0
 lgkiL: lDetY[21,14] - lDetZ[21,14] = 0
 lgkiLa: ldetYa[21,14] - ldetZa[21,14] = 0
 lgkiL: lDetY[21,15] - lDetZ[21,15] = 0
 lgkiLa: ldetYa[21,15] - ldetZa[21,15] = 0
 lgkiU: - lGstk[49,0] + lDetY[22,0] <= 0
 lgkiU: lDetY[22,0] - lDetZ[22,0] <= 0
 lgkiU: - lGstk[49,0] + lDetY[22,0] - lDetZ[22,0] >= -1
 lgkiUa: - lgstka[49,0] + ldetYa[22,0] <= 0
 lgkiUa: ldetYa[22,0] - ldetZa[22,0] <= 0
 lgkiUa: - lgstka[49,0] + ldetYa[22,0] - ldetZa[22,0] >= -1
 lgkiU: - lGstk[49,1] + lDetY[22,1] <= 0
 lgkiU: lDetY[22,1] - lDetZ[22,1] <= 0
 lgkiU: - lGstk[49,1] + lDetY[22,1] - lDetZ[22,1] >= -1
 lgkiUa: - lgstka[49,1] + ldetYa[22,1] <= 0
 lgkiUa: ldetYa[22,1] - ldetZa[22,1] <= 0
 lgkiUa: - lgstka[49,1] + ldetYa[22,1] - ldetZa[22,1] >= -1
 lgkiU: - lGstk[49,2] + lDetY[22,2] <= 0
 lgkiU: lDetY[22,2] - lDetZ[22,2] <= 0
 lgkiU: - lGstk[49,2] + lDetY[22,2] - lDetZ[22,2] >= -1
 lgkiUa: - lgstka[49,2] + ldetYa[22,2] <= 0
 lgkiUa: ldetYa[22,2] - ldetZa[22,2] <= 0
 lgkiUa: - lgstka[49,2] + ldetYa[22,2] - ldetZa[22,2] >= -1
 lgkiU: - lGstk[49,3] + lDetY[22,3] <= 0
 lgkiU: lDetY[22,3] - lDetZ[22,3] <= 0
 lgkiU: - lGstk[49,3] + lDetY[22,3] - lDetZ[22,3] >= -1
 lgkiUa: - lgstka[49,3] + ldetYa[22,3] <= 0
 lgkiUa: ldetYa[22,3] - ldetZa[22,3] <= 0
 lgkiUa: - lgstka[49,3] + ldetYa[22,3] - ldetZa[22,3] >= -1
 lgkiU: - lGstk[49,4] + lDetY[22,4] <= 0
 lgkiU: lDetY[22,4] - lDetZ[22,4] <= 0
 lgkiU: - lGstk[49,4] + lDetY[22,4] - lDetZ[22,4] >= -1
 lgkiUa: - lgstka[49,4] + ldetYa[22,4] <= 0
 lgkiUa: ldetYa[22,4] - ldetZa[22,4] <= 0
 lgkiUa: - lgstka[49,4] + ldetYa[22,4] - ldetZa[22,4] >= -1
 lgkiU: - lGstk[49,5] + lDetY[22,5] <= 0
 lgkiU: lDetY[22,5] - lDetZ[22,5] <= 0
 lgkiU: - lGstk[49,5] + lDetY[22,5] - lDetZ[22,5] >= -1
 lgkiUa: - lgstka[49,5] + ldetYa[22,5] <= 0
 lgkiUa: ldetYa[22,5] - ldetZa[22,5] <= 0
 lgkiUa: - lgstka[49,5] + ldetYa[22,5] - ldetZa[22,5] >= -1
 lgkiU: - lGstk[49,6] + lDetY[22,6] <= 0
 lgkiU: lDetY[22,6] - lDetZ[22,6] <= 0
 lgkiU: - lGstk[49,6] + lDetY[22,6] - lDetZ[22,6] >= -1
 lgkiUa: - lgstka[49,6] + ldetYa[22,6] <= 0
 lgkiUa: ldetYa[22,6] - ldetZa[22,6] <= 0
 lgkiUa: - lgstka[49,6] + ldetYa[22,6] - ldetZa[22,6] >= -1
 lgkiU: - lGstk[49,7] + lDetY[22,7] <= 0
 lgkiU: lDetY[22,7] - lDetZ[22,7] <= 0
 lgkiU: - lGstk[49,7] + lDetY[22,7] - lDetZ[22,7] >= -1
 lgkiUa: - lgstka[49,7] + ldetYa[22,7] <= 0
 lgkiUa: ldetYa[22,7] - ldetZa[22,7] <= 0
 lgkiUa: - lgstka[49,7] + ldetYa[22,7] - ldetZa[22,7] >= -1
 lgkiL: lDetY[22,8] - lDetZ[22,8] = 0
 lgkiLa: ldetYa[22,8] - ldetZa[22,8] = 0
 lgkiL: lDetY[22,9] - lDetZ[22,9] = 0
 lgkiLa: ldetYa[22,9] - ldetZa[22,9] = 0
 lgkiL: lDetY[22,10] - lDetZ[22,10] = 0
 lgkiLa: ldetYa[22,10] - ldetZa[22,10] = 0
 lgkiL: lDetY[22,11] - lDetZ[22,11] = 0
 lgkiLa: ldetYa[22,11] - ldetZa[22,11] = 0
 lgkiL: lDetY[22,12] - lDetZ[22,12] = 0
 lgkiLa: ldetYa[22,12] - ldetZa[22,12] = 0
 lgkiL: lDetY[22,13] - lDetZ[22,13] = 0
 lgkiLa: ldetYa[22,13] - ldetZa[22,13] = 0
 lgkiL: lDetY[22,14] - lDetZ[22,14] = 0
 lgkiLa: ldetYa[22,14] - ldetZa[22,14] = 0
 lgkiL: lDetY[22,15] - lDetZ[22,15] = 0
 lgkiLa: ldetYa[22,15] - ldetZa[22,15] = 0
 lgkiU: - lGstk[50,0] + lDetY[23,0] <= 0
 lgkiU: lDetY[23,0] - lDetZ[23,0] <= 0
 lgkiU: - lGstk[50,0] + lDetY[23,0] - lDetZ[23,0] >= -1
 lgkiUa: - lgstka[50,0] + ldetYa[23,0] <= 0
 lgkiUa: ldetYa[23,0] - ldetZa[23,0] <= 0
 lgkiUa: - lgstka[50,0] + ldetYa[23,0] - ldetZa[23,0] >= -1
 lgkiU: - lGstk[50,1] + lDetY[23,1] <= 0
 lgkiU: lDetY[23,1] - lDetZ[23,1] <= 0
 lgkiU: - lGstk[50,1] + lDetY[23,1] - lDetZ[23,1] >= -1
 lgkiUa: - lgstka[50,1] + ldetYa[23,1] <= 0
 lgkiUa: ldetYa[23,1] - ldetZa[23,1] <= 0
 lgkiUa: - lgstka[50,1] + ldetYa[23,1] - ldetZa[23,1] >= -1
 lgkiU: - lGstk[50,2] + lDetY[23,2] <= 0
 lgkiU: lDetY[23,2] - lDetZ[23,2] <= 0
 lgkiU: - lGstk[50,2] + lDetY[23,2] - lDetZ[23,2] >= -1
 lgkiUa: - lgstka[50,2] + ldetYa[23,2] <= 0
 lgkiUa: ldetYa[23,2] - ldetZa[23,2] <= 0
 lgkiUa: - lgstka[50,2] + ldetYa[23,2] - ldetZa[23,2] >= -1
 lgkiU: - lGstk[50,3] + lDetY[23,3] <= 0
 lgkiU: lDetY[23,3] - lDetZ[23,3] <= 0
 lgkiU: - lGstk[50,3] + lDetY[23,3] - lDetZ[23,3] >= -1
 lgkiUa: - lgstka[50,3] + ldetYa[23,3] <= 0
 lgkiUa: ldetYa[23,3] - ldetZa[23,3] <= 0
 lgkiUa: - lgstka[50,3] + ldetYa[23,3] - ldetZa[23,3] >= -1
 lgkiU: - lGstk[50,4] + lDetY[23,4] <= 0
 lgkiU: lDetY[23,4] - lDetZ[23,4] <= 0
 lgkiU: - lGstk[50,4] + lDetY[23,4] - lDetZ[23,4] >= -1
 lgkiUa: - lgstka[50,4] + ldetYa[23,4] <= 0
 lgkiUa: ldetYa[23,4] - ldetZa[23,4] <= 0
 lgkiUa: - lgstka[50,4] + ldetYa[23,4] - ldetZa[23,4] >= -1
 lgkiU: - lGstk[50,5] + lDetY[23,5] <= 0
 lgkiU: lDetY[23,5] - lDetZ[23,5] <= 0
 lgkiU: - lGstk[50,5] + lDetY[23,5] - lDetZ[23,5] >= -1
 lgkiUa: - lgstka[50,5] + ldetYa[23,5] <= 0
 lgkiUa: ldetYa[23,5] - ldetZa[23,5] <= 0
 lgkiUa: - lgstka[50,5] + ldetYa[23,5] - ldetZa[23,5] >= -1
 lgkiU: - lGstk[50,6] + lDetY[23,6] <= 0
 lgkiU: lDetY[23,6] - lDetZ[23,6] <= 0
 lgkiU: - lGstk[50,6] + lDetY[23,6] - lDetZ[23,6] >= -1
 lgkiUa: - lgstka[50,6] + ldetYa[23,6] <= 0
 lgkiUa: ldetYa[23,6] - ldetZa[23,6] <= 0
 lgkiUa: - lgstka[50,6] + ldetYa[23,6] - ldetZa[23,6] >= -1
 lgkiU: - lGstk[50,7] + lDetY[23,7] <= 0
 lgkiU: lDetY[23,7] - lDetZ[23,7] <= 0
 lgkiU: - lGstk[50,7] + lDetY[23,7] - lDetZ[23,7] >= -1
 lgkiUa: - lgstka[50,7] + ldetYa[23,7] <= 0
 lgkiUa: ldetYa[23,7] - ldetZa[23,7] <= 0
 lgkiUa: - lgstka[50,7] + ldetYa[23,7] - ldetZa[23,7] >= -1
 lgkiL: lDetY[23,8] - lDetZ[23,8] = 0
 lgkiLa: ldetYa[23,8] - ldetZa[23,8] = 0
 lgkiL: lDetY[23,9] - lDetZ[23,9] = 0
 lgkiLa: ldetYa[23,9] - ldetZa[23,9] = 0
 lgkiL: lDetY[23,10] - lDetZ[23,10] = 0
 lgkiLa: ldetYa[23,10] - ldetZa[23,10] = 0
 lgkiL: lDetY[23,11] - lDetZ[23,11] = 0
 lgkiLa: ldetYa[23,11] - ldetZa[23,11] = 0
 lgkiL: lDetY[23,12] - lDetZ[23,12] = 0
 lgkiLa: ldetYa[23,12] - ldetZa[23,12] = 0
 lgkiL: lDetY[23,13] - lDetZ[23,13] = 0
 lgkiLa: ldetYa[23,13] - ldetZa[23,13] = 0
 lgkiL: lDetY[23,14] - lDetZ[23,14] = 0
 lgkiLa: ldetYa[23,14] - ldetZa[23,14] = 0
 lgkiL: lDetY[23,15] - lDetZ[23,15] = 0
 lgkiLa: ldetYa[23,15] - ldetZa[23,15] = 0
 ldetSR: lDetZ[19,0] - lDetW[19,0] = 0
 ldetSRa: ldetZa[19,0] - ldetWa[19,0] = 0
 ldetSR: lDetZ[19,1] - lDetW[19,1] = 0
 ldetSRa: ldetZa[19,1] - ldetWa[19,1] = 0
 ldetSR: lDetZ[19,2] - lDetW[19,2] = 0
 ldetSRa: ldetZa[19,2] - ldetWa[19,2] = 0
 ldetSR: lDetZ[19,3] - lDetW[19,3] = 0
 ldetSRa: ldetZa[19,3] - ldetWa[19,3] = 0
 ldetSR: lDetZ[19,4] - lDetW[19,5] = 0
 ldetSRa: ldetZa[19,4] - ldetWa[19,5] = 0
 ldetSR: lDetZ[19,5] - lDetW[19,6] = 0
 ldetSRa: ldetZa[19,5] - ldetWa[19,6] = 0
 ldetSR: lDetZ[19,6] - lDetW[19,7] = 0
 ldetSRa: ldetZa[19,6] - ldetWa[19,7] = 0
 ldetSR: lDetZ[19,7] - lDetW[19,4] = 0
 ldetSRa: ldetZa[19,7] - ldetWa[19,4] = 0
 ldetSR: lDetZ[19,8] - lDetW[19,10] = 0
 ldetSRa: ldetZa[19,8] - ldetWa[19,10] = 0
 ldetSR: lDetZ[19,9] - lDetW[19,11] = 0
 ldetSRa: ldetZa[19,9] - ldetWa[19,11] = 0
 ldetSR: lDetZ[19,10] - lDetW[19,8] = 0
 ldetSRa: ldetZa[19,10] - ldetWa[19,8] = 0
 ldetSR: lDetZ[19,11] - lDetW[19,9] = 0
 ldetSRa: ldetZa[19,11] - ldetWa[19,9] = 0
 ldetSR: lDetZ[19,12] - lDetW[19,15] = 0
 ldetSRa: ldetZa[19,12] - ldetWa[19,15] = 0
 ldetSR: lDetZ[19,13] - lDetW[19,12] = 0
 ldetSRa: ldetZa[19,13] - ldetWa[19,12] = 0
 ldetSR: lDetZ[19,14] - lDetW[19,13] = 0
 ldetSRa: ldetZa[19,14] - ldetWa[19,13] = 0
 ldetSR: lDetZ[19,15] - lDetW[19,14] = 0
 ldetSRa: ldetZa[19,15] - ldetWa[19,14] = 0
 ldetSR: lDetZ[20,0] - lDetW[20,0] = 0
 ldetSRa: ldetZa[20,0] - ldetWa[20,0] = 0
 ldetSR: lDetZ[20,1] - lDetW[20,1] = 0
 ldetSRa: ldetZa[20,1] - ldetWa[20,1] = 0
 ldetSR: lDetZ[20,2] - lDetW[20,2] = 0
 ldetSRa: ldetZa[20,2] - ldetWa[20,2] = 0
 ldetSR: lDetZ[20,3] - lDetW[20,3] = 0
 ldetSRa: ldetZa[20,3] - ldetWa[20,3] = 0
 ldetSR: lDetZ[20,4] - lDetW[20,5] = 0
 ldetSRa: ldetZa[20,4] - ldetWa[20,5] = 0
 ldetSR: lDetZ[20,5] - lDetW[20,6] = 0
 ldetSRa: ldetZa[20,5] - ldetWa[20,6] = 0
 ldetSR: lDetZ[20,6] - lDetW[20,7] = 0
 ldetSRa: ldetZa[20,6] - ldetWa[20,7] = 0
 ldetSR: lDetZ[20,7] - lDetW[20,4] = 0
 ldetSRa: ldetZa[20,7] - ldetWa[20,4] = 0
 ldetSR: lDetZ[20,8] - lDetW[20,10] = 0
 ldetSRa: ldetZa[20,8] - ldetWa[20,10] = 0
 ldetSR: lDetZ[20,9] - lDetW[20,11] = 0
 ldetSRa: ldetZa[20,9] - ldetWa[20,11] = 0
 ldetSR: lDetZ[20,10] - lDetW[20,8] = 0
 ldetSRa: ldetZa[20,10] - ldetWa[20,8] = 0
 ldetSR: lDetZ[20,11] - lDetW[20,9] = 0
 ldetSRa: ldetZa[20,11] - ldetWa[20,9] = 0
 ldetSR: lDetZ[20,12] - lDetW[20,15] = 0
 ldetSRa: ldetZa[20,12] - ldetWa[20,15] = 0
 ldetSR: lDetZ[20,13] - lDetW[20,12] = 0
 ldetSRa: ldetZa[20,13] - ldetWa[20,12] = 0
 ldetSR: lDetZ[20,14] - lDetW[20,13] = 0
 ldetSRa: ldetZa[20,14] - ldetWa[20,13] = 0
 ldetSR: lDetZ[20,15] - lDetW[20,14] = 0
 ldetSRa: ldetZa[20,15] - ldetWa[20,14] = 0
 ldetSR: lDetZ[21,0] - lDetW[21,0] = 0
 ldetSRa: ldetZa[21,0] - ldetWa[21,0] = 0
 ldetSR: lDetZ[21,1] - lDetW[21,1] = 0
 ldetSRa: ldetZa[21,1] - ldetWa[21,1] = 0
 ldetSR: lDetZ[21,2] - lDetW[21,2] = 0
 ldetSRa: ldetZa[21,2] - ldetWa[21,2] = 0
 ldetSR: lDetZ[21,3] - lDetW[21,3] = 0
 ldetSRa: ldetZa[21,3] - ldetWa[21,3] = 0
 ldetSR: lDetZ[21,4] - lDetW[21,5] = 0
 ldetSRa: ldetZa[21,4] - ldetWa[21,5] = 0
 ldetSR: lDetZ[21,5] - lDetW[21,6] = 0
 ldetSRa: ldetZa[21,5] - ldetWa[21,6] = 0
 ldetSR: lDetZ[21,6] - lDetW[21,7] = 0
 ldetSRa: ldetZa[21,6] - ldetWa[21,7] = 0
 ldetSR: lDetZ[21,7] - lDetW[21,4] = 0
 ldetSRa: ldetZa[21,7] - ldetWa[21,4] = 0
 ldetSR: lDetZ[21,8] - lDetW[21,10] = 0
 ldetSRa: ldetZa[21,8] - ldetWa[21,10] = 0
 ldetSR: lDetZ[21,9] - lDetW[21,11] = 0
 ldetSRa: ldetZa[21,9] - ldetWa[21,11] = 0
 ldetSR: lDetZ[21,10] - lDetW[21,8] = 0
 ldetSRa: ldetZa[21,10] - ldetWa[21,8] = 0
 ldetSR: lDetZ[21,11] - lDetW[21,9] = 0
 ldetSRa: ldetZa[21,11] - ldetWa[21,9] = 0
 ldetSR: lDetZ[21,12] - lDetW[21,15] = 0
 ldetSRa: ldetZa[21,12] - ldetWa[21,15] = 0
 ldetSR: lDetZ[21,13] - lDetW[21,12] = 0
 ldetSRa: ldetZa[21,13] - ldetWa[21,12] = 0
 ldetSR: lDetZ[21,14] - lDetW[21,13] = 0
 ldetSRa: ldetZa[21,14] - ldetWa[21,13] = 0
 ldetSR: lDetZ[21,15] - lDetW[21,14] = 0
 ldetSRa: ldetZa[21,15] - ldetWa[21,14] = 0
 ldetSR: lDetZ[22,0] - lDetW[22,0] = 0
 ldetSRa: ldetZa[22,0] - ldetWa[22,0] = 0
 ldetSR: lDetZ[22,1] - lDetW[22,1] = 0
 ldetSRa: ldetZa[22,1] - ldetWa[22,1] = 0
 ldetSR: lDetZ[22,2] - lDetW[22,2] = 0
 ldetSRa: ldetZa[22,2] - ldetWa[22,2] = 0
 ldetSR: lDetZ[22,3] - lDetW[22,3] = 0
 ldetSRa: ldetZa[22,3] - ldetWa[22,3] = 0
 ldetSR: lDetZ[22,4] - lDetW[22,5] = 0
 ldetSRa: ldetZa[22,4] - ldetWa[22,5] = 0
 ldetSR: lDetZ[22,5] - lDetW[22,6] = 0
 ldetSRa: ldetZa[22,5] - ldetWa[22,6] = 0
 ldetSR: lDetZ[22,6] - lDetW[22,7] = 0
 ldetSRa: ldetZa[22,6] - ldetWa[22,7] = 0
 ldetSR: lDetZ[22,7] - lDetW[22,4] = 0
 ldetSRa: ldetZa[22,7] - ldetWa[22,4] = 0
 ldetSR: lDetZ[22,8] - lDetW[22,10] = 0
 ldetSRa: ldetZa[22,8] - ldetWa[22,10] = 0
 ldetSR: lDetZ[22,9] - lDetW[22,11] = 0
 ldetSRa: ldetZa[22,9] - ldetWa[22,11] = 0
 ldetSR: lDetZ[22,10] - lDetW[22,8] = 0
 ldetSRa: ldetZa[22,10] - ldetWa[22,8] = 0
 ldetSR: lDetZ[22,11] - lDetW[22,9] = 0
 ldetSRa: ldetZa[22,11] - ldetWa[22,9] = 0
 ldetSR: lDetZ[22,12] - lDetW[22,15] = 0
 ldetSRa: ldetZa[22,12] - ldetWa[22,15] = 0
 ldetSR: lDetZ[22,13] - lDetW[22,12] = 0
 ldetSRa: ldetZa[22,13] - ldetWa[22,12] = 0
 ldetSR: lDetZ[22,14] - lDetW[22,13] = 0
 ldetSRa: ldetZa[22,14] - ldetWa[22,13] = 0
 ldetSR: lDetZ[22,15] - lDetW[22,14] = 0
 ldetSRa: ldetZa[22,15] - ldetWa[22,14] = 0
 ldetMC0: - lDetX[20,4] + lDetW[19,0] = 0
 ldetMC1: lDetW[19,4] - lDetW[19,8] <= 0
 ldetMC1: - lDetX[20,8] + lDetW[19,4] <= 0
 ldetMC1: - lDetX[20,8] + lDetW[19,4] - lDetW[19,8] >= -1
 ldetMC2: - lDetX[20,4] + lDetW[19,8] <= 0
 ldetMC2: - lDetX[20,12] + lDetW[19,8] <= 0
 ldetMC2: - lDetX[20,4] - lDetX[20,12] + lDetW[19,8] >= -1
 ldetMC3: - lDetX[20,0] + lDetW[19,12] <= 0
 ldetMC3: - lDetX[20,12] + lDetW[19,12] <= 0
 ldetMC3: - lDetX[20,0] - lDetX[20,12] + lDetW[19,12] >= -1
 ldetMC0a: - ldetXa[20,4] + ldetWa[19,0] = 0
 ldetMC1a: ldetWa[19,4] - ldetWa[19,8] <= 0
 ldetMC1a: - ldetXa[20,8] + ldetWa[19,4] <= 0
 ldetMC1a: - ldetXa[20,8] + ldetWa[19,4] - ldetWa[19,8] >= -1
 ldetMC2a: - ldetXa[20,4] + ldetWa[19,8] <= 0
 ldetMC2a: - ldetXa[20,12] + ldetWa[19,8] <= 0
 ldetMC2a: - ldetXa[20,4] - ldetXa[20,12] + ldetWa[19,8] >= -1
 ldetMC3a: - ldetXa[20,0] + ldetWa[19,12] <= 0
 ldetMC3a: - ldetXa[20,12] + ldetWa[19,12] <= 0
 ldetMC3a: - ldetXa[20,0] - ldetXa[20,12] + ldetWa[19,12] >= -1
 ldetMC0: - lDetX[20,5] + lDetW[19,1] = 0
 ldetMC1: lDetW[19,5] - lDetW[19,9] <= 0
 ldetMC1: - lDetX[20,9] + lDetW[19,5] <= 0
 ldetMC1: - lDetX[20,9] + lDetW[19,5] - lDetW[19,9] >= -1
 ldetMC2: - lDetX[20,5] + lDetW[19,9] <= 0
 ldetMC2: - lDetX[20,13] + lDetW[19,9] <= 0
 ldetMC2: - lDetX[20,5] - lDetX[20,13] + lDetW[19,9] >= -1
 ldetMC3: - lDetX[20,1] + lDetW[19,13] <= 0
 ldetMC3: - lDetX[20,13] + lDetW[19,13] <= 0
 ldetMC3: - lDetX[20,1] - lDetX[20,13] + lDetW[19,13] >= -1
 ldetMC0a: - ldetXa[20,5] + ldetWa[19,1] = 0
 ldetMC1a: ldetWa[19,5] - ldetWa[19,9] <= 0
 ldetMC1a: - ldetXa[20,9] + ldetWa[19,5] <= 0
 ldetMC1a: - ldetXa[20,9] + ldetWa[19,5] - ldetWa[19,9] >= -1
 ldetMC2a: - ldetXa[20,5] + ldetWa[19,9] <= 0
 ldetMC2a: - ldetXa[20,13] + ldetWa[19,9] <= 0
 ldetMC2a: - ldetXa[20,5] - ldetXa[20,13] + ldetWa[19,9] >= -1
 ldetMC3a: - ldetXa[20,1] + ldetWa[19,13] <= 0
 ldetMC3a: - ldetXa[20,13] + ldetWa[19,13] <= 0
 ldetMC3a: - ldetXa[20,1] - ldetXa[20,13] + ldetWa[19,13] >= -1
 ldetMC0: - lDetX[20,6] + lDetW[19,2] = 0
 ldetMC1: lDetW[19,6] - lDetW[19,10] <= 0
 ldetMC1: - lDetX[20,10] + lDetW[19,6] <= 0
 ldetMC1: - lDetX[20,10] + lDetW[19,6] - lDetW[19,10] >= -1
 ldetMC2: - lDetX[20,6] + lDetW[19,10] <= 0
 ldetMC2: - lDetX[20,14] + lDetW[19,10] <= 0
 ldetMC2: - lDetX[20,6] - lDetX[20,14] + lDetW[19,10] >= -1
 ldetMC3: - lDetX[20,2] + lDetW[19,14] <= 0
 ldetMC3: - lDetX[20,14] + lDetW[19,14] <= 0
 ldetMC3: - lDetX[20,2] - lDetX[20,14] + lDetW[19,14] >= -1
 ldetMC0a: - ldetXa[20,6] + ldetWa[19,2] = 0
 ldetMC1a: ldetWa[19,6] - ldetWa[19,10] <= 0
 ldetMC1a: - ldetXa[20,10] + ldetWa[19,6] <= 0
 ldetMC1a: - ldetXa[20,10] + ldetWa[19,6] - ldetWa[19,10] >= -1
 ldetMC2a: - ldetXa[20,6] + ldetWa[19,10] <= 0
 ldetMC2a: - ldetXa[20,14] + ldetWa[19,10] <= 0
 ldetMC2a: - ldetXa[20,6] - ldetXa[20,14] + ldetWa[19,10] >= -1
 ldetMC3a: - ldetXa[20,2] + ldetWa[19,14] <= 0
 ldetMC3a: - ldetXa[20,14] + ldetWa[19,14] <= 0
 ldetMC3a: - ldetXa[20,2] - ldetXa[20,14] + ldetWa[19,14] >= -1
 ldetMC0: - lDetX[20,7] + lDetW[19,3] = 0
 ldetMC1: lDetW[19,7] - lDetW[19,11] <= 0
 ldetMC1: - lDetX[20,11] + lDetW[19,7] <= 0
 ldetMC1: - lDetX[20,11] + lDetW[19,7] - lDetW[19,11] >= -1
 ldetMC2: - lDetX[20,7] + lDetW[19,11] <= 0
 ldetMC2: - lDetX[20,15] + lDetW[19,11] <= 0
 ldetMC2: - lDetX[20,7] - lDetX[20,15] + lDetW[19,11] >= -1
 ldetMC3: - lDetX[20,3] + lDetW[19,15] <= 0
 ldetMC3: - lDetX[20,15] + lDetW[19,15] <= 0
 ldetMC3: - lDetX[20,3] - lDetX[20,15] + lDetW[19,15] >= -1
 ldetMC0a: - ldetXa[20,7] + ldetWa[19,3] = 0
 ldetMC1a: ldetWa[19,7] - ldetWa[19,11] <= 0
 ldetMC1a: - ldetXa[20,11] + ldetWa[19,7] <= 0
 ldetMC1a: - ldetXa[20,11] + ldetWa[19,7] - ldetWa[19,11] >= -1
 ldetMC2a: - ldetXa[20,7] + ldetWa[19,11] <= 0
 ldetMC2a: - ldetXa[20,15] + ldetWa[19,11] <= 0
 ldetMC2a: - ldetXa[20,7] - ldetXa[20,15] + ldetWa[19,11] >= -1
 ldetMC3a: - ldetXa[20,3] + ldetWa[19,15] <= 0
 ldetMC3a: - ldetXa[20,15] + ldetWa[19,15] <= 0
 ldetMC3a: - ldetXa[20,3] - ldetXa[20,15] + ldetWa[19,15] >= -1
 ldetMC0: - lDetX[21,4] + lDetW[20,0] = 0
 ldetMC1: lDetW[20,4] - lDetW[20,8] <= 0
 ldetMC1: - lDetX[21,8] + lDetW[20,4] <= 0
 ldetMC1: - lDetX[21,8] + lDetW[20,4] - lDetW[20,8] >= -1
 ldetMC2: - lDetX[21,4] + lDetW[20,8] <= 0
 ldetMC2: - lDetX[21,12] + lDetW[20,8] <= 0
 ldetMC2: - lDetX[21,4] - lDetX[21,12] + lDetW[20,8] >= -1
 ldetMC3: - lDetX[21,0] + lDetW[20,12] <= 0
 ldetMC3: - lDetX[21,12] + lDetW[20,12] <= 0
 ldetMC3: - lDetX[21,0] - lDetX[21,12] + lDetW[20,12] >= -1
 ldetMC0a: - ldetXa[21,4] + ldetWa[20,0] = 0
 ldetMC1a: ldetWa[20,4] - ldetWa[20,8] <= 0
 ldetMC1a: - ldetXa[21,8] + ldetWa[20,4] <= 0
 ldetMC1a: - ldetXa[21,8] + ldetWa[20,4] - ldetWa[20,8] >= -1
 ldetMC2a: - ldetXa[21,4] + ldetWa[20,8] <= 0
 ldetMC2a: - ldetXa[21,12] + ldetWa[20,8] <= 0
 ldetMC2a: - ldetXa[21,4] - ldetXa[21,12] + ldetWa[20,8] >= -1
 ldetMC3a: - ldetXa[21,0] + ldetWa[20,12] <= 0
 ldetMC3a: - ldetXa[21,12] + ldetWa[20,12] <= 0
 ldetMC3a: - ldetXa[21,0] - ldetXa[21,12] + ldetWa[20,12] >= -1
 ldetMC0: - lDetX[21,5] + lDetW[20,1] = 0
 ldetMC1: lDetW[20,5] - lDetW[20,9] <= 0
 ldetMC1: - lDetX[21,9] + lDetW[20,5] <= 0
 ldetMC1: - lDetX[21,9] + lDetW[20,5] - lDetW[20,9] >= -1
 ldetMC2: - lDetX[21,5] + lDetW[20,9] <= 0
 ldetMC2: - lDetX[21,13] + lDetW[20,9] <= 0
 ldetMC2: - lDetX[21,5] - lDetX[21,13] + lDetW[20,9] >= -1
 ldetMC3: - lDetX[21,1] + lDetW[20,13] <= 0
 ldetMC3: - lDetX[21,13] + lDetW[20,13] <= 0
 ldetMC3: - lDetX[21,1] - lDetX[21,13] + lDetW[20,13] >= -1
 ldetMC0a: - ldetXa[21,5] + ldetWa[20,1] = 0
 ldetMC1a: ldetWa[20,5] - ldetWa[20,9] <= 0
 ldetMC1a: - ldetXa[21,9] + ldetWa[20,5] <= 0
 ldetMC1a: - ldetXa[21,9] + ldetWa[20,5] - ldetWa[20,9] >= -1
 ldetMC2a: - ldetXa[21,5] + ldetWa[20,9] <= 0
 ldetMC2a: - ldetXa[21,13] + ldetWa[20,9] <= 0
 ldetMC2a: - ldetXa[21,5] - ldetXa[21,13] + ldetWa[20,9] >= -1
 ldetMC3a: - ldetXa[21,1] + ldetWa[20,13] <= 0
 ldetMC3a: - ldetXa[21,13] + ldetWa[20,13] <= 0
 ldetMC3a: - ldetXa[21,1] - ldetXa[21,13] + ldetWa[20,13] >= -1
 ldetMC0: - lDetX[21,6] + lDetW[20,2] = 0
 ldetMC1: lDetW[20,6] - lDetW[20,10] <= 0
 ldetMC1: - lDetX[21,10] + lDetW[20,6] <= 0
 ldetMC1: - lDetX[21,10] + lDetW[20,6] - lDetW[20,10] >= -1
 ldetMC2: - lDetX[21,6] + lDetW[20,10] <= 0
 ldetMC2: - lDetX[21,14] + lDetW[20,10] <= 0
 ldetMC2: - lDetX[21,6] - lDetX[21,14] + lDetW[20,10] >= -1
 ldetMC3: - lDetX[21,2] + lDetW[20,14] <= 0
 ldetMC3: - lDetX[21,14] + lDetW[20,14] <= 0
 ldetMC3: - lDetX[21,2] - lDetX[21,14] + lDetW[20,14] >= -1
 ldetMC0a: - ldetXa[21,6] + ldetWa[20,2] = 0
 ldetMC1a: ldetWa[20,6] - ldetWa[20,10] <= 0
 ldetMC1a: - ldetXa[21,10] + ldetWa[20,6] <= 0
 ldetMC1a: - ldetXa[21,10] + ldetWa[20,6] - ldetWa[20,10] >= -1
 ldetMC2a: - ldetXa[21,6] + ldetWa[20,10] <= 0
 ldetMC2a: - ldetXa[21,14] + ldetWa[20,10] <= 0
 ldetMC2a: - ldetXa[21,6] - ldetXa[21,14] + ldetWa[20,10] >= -1
 ldetMC3a: - ldetXa[21,2] + ldetWa[20,14] <= 0
 ldetMC3a: - ldetXa[21,14] + ldetWa[20,14] <= 0
 ldetMC3a: - ldetXa[21,2] - ldetXa[21,14] + ldetWa[20,14] >= -1
 ldetMC0: - lDetX[21,7] + lDetW[20,3] = 0
 ldetMC1: lDetW[20,7] - lDetW[20,11] <= 0
 ldetMC1: - lDetX[21,11] + lDetW[20,7] <= 0
 ldetMC1: - lDetX[21,11] + lDetW[20,7] - lDetW[20,11] >= -1
 ldetMC2: - lDetX[21,7] + lDetW[20,11] <= 0
 ldetMC2: - lDetX[21,15] + lDetW[20,11] <= 0
 ldetMC2: - lDetX[21,7] - lDetX[21,15] + lDetW[20,11] >= -1
 ldetMC3: - lDetX[21,3] + lDetW[20,15] <= 0
 ldetMC3: - lDetX[21,15] + lDetW[20,15] <= 0
 ldetMC3: - lDetX[21,3] - lDetX[21,15] + lDetW[20,15] >= -1
 ldetMC0a: - ldetXa[21,7] + ldetWa[20,3] = 0
 ldetMC1a: ldetWa[20,7] - ldetWa[20,11] <= 0
 ldetMC1a: - ldetXa[21,11] + ldetWa[20,7] <= 0
 ldetMC1a: - ldetXa[21,11] + ldetWa[20,7] - ldetWa[20,11] >= -1
 ldetMC2a: - ldetXa[21,7] + ldetWa[20,11] <= 0
 ldetMC2a: - ldetXa[21,15] + ldetWa[20,11] <= 0
 ldetMC2a: - ldetXa[21,7] - ldetXa[21,15] + ldetWa[20,11] >= -1
 ldetMC3a: - ldetXa[21,3] + ldetWa[20,15] <= 0
 ldetMC3a: - ldetXa[21,15] + ldetWa[20,15] <= 0
 ldetMC3a: - ldetXa[21,3] - ldetXa[21,15] + ldetWa[20,15] >= -1
 ldetMC0: - lDetX[22,4] + lDetW[21,0] = 0
 ldetMC1: lDetW[21,4] - lDetW[21,8] <= 0
 ldetMC1: - lDetX[22,8] + lDetW[21,4] <= 0
 ldetMC1: - lDetX[22,8] + lDetW[21,4] - lDetW[21,8] >= -1
 ldetMC2: - lDetX[22,4] + lDetW[21,8] <= 0
 ldetMC2: - lDetX[22,12] + lDetW[21,8] <= 0
 ldetMC2: - lDetX[22,4] - lDetX[22,12] + lDetW[21,8] >= -1
 ldetMC3: - lDetX[22,0] + lDetW[21,12] <= 0
 ldetMC3: - lDetX[22,12] + lDetW[21,12] <= 0
 ldetMC3: - lDetX[22,0] - lDetX[22,12] + lDetW[21,12] >= -1
 ldetMC0a: - ldetXa[22,4] + ldetWa[21,0] = 0
 ldetMC1a: ldetWa[21,4] - ldetWa[21,8] <= 0
 ldetMC1a: - ldetXa[22,8] + ldetWa[21,4] <= 0
 ldetMC1a: - ldetXa[22,8] + ldetWa[21,4] - ldetWa[21,8] >= -1
 ldetMC2a: - ldetXa[22,4] + ldetWa[21,8] <= 0
 ldetMC2a: - ldetXa[22,12] + ldetWa[21,8] <= 0
 ldetMC2a: - ldetXa[22,4] - ldetXa[22,12] + ldetWa[21,8] >= -1
 ldetMC3a: - ldetXa[22,0] + ldetWa[21,12] <= 0
 ldetMC3a: - ldetXa[22,12] + ldetWa[21,12] <= 0
 ldetMC3a: - ldetXa[22,0] - ldetXa[22,12] + ldetWa[21,12] >= -1
 ldetMC0: - lDetX[22,5] + lDetW[21,1] = 0
 ldetMC1: lDetW[21,5] - lDetW[21,9] <= 0
 ldetMC1: - lDetX[22,9] + lDetW[21,5] <= 0
 ldetMC1: - lDetX[22,9] + lDetW[21,5] - lDetW[21,9] >= -1
 ldetMC2: - lDetX[22,5] + lDetW[21,9] <= 0
 ldetMC2: - lDetX[22,13] + lDetW[21,9] <= 0
 ldetMC2: - lDetX[22,5] - lDetX[22,13] + lDetW[21,9] >= -1
 ldetMC3: - lDetX[22,1] + lDetW[21,13] <= 0
 ldetMC3: - lDetX[22,13] + lDetW[21,13] <= 0
 ldetMC3: - lDetX[22,1] - lDetX[22,13] + lDetW[21,13] >= -1
 ldetMC0a: - ldetXa[22,5] + ldetWa[21,1] = 0
 ldetMC1a: ldetWa[21,5] - ldetWa[21,9] <= 0
 ldetMC1a: - ldetXa[22,9] + ldetWa[21,5] <= 0
 ldetMC1a: - ldetXa[22,9] + ldetWa[21,5] - ldetWa[21,9] >= -1
 ldetMC2a: - ldetXa[22,5] + ldetWa[21,9] <= 0
 ldetMC2a: - ldetXa[22,13] + ldetWa[21,9] <= 0
 ldetMC2a: - ldetXa[22,5] - ldetXa[22,13] + ldetWa[21,9] >= -1
 ldetMC3a: - ldetXa[22,1] + ldetWa[21,13] <= 0
 ldetMC3a: - ldetXa[22,13] + ldetWa[21,13] <= 0
 ldetMC3a: - ldetXa[22,1] - ldetXa[22,13] + ldetWa[21,13] >= -1
 ldetMC0: - lDetX[22,6] + lDetW[21,2] = 0
 ldetMC1: lDetW[21,6] - lDetW[21,10] <= 0
 ldetMC1: - lDetX[22,10] + lDetW[21,6] <= 0
 ldetMC1: - lDetX[22,10] + lDetW[21,6] - lDetW[21,10] >= -1
 ldetMC2: - lDetX[22,6] + lDetW[21,10] <= 0
 ldetMC2: - lDetX[22,14] + lDetW[21,10] <= 0
 ldetMC2: - lDetX[22,6] - lDetX[22,14] + lDetW[21,10] >= -1
 ldetMC3: - lDetX[22,2] + lDetW[21,14] <= 0
 ldetMC3: - lDetX[22,14] + lDetW[21,14] <= 0
 ldetMC3: - lDetX[22,2] - lDetX[22,14] + lDetW[21,14] >= -1
 ldetMC0a: - ldetXa[22,6] + ldetWa[21,2] = 0
 ldetMC1a: ldetWa[21,6] - ldetWa[21,10] <= 0
 ldetMC1a: - ldetXa[22,10] + ldetWa[21,6] <= 0
 ldetMC1a: - ldetXa[22,10] + ldetWa[21,6] - ldetWa[21,10] >= -1
 ldetMC2a: - ldetXa[22,6] + ldetWa[21,10] <= 0
 ldetMC2a: - ldetXa[22,14] + ldetWa[21,10] <= 0
 ldetMC2a: - ldetXa[22,6] - ldetXa[22,14] + ldetWa[21,10] >= -1
 ldetMC3a: - ldetXa[22,2] + ldetWa[21,14] <= 0
 ldetMC3a: - ldetXa[22,14] + ldetWa[21,14] <= 0
 ldetMC3a: - ldetXa[22,2] - ldetXa[22,14] + ldetWa[21,14] >= -1
 ldetMC0: - lDetX[22,7] + lDetW[21,3] = 0
 ldetMC1: lDetW[21,7] - lDetW[21,11] <= 0
 ldetMC1: - lDetX[22,11] + lDetW[21,7] <= 0
 ldetMC1: - lDetX[22,11] + lDetW[21,7] - lDetW[21,11] >= -1
 ldetMC2: - lDetX[22,7] + lDetW[21,11] <= 0
 ldetMC2: - lDetX[22,15] + lDetW[21,11] <= 0
 ldetMC2: - lDetX[22,7] - lDetX[22,15] + lDetW[21,11] >= -1
 ldetMC3: - lDetX[22,3] + lDetW[21,15] <= 0
 ldetMC3: - lDetX[22,15] + lDetW[21,15] <= 0
 ldetMC3: - lDetX[22,3] - lDetX[22,15] + lDetW[21,15] >= -1
 ldetMC0a: - ldetXa[22,7] + ldetWa[21,3] = 0
 ldetMC1a: ldetWa[21,7] - ldetWa[21,11] <= 0
 ldetMC1a: - ldetXa[22,11] + ldetWa[21,7] <= 0
 ldetMC1a: - ldetXa[22,11] + ldetWa[21,7] - ldetWa[21,11] >= -1
 ldetMC2a: - ldetXa[22,7] + ldetWa[21,11] <= 0
 ldetMC2a: - ldetXa[22,15] + ldetWa[21,11] <= 0
 ldetMC2a: - ldetXa[22,7] - ldetXa[22,15] + ldetWa[21,11] >= -1
 ldetMC3a: - ldetXa[22,3] + ldetWa[21,15] <= 0
 ldetMC3a: - ldetXa[22,15] + ldetWa[21,15] <= 0
 ldetMC3a: - ldetXa[22,3] - ldetXa[22,15] + ldetWa[21,15] >= -1
 ldetMC0: - lDetX[23,4] + lDetW[22,0] = 0
 ldetMC1: lDetW[22,4] - lDetW[22,8] <= 0
 ldetMC1: - lDetX[23,8] + lDetW[22,4] <= 0
 ldetMC1: - lDetX[23,8] + lDetW[22,4] - lDetW[22,8] >= -1
 ldetMC2: - lDetX[23,4] + lDetW[22,8] <= 0
 ldetMC2: - lDetX[23,12] + lDetW[22,8] <= 0
 ldetMC2: - lDetX[23,4] - lDetX[23,12] + lDetW[22,8] >= -1
 ldetMC3: - lDetX[23,0] + lDetW[22,12] <= 0
 ldetMC3: - lDetX[23,12] + lDetW[22,12] <= 0
 ldetMC3: - lDetX[23,0] - lDetX[23,12] + lDetW[22,12] >= -1
 ldetMC0a: - ldetXa[23,4] + ldetWa[22,0] = 0
 ldetMC1a: ldetWa[22,4] - ldetWa[22,8] <= 0
 ldetMC1a: - ldetXa[23,8] + ldetWa[22,4] <= 0
 ldetMC1a: - ldetXa[23,8] + ldetWa[22,4] - ldetWa[22,8] >= -1
 ldetMC2a: - ldetXa[23,4] + ldetWa[22,8] <= 0
 ldetMC2a: - ldetXa[23,12] + ldetWa[22,8] <= 0
 ldetMC2a: - ldetXa[23,4] - ldetXa[23,12] + ldetWa[22,8] >= -1
 ldetMC3a: - ldetXa[23,0] + ldetWa[22,12] <= 0
 ldetMC3a: - ldetXa[23,12] + ldetWa[22,12] <= 0
 ldetMC3a: - ldetXa[23,0] - ldetXa[23,12] + ldetWa[22,12] >= -1
 ldetMC0: - lDetX[23,5] + lDetW[22,1] = 0
 ldetMC1: lDetW[22,5] - lDetW[22,9] <= 0
 ldetMC1: - lDetX[23,9] + lDetW[22,5] <= 0
 ldetMC1: - lDetX[23,9] + lDetW[22,5] - lDetW[22,9] >= -1
 ldetMC2: - lDetX[23,5] + lDetW[22,9] <= 0
 ldetMC2: - lDetX[23,13] + lDetW[22,9] <= 0
 ldetMC2: - lDetX[23,5] - lDetX[23,13] + lDetW[22,9] >= -1
 ldetMC3: - lDetX[23,1] + lDetW[22,13] <= 0
 ldetMC3: - lDetX[23,13] + lDetW[22,13] <= 0
 ldetMC3: - lDetX[23,1] - lDetX[23,13] + lDetW[22,13] >= -1
 ldetMC0a: - ldetXa[23,5] + ldetWa[22,1] = 0
 ldetMC1a: ldetWa[22,5] - ldetWa[22,9] <= 0
 ldetMC1a: - ldetXa[23,9] + ldetWa[22,5] <= 0
 ldetMC1a: - ldetXa[23,9] + ldetWa[22,5] - ldetWa[22,9] >= -1
 ldetMC2a: - ldetXa[23,5] + ldetWa[22,9] <= 0
 ldetMC2a: - ldetXa[23,13] + ldetWa[22,9] <= 0
 ldetMC2a: - ldetXa[23,5] - ldetXa[23,13] + ldetWa[22,9] >= -1
 ldetMC3a: - ldetXa[23,1] + ldetWa[22,13] <= 0
 ldetMC3a: - ldetXa[23,13] + ldetWa[22,13] <= 0
 ldetMC3a: - ldetXa[23,1] - ldetXa[23,13] + ldetWa[22,13] >= -1
 ldetMC0: - lDetX[23,6] + lDetW[22,2] = 0
 ldetMC1: lDetW[22,6] - lDetW[22,10] <= 0
 ldetMC1: - lDetX[23,10] + lDetW[22,6] <= 0
 ldetMC1: - lDetX[23,10] + lDetW[22,6] - lDetW[22,10] >= -1
 ldetMC2: - lDetX[23,6] + lDetW[22,10] <= 0
 ldetMC2: - lDetX[23,14] + lDetW[22,10] <= 0
 ldetMC2: - lDetX[23,6] - lDetX[23,14] + lDetW[22,10] >= -1
 ldetMC3: - lDetX[23,2] + lDetW[22,14] <= 0
 ldetMC3: - lDetX[23,14] + lDetW[22,14] <= 0
 ldetMC3: - lDetX[23,2] - lDetX[23,14] + lDetW[22,14] >= -1
 ldetMC0a: - ldetXa[23,6] + ldetWa[22,2] = 0
 ldetMC1a: ldetWa[22,6] - ldetWa[22,10] <= 0
 ldetMC1a: - ldetXa[23,10] + ldetWa[22,6] <= 0
 ldetMC1a: - ldetXa[23,10] + ldetWa[22,6] - ldetWa[22,10] >= -1
 ldetMC2a: - ldetXa[23,6] + ldetWa[22,10] <= 0
 ldetMC2a: - ldetXa[23,14] + ldetWa[22,10] <= 0
 ldetMC2a: - ldetXa[23,6] - ldetXa[23,14] + ldetWa[22,10] >= -1
 ldetMC3a: - ldetXa[23,2] + ldetWa[22,14] <= 0
 ldetMC3a: - ldetXa[23,14] + ldetWa[22,14] <= 0
 ldetMC3a: - ldetXa[23,2] - ldetXa[23,14] + ldetWa[22,14] >= -1
 ldetMC0: - lDetX[23,7] + lDetW[22,3] = 0
 ldetMC1: lDetW[22,7] - lDetW[22,11] <= 0
 ldetMC1: - lDetX[23,11] + lDetW[22,7] <= 0
 ldetMC1: - lDetX[23,11] + lDetW[22,7] - lDetW[22,11] >= -1
 ldetMC2: - lDetX[23,7] + lDetW[22,11] <= 0
 ldetMC2: - lDetX[23,15] + lDetW[22,11] <= 0
 ldetMC2: - lDetX[23,7] - lDetX[23,15] + lDetW[22,11] >= -1
 ldetMC3: - lDetX[23,3] + lDetW[22,15] <= 0
 ldetMC3: - lDetX[23,15] + lDetW[22,15] <= 0
 ldetMC3: - lDetX[23,3] - lDetX[23,15] + lDetW[22,15] >= -1
 ldetMC0a: - ldetXa[23,7] + ldetWa[22,3] = 0
 ldetMC1a: ldetWa[22,7] - ldetWa[22,11] <= 0
 ldetMC1a: - ldetXa[23,11] + ldetWa[22,7] <= 0
 ldetMC1a: - ldetXa[23,11] + ldetWa[22,7] - ldetWa[22,11] >= -1
 ldetMC2a: - ldetXa[23,7] + ldetWa[22,11] <= 0
 ldetMC2a: - ldetXa[23,15] + ldetWa[22,11] <= 0
 ldetMC2a: - ldetXa[23,7] - ldetXa[23,15] + ldetWa[22,11] >= -1
 ldetMC3a: - ldetXa[23,3] + ldetWa[22,15] <= 0
 ldetMC3a: - ldetXa[23,15] + ldetWa[22,15] <= 0
 ldetMC3a: - ldetXa[23,3] - ldetXa[23,15] + ldetWa[22,15] >= -1
 lfrMC0: lFrMC[19,0] = 0
 lfrMC1: - lDetW[19,4] + lFrMC[19,4] <= 0
 lfrMC1: - ldw[19,8] + lFrMC[19,4] <= 0
 lfrMC1: - ldx[20,8] + lFrMC[19,4] <= 0
 lfrMC1: ldw[19,4] + lFrMC[19,4] <= 1
 lfrMC1: - ldx[20,8] + ldw[19,4] - ldw[19,8] - lDetW[19,4] + lFrMC[19,4]
   >= -2
 lfrMC2: - lDetW[19,8] + lFrMC[19,8] <= 0
 lfrMC2: - ldx[20,4] + lFrMC[19,8] <= 0
 lfrMC2: - ldx[20,12] + lFrMC[19,8] <= 0
 lfrMC2: ldw[19,8] + lFrMC[19,8] <= 1
 lfrMC2: - ldx[20,4] - ldx[20,12] + ldw[19,8] - lDetW[19,8] + lFrMC[19,8]
   >= -2
 lfrMC3: - lDetW[19,12] + lFrMC[19,12] <= 0
 lfrMC3: - ldx[20,0] + lFrMC[19,12] <= 0
 lfrMC3: - ldx[20,12] + lFrMC[19,12] <= 0
 lfrMC3: ldw[19,12] + lFrMC[19,12] <= 1
 lfrMC3: - ldx[20,0] - ldx[20,12] + ldw[19,12] - lDetW[19,12]
   + lFrMC[19,12] >= -2
 lfrMC0a: lfrMCa[19,0] = 0
 lfrMC1a: - ldetWa[19,4] + lfrMCa[19,4] <= 0
 lfrMC1a: - ldw[19,8] + lfrMCa[19,4] <= 0
 lfrMC1a: - ldx[20,8] + lfrMCa[19,4] <= 0
 lfrMC1a: ldw[19,4] + lfrMCa[19,4] <= 1
 lfrMC1a: - ldx[20,8] + ldw[19,4] - ldw[19,8] - ldetWa[19,4] + lfrMCa[19,4]
   >= -2
 lfrMC2a: - ldetWa[19,8] + lfrMCa[19,8] <= 0
 lfrMC2a: - ldx[20,4] + lfrMCa[19,8] <= 0
 lfrMC2a: - ldx[20,12] + lfrMCa[19,8] <= 0
 lfrMC2a: ldw[19,8] + lfrMCa[19,8] <= 1
 lfrMC2a: - ldx[20,4] - ldx[20,12] + ldw[19,8] - ldetWa[19,8]
   + lfrMCa[19,8] >= -2
 lfrMC3a: - ldetWa[19,12] + lfrMCa[19,12] <= 0
 lfrMC3a: - ldx[20,0] + lfrMCa[19,12] <= 0
 lfrMC3a: - ldx[20,12] + lfrMCa[19,12] <= 0
 lfrMC3a: ldw[19,12] + lfrMCa[19,12] <= 1
 lfrMC3a: - ldx[20,0] - ldx[20,12] + ldw[19,12] - ldetWa[19,12]
   + lfrMCa[19,12] >= -2
 lfrMC0: lFrMC[19,1] = 0
 lfrMC1: - lDetW[19,5] + lFrMC[19,5] <= 0
 lfrMC1: - ldw[19,9] + lFrMC[19,5] <= 0
 lfrMC1: - ldx[20,9] + lFrMC[19,5] <= 0
 lfrMC1: ldw[19,5] + lFrMC[19,5] <= 1
 lfrMC1: - ldx[20,9] + ldw[19,5] - ldw[19,9] - lDetW[19,5] + lFrMC[19,5]
   >= -2
 lfrMC2: - lDetW[19,9] + lFrMC[19,9] <= 0
 lfrMC2: - ldx[20,5] + lFrMC[19,9] <= 0
 lfrMC2: - ldx[20,13] + lFrMC[19,9] <= 0
 lfrMC2: ldw[19,9] + lFrMC[19,9] <= 1
 lfrMC2: - ldx[20,5] - ldx[20,13] + ldw[19,9] - lDetW[19,9] + lFrMC[19,9]
   >= -2
 lfrMC3: - lDetW[19,13] + lFrMC[19,13] <= 0
 lfrMC3: - ldx[20,1] + lFrMC[19,13] <= 0
 lfrMC3: - ldx[20,13] + lFrMC[19,13] <= 0
 lfrMC3: ldw[19,13] + lFrMC[19,13] <= 1
 lfrMC3: - ldx[20,1] - ldx[20,13] + ldw[19,13] - lDetW[19,13]
   + lFrMC[19,13] >= -2
 lfrMC0a: lfrMCa[19,1] = 0
 lfrMC1a: - ldetWa[19,5] + lfrMCa[19,5] <= 0
 lfrMC1a: - ldw[19,9] + lfrMCa[19,5] <= 0
 lfrMC1a: - ldx[20,9] + lfrMCa[19,5] <= 0
 lfrMC1a: ldw[19,5] + lfrMCa[19,5] <= 1
 lfrMC1a: - ldx[20,9] + ldw[19,5] - ldw[19,9] - ldetWa[19,5] + lfrMCa[19,5]
   >= -2
 lfrMC2a: - ldetWa[19,9] + lfrMCa[19,9] <= 0
 lfrMC2a: - ldx[20,5] + lfrMCa[19,9] <= 0
 lfrMC2a: - ldx[20,13] + lfrMCa[19,9] <= 0
 lfrMC2a: ldw[19,9] + lfrMCa[19,9] <= 1
 lfrMC2a: - ldx[20,5] - ldx[20,13] + ldw[19,9] - ldetWa[19,9]
   + lfrMCa[19,9] >= -2
 lfrMC3a: - ldetWa[19,13] + lfrMCa[19,13] <= 0
 lfrMC3a: - ldx[20,1] + lfrMCa[19,13] <= 0
 lfrMC3a: - ldx[20,13] + lfrMCa[19,13] <= 0
 lfrMC3a: ldw[19,13] + lfrMCa[19,13] <= 1
 lfrMC3a: - ldx[20,1] - ldx[20,13] + ldw[19,13] - ldetWa[19,13]
   + lfrMCa[19,13] >= -2
 lfrMC0: lFrMC[19,2] = 0
 lfrMC1: - lDetW[19,6] + lFrMC[19,6] <= 0
 lfrMC1: - ldw[19,10] + lFrMC[19,6] <= 0
 lfrMC1: - ldx[20,10] + lFrMC[19,6] <= 0
 lfrMC1: ldw[19,6] + lFrMC[19,6] <= 1
 lfrMC1: - ldx[20,10] + ldw[19,6] - ldw[19,10] - lDetW[19,6] + lFrMC[19,6]
   >= -2
 lfrMC2: - lDetW[19,10] + lFrMC[19,10] <= 0
 lfrMC2: - ldx[20,6] + lFrMC[19,10] <= 0
 lfrMC2: - ldx[20,14] + lFrMC[19,10] <= 0
 lfrMC2: ldw[19,10] + lFrMC[19,10] <= 1
 lfrMC2: - ldx[20,6] - ldx[20,14] + ldw[19,10] - lDetW[19,10]
   + lFrMC[19,10] >= -2
 lfrMC3: - lDetW[19,14] + lFrMC[19,14] <= 0
 lfrMC3: - ldx[20,2] + lFrMC[19,14] <= 0
 lfrMC3: - ldx[20,14] + lFrMC[19,14] <= 0
 lfrMC3: ldw[19,14] + lFrMC[19,14] <= 1
 lfrMC3: - ldx[20,2] - ldx[20,14] + ldw[19,14] - lDetW[19,14]
   + lFrMC[19,14] >= -2
 lfrMC0a: lfrMCa[19,2] = 0
 lfrMC1a: - ldetWa[19,6] + lfrMCa[19,6] <= 0
 lfrMC1a: - ldw[19,10] + lfrMCa[19,6] <= 0
 lfrMC1a: - ldx[20,10] + lfrMCa[19,6] <= 0
 lfrMC1a: ldw[19,6] + lfrMCa[19,6] <= 1
 lfrMC1a: - ldx[20,10] + ldw[19,6] - ldw[19,10] - ldetWa[19,6]
   + lfrMCa[19,6] >= -2
 lfrMC2a: - ldetWa[19,10] + lfrMCa[19,10] <= 0
 lfrMC2a: - ldx[20,6] + lfrMCa[19,10] <= 0
 lfrMC2a: - ldx[20,14] + lfrMCa[19,10] <= 0
 lfrMC2a: ldw[19,10] + lfrMCa[19,10] <= 1
 lfrMC2a: - ldx[20,6] - ldx[20,14] + ldw[19,10] - ldetWa[19,10]
   + lfrMCa[19,10] >= -2
 lfrMC3a: - ldetWa[19,14] + lfrMCa[19,14] <= 0
 lfrMC3a: - ldx[20,2] + lfrMCa[19,14] <= 0
 lfrMC3a: - ldx[20,14] + lfrMCa[19,14] <= 0
 lfrMC3a: ldw[19,14] + lfrMCa[19,14] <= 1
 lfrMC3a: - ldx[20,2] - ldx[20,14] + ldw[19,14] - ldetWa[19,14]
   + lfrMCa[19,14] >= -2
 lfrMC0: lFrMC[19,3] = 0
 lfrMC1: - lDetW[19,7] + lFrMC[19,7] <= 0
 lfrMC1: - ldw[19,11] + lFrMC[19,7] <= 0
 lfrMC1: - ldx[20,11] + lFrMC[19,7] <= 0
 lfrMC1: ldw[19,7] + lFrMC[19,7] <= 1
 lfrMC1: - ldx[20,11] + ldw[19,7] - ldw[19,11] - lDetW[19,7] + lFrMC[19,7]
   >= -2
 lfrMC2: - lDetW[19,11] + lFrMC[19,11] <= 0
 lfrMC2: - ldx[20,7] + lFrMC[19,11] <= 0
 lfrMC2: - ldx[20,15] + lFrMC[19,11] <= 0
 lfrMC2: ldw[19,11] + lFrMC[19,11] <= 1
 lfrMC2: - ldx[20,7] - ldx[20,15] + ldw[19,11] - lDetW[19,11]
   + lFrMC[19,11] >= -2
 lfrMC3: - lDetW[19,15] + lFrMC[19,15] <= 0
 lfrMC3: - ldx[20,3] + lFrMC[19,15] <= 0
 lfrMC3: - ldx[20,15] + lFrMC[19,15] <= 0
 lfrMC3: ldw[19,15] + lFrMC[19,15] <= 1
 lfrMC3: - ldx[20,3] - ldx[20,15] + ldw[19,15] - lDetW[19,15]
   + lFrMC[19,15] >= -2
 lfrMC0a: lfrMCa[19,3] = 0
 lfrMC1a: - ldetWa[19,7] + lfrMCa[19,7] <= 0
 lfrMC1a: - ldw[19,11] + lfrMCa[19,7] <= 0
 lfrMC1a: - ldx[20,11] + lfrMCa[19,7] <= 0
 lfrMC1a: ldw[19,7] + lfrMCa[19,7] <= 1
 lfrMC1a: - ldx[20,11] + ldw[19,7] - ldw[19,11] - ldetWa[19,7]
   + lfrMCa[19,7] >= -2
 lfrMC2a: - ldetWa[19,11] + lfrMCa[19,11] <= 0
 lfrMC2a: - ldx[20,7] + lfrMCa[19,11] <= 0
 lfrMC2a: - ldx[20,15] + lfrMCa[19,11] <= 0
 lfrMC2a: ldw[19,11] + lfrMCa[19,11] <= 1
 lfrMC2a: - ldx[20,7] - ldx[20,15] + ldw[19,11] - ldetWa[19,11]
   + lfrMCa[19,11] >= -2
 lfrMC3a: - ldetWa[19,15] + lfrMCa[19,15] <= 0
 lfrMC3a: - ldx[20,3] + lfrMCa[19,15] <= 0
 lfrMC3a: - ldx[20,15] + lfrMCa[19,15] <= 0
 lfrMC3a: ldw[19,15] + lfrMCa[19,15] <= 1
 lfrMC3a: - ldx[20,3] - ldx[20,15] + ldw[19,15] - ldetWa[19,15]
   + lfrMCa[19,15] >= -2
 lfrMC0: lFrMC[20,0] = 0
 lfrMC1: - lDetW[20,4] + lFrMC[20,4] <= 0
 lfrMC1: - ldw[20,8] + lFrMC[20,4] <= 0
 lfrMC1: - ldx[21,8] + lFrMC[20,4] <= 0
 lfrMC1: ldw[20,4] + lFrMC[20,4] <= 1
 lfrMC1: - ldx[21,8] + ldw[20,4] - ldw[20,8] - lDetW[20,4] + lFrMC[20,4]
   >= -2
 lfrMC2: - lDetW[20,8] + lFrMC[20,8] <= 0
 lfrMC2: - ldx[21,4] + lFrMC[20,8] <= 0
 lfrMC2: - ldx[21,12] + lFrMC[20,8] <= 0
 lfrMC2: ldw[20,8] + lFrMC[20,8] <= 1
 lfrMC2: - ldx[21,4] - ldx[21,12] + ldw[20,8] - lDetW[20,8] + lFrMC[20,8]
   >= -2
 lfrMC3: - lDetW[20,12] + lFrMC[20,12] <= 0
 lfrMC3: - ldx[21,0] + lFrMC[20,12] <= 0
 lfrMC3: - ldx[21,12] + lFrMC[20,12] <= 0
 lfrMC3: ldw[20,12] + lFrMC[20,12] <= 1
 lfrMC3: - ldx[21,0] - ldx[21,12] + ldw[20,12] - lDetW[20,12]
   + lFrMC[20,12] >= -2
 lfrMC0a: lfrMCa[20,0] = 0
 lfrMC1a: - ldetWa[20,4] + lfrMCa[20,4] <= 0
 lfrMC1a: - ldw[20,8] + lfrMCa[20,4] <= 0
 lfrMC1a: - ldx[21,8] + lfrMCa[20,4] <= 0
 lfrMC1a: ldw[20,4] + lfrMCa[20,4] <= 1
 lfrMC1a: - ldx[21,8] + ldw[20,4] - ldw[20,8] - ldetWa[20,4] + lfrMCa[20,4]
   >= -2
 lfrMC2a: - ldetWa[20,8] + lfrMCa[20,8] <= 0
 lfrMC2a: - ldx[21,4] + lfrMCa[20,8] <= 0
 lfrMC2a: - ldx[21,12] + lfrMCa[20,8] <= 0
 lfrMC2a: ldw[20,8] + lfrMCa[20,8] <= 1
 lfrMC2a: - ldx[21,4] - ldx[21,12] + ldw[20,8] - ldetWa[20,8]
   + lfrMCa[20,8] >= -2
 lfrMC3a: - ldetWa[20,12] + lfrMCa[20,12] <= 0
 lfrMC3a: - ldx[21,0] + lfrMCa[20,12] <= 0
 lfrMC3a: - ldx[21,12] + lfrMCa[20,12] <= 0
 lfrMC3a: ldw[20,12] + lfrMCa[20,12] <= 1
 lfrMC3a: - ldx[21,0] - ldx[21,12] + ldw[20,12] - ldetWa[20,12]
   + lfrMCa[20,12] >= -2
 lfrMC0: lFrMC[20,1] = 0
 lfrMC1: - lDetW[20,5] + lFrMC[20,5] <= 0
 lfrMC1: - ldw[20,9] + lFrMC[20,5] <= 0
 lfrMC1: - ldx[21,9] + lFrMC[20,5] <= 0
 lfrMC1: ldw[20,5] + lFrMC[20,5] <= 1
 lfrMC1: - ldx[21,9] + ldw[20,5] - ldw[20,9] - lDetW[20,5] + lFrMC[20,5]
   >= -2
 lfrMC2: - lDetW[20,9] + lFrMC[20,9] <= 0
 lfrMC2: - ldx[21,5] + lFrMC[20,9] <= 0
 lfrMC2: - ldx[21,13] + lFrMC[20,9] <= 0
 lfrMC2: ldw[20,9] + lFrMC[20,9] <= 1
 lfrMC2: - ldx[21,5] - ldx[21,13] + ldw[20,9] - lDetW[20,9] + lFrMC[20,9]
   >= -2
 lfrMC3: - lDetW[20,13] + lFrMC[20,13] <= 0
 lfrMC3: - ldx[21,1] + lFrMC[20,13] <= 0
 lfrMC3: - ldx[21,13] + lFrMC[20,13] <= 0
 lfrMC3: ldw[20,13] + lFrMC[20,13] <= 1
 lfrMC3: - ldx[21,1] - ldx[21,13] + ldw[20,13] - lDetW[20,13]
   + lFrMC[20,13] >= -2
 lfrMC0a: lfrMCa[20,1] = 0
 lfrMC1a: - ldetWa[20,5] + lfrMCa[20,5] <= 0
 lfrMC1a: - ldw[20,9] + lfrMCa[20,5] <= 0
 lfrMC1a: - ldx[21,9] + lfrMCa[20,5] <= 0
 lfrMC1a: ldw[20,5] + lfrMCa[20,5] <= 1
 lfrMC1a: - ldx[21,9] + ldw[20,5] - ldw[20,9] - ldetWa[20,5] + lfrMCa[20,5]
   >= -2
 lfrMC2a: - ldetWa[20,9] + lfrMCa[20,9] <= 0
 lfrMC2a: - ldx[21,5] + lfrMCa[20,9] <= 0
 lfrMC2a: - ldx[21,13] + lfrMCa[20,9] <= 0
 lfrMC2a: ldw[20,9] + lfrMCa[20,9] <= 1
 lfrMC2a: - ldx[21,5] - ldx[21,13] + ldw[20,9] - ldetWa[20,9]
   + lfrMCa[20,9] >= -2
 lfrMC3a: - ldetWa[20,13] + lfrMCa[20,13] <= 0
 lfrMC3a: - ldx[21,1] + lfrMCa[20,13] <= 0
 lfrMC3a: - ldx[21,13] + lfrMCa[20,13] <= 0
 lfrMC3a: ldw[20,13] + lfrMCa[20,13] <= 1
 lfrMC3a: - ldx[21,1] - ldx[21,13] + ldw[20,13] - ldetWa[20,13]
   + lfrMCa[20,13] >= -2
 lfrMC0: lFrMC[20,2] = 0
 lfrMC1: - lDetW[20,6] + lFrMC[20,6] <= 0
 lfrMC1: - ldw[20,10] + lFrMC[20,6] <= 0
 lfrMC1: - ldx[21,10] + lFrMC[20,6] <= 0
 lfrMC1: ldw[20,6] + lFrMC[20,6] <= 1
 lfrMC1: - ldx[21,10] + ldw[20,6] - ldw[20,10] - lDetW[20,6] + lFrMC[20,6]
   >= -2
 lfrMC2: - lDetW[20,10] + lFrMC[20,10] <= 0
 lfrMC2: - ldx[21,6] + lFrMC[20,10] <= 0
 lfrMC2: - ldx[21,14] + lFrMC[20,10] <= 0
 lfrMC2: ldw[20,10] + lFrMC[20,10] <= 1
 lfrMC2: - ldx[21,6] - ldx[21,14] + ldw[20,10] - lDetW[20,10]
   + lFrMC[20,10] >= -2
 lfrMC3: - lDetW[20,14] + lFrMC[20,14] <= 0
 lfrMC3: - ldx[21,2] + lFrMC[20,14] <= 0
 lfrMC3: - ldx[21,14] + lFrMC[20,14] <= 0
 lfrMC3: ldw[20,14] + lFrMC[20,14] <= 1
 lfrMC3: - ldx[21,2] - ldx[21,14] + ldw[20,14] - lDetW[20,14]
   + lFrMC[20,14] >= -2
 lfrMC0a: lfrMCa[20,2] = 0
 lfrMC1a: - ldetWa[20,6] + lfrMCa[20,6] <= 0
 lfrMC1a: - ldw[20,10] + lfrMCa[20,6] <= 0
 lfrMC1a: - ldx[21,10] + lfrMCa[20,6] <= 0
 lfrMC1a: ldw[20,6] + lfrMCa[20,6] <= 1
 lfrMC1a: - ldx[21,10] + ldw[20,6] - ldw[20,10] - ldetWa[20,6]
   + lfrMCa[20,6] >= -2
 lfrMC2a: - ldetWa[20,10] + lfrMCa[20,10] <= 0
 lfrMC2a: - ldx[21,6] + lfrMCa[20,10] <= 0
 lfrMC2a: - ldx[21,14] + lfrMCa[20,10] <= 0
 lfrMC2a: ldw[20,10] + lfrMCa[20,10] <= 1
 lfrMC2a: - ldx[21,6] - ldx[21,14] + ldw[20,10] - ldetWa[20,10]
   + lfrMCa[20,10] >= -2
 lfrMC3a: - ldetWa[20,14] + lfrMCa[20,14] <= 0
 lfrMC3a: - ldx[21,2] + lfrMCa[20,14] <= 0
 lfrMC3a: - ldx[21,14] + lfrMCa[20,14] <= 0
 lfrMC3a: ldw[20,14] + lfrMCa[20,14] <= 1
 lfrMC3a: - ldx[21,2] - ldx[21,14] + ldw[20,14] - ldetWa[20,14]
   + lfrMCa[20,14] >= -2
 lfrMC0: lFrMC[20,3] = 0
 lfrMC1: - lDetW[20,7] + lFrMC[20,7] <= 0
 lfrMC1: - ldw[20,11] + lFrMC[20,7] <= 0
 lfrMC1: - ldx[21,11] + lFrMC[20,7] <= 0
 lfrMC1: ldw[20,7] + lFrMC[20,7] <= 1
 lfrMC1: - ldx[21,11] + ldw[20,7] - ldw[20,11] - lDetW[20,7] + lFrMC[20,7]
   >= -2
 lfrMC2: - lDetW[20,11] + lFrMC[20,11] <= 0
 lfrMC2: - ldx[21,7] + lFrMC[20,11] <= 0
 lfrMC2: - ldx[21,15] + lFrMC[20,11] <= 0
 lfrMC2: ldw[20,11] + lFrMC[20,11] <= 1
 lfrMC2: - ldx[21,7] - ldx[21,15] + ldw[20,11] - lDetW[20,11]
   + lFrMC[20,11] >= -2
 lfrMC3: - lDetW[20,15] + lFrMC[20,15] <= 0
 lfrMC3: - ldx[21,3] + lFrMC[20,15] <= 0
 lfrMC3: - ldx[21,15] + lFrMC[20,15] <= 0
 lfrMC3: ldw[20,15] + lFrMC[20,15] <= 1
 lfrMC3: - ldx[21,3] - ldx[21,15] + ldw[20,15] - lDetW[20,15]
   + lFrMC[20,15] >= -2
 lfrMC0a: lfrMCa[20,3] = 0
 lfrMC1a: - ldetWa[20,7] + lfrMCa[20,7] <= 0
 lfrMC1a: - ldw[20,11] + lfrMCa[20,7] <= 0
 lfrMC1a: - ldx[21,11] + lfrMCa[20,7] <= 0
 lfrMC1a: ldw[20,7] + lfrMCa[20,7] <= 1
 lfrMC1a: - ldx[21,11] + ldw[20,7] - ldw[20,11] - ldetWa[20,7]
   + lfrMCa[20,7] >= -2
 lfrMC2a: - ldetWa[20,11] + lfrMCa[20,11] <= 0
 lfrMC2a: - ldx[21,7] + lfrMCa[20,11] <= 0
 lfrMC2a: - ldx[21,15] + lfrMCa[20,11] <= 0
 lfrMC2a: ldw[20,11] + lfrMCa[20,11] <= 1
 lfrMC2a: - ldx[21,7] - ldx[21,15] + ldw[20,11] - ldetWa[20,11]
   + lfrMCa[20,11] >= -2
 lfrMC3a: - ldetWa[20,15] + lfrMCa[20,15] <= 0
 lfrMC3a: - ldx[21,3] + lfrMCa[20,15] <= 0
 lfrMC3a: - ldx[21,15] + lfrMCa[20,15] <= 0
 lfrMC3a: ldw[20,15] + lfrMCa[20,15] <= 1
 lfrMC3a: - ldx[21,3] - ldx[21,15] + ldw[20,15] - ldetWa[20,15]
   + lfrMCa[20,15] >= -2
 lfrMC0: lFrMC[21,0] = 0
 lfrMC1: - lDetW[21,4] + lFrMC[21,4] <= 0
 lfrMC1: - ldw[21,8] + lFrMC[21,4] <= 0
 lfrMC1: - ldx[22,8] + lFrMC[21,4] <= 0
 lfrMC1: ldw[21,4] + lFrMC[21,4] <= 1
 lfrMC1: - ldx[22,8] + ldw[21,4] - ldw[21,8] - lDetW[21,4] + lFrMC[21,4]
   >= -2
 lfrMC2: - lDetW[21,8] + lFrMC[21,8] <= 0
 lfrMC2: - ldx[22,4] + lFrMC[21,8] <= 0
 lfrMC2: - ldx[22,12] + lFrMC[21,8] <= 0
 lfrMC2: ldw[21,8] + lFrMC[21,8] <= 1
 lfrMC2: - ldx[22,4] - ldx[22,12] + ldw[21,8] - lDetW[21,8] + lFrMC[21,8]
   >= -2
 lfrMC3: - lDetW[21,12] + lFrMC[21,12] <= 0
 lfrMC3: - ldx[22,0] + lFrMC[21,12] <= 0
 lfrMC3: - ldx[22,12] + lFrMC[21,12] <= 0
 lfrMC3: ldw[21,12] + lFrMC[21,12] <= 1
 lfrMC3: - ldx[22,0] - ldx[22,12] + ldw[21,12] - lDetW[21,12]
   + lFrMC[21,12] >= -2
 lfrMC0a: lfrMCa[21,0] = 0
 lfrMC1a: - ldetWa[21,4] + lfrMCa[21,4] <= 0
 lfrMC1a: - ldw[21,8] + lfrMCa[21,4] <= 0
 lfrMC1a: - ldx[22,8] + lfrMCa[21,4] <= 0
 lfrMC1a: ldw[21,4] + lfrMCa[21,4] <= 1
 lfrMC1a: - ldx[22,8] + ldw[21,4] - ldw[21,8] - ldetWa[21,4] + lfrMCa[21,4]
   >= -2
 lfrMC2a: - ldetWa[21,8] + lfrMCa[21,8] <= 0
 lfrMC2a: - ldx[22,4] + lfrMCa[21,8] <= 0
 lfrMC2a: - ldx[22,12] + lfrMCa[21,8] <= 0
 lfrMC2a: ldw[21,8] + lfrMCa[21,8] <= 1
 lfrMC2a: - ldx[22,4] - ldx[22,12] + ldw[21,8] - ldetWa[21,8]
   + lfrMCa[21,8] >= -2
 lfrMC3a: - ldetWa[21,12] + lfrMCa[21,12] <= 0
 lfrMC3a: - ldx[22,0] + lfrMCa[21,12] <= 0
 lfrMC3a: - ldx[22,12] + lfrMCa[21,12] <= 0
 lfrMC3a: ldw[21,12] + lfrMCa[21,12] <= 1
 lfrMC3a: - ldx[22,0] - ldx[22,12] + ldw[21,12] - ldetWa[21,12]
   + lfrMCa[21,12] >= -2
 lfrMC0: lFrMC[21,1] = 0
 lfrMC1: - lDetW[21,5] + lFrMC[21,5] <= 0
 lfrMC1: - ldw[21,9] + lFrMC[21,5] <= 0
 lfrMC1: - ldx[22,9] + lFrMC[21,5] <= 0
 lfrMC1: ldw[21,5] + lFrMC[21,5] <= 1
 lfrMC1: - ldx[22,9] + ldw[21,5] - ldw[21,9] - lDetW[21,5] + lFrMC[21,5]
   >= -2
 lfrMC2: - lDetW[21,9] + lFrMC[21,9] <= 0
 lfrMC2: - ldx[22,5] + lFrMC[21,9] <= 0
 lfrMC2: - ldx[22,13] + lFrMC[21,9] <= 0
 lfrMC2: ldw[21,9] + lFrMC[21,9] <= 1
 lfrMC2: - ldx[22,5] - ldx[22,13] + ldw[21,9] - lDetW[21,9] + lFrMC[21,9]
   >= -2
 lfrMC3: - lDetW[21,13] + lFrMC[21,13] <= 0
 lfrMC3: - ldx[22,1] + lFrMC[21,13] <= 0
 lfrMC3: - ldx[22,13] + lFrMC[21,13] <= 0
 lfrMC3: ldw[21,13] + lFrMC[21,13] <= 1
 lfrMC3: - ldx[22,1] - ldx[22,13] + ldw[21,13] - lDetW[21,13]
   + lFrMC[21,13] >= -2
 lfrMC0a: lfrMCa[21,1] = 0
 lfrMC1a: - ldetWa[21,5] + lfrMCa[21,5] <= 0
 lfrMC1a: - ldw[21,9] + lfrMCa[21,5] <= 0
 lfrMC1a: - ldx[22,9] + lfrMCa[21,5] <= 0
 lfrMC1a: ldw[21,5] + lfrMCa[21,5] <= 1
 lfrMC1a: - ldx[22,9] + ldw[21,5] - ldw[21,9] - ldetWa[21,5] + lfrMCa[21,5]
   >= -2
 lfrMC2a: - ldetWa[21,9] + lfrMCa[21,9] <= 0
 lfrMC2a: - ldx[22,5] + lfrMCa[21,9] <= 0
 lfrMC2a: - ldx[22,13] + lfrMCa[21,9] <= 0
 lfrMC2a: ldw[21,9] + lfrMCa[21,9] <= 1
 lfrMC2a: - ldx[22,5] - ldx[22,13] + ldw[21,9] - ldetWa[21,9]
   + lfrMCa[21,9] >= -2
 lfrMC3a: - ldetWa[21,13] + lfrMCa[21,13] <= 0
 lfrMC3a: - ldx[22,1] + lfrMCa[21,13] <= 0
 lfrMC3a: - ldx[22,13] + lfrMCa[21,13] <= 0
 lfrMC3a: ldw[21,13] + lfrMCa[21,13] <= 1
 lfrMC3a: - ldx[22,1] - ldx[22,13] + ldw[21,13] - ldetWa[21,13]
   + lfrMCa[21,13] >= -2
 lfrMC0: lFrMC[21,2] = 0
 lfrMC1: - lDetW[21,6] + lFrMC[21,6] <= 0
 lfrMC1: - ldw[21,10] + lFrMC[21,6] <= 0
 lfrMC1: - ldx[22,10] + lFrMC[21,6] <= 0
 lfrMC1: ldw[21,6] + lFrMC[21,6] <= 1
 lfrMC1: - ldx[22,10] + ldw[21,6] - ldw[21,10] - lDetW[21,6] + lFrMC[21,6]
   >= -2
 lfrMC2: - lDetW[21,10] + lFrMC[21,10] <= 0
 lfrMC2: - ldx[22,6] + lFrMC[21,10] <= 0
 lfrMC2: - ldx[22,14] + lFrMC[21,10] <= 0
 lfrMC2: ldw[21,10] + lFrMC[21,10] <= 1
 lfrMC2: - ldx[22,6] - ldx[22,14] + ldw[21,10] - lDetW[21,10]
   + lFrMC[21,10] >= -2
 lfrMC3: - lDetW[21,14] + lFrMC[21,14] <= 0
 lfrMC3: - ldx[22,2] + lFrMC[21,14] <= 0
 lfrMC3: - ldx[22,14] + lFrMC[21,14] <= 0
 lfrMC3: ldw[21,14] + lFrMC[21,14] <= 1
 lfrMC3: - ldx[22,2] - ldx[22,14] + ldw[21,14] - lDetW[21,14]
   + lFrMC[21,14] >= -2
 lfrMC0a: lfrMCa[21,2] = 0
 lfrMC1a: - ldetWa[21,6] + lfrMCa[21,6] <= 0
 lfrMC1a: - ldw[21,10] + lfrMCa[21,6] <= 0
 lfrMC1a: - ldx[22,10] + lfrMCa[21,6] <= 0
 lfrMC1a: ldw[21,6] + lfrMCa[21,6] <= 1
 lfrMC1a: - ldx[22,10] + ldw[21,6] - ldw[21,10] - ldetWa[21,6]
   + lfrMCa[21,6] >= -2
 lfrMC2a: - ldetWa[21,10] + lfrMCa[21,10] <= 0
 lfrMC2a: - ldx[22,6] + lfrMCa[21,10] <= 0
 lfrMC2a: - ldx[22,14] + lfrMCa[21,10] <= 0
 lfrMC2a: ldw[21,10] + lfrMCa[21,10] <= 1
 lfrMC2a: - ldx[22,6] - ldx[22,14] + ldw[21,10] - ldetWa[21,10]
   + lfrMCa[21,10] >= -2
 lfrMC3a: - ldetWa[21,14] + lfrMCa[21,14] <= 0
 lfrMC3a: - ldx[22,2] + lfrMCa[21,14] <= 0
 lfrMC3a: - ldx[22,14] + lfrMCa[21,14] <= 0
 lfrMC3a: ldw[21,14] + lfrMCa[21,14] <= 1
 lfrMC3a: - ldx[22,2] - ldx[22,14] + ldw[21,14] - ldetWa[21,14]
   + lfrMCa[21,14] >= -2
 lfrMC0: lFrMC[21,3] = 0
 lfrMC1: - lDetW[21,7] + lFrMC[21,7] <= 0
 lfrMC1: - ldw[21,11] + lFrMC[21,7] <= 0
 lfrMC1: - ldx[22,11] + lFrMC[21,7] <= 0
 lfrMC1: ldw[21,7] + lFrMC[21,7] <= 1
 lfrMC1: - ldx[22,11] + ldw[21,7] - ldw[21,11] - lDetW[21,7] + lFrMC[21,7]
   >= -2
 lfrMC2: - lDetW[21,11] + lFrMC[21,11] <= 0
 lfrMC2: - ldx[22,7] + lFrMC[21,11] <= 0
 lfrMC2: - ldx[22,15] + lFrMC[21,11] <= 0
 lfrMC2: ldw[21,11] + lFrMC[21,11] <= 1
 lfrMC2: - ldx[22,7] - ldx[22,15] + ldw[21,11] - lDetW[21,11]
   + lFrMC[21,11] >= -2
 lfrMC3: - lDetW[21,15] + lFrMC[21,15] <= 0
 lfrMC3: - ldx[22,3] + lFrMC[21,15] <= 0
 lfrMC3: - ldx[22,15] + lFrMC[21,15] <= 0
 lfrMC3: ldw[21,15] + lFrMC[21,15] <= 1
 lfrMC3: - ldx[22,3] - ldx[22,15] + ldw[21,15] - lDetW[21,15]
   + lFrMC[21,15] >= -2
 lfrMC0a: lfrMCa[21,3] = 0
 lfrMC1a: - ldetWa[21,7] + lfrMCa[21,7] <= 0
 lfrMC1a: - ldw[21,11] + lfrMCa[21,7] <= 0
 lfrMC1a: - ldx[22,11] + lfrMCa[21,7] <= 0
 lfrMC1a: ldw[21,7] + lfrMCa[21,7] <= 1
 lfrMC1a: - ldx[22,11] + ldw[21,7] - ldw[21,11] - ldetWa[21,7]
   + lfrMCa[21,7] >= -2
 lfrMC2a: - ldetWa[21,11] + lfrMCa[21,11] <= 0
 lfrMC2a: - ldx[22,7] + lfrMCa[21,11] <= 0
 lfrMC2a: - ldx[22,15] + lfrMCa[21,11] <= 0
 lfrMC2a: ldw[21,11] + lfrMCa[21,11] <= 1
 lfrMC2a: - ldx[22,7] - ldx[22,15] + ldw[21,11] - ldetWa[21,11]
   + lfrMCa[21,11] >= -2
 lfrMC3a: - ldetWa[21,15] + lfrMCa[21,15] <= 0
 lfrMC3a: - ldx[22,3] + lfrMCa[21,15] <= 0
 lfrMC3a: - ldx[22,15] + lfrMCa[21,15] <= 0
 lfrMC3a: ldw[21,15] + lfrMCa[21,15] <= 1
 lfrMC3a: - ldx[22,3] - ldx[22,15] + ldw[21,15] - ldetWa[21,15]
   + lfrMCa[21,15] >= -2
 lfrMC0: lFrMC[22,0] = 0
 lfrMC1: - lDetW[22,4] + lFrMC[22,4] <= 0
 lfrMC1: - ldw[22,8] + lFrMC[22,4] <= 0
 lfrMC1: - ldx[23,8] + lFrMC[22,4] <= 0
 lfrMC1: ldw[22,4] + lFrMC[22,4] <= 1
 lfrMC1: - ldx[23,8] + ldw[22,4] - ldw[22,8] - lDetW[22,4] + lFrMC[22,4]
   >= -2
 lfrMC2: - lDetW[22,8] + lFrMC[22,8] <= 0
 lfrMC2: - ldx[23,4] + lFrMC[22,8] <= 0
 lfrMC2: - ldx[23,12] + lFrMC[22,8] <= 0
 lfrMC2: ldw[22,8] + lFrMC[22,8] <= 1
 lfrMC2: - ldx[23,4] - ldx[23,12] + ldw[22,8] - lDetW[22,8] + lFrMC[22,8]
   >= -2
 lfrMC3: - lDetW[22,12] + lFrMC[22,12] <= 0
 lfrMC3: - ldx[23,0] + lFrMC[22,12] <= 0
 lfrMC3: - ldx[23,12] + lFrMC[22,12] <= 0
 lfrMC3: ldw[22,12] + lFrMC[22,12] <= 1
 lfrMC3: - ldx[23,0] - ldx[23,12] + ldw[22,12] - lDetW[22,12]
   + lFrMC[22,12] >= -2
 lfrMC0a: lfrMCa[22,0] = 0
 lfrMC1a: - ldetWa[22,4] + lfrMCa[22,4] <= 0
 lfrMC1a: - ldw[22,8] + lfrMCa[22,4] <= 0
 lfrMC1a: - ldx[23,8] + lfrMCa[22,4] <= 0
 lfrMC1a: ldw[22,4] + lfrMCa[22,4] <= 1
 lfrMC1a: - ldx[23,8] + ldw[22,4] - ldw[22,8] - ldetWa[22,4] + lfrMCa[22,4]
   >= -2
 lfrMC2a: - ldetWa[22,8] + lfrMCa[22,8] <= 0
 lfrMC2a: - ldx[23,4] + lfrMCa[22,8] <= 0
 lfrMC2a: - ldx[23,12] + lfrMCa[22,8] <= 0
 lfrMC2a: ldw[22,8] + lfrMCa[22,8] <= 1
 lfrMC2a: - ldx[23,4] - ldx[23,12] + ldw[22,8] - ldetWa[22,8]
   + lfrMCa[22,8] >= -2
 lfrMC3a: - ldetWa[22,12] + lfrMCa[22,12] <= 0
 lfrMC3a: - ldx[23,0] + lfrMCa[22,12] <= 0
 lfrMC3a: - ldx[23,12] + lfrMCa[22,12] <= 0
 lfrMC3a: ldw[22,12] + lfrMCa[22,12] <= 1
 lfrMC3a: - ldx[23,0] - ldx[23,12] + ldw[22,12] - ldetWa[22,12]
   + lfrMCa[22,12] >= -2
 lfrMC0: lFrMC[22,1] = 0
 lfrMC1: - lDetW[22,5] + lFrMC[22,5] <= 0
 lfrMC1: - ldw[22,9] + lFrMC[22,5] <= 0
 lfrMC1: - ldx[23,9] + lFrMC[22,5] <= 0
 lfrMC1: ldw[22,5] + lFrMC[22,5] <= 1
 lfrMC1: - ldx[23,9] + ldw[22,5] - ldw[22,9] - lDetW[22,5] + lFrMC[22,5]
   >= -2
 lfrMC2: - lDetW[22,9] + lFrMC[22,9] <= 0
 lfrMC2: - ldx[23,5] + lFrMC[22,9] <= 0
 lfrMC2: - ldx[23,13] + lFrMC[22,9] <= 0
 lfrMC2: ldw[22,9] + lFrMC[22,9] <= 1
 lfrMC2: - ldx[23,5] - ldx[23,13] + ldw[22,9] - lDetW[22,9] + lFrMC[22,9]
   >= -2
 lfrMC3: - lDetW[22,13] + lFrMC[22,13] <= 0
 lfrMC3: - ldx[23,1] + lFrMC[22,13] <= 0
 lfrMC3: - ldx[23,13] + lFrMC[22,13] <= 0
 lfrMC3: ldw[22,13] + lFrMC[22,13] <= 1
 lfrMC3: - ldx[23,1] - ldx[23,13] + ldw[22,13] - lDetW[22,13]
   + lFrMC[22,13] >= -2
 lfrMC0a: lfrMCa[22,1] = 0
 lfrMC1a: - ldetWa[22,5] + lfrMCa[22,5] <= 0
 lfrMC1a: - ldw[22,9] + lfrMCa[22,5] <= 0
 lfrMC1a: - ldx[23,9] + lfrMCa[22,5] <= 0
 lfrMC1a: ldw[22,5] + lfrMCa[22,5] <= 1
 lfrMC1a: - ldx[23,9] + ldw[22,5] - ldw[22,9] - ldetWa[22,5] + lfrMCa[22,5]
   >= -2
 lfrMC2a: - ldetWa[22,9] + lfrMCa[22,9] <= 0
 lfrMC2a: - ldx[23,5] + lfrMCa[22,9] <= 0
 lfrMC2a: - ldx[23,13] + lfrMCa[22,9] <= 0
 lfrMC2a: ldw[22,9] + lfrMCa[22,9] <= 1
 lfrMC2a: - ldx[23,5] - ldx[23,13] + ldw[22,9] - ldetWa[22,9]
   + lfrMCa[22,9] >= -2
 lfrMC3a: - ldetWa[22,13] + lfrMCa[22,13] <= 0
 lfrMC3a: - ldx[23,1] + lfrMCa[22,13] <= 0
 lfrMC3a: - ldx[23,13] + lfrMCa[22,13] <= 0
 lfrMC3a: ldw[22,13] + lfrMCa[22,13] <= 1
 lfrMC3a: - ldx[23,1] - ldx[23,13] + ldw[22,13] - ldetWa[22,13]
   + lfrMCa[22,13] >= -2
 lfrMC0: lFrMC[22,2] = 0
 lfrMC1: - lDetW[22,6] + lFrMC[22,6] <= 0
 lfrMC1: - ldw[22,10] + lFrMC[22,6] <= 0
 lfrMC1: - ldx[23,10] + lFrMC[22,6] <= 0
 lfrMC1: ldw[22,6] + lFrMC[22,6] <= 1
 lfrMC1: - ldx[23,10] + ldw[22,6] - ldw[22,10] - lDetW[22,6] + lFrMC[22,6]
   >= -2
 lfrMC2: - lDetW[22,10] + lFrMC[22,10] <= 0
 lfrMC2: - ldx[23,6] + lFrMC[22,10] <= 0
 lfrMC2: - ldx[23,14] + lFrMC[22,10] <= 0
 lfrMC2: ldw[22,10] + lFrMC[22,10] <= 1
 lfrMC2: - ldx[23,6] - ldx[23,14] + ldw[22,10] - lDetW[22,10]
   + lFrMC[22,10] >= -2
 lfrMC3: - lDetW[22,14] + lFrMC[22,14] <= 0
 lfrMC3: - ldx[23,2] + lFrMC[22,14] <= 0
 lfrMC3: - ldx[23,14] + lFrMC[22,14] <= 0
 lfrMC3: ldw[22,14] + lFrMC[22,14] <= 1
 lfrMC3: - ldx[23,2] - ldx[23,14] + ldw[22,14] - lDetW[22,14]
   + lFrMC[22,14] >= -2
 lfrMC0a: lfrMCa[22,2] = 0
 lfrMC1a: - ldetWa[22,6] + lfrMCa[22,6] <= 0
 lfrMC1a: - ldw[22,10] + lfrMCa[22,6] <= 0
 lfrMC1a: - ldx[23,10] + lfrMCa[22,6] <= 0
 lfrMC1a: ldw[22,6] + lfrMCa[22,6] <= 1
 lfrMC1a: - ldx[23,10] + ldw[22,6] - ldw[22,10] - ldetWa[22,6]
   + lfrMCa[22,6] >= -2
 lfrMC2a: - ldetWa[22,10] + lfrMCa[22,10] <= 0
 lfrMC2a: - ldx[23,6] + lfrMCa[22,10] <= 0
 lfrMC2a: - ldx[23,14] + lfrMCa[22,10] <= 0
 lfrMC2a: ldw[22,10] + lfrMCa[22,10] <= 1
 lfrMC2a: - ldx[23,6] - ldx[23,14] + ldw[22,10] - ldetWa[22,10]
   + lfrMCa[22,10] >= -2
 lfrMC3a: - ldetWa[22,14] + lfrMCa[22,14] <= 0
 lfrMC3a: - ldx[23,2] + lfrMCa[22,14] <= 0
 lfrMC3a: - ldx[23,14] + lfrMCa[22,14] <= 0
 lfrMC3a: ldw[22,14] + lfrMCa[22,14] <= 1
 lfrMC3a: - ldx[23,2] - ldx[23,14] + ldw[22,14] - ldetWa[22,14]
   + lfrMCa[22,14] >= -2
 lfrMC0: lFrMC[22,3] = 0
 lfrMC1: - lDetW[22,7] + lFrMC[22,7] <= 0
 lfrMC1: - ldw[22,11] + lFrMC[22,7] <= 0
 lfrMC1: - ldx[23,11] + lFrMC[22,7] <= 0
 lfrMC1: ldw[22,7] + lFrMC[22,7] <= 1
 lfrMC1: - ldx[23,11] + ldw[22,7] - ldw[22,11] - lDetW[22,7] + lFrMC[22,7]
   >= -2
 lfrMC2: - lDetW[22,11] + lFrMC[22,11] <= 0
 lfrMC2: - ldx[23,7] + lFrMC[22,11] <= 0
 lfrMC2: - ldx[23,15] + lFrMC[22,11] <= 0
 lfrMC2: ldw[22,11] + lFrMC[22,11] <= 1
 lfrMC2: - ldx[23,7] - ldx[23,15] + ldw[22,11] - lDetW[22,11]
   + lFrMC[22,11] >= -2
 lfrMC3: - lDetW[22,15] + lFrMC[22,15] <= 0
 lfrMC3: - ldx[23,3] + lFrMC[22,15] <= 0
 lfrMC3: - ldx[23,15] + lFrMC[22,15] <= 0
 lfrMC3: ldw[22,15] + lFrMC[22,15] <= 1
 lfrMC3: - ldx[23,3] - ldx[23,15] + ldw[22,15] - lDetW[22,15]
   + lFrMC[22,15] >= -2
 lfrMC0a: lfrMCa[22,3] = 0
 lfrMC1a: - ldetWa[22,7] + lfrMCa[22,7] <= 0
 lfrMC1a: - ldw[22,11] + lfrMCa[22,7] <= 0
 lfrMC1a: - ldx[23,11] + lfrMCa[22,7] <= 0
 lfrMC1a: ldw[22,7] + lfrMCa[22,7] <= 1
 lfrMC1a: - ldx[23,11] + ldw[22,7] - ldw[22,11] - ldetWa[22,7]
   + lfrMCa[22,7] >= -2
 lfrMC2a: - ldetWa[22,11] + lfrMCa[22,11] <= 0
 lfrMC2a: - ldx[23,7] + lfrMCa[22,11] <= 0
 lfrMC2a: - ldx[23,15] + lfrMCa[22,11] <= 0
 lfrMC2a: ldw[22,11] + lfrMCa[22,11] <= 1
 lfrMC2a: - ldx[23,7] - ldx[23,15] + ldw[22,11] - ldetWa[22,11]
   + lfrMCa[22,11] >= -2
 lfrMC3a: - ldetWa[22,15] + lfrMCa[22,15] <= 0
 lfrMC3a: - ldx[23,3] + lfrMCa[22,15] <= 0
 lfrMC3a: - ldx[23,15] + lfrMCa[22,15] <= 0
 lfrMC3a: ldw[22,15] + lfrMCa[22,15] <= 1
 lfrMC3a: - ldx[23,3] - ldx[23,15] + ldw[22,15] - ldetWa[22,15]
   + lfrMCa[22,15] >= -2
 rfLfr: ldz[23,0] + ldz[23,1] + ldz[23,2] + ldz[23,3] + ldz[23,4]
   + ldz[23,5] + ldz[23,6] + ldz[23,7] + ldz[23,8] + ldz[23,9] + ldz[23,10]
   + ldz[23,11] + ldz[23,12] + ldz[23,13] + ldz[23,14] + ldz[23,15]
   - lfrSCa[19,0] - lfrSCa[19,1] - lfrSCa[19,2] - lfrSCa[19,3]
   - lfrSCa[19,4] - lfrSCa[19,5] - lfrSCa[19,6] - lfrSCa[19,7]
   - lfrSCa[19,8] - lfrSCa[19,9] - lfrSCa[19,10] - lfrSCa[19,11]
   - lfrSCa[19,12] - lfrSCa[19,13] - lfrSCa[19,14] - lfrSCa[19,15]
   - lfrSCa[20,0] - lfrSCa[20,1] - lfrSCa[20,2] - lfrSCa[20,3]
   - lfrSCa[20,4] - lfrSCa[20,5] - lfrSCa[20,6] - lfrSCa[20,7]
   - lfrSCa[20,8] - lfrSCa[20,9] - lfrSCa[20,10] - lfrSCa[20,11]
   - lfrSCa[20,12] - lfrSCa[20,13] - lfrSCa[20,14] - lfrSCa[20,15]
   - lfrSCa[21,0] - lfrSCa[21,1] - lfrSCa[21,2] - lfrSCa[21,3]
   - lfrSCa[21,4] - lfrSCa[21,5] - lfrSCa[21,6] - lfrSCa[21,7]
   - lfrSCa[21,8] - lfrSCa[21,9] - lfrSCa[21,10] - lfrSCa[21,11]
   - lfrSCa[21,12] - lfrSCa[21,13] - lfrSCa[21,14] - lfrSCa[21,15]
   - lfrSCa[22,0] - lfrSCa[22,1] - lfrSCa[22,2] - lfrSCa[22,3]
   - lfrSCa[22,4] - lfrSCa[22,5] - lfrSCa[22,6] - lfrSCa[22,7]
   - lfrSCa[22,8] - lfrSCa[22,9] - lfrSCa[22,10] - lfrSCa[22,11]
   - lfrSCa[22,12] - lfrSCa[22,13] - lfrSCa[22,14] - lfrSCa[22,15]
   - lfrSCa[23,0] - lfrSCa[23,1] - lfrSCa[23,2] - lfrSCa[23,3]
   - lfrSCa[23,4] - lfrSCa[23,5] - lfrSCa[23,6] - lfrSCa[23,7]
   - lfrSCa[23,8] - lfrSCa[23,9] - lfrSCa[23,10] - lfrSCa[23,11]
   - lfrSCa[23,12] - lfrSCa[23,13] - lfrSCa[23,14] - lfrSCa[23,15]
   - lfrMCa[19,0] - lfrMCa[19,1] - lfrMCa[19,2] - lfrMCa[19,3]
   - lfrMCa[19,4] - lfrMCa[19,5] - lfrMCa[19,6] - lfrMCa[19,7]
   - lfrMCa[19,8] - lfrMCa[19,9] - lfrMCa[19,10] - lfrMCa[19,11]
   - lfrMCa[19,12] - lfrMCa[19,13] - lfrMCa[19,14] - lfrMCa[19,15]
   - lfrMCa[20,0] - lfrMCa[20,1] - lfrMCa[20,2] - lfrMCa[20,3]
   - lfrMCa[20,4] - lfrMCa[20,5] - lfrMCa[20,6] - lfrMCa[20,7]
   - lfrMCa[20,8] - lfrMCa[20,9] - lfrMCa[20,10] - lfrMCa[20,11]
   - lfrMCa[20,12] - lfrMCa[20,13] - lfrMCa[20,14] - lfrMCa[20,15]
   - lfrMCa[21,0] - lfrMCa[21,1] - lfrMCa[21,2] - lfrMCa[21,3]
   - lfrMCa[21,4] - lfrMCa[21,5] - lfrMCa[21,6] - lfrMCa[21,7]
   - lfrMCa[21,8] - lfrMCa[21,9] - lfrMCa[21,10] - lfrMCa[21,11]
   - lfrMCa[21,12] - lfrMCa[21,13] - lfrMCa[21,14] - lfrMCa[21,15]
   - lfrMCa[22,0] - lfrMCa[22,1] - lfrMCa[22,2] - lfrMCa[22,3]
   - lfrMCa[22,4] - lfrMCa[22,5] - lfrMCa[22,6] - lfrMCa[22,7]
   - lfrMCa[22,8] - lfrMCa[22,9] - lfrMCa[22,10] - lfrMCa[22,11]
   - lfrMCa[22,12] - lfrMCa[22,13] - lfrMCa[22,14] - lfrMCa[22,15] = 0
 KeyBgeI: - ugstka[0,0] - ugstka[1,8] - ugstka[2,2] - ugstka[3,10]
   - lgstka[46,1] - lgstka[47,9] - lgstka[48,0] - lgstka[49,8]
   - lgstka[50,2] + iLANE[0] = 0
 KeyBgeI: - ugstka[0,1] - ugstka[1,9] - ugstka[2,0] - ugstka[3,8]
   - lgstka[46,7] - lgstka[47,15] - lgstka[48,1] - lgstka[49,9]
   - lgstka[50,0] + iLANE[1] = 0
 KeyBgeI: - ugstka[0,2] - ugstka[1,10] - ugstka[2,4] - ugstka[3,12]
   - lgstka[46,0] - lgstka[47,8] - lgstka[48,2] - lgstka[49,10]
   - lgstka[50,4] + iLANE[2] = 0
 KeyBgeI: - ugstka[0,3] - ugstka[1,11] - ugstka[2,7] - ugstka[3,15]
   - lgstka[46,5] - lgstka[47,13] - lgstka[48,3] - lgstka[49,11]
   - lgstka[50,7] + iLANE[3] = 0
 KeyBgeI: - ugstka[0,4] - ugstka[1,12] - ugstka[2,6] - ugstka[3,14]
   - lgstka[46,2] - lgstka[47,10] - lgstka[48,4] - lgstka[49,12]
   - lgstka[50,6] + iLANE[4] = 0
 KeyBgeI: - ugstka[0,5] - ugstka[1,13] - ugstka[2,3] - ugstka[3,11]
   - lgstka[46,6] - lgstka[47,14] - lgstka[48,5] - lgstka[49,13]
   - lgstka[50,3] + iLANE[5] = 0
 KeyBgeI: - ugstka[0,6] - ugstka[1,14] - ugstka[2,5] - ugstka[3,13]
   - lgstka[46,4] - lgstka[47,12] - lgstka[48,6] - lgstka[49,14]
   - lgstka[50,5] + iLANE[6] = 0
 KeyBgeI: - ugstka[0,7] - ugstka[1,15] - ugstka[2,1] - ugstka[3,9]
   - lgstka[46,3] - lgstka[47,11] - lgstka[48,7] - lgstka[49,15]
   - lgstka[50,1] + iLANE[7] = 0
 KeyBgeI: - ugstka[0,8] - ugstka[1,2] - ugstka[2,10] - ugstka[3,4]
   - lgstka[46,9] - lgstka[47,0] - lgstka[48,8] - lgstka[49,2]
   - lgstka[50,10] + iLANE[8] = 0
 KeyBgeI: - ugstka[0,9] - ugstka[1,0] - ugstka[2,8] - ugstka[3,2]
   - lgstka[46,15] - lgstka[47,1] - lgstka[48,9] - lgstka[49,0]
   - lgstka[50,8] + iLANE[9] = 0
 KeyBgeI: - ugstka[0,10] - ugstka[1,4] - ugstka[2,12] - ugstka[3,6]
   - lgstka[46,8] - lgstka[47,2] - lgstka[48,10] - lgstka[49,4]
   - lgstka[50,12] + iLANE[10] = 0
 KeyBgeI: - ugstka[0,11] - ugstka[1,7] - ugstka[2,15] - ugstka[3,1]
   - lgstka[46,13] - lgstka[47,3] - lgstka[48,11] - lgstka[49,7]
   - lgstka[50,15] + iLANE[11] = 0
 KeyBgeI: - ugstka[0,12] - ugstka[1,6] - ugstka[2,14] - ugstka[3,5]
   - lgstka[46,10] - lgstka[47,4] - lgstka[48,12] - lgstka[49,6]
   - lgstka[50,14] + iLANE[12] = 0
 KeyBgeI: - ugstka[0,13] - ugstka[1,3] - ugstka[2,11] - ugstka[3,7]
   - lgstka[46,14] - lgstka[47,5] - lgstka[48,13] - lgstka[49,3]
   - lgstka[50,11] + iLANE[13] = 0
 KeyBgeI: - ugstka[0,14] - ugstka[1,5] - ugstka[2,13] - ugstka[3,3]
   - lgstka[46,12] - lgstka[47,6] - lgstka[48,14] - lgstka[49,5]
   - lgstka[50,13] + iLANE[14] = 0
 KeyBgeI: - ugstka[0,15] - ugstka[1,1] - ugstka[2,9] - ugstka[3,0]
   - lgstka[46,11] - lgstka[47,7] - lgstka[48,15] - lgstka[49,1]
   - lgstka[50,9] + iLANE[15] = 0
 KeyBgeI: iLANE[0] - 1000 ikz[0] <= 2
 KeyBgeI: iLANE[0] - 1000 ikz[0] >= -997
 KeyBgeI: - iLANE[0] + ikSUM[0] - 1000 ikz[0] <= 0
 KeyBgeI: - iLANE[0] + ikSUM[0] + 1000 ikz[0] >= 0
 KeyBgeI: ikSUM[0] + 1000 ikz[0] <= 1002
 KeyBgeI: ikSUM[0] - 1000 ikz[0] >= -998
 KeyBgeI: iLANE[1] - 1000 ikz[1] <= 2
 KeyBgeI: iLANE[1] - 1000 ikz[1] >= -997
 KeyBgeI: - iLANE[1] + ikSUM[1] - 1000 ikz[1] <= 0
 KeyBgeI: - iLANE[1] + ikSUM[1] + 1000 ikz[1] >= 0
 KeyBgeI: ikSUM[1] + 1000 ikz[1] <= 1002
 KeyBgeI: ikSUM[1] - 1000 ikz[1] >= -998
 KeyBgeI: iLANE[2] - 1000 ikz[2] <= 2
 KeyBgeI: iLANE[2] - 1000 ikz[2] >= -997
 KeyBgeI: - iLANE[2] + ikSUM[2] - 1000 ikz[2] <= 0
 KeyBgeI: - iLANE[2] + ikSUM[2] + 1000 ikz[2] >= 0
 KeyBgeI: ikSUM[2] + 1000 ikz[2] <= 1002
 KeyBgeI: ikSUM[2] - 1000 ikz[2] >= -998
 KeyBgeI: iLANE[3] - 1000 ikz[3] <= 2
 KeyBgeI: iLANE[3] - 1000 ikz[3] >= -997
 KeyBgeI: - iLANE[3] + ikSUM[3] - 1000 ikz[3] <= 0
 KeyBgeI: - iLANE[3] + ikSUM[3] + 1000 ikz[3] >= 0
 KeyBgeI: ikSUM[3] + 1000 ikz[3] <= 1002
 KeyBgeI: ikSUM[3] - 1000 ikz[3] >= -998
 KeyBgeI: iLANE[4] - 1000 ikz[4] <= 2
 KeyBgeI: iLANE[4] - 1000 ikz[4] >= -997
 KeyBgeI: - iLANE[4] + ikSUM[4] - 1000 ikz[4] <= 0
 KeyBgeI: - iLANE[4] + ikSUM[4] + 1000 ikz[4] >= 0
 KeyBgeI: ikSUM[4] + 1000 ikz[4] <= 1002
 KeyBgeI: ikSUM[4] - 1000 ikz[4] >= -998
 KeyBgeI: iLANE[5] - 1000 ikz[5] <= 2
 KeyBgeI: iLANE[5] - 1000 ikz[5] >= -997
 KeyBgeI: - iLANE[5] + ikSUM[5] - 1000 ikz[5] <= 0
 KeyBgeI: - iLANE[5] + ikSUM[5] + 1000 ikz[5] >= 0
 KeyBgeI: ikSUM[5] + 1000 ikz[5] <= 1002
 KeyBgeI: ikSUM[5] - 1000 ikz[5] >= -998
 KeyBgeI: iLANE[6] - 1000 ikz[6] <= 2
 KeyBgeI: iLANE[6] - 1000 ikz[6] >= -997
 KeyBgeI: - iLANE[6] + ikSUM[6] - 1000 ikz[6] <= 0
 KeyBgeI: - iLANE[6] + ikSUM[6] + 1000 ikz[6] >= 0
 KeyBgeI: ikSUM[6] + 1000 ikz[6] <= 1002
 KeyBgeI: ikSUM[6] - 1000 ikz[6] >= -998
 KeyBgeI: iLANE[7] - 1000 ikz[7] <= 2
 KeyBgeI: iLANE[7] - 1000 ikz[7] >= -997
 KeyBgeI: - iLANE[7] + ikSUM[7] - 1000 ikz[7] <= 0
 KeyBgeI: - iLANE[7] + ikSUM[7] + 1000 ikz[7] >= 0
 KeyBgeI: ikSUM[7] + 1000 ikz[7] <= 1002
 KeyBgeI: ikSUM[7] - 1000 ikz[7] >= -998
 KeyBgeI: iLANE[8] - 1000 ikz[8] <= 2
 KeyBgeI: iLANE[8] - 1000 ikz[8] >= -997
 KeyBgeI: - iLANE[8] + ikSUM[8] - 1000 ikz[8] <= 0
 KeyBgeI: - iLANE[8] + ikSUM[8] + 1000 ikz[8] >= 0
 KeyBgeI: ikSUM[8] + 1000 ikz[8] <= 1002
 KeyBgeI: ikSUM[8] - 1000 ikz[8] >= -998
 KeyBgeI: iLANE[9] - 1000 ikz[9] <= 2
 KeyBgeI: iLANE[9] - 1000 ikz[9] >= -997
 KeyBgeI: - iLANE[9] + ikSUM[9] - 1000 ikz[9] <= 0
 KeyBgeI: - iLANE[9] + ikSUM[9] + 1000 ikz[9] >= 0
 KeyBgeI: ikSUM[9] + 1000 ikz[9] <= 1002
 KeyBgeI: ikSUM[9] - 1000 ikz[9] >= -998
 KeyBgeI: iLANE[10] - 1000 ikz[10] <= 2
 KeyBgeI: iLANE[10] - 1000 ikz[10] >= -997
 KeyBgeI: - iLANE[10] + ikSUM[10] - 1000 ikz[10] <= 0
 KeyBgeI: - iLANE[10] + ikSUM[10] + 1000 ikz[10] >= 0
 KeyBgeI: ikSUM[10] + 1000 ikz[10] <= 1002
 KeyBgeI: ikSUM[10] - 1000 ikz[10] >= -998
 KeyBgeI: iLANE[11] - 1000 ikz[11] <= 2
 KeyBgeI: iLANE[11] - 1000 ikz[11] >= -997
 KeyBgeI: - iLANE[11] + ikSUM[11] - 1000 ikz[11] <= 0
 KeyBgeI: - iLANE[11] + ikSUM[11] + 1000 ikz[11] >= 0
 KeyBgeI: ikSUM[11] + 1000 ikz[11] <= 1002
 KeyBgeI: ikSUM[11] - 1000 ikz[11] >= -998
 KeyBgeI: iLANE[12] - 1000 ikz[12] <= 2
 KeyBgeI: iLANE[12] - 1000 ikz[12] >= -997
 KeyBgeI: - iLANE[12] + ikSUM[12] - 1000 ikz[12] <= 0
 KeyBgeI: - iLANE[12] + ikSUM[12] + 1000 ikz[12] >= 0
 KeyBgeI: ikSUM[12] + 1000 ikz[12] <= 1002
 KeyBgeI: ikSUM[12] - 1000 ikz[12] >= -998
 KeyBgeI: iLANE[13] - 1000 ikz[13] <= 2
 KeyBgeI: iLANE[13] - 1000 ikz[13] >= -997
 KeyBgeI: - iLANE[13] + ikSUM[13] - 1000 ikz[13] <= 0
 KeyBgeI: - iLANE[13] + ikSUM[13] + 1000 ikz[13] >= 0
 KeyBgeI: ikSUM[13] + 1000 ikz[13] <= 1002
 KeyBgeI: ikSUM[13] - 1000 ikz[13] >= -998
 KeyBgeI: iLANE[14] - 1000 ikz[14] <= 2
 KeyBgeI: iLANE[14] - 1000 ikz[14] >= -997
 KeyBgeI: - iLANE[14] + ikSUM[14] - 1000 ikz[14] <= 0
 KeyBgeI: - iLANE[14] + ikSUM[14] + 1000 ikz[14] >= 0
 KeyBgeI: ikSUM[14] + 1000 ikz[14] <= 1002
 KeyBgeI: ikSUM[14] - 1000 ikz[14] >= -998
 KeyBgeI: iLANE[15] - 1000 ikz[15] <= 2
 KeyBgeI: iLANE[15] - 1000 ikz[15] >= -997
 KeyBgeI: - iLANE[15] + ikSUM[15] - 1000 ikz[15] <= 0
 KeyBgeI: - iLANE[15] + ikSUM[15] + 1000 ikz[15] >= 0
 KeyBgeI: ikSUM[15] + 1000 ikz[15] <= 1002
 KeyBgeI: ikSUM[15] - 1000 ikz[15] >= -998
 KeyBgeG: - uGstk[0,0] - uGstk[1,8] - uGstk[2,2] - uGstk[3,10]
   - lGstk[46,1] - lGstk[47,9] - lGstk[48,0] - lGstk[49,8] - lGstk[50,2]
   + gLANE[0] = 0
 KeyBgeG: - uGstk[0,1] - uGstk[1,9] - uGstk[2,0] - uGstk[3,8] - lGstk[46,7]
   - lGstk[47,15] - lGstk[48,1] - lGstk[49,9] - lGstk[50,0] + gLANE[1] = 0
 KeyBgeG: - uGstk[0,2] - uGstk[1,10] - uGstk[2,4] - uGstk[3,12]
   - lGstk[46,0] - lGstk[47,8] - lGstk[48,2] - lGstk[49,10] - lGstk[50,4]
   + gLANE[2] = 0
 KeyBgeG: - uGstk[0,3] - uGstk[1,11] - uGstk[2,7] - uGstk[3,15]
   - lGstk[46,5] - lGstk[47,13] - lGstk[48,3] - lGstk[49,11] - lGstk[50,7]
   + gLANE[3] = 0
 KeyBgeG: - uGstk[0,4] - uGstk[1,12] - uGstk[2,6] - uGstk[3,14]
   - lGstk[46,2] - lGstk[47,10] - lGstk[48,4] - lGstk[49,12] - lGstk[50,6]
   + gLANE[4] = 0
 KeyBgeG: - uGstk[0,5] - uGstk[1,13] - uGstk[2,3] - uGstk[3,11]
   - lGstk[46,6] - lGstk[47,14] - lGstk[48,5] - lGstk[49,13] - lGstk[50,3]
   + gLANE[5] = 0
 KeyBgeG: - uGstk[0,6] - uGstk[1,14] - uGstk[2,5] - uGstk[3,13]
   - lGstk[46,4] - lGstk[47,12] - lGstk[48,6] - lGstk[49,14] - lGstk[50,5]
   + gLANE[6] = 0
 KeyBgeG: - uGstk[0,7] - uGstk[1,15] - uGstk[2,1] - uGstk[3,9]
   - lGstk[46,3] - lGstk[47,11] - lGstk[48,7] - lGstk[49,15] - lGstk[50,1]
   + gLANE[7] = 0
 KeyBgeG: - uGstk[0,8] - uGstk[1,2] - uGstk[2,10] - uGstk[3,4]
   - lGstk[46,9] - lGstk[47,0] - lGstk[48,8] - lGstk[49,2] - lGstk[50,10]
   + gLANE[8] = 0
 KeyBgeG: - uGstk[0,9] - uGstk[1,0] - uGstk[2,8] - uGstk[3,2]
   - lGstk[46,15] - lGstk[47,1] - lGstk[48,9] - lGstk[49,0] - lGstk[50,8]
   + gLANE[9] = 0
 KeyBgeG: - uGstk[0,10] - uGstk[1,4] - uGstk[2,12] - uGstk[3,6]
   - lGstk[46,8] - lGstk[47,2] - lGstk[48,10] - lGstk[49,4] - lGstk[50,12]
   + gLANE[10] = 0
 KeyBgeG: - uGstk[0,11] - uGstk[1,7] - uGstk[2,15] - uGstk[3,1]
   - lGstk[46,13] - lGstk[47,3] - lGstk[48,11] - lGstk[49,7] - lGstk[50,15]
   + gLANE[11] = 0
 KeyBgeG: - uGstk[0,12] - uGstk[1,6] - uGstk[2,14] - uGstk[3,5]
   - lGstk[46,10] - lGstk[47,4] - lGstk[48,12] - lGstk[49,6] - lGstk[50,14]
   + gLANE[12] = 0
 KeyBgeG: - uGstk[0,13] - uGstk[1,3] - uGstk[2,11] - uGstk[3,7]
   - lGstk[46,14] - lGstk[47,5] - lGstk[48,13] - lGstk[49,3] - lGstk[50,11]
   + gLANE[13] = 0
 KeyBgeG: - uGstk[0,14] - uGstk[1,5] - uGstk[2,13] - uGstk[3,3]
   - lGstk[46,12] - lGstk[47,6] - lGstk[48,14] - lGstk[49,5] - lGstk[50,13]
   + gLANE[14] = 0
 KeyBgeG: - uGstk[0,15] - uGstk[1,1] - uGstk[2,9] - uGstk[3,0]
   - lGstk[46,11] - lGstk[47,7] - lGstk[48,15] - lGstk[49,1] - lGstk[50,9]
   + gLANE[15] = 0
 KeyBgeG: gLANE[0] - 1000 gkz[0] <= 2
 KeyBgeG: gLANE[0] - 1000 gkz[0] >= -997
 KeyBgeG: - gLANE[0] + gkSUM[0] - 1000 gkz[0] <= 0
 KeyBgeG: - gLANE[0] + gkSUM[0] + 1000 gkz[0] >= 0
 KeyBgeG: gkSUM[0] + 1000 gkz[0] <= 1002
 KeyBgeG: gkSUM[0] - 1000 gkz[0] >= -998
 KeyBgeG: gLANE[1] - 1000 gkz[1] <= 2
 KeyBgeG: gLANE[1] - 1000 gkz[1] >= -997
 KeyBgeG: - gLANE[1] + gkSUM[1] - 1000 gkz[1] <= 0
 KeyBgeG: - gLANE[1] + gkSUM[1] + 1000 gkz[1] >= 0
 KeyBgeG: gkSUM[1] + 1000 gkz[1] <= 1002
 KeyBgeG: gkSUM[1] - 1000 gkz[1] >= -998
 KeyBgeG: gLANE[2] - 1000 gkz[2] <= 2
 KeyBgeG: gLANE[2] - 1000 gkz[2] >= -997
 KeyBgeG: - gLANE[2] + gkSUM[2] - 1000 gkz[2] <= 0
 KeyBgeG: - gLANE[2] + gkSUM[2] + 1000 gkz[2] >= 0
 KeyBgeG: gkSUM[2] + 1000 gkz[2] <= 1002
 KeyBgeG: gkSUM[2] - 1000 gkz[2] >= -998
 KeyBgeG: gLANE[3] - 1000 gkz[3] <= 2
 KeyBgeG: gLANE[3] - 1000 gkz[3] >= -997
 KeyBgeG: - gLANE[3] + gkSUM[3] - 1000 gkz[3] <= 0
 KeyBgeG: - gLANE[3] + gkSUM[3] + 1000 gkz[3] >= 0
 KeyBgeG: gkSUM[3] + 1000 gkz[3] <= 1002
 KeyBgeG: gkSUM[3] - 1000 gkz[3] >= -998
 KeyBgeG: gLANE[4] - 1000 gkz[4] <= 2
 KeyBgeG: gLANE[4] - 1000 gkz[4] >= -997
 KeyBgeG: - gLANE[4] + gkSUM[4] - 1000 gkz[4] <= 0
 KeyBgeG: - gLANE[4] + gkSUM[4] + 1000 gkz[4] >= 0
 KeyBgeG: gkSUM[4] + 1000 gkz[4] <= 1002
 KeyBgeG: gkSUM[4] - 1000 gkz[4] >= -998
 KeyBgeG: gLANE[5] - 1000 gkz[5] <= 2
 KeyBgeG: gLANE[5] - 1000 gkz[5] >= -997
 KeyBgeG: - gLANE[5] + gkSUM[5] - 1000 gkz[5] <= 0
 KeyBgeG: - gLANE[5] + gkSUM[5] + 1000 gkz[5] >= 0
 KeyBgeG: gkSUM[5] + 1000 gkz[5] <= 1002
 KeyBgeG: gkSUM[5] - 1000 gkz[5] >= -998
 KeyBgeG: gLANE[6] - 1000 gkz[6] <= 2
 KeyBgeG: gLANE[6] - 1000 gkz[6] >= -997
 KeyBgeG: - gLANE[6] + gkSUM[6] - 1000 gkz[6] <= 0
 KeyBgeG: - gLANE[6] + gkSUM[6] + 1000 gkz[6] >= 0
 KeyBgeG: gkSUM[6] + 1000 gkz[6] <= 1002
 KeyBgeG: gkSUM[6] - 1000 gkz[6] >= -998
 KeyBgeG: gLANE[7] - 1000 gkz[7] <= 2
 KeyBgeG: gLANE[7] - 1000 gkz[7] >= -997
 KeyBgeG: - gLANE[7] + gkSUM[7] - 1000 gkz[7] <= 0
 KeyBgeG: - gLANE[7] + gkSUM[7] + 1000 gkz[7] >= 0
 KeyBgeG: gkSUM[7] + 1000 gkz[7] <= 1002
 KeyBgeG: gkSUM[7] - 1000 gkz[7] >= -998
 KeyBgeG: gLANE[8] - 1000 gkz[8] <= 2
 KeyBgeG: gLANE[8] - 1000 gkz[8] >= -997
 KeyBgeG: - gLANE[8] + gkSUM[8] - 1000 gkz[8] <= 0
 KeyBgeG: - gLANE[8] + gkSUM[8] + 1000 gkz[8] >= 0
 KeyBgeG: gkSUM[8] + 1000 gkz[8] <= 1002
 KeyBgeG: gkSUM[8] - 1000 gkz[8] >= -998
 KeyBgeG: gLANE[9] - 1000 gkz[9] <= 2
 KeyBgeG: gLANE[9] - 1000 gkz[9] >= -997
 KeyBgeG: - gLANE[9] + gkSUM[9] - 1000 gkz[9] <= 0
 KeyBgeG: - gLANE[9] + gkSUM[9] + 1000 gkz[9] >= 0
 KeyBgeG: gkSUM[9] + 1000 gkz[9] <= 1002
 KeyBgeG: gkSUM[9] - 1000 gkz[9] >= -998
 KeyBgeG: gLANE[10] - 1000 gkz[10] <= 2
 KeyBgeG: gLANE[10] - 1000 gkz[10] >= -997
 KeyBgeG: - gLANE[10] + gkSUM[10] - 1000 gkz[10] <= 0
 KeyBgeG: - gLANE[10] + gkSUM[10] + 1000 gkz[10] >= 0
 KeyBgeG: gkSUM[10] + 1000 gkz[10] <= 1002
 KeyBgeG: gkSUM[10] - 1000 gkz[10] >= -998
 KeyBgeG: gLANE[11] - 1000 gkz[11] <= 2
 KeyBgeG: gLANE[11] - 1000 gkz[11] >= -997
 KeyBgeG: - gLANE[11] + gkSUM[11] - 1000 gkz[11] <= 0
 KeyBgeG: - gLANE[11] + gkSUM[11] + 1000 gkz[11] >= 0
 KeyBgeG: gkSUM[11] + 1000 gkz[11] <= 1002
 KeyBgeG: gkSUM[11] - 1000 gkz[11] >= -998
 KeyBgeG: gLANE[12] - 1000 gkz[12] <= 2
 KeyBgeG: gLANE[12] - 1000 gkz[12] >= -997
 KeyBgeG: - gLANE[12] + gkSUM[12] - 1000 gkz[12] <= 0
 KeyBgeG: - gLANE[12] + gkSUM[12] + 1000 gkz[12] >= 0
 KeyBgeG: gkSUM[12] + 1000 gkz[12] <= 1002
 KeyBgeG: gkSUM[12] - 1000 gkz[12] >= -998
 KeyBgeG: gLANE[13] - 1000 gkz[13] <= 2
 KeyBgeG: gLANE[13] - 1000 gkz[13] >= -997
 KeyBgeG: - gLANE[13] + gkSUM[13] - 1000 gkz[13] <= 0
 KeyBgeG: - gLANE[13] + gkSUM[13] + 1000 gkz[13] >= 0
 KeyBgeG: gkSUM[13] + 1000 gkz[13] <= 1002
 KeyBgeG: gkSUM[13] - 1000 gkz[13] >= -998
 KeyBgeG: gLANE[14] - 1000 gkz[14] <= 2
 KeyBgeG: gLANE[14] - 1000 gkz[14] >= -997
 KeyBgeG: - gLANE[14] + gkSUM[14] - 1000 gkz[14] <= 0
 KeyBgeG: - gLANE[14] + gkSUM[14] + 1000 gkz[14] >= 0
 KeyBgeG: gkSUM[14] + 1000 gkz[14] <= 1002
 KeyBgeG: gkSUM[14] - 1000 gkz[14] >= -998
 KeyBgeG: gLANE[15] - 1000 gkz[15] <= 2
 KeyBgeG: gLANE[15] - 1000 gkz[15] >= -997
 KeyBgeG: - gLANE[15] + gkSUM[15] - 1000 gkz[15] <= 0
 KeyBgeG: - gLANE[15] + gkSUM[15] + 1000 gkz[15] >= 0
 KeyBgeG: gkSUM[15] + 1000 gkz[15] <= 1002
 KeyBgeG: gkSUM[15] - 1000 gkz[15] >= -998
 IKS: - 8 ikSUM[0] - 8 ikSUM[1] - 8 ikSUM[2] - 8 ikSUM[3] - 8 ikSUM[4]
   - 8 ikSUM[5] - 8 ikSUM[6] - 8 ikSUM[7] - 8 ikSUM[8] - 8 ikSUM[9]
   - 8 ikSUM[10] - 8 ikSUM[11] - 8 ikSUM[12] - 8 ikSUM[13] - 8 ikSUM[14]
   - 8 ikSUM[15] + iks = 0
 GKS: - 8 gkSUM[0] - 8 gkSUM[1] - 8 gkSUM[2] - 8 gkSUM[3] - 8 gkSUM[4]
   - 8 gkSUM[5] - 8 gkSUM[6] - 8 gkSUM[7] - 8 gkSUM[8] - 8 gkSUM[9]
   - 8 gkSUM[10] - 8 gkSUM[11] - 8 gkSUM[12] - 8 gkSUM[13] - 8 gkSUM[14]
   - 8 gkSUM[15] + gks = 0
 rb: - 8 udw[0,0] - 8 udw[0,1] - 8 udw[0,2] - 8 udw[0,3] - 8 udw[0,4]
   - 8 udw[0,5] - 8 udw[0,6] - 8 udw[0,7] - 8 udw[0,8] - 8 udw[0,9]
   - 8 udw[0,10] - 8 udw[0,11] - 8 udw[0,12] - 8 udw[0,13] - 8 udw[0,14]
   - 8 udw[0,15] + rb = 0
 rf: - 8 ldz[23,0] - 8 ldz[23,1] - 8 ldz[23,2] - 8 ldz[23,3] - 8 ldz[23,4]
   - 8 ldz[23,5] - 8 ldz[23,6] - 8 ldz[23,7] - 8 ldz[23,8] - 8 ldz[23,9]
   - 8 ldz[23,10] - 8 ldz[23,11] - 8 ldz[23,12] - 8 ldz[23,13]
   - 8 ldz[23,14] - 8 ldz[23,15] + rf = 0
 cb: 8 udy[4,0] + 8 udy[4,1] + 8 udy[4,2] + 8 udy[4,3] + 8 udy[4,4]
   + 8 udy[4,5] + 8 udy[4,6] + 8 udy[4,7] + 8 udy[4,8] + 8 udy[4,9]
   + 8 udy[4,10] + 8 udy[4,11] + 8 udy[4,12] + 8 udy[4,13] + 8 udy[4,14]
   + 8 udy[4,15] - rb + cb = 0
 cf: 8 ldx[19,0] + 8 ldx[19,1] + 8 ldx[19,2] + 8 ldx[19,3] + 8 ldx[19,4]
   + 8 ldx[19,5] + 8 ldx[19,6] + 8 ldx[19,7] + 8 ldx[19,8] + 8 ldx[19,9]
   + 8 ldx[19,10] + 8 ldx[19,11] + 8 ldx[19,12] + 8 ldx[19,13]
   + 8 ldx[19,14] + 8 ldx[19,15] - rf + cf = 0
 cbp: - 8 uFrSC[1,0] - 8 uFrSC[1,1] - 8 uFrSC[1,2] - 8 uFrSC[1,3]
   - 8 uFrSC[1,4] - 8 uFrSC[1,5] - 8 uFrSC[1,6] - 8 uFrSC[1,7]
   - 8 uFrSC[1,8] - 8 uFrSC[1,9] - 8 uFrSC[1,10] - 8 uFrSC[1,11]
   - 8 uFrSC[1,12] - 8 uFrSC[1,13] - 8 uFrSC[1,14] - 8 uFrSC[1,15]
   - 8 uFrSC[2,0] - 8 uFrSC[2,1] - 8 uFrSC[2,2] - 8 uFrSC[2,3]
   - 8 uFrSC[2,4] - 8 uFrSC[2,5] - 8 uFrSC[2,6] - 8 uFrSC[2,7]
   - 8 uFrSC[2,8] - 8 uFrSC[2,9] - 8 uFrSC[2,10] - 8 uFrSC[2,11]
   - 8 uFrSC[2,12] - 8 uFrSC[2,13] - 8 uFrSC[2,14] - 8 uFrSC[2,15]
   - 8 uFrSC[3,0] - 8 uFrSC[3,1] - 8 uFrSC[3,2] - 8 uFrSC[3,3]
   - 8 uFrSC[3,4] - 8 uFrSC[3,5] - 8 uFrSC[3,6] - 8 uFrSC[3,7]
   - 8 uFrSC[3,8] - 8 uFrSC[3,9] - 8 uFrSC[3,10] - 8 uFrSC[3,11]
   - 8 uFrSC[3,12] - 8 uFrSC[3,13] - 8 uFrSC[3,14] - 8 uFrSC[3,15]
   - 8 uFrSC[4,0] - 8 uFrSC[4,1] - 8 uFrSC[4,2] - 8 uFrSC[4,3]
   - 8 uFrSC[4,4] - 8 uFrSC[4,5] - 8 uFrSC[4,6] - 8 uFrSC[4,7]
   - 8 uFrSC[4,8] - 8 uFrSC[4,9] - 8 uFrSC[4,10] - 8 uFrSC[4,11]
   - 8 uFrSC[4,12] - 8 uFrSC[4,13] - 8 uFrSC[4,14] - 8 uFrSC[4,15]
   - 8 uFrMC[2,0] - 8 uFrMC[2,1] - 8 uFrMC[2,2] - 8 uFrMC[2,3]
   - 8 uFrMC[2,4] - 8 uFrMC[2,5] - 8 uFrMC[2,6] - 8 uFrMC[2,7]
   - 8 uFrMC[2,8] - 8 uFrMC[2,9] - 8 uFrMC[2,10] - 8 uFrMC[2,11]
   - 8 uFrMC[2,12] - 8 uFrMC[2,13] - 8 uFrMC[2,14] - 8 uFrMC[2,15]
   - 8 uFrMC[3,0] - 8 uFrMC[3,1] - 8 uFrMC[3,2] - 8 uFrMC[3,3]
   - 8 uFrMC[3,4] - 8 uFrMC[3,5] - 8 uFrMC[3,6] - 8 uFrMC[3,7]
   - 8 uFrMC[3,8] - 8 uFrMC[3,9] - 8 uFrMC[3,10] - 8 uFrMC[3,11]
   - 8 uFrMC[3,12] - 8 uFrMC[3,13] - 8 uFrMC[3,14] - 8 uFrMC[3,15]
   - 8 uFrMC[4,0] - 8 uFrMC[4,1] - 8 uFrMC[4,2] - 8 uFrMC[4,3]
   - 8 uFrMC[4,4] - 8 uFrMC[4,5] - 8 uFrMC[4,6] - 8 uFrMC[4,7]
   - 8 uFrMC[4,8] - 8 uFrMC[4,9] - 8 uFrMC[4,10] - 8 uFrMC[4,11]
   - 8 uFrMC[4,12] - 8 uFrMC[4,13] - 8 uFrMC[4,14] - 8 uFrMC[4,15] + cbp
   = 0
 cfp: - 8 lFrSC[19,0] - 8 lFrSC[19,1] - 8 lFrSC[19,2] - 8 lFrSC[19,3]
   - 8 lFrSC[19,4] - 8 lFrSC[19,5] - 8 lFrSC[19,6] - 8 lFrSC[19,7]
   - 8 lFrSC[19,8] - 8 lFrSC[19,9] - 8 lFrSC[19,10] - 8 lFrSC[19,11]
   - 8 lFrSC[19,12] - 8 lFrSC[19,13] - 8 lFrSC[19,14] - 8 lFrSC[19,15]
   - 8 lFrSC[20,0] - 8 lFrSC[20,1] - 8 lFrSC[20,2] - 8 lFrSC[20,3]
   - 8 lFrSC[20,4] - 8 lFrSC[20,5] - 8 lFrSC[20,6] - 8 lFrSC[20,7]
   - 8 lFrSC[20,8] - 8 lFrSC[20,9] - 8 lFrSC[20,10] - 8 lFrSC[20,11]
   - 8 lFrSC[20,12] - 8 lFrSC[20,13] - 8 lFrSC[20,14] - 8 lFrSC[20,15]
   - 8 lFrSC[21,0] - 8 lFrSC[21,1] - 8 lFrSC[21,2] - 8 lFrSC[21,3]
   - 8 lFrSC[21,4] - 8 lFrSC[21,5] - 8 lFrSC[21,6] - 8 lFrSC[21,7]
   - 8 lFrSC[21,8] - 8 lFrSC[21,9] - 8 lFrSC[21,10] - 8 lFrSC[21,11]
   - 8 lFrSC[21,12] - 8 lFrSC[21,13] - 8 lFrSC[21,14] - 8 lFrSC[21,15]
   - 8 lFrSC[22,0] - 8 lFrSC[22,1] - 8 lFrSC[22,2] - 8 lFrSC[22,3]
   - 8 lFrSC[22,4] - 8 lFrSC[22,5] - 8 lFrSC[22,6] - 8 lFrSC[22,7]
   - 8 lFrSC[22,8] - 8 lFrSC[22,9] - 8 lFrSC[22,10] - 8 lFrSC[22,11]
   - 8 lFrSC[22,12] - 8 lFrSC[22,13] - 8 lFrSC[22,14] - 8 lFrSC[22,15]
   - 8 lFrSC[23,0] - 8 lFrSC[23,1] - 8 lFrSC[23,2] - 8 lFrSC[23,3]
   - 8 lFrSC[23,4] - 8 lFrSC[23,5] - 8 lFrSC[23,6] - 8 lFrSC[23,7]
   - 8 lFrSC[23,8] - 8 lFrSC[23,9] - 8 lFrSC[23,10] - 8 lFrSC[23,11]
   - 8 lFrSC[23,12] - 8 lFrSC[23,13] - 8 lFrSC[23,14] - 8 lFrSC[23,15]
   - 8 lFrMC[19,0] - 8 lFrMC[19,1] - 8 lFrMC[19,2] - 8 lFrMC[19,3]
   - 8 lFrMC[19,4] - 8 lFrMC[19,5] - 8 lFrMC[19,6] - 8 lFrMC[19,7]
   - 8 lFrMC[19,8] - 8 lFrMC[19,9] - 8 lFrMC[19,10] - 8 lFrMC[19,11]
   - 8 lFrMC[19,12] - 8 lFrMC[19,13] - 8 lFrMC[19,14] - 8 lFrMC[19,15]
   - 8 lFrMC[20,0] - 8 lFrMC[20,1] - 8 lFrMC[20,2] - 8 lFrMC[20,3]
   - 8 lFrMC[20,4] - 8 lFrMC[20,5] - 8 lFrMC[20,6] - 8 lFrMC[20,7]
   - 8 lFrMC[20,8] - 8 lFrMC[20,9] - 8 lFrMC[20,10] - 8 lFrMC[20,11]
   - 8 lFrMC[20,12] - 8 lFrMC[20,13] - 8 lFrMC[20,14] - 8 lFrMC[20,15]
   - 8 lFrMC[21,0] - 8 lFrMC[21,1] - 8 lFrMC[21,2] - 8 lFrMC[21,3]
   - 8 lFrMC[21,4] - 8 lFrMC[21,5] - 8 lFrMC[21,6] - 8 lFrMC[21,7]
   - 8 lFrMC[21,8] - 8 lFrMC[21,9] - 8 lFrMC[21,10] - 8 lFrMC[21,11]
   - 8 lFrMC[21,12] - 8 lFrMC[21,13] - 8 lFrMC[21,14] - 8 lFrMC[21,15]
   - 8 lFrMC[22,0] - 8 lFrMC[22,1] - 8 lFrMC[22,2] - 8 lFrMC[22,3]
   - 8 lFrMC[22,4] - 8 lFrMC[22,5] - 8 lFrMC[22,6] - 8 lFrMC[22,7]
   - 8 lFrMC[22,8] - 8 lFrMC[22,9] - 8 lFrMC[22,10] - 8 lFrMC[22,11]
   - 8 lFrMC[22,12] - 8 lFrMC[22,13] - 8 lFrMC[22,14] - 8 lFrMC[22,15]
   + cfp = 0
 D: D = 130.55
 Dc: - D + Dc = 2
 Qc: - 2 cb - 2 cf + Qc = 5.1
 Mc0: - Dc + Mc >= 2
 Mc1: gks + Mc - T2 >= 0
 Mc2: 2 cbp + 2 cfp - Qc + Mc >= 0
 Mc3: - iks + Mc >= 0
 T0: - Dc + T0 = 0
 T1: - gks - Dc + T1 = 2
 T2p: - rb - gks + cbp - Dc + T2 = 1
 T31: - gks + 2 cbp + 2 cfp - Qc + T31 = 0
 T32: - iks + T32 = 5.1
 T4: T4 = 208
 Tc: - T0 + Tc >= 0
 Tc: - T1 + Tc >= 0
 Tc: - T2 + Tc >= 0
 Tc: - T31 + Tc >= 0
 Tc: - T32 + Tc >= 0
 Tc: - T4 + Tc >= 0
 R11540: Tc <= 255
 iBCT: [ uDX[12,0] * lDY[12,0] - uDX[12,0] * ldy[12,0]
   + uDX[12,1] * lDY[12,1] - uDX[12,1] * ldy[12,1] + uDX[12,2] * lDY[12,2]
   - uDX[12,2] * ldy[12,2] + uDX[12,3] * lDY[12,3] - uDX[12,3] * ldy[12,3]
   + uDX[12,4] * lDY[12,4] - uDX[12,4] * ldy[12,4] + uDX[12,5] * lDY[12,5]
   - uDX[12,5] * ldy[12,5] + uDX[12,6] * lDY[12,6] - uDX[12,6] * ldy[12,6]
   + uDX[12,7] * lDY[12,7] - uDX[12,7] * ldy[12,7] + uDX[12,8] * lDY[12,8]
   - uDX[12,8] * ldy[12,8] + uDX[12,9] * lDY[12,9] - uDX[12,9] * ldy[12,9]
   + uDX[12,10] * lDY[12,10] - uDX[12,10] * ldy[12,10]
   + uDX[12,11] * lDY[12,11] - uDX[12,11] * ldy[12,11]
   + uDX[12,12] * lDY[12,12] - uDX[12,12] * ldy[12,12]
   + uDX[12,13] * lDY[12,13] - uDX[12,13] * ldy[12,13]
   + uDX[12,14] * lDY[12,14] - uDX[12,14] * ldy[12,14]
   + uDX[12,15] * lDY[12,15] - uDX[12,15] * ldy[12,15]
   - udx[12,0] * lDY[12,0] + udx[12,0] * ldy[12,0] - udx[12,1] * lDY[12,1]
   + udx[12,1] * ldy[12,1] - udx[12,2] * lDY[12,2] + udx[12,2] * ldy[12,2]
   - udx[12,3] * lDY[12,3] + udx[12,3] * ldy[12,3] - udx[12,4] * lDY[12,4]
   + udx[12,4] * ldy[12,4] - udx[12,5] * lDY[12,5] + udx[12,5] * ldy[12,5]
   - udx[12,6] * lDY[12,6] + udx[12,6] * ldy[12,6] - udx[12,7] * lDY[12,7]
   + udx[12,7] * ldy[12,7] - udx[12,8] * lDY[12,8] + udx[12,8] * ldy[12,8]
   - udx[12,9] * lDY[12,9] + udx[12,9] * ldy[12,9]
   - udx[12,10] * lDY[12,10] + udx[12,10] * ldy[12,10]
   - udx[12,11] * lDY[12,11] + udx[12,11] * ldy[12,11]
   - udx[12,12] * lDY[12,12] + udx[12,12] * ldy[12,12]
   - udx[12,13] * lDY[12,13] + udx[12,13] * ldy[12,13]
   - udx[12,14] * lDY[12,14] + udx[12,14] * ldy[12,14]
   - udx[12,15] * lDY[12,15] + udx[12,15] * ldy[12,15] ] >= 1
Bounds
 ikSUM[0] <= 2
 ikSUM[1] <= 2
 ikSUM[2] <= 2
 ikSUM[3] <= 2
 ikSUM[4] <= 2
 ikSUM[5] <= 2
 ikSUM[6] <= 2
 ikSUM[7] <= 2
 ikSUM[8] <= 2
 ikSUM[9] <= 2
 ikSUM[10] <= 2
 ikSUM[11] <= 2
 ikSUM[12] <= 2
 ikSUM[13] <= 2
 ikSUM[14] <= 2
 ikSUM[15] <= 2
 gkSUM[0] <= 2
 gkSUM[1] <= 2
 gkSUM[2] <= 2
 gkSUM[3] <= 2
 gkSUM[4] <= 2
 gkSUM[5] <= 2
 gkSUM[6] <= 2
 gkSUM[7] <= 2
 gkSUM[8] <= 2
 gkSUM[9] <= 2
 gkSUM[10] <= 2
 gkSUM[11] <= 2
 gkSUM[12] <= 2
 gkSUM[13] <= 2
 gkSUM[14] <= 2
 gkSUM[15] <= 2
Binaries
 uDX[1,0] uDX[1,1] uDX[1,2] uDX[1,3] uDX[1,4] uDX[1,5] uDX[1,6] uDX[1,7]
 uDX[1,8] uDX[1,9] uDX[1,10] uDX[1,11] uDX[1,12] uDX[1,13] uDX[1,14]
 uDX[1,15] uDX[2,0] uDX[2,1] uDX[2,2] uDX[2,3] uDX[2,4] uDX[2,5] uDX[2,6]
 uDX[2,7] uDX[2,8] uDX[2,9] uDX[2,10] uDX[2,11] uDX[2,12] uDX[2,13]
 uDX[2,14] uDX[2,15] uDX[3,0] uDX[3,1] uDX[3,2] uDX[3,3] uDX[3,4] uDX[3,5]
 uDX[3,6] uDX[3,7] uDX[3,8] uDX[3,9] uDX[3,10] uDX[3,11] uDX[3,12]
 uDX[3,13] uDX[3,14] uDX[3,15] uDX[4,0] uDX[4,1] uDX[4,2] uDX[4,3] uDX[4,4]
 uDX[4,5] uDX[4,6] uDX[4,7] uDX[4,8] uDX[4,9] uDX[4,10] uDX[4,11] uDX[4,12]
 uDX[4,13] uDX[4,14] uDX[4,15] uDX[5,0] uDX[5,1] uDX[5,2] uDX[5,3] uDX[5,4]
 uDX[5,5] uDX[5,6] uDX[5,7] uDX[5,8] uDX[5,9] uDX[5,10] uDX[5,11] uDX[5,12]
 uDX[5,13] uDX[5,14] uDX[5,15] uDX[6,0] uDX[6,1] uDX[6,2] uDX[6,3] uDX[6,4]
 uDX[6,5] uDX[6,6] uDX[6,7] uDX[6,8] uDX[6,9] uDX[6,10] uDX[6,11] uDX[6,12]
 uDX[6,13] uDX[6,14] uDX[6,15] uDX[7,0] uDX[7,1] uDX[7,2] uDX[7,3] uDX[7,4]
 uDX[7,5] uDX[7,6] uDX[7,7] uDX[7,8] uDX[7,9] uDX[7,10] uDX[7,11] uDX[7,12]
 uDX[7,13] uDX[7,14] uDX[7,15] uDX[8,0] uDX[8,1] uDX[8,2] uDX[8,3] uDX[8,4]
 uDX[8,5] uDX[8,6] uDX[8,7] uDX[8,8] uDX[8,9] uDX[8,10] uDX[8,11] uDX[8,12]
 uDX[8,13] uDX[8,14] uDX[8,15] uDX[9,0] uDX[9,1] uDX[9,2] uDX[9,3] uDX[9,4]
 uDX[9,5] uDX[9,6] uDX[9,7] uDX[9,8] uDX[9,9] uDX[9,10] uDX[9,11] uDX[9,12]
 uDX[9,13] uDX[9,14] uDX[9,15] uDX[10,0] uDX[10,1] uDX[10,2] uDX[10,3]
 uDX[10,4] uDX[10,5] uDX[10,6] uDX[10,7] uDX[10,8] uDX[10,9] uDX[10,10]
 uDX[10,11] uDX[10,12] uDX[10,13] uDX[10,14] uDX[10,15] uDX[11,0] uDX[11,1]
 uDX[11,2] uDX[11,3] uDX[11,4] uDX[11,5] uDX[11,6] uDX[11,7] uDX[11,8]
 uDX[11,9] uDX[11,10] uDX[11,11] uDX[11,12] uDX[11,13] uDX[11,14]
 uDX[11,15] uDX[12,0] uDX[12,1] uDX[12,2] uDX[12,3] uDX[12,4] uDX[12,5]
 uDX[12,6] uDX[12,7] uDX[12,8] uDX[12,9] uDX[12,10] uDX[12,11] uDX[12,12]
 uDX[12,13] uDX[12,14] uDX[12,15] udx[1,0] udx[1,1] udx[1,2] udx[1,3]
 udx[1,4] udx[1,5] udx[1,6] udx[1,7] udx[1,8] udx[1,9] udx[1,10] udx[1,11]
 udx[1,12] udx[1,13] udx[1,14] udx[1,15] udx[2,0] udx[2,1] udx[2,2]
 udx[2,3] udx[2,4] udx[2,5] udx[2,6] udx[2,7] udx[2,8] udx[2,9] udx[2,10]
 udx[2,11] udx[2,12] udx[2,13] udx[2,14] udx[2,15] udx[3,0] udx[3,1]
 udx[3,2] udx[3,3] udx[3,4] udx[3,5] udx[3,6] udx[3,7] udx[3,8] udx[3,9]
 udx[3,10] udx[3,11] udx[3,12] udx[3,13] udx[3,14] udx[3,15] udx[4,0]
 udx[4,1] udx[4,2] udx[4,3] udx[4,4] udx[4,5] udx[4,6] udx[4,7] udx[4,8]
 udx[4,9] udx[4,10] udx[4,11] udx[4,12] udx[4,13] udx[4,14] udx[4,15]
 udx[5,0] udx[5,1] udx[5,2] udx[5,3] udx[5,4] udx[5,5] udx[5,6] udx[5,7]
 udx[5,8] udx[5,9] udx[5,10] udx[5,11] udx[5,12] udx[5,13] udx[5,14]
 udx[5,15] udx[6,0] udx[6,1] udx[6,2] udx[6,3] udx[6,4] udx[6,5] udx[6,6]
 udx[6,7] udx[6,8] udx[6,9] udx[6,10] udx[6,11] udx[6,12] udx[6,13]
 udx[6,14] udx[6,15] udx[7,0] udx[7,1] udx[7,2] udx[7,3] udx[7,4] udx[7,5]
 udx[7,6] udx[7,7] udx[7,8] udx[7,9] udx[7,10] udx[7,11] udx[7,12]
 udx[7,13] udx[7,14] udx[7,15] udx[8,0] udx[8,1] udx[8,2] udx[8,3] udx[8,4]
 udx[8,5] udx[8,6] udx[8,7] udx[8,8] udx[8,9] udx[8,10] udx[8,11] udx[8,12]
 udx[8,13] udx[8,14] udx[8,15] udx[9,0] udx[9,1] udx[9,2] udx[9,3] udx[9,4]
 udx[9,5] udx[9,6] udx[9,7] udx[9,8] udx[9,9] udx[9,10] udx[9,11] udx[9,12]
 udx[9,13] udx[9,14] udx[9,15] udx[10,0] udx[10,1] udx[10,2] udx[10,3]
 udx[10,4] udx[10,5] udx[10,6] udx[10,7] udx[10,8] udx[10,9] udx[10,10]
 udx[10,11] udx[10,12] udx[10,13] udx[10,14] udx[10,15] udx[11,0] udx[11,1]
 udx[11,2] udx[11,3] udx[11,4] udx[11,5] udx[11,6] udx[11,7] udx[11,8]
 udx[11,9] udx[11,10] udx[11,11] udx[11,12] udx[11,13] udx[11,14]
 udx[11,15] udx[12,0] udx[12,1] udx[12,2] udx[12,3] udx[12,4] udx[12,5]
 udx[12,6] udx[12,7] udx[12,8] udx[12,9] udx[12,10] udx[12,11] udx[12,12]
 udx[12,13] udx[12,14] udx[12,15] uDY[1,0] uDY[1,1] uDY[1,2] uDY[1,3]
 uDY[1,4] uDY[1,5] uDY[1,6] uDY[1,7] uDY[1,8] uDY[1,9] uDY[1,10] uDY[1,11]
 uDY[1,12] uDY[1,13] uDY[1,14] uDY[1,15] uDY[2,0] uDY[2,1] uDY[2,2]
 uDY[2,3] uDY[2,4] uDY[2,5] uDY[2,6] uDY[2,7] uDY[2,8] uDY[2,9] uDY[2,10]
 uDY[2,11] uDY[2,12] uDY[2,13] uDY[2,14] uDY[2,15] uDY[3,0] uDY[3,1]
 uDY[3,2] uDY[3,3] uDY[3,4] uDY[3,5] uDY[3,6] uDY[3,7] uDY[3,8] uDY[3,9]
 uDY[3,10] uDY[3,11] uDY[3,12] uDY[3,13] uDY[3,14] uDY[3,15] uDY[4,0]
 uDY[4,1] uDY[4,2] uDY[4,3] uDY[4,4] uDY[4,5] uDY[4,6] uDY[4,7] uDY[4,8]
 uDY[4,9] uDY[4,10] uDY[4,11] uDY[4,12] uDY[4,13] uDY[4,14] uDY[4,15]
 uDY[5,0] uDY[5,1] uDY[5,2] uDY[5,3] uDY[5,4] uDY[5,5] uDY[5,6] uDY[5,7]
 uDY[5,8] uDY[5,9] uDY[5,10] uDY[5,11] uDY[5,12] uDY[5,13] uDY[5,14]
 uDY[5,15] uDY[6,0] uDY[6,1] uDY[6,2] uDY[6,3] uDY[6,4] uDY[6,5] uDY[6,6]
 uDY[6,7] uDY[6,8] uDY[6,9] uDY[6,10] uDY[6,11] uDY[6,12] uDY[6,13]
 uDY[6,14] uDY[6,15] uDY[7,0] uDY[7,1] uDY[7,2] uDY[7,3] uDY[7,4] uDY[7,5]
 uDY[7,6] uDY[7,7] uDY[7,8] uDY[7,9] uDY[7,10] uDY[7,11] uDY[7,12]
 uDY[7,13] uDY[7,14] uDY[7,15] uDY[8,0] uDY[8,1] uDY[8,2] uDY[8,3] uDY[8,4]
 uDY[8,5] uDY[8,6] uDY[8,7] uDY[8,8] uDY[8,9] uDY[8,10] uDY[8,11] uDY[8,12]
 uDY[8,13] uDY[8,14] uDY[8,15] uDY[9,0] uDY[9,1] uDY[9,2] uDY[9,3] uDY[9,4]
 uDY[9,5] uDY[9,6] uDY[9,7] uDY[9,8] uDY[9,9] uDY[9,10] uDY[9,11] uDY[9,12]
 uDY[9,13] uDY[9,14] uDY[9,15] uDY[10,0] uDY[10,1] uDY[10,2] uDY[10,3]
 uDY[10,4] uDY[10,5] uDY[10,6] uDY[10,7] uDY[10,8] uDY[10,9] uDY[10,10]
 uDY[10,11] uDY[10,12] uDY[10,13] uDY[10,14] uDY[10,15] uDY[11,0] uDY[11,1]
 uDY[11,2] uDY[11,3] uDY[11,4] uDY[11,5] uDY[11,6] uDY[11,7] uDY[11,8]
 uDY[11,9] uDY[11,10] uDY[11,11] uDY[11,12] uDY[11,13] uDY[11,14]
 uDY[11,15] udy[1,0] udy[1,1] udy[1,2] udy[1,3] udy[1,4] udy[1,5] udy[1,6]
 udy[1,7] udy[1,8] udy[1,9] udy[1,10] udy[1,11] udy[1,12] udy[1,13]
 udy[1,14] udy[1,15] udy[2,0] udy[2,1] udy[2,2] udy[2,3] udy[2,4] udy[2,5]
 udy[2,6] udy[2,7] udy[2,8] udy[2,9] udy[2,10] udy[2,11] udy[2,12]
 udy[2,13] udy[2,14] udy[2,15] udy[3,0] udy[3,1] udy[3,2] udy[3,3] udy[3,4]
 udy[3,5] udy[3,6] udy[3,7] udy[3,8] udy[3,9] udy[3,10] udy[3,11] udy[3,12]
 udy[3,13] udy[3,14] udy[3,15] udy[4,0] udy[4,1] udy[4,2] udy[4,3] udy[4,4]
 udy[4,5] udy[4,6] udy[4,7] udy[4,8] udy[4,9] udy[4,10] udy[4,11] udy[4,12]
 udy[4,13] udy[4,14] udy[4,15] udy[5,0] udy[5,1] udy[5,2] udy[5,3] udy[5,4]
 udy[5,5] udy[5,6] udy[5,7] udy[5,8] udy[5,9] udy[5,10] udy[5,11] udy[5,12]
 udy[5,13] udy[5,14] udy[5,15] udy[6,0] udy[6,1] udy[6,2] udy[6,3] udy[6,4]
 udy[6,5] udy[6,6] udy[6,7] udy[6,8] udy[6,9] udy[6,10] udy[6,11] udy[6,12]
 udy[6,13] udy[6,14] udy[6,15] udy[7,0] udy[7,1] udy[7,2] udy[7,3] udy[7,4]
 udy[7,5] udy[7,6] udy[7,7] udy[7,8] udy[7,9] udy[7,10] udy[7,11] udy[7,12]
 udy[7,13] udy[7,14] udy[7,15] udy[8,0] udy[8,1] udy[8,2] udy[8,3] udy[8,4]
 udy[8,5] udy[8,6] udy[8,7] udy[8,8] udy[8,9] udy[8,10] udy[8,11] udy[8,12]
 udy[8,13] udy[8,14] udy[8,15] udy[9,0] udy[9,1] udy[9,2] udy[9,3] udy[9,4]
 udy[9,5] udy[9,6] udy[9,7] udy[9,8] udy[9,9] udy[9,10] udy[9,11] udy[9,12]
 udy[9,13] udy[9,14] udy[9,15] udy[10,0] udy[10,1] udy[10,2] udy[10,3]
 udy[10,4] udy[10,5] udy[10,6] udy[10,7] udy[10,8] udy[10,9] udy[10,10]
 udy[10,11] udy[10,12] udy[10,13] udy[10,14] udy[10,15] udy[11,0] udy[11,1]
 udy[11,2] udy[11,3] udy[11,4] udy[11,5] udy[11,6] udy[11,7] udy[11,8]
 udy[11,9] udy[11,10] udy[11,11] udy[11,12] udy[11,13] udy[11,14]
 udy[11,15] uDZ[1,0] uDZ[1,1] uDZ[1,2] uDZ[1,3] uDZ[1,4] uDZ[1,5] uDZ[1,6]
 uDZ[1,7] uDZ[1,8] uDZ[1,9] uDZ[1,10] uDZ[1,11] uDZ[1,12] uDZ[1,13]
 uDZ[1,14] uDZ[1,15] uDZ[2,0] uDZ[2,1] uDZ[2,2] uDZ[2,3] uDZ[2,4] uDZ[2,5]
 uDZ[2,6] uDZ[2,7] uDZ[2,8] uDZ[2,9] uDZ[2,10] uDZ[2,11] uDZ[2,12]
 uDZ[2,13] uDZ[2,14] uDZ[2,15] uDZ[3,0] uDZ[3,1] uDZ[3,2] uDZ[3,3] uDZ[3,4]
 uDZ[3,5] uDZ[3,6] uDZ[3,7] uDZ[3,8] uDZ[3,9] uDZ[3,10] uDZ[3,11] uDZ[3,12]
 uDZ[3,13] uDZ[3,14] uDZ[3,15] uDZ[4,0] uDZ[4,1] uDZ[4,2] uDZ[4,3] uDZ[4,4]
 uDZ[4,5] uDZ[4,6] uDZ[4,7] uDZ[4,8] uDZ[4,9] uDZ[4,10] uDZ[4,11] uDZ[4,12]
 uDZ[4,13] uDZ[4,14] uDZ[4,15] uDZ[5,0] uDZ[5,1] uDZ[5,2] uDZ[5,3] uDZ[5,4]
 uDZ[5,5] uDZ[5,6] uDZ[5,7] uDZ[5,8] uDZ[5,9] uDZ[5,10] uDZ[5,11] uDZ[5,12]
 uDZ[5,13] uDZ[5,14] uDZ[5,15] uDZ[6,0] uDZ[6,1] uDZ[6,2] uDZ[6,3] uDZ[6,4]
 uDZ[6,5] uDZ[6,6] uDZ[6,7] uDZ[6,8] uDZ[6,9] uDZ[6,10] uDZ[6,11] uDZ[6,12]
 uDZ[6,13] uDZ[6,14] uDZ[6,15] uDZ[7,0] uDZ[7,1] uDZ[7,2] uDZ[7,3] uDZ[7,4]
 uDZ[7,5] uDZ[7,6] uDZ[7,7] uDZ[7,8] uDZ[7,9] uDZ[7,10] uDZ[7,11] uDZ[7,12]
 uDZ[7,13] uDZ[7,14] uDZ[7,15] uDZ[8,0] uDZ[8,1] uDZ[8,2] uDZ[8,3] uDZ[8,4]
 uDZ[8,5] uDZ[8,6] uDZ[8,7] uDZ[8,8] uDZ[8,9] uDZ[8,10] uDZ[8,11] uDZ[8,12]
 uDZ[8,13] uDZ[8,14] uDZ[8,15] uDZ[9,0] uDZ[9,1] uDZ[9,2] uDZ[9,3] uDZ[9,4]
 uDZ[9,5] uDZ[9,6] uDZ[9,7] uDZ[9,8] uDZ[9,9] uDZ[9,10] uDZ[9,11] uDZ[9,12]
 uDZ[9,13] uDZ[9,14] uDZ[9,15] uDZ[10,0] uDZ[10,1] uDZ[10,2] uDZ[10,3]
 uDZ[10,4] uDZ[10,5] uDZ[10,6] uDZ[10,7] uDZ[10,8] uDZ[10,9] uDZ[10,10]
 uDZ[10,11] uDZ[10,12] uDZ[10,13] uDZ[10,14] uDZ[10,15] uDZ[11,0] uDZ[11,1]
 uDZ[11,2] uDZ[11,3] uDZ[11,4] uDZ[11,5] uDZ[11,6] uDZ[11,7] uDZ[11,8]
 uDZ[11,9] uDZ[11,10] uDZ[11,11] uDZ[11,12] uDZ[11,13] uDZ[11,14]
 uDZ[11,15] udz[1,0] udz[1,1] udz[1,2] udz[1,3] udz[1,4] udz[1,5] udz[1,6]
 udz[1,7] udz[1,8] udz[1,9] udz[1,10] udz[1,11] udz[1,12] udz[1,13]
 udz[1,14] udz[1,15] udz[2,0] udz[2,1] udz[2,2] udz[2,3] udz[2,4] udz[2,5]
 udz[2,6] udz[2,7] udz[2,8] udz[2,9] udz[2,10] udz[2,11] udz[2,12]
 udz[2,13] udz[2,14] udz[2,15] udz[3,0] udz[3,1] udz[3,2] udz[3,3] udz[3,4]
 udz[3,5] udz[3,6] udz[3,7] udz[3,8] udz[3,9] udz[3,10] udz[3,11] udz[3,12]
 udz[3,13] udz[3,14] udz[3,15] udz[4,0] udz[4,1] udz[4,2] udz[4,3] udz[4,4]
 udz[4,5] udz[4,6] udz[4,7] udz[4,8] udz[4,9] udz[4,10] udz[4,11] udz[4,12]
 udz[4,13] udz[4,14] udz[4,15] udz[5,0] udz[5,1] udz[5,2] udz[5,3] udz[5,4]
 udz[5,5] udz[5,6] udz[5,7] udz[5,8] udz[5,9] udz[5,10] udz[5,11] udz[5,12]
 udz[5,13] udz[5,14] udz[5,15] udz[6,0] udz[6,1] udz[6,2] udz[6,3] udz[6,4]
 udz[6,5] udz[6,6] udz[6,7] udz[6,8] udz[6,9] udz[6,10] udz[6,11] udz[6,12]
 udz[6,13] udz[6,14] udz[6,15] udz[7,0] udz[7,1] udz[7,2] udz[7,3] udz[7,4]
 udz[7,5] udz[7,6] udz[7,7] udz[7,8] udz[7,9] udz[7,10] udz[7,11] udz[7,12]
 udz[7,13] udz[7,14] udz[7,15] udz[8,0] udz[8,1] udz[8,2] udz[8,3] udz[8,4]
 udz[8,5] udz[8,6] udz[8,7] udz[8,8] udz[8,9] udz[8,10] udz[8,11] udz[8,12]
 udz[8,13] udz[8,14] udz[8,15] udz[9,0] udz[9,1] udz[9,2] udz[9,3] udz[9,4]
 udz[9,5] udz[9,6] udz[9,7] udz[9,8] udz[9,9] udz[9,10] udz[9,11] udz[9,12]
 udz[9,13] udz[9,14] udz[9,15] udz[10,0] udz[10,1] udz[10,2] udz[10,3]
 udz[10,4] udz[10,5] udz[10,6] udz[10,7] udz[10,8] udz[10,9] udz[10,10]
 udz[10,11] udz[10,12] udz[10,13] udz[10,14] udz[10,15] udz[11,0] udz[11,1]
 udz[11,2] udz[11,3] udz[11,4] udz[11,5] udz[11,6] udz[11,7] udz[11,8]
 udz[11,9] udz[11,10] udz[11,11] udz[11,12] udz[11,13] udz[11,14]
 udz[11,15] uDW[0,0] uDW[0,1] uDW[0,2] uDW[0,3] uDW[0,4] uDW[0,5] uDW[0,6]
 uDW[0,7] uDW[0,8] uDW[0,9] uDW[0,10] uDW[0,11] uDW[0,12] uDW[0,13]
 uDW[0,14] uDW[0,15] uDW[1,0] uDW[1,1] uDW[1,2] uDW[1,3] uDW[1,4] uDW[1,5]
 uDW[1,6] uDW[1,7] uDW[1,8] uDW[1,9] uDW[1,10] uDW[1,11] uDW[1,12]
 uDW[1,13] uDW[1,14] uDW[1,15] uDW[2,0] uDW[2,1] uDW[2,2] uDW[2,3] uDW[2,4]
 uDW[2,5] uDW[2,6] uDW[2,7] uDW[2,8] uDW[2,9] uDW[2,10] uDW[2,11] uDW[2,12]
 uDW[2,13] uDW[2,14] uDW[2,15] uDW[3,0] uDW[3,1] uDW[3,2] uDW[3,3] uDW[3,4]
 uDW[3,5] uDW[3,6] uDW[3,7] uDW[3,8] uDW[3,9] uDW[3,10] uDW[3,11] uDW[3,12]
 uDW[3,13] uDW[3,14] uDW[3,15] uDW[4,0] uDW[4,1] uDW[4,2] uDW[4,3] uDW[4,4]
 uDW[4,5] uDW[4,6] uDW[4,7] uDW[4,8] uDW[4,9] uDW[4,10] uDW[4,11] uDW[4,12]
 uDW[4,13] uDW[4,14] uDW[4,15] uDW[5,0] uDW[5,1] uDW[5,2] uDW[5,3] uDW[5,4]
 uDW[5,5] uDW[5,6] uDW[5,7] uDW[5,8] uDW[5,9] uDW[5,10] uDW[5,11] uDW[5,12]
 uDW[5,13] uDW[5,14] uDW[5,15] uDW[6,0] uDW[6,1] uDW[6,2] uDW[6,3] uDW[6,4]
 uDW[6,5] uDW[6,6] uDW[6,7] uDW[6,8] uDW[6,9] uDW[6,10] uDW[6,11] uDW[6,12]
 uDW[6,13] uDW[6,14] uDW[6,15] uDW[7,0] uDW[7,1] uDW[7,2] uDW[7,3] uDW[7,4]
 uDW[7,5] uDW[7,6] uDW[7,7] uDW[7,8] uDW[7,9] uDW[7,10] uDW[7,11] uDW[7,12]
 uDW[7,13] uDW[7,14] uDW[7,15] uDW[8,0] uDW[8,1] uDW[8,2] uDW[8,3] uDW[8,4]
 uDW[8,5] uDW[8,6] uDW[8,7] uDW[8,8] uDW[8,9] uDW[8,10] uDW[8,11] uDW[8,12]
 uDW[8,13] uDW[8,14] uDW[8,15] uDW[9,0] uDW[9,1] uDW[9,2] uDW[9,3] uDW[9,4]
 uDW[9,5] uDW[9,6] uDW[9,7] uDW[9,8] uDW[9,9] uDW[9,10] uDW[9,11] uDW[9,12]
 uDW[9,13] uDW[9,14] uDW[9,15] uDW[10,0] uDW[10,1] uDW[10,2] uDW[10,3]
 uDW[10,4] uDW[10,5] uDW[10,6] uDW[10,7] uDW[10,8] uDW[10,9] uDW[10,10]
 uDW[10,11] uDW[10,12] uDW[10,13] uDW[10,14] uDW[10,15] uDW[11,0] uDW[11,1]
 uDW[11,2] uDW[11,3] uDW[11,4] uDW[11,5] uDW[11,6] uDW[11,7] uDW[11,8]
 uDW[11,9] uDW[11,10] uDW[11,11] uDW[11,12] uDW[11,13] uDW[11,14]
 uDW[11,15] udw[0,0] udw[0,1] udw[0,2] udw[0,3] udw[0,4] udw[0,5] udw[0,6]
 udw[0,7] udw[0,8] udw[0,9] udw[0,10] udw[0,11] udw[0,12] udw[0,13]
 udw[0,14] udw[0,15] udw[1,0] udw[1,1] udw[1,2] udw[1,3] udw[1,4] udw[1,5]
 udw[1,6] udw[1,7] udw[1,8] udw[1,9] udw[1,10] udw[1,11] udw[1,12]
 udw[1,13] udw[1,14] udw[1,15] udw[2,0] udw[2,1] udw[2,2] udw[2,3] udw[2,4]
 udw[2,5] udw[2,6] udw[2,7] udw[2,8] udw[2,9] udw[2,10] udw[2,11] udw[2,12]
 udw[2,13] udw[2,14] udw[2,15] udw[3,0] udw[3,1] udw[3,2] udw[3,3] udw[3,4]
 udw[3,5] udw[3,6] udw[3,7] udw[3,8] udw[3,9] udw[3,10] udw[3,11] udw[3,12]
 udw[3,13] udw[3,14] udw[3,15] udw[4,0] udw[4,1] udw[4,2] udw[4,3] udw[4,4]
 udw[4,5] udw[4,6] udw[4,7] udw[4,8] udw[4,9] udw[4,10] udw[4,11] udw[4,12]
 udw[4,13] udw[4,14] udw[4,15] udw[5,0] udw[5,1] udw[5,2] udw[5,3] udw[5,4]
 udw[5,5] udw[5,6] udw[5,7] udw[5,8] udw[5,9] udw[5,10] udw[5,11] udw[5,12]
 udw[5,13] udw[5,14] udw[5,15] udw[6,0] udw[6,1] udw[6,2] udw[6,3] udw[6,4]
 udw[6,5] udw[6,6] udw[6,7] udw[6,8] udw[6,9] udw[6,10] udw[6,11] udw[6,12]
 udw[6,13] udw[6,14] udw[6,15] udw[7,0] udw[7,1] udw[7,2] udw[7,3] udw[7,4]
 udw[7,5] udw[7,6] udw[7,7] udw[7,8] udw[7,9] udw[7,10] udw[7,11] udw[7,12]
 udw[7,13] udw[7,14] udw[7,15] udw[8,0] udw[8,1] udw[8,2] udw[8,3] udw[8,4]
 udw[8,5] udw[8,6] udw[8,7] udw[8,8] udw[8,9] udw[8,10] udw[8,11] udw[8,12]
 udw[8,13] udw[8,14] udw[8,15] udw[9,0] udw[9,1] udw[9,2] udw[9,3] udw[9,4]
 udw[9,5] udw[9,6] udw[9,7] udw[9,8] udw[9,9] udw[9,10] udw[9,11] udw[9,12]
 udw[9,13] udw[9,14] udw[9,15] udw[10,0] udw[10,1] udw[10,2] udw[10,3]
 udw[10,4] udw[10,5] udw[10,6] udw[10,7] udw[10,8] udw[10,9] udw[10,10]
 udw[10,11] udw[10,12] udw[10,13] udw[10,14] udw[10,15] udw[11,0] udw[11,1]
 udw[11,2] udw[11,3] udw[11,4] udw[11,5] udw[11,6] udw[11,7] udw[11,8]
 udw[11,9] udw[11,10] udw[11,11] udw[11,12] udw[11,13] udw[11,14]
 udw[11,15] ustk[0,0] ustk[0,1] ustk[0,2] ustk[0,3] ustk[0,4] ustk[0,5]
 ustk[0,6] ustk[0,7] ustk[0,8] ustk[0,9] ustk[0,10] ustk[0,11] ustk[0,12]
 ustk[0,13] ustk[0,14] ustk[0,15] ustk[1,0] ustk[1,1] ustk[1,2] ustk[1,3]
 ustk[1,4] ustk[1,5] ustk[1,6] ustk[1,7] ustk[1,8] ustk[1,9] ustk[1,10]
 ustk[1,11] ustk[1,12] ustk[1,13] ustk[1,14] ustk[1,15] ustk[2,0] ustk[2,1]
 ustk[2,2] ustk[2,3] ustk[2,4] ustk[2,5] ustk[2,6] ustk[2,7] ustk[2,8]
 ustk[2,9] ustk[2,10] ustk[2,11] ustk[2,12] ustk[2,13] ustk[2,14]
 ustk[2,15] ustk[3,0] ustk[3,1] ustk[3,2] ustk[3,3] ustk[3,4] ustk[3,5]
 ustk[3,6] ustk[3,7] ustk[3,8] ustk[3,9] ustk[3,10] ustk[3,11] ustk[3,12]
 ustk[3,13] ustk[3,14] ustk[3,15] ustk[4,0] ustk[4,1] ustk[4,2] ustk[4,3]
 ustk[4,4] ustk[4,5] ustk[4,6] ustk[4,7] ustk[4,8] ustk[4,9] ustk[4,10]
 ustk[4,11] ustk[4,12] ustk[4,13] ustk[4,14] ustk[4,15] ustk[5,0] ustk[5,1]
 ustk[5,2] ustk[5,3] ustk[5,4] ustk[5,5] ustk[5,6] ustk[5,7] ustk[5,8]
 ustk[5,9] ustk[5,10] ustk[5,11] ustk[5,12] ustk[5,13] ustk[5,14]
 ustk[5,15] ustk[6,0] ustk[6,1] ustk[6,2] ustk[6,3] ustk[6,4] ustk[6,5]
 ustk[6,6] ustk[6,7] ustk[6,8] ustk[6,9] ustk[6,10] ustk[6,11] ustk[6,12]
 ustk[6,13] ustk[6,14] ustk[6,15] ustk[7,0] ustk[7,1] ustk[7,2] ustk[7,3]
 ustk[7,4] ustk[7,5] ustk[7,6] ustk[7,7] ustk[7,8] ustk[7,9] ustk[7,10]
 ustk[7,11] ustk[7,12] ustk[7,13] ustk[7,14] ustk[7,15] ustk[8,0] ustk[8,1]
 ustk[8,2] ustk[8,3] ustk[8,4] ustk[8,5] ustk[8,6] ustk[8,7] ustk[8,8]
 ustk[8,9] ustk[8,10] ustk[8,11] ustk[8,12] ustk[8,13] ustk[8,14]
 ustk[8,15] ustk[9,0] ustk[9,1] ustk[9,2] ustk[9,3] ustk[9,4] ustk[9,5]
 ustk[9,6] ustk[9,7] ustk[9,8] ustk[9,9] ustk[9,10] ustk[9,11] ustk[9,12]
 ustk[9,13] ustk[9,14] ustk[9,15] ustk[10,0] ustk[10,1] ustk[10,2]
 ustk[10,3] ustk[10,4] ustk[10,5] ustk[10,6] ustk[10,7] ustk[10,8]
 ustk[10,9] ustk[10,10] ustk[10,11] ustk[10,12] ustk[10,13] ustk[10,14]
 ustk[10,15] ustk[11,0] ustk[11,1] ustk[11,2] ustk[11,3] ustk[11,4]
 ustk[11,5] ustk[11,6] ustk[11,7] ustk[11,8] ustk[11,9] ustk[11,10]
 ustk[11,11] ustk[11,12] ustk[11,13] ustk[11,14] ustk[11,15] ustk[12,0]
 ustk[12,1] ustk[12,2] ustk[12,3] ustk[12,4] ustk[12,5] ustk[12,6]
 ustk[12,7] ustk[12,8] ustk[12,9] ustk[12,10] ustk[12,11] ustk[12,12]
 ustk[12,13] ustk[12,14] ustk[12,15] ustk[13,0] ustk[13,1] ustk[13,2]
 ustk[13,3] ustk[13,4] ustk[13,5] ustk[13,6] ustk[13,7] ustk[13,8]
 ustk[13,9] ustk[13,10] ustk[13,11] ustk[13,12] ustk[13,13] ustk[13,14]
 ustk[13,15] ustk[14,0] ustk[14,1] ustk[14,2] ustk[14,3] ustk[14,4]
 ustk[14,5] ustk[14,6] ustk[14,7] ustk[14,8] ustk[14,9] ustk[14,10]
 ustk[14,11] ustk[14,12] ustk[14,13] ustk[14,14] ustk[14,15] ustk[15,0]
 ustk[15,1] ustk[15,2] ustk[15,3] ustk[15,4] ustk[15,5] ustk[15,6]
 ustk[15,7] ustk[15,8] ustk[15,9] ustk[15,10] ustk[15,11] ustk[15,12]
 ustk[15,13] ustk[15,14] ustk[15,15] ustk[16,0] ustk[16,1] ustk[16,2]
 ustk[16,3] ustk[16,4] ustk[16,5] ustk[16,6] ustk[16,7] ustk[16,8]
 ustk[16,9] ustk[16,10] ustk[16,11] ustk[16,12] ustk[16,13] ustk[16,14]
 ustk[16,15] ustk[17,0] ustk[17,1] ustk[17,2] ustk[17,3] ustk[17,4]
 ustk[17,5] ustk[17,6] ustk[17,7] ustk[17,8] ustk[17,9] ustk[17,10]
 ustk[17,11] ustk[17,12] ustk[17,13] ustk[17,14] ustk[17,15] ustk[18,0]
 ustk[18,1] ustk[18,2] ustk[18,3] ustk[18,4] ustk[18,5] ustk[18,6]
 ustk[18,7] ustk[18,8] ustk[18,9] ustk[18,10] ustk[18,11] ustk[18,12]
 ustk[18,13] ustk[18,14] ustk[18,15] ustk[19,0] ustk[19,1] ustk[19,2]
 ustk[19,3] ustk[19,4] ustk[19,5] ustk[19,6] ustk[19,7] ustk[19,8]
 ustk[19,9] ustk[19,10] ustk[19,11] ustk[19,12] ustk[19,13] ustk[19,14]
 ustk[19,15] ustk[20,0] ustk[20,1] ustk[20,2] ustk[20,3] ustk[20,4]
 ustk[20,5] ustk[20,6] ustk[20,7] ustk[20,8] ustk[20,9] ustk[20,10]
 ustk[20,11] ustk[20,12] ustk[20,13] ustk[20,14] ustk[20,15] ustk[21,0]
 ustk[21,1] ustk[21,2] ustk[21,3] ustk[21,4] ustk[21,5] ustk[21,6]
 ustk[21,7] ustk[21,8] ustk[21,9] ustk[21,10] ustk[21,11] ustk[21,12]
 ustk[21,13] ustk[21,14] ustk[21,15] ustk[22,0] ustk[22,1] ustk[22,2]
 ustk[22,3] ustk[22,4] ustk[22,5] ustk[22,6] ustk[22,7] ustk[22,8]
 ustk[22,9] ustk[22,10] ustk[22,11] ustk[22,12] ustk[22,13] ustk[22,14]
 ustk[22,15] ustk[23,0] ustk[23,1] ustk[23,2] ustk[23,3] ustk[23,4]
 ustk[23,5] ustk[23,6] ustk[23,7] ustk[23,8] ustk[23,9] ustk[23,10]
 ustk[23,11] ustk[23,12] ustk[23,13] ustk[23,14] ustk[23,15] ustk[24,0]
 ustk[24,1] ustk[24,2] ustk[24,3] ustk[24,4] ustk[24,5] ustk[24,6]
 ustk[24,7] ustk[24,8] ustk[24,9] ustk[24,10] ustk[24,11] ustk[24,12]
 ustk[24,13] ustk[24,14] ustk[24,15] ustk[25,0] ustk[25,1] ustk[25,2]
 ustk[25,3] ustk[25,4] ustk[25,5] ustk[25,6] ustk[25,7] ustk[25,8]
 ustk[25,9] ustk[25,10] ustk[25,11] ustk[25,12] ustk[25,13] ustk[25,14]
 ustk[25,15] ustk[26,0] ustk[26,1] ustk[26,2] ustk[26,3] ustk[26,4]
 ustk[26,5] ustk[26,6] ustk[26,7] ustk[26,8] ustk[26,9] ustk[26,10]
 ustk[26,11] ustk[26,12] ustk[26,13] ustk[26,14] ustk[26,15] ustk[27,0]
 ustk[27,1] ustk[27,2] ustk[27,3] ustk[27,4] ustk[27,5] ustk[27,6]
 ustk[27,7] ustk[27,8] ustk[27,9] ustk[27,10] ustk[27,11] ustk[27,12]
 ustk[27,13] ustk[27,14] ustk[27,15] ustk[28,0] ustk[28,1] ustk[28,2]
 ustk[28,3] ustk[28,4] ustk[28,5] ustk[28,6] ustk[28,7] ustk[28,8]
 ustk[28,9] ustk[28,10] ustk[28,11] ustk[28,12] ustk[28,13] ustk[28,14]
 ustk[28,15] ustk[29,0] ustk[29,1] ustk[29,2] ustk[29,3] ustk[29,4]
 ustk[29,5] ustk[29,6] ustk[29,7] ustk[29,8] ustk[29,9] ustk[29,10]
 ustk[29,11] ustk[29,12] ustk[29,13] ustk[29,14] ustk[29,15] ustk[30,0]
 ustk[30,1] ustk[30,2] ustk[30,3] ustk[30,4] ustk[30,5] ustk[30,6]
 ustk[30,7] ustk[30,8] ustk[30,9] ustk[30,10] ustk[30,11] ustk[30,12]
 ustk[30,13] ustk[30,14] ustk[30,15] ustk[31,0] ustk[31,1] ustk[31,2]
 ustk[31,3] ustk[31,4] ustk[31,5] ustk[31,6] ustk[31,7] ustk[31,8]
 ustk[31,9] ustk[31,10] ustk[31,11] ustk[31,12] ustk[31,13] ustk[31,14]
 ustk[31,15] ustk[32,0] ustk[32,1] ustk[32,2] ustk[32,3] ustk[32,4]
 ustk[32,5] ustk[32,6] ustk[32,7] ustk[32,8] ustk[32,9] ustk[32,10]
 ustk[32,11] ustk[32,12] ustk[32,13] ustk[32,14] ustk[32,15] ustk[33,0]
 ustk[33,1] ustk[33,2] ustk[33,3] ustk[33,4] ustk[33,5] ustk[33,6]
 ustk[33,7] ustk[33,8] ustk[33,9] ustk[33,10] ustk[33,11] ustk[33,12]
 ustk[33,13] ustk[33,14] ustk[33,15] ustk[34,0] ustk[34,1] ustk[34,2]
 ustk[34,3] ustk[34,4] ustk[34,5] ustk[34,6] ustk[34,7] ustk[34,8]
 ustk[34,9] ustk[34,10] ustk[34,11] ustk[34,12] ustk[34,13] ustk[34,14]
 ustk[34,15] ustk[35,0] ustk[35,1] ustk[35,2] ustk[35,3] ustk[35,4]
 ustk[35,5] ustk[35,6] ustk[35,7] ustk[35,8] ustk[35,9] ustk[35,10]
 ustk[35,11] ustk[35,12] ustk[35,13] ustk[35,14] ustk[35,15] ustk[36,0]
 ustk[36,1] ustk[36,2] ustk[36,3] ustk[36,4] ustk[36,5] ustk[36,6]
 ustk[36,7] ustk[36,8] ustk[36,9] ustk[36,10] ustk[36,11] ustk[36,12]
 ustk[36,13] ustk[36,14] ustk[36,15] ustk[37,0] ustk[37,1] ustk[37,2]
 ustk[37,3] ustk[37,4] ustk[37,5] ustk[37,6] ustk[37,7] ustk[37,8]
 ustk[37,9] ustk[37,10] ustk[37,11] ustk[37,12] ustk[37,13] ustk[37,14]
 ustk[37,15] ustk[38,0] ustk[38,1] ustk[38,2] ustk[38,3] ustk[38,4]
 ustk[38,5] ustk[38,6] ustk[38,7] ustk[38,8] ustk[38,9] ustk[38,10]
 ustk[38,11] ustk[38,12] ustk[38,13] ustk[38,14] ustk[38,15] eqk[0] eqk[1]
 eqk[2] eqk[3] eqk[4] eqk[5] eqk[6] eqk[7] eqk[8] eqk[9] eqk[10] eqk[11]
 eqk[12] eqk[13] eqk[14] eqk[15] uLANE[0] uLANE[1] uLANE[2] uLANE[3]
 uLANE[4] uLANE[5] uLANE[6] uLANE[7] uLANE[8] uLANE[9] uLANE[10] uLANE[11]
 uLANE[12] uLANE[13] uLANE[14] uLANE[15] lDX[13,0] lDX[13,1] lDX[13,2]
 lDX[13,3] lDX[13,4] lDX[13,5] lDX[13,6] lDX[13,7] lDX[13,8] lDX[13,9]
 lDX[13,10] lDX[13,11] lDX[13,12] lDX[13,13] lDX[13,14] lDX[13,15]
 lDX[14,0] lDX[14,1] lDX[14,2] lDX[14,3] lDX[14,4] lDX[14,5] lDX[14,6]
 lDX[14,7] lDX[14,8] lDX[14,9] lDX[14,10] lDX[14,11] lDX[14,12] lDX[14,13]
 lDX[14,14] lDX[14,15] lDX[15,0] lDX[15,1] lDX[15,2] lDX[15,3] lDX[15,4]
 lDX[15,5] lDX[15,6] lDX[15,7] lDX[15,8] lDX[15,9] lDX[15,10] lDX[15,11]
 lDX[15,12] lDX[15,13] lDX[15,14] lDX[15,15] lDX[16,0] lDX[16,1] lDX[16,2]
 lDX[16,3] lDX[16,4] lDX[16,5] lDX[16,6] lDX[16,7] lDX[16,8] lDX[16,9]
 lDX[16,10] lDX[16,11] lDX[16,12] lDX[16,13] lDX[16,14] lDX[16,15]
 lDX[17,0] lDX[17,1] lDX[17,2] lDX[17,3] lDX[17,4] lDX[17,5] lDX[17,6]
 lDX[17,7] lDX[17,8] lDX[17,9] lDX[17,10] lDX[17,11] lDX[17,12] lDX[17,13]
 lDX[17,14] lDX[17,15] lDX[18,0] lDX[18,1] lDX[18,2] lDX[18,3] lDX[18,4]
 lDX[18,5] lDX[18,6] lDX[18,7] lDX[18,8] lDX[18,9] lDX[18,10] lDX[18,11]
 lDX[18,12] lDX[18,13] lDX[18,14] lDX[18,15] lDX[19,0] lDX[19,1] lDX[19,2]
 lDX[19,3] lDX[19,4] lDX[19,5] lDX[19,6] lDX[19,7] lDX[19,8] lDX[19,9]
 lDX[19,10] lDX[19,11] lDX[19,12] lDX[19,13] lDX[19,14] lDX[19,15]
 lDX[20,0] lDX[20,1] lDX[20,2] lDX[20,3] lDX[20,4] lDX[20,5] lDX[20,6]
 lDX[20,7] lDX[20,8] lDX[20,9] lDX[20,10] lDX[20,11] lDX[20,12] lDX[20,13]
 lDX[20,14] lDX[20,15] lDX[21,0] lDX[21,1] lDX[21,2] lDX[21,3] lDX[21,4]
 lDX[21,5] lDX[21,6] lDX[21,7] lDX[21,8] lDX[21,9] lDX[21,10] lDX[21,11]
 lDX[21,12] lDX[21,13] lDX[21,14] lDX[21,15] lDX[22,0] lDX[22,1] lDX[22,2]
 lDX[22,3] lDX[22,4] lDX[22,5] lDX[22,6] lDX[22,7] lDX[22,8] lDX[22,9]
 lDX[22,10] lDX[22,11] lDX[22,12] lDX[22,13] lDX[22,14] lDX[22,15]
 lDX[23,0] lDX[23,1] lDX[23,2] lDX[23,3] lDX[23,4] lDX[23,5] lDX[23,6]
 lDX[23,7] lDX[23,8] lDX[23,9] lDX[23,10] lDX[23,11] lDX[23,12] lDX[23,13]
 lDX[23,14] lDX[23,15] lDX[24,0] lDX[24,1] lDX[24,2] lDX[24,3] lDX[24,4]
 lDX[24,5] lDX[24,6] lDX[24,7] lDX[24,8] lDX[24,9] lDX[24,10] lDX[24,11]
 lDX[24,12] lDX[24,13] lDX[24,14] lDX[24,15] ldx[13,0] ldx[13,1] ldx[13,2]
 ldx[13,3] ldx[13,4] ldx[13,5] ldx[13,6] ldx[13,7] ldx[13,8] ldx[13,9]
 ldx[13,10] ldx[13,11] ldx[13,12] ldx[13,13] ldx[13,14] ldx[13,15]
 ldx[14,0] ldx[14,1] ldx[14,2] ldx[14,3] ldx[14,4] ldx[14,5] ldx[14,6]
 ldx[14,7] ldx[14,8] ldx[14,9] ldx[14,10] ldx[14,11] ldx[14,12] ldx[14,13]
 ldx[14,14] ldx[14,15] ldx[15,0] ldx[15,1] ldx[15,2] ldx[15,3] ldx[15,4]
 ldx[15,5] ldx[15,6] ldx[15,7] ldx[15,8] ldx[15,9] ldx[15,10] ldx[15,11]
 ldx[15,12] ldx[15,13] ldx[15,14] ldx[15,15] ldx[16,0] ldx[16,1] ldx[16,2]
 ldx[16,3] ldx[16,4] ldx[16,5] ldx[16,6] ldx[16,7] ldx[16,8] ldx[16,9]
 ldx[16,10] ldx[16,11] ldx[16,12] ldx[16,13] ldx[16,14] ldx[16,15]
 ldx[17,0] ldx[17,1] ldx[17,2] ldx[17,3] ldx[17,4] ldx[17,5] ldx[17,6]
 ldx[17,7] ldx[17,8] ldx[17,9] ldx[17,10] ldx[17,11] ldx[17,12] ldx[17,13]
 ldx[17,14] ldx[17,15] ldx[18,0] ldx[18,1] ldx[18,2] ldx[18,3] ldx[18,4]
 ldx[18,5] ldx[18,6] ldx[18,7] ldx[18,8] ldx[18,9] ldx[18,10] ldx[18,11]
 ldx[18,12] ldx[18,13] ldx[18,14] ldx[18,15] ldx[19,0] ldx[19,1] ldx[19,2]
 ldx[19,3] ldx[19,4] ldx[19,5] ldx[19,6] ldx[19,7] ldx[19,8] ldx[19,9]
 ldx[19,10] ldx[19,11] ldx[19,12] ldx[19,13] ldx[19,14] ldx[19,15]
 ldx[20,0] ldx[20,1] ldx[20,2] ldx[20,3] ldx[20,4] ldx[20,5] ldx[20,6]
 ldx[20,7] ldx[20,8] ldx[20,9] ldx[20,10] ldx[20,11] ldx[20,12] ldx[20,13]
 ldx[20,14] ldx[20,15] ldx[21,0] ldx[21,1] ldx[21,2] ldx[21,3] ldx[21,4]
 ldx[21,5] ldx[21,6] ldx[21,7] ldx[21,8] ldx[21,9] ldx[21,10] ldx[21,11]
 ldx[21,12] ldx[21,13] ldx[21,14] ldx[21,15] ldx[22,0] ldx[22,1] ldx[22,2]
 ldx[22,3] ldx[22,4] ldx[22,5] ldx[22,6] ldx[22,7] ldx[22,8] ldx[22,9]
 ldx[22,10] ldx[22,11] ldx[22,12] ldx[22,13] ldx[22,14] ldx[22,15]
 ldx[23,0] ldx[23,1] ldx[23,2] ldx[23,3] ldx[23,4] ldx[23,5] ldx[23,6]
 ldx[23,7] ldx[23,8] ldx[23,9] ldx[23,10] ldx[23,11] ldx[23,12] ldx[23,13]
 ldx[23,14] ldx[23,15] ldx[24,0] ldx[24,1] ldx[24,2] ldx[24,3] ldx[24,4]
 ldx[24,5] ldx[24,6] ldx[24,7] ldx[24,8] ldx[24,9] ldx[24,10] ldx[24,11]
 ldx[24,12] ldx[24,13] ldx[24,14] ldx[24,15] lDY[12,0] lDY[12,1] lDY[12,2]
 lDY[12,3] lDY[12,4] lDY[12,5] lDY[12,6] lDY[12,7] lDY[12,8] lDY[12,9]
 lDY[12,10] lDY[12,11] lDY[12,12] lDY[12,13] lDY[12,14] lDY[12,15]
 lDY[13,0] lDY[13,1] lDY[13,2] lDY[13,3] lDY[13,4] lDY[13,5] lDY[13,6]
 lDY[13,7] lDY[13,8] lDY[13,9] lDY[13,10] lDY[13,11] lDY[13,12] lDY[13,13]
 lDY[13,14] lDY[13,15] lDY[14,0] lDY[14,1] lDY[14,2] lDY[14,3] lDY[14,4]
 lDY[14,5] lDY[14,6] lDY[14,7] lDY[14,8] lDY[14,9] lDY[14,10] lDY[14,11]
 lDY[14,12] lDY[14,13] lDY[14,14] lDY[14,15] lDY[15,0] lDY[15,1] lDY[15,2]
 lDY[15,3] lDY[15,4] lDY[15,5] lDY[15,6] lDY[15,7] lDY[15,8] lDY[15,9]
 lDY[15,10] lDY[15,11] lDY[15,12] lDY[15,13] lDY[15,14] lDY[15,15]
 lDY[16,0] lDY[16,1] lDY[16,2] lDY[16,3] lDY[16,4] lDY[16,5] lDY[16,6]
 lDY[16,7] lDY[16,8] lDY[16,9] lDY[16,10] lDY[16,11] lDY[16,12] lDY[16,13]
 lDY[16,14] lDY[16,15] lDY[17,0] lDY[17,1] lDY[17,2] lDY[17,3] lDY[17,4]
 lDY[17,5] lDY[17,6] lDY[17,7] lDY[17,8] lDY[17,9] lDY[17,10] lDY[17,11]
 lDY[17,12] lDY[17,13] lDY[17,14] lDY[17,15] lDY[18,0] lDY[18,1] lDY[18,2]
 lDY[18,3] lDY[18,4] lDY[18,5] lDY[18,6] lDY[18,7] lDY[18,8] lDY[18,9]
 lDY[18,10] lDY[18,11] lDY[18,12] lDY[18,13] lDY[18,14] lDY[18,15]
 lDY[19,0] lDY[19,1] lDY[19,2] lDY[19,3] lDY[19,4] lDY[19,5] lDY[19,6]
 lDY[19,7] lDY[19,8] lDY[19,9] lDY[19,10] lDY[19,11] lDY[19,12] lDY[19,13]
 lDY[19,14] lDY[19,15] lDY[20,0] lDY[20,1] lDY[20,2] lDY[20,3] lDY[20,4]
 lDY[20,5] lDY[20,6] lDY[20,7] lDY[20,8] lDY[20,9] lDY[20,10] lDY[20,11]
 lDY[20,12] lDY[20,13] lDY[20,14] lDY[20,15] lDY[21,0] lDY[21,1] lDY[21,2]
 lDY[21,3] lDY[21,4] lDY[21,5] lDY[21,6] lDY[21,7] lDY[21,8] lDY[21,9]
 lDY[21,10] lDY[21,11] lDY[21,12] lDY[21,13] lDY[21,14] lDY[21,15]
 lDY[22,0] lDY[22,1] lDY[22,2] lDY[22,3] lDY[22,4] lDY[22,5] lDY[22,6]
 lDY[22,7] lDY[22,8] lDY[22,9] lDY[22,10] lDY[22,11] lDY[22,12] lDY[22,13]
 lDY[22,14] lDY[22,15] lDY[23,0] lDY[23,1] lDY[23,2] lDY[23,3] lDY[23,4]
 lDY[23,5] lDY[23,6] lDY[23,7] lDY[23,8] lDY[23,9] lDY[23,10] lDY[23,11]
 lDY[23,12] lDY[23,13] lDY[23,14] lDY[23,15] ldy[12,0] ldy[12,1] ldy[12,2]
 ldy[12,3] ldy[12,4] ldy[12,5] ldy[12,6] ldy[12,7] ldy[12,8] ldy[12,9]
 ldy[12,10] ldy[12,11] ldy[12,12] ldy[12,13] ldy[12,14] ldy[12,15]
 ldy[13,0] ldy[13,1] ldy[13,2] ldy[13,3] ldy[13,4] ldy[13,5] ldy[13,6]
 ldy[13,7] ldy[13,8] ldy[13,9] ldy[13,10] ldy[13,11] ldy[13,12] ldy[13,13]
 ldy[13,14] ldy[13,15] ldy[14,0] ldy[14,1] ldy[14,2] ldy[14,3] ldy[14,4]
 ldy[14,5] ldy[14,6] ldy[14,7] ldy[14,8] ldy[14,9] ldy[14,10] ldy[14,11]
 ldy[14,12] ldy[14,13] ldy[14,14] ldy[14,15] ldy[15,0] ldy[15,1] ldy[15,2]
 ldy[15,3] ldy[15,4] ldy[15,5] ldy[15,6] ldy[15,7] ldy[15,8] ldy[15,9]
 ldy[15,10] ldy[15,11] ldy[15,12] ldy[15,13] ldy[15,14] ldy[15,15]
 ldy[16,0] ldy[16,1] ldy[16,2] ldy[16,3] ldy[16,4] ldy[16,5] ldy[16,6]
 ldy[16,7] ldy[16,8] ldy[16,9] ldy[16,10] ldy[16,11] ldy[16,12] ldy[16,13]
 ldy[16,14] ldy[16,15] ldy[17,0] ldy[17,1] ldy[17,2] ldy[17,3] ldy[17,4]
 ldy[17,5] ldy[17,6] ldy[17,7] ldy[17,8] ldy[17,9] ldy[17,10] ldy[17,11]
 ldy[17,12] ldy[17,13] ldy[17,14] ldy[17,15] ldy[18,0] ldy[18,1] ldy[18,2]
 ldy[18,3] ldy[18,4] ldy[18,5] ldy[18,6] ldy[18,7] ldy[18,8] ldy[18,9]
 ldy[18,10] ldy[18,11] ldy[18,12] ldy[18,13] ldy[18,14] ldy[18,15]
 ldy[19,0] ldy[19,1] ldy[19,2] ldy[19,3] ldy[19,4] ldy[19,5] ldy[19,6]
 ldy[19,7] ldy[19,8] ldy[19,9] ldy[19,10] ldy[19,11] ldy[19,12] ldy[19,13]
 ldy[19,14] ldy[19,15] ldy[20,0] ldy[20,1] ldy[20,2] ldy[20,3] ldy[20,4]
 ldy[20,5] ldy[20,6] ldy[20,7] ldy[20,8] ldy[20,9] ldy[20,10] ldy[20,11]
 ldy[20,12] ldy[20,13] ldy[20,14] ldy[20,15] ldy[21,0] ldy[21,1] ldy[21,2]
 ldy[21,3] ldy[21,4] ldy[21,5] ldy[21,6] ldy[21,7] ldy[21,8] ldy[21,9]
 ldy[21,10] ldy[21,11] ldy[21,12] ldy[21,13] ldy[21,14] ldy[21,15]
 ldy[22,0] ldy[22,1] ldy[22,2] ldy[22,3] ldy[22,4] ldy[22,5] ldy[22,6]
 ldy[22,7] ldy[22,8] ldy[22,9] ldy[22,10] ldy[22,11] ldy[22,12] ldy[22,13]
 ldy[22,14] ldy[22,15] ldy[23,0] ldy[23,1] ldy[23,2] ldy[23,3] ldy[23,4]
 ldy[23,5] ldy[23,6] ldy[23,7] ldy[23,8] ldy[23,9] ldy[23,10] ldy[23,11]
 ldy[23,12] ldy[23,13] ldy[23,14] ldy[23,15] lDZ[12,0] lDZ[12,1] lDZ[12,2]
 lDZ[12,3] lDZ[12,4] lDZ[12,5] lDZ[12,6] lDZ[12,7] lDZ[12,8] lDZ[12,9]
 lDZ[12,10] lDZ[12,11] lDZ[12,12] lDZ[12,13] lDZ[12,14] lDZ[12,15]
 lDZ[13,0] lDZ[13,1] lDZ[13,2] lDZ[13,3] lDZ[13,4] lDZ[13,5] lDZ[13,6]
 lDZ[13,7] lDZ[13,8] lDZ[13,9] lDZ[13,10] lDZ[13,11] lDZ[13,12] lDZ[13,13]
 lDZ[13,14] lDZ[13,15] lDZ[14,0] lDZ[14,1] lDZ[14,2] lDZ[14,3] lDZ[14,4]
 lDZ[14,5] lDZ[14,6] lDZ[14,7] lDZ[14,8] lDZ[14,9] lDZ[14,10] lDZ[14,11]
 lDZ[14,12] lDZ[14,13] lDZ[14,14] lDZ[14,15] lDZ[15,0] lDZ[15,1] lDZ[15,2]
 lDZ[15,3] lDZ[15,4] lDZ[15,5] lDZ[15,6] lDZ[15,7] lDZ[15,8] lDZ[15,9]
 lDZ[15,10] lDZ[15,11] lDZ[15,12] lDZ[15,13] lDZ[15,14] lDZ[15,15]
 lDZ[16,0] lDZ[16,1] lDZ[16,2] lDZ[16,3] lDZ[16,4] lDZ[16,5] lDZ[16,6]
 lDZ[16,7] lDZ[16,8] lDZ[16,9] lDZ[16,10] lDZ[16,11] lDZ[16,12] lDZ[16,13]
 lDZ[16,14] lDZ[16,15] lDZ[17,0] lDZ[17,1] lDZ[17,2] lDZ[17,3] lDZ[17,4]
 lDZ[17,5] lDZ[17,6] lDZ[17,7] lDZ[17,8] lDZ[17,9] lDZ[17,10] lDZ[17,11]
 lDZ[17,12] lDZ[17,13] lDZ[17,14] lDZ[17,15] lDZ[18,0] lDZ[18,1] lDZ[18,2]
 lDZ[18,3] lDZ[18,4] lDZ[18,5] lDZ[18,6] lDZ[18,7] lDZ[18,8] lDZ[18,9]
 lDZ[18,10] lDZ[18,11] lDZ[18,12] lDZ[18,13] lDZ[18,14] lDZ[18,15]
 lDZ[19,0] lDZ[19,1] lDZ[19,2] lDZ[19,3] lDZ[19,4] lDZ[19,5] lDZ[19,6]
 lDZ[19,7] lDZ[19,8] lDZ[19,9] lDZ[19,10] lDZ[19,11] lDZ[19,12] lDZ[19,13]
 lDZ[19,14] lDZ[19,15] lDZ[20,0] lDZ[20,1] lDZ[20,2] lDZ[20,3] lDZ[20,4]
 lDZ[20,5] lDZ[20,6] lDZ[20,7] lDZ[20,8] lDZ[20,9] lDZ[20,10] lDZ[20,11]
 lDZ[20,12] lDZ[20,13] lDZ[20,14] lDZ[20,15] lDZ[21,0] lDZ[21,1] lDZ[21,2]
 lDZ[21,3] lDZ[21,4] lDZ[21,5] lDZ[21,6] lDZ[21,7] lDZ[21,8] lDZ[21,9]
 lDZ[21,10] lDZ[21,11] lDZ[21,12] lDZ[21,13] lDZ[21,14] lDZ[21,15]
 lDZ[22,0] lDZ[22,1] lDZ[22,2] lDZ[22,3] lDZ[22,4] lDZ[22,5] lDZ[22,6]
 lDZ[22,7] lDZ[22,8] lDZ[22,9] lDZ[22,10] lDZ[22,11] lDZ[22,12] lDZ[22,13]
 lDZ[22,14] lDZ[22,15] lDZ[23,0] lDZ[23,1] lDZ[23,2] lDZ[23,3] lDZ[23,4]
 lDZ[23,5] lDZ[23,6] lDZ[23,7] lDZ[23,8] lDZ[23,9] lDZ[23,10] lDZ[23,11]
 lDZ[23,12] lDZ[23,13] lDZ[23,14] lDZ[23,15] ldz[12,0] ldz[12,1] ldz[12,2]
 ldz[12,3] ldz[12,4] ldz[12,5] ldz[12,6] ldz[12,7] ldz[12,8] ldz[12,9]
 ldz[12,10] ldz[12,11] ldz[12,12] ldz[12,13] ldz[12,14] ldz[12,15]
 ldz[13,0] ldz[13,1] ldz[13,2] ldz[13,3] ldz[13,4] ldz[13,5] ldz[13,6]
 ldz[13,7] ldz[13,8] ldz[13,9] ldz[13,10] ldz[13,11] ldz[13,12] ldz[13,13]
 ldz[13,14] ldz[13,15] ldz[14,0] ldz[14,1] ldz[14,2] ldz[14,3] ldz[14,4]
 ldz[14,5] ldz[14,6] ldz[14,7] ldz[14,8] ldz[14,9] ldz[14,10] ldz[14,11]
 ldz[14,12] ldz[14,13] ldz[14,14] ldz[14,15] ldz[15,0] ldz[15,1] ldz[15,2]
 ldz[15,3] ldz[15,4] ldz[15,5] ldz[15,6] ldz[15,7] ldz[15,8] ldz[15,9]
 ldz[15,10] ldz[15,11] ldz[15,12] ldz[15,13] ldz[15,14] ldz[15,15]
 ldz[16,0] ldz[16,1] ldz[16,2] ldz[16,3] ldz[16,4] ldz[16,5] ldz[16,6]
 ldz[16,7] ldz[16,8] ldz[16,9] ldz[16,10] ldz[16,11] ldz[16,12] ldz[16,13]
 ldz[16,14] ldz[16,15] ldz[17,0] ldz[17,1] ldz[17,2] ldz[17,3] ldz[17,4]
 ldz[17,5] ldz[17,6] ldz[17,7] ldz[17,8] ldz[17,9] ldz[17,10] ldz[17,11]
 ldz[17,12] ldz[17,13] ldz[17,14] ldz[17,15] ldz[18,0] ldz[18,1] ldz[18,2]
 ldz[18,3] ldz[18,4] ldz[18,5] ldz[18,6] ldz[18,7] ldz[18,8] ldz[18,9]
 ldz[18,10] ldz[18,11] ldz[18,12] ldz[18,13] ldz[18,14] ldz[18,15]
 ldz[19,0] ldz[19,1] ldz[19,2] ldz[19,3] ldz[19,4] ldz[19,5] ldz[19,6]
 ldz[19,7] ldz[19,8] ldz[19,9] ldz[19,10] ldz[19,11] ldz[19,12] ldz[19,13]
 ldz[19,14] ldz[19,15] ldz[20,0] ldz[20,1] ldz[20,2] ldz[20,3] ldz[20,4]
 ldz[20,5] ldz[20,6] ldz[20,7] ldz[20,8] ldz[20,9] ldz[20,10] ldz[20,11]
 ldz[20,12] ldz[20,13] ldz[20,14] ldz[20,15] ldz[21,0] ldz[21,1] ldz[21,2]
 ldz[21,3] ldz[21,4] ldz[21,5] ldz[21,6] ldz[21,7] ldz[21,8] ldz[21,9]
 ldz[21,10] ldz[21,11] ldz[21,12] ldz[21,13] ldz[21,14] ldz[21,15]
 ldz[22,0] ldz[22,1] ldz[22,2] ldz[22,3] ldz[22,4] ldz[22,5] ldz[22,6]
 ldz[22,7] ldz[22,8] ldz[22,9] ldz[22,10] ldz[22,11] ldz[22,12] ldz[22,13]
 ldz[22,14] ldz[22,15] ldz[23,0] ldz[23,1] ldz[23,2] ldz[23,3] ldz[23,4]
 ldz[23,5] ldz[23,6] ldz[23,7] ldz[23,8] ldz[23,9] ldz[23,10] ldz[23,11]
 ldz[23,12] ldz[23,13] ldz[23,14] ldz[23,15] lDW[12,0] lDW[12,1] lDW[12,2]
 lDW[12,3] lDW[12,4] lDW[12,5] lDW[12,6] lDW[12,7] lDW[12,8] lDW[12,9]
 lDW[12,10] lDW[12,11] lDW[12,12] lDW[12,13] lDW[12,14] lDW[12,15]
 lDW[13,0] lDW[13,1] lDW[13,2] lDW[13,3] lDW[13,4] lDW[13,5] lDW[13,6]
 lDW[13,7] lDW[13,8] lDW[13,9] lDW[13,10] lDW[13,11] lDW[13,12] lDW[13,13]
 lDW[13,14] lDW[13,15] lDW[14,0] lDW[14,1] lDW[14,2] lDW[14,3] lDW[14,4]
 lDW[14,5] lDW[14,6] lDW[14,7] lDW[14,8] lDW[14,9] lDW[14,10] lDW[14,11]
 lDW[14,12] lDW[14,13] lDW[14,14] lDW[14,15] lDW[15,0] lDW[15,1] lDW[15,2]
 lDW[15,3] lDW[15,4] lDW[15,5] lDW[15,6] lDW[15,7] lDW[15,8] lDW[15,9]
 lDW[15,10] lDW[15,11] lDW[15,12] lDW[15,13] lDW[15,14] lDW[15,15]
 lDW[16,0] lDW[16,1] lDW[16,2] lDW[16,3] lDW[16,4] lDW[16,5] lDW[16,6]
 lDW[16,7] lDW[16,8] lDW[16,9] lDW[16,10] lDW[16,11] lDW[16,12] lDW[16,13]
 lDW[16,14] lDW[16,15] lDW[17,0] lDW[17,1] lDW[17,2] lDW[17,3] lDW[17,4]
 lDW[17,5] lDW[17,6] lDW[17,7] lDW[17,8] lDW[17,9] lDW[17,10] lDW[17,11]
 lDW[17,12] lDW[17,13] lDW[17,14] lDW[17,15] lDW[18,0] lDW[18,1] lDW[18,2]
 lDW[18,3] lDW[18,4] lDW[18,5] lDW[18,6] lDW[18,7] lDW[18,8] lDW[18,9]
 lDW[18,10] lDW[18,11] lDW[18,12] lDW[18,13] lDW[18,14] lDW[18,15]
 lDW[19,0] lDW[19,1] lDW[19,2] lDW[19,3] lDW[19,4] lDW[19,5] lDW[19,6]
 lDW[19,7] lDW[19,8] lDW[19,9] lDW[19,10] lDW[19,11] lDW[19,12] lDW[19,13]
 lDW[19,14] lDW[19,15] lDW[20,0] lDW[20,1] lDW[20,2] lDW[20,3] lDW[20,4]
 lDW[20,5] lDW[20,6] lDW[20,7] lDW[20,8] lDW[20,9] lDW[20,10] lDW[20,11]
 lDW[20,12] lDW[20,13] lDW[20,14] lDW[20,15] lDW[21,0] lDW[21,1] lDW[21,2]
 lDW[21,3] lDW[21,4] lDW[21,5] lDW[21,6] lDW[21,7] lDW[21,8] lDW[21,9]
 lDW[21,10] lDW[21,11] lDW[21,12] lDW[21,13] lDW[21,14] lDW[21,15]
 lDW[22,0] lDW[22,1] lDW[22,2] lDW[22,3] lDW[22,4] lDW[22,5] lDW[22,6]
 lDW[22,7] lDW[22,8] lDW[22,9] lDW[22,10] lDW[22,11] lDW[22,12] lDW[22,13]
 lDW[22,14] lDW[22,15] lDW[23,0] lDW[23,1] lDW[23,2] lDW[23,3] lDW[23,4]
 lDW[23,5] lDW[23,6] lDW[23,7] lDW[23,8] lDW[23,9] lDW[23,10] lDW[23,11]
 lDW[23,12] lDW[23,13] lDW[23,14] lDW[23,15] ldw[12,0] ldw[12,1] ldw[12,2]
 ldw[12,3] ldw[12,4] ldw[12,5] ldw[12,6] ldw[12,7] ldw[12,8] ldw[12,9]
 ldw[12,10] ldw[12,11] ldw[12,12] ldw[12,13] ldw[12,14] ldw[12,15]
 ldw[13,0] ldw[13,1] ldw[13,2] ldw[13,3] ldw[13,4] ldw[13,5] ldw[13,6]
 ldw[13,7] ldw[13,8] ldw[13,9] ldw[13,10] ldw[13,11] ldw[13,12] ldw[13,13]
 ldw[13,14] ldw[13,15] ldw[14,0] ldw[14,1] ldw[14,2] ldw[14,3] ldw[14,4]
 ldw[14,5] ldw[14,6] ldw[14,7] ldw[14,8] ldw[14,9] ldw[14,10] ldw[14,11]
 ldw[14,12] ldw[14,13] ldw[14,14] ldw[14,15] ldw[15,0] ldw[15,1] ldw[15,2]
 ldw[15,3] ldw[15,4] ldw[15,5] ldw[15,6] ldw[15,7] ldw[15,8] ldw[15,9]
 ldw[15,10] ldw[15,11] ldw[15,12] ldw[15,13] ldw[15,14] ldw[15,15]
 ldw[16,0] ldw[16,1] ldw[16,2] ldw[16,3] ldw[16,4] ldw[16,5] ldw[16,6]
 ldw[16,7] ldw[16,8] ldw[16,9] ldw[16,10] ldw[16,11] ldw[16,12] ldw[16,13]
 ldw[16,14] ldw[16,15] ldw[17,0] ldw[17,1] ldw[17,2] ldw[17,3] ldw[17,4]
 ldw[17,5] ldw[17,6] ldw[17,7] ldw[17,8] ldw[17,9] ldw[17,10] ldw[17,11]
 ldw[17,12] ldw[17,13] ldw[17,14] ldw[17,15] ldw[18,0] ldw[18,1] ldw[18,2]
 ldw[18,3] ldw[18,4] ldw[18,5] ldw[18,6] ldw[18,7] ldw[18,8] ldw[18,9]
 ldw[18,10] ldw[18,11] ldw[18,12] ldw[18,13] ldw[18,14] ldw[18,15]
 ldw[19,0] ldw[19,1] ldw[19,2] ldw[19,3] ldw[19,4] ldw[19,5] ldw[19,6]
 ldw[19,7] ldw[19,8] ldw[19,9] ldw[19,10] ldw[19,11] ldw[19,12] ldw[19,13]
 ldw[19,14] ldw[19,15] ldw[20,0] ldw[20,1] ldw[20,2] ldw[20,3] ldw[20,4]
 ldw[20,5] ldw[20,6] ldw[20,7] ldw[20,8] ldw[20,9] ldw[20,10] ldw[20,11]
 ldw[20,12] ldw[20,13] ldw[20,14] ldw[20,15] ldw[21,0] ldw[21,1] ldw[21,2]
 ldw[21,3] ldw[21,4] ldw[21,5] ldw[21,6] ldw[21,7] ldw[21,8] ldw[21,9]
 ldw[21,10] ldw[21,11] ldw[21,12] ldw[21,13] ldw[21,14] ldw[21,15]
 ldw[22,0] ldw[22,1] ldw[22,2] ldw[22,3] ldw[22,4] ldw[22,5] ldw[22,6]
 ldw[22,7] ldw[22,8] ldw[22,9] ldw[22,10] ldw[22,11] ldw[22,12] ldw[22,13]
 ldw[22,14] ldw[22,15] ldw[23,0] ldw[23,1] ldw[23,2] ldw[23,3] ldw[23,4]
 ldw[23,5] ldw[23,6] ldw[23,7] ldw[23,8] ldw[23,9] ldw[23,10] ldw[23,11]
 ldw[23,12] ldw[23,13] ldw[23,14] ldw[23,15] lstk[39,0] lstk[39,1]
 lstk[39,2] lstk[39,3] lstk[39,4] lstk[39,5] lstk[39,6] lstk[39,7]
 lstk[39,8] lstk[39,9] lstk[39,10] lstk[39,11] lstk[39,12] lstk[39,13]
 lstk[39,14] lstk[39,15] lstk[40,0] lstk[40,1] lstk[40,2] lstk[40,3]
 lstk[40,4] lstk[40,5] lstk[40,6] lstk[40,7] lstk[40,8] lstk[40,9]
 lstk[40,10] lstk[40,11] lstk[40,12] lstk[40,13] lstk[40,14] lstk[40,15]
 lstk[41,0] lstk[41,1] lstk[41,2] lstk[41,3] lstk[41,4] lstk[41,5]
 lstk[41,6] lstk[41,7] lstk[41,8] lstk[41,9] lstk[41,10] lstk[41,11]
 lstk[41,12] lstk[41,13] lstk[41,14] lstk[41,15] lstk[42,0] lstk[42,1]
 lstk[42,2] lstk[42,3] lstk[42,4] lstk[42,5] lstk[42,6] lstk[42,7]
 lstk[42,8] lstk[42,9] lstk[42,10] lstk[42,11] lstk[42,12] lstk[42,13]
 lstk[42,14] lstk[42,15] lstk[43,0] lstk[43,1] lstk[43,2] lstk[43,3]
 lstk[43,4] lstk[43,5] lstk[43,6] lstk[43,7] lstk[43,8] lstk[43,9]
 lstk[43,10] lstk[43,11] lstk[43,12] lstk[43,13] lstk[43,14] lstk[43,15]
 lstk[44,0] lstk[44,1] lstk[44,2] lstk[44,3] lstk[44,4] lstk[44,5]
 lstk[44,6] lstk[44,7] lstk[44,8] lstk[44,9] lstk[44,10] lstk[44,11]
 lstk[44,12] lstk[44,13] lstk[44,14] lstk[44,15] lstk[45,0] lstk[45,1]
 lstk[45,2] lstk[45,3] lstk[45,4] lstk[45,5] lstk[45,6] lstk[45,7]
 lstk[45,8] lstk[45,9] lstk[45,10] lstk[45,11] lstk[45,12] lstk[45,13]
 lstk[45,14] lstk[45,15] lstk[46,0] lstk[46,1] lstk[46,2] lstk[46,3]
 lstk[46,4] lstk[46,5] lstk[46,6] lstk[46,7] lstk[46,8] lstk[46,9]
 lstk[46,10] lstk[46,11] lstk[46,12] lstk[46,13] lstk[46,14] lstk[46,15]
 lstk[47,0] lstk[47,1] lstk[47,2] lstk[47,3] lstk[47,4] lstk[47,5]
 lstk[47,6] lstk[47,7] lstk[47,8] lstk[47,9] lstk[47,10] lstk[47,11]
 lstk[47,12] lstk[47,13] lstk[47,14] lstk[47,15] lstk[48,0] lstk[48,1]
 lstk[48,2] lstk[48,3] lstk[48,4] lstk[48,5] lstk[48,6] lstk[48,7]
 lstk[48,8] lstk[48,9] lstk[48,10] lstk[48,11] lstk[48,12] lstk[48,13]
 lstk[48,14] lstk[48,15] lstk[49,0] lstk[49,1] lstk[49,2] lstk[49,3]
 lstk[49,4] lstk[49,5] lstk[49,6] lstk[49,7] lstk[49,8] lstk[49,9]
 lstk[49,10] lstk[49,11] lstk[49,12] lstk[49,13] lstk[49,14] lstk[49,15]
 lstk[50,0] lstk[50,1] lstk[50,2] lstk[50,3] lstk[50,4] lstk[50,5]
 lstk[50,6] lstk[50,7] lstk[50,8] lstk[50,9] lstk[50,10] lstk[50,11]
 lstk[50,12] lstk[50,13] lstk[50,14] lstk[50,15] lLANE[0] lLANE[1] lLANE[2]
 lLANE[3] lLANE[4] lLANE[5] lLANE[6] lLANE[7] lLANE[8] lLANE[9] lLANE[10]
 lLANE[11] lLANE[12] lLANE[13] lLANE[14] lLANE[15] ugEQK[0,0] ugEQK[0,1]
 ugEQK[0,2] ugEQK[0,3] ugEQK[0,4] ugEQK[0,5] ugEQK[0,6] ugEQK[0,7]
 ugEQK[0,8] ugEQK[0,9] ugEQK[0,10] ugEQK[0,11] ugEQK[0,12] ugEQK[0,13]
 ugEQK[0,14] ugEQK[0,15] ugEQK[1,0] ugEQK[1,1] ugEQK[1,2] ugEQK[1,3]
 ugEQK[1,4] ugEQK[1,5] ugEQK[1,6] ugEQK[1,7] ugEQK[1,8] ugEQK[1,9]
 ugEQK[1,10] ugEQK[1,11] ugEQK[1,12] ugEQK[1,13] ugEQK[1,14] ugEQK[1,15]
 ugEQK[2,0] ugEQK[2,1] ugEQK[2,2] ugEQK[2,3] ugEQK[2,4] ugEQK[2,5]
 ugEQK[2,6] ugEQK[2,7] ugEQK[2,8] ugEQK[2,9] ugEQK[2,10] ugEQK[2,11]
 ugEQK[2,12] ugEQK[2,13] ugEQK[2,14] ugEQK[2,15] ugEQK[3,0] ugEQK[3,1]
 ugEQK[3,2] ugEQK[3,3] ugEQK[3,4] ugEQK[3,5] ugEQK[3,6] ugEQK[3,7]
 ugEQK[3,8] ugEQK[3,9] ugEQK[3,10] ugEQK[3,11] ugEQK[3,12] ugEQK[3,13]
 ugEQK[3,14] ugEQK[3,15] uGstk[0,0] uGstk[0,1] uGstk[0,2] uGstk[0,3]
 uGstk[0,4] uGstk[0,5] uGstk[0,6] uGstk[0,7] uGstk[0,8] uGstk[0,9]
 uGstk[0,10] uGstk[0,11] uGstk[0,12] uGstk[0,13] uGstk[0,14] uGstk[0,15]
 uGstk[1,0] uGstk[1,1] uGstk[1,2] uGstk[1,3] uGstk[1,4] uGstk[1,5]
 uGstk[1,6] uGstk[1,7] uGstk[1,8] uGstk[1,9] uGstk[1,10] uGstk[1,11]
 uGstk[1,12] uGstk[1,13] uGstk[1,14] uGstk[1,15] uGstk[2,0] uGstk[2,1]
 uGstk[2,2] uGstk[2,3] uGstk[2,4] uGstk[2,5] uGstk[2,6] uGstk[2,7]
 uGstk[2,8] uGstk[2,9] uGstk[2,10] uGstk[2,11] uGstk[2,12] uGstk[2,13]
 uGstk[2,14] uGstk[2,15] uGstk[3,0] uGstk[3,1] uGstk[3,2] uGstk[3,3]
 uGstk[3,4] uGstk[3,5] uGstk[3,6] uGstk[3,7] uGstk[3,8] uGstk[3,9]
 uGstk[3,10] uGstk[3,11] uGstk[3,12] uGstk[3,13] uGstk[3,14] uGstk[3,15]
 uDetX[1,0] uDetX[1,1] uDetX[1,2] uDetX[1,3] uDetX[1,4] uDetX[1,5]
 uDetX[1,6] uDetX[1,7] uDetX[1,8] uDetX[1,9] uDetX[1,10] uDetX[1,11]
 uDetX[1,12] uDetX[1,13] uDetX[1,14] uDetX[1,15] uDetX[2,0] uDetX[2,1]
 uDetX[2,2] uDetX[2,3] uDetX[2,4] uDetX[2,5] uDetX[2,6] uDetX[2,7]
 uDetX[2,8] uDetX[2,9] uDetX[2,10] uDetX[2,11] uDetX[2,12] uDetX[2,13]
 uDetX[2,14] uDetX[2,15] uDetX[3,0] uDetX[3,1] uDetX[3,2] uDetX[3,3]
 uDetX[3,4] uDetX[3,5] uDetX[3,6] uDetX[3,7] uDetX[3,8] uDetX[3,9]
 uDetX[3,10] uDetX[3,11] uDetX[3,12] uDetX[3,13] uDetX[3,14] uDetX[3,15]
 uDetX[4,0] uDetX[4,1] uDetX[4,2] uDetX[4,3] uDetX[4,4] uDetX[4,5]
 uDetX[4,6] uDetX[4,7] uDetX[4,8] uDetX[4,9] uDetX[4,10] uDetX[4,11]
 uDetX[4,12] uDetX[4,13] uDetX[4,14] uDetX[4,15] uDetY[1,0] uDetY[1,1]
 uDetY[1,2] uDetY[1,3] uDetY[1,4] uDetY[1,5] uDetY[1,6] uDetY[1,7]
 uDetY[1,8] uDetY[1,9] uDetY[1,10] uDetY[1,11] uDetY[1,12] uDetY[1,13]
 uDetY[1,14] uDetY[1,15] uDetY[2,0] uDetY[2,1] uDetY[2,2] uDetY[2,3]
 uDetY[2,4] uDetY[2,5] uDetY[2,6] uDetY[2,7] uDetY[2,8] uDetY[2,9]
 uDetY[2,10] uDetY[2,11] uDetY[2,12] uDetY[2,13] uDetY[2,14] uDetY[2,15]
 uDetY[3,0] uDetY[3,1] uDetY[3,2] uDetY[3,3] uDetY[3,4] uDetY[3,5]
 uDetY[3,6] uDetY[3,7] uDetY[3,8] uDetY[3,9] uDetY[3,10] uDetY[3,11]
 uDetY[3,12] uDetY[3,13] uDetY[3,14] uDetY[3,15] uDetY[4,0] uDetY[4,1]
 uDetY[4,2] uDetY[4,3] uDetY[4,4] uDetY[4,5] uDetY[4,6] uDetY[4,7]
 uDetY[4,8] uDetY[4,9] uDetY[4,10] uDetY[4,11] uDetY[4,12] uDetY[4,13]
 uDetY[4,14] uDetY[4,15] uDetZ[1,0] uDetZ[1,1] uDetZ[1,2] uDetZ[1,3]
 uDetZ[1,4] uDetZ[1,5] uDetZ[1,6] uDetZ[1,7] uDetZ[1,8] uDetZ[1,9]
 uDetZ[1,10] uDetZ[1,11] uDetZ[1,12] uDetZ[1,13] uDetZ[1,14] uDetZ[1,15]
 uDetZ[2,0] uDetZ[2,1] uDetZ[2,2] uDetZ[2,3] uDetZ[2,4] uDetZ[2,5]
 uDetZ[2,6] uDetZ[2,7] uDetZ[2,8] uDetZ[2,9] uDetZ[2,10] uDetZ[2,11]
 uDetZ[2,12] uDetZ[2,13] uDetZ[2,14] uDetZ[2,15] uDetZ[3,0] uDetZ[3,1]
 uDetZ[3,2] uDetZ[3,3] uDetZ[3,4] uDetZ[3,5] uDetZ[3,6] uDetZ[3,7]
 uDetZ[3,8] uDetZ[3,9] uDetZ[3,10] uDetZ[3,11] uDetZ[3,12] uDetZ[3,13]
 uDetZ[3,14] uDetZ[3,15] uDetW[0,0] uDetW[0,1] uDetW[0,2] uDetW[0,3]
 uDetW[0,4] uDetW[0,5] uDetW[0,6] uDetW[0,7] uDetW[0,8] uDetW[0,9]
 uDetW[0,10] uDetW[0,11] uDetW[0,12] uDetW[0,13] uDetW[0,14] uDetW[0,15]
 uDetW[1,0] uDetW[1,1] uDetW[1,2] uDetW[1,3] uDetW[1,4] uDetW[1,5]
 uDetW[1,6] uDetW[1,7] uDetW[1,8] uDetW[1,9] uDetW[1,10] uDetW[1,11]
 uDetW[1,12] uDetW[1,13] uDetW[1,14] uDetW[1,15] uDetW[2,0] uDetW[2,1]
 uDetW[2,2] uDetW[2,3] uDetW[2,4] uDetW[2,5] uDetW[2,6] uDetW[2,7]
 uDetW[2,8] uDetW[2,9] uDetW[2,10] uDetW[2,11] uDetW[2,12] uDetW[2,13]
 uDetW[2,14] uDetW[2,15] uDetW[3,0] uDetW[3,1] uDetW[3,2] uDetW[3,3]
 uDetW[3,4] uDetW[3,5] uDetW[3,6] uDetW[3,7] uDetW[3,8] uDetW[3,9]
 uDetW[3,10] uDetW[3,11] uDetW[3,12] uDetW[3,13] uDetW[3,14] uDetW[3,15]
 uFrSC[1,0] uFrSC[1,1] uFrSC[1,2] uFrSC[1,3] uFrSC[1,4] uFrSC[1,5]
 uFrSC[1,6] uFrSC[1,7] uFrSC[1,8] uFrSC[1,9] uFrSC[1,10] uFrSC[1,11]
 uFrSC[1,12] uFrSC[1,13] uFrSC[1,14] uFrSC[1,15] uFrSC[2,0] uFrSC[2,1]
 uFrSC[2,2] uFrSC[2,3] uFrSC[2,4] uFrSC[2,5] uFrSC[2,6] uFrSC[2,7]
 uFrSC[2,8] uFrSC[2,9] uFrSC[2,10] uFrSC[2,11] uFrSC[2,12] uFrSC[2,13]
 uFrSC[2,14] uFrSC[2,15] uFrSC[3,0] uFrSC[3,1] uFrSC[3,2] uFrSC[3,3]
 uFrSC[3,4] uFrSC[3,5] uFrSC[3,6] uFrSC[3,7] uFrSC[3,8] uFrSC[3,9]
 uFrSC[3,10] uFrSC[3,11] uFrSC[3,12] uFrSC[3,13] uFrSC[3,14] uFrSC[3,15]
 uFrSC[4,0] uFrSC[4,1] uFrSC[4,2] uFrSC[4,3] uFrSC[4,4] uFrSC[4,5]
 uFrSC[4,6] uFrSC[4,7] uFrSC[4,8] uFrSC[4,9] uFrSC[4,10] uFrSC[4,11]
 uFrSC[4,12] uFrSC[4,13] uFrSC[4,14] uFrSC[4,15] uFrMC[2,0] uFrMC[2,1]
 uFrMC[2,2] uFrMC[2,3] uFrMC[2,4] uFrMC[2,5] uFrMC[2,6] uFrMC[2,7]
 uFrMC[2,8] uFrMC[2,9] uFrMC[2,10] uFrMC[2,11] uFrMC[2,12] uFrMC[2,13]
 uFrMC[2,14] uFrMC[2,15] uFrMC[3,0] uFrMC[3,1] uFrMC[3,2] uFrMC[3,3]
 uFrMC[3,4] uFrMC[3,5] uFrMC[3,6] uFrMC[3,7] uFrMC[3,8] uFrMC[3,9]
 uFrMC[3,10] uFrMC[3,11] uFrMC[3,12] uFrMC[3,13] uFrMC[3,14] uFrMC[3,15]
 uFrMC[4,0] uFrMC[4,1] uFrMC[4,2] uFrMC[4,3] uFrMC[4,4] uFrMC[4,5]
 uFrMC[4,6] uFrMC[4,7] uFrMC[4,8] uFrMC[4,9] uFrMC[4,10] uFrMC[4,11]
 uFrMC[4,12] uFrMC[4,13] uFrMC[4,14] uFrMC[4,15] ugeqka[0,0] ugeqka[0,1]
 ugeqka[0,2] ugeqka[0,3] ugeqka[0,4] ugeqka[0,5] ugeqka[0,6] ugeqka[0,7]
 ugeqka[0,8] ugeqka[0,9] ugeqka[0,10] ugeqka[0,11] ugeqka[0,12]
 ugeqka[0,13] ugeqka[0,14] ugeqka[0,15] ugeqka[1,0] ugeqka[1,1] ugeqka[1,2]
 ugeqka[1,3] ugeqka[1,4] ugeqka[1,5] ugeqka[1,6] ugeqka[1,7] ugeqka[1,8]
 ugeqka[1,9] ugeqka[1,10] ugeqka[1,11] ugeqka[1,12] ugeqka[1,13]
 ugeqka[1,14] ugeqka[1,15] ugeqka[2,0] ugeqka[2,1] ugeqka[2,2] ugeqka[2,3]
 ugeqka[2,4] ugeqka[2,5] ugeqka[2,6] ugeqka[2,7] ugeqka[2,8] ugeqka[2,9]
 ugeqka[2,10] ugeqka[2,11] ugeqka[2,12] ugeqka[2,13] ugeqka[2,14]
 ugeqka[2,15] ugeqka[3,0] ugeqka[3,1] ugeqka[3,2] ugeqka[3,3] ugeqka[3,4]
 ugeqka[3,5] ugeqka[3,6] ugeqka[3,7] ugeqka[3,8] ugeqka[3,9] ugeqka[3,10]
 ugeqka[3,11] ugeqka[3,12] ugeqka[3,13] ugeqka[3,14] ugeqka[3,15]
 ugstka[0,0] ugstka[0,1] ugstka[0,2] ugstka[0,3] ugstka[0,4] ugstka[0,5]
 ugstka[0,6] ugstka[0,7] ugstka[0,8] ugstka[0,9] ugstka[0,10] ugstka[0,11]
 ugstka[0,12] ugstka[0,13] ugstka[0,14] ugstka[0,15] ugstka[1,0]
 ugstka[1,1] ugstka[1,2] ugstka[1,3] ugstka[1,4] ugstka[1,5] ugstka[1,6]
 ugstka[1,7] ugstka[1,8] ugstka[1,9] ugstka[1,10] ugstka[1,11] ugstka[1,12]
 ugstka[1,13] ugstka[1,14] ugstka[1,15] ugstka[2,0] ugstka[2,1] ugstka[2,2]
 ugstka[2,3] ugstka[2,4] ugstka[2,5] ugstka[2,6] ugstka[2,7] ugstka[2,8]
 ugstka[2,9] ugstka[2,10] ugstka[2,11] ugstka[2,12] ugstka[2,13]
 ugstka[2,14] ugstka[2,15] ugstka[3,0] ugstka[3,1] ugstka[3,2] ugstka[3,3]
 ugstka[3,4] ugstka[3,5] ugstka[3,6] ugstka[3,7] ugstka[3,8] ugstka[3,9]
 ugstka[3,10] ugstka[3,11] ugstka[3,12] ugstka[3,13] ugstka[3,14]
 ugstka[3,15] udetXa[1,0] udetXa[1,1] udetXa[1,2] udetXa[1,3] udetXa[1,4]
 udetXa[1,5] udetXa[1,6] udetXa[1,7] udetXa[1,8] udetXa[1,9] udetXa[1,10]
 udetXa[1,11] udetXa[1,12] udetXa[1,13] udetXa[1,14] udetXa[1,15]
 udetXa[2,0] udetXa[2,1] udetXa[2,2] udetXa[2,3] udetXa[2,4] udetXa[2,5]
 udetXa[2,6] udetXa[2,7] udetXa[2,8] udetXa[2,9] udetXa[2,10] udetXa[2,11]
 udetXa[2,12] udetXa[2,13] udetXa[2,14] udetXa[2,15] udetXa[3,0]
 udetXa[3,1] udetXa[3,2] udetXa[3,3] udetXa[3,4] udetXa[3,5] udetXa[3,6]
 udetXa[3,7] udetXa[3,8] udetXa[3,9] udetXa[3,10] udetXa[3,11] udetXa[3,12]
 udetXa[3,13] udetXa[3,14] udetXa[3,15] udetXa[4,0] udetXa[4,1] udetXa[4,2]
 udetXa[4,3] udetXa[4,4] udetXa[4,5] udetXa[4,6] udetXa[4,7] udetXa[4,8]
 udetXa[4,9] udetXa[4,10] udetXa[4,11] udetXa[4,12] udetXa[4,13]
 udetXa[4,14] udetXa[4,15] udetYa[1,0] udetYa[1,1] udetYa[1,2] udetYa[1,3]
 udetYa[1,4] udetYa[1,5] udetYa[1,6] udetYa[1,7] udetYa[1,8] udetYa[1,9]
 udetYa[1,10] udetYa[1,11] udetYa[1,12] udetYa[1,13] udetYa[1,14]
 udetYa[1,15] udetYa[2,0] udetYa[2,1] udetYa[2,2] udetYa[2,3] udetYa[2,4]
 udetYa[2,5] udetYa[2,6] udetYa[2,7] udetYa[2,8] udetYa[2,9] udetYa[2,10]
 udetYa[2,11] udetYa[2,12] udetYa[2,13] udetYa[2,14] udetYa[2,15]
 udetYa[3,0] udetYa[3,1] udetYa[3,2] udetYa[3,3] udetYa[3,4] udetYa[3,5]
 udetYa[3,6] udetYa[3,7] udetYa[3,8] udetYa[3,9] udetYa[3,10] udetYa[3,11]
 udetYa[3,12] udetYa[3,13] udetYa[3,14] udetYa[3,15] udetYa[4,0]
 udetYa[4,1] udetYa[4,2] udetYa[4,3] udetYa[4,4] udetYa[4,5] udetYa[4,6]
 udetYa[4,7] udetYa[4,8] udetYa[4,9] udetYa[4,10] udetYa[4,11] udetYa[4,12]
 udetYa[4,13] udetYa[4,14] udetYa[4,15] udetZa[1,0] udetZa[1,1] udetZa[1,2]
 udetZa[1,3] udetZa[1,4] udetZa[1,5] udetZa[1,6] udetZa[1,7] udetZa[1,8]
 udetZa[1,9] udetZa[1,10] udetZa[1,11] udetZa[1,12] udetZa[1,13]
 udetZa[1,14] udetZa[1,15] udetZa[2,0] udetZa[2,1] udetZa[2,2] udetZa[2,3]
 udetZa[2,4] udetZa[2,5] udetZa[2,6] udetZa[2,7] udetZa[2,8] udetZa[2,9]
 udetZa[2,10] udetZa[2,11] udetZa[2,12] udetZa[2,13] udetZa[2,14]
 udetZa[2,15] udetZa[3,0] udetZa[3,1] udetZa[3,2] udetZa[3,3] udetZa[3,4]
 udetZa[3,5] udetZa[3,6] udetZa[3,7] udetZa[3,8] udetZa[3,9] udetZa[3,10]
 udetZa[3,11] udetZa[3,12] udetZa[3,13] udetZa[3,14] udetZa[3,15]
 udetWa[0,0] udetWa[0,1] udetWa[0,2] udetWa[0,3] udetWa[0,4] udetWa[0,5]
 udetWa[0,6] udetWa[0,7] udetWa[0,8] udetWa[0,9] udetWa[0,10] udetWa[0,11]
 udetWa[0,12] udetWa[0,13] udetWa[0,14] udetWa[0,15] udetWa[1,0]
 udetWa[1,1] udetWa[1,2] udetWa[1,3] udetWa[1,4] udetWa[1,5] udetWa[1,6]
 udetWa[1,7] udetWa[1,8] udetWa[1,9] udetWa[1,10] udetWa[1,11] udetWa[1,12]
 udetWa[1,13] udetWa[1,14] udetWa[1,15] udetWa[2,0] udetWa[2,1] udetWa[2,2]
 udetWa[2,3] udetWa[2,4] udetWa[2,5] udetWa[2,6] udetWa[2,7] udetWa[2,8]
 udetWa[2,9] udetWa[2,10] udetWa[2,11] udetWa[2,12] udetWa[2,13]
 udetWa[2,14] udetWa[2,15] udetWa[3,0] udetWa[3,1] udetWa[3,2] udetWa[3,3]
 udetWa[3,4] udetWa[3,5] udetWa[3,6] udetWa[3,7] udetWa[3,8] udetWa[3,9]
 udetWa[3,10] udetWa[3,11] udetWa[3,12] udetWa[3,13] udetWa[3,14]
 udetWa[3,15] ufrSCa[1,0] ufrSCa[1,1] ufrSCa[1,2] ufrSCa[1,3] ufrSCa[1,4]
 ufrSCa[1,5] ufrSCa[1,6] ufrSCa[1,7] ufrSCa[1,8] ufrSCa[1,9] ufrSCa[1,10]
 ufrSCa[1,11] ufrSCa[1,12] ufrSCa[1,13] ufrSCa[1,14] ufrSCa[1,15]
 ufrSCa[2,0] ufrSCa[2,1] ufrSCa[2,2] ufrSCa[2,3] ufrSCa[2,4] ufrSCa[2,5]
 ufrSCa[2,6] ufrSCa[2,7] ufrSCa[2,8] ufrSCa[2,9] ufrSCa[2,10] ufrSCa[2,11]
 ufrSCa[2,12] ufrSCa[2,13] ufrSCa[2,14] ufrSCa[2,15] ufrSCa[3,0]
 ufrSCa[3,1] ufrSCa[3,2] ufrSCa[3,3] ufrSCa[3,4] ufrSCa[3,5] ufrSCa[3,6]
 ufrSCa[3,7] ufrSCa[3,8] ufrSCa[3,9] ufrSCa[3,10] ufrSCa[3,11] ufrSCa[3,12]
 ufrSCa[3,13] ufrSCa[3,14] ufrSCa[3,15] ufrSCa[4,0] ufrSCa[4,1] ufrSCa[4,2]
 ufrSCa[4,3] ufrSCa[4,4] ufrSCa[4,5] ufrSCa[4,6] ufrSCa[4,7] ufrSCa[4,8]
 ufrSCa[4,9] ufrSCa[4,10] ufrSCa[4,11] ufrSCa[4,12] ufrSCa[4,13]
 ufrSCa[4,14] ufrSCa[4,15] ufrMCa[2,0] ufrMCa[2,1] ufrMCa[2,2] ufrMCa[2,3]
 ufrMCa[2,4] ufrMCa[2,5] ufrMCa[2,6] ufrMCa[2,7] ufrMCa[2,8] ufrMCa[2,9]
 ufrMCa[2,10] ufrMCa[2,11] ufrMCa[2,12] ufrMCa[2,13] ufrMCa[2,14]
 ufrMCa[2,15] ufrMCa[3,0] ufrMCa[3,1] ufrMCa[3,2] ufrMCa[3,3] ufrMCa[3,4]
 ufrMCa[3,5] ufrMCa[3,6] ufrMCa[3,7] ufrMCa[3,8] ufrMCa[3,9] ufrMCa[3,10]
 ufrMCa[3,11] ufrMCa[3,12] ufrMCa[3,13] ufrMCa[3,14] ufrMCa[3,15]
 ufrMCa[4,0] ufrMCa[4,1] ufrMCa[4,2] ufrMCa[4,3] ufrMCa[4,4] ufrMCa[4,5]
 ufrMCa[4,6] ufrMCa[4,7] ufrMCa[4,8] ufrMCa[4,9] ufrMCa[4,10] ufrMCa[4,11]
 ufrMCa[4,12] ufrMCa[4,13] ufrMCa[4,14] ufrMCa[4,15] lGstk[46,0]
 lGstk[46,1] lGstk[46,2] lGstk[46,3] lGstk[46,4] lGstk[46,5] lGstk[46,6]
 lGstk[46,7] lGstk[46,8] lGstk[46,9] lGstk[46,10] lGstk[46,11] lGstk[46,12]
 lGstk[46,13] lGstk[46,14] lGstk[46,15] lGstk[47,0] lGstk[47,1] lGstk[47,2]
 lGstk[47,3] lGstk[47,4] lGstk[47,5] lGstk[47,6] lGstk[47,7] lGstk[47,8]
 lGstk[47,9] lGstk[47,10] lGstk[47,11] lGstk[47,12] lGstk[47,13]
 lGstk[47,14] lGstk[47,15] lGstk[48,0] lGstk[48,1] lGstk[48,2] lGstk[48,3]
 lGstk[48,4] lGstk[48,5] lGstk[48,6] lGstk[48,7] lGstk[48,8] lGstk[48,9]
 lGstk[48,10] lGstk[48,11] lGstk[48,12] lGstk[48,13] lGstk[48,14]
 lGstk[48,15] lGstk[49,0] lGstk[49,1] lGstk[49,2] lGstk[49,3] lGstk[49,4]
 lGstk[49,5] lGstk[49,6] lGstk[49,7] lGstk[49,8] lGstk[49,9] lGstk[49,10]
 lGstk[49,11] lGstk[49,12] lGstk[49,13] lGstk[49,14] lGstk[49,15]
 lGstk[50,0] lGstk[50,1] lGstk[50,2] lGstk[50,3] lGstk[50,4] lGstk[50,5]
 lGstk[50,6] lGstk[50,7] lGstk[50,8] lGstk[50,9] lGstk[50,10] lGstk[50,11]
 lGstk[50,12] lGstk[50,13] lGstk[50,14] lGstk[50,15] lDetX[19,0]
 lDetX[19,1] lDetX[19,2] lDetX[19,3] lDetX[19,4] lDetX[19,5] lDetX[19,6]
 lDetX[19,7] lDetX[19,8] lDetX[19,9] lDetX[19,10] lDetX[19,11] lDetX[19,12]
 lDetX[19,13] lDetX[19,14] lDetX[19,15] lDetX[20,0] lDetX[20,1] lDetX[20,2]
 lDetX[20,3] lDetX[20,4] lDetX[20,5] lDetX[20,6] lDetX[20,7] lDetX[20,8]
 lDetX[20,9] lDetX[20,10] lDetX[20,11] lDetX[20,12] lDetX[20,13]
 lDetX[20,14] lDetX[20,15] lDetX[21,0] lDetX[21,1] lDetX[21,2] lDetX[21,3]
 lDetX[21,4] lDetX[21,5] lDetX[21,6] lDetX[21,7] lDetX[21,8] lDetX[21,9]
 lDetX[21,10] lDetX[21,11] lDetX[21,12] lDetX[21,13] lDetX[21,14]
 lDetX[21,15] lDetX[22,0] lDetX[22,1] lDetX[22,2] lDetX[22,3] lDetX[22,4]
 lDetX[22,5] lDetX[22,6] lDetX[22,7] lDetX[22,8] lDetX[22,9] lDetX[22,10]
 lDetX[22,11] lDetX[22,12] lDetX[22,13] lDetX[22,14] lDetX[22,15]
 lDetX[23,0] lDetX[23,1] lDetX[23,2] lDetX[23,3] lDetX[23,4] lDetX[23,5]
 lDetX[23,6] lDetX[23,7] lDetX[23,8] lDetX[23,9] lDetX[23,10] lDetX[23,11]
 lDetX[23,12] lDetX[23,13] lDetX[23,14] lDetX[23,15] lDetY[19,0]
 lDetY[19,1] lDetY[19,2] lDetY[19,3] lDetY[19,4] lDetY[19,5] lDetY[19,6]
 lDetY[19,7] lDetY[19,8] lDetY[19,9] lDetY[19,10] lDetY[19,11] lDetY[19,12]
 lDetY[19,13] lDetY[19,14] lDetY[19,15] lDetY[20,0] lDetY[20,1] lDetY[20,2]
 lDetY[20,3] lDetY[20,4] lDetY[20,5] lDetY[20,6] lDetY[20,7] lDetY[20,8]
 lDetY[20,9] lDetY[20,10] lDetY[20,11] lDetY[20,12] lDetY[20,13]
 lDetY[20,14] lDetY[20,15] lDetY[21,0] lDetY[21,1] lDetY[21,2] lDetY[21,3]
 lDetY[21,4] lDetY[21,5] lDetY[21,6] lDetY[21,7] lDetY[21,8] lDetY[21,9]
 lDetY[21,10] lDetY[21,11] lDetY[21,12] lDetY[21,13] lDetY[21,14]
 lDetY[21,15] lDetY[22,0] lDetY[22,1] lDetY[22,2] lDetY[22,3] lDetY[22,4]
 lDetY[22,5] lDetY[22,6] lDetY[22,7] lDetY[22,8] lDetY[22,9] lDetY[22,10]
 lDetY[22,11] lDetY[22,12] lDetY[22,13] lDetY[22,14] lDetY[22,15]
 lDetY[23,0] lDetY[23,1] lDetY[23,2] lDetY[23,3] lDetY[23,4] lDetY[23,5]
 lDetY[23,6] lDetY[23,7] lDetY[23,8] lDetY[23,9] lDetY[23,10] lDetY[23,11]
 lDetY[23,12] lDetY[23,13] lDetY[23,14] lDetY[23,15] lDetZ[19,0]
 lDetZ[19,1] lDetZ[19,2] lDetZ[19,3] lDetZ[19,4] lDetZ[19,5] lDetZ[19,6]
 lDetZ[19,7] lDetZ[19,8] lDetZ[19,9] lDetZ[19,10] lDetZ[19,11] lDetZ[19,12]
 lDetZ[19,13] lDetZ[19,14] lDetZ[19,15] lDetZ[20,0] lDetZ[20,1] lDetZ[20,2]
 lDetZ[20,3] lDetZ[20,4] lDetZ[20,5] lDetZ[20,6] lDetZ[20,7] lDetZ[20,8]
 lDetZ[20,9] lDetZ[20,10] lDetZ[20,11] lDetZ[20,12] lDetZ[20,13]
 lDetZ[20,14] lDetZ[20,15] lDetZ[21,0] lDetZ[21,1] lDetZ[21,2] lDetZ[21,3]
 lDetZ[21,4] lDetZ[21,5] lDetZ[21,6] lDetZ[21,7] lDetZ[21,8] lDetZ[21,9]
 lDetZ[21,10] lDetZ[21,11] lDetZ[21,12] lDetZ[21,13] lDetZ[21,14]
 lDetZ[21,15] lDetZ[22,0] lDetZ[22,1] lDetZ[22,2] lDetZ[22,3] lDetZ[22,4]
 lDetZ[22,5] lDetZ[22,6] lDetZ[22,7] lDetZ[22,8] lDetZ[22,9] lDetZ[22,10]
 lDetZ[22,11] lDetZ[22,12] lDetZ[22,13] lDetZ[22,14] lDetZ[22,15]
 lDetZ[23,0] lDetZ[23,1] lDetZ[23,2] lDetZ[23,3] lDetZ[23,4] lDetZ[23,5]
 lDetZ[23,6] lDetZ[23,7] lDetZ[23,8] lDetZ[23,9] lDetZ[23,10] lDetZ[23,11]
 lDetZ[23,12] lDetZ[23,13] lDetZ[23,14] lDetZ[23,15] lDetW[19,0]
 lDetW[19,1] lDetW[19,2] lDetW[19,3] lDetW[19,4] lDetW[19,5] lDetW[19,6]
 lDetW[19,7] lDetW[19,8] lDetW[19,9] lDetW[19,10] lDetW[19,11] lDetW[19,12]
 lDetW[19,13] lDetW[19,14] lDetW[19,15] lDetW[20,0] lDetW[20,1] lDetW[20,2]
 lDetW[20,3] lDetW[20,4] lDetW[20,5] lDetW[20,6] lDetW[20,7] lDetW[20,8]
 lDetW[20,9] lDetW[20,10] lDetW[20,11] lDetW[20,12] lDetW[20,13]
 lDetW[20,14] lDetW[20,15] lDetW[21,0] lDetW[21,1] lDetW[21,2] lDetW[21,3]
 lDetW[21,4] lDetW[21,5] lDetW[21,6] lDetW[21,7] lDetW[21,8] lDetW[21,9]
 lDetW[21,10] lDetW[21,11] lDetW[21,12] lDetW[21,13] lDetW[21,14]
 lDetW[21,15] lDetW[22,0] lDetW[22,1] lDetW[22,2] lDetW[22,3] lDetW[22,4]
 lDetW[22,5] lDetW[22,6] lDetW[22,7] lDetW[22,8] lDetW[22,9] lDetW[22,10]
 lDetW[22,11] lDetW[22,12] lDetW[22,13] lDetW[22,14] lDetW[22,15]
 lFrSC[19,0] lFrSC[19,1] lFrSC[19,2] lFrSC[19,3] lFrSC[19,4] lFrSC[19,5]
 lFrSC[19,6] lFrSC[19,7] lFrSC[19,8] lFrSC[19,9] lFrSC[19,10] lFrSC[19,11]
 lFrSC[19,12] lFrSC[19,13] lFrSC[19,14] lFrSC[19,15] lFrSC[20,0]
 lFrSC[20,1] lFrSC[20,2] lFrSC[20,3] lFrSC[20,4] lFrSC[20,5] lFrSC[20,6]
 lFrSC[20,7] lFrSC[20,8] lFrSC[20,9] lFrSC[20,10] lFrSC[20,11] lFrSC[20,12]
 lFrSC[20,13] lFrSC[20,14] lFrSC[20,15] lFrSC[21,0] lFrSC[21,1] lFrSC[21,2]
 lFrSC[21,3] lFrSC[21,4] lFrSC[21,5] lFrSC[21,6] lFrSC[21,7] lFrSC[21,8]
 lFrSC[21,9] lFrSC[21,10] lFrSC[21,11] lFrSC[21,12] lFrSC[21,13]
 lFrSC[21,14] lFrSC[21,15] lFrSC[22,0] lFrSC[22,1] lFrSC[22,2] lFrSC[22,3]
 lFrSC[22,4] lFrSC[22,5] lFrSC[22,6] lFrSC[22,7] lFrSC[22,8] lFrSC[22,9]
 lFrSC[22,10] lFrSC[22,11] lFrSC[22,12] lFrSC[22,13] lFrSC[22,14]
 lFrSC[22,15] lFrSC[23,0] lFrSC[23,1] lFrSC[23,2] lFrSC[23,3] lFrSC[23,4]
 lFrSC[23,5] lFrSC[23,6] lFrSC[23,7] lFrSC[23,8] lFrSC[23,9] lFrSC[23,10]
 lFrSC[23,11] lFrSC[23,12] lFrSC[23,13] lFrSC[23,14] lFrSC[23,15]
 lFrMC[19,0] lFrMC[19,1] lFrMC[19,2] lFrMC[19,3] lFrMC[19,4] lFrMC[19,5]
 lFrMC[19,6] lFrMC[19,7] lFrMC[19,8] lFrMC[19,9] lFrMC[19,10] lFrMC[19,11]
 lFrMC[19,12] lFrMC[19,13] lFrMC[19,14] lFrMC[19,15] lFrMC[20,0]
 lFrMC[20,1] lFrMC[20,2] lFrMC[20,3] lFrMC[20,4] lFrMC[20,5] lFrMC[20,6]
 lFrMC[20,7] lFrMC[20,8] lFrMC[20,9] lFrMC[20,10] lFrMC[20,11] lFrMC[20,12]
 lFrMC[20,13] lFrMC[20,14] lFrMC[20,15] lFrMC[21,0] lFrMC[21,1] lFrMC[21,2]
 lFrMC[21,3] lFrMC[21,4] lFrMC[21,5] lFrMC[21,6] lFrMC[21,7] lFrMC[21,8]
 lFrMC[21,9] lFrMC[21,10] lFrMC[21,11] lFrMC[21,12] lFrMC[21,13]
 lFrMC[21,14] lFrMC[21,15] lFrMC[22,0] lFrMC[22,1] lFrMC[22,2] lFrMC[22,3]
 lFrMC[22,4] lFrMC[22,5] lFrMC[22,6] lFrMC[22,7] lFrMC[22,8] lFrMC[22,9]
 lFrMC[22,10] lFrMC[22,11] lFrMC[22,12] lFrMC[22,13] lFrMC[22,14]
 lFrMC[22,15] lgstka[46,0] lgstka[46,1] lgstka[46,2] lgstka[46,3]
 lgstka[46,4] lgstka[46,5] lgstka[46,6] lgstka[46,7] lgstka[46,8]
 lgstka[46,9] lgstka[46,10] lgstka[46,11] lgstka[46,12] lgstka[46,13]
 lgstka[46,14] lgstka[46,15] lgstka[47,0] lgstka[47,1] lgstka[47,2]
 lgstka[47,3] lgstka[47,4] lgstka[47,5] lgstka[47,6] lgstka[47,7]
 lgstka[47,8] lgstka[47,9] lgstka[47,10] lgstka[47,11] lgstka[47,12]
 lgstka[47,13] lgstka[47,14] lgstka[47,15] lgstka[48,0] lgstka[48,1]
 lgstka[48,2] lgstka[48,3] lgstka[48,4] lgstka[48,5] lgstka[48,6]
 lgstka[48,7] lgstka[48,8] lgstka[48,9] lgstka[48,10] lgstka[48,11]
 lgstka[48,12] lgstka[48,13] lgstka[48,14] lgstka[48,15] lgstka[49,0]
 lgstka[49,1] lgstka[49,2] lgstka[49,3] lgstka[49,4] lgstka[49,5]
 lgstka[49,6] lgstka[49,7] lgstka[49,8] lgstka[49,9] lgstka[49,10]
 lgstka[49,11] lgstka[49,12] lgstka[49,13] lgstka[49,14] lgstka[49,15]
 lgstka[50,0] lgstka[50,1] lgstka[50,2] lgstka[50,3] lgstka[50,4]
 lgstka[50,5] lgstka[50,6] lgstka[50,7] lgstka[50,8] lgstka[50,9]
 lgstka[50,10] lgstka[50,11] lgstka[50,12] lgstka[50,13] lgstka[50,14]
 lgstka[50,15] ldetXa[19,0] ldetXa[19,1] ldetXa[19,2] ldetXa[19,3]
 ldetXa[19,4] ldetXa[19,5] ldetXa[19,6] ldetXa[19,7] ldetXa[19,8]
 ldetXa[19,9] ldetXa[19,10] ldetXa[19,11] ldetXa[19,12] ldetXa[19,13]
 ldetXa[19,14] ldetXa[19,15] ldetXa[20,0] ldetXa[20,1] ldetXa[20,2]
 ldetXa[20,3] ldetXa[20,4] ldetXa[20,5] ldetXa[20,6] ldetXa[20,7]
 ldetXa[20,8] ldetXa[20,9] ldetXa[20,10] ldetXa[20,11] ldetXa[20,12]
 ldetXa[20,13] ldetXa[20,14] ldetXa[20,15] ldetXa[21,0] ldetXa[21,1]
 ldetXa[21,2] ldetXa[21,3] ldetXa[21,4] ldetXa[21,5] ldetXa[21,6]
 ldetXa[21,7] ldetXa[21,8] ldetXa[21,9] ldetXa[21,10] ldetXa[21,11]
 ldetXa[21,12] ldetXa[21,13] ldetXa[21,14] ldetXa[21,15] ldetXa[22,0]
 ldetXa[22,1] ldetXa[22,2] ldetXa[22,3] ldetXa[22,4] ldetXa[22,5]
 ldetXa[22,6] ldetXa[22,7] ldetXa[22,8] ldetXa[22,9] ldetXa[22,10]
 ldetXa[22,11] ldetXa[22,12] ldetXa[22,13] ldetXa[22,14] ldetXa[22,15]
 ldetXa[23,0] ldetXa[23,1] ldetXa[23,2] ldetXa[23,3] ldetXa[23,4]
 ldetXa[23,5] ldetXa[23,6] ldetXa[23,7] ldetXa[23,8] ldetXa[23,9]
 ldetXa[23,10] ldetXa[23,11] ldetXa[23,12] ldetXa[23,13] ldetXa[23,14]
 ldetXa[23,15] ldetYa[19,0] ldetYa[19,1] ldetYa[19,2] ldetYa[19,3]
 ldetYa[19,4] ldetYa[19,5] ldetYa[19,6] ldetYa[19,7] ldetYa[19,8]
 ldetYa[19,9] ldetYa[19,10] ldetYa[19,11] ldetYa[19,12] ldetYa[19,13]
 ldetYa[19,14] ldetYa[19,15] ldetYa[20,0] ldetYa[20,1] ldetYa[20,2]
 ldetYa[20,3] ldetYa[20,4] ldetYa[20,5] ldetYa[20,6] ldetYa[20,7]
 ldetYa[20,8] ldetYa[20,9] ldetYa[20,10] ldetYa[20,11] ldetYa[20,12]
 ldetYa[20,13] ldetYa[20,14] ldetYa[20,15] ldetYa[21,0] ldetYa[21,1]
 ldetYa[21,2] ldetYa[21,3] ldetYa[21,4] ldetYa[21,5] ldetYa[21,6]
 ldetYa[21,7] ldetYa[21,8] ldetYa[21,9] ldetYa[21,10] ldetYa[21,11]
 ldetYa[21,12] ldetYa[21,13] ldetYa[21,14] ldetYa[21,15] ldetYa[22,0]
 ldetYa[22,1] ldetYa[22,2] ldetYa[22,3] ldetYa[22,4] ldetYa[22,5]
 ldetYa[22,6] ldetYa[22,7] ldetYa[22,8] ldetYa[22,9] ldetYa[22,10]
 ldetYa[22,11] ldetYa[22,12] ldetYa[22,13] ldetYa[22,14] ldetYa[22,15]
 ldetYa[23,0] ldetYa[23,1] ldetYa[23,2] ldetYa[23,3] ldetYa[23,4]
 ldetYa[23,5] ldetYa[23,6] ldetYa[23,7] ldetYa[23,8] ldetYa[23,9]
 ldetYa[23,10] ldetYa[23,11] ldetYa[23,12] ldetYa[23,13] ldetYa[23,14]
 ldetYa[23,15] ldetZa[19,0] ldetZa[19,1] ldetZa[19,2] ldetZa[19,3]
 ldetZa[19,4] ldetZa[19,5] ldetZa[19,6] ldetZa[19,7] ldetZa[19,8]
 ldetZa[19,9] ldetZa[19,10] ldetZa[19,11] ldetZa[19,12] ldetZa[19,13]
 ldetZa[19,14] ldetZa[19,15] ldetZa[20,0] ldetZa[20,1] ldetZa[20,2]
 ldetZa[20,3] ldetZa[20,4] ldetZa[20,5] ldetZa[20,6] ldetZa[20,7]
 ldetZa[20,8] ldetZa[20,9] ldetZa[20,10] ldetZa[20,11] ldetZa[20,12]
 ldetZa[20,13] ldetZa[20,14] ldetZa[20,15] ldetZa[21,0] ldetZa[21,1]
 ldetZa[21,2] ldetZa[21,3] ldetZa[21,4] ldetZa[21,5] ldetZa[21,6]
 ldetZa[21,7] ldetZa[21,8] ldetZa[21,9] ldetZa[21,10] ldetZa[21,11]
 ldetZa[21,12] ldetZa[21,13] ldetZa[21,14] ldetZa[21,15] ldetZa[22,0]
 ldetZa[22,1] ldetZa[22,2] ldetZa[22,3] ldetZa[22,4] ldetZa[22,5]
 ldetZa[22,6] ldetZa[22,7] ldetZa[22,8] ldetZa[22,9] ldetZa[22,10]
 ldetZa[22,11] ldetZa[22,12] ldetZa[22,13] ldetZa[22,14] ldetZa[22,15]
 ldetZa[23,0] ldetZa[23,1] ldetZa[23,2] ldetZa[23,3] ldetZa[23,4]
 ldetZa[23,5] ldetZa[23,6] ldetZa[23,7] ldetZa[23,8] ldetZa[23,9]
 ldetZa[23,10] ldetZa[23,11] ldetZa[23,12] ldetZa[23,13] ldetZa[23,14]
 ldetZa[23,15] ldetWa[19,0] ldetWa[19,1] ldetWa[19,2] ldetWa[19,3]
 ldetWa[19,4] ldetWa[19,5] ldetWa[19,6] ldetWa[19,7] ldetWa[19,8]
 ldetWa[19,9] ldetWa[19,10] ldetWa[19,11] ldetWa[19,12] ldetWa[19,13]
 ldetWa[19,14] ldetWa[19,15] ldetWa[20,0] ldetWa[20,1] ldetWa[20,2]
 ldetWa[20,3] ldetWa[20,4] ldetWa[20,5] ldetWa[20,6] ldetWa[20,7]
 ldetWa[20,8] ldetWa[20,9] ldetWa[20,10] ldetWa[20,11] ldetWa[20,12]
 ldetWa[20,13] ldetWa[20,14] ldetWa[20,15] ldetWa[21,0] ldetWa[21,1]
 ldetWa[21,2] ldetWa[21,3] ldetWa[21,4] ldetWa[21,5] ldetWa[21,6]
 ldetWa[21,7] ldetWa[21,8] ldetWa[21,9] ldetWa[21,10] ldetWa[21,11]
 ldetWa[21,12] ldetWa[21,13] ldetWa[21,14] ldetWa[21,15] ldetWa[22,0]
 ldetWa[22,1] ldetWa[22,2] ldetWa[22,3] ldetWa[22,4] ldetWa[22,5]
 ldetWa[22,6] ldetWa[22,7] ldetWa[22,8] ldetWa[22,9] ldetWa[22,10]
 ldetWa[22,11] ldetWa[22,12] ldetWa[22,13] ldetWa[22,14] ldetWa[22,15]
 lfrSCa[19,0] lfrSCa[19,1] lfrSCa[19,2] lfrSCa[19,3] lfrSCa[19,4]
 lfrSCa[19,5] lfrSCa[19,6] lfrSCa[19,7] lfrSCa[19,8] lfrSCa[19,9]
 lfrSCa[19,10] lfrSCa[19,11] lfrSCa[19,12] lfrSCa[19,13] lfrSCa[19,14]
 lfrSCa[19,15] lfrSCa[20,0] lfrSCa[20,1] lfrSCa[20,2] lfrSCa[20,3]
 lfrSCa[20,4] lfrSCa[20,5] lfrSCa[20,6] lfrSCa[20,7] lfrSCa[20,8]
 lfrSCa[20,9] lfrSCa[20,10] lfrSCa[20,11] lfrSCa[20,12] lfrSCa[20,13]
 lfrSCa[20,14] lfrSCa[20,15] lfrSCa[21,0] lfrSCa[21,1] lfrSCa[21,2]
 lfrSCa[21,3] lfrSCa[21,4] lfrSCa[21,5] lfrSCa[21,6] lfrSCa[21,7]
 lfrSCa[21,8] lfrSCa[21,9] lfrSCa[21,10] lfrSCa[21,11] lfrSCa[21,12]
 lfrSCa[21,13] lfrSCa[21,14] lfrSCa[21,15] lfrSCa[22,0] lfrSCa[22,1]
 lfrSCa[22,2] lfrSCa[22,3] lfrSCa[22,4] lfrSCa[22,5] lfrSCa[22,6]
 lfrSCa[22,7] lfrSCa[22,8] lfrSCa[22,9] lfrSCa[22,10] lfrSCa[22,11]
 lfrSCa[22,12] lfrSCa[22,13] lfrSCa[22,14] lfrSCa[22,15] lfrSCa[23,0]
 lfrSCa[23,1] lfrSCa[23,2] lfrSCa[23,3] lfrSCa[23,4] lfrSCa[23,5]
 lfrSCa[23,6] lfrSCa[23,7] lfrSCa[23,8] lfrSCa[23,9] lfrSCa[23,10]
 lfrSCa[23,11] lfrSCa[23,12] lfrSCa[23,13] lfrSCa[23,14] lfrSCa[23,15]
 lfrMCa[19,0] lfrMCa[19,1] lfrMCa[19,2] lfrMCa[19,3] lfrMCa[19,4]
 lfrMCa[19,5] lfrMCa[19,6] lfrMCa[19,7] lfrMCa[19,8] lfrMCa[19,9]
 lfrMCa[19,10] lfrMCa[19,11] lfrMCa[19,12] lfrMCa[19,13] lfrMCa[19,14]
 lfrMCa[19,15] lfrMCa[20,0] lfrMCa[20,1] lfrMCa[20,2] lfrMCa[20,3]
 lfrMCa[20,4] lfrMCa[20,5] lfrMCa[20,6] lfrMCa[20,7] lfrMCa[20,8]
 lfrMCa[20,9] lfrMCa[20,10] lfrMCa[20,11] lfrMCa[20,12] lfrMCa[20,13]
 lfrMCa[20,14] lfrMCa[20,15] lfrMCa[21,0] lfrMCa[21,1] lfrMCa[21,2]
 lfrMCa[21,3] lfrMCa[21,4] lfrMCa[21,5] lfrMCa[21,6] lfrMCa[21,7]
 lfrMCa[21,8] lfrMCa[21,9] lfrMCa[21,10] lfrMCa[21,11] lfrMCa[21,12]
 lfrMCa[21,13] lfrMCa[21,14] lfrMCa[21,15] lfrMCa[22,0] lfrMCa[22,1]
 lfrMCa[22,2] lfrMCa[22,3] lfrMCa[22,4] lfrMCa[22,5] lfrMCa[22,6]
 lfrMCa[22,7] lfrMCa[22,8] lfrMCa[22,9] lfrMCa[22,10] lfrMCa[22,11]
 lfrMCa[22,12] lfrMCa[22,13] lfrMCa[22,14] lfrMCa[22,15] ikz[0] ikz[1]
 ikz[2] ikz[3] ikz[4] ikz[5] ikz[6] ikz[7] ikz[8] ikz[9] ikz[10] ikz[11]
 ikz[12] ikz[13] ikz[14] ikz[15] gkz[0] gkz[1] gkz[2] gkz[3] gkz[4] gkz[5]
 gkz[6] gkz[7] gkz[8] gkz[9] gkz[10] gkz[11] gkz[12] gkz[13] gkz[14]
 gkz[15]
Generals
 iLANE[0] iLANE[1] iLANE[2] iLANE[3] iLANE[4] iLANE[5] iLANE[6] iLANE[7]
 iLANE[8] iLANE[9] iLANE[10] iLANE[11] iLANE[12] iLANE[13] iLANE[14]
 iLANE[15] ikSUM[0] ikSUM[1] ikSUM[2] ikSUM[3] ikSUM[4] ikSUM[5] ikSUM[6]
 ikSUM[7] ikSUM[8] ikSUM[9] ikSUM[10] ikSUM[11] ikSUM[12] ikSUM[13]
 ikSUM[14] ikSUM[15] gLANE[0] gLANE[1] gLANE[2] gLANE[3] gLANE[4] gLANE[5]
 gLANE[6] gLANE[7] gLANE[8] gLANE[9] gLANE[10] gLANE[11] gLANE[12]
 gLANE[13] gLANE[14] gLANE[15] gkSUM[0] gkSUM[1] gkSUM[2] gkSUM[3] gkSUM[4]
 gkSUM[5] gkSUM[6] gkSUM[7] gkSUM[8] gkSUM[9] gkSUM[10] gkSUM[11] gkSUM[12]
 gkSUM[13] gkSUM[14] gkSUM[15] rb rf cb cf iks gks cbp cfp
End
