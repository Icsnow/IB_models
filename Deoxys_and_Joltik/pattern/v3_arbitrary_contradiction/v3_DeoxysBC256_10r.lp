\ Model ./v3_DeoxysBC256_10r
\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize multi-objectives
  Tc: Priority=2 Weight=1 AbsTol=1e-06 RelTol=0
   Tc
  T32: Priority=1 Weight=1 AbsTol=1e-06 RelTol=0
   T32
Subject To
 uBas: uDW[-1,0] - udw[-1,0] >= 0
 uBas: uDW[0,0] - udw[0,0] >= 0
 uBas: uDW[1,0] - udw[1,0] >= 0
 uBas: uDW[2,0] - udw[2,0] >= 0
 uBas: uDW[3,0] - udw[3,0] >= 0
 uBas: uDW[4,0] - udw[4,0] >= 0
 uBas: uDW[5,0] - udw[5,0] >= 0
 uBas: uDW[6,0] - udw[6,0] >= 0
 uBas: uDW[7,0] - udw[7,0] >= 0
 uBas: uDX[0,0] - udx[0,0] >= 0
 uBas: uDY[0,0] - udy[0,0] >= 0
 uBas: uDZ[0,0] - udz[0,0] >= 0
 uBas: uDX[1,0] - udx[1,0] >= 0
 uBas: uDY[1,0] - udy[1,0] >= 0
 uBas: uDZ[1,0] - udz[1,0] >= 0
 uBas: uDX[2,0] - udx[2,0] >= 0
 uBas: uDY[2,0] - udy[2,0] >= 0
 uBas: uDZ[2,0] - udz[2,0] >= 0
 uBas: uDX[3,0] - udx[3,0] >= 0
 uBas: uDY[3,0] - udy[3,0] >= 0
 uBas: uDZ[3,0] - udz[3,0] >= 0
 uBas: uDX[4,0] - udx[4,0] >= 0
 uBas: uDY[4,0] - udy[4,0] >= 0
 uBas: uDZ[4,0] - udz[4,0] >= 0
 uBas: uDX[5,0] - udx[5,0] >= 0
 uBas: uDY[5,0] - udy[5,0] >= 0
 uBas: uDZ[5,0] - udz[5,0] >= 0
 uBas: uDX[6,0] - udx[6,0] >= 0
 uBas: uDY[6,0] - udy[6,0] >= 0
 uBas: uDZ[6,0] - udz[6,0] >= 0
 uBas: uDX[7,0] - udx[7,0] >= 0
 uBas: uDY[7,0] - udy[7,0] >= 0
 uBas: uDZ[7,0] - udz[7,0] >= 0
 uBas: uDW[-1,1] - udw[-1,1] >= 0
 uBas: uDW[0,1] - udw[0,1] >= 0
 uBas: uDW[1,1] - udw[1,1] >= 0
 uBas: uDW[2,1] - udw[2,1] >= 0
 uBas: uDW[3,1] - udw[3,1] >= 0
 uBas: uDW[4,1] - udw[4,1] >= 0
 uBas: uDW[5,1] - udw[5,1] >= 0
 uBas: uDW[6,1] - udw[6,1] >= 0
 uBas: uDW[7,1] - udw[7,1] >= 0
 uBas: uDX[0,1] - udx[0,1] >= 0
 uBas: uDY[0,1] - udy[0,1] >= 0
 uBas: uDZ[0,1] - udz[0,1] >= 0
 uBas: uDX[1,1] - udx[1,1] >= 0
 uBas: uDY[1,1] - udy[1,1] >= 0
 uBas: uDZ[1,1] - udz[1,1] >= 0
 uBas: uDX[2,1] - udx[2,1] >= 0
 uBas: uDY[2,1] - udy[2,1] >= 0
 uBas: uDZ[2,1] - udz[2,1] >= 0
 uBas: uDX[3,1] - udx[3,1] >= 0
 uBas: uDY[3,1] - udy[3,1] >= 0
 uBas: uDZ[3,1] - udz[3,1] >= 0
 uBas: uDX[4,1] - udx[4,1] >= 0
 uBas: uDY[4,1] - udy[4,1] >= 0
 uBas: uDZ[4,1] - udz[4,1] >= 0
 uBas: uDX[5,1] - udx[5,1] >= 0
 uBas: uDY[5,1] - udy[5,1] >= 0
 uBas: uDZ[5,1] - udz[5,1] >= 0
 uBas: uDX[6,1] - udx[6,1] >= 0
 uBas: uDY[6,1] - udy[6,1] >= 0
 uBas: uDZ[6,1] - udz[6,1] >= 0
 uBas: uDX[7,1] - udx[7,1] >= 0
 uBas: uDY[7,1] - udy[7,1] >= 0
 uBas: uDZ[7,1] - udz[7,1] >= 0
 uBas: uDW[-1,2] - udw[-1,2] >= 0
 uBas: uDW[0,2] - udw[0,2] >= 0
 uBas: uDW[1,2] - udw[1,2] >= 0
 uBas: uDW[2,2] - udw[2,2] >= 0
 uBas: uDW[3,2] - udw[3,2] >= 0
 uBas: uDW[4,2] - udw[4,2] >= 0
 uBas: uDW[5,2] - udw[5,2] >= 0
 uBas: uDW[6,2] - udw[6,2] >= 0
 uBas: uDW[7,2] - udw[7,2] >= 0
 uBas: uDX[0,2] - udx[0,2] >= 0
 uBas: uDY[0,2] - udy[0,2] >= 0
 uBas: uDZ[0,2] - udz[0,2] >= 0
 uBas: uDX[1,2] - udx[1,2] >= 0
 uBas: uDY[1,2] - udy[1,2] >= 0
 uBas: uDZ[1,2] - udz[1,2] >= 0
 uBas: uDX[2,2] - udx[2,2] >= 0
 uBas: uDY[2,2] - udy[2,2] >= 0
 uBas: uDZ[2,2] - udz[2,2] >= 0
 uBas: uDX[3,2] - udx[3,2] >= 0
 uBas: uDY[3,2] - udy[3,2] >= 0
 uBas: uDZ[3,2] - udz[3,2] >= 0
 uBas: uDX[4,2] - udx[4,2] >= 0
 uBas: uDY[4,2] - udy[4,2] >= 0
 uBas: uDZ[4,2] - udz[4,2] >= 0
 uBas: uDX[5,2] - udx[5,2] >= 0
 uBas: uDY[5,2] - udy[5,2] >= 0
 uBas: uDZ[5,2] - udz[5,2] >= 0
 uBas: uDX[6,2] - udx[6,2] >= 0
 uBas: uDY[6,2] - udy[6,2] >= 0
 uBas: uDZ[6,2] - udz[6,2] >= 0
 uBas: uDX[7,2] - udx[7,2] >= 0
 uBas: uDY[7,2] - udy[7,2] >= 0
 uBas: uDZ[7,2] - udz[7,2] >= 0
 uBas: uDW[-1,3] - udw[-1,3] >= 0
 uBas: uDW[0,3] - udw[0,3] >= 0
 uBas: uDW[1,3] - udw[1,3] >= 0
 uBas: uDW[2,3] - udw[2,3] >= 0
 uBas: uDW[3,3] - udw[3,3] >= 0
 uBas: uDW[4,3] - udw[4,3] >= 0
 uBas: uDW[5,3] - udw[5,3] >= 0
 uBas: uDW[6,3] - udw[6,3] >= 0
 uBas: uDW[7,3] - udw[7,3] >= 0
 uBas: uDX[0,3] - udx[0,3] >= 0
 uBas: uDY[0,3] - udy[0,3] >= 0
 uBas: uDZ[0,3] - udz[0,3] >= 0
 uBas: uDX[1,3] - udx[1,3] >= 0
 uBas: uDY[1,3] - udy[1,3] >= 0
 uBas: uDZ[1,3] - udz[1,3] >= 0
 uBas: uDX[2,3] - udx[2,3] >= 0
 uBas: uDY[2,3] - udy[2,3] >= 0
 uBas: uDZ[2,3] - udz[2,3] >= 0
 uBas: uDX[3,3] - udx[3,3] >= 0
 uBas: uDY[3,3] - udy[3,3] >= 0
 uBas: uDZ[3,3] - udz[3,3] >= 0
 uBas: uDX[4,3] - udx[4,3] >= 0
 uBas: uDY[4,3] - udy[4,3] >= 0
 uBas: uDZ[4,3] - udz[4,3] >= 0
 uBas: uDX[5,3] - udx[5,3] >= 0
 uBas: uDY[5,3] - udy[5,3] >= 0
 uBas: uDZ[5,3] - udz[5,3] >= 0
 uBas: uDX[6,3] - udx[6,3] >= 0
 uBas: uDY[6,3] - udy[6,3] >= 0
 uBas: uDZ[6,3] - udz[6,3] >= 0
 uBas: uDX[7,3] - udx[7,3] >= 0
 uBas: uDY[7,3] - udy[7,3] >= 0
 uBas: uDZ[7,3] - udz[7,3] >= 0
 uBas: uDW[-1,4] - udw[-1,4] >= 0
 uBas: uDW[0,4] - udw[0,4] >= 0
 uBas: uDW[1,4] - udw[1,4] >= 0
 uBas: uDW[2,4] - udw[2,4] >= 0
 uBas: uDW[3,4] - udw[3,4] >= 0
 uBas: uDW[4,4] - udw[4,4] >= 0
 uBas: uDW[5,4] - udw[5,4] >= 0
 uBas: uDW[6,4] - udw[6,4] >= 0
 uBas: uDW[7,4] - udw[7,4] >= 0
 uBas: uDX[0,4] - udx[0,4] >= 0
 uBas: uDY[0,4] - udy[0,4] >= 0
 uBas: uDZ[0,4] - udz[0,4] >= 0
 uBas: uDX[1,4] - udx[1,4] >= 0
 uBas: uDY[1,4] - udy[1,4] >= 0
 uBas: uDZ[1,4] - udz[1,4] >= 0
 uBas: uDX[2,4] - udx[2,4] >= 0
 uBas: uDY[2,4] - udy[2,4] >= 0
 uBas: uDZ[2,4] - udz[2,4] >= 0
 uBas: uDX[3,4] - udx[3,4] >= 0
 uBas: uDY[3,4] - udy[3,4] >= 0
 uBas: uDZ[3,4] - udz[3,4] >= 0
 uBas: uDX[4,4] - udx[4,4] >= 0
 uBas: uDY[4,4] - udy[4,4] >= 0
 uBas: uDZ[4,4] - udz[4,4] >= 0
 uBas: uDX[5,4] - udx[5,4] >= 0
 uBas: uDY[5,4] - udy[5,4] >= 0
 uBas: uDZ[5,4] - udz[5,4] >= 0
 uBas: uDX[6,4] - udx[6,4] >= 0
 uBas: uDY[6,4] - udy[6,4] >= 0
 uBas: uDZ[6,4] - udz[6,4] >= 0
 uBas: uDX[7,4] - udx[7,4] >= 0
 uBas: uDY[7,4] - udy[7,4] >= 0
 uBas: uDZ[7,4] - udz[7,4] >= 0
 uBas: uDW[-1,5] - udw[-1,5] >= 0
 uBas: uDW[0,5] - udw[0,5] >= 0
 uBas: uDW[1,5] - udw[1,5] >= 0
 uBas: uDW[2,5] - udw[2,5] >= 0
 uBas: uDW[3,5] - udw[3,5] >= 0
 uBas: uDW[4,5] - udw[4,5] >= 0
 uBas: uDW[5,5] - udw[5,5] >= 0
 uBas: uDW[6,5] - udw[6,5] >= 0
 uBas: uDW[7,5] - udw[7,5] >= 0
 uBas: uDX[0,5] - udx[0,5] >= 0
 uBas: uDY[0,5] - udy[0,5] >= 0
 uBas: uDZ[0,5] - udz[0,5] >= 0
 uBas: uDX[1,5] - udx[1,5] >= 0
 uBas: uDY[1,5] - udy[1,5] >= 0
 uBas: uDZ[1,5] - udz[1,5] >= 0
 uBas: uDX[2,5] - udx[2,5] >= 0
 uBas: uDY[2,5] - udy[2,5] >= 0
 uBas: uDZ[2,5] - udz[2,5] >= 0
 uBas: uDX[3,5] - udx[3,5] >= 0
 uBas: uDY[3,5] - udy[3,5] >= 0
 uBas: uDZ[3,5] - udz[3,5] >= 0
 uBas: uDX[4,5] - udx[4,5] >= 0
 uBas: uDY[4,5] - udy[4,5] >= 0
 uBas: uDZ[4,5] - udz[4,5] >= 0
 uBas: uDX[5,5] - udx[5,5] >= 0
 uBas: uDY[5,5] - udy[5,5] >= 0
 uBas: uDZ[5,5] - udz[5,5] >= 0
 uBas: uDX[6,5] - udx[6,5] >= 0
 uBas: uDY[6,5] - udy[6,5] >= 0
 uBas: uDZ[6,5] - udz[6,5] >= 0
 uBas: uDX[7,5] - udx[7,5] >= 0
 uBas: uDY[7,5] - udy[7,5] >= 0
 uBas: uDZ[7,5] - udz[7,5] >= 0
 uBas: uDW[-1,6] - udw[-1,6] >= 0
 uBas: uDW[0,6] - udw[0,6] >= 0
 uBas: uDW[1,6] - udw[1,6] >= 0
 uBas: uDW[2,6] - udw[2,6] >= 0
 uBas: uDW[3,6] - udw[3,6] >= 0
 uBas: uDW[4,6] - udw[4,6] >= 0
 uBas: uDW[5,6] - udw[5,6] >= 0
 uBas: uDW[6,6] - udw[6,6] >= 0
 uBas: uDW[7,6] - udw[7,6] >= 0
 uBas: uDX[0,6] - udx[0,6] >= 0
 uBas: uDY[0,6] - udy[0,6] >= 0
 uBas: uDZ[0,6] - udz[0,6] >= 0
 uBas: uDX[1,6] - udx[1,6] >= 0
 uBas: uDY[1,6] - udy[1,6] >= 0
 uBas: uDZ[1,6] - udz[1,6] >= 0
 uBas: uDX[2,6] - udx[2,6] >= 0
 uBas: uDY[2,6] - udy[2,6] >= 0
 uBas: uDZ[2,6] - udz[2,6] >= 0
 uBas: uDX[3,6] - udx[3,6] >= 0
 uBas: uDY[3,6] - udy[3,6] >= 0
 uBas: uDZ[3,6] - udz[3,6] >= 0
 uBas: uDX[4,6] - udx[4,6] >= 0
 uBas: uDY[4,6] - udy[4,6] >= 0
 uBas: uDZ[4,6] - udz[4,6] >= 0
 uBas: uDX[5,6] - udx[5,6] >= 0
 uBas: uDY[5,6] - udy[5,6] >= 0
 uBas: uDZ[5,6] - udz[5,6] >= 0
 uBas: uDX[6,6] - udx[6,6] >= 0
 uBas: uDY[6,6] - udy[6,6] >= 0
 uBas: uDZ[6,6] - udz[6,6] >= 0
 uBas: uDX[7,6] - udx[7,6] >= 0
 uBas: uDY[7,6] - udy[7,6] >= 0
 uBas: uDZ[7,6] - udz[7,6] >= 0
 uBas: uDW[-1,7] - udw[-1,7] >= 0
 uBas: uDW[0,7] - udw[0,7] >= 0
 uBas: uDW[1,7] - udw[1,7] >= 0
 uBas: uDW[2,7] - udw[2,7] >= 0
 uBas: uDW[3,7] - udw[3,7] >= 0
 uBas: uDW[4,7] - udw[4,7] >= 0
 uBas: uDW[5,7] - udw[5,7] >= 0
 uBas: uDW[6,7] - udw[6,7] >= 0
 uBas: uDW[7,7] - udw[7,7] >= 0
 uBas: uDX[0,7] - udx[0,7] >= 0
 uBas: uDY[0,7] - udy[0,7] >= 0
 uBas: uDZ[0,7] - udz[0,7] >= 0
 uBas: uDX[1,7] - udx[1,7] >= 0
 uBas: uDY[1,7] - udy[1,7] >= 0
 uBas: uDZ[1,7] - udz[1,7] >= 0
 uBas: uDX[2,7] - udx[2,7] >= 0
 uBas: uDY[2,7] - udy[2,7] >= 0
 uBas: uDZ[2,7] - udz[2,7] >= 0
 uBas: uDX[3,7] - udx[3,7] >= 0
 uBas: uDY[3,7] - udy[3,7] >= 0
 uBas: uDZ[3,7] - udz[3,7] >= 0
 uBas: uDX[4,7] - udx[4,7] >= 0
 uBas: uDY[4,7] - udy[4,7] >= 0
 uBas: uDZ[4,7] - udz[4,7] >= 0
 uBas: uDX[5,7] - udx[5,7] >= 0
 uBas: uDY[5,7] - udy[5,7] >= 0
 uBas: uDZ[5,7] - udz[5,7] >= 0
 uBas: uDX[6,7] - udx[6,7] >= 0
 uBas: uDY[6,7] - udy[6,7] >= 0
 uBas: uDZ[6,7] - udz[6,7] >= 0
 uBas: uDX[7,7] - udx[7,7] >= 0
 uBas: uDY[7,7] - udy[7,7] >= 0
 uBas: uDZ[7,7] - udz[7,7] >= 0
 uBas: uDW[-1,8] - udw[-1,8] >= 0
 uBas: uDW[0,8] - udw[0,8] >= 0
 uBas: uDW[1,8] - udw[1,8] >= 0
 uBas: uDW[2,8] - udw[2,8] >= 0
 uBas: uDW[3,8] - udw[3,8] >= 0
 uBas: uDW[4,8] - udw[4,8] >= 0
 uBas: uDW[5,8] - udw[5,8] >= 0
 uBas: uDW[6,8] - udw[6,8] >= 0
 uBas: uDW[7,8] - udw[7,8] >= 0
 uBas: uDX[0,8] - udx[0,8] >= 0
 uBas: uDY[0,8] - udy[0,8] >= 0
 uBas: uDZ[0,8] - udz[0,8] >= 0
 uBas: uDX[1,8] - udx[1,8] >= 0
 uBas: uDY[1,8] - udy[1,8] >= 0
 uBas: uDZ[1,8] - udz[1,8] >= 0
 uBas: uDX[2,8] - udx[2,8] >= 0
 uBas: uDY[2,8] - udy[2,8] >= 0
 uBas: uDZ[2,8] - udz[2,8] >= 0
 uBas: uDX[3,8] - udx[3,8] >= 0
 uBas: uDY[3,8] - udy[3,8] >= 0
 uBas: uDZ[3,8] - udz[3,8] >= 0
 uBas: uDX[4,8] - udx[4,8] >= 0
 uBas: uDY[4,8] - udy[4,8] >= 0
 uBas: uDZ[4,8] - udz[4,8] >= 0
 uBas: uDX[5,8] - udx[5,8] >= 0
 uBas: uDY[5,8] - udy[5,8] >= 0
 uBas: uDZ[5,8] - udz[5,8] >= 0
 uBas: uDX[6,8] - udx[6,8] >= 0
 uBas: uDY[6,8] - udy[6,8] >= 0
 uBas: uDZ[6,8] - udz[6,8] >= 0
 uBas: uDX[7,8] - udx[7,8] >= 0
 uBas: uDY[7,8] - udy[7,8] >= 0
 uBas: uDZ[7,8] - udz[7,8] >= 0
 uBas: uDW[-1,9] - udw[-1,9] >= 0
 uBas: uDW[0,9] - udw[0,9] >= 0
 uBas: uDW[1,9] - udw[1,9] >= 0
 uBas: uDW[2,9] - udw[2,9] >= 0
 uBas: uDW[3,9] - udw[3,9] >= 0
 uBas: uDW[4,9] - udw[4,9] >= 0
 uBas: uDW[5,9] - udw[5,9] >= 0
 uBas: uDW[6,9] - udw[6,9] >= 0
 uBas: uDW[7,9] - udw[7,9] >= 0
 uBas: uDX[0,9] - udx[0,9] >= 0
 uBas: uDY[0,9] - udy[0,9] >= 0
 uBas: uDZ[0,9] - udz[0,9] >= 0
 uBas: uDX[1,9] - udx[1,9] >= 0
 uBas: uDY[1,9] - udy[1,9] >= 0
 uBas: uDZ[1,9] - udz[1,9] >= 0
 uBas: uDX[2,9] - udx[2,9] >= 0
 uBas: uDY[2,9] - udy[2,9] >= 0
 uBas: uDZ[2,9] - udz[2,9] >= 0
 uBas: uDX[3,9] - udx[3,9] >= 0
 uBas: uDY[3,9] - udy[3,9] >= 0
 uBas: uDZ[3,9] - udz[3,9] >= 0
 uBas: uDX[4,9] - udx[4,9] >= 0
 uBas: uDY[4,9] - udy[4,9] >= 0
 uBas: uDZ[4,9] - udz[4,9] >= 0
 uBas: uDX[5,9] - udx[5,9] >= 0
 uBas: uDY[5,9] - udy[5,9] >= 0
 uBas: uDZ[5,9] - udz[5,9] >= 0
 uBas: uDX[6,9] - udx[6,9] >= 0
 uBas: uDY[6,9] - udy[6,9] >= 0
 uBas: uDZ[6,9] - udz[6,9] >= 0
 uBas: uDX[7,9] - udx[7,9] >= 0
 uBas: uDY[7,9] - udy[7,9] >= 0
 uBas: uDZ[7,9] - udz[7,9] >= 0
 uBas: uDW[-1,10] - udw[-1,10] >= 0
 uBas: uDW[0,10] - udw[0,10] >= 0
 uBas: uDW[1,10] - udw[1,10] >= 0
 uBas: uDW[2,10] - udw[2,10] >= 0
 uBas: uDW[3,10] - udw[3,10] >= 0
 uBas: uDW[4,10] - udw[4,10] >= 0
 uBas: uDW[5,10] - udw[5,10] >= 0
 uBas: uDW[6,10] - udw[6,10] >= 0
 uBas: uDW[7,10] - udw[7,10] >= 0
 uBas: uDX[0,10] - udx[0,10] >= 0
 uBas: uDY[0,10] - udy[0,10] >= 0
 uBas: uDZ[0,10] - udz[0,10] >= 0
 uBas: uDX[1,10] - udx[1,10] >= 0
 uBas: uDY[1,10] - udy[1,10] >= 0
 uBas: uDZ[1,10] - udz[1,10] >= 0
 uBas: uDX[2,10] - udx[2,10] >= 0
 uBas: uDY[2,10] - udy[2,10] >= 0
 uBas: uDZ[2,10] - udz[2,10] >= 0
 uBas: uDX[3,10] - udx[3,10] >= 0
 uBas: uDY[3,10] - udy[3,10] >= 0
 uBas: uDZ[3,10] - udz[3,10] >= 0
 uBas: uDX[4,10] - udx[4,10] >= 0
 uBas: uDY[4,10] - udy[4,10] >= 0
 uBas: uDZ[4,10] - udz[4,10] >= 0
 uBas: uDX[5,10] - udx[5,10] >= 0
 uBas: uDY[5,10] - udy[5,10] >= 0
 uBas: uDZ[5,10] - udz[5,10] >= 0
 uBas: uDX[6,10] - udx[6,10] >= 0
 uBas: uDY[6,10] - udy[6,10] >= 0
 uBas: uDZ[6,10] - udz[6,10] >= 0
 uBas: uDX[7,10] - udx[7,10] >= 0
 uBas: uDY[7,10] - udy[7,10] >= 0
 uBas: uDZ[7,10] - udz[7,10] >= 0
 uBas: uDW[-1,11] - udw[-1,11] >= 0
 uBas: uDW[0,11] - udw[0,11] >= 0
 uBas: uDW[1,11] - udw[1,11] >= 0
 uBas: uDW[2,11] - udw[2,11] >= 0
 uBas: uDW[3,11] - udw[3,11] >= 0
 uBas: uDW[4,11] - udw[4,11] >= 0
 uBas: uDW[5,11] - udw[5,11] >= 0
 uBas: uDW[6,11] - udw[6,11] >= 0
 uBas: uDW[7,11] - udw[7,11] >= 0
 uBas: uDX[0,11] - udx[0,11] >= 0
 uBas: uDY[0,11] - udy[0,11] >= 0
 uBas: uDZ[0,11] - udz[0,11] >= 0
 uBas: uDX[1,11] - udx[1,11] >= 0
 uBas: uDY[1,11] - udy[1,11] >= 0
 uBas: uDZ[1,11] - udz[1,11] >= 0
 uBas: uDX[2,11] - udx[2,11] >= 0
 uBas: uDY[2,11] - udy[2,11] >= 0
 uBas: uDZ[2,11] - udz[2,11] >= 0
 uBas: uDX[3,11] - udx[3,11] >= 0
 uBas: uDY[3,11] - udy[3,11] >= 0
 uBas: uDZ[3,11] - udz[3,11] >= 0
 uBas: uDX[4,11] - udx[4,11] >= 0
 uBas: uDY[4,11] - udy[4,11] >= 0
 uBas: uDZ[4,11] - udz[4,11] >= 0
 uBas: uDX[5,11] - udx[5,11] >= 0
 uBas: uDY[5,11] - udy[5,11] >= 0
 uBas: uDZ[5,11] - udz[5,11] >= 0
 uBas: uDX[6,11] - udx[6,11] >= 0
 uBas: uDY[6,11] - udy[6,11] >= 0
 uBas: uDZ[6,11] - udz[6,11] >= 0
 uBas: uDX[7,11] - udx[7,11] >= 0
 uBas: uDY[7,11] - udy[7,11] >= 0
 uBas: uDZ[7,11] - udz[7,11] >= 0
 uBas: uDW[-1,12] - udw[-1,12] >= 0
 uBas: uDW[0,12] - udw[0,12] >= 0
 uBas: uDW[1,12] - udw[1,12] >= 0
 uBas: uDW[2,12] - udw[2,12] >= 0
 uBas: uDW[3,12] - udw[3,12] >= 0
 uBas: uDW[4,12] - udw[4,12] >= 0
 uBas: uDW[5,12] - udw[5,12] >= 0
 uBas: uDW[6,12] - udw[6,12] >= 0
 uBas: uDW[7,12] - udw[7,12] >= 0
 uBas: uDX[0,12] - udx[0,12] >= 0
 uBas: uDY[0,12] - udy[0,12] >= 0
 uBas: uDZ[0,12] - udz[0,12] >= 0
 uBas: uDX[1,12] - udx[1,12] >= 0
 uBas: uDY[1,12] - udy[1,12] >= 0
 uBas: uDZ[1,12] - udz[1,12] >= 0
 uBas: uDX[2,12] - udx[2,12] >= 0
 uBas: uDY[2,12] - udy[2,12] >= 0
 uBas: uDZ[2,12] - udz[2,12] >= 0
 uBas: uDX[3,12] - udx[3,12] >= 0
 uBas: uDY[3,12] - udy[3,12] >= 0
 uBas: uDZ[3,12] - udz[3,12] >= 0
 uBas: uDX[4,12] - udx[4,12] >= 0
 uBas: uDY[4,12] - udy[4,12] >= 0
 uBas: uDZ[4,12] - udz[4,12] >= 0
 uBas: uDX[5,12] - udx[5,12] >= 0
 uBas: uDY[5,12] - udy[5,12] >= 0
 uBas: uDZ[5,12] - udz[5,12] >= 0
 uBas: uDX[6,12] - udx[6,12] >= 0
 uBas: uDY[6,12] - udy[6,12] >= 0
 uBas: uDZ[6,12] - udz[6,12] >= 0
 uBas: uDX[7,12] - udx[7,12] >= 0
 uBas: uDY[7,12] - udy[7,12] >= 0
 uBas: uDZ[7,12] - udz[7,12] >= 0
 uBas: uDW[-1,13] - udw[-1,13] >= 0
 uBas: uDW[0,13] - udw[0,13] >= 0
 uBas: uDW[1,13] - udw[1,13] >= 0
 uBas: uDW[2,13] - udw[2,13] >= 0
 uBas: uDW[3,13] - udw[3,13] >= 0
 uBas: uDW[4,13] - udw[4,13] >= 0
 uBas: uDW[5,13] - udw[5,13] >= 0
 uBas: uDW[6,13] - udw[6,13] >= 0
 uBas: uDW[7,13] - udw[7,13] >= 0
 uBas: uDX[0,13] - udx[0,13] >= 0
 uBas: uDY[0,13] - udy[0,13] >= 0
 uBas: uDZ[0,13] - udz[0,13] >= 0
 uBas: uDX[1,13] - udx[1,13] >= 0
 uBas: uDY[1,13] - udy[1,13] >= 0
 uBas: uDZ[1,13] - udz[1,13] >= 0
 uBas: uDX[2,13] - udx[2,13] >= 0
 uBas: uDY[2,13] - udy[2,13] >= 0
 uBas: uDZ[2,13] - udz[2,13] >= 0
 uBas: uDX[3,13] - udx[3,13] >= 0
 uBas: uDY[3,13] - udy[3,13] >= 0
 uBas: uDZ[3,13] - udz[3,13] >= 0
 uBas: uDX[4,13] - udx[4,13] >= 0
 uBas: uDY[4,13] - udy[4,13] >= 0
 uBas: uDZ[4,13] - udz[4,13] >= 0
 uBas: uDX[5,13] - udx[5,13] >= 0
 uBas: uDY[5,13] - udy[5,13] >= 0
 uBas: uDZ[5,13] - udz[5,13] >= 0
 uBas: uDX[6,13] - udx[6,13] >= 0
 uBas: uDY[6,13] - udy[6,13] >= 0
 uBas: uDZ[6,13] - udz[6,13] >= 0
 uBas: uDX[7,13] - udx[7,13] >= 0
 uBas: uDY[7,13] - udy[7,13] >= 0
 uBas: uDZ[7,13] - udz[7,13] >= 0
 uBas: uDW[-1,14] - udw[-1,14] >= 0
 uBas: uDW[0,14] - udw[0,14] >= 0
 uBas: uDW[1,14] - udw[1,14] >= 0
 uBas: uDW[2,14] - udw[2,14] >= 0
 uBas: uDW[3,14] - udw[3,14] >= 0
 uBas: uDW[4,14] - udw[4,14] >= 0
 uBas: uDW[5,14] - udw[5,14] >= 0
 uBas: uDW[6,14] - udw[6,14] >= 0
 uBas: uDW[7,14] - udw[7,14] >= 0
 uBas: uDX[0,14] - udx[0,14] >= 0
 uBas: uDY[0,14] - udy[0,14] >= 0
 uBas: uDZ[0,14] - udz[0,14] >= 0
 uBas: uDX[1,14] - udx[1,14] >= 0
 uBas: uDY[1,14] - udy[1,14] >= 0
 uBas: uDZ[1,14] - udz[1,14] >= 0
 uBas: uDX[2,14] - udx[2,14] >= 0
 uBas: uDY[2,14] - udy[2,14] >= 0
 uBas: uDZ[2,14] - udz[2,14] >= 0
 uBas: uDX[3,14] - udx[3,14] >= 0
 uBas: uDY[3,14] - udy[3,14] >= 0
 uBas: uDZ[3,14] - udz[3,14] >= 0
 uBas: uDX[4,14] - udx[4,14] >= 0
 uBas: uDY[4,14] - udy[4,14] >= 0
 uBas: uDZ[4,14] - udz[4,14] >= 0
 uBas: uDX[5,14] - udx[5,14] >= 0
 uBas: uDY[5,14] - udy[5,14] >= 0
 uBas: uDZ[5,14] - udz[5,14] >= 0
 uBas: uDX[6,14] - udx[6,14] >= 0
 uBas: uDY[6,14] - udy[6,14] >= 0
 uBas: uDZ[6,14] - udz[6,14] >= 0
 uBas: uDX[7,14] - udx[7,14] >= 0
 uBas: uDY[7,14] - udy[7,14] >= 0
 uBas: uDZ[7,14] - udz[7,14] >= 0
 uBas: uDW[-1,15] - udw[-1,15] >= 0
 uBas: uDW[0,15] - udw[0,15] >= 0
 uBas: uDW[1,15] - udw[1,15] >= 0
 uBas: uDW[2,15] - udw[2,15] >= 0
 uBas: uDW[3,15] - udw[3,15] >= 0
 uBas: uDW[4,15] - udw[4,15] >= 0
 uBas: uDW[5,15] - udw[5,15] >= 0
 uBas: uDW[6,15] - udw[6,15] >= 0
 uBas: uDW[7,15] - udw[7,15] >= 0
 uBas: uDX[0,15] - udx[0,15] >= 0
 uBas: uDY[0,15] - udy[0,15] >= 0
 uBas: uDZ[0,15] - udz[0,15] >= 0
 uBas: uDX[1,15] - udx[1,15] >= 0
 uBas: uDY[1,15] - udy[1,15] >= 0
 uBas: uDZ[1,15] - udz[1,15] >= 0
 uBas: uDX[2,15] - udx[2,15] >= 0
 uBas: uDY[2,15] - udy[2,15] >= 0
 uBas: uDZ[2,15] - udz[2,15] >= 0
 uBas: uDX[3,15] - udx[3,15] >= 0
 uBas: uDY[3,15] - udy[3,15] >= 0
 uBas: uDZ[3,15] - udz[3,15] >= 0
 uBas: uDX[4,15] - udx[4,15] >= 0
 uBas: uDY[4,15] - udy[4,15] >= 0
 uBas: uDZ[4,15] - udz[4,15] >= 0
 uBas: uDX[5,15] - udx[5,15] >= 0
 uBas: uDY[5,15] - udy[5,15] >= 0
 uBas: uDZ[5,15] - udz[5,15] >= 0
 uBas: uDX[6,15] - udx[6,15] >= 0
 uBas: uDY[6,15] - udy[6,15] >= 0
 uBas: uDZ[6,15] - udz[6,15] >= 0
 uBas: uDX[7,15] - udx[7,15] >= 0
 uBas: uDY[7,15] - udy[7,15] >= 0
 uBas: uDZ[7,15] - udz[7,15] >= 0
 uATK: uDW[-1,0] - udx[0,0] - ucan[0,0] >= 0
 uATK: - uDW[-1,0] + uDX[0,0] + ucan[0,0] >= 0
 uATK: uDW[-1,0] - uDX[0,0] + ustk[0,0] - 2 ucan[0,0] >= 0
 uATK: uDX[0,0] - ustk[0,0] + ucan[0,0] >= 0
 uATK: udw[-1,0] - udx[0,0] = 0
 uATK: uDW[-1,1] - udx[0,1] - ucan[0,1] >= 0
 uATK: - uDW[-1,1] + uDX[0,1] + ucan[0,1] >= 0
 uATK: uDW[-1,1] - uDX[0,1] + ustk[0,1] - 2 ucan[0,1] >= 0
 uATK: uDX[0,1] - ustk[0,1] + ucan[0,1] >= 0
 uATK: udw[-1,1] - udx[0,1] = 0
 uATK: uDW[-1,2] - udx[0,2] - ucan[0,2] >= 0
 uATK: - uDW[-1,2] + uDX[0,2] + ucan[0,2] >= 0
 uATK: uDW[-1,2] - uDX[0,2] + ustk[0,2] - 2 ucan[0,2] >= 0
 uATK: uDX[0,2] - ustk[0,2] + ucan[0,2] >= 0
 uATK: udw[-1,2] - udx[0,2] = 0
 uATK: uDW[-1,3] - udx[0,3] - ucan[0,3] >= 0
 uATK: - uDW[-1,3] + uDX[0,3] + ucan[0,3] >= 0
 uATK: uDW[-1,3] - uDX[0,3] + ustk[0,3] - 2 ucan[0,3] >= 0
 uATK: uDX[0,3] - ustk[0,3] + ucan[0,3] >= 0
 uATK: udw[-1,3] - udx[0,3] = 0
 uATK: uDW[-1,4] - udx[0,4] - ucan[0,4] >= 0
 uATK: - uDW[-1,4] + uDX[0,4] + ucan[0,4] >= 0
 uATK: uDW[-1,4] - uDX[0,4] + ustk[0,4] - 2 ucan[0,4] >= 0
 uATK: uDX[0,4] - ustk[0,4] + ucan[0,4] >= 0
 uATK: udw[-1,4] - udx[0,4] = 0
 uATK: uDW[-1,5] - udx[0,5] - ucan[0,5] >= 0
 uATK: - uDW[-1,5] + uDX[0,5] + ucan[0,5] >= 0
 uATK: uDW[-1,5] - uDX[0,5] + ustk[0,5] - 2 ucan[0,5] >= 0
 uATK: uDX[0,5] - ustk[0,5] + ucan[0,5] >= 0
 uATK: udw[-1,5] - udx[0,5] = 0
 uATK: uDW[-1,6] - udx[0,6] - ucan[0,6] >= 0
 uATK: - uDW[-1,6] + uDX[0,6] + ucan[0,6] >= 0
 uATK: uDW[-1,6] - uDX[0,6] + ustk[0,6] - 2 ucan[0,6] >= 0
 uATK: uDX[0,6] - ustk[0,6] + ucan[0,6] >= 0
 uATK: udw[-1,6] - udx[0,6] = 0
 uATK: uDW[-1,7] - udx[0,7] - ucan[0,7] >= 0
 uATK: - uDW[-1,7] + uDX[0,7] + ucan[0,7] >= 0
 uATK: uDW[-1,7] - uDX[0,7] + ustk[0,7] - 2 ucan[0,7] >= 0
 uATK: uDX[0,7] - ustk[0,7] + ucan[0,7] >= 0
 uATK: udw[-1,7] - udx[0,7] = 0
 uATK: uDW[-1,8] - udx[0,8] - ucan[0,8] >= 0
 uATK: - uDW[-1,8] + uDX[0,8] + ucan[0,8] >= 0
 uATK: uDW[-1,8] - uDX[0,8] + ustk[0,8] - 2 ucan[0,8] >= 0
 uATK: uDX[0,8] - ustk[0,8] + ucan[0,8] >= 0
 uATK: udw[-1,8] - udx[0,8] = 0
 uATK: uDW[-1,9] - udx[0,9] - ucan[0,9] >= 0
 uATK: - uDW[-1,9] + uDX[0,9] + ucan[0,9] >= 0
 uATK: uDW[-1,9] - uDX[0,9] + ustk[0,9] - 2 ucan[0,9] >= 0
 uATK: uDX[0,9] - ustk[0,9] + ucan[0,9] >= 0
 uATK: udw[-1,9] - udx[0,9] = 0
 uATK: uDW[-1,10] - udx[0,10] - ucan[0,10] >= 0
 uATK: - uDW[-1,10] + uDX[0,10] + ucan[0,10] >= 0
 uATK: uDW[-1,10] - uDX[0,10] + ustk[0,10] - 2 ucan[0,10] >= 0
 uATK: uDX[0,10] - ustk[0,10] + ucan[0,10] >= 0
 uATK: udw[-1,10] - udx[0,10] = 0
 uATK: uDW[-1,11] - udx[0,11] - ucan[0,11] >= 0
 uATK: - uDW[-1,11] + uDX[0,11] + ucan[0,11] >= 0
 uATK: uDW[-1,11] - uDX[0,11] + ustk[0,11] - 2 ucan[0,11] >= 0
 uATK: uDX[0,11] - ustk[0,11] + ucan[0,11] >= 0
 uATK: udw[-1,11] - udx[0,11] = 0
 uATK: uDW[-1,12] - udx[0,12] - ucan[0,12] >= 0
 uATK: - uDW[-1,12] + uDX[0,12] + ucan[0,12] >= 0
 uATK: uDW[-1,12] - uDX[0,12] + ustk[0,12] - 2 ucan[0,12] >= 0
 uATK: uDX[0,12] - ustk[0,12] + ucan[0,12] >= 0
 uATK: udw[-1,12] - udx[0,12] = 0
 uATK: uDW[-1,13] - udx[0,13] - ucan[0,13] >= 0
 uATK: - uDW[-1,13] + uDX[0,13] + ucan[0,13] >= 0
 uATK: uDW[-1,13] - uDX[0,13] + ustk[0,13] - 2 ucan[0,13] >= 0
 uATK: uDX[0,13] - ustk[0,13] + ucan[0,13] >= 0
 uATK: udw[-1,13] - udx[0,13] = 0
 uATK: uDW[-1,14] - udx[0,14] - ucan[0,14] >= 0
 uATK: - uDW[-1,14] + uDX[0,14] + ucan[0,14] >= 0
 uATK: uDW[-1,14] - uDX[0,14] + ustk[0,14] - 2 ucan[0,14] >= 0
 uATK: uDX[0,14] - ustk[0,14] + ucan[0,14] >= 0
 uATK: udw[-1,14] - udx[0,14] = 0
 uATK: uDW[-1,15] - udx[0,15] - ucan[0,15] >= 0
 uATK: - uDW[-1,15] + uDX[0,15] + ucan[0,15] >= 0
 uATK: uDW[-1,15] - uDX[0,15] + ustk[0,15] - 2 ucan[0,15] >= 0
 uATK: uDX[0,15] - ustk[0,15] + ucan[0,15] >= 0
 uATK: udw[-1,15] - udx[0,15] = 0
 uATK: uDW[0,0] - udx[1,0] - ucan[1,0] >= 0
 uATK: - uDW[0,0] + uDX[1,0] + ucan[1,0] >= 0
 uATK: uDW[0,0] - uDX[1,0] + ustk[1,0] - 2 ucan[1,0] >= 0
 uATK: uDX[1,0] - ustk[1,0] + ucan[1,0] >= 0
 uATK: udw[0,0] - udx[1,0] = 0
 uATK: uDW[0,1] - udx[1,1] - ucan[1,1] >= 0
 uATK: - uDW[0,1] + uDX[1,1] + ucan[1,1] >= 0
 uATK: uDW[0,1] - uDX[1,1] + ustk[1,1] - 2 ucan[1,1] >= 0
 uATK: uDX[1,1] - ustk[1,1] + ucan[1,1] >= 0
 uATK: udw[0,1] - udx[1,1] = 0
 uATK: uDW[0,2] - udx[1,2] - ucan[1,2] >= 0
 uATK: - uDW[0,2] + uDX[1,2] + ucan[1,2] >= 0
 uATK: uDW[0,2] - uDX[1,2] + ustk[1,2] - 2 ucan[1,2] >= 0
 uATK: uDX[1,2] - ustk[1,2] + ucan[1,2] >= 0
 uATK: udw[0,2] - udx[1,2] = 0
 uATK: uDW[0,3] - udx[1,3] - ucan[1,3] >= 0
 uATK: - uDW[0,3] + uDX[1,3] + ucan[1,3] >= 0
 uATK: uDW[0,3] - uDX[1,3] + ustk[1,3] - 2 ucan[1,3] >= 0
 uATK: uDX[1,3] - ustk[1,3] + ucan[1,3] >= 0
 uATK: udw[0,3] - udx[1,3] = 0
 uATK: uDW[0,4] - udx[1,4] - ucan[1,4] >= 0
 uATK: - uDW[0,4] + uDX[1,4] + ucan[1,4] >= 0
 uATK: uDW[0,4] - uDX[1,4] + ustk[1,4] - 2 ucan[1,4] >= 0
 uATK: uDX[1,4] - ustk[1,4] + ucan[1,4] >= 0
 uATK: udw[0,4] - udx[1,4] = 0
 uATK: uDW[0,5] - udx[1,5] - ucan[1,5] >= 0
 uATK: - uDW[0,5] + uDX[1,5] + ucan[1,5] >= 0
 uATK: uDW[0,5] - uDX[1,5] + ustk[1,5] - 2 ucan[1,5] >= 0
 uATK: uDX[1,5] - ustk[1,5] + ucan[1,5] >= 0
 uATK: udw[0,5] - udx[1,5] = 0
 uATK: uDW[0,6] - udx[1,6] - ucan[1,6] >= 0
 uATK: - uDW[0,6] + uDX[1,6] + ucan[1,6] >= 0
 uATK: uDW[0,6] - uDX[1,6] + ustk[1,6] - 2 ucan[1,6] >= 0
 uATK: uDX[1,6] - ustk[1,6] + ucan[1,6] >= 0
 uATK: udw[0,6] - udx[1,6] = 0
 uATK: uDW[0,7] - udx[1,7] - ucan[1,7] >= 0
 uATK: - uDW[0,7] + uDX[1,7] + ucan[1,7] >= 0
 uATK: uDW[0,7] - uDX[1,7] + ustk[1,7] - 2 ucan[1,7] >= 0
 uATK: uDX[1,7] - ustk[1,7] + ucan[1,7] >= 0
 uATK: udw[0,7] - udx[1,7] = 0
 uATK: uDW[0,8] - udx[1,8] - ucan[1,8] >= 0
 uATK: - uDW[0,8] + uDX[1,8] + ucan[1,8] >= 0
 uATK: uDW[0,8] - uDX[1,8] + ustk[1,8] - 2 ucan[1,8] >= 0
 uATK: uDX[1,8] - ustk[1,8] + ucan[1,8] >= 0
 uATK: udw[0,8] - udx[1,8] = 0
 uATK: uDW[0,9] - udx[1,9] - ucan[1,9] >= 0
 uATK: - uDW[0,9] + uDX[1,9] + ucan[1,9] >= 0
 uATK: uDW[0,9] - uDX[1,9] + ustk[1,9] - 2 ucan[1,9] >= 0
 uATK: uDX[1,9] - ustk[1,9] + ucan[1,9] >= 0
 uATK: udw[0,9] - udx[1,9] = 0
 uATK: uDW[0,10] - udx[1,10] - ucan[1,10] >= 0
 uATK: - uDW[0,10] + uDX[1,10] + ucan[1,10] >= 0
 uATK: uDW[0,10] - uDX[1,10] + ustk[1,10] - 2 ucan[1,10] >= 0
 uATK: uDX[1,10] - ustk[1,10] + ucan[1,10] >= 0
 uATK: udw[0,10] - udx[1,10] = 0
 uATK: uDW[0,11] - udx[1,11] - ucan[1,11] >= 0
 uATK: - uDW[0,11] + uDX[1,11] + ucan[1,11] >= 0
 uATK: uDW[0,11] - uDX[1,11] + ustk[1,11] - 2 ucan[1,11] >= 0
 uATK: uDX[1,11] - ustk[1,11] + ucan[1,11] >= 0
 uATK: udw[0,11] - udx[1,11] = 0
 uATK: uDW[0,12] - udx[1,12] - ucan[1,12] >= 0
 uATK: - uDW[0,12] + uDX[1,12] + ucan[1,12] >= 0
 uATK: uDW[0,12] - uDX[1,12] + ustk[1,12] - 2 ucan[1,12] >= 0
 uATK: uDX[1,12] - ustk[1,12] + ucan[1,12] >= 0
 uATK: udw[0,12] - udx[1,12] = 0
 uATK: uDW[0,13] - udx[1,13] - ucan[1,13] >= 0
 uATK: - uDW[0,13] + uDX[1,13] + ucan[1,13] >= 0
 uATK: uDW[0,13] - uDX[1,13] + ustk[1,13] - 2 ucan[1,13] >= 0
 uATK: uDX[1,13] - ustk[1,13] + ucan[1,13] >= 0
 uATK: udw[0,13] - udx[1,13] = 0
 uATK: uDW[0,14] - udx[1,14] - ucan[1,14] >= 0
 uATK: - uDW[0,14] + uDX[1,14] + ucan[1,14] >= 0
 uATK: uDW[0,14] - uDX[1,14] + ustk[1,14] - 2 ucan[1,14] >= 0
 uATK: uDX[1,14] - ustk[1,14] + ucan[1,14] >= 0
 uATK: udw[0,14] - udx[1,14] = 0
 uATK: uDW[0,15] - udx[1,15] - ucan[1,15] >= 0
 uATK: - uDW[0,15] + uDX[1,15] + ucan[1,15] >= 0
 uATK: uDW[0,15] - uDX[1,15] + ustk[1,15] - 2 ucan[1,15] >= 0
 uATK: uDX[1,15] - ustk[1,15] + ucan[1,15] >= 0
 uATK: udw[0,15] - udx[1,15] = 0
 uATK: uDW[1,0] - udx[2,0] - ucan[2,0] >= 0
 uATK: - uDW[1,0] + uDX[2,0] + ucan[2,0] >= 0
 uATK: uDW[1,0] - uDX[2,0] + ustk[2,0] - 2 ucan[2,0] >= 0
 uATK: uDX[2,0] - ustk[2,0] + ucan[2,0] >= 0
 uATK: udw[1,0] - udx[2,0] = 0
 uATK: uDW[1,1] - udx[2,1] - ucan[2,1] >= 0
 uATK: - uDW[1,1] + uDX[2,1] + ucan[2,1] >= 0
 uATK: uDW[1,1] - uDX[2,1] + ustk[2,1] - 2 ucan[2,1] >= 0
 uATK: uDX[2,1] - ustk[2,1] + ucan[2,1] >= 0
 uATK: udw[1,1] - udx[2,1] = 0
 uATK: uDW[1,2] - udx[2,2] - ucan[2,2] >= 0
 uATK: - uDW[1,2] + uDX[2,2] + ucan[2,2] >= 0
 uATK: uDW[1,2] - uDX[2,2] + ustk[2,2] - 2 ucan[2,2] >= 0
 uATK: uDX[2,2] - ustk[2,2] + ucan[2,2] >= 0
 uATK: udw[1,2] - udx[2,2] = 0
 uATK: uDW[1,3] - udx[2,3] - ucan[2,3] >= 0
 uATK: - uDW[1,3] + uDX[2,3] + ucan[2,3] >= 0
 uATK: uDW[1,3] - uDX[2,3] + ustk[2,3] - 2 ucan[2,3] >= 0
 uATK: uDX[2,3] - ustk[2,3] + ucan[2,3] >= 0
 uATK: udw[1,3] - udx[2,3] = 0
 uATK: uDW[1,4] - udx[2,4] - ucan[2,4] >= 0
 uATK: - uDW[1,4] + uDX[2,4] + ucan[2,4] >= 0
 uATK: uDW[1,4] - uDX[2,4] + ustk[2,4] - 2 ucan[2,4] >= 0
 uATK: uDX[2,4] - ustk[2,4] + ucan[2,4] >= 0
 uATK: udw[1,4] - udx[2,4] = 0
 uATK: uDW[1,5] - udx[2,5] - ucan[2,5] >= 0
 uATK: - uDW[1,5] + uDX[2,5] + ucan[2,5] >= 0
 uATK: uDW[1,5] - uDX[2,5] + ustk[2,5] - 2 ucan[2,5] >= 0
 uATK: uDX[2,5] - ustk[2,5] + ucan[2,5] >= 0
 uATK: udw[1,5] - udx[2,5] = 0
 uATK: uDW[1,6] - udx[2,6] - ucan[2,6] >= 0
 uATK: - uDW[1,6] + uDX[2,6] + ucan[2,6] >= 0
 uATK: uDW[1,6] - uDX[2,6] + ustk[2,6] - 2 ucan[2,6] >= 0
 uATK: uDX[2,6] - ustk[2,6] + ucan[2,6] >= 0
 uATK: udw[1,6] - udx[2,6] = 0
 uATK: uDW[1,7] - udx[2,7] - ucan[2,7] >= 0
 uATK: - uDW[1,7] + uDX[2,7] + ucan[2,7] >= 0
 uATK: uDW[1,7] - uDX[2,7] + ustk[2,7] - 2 ucan[2,7] >= 0
 uATK: uDX[2,7] - ustk[2,7] + ucan[2,7] >= 0
 uATK: udw[1,7] - udx[2,7] = 0
 uATK: uDW[1,8] - udx[2,8] - ucan[2,8] >= 0
 uATK: - uDW[1,8] + uDX[2,8] + ucan[2,8] >= 0
 uATK: uDW[1,8] - uDX[2,8] + ustk[2,8] - 2 ucan[2,8] >= 0
 uATK: uDX[2,8] - ustk[2,8] + ucan[2,8] >= 0
 uATK: udw[1,8] - udx[2,8] = 0
 uATK: uDW[1,9] - udx[2,9] - ucan[2,9] >= 0
 uATK: - uDW[1,9] + uDX[2,9] + ucan[2,9] >= 0
 uATK: uDW[1,9] - uDX[2,9] + ustk[2,9] - 2 ucan[2,9] >= 0
 uATK: uDX[2,9] - ustk[2,9] + ucan[2,9] >= 0
 uATK: udw[1,9] - udx[2,9] = 0
 uATK: uDW[1,10] - udx[2,10] - ucan[2,10] >= 0
 uATK: - uDW[1,10] + uDX[2,10] + ucan[2,10] >= 0
 uATK: uDW[1,10] - uDX[2,10] + ustk[2,10] - 2 ucan[2,10] >= 0
 uATK: uDX[2,10] - ustk[2,10] + ucan[2,10] >= 0
 uATK: udw[1,10] - udx[2,10] = 0
 uATK: uDW[1,11] - udx[2,11] - ucan[2,11] >= 0
 uATK: - uDW[1,11] + uDX[2,11] + ucan[2,11] >= 0
 uATK: uDW[1,11] - uDX[2,11] + ustk[2,11] - 2 ucan[2,11] >= 0
 uATK: uDX[2,11] - ustk[2,11] + ucan[2,11] >= 0
 uATK: udw[1,11] - udx[2,11] = 0
 uATK: uDW[1,12] - udx[2,12] - ucan[2,12] >= 0
 uATK: - uDW[1,12] + uDX[2,12] + ucan[2,12] >= 0
 uATK: uDW[1,12] - uDX[2,12] + ustk[2,12] - 2 ucan[2,12] >= 0
 uATK: uDX[2,12] - ustk[2,12] + ucan[2,12] >= 0
 uATK: udw[1,12] - udx[2,12] = 0
 uATK: uDW[1,13] - udx[2,13] - ucan[2,13] >= 0
 uATK: - uDW[1,13] + uDX[2,13] + ucan[2,13] >= 0
 uATK: uDW[1,13] - uDX[2,13] + ustk[2,13] - 2 ucan[2,13] >= 0
 uATK: uDX[2,13] - ustk[2,13] + ucan[2,13] >= 0
 uATK: udw[1,13] - udx[2,13] = 0
 uATK: uDW[1,14] - udx[2,14] - ucan[2,14] >= 0
 uATK: - uDW[1,14] + uDX[2,14] + ucan[2,14] >= 0
 uATK: uDW[1,14] - uDX[2,14] + ustk[2,14] - 2 ucan[2,14] >= 0
 uATK: uDX[2,14] - ustk[2,14] + ucan[2,14] >= 0
 uATK: udw[1,14] - udx[2,14] = 0
 uATK: uDW[1,15] - udx[2,15] - ucan[2,15] >= 0
 uATK: - uDW[1,15] + uDX[2,15] + ucan[2,15] >= 0
 uATK: uDW[1,15] - uDX[2,15] + ustk[2,15] - 2 ucan[2,15] >= 0
 uATK: uDX[2,15] - ustk[2,15] + ucan[2,15] >= 0
 uATK: udw[1,15] - udx[2,15] = 0
 uATK: uDW[2,0] - udx[3,0] - ucan[3,0] >= 0
 uATK: - uDW[2,0] + uDX[3,0] + ucan[3,0] >= 0
 uATK: uDW[2,0] - uDX[3,0] + ustk[3,0] - 2 ucan[3,0] >= 0
 uATK: uDX[3,0] - ustk[3,0] + ucan[3,0] >= 0
 uATK: udw[2,0] - udx[3,0] = 0
 uATK: uDW[2,1] - udx[3,1] - ucan[3,1] >= 0
 uATK: - uDW[2,1] + uDX[3,1] + ucan[3,1] >= 0
 uATK: uDW[2,1] - uDX[3,1] + ustk[3,1] - 2 ucan[3,1] >= 0
 uATK: uDX[3,1] - ustk[3,1] + ucan[3,1] >= 0
 uATK: udw[2,1] - udx[3,1] = 0
 uATK: uDW[2,2] - udx[3,2] - ucan[3,2] >= 0
 uATK: - uDW[2,2] + uDX[3,2] + ucan[3,2] >= 0
 uATK: uDW[2,2] - uDX[3,2] + ustk[3,2] - 2 ucan[3,2] >= 0
 uATK: uDX[3,2] - ustk[3,2] + ucan[3,2] >= 0
 uATK: udw[2,2] - udx[3,2] = 0
 uATK: uDW[2,3] - udx[3,3] - ucan[3,3] >= 0
 uATK: - uDW[2,3] + uDX[3,3] + ucan[3,3] >= 0
 uATK: uDW[2,3] - uDX[3,3] + ustk[3,3] - 2 ucan[3,3] >= 0
 uATK: uDX[3,3] - ustk[3,3] + ucan[3,3] >= 0
 uATK: udw[2,3] - udx[3,3] = 0
 uATK: uDW[2,4] - udx[3,4] - ucan[3,4] >= 0
 uATK: - uDW[2,4] + uDX[3,4] + ucan[3,4] >= 0
 uATK: uDW[2,4] - uDX[3,4] + ustk[3,4] - 2 ucan[3,4] >= 0
 uATK: uDX[3,4] - ustk[3,4] + ucan[3,4] >= 0
 uATK: udw[2,4] - udx[3,4] = 0
 uATK: uDW[2,5] - udx[3,5] - ucan[3,5] >= 0
 uATK: - uDW[2,5] + uDX[3,5] + ucan[3,5] >= 0
 uATK: uDW[2,5] - uDX[3,5] + ustk[3,5] - 2 ucan[3,5] >= 0
 uATK: uDX[3,5] - ustk[3,5] + ucan[3,5] >= 0
 uATK: udw[2,5] - udx[3,5] = 0
 uATK: uDW[2,6] - udx[3,6] - ucan[3,6] >= 0
 uATK: - uDW[2,6] + uDX[3,6] + ucan[3,6] >= 0
 uATK: uDW[2,6] - uDX[3,6] + ustk[3,6] - 2 ucan[3,6] >= 0
 uATK: uDX[3,6] - ustk[3,6] + ucan[3,6] >= 0
 uATK: udw[2,6] - udx[3,6] = 0
 uATK: uDW[2,7] - udx[3,7] - ucan[3,7] >= 0
 uATK: - uDW[2,7] + uDX[3,7] + ucan[3,7] >= 0
 uATK: uDW[2,7] - uDX[3,7] + ustk[3,7] - 2 ucan[3,7] >= 0
 uATK: uDX[3,7] - ustk[3,7] + ucan[3,7] >= 0
 uATK: udw[2,7] - udx[3,7] = 0
 uATK: uDW[2,8] - udx[3,8] - ucan[3,8] >= 0
 uATK: - uDW[2,8] + uDX[3,8] + ucan[3,8] >= 0
 uATK: uDW[2,8] - uDX[3,8] + ustk[3,8] - 2 ucan[3,8] >= 0
 uATK: uDX[3,8] - ustk[3,8] + ucan[3,8] >= 0
 uATK: udw[2,8] - udx[3,8] = 0
 uATK: uDW[2,9] - udx[3,9] - ucan[3,9] >= 0
 uATK: - uDW[2,9] + uDX[3,9] + ucan[3,9] >= 0
 uATK: uDW[2,9] - uDX[3,9] + ustk[3,9] - 2 ucan[3,9] >= 0
 uATK: uDX[3,9] - ustk[3,9] + ucan[3,9] >= 0
 uATK: udw[2,9] - udx[3,9] = 0
 uATK: uDW[2,10] - udx[3,10] - ucan[3,10] >= 0
 uATK: - uDW[2,10] + uDX[3,10] + ucan[3,10] >= 0
 uATK: uDW[2,10] - uDX[3,10] + ustk[3,10] - 2 ucan[3,10] >= 0
 uATK: uDX[3,10] - ustk[3,10] + ucan[3,10] >= 0
 uATK: udw[2,10] - udx[3,10] = 0
 uATK: uDW[2,11] - udx[3,11] - ucan[3,11] >= 0
 uATK: - uDW[2,11] + uDX[3,11] + ucan[3,11] >= 0
 uATK: uDW[2,11] - uDX[3,11] + ustk[3,11] - 2 ucan[3,11] >= 0
 uATK: uDX[3,11] - ustk[3,11] + ucan[3,11] >= 0
 uATK: udw[2,11] - udx[3,11] = 0
 uATK: uDW[2,12] - udx[3,12] - ucan[3,12] >= 0
 uATK: - uDW[2,12] + uDX[3,12] + ucan[3,12] >= 0
 uATK: uDW[2,12] - uDX[3,12] + ustk[3,12] - 2 ucan[3,12] >= 0
 uATK: uDX[3,12] - ustk[3,12] + ucan[3,12] >= 0
 uATK: udw[2,12] - udx[3,12] = 0
 uATK: uDW[2,13] - udx[3,13] - ucan[3,13] >= 0
 uATK: - uDW[2,13] + uDX[3,13] + ucan[3,13] >= 0
 uATK: uDW[2,13] - uDX[3,13] + ustk[3,13] - 2 ucan[3,13] >= 0
 uATK: uDX[3,13] - ustk[3,13] + ucan[3,13] >= 0
 uATK: udw[2,13] - udx[3,13] = 0
 uATK: uDW[2,14] - udx[3,14] - ucan[3,14] >= 0
 uATK: - uDW[2,14] + uDX[3,14] + ucan[3,14] >= 0
 uATK: uDW[2,14] - uDX[3,14] + ustk[3,14] - 2 ucan[3,14] >= 0
 uATK: uDX[3,14] - ustk[3,14] + ucan[3,14] >= 0
 uATK: udw[2,14] - udx[3,14] = 0
 uATK: uDW[2,15] - udx[3,15] - ucan[3,15] >= 0
 uATK: - uDW[2,15] + uDX[3,15] + ucan[3,15] >= 0
 uATK: uDW[2,15] - uDX[3,15] + ustk[3,15] - 2 ucan[3,15] >= 0
 uATK: uDX[3,15] - ustk[3,15] + ucan[3,15] >= 0
 uATK: udw[2,15] - udx[3,15] = 0
 uATK: uDW[3,0] - udx[4,0] - ucan[4,0] >= 0
 uATK: - uDW[3,0] + uDX[4,0] + ucan[4,0] >= 0
 uATK: uDW[3,0] - uDX[4,0] + ustk[4,0] - 2 ucan[4,0] >= 0
 uATK: uDX[4,0] - ustk[4,0] + ucan[4,0] >= 0
 uATK: udw[3,0] - udx[4,0] = 0
 uATK: uDW[3,1] - udx[4,1] - ucan[4,1] >= 0
 uATK: - uDW[3,1] + uDX[4,1] + ucan[4,1] >= 0
 uATK: uDW[3,1] - uDX[4,1] + ustk[4,1] - 2 ucan[4,1] >= 0
 uATK: uDX[4,1] - ustk[4,1] + ucan[4,1] >= 0
 uATK: udw[3,1] - udx[4,1] = 0
 uATK: uDW[3,2] - udx[4,2] - ucan[4,2] >= 0
 uATK: - uDW[3,2] + uDX[4,2] + ucan[4,2] >= 0
 uATK: uDW[3,2] - uDX[4,2] + ustk[4,2] - 2 ucan[4,2] >= 0
 uATK: uDX[4,2] - ustk[4,2] + ucan[4,2] >= 0
 uATK: udw[3,2] - udx[4,2] = 0
 uATK: uDW[3,3] - udx[4,3] - ucan[4,3] >= 0
 uATK: - uDW[3,3] + uDX[4,3] + ucan[4,3] >= 0
 uATK: uDW[3,3] - uDX[4,3] + ustk[4,3] - 2 ucan[4,3] >= 0
 uATK: uDX[4,3] - ustk[4,3] + ucan[4,3] >= 0
 uATK: udw[3,3] - udx[4,3] = 0
 uATK: uDW[3,4] - udx[4,4] - ucan[4,4] >= 0
 uATK: - uDW[3,4] + uDX[4,4] + ucan[4,4] >= 0
 uATK: uDW[3,4] - uDX[4,4] + ustk[4,4] - 2 ucan[4,4] >= 0
 uATK: uDX[4,4] - ustk[4,4] + ucan[4,4] >= 0
 uATK: udw[3,4] - udx[4,4] = 0
 uATK: uDW[3,5] - udx[4,5] - ucan[4,5] >= 0
 uATK: - uDW[3,5] + uDX[4,5] + ucan[4,5] >= 0
 uATK: uDW[3,5] - uDX[4,5] + ustk[4,5] - 2 ucan[4,5] >= 0
 uATK: uDX[4,5] - ustk[4,5] + ucan[4,5] >= 0
 uATK: udw[3,5] - udx[4,5] = 0
 uATK: uDW[3,6] - udx[4,6] - ucan[4,6] >= 0
 uATK: - uDW[3,6] + uDX[4,6] + ucan[4,6] >= 0
 uATK: uDW[3,6] - uDX[4,6] + ustk[4,6] - 2 ucan[4,6] >= 0
 uATK: uDX[4,6] - ustk[4,6] + ucan[4,6] >= 0
 uATK: udw[3,6] - udx[4,6] = 0
 uATK: uDW[3,7] - udx[4,7] - ucan[4,7] >= 0
 uATK: - uDW[3,7] + uDX[4,7] + ucan[4,7] >= 0
 uATK: uDW[3,7] - uDX[4,7] + ustk[4,7] - 2 ucan[4,7] >= 0
 uATK: uDX[4,7] - ustk[4,7] + ucan[4,7] >= 0
 uATK: udw[3,7] - udx[4,7] = 0
 uATK: uDW[3,8] - udx[4,8] - ucan[4,8] >= 0
 uATK: - uDW[3,8] + uDX[4,8] + ucan[4,8] >= 0
 uATK: uDW[3,8] - uDX[4,8] + ustk[4,8] - 2 ucan[4,8] >= 0
 uATK: uDX[4,8] - ustk[4,8] + ucan[4,8] >= 0
 uATK: udw[3,8] - udx[4,8] = 0
 uATK: uDW[3,9] - udx[4,9] - ucan[4,9] >= 0
 uATK: - uDW[3,9] + uDX[4,9] + ucan[4,9] >= 0
 uATK: uDW[3,9] - uDX[4,9] + ustk[4,9] - 2 ucan[4,9] >= 0
 uATK: uDX[4,9] - ustk[4,9] + ucan[4,9] >= 0
 uATK: udw[3,9] - udx[4,9] = 0
 uATK: uDW[3,10] - udx[4,10] - ucan[4,10] >= 0
 uATK: - uDW[3,10] + uDX[4,10] + ucan[4,10] >= 0
 uATK: uDW[3,10] - uDX[4,10] + ustk[4,10] - 2 ucan[4,10] >= 0
 uATK: uDX[4,10] - ustk[4,10] + ucan[4,10] >= 0
 uATK: udw[3,10] - udx[4,10] = 0
 uATK: uDW[3,11] - udx[4,11] - ucan[4,11] >= 0
 uATK: - uDW[3,11] + uDX[4,11] + ucan[4,11] >= 0
 uATK: uDW[3,11] - uDX[4,11] + ustk[4,11] - 2 ucan[4,11] >= 0
 uATK: uDX[4,11] - ustk[4,11] + ucan[4,11] >= 0
 uATK: udw[3,11] - udx[4,11] = 0
 uATK: uDW[3,12] - udx[4,12] - ucan[4,12] >= 0
 uATK: - uDW[3,12] + uDX[4,12] + ucan[4,12] >= 0
 uATK: uDW[3,12] - uDX[4,12] + ustk[4,12] - 2 ucan[4,12] >= 0
 uATK: uDX[4,12] - ustk[4,12] + ucan[4,12] >= 0
 uATK: udw[3,12] - udx[4,12] = 0
 uATK: uDW[3,13] - udx[4,13] - ucan[4,13] >= 0
 uATK: - uDW[3,13] + uDX[4,13] + ucan[4,13] >= 0
 uATK: uDW[3,13] - uDX[4,13] + ustk[4,13] - 2 ucan[4,13] >= 0
 uATK: uDX[4,13] - ustk[4,13] + ucan[4,13] >= 0
 uATK: udw[3,13] - udx[4,13] = 0
 uATK: uDW[3,14] - udx[4,14] - ucan[4,14] >= 0
 uATK: - uDW[3,14] + uDX[4,14] + ucan[4,14] >= 0
 uATK: uDW[3,14] - uDX[4,14] + ustk[4,14] - 2 ucan[4,14] >= 0
 uATK: uDX[4,14] - ustk[4,14] + ucan[4,14] >= 0
 uATK: udw[3,14] - udx[4,14] = 0
 uATK: uDW[3,15] - udx[4,15] - ucan[4,15] >= 0
 uATK: - uDW[3,15] + uDX[4,15] + ucan[4,15] >= 0
 uATK: uDW[3,15] - uDX[4,15] + ustk[4,15] - 2 ucan[4,15] >= 0
 uATK: uDX[4,15] - ustk[4,15] + ucan[4,15] >= 0
 uATK: udw[3,15] - udx[4,15] = 0
 uATK: uDW[4,0] - udx[5,0] - ucan[5,0] >= 0
 uATK: - uDW[4,0] + uDX[5,0] + ucan[5,0] >= 0
 uATK: uDW[4,0] - uDX[5,0] + ustk[5,0] - 2 ucan[5,0] >= 0
 uATK: uDX[5,0] - ustk[5,0] + ucan[5,0] >= 0
 uATK: udw[4,0] - udx[5,0] = 0
 uATK: uDW[4,1] - udx[5,1] - ucan[5,1] >= 0
 uATK: - uDW[4,1] + uDX[5,1] + ucan[5,1] >= 0
 uATK: uDW[4,1] - uDX[5,1] + ustk[5,1] - 2 ucan[5,1] >= 0
 uATK: uDX[5,1] - ustk[5,1] + ucan[5,1] >= 0
 uATK: udw[4,1] - udx[5,1] = 0
 uATK: uDW[4,2] - udx[5,2] - ucan[5,2] >= 0
 uATK: - uDW[4,2] + uDX[5,2] + ucan[5,2] >= 0
 uATK: uDW[4,2] - uDX[5,2] + ustk[5,2] - 2 ucan[5,2] >= 0
 uATK: uDX[5,2] - ustk[5,2] + ucan[5,2] >= 0
 uATK: udw[4,2] - udx[5,2] = 0
 uATK: uDW[4,3] - udx[5,3] - ucan[5,3] >= 0
 uATK: - uDW[4,3] + uDX[5,3] + ucan[5,3] >= 0
 uATK: uDW[4,3] - uDX[5,3] + ustk[5,3] - 2 ucan[5,3] >= 0
 uATK: uDX[5,3] - ustk[5,3] + ucan[5,3] >= 0
 uATK: udw[4,3] - udx[5,3] = 0
 uATK: uDW[4,4] - udx[5,4] - ucan[5,4] >= 0
 uATK: - uDW[4,4] + uDX[5,4] + ucan[5,4] >= 0
 uATK: uDW[4,4] - uDX[5,4] + ustk[5,4] - 2 ucan[5,4] >= 0
 uATK: uDX[5,4] - ustk[5,4] + ucan[5,4] >= 0
 uATK: udw[4,4] - udx[5,4] = 0
 uATK: uDW[4,5] - udx[5,5] - ucan[5,5] >= 0
 uATK: - uDW[4,5] + uDX[5,5] + ucan[5,5] >= 0
 uATK: uDW[4,5] - uDX[5,5] + ustk[5,5] - 2 ucan[5,5] >= 0
 uATK: uDX[5,5] - ustk[5,5] + ucan[5,5] >= 0
 uATK: udw[4,5] - udx[5,5] = 0
 uATK: uDW[4,6] - udx[5,6] - ucan[5,6] >= 0
 uATK: - uDW[4,6] + uDX[5,6] + ucan[5,6] >= 0
 uATK: uDW[4,6] - uDX[5,6] + ustk[5,6] - 2 ucan[5,6] >= 0
 uATK: uDX[5,6] - ustk[5,6] + ucan[5,6] >= 0
 uATK: udw[4,6] - udx[5,6] = 0
 uATK: uDW[4,7] - udx[5,7] - ucan[5,7] >= 0
 uATK: - uDW[4,7] + uDX[5,7] + ucan[5,7] >= 0
 uATK: uDW[4,7] - uDX[5,7] + ustk[5,7] - 2 ucan[5,7] >= 0
 uATK: uDX[5,7] - ustk[5,7] + ucan[5,7] >= 0
 uATK: udw[4,7] - udx[5,7] = 0
 uATK: uDW[4,8] - udx[5,8] - ucan[5,8] >= 0
 uATK: - uDW[4,8] + uDX[5,8] + ucan[5,8] >= 0
 uATK: uDW[4,8] - uDX[5,8] + ustk[5,8] - 2 ucan[5,8] >= 0
 uATK: uDX[5,8] - ustk[5,8] + ucan[5,8] >= 0
 uATK: udw[4,8] - udx[5,8] = 0
 uATK: uDW[4,9] - udx[5,9] - ucan[5,9] >= 0
 uATK: - uDW[4,9] + uDX[5,9] + ucan[5,9] >= 0
 uATK: uDW[4,9] - uDX[5,9] + ustk[5,9] - 2 ucan[5,9] >= 0
 uATK: uDX[5,9] - ustk[5,9] + ucan[5,9] >= 0
 uATK: udw[4,9] - udx[5,9] = 0
 uATK: uDW[4,10] - udx[5,10] - ucan[5,10] >= 0
 uATK: - uDW[4,10] + uDX[5,10] + ucan[5,10] >= 0
 uATK: uDW[4,10] - uDX[5,10] + ustk[5,10] - 2 ucan[5,10] >= 0
 uATK: uDX[5,10] - ustk[5,10] + ucan[5,10] >= 0
 uATK: udw[4,10] - udx[5,10] = 0
 uATK: uDW[4,11] - udx[5,11] - ucan[5,11] >= 0
 uATK: - uDW[4,11] + uDX[5,11] + ucan[5,11] >= 0
 uATK: uDW[4,11] - uDX[5,11] + ustk[5,11] - 2 ucan[5,11] >= 0
 uATK: uDX[5,11] - ustk[5,11] + ucan[5,11] >= 0
 uATK: udw[4,11] - udx[5,11] = 0
 uATK: uDW[4,12] - udx[5,12] - ucan[5,12] >= 0
 uATK: - uDW[4,12] + uDX[5,12] + ucan[5,12] >= 0
 uATK: uDW[4,12] - uDX[5,12] + ustk[5,12] - 2 ucan[5,12] >= 0
 uATK: uDX[5,12] - ustk[5,12] + ucan[5,12] >= 0
 uATK: udw[4,12] - udx[5,12] = 0
 uATK: uDW[4,13] - udx[5,13] - ucan[5,13] >= 0
 uATK: - uDW[4,13] + uDX[5,13] + ucan[5,13] >= 0
 uATK: uDW[4,13] - uDX[5,13] + ustk[5,13] - 2 ucan[5,13] >= 0
 uATK: uDX[5,13] - ustk[5,13] + ucan[5,13] >= 0
 uATK: udw[4,13] - udx[5,13] = 0
 uATK: uDW[4,14] - udx[5,14] - ucan[5,14] >= 0
 uATK: - uDW[4,14] + uDX[5,14] + ucan[5,14] >= 0
 uATK: uDW[4,14] - uDX[5,14] + ustk[5,14] - 2 ucan[5,14] >= 0
 uATK: uDX[5,14] - ustk[5,14] + ucan[5,14] >= 0
 uATK: udw[4,14] - udx[5,14] = 0
 uATK: uDW[4,15] - udx[5,15] - ucan[5,15] >= 0
 uATK: - uDW[4,15] + uDX[5,15] + ucan[5,15] >= 0
 uATK: uDW[4,15] - uDX[5,15] + ustk[5,15] - 2 ucan[5,15] >= 0
 uATK: uDX[5,15] - ustk[5,15] + ucan[5,15] >= 0
 uATK: udw[4,15] - udx[5,15] = 0
 uATK: uDW[5,0] - udx[6,0] - ucan[6,0] >= 0
 uATK: - uDW[5,0] + uDX[6,0] + ucan[6,0] >= 0
 uATK: uDW[5,0] - uDX[6,0] + ustk[6,0] - 2 ucan[6,0] >= 0
 uATK: uDX[6,0] - ustk[6,0] + ucan[6,0] >= 0
 uATK: udw[5,0] - udx[6,0] = 0
 uATK: uDW[5,1] - udx[6,1] - ucan[6,1] >= 0
 uATK: - uDW[5,1] + uDX[6,1] + ucan[6,1] >= 0
 uATK: uDW[5,1] - uDX[6,1] + ustk[6,1] - 2 ucan[6,1] >= 0
 uATK: uDX[6,1] - ustk[6,1] + ucan[6,1] >= 0
 uATK: udw[5,1] - udx[6,1] = 0
 uATK: uDW[5,2] - udx[6,2] - ucan[6,2] >= 0
 uATK: - uDW[5,2] + uDX[6,2] + ucan[6,2] >= 0
 uATK: uDW[5,2] - uDX[6,2] + ustk[6,2] - 2 ucan[6,2] >= 0
 uATK: uDX[6,2] - ustk[6,2] + ucan[6,2] >= 0
 uATK: udw[5,2] - udx[6,2] = 0
 uATK: uDW[5,3] - udx[6,3] - ucan[6,3] >= 0
 uATK: - uDW[5,3] + uDX[6,3] + ucan[6,3] >= 0
 uATK: uDW[5,3] - uDX[6,3] + ustk[6,3] - 2 ucan[6,3] >= 0
 uATK: uDX[6,3] - ustk[6,3] + ucan[6,3] >= 0
 uATK: udw[5,3] - udx[6,3] = 0
 uATK: uDW[5,4] - udx[6,4] - ucan[6,4] >= 0
 uATK: - uDW[5,4] + uDX[6,4] + ucan[6,4] >= 0
 uATK: uDW[5,4] - uDX[6,4] + ustk[6,4] - 2 ucan[6,4] >= 0
 uATK: uDX[6,4] - ustk[6,4] + ucan[6,4] >= 0
 uATK: udw[5,4] - udx[6,4] = 0
 uATK: uDW[5,5] - udx[6,5] - ucan[6,5] >= 0
 uATK: - uDW[5,5] + uDX[6,5] + ucan[6,5] >= 0
 uATK: uDW[5,5] - uDX[6,5] + ustk[6,5] - 2 ucan[6,5] >= 0
 uATK: uDX[6,5] - ustk[6,5] + ucan[6,5] >= 0
 uATK: udw[5,5] - udx[6,5] = 0
 uATK: uDW[5,6] - udx[6,6] - ucan[6,6] >= 0
 uATK: - uDW[5,6] + uDX[6,6] + ucan[6,6] >= 0
 uATK: uDW[5,6] - uDX[6,6] + ustk[6,6] - 2 ucan[6,6] >= 0
 uATK: uDX[6,6] - ustk[6,6] + ucan[6,6] >= 0
 uATK: udw[5,6] - udx[6,6] = 0
 uATK: uDW[5,7] - udx[6,7] - ucan[6,7] >= 0
 uATK: - uDW[5,7] + uDX[6,7] + ucan[6,7] >= 0
 uATK: uDW[5,7] - uDX[6,7] + ustk[6,7] - 2 ucan[6,7] >= 0
 uATK: uDX[6,7] - ustk[6,7] + ucan[6,7] >= 0
 uATK: udw[5,7] - udx[6,7] = 0
 uATK: uDW[5,8] - udx[6,8] - ucan[6,8] >= 0
 uATK: - uDW[5,8] + uDX[6,8] + ucan[6,8] >= 0
 uATK: uDW[5,8] - uDX[6,8] + ustk[6,8] - 2 ucan[6,8] >= 0
 uATK: uDX[6,8] - ustk[6,8] + ucan[6,8] >= 0
 uATK: udw[5,8] - udx[6,8] = 0
 uATK: uDW[5,9] - udx[6,9] - ucan[6,9] >= 0
 uATK: - uDW[5,9] + uDX[6,9] + ucan[6,9] >= 0
 uATK: uDW[5,9] - uDX[6,9] + ustk[6,9] - 2 ucan[6,9] >= 0
 uATK: uDX[6,9] - ustk[6,9] + ucan[6,9] >= 0
 uATK: udw[5,9] - udx[6,9] = 0
 uATK: uDW[5,10] - udx[6,10] - ucan[6,10] >= 0
 uATK: - uDW[5,10] + uDX[6,10] + ucan[6,10] >= 0
 uATK: uDW[5,10] - uDX[6,10] + ustk[6,10] - 2 ucan[6,10] >= 0
 uATK: uDX[6,10] - ustk[6,10] + ucan[6,10] >= 0
 uATK: udw[5,10] - udx[6,10] = 0
 uATK: uDW[5,11] - udx[6,11] - ucan[6,11] >= 0
 uATK: - uDW[5,11] + uDX[6,11] + ucan[6,11] >= 0
 uATK: uDW[5,11] - uDX[6,11] + ustk[6,11] - 2 ucan[6,11] >= 0
 uATK: uDX[6,11] - ustk[6,11] + ucan[6,11] >= 0
 uATK: udw[5,11] - udx[6,11] = 0
 uATK: uDW[5,12] - udx[6,12] - ucan[6,12] >= 0
 uATK: - uDW[5,12] + uDX[6,12] + ucan[6,12] >= 0
 uATK: uDW[5,12] - uDX[6,12] + ustk[6,12] - 2 ucan[6,12] >= 0
 uATK: uDX[6,12] - ustk[6,12] + ucan[6,12] >= 0
 uATK: udw[5,12] - udx[6,12] = 0
 uATK: uDW[5,13] - udx[6,13] - ucan[6,13] >= 0
 uATK: - uDW[5,13] + uDX[6,13] + ucan[6,13] >= 0
 uATK: uDW[5,13] - uDX[6,13] + ustk[6,13] - 2 ucan[6,13] >= 0
 uATK: uDX[6,13] - ustk[6,13] + ucan[6,13] >= 0
 uATK: udw[5,13] - udx[6,13] = 0
 uATK: uDW[5,14] - udx[6,14] - ucan[6,14] >= 0
 uATK: - uDW[5,14] + uDX[6,14] + ucan[6,14] >= 0
 uATK: uDW[5,14] - uDX[6,14] + ustk[6,14] - 2 ucan[6,14] >= 0
 uATK: uDX[6,14] - ustk[6,14] + ucan[6,14] >= 0
 uATK: udw[5,14] - udx[6,14] = 0
 uATK: uDW[5,15] - udx[6,15] - ucan[6,15] >= 0
 uATK: - uDW[5,15] + uDX[6,15] + ucan[6,15] >= 0
 uATK: uDW[5,15] - uDX[6,15] + ustk[6,15] - 2 ucan[6,15] >= 0
 uATK: uDX[6,15] - ustk[6,15] + ucan[6,15] >= 0
 uATK: udw[5,15] - udx[6,15] = 0
 uATK: uDW[6,0] - udx[7,0] - ucan[7,0] >= 0
 uATK: - uDW[6,0] + uDX[7,0] + ucan[7,0] >= 0
 uATK: uDW[6,0] - uDX[7,0] + ustk[7,0] - 2 ucan[7,0] >= 0
 uATK: uDX[7,0] - ustk[7,0] + ucan[7,0] >= 0
 uATK: udw[6,0] - udx[7,0] = 0
 uATK: uDW[6,1] - udx[7,1] - ucan[7,1] >= 0
 uATK: - uDW[6,1] + uDX[7,1] + ucan[7,1] >= 0
 uATK: uDW[6,1] - uDX[7,1] + ustk[7,1] - 2 ucan[7,1] >= 0
 uATK: uDX[7,1] - ustk[7,1] + ucan[7,1] >= 0
 uATK: udw[6,1] - udx[7,1] = 0
 uATK: uDW[6,2] - udx[7,2] - ucan[7,2] >= 0
 uATK: - uDW[6,2] + uDX[7,2] + ucan[7,2] >= 0
 uATK: uDW[6,2] - uDX[7,2] + ustk[7,2] - 2 ucan[7,2] >= 0
 uATK: uDX[7,2] - ustk[7,2] + ucan[7,2] >= 0
 uATK: udw[6,2] - udx[7,2] = 0
 uATK: uDW[6,3] - udx[7,3] - ucan[7,3] >= 0
 uATK: - uDW[6,3] + uDX[7,3] + ucan[7,3] >= 0
 uATK: uDW[6,3] - uDX[7,3] + ustk[7,3] - 2 ucan[7,3] >= 0
 uATK: uDX[7,3] - ustk[7,3] + ucan[7,3] >= 0
 uATK: udw[6,3] - udx[7,3] = 0
 uATK: uDW[6,4] - udx[7,4] - ucan[7,4] >= 0
 uATK: - uDW[6,4] + uDX[7,4] + ucan[7,4] >= 0
 uATK: uDW[6,4] - uDX[7,4] + ustk[7,4] - 2 ucan[7,4] >= 0
 uATK: uDX[7,4] - ustk[7,4] + ucan[7,4] >= 0
 uATK: udw[6,4] - udx[7,4] = 0
 uATK: uDW[6,5] - udx[7,5] - ucan[7,5] >= 0
 uATK: - uDW[6,5] + uDX[7,5] + ucan[7,5] >= 0
 uATK: uDW[6,5] - uDX[7,5] + ustk[7,5] - 2 ucan[7,5] >= 0
 uATK: uDX[7,5] - ustk[7,5] + ucan[7,5] >= 0
 uATK: udw[6,5] - udx[7,5] = 0
 uATK: uDW[6,6] - udx[7,6] - ucan[7,6] >= 0
 uATK: - uDW[6,6] + uDX[7,6] + ucan[7,6] >= 0
 uATK: uDW[6,6] - uDX[7,6] + ustk[7,6] - 2 ucan[7,6] >= 0
 uATK: uDX[7,6] - ustk[7,6] + ucan[7,6] >= 0
 uATK: udw[6,6] - udx[7,6] = 0
 uATK: uDW[6,7] - udx[7,7] - ucan[7,7] >= 0
 uATK: - uDW[6,7] + uDX[7,7] + ucan[7,7] >= 0
 uATK: uDW[6,7] - uDX[7,7] + ustk[7,7] - 2 ucan[7,7] >= 0
 uATK: uDX[7,7] - ustk[7,7] + ucan[7,7] >= 0
 uATK: udw[6,7] - udx[7,7] = 0
 uATK: uDW[6,8] - udx[7,8] - ucan[7,8] >= 0
 uATK: - uDW[6,8] + uDX[7,8] + ucan[7,8] >= 0
 uATK: uDW[6,8] - uDX[7,8] + ustk[7,8] - 2 ucan[7,8] >= 0
 uATK: uDX[7,8] - ustk[7,8] + ucan[7,8] >= 0
 uATK: udw[6,8] - udx[7,8] = 0
 uATK: uDW[6,9] - udx[7,9] - ucan[7,9] >= 0
 uATK: - uDW[6,9] + uDX[7,9] + ucan[7,9] >= 0
 uATK: uDW[6,9] - uDX[7,9] + ustk[7,9] - 2 ucan[7,9] >= 0
 uATK: uDX[7,9] - ustk[7,9] + ucan[7,9] >= 0
 uATK: udw[6,9] - udx[7,9] = 0
 uATK: uDW[6,10] - udx[7,10] - ucan[7,10] >= 0
 uATK: - uDW[6,10] + uDX[7,10] + ucan[7,10] >= 0
 uATK: uDW[6,10] - uDX[7,10] + ustk[7,10] - 2 ucan[7,10] >= 0
 uATK: uDX[7,10] - ustk[7,10] + ucan[7,10] >= 0
 uATK: udw[6,10] - udx[7,10] = 0
 uATK: uDW[6,11] - udx[7,11] - ucan[7,11] >= 0
 uATK: - uDW[6,11] + uDX[7,11] + ucan[7,11] >= 0
 uATK: uDW[6,11] - uDX[7,11] + ustk[7,11] - 2 ucan[7,11] >= 0
 uATK: uDX[7,11] - ustk[7,11] + ucan[7,11] >= 0
 uATK: udw[6,11] - udx[7,11] = 0
 uATK: uDW[6,12] - udx[7,12] - ucan[7,12] >= 0
 uATK: - uDW[6,12] + uDX[7,12] + ucan[7,12] >= 0
 uATK: uDW[6,12] - uDX[7,12] + ustk[7,12] - 2 ucan[7,12] >= 0
 uATK: uDX[7,12] - ustk[7,12] + ucan[7,12] >= 0
 uATK: udw[6,12] - udx[7,12] = 0
 uATK: uDW[6,13] - udx[7,13] - ucan[7,13] >= 0
 uATK: - uDW[6,13] + uDX[7,13] + ucan[7,13] >= 0
 uATK: uDW[6,13] - uDX[7,13] + ustk[7,13] - 2 ucan[7,13] >= 0
 uATK: uDX[7,13] - ustk[7,13] + ucan[7,13] >= 0
 uATK: udw[6,13] - udx[7,13] = 0
 uATK: uDW[6,14] - udx[7,14] - ucan[7,14] >= 0
 uATK: - uDW[6,14] + uDX[7,14] + ucan[7,14] >= 0
 uATK: uDW[6,14] - uDX[7,14] + ustk[7,14] - 2 ucan[7,14] >= 0
 uATK: uDX[7,14] - ustk[7,14] + ucan[7,14] >= 0
 uATK: udw[6,14] - udx[7,14] = 0
 uATK: uDW[6,15] - udx[7,15] - ucan[7,15] >= 0
 uATK: - uDW[6,15] + uDX[7,15] + ucan[7,15] >= 0
 uATK: uDW[6,15] - uDX[7,15] + ustk[7,15] - 2 ucan[7,15] >= 0
 uATK: uDX[7,15] - ustk[7,15] + ucan[7,15] >= 0
 uATK: udw[6,15] - udx[7,15] = 0
 bSC: uDX[0,0] - uDY[0,0] = 0
 bSC: udx[0,0] - uDY[0,0] = 0
 bSC: uDX[0,1] - uDY[0,1] = 0
 bSC: udx[0,1] - uDY[0,1] = 0
 bSC: uDX[0,2] - uDY[0,2] = 0
 bSC: udx[0,2] - uDY[0,2] = 0
 bSC: uDX[0,3] - uDY[0,3] = 0
 bSC: udx[0,3] - uDY[0,3] = 0
 bSC: uDX[0,4] - uDY[0,4] = 0
 bSC: udx[0,4] - uDY[0,4] = 0
 bSC: uDX[0,5] - uDY[0,5] = 0
 bSC: udx[0,5] - uDY[0,5] = 0
 bSC: uDX[0,6] - uDY[0,6] = 0
 bSC: udx[0,6] - uDY[0,6] = 0
 bSC: uDX[0,7] - uDY[0,7] = 0
 bSC: udx[0,7] - uDY[0,7] = 0
 bSC: uDX[0,8] - uDY[0,8] = 0
 bSC: udx[0,8] - uDY[0,8] = 0
 bSC: uDX[0,9] - uDY[0,9] = 0
 bSC: udx[0,9] - uDY[0,9] = 0
 bSC: uDX[0,10] - uDY[0,10] = 0
 bSC: udx[0,10] - uDY[0,10] = 0
 bSC: uDX[0,11] - uDY[0,11] = 0
 bSC: udx[0,11] - uDY[0,11] = 0
 bSC: uDX[0,12] - uDY[0,12] = 0
 bSC: udx[0,12] - uDY[0,12] = 0
 bSC: uDX[0,13] - uDY[0,13] = 0
 bSC: udx[0,13] - uDY[0,13] = 0
 bSC: uDX[0,14] - uDY[0,14] = 0
 bSC: udx[0,14] - uDY[0,14] = 0
 bSC: uDX[0,15] - uDY[0,15] = 0
 bSC: udx[0,15] - uDY[0,15] = 0
 uSC: - uDX[1,0] + uDY[1,0] = 0
 uSC: - uDX[1,0] + udy[1,0] = 0
 uSC: - uDX[1,1] + uDY[1,1] = 0
 uSC: - uDX[1,1] + udy[1,1] = 0
 uSC: - uDX[1,2] + uDY[1,2] = 0
 uSC: - uDX[1,2] + udy[1,2] = 0
 uSC: - uDX[1,3] + uDY[1,3] = 0
 uSC: - uDX[1,3] + udy[1,3] = 0
 uSC: - uDX[1,4] + uDY[1,4] = 0
 uSC: - uDX[1,4] + udy[1,4] = 0
 uSC: - uDX[1,5] + uDY[1,5] = 0
 uSC: - uDX[1,5] + udy[1,5] = 0
 uSC: - uDX[1,6] + uDY[1,6] = 0
 uSC: - uDX[1,6] + udy[1,6] = 0
 uSC: - uDX[1,7] + uDY[1,7] = 0
 uSC: - uDX[1,7] + udy[1,7] = 0
 uSC: - uDX[1,8] + uDY[1,8] = 0
 uSC: - uDX[1,8] + udy[1,8] = 0
 uSC: - uDX[1,9] + uDY[1,9] = 0
 uSC: - uDX[1,9] + udy[1,9] = 0
 uSC: - uDX[1,10] + uDY[1,10] = 0
 uSC: - uDX[1,10] + udy[1,10] = 0
 uSC: - uDX[1,11] + uDY[1,11] = 0
 uSC: - uDX[1,11] + udy[1,11] = 0
 uSC: - uDX[1,12] + uDY[1,12] = 0
 uSC: - uDX[1,12] + udy[1,12] = 0
 uSC: - uDX[1,13] + uDY[1,13] = 0
 uSC: - uDX[1,13] + udy[1,13] = 0
 uSC: - uDX[1,14] + uDY[1,14] = 0
 uSC: - uDX[1,14] + udy[1,14] = 0
 uSC: - uDX[1,15] + uDY[1,15] = 0
 uSC: - uDX[1,15] + udy[1,15] = 0
 uSC: - uDX[2,0] + uDY[2,0] = 0
 uSC: - uDX[2,0] + udy[2,0] = 0
 uSC: - uDX[2,1] + uDY[2,1] = 0
 uSC: - uDX[2,1] + udy[2,1] = 0
 uSC: - uDX[2,2] + uDY[2,2] = 0
 uSC: - uDX[2,2] + udy[2,2] = 0
 uSC: - uDX[2,3] + uDY[2,3] = 0
 uSC: - uDX[2,3] + udy[2,3] = 0
 uSC: - uDX[2,4] + uDY[2,4] = 0
 uSC: - uDX[2,4] + udy[2,4] = 0
 uSC: - uDX[2,5] + uDY[2,5] = 0
 uSC: - uDX[2,5] + udy[2,5] = 0
 uSC: - uDX[2,6] + uDY[2,6] = 0
 uSC: - uDX[2,6] + udy[2,6] = 0
 uSC: - uDX[2,7] + uDY[2,7] = 0
 uSC: - uDX[2,7] + udy[2,7] = 0
 uSC: - uDX[2,8] + uDY[2,8] = 0
 uSC: - uDX[2,8] + udy[2,8] = 0
 uSC: - uDX[2,9] + uDY[2,9] = 0
 uSC: - uDX[2,9] + udy[2,9] = 0
 uSC: - uDX[2,10] + uDY[2,10] = 0
 uSC: - uDX[2,10] + udy[2,10] = 0
 uSC: - uDX[2,11] + uDY[2,11] = 0
 uSC: - uDX[2,11] + udy[2,11] = 0
 uSC: - uDX[2,12] + uDY[2,12] = 0
 uSC: - uDX[2,12] + udy[2,12] = 0
 uSC: - uDX[2,13] + uDY[2,13] = 0
 uSC: - uDX[2,13] + udy[2,13] = 0
 uSC: - uDX[2,14] + uDY[2,14] = 0
 uSC: - uDX[2,14] + udy[2,14] = 0
 uSC: - uDX[2,15] + uDY[2,15] = 0
 uSC: - uDX[2,15] + udy[2,15] = 0
 uSC: - uDX[3,0] + uDY[3,0] = 0
 uSC: - uDX[3,0] + udy[3,0] = 0
 uSC: - uDX[3,1] + uDY[3,1] = 0
 uSC: - uDX[3,1] + udy[3,1] = 0
 uSC: - uDX[3,2] + uDY[3,2] = 0
 uSC: - uDX[3,2] + udy[3,2] = 0
 uSC: - uDX[3,3] + uDY[3,3] = 0
 uSC: - uDX[3,3] + udy[3,3] = 0
 uSC: - uDX[3,4] + uDY[3,4] = 0
 uSC: - uDX[3,4] + udy[3,4] = 0
 uSC: - uDX[3,5] + uDY[3,5] = 0
 uSC: - uDX[3,5] + udy[3,5] = 0
 uSC: - uDX[3,6] + uDY[3,6] = 0
 uSC: - uDX[3,6] + udy[3,6] = 0
 uSC: - uDX[3,7] + uDY[3,7] = 0
 uSC: - uDX[3,7] + udy[3,7] = 0
 uSC: - uDX[3,8] + uDY[3,8] = 0
 uSC: - uDX[3,8] + udy[3,8] = 0
 uSC: - uDX[3,9] + uDY[3,9] = 0
 uSC: - uDX[3,9] + udy[3,9] = 0
 uSC: - uDX[3,10] + uDY[3,10] = 0
 uSC: - uDX[3,10] + udy[3,10] = 0
 uSC: - uDX[3,11] + uDY[3,11] = 0
 uSC: - uDX[3,11] + udy[3,11] = 0
 uSC: - uDX[3,12] + uDY[3,12] = 0
 uSC: - uDX[3,12] + udy[3,12] = 0
 uSC: - uDX[3,13] + uDY[3,13] = 0
 uSC: - uDX[3,13] + udy[3,13] = 0
 uSC: - uDX[3,14] + uDY[3,14] = 0
 uSC: - uDX[3,14] + udy[3,14] = 0
 uSC: - uDX[3,15] + uDY[3,15] = 0
 uSC: - uDX[3,15] + udy[3,15] = 0
 uSC: - uDX[4,0] + uDY[4,0] = 0
 uSC: - uDX[4,0] + udy[4,0] = 0
 uSC: - uDX[4,1] + uDY[4,1] = 0
 uSC: - uDX[4,1] + udy[4,1] = 0
 uSC: - uDX[4,2] + uDY[4,2] = 0
 uSC: - uDX[4,2] + udy[4,2] = 0
 uSC: - uDX[4,3] + uDY[4,3] = 0
 uSC: - uDX[4,3] + udy[4,3] = 0
 uSC: - uDX[4,4] + uDY[4,4] = 0
 uSC: - uDX[4,4] + udy[4,4] = 0
 uSC: - uDX[4,5] + uDY[4,5] = 0
 uSC: - uDX[4,5] + udy[4,5] = 0
 uSC: - uDX[4,6] + uDY[4,6] = 0
 uSC: - uDX[4,6] + udy[4,6] = 0
 uSC: - uDX[4,7] + uDY[4,7] = 0
 uSC: - uDX[4,7] + udy[4,7] = 0
 uSC: - uDX[4,8] + uDY[4,8] = 0
 uSC: - uDX[4,8] + udy[4,8] = 0
 uSC: - uDX[4,9] + uDY[4,9] = 0
 uSC: - uDX[4,9] + udy[4,9] = 0
 uSC: - uDX[4,10] + uDY[4,10] = 0
 uSC: - uDX[4,10] + udy[4,10] = 0
 uSC: - uDX[4,11] + uDY[4,11] = 0
 uSC: - uDX[4,11] + udy[4,11] = 0
 uSC: - uDX[4,12] + uDY[4,12] = 0
 uSC: - uDX[4,12] + udy[4,12] = 0
 uSC: - uDX[4,13] + uDY[4,13] = 0
 uSC: - uDX[4,13] + udy[4,13] = 0
 uSC: - uDX[4,14] + uDY[4,14] = 0
 uSC: - uDX[4,14] + udy[4,14] = 0
 uSC: - uDX[4,15] + uDY[4,15] = 0
 uSC: - uDX[4,15] + udy[4,15] = 0
 uSC: - uDX[5,0] + uDY[5,0] = 0
 uSC: - uDX[5,0] + udy[5,0] = 0
 uSC: - uDX[5,1] + uDY[5,1] = 0
 uSC: - uDX[5,1] + udy[5,1] = 0
 uSC: - uDX[5,2] + uDY[5,2] = 0
 uSC: - uDX[5,2] + udy[5,2] = 0
 uSC: - uDX[5,3] + uDY[5,3] = 0
 uSC: - uDX[5,3] + udy[5,3] = 0
 uSC: - uDX[5,4] + uDY[5,4] = 0
 uSC: - uDX[5,4] + udy[5,4] = 0
 uSC: - uDX[5,5] + uDY[5,5] = 0
 uSC: - uDX[5,5] + udy[5,5] = 0
 uSC: - uDX[5,6] + uDY[5,6] = 0
 uSC: - uDX[5,6] + udy[5,6] = 0
 uSC: - uDX[5,7] + uDY[5,7] = 0
 uSC: - uDX[5,7] + udy[5,7] = 0
 uSC: - uDX[5,8] + uDY[5,8] = 0
 uSC: - uDX[5,8] + udy[5,8] = 0
 uSC: - uDX[5,9] + uDY[5,9] = 0
 uSC: - uDX[5,9] + udy[5,9] = 0
 uSC: - uDX[5,10] + uDY[5,10] = 0
 uSC: - uDX[5,10] + udy[5,10] = 0
 uSC: - uDX[5,11] + uDY[5,11] = 0
 uSC: - uDX[5,11] + udy[5,11] = 0
 uSC: - uDX[5,12] + uDY[5,12] = 0
 uSC: - uDX[5,12] + udy[5,12] = 0
 uSC: - uDX[5,13] + uDY[5,13] = 0
 uSC: - uDX[5,13] + udy[5,13] = 0
 uSC: - uDX[5,14] + uDY[5,14] = 0
 uSC: - uDX[5,14] + udy[5,14] = 0
 uSC: - uDX[5,15] + uDY[5,15] = 0
 uSC: - uDX[5,15] + udy[5,15] = 0
 uSC: - uDX[6,0] + uDY[6,0] = 0
 uSC: - uDX[6,0] + udy[6,0] = 0
 uSC: - uDX[6,1] + uDY[6,1] = 0
 uSC: - uDX[6,1] + udy[6,1] = 0
 uSC: - uDX[6,2] + uDY[6,2] = 0
 uSC: - uDX[6,2] + udy[6,2] = 0
 uSC: - uDX[6,3] + uDY[6,3] = 0
 uSC: - uDX[6,3] + udy[6,3] = 0
 uSC: - uDX[6,4] + uDY[6,4] = 0
 uSC: - uDX[6,4] + udy[6,4] = 0
 uSC: - uDX[6,5] + uDY[6,5] = 0
 uSC: - uDX[6,5] + udy[6,5] = 0
 uSC: - uDX[6,6] + uDY[6,6] = 0
 uSC: - uDX[6,6] + udy[6,6] = 0
 uSC: - uDX[6,7] + uDY[6,7] = 0
 uSC: - uDX[6,7] + udy[6,7] = 0
 uSC: - uDX[6,8] + uDY[6,8] = 0
 uSC: - uDX[6,8] + udy[6,8] = 0
 uSC: - uDX[6,9] + uDY[6,9] = 0
 uSC: - uDX[6,9] + udy[6,9] = 0
 uSC: - uDX[6,10] + uDY[6,10] = 0
 uSC: - uDX[6,10] + udy[6,10] = 0
 uSC: - uDX[6,11] + uDY[6,11] = 0
 uSC: - uDX[6,11] + udy[6,11] = 0
 uSC: - uDX[6,12] + uDY[6,12] = 0
 uSC: - uDX[6,12] + udy[6,12] = 0
 uSC: - uDX[6,13] + uDY[6,13] = 0
 uSC: - uDX[6,13] + udy[6,13] = 0
 uSC: - uDX[6,14] + uDY[6,14] = 0
 uSC: - uDX[6,14] + udy[6,14] = 0
 uSC: - uDX[6,15] + uDY[6,15] = 0
 uSC: - uDX[6,15] + udy[6,15] = 0
 uSC: - uDX[7,0] + uDY[7,0] = 0
 uSC: - uDX[7,0] + udy[7,0] = 0
 uSC: - uDX[7,1] + uDY[7,1] = 0
 uSC: - uDX[7,1] + udy[7,1] = 0
 uSC: - uDX[7,2] + uDY[7,2] = 0
 uSC: - uDX[7,2] + udy[7,2] = 0
 uSC: - uDX[7,3] + uDY[7,3] = 0
 uSC: - uDX[7,3] + udy[7,3] = 0
 uSC: - uDX[7,4] + uDY[7,4] = 0
 uSC: - uDX[7,4] + udy[7,4] = 0
 uSC: - uDX[7,5] + uDY[7,5] = 0
 uSC: - uDX[7,5] + udy[7,5] = 0
 uSC: - uDX[7,6] + uDY[7,6] = 0
 uSC: - uDX[7,6] + udy[7,6] = 0
 uSC: - uDX[7,7] + uDY[7,7] = 0
 uSC: - uDX[7,7] + udy[7,7] = 0
 uSC: - uDX[7,8] + uDY[7,8] = 0
 uSC: - uDX[7,8] + udy[7,8] = 0
 uSC: - uDX[7,9] + uDY[7,9] = 0
 uSC: - uDX[7,9] + udy[7,9] = 0
 uSC: - uDX[7,10] + uDY[7,10] = 0
 uSC: - uDX[7,10] + udy[7,10] = 0
 uSC: - uDX[7,11] + uDY[7,11] = 0
 uSC: - uDX[7,11] + udy[7,11] = 0
 uSC: - uDX[7,12] + uDY[7,12] = 0
 uSC: - uDX[7,12] + udy[7,12] = 0
 uSC: - uDX[7,13] + uDY[7,13] = 0
 uSC: - uDX[7,13] + udy[7,13] = 0
 uSC: - uDX[7,14] + uDY[7,14] = 0
 uSC: - uDX[7,14] + udy[7,14] = 0
 uSC: - uDX[7,15] + uDY[7,15] = 0
 uSC: - uDX[7,15] + udy[7,15] = 0
 uSR: uDY[0,0] - uDZ[0,0] = 0
 uSR: udy[0,0] - udz[0,0] = 0
 uSR: uDY[0,5] - uDZ[0,1] = 0
 uSR: udy[0,5] - udz[0,1] = 0
 uSR: uDY[0,10] - uDZ[0,2] = 0
 uSR: udy[0,10] - udz[0,2] = 0
 uSR: uDY[0,15] - uDZ[0,3] = 0
 uSR: udy[0,15] - udz[0,3] = 0
 uSR: uDY[0,4] - uDZ[0,4] = 0
 uSR: udy[0,4] - udz[0,4] = 0
 uSR: uDY[0,9] - uDZ[0,5] = 0
 uSR: udy[0,9] - udz[0,5] = 0
 uSR: uDY[0,14] - uDZ[0,6] = 0
 uSR: udy[0,14] - udz[0,6] = 0
 uSR: uDY[0,3] - uDZ[0,7] = 0
 uSR: udy[0,3] - udz[0,7] = 0
 uSR: uDY[0,8] - uDZ[0,8] = 0
 uSR: udy[0,8] - udz[0,8] = 0
 uSR: uDY[0,13] - uDZ[0,9] = 0
 uSR: udy[0,13] - udz[0,9] = 0
 uSR: uDY[0,2] - uDZ[0,10] = 0
 uSR: udy[0,2] - udz[0,10] = 0
 uSR: uDY[0,7] - uDZ[0,11] = 0
 uSR: udy[0,7] - udz[0,11] = 0
 uSR: uDY[0,12] - uDZ[0,12] = 0
 uSR: udy[0,12] - udz[0,12] = 0
 uSR: uDY[0,1] - uDZ[0,13] = 0
 uSR: udy[0,1] - udz[0,13] = 0
 uSR: uDY[0,6] - uDZ[0,14] = 0
 uSR: udy[0,6] - udz[0,14] = 0
 uSR: uDY[0,11] - uDZ[0,15] = 0
 uSR: udy[0,11] - udz[0,15] = 0
 uSR: uDY[1,0] - uDZ[1,0] = 0
 uSR: udy[1,0] - udz[1,0] = 0
 uSR: uDY[1,5] - uDZ[1,1] = 0
 uSR: udy[1,5] - udz[1,1] = 0
 uSR: uDY[1,10] - uDZ[1,2] = 0
 uSR: udy[1,10] - udz[1,2] = 0
 uSR: uDY[1,15] - uDZ[1,3] = 0
 uSR: udy[1,15] - udz[1,3] = 0
 uSR: uDY[1,4] - uDZ[1,4] = 0
 uSR: udy[1,4] - udz[1,4] = 0
 uSR: uDY[1,9] - uDZ[1,5] = 0
 uSR: udy[1,9] - udz[1,5] = 0
 uSR: uDY[1,14] - uDZ[1,6] = 0
 uSR: udy[1,14] - udz[1,6] = 0
 uSR: uDY[1,3] - uDZ[1,7] = 0
 uSR: udy[1,3] - udz[1,7] = 0
 uSR: uDY[1,8] - uDZ[1,8] = 0
 uSR: udy[1,8] - udz[1,8] = 0
 uSR: uDY[1,13] - uDZ[1,9] = 0
 uSR: udy[1,13] - udz[1,9] = 0
 uSR: uDY[1,2] - uDZ[1,10] = 0
 uSR: udy[1,2] - udz[1,10] = 0
 uSR: uDY[1,7] - uDZ[1,11] = 0
 uSR: udy[1,7] - udz[1,11] = 0
 uSR: uDY[1,12] - uDZ[1,12] = 0
 uSR: udy[1,12] - udz[1,12] = 0
 uSR: uDY[1,1] - uDZ[1,13] = 0
 uSR: udy[1,1] - udz[1,13] = 0
 uSR: uDY[1,6] - uDZ[1,14] = 0
 uSR: udy[1,6] - udz[1,14] = 0
 uSR: uDY[1,11] - uDZ[1,15] = 0
 uSR: udy[1,11] - udz[1,15] = 0
 uSR: uDY[2,0] - uDZ[2,0] = 0
 uSR: udy[2,0] - udz[2,0] = 0
 uSR: uDY[2,5] - uDZ[2,1] = 0
 uSR: udy[2,5] - udz[2,1] = 0
 uSR: uDY[2,10] - uDZ[2,2] = 0
 uSR: udy[2,10] - udz[2,2] = 0
 uSR: uDY[2,15] - uDZ[2,3] = 0
 uSR: udy[2,15] - udz[2,3] = 0
 uSR: uDY[2,4] - uDZ[2,4] = 0
 uSR: udy[2,4] - udz[2,4] = 0
 uSR: uDY[2,9] - uDZ[2,5] = 0
 uSR: udy[2,9] - udz[2,5] = 0
 uSR: uDY[2,14] - uDZ[2,6] = 0
 uSR: udy[2,14] - udz[2,6] = 0
 uSR: uDY[2,3] - uDZ[2,7] = 0
 uSR: udy[2,3] - udz[2,7] = 0
 uSR: uDY[2,8] - uDZ[2,8] = 0
 uSR: udy[2,8] - udz[2,8] = 0
 uSR: uDY[2,13] - uDZ[2,9] = 0
 uSR: udy[2,13] - udz[2,9] = 0
 uSR: uDY[2,2] - uDZ[2,10] = 0
 uSR: udy[2,2] - udz[2,10] = 0
 uSR: uDY[2,7] - uDZ[2,11] = 0
 uSR: udy[2,7] - udz[2,11] = 0
 uSR: uDY[2,12] - uDZ[2,12] = 0
 uSR: udy[2,12] - udz[2,12] = 0
 uSR: uDY[2,1] - uDZ[2,13] = 0
 uSR: udy[2,1] - udz[2,13] = 0
 uSR: uDY[2,6] - uDZ[2,14] = 0
 uSR: udy[2,6] - udz[2,14] = 0
 uSR: uDY[2,11] - uDZ[2,15] = 0
 uSR: udy[2,11] - udz[2,15] = 0
 uSR: uDY[3,0] - uDZ[3,0] = 0
 uSR: udy[3,0] - udz[3,0] = 0
 uSR: uDY[3,5] - uDZ[3,1] = 0
 uSR: udy[3,5] - udz[3,1] = 0
 uSR: uDY[3,10] - uDZ[3,2] = 0
 uSR: udy[3,10] - udz[3,2] = 0
 uSR: uDY[3,15] - uDZ[3,3] = 0
 uSR: udy[3,15] - udz[3,3] = 0
 uSR: uDY[3,4] - uDZ[3,4] = 0
 uSR: udy[3,4] - udz[3,4] = 0
 uSR: uDY[3,9] - uDZ[3,5] = 0
 uSR: udy[3,9] - udz[3,5] = 0
 uSR: uDY[3,14] - uDZ[3,6] = 0
 uSR: udy[3,14] - udz[3,6] = 0
 uSR: uDY[3,3] - uDZ[3,7] = 0
 uSR: udy[3,3] - udz[3,7] = 0
 uSR: uDY[3,8] - uDZ[3,8] = 0
 uSR: udy[3,8] - udz[3,8] = 0
 uSR: uDY[3,13] - uDZ[3,9] = 0
 uSR: udy[3,13] - udz[3,9] = 0
 uSR: uDY[3,2] - uDZ[3,10] = 0
 uSR: udy[3,2] - udz[3,10] = 0
 uSR: uDY[3,7] - uDZ[3,11] = 0
 uSR: udy[3,7] - udz[3,11] = 0
 uSR: uDY[3,12] - uDZ[3,12] = 0
 uSR: udy[3,12] - udz[3,12] = 0
 uSR: uDY[3,1] - uDZ[3,13] = 0
 uSR: udy[3,1] - udz[3,13] = 0
 uSR: uDY[3,6] - uDZ[3,14] = 0
 uSR: udy[3,6] - udz[3,14] = 0
 uSR: uDY[3,11] - uDZ[3,15] = 0
 uSR: udy[3,11] - udz[3,15] = 0
 uSR: uDY[4,0] - uDZ[4,0] = 0
 uSR: udy[4,0] - udz[4,0] = 0
 uSR: uDY[4,5] - uDZ[4,1] = 0
 uSR: udy[4,5] - udz[4,1] = 0
 uSR: uDY[4,10] - uDZ[4,2] = 0
 uSR: udy[4,10] - udz[4,2] = 0
 uSR: uDY[4,15] - uDZ[4,3] = 0
 uSR: udy[4,15] - udz[4,3] = 0
 uSR: uDY[4,4] - uDZ[4,4] = 0
 uSR: udy[4,4] - udz[4,4] = 0
 uSR: uDY[4,9] - uDZ[4,5] = 0
 uSR: udy[4,9] - udz[4,5] = 0
 uSR: uDY[4,14] - uDZ[4,6] = 0
 uSR: udy[4,14] - udz[4,6] = 0
 uSR: uDY[4,3] - uDZ[4,7] = 0
 uSR: udy[4,3] - udz[4,7] = 0
 uSR: uDY[4,8] - uDZ[4,8] = 0
 uSR: udy[4,8] - udz[4,8] = 0
 uSR: uDY[4,13] - uDZ[4,9] = 0
 uSR: udy[4,13] - udz[4,9] = 0
 uSR: uDY[4,2] - uDZ[4,10] = 0
 uSR: udy[4,2] - udz[4,10] = 0
 uSR: uDY[4,7] - uDZ[4,11] = 0
 uSR: udy[4,7] - udz[4,11] = 0
 uSR: uDY[4,12] - uDZ[4,12] = 0
 uSR: udy[4,12] - udz[4,12] = 0
 uSR: uDY[4,1] - uDZ[4,13] = 0
 uSR: udy[4,1] - udz[4,13] = 0
 uSR: uDY[4,6] - uDZ[4,14] = 0
 uSR: udy[4,6] - udz[4,14] = 0
 uSR: uDY[4,11] - uDZ[4,15] = 0
 uSR: udy[4,11] - udz[4,15] = 0
 uSR: uDY[5,0] - uDZ[5,0] = 0
 uSR: udy[5,0] - udz[5,0] = 0
 uSR: uDY[5,5] - uDZ[5,1] = 0
 uSR: udy[5,5] - udz[5,1] = 0
 uSR: uDY[5,10] - uDZ[5,2] = 0
 uSR: udy[5,10] - udz[5,2] = 0
 uSR: uDY[5,15] - uDZ[5,3] = 0
 uSR: udy[5,15] - udz[5,3] = 0
 uSR: uDY[5,4] - uDZ[5,4] = 0
 uSR: udy[5,4] - udz[5,4] = 0
 uSR: uDY[5,9] - uDZ[5,5] = 0
 uSR: udy[5,9] - udz[5,5] = 0
 uSR: uDY[5,14] - uDZ[5,6] = 0
 uSR: udy[5,14] - udz[5,6] = 0
 uSR: uDY[5,3] - uDZ[5,7] = 0
 uSR: udy[5,3] - udz[5,7] = 0
 uSR: uDY[5,8] - uDZ[5,8] = 0
 uSR: udy[5,8] - udz[5,8] = 0
 uSR: uDY[5,13] - uDZ[5,9] = 0
 uSR: udy[5,13] - udz[5,9] = 0
 uSR: uDY[5,2] - uDZ[5,10] = 0
 uSR: udy[5,2] - udz[5,10] = 0
 uSR: uDY[5,7] - uDZ[5,11] = 0
 uSR: udy[5,7] - udz[5,11] = 0
 uSR: uDY[5,12] - uDZ[5,12] = 0
 uSR: udy[5,12] - udz[5,12] = 0
 uSR: uDY[5,1] - uDZ[5,13] = 0
 uSR: udy[5,1] - udz[5,13] = 0
 uSR: uDY[5,6] - uDZ[5,14] = 0
 uSR: udy[5,6] - udz[5,14] = 0
 uSR: uDY[5,11] - uDZ[5,15] = 0
 uSR: udy[5,11] - udz[5,15] = 0
 uSR: uDY[6,0] - uDZ[6,0] = 0
 uSR: udy[6,0] - udz[6,0] = 0
 uSR: uDY[6,5] - uDZ[6,1] = 0
 uSR: udy[6,5] - udz[6,1] = 0
 uSR: uDY[6,10] - uDZ[6,2] = 0
 uSR: udy[6,10] - udz[6,2] = 0
 uSR: uDY[6,15] - uDZ[6,3] = 0
 uSR: udy[6,15] - udz[6,3] = 0
 uSR: uDY[6,4] - uDZ[6,4] = 0
 uSR: udy[6,4] - udz[6,4] = 0
 uSR: uDY[6,9] - uDZ[6,5] = 0
 uSR: udy[6,9] - udz[6,5] = 0
 uSR: uDY[6,14] - uDZ[6,6] = 0
 uSR: udy[6,14] - udz[6,6] = 0
 uSR: uDY[6,3] - uDZ[6,7] = 0
 uSR: udy[6,3] - udz[6,7] = 0
 uSR: uDY[6,8] - uDZ[6,8] = 0
 uSR: udy[6,8] - udz[6,8] = 0
 uSR: uDY[6,13] - uDZ[6,9] = 0
 uSR: udy[6,13] - udz[6,9] = 0
 uSR: uDY[6,2] - uDZ[6,10] = 0
 uSR: udy[6,2] - udz[6,10] = 0
 uSR: uDY[6,7] - uDZ[6,11] = 0
 uSR: udy[6,7] - udz[6,11] = 0
 uSR: uDY[6,12] - uDZ[6,12] = 0
 uSR: udy[6,12] - udz[6,12] = 0
 uSR: uDY[6,1] - uDZ[6,13] = 0
 uSR: udy[6,1] - udz[6,13] = 0
 uSR: uDY[6,6] - uDZ[6,14] = 0
 uSR: udy[6,6] - udz[6,14] = 0
 uSR: uDY[6,11] - uDZ[6,15] = 0
 uSR: udy[6,11] - udz[6,15] = 0
 uSR: uDY[7,0] - uDZ[7,0] = 0
 uSR: udy[7,0] - udz[7,0] = 0
 uSR: uDY[7,5] - uDZ[7,1] = 0
 uSR: udy[7,5] - udz[7,1] = 0
 uSR: uDY[7,10] - uDZ[7,2] = 0
 uSR: udy[7,10] - udz[7,2] = 0
 uSR: uDY[7,15] - uDZ[7,3] = 0
 uSR: udy[7,15] - udz[7,3] = 0
 uSR: uDY[7,4] - uDZ[7,4] = 0
 uSR: udy[7,4] - udz[7,4] = 0
 uSR: uDY[7,9] - uDZ[7,5] = 0
 uSR: udy[7,9] - udz[7,5] = 0
 uSR: uDY[7,14] - uDZ[7,6] = 0
 uSR: udy[7,14] - udz[7,6] = 0
 uSR: uDY[7,3] - uDZ[7,7] = 0
 uSR: udy[7,3] - udz[7,7] = 0
 uSR: uDY[7,8] - uDZ[7,8] = 0
 uSR: udy[7,8] - udz[7,8] = 0
 uSR: uDY[7,13] - uDZ[7,9] = 0
 uSR: udy[7,13] - udz[7,9] = 0
 uSR: uDY[7,2] - uDZ[7,10] = 0
 uSR: udy[7,2] - udz[7,10] = 0
 uSR: uDY[7,7] - uDZ[7,11] = 0
 uSR: udy[7,7] - udz[7,11] = 0
 uSR: uDY[7,12] - uDZ[7,12] = 0
 uSR: udy[7,12] - udz[7,12] = 0
 uSR: uDY[7,1] - uDZ[7,13] = 0
 uSR: udy[7,1] - udz[7,13] = 0
 uSR: uDY[7,6] - uDZ[7,14] = 0
 uSR: udy[7,6] - udz[7,14] = 0
 uSR: uDY[7,11] - uDZ[7,15] = 0
 uSR: udy[7,11] - udz[7,15] = 0
 bMC: - udw[0,0] + uTC[0,0] >= 0
 bMC: - udw[0,1] + uTC[0,0] >= 0
 bMC: - udw[0,2] + uTC[0,0] >= 0
 bMC: - udw[0,3] + uTC[0,0] >= 0
 bMC: - udw[0,0] - udw[0,1] - udw[0,2] - udw[0,3] + uTC[0,0] <= 0
 bMC: - udz[0,0] + uTC[0,0] <= 0
 bMC: - udz[0,1] + uTC[0,0] <= 0
 bMC: - udz[0,2] + uTC[0,0] <= 0
 bMC: - udz[0,3] + uTC[0,0] <= 0
 bMC: - udz[0,0] - udz[0,1] - udz[0,2] - udz[0,3] + 4 uTC[0,0] >= 0
 bMC: - uDW[0,0] + uAC[0,0] >= 0
 bMC: - uDW[0,1] + uAC[0,0] >= 0
 bMC: - uDW[0,2] + uAC[0,0] >= 0
 bMC: - uDW[0,3] + uAC[0,0] >= 0
 bMC: - uDW[0,0] - uDW[0,1] - uDW[0,2] - uDW[0,3] + uAC[0,0] <= 0
 bMC: uDW[0,0] + uDW[0,1] + uDW[0,2] + uDW[0,3] + uDZ[0,0] + uDZ[0,1]
   + uDZ[0,2] + uDZ[0,3] - 5 uAC[0,0] >= 0
 bMC: uDW[0,0] + uDW[0,1] + uDW[0,2] + uDW[0,3] + uDZ[0,0] + uDZ[0,1]
   + uDZ[0,2] + uDZ[0,3] - 8 uAC[0,0] <= 0
 bMC: - udw[0,4] + uTC[0,1] >= 0
 bMC: - udw[0,5] + uTC[0,1] >= 0
 bMC: - udw[0,6] + uTC[0,1] >= 0
 bMC: - udw[0,7] + uTC[0,1] >= 0
 bMC: - udw[0,4] - udw[0,5] - udw[0,6] - udw[0,7] + uTC[0,1] <= 0
 bMC: - udz[0,4] + uTC[0,1] <= 0
 bMC: - udz[0,5] + uTC[0,1] <= 0
 bMC: - udz[0,6] + uTC[0,1] <= 0
 bMC: - udz[0,7] + uTC[0,1] <= 0
 bMC: - udz[0,4] - udz[0,5] - udz[0,6] - udz[0,7] + 4 uTC[0,1] >= 0
 bMC: - uDW[0,4] + uAC[0,1] >= 0
 bMC: - uDW[0,5] + uAC[0,1] >= 0
 bMC: - uDW[0,6] + uAC[0,1] >= 0
 bMC: - uDW[0,7] + uAC[0,1] >= 0
 bMC: - uDW[0,4] - uDW[0,5] - uDW[0,6] - uDW[0,7] + uAC[0,1] <= 0
 bMC: uDW[0,4] + uDW[0,5] + uDW[0,6] + uDW[0,7] + uDZ[0,4] + uDZ[0,5]
   + uDZ[0,6] + uDZ[0,7] - 5 uAC[0,1] >= 0
 bMC: uDW[0,4] + uDW[0,5] + uDW[0,6] + uDW[0,7] + uDZ[0,4] + uDZ[0,5]
   + uDZ[0,6] + uDZ[0,7] - 8 uAC[0,1] <= 0
 bMC: - udw[0,8] + uTC[0,2] >= 0
 bMC: - udw[0,9] + uTC[0,2] >= 0
 bMC: - udw[0,10] + uTC[0,2] >= 0
 bMC: - udw[0,11] + uTC[0,2] >= 0
 bMC: - udw[0,8] - udw[0,9] - udw[0,10] - udw[0,11] + uTC[0,2] <= 0
 bMC: - udz[0,8] + uTC[0,2] <= 0
 bMC: - udz[0,9] + uTC[0,2] <= 0
 bMC: - udz[0,10] + uTC[0,2] <= 0
 bMC: - udz[0,11] + uTC[0,2] <= 0
 bMC: - udz[0,8] - udz[0,9] - udz[0,10] - udz[0,11] + 4 uTC[0,2] >= 0
 bMC: - uDW[0,8] + uAC[0,2] >= 0
 bMC: - uDW[0,9] + uAC[0,2] >= 0
 bMC: - uDW[0,10] + uAC[0,2] >= 0
 bMC: - uDW[0,11] + uAC[0,2] >= 0
 bMC: - uDW[0,8] - uDW[0,9] - uDW[0,10] - uDW[0,11] + uAC[0,2] <= 0
 bMC: uDW[0,8] + uDW[0,9] + uDW[0,10] + uDW[0,11] + uDZ[0,8] + uDZ[0,9]
   + uDZ[0,10] + uDZ[0,11] - 5 uAC[0,2] >= 0
 bMC: uDW[0,8] + uDW[0,9] + uDW[0,10] + uDW[0,11] + uDZ[0,8] + uDZ[0,9]
   + uDZ[0,10] + uDZ[0,11] - 8 uAC[0,2] <= 0
 bMC: - udw[0,12] + uTC[0,3] >= 0
 bMC: - udw[0,13] + uTC[0,3] >= 0
 bMC: - udw[0,14] + uTC[0,3] >= 0
 bMC: - udw[0,15] + uTC[0,3] >= 0
 bMC: - udw[0,12] - udw[0,13] - udw[0,14] - udw[0,15] + uTC[0,3] <= 0
 bMC: - udz[0,12] + uTC[0,3] <= 0
 bMC: - udz[0,13] + uTC[0,3] <= 0
 bMC: - udz[0,14] + uTC[0,3] <= 0
 bMC: - udz[0,15] + uTC[0,3] <= 0
 bMC: - udz[0,12] - udz[0,13] - udz[0,14] - udz[0,15] + 4 uTC[0,3] >= 0
 bMC: - uDW[0,12] + uAC[0,3] >= 0
 bMC: - uDW[0,13] + uAC[0,3] >= 0
 bMC: - uDW[0,14] + uAC[0,3] >= 0
 bMC: - uDW[0,15] + uAC[0,3] >= 0
 bMC: - uDW[0,12] - uDW[0,13] - uDW[0,14] - uDW[0,15] + uAC[0,3] <= 0
 bMC: uDW[0,12] + uDW[0,13] + uDW[0,14] + uDW[0,15] + uDZ[0,12] + uDZ[0,13]
   + uDZ[0,14] + uDZ[0,15] - 5 uAC[0,3] >= 0
 bMC: uDW[0,12] + uDW[0,13] + uDW[0,14] + uDW[0,15] + uDZ[0,12] + uDZ[0,13]
   + uDZ[0,14] + uDZ[0,15] - 8 uAC[0,3] <= 0
 uMC: - udz[1,0] + uTC[1,0] >= 0
 uMC: - udz[1,1] + uTC[1,0] >= 0
 uMC: - udz[1,2] + uTC[1,0] >= 0
 uMC: - udz[1,3] + uTC[1,0] >= 0
 uMC: - udz[1,0] - udz[1,1] - udz[1,2] - udz[1,3] + uTC[1,0] <= 0
 uMC: - udw[1,0] + uTC[1,0] <= 0
 uMC: - udw[1,1] + uTC[1,0] <= 0
 uMC: - udw[1,2] + uTC[1,0] <= 0
 uMC: - udw[1,3] + uTC[1,0] <= 0
 uMC: - udw[1,0] - udw[1,1] - udw[1,2] - udw[1,3] + 4 uTC[1,0] >= 0
 uMC: - uDZ[1,0] + uAC[1,0] >= 0
 uMC: - uDZ[1,1] + uAC[1,0] >= 0
 uMC: - uDZ[1,2] + uAC[1,0] >= 0
 uMC: - uDZ[1,3] + uAC[1,0] >= 0
 uMC: - uDZ[1,0] - uDZ[1,1] - uDZ[1,2] - uDZ[1,3] + uAC[1,0] <= 0
 uMC: uDW[1,0] + uDW[1,1] + uDW[1,2] + uDW[1,3] + uDZ[1,0] + uDZ[1,1]
   + uDZ[1,2] + uDZ[1,3] - 5 uAC[1,0] >= 0
 uMC: uDW[1,0] + uDW[1,1] + uDW[1,2] + uDW[1,3] + uDZ[1,0] + uDZ[1,1]
   + uDZ[1,2] + uDZ[1,3] - 8 uAC[1,0] <= 0
 uMC: - udz[1,4] + uTC[1,1] >= 0
 uMC: - udz[1,5] + uTC[1,1] >= 0
 uMC: - udz[1,6] + uTC[1,1] >= 0
 uMC: - udz[1,7] + uTC[1,1] >= 0
 uMC: - udz[1,4] - udz[1,5] - udz[1,6] - udz[1,7] + uTC[1,1] <= 0
 uMC: - udw[1,4] + uTC[1,1] <= 0
 uMC: - udw[1,5] + uTC[1,1] <= 0
 uMC: - udw[1,6] + uTC[1,1] <= 0
 uMC: - udw[1,7] + uTC[1,1] <= 0
 uMC: - udw[1,4] - udw[1,5] - udw[1,6] - udw[1,7] + 4 uTC[1,1] >= 0
 uMC: - uDZ[1,4] + uAC[1,1] >= 0
 uMC: - uDZ[1,5] + uAC[1,1] >= 0
 uMC: - uDZ[1,6] + uAC[1,1] >= 0
 uMC: - uDZ[1,7] + uAC[1,1] >= 0
 uMC: - uDZ[1,4] - uDZ[1,5] - uDZ[1,6] - uDZ[1,7] + uAC[1,1] <= 0
 uMC: uDW[1,4] + uDW[1,5] + uDW[1,6] + uDW[1,7] + uDZ[1,4] + uDZ[1,5]
   + uDZ[1,6] + uDZ[1,7] - 5 uAC[1,1] >= 0
 uMC: uDW[1,4] + uDW[1,5] + uDW[1,6] + uDW[1,7] + uDZ[1,4] + uDZ[1,5]
   + uDZ[1,6] + uDZ[1,7] - 8 uAC[1,1] <= 0
 uMC: - udz[1,8] + uTC[1,2] >= 0
 uMC: - udz[1,9] + uTC[1,2] >= 0
 uMC: - udz[1,10] + uTC[1,2] >= 0
 uMC: - udz[1,11] + uTC[1,2] >= 0
 uMC: - udz[1,8] - udz[1,9] - udz[1,10] - udz[1,11] + uTC[1,2] <= 0
 uMC: - udw[1,8] + uTC[1,2] <= 0
 uMC: - udw[1,9] + uTC[1,2] <= 0
 uMC: - udw[1,10] + uTC[1,2] <= 0
 uMC: - udw[1,11] + uTC[1,2] <= 0
 uMC: - udw[1,8] - udw[1,9] - udw[1,10] - udw[1,11] + 4 uTC[1,2] >= 0
 uMC: - uDZ[1,8] + uAC[1,2] >= 0
 uMC: - uDZ[1,9] + uAC[1,2] >= 0
 uMC: - uDZ[1,10] + uAC[1,2] >= 0
 uMC: - uDZ[1,11] + uAC[1,2] >= 0
 uMC: - uDZ[1,8] - uDZ[1,9] - uDZ[1,10] - uDZ[1,11] + uAC[1,2] <= 0
 uMC: uDW[1,8] + uDW[1,9] + uDW[1,10] + uDW[1,11] + uDZ[1,8] + uDZ[1,9]
   + uDZ[1,10] + uDZ[1,11] - 5 uAC[1,2] >= 0
 uMC: uDW[1,8] + uDW[1,9] + uDW[1,10] + uDW[1,11] + uDZ[1,8] + uDZ[1,9]
   + uDZ[1,10] + uDZ[1,11] - 8 uAC[1,2] <= 0
 uMC: - udz[1,12] + uTC[1,3] >= 0
 uMC: - udz[1,13] + uTC[1,3] >= 0
 uMC: - udz[1,14] + uTC[1,3] >= 0
 uMC: - udz[1,15] + uTC[1,3] >= 0
 uMC: - udz[1,12] - udz[1,13] - udz[1,14] - udz[1,15] + uTC[1,3] <= 0
 uMC: - udw[1,12] + uTC[1,3] <= 0
 uMC: - udw[1,13] + uTC[1,3] <= 0
 uMC: - udw[1,14] + uTC[1,3] <= 0
 uMC: - udw[1,15] + uTC[1,3] <= 0
 uMC: - udw[1,12] - udw[1,13] - udw[1,14] - udw[1,15] + 4 uTC[1,3] >= 0
 uMC: - uDZ[1,12] + uAC[1,3] >= 0
 uMC: - uDZ[1,13] + uAC[1,3] >= 0
 uMC: - uDZ[1,14] + uAC[1,3] >= 0
 uMC: - uDZ[1,15] + uAC[1,3] >= 0
 uMC: - uDZ[1,12] - uDZ[1,13] - uDZ[1,14] - uDZ[1,15] + uAC[1,3] <= 0
 uMC: uDW[1,12] + uDW[1,13] + uDW[1,14] + uDW[1,15] + uDZ[1,12] + uDZ[1,13]
   + uDZ[1,14] + uDZ[1,15] - 5 uAC[1,3] >= 0
 uMC: uDW[1,12] + uDW[1,13] + uDW[1,14] + uDW[1,15] + uDZ[1,12] + uDZ[1,13]
   + uDZ[1,14] + uDZ[1,15] - 8 uAC[1,3] <= 0
 uMC: - udz[2,0] + uTC[2,0] >= 0
 uMC: - udz[2,1] + uTC[2,0] >= 0
 uMC: - udz[2,2] + uTC[2,0] >= 0
 uMC: - udz[2,3] + uTC[2,0] >= 0
 uMC: - udz[2,0] - udz[2,1] - udz[2,2] - udz[2,3] + uTC[2,0] <= 0
 uMC: - udw[2,0] + uTC[2,0] <= 0
 uMC: - udw[2,1] + uTC[2,0] <= 0
 uMC: - udw[2,2] + uTC[2,0] <= 0
 uMC: - udw[2,3] + uTC[2,0] <= 0
 uMC: - udw[2,0] - udw[2,1] - udw[2,2] - udw[2,3] + 4 uTC[2,0] >= 0
 uMC: - uDZ[2,0] + uAC[2,0] >= 0
 uMC: - uDZ[2,1] + uAC[2,0] >= 0
 uMC: - uDZ[2,2] + uAC[2,0] >= 0
 uMC: - uDZ[2,3] + uAC[2,0] >= 0
 uMC: - uDZ[2,0] - uDZ[2,1] - uDZ[2,2] - uDZ[2,3] + uAC[2,0] <= 0
 uMC: uDW[2,0] + uDW[2,1] + uDW[2,2] + uDW[2,3] + uDZ[2,0] + uDZ[2,1]
   + uDZ[2,2] + uDZ[2,3] - 5 uAC[2,0] >= 0
 uMC: uDW[2,0] + uDW[2,1] + uDW[2,2] + uDW[2,3] + uDZ[2,0] + uDZ[2,1]
   + uDZ[2,2] + uDZ[2,3] - 8 uAC[2,0] <= 0
 uMC: - udz[2,4] + uTC[2,1] >= 0
 uMC: - udz[2,5] + uTC[2,1] >= 0
 uMC: - udz[2,6] + uTC[2,1] >= 0
 uMC: - udz[2,7] + uTC[2,1] >= 0
 uMC: - udz[2,4] - udz[2,5] - udz[2,6] - udz[2,7] + uTC[2,1] <= 0
 uMC: - udw[2,4] + uTC[2,1] <= 0
 uMC: - udw[2,5] + uTC[2,1] <= 0
 uMC: - udw[2,6] + uTC[2,1] <= 0
 uMC: - udw[2,7] + uTC[2,1] <= 0
 uMC: - udw[2,4] - udw[2,5] - udw[2,6] - udw[2,7] + 4 uTC[2,1] >= 0
 uMC: - uDZ[2,4] + uAC[2,1] >= 0
 uMC: - uDZ[2,5] + uAC[2,1] >= 0
 uMC: - uDZ[2,6] + uAC[2,1] >= 0
 uMC: - uDZ[2,7] + uAC[2,1] >= 0
 uMC: - uDZ[2,4] - uDZ[2,5] - uDZ[2,6] - uDZ[2,7] + uAC[2,1] <= 0
 uMC: uDW[2,4] + uDW[2,5] + uDW[2,6] + uDW[2,7] + uDZ[2,4] + uDZ[2,5]
   + uDZ[2,6] + uDZ[2,7] - 5 uAC[2,1] >= 0
 uMC: uDW[2,4] + uDW[2,5] + uDW[2,6] + uDW[2,7] + uDZ[2,4] + uDZ[2,5]
   + uDZ[2,6] + uDZ[2,7] - 8 uAC[2,1] <= 0
 uMC: - udz[2,8] + uTC[2,2] >= 0
 uMC: - udz[2,9] + uTC[2,2] >= 0
 uMC: - udz[2,10] + uTC[2,2] >= 0
 uMC: - udz[2,11] + uTC[2,2] >= 0
 uMC: - udz[2,8] - udz[2,9] - udz[2,10] - udz[2,11] + uTC[2,2] <= 0
 uMC: - udw[2,8] + uTC[2,2] <= 0
 uMC: - udw[2,9] + uTC[2,2] <= 0
 uMC: - udw[2,10] + uTC[2,2] <= 0
 uMC: - udw[2,11] + uTC[2,2] <= 0
 uMC: - udw[2,8] - udw[2,9] - udw[2,10] - udw[2,11] + 4 uTC[2,2] >= 0
 uMC: - uDZ[2,8] + uAC[2,2] >= 0
 uMC: - uDZ[2,9] + uAC[2,2] >= 0
 uMC: - uDZ[2,10] + uAC[2,2] >= 0
 uMC: - uDZ[2,11] + uAC[2,2] >= 0
 uMC: - uDZ[2,8] - uDZ[2,9] - uDZ[2,10] - uDZ[2,11] + uAC[2,2] <= 0
 uMC: uDW[2,8] + uDW[2,9] + uDW[2,10] + uDW[2,11] + uDZ[2,8] + uDZ[2,9]
   + uDZ[2,10] + uDZ[2,11] - 5 uAC[2,2] >= 0
 uMC: uDW[2,8] + uDW[2,9] + uDW[2,10] + uDW[2,11] + uDZ[2,8] + uDZ[2,9]
   + uDZ[2,10] + uDZ[2,11] - 8 uAC[2,2] <= 0
 uMC: - udz[2,12] + uTC[2,3] >= 0
 uMC: - udz[2,13] + uTC[2,3] >= 0
 uMC: - udz[2,14] + uTC[2,3] >= 0
 uMC: - udz[2,15] + uTC[2,3] >= 0
 uMC: - udz[2,12] - udz[2,13] - udz[2,14] - udz[2,15] + uTC[2,3] <= 0
 uMC: - udw[2,12] + uTC[2,3] <= 0
 uMC: - udw[2,13] + uTC[2,3] <= 0
 uMC: - udw[2,14] + uTC[2,3] <= 0
 uMC: - udw[2,15] + uTC[2,3] <= 0
 uMC: - udw[2,12] - udw[2,13] - udw[2,14] - udw[2,15] + 4 uTC[2,3] >= 0
 uMC: - uDZ[2,12] + uAC[2,3] >= 0
 uMC: - uDZ[2,13] + uAC[2,3] >= 0
 uMC: - uDZ[2,14] + uAC[2,3] >= 0
 uMC: - uDZ[2,15] + uAC[2,3] >= 0
 uMC: - uDZ[2,12] - uDZ[2,13] - uDZ[2,14] - uDZ[2,15] + uAC[2,3] <= 0
 uMC: uDW[2,12] + uDW[2,13] + uDW[2,14] + uDW[2,15] + uDZ[2,12] + uDZ[2,13]
   + uDZ[2,14] + uDZ[2,15] - 5 uAC[2,3] >= 0
 uMC: uDW[2,12] + uDW[2,13] + uDW[2,14] + uDW[2,15] + uDZ[2,12] + uDZ[2,13]
   + uDZ[2,14] + uDZ[2,15] - 8 uAC[2,3] <= 0
 uMC: - udz[3,0] + uTC[3,0] >= 0
 uMC: - udz[3,1] + uTC[3,0] >= 0
 uMC: - udz[3,2] + uTC[3,0] >= 0
 uMC: - udz[3,3] + uTC[3,0] >= 0
 uMC: - udz[3,0] - udz[3,1] - udz[3,2] - udz[3,3] + uTC[3,0] <= 0
 uMC: - udw[3,0] + uTC[3,0] <= 0
 uMC: - udw[3,1] + uTC[3,0] <= 0
 uMC: - udw[3,2] + uTC[3,0] <= 0
 uMC: - udw[3,3] + uTC[3,0] <= 0
 uMC: - udw[3,0] - udw[3,1] - udw[3,2] - udw[3,3] + 4 uTC[3,0] >= 0
 uMC: - uDZ[3,0] + uAC[3,0] >= 0
 uMC: - uDZ[3,1] + uAC[3,0] >= 0
 uMC: - uDZ[3,2] + uAC[3,0] >= 0
 uMC: - uDZ[3,3] + uAC[3,0] >= 0
 uMC: - uDZ[3,0] - uDZ[3,1] - uDZ[3,2] - uDZ[3,3] + uAC[3,0] <= 0
 uMC: uDW[3,0] + uDW[3,1] + uDW[3,2] + uDW[3,3] + uDZ[3,0] + uDZ[3,1]
   + uDZ[3,2] + uDZ[3,3] - 5 uAC[3,0] >= 0
 uMC: uDW[3,0] + uDW[3,1] + uDW[3,2] + uDW[3,3] + uDZ[3,0] + uDZ[3,1]
   + uDZ[3,2] + uDZ[3,3] - 8 uAC[3,0] <= 0
 uMC: - udz[3,4] + uTC[3,1] >= 0
 uMC: - udz[3,5] + uTC[3,1] >= 0
 uMC: - udz[3,6] + uTC[3,1] >= 0
 uMC: - udz[3,7] + uTC[3,1] >= 0
 uMC: - udz[3,4] - udz[3,5] - udz[3,6] - udz[3,7] + uTC[3,1] <= 0
 uMC: - udw[3,4] + uTC[3,1] <= 0
 uMC: - udw[3,5] + uTC[3,1] <= 0
 uMC: - udw[3,6] + uTC[3,1] <= 0
 uMC: - udw[3,7] + uTC[3,1] <= 0
 uMC: - udw[3,4] - udw[3,5] - udw[3,6] - udw[3,7] + 4 uTC[3,1] >= 0
 uMC: - uDZ[3,4] + uAC[3,1] >= 0
 uMC: - uDZ[3,5] + uAC[3,1] >= 0
 uMC: - uDZ[3,6] + uAC[3,1] >= 0
 uMC: - uDZ[3,7] + uAC[3,1] >= 0
 uMC: - uDZ[3,4] - uDZ[3,5] - uDZ[3,6] - uDZ[3,7] + uAC[3,1] <= 0
 uMC: uDW[3,4] + uDW[3,5] + uDW[3,6] + uDW[3,7] + uDZ[3,4] + uDZ[3,5]
   + uDZ[3,6] + uDZ[3,7] - 5 uAC[3,1] >= 0
 uMC: uDW[3,4] + uDW[3,5] + uDW[3,6] + uDW[3,7] + uDZ[3,4] + uDZ[3,5]
   + uDZ[3,6] + uDZ[3,7] - 8 uAC[3,1] <= 0
 uMC: - udz[3,8] + uTC[3,2] >= 0
 uMC: - udz[3,9] + uTC[3,2] >= 0
 uMC: - udz[3,10] + uTC[3,2] >= 0
 uMC: - udz[3,11] + uTC[3,2] >= 0
 uMC: - udz[3,8] - udz[3,9] - udz[3,10] - udz[3,11] + uTC[3,2] <= 0
 uMC: - udw[3,8] + uTC[3,2] <= 0
 uMC: - udw[3,9] + uTC[3,2] <= 0
 uMC: - udw[3,10] + uTC[3,2] <= 0
 uMC: - udw[3,11] + uTC[3,2] <= 0
 uMC: - udw[3,8] - udw[3,9] - udw[3,10] - udw[3,11] + 4 uTC[3,2] >= 0
 uMC: - uDZ[3,8] + uAC[3,2] >= 0
 uMC: - uDZ[3,9] + uAC[3,2] >= 0
 uMC: - uDZ[3,10] + uAC[3,2] >= 0
 uMC: - uDZ[3,11] + uAC[3,2] >= 0
 uMC: - uDZ[3,8] - uDZ[3,9] - uDZ[3,10] - uDZ[3,11] + uAC[3,2] <= 0
 uMC: uDW[3,8] + uDW[3,9] + uDW[3,10] + uDW[3,11] + uDZ[3,8] + uDZ[3,9]
   + uDZ[3,10] + uDZ[3,11] - 5 uAC[3,2] >= 0
 uMC: uDW[3,8] + uDW[3,9] + uDW[3,10] + uDW[3,11] + uDZ[3,8] + uDZ[3,9]
   + uDZ[3,10] + uDZ[3,11] - 8 uAC[3,2] <= 0
 uMC: - udz[3,12] + uTC[3,3] >= 0
 uMC: - udz[3,13] + uTC[3,3] >= 0
 uMC: - udz[3,14] + uTC[3,3] >= 0
 uMC: - udz[3,15] + uTC[3,3] >= 0
 uMC: - udz[3,12] - udz[3,13] - udz[3,14] - udz[3,15] + uTC[3,3] <= 0
 uMC: - udw[3,12] + uTC[3,3] <= 0
 uMC: - udw[3,13] + uTC[3,3] <= 0
 uMC: - udw[3,14] + uTC[3,3] <= 0
 uMC: - udw[3,15] + uTC[3,3] <= 0
 uMC: - udw[3,12] - udw[3,13] - udw[3,14] - udw[3,15] + 4 uTC[3,3] >= 0
 uMC: - uDZ[3,12] + uAC[3,3] >= 0
 uMC: - uDZ[3,13] + uAC[3,3] >= 0
 uMC: - uDZ[3,14] + uAC[3,3] >= 0
 uMC: - uDZ[3,15] + uAC[3,3] >= 0
 uMC: - uDZ[3,12] - uDZ[3,13] - uDZ[3,14] - uDZ[3,15] + uAC[3,3] <= 0
 uMC: uDW[3,12] + uDW[3,13] + uDW[3,14] + uDW[3,15] + uDZ[3,12] + uDZ[3,13]
   + uDZ[3,14] + uDZ[3,15] - 5 uAC[3,3] >= 0
 uMC: uDW[3,12] + uDW[3,13] + uDW[3,14] + uDW[3,15] + uDZ[3,12] + uDZ[3,13]
   + uDZ[3,14] + uDZ[3,15] - 8 uAC[3,3] <= 0
 uMC: - udz[4,0] + uTC[4,0] >= 0
 uMC: - udz[4,1] + uTC[4,0] >= 0
 uMC: - udz[4,2] + uTC[4,0] >= 0
 uMC: - udz[4,3] + uTC[4,0] >= 0
 uMC: - udz[4,0] - udz[4,1] - udz[4,2] - udz[4,3] + uTC[4,0] <= 0
 uMC: - udw[4,0] + uTC[4,0] <= 0
 uMC: - udw[4,1] + uTC[4,0] <= 0
 uMC: - udw[4,2] + uTC[4,0] <= 0
 uMC: - udw[4,3] + uTC[4,0] <= 0
 uMC: - udw[4,0] - udw[4,1] - udw[4,2] - udw[4,3] + 4 uTC[4,0] >= 0
 uMC: - uDZ[4,0] + uAC[4,0] >= 0
 uMC: - uDZ[4,1] + uAC[4,0] >= 0
 uMC: - uDZ[4,2] + uAC[4,0] >= 0
 uMC: - uDZ[4,3] + uAC[4,0] >= 0
 uMC: - uDZ[4,0] - uDZ[4,1] - uDZ[4,2] - uDZ[4,3] + uAC[4,0] <= 0
 uMC: uDW[4,0] + uDW[4,1] + uDW[4,2] + uDW[4,3] + uDZ[4,0] + uDZ[4,1]
   + uDZ[4,2] + uDZ[4,3] - 5 uAC[4,0] >= 0
 uMC: uDW[4,0] + uDW[4,1] + uDW[4,2] + uDW[4,3] + uDZ[4,0] + uDZ[4,1]
   + uDZ[4,2] + uDZ[4,3] - 8 uAC[4,0] <= 0
 uMC: - udz[4,4] + uTC[4,1] >= 0
 uMC: - udz[4,5] + uTC[4,1] >= 0
 uMC: - udz[4,6] + uTC[4,1] >= 0
 uMC: - udz[4,7] + uTC[4,1] >= 0
 uMC: - udz[4,4] - udz[4,5] - udz[4,6] - udz[4,7] + uTC[4,1] <= 0
 uMC: - udw[4,4] + uTC[4,1] <= 0
 uMC: - udw[4,5] + uTC[4,1] <= 0
 uMC: - udw[4,6] + uTC[4,1] <= 0
 uMC: - udw[4,7] + uTC[4,1] <= 0
 uMC: - udw[4,4] - udw[4,5] - udw[4,6] - udw[4,7] + 4 uTC[4,1] >= 0
 uMC: - uDZ[4,4] + uAC[4,1] >= 0
 uMC: - uDZ[4,5] + uAC[4,1] >= 0
 uMC: - uDZ[4,6] + uAC[4,1] >= 0
 uMC: - uDZ[4,7] + uAC[4,1] >= 0
 uMC: - uDZ[4,4] - uDZ[4,5] - uDZ[4,6] - uDZ[4,7] + uAC[4,1] <= 0
 uMC: uDW[4,4] + uDW[4,5] + uDW[4,6] + uDW[4,7] + uDZ[4,4] + uDZ[4,5]
   + uDZ[4,6] + uDZ[4,7] - 5 uAC[4,1] >= 0
 uMC: uDW[4,4] + uDW[4,5] + uDW[4,6] + uDW[4,7] + uDZ[4,4] + uDZ[4,5]
   + uDZ[4,6] + uDZ[4,7] - 8 uAC[4,1] <= 0
 uMC: - udz[4,8] + uTC[4,2] >= 0
 uMC: - udz[4,9] + uTC[4,2] >= 0
 uMC: - udz[4,10] + uTC[4,2] >= 0
 uMC: - udz[4,11] + uTC[4,2] >= 0
 uMC: - udz[4,8] - udz[4,9] - udz[4,10] - udz[4,11] + uTC[4,2] <= 0
 uMC: - udw[4,8] + uTC[4,2] <= 0
 uMC: - udw[4,9] + uTC[4,2] <= 0
 uMC: - udw[4,10] + uTC[4,2] <= 0
 uMC: - udw[4,11] + uTC[4,2] <= 0
 uMC: - udw[4,8] - udw[4,9] - udw[4,10] - udw[4,11] + 4 uTC[4,2] >= 0
 uMC: - uDZ[4,8] + uAC[4,2] >= 0
 uMC: - uDZ[4,9] + uAC[4,2] >= 0
 uMC: - uDZ[4,10] + uAC[4,2] >= 0
 uMC: - uDZ[4,11] + uAC[4,2] >= 0
 uMC: - uDZ[4,8] - uDZ[4,9] - uDZ[4,10] - uDZ[4,11] + uAC[4,2] <= 0
 uMC: uDW[4,8] + uDW[4,9] + uDW[4,10] + uDW[4,11] + uDZ[4,8] + uDZ[4,9]
   + uDZ[4,10] + uDZ[4,11] - 5 uAC[4,2] >= 0
 uMC: uDW[4,8] + uDW[4,9] + uDW[4,10] + uDW[4,11] + uDZ[4,8] + uDZ[4,9]
   + uDZ[4,10] + uDZ[4,11] - 8 uAC[4,2] <= 0
 uMC: - udz[4,12] + uTC[4,3] >= 0
 uMC: - udz[4,13] + uTC[4,3] >= 0
 uMC: - udz[4,14] + uTC[4,3] >= 0
 uMC: - udz[4,15] + uTC[4,3] >= 0
 uMC: - udz[4,12] - udz[4,13] - udz[4,14] - udz[4,15] + uTC[4,3] <= 0
 uMC: - udw[4,12] + uTC[4,3] <= 0
 uMC: - udw[4,13] + uTC[4,3] <= 0
 uMC: - udw[4,14] + uTC[4,3] <= 0
 uMC: - udw[4,15] + uTC[4,3] <= 0
 uMC: - udw[4,12] - udw[4,13] - udw[4,14] - udw[4,15] + 4 uTC[4,3] >= 0
 uMC: - uDZ[4,12] + uAC[4,3] >= 0
 uMC: - uDZ[4,13] + uAC[4,3] >= 0
 uMC: - uDZ[4,14] + uAC[4,3] >= 0
 uMC: - uDZ[4,15] + uAC[4,3] >= 0
 uMC: - uDZ[4,12] - uDZ[4,13] - uDZ[4,14] - uDZ[4,15] + uAC[4,3] <= 0
 uMC: uDW[4,12] + uDW[4,13] + uDW[4,14] + uDW[4,15] + uDZ[4,12] + uDZ[4,13]
   + uDZ[4,14] + uDZ[4,15] - 5 uAC[4,3] >= 0
 uMC: uDW[4,12] + uDW[4,13] + uDW[4,14] + uDW[4,15] + uDZ[4,12] + uDZ[4,13]
   + uDZ[4,14] + uDZ[4,15] - 8 uAC[4,3] <= 0
 uMC: - udz[5,0] + uTC[5,0] >= 0
 uMC: - udz[5,1] + uTC[5,0] >= 0
 uMC: - udz[5,2] + uTC[5,0] >= 0
 uMC: - udz[5,3] + uTC[5,0] >= 0
 uMC: - udz[5,0] - udz[5,1] - udz[5,2] - udz[5,3] + uTC[5,0] <= 0
 uMC: - udw[5,0] + uTC[5,0] <= 0
 uMC: - udw[5,1] + uTC[5,0] <= 0
 uMC: - udw[5,2] + uTC[5,0] <= 0
 uMC: - udw[5,3] + uTC[5,0] <= 0
 uMC: - udw[5,0] - udw[5,1] - udw[5,2] - udw[5,3] + 4 uTC[5,0] >= 0
 uMC: - uDZ[5,0] + uAC[5,0] >= 0
 uMC: - uDZ[5,1] + uAC[5,0] >= 0
 uMC: - uDZ[5,2] + uAC[5,0] >= 0
 uMC: - uDZ[5,3] + uAC[5,0] >= 0
 uMC: - uDZ[5,0] - uDZ[5,1] - uDZ[5,2] - uDZ[5,3] + uAC[5,0] <= 0
 uMC: uDW[5,0] + uDW[5,1] + uDW[5,2] + uDW[5,3] + uDZ[5,0] + uDZ[5,1]
   + uDZ[5,2] + uDZ[5,3] - 5 uAC[5,0] >= 0
 uMC: uDW[5,0] + uDW[5,1] + uDW[5,2] + uDW[5,3] + uDZ[5,0] + uDZ[5,1]
   + uDZ[5,2] + uDZ[5,3] - 8 uAC[5,0] <= 0
 uMC: - udz[5,4] + uTC[5,1] >= 0
 uMC: - udz[5,5] + uTC[5,1] >= 0
 uMC: - udz[5,6] + uTC[5,1] >= 0
 uMC: - udz[5,7] + uTC[5,1] >= 0
 uMC: - udz[5,4] - udz[5,5] - udz[5,6] - udz[5,7] + uTC[5,1] <= 0
 uMC: - udw[5,4] + uTC[5,1] <= 0
 uMC: - udw[5,5] + uTC[5,1] <= 0
 uMC: - udw[5,6] + uTC[5,1] <= 0
 uMC: - udw[5,7] + uTC[5,1] <= 0
 uMC: - udw[5,4] - udw[5,5] - udw[5,6] - udw[5,7] + 4 uTC[5,1] >= 0
 uMC: - uDZ[5,4] + uAC[5,1] >= 0
 uMC: - uDZ[5,5] + uAC[5,1] >= 0
 uMC: - uDZ[5,6] + uAC[5,1] >= 0
 uMC: - uDZ[5,7] + uAC[5,1] >= 0
 uMC: - uDZ[5,4] - uDZ[5,5] - uDZ[5,6] - uDZ[5,7] + uAC[5,1] <= 0
 uMC: uDW[5,4] + uDW[5,5] + uDW[5,6] + uDW[5,7] + uDZ[5,4] + uDZ[5,5]
   + uDZ[5,6] + uDZ[5,7] - 5 uAC[5,1] >= 0
 uMC: uDW[5,4] + uDW[5,5] + uDW[5,6] + uDW[5,7] + uDZ[5,4] + uDZ[5,5]
   + uDZ[5,6] + uDZ[5,7] - 8 uAC[5,1] <= 0
 uMC: - udz[5,8] + uTC[5,2] >= 0
 uMC: - udz[5,9] + uTC[5,2] >= 0
 uMC: - udz[5,10] + uTC[5,2] >= 0
 uMC: - udz[5,11] + uTC[5,2] >= 0
 uMC: - udz[5,8] - udz[5,9] - udz[5,10] - udz[5,11] + uTC[5,2] <= 0
 uMC: - udw[5,8] + uTC[5,2] <= 0
 uMC: - udw[5,9] + uTC[5,2] <= 0
 uMC: - udw[5,10] + uTC[5,2] <= 0
 uMC: - udw[5,11] + uTC[5,2] <= 0
 uMC: - udw[5,8] - udw[5,9] - udw[5,10] - udw[5,11] + 4 uTC[5,2] >= 0
 uMC: - uDZ[5,8] + uAC[5,2] >= 0
 uMC: - uDZ[5,9] + uAC[5,2] >= 0
 uMC: - uDZ[5,10] + uAC[5,2] >= 0
 uMC: - uDZ[5,11] + uAC[5,2] >= 0
 uMC: - uDZ[5,8] - uDZ[5,9] - uDZ[5,10] - uDZ[5,11] + uAC[5,2] <= 0
 uMC: uDW[5,8] + uDW[5,9] + uDW[5,10] + uDW[5,11] + uDZ[5,8] + uDZ[5,9]
   + uDZ[5,10] + uDZ[5,11] - 5 uAC[5,2] >= 0
 uMC: uDW[5,8] + uDW[5,9] + uDW[5,10] + uDW[5,11] + uDZ[5,8] + uDZ[5,9]
   + uDZ[5,10] + uDZ[5,11] - 8 uAC[5,2] <= 0
 uMC: - udz[5,12] + uTC[5,3] >= 0
 uMC: - udz[5,13] + uTC[5,3] >= 0
 uMC: - udz[5,14] + uTC[5,3] >= 0
 uMC: - udz[5,15] + uTC[5,3] >= 0
 uMC: - udz[5,12] - udz[5,13] - udz[5,14] - udz[5,15] + uTC[5,3] <= 0
 uMC: - udw[5,12] + uTC[5,3] <= 0
 uMC: - udw[5,13] + uTC[5,3] <= 0
 uMC: - udw[5,14] + uTC[5,3] <= 0
 uMC: - udw[5,15] + uTC[5,3] <= 0
 uMC: - udw[5,12] - udw[5,13] - udw[5,14] - udw[5,15] + 4 uTC[5,3] >= 0
 uMC: - uDZ[5,12] + uAC[5,3] >= 0
 uMC: - uDZ[5,13] + uAC[5,3] >= 0
 uMC: - uDZ[5,14] + uAC[5,3] >= 0
 uMC: - uDZ[5,15] + uAC[5,3] >= 0
 uMC: - uDZ[5,12] - uDZ[5,13] - uDZ[5,14] - uDZ[5,15] + uAC[5,3] <= 0
 uMC: uDW[5,12] + uDW[5,13] + uDW[5,14] + uDW[5,15] + uDZ[5,12] + uDZ[5,13]
   + uDZ[5,14] + uDZ[5,15] - 5 uAC[5,3] >= 0
 uMC: uDW[5,12] + uDW[5,13] + uDW[5,14] + uDW[5,15] + uDZ[5,12] + uDZ[5,13]
   + uDZ[5,14] + uDZ[5,15] - 8 uAC[5,3] <= 0
 uMC: - udz[6,0] + uTC[6,0] >= 0
 uMC: - udz[6,1] + uTC[6,0] >= 0
 uMC: - udz[6,2] + uTC[6,0] >= 0
 uMC: - udz[6,3] + uTC[6,0] >= 0
 uMC: - udz[6,0] - udz[6,1] - udz[6,2] - udz[6,3] + uTC[6,0] <= 0
 uMC: - udw[6,0] + uTC[6,0] <= 0
 uMC: - udw[6,1] + uTC[6,0] <= 0
 uMC: - udw[6,2] + uTC[6,0] <= 0
 uMC: - udw[6,3] + uTC[6,0] <= 0
 uMC: - udw[6,0] - udw[6,1] - udw[6,2] - udw[6,3] + 4 uTC[6,0] >= 0
 uMC: - uDZ[6,0] + uAC[6,0] >= 0
 uMC: - uDZ[6,1] + uAC[6,0] >= 0
 uMC: - uDZ[6,2] + uAC[6,0] >= 0
 uMC: - uDZ[6,3] + uAC[6,0] >= 0
 uMC: - uDZ[6,0] - uDZ[6,1] - uDZ[6,2] - uDZ[6,3] + uAC[6,0] <= 0
 uMC: uDW[6,0] + uDW[6,1] + uDW[6,2] + uDW[6,3] + uDZ[6,0] + uDZ[6,1]
   + uDZ[6,2] + uDZ[6,3] - 5 uAC[6,0] >= 0
 uMC: uDW[6,0] + uDW[6,1] + uDW[6,2] + uDW[6,3] + uDZ[6,0] + uDZ[6,1]
   + uDZ[6,2] + uDZ[6,3] - 8 uAC[6,0] <= 0
 uMC: - udz[6,4] + uTC[6,1] >= 0
 uMC: - udz[6,5] + uTC[6,1] >= 0
 uMC: - udz[6,6] + uTC[6,1] >= 0
 uMC: - udz[6,7] + uTC[6,1] >= 0
 uMC: - udz[6,4] - udz[6,5] - udz[6,6] - udz[6,7] + uTC[6,1] <= 0
 uMC: - udw[6,4] + uTC[6,1] <= 0
 uMC: - udw[6,5] + uTC[6,1] <= 0
 uMC: - udw[6,6] + uTC[6,1] <= 0
 uMC: - udw[6,7] + uTC[6,1] <= 0
 uMC: - udw[6,4] - udw[6,5] - udw[6,6] - udw[6,7] + 4 uTC[6,1] >= 0
 uMC: - uDZ[6,4] + uAC[6,1] >= 0
 uMC: - uDZ[6,5] + uAC[6,1] >= 0
 uMC: - uDZ[6,6] + uAC[6,1] >= 0
 uMC: - uDZ[6,7] + uAC[6,1] >= 0
 uMC: - uDZ[6,4] - uDZ[6,5] - uDZ[6,6] - uDZ[6,7] + uAC[6,1] <= 0
 uMC: uDW[6,4] + uDW[6,5] + uDW[6,6] + uDW[6,7] + uDZ[6,4] + uDZ[6,5]
   + uDZ[6,6] + uDZ[6,7] - 5 uAC[6,1] >= 0
 uMC: uDW[6,4] + uDW[6,5] + uDW[6,6] + uDW[6,7] + uDZ[6,4] + uDZ[6,5]
   + uDZ[6,6] + uDZ[6,7] - 8 uAC[6,1] <= 0
 uMC: - udz[6,8] + uTC[6,2] >= 0
 uMC: - udz[6,9] + uTC[6,2] >= 0
 uMC: - udz[6,10] + uTC[6,2] >= 0
 uMC: - udz[6,11] + uTC[6,2] >= 0
 uMC: - udz[6,8] - udz[6,9] - udz[6,10] - udz[6,11] + uTC[6,2] <= 0
 uMC: - udw[6,8] + uTC[6,2] <= 0
 uMC: - udw[6,9] + uTC[6,2] <= 0
 uMC: - udw[6,10] + uTC[6,2] <= 0
 uMC: - udw[6,11] + uTC[6,2] <= 0
 uMC: - udw[6,8] - udw[6,9] - udw[6,10] - udw[6,11] + 4 uTC[6,2] >= 0
 uMC: - uDZ[6,8] + uAC[6,2] >= 0
 uMC: - uDZ[6,9] + uAC[6,2] >= 0
 uMC: - uDZ[6,10] + uAC[6,2] >= 0
 uMC: - uDZ[6,11] + uAC[6,2] >= 0
 uMC: - uDZ[6,8] - uDZ[6,9] - uDZ[6,10] - uDZ[6,11] + uAC[6,2] <= 0
 uMC: uDW[6,8] + uDW[6,9] + uDW[6,10] + uDW[6,11] + uDZ[6,8] + uDZ[6,9]
   + uDZ[6,10] + uDZ[6,11] - 5 uAC[6,2] >= 0
 uMC: uDW[6,8] + uDW[6,9] + uDW[6,10] + uDW[6,11] + uDZ[6,8] + uDZ[6,9]
   + uDZ[6,10] + uDZ[6,11] - 8 uAC[6,2] <= 0
 uMC: - udz[6,12] + uTC[6,3] >= 0
 uMC: - udz[6,13] + uTC[6,3] >= 0
 uMC: - udz[6,14] + uTC[6,3] >= 0
 uMC: - udz[6,15] + uTC[6,3] >= 0
 uMC: - udz[6,12] - udz[6,13] - udz[6,14] - udz[6,15] + uTC[6,3] <= 0
 uMC: - udw[6,12] + uTC[6,3] <= 0
 uMC: - udw[6,13] + uTC[6,3] <= 0
 uMC: - udw[6,14] + uTC[6,3] <= 0
 uMC: - udw[6,15] + uTC[6,3] <= 0
 uMC: - udw[6,12] - udw[6,13] - udw[6,14] - udw[6,15] + 4 uTC[6,3] >= 0
 uMC: - uDZ[6,12] + uAC[6,3] >= 0
 uMC: - uDZ[6,13] + uAC[6,3] >= 0
 uMC: - uDZ[6,14] + uAC[6,3] >= 0
 uMC: - uDZ[6,15] + uAC[6,3] >= 0
 uMC: - uDZ[6,12] - uDZ[6,13] - uDZ[6,14] - uDZ[6,15] + uAC[6,3] <= 0
 uMC: uDW[6,12] + uDW[6,13] + uDW[6,14] + uDW[6,15] + uDZ[6,12] + uDZ[6,13]
   + uDZ[6,14] + uDZ[6,15] - 5 uAC[6,3] >= 0
 uMC: uDW[6,12] + uDW[6,13] + uDW[6,14] + uDW[6,15] + uDZ[6,12] + uDZ[6,13]
   + uDZ[6,14] + uDZ[6,15] - 8 uAC[6,3] <= 0
 uMC: - udz[7,0] + uTC[7,0] >= 0
 uMC: - udz[7,1] + uTC[7,0] >= 0
 uMC: - udz[7,2] + uTC[7,0] >= 0
 uMC: - udz[7,3] + uTC[7,0] >= 0
 uMC: - udz[7,0] - udz[7,1] - udz[7,2] - udz[7,3] + uTC[7,0] <= 0
 uMC: - udw[7,0] + uTC[7,0] <= 0
 uMC: - udw[7,1] + uTC[7,0] <= 0
 uMC: - udw[7,2] + uTC[7,0] <= 0
 uMC: - udw[7,3] + uTC[7,0] <= 0
 uMC: - udw[7,0] - udw[7,1] - udw[7,2] - udw[7,3] + 4 uTC[7,0] >= 0
 uMC: - uDZ[7,0] + uAC[7,0] >= 0
 uMC: - uDZ[7,1] + uAC[7,0] >= 0
 uMC: - uDZ[7,2] + uAC[7,0] >= 0
 uMC: - uDZ[7,3] + uAC[7,0] >= 0
 uMC: - uDZ[7,0] - uDZ[7,1] - uDZ[7,2] - uDZ[7,3] + uAC[7,0] <= 0
 uMC: uDW[7,0] + uDW[7,1] + uDW[7,2] + uDW[7,3] + uDZ[7,0] + uDZ[7,1]
   + uDZ[7,2] + uDZ[7,3] - 5 uAC[7,0] >= 0
 uMC: uDW[7,0] + uDW[7,1] + uDW[7,2] + uDW[7,3] + uDZ[7,0] + uDZ[7,1]
   + uDZ[7,2] + uDZ[7,3] - 8 uAC[7,0] <= 0
 uMC: - udz[7,4] + uTC[7,1] >= 0
 uMC: - udz[7,5] + uTC[7,1] >= 0
 uMC: - udz[7,6] + uTC[7,1] >= 0
 uMC: - udz[7,7] + uTC[7,1] >= 0
 uMC: - udz[7,4] - udz[7,5] - udz[7,6] - udz[7,7] + uTC[7,1] <= 0
 uMC: - udw[7,4] + uTC[7,1] <= 0
 uMC: - udw[7,5] + uTC[7,1] <= 0
 uMC: - udw[7,6] + uTC[7,1] <= 0
 uMC: - udw[7,7] + uTC[7,1] <= 0
 uMC: - udw[7,4] - udw[7,5] - udw[7,6] - udw[7,7] + 4 uTC[7,1] >= 0
 uMC: - uDZ[7,4] + uAC[7,1] >= 0
 uMC: - uDZ[7,5] + uAC[7,1] >= 0
 uMC: - uDZ[7,6] + uAC[7,1] >= 0
 uMC: - uDZ[7,7] + uAC[7,1] >= 0
 uMC: - uDZ[7,4] - uDZ[7,5] - uDZ[7,6] - uDZ[7,7] + uAC[7,1] <= 0
 uMC: uDW[7,4] + uDW[7,5] + uDW[7,6] + uDW[7,7] + uDZ[7,4] + uDZ[7,5]
   + uDZ[7,6] + uDZ[7,7] - 5 uAC[7,1] >= 0
 uMC: uDW[7,4] + uDW[7,5] + uDW[7,6] + uDW[7,7] + uDZ[7,4] + uDZ[7,5]
   + uDZ[7,6] + uDZ[7,7] - 8 uAC[7,1] <= 0
 uMC: - udz[7,8] + uTC[7,2] >= 0
 uMC: - udz[7,9] + uTC[7,2] >= 0
 uMC: - udz[7,10] + uTC[7,2] >= 0
 uMC: - udz[7,11] + uTC[7,2] >= 0
 uMC: - udz[7,8] - udz[7,9] - udz[7,10] - udz[7,11] + uTC[7,2] <= 0
 uMC: - udw[7,8] + uTC[7,2] <= 0
 uMC: - udw[7,9] + uTC[7,2] <= 0
 uMC: - udw[7,10] + uTC[7,2] <= 0
 uMC: - udw[7,11] + uTC[7,2] <= 0
 uMC: - udw[7,8] - udw[7,9] - udw[7,10] - udw[7,11] + 4 uTC[7,2] >= 0
 uMC: - uDZ[7,8] + uAC[7,2] >= 0
 uMC: - uDZ[7,9] + uAC[7,2] >= 0
 uMC: - uDZ[7,10] + uAC[7,2] >= 0
 uMC: - uDZ[7,11] + uAC[7,2] >= 0
 uMC: - uDZ[7,8] - uDZ[7,9] - uDZ[7,10] - uDZ[7,11] + uAC[7,2] <= 0
 uMC: uDW[7,8] + uDW[7,9] + uDW[7,10] + uDW[7,11] + uDZ[7,8] + uDZ[7,9]
   + uDZ[7,10] + uDZ[7,11] - 5 uAC[7,2] >= 0
 uMC: uDW[7,8] + uDW[7,9] + uDW[7,10] + uDW[7,11] + uDZ[7,8] + uDZ[7,9]
   + uDZ[7,10] + uDZ[7,11] - 8 uAC[7,2] <= 0
 uMC: - udz[7,12] + uTC[7,3] >= 0
 uMC: - udz[7,13] + uTC[7,3] >= 0
 uMC: - udz[7,14] + uTC[7,3] >= 0
 uMC: - udz[7,15] + uTC[7,3] >= 0
 uMC: - udz[7,12] - udz[7,13] - udz[7,14] - udz[7,15] + uTC[7,3] <= 0
 uMC: - udw[7,12] + uTC[7,3] <= 0
 uMC: - udw[7,13] + uTC[7,3] <= 0
 uMC: - udw[7,14] + uTC[7,3] <= 0
 uMC: - udw[7,15] + uTC[7,3] <= 0
 uMC: - udw[7,12] - udw[7,13] - udw[7,14] - udw[7,15] + 4 uTC[7,3] >= 0
 uMC: - uDZ[7,12] + uAC[7,3] >= 0
 uMC: - uDZ[7,13] + uAC[7,3] >= 0
 uMC: - uDZ[7,14] + uAC[7,3] >= 0
 uMC: - uDZ[7,15] + uAC[7,3] >= 0
 uMC: - uDZ[7,12] - uDZ[7,13] - uDZ[7,14] - uDZ[7,15] + uAC[7,3] <= 0
 uMC: uDW[7,12] + uDW[7,13] + uDW[7,14] + uDW[7,15] + uDZ[7,12] + uDZ[7,13]
   + uDZ[7,14] + uDZ[7,15] - 5 uAC[7,3] >= 0
 uMC: uDW[7,12] + uDW[7,13] + uDW[7,14] + uDW[7,15] + uDZ[7,12] + uDZ[7,13]
   + uDZ[7,14] + uDZ[7,15] - 8 uAC[7,3] <= 0
 R2224: uLANE[0] + uLANE[1] + uLANE[2] + uLANE[3] + uLANE[4] + uLANE[5]
   + uLANE[6] + uLANE[7] + uLANE[8] + uLANE[9] + uLANE[10] + uLANE[11]
   + uLANE[12] + uLANE[13] + uLANE[14] + uLANE[15] >= 1
 uMLANE: - ustk[0,0] + uLANE[0] >= 0
 uMLANE: - ustk[1,7] + uLANE[0] >= 0
 uMLANE: - ustk[2,14] + uLANE[0] >= 0
 uMLANE: - ustk[3,9] + uLANE[0] >= 0
 uMLANE: - ustk[4,8] + uLANE[0] >= 0
 uMLANE: - ustk[5,15] + uLANE[0] >= 0
 uMLANE: - ustk[6,6] + uLANE[0] >= 0
 uMLANE: - ustk[7,1] + uLANE[0] >= 0
 uMLANE: - ustk[8,0] + uLANE[0] >= 0
 uMLANE: - ustk[0,1] + uLANE[1] >= 0
 uMLANE: - ustk[1,0] + uLANE[1] >= 0
 uMLANE: - ustk[2,7] + uLANE[1] >= 0
 uMLANE: - ustk[3,14] + uLANE[1] >= 0
 uMLANE: - ustk[4,9] + uLANE[1] >= 0
 uMLANE: - ustk[5,8] + uLANE[1] >= 0
 uMLANE: - ustk[6,15] + uLANE[1] >= 0
 uMLANE: - ustk[7,6] + uLANE[1] >= 0
 uMLANE: - ustk[8,1] + uLANE[1] >= 0
 uMLANE: - ustk[0,2] + uLANE[2] >= 0
 uMLANE: - ustk[1,13] + uLANE[2] >= 0
 uMLANE: - ustk[2,12] + uLANE[2] >= 0
 uMLANE: - ustk[3,3] + uLANE[2] >= 0
 uMLANE: - ustk[4,10] + uLANE[2] >= 0
 uMLANE: - ustk[5,5] + uLANE[2] >= 0
 uMLANE: - ustk[6,4] + uLANE[2] >= 0
 uMLANE: - ustk[7,11] + uLANE[2] >= 0
 uMLANE: - ustk[8,2] + uLANE[2] >= 0
 uMLANE: - ustk[0,3] + uLANE[3] >= 0
 uMLANE: - ustk[1,10] + uLANE[3] >= 0
 uMLANE: - ustk[2,5] + uLANE[3] >= 0
 uMLANE: - ustk[3,4] + uLANE[3] >= 0
 uMLANE: - ustk[4,11] + uLANE[3] >= 0
 uMLANE: - ustk[5,2] + uLANE[3] >= 0
 uMLANE: - ustk[6,13] + uLANE[3] >= 0
 uMLANE: - ustk[7,12] + uLANE[3] >= 0
 uMLANE: - ustk[8,3] + uLANE[3] >= 0
 uMLANE: - ustk[0,4] + uLANE[4] >= 0
 uMLANE: - ustk[1,11] + uLANE[4] >= 0
 uMLANE: - ustk[2,2] + uLANE[4] >= 0
 uMLANE: - ustk[3,13] + uLANE[4] >= 0
 uMLANE: - ustk[4,12] + uLANE[4] >= 0
 uMLANE: - ustk[5,3] + uLANE[4] >= 0
 uMLANE: - ustk[6,10] + uLANE[4] >= 0
 uMLANE: - ustk[7,5] + uLANE[4] >= 0
 uMLANE: - ustk[8,4] + uLANE[4] >= 0
 uMLANE: - ustk[0,5] + uLANE[5] >= 0
 uMLANE: - ustk[1,4] + uLANE[5] >= 0
 uMLANE: - ustk[2,11] + uLANE[5] >= 0
 uMLANE: - ustk[3,2] + uLANE[5] >= 0
 uMLANE: - ustk[4,13] + uLANE[5] >= 0
 uMLANE: - ustk[5,12] + uLANE[5] >= 0
 uMLANE: - ustk[6,3] + uLANE[5] >= 0
 uMLANE: - ustk[7,10] + uLANE[5] >= 0
 uMLANE: - ustk[8,5] + uLANE[5] >= 0
 uMLANE: - ustk[0,6] + uLANE[6] >= 0
 uMLANE: - ustk[1,1] + uLANE[6] >= 0
 uMLANE: - ustk[2,0] + uLANE[6] >= 0
 uMLANE: - ustk[3,7] + uLANE[6] >= 0
 uMLANE: - ustk[4,14] + uLANE[6] >= 0
 uMLANE: - ustk[5,9] + uLANE[6] >= 0
 uMLANE: - ustk[6,8] + uLANE[6] >= 0
 uMLANE: - ustk[7,15] + uLANE[6] >= 0
 uMLANE: - ustk[8,6] + uLANE[6] >= 0
 uMLANE: - ustk[0,7] + uLANE[7] >= 0
 uMLANE: - ustk[1,14] + uLANE[7] >= 0
 uMLANE: - ustk[2,9] + uLANE[7] >= 0
 uMLANE: - ustk[3,8] + uLANE[7] >= 0
 uMLANE: - ustk[4,15] + uLANE[7] >= 0
 uMLANE: - ustk[5,6] + uLANE[7] >= 0
 uMLANE: - ustk[6,1] + uLANE[7] >= 0
 uMLANE: - ustk[7,0] + uLANE[7] >= 0
 uMLANE: - ustk[8,7] + uLANE[7] >= 0
 uMLANE: - ustk[0,8] + uLANE[8] >= 0
 uMLANE: - ustk[1,15] + uLANE[8] >= 0
 uMLANE: - ustk[2,6] + uLANE[8] >= 0
 uMLANE: - ustk[3,1] + uLANE[8] >= 0
 uMLANE: - ustk[4,0] + uLANE[8] >= 0
 uMLANE: - ustk[5,7] + uLANE[8] >= 0
 uMLANE: - ustk[6,14] + uLANE[8] >= 0
 uMLANE: - ustk[7,9] + uLANE[8] >= 0
 uMLANE: - ustk[8,8] + uLANE[8] >= 0
 uMLANE: - ustk[0,9] + uLANE[9] >= 0
 uMLANE: - ustk[1,8] + uLANE[9] >= 0
 uMLANE: - ustk[2,15] + uLANE[9] >= 0
 uMLANE: - ustk[3,6] + uLANE[9] >= 0
 uMLANE: - ustk[4,1] + uLANE[9] >= 0
 uMLANE: - ustk[5,0] + uLANE[9] >= 0
 uMLANE: - ustk[6,7] + uLANE[9] >= 0
 uMLANE: - ustk[7,14] + uLANE[9] >= 0
 uMLANE: - ustk[8,9] + uLANE[9] >= 0
 uMLANE: - ustk[0,10] + uLANE[10] >= 0
 uMLANE: - ustk[1,5] + uLANE[10] >= 0
 uMLANE: - ustk[2,4] + uLANE[10] >= 0
 uMLANE: - ustk[3,11] + uLANE[10] >= 0
 uMLANE: - ustk[4,2] + uLANE[10] >= 0
 uMLANE: - ustk[5,13] + uLANE[10] >= 0
 uMLANE: - ustk[6,12] + uLANE[10] >= 0
 uMLANE: - ustk[7,3] + uLANE[10] >= 0
 uMLANE: - ustk[8,10] + uLANE[10] >= 0
 uMLANE: - ustk[0,11] + uLANE[11] >= 0
 uMLANE: - ustk[1,2] + uLANE[11] >= 0
 uMLANE: - ustk[2,13] + uLANE[11] >= 0
 uMLANE: - ustk[3,12] + uLANE[11] >= 0
 uMLANE: - ustk[4,3] + uLANE[11] >= 0
 uMLANE: - ustk[5,10] + uLANE[11] >= 0
 uMLANE: - ustk[6,5] + uLANE[11] >= 0
 uMLANE: - ustk[7,4] + uLANE[11] >= 0
 uMLANE: - ustk[8,11] + uLANE[11] >= 0
 uMLANE: - ustk[0,12] + uLANE[12] >= 0
 uMLANE: - ustk[1,3] + uLANE[12] >= 0
 uMLANE: - ustk[2,10] + uLANE[12] >= 0
 uMLANE: - ustk[3,5] + uLANE[12] >= 0
 uMLANE: - ustk[4,4] + uLANE[12] >= 0
 uMLANE: - ustk[5,11] + uLANE[12] >= 0
 uMLANE: - ustk[6,2] + uLANE[12] >= 0
 uMLANE: - ustk[7,13] + uLANE[12] >= 0
 uMLANE: - ustk[8,12] + uLANE[12] >= 0
 uMLANE: - ustk[0,13] + uLANE[13] >= 0
 uMLANE: - ustk[1,12] + uLANE[13] >= 0
 uMLANE: - ustk[2,3] + uLANE[13] >= 0
 uMLANE: - ustk[3,10] + uLANE[13] >= 0
 uMLANE: - ustk[4,5] + uLANE[13] >= 0
 uMLANE: - ustk[5,4] + uLANE[13] >= 0
 uMLANE: - ustk[6,11] + uLANE[13] >= 0
 uMLANE: - ustk[7,2] + uLANE[13] >= 0
 uMLANE: - ustk[8,13] + uLANE[13] >= 0
 uMLANE: - ustk[0,14] + uLANE[14] >= 0
 uMLANE: - ustk[1,9] + uLANE[14] >= 0
 uMLANE: - ustk[2,8] + uLANE[14] >= 0
 uMLANE: - ustk[3,15] + uLANE[14] >= 0
 uMLANE: - ustk[4,6] + uLANE[14] >= 0
 uMLANE: - ustk[5,1] + uLANE[14] >= 0
 uMLANE: - ustk[6,0] + uLANE[14] >= 0
 uMLANE: - ustk[7,7] + uLANE[14] >= 0
 uMLANE: - ustk[8,14] + uLANE[14] >= 0
 uMLANE: - ustk[0,15] + uLANE[15] >= 0
 uMLANE: - ustk[1,6] + uLANE[15] >= 0
 uMLANE: - ustk[2,1] + uLANE[15] >= 0
 uMLANE: - ustk[3,0] + uLANE[15] >= 0
 uMLANE: - ustk[4,7] + uLANE[15] >= 0
 uMLANE: - ustk[5,14] + uLANE[15] >= 0
 uMLANE: - ustk[6,9] + uLANE[15] >= 0
 uMLANE: - ustk[7,8] + uLANE[15] >= 0
 uMLANE: - ustk[8,15] + uLANE[15] >= 0
 uType1: ustk[0,0] + ustk[1,7] + ustk[2,14] + ustk[3,9] + ustk[4,8]
   + ustk[5,15] + ustk[6,6] + ustk[7,1] + ustk[8,0] - 9 uLANE[0]
   + uType1[0] = 0
 uType1: uType1[0] <= 1
 uType1: ustk[0,1] + ustk[1,0] + ustk[2,7] + ustk[3,14] + ustk[4,9]
   + ustk[5,8] + ustk[6,15] + ustk[7,6] + ustk[8,1] - 9 uLANE[1]
   + uType1[1] = 0
 uType1: uType1[1] <= 1
 uType1: ustk[0,2] + ustk[1,13] + ustk[2,12] + ustk[3,3] + ustk[4,10]
   + ustk[5,5] + ustk[6,4] + ustk[7,11] + ustk[8,2] - 9 uLANE[2]
   + uType1[2] = 0
 uType1: uType1[2] <= 1
 uType1: ustk[0,3] + ustk[1,10] + ustk[2,5] + ustk[3,4] + ustk[4,11]
   + ustk[5,2] + ustk[6,13] + ustk[7,12] + ustk[8,3] - 9 uLANE[3]
   + uType1[3] = 0
 uType1: uType1[3] <= 1
 uType1: ustk[0,4] + ustk[1,11] + ustk[2,2] + ustk[3,13] + ustk[4,12]
   + ustk[5,3] + ustk[6,10] + ustk[7,5] + ustk[8,4] - 9 uLANE[4]
   + uType1[4] = 0
 uType1: uType1[4] <= 1
 uType1: ustk[0,5] + ustk[1,4] + ustk[2,11] + ustk[3,2] + ustk[4,13]
   + ustk[5,12] + ustk[6,3] + ustk[7,10] + ustk[8,5] - 9 uLANE[5]
   + uType1[5] = 0
 uType1: uType1[5] <= 1
 uType1: ustk[0,6] + ustk[1,1] + ustk[2,0] + ustk[3,7] + ustk[4,14]
   + ustk[5,9] + ustk[6,8] + ustk[7,15] + ustk[8,6] - 9 uLANE[6]
   + uType1[6] = 0
 uType1: uType1[6] <= 1
 uType1: ustk[0,7] + ustk[1,14] + ustk[2,9] + ustk[3,8] + ustk[4,15]
   + ustk[5,6] + ustk[6,1] + ustk[7,0] + ustk[8,7] - 9 uLANE[7] + uType1[7]
   = 0
 uType1: uType1[7] <= 1
 uType1: ustk[0,8] + ustk[1,15] + ustk[2,6] + ustk[3,1] + ustk[4,0]
   + ustk[5,7] + ustk[6,14] + ustk[7,9] + ustk[8,8] - 9 uLANE[8]
   + uType1[8] = 0
 uType1: uType1[8] <= 1
 uType1: ustk[0,9] + ustk[1,8] + ustk[2,15] + ustk[3,6] + ustk[4,1]
   + ustk[5,0] + ustk[6,7] + ustk[7,14] + ustk[8,9] - 9 uLANE[9]
   + uType1[9] = 0
 uType1: uType1[9] <= 1
 uType1: ustk[0,10] + ustk[1,5] + ustk[2,4] + ustk[3,11] + ustk[4,2]
   + ustk[5,13] + ustk[6,12] + ustk[7,3] + ustk[8,10] - 9 uLANE[10]
   + uType1[10] = 0
 uType1: uType1[10] <= 1
 uType1: ustk[0,11] + ustk[1,2] + ustk[2,13] + ustk[3,12] + ustk[4,3]
   + ustk[5,10] + ustk[6,5] + ustk[7,4] + ustk[8,11] - 9 uLANE[11]
   + uType1[11] = 0
 uType1: uType1[11] <= 1
 uType1: ustk[0,12] + ustk[1,3] + ustk[2,10] + ustk[3,5] + ustk[4,4]
   + ustk[5,11] + ustk[6,2] + ustk[7,13] + ustk[8,12] - 9 uLANE[12]
   + uType1[12] = 0
 uType1: uType1[12] <= 1
 uType1: ustk[0,13] + ustk[1,12] + ustk[2,3] + ustk[3,10] + ustk[4,5]
   + ustk[5,4] + ustk[6,11] + ustk[7,2] + ustk[8,13] - 9 uLANE[13]
   + uType1[13] = 0
 uType1: uType1[13] <= 1
 uType1: ustk[0,14] + ustk[1,9] + ustk[2,8] + ustk[3,15] + ustk[4,6]
   + ustk[5,1] + ustk[6,0] + ustk[7,7] + ustk[8,14] - 9 uLANE[14]
   + uType1[14] = 0
 uType1: uType1[14] <= 1
 uType1: ustk[0,15] + ustk[1,6] + ustk[2,1] + ustk[3,0] + ustk[4,7]
   + ustk[5,14] + ustk[6,9] + ustk[7,8] + ustk[8,15] - 9 uLANE[15]
   + uType1[15] = 0
 uType1: uType1[15] <= 1
 uType2A: - uDZ[0,0] - uDZ[0,1] - uDZ[0,2] - uDZ[0,3] + uT2CanA[0,0] = 0
 uType2B: uDW[0,0] + uDW[0,1] + uDW[0,2] + uDW[0,3] - 4 uAC[0,0]
   + uT2CanB[0,0] = 0
 uType2A: - uDZ[0,4] - uDZ[0,5] - uDZ[0,6] - uDZ[0,7] + uT2CanA[0,1] = 0
 uType2B: uDW[0,4] + uDW[0,5] + uDW[0,6] + uDW[0,7] - 4 uAC[0,1]
   + uT2CanB[0,1] = 0
 uType2A: - uDZ[0,8] - uDZ[0,9] - uDZ[0,10] - uDZ[0,11] + uT2CanA[0,2] = 0
 uType2B: uDW[0,8] + uDW[0,9] + uDW[0,10] + uDW[0,11] - 4 uAC[0,2]
   + uT2CanB[0,2] = 0
 uType2A: - uDZ[0,12] - uDZ[0,13] - uDZ[0,14] - uDZ[0,15] + uT2CanA[0,3]
   = 0
 uType2B: uDW[0,12] + uDW[0,13] + uDW[0,14] + uDW[0,15] - 4 uAC[0,3]
   + uT2CanB[0,3] = 0
 uType2A: - uDZ[1,0] - uDZ[1,1] - uDZ[1,2] - uDZ[1,3] + uT2CanA[1,0] = 0
 uType2B: uDW[1,0] + uDW[1,1] + uDW[1,2] + uDW[1,3] - 4 uAC[1,0]
   + uT2CanB[1,0] = 0
 uType2A: - uDZ[1,4] - uDZ[1,5] - uDZ[1,6] - uDZ[1,7] + uT2CanA[1,1] = 0
 uType2B: uDW[1,4] + uDW[1,5] + uDW[1,6] + uDW[1,7] - 4 uAC[1,1]
   + uT2CanB[1,1] = 0
 uType2A: - uDZ[1,8] - uDZ[1,9] - uDZ[1,10] - uDZ[1,11] + uT2CanA[1,2] = 0
 uType2B: uDW[1,8] + uDW[1,9] + uDW[1,10] + uDW[1,11] - 4 uAC[1,2]
   + uT2CanB[1,2] = 0
 uType2A: - uDZ[1,12] - uDZ[1,13] - uDZ[1,14] - uDZ[1,15] + uT2CanA[1,3]
   = 0
 uType2B: uDW[1,12] + uDW[1,13] + uDW[1,14] + uDW[1,15] - 4 uAC[1,3]
   + uT2CanB[1,3] = 0
 uType2A: - uDZ[2,0] - uDZ[2,1] - uDZ[2,2] - uDZ[2,3] + uT2CanA[2,0] = 0
 uType2B: uDW[2,0] + uDW[2,1] + uDW[2,2] + uDW[2,3] - 4 uAC[2,0]
   + uT2CanB[2,0] = 0
 uType2A: - uDZ[2,4] - uDZ[2,5] - uDZ[2,6] - uDZ[2,7] + uT2CanA[2,1] = 0
 uType2B: uDW[2,4] + uDW[2,5] + uDW[2,6] + uDW[2,7] - 4 uAC[2,1]
   + uT2CanB[2,1] = 0
 uType2A: - uDZ[2,8] - uDZ[2,9] - uDZ[2,10] - uDZ[2,11] + uT2CanA[2,2] = 0
 uType2B: uDW[2,8] + uDW[2,9] + uDW[2,10] + uDW[2,11] - 4 uAC[2,2]
   + uT2CanB[2,2] = 0
 uType2A: - uDZ[2,12] - uDZ[2,13] - uDZ[2,14] - uDZ[2,15] + uT2CanA[2,3]
   = 0
 uType2B: uDW[2,12] + uDW[2,13] + uDW[2,14] + uDW[2,15] - 4 uAC[2,3]
   + uT2CanB[2,3] = 0
 uType2A: - uDZ[3,0] - uDZ[3,1] - uDZ[3,2] - uDZ[3,3] + uT2CanA[3,0] = 0
 uType2B: uDW[3,0] + uDW[3,1] + uDW[3,2] + uDW[3,3] - 4 uAC[3,0]
   + uT2CanB[3,0] = 0
 uType2A: - uDZ[3,4] - uDZ[3,5] - uDZ[3,6] - uDZ[3,7] + uT2CanA[3,1] = 0
 uType2B: uDW[3,4] + uDW[3,5] + uDW[3,6] + uDW[3,7] - 4 uAC[3,1]
   + uT2CanB[3,1] = 0
 uType2A: - uDZ[3,8] - uDZ[3,9] - uDZ[3,10] - uDZ[3,11] + uT2CanA[3,2] = 0
 uType2B: uDW[3,8] + uDW[3,9] + uDW[3,10] + uDW[3,11] - 4 uAC[3,2]
   + uT2CanB[3,2] = 0
 uType2A: - uDZ[3,12] - uDZ[3,13] - uDZ[3,14] - uDZ[3,15] + uT2CanA[3,3]
   = 0
 uType2B: uDW[3,12] + uDW[3,13] + uDW[3,14] + uDW[3,15] - 4 uAC[3,3]
   + uT2CanB[3,3] = 0
 uType2A: - uDZ[4,0] - uDZ[4,1] - uDZ[4,2] - uDZ[4,3] + uT2CanA[4,0] = 0
 uType2B: uDW[4,0] + uDW[4,1] + uDW[4,2] + uDW[4,3] - 4 uAC[4,0]
   + uT2CanB[4,0] = 0
 uType2A: - uDZ[4,4] - uDZ[4,5] - uDZ[4,6] - uDZ[4,7] + uT2CanA[4,1] = 0
 uType2B: uDW[4,4] + uDW[4,5] + uDW[4,6] + uDW[4,7] - 4 uAC[4,1]
   + uT2CanB[4,1] = 0
 uType2A: - uDZ[4,8] - uDZ[4,9] - uDZ[4,10] - uDZ[4,11] + uT2CanA[4,2] = 0
 uType2B: uDW[4,8] + uDW[4,9] + uDW[4,10] + uDW[4,11] - 4 uAC[4,2]
   + uT2CanB[4,2] = 0
 uType2A: - uDZ[4,12] - uDZ[4,13] - uDZ[4,14] - uDZ[4,15] + uT2CanA[4,3]
   = 0
 uType2B: uDW[4,12] + uDW[4,13] + uDW[4,14] + uDW[4,15] - 4 uAC[4,3]
   + uT2CanB[4,3] = 0
 uType2A: - uDZ[5,0] - uDZ[5,1] - uDZ[5,2] - uDZ[5,3] + uT2CanA[5,0] = 0
 uType2B: uDW[5,0] + uDW[5,1] + uDW[5,2] + uDW[5,3] - 4 uAC[5,0]
   + uT2CanB[5,0] = 0
 uType2A: - uDZ[5,4] - uDZ[5,5] - uDZ[5,6] - uDZ[5,7] + uT2CanA[5,1] = 0
 uType2B: uDW[5,4] + uDW[5,5] + uDW[5,6] + uDW[5,7] - 4 uAC[5,1]
   + uT2CanB[5,1] = 0
 uType2A: - uDZ[5,8] - uDZ[5,9] - uDZ[5,10] - uDZ[5,11] + uT2CanA[5,2] = 0
 uType2B: uDW[5,8] + uDW[5,9] + uDW[5,10] + uDW[5,11] - 4 uAC[5,2]
   + uT2CanB[5,2] = 0
 uType2A: - uDZ[5,12] - uDZ[5,13] - uDZ[5,14] - uDZ[5,15] + uT2CanA[5,3]
   = 0
 uType2B: uDW[5,12] + uDW[5,13] + uDW[5,14] + uDW[5,15] - 4 uAC[5,3]
   + uT2CanB[5,3] = 0
 uType2A: - uDZ[6,0] - uDZ[6,1] - uDZ[6,2] - uDZ[6,3] + uT2CanA[6,0] = 0
 uType2B: uDW[6,0] + uDW[6,1] + uDW[6,2] + uDW[6,3] - 4 uAC[6,0]
   + uT2CanB[6,0] = 0
 uType2A: - uDZ[6,4] - uDZ[6,5] - uDZ[6,6] - uDZ[6,7] + uT2CanA[6,1] = 0
 uType2B: uDW[6,4] + uDW[6,5] + uDW[6,6] + uDW[6,7] - 4 uAC[6,1]
   + uT2CanB[6,1] = 0
 uType2A: - uDZ[6,8] - uDZ[6,9] - uDZ[6,10] - uDZ[6,11] + uT2CanA[6,2] = 0
 uType2B: uDW[6,8] + uDW[6,9] + uDW[6,10] + uDW[6,11] - 4 uAC[6,2]
   + uT2CanB[6,2] = 0
 uType2A: - uDZ[6,12] - uDZ[6,13] - uDZ[6,14] - uDZ[6,15] + uT2CanA[6,3]
   = 0
 uType2B: uDW[6,12] + uDW[6,13] + uDW[6,14] + uDW[6,15] - 4 uAC[6,3]
   + uT2CanB[6,3] = 0
 uType2A: - uDZ[7,0] - uDZ[7,1] - uDZ[7,2] - uDZ[7,3] + uT2CanA[7,0] = 0
 uType2B: uDW[7,0] + uDW[7,1] + uDW[7,2] + uDW[7,3] - 4 uAC[7,0]
   + uT2CanB[7,0] = 0
 uType2A: - uDZ[7,4] - uDZ[7,5] - uDZ[7,6] - uDZ[7,7] + uT2CanA[7,1] = 0
 uType2B: uDW[7,4] + uDW[7,5] + uDW[7,6] + uDW[7,7] - 4 uAC[7,1]
   + uT2CanB[7,1] = 0
 uType2A: - uDZ[7,8] - uDZ[7,9] - uDZ[7,10] - uDZ[7,11] + uT2CanA[7,2] = 0
 uType2B: uDW[7,8] + uDW[7,9] + uDW[7,10] + uDW[7,11] - 4 uAC[7,2]
   + uT2CanB[7,2] = 0
 uType2A: - uDZ[7,12] - uDZ[7,13] - uDZ[7,14] - uDZ[7,15] + uT2CanA[7,3]
   = 0
 uType2B: uDW[7,12] + uDW[7,13] + uDW[7,14] + uDW[7,15] - 4 uAC[7,3]
   + uT2CanB[7,3] = 0
 uType2C: - ucan[1,0] - ucan[1,1] - ucan[1,2] - ucan[1,3] + uT2CanC[1,0]
   = 0
 uType2C: - ucan[1,4] - ucan[1,5] - ucan[1,6] - ucan[1,7] + uT2CanC[1,1]
   = 0
 uType2C: - ucan[1,8] - ucan[1,9] - ucan[1,10] - ucan[1,11] + uT2CanC[1,2]
   = 0
 uType2C: - ucan[1,12] - ucan[1,13] - ucan[1,14] - ucan[1,15]
   + uT2CanC[1,3] = 0
 uType2C: - ucan[2,0] - ucan[2,1] - ucan[2,2] - ucan[2,3] + uT2CanC[2,0]
   = 0
 uType2C: - ucan[2,4] - ucan[2,5] - ucan[2,6] - ucan[2,7] + uT2CanC[2,1]
   = 0
 uType2C: - ucan[2,8] - ucan[2,9] - ucan[2,10] - ucan[2,11] + uT2CanC[2,2]
   = 0
 uType2C: - ucan[2,12] - ucan[2,13] - ucan[2,14] - ucan[2,15]
   + uT2CanC[2,3] = 0
 uType2C: - ucan[3,0] - ucan[3,1] - ucan[3,2] - ucan[3,3] + uT2CanC[3,0]
   = 0
 uType2C: - ucan[3,4] - ucan[3,5] - ucan[3,6] - ucan[3,7] + uT2CanC[3,1]
   = 0
 uType2C: - ucan[3,8] - ucan[3,9] - ucan[3,10] - ucan[3,11] + uT2CanC[3,2]
   = 0
 uType2C: - ucan[3,12] - ucan[3,13] - ucan[3,14] - ucan[3,15]
   + uT2CanC[3,3] = 0
 uType2C: - ucan[4,0] - ucan[4,1] - ucan[4,2] - ucan[4,3] + uT2CanC[4,0]
   = 0
 uType2C: - ucan[4,4] - ucan[4,5] - ucan[4,6] - ucan[4,7] + uT2CanC[4,1]
   = 0
 uType2C: - ucan[4,8] - ucan[4,9] - ucan[4,10] - ucan[4,11] + uT2CanC[4,2]
   = 0
 uType2C: - ucan[4,12] - ucan[4,13] - ucan[4,14] - ucan[4,15]
   + uT2CanC[4,3] = 0
 uType2C: - ucan[5,0] - ucan[5,1] - ucan[5,2] - ucan[5,3] + uT2CanC[5,0]
   = 0
 uType2C: - ucan[5,4] - ucan[5,5] - ucan[5,6] - ucan[5,7] + uT2CanC[5,1]
   = 0
 uType2C: - ucan[5,8] - ucan[5,9] - ucan[5,10] - ucan[5,11] + uT2CanC[5,2]
   = 0
 uType2C: - ucan[5,12] - ucan[5,13] - ucan[5,14] - ucan[5,15]
   + uT2CanC[5,3] = 0
 uType2C: - ucan[6,0] - ucan[6,1] - ucan[6,2] - ucan[6,3] + uT2CanC[6,0]
   = 0
 uType2C: - ucan[6,4] - ucan[6,5] - ucan[6,6] - ucan[6,7] + uT2CanC[6,1]
   = 0
 uType2C: - ucan[6,8] - ucan[6,9] - ucan[6,10] - ucan[6,11] + uT2CanC[6,2]
   = 0
 uType2C: - ucan[6,12] - ucan[6,13] - ucan[6,14] - ucan[6,15]
   + uT2CanC[6,3] = 0
 uType2C: - ucan[7,0] - ucan[7,1] - ucan[7,2] - ucan[7,3] + uT2CanC[7,0]
   = 0
 uType2C: - ucan[7,4] - ucan[7,5] - ucan[7,6] - ucan[7,7] + uT2CanC[7,1]
   = 0
 uType2C: - ucan[7,8] - ucan[7,9] - ucan[7,10] - ucan[7,11] + uT2CanC[7,2]
   = 0
 uType2C: - ucan[7,12] - ucan[7,13] - ucan[7,14] - ucan[7,15]
   + uT2CanC[7,3] = 0
 uType2: uT2CanA[0,0] - uT2CanB[0,0] - uT2CanC[1,0] + uType2[0,0] >= 0
 uType2: uType2[0,0] >= 0
 uType2: uT2CanA[0,1] - uT2CanB[0,1] - uT2CanC[1,1] + uType2[0,1] >= 0
 uType2: uType2[0,1] >= 0
 uType2: uT2CanA[0,2] - uT2CanB[0,2] - uT2CanC[1,2] + uType2[0,2] >= 0
 uType2: uType2[0,2] >= 0
 uType2: uT2CanA[0,3] - uT2CanB[0,3] - uT2CanC[1,3] + uType2[0,3] >= 0
 uType2: uType2[0,3] >= 0
 uType2: uT2CanA[1,0] - uT2CanB[1,0] - uT2CanC[2,0] + uType2[1,0] >= 0
 uType2: uType2[1,0] >= 0
 uType2: uT2CanA[1,1] - uT2CanB[1,1] - uT2CanC[2,1] + uType2[1,1] >= 0
 uType2: uType2[1,1] >= 0
 uType2: uT2CanA[1,2] - uT2CanB[1,2] - uT2CanC[2,2] + uType2[1,2] >= 0
 uType2: uType2[1,2] >= 0
 uType2: uT2CanA[1,3] - uT2CanB[1,3] - uT2CanC[2,3] + uType2[1,3] >= 0
 uType2: uType2[1,3] >= 0
 uType2: uT2CanA[2,0] - uT2CanB[2,0] - uT2CanC[3,0] + uType2[2,0] >= 0
 uType2: uType2[2,0] >= 0
 uType2: uT2CanA[2,1] - uT2CanB[2,1] - uT2CanC[3,1] + uType2[2,1] >= 0
 uType2: uType2[2,1] >= 0
 uType2: uT2CanA[2,2] - uT2CanB[2,2] - uT2CanC[3,2] + uType2[2,2] >= 0
 uType2: uType2[2,2] >= 0
 uType2: uT2CanA[2,3] - uT2CanB[2,3] - uT2CanC[3,3] + uType2[2,3] >= 0
 uType2: uType2[2,3] >= 0
 uType2: uT2CanA[3,0] - uT2CanB[3,0] - uT2CanC[4,0] + uType2[3,0] >= 0
 uType2: uType2[3,0] >= 0
 uType2: uT2CanA[3,1] - uT2CanB[3,1] - uT2CanC[4,1] + uType2[3,1] >= 0
 uType2: uType2[3,1] >= 0
 uType2: uT2CanA[3,2] - uT2CanB[3,2] - uT2CanC[4,2] + uType2[3,2] >= 0
 uType2: uType2[3,2] >= 0
 uType2: uT2CanA[3,3] - uT2CanB[3,3] - uT2CanC[4,3] + uType2[3,3] >= 0
 uType2: uType2[3,3] >= 0
 uType2: uT2CanA[4,0] - uT2CanB[4,0] - uT2CanC[5,0] + uType2[4,0] >= 0
 uType2: uType2[4,0] >= 0
 uType2: uT2CanA[4,1] - uT2CanB[4,1] - uT2CanC[5,1] + uType2[4,1] >= 0
 uType2: uType2[4,1] >= 0
 uType2: uT2CanA[4,2] - uT2CanB[4,2] - uT2CanC[5,2] + uType2[4,2] >= 0
 uType2: uType2[4,2] >= 0
 uType2: uT2CanA[4,3] - uT2CanB[4,3] - uT2CanC[5,3] + uType2[4,3] >= 0
 uType2: uType2[4,3] >= 0
 uType2: uT2CanA[5,0] - uT2CanB[5,0] - uT2CanC[6,0] + uType2[5,0] >= 0
 uType2: uType2[5,0] >= 0
 uType2: uT2CanA[5,1] - uT2CanB[5,1] - uT2CanC[6,1] + uType2[5,1] >= 0
 uType2: uType2[5,1] >= 0
 uType2: uT2CanA[5,2] - uT2CanB[5,2] - uT2CanC[6,2] + uType2[5,2] >= 0
 uType2: uType2[5,2] >= 0
 uType2: uT2CanA[5,3] - uT2CanB[5,3] - uT2CanC[6,3] + uType2[5,3] >= 0
 uType2: uType2[5,3] >= 0
 uType2: uT2CanA[6,0] - uT2CanB[6,0] - uT2CanC[7,0] + uType2[6,0] >= 0
 uType2: uType2[6,0] >= 0
 uType2: uT2CanA[6,1] - uT2CanB[6,1] - uT2CanC[7,1] + uType2[6,1] >= 0
 uType2: uType2[6,1] >= 0
 uType2: uT2CanA[6,2] - uT2CanB[6,2] - uT2CanC[7,2] + uType2[6,2] >= 0
 uType2: uType2[6,2] >= 0
 uType2: uT2CanA[6,3] - uT2CanB[6,3] - uT2CanC[7,3] + uType2[6,3] >= 0
 uType2: uType2[6,3] >= 0
 uCAN: 2 uLANE[0] + 2 uLANE[1] + 2 uLANE[2] + 2 uLANE[3] + 2 uLANE[4]
   + 2 uLANE[5] + 2 uLANE[6] + 2 uLANE[7] + 2 uLANE[8] + 2 uLANE[9]
   + 2 uLANE[10] + 2 uLANE[11] + 2 uLANE[12] + 2 uLANE[13] + 2 uLANE[14]
   + 2 uLANE[15] - uType1[0] - uType1[1] - uType1[2] - uType1[3]
   - uType1[4] - uType1[5] - uType1[6] - uType1[7] - uType1[8] - uType1[9]
   - uType1[10] - uType1[11] - uType1[12] - uType1[13] - uType1[14]
   - uType1[15] - uType2[0,0] - uType2[0,1] - uType2[0,2] - uType2[0,3]
   - uType2[1,0] - uType2[1,1] - uType2[1,2] - uType2[1,3] - uType2[2,0]
   - uType2[2,1] - uType2[2,2] - uType2[2,3] - uType2[3,0] - uType2[3,1]
   - uType2[3,2] - uType2[3,3] - uType2[4,0] - uType2[4,1] - uType2[4,2]
   - uType2[4,3] - uType2[5,0] - uType2[5,1] - uType2[5,2] - uType2[5,3]
   - uType2[6,0] - uType2[6,1] - uType2[6,2] - uType2[6,3] - uType2[7,0]
   - uType2[7,1] - uType2[7,2] - uType2[7,3] >= 1
 bDet_ATK: uGstk[0,0] - uDetX[0,0] >= 0
 bDet_ATK: - uDetX[0,0] + uDetW[-1,0] >= 0
 bDet_ATK: - uGstk[0,0] + uDetX[0,0] - uDetW[-1,0] >= -1
 bDet_ATK: uGstk[0,1] - uDetX[0,1] >= 0
 bDet_ATK: - uDetX[0,1] + uDetW[-1,1] >= 0
 bDet_ATK: - uGstk[0,1] + uDetX[0,1] - uDetW[-1,1] >= -1
 bDet_ATK: uGstk[0,2] - uDetX[0,2] >= 0
 bDet_ATK: - uDetX[0,2] + uDetW[-1,2] >= 0
 bDet_ATK: - uGstk[0,2] + uDetX[0,2] - uDetW[-1,2] >= -1
 bDet_ATK: uGstk[0,3] - uDetX[0,3] >= 0
 bDet_ATK: - uDetX[0,3] + uDetW[-1,3] >= 0
 bDet_ATK: - uGstk[0,3] + uDetX[0,3] - uDetW[-1,3] >= -1
 bDet_ATK: uGstk[0,4] - uDetX[0,4] >= 0
 bDet_ATK: - uDetX[0,4] + uDetW[-1,4] >= 0
 bDet_ATK: - uGstk[0,4] + uDetX[0,4] - uDetW[-1,4] >= -1
 bDet_ATK: uGstk[0,5] - uDetX[0,5] >= 0
 bDet_ATK: - uDetX[0,5] + uDetW[-1,5] >= 0
 bDet_ATK: - uGstk[0,5] + uDetX[0,5] - uDetW[-1,5] >= -1
 bDet_ATK: uGstk[0,6] - uDetX[0,6] >= 0
 bDet_ATK: - uDetX[0,6] + uDetW[-1,6] >= 0
 bDet_ATK: - uGstk[0,6] + uDetX[0,6] - uDetW[-1,6] >= -1
 bDet_ATK: uGstk[0,7] - uDetX[0,7] >= 0
 bDet_ATK: - uDetX[0,7] + uDetW[-1,7] >= 0
 bDet_ATK: - uGstk[0,7] + uDetX[0,7] - uDetW[-1,7] >= -1
 bDet_ATK: uGstk[0,8] - uDetX[0,8] >= 0
 bDet_ATK: - uDetX[0,8] + uDetW[-1,8] >= 0
 bDet_ATK: - uGstk[0,8] + uDetX[0,8] - uDetW[-1,8] >= -1
 bDet_ATK: uGstk[0,9] - uDetX[0,9] >= 0
 bDet_ATK: - uDetX[0,9] + uDetW[-1,9] >= 0
 bDet_ATK: - uGstk[0,9] + uDetX[0,9] - uDetW[-1,9] >= -1
 bDet_ATK: uGstk[0,10] - uDetX[0,10] >= 0
 bDet_ATK: - uDetX[0,10] + uDetW[-1,10] >= 0
 bDet_ATK: - uGstk[0,10] + uDetX[0,10] - uDetW[-1,10] >= -1
 bDet_ATK: uGstk[0,11] - uDetX[0,11] >= 0
 bDet_ATK: - uDetX[0,11] + uDetW[-1,11] >= 0
 bDet_ATK: - uGstk[0,11] + uDetX[0,11] - uDetW[-1,11] >= -1
 bDet_ATK: uGstk[0,12] - uDetX[0,12] >= 0
 bDet_ATK: - uDetX[0,12] + uDetW[-1,12] >= 0
 bDet_ATK: - uGstk[0,12] + uDetX[0,12] - uDetW[-1,12] >= -1
 bDet_ATK: uGstk[0,13] - uDetX[0,13] >= 0
 bDet_ATK: - uDetX[0,13] + uDetW[-1,13] >= 0
 bDet_ATK: - uGstk[0,13] + uDetX[0,13] - uDetW[-1,13] >= -1
 bDet_ATK: uGstk[0,14] - uDetX[0,14] >= 0
 bDet_ATK: - uDetX[0,14] + uDetW[-1,14] >= 0
 bDet_ATK: - uGstk[0,14] + uDetX[0,14] - uDetW[-1,14] >= -1
 bDet_ATK: uGstk[0,15] - uDetX[0,15] >= 0
 bDet_ATK: - uDetX[0,15] + uDetW[-1,15] >= 0
 bDet_ATK: - uGstk[0,15] + uDetX[0,15] - uDetW[-1,15] >= -1
 bDet_SC: - uDetX[0,0] + uDetY[0,0] = 0
 bFr_SC: uDY[0,0] - udy[0,0] - uFrSB[0,0] >= 0
 bFr_SC: - uFrSB[0,0] + uDetY[0,0] >= 0
 bFr_SC: - uDY[0,0] + udy[0,0] + uFrSB[0,0] - uDetY[0,0] >= -1
 bDet_SC: - uDetX[0,1] + uDetY[0,1] = 0
 bFr_SC: uDY[0,1] - udy[0,1] - uFrSB[0,1] >= 0
 bFr_SC: - uFrSB[0,1] + uDetY[0,1] >= 0
 bFr_SC: - uDY[0,1] + udy[0,1] + uFrSB[0,1] - uDetY[0,1] >= -1
 bDet_SC: - uDetX[0,2] + uDetY[0,2] = 0
 bFr_SC: uDY[0,2] - udy[0,2] - uFrSB[0,2] >= 0
 bFr_SC: - uFrSB[0,2] + uDetY[0,2] >= 0
 bFr_SC: - uDY[0,2] + udy[0,2] + uFrSB[0,2] - uDetY[0,2] >= -1
 bDet_SC: - uDetX[0,3] + uDetY[0,3] = 0
 bFr_SC: uDY[0,3] - udy[0,3] - uFrSB[0,3] >= 0
 bFr_SC: - uFrSB[0,3] + uDetY[0,3] >= 0
 bFr_SC: - uDY[0,3] + udy[0,3] + uFrSB[0,3] - uDetY[0,3] >= -1
 bDet_SC: - uDetX[0,4] + uDetY[0,4] = 0
 bFr_SC: uDY[0,4] - udy[0,4] - uFrSB[0,4] >= 0
 bFr_SC: - uFrSB[0,4] + uDetY[0,4] >= 0
 bFr_SC: - uDY[0,4] + udy[0,4] + uFrSB[0,4] - uDetY[0,4] >= -1
 bDet_SC: - uDetX[0,5] + uDetY[0,5] = 0
 bFr_SC: uDY[0,5] - udy[0,5] - uFrSB[0,5] >= 0
 bFr_SC: - uFrSB[0,5] + uDetY[0,5] >= 0
 bFr_SC: - uDY[0,5] + udy[0,5] + uFrSB[0,5] - uDetY[0,5] >= -1
 bDet_SC: - uDetX[0,6] + uDetY[0,6] = 0
 bFr_SC: uDY[0,6] - udy[0,6] - uFrSB[0,6] >= 0
 bFr_SC: - uFrSB[0,6] + uDetY[0,6] >= 0
 bFr_SC: - uDY[0,6] + udy[0,6] + uFrSB[0,6] - uDetY[0,6] >= -1
 bDet_SC: - uDetX[0,7] + uDetY[0,7] = 0
 bFr_SC: uDY[0,7] - udy[0,7] - uFrSB[0,7] >= 0
 bFr_SC: - uFrSB[0,7] + uDetY[0,7] >= 0
 bFr_SC: - uDY[0,7] + udy[0,7] + uFrSB[0,7] - uDetY[0,7] >= -1
 bDet_SC: - uDetX[0,8] + uDetY[0,8] = 0
 bFr_SC: uDY[0,8] - udy[0,8] - uFrSB[0,8] >= 0
 bFr_SC: - uFrSB[0,8] + uDetY[0,8] >= 0
 bFr_SC: - uDY[0,8] + udy[0,8] + uFrSB[0,8] - uDetY[0,8] >= -1
 bDet_SC: - uDetX[0,9] + uDetY[0,9] = 0
 bFr_SC: uDY[0,9] - udy[0,9] - uFrSB[0,9] >= 0
 bFr_SC: - uFrSB[0,9] + uDetY[0,9] >= 0
 bFr_SC: - uDY[0,9] + udy[0,9] + uFrSB[0,9] - uDetY[0,9] >= -1
 bDet_SC: - uDetX[0,10] + uDetY[0,10] = 0
 bFr_SC: uDY[0,10] - udy[0,10] - uFrSB[0,10] >= 0
 bFr_SC: - uFrSB[0,10] + uDetY[0,10] >= 0
 bFr_SC: - uDY[0,10] + udy[0,10] + uFrSB[0,10] - uDetY[0,10] >= -1
 bDet_SC: - uDetX[0,11] + uDetY[0,11] = 0
 bFr_SC: uDY[0,11] - udy[0,11] - uFrSB[0,11] >= 0
 bFr_SC: - uFrSB[0,11] + uDetY[0,11] >= 0
 bFr_SC: - uDY[0,11] + udy[0,11] + uFrSB[0,11] - uDetY[0,11] >= -1
 bDet_SC: - uDetX[0,12] + uDetY[0,12] = 0
 bFr_SC: uDY[0,12] - udy[0,12] - uFrSB[0,12] >= 0
 bFr_SC: - uFrSB[0,12] + uDetY[0,12] >= 0
 bFr_SC: - uDY[0,12] + udy[0,12] + uFrSB[0,12] - uDetY[0,12] >= -1
 bDet_SC: - uDetX[0,13] + uDetY[0,13] = 0
 bFr_SC: uDY[0,13] - udy[0,13] - uFrSB[0,13] >= 0
 bFr_SC: - uFrSB[0,13] + uDetY[0,13] >= 0
 bFr_SC: - uDY[0,13] + udy[0,13] + uFrSB[0,13] - uDetY[0,13] >= -1
 bDet_SC: - uDetX[0,14] + uDetY[0,14] = 0
 bFr_SC: uDY[0,14] - udy[0,14] - uFrSB[0,14] >= 0
 bFr_SC: - uFrSB[0,14] + uDetY[0,14] >= 0
 bFr_SC: - uDY[0,14] + udy[0,14] + uFrSB[0,14] - uDetY[0,14] >= -1
 bDet_SC: - uDetX[0,15] + uDetY[0,15] = 0
 bFr_SC: uDY[0,15] - udy[0,15] - uFrSB[0,15] >= 0
 bFr_SC: - uFrSB[0,15] + uDetY[0,15] >= 0
 bFr_SC: - uDY[0,15] + udy[0,15] + uFrSB[0,15] - uDetY[0,15] >= -1
 bDet_SR: - uDetY[0,0] + uDetZ[0,0] = 0
 bDet_SR: - uDetY[0,5] + uDetZ[0,1] = 0
 bDet_SR: - uDetY[0,10] + uDetZ[0,2] = 0
 bDet_SR: - uDetY[0,15] + uDetZ[0,3] = 0
 bDet_SR: - uDetY[0,4] + uDetZ[0,4] = 0
 bDet_SR: - uDetY[0,9] + uDetZ[0,5] = 0
 bDet_SR: - uDetY[0,14] + uDetZ[0,6] = 0
 bDet_SR: - uDetY[0,3] + uDetZ[0,7] = 0
 bDet_SR: - uDetY[0,8] + uDetZ[0,8] = 0
 bDet_SR: - uDetY[0,13] + uDetZ[0,9] = 0
 bDet_SR: - uDetY[0,2] + uDetZ[0,10] = 0
 bDet_SR: - uDetY[0,7] + uDetZ[0,11] = 0
 bDet_SR: - uDetY[0,12] + uDetZ[0,12] = 0
 bDet_SR: - uDetY[0,1] + uDetZ[0,13] = 0
 bDet_SR: - uDetY[0,6] + uDetZ[0,14] = 0
 bDet_SR: - uDetY[0,11] + uDetZ[0,15] = 0
 bDet_MC: - uDetZ[0,0] + uDetC[0,0] <= 0
 bDet_MC: - uDetZ[0,1] + uDetC[0,0] <= 0
 bDet_MC: - uDetZ[0,2] + uDetC[0,0] <= 0
 bDet_MC: - uDetZ[0,3] + uDetC[0,0] <= 0
 bDet_MC: uDetZ[0,0] + uDetZ[0,1] + uDetZ[0,2] + uDetZ[0,3] - uDetC[0,0]
   <= 3
 bDet_MC: - uDetW[0,0] + uDetC[0,0] = 0
 bDet_MC: - uDetW[0,1] + uDetC[0,0] = 0
 bDet_MC: - uDetW[0,2] + uDetC[0,0] = 0
 bDet_MC: - uDetW[0,3] + uDetC[0,0] = 0
 bDet_MC: - uDetZ[0,4] + uDetC[0,1] <= 0
 bDet_MC: - uDetZ[0,5] + uDetC[0,1] <= 0
 bDet_MC: - uDetZ[0,6] + uDetC[0,1] <= 0
 bDet_MC: - uDetZ[0,7] + uDetC[0,1] <= 0
 bDet_MC: uDetZ[0,4] + uDetZ[0,5] + uDetZ[0,6] + uDetZ[0,7] - uDetC[0,1]
   <= 3
 bDet_MC: - uDetW[0,4] + uDetC[0,1] = 0
 bDet_MC: - uDetW[0,5] + uDetC[0,1] = 0
 bDet_MC: - uDetW[0,6] + uDetC[0,1] = 0
 bDet_MC: - uDetW[0,7] + uDetC[0,1] = 0
 bDet_MC: - uDetZ[0,8] + uDetC[0,2] <= 0
 bDet_MC: - uDetZ[0,9] + uDetC[0,2] <= 0
 bDet_MC: - uDetZ[0,10] + uDetC[0,2] <= 0
 bDet_MC: - uDetZ[0,11] + uDetC[0,2] <= 0
 bDet_MC: uDetZ[0,8] + uDetZ[0,9] + uDetZ[0,10] + uDetZ[0,11] - uDetC[0,2]
   <= 3
 bDet_MC: - uDetW[0,8] + uDetC[0,2] = 0
 bDet_MC: - uDetW[0,9] + uDetC[0,2] = 0
 bDet_MC: - uDetW[0,10] + uDetC[0,2] = 0
 bDet_MC: - uDetW[0,11] + uDetC[0,2] = 0
 bDet_MC: - uDetZ[0,12] + uDetC[0,3] <= 0
 bDet_MC: - uDetZ[0,13] + uDetC[0,3] <= 0
 bDet_MC: - uDetZ[0,14] + uDetC[0,3] <= 0
 bDet_MC: - uDetZ[0,15] + uDetC[0,3] <= 0
 bDet_MC: uDetZ[0,12] + uDetZ[0,13] + uDetZ[0,14] + uDetZ[0,15]
   - uDetC[0,3] <= 3
 bDet_MC: - uDetW[0,12] + uDetC[0,3] = 0
 bDet_MC: - uDetW[0,13] + uDetC[0,3] = 0
 bDet_MC: - uDetW[0,14] + uDetC[0,3] = 0
 bDet_MC: - uDetW[0,15] + uDetC[0,3] = 0
 bDet_MC: - uDetZ[0,0] + uDiffDetMC[0,0] >= 0
 bDet_MC: - uDZ[0,0] + uDetZ[0,0] - uDiffDetMC[0,0] >= -1
 bDet_MC: uDZ[0,0] + uDiffDetMC[0,0] >= 1
 bDet_MC: - uDetZ[0,1] + uDiffDetMC[0,1] >= 0
 bDet_MC: - uDZ[0,1] + uDetZ[0,1] - uDiffDetMC[0,1] >= -1
 bDet_MC: uDZ[0,1] + uDiffDetMC[0,1] >= 1
 bDet_MC: - uDetZ[0,2] + uDiffDetMC[0,2] >= 0
 bDet_MC: - uDZ[0,2] + uDetZ[0,2] - uDiffDetMC[0,2] >= -1
 bDet_MC: uDZ[0,2] + uDiffDetMC[0,2] >= 1
 bDet_MC: - uDetZ[0,3] + uDiffDetMC[0,3] >= 0
 bDet_MC: - uDZ[0,3] + uDetZ[0,3] - uDiffDetMC[0,3] >= -1
 bDet_MC: uDZ[0,3] + uDiffDetMC[0,3] >= 1
 bDet_MC: - uDetZ[0,4] + uDiffDetMC[0,4] >= 0
 bDet_MC: - uDZ[0,4] + uDetZ[0,4] - uDiffDetMC[0,4] >= -1
 bDet_MC: uDZ[0,4] + uDiffDetMC[0,4] >= 1
 bDet_MC: - uDetZ[0,5] + uDiffDetMC[0,5] >= 0
 bDet_MC: - uDZ[0,5] + uDetZ[0,5] - uDiffDetMC[0,5] >= -1
 bDet_MC: uDZ[0,5] + uDiffDetMC[0,5] >= 1
 bDet_MC: - uDetZ[0,6] + uDiffDetMC[0,6] >= 0
 bDet_MC: - uDZ[0,6] + uDetZ[0,6] - uDiffDetMC[0,6] >= -1
 bDet_MC: uDZ[0,6] + uDiffDetMC[0,6] >= 1
 bDet_MC: - uDetZ[0,7] + uDiffDetMC[0,7] >= 0
 bDet_MC: - uDZ[0,7] + uDetZ[0,7] - uDiffDetMC[0,7] >= -1
 bDet_MC: uDZ[0,7] + uDiffDetMC[0,7] >= 1
 bDet_MC: - uDetZ[0,8] + uDiffDetMC[0,8] >= 0
 bDet_MC: - uDZ[0,8] + uDetZ[0,8] - uDiffDetMC[0,8] >= -1
 bDet_MC: uDZ[0,8] + uDiffDetMC[0,8] >= 1
 bDet_MC: - uDetZ[0,9] + uDiffDetMC[0,9] >= 0
 bDet_MC: - uDZ[0,9] + uDetZ[0,9] - uDiffDetMC[0,9] >= -1
 bDet_MC: uDZ[0,9] + uDiffDetMC[0,9] >= 1
 bDet_MC: - uDetZ[0,10] + uDiffDetMC[0,10] >= 0
 bDet_MC: - uDZ[0,10] + uDetZ[0,10] - uDiffDetMC[0,10] >= -1
 bDet_MC: uDZ[0,10] + uDiffDetMC[0,10] >= 1
 bDet_MC: - uDetZ[0,11] + uDiffDetMC[0,11] >= 0
 bDet_MC: - uDZ[0,11] + uDetZ[0,11] - uDiffDetMC[0,11] >= -1
 bDet_MC: uDZ[0,11] + uDiffDetMC[0,11] >= 1
 bDet_MC: - uDetZ[0,12] + uDiffDetMC[0,12] >= 0
 bDet_MC: - uDZ[0,12] + uDetZ[0,12] - uDiffDetMC[0,12] >= -1
 bDet_MC: uDZ[0,12] + uDiffDetMC[0,12] >= 1
 bDet_MC: - uDetZ[0,13] + uDiffDetMC[0,13] >= 0
 bDet_MC: - uDZ[0,13] + uDetZ[0,13] - uDiffDetMC[0,13] >= -1
 bDet_MC: uDZ[0,13] + uDiffDetMC[0,13] >= 1
 bDet_MC: - uDetZ[0,14] + uDiffDetMC[0,14] >= 0
 bDet_MC: - uDZ[0,14] + uDetZ[0,14] - uDiffDetMC[0,14] >= -1
 bDet_MC: uDZ[0,14] + uDiffDetMC[0,14] >= 1
 bDet_MC: - uDetZ[0,15] + uDiffDetMC[0,15] >= 0
 bDet_MC: - uDZ[0,15] + uDetZ[0,15] - uDiffDetMC[0,15] >= -1
 bDet_MC: uDZ[0,15] + uDiffDetMC[0,15] >= 1
 bDet_MC: - uDiffDetMC[0,0] + uDiffDetC[0,0] <= 0
 bDet_MC: - uDiffDetMC[0,1] + uDiffDetC[0,0] <= 0
 bDet_MC: - uDiffDetMC[0,2] + uDiffDetC[0,0] <= 0
 bDet_MC: - uDiffDetMC[0,3] + uDiffDetC[0,0] <= 0
 bDet_MC: uDiffDetMC[0,0] + uDiffDetMC[0,1] + uDiffDetMC[0,2]
   + uDiffDetMC[0,3] - uDiffDetC[0,0] <= 3
 bDet_MC: - uDiffDetMC[0,4] + uDiffDetC[0,1] <= 0
 bDet_MC: - uDiffDetMC[0,5] + uDiffDetC[0,1] <= 0
 bDet_MC: - uDiffDetMC[0,6] + uDiffDetC[0,1] <= 0
 bDet_MC: - uDiffDetMC[0,7] + uDiffDetC[0,1] <= 0
 bDet_MC: uDiffDetMC[0,4] + uDiffDetMC[0,5] + uDiffDetMC[0,6]
   + uDiffDetMC[0,7] - uDiffDetC[0,1] <= 3
 bDet_MC: - uDiffDetMC[0,8] + uDiffDetC[0,2] <= 0
 bDet_MC: - uDiffDetMC[0,9] + uDiffDetC[0,2] <= 0
 bDet_MC: - uDiffDetMC[0,10] + uDiffDetC[0,2] <= 0
 bDet_MC: - uDiffDetMC[0,11] + uDiffDetC[0,2] <= 0
 bDet_MC: uDiffDetMC[0,8] + uDiffDetMC[0,9] + uDiffDetMC[0,10]
   + uDiffDetMC[0,11] - uDiffDetC[0,2] <= 3
 bDet_MC: - uDiffDetMC[0,12] + uDiffDetC[0,3] <= 0
 bDet_MC: - uDiffDetMC[0,13] + uDiffDetC[0,3] <= 0
 bDet_MC: - uDiffDetMC[0,14] + uDiffDetC[0,3] <= 0
 bDet_MC: - uDiffDetMC[0,15] + uDiffDetC[0,3] <= 0
 bDet_MC: uDiffDetMC[0,12] + uDiffDetMC[0,13] + uDiffDetMC[0,14]
   + uDiffDetMC[0,15] - uDiffDetC[0,3] <= 3
 bFr_MC: uTC[0,0] - uFrMC[0,0] >= 0
 bFr_MC: udw[0,0] + uFrMC[0,0] <= 1
 bFr_MC: - udw[0,0] + uTC[0,0] + uDiffDetC[0,0] - uFrMC[0,0] <= 1
 bFr_MC: uDiffDetC[0,0] - uFrMC[0,0] >= 0
 bFr_MC: uTC[0,0] - uFrMC[0,1] >= 0
 bFr_MC: udw[0,1] + uFrMC[0,1] <= 1
 bFr_MC: - udw[0,1] + uTC[0,0] + uDiffDetC[0,0] - uFrMC[0,1] <= 1
 bFr_MC: uDiffDetC[0,0] - uFrMC[0,1] >= 0
 bFr_MC: uTC[0,0] - uFrMC[0,2] >= 0
 bFr_MC: udw[0,2] + uFrMC[0,2] <= 1
 bFr_MC: - udw[0,2] + uTC[0,0] + uDiffDetC[0,0] - uFrMC[0,2] <= 1
 bFr_MC: uDiffDetC[0,0] - uFrMC[0,2] >= 0
 bFr_MC: uTC[0,0] - uFrMC[0,3] >= 0
 bFr_MC: udw[0,3] + uFrMC[0,3] <= 1
 bFr_MC: - udw[0,3] + uTC[0,0] + uDiffDetC[0,0] - uFrMC[0,3] <= 1
 bFr_MC: uDiffDetC[0,0] - uFrMC[0,3] >= 0
 bFr_MC: uTC[0,1] - uFrMC[0,4] >= 0
 bFr_MC: udw[0,4] + uFrMC[0,4] <= 1
 bFr_MC: - udw[0,4] + uTC[0,1] + uDiffDetC[0,1] - uFrMC[0,4] <= 1
 bFr_MC: uDiffDetC[0,1] - uFrMC[0,4] >= 0
 bFr_MC: uTC[0,1] - uFrMC[0,5] >= 0
 bFr_MC: udw[0,5] + uFrMC[0,5] <= 1
 bFr_MC: - udw[0,5] + uTC[0,1] + uDiffDetC[0,1] - uFrMC[0,5] <= 1
 bFr_MC: uDiffDetC[0,1] - uFrMC[0,5] >= 0
 bFr_MC: uTC[0,1] - uFrMC[0,6] >= 0
 bFr_MC: udw[0,6] + uFrMC[0,6] <= 1
 bFr_MC: - udw[0,6] + uTC[0,1] + uDiffDetC[0,1] - uFrMC[0,6] <= 1
 bFr_MC: uDiffDetC[0,1] - uFrMC[0,6] >= 0
 bFr_MC: uTC[0,1] - uFrMC[0,7] >= 0
 bFr_MC: udw[0,7] + uFrMC[0,7] <= 1
 bFr_MC: - udw[0,7] + uTC[0,1] + uDiffDetC[0,1] - uFrMC[0,7] <= 1
 bFr_MC: uDiffDetC[0,1] - uFrMC[0,7] >= 0
 bFr_MC: uTC[0,2] - uFrMC[0,8] >= 0
 bFr_MC: udw[0,8] + uFrMC[0,8] <= 1
 bFr_MC: - udw[0,8] + uTC[0,2] + uDiffDetC[0,2] - uFrMC[0,8] <= 1
 bFr_MC: uDiffDetC[0,2] - uFrMC[0,8] >= 0
 bFr_MC: uTC[0,2] - uFrMC[0,9] >= 0
 bFr_MC: udw[0,9] + uFrMC[0,9] <= 1
 bFr_MC: - udw[0,9] + uTC[0,2] + uDiffDetC[0,2] - uFrMC[0,9] <= 1
 bFr_MC: uDiffDetC[0,2] - uFrMC[0,9] >= 0
 bFr_MC: uTC[0,2] - uFrMC[0,10] >= 0
 bFr_MC: udw[0,10] + uFrMC[0,10] <= 1
 bFr_MC: - udw[0,10] + uTC[0,2] + uDiffDetC[0,2] - uFrMC[0,10] <= 1
 bFr_MC: uDiffDetC[0,2] - uFrMC[0,10] >= 0
 bFr_MC: uTC[0,2] - uFrMC[0,11] >= 0
 bFr_MC: udw[0,11] + uFrMC[0,11] <= 1
 bFr_MC: - udw[0,11] + uTC[0,2] + uDiffDetC[0,2] - uFrMC[0,11] <= 1
 bFr_MC: uDiffDetC[0,2] - uFrMC[0,11] >= 0
 bFr_MC: uTC[0,3] - uFrMC[0,12] >= 0
 bFr_MC: udw[0,12] + uFrMC[0,12] <= 1
 bFr_MC: - udw[0,12] + uTC[0,3] + uDiffDetC[0,3] - uFrMC[0,12] <= 1
 bFr_MC: uDiffDetC[0,3] - uFrMC[0,12] >= 0
 bFr_MC: uTC[0,3] - uFrMC[0,13] >= 0
 bFr_MC: udw[0,13] + uFrMC[0,13] <= 1
 bFr_MC: - udw[0,13] + uTC[0,3] + uDiffDetC[0,3] - uFrMC[0,13] <= 1
 bFr_MC: uDiffDetC[0,3] - uFrMC[0,13] >= 0
 bFr_MC: uTC[0,3] - uFrMC[0,14] >= 0
 bFr_MC: udw[0,14] + uFrMC[0,14] <= 1
 bFr_MC: - udw[0,14] + uTC[0,3] + uDiffDetC[0,3] - uFrMC[0,14] <= 1
 bFr_MC: uDiffDetC[0,3] - uFrMC[0,14] >= 0
 bFr_MC: uTC[0,3] - uFrMC[0,15] >= 0
 bFr_MC: udw[0,15] + uFrMC[0,15] <= 1
 bFr_MC: - udw[0,15] + uTC[0,3] + uDiffDetC[0,3] - uFrMC[0,15] <= 1
 bFr_MC: uDiffDetC[0,3] - uFrMC[0,15] >= 0
 lBas: lDX[2,0] - ldx[2,0] >= 0
 lBas: lDY[1,0] - ldy[1,0] >= 0
 lBas: lDZ[1,0] - ldz[1,0] >= 0
 lBas: lDW[1,0] - ldw[1,0] >= 0
 lBas: lDX[3,0] - ldx[3,0] >= 0
 lBas: lDY[2,0] - ldy[2,0] >= 0
 lBas: lDZ[2,0] - ldz[2,0] >= 0
 lBas: lDW[2,0] - ldw[2,0] >= 0
 lBas: lDX[4,0] - ldx[4,0] >= 0
 lBas: lDY[3,0] - ldy[3,0] >= 0
 lBas: lDZ[3,0] - ldz[3,0] >= 0
 lBas: lDW[3,0] - ldw[3,0] >= 0
 lBas: lDX[5,0] - ldx[5,0] >= 0
 lBas: lDY[4,0] - ldy[4,0] >= 0
 lBas: lDZ[4,0] - ldz[4,0] >= 0
 lBas: lDW[4,0] - ldw[4,0] >= 0
 lBas: lDX[6,0] - ldx[6,0] >= 0
 lBas: lDY[5,0] - ldy[5,0] >= 0
 lBas: lDZ[5,0] - ldz[5,0] >= 0
 lBas: lDW[5,0] - ldw[5,0] >= 0
 lBas: lDX[7,0] - ldx[7,0] >= 0
 lBas: lDY[6,0] - ldy[6,0] >= 0
 lBas: lDZ[6,0] - ldz[6,0] >= 0
 lBas: lDW[6,0] - ldw[6,0] >= 0
 lBas: lDX[8,0] - ldx[8,0] >= 0
 lBas: lDY[7,0] - ldy[7,0] >= 0
 lBas: lDZ[7,0] - ldz[7,0] >= 0
 lBas: lDW[7,0] - ldw[7,0] >= 0
 lBas: lDX[2,1] - ldx[2,1] >= 0
 lBas: lDY[1,1] - ldy[1,1] >= 0
 lBas: lDZ[1,1] - ldz[1,1] >= 0
 lBas: lDW[1,1] - ldw[1,1] >= 0
 lBas: lDX[3,1] - ldx[3,1] >= 0
 lBas: lDY[2,1] - ldy[2,1] >= 0
 lBas: lDZ[2,1] - ldz[2,1] >= 0
 lBas: lDW[2,1] - ldw[2,1] >= 0
 lBas: lDX[4,1] - ldx[4,1] >= 0
 lBas: lDY[3,1] - ldy[3,1] >= 0
 lBas: lDZ[3,1] - ldz[3,1] >= 0
 lBas: lDW[3,1] - ldw[3,1] >= 0
 lBas: lDX[5,1] - ldx[5,1] >= 0
 lBas: lDY[4,1] - ldy[4,1] >= 0
 lBas: lDZ[4,1] - ldz[4,1] >= 0
 lBas: lDW[4,1] - ldw[4,1] >= 0
 lBas: lDX[6,1] - ldx[6,1] >= 0
 lBas: lDY[5,1] - ldy[5,1] >= 0
 lBas: lDZ[5,1] - ldz[5,1] >= 0
 lBas: lDW[5,1] - ldw[5,1] >= 0
 lBas: lDX[7,1] - ldx[7,1] >= 0
 lBas: lDY[6,1] - ldy[6,1] >= 0
 lBas: lDZ[6,1] - ldz[6,1] >= 0
 lBas: lDW[6,1] - ldw[6,1] >= 0
 lBas: lDX[8,1] - ldx[8,1] >= 0
 lBas: lDY[7,1] - ldy[7,1] >= 0
 lBas: lDZ[7,1] - ldz[7,1] >= 0
 lBas: lDW[7,1] - ldw[7,1] >= 0
 lBas: lDX[2,2] - ldx[2,2] >= 0
 lBas: lDY[1,2] - ldy[1,2] >= 0
 lBas: lDZ[1,2] - ldz[1,2] >= 0
 lBas: lDW[1,2] - ldw[1,2] >= 0
 lBas: lDX[3,2] - ldx[3,2] >= 0
 lBas: lDY[2,2] - ldy[2,2] >= 0
 lBas: lDZ[2,2] - ldz[2,2] >= 0
 lBas: lDW[2,2] - ldw[2,2] >= 0
 lBas: lDX[4,2] - ldx[4,2] >= 0
 lBas: lDY[3,2] - ldy[3,2] >= 0
 lBas: lDZ[3,2] - ldz[3,2] >= 0
 lBas: lDW[3,2] - ldw[3,2] >= 0
 lBas: lDX[5,2] - ldx[5,2] >= 0
 lBas: lDY[4,2] - ldy[4,2] >= 0
 lBas: lDZ[4,2] - ldz[4,2] >= 0
 lBas: lDW[4,2] - ldw[4,2] >= 0
 lBas: lDX[6,2] - ldx[6,2] >= 0
 lBas: lDY[5,2] - ldy[5,2] >= 0
 lBas: lDZ[5,2] - ldz[5,2] >= 0
 lBas: lDW[5,2] - ldw[5,2] >= 0
 lBas: lDX[7,2] - ldx[7,2] >= 0
 lBas: lDY[6,2] - ldy[6,2] >= 0
 lBas: lDZ[6,2] - ldz[6,2] >= 0
 lBas: lDW[6,2] - ldw[6,2] >= 0
 lBas: lDX[8,2] - ldx[8,2] >= 0
 lBas: lDY[7,2] - ldy[7,2] >= 0
 lBas: lDZ[7,2] - ldz[7,2] >= 0
 lBas: lDW[7,2] - ldw[7,2] >= 0
 lBas: lDX[2,3] - ldx[2,3] >= 0
 lBas: lDY[1,3] - ldy[1,3] >= 0
 lBas: lDZ[1,3] - ldz[1,3] >= 0
 lBas: lDW[1,3] - ldw[1,3] >= 0
 lBas: lDX[3,3] - ldx[3,3] >= 0
 lBas: lDY[2,3] - ldy[2,3] >= 0
 lBas: lDZ[2,3] - ldz[2,3] >= 0
 lBas: lDW[2,3] - ldw[2,3] >= 0
 lBas: lDX[4,3] - ldx[4,3] >= 0
 lBas: lDY[3,3] - ldy[3,3] >= 0
 lBas: lDZ[3,3] - ldz[3,3] >= 0
 lBas: lDW[3,3] - ldw[3,3] >= 0
 lBas: lDX[5,3] - ldx[5,3] >= 0
 lBas: lDY[4,3] - ldy[4,3] >= 0
 lBas: lDZ[4,3] - ldz[4,3] >= 0
 lBas: lDW[4,3] - ldw[4,3] >= 0
 lBas: lDX[6,3] - ldx[6,3] >= 0
 lBas: lDY[5,3] - ldy[5,3] >= 0
 lBas: lDZ[5,3] - ldz[5,3] >= 0
 lBas: lDW[5,3] - ldw[5,3] >= 0
 lBas: lDX[7,3] - ldx[7,3] >= 0
 lBas: lDY[6,3] - ldy[6,3] >= 0
 lBas: lDZ[6,3] - ldz[6,3] >= 0
 lBas: lDW[6,3] - ldw[6,3] >= 0
 lBas: lDX[8,3] - ldx[8,3] >= 0
 lBas: lDY[7,3] - ldy[7,3] >= 0
 lBas: lDZ[7,3] - ldz[7,3] >= 0
 lBas: lDW[7,3] - ldw[7,3] >= 0
 lBas: lDX[2,4] - ldx[2,4] >= 0
 lBas: lDY[1,4] - ldy[1,4] >= 0
 lBas: lDZ[1,4] - ldz[1,4] >= 0
 lBas: lDW[1,4] - ldw[1,4] >= 0
 lBas: lDX[3,4] - ldx[3,4] >= 0
 lBas: lDY[2,4] - ldy[2,4] >= 0
 lBas: lDZ[2,4] - ldz[2,4] >= 0
 lBas: lDW[2,4] - ldw[2,4] >= 0
 lBas: lDX[4,4] - ldx[4,4] >= 0
 lBas: lDY[3,4] - ldy[3,4] >= 0
 lBas: lDZ[3,4] - ldz[3,4] >= 0
 lBas: lDW[3,4] - ldw[3,4] >= 0
 lBas: lDX[5,4] - ldx[5,4] >= 0
 lBas: lDY[4,4] - ldy[4,4] >= 0
 lBas: lDZ[4,4] - ldz[4,4] >= 0
 lBas: lDW[4,4] - ldw[4,4] >= 0
 lBas: lDX[6,4] - ldx[6,4] >= 0
 lBas: lDY[5,4] - ldy[5,4] >= 0
 lBas: lDZ[5,4] - ldz[5,4] >= 0
 lBas: lDW[5,4] - ldw[5,4] >= 0
 lBas: lDX[7,4] - ldx[7,4] >= 0
 lBas: lDY[6,4] - ldy[6,4] >= 0
 lBas: lDZ[6,4] - ldz[6,4] >= 0
 lBas: lDW[6,4] - ldw[6,4] >= 0
 lBas: lDX[8,4] - ldx[8,4] >= 0
 lBas: lDY[7,4] - ldy[7,4] >= 0
 lBas: lDZ[7,4] - ldz[7,4] >= 0
 lBas: lDW[7,4] - ldw[7,4] >= 0
 lBas: lDX[2,5] - ldx[2,5] >= 0
 lBas: lDY[1,5] - ldy[1,5] >= 0
 lBas: lDZ[1,5] - ldz[1,5] >= 0
 lBas: lDW[1,5] - ldw[1,5] >= 0
 lBas: lDX[3,5] - ldx[3,5] >= 0
 lBas: lDY[2,5] - ldy[2,5] >= 0
 lBas: lDZ[2,5] - ldz[2,5] >= 0
 lBas: lDW[2,5] - ldw[2,5] >= 0
 lBas: lDX[4,5] - ldx[4,5] >= 0
 lBas: lDY[3,5] - ldy[3,5] >= 0
 lBas: lDZ[3,5] - ldz[3,5] >= 0
 lBas: lDW[3,5] - ldw[3,5] >= 0
 lBas: lDX[5,5] - ldx[5,5] >= 0
 lBas: lDY[4,5] - ldy[4,5] >= 0
 lBas: lDZ[4,5] - ldz[4,5] >= 0
 lBas: lDW[4,5] - ldw[4,5] >= 0
 lBas: lDX[6,5] - ldx[6,5] >= 0
 lBas: lDY[5,5] - ldy[5,5] >= 0
 lBas: lDZ[5,5] - ldz[5,5] >= 0
 lBas: lDW[5,5] - ldw[5,5] >= 0
 lBas: lDX[7,5] - ldx[7,5] >= 0
 lBas: lDY[6,5] - ldy[6,5] >= 0
 lBas: lDZ[6,5] - ldz[6,5] >= 0
 lBas: lDW[6,5] - ldw[6,5] >= 0
 lBas: lDX[8,5] - ldx[8,5] >= 0
 lBas: lDY[7,5] - ldy[7,5] >= 0
 lBas: lDZ[7,5] - ldz[7,5] >= 0
 lBas: lDW[7,5] - ldw[7,5] >= 0
 lBas: lDX[2,6] - ldx[2,6] >= 0
 lBas: lDY[1,6] - ldy[1,6] >= 0
 lBas: lDZ[1,6] - ldz[1,6] >= 0
 lBas: lDW[1,6] - ldw[1,6] >= 0
 lBas: lDX[3,6] - ldx[3,6] >= 0
 lBas: lDY[2,6] - ldy[2,6] >= 0
 lBas: lDZ[2,6] - ldz[2,6] >= 0
 lBas: lDW[2,6] - ldw[2,6] >= 0
 lBas: lDX[4,6] - ldx[4,6] >= 0
 lBas: lDY[3,6] - ldy[3,6] >= 0
 lBas: lDZ[3,6] - ldz[3,6] >= 0
 lBas: lDW[3,6] - ldw[3,6] >= 0
 lBas: lDX[5,6] - ldx[5,6] >= 0
 lBas: lDY[4,6] - ldy[4,6] >= 0
 lBas: lDZ[4,6] - ldz[4,6] >= 0
 lBas: lDW[4,6] - ldw[4,6] >= 0
 lBas: lDX[6,6] - ldx[6,6] >= 0
 lBas: lDY[5,6] - ldy[5,6] >= 0
 lBas: lDZ[5,6] - ldz[5,6] >= 0
 lBas: lDW[5,6] - ldw[5,6] >= 0
 lBas: lDX[7,6] - ldx[7,6] >= 0
 lBas: lDY[6,6] - ldy[6,6] >= 0
 lBas: lDZ[6,6] - ldz[6,6] >= 0
 lBas: lDW[6,6] - ldw[6,6] >= 0
 lBas: lDX[8,6] - ldx[8,6] >= 0
 lBas: lDY[7,6] - ldy[7,6] >= 0
 lBas: lDZ[7,6] - ldz[7,6] >= 0
 lBas: lDW[7,6] - ldw[7,6] >= 0
 lBas: lDX[2,7] - ldx[2,7] >= 0
 lBas: lDY[1,7] - ldy[1,7] >= 0
 lBas: lDZ[1,7] - ldz[1,7] >= 0
 lBas: lDW[1,7] - ldw[1,7] >= 0
 lBas: lDX[3,7] - ldx[3,7] >= 0
 lBas: lDY[2,7] - ldy[2,7] >= 0
 lBas: lDZ[2,7] - ldz[2,7] >= 0
 lBas: lDW[2,7] - ldw[2,7] >= 0
 lBas: lDX[4,7] - ldx[4,7] >= 0
 lBas: lDY[3,7] - ldy[3,7] >= 0
 lBas: lDZ[3,7] - ldz[3,7] >= 0
 lBas: lDW[3,7] - ldw[3,7] >= 0
 lBas: lDX[5,7] - ldx[5,7] >= 0
 lBas: lDY[4,7] - ldy[4,7] >= 0
 lBas: lDZ[4,7] - ldz[4,7] >= 0
 lBas: lDW[4,7] - ldw[4,7] >= 0
 lBas: lDX[6,7] - ldx[6,7] >= 0
 lBas: lDY[5,7] - ldy[5,7] >= 0
 lBas: lDZ[5,7] - ldz[5,7] >= 0
 lBas: lDW[5,7] - ldw[5,7] >= 0
 lBas: lDX[7,7] - ldx[7,7] >= 0
 lBas: lDY[6,7] - ldy[6,7] >= 0
 lBas: lDZ[6,7] - ldz[6,7] >= 0
 lBas: lDW[6,7] - ldw[6,7] >= 0
 lBas: lDX[8,7] - ldx[8,7] >= 0
 lBas: lDY[7,7] - ldy[7,7] >= 0
 lBas: lDZ[7,7] - ldz[7,7] >= 0
 lBas: lDW[7,7] - ldw[7,7] >= 0
 lBas: lDX[2,8] - ldx[2,8] >= 0
 lBas: lDY[1,8] - ldy[1,8] >= 0
 lBas: lDZ[1,8] - ldz[1,8] >= 0
 lBas: lDW[1,8] - ldw[1,8] >= 0
 lBas: lDX[3,8] - ldx[3,8] >= 0
 lBas: lDY[2,8] - ldy[2,8] >= 0
 lBas: lDZ[2,8] - ldz[2,8] >= 0
 lBas: lDW[2,8] - ldw[2,8] >= 0
 lBas: lDX[4,8] - ldx[4,8] >= 0
 lBas: lDY[3,8] - ldy[3,8] >= 0
 lBas: lDZ[3,8] - ldz[3,8] >= 0
 lBas: lDW[3,8] - ldw[3,8] >= 0
 lBas: lDX[5,8] - ldx[5,8] >= 0
 lBas: lDY[4,8] - ldy[4,8] >= 0
 lBas: lDZ[4,8] - ldz[4,8] >= 0
 lBas: lDW[4,8] - ldw[4,8] >= 0
 lBas: lDX[6,8] - ldx[6,8] >= 0
 lBas: lDY[5,8] - ldy[5,8] >= 0
 lBas: lDZ[5,8] - ldz[5,8] >= 0
 lBas: lDW[5,8] - ldw[5,8] >= 0
 lBas: lDX[7,8] - ldx[7,8] >= 0
 lBas: lDY[6,8] - ldy[6,8] >= 0
 lBas: lDZ[6,8] - ldz[6,8] >= 0
 lBas: lDW[6,8] - ldw[6,8] >= 0
 lBas: lDX[8,8] - ldx[8,8] >= 0
 lBas: lDY[7,8] - ldy[7,8] >= 0
 lBas: lDZ[7,8] - ldz[7,8] >= 0
 lBas: lDW[7,8] - ldw[7,8] >= 0
 lBas: lDX[2,9] - ldx[2,9] >= 0
 lBas: lDY[1,9] - ldy[1,9] >= 0
 lBas: lDZ[1,9] - ldz[1,9] >= 0
 lBas: lDW[1,9] - ldw[1,9] >= 0
 lBas: lDX[3,9] - ldx[3,9] >= 0
 lBas: lDY[2,9] - ldy[2,9] >= 0
 lBas: lDZ[2,9] - ldz[2,9] >= 0
 lBas: lDW[2,9] - ldw[2,9] >= 0
 lBas: lDX[4,9] - ldx[4,9] >= 0
 lBas: lDY[3,9] - ldy[3,9] >= 0
 lBas: lDZ[3,9] - ldz[3,9] >= 0
 lBas: lDW[3,9] - ldw[3,9] >= 0
 lBas: lDX[5,9] - ldx[5,9] >= 0
 lBas: lDY[4,9] - ldy[4,9] >= 0
 lBas: lDZ[4,9] - ldz[4,9] >= 0
 lBas: lDW[4,9] - ldw[4,9] >= 0
 lBas: lDX[6,9] - ldx[6,9] >= 0
 lBas: lDY[5,9] - ldy[5,9] >= 0
 lBas: lDZ[5,9] - ldz[5,9] >= 0
 lBas: lDW[5,9] - ldw[5,9] >= 0
 lBas: lDX[7,9] - ldx[7,9] >= 0
 lBas: lDY[6,9] - ldy[6,9] >= 0
 lBas: lDZ[6,9] - ldz[6,9] >= 0
 lBas: lDW[6,9] - ldw[6,9] >= 0
 lBas: lDX[8,9] - ldx[8,9] >= 0
 lBas: lDY[7,9] - ldy[7,9] >= 0
 lBas: lDZ[7,9] - ldz[7,9] >= 0
 lBas: lDW[7,9] - ldw[7,9] >= 0
 lBas: lDX[2,10] - ldx[2,10] >= 0
 lBas: lDY[1,10] - ldy[1,10] >= 0
 lBas: lDZ[1,10] - ldz[1,10] >= 0
 lBas: lDW[1,10] - ldw[1,10] >= 0
 lBas: lDX[3,10] - ldx[3,10] >= 0
 lBas: lDY[2,10] - ldy[2,10] >= 0
 lBas: lDZ[2,10] - ldz[2,10] >= 0
 lBas: lDW[2,10] - ldw[2,10] >= 0
 lBas: lDX[4,10] - ldx[4,10] >= 0
 lBas: lDY[3,10] - ldy[3,10] >= 0
 lBas: lDZ[3,10] - ldz[3,10] >= 0
 lBas: lDW[3,10] - ldw[3,10] >= 0
 lBas: lDX[5,10] - ldx[5,10] >= 0
 lBas: lDY[4,10] - ldy[4,10] >= 0
 lBas: lDZ[4,10] - ldz[4,10] >= 0
 lBas: lDW[4,10] - ldw[4,10] >= 0
 lBas: lDX[6,10] - ldx[6,10] >= 0
 lBas: lDY[5,10] - ldy[5,10] >= 0
 lBas: lDZ[5,10] - ldz[5,10] >= 0
 lBas: lDW[5,10] - ldw[5,10] >= 0
 lBas: lDX[7,10] - ldx[7,10] >= 0
 lBas: lDY[6,10] - ldy[6,10] >= 0
 lBas: lDZ[6,10] - ldz[6,10] >= 0
 lBas: lDW[6,10] - ldw[6,10] >= 0
 lBas: lDX[8,10] - ldx[8,10] >= 0
 lBas: lDY[7,10] - ldy[7,10] >= 0
 lBas: lDZ[7,10] - ldz[7,10] >= 0
 lBas: lDW[7,10] - ldw[7,10] >= 0
 lBas: lDX[2,11] - ldx[2,11] >= 0
 lBas: lDY[1,11] - ldy[1,11] >= 0
 lBas: lDZ[1,11] - ldz[1,11] >= 0
 lBas: lDW[1,11] - ldw[1,11] >= 0
 lBas: lDX[3,11] - ldx[3,11] >= 0
 lBas: lDY[2,11] - ldy[2,11] >= 0
 lBas: lDZ[2,11] - ldz[2,11] >= 0
 lBas: lDW[2,11] - ldw[2,11] >= 0
 lBas: lDX[4,11] - ldx[4,11] >= 0
 lBas: lDY[3,11] - ldy[3,11] >= 0
 lBas: lDZ[3,11] - ldz[3,11] >= 0
 lBas: lDW[3,11] - ldw[3,11] >= 0
 lBas: lDX[5,11] - ldx[5,11] >= 0
 lBas: lDY[4,11] - ldy[4,11] >= 0
 lBas: lDZ[4,11] - ldz[4,11] >= 0
 lBas: lDW[4,11] - ldw[4,11] >= 0
 lBas: lDX[6,11] - ldx[6,11] >= 0
 lBas: lDY[5,11] - ldy[5,11] >= 0
 lBas: lDZ[5,11] - ldz[5,11] >= 0
 lBas: lDW[5,11] - ldw[5,11] >= 0
 lBas: lDX[7,11] - ldx[7,11] >= 0
 lBas: lDY[6,11] - ldy[6,11] >= 0
 lBas: lDZ[6,11] - ldz[6,11] >= 0
 lBas: lDW[6,11] - ldw[6,11] >= 0
 lBas: lDX[8,11] - ldx[8,11] >= 0
 lBas: lDY[7,11] - ldy[7,11] >= 0
 lBas: lDZ[7,11] - ldz[7,11] >= 0
 lBas: lDW[7,11] - ldw[7,11] >= 0
 lBas: lDX[2,12] - ldx[2,12] >= 0
 lBas: lDY[1,12] - ldy[1,12] >= 0
 lBas: lDZ[1,12] - ldz[1,12] >= 0
 lBas: lDW[1,12] - ldw[1,12] >= 0
 lBas: lDX[3,12] - ldx[3,12] >= 0
 lBas: lDY[2,12] - ldy[2,12] >= 0
 lBas: lDZ[2,12] - ldz[2,12] >= 0
 lBas: lDW[2,12] - ldw[2,12] >= 0
 lBas: lDX[4,12] - ldx[4,12] >= 0
 lBas: lDY[3,12] - ldy[3,12] >= 0
 lBas: lDZ[3,12] - ldz[3,12] >= 0
 lBas: lDW[3,12] - ldw[3,12] >= 0
 lBas: lDX[5,12] - ldx[5,12] >= 0
 lBas: lDY[4,12] - ldy[4,12] >= 0
 lBas: lDZ[4,12] - ldz[4,12] >= 0
 lBas: lDW[4,12] - ldw[4,12] >= 0
 lBas: lDX[6,12] - ldx[6,12] >= 0
 lBas: lDY[5,12] - ldy[5,12] >= 0
 lBas: lDZ[5,12] - ldz[5,12] >= 0
 lBas: lDW[5,12] - ldw[5,12] >= 0
 lBas: lDX[7,12] - ldx[7,12] >= 0
 lBas: lDY[6,12] - ldy[6,12] >= 0
 lBas: lDZ[6,12] - ldz[6,12] >= 0
 lBas: lDW[6,12] - ldw[6,12] >= 0
 lBas: lDX[8,12] - ldx[8,12] >= 0
 lBas: lDY[7,12] - ldy[7,12] >= 0
 lBas: lDZ[7,12] - ldz[7,12] >= 0
 lBas: lDW[7,12] - ldw[7,12] >= 0
 lBas: lDX[2,13] - ldx[2,13] >= 0
 lBas: lDY[1,13] - ldy[1,13] >= 0
 lBas: lDZ[1,13] - ldz[1,13] >= 0
 lBas: lDW[1,13] - ldw[1,13] >= 0
 lBas: lDX[3,13] - ldx[3,13] >= 0
 lBas: lDY[2,13] - ldy[2,13] >= 0
 lBas: lDZ[2,13] - ldz[2,13] >= 0
 lBas: lDW[2,13] - ldw[2,13] >= 0
 lBas: lDX[4,13] - ldx[4,13] >= 0
 lBas: lDY[3,13] - ldy[3,13] >= 0
 lBas: lDZ[3,13] - ldz[3,13] >= 0
 lBas: lDW[3,13] - ldw[3,13] >= 0
 lBas: lDX[5,13] - ldx[5,13] >= 0
 lBas: lDY[4,13] - ldy[4,13] >= 0
 lBas: lDZ[4,13] - ldz[4,13] >= 0
 lBas: lDW[4,13] - ldw[4,13] >= 0
 lBas: lDX[6,13] - ldx[6,13] >= 0
 lBas: lDY[5,13] - ldy[5,13] >= 0
 lBas: lDZ[5,13] - ldz[5,13] >= 0
 lBas: lDW[5,13] - ldw[5,13] >= 0
 lBas: lDX[7,13] - ldx[7,13] >= 0
 lBas: lDY[6,13] - ldy[6,13] >= 0
 lBas: lDZ[6,13] - ldz[6,13] >= 0
 lBas: lDW[6,13] - ldw[6,13] >= 0
 lBas: lDX[8,13] - ldx[8,13] >= 0
 lBas: lDY[7,13] - ldy[7,13] >= 0
 lBas: lDZ[7,13] - ldz[7,13] >= 0
 lBas: lDW[7,13] - ldw[7,13] >= 0
 lBas: lDX[2,14] - ldx[2,14] >= 0
 lBas: lDY[1,14] - ldy[1,14] >= 0
 lBas: lDZ[1,14] - ldz[1,14] >= 0
 lBas: lDW[1,14] - ldw[1,14] >= 0
 lBas: lDX[3,14] - ldx[3,14] >= 0
 lBas: lDY[2,14] - ldy[2,14] >= 0
 lBas: lDZ[2,14] - ldz[2,14] >= 0
 lBas: lDW[2,14] - ldw[2,14] >= 0
 lBas: lDX[4,14] - ldx[4,14] >= 0
 lBas: lDY[3,14] - ldy[3,14] >= 0
 lBas: lDZ[3,14] - ldz[3,14] >= 0
 lBas: lDW[3,14] - ldw[3,14] >= 0
 lBas: lDX[5,14] - ldx[5,14] >= 0
 lBas: lDY[4,14] - ldy[4,14] >= 0
 lBas: lDZ[4,14] - ldz[4,14] >= 0
 lBas: lDW[4,14] - ldw[4,14] >= 0
 lBas: lDX[6,14] - ldx[6,14] >= 0
 lBas: lDY[5,14] - ldy[5,14] >= 0
 lBas: lDZ[5,14] - ldz[5,14] >= 0
 lBas: lDW[5,14] - ldw[5,14] >= 0
 lBas: lDX[7,14] - ldx[7,14] >= 0
 lBas: lDY[6,14] - ldy[6,14] >= 0
 lBas: lDZ[6,14] - ldz[6,14] >= 0
 lBas: lDW[6,14] - ldw[6,14] >= 0
 lBas: lDX[8,14] - ldx[8,14] >= 0
 lBas: lDY[7,14] - ldy[7,14] >= 0
 lBas: lDZ[7,14] - ldz[7,14] >= 0
 lBas: lDW[7,14] - ldw[7,14] >= 0
 lBas: lDX[2,15] - ldx[2,15] >= 0
 lBas: lDY[1,15] - ldy[1,15] >= 0
 lBas: lDZ[1,15] - ldz[1,15] >= 0
 lBas: lDW[1,15] - ldw[1,15] >= 0
 lBas: lDX[3,15] - ldx[3,15] >= 0
 lBas: lDY[2,15] - ldy[2,15] >= 0
 lBas: lDZ[2,15] - ldz[2,15] >= 0
 lBas: lDW[2,15] - ldw[2,15] >= 0
 lBas: lDX[4,15] - ldx[4,15] >= 0
 lBas: lDY[3,15] - ldy[3,15] >= 0
 lBas: lDZ[3,15] - ldz[3,15] >= 0
 lBas: lDW[3,15] - ldw[3,15] >= 0
 lBas: lDX[5,15] - ldx[5,15] >= 0
 lBas: lDY[4,15] - ldy[4,15] >= 0
 lBas: lDZ[4,15] - ldz[4,15] >= 0
 lBas: lDW[4,15] - ldw[4,15] >= 0
 lBas: lDX[6,15] - ldx[6,15] >= 0
 lBas: lDY[5,15] - ldy[5,15] >= 0
 lBas: lDZ[5,15] - ldz[5,15] >= 0
 lBas: lDW[5,15] - ldw[5,15] >= 0
 lBas: lDX[7,15] - ldx[7,15] >= 0
 lBas: lDY[6,15] - ldy[6,15] >= 0
 lBas: lDZ[6,15] - ldz[6,15] >= 0
 lBas: lDW[6,15] - ldw[6,15] >= 0
 lBas: lDX[8,15] - ldx[8,15] >= 0
 lBas: lDY[7,15] - ldy[7,15] >= 0
 lBas: lDZ[7,15] - ldz[7,15] >= 0
 lBas: lDW[7,15] - ldw[7,15] >= 0
 lATK: lDW[0,0] - ldx[1,0] - lcan[1,0] >= 0
 lATK: - lDW[0,0] + lDX[1,0] + lcan[1,0] >= 0
 lATK: lDW[0,0] - lDX[1,0] + lstk[1,0] - 2 lcan[1,0] >= 0
 lATK: lDX[1,0] - lstk[1,0] + lcan[1,0] >= 0
 lATK: ldw[0,0] - ldx[1,0] = 0
 lATK: lDW[0,1] - ldx[1,1] - lcan[1,1] >= 0
 lATK: - lDW[0,1] + lDX[1,1] + lcan[1,1] >= 0
 lATK: lDW[0,1] - lDX[1,1] + lstk[1,1] - 2 lcan[1,1] >= 0
 lATK: lDX[1,1] - lstk[1,1] + lcan[1,1] >= 0
 lATK: ldw[0,1] - ldx[1,1] = 0
 lATK: lDW[0,2] - ldx[1,2] - lcan[1,2] >= 0
 lATK: - lDW[0,2] + lDX[1,2] + lcan[1,2] >= 0
 lATK: lDW[0,2] - lDX[1,2] + lstk[1,2] - 2 lcan[1,2] >= 0
 lATK: lDX[1,2] - lstk[1,2] + lcan[1,2] >= 0
 lATK: ldw[0,2] - ldx[1,2] = 0
 lATK: lDW[0,3] - ldx[1,3] - lcan[1,3] >= 0
 lATK: - lDW[0,3] + lDX[1,3] + lcan[1,3] >= 0
 lATK: lDW[0,3] - lDX[1,3] + lstk[1,3] - 2 lcan[1,3] >= 0
 lATK: lDX[1,3] - lstk[1,3] + lcan[1,3] >= 0
 lATK: ldw[0,3] - ldx[1,3] = 0
 lATK: lDW[0,4] - ldx[1,4] - lcan[1,4] >= 0
 lATK: - lDW[0,4] + lDX[1,4] + lcan[1,4] >= 0
 lATK: lDW[0,4] - lDX[1,4] + lstk[1,4] - 2 lcan[1,4] >= 0
 lATK: lDX[1,4] - lstk[1,4] + lcan[1,4] >= 0
 lATK: ldw[0,4] - ldx[1,4] = 0
 lATK: lDW[0,5] - ldx[1,5] - lcan[1,5] >= 0
 lATK: - lDW[0,5] + lDX[1,5] + lcan[1,5] >= 0
 lATK: lDW[0,5] - lDX[1,5] + lstk[1,5] - 2 lcan[1,5] >= 0
 lATK: lDX[1,5] - lstk[1,5] + lcan[1,5] >= 0
 lATK: ldw[0,5] - ldx[1,5] = 0
 lATK: lDW[0,6] - ldx[1,6] - lcan[1,6] >= 0
 lATK: - lDW[0,6] + lDX[1,6] + lcan[1,6] >= 0
 lATK: lDW[0,6] - lDX[1,6] + lstk[1,6] - 2 lcan[1,6] >= 0
 lATK: lDX[1,6] - lstk[1,6] + lcan[1,6] >= 0
 lATK: ldw[0,6] - ldx[1,6] = 0
 lATK: lDW[0,7] - ldx[1,7] - lcan[1,7] >= 0
 lATK: - lDW[0,7] + lDX[1,7] + lcan[1,7] >= 0
 lATK: lDW[0,7] - lDX[1,7] + lstk[1,7] - 2 lcan[1,7] >= 0
 lATK: lDX[1,7] - lstk[1,7] + lcan[1,7] >= 0
 lATK: ldw[0,7] - ldx[1,7] = 0
 lATK: lDW[0,8] - ldx[1,8] - lcan[1,8] >= 0
 lATK: - lDW[0,8] + lDX[1,8] + lcan[1,8] >= 0
 lATK: lDW[0,8] - lDX[1,8] + lstk[1,8] - 2 lcan[1,8] >= 0
 lATK: lDX[1,8] - lstk[1,8] + lcan[1,8] >= 0
 lATK: ldw[0,8] - ldx[1,8] = 0
 lATK: lDW[0,9] - ldx[1,9] - lcan[1,9] >= 0
 lATK: - lDW[0,9] + lDX[1,9] + lcan[1,9] >= 0
 lATK: lDW[0,9] - lDX[1,9] + lstk[1,9] - 2 lcan[1,9] >= 0
 lATK: lDX[1,9] - lstk[1,9] + lcan[1,9] >= 0
 lATK: ldw[0,9] - ldx[1,9] = 0
 lATK: lDW[0,10] - ldx[1,10] - lcan[1,10] >= 0
 lATK: - lDW[0,10] + lDX[1,10] + lcan[1,10] >= 0
 lATK: lDW[0,10] - lDX[1,10] + lstk[1,10] - 2 lcan[1,10] >= 0
 lATK: lDX[1,10] - lstk[1,10] + lcan[1,10] >= 0
 lATK: ldw[0,10] - ldx[1,10] = 0
 lATK: lDW[0,11] - ldx[1,11] - lcan[1,11] >= 0
 lATK: - lDW[0,11] + lDX[1,11] + lcan[1,11] >= 0
 lATK: lDW[0,11] - lDX[1,11] + lstk[1,11] - 2 lcan[1,11] >= 0
 lATK: lDX[1,11] - lstk[1,11] + lcan[1,11] >= 0
 lATK: ldw[0,11] - ldx[1,11] = 0
 lATK: lDW[0,12] - ldx[1,12] - lcan[1,12] >= 0
 lATK: - lDW[0,12] + lDX[1,12] + lcan[1,12] >= 0
 lATK: lDW[0,12] - lDX[1,12] + lstk[1,12] - 2 lcan[1,12] >= 0
 lATK: lDX[1,12] - lstk[1,12] + lcan[1,12] >= 0
 lATK: ldw[0,12] - ldx[1,12] = 0
 lATK: lDW[0,13] - ldx[1,13] - lcan[1,13] >= 0
 lATK: - lDW[0,13] + lDX[1,13] + lcan[1,13] >= 0
 lATK: lDW[0,13] - lDX[1,13] + lstk[1,13] - 2 lcan[1,13] >= 0
 lATK: lDX[1,13] - lstk[1,13] + lcan[1,13] >= 0
 lATK: ldw[0,13] - ldx[1,13] = 0
 lATK: lDW[0,14] - ldx[1,14] - lcan[1,14] >= 0
 lATK: - lDW[0,14] + lDX[1,14] + lcan[1,14] >= 0
 lATK: lDW[0,14] - lDX[1,14] + lstk[1,14] - 2 lcan[1,14] >= 0
 lATK: lDX[1,14] - lstk[1,14] + lcan[1,14] >= 0
 lATK: ldw[0,14] - ldx[1,14] = 0
 lATK: lDW[0,15] - ldx[1,15] - lcan[1,15] >= 0
 lATK: - lDW[0,15] + lDX[1,15] + lcan[1,15] >= 0
 lATK: lDW[0,15] - lDX[1,15] + lstk[1,15] - 2 lcan[1,15] >= 0
 lATK: lDX[1,15] - lstk[1,15] + lcan[1,15] >= 0
 lATK: ldw[0,15] - ldx[1,15] = 0
 lATK: lDW[1,0] - ldx[2,0] - lcan[2,0] >= 0
 lATK: - lDW[1,0] + lDX[2,0] + lcan[2,0] >= 0
 lATK: lDW[1,0] - lDX[2,0] + lstk[2,0] - 2 lcan[2,0] >= 0
 lATK: lDX[2,0] - lstk[2,0] + lcan[2,0] >= 0
 lATK: ldw[1,0] - ldx[2,0] = 0
 lATK: lDW[1,1] - ldx[2,1] - lcan[2,1] >= 0
 lATK: - lDW[1,1] + lDX[2,1] + lcan[2,1] >= 0
 lATK: lDW[1,1] - lDX[2,1] + lstk[2,1] - 2 lcan[2,1] >= 0
 lATK: lDX[2,1] - lstk[2,1] + lcan[2,1] >= 0
 lATK: ldw[1,1] - ldx[2,1] = 0
 lATK: lDW[1,2] - ldx[2,2] - lcan[2,2] >= 0
 lATK: - lDW[1,2] + lDX[2,2] + lcan[2,2] >= 0
 lATK: lDW[1,2] - lDX[2,2] + lstk[2,2] - 2 lcan[2,2] >= 0
 lATK: lDX[2,2] - lstk[2,2] + lcan[2,2] >= 0
 lATK: ldw[1,2] - ldx[2,2] = 0
 lATK: lDW[1,3] - ldx[2,3] - lcan[2,3] >= 0
 lATK: - lDW[1,3] + lDX[2,3] + lcan[2,3] >= 0
 lATK: lDW[1,3] - lDX[2,3] + lstk[2,3] - 2 lcan[2,3] >= 0
 lATK: lDX[2,3] - lstk[2,3] + lcan[2,3] >= 0
 lATK: ldw[1,3] - ldx[2,3] = 0
 lATK: lDW[1,4] - ldx[2,4] - lcan[2,4] >= 0
 lATK: - lDW[1,4] + lDX[2,4] + lcan[2,4] >= 0
 lATK: lDW[1,4] - lDX[2,4] + lstk[2,4] - 2 lcan[2,4] >= 0
 lATK: lDX[2,4] - lstk[2,4] + lcan[2,4] >= 0
 lATK: ldw[1,4] - ldx[2,4] = 0
 lATK: lDW[1,5] - ldx[2,5] - lcan[2,5] >= 0
 lATK: - lDW[1,5] + lDX[2,5] + lcan[2,5] >= 0
 lATK: lDW[1,5] - lDX[2,5] + lstk[2,5] - 2 lcan[2,5] >= 0
 lATK: lDX[2,5] - lstk[2,5] + lcan[2,5] >= 0
 lATK: ldw[1,5] - ldx[2,5] = 0
 lATK: lDW[1,6] - ldx[2,6] - lcan[2,6] >= 0
 lATK: - lDW[1,6] + lDX[2,6] + lcan[2,6] >= 0
 lATK: lDW[1,6] - lDX[2,6] + lstk[2,6] - 2 lcan[2,6] >= 0
 lATK: lDX[2,6] - lstk[2,6] + lcan[2,6] >= 0
 lATK: ldw[1,6] - ldx[2,6] = 0
 lATK: lDW[1,7] - ldx[2,7] - lcan[2,7] >= 0
 lATK: - lDW[1,7] + lDX[2,7] + lcan[2,7] >= 0
 lATK: lDW[1,7] - lDX[2,7] + lstk[2,7] - 2 lcan[2,7] >= 0
 lATK: lDX[2,7] - lstk[2,7] + lcan[2,7] >= 0
 lATK: ldw[1,7] - ldx[2,7] = 0
 lATK: lDW[1,8] - ldx[2,8] - lcan[2,8] >= 0
 lATK: - lDW[1,8] + lDX[2,8] + lcan[2,8] >= 0
 lATK: lDW[1,8] - lDX[2,8] + lstk[2,8] - 2 lcan[2,8] >= 0
 lATK: lDX[2,8] - lstk[2,8] + lcan[2,8] >= 0
 lATK: ldw[1,8] - ldx[2,8] = 0
 lATK: lDW[1,9] - ldx[2,9] - lcan[2,9] >= 0
 lATK: - lDW[1,9] + lDX[2,9] + lcan[2,9] >= 0
 lATK: lDW[1,9] - lDX[2,9] + lstk[2,9] - 2 lcan[2,9] >= 0
 lATK: lDX[2,9] - lstk[2,9] + lcan[2,9] >= 0
 lATK: ldw[1,9] - ldx[2,9] = 0
 lATK: lDW[1,10] - ldx[2,10] - lcan[2,10] >= 0
 lATK: - lDW[1,10] + lDX[2,10] + lcan[2,10] >= 0
 lATK: lDW[1,10] - lDX[2,10] + lstk[2,10] - 2 lcan[2,10] >= 0
 lATK: lDX[2,10] - lstk[2,10] + lcan[2,10] >= 0
 lATK: ldw[1,10] - ldx[2,10] = 0
 lATK: lDW[1,11] - ldx[2,11] - lcan[2,11] >= 0
 lATK: - lDW[1,11] + lDX[2,11] + lcan[2,11] >= 0
 lATK: lDW[1,11] - lDX[2,11] + lstk[2,11] - 2 lcan[2,11] >= 0
 lATK: lDX[2,11] - lstk[2,11] + lcan[2,11] >= 0
 lATK: ldw[1,11] - ldx[2,11] = 0
 lATK: lDW[1,12] - ldx[2,12] - lcan[2,12] >= 0
 lATK: - lDW[1,12] + lDX[2,12] + lcan[2,12] >= 0
 lATK: lDW[1,12] - lDX[2,12] + lstk[2,12] - 2 lcan[2,12] >= 0
 lATK: lDX[2,12] - lstk[2,12] + lcan[2,12] >= 0
 lATK: ldw[1,12] - ldx[2,12] = 0
 lATK: lDW[1,13] - ldx[2,13] - lcan[2,13] >= 0
 lATK: - lDW[1,13] + lDX[2,13] + lcan[2,13] >= 0
 lATK: lDW[1,13] - lDX[2,13] + lstk[2,13] - 2 lcan[2,13] >= 0
 lATK: lDX[2,13] - lstk[2,13] + lcan[2,13] >= 0
 lATK: ldw[1,13] - ldx[2,13] = 0
 lATK: lDW[1,14] - ldx[2,14] - lcan[2,14] >= 0
 lATK: - lDW[1,14] + lDX[2,14] + lcan[2,14] >= 0
 lATK: lDW[1,14] - lDX[2,14] + lstk[2,14] - 2 lcan[2,14] >= 0
 lATK: lDX[2,14] - lstk[2,14] + lcan[2,14] >= 0
 lATK: ldw[1,14] - ldx[2,14] = 0
 lATK: lDW[1,15] - ldx[2,15] - lcan[2,15] >= 0
 lATK: - lDW[1,15] + lDX[2,15] + lcan[2,15] >= 0
 lATK: lDW[1,15] - lDX[2,15] + lstk[2,15] - 2 lcan[2,15] >= 0
 lATK: lDX[2,15] - lstk[2,15] + lcan[2,15] >= 0
 lATK: ldw[1,15] - ldx[2,15] = 0
 lATK: lDW[2,0] - ldx[3,0] - lcan[3,0] >= 0
 lATK: - lDW[2,0] + lDX[3,0] + lcan[3,0] >= 0
 lATK: lDW[2,0] - lDX[3,0] + lstk[3,0] - 2 lcan[3,0] >= 0
 lATK: lDX[3,0] - lstk[3,0] + lcan[3,0] >= 0
 lATK: ldw[2,0] - ldx[3,0] = 0
 lATK: lDW[2,1] - ldx[3,1] - lcan[3,1] >= 0
 lATK: - lDW[2,1] + lDX[3,1] + lcan[3,1] >= 0
 lATK: lDW[2,1] - lDX[3,1] + lstk[3,1] - 2 lcan[3,1] >= 0
 lATK: lDX[3,1] - lstk[3,1] + lcan[3,1] >= 0
 lATK: ldw[2,1] - ldx[3,1] = 0
 lATK: lDW[2,2] - ldx[3,2] - lcan[3,2] >= 0
 lATK: - lDW[2,2] + lDX[3,2] + lcan[3,2] >= 0
 lATK: lDW[2,2] - lDX[3,2] + lstk[3,2] - 2 lcan[3,2] >= 0
 lATK: lDX[3,2] - lstk[3,2] + lcan[3,2] >= 0
 lATK: ldw[2,2] - ldx[3,2] = 0
 lATK: lDW[2,3] - ldx[3,3] - lcan[3,3] >= 0
 lATK: - lDW[2,3] + lDX[3,3] + lcan[3,3] >= 0
 lATK: lDW[2,3] - lDX[3,3] + lstk[3,3] - 2 lcan[3,3] >= 0
 lATK: lDX[3,3] - lstk[3,3] + lcan[3,3] >= 0
 lATK: ldw[2,3] - ldx[3,3] = 0
 lATK: lDW[2,4] - ldx[3,4] - lcan[3,4] >= 0
 lATK: - lDW[2,4] + lDX[3,4] + lcan[3,4] >= 0
 lATK: lDW[2,4] - lDX[3,4] + lstk[3,4] - 2 lcan[3,4] >= 0
 lATK: lDX[3,4] - lstk[3,4] + lcan[3,4] >= 0
 lATK: ldw[2,4] - ldx[3,4] = 0
 lATK: lDW[2,5] - ldx[3,5] - lcan[3,5] >= 0
 lATK: - lDW[2,5] + lDX[3,5] + lcan[3,5] >= 0
 lATK: lDW[2,5] - lDX[3,5] + lstk[3,5] - 2 lcan[3,5] >= 0
 lATK: lDX[3,5] - lstk[3,5] + lcan[3,5] >= 0
 lATK: ldw[2,5] - ldx[3,5] = 0
 lATK: lDW[2,6] - ldx[3,6] - lcan[3,6] >= 0
 lATK: - lDW[2,6] + lDX[3,6] + lcan[3,6] >= 0
 lATK: lDW[2,6] - lDX[3,6] + lstk[3,6] - 2 lcan[3,6] >= 0
 lATK: lDX[3,6] - lstk[3,6] + lcan[3,6] >= 0
 lATK: ldw[2,6] - ldx[3,6] = 0
 lATK: lDW[2,7] - ldx[3,7] - lcan[3,7] >= 0
 lATK: - lDW[2,7] + lDX[3,7] + lcan[3,7] >= 0
 lATK: lDW[2,7] - lDX[3,7] + lstk[3,7] - 2 lcan[3,7] >= 0
 lATK: lDX[3,7] - lstk[3,7] + lcan[3,7] >= 0
 lATK: ldw[2,7] - ldx[3,7] = 0
 lATK: lDW[2,8] - ldx[3,8] - lcan[3,8] >= 0
 lATK: - lDW[2,8] + lDX[3,8] + lcan[3,8] >= 0
 lATK: lDW[2,8] - lDX[3,8] + lstk[3,8] - 2 lcan[3,8] >= 0
 lATK: lDX[3,8] - lstk[3,8] + lcan[3,8] >= 0
 lATK: ldw[2,8] - ldx[3,8] = 0
 lATK: lDW[2,9] - ldx[3,9] - lcan[3,9] >= 0
 lATK: - lDW[2,9] + lDX[3,9] + lcan[3,9] >= 0
 lATK: lDW[2,9] - lDX[3,9] + lstk[3,9] - 2 lcan[3,9] >= 0
 lATK: lDX[3,9] - lstk[3,9] + lcan[3,9] >= 0
 lATK: ldw[2,9] - ldx[3,9] = 0
 lATK: lDW[2,10] - ldx[3,10] - lcan[3,10] >= 0
 lATK: - lDW[2,10] + lDX[3,10] + lcan[3,10] >= 0
 lATK: lDW[2,10] - lDX[3,10] + lstk[3,10] - 2 lcan[3,10] >= 0
 lATK: lDX[3,10] - lstk[3,10] + lcan[3,10] >= 0
 lATK: ldw[2,10] - ldx[3,10] = 0
 lATK: lDW[2,11] - ldx[3,11] - lcan[3,11] >= 0
 lATK: - lDW[2,11] + lDX[3,11] + lcan[3,11] >= 0
 lATK: lDW[2,11] - lDX[3,11] + lstk[3,11] - 2 lcan[3,11] >= 0
 lATK: lDX[3,11] - lstk[3,11] + lcan[3,11] >= 0
 lATK: ldw[2,11] - ldx[3,11] = 0
 lATK: lDW[2,12] - ldx[3,12] - lcan[3,12] >= 0
 lATK: - lDW[2,12] + lDX[3,12] + lcan[3,12] >= 0
 lATK: lDW[2,12] - lDX[3,12] + lstk[3,12] - 2 lcan[3,12] >= 0
 lATK: lDX[3,12] - lstk[3,12] + lcan[3,12] >= 0
 lATK: ldw[2,12] - ldx[3,12] = 0
 lATK: lDW[2,13] - ldx[3,13] - lcan[3,13] >= 0
 lATK: - lDW[2,13] + lDX[3,13] + lcan[3,13] >= 0
 lATK: lDW[2,13] - lDX[3,13] + lstk[3,13] - 2 lcan[3,13] >= 0
 lATK: lDX[3,13] - lstk[3,13] + lcan[3,13] >= 0
 lATK: ldw[2,13] - ldx[3,13] = 0
 lATK: lDW[2,14] - ldx[3,14] - lcan[3,14] >= 0
 lATK: - lDW[2,14] + lDX[3,14] + lcan[3,14] >= 0
 lATK: lDW[2,14] - lDX[3,14] + lstk[3,14] - 2 lcan[3,14] >= 0
 lATK: lDX[3,14] - lstk[3,14] + lcan[3,14] >= 0
 lATK: ldw[2,14] - ldx[3,14] = 0
 lATK: lDW[2,15] - ldx[3,15] - lcan[3,15] >= 0
 lATK: - lDW[2,15] + lDX[3,15] + lcan[3,15] >= 0
 lATK: lDW[2,15] - lDX[3,15] + lstk[3,15] - 2 lcan[3,15] >= 0
 lATK: lDX[3,15] - lstk[3,15] + lcan[3,15] >= 0
 lATK: ldw[2,15] - ldx[3,15] = 0
 lATK: lDW[3,0] - ldx[4,0] - lcan[4,0] >= 0
 lATK: - lDW[3,0] + lDX[4,0] + lcan[4,0] >= 0
 lATK: lDW[3,0] - lDX[4,0] + lstk[4,0] - 2 lcan[4,0] >= 0
 lATK: lDX[4,0] - lstk[4,0] + lcan[4,0] >= 0
 lATK: ldw[3,0] - ldx[4,0] = 0
 lATK: lDW[3,1] - ldx[4,1] - lcan[4,1] >= 0
 lATK: - lDW[3,1] + lDX[4,1] + lcan[4,1] >= 0
 lATK: lDW[3,1] - lDX[4,1] + lstk[4,1] - 2 lcan[4,1] >= 0
 lATK: lDX[4,1] - lstk[4,1] + lcan[4,1] >= 0
 lATK: ldw[3,1] - ldx[4,1] = 0
 lATK: lDW[3,2] - ldx[4,2] - lcan[4,2] >= 0
 lATK: - lDW[3,2] + lDX[4,2] + lcan[4,2] >= 0
 lATK: lDW[3,2] - lDX[4,2] + lstk[4,2] - 2 lcan[4,2] >= 0
 lATK: lDX[4,2] - lstk[4,2] + lcan[4,2] >= 0
 lATK: ldw[3,2] - ldx[4,2] = 0
 lATK: lDW[3,3] - ldx[4,3] - lcan[4,3] >= 0
 lATK: - lDW[3,3] + lDX[4,3] + lcan[4,3] >= 0
 lATK: lDW[3,3] - lDX[4,3] + lstk[4,3] - 2 lcan[4,3] >= 0
 lATK: lDX[4,3] - lstk[4,3] + lcan[4,3] >= 0
 lATK: ldw[3,3] - ldx[4,3] = 0
 lATK: lDW[3,4] - ldx[4,4] - lcan[4,4] >= 0
 lATK: - lDW[3,4] + lDX[4,4] + lcan[4,4] >= 0
 lATK: lDW[3,4] - lDX[4,4] + lstk[4,4] - 2 lcan[4,4] >= 0
 lATK: lDX[4,4] - lstk[4,4] + lcan[4,4] >= 0
 lATK: ldw[3,4] - ldx[4,4] = 0
 lATK: lDW[3,5] - ldx[4,5] - lcan[4,5] >= 0
 lATK: - lDW[3,5] + lDX[4,5] + lcan[4,5] >= 0
 lATK: lDW[3,5] - lDX[4,5] + lstk[4,5] - 2 lcan[4,5] >= 0
 lATK: lDX[4,5] - lstk[4,5] + lcan[4,5] >= 0
 lATK: ldw[3,5] - ldx[4,5] = 0
 lATK: lDW[3,6] - ldx[4,6] - lcan[4,6] >= 0
 lATK: - lDW[3,6] + lDX[4,6] + lcan[4,6] >= 0
 lATK: lDW[3,6] - lDX[4,6] + lstk[4,6] - 2 lcan[4,6] >= 0
 lATK: lDX[4,6] - lstk[4,6] + lcan[4,6] >= 0
 lATK: ldw[3,6] - ldx[4,6] = 0
 lATK: lDW[3,7] - ldx[4,7] - lcan[4,7] >= 0
 lATK: - lDW[3,7] + lDX[4,7] + lcan[4,7] >= 0
 lATK: lDW[3,7] - lDX[4,7] + lstk[4,7] - 2 lcan[4,7] >= 0
 lATK: lDX[4,7] - lstk[4,7] + lcan[4,7] >= 0
 lATK: ldw[3,7] - ldx[4,7] = 0
 lATK: lDW[3,8] - ldx[4,8] - lcan[4,8] >= 0
 lATK: - lDW[3,8] + lDX[4,8] + lcan[4,8] >= 0
 lATK: lDW[3,8] - lDX[4,8] + lstk[4,8] - 2 lcan[4,8] >= 0
 lATK: lDX[4,8] - lstk[4,8] + lcan[4,8] >= 0
 lATK: ldw[3,8] - ldx[4,8] = 0
 lATK: lDW[3,9] - ldx[4,9] - lcan[4,9] >= 0
 lATK: - lDW[3,9] + lDX[4,9] + lcan[4,9] >= 0
 lATK: lDW[3,9] - lDX[4,9] + lstk[4,9] - 2 lcan[4,9] >= 0
 lATK: lDX[4,9] - lstk[4,9] + lcan[4,9] >= 0
 lATK: ldw[3,9] - ldx[4,9] = 0
 lATK: lDW[3,10] - ldx[4,10] - lcan[4,10] >= 0
 lATK: - lDW[3,10] + lDX[4,10] + lcan[4,10] >= 0
 lATK: lDW[3,10] - lDX[4,10] + lstk[4,10] - 2 lcan[4,10] >= 0
 lATK: lDX[4,10] - lstk[4,10] + lcan[4,10] >= 0
 lATK: ldw[3,10] - ldx[4,10] = 0
 lATK: lDW[3,11] - ldx[4,11] - lcan[4,11] >= 0
 lATK: - lDW[3,11] + lDX[4,11] + lcan[4,11] >= 0
 lATK: lDW[3,11] - lDX[4,11] + lstk[4,11] - 2 lcan[4,11] >= 0
 lATK: lDX[4,11] - lstk[4,11] + lcan[4,11] >= 0
 lATK: ldw[3,11] - ldx[4,11] = 0
 lATK: lDW[3,12] - ldx[4,12] - lcan[4,12] >= 0
 lATK: - lDW[3,12] + lDX[4,12] + lcan[4,12] >= 0
 lATK: lDW[3,12] - lDX[4,12] + lstk[4,12] - 2 lcan[4,12] >= 0
 lATK: lDX[4,12] - lstk[4,12] + lcan[4,12] >= 0
 lATK: ldw[3,12] - ldx[4,12] = 0
 lATK: lDW[3,13] - ldx[4,13] - lcan[4,13] >= 0
 lATK: - lDW[3,13] + lDX[4,13] + lcan[4,13] >= 0
 lATK: lDW[3,13] - lDX[4,13] + lstk[4,13] - 2 lcan[4,13] >= 0
 lATK: lDX[4,13] - lstk[4,13] + lcan[4,13] >= 0
 lATK: ldw[3,13] - ldx[4,13] = 0
 lATK: lDW[3,14] - ldx[4,14] - lcan[4,14] >= 0
 lATK: - lDW[3,14] + lDX[4,14] + lcan[4,14] >= 0
 lATK: lDW[3,14] - lDX[4,14] + lstk[4,14] - 2 lcan[4,14] >= 0
 lATK: lDX[4,14] - lstk[4,14] + lcan[4,14] >= 0
 lATK: ldw[3,14] - ldx[4,14] = 0
 lATK: lDW[3,15] - ldx[4,15] - lcan[4,15] >= 0
 lATK: - lDW[3,15] + lDX[4,15] + lcan[4,15] >= 0
 lATK: lDW[3,15] - lDX[4,15] + lstk[4,15] - 2 lcan[4,15] >= 0
 lATK: lDX[4,15] - lstk[4,15] + lcan[4,15] >= 0
 lATK: ldw[3,15] - ldx[4,15] = 0
 lATK: lDW[4,0] - ldx[5,0] - lcan[5,0] >= 0
 lATK: - lDW[4,0] + lDX[5,0] + lcan[5,0] >= 0
 lATK: lDW[4,0] - lDX[5,0] + lstk[5,0] - 2 lcan[5,0] >= 0
 lATK: lDX[5,0] - lstk[5,0] + lcan[5,0] >= 0
 lATK: ldw[4,0] - ldx[5,0] = 0
 lATK: lDW[4,1] - ldx[5,1] - lcan[5,1] >= 0
 lATK: - lDW[4,1] + lDX[5,1] + lcan[5,1] >= 0
 lATK: lDW[4,1] - lDX[5,1] + lstk[5,1] - 2 lcan[5,1] >= 0
 lATK: lDX[5,1] - lstk[5,1] + lcan[5,1] >= 0
 lATK: ldw[4,1] - ldx[5,1] = 0
 lATK: lDW[4,2] - ldx[5,2] - lcan[5,2] >= 0
 lATK: - lDW[4,2] + lDX[5,2] + lcan[5,2] >= 0
 lATK: lDW[4,2] - lDX[5,2] + lstk[5,2] - 2 lcan[5,2] >= 0
 lATK: lDX[5,2] - lstk[5,2] + lcan[5,2] >= 0
 lATK: ldw[4,2] - ldx[5,2] = 0
 lATK: lDW[4,3] - ldx[5,3] - lcan[5,3] >= 0
 lATK: - lDW[4,3] + lDX[5,3] + lcan[5,3] >= 0
 lATK: lDW[4,3] - lDX[5,3] + lstk[5,3] - 2 lcan[5,3] >= 0
 lATK: lDX[5,3] - lstk[5,3] + lcan[5,3] >= 0
 lATK: ldw[4,3] - ldx[5,3] = 0
 lATK: lDW[4,4] - ldx[5,4] - lcan[5,4] >= 0
 lATK: - lDW[4,4] + lDX[5,4] + lcan[5,4] >= 0
 lATK: lDW[4,4] - lDX[5,4] + lstk[5,4] - 2 lcan[5,4] >= 0
 lATK: lDX[5,4] - lstk[5,4] + lcan[5,4] >= 0
 lATK: ldw[4,4] - ldx[5,4] = 0
 lATK: lDW[4,5] - ldx[5,5] - lcan[5,5] >= 0
 lATK: - lDW[4,5] + lDX[5,5] + lcan[5,5] >= 0
 lATK: lDW[4,5] - lDX[5,5] + lstk[5,5] - 2 lcan[5,5] >= 0
 lATK: lDX[5,5] - lstk[5,5] + lcan[5,5] >= 0
 lATK: ldw[4,5] - ldx[5,5] = 0
 lATK: lDW[4,6] - ldx[5,6] - lcan[5,6] >= 0
 lATK: - lDW[4,6] + lDX[5,6] + lcan[5,6] >= 0
 lATK: lDW[4,6] - lDX[5,6] + lstk[5,6] - 2 lcan[5,6] >= 0
 lATK: lDX[5,6] - lstk[5,6] + lcan[5,6] >= 0
 lATK: ldw[4,6] - ldx[5,6] = 0
 lATK: lDW[4,7] - ldx[5,7] - lcan[5,7] >= 0
 lATK: - lDW[4,7] + lDX[5,7] + lcan[5,7] >= 0
 lATK: lDW[4,7] - lDX[5,7] + lstk[5,7] - 2 lcan[5,7] >= 0
 lATK: lDX[5,7] - lstk[5,7] + lcan[5,7] >= 0
 lATK: ldw[4,7] - ldx[5,7] = 0
 lATK: lDW[4,8] - ldx[5,8] - lcan[5,8] >= 0
 lATK: - lDW[4,8] + lDX[5,8] + lcan[5,8] >= 0
 lATK: lDW[4,8] - lDX[5,8] + lstk[5,8] - 2 lcan[5,8] >= 0
 lATK: lDX[5,8] - lstk[5,8] + lcan[5,8] >= 0
 lATK: ldw[4,8] - ldx[5,8] = 0
 lATK: lDW[4,9] - ldx[5,9] - lcan[5,9] >= 0
 lATK: - lDW[4,9] + lDX[5,9] + lcan[5,9] >= 0
 lATK: lDW[4,9] - lDX[5,9] + lstk[5,9] - 2 lcan[5,9] >= 0
 lATK: lDX[5,9] - lstk[5,9] + lcan[5,9] >= 0
 lATK: ldw[4,9] - ldx[5,9] = 0
 lATK: lDW[4,10] - ldx[5,10] - lcan[5,10] >= 0
 lATK: - lDW[4,10] + lDX[5,10] + lcan[5,10] >= 0
 lATK: lDW[4,10] - lDX[5,10] + lstk[5,10] - 2 lcan[5,10] >= 0
 lATK: lDX[5,10] - lstk[5,10] + lcan[5,10] >= 0
 lATK: ldw[4,10] - ldx[5,10] = 0
 lATK: lDW[4,11] - ldx[5,11] - lcan[5,11] >= 0
 lATK: - lDW[4,11] + lDX[5,11] + lcan[5,11] >= 0
 lATK: lDW[4,11] - lDX[5,11] + lstk[5,11] - 2 lcan[5,11] >= 0
 lATK: lDX[5,11] - lstk[5,11] + lcan[5,11] >= 0
 lATK: ldw[4,11] - ldx[5,11] = 0
 lATK: lDW[4,12] - ldx[5,12] - lcan[5,12] >= 0
 lATK: - lDW[4,12] + lDX[5,12] + lcan[5,12] >= 0
 lATK: lDW[4,12] - lDX[5,12] + lstk[5,12] - 2 lcan[5,12] >= 0
 lATK: lDX[5,12] - lstk[5,12] + lcan[5,12] >= 0
 lATK: ldw[4,12] - ldx[5,12] = 0
 lATK: lDW[4,13] - ldx[5,13] - lcan[5,13] >= 0
 lATK: - lDW[4,13] + lDX[5,13] + lcan[5,13] >= 0
 lATK: lDW[4,13] - lDX[5,13] + lstk[5,13] - 2 lcan[5,13] >= 0
 lATK: lDX[5,13] - lstk[5,13] + lcan[5,13] >= 0
 lATK: ldw[4,13] - ldx[5,13] = 0
 lATK: lDW[4,14] - ldx[5,14] - lcan[5,14] >= 0
 lATK: - lDW[4,14] + lDX[5,14] + lcan[5,14] >= 0
 lATK: lDW[4,14] - lDX[5,14] + lstk[5,14] - 2 lcan[5,14] >= 0
 lATK: lDX[5,14] - lstk[5,14] + lcan[5,14] >= 0
 lATK: ldw[4,14] - ldx[5,14] = 0
 lATK: lDW[4,15] - ldx[5,15] - lcan[5,15] >= 0
 lATK: - lDW[4,15] + lDX[5,15] + lcan[5,15] >= 0
 lATK: lDW[4,15] - lDX[5,15] + lstk[5,15] - 2 lcan[5,15] >= 0
 lATK: lDX[5,15] - lstk[5,15] + lcan[5,15] >= 0
 lATK: ldw[4,15] - ldx[5,15] = 0
 lATK: lDW[5,0] - ldx[6,0] - lcan[6,0] >= 0
 lATK: - lDW[5,0] + lDX[6,0] + lcan[6,0] >= 0
 lATK: lDW[5,0] - lDX[6,0] + lstk[6,0] - 2 lcan[6,0] >= 0
 lATK: lDX[6,0] - lstk[6,0] + lcan[6,0] >= 0
 lATK: ldw[5,0] - ldx[6,0] = 0
 lATK: lDW[5,1] - ldx[6,1] - lcan[6,1] >= 0
 lATK: - lDW[5,1] + lDX[6,1] + lcan[6,1] >= 0
 lATK: lDW[5,1] - lDX[6,1] + lstk[6,1] - 2 lcan[6,1] >= 0
 lATK: lDX[6,1] - lstk[6,1] + lcan[6,1] >= 0
 lATK: ldw[5,1] - ldx[6,1] = 0
 lATK: lDW[5,2] - ldx[6,2] - lcan[6,2] >= 0
 lATK: - lDW[5,2] + lDX[6,2] + lcan[6,2] >= 0
 lATK: lDW[5,2] - lDX[6,2] + lstk[6,2] - 2 lcan[6,2] >= 0
 lATK: lDX[6,2] - lstk[6,2] + lcan[6,2] >= 0
 lATK: ldw[5,2] - ldx[6,2] = 0
 lATK: lDW[5,3] - ldx[6,3] - lcan[6,3] >= 0
 lATK: - lDW[5,3] + lDX[6,3] + lcan[6,3] >= 0
 lATK: lDW[5,3] - lDX[6,3] + lstk[6,3] - 2 lcan[6,3] >= 0
 lATK: lDX[6,3] - lstk[6,3] + lcan[6,3] >= 0
 lATK: ldw[5,3] - ldx[6,3] = 0
 lATK: lDW[5,4] - ldx[6,4] - lcan[6,4] >= 0
 lATK: - lDW[5,4] + lDX[6,4] + lcan[6,4] >= 0
 lATK: lDW[5,4] - lDX[6,4] + lstk[6,4] - 2 lcan[6,4] >= 0
 lATK: lDX[6,4] - lstk[6,4] + lcan[6,4] >= 0
 lATK: ldw[5,4] - ldx[6,4] = 0
 lATK: lDW[5,5] - ldx[6,5] - lcan[6,5] >= 0
 lATK: - lDW[5,5] + lDX[6,5] + lcan[6,5] >= 0
 lATK: lDW[5,5] - lDX[6,5] + lstk[6,5] - 2 lcan[6,5] >= 0
 lATK: lDX[6,5] - lstk[6,5] + lcan[6,5] >= 0
 lATK: ldw[5,5] - ldx[6,5] = 0
 lATK: lDW[5,6] - ldx[6,6] - lcan[6,6] >= 0
 lATK: - lDW[5,6] + lDX[6,6] + lcan[6,6] >= 0
 lATK: lDW[5,6] - lDX[6,6] + lstk[6,6] - 2 lcan[6,6] >= 0
 lATK: lDX[6,6] - lstk[6,6] + lcan[6,6] >= 0
 lATK: ldw[5,6] - ldx[6,6] = 0
 lATK: lDW[5,7] - ldx[6,7] - lcan[6,7] >= 0
 lATK: - lDW[5,7] + lDX[6,7] + lcan[6,7] >= 0
 lATK: lDW[5,7] - lDX[6,7] + lstk[6,7] - 2 lcan[6,7] >= 0
 lATK: lDX[6,7] - lstk[6,7] + lcan[6,7] >= 0
 lATK: ldw[5,7] - ldx[6,7] = 0
 lATK: lDW[5,8] - ldx[6,8] - lcan[6,8] >= 0
 lATK: - lDW[5,8] + lDX[6,8] + lcan[6,8] >= 0
 lATK: lDW[5,8] - lDX[6,8] + lstk[6,8] - 2 lcan[6,8] >= 0
 lATK: lDX[6,8] - lstk[6,8] + lcan[6,8] >= 0
 lATK: ldw[5,8] - ldx[6,8] = 0
 lATK: lDW[5,9] - ldx[6,9] - lcan[6,9] >= 0
 lATK: - lDW[5,9] + lDX[6,9] + lcan[6,9] >= 0
 lATK: lDW[5,9] - lDX[6,9] + lstk[6,9] - 2 lcan[6,9] >= 0
 lATK: lDX[6,9] - lstk[6,9] + lcan[6,9] >= 0
 lATK: ldw[5,9] - ldx[6,9] = 0
 lATK: lDW[5,10] - ldx[6,10] - lcan[6,10] >= 0
 lATK: - lDW[5,10] + lDX[6,10] + lcan[6,10] >= 0
 lATK: lDW[5,10] - lDX[6,10] + lstk[6,10] - 2 lcan[6,10] >= 0
 lATK: lDX[6,10] - lstk[6,10] + lcan[6,10] >= 0
 lATK: ldw[5,10] - ldx[6,10] = 0
 lATK: lDW[5,11] - ldx[6,11] - lcan[6,11] >= 0
 lATK: - lDW[5,11] + lDX[6,11] + lcan[6,11] >= 0
 lATK: lDW[5,11] - lDX[6,11] + lstk[6,11] - 2 lcan[6,11] >= 0
 lATK: lDX[6,11] - lstk[6,11] + lcan[6,11] >= 0
 lATK: ldw[5,11] - ldx[6,11] = 0
 lATK: lDW[5,12] - ldx[6,12] - lcan[6,12] >= 0
 lATK: - lDW[5,12] + lDX[6,12] + lcan[6,12] >= 0
 lATK: lDW[5,12] - lDX[6,12] + lstk[6,12] - 2 lcan[6,12] >= 0
 lATK: lDX[6,12] - lstk[6,12] + lcan[6,12] >= 0
 lATK: ldw[5,12] - ldx[6,12] = 0
 lATK: lDW[5,13] - ldx[6,13] - lcan[6,13] >= 0
 lATK: - lDW[5,13] + lDX[6,13] + lcan[6,13] >= 0
 lATK: lDW[5,13] - lDX[6,13] + lstk[6,13] - 2 lcan[6,13] >= 0
 lATK: lDX[6,13] - lstk[6,13] + lcan[6,13] >= 0
 lATK: ldw[5,13] - ldx[6,13] = 0
 lATK: lDW[5,14] - ldx[6,14] - lcan[6,14] >= 0
 lATK: - lDW[5,14] + lDX[6,14] + lcan[6,14] >= 0
 lATK: lDW[5,14] - lDX[6,14] + lstk[6,14] - 2 lcan[6,14] >= 0
 lATK: lDX[6,14] - lstk[6,14] + lcan[6,14] >= 0
 lATK: ldw[5,14] - ldx[6,14] = 0
 lATK: lDW[5,15] - ldx[6,15] - lcan[6,15] >= 0
 lATK: - lDW[5,15] + lDX[6,15] + lcan[6,15] >= 0
 lATK: lDW[5,15] - lDX[6,15] + lstk[6,15] - 2 lcan[6,15] >= 0
 lATK: lDX[6,15] - lstk[6,15] + lcan[6,15] >= 0
 lATK: ldw[5,15] - ldx[6,15] = 0
 lATK: lDW[6,0] - ldx[7,0] - lcan[7,0] >= 0
 lATK: - lDW[6,0] + lDX[7,0] + lcan[7,0] >= 0
 lATK: lDW[6,0] - lDX[7,0] + lstk[7,0] - 2 lcan[7,0] >= 0
 lATK: lDX[7,0] - lstk[7,0] + lcan[7,0] >= 0
 lATK: ldw[6,0] - ldx[7,0] = 0
 lATK: lDW[6,1] - ldx[7,1] - lcan[7,1] >= 0
 lATK: - lDW[6,1] + lDX[7,1] + lcan[7,1] >= 0
 lATK: lDW[6,1] - lDX[7,1] + lstk[7,1] - 2 lcan[7,1] >= 0
 lATK: lDX[7,1] - lstk[7,1] + lcan[7,1] >= 0
 lATK: ldw[6,1] - ldx[7,1] = 0
 lATK: lDW[6,2] - ldx[7,2] - lcan[7,2] >= 0
 lATK: - lDW[6,2] + lDX[7,2] + lcan[7,2] >= 0
 lATK: lDW[6,2] - lDX[7,2] + lstk[7,2] - 2 lcan[7,2] >= 0
 lATK: lDX[7,2] - lstk[7,2] + lcan[7,2] >= 0
 lATK: ldw[6,2] - ldx[7,2] = 0
 lATK: lDW[6,3] - ldx[7,3] - lcan[7,3] >= 0
 lATK: - lDW[6,3] + lDX[7,3] + lcan[7,3] >= 0
 lATK: lDW[6,3] - lDX[7,3] + lstk[7,3] - 2 lcan[7,3] >= 0
 lATK: lDX[7,3] - lstk[7,3] + lcan[7,3] >= 0
 lATK: ldw[6,3] - ldx[7,3] = 0
 lATK: lDW[6,4] - ldx[7,4] - lcan[7,4] >= 0
 lATK: - lDW[6,4] + lDX[7,4] + lcan[7,4] >= 0
 lATK: lDW[6,4] - lDX[7,4] + lstk[7,4] - 2 lcan[7,4] >= 0
 lATK: lDX[7,4] - lstk[7,4] + lcan[7,4] >= 0
 lATK: ldw[6,4] - ldx[7,4] = 0
 lATK: lDW[6,5] - ldx[7,5] - lcan[7,5] >= 0
 lATK: - lDW[6,5] + lDX[7,5] + lcan[7,5] >= 0
 lATK: lDW[6,5] - lDX[7,5] + lstk[7,5] - 2 lcan[7,5] >= 0
 lATK: lDX[7,5] - lstk[7,5] + lcan[7,5] >= 0
 lATK: ldw[6,5] - ldx[7,5] = 0
 lATK: lDW[6,6] - ldx[7,6] - lcan[7,6] >= 0
 lATK: - lDW[6,6] + lDX[7,6] + lcan[7,6] >= 0
 lATK: lDW[6,6] - lDX[7,6] + lstk[7,6] - 2 lcan[7,6] >= 0
 lATK: lDX[7,6] - lstk[7,6] + lcan[7,6] >= 0
 lATK: ldw[6,6] - ldx[7,6] = 0
 lATK: lDW[6,7] - ldx[7,7] - lcan[7,7] >= 0
 lATK: - lDW[6,7] + lDX[7,7] + lcan[7,7] >= 0
 lATK: lDW[6,7] - lDX[7,7] + lstk[7,7] - 2 lcan[7,7] >= 0
 lATK: lDX[7,7] - lstk[7,7] + lcan[7,7] >= 0
 lATK: ldw[6,7] - ldx[7,7] = 0
 lATK: lDW[6,8] - ldx[7,8] - lcan[7,8] >= 0
 lATK: - lDW[6,8] + lDX[7,8] + lcan[7,8] >= 0
 lATK: lDW[6,8] - lDX[7,8] + lstk[7,8] - 2 lcan[7,8] >= 0
 lATK: lDX[7,8] - lstk[7,8] + lcan[7,8] >= 0
 lATK: ldw[6,8] - ldx[7,8] = 0
 lATK: lDW[6,9] - ldx[7,9] - lcan[7,9] >= 0
 lATK: - lDW[6,9] + lDX[7,9] + lcan[7,9] >= 0
 lATK: lDW[6,9] - lDX[7,9] + lstk[7,9] - 2 lcan[7,9] >= 0
 lATK: lDX[7,9] - lstk[7,9] + lcan[7,9] >= 0
 lATK: ldw[6,9] - ldx[7,9] = 0
 lATK: lDW[6,10] - ldx[7,10] - lcan[7,10] >= 0
 lATK: - lDW[6,10] + lDX[7,10] + lcan[7,10] >= 0
 lATK: lDW[6,10] - lDX[7,10] + lstk[7,10] - 2 lcan[7,10] >= 0
 lATK: lDX[7,10] - lstk[7,10] + lcan[7,10] >= 0
 lATK: ldw[6,10] - ldx[7,10] = 0
 lATK: lDW[6,11] - ldx[7,11] - lcan[7,11] >= 0
 lATK: - lDW[6,11] + lDX[7,11] + lcan[7,11] >= 0
 lATK: lDW[6,11] - lDX[7,11] + lstk[7,11] - 2 lcan[7,11] >= 0
 lATK: lDX[7,11] - lstk[7,11] + lcan[7,11] >= 0
 lATK: ldw[6,11] - ldx[7,11] = 0
 lATK: lDW[6,12] - ldx[7,12] - lcan[7,12] >= 0
 lATK: - lDW[6,12] + lDX[7,12] + lcan[7,12] >= 0
 lATK: lDW[6,12] - lDX[7,12] + lstk[7,12] - 2 lcan[7,12] >= 0
 lATK: lDX[7,12] - lstk[7,12] + lcan[7,12] >= 0
 lATK: ldw[6,12] - ldx[7,12] = 0
 lATK: lDW[6,13] - ldx[7,13] - lcan[7,13] >= 0
 lATK: - lDW[6,13] + lDX[7,13] + lcan[7,13] >= 0
 lATK: lDW[6,13] - lDX[7,13] + lstk[7,13] - 2 lcan[7,13] >= 0
 lATK: lDX[7,13] - lstk[7,13] + lcan[7,13] >= 0
 lATK: ldw[6,13] - ldx[7,13] = 0
 lATK: lDW[6,14] - ldx[7,14] - lcan[7,14] >= 0
 lATK: - lDW[6,14] + lDX[7,14] + lcan[7,14] >= 0
 lATK: lDW[6,14] - lDX[7,14] + lstk[7,14] - 2 lcan[7,14] >= 0
 lATK: lDX[7,14] - lstk[7,14] + lcan[7,14] >= 0
 lATK: ldw[6,14] - ldx[7,14] = 0
 lATK: lDW[6,15] - ldx[7,15] - lcan[7,15] >= 0
 lATK: - lDW[6,15] + lDX[7,15] + lcan[7,15] >= 0
 lATK: lDW[6,15] - lDX[7,15] + lstk[7,15] - 2 lcan[7,15] >= 0
 lATK: lDX[7,15] - lstk[7,15] + lcan[7,15] >= 0
 lATK: ldw[6,15] - ldx[7,15] = 0
 lATK: lDW[7,0] - ldx[8,0] - lcan[8,0] >= 0
 lATK: - lDW[7,0] + lDX[8,0] + lcan[8,0] >= 0
 lATK: lDW[7,0] - lDX[8,0] + lstk[8,0] - 2 lcan[8,0] >= 0
 lATK: lDX[8,0] - lstk[8,0] + lcan[8,0] >= 0
 lATK: ldw[7,0] - ldx[8,0] = 0
 lATK: lDW[7,1] - ldx[8,1] - lcan[8,1] >= 0
 lATK: - lDW[7,1] + lDX[8,1] + lcan[8,1] >= 0
 lATK: lDW[7,1] - lDX[8,1] + lstk[8,1] - 2 lcan[8,1] >= 0
 lATK: lDX[8,1] - lstk[8,1] + lcan[8,1] >= 0
 lATK: ldw[7,1] - ldx[8,1] = 0
 lATK: lDW[7,2] - ldx[8,2] - lcan[8,2] >= 0
 lATK: - lDW[7,2] + lDX[8,2] + lcan[8,2] >= 0
 lATK: lDW[7,2] - lDX[8,2] + lstk[8,2] - 2 lcan[8,2] >= 0
 lATK: lDX[8,2] - lstk[8,2] + lcan[8,2] >= 0
 lATK: ldw[7,2] - ldx[8,2] = 0
 lATK: lDW[7,3] - ldx[8,3] - lcan[8,3] >= 0
 lATK: - lDW[7,3] + lDX[8,3] + lcan[8,3] >= 0
 lATK: lDW[7,3] - lDX[8,3] + lstk[8,3] - 2 lcan[8,3] >= 0
 lATK: lDX[8,3] - lstk[8,3] + lcan[8,3] >= 0
 lATK: ldw[7,3] - ldx[8,3] = 0
 lATK: lDW[7,4] - ldx[8,4] - lcan[8,4] >= 0
 lATK: - lDW[7,4] + lDX[8,4] + lcan[8,4] >= 0
 lATK: lDW[7,4] - lDX[8,4] + lstk[8,4] - 2 lcan[8,4] >= 0
 lATK: lDX[8,4] - lstk[8,4] + lcan[8,4] >= 0
 lATK: ldw[7,4] - ldx[8,4] = 0
 lATK: lDW[7,5] - ldx[8,5] - lcan[8,5] >= 0
 lATK: - lDW[7,5] + lDX[8,5] + lcan[8,5] >= 0
 lATK: lDW[7,5] - lDX[8,5] + lstk[8,5] - 2 lcan[8,5] >= 0
 lATK: lDX[8,5] - lstk[8,5] + lcan[8,5] >= 0
 lATK: ldw[7,5] - ldx[8,5] = 0
 lATK: lDW[7,6] - ldx[8,6] - lcan[8,6] >= 0
 lATK: - lDW[7,6] + lDX[8,6] + lcan[8,6] >= 0
 lATK: lDW[7,6] - lDX[8,6] + lstk[8,6] - 2 lcan[8,6] >= 0
 lATK: lDX[8,6] - lstk[8,6] + lcan[8,6] >= 0
 lATK: ldw[7,6] - ldx[8,6] = 0
 lATK: lDW[7,7] - ldx[8,7] - lcan[8,7] >= 0
 lATK: - lDW[7,7] + lDX[8,7] + lcan[8,7] >= 0
 lATK: lDW[7,7] - lDX[8,7] + lstk[8,7] - 2 lcan[8,7] >= 0
 lATK: lDX[8,7] - lstk[8,7] + lcan[8,7] >= 0
 lATK: ldw[7,7] - ldx[8,7] = 0
 lATK: lDW[7,8] - ldx[8,8] - lcan[8,8] >= 0
 lATK: - lDW[7,8] + lDX[8,8] + lcan[8,8] >= 0
 lATK: lDW[7,8] - lDX[8,8] + lstk[8,8] - 2 lcan[8,8] >= 0
 lATK: lDX[8,8] - lstk[8,8] + lcan[8,8] >= 0
 lATK: ldw[7,8] - ldx[8,8] = 0
 lATK: lDW[7,9] - ldx[8,9] - lcan[8,9] >= 0
 lATK: - lDW[7,9] + lDX[8,9] + lcan[8,9] >= 0
 lATK: lDW[7,9] - lDX[8,9] + lstk[8,9] - 2 lcan[8,9] >= 0
 lATK: lDX[8,9] - lstk[8,9] + lcan[8,9] >= 0
 lATK: ldw[7,9] - ldx[8,9] = 0
 lATK: lDW[7,10] - ldx[8,10] - lcan[8,10] >= 0
 lATK: - lDW[7,10] + lDX[8,10] + lcan[8,10] >= 0
 lATK: lDW[7,10] - lDX[8,10] + lstk[8,10] - 2 lcan[8,10] >= 0
 lATK: lDX[8,10] - lstk[8,10] + lcan[8,10] >= 0
 lATK: ldw[7,10] - ldx[8,10] = 0
 lATK: lDW[7,11] - ldx[8,11] - lcan[8,11] >= 0
 lATK: - lDW[7,11] + lDX[8,11] + lcan[8,11] >= 0
 lATK: lDW[7,11] - lDX[8,11] + lstk[8,11] - 2 lcan[8,11] >= 0
 lATK: lDX[8,11] - lstk[8,11] + lcan[8,11] >= 0
 lATK: ldw[7,11] - ldx[8,11] = 0
 lATK: lDW[7,12] - ldx[8,12] - lcan[8,12] >= 0
 lATK: - lDW[7,12] + lDX[8,12] + lcan[8,12] >= 0
 lATK: lDW[7,12] - lDX[8,12] + lstk[8,12] - 2 lcan[8,12] >= 0
 lATK: lDX[8,12] - lstk[8,12] + lcan[8,12] >= 0
 lATK: ldw[7,12] - ldx[8,12] = 0
 lATK: lDW[7,13] - ldx[8,13] - lcan[8,13] >= 0
 lATK: - lDW[7,13] + lDX[8,13] + lcan[8,13] >= 0
 lATK: lDW[7,13] - lDX[8,13] + lstk[8,13] - 2 lcan[8,13] >= 0
 lATK: lDX[8,13] - lstk[8,13] + lcan[8,13] >= 0
 lATK: ldw[7,13] - ldx[8,13] = 0
 lATK: lDW[7,14] - ldx[8,14] - lcan[8,14] >= 0
 lATK: - lDW[7,14] + lDX[8,14] + lcan[8,14] >= 0
 lATK: lDW[7,14] - lDX[8,14] + lstk[8,14] - 2 lcan[8,14] >= 0
 lATK: lDX[8,14] - lstk[8,14] + lcan[8,14] >= 0
 lATK: ldw[7,14] - ldx[8,14] = 0
 lATK: lDW[7,15] - ldx[8,15] - lcan[8,15] >= 0
 lATK: - lDW[7,15] + lDX[8,15] + lcan[8,15] >= 0
 lATK: lDW[7,15] - lDX[8,15] + lstk[8,15] - 2 lcan[8,15] >= 0
 lATK: lDX[8,15] - lstk[8,15] + lcan[8,15] >= 0
 lATK: ldw[7,15] - ldx[8,15] = 0
 lSC: lDX[1,0] - lDY[1,0] = 0
 lSC: ldx[1,0] - lDY[1,0] = 0
 lSC: lDX[1,1] - lDY[1,1] = 0
 lSC: ldx[1,1] - lDY[1,1] = 0
 lSC: lDX[1,2] - lDY[1,2] = 0
 lSC: ldx[1,2] - lDY[1,2] = 0
 lSC: lDX[1,3] - lDY[1,3] = 0
 lSC: ldx[1,3] - lDY[1,3] = 0
 lSC: lDX[1,4] - lDY[1,4] = 0
 lSC: ldx[1,4] - lDY[1,4] = 0
 lSC: lDX[1,5] - lDY[1,5] = 0
 lSC: ldx[1,5] - lDY[1,5] = 0
 lSC: lDX[1,6] - lDY[1,6] = 0
 lSC: ldx[1,6] - lDY[1,6] = 0
 lSC: lDX[1,7] - lDY[1,7] = 0
 lSC: ldx[1,7] - lDY[1,7] = 0
 lSC: lDX[1,8] - lDY[1,8] = 0
 lSC: ldx[1,8] - lDY[1,8] = 0
 lSC: lDX[1,9] - lDY[1,9] = 0
 lSC: ldx[1,9] - lDY[1,9] = 0
 lSC: lDX[1,10] - lDY[1,10] = 0
 lSC: ldx[1,10] - lDY[1,10] = 0
 lSC: lDX[1,11] - lDY[1,11] = 0
 lSC: ldx[1,11] - lDY[1,11] = 0
 lSC: lDX[1,12] - lDY[1,12] = 0
 lSC: ldx[1,12] - lDY[1,12] = 0
 lSC: lDX[1,13] - lDY[1,13] = 0
 lSC: ldx[1,13] - lDY[1,13] = 0
 lSC: lDX[1,14] - lDY[1,14] = 0
 lSC: ldx[1,14] - lDY[1,14] = 0
 lSC: lDX[1,15] - lDY[1,15] = 0
 lSC: ldx[1,15] - lDY[1,15] = 0
 lSC: lDX[2,0] - lDY[2,0] = 0
 lSC: ldx[2,0] - lDY[2,0] = 0
 lSC: lDX[2,1] - lDY[2,1] = 0
 lSC: ldx[2,1] - lDY[2,1] = 0
 lSC: lDX[2,2] - lDY[2,2] = 0
 lSC: ldx[2,2] - lDY[2,2] = 0
 lSC: lDX[2,3] - lDY[2,3] = 0
 lSC: ldx[2,3] - lDY[2,3] = 0
 lSC: lDX[2,4] - lDY[2,4] = 0
 lSC: ldx[2,4] - lDY[2,4] = 0
 lSC: lDX[2,5] - lDY[2,5] = 0
 lSC: ldx[2,5] - lDY[2,5] = 0
 lSC: lDX[2,6] - lDY[2,6] = 0
 lSC: ldx[2,6] - lDY[2,6] = 0
 lSC: lDX[2,7] - lDY[2,7] = 0
 lSC: ldx[2,7] - lDY[2,7] = 0
 lSC: lDX[2,8] - lDY[2,8] = 0
 lSC: ldx[2,8] - lDY[2,8] = 0
 lSC: lDX[2,9] - lDY[2,9] = 0
 lSC: ldx[2,9] - lDY[2,9] = 0
 lSC: lDX[2,10] - lDY[2,10] = 0
 lSC: ldx[2,10] - lDY[2,10] = 0
 lSC: lDX[2,11] - lDY[2,11] = 0
 lSC: ldx[2,11] - lDY[2,11] = 0
 lSC: lDX[2,12] - lDY[2,12] = 0
 lSC: ldx[2,12] - lDY[2,12] = 0
 lSC: lDX[2,13] - lDY[2,13] = 0
 lSC: ldx[2,13] - lDY[2,13] = 0
 lSC: lDX[2,14] - lDY[2,14] = 0
 lSC: ldx[2,14] - lDY[2,14] = 0
 lSC: lDX[2,15] - lDY[2,15] = 0
 lSC: ldx[2,15] - lDY[2,15] = 0
 lSC: lDX[3,0] - lDY[3,0] = 0
 lSC: ldx[3,0] - lDY[3,0] = 0
 lSC: lDX[3,1] - lDY[3,1] = 0
 lSC: ldx[3,1] - lDY[3,1] = 0
 lSC: lDX[3,2] - lDY[3,2] = 0
 lSC: ldx[3,2] - lDY[3,2] = 0
 lSC: lDX[3,3] - lDY[3,3] = 0
 lSC: ldx[3,3] - lDY[3,3] = 0
 lSC: lDX[3,4] - lDY[3,4] = 0
 lSC: ldx[3,4] - lDY[3,4] = 0
 lSC: lDX[3,5] - lDY[3,5] = 0
 lSC: ldx[3,5] - lDY[3,5] = 0
 lSC: lDX[3,6] - lDY[3,6] = 0
 lSC: ldx[3,6] - lDY[3,6] = 0
 lSC: lDX[3,7] - lDY[3,7] = 0
 lSC: ldx[3,7] - lDY[3,7] = 0
 lSC: lDX[3,8] - lDY[3,8] = 0
 lSC: ldx[3,8] - lDY[3,8] = 0
 lSC: lDX[3,9] - lDY[3,9] = 0
 lSC: ldx[3,9] - lDY[3,9] = 0
 lSC: lDX[3,10] - lDY[3,10] = 0
 lSC: ldx[3,10] - lDY[3,10] = 0
 lSC: lDX[3,11] - lDY[3,11] = 0
 lSC: ldx[3,11] - lDY[3,11] = 0
 lSC: lDX[3,12] - lDY[3,12] = 0
 lSC: ldx[3,12] - lDY[3,12] = 0
 lSC: lDX[3,13] - lDY[3,13] = 0
 lSC: ldx[3,13] - lDY[3,13] = 0
 lSC: lDX[3,14] - lDY[3,14] = 0
 lSC: ldx[3,14] - lDY[3,14] = 0
 lSC: lDX[3,15] - lDY[3,15] = 0
 lSC: ldx[3,15] - lDY[3,15] = 0
 lSC: lDX[4,0] - lDY[4,0] = 0
 lSC: ldx[4,0] - lDY[4,0] = 0
 lSC: lDX[4,1] - lDY[4,1] = 0
 lSC: ldx[4,1] - lDY[4,1] = 0
 lSC: lDX[4,2] - lDY[4,2] = 0
 lSC: ldx[4,2] - lDY[4,2] = 0
 lSC: lDX[4,3] - lDY[4,3] = 0
 lSC: ldx[4,3] - lDY[4,3] = 0
 lSC: lDX[4,4] - lDY[4,4] = 0
 lSC: ldx[4,4] - lDY[4,4] = 0
 lSC: lDX[4,5] - lDY[4,5] = 0
 lSC: ldx[4,5] - lDY[4,5] = 0
 lSC: lDX[4,6] - lDY[4,6] = 0
 lSC: ldx[4,6] - lDY[4,6] = 0
 lSC: lDX[4,7] - lDY[4,7] = 0
 lSC: ldx[4,7] - lDY[4,7] = 0
 lSC: lDX[4,8] - lDY[4,8] = 0
 lSC: ldx[4,8] - lDY[4,8] = 0
 lSC: lDX[4,9] - lDY[4,9] = 0
 lSC: ldx[4,9] - lDY[4,9] = 0
 lSC: lDX[4,10] - lDY[4,10] = 0
 lSC: ldx[4,10] - lDY[4,10] = 0
 lSC: lDX[4,11] - lDY[4,11] = 0
 lSC: ldx[4,11] - lDY[4,11] = 0
 lSC: lDX[4,12] - lDY[4,12] = 0
 lSC: ldx[4,12] - lDY[4,12] = 0
 lSC: lDX[4,13] - lDY[4,13] = 0
 lSC: ldx[4,13] - lDY[4,13] = 0
 lSC: lDX[4,14] - lDY[4,14] = 0
 lSC: ldx[4,14] - lDY[4,14] = 0
 lSC: lDX[4,15] - lDY[4,15] = 0
 lSC: ldx[4,15] - lDY[4,15] = 0
 lSC: lDX[5,0] - lDY[5,0] = 0
 lSC: ldx[5,0] - lDY[5,0] = 0
 lSC: lDX[5,1] - lDY[5,1] = 0
 lSC: ldx[5,1] - lDY[5,1] = 0
 lSC: lDX[5,2] - lDY[5,2] = 0
 lSC: ldx[5,2] - lDY[5,2] = 0
 lSC: lDX[5,3] - lDY[5,3] = 0
 lSC: ldx[5,3] - lDY[5,3] = 0
 lSC: lDX[5,4] - lDY[5,4] = 0
 lSC: ldx[5,4] - lDY[5,4] = 0
 lSC: lDX[5,5] - lDY[5,5] = 0
 lSC: ldx[5,5] - lDY[5,5] = 0
 lSC: lDX[5,6] - lDY[5,6] = 0
 lSC: ldx[5,6] - lDY[5,6] = 0
 lSC: lDX[5,7] - lDY[5,7] = 0
 lSC: ldx[5,7] - lDY[5,7] = 0
 lSC: lDX[5,8] - lDY[5,8] = 0
 lSC: ldx[5,8] - lDY[5,8] = 0
 lSC: lDX[5,9] - lDY[5,9] = 0
 lSC: ldx[5,9] - lDY[5,9] = 0
 lSC: lDX[5,10] - lDY[5,10] = 0
 lSC: ldx[5,10] - lDY[5,10] = 0
 lSC: lDX[5,11] - lDY[5,11] = 0
 lSC: ldx[5,11] - lDY[5,11] = 0
 lSC: lDX[5,12] - lDY[5,12] = 0
 lSC: ldx[5,12] - lDY[5,12] = 0
 lSC: lDX[5,13] - lDY[5,13] = 0
 lSC: ldx[5,13] - lDY[5,13] = 0
 lSC: lDX[5,14] - lDY[5,14] = 0
 lSC: ldx[5,14] - lDY[5,14] = 0
 lSC: lDX[5,15] - lDY[5,15] = 0
 lSC: ldx[5,15] - lDY[5,15] = 0
 lSC: lDX[6,0] - lDY[6,0] = 0
 lSC: ldx[6,0] - lDY[6,0] = 0
 lSC: lDX[6,1] - lDY[6,1] = 0
 lSC: ldx[6,1] - lDY[6,1] = 0
 lSC: lDX[6,2] - lDY[6,2] = 0
 lSC: ldx[6,2] - lDY[6,2] = 0
 lSC: lDX[6,3] - lDY[6,3] = 0
 lSC: ldx[6,3] - lDY[6,3] = 0
 lSC: lDX[6,4] - lDY[6,4] = 0
 lSC: ldx[6,4] - lDY[6,4] = 0
 lSC: lDX[6,5] - lDY[6,5] = 0
 lSC: ldx[6,5] - lDY[6,5] = 0
 lSC: lDX[6,6] - lDY[6,6] = 0
 lSC: ldx[6,6] - lDY[6,6] = 0
 lSC: lDX[6,7] - lDY[6,7] = 0
 lSC: ldx[6,7] - lDY[6,7] = 0
 lSC: lDX[6,8] - lDY[6,8] = 0
 lSC: ldx[6,8] - lDY[6,8] = 0
 lSC: lDX[6,9] - lDY[6,9] = 0
 lSC: ldx[6,9] - lDY[6,9] = 0
 lSC: lDX[6,10] - lDY[6,10] = 0
 lSC: ldx[6,10] - lDY[6,10] = 0
 lSC: lDX[6,11] - lDY[6,11] = 0
 lSC: ldx[6,11] - lDY[6,11] = 0
 lSC: lDX[6,12] - lDY[6,12] = 0
 lSC: ldx[6,12] - lDY[6,12] = 0
 lSC: lDX[6,13] - lDY[6,13] = 0
 lSC: ldx[6,13] - lDY[6,13] = 0
 lSC: lDX[6,14] - lDY[6,14] = 0
 lSC: ldx[6,14] - lDY[6,14] = 0
 lSC: lDX[6,15] - lDY[6,15] = 0
 lSC: ldx[6,15] - lDY[6,15] = 0
 lSC: lDX[7,0] - lDY[7,0] = 0
 lSC: ldx[7,0] - lDY[7,0] = 0
 lSC: lDX[7,1] - lDY[7,1] = 0
 lSC: ldx[7,1] - lDY[7,1] = 0
 lSC: lDX[7,2] - lDY[7,2] = 0
 lSC: ldx[7,2] - lDY[7,2] = 0
 lSC: lDX[7,3] - lDY[7,3] = 0
 lSC: ldx[7,3] - lDY[7,3] = 0
 lSC: lDX[7,4] - lDY[7,4] = 0
 lSC: ldx[7,4] - lDY[7,4] = 0
 lSC: lDX[7,5] - lDY[7,5] = 0
 lSC: ldx[7,5] - lDY[7,5] = 0
 lSC: lDX[7,6] - lDY[7,6] = 0
 lSC: ldx[7,6] - lDY[7,6] = 0
 lSC: lDX[7,7] - lDY[7,7] = 0
 lSC: ldx[7,7] - lDY[7,7] = 0
 lSC: lDX[7,8] - lDY[7,8] = 0
 lSC: ldx[7,8] - lDY[7,8] = 0
 lSC: lDX[7,9] - lDY[7,9] = 0
 lSC: ldx[7,9] - lDY[7,9] = 0
 lSC: lDX[7,10] - lDY[7,10] = 0
 lSC: ldx[7,10] - lDY[7,10] = 0
 lSC: lDX[7,11] - lDY[7,11] = 0
 lSC: ldx[7,11] - lDY[7,11] = 0
 lSC: lDX[7,12] - lDY[7,12] = 0
 lSC: ldx[7,12] - lDY[7,12] = 0
 lSC: lDX[7,13] - lDY[7,13] = 0
 lSC: ldx[7,13] - lDY[7,13] = 0
 lSC: lDX[7,14] - lDY[7,14] = 0
 lSC: ldx[7,14] - lDY[7,14] = 0
 lSC: lDX[7,15] - lDY[7,15] = 0
 lSC: ldx[7,15] - lDY[7,15] = 0
 lSC: lDY[1,0] - lDZ[1,0] = 0
 lSC: ldy[1,0] - ldz[1,0] = 0
 lSC: lDY[1,5] - lDZ[1,1] = 0
 lSC: ldy[1,5] - ldz[1,1] = 0
 lSC: lDY[1,10] - lDZ[1,2] = 0
 lSC: ldy[1,10] - ldz[1,2] = 0
 lSC: lDY[1,15] - lDZ[1,3] = 0
 lSC: ldy[1,15] - ldz[1,3] = 0
 lSC: lDY[1,4] - lDZ[1,4] = 0
 lSC: ldy[1,4] - ldz[1,4] = 0
 lSC: lDY[1,9] - lDZ[1,5] = 0
 lSC: ldy[1,9] - ldz[1,5] = 0
 lSC: lDY[1,14] - lDZ[1,6] = 0
 lSC: ldy[1,14] - ldz[1,6] = 0
 lSC: lDY[1,3] - lDZ[1,7] = 0
 lSC: ldy[1,3] - ldz[1,7] = 0
 lSC: lDY[1,8] - lDZ[1,8] = 0
 lSC: ldy[1,8] - ldz[1,8] = 0
 lSC: lDY[1,13] - lDZ[1,9] = 0
 lSC: ldy[1,13] - ldz[1,9] = 0
 lSC: lDY[1,2] - lDZ[1,10] = 0
 lSC: ldy[1,2] - ldz[1,10] = 0
 lSC: lDY[1,7] - lDZ[1,11] = 0
 lSC: ldy[1,7] - ldz[1,11] = 0
 lSC: lDY[1,12] - lDZ[1,12] = 0
 lSC: ldy[1,12] - ldz[1,12] = 0
 lSC: lDY[1,1] - lDZ[1,13] = 0
 lSC: ldy[1,1] - ldz[1,13] = 0
 lSC: lDY[1,6] - lDZ[1,14] = 0
 lSC: ldy[1,6] - ldz[1,14] = 0
 lSC: lDY[1,11] - lDZ[1,15] = 0
 lSC: ldy[1,11] - ldz[1,15] = 0
 lSC: lDY[2,0] - lDZ[2,0] = 0
 lSC: ldy[2,0] - ldz[2,0] = 0
 lSC: lDY[2,5] - lDZ[2,1] = 0
 lSC: ldy[2,5] - ldz[2,1] = 0
 lSC: lDY[2,10] - lDZ[2,2] = 0
 lSC: ldy[2,10] - ldz[2,2] = 0
 lSC: lDY[2,15] - lDZ[2,3] = 0
 lSC: ldy[2,15] - ldz[2,3] = 0
 lSC: lDY[2,4] - lDZ[2,4] = 0
 lSC: ldy[2,4] - ldz[2,4] = 0
 lSC: lDY[2,9] - lDZ[2,5] = 0
 lSC: ldy[2,9] - ldz[2,5] = 0
 lSC: lDY[2,14] - lDZ[2,6] = 0
 lSC: ldy[2,14] - ldz[2,6] = 0
 lSC: lDY[2,3] - lDZ[2,7] = 0
 lSC: ldy[2,3] - ldz[2,7] = 0
 lSC: lDY[2,8] - lDZ[2,8] = 0
 lSC: ldy[2,8] - ldz[2,8] = 0
 lSC: lDY[2,13] - lDZ[2,9] = 0
 lSC: ldy[2,13] - ldz[2,9] = 0
 lSC: lDY[2,2] - lDZ[2,10] = 0
 lSC: ldy[2,2] - ldz[2,10] = 0
 lSC: lDY[2,7] - lDZ[2,11] = 0
 lSC: ldy[2,7] - ldz[2,11] = 0
 lSC: lDY[2,12] - lDZ[2,12] = 0
 lSC: ldy[2,12] - ldz[2,12] = 0
 lSC: lDY[2,1] - lDZ[2,13] = 0
 lSC: ldy[2,1] - ldz[2,13] = 0
 lSC: lDY[2,6] - lDZ[2,14] = 0
 lSC: ldy[2,6] - ldz[2,14] = 0
 lSC: lDY[2,11] - lDZ[2,15] = 0
 lSC: ldy[2,11] - ldz[2,15] = 0
 lSC: lDY[3,0] - lDZ[3,0] = 0
 lSC: ldy[3,0] - ldz[3,0] = 0
 lSC: lDY[3,5] - lDZ[3,1] = 0
 lSC: ldy[3,5] - ldz[3,1] = 0
 lSC: lDY[3,10] - lDZ[3,2] = 0
 lSC: ldy[3,10] - ldz[3,2] = 0
 lSC: lDY[3,15] - lDZ[3,3] = 0
 lSC: ldy[3,15] - ldz[3,3] = 0
 lSC: lDY[3,4] - lDZ[3,4] = 0
 lSC: ldy[3,4] - ldz[3,4] = 0
 lSC: lDY[3,9] - lDZ[3,5] = 0
 lSC: ldy[3,9] - ldz[3,5] = 0
 lSC: lDY[3,14] - lDZ[3,6] = 0
 lSC: ldy[3,14] - ldz[3,6] = 0
 lSC: lDY[3,3] - lDZ[3,7] = 0
 lSC: ldy[3,3] - ldz[3,7] = 0
 lSC: lDY[3,8] - lDZ[3,8] = 0
 lSC: ldy[3,8] - ldz[3,8] = 0
 lSC: lDY[3,13] - lDZ[3,9] = 0
 lSC: ldy[3,13] - ldz[3,9] = 0
 lSC: lDY[3,2] - lDZ[3,10] = 0
 lSC: ldy[3,2] - ldz[3,10] = 0
 lSC: lDY[3,7] - lDZ[3,11] = 0
 lSC: ldy[3,7] - ldz[3,11] = 0
 lSC: lDY[3,12] - lDZ[3,12] = 0
 lSC: ldy[3,12] - ldz[3,12] = 0
 lSC: lDY[3,1] - lDZ[3,13] = 0
 lSC: ldy[3,1] - ldz[3,13] = 0
 lSC: lDY[3,6] - lDZ[3,14] = 0
 lSC: ldy[3,6] - ldz[3,14] = 0
 lSC: lDY[3,11] - lDZ[3,15] = 0
 lSC: ldy[3,11] - ldz[3,15] = 0
 lSC: lDY[4,0] - lDZ[4,0] = 0
 lSC: ldy[4,0] - ldz[4,0] = 0
 lSC: lDY[4,5] - lDZ[4,1] = 0
 lSC: ldy[4,5] - ldz[4,1] = 0
 lSC: lDY[4,10] - lDZ[4,2] = 0
 lSC: ldy[4,10] - ldz[4,2] = 0
 lSC: lDY[4,15] - lDZ[4,3] = 0
 lSC: ldy[4,15] - ldz[4,3] = 0
 lSC: lDY[4,4] - lDZ[4,4] = 0
 lSC: ldy[4,4] - ldz[4,4] = 0
 lSC: lDY[4,9] - lDZ[4,5] = 0
 lSC: ldy[4,9] - ldz[4,5] = 0
 lSC: lDY[4,14] - lDZ[4,6] = 0
 lSC: ldy[4,14] - ldz[4,6] = 0
 lSC: lDY[4,3] - lDZ[4,7] = 0
 lSC: ldy[4,3] - ldz[4,7] = 0
 lSC: lDY[4,8] - lDZ[4,8] = 0
 lSC: ldy[4,8] - ldz[4,8] = 0
 lSC: lDY[4,13] - lDZ[4,9] = 0
 lSC: ldy[4,13] - ldz[4,9] = 0
 lSC: lDY[4,2] - lDZ[4,10] = 0
 lSC: ldy[4,2] - ldz[4,10] = 0
 lSC: lDY[4,7] - lDZ[4,11] = 0
 lSC: ldy[4,7] - ldz[4,11] = 0
 lSC: lDY[4,12] - lDZ[4,12] = 0
 lSC: ldy[4,12] - ldz[4,12] = 0
 lSC: lDY[4,1] - lDZ[4,13] = 0
 lSC: ldy[4,1] - ldz[4,13] = 0
 lSC: lDY[4,6] - lDZ[4,14] = 0
 lSC: ldy[4,6] - ldz[4,14] = 0
 lSC: lDY[4,11] - lDZ[4,15] = 0
 lSC: ldy[4,11] - ldz[4,15] = 0
 lSC: lDY[5,0] - lDZ[5,0] = 0
 lSC: ldy[5,0] - ldz[5,0] = 0
 lSC: lDY[5,5] - lDZ[5,1] = 0
 lSC: ldy[5,5] - ldz[5,1] = 0
 lSC: lDY[5,10] - lDZ[5,2] = 0
 lSC: ldy[5,10] - ldz[5,2] = 0
 lSC: lDY[5,15] - lDZ[5,3] = 0
 lSC: ldy[5,15] - ldz[5,3] = 0
 lSC: lDY[5,4] - lDZ[5,4] = 0
 lSC: ldy[5,4] - ldz[5,4] = 0
 lSC: lDY[5,9] - lDZ[5,5] = 0
 lSC: ldy[5,9] - ldz[5,5] = 0
 lSC: lDY[5,14] - lDZ[5,6] = 0
 lSC: ldy[5,14] - ldz[5,6] = 0
 lSC: lDY[5,3] - lDZ[5,7] = 0
 lSC: ldy[5,3] - ldz[5,7] = 0
 lSC: lDY[5,8] - lDZ[5,8] = 0
 lSC: ldy[5,8] - ldz[5,8] = 0
 lSC: lDY[5,13] - lDZ[5,9] = 0
 lSC: ldy[5,13] - ldz[5,9] = 0
 lSC: lDY[5,2] - lDZ[5,10] = 0
 lSC: ldy[5,2] - ldz[5,10] = 0
 lSC: lDY[5,7] - lDZ[5,11] = 0
 lSC: ldy[5,7] - ldz[5,11] = 0
 lSC: lDY[5,12] - lDZ[5,12] = 0
 lSC: ldy[5,12] - ldz[5,12] = 0
 lSC: lDY[5,1] - lDZ[5,13] = 0
 lSC: ldy[5,1] - ldz[5,13] = 0
 lSC: lDY[5,6] - lDZ[5,14] = 0
 lSC: ldy[5,6] - ldz[5,14] = 0
 lSC: lDY[5,11] - lDZ[5,15] = 0
 lSC: ldy[5,11] - ldz[5,15] = 0
 lSC: lDY[6,0] - lDZ[6,0] = 0
 lSC: ldy[6,0] - ldz[6,0] = 0
 lSC: lDY[6,5] - lDZ[6,1] = 0
 lSC: ldy[6,5] - ldz[6,1] = 0
 lSC: lDY[6,10] - lDZ[6,2] = 0
 lSC: ldy[6,10] - ldz[6,2] = 0
 lSC: lDY[6,15] - lDZ[6,3] = 0
 lSC: ldy[6,15] - ldz[6,3] = 0
 lSC: lDY[6,4] - lDZ[6,4] = 0
 lSC: ldy[6,4] - ldz[6,4] = 0
 lSC: lDY[6,9] - lDZ[6,5] = 0
 lSC: ldy[6,9] - ldz[6,5] = 0
 lSC: lDY[6,14] - lDZ[6,6] = 0
 lSC: ldy[6,14] - ldz[6,6] = 0
 lSC: lDY[6,3] - lDZ[6,7] = 0
 lSC: ldy[6,3] - ldz[6,7] = 0
 lSC: lDY[6,8] - lDZ[6,8] = 0
 lSC: ldy[6,8] - ldz[6,8] = 0
 lSC: lDY[6,13] - lDZ[6,9] = 0
 lSC: ldy[6,13] - ldz[6,9] = 0
 lSC: lDY[6,2] - lDZ[6,10] = 0
 lSC: ldy[6,2] - ldz[6,10] = 0
 lSC: lDY[6,7] - lDZ[6,11] = 0
 lSC: ldy[6,7] - ldz[6,11] = 0
 lSC: lDY[6,12] - lDZ[6,12] = 0
 lSC: ldy[6,12] - ldz[6,12] = 0
 lSC: lDY[6,1] - lDZ[6,13] = 0
 lSC: ldy[6,1] - ldz[6,13] = 0
 lSC: lDY[6,6] - lDZ[6,14] = 0
 lSC: ldy[6,6] - ldz[6,14] = 0
 lSC: lDY[6,11] - lDZ[6,15] = 0
 lSC: ldy[6,11] - ldz[6,15] = 0
 lSC: lDY[7,0] - lDZ[7,0] = 0
 lSC: ldy[7,0] - ldz[7,0] = 0
 lSC: lDY[7,5] - lDZ[7,1] = 0
 lSC: ldy[7,5] - ldz[7,1] = 0
 lSC: lDY[7,10] - lDZ[7,2] = 0
 lSC: ldy[7,10] - ldz[7,2] = 0
 lSC: lDY[7,15] - lDZ[7,3] = 0
 lSC: ldy[7,15] - ldz[7,3] = 0
 lSC: lDY[7,4] - lDZ[7,4] = 0
 lSC: ldy[7,4] - ldz[7,4] = 0
 lSC: lDY[7,9] - lDZ[7,5] = 0
 lSC: ldy[7,9] - ldz[7,5] = 0
 lSC: lDY[7,14] - lDZ[7,6] = 0
 lSC: ldy[7,14] - ldz[7,6] = 0
 lSC: lDY[7,3] - lDZ[7,7] = 0
 lSC: ldy[7,3] - ldz[7,7] = 0
 lSC: lDY[7,8] - lDZ[7,8] = 0
 lSC: ldy[7,8] - ldz[7,8] = 0
 lSC: lDY[7,13] - lDZ[7,9] = 0
 lSC: ldy[7,13] - ldz[7,9] = 0
 lSC: lDY[7,2] - lDZ[7,10] = 0
 lSC: ldy[7,2] - ldz[7,10] = 0
 lSC: lDY[7,7] - lDZ[7,11] = 0
 lSC: ldy[7,7] - ldz[7,11] = 0
 lSC: lDY[7,12] - lDZ[7,12] = 0
 lSC: ldy[7,12] - ldz[7,12] = 0
 lSC: lDY[7,1] - lDZ[7,13] = 0
 lSC: ldy[7,1] - ldz[7,13] = 0
 lSC: lDY[7,6] - lDZ[7,14] = 0
 lSC: ldy[7,6] - ldz[7,14] = 0
 lSC: lDY[7,11] - lDZ[7,15] = 0
 lSC: ldy[7,11] - ldz[7,15] = 0
 lMC: - ldw[1,0] + lTC[1,0] >= 0
 lMC: - ldw[1,1] + lTC[1,0] >= 0
 lMC: - ldw[1,2] + lTC[1,0] >= 0
 lMC: - ldw[1,3] + lTC[1,0] >= 0
 lMC: - ldw[1,0] - ldw[1,1] - ldw[1,2] - ldw[1,3] + lTC[1,0] <= 0
 lMC: - ldz[1,0] + lTC[1,0] <= 0
 lMC: - ldz[1,1] + lTC[1,0] <= 0
 lMC: - ldz[1,2] + lTC[1,0] <= 0
 lMC: - ldz[1,3] + lTC[1,0] <= 0
 lMC: - ldz[1,0] - ldz[1,1] - ldz[1,2] - ldz[1,3] + 4 lTC[1,0] >= 0
 lMC: - lDW[1,0] + lAC[1,0] >= 0
 lMC: - lDW[1,1] + lAC[1,0] >= 0
 lMC: - lDW[1,2] + lAC[1,0] >= 0
 lMC: - lDW[1,3] + lAC[1,0] >= 0
 lMC: - lDW[1,0] - lDW[1,1] - lDW[1,2] - lDW[1,3] + lAC[1,0] <= 0
 lMC: lDW[1,0] + lDW[1,1] + lDW[1,2] + lDW[1,3] + lDZ[1,0] + lDZ[1,1]
   + lDZ[1,2] + lDZ[1,3] - 5 lAC[1,0] >= 0
 lMC: lDW[1,0] + lDW[1,1] + lDW[1,2] + lDW[1,3] + lDZ[1,0] + lDZ[1,1]
   + lDZ[1,2] + lDZ[1,3] - 8 lAC[1,0] <= 0
 lMC: - ldw[1,4] + lTC[1,1] >= 0
 lMC: - ldw[1,5] + lTC[1,1] >= 0
 lMC: - ldw[1,6] + lTC[1,1] >= 0
 lMC: - ldw[1,7] + lTC[1,1] >= 0
 lMC: - ldw[1,4] - ldw[1,5] - ldw[1,6] - ldw[1,7] + lTC[1,1] <= 0
 lMC: - ldz[1,4] + lTC[1,1] <= 0
 lMC: - ldz[1,5] + lTC[1,1] <= 0
 lMC: - ldz[1,6] + lTC[1,1] <= 0
 lMC: - ldz[1,7] + lTC[1,1] <= 0
 lMC: - ldz[1,4] - ldz[1,5] - ldz[1,6] - ldz[1,7] + 4 lTC[1,1] >= 0
 lMC: - lDW[1,4] + lAC[1,1] >= 0
 lMC: - lDW[1,5] + lAC[1,1] >= 0
 lMC: - lDW[1,6] + lAC[1,1] >= 0
 lMC: - lDW[1,7] + lAC[1,1] >= 0
 lMC: - lDW[1,4] - lDW[1,5] - lDW[1,6] - lDW[1,7] + lAC[1,1] <= 0
 lMC: lDW[1,4] + lDW[1,5] + lDW[1,6] + lDW[1,7] + lDZ[1,4] + lDZ[1,5]
   + lDZ[1,6] + lDZ[1,7] - 5 lAC[1,1] >= 0
 lMC: lDW[1,4] + lDW[1,5] + lDW[1,6] + lDW[1,7] + lDZ[1,4] + lDZ[1,5]
   + lDZ[1,6] + lDZ[1,7] - 8 lAC[1,1] <= 0
 lMC: - ldw[1,8] + lTC[1,2] >= 0
 lMC: - ldw[1,9] + lTC[1,2] >= 0
 lMC: - ldw[1,10] + lTC[1,2] >= 0
 lMC: - ldw[1,11] + lTC[1,2] >= 0
 lMC: - ldw[1,8] - ldw[1,9] - ldw[1,10] - ldw[1,11] + lTC[1,2] <= 0
 lMC: - ldz[1,8] + lTC[1,2] <= 0
 lMC: - ldz[1,9] + lTC[1,2] <= 0
 lMC: - ldz[1,10] + lTC[1,2] <= 0
 lMC: - ldz[1,11] + lTC[1,2] <= 0
 lMC: - ldz[1,8] - ldz[1,9] - ldz[1,10] - ldz[1,11] + 4 lTC[1,2] >= 0
 lMC: - lDW[1,8] + lAC[1,2] >= 0
 lMC: - lDW[1,9] + lAC[1,2] >= 0
 lMC: - lDW[1,10] + lAC[1,2] >= 0
 lMC: - lDW[1,11] + lAC[1,2] >= 0
 lMC: - lDW[1,8] - lDW[1,9] - lDW[1,10] - lDW[1,11] + lAC[1,2] <= 0
 lMC: lDW[1,8] + lDW[1,9] + lDW[1,10] + lDW[1,11] + lDZ[1,8] + lDZ[1,9]
   + lDZ[1,10] + lDZ[1,11] - 5 lAC[1,2] >= 0
 lMC: lDW[1,8] + lDW[1,9] + lDW[1,10] + lDW[1,11] + lDZ[1,8] + lDZ[1,9]
   + lDZ[1,10] + lDZ[1,11] - 8 lAC[1,2] <= 0
 lMC: - ldw[1,12] + lTC[1,3] >= 0
 lMC: - ldw[1,13] + lTC[1,3] >= 0
 lMC: - ldw[1,14] + lTC[1,3] >= 0
 lMC: - ldw[1,15] + lTC[1,3] >= 0
 lMC: - ldw[1,12] - ldw[1,13] - ldw[1,14] - ldw[1,15] + lTC[1,3] <= 0
 lMC: - ldz[1,12] + lTC[1,3] <= 0
 lMC: - ldz[1,13] + lTC[1,3] <= 0
 lMC: - ldz[1,14] + lTC[1,3] <= 0
 lMC: - ldz[1,15] + lTC[1,3] <= 0
 lMC: - ldz[1,12] - ldz[1,13] - ldz[1,14] - ldz[1,15] + 4 lTC[1,3] >= 0
 lMC: - lDW[1,12] + lAC[1,3] >= 0
 lMC: - lDW[1,13] + lAC[1,3] >= 0
 lMC: - lDW[1,14] + lAC[1,3] >= 0
 lMC: - lDW[1,15] + lAC[1,3] >= 0
 lMC: - lDW[1,12] - lDW[1,13] - lDW[1,14] - lDW[1,15] + lAC[1,3] <= 0
 lMC: lDW[1,12] + lDW[1,13] + lDW[1,14] + lDW[1,15] + lDZ[1,12] + lDZ[1,13]
   + lDZ[1,14] + lDZ[1,15] - 5 lAC[1,3] >= 0
 lMC: lDW[1,12] + lDW[1,13] + lDW[1,14] + lDW[1,15] + lDZ[1,12] + lDZ[1,13]
   + lDZ[1,14] + lDZ[1,15] - 8 lAC[1,3] <= 0
 lMC: - ldw[2,0] + lTC[2,0] >= 0
 lMC: - ldw[2,1] + lTC[2,0] >= 0
 lMC: - ldw[2,2] + lTC[2,0] >= 0
 lMC: - ldw[2,3] + lTC[2,0] >= 0
 lMC: - ldw[2,0] - ldw[2,1] - ldw[2,2] - ldw[2,3] + lTC[2,0] <= 0
 lMC: - ldz[2,0] + lTC[2,0] <= 0
 lMC: - ldz[2,1] + lTC[2,0] <= 0
 lMC: - ldz[2,2] + lTC[2,0] <= 0
 lMC: - ldz[2,3] + lTC[2,0] <= 0
 lMC: - ldz[2,0] - ldz[2,1] - ldz[2,2] - ldz[2,3] + 4 lTC[2,0] >= 0
 lMC: - lDW[2,0] + lAC[2,0] >= 0
 lMC: - lDW[2,1] + lAC[2,0] >= 0
 lMC: - lDW[2,2] + lAC[2,0] >= 0
 lMC: - lDW[2,3] + lAC[2,0] >= 0
 lMC: - lDW[2,0] - lDW[2,1] - lDW[2,2] - lDW[2,3] + lAC[2,0] <= 0
 lMC: lDW[2,0] + lDW[2,1] + lDW[2,2] + lDW[2,3] + lDZ[2,0] + lDZ[2,1]
   + lDZ[2,2] + lDZ[2,3] - 5 lAC[2,0] >= 0
 lMC: lDW[2,0] + lDW[2,1] + lDW[2,2] + lDW[2,3] + lDZ[2,0] + lDZ[2,1]
   + lDZ[2,2] + lDZ[2,3] - 8 lAC[2,0] <= 0
 lMC: - ldw[2,4] + lTC[2,1] >= 0
 lMC: - ldw[2,5] + lTC[2,1] >= 0
 lMC: - ldw[2,6] + lTC[2,1] >= 0
 lMC: - ldw[2,7] + lTC[2,1] >= 0
 lMC: - ldw[2,4] - ldw[2,5] - ldw[2,6] - ldw[2,7] + lTC[2,1] <= 0
 lMC: - ldz[2,4] + lTC[2,1] <= 0
 lMC: - ldz[2,5] + lTC[2,1] <= 0
 lMC: - ldz[2,6] + lTC[2,1] <= 0
 lMC: - ldz[2,7] + lTC[2,1] <= 0
 lMC: - ldz[2,4] - ldz[2,5] - ldz[2,6] - ldz[2,7] + 4 lTC[2,1] >= 0
 lMC: - lDW[2,4] + lAC[2,1] >= 0
 lMC: - lDW[2,5] + lAC[2,1] >= 0
 lMC: - lDW[2,6] + lAC[2,1] >= 0
 lMC: - lDW[2,7] + lAC[2,1] >= 0
 lMC: - lDW[2,4] - lDW[2,5] - lDW[2,6] - lDW[2,7] + lAC[2,1] <= 0
 lMC: lDW[2,4] + lDW[2,5] + lDW[2,6] + lDW[2,7] + lDZ[2,4] + lDZ[2,5]
   + lDZ[2,6] + lDZ[2,7] - 5 lAC[2,1] >= 0
 lMC: lDW[2,4] + lDW[2,5] + lDW[2,6] + lDW[2,7] + lDZ[2,4] + lDZ[2,5]
   + lDZ[2,6] + lDZ[2,7] - 8 lAC[2,1] <= 0
 lMC: - ldw[2,8] + lTC[2,2] >= 0
 lMC: - ldw[2,9] + lTC[2,2] >= 0
 lMC: - ldw[2,10] + lTC[2,2] >= 0
 lMC: - ldw[2,11] + lTC[2,2] >= 0
 lMC: - ldw[2,8] - ldw[2,9] - ldw[2,10] - ldw[2,11] + lTC[2,2] <= 0
 lMC: - ldz[2,8] + lTC[2,2] <= 0
 lMC: - ldz[2,9] + lTC[2,2] <= 0
 lMC: - ldz[2,10] + lTC[2,2] <= 0
 lMC: - ldz[2,11] + lTC[2,2] <= 0
 lMC: - ldz[2,8] - ldz[2,9] - ldz[2,10] - ldz[2,11] + 4 lTC[2,2] >= 0
 lMC: - lDW[2,8] + lAC[2,2] >= 0
 lMC: - lDW[2,9] + lAC[2,2] >= 0
 lMC: - lDW[2,10] + lAC[2,2] >= 0
 lMC: - lDW[2,11] + lAC[2,2] >= 0
 lMC: - lDW[2,8] - lDW[2,9] - lDW[2,10] - lDW[2,11] + lAC[2,2] <= 0
 lMC: lDW[2,8] + lDW[2,9] + lDW[2,10] + lDW[2,11] + lDZ[2,8] + lDZ[2,9]
   + lDZ[2,10] + lDZ[2,11] - 5 lAC[2,2] >= 0
 lMC: lDW[2,8] + lDW[2,9] + lDW[2,10] + lDW[2,11] + lDZ[2,8] + lDZ[2,9]
   + lDZ[2,10] + lDZ[2,11] - 8 lAC[2,2] <= 0
 lMC: - ldw[2,12] + lTC[2,3] >= 0
 lMC: - ldw[2,13] + lTC[2,3] >= 0
 lMC: - ldw[2,14] + lTC[2,3] >= 0
 lMC: - ldw[2,15] + lTC[2,3] >= 0
 lMC: - ldw[2,12] - ldw[2,13] - ldw[2,14] - ldw[2,15] + lTC[2,3] <= 0
 lMC: - ldz[2,12] + lTC[2,3] <= 0
 lMC: - ldz[2,13] + lTC[2,3] <= 0
 lMC: - ldz[2,14] + lTC[2,3] <= 0
 lMC: - ldz[2,15] + lTC[2,3] <= 0
 lMC: - ldz[2,12] - ldz[2,13] - ldz[2,14] - ldz[2,15] + 4 lTC[2,3] >= 0
 lMC: - lDW[2,12] + lAC[2,3] >= 0
 lMC: - lDW[2,13] + lAC[2,3] >= 0
 lMC: - lDW[2,14] + lAC[2,3] >= 0
 lMC: - lDW[2,15] + lAC[2,3] >= 0
 lMC: - lDW[2,12] - lDW[2,13] - lDW[2,14] - lDW[2,15] + lAC[2,3] <= 0
 lMC: lDW[2,12] + lDW[2,13] + lDW[2,14] + lDW[2,15] + lDZ[2,12] + lDZ[2,13]
   + lDZ[2,14] + lDZ[2,15] - 5 lAC[2,3] >= 0
 lMC: lDW[2,12] + lDW[2,13] + lDW[2,14] + lDW[2,15] + lDZ[2,12] + lDZ[2,13]
   + lDZ[2,14] + lDZ[2,15] - 8 lAC[2,3] <= 0
 lMC: - ldw[3,0] + lTC[3,0] >= 0
 lMC: - ldw[3,1] + lTC[3,0] >= 0
 lMC: - ldw[3,2] + lTC[3,0] >= 0
 lMC: - ldw[3,3] + lTC[3,0] >= 0
 lMC: - ldw[3,0] - ldw[3,1] - ldw[3,2] - ldw[3,3] + lTC[3,0] <= 0
 lMC: - ldz[3,0] + lTC[3,0] <= 0
 lMC: - ldz[3,1] + lTC[3,0] <= 0
 lMC: - ldz[3,2] + lTC[3,0] <= 0
 lMC: - ldz[3,3] + lTC[3,0] <= 0
 lMC: - ldz[3,0] - ldz[3,1] - ldz[3,2] - ldz[3,3] + 4 lTC[3,0] >= 0
 lMC: - lDW[3,0] + lAC[3,0] >= 0
 lMC: - lDW[3,1] + lAC[3,0] >= 0
 lMC: - lDW[3,2] + lAC[3,0] >= 0
 lMC: - lDW[3,3] + lAC[3,0] >= 0
 lMC: - lDW[3,0] - lDW[3,1] - lDW[3,2] - lDW[3,3] + lAC[3,0] <= 0
 lMC: lDW[3,0] + lDW[3,1] + lDW[3,2] + lDW[3,3] + lDZ[3,0] + lDZ[3,1]
   + lDZ[3,2] + lDZ[3,3] - 5 lAC[3,0] >= 0
 lMC: lDW[3,0] + lDW[3,1] + lDW[3,2] + lDW[3,3] + lDZ[3,0] + lDZ[3,1]
   + lDZ[3,2] + lDZ[3,3] - 8 lAC[3,0] <= 0
 lMC: - ldw[3,4] + lTC[3,1] >= 0
 lMC: - ldw[3,5] + lTC[3,1] >= 0
 lMC: - ldw[3,6] + lTC[3,1] >= 0
 lMC: - ldw[3,7] + lTC[3,1] >= 0
 lMC: - ldw[3,4] - ldw[3,5] - ldw[3,6] - ldw[3,7] + lTC[3,1] <= 0
 lMC: - ldz[3,4] + lTC[3,1] <= 0
 lMC: - ldz[3,5] + lTC[3,1] <= 0
 lMC: - ldz[3,6] + lTC[3,1] <= 0
 lMC: - ldz[3,7] + lTC[3,1] <= 0
 lMC: - ldz[3,4] - ldz[3,5] - ldz[3,6] - ldz[3,7] + 4 lTC[3,1] >= 0
 lMC: - lDW[3,4] + lAC[3,1] >= 0
 lMC: - lDW[3,5] + lAC[3,1] >= 0
 lMC: - lDW[3,6] + lAC[3,1] >= 0
 lMC: - lDW[3,7] + lAC[3,1] >= 0
 lMC: - lDW[3,4] - lDW[3,5] - lDW[3,6] - lDW[3,7] + lAC[3,1] <= 0
 lMC: lDW[3,4] + lDW[3,5] + lDW[3,6] + lDW[3,7] + lDZ[3,4] + lDZ[3,5]
   + lDZ[3,6] + lDZ[3,7] - 5 lAC[3,1] >= 0
 lMC: lDW[3,4] + lDW[3,5] + lDW[3,6] + lDW[3,7] + lDZ[3,4] + lDZ[3,5]
   + lDZ[3,6] + lDZ[3,7] - 8 lAC[3,1] <= 0
 lMC: - ldw[3,8] + lTC[3,2] >= 0
 lMC: - ldw[3,9] + lTC[3,2] >= 0
 lMC: - ldw[3,10] + lTC[3,2] >= 0
 lMC: - ldw[3,11] + lTC[3,2] >= 0
 lMC: - ldw[3,8] - ldw[3,9] - ldw[3,10] - ldw[3,11] + lTC[3,2] <= 0
 lMC: - ldz[3,8] + lTC[3,2] <= 0
 lMC: - ldz[3,9] + lTC[3,2] <= 0
 lMC: - ldz[3,10] + lTC[3,2] <= 0
 lMC: - ldz[3,11] + lTC[3,2] <= 0
 lMC: - ldz[3,8] - ldz[3,9] - ldz[3,10] - ldz[3,11] + 4 lTC[3,2] >= 0
 lMC: - lDW[3,8] + lAC[3,2] >= 0
 lMC: - lDW[3,9] + lAC[3,2] >= 0
 lMC: - lDW[3,10] + lAC[3,2] >= 0
 lMC: - lDW[3,11] + lAC[3,2] >= 0
 lMC: - lDW[3,8] - lDW[3,9] - lDW[3,10] - lDW[3,11] + lAC[3,2] <= 0
 lMC: lDW[3,8] + lDW[3,9] + lDW[3,10] + lDW[3,11] + lDZ[3,8] + lDZ[3,9]
   + lDZ[3,10] + lDZ[3,11] - 5 lAC[3,2] >= 0
 lMC: lDW[3,8] + lDW[3,9] + lDW[3,10] + lDW[3,11] + lDZ[3,8] + lDZ[3,9]
   + lDZ[3,10] + lDZ[3,11] - 8 lAC[3,2] <= 0
 lMC: - ldw[3,12] + lTC[3,3] >= 0
 lMC: - ldw[3,13] + lTC[3,3] >= 0
 lMC: - ldw[3,14] + lTC[3,3] >= 0
 lMC: - ldw[3,15] + lTC[3,3] >= 0
 lMC: - ldw[3,12] - ldw[3,13] - ldw[3,14] - ldw[3,15] + lTC[3,3] <= 0
 lMC: - ldz[3,12] + lTC[3,3] <= 0
 lMC: - ldz[3,13] + lTC[3,3] <= 0
 lMC: - ldz[3,14] + lTC[3,3] <= 0
 lMC: - ldz[3,15] + lTC[3,3] <= 0
 lMC: - ldz[3,12] - ldz[3,13] - ldz[3,14] - ldz[3,15] + 4 lTC[3,3] >= 0
 lMC: - lDW[3,12] + lAC[3,3] >= 0
 lMC: - lDW[3,13] + lAC[3,3] >= 0
 lMC: - lDW[3,14] + lAC[3,3] >= 0
 lMC: - lDW[3,15] + lAC[3,3] >= 0
 lMC: - lDW[3,12] - lDW[3,13] - lDW[3,14] - lDW[3,15] + lAC[3,3] <= 0
 lMC: lDW[3,12] + lDW[3,13] + lDW[3,14] + lDW[3,15] + lDZ[3,12] + lDZ[3,13]
   + lDZ[3,14] + lDZ[3,15] - 5 lAC[3,3] >= 0
 lMC: lDW[3,12] + lDW[3,13] + lDW[3,14] + lDW[3,15] + lDZ[3,12] + lDZ[3,13]
   + lDZ[3,14] + lDZ[3,15] - 8 lAC[3,3] <= 0
 lMC: - ldw[4,0] + lTC[4,0] >= 0
 lMC: - ldw[4,1] + lTC[4,0] >= 0
 lMC: - ldw[4,2] + lTC[4,0] >= 0
 lMC: - ldw[4,3] + lTC[4,0] >= 0
 lMC: - ldw[4,0] - ldw[4,1] - ldw[4,2] - ldw[4,3] + lTC[4,0] <= 0
 lMC: - ldz[4,0] + lTC[4,0] <= 0
 lMC: - ldz[4,1] + lTC[4,0] <= 0
 lMC: - ldz[4,2] + lTC[4,0] <= 0
 lMC: - ldz[4,3] + lTC[4,0] <= 0
 lMC: - ldz[4,0] - ldz[4,1] - ldz[4,2] - ldz[4,3] + 4 lTC[4,0] >= 0
 lMC: - lDW[4,0] + lAC[4,0] >= 0
 lMC: - lDW[4,1] + lAC[4,0] >= 0
 lMC: - lDW[4,2] + lAC[4,0] >= 0
 lMC: - lDW[4,3] + lAC[4,0] >= 0
 lMC: - lDW[4,0] - lDW[4,1] - lDW[4,2] - lDW[4,3] + lAC[4,0] <= 0
 lMC: lDW[4,0] + lDW[4,1] + lDW[4,2] + lDW[4,3] + lDZ[4,0] + lDZ[4,1]
   + lDZ[4,2] + lDZ[4,3] - 5 lAC[4,0] >= 0
 lMC: lDW[4,0] + lDW[4,1] + lDW[4,2] + lDW[4,3] + lDZ[4,0] + lDZ[4,1]
   + lDZ[4,2] + lDZ[4,3] - 8 lAC[4,0] <= 0
 lMC: - ldw[4,4] + lTC[4,1] >= 0
 lMC: - ldw[4,5] + lTC[4,1] >= 0
 lMC: - ldw[4,6] + lTC[4,1] >= 0
 lMC: - ldw[4,7] + lTC[4,1] >= 0
 lMC: - ldw[4,4] - ldw[4,5] - ldw[4,6] - ldw[4,7] + lTC[4,1] <= 0
 lMC: - ldz[4,4] + lTC[4,1] <= 0
 lMC: - ldz[4,5] + lTC[4,1] <= 0
 lMC: - ldz[4,6] + lTC[4,1] <= 0
 lMC: - ldz[4,7] + lTC[4,1] <= 0
 lMC: - ldz[4,4] - ldz[4,5] - ldz[4,6] - ldz[4,7] + 4 lTC[4,1] >= 0
 lMC: - lDW[4,4] + lAC[4,1] >= 0
 lMC: - lDW[4,5] + lAC[4,1] >= 0
 lMC: - lDW[4,6] + lAC[4,1] >= 0
 lMC: - lDW[4,7] + lAC[4,1] >= 0
 lMC: - lDW[4,4] - lDW[4,5] - lDW[4,6] - lDW[4,7] + lAC[4,1] <= 0
 lMC: lDW[4,4] + lDW[4,5] + lDW[4,6] + lDW[4,7] + lDZ[4,4] + lDZ[4,5]
   + lDZ[4,6] + lDZ[4,7] - 5 lAC[4,1] >= 0
 lMC: lDW[4,4] + lDW[4,5] + lDW[4,6] + lDW[4,7] + lDZ[4,4] + lDZ[4,5]
   + lDZ[4,6] + lDZ[4,7] - 8 lAC[4,1] <= 0
 lMC: - ldw[4,8] + lTC[4,2] >= 0
 lMC: - ldw[4,9] + lTC[4,2] >= 0
 lMC: - ldw[4,10] + lTC[4,2] >= 0
 lMC: - ldw[4,11] + lTC[4,2] >= 0
 lMC: - ldw[4,8] - ldw[4,9] - ldw[4,10] - ldw[4,11] + lTC[4,2] <= 0
 lMC: - ldz[4,8] + lTC[4,2] <= 0
 lMC: - ldz[4,9] + lTC[4,2] <= 0
 lMC: - ldz[4,10] + lTC[4,2] <= 0
 lMC: - ldz[4,11] + lTC[4,2] <= 0
 lMC: - ldz[4,8] - ldz[4,9] - ldz[4,10] - ldz[4,11] + 4 lTC[4,2] >= 0
 lMC: - lDW[4,8] + lAC[4,2] >= 0
 lMC: - lDW[4,9] + lAC[4,2] >= 0
 lMC: - lDW[4,10] + lAC[4,2] >= 0
 lMC: - lDW[4,11] + lAC[4,2] >= 0
 lMC: - lDW[4,8] - lDW[4,9] - lDW[4,10] - lDW[4,11] + lAC[4,2] <= 0
 lMC: lDW[4,8] + lDW[4,9] + lDW[4,10] + lDW[4,11] + lDZ[4,8] + lDZ[4,9]
   + lDZ[4,10] + lDZ[4,11] - 5 lAC[4,2] >= 0
 lMC: lDW[4,8] + lDW[4,9] + lDW[4,10] + lDW[4,11] + lDZ[4,8] + lDZ[4,9]
   + lDZ[4,10] + lDZ[4,11] - 8 lAC[4,2] <= 0
 lMC: - ldw[4,12] + lTC[4,3] >= 0
 lMC: - ldw[4,13] + lTC[4,3] >= 0
 lMC: - ldw[4,14] + lTC[4,3] >= 0
 lMC: - ldw[4,15] + lTC[4,3] >= 0
 lMC: - ldw[4,12] - ldw[4,13] - ldw[4,14] - ldw[4,15] + lTC[4,3] <= 0
 lMC: - ldz[4,12] + lTC[4,3] <= 0
 lMC: - ldz[4,13] + lTC[4,3] <= 0
 lMC: - ldz[4,14] + lTC[4,3] <= 0
 lMC: - ldz[4,15] + lTC[4,3] <= 0
 lMC: - ldz[4,12] - ldz[4,13] - ldz[4,14] - ldz[4,15] + 4 lTC[4,3] >= 0
 lMC: - lDW[4,12] + lAC[4,3] >= 0
 lMC: - lDW[4,13] + lAC[4,3] >= 0
 lMC: - lDW[4,14] + lAC[4,3] >= 0
 lMC: - lDW[4,15] + lAC[4,3] >= 0
 lMC: - lDW[4,12] - lDW[4,13] - lDW[4,14] - lDW[4,15] + lAC[4,3] <= 0
 lMC: lDW[4,12] + lDW[4,13] + lDW[4,14] + lDW[4,15] + lDZ[4,12] + lDZ[4,13]
   + lDZ[4,14] + lDZ[4,15] - 5 lAC[4,3] >= 0
 lMC: lDW[4,12] + lDW[4,13] + lDW[4,14] + lDW[4,15] + lDZ[4,12] + lDZ[4,13]
   + lDZ[4,14] + lDZ[4,15] - 8 lAC[4,3] <= 0
 lMC: - ldw[5,0] + lTC[5,0] >= 0
 lMC: - ldw[5,1] + lTC[5,0] >= 0
 lMC: - ldw[5,2] + lTC[5,0] >= 0
 lMC: - ldw[5,3] + lTC[5,0] >= 0
 lMC: - ldw[5,0] - ldw[5,1] - ldw[5,2] - ldw[5,3] + lTC[5,0] <= 0
 lMC: - ldz[5,0] + lTC[5,0] <= 0
 lMC: - ldz[5,1] + lTC[5,0] <= 0
 lMC: - ldz[5,2] + lTC[5,0] <= 0
 lMC: - ldz[5,3] + lTC[5,0] <= 0
 lMC: - ldz[5,0] - ldz[5,1] - ldz[5,2] - ldz[5,3] + 4 lTC[5,0] >= 0
 lMC: - lDW[5,0] + lAC[5,0] >= 0
 lMC: - lDW[5,1] + lAC[5,0] >= 0
 lMC: - lDW[5,2] + lAC[5,0] >= 0
 lMC: - lDW[5,3] + lAC[5,0] >= 0
 lMC: - lDW[5,0] - lDW[5,1] - lDW[5,2] - lDW[5,3] + lAC[5,0] <= 0
 lMC: lDW[5,0] + lDW[5,1] + lDW[5,2] + lDW[5,3] + lDZ[5,0] + lDZ[5,1]
   + lDZ[5,2] + lDZ[5,3] - 5 lAC[5,0] >= 0
 lMC: lDW[5,0] + lDW[5,1] + lDW[5,2] + lDW[5,3] + lDZ[5,0] + lDZ[5,1]
   + lDZ[5,2] + lDZ[5,3] - 8 lAC[5,0] <= 0
 lMC: - ldw[5,4] + lTC[5,1] >= 0
 lMC: - ldw[5,5] + lTC[5,1] >= 0
 lMC: - ldw[5,6] + lTC[5,1] >= 0
 lMC: - ldw[5,7] + lTC[5,1] >= 0
 lMC: - ldw[5,4] - ldw[5,5] - ldw[5,6] - ldw[5,7] + lTC[5,1] <= 0
 lMC: - ldz[5,4] + lTC[5,1] <= 0
 lMC: - ldz[5,5] + lTC[5,1] <= 0
 lMC: - ldz[5,6] + lTC[5,1] <= 0
 lMC: - ldz[5,7] + lTC[5,1] <= 0
 lMC: - ldz[5,4] - ldz[5,5] - ldz[5,6] - ldz[5,7] + 4 lTC[5,1] >= 0
 lMC: - lDW[5,4] + lAC[5,1] >= 0
 lMC: - lDW[5,5] + lAC[5,1] >= 0
 lMC: - lDW[5,6] + lAC[5,1] >= 0
 lMC: - lDW[5,7] + lAC[5,1] >= 0
 lMC: - lDW[5,4] - lDW[5,5] - lDW[5,6] - lDW[5,7] + lAC[5,1] <= 0
 lMC: lDW[5,4] + lDW[5,5] + lDW[5,6] + lDW[5,7] + lDZ[5,4] + lDZ[5,5]
   + lDZ[5,6] + lDZ[5,7] - 5 lAC[5,1] >= 0
 lMC: lDW[5,4] + lDW[5,5] + lDW[5,6] + lDW[5,7] + lDZ[5,4] + lDZ[5,5]
   + lDZ[5,6] + lDZ[5,7] - 8 lAC[5,1] <= 0
 lMC: - ldw[5,8] + lTC[5,2] >= 0
 lMC: - ldw[5,9] + lTC[5,2] >= 0
 lMC: - ldw[5,10] + lTC[5,2] >= 0
 lMC: - ldw[5,11] + lTC[5,2] >= 0
 lMC: - ldw[5,8] - ldw[5,9] - ldw[5,10] - ldw[5,11] + lTC[5,2] <= 0
 lMC: - ldz[5,8] + lTC[5,2] <= 0
 lMC: - ldz[5,9] + lTC[5,2] <= 0
 lMC: - ldz[5,10] + lTC[5,2] <= 0
 lMC: - ldz[5,11] + lTC[5,2] <= 0
 lMC: - ldz[5,8] - ldz[5,9] - ldz[5,10] - ldz[5,11] + 4 lTC[5,2] >= 0
 lMC: - lDW[5,8] + lAC[5,2] >= 0
 lMC: - lDW[5,9] + lAC[5,2] >= 0
 lMC: - lDW[5,10] + lAC[5,2] >= 0
 lMC: - lDW[5,11] + lAC[5,2] >= 0
 lMC: - lDW[5,8] - lDW[5,9] - lDW[5,10] - lDW[5,11] + lAC[5,2] <= 0
 lMC: lDW[5,8] + lDW[5,9] + lDW[5,10] + lDW[5,11] + lDZ[5,8] + lDZ[5,9]
   + lDZ[5,10] + lDZ[5,11] - 5 lAC[5,2] >= 0
 lMC: lDW[5,8] + lDW[5,9] + lDW[5,10] + lDW[5,11] + lDZ[5,8] + lDZ[5,9]
   + lDZ[5,10] + lDZ[5,11] - 8 lAC[5,2] <= 0
 lMC: - ldw[5,12] + lTC[5,3] >= 0
 lMC: - ldw[5,13] + lTC[5,3] >= 0
 lMC: - ldw[5,14] + lTC[5,3] >= 0
 lMC: - ldw[5,15] + lTC[5,3] >= 0
 lMC: - ldw[5,12] - ldw[5,13] - ldw[5,14] - ldw[5,15] + lTC[5,3] <= 0
 lMC: - ldz[5,12] + lTC[5,3] <= 0
 lMC: - ldz[5,13] + lTC[5,3] <= 0
 lMC: - ldz[5,14] + lTC[5,3] <= 0
 lMC: - ldz[5,15] + lTC[5,3] <= 0
 lMC: - ldz[5,12] - ldz[5,13] - ldz[5,14] - ldz[5,15] + 4 lTC[5,3] >= 0
 lMC: - lDW[5,12] + lAC[5,3] >= 0
 lMC: - lDW[5,13] + lAC[5,3] >= 0
 lMC: - lDW[5,14] + lAC[5,3] >= 0
 lMC: - lDW[5,15] + lAC[5,3] >= 0
 lMC: - lDW[5,12] - lDW[5,13] - lDW[5,14] - lDW[5,15] + lAC[5,3] <= 0
 lMC: lDW[5,12] + lDW[5,13] + lDW[5,14] + lDW[5,15] + lDZ[5,12] + lDZ[5,13]
   + lDZ[5,14] + lDZ[5,15] - 5 lAC[5,3] >= 0
 lMC: lDW[5,12] + lDW[5,13] + lDW[5,14] + lDW[5,15] + lDZ[5,12] + lDZ[5,13]
   + lDZ[5,14] + lDZ[5,15] - 8 lAC[5,3] <= 0
 lMC: - ldw[6,0] + lTC[6,0] >= 0
 lMC: - ldw[6,1] + lTC[6,0] >= 0
 lMC: - ldw[6,2] + lTC[6,0] >= 0
 lMC: - ldw[6,3] + lTC[6,0] >= 0
 lMC: - ldw[6,0] - ldw[6,1] - ldw[6,2] - ldw[6,3] + lTC[6,0] <= 0
 lMC: - ldz[6,0] + lTC[6,0] <= 0
 lMC: - ldz[6,1] + lTC[6,0] <= 0
 lMC: - ldz[6,2] + lTC[6,0] <= 0
 lMC: - ldz[6,3] + lTC[6,0] <= 0
 lMC: - ldz[6,0] - ldz[6,1] - ldz[6,2] - ldz[6,3] + 4 lTC[6,0] >= 0
 lMC: - lDW[6,0] + lAC[6,0] >= 0
 lMC: - lDW[6,1] + lAC[6,0] >= 0
 lMC: - lDW[6,2] + lAC[6,0] >= 0
 lMC: - lDW[6,3] + lAC[6,0] >= 0
 lMC: - lDW[6,0] - lDW[6,1] - lDW[6,2] - lDW[6,3] + lAC[6,0] <= 0
 lMC: lDW[6,0] + lDW[6,1] + lDW[6,2] + lDW[6,3] + lDZ[6,0] + lDZ[6,1]
   + lDZ[6,2] + lDZ[6,3] - 5 lAC[6,0] >= 0
 lMC: lDW[6,0] + lDW[6,1] + lDW[6,2] + lDW[6,3] + lDZ[6,0] + lDZ[6,1]
   + lDZ[6,2] + lDZ[6,3] - 8 lAC[6,0] <= 0
 lMC: - ldw[6,4] + lTC[6,1] >= 0
 lMC: - ldw[6,5] + lTC[6,1] >= 0
 lMC: - ldw[6,6] + lTC[6,1] >= 0
 lMC: - ldw[6,7] + lTC[6,1] >= 0
 lMC: - ldw[6,4] - ldw[6,5] - ldw[6,6] - ldw[6,7] + lTC[6,1] <= 0
 lMC: - ldz[6,4] + lTC[6,1] <= 0
 lMC: - ldz[6,5] + lTC[6,1] <= 0
 lMC: - ldz[6,6] + lTC[6,1] <= 0
 lMC: - ldz[6,7] + lTC[6,1] <= 0
 lMC: - ldz[6,4] - ldz[6,5] - ldz[6,6] - ldz[6,7] + 4 lTC[6,1] >= 0
 lMC: - lDW[6,4] + lAC[6,1] >= 0
 lMC: - lDW[6,5] + lAC[6,1] >= 0
 lMC: - lDW[6,6] + lAC[6,1] >= 0
 lMC: - lDW[6,7] + lAC[6,1] >= 0
 lMC: - lDW[6,4] - lDW[6,5] - lDW[6,6] - lDW[6,7] + lAC[6,1] <= 0
 lMC: lDW[6,4] + lDW[6,5] + lDW[6,6] + lDW[6,7] + lDZ[6,4] + lDZ[6,5]
   + lDZ[6,6] + lDZ[6,7] - 5 lAC[6,1] >= 0
 lMC: lDW[6,4] + lDW[6,5] + lDW[6,6] + lDW[6,7] + lDZ[6,4] + lDZ[6,5]
   + lDZ[6,6] + lDZ[6,7] - 8 lAC[6,1] <= 0
 lMC: - ldw[6,8] + lTC[6,2] >= 0
 lMC: - ldw[6,9] + lTC[6,2] >= 0
 lMC: - ldw[6,10] + lTC[6,2] >= 0
 lMC: - ldw[6,11] + lTC[6,2] >= 0
 lMC: - ldw[6,8] - ldw[6,9] - ldw[6,10] - ldw[6,11] + lTC[6,2] <= 0
 lMC: - ldz[6,8] + lTC[6,2] <= 0
 lMC: - ldz[6,9] + lTC[6,2] <= 0
 lMC: - ldz[6,10] + lTC[6,2] <= 0
 lMC: - ldz[6,11] + lTC[6,2] <= 0
 lMC: - ldz[6,8] - ldz[6,9] - ldz[6,10] - ldz[6,11] + 4 lTC[6,2] >= 0
 lMC: - lDW[6,8] + lAC[6,2] >= 0
 lMC: - lDW[6,9] + lAC[6,2] >= 0
 lMC: - lDW[6,10] + lAC[6,2] >= 0
 lMC: - lDW[6,11] + lAC[6,2] >= 0
 lMC: - lDW[6,8] - lDW[6,9] - lDW[6,10] - lDW[6,11] + lAC[6,2] <= 0
 lMC: lDW[6,8] + lDW[6,9] + lDW[6,10] + lDW[6,11] + lDZ[6,8] + lDZ[6,9]
   + lDZ[6,10] + lDZ[6,11] - 5 lAC[6,2] >= 0
 lMC: lDW[6,8] + lDW[6,9] + lDW[6,10] + lDW[6,11] + lDZ[6,8] + lDZ[6,9]
   + lDZ[6,10] + lDZ[6,11] - 8 lAC[6,2] <= 0
 lMC: - ldw[6,12] + lTC[6,3] >= 0
 lMC: - ldw[6,13] + lTC[6,3] >= 0
 lMC: - ldw[6,14] + lTC[6,3] >= 0
 lMC: - ldw[6,15] + lTC[6,3] >= 0
 lMC: - ldw[6,12] - ldw[6,13] - ldw[6,14] - ldw[6,15] + lTC[6,3] <= 0
 lMC: - ldz[6,12] + lTC[6,3] <= 0
 lMC: - ldz[6,13] + lTC[6,3] <= 0
 lMC: - ldz[6,14] + lTC[6,3] <= 0
 lMC: - ldz[6,15] + lTC[6,3] <= 0
 lMC: - ldz[6,12] - ldz[6,13] - ldz[6,14] - ldz[6,15] + 4 lTC[6,3] >= 0
 lMC: - lDW[6,12] + lAC[6,3] >= 0
 lMC: - lDW[6,13] + lAC[6,3] >= 0
 lMC: - lDW[6,14] + lAC[6,3] >= 0
 lMC: - lDW[6,15] + lAC[6,3] >= 0
 lMC: - lDW[6,12] - lDW[6,13] - lDW[6,14] - lDW[6,15] + lAC[6,3] <= 0
 lMC: lDW[6,12] + lDW[6,13] + lDW[6,14] + lDW[6,15] + lDZ[6,12] + lDZ[6,13]
   + lDZ[6,14] + lDZ[6,15] - 5 lAC[6,3] >= 0
 lMC: lDW[6,12] + lDW[6,13] + lDW[6,14] + lDW[6,15] + lDZ[6,12] + lDZ[6,13]
   + lDZ[6,14] + lDZ[6,15] - 8 lAC[6,3] <= 0
 lMC: - ldw[7,0] + lTC[7,0] >= 0
 lMC: - ldw[7,1] + lTC[7,0] >= 0
 lMC: - ldw[7,2] + lTC[7,0] >= 0
 lMC: - ldw[7,3] + lTC[7,0] >= 0
 lMC: - ldw[7,0] - ldw[7,1] - ldw[7,2] - ldw[7,3] + lTC[7,0] <= 0
 lMC: - ldz[7,0] + lTC[7,0] <= 0
 lMC: - ldz[7,1] + lTC[7,0] <= 0
 lMC: - ldz[7,2] + lTC[7,0] <= 0
 lMC: - ldz[7,3] + lTC[7,0] <= 0
 lMC: - ldz[7,0] - ldz[7,1] - ldz[7,2] - ldz[7,3] + 4 lTC[7,0] >= 0
 lMC: - lDW[7,0] + lAC[7,0] >= 0
 lMC: - lDW[7,1] + lAC[7,0] >= 0
 lMC: - lDW[7,2] + lAC[7,0] >= 0
 lMC: - lDW[7,3] + lAC[7,0] >= 0
 lMC: - lDW[7,0] - lDW[7,1] - lDW[7,2] - lDW[7,3] + lAC[7,0] <= 0
 lMC: lDW[7,0] + lDW[7,1] + lDW[7,2] + lDW[7,3] + lDZ[7,0] + lDZ[7,1]
   + lDZ[7,2] + lDZ[7,3] - 5 lAC[7,0] >= 0
 lMC: lDW[7,0] + lDW[7,1] + lDW[7,2] + lDW[7,3] + lDZ[7,0] + lDZ[7,1]
   + lDZ[7,2] + lDZ[7,3] - 8 lAC[7,0] <= 0
 lMC: - ldw[7,4] + lTC[7,1] >= 0
 lMC: - ldw[7,5] + lTC[7,1] >= 0
 lMC: - ldw[7,6] + lTC[7,1] >= 0
 lMC: - ldw[7,7] + lTC[7,1] >= 0
 lMC: - ldw[7,4] - ldw[7,5] - ldw[7,6] - ldw[7,7] + lTC[7,1] <= 0
 lMC: - ldz[7,4] + lTC[7,1] <= 0
 lMC: - ldz[7,5] + lTC[7,1] <= 0
 lMC: - ldz[7,6] + lTC[7,1] <= 0
 lMC: - ldz[7,7] + lTC[7,1] <= 0
 lMC: - ldz[7,4] - ldz[7,5] - ldz[7,6] - ldz[7,7] + 4 lTC[7,1] >= 0
 lMC: - lDW[7,4] + lAC[7,1] >= 0
 lMC: - lDW[7,5] + lAC[7,1] >= 0
 lMC: - lDW[7,6] + lAC[7,1] >= 0
 lMC: - lDW[7,7] + lAC[7,1] >= 0
 lMC: - lDW[7,4] - lDW[7,5] - lDW[7,6] - lDW[7,7] + lAC[7,1] <= 0
 lMC: lDW[7,4] + lDW[7,5] + lDW[7,6] + lDW[7,7] + lDZ[7,4] + lDZ[7,5]
   + lDZ[7,6] + lDZ[7,7] - 5 lAC[7,1] >= 0
 lMC: lDW[7,4] + lDW[7,5] + lDW[7,6] + lDW[7,7] + lDZ[7,4] + lDZ[7,5]
   + lDZ[7,6] + lDZ[7,7] - 8 lAC[7,1] <= 0
 lMC: - ldw[7,8] + lTC[7,2] >= 0
 lMC: - ldw[7,9] + lTC[7,2] >= 0
 lMC: - ldw[7,10] + lTC[7,2] >= 0
 lMC: - ldw[7,11] + lTC[7,2] >= 0
 lMC: - ldw[7,8] - ldw[7,9] - ldw[7,10] - ldw[7,11] + lTC[7,2] <= 0
 lMC: - ldz[7,8] + lTC[7,2] <= 0
 lMC: - ldz[7,9] + lTC[7,2] <= 0
 lMC: - ldz[7,10] + lTC[7,2] <= 0
 lMC: - ldz[7,11] + lTC[7,2] <= 0
 lMC: - ldz[7,8] - ldz[7,9] - ldz[7,10] - ldz[7,11] + 4 lTC[7,2] >= 0
 lMC: - lDW[7,8] + lAC[7,2] >= 0
 lMC: - lDW[7,9] + lAC[7,2] >= 0
 lMC: - lDW[7,10] + lAC[7,2] >= 0
 lMC: - lDW[7,11] + lAC[7,2] >= 0
 lMC: - lDW[7,8] - lDW[7,9] - lDW[7,10] - lDW[7,11] + lAC[7,2] <= 0
 lMC: lDW[7,8] + lDW[7,9] + lDW[7,10] + lDW[7,11] + lDZ[7,8] + lDZ[7,9]
   + lDZ[7,10] + lDZ[7,11] - 5 lAC[7,2] >= 0
 lMC: lDW[7,8] + lDW[7,9] + lDW[7,10] + lDW[7,11] + lDZ[7,8] + lDZ[7,9]
   + lDZ[7,10] + lDZ[7,11] - 8 lAC[7,2] <= 0
 lMC: - ldw[7,12] + lTC[7,3] >= 0
 lMC: - ldw[7,13] + lTC[7,3] >= 0
 lMC: - ldw[7,14] + lTC[7,3] >= 0
 lMC: - ldw[7,15] + lTC[7,3] >= 0
 lMC: - ldw[7,12] - ldw[7,13] - ldw[7,14] - ldw[7,15] + lTC[7,3] <= 0
 lMC: - ldz[7,12] + lTC[7,3] <= 0
 lMC: - ldz[7,13] + lTC[7,3] <= 0
 lMC: - ldz[7,14] + lTC[7,3] <= 0
 lMC: - ldz[7,15] + lTC[7,3] <= 0
 lMC: - ldz[7,12] - ldz[7,13] - ldz[7,14] - ldz[7,15] + 4 lTC[7,3] >= 0
 lMC: - lDW[7,12] + lAC[7,3] >= 0
 lMC: - lDW[7,13] + lAC[7,3] >= 0
 lMC: - lDW[7,14] + lAC[7,3] >= 0
 lMC: - lDW[7,15] + lAC[7,3] >= 0
 lMC: - lDW[7,12] - lDW[7,13] - lDW[7,14] - lDW[7,15] + lAC[7,3] <= 0
 lMC: lDW[7,12] + lDW[7,13] + lDW[7,14] + lDW[7,15] + lDZ[7,12] + lDZ[7,13]
   + lDZ[7,14] + lDZ[7,15] - 5 lAC[7,3] >= 0
 lMC: lDW[7,12] + lDW[7,13] + lDW[7,14] + lDW[7,15] + lDZ[7,12] + lDZ[7,13]
   + lDZ[7,14] + lDZ[7,15] - 8 lAC[7,3] <= 0
 lBaseq: leqDY[8,0] - leqdy[8,0] >= 0
 lBaseq: leqDZ[8,0] - leqdz[8,0] >= 0
 lBaseq: leqDY[9,0] - leqdy[9,0] >= 0
 lBaseq: leqDZ[9,0] - leqdz[9,0] >= 0
 lBaseq: leqDX[9,0] - leqdx[9,0] >= 0
 lBaseq: leqDW[8,0] - leqdw[8,0] >= 0
 lBaseq: leqDY[8,1] - leqdy[8,1] >= 0
 lBaseq: leqDZ[8,1] - leqdz[8,1] >= 0
 lBaseq: leqDY[9,1] - leqdy[9,1] >= 0
 lBaseq: leqDZ[9,1] - leqdz[9,1] >= 0
 lBaseq: leqDX[9,1] - leqdx[9,1] >= 0
 lBaseq: leqDW[8,1] - leqdw[8,1] >= 0
 lBaseq: leqDY[8,2] - leqdy[8,2] >= 0
 lBaseq: leqDZ[8,2] - leqdz[8,2] >= 0
 lBaseq: leqDY[9,2] - leqdy[9,2] >= 0
 lBaseq: leqDZ[9,2] - leqdz[9,2] >= 0
 lBaseq: leqDX[9,2] - leqdx[9,2] >= 0
 lBaseq: leqDW[8,2] - leqdw[8,2] >= 0
 lBaseq: leqDY[8,3] - leqdy[8,3] >= 0
 lBaseq: leqDZ[8,3] - leqdz[8,3] >= 0
 lBaseq: leqDY[9,3] - leqdy[9,3] >= 0
 lBaseq: leqDZ[9,3] - leqdz[9,3] >= 0
 lBaseq: leqDX[9,3] - leqdx[9,3] >= 0
 lBaseq: leqDW[8,3] - leqdw[8,3] >= 0
 lBaseq: leqDY[8,4] - leqdy[8,4] >= 0
 lBaseq: leqDZ[8,4] - leqdz[8,4] >= 0
 lBaseq: leqDY[9,4] - leqdy[9,4] >= 0
 lBaseq: leqDZ[9,4] - leqdz[9,4] >= 0
 lBaseq: leqDX[9,4] - leqdx[9,4] >= 0
 lBaseq: leqDW[8,4] - leqdw[8,4] >= 0
 lBaseq: leqDY[8,5] - leqdy[8,5] >= 0
 lBaseq: leqDZ[8,5] - leqdz[8,5] >= 0
 lBaseq: leqDY[9,5] - leqdy[9,5] >= 0
 lBaseq: leqDZ[9,5] - leqdz[9,5] >= 0
 lBaseq: leqDX[9,5] - leqdx[9,5] >= 0
 lBaseq: leqDW[8,5] - leqdw[8,5] >= 0
 lBaseq: leqDY[8,6] - leqdy[8,6] >= 0
 lBaseq: leqDZ[8,6] - leqdz[8,6] >= 0
 lBaseq: leqDY[9,6] - leqdy[9,6] >= 0
 lBaseq: leqDZ[9,6] - leqdz[9,6] >= 0
 lBaseq: leqDX[9,6] - leqdx[9,6] >= 0
 lBaseq: leqDW[8,6] - leqdw[8,6] >= 0
 lBaseq: leqDY[8,7] - leqdy[8,7] >= 0
 lBaseq: leqDZ[8,7] - leqdz[8,7] >= 0
 lBaseq: leqDY[9,7] - leqdy[9,7] >= 0
 lBaseq: leqDZ[9,7] - leqdz[9,7] >= 0
 lBaseq: leqDX[9,7] - leqdx[9,7] >= 0
 lBaseq: leqDW[8,7] - leqdw[8,7] >= 0
 lBaseq: leqDY[8,8] - leqdy[8,8] >= 0
 lBaseq: leqDZ[8,8] - leqdz[8,8] >= 0
 lBaseq: leqDY[9,8] - leqdy[9,8] >= 0
 lBaseq: leqDZ[9,8] - leqdz[9,8] >= 0
 lBaseq: leqDX[9,8] - leqdx[9,8] >= 0
 lBaseq: leqDW[8,8] - leqdw[8,8] >= 0
 lBaseq: leqDY[8,9] - leqdy[8,9] >= 0
 lBaseq: leqDZ[8,9] - leqdz[8,9] >= 0
 lBaseq: leqDY[9,9] - leqdy[9,9] >= 0
 lBaseq: leqDZ[9,9] - leqdz[9,9] >= 0
 lBaseq: leqDX[9,9] - leqdx[9,9] >= 0
 lBaseq: leqDW[8,9] - leqdw[8,9] >= 0
 lBaseq: leqDY[8,10] - leqdy[8,10] >= 0
 lBaseq: leqDZ[8,10] - leqdz[8,10] >= 0
 lBaseq: leqDY[9,10] - leqdy[9,10] >= 0
 lBaseq: leqDZ[9,10] - leqdz[9,10] >= 0
 lBaseq: leqDX[9,10] - leqdx[9,10] >= 0
 lBaseq: leqDW[8,10] - leqdw[8,10] >= 0
 lBaseq: leqDY[8,11] - leqdy[8,11] >= 0
 lBaseq: leqDZ[8,11] - leqdz[8,11] >= 0
 lBaseq: leqDY[9,11] - leqdy[9,11] >= 0
 lBaseq: leqDZ[9,11] - leqdz[9,11] >= 0
 lBaseq: leqDX[9,11] - leqdx[9,11] >= 0
 lBaseq: leqDW[8,11] - leqdw[8,11] >= 0
 lBaseq: leqDY[8,12] - leqdy[8,12] >= 0
 lBaseq: leqDZ[8,12] - leqdz[8,12] >= 0
 lBaseq: leqDY[9,12] - leqdy[9,12] >= 0
 lBaseq: leqDZ[9,12] - leqdz[9,12] >= 0
 lBaseq: leqDX[9,12] - leqdx[9,12] >= 0
 lBaseq: leqDW[8,12] - leqdw[8,12] >= 0
 lBaseq: leqDY[8,13] - leqdy[8,13] >= 0
 lBaseq: leqDZ[8,13] - leqdz[8,13] >= 0
 lBaseq: leqDY[9,13] - leqdy[9,13] >= 0
 lBaseq: leqDZ[9,13] - leqdz[9,13] >= 0
 lBaseq: leqDX[9,13] - leqdx[9,13] >= 0
 lBaseq: leqDW[8,13] - leqdw[8,13] >= 0
 lBaseq: leqDY[8,14] - leqdy[8,14] >= 0
 lBaseq: leqDZ[8,14] - leqdz[8,14] >= 0
 lBaseq: leqDY[9,14] - leqdy[9,14] >= 0
 lBaseq: leqDZ[9,14] - leqdz[9,14] >= 0
 lBaseq: leqDX[9,14] - leqdx[9,14] >= 0
 lBaseq: leqDW[8,14] - leqdw[8,14] >= 0
 lBaseq: leqDY[8,15] - leqdy[8,15] >= 0
 lBaseq: leqDZ[8,15] - leqdz[8,15] >= 0
 lBaseq: leqDY[9,15] - leqdy[9,15] >= 0
 lBaseq: leqDZ[9,15] - leqdz[9,15] >= 0
 lBaseq: leqDX[9,15] - leqdx[9,15] >= 0
 lBaseq: leqDW[8,15] - leqdw[8,15] >= 0
 l2eq: - lstk[9,0] + leqT[9,0] >= 0
 l2eq: - lstk[9,1] + leqT[9,0] >= 0
 l2eq: - lstk[9,2] + leqT[9,0] >= 0
 l2eq: - lstk[9,3] + leqT[9,0] >= 0
 l2eq: - lstk[9,0] - lstk[9,1] - lstk[9,2] - lstk[9,3] + leqT[9,0] <= 0
 l2eq: lstk[9,0] + lstk[9,1] + lstk[9,2] + lstk[9,3] - 5 leqT[9,0]
   + leqk[9,0] + leqk[9,5] + leqk[9,10] + leqk[9,15] >= 0
 l2eq: lstk[9,0] + lstk[9,1] + lstk[9,2] + lstk[9,3] - 8 leqT[9,0]
   + leqk[9,0] + leqk[9,5] + leqk[9,10] + leqk[9,15] <= 0
 l2eq: - lstk[9,4] + leqT[9,1] >= 0
 l2eq: - lstk[9,5] + leqT[9,1] >= 0
 l2eq: - lstk[9,6] + leqT[9,1] >= 0
 l2eq: - lstk[9,7] + leqT[9,1] >= 0
 l2eq: - lstk[9,4] - lstk[9,5] - lstk[9,6] - lstk[9,7] + leqT[9,1] <= 0
 l2eq: lstk[9,4] + lstk[9,5] + lstk[9,6] + lstk[9,7] - 5 leqT[9,1]
   + leqk[9,3] + leqk[9,4] + leqk[9,9] + leqk[9,14] >= 0
 l2eq: lstk[9,4] + lstk[9,5] + lstk[9,6] + lstk[9,7] - 8 leqT[9,1]
   + leqk[9,3] + leqk[9,4] + leqk[9,9] + leqk[9,14] <= 0
 l2eq: - lstk[9,8] + leqT[9,2] >= 0
 l2eq: - lstk[9,9] + leqT[9,2] >= 0
 l2eq: - lstk[9,10] + leqT[9,2] >= 0
 l2eq: - lstk[9,11] + leqT[9,2] >= 0
 l2eq: - lstk[9,8] - lstk[9,9] - lstk[9,10] - lstk[9,11] + leqT[9,2] <= 0
 l2eq: lstk[9,8] + lstk[9,9] + lstk[9,10] + lstk[9,11] - 5 leqT[9,2]
   + leqk[9,2] + leqk[9,7] + leqk[9,8] + leqk[9,13] >= 0
 l2eq: lstk[9,8] + lstk[9,9] + lstk[9,10] + lstk[9,11] - 8 leqT[9,2]
   + leqk[9,2] + leqk[9,7] + leqk[9,8] + leqk[9,13] <= 0
 l2eq: - lstk[9,12] + leqT[9,3] >= 0
 l2eq: - lstk[9,13] + leqT[9,3] >= 0
 l2eq: - lstk[9,14] + leqT[9,3] >= 0
 l2eq: - lstk[9,15] + leqT[9,3] >= 0
 l2eq: - lstk[9,12] - lstk[9,13] - lstk[9,14] - lstk[9,15] + leqT[9,3]
   <= 0
 l2eq: lstk[9,12] + lstk[9,13] + lstk[9,14] + lstk[9,15] - 5 leqT[9,3]
   + leqk[9,1] + leqk[9,6] + leqk[9,11] + leqk[9,12] >= 0
 l2eq: lstk[9,12] + lstk[9,13] + lstk[9,14] + lstk[9,15] - 8 leqT[9,3]
   + leqk[9,1] + leqk[9,6] + leqk[9,11] + leqk[9,12] <= 0
 l2eq: - lstk[10,0] + leqT[10,0] >= 0
 l2eq: - lstk[10,1] + leqT[10,0] >= 0
 l2eq: - lstk[10,2] + leqT[10,0] >= 0
 l2eq: - lstk[10,3] + leqT[10,0] >= 0
 l2eq: - lstk[10,0] - lstk[10,1] - lstk[10,2] - lstk[10,3] + leqT[10,0]
   <= 0
 l2eq: lstk[10,0] + lstk[10,1] + lstk[10,2] + lstk[10,3] - 5 leqT[10,0]
   + leqk[10,0] + leqk[10,5] + leqk[10,10] + leqk[10,15] >= 0
 l2eq: lstk[10,0] + lstk[10,1] + lstk[10,2] + lstk[10,3] - 8 leqT[10,0]
   + leqk[10,0] + leqk[10,5] + leqk[10,10] + leqk[10,15] <= 0
 l2eq: - lstk[10,4] + leqT[10,1] >= 0
 l2eq: - lstk[10,5] + leqT[10,1] >= 0
 l2eq: - lstk[10,6] + leqT[10,1] >= 0
 l2eq: - lstk[10,7] + leqT[10,1] >= 0
 l2eq: - lstk[10,4] - lstk[10,5] - lstk[10,6] - lstk[10,7] + leqT[10,1]
   <= 0
 l2eq: lstk[10,4] + lstk[10,5] + lstk[10,6] + lstk[10,7] - 5 leqT[10,1]
   + leqk[10,3] + leqk[10,4] + leqk[10,9] + leqk[10,14] >= 0
 l2eq: lstk[10,4] + lstk[10,5] + lstk[10,6] + lstk[10,7] - 8 leqT[10,1]
   + leqk[10,3] + leqk[10,4] + leqk[10,9] + leqk[10,14] <= 0
 l2eq: - lstk[10,8] + leqT[10,2] >= 0
 l2eq: - lstk[10,9] + leqT[10,2] >= 0
 l2eq: - lstk[10,10] + leqT[10,2] >= 0
 l2eq: - lstk[10,11] + leqT[10,2] >= 0
 l2eq: - lstk[10,8] - lstk[10,9] - lstk[10,10] - lstk[10,11] + leqT[10,2]
   <= 0
 l2eq: lstk[10,8] + lstk[10,9] + lstk[10,10] + lstk[10,11] - 5 leqT[10,2]
   + leqk[10,2] + leqk[10,7] + leqk[10,8] + leqk[10,13] >= 0
 l2eq: lstk[10,8] + lstk[10,9] + lstk[10,10] + lstk[10,11] - 8 leqT[10,2]
   + leqk[10,2] + leqk[10,7] + leqk[10,8] + leqk[10,13] <= 0
 l2eq: - lstk[10,12] + leqT[10,3] >= 0
 l2eq: - lstk[10,13] + leqT[10,3] >= 0
 l2eq: - lstk[10,14] + leqT[10,3] >= 0
 l2eq: - lstk[10,15] + leqT[10,3] >= 0
 l2eq: - lstk[10,12] - lstk[10,13] - lstk[10,14] - lstk[10,15] + leqT[10,3]
   <= 0
 l2eq: lstk[10,12] + lstk[10,13] + lstk[10,14] + lstk[10,15] - 5 leqT[10,3]
   + leqk[10,1] + leqk[10,6] + leqk[10,11] + leqk[10,12] >= 0
 l2eq: lstk[10,12] + lstk[10,13] + lstk[10,14] + lstk[10,15] - 8 leqT[10,3]
   + leqk[10,1] + leqk[10,6] + leqk[10,11] + leqk[10,12] <= 0
 lfSC: - lDX[8,0] + leqDY[8,0] = 0
 lfSC: - lDX[8,0] + leqdy[8,0] = 0
 lfSC: - lDX[8,1] + leqDY[8,1] = 0
 lfSC: - lDX[8,1] + leqdy[8,1] = 0
 lfSC: - lDX[8,2] + leqDY[8,2] = 0
 lfSC: - lDX[8,2] + leqdy[8,2] = 0
 lfSC: - lDX[8,3] + leqDY[8,3] = 0
 lfSC: - lDX[8,3] + leqdy[8,3] = 0
 lfSC: - lDX[8,4] + leqDY[8,4] = 0
 lfSC: - lDX[8,4] + leqdy[8,4] = 0
 lfSC: - lDX[8,5] + leqDY[8,5] = 0
 lfSC: - lDX[8,5] + leqdy[8,5] = 0
 lfSC: - lDX[8,6] + leqDY[8,6] = 0
 lfSC: - lDX[8,6] + leqdy[8,6] = 0
 lfSC: - lDX[8,7] + leqDY[8,7] = 0
 lfSC: - lDX[8,7] + leqdy[8,7] = 0
 lfSC: - lDX[8,8] + leqDY[8,8] = 0
 lfSC: - lDX[8,8] + leqdy[8,8] = 0
 lfSC: - lDX[8,9] + leqDY[8,9] = 0
 lfSC: - lDX[8,9] + leqdy[8,9] = 0
 lfSC: - lDX[8,10] + leqDY[8,10] = 0
 lfSC: - lDX[8,10] + leqdy[8,10] = 0
 lfSC: - lDX[8,11] + leqDY[8,11] = 0
 lfSC: - lDX[8,11] + leqdy[8,11] = 0
 lfSC: - lDX[8,12] + leqDY[8,12] = 0
 lfSC: - lDX[8,12] + leqdy[8,12] = 0
 lfSC: - lDX[8,13] + leqDY[8,13] = 0
 lfSC: - lDX[8,13] + leqdy[8,13] = 0
 lfSC: - lDX[8,14] + leqDY[8,14] = 0
 lfSC: - lDX[8,14] + leqdy[8,14] = 0
 lfSC: - lDX[8,15] + leqDY[8,15] = 0
 lfSC: - lDX[8,15] + leqdy[8,15] = 0
 lfSC: - leqDX[9,0] + leqDY[9,0] = 0
 lfSC: - leqDX[9,0] + leqdy[9,0] = 0
 lfSC: - leqDX[9,1] + leqDY[9,1] = 0
 lfSC: - leqDX[9,1] + leqdy[9,1] = 0
 lfSC: - leqDX[9,2] + leqDY[9,2] = 0
 lfSC: - leqDX[9,2] + leqdy[9,2] = 0
 lfSC: - leqDX[9,3] + leqDY[9,3] = 0
 lfSC: - leqDX[9,3] + leqdy[9,3] = 0
 lfSC: - leqDX[9,4] + leqDY[9,4] = 0
 lfSC: - leqDX[9,4] + leqdy[9,4] = 0
 lfSC: - leqDX[9,5] + leqDY[9,5] = 0
 lfSC: - leqDX[9,5] + leqdy[9,5] = 0
 lfSC: - leqDX[9,6] + leqDY[9,6] = 0
 lfSC: - leqDX[9,6] + leqdy[9,6] = 0
 lfSC: - leqDX[9,7] + leqDY[9,7] = 0
 lfSC: - leqDX[9,7] + leqdy[9,7] = 0
 lfSC: - leqDX[9,8] + leqDY[9,8] = 0
 lfSC: - leqDX[9,8] + leqdy[9,8] = 0
 lfSC: - leqDX[9,9] + leqDY[9,9] = 0
 lfSC: - leqDX[9,9] + leqdy[9,9] = 0
 lfSC: - leqDX[9,10] + leqDY[9,10] = 0
 lfSC: - leqDX[9,10] + leqdy[9,10] = 0
 lfSC: - leqDX[9,11] + leqDY[9,11] = 0
 lfSC: - leqDX[9,11] + leqdy[9,11] = 0
 lfSC: - leqDX[9,12] + leqDY[9,12] = 0
 lfSC: - leqDX[9,12] + leqdy[9,12] = 0
 lfSC: - leqDX[9,13] + leqDY[9,13] = 0
 lfSC: - leqDX[9,13] + leqdy[9,13] = 0
 lfSC: - leqDX[9,14] + leqDY[9,14] = 0
 lfSC: - leqDX[9,14] + leqdy[9,14] = 0
 lfSC: - leqDX[9,15] + leqDY[9,15] = 0
 lfSC: - leqDX[9,15] + leqdy[9,15] = 0
 lfATK: leqdy[8,0] - leqdz[8,0] = 0
 lfATK: - leqDZ[8,0] + leqdz[8,0] + leqk[9,0] >= 0
 lfATK: leqDZ[8,0] - leqk[9,0] >= 0
 lfATK: leqdy[8,1] - leqdz[8,1] = 0
 lfATK: - leqDZ[8,1] + leqdz[8,1] + leqk[9,1] >= 0
 lfATK: leqDZ[8,1] - leqk[9,1] >= 0
 lfATK: leqdy[8,2] - leqdz[8,2] = 0
 lfATK: - leqDZ[8,2] + leqdz[8,2] + leqk[9,2] >= 0
 lfATK: leqDZ[8,2] - leqk[9,2] >= 0
 lfATK: leqdy[8,3] - leqdz[8,3] = 0
 lfATK: - leqDZ[8,3] + leqdz[8,3] + leqk[9,3] >= 0
 lfATK: leqDZ[8,3] - leqk[9,3] >= 0
 lfATK: leqdy[8,4] - leqdz[8,4] = 0
 lfATK: - leqDZ[8,4] + leqdz[8,4] + leqk[9,4] >= 0
 lfATK: leqDZ[8,4] - leqk[9,4] >= 0
 lfATK: leqdy[8,5] - leqdz[8,5] = 0
 lfATK: - leqDZ[8,5] + leqdz[8,5] + leqk[9,5] >= 0
 lfATK: leqDZ[8,5] - leqk[9,5] >= 0
 lfATK: leqdy[8,6] - leqdz[8,6] = 0
 lfATK: - leqDZ[8,6] + leqdz[8,6] + leqk[9,6] >= 0
 lfATK: leqDZ[8,6] - leqk[9,6] >= 0
 lfATK: leqdy[8,7] - leqdz[8,7] = 0
 lfATK: - leqDZ[8,7] + leqdz[8,7] + leqk[9,7] >= 0
 lfATK: leqDZ[8,7] - leqk[9,7] >= 0
 lfATK: leqdy[8,8] - leqdz[8,8] = 0
 lfATK: - leqDZ[8,8] + leqdz[8,8] + leqk[9,8] >= 0
 lfATK: leqDZ[8,8] - leqk[9,8] >= 0
 lfATK: leqdy[8,9] - leqdz[8,9] = 0
 lfATK: - leqDZ[8,9] + leqdz[8,9] + leqk[9,9] >= 0
 lfATK: leqDZ[8,9] - leqk[9,9] >= 0
 lfATK: leqdy[8,10] - leqdz[8,10] = 0
 lfATK: - leqDZ[8,10] + leqdz[8,10] + leqk[9,10] >= 0
 lfATK: leqDZ[8,10] - leqk[9,10] >= 0
 lfATK: leqdy[8,11] - leqdz[8,11] = 0
 lfATK: - leqDZ[8,11] + leqdz[8,11] + leqk[9,11] >= 0
 lfATK: leqDZ[8,11] - leqk[9,11] >= 0
 lfATK: leqdy[8,12] - leqdz[8,12] = 0
 lfATK: - leqDZ[8,12] + leqdz[8,12] + leqk[9,12] >= 0
 lfATK: leqDZ[8,12] - leqk[9,12] >= 0
 lfATK: leqdy[8,13] - leqdz[8,13] = 0
 lfATK: - leqDZ[8,13] + leqdz[8,13] + leqk[9,13] >= 0
 lfATK: leqDZ[8,13] - leqk[9,13] >= 0
 lfATK: leqdy[8,14] - leqdz[8,14] = 0
 lfATK: - leqDZ[8,14] + leqdz[8,14] + leqk[9,14] >= 0
 lfATK: leqDZ[8,14] - leqk[9,14] >= 0
 lfATK: leqdy[8,15] - leqdz[8,15] = 0
 lfATK: - leqDZ[8,15] + leqdz[8,15] + leqk[9,15] >= 0
 lfATK: leqDZ[8,15] - leqk[9,15] >= 0
 lfATK: leqdy[9,0] - leqdz[9,0] = 0
 lfATK: - leqDZ[9,0] + leqdz[9,0] + leqk[10,0] >= 0
 lfATK: leqDZ[9,0] - leqk[10,0] >= 0
 lfATK: leqdy[9,1] - leqdz[9,1] = 0
 lfATK: - leqDZ[9,1] + leqdz[9,1] + leqk[10,1] >= 0
 lfATK: leqDZ[9,1] - leqk[10,1] >= 0
 lfATK: leqdy[9,2] - leqdz[9,2] = 0
 lfATK: - leqDZ[9,2] + leqdz[9,2] + leqk[10,2] >= 0
 lfATK: leqDZ[9,2] - leqk[10,2] >= 0
 lfATK: leqdy[9,3] - leqdz[9,3] = 0
 lfATK: - leqDZ[9,3] + leqdz[9,3] + leqk[10,3] >= 0
 lfATK: leqDZ[9,3] - leqk[10,3] >= 0
 lfATK: leqdy[9,4] - leqdz[9,4] = 0
 lfATK: - leqDZ[9,4] + leqdz[9,4] + leqk[10,4] >= 0
 lfATK: leqDZ[9,4] - leqk[10,4] >= 0
 lfATK: leqdy[9,5] - leqdz[9,5] = 0
 lfATK: - leqDZ[9,5] + leqdz[9,5] + leqk[10,5] >= 0
 lfATK: leqDZ[9,5] - leqk[10,5] >= 0
 lfATK: leqdy[9,6] - leqdz[9,6] = 0
 lfATK: - leqDZ[9,6] + leqdz[9,6] + leqk[10,6] >= 0
 lfATK: leqDZ[9,6] - leqk[10,6] >= 0
 lfATK: leqdy[9,7] - leqdz[9,7] = 0
 lfATK: - leqDZ[9,7] + leqdz[9,7] + leqk[10,7] >= 0
 lfATK: leqDZ[9,7] - leqk[10,7] >= 0
 lfATK: leqdy[9,8] - leqdz[9,8] = 0
 lfATK: - leqDZ[9,8] + leqdz[9,8] + leqk[10,8] >= 0
 lfATK: leqDZ[9,8] - leqk[10,8] >= 0
 lfATK: leqdy[9,9] - leqdz[9,9] = 0
 lfATK: - leqDZ[9,9] + leqdz[9,9] + leqk[10,9] >= 0
 lfATK: leqDZ[9,9] - leqk[10,9] >= 0
 lfATK: leqdy[9,10] - leqdz[9,10] = 0
 lfATK: - leqDZ[9,10] + leqdz[9,10] + leqk[10,10] >= 0
 lfATK: leqDZ[9,10] - leqk[10,10] >= 0
 lfATK: leqdy[9,11] - leqdz[9,11] = 0
 lfATK: - leqDZ[9,11] + leqdz[9,11] + leqk[10,11] >= 0
 lfATK: leqDZ[9,11] - leqk[10,11] >= 0
 lfATK: leqdy[9,12] - leqdz[9,12] = 0
 lfATK: - leqDZ[9,12] + leqdz[9,12] + leqk[10,12] >= 0
 lfATK: leqDZ[9,12] - leqk[10,12] >= 0
 lfATK: leqdy[9,13] - leqdz[9,13] = 0
 lfATK: - leqDZ[9,13] + leqdz[9,13] + leqk[10,13] >= 0
 lfATK: leqDZ[9,13] - leqk[10,13] >= 0
 lfATK: leqdy[9,14] - leqdz[9,14] = 0
 lfATK: - leqDZ[9,14] + leqdz[9,14] + leqk[10,14] >= 0
 lfATK: leqDZ[9,14] - leqk[10,14] >= 0
 lfATK: leqdy[9,15] - leqdz[9,15] = 0
 lfATK: - leqDZ[9,15] + leqdz[9,15] + leqk[10,15] >= 0
 lfATK: leqDZ[9,15] - leqk[10,15] >= 0
 lfSR: - leqDW[8,0] + leqDZ[8,0] = 0
 lfSR: - leqdw[8,0] + leqdz[8,0] = 0
 lfSR: - leqDW[8,1] + leqDZ[8,5] = 0
 lfSR: - leqdw[8,1] + leqdz[8,5] = 0
 lfSR: - leqDW[8,2] + leqDZ[8,10] = 0
 lfSR: - leqdw[8,2] + leqdz[8,10] = 0
 lfSR: - leqDW[8,3] + leqDZ[8,15] = 0
 lfSR: - leqdw[8,3] + leqdz[8,15] = 0
 lfSR: - leqDW[8,4] + leqDZ[8,4] = 0
 lfSR: - leqdw[8,4] + leqdz[8,4] = 0
 lfSR: - leqDW[8,5] + leqDZ[8,9] = 0
 lfSR: - leqdw[8,5] + leqdz[8,9] = 0
 lfSR: - leqDW[8,6] + leqDZ[8,14] = 0
 lfSR: - leqdw[8,6] + leqdz[8,14] = 0
 lfSR: - leqDW[8,7] + leqDZ[8,3] = 0
 lfSR: - leqdw[8,7] + leqdz[8,3] = 0
 lfSR: - leqDW[8,8] + leqDZ[8,8] = 0
 lfSR: - leqdw[8,8] + leqdz[8,8] = 0
 lfSR: - leqDW[8,9] + leqDZ[8,13] = 0
 lfSR: - leqdw[8,9] + leqdz[8,13] = 0
 lfSR: - leqDW[8,10] + leqDZ[8,2] = 0
 lfSR: - leqdw[8,10] + leqdz[8,2] = 0
 lfSR: - leqDW[8,11] + leqDZ[8,7] = 0
 lfSR: - leqdw[8,11] + leqdz[8,7] = 0
 lfSR: - leqDW[8,12] + leqDZ[8,12] = 0
 lfSR: - leqdw[8,12] + leqdz[8,12] = 0
 lfSR: - leqDW[8,13] + leqDZ[8,1] = 0
 lfSR: - leqdw[8,13] + leqdz[8,1] = 0
 lfSR: - leqDW[8,14] + leqDZ[8,6] = 0
 lfSR: - leqdw[8,14] + leqdz[8,6] = 0
 lfSR: - leqDW[8,15] + leqDZ[8,11] = 0
 lfSR: - leqdw[8,15] + leqdz[8,11] = 0
 lfMC: lTC[8,0] - leqdw[8,0] >= 0
 lfMC: lTC[8,0] - leqdw[8,1] >= 0
 lfMC: lTC[8,0] - leqdw[8,2] >= 0
 lfMC: lTC[8,0] - leqdw[8,3] >= 0
 lfMC: lTC[8,0] - leqdw[8,0] - leqdw[8,1] - leqdw[8,2] - leqdw[8,3] <= 0
 lfMC: lTC[8,0] - leqdx[9,0] <= 0
 lfMC: lTC[8,0] - leqdx[9,1] <= 0
 lfMC: lTC[8,0] - leqdx[9,2] <= 0
 lfMC: lTC[8,0] - leqdx[9,3] <= 0
 lfMC: 4 lTC[8,0] - leqdx[9,0] - leqdx[9,1] - leqdx[9,2] - leqdx[9,3] >= 0
 lfMC: lAC[8,0] - leqDW[8,0] >= 0
 lfMC: lAC[8,0] - leqDW[8,1] >= 0
 lfMC: lAC[8,0] - leqDW[8,2] >= 0
 lfMC: lAC[8,0] - leqDW[8,3] >= 0
 lfMC: lAC[8,0] - leqDW[8,0] - leqDW[8,1] - leqDW[8,2] - leqDW[8,3] <= 0
 lfMC: - 5 lAC[8,0] + leqDW[8,0] + leqDW[8,1] + leqDW[8,2] + leqDW[8,3]
   + leqDX[9,0] + leqDX[9,1] + leqDX[9,2] + leqDX[9,3] >= 0
 lfMC: - 8 lAC[8,0] + leqDW[8,0] + leqDW[8,1] + leqDW[8,2] + leqDW[8,3]
   + leqDX[9,0] + leqDX[9,1] + leqDX[9,2] + leqDX[9,3] <= 0
 lfMC: lTC[8,1] - leqdw[8,4] >= 0
 lfMC: lTC[8,1] - leqdw[8,5] >= 0
 lfMC: lTC[8,1] - leqdw[8,6] >= 0
 lfMC: lTC[8,1] - leqdw[8,7] >= 0
 lfMC: lTC[8,1] - leqdw[8,4] - leqdw[8,5] - leqdw[8,6] - leqdw[8,7] <= 0
 lfMC: lTC[8,1] - leqdx[9,4] <= 0
 lfMC: lTC[8,1] - leqdx[9,5] <= 0
 lfMC: lTC[8,1] - leqdx[9,6] <= 0
 lfMC: lTC[8,1] - leqdx[9,7] <= 0
 lfMC: 4 lTC[8,1] - leqdx[9,4] - leqdx[9,5] - leqdx[9,6] - leqdx[9,7] >= 0
 lfMC: lAC[8,1] - leqDW[8,4] >= 0
 lfMC: lAC[8,1] - leqDW[8,5] >= 0
 lfMC: lAC[8,1] - leqDW[8,6] >= 0
 lfMC: lAC[8,1] - leqDW[8,7] >= 0
 lfMC: lAC[8,1] - leqDW[8,4] - leqDW[8,5] - leqDW[8,6] - leqDW[8,7] <= 0
 lfMC: - 5 lAC[8,1] + leqDW[8,4] + leqDW[8,5] + leqDW[8,6] + leqDW[8,7]
   + leqDX[9,4] + leqDX[9,5] + leqDX[9,6] + leqDX[9,7] >= 0
 lfMC: - 8 lAC[8,1] + leqDW[8,4] + leqDW[8,5] + leqDW[8,6] + leqDW[8,7]
   + leqDX[9,4] + leqDX[9,5] + leqDX[9,6] + leqDX[9,7] <= 0
 lfMC: lTC[8,2] - leqdw[8,8] >= 0
 lfMC: lTC[8,2] - leqdw[8,9] >= 0
 lfMC: lTC[8,2] - leqdw[8,10] >= 0
 lfMC: lTC[8,2] - leqdw[8,11] >= 0
 lfMC: lTC[8,2] - leqdw[8,8] - leqdw[8,9] - leqdw[8,10] - leqdw[8,11] <= 0
 lfMC: lTC[8,2] - leqdx[9,8] <= 0
 lfMC: lTC[8,2] - leqdx[9,9] <= 0
 lfMC: lTC[8,2] - leqdx[9,10] <= 0
 lfMC: lTC[8,2] - leqdx[9,11] <= 0
 lfMC: 4 lTC[8,2] - leqdx[9,8] - leqdx[9,9] - leqdx[9,10] - leqdx[9,11]
   >= 0
 lfMC: lAC[8,2] - leqDW[8,8] >= 0
 lfMC: lAC[8,2] - leqDW[8,9] >= 0
 lfMC: lAC[8,2] - leqDW[8,10] >= 0
 lfMC: lAC[8,2] - leqDW[8,11] >= 0
 lfMC: lAC[8,2] - leqDW[8,8] - leqDW[8,9] - leqDW[8,10] - leqDW[8,11] <= 0
 lfMC: - 5 lAC[8,2] + leqDW[8,8] + leqDW[8,9] + leqDW[8,10] + leqDW[8,11]
   + leqDX[9,8] + leqDX[9,9] + leqDX[9,10] + leqDX[9,11] >= 0
 lfMC: - 8 lAC[8,2] + leqDW[8,8] + leqDW[8,9] + leqDW[8,10] + leqDW[8,11]
   + leqDX[9,8] + leqDX[9,9] + leqDX[9,10] + leqDX[9,11] <= 0
 lfMC: lTC[8,3] - leqdw[8,12] >= 0
 lfMC: lTC[8,3] - leqdw[8,13] >= 0
 lfMC: lTC[8,3] - leqdw[8,14] >= 0
 lfMC: lTC[8,3] - leqdw[8,15] >= 0
 lfMC: lTC[8,3] - leqdw[8,12] - leqdw[8,13] - leqdw[8,14] - leqdw[8,15]
   <= 0
 lfMC: lTC[8,3] - leqdx[9,12] <= 0
 lfMC: lTC[8,3] - leqdx[9,13] <= 0
 lfMC: lTC[8,3] - leqdx[9,14] <= 0
 lfMC: lTC[8,3] - leqdx[9,15] <= 0
 lfMC: 4 lTC[8,3] - leqdx[9,12] - leqdx[9,13] - leqdx[9,14] - leqdx[9,15]
   >= 0
 lfMC: lAC[8,3] - leqDW[8,12] >= 0
 lfMC: lAC[8,3] - leqDW[8,13] >= 0
 lfMC: lAC[8,3] - leqDW[8,14] >= 0
 lfMC: lAC[8,3] - leqDW[8,15] >= 0
 lfMC: lAC[8,3] - leqDW[8,12] - leqDW[8,13] - leqDW[8,14] - leqDW[8,15]
   <= 0
 lfMC: - 5 lAC[8,3] + leqDW[8,12] + leqDW[8,13] + leqDW[8,14] + leqDW[8,15]
   + leqDX[9,12] + leqDX[9,13] + leqDX[9,14] + leqDX[9,15] >= 0
 lfMC: - 8 lAC[8,3] + leqDW[8,12] + leqDW[8,13] + leqDW[8,14] + leqDW[8,15]
   + leqDX[9,12] + leqDX[9,13] + leqDX[9,14] + leqDX[9,15] <= 0
 R5270: lLANE[0] + lLANE[1] + lLANE[2] + lLANE[3] + lLANE[4] + lLANE[5]
   + lLANE[6] + lLANE[7] + lLANE[8] + lLANE[9] + lLANE[10] + lLANE[11]
   + lLANE[12] + lLANE[13] + lLANE[14] + lLANE[15] >= 1
 lMLANE: - lstk[0,0] + lLANE[0] >= 0
 lMLANE: - lstk[1,7] + lLANE[0] >= 0
 lMLANE: - lstk[2,14] + lLANE[0] >= 0
 lMLANE: - lstk[3,9] + lLANE[0] >= 0
 lMLANE: - lstk[4,8] + lLANE[0] >= 0
 lMLANE: - lstk[5,15] + lLANE[0] >= 0
 lMLANE: - lstk[6,6] + lLANE[0] >= 0
 lMLANE: - lstk[7,1] + lLANE[0] >= 0
 lMLANE: - lstk[8,0] + lLANE[0] >= 0
 lMLANE: - lstk[9,7] + lLANE[0] >= 0
 lMLANE: - lstk[10,14] + lLANE[0] >= 0
 lMLANE: - lstk[0,1] + lLANE[1] >= 0
 lMLANE: - lstk[1,0] + lLANE[1] >= 0
 lMLANE: - lstk[2,7] + lLANE[1] >= 0
 lMLANE: - lstk[3,14] + lLANE[1] >= 0
 lMLANE: - lstk[4,9] + lLANE[1] >= 0
 lMLANE: - lstk[5,8] + lLANE[1] >= 0
 lMLANE: - lstk[6,15] + lLANE[1] >= 0
 lMLANE: - lstk[7,6] + lLANE[1] >= 0
 lMLANE: - lstk[8,1] + lLANE[1] >= 0
 lMLANE: - lstk[9,0] + lLANE[1] >= 0
 lMLANE: - lstk[10,7] + lLANE[1] >= 0
 lMLANE: - lstk[0,2] + lLANE[2] >= 0
 lMLANE: - lstk[1,13] + lLANE[2] >= 0
 lMLANE: - lstk[2,12] + lLANE[2] >= 0
 lMLANE: - lstk[3,3] + lLANE[2] >= 0
 lMLANE: - lstk[4,10] + lLANE[2] >= 0
 lMLANE: - lstk[5,5] + lLANE[2] >= 0
 lMLANE: - lstk[6,4] + lLANE[2] >= 0
 lMLANE: - lstk[7,11] + lLANE[2] >= 0
 lMLANE: - lstk[8,2] + lLANE[2] >= 0
 lMLANE: - lstk[9,13] + lLANE[2] >= 0
 lMLANE: - lstk[10,12] + lLANE[2] >= 0
 lMLANE: - lstk[0,3] + lLANE[3] >= 0
 lMLANE: - lstk[1,10] + lLANE[3] >= 0
 lMLANE: - lstk[2,5] + lLANE[3] >= 0
 lMLANE: - lstk[3,4] + lLANE[3] >= 0
 lMLANE: - lstk[4,11] + lLANE[3] >= 0
 lMLANE: - lstk[5,2] + lLANE[3] >= 0
 lMLANE: - lstk[6,13] + lLANE[3] >= 0
 lMLANE: - lstk[7,12] + lLANE[3] >= 0
 lMLANE: - lstk[8,3] + lLANE[3] >= 0
 lMLANE: - lstk[9,10] + lLANE[3] >= 0
 lMLANE: - lstk[10,5] + lLANE[3] >= 0
 lMLANE: - lstk[0,4] + lLANE[4] >= 0
 lMLANE: - lstk[1,11] + lLANE[4] >= 0
 lMLANE: - lstk[2,2] + lLANE[4] >= 0
 lMLANE: - lstk[3,13] + lLANE[4] >= 0
 lMLANE: - lstk[4,12] + lLANE[4] >= 0
 lMLANE: - lstk[5,3] + lLANE[4] >= 0
 lMLANE: - lstk[6,10] + lLANE[4] >= 0
 lMLANE: - lstk[7,5] + lLANE[4] >= 0
 lMLANE: - lstk[8,4] + lLANE[4] >= 0
 lMLANE: - lstk[9,11] + lLANE[4] >= 0
 lMLANE: - lstk[10,2] + lLANE[4] >= 0
 lMLANE: - lstk[0,5] + lLANE[5] >= 0
 lMLANE: - lstk[1,4] + lLANE[5] >= 0
 lMLANE: - lstk[2,11] + lLANE[5] >= 0
 lMLANE: - lstk[3,2] + lLANE[5] >= 0
 lMLANE: - lstk[4,13] + lLANE[5] >= 0
 lMLANE: - lstk[5,12] + lLANE[5] >= 0
 lMLANE: - lstk[6,3] + lLANE[5] >= 0
 lMLANE: - lstk[7,10] + lLANE[5] >= 0
 lMLANE: - lstk[8,5] + lLANE[5] >= 0
 lMLANE: - lstk[9,4] + lLANE[5] >= 0
 lMLANE: - lstk[10,11] + lLANE[5] >= 0
 lMLANE: - lstk[0,6] + lLANE[6] >= 0
 lMLANE: - lstk[1,1] + lLANE[6] >= 0
 lMLANE: - lstk[2,0] + lLANE[6] >= 0
 lMLANE: - lstk[3,7] + lLANE[6] >= 0
 lMLANE: - lstk[4,14] + lLANE[6] >= 0
 lMLANE: - lstk[5,9] + lLANE[6] >= 0
 lMLANE: - lstk[6,8] + lLANE[6] >= 0
 lMLANE: - lstk[7,15] + lLANE[6] >= 0
 lMLANE: - lstk[8,6] + lLANE[6] >= 0
 lMLANE: - lstk[9,1] + lLANE[6] >= 0
 lMLANE: - lstk[10,0] + lLANE[6] >= 0
 lMLANE: - lstk[0,7] + lLANE[7] >= 0
 lMLANE: - lstk[1,14] + lLANE[7] >= 0
 lMLANE: - lstk[2,9] + lLANE[7] >= 0
 lMLANE: - lstk[3,8] + lLANE[7] >= 0
 lMLANE: - lstk[4,15] + lLANE[7] >= 0
 lMLANE: - lstk[5,6] + lLANE[7] >= 0
 lMLANE: - lstk[6,1] + lLANE[7] >= 0
 lMLANE: - lstk[7,0] + lLANE[7] >= 0
 lMLANE: - lstk[8,7] + lLANE[7] >= 0
 lMLANE: - lstk[9,14] + lLANE[7] >= 0
 lMLANE: - lstk[10,9] + lLANE[7] >= 0
 lMLANE: - lstk[0,8] + lLANE[8] >= 0
 lMLANE: - lstk[1,15] + lLANE[8] >= 0
 lMLANE: - lstk[2,6] + lLANE[8] >= 0
 lMLANE: - lstk[3,1] + lLANE[8] >= 0
 lMLANE: - lstk[4,0] + lLANE[8] >= 0
 lMLANE: - lstk[5,7] + lLANE[8] >= 0
 lMLANE: - lstk[6,14] + lLANE[8] >= 0
 lMLANE: - lstk[7,9] + lLANE[8] >= 0
 lMLANE: - lstk[8,8] + lLANE[8] >= 0
 lMLANE: - lstk[9,15] + lLANE[8] >= 0
 lMLANE: - lstk[10,6] + lLANE[8] >= 0
 lMLANE: - lstk[0,9] + lLANE[9] >= 0
 lMLANE: - lstk[1,8] + lLANE[9] >= 0
 lMLANE: - lstk[2,15] + lLANE[9] >= 0
 lMLANE: - lstk[3,6] + lLANE[9] >= 0
 lMLANE: - lstk[4,1] + lLANE[9] >= 0
 lMLANE: - lstk[5,0] + lLANE[9] >= 0
 lMLANE: - lstk[6,7] + lLANE[9] >= 0
 lMLANE: - lstk[7,14] + lLANE[9] >= 0
 lMLANE: - lstk[8,9] + lLANE[9] >= 0
 lMLANE: - lstk[9,8] + lLANE[9] >= 0
 lMLANE: - lstk[10,15] + lLANE[9] >= 0
 lMLANE: - lstk[0,10] + lLANE[10] >= 0
 lMLANE: - lstk[1,5] + lLANE[10] >= 0
 lMLANE: - lstk[2,4] + lLANE[10] >= 0
 lMLANE: - lstk[3,11] + lLANE[10] >= 0
 lMLANE: - lstk[4,2] + lLANE[10] >= 0
 lMLANE: - lstk[5,13] + lLANE[10] >= 0
 lMLANE: - lstk[6,12] + lLANE[10] >= 0
 lMLANE: - lstk[7,3] + lLANE[10] >= 0
 lMLANE: - lstk[8,10] + lLANE[10] >= 0
 lMLANE: - lstk[9,5] + lLANE[10] >= 0
 lMLANE: - lstk[10,4] + lLANE[10] >= 0
 lMLANE: - lstk[0,11] + lLANE[11] >= 0
 lMLANE: - lstk[1,2] + lLANE[11] >= 0
 lMLANE: - lstk[2,13] + lLANE[11] >= 0
 lMLANE: - lstk[3,12] + lLANE[11] >= 0
 lMLANE: - lstk[4,3] + lLANE[11] >= 0
 lMLANE: - lstk[5,10] + lLANE[11] >= 0
 lMLANE: - lstk[6,5] + lLANE[11] >= 0
 lMLANE: - lstk[7,4] + lLANE[11] >= 0
 lMLANE: - lstk[8,11] + lLANE[11] >= 0
 lMLANE: - lstk[9,2] + lLANE[11] >= 0
 lMLANE: - lstk[10,13] + lLANE[11] >= 0
 lMLANE: - lstk[0,12] + lLANE[12] >= 0
 lMLANE: - lstk[1,3] + lLANE[12] >= 0
 lMLANE: - lstk[2,10] + lLANE[12] >= 0
 lMLANE: - lstk[3,5] + lLANE[12] >= 0
 lMLANE: - lstk[4,4] + lLANE[12] >= 0
 lMLANE: - lstk[5,11] + lLANE[12] >= 0
 lMLANE: - lstk[6,2] + lLANE[12] >= 0
 lMLANE: - lstk[7,13] + lLANE[12] >= 0
 lMLANE: - lstk[8,12] + lLANE[12] >= 0
 lMLANE: - lstk[9,3] + lLANE[12] >= 0
 lMLANE: - lstk[10,10] + lLANE[12] >= 0
 lMLANE: - lstk[0,13] + lLANE[13] >= 0
 lMLANE: - lstk[1,12] + lLANE[13] >= 0
 lMLANE: - lstk[2,3] + lLANE[13] >= 0
 lMLANE: - lstk[3,10] + lLANE[13] >= 0
 lMLANE: - lstk[4,5] + lLANE[13] >= 0
 lMLANE: - lstk[5,4] + lLANE[13] >= 0
 lMLANE: - lstk[6,11] + lLANE[13] >= 0
 lMLANE: - lstk[7,2] + lLANE[13] >= 0
 lMLANE: - lstk[8,13] + lLANE[13] >= 0
 lMLANE: - lstk[9,12] + lLANE[13] >= 0
 lMLANE: - lstk[10,3] + lLANE[13] >= 0
 lMLANE: - lstk[0,14] + lLANE[14] >= 0
 lMLANE: - lstk[1,9] + lLANE[14] >= 0
 lMLANE: - lstk[2,8] + lLANE[14] >= 0
 lMLANE: - lstk[3,15] + lLANE[14] >= 0
 lMLANE: - lstk[4,6] + lLANE[14] >= 0
 lMLANE: - lstk[5,1] + lLANE[14] >= 0
 lMLANE: - lstk[6,0] + lLANE[14] >= 0
 lMLANE: - lstk[7,7] + lLANE[14] >= 0
 lMLANE: - lstk[8,14] + lLANE[14] >= 0
 lMLANE: - lstk[9,9] + lLANE[14] >= 0
 lMLANE: - lstk[10,8] + lLANE[14] >= 0
 lMLANE: - lstk[0,15] + lLANE[15] >= 0
 lMLANE: - lstk[1,6] + lLANE[15] >= 0
 lMLANE: - lstk[2,1] + lLANE[15] >= 0
 lMLANE: - lstk[3,0] + lLANE[15] >= 0
 lMLANE: - lstk[4,7] + lLANE[15] >= 0
 lMLANE: - lstk[5,14] + lLANE[15] >= 0
 lMLANE: - lstk[6,9] + lLANE[15] >= 0
 lMLANE: - lstk[7,8] + lLANE[15] >= 0
 lMLANE: - lstk[8,15] + lLANE[15] >= 0
 lMLANE: - lstk[9,6] + lLANE[15] >= 0
 lMLANE: - lstk[10,1] + lLANE[15] >= 0
 lType1C: lstk[0,0] + lstk[1,7] + lstk[2,14] + lstk[3,9] + lstk[4,8]
   + lstk[5,15] + lstk[6,6] + lstk[7,1] + lstk[8,0] + lstk[9,7]
   + lstk[10,14] - 11 lLANE[0] + lType1[0] = 0
 lType1C: lType1[0] <= 1
 lType1C: lstk[0,1] + lstk[1,0] + lstk[2,7] + lstk[3,14] + lstk[4,9]
   + lstk[5,8] + lstk[6,15] + lstk[7,6] + lstk[8,1] + lstk[9,0]
   + lstk[10,7] - 11 lLANE[1] + lType1[1] = 0
 lType1C: lType1[1] <= 1
 lType1C: lstk[0,2] + lstk[1,13] + lstk[2,12] + lstk[3,3] + lstk[4,10]
   + lstk[5,5] + lstk[6,4] + lstk[7,11] + lstk[8,2] + lstk[9,13]
   + lstk[10,12] - 11 lLANE[2] + lType1[2] = 0
 lType1C: lType1[2] <= 1
 lType1C: lstk[0,3] + lstk[1,10] + lstk[2,5] + lstk[3,4] + lstk[4,11]
   + lstk[5,2] + lstk[6,13] + lstk[7,12] + lstk[8,3] + lstk[9,10]
   + lstk[10,5] - 11 lLANE[3] + lType1[3] = 0
 lType1C: lType1[3] <= 1
 lType1C: lstk[0,4] + lstk[1,11] + lstk[2,2] + lstk[3,13] + lstk[4,12]
   + lstk[5,3] + lstk[6,10] + lstk[7,5] + lstk[8,4] + lstk[9,11]
   + lstk[10,2] - 11 lLANE[4] + lType1[4] = 0
 lType1C: lType1[4] <= 1
 lType1C: lstk[0,5] + lstk[1,4] + lstk[2,11] + lstk[3,2] + lstk[4,13]
   + lstk[5,12] + lstk[6,3] + lstk[7,10] + lstk[8,5] + lstk[9,4]
   + lstk[10,11] - 11 lLANE[5] + lType1[5] = 0
 lType1C: lType1[5] <= 1
 lType1C: lstk[0,6] + lstk[1,1] + lstk[2,0] + lstk[3,7] + lstk[4,14]
   + lstk[5,9] + lstk[6,8] + lstk[7,15] + lstk[8,6] + lstk[9,1]
   + lstk[10,0] - 11 lLANE[6] + lType1[6] = 0
 lType1C: lType1[6] <= 1
 lType1C: lstk[0,7] + lstk[1,14] + lstk[2,9] + lstk[3,8] + lstk[4,15]
   + lstk[5,6] + lstk[6,1] + lstk[7,0] + lstk[8,7] + lstk[9,14]
   + lstk[10,9] - 11 lLANE[7] + lType1[7] = 0
 lType1C: lType1[7] <= 1
 lType1C: lstk[0,8] + lstk[1,15] + lstk[2,6] + lstk[3,1] + lstk[4,0]
   + lstk[5,7] + lstk[6,14] + lstk[7,9] + lstk[8,8] + lstk[9,15]
   + lstk[10,6] - 11 lLANE[8] + lType1[8] = 0
 lType1C: lType1[8] <= 1
 lType1C: lstk[0,9] + lstk[1,8] + lstk[2,15] + lstk[3,6] + lstk[4,1]
   + lstk[5,0] + lstk[6,7] + lstk[7,14] + lstk[8,9] + lstk[9,8]
   + lstk[10,15] - 11 lLANE[9] + lType1[9] = 0
 lType1C: lType1[9] <= 1
 lType1C: lstk[0,10] + lstk[1,5] + lstk[2,4] + lstk[3,11] + lstk[4,2]
   + lstk[5,13] + lstk[6,12] + lstk[7,3] + lstk[8,10] + lstk[9,5]
   + lstk[10,4] - 11 lLANE[10] + lType1[10] = 0
 lType1C: lType1[10] <= 1
 lType1C: lstk[0,11] + lstk[1,2] + lstk[2,13] + lstk[3,12] + lstk[4,3]
   + lstk[5,10] + lstk[6,5] + lstk[7,4] + lstk[8,11] + lstk[9,2]
   + lstk[10,13] - 11 lLANE[11] + lType1[11] = 0
 lType1C: lType1[11] <= 1
 lType1C: lstk[0,12] + lstk[1,3] + lstk[2,10] + lstk[3,5] + lstk[4,4]
   + lstk[5,11] + lstk[6,2] + lstk[7,13] + lstk[8,12] + lstk[9,3]
   + lstk[10,10] - 11 lLANE[12] + lType1[12] = 0
 lType1C: lType1[12] <= 1
 lType1C: lstk[0,13] + lstk[1,12] + lstk[2,3] + lstk[3,10] + lstk[4,5]
   + lstk[5,4] + lstk[6,11] + lstk[7,2] + lstk[8,13] + lstk[9,12]
   + lstk[10,3] - 11 lLANE[13] + lType1[13] = 0
 lType1C: lType1[13] <= 1
 lType1C: lstk[0,14] + lstk[1,9] + lstk[2,8] + lstk[3,15] + lstk[4,6]
   + lstk[5,1] + lstk[6,0] + lstk[7,7] + lstk[8,14] + lstk[9,9]
   + lstk[10,8] - 11 lLANE[14] + lType1[14] = 0
 lType1C: lType1[14] <= 1
 lType1C: lstk[0,15] + lstk[1,6] + lstk[2,1] + lstk[3,0] + lstk[4,7]
   + lstk[5,14] + lstk[6,9] + lstk[7,8] + lstk[8,15] + lstk[9,6]
   + lstk[10,1] - 11 lLANE[15] + lType1[15] = 0
 lType1C: lType1[15] <= 1
 lType2A: - lDZ[1,0] - lDZ[1,1] - lDZ[1,2] - lDZ[1,3] + lT2CanA[1,0] = 0
 lType1B: lDW[1,0] + lDW[1,1] + lDW[1,2] + lDW[1,3] - 4 lAC[1,0]
   + lT2CanB[1,0] = 0
 lType2C: - lcan[2,0] - lcan[2,1] - lcan[2,2] - lcan[2,3] + lT2CanC[2,0]
   = 0
 lType2A: - lDZ[1,4] - lDZ[1,5] - lDZ[1,6] - lDZ[1,7] + lT2CanA[1,1] = 0
 lType1B: lDW[1,4] + lDW[1,5] + lDW[1,6] + lDW[1,7] - 4 lAC[1,1]
   + lT2CanB[1,1] = 0
 lType2C: - lcan[2,4] - lcan[2,5] - lcan[2,6] - lcan[2,7] + lT2CanC[2,1]
   = 0
 lType2A: - lDZ[1,8] - lDZ[1,9] - lDZ[1,10] - lDZ[1,11] + lT2CanA[1,2] = 0
 lType1B: lDW[1,8] + lDW[1,9] + lDW[1,10] + lDW[1,11] - 4 lAC[1,2]
   + lT2CanB[1,2] = 0
 lType2C: - lcan[2,8] - lcan[2,9] - lcan[2,10] - lcan[2,11] + lT2CanC[2,2]
   = 0
 lType2A: - lDZ[1,12] - lDZ[1,13] - lDZ[1,14] - lDZ[1,15] + lT2CanA[1,3]
   = 0
 lType1B: lDW[1,12] + lDW[1,13] + lDW[1,14] + lDW[1,15] - 4 lAC[1,3]
   + lT2CanB[1,3] = 0
 lType2C: - lcan[2,12] - lcan[2,13] - lcan[2,14] - lcan[2,15]
   + lT2CanC[2,3] = 0
 lType2A: - lDZ[2,0] - lDZ[2,1] - lDZ[2,2] - lDZ[2,3] + lT2CanA[2,0] = 0
 lType1B: lDW[2,0] + lDW[2,1] + lDW[2,2] + lDW[2,3] - 4 lAC[2,0]
   + lT2CanB[2,0] = 0
 lType2C: - lcan[3,0] - lcan[3,1] - lcan[3,2] - lcan[3,3] + lT2CanC[3,0]
   = 0
 lType2A: - lDZ[2,4] - lDZ[2,5] - lDZ[2,6] - lDZ[2,7] + lT2CanA[2,1] = 0
 lType1B: lDW[2,4] + lDW[2,5] + lDW[2,6] + lDW[2,7] - 4 lAC[2,1]
   + lT2CanB[2,1] = 0
 lType2C: - lcan[3,4] - lcan[3,5] - lcan[3,6] - lcan[3,7] + lT2CanC[3,1]
   = 0
 lType2A: - lDZ[2,8] - lDZ[2,9] - lDZ[2,10] - lDZ[2,11] + lT2CanA[2,2] = 0
 lType1B: lDW[2,8] + lDW[2,9] + lDW[2,10] + lDW[2,11] - 4 lAC[2,2]
   + lT2CanB[2,2] = 0
 lType2C: - lcan[3,8] - lcan[3,9] - lcan[3,10] - lcan[3,11] + lT2CanC[3,2]
   = 0
 lType2A: - lDZ[2,12] - lDZ[2,13] - lDZ[2,14] - lDZ[2,15] + lT2CanA[2,3]
   = 0
 lType1B: lDW[2,12] + lDW[2,13] + lDW[2,14] + lDW[2,15] - 4 lAC[2,3]
   + lT2CanB[2,3] = 0
 lType2C: - lcan[3,12] - lcan[3,13] - lcan[3,14] - lcan[3,15]
   + lT2CanC[3,3] = 0
 lType2A: - lDZ[3,0] - lDZ[3,1] - lDZ[3,2] - lDZ[3,3] + lT2CanA[3,0] = 0
 lType1B: lDW[3,0] + lDW[3,1] + lDW[3,2] + lDW[3,3] - 4 lAC[3,0]
   + lT2CanB[3,0] = 0
 lType2C: - lcan[4,0] - lcan[4,1] - lcan[4,2] - lcan[4,3] + lT2CanC[4,0]
   = 0
 lType2A: - lDZ[3,4] - lDZ[3,5] - lDZ[3,6] - lDZ[3,7] + lT2CanA[3,1] = 0
 lType1B: lDW[3,4] + lDW[3,5] + lDW[3,6] + lDW[3,7] - 4 lAC[3,1]
   + lT2CanB[3,1] = 0
 lType2C: - lcan[4,4] - lcan[4,5] - lcan[4,6] - lcan[4,7] + lT2CanC[4,1]
   = 0
 lType2A: - lDZ[3,8] - lDZ[3,9] - lDZ[3,10] - lDZ[3,11] + lT2CanA[3,2] = 0
 lType1B: lDW[3,8] + lDW[3,9] + lDW[3,10] + lDW[3,11] - 4 lAC[3,2]
   + lT2CanB[3,2] = 0
 lType2C: - lcan[4,8] - lcan[4,9] - lcan[4,10] - lcan[4,11] + lT2CanC[4,2]
   = 0
 lType2A: - lDZ[3,12] - lDZ[3,13] - lDZ[3,14] - lDZ[3,15] + lT2CanA[3,3]
   = 0
 lType1B: lDW[3,12] + lDW[3,13] + lDW[3,14] + lDW[3,15] - 4 lAC[3,3]
   + lT2CanB[3,3] = 0
 lType2C: - lcan[4,12] - lcan[4,13] - lcan[4,14] - lcan[4,15]
   + lT2CanC[4,3] = 0
 lType2A: - lDZ[4,0] - lDZ[4,1] - lDZ[4,2] - lDZ[4,3] + lT2CanA[4,0] = 0
 lType1B: lDW[4,0] + lDW[4,1] + lDW[4,2] + lDW[4,3] - 4 lAC[4,0]
   + lT2CanB[4,0] = 0
 lType2C: - lcan[5,0] - lcan[5,1] - lcan[5,2] - lcan[5,3] + lT2CanC[5,0]
   = 0
 lType2A: - lDZ[4,4] - lDZ[4,5] - lDZ[4,6] - lDZ[4,7] + lT2CanA[4,1] = 0
 lType1B: lDW[4,4] + lDW[4,5] + lDW[4,6] + lDW[4,7] - 4 lAC[4,1]
   + lT2CanB[4,1] = 0
 lType2C: - lcan[5,4] - lcan[5,5] - lcan[5,6] - lcan[5,7] + lT2CanC[5,1]
   = 0
 lType2A: - lDZ[4,8] - lDZ[4,9] - lDZ[4,10] - lDZ[4,11] + lT2CanA[4,2] = 0
 lType1B: lDW[4,8] + lDW[4,9] + lDW[4,10] + lDW[4,11] - 4 lAC[4,2]
   + lT2CanB[4,2] = 0
 lType2C: - lcan[5,8] - lcan[5,9] - lcan[5,10] - lcan[5,11] + lT2CanC[5,2]
   = 0
 lType2A: - lDZ[4,12] - lDZ[4,13] - lDZ[4,14] - lDZ[4,15] + lT2CanA[4,3]
   = 0
 lType1B: lDW[4,12] + lDW[4,13] + lDW[4,14] + lDW[4,15] - 4 lAC[4,3]
   + lT2CanB[4,3] = 0
 lType2C: - lcan[5,12] - lcan[5,13] - lcan[5,14] - lcan[5,15]
   + lT2CanC[5,3] = 0
 lType2A: - lDZ[5,0] - lDZ[5,1] - lDZ[5,2] - lDZ[5,3] + lT2CanA[5,0] = 0
 lType1B: lDW[5,0] + lDW[5,1] + lDW[5,2] + lDW[5,3] - 4 lAC[5,0]
   + lT2CanB[5,0] = 0
 lType2C: - lcan[6,0] - lcan[6,1] - lcan[6,2] - lcan[6,3] + lT2CanC[6,0]
   = 0
 lType2A: - lDZ[5,4] - lDZ[5,5] - lDZ[5,6] - lDZ[5,7] + lT2CanA[5,1] = 0
 lType1B: lDW[5,4] + lDW[5,5] + lDW[5,6] + lDW[5,7] - 4 lAC[5,1]
   + lT2CanB[5,1] = 0
 lType2C: - lcan[6,4] - lcan[6,5] - lcan[6,6] - lcan[6,7] + lT2CanC[6,1]
   = 0
 lType2A: - lDZ[5,8] - lDZ[5,9] - lDZ[5,10] - lDZ[5,11] + lT2CanA[5,2] = 0
 lType1B: lDW[5,8] + lDW[5,9] + lDW[5,10] + lDW[5,11] - 4 lAC[5,2]
   + lT2CanB[5,2] = 0
 lType2C: - lcan[6,8] - lcan[6,9] - lcan[6,10] - lcan[6,11] + lT2CanC[6,2]
   = 0
 lType2A: - lDZ[5,12] - lDZ[5,13] - lDZ[5,14] - lDZ[5,15] + lT2CanA[5,3]
   = 0
 lType1B: lDW[5,12] + lDW[5,13] + lDW[5,14] + lDW[5,15] - 4 lAC[5,3]
   + lT2CanB[5,3] = 0
 lType2C: - lcan[6,12] - lcan[6,13] - lcan[6,14] - lcan[6,15]
   + lT2CanC[6,3] = 0
 lType2A: - lDZ[6,0] - lDZ[6,1] - lDZ[6,2] - lDZ[6,3] + lT2CanA[6,0] = 0
 lType1B: lDW[6,0] + lDW[6,1] + lDW[6,2] + lDW[6,3] - 4 lAC[6,0]
   + lT2CanB[6,0] = 0
 lType2C: - lcan[7,0] - lcan[7,1] - lcan[7,2] - lcan[7,3] + lT2CanC[7,0]
   = 0
 lType2A: - lDZ[6,4] - lDZ[6,5] - lDZ[6,6] - lDZ[6,7] + lT2CanA[6,1] = 0
 lType1B: lDW[6,4] + lDW[6,5] + lDW[6,6] + lDW[6,7] - 4 lAC[6,1]
   + lT2CanB[6,1] = 0
 lType2C: - lcan[7,4] - lcan[7,5] - lcan[7,6] - lcan[7,7] + lT2CanC[7,1]
   = 0
 lType2A: - lDZ[6,8] - lDZ[6,9] - lDZ[6,10] - lDZ[6,11] + lT2CanA[6,2] = 0
 lType1B: lDW[6,8] + lDW[6,9] + lDW[6,10] + lDW[6,11] - 4 lAC[6,2]
   + lT2CanB[6,2] = 0
 lType2C: - lcan[7,8] - lcan[7,9] - lcan[7,10] - lcan[7,11] + lT2CanC[7,2]
   = 0
 lType2A: - lDZ[6,12] - lDZ[6,13] - lDZ[6,14] - lDZ[6,15] + lT2CanA[6,3]
   = 0
 lType1B: lDW[6,12] + lDW[6,13] + lDW[6,14] + lDW[6,15] - 4 lAC[6,3]
   + lT2CanB[6,3] = 0
 lType2C: - lcan[7,12] - lcan[7,13] - lcan[7,14] - lcan[7,15]
   + lT2CanC[7,3] = 0
 lType2A: - lDZ[7,0] - lDZ[7,1] - lDZ[7,2] - lDZ[7,3] + lT2CanA[7,0] = 0
 lType1B: lDW[7,0] + lDW[7,1] + lDW[7,2] + lDW[7,3] - 4 lAC[7,0]
   + lT2CanB[7,0] = 0
 lType2C: - lcan[8,0] - lcan[8,1] - lcan[8,2] - lcan[8,3] + lT2CanC[8,0]
   = 0
 lType2A: - lDZ[7,4] - lDZ[7,5] - lDZ[7,6] - lDZ[7,7] + lT2CanA[7,1] = 0
 lType1B: lDW[7,4] + lDW[7,5] + lDW[7,6] + lDW[7,7] - 4 lAC[7,1]
   + lT2CanB[7,1] = 0
 lType2C: - lcan[8,4] - lcan[8,5] - lcan[8,6] - lcan[8,7] + lT2CanC[8,1]
   = 0
 lType2A: - lDZ[7,8] - lDZ[7,9] - lDZ[7,10] - lDZ[7,11] + lT2CanA[7,2] = 0
 lType1B: lDW[7,8] + lDW[7,9] + lDW[7,10] + lDW[7,11] - 4 lAC[7,2]
   + lT2CanB[7,2] = 0
 lType2C: - lcan[8,8] - lcan[8,9] - lcan[8,10] - lcan[8,11] + lT2CanC[8,2]
   = 0
 lType2A: - lDZ[7,12] - lDZ[7,13] - lDZ[7,14] - lDZ[7,15] + lT2CanA[7,3]
   = 0
 lType1B: lDW[7,12] + lDW[7,13] + lDW[7,14] + lDW[7,15] - 4 lAC[7,3]
   + lT2CanB[7,3] = 0
 lType2C: - lcan[8,12] - lcan[8,13] - lcan[8,14] - lcan[8,15]
   + lT2CanC[8,3] = 0
 lType2CAN: lT2CanA[0,0] - lT2CanB[0,0] - lT2CanC[1,0] + lType2[0,0] >= 0
 lType2CAN: lT2CanA[0,1] - lT2CanB[0,1] - lT2CanC[1,1] + lType2[0,1] >= 0
 lType2CAN: lT2CanA[0,2] - lT2CanB[0,2] - lT2CanC[1,2] + lType2[0,2] >= 0
 lType2CAN: lT2CanA[0,3] - lT2CanB[0,3] - lT2CanC[1,3] + lType2[0,3] >= 0
 lType2CAN: lT2CanA[1,0] - lT2CanB[1,0] - lT2CanC[2,0] + lType2[1,0] >= 0
 lType2CAN: lT2CanA[1,1] - lT2CanB[1,1] - lT2CanC[2,1] + lType2[1,1] >= 0
 lType2CAN: lT2CanA[1,2] - lT2CanB[1,2] - lT2CanC[2,2] + lType2[1,2] >= 0
 lType2CAN: lT2CanA[1,3] - lT2CanB[1,3] - lT2CanC[2,3] + lType2[1,3] >= 0
 lType2CAN: lT2CanA[2,0] - lT2CanB[2,0] - lT2CanC[3,0] + lType2[2,0] >= 0
 lType2CAN: lT2CanA[2,1] - lT2CanB[2,1] - lT2CanC[3,1] + lType2[2,1] >= 0
 lType2CAN: lT2CanA[2,2] - lT2CanB[2,2] - lT2CanC[3,2] + lType2[2,2] >= 0
 lType2CAN: lT2CanA[2,3] - lT2CanB[2,3] - lT2CanC[3,3] + lType2[2,3] >= 0
 lType2CAN: lT2CanA[3,0] - lT2CanB[3,0] - lT2CanC[4,0] + lType2[3,0] >= 0
 lType2CAN: lT2CanA[3,1] - lT2CanB[3,1] - lT2CanC[4,1] + lType2[3,1] >= 0
 lType2CAN: lT2CanA[3,2] - lT2CanB[3,2] - lT2CanC[4,2] + lType2[3,2] >= 0
 lType2CAN: lT2CanA[3,3] - lT2CanB[3,3] - lT2CanC[4,3] + lType2[3,3] >= 0
 lType2CAN: lT2CanA[4,0] - lT2CanB[4,0] - lT2CanC[5,0] + lType2[4,0] >= 0
 lType2CAN: lT2CanA[4,1] - lT2CanB[4,1] - lT2CanC[5,1] + lType2[4,1] >= 0
 lType2CAN: lT2CanA[4,2] - lT2CanB[4,2] - lT2CanC[5,2] + lType2[4,2] >= 0
 lType2CAN: lT2CanA[4,3] - lT2CanB[4,3] - lT2CanC[5,3] + lType2[4,3] >= 0
 lType2CAN: lT2CanA[5,0] - lT2CanB[5,0] - lT2CanC[6,0] + lType2[5,0] >= 0
 lType2CAN: lT2CanA[5,1] - lT2CanB[5,1] - lT2CanC[6,1] + lType2[5,1] >= 0
 lType2CAN: lT2CanA[5,2] - lT2CanB[5,2] - lT2CanC[6,2] + lType2[5,2] >= 0
 lType2CAN: lT2CanA[5,3] - lT2CanB[5,3] - lT2CanC[6,3] + lType2[5,3] >= 0
 lType2CAN: lT2CanA[6,0] - lT2CanB[6,0] - lT2CanC[7,0] + lType2[6,0] >= 0
 lType2CAN: lT2CanA[6,1] - lT2CanB[6,1] - lT2CanC[7,1] + lType2[6,1] >= 0
 lType2CAN: lT2CanA[6,2] - lT2CanB[6,2] - lT2CanC[7,2] + lType2[6,2] >= 0
 lType2CAN: lT2CanA[6,3] - lT2CanB[6,3] - lT2CanC[7,3] + lType2[6,3] >= 0
 lType2CAN: lT2CanA[7,0] - lT2CanB[7,0] - lT2CanC[8,0] + lType2[7,0] >= 0
 lType2CAN: lT2CanA[7,1] - lT2CanB[7,1] - lT2CanC[8,1] + lType2[7,1] >= 0
 lType2CAN: lT2CanA[7,2] - lT2CanB[7,2] - lT2CanC[8,2] + lType2[7,2] >= 0
 lType2CAN: lT2CanA[7,3] - lT2CanB[7,3] - lT2CanC[8,3] + lType2[7,3] >= 0
 lCAN: 2 lLANE[0] + 2 lLANE[1] + 2 lLANE[2] + 2 lLANE[3] + 2 lLANE[4]
   + 2 lLANE[5] + 2 lLANE[6] + 2 lLANE[7] + 2 lLANE[8] + 2 lLANE[9]
   + 2 lLANE[10] + 2 lLANE[11] + 2 lLANE[12] + 2 lLANE[13] + 2 lLANE[14]
   + 2 lLANE[15] - lType1[0] - lType1[1] - lType1[2] - lType1[3]
   - lType1[4] - lType1[5] - lType1[6] - lType1[7] - lType1[8] - lType1[9]
   - lType1[10] - lType1[11] - lType1[12] - lType1[13] - lType1[14]
   - lType1[15] - lType2[0,0] - lType2[0,1] - lType2[0,2] - lType2[0,3]
   - lType2[1,0] - lType2[1,1] - lType2[1,2] - lType2[1,3] - lType2[2,0]
   - lType2[2,1] - lType2[2,2] - lType2[2,3] - lType2[3,0] - lType2[3,1]
   - lType2[3,2] - lType2[3,3] - lType2[4,0] - lType2[4,1] - lType2[4,2]
   - lType2[4,3] - lType2[5,0] - lType2[5,1] - lType2[5,2] - lType2[5,3]
   - lType2[6,0] - lType2[6,1] - lType2[6,2] - lType2[6,3] - lType2[7,0]
   - lType2[7,1] - lType2[7,2] - lType2[7,3] >= 1
 R5596: - lstk[9,0] + leqDX[9,0] >= 0
 R5597: - lstk[9,1] + leqDX[9,1] >= 0
 R5598: - lstk[9,2] + leqDX[9,2] >= 0
 R5599: - lstk[9,3] + leqDX[9,3] >= 0
 R5600: - lstk[9,4] + leqDX[9,4] >= 0
 R5601: - lstk[9,5] + leqDX[9,5] >= 0
 R5602: - lstk[9,6] + leqDX[9,6] >= 0
 R5603: - lstk[9,7] + leqDX[9,7] >= 0
 R5604: - lstk[9,8] + leqDX[9,8] >= 0
 R5605: - lstk[9,9] + leqDX[9,9] >= 0
 R5606: - lstk[9,10] + leqDX[9,10] >= 0
 R5607: - lstk[9,11] + leqDX[9,11] >= 0
 R5608: - lstk[9,12] + leqDX[9,12] >= 0
 R5609: - lstk[9,13] + leqDX[9,13] >= 0
 R5610: - lstk[9,14] + leqDX[9,14] >= 0
 R5611: - lstk[9,15] + leqDX[9,15] >= 0
 lDetAK: lGstk[9,0] - ldetEQY[8,0] >= 0
 lDetAK: - ldetEQY[8,0] + ldetEQZ[8,0] >= 0
 lDetAK: - lGstk[9,0] + ldetEQY[8,0] - ldetEQZ[8,0] >= -1
 lDetAK: lGstk[9,1] - ldetEQY[8,1] >= 0
 lDetAK: - ldetEQY[8,1] + ldetEQZ[8,1] >= 0
 lDetAK: - lGstk[9,1] + ldetEQY[8,1] - ldetEQZ[8,1] >= -1
 lDetAK: lGstk[9,2] - ldetEQY[8,2] >= 0
 lDetAK: - ldetEQY[8,2] + ldetEQZ[8,2] >= 0
 lDetAK: - lGstk[9,2] + ldetEQY[8,2] - ldetEQZ[8,2] >= -1
 lDetAK: lGstk[9,3] - ldetEQY[8,3] >= 0
 lDetAK: - ldetEQY[8,3] + ldetEQZ[8,3] >= 0
 lDetAK: - lGstk[9,3] + ldetEQY[8,3] - ldetEQZ[8,3] >= -1
 lDetAK: lGstk[9,4] - ldetEQY[8,4] >= 0
 lDetAK: - ldetEQY[8,4] + ldetEQZ[8,4] >= 0
 lDetAK: - lGstk[9,4] + ldetEQY[8,4] - ldetEQZ[8,4] >= -1
 lDetAK: lGstk[9,5] - ldetEQY[8,5] >= 0
 lDetAK: - ldetEQY[8,5] + ldetEQZ[8,5] >= 0
 lDetAK: - lGstk[9,5] + ldetEQY[8,5] - ldetEQZ[8,5] >= -1
 lDetAK: lGstk[9,6] - ldetEQY[8,6] >= 0
 lDetAK: - ldetEQY[8,6] + ldetEQZ[8,6] >= 0
 lDetAK: - lGstk[9,6] + ldetEQY[8,6] - ldetEQZ[8,6] >= -1
 lDetAK: lGstk[9,7] - ldetEQY[8,7] >= 0
 lDetAK: - ldetEQY[8,7] + ldetEQZ[8,7] >= 0
 lDetAK: - lGstk[9,7] + ldetEQY[8,7] - ldetEQZ[8,7] >= -1
 lDetAK: lGstk[9,8] - ldetEQY[8,8] >= 0
 lDetAK: - ldetEQY[8,8] + ldetEQZ[8,8] >= 0
 lDetAK: - lGstk[9,8] + ldetEQY[8,8] - ldetEQZ[8,8] >= -1
 lDetAK: lGstk[9,9] - ldetEQY[8,9] >= 0
 lDetAK: - ldetEQY[8,9] + ldetEQZ[8,9] >= 0
 lDetAK: - lGstk[9,9] + ldetEQY[8,9] - ldetEQZ[8,9] >= -1
 lDetAK: lGstk[9,10] - ldetEQY[8,10] >= 0
 lDetAK: - ldetEQY[8,10] + ldetEQZ[8,10] >= 0
 lDetAK: - lGstk[9,10] + ldetEQY[8,10] - ldetEQZ[8,10] >= -1
 lDetAK: lGstk[9,11] - ldetEQY[8,11] >= 0
 lDetAK: - ldetEQY[8,11] + ldetEQZ[8,11] >= 0
 lDetAK: - lGstk[9,11] + ldetEQY[8,11] - ldetEQZ[8,11] >= -1
 lDetAK: lGstk[9,12] - ldetEQY[8,12] >= 0
 lDetAK: - ldetEQY[8,12] + ldetEQZ[8,12] >= 0
 lDetAK: - lGstk[9,12] + ldetEQY[8,12] - ldetEQZ[8,12] >= -1
 lDetAK: lGstk[9,13] - ldetEQY[8,13] >= 0
 lDetAK: - ldetEQY[8,13] + ldetEQZ[8,13] >= 0
 lDetAK: - lGstk[9,13] + ldetEQY[8,13] - ldetEQZ[8,13] >= -1
 lDetAK: lGstk[9,14] - ldetEQY[8,14] >= 0
 lDetAK: - ldetEQY[8,14] + ldetEQZ[8,14] >= 0
 lDetAK: - lGstk[9,14] + ldetEQY[8,14] - ldetEQZ[8,14] >= -1
 lDetAK: lGstk[9,15] - ldetEQY[8,15] >= 0
 lDetAK: - ldetEQY[8,15] + ldetEQZ[8,15] >= 0
 lDetAK: - lGstk[9,15] + ldetEQY[8,15] - ldetEQZ[8,15] >= -1
 lDetAK: lGstk[10,0] - ldetEQY[9,0] >= 0
 lDetAK: - ldetEQY[9,0] + ldetEQZ[9,0] >= 0
 lDetAK: - lGstk[10,0] + ldetEQY[9,0] - ldetEQZ[9,0] >= -1
 lDetAK: lGstk[10,1] - ldetEQY[9,1] >= 0
 lDetAK: - ldetEQY[9,1] + ldetEQZ[9,1] >= 0
 lDetAK: - lGstk[10,1] + ldetEQY[9,1] - ldetEQZ[9,1] >= -1
 lDetAK: lGstk[10,2] - ldetEQY[9,2] >= 0
 lDetAK: - ldetEQY[9,2] + ldetEQZ[9,2] >= 0
 lDetAK: - lGstk[10,2] + ldetEQY[9,2] - ldetEQZ[9,2] >= -1
 lDetAK: lGstk[10,3] - ldetEQY[9,3] >= 0
 lDetAK: - ldetEQY[9,3] + ldetEQZ[9,3] >= 0
 lDetAK: - lGstk[10,3] + ldetEQY[9,3] - ldetEQZ[9,3] >= -1
 lDetAK: lGstk[10,4] - ldetEQY[9,4] >= 0
 lDetAK: - ldetEQY[9,4] + ldetEQZ[9,4] >= 0
 lDetAK: - lGstk[10,4] + ldetEQY[9,4] - ldetEQZ[9,4] >= -1
 lDetAK: lGstk[10,5] - ldetEQY[9,5] >= 0
 lDetAK: - ldetEQY[9,5] + ldetEQZ[9,5] >= 0
 lDetAK: - lGstk[10,5] + ldetEQY[9,5] - ldetEQZ[9,5] >= -1
 lDetAK: lGstk[10,6] - ldetEQY[9,6] >= 0
 lDetAK: - ldetEQY[9,6] + ldetEQZ[9,6] >= 0
 lDetAK: - lGstk[10,6] + ldetEQY[9,6] - ldetEQZ[9,6] >= -1
 lDetAK: lGstk[10,7] - ldetEQY[9,7] >= 0
 lDetAK: - ldetEQY[9,7] + ldetEQZ[9,7] >= 0
 lDetAK: - lGstk[10,7] + ldetEQY[9,7] - ldetEQZ[9,7] >= -1
 lDetAK: lGstk[10,8] - ldetEQY[9,8] >= 0
 lDetAK: - ldetEQY[9,8] + ldetEQZ[9,8] >= 0
 lDetAK: - lGstk[10,8] + ldetEQY[9,8] - ldetEQZ[9,8] >= -1
 lDetAK: lGstk[10,9] - ldetEQY[9,9] >= 0
 lDetAK: - ldetEQY[9,9] + ldetEQZ[9,9] >= 0
 lDetAK: - lGstk[10,9] + ldetEQY[9,9] - ldetEQZ[9,9] >= -1
 lDetAK: lGstk[10,10] - ldetEQY[9,10] >= 0
 lDetAK: - ldetEQY[9,10] + ldetEQZ[9,10] >= 0
 lDetAK: - lGstk[10,10] + ldetEQY[9,10] - ldetEQZ[9,10] >= -1
 lDetAK: lGstk[10,11] - ldetEQY[9,11] >= 0
 lDetAK: - ldetEQY[9,11] + ldetEQZ[9,11] >= 0
 lDetAK: - lGstk[10,11] + ldetEQY[9,11] - ldetEQZ[9,11] >= -1
 lDetAK: lGstk[10,12] - ldetEQY[9,12] >= 0
 lDetAK: - ldetEQY[9,12] + ldetEQZ[9,12] >= 0
 lDetAK: - lGstk[10,12] + ldetEQY[9,12] - ldetEQZ[9,12] >= -1
 lDetAK: lGstk[10,13] - ldetEQY[9,13] >= 0
 lDetAK: - ldetEQY[9,13] + ldetEQZ[9,13] >= 0
 lDetAK: - lGstk[10,13] + ldetEQY[9,13] - ldetEQZ[9,13] >= -1
 lDetAK: lGstk[10,14] - ldetEQY[9,14] >= 0
 lDetAK: - ldetEQY[9,14] + ldetEQZ[9,14] >= 0
 lDetAK: - lGstk[10,14] + ldetEQY[9,14] - ldetEQZ[9,14] >= -1
 lDetAK: lGstk[10,15] - ldetEQY[9,15] >= 0
 lDetAK: - ldetEQY[9,15] + ldetEQZ[9,15] >= 0
 lDetAK: - lGstk[10,15] + ldetEQY[9,15] - ldetEQZ[9,15] >= -1
 lDetSC: ldetEQX[8,0] - ldetEQY[8,0] = 0
 lDetSC: ldetEQX[8,1] - ldetEQY[8,1] = 0
 lDetSC: ldetEQX[8,2] - ldetEQY[8,2] = 0
 lDetSC: ldetEQX[8,3] - ldetEQY[8,3] = 0
 lDetSC: ldetEQX[8,4] - ldetEQY[8,4] = 0
 lDetSC: ldetEQX[8,5] - ldetEQY[8,5] = 0
 lDetSC: ldetEQX[8,6] - ldetEQY[8,6] = 0
 lDetSC: ldetEQX[8,7] - ldetEQY[8,7] = 0
 lDetSC: ldetEQX[8,8] - ldetEQY[8,8] = 0
 lDetSC: ldetEQX[8,9] - ldetEQY[8,9] = 0
 lDetSC: ldetEQX[8,10] - ldetEQY[8,10] = 0
 lDetSC: ldetEQX[8,11] - ldetEQY[8,11] = 0
 lDetSC: ldetEQX[8,12] - ldetEQY[8,12] = 0
 lDetSC: ldetEQX[8,13] - ldetEQY[8,13] = 0
 lDetSC: ldetEQX[8,14] - ldetEQY[8,14] = 0
 lDetSC: ldetEQX[8,15] - ldetEQY[8,15] = 0
 lDetSC: ldetEQX[9,0] - ldetEQY[9,0] = 0
 lDetSC: ldetEQX[9,1] - ldetEQY[9,1] = 0
 lDetSC: ldetEQX[9,2] - ldetEQY[9,2] = 0
 lDetSC: ldetEQX[9,3] - ldetEQY[9,3] = 0
 lDetSC: ldetEQX[9,4] - ldetEQY[9,4] = 0
 lDetSC: ldetEQX[9,5] - ldetEQY[9,5] = 0
 lDetSC: ldetEQX[9,6] - ldetEQY[9,6] = 0
 lDetSC: ldetEQX[9,7] - ldetEQY[9,7] = 0
 lDetSC: ldetEQX[9,8] - ldetEQY[9,8] = 0
 lDetSC: ldetEQX[9,9] - ldetEQY[9,9] = 0
 lDetSC: ldetEQX[9,10] - ldetEQY[9,10] = 0
 lDetSC: ldetEQX[9,11] - ldetEQY[9,11] = 0
 lDetSC: ldetEQX[9,12] - ldetEQY[9,12] = 0
 lDetSC: ldetEQX[9,13] - ldetEQY[9,13] = 0
 lDetSC: ldetEQX[9,14] - ldetEQY[9,14] = 0
 lDetSC: ldetEQX[9,15] - ldetEQY[9,15] = 0
 lDetMC: - ldetEQX[9,0] + lDetC[9,0] <= 0
 lDetMC: - ldetEQX[9,1] + lDetC[9,0] <= 0
 lDetMC: - ldetEQX[9,2] + lDetC[9,0] <= 0
 lDetMC: - ldetEQX[9,3] + lDetC[9,0] <= 0
 lDetMC: ldetEQX[9,0] + ldetEQX[9,1] + ldetEQX[9,2] + ldetEQX[9,3]
   - lDetC[9,0] <= 3
 lDetMC: - ldetEQW[8,0] + lDetC[9,0] = 0
 lDetMC: - ldetEQW[8,1] + lDetC[9,0] = 0
 lDetMC: - ldetEQW[8,2] + lDetC[9,0] = 0
 lDetMC: - ldetEQW[8,3] + lDetC[9,0] = 0
 lDetMC: - ldetEQX[9,4] + lDetC[9,1] <= 0
 lDetMC: - ldetEQX[9,5] + lDetC[9,1] <= 0
 lDetMC: - ldetEQX[9,6] + lDetC[9,1] <= 0
 lDetMC: - ldetEQX[9,7] + lDetC[9,1] <= 0
 lDetMC: ldetEQX[9,4] + ldetEQX[9,5] + ldetEQX[9,6] + ldetEQX[9,7]
   - lDetC[9,1] <= 3
 lDetMC: - ldetEQW[8,4] + lDetC[9,1] = 0
 lDetMC: - ldetEQW[8,5] + lDetC[9,1] = 0
 lDetMC: - ldetEQW[8,6] + lDetC[9,1] = 0
 lDetMC: - ldetEQW[8,7] + lDetC[9,1] = 0
 lDetMC: - ldetEQX[9,8] + lDetC[9,2] <= 0
 lDetMC: - ldetEQX[9,9] + lDetC[9,2] <= 0
 lDetMC: - ldetEQX[9,10] + lDetC[9,2] <= 0
 lDetMC: - ldetEQX[9,11] + lDetC[9,2] <= 0
 lDetMC: ldetEQX[9,8] + ldetEQX[9,9] + ldetEQX[9,10] + ldetEQX[9,11]
   - lDetC[9,2] <= 3
 lDetMC: - ldetEQW[8,8] + lDetC[9,2] = 0
 lDetMC: - ldetEQW[8,9] + lDetC[9,2] = 0
 lDetMC: - ldetEQW[8,10] + lDetC[9,2] = 0
 lDetMC: - ldetEQW[8,11] + lDetC[9,2] = 0
 lDetMC: - ldetEQX[9,12] + lDetC[9,3] <= 0
 lDetMC: - ldetEQX[9,13] + lDetC[9,3] <= 0
 lDetMC: - ldetEQX[9,14] + lDetC[9,3] <= 0
 lDetMC: - ldetEQX[9,15] + lDetC[9,3] <= 0
 lDetMC: ldetEQX[9,12] + ldetEQX[9,13] + ldetEQX[9,14] + ldetEQX[9,15]
   - lDetC[9,3] <= 3
 lDetMC: - ldetEQW[8,12] + lDetC[9,3] = 0
 lDetMC: - ldetEQW[8,13] + lDetC[9,3] = 0
 lDetMC: - ldetEQW[8,14] + lDetC[9,3] = 0
 lDetMC: - ldetEQW[8,15] + lDetC[9,3] = 0
 lDetSR: ldetEQZ[8,0] - ldetEQW[8,0] = 0
 lDetSR: ldetEQZ[8,5] - ldetEQW[8,1] = 0
 lDetSR: ldetEQZ[8,10] - ldetEQW[8,2] = 0
 lDetSR: ldetEQZ[8,15] - ldetEQW[8,3] = 0
 lDetSR: ldetEQZ[8,4] - ldetEQW[8,4] = 0
 lDetSR: ldetEQZ[8,9] - ldetEQW[8,5] = 0
 lDetSR: ldetEQZ[8,14] - ldetEQW[8,6] = 0
 lDetSR: ldetEQZ[8,3] - ldetEQW[8,7] = 0
 lDetSR: ldetEQZ[8,8] - ldetEQW[8,8] = 0
 lDetSR: ldetEQZ[8,13] - ldetEQW[8,9] = 0
 lDetSR: ldetEQZ[8,2] - ldetEQW[8,10] = 0
 lDetSR: ldetEQZ[8,7] - ldetEQW[8,11] = 0
 lDetSR: ldetEQZ[8,12] - ldetEQW[8,12] = 0
 lDetSR: ldetEQZ[8,1] - ldetEQW[8,13] = 0
 lDetSR: ldetEQZ[8,6] - ldetEQW[8,14] = 0
 lDetSR: ldetEQZ[8,11] - ldetEQW[8,15] = 0
 lFrSC: lDX[8,0] - ldx[8,0] - lFrSB[8,0] >= 0
 lFrSC: ldetEQX[8,0] - lFrSB[8,0] >= 0
 lFrSC: - lDX[8,0] + ldx[8,0] - ldetEQX[8,0] + lFrSB[8,0] >= -1
 lFrSC: lDX[8,1] - ldx[8,1] - lFrSB[8,1] >= 0
 lFrSC: ldetEQX[8,1] - lFrSB[8,1] >= 0
 lFrSC: - lDX[8,1] + ldx[8,1] - ldetEQX[8,1] + lFrSB[8,1] >= -1
 lFrSC: lDX[8,2] - ldx[8,2] - lFrSB[8,2] >= 0
 lFrSC: ldetEQX[8,2] - lFrSB[8,2] >= 0
 lFrSC: - lDX[8,2] + ldx[8,2] - ldetEQX[8,2] + lFrSB[8,2] >= -1
 lFrSC: lDX[8,3] - ldx[8,3] - lFrSB[8,3] >= 0
 lFrSC: ldetEQX[8,3] - lFrSB[8,3] >= 0
 lFrSC: - lDX[8,3] + ldx[8,3] - ldetEQX[8,3] + lFrSB[8,3] >= -1
 lFrSC: lDX[8,4] - ldx[8,4] - lFrSB[8,4] >= 0
 lFrSC: ldetEQX[8,4] - lFrSB[8,4] >= 0
 lFrSC: - lDX[8,4] + ldx[8,4] - ldetEQX[8,4] + lFrSB[8,4] >= -1
 lFrSC: lDX[8,5] - ldx[8,5] - lFrSB[8,5] >= 0
 lFrSC: ldetEQX[8,5] - lFrSB[8,5] >= 0
 lFrSC: - lDX[8,5] + ldx[8,5] - ldetEQX[8,5] + lFrSB[8,5] >= -1
 lFrSC: lDX[8,6] - ldx[8,6] - lFrSB[8,6] >= 0
 lFrSC: ldetEQX[8,6] - lFrSB[8,6] >= 0
 lFrSC: - lDX[8,6] + ldx[8,6] - ldetEQX[8,6] + lFrSB[8,6] >= -1
 lFrSC: lDX[8,7] - ldx[8,7] - lFrSB[8,7] >= 0
 lFrSC: ldetEQX[8,7] - lFrSB[8,7] >= 0
 lFrSC: - lDX[8,7] + ldx[8,7] - ldetEQX[8,7] + lFrSB[8,7] >= -1
 lFrSC: lDX[8,8] - ldx[8,8] - lFrSB[8,8] >= 0
 lFrSC: ldetEQX[8,8] - lFrSB[8,8] >= 0
 lFrSC: - lDX[8,8] + ldx[8,8] - ldetEQX[8,8] + lFrSB[8,8] >= -1
 lFrSC: lDX[8,9] - ldx[8,9] - lFrSB[8,9] >= 0
 lFrSC: ldetEQX[8,9] - lFrSB[8,9] >= 0
 lFrSC: - lDX[8,9] + ldx[8,9] - ldetEQX[8,9] + lFrSB[8,9] >= -1
 lFrSC: lDX[8,10] - ldx[8,10] - lFrSB[8,10] >= 0
 lFrSC: ldetEQX[8,10] - lFrSB[8,10] >= 0
 lFrSC: - lDX[8,10] + ldx[8,10] - ldetEQX[8,10] + lFrSB[8,10] >= -1
 lFrSC: lDX[8,11] - ldx[8,11] - lFrSB[8,11] >= 0
 lFrSC: ldetEQX[8,11] - lFrSB[8,11] >= 0
 lFrSC: - lDX[8,11] + ldx[8,11] - ldetEQX[8,11] + lFrSB[8,11] >= -1
 lFrSC: lDX[8,12] - ldx[8,12] - lFrSB[8,12] >= 0
 lFrSC: ldetEQX[8,12] - lFrSB[8,12] >= 0
 lFrSC: - lDX[8,12] + ldx[8,12] - ldetEQX[8,12] + lFrSB[8,12] >= -1
 lFrSC: lDX[8,13] - ldx[8,13] - lFrSB[8,13] >= 0
 lFrSC: ldetEQX[8,13] - lFrSB[8,13] >= 0
 lFrSC: - lDX[8,13] + ldx[8,13] - ldetEQX[8,13] + lFrSB[8,13] >= -1
 lFrSC: lDX[8,14] - ldx[8,14] - lFrSB[8,14] >= 0
 lFrSC: ldetEQX[8,14] - lFrSB[8,14] >= 0
 lFrSC: - lDX[8,14] + ldx[8,14] - ldetEQX[8,14] + lFrSB[8,14] >= -1
 lFrSC: lDX[8,15] - ldx[8,15] - lFrSB[8,15] >= 0
 lFrSC: ldetEQX[8,15] - lFrSB[8,15] >= 0
 lFrSC: - lDX[8,15] + ldx[8,15] - ldetEQX[8,15] + lFrSB[8,15] >= -1
 lFrSC: leqDX[9,0] - leqdx[9,0] - lFrSB[9,0] >= 0
 lFrSC: ldetEQX[9,0] - lFrSB[9,0] >= 0
 lFrSC: - leqDX[9,0] + leqdx[9,0] - ldetEQX[9,0] + lFrSB[9,0] >= -1
 lFrSC: leqDX[9,1] - leqdx[9,1] - lFrSB[9,1] >= 0
 lFrSC: ldetEQX[9,1] - lFrSB[9,1] >= 0
 lFrSC: - leqDX[9,1] + leqdx[9,1] - ldetEQX[9,1] + lFrSB[9,1] >= -1
 lFrSC: leqDX[9,2] - leqdx[9,2] - lFrSB[9,2] >= 0
 lFrSC: ldetEQX[9,2] - lFrSB[9,2] >= 0
 lFrSC: - leqDX[9,2] + leqdx[9,2] - ldetEQX[9,2] + lFrSB[9,2] >= -1
 lFrSC: leqDX[9,3] - leqdx[9,3] - lFrSB[9,3] >= 0
 lFrSC: ldetEQX[9,3] - lFrSB[9,3] >= 0
 lFrSC: - leqDX[9,3] + leqdx[9,3] - ldetEQX[9,3] + lFrSB[9,3] >= -1
 lFrSC: leqDX[9,4] - leqdx[9,4] - lFrSB[9,4] >= 0
 lFrSC: ldetEQX[9,4] - lFrSB[9,4] >= 0
 lFrSC: - leqDX[9,4] + leqdx[9,4] - ldetEQX[9,4] + lFrSB[9,4] >= -1
 lFrSC: leqDX[9,5] - leqdx[9,5] - lFrSB[9,5] >= 0
 lFrSC: ldetEQX[9,5] - lFrSB[9,5] >= 0
 lFrSC: - leqDX[9,5] + leqdx[9,5] - ldetEQX[9,5] + lFrSB[9,5] >= -1
 lFrSC: leqDX[9,6] - leqdx[9,6] - lFrSB[9,6] >= 0
 lFrSC: ldetEQX[9,6] - lFrSB[9,6] >= 0
 lFrSC: - leqDX[9,6] + leqdx[9,6] - ldetEQX[9,6] + lFrSB[9,6] >= -1
 lFrSC: leqDX[9,7] - leqdx[9,7] - lFrSB[9,7] >= 0
 lFrSC: ldetEQX[9,7] - lFrSB[9,7] >= 0
 lFrSC: - leqDX[9,7] + leqdx[9,7] - ldetEQX[9,7] + lFrSB[9,7] >= -1
 lFrSC: leqDX[9,8] - leqdx[9,8] - lFrSB[9,8] >= 0
 lFrSC: ldetEQX[9,8] - lFrSB[9,8] >= 0
 lFrSC: - leqDX[9,8] + leqdx[9,8] - ldetEQX[9,8] + lFrSB[9,8] >= -1
 lFrSC: leqDX[9,9] - leqdx[9,9] - lFrSB[9,9] >= 0
 lFrSC: ldetEQX[9,9] - lFrSB[9,9] >= 0
 lFrSC: - leqDX[9,9] + leqdx[9,9] - ldetEQX[9,9] + lFrSB[9,9] >= -1
 lFrSC: leqDX[9,10] - leqdx[9,10] - lFrSB[9,10] >= 0
 lFrSC: ldetEQX[9,10] - lFrSB[9,10] >= 0
 lFrSC: - leqDX[9,10] + leqdx[9,10] - ldetEQX[9,10] + lFrSB[9,10] >= -1
 lFrSC: leqDX[9,11] - leqdx[9,11] - lFrSB[9,11] >= 0
 lFrSC: ldetEQX[9,11] - lFrSB[9,11] >= 0
 lFrSC: - leqDX[9,11] + leqdx[9,11] - ldetEQX[9,11] + lFrSB[9,11] >= -1
 lFrSC: leqDX[9,12] - leqdx[9,12] - lFrSB[9,12] >= 0
 lFrSC: ldetEQX[9,12] - lFrSB[9,12] >= 0
 lFrSC: - leqDX[9,12] + leqdx[9,12] - ldetEQX[9,12] + lFrSB[9,12] >= -1
 lFrSC: leqDX[9,13] - leqdx[9,13] - lFrSB[9,13] >= 0
 lFrSC: ldetEQX[9,13] - lFrSB[9,13] >= 0
 lFrSC: - leqDX[9,13] + leqdx[9,13] - ldetEQX[9,13] + lFrSB[9,13] >= -1
 lFrSC: leqDX[9,14] - leqdx[9,14] - lFrSB[9,14] >= 0
 lFrSC: ldetEQX[9,14] - lFrSB[9,14] >= 0
 lFrSC: - leqDX[9,14] + leqdx[9,14] - ldetEQX[9,14] + lFrSB[9,14] >= -1
 lFrSC: leqDX[9,15] - leqdx[9,15] - lFrSB[9,15] >= 0
 lFrSC: ldetEQX[9,15] - lFrSB[9,15] >= 0
 lFrSC: - leqDX[9,15] + leqdx[9,15] - ldetEQX[9,15] + lFrSB[9,15] >= -1
 lFrMC: - ldetEQX[9,0] + lDiffDetMC[8,0] >= 0
 lFrMC: - leqDX[9,0] + ldetEQX[9,0] - lDiffDetMC[8,0] >= -1
 lFrMC: leqDX[9,0] + lDiffDetMC[8,0] >= 1
 lFrMC: - ldetEQX[9,1] + lDiffDetMC[8,1] >= 0
 lFrMC: - leqDX[9,1] + ldetEQX[9,1] - lDiffDetMC[8,1] >= -1
 lFrMC: leqDX[9,1] + lDiffDetMC[8,1] >= 1
 lFrMC: - ldetEQX[9,2] + lDiffDetMC[8,2] >= 0
 lFrMC: - leqDX[9,2] + ldetEQX[9,2] - lDiffDetMC[8,2] >= -1
 lFrMC: leqDX[9,2] + lDiffDetMC[8,2] >= 1
 lFrMC: - ldetEQX[9,3] + lDiffDetMC[8,3] >= 0
 lFrMC: - leqDX[9,3] + ldetEQX[9,3] - lDiffDetMC[8,3] >= -1
 lFrMC: leqDX[9,3] + lDiffDetMC[8,3] >= 1
 lFrMC: - ldetEQX[9,4] + lDiffDetMC[8,4] >= 0
 lFrMC: - leqDX[9,4] + ldetEQX[9,4] - lDiffDetMC[8,4] >= -1
 lFrMC: leqDX[9,4] + lDiffDetMC[8,4] >= 1
 lFrMC: - ldetEQX[9,5] + lDiffDetMC[8,5] >= 0
 lFrMC: - leqDX[9,5] + ldetEQX[9,5] - lDiffDetMC[8,5] >= -1
 lFrMC: leqDX[9,5] + lDiffDetMC[8,5] >= 1
 lFrMC: - ldetEQX[9,6] + lDiffDetMC[8,6] >= 0
 lFrMC: - leqDX[9,6] + ldetEQX[9,6] - lDiffDetMC[8,6] >= -1
 lFrMC: leqDX[9,6] + lDiffDetMC[8,6] >= 1
 lFrMC: - ldetEQX[9,7] + lDiffDetMC[8,7] >= 0
 lFrMC: - leqDX[9,7] + ldetEQX[9,7] - lDiffDetMC[8,7] >= -1
 lFrMC: leqDX[9,7] + lDiffDetMC[8,7] >= 1
 lFrMC: - ldetEQX[9,8] + lDiffDetMC[8,8] >= 0
 lFrMC: - leqDX[9,8] + ldetEQX[9,8] - lDiffDetMC[8,8] >= -1
 lFrMC: leqDX[9,8] + lDiffDetMC[8,8] >= 1
 lFrMC: - ldetEQX[9,9] + lDiffDetMC[8,9] >= 0
 lFrMC: - leqDX[9,9] + ldetEQX[9,9] - lDiffDetMC[8,9] >= -1
 lFrMC: leqDX[9,9] + lDiffDetMC[8,9] >= 1
 lFrMC: - ldetEQX[9,10] + lDiffDetMC[8,10] >= 0
 lFrMC: - leqDX[9,10] + ldetEQX[9,10] - lDiffDetMC[8,10] >= -1
 lFrMC: leqDX[9,10] + lDiffDetMC[8,10] >= 1
 lFrMC: - ldetEQX[9,11] + lDiffDetMC[8,11] >= 0
 lFrMC: - leqDX[9,11] + ldetEQX[9,11] - lDiffDetMC[8,11] >= -1
 lFrMC: leqDX[9,11] + lDiffDetMC[8,11] >= 1
 lFrMC: - ldetEQX[9,12] + lDiffDetMC[8,12] >= 0
 lFrMC: - leqDX[9,12] + ldetEQX[9,12] - lDiffDetMC[8,12] >= -1
 lFrMC: leqDX[9,12] + lDiffDetMC[8,12] >= 1
 lFrMC: - ldetEQX[9,13] + lDiffDetMC[8,13] >= 0
 lFrMC: - leqDX[9,13] + ldetEQX[9,13] - lDiffDetMC[8,13] >= -1
 lFrMC: leqDX[9,13] + lDiffDetMC[8,13] >= 1
 lFrMC: - ldetEQX[9,14] + lDiffDetMC[8,14] >= 0
 lFrMC: - leqDX[9,14] + ldetEQX[9,14] - lDiffDetMC[8,14] >= -1
 lFrMC: leqDX[9,14] + lDiffDetMC[8,14] >= 1
 lFrMC: - ldetEQX[9,15] + lDiffDetMC[8,15] >= 0
 lFrMC: - leqDX[9,15] + ldetEQX[9,15] - lDiffDetMC[8,15] >= -1
 lFrMC: leqDX[9,15] + lDiffDetMC[8,15] >= 1
 lFrMC: - lDiffDetMC[8,0] + lDiffDetC[8,0] <= 0
 lFrMC: - lDiffDetMC[8,1] + lDiffDetC[8,0] <= 0
 lFrMC: - lDiffDetMC[8,2] + lDiffDetC[8,0] <= 0
 lFrMC: - lDiffDetMC[8,3] + lDiffDetC[8,0] <= 0
 lFrMC: lDiffDetMC[8,0] + lDiffDetMC[8,1] + lDiffDetMC[8,2]
   + lDiffDetMC[8,3] - lDiffDetC[8,0] <= 3
 lFrMC: - lDiffDetMC[8,4] + lDiffDetC[8,1] <= 0
 lFrMC: - lDiffDetMC[8,5] + lDiffDetC[8,1] <= 0
 lFrMC: - lDiffDetMC[8,6] + lDiffDetC[8,1] <= 0
 lFrMC: - lDiffDetMC[8,7] + lDiffDetC[8,1] <= 0
 lFrMC: lDiffDetMC[8,4] + lDiffDetMC[8,5] + lDiffDetMC[8,6]
   + lDiffDetMC[8,7] - lDiffDetC[8,1] <= 3
 lFrMC: - lDiffDetMC[8,8] + lDiffDetC[8,2] <= 0
 lFrMC: - lDiffDetMC[8,9] + lDiffDetC[8,2] <= 0
 lFrMC: - lDiffDetMC[8,10] + lDiffDetC[8,2] <= 0
 lFrMC: - lDiffDetMC[8,11] + lDiffDetC[8,2] <= 0
 lFrMC: lDiffDetMC[8,8] + lDiffDetMC[8,9] + lDiffDetMC[8,10]
   + lDiffDetMC[8,11] - lDiffDetC[8,2] <= 3
 lFrMC: - lDiffDetMC[8,12] + lDiffDetC[8,3] <= 0
 lFrMC: - lDiffDetMC[8,13] + lDiffDetC[8,3] <= 0
 lFrMC: - lDiffDetMC[8,14] + lDiffDetC[8,3] <= 0
 lFrMC: - lDiffDetMC[8,15] + lDiffDetC[8,3] <= 0
 lFrMC: lDiffDetMC[8,12] + lDiffDetMC[8,13] + lDiffDetMC[8,14]
   + lDiffDetMC[8,15] - lDiffDetC[8,3] <= 3
 lFreqMC: lTC[8,0] - lFrMC[8,0] >= 0
 lFreqMC: leqdw[8,0] + lFrMC[8,0] <= 1
 lFreqMC: lTC[8,0] - leqdw[8,0] + lDiffDetC[8,0] - lFrMC[8,0] <= 1
 lFreqMC: lDiffDetC[8,0] - lFrMC[8,0] >= 0
 lFreqMC: lTC[8,0] - lFrMC[8,1] >= 0
 lFreqMC: leqdw[8,1] + lFrMC[8,1] <= 1
 lFreqMC: lTC[8,0] - leqdw[8,1] + lDiffDetC[8,0] - lFrMC[8,1] <= 1
 lFreqMC: lDiffDetC[8,0] - lFrMC[8,1] >= 0
 lFreqMC: lTC[8,0] - lFrMC[8,2] >= 0
 lFreqMC: leqdw[8,2] + lFrMC[8,2] <= 1
 lFreqMC: lTC[8,0] - leqdw[8,2] + lDiffDetC[8,0] - lFrMC[8,2] <= 1
 lFreqMC: lDiffDetC[8,0] - lFrMC[8,2] >= 0
 lFreqMC: lTC[8,0] - lFrMC[8,3] >= 0
 lFreqMC: leqdw[8,3] + lFrMC[8,3] <= 1
 lFreqMC: lTC[8,0] - leqdw[8,3] + lDiffDetC[8,0] - lFrMC[8,3] <= 1
 lFreqMC: lDiffDetC[8,0] - lFrMC[8,3] >= 0
 lFreqMC: lTC[8,1] - lFrMC[8,4] >= 0
 lFreqMC: leqdw[8,4] + lFrMC[8,4] <= 1
 lFreqMC: lTC[8,1] - leqdw[8,4] + lDiffDetC[8,1] - lFrMC[8,4] <= 1
 lFreqMC: lDiffDetC[8,1] - lFrMC[8,4] >= 0
 lFreqMC: lTC[8,1] - lFrMC[8,5] >= 0
 lFreqMC: leqdw[8,5] + lFrMC[8,5] <= 1
 lFreqMC: lTC[8,1] - leqdw[8,5] + lDiffDetC[8,1] - lFrMC[8,5] <= 1
 lFreqMC: lDiffDetC[8,1] - lFrMC[8,5] >= 0
 lFreqMC: lTC[8,1] - lFrMC[8,6] >= 0
 lFreqMC: leqdw[8,6] + lFrMC[8,6] <= 1
 lFreqMC: lTC[8,1] - leqdw[8,6] + lDiffDetC[8,1] - lFrMC[8,6] <= 1
 lFreqMC: lDiffDetC[8,1] - lFrMC[8,6] >= 0
 lFreqMC: lTC[8,1] - lFrMC[8,7] >= 0
 lFreqMC: leqdw[8,7] + lFrMC[8,7] <= 1
 lFreqMC: lTC[8,1] - leqdw[8,7] + lDiffDetC[8,1] - lFrMC[8,7] <= 1
 lFreqMC: lDiffDetC[8,1] - lFrMC[8,7] >= 0
 lFreqMC: lTC[8,2] - lFrMC[8,8] >= 0
 lFreqMC: leqdw[8,8] + lFrMC[8,8] <= 1
 lFreqMC: lTC[8,2] - leqdw[8,8] + lDiffDetC[8,2] - lFrMC[8,8] <= 1
 lFreqMC: lDiffDetC[8,2] - lFrMC[8,8] >= 0
 lFreqMC: lTC[8,2] - lFrMC[8,9] >= 0
 lFreqMC: leqdw[8,9] + lFrMC[8,9] <= 1
 lFreqMC: lTC[8,2] - leqdw[8,9] + lDiffDetC[8,2] - lFrMC[8,9] <= 1
 lFreqMC: lDiffDetC[8,2] - lFrMC[8,9] >= 0
 lFreqMC: lTC[8,2] - lFrMC[8,10] >= 0
 lFreqMC: leqdw[8,10] + lFrMC[8,10] <= 1
 lFreqMC: lTC[8,2] - leqdw[8,10] + lDiffDetC[8,2] - lFrMC[8,10] <= 1
 lFreqMC: lDiffDetC[8,2] - lFrMC[8,10] >= 0
 lFreqMC: lTC[8,2] - lFrMC[8,11] >= 0
 lFreqMC: leqdw[8,11] + lFrMC[8,11] <= 1
 lFreqMC: lTC[8,2] - leqdw[8,11] + lDiffDetC[8,2] - lFrMC[8,11] <= 1
 lFreqMC: lDiffDetC[8,2] - lFrMC[8,11] >= 0
 lFreqMC: lTC[8,3] - lFrMC[8,12] >= 0
 lFreqMC: leqdw[8,12] + lFrMC[8,12] <= 1
 lFreqMC: lTC[8,3] - leqdw[8,12] + lDiffDetC[8,3] - lFrMC[8,12] <= 1
 lFreqMC: lDiffDetC[8,3] - lFrMC[8,12] >= 0
 lFreqMC: lTC[8,3] - lFrMC[8,13] >= 0
 lFreqMC: leqdw[8,13] + lFrMC[8,13] <= 1
 lFreqMC: lTC[8,3] - leqdw[8,13] + lDiffDetC[8,3] - lFrMC[8,13] <= 1
 lFreqMC: lDiffDetC[8,3] - lFrMC[8,13] >= 0
 lFreqMC: lTC[8,3] - lFrMC[8,14] >= 0
 lFreqMC: leqdw[8,14] + lFrMC[8,14] <= 1
 lFreqMC: lTC[8,3] - leqdw[8,14] + lDiffDetC[8,3] - lFrMC[8,14] <= 1
 lFreqMC: lDiffDetC[8,3] - lFrMC[8,14] >= 0
 lFreqMC: lTC[8,3] - lFrMC[8,15] >= 0
 lFreqMC: leqdw[8,15] + lFrMC[8,15] <= 1
 lFreqMC: lTC[8,3] - leqdw[8,15] + lDiffDetC[8,3] - lFrMC[8,15] <= 1
 lFreqMC: lDiffDetC[8,3] - lFrMC[8,15] >= 0
 R6020: uDetW[-1,0] = 1
 R6021: ldetEQZ[9,0] = 1
 R6022: uDetW[-1,1] = 1
 R6023: ldetEQZ[9,1] = 1
 R6024: uDetW[-1,2] = 1
 R6025: ldetEQZ[9,2] = 1
 R6026: uDetW[-1,3] = 1
 R6027: ldetEQZ[9,3] = 1
 R6028: uDetW[-1,4] = 1
 R6029: ldetEQZ[9,4] = 1
 R6030: uDetW[-1,5] = 1
 R6031: ldetEQZ[9,5] = 1
 R6032: uDetW[-1,6] = 1
 R6033: ldetEQZ[9,6] = 1
 R6034: uDetW[-1,7] = 1
 R6035: ldetEQZ[9,7] = 1
 R6036: uDetW[-1,8] = 1
 R6037: ldetEQZ[9,8] = 1
 R6038: uDetW[-1,9] = 1
 R6039: ldetEQZ[9,9] = 1
 R6040: uDetW[-1,10] = 1
 R6041: ldetEQZ[9,10] = 1
 R6042: uDetW[-1,11] = 1
 R6043: ldetEQZ[9,11] = 1
 R6044: uDetW[-1,12] = 1
 R6045: ldetEQZ[9,12] = 1
 R6046: uDetW[-1,13] = 1
 R6047: ldetEQZ[9,13] = 1
 R6048: uDetW[-1,14] = 1
 R6049: ldetEQZ[9,14] = 1
 R6050: uDetW[-1,15] = 1
 R6051: ldetEQZ[9,15] = 1
 R6052: - 8 udw[-1,0] - 8 udw[-1,1] - 8 udw[-1,2] - 8 udw[-1,3]
   - 8 udw[-1,4] - 8 udw[-1,5] - 8 udw[-1,6] - 8 udw[-1,7] - 8 udw[-1,8]
   - 8 udw[-1,9] - 8 udw[-1,10] - 8 udw[-1,11] - 8 udw[-1,12]
   - 8 udw[-1,13] - 8 udw[-1,14] - 8 udw[-1,15] + rb = 0
 R6053: 8 udx[1,0] + 8 udx[1,1] + 8 udx[1,2] + 8 udx[1,3] + 8 udx[1,4]
   + 8 udx[1,5] + 8 udx[1,6] + 8 udx[1,7] + 8 udx[1,8] + 8 udx[1,9]
   + 8 udx[1,10] + 8 udx[1,11] + 8 udx[1,12] + 8 udx[1,13] + 8 udx[1,14]
   + 8 udx[1,15] - rb + cb = 0
 R6054: - 8 udx[0,0] - 8 udx[0,1] - 8 udx[0,2] - 8 udx[0,3] - 8 udx[0,4]
   - 8 udx[0,5] - 8 udx[0,6] - 8 udx[0,7] - 8 udx[0,8] - 8 udx[0,9]
   - 8 udx[0,10] - 8 udx[0,11] - 8 udx[0,12] - 8 udx[0,13] - 8 udx[0,14]
   - 8 udx[0,15] + mb = 0
 R6055: - 8 uGstk[0,0] - 8 uGstk[0,1] - 8 uGstk[0,2] - 8 uGstk[0,3]
   - 8 uGstk[0,4] - 8 uGstk[0,5] - 8 uGstk[0,6] - 8 uGstk[0,7]
   - 8 uGstk[0,8] - 8 uGstk[0,9] - 8 uGstk[0,10] - 8 uGstk[0,11]
   - 8 uGstk[0,12] - 8 uGstk[0,13] - 8 uGstk[0,14] - 8 uGstk[0,15] + mb'
   = 0
 R6056: - 8 uFrSB[0,0] - 8 uFrSB[0,1] - 8 uFrSB[0,2] - 8 uFrSB[0,3]
   - 8 uFrSB[0,4] - 8 uFrSB[0,5] - 8 uFrSB[0,6] - 8 uFrSB[0,7]
   - 8 uFrSB[0,8] - 8 uFrSB[0,9] - 8 uFrSB[0,10] - 8 uFrSB[0,11]
   - 8 uFrSB[0,12] - 8 uFrSB[0,13] - 8 uFrSB[0,14] - 8 uFrSB[0,15]
   - 8 uFrMC[0,0] - 8 uFrMC[0,1] - 8 uFrMC[0,2] - 8 uFrMC[0,3]
   - 8 uFrMC[0,4] - 8 uFrMC[0,5] - 8 uFrMC[0,6] - 8 uFrMC[0,7]
   - 8 uFrMC[0,8] - 8 uFrMC[0,9] - 8 uFrMC[0,10] - 8 uFrMC[0,11]
   - 8 uFrMC[0,12] - 8 uFrMC[0,13] - 8 uFrMC[0,14] - 8 uFrMC[0,15] + cb'
   = 0
 R6057: - 8 leqdz[9,0] - 8 leqdz[9,1] - 8 leqdz[9,2] - 8 leqdz[9,3]
   - 8 leqdz[9,4] - 8 leqdz[9,5] - 8 leqdz[9,6] - 8 leqdz[9,7]
   - 8 leqdz[9,8] - 8 leqdz[9,9] - 8 leqdz[9,10] - 8 leqdz[9,11]
   - 8 leqdz[9,12] - 8 leqdz[9,13] - 8 leqdz[9,14] - 8 leqdz[9,15] + rf
   = 0
 R6058: 8 ldx[8,0] + 8 ldx[8,1] + 8 ldx[8,2] + 8 ldx[8,3] + 8 ldx[8,4]
   + 8 ldx[8,5] + 8 ldx[8,6] + 8 ldx[8,7] + 8 ldx[8,8] + 8 ldx[8,9]
   + 8 ldx[8,10] + 8 ldx[8,11] + 8 ldx[8,12] + 8 ldx[8,13] + 8 ldx[8,14]
   + 8 ldx[8,15] - rf + cf = 0
 R6059: - 8 leqdz[8,0] - 8 leqdz[8,1] - 8 leqdz[8,2] - 8 leqdz[8,3]
   - 8 leqdz[8,4] - 8 leqdz[8,5] - 8 leqdz[8,6] - 8 leqdz[8,7]
   - 8 leqdz[8,8] - 8 leqdz[8,9] - 8 leqdz[8,10] - 8 leqdz[8,11]
   - 8 leqdz[8,12] - 8 leqdz[8,13] - 8 leqdz[8,14] - 8 leqdz[8,15]
   - 8 leqdz[9,0] - 8 leqdz[9,1] - 8 leqdz[9,2] - 8 leqdz[9,3]
   - 8 leqdz[9,4] - 8 leqdz[9,5] - 8 leqdz[9,6] - 8 leqdz[9,7]
   - 8 leqdz[9,8] - 8 leqdz[9,9] - 8 leqdz[9,10] - 8 leqdz[9,11]
   - 8 leqdz[9,12] - 8 leqdz[9,13] - 8 leqdz[9,14] - 8 leqdz[9,15] + mf
   = 0
 R6060: - 8 lGstk[9,0] - 8 lGstk[9,1] - 8 lGstk[9,2] - 8 lGstk[9,3]
   - 8 lGstk[9,4] - 8 lGstk[9,5] - 8 lGstk[9,6] - 8 lGstk[9,7]
   - 8 lGstk[9,8] - 8 lGstk[9,9] - 8 lGstk[9,10] - 8 lGstk[9,11]
   - 8 lGstk[9,12] - 8 lGstk[9,13] - 8 lGstk[9,14] - 8 lGstk[9,15]
   - 8 lGstk[10,0] - 8 lGstk[10,1] - 8 lGstk[10,2] - 8 lGstk[10,3]
   - 8 lGstk[10,4] - 8 lGstk[10,5] - 8 lGstk[10,6] - 8 lGstk[10,7]
   - 8 lGstk[10,8] - 8 lGstk[10,9] - 8 lGstk[10,10] - 8 lGstk[10,11]
   - 8 lGstk[10,12] - 8 lGstk[10,13] - 8 lGstk[10,14] - 8 lGstk[10,15]
   + mf' = 0
 R6061: - 8 lFrSB[8,0] - 8 lFrSB[8,1] - 8 lFrSB[8,2] - 8 lFrSB[8,3]
   - 8 lFrSB[8,4] - 8 lFrSB[8,5] - 8 lFrSB[8,6] - 8 lFrSB[8,7]
   - 8 lFrSB[8,8] - 8 lFrSB[8,9] - 8 lFrSB[8,10] - 8 lFrSB[8,11]
   - 8 lFrSB[8,12] - 8 lFrSB[8,13] - 8 lFrSB[8,14] - 8 lFrSB[8,15]
   - 8 lFrSB[9,0] - 8 lFrSB[9,1] - 8 lFrSB[9,2] - 8 lFrSB[9,3]
   - 8 lFrSB[9,4] - 8 lFrSB[9,5] - 8 lFrSB[9,6] - 8 lFrSB[9,7]
   - 8 lFrSB[9,8] - 8 lFrSB[9,9] - 8 lFrSB[9,10] - 8 lFrSB[9,11]
   - 8 lFrSB[9,12] - 8 lFrSB[9,13] - 8 lFrSB[9,14] - 8 lFrSB[9,15]
   - 8 lFrMC[8,0] - 8 lFrMC[8,1] - 8 lFrMC[8,2] - 8 lFrMC[8,3]
   - 8 lFrMC[8,4] - 8 lFrMC[8,5] - 8 lFrMC[8,6] - 8 lFrMC[8,7]
   - 8 lFrMC[8,8] - 8 lFrMC[8,9] - 8 lFrMC[8,10] - 8 lFrMC[8,11]
   - 8 lFrMC[8,12] - 8 lFrMC[8,13] - 8 lFrMC[8,14] - 8 lFrMC[8,15] + cf'
   = 0
 R6062: Dc = 131.1
 R6063: - 2 cb - 2 cf + Qc = 6.2
 R6064: - Dc + T0 = 2
 R6065: - mb' - mf' - Dc + T1 = 2
 R6066: - rb - mb' - mf' + cb' - Dc + T2 = 1
 R6067: - mb' - mf' + 2 cb' + 2 cf' - Qc + T31 = 0
 R6068: - mb - mf + T32 - epsilon = 6.2
 R6069: T4 = 152
 R6070: - T0 + Tc >= 0
 R6071: - T1 + Tc >= 0
 R6072: - T2 + Tc >= 0
 R6073: - T31 + Tc >= 0
 R6074: - T32 + Tc >= 0
 R6075: - T4 + Tc >= 0
 qc0: [ uDX[1,0] * lDY[1,0] - uDX[1,0] * ldy[1,0] + uDX[1,1] * lDY[1,1]
   - uDX[1,1] * ldy[1,1] + uDX[1,2] * lDY[1,2] - uDX[1,2] * ldy[1,2]
   + uDX[1,3] * lDY[1,3] - uDX[1,3] * ldy[1,3] + uDX[1,4] * lDY[1,4]
   - uDX[1,4] * ldy[1,4] + uDX[1,5] * lDY[1,5] - uDX[1,5] * ldy[1,5]
   + uDX[1,6] * lDY[1,6] - uDX[1,6] * ldy[1,6] + uDX[1,7] * lDY[1,7]
   - uDX[1,7] * ldy[1,7] + uDX[1,8] * lDY[1,8] - uDX[1,8] * ldy[1,8]
   + uDX[1,9] * lDY[1,9] - uDX[1,9] * ldy[1,9] + uDX[1,10] * lDY[1,10]
   - uDX[1,10] * ldy[1,10] + uDX[1,11] * lDY[1,11] - uDX[1,11] * ldy[1,11]
   + uDX[1,12] * lDY[1,12] - uDX[1,12] * ldy[1,12] + uDX[1,13] * lDY[1,13]
   - uDX[1,13] * ldy[1,13] + uDX[1,14] * lDY[1,14] - uDX[1,14] * ldy[1,14]
   + uDX[1,15] * lDY[1,15] - uDX[1,15] * ldy[1,15] + uDX[2,0] * lDY[2,0]
   - uDX[2,0] * ldy[2,0] + uDX[2,1] * lDY[2,1] - uDX[2,1] * ldy[2,1]
   + uDX[2,2] * lDY[2,2] - uDX[2,2] * ldy[2,2] + uDX[2,3] * lDY[2,3]
   - uDX[2,3] * ldy[2,3] + uDX[2,4] * lDY[2,4] - uDX[2,4] * ldy[2,4]
   + uDX[2,5] * lDY[2,5] - uDX[2,5] * ldy[2,5] + uDX[2,6] * lDY[2,6]
   - uDX[2,6] * ldy[2,6] + uDX[2,7] * lDY[2,7] - uDX[2,7] * ldy[2,7]
   + uDX[2,8] * lDY[2,8] - uDX[2,8] * ldy[2,8] + uDX[2,9] * lDY[2,9]
   - uDX[2,9] * ldy[2,9] + uDX[2,10] * lDY[2,10] - uDX[2,10] * ldy[2,10]
   + uDX[2,11] * lDY[2,11] - uDX[2,11] * ldy[2,11] + uDX[2,12] * lDY[2,12]
   - uDX[2,12] * ldy[2,12] + uDX[2,13] * lDY[2,13] - uDX[2,13] * ldy[2,13]
   + uDX[2,14] * lDY[2,14] - uDX[2,14] * ldy[2,14] + uDX[2,15] * lDY[2,15]
   - uDX[2,15] * ldy[2,15] + uDX[3,0] * lDY[3,0] - uDX[3,0] * ldy[3,0]
   + uDX[3,1] * lDY[3,1] - uDX[3,1] * ldy[3,1] + uDX[3,2] * lDY[3,2]
   - uDX[3,2] * ldy[3,2] + uDX[3,3] * lDY[3,3] - uDX[3,3] * ldy[3,3]
   + uDX[3,4] * lDY[3,4] - uDX[3,4] * ldy[3,4] + uDX[3,5] * lDY[3,5]
   - uDX[3,5] * ldy[3,5] + uDX[3,6] * lDY[3,6] - uDX[3,6] * ldy[3,6]
   + uDX[3,7] * lDY[3,7] - uDX[3,7] * ldy[3,7] + uDX[3,8] * lDY[3,8]
   - uDX[3,8] * ldy[3,8] + uDX[3,9] * lDY[3,9] - uDX[3,9] * ldy[3,9]
   + uDX[3,10] * lDY[3,10] - uDX[3,10] * ldy[3,10] + uDX[3,11] * lDY[3,11]
   - uDX[3,11] * ldy[3,11] + uDX[3,12] * lDY[3,12] - uDX[3,12] * ldy[3,12]
   + uDX[3,13] * lDY[3,13] - uDX[3,13] * ldy[3,13] + uDX[3,14] * lDY[3,14]
   - uDX[3,14] * ldy[3,14] + uDX[3,15] * lDY[3,15] - uDX[3,15] * ldy[3,15]
   + uDX[4,0] * lDY[4,0] - uDX[4,0] * ldy[4,0] + uDX[4,1] * lDY[4,1]
   - uDX[4,1] * ldy[4,1] + uDX[4,2] * lDY[4,2] - uDX[4,2] * ldy[4,2]
   + uDX[4,3] * lDY[4,3] - uDX[4,3] * ldy[4,3] + uDX[4,4] * lDY[4,4]
   - uDX[4,4] * ldy[4,4] + uDX[4,5] * lDY[4,5] - uDX[4,5] * ldy[4,5]
   + uDX[4,6] * lDY[4,6] - uDX[4,6] * ldy[4,6] + uDX[4,7] * lDY[4,7]
   - uDX[4,7] * ldy[4,7] + uDX[4,8] * lDY[4,8] - uDX[4,8] * ldy[4,8]
   + uDX[4,9] * lDY[4,9] - uDX[4,9] * ldy[4,9] + uDX[4,10] * lDY[4,10]
   - uDX[4,10] * ldy[4,10] + uDX[4,11] * lDY[4,11] - uDX[4,11] * ldy[4,11]
   + uDX[4,12] * lDY[4,12] - uDX[4,12] * ldy[4,12] + uDX[4,13] * lDY[4,13]
   - uDX[4,13] * ldy[4,13] + uDX[4,14] * lDY[4,14] - uDX[4,14] * ldy[4,14]
   + uDX[4,15] * lDY[4,15] - uDX[4,15] * ldy[4,15] + uDX[5,0] * lDY[5,0]
   - uDX[5,0] * ldy[5,0] + uDX[5,1] * lDY[5,1] - uDX[5,1] * ldy[5,1]
   + uDX[5,2] * lDY[5,2] - uDX[5,2] * ldy[5,2] + uDX[5,3] * lDY[5,3]
   - uDX[5,3] * ldy[5,3] + uDX[5,4] * lDY[5,4] - uDX[5,4] * ldy[5,4]
   + uDX[5,5] * lDY[5,5] - uDX[5,5] * ldy[5,5] + uDX[5,6] * lDY[5,6]
   - uDX[5,6] * ldy[5,6] + uDX[5,7] * lDY[5,7] - uDX[5,7] * ldy[5,7]
   + uDX[5,8] * lDY[5,8] - uDX[5,8] * ldy[5,8] + uDX[5,9] * lDY[5,9]
   - uDX[5,9] * ldy[5,9] + uDX[5,10] * lDY[5,10] - uDX[5,10] * ldy[5,10]
   + uDX[5,11] * lDY[5,11] - uDX[5,11] * ldy[5,11] + uDX[5,12] * lDY[5,12]
   - uDX[5,12] * ldy[5,12] + uDX[5,13] * lDY[5,13] - uDX[5,13] * ldy[5,13]
   + uDX[5,14] * lDY[5,14] - uDX[5,14] * ldy[5,14] + uDX[5,15] * lDY[5,15]
   - uDX[5,15] * ldy[5,15] + uDX[6,0] * lDY[6,0] - uDX[6,0] * ldy[6,0]
   + uDX[6,1] * lDY[6,1] - uDX[6,1] * ldy[6,1] + uDX[6,2] * lDY[6,2]
   - uDX[6,2] * ldy[6,2] + uDX[6,3] * lDY[6,3] - uDX[6,3] * ldy[6,3]
   + uDX[6,4] * lDY[6,4] - uDX[6,4] * ldy[6,4] + uDX[6,5] * lDY[6,5]
   - uDX[6,5] * ldy[6,5] + uDX[6,6] * lDY[6,6] - uDX[6,6] * ldy[6,6]
   + uDX[6,7] * lDY[6,7] - uDX[6,7] * ldy[6,7] + uDX[6,8] * lDY[6,8]
   - uDX[6,8] * ldy[6,8] + uDX[6,9] * lDY[6,9] - uDX[6,9] * ldy[6,9]
   + uDX[6,10] * lDY[6,10] - uDX[6,10] * ldy[6,10] + uDX[6,11] * lDY[6,11]
   - uDX[6,11] * ldy[6,11] + uDX[6,12] * lDY[6,12] - uDX[6,12] * ldy[6,12]
   + uDX[6,13] * lDY[6,13] - uDX[6,13] * ldy[6,13] + uDX[6,14] * lDY[6,14]
   - uDX[6,14] * ldy[6,14] + uDX[6,15] * lDY[6,15] - uDX[6,15] * ldy[6,15]
   + uDX[7,0] * lDY[7,0] - uDX[7,0] * ldy[7,0] + uDX[7,1] * lDY[7,1]
   - uDX[7,1] * ldy[7,1] + uDX[7,2] * lDY[7,2] - uDX[7,2] * ldy[7,2]
   + uDX[7,3] * lDY[7,3] - uDX[7,3] * ldy[7,3] + uDX[7,4] * lDY[7,4]
   - uDX[7,4] * ldy[7,4] + uDX[7,5] * lDY[7,5] - uDX[7,5] * ldy[7,5]
   + uDX[7,6] * lDY[7,6] - uDX[7,6] * ldy[7,6] + uDX[7,7] * lDY[7,7]
   - uDX[7,7] * ldy[7,7] + uDX[7,8] * lDY[7,8] - uDX[7,8] * ldy[7,8]
   + uDX[7,9] * lDY[7,9] - uDX[7,9] * ldy[7,9] + uDX[7,10] * lDY[7,10]
   - uDX[7,10] * ldy[7,10] + uDX[7,11] * lDY[7,11] - uDX[7,11] * ldy[7,11]
   + uDX[7,12] * lDY[7,12] - uDX[7,12] * ldy[7,12] + uDX[7,13] * lDY[7,13]
   - uDX[7,13] * ldy[7,13] + uDX[7,14] * lDY[7,14] - uDX[7,14] * ldy[7,14]
   + uDX[7,15] * lDY[7,15] - uDX[7,15] * ldy[7,15] - udx[1,0] * lDY[1,0]
   + udx[1,0] * ldy[1,0] - udx[1,1] * lDY[1,1] + udx[1,1] * ldy[1,1]
   - udx[1,2] * lDY[1,2] + udx[1,2] * ldy[1,2] - udx[1,3] * lDY[1,3]
   + udx[1,3] * ldy[1,3] - udx[1,4] * lDY[1,4] + udx[1,4] * ldy[1,4]
   - udx[1,5] * lDY[1,5] + udx[1,5] * ldy[1,5] - udx[1,6] * lDY[1,6]
   + udx[1,6] * ldy[1,6] - udx[1,7] * lDY[1,7] + udx[1,7] * ldy[1,7]
   - udx[1,8] * lDY[1,8] + udx[1,8] * ldy[1,8] - udx[1,9] * lDY[1,9]
   + udx[1,9] * ldy[1,9] - udx[1,10] * lDY[1,10] + udx[1,10] * ldy[1,10]
   - udx[1,11] * lDY[1,11] + udx[1,11] * ldy[1,11] - udx[1,12] * lDY[1,12]
   + udx[1,12] * ldy[1,12] - udx[1,13] * lDY[1,13] + udx[1,13] * ldy[1,13]
   - udx[1,14] * lDY[1,14] + udx[1,14] * ldy[1,14] - udx[1,15] * lDY[1,15]
   + udx[1,15] * ldy[1,15] - udx[2,0] * lDY[2,0] + udx[2,0] * ldy[2,0]
   - udx[2,1] * lDY[2,1] + udx[2,1] * ldy[2,1] - udx[2,2] * lDY[2,2]
   + udx[2,2] * ldy[2,2] - udx[2,3] * lDY[2,3] + udx[2,3] * ldy[2,3]
   - udx[2,4] * lDY[2,4] + udx[2,4] * ldy[2,4] - udx[2,5] * lDY[2,5]
   + udx[2,5] * ldy[2,5] - udx[2,6] * lDY[2,6] + udx[2,6] * ldy[2,6]
   - udx[2,7] * lDY[2,7] + udx[2,7] * ldy[2,7] - udx[2,8] * lDY[2,8]
   + udx[2,8] * ldy[2,8] - udx[2,9] * lDY[2,9] + udx[2,9] * ldy[2,9]
   - udx[2,10] * lDY[2,10] + udx[2,10] * ldy[2,10] - udx[2,11] * lDY[2,11]
   + udx[2,11] * ldy[2,11] - udx[2,12] * lDY[2,12] + udx[2,12] * ldy[2,12]
   - udx[2,13] * lDY[2,13] + udx[2,13] * ldy[2,13] - udx[2,14] * lDY[2,14]
   + udx[2,14] * ldy[2,14] - udx[2,15] * lDY[2,15] + udx[2,15] * ldy[2,15]
   - udx[3,0] * lDY[3,0] + udx[3,0] * ldy[3,0] - udx[3,1] * lDY[3,1]
   + udx[3,1] * ldy[3,1] - udx[3,2] * lDY[3,2] + udx[3,2] * ldy[3,2]
   - udx[3,3] * lDY[3,3] + udx[3,3] * ldy[3,3] - udx[3,4] * lDY[3,4]
   + udx[3,4] * ldy[3,4] - udx[3,5] * lDY[3,5] + udx[3,5] * ldy[3,5]
   - udx[3,6] * lDY[3,6] + udx[3,6] * ldy[3,6] - udx[3,7] * lDY[3,7]
   + udx[3,7] * ldy[3,7] - udx[3,8] * lDY[3,8] + udx[3,8] * ldy[3,8]
   - udx[3,9] * lDY[3,9] + udx[3,9] * ldy[3,9] - udx[3,10] * lDY[3,10]
   + udx[3,10] * ldy[3,10] - udx[3,11] * lDY[3,11] + udx[3,11] * ldy[3,11]
   - udx[3,12] * lDY[3,12] + udx[3,12] * ldy[3,12] - udx[3,13] * lDY[3,13]
   + udx[3,13] * ldy[3,13] - udx[3,14] * lDY[3,14] + udx[3,14] * ldy[3,14]
   - udx[3,15] * lDY[3,15] + udx[3,15] * ldy[3,15] - udx[4,0] * lDY[4,0]
   + udx[4,0] * ldy[4,0] - udx[4,1] * lDY[4,1] + udx[4,1] * ldy[4,1]
   - udx[4,2] * lDY[4,2] + udx[4,2] * ldy[4,2] - udx[4,3] * lDY[4,3]
   + udx[4,3] * ldy[4,3] - udx[4,4] * lDY[4,4] + udx[4,4] * ldy[4,4]
   - udx[4,5] * lDY[4,5] + udx[4,5] * ldy[4,5] - udx[4,6] * lDY[4,6]
   + udx[4,6] * ldy[4,6] - udx[4,7] * lDY[4,7] + udx[4,7] * ldy[4,7]
   - udx[4,8] * lDY[4,8] + udx[4,8] * ldy[4,8] - udx[4,9] * lDY[4,9]
   + udx[4,9] * ldy[4,9] - udx[4,10] * lDY[4,10] + udx[4,10] * ldy[4,10]
   - udx[4,11] * lDY[4,11] + udx[4,11] * ldy[4,11] - udx[4,12] * lDY[4,12]
   + udx[4,12] * ldy[4,12] - udx[4,13] * lDY[4,13] + udx[4,13] * ldy[4,13]
   - udx[4,14] * lDY[4,14] + udx[4,14] * ldy[4,14] - udx[4,15] * lDY[4,15]
   + udx[4,15] * ldy[4,15] - udx[5,0] * lDY[5,0] + udx[5,0] * ldy[5,0]
   - udx[5,1] * lDY[5,1] + udx[5,1] * ldy[5,1] - udx[5,2] * lDY[5,2]
   + udx[5,2] * ldy[5,2] - udx[5,3] * lDY[5,3] + udx[5,3] * ldy[5,3]
   - udx[5,4] * lDY[5,4] + udx[5,4] * ldy[5,4] - udx[5,5] * lDY[5,5]
   + udx[5,5] * ldy[5,5] - udx[5,6] * lDY[5,6] + udx[5,6] * ldy[5,6]
   - udx[5,7] * lDY[5,7] + udx[5,7] * ldy[5,7] - udx[5,8] * lDY[5,8]
   + udx[5,8] * ldy[5,8] - udx[5,9] * lDY[5,9] + udx[5,9] * ldy[5,9]
   - udx[5,10] * lDY[5,10] + udx[5,10] * ldy[5,10] - udx[5,11] * lDY[5,11]
   + udx[5,11] * ldy[5,11] - udx[5,12] * lDY[5,12] + udx[5,12] * ldy[5,12]
   - udx[5,13] * lDY[5,13] + udx[5,13] * ldy[5,13] - udx[5,14] * lDY[5,14]
   + udx[5,14] * ldy[5,14] - udx[5,15] * lDY[5,15] + udx[5,15] * ldy[5,15]
   - udx[6,0] * lDY[6,0] + udx[6,0] * ldy[6,0] - udx[6,1] * lDY[6,1]
   + udx[6,1] * ldy[6,1] - udx[6,2] * lDY[6,2] + udx[6,2] * ldy[6,2]
   - udx[6,3] * lDY[6,3] + udx[6,3] * ldy[6,3] - udx[6,4] * lDY[6,4]
   + udx[6,4] * ldy[6,4] - udx[6,5] * lDY[6,5] + udx[6,5] * ldy[6,5]
   - udx[6,6] * lDY[6,6] + udx[6,6] * ldy[6,6] - udx[6,7] * lDY[6,7]
   + udx[6,7] * ldy[6,7] - udx[6,8] * lDY[6,8] + udx[6,8] * ldy[6,8]
   - udx[6,9] * lDY[6,9] + udx[6,9] * ldy[6,9] - udx[6,10] * lDY[6,10]
   + udx[6,10] * ldy[6,10] - udx[6,11] * lDY[6,11] + udx[6,11] * ldy[6,11]
   - udx[6,12] * lDY[6,12] + udx[6,12] * ldy[6,12] - udx[6,13] * lDY[6,13]
   + udx[6,13] * ldy[6,13] - udx[6,14] * lDY[6,14] + udx[6,14] * ldy[6,14]
   - udx[6,15] * lDY[6,15] + udx[6,15] * ldy[6,15] - udx[7,0] * lDY[7,0]
   + udx[7,0] * ldy[7,0] - udx[7,1] * lDY[7,1] + udx[7,1] * ldy[7,1]
   - udx[7,2] * lDY[7,2] + udx[7,2] * ldy[7,2] - udx[7,3] * lDY[7,3]
   + udx[7,3] * ldy[7,3] - udx[7,4] * lDY[7,4] + udx[7,4] * ldy[7,4]
   - udx[7,5] * lDY[7,5] + udx[7,5] * ldy[7,5] - udx[7,6] * lDY[7,6]
   + udx[7,6] * ldy[7,6] - udx[7,7] * lDY[7,7] + udx[7,7] * ldy[7,7]
   - udx[7,8] * lDY[7,8] + udx[7,8] * ldy[7,8] - udx[7,9] * lDY[7,9]
   + udx[7,9] * ldy[7,9] - udx[7,10] * lDY[7,10] + udx[7,10] * ldy[7,10]
   - udx[7,11] * lDY[7,11] + udx[7,11] * ldy[7,11] - udx[7,12] * lDY[7,12]
   + udx[7,12] * ldy[7,12] - udx[7,13] * lDY[7,13] + udx[7,13] * ldy[7,13]
   - udx[7,14] * lDY[7,14] + udx[7,14] * ldy[7,14] - udx[7,15] * lDY[7,15]
   + udx[7,15] * ldy[7,15] ] >= 1
Bounds
Binaries
 uDW[-1,0] uDW[-1,1] uDW[-1,2] uDW[-1,3] uDW[-1,4] uDW[-1,5] uDW[-1,6]
 uDW[-1,7] uDW[-1,8] uDW[-1,9] uDW[-1,10] uDW[-1,11] uDW[-1,12] uDW[-1,13]
 uDW[-1,14] uDW[-1,15] uDW[0,0] uDW[0,1] uDW[0,2] uDW[0,3] uDW[0,4]
 uDW[0,5] uDW[0,6] uDW[0,7] uDW[0,8] uDW[0,9] uDW[0,10] uDW[0,11] uDW[0,12]
 uDW[0,13] uDW[0,14] uDW[0,15] uDW[1,0] uDW[1,1] uDW[1,2] uDW[1,3] uDW[1,4]
 uDW[1,5] uDW[1,6] uDW[1,7] uDW[1,8] uDW[1,9] uDW[1,10] uDW[1,11] uDW[1,12]
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
 uDW[7,13] uDW[7,14] uDW[7,15] udw[-1,0] udw[-1,1] udw[-1,2] udw[-1,3]
 udw[-1,4] udw[-1,5] udw[-1,6] udw[-1,7] udw[-1,8] udw[-1,9] udw[-1,10]
 udw[-1,11] udw[-1,12] udw[-1,13] udw[-1,14] udw[-1,15] udw[0,0] udw[0,1]
 udw[0,2] udw[0,3] udw[0,4] udw[0,5] udw[0,6] udw[0,7] udw[0,8] udw[0,9]
 udw[0,10] udw[0,11] udw[0,12] udw[0,13] udw[0,14] udw[0,15] udw[1,0]
 udw[1,1] udw[1,2] udw[1,3] udw[1,4] udw[1,5] udw[1,6] udw[1,7] udw[1,8]
 udw[1,9] udw[1,10] udw[1,11] udw[1,12] udw[1,13] udw[1,14] udw[1,15]
 udw[2,0] udw[2,1] udw[2,2] udw[2,3] udw[2,4] udw[2,5] udw[2,6] udw[2,7]
 udw[2,8] udw[2,9] udw[2,10] udw[2,11] udw[2,12] udw[2,13] udw[2,14]
 udw[2,15] udw[3,0] udw[3,1] udw[3,2] udw[3,3] udw[3,4] udw[3,5] udw[3,6]
 udw[3,7] udw[3,8] udw[3,9] udw[3,10] udw[3,11] udw[3,12] udw[3,13]
 udw[3,14] udw[3,15] udw[4,0] udw[4,1] udw[4,2] udw[4,3] udw[4,4] udw[4,5]
 udw[4,6] udw[4,7] udw[4,8] udw[4,9] udw[4,10] udw[4,11] udw[4,12]
 udw[4,13] udw[4,14] udw[4,15] udw[5,0] udw[5,1] udw[5,2] udw[5,3] udw[5,4]
 udw[5,5] udw[5,6] udw[5,7] udw[5,8] udw[5,9] udw[5,10] udw[5,11] udw[5,12]
 udw[5,13] udw[5,14] udw[5,15] udw[6,0] udw[6,1] udw[6,2] udw[6,3] udw[6,4]
 udw[6,5] udw[6,6] udw[6,7] udw[6,8] udw[6,9] udw[6,10] udw[6,11] udw[6,12]
 udw[6,13] udw[6,14] udw[6,15] udw[7,0] udw[7,1] udw[7,2] udw[7,3] udw[7,4]
 udw[7,5] udw[7,6] udw[7,7] udw[7,8] udw[7,9] udw[7,10] udw[7,11] udw[7,12]
 udw[7,13] udw[7,14] udw[7,15] uDX[0,0] uDX[0,1] uDX[0,2] uDX[0,3] uDX[0,4]
 uDX[0,5] uDX[0,6] uDX[0,7] uDX[0,8] uDX[0,9] uDX[0,10] uDX[0,11] uDX[0,12]
 uDX[0,13] uDX[0,14] uDX[0,15] uDX[1,0] uDX[1,1] uDX[1,2] uDX[1,3] uDX[1,4]
 uDX[1,5] uDX[1,6] uDX[1,7] uDX[1,8] uDX[1,9] uDX[1,10] uDX[1,11] uDX[1,12]
 uDX[1,13] uDX[1,14] uDX[1,15] uDX[2,0] uDX[2,1] uDX[2,2] uDX[2,3] uDX[2,4]
 uDX[2,5] uDX[2,6] uDX[2,7] uDX[2,8] uDX[2,9] uDX[2,10] uDX[2,11] uDX[2,12]
 uDX[2,13] uDX[2,14] uDX[2,15] uDX[3,0] uDX[3,1] uDX[3,2] uDX[3,3] uDX[3,4]
 uDX[3,5] uDX[3,6] uDX[3,7] uDX[3,8] uDX[3,9] uDX[3,10] uDX[3,11] uDX[3,12]
 uDX[3,13] uDX[3,14] uDX[3,15] uDX[4,0] uDX[4,1] uDX[4,2] uDX[4,3] uDX[4,4]
 uDX[4,5] uDX[4,6] uDX[4,7] uDX[4,8] uDX[4,9] uDX[4,10] uDX[4,11] uDX[4,12]
 uDX[4,13] uDX[4,14] uDX[4,15] uDX[5,0] uDX[5,1] uDX[5,2] uDX[5,3] uDX[5,4]
 uDX[5,5] uDX[5,6] uDX[5,7] uDX[5,8] uDX[5,9] uDX[5,10] uDX[5,11] uDX[5,12]
 uDX[5,13] uDX[5,14] uDX[5,15] uDX[6,0] uDX[6,1] uDX[6,2] uDX[6,3] uDX[6,4]
 uDX[6,5] uDX[6,6] uDX[6,7] uDX[6,8] uDX[6,9] uDX[6,10] uDX[6,11] uDX[6,12]
 uDX[6,13] uDX[6,14] uDX[6,15] uDX[7,0] uDX[7,1] uDX[7,2] uDX[7,3] uDX[7,4]
 uDX[7,5] uDX[7,6] uDX[7,7] uDX[7,8] uDX[7,9] uDX[7,10] uDX[7,11] uDX[7,12]
 uDX[7,13] uDX[7,14] uDX[7,15] udx[0,0] udx[0,1] udx[0,2] udx[0,3] udx[0,4]
 udx[0,5] udx[0,6] udx[0,7] udx[0,8] udx[0,9] udx[0,10] udx[0,11] udx[0,12]
 udx[0,13] udx[0,14] udx[0,15] udx[1,0] udx[1,1] udx[1,2] udx[1,3] udx[1,4]
 udx[1,5] udx[1,6] udx[1,7] udx[1,8] udx[1,9] udx[1,10] udx[1,11] udx[1,12]
 udx[1,13] udx[1,14] udx[1,15] udx[2,0] udx[2,1] udx[2,2] udx[2,3] udx[2,4]
 udx[2,5] udx[2,6] udx[2,7] udx[2,8] udx[2,9] udx[2,10] udx[2,11] udx[2,12]
 udx[2,13] udx[2,14] udx[2,15] udx[3,0] udx[3,1] udx[3,2] udx[3,3] udx[3,4]
 udx[3,5] udx[3,6] udx[3,7] udx[3,8] udx[3,9] udx[3,10] udx[3,11] udx[3,12]
 udx[3,13] udx[3,14] udx[3,15] udx[4,0] udx[4,1] udx[4,2] udx[4,3] udx[4,4]
 udx[4,5] udx[4,6] udx[4,7] udx[4,8] udx[4,9] udx[4,10] udx[4,11] udx[4,12]
 udx[4,13] udx[4,14] udx[4,15] udx[5,0] udx[5,1] udx[5,2] udx[5,3] udx[5,4]
 udx[5,5] udx[5,6] udx[5,7] udx[5,8] udx[5,9] udx[5,10] udx[5,11] udx[5,12]
 udx[5,13] udx[5,14] udx[5,15] udx[6,0] udx[6,1] udx[6,2] udx[6,3] udx[6,4]
 udx[6,5] udx[6,6] udx[6,7] udx[6,8] udx[6,9] udx[6,10] udx[6,11] udx[6,12]
 udx[6,13] udx[6,14] udx[6,15] udx[7,0] udx[7,1] udx[7,2] udx[7,3] udx[7,4]
 udx[7,5] udx[7,6] udx[7,7] udx[7,8] udx[7,9] udx[7,10] udx[7,11] udx[7,12]
 udx[7,13] udx[7,14] udx[7,15] uDY[0,0] uDY[0,1] uDY[0,2] uDY[0,3] uDY[0,4]
 uDY[0,5] uDY[0,6] uDY[0,7] uDY[0,8] uDY[0,9] uDY[0,10] uDY[0,11] uDY[0,12]
 uDY[0,13] uDY[0,14] uDY[0,15] uDY[1,0] uDY[1,1] uDY[1,2] uDY[1,3] uDY[1,4]
 uDY[1,5] uDY[1,6] uDY[1,7] uDY[1,8] uDY[1,9] uDY[1,10] uDY[1,11] uDY[1,12]
 uDY[1,13] uDY[1,14] uDY[1,15] uDY[2,0] uDY[2,1] uDY[2,2] uDY[2,3] uDY[2,4]
 uDY[2,5] uDY[2,6] uDY[2,7] uDY[2,8] uDY[2,9] uDY[2,10] uDY[2,11] uDY[2,12]
 uDY[2,13] uDY[2,14] uDY[2,15] uDY[3,0] uDY[3,1] uDY[3,2] uDY[3,3] uDY[3,4]
 uDY[3,5] uDY[3,6] uDY[3,7] uDY[3,8] uDY[3,9] uDY[3,10] uDY[3,11] uDY[3,12]
 uDY[3,13] uDY[3,14] uDY[3,15] uDY[4,0] uDY[4,1] uDY[4,2] uDY[4,3] uDY[4,4]
 uDY[4,5] uDY[4,6] uDY[4,7] uDY[4,8] uDY[4,9] uDY[4,10] uDY[4,11] uDY[4,12]
 uDY[4,13] uDY[4,14] uDY[4,15] uDY[5,0] uDY[5,1] uDY[5,2] uDY[5,3] uDY[5,4]
 uDY[5,5] uDY[5,6] uDY[5,7] uDY[5,8] uDY[5,9] uDY[5,10] uDY[5,11] uDY[5,12]
 uDY[5,13] uDY[5,14] uDY[5,15] uDY[6,0] uDY[6,1] uDY[6,2] uDY[6,3] uDY[6,4]
 uDY[6,5] uDY[6,6] uDY[6,7] uDY[6,8] uDY[6,9] uDY[6,10] uDY[6,11] uDY[6,12]
 uDY[6,13] uDY[6,14] uDY[6,15] uDY[7,0] uDY[7,1] uDY[7,2] uDY[7,3] uDY[7,4]
 uDY[7,5] uDY[7,6] uDY[7,7] uDY[7,8] uDY[7,9] uDY[7,10] uDY[7,11] uDY[7,12]
 uDY[7,13] uDY[7,14] uDY[7,15] udy[0,0] udy[0,1] udy[0,2] udy[0,3] udy[0,4]
 udy[0,5] udy[0,6] udy[0,7] udy[0,8] udy[0,9] udy[0,10] udy[0,11] udy[0,12]
 udy[0,13] udy[0,14] udy[0,15] udy[1,0] udy[1,1] udy[1,2] udy[1,3] udy[1,4]
 udy[1,5] udy[1,6] udy[1,7] udy[1,8] udy[1,9] udy[1,10] udy[1,11] udy[1,12]
 udy[1,13] udy[1,14] udy[1,15] udy[2,0] udy[2,1] udy[2,2] udy[2,3] udy[2,4]
 udy[2,5] udy[2,6] udy[2,7] udy[2,8] udy[2,9] udy[2,10] udy[2,11] udy[2,12]
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
 udy[7,13] udy[7,14] udy[7,15] uDZ[0,0] uDZ[0,1] uDZ[0,2] uDZ[0,3] uDZ[0,4]
 uDZ[0,5] uDZ[0,6] uDZ[0,7] uDZ[0,8] uDZ[0,9] uDZ[0,10] uDZ[0,11] uDZ[0,12]
 uDZ[0,13] uDZ[0,14] uDZ[0,15] uDZ[1,0] uDZ[1,1] uDZ[1,2] uDZ[1,3] uDZ[1,4]
 uDZ[1,5] uDZ[1,6] uDZ[1,7] uDZ[1,8] uDZ[1,9] uDZ[1,10] uDZ[1,11] uDZ[1,12]
 uDZ[1,13] uDZ[1,14] uDZ[1,15] uDZ[2,0] uDZ[2,1] uDZ[2,2] uDZ[2,3] uDZ[2,4]
 uDZ[2,5] uDZ[2,6] uDZ[2,7] uDZ[2,8] uDZ[2,9] uDZ[2,10] uDZ[2,11] uDZ[2,12]
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
 uDZ[7,13] uDZ[7,14] uDZ[7,15] udz[0,0] udz[0,1] udz[0,2] udz[0,3] udz[0,4]
 udz[0,5] udz[0,6] udz[0,7] udz[0,8] udz[0,9] udz[0,10] udz[0,11] udz[0,12]
 udz[0,13] udz[0,14] udz[0,15] udz[1,0] udz[1,1] udz[1,2] udz[1,3] udz[1,4]
 udz[1,5] udz[1,6] udz[1,7] udz[1,8] udz[1,9] udz[1,10] udz[1,11] udz[1,12]
 udz[1,13] udz[1,14] udz[1,15] udz[2,0] udz[2,1] udz[2,2] udz[2,3] udz[2,4]
 udz[2,5] udz[2,6] udz[2,7] udz[2,8] udz[2,9] udz[2,10] udz[2,11] udz[2,12]
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
 udz[7,13] udz[7,14] udz[7,15] ustk[0,0] ustk[0,1] ustk[0,2] ustk[0,3]
 ustk[0,4] ustk[0,5] ustk[0,6] ustk[0,7] ustk[0,8] ustk[0,9] ustk[0,10]
 ustk[0,11] ustk[0,12] ustk[0,13] ustk[0,14] ustk[0,15] ustk[1,0] ustk[1,1]
 ustk[1,2] ustk[1,3] ustk[1,4] ustk[1,5] ustk[1,6] ustk[1,7] ustk[1,8]
 ustk[1,9] ustk[1,10] ustk[1,11] ustk[1,12] ustk[1,13] ustk[1,14]
 ustk[1,15] ustk[2,0] ustk[2,1] ustk[2,2] ustk[2,3] ustk[2,4] ustk[2,5]
 ustk[2,6] ustk[2,7] ustk[2,8] ustk[2,9] ustk[2,10] ustk[2,11] ustk[2,12]
 ustk[2,13] ustk[2,14] ustk[2,15] ustk[3,0] ustk[3,1] ustk[3,2] ustk[3,3]
 ustk[3,4] ustk[3,5] ustk[3,6] ustk[3,7] ustk[3,8] ustk[3,9] ustk[3,10]
 ustk[3,11] ustk[3,12] ustk[3,13] ustk[3,14] ustk[3,15] ustk[4,0] ustk[4,1]
 ustk[4,2] ustk[4,3] ustk[4,4] ustk[4,5] ustk[4,6] ustk[4,7] ustk[4,8]
 ustk[4,9] ustk[4,10] ustk[4,11] ustk[4,12] ustk[4,13] ustk[4,14]
 ustk[4,15] ustk[5,0] ustk[5,1] ustk[5,2] ustk[5,3] ustk[5,4] ustk[5,5]
 ustk[5,6] ustk[5,7] ustk[5,8] ustk[5,9] ustk[5,10] ustk[5,11] ustk[5,12]
 ustk[5,13] ustk[5,14] ustk[5,15] ustk[6,0] ustk[6,1] ustk[6,2] ustk[6,3]
 ustk[6,4] ustk[6,5] ustk[6,6] ustk[6,7] ustk[6,8] ustk[6,9] ustk[6,10]
 ustk[6,11] ustk[6,12] ustk[6,13] ustk[6,14] ustk[6,15] ustk[7,0] ustk[7,1]
 ustk[7,2] ustk[7,3] ustk[7,4] ustk[7,5] ustk[7,6] ustk[7,7] ustk[7,8]
 ustk[7,9] ustk[7,10] ustk[7,11] ustk[7,12] ustk[7,13] ustk[7,14]
 ustk[7,15] ustk[8,0] ustk[8,1] ustk[8,2] ustk[8,3] ustk[8,4] ustk[8,5]
 ustk[8,6] ustk[8,7] ustk[8,8] ustk[8,9] ustk[8,10] ustk[8,11] ustk[8,12]
 ustk[8,13] ustk[8,14] ustk[8,15] ucan[0,0] ucan[0,1] ucan[0,2] ucan[0,3]
 ucan[0,4] ucan[0,5] ucan[0,6] ucan[0,7] ucan[0,8] ucan[0,9] ucan[0,10]
 ucan[0,11] ucan[0,12] ucan[0,13] ucan[0,14] ucan[0,15] ucan[1,0] ucan[1,1]
 ucan[1,2] ucan[1,3] ucan[1,4] ucan[1,5] ucan[1,6] ucan[1,7] ucan[1,8]
 ucan[1,9] ucan[1,10] ucan[1,11] ucan[1,12] ucan[1,13] ucan[1,14]
 ucan[1,15] ucan[2,0] ucan[2,1] ucan[2,2] ucan[2,3] ucan[2,4] ucan[2,5]
 ucan[2,6] ucan[2,7] ucan[2,8] ucan[2,9] ucan[2,10] ucan[2,11] ucan[2,12]
 ucan[2,13] ucan[2,14] ucan[2,15] ucan[3,0] ucan[3,1] ucan[3,2] ucan[3,3]
 ucan[3,4] ucan[3,5] ucan[3,6] ucan[3,7] ucan[3,8] ucan[3,9] ucan[3,10]
 ucan[3,11] ucan[3,12] ucan[3,13] ucan[3,14] ucan[3,15] ucan[4,0] ucan[4,1]
 ucan[4,2] ucan[4,3] ucan[4,4] ucan[4,5] ucan[4,6] ucan[4,7] ucan[4,8]
 ucan[4,9] ucan[4,10] ucan[4,11] ucan[4,12] ucan[4,13] ucan[4,14]
 ucan[4,15] ucan[5,0] ucan[5,1] ucan[5,2] ucan[5,3] ucan[5,4] ucan[5,5]
 ucan[5,6] ucan[5,7] ucan[5,8] ucan[5,9] ucan[5,10] ucan[5,11] ucan[5,12]
 ucan[5,13] ucan[5,14] ucan[5,15] ucan[6,0] ucan[6,1] ucan[6,2] ucan[6,3]
 ucan[6,4] ucan[6,5] ucan[6,6] ucan[6,7] ucan[6,8] ucan[6,9] ucan[6,10]
 ucan[6,11] ucan[6,12] ucan[6,13] ucan[6,14] ucan[6,15] ucan[7,0] ucan[7,1]
 ucan[7,2] ucan[7,3] ucan[7,4] ucan[7,5] ucan[7,6] ucan[7,7] ucan[7,8]
 ucan[7,9] ucan[7,10] ucan[7,11] ucan[7,12] ucan[7,13] ucan[7,14]
 ucan[7,15] uAC[0,0] uAC[0,1] uAC[0,2] uAC[0,3] uAC[1,0] uAC[1,1] uAC[1,2]
 uAC[1,3] uAC[2,0] uAC[2,1] uAC[2,2] uAC[2,3] uAC[3,0] uAC[3,1] uAC[3,2]
 uAC[3,3] uAC[4,0] uAC[4,1] uAC[4,2] uAC[4,3] uAC[5,0] uAC[5,1] uAC[5,2]
 uAC[5,3] uAC[6,0] uAC[6,1] uAC[6,2] uAC[6,3] uAC[7,0] uAC[7,1] uAC[7,2]
 uAC[7,3] uTC[0,0] uTC[0,1] uTC[0,2] uTC[0,3] uTC[1,0] uTC[1,1] uTC[1,2]
 uTC[1,3] uTC[2,0] uTC[2,1] uTC[2,2] uTC[2,3] uTC[3,0] uTC[3,1] uTC[3,2]
 uTC[3,3] uTC[4,0] uTC[4,1] uTC[4,2] uTC[4,3] uTC[5,0] uTC[5,1] uTC[5,2]
 uTC[5,3] uTC[6,0] uTC[6,1] uTC[6,2] uTC[6,3] uTC[7,0] uTC[7,1] uTC[7,2]
 uTC[7,3] uLANE[0] uLANE[1] uLANE[2] uLANE[3] uLANE[4] uLANE[5] uLANE[6]
 uLANE[7] uLANE[8] uLANE[9] uLANE[10] uLANE[11] uLANE[12] uLANE[13]
 uLANE[14] uLANE[15] uGstk[0,0] uGstk[0,1] uGstk[0,2] uGstk[0,3] uGstk[0,4]
 uGstk[0,5] uGstk[0,6] uGstk[0,7] uGstk[0,8] uGstk[0,9] uGstk[0,10]
 uGstk[0,11] uGstk[0,12] uGstk[0,13] uGstk[0,14] uGstk[0,15] uDetX[0,0]
 uDetX[0,1] uDetX[0,2] uDetX[0,3] uDetX[0,4] uDetX[0,5] uDetX[0,6]
 uDetX[0,7] uDetX[0,8] uDetX[0,9] uDetX[0,10] uDetX[0,11] uDetX[0,12]
 uDetX[0,13] uDetX[0,14] uDetX[0,15] uFrSB[0,0] uFrSB[0,1] uFrSB[0,2]
 uFrSB[0,3] uFrSB[0,4] uFrSB[0,5] uFrSB[0,6] uFrSB[0,7] uFrSB[0,8]
 uFrSB[0,9] uFrSB[0,10] uFrSB[0,11] uFrSB[0,12] uFrSB[0,13] uFrSB[0,14]
 uFrSB[0,15] uDetY[0,0] uDetY[0,1] uDetY[0,2] uDetY[0,3] uDetY[0,4]
 uDetY[0,5] uDetY[0,6] uDetY[0,7] uDetY[0,8] uDetY[0,9] uDetY[0,10]
 uDetY[0,11] uDetY[0,12] uDetY[0,13] uDetY[0,14] uDetY[0,15] uDetZ[0,0]
 uDetZ[0,1] uDetZ[0,2] uDetZ[0,3] uDetZ[0,4] uDetZ[0,5] uDetZ[0,6]
 uDetZ[0,7] uDetZ[0,8] uDetZ[0,9] uDetZ[0,10] uDetZ[0,11] uDetZ[0,12]
 uDetZ[0,13] uDetZ[0,14] uDetZ[0,15] uDetW[-1,0] uDetW[-1,1] uDetW[-1,2]
 uDetW[-1,3] uDetW[-1,4] uDetW[-1,5] uDetW[-1,6] uDetW[-1,7] uDetW[-1,8]
 uDetW[-1,9] uDetW[-1,10] uDetW[-1,11] uDetW[-1,12] uDetW[-1,13]
 uDetW[-1,14] uDetW[-1,15] uDetW[0,0] uDetW[0,1] uDetW[0,2] uDetW[0,3]
 uDetW[0,4] uDetW[0,5] uDetW[0,6] uDetW[0,7] uDetW[0,8] uDetW[0,9]
 uDetW[0,10] uDetW[0,11] uDetW[0,12] uDetW[0,13] uDetW[0,14] uDetW[0,15]
 uDetC[0,0] uDetC[0,1] uDetC[0,2] uDetC[0,3] uDiffDetMC[0,0]
 uDiffDetMC[0,1] uDiffDetMC[0,2] uDiffDetMC[0,3] uDiffDetMC[0,4]
 uDiffDetMC[0,5] uDiffDetMC[0,6] uDiffDetMC[0,7] uDiffDetMC[0,8]
 uDiffDetMC[0,9] uDiffDetMC[0,10] uDiffDetMC[0,11] uDiffDetMC[0,12]
 uDiffDetMC[0,13] uDiffDetMC[0,14] uDiffDetMC[0,15] uDiffDetC[0,0]
 uDiffDetC[0,1] uDiffDetC[0,2] uDiffDetC[0,3] uFrMC[0,0] uFrMC[0,1]
 uFrMC[0,2] uFrMC[0,3] uFrMC[0,4] uFrMC[0,5] uFrMC[0,6] uFrMC[0,7]
 uFrMC[0,8] uFrMC[0,9] uFrMC[0,10] uFrMC[0,11] uFrMC[0,12] uFrMC[0,13]
 uFrMC[0,14] uFrMC[0,15] lDW[0,0] lDW[0,1] lDW[0,2] lDW[0,3] lDW[0,4]
 lDW[0,5] lDW[0,6] lDW[0,7] lDW[0,8] lDW[0,9] lDW[0,10] lDW[0,11] lDW[0,12]
 lDW[0,13] lDW[0,14] lDW[0,15] lDW[1,0] lDW[1,1] lDW[1,2] lDW[1,3] lDW[1,4]
 lDW[1,5] lDW[1,6] lDW[1,7] lDW[1,8] lDW[1,9] lDW[1,10] lDW[1,11] lDW[1,12]
 lDW[1,13] lDW[1,14] lDW[1,15] lDW[2,0] lDW[2,1] lDW[2,2] lDW[2,3] lDW[2,4]
 lDW[2,5] lDW[2,6] lDW[2,7] lDW[2,8] lDW[2,9] lDW[2,10] lDW[2,11] lDW[2,12]
 lDW[2,13] lDW[2,14] lDW[2,15] lDW[3,0] lDW[3,1] lDW[3,2] lDW[3,3] lDW[3,4]
 lDW[3,5] lDW[3,6] lDW[3,7] lDW[3,8] lDW[3,9] lDW[3,10] lDW[3,11] lDW[3,12]
 lDW[3,13] lDW[3,14] lDW[3,15] lDW[4,0] lDW[4,1] lDW[4,2] lDW[4,3] lDW[4,4]
 lDW[4,5] lDW[4,6] lDW[4,7] lDW[4,8] lDW[4,9] lDW[4,10] lDW[4,11] lDW[4,12]
 lDW[4,13] lDW[4,14] lDW[4,15] lDW[5,0] lDW[5,1] lDW[5,2] lDW[5,3] lDW[5,4]
 lDW[5,5] lDW[5,6] lDW[5,7] lDW[5,8] lDW[5,9] lDW[5,10] lDW[5,11] lDW[5,12]
 lDW[5,13] lDW[5,14] lDW[5,15] lDW[6,0] lDW[6,1] lDW[6,2] lDW[6,3] lDW[6,4]
 lDW[6,5] lDW[6,6] lDW[6,7] lDW[6,8] lDW[6,9] lDW[6,10] lDW[6,11] lDW[6,12]
 lDW[6,13] lDW[6,14] lDW[6,15] lDW[7,0] lDW[7,1] lDW[7,2] lDW[7,3] lDW[7,4]
 lDW[7,5] lDW[7,6] lDW[7,7] lDW[7,8] lDW[7,9] lDW[7,10] lDW[7,11] lDW[7,12]
 lDW[7,13] lDW[7,14] lDW[7,15] ldw[0,0] ldw[0,1] ldw[0,2] ldw[0,3] ldw[0,4]
 ldw[0,5] ldw[0,6] ldw[0,7] ldw[0,8] ldw[0,9] ldw[0,10] ldw[0,11] ldw[0,12]
 ldw[0,13] ldw[0,14] ldw[0,15] ldw[1,0] ldw[1,1] ldw[1,2] ldw[1,3] ldw[1,4]
 ldw[1,5] ldw[1,6] ldw[1,7] ldw[1,8] ldw[1,9] ldw[1,10] ldw[1,11] ldw[1,12]
 ldw[1,13] ldw[1,14] ldw[1,15] ldw[2,0] ldw[2,1] ldw[2,2] ldw[2,3] ldw[2,4]
 ldw[2,5] ldw[2,6] ldw[2,7] ldw[2,8] ldw[2,9] ldw[2,10] ldw[2,11] ldw[2,12]
 ldw[2,13] ldw[2,14] ldw[2,15] ldw[3,0] ldw[3,1] ldw[3,2] ldw[3,3] ldw[3,4]
 ldw[3,5] ldw[3,6] ldw[3,7] ldw[3,8] ldw[3,9] ldw[3,10] ldw[3,11] ldw[3,12]
 ldw[3,13] ldw[3,14] ldw[3,15] ldw[4,0] ldw[4,1] ldw[4,2] ldw[4,3] ldw[4,4]
 ldw[4,5] ldw[4,6] ldw[4,7] ldw[4,8] ldw[4,9] ldw[4,10] ldw[4,11] ldw[4,12]
 ldw[4,13] ldw[4,14] ldw[4,15] ldw[5,0] ldw[5,1] ldw[5,2] ldw[5,3] ldw[5,4]
 ldw[5,5] ldw[5,6] ldw[5,7] ldw[5,8] ldw[5,9] ldw[5,10] ldw[5,11] ldw[5,12]
 ldw[5,13] ldw[5,14] ldw[5,15] ldw[6,0] ldw[6,1] ldw[6,2] ldw[6,3] ldw[6,4]
 ldw[6,5] ldw[6,6] ldw[6,7] ldw[6,8] ldw[6,9] ldw[6,10] ldw[6,11] ldw[6,12]
 ldw[6,13] ldw[6,14] ldw[6,15] ldw[7,0] ldw[7,1] ldw[7,2] ldw[7,3] ldw[7,4]
 ldw[7,5] ldw[7,6] ldw[7,7] ldw[7,8] ldw[7,9] ldw[7,10] ldw[7,11] ldw[7,12]
 ldw[7,13] ldw[7,14] ldw[7,15] lDX[1,0] lDX[1,1] lDX[1,2] lDX[1,3] lDX[1,4]
 lDX[1,5] lDX[1,6] lDX[1,7] lDX[1,8] lDX[1,9] lDX[1,10] lDX[1,11] lDX[1,12]
 lDX[1,13] lDX[1,14] lDX[1,15] lDX[2,0] lDX[2,1] lDX[2,2] lDX[2,3] lDX[2,4]
 lDX[2,5] lDX[2,6] lDX[2,7] lDX[2,8] lDX[2,9] lDX[2,10] lDX[2,11] lDX[2,12]
 lDX[2,13] lDX[2,14] lDX[2,15] lDX[3,0] lDX[3,1] lDX[3,2] lDX[3,3] lDX[3,4]
 lDX[3,5] lDX[3,6] lDX[3,7] lDX[3,8] lDX[3,9] lDX[3,10] lDX[3,11] lDX[3,12]
 lDX[3,13] lDX[3,14] lDX[3,15] lDX[4,0] lDX[4,1] lDX[4,2] lDX[4,3] lDX[4,4]
 lDX[4,5] lDX[4,6] lDX[4,7] lDX[4,8] lDX[4,9] lDX[4,10] lDX[4,11] lDX[4,12]
 lDX[4,13] lDX[4,14] lDX[4,15] lDX[5,0] lDX[5,1] lDX[5,2] lDX[5,3] lDX[5,4]
 lDX[5,5] lDX[5,6] lDX[5,7] lDX[5,8] lDX[5,9] lDX[5,10] lDX[5,11] lDX[5,12]
 lDX[5,13] lDX[5,14] lDX[5,15] lDX[6,0] lDX[6,1] lDX[6,2] lDX[6,3] lDX[6,4]
 lDX[6,5] lDX[6,6] lDX[6,7] lDX[6,8] lDX[6,9] lDX[6,10] lDX[6,11] lDX[6,12]
 lDX[6,13] lDX[6,14] lDX[6,15] lDX[7,0] lDX[7,1] lDX[7,2] lDX[7,3] lDX[7,4]
 lDX[7,5] lDX[7,6] lDX[7,7] lDX[7,8] lDX[7,9] lDX[7,10] lDX[7,11] lDX[7,12]
 lDX[7,13] lDX[7,14] lDX[7,15] lDX[8,0] lDX[8,1] lDX[8,2] lDX[8,3] lDX[8,4]
 lDX[8,5] lDX[8,6] lDX[8,7] lDX[8,8] lDX[8,9] lDX[8,10] lDX[8,11] lDX[8,12]
 lDX[8,13] lDX[8,14] lDX[8,15] ldx[1,0] ldx[1,1] ldx[1,2] ldx[1,3] ldx[1,4]
 ldx[1,5] ldx[1,6] ldx[1,7] ldx[1,8] ldx[1,9] ldx[1,10] ldx[1,11] ldx[1,12]
 ldx[1,13] ldx[1,14] ldx[1,15] ldx[2,0] ldx[2,1] ldx[2,2] ldx[2,3] ldx[2,4]
 ldx[2,5] ldx[2,6] ldx[2,7] ldx[2,8] ldx[2,9] ldx[2,10] ldx[2,11] ldx[2,12]
 ldx[2,13] ldx[2,14] ldx[2,15] ldx[3,0] ldx[3,1] ldx[3,2] ldx[3,3] ldx[3,4]
 ldx[3,5] ldx[3,6] ldx[3,7] ldx[3,8] ldx[3,9] ldx[3,10] ldx[3,11] ldx[3,12]
 ldx[3,13] ldx[3,14] ldx[3,15] ldx[4,0] ldx[4,1] ldx[4,2] ldx[4,3] ldx[4,4]
 ldx[4,5] ldx[4,6] ldx[4,7] ldx[4,8] ldx[4,9] ldx[4,10] ldx[4,11] ldx[4,12]
 ldx[4,13] ldx[4,14] ldx[4,15] ldx[5,0] ldx[5,1] ldx[5,2] ldx[5,3] ldx[5,4]
 ldx[5,5] ldx[5,6] ldx[5,7] ldx[5,8] ldx[5,9] ldx[5,10] ldx[5,11] ldx[5,12]
 ldx[5,13] ldx[5,14] ldx[5,15] ldx[6,0] ldx[6,1] ldx[6,2] ldx[6,3] ldx[6,4]
 ldx[6,5] ldx[6,6] ldx[6,7] ldx[6,8] ldx[6,9] ldx[6,10] ldx[6,11] ldx[6,12]
 ldx[6,13] ldx[6,14] ldx[6,15] ldx[7,0] ldx[7,1] ldx[7,2] ldx[7,3] ldx[7,4]
 ldx[7,5] ldx[7,6] ldx[7,7] ldx[7,8] ldx[7,9] ldx[7,10] ldx[7,11] ldx[7,12]
 ldx[7,13] ldx[7,14] ldx[7,15] ldx[8,0] ldx[8,1] ldx[8,2] ldx[8,3] ldx[8,4]
 ldx[8,5] ldx[8,6] ldx[8,7] ldx[8,8] ldx[8,9] ldx[8,10] ldx[8,11] ldx[8,12]
 ldx[8,13] ldx[8,14] ldx[8,15] lDY[1,0] lDY[1,1] lDY[1,2] lDY[1,3] lDY[1,4]
 lDY[1,5] lDY[1,6] lDY[1,7] lDY[1,8] lDY[1,9] lDY[1,10] lDY[1,11] lDY[1,12]
 lDY[1,13] lDY[1,14] lDY[1,15] lDY[2,0] lDY[2,1] lDY[2,2] lDY[2,3] lDY[2,4]
 lDY[2,5] lDY[2,6] lDY[2,7] lDY[2,8] lDY[2,9] lDY[2,10] lDY[2,11] lDY[2,12]
 lDY[2,13] lDY[2,14] lDY[2,15] lDY[3,0] lDY[3,1] lDY[3,2] lDY[3,3] lDY[3,4]
 lDY[3,5] lDY[3,6] lDY[3,7] lDY[3,8] lDY[3,9] lDY[3,10] lDY[3,11] lDY[3,12]
 lDY[3,13] lDY[3,14] lDY[3,15] lDY[4,0] lDY[4,1] lDY[4,2] lDY[4,3] lDY[4,4]
 lDY[4,5] lDY[4,6] lDY[4,7] lDY[4,8] lDY[4,9] lDY[4,10] lDY[4,11] lDY[4,12]
 lDY[4,13] lDY[4,14] lDY[4,15] lDY[5,0] lDY[5,1] lDY[5,2] lDY[5,3] lDY[5,4]
 lDY[5,5] lDY[5,6] lDY[5,7] lDY[5,8] lDY[5,9] lDY[5,10] lDY[5,11] lDY[5,12]
 lDY[5,13] lDY[5,14] lDY[5,15] lDY[6,0] lDY[6,1] lDY[6,2] lDY[6,3] lDY[6,4]
 lDY[6,5] lDY[6,6] lDY[6,7] lDY[6,8] lDY[6,9] lDY[6,10] lDY[6,11] lDY[6,12]
 lDY[6,13] lDY[6,14] lDY[6,15] lDY[7,0] lDY[7,1] lDY[7,2] lDY[7,3] lDY[7,4]
 lDY[7,5] lDY[7,6] lDY[7,7] lDY[7,8] lDY[7,9] lDY[7,10] lDY[7,11] lDY[7,12]
 lDY[7,13] lDY[7,14] lDY[7,15] ldy[1,0] ldy[1,1] ldy[1,2] ldy[1,3] ldy[1,4]
 ldy[1,5] ldy[1,6] ldy[1,7] ldy[1,8] ldy[1,9] ldy[1,10] ldy[1,11] ldy[1,12]
 ldy[1,13] ldy[1,14] ldy[1,15] ldy[2,0] ldy[2,1] ldy[2,2] ldy[2,3] ldy[2,4]
 ldy[2,5] ldy[2,6] ldy[2,7] ldy[2,8] ldy[2,9] ldy[2,10] ldy[2,11] ldy[2,12]
 ldy[2,13] ldy[2,14] ldy[2,15] ldy[3,0] ldy[3,1] ldy[3,2] ldy[3,3] ldy[3,4]
 ldy[3,5] ldy[3,6] ldy[3,7] ldy[3,8] ldy[3,9] ldy[3,10] ldy[3,11] ldy[3,12]
 ldy[3,13] ldy[3,14] ldy[3,15] ldy[4,0] ldy[4,1] ldy[4,2] ldy[4,3] ldy[4,4]
 ldy[4,5] ldy[4,6] ldy[4,7] ldy[4,8] ldy[4,9] ldy[4,10] ldy[4,11] ldy[4,12]
 ldy[4,13] ldy[4,14] ldy[4,15] ldy[5,0] ldy[5,1] ldy[5,2] ldy[5,3] ldy[5,4]
 ldy[5,5] ldy[5,6] ldy[5,7] ldy[5,8] ldy[5,9] ldy[5,10] ldy[5,11] ldy[5,12]
 ldy[5,13] ldy[5,14] ldy[5,15] ldy[6,0] ldy[6,1] ldy[6,2] ldy[6,3] ldy[6,4]
 ldy[6,5] ldy[6,6] ldy[6,7] ldy[6,8] ldy[6,9] ldy[6,10] ldy[6,11] ldy[6,12]
 ldy[6,13] ldy[6,14] ldy[6,15] ldy[7,0] ldy[7,1] ldy[7,2] ldy[7,3] ldy[7,4]
 ldy[7,5] ldy[7,6] ldy[7,7] ldy[7,8] ldy[7,9] ldy[7,10] ldy[7,11] ldy[7,12]
 ldy[7,13] ldy[7,14] ldy[7,15] lDZ[1,0] lDZ[1,1] lDZ[1,2] lDZ[1,3] lDZ[1,4]
 lDZ[1,5] lDZ[1,6] lDZ[1,7] lDZ[1,8] lDZ[1,9] lDZ[1,10] lDZ[1,11] lDZ[1,12]
 lDZ[1,13] lDZ[1,14] lDZ[1,15] lDZ[2,0] lDZ[2,1] lDZ[2,2] lDZ[2,3] lDZ[2,4]
 lDZ[2,5] lDZ[2,6] lDZ[2,7] lDZ[2,8] lDZ[2,9] lDZ[2,10] lDZ[2,11] lDZ[2,12]
 lDZ[2,13] lDZ[2,14] lDZ[2,15] lDZ[3,0] lDZ[3,1] lDZ[3,2] lDZ[3,3] lDZ[3,4]
 lDZ[3,5] lDZ[3,6] lDZ[3,7] lDZ[3,8] lDZ[3,9] lDZ[3,10] lDZ[3,11] lDZ[3,12]
 lDZ[3,13] lDZ[3,14] lDZ[3,15] lDZ[4,0] lDZ[4,1] lDZ[4,2] lDZ[4,3] lDZ[4,4]
 lDZ[4,5] lDZ[4,6] lDZ[4,7] lDZ[4,8] lDZ[4,9] lDZ[4,10] lDZ[4,11] lDZ[4,12]
 lDZ[4,13] lDZ[4,14] lDZ[4,15] lDZ[5,0] lDZ[5,1] lDZ[5,2] lDZ[5,3] lDZ[5,4]
 lDZ[5,5] lDZ[5,6] lDZ[5,7] lDZ[5,8] lDZ[5,9] lDZ[5,10] lDZ[5,11] lDZ[5,12]
 lDZ[5,13] lDZ[5,14] lDZ[5,15] lDZ[6,0] lDZ[6,1] lDZ[6,2] lDZ[6,3] lDZ[6,4]
 lDZ[6,5] lDZ[6,6] lDZ[6,7] lDZ[6,8] lDZ[6,9] lDZ[6,10] lDZ[6,11] lDZ[6,12]
 lDZ[6,13] lDZ[6,14] lDZ[6,15] lDZ[7,0] lDZ[7,1] lDZ[7,2] lDZ[7,3] lDZ[7,4]
 lDZ[7,5] lDZ[7,6] lDZ[7,7] lDZ[7,8] lDZ[7,9] lDZ[7,10] lDZ[7,11] lDZ[7,12]
 lDZ[7,13] lDZ[7,14] lDZ[7,15] ldz[1,0] ldz[1,1] ldz[1,2] ldz[1,3] ldz[1,4]
 ldz[1,5] ldz[1,6] ldz[1,7] ldz[1,8] ldz[1,9] ldz[1,10] ldz[1,11] ldz[1,12]
 ldz[1,13] ldz[1,14] ldz[1,15] ldz[2,0] ldz[2,1] ldz[2,2] ldz[2,3] ldz[2,4]
 ldz[2,5] ldz[2,6] ldz[2,7] ldz[2,8] ldz[2,9] ldz[2,10] ldz[2,11] ldz[2,12]
 ldz[2,13] ldz[2,14] ldz[2,15] ldz[3,0] ldz[3,1] ldz[3,2] ldz[3,3] ldz[3,4]
 ldz[3,5] ldz[3,6] ldz[3,7] ldz[3,8] ldz[3,9] ldz[3,10] ldz[3,11] ldz[3,12]
 ldz[3,13] ldz[3,14] ldz[3,15] ldz[4,0] ldz[4,1] ldz[4,2] ldz[4,3] ldz[4,4]
 ldz[4,5] ldz[4,6] ldz[4,7] ldz[4,8] ldz[4,9] ldz[4,10] ldz[4,11] ldz[4,12]
 ldz[4,13] ldz[4,14] ldz[4,15] ldz[5,0] ldz[5,1] ldz[5,2] ldz[5,3] ldz[5,4]
 ldz[5,5] ldz[5,6] ldz[5,7] ldz[5,8] ldz[5,9] ldz[5,10] ldz[5,11] ldz[5,12]
 ldz[5,13] ldz[5,14] ldz[5,15] ldz[6,0] ldz[6,1] ldz[6,2] ldz[6,3] ldz[6,4]
 ldz[6,5] ldz[6,6] ldz[6,7] ldz[6,8] ldz[6,9] ldz[6,10] ldz[6,11] ldz[6,12]
 ldz[6,13] ldz[6,14] ldz[6,15] ldz[7,0] ldz[7,1] ldz[7,2] ldz[7,3] ldz[7,4]
 ldz[7,5] ldz[7,6] ldz[7,7] ldz[7,8] ldz[7,9] ldz[7,10] ldz[7,11] ldz[7,12]
 ldz[7,13] ldz[7,14] ldz[7,15] lstk[0,0] lstk[0,1] lstk[0,2] lstk[0,3]
 lstk[0,4] lstk[0,5] lstk[0,6] lstk[0,7] lstk[0,8] lstk[0,9] lstk[0,10]
 lstk[0,11] lstk[0,12] lstk[0,13] lstk[0,14] lstk[0,15] lstk[1,0] lstk[1,1]
 lstk[1,2] lstk[1,3] lstk[1,4] lstk[1,5] lstk[1,6] lstk[1,7] lstk[1,8]
 lstk[1,9] lstk[1,10] lstk[1,11] lstk[1,12] lstk[1,13] lstk[1,14]
 lstk[1,15] lstk[2,0] lstk[2,1] lstk[2,2] lstk[2,3] lstk[2,4] lstk[2,5]
 lstk[2,6] lstk[2,7] lstk[2,8] lstk[2,9] lstk[2,10] lstk[2,11] lstk[2,12]
 lstk[2,13] lstk[2,14] lstk[2,15] lstk[3,0] lstk[3,1] lstk[3,2] lstk[3,3]
 lstk[3,4] lstk[3,5] lstk[3,6] lstk[3,7] lstk[3,8] lstk[3,9] lstk[3,10]
 lstk[3,11] lstk[3,12] lstk[3,13] lstk[3,14] lstk[3,15] lstk[4,0] lstk[4,1]
 lstk[4,2] lstk[4,3] lstk[4,4] lstk[4,5] lstk[4,6] lstk[4,7] lstk[4,8]
 lstk[4,9] lstk[4,10] lstk[4,11] lstk[4,12] lstk[4,13] lstk[4,14]
 lstk[4,15] lstk[5,0] lstk[5,1] lstk[5,2] lstk[5,3] lstk[5,4] lstk[5,5]
 lstk[5,6] lstk[5,7] lstk[5,8] lstk[5,9] lstk[5,10] lstk[5,11] lstk[5,12]
 lstk[5,13] lstk[5,14] lstk[5,15] lstk[6,0] lstk[6,1] lstk[6,2] lstk[6,3]
 lstk[6,4] lstk[6,5] lstk[6,6] lstk[6,7] lstk[6,8] lstk[6,9] lstk[6,10]
 lstk[6,11] lstk[6,12] lstk[6,13] lstk[6,14] lstk[6,15] lstk[7,0] lstk[7,1]
 lstk[7,2] lstk[7,3] lstk[7,4] lstk[7,5] lstk[7,6] lstk[7,7] lstk[7,8]
 lstk[7,9] lstk[7,10] lstk[7,11] lstk[7,12] lstk[7,13] lstk[7,14]
 lstk[7,15] lstk[8,0] lstk[8,1] lstk[8,2] lstk[8,3] lstk[8,4] lstk[8,5]
 lstk[8,6] lstk[8,7] lstk[8,8] lstk[8,9] lstk[8,10] lstk[8,11] lstk[8,12]
 lstk[8,13] lstk[8,14] lstk[8,15] lstk[9,0] lstk[9,1] lstk[9,2] lstk[9,3]
 lstk[9,4] lstk[9,5] lstk[9,6] lstk[9,7] lstk[9,8] lstk[9,9] lstk[9,10]
 lstk[9,11] lstk[9,12] lstk[9,13] lstk[9,14] lstk[9,15] lstk[10,0]
 lstk[10,1] lstk[10,2] lstk[10,3] lstk[10,4] lstk[10,5] lstk[10,6]
 lstk[10,7] lstk[10,8] lstk[10,9] lstk[10,10] lstk[10,11] lstk[10,12]
 lstk[10,13] lstk[10,14] lstk[10,15] lcan[1,0] lcan[1,1] lcan[1,2]
 lcan[1,3] lcan[1,4] lcan[1,5] lcan[1,6] lcan[1,7] lcan[1,8] lcan[1,9]
 lcan[1,10] lcan[1,11] lcan[1,12] lcan[1,13] lcan[1,14] lcan[1,15]
 lcan[2,0] lcan[2,1] lcan[2,2] lcan[2,3] lcan[2,4] lcan[2,5] lcan[2,6]
 lcan[2,7] lcan[2,8] lcan[2,9] lcan[2,10] lcan[2,11] lcan[2,12] lcan[2,13]
 lcan[2,14] lcan[2,15] lcan[3,0] lcan[3,1] lcan[3,2] lcan[3,3] lcan[3,4]
 lcan[3,5] lcan[3,6] lcan[3,7] lcan[3,8] lcan[3,9] lcan[3,10] lcan[3,11]
 lcan[3,12] lcan[3,13] lcan[3,14] lcan[3,15] lcan[4,0] lcan[4,1] lcan[4,2]
 lcan[4,3] lcan[4,4] lcan[4,5] lcan[4,6] lcan[4,7] lcan[4,8] lcan[4,9]
 lcan[4,10] lcan[4,11] lcan[4,12] lcan[4,13] lcan[4,14] lcan[4,15]
 lcan[5,0] lcan[5,1] lcan[5,2] lcan[5,3] lcan[5,4] lcan[5,5] lcan[5,6]
 lcan[5,7] lcan[5,8] lcan[5,9] lcan[5,10] lcan[5,11] lcan[5,12] lcan[5,13]
 lcan[5,14] lcan[5,15] lcan[6,0] lcan[6,1] lcan[6,2] lcan[6,3] lcan[6,4]
 lcan[6,5] lcan[6,6] lcan[6,7] lcan[6,8] lcan[6,9] lcan[6,10] lcan[6,11]
 lcan[6,12] lcan[6,13] lcan[6,14] lcan[6,15] lcan[7,0] lcan[7,1] lcan[7,2]
 lcan[7,3] lcan[7,4] lcan[7,5] lcan[7,6] lcan[7,7] lcan[7,8] lcan[7,9]
 lcan[7,10] lcan[7,11] lcan[7,12] lcan[7,13] lcan[7,14] lcan[7,15]
 lcan[8,0] lcan[8,1] lcan[8,2] lcan[8,3] lcan[8,4] lcan[8,5] lcan[8,6]
 lcan[8,7] lcan[8,8] lcan[8,9] lcan[8,10] lcan[8,11] lcan[8,12] lcan[8,13]
 lcan[8,14] lcan[8,15] lAC[1,0] lAC[1,1] lAC[1,2] lAC[1,3] lAC[2,0]
 lAC[2,1] lAC[2,2] lAC[2,3] lAC[3,0] lAC[3,1] lAC[3,2] lAC[3,3] lAC[4,0]
 lAC[4,1] lAC[4,2] lAC[4,3] lAC[5,0] lAC[5,1] lAC[5,2] lAC[5,3] lAC[6,0]
 lAC[6,1] lAC[6,2] lAC[6,3] lAC[7,0] lAC[7,1] lAC[7,2] lAC[7,3] lAC[8,0]
 lAC[8,1] lAC[8,2] lAC[8,3] lTC[1,0] lTC[1,1] lTC[1,2] lTC[1,3] lTC[2,0]
 lTC[2,1] lTC[2,2] lTC[2,3] lTC[3,0] lTC[3,1] lTC[3,2] lTC[3,3] lTC[4,0]
 lTC[4,1] lTC[4,2] lTC[4,3] lTC[5,0] lTC[5,1] lTC[5,2] lTC[5,3] lTC[6,0]
 lTC[6,1] lTC[6,2] lTC[6,3] lTC[7,0] lTC[7,1] lTC[7,2] lTC[7,3] lTC[8,0]
 lTC[8,1] lTC[8,2] lTC[8,3] leqDW[8,0] leqDW[8,1] leqDW[8,2] leqDW[8,3]
 leqDW[8,4] leqDW[8,5] leqDW[8,6] leqDW[8,7] leqDW[8,8] leqDW[8,9]
 leqDW[8,10] leqDW[8,11] leqDW[8,12] leqDW[8,13] leqDW[8,14] leqDW[8,15]
 leqdw[8,0] leqdw[8,1] leqdw[8,2] leqdw[8,3] leqdw[8,4] leqdw[8,5]
 leqdw[8,6] leqdw[8,7] leqdw[8,8] leqdw[8,9] leqdw[8,10] leqdw[8,11]
 leqdw[8,12] leqdw[8,13] leqdw[8,14] leqdw[8,15] leqDX[9,0] leqDX[9,1]
 leqDX[9,2] leqDX[9,3] leqDX[9,4] leqDX[9,5] leqDX[9,6] leqDX[9,7]
 leqDX[9,8] leqDX[9,9] leqDX[9,10] leqDX[9,11] leqDX[9,12] leqDX[9,13]
 leqDX[9,14] leqDX[9,15] leqdx[9,0] leqdx[9,1] leqdx[9,2] leqdx[9,3]
 leqdx[9,4] leqdx[9,5] leqdx[9,6] leqdx[9,7] leqdx[9,8] leqdx[9,9]
 leqdx[9,10] leqdx[9,11] leqdx[9,12] leqdx[9,13] leqdx[9,14] leqdx[9,15]
 leqDY[8,0] leqDY[8,1] leqDY[8,2] leqDY[8,3] leqDY[8,4] leqDY[8,5]
 leqDY[8,6] leqDY[8,7] leqDY[8,8] leqDY[8,9] leqDY[8,10] leqDY[8,11]
 leqDY[8,12] leqDY[8,13] leqDY[8,14] leqDY[8,15] leqDY[9,0] leqDY[9,1]
 leqDY[9,2] leqDY[9,3] leqDY[9,4] leqDY[9,5] leqDY[9,6] leqDY[9,7]
 leqDY[9,8] leqDY[9,9] leqDY[9,10] leqDY[9,11] leqDY[9,12] leqDY[9,13]
 leqDY[9,14] leqDY[9,15] leqdy[8,0] leqdy[8,1] leqdy[8,2] leqdy[8,3]
 leqdy[8,4] leqdy[8,5] leqdy[8,6] leqdy[8,7] leqdy[8,8] leqdy[8,9]
 leqdy[8,10] leqdy[8,11] leqdy[8,12] leqdy[8,13] leqdy[8,14] leqdy[8,15]
 leqdy[9,0] leqdy[9,1] leqdy[9,2] leqdy[9,3] leqdy[9,4] leqdy[9,5]
 leqdy[9,6] leqdy[9,7] leqdy[9,8] leqdy[9,9] leqdy[9,10] leqdy[9,11]
 leqdy[9,12] leqdy[9,13] leqdy[9,14] leqdy[9,15] leqDZ[8,0] leqDZ[8,1]
 leqDZ[8,2] leqDZ[8,3] leqDZ[8,4] leqDZ[8,5] leqDZ[8,6] leqDZ[8,7]
 leqDZ[8,8] leqDZ[8,9] leqDZ[8,10] leqDZ[8,11] leqDZ[8,12] leqDZ[8,13]
 leqDZ[8,14] leqDZ[8,15] leqDZ[9,0] leqDZ[9,1] leqDZ[9,2] leqDZ[9,3]
 leqDZ[9,4] leqDZ[9,5] leqDZ[9,6] leqDZ[9,7] leqDZ[9,8] leqDZ[9,9]
 leqDZ[9,10] leqDZ[9,11] leqDZ[9,12] leqDZ[9,13] leqDZ[9,14] leqDZ[9,15]
 leqdz[8,0] leqdz[8,1] leqdz[8,2] leqdz[8,3] leqdz[8,4] leqdz[8,5]
 leqdz[8,6] leqdz[8,7] leqdz[8,8] leqdz[8,9] leqdz[8,10] leqdz[8,11]
 leqdz[8,12] leqdz[8,13] leqdz[8,14] leqdz[8,15] leqdz[9,0] leqdz[9,1]
 leqdz[9,2] leqdz[9,3] leqdz[9,4] leqdz[9,5] leqdz[9,6] leqdz[9,7]
 leqdz[9,8] leqdz[9,9] leqdz[9,10] leqdz[9,11] leqdz[9,12] leqdz[9,13]
 leqdz[9,14] leqdz[9,15] leqT[9,0] leqT[9,1] leqT[9,2] leqT[9,3] leqT[10,0]
 leqT[10,1] leqT[10,2] leqT[10,3] leqk[9,0] leqk[9,1] leqk[9,2] leqk[9,3]
 leqk[9,4] leqk[9,5] leqk[9,6] leqk[9,7] leqk[9,8] leqk[9,9] leqk[9,10]
 leqk[9,11] leqk[9,12] leqk[9,13] leqk[9,14] leqk[9,15] leqk[10,0]
 leqk[10,1] leqk[10,2] leqk[10,3] leqk[10,4] leqk[10,5] leqk[10,6]
 leqk[10,7] leqk[10,8] leqk[10,9] leqk[10,10] leqk[10,11] leqk[10,12]
 leqk[10,13] leqk[10,14] leqk[10,15] lLANE[0] lLANE[1] lLANE[2] lLANE[3]
 lLANE[4] lLANE[5] lLANE[6] lLANE[7] lLANE[8] lLANE[9] lLANE[10] lLANE[11]
 lLANE[12] lLANE[13] lLANE[14] lLANE[15] lGstk[9,0] lGstk[9,1] lGstk[9,2]
 lGstk[9,3] lGstk[9,4] lGstk[9,5] lGstk[9,6] lGstk[9,7] lGstk[9,8]
 lGstk[9,9] lGstk[9,10] lGstk[9,11] lGstk[9,12] lGstk[9,13] lGstk[9,14]
 lGstk[9,15] lGstk[10,0] lGstk[10,1] lGstk[10,2] lGstk[10,3] lGstk[10,4]
 lGstk[10,5] lGstk[10,6] lGstk[10,7] lGstk[10,8] lGstk[10,9] lGstk[10,10]
 lGstk[10,11] lGstk[10,12] lGstk[10,13] lGstk[10,14] lGstk[10,15]
 ldetEQX[8,0] ldetEQX[8,1] ldetEQX[8,2] ldetEQX[8,3] ldetEQX[8,4]
 ldetEQX[8,5] ldetEQX[8,6] ldetEQX[8,7] ldetEQX[8,8] ldetEQX[8,9]
 ldetEQX[8,10] ldetEQX[8,11] ldetEQX[8,12] ldetEQX[8,13] ldetEQX[8,14]
 ldetEQX[8,15] ldetEQX[9,0] ldetEQX[9,1] ldetEQX[9,2] ldetEQX[9,3]
 ldetEQX[9,4] ldetEQX[9,5] ldetEQX[9,6] ldetEQX[9,7] ldetEQX[9,8]
 ldetEQX[9,9] ldetEQX[9,10] ldetEQX[9,11] ldetEQX[9,12] ldetEQX[9,13]
 ldetEQX[9,14] ldetEQX[9,15] ldetEQY[8,0] ldetEQY[8,1] ldetEQY[8,2]
 ldetEQY[8,3] ldetEQY[8,4] ldetEQY[8,5] ldetEQY[8,6] ldetEQY[8,7]
 ldetEQY[8,8] ldetEQY[8,9] ldetEQY[8,10] ldetEQY[8,11] ldetEQY[8,12]
 ldetEQY[8,13] ldetEQY[8,14] ldetEQY[8,15] ldetEQY[9,0] ldetEQY[9,1]
 ldetEQY[9,2] ldetEQY[9,3] ldetEQY[9,4] ldetEQY[9,5] ldetEQY[9,6]
 ldetEQY[9,7] ldetEQY[9,8] ldetEQY[9,9] ldetEQY[9,10] ldetEQY[9,11]
 ldetEQY[9,12] ldetEQY[9,13] ldetEQY[9,14] ldetEQY[9,15] ldetEQZ[8,0]
 ldetEQZ[8,1] ldetEQZ[8,2] ldetEQZ[8,3] ldetEQZ[8,4] ldetEQZ[8,5]
 ldetEQZ[8,6] ldetEQZ[8,7] ldetEQZ[8,8] ldetEQZ[8,9] ldetEQZ[8,10]
 ldetEQZ[8,11] ldetEQZ[8,12] ldetEQZ[8,13] ldetEQZ[8,14] ldetEQZ[8,15]
 ldetEQZ[9,0] ldetEQZ[9,1] ldetEQZ[9,2] ldetEQZ[9,3] ldetEQZ[9,4]
 ldetEQZ[9,5] ldetEQZ[9,6] ldetEQZ[9,7] ldetEQZ[9,8] ldetEQZ[9,9]
 ldetEQZ[9,10] ldetEQZ[9,11] ldetEQZ[9,12] ldetEQZ[9,13] ldetEQZ[9,14]
 ldetEQZ[9,15] ldetEQW[8,0] ldetEQW[8,1] ldetEQW[8,2] ldetEQW[8,3]
 ldetEQW[8,4] ldetEQW[8,5] ldetEQW[8,6] ldetEQW[8,7] ldetEQW[8,8]
 ldetEQW[8,9] ldetEQW[8,10] ldetEQW[8,11] ldetEQW[8,12] ldetEQW[8,13]
 ldetEQW[8,14] ldetEQW[8,15] lDetC[9,0] lDetC[9,1] lDetC[9,2] lDetC[9,3]
 lFrSB[8,0] lFrSB[8,1] lFrSB[8,2] lFrSB[8,3] lFrSB[8,4] lFrSB[8,5]
 lFrSB[8,6] lFrSB[8,7] lFrSB[8,8] lFrSB[8,9] lFrSB[8,10] lFrSB[8,11]
 lFrSB[8,12] lFrSB[8,13] lFrSB[8,14] lFrSB[8,15] lFrSB[9,0] lFrSB[9,1]
 lFrSB[9,2] lFrSB[9,3] lFrSB[9,4] lFrSB[9,5] lFrSB[9,6] lFrSB[9,7]
 lFrSB[9,8] lFrSB[9,9] lFrSB[9,10] lFrSB[9,11] lFrSB[9,12] lFrSB[9,13]
 lFrSB[9,14] lFrSB[9,15] lDiffDetMC[8,0] lDiffDetMC[8,1] lDiffDetMC[8,2]
 lDiffDetMC[8,3] lDiffDetMC[8,4] lDiffDetMC[8,5] lDiffDetMC[8,6]
 lDiffDetMC[8,7] lDiffDetMC[8,8] lDiffDetMC[8,9] lDiffDetMC[8,10]
 lDiffDetMC[8,11] lDiffDetMC[8,12] lDiffDetMC[8,13] lDiffDetMC[8,14]
 lDiffDetMC[8,15] lDiffDetC[8,0] lDiffDetC[8,1] lDiffDetC[8,2]
 lDiffDetC[8,3] lFrMC[8,0] lFrMC[8,1] lFrMC[8,2] lFrMC[8,3] lFrMC[8,4]
 lFrMC[8,5] lFrMC[8,6] lFrMC[8,7] lFrMC[8,8] lFrMC[8,9] lFrMC[8,10]
 lFrMC[8,11] lFrMC[8,12] lFrMC[8,13] lFrMC[8,14] lFrMC[8,15]
Generals
 uT2CanA[0,0] uT2CanA[0,1] uT2CanA[0,2] uT2CanA[0,3] uT2CanA[1,0]
 uT2CanA[1,1] uT2CanA[1,2] uT2CanA[1,3] uT2CanA[2,0] uT2CanA[2,1]
 uT2CanA[2,2] uT2CanA[2,3] uT2CanA[3,0] uT2CanA[3,1] uT2CanA[3,2]
 uT2CanA[3,3] uT2CanA[4,0] uT2CanA[4,1] uT2CanA[4,2] uT2CanA[4,3]
 uT2CanA[5,0] uT2CanA[5,1] uT2CanA[5,2] uT2CanA[5,3] uT2CanA[6,0]
 uT2CanA[6,1] uT2CanA[6,2] uT2CanA[6,3] uT2CanA[7,0] uT2CanA[7,1]
 uT2CanA[7,2] uT2CanA[7,3] uT2CanB[0,0] uT2CanB[0,1] uT2CanB[0,2]
 uT2CanB[0,3] uT2CanB[1,0] uT2CanB[1,1] uT2CanB[1,2] uT2CanB[1,3]
 uT2CanB[2,0] uT2CanB[2,1] uT2CanB[2,2] uT2CanB[2,3] uT2CanB[3,0]
 uT2CanB[3,1] uT2CanB[3,2] uT2CanB[3,3] uT2CanB[4,0] uT2CanB[4,1]
 uT2CanB[4,2] uT2CanB[4,3] uT2CanB[5,0] uT2CanB[5,1] uT2CanB[5,2]
 uT2CanB[5,3] uT2CanB[6,0] uT2CanB[6,1] uT2CanB[6,2] uT2CanB[6,3]
 uT2CanB[7,0] uT2CanB[7,1] uT2CanB[7,2] uT2CanB[7,3] uT2CanC[1,0]
 uT2CanC[1,1] uT2CanC[1,2] uT2CanC[1,3] uT2CanC[2,0] uT2CanC[2,1]
 uT2CanC[2,2] uT2CanC[2,3] uT2CanC[3,0] uT2CanC[3,1] uT2CanC[3,2]
 uT2CanC[3,3] uT2CanC[4,0] uT2CanC[4,1] uT2CanC[4,2] uT2CanC[4,3]
 uT2CanC[5,0] uT2CanC[5,1] uT2CanC[5,2] uT2CanC[5,3] uT2CanC[6,0]
 uT2CanC[6,1] uT2CanC[6,2] uT2CanC[6,3] uT2CanC[7,0] uT2CanC[7,1]
 uT2CanC[7,2] uT2CanC[7,3] uType1[0] uType1[1] uType1[2] uType1[3]
 uType1[4] uType1[5] uType1[6] uType1[7] uType1[8] uType1[9] uType1[10]
 uType1[11] uType1[12] uType1[13] uType1[14] uType1[15] uType2[0,0]
 uType2[0,1] uType2[0,2] uType2[0,3] uType2[1,0] uType2[1,1] uType2[1,2]
 uType2[1,3] uType2[2,0] uType2[2,1] uType2[2,2] uType2[2,3] uType2[3,0]
 uType2[3,1] uType2[3,2] uType2[3,3] uType2[4,0] uType2[4,1] uType2[4,2]
 uType2[4,3] uType2[5,0] uType2[5,1] uType2[5,2] uType2[5,3] uType2[6,0]
 uType2[6,1] uType2[6,2] uType2[6,3] uType2[7,0] uType2[7,1] uType2[7,2]
 uType2[7,3] lT2CanA[0,0] lT2CanA[0,1] lT2CanA[0,2] lT2CanA[0,3]
 lT2CanA[1,0] lT2CanA[1,1] lT2CanA[1,2] lT2CanA[1,3] lT2CanA[2,0]
 lT2CanA[2,1] lT2CanA[2,2] lT2CanA[2,3] lT2CanA[3,0] lT2CanA[3,1]
 lT2CanA[3,2] lT2CanA[3,3] lT2CanA[4,0] lT2CanA[4,1] lT2CanA[4,2]
 lT2CanA[4,3] lT2CanA[5,0] lT2CanA[5,1] lT2CanA[5,2] lT2CanA[5,3]
 lT2CanA[6,0] lT2CanA[6,1] lT2CanA[6,2] lT2CanA[6,3] lT2CanA[7,0]
 lT2CanA[7,1] lT2CanA[7,2] lT2CanA[7,3] lT2CanB[0,0] lT2CanB[0,1]
 lT2CanB[0,2] lT2CanB[0,3] lT2CanB[1,0] lT2CanB[1,1] lT2CanB[1,2]
 lT2CanB[1,3] lT2CanB[2,0] lT2CanB[2,1] lT2CanB[2,2] lT2CanB[2,3]
 lT2CanB[3,0] lT2CanB[3,1] lT2CanB[3,2] lT2CanB[3,3] lT2CanB[4,0]
 lT2CanB[4,1] lT2CanB[4,2] lT2CanB[4,3] lT2CanB[5,0] lT2CanB[5,1]
 lT2CanB[5,2] lT2CanB[5,3] lT2CanB[6,0] lT2CanB[6,1] lT2CanB[6,2]
 lT2CanB[6,3] lT2CanB[7,0] lT2CanB[7,1] lT2CanB[7,2] lT2CanB[7,3]
 lT2CanC[1,0] lT2CanC[1,1] lT2CanC[1,2] lT2CanC[1,3] lT2CanC[2,0]
 lT2CanC[2,1] lT2CanC[2,2] lT2CanC[2,3] lT2CanC[3,0] lT2CanC[3,1]
 lT2CanC[3,2] lT2CanC[3,3] lT2CanC[4,0] lT2CanC[4,1] lT2CanC[4,2]
 lT2CanC[4,3] lT2CanC[5,0] lT2CanC[5,1] lT2CanC[5,2] lT2CanC[5,3]
 lT2CanC[6,0] lT2CanC[6,1] lT2CanC[6,2] lT2CanC[6,3] lT2CanC[7,0]
 lT2CanC[7,1] lT2CanC[7,2] lT2CanC[7,3] lT2CanC[8,0] lT2CanC[8,1]
 lT2CanC[8,2] lT2CanC[8,3] lType1[0] lType1[1] lType1[2] lType1[3]
 lType1[4] lType1[5] lType1[6] lType1[7] lType1[8] lType1[9] lType1[10]
 lType1[11] lType1[12] lType1[13] lType1[14] lType1[15] lType2[0,0]
 lType2[0,1] lType2[0,2] lType2[0,3] lType2[1,0] lType2[1,1] lType2[1,2]
 lType2[1,3] lType2[2,0] lType2[2,1] lType2[2,2] lType2[2,3] lType2[3,0]
 lType2[3,1] lType2[3,2] lType2[3,3] lType2[4,0] lType2[4,1] lType2[4,2]
 lType2[4,3] lType2[5,0] lType2[5,1] lType2[5,2] lType2[5,3] lType2[6,0]
 lType2[6,1] lType2[6,2] lType2[6,3] lType2[7,0] lType2[7,1] lType2[7,2]
 lType2[7,3] rb rf cb cf mb mf mb' mf' cb' cf'
End
