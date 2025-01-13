\ Model ./SKINNYe-64-256_33-6+11+11+5_False
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
  OBJ4: Priority=0 Weight=1 AbsTol=1e-06 RelTol=0
  
  OBJ5: Priority=0 Weight=1 AbsTol=1e-06 RelTol=0
  
  OBJ6: Priority=0 Weight=1 AbsTol=1e-06 RelTol=0
  
  OBJ7: Priority=0 Weight=1 AbsTol=1e-06 RelTol=0
  
  OBJ8: Priority=0 Weight=1 AbsTol=1e-06 RelTol=0
  
  OBJ9: Priority=0 Weight=1 AbsTol=1e-06 RelTol=0
  
  mbi: Priority=0 Weight=1 AbsTol=1e-06 RelTol=0
   mbi
  mfi: Priority=0 Weight=1 AbsTol=1e-06 RelTol=0
   mfi
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
 uBasic: uDX[13,0] - udx[13,0] >= 0
 uBasic: uDX[14,0] - udx[14,0] >= 0
 uBasic: uDX[15,0] - udx[15,0] >= 0
 uBasic: uDX[16,0] - udx[16,0] >= 0
 uBasic: uDX[17,0] - udx[17,0] >= 0
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
 uBasic: uDY[12,0] - udy[12,0] >= 0
 uBasic: uDZ[12,0] - udz[12,0] >= 0
 uBasic: uDY[13,0] - udy[13,0] >= 0
 uBasic: uDZ[13,0] - udz[13,0] >= 0
 uBasic: uDY[14,0] - udy[14,0] >= 0
 uBasic: uDZ[14,0] - udz[14,0] >= 0
 uBasic: uDY[15,0] - udy[15,0] >= 0
 uBasic: uDZ[15,0] - udz[15,0] >= 0
 uBasic: uDY[16,0] - udy[16,0] >= 0
 uBasic: uDZ[16,0] - udz[16,0] >= 0
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
 uBasic: uDW[12,0] - udw[12,0] >= 0
 uBasic: uDW[13,0] - udw[13,0] >= 0
 uBasic: uDW[14,0] - udw[14,0] >= 0
 uBasic: uDW[15,0] - udw[15,0] >= 0
 uBasic: uDW[16,0] - udw[16,0] >= 0
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
 uBasic: uDX[13,1] - udx[13,1] >= 0
 uBasic: uDX[14,1] - udx[14,1] >= 0
 uBasic: uDX[15,1] - udx[15,1] >= 0
 uBasic: uDX[16,1] - udx[16,1] >= 0
 uBasic: uDX[17,1] - udx[17,1] >= 0
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
 uBasic: uDY[12,1] - udy[12,1] >= 0
 uBasic: uDZ[12,1] - udz[12,1] >= 0
 uBasic: uDY[13,1] - udy[13,1] >= 0
 uBasic: uDZ[13,1] - udz[13,1] >= 0
 uBasic: uDY[14,1] - udy[14,1] >= 0
 uBasic: uDZ[14,1] - udz[14,1] >= 0
 uBasic: uDY[15,1] - udy[15,1] >= 0
 uBasic: uDZ[15,1] - udz[15,1] >= 0
 uBasic: uDY[16,1] - udy[16,1] >= 0
 uBasic: uDZ[16,1] - udz[16,1] >= 0
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
 uBasic: uDW[12,1] - udw[12,1] >= 0
 uBasic: uDW[13,1] - udw[13,1] >= 0
 uBasic: uDW[14,1] - udw[14,1] >= 0
 uBasic: uDW[15,1] - udw[15,1] >= 0
 uBasic: uDW[16,1] - udw[16,1] >= 0
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
 uBasic: uDX[13,2] - udx[13,2] >= 0
 uBasic: uDX[14,2] - udx[14,2] >= 0
 uBasic: uDX[15,2] - udx[15,2] >= 0
 uBasic: uDX[16,2] - udx[16,2] >= 0
 uBasic: uDX[17,2] - udx[17,2] >= 0
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
 uBasic: uDY[12,2] - udy[12,2] >= 0
 uBasic: uDZ[12,2] - udz[12,2] >= 0
 uBasic: uDY[13,2] - udy[13,2] >= 0
 uBasic: uDZ[13,2] - udz[13,2] >= 0
 uBasic: uDY[14,2] - udy[14,2] >= 0
 uBasic: uDZ[14,2] - udz[14,2] >= 0
 uBasic: uDY[15,2] - udy[15,2] >= 0
 uBasic: uDZ[15,2] - udz[15,2] >= 0
 uBasic: uDY[16,2] - udy[16,2] >= 0
 uBasic: uDZ[16,2] - udz[16,2] >= 0
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
 uBasic: uDW[12,2] - udw[12,2] >= 0
 uBasic: uDW[13,2] - udw[13,2] >= 0
 uBasic: uDW[14,2] - udw[14,2] >= 0
 uBasic: uDW[15,2] - udw[15,2] >= 0
 uBasic: uDW[16,2] - udw[16,2] >= 0
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
 uBasic: uDX[13,3] - udx[13,3] >= 0
 uBasic: uDX[14,3] - udx[14,3] >= 0
 uBasic: uDX[15,3] - udx[15,3] >= 0
 uBasic: uDX[16,3] - udx[16,3] >= 0
 uBasic: uDX[17,3] - udx[17,3] >= 0
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
 uBasic: uDY[12,3] - udy[12,3] >= 0
 uBasic: uDZ[12,3] - udz[12,3] >= 0
 uBasic: uDY[13,3] - udy[13,3] >= 0
 uBasic: uDZ[13,3] - udz[13,3] >= 0
 uBasic: uDY[14,3] - udy[14,3] >= 0
 uBasic: uDZ[14,3] - udz[14,3] >= 0
 uBasic: uDY[15,3] - udy[15,3] >= 0
 uBasic: uDZ[15,3] - udz[15,3] >= 0
 uBasic: uDY[16,3] - udy[16,3] >= 0
 uBasic: uDZ[16,3] - udz[16,3] >= 0
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
 uBasic: uDW[12,3] - udw[12,3] >= 0
 uBasic: uDW[13,3] - udw[13,3] >= 0
 uBasic: uDW[14,3] - udw[14,3] >= 0
 uBasic: uDW[15,3] - udw[15,3] >= 0
 uBasic: uDW[16,3] - udw[16,3] >= 0
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
 uBasic: uDX[13,4] - udx[13,4] >= 0
 uBasic: uDX[14,4] - udx[14,4] >= 0
 uBasic: uDX[15,4] - udx[15,4] >= 0
 uBasic: uDX[16,4] - udx[16,4] >= 0
 uBasic: uDX[17,4] - udx[17,4] >= 0
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
 uBasic: uDY[12,4] - udy[12,4] >= 0
 uBasic: uDZ[12,4] - udz[12,4] >= 0
 uBasic: uDY[13,4] - udy[13,4] >= 0
 uBasic: uDZ[13,4] - udz[13,4] >= 0
 uBasic: uDY[14,4] - udy[14,4] >= 0
 uBasic: uDZ[14,4] - udz[14,4] >= 0
 uBasic: uDY[15,4] - udy[15,4] >= 0
 uBasic: uDZ[15,4] - udz[15,4] >= 0
 uBasic: uDY[16,4] - udy[16,4] >= 0
 uBasic: uDZ[16,4] - udz[16,4] >= 0
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
 uBasic: uDW[12,4] - udw[12,4] >= 0
 uBasic: uDW[13,4] - udw[13,4] >= 0
 uBasic: uDW[14,4] - udw[14,4] >= 0
 uBasic: uDW[15,4] - udw[15,4] >= 0
 uBasic: uDW[16,4] - udw[16,4] >= 0
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
 uBasic: uDX[13,5] - udx[13,5] >= 0
 uBasic: uDX[14,5] - udx[14,5] >= 0
 uBasic: uDX[15,5] - udx[15,5] >= 0
 uBasic: uDX[16,5] - udx[16,5] >= 0
 uBasic: uDX[17,5] - udx[17,5] >= 0
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
 uBasic: uDY[12,5] - udy[12,5] >= 0
 uBasic: uDZ[12,5] - udz[12,5] >= 0
 uBasic: uDY[13,5] - udy[13,5] >= 0
 uBasic: uDZ[13,5] - udz[13,5] >= 0
 uBasic: uDY[14,5] - udy[14,5] >= 0
 uBasic: uDZ[14,5] - udz[14,5] >= 0
 uBasic: uDY[15,5] - udy[15,5] >= 0
 uBasic: uDZ[15,5] - udz[15,5] >= 0
 uBasic: uDY[16,5] - udy[16,5] >= 0
 uBasic: uDZ[16,5] - udz[16,5] >= 0
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
 uBasic: uDW[12,5] - udw[12,5] >= 0
 uBasic: uDW[13,5] - udw[13,5] >= 0
 uBasic: uDW[14,5] - udw[14,5] >= 0
 uBasic: uDW[15,5] - udw[15,5] >= 0
 uBasic: uDW[16,5] - udw[16,5] >= 0
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
 uBasic: uDX[13,6] - udx[13,6] >= 0
 uBasic: uDX[14,6] - udx[14,6] >= 0
 uBasic: uDX[15,6] - udx[15,6] >= 0
 uBasic: uDX[16,6] - udx[16,6] >= 0
 uBasic: uDX[17,6] - udx[17,6] >= 0
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
 uBasic: uDY[12,6] - udy[12,6] >= 0
 uBasic: uDZ[12,6] - udz[12,6] >= 0
 uBasic: uDY[13,6] - udy[13,6] >= 0
 uBasic: uDZ[13,6] - udz[13,6] >= 0
 uBasic: uDY[14,6] - udy[14,6] >= 0
 uBasic: uDZ[14,6] - udz[14,6] >= 0
 uBasic: uDY[15,6] - udy[15,6] >= 0
 uBasic: uDZ[15,6] - udz[15,6] >= 0
 uBasic: uDY[16,6] - udy[16,6] >= 0
 uBasic: uDZ[16,6] - udz[16,6] >= 0
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
 uBasic: uDW[12,6] - udw[12,6] >= 0
 uBasic: uDW[13,6] - udw[13,6] >= 0
 uBasic: uDW[14,6] - udw[14,6] >= 0
 uBasic: uDW[15,6] - udw[15,6] >= 0
 uBasic: uDW[16,6] - udw[16,6] >= 0
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
 uBasic: uDX[13,7] - udx[13,7] >= 0
 uBasic: uDX[14,7] - udx[14,7] >= 0
 uBasic: uDX[15,7] - udx[15,7] >= 0
 uBasic: uDX[16,7] - udx[16,7] >= 0
 uBasic: uDX[17,7] - udx[17,7] >= 0
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
 uBasic: uDY[12,7] - udy[12,7] >= 0
 uBasic: uDZ[12,7] - udz[12,7] >= 0
 uBasic: uDY[13,7] - udy[13,7] >= 0
 uBasic: uDZ[13,7] - udz[13,7] >= 0
 uBasic: uDY[14,7] - udy[14,7] >= 0
 uBasic: uDZ[14,7] - udz[14,7] >= 0
 uBasic: uDY[15,7] - udy[15,7] >= 0
 uBasic: uDZ[15,7] - udz[15,7] >= 0
 uBasic: uDY[16,7] - udy[16,7] >= 0
 uBasic: uDZ[16,7] - udz[16,7] >= 0
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
 uBasic: uDW[12,7] - udw[12,7] >= 0
 uBasic: uDW[13,7] - udw[13,7] >= 0
 uBasic: uDW[14,7] - udw[14,7] >= 0
 uBasic: uDW[15,7] - udw[15,7] >= 0
 uBasic: uDW[16,7] - udw[16,7] >= 0
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
 uBasic: uDX[13,8] - udx[13,8] >= 0
 uBasic: uDX[14,8] - udx[14,8] >= 0
 uBasic: uDX[15,8] - udx[15,8] >= 0
 uBasic: uDX[16,8] - udx[16,8] >= 0
 uBasic: uDX[17,8] - udx[17,8] >= 0
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
 uBasic: uDY[12,8] - udy[12,8] >= 0
 uBasic: uDZ[12,8] - udz[12,8] >= 0
 uBasic: uDY[13,8] - udy[13,8] >= 0
 uBasic: uDZ[13,8] - udz[13,8] >= 0
 uBasic: uDY[14,8] - udy[14,8] >= 0
 uBasic: uDZ[14,8] - udz[14,8] >= 0
 uBasic: uDY[15,8] - udy[15,8] >= 0
 uBasic: uDZ[15,8] - udz[15,8] >= 0
 uBasic: uDY[16,8] - udy[16,8] >= 0
 uBasic: uDZ[16,8] - udz[16,8] >= 0
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
 uBasic: uDW[12,8] - udw[12,8] >= 0
 uBasic: uDW[13,8] - udw[13,8] >= 0
 uBasic: uDW[14,8] - udw[14,8] >= 0
 uBasic: uDW[15,8] - udw[15,8] >= 0
 uBasic: uDW[16,8] - udw[16,8] >= 0
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
 uBasic: uDX[13,9] - udx[13,9] >= 0
 uBasic: uDX[14,9] - udx[14,9] >= 0
 uBasic: uDX[15,9] - udx[15,9] >= 0
 uBasic: uDX[16,9] - udx[16,9] >= 0
 uBasic: uDX[17,9] - udx[17,9] >= 0
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
 uBasic: uDY[12,9] - udy[12,9] >= 0
 uBasic: uDZ[12,9] - udz[12,9] >= 0
 uBasic: uDY[13,9] - udy[13,9] >= 0
 uBasic: uDZ[13,9] - udz[13,9] >= 0
 uBasic: uDY[14,9] - udy[14,9] >= 0
 uBasic: uDZ[14,9] - udz[14,9] >= 0
 uBasic: uDY[15,9] - udy[15,9] >= 0
 uBasic: uDZ[15,9] - udz[15,9] >= 0
 uBasic: uDY[16,9] - udy[16,9] >= 0
 uBasic: uDZ[16,9] - udz[16,9] >= 0
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
 uBasic: uDW[12,9] - udw[12,9] >= 0
 uBasic: uDW[13,9] - udw[13,9] >= 0
 uBasic: uDW[14,9] - udw[14,9] >= 0
 uBasic: uDW[15,9] - udw[15,9] >= 0
 uBasic: uDW[16,9] - udw[16,9] >= 0
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
 uBasic: uDX[13,10] - udx[13,10] >= 0
 uBasic: uDX[14,10] - udx[14,10] >= 0
 uBasic: uDX[15,10] - udx[15,10] >= 0
 uBasic: uDX[16,10] - udx[16,10] >= 0
 uBasic: uDX[17,10] - udx[17,10] >= 0
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
 uBasic: uDY[12,10] - udy[12,10] >= 0
 uBasic: uDZ[12,10] - udz[12,10] >= 0
 uBasic: uDY[13,10] - udy[13,10] >= 0
 uBasic: uDZ[13,10] - udz[13,10] >= 0
 uBasic: uDY[14,10] - udy[14,10] >= 0
 uBasic: uDZ[14,10] - udz[14,10] >= 0
 uBasic: uDY[15,10] - udy[15,10] >= 0
 uBasic: uDZ[15,10] - udz[15,10] >= 0
 uBasic: uDY[16,10] - udy[16,10] >= 0
 uBasic: uDZ[16,10] - udz[16,10] >= 0
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
 uBasic: uDW[12,10] - udw[12,10] >= 0
 uBasic: uDW[13,10] - udw[13,10] >= 0
 uBasic: uDW[14,10] - udw[14,10] >= 0
 uBasic: uDW[15,10] - udw[15,10] >= 0
 uBasic: uDW[16,10] - udw[16,10] >= 0
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
 uBasic: uDX[13,11] - udx[13,11] >= 0
 uBasic: uDX[14,11] - udx[14,11] >= 0
 uBasic: uDX[15,11] - udx[15,11] >= 0
 uBasic: uDX[16,11] - udx[16,11] >= 0
 uBasic: uDX[17,11] - udx[17,11] >= 0
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
 uBasic: uDY[12,11] - udy[12,11] >= 0
 uBasic: uDZ[12,11] - udz[12,11] >= 0
 uBasic: uDY[13,11] - udy[13,11] >= 0
 uBasic: uDZ[13,11] - udz[13,11] >= 0
 uBasic: uDY[14,11] - udy[14,11] >= 0
 uBasic: uDZ[14,11] - udz[14,11] >= 0
 uBasic: uDY[15,11] - udy[15,11] >= 0
 uBasic: uDZ[15,11] - udz[15,11] >= 0
 uBasic: uDY[16,11] - udy[16,11] >= 0
 uBasic: uDZ[16,11] - udz[16,11] >= 0
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
 uBasic: uDW[12,11] - udw[12,11] >= 0
 uBasic: uDW[13,11] - udw[13,11] >= 0
 uBasic: uDW[14,11] - udw[14,11] >= 0
 uBasic: uDW[15,11] - udw[15,11] >= 0
 uBasic: uDW[16,11] - udw[16,11] >= 0
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
 uBasic: uDX[13,12] - udx[13,12] >= 0
 uBasic: uDX[14,12] - udx[14,12] >= 0
 uBasic: uDX[15,12] - udx[15,12] >= 0
 uBasic: uDX[16,12] - udx[16,12] >= 0
 uBasic: uDX[17,12] - udx[17,12] >= 0
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
 uBasic: uDY[12,12] - udy[12,12] >= 0
 uBasic: uDZ[12,12] - udz[12,12] >= 0
 uBasic: uDY[13,12] - udy[13,12] >= 0
 uBasic: uDZ[13,12] - udz[13,12] >= 0
 uBasic: uDY[14,12] - udy[14,12] >= 0
 uBasic: uDZ[14,12] - udz[14,12] >= 0
 uBasic: uDY[15,12] - udy[15,12] >= 0
 uBasic: uDZ[15,12] - udz[15,12] >= 0
 uBasic: uDY[16,12] - udy[16,12] >= 0
 uBasic: uDZ[16,12] - udz[16,12] >= 0
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
 uBasic: uDW[12,12] - udw[12,12] >= 0
 uBasic: uDW[13,12] - udw[13,12] >= 0
 uBasic: uDW[14,12] - udw[14,12] >= 0
 uBasic: uDW[15,12] - udw[15,12] >= 0
 uBasic: uDW[16,12] - udw[16,12] >= 0
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
 uBasic: uDX[13,13] - udx[13,13] >= 0
 uBasic: uDX[14,13] - udx[14,13] >= 0
 uBasic: uDX[15,13] - udx[15,13] >= 0
 uBasic: uDX[16,13] - udx[16,13] >= 0
 uBasic: uDX[17,13] - udx[17,13] >= 0
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
 uBasic: uDY[12,13] - udy[12,13] >= 0
 uBasic: uDZ[12,13] - udz[12,13] >= 0
 uBasic: uDY[13,13] - udy[13,13] >= 0
 uBasic: uDZ[13,13] - udz[13,13] >= 0
 uBasic: uDY[14,13] - udy[14,13] >= 0
 uBasic: uDZ[14,13] - udz[14,13] >= 0
 uBasic: uDY[15,13] - udy[15,13] >= 0
 uBasic: uDZ[15,13] - udz[15,13] >= 0
 uBasic: uDY[16,13] - udy[16,13] >= 0
 uBasic: uDZ[16,13] - udz[16,13] >= 0
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
 uBasic: uDW[12,13] - udw[12,13] >= 0
 uBasic: uDW[13,13] - udw[13,13] >= 0
 uBasic: uDW[14,13] - udw[14,13] >= 0
 uBasic: uDW[15,13] - udw[15,13] >= 0
 uBasic: uDW[16,13] - udw[16,13] >= 0
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
 uBasic: uDX[13,14] - udx[13,14] >= 0
 uBasic: uDX[14,14] - udx[14,14] >= 0
 uBasic: uDX[15,14] - udx[15,14] >= 0
 uBasic: uDX[16,14] - udx[16,14] >= 0
 uBasic: uDX[17,14] - udx[17,14] >= 0
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
 uBasic: uDY[12,14] - udy[12,14] >= 0
 uBasic: uDZ[12,14] - udz[12,14] >= 0
 uBasic: uDY[13,14] - udy[13,14] >= 0
 uBasic: uDZ[13,14] - udz[13,14] >= 0
 uBasic: uDY[14,14] - udy[14,14] >= 0
 uBasic: uDZ[14,14] - udz[14,14] >= 0
 uBasic: uDY[15,14] - udy[15,14] >= 0
 uBasic: uDZ[15,14] - udz[15,14] >= 0
 uBasic: uDY[16,14] - udy[16,14] >= 0
 uBasic: uDZ[16,14] - udz[16,14] >= 0
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
 uBasic: uDW[12,14] - udw[12,14] >= 0
 uBasic: uDW[13,14] - udw[13,14] >= 0
 uBasic: uDW[14,14] - udw[14,14] >= 0
 uBasic: uDW[15,14] - udw[15,14] >= 0
 uBasic: uDW[16,14] - udw[16,14] >= 0
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
 uBasic: uDX[13,15] - udx[13,15] >= 0
 uBasic: uDX[14,15] - udx[14,15] >= 0
 uBasic: uDX[15,15] - udx[15,15] >= 0
 uBasic: uDX[16,15] - udx[16,15] >= 0
 uBasic: uDX[17,15] - udx[17,15] >= 0
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
 uBasic: uDY[12,15] - udy[12,15] >= 0
 uBasic: uDZ[12,15] - udz[12,15] >= 0
 uBasic: uDY[13,15] - udy[13,15] >= 0
 uBasic: uDZ[13,15] - udz[13,15] >= 0
 uBasic: uDY[14,15] - udy[14,15] >= 0
 uBasic: uDZ[14,15] - udz[14,15] >= 0
 uBasic: uDY[15,15] - udy[15,15] >= 0
 uBasic: uDZ[15,15] - udz[15,15] >= 0
 uBasic: uDY[16,15] - udy[16,15] >= 0
 uBasic: uDZ[16,15] - udz[16,15] >= 0
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
 uBasic: uDW[12,15] - udw[12,15] >= 0
 uBasic: uDW[13,15] - udw[13,15] >= 0
 uBasic: uDW[14,15] - udw[14,15] >= 0
 uBasic: uDW[15,15] - udw[15,15] >= 0
 uBasic: uDW[16,15] - udw[16,15] >= 0
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
 ueSC: uDX[5,0] - uDY[5,0] = 0
 ueSC: udx[5,0] - uDY[5,0] = 0
 ueSC: uDX[5,1] - uDY[5,1] = 0
 ueSC: udx[5,1] - uDY[5,1] = 0
 ueSC: uDX[5,2] - uDY[5,2] = 0
 ueSC: udx[5,2] - uDY[5,2] = 0
 ueSC: uDX[5,3] - uDY[5,3] = 0
 ueSC: udx[5,3] - uDY[5,3] = 0
 ueSC: uDX[5,4] - uDY[5,4] = 0
 ueSC: udx[5,4] - uDY[5,4] = 0
 ueSC: uDX[5,5] - uDY[5,5] = 0
 ueSC: udx[5,5] - uDY[5,5] = 0
 ueSC: uDX[5,6] - uDY[5,6] = 0
 ueSC: udx[5,6] - uDY[5,6] = 0
 ueSC: uDX[5,7] - uDY[5,7] = 0
 ueSC: udx[5,7] - uDY[5,7] = 0
 ueSC: uDX[5,8] - uDY[5,8] = 0
 ueSC: udx[5,8] - uDY[5,8] = 0
 ueSC: uDX[5,9] - uDY[5,9] = 0
 ueSC: udx[5,9] - uDY[5,9] = 0
 ueSC: uDX[5,10] - uDY[5,10] = 0
 ueSC: udx[5,10] - uDY[5,10] = 0
 ueSC: uDX[5,11] - uDY[5,11] = 0
 ueSC: udx[5,11] - uDY[5,11] = 0
 ueSC: uDX[5,12] - uDY[5,12] = 0
 ueSC: udx[5,12] - uDY[5,12] = 0
 ueSC: uDX[5,13] - uDY[5,13] = 0
 ueSC: udx[5,13] - uDY[5,13] = 0
 ueSC: uDX[5,14] - uDY[5,14] = 0
 ueSC: udx[5,14] - uDY[5,14] = 0
 ueSC: uDX[5,15] - uDY[5,15] = 0
 ueSC: udx[5,15] - uDY[5,15] = 0
 ueSC: uDX[6,0] - uDY[6,0] = 0
 ueSC: udx[6,0] - uDY[6,0] = 0
 ueSC: uDX[6,1] - uDY[6,1] = 0
 ueSC: udx[6,1] - uDY[6,1] = 0
 ueSC: uDX[6,2] - uDY[6,2] = 0
 ueSC: udx[6,2] - uDY[6,2] = 0
 ueSC: uDX[6,3] - uDY[6,3] = 0
 ueSC: udx[6,3] - uDY[6,3] = 0
 ueSC: uDX[6,4] - uDY[6,4] = 0
 ueSC: udx[6,4] - uDY[6,4] = 0
 ueSC: uDX[6,5] - uDY[6,5] = 0
 ueSC: udx[6,5] - uDY[6,5] = 0
 ueSC: uDX[6,6] - uDY[6,6] = 0
 ueSC: udx[6,6] - uDY[6,6] = 0
 ueSC: uDX[6,7] - uDY[6,7] = 0
 ueSC: udx[6,7] - uDY[6,7] = 0
 ueSC: uDX[6,8] - uDY[6,8] = 0
 ueSC: udx[6,8] - uDY[6,8] = 0
 ueSC: uDX[6,9] - uDY[6,9] = 0
 ueSC: udx[6,9] - uDY[6,9] = 0
 ueSC: uDX[6,10] - uDY[6,10] = 0
 ueSC: udx[6,10] - uDY[6,10] = 0
 ueSC: uDX[6,11] - uDY[6,11] = 0
 ueSC: udx[6,11] - uDY[6,11] = 0
 ueSC: uDX[6,12] - uDY[6,12] = 0
 ueSC: udx[6,12] - uDY[6,12] = 0
 ueSC: uDX[6,13] - uDY[6,13] = 0
 ueSC: udx[6,13] - uDY[6,13] = 0
 ueSC: uDX[6,14] - uDY[6,14] = 0
 ueSC: udx[6,14] - uDY[6,14] = 0
 ueSC: uDX[6,15] - uDY[6,15] = 0
 ueSC: udx[6,15] - uDY[6,15] = 0
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
 udSC: - uDX[12,0] + uDY[12,0] = 0
 udSC: - uDX[12,0] + udy[12,0] = 0
 udSC: - uDX[12,1] + uDY[12,1] = 0
 udSC: - uDX[12,1] + udy[12,1] = 0
 udSC: - uDX[12,2] + uDY[12,2] = 0
 udSC: - uDX[12,2] + udy[12,2] = 0
 udSC: - uDX[12,3] + uDY[12,3] = 0
 udSC: - uDX[12,3] + udy[12,3] = 0
 udSC: - uDX[12,4] + uDY[12,4] = 0
 udSC: - uDX[12,4] + udy[12,4] = 0
 udSC: - uDX[12,5] + uDY[12,5] = 0
 udSC: - uDX[12,5] + udy[12,5] = 0
 udSC: - uDX[12,6] + uDY[12,6] = 0
 udSC: - uDX[12,6] + udy[12,6] = 0
 udSC: - uDX[12,7] + uDY[12,7] = 0
 udSC: - uDX[12,7] + udy[12,7] = 0
 udSC: - uDX[12,8] + uDY[12,8] = 0
 udSC: - uDX[12,8] + udy[12,8] = 0
 udSC: - uDX[12,9] + uDY[12,9] = 0
 udSC: - uDX[12,9] + udy[12,9] = 0
 udSC: - uDX[12,10] + uDY[12,10] = 0
 udSC: - uDX[12,10] + udy[12,10] = 0
 udSC: - uDX[12,11] + uDY[12,11] = 0
 udSC: - uDX[12,11] + udy[12,11] = 0
 udSC: - uDX[12,12] + uDY[12,12] = 0
 udSC: - uDX[12,12] + udy[12,12] = 0
 udSC: - uDX[12,13] + uDY[12,13] = 0
 udSC: - uDX[12,13] + udy[12,13] = 0
 udSC: - uDX[12,14] + uDY[12,14] = 0
 udSC: - uDX[12,14] + udy[12,14] = 0
 udSC: - uDX[12,15] + uDY[12,15] = 0
 udSC: - uDX[12,15] + udy[12,15] = 0
 udSC: - uDX[13,0] + uDY[13,0] = 0
 udSC: - uDX[13,0] + udy[13,0] = 0
 udSC: - uDX[13,1] + uDY[13,1] = 0
 udSC: - uDX[13,1] + udy[13,1] = 0
 udSC: - uDX[13,2] + uDY[13,2] = 0
 udSC: - uDX[13,2] + udy[13,2] = 0
 udSC: - uDX[13,3] + uDY[13,3] = 0
 udSC: - uDX[13,3] + udy[13,3] = 0
 udSC: - uDX[13,4] + uDY[13,4] = 0
 udSC: - uDX[13,4] + udy[13,4] = 0
 udSC: - uDX[13,5] + uDY[13,5] = 0
 udSC: - uDX[13,5] + udy[13,5] = 0
 udSC: - uDX[13,6] + uDY[13,6] = 0
 udSC: - uDX[13,6] + udy[13,6] = 0
 udSC: - uDX[13,7] + uDY[13,7] = 0
 udSC: - uDX[13,7] + udy[13,7] = 0
 udSC: - uDX[13,8] + uDY[13,8] = 0
 udSC: - uDX[13,8] + udy[13,8] = 0
 udSC: - uDX[13,9] + uDY[13,9] = 0
 udSC: - uDX[13,9] + udy[13,9] = 0
 udSC: - uDX[13,10] + uDY[13,10] = 0
 udSC: - uDX[13,10] + udy[13,10] = 0
 udSC: - uDX[13,11] + uDY[13,11] = 0
 udSC: - uDX[13,11] + udy[13,11] = 0
 udSC: - uDX[13,12] + uDY[13,12] = 0
 udSC: - uDX[13,12] + udy[13,12] = 0
 udSC: - uDX[13,13] + uDY[13,13] = 0
 udSC: - uDX[13,13] + udy[13,13] = 0
 udSC: - uDX[13,14] + uDY[13,14] = 0
 udSC: - uDX[13,14] + udy[13,14] = 0
 udSC: - uDX[13,15] + uDY[13,15] = 0
 udSC: - uDX[13,15] + udy[13,15] = 0
 udSC: - uDX[14,0] + uDY[14,0] = 0
 udSC: - uDX[14,0] + udy[14,0] = 0
 udSC: - uDX[14,1] + uDY[14,1] = 0
 udSC: - uDX[14,1] + udy[14,1] = 0
 udSC: - uDX[14,2] + uDY[14,2] = 0
 udSC: - uDX[14,2] + udy[14,2] = 0
 udSC: - uDX[14,3] + uDY[14,3] = 0
 udSC: - uDX[14,3] + udy[14,3] = 0
 udSC: - uDX[14,4] + uDY[14,4] = 0
 udSC: - uDX[14,4] + udy[14,4] = 0
 udSC: - uDX[14,5] + uDY[14,5] = 0
 udSC: - uDX[14,5] + udy[14,5] = 0
 udSC: - uDX[14,6] + uDY[14,6] = 0
 udSC: - uDX[14,6] + udy[14,6] = 0
 udSC: - uDX[14,7] + uDY[14,7] = 0
 udSC: - uDX[14,7] + udy[14,7] = 0
 udSC: - uDX[14,8] + uDY[14,8] = 0
 udSC: - uDX[14,8] + udy[14,8] = 0
 udSC: - uDX[14,9] + uDY[14,9] = 0
 udSC: - uDX[14,9] + udy[14,9] = 0
 udSC: - uDX[14,10] + uDY[14,10] = 0
 udSC: - uDX[14,10] + udy[14,10] = 0
 udSC: - uDX[14,11] + uDY[14,11] = 0
 udSC: - uDX[14,11] + udy[14,11] = 0
 udSC: - uDX[14,12] + uDY[14,12] = 0
 udSC: - uDX[14,12] + udy[14,12] = 0
 udSC: - uDX[14,13] + uDY[14,13] = 0
 udSC: - uDX[14,13] + udy[14,13] = 0
 udSC: - uDX[14,14] + uDY[14,14] = 0
 udSC: - uDX[14,14] + udy[14,14] = 0
 udSC: - uDX[14,15] + uDY[14,15] = 0
 udSC: - uDX[14,15] + udy[14,15] = 0
 udSC: - uDX[15,0] + uDY[15,0] = 0
 udSC: - uDX[15,0] + udy[15,0] = 0
 udSC: - uDX[15,1] + uDY[15,1] = 0
 udSC: - uDX[15,1] + udy[15,1] = 0
 udSC: - uDX[15,2] + uDY[15,2] = 0
 udSC: - uDX[15,2] + udy[15,2] = 0
 udSC: - uDX[15,3] + uDY[15,3] = 0
 udSC: - uDX[15,3] + udy[15,3] = 0
 udSC: - uDX[15,4] + uDY[15,4] = 0
 udSC: - uDX[15,4] + udy[15,4] = 0
 udSC: - uDX[15,5] + uDY[15,5] = 0
 udSC: - uDX[15,5] + udy[15,5] = 0
 udSC: - uDX[15,6] + uDY[15,6] = 0
 udSC: - uDX[15,6] + udy[15,6] = 0
 udSC: - uDX[15,7] + uDY[15,7] = 0
 udSC: - uDX[15,7] + udy[15,7] = 0
 udSC: - uDX[15,8] + uDY[15,8] = 0
 udSC: - uDX[15,8] + udy[15,8] = 0
 udSC: - uDX[15,9] + uDY[15,9] = 0
 udSC: - uDX[15,9] + udy[15,9] = 0
 udSC: - uDX[15,10] + uDY[15,10] = 0
 udSC: - uDX[15,10] + udy[15,10] = 0
 udSC: - uDX[15,11] + uDY[15,11] = 0
 udSC: - uDX[15,11] + udy[15,11] = 0
 udSC: - uDX[15,12] + uDY[15,12] = 0
 udSC: - uDX[15,12] + udy[15,12] = 0
 udSC: - uDX[15,13] + uDY[15,13] = 0
 udSC: - uDX[15,13] + udy[15,13] = 0
 udSC: - uDX[15,14] + uDY[15,14] = 0
 udSC: - uDX[15,14] + udy[15,14] = 0
 udSC: - uDX[15,15] + uDY[15,15] = 0
 udSC: - uDX[15,15] + udy[15,15] = 0
 udSC: - uDX[16,0] + uDY[16,0] = 0
 udSC: - uDX[16,0] + udy[16,0] = 0
 udSC: - uDX[16,1] + uDY[16,1] = 0
 udSC: - uDX[16,1] + udy[16,1] = 0
 udSC: - uDX[16,2] + uDY[16,2] = 0
 udSC: - uDX[16,2] + udy[16,2] = 0
 udSC: - uDX[16,3] + uDY[16,3] = 0
 udSC: - uDX[16,3] + udy[16,3] = 0
 udSC: - uDX[16,4] + uDY[16,4] = 0
 udSC: - uDX[16,4] + udy[16,4] = 0
 udSC: - uDX[16,5] + uDY[16,5] = 0
 udSC: - uDX[16,5] + udy[16,5] = 0
 udSC: - uDX[16,6] + uDY[16,6] = 0
 udSC: - uDX[16,6] + udy[16,6] = 0
 udSC: - uDX[16,7] + uDY[16,7] = 0
 udSC: - uDX[16,7] + udy[16,7] = 0
 udSC: - uDX[16,8] + uDY[16,8] = 0
 udSC: - uDX[16,8] + udy[16,8] = 0
 udSC: - uDX[16,9] + uDY[16,9] = 0
 udSC: - uDX[16,9] + udy[16,9] = 0
 udSC: - uDX[16,10] + uDY[16,10] = 0
 udSC: - uDX[16,10] + udy[16,10] = 0
 udSC: - uDX[16,11] + uDY[16,11] = 0
 udSC: - uDX[16,11] + udy[16,11] = 0
 udSC: - uDX[16,12] + uDY[16,12] = 0
 udSC: - uDX[16,12] + udy[16,12] = 0
 udSC: - uDX[16,13] + uDY[16,13] = 0
 udSC: - uDX[16,13] + udy[16,13] = 0
 udSC: - uDX[16,14] + uDY[16,14] = 0
 udSC: - uDX[16,14] + udy[16,14] = 0
 udSC: - uDX[16,15] + uDY[16,15] = 0
 udSC: - uDX[16,15] + udy[16,15] = 0
 eqk0: - ustk[0,0] + eqk[0] = 0
 eqk1: - ustk[0,0] + eqk[4] = 0
 eqk2: - ustk[0,7] + eqk[8] = 0
 eqk2: - ustk[0,0] + eqk[12] = 0
 eqk0: - ustk[0,1] + eqk[1] = 0
 eqk1: - ustk[0,1] + eqk[5] = 0
 eqk2: - ustk[0,4] + eqk[9] = 0
 eqk2: - ustk[0,1] + eqk[13] = 0
 eqk0: - ustk[0,2] + eqk[2] = 0
 eqk1: - ustk[0,2] + eqk[6] = 0
 eqk2: - ustk[0,5] + eqk[10] = 0
 eqk2: - ustk[0,2] + eqk[14] = 0
 eqk0: - ustk[0,3] + eqk[3] = 0
 eqk1: - ustk[0,3] + eqk[7] = 0
 eqk2: - ustk[0,6] + eqk[11] = 0
 eqk2: - ustk[0,3] + eqk[15] = 0
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
 uATKri: udy[6,0] - udz[6,0] = 0
 uATKri: uDY[6,0] - uDZ[6,0] + ustk[6,0] >= 0
 uATKri: - uDY[6,0] + uDZ[6,0] + ustk[6,0] >= 0
 uATKri: uDY[6,0] + uDZ[6,0] - udz[6,0] - ustk[6,0] >= 0
 uATKri: uDY[6,0] - udz[6,0] >= 0
 uATKri: udy[6,1] - udz[6,1] = 0
 uATKri: uDY[6,1] - uDZ[6,1] + ustk[6,1] >= 0
 uATKri: - uDY[6,1] + uDZ[6,1] + ustk[6,1] >= 0
 uATKri: uDY[6,1] + uDZ[6,1] - udz[6,1] - ustk[6,1] >= 0
 uATKri: uDY[6,1] - udz[6,1] >= 0
 uATKri: udy[6,2] - udz[6,2] = 0
 uATKri: uDY[6,2] - uDZ[6,2] + ustk[6,2] >= 0
 uATKri: - uDY[6,2] + uDZ[6,2] + ustk[6,2] >= 0
 uATKri: uDY[6,2] + uDZ[6,2] - udz[6,2] - ustk[6,2] >= 0
 uATKri: uDY[6,2] - udz[6,2] >= 0
 uATKri: udy[6,3] - udz[6,3] = 0
 uATKri: uDY[6,3] - uDZ[6,3] + ustk[6,3] >= 0
 uATKri: - uDY[6,3] + uDZ[6,3] + ustk[6,3] >= 0
 uATKri: uDY[6,3] + uDZ[6,3] - udz[6,3] - ustk[6,3] >= 0
 uATKri: uDY[6,3] - udz[6,3] >= 0
 uATKri: udy[6,4] - udz[6,4] = 0
 uATKri: uDY[6,4] - uDZ[6,4] + ustk[6,4] >= 0
 uATKri: - uDY[6,4] + uDZ[6,4] + ustk[6,4] >= 0
 uATKri: uDY[6,4] + uDZ[6,4] - udz[6,4] - ustk[6,4] >= 0
 uATKri: uDY[6,4] - udz[6,4] >= 0
 uATKri: udy[6,5] - udz[6,5] = 0
 uATKri: uDY[6,5] - uDZ[6,5] + ustk[6,5] >= 0
 uATKri: - uDY[6,5] + uDZ[6,5] + ustk[6,5] >= 0
 uATKri: uDY[6,5] + uDZ[6,5] - udz[6,5] - ustk[6,5] >= 0
 uATKri: uDY[6,5] - udz[6,5] >= 0
 uATKri: udy[6,6] - udz[6,6] = 0
 uATKri: uDY[6,6] - uDZ[6,6] + ustk[6,6] >= 0
 uATKri: - uDY[6,6] + uDZ[6,6] + ustk[6,6] >= 0
 uATKri: uDY[6,6] + uDZ[6,6] - udz[6,6] - ustk[6,6] >= 0
 uATKri: uDY[6,6] - udz[6,6] >= 0
 uATKri: udy[6,7] - udz[6,7] = 0
 uATKri: uDY[6,7] - uDZ[6,7] + ustk[6,7] >= 0
 uATKri: - uDY[6,7] + uDZ[6,7] + ustk[6,7] >= 0
 uATKri: uDY[6,7] + uDZ[6,7] - udz[6,7] - ustk[6,7] >= 0
 uATKri: uDY[6,7] - udz[6,7] >= 0
 uATKri: uDY[6,8] - uDZ[6,8] = 0
 uATKri: udy[6,8] - udz[6,8] = 0
 uATKri: uDY[6,9] - uDZ[6,9] = 0
 uATKri: udy[6,9] - udz[6,9] = 0
 uATKri: uDY[6,10] - uDZ[6,10] = 0
 uATKri: udy[6,10] - udz[6,10] = 0
 uATKri: uDY[6,11] - uDZ[6,11] = 0
 uATKri: udy[6,11] - udz[6,11] = 0
 uATKri: uDY[6,12] - uDZ[6,12] = 0
 uATKri: udy[6,12] - udz[6,12] = 0
 uATKri: uDY[6,13] - uDZ[6,13] = 0
 uATKri: udy[6,13] - udz[6,13] = 0
 uATKri: uDY[6,14] - uDZ[6,14] = 0
 uATKri: udy[6,14] - udz[6,14] = 0
 uATKri: uDY[6,15] - uDZ[6,15] = 0
 uATKri: udy[6,15] - udz[6,15] = 0
 uATKr1: udy[7,0] - udz[7,0] = 0
 uATKr1: uDY[7,0] - uDZ[7,0] + ustk[7,0] >= 0
 uATKr1: - uDY[7,0] + uDZ[7,0] + ustk[7,0] >= 0
 uATKr1: uDY[7,0] + uDZ[7,0] - udz[7,0] - ustk[7,0] >= 0
 uATKr1: uDY[7,0] - udz[7,0] >= 0
 uATKr1: udy[7,1] - udz[7,1] = 0
 uATKr1: uDY[7,1] - uDZ[7,1] + ustk[7,1] >= 0
 uATKr1: - uDY[7,1] + uDZ[7,1] + ustk[7,1] >= 0
 uATKr1: uDY[7,1] + uDZ[7,1] - udz[7,1] - ustk[7,1] >= 0
 uATKr1: uDY[7,1] - udz[7,1] >= 0
 uATKr1: udy[7,2] - udz[7,2] = 0
 uATKr1: uDY[7,2] - uDZ[7,2] + ustk[7,2] >= 0
 uATKr1: - uDY[7,2] + uDZ[7,2] + ustk[7,2] >= 0
 uATKr1: uDY[7,2] + uDZ[7,2] - udz[7,2] - ustk[7,2] >= 0
 uATKr1: uDY[7,2] - udz[7,2] >= 0
 uATKr1: udy[7,3] - udz[7,3] = 0
 uATKr1: uDY[7,3] - uDZ[7,3] + ustk[7,3] >= 0
 uATKr1: - uDY[7,3] + uDZ[7,3] + ustk[7,3] >= 0
 uATKr1: uDY[7,3] + uDZ[7,3] - udz[7,3] - ustk[7,3] >= 0
 uATKr1: uDY[7,3] - udz[7,3] >= 0
 uATKr1: udy[7,4] - udz[7,4] = 0
 uATKr1: uDY[7,4] - uDZ[7,4] + ustk[7,4] >= 0
 uATKr1: - uDY[7,4] + uDZ[7,4] + ustk[7,4] >= 0
 uATKr1: uDY[7,4] + uDZ[7,4] - udz[7,4] - ustk[7,4] >= 0
 uATKr1: uDY[7,4] - udz[7,4] >= 0
 uATKr1: udy[7,5] - udz[7,5] = 0
 uATKr1: uDY[7,5] - uDZ[7,5] + ustk[7,5] >= 0
 uATKr1: - uDY[7,5] + uDZ[7,5] + ustk[7,5] >= 0
 uATKr1: uDY[7,5] + uDZ[7,5] - udz[7,5] - ustk[7,5] >= 0
 uATKr1: uDY[7,5] - udz[7,5] >= 0
 uATKr1: udy[7,6] - udz[7,6] = 0
 uATKr1: uDY[7,6] - uDZ[7,6] + ustk[7,6] >= 0
 uATKr1: - uDY[7,6] + uDZ[7,6] + ustk[7,6] >= 0
 uATKr1: uDY[7,6] + uDZ[7,6] - udz[7,6] - ustk[7,6] >= 0
 uATKr1: uDY[7,6] - udz[7,6] >= 0
 uATKr1: udy[7,7] - udz[7,7] = 0
 uATKr1: uDY[7,7] - uDZ[7,7] + ustk[7,7] >= 0
 uATKr1: - uDY[7,7] + uDZ[7,7] + ustk[7,7] >= 0
 uATKr1: uDY[7,7] + uDZ[7,7] - udz[7,7] - ustk[7,7] >= 0
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
 uATKr1: uDY[8,0] - uDZ[8,0] + ustk[8,0] >= 0
 uATKr1: - uDY[8,0] + uDZ[8,0] + ustk[8,0] >= 0
 uATKr1: uDY[8,0] + uDZ[8,0] - udz[8,0] - ustk[8,0] >= 0
 uATKr1: uDY[8,0] - udz[8,0] >= 0
 uATKr1: udy[8,1] - udz[8,1] = 0
 uATKr1: uDY[8,1] - uDZ[8,1] + ustk[8,1] >= 0
 uATKr1: - uDY[8,1] + uDZ[8,1] + ustk[8,1] >= 0
 uATKr1: uDY[8,1] + uDZ[8,1] - udz[8,1] - ustk[8,1] >= 0
 uATKr1: uDY[8,1] - udz[8,1] >= 0
 uATKr1: udy[8,2] - udz[8,2] = 0
 uATKr1: uDY[8,2] - uDZ[8,2] + ustk[8,2] >= 0
 uATKr1: - uDY[8,2] + uDZ[8,2] + ustk[8,2] >= 0
 uATKr1: uDY[8,2] + uDZ[8,2] - udz[8,2] - ustk[8,2] >= 0
 uATKr1: uDY[8,2] - udz[8,2] >= 0
 uATKr1: udy[8,3] - udz[8,3] = 0
 uATKr1: uDY[8,3] - uDZ[8,3] + ustk[8,3] >= 0
 uATKr1: - uDY[8,3] + uDZ[8,3] + ustk[8,3] >= 0
 uATKr1: uDY[8,3] + uDZ[8,3] - udz[8,3] - ustk[8,3] >= 0
 uATKr1: uDY[8,3] - udz[8,3] >= 0
 uATKr1: udy[8,4] - udz[8,4] = 0
 uATKr1: uDY[8,4] - uDZ[8,4] + ustk[8,4] >= 0
 uATKr1: - uDY[8,4] + uDZ[8,4] + ustk[8,4] >= 0
 uATKr1: uDY[8,4] + uDZ[8,4] - udz[8,4] - ustk[8,4] >= 0
 uATKr1: uDY[8,4] - udz[8,4] >= 0
 uATKr1: udy[8,5] - udz[8,5] = 0
 uATKr1: uDY[8,5] - uDZ[8,5] + ustk[8,5] >= 0
 uATKr1: - uDY[8,5] + uDZ[8,5] + ustk[8,5] >= 0
 uATKr1: uDY[8,5] + uDZ[8,5] - udz[8,5] - ustk[8,5] >= 0
 uATKr1: uDY[8,5] - udz[8,5] >= 0
 uATKr1: udy[8,6] - udz[8,6] = 0
 uATKr1: uDY[8,6] - uDZ[8,6] + ustk[8,6] >= 0
 uATKr1: - uDY[8,6] + uDZ[8,6] + ustk[8,6] >= 0
 uATKr1: uDY[8,6] + uDZ[8,6] - udz[8,6] - ustk[8,6] >= 0
 uATKr1: uDY[8,6] - udz[8,6] >= 0
 uATKr1: udy[8,7] - udz[8,7] = 0
 uATKr1: uDY[8,7] - uDZ[8,7] + ustk[8,7] >= 0
 uATKr1: - uDY[8,7] + uDZ[8,7] + ustk[8,7] >= 0
 uATKr1: uDY[8,7] + uDZ[8,7] - udz[8,7] - ustk[8,7] >= 0
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
 uATKr1: uDY[9,0] - uDZ[9,0] + ustk[9,0] >= 0
 uATKr1: - uDY[9,0] + uDZ[9,0] + ustk[9,0] >= 0
 uATKr1: uDY[9,0] + uDZ[9,0] - udz[9,0] - ustk[9,0] >= 0
 uATKr1: uDY[9,0] - udz[9,0] >= 0
 uATKr1: udy[9,1] - udz[9,1] = 0
 uATKr1: uDY[9,1] - uDZ[9,1] + ustk[9,1] >= 0
 uATKr1: - uDY[9,1] + uDZ[9,1] + ustk[9,1] >= 0
 uATKr1: uDY[9,1] + uDZ[9,1] - udz[9,1] - ustk[9,1] >= 0
 uATKr1: uDY[9,1] - udz[9,1] >= 0
 uATKr1: udy[9,2] - udz[9,2] = 0
 uATKr1: uDY[9,2] - uDZ[9,2] + ustk[9,2] >= 0
 uATKr1: - uDY[9,2] + uDZ[9,2] + ustk[9,2] >= 0
 uATKr1: uDY[9,2] + uDZ[9,2] - udz[9,2] - ustk[9,2] >= 0
 uATKr1: uDY[9,2] - udz[9,2] >= 0
 uATKr1: udy[9,3] - udz[9,3] = 0
 uATKr1: uDY[9,3] - uDZ[9,3] + ustk[9,3] >= 0
 uATKr1: - uDY[9,3] + uDZ[9,3] + ustk[9,3] >= 0
 uATKr1: uDY[9,3] + uDZ[9,3] - udz[9,3] - ustk[9,3] >= 0
 uATKr1: uDY[9,3] - udz[9,3] >= 0
 uATKr1: udy[9,4] - udz[9,4] = 0
 uATKr1: uDY[9,4] - uDZ[9,4] + ustk[9,4] >= 0
 uATKr1: - uDY[9,4] + uDZ[9,4] + ustk[9,4] >= 0
 uATKr1: uDY[9,4] + uDZ[9,4] - udz[9,4] - ustk[9,4] >= 0
 uATKr1: uDY[9,4] - udz[9,4] >= 0
 uATKr1: udy[9,5] - udz[9,5] = 0
 uATKr1: uDY[9,5] - uDZ[9,5] + ustk[9,5] >= 0
 uATKr1: - uDY[9,5] + uDZ[9,5] + ustk[9,5] >= 0
 uATKr1: uDY[9,5] + uDZ[9,5] - udz[9,5] - ustk[9,5] >= 0
 uATKr1: uDY[9,5] - udz[9,5] >= 0
 uATKr1: udy[9,6] - udz[9,6] = 0
 uATKr1: uDY[9,6] - uDZ[9,6] + ustk[9,6] >= 0
 uATKr1: - uDY[9,6] + uDZ[9,6] + ustk[9,6] >= 0
 uATKr1: uDY[9,6] + uDZ[9,6] - udz[9,6] - ustk[9,6] >= 0
 uATKr1: uDY[9,6] - udz[9,6] >= 0
 uATKr1: udy[9,7] - udz[9,7] = 0
 uATKr1: uDY[9,7] - uDZ[9,7] + ustk[9,7] >= 0
 uATKr1: - uDY[9,7] + uDZ[9,7] + ustk[9,7] >= 0
 uATKr1: uDY[9,7] + uDZ[9,7] - udz[9,7] - ustk[9,7] >= 0
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
 uATKr1: uDY[10,0] - uDZ[10,0] + ustk[10,0] >= 0
 uATKr1: - uDY[10,0] + uDZ[10,0] + ustk[10,0] >= 0
 uATKr1: uDY[10,0] + uDZ[10,0] - udz[10,0] - ustk[10,0] >= 0
 uATKr1: uDY[10,0] - udz[10,0] >= 0
 uATKr1: udy[10,1] - udz[10,1] = 0
 uATKr1: uDY[10,1] - uDZ[10,1] + ustk[10,1] >= 0
 uATKr1: - uDY[10,1] + uDZ[10,1] + ustk[10,1] >= 0
 uATKr1: uDY[10,1] + uDZ[10,1] - udz[10,1] - ustk[10,1] >= 0
 uATKr1: uDY[10,1] - udz[10,1] >= 0
 uATKr1: udy[10,2] - udz[10,2] = 0
 uATKr1: uDY[10,2] - uDZ[10,2] + ustk[10,2] >= 0
 uATKr1: - uDY[10,2] + uDZ[10,2] + ustk[10,2] >= 0
 uATKr1: uDY[10,2] + uDZ[10,2] - udz[10,2] - ustk[10,2] >= 0
 uATKr1: uDY[10,2] - udz[10,2] >= 0
 uATKr1: udy[10,3] - udz[10,3] = 0
 uATKr1: uDY[10,3] - uDZ[10,3] + ustk[10,3] >= 0
 uATKr1: - uDY[10,3] + uDZ[10,3] + ustk[10,3] >= 0
 uATKr1: uDY[10,3] + uDZ[10,3] - udz[10,3] - ustk[10,3] >= 0
 uATKr1: uDY[10,3] - udz[10,3] >= 0
 uATKr1: udy[10,4] - udz[10,4] = 0
 uATKr1: uDY[10,4] - uDZ[10,4] + ustk[10,4] >= 0
 uATKr1: - uDY[10,4] + uDZ[10,4] + ustk[10,4] >= 0
 uATKr1: uDY[10,4] + uDZ[10,4] - udz[10,4] - ustk[10,4] >= 0
 uATKr1: uDY[10,4] - udz[10,4] >= 0
 uATKr1: udy[10,5] - udz[10,5] = 0
 uATKr1: uDY[10,5] - uDZ[10,5] + ustk[10,5] >= 0
 uATKr1: - uDY[10,5] + uDZ[10,5] + ustk[10,5] >= 0
 uATKr1: uDY[10,5] + uDZ[10,5] - udz[10,5] - ustk[10,5] >= 0
 uATKr1: uDY[10,5] - udz[10,5] >= 0
 uATKr1: udy[10,6] - udz[10,6] = 0
 uATKr1: uDY[10,6] - uDZ[10,6] + ustk[10,6] >= 0
 uATKr1: - uDY[10,6] + uDZ[10,6] + ustk[10,6] >= 0
 uATKr1: uDY[10,6] + uDZ[10,6] - udz[10,6] - ustk[10,6] >= 0
 uATKr1: uDY[10,6] - udz[10,6] >= 0
 uATKr1: udy[10,7] - udz[10,7] = 0
 uATKr1: uDY[10,7] - uDZ[10,7] + ustk[10,7] >= 0
 uATKr1: - uDY[10,7] + uDZ[10,7] + ustk[10,7] >= 0
 uATKr1: uDY[10,7] + uDZ[10,7] - udz[10,7] - ustk[10,7] >= 0
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
 uATKr1: uDY[11,0] - uDZ[11,0] + ustk[11,0] >= 0
 uATKr1: - uDY[11,0] + uDZ[11,0] + ustk[11,0] >= 0
 uATKr1: uDY[11,0] + uDZ[11,0] - udz[11,0] - ustk[11,0] >= 0
 uATKr1: uDY[11,0] - udz[11,0] >= 0
 uATKr1: udy[11,1] - udz[11,1] = 0
 uATKr1: uDY[11,1] - uDZ[11,1] + ustk[11,1] >= 0
 uATKr1: - uDY[11,1] + uDZ[11,1] + ustk[11,1] >= 0
 uATKr1: uDY[11,1] + uDZ[11,1] - udz[11,1] - ustk[11,1] >= 0
 uATKr1: uDY[11,1] - udz[11,1] >= 0
 uATKr1: udy[11,2] - udz[11,2] = 0
 uATKr1: uDY[11,2] - uDZ[11,2] + ustk[11,2] >= 0
 uATKr1: - uDY[11,2] + uDZ[11,2] + ustk[11,2] >= 0
 uATKr1: uDY[11,2] + uDZ[11,2] - udz[11,2] - ustk[11,2] >= 0
 uATKr1: uDY[11,2] - udz[11,2] >= 0
 uATKr1: udy[11,3] - udz[11,3] = 0
 uATKr1: uDY[11,3] - uDZ[11,3] + ustk[11,3] >= 0
 uATKr1: - uDY[11,3] + uDZ[11,3] + ustk[11,3] >= 0
 uATKr1: uDY[11,3] + uDZ[11,3] - udz[11,3] - ustk[11,3] >= 0
 uATKr1: uDY[11,3] - udz[11,3] >= 0
 uATKr1: udy[11,4] - udz[11,4] = 0
 uATKr1: uDY[11,4] - uDZ[11,4] + ustk[11,4] >= 0
 uATKr1: - uDY[11,4] + uDZ[11,4] + ustk[11,4] >= 0
 uATKr1: uDY[11,4] + uDZ[11,4] - udz[11,4] - ustk[11,4] >= 0
 uATKr1: uDY[11,4] - udz[11,4] >= 0
 uATKr1: udy[11,5] - udz[11,5] = 0
 uATKr1: uDY[11,5] - uDZ[11,5] + ustk[11,5] >= 0
 uATKr1: - uDY[11,5] + uDZ[11,5] + ustk[11,5] >= 0
 uATKr1: uDY[11,5] + uDZ[11,5] - udz[11,5] - ustk[11,5] >= 0
 uATKr1: uDY[11,5] - udz[11,5] >= 0
 uATKr1: udy[11,6] - udz[11,6] = 0
 uATKr1: uDY[11,6] - uDZ[11,6] + ustk[11,6] >= 0
 uATKr1: - uDY[11,6] + uDZ[11,6] + ustk[11,6] >= 0
 uATKr1: uDY[11,6] + uDZ[11,6] - udz[11,6] - ustk[11,6] >= 0
 uATKr1: uDY[11,6] - udz[11,6] >= 0
 uATKr1: udy[11,7] - udz[11,7] = 0
 uATKr1: uDY[11,7] - uDZ[11,7] + ustk[11,7] >= 0
 uATKr1: - uDY[11,7] + uDZ[11,7] + ustk[11,7] >= 0
 uATKr1: uDY[11,7] + uDZ[11,7] - udz[11,7] - ustk[11,7] >= 0
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
 uATKr1: udy[12,0] - udz[12,0] = 0
 uATKr1: uDY[12,0] - uDZ[12,0] + ustk[12,0] >= 0
 uATKr1: - uDY[12,0] + uDZ[12,0] + ustk[12,0] >= 0
 uATKr1: uDY[12,0] + uDZ[12,0] - udz[12,0] - ustk[12,0] >= 0
 uATKr1: uDY[12,0] - udz[12,0] >= 0
 uATKr1: udy[12,1] - udz[12,1] = 0
 uATKr1: uDY[12,1] - uDZ[12,1] + ustk[12,1] >= 0
 uATKr1: - uDY[12,1] + uDZ[12,1] + ustk[12,1] >= 0
 uATKr1: uDY[12,1] + uDZ[12,1] - udz[12,1] - ustk[12,1] >= 0
 uATKr1: uDY[12,1] - udz[12,1] >= 0
 uATKr1: udy[12,2] - udz[12,2] = 0
 uATKr1: uDY[12,2] - uDZ[12,2] + ustk[12,2] >= 0
 uATKr1: - uDY[12,2] + uDZ[12,2] + ustk[12,2] >= 0
 uATKr1: uDY[12,2] + uDZ[12,2] - udz[12,2] - ustk[12,2] >= 0
 uATKr1: uDY[12,2] - udz[12,2] >= 0
 uATKr1: udy[12,3] - udz[12,3] = 0
 uATKr1: uDY[12,3] - uDZ[12,3] + ustk[12,3] >= 0
 uATKr1: - uDY[12,3] + uDZ[12,3] + ustk[12,3] >= 0
 uATKr1: uDY[12,3] + uDZ[12,3] - udz[12,3] - ustk[12,3] >= 0
 uATKr1: uDY[12,3] - udz[12,3] >= 0
 uATKr1: udy[12,4] - udz[12,4] = 0
 uATKr1: uDY[12,4] - uDZ[12,4] + ustk[12,4] >= 0
 uATKr1: - uDY[12,4] + uDZ[12,4] + ustk[12,4] >= 0
 uATKr1: uDY[12,4] + uDZ[12,4] - udz[12,4] - ustk[12,4] >= 0
 uATKr1: uDY[12,4] - udz[12,4] >= 0
 uATKr1: udy[12,5] - udz[12,5] = 0
 uATKr1: uDY[12,5] - uDZ[12,5] + ustk[12,5] >= 0
 uATKr1: - uDY[12,5] + uDZ[12,5] + ustk[12,5] >= 0
 uATKr1: uDY[12,5] + uDZ[12,5] - udz[12,5] - ustk[12,5] >= 0
 uATKr1: uDY[12,5] - udz[12,5] >= 0
 uATKr1: udy[12,6] - udz[12,6] = 0
 uATKr1: uDY[12,6] - uDZ[12,6] + ustk[12,6] >= 0
 uATKr1: - uDY[12,6] + uDZ[12,6] + ustk[12,6] >= 0
 uATKr1: uDY[12,6] + uDZ[12,6] - udz[12,6] - ustk[12,6] >= 0
 uATKr1: uDY[12,6] - udz[12,6] >= 0
 uATKr1: udy[12,7] - udz[12,7] = 0
 uATKr1: uDY[12,7] - uDZ[12,7] + ustk[12,7] >= 0
 uATKr1: - uDY[12,7] + uDZ[12,7] + ustk[12,7] >= 0
 uATKr1: uDY[12,7] + uDZ[12,7] - udz[12,7] - ustk[12,7] >= 0
 uATKr1: uDY[12,7] - udz[12,7] >= 0
 uATKr1: uDY[12,8] - uDZ[12,8] = 0
 uATKr1: udy[12,8] - udz[12,8] = 0
 uATKr1: uDY[12,9] - uDZ[12,9] = 0
 uATKr1: udy[12,9] - udz[12,9] = 0
 uATKr1: uDY[12,10] - uDZ[12,10] = 0
 uATKr1: udy[12,10] - udz[12,10] = 0
 uATKr1: uDY[12,11] - uDZ[12,11] = 0
 uATKr1: udy[12,11] - udz[12,11] = 0
 uATKr1: uDY[12,12] - uDZ[12,12] = 0
 uATKr1: udy[12,12] - udz[12,12] = 0
 uATKr1: uDY[12,13] - uDZ[12,13] = 0
 uATKr1: udy[12,13] - udz[12,13] = 0
 uATKr1: uDY[12,14] - uDZ[12,14] = 0
 uATKr1: udy[12,14] - udz[12,14] = 0
 uATKr1: uDY[12,15] - uDZ[12,15] = 0
 uATKr1: udy[12,15] - udz[12,15] = 0
 uATKr1: udy[13,0] - udz[13,0] = 0
 uATKr1: uDY[13,0] - uDZ[13,0] + ustk[13,0] >= 0
 uATKr1: - uDY[13,0] + uDZ[13,0] + ustk[13,0] >= 0
 uATKr1: uDY[13,0] + uDZ[13,0] - udz[13,0] - ustk[13,0] >= 0
 uATKr1: uDY[13,0] - udz[13,0] >= 0
 uATKr1: udy[13,1] - udz[13,1] = 0
 uATKr1: uDY[13,1] - uDZ[13,1] + ustk[13,1] >= 0
 uATKr1: - uDY[13,1] + uDZ[13,1] + ustk[13,1] >= 0
 uATKr1: uDY[13,1] + uDZ[13,1] - udz[13,1] - ustk[13,1] >= 0
 uATKr1: uDY[13,1] - udz[13,1] >= 0
 uATKr1: udy[13,2] - udz[13,2] = 0
 uATKr1: uDY[13,2] - uDZ[13,2] + ustk[13,2] >= 0
 uATKr1: - uDY[13,2] + uDZ[13,2] + ustk[13,2] >= 0
 uATKr1: uDY[13,2] + uDZ[13,2] - udz[13,2] - ustk[13,2] >= 0
 uATKr1: uDY[13,2] - udz[13,2] >= 0
 uATKr1: udy[13,3] - udz[13,3] = 0
 uATKr1: uDY[13,3] - uDZ[13,3] + ustk[13,3] >= 0
 uATKr1: - uDY[13,3] + uDZ[13,3] + ustk[13,3] >= 0
 uATKr1: uDY[13,3] + uDZ[13,3] - udz[13,3] - ustk[13,3] >= 0
 uATKr1: uDY[13,3] - udz[13,3] >= 0
 uATKr1: udy[13,4] - udz[13,4] = 0
 uATKr1: uDY[13,4] - uDZ[13,4] + ustk[13,4] >= 0
 uATKr1: - uDY[13,4] + uDZ[13,4] + ustk[13,4] >= 0
 uATKr1: uDY[13,4] + uDZ[13,4] - udz[13,4] - ustk[13,4] >= 0
 uATKr1: uDY[13,4] - udz[13,4] >= 0
 uATKr1: udy[13,5] - udz[13,5] = 0
 uATKr1: uDY[13,5] - uDZ[13,5] + ustk[13,5] >= 0
 uATKr1: - uDY[13,5] + uDZ[13,5] + ustk[13,5] >= 0
 uATKr1: uDY[13,5] + uDZ[13,5] - udz[13,5] - ustk[13,5] >= 0
 uATKr1: uDY[13,5] - udz[13,5] >= 0
 uATKr1: udy[13,6] - udz[13,6] = 0
 uATKr1: uDY[13,6] - uDZ[13,6] + ustk[13,6] >= 0
 uATKr1: - uDY[13,6] + uDZ[13,6] + ustk[13,6] >= 0
 uATKr1: uDY[13,6] + uDZ[13,6] - udz[13,6] - ustk[13,6] >= 0
 uATKr1: uDY[13,6] - udz[13,6] >= 0
 uATKr1: udy[13,7] - udz[13,7] = 0
 uATKr1: uDY[13,7] - uDZ[13,7] + ustk[13,7] >= 0
 uATKr1: - uDY[13,7] + uDZ[13,7] + ustk[13,7] >= 0
 uATKr1: uDY[13,7] + uDZ[13,7] - udz[13,7] - ustk[13,7] >= 0
 uATKr1: uDY[13,7] - udz[13,7] >= 0
 uATKr1: uDY[13,8] - uDZ[13,8] = 0
 uATKr1: udy[13,8] - udz[13,8] = 0
 uATKr1: uDY[13,9] - uDZ[13,9] = 0
 uATKr1: udy[13,9] - udz[13,9] = 0
 uATKr1: uDY[13,10] - uDZ[13,10] = 0
 uATKr1: udy[13,10] - udz[13,10] = 0
 uATKr1: uDY[13,11] - uDZ[13,11] = 0
 uATKr1: udy[13,11] - udz[13,11] = 0
 uATKr1: uDY[13,12] - uDZ[13,12] = 0
 uATKr1: udy[13,12] - udz[13,12] = 0
 uATKr1: uDY[13,13] - uDZ[13,13] = 0
 uATKr1: udy[13,13] - udz[13,13] = 0
 uATKr1: uDY[13,14] - uDZ[13,14] = 0
 uATKr1: udy[13,14] - udz[13,14] = 0
 uATKr1: uDY[13,15] - uDZ[13,15] = 0
 uATKr1: udy[13,15] - udz[13,15] = 0
 uATKr1: udy[14,0] - udz[14,0] = 0
 uATKr1: uDY[14,0] - uDZ[14,0] + ustk[14,0] >= 0
 uATKr1: - uDY[14,0] + uDZ[14,0] + ustk[14,0] >= 0
 uATKr1: uDY[14,0] + uDZ[14,0] - udz[14,0] - ustk[14,0] >= 0
 uATKr1: uDY[14,0] - udz[14,0] >= 0
 uATKr1: udy[14,1] - udz[14,1] = 0
 uATKr1: uDY[14,1] - uDZ[14,1] + ustk[14,1] >= 0
 uATKr1: - uDY[14,1] + uDZ[14,1] + ustk[14,1] >= 0
 uATKr1: uDY[14,1] + uDZ[14,1] - udz[14,1] - ustk[14,1] >= 0
 uATKr1: uDY[14,1] - udz[14,1] >= 0
 uATKr1: udy[14,2] - udz[14,2] = 0
 uATKr1: uDY[14,2] - uDZ[14,2] + ustk[14,2] >= 0
 uATKr1: - uDY[14,2] + uDZ[14,2] + ustk[14,2] >= 0
 uATKr1: uDY[14,2] + uDZ[14,2] - udz[14,2] - ustk[14,2] >= 0
 uATKr1: uDY[14,2] - udz[14,2] >= 0
 uATKr1: udy[14,3] - udz[14,3] = 0
 uATKr1: uDY[14,3] - uDZ[14,3] + ustk[14,3] >= 0
 uATKr1: - uDY[14,3] + uDZ[14,3] + ustk[14,3] >= 0
 uATKr1: uDY[14,3] + uDZ[14,3] - udz[14,3] - ustk[14,3] >= 0
 uATKr1: uDY[14,3] - udz[14,3] >= 0
 uATKr1: udy[14,4] - udz[14,4] = 0
 uATKr1: uDY[14,4] - uDZ[14,4] + ustk[14,4] >= 0
 uATKr1: - uDY[14,4] + uDZ[14,4] + ustk[14,4] >= 0
 uATKr1: uDY[14,4] + uDZ[14,4] - udz[14,4] - ustk[14,4] >= 0
 uATKr1: uDY[14,4] - udz[14,4] >= 0
 uATKr1: udy[14,5] - udz[14,5] = 0
 uATKr1: uDY[14,5] - uDZ[14,5] + ustk[14,5] >= 0
 uATKr1: - uDY[14,5] + uDZ[14,5] + ustk[14,5] >= 0
 uATKr1: uDY[14,5] + uDZ[14,5] - udz[14,5] - ustk[14,5] >= 0
 uATKr1: uDY[14,5] - udz[14,5] >= 0
 uATKr1: udy[14,6] - udz[14,6] = 0
 uATKr1: uDY[14,6] - uDZ[14,6] + ustk[14,6] >= 0
 uATKr1: - uDY[14,6] + uDZ[14,6] + ustk[14,6] >= 0
 uATKr1: uDY[14,6] + uDZ[14,6] - udz[14,6] - ustk[14,6] >= 0
 uATKr1: uDY[14,6] - udz[14,6] >= 0
 uATKr1: udy[14,7] - udz[14,7] = 0
 uATKr1: uDY[14,7] - uDZ[14,7] + ustk[14,7] >= 0
 uATKr1: - uDY[14,7] + uDZ[14,7] + ustk[14,7] >= 0
 uATKr1: uDY[14,7] + uDZ[14,7] - udz[14,7] - ustk[14,7] >= 0
 uATKr1: uDY[14,7] - udz[14,7] >= 0
 uATKr1: uDY[14,8] - uDZ[14,8] = 0
 uATKr1: udy[14,8] - udz[14,8] = 0
 uATKr1: uDY[14,9] - uDZ[14,9] = 0
 uATKr1: udy[14,9] - udz[14,9] = 0
 uATKr1: uDY[14,10] - uDZ[14,10] = 0
 uATKr1: udy[14,10] - udz[14,10] = 0
 uATKr1: uDY[14,11] - uDZ[14,11] = 0
 uATKr1: udy[14,11] - udz[14,11] = 0
 uATKr1: uDY[14,12] - uDZ[14,12] = 0
 uATKr1: udy[14,12] - udz[14,12] = 0
 uATKr1: uDY[14,13] - uDZ[14,13] = 0
 uATKr1: udy[14,13] - udz[14,13] = 0
 uATKr1: uDY[14,14] - uDZ[14,14] = 0
 uATKr1: udy[14,14] - udz[14,14] = 0
 uATKr1: uDY[14,15] - uDZ[14,15] = 0
 uATKr1: udy[14,15] - udz[14,15] = 0
 uATKr1: udy[15,0] - udz[15,0] = 0
 uATKr1: uDY[15,0] - uDZ[15,0] + ustk[15,0] >= 0
 uATKr1: - uDY[15,0] + uDZ[15,0] + ustk[15,0] >= 0
 uATKr1: uDY[15,0] + uDZ[15,0] - udz[15,0] - ustk[15,0] >= 0
 uATKr1: uDY[15,0] - udz[15,0] >= 0
 uATKr1: udy[15,1] - udz[15,1] = 0
 uATKr1: uDY[15,1] - uDZ[15,1] + ustk[15,1] >= 0
 uATKr1: - uDY[15,1] + uDZ[15,1] + ustk[15,1] >= 0
 uATKr1: uDY[15,1] + uDZ[15,1] - udz[15,1] - ustk[15,1] >= 0
 uATKr1: uDY[15,1] - udz[15,1] >= 0
 uATKr1: udy[15,2] - udz[15,2] = 0
 uATKr1: uDY[15,2] - uDZ[15,2] + ustk[15,2] >= 0
 uATKr1: - uDY[15,2] + uDZ[15,2] + ustk[15,2] >= 0
 uATKr1: uDY[15,2] + uDZ[15,2] - udz[15,2] - ustk[15,2] >= 0
 uATKr1: uDY[15,2] - udz[15,2] >= 0
 uATKr1: udy[15,3] - udz[15,3] = 0
 uATKr1: uDY[15,3] - uDZ[15,3] + ustk[15,3] >= 0
 uATKr1: - uDY[15,3] + uDZ[15,3] + ustk[15,3] >= 0
 uATKr1: uDY[15,3] + uDZ[15,3] - udz[15,3] - ustk[15,3] >= 0
 uATKr1: uDY[15,3] - udz[15,3] >= 0
 uATKr1: udy[15,4] - udz[15,4] = 0
 uATKr1: uDY[15,4] - uDZ[15,4] + ustk[15,4] >= 0
 uATKr1: - uDY[15,4] + uDZ[15,4] + ustk[15,4] >= 0
 uATKr1: uDY[15,4] + uDZ[15,4] - udz[15,4] - ustk[15,4] >= 0
 uATKr1: uDY[15,4] - udz[15,4] >= 0
 uATKr1: udy[15,5] - udz[15,5] = 0
 uATKr1: uDY[15,5] - uDZ[15,5] + ustk[15,5] >= 0
 uATKr1: - uDY[15,5] + uDZ[15,5] + ustk[15,5] >= 0
 uATKr1: uDY[15,5] + uDZ[15,5] - udz[15,5] - ustk[15,5] >= 0
 uATKr1: uDY[15,5] - udz[15,5] >= 0
 uATKr1: udy[15,6] - udz[15,6] = 0
 uATKr1: uDY[15,6] - uDZ[15,6] + ustk[15,6] >= 0
 uATKr1: - uDY[15,6] + uDZ[15,6] + ustk[15,6] >= 0
 uATKr1: uDY[15,6] + uDZ[15,6] - udz[15,6] - ustk[15,6] >= 0
 uATKr1: uDY[15,6] - udz[15,6] >= 0
 uATKr1: udy[15,7] - udz[15,7] = 0
 uATKr1: uDY[15,7] - uDZ[15,7] + ustk[15,7] >= 0
 uATKr1: - uDY[15,7] + uDZ[15,7] + ustk[15,7] >= 0
 uATKr1: uDY[15,7] + uDZ[15,7] - udz[15,7] - ustk[15,7] >= 0
 uATKr1: uDY[15,7] - udz[15,7] >= 0
 uATKr1: uDY[15,8] - uDZ[15,8] = 0
 uATKr1: udy[15,8] - udz[15,8] = 0
 uATKr1: uDY[15,9] - uDZ[15,9] = 0
 uATKr1: udy[15,9] - udz[15,9] = 0
 uATKr1: uDY[15,10] - uDZ[15,10] = 0
 uATKr1: udy[15,10] - udz[15,10] = 0
 uATKr1: uDY[15,11] - uDZ[15,11] = 0
 uATKr1: udy[15,11] - udz[15,11] = 0
 uATKr1: uDY[15,12] - uDZ[15,12] = 0
 uATKr1: udy[15,12] - udz[15,12] = 0
 uATKr1: uDY[15,13] - uDZ[15,13] = 0
 uATKr1: udy[15,13] - udz[15,13] = 0
 uATKr1: uDY[15,14] - uDZ[15,14] = 0
 uATKr1: udy[15,14] - udz[15,14] = 0
 uATKr1: uDY[15,15] - uDZ[15,15] = 0
 uATKr1: udy[15,15] - udz[15,15] = 0
 uATKr1: udy[16,0] - udz[16,0] = 0
 uATKr1: uDY[16,0] - uDZ[16,0] + ustk[16,0] >= 0
 uATKr1: - uDY[16,0] + uDZ[16,0] + ustk[16,0] >= 0
 uATKr1: uDY[16,0] + uDZ[16,0] - udz[16,0] - ustk[16,0] >= 0
 uATKr1: uDY[16,0] - udz[16,0] >= 0
 uATKr1: udy[16,1] - udz[16,1] = 0
 uATKr1: uDY[16,1] - uDZ[16,1] + ustk[16,1] >= 0
 uATKr1: - uDY[16,1] + uDZ[16,1] + ustk[16,1] >= 0
 uATKr1: uDY[16,1] + uDZ[16,1] - udz[16,1] - ustk[16,1] >= 0
 uATKr1: uDY[16,1] - udz[16,1] >= 0
 uATKr1: udy[16,2] - udz[16,2] = 0
 uATKr1: uDY[16,2] - uDZ[16,2] + ustk[16,2] >= 0
 uATKr1: - uDY[16,2] + uDZ[16,2] + ustk[16,2] >= 0
 uATKr1: uDY[16,2] + uDZ[16,2] - udz[16,2] - ustk[16,2] >= 0
 uATKr1: uDY[16,2] - udz[16,2] >= 0
 uATKr1: udy[16,3] - udz[16,3] = 0
 uATKr1: uDY[16,3] - uDZ[16,3] + ustk[16,3] >= 0
 uATKr1: - uDY[16,3] + uDZ[16,3] + ustk[16,3] >= 0
 uATKr1: uDY[16,3] + uDZ[16,3] - udz[16,3] - ustk[16,3] >= 0
 uATKr1: uDY[16,3] - udz[16,3] >= 0
 uATKr1: udy[16,4] - udz[16,4] = 0
 uATKr1: uDY[16,4] - uDZ[16,4] + ustk[16,4] >= 0
 uATKr1: - uDY[16,4] + uDZ[16,4] + ustk[16,4] >= 0
 uATKr1: uDY[16,4] + uDZ[16,4] - udz[16,4] - ustk[16,4] >= 0
 uATKr1: uDY[16,4] - udz[16,4] >= 0
 uATKr1: udy[16,5] - udz[16,5] = 0
 uATKr1: uDY[16,5] - uDZ[16,5] + ustk[16,5] >= 0
 uATKr1: - uDY[16,5] + uDZ[16,5] + ustk[16,5] >= 0
 uATKr1: uDY[16,5] + uDZ[16,5] - udz[16,5] - ustk[16,5] >= 0
 uATKr1: uDY[16,5] - udz[16,5] >= 0
 uATKr1: udy[16,6] - udz[16,6] = 0
 uATKr1: uDY[16,6] - uDZ[16,6] + ustk[16,6] >= 0
 uATKr1: - uDY[16,6] + uDZ[16,6] + ustk[16,6] >= 0
 uATKr1: uDY[16,6] + uDZ[16,6] - udz[16,6] - ustk[16,6] >= 0
 uATKr1: uDY[16,6] - udz[16,6] >= 0
 uATKr1: udy[16,7] - udz[16,7] = 0
 uATKr1: uDY[16,7] - uDZ[16,7] + ustk[16,7] >= 0
 uATKr1: - uDY[16,7] + uDZ[16,7] + ustk[16,7] >= 0
 uATKr1: uDY[16,7] + uDZ[16,7] - udz[16,7] - ustk[16,7] >= 0
 uATKr1: uDY[16,7] - udz[16,7] >= 0
 uATKr1: uDY[16,8] - uDZ[16,8] = 0
 uATKr1: udy[16,8] - udz[16,8] = 0
 uATKr1: uDY[16,9] - uDZ[16,9] = 0
 uATKr1: udy[16,9] - udz[16,9] = 0
 uATKr1: uDY[16,10] - uDZ[16,10] = 0
 uATKr1: udy[16,10] - udz[16,10] = 0
 uATKr1: uDY[16,11] - uDZ[16,11] = 0
 uATKr1: udy[16,11] - udz[16,11] = 0
 uATKr1: uDY[16,12] - uDZ[16,12] = 0
 uATKr1: udy[16,12] - udz[16,12] = 0
 uATKr1: uDY[16,13] - uDZ[16,13] = 0
 uATKr1: udy[16,13] - udz[16,13] = 0
 uATKr1: uDY[16,14] - uDZ[16,14] = 0
 uATKr1: udy[16,14] - udz[16,14] = 0
 uATKr1: uDY[16,15] - uDZ[16,15] = 0
 uATKr1: udy[16,15] - udz[16,15] = 0
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
 uSR: uDZ[12,0] - uDW[12,0] = 0
 uSR: udz[12,0] - udw[12,0] = 0
 uSR: uDZ[12,1] - uDW[12,1] = 0
 uSR: udz[12,1] - udw[12,1] = 0
 uSR: uDZ[12,2] - uDW[12,2] = 0
 uSR: udz[12,2] - udw[12,2] = 0
 uSR: uDZ[12,3] - uDW[12,3] = 0
 uSR: udz[12,3] - udw[12,3] = 0
 uSR: uDZ[12,4] - uDW[12,5] = 0
 uSR: udz[12,4] - udw[12,5] = 0
 uSR: uDZ[12,5] - uDW[12,6] = 0
 uSR: udz[12,5] - udw[12,6] = 0
 uSR: uDZ[12,6] - uDW[12,7] = 0
 uSR: udz[12,6] - udw[12,7] = 0
 uSR: uDZ[12,7] - uDW[12,4] = 0
 uSR: udz[12,7] - udw[12,4] = 0
 uSR: uDZ[12,8] - uDW[12,10] = 0
 uSR: udz[12,8] - udw[12,10] = 0
 uSR: uDZ[12,9] - uDW[12,11] = 0
 uSR: udz[12,9] - udw[12,11] = 0
 uSR: uDZ[12,10] - uDW[12,8] = 0
 uSR: udz[12,10] - udw[12,8] = 0
 uSR: uDZ[12,11] - uDW[12,9] = 0
 uSR: udz[12,11] - udw[12,9] = 0
 uSR: uDZ[12,12] - uDW[12,15] = 0
 uSR: udz[12,12] - udw[12,15] = 0
 uSR: uDZ[12,13] - uDW[12,12] = 0
 uSR: udz[12,13] - udw[12,12] = 0
 uSR: uDZ[12,14] - uDW[12,13] = 0
 uSR: udz[12,14] - udw[12,13] = 0
 uSR: uDZ[12,15] - uDW[12,14] = 0
 uSR: udz[12,15] - udw[12,14] = 0
 uSR: uDZ[13,0] - uDW[13,0] = 0
 uSR: udz[13,0] - udw[13,0] = 0
 uSR: uDZ[13,1] - uDW[13,1] = 0
 uSR: udz[13,1] - udw[13,1] = 0
 uSR: uDZ[13,2] - uDW[13,2] = 0
 uSR: udz[13,2] - udw[13,2] = 0
 uSR: uDZ[13,3] - uDW[13,3] = 0
 uSR: udz[13,3] - udw[13,3] = 0
 uSR: uDZ[13,4] - uDW[13,5] = 0
 uSR: udz[13,4] - udw[13,5] = 0
 uSR: uDZ[13,5] - uDW[13,6] = 0
 uSR: udz[13,5] - udw[13,6] = 0
 uSR: uDZ[13,6] - uDW[13,7] = 0
 uSR: udz[13,6] - udw[13,7] = 0
 uSR: uDZ[13,7] - uDW[13,4] = 0
 uSR: udz[13,7] - udw[13,4] = 0
 uSR: uDZ[13,8] - uDW[13,10] = 0
 uSR: udz[13,8] - udw[13,10] = 0
 uSR: uDZ[13,9] - uDW[13,11] = 0
 uSR: udz[13,9] - udw[13,11] = 0
 uSR: uDZ[13,10] - uDW[13,8] = 0
 uSR: udz[13,10] - udw[13,8] = 0
 uSR: uDZ[13,11] - uDW[13,9] = 0
 uSR: udz[13,11] - udw[13,9] = 0
 uSR: uDZ[13,12] - uDW[13,15] = 0
 uSR: udz[13,12] - udw[13,15] = 0
 uSR: uDZ[13,13] - uDW[13,12] = 0
 uSR: udz[13,13] - udw[13,12] = 0
 uSR: uDZ[13,14] - uDW[13,13] = 0
 uSR: udz[13,14] - udw[13,13] = 0
 uSR: uDZ[13,15] - uDW[13,14] = 0
 uSR: udz[13,15] - udw[13,14] = 0
 uSR: uDZ[14,0] - uDW[14,0] = 0
 uSR: udz[14,0] - udw[14,0] = 0
 uSR: uDZ[14,1] - uDW[14,1] = 0
 uSR: udz[14,1] - udw[14,1] = 0
 uSR: uDZ[14,2] - uDW[14,2] = 0
 uSR: udz[14,2] - udw[14,2] = 0
 uSR: uDZ[14,3] - uDW[14,3] = 0
 uSR: udz[14,3] - udw[14,3] = 0
 uSR: uDZ[14,4] - uDW[14,5] = 0
 uSR: udz[14,4] - udw[14,5] = 0
 uSR: uDZ[14,5] - uDW[14,6] = 0
 uSR: udz[14,5] - udw[14,6] = 0
 uSR: uDZ[14,6] - uDW[14,7] = 0
 uSR: udz[14,6] - udw[14,7] = 0
 uSR: uDZ[14,7] - uDW[14,4] = 0
 uSR: udz[14,7] - udw[14,4] = 0
 uSR: uDZ[14,8] - uDW[14,10] = 0
 uSR: udz[14,8] - udw[14,10] = 0
 uSR: uDZ[14,9] - uDW[14,11] = 0
 uSR: udz[14,9] - udw[14,11] = 0
 uSR: uDZ[14,10] - uDW[14,8] = 0
 uSR: udz[14,10] - udw[14,8] = 0
 uSR: uDZ[14,11] - uDW[14,9] = 0
 uSR: udz[14,11] - udw[14,9] = 0
 uSR: uDZ[14,12] - uDW[14,15] = 0
 uSR: udz[14,12] - udw[14,15] = 0
 uSR: uDZ[14,13] - uDW[14,12] = 0
 uSR: udz[14,13] - udw[14,12] = 0
 uSR: uDZ[14,14] - uDW[14,13] = 0
 uSR: udz[14,14] - udw[14,13] = 0
 uSR: uDZ[14,15] - uDW[14,14] = 0
 uSR: udz[14,15] - udw[14,14] = 0
 uSR: uDZ[15,0] - uDW[15,0] = 0
 uSR: udz[15,0] - udw[15,0] = 0
 uSR: uDZ[15,1] - uDW[15,1] = 0
 uSR: udz[15,1] - udw[15,1] = 0
 uSR: uDZ[15,2] - uDW[15,2] = 0
 uSR: udz[15,2] - udw[15,2] = 0
 uSR: uDZ[15,3] - uDW[15,3] = 0
 uSR: udz[15,3] - udw[15,3] = 0
 uSR: uDZ[15,4] - uDW[15,5] = 0
 uSR: udz[15,4] - udw[15,5] = 0
 uSR: uDZ[15,5] - uDW[15,6] = 0
 uSR: udz[15,5] - udw[15,6] = 0
 uSR: uDZ[15,6] - uDW[15,7] = 0
 uSR: udz[15,6] - udw[15,7] = 0
 uSR: uDZ[15,7] - uDW[15,4] = 0
 uSR: udz[15,7] - udw[15,4] = 0
 uSR: uDZ[15,8] - uDW[15,10] = 0
 uSR: udz[15,8] - udw[15,10] = 0
 uSR: uDZ[15,9] - uDW[15,11] = 0
 uSR: udz[15,9] - udw[15,11] = 0
 uSR: uDZ[15,10] - uDW[15,8] = 0
 uSR: udz[15,10] - udw[15,8] = 0
 uSR: uDZ[15,11] - uDW[15,9] = 0
 uSR: udz[15,11] - udw[15,9] = 0
 uSR: uDZ[15,12] - uDW[15,15] = 0
 uSR: udz[15,12] - udw[15,15] = 0
 uSR: uDZ[15,13] - uDW[15,12] = 0
 uSR: udz[15,13] - udw[15,12] = 0
 uSR: uDZ[15,14] - uDW[15,13] = 0
 uSR: udz[15,14] - udw[15,13] = 0
 uSR: uDZ[15,15] - uDW[15,14] = 0
 uSR: udz[15,15] - udw[15,14] = 0
 uSR: uDZ[16,0] - uDW[16,0] = 0
 uSR: udz[16,0] - udw[16,0] = 0
 uSR: uDZ[16,1] - uDW[16,1] = 0
 uSR: udz[16,1] - udw[16,1] = 0
 uSR: uDZ[16,2] - uDW[16,2] = 0
 uSR: udz[16,2] - udw[16,2] = 0
 uSR: uDZ[16,3] - uDW[16,3] = 0
 uSR: udz[16,3] - udw[16,3] = 0
 uSR: uDZ[16,4] - uDW[16,5] = 0
 uSR: udz[16,4] - udw[16,5] = 0
 uSR: uDZ[16,5] - uDW[16,6] = 0
 uSR: udz[16,5] - udw[16,6] = 0
 uSR: uDZ[16,6] - uDW[16,7] = 0
 uSR: udz[16,6] - udw[16,7] = 0
 uSR: uDZ[16,7] - uDW[16,4] = 0
 uSR: udz[16,7] - udw[16,4] = 0
 uSR: uDZ[16,8] - uDW[16,10] = 0
 uSR: udz[16,8] - udw[16,10] = 0
 uSR: uDZ[16,9] - uDW[16,11] = 0
 uSR: udz[16,9] - udw[16,11] = 0
 uSR: uDZ[16,10] - uDW[16,8] = 0
 uSR: udz[16,10] - udw[16,8] = 0
 uSR: uDZ[16,11] - uDW[16,9] = 0
 uSR: udz[16,11] - udw[16,9] = 0
 uSR: uDZ[16,12] - uDW[16,15] = 0
 uSR: udz[16,12] - udw[16,15] = 0
 uSR: uDZ[16,13] - uDW[16,12] = 0
 uSR: udz[16,13] - udw[16,12] = 0
 uSR: uDZ[16,14] - uDW[16,13] = 0
 uSR: udz[16,14] - udw[16,13] = 0
 uSR: uDZ[16,15] - uDW[16,14] = 0
 uSR: udz[16,15] - udw[16,14] = 0
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
 ueMC0: - udx[5,4] + udw[4,0] = 0
 ueMC0: - uDX[5,4] + uDW[4,0] = 0
 ueMC1: udw[4,4] - udw[4,8] >= 0
 ueMC1: - udx[5,8] + udw[4,4] >= 0
 ueMC1: - udx[5,8] + udw[4,4] - udw[4,8] <= 0
 ueMC1: uDX[5,8] - uDW[4,4] + uDW[4,8] >= 0
 ueMC1: uDX[5,8] + uDW[4,4] - uDW[4,8] >= 0
 ueMC1: - uDX[5,8] + uDW[4,4] + uDW[4,8] >= 0
 ueMC1: - uDX[5,8] + uDW[4,4] - uDW[4,8] <= 0
 ueMC2: - udx[5,4] + udw[4,8] >= 0
 ueMC2: - udx[5,12] + udw[4,8] >= 0
 ueMC2: - udx[5,4] - udx[5,12] + udw[4,8] <= 0
 ueMC2: uDX[5,4] + uDX[5,12] - uDW[4,8] >= 0
 ueMC2: - uDX[5,4] + uDX[5,12] + uDW[4,8] >= 0
 ueMC2: uDX[5,4] - uDX[5,12] + uDW[4,8] >= 0
 ueMC2: - uDX[5,4] - uDX[5,12] + uDW[4,8] <= 0
 ueMC3: - udx[5,0] + udw[4,12] >= 0
 ueMC3: - udx[5,12] + udw[4,12] >= 0
 ueMC3: - udx[5,0] - udx[5,12] + udw[4,12] <= 0
 ueMC3: uDX[5,0] + uDX[5,12] - uDW[4,12] >= 0
 ueMC3: - uDX[5,0] + uDX[5,12] + uDW[4,12] >= 0
 ueMC3: uDX[5,0] - uDX[5,12] + uDW[4,12] >= 0
 ueMC3: - uDX[5,0] - uDX[5,12] + uDW[4,12] <= 0
 ueMC0: - udx[5,5] + udw[4,1] = 0
 ueMC0: - uDX[5,5] + uDW[4,1] = 0
 ueMC1: udw[4,5] - udw[4,9] >= 0
 ueMC1: - udx[5,9] + udw[4,5] >= 0
 ueMC1: - udx[5,9] + udw[4,5] - udw[4,9] <= 0
 ueMC1: uDX[5,9] - uDW[4,5] + uDW[4,9] >= 0
 ueMC1: uDX[5,9] + uDW[4,5] - uDW[4,9] >= 0
 ueMC1: - uDX[5,9] + uDW[4,5] + uDW[4,9] >= 0
 ueMC1: - uDX[5,9] + uDW[4,5] - uDW[4,9] <= 0
 ueMC2: - udx[5,5] + udw[4,9] >= 0
 ueMC2: - udx[5,13] + udw[4,9] >= 0
 ueMC2: - udx[5,5] - udx[5,13] + udw[4,9] <= 0
 ueMC2: uDX[5,5] + uDX[5,13] - uDW[4,9] >= 0
 ueMC2: - uDX[5,5] + uDX[5,13] + uDW[4,9] >= 0
 ueMC2: uDX[5,5] - uDX[5,13] + uDW[4,9] >= 0
 ueMC2: - uDX[5,5] - uDX[5,13] + uDW[4,9] <= 0
 ueMC3: - udx[5,1] + udw[4,13] >= 0
 ueMC3: - udx[5,13] + udw[4,13] >= 0
 ueMC3: - udx[5,1] - udx[5,13] + udw[4,13] <= 0
 ueMC3: uDX[5,1] + uDX[5,13] - uDW[4,13] >= 0
 ueMC3: - uDX[5,1] + uDX[5,13] + uDW[4,13] >= 0
 ueMC3: uDX[5,1] - uDX[5,13] + uDW[4,13] >= 0
 ueMC3: - uDX[5,1] - uDX[5,13] + uDW[4,13] <= 0
 ueMC0: - udx[5,6] + udw[4,2] = 0
 ueMC0: - uDX[5,6] + uDW[4,2] = 0
 ueMC1: udw[4,6] - udw[4,10] >= 0
 ueMC1: - udx[5,10] + udw[4,6] >= 0
 ueMC1: - udx[5,10] + udw[4,6] - udw[4,10] <= 0
 ueMC1: uDX[5,10] - uDW[4,6] + uDW[4,10] >= 0
 ueMC1: uDX[5,10] + uDW[4,6] - uDW[4,10] >= 0
 ueMC1: - uDX[5,10] + uDW[4,6] + uDW[4,10] >= 0
 ueMC1: - uDX[5,10] + uDW[4,6] - uDW[4,10] <= 0
 ueMC2: - udx[5,6] + udw[4,10] >= 0
 ueMC2: - udx[5,14] + udw[4,10] >= 0
 ueMC2: - udx[5,6] - udx[5,14] + udw[4,10] <= 0
 ueMC2: uDX[5,6] + uDX[5,14] - uDW[4,10] >= 0
 ueMC2: - uDX[5,6] + uDX[5,14] + uDW[4,10] >= 0
 ueMC2: uDX[5,6] - uDX[5,14] + uDW[4,10] >= 0
 ueMC2: - uDX[5,6] - uDX[5,14] + uDW[4,10] <= 0
 ueMC3: - udx[5,2] + udw[4,14] >= 0
 ueMC3: - udx[5,14] + udw[4,14] >= 0
 ueMC3: - udx[5,2] - udx[5,14] + udw[4,14] <= 0
 ueMC3: uDX[5,2] + uDX[5,14] - uDW[4,14] >= 0
 ueMC3: - uDX[5,2] + uDX[5,14] + uDW[4,14] >= 0
 ueMC3: uDX[5,2] - uDX[5,14] + uDW[4,14] >= 0
 ueMC3: - uDX[5,2] - uDX[5,14] + uDW[4,14] <= 0
 ueMC0: - udx[5,7] + udw[4,3] = 0
 ueMC0: - uDX[5,7] + uDW[4,3] = 0
 ueMC1: udw[4,7] - udw[4,11] >= 0
 ueMC1: - udx[5,11] + udw[4,7] >= 0
 ueMC1: - udx[5,11] + udw[4,7] - udw[4,11] <= 0
 ueMC1: uDX[5,11] - uDW[4,7] + uDW[4,11] >= 0
 ueMC1: uDX[5,11] + uDW[4,7] - uDW[4,11] >= 0
 ueMC1: - uDX[5,11] + uDW[4,7] + uDW[4,11] >= 0
 ueMC1: - uDX[5,11] + uDW[4,7] - uDW[4,11] <= 0
 ueMC2: - udx[5,7] + udw[4,11] >= 0
 ueMC2: - udx[5,15] + udw[4,11] >= 0
 ueMC2: - udx[5,7] - udx[5,15] + udw[4,11] <= 0
 ueMC2: uDX[5,7] + uDX[5,15] - uDW[4,11] >= 0
 ueMC2: - uDX[5,7] + uDX[5,15] + uDW[4,11] >= 0
 ueMC2: uDX[5,7] - uDX[5,15] + uDW[4,11] >= 0
 ueMC2: - uDX[5,7] - uDX[5,15] + uDW[4,11] <= 0
 ueMC3: - udx[5,3] + udw[4,15] >= 0
 ueMC3: - udx[5,15] + udw[4,15] >= 0
 ueMC3: - udx[5,3] - udx[5,15] + udw[4,15] <= 0
 ueMC3: uDX[5,3] + uDX[5,15] - uDW[4,15] >= 0
 ueMC3: - uDX[5,3] + uDX[5,15] + uDW[4,15] >= 0
 ueMC3: uDX[5,3] - uDX[5,15] + uDW[4,15] >= 0
 ueMC3: - uDX[5,3] - uDX[5,15] + uDW[4,15] <= 0
 ueMC0: - udx[6,4] + udw[5,0] = 0
 ueMC0: - uDX[6,4] + uDW[5,0] = 0
 ueMC1: udw[5,4] - udw[5,8] >= 0
 ueMC1: - udx[6,8] + udw[5,4] >= 0
 ueMC1: - udx[6,8] + udw[5,4] - udw[5,8] <= 0
 ueMC1: uDX[6,8] - uDW[5,4] + uDW[5,8] >= 0
 ueMC1: uDX[6,8] + uDW[5,4] - uDW[5,8] >= 0
 ueMC1: - uDX[6,8] + uDW[5,4] + uDW[5,8] >= 0
 ueMC1: - uDX[6,8] + uDW[5,4] - uDW[5,8] <= 0
 ueMC2: - udx[6,4] + udw[5,8] >= 0
 ueMC2: - udx[6,12] + udw[5,8] >= 0
 ueMC2: - udx[6,4] - udx[6,12] + udw[5,8] <= 0
 ueMC2: uDX[6,4] + uDX[6,12] - uDW[5,8] >= 0
 ueMC2: - uDX[6,4] + uDX[6,12] + uDW[5,8] >= 0
 ueMC2: uDX[6,4] - uDX[6,12] + uDW[5,8] >= 0
 ueMC2: - uDX[6,4] - uDX[6,12] + uDW[5,8] <= 0
 ueMC3: - udx[6,0] + udw[5,12] >= 0
 ueMC3: - udx[6,12] + udw[5,12] >= 0
 ueMC3: - udx[6,0] - udx[6,12] + udw[5,12] <= 0
 ueMC3: uDX[6,0] + uDX[6,12] - uDW[5,12] >= 0
 ueMC3: - uDX[6,0] + uDX[6,12] + uDW[5,12] >= 0
 ueMC3: uDX[6,0] - uDX[6,12] + uDW[5,12] >= 0
 ueMC3: - uDX[6,0] - uDX[6,12] + uDW[5,12] <= 0
 ueMC0: - udx[6,5] + udw[5,1] = 0
 ueMC0: - uDX[6,5] + uDW[5,1] = 0
 ueMC1: udw[5,5] - udw[5,9] >= 0
 ueMC1: - udx[6,9] + udw[5,5] >= 0
 ueMC1: - udx[6,9] + udw[5,5] - udw[5,9] <= 0
 ueMC1: uDX[6,9] - uDW[5,5] + uDW[5,9] >= 0
 ueMC1: uDX[6,9] + uDW[5,5] - uDW[5,9] >= 0
 ueMC1: - uDX[6,9] + uDW[5,5] + uDW[5,9] >= 0
 ueMC1: - uDX[6,9] + uDW[5,5] - uDW[5,9] <= 0
 ueMC2: - udx[6,5] + udw[5,9] >= 0
 ueMC2: - udx[6,13] + udw[5,9] >= 0
 ueMC2: - udx[6,5] - udx[6,13] + udw[5,9] <= 0
 ueMC2: uDX[6,5] + uDX[6,13] - uDW[5,9] >= 0
 ueMC2: - uDX[6,5] + uDX[6,13] + uDW[5,9] >= 0
 ueMC2: uDX[6,5] - uDX[6,13] + uDW[5,9] >= 0
 ueMC2: - uDX[6,5] - uDX[6,13] + uDW[5,9] <= 0
 ueMC3: - udx[6,1] + udw[5,13] >= 0
 ueMC3: - udx[6,13] + udw[5,13] >= 0
 ueMC3: - udx[6,1] - udx[6,13] + udw[5,13] <= 0
 ueMC3: uDX[6,1] + uDX[6,13] - uDW[5,13] >= 0
 ueMC3: - uDX[6,1] + uDX[6,13] + uDW[5,13] >= 0
 ueMC3: uDX[6,1] - uDX[6,13] + uDW[5,13] >= 0
 ueMC3: - uDX[6,1] - uDX[6,13] + uDW[5,13] <= 0
 ueMC0: - udx[6,6] + udw[5,2] = 0
 ueMC0: - uDX[6,6] + uDW[5,2] = 0
 ueMC1: udw[5,6] - udw[5,10] >= 0
 ueMC1: - udx[6,10] + udw[5,6] >= 0
 ueMC1: - udx[6,10] + udw[5,6] - udw[5,10] <= 0
 ueMC1: uDX[6,10] - uDW[5,6] + uDW[5,10] >= 0
 ueMC1: uDX[6,10] + uDW[5,6] - uDW[5,10] >= 0
 ueMC1: - uDX[6,10] + uDW[5,6] + uDW[5,10] >= 0
 ueMC1: - uDX[6,10] + uDW[5,6] - uDW[5,10] <= 0
 ueMC2: - udx[6,6] + udw[5,10] >= 0
 ueMC2: - udx[6,14] + udw[5,10] >= 0
 ueMC2: - udx[6,6] - udx[6,14] + udw[5,10] <= 0
 ueMC2: uDX[6,6] + uDX[6,14] - uDW[5,10] >= 0
 ueMC2: - uDX[6,6] + uDX[6,14] + uDW[5,10] >= 0
 ueMC2: uDX[6,6] - uDX[6,14] + uDW[5,10] >= 0
 ueMC2: - uDX[6,6] - uDX[6,14] + uDW[5,10] <= 0
 ueMC3: - udx[6,2] + udw[5,14] >= 0
 ueMC3: - udx[6,14] + udw[5,14] >= 0
 ueMC3: - udx[6,2] - udx[6,14] + udw[5,14] <= 0
 ueMC3: uDX[6,2] + uDX[6,14] - uDW[5,14] >= 0
 ueMC3: - uDX[6,2] + uDX[6,14] + uDW[5,14] >= 0
 ueMC3: uDX[6,2] - uDX[6,14] + uDW[5,14] >= 0
 ueMC3: - uDX[6,2] - uDX[6,14] + uDW[5,14] <= 0
 ueMC0: - udx[6,7] + udw[5,3] = 0
 ueMC0: - uDX[6,7] + uDW[5,3] = 0
 ueMC1: udw[5,7] - udw[5,11] >= 0
 ueMC1: - udx[6,11] + udw[5,7] >= 0
 ueMC1: - udx[6,11] + udw[5,7] - udw[5,11] <= 0
 ueMC1: uDX[6,11] - uDW[5,7] + uDW[5,11] >= 0
 ueMC1: uDX[6,11] + uDW[5,7] - uDW[5,11] >= 0
 ueMC1: - uDX[6,11] + uDW[5,7] + uDW[5,11] >= 0
 ueMC1: - uDX[6,11] + uDW[5,7] - uDW[5,11] <= 0
 ueMC2: - udx[6,7] + udw[5,11] >= 0
 ueMC2: - udx[6,15] + udw[5,11] >= 0
 ueMC2: - udx[6,7] - udx[6,15] + udw[5,11] <= 0
 ueMC2: uDX[6,7] + uDX[6,15] - uDW[5,11] >= 0
 ueMC2: - uDX[6,7] + uDX[6,15] + uDW[5,11] >= 0
 ueMC2: uDX[6,7] - uDX[6,15] + uDW[5,11] >= 0
 ueMC2: - uDX[6,7] - uDX[6,15] + uDW[5,11] <= 0
 ueMC3: - udx[6,3] + udw[5,15] >= 0
 ueMC3: - udx[6,15] + udw[5,15] >= 0
 ueMC3: - udx[6,3] - udx[6,15] + udw[5,15] <= 0
 ueMC3: uDX[6,3] + uDX[6,15] - uDW[5,15] >= 0
 ueMC3: - uDX[6,3] + uDX[6,15] + uDW[5,15] >= 0
 ueMC3: uDX[6,3] - uDX[6,15] + uDW[5,15] >= 0
 ueMC3: - uDX[6,3] - uDX[6,15] + uDW[5,15] <= 0
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
 udMC0: udx[13,0] - udw[12,12] >= 0
 udMC0: udx[13,0] - udx[13,12] >= 0
 udMC0: udx[13,0] - udx[13,12] - udw[12,12] <= 0
 udMC0: - uDX[13,0] + uDX[13,12] + uDW[12,12] >= 0
 udMC0: uDX[13,0] + uDX[13,12] - uDW[12,12] >= 0
 udMC0: uDX[13,0] - uDX[13,12] + uDW[12,12] >= 0
 udMC0: uDX[13,0] - uDX[13,12] - uDW[12,12] <= 0
 udMC1: udx[13,4] - udw[12,0] = 0
 udMC1: uDX[13,4] - uDW[12,0] = 0
 udMC2: udx[13,8] - udw[12,4] >= 0
 udMC2: udx[13,8] - udw[12,8] >= 0
 udMC2: udx[13,8] - udw[12,4] - udw[12,8] <= 0
 udMC2: - uDX[13,8] + uDW[12,4] + uDW[12,8] >= 0
 udMC2: uDX[13,8] - uDW[12,4] + uDW[12,8] >= 0
 udMC2: uDX[13,8] + uDW[12,4] - uDW[12,8] >= 0
 udMC2: uDX[13,8] - uDW[12,4] - uDW[12,8] <= 0
 udMC3: udx[13,12] - udw[12,0] >= 0
 udMC3: udx[13,12] - udw[12,8] >= 0
 udMC3: udx[13,12] - udw[12,0] - udw[12,8] <= 0
 udMC3: - uDX[13,12] + uDW[12,0] + uDW[12,8] >= 0
 udMC3: uDX[13,12] - uDW[12,0] + uDW[12,8] >= 0
 udMC3: uDX[13,12] + uDW[12,0] - uDW[12,8] >= 0
 udMC3: uDX[13,12] - uDW[12,0] - uDW[12,8] <= 0
 udMC0: udx[13,1] - udw[12,13] >= 0
 udMC0: udx[13,1] - udx[13,13] >= 0
 udMC0: udx[13,1] - udx[13,13] - udw[12,13] <= 0
 udMC0: - uDX[13,1] + uDX[13,13] + uDW[12,13] >= 0
 udMC0: uDX[13,1] + uDX[13,13] - uDW[12,13] >= 0
 udMC0: uDX[13,1] - uDX[13,13] + uDW[12,13] >= 0
 udMC0: uDX[13,1] - uDX[13,13] - uDW[12,13] <= 0
 udMC1: udx[13,5] - udw[12,1] = 0
 udMC1: uDX[13,5] - uDW[12,1] = 0
 udMC2: udx[13,9] - udw[12,5] >= 0
 udMC2: udx[13,9] - udw[12,9] >= 0
 udMC2: udx[13,9] - udw[12,5] - udw[12,9] <= 0
 udMC2: - uDX[13,9] + uDW[12,5] + uDW[12,9] >= 0
 udMC2: uDX[13,9] - uDW[12,5] + uDW[12,9] >= 0
 udMC2: uDX[13,9] + uDW[12,5] - uDW[12,9] >= 0
 udMC2: uDX[13,9] - uDW[12,5] - uDW[12,9] <= 0
 udMC3: udx[13,13] - udw[12,1] >= 0
 udMC3: udx[13,13] - udw[12,9] >= 0
 udMC3: udx[13,13] - udw[12,1] - udw[12,9] <= 0
 udMC3: - uDX[13,13] + uDW[12,1] + uDW[12,9] >= 0
 udMC3: uDX[13,13] - uDW[12,1] + uDW[12,9] >= 0
 udMC3: uDX[13,13] + uDW[12,1] - uDW[12,9] >= 0
 udMC3: uDX[13,13] - uDW[12,1] - uDW[12,9] <= 0
 udMC0: udx[13,2] - udw[12,14] >= 0
 udMC0: udx[13,2] - udx[13,14] >= 0
 udMC0: udx[13,2] - udx[13,14] - udw[12,14] <= 0
 udMC0: - uDX[13,2] + uDX[13,14] + uDW[12,14] >= 0
 udMC0: uDX[13,2] + uDX[13,14] - uDW[12,14] >= 0
 udMC0: uDX[13,2] - uDX[13,14] + uDW[12,14] >= 0
 udMC0: uDX[13,2] - uDX[13,14] - uDW[12,14] <= 0
 udMC1: udx[13,6] - udw[12,2] = 0
 udMC1: uDX[13,6] - uDW[12,2] = 0
 udMC2: udx[13,10] - udw[12,6] >= 0
 udMC2: udx[13,10] - udw[12,10] >= 0
 udMC2: udx[13,10] - udw[12,6] - udw[12,10] <= 0
 udMC2: - uDX[13,10] + uDW[12,6] + uDW[12,10] >= 0
 udMC2: uDX[13,10] - uDW[12,6] + uDW[12,10] >= 0
 udMC2: uDX[13,10] + uDW[12,6] - uDW[12,10] >= 0
 udMC2: uDX[13,10] - uDW[12,6] - uDW[12,10] <= 0
 udMC3: udx[13,14] - udw[12,2] >= 0
 udMC3: udx[13,14] - udw[12,10] >= 0
 udMC3: udx[13,14] - udw[12,2] - udw[12,10] <= 0
 udMC3: - uDX[13,14] + uDW[12,2] + uDW[12,10] >= 0
 udMC3: uDX[13,14] - uDW[12,2] + uDW[12,10] >= 0
 udMC3: uDX[13,14] + uDW[12,2] - uDW[12,10] >= 0
 udMC3: uDX[13,14] - uDW[12,2] - uDW[12,10] <= 0
 udMC0: udx[13,3] - udw[12,15] >= 0
 udMC0: udx[13,3] - udx[13,15] >= 0
 udMC0: udx[13,3] - udx[13,15] - udw[12,15] <= 0
 udMC0: - uDX[13,3] + uDX[13,15] + uDW[12,15] >= 0
 udMC0: uDX[13,3] + uDX[13,15] - uDW[12,15] >= 0
 udMC0: uDX[13,3] - uDX[13,15] + uDW[12,15] >= 0
 udMC0: uDX[13,3] - uDX[13,15] - uDW[12,15] <= 0
 udMC1: udx[13,7] - udw[12,3] = 0
 udMC1: uDX[13,7] - uDW[12,3] = 0
 udMC2: udx[13,11] - udw[12,7] >= 0
 udMC2: udx[13,11] - udw[12,11] >= 0
 udMC2: udx[13,11] - udw[12,7] - udw[12,11] <= 0
 udMC2: - uDX[13,11] + uDW[12,7] + uDW[12,11] >= 0
 udMC2: uDX[13,11] - uDW[12,7] + uDW[12,11] >= 0
 udMC2: uDX[13,11] + uDW[12,7] - uDW[12,11] >= 0
 udMC2: uDX[13,11] - uDW[12,7] - uDW[12,11] <= 0
 udMC3: udx[13,15] - udw[12,3] >= 0
 udMC3: udx[13,15] - udw[12,11] >= 0
 udMC3: udx[13,15] - udw[12,3] - udw[12,11] <= 0
 udMC3: - uDX[13,15] + uDW[12,3] + uDW[12,11] >= 0
 udMC3: uDX[13,15] - uDW[12,3] + uDW[12,11] >= 0
 udMC3: uDX[13,15] + uDW[12,3] - uDW[12,11] >= 0
 udMC3: uDX[13,15] - uDW[12,3] - uDW[12,11] <= 0
 udMC0: udx[14,0] - udw[13,12] >= 0
 udMC0: udx[14,0] - udx[14,12] >= 0
 udMC0: udx[14,0] - udx[14,12] - udw[13,12] <= 0
 udMC0: - uDX[14,0] + uDX[14,12] + uDW[13,12] >= 0
 udMC0: uDX[14,0] + uDX[14,12] - uDW[13,12] >= 0
 udMC0: uDX[14,0] - uDX[14,12] + uDW[13,12] >= 0
 udMC0: uDX[14,0] - uDX[14,12] - uDW[13,12] <= 0
 udMC1: udx[14,4] - udw[13,0] = 0
 udMC1: uDX[14,4] - uDW[13,0] = 0
 udMC2: udx[14,8] - udw[13,4] >= 0
 udMC2: udx[14,8] - udw[13,8] >= 0
 udMC2: udx[14,8] - udw[13,4] - udw[13,8] <= 0
 udMC2: - uDX[14,8] + uDW[13,4] + uDW[13,8] >= 0
 udMC2: uDX[14,8] - uDW[13,4] + uDW[13,8] >= 0
 udMC2: uDX[14,8] + uDW[13,4] - uDW[13,8] >= 0
 udMC2: uDX[14,8] - uDW[13,4] - uDW[13,8] <= 0
 udMC3: udx[14,12] - udw[13,0] >= 0
 udMC3: udx[14,12] - udw[13,8] >= 0
 udMC3: udx[14,12] - udw[13,0] - udw[13,8] <= 0
 udMC3: - uDX[14,12] + uDW[13,0] + uDW[13,8] >= 0
 udMC3: uDX[14,12] - uDW[13,0] + uDW[13,8] >= 0
 udMC3: uDX[14,12] + uDW[13,0] - uDW[13,8] >= 0
 udMC3: uDX[14,12] - uDW[13,0] - uDW[13,8] <= 0
 udMC0: udx[14,1] - udw[13,13] >= 0
 udMC0: udx[14,1] - udx[14,13] >= 0
 udMC0: udx[14,1] - udx[14,13] - udw[13,13] <= 0
 udMC0: - uDX[14,1] + uDX[14,13] + uDW[13,13] >= 0
 udMC0: uDX[14,1] + uDX[14,13] - uDW[13,13] >= 0
 udMC0: uDX[14,1] - uDX[14,13] + uDW[13,13] >= 0
 udMC0: uDX[14,1] - uDX[14,13] - uDW[13,13] <= 0
 udMC1: udx[14,5] - udw[13,1] = 0
 udMC1: uDX[14,5] - uDW[13,1] = 0
 udMC2: udx[14,9] - udw[13,5] >= 0
 udMC2: udx[14,9] - udw[13,9] >= 0
 udMC2: udx[14,9] - udw[13,5] - udw[13,9] <= 0
 udMC2: - uDX[14,9] + uDW[13,5] + uDW[13,9] >= 0
 udMC2: uDX[14,9] - uDW[13,5] + uDW[13,9] >= 0
 udMC2: uDX[14,9] + uDW[13,5] - uDW[13,9] >= 0
 udMC2: uDX[14,9] - uDW[13,5] - uDW[13,9] <= 0
 udMC3: udx[14,13] - udw[13,1] >= 0
 udMC3: udx[14,13] - udw[13,9] >= 0
 udMC3: udx[14,13] - udw[13,1] - udw[13,9] <= 0
 udMC3: - uDX[14,13] + uDW[13,1] + uDW[13,9] >= 0
 udMC3: uDX[14,13] - uDW[13,1] + uDW[13,9] >= 0
 udMC3: uDX[14,13] + uDW[13,1] - uDW[13,9] >= 0
 udMC3: uDX[14,13] - uDW[13,1] - uDW[13,9] <= 0
 udMC0: udx[14,2] - udw[13,14] >= 0
 udMC0: udx[14,2] - udx[14,14] >= 0
 udMC0: udx[14,2] - udx[14,14] - udw[13,14] <= 0
 udMC0: - uDX[14,2] + uDX[14,14] + uDW[13,14] >= 0
 udMC0: uDX[14,2] + uDX[14,14] - uDW[13,14] >= 0
 udMC0: uDX[14,2] - uDX[14,14] + uDW[13,14] >= 0
 udMC0: uDX[14,2] - uDX[14,14] - uDW[13,14] <= 0
 udMC1: udx[14,6] - udw[13,2] = 0
 udMC1: uDX[14,6] - uDW[13,2] = 0
 udMC2: udx[14,10] - udw[13,6] >= 0
 udMC2: udx[14,10] - udw[13,10] >= 0
 udMC2: udx[14,10] - udw[13,6] - udw[13,10] <= 0
 udMC2: - uDX[14,10] + uDW[13,6] + uDW[13,10] >= 0
 udMC2: uDX[14,10] - uDW[13,6] + uDW[13,10] >= 0
 udMC2: uDX[14,10] + uDW[13,6] - uDW[13,10] >= 0
 udMC2: uDX[14,10] - uDW[13,6] - uDW[13,10] <= 0
 udMC3: udx[14,14] - udw[13,2] >= 0
 udMC3: udx[14,14] - udw[13,10] >= 0
 udMC3: udx[14,14] - udw[13,2] - udw[13,10] <= 0
 udMC3: - uDX[14,14] + uDW[13,2] + uDW[13,10] >= 0
 udMC3: uDX[14,14] - uDW[13,2] + uDW[13,10] >= 0
 udMC3: uDX[14,14] + uDW[13,2] - uDW[13,10] >= 0
 udMC3: uDX[14,14] - uDW[13,2] - uDW[13,10] <= 0
 udMC0: udx[14,3] - udw[13,15] >= 0
 udMC0: udx[14,3] - udx[14,15] >= 0
 udMC0: udx[14,3] - udx[14,15] - udw[13,15] <= 0
 udMC0: - uDX[14,3] + uDX[14,15] + uDW[13,15] >= 0
 udMC0: uDX[14,3] + uDX[14,15] - uDW[13,15] >= 0
 udMC0: uDX[14,3] - uDX[14,15] + uDW[13,15] >= 0
 udMC0: uDX[14,3] - uDX[14,15] - uDW[13,15] <= 0
 udMC1: udx[14,7] - udw[13,3] = 0
 udMC1: uDX[14,7] - uDW[13,3] = 0
 udMC2: udx[14,11] - udw[13,7] >= 0
 udMC2: udx[14,11] - udw[13,11] >= 0
 udMC2: udx[14,11] - udw[13,7] - udw[13,11] <= 0
 udMC2: - uDX[14,11] + uDW[13,7] + uDW[13,11] >= 0
 udMC2: uDX[14,11] - uDW[13,7] + uDW[13,11] >= 0
 udMC2: uDX[14,11] + uDW[13,7] - uDW[13,11] >= 0
 udMC2: uDX[14,11] - uDW[13,7] - uDW[13,11] <= 0
 udMC3: udx[14,15] - udw[13,3] >= 0
 udMC3: udx[14,15] - udw[13,11] >= 0
 udMC3: udx[14,15] - udw[13,3] - udw[13,11] <= 0
 udMC3: - uDX[14,15] + uDW[13,3] + uDW[13,11] >= 0
 udMC3: uDX[14,15] - uDW[13,3] + uDW[13,11] >= 0
 udMC3: uDX[14,15] + uDW[13,3] - uDW[13,11] >= 0
 udMC3: uDX[14,15] - uDW[13,3] - uDW[13,11] <= 0
 udMC0: udx[15,0] - udw[14,12] >= 0
 udMC0: udx[15,0] - udx[15,12] >= 0
 udMC0: udx[15,0] - udx[15,12] - udw[14,12] <= 0
 udMC0: - uDX[15,0] + uDX[15,12] + uDW[14,12] >= 0
 udMC0: uDX[15,0] + uDX[15,12] - uDW[14,12] >= 0
 udMC0: uDX[15,0] - uDX[15,12] + uDW[14,12] >= 0
 udMC0: uDX[15,0] - uDX[15,12] - uDW[14,12] <= 0
 udMC1: udx[15,4] - udw[14,0] = 0
 udMC1: uDX[15,4] - uDW[14,0] = 0
 udMC2: udx[15,8] - udw[14,4] >= 0
 udMC2: udx[15,8] - udw[14,8] >= 0
 udMC2: udx[15,8] - udw[14,4] - udw[14,8] <= 0
 udMC2: - uDX[15,8] + uDW[14,4] + uDW[14,8] >= 0
 udMC2: uDX[15,8] - uDW[14,4] + uDW[14,8] >= 0
 udMC2: uDX[15,8] + uDW[14,4] - uDW[14,8] >= 0
 udMC2: uDX[15,8] - uDW[14,4] - uDW[14,8] <= 0
 udMC3: udx[15,12] - udw[14,0] >= 0
 udMC3: udx[15,12] - udw[14,8] >= 0
 udMC3: udx[15,12] - udw[14,0] - udw[14,8] <= 0
 udMC3: - uDX[15,12] + uDW[14,0] + uDW[14,8] >= 0
 udMC3: uDX[15,12] - uDW[14,0] + uDW[14,8] >= 0
 udMC3: uDX[15,12] + uDW[14,0] - uDW[14,8] >= 0
 udMC3: uDX[15,12] - uDW[14,0] - uDW[14,8] <= 0
 udMC0: udx[15,1] - udw[14,13] >= 0
 udMC0: udx[15,1] - udx[15,13] >= 0
 udMC0: udx[15,1] - udx[15,13] - udw[14,13] <= 0
 udMC0: - uDX[15,1] + uDX[15,13] + uDW[14,13] >= 0
 udMC0: uDX[15,1] + uDX[15,13] - uDW[14,13] >= 0
 udMC0: uDX[15,1] - uDX[15,13] + uDW[14,13] >= 0
 udMC0: uDX[15,1] - uDX[15,13] - uDW[14,13] <= 0
 udMC1: udx[15,5] - udw[14,1] = 0
 udMC1: uDX[15,5] - uDW[14,1] = 0
 udMC2: udx[15,9] - udw[14,5] >= 0
 udMC2: udx[15,9] - udw[14,9] >= 0
 udMC2: udx[15,9] - udw[14,5] - udw[14,9] <= 0
 udMC2: - uDX[15,9] + uDW[14,5] + uDW[14,9] >= 0
 udMC2: uDX[15,9] - uDW[14,5] + uDW[14,9] >= 0
 udMC2: uDX[15,9] + uDW[14,5] - uDW[14,9] >= 0
 udMC2: uDX[15,9] - uDW[14,5] - uDW[14,9] <= 0
 udMC3: udx[15,13] - udw[14,1] >= 0
 udMC3: udx[15,13] - udw[14,9] >= 0
 udMC3: udx[15,13] - udw[14,1] - udw[14,9] <= 0
 udMC3: - uDX[15,13] + uDW[14,1] + uDW[14,9] >= 0
 udMC3: uDX[15,13] - uDW[14,1] + uDW[14,9] >= 0
 udMC3: uDX[15,13] + uDW[14,1] - uDW[14,9] >= 0
 udMC3: uDX[15,13] - uDW[14,1] - uDW[14,9] <= 0
 udMC0: udx[15,2] - udw[14,14] >= 0
 udMC0: udx[15,2] - udx[15,14] >= 0
 udMC0: udx[15,2] - udx[15,14] - udw[14,14] <= 0
 udMC0: - uDX[15,2] + uDX[15,14] + uDW[14,14] >= 0
 udMC0: uDX[15,2] + uDX[15,14] - uDW[14,14] >= 0
 udMC0: uDX[15,2] - uDX[15,14] + uDW[14,14] >= 0
 udMC0: uDX[15,2] - uDX[15,14] - uDW[14,14] <= 0
 udMC1: udx[15,6] - udw[14,2] = 0
 udMC1: uDX[15,6] - uDW[14,2] = 0
 udMC2: udx[15,10] - udw[14,6] >= 0
 udMC2: udx[15,10] - udw[14,10] >= 0
 udMC2: udx[15,10] - udw[14,6] - udw[14,10] <= 0
 udMC2: - uDX[15,10] + uDW[14,6] + uDW[14,10] >= 0
 udMC2: uDX[15,10] - uDW[14,6] + uDW[14,10] >= 0
 udMC2: uDX[15,10] + uDW[14,6] - uDW[14,10] >= 0
 udMC2: uDX[15,10] - uDW[14,6] - uDW[14,10] <= 0
 udMC3: udx[15,14] - udw[14,2] >= 0
 udMC3: udx[15,14] - udw[14,10] >= 0
 udMC3: udx[15,14] - udw[14,2] - udw[14,10] <= 0
 udMC3: - uDX[15,14] + uDW[14,2] + uDW[14,10] >= 0
 udMC3: uDX[15,14] - uDW[14,2] + uDW[14,10] >= 0
 udMC3: uDX[15,14] + uDW[14,2] - uDW[14,10] >= 0
 udMC3: uDX[15,14] - uDW[14,2] - uDW[14,10] <= 0
 udMC0: udx[15,3] - udw[14,15] >= 0
 udMC0: udx[15,3] - udx[15,15] >= 0
 udMC0: udx[15,3] - udx[15,15] - udw[14,15] <= 0
 udMC0: - uDX[15,3] + uDX[15,15] + uDW[14,15] >= 0
 udMC0: uDX[15,3] + uDX[15,15] - uDW[14,15] >= 0
 udMC0: uDX[15,3] - uDX[15,15] + uDW[14,15] >= 0
 udMC0: uDX[15,3] - uDX[15,15] - uDW[14,15] <= 0
 udMC1: udx[15,7] - udw[14,3] = 0
 udMC1: uDX[15,7] - uDW[14,3] = 0
 udMC2: udx[15,11] - udw[14,7] >= 0
 udMC2: udx[15,11] - udw[14,11] >= 0
 udMC2: udx[15,11] - udw[14,7] - udw[14,11] <= 0
 udMC2: - uDX[15,11] + uDW[14,7] + uDW[14,11] >= 0
 udMC2: uDX[15,11] - uDW[14,7] + uDW[14,11] >= 0
 udMC2: uDX[15,11] + uDW[14,7] - uDW[14,11] >= 0
 udMC2: uDX[15,11] - uDW[14,7] - uDW[14,11] <= 0
 udMC3: udx[15,15] - udw[14,3] >= 0
 udMC3: udx[15,15] - udw[14,11] >= 0
 udMC3: udx[15,15] - udw[14,3] - udw[14,11] <= 0
 udMC3: - uDX[15,15] + uDW[14,3] + uDW[14,11] >= 0
 udMC3: uDX[15,15] - uDW[14,3] + uDW[14,11] >= 0
 udMC3: uDX[15,15] + uDW[14,3] - uDW[14,11] >= 0
 udMC3: uDX[15,15] - uDW[14,3] - uDW[14,11] <= 0
 udMC0: udx[16,0] - udw[15,12] >= 0
 udMC0: udx[16,0] - udx[16,12] >= 0
 udMC0: udx[16,0] - udx[16,12] - udw[15,12] <= 0
 udMC0: - uDX[16,0] + uDX[16,12] + uDW[15,12] >= 0
 udMC0: uDX[16,0] + uDX[16,12] - uDW[15,12] >= 0
 udMC0: uDX[16,0] - uDX[16,12] + uDW[15,12] >= 0
 udMC0: uDX[16,0] - uDX[16,12] - uDW[15,12] <= 0
 udMC1: udx[16,4] - udw[15,0] = 0
 udMC1: uDX[16,4] - uDW[15,0] = 0
 udMC2: udx[16,8] - udw[15,4] >= 0
 udMC2: udx[16,8] - udw[15,8] >= 0
 udMC2: udx[16,8] - udw[15,4] - udw[15,8] <= 0
 udMC2: - uDX[16,8] + uDW[15,4] + uDW[15,8] >= 0
 udMC2: uDX[16,8] - uDW[15,4] + uDW[15,8] >= 0
 udMC2: uDX[16,8] + uDW[15,4] - uDW[15,8] >= 0
 udMC2: uDX[16,8] - uDW[15,4] - uDW[15,8] <= 0
 udMC3: udx[16,12] - udw[15,0] >= 0
 udMC3: udx[16,12] - udw[15,8] >= 0
 udMC3: udx[16,12] - udw[15,0] - udw[15,8] <= 0
 udMC3: - uDX[16,12] + uDW[15,0] + uDW[15,8] >= 0
 udMC3: uDX[16,12] - uDW[15,0] + uDW[15,8] >= 0
 udMC3: uDX[16,12] + uDW[15,0] - uDW[15,8] >= 0
 udMC3: uDX[16,12] - uDW[15,0] - uDW[15,8] <= 0
 udMC0: udx[16,1] - udw[15,13] >= 0
 udMC0: udx[16,1] - udx[16,13] >= 0
 udMC0: udx[16,1] - udx[16,13] - udw[15,13] <= 0
 udMC0: - uDX[16,1] + uDX[16,13] + uDW[15,13] >= 0
 udMC0: uDX[16,1] + uDX[16,13] - uDW[15,13] >= 0
 udMC0: uDX[16,1] - uDX[16,13] + uDW[15,13] >= 0
 udMC0: uDX[16,1] - uDX[16,13] - uDW[15,13] <= 0
 udMC1: udx[16,5] - udw[15,1] = 0
 udMC1: uDX[16,5] - uDW[15,1] = 0
 udMC2: udx[16,9] - udw[15,5] >= 0
 udMC2: udx[16,9] - udw[15,9] >= 0
 udMC2: udx[16,9] - udw[15,5] - udw[15,9] <= 0
 udMC2: - uDX[16,9] + uDW[15,5] + uDW[15,9] >= 0
 udMC2: uDX[16,9] - uDW[15,5] + uDW[15,9] >= 0
 udMC2: uDX[16,9] + uDW[15,5] - uDW[15,9] >= 0
 udMC2: uDX[16,9] - uDW[15,5] - uDW[15,9] <= 0
 udMC3: udx[16,13] - udw[15,1] >= 0
 udMC3: udx[16,13] - udw[15,9] >= 0
 udMC3: udx[16,13] - udw[15,1] - udw[15,9] <= 0
 udMC3: - uDX[16,13] + uDW[15,1] + uDW[15,9] >= 0
 udMC3: uDX[16,13] - uDW[15,1] + uDW[15,9] >= 0
 udMC3: uDX[16,13] + uDW[15,1] - uDW[15,9] >= 0
 udMC3: uDX[16,13] - uDW[15,1] - uDW[15,9] <= 0
 udMC0: udx[16,2] - udw[15,14] >= 0
 udMC0: udx[16,2] - udx[16,14] >= 0
 udMC0: udx[16,2] - udx[16,14] - udw[15,14] <= 0
 udMC0: - uDX[16,2] + uDX[16,14] + uDW[15,14] >= 0
 udMC0: uDX[16,2] + uDX[16,14] - uDW[15,14] >= 0
 udMC0: uDX[16,2] - uDX[16,14] + uDW[15,14] >= 0
 udMC0: uDX[16,2] - uDX[16,14] - uDW[15,14] <= 0
 udMC1: udx[16,6] - udw[15,2] = 0
 udMC1: uDX[16,6] - uDW[15,2] = 0
 udMC2: udx[16,10] - udw[15,6] >= 0
 udMC2: udx[16,10] - udw[15,10] >= 0
 udMC2: udx[16,10] - udw[15,6] - udw[15,10] <= 0
 udMC2: - uDX[16,10] + uDW[15,6] + uDW[15,10] >= 0
 udMC2: uDX[16,10] - uDW[15,6] + uDW[15,10] >= 0
 udMC2: uDX[16,10] + uDW[15,6] - uDW[15,10] >= 0
 udMC2: uDX[16,10] - uDW[15,6] - uDW[15,10] <= 0
 udMC3: udx[16,14] - udw[15,2] >= 0
 udMC3: udx[16,14] - udw[15,10] >= 0
 udMC3: udx[16,14] - udw[15,2] - udw[15,10] <= 0
 udMC3: - uDX[16,14] + uDW[15,2] + uDW[15,10] >= 0
 udMC3: uDX[16,14] - uDW[15,2] + uDW[15,10] >= 0
 udMC3: uDX[16,14] + uDW[15,2] - uDW[15,10] >= 0
 udMC3: uDX[16,14] - uDW[15,2] - uDW[15,10] <= 0
 udMC0: udx[16,3] - udw[15,15] >= 0
 udMC0: udx[16,3] - udx[16,15] >= 0
 udMC0: udx[16,3] - udx[16,15] - udw[15,15] <= 0
 udMC0: - uDX[16,3] + uDX[16,15] + uDW[15,15] >= 0
 udMC0: uDX[16,3] + uDX[16,15] - uDW[15,15] >= 0
 udMC0: uDX[16,3] - uDX[16,15] + uDW[15,15] >= 0
 udMC0: uDX[16,3] - uDX[16,15] - uDW[15,15] <= 0
 udMC1: udx[16,7] - udw[15,3] = 0
 udMC1: uDX[16,7] - uDW[15,3] = 0
 udMC2: udx[16,11] - udw[15,7] >= 0
 udMC2: udx[16,11] - udw[15,11] >= 0
 udMC2: udx[16,11] - udw[15,7] - udw[15,11] <= 0
 udMC2: - uDX[16,11] + uDW[15,7] + uDW[15,11] >= 0
 udMC2: uDX[16,11] - uDW[15,7] + uDW[15,11] >= 0
 udMC2: uDX[16,11] + uDW[15,7] - uDW[15,11] >= 0
 udMC2: uDX[16,11] - uDW[15,7] - uDW[15,11] <= 0
 udMC3: udx[16,15] - udw[15,3] >= 0
 udMC3: udx[16,15] - udw[15,11] >= 0
 udMC3: udx[16,15] - udw[15,3] - udw[15,11] <= 0
 udMC3: - uDX[16,15] + uDW[15,3] + uDW[15,11] >= 0
 udMC3: uDX[16,15] - uDW[15,3] + uDW[15,11] >= 0
 udMC3: uDX[16,15] + uDW[15,3] - uDW[15,11] >= 0
 udMC3: uDX[16,15] - uDW[15,3] - uDW[15,11] <= 0
 udMC0: udx[17,0] - udw[16,12] >= 0
 udMC0: udx[17,0] - udx[17,12] >= 0
 udMC0: udx[17,0] - udx[17,12] - udw[16,12] <= 0
 udMC0: - uDX[17,0] + uDX[17,12] + uDW[16,12] >= 0
 udMC0: uDX[17,0] + uDX[17,12] - uDW[16,12] >= 0
 udMC0: uDX[17,0] - uDX[17,12] + uDW[16,12] >= 0
 udMC0: uDX[17,0] - uDX[17,12] - uDW[16,12] <= 0
 udMC1: udx[17,4] - udw[16,0] = 0
 udMC1: uDX[17,4] - uDW[16,0] = 0
 udMC2: udx[17,8] - udw[16,4] >= 0
 udMC2: udx[17,8] - udw[16,8] >= 0
 udMC2: udx[17,8] - udw[16,4] - udw[16,8] <= 0
 udMC2: - uDX[17,8] + uDW[16,4] + uDW[16,8] >= 0
 udMC2: uDX[17,8] - uDW[16,4] + uDW[16,8] >= 0
 udMC2: uDX[17,8] + uDW[16,4] - uDW[16,8] >= 0
 udMC2: uDX[17,8] - uDW[16,4] - uDW[16,8] <= 0
 udMC3: udx[17,12] - udw[16,0] >= 0
 udMC3: udx[17,12] - udw[16,8] >= 0
 udMC3: udx[17,12] - udw[16,0] - udw[16,8] <= 0
 udMC3: - uDX[17,12] + uDW[16,0] + uDW[16,8] >= 0
 udMC3: uDX[17,12] - uDW[16,0] + uDW[16,8] >= 0
 udMC3: uDX[17,12] + uDW[16,0] - uDW[16,8] >= 0
 udMC3: uDX[17,12] - uDW[16,0] - uDW[16,8] <= 0
 udMC0: udx[17,1] - udw[16,13] >= 0
 udMC0: udx[17,1] - udx[17,13] >= 0
 udMC0: udx[17,1] - udx[17,13] - udw[16,13] <= 0
 udMC0: - uDX[17,1] + uDX[17,13] + uDW[16,13] >= 0
 udMC0: uDX[17,1] + uDX[17,13] - uDW[16,13] >= 0
 udMC0: uDX[17,1] - uDX[17,13] + uDW[16,13] >= 0
 udMC0: uDX[17,1] - uDX[17,13] - uDW[16,13] <= 0
 udMC1: udx[17,5] - udw[16,1] = 0
 udMC1: uDX[17,5] - uDW[16,1] = 0
 udMC2: udx[17,9] - udw[16,5] >= 0
 udMC2: udx[17,9] - udw[16,9] >= 0
 udMC2: udx[17,9] - udw[16,5] - udw[16,9] <= 0
 udMC2: - uDX[17,9] + uDW[16,5] + uDW[16,9] >= 0
 udMC2: uDX[17,9] - uDW[16,5] + uDW[16,9] >= 0
 udMC2: uDX[17,9] + uDW[16,5] - uDW[16,9] >= 0
 udMC2: uDX[17,9] - uDW[16,5] - uDW[16,9] <= 0
 udMC3: udx[17,13] - udw[16,1] >= 0
 udMC3: udx[17,13] - udw[16,9] >= 0
 udMC3: udx[17,13] - udw[16,1] - udw[16,9] <= 0
 udMC3: - uDX[17,13] + uDW[16,1] + uDW[16,9] >= 0
 udMC3: uDX[17,13] - uDW[16,1] + uDW[16,9] >= 0
 udMC3: uDX[17,13] + uDW[16,1] - uDW[16,9] >= 0
 udMC3: uDX[17,13] - uDW[16,1] - uDW[16,9] <= 0
 udMC0: udx[17,2] - udw[16,14] >= 0
 udMC0: udx[17,2] - udx[17,14] >= 0
 udMC0: udx[17,2] - udx[17,14] - udw[16,14] <= 0
 udMC0: - uDX[17,2] + uDX[17,14] + uDW[16,14] >= 0
 udMC0: uDX[17,2] + uDX[17,14] - uDW[16,14] >= 0
 udMC0: uDX[17,2] - uDX[17,14] + uDW[16,14] >= 0
 udMC0: uDX[17,2] - uDX[17,14] - uDW[16,14] <= 0
 udMC1: udx[17,6] - udw[16,2] = 0
 udMC1: uDX[17,6] - uDW[16,2] = 0
 udMC2: udx[17,10] - udw[16,6] >= 0
 udMC2: udx[17,10] - udw[16,10] >= 0
 udMC2: udx[17,10] - udw[16,6] - udw[16,10] <= 0
 udMC2: - uDX[17,10] + uDW[16,6] + uDW[16,10] >= 0
 udMC2: uDX[17,10] - uDW[16,6] + uDW[16,10] >= 0
 udMC2: uDX[17,10] + uDW[16,6] - uDW[16,10] >= 0
 udMC2: uDX[17,10] - uDW[16,6] - uDW[16,10] <= 0
 udMC3: udx[17,14] - udw[16,2] >= 0
 udMC3: udx[17,14] - udw[16,10] >= 0
 udMC3: udx[17,14] - udw[16,2] - udw[16,10] <= 0
 udMC3: - uDX[17,14] + uDW[16,2] + uDW[16,10] >= 0
 udMC3: uDX[17,14] - uDW[16,2] + uDW[16,10] >= 0
 udMC3: uDX[17,14] + uDW[16,2] - uDW[16,10] >= 0
 udMC3: uDX[17,14] - uDW[16,2] - uDW[16,10] <= 0
 udMC0: udx[17,3] - udw[16,15] >= 0
 udMC0: udx[17,3] - udx[17,15] >= 0
 udMC0: udx[17,3] - udx[17,15] - udw[16,15] <= 0
 udMC0: - uDX[17,3] + uDX[17,15] + uDW[16,15] >= 0
 udMC0: uDX[17,3] + uDX[17,15] - uDW[16,15] >= 0
 udMC0: uDX[17,3] - uDX[17,15] + uDW[16,15] >= 0
 udMC0: uDX[17,3] - uDX[17,15] - uDW[16,15] <= 0
 udMC1: udx[17,7] - udw[16,3] = 0
 udMC1: uDX[17,7] - uDW[16,3] = 0
 udMC2: udx[17,11] - udw[16,7] >= 0
 udMC2: udx[17,11] - udw[16,11] >= 0
 udMC2: udx[17,11] - udw[16,7] - udw[16,11] <= 0
 udMC2: - uDX[17,11] + uDW[16,7] + uDW[16,11] >= 0
 udMC2: uDX[17,11] - uDW[16,7] + uDW[16,11] >= 0
 udMC2: uDX[17,11] + uDW[16,7] - uDW[16,11] >= 0
 udMC2: uDX[17,11] - uDW[16,7] - uDW[16,11] <= 0
 udMC3: udx[17,15] - udw[16,3] >= 0
 udMC3: udx[17,15] - udw[16,11] >= 0
 udMC3: udx[17,15] - udw[16,3] - udw[16,11] <= 0
 udMC3: - uDX[17,15] + uDW[16,3] + uDW[16,11] >= 0
 udMC3: uDX[17,15] - uDW[16,3] + uDW[16,11] >= 0
 udMC3: uDX[17,15] + uDW[16,3] - uDW[16,11] >= 0
 udMC3: uDX[17,15] - uDW[16,3] - uDW[16,11] <= 0
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
 uKS: - ustk[0,0] - ustk[1,8] - ustk[2,2] - ustk[3,10] - ustk[4,4]
   - ustk[5,12] - ustk[6,6] - ustk[7,14] - ustk[8,5] - ustk[9,13]
   - ustk[10,3] - ustk[11,11] - ustk[12,7] - ustk[13,15] - ustk[14,1]
   - ustk[15,9] - ustk[16,0] + 17 uLANE[0] <= 3
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
 uKS: - ustk[0,1] - ustk[1,9] - ustk[2,0] - ustk[3,8] - ustk[4,2]
   - ustk[5,10] - ustk[6,4] - ustk[7,12] - ustk[8,6] - ustk[9,14]
   - ustk[10,5] - ustk[11,13] - ustk[12,3] - ustk[13,11] - ustk[14,7]
   - ustk[15,15] - ustk[16,1] + 17 uLANE[1] <= 3
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
 uKS: - ustk[0,2] - ustk[1,10] - ustk[2,4] - ustk[3,12] - ustk[4,6]
   - ustk[5,14] - ustk[6,5] - ustk[7,13] - ustk[8,3] - ustk[9,11]
   - ustk[10,7] - ustk[11,15] - ustk[12,1] - ustk[13,9] - ustk[14,0]
   - ustk[15,8] - ustk[16,2] + 17 uLANE[2] <= 3
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
 uKS: - ustk[0,3] - ustk[1,11] - ustk[2,7] - ustk[3,15] - ustk[4,1]
   - ustk[5,9] - ustk[6,0] - ustk[7,8] - ustk[8,2] - ustk[9,10]
   - ustk[10,4] - ustk[11,12] - ustk[12,6] - ustk[13,14] - ustk[14,5]
   - ustk[15,13] - ustk[16,3] + 17 uLANE[3] <= 3
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
 uKS: - ustk[0,4] - ustk[1,12] - ustk[2,6] - ustk[3,14] - ustk[4,5]
   - ustk[5,13] - ustk[6,3] - ustk[7,11] - ustk[8,7] - ustk[9,15]
   - ustk[10,1] - ustk[11,9] - ustk[12,0] - ustk[13,8] - ustk[14,2]
   - ustk[15,10] - ustk[16,4] + 17 uLANE[4] <= 3
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
 uKS: - ustk[0,5] - ustk[1,13] - ustk[2,3] - ustk[3,11] - ustk[4,7]
   - ustk[5,15] - ustk[6,1] - ustk[7,9] - ustk[8,0] - ustk[9,8]
   - ustk[10,2] - ustk[11,10] - ustk[12,4] - ustk[13,12] - ustk[14,6]
   - ustk[15,14] - ustk[16,5] + 17 uLANE[5] <= 3
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
 uKS: - ustk[0,6] - ustk[1,14] - ustk[2,5] - ustk[3,13] - ustk[4,3]
   - ustk[5,11] - ustk[6,7] - ustk[7,15] - ustk[8,1] - ustk[9,9]
   - ustk[10,0] - ustk[11,8] - ustk[12,2] - ustk[13,10] - ustk[14,4]
   - ustk[15,12] - ustk[16,6] + 17 uLANE[6] <= 3
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
 uKS: - ustk[0,7] - ustk[1,15] - ustk[2,1] - ustk[3,9] - ustk[4,0]
   - ustk[5,8] - ustk[6,2] - ustk[7,10] - ustk[8,4] - ustk[9,12]
   - ustk[10,6] - ustk[11,14] - ustk[12,5] - ustk[13,13] - ustk[14,3]
   - ustk[15,11] - ustk[16,7] + 17 uLANE[7] <= 3
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
 uKS: - ustk[0,8] - ustk[1,2] - ustk[2,10] - ustk[3,4] - ustk[4,12]
   - ustk[5,6] - ustk[6,14] - ustk[7,5] - ustk[8,13] - ustk[9,3]
   - ustk[10,11] - ustk[11,7] - ustk[12,15] - ustk[13,1] - ustk[14,9]
   - ustk[15,0] - ustk[16,8] + 17 uLANE[8] <= 3
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
 uKS: - ustk[0,9] - ustk[1,0] - ustk[2,8] - ustk[3,2] - ustk[4,10]
   - ustk[5,4] - ustk[6,12] - ustk[7,6] - ustk[8,14] - ustk[9,5]
   - ustk[10,13] - ustk[11,3] - ustk[12,11] - ustk[13,7] - ustk[14,15]
   - ustk[15,1] - ustk[16,9] + 17 uLANE[9] <= 3
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
 uKS: - ustk[0,10] - ustk[1,4] - ustk[2,12] - ustk[3,6] - ustk[4,14]
   - ustk[5,5] - ustk[6,13] - ustk[7,3] - ustk[8,11] - ustk[9,7]
   - ustk[10,15] - ustk[11,1] - ustk[12,9] - ustk[13,0] - ustk[14,8]
   - ustk[15,2] - ustk[16,10] + 17 uLANE[10] <= 3
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
 uKS: - ustk[0,11] - ustk[1,7] - ustk[2,15] - ustk[3,1] - ustk[4,9]
   - ustk[5,0] - ustk[6,8] - ustk[7,2] - ustk[8,10] - ustk[9,4]
   - ustk[10,12] - ustk[11,6] - ustk[12,14] - ustk[13,5] - ustk[14,13]
   - ustk[15,3] - ustk[16,11] + 17 uLANE[11] <= 3
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
 uKS: - ustk[0,12] - ustk[1,6] - ustk[2,14] - ustk[3,5] - ustk[4,13]
   - ustk[5,3] - ustk[6,11] - ustk[7,7] - ustk[8,15] - ustk[9,1]
   - ustk[10,9] - ustk[11,0] - ustk[12,8] - ustk[13,2] - ustk[14,10]
   - ustk[15,4] - ustk[16,12] + 17 uLANE[12] <= 3
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
 uKS: - ustk[0,13] - ustk[1,3] - ustk[2,11] - ustk[3,7] - ustk[4,15]
   - ustk[5,1] - ustk[6,9] - ustk[7,0] - ustk[8,8] - ustk[9,2]
   - ustk[10,10] - ustk[11,4] - ustk[12,12] - ustk[13,6] - ustk[14,14]
   - ustk[15,5] - ustk[16,13] + 17 uLANE[13] <= 3
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
 uKS: - ustk[0,14] - ustk[1,5] - ustk[2,13] - ustk[3,3] - ustk[4,11]
   - ustk[5,7] - ustk[6,15] - ustk[7,1] - ustk[8,9] - ustk[9,0]
   - ustk[10,8] - ustk[11,2] - ustk[12,10] - ustk[13,4] - ustk[14,12]
   - ustk[15,6] - ustk[16,14] + 17 uLANE[14] <= 3
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
 uKS: - ustk[0,15] - ustk[1,1] - ustk[2,9] - ustk[3,0] - ustk[4,8]
   - ustk[5,2] - ustk[6,10] - ustk[7,4] - ustk[8,12] - ustk[9,6]
   - ustk[10,14] - ustk[11,5] - ustk[12,13] - ustk[13,3] - ustk[14,11]
   - ustk[15,7] - ustk[16,15] + 17 uLANE[15] <= 3
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
 lBasic: lDX[24,0] - ldx[24,0] >= 0
 lBasic: lDY[24,0] - ldy[24,0] >= 0
 lBasic: lDZ[24,0] - ldz[24,0] >= 0
 lBasic: lDW[24,0] - ldw[24,0] >= 0
 lBasic: lDX[24,1] - ldx[24,1] >= 0
 lBasic: lDY[24,1] - ldy[24,1] >= 0
 lBasic: lDZ[24,1] - ldz[24,1] >= 0
 lBasic: lDW[24,1] - ldw[24,1] >= 0
 lBasic: lDX[24,2] - ldx[24,2] >= 0
 lBasic: lDY[24,2] - ldy[24,2] >= 0
 lBasic: lDZ[24,2] - ldz[24,2] >= 0
 lBasic: lDW[24,2] - ldw[24,2] >= 0
 lBasic: lDX[24,3] - ldx[24,3] >= 0
 lBasic: lDY[24,3] - ldy[24,3] >= 0
 lBasic: lDZ[24,3] - ldz[24,3] >= 0
 lBasic: lDW[24,3] - ldw[24,3] >= 0
 lBasic: lDX[24,4] - ldx[24,4] >= 0
 lBasic: lDY[24,4] - ldy[24,4] >= 0
 lBasic: lDZ[24,4] - ldz[24,4] >= 0
 lBasic: lDW[24,4] - ldw[24,4] >= 0
 lBasic: lDX[24,5] - ldx[24,5] >= 0
 lBasic: lDY[24,5] - ldy[24,5] >= 0
 lBasic: lDZ[24,5] - ldz[24,5] >= 0
 lBasic: lDW[24,5] - ldw[24,5] >= 0
 lBasic: lDX[24,6] - ldx[24,6] >= 0
 lBasic: lDY[24,6] - ldy[24,6] >= 0
 lBasic: lDZ[24,6] - ldz[24,6] >= 0
 lBasic: lDW[24,6] - ldw[24,6] >= 0
 lBasic: lDX[24,7] - ldx[24,7] >= 0
 lBasic: lDY[24,7] - ldy[24,7] >= 0
 lBasic: lDZ[24,7] - ldz[24,7] >= 0
 lBasic: lDW[24,7] - ldw[24,7] >= 0
 lBasic: lDX[24,8] - ldx[24,8] >= 0
 lBasic: lDY[24,8] - ldy[24,8] >= 0
 lBasic: lDZ[24,8] - ldz[24,8] >= 0
 lBasic: lDW[24,8] - ldw[24,8] >= 0
 lBasic: lDX[24,9] - ldx[24,9] >= 0
 lBasic: lDY[24,9] - ldy[24,9] >= 0
 lBasic: lDZ[24,9] - ldz[24,9] >= 0
 lBasic: lDW[24,9] - ldw[24,9] >= 0
 lBasic: lDX[24,10] - ldx[24,10] >= 0
 lBasic: lDY[24,10] - ldy[24,10] >= 0
 lBasic: lDZ[24,10] - ldz[24,10] >= 0
 lBasic: lDW[24,10] - ldw[24,10] >= 0
 lBasic: lDX[24,11] - ldx[24,11] >= 0
 lBasic: lDY[24,11] - ldy[24,11] >= 0
 lBasic: lDZ[24,11] - ldz[24,11] >= 0
 lBasic: lDW[24,11] - ldw[24,11] >= 0
 lBasic: lDX[24,12] - ldx[24,12] >= 0
 lBasic: lDY[24,12] - ldy[24,12] >= 0
 lBasic: lDZ[24,12] - ldz[24,12] >= 0
 lBasic: lDW[24,12] - ldw[24,12] >= 0
 lBasic: lDX[24,13] - ldx[24,13] >= 0
 lBasic: lDY[24,13] - ldy[24,13] >= 0
 lBasic: lDZ[24,13] - ldz[24,13] >= 0
 lBasic: lDW[24,13] - ldw[24,13] >= 0
 lBasic: lDX[24,14] - ldx[24,14] >= 0
 lBasic: lDY[24,14] - ldy[24,14] >= 0
 lBasic: lDZ[24,14] - ldz[24,14] >= 0
 lBasic: lDW[24,14] - ldw[24,14] >= 0
 lBasic: lDX[24,15] - ldx[24,15] >= 0
 lBasic: lDY[24,15] - ldy[24,15] >= 0
 lBasic: lDZ[24,15] - ldz[24,15] >= 0
 lBasic: lDW[24,15] - ldw[24,15] >= 0
 lBasic: lDX[25,0] - ldx[25,0] >= 0
 lBasic: lDY[25,0] - ldy[25,0] >= 0
 lBasic: lDZ[25,0] - ldz[25,0] >= 0
 lBasic: lDW[25,0] - ldw[25,0] >= 0
 lBasic: lDX[25,1] - ldx[25,1] >= 0
 lBasic: lDY[25,1] - ldy[25,1] >= 0
 lBasic: lDZ[25,1] - ldz[25,1] >= 0
 lBasic: lDW[25,1] - ldw[25,1] >= 0
 lBasic: lDX[25,2] - ldx[25,2] >= 0
 lBasic: lDY[25,2] - ldy[25,2] >= 0
 lBasic: lDZ[25,2] - ldz[25,2] >= 0
 lBasic: lDW[25,2] - ldw[25,2] >= 0
 lBasic: lDX[25,3] - ldx[25,3] >= 0
 lBasic: lDY[25,3] - ldy[25,3] >= 0
 lBasic: lDZ[25,3] - ldz[25,3] >= 0
 lBasic: lDW[25,3] - ldw[25,3] >= 0
 lBasic: lDX[25,4] - ldx[25,4] >= 0
 lBasic: lDY[25,4] - ldy[25,4] >= 0
 lBasic: lDZ[25,4] - ldz[25,4] >= 0
 lBasic: lDW[25,4] - ldw[25,4] >= 0
 lBasic: lDX[25,5] - ldx[25,5] >= 0
 lBasic: lDY[25,5] - ldy[25,5] >= 0
 lBasic: lDZ[25,5] - ldz[25,5] >= 0
 lBasic: lDW[25,5] - ldw[25,5] >= 0
 lBasic: lDX[25,6] - ldx[25,6] >= 0
 lBasic: lDY[25,6] - ldy[25,6] >= 0
 lBasic: lDZ[25,6] - ldz[25,6] >= 0
 lBasic: lDW[25,6] - ldw[25,6] >= 0
 lBasic: lDX[25,7] - ldx[25,7] >= 0
 lBasic: lDY[25,7] - ldy[25,7] >= 0
 lBasic: lDZ[25,7] - ldz[25,7] >= 0
 lBasic: lDW[25,7] - ldw[25,7] >= 0
 lBasic: lDX[25,8] - ldx[25,8] >= 0
 lBasic: lDY[25,8] - ldy[25,8] >= 0
 lBasic: lDZ[25,8] - ldz[25,8] >= 0
 lBasic: lDW[25,8] - ldw[25,8] >= 0
 lBasic: lDX[25,9] - ldx[25,9] >= 0
 lBasic: lDY[25,9] - ldy[25,9] >= 0
 lBasic: lDZ[25,9] - ldz[25,9] >= 0
 lBasic: lDW[25,9] - ldw[25,9] >= 0
 lBasic: lDX[25,10] - ldx[25,10] >= 0
 lBasic: lDY[25,10] - ldy[25,10] >= 0
 lBasic: lDZ[25,10] - ldz[25,10] >= 0
 lBasic: lDW[25,10] - ldw[25,10] >= 0
 lBasic: lDX[25,11] - ldx[25,11] >= 0
 lBasic: lDY[25,11] - ldy[25,11] >= 0
 lBasic: lDZ[25,11] - ldz[25,11] >= 0
 lBasic: lDW[25,11] - ldw[25,11] >= 0
 lBasic: lDX[25,12] - ldx[25,12] >= 0
 lBasic: lDY[25,12] - ldy[25,12] >= 0
 lBasic: lDZ[25,12] - ldz[25,12] >= 0
 lBasic: lDW[25,12] - ldw[25,12] >= 0
 lBasic: lDX[25,13] - ldx[25,13] >= 0
 lBasic: lDY[25,13] - ldy[25,13] >= 0
 lBasic: lDZ[25,13] - ldz[25,13] >= 0
 lBasic: lDW[25,13] - ldw[25,13] >= 0
 lBasic: lDX[25,14] - ldx[25,14] >= 0
 lBasic: lDY[25,14] - ldy[25,14] >= 0
 lBasic: lDZ[25,14] - ldz[25,14] >= 0
 lBasic: lDW[25,14] - ldw[25,14] >= 0
 lBasic: lDX[25,15] - ldx[25,15] >= 0
 lBasic: lDY[25,15] - ldy[25,15] >= 0
 lBasic: lDZ[25,15] - ldz[25,15] >= 0
 lBasic: lDW[25,15] - ldw[25,15] >= 0
 lBasic: lDX[26,0] - ldx[26,0] >= 0
 lBasic: lDY[26,0] - ldy[26,0] >= 0
 lBasic: lDZ[26,0] - ldz[26,0] >= 0
 lBasic: lDW[26,0] - ldw[26,0] >= 0
 lBasic: lDX[26,1] - ldx[26,1] >= 0
 lBasic: lDY[26,1] - ldy[26,1] >= 0
 lBasic: lDZ[26,1] - ldz[26,1] >= 0
 lBasic: lDW[26,1] - ldw[26,1] >= 0
 lBasic: lDX[26,2] - ldx[26,2] >= 0
 lBasic: lDY[26,2] - ldy[26,2] >= 0
 lBasic: lDZ[26,2] - ldz[26,2] >= 0
 lBasic: lDW[26,2] - ldw[26,2] >= 0
 lBasic: lDX[26,3] - ldx[26,3] >= 0
 lBasic: lDY[26,3] - ldy[26,3] >= 0
 lBasic: lDZ[26,3] - ldz[26,3] >= 0
 lBasic: lDW[26,3] - ldw[26,3] >= 0
 lBasic: lDX[26,4] - ldx[26,4] >= 0
 lBasic: lDY[26,4] - ldy[26,4] >= 0
 lBasic: lDZ[26,4] - ldz[26,4] >= 0
 lBasic: lDW[26,4] - ldw[26,4] >= 0
 lBasic: lDX[26,5] - ldx[26,5] >= 0
 lBasic: lDY[26,5] - ldy[26,5] >= 0
 lBasic: lDZ[26,5] - ldz[26,5] >= 0
 lBasic: lDW[26,5] - ldw[26,5] >= 0
 lBasic: lDX[26,6] - ldx[26,6] >= 0
 lBasic: lDY[26,6] - ldy[26,6] >= 0
 lBasic: lDZ[26,6] - ldz[26,6] >= 0
 lBasic: lDW[26,6] - ldw[26,6] >= 0
 lBasic: lDX[26,7] - ldx[26,7] >= 0
 lBasic: lDY[26,7] - ldy[26,7] >= 0
 lBasic: lDZ[26,7] - ldz[26,7] >= 0
 lBasic: lDW[26,7] - ldw[26,7] >= 0
 lBasic: lDX[26,8] - ldx[26,8] >= 0
 lBasic: lDY[26,8] - ldy[26,8] >= 0
 lBasic: lDZ[26,8] - ldz[26,8] >= 0
 lBasic: lDW[26,8] - ldw[26,8] >= 0
 lBasic: lDX[26,9] - ldx[26,9] >= 0
 lBasic: lDY[26,9] - ldy[26,9] >= 0
 lBasic: lDZ[26,9] - ldz[26,9] >= 0
 lBasic: lDW[26,9] - ldw[26,9] >= 0
 lBasic: lDX[26,10] - ldx[26,10] >= 0
 lBasic: lDY[26,10] - ldy[26,10] >= 0
 lBasic: lDZ[26,10] - ldz[26,10] >= 0
 lBasic: lDW[26,10] - ldw[26,10] >= 0
 lBasic: lDX[26,11] - ldx[26,11] >= 0
 lBasic: lDY[26,11] - ldy[26,11] >= 0
 lBasic: lDZ[26,11] - ldz[26,11] >= 0
 lBasic: lDW[26,11] - ldw[26,11] >= 0
 lBasic: lDX[26,12] - ldx[26,12] >= 0
 lBasic: lDY[26,12] - ldy[26,12] >= 0
 lBasic: lDZ[26,12] - ldz[26,12] >= 0
 lBasic: lDW[26,12] - ldw[26,12] >= 0
 lBasic: lDX[26,13] - ldx[26,13] >= 0
 lBasic: lDY[26,13] - ldy[26,13] >= 0
 lBasic: lDZ[26,13] - ldz[26,13] >= 0
 lBasic: lDW[26,13] - ldw[26,13] >= 0
 lBasic: lDX[26,14] - ldx[26,14] >= 0
 lBasic: lDY[26,14] - ldy[26,14] >= 0
 lBasic: lDZ[26,14] - ldz[26,14] >= 0
 lBasic: lDW[26,14] - ldw[26,14] >= 0
 lBasic: lDX[26,15] - ldx[26,15] >= 0
 lBasic: lDY[26,15] - ldy[26,15] >= 0
 lBasic: lDZ[26,15] - ldz[26,15] >= 0
 lBasic: lDW[26,15] - ldw[26,15] >= 0
 lBasic: lDX[27,0] - ldx[27,0] >= 0
 lBasic: lDY[27,0] - ldy[27,0] >= 0
 lBasic: lDZ[27,0] - ldz[27,0] >= 0
 lBasic: lDW[27,0] - ldw[27,0] >= 0
 lBasic: lDX[27,1] - ldx[27,1] >= 0
 lBasic: lDY[27,1] - ldy[27,1] >= 0
 lBasic: lDZ[27,1] - ldz[27,1] >= 0
 lBasic: lDW[27,1] - ldw[27,1] >= 0
 lBasic: lDX[27,2] - ldx[27,2] >= 0
 lBasic: lDY[27,2] - ldy[27,2] >= 0
 lBasic: lDZ[27,2] - ldz[27,2] >= 0
 lBasic: lDW[27,2] - ldw[27,2] >= 0
 lBasic: lDX[27,3] - ldx[27,3] >= 0
 lBasic: lDY[27,3] - ldy[27,3] >= 0
 lBasic: lDZ[27,3] - ldz[27,3] >= 0
 lBasic: lDW[27,3] - ldw[27,3] >= 0
 lBasic: lDX[27,4] - ldx[27,4] >= 0
 lBasic: lDY[27,4] - ldy[27,4] >= 0
 lBasic: lDZ[27,4] - ldz[27,4] >= 0
 lBasic: lDW[27,4] - ldw[27,4] >= 0
 lBasic: lDX[27,5] - ldx[27,5] >= 0
 lBasic: lDY[27,5] - ldy[27,5] >= 0
 lBasic: lDZ[27,5] - ldz[27,5] >= 0
 lBasic: lDW[27,5] - ldw[27,5] >= 0
 lBasic: lDX[27,6] - ldx[27,6] >= 0
 lBasic: lDY[27,6] - ldy[27,6] >= 0
 lBasic: lDZ[27,6] - ldz[27,6] >= 0
 lBasic: lDW[27,6] - ldw[27,6] >= 0
 lBasic: lDX[27,7] - ldx[27,7] >= 0
 lBasic: lDY[27,7] - ldy[27,7] >= 0
 lBasic: lDZ[27,7] - ldz[27,7] >= 0
 lBasic: lDW[27,7] - ldw[27,7] >= 0
 lBasic: lDX[27,8] - ldx[27,8] >= 0
 lBasic: lDY[27,8] - ldy[27,8] >= 0
 lBasic: lDZ[27,8] - ldz[27,8] >= 0
 lBasic: lDW[27,8] - ldw[27,8] >= 0
 lBasic: lDX[27,9] - ldx[27,9] >= 0
 lBasic: lDY[27,9] - ldy[27,9] >= 0
 lBasic: lDZ[27,9] - ldz[27,9] >= 0
 lBasic: lDW[27,9] - ldw[27,9] >= 0
 lBasic: lDX[27,10] - ldx[27,10] >= 0
 lBasic: lDY[27,10] - ldy[27,10] >= 0
 lBasic: lDZ[27,10] - ldz[27,10] >= 0
 lBasic: lDW[27,10] - ldw[27,10] >= 0
 lBasic: lDX[27,11] - ldx[27,11] >= 0
 lBasic: lDY[27,11] - ldy[27,11] >= 0
 lBasic: lDZ[27,11] - ldz[27,11] >= 0
 lBasic: lDW[27,11] - ldw[27,11] >= 0
 lBasic: lDX[27,12] - ldx[27,12] >= 0
 lBasic: lDY[27,12] - ldy[27,12] >= 0
 lBasic: lDZ[27,12] - ldz[27,12] >= 0
 lBasic: lDW[27,12] - ldw[27,12] >= 0
 lBasic: lDX[27,13] - ldx[27,13] >= 0
 lBasic: lDY[27,13] - ldy[27,13] >= 0
 lBasic: lDZ[27,13] - ldz[27,13] >= 0
 lBasic: lDW[27,13] - ldw[27,13] >= 0
 lBasic: lDX[27,14] - ldx[27,14] >= 0
 lBasic: lDY[27,14] - ldy[27,14] >= 0
 lBasic: lDZ[27,14] - ldz[27,14] >= 0
 lBasic: lDW[27,14] - ldw[27,14] >= 0
 lBasic: lDX[27,15] - ldx[27,15] >= 0
 lBasic: lDY[27,15] - ldy[27,15] >= 0
 lBasic: lDZ[27,15] - ldz[27,15] >= 0
 lBasic: lDW[27,15] - ldw[27,15] >= 0
 lBasic: lDX[28,0] - ldx[28,0] >= 0
 lBasic: lDY[28,0] - ldy[28,0] >= 0
 lBasic: lDZ[28,0] - ldz[28,0] >= 0
 lBasic: lDW[28,0] - ldw[28,0] >= 0
 lBasic: lDX[28,1] - ldx[28,1] >= 0
 lBasic: lDY[28,1] - ldy[28,1] >= 0
 lBasic: lDZ[28,1] - ldz[28,1] >= 0
 lBasic: lDW[28,1] - ldw[28,1] >= 0
 lBasic: lDX[28,2] - ldx[28,2] >= 0
 lBasic: lDY[28,2] - ldy[28,2] >= 0
 lBasic: lDZ[28,2] - ldz[28,2] >= 0
 lBasic: lDW[28,2] - ldw[28,2] >= 0
 lBasic: lDX[28,3] - ldx[28,3] >= 0
 lBasic: lDY[28,3] - ldy[28,3] >= 0
 lBasic: lDZ[28,3] - ldz[28,3] >= 0
 lBasic: lDW[28,3] - ldw[28,3] >= 0
 lBasic: lDX[28,4] - ldx[28,4] >= 0
 lBasic: lDY[28,4] - ldy[28,4] >= 0
 lBasic: lDZ[28,4] - ldz[28,4] >= 0
 lBasic: lDW[28,4] - ldw[28,4] >= 0
 lBasic: lDX[28,5] - ldx[28,5] >= 0
 lBasic: lDY[28,5] - ldy[28,5] >= 0
 lBasic: lDZ[28,5] - ldz[28,5] >= 0
 lBasic: lDW[28,5] - ldw[28,5] >= 0
 lBasic: lDX[28,6] - ldx[28,6] >= 0
 lBasic: lDY[28,6] - ldy[28,6] >= 0
 lBasic: lDZ[28,6] - ldz[28,6] >= 0
 lBasic: lDW[28,6] - ldw[28,6] >= 0
 lBasic: lDX[28,7] - ldx[28,7] >= 0
 lBasic: lDY[28,7] - ldy[28,7] >= 0
 lBasic: lDZ[28,7] - ldz[28,7] >= 0
 lBasic: lDW[28,7] - ldw[28,7] >= 0
 lBasic: lDX[28,8] - ldx[28,8] >= 0
 lBasic: lDY[28,8] - ldy[28,8] >= 0
 lBasic: lDZ[28,8] - ldz[28,8] >= 0
 lBasic: lDW[28,8] - ldw[28,8] >= 0
 lBasic: lDX[28,9] - ldx[28,9] >= 0
 lBasic: lDY[28,9] - ldy[28,9] >= 0
 lBasic: lDZ[28,9] - ldz[28,9] >= 0
 lBasic: lDW[28,9] - ldw[28,9] >= 0
 lBasic: lDX[28,10] - ldx[28,10] >= 0
 lBasic: lDY[28,10] - ldy[28,10] >= 0
 lBasic: lDZ[28,10] - ldz[28,10] >= 0
 lBasic: lDW[28,10] - ldw[28,10] >= 0
 lBasic: lDX[28,11] - ldx[28,11] >= 0
 lBasic: lDY[28,11] - ldy[28,11] >= 0
 lBasic: lDZ[28,11] - ldz[28,11] >= 0
 lBasic: lDW[28,11] - ldw[28,11] >= 0
 lBasic: lDX[28,12] - ldx[28,12] >= 0
 lBasic: lDY[28,12] - ldy[28,12] >= 0
 lBasic: lDZ[28,12] - ldz[28,12] >= 0
 lBasic: lDW[28,12] - ldw[28,12] >= 0
 lBasic: lDX[28,13] - ldx[28,13] >= 0
 lBasic: lDY[28,13] - ldy[28,13] >= 0
 lBasic: lDZ[28,13] - ldz[28,13] >= 0
 lBasic: lDW[28,13] - ldw[28,13] >= 0
 lBasic: lDX[28,14] - ldx[28,14] >= 0
 lBasic: lDY[28,14] - ldy[28,14] >= 0
 lBasic: lDZ[28,14] - ldz[28,14] >= 0
 lBasic: lDW[28,14] - ldw[28,14] >= 0
 lBasic: lDX[28,15] - ldx[28,15] >= 0
 lBasic: lDY[28,15] - ldy[28,15] >= 0
 lBasic: lDZ[28,15] - ldz[28,15] >= 0
 lBasic: lDW[28,15] - ldw[28,15] >= 0
 lBasic: lDX[29,0] - ldx[29,0] >= 0
 lBasic: lDY[29,0] - ldy[29,0] >= 0
 lBasic: lDZ[29,0] - ldz[29,0] >= 0
 lBasic: lDW[29,0] - ldw[29,0] >= 0
 lBasic: lDX[29,1] - ldx[29,1] >= 0
 lBasic: lDY[29,1] - ldy[29,1] >= 0
 lBasic: lDZ[29,1] - ldz[29,1] >= 0
 lBasic: lDW[29,1] - ldw[29,1] >= 0
 lBasic: lDX[29,2] - ldx[29,2] >= 0
 lBasic: lDY[29,2] - ldy[29,2] >= 0
 lBasic: lDZ[29,2] - ldz[29,2] >= 0
 lBasic: lDW[29,2] - ldw[29,2] >= 0
 lBasic: lDX[29,3] - ldx[29,3] >= 0
 lBasic: lDY[29,3] - ldy[29,3] >= 0
 lBasic: lDZ[29,3] - ldz[29,3] >= 0
 lBasic: lDW[29,3] - ldw[29,3] >= 0
 lBasic: lDX[29,4] - ldx[29,4] >= 0
 lBasic: lDY[29,4] - ldy[29,4] >= 0
 lBasic: lDZ[29,4] - ldz[29,4] >= 0
 lBasic: lDW[29,4] - ldw[29,4] >= 0
 lBasic: lDX[29,5] - ldx[29,5] >= 0
 lBasic: lDY[29,5] - ldy[29,5] >= 0
 lBasic: lDZ[29,5] - ldz[29,5] >= 0
 lBasic: lDW[29,5] - ldw[29,5] >= 0
 lBasic: lDX[29,6] - ldx[29,6] >= 0
 lBasic: lDY[29,6] - ldy[29,6] >= 0
 lBasic: lDZ[29,6] - ldz[29,6] >= 0
 lBasic: lDW[29,6] - ldw[29,6] >= 0
 lBasic: lDX[29,7] - ldx[29,7] >= 0
 lBasic: lDY[29,7] - ldy[29,7] >= 0
 lBasic: lDZ[29,7] - ldz[29,7] >= 0
 lBasic: lDW[29,7] - ldw[29,7] >= 0
 lBasic: lDX[29,8] - ldx[29,8] >= 0
 lBasic: lDY[29,8] - ldy[29,8] >= 0
 lBasic: lDZ[29,8] - ldz[29,8] >= 0
 lBasic: lDW[29,8] - ldw[29,8] >= 0
 lBasic: lDX[29,9] - ldx[29,9] >= 0
 lBasic: lDY[29,9] - ldy[29,9] >= 0
 lBasic: lDZ[29,9] - ldz[29,9] >= 0
 lBasic: lDW[29,9] - ldw[29,9] >= 0
 lBasic: lDX[29,10] - ldx[29,10] >= 0
 lBasic: lDY[29,10] - ldy[29,10] >= 0
 lBasic: lDZ[29,10] - ldz[29,10] >= 0
 lBasic: lDW[29,10] - ldw[29,10] >= 0
 lBasic: lDX[29,11] - ldx[29,11] >= 0
 lBasic: lDY[29,11] - ldy[29,11] >= 0
 lBasic: lDZ[29,11] - ldz[29,11] >= 0
 lBasic: lDW[29,11] - ldw[29,11] >= 0
 lBasic: lDX[29,12] - ldx[29,12] >= 0
 lBasic: lDY[29,12] - ldy[29,12] >= 0
 lBasic: lDZ[29,12] - ldz[29,12] >= 0
 lBasic: lDW[29,12] - ldw[29,12] >= 0
 lBasic: lDX[29,13] - ldx[29,13] >= 0
 lBasic: lDY[29,13] - ldy[29,13] >= 0
 lBasic: lDZ[29,13] - ldz[29,13] >= 0
 lBasic: lDW[29,13] - ldw[29,13] >= 0
 lBasic: lDX[29,14] - ldx[29,14] >= 0
 lBasic: lDY[29,14] - ldy[29,14] >= 0
 lBasic: lDZ[29,14] - ldz[29,14] >= 0
 lBasic: lDW[29,14] - ldw[29,14] >= 0
 lBasic: lDX[29,15] - ldx[29,15] >= 0
 lBasic: lDY[29,15] - ldy[29,15] >= 0
 lBasic: lDZ[29,15] - ldz[29,15] >= 0
 lBasic: lDW[29,15] - ldw[29,15] >= 0
 lBasic: lDX[30,0] - ldx[30,0] >= 0
 lBasic: lDY[30,0] - ldy[30,0] >= 0
 lBasic: lDZ[30,0] - ldz[30,0] >= 0
 lBasic: lDW[30,0] - ldw[30,0] >= 0
 lBasic: lDX[30,1] - ldx[30,1] >= 0
 lBasic: lDY[30,1] - ldy[30,1] >= 0
 lBasic: lDZ[30,1] - ldz[30,1] >= 0
 lBasic: lDW[30,1] - ldw[30,1] >= 0
 lBasic: lDX[30,2] - ldx[30,2] >= 0
 lBasic: lDY[30,2] - ldy[30,2] >= 0
 lBasic: lDZ[30,2] - ldz[30,2] >= 0
 lBasic: lDW[30,2] - ldw[30,2] >= 0
 lBasic: lDX[30,3] - ldx[30,3] >= 0
 lBasic: lDY[30,3] - ldy[30,3] >= 0
 lBasic: lDZ[30,3] - ldz[30,3] >= 0
 lBasic: lDW[30,3] - ldw[30,3] >= 0
 lBasic: lDX[30,4] - ldx[30,4] >= 0
 lBasic: lDY[30,4] - ldy[30,4] >= 0
 lBasic: lDZ[30,4] - ldz[30,4] >= 0
 lBasic: lDW[30,4] - ldw[30,4] >= 0
 lBasic: lDX[30,5] - ldx[30,5] >= 0
 lBasic: lDY[30,5] - ldy[30,5] >= 0
 lBasic: lDZ[30,5] - ldz[30,5] >= 0
 lBasic: lDW[30,5] - ldw[30,5] >= 0
 lBasic: lDX[30,6] - ldx[30,6] >= 0
 lBasic: lDY[30,6] - ldy[30,6] >= 0
 lBasic: lDZ[30,6] - ldz[30,6] >= 0
 lBasic: lDW[30,6] - ldw[30,6] >= 0
 lBasic: lDX[30,7] - ldx[30,7] >= 0
 lBasic: lDY[30,7] - ldy[30,7] >= 0
 lBasic: lDZ[30,7] - ldz[30,7] >= 0
 lBasic: lDW[30,7] - ldw[30,7] >= 0
 lBasic: lDX[30,8] - ldx[30,8] >= 0
 lBasic: lDY[30,8] - ldy[30,8] >= 0
 lBasic: lDZ[30,8] - ldz[30,8] >= 0
 lBasic: lDW[30,8] - ldw[30,8] >= 0
 lBasic: lDX[30,9] - ldx[30,9] >= 0
 lBasic: lDY[30,9] - ldy[30,9] >= 0
 lBasic: lDZ[30,9] - ldz[30,9] >= 0
 lBasic: lDW[30,9] - ldw[30,9] >= 0
 lBasic: lDX[30,10] - ldx[30,10] >= 0
 lBasic: lDY[30,10] - ldy[30,10] >= 0
 lBasic: lDZ[30,10] - ldz[30,10] >= 0
 lBasic: lDW[30,10] - ldw[30,10] >= 0
 lBasic: lDX[30,11] - ldx[30,11] >= 0
 lBasic: lDY[30,11] - ldy[30,11] >= 0
 lBasic: lDZ[30,11] - ldz[30,11] >= 0
 lBasic: lDW[30,11] - ldw[30,11] >= 0
 lBasic: lDX[30,12] - ldx[30,12] >= 0
 lBasic: lDY[30,12] - ldy[30,12] >= 0
 lBasic: lDZ[30,12] - ldz[30,12] >= 0
 lBasic: lDW[30,12] - ldw[30,12] >= 0
 lBasic: lDX[30,13] - ldx[30,13] >= 0
 lBasic: lDY[30,13] - ldy[30,13] >= 0
 lBasic: lDZ[30,13] - ldz[30,13] >= 0
 lBasic: lDW[30,13] - ldw[30,13] >= 0
 lBasic: lDX[30,14] - ldx[30,14] >= 0
 lBasic: lDY[30,14] - ldy[30,14] >= 0
 lBasic: lDZ[30,14] - ldz[30,14] >= 0
 lBasic: lDW[30,14] - ldw[30,14] >= 0
 lBasic: lDX[30,15] - ldx[30,15] >= 0
 lBasic: lDY[30,15] - ldy[30,15] >= 0
 lBasic: lDZ[30,15] - ldz[30,15] >= 0
 lBasic: lDW[30,15] - ldw[30,15] >= 0
 lBasic: lDX[31,0] - ldx[31,0] >= 0
 lBasic: lDY[31,0] - ldy[31,0] >= 0
 lBasic: lDZ[31,0] - ldz[31,0] >= 0
 lBasic: lDW[31,0] - ldw[31,0] >= 0
 lBasic: lDX[31,1] - ldx[31,1] >= 0
 lBasic: lDY[31,1] - ldy[31,1] >= 0
 lBasic: lDZ[31,1] - ldz[31,1] >= 0
 lBasic: lDW[31,1] - ldw[31,1] >= 0
 lBasic: lDX[31,2] - ldx[31,2] >= 0
 lBasic: lDY[31,2] - ldy[31,2] >= 0
 lBasic: lDZ[31,2] - ldz[31,2] >= 0
 lBasic: lDW[31,2] - ldw[31,2] >= 0
 lBasic: lDX[31,3] - ldx[31,3] >= 0
 lBasic: lDY[31,3] - ldy[31,3] >= 0
 lBasic: lDZ[31,3] - ldz[31,3] >= 0
 lBasic: lDW[31,3] - ldw[31,3] >= 0
 lBasic: lDX[31,4] - ldx[31,4] >= 0
 lBasic: lDY[31,4] - ldy[31,4] >= 0
 lBasic: lDZ[31,4] - ldz[31,4] >= 0
 lBasic: lDW[31,4] - ldw[31,4] >= 0
 lBasic: lDX[31,5] - ldx[31,5] >= 0
 lBasic: lDY[31,5] - ldy[31,5] >= 0
 lBasic: lDZ[31,5] - ldz[31,5] >= 0
 lBasic: lDW[31,5] - ldw[31,5] >= 0
 lBasic: lDX[31,6] - ldx[31,6] >= 0
 lBasic: lDY[31,6] - ldy[31,6] >= 0
 lBasic: lDZ[31,6] - ldz[31,6] >= 0
 lBasic: lDW[31,6] - ldw[31,6] >= 0
 lBasic: lDX[31,7] - ldx[31,7] >= 0
 lBasic: lDY[31,7] - ldy[31,7] >= 0
 lBasic: lDZ[31,7] - ldz[31,7] >= 0
 lBasic: lDW[31,7] - ldw[31,7] >= 0
 lBasic: lDX[31,8] - ldx[31,8] >= 0
 lBasic: lDY[31,8] - ldy[31,8] >= 0
 lBasic: lDZ[31,8] - ldz[31,8] >= 0
 lBasic: lDW[31,8] - ldw[31,8] >= 0
 lBasic: lDX[31,9] - ldx[31,9] >= 0
 lBasic: lDY[31,9] - ldy[31,9] >= 0
 lBasic: lDZ[31,9] - ldz[31,9] >= 0
 lBasic: lDW[31,9] - ldw[31,9] >= 0
 lBasic: lDX[31,10] - ldx[31,10] >= 0
 lBasic: lDY[31,10] - ldy[31,10] >= 0
 lBasic: lDZ[31,10] - ldz[31,10] >= 0
 lBasic: lDW[31,10] - ldw[31,10] >= 0
 lBasic: lDX[31,11] - ldx[31,11] >= 0
 lBasic: lDY[31,11] - ldy[31,11] >= 0
 lBasic: lDZ[31,11] - ldz[31,11] >= 0
 lBasic: lDW[31,11] - ldw[31,11] >= 0
 lBasic: lDX[31,12] - ldx[31,12] >= 0
 lBasic: lDY[31,12] - ldy[31,12] >= 0
 lBasic: lDZ[31,12] - ldz[31,12] >= 0
 lBasic: lDW[31,12] - ldw[31,12] >= 0
 lBasic: lDX[31,13] - ldx[31,13] >= 0
 lBasic: lDY[31,13] - ldy[31,13] >= 0
 lBasic: lDZ[31,13] - ldz[31,13] >= 0
 lBasic: lDW[31,13] - ldw[31,13] >= 0
 lBasic: lDX[31,14] - ldx[31,14] >= 0
 lBasic: lDY[31,14] - ldy[31,14] >= 0
 lBasic: lDZ[31,14] - ldz[31,14] >= 0
 lBasic: lDW[31,14] - ldw[31,14] >= 0
 lBasic: lDX[31,15] - ldx[31,15] >= 0
 lBasic: lDY[31,15] - ldy[31,15] >= 0
 lBasic: lDZ[31,15] - ldz[31,15] >= 0
 lBasic: lDW[31,15] - ldw[31,15] >= 0
 lBasic: lDX[32,0] - ldx[32,0] >= 0
 lBasic: lDY[32,0] - ldy[32,0] >= 0
 lBasic: lDZ[32,0] - ldz[32,0] >= 0
 lBasic: lDW[32,0] - ldw[32,0] >= 0
 lBasic: lDX[32,1] - ldx[32,1] >= 0
 lBasic: lDY[32,1] - ldy[32,1] >= 0
 lBasic: lDZ[32,1] - ldz[32,1] >= 0
 lBasic: lDW[32,1] - ldw[32,1] >= 0
 lBasic: lDX[32,2] - ldx[32,2] >= 0
 lBasic: lDY[32,2] - ldy[32,2] >= 0
 lBasic: lDZ[32,2] - ldz[32,2] >= 0
 lBasic: lDW[32,2] - ldw[32,2] >= 0
 lBasic: lDX[32,3] - ldx[32,3] >= 0
 lBasic: lDY[32,3] - ldy[32,3] >= 0
 lBasic: lDZ[32,3] - ldz[32,3] >= 0
 lBasic: lDW[32,3] - ldw[32,3] >= 0
 lBasic: lDX[32,4] - ldx[32,4] >= 0
 lBasic: lDY[32,4] - ldy[32,4] >= 0
 lBasic: lDZ[32,4] - ldz[32,4] >= 0
 lBasic: lDW[32,4] - ldw[32,4] >= 0
 lBasic: lDX[32,5] - ldx[32,5] >= 0
 lBasic: lDY[32,5] - ldy[32,5] >= 0
 lBasic: lDZ[32,5] - ldz[32,5] >= 0
 lBasic: lDW[32,5] - ldw[32,5] >= 0
 lBasic: lDX[32,6] - ldx[32,6] >= 0
 lBasic: lDY[32,6] - ldy[32,6] >= 0
 lBasic: lDZ[32,6] - ldz[32,6] >= 0
 lBasic: lDW[32,6] - ldw[32,6] >= 0
 lBasic: lDX[32,7] - ldx[32,7] >= 0
 lBasic: lDY[32,7] - ldy[32,7] >= 0
 lBasic: lDZ[32,7] - ldz[32,7] >= 0
 lBasic: lDW[32,7] - ldw[32,7] >= 0
 lBasic: lDX[32,8] - ldx[32,8] >= 0
 lBasic: lDY[32,8] - ldy[32,8] >= 0
 lBasic: lDZ[32,8] - ldz[32,8] >= 0
 lBasic: lDW[32,8] - ldw[32,8] >= 0
 lBasic: lDX[32,9] - ldx[32,9] >= 0
 lBasic: lDY[32,9] - ldy[32,9] >= 0
 lBasic: lDZ[32,9] - ldz[32,9] >= 0
 lBasic: lDW[32,9] - ldw[32,9] >= 0
 lBasic: lDX[32,10] - ldx[32,10] >= 0
 lBasic: lDY[32,10] - ldy[32,10] >= 0
 lBasic: lDZ[32,10] - ldz[32,10] >= 0
 lBasic: lDW[32,10] - ldw[32,10] >= 0
 lBasic: lDX[32,11] - ldx[32,11] >= 0
 lBasic: lDY[32,11] - ldy[32,11] >= 0
 lBasic: lDZ[32,11] - ldz[32,11] >= 0
 lBasic: lDW[32,11] - ldw[32,11] >= 0
 lBasic: lDX[32,12] - ldx[32,12] >= 0
 lBasic: lDY[32,12] - ldy[32,12] >= 0
 lBasic: lDZ[32,12] - ldz[32,12] >= 0
 lBasic: lDW[32,12] - ldw[32,12] >= 0
 lBasic: lDX[32,13] - ldx[32,13] >= 0
 lBasic: lDY[32,13] - ldy[32,13] >= 0
 lBasic: lDZ[32,13] - ldz[32,13] >= 0
 lBasic: lDW[32,13] - ldw[32,13] >= 0
 lBasic: lDX[32,14] - ldx[32,14] >= 0
 lBasic: lDY[32,14] - ldy[32,14] >= 0
 lBasic: lDZ[32,14] - ldz[32,14] >= 0
 lBasic: lDW[32,14] - ldw[32,14] >= 0
 lBasic: lDX[32,15] - ldx[32,15] >= 0
 lBasic: lDY[32,15] - ldy[32,15] >= 0
 lBasic: lDZ[32,15] - ldz[32,15] >= 0
 lBasic: lDW[32,15] - ldw[32,15] >= 0
 lBasic: lDY[17,0] - ldy[17,0] >= 0
 lBasic: lDZ[17,0] - ldz[17,0] >= 0
 lBasic: lDW[17,0] - ldw[17,0] >= 0
 lBasic: lDY[17,1] - ldy[17,1] >= 0
 lBasic: lDZ[17,1] - ldz[17,1] >= 0
 lBasic: lDW[17,1] - ldw[17,1] >= 0
 lBasic: lDY[17,2] - ldy[17,2] >= 0
 lBasic: lDZ[17,2] - ldz[17,2] >= 0
 lBasic: lDW[17,2] - ldw[17,2] >= 0
 lBasic: lDY[17,3] - ldy[17,3] >= 0
 lBasic: lDZ[17,3] - ldz[17,3] >= 0
 lBasic: lDW[17,3] - ldw[17,3] >= 0
 lBasic: lDY[17,4] - ldy[17,4] >= 0
 lBasic: lDZ[17,4] - ldz[17,4] >= 0
 lBasic: lDW[17,4] - ldw[17,4] >= 0
 lBasic: lDY[17,5] - ldy[17,5] >= 0
 lBasic: lDZ[17,5] - ldz[17,5] >= 0
 lBasic: lDW[17,5] - ldw[17,5] >= 0
 lBasic: lDY[17,6] - ldy[17,6] >= 0
 lBasic: lDZ[17,6] - ldz[17,6] >= 0
 lBasic: lDW[17,6] - ldw[17,6] >= 0
 lBasic: lDY[17,7] - ldy[17,7] >= 0
 lBasic: lDZ[17,7] - ldz[17,7] >= 0
 lBasic: lDW[17,7] - ldw[17,7] >= 0
 lBasic: lDY[17,8] - ldy[17,8] >= 0
 lBasic: lDZ[17,8] - ldz[17,8] >= 0
 lBasic: lDW[17,8] - ldw[17,8] >= 0
 lBasic: lDY[17,9] - ldy[17,9] >= 0
 lBasic: lDZ[17,9] - ldz[17,9] >= 0
 lBasic: lDW[17,9] - ldw[17,9] >= 0
 lBasic: lDY[17,10] - ldy[17,10] >= 0
 lBasic: lDZ[17,10] - ldz[17,10] >= 0
 lBasic: lDW[17,10] - ldw[17,10] >= 0
 lBasic: lDY[17,11] - ldy[17,11] >= 0
 lBasic: lDZ[17,11] - ldz[17,11] >= 0
 lBasic: lDW[17,11] - ldw[17,11] >= 0
 lBasic: lDY[17,12] - ldy[17,12] >= 0
 lBasic: lDZ[17,12] - ldz[17,12] >= 0
 lBasic: lDW[17,12] - ldw[17,12] >= 0
 lBasic: lDY[17,13] - ldy[17,13] >= 0
 lBasic: lDZ[17,13] - ldz[17,13] >= 0
 lBasic: lDW[17,13] - ldw[17,13] >= 0
 lBasic: lDY[17,14] - ldy[17,14] >= 0
 lBasic: lDZ[17,14] - ldz[17,14] >= 0
 lBasic: lDW[17,14] - ldw[17,14] >= 0
 lBasic: lDY[17,15] - ldy[17,15] >= 0
 lBasic: lDZ[17,15] - ldz[17,15] >= 0
 lBasic: lDW[17,15] - ldw[17,15] >= 0
 lBasic: lDX[33,0] - ldx[33,0] >= 0
 lBasic: lDX[33,1] - ldx[33,1] >= 0
 lBasic: lDX[33,2] - ldx[33,2] >= 0
 lBasic: lDX[33,3] - ldx[33,3] >= 0
 lBasic: lDX[33,4] - ldx[33,4] >= 0
 lBasic: lDX[33,5] - ldx[33,5] >= 0
 lBasic: lDX[33,6] - ldx[33,6] >= 0
 lBasic: lDX[33,7] - ldx[33,7] >= 0
 lBasic: lDX[33,8] - ldx[33,8] >= 0
 lBasic: lDX[33,9] - ldx[33,9] >= 0
 lBasic: lDX[33,10] - ldx[33,10] >= 0
 lBasic: lDX[33,11] - ldx[33,11] >= 0
 lBasic: lDX[33,12] - ldx[33,12] >= 0
 lBasic: lDX[33,13] - ldx[33,13] >= 0
 lBasic: lDX[33,14] - ldx[33,14] >= 0
 lBasic: lDX[33,15] - ldx[33,15] >= 0
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
 ldSC: lDX[19,0] - lDY[19,0] = 0
 ldSC: ldx[19,0] - lDY[19,0] = 0
 ldSC: lDX[19,1] - lDY[19,1] = 0
 ldSC: ldx[19,1] - lDY[19,1] = 0
 ldSC: lDX[19,2] - lDY[19,2] = 0
 ldSC: ldx[19,2] - lDY[19,2] = 0
 ldSC: lDX[19,3] - lDY[19,3] = 0
 ldSC: ldx[19,3] - lDY[19,3] = 0
 ldSC: lDX[19,4] - lDY[19,4] = 0
 ldSC: ldx[19,4] - lDY[19,4] = 0
 ldSC: lDX[19,5] - lDY[19,5] = 0
 ldSC: ldx[19,5] - lDY[19,5] = 0
 ldSC: lDX[19,6] - lDY[19,6] = 0
 ldSC: ldx[19,6] - lDY[19,6] = 0
 ldSC: lDX[19,7] - lDY[19,7] = 0
 ldSC: ldx[19,7] - lDY[19,7] = 0
 ldSC: lDX[19,8] - lDY[19,8] = 0
 ldSC: ldx[19,8] - lDY[19,8] = 0
 ldSC: lDX[19,9] - lDY[19,9] = 0
 ldSC: ldx[19,9] - lDY[19,9] = 0
 ldSC: lDX[19,10] - lDY[19,10] = 0
 ldSC: ldx[19,10] - lDY[19,10] = 0
 ldSC: lDX[19,11] - lDY[19,11] = 0
 ldSC: ldx[19,11] - lDY[19,11] = 0
 ldSC: lDX[19,12] - lDY[19,12] = 0
 ldSC: ldx[19,12] - lDY[19,12] = 0
 ldSC: lDX[19,13] - lDY[19,13] = 0
 ldSC: ldx[19,13] - lDY[19,13] = 0
 ldSC: lDX[19,14] - lDY[19,14] = 0
 ldSC: ldx[19,14] - lDY[19,14] = 0
 ldSC: lDX[19,15] - lDY[19,15] = 0
 ldSC: ldx[19,15] - lDY[19,15] = 0
 ldSC: lDX[20,0] - lDY[20,0] = 0
 ldSC: ldx[20,0] - lDY[20,0] = 0
 ldSC: lDX[20,1] - lDY[20,1] = 0
 ldSC: ldx[20,1] - lDY[20,1] = 0
 ldSC: lDX[20,2] - lDY[20,2] = 0
 ldSC: ldx[20,2] - lDY[20,2] = 0
 ldSC: lDX[20,3] - lDY[20,3] = 0
 ldSC: ldx[20,3] - lDY[20,3] = 0
 ldSC: lDX[20,4] - lDY[20,4] = 0
 ldSC: ldx[20,4] - lDY[20,4] = 0
 ldSC: lDX[20,5] - lDY[20,5] = 0
 ldSC: ldx[20,5] - lDY[20,5] = 0
 ldSC: lDX[20,6] - lDY[20,6] = 0
 ldSC: ldx[20,6] - lDY[20,6] = 0
 ldSC: lDX[20,7] - lDY[20,7] = 0
 ldSC: ldx[20,7] - lDY[20,7] = 0
 ldSC: lDX[20,8] - lDY[20,8] = 0
 ldSC: ldx[20,8] - lDY[20,8] = 0
 ldSC: lDX[20,9] - lDY[20,9] = 0
 ldSC: ldx[20,9] - lDY[20,9] = 0
 ldSC: lDX[20,10] - lDY[20,10] = 0
 ldSC: ldx[20,10] - lDY[20,10] = 0
 ldSC: lDX[20,11] - lDY[20,11] = 0
 ldSC: ldx[20,11] - lDY[20,11] = 0
 ldSC: lDX[20,12] - lDY[20,12] = 0
 ldSC: ldx[20,12] - lDY[20,12] = 0
 ldSC: lDX[20,13] - lDY[20,13] = 0
 ldSC: ldx[20,13] - lDY[20,13] = 0
 ldSC: lDX[20,14] - lDY[20,14] = 0
 ldSC: ldx[20,14] - lDY[20,14] = 0
 ldSC: lDX[20,15] - lDY[20,15] = 0
 ldSC: ldx[20,15] - lDY[20,15] = 0
 ldSC: lDX[21,0] - lDY[21,0] = 0
 ldSC: ldx[21,0] - lDY[21,0] = 0
 ldSC: lDX[21,1] - lDY[21,1] = 0
 ldSC: ldx[21,1] - lDY[21,1] = 0
 ldSC: lDX[21,2] - lDY[21,2] = 0
 ldSC: ldx[21,2] - lDY[21,2] = 0
 ldSC: lDX[21,3] - lDY[21,3] = 0
 ldSC: ldx[21,3] - lDY[21,3] = 0
 ldSC: lDX[21,4] - lDY[21,4] = 0
 ldSC: ldx[21,4] - lDY[21,4] = 0
 ldSC: lDX[21,5] - lDY[21,5] = 0
 ldSC: ldx[21,5] - lDY[21,5] = 0
 ldSC: lDX[21,6] - lDY[21,6] = 0
 ldSC: ldx[21,6] - lDY[21,6] = 0
 ldSC: lDX[21,7] - lDY[21,7] = 0
 ldSC: ldx[21,7] - lDY[21,7] = 0
 ldSC: lDX[21,8] - lDY[21,8] = 0
 ldSC: ldx[21,8] - lDY[21,8] = 0
 ldSC: lDX[21,9] - lDY[21,9] = 0
 ldSC: ldx[21,9] - lDY[21,9] = 0
 ldSC: lDX[21,10] - lDY[21,10] = 0
 ldSC: ldx[21,10] - lDY[21,10] = 0
 ldSC: lDX[21,11] - lDY[21,11] = 0
 ldSC: ldx[21,11] - lDY[21,11] = 0
 ldSC: lDX[21,12] - lDY[21,12] = 0
 ldSC: ldx[21,12] - lDY[21,12] = 0
 ldSC: lDX[21,13] - lDY[21,13] = 0
 ldSC: ldx[21,13] - lDY[21,13] = 0
 ldSC: lDX[21,14] - lDY[21,14] = 0
 ldSC: ldx[21,14] - lDY[21,14] = 0
 ldSC: lDX[21,15] - lDY[21,15] = 0
 ldSC: ldx[21,15] - lDY[21,15] = 0
 ldSC: lDX[22,0] - lDY[22,0] = 0
 ldSC: ldx[22,0] - lDY[22,0] = 0
 ldSC: lDX[22,1] - lDY[22,1] = 0
 ldSC: ldx[22,1] - lDY[22,1] = 0
 ldSC: lDX[22,2] - lDY[22,2] = 0
 ldSC: ldx[22,2] - lDY[22,2] = 0
 ldSC: lDX[22,3] - lDY[22,3] = 0
 ldSC: ldx[22,3] - lDY[22,3] = 0
 ldSC: lDX[22,4] - lDY[22,4] = 0
 ldSC: ldx[22,4] - lDY[22,4] = 0
 ldSC: lDX[22,5] - lDY[22,5] = 0
 ldSC: ldx[22,5] - lDY[22,5] = 0
 ldSC: lDX[22,6] - lDY[22,6] = 0
 ldSC: ldx[22,6] - lDY[22,6] = 0
 ldSC: lDX[22,7] - lDY[22,7] = 0
 ldSC: ldx[22,7] - lDY[22,7] = 0
 ldSC: lDX[22,8] - lDY[22,8] = 0
 ldSC: ldx[22,8] - lDY[22,8] = 0
 ldSC: lDX[22,9] - lDY[22,9] = 0
 ldSC: ldx[22,9] - lDY[22,9] = 0
 ldSC: lDX[22,10] - lDY[22,10] = 0
 ldSC: ldx[22,10] - lDY[22,10] = 0
 ldSC: lDX[22,11] - lDY[22,11] = 0
 ldSC: ldx[22,11] - lDY[22,11] = 0
 ldSC: lDX[22,12] - lDY[22,12] = 0
 ldSC: ldx[22,12] - lDY[22,12] = 0
 ldSC: lDX[22,13] - lDY[22,13] = 0
 ldSC: ldx[22,13] - lDY[22,13] = 0
 ldSC: lDX[22,14] - lDY[22,14] = 0
 ldSC: ldx[22,14] - lDY[22,14] = 0
 ldSC: lDX[22,15] - lDY[22,15] = 0
 ldSC: ldx[22,15] - lDY[22,15] = 0
 ldSC: lDX[23,0] - lDY[23,0] = 0
 ldSC: ldx[23,0] - lDY[23,0] = 0
 ldSC: lDX[23,1] - lDY[23,1] = 0
 ldSC: ldx[23,1] - lDY[23,1] = 0
 ldSC: lDX[23,2] - lDY[23,2] = 0
 ldSC: ldx[23,2] - lDY[23,2] = 0
 ldSC: lDX[23,3] - lDY[23,3] = 0
 ldSC: ldx[23,3] - lDY[23,3] = 0
 ldSC: lDX[23,4] - lDY[23,4] = 0
 ldSC: ldx[23,4] - lDY[23,4] = 0
 ldSC: lDX[23,5] - lDY[23,5] = 0
 ldSC: ldx[23,5] - lDY[23,5] = 0
 ldSC: lDX[23,6] - lDY[23,6] = 0
 ldSC: ldx[23,6] - lDY[23,6] = 0
 ldSC: lDX[23,7] - lDY[23,7] = 0
 ldSC: ldx[23,7] - lDY[23,7] = 0
 ldSC: lDX[23,8] - lDY[23,8] = 0
 ldSC: ldx[23,8] - lDY[23,8] = 0
 ldSC: lDX[23,9] - lDY[23,9] = 0
 ldSC: ldx[23,9] - lDY[23,9] = 0
 ldSC: lDX[23,10] - lDY[23,10] = 0
 ldSC: ldx[23,10] - lDY[23,10] = 0
 ldSC: lDX[23,11] - lDY[23,11] = 0
 ldSC: ldx[23,11] - lDY[23,11] = 0
 ldSC: lDX[23,12] - lDY[23,12] = 0
 ldSC: ldx[23,12] - lDY[23,12] = 0
 ldSC: lDX[23,13] - lDY[23,13] = 0
 ldSC: ldx[23,13] - lDY[23,13] = 0
 ldSC: lDX[23,14] - lDY[23,14] = 0
 ldSC: ldx[23,14] - lDY[23,14] = 0
 ldSC: lDX[23,15] - lDY[23,15] = 0
 ldSC: ldx[23,15] - lDY[23,15] = 0
 ldSC: lDX[24,0] - lDY[24,0] = 0
 ldSC: ldx[24,0] - lDY[24,0] = 0
 ldSC: lDX[24,1] - lDY[24,1] = 0
 ldSC: ldx[24,1] - lDY[24,1] = 0
 ldSC: lDX[24,2] - lDY[24,2] = 0
 ldSC: ldx[24,2] - lDY[24,2] = 0
 ldSC: lDX[24,3] - lDY[24,3] = 0
 ldSC: ldx[24,3] - lDY[24,3] = 0
 ldSC: lDX[24,4] - lDY[24,4] = 0
 ldSC: ldx[24,4] - lDY[24,4] = 0
 ldSC: lDX[24,5] - lDY[24,5] = 0
 ldSC: ldx[24,5] - lDY[24,5] = 0
 ldSC: lDX[24,6] - lDY[24,6] = 0
 ldSC: ldx[24,6] - lDY[24,6] = 0
 ldSC: lDX[24,7] - lDY[24,7] = 0
 ldSC: ldx[24,7] - lDY[24,7] = 0
 ldSC: lDX[24,8] - lDY[24,8] = 0
 ldSC: ldx[24,8] - lDY[24,8] = 0
 ldSC: lDX[24,9] - lDY[24,9] = 0
 ldSC: ldx[24,9] - lDY[24,9] = 0
 ldSC: lDX[24,10] - lDY[24,10] = 0
 ldSC: ldx[24,10] - lDY[24,10] = 0
 ldSC: lDX[24,11] - lDY[24,11] = 0
 ldSC: ldx[24,11] - lDY[24,11] = 0
 ldSC: lDX[24,12] - lDY[24,12] = 0
 ldSC: ldx[24,12] - lDY[24,12] = 0
 ldSC: lDX[24,13] - lDY[24,13] = 0
 ldSC: ldx[24,13] - lDY[24,13] = 0
 ldSC: lDX[24,14] - lDY[24,14] = 0
 ldSC: ldx[24,14] - lDY[24,14] = 0
 ldSC: lDX[24,15] - lDY[24,15] = 0
 ldSC: ldx[24,15] - lDY[24,15] = 0
 ldSC: lDX[25,0] - lDY[25,0] = 0
 ldSC: ldx[25,0] - lDY[25,0] = 0
 ldSC: lDX[25,1] - lDY[25,1] = 0
 ldSC: ldx[25,1] - lDY[25,1] = 0
 ldSC: lDX[25,2] - lDY[25,2] = 0
 ldSC: ldx[25,2] - lDY[25,2] = 0
 ldSC: lDX[25,3] - lDY[25,3] = 0
 ldSC: ldx[25,3] - lDY[25,3] = 0
 ldSC: lDX[25,4] - lDY[25,4] = 0
 ldSC: ldx[25,4] - lDY[25,4] = 0
 ldSC: lDX[25,5] - lDY[25,5] = 0
 ldSC: ldx[25,5] - lDY[25,5] = 0
 ldSC: lDX[25,6] - lDY[25,6] = 0
 ldSC: ldx[25,6] - lDY[25,6] = 0
 ldSC: lDX[25,7] - lDY[25,7] = 0
 ldSC: ldx[25,7] - lDY[25,7] = 0
 ldSC: lDX[25,8] - lDY[25,8] = 0
 ldSC: ldx[25,8] - lDY[25,8] = 0
 ldSC: lDX[25,9] - lDY[25,9] = 0
 ldSC: ldx[25,9] - lDY[25,9] = 0
 ldSC: lDX[25,10] - lDY[25,10] = 0
 ldSC: ldx[25,10] - lDY[25,10] = 0
 ldSC: lDX[25,11] - lDY[25,11] = 0
 ldSC: ldx[25,11] - lDY[25,11] = 0
 ldSC: lDX[25,12] - lDY[25,12] = 0
 ldSC: ldx[25,12] - lDY[25,12] = 0
 ldSC: lDX[25,13] - lDY[25,13] = 0
 ldSC: ldx[25,13] - lDY[25,13] = 0
 ldSC: lDX[25,14] - lDY[25,14] = 0
 ldSC: ldx[25,14] - lDY[25,14] = 0
 ldSC: lDX[25,15] - lDY[25,15] = 0
 ldSC: ldx[25,15] - lDY[25,15] = 0
 ldSC: lDX[26,0] - lDY[26,0] = 0
 ldSC: ldx[26,0] - lDY[26,0] = 0
 ldSC: lDX[26,1] - lDY[26,1] = 0
 ldSC: ldx[26,1] - lDY[26,1] = 0
 ldSC: lDX[26,2] - lDY[26,2] = 0
 ldSC: ldx[26,2] - lDY[26,2] = 0
 ldSC: lDX[26,3] - lDY[26,3] = 0
 ldSC: ldx[26,3] - lDY[26,3] = 0
 ldSC: lDX[26,4] - lDY[26,4] = 0
 ldSC: ldx[26,4] - lDY[26,4] = 0
 ldSC: lDX[26,5] - lDY[26,5] = 0
 ldSC: ldx[26,5] - lDY[26,5] = 0
 ldSC: lDX[26,6] - lDY[26,6] = 0
 ldSC: ldx[26,6] - lDY[26,6] = 0
 ldSC: lDX[26,7] - lDY[26,7] = 0
 ldSC: ldx[26,7] - lDY[26,7] = 0
 ldSC: lDX[26,8] - lDY[26,8] = 0
 ldSC: ldx[26,8] - lDY[26,8] = 0
 ldSC: lDX[26,9] - lDY[26,9] = 0
 ldSC: ldx[26,9] - lDY[26,9] = 0
 ldSC: lDX[26,10] - lDY[26,10] = 0
 ldSC: ldx[26,10] - lDY[26,10] = 0
 ldSC: lDX[26,11] - lDY[26,11] = 0
 ldSC: ldx[26,11] - lDY[26,11] = 0
 ldSC: lDX[26,12] - lDY[26,12] = 0
 ldSC: ldx[26,12] - lDY[26,12] = 0
 ldSC: lDX[26,13] - lDY[26,13] = 0
 ldSC: ldx[26,13] - lDY[26,13] = 0
 ldSC: lDX[26,14] - lDY[26,14] = 0
 ldSC: ldx[26,14] - lDY[26,14] = 0
 ldSC: lDX[26,15] - lDY[26,15] = 0
 ldSC: ldx[26,15] - lDY[26,15] = 0
 ldSC: lDX[27,0] - lDY[27,0] = 0
 ldSC: ldx[27,0] - lDY[27,0] = 0
 ldSC: lDX[27,1] - lDY[27,1] = 0
 ldSC: ldx[27,1] - lDY[27,1] = 0
 ldSC: lDX[27,2] - lDY[27,2] = 0
 ldSC: ldx[27,2] - lDY[27,2] = 0
 ldSC: lDX[27,3] - lDY[27,3] = 0
 ldSC: ldx[27,3] - lDY[27,3] = 0
 ldSC: lDX[27,4] - lDY[27,4] = 0
 ldSC: ldx[27,4] - lDY[27,4] = 0
 ldSC: lDX[27,5] - lDY[27,5] = 0
 ldSC: ldx[27,5] - lDY[27,5] = 0
 ldSC: lDX[27,6] - lDY[27,6] = 0
 ldSC: ldx[27,6] - lDY[27,6] = 0
 ldSC: lDX[27,7] - lDY[27,7] = 0
 ldSC: ldx[27,7] - lDY[27,7] = 0
 ldSC: lDX[27,8] - lDY[27,8] = 0
 ldSC: ldx[27,8] - lDY[27,8] = 0
 ldSC: lDX[27,9] - lDY[27,9] = 0
 ldSC: ldx[27,9] - lDY[27,9] = 0
 ldSC: lDX[27,10] - lDY[27,10] = 0
 ldSC: ldx[27,10] - lDY[27,10] = 0
 ldSC: lDX[27,11] - lDY[27,11] = 0
 ldSC: ldx[27,11] - lDY[27,11] = 0
 ldSC: lDX[27,12] - lDY[27,12] = 0
 ldSC: ldx[27,12] - lDY[27,12] = 0
 ldSC: lDX[27,13] - lDY[27,13] = 0
 ldSC: ldx[27,13] - lDY[27,13] = 0
 ldSC: lDX[27,14] - lDY[27,14] = 0
 ldSC: ldx[27,14] - lDY[27,14] = 0
 ldSC: lDX[27,15] - lDY[27,15] = 0
 ldSC: ldx[27,15] - lDY[27,15] = 0
 leSC: - lDX[28,0] + lDY[28,0] = 0
 leSC: - lDX[28,0] + ldy[28,0] = 0
 leSC: - lDX[28,1] + lDY[28,1] = 0
 leSC: - lDX[28,1] + ldy[28,1] = 0
 leSC: - lDX[28,2] + lDY[28,2] = 0
 leSC: - lDX[28,2] + ldy[28,2] = 0
 leSC: - lDX[28,3] + lDY[28,3] = 0
 leSC: - lDX[28,3] + ldy[28,3] = 0
 leSC: - lDX[28,4] + lDY[28,4] = 0
 leSC: - lDX[28,4] + ldy[28,4] = 0
 leSC: - lDX[28,5] + lDY[28,5] = 0
 leSC: - lDX[28,5] + ldy[28,5] = 0
 leSC: - lDX[28,6] + lDY[28,6] = 0
 leSC: - lDX[28,6] + ldy[28,6] = 0
 leSC: - lDX[28,7] + lDY[28,7] = 0
 leSC: - lDX[28,7] + ldy[28,7] = 0
 leSC: - lDX[28,8] + lDY[28,8] = 0
 leSC: - lDX[28,8] + ldy[28,8] = 0
 leSC: - lDX[28,9] + lDY[28,9] = 0
 leSC: - lDX[28,9] + ldy[28,9] = 0
 leSC: - lDX[28,10] + lDY[28,10] = 0
 leSC: - lDX[28,10] + ldy[28,10] = 0
 leSC: - lDX[28,11] + lDY[28,11] = 0
 leSC: - lDX[28,11] + ldy[28,11] = 0
 leSC: - lDX[28,12] + lDY[28,12] = 0
 leSC: - lDX[28,12] + ldy[28,12] = 0
 leSC: - lDX[28,13] + lDY[28,13] = 0
 leSC: - lDX[28,13] + ldy[28,13] = 0
 leSC: - lDX[28,14] + lDY[28,14] = 0
 leSC: - lDX[28,14] + ldy[28,14] = 0
 leSC: - lDX[28,15] + lDY[28,15] = 0
 leSC: - lDX[28,15] + ldy[28,15] = 0
 leSC: - lDX[29,0] + lDY[29,0] = 0
 leSC: - lDX[29,0] + ldy[29,0] = 0
 leSC: - lDX[29,1] + lDY[29,1] = 0
 leSC: - lDX[29,1] + ldy[29,1] = 0
 leSC: - lDX[29,2] + lDY[29,2] = 0
 leSC: - lDX[29,2] + ldy[29,2] = 0
 leSC: - lDX[29,3] + lDY[29,3] = 0
 leSC: - lDX[29,3] + ldy[29,3] = 0
 leSC: - lDX[29,4] + lDY[29,4] = 0
 leSC: - lDX[29,4] + ldy[29,4] = 0
 leSC: - lDX[29,5] + lDY[29,5] = 0
 leSC: - lDX[29,5] + ldy[29,5] = 0
 leSC: - lDX[29,6] + lDY[29,6] = 0
 leSC: - lDX[29,6] + ldy[29,6] = 0
 leSC: - lDX[29,7] + lDY[29,7] = 0
 leSC: - lDX[29,7] + ldy[29,7] = 0
 leSC: - lDX[29,8] + lDY[29,8] = 0
 leSC: - lDX[29,8] + ldy[29,8] = 0
 leSC: - lDX[29,9] + lDY[29,9] = 0
 leSC: - lDX[29,9] + ldy[29,9] = 0
 leSC: - lDX[29,10] + lDY[29,10] = 0
 leSC: - lDX[29,10] + ldy[29,10] = 0
 leSC: - lDX[29,11] + lDY[29,11] = 0
 leSC: - lDX[29,11] + ldy[29,11] = 0
 leSC: - lDX[29,12] + lDY[29,12] = 0
 leSC: - lDX[29,12] + ldy[29,12] = 0
 leSC: - lDX[29,13] + lDY[29,13] = 0
 leSC: - lDX[29,13] + ldy[29,13] = 0
 leSC: - lDX[29,14] + lDY[29,14] = 0
 leSC: - lDX[29,14] + ldy[29,14] = 0
 leSC: - lDX[29,15] + lDY[29,15] = 0
 leSC: - lDX[29,15] + ldy[29,15] = 0
 leSC: - lDX[30,0] + lDY[30,0] = 0
 leSC: - lDX[30,0] + ldy[30,0] = 0
 leSC: - lDX[30,1] + lDY[30,1] = 0
 leSC: - lDX[30,1] + ldy[30,1] = 0
 leSC: - lDX[30,2] + lDY[30,2] = 0
 leSC: - lDX[30,2] + ldy[30,2] = 0
 leSC: - lDX[30,3] + lDY[30,3] = 0
 leSC: - lDX[30,3] + ldy[30,3] = 0
 leSC: - lDX[30,4] + lDY[30,4] = 0
 leSC: - lDX[30,4] + ldy[30,4] = 0
 leSC: - lDX[30,5] + lDY[30,5] = 0
 leSC: - lDX[30,5] + ldy[30,5] = 0
 leSC: - lDX[30,6] + lDY[30,6] = 0
 leSC: - lDX[30,6] + ldy[30,6] = 0
 leSC: - lDX[30,7] + lDY[30,7] = 0
 leSC: - lDX[30,7] + ldy[30,7] = 0
 leSC: - lDX[30,8] + lDY[30,8] = 0
 leSC: - lDX[30,8] + ldy[30,8] = 0
 leSC: - lDX[30,9] + lDY[30,9] = 0
 leSC: - lDX[30,9] + ldy[30,9] = 0
 leSC: - lDX[30,10] + lDY[30,10] = 0
 leSC: - lDX[30,10] + ldy[30,10] = 0
 leSC: - lDX[30,11] + lDY[30,11] = 0
 leSC: - lDX[30,11] + ldy[30,11] = 0
 leSC: - lDX[30,12] + lDY[30,12] = 0
 leSC: - lDX[30,12] + ldy[30,12] = 0
 leSC: - lDX[30,13] + lDY[30,13] = 0
 leSC: - lDX[30,13] + ldy[30,13] = 0
 leSC: - lDX[30,14] + lDY[30,14] = 0
 leSC: - lDX[30,14] + ldy[30,14] = 0
 leSC: - lDX[30,15] + lDY[30,15] = 0
 leSC: - lDX[30,15] + ldy[30,15] = 0
 leSC: - lDX[31,0] + lDY[31,0] = 0
 leSC: - lDX[31,0] + ldy[31,0] = 0
 leSC: - lDX[31,1] + lDY[31,1] = 0
 leSC: - lDX[31,1] + ldy[31,1] = 0
 leSC: - lDX[31,2] + lDY[31,2] = 0
 leSC: - lDX[31,2] + ldy[31,2] = 0
 leSC: - lDX[31,3] + lDY[31,3] = 0
 leSC: - lDX[31,3] + ldy[31,3] = 0
 leSC: - lDX[31,4] + lDY[31,4] = 0
 leSC: - lDX[31,4] + ldy[31,4] = 0
 leSC: - lDX[31,5] + lDY[31,5] = 0
 leSC: - lDX[31,5] + ldy[31,5] = 0
 leSC: - lDX[31,6] + lDY[31,6] = 0
 leSC: - lDX[31,6] + ldy[31,6] = 0
 leSC: - lDX[31,7] + lDY[31,7] = 0
 leSC: - lDX[31,7] + ldy[31,7] = 0
 leSC: - lDX[31,8] + lDY[31,8] = 0
 leSC: - lDX[31,8] + ldy[31,8] = 0
 leSC: - lDX[31,9] + lDY[31,9] = 0
 leSC: - lDX[31,9] + ldy[31,9] = 0
 leSC: - lDX[31,10] + lDY[31,10] = 0
 leSC: - lDX[31,10] + ldy[31,10] = 0
 leSC: - lDX[31,11] + lDY[31,11] = 0
 leSC: - lDX[31,11] + ldy[31,11] = 0
 leSC: - lDX[31,12] + lDY[31,12] = 0
 leSC: - lDX[31,12] + ldy[31,12] = 0
 leSC: - lDX[31,13] + lDY[31,13] = 0
 leSC: - lDX[31,13] + ldy[31,13] = 0
 leSC: - lDX[31,14] + lDY[31,14] = 0
 leSC: - lDX[31,14] + ldy[31,14] = 0
 leSC: - lDX[31,15] + lDY[31,15] = 0
 leSC: - lDX[31,15] + ldy[31,15] = 0
 leSC: - lDX[32,0] + lDY[32,0] = 0
 leSC: - lDX[32,0] + ldy[32,0] = 0
 leSC: - lDX[32,1] + lDY[32,1] = 0
 leSC: - lDX[32,1] + ldy[32,1] = 0
 leSC: - lDX[32,2] + lDY[32,2] = 0
 leSC: - lDX[32,2] + ldy[32,2] = 0
 leSC: - lDX[32,3] + lDY[32,3] = 0
 leSC: - lDX[32,3] + ldy[32,3] = 0
 leSC: - lDX[32,4] + lDY[32,4] = 0
 leSC: - lDX[32,4] + ldy[32,4] = 0
 leSC: - lDX[32,5] + lDY[32,5] = 0
 leSC: - lDX[32,5] + ldy[32,5] = 0
 leSC: - lDX[32,6] + lDY[32,6] = 0
 leSC: - lDX[32,6] + ldy[32,6] = 0
 leSC: - lDX[32,7] + lDY[32,7] = 0
 leSC: - lDX[32,7] + ldy[32,7] = 0
 leSC: - lDX[32,8] + lDY[32,8] = 0
 leSC: - lDX[32,8] + ldy[32,8] = 0
 leSC: - lDX[32,9] + lDY[32,9] = 0
 leSC: - lDX[32,9] + ldy[32,9] = 0
 leSC: - lDX[32,10] + lDY[32,10] = 0
 leSC: - lDX[32,10] + ldy[32,10] = 0
 leSC: - lDX[32,11] + lDY[32,11] = 0
 leSC: - lDX[32,11] + ldy[32,11] = 0
 leSC: - lDX[32,12] + lDY[32,12] = 0
 leSC: - lDX[32,12] + ldy[32,12] = 0
 leSC: - lDX[32,13] + lDY[32,13] = 0
 leSC: - lDX[32,13] + ldy[32,13] = 0
 leSC: - lDX[32,14] + lDY[32,14] = 0
 leSC: - lDX[32,14] + ldy[32,14] = 0
 leSC: - lDX[32,15] + lDY[32,15] = 0
 leSC: - lDX[32,15] + ldy[32,15] = 0
 lATK: ldy[17,0] - ldz[17,0] = 0
 lATK: lDY[17,0] - lDZ[17,0] + lstk[17,0] >= 0
 lATK: - lDY[17,0] + lDZ[17,0] + lstk[17,0] >= 0
 lATK: lDY[17,0] + lDZ[17,0] - ldz[17,0] - lstk[17,0] >= 0
 lATK: lDY[17,0] - ldz[17,0] >= 0
 lATK: ldy[17,1] - ldz[17,1] = 0
 lATK: lDY[17,1] - lDZ[17,1] + lstk[17,1] >= 0
 lATK: - lDY[17,1] + lDZ[17,1] + lstk[17,1] >= 0
 lATK: lDY[17,1] + lDZ[17,1] - ldz[17,1] - lstk[17,1] >= 0
 lATK: lDY[17,1] - ldz[17,1] >= 0
 lATK: ldy[17,2] - ldz[17,2] = 0
 lATK: lDY[17,2] - lDZ[17,2] + lstk[17,2] >= 0
 lATK: - lDY[17,2] + lDZ[17,2] + lstk[17,2] >= 0
 lATK: lDY[17,2] + lDZ[17,2] - ldz[17,2] - lstk[17,2] >= 0
 lATK: lDY[17,2] - ldz[17,2] >= 0
 lATK: ldy[17,3] - ldz[17,3] = 0
 lATK: lDY[17,3] - lDZ[17,3] + lstk[17,3] >= 0
 lATK: - lDY[17,3] + lDZ[17,3] + lstk[17,3] >= 0
 lATK: lDY[17,3] + lDZ[17,3] - ldz[17,3] - lstk[17,3] >= 0
 lATK: lDY[17,3] - ldz[17,3] >= 0
 lATK: ldy[17,4] - ldz[17,4] = 0
 lATK: lDY[17,4] - lDZ[17,4] + lstk[17,4] >= 0
 lATK: - lDY[17,4] + lDZ[17,4] + lstk[17,4] >= 0
 lATK: lDY[17,4] + lDZ[17,4] - ldz[17,4] - lstk[17,4] >= 0
 lATK: lDY[17,4] - ldz[17,4] >= 0
 lATK: ldy[17,5] - ldz[17,5] = 0
 lATK: lDY[17,5] - lDZ[17,5] + lstk[17,5] >= 0
 lATK: - lDY[17,5] + lDZ[17,5] + lstk[17,5] >= 0
 lATK: lDY[17,5] + lDZ[17,5] - ldz[17,5] - lstk[17,5] >= 0
 lATK: lDY[17,5] - ldz[17,5] >= 0
 lATK: ldy[17,6] - ldz[17,6] = 0
 lATK: lDY[17,6] - lDZ[17,6] + lstk[17,6] >= 0
 lATK: - lDY[17,6] + lDZ[17,6] + lstk[17,6] >= 0
 lATK: lDY[17,6] + lDZ[17,6] - ldz[17,6] - lstk[17,6] >= 0
 lATK: lDY[17,6] - ldz[17,6] >= 0
 lATK: ldy[17,7] - ldz[17,7] = 0
 lATK: lDY[17,7] - lDZ[17,7] + lstk[17,7] >= 0
 lATK: - lDY[17,7] + lDZ[17,7] + lstk[17,7] >= 0
 lATK: lDY[17,7] + lDZ[17,7] - ldz[17,7] - lstk[17,7] >= 0
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
 lATK: lDY[18,0] - lDZ[18,0] + lstk[18,0] >= 0
 lATK: - lDY[18,0] + lDZ[18,0] + lstk[18,0] >= 0
 lATK: lDY[18,0] + lDZ[18,0] - ldz[18,0] - lstk[18,0] >= 0
 lATK: lDY[18,0] - ldz[18,0] >= 0
 lATK: ldy[18,1] - ldz[18,1] = 0
 lATK: lDY[18,1] - lDZ[18,1] + lstk[18,1] >= 0
 lATK: - lDY[18,1] + lDZ[18,1] + lstk[18,1] >= 0
 lATK: lDY[18,1] + lDZ[18,1] - ldz[18,1] - lstk[18,1] >= 0
 lATK: lDY[18,1] - ldz[18,1] >= 0
 lATK: ldy[18,2] - ldz[18,2] = 0
 lATK: lDY[18,2] - lDZ[18,2] + lstk[18,2] >= 0
 lATK: - lDY[18,2] + lDZ[18,2] + lstk[18,2] >= 0
 lATK: lDY[18,2] + lDZ[18,2] - ldz[18,2] - lstk[18,2] >= 0
 lATK: lDY[18,2] - ldz[18,2] >= 0
 lATK: ldy[18,3] - ldz[18,3] = 0
 lATK: lDY[18,3] - lDZ[18,3] + lstk[18,3] >= 0
 lATK: - lDY[18,3] + lDZ[18,3] + lstk[18,3] >= 0
 lATK: lDY[18,3] + lDZ[18,3] - ldz[18,3] - lstk[18,3] >= 0
 lATK: lDY[18,3] - ldz[18,3] >= 0
 lATK: ldy[18,4] - ldz[18,4] = 0
 lATK: lDY[18,4] - lDZ[18,4] + lstk[18,4] >= 0
 lATK: - lDY[18,4] + lDZ[18,4] + lstk[18,4] >= 0
 lATK: lDY[18,4] + lDZ[18,4] - ldz[18,4] - lstk[18,4] >= 0
 lATK: lDY[18,4] - ldz[18,4] >= 0
 lATK: ldy[18,5] - ldz[18,5] = 0
 lATK: lDY[18,5] - lDZ[18,5] + lstk[18,5] >= 0
 lATK: - lDY[18,5] + lDZ[18,5] + lstk[18,5] >= 0
 lATK: lDY[18,5] + lDZ[18,5] - ldz[18,5] - lstk[18,5] >= 0
 lATK: lDY[18,5] - ldz[18,5] >= 0
 lATK: ldy[18,6] - ldz[18,6] = 0
 lATK: lDY[18,6] - lDZ[18,6] + lstk[18,6] >= 0
 lATK: - lDY[18,6] + lDZ[18,6] + lstk[18,6] >= 0
 lATK: lDY[18,6] + lDZ[18,6] - ldz[18,6] - lstk[18,6] >= 0
 lATK: lDY[18,6] - ldz[18,6] >= 0
 lATK: ldy[18,7] - ldz[18,7] = 0
 lATK: lDY[18,7] - lDZ[18,7] + lstk[18,7] >= 0
 lATK: - lDY[18,7] + lDZ[18,7] + lstk[18,7] >= 0
 lATK: lDY[18,7] + lDZ[18,7] - ldz[18,7] - lstk[18,7] >= 0
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
 lATK: lDY[19,0] - lDZ[19,0] + lstk[19,0] >= 0
 lATK: - lDY[19,0] + lDZ[19,0] + lstk[19,0] >= 0
 lATK: lDY[19,0] + lDZ[19,0] - ldz[19,0] - lstk[19,0] >= 0
 lATK: lDY[19,0] - ldz[19,0] >= 0
 lATK: ldy[19,1] - ldz[19,1] = 0
 lATK: lDY[19,1] - lDZ[19,1] + lstk[19,1] >= 0
 lATK: - lDY[19,1] + lDZ[19,1] + lstk[19,1] >= 0
 lATK: lDY[19,1] + lDZ[19,1] - ldz[19,1] - lstk[19,1] >= 0
 lATK: lDY[19,1] - ldz[19,1] >= 0
 lATK: ldy[19,2] - ldz[19,2] = 0
 lATK: lDY[19,2] - lDZ[19,2] + lstk[19,2] >= 0
 lATK: - lDY[19,2] + lDZ[19,2] + lstk[19,2] >= 0
 lATK: lDY[19,2] + lDZ[19,2] - ldz[19,2] - lstk[19,2] >= 0
 lATK: lDY[19,2] - ldz[19,2] >= 0
 lATK: ldy[19,3] - ldz[19,3] = 0
 lATK: lDY[19,3] - lDZ[19,3] + lstk[19,3] >= 0
 lATK: - lDY[19,3] + lDZ[19,3] + lstk[19,3] >= 0
 lATK: lDY[19,3] + lDZ[19,3] - ldz[19,3] - lstk[19,3] >= 0
 lATK: lDY[19,3] - ldz[19,3] >= 0
 lATK: ldy[19,4] - ldz[19,4] = 0
 lATK: lDY[19,4] - lDZ[19,4] + lstk[19,4] >= 0
 lATK: - lDY[19,4] + lDZ[19,4] + lstk[19,4] >= 0
 lATK: lDY[19,4] + lDZ[19,4] - ldz[19,4] - lstk[19,4] >= 0
 lATK: lDY[19,4] - ldz[19,4] >= 0
 lATK: ldy[19,5] - ldz[19,5] = 0
 lATK: lDY[19,5] - lDZ[19,5] + lstk[19,5] >= 0
 lATK: - lDY[19,5] + lDZ[19,5] + lstk[19,5] >= 0
 lATK: lDY[19,5] + lDZ[19,5] - ldz[19,5] - lstk[19,5] >= 0
 lATK: lDY[19,5] - ldz[19,5] >= 0
 lATK: ldy[19,6] - ldz[19,6] = 0
 lATK: lDY[19,6] - lDZ[19,6] + lstk[19,6] >= 0
 lATK: - lDY[19,6] + lDZ[19,6] + lstk[19,6] >= 0
 lATK: lDY[19,6] + lDZ[19,6] - ldz[19,6] - lstk[19,6] >= 0
 lATK: lDY[19,6] - ldz[19,6] >= 0
 lATK: ldy[19,7] - ldz[19,7] = 0
 lATK: lDY[19,7] - lDZ[19,7] + lstk[19,7] >= 0
 lATK: - lDY[19,7] + lDZ[19,7] + lstk[19,7] >= 0
 lATK: lDY[19,7] + lDZ[19,7] - ldz[19,7] - lstk[19,7] >= 0
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
 lATK: lDY[20,0] - lDZ[20,0] + lstk[20,0] >= 0
 lATK: - lDY[20,0] + lDZ[20,0] + lstk[20,0] >= 0
 lATK: lDY[20,0] + lDZ[20,0] - ldz[20,0] - lstk[20,0] >= 0
 lATK: lDY[20,0] - ldz[20,0] >= 0
 lATK: ldy[20,1] - ldz[20,1] = 0
 lATK: lDY[20,1] - lDZ[20,1] + lstk[20,1] >= 0
 lATK: - lDY[20,1] + lDZ[20,1] + lstk[20,1] >= 0
 lATK: lDY[20,1] + lDZ[20,1] - ldz[20,1] - lstk[20,1] >= 0
 lATK: lDY[20,1] - ldz[20,1] >= 0
 lATK: ldy[20,2] - ldz[20,2] = 0
 lATK: lDY[20,2] - lDZ[20,2] + lstk[20,2] >= 0
 lATK: - lDY[20,2] + lDZ[20,2] + lstk[20,2] >= 0
 lATK: lDY[20,2] + lDZ[20,2] - ldz[20,2] - lstk[20,2] >= 0
 lATK: lDY[20,2] - ldz[20,2] >= 0
 lATK: ldy[20,3] - ldz[20,3] = 0
 lATK: lDY[20,3] - lDZ[20,3] + lstk[20,3] >= 0
 lATK: - lDY[20,3] + lDZ[20,3] + lstk[20,3] >= 0
 lATK: lDY[20,3] + lDZ[20,3] - ldz[20,3] - lstk[20,3] >= 0
 lATK: lDY[20,3] - ldz[20,3] >= 0
 lATK: ldy[20,4] - ldz[20,4] = 0
 lATK: lDY[20,4] - lDZ[20,4] + lstk[20,4] >= 0
 lATK: - lDY[20,4] + lDZ[20,4] + lstk[20,4] >= 0
 lATK: lDY[20,4] + lDZ[20,4] - ldz[20,4] - lstk[20,4] >= 0
 lATK: lDY[20,4] - ldz[20,4] >= 0
 lATK: ldy[20,5] - ldz[20,5] = 0
 lATK: lDY[20,5] - lDZ[20,5] + lstk[20,5] >= 0
 lATK: - lDY[20,5] + lDZ[20,5] + lstk[20,5] >= 0
 lATK: lDY[20,5] + lDZ[20,5] - ldz[20,5] - lstk[20,5] >= 0
 lATK: lDY[20,5] - ldz[20,5] >= 0
 lATK: ldy[20,6] - ldz[20,6] = 0
 lATK: lDY[20,6] - lDZ[20,6] + lstk[20,6] >= 0
 lATK: - lDY[20,6] + lDZ[20,6] + lstk[20,6] >= 0
 lATK: lDY[20,6] + lDZ[20,6] - ldz[20,6] - lstk[20,6] >= 0
 lATK: lDY[20,6] - ldz[20,6] >= 0
 lATK: ldy[20,7] - ldz[20,7] = 0
 lATK: lDY[20,7] - lDZ[20,7] + lstk[20,7] >= 0
 lATK: - lDY[20,7] + lDZ[20,7] + lstk[20,7] >= 0
 lATK: lDY[20,7] + lDZ[20,7] - ldz[20,7] - lstk[20,7] >= 0
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
 lATK: lDY[21,0] - lDZ[21,0] + lstk[21,0] >= 0
 lATK: - lDY[21,0] + lDZ[21,0] + lstk[21,0] >= 0
 lATK: lDY[21,0] + lDZ[21,0] - ldz[21,0] - lstk[21,0] >= 0
 lATK: lDY[21,0] - ldz[21,0] >= 0
 lATK: ldy[21,1] - ldz[21,1] = 0
 lATK: lDY[21,1] - lDZ[21,1] + lstk[21,1] >= 0
 lATK: - lDY[21,1] + lDZ[21,1] + lstk[21,1] >= 0
 lATK: lDY[21,1] + lDZ[21,1] - ldz[21,1] - lstk[21,1] >= 0
 lATK: lDY[21,1] - ldz[21,1] >= 0
 lATK: ldy[21,2] - ldz[21,2] = 0
 lATK: lDY[21,2] - lDZ[21,2] + lstk[21,2] >= 0
 lATK: - lDY[21,2] + lDZ[21,2] + lstk[21,2] >= 0
 lATK: lDY[21,2] + lDZ[21,2] - ldz[21,2] - lstk[21,2] >= 0
 lATK: lDY[21,2] - ldz[21,2] >= 0
 lATK: ldy[21,3] - ldz[21,3] = 0
 lATK: lDY[21,3] - lDZ[21,3] + lstk[21,3] >= 0
 lATK: - lDY[21,3] + lDZ[21,3] + lstk[21,3] >= 0
 lATK: lDY[21,3] + lDZ[21,3] - ldz[21,3] - lstk[21,3] >= 0
 lATK: lDY[21,3] - ldz[21,3] >= 0
 lATK: ldy[21,4] - ldz[21,4] = 0
 lATK: lDY[21,4] - lDZ[21,4] + lstk[21,4] >= 0
 lATK: - lDY[21,4] + lDZ[21,4] + lstk[21,4] >= 0
 lATK: lDY[21,4] + lDZ[21,4] - ldz[21,4] - lstk[21,4] >= 0
 lATK: lDY[21,4] - ldz[21,4] >= 0
 lATK: ldy[21,5] - ldz[21,5] = 0
 lATK: lDY[21,5] - lDZ[21,5] + lstk[21,5] >= 0
 lATK: - lDY[21,5] + lDZ[21,5] + lstk[21,5] >= 0
 lATK: lDY[21,5] + lDZ[21,5] - ldz[21,5] - lstk[21,5] >= 0
 lATK: lDY[21,5] - ldz[21,5] >= 0
 lATK: ldy[21,6] - ldz[21,6] = 0
 lATK: lDY[21,6] - lDZ[21,6] + lstk[21,6] >= 0
 lATK: - lDY[21,6] + lDZ[21,6] + lstk[21,6] >= 0
 lATK: lDY[21,6] + lDZ[21,6] - ldz[21,6] - lstk[21,6] >= 0
 lATK: lDY[21,6] - ldz[21,6] >= 0
 lATK: ldy[21,7] - ldz[21,7] = 0
 lATK: lDY[21,7] - lDZ[21,7] + lstk[21,7] >= 0
 lATK: - lDY[21,7] + lDZ[21,7] + lstk[21,7] >= 0
 lATK: lDY[21,7] + lDZ[21,7] - ldz[21,7] - lstk[21,7] >= 0
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
 lATK: lDY[22,0] - lDZ[22,0] + lstk[22,0] >= 0
 lATK: - lDY[22,0] + lDZ[22,0] + lstk[22,0] >= 0
 lATK: lDY[22,0] + lDZ[22,0] - ldz[22,0] - lstk[22,0] >= 0
 lATK: lDY[22,0] - ldz[22,0] >= 0
 lATK: ldy[22,1] - ldz[22,1] = 0
 lATK: lDY[22,1] - lDZ[22,1] + lstk[22,1] >= 0
 lATK: - lDY[22,1] + lDZ[22,1] + lstk[22,1] >= 0
 lATK: lDY[22,1] + lDZ[22,1] - ldz[22,1] - lstk[22,1] >= 0
 lATK: lDY[22,1] - ldz[22,1] >= 0
 lATK: ldy[22,2] - ldz[22,2] = 0
 lATK: lDY[22,2] - lDZ[22,2] + lstk[22,2] >= 0
 lATK: - lDY[22,2] + lDZ[22,2] + lstk[22,2] >= 0
 lATK: lDY[22,2] + lDZ[22,2] - ldz[22,2] - lstk[22,2] >= 0
 lATK: lDY[22,2] - ldz[22,2] >= 0
 lATK: ldy[22,3] - ldz[22,3] = 0
 lATK: lDY[22,3] - lDZ[22,3] + lstk[22,3] >= 0
 lATK: - lDY[22,3] + lDZ[22,3] + lstk[22,3] >= 0
 lATK: lDY[22,3] + lDZ[22,3] - ldz[22,3] - lstk[22,3] >= 0
 lATK: lDY[22,3] - ldz[22,3] >= 0
 lATK: ldy[22,4] - ldz[22,4] = 0
 lATK: lDY[22,4] - lDZ[22,4] + lstk[22,4] >= 0
 lATK: - lDY[22,4] + lDZ[22,4] + lstk[22,4] >= 0
 lATK: lDY[22,4] + lDZ[22,4] - ldz[22,4] - lstk[22,4] >= 0
 lATK: lDY[22,4] - ldz[22,4] >= 0
 lATK: ldy[22,5] - ldz[22,5] = 0
 lATK: lDY[22,5] - lDZ[22,5] + lstk[22,5] >= 0
 lATK: - lDY[22,5] + lDZ[22,5] + lstk[22,5] >= 0
 lATK: lDY[22,5] + lDZ[22,5] - ldz[22,5] - lstk[22,5] >= 0
 lATK: lDY[22,5] - ldz[22,5] >= 0
 lATK: ldy[22,6] - ldz[22,6] = 0
 lATK: lDY[22,6] - lDZ[22,6] + lstk[22,6] >= 0
 lATK: - lDY[22,6] + lDZ[22,6] + lstk[22,6] >= 0
 lATK: lDY[22,6] + lDZ[22,6] - ldz[22,6] - lstk[22,6] >= 0
 lATK: lDY[22,6] - ldz[22,6] >= 0
 lATK: ldy[22,7] - ldz[22,7] = 0
 lATK: lDY[22,7] - lDZ[22,7] + lstk[22,7] >= 0
 lATK: - lDY[22,7] + lDZ[22,7] + lstk[22,7] >= 0
 lATK: lDY[22,7] + lDZ[22,7] - ldz[22,7] - lstk[22,7] >= 0
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
 lATK: lDY[23,0] - lDZ[23,0] + lstk[23,0] >= 0
 lATK: - lDY[23,0] + lDZ[23,0] + lstk[23,0] >= 0
 lATK: lDY[23,0] + lDZ[23,0] - ldz[23,0] - lstk[23,0] >= 0
 lATK: lDY[23,0] - ldz[23,0] >= 0
 lATK: ldy[23,1] - ldz[23,1] = 0
 lATK: lDY[23,1] - lDZ[23,1] + lstk[23,1] >= 0
 lATK: - lDY[23,1] + lDZ[23,1] + lstk[23,1] >= 0
 lATK: lDY[23,1] + lDZ[23,1] - ldz[23,1] - lstk[23,1] >= 0
 lATK: lDY[23,1] - ldz[23,1] >= 0
 lATK: ldy[23,2] - ldz[23,2] = 0
 lATK: lDY[23,2] - lDZ[23,2] + lstk[23,2] >= 0
 lATK: - lDY[23,2] + lDZ[23,2] + lstk[23,2] >= 0
 lATK: lDY[23,2] + lDZ[23,2] - ldz[23,2] - lstk[23,2] >= 0
 lATK: lDY[23,2] - ldz[23,2] >= 0
 lATK: ldy[23,3] - ldz[23,3] = 0
 lATK: lDY[23,3] - lDZ[23,3] + lstk[23,3] >= 0
 lATK: - lDY[23,3] + lDZ[23,3] + lstk[23,3] >= 0
 lATK: lDY[23,3] + lDZ[23,3] - ldz[23,3] - lstk[23,3] >= 0
 lATK: lDY[23,3] - ldz[23,3] >= 0
 lATK: ldy[23,4] - ldz[23,4] = 0
 lATK: lDY[23,4] - lDZ[23,4] + lstk[23,4] >= 0
 lATK: - lDY[23,4] + lDZ[23,4] + lstk[23,4] >= 0
 lATK: lDY[23,4] + lDZ[23,4] - ldz[23,4] - lstk[23,4] >= 0
 lATK: lDY[23,4] - ldz[23,4] >= 0
 lATK: ldy[23,5] - ldz[23,5] = 0
 lATK: lDY[23,5] - lDZ[23,5] + lstk[23,5] >= 0
 lATK: - lDY[23,5] + lDZ[23,5] + lstk[23,5] >= 0
 lATK: lDY[23,5] + lDZ[23,5] - ldz[23,5] - lstk[23,5] >= 0
 lATK: lDY[23,5] - ldz[23,5] >= 0
 lATK: ldy[23,6] - ldz[23,6] = 0
 lATK: lDY[23,6] - lDZ[23,6] + lstk[23,6] >= 0
 lATK: - lDY[23,6] + lDZ[23,6] + lstk[23,6] >= 0
 lATK: lDY[23,6] + lDZ[23,6] - ldz[23,6] - lstk[23,6] >= 0
 lATK: lDY[23,6] - ldz[23,6] >= 0
 lATK: ldy[23,7] - ldz[23,7] = 0
 lATK: lDY[23,7] - lDZ[23,7] + lstk[23,7] >= 0
 lATK: - lDY[23,7] + lDZ[23,7] + lstk[23,7] >= 0
 lATK: lDY[23,7] + lDZ[23,7] - ldz[23,7] - lstk[23,7] >= 0
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
 lATK: ldy[24,0] - ldz[24,0] = 0
 lATK: lDY[24,0] - lDZ[24,0] + lstk[24,0] >= 0
 lATK: - lDY[24,0] + lDZ[24,0] + lstk[24,0] >= 0
 lATK: lDY[24,0] + lDZ[24,0] - ldz[24,0] - lstk[24,0] >= 0
 lATK: lDY[24,0] - ldz[24,0] >= 0
 lATK: ldy[24,1] - ldz[24,1] = 0
 lATK: lDY[24,1] - lDZ[24,1] + lstk[24,1] >= 0
 lATK: - lDY[24,1] + lDZ[24,1] + lstk[24,1] >= 0
 lATK: lDY[24,1] + lDZ[24,1] - ldz[24,1] - lstk[24,1] >= 0
 lATK: lDY[24,1] - ldz[24,1] >= 0
 lATK: ldy[24,2] - ldz[24,2] = 0
 lATK: lDY[24,2] - lDZ[24,2] + lstk[24,2] >= 0
 lATK: - lDY[24,2] + lDZ[24,2] + lstk[24,2] >= 0
 lATK: lDY[24,2] + lDZ[24,2] - ldz[24,2] - lstk[24,2] >= 0
 lATK: lDY[24,2] - ldz[24,2] >= 0
 lATK: ldy[24,3] - ldz[24,3] = 0
 lATK: lDY[24,3] - lDZ[24,3] + lstk[24,3] >= 0
 lATK: - lDY[24,3] + lDZ[24,3] + lstk[24,3] >= 0
 lATK: lDY[24,3] + lDZ[24,3] - ldz[24,3] - lstk[24,3] >= 0
 lATK: lDY[24,3] - ldz[24,3] >= 0
 lATK: ldy[24,4] - ldz[24,4] = 0
 lATK: lDY[24,4] - lDZ[24,4] + lstk[24,4] >= 0
 lATK: - lDY[24,4] + lDZ[24,4] + lstk[24,4] >= 0
 lATK: lDY[24,4] + lDZ[24,4] - ldz[24,4] - lstk[24,4] >= 0
 lATK: lDY[24,4] - ldz[24,4] >= 0
 lATK: ldy[24,5] - ldz[24,5] = 0
 lATK: lDY[24,5] - lDZ[24,5] + lstk[24,5] >= 0
 lATK: - lDY[24,5] + lDZ[24,5] + lstk[24,5] >= 0
 lATK: lDY[24,5] + lDZ[24,5] - ldz[24,5] - lstk[24,5] >= 0
 lATK: lDY[24,5] - ldz[24,5] >= 0
 lATK: ldy[24,6] - ldz[24,6] = 0
 lATK: lDY[24,6] - lDZ[24,6] + lstk[24,6] >= 0
 lATK: - lDY[24,6] + lDZ[24,6] + lstk[24,6] >= 0
 lATK: lDY[24,6] + lDZ[24,6] - ldz[24,6] - lstk[24,6] >= 0
 lATK: lDY[24,6] - ldz[24,6] >= 0
 lATK: ldy[24,7] - ldz[24,7] = 0
 lATK: lDY[24,7] - lDZ[24,7] + lstk[24,7] >= 0
 lATK: - lDY[24,7] + lDZ[24,7] + lstk[24,7] >= 0
 lATK: lDY[24,7] + lDZ[24,7] - ldz[24,7] - lstk[24,7] >= 0
 lATK: lDY[24,7] - ldz[24,7] >= 0
 lATK: lDY[24,8] - lDZ[24,8] = 0
 lATK: ldy[24,8] - ldz[24,8] = 0
 lATK: lDY[24,9] - lDZ[24,9] = 0
 lATK: ldy[24,9] - ldz[24,9] = 0
 lATK: lDY[24,10] - lDZ[24,10] = 0
 lATK: ldy[24,10] - ldz[24,10] = 0
 lATK: lDY[24,11] - lDZ[24,11] = 0
 lATK: ldy[24,11] - ldz[24,11] = 0
 lATK: lDY[24,12] - lDZ[24,12] = 0
 lATK: ldy[24,12] - ldz[24,12] = 0
 lATK: lDY[24,13] - lDZ[24,13] = 0
 lATK: ldy[24,13] - ldz[24,13] = 0
 lATK: lDY[24,14] - lDZ[24,14] = 0
 lATK: ldy[24,14] - ldz[24,14] = 0
 lATK: lDY[24,15] - lDZ[24,15] = 0
 lATK: ldy[24,15] - ldz[24,15] = 0
 lATK: ldy[25,0] - ldz[25,0] = 0
 lATK: lDY[25,0] - lDZ[25,0] + lstk[25,0] >= 0
 lATK: - lDY[25,0] + lDZ[25,0] + lstk[25,0] >= 0
 lATK: lDY[25,0] + lDZ[25,0] - ldz[25,0] - lstk[25,0] >= 0
 lATK: lDY[25,0] - ldz[25,0] >= 0
 lATK: ldy[25,1] - ldz[25,1] = 0
 lATK: lDY[25,1] - lDZ[25,1] + lstk[25,1] >= 0
 lATK: - lDY[25,1] + lDZ[25,1] + lstk[25,1] >= 0
 lATK: lDY[25,1] + lDZ[25,1] - ldz[25,1] - lstk[25,1] >= 0
 lATK: lDY[25,1] - ldz[25,1] >= 0
 lATK: ldy[25,2] - ldz[25,2] = 0
 lATK: lDY[25,2] - lDZ[25,2] + lstk[25,2] >= 0
 lATK: - lDY[25,2] + lDZ[25,2] + lstk[25,2] >= 0
 lATK: lDY[25,2] + lDZ[25,2] - ldz[25,2] - lstk[25,2] >= 0
 lATK: lDY[25,2] - ldz[25,2] >= 0
 lATK: ldy[25,3] - ldz[25,3] = 0
 lATK: lDY[25,3] - lDZ[25,3] + lstk[25,3] >= 0
 lATK: - lDY[25,3] + lDZ[25,3] + lstk[25,3] >= 0
 lATK: lDY[25,3] + lDZ[25,3] - ldz[25,3] - lstk[25,3] >= 0
 lATK: lDY[25,3] - ldz[25,3] >= 0
 lATK: ldy[25,4] - ldz[25,4] = 0
 lATK: lDY[25,4] - lDZ[25,4] + lstk[25,4] >= 0
 lATK: - lDY[25,4] + lDZ[25,4] + lstk[25,4] >= 0
 lATK: lDY[25,4] + lDZ[25,4] - ldz[25,4] - lstk[25,4] >= 0
 lATK: lDY[25,4] - ldz[25,4] >= 0
 lATK: ldy[25,5] - ldz[25,5] = 0
 lATK: lDY[25,5] - lDZ[25,5] + lstk[25,5] >= 0
 lATK: - lDY[25,5] + lDZ[25,5] + lstk[25,5] >= 0
 lATK: lDY[25,5] + lDZ[25,5] - ldz[25,5] - lstk[25,5] >= 0
 lATK: lDY[25,5] - ldz[25,5] >= 0
 lATK: ldy[25,6] - ldz[25,6] = 0
 lATK: lDY[25,6] - lDZ[25,6] + lstk[25,6] >= 0
 lATK: - lDY[25,6] + lDZ[25,6] + lstk[25,6] >= 0
 lATK: lDY[25,6] + lDZ[25,6] - ldz[25,6] - lstk[25,6] >= 0
 lATK: lDY[25,6] - ldz[25,6] >= 0
 lATK: ldy[25,7] - ldz[25,7] = 0
 lATK: lDY[25,7] - lDZ[25,7] + lstk[25,7] >= 0
 lATK: - lDY[25,7] + lDZ[25,7] + lstk[25,7] >= 0
 lATK: lDY[25,7] + lDZ[25,7] - ldz[25,7] - lstk[25,7] >= 0
 lATK: lDY[25,7] - ldz[25,7] >= 0
 lATK: lDY[25,8] - lDZ[25,8] = 0
 lATK: ldy[25,8] - ldz[25,8] = 0
 lATK: lDY[25,9] - lDZ[25,9] = 0
 lATK: ldy[25,9] - ldz[25,9] = 0
 lATK: lDY[25,10] - lDZ[25,10] = 0
 lATK: ldy[25,10] - ldz[25,10] = 0
 lATK: lDY[25,11] - lDZ[25,11] = 0
 lATK: ldy[25,11] - ldz[25,11] = 0
 lATK: lDY[25,12] - lDZ[25,12] = 0
 lATK: ldy[25,12] - ldz[25,12] = 0
 lATK: lDY[25,13] - lDZ[25,13] = 0
 lATK: ldy[25,13] - ldz[25,13] = 0
 lATK: lDY[25,14] - lDZ[25,14] = 0
 lATK: ldy[25,14] - ldz[25,14] = 0
 lATK: lDY[25,15] - lDZ[25,15] = 0
 lATK: ldy[25,15] - ldz[25,15] = 0
 lATK: ldy[26,0] - ldz[26,0] = 0
 lATK: lDY[26,0] - lDZ[26,0] + lstk[26,0] >= 0
 lATK: - lDY[26,0] + lDZ[26,0] + lstk[26,0] >= 0
 lATK: lDY[26,0] + lDZ[26,0] - ldz[26,0] - lstk[26,0] >= 0
 lATK: lDY[26,0] - ldz[26,0] >= 0
 lATK: ldy[26,1] - ldz[26,1] = 0
 lATK: lDY[26,1] - lDZ[26,1] + lstk[26,1] >= 0
 lATK: - lDY[26,1] + lDZ[26,1] + lstk[26,1] >= 0
 lATK: lDY[26,1] + lDZ[26,1] - ldz[26,1] - lstk[26,1] >= 0
 lATK: lDY[26,1] - ldz[26,1] >= 0
 lATK: ldy[26,2] - ldz[26,2] = 0
 lATK: lDY[26,2] - lDZ[26,2] + lstk[26,2] >= 0
 lATK: - lDY[26,2] + lDZ[26,2] + lstk[26,2] >= 0
 lATK: lDY[26,2] + lDZ[26,2] - ldz[26,2] - lstk[26,2] >= 0
 lATK: lDY[26,2] - ldz[26,2] >= 0
 lATK: ldy[26,3] - ldz[26,3] = 0
 lATK: lDY[26,3] - lDZ[26,3] + lstk[26,3] >= 0
 lATK: - lDY[26,3] + lDZ[26,3] + lstk[26,3] >= 0
 lATK: lDY[26,3] + lDZ[26,3] - ldz[26,3] - lstk[26,3] >= 0
 lATK: lDY[26,3] - ldz[26,3] >= 0
 lATK: ldy[26,4] - ldz[26,4] = 0
 lATK: lDY[26,4] - lDZ[26,4] + lstk[26,4] >= 0
 lATK: - lDY[26,4] + lDZ[26,4] + lstk[26,4] >= 0
 lATK: lDY[26,4] + lDZ[26,4] - ldz[26,4] - lstk[26,4] >= 0
 lATK: lDY[26,4] - ldz[26,4] >= 0
 lATK: ldy[26,5] - ldz[26,5] = 0
 lATK: lDY[26,5] - lDZ[26,5] + lstk[26,5] >= 0
 lATK: - lDY[26,5] + lDZ[26,5] + lstk[26,5] >= 0
 lATK: lDY[26,5] + lDZ[26,5] - ldz[26,5] - lstk[26,5] >= 0
 lATK: lDY[26,5] - ldz[26,5] >= 0
 lATK: ldy[26,6] - ldz[26,6] = 0
 lATK: lDY[26,6] - lDZ[26,6] + lstk[26,6] >= 0
 lATK: - lDY[26,6] + lDZ[26,6] + lstk[26,6] >= 0
 lATK: lDY[26,6] + lDZ[26,6] - ldz[26,6] - lstk[26,6] >= 0
 lATK: lDY[26,6] - ldz[26,6] >= 0
 lATK: ldy[26,7] - ldz[26,7] = 0
 lATK: lDY[26,7] - lDZ[26,7] + lstk[26,7] >= 0
 lATK: - lDY[26,7] + lDZ[26,7] + lstk[26,7] >= 0
 lATK: lDY[26,7] + lDZ[26,7] - ldz[26,7] - lstk[26,7] >= 0
 lATK: lDY[26,7] - ldz[26,7] >= 0
 lATK: lDY[26,8] - lDZ[26,8] = 0
 lATK: ldy[26,8] - ldz[26,8] = 0
 lATK: lDY[26,9] - lDZ[26,9] = 0
 lATK: ldy[26,9] - ldz[26,9] = 0
 lATK: lDY[26,10] - lDZ[26,10] = 0
 lATK: ldy[26,10] - ldz[26,10] = 0
 lATK: lDY[26,11] - lDZ[26,11] = 0
 lATK: ldy[26,11] - ldz[26,11] = 0
 lATK: lDY[26,12] - lDZ[26,12] = 0
 lATK: ldy[26,12] - ldz[26,12] = 0
 lATK: lDY[26,13] - lDZ[26,13] = 0
 lATK: ldy[26,13] - ldz[26,13] = 0
 lATK: lDY[26,14] - lDZ[26,14] = 0
 lATK: ldy[26,14] - ldz[26,14] = 0
 lATK: lDY[26,15] - lDZ[26,15] = 0
 lATK: ldy[26,15] - ldz[26,15] = 0
 lATK: ldy[27,0] - ldz[27,0] = 0
 lATK: lDY[27,0] - lDZ[27,0] + lstk[27,0] >= 0
 lATK: - lDY[27,0] + lDZ[27,0] + lstk[27,0] >= 0
 lATK: lDY[27,0] + lDZ[27,0] - ldz[27,0] - lstk[27,0] >= 0
 lATK: lDY[27,0] - ldz[27,0] >= 0
 lATK: ldy[27,1] - ldz[27,1] = 0
 lATK: lDY[27,1] - lDZ[27,1] + lstk[27,1] >= 0
 lATK: - lDY[27,1] + lDZ[27,1] + lstk[27,1] >= 0
 lATK: lDY[27,1] + lDZ[27,1] - ldz[27,1] - lstk[27,1] >= 0
 lATK: lDY[27,1] - ldz[27,1] >= 0
 lATK: ldy[27,2] - ldz[27,2] = 0
 lATK: lDY[27,2] - lDZ[27,2] + lstk[27,2] >= 0
 lATK: - lDY[27,2] + lDZ[27,2] + lstk[27,2] >= 0
 lATK: lDY[27,2] + lDZ[27,2] - ldz[27,2] - lstk[27,2] >= 0
 lATK: lDY[27,2] - ldz[27,2] >= 0
 lATK: ldy[27,3] - ldz[27,3] = 0
 lATK: lDY[27,3] - lDZ[27,3] + lstk[27,3] >= 0
 lATK: - lDY[27,3] + lDZ[27,3] + lstk[27,3] >= 0
 lATK: lDY[27,3] + lDZ[27,3] - ldz[27,3] - lstk[27,3] >= 0
 lATK: lDY[27,3] - ldz[27,3] >= 0
 lATK: ldy[27,4] - ldz[27,4] = 0
 lATK: lDY[27,4] - lDZ[27,4] + lstk[27,4] >= 0
 lATK: - lDY[27,4] + lDZ[27,4] + lstk[27,4] >= 0
 lATK: lDY[27,4] + lDZ[27,4] - ldz[27,4] - lstk[27,4] >= 0
 lATK: lDY[27,4] - ldz[27,4] >= 0
 lATK: ldy[27,5] - ldz[27,5] = 0
 lATK: lDY[27,5] - lDZ[27,5] + lstk[27,5] >= 0
 lATK: - lDY[27,5] + lDZ[27,5] + lstk[27,5] >= 0
 lATK: lDY[27,5] + lDZ[27,5] - ldz[27,5] - lstk[27,5] >= 0
 lATK: lDY[27,5] - ldz[27,5] >= 0
 lATK: ldy[27,6] - ldz[27,6] = 0
 lATK: lDY[27,6] - lDZ[27,6] + lstk[27,6] >= 0
 lATK: - lDY[27,6] + lDZ[27,6] + lstk[27,6] >= 0
 lATK: lDY[27,6] + lDZ[27,6] - ldz[27,6] - lstk[27,6] >= 0
 lATK: lDY[27,6] - ldz[27,6] >= 0
 lATK: ldy[27,7] - ldz[27,7] = 0
 lATK: lDY[27,7] - lDZ[27,7] + lstk[27,7] >= 0
 lATK: - lDY[27,7] + lDZ[27,7] + lstk[27,7] >= 0
 lATK: lDY[27,7] + lDZ[27,7] - ldz[27,7] - lstk[27,7] >= 0
 lATK: lDY[27,7] - ldz[27,7] >= 0
 lATK: lDY[27,8] - lDZ[27,8] = 0
 lATK: ldy[27,8] - ldz[27,8] = 0
 lATK: lDY[27,9] - lDZ[27,9] = 0
 lATK: ldy[27,9] - ldz[27,9] = 0
 lATK: lDY[27,10] - lDZ[27,10] = 0
 lATK: ldy[27,10] - ldz[27,10] = 0
 lATK: lDY[27,11] - lDZ[27,11] = 0
 lATK: ldy[27,11] - ldz[27,11] = 0
 lATK: lDY[27,12] - lDZ[27,12] = 0
 lATK: ldy[27,12] - ldz[27,12] = 0
 lATK: lDY[27,13] - lDZ[27,13] = 0
 lATK: ldy[27,13] - ldz[27,13] = 0
 lATK: lDY[27,14] - lDZ[27,14] = 0
 lATK: ldy[27,14] - ldz[27,14] = 0
 lATK: lDY[27,15] - lDZ[27,15] = 0
 lATK: ldy[27,15] - ldz[27,15] = 0
 lATK: ldy[28,0] - ldz[28,0] = 0
 lATK: lDY[28,0] - lDZ[28,0] + lstk[28,0] >= 0
 lATK: - lDY[28,0] + lDZ[28,0] + lstk[28,0] >= 0
 lATK: lDY[28,0] + lDZ[28,0] - ldz[28,0] - lstk[28,0] >= 0
 lATK: lDY[28,0] - ldz[28,0] >= 0
 lATK: ldy[28,1] - ldz[28,1] = 0
 lATK: lDY[28,1] - lDZ[28,1] + lstk[28,1] >= 0
 lATK: - lDY[28,1] + lDZ[28,1] + lstk[28,1] >= 0
 lATK: lDY[28,1] + lDZ[28,1] - ldz[28,1] - lstk[28,1] >= 0
 lATK: lDY[28,1] - ldz[28,1] >= 0
 lATK: ldy[28,2] - ldz[28,2] = 0
 lATK: lDY[28,2] - lDZ[28,2] + lstk[28,2] >= 0
 lATK: - lDY[28,2] + lDZ[28,2] + lstk[28,2] >= 0
 lATK: lDY[28,2] + lDZ[28,2] - ldz[28,2] - lstk[28,2] >= 0
 lATK: lDY[28,2] - ldz[28,2] >= 0
 lATK: ldy[28,3] - ldz[28,3] = 0
 lATK: lDY[28,3] - lDZ[28,3] + lstk[28,3] >= 0
 lATK: - lDY[28,3] + lDZ[28,3] + lstk[28,3] >= 0
 lATK: lDY[28,3] + lDZ[28,3] - ldz[28,3] - lstk[28,3] >= 0
 lATK: lDY[28,3] - ldz[28,3] >= 0
 lATK: ldy[28,4] - ldz[28,4] = 0
 lATK: lDY[28,4] - lDZ[28,4] + lstk[28,4] >= 0
 lATK: - lDY[28,4] + lDZ[28,4] + lstk[28,4] >= 0
 lATK: lDY[28,4] + lDZ[28,4] - ldz[28,4] - lstk[28,4] >= 0
 lATK: lDY[28,4] - ldz[28,4] >= 0
 lATK: ldy[28,5] - ldz[28,5] = 0
 lATK: lDY[28,5] - lDZ[28,5] + lstk[28,5] >= 0
 lATK: - lDY[28,5] + lDZ[28,5] + lstk[28,5] >= 0
 lATK: lDY[28,5] + lDZ[28,5] - ldz[28,5] - lstk[28,5] >= 0
 lATK: lDY[28,5] - ldz[28,5] >= 0
 lATK: ldy[28,6] - ldz[28,6] = 0
 lATK: lDY[28,6] - lDZ[28,6] + lstk[28,6] >= 0
 lATK: - lDY[28,6] + lDZ[28,6] + lstk[28,6] >= 0
 lATK: lDY[28,6] + lDZ[28,6] - ldz[28,6] - lstk[28,6] >= 0
 lATK: lDY[28,6] - ldz[28,6] >= 0
 lATK: ldy[28,7] - ldz[28,7] = 0
 lATK: lDY[28,7] - lDZ[28,7] + lstk[28,7] >= 0
 lATK: - lDY[28,7] + lDZ[28,7] + lstk[28,7] >= 0
 lATK: lDY[28,7] + lDZ[28,7] - ldz[28,7] - lstk[28,7] >= 0
 lATK: lDY[28,7] - ldz[28,7] >= 0
 lATK: lDY[28,8] - lDZ[28,8] = 0
 lATK: ldy[28,8] - ldz[28,8] = 0
 lATK: lDY[28,9] - lDZ[28,9] = 0
 lATK: ldy[28,9] - ldz[28,9] = 0
 lATK: lDY[28,10] - lDZ[28,10] = 0
 lATK: ldy[28,10] - ldz[28,10] = 0
 lATK: lDY[28,11] - lDZ[28,11] = 0
 lATK: ldy[28,11] - ldz[28,11] = 0
 lATK: lDY[28,12] - lDZ[28,12] = 0
 lATK: ldy[28,12] - ldz[28,12] = 0
 lATK: lDY[28,13] - lDZ[28,13] = 0
 lATK: ldy[28,13] - ldz[28,13] = 0
 lATK: lDY[28,14] - lDZ[28,14] = 0
 lATK: ldy[28,14] - ldz[28,14] = 0
 lATK: lDY[28,15] - lDZ[28,15] = 0
 lATK: ldy[28,15] - ldz[28,15] = 0
 lATK: ldy[29,0] - ldz[29,0] = 0
 lATK: lDY[29,0] - lDZ[29,0] + lstk[29,0] >= 0
 lATK: - lDY[29,0] + lDZ[29,0] + lstk[29,0] >= 0
 lATK: lDY[29,0] + lDZ[29,0] - ldz[29,0] - lstk[29,0] >= 0
 lATK: lDY[29,0] - ldz[29,0] >= 0
 lATK: ldy[29,1] - ldz[29,1] = 0
 lATK: lDY[29,1] - lDZ[29,1] + lstk[29,1] >= 0
 lATK: - lDY[29,1] + lDZ[29,1] + lstk[29,1] >= 0
 lATK: lDY[29,1] + lDZ[29,1] - ldz[29,1] - lstk[29,1] >= 0
 lATK: lDY[29,1] - ldz[29,1] >= 0
 lATK: ldy[29,2] - ldz[29,2] = 0
 lATK: lDY[29,2] - lDZ[29,2] + lstk[29,2] >= 0
 lATK: - lDY[29,2] + lDZ[29,2] + lstk[29,2] >= 0
 lATK: lDY[29,2] + lDZ[29,2] - ldz[29,2] - lstk[29,2] >= 0
 lATK: lDY[29,2] - ldz[29,2] >= 0
 lATK: ldy[29,3] - ldz[29,3] = 0
 lATK: lDY[29,3] - lDZ[29,3] + lstk[29,3] >= 0
 lATK: - lDY[29,3] + lDZ[29,3] + lstk[29,3] >= 0
 lATK: lDY[29,3] + lDZ[29,3] - ldz[29,3] - lstk[29,3] >= 0
 lATK: lDY[29,3] - ldz[29,3] >= 0
 lATK: ldy[29,4] - ldz[29,4] = 0
 lATK: lDY[29,4] - lDZ[29,4] + lstk[29,4] >= 0
 lATK: - lDY[29,4] + lDZ[29,4] + lstk[29,4] >= 0
 lATK: lDY[29,4] + lDZ[29,4] - ldz[29,4] - lstk[29,4] >= 0
 lATK: lDY[29,4] - ldz[29,4] >= 0
 lATK: ldy[29,5] - ldz[29,5] = 0
 lATK: lDY[29,5] - lDZ[29,5] + lstk[29,5] >= 0
 lATK: - lDY[29,5] + lDZ[29,5] + lstk[29,5] >= 0
 lATK: lDY[29,5] + lDZ[29,5] - ldz[29,5] - lstk[29,5] >= 0
 lATK: lDY[29,5] - ldz[29,5] >= 0
 lATK: ldy[29,6] - ldz[29,6] = 0
 lATK: lDY[29,6] - lDZ[29,6] + lstk[29,6] >= 0
 lATK: - lDY[29,6] + lDZ[29,6] + lstk[29,6] >= 0
 lATK: lDY[29,6] + lDZ[29,6] - ldz[29,6] - lstk[29,6] >= 0
 lATK: lDY[29,6] - ldz[29,6] >= 0
 lATK: ldy[29,7] - ldz[29,7] = 0
 lATK: lDY[29,7] - lDZ[29,7] + lstk[29,7] >= 0
 lATK: - lDY[29,7] + lDZ[29,7] + lstk[29,7] >= 0
 lATK: lDY[29,7] + lDZ[29,7] - ldz[29,7] - lstk[29,7] >= 0
 lATK: lDY[29,7] - ldz[29,7] >= 0
 lATK: lDY[29,8] - lDZ[29,8] = 0
 lATK: ldy[29,8] - ldz[29,8] = 0
 lATK: lDY[29,9] - lDZ[29,9] = 0
 lATK: ldy[29,9] - ldz[29,9] = 0
 lATK: lDY[29,10] - lDZ[29,10] = 0
 lATK: ldy[29,10] - ldz[29,10] = 0
 lATK: lDY[29,11] - lDZ[29,11] = 0
 lATK: ldy[29,11] - ldz[29,11] = 0
 lATK: lDY[29,12] - lDZ[29,12] = 0
 lATK: ldy[29,12] - ldz[29,12] = 0
 lATK: lDY[29,13] - lDZ[29,13] = 0
 lATK: ldy[29,13] - ldz[29,13] = 0
 lATK: lDY[29,14] - lDZ[29,14] = 0
 lATK: ldy[29,14] - ldz[29,14] = 0
 lATK: lDY[29,15] - lDZ[29,15] = 0
 lATK: ldy[29,15] - ldz[29,15] = 0
 lATK: ldy[30,0] - ldz[30,0] = 0
 lATK: lDY[30,0] - lDZ[30,0] + lstk[30,0] >= 0
 lATK: - lDY[30,0] + lDZ[30,0] + lstk[30,0] >= 0
 lATK: lDY[30,0] + lDZ[30,0] - ldz[30,0] - lstk[30,0] >= 0
 lATK: lDY[30,0] - ldz[30,0] >= 0
 lATK: ldy[30,1] - ldz[30,1] = 0
 lATK: lDY[30,1] - lDZ[30,1] + lstk[30,1] >= 0
 lATK: - lDY[30,1] + lDZ[30,1] + lstk[30,1] >= 0
 lATK: lDY[30,1] + lDZ[30,1] - ldz[30,1] - lstk[30,1] >= 0
 lATK: lDY[30,1] - ldz[30,1] >= 0
 lATK: ldy[30,2] - ldz[30,2] = 0
 lATK: lDY[30,2] - lDZ[30,2] + lstk[30,2] >= 0
 lATK: - lDY[30,2] + lDZ[30,2] + lstk[30,2] >= 0
 lATK: lDY[30,2] + lDZ[30,2] - ldz[30,2] - lstk[30,2] >= 0
 lATK: lDY[30,2] - ldz[30,2] >= 0
 lATK: ldy[30,3] - ldz[30,3] = 0
 lATK: lDY[30,3] - lDZ[30,3] + lstk[30,3] >= 0
 lATK: - lDY[30,3] + lDZ[30,3] + lstk[30,3] >= 0
 lATK: lDY[30,3] + lDZ[30,3] - ldz[30,3] - lstk[30,3] >= 0
 lATK: lDY[30,3] - ldz[30,3] >= 0
 lATK: ldy[30,4] - ldz[30,4] = 0
 lATK: lDY[30,4] - lDZ[30,4] + lstk[30,4] >= 0
 lATK: - lDY[30,4] + lDZ[30,4] + lstk[30,4] >= 0
 lATK: lDY[30,4] + lDZ[30,4] - ldz[30,4] - lstk[30,4] >= 0
 lATK: lDY[30,4] - ldz[30,4] >= 0
 lATK: ldy[30,5] - ldz[30,5] = 0
 lATK: lDY[30,5] - lDZ[30,5] + lstk[30,5] >= 0
 lATK: - lDY[30,5] + lDZ[30,5] + lstk[30,5] >= 0
 lATK: lDY[30,5] + lDZ[30,5] - ldz[30,5] - lstk[30,5] >= 0
 lATK: lDY[30,5] - ldz[30,5] >= 0
 lATK: ldy[30,6] - ldz[30,6] = 0
 lATK: lDY[30,6] - lDZ[30,6] + lstk[30,6] >= 0
 lATK: - lDY[30,6] + lDZ[30,6] + lstk[30,6] >= 0
 lATK: lDY[30,6] + lDZ[30,6] - ldz[30,6] - lstk[30,6] >= 0
 lATK: lDY[30,6] - ldz[30,6] >= 0
 lATK: ldy[30,7] - ldz[30,7] = 0
 lATK: lDY[30,7] - lDZ[30,7] + lstk[30,7] >= 0
 lATK: - lDY[30,7] + lDZ[30,7] + lstk[30,7] >= 0
 lATK: lDY[30,7] + lDZ[30,7] - ldz[30,7] - lstk[30,7] >= 0
 lATK: lDY[30,7] - ldz[30,7] >= 0
 lATK: lDY[30,8] - lDZ[30,8] = 0
 lATK: ldy[30,8] - ldz[30,8] = 0
 lATK: lDY[30,9] - lDZ[30,9] = 0
 lATK: ldy[30,9] - ldz[30,9] = 0
 lATK: lDY[30,10] - lDZ[30,10] = 0
 lATK: ldy[30,10] - ldz[30,10] = 0
 lATK: lDY[30,11] - lDZ[30,11] = 0
 lATK: ldy[30,11] - ldz[30,11] = 0
 lATK: lDY[30,12] - lDZ[30,12] = 0
 lATK: ldy[30,12] - ldz[30,12] = 0
 lATK: lDY[30,13] - lDZ[30,13] = 0
 lATK: ldy[30,13] - ldz[30,13] = 0
 lATK: lDY[30,14] - lDZ[30,14] = 0
 lATK: ldy[30,14] - ldz[30,14] = 0
 lATK: lDY[30,15] - lDZ[30,15] = 0
 lATK: ldy[30,15] - ldz[30,15] = 0
 lATK: ldy[31,0] - ldz[31,0] = 0
 lATK: lDY[31,0] - lDZ[31,0] + lstk[31,0] >= 0
 lATK: - lDY[31,0] + lDZ[31,0] + lstk[31,0] >= 0
 lATK: lDY[31,0] + lDZ[31,0] - ldz[31,0] - lstk[31,0] >= 0
 lATK: lDY[31,0] - ldz[31,0] >= 0
 lATK: ldy[31,1] - ldz[31,1] = 0
 lATK: lDY[31,1] - lDZ[31,1] + lstk[31,1] >= 0
 lATK: - lDY[31,1] + lDZ[31,1] + lstk[31,1] >= 0
 lATK: lDY[31,1] + lDZ[31,1] - ldz[31,1] - lstk[31,1] >= 0
 lATK: lDY[31,1] - ldz[31,1] >= 0
 lATK: ldy[31,2] - ldz[31,2] = 0
 lATK: lDY[31,2] - lDZ[31,2] + lstk[31,2] >= 0
 lATK: - lDY[31,2] + lDZ[31,2] + lstk[31,2] >= 0
 lATK: lDY[31,2] + lDZ[31,2] - ldz[31,2] - lstk[31,2] >= 0
 lATK: lDY[31,2] - ldz[31,2] >= 0
 lATK: ldy[31,3] - ldz[31,3] = 0
 lATK: lDY[31,3] - lDZ[31,3] + lstk[31,3] >= 0
 lATK: - lDY[31,3] + lDZ[31,3] + lstk[31,3] >= 0
 lATK: lDY[31,3] + lDZ[31,3] - ldz[31,3] - lstk[31,3] >= 0
 lATK: lDY[31,3] - ldz[31,3] >= 0
 lATK: ldy[31,4] - ldz[31,4] = 0
 lATK: lDY[31,4] - lDZ[31,4] + lstk[31,4] >= 0
 lATK: - lDY[31,4] + lDZ[31,4] + lstk[31,4] >= 0
 lATK: lDY[31,4] + lDZ[31,4] - ldz[31,4] - lstk[31,4] >= 0
 lATK: lDY[31,4] - ldz[31,4] >= 0
 lATK: ldy[31,5] - ldz[31,5] = 0
 lATK: lDY[31,5] - lDZ[31,5] + lstk[31,5] >= 0
 lATK: - lDY[31,5] + lDZ[31,5] + lstk[31,5] >= 0
 lATK: lDY[31,5] + lDZ[31,5] - ldz[31,5] - lstk[31,5] >= 0
 lATK: lDY[31,5] - ldz[31,5] >= 0
 lATK: ldy[31,6] - ldz[31,6] = 0
 lATK: lDY[31,6] - lDZ[31,6] + lstk[31,6] >= 0
 lATK: - lDY[31,6] + lDZ[31,6] + lstk[31,6] >= 0
 lATK: lDY[31,6] + lDZ[31,6] - ldz[31,6] - lstk[31,6] >= 0
 lATK: lDY[31,6] - ldz[31,6] >= 0
 lATK: ldy[31,7] - ldz[31,7] = 0
 lATK: lDY[31,7] - lDZ[31,7] + lstk[31,7] >= 0
 lATK: - lDY[31,7] + lDZ[31,7] + lstk[31,7] >= 0
 lATK: lDY[31,7] + lDZ[31,7] - ldz[31,7] - lstk[31,7] >= 0
 lATK: lDY[31,7] - ldz[31,7] >= 0
 lATK: lDY[31,8] - lDZ[31,8] = 0
 lATK: ldy[31,8] - ldz[31,8] = 0
 lATK: lDY[31,9] - lDZ[31,9] = 0
 lATK: ldy[31,9] - ldz[31,9] = 0
 lATK: lDY[31,10] - lDZ[31,10] = 0
 lATK: ldy[31,10] - ldz[31,10] = 0
 lATK: lDY[31,11] - lDZ[31,11] = 0
 lATK: ldy[31,11] - ldz[31,11] = 0
 lATK: lDY[31,12] - lDZ[31,12] = 0
 lATK: ldy[31,12] - ldz[31,12] = 0
 lATK: lDY[31,13] - lDZ[31,13] = 0
 lATK: ldy[31,13] - ldz[31,13] = 0
 lATK: lDY[31,14] - lDZ[31,14] = 0
 lATK: ldy[31,14] - ldz[31,14] = 0
 lATK: lDY[31,15] - lDZ[31,15] = 0
 lATK: ldy[31,15] - ldz[31,15] = 0
 lATK: ldy[32,0] - ldz[32,0] = 0
 lATK: lDY[32,0] - lDZ[32,0] + lstk[32,0] >= 0
 lATK: - lDY[32,0] + lDZ[32,0] + lstk[32,0] >= 0
 lATK: lDY[32,0] + lDZ[32,0] - ldz[32,0] - lstk[32,0] >= 0
 lATK: lDY[32,0] - ldz[32,0] >= 0
 lATK: ldy[32,1] - ldz[32,1] = 0
 lATK: lDY[32,1] - lDZ[32,1] + lstk[32,1] >= 0
 lATK: - lDY[32,1] + lDZ[32,1] + lstk[32,1] >= 0
 lATK: lDY[32,1] + lDZ[32,1] - ldz[32,1] - lstk[32,1] >= 0
 lATK: lDY[32,1] - ldz[32,1] >= 0
 lATK: ldy[32,2] - ldz[32,2] = 0
 lATK: lDY[32,2] - lDZ[32,2] + lstk[32,2] >= 0
 lATK: - lDY[32,2] + lDZ[32,2] + lstk[32,2] >= 0
 lATK: lDY[32,2] + lDZ[32,2] - ldz[32,2] - lstk[32,2] >= 0
 lATK: lDY[32,2] - ldz[32,2] >= 0
 lATK: ldy[32,3] - ldz[32,3] = 0
 lATK: lDY[32,3] - lDZ[32,3] + lstk[32,3] >= 0
 lATK: - lDY[32,3] + lDZ[32,3] + lstk[32,3] >= 0
 lATK: lDY[32,3] + lDZ[32,3] - ldz[32,3] - lstk[32,3] >= 0
 lATK: lDY[32,3] - ldz[32,3] >= 0
 lATK: ldy[32,4] - ldz[32,4] = 0
 lATK: lDY[32,4] - lDZ[32,4] + lstk[32,4] >= 0
 lATK: - lDY[32,4] + lDZ[32,4] + lstk[32,4] >= 0
 lATK: lDY[32,4] + lDZ[32,4] - ldz[32,4] - lstk[32,4] >= 0
 lATK: lDY[32,4] - ldz[32,4] >= 0
 lATK: ldy[32,5] - ldz[32,5] = 0
 lATK: lDY[32,5] - lDZ[32,5] + lstk[32,5] >= 0
 lATK: - lDY[32,5] + lDZ[32,5] + lstk[32,5] >= 0
 lATK: lDY[32,5] + lDZ[32,5] - ldz[32,5] - lstk[32,5] >= 0
 lATK: lDY[32,5] - ldz[32,5] >= 0
 lATK: ldy[32,6] - ldz[32,6] = 0
 lATK: lDY[32,6] - lDZ[32,6] + lstk[32,6] >= 0
 lATK: - lDY[32,6] + lDZ[32,6] + lstk[32,6] >= 0
 lATK: lDY[32,6] + lDZ[32,6] - ldz[32,6] - lstk[32,6] >= 0
 lATK: lDY[32,6] - ldz[32,6] >= 0
 lATK: ldy[32,7] - ldz[32,7] = 0
 lATK: lDY[32,7] - lDZ[32,7] + lstk[32,7] >= 0
 lATK: - lDY[32,7] + lDZ[32,7] + lstk[32,7] >= 0
 lATK: lDY[32,7] + lDZ[32,7] - ldz[32,7] - lstk[32,7] >= 0
 lATK: lDY[32,7] - ldz[32,7] >= 0
 lATK: lDY[32,8] - lDZ[32,8] = 0
 lATK: ldy[32,8] - ldz[32,8] = 0
 lATK: lDY[32,9] - lDZ[32,9] = 0
 lATK: ldy[32,9] - ldz[32,9] = 0
 lATK: lDY[32,10] - lDZ[32,10] = 0
 lATK: ldy[32,10] - ldz[32,10] = 0
 lATK: lDY[32,11] - lDZ[32,11] = 0
 lATK: ldy[32,11] - ldz[32,11] = 0
 lATK: lDY[32,12] - lDZ[32,12] = 0
 lATK: ldy[32,12] - ldz[32,12] = 0
 lATK: lDY[32,13] - lDZ[32,13] = 0
 lATK: ldy[32,13] - ldz[32,13] = 0
 lATK: lDY[32,14] - lDZ[32,14] = 0
 lATK: ldy[32,14] - ldz[32,14] = 0
 lATK: lDY[32,15] - lDZ[32,15] = 0
 lATK: ldy[32,15] - ldz[32,15] = 0
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
 lSR: lDZ[24,0] - lDW[24,0] = 0
 lSR: ldz[24,0] - ldw[24,0] = 0
 lSR: lDZ[24,1] - lDW[24,1] = 0
 lSR: ldz[24,1] - ldw[24,1] = 0
 lSR: lDZ[24,2] - lDW[24,2] = 0
 lSR: ldz[24,2] - ldw[24,2] = 0
 lSR: lDZ[24,3] - lDW[24,3] = 0
 lSR: ldz[24,3] - ldw[24,3] = 0
 lSR: lDZ[24,4] - lDW[24,5] = 0
 lSR: ldz[24,4] - ldw[24,5] = 0
 lSR: lDZ[24,5] - lDW[24,6] = 0
 lSR: ldz[24,5] - ldw[24,6] = 0
 lSR: lDZ[24,6] - lDW[24,7] = 0
 lSR: ldz[24,6] - ldw[24,7] = 0
 lSR: lDZ[24,7] - lDW[24,4] = 0
 lSR: ldz[24,7] - ldw[24,4] = 0
 lSR: lDZ[24,8] - lDW[24,10] = 0
 lSR: ldz[24,8] - ldw[24,10] = 0
 lSR: lDZ[24,9] - lDW[24,11] = 0
 lSR: ldz[24,9] - ldw[24,11] = 0
 lSR: lDZ[24,10] - lDW[24,8] = 0
 lSR: ldz[24,10] - ldw[24,8] = 0
 lSR: lDZ[24,11] - lDW[24,9] = 0
 lSR: ldz[24,11] - ldw[24,9] = 0
 lSR: lDZ[24,12] - lDW[24,15] = 0
 lSR: ldz[24,12] - ldw[24,15] = 0
 lSR: lDZ[24,13] - lDW[24,12] = 0
 lSR: ldz[24,13] - ldw[24,12] = 0
 lSR: lDZ[24,14] - lDW[24,13] = 0
 lSR: ldz[24,14] - ldw[24,13] = 0
 lSR: lDZ[24,15] - lDW[24,14] = 0
 lSR: ldz[24,15] - ldw[24,14] = 0
 lSR: lDZ[25,0] - lDW[25,0] = 0
 lSR: ldz[25,0] - ldw[25,0] = 0
 lSR: lDZ[25,1] - lDW[25,1] = 0
 lSR: ldz[25,1] - ldw[25,1] = 0
 lSR: lDZ[25,2] - lDW[25,2] = 0
 lSR: ldz[25,2] - ldw[25,2] = 0
 lSR: lDZ[25,3] - lDW[25,3] = 0
 lSR: ldz[25,3] - ldw[25,3] = 0
 lSR: lDZ[25,4] - lDW[25,5] = 0
 lSR: ldz[25,4] - ldw[25,5] = 0
 lSR: lDZ[25,5] - lDW[25,6] = 0
 lSR: ldz[25,5] - ldw[25,6] = 0
 lSR: lDZ[25,6] - lDW[25,7] = 0
 lSR: ldz[25,6] - ldw[25,7] = 0
 lSR: lDZ[25,7] - lDW[25,4] = 0
 lSR: ldz[25,7] - ldw[25,4] = 0
 lSR: lDZ[25,8] - lDW[25,10] = 0
 lSR: ldz[25,8] - ldw[25,10] = 0
 lSR: lDZ[25,9] - lDW[25,11] = 0
 lSR: ldz[25,9] - ldw[25,11] = 0
 lSR: lDZ[25,10] - lDW[25,8] = 0
 lSR: ldz[25,10] - ldw[25,8] = 0
 lSR: lDZ[25,11] - lDW[25,9] = 0
 lSR: ldz[25,11] - ldw[25,9] = 0
 lSR: lDZ[25,12] - lDW[25,15] = 0
 lSR: ldz[25,12] - ldw[25,15] = 0
 lSR: lDZ[25,13] - lDW[25,12] = 0
 lSR: ldz[25,13] - ldw[25,12] = 0
 lSR: lDZ[25,14] - lDW[25,13] = 0
 lSR: ldz[25,14] - ldw[25,13] = 0
 lSR: lDZ[25,15] - lDW[25,14] = 0
 lSR: ldz[25,15] - ldw[25,14] = 0
 lSR: lDZ[26,0] - lDW[26,0] = 0
 lSR: ldz[26,0] - ldw[26,0] = 0
 lSR: lDZ[26,1] - lDW[26,1] = 0
 lSR: ldz[26,1] - ldw[26,1] = 0
 lSR: lDZ[26,2] - lDW[26,2] = 0
 lSR: ldz[26,2] - ldw[26,2] = 0
 lSR: lDZ[26,3] - lDW[26,3] = 0
 lSR: ldz[26,3] - ldw[26,3] = 0
 lSR: lDZ[26,4] - lDW[26,5] = 0
 lSR: ldz[26,4] - ldw[26,5] = 0
 lSR: lDZ[26,5] - lDW[26,6] = 0
 lSR: ldz[26,5] - ldw[26,6] = 0
 lSR: lDZ[26,6] - lDW[26,7] = 0
 lSR: ldz[26,6] - ldw[26,7] = 0
 lSR: lDZ[26,7] - lDW[26,4] = 0
 lSR: ldz[26,7] - ldw[26,4] = 0
 lSR: lDZ[26,8] - lDW[26,10] = 0
 lSR: ldz[26,8] - ldw[26,10] = 0
 lSR: lDZ[26,9] - lDW[26,11] = 0
 lSR: ldz[26,9] - ldw[26,11] = 0
 lSR: lDZ[26,10] - lDW[26,8] = 0
 lSR: ldz[26,10] - ldw[26,8] = 0
 lSR: lDZ[26,11] - lDW[26,9] = 0
 lSR: ldz[26,11] - ldw[26,9] = 0
 lSR: lDZ[26,12] - lDW[26,15] = 0
 lSR: ldz[26,12] - ldw[26,15] = 0
 lSR: lDZ[26,13] - lDW[26,12] = 0
 lSR: ldz[26,13] - ldw[26,12] = 0
 lSR: lDZ[26,14] - lDW[26,13] = 0
 lSR: ldz[26,14] - ldw[26,13] = 0
 lSR: lDZ[26,15] - lDW[26,14] = 0
 lSR: ldz[26,15] - ldw[26,14] = 0
 lSR: lDZ[27,0] - lDW[27,0] = 0
 lSR: ldz[27,0] - ldw[27,0] = 0
 lSR: lDZ[27,1] - lDW[27,1] = 0
 lSR: ldz[27,1] - ldw[27,1] = 0
 lSR: lDZ[27,2] - lDW[27,2] = 0
 lSR: ldz[27,2] - ldw[27,2] = 0
 lSR: lDZ[27,3] - lDW[27,3] = 0
 lSR: ldz[27,3] - ldw[27,3] = 0
 lSR: lDZ[27,4] - lDW[27,5] = 0
 lSR: ldz[27,4] - ldw[27,5] = 0
 lSR: lDZ[27,5] - lDW[27,6] = 0
 lSR: ldz[27,5] - ldw[27,6] = 0
 lSR: lDZ[27,6] - lDW[27,7] = 0
 lSR: ldz[27,6] - ldw[27,7] = 0
 lSR: lDZ[27,7] - lDW[27,4] = 0
 lSR: ldz[27,7] - ldw[27,4] = 0
 lSR: lDZ[27,8] - lDW[27,10] = 0
 lSR: ldz[27,8] - ldw[27,10] = 0
 lSR: lDZ[27,9] - lDW[27,11] = 0
 lSR: ldz[27,9] - ldw[27,11] = 0
 lSR: lDZ[27,10] - lDW[27,8] = 0
 lSR: ldz[27,10] - ldw[27,8] = 0
 lSR: lDZ[27,11] - lDW[27,9] = 0
 lSR: ldz[27,11] - ldw[27,9] = 0
 lSR: lDZ[27,12] - lDW[27,15] = 0
 lSR: ldz[27,12] - ldw[27,15] = 0
 lSR: lDZ[27,13] - lDW[27,12] = 0
 lSR: ldz[27,13] - ldw[27,12] = 0
 lSR: lDZ[27,14] - lDW[27,13] = 0
 lSR: ldz[27,14] - ldw[27,13] = 0
 lSR: lDZ[27,15] - lDW[27,14] = 0
 lSR: ldz[27,15] - ldw[27,14] = 0
 lSR: lDZ[28,0] - lDW[28,0] = 0
 lSR: ldz[28,0] - ldw[28,0] = 0
 lSR: lDZ[28,1] - lDW[28,1] = 0
 lSR: ldz[28,1] - ldw[28,1] = 0
 lSR: lDZ[28,2] - lDW[28,2] = 0
 lSR: ldz[28,2] - ldw[28,2] = 0
 lSR: lDZ[28,3] - lDW[28,3] = 0
 lSR: ldz[28,3] - ldw[28,3] = 0
 lSR: lDZ[28,4] - lDW[28,5] = 0
 lSR: ldz[28,4] - ldw[28,5] = 0
 lSR: lDZ[28,5] - lDW[28,6] = 0
 lSR: ldz[28,5] - ldw[28,6] = 0
 lSR: lDZ[28,6] - lDW[28,7] = 0
 lSR: ldz[28,6] - ldw[28,7] = 0
 lSR: lDZ[28,7] - lDW[28,4] = 0
 lSR: ldz[28,7] - ldw[28,4] = 0
 lSR: lDZ[28,8] - lDW[28,10] = 0
 lSR: ldz[28,8] - ldw[28,10] = 0
 lSR: lDZ[28,9] - lDW[28,11] = 0
 lSR: ldz[28,9] - ldw[28,11] = 0
 lSR: lDZ[28,10] - lDW[28,8] = 0
 lSR: ldz[28,10] - ldw[28,8] = 0
 lSR: lDZ[28,11] - lDW[28,9] = 0
 lSR: ldz[28,11] - ldw[28,9] = 0
 lSR: lDZ[28,12] - lDW[28,15] = 0
 lSR: ldz[28,12] - ldw[28,15] = 0
 lSR: lDZ[28,13] - lDW[28,12] = 0
 lSR: ldz[28,13] - ldw[28,12] = 0
 lSR: lDZ[28,14] - lDW[28,13] = 0
 lSR: ldz[28,14] - ldw[28,13] = 0
 lSR: lDZ[28,15] - lDW[28,14] = 0
 lSR: ldz[28,15] - ldw[28,14] = 0
 lSR: lDZ[29,0] - lDW[29,0] = 0
 lSR: ldz[29,0] - ldw[29,0] = 0
 lSR: lDZ[29,1] - lDW[29,1] = 0
 lSR: ldz[29,1] - ldw[29,1] = 0
 lSR: lDZ[29,2] - lDW[29,2] = 0
 lSR: ldz[29,2] - ldw[29,2] = 0
 lSR: lDZ[29,3] - lDW[29,3] = 0
 lSR: ldz[29,3] - ldw[29,3] = 0
 lSR: lDZ[29,4] - lDW[29,5] = 0
 lSR: ldz[29,4] - ldw[29,5] = 0
 lSR: lDZ[29,5] - lDW[29,6] = 0
 lSR: ldz[29,5] - ldw[29,6] = 0
 lSR: lDZ[29,6] - lDW[29,7] = 0
 lSR: ldz[29,6] - ldw[29,7] = 0
 lSR: lDZ[29,7] - lDW[29,4] = 0
 lSR: ldz[29,7] - ldw[29,4] = 0
 lSR: lDZ[29,8] - lDW[29,10] = 0
 lSR: ldz[29,8] - ldw[29,10] = 0
 lSR: lDZ[29,9] - lDW[29,11] = 0
 lSR: ldz[29,9] - ldw[29,11] = 0
 lSR: lDZ[29,10] - lDW[29,8] = 0
 lSR: ldz[29,10] - ldw[29,8] = 0
 lSR: lDZ[29,11] - lDW[29,9] = 0
 lSR: ldz[29,11] - ldw[29,9] = 0
 lSR: lDZ[29,12] - lDW[29,15] = 0
 lSR: ldz[29,12] - ldw[29,15] = 0
 lSR: lDZ[29,13] - lDW[29,12] = 0
 lSR: ldz[29,13] - ldw[29,12] = 0
 lSR: lDZ[29,14] - lDW[29,13] = 0
 lSR: ldz[29,14] - ldw[29,13] = 0
 lSR: lDZ[29,15] - lDW[29,14] = 0
 lSR: ldz[29,15] - ldw[29,14] = 0
 lSR: lDZ[30,0] - lDW[30,0] = 0
 lSR: ldz[30,0] - ldw[30,0] = 0
 lSR: lDZ[30,1] - lDW[30,1] = 0
 lSR: ldz[30,1] - ldw[30,1] = 0
 lSR: lDZ[30,2] - lDW[30,2] = 0
 lSR: ldz[30,2] - ldw[30,2] = 0
 lSR: lDZ[30,3] - lDW[30,3] = 0
 lSR: ldz[30,3] - ldw[30,3] = 0
 lSR: lDZ[30,4] - lDW[30,5] = 0
 lSR: ldz[30,4] - ldw[30,5] = 0
 lSR: lDZ[30,5] - lDW[30,6] = 0
 lSR: ldz[30,5] - ldw[30,6] = 0
 lSR: lDZ[30,6] - lDW[30,7] = 0
 lSR: ldz[30,6] - ldw[30,7] = 0
 lSR: lDZ[30,7] - lDW[30,4] = 0
 lSR: ldz[30,7] - ldw[30,4] = 0
 lSR: lDZ[30,8] - lDW[30,10] = 0
 lSR: ldz[30,8] - ldw[30,10] = 0
 lSR: lDZ[30,9] - lDW[30,11] = 0
 lSR: ldz[30,9] - ldw[30,11] = 0
 lSR: lDZ[30,10] - lDW[30,8] = 0
 lSR: ldz[30,10] - ldw[30,8] = 0
 lSR: lDZ[30,11] - lDW[30,9] = 0
 lSR: ldz[30,11] - ldw[30,9] = 0
 lSR: lDZ[30,12] - lDW[30,15] = 0
 lSR: ldz[30,12] - ldw[30,15] = 0
 lSR: lDZ[30,13] - lDW[30,12] = 0
 lSR: ldz[30,13] - ldw[30,12] = 0
 lSR: lDZ[30,14] - lDW[30,13] = 0
 lSR: ldz[30,14] - ldw[30,13] = 0
 lSR: lDZ[30,15] - lDW[30,14] = 0
 lSR: ldz[30,15] - ldw[30,14] = 0
 lSR: lDZ[31,0] - lDW[31,0] = 0
 lSR: ldz[31,0] - ldw[31,0] = 0
 lSR: lDZ[31,1] - lDW[31,1] = 0
 lSR: ldz[31,1] - ldw[31,1] = 0
 lSR: lDZ[31,2] - lDW[31,2] = 0
 lSR: ldz[31,2] - ldw[31,2] = 0
 lSR: lDZ[31,3] - lDW[31,3] = 0
 lSR: ldz[31,3] - ldw[31,3] = 0
 lSR: lDZ[31,4] - lDW[31,5] = 0
 lSR: ldz[31,4] - ldw[31,5] = 0
 lSR: lDZ[31,5] - lDW[31,6] = 0
 lSR: ldz[31,5] - ldw[31,6] = 0
 lSR: lDZ[31,6] - lDW[31,7] = 0
 lSR: ldz[31,6] - ldw[31,7] = 0
 lSR: lDZ[31,7] - lDW[31,4] = 0
 lSR: ldz[31,7] - ldw[31,4] = 0
 lSR: lDZ[31,8] - lDW[31,10] = 0
 lSR: ldz[31,8] - ldw[31,10] = 0
 lSR: lDZ[31,9] - lDW[31,11] = 0
 lSR: ldz[31,9] - ldw[31,11] = 0
 lSR: lDZ[31,10] - lDW[31,8] = 0
 lSR: ldz[31,10] - ldw[31,8] = 0
 lSR: lDZ[31,11] - lDW[31,9] = 0
 lSR: ldz[31,11] - ldw[31,9] = 0
 lSR: lDZ[31,12] - lDW[31,15] = 0
 lSR: ldz[31,12] - ldw[31,15] = 0
 lSR: lDZ[31,13] - lDW[31,12] = 0
 lSR: ldz[31,13] - ldw[31,12] = 0
 lSR: lDZ[31,14] - lDW[31,13] = 0
 lSR: ldz[31,14] - ldw[31,13] = 0
 lSR: lDZ[31,15] - lDW[31,14] = 0
 lSR: ldz[31,15] - ldw[31,14] = 0
 lSR: lDZ[32,0] - lDW[32,0] = 0
 lSR: ldz[32,0] - ldw[32,0] = 0
 lSR: lDZ[32,1] - lDW[32,1] = 0
 lSR: ldz[32,1] - ldw[32,1] = 0
 lSR: lDZ[32,2] - lDW[32,2] = 0
 lSR: ldz[32,2] - ldw[32,2] = 0
 lSR: lDZ[32,3] - lDW[32,3] = 0
 lSR: ldz[32,3] - ldw[32,3] = 0
 lSR: lDZ[32,4] - lDW[32,5] = 0
 lSR: ldz[32,4] - ldw[32,5] = 0
 lSR: lDZ[32,5] - lDW[32,6] = 0
 lSR: ldz[32,5] - ldw[32,6] = 0
 lSR: lDZ[32,6] - lDW[32,7] = 0
 lSR: ldz[32,6] - ldw[32,7] = 0
 lSR: lDZ[32,7] - lDW[32,4] = 0
 lSR: ldz[32,7] - ldw[32,4] = 0
 lSR: lDZ[32,8] - lDW[32,10] = 0
 lSR: ldz[32,8] - ldw[32,10] = 0
 lSR: lDZ[32,9] - lDW[32,11] = 0
 lSR: ldz[32,9] - ldw[32,11] = 0
 lSR: lDZ[32,10] - lDW[32,8] = 0
 lSR: ldz[32,10] - ldw[32,8] = 0
 lSR: lDZ[32,11] - lDW[32,9] = 0
 lSR: ldz[32,11] - ldw[32,9] = 0
 lSR: lDZ[32,12] - lDW[32,15] = 0
 lSR: ldz[32,12] - ldw[32,15] = 0
 lSR: lDZ[32,13] - lDW[32,12] = 0
 lSR: ldz[32,13] - ldw[32,12] = 0
 lSR: lDZ[32,14] - lDW[32,13] = 0
 lSR: ldz[32,14] - ldw[32,13] = 0
 lSR: lDZ[32,15] - lDW[32,14] = 0
 lSR: ldz[32,15] - ldw[32,14] = 0
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
 ldMC0: - ldx[20,4] + ldw[19,0] = 0
 ldMC0: - lDX[20,4] + lDW[19,0] = 0
 ldMC1: ldw[19,4] - ldw[19,8] >= 0
 ldMC1: - ldx[20,8] + ldw[19,4] >= 0
 ldMC1: - ldx[20,8] + ldw[19,4] - ldw[19,8] <= 0
 ldMC1: lDX[20,8] - lDW[19,4] + lDW[19,8] >= 0
 ldMC1: lDX[20,8] + lDW[19,4] - lDW[19,8] >= 0
 ldMC1: - lDX[20,8] + lDW[19,4] + lDW[19,8] >= 0
 ldMC1: - lDX[20,8] + lDW[19,4] - lDW[19,8] <= 0
 ldMC2: - ldx[20,4] + ldw[19,8] >= 0
 ldMC2: - ldx[20,12] + ldw[19,8] >= 0
 ldMC2: - ldx[20,4] - ldx[20,12] + ldw[19,8] <= 0
 ldMC2: lDX[20,4] + lDX[20,12] - lDW[19,8] >= 0
 ldMC2: - lDX[20,4] + lDX[20,12] + lDW[19,8] >= 0
 ldMC2: lDX[20,4] - lDX[20,12] + lDW[19,8] >= 0
 ldMC2: - lDX[20,4] - lDX[20,12] + lDW[19,8] <= 0
 ldMC3: - ldx[20,0] + ldw[19,12] >= 0
 ldMC3: - ldx[20,12] + ldw[19,12] >= 0
 ldMC3: - ldx[20,0] - ldx[20,12] + ldw[19,12] <= 0
 ldMC3: lDX[20,0] + lDX[20,12] - lDW[19,12] >= 0
 ldMC3: - lDX[20,0] + lDX[20,12] + lDW[19,12] >= 0
 ldMC3: lDX[20,0] - lDX[20,12] + lDW[19,12] >= 0
 ldMC3: - lDX[20,0] - lDX[20,12] + lDW[19,12] <= 0
 ldMC0: - ldx[20,5] + ldw[19,1] = 0
 ldMC0: - lDX[20,5] + lDW[19,1] = 0
 ldMC1: ldw[19,5] - ldw[19,9] >= 0
 ldMC1: - ldx[20,9] + ldw[19,5] >= 0
 ldMC1: - ldx[20,9] + ldw[19,5] - ldw[19,9] <= 0
 ldMC1: lDX[20,9] - lDW[19,5] + lDW[19,9] >= 0
 ldMC1: lDX[20,9] + lDW[19,5] - lDW[19,9] >= 0
 ldMC1: - lDX[20,9] + lDW[19,5] + lDW[19,9] >= 0
 ldMC1: - lDX[20,9] + lDW[19,5] - lDW[19,9] <= 0
 ldMC2: - ldx[20,5] + ldw[19,9] >= 0
 ldMC2: - ldx[20,13] + ldw[19,9] >= 0
 ldMC2: - ldx[20,5] - ldx[20,13] + ldw[19,9] <= 0
 ldMC2: lDX[20,5] + lDX[20,13] - lDW[19,9] >= 0
 ldMC2: - lDX[20,5] + lDX[20,13] + lDW[19,9] >= 0
 ldMC2: lDX[20,5] - lDX[20,13] + lDW[19,9] >= 0
 ldMC2: - lDX[20,5] - lDX[20,13] + lDW[19,9] <= 0
 ldMC3: - ldx[20,1] + ldw[19,13] >= 0
 ldMC3: - ldx[20,13] + ldw[19,13] >= 0
 ldMC3: - ldx[20,1] - ldx[20,13] + ldw[19,13] <= 0
 ldMC3: lDX[20,1] + lDX[20,13] - lDW[19,13] >= 0
 ldMC3: - lDX[20,1] + lDX[20,13] + lDW[19,13] >= 0
 ldMC3: lDX[20,1] - lDX[20,13] + lDW[19,13] >= 0
 ldMC3: - lDX[20,1] - lDX[20,13] + lDW[19,13] <= 0
 ldMC0: - ldx[20,6] + ldw[19,2] = 0
 ldMC0: - lDX[20,6] + lDW[19,2] = 0
 ldMC1: ldw[19,6] - ldw[19,10] >= 0
 ldMC1: - ldx[20,10] + ldw[19,6] >= 0
 ldMC1: - ldx[20,10] + ldw[19,6] - ldw[19,10] <= 0
 ldMC1: lDX[20,10] - lDW[19,6] + lDW[19,10] >= 0
 ldMC1: lDX[20,10] + lDW[19,6] - lDW[19,10] >= 0
 ldMC1: - lDX[20,10] + lDW[19,6] + lDW[19,10] >= 0
 ldMC1: - lDX[20,10] + lDW[19,6] - lDW[19,10] <= 0
 ldMC2: - ldx[20,6] + ldw[19,10] >= 0
 ldMC2: - ldx[20,14] + ldw[19,10] >= 0
 ldMC2: - ldx[20,6] - ldx[20,14] + ldw[19,10] <= 0
 ldMC2: lDX[20,6] + lDX[20,14] - lDW[19,10] >= 0
 ldMC2: - lDX[20,6] + lDX[20,14] + lDW[19,10] >= 0
 ldMC2: lDX[20,6] - lDX[20,14] + lDW[19,10] >= 0
 ldMC2: - lDX[20,6] - lDX[20,14] + lDW[19,10] <= 0
 ldMC3: - ldx[20,2] + ldw[19,14] >= 0
 ldMC3: - ldx[20,14] + ldw[19,14] >= 0
 ldMC3: - ldx[20,2] - ldx[20,14] + ldw[19,14] <= 0
 ldMC3: lDX[20,2] + lDX[20,14] - lDW[19,14] >= 0
 ldMC3: - lDX[20,2] + lDX[20,14] + lDW[19,14] >= 0
 ldMC3: lDX[20,2] - lDX[20,14] + lDW[19,14] >= 0
 ldMC3: - lDX[20,2] - lDX[20,14] + lDW[19,14] <= 0
 ldMC0: - ldx[20,7] + ldw[19,3] = 0
 ldMC0: - lDX[20,7] + lDW[19,3] = 0
 ldMC1: ldw[19,7] - ldw[19,11] >= 0
 ldMC1: - ldx[20,11] + ldw[19,7] >= 0
 ldMC1: - ldx[20,11] + ldw[19,7] - ldw[19,11] <= 0
 ldMC1: lDX[20,11] - lDW[19,7] + lDW[19,11] >= 0
 ldMC1: lDX[20,11] + lDW[19,7] - lDW[19,11] >= 0
 ldMC1: - lDX[20,11] + lDW[19,7] + lDW[19,11] >= 0
 ldMC1: - lDX[20,11] + lDW[19,7] - lDW[19,11] <= 0
 ldMC2: - ldx[20,7] + ldw[19,11] >= 0
 ldMC2: - ldx[20,15] + ldw[19,11] >= 0
 ldMC2: - ldx[20,7] - ldx[20,15] + ldw[19,11] <= 0
 ldMC2: lDX[20,7] + lDX[20,15] - lDW[19,11] >= 0
 ldMC2: - lDX[20,7] + lDX[20,15] + lDW[19,11] >= 0
 ldMC2: lDX[20,7] - lDX[20,15] + lDW[19,11] >= 0
 ldMC2: - lDX[20,7] - lDX[20,15] + lDW[19,11] <= 0
 ldMC3: - ldx[20,3] + ldw[19,15] >= 0
 ldMC3: - ldx[20,15] + ldw[19,15] >= 0
 ldMC3: - ldx[20,3] - ldx[20,15] + ldw[19,15] <= 0
 ldMC3: lDX[20,3] + lDX[20,15] - lDW[19,15] >= 0
 ldMC3: - lDX[20,3] + lDX[20,15] + lDW[19,15] >= 0
 ldMC3: lDX[20,3] - lDX[20,15] + lDW[19,15] >= 0
 ldMC3: - lDX[20,3] - lDX[20,15] + lDW[19,15] <= 0
 ldMC0: - ldx[21,4] + ldw[20,0] = 0
 ldMC0: - lDX[21,4] + lDW[20,0] = 0
 ldMC1: ldw[20,4] - ldw[20,8] >= 0
 ldMC1: - ldx[21,8] + ldw[20,4] >= 0
 ldMC1: - ldx[21,8] + ldw[20,4] - ldw[20,8] <= 0
 ldMC1: lDX[21,8] - lDW[20,4] + lDW[20,8] >= 0
 ldMC1: lDX[21,8] + lDW[20,4] - lDW[20,8] >= 0
 ldMC1: - lDX[21,8] + lDW[20,4] + lDW[20,8] >= 0
 ldMC1: - lDX[21,8] + lDW[20,4] - lDW[20,8] <= 0
 ldMC2: - ldx[21,4] + ldw[20,8] >= 0
 ldMC2: - ldx[21,12] + ldw[20,8] >= 0
 ldMC2: - ldx[21,4] - ldx[21,12] + ldw[20,8] <= 0
 ldMC2: lDX[21,4] + lDX[21,12] - lDW[20,8] >= 0
 ldMC2: - lDX[21,4] + lDX[21,12] + lDW[20,8] >= 0
 ldMC2: lDX[21,4] - lDX[21,12] + lDW[20,8] >= 0
 ldMC2: - lDX[21,4] - lDX[21,12] + lDW[20,8] <= 0
 ldMC3: - ldx[21,0] + ldw[20,12] >= 0
 ldMC3: - ldx[21,12] + ldw[20,12] >= 0
 ldMC3: - ldx[21,0] - ldx[21,12] + ldw[20,12] <= 0
 ldMC3: lDX[21,0] + lDX[21,12] - lDW[20,12] >= 0
 ldMC3: - lDX[21,0] + lDX[21,12] + lDW[20,12] >= 0
 ldMC3: lDX[21,0] - lDX[21,12] + lDW[20,12] >= 0
 ldMC3: - lDX[21,0] - lDX[21,12] + lDW[20,12] <= 0
 ldMC0: - ldx[21,5] + ldw[20,1] = 0
 ldMC0: - lDX[21,5] + lDW[20,1] = 0
 ldMC1: ldw[20,5] - ldw[20,9] >= 0
 ldMC1: - ldx[21,9] + ldw[20,5] >= 0
 ldMC1: - ldx[21,9] + ldw[20,5] - ldw[20,9] <= 0
 ldMC1: lDX[21,9] - lDW[20,5] + lDW[20,9] >= 0
 ldMC1: lDX[21,9] + lDW[20,5] - lDW[20,9] >= 0
 ldMC1: - lDX[21,9] + lDW[20,5] + lDW[20,9] >= 0
 ldMC1: - lDX[21,9] + lDW[20,5] - lDW[20,9] <= 0
 ldMC2: - ldx[21,5] + ldw[20,9] >= 0
 ldMC2: - ldx[21,13] + ldw[20,9] >= 0
 ldMC2: - ldx[21,5] - ldx[21,13] + ldw[20,9] <= 0
 ldMC2: lDX[21,5] + lDX[21,13] - lDW[20,9] >= 0
 ldMC2: - lDX[21,5] + lDX[21,13] + lDW[20,9] >= 0
 ldMC2: lDX[21,5] - lDX[21,13] + lDW[20,9] >= 0
 ldMC2: - lDX[21,5] - lDX[21,13] + lDW[20,9] <= 0
 ldMC3: - ldx[21,1] + ldw[20,13] >= 0
 ldMC3: - ldx[21,13] + ldw[20,13] >= 0
 ldMC3: - ldx[21,1] - ldx[21,13] + ldw[20,13] <= 0
 ldMC3: lDX[21,1] + lDX[21,13] - lDW[20,13] >= 0
 ldMC3: - lDX[21,1] + lDX[21,13] + lDW[20,13] >= 0
 ldMC3: lDX[21,1] - lDX[21,13] + lDW[20,13] >= 0
 ldMC3: - lDX[21,1] - lDX[21,13] + lDW[20,13] <= 0
 ldMC0: - ldx[21,6] + ldw[20,2] = 0
 ldMC0: - lDX[21,6] + lDW[20,2] = 0
 ldMC1: ldw[20,6] - ldw[20,10] >= 0
 ldMC1: - ldx[21,10] + ldw[20,6] >= 0
 ldMC1: - ldx[21,10] + ldw[20,6] - ldw[20,10] <= 0
 ldMC1: lDX[21,10] - lDW[20,6] + lDW[20,10] >= 0
 ldMC1: lDX[21,10] + lDW[20,6] - lDW[20,10] >= 0
 ldMC1: - lDX[21,10] + lDW[20,6] + lDW[20,10] >= 0
 ldMC1: - lDX[21,10] + lDW[20,6] - lDW[20,10] <= 0
 ldMC2: - ldx[21,6] + ldw[20,10] >= 0
 ldMC2: - ldx[21,14] + ldw[20,10] >= 0
 ldMC2: - ldx[21,6] - ldx[21,14] + ldw[20,10] <= 0
 ldMC2: lDX[21,6] + lDX[21,14] - lDW[20,10] >= 0
 ldMC2: - lDX[21,6] + lDX[21,14] + lDW[20,10] >= 0
 ldMC2: lDX[21,6] - lDX[21,14] + lDW[20,10] >= 0
 ldMC2: - lDX[21,6] - lDX[21,14] + lDW[20,10] <= 0
 ldMC3: - ldx[21,2] + ldw[20,14] >= 0
 ldMC3: - ldx[21,14] + ldw[20,14] >= 0
 ldMC3: - ldx[21,2] - ldx[21,14] + ldw[20,14] <= 0
 ldMC3: lDX[21,2] + lDX[21,14] - lDW[20,14] >= 0
 ldMC3: - lDX[21,2] + lDX[21,14] + lDW[20,14] >= 0
 ldMC3: lDX[21,2] - lDX[21,14] + lDW[20,14] >= 0
 ldMC3: - lDX[21,2] - lDX[21,14] + lDW[20,14] <= 0
 ldMC0: - ldx[21,7] + ldw[20,3] = 0
 ldMC0: - lDX[21,7] + lDW[20,3] = 0
 ldMC1: ldw[20,7] - ldw[20,11] >= 0
 ldMC1: - ldx[21,11] + ldw[20,7] >= 0
 ldMC1: - ldx[21,11] + ldw[20,7] - ldw[20,11] <= 0
 ldMC1: lDX[21,11] - lDW[20,7] + lDW[20,11] >= 0
 ldMC1: lDX[21,11] + lDW[20,7] - lDW[20,11] >= 0
 ldMC1: - lDX[21,11] + lDW[20,7] + lDW[20,11] >= 0
 ldMC1: - lDX[21,11] + lDW[20,7] - lDW[20,11] <= 0
 ldMC2: - ldx[21,7] + ldw[20,11] >= 0
 ldMC2: - ldx[21,15] + ldw[20,11] >= 0
 ldMC2: - ldx[21,7] - ldx[21,15] + ldw[20,11] <= 0
 ldMC2: lDX[21,7] + lDX[21,15] - lDW[20,11] >= 0
 ldMC2: - lDX[21,7] + lDX[21,15] + lDW[20,11] >= 0
 ldMC2: lDX[21,7] - lDX[21,15] + lDW[20,11] >= 0
 ldMC2: - lDX[21,7] - lDX[21,15] + lDW[20,11] <= 0
 ldMC3: - ldx[21,3] + ldw[20,15] >= 0
 ldMC3: - ldx[21,15] + ldw[20,15] >= 0
 ldMC3: - ldx[21,3] - ldx[21,15] + ldw[20,15] <= 0
 ldMC3: lDX[21,3] + lDX[21,15] - lDW[20,15] >= 0
 ldMC3: - lDX[21,3] + lDX[21,15] + lDW[20,15] >= 0
 ldMC3: lDX[21,3] - lDX[21,15] + lDW[20,15] >= 0
 ldMC3: - lDX[21,3] - lDX[21,15] + lDW[20,15] <= 0
 ldMC0: - ldx[22,4] + ldw[21,0] = 0
 ldMC0: - lDX[22,4] + lDW[21,0] = 0
 ldMC1: ldw[21,4] - ldw[21,8] >= 0
 ldMC1: - ldx[22,8] + ldw[21,4] >= 0
 ldMC1: - ldx[22,8] + ldw[21,4] - ldw[21,8] <= 0
 ldMC1: lDX[22,8] - lDW[21,4] + lDW[21,8] >= 0
 ldMC1: lDX[22,8] + lDW[21,4] - lDW[21,8] >= 0
 ldMC1: - lDX[22,8] + lDW[21,4] + lDW[21,8] >= 0
 ldMC1: - lDX[22,8] + lDW[21,4] - lDW[21,8] <= 0
 ldMC2: - ldx[22,4] + ldw[21,8] >= 0
 ldMC2: - ldx[22,12] + ldw[21,8] >= 0
 ldMC2: - ldx[22,4] - ldx[22,12] + ldw[21,8] <= 0
 ldMC2: lDX[22,4] + lDX[22,12] - lDW[21,8] >= 0
 ldMC2: - lDX[22,4] + lDX[22,12] + lDW[21,8] >= 0
 ldMC2: lDX[22,4] - lDX[22,12] + lDW[21,8] >= 0
 ldMC2: - lDX[22,4] - lDX[22,12] + lDW[21,8] <= 0
 ldMC3: - ldx[22,0] + ldw[21,12] >= 0
 ldMC3: - ldx[22,12] + ldw[21,12] >= 0
 ldMC3: - ldx[22,0] - ldx[22,12] + ldw[21,12] <= 0
 ldMC3: lDX[22,0] + lDX[22,12] - lDW[21,12] >= 0
 ldMC3: - lDX[22,0] + lDX[22,12] + lDW[21,12] >= 0
 ldMC3: lDX[22,0] - lDX[22,12] + lDW[21,12] >= 0
 ldMC3: - lDX[22,0] - lDX[22,12] + lDW[21,12] <= 0
 ldMC0: - ldx[22,5] + ldw[21,1] = 0
 ldMC0: - lDX[22,5] + lDW[21,1] = 0
 ldMC1: ldw[21,5] - ldw[21,9] >= 0
 ldMC1: - ldx[22,9] + ldw[21,5] >= 0
 ldMC1: - ldx[22,9] + ldw[21,5] - ldw[21,9] <= 0
 ldMC1: lDX[22,9] - lDW[21,5] + lDW[21,9] >= 0
 ldMC1: lDX[22,9] + lDW[21,5] - lDW[21,9] >= 0
 ldMC1: - lDX[22,9] + lDW[21,5] + lDW[21,9] >= 0
 ldMC1: - lDX[22,9] + lDW[21,5] - lDW[21,9] <= 0
 ldMC2: - ldx[22,5] + ldw[21,9] >= 0
 ldMC2: - ldx[22,13] + ldw[21,9] >= 0
 ldMC2: - ldx[22,5] - ldx[22,13] + ldw[21,9] <= 0
 ldMC2: lDX[22,5] + lDX[22,13] - lDW[21,9] >= 0
 ldMC2: - lDX[22,5] + lDX[22,13] + lDW[21,9] >= 0
 ldMC2: lDX[22,5] - lDX[22,13] + lDW[21,9] >= 0
 ldMC2: - lDX[22,5] - lDX[22,13] + lDW[21,9] <= 0
 ldMC3: - ldx[22,1] + ldw[21,13] >= 0
 ldMC3: - ldx[22,13] + ldw[21,13] >= 0
 ldMC3: - ldx[22,1] - ldx[22,13] + ldw[21,13] <= 0
 ldMC3: lDX[22,1] + lDX[22,13] - lDW[21,13] >= 0
 ldMC3: - lDX[22,1] + lDX[22,13] + lDW[21,13] >= 0
 ldMC3: lDX[22,1] - lDX[22,13] + lDW[21,13] >= 0
 ldMC3: - lDX[22,1] - lDX[22,13] + lDW[21,13] <= 0
 ldMC0: - ldx[22,6] + ldw[21,2] = 0
 ldMC0: - lDX[22,6] + lDW[21,2] = 0
 ldMC1: ldw[21,6] - ldw[21,10] >= 0
 ldMC1: - ldx[22,10] + ldw[21,6] >= 0
 ldMC1: - ldx[22,10] + ldw[21,6] - ldw[21,10] <= 0
 ldMC1: lDX[22,10] - lDW[21,6] + lDW[21,10] >= 0
 ldMC1: lDX[22,10] + lDW[21,6] - lDW[21,10] >= 0
 ldMC1: - lDX[22,10] + lDW[21,6] + lDW[21,10] >= 0
 ldMC1: - lDX[22,10] + lDW[21,6] - lDW[21,10] <= 0
 ldMC2: - ldx[22,6] + ldw[21,10] >= 0
 ldMC2: - ldx[22,14] + ldw[21,10] >= 0
 ldMC2: - ldx[22,6] - ldx[22,14] + ldw[21,10] <= 0
 ldMC2: lDX[22,6] + lDX[22,14] - lDW[21,10] >= 0
 ldMC2: - lDX[22,6] + lDX[22,14] + lDW[21,10] >= 0
 ldMC2: lDX[22,6] - lDX[22,14] + lDW[21,10] >= 0
 ldMC2: - lDX[22,6] - lDX[22,14] + lDW[21,10] <= 0
 ldMC3: - ldx[22,2] + ldw[21,14] >= 0
 ldMC3: - ldx[22,14] + ldw[21,14] >= 0
 ldMC3: - ldx[22,2] - ldx[22,14] + ldw[21,14] <= 0
 ldMC3: lDX[22,2] + lDX[22,14] - lDW[21,14] >= 0
 ldMC3: - lDX[22,2] + lDX[22,14] + lDW[21,14] >= 0
 ldMC3: lDX[22,2] - lDX[22,14] + lDW[21,14] >= 0
 ldMC3: - lDX[22,2] - lDX[22,14] + lDW[21,14] <= 0
 ldMC0: - ldx[22,7] + ldw[21,3] = 0
 ldMC0: - lDX[22,7] + lDW[21,3] = 0
 ldMC1: ldw[21,7] - ldw[21,11] >= 0
 ldMC1: - ldx[22,11] + ldw[21,7] >= 0
 ldMC1: - ldx[22,11] + ldw[21,7] - ldw[21,11] <= 0
 ldMC1: lDX[22,11] - lDW[21,7] + lDW[21,11] >= 0
 ldMC1: lDX[22,11] + lDW[21,7] - lDW[21,11] >= 0
 ldMC1: - lDX[22,11] + lDW[21,7] + lDW[21,11] >= 0
 ldMC1: - lDX[22,11] + lDW[21,7] - lDW[21,11] <= 0
 ldMC2: - ldx[22,7] + ldw[21,11] >= 0
 ldMC2: - ldx[22,15] + ldw[21,11] >= 0
 ldMC2: - ldx[22,7] - ldx[22,15] + ldw[21,11] <= 0
 ldMC2: lDX[22,7] + lDX[22,15] - lDW[21,11] >= 0
 ldMC2: - lDX[22,7] + lDX[22,15] + lDW[21,11] >= 0
 ldMC2: lDX[22,7] - lDX[22,15] + lDW[21,11] >= 0
 ldMC2: - lDX[22,7] - lDX[22,15] + lDW[21,11] <= 0
 ldMC3: - ldx[22,3] + ldw[21,15] >= 0
 ldMC3: - ldx[22,15] + ldw[21,15] >= 0
 ldMC3: - ldx[22,3] - ldx[22,15] + ldw[21,15] <= 0
 ldMC3: lDX[22,3] + lDX[22,15] - lDW[21,15] >= 0
 ldMC3: - lDX[22,3] + lDX[22,15] + lDW[21,15] >= 0
 ldMC3: lDX[22,3] - lDX[22,15] + lDW[21,15] >= 0
 ldMC3: - lDX[22,3] - lDX[22,15] + lDW[21,15] <= 0
 ldMC0: - ldx[23,4] + ldw[22,0] = 0
 ldMC0: - lDX[23,4] + lDW[22,0] = 0
 ldMC1: ldw[22,4] - ldw[22,8] >= 0
 ldMC1: - ldx[23,8] + ldw[22,4] >= 0
 ldMC1: - ldx[23,8] + ldw[22,4] - ldw[22,8] <= 0
 ldMC1: lDX[23,8] - lDW[22,4] + lDW[22,8] >= 0
 ldMC1: lDX[23,8] + lDW[22,4] - lDW[22,8] >= 0
 ldMC1: - lDX[23,8] + lDW[22,4] + lDW[22,8] >= 0
 ldMC1: - lDX[23,8] + lDW[22,4] - lDW[22,8] <= 0
 ldMC2: - ldx[23,4] + ldw[22,8] >= 0
 ldMC2: - ldx[23,12] + ldw[22,8] >= 0
 ldMC2: - ldx[23,4] - ldx[23,12] + ldw[22,8] <= 0
 ldMC2: lDX[23,4] + lDX[23,12] - lDW[22,8] >= 0
 ldMC2: - lDX[23,4] + lDX[23,12] + lDW[22,8] >= 0
 ldMC2: lDX[23,4] - lDX[23,12] + lDW[22,8] >= 0
 ldMC2: - lDX[23,4] - lDX[23,12] + lDW[22,8] <= 0
 ldMC3: - ldx[23,0] + ldw[22,12] >= 0
 ldMC3: - ldx[23,12] + ldw[22,12] >= 0
 ldMC3: - ldx[23,0] - ldx[23,12] + ldw[22,12] <= 0
 ldMC3: lDX[23,0] + lDX[23,12] - lDW[22,12] >= 0
 ldMC3: - lDX[23,0] + lDX[23,12] + lDW[22,12] >= 0
 ldMC3: lDX[23,0] - lDX[23,12] + lDW[22,12] >= 0
 ldMC3: - lDX[23,0] - lDX[23,12] + lDW[22,12] <= 0
 ldMC0: - ldx[23,5] + ldw[22,1] = 0
 ldMC0: - lDX[23,5] + lDW[22,1] = 0
 ldMC1: ldw[22,5] - ldw[22,9] >= 0
 ldMC1: - ldx[23,9] + ldw[22,5] >= 0
 ldMC1: - ldx[23,9] + ldw[22,5] - ldw[22,9] <= 0
 ldMC1: lDX[23,9] - lDW[22,5] + lDW[22,9] >= 0
 ldMC1: lDX[23,9] + lDW[22,5] - lDW[22,9] >= 0
 ldMC1: - lDX[23,9] + lDW[22,5] + lDW[22,9] >= 0
 ldMC1: - lDX[23,9] + lDW[22,5] - lDW[22,9] <= 0
 ldMC2: - ldx[23,5] + ldw[22,9] >= 0
 ldMC2: - ldx[23,13] + ldw[22,9] >= 0
 ldMC2: - ldx[23,5] - ldx[23,13] + ldw[22,9] <= 0
 ldMC2: lDX[23,5] + lDX[23,13] - lDW[22,9] >= 0
 ldMC2: - lDX[23,5] + lDX[23,13] + lDW[22,9] >= 0
 ldMC2: lDX[23,5] - lDX[23,13] + lDW[22,9] >= 0
 ldMC2: - lDX[23,5] - lDX[23,13] + lDW[22,9] <= 0
 ldMC3: - ldx[23,1] + ldw[22,13] >= 0
 ldMC3: - ldx[23,13] + ldw[22,13] >= 0
 ldMC3: - ldx[23,1] - ldx[23,13] + ldw[22,13] <= 0
 ldMC3: lDX[23,1] + lDX[23,13] - lDW[22,13] >= 0
 ldMC3: - lDX[23,1] + lDX[23,13] + lDW[22,13] >= 0
 ldMC3: lDX[23,1] - lDX[23,13] + lDW[22,13] >= 0
 ldMC3: - lDX[23,1] - lDX[23,13] + lDW[22,13] <= 0
 ldMC0: - ldx[23,6] + ldw[22,2] = 0
 ldMC0: - lDX[23,6] + lDW[22,2] = 0
 ldMC1: ldw[22,6] - ldw[22,10] >= 0
 ldMC1: - ldx[23,10] + ldw[22,6] >= 0
 ldMC1: - ldx[23,10] + ldw[22,6] - ldw[22,10] <= 0
 ldMC1: lDX[23,10] - lDW[22,6] + lDW[22,10] >= 0
 ldMC1: lDX[23,10] + lDW[22,6] - lDW[22,10] >= 0
 ldMC1: - lDX[23,10] + lDW[22,6] + lDW[22,10] >= 0
 ldMC1: - lDX[23,10] + lDW[22,6] - lDW[22,10] <= 0
 ldMC2: - ldx[23,6] + ldw[22,10] >= 0
 ldMC2: - ldx[23,14] + ldw[22,10] >= 0
 ldMC2: - ldx[23,6] - ldx[23,14] + ldw[22,10] <= 0
 ldMC2: lDX[23,6] + lDX[23,14] - lDW[22,10] >= 0
 ldMC2: - lDX[23,6] + lDX[23,14] + lDW[22,10] >= 0
 ldMC2: lDX[23,6] - lDX[23,14] + lDW[22,10] >= 0
 ldMC2: - lDX[23,6] - lDX[23,14] + lDW[22,10] <= 0
 ldMC3: - ldx[23,2] + ldw[22,14] >= 0
 ldMC3: - ldx[23,14] + ldw[22,14] >= 0
 ldMC3: - ldx[23,2] - ldx[23,14] + ldw[22,14] <= 0
 ldMC3: lDX[23,2] + lDX[23,14] - lDW[22,14] >= 0
 ldMC3: - lDX[23,2] + lDX[23,14] + lDW[22,14] >= 0
 ldMC3: lDX[23,2] - lDX[23,14] + lDW[22,14] >= 0
 ldMC3: - lDX[23,2] - lDX[23,14] + lDW[22,14] <= 0
 ldMC0: - ldx[23,7] + ldw[22,3] = 0
 ldMC0: - lDX[23,7] + lDW[22,3] = 0
 ldMC1: ldw[22,7] - ldw[22,11] >= 0
 ldMC1: - ldx[23,11] + ldw[22,7] >= 0
 ldMC1: - ldx[23,11] + ldw[22,7] - ldw[22,11] <= 0
 ldMC1: lDX[23,11] - lDW[22,7] + lDW[22,11] >= 0
 ldMC1: lDX[23,11] + lDW[22,7] - lDW[22,11] >= 0
 ldMC1: - lDX[23,11] + lDW[22,7] + lDW[22,11] >= 0
 ldMC1: - lDX[23,11] + lDW[22,7] - lDW[22,11] <= 0
 ldMC2: - ldx[23,7] + ldw[22,11] >= 0
 ldMC2: - ldx[23,15] + ldw[22,11] >= 0
 ldMC2: - ldx[23,7] - ldx[23,15] + ldw[22,11] <= 0
 ldMC2: lDX[23,7] + lDX[23,15] - lDW[22,11] >= 0
 ldMC2: - lDX[23,7] + lDX[23,15] + lDW[22,11] >= 0
 ldMC2: lDX[23,7] - lDX[23,15] + lDW[22,11] >= 0
 ldMC2: - lDX[23,7] - lDX[23,15] + lDW[22,11] <= 0
 ldMC3: - ldx[23,3] + ldw[22,15] >= 0
 ldMC3: - ldx[23,15] + ldw[22,15] >= 0
 ldMC3: - ldx[23,3] - ldx[23,15] + ldw[22,15] <= 0
 ldMC3: lDX[23,3] + lDX[23,15] - lDW[22,15] >= 0
 ldMC3: - lDX[23,3] + lDX[23,15] + lDW[22,15] >= 0
 ldMC3: lDX[23,3] - lDX[23,15] + lDW[22,15] >= 0
 ldMC3: - lDX[23,3] - lDX[23,15] + lDW[22,15] <= 0
 ldMC0: - ldx[24,4] + ldw[23,0] = 0
 ldMC0: - lDX[24,4] + lDW[23,0] = 0
 ldMC1: ldw[23,4] - ldw[23,8] >= 0
 ldMC1: - ldx[24,8] + ldw[23,4] >= 0
 ldMC1: - ldx[24,8] + ldw[23,4] - ldw[23,8] <= 0
 ldMC1: lDX[24,8] - lDW[23,4] + lDW[23,8] >= 0
 ldMC1: lDX[24,8] + lDW[23,4] - lDW[23,8] >= 0
 ldMC1: - lDX[24,8] + lDW[23,4] + lDW[23,8] >= 0
 ldMC1: - lDX[24,8] + lDW[23,4] - lDW[23,8] <= 0
 ldMC2: - ldx[24,4] + ldw[23,8] >= 0
 ldMC2: - ldx[24,12] + ldw[23,8] >= 0
 ldMC2: - ldx[24,4] - ldx[24,12] + ldw[23,8] <= 0
 ldMC2: lDX[24,4] + lDX[24,12] - lDW[23,8] >= 0
 ldMC2: - lDX[24,4] + lDX[24,12] + lDW[23,8] >= 0
 ldMC2: lDX[24,4] - lDX[24,12] + lDW[23,8] >= 0
 ldMC2: - lDX[24,4] - lDX[24,12] + lDW[23,8] <= 0
 ldMC3: - ldx[24,0] + ldw[23,12] >= 0
 ldMC3: - ldx[24,12] + ldw[23,12] >= 0
 ldMC3: - ldx[24,0] - ldx[24,12] + ldw[23,12] <= 0
 ldMC3: lDX[24,0] + lDX[24,12] - lDW[23,12] >= 0
 ldMC3: - lDX[24,0] + lDX[24,12] + lDW[23,12] >= 0
 ldMC3: lDX[24,0] - lDX[24,12] + lDW[23,12] >= 0
 ldMC3: - lDX[24,0] - lDX[24,12] + lDW[23,12] <= 0
 ldMC0: - ldx[24,5] + ldw[23,1] = 0
 ldMC0: - lDX[24,5] + lDW[23,1] = 0
 ldMC1: ldw[23,5] - ldw[23,9] >= 0
 ldMC1: - ldx[24,9] + ldw[23,5] >= 0
 ldMC1: - ldx[24,9] + ldw[23,5] - ldw[23,9] <= 0
 ldMC1: lDX[24,9] - lDW[23,5] + lDW[23,9] >= 0
 ldMC1: lDX[24,9] + lDW[23,5] - lDW[23,9] >= 0
 ldMC1: - lDX[24,9] + lDW[23,5] + lDW[23,9] >= 0
 ldMC1: - lDX[24,9] + lDW[23,5] - lDW[23,9] <= 0
 ldMC2: - ldx[24,5] + ldw[23,9] >= 0
 ldMC2: - ldx[24,13] + ldw[23,9] >= 0
 ldMC2: - ldx[24,5] - ldx[24,13] + ldw[23,9] <= 0
 ldMC2: lDX[24,5] + lDX[24,13] - lDW[23,9] >= 0
 ldMC2: - lDX[24,5] + lDX[24,13] + lDW[23,9] >= 0
 ldMC2: lDX[24,5] - lDX[24,13] + lDW[23,9] >= 0
 ldMC2: - lDX[24,5] - lDX[24,13] + lDW[23,9] <= 0
 ldMC3: - ldx[24,1] + ldw[23,13] >= 0
 ldMC3: - ldx[24,13] + ldw[23,13] >= 0
 ldMC3: - ldx[24,1] - ldx[24,13] + ldw[23,13] <= 0
 ldMC3: lDX[24,1] + lDX[24,13] - lDW[23,13] >= 0
 ldMC3: - lDX[24,1] + lDX[24,13] + lDW[23,13] >= 0
 ldMC3: lDX[24,1] - lDX[24,13] + lDW[23,13] >= 0
 ldMC3: - lDX[24,1] - lDX[24,13] + lDW[23,13] <= 0
 ldMC0: - ldx[24,6] + ldw[23,2] = 0
 ldMC0: - lDX[24,6] + lDW[23,2] = 0
 ldMC1: ldw[23,6] - ldw[23,10] >= 0
 ldMC1: - ldx[24,10] + ldw[23,6] >= 0
 ldMC1: - ldx[24,10] + ldw[23,6] - ldw[23,10] <= 0
 ldMC1: lDX[24,10] - lDW[23,6] + lDW[23,10] >= 0
 ldMC1: lDX[24,10] + lDW[23,6] - lDW[23,10] >= 0
 ldMC1: - lDX[24,10] + lDW[23,6] + lDW[23,10] >= 0
 ldMC1: - lDX[24,10] + lDW[23,6] - lDW[23,10] <= 0
 ldMC2: - ldx[24,6] + ldw[23,10] >= 0
 ldMC2: - ldx[24,14] + ldw[23,10] >= 0
 ldMC2: - ldx[24,6] - ldx[24,14] + ldw[23,10] <= 0
 ldMC2: lDX[24,6] + lDX[24,14] - lDW[23,10] >= 0
 ldMC2: - lDX[24,6] + lDX[24,14] + lDW[23,10] >= 0
 ldMC2: lDX[24,6] - lDX[24,14] + lDW[23,10] >= 0
 ldMC2: - lDX[24,6] - lDX[24,14] + lDW[23,10] <= 0
 ldMC3: - ldx[24,2] + ldw[23,14] >= 0
 ldMC3: - ldx[24,14] + ldw[23,14] >= 0
 ldMC3: - ldx[24,2] - ldx[24,14] + ldw[23,14] <= 0
 ldMC3: lDX[24,2] + lDX[24,14] - lDW[23,14] >= 0
 ldMC3: - lDX[24,2] + lDX[24,14] + lDW[23,14] >= 0
 ldMC3: lDX[24,2] - lDX[24,14] + lDW[23,14] >= 0
 ldMC3: - lDX[24,2] - lDX[24,14] + lDW[23,14] <= 0
 ldMC0: - ldx[24,7] + ldw[23,3] = 0
 ldMC0: - lDX[24,7] + lDW[23,3] = 0
 ldMC1: ldw[23,7] - ldw[23,11] >= 0
 ldMC1: - ldx[24,11] + ldw[23,7] >= 0
 ldMC1: - ldx[24,11] + ldw[23,7] - ldw[23,11] <= 0
 ldMC1: lDX[24,11] - lDW[23,7] + lDW[23,11] >= 0
 ldMC1: lDX[24,11] + lDW[23,7] - lDW[23,11] >= 0
 ldMC1: - lDX[24,11] + lDW[23,7] + lDW[23,11] >= 0
 ldMC1: - lDX[24,11] + lDW[23,7] - lDW[23,11] <= 0
 ldMC2: - ldx[24,7] + ldw[23,11] >= 0
 ldMC2: - ldx[24,15] + ldw[23,11] >= 0
 ldMC2: - ldx[24,7] - ldx[24,15] + ldw[23,11] <= 0
 ldMC2: lDX[24,7] + lDX[24,15] - lDW[23,11] >= 0
 ldMC2: - lDX[24,7] + lDX[24,15] + lDW[23,11] >= 0
 ldMC2: lDX[24,7] - lDX[24,15] + lDW[23,11] >= 0
 ldMC2: - lDX[24,7] - lDX[24,15] + lDW[23,11] <= 0
 ldMC3: - ldx[24,3] + ldw[23,15] >= 0
 ldMC3: - ldx[24,15] + ldw[23,15] >= 0
 ldMC3: - ldx[24,3] - ldx[24,15] + ldw[23,15] <= 0
 ldMC3: lDX[24,3] + lDX[24,15] - lDW[23,15] >= 0
 ldMC3: - lDX[24,3] + lDX[24,15] + lDW[23,15] >= 0
 ldMC3: lDX[24,3] - lDX[24,15] + lDW[23,15] >= 0
 ldMC3: - lDX[24,3] - lDX[24,15] + lDW[23,15] <= 0
 ldMC0: - ldx[25,4] + ldw[24,0] = 0
 ldMC0: - lDX[25,4] + lDW[24,0] = 0
 ldMC1: ldw[24,4] - ldw[24,8] >= 0
 ldMC1: - ldx[25,8] + ldw[24,4] >= 0
 ldMC1: - ldx[25,8] + ldw[24,4] - ldw[24,8] <= 0
 ldMC1: lDX[25,8] - lDW[24,4] + lDW[24,8] >= 0
 ldMC1: lDX[25,8] + lDW[24,4] - lDW[24,8] >= 0
 ldMC1: - lDX[25,8] + lDW[24,4] + lDW[24,8] >= 0
 ldMC1: - lDX[25,8] + lDW[24,4] - lDW[24,8] <= 0
 ldMC2: - ldx[25,4] + ldw[24,8] >= 0
 ldMC2: - ldx[25,12] + ldw[24,8] >= 0
 ldMC2: - ldx[25,4] - ldx[25,12] + ldw[24,8] <= 0
 ldMC2: lDX[25,4] + lDX[25,12] - lDW[24,8] >= 0
 ldMC2: - lDX[25,4] + lDX[25,12] + lDW[24,8] >= 0
 ldMC2: lDX[25,4] - lDX[25,12] + lDW[24,8] >= 0
 ldMC2: - lDX[25,4] - lDX[25,12] + lDW[24,8] <= 0
 ldMC3: - ldx[25,0] + ldw[24,12] >= 0
 ldMC3: - ldx[25,12] + ldw[24,12] >= 0
 ldMC3: - ldx[25,0] - ldx[25,12] + ldw[24,12] <= 0
 ldMC3: lDX[25,0] + lDX[25,12] - lDW[24,12] >= 0
 ldMC3: - lDX[25,0] + lDX[25,12] + lDW[24,12] >= 0
 ldMC3: lDX[25,0] - lDX[25,12] + lDW[24,12] >= 0
 ldMC3: - lDX[25,0] - lDX[25,12] + lDW[24,12] <= 0
 ldMC0: - ldx[25,5] + ldw[24,1] = 0
 ldMC0: - lDX[25,5] + lDW[24,1] = 0
 ldMC1: ldw[24,5] - ldw[24,9] >= 0
 ldMC1: - ldx[25,9] + ldw[24,5] >= 0
 ldMC1: - ldx[25,9] + ldw[24,5] - ldw[24,9] <= 0
 ldMC1: lDX[25,9] - lDW[24,5] + lDW[24,9] >= 0
 ldMC1: lDX[25,9] + lDW[24,5] - lDW[24,9] >= 0
 ldMC1: - lDX[25,9] + lDW[24,5] + lDW[24,9] >= 0
 ldMC1: - lDX[25,9] + lDW[24,5] - lDW[24,9] <= 0
 ldMC2: - ldx[25,5] + ldw[24,9] >= 0
 ldMC2: - ldx[25,13] + ldw[24,9] >= 0
 ldMC2: - ldx[25,5] - ldx[25,13] + ldw[24,9] <= 0
 ldMC2: lDX[25,5] + lDX[25,13] - lDW[24,9] >= 0
 ldMC2: - lDX[25,5] + lDX[25,13] + lDW[24,9] >= 0
 ldMC2: lDX[25,5] - lDX[25,13] + lDW[24,9] >= 0
 ldMC2: - lDX[25,5] - lDX[25,13] + lDW[24,9] <= 0
 ldMC3: - ldx[25,1] + ldw[24,13] >= 0
 ldMC3: - ldx[25,13] + ldw[24,13] >= 0
 ldMC3: - ldx[25,1] - ldx[25,13] + ldw[24,13] <= 0
 ldMC3: lDX[25,1] + lDX[25,13] - lDW[24,13] >= 0
 ldMC3: - lDX[25,1] + lDX[25,13] + lDW[24,13] >= 0
 ldMC3: lDX[25,1] - lDX[25,13] + lDW[24,13] >= 0
 ldMC3: - lDX[25,1] - lDX[25,13] + lDW[24,13] <= 0
 ldMC0: - ldx[25,6] + ldw[24,2] = 0
 ldMC0: - lDX[25,6] + lDW[24,2] = 0
 ldMC1: ldw[24,6] - ldw[24,10] >= 0
 ldMC1: - ldx[25,10] + ldw[24,6] >= 0
 ldMC1: - ldx[25,10] + ldw[24,6] - ldw[24,10] <= 0
 ldMC1: lDX[25,10] - lDW[24,6] + lDW[24,10] >= 0
 ldMC1: lDX[25,10] + lDW[24,6] - lDW[24,10] >= 0
 ldMC1: - lDX[25,10] + lDW[24,6] + lDW[24,10] >= 0
 ldMC1: - lDX[25,10] + lDW[24,6] - lDW[24,10] <= 0
 ldMC2: - ldx[25,6] + ldw[24,10] >= 0
 ldMC2: - ldx[25,14] + ldw[24,10] >= 0
 ldMC2: - ldx[25,6] - ldx[25,14] + ldw[24,10] <= 0
 ldMC2: lDX[25,6] + lDX[25,14] - lDW[24,10] >= 0
 ldMC2: - lDX[25,6] + lDX[25,14] + lDW[24,10] >= 0
 ldMC2: lDX[25,6] - lDX[25,14] + lDW[24,10] >= 0
 ldMC2: - lDX[25,6] - lDX[25,14] + lDW[24,10] <= 0
 ldMC3: - ldx[25,2] + ldw[24,14] >= 0
 ldMC3: - ldx[25,14] + ldw[24,14] >= 0
 ldMC3: - ldx[25,2] - ldx[25,14] + ldw[24,14] <= 0
 ldMC3: lDX[25,2] + lDX[25,14] - lDW[24,14] >= 0
 ldMC3: - lDX[25,2] + lDX[25,14] + lDW[24,14] >= 0
 ldMC3: lDX[25,2] - lDX[25,14] + lDW[24,14] >= 0
 ldMC3: - lDX[25,2] - lDX[25,14] + lDW[24,14] <= 0
 ldMC0: - ldx[25,7] + ldw[24,3] = 0
 ldMC0: - lDX[25,7] + lDW[24,3] = 0
 ldMC1: ldw[24,7] - ldw[24,11] >= 0
 ldMC1: - ldx[25,11] + ldw[24,7] >= 0
 ldMC1: - ldx[25,11] + ldw[24,7] - ldw[24,11] <= 0
 ldMC1: lDX[25,11] - lDW[24,7] + lDW[24,11] >= 0
 ldMC1: lDX[25,11] + lDW[24,7] - lDW[24,11] >= 0
 ldMC1: - lDX[25,11] + lDW[24,7] + lDW[24,11] >= 0
 ldMC1: - lDX[25,11] + lDW[24,7] - lDW[24,11] <= 0
 ldMC2: - ldx[25,7] + ldw[24,11] >= 0
 ldMC2: - ldx[25,15] + ldw[24,11] >= 0
 ldMC2: - ldx[25,7] - ldx[25,15] + ldw[24,11] <= 0
 ldMC2: lDX[25,7] + lDX[25,15] - lDW[24,11] >= 0
 ldMC2: - lDX[25,7] + lDX[25,15] + lDW[24,11] >= 0
 ldMC2: lDX[25,7] - lDX[25,15] + lDW[24,11] >= 0
 ldMC2: - lDX[25,7] - lDX[25,15] + lDW[24,11] <= 0
 ldMC3: - ldx[25,3] + ldw[24,15] >= 0
 ldMC3: - ldx[25,15] + ldw[24,15] >= 0
 ldMC3: - ldx[25,3] - ldx[25,15] + ldw[24,15] <= 0
 ldMC3: lDX[25,3] + lDX[25,15] - lDW[24,15] >= 0
 ldMC3: - lDX[25,3] + lDX[25,15] + lDW[24,15] >= 0
 ldMC3: lDX[25,3] - lDX[25,15] + lDW[24,15] >= 0
 ldMC3: - lDX[25,3] - lDX[25,15] + lDW[24,15] <= 0
 ldMC0: - ldx[26,4] + ldw[25,0] = 0
 ldMC0: - lDX[26,4] + lDW[25,0] = 0
 ldMC1: ldw[25,4] - ldw[25,8] >= 0
 ldMC1: - ldx[26,8] + ldw[25,4] >= 0
 ldMC1: - ldx[26,8] + ldw[25,4] - ldw[25,8] <= 0
 ldMC1: lDX[26,8] - lDW[25,4] + lDW[25,8] >= 0
 ldMC1: lDX[26,8] + lDW[25,4] - lDW[25,8] >= 0
 ldMC1: - lDX[26,8] + lDW[25,4] + lDW[25,8] >= 0
 ldMC1: - lDX[26,8] + lDW[25,4] - lDW[25,8] <= 0
 ldMC2: - ldx[26,4] + ldw[25,8] >= 0
 ldMC2: - ldx[26,12] + ldw[25,8] >= 0
 ldMC2: - ldx[26,4] - ldx[26,12] + ldw[25,8] <= 0
 ldMC2: lDX[26,4] + lDX[26,12] - lDW[25,8] >= 0
 ldMC2: - lDX[26,4] + lDX[26,12] + lDW[25,8] >= 0
 ldMC2: lDX[26,4] - lDX[26,12] + lDW[25,8] >= 0
 ldMC2: - lDX[26,4] - lDX[26,12] + lDW[25,8] <= 0
 ldMC3: - ldx[26,0] + ldw[25,12] >= 0
 ldMC3: - ldx[26,12] + ldw[25,12] >= 0
 ldMC3: - ldx[26,0] - ldx[26,12] + ldw[25,12] <= 0
 ldMC3: lDX[26,0] + lDX[26,12] - lDW[25,12] >= 0
 ldMC3: - lDX[26,0] + lDX[26,12] + lDW[25,12] >= 0
 ldMC3: lDX[26,0] - lDX[26,12] + lDW[25,12] >= 0
 ldMC3: - lDX[26,0] - lDX[26,12] + lDW[25,12] <= 0
 ldMC0: - ldx[26,5] + ldw[25,1] = 0
 ldMC0: - lDX[26,5] + lDW[25,1] = 0
 ldMC1: ldw[25,5] - ldw[25,9] >= 0
 ldMC1: - ldx[26,9] + ldw[25,5] >= 0
 ldMC1: - ldx[26,9] + ldw[25,5] - ldw[25,9] <= 0
 ldMC1: lDX[26,9] - lDW[25,5] + lDW[25,9] >= 0
 ldMC1: lDX[26,9] + lDW[25,5] - lDW[25,9] >= 0
 ldMC1: - lDX[26,9] + lDW[25,5] + lDW[25,9] >= 0
 ldMC1: - lDX[26,9] + lDW[25,5] - lDW[25,9] <= 0
 ldMC2: - ldx[26,5] + ldw[25,9] >= 0
 ldMC2: - ldx[26,13] + ldw[25,9] >= 0
 ldMC2: - ldx[26,5] - ldx[26,13] + ldw[25,9] <= 0
 ldMC2: lDX[26,5] + lDX[26,13] - lDW[25,9] >= 0
 ldMC2: - lDX[26,5] + lDX[26,13] + lDW[25,9] >= 0
 ldMC2: lDX[26,5] - lDX[26,13] + lDW[25,9] >= 0
 ldMC2: - lDX[26,5] - lDX[26,13] + lDW[25,9] <= 0
 ldMC3: - ldx[26,1] + ldw[25,13] >= 0
 ldMC3: - ldx[26,13] + ldw[25,13] >= 0
 ldMC3: - ldx[26,1] - ldx[26,13] + ldw[25,13] <= 0
 ldMC3: lDX[26,1] + lDX[26,13] - lDW[25,13] >= 0
 ldMC3: - lDX[26,1] + lDX[26,13] + lDW[25,13] >= 0
 ldMC3: lDX[26,1] - lDX[26,13] + lDW[25,13] >= 0
 ldMC3: - lDX[26,1] - lDX[26,13] + lDW[25,13] <= 0
 ldMC0: - ldx[26,6] + ldw[25,2] = 0
 ldMC0: - lDX[26,6] + lDW[25,2] = 0
 ldMC1: ldw[25,6] - ldw[25,10] >= 0
 ldMC1: - ldx[26,10] + ldw[25,6] >= 0
 ldMC1: - ldx[26,10] + ldw[25,6] - ldw[25,10] <= 0
 ldMC1: lDX[26,10] - lDW[25,6] + lDW[25,10] >= 0
 ldMC1: lDX[26,10] + lDW[25,6] - lDW[25,10] >= 0
 ldMC1: - lDX[26,10] + lDW[25,6] + lDW[25,10] >= 0
 ldMC1: - lDX[26,10] + lDW[25,6] - lDW[25,10] <= 0
 ldMC2: - ldx[26,6] + ldw[25,10] >= 0
 ldMC2: - ldx[26,14] + ldw[25,10] >= 0
 ldMC2: - ldx[26,6] - ldx[26,14] + ldw[25,10] <= 0
 ldMC2: lDX[26,6] + lDX[26,14] - lDW[25,10] >= 0
 ldMC2: - lDX[26,6] + lDX[26,14] + lDW[25,10] >= 0
 ldMC2: lDX[26,6] - lDX[26,14] + lDW[25,10] >= 0
 ldMC2: - lDX[26,6] - lDX[26,14] + lDW[25,10] <= 0
 ldMC3: - ldx[26,2] + ldw[25,14] >= 0
 ldMC3: - ldx[26,14] + ldw[25,14] >= 0
 ldMC3: - ldx[26,2] - ldx[26,14] + ldw[25,14] <= 0
 ldMC3: lDX[26,2] + lDX[26,14] - lDW[25,14] >= 0
 ldMC3: - lDX[26,2] + lDX[26,14] + lDW[25,14] >= 0
 ldMC3: lDX[26,2] - lDX[26,14] + lDW[25,14] >= 0
 ldMC3: - lDX[26,2] - lDX[26,14] + lDW[25,14] <= 0
 ldMC0: - ldx[26,7] + ldw[25,3] = 0
 ldMC0: - lDX[26,7] + lDW[25,3] = 0
 ldMC1: ldw[25,7] - ldw[25,11] >= 0
 ldMC1: - ldx[26,11] + ldw[25,7] >= 0
 ldMC1: - ldx[26,11] + ldw[25,7] - ldw[25,11] <= 0
 ldMC1: lDX[26,11] - lDW[25,7] + lDW[25,11] >= 0
 ldMC1: lDX[26,11] + lDW[25,7] - lDW[25,11] >= 0
 ldMC1: - lDX[26,11] + lDW[25,7] + lDW[25,11] >= 0
 ldMC1: - lDX[26,11] + lDW[25,7] - lDW[25,11] <= 0
 ldMC2: - ldx[26,7] + ldw[25,11] >= 0
 ldMC2: - ldx[26,15] + ldw[25,11] >= 0
 ldMC2: - ldx[26,7] - ldx[26,15] + ldw[25,11] <= 0
 ldMC2: lDX[26,7] + lDX[26,15] - lDW[25,11] >= 0
 ldMC2: - lDX[26,7] + lDX[26,15] + lDW[25,11] >= 0
 ldMC2: lDX[26,7] - lDX[26,15] + lDW[25,11] >= 0
 ldMC2: - lDX[26,7] - lDX[26,15] + lDW[25,11] <= 0
 ldMC3: - ldx[26,3] + ldw[25,15] >= 0
 ldMC3: - ldx[26,15] + ldw[25,15] >= 0
 ldMC3: - ldx[26,3] - ldx[26,15] + ldw[25,15] <= 0
 ldMC3: lDX[26,3] + lDX[26,15] - lDW[25,15] >= 0
 ldMC3: - lDX[26,3] + lDX[26,15] + lDW[25,15] >= 0
 ldMC3: lDX[26,3] - lDX[26,15] + lDW[25,15] >= 0
 ldMC3: - lDX[26,3] - lDX[26,15] + lDW[25,15] <= 0
 ldMC0: - ldx[27,4] + ldw[26,0] = 0
 ldMC0: - lDX[27,4] + lDW[26,0] = 0
 ldMC1: ldw[26,4] - ldw[26,8] >= 0
 ldMC1: - ldx[27,8] + ldw[26,4] >= 0
 ldMC1: - ldx[27,8] + ldw[26,4] - ldw[26,8] <= 0
 ldMC1: lDX[27,8] - lDW[26,4] + lDW[26,8] >= 0
 ldMC1: lDX[27,8] + lDW[26,4] - lDW[26,8] >= 0
 ldMC1: - lDX[27,8] + lDW[26,4] + lDW[26,8] >= 0
 ldMC1: - lDX[27,8] + lDW[26,4] - lDW[26,8] <= 0
 ldMC2: - ldx[27,4] + ldw[26,8] >= 0
 ldMC2: - ldx[27,12] + ldw[26,8] >= 0
 ldMC2: - ldx[27,4] - ldx[27,12] + ldw[26,8] <= 0
 ldMC2: lDX[27,4] + lDX[27,12] - lDW[26,8] >= 0
 ldMC2: - lDX[27,4] + lDX[27,12] + lDW[26,8] >= 0
 ldMC2: lDX[27,4] - lDX[27,12] + lDW[26,8] >= 0
 ldMC2: - lDX[27,4] - lDX[27,12] + lDW[26,8] <= 0
 ldMC3: - ldx[27,0] + ldw[26,12] >= 0
 ldMC3: - ldx[27,12] + ldw[26,12] >= 0
 ldMC3: - ldx[27,0] - ldx[27,12] + ldw[26,12] <= 0
 ldMC3: lDX[27,0] + lDX[27,12] - lDW[26,12] >= 0
 ldMC3: - lDX[27,0] + lDX[27,12] + lDW[26,12] >= 0
 ldMC3: lDX[27,0] - lDX[27,12] + lDW[26,12] >= 0
 ldMC3: - lDX[27,0] - lDX[27,12] + lDW[26,12] <= 0
 ldMC0: - ldx[27,5] + ldw[26,1] = 0
 ldMC0: - lDX[27,5] + lDW[26,1] = 0
 ldMC1: ldw[26,5] - ldw[26,9] >= 0
 ldMC1: - ldx[27,9] + ldw[26,5] >= 0
 ldMC1: - ldx[27,9] + ldw[26,5] - ldw[26,9] <= 0
 ldMC1: lDX[27,9] - lDW[26,5] + lDW[26,9] >= 0
 ldMC1: lDX[27,9] + lDW[26,5] - lDW[26,9] >= 0
 ldMC1: - lDX[27,9] + lDW[26,5] + lDW[26,9] >= 0
 ldMC1: - lDX[27,9] + lDW[26,5] - lDW[26,9] <= 0
 ldMC2: - ldx[27,5] + ldw[26,9] >= 0
 ldMC2: - ldx[27,13] + ldw[26,9] >= 0
 ldMC2: - ldx[27,5] - ldx[27,13] + ldw[26,9] <= 0
 ldMC2: lDX[27,5] + lDX[27,13] - lDW[26,9] >= 0
 ldMC2: - lDX[27,5] + lDX[27,13] + lDW[26,9] >= 0
 ldMC2: lDX[27,5] - lDX[27,13] + lDW[26,9] >= 0
 ldMC2: - lDX[27,5] - lDX[27,13] + lDW[26,9] <= 0
 ldMC3: - ldx[27,1] + ldw[26,13] >= 0
 ldMC3: - ldx[27,13] + ldw[26,13] >= 0
 ldMC3: - ldx[27,1] - ldx[27,13] + ldw[26,13] <= 0
 ldMC3: lDX[27,1] + lDX[27,13] - lDW[26,13] >= 0
 ldMC3: - lDX[27,1] + lDX[27,13] + lDW[26,13] >= 0
 ldMC3: lDX[27,1] - lDX[27,13] + lDW[26,13] >= 0
 ldMC3: - lDX[27,1] - lDX[27,13] + lDW[26,13] <= 0
 ldMC0: - ldx[27,6] + ldw[26,2] = 0
 ldMC0: - lDX[27,6] + lDW[26,2] = 0
 ldMC1: ldw[26,6] - ldw[26,10] >= 0
 ldMC1: - ldx[27,10] + ldw[26,6] >= 0
 ldMC1: - ldx[27,10] + ldw[26,6] - ldw[26,10] <= 0
 ldMC1: lDX[27,10] - lDW[26,6] + lDW[26,10] >= 0
 ldMC1: lDX[27,10] + lDW[26,6] - lDW[26,10] >= 0
 ldMC1: - lDX[27,10] + lDW[26,6] + lDW[26,10] >= 0
 ldMC1: - lDX[27,10] + lDW[26,6] - lDW[26,10] <= 0
 ldMC2: - ldx[27,6] + ldw[26,10] >= 0
 ldMC2: - ldx[27,14] + ldw[26,10] >= 0
 ldMC2: - ldx[27,6] - ldx[27,14] + ldw[26,10] <= 0
 ldMC2: lDX[27,6] + lDX[27,14] - lDW[26,10] >= 0
 ldMC2: - lDX[27,6] + lDX[27,14] + lDW[26,10] >= 0
 ldMC2: lDX[27,6] - lDX[27,14] + lDW[26,10] >= 0
 ldMC2: - lDX[27,6] - lDX[27,14] + lDW[26,10] <= 0
 ldMC3: - ldx[27,2] + ldw[26,14] >= 0
 ldMC3: - ldx[27,14] + ldw[26,14] >= 0
 ldMC3: - ldx[27,2] - ldx[27,14] + ldw[26,14] <= 0
 ldMC3: lDX[27,2] + lDX[27,14] - lDW[26,14] >= 0
 ldMC3: - lDX[27,2] + lDX[27,14] + lDW[26,14] >= 0
 ldMC3: lDX[27,2] - lDX[27,14] + lDW[26,14] >= 0
 ldMC3: - lDX[27,2] - lDX[27,14] + lDW[26,14] <= 0
 ldMC0: - ldx[27,7] + ldw[26,3] = 0
 ldMC0: - lDX[27,7] + lDW[26,3] = 0
 ldMC1: ldw[26,7] - ldw[26,11] >= 0
 ldMC1: - ldx[27,11] + ldw[26,7] >= 0
 ldMC1: - ldx[27,11] + ldw[26,7] - ldw[26,11] <= 0
 ldMC1: lDX[27,11] - lDW[26,7] + lDW[26,11] >= 0
 ldMC1: lDX[27,11] + lDW[26,7] - lDW[26,11] >= 0
 ldMC1: - lDX[27,11] + lDW[26,7] + lDW[26,11] >= 0
 ldMC1: - lDX[27,11] + lDW[26,7] - lDW[26,11] <= 0
 ldMC2: - ldx[27,7] + ldw[26,11] >= 0
 ldMC2: - ldx[27,15] + ldw[26,11] >= 0
 ldMC2: - ldx[27,7] - ldx[27,15] + ldw[26,11] <= 0
 ldMC2: lDX[27,7] + lDX[27,15] - lDW[26,11] >= 0
 ldMC2: - lDX[27,7] + lDX[27,15] + lDW[26,11] >= 0
 ldMC2: lDX[27,7] - lDX[27,15] + lDW[26,11] >= 0
 ldMC2: - lDX[27,7] - lDX[27,15] + lDW[26,11] <= 0
 ldMC3: - ldx[27,3] + ldw[26,15] >= 0
 ldMC3: - ldx[27,15] + ldw[26,15] >= 0
 ldMC3: - ldx[27,3] - ldx[27,15] + ldw[26,15] <= 0
 ldMC3: lDX[27,3] + lDX[27,15] - lDW[26,15] >= 0
 ldMC3: - lDX[27,3] + lDX[27,15] + lDW[26,15] >= 0
 ldMC3: lDX[27,3] - lDX[27,15] + lDW[26,15] >= 0
 ldMC3: - lDX[27,3] - lDX[27,15] + lDW[26,15] <= 0
 ldMC0: - ldx[28,4] + ldw[27,0] = 0
 ldMC0: - lDX[28,4] + lDW[27,0] = 0
 ldMC1: ldw[27,4] - ldw[27,8] >= 0
 ldMC1: - ldx[28,8] + ldw[27,4] >= 0
 ldMC1: - ldx[28,8] + ldw[27,4] - ldw[27,8] <= 0
 ldMC1: lDX[28,8] - lDW[27,4] + lDW[27,8] >= 0
 ldMC1: lDX[28,8] + lDW[27,4] - lDW[27,8] >= 0
 ldMC1: - lDX[28,8] + lDW[27,4] + lDW[27,8] >= 0
 ldMC1: - lDX[28,8] + lDW[27,4] - lDW[27,8] <= 0
 ldMC2: - ldx[28,4] + ldw[27,8] >= 0
 ldMC2: - ldx[28,12] + ldw[27,8] >= 0
 ldMC2: - ldx[28,4] - ldx[28,12] + ldw[27,8] <= 0
 ldMC2: lDX[28,4] + lDX[28,12] - lDW[27,8] >= 0
 ldMC2: - lDX[28,4] + lDX[28,12] + lDW[27,8] >= 0
 ldMC2: lDX[28,4] - lDX[28,12] + lDW[27,8] >= 0
 ldMC2: - lDX[28,4] - lDX[28,12] + lDW[27,8] <= 0
 ldMC3: - ldx[28,0] + ldw[27,12] >= 0
 ldMC3: - ldx[28,12] + ldw[27,12] >= 0
 ldMC3: - ldx[28,0] - ldx[28,12] + ldw[27,12] <= 0
 ldMC3: lDX[28,0] + lDX[28,12] - lDW[27,12] >= 0
 ldMC3: - lDX[28,0] + lDX[28,12] + lDW[27,12] >= 0
 ldMC3: lDX[28,0] - lDX[28,12] + lDW[27,12] >= 0
 ldMC3: - lDX[28,0] - lDX[28,12] + lDW[27,12] <= 0
 ldMC0: - ldx[28,5] + ldw[27,1] = 0
 ldMC0: - lDX[28,5] + lDW[27,1] = 0
 ldMC1: ldw[27,5] - ldw[27,9] >= 0
 ldMC1: - ldx[28,9] + ldw[27,5] >= 0
 ldMC1: - ldx[28,9] + ldw[27,5] - ldw[27,9] <= 0
 ldMC1: lDX[28,9] - lDW[27,5] + lDW[27,9] >= 0
 ldMC1: lDX[28,9] + lDW[27,5] - lDW[27,9] >= 0
 ldMC1: - lDX[28,9] + lDW[27,5] + lDW[27,9] >= 0
 ldMC1: - lDX[28,9] + lDW[27,5] - lDW[27,9] <= 0
 ldMC2: - ldx[28,5] + ldw[27,9] >= 0
 ldMC2: - ldx[28,13] + ldw[27,9] >= 0
 ldMC2: - ldx[28,5] - ldx[28,13] + ldw[27,9] <= 0
 ldMC2: lDX[28,5] + lDX[28,13] - lDW[27,9] >= 0
 ldMC2: - lDX[28,5] + lDX[28,13] + lDW[27,9] >= 0
 ldMC2: lDX[28,5] - lDX[28,13] + lDW[27,9] >= 0
 ldMC2: - lDX[28,5] - lDX[28,13] + lDW[27,9] <= 0
 ldMC3: - ldx[28,1] + ldw[27,13] >= 0
 ldMC3: - ldx[28,13] + ldw[27,13] >= 0
 ldMC3: - ldx[28,1] - ldx[28,13] + ldw[27,13] <= 0
 ldMC3: lDX[28,1] + lDX[28,13] - lDW[27,13] >= 0
 ldMC3: - lDX[28,1] + lDX[28,13] + lDW[27,13] >= 0
 ldMC3: lDX[28,1] - lDX[28,13] + lDW[27,13] >= 0
 ldMC3: - lDX[28,1] - lDX[28,13] + lDW[27,13] <= 0
 ldMC0: - ldx[28,6] + ldw[27,2] = 0
 ldMC0: - lDX[28,6] + lDW[27,2] = 0
 ldMC1: ldw[27,6] - ldw[27,10] >= 0
 ldMC1: - ldx[28,10] + ldw[27,6] >= 0
 ldMC1: - ldx[28,10] + ldw[27,6] - ldw[27,10] <= 0
 ldMC1: lDX[28,10] - lDW[27,6] + lDW[27,10] >= 0
 ldMC1: lDX[28,10] + lDW[27,6] - lDW[27,10] >= 0
 ldMC1: - lDX[28,10] + lDW[27,6] + lDW[27,10] >= 0
 ldMC1: - lDX[28,10] + lDW[27,6] - lDW[27,10] <= 0
 ldMC2: - ldx[28,6] + ldw[27,10] >= 0
 ldMC2: - ldx[28,14] + ldw[27,10] >= 0
 ldMC2: - ldx[28,6] - ldx[28,14] + ldw[27,10] <= 0
 ldMC2: lDX[28,6] + lDX[28,14] - lDW[27,10] >= 0
 ldMC2: - lDX[28,6] + lDX[28,14] + lDW[27,10] >= 0
 ldMC2: lDX[28,6] - lDX[28,14] + lDW[27,10] >= 0
 ldMC2: - lDX[28,6] - lDX[28,14] + lDW[27,10] <= 0
 ldMC3: - ldx[28,2] + ldw[27,14] >= 0
 ldMC3: - ldx[28,14] + ldw[27,14] >= 0
 ldMC3: - ldx[28,2] - ldx[28,14] + ldw[27,14] <= 0
 ldMC3: lDX[28,2] + lDX[28,14] - lDW[27,14] >= 0
 ldMC3: - lDX[28,2] + lDX[28,14] + lDW[27,14] >= 0
 ldMC3: lDX[28,2] - lDX[28,14] + lDW[27,14] >= 0
 ldMC3: - lDX[28,2] - lDX[28,14] + lDW[27,14] <= 0
 ldMC0: - ldx[28,7] + ldw[27,3] = 0
 ldMC0: - lDX[28,7] + lDW[27,3] = 0
 ldMC1: ldw[27,7] - ldw[27,11] >= 0
 ldMC1: - ldx[28,11] + ldw[27,7] >= 0
 ldMC1: - ldx[28,11] + ldw[27,7] - ldw[27,11] <= 0
 ldMC1: lDX[28,11] - lDW[27,7] + lDW[27,11] >= 0
 ldMC1: lDX[28,11] + lDW[27,7] - lDW[27,11] >= 0
 ldMC1: - lDX[28,11] + lDW[27,7] + lDW[27,11] >= 0
 ldMC1: - lDX[28,11] + lDW[27,7] - lDW[27,11] <= 0
 ldMC2: - ldx[28,7] + ldw[27,11] >= 0
 ldMC2: - ldx[28,15] + ldw[27,11] >= 0
 ldMC2: - ldx[28,7] - ldx[28,15] + ldw[27,11] <= 0
 ldMC2: lDX[28,7] + lDX[28,15] - lDW[27,11] >= 0
 ldMC2: - lDX[28,7] + lDX[28,15] + lDW[27,11] >= 0
 ldMC2: lDX[28,7] - lDX[28,15] + lDW[27,11] >= 0
 ldMC2: - lDX[28,7] - lDX[28,15] + lDW[27,11] <= 0
 ldMC3: - ldx[28,3] + ldw[27,15] >= 0
 ldMC3: - ldx[28,15] + ldw[27,15] >= 0
 ldMC3: - ldx[28,3] - ldx[28,15] + ldw[27,15] <= 0
 ldMC3: lDX[28,3] + lDX[28,15] - lDW[27,15] >= 0
 ldMC3: - lDX[28,3] + lDX[28,15] + lDW[27,15] >= 0
 ldMC3: lDX[28,3] - lDX[28,15] + lDW[27,15] >= 0
 ldMC3: - lDX[28,3] - lDX[28,15] + lDW[27,15] <= 0
 leMC0: ldx[29,0] - ldw[28,12] >= 0
 leMC0: ldx[29,0] - ldx[29,12] >= 0
 leMC0: ldx[29,0] - ldx[29,12] - ldw[28,12] <= 0
 leMC0: - lDX[29,0] + lDX[29,12] + lDW[28,12] >= 0
 leMC0: lDX[29,0] + lDX[29,12] - lDW[28,12] >= 0
 leMC0: lDX[29,0] - lDX[29,12] + lDW[28,12] >= 0
 leMC0: lDX[29,0] - lDX[29,12] - lDW[28,12] <= 0
 leMC1: ldx[29,4] - ldw[28,0] = 0
 leMC1: lDX[29,4] - lDW[28,0] = 0
 leMC2: ldx[29,8] - ldw[28,4] >= 0
 leMC2: ldx[29,8] - ldw[28,8] >= 0
 leMC2: ldx[29,8] - ldw[28,4] - ldw[28,8] <= 0
 leMC2: - lDX[29,8] + lDW[28,4] + lDW[28,8] >= 0
 leMC2: lDX[29,8] - lDW[28,4] + lDW[28,8] >= 0
 leMC2: lDX[29,8] + lDW[28,4] - lDW[28,8] >= 0
 leMC2: lDX[29,8] - lDW[28,4] - lDW[28,8] <= 0
 leMC3: ldx[29,12] - ldw[28,0] >= 0
 leMC3: ldx[29,12] - ldw[28,8] >= 0
 leMC3: ldx[29,12] - ldw[28,0] - ldw[28,8] <= 0
 leMC3: - lDX[29,12] + lDW[28,0] + lDW[28,8] >= 0
 leMC3: lDX[29,12] - lDW[28,0] + lDW[28,8] >= 0
 leMC3: lDX[29,12] + lDW[28,0] - lDW[28,8] >= 0
 leMC3: lDX[29,12] - lDW[28,0] - lDW[28,8] <= 0
 leMC0: ldx[29,1] - ldw[28,13] >= 0
 leMC0: ldx[29,1] - ldx[29,13] >= 0
 leMC0: ldx[29,1] - ldx[29,13] - ldw[28,13] <= 0
 leMC0: - lDX[29,1] + lDX[29,13] + lDW[28,13] >= 0
 leMC0: lDX[29,1] + lDX[29,13] - lDW[28,13] >= 0
 leMC0: lDX[29,1] - lDX[29,13] + lDW[28,13] >= 0
 leMC0: lDX[29,1] - lDX[29,13] - lDW[28,13] <= 0
 leMC1: ldx[29,5] - ldw[28,1] = 0
 leMC1: lDX[29,5] - lDW[28,1] = 0
 leMC2: ldx[29,9] - ldw[28,5] >= 0
 leMC2: ldx[29,9] - ldw[28,9] >= 0
 leMC2: ldx[29,9] - ldw[28,5] - ldw[28,9] <= 0
 leMC2: - lDX[29,9] + lDW[28,5] + lDW[28,9] >= 0
 leMC2: lDX[29,9] - lDW[28,5] + lDW[28,9] >= 0
 leMC2: lDX[29,9] + lDW[28,5] - lDW[28,9] >= 0
 leMC2: lDX[29,9] - lDW[28,5] - lDW[28,9] <= 0
 leMC3: ldx[29,13] - ldw[28,1] >= 0
 leMC3: ldx[29,13] - ldw[28,9] >= 0
 leMC3: ldx[29,13] - ldw[28,1] - ldw[28,9] <= 0
 leMC3: - lDX[29,13] + lDW[28,1] + lDW[28,9] >= 0
 leMC3: lDX[29,13] - lDW[28,1] + lDW[28,9] >= 0
 leMC3: lDX[29,13] + lDW[28,1] - lDW[28,9] >= 0
 leMC3: lDX[29,13] - lDW[28,1] - lDW[28,9] <= 0
 leMC0: ldx[29,2] - ldw[28,14] >= 0
 leMC0: ldx[29,2] - ldx[29,14] >= 0
 leMC0: ldx[29,2] - ldx[29,14] - ldw[28,14] <= 0
 leMC0: - lDX[29,2] + lDX[29,14] + lDW[28,14] >= 0
 leMC0: lDX[29,2] + lDX[29,14] - lDW[28,14] >= 0
 leMC0: lDX[29,2] - lDX[29,14] + lDW[28,14] >= 0
 leMC0: lDX[29,2] - lDX[29,14] - lDW[28,14] <= 0
 leMC1: ldx[29,6] - ldw[28,2] = 0
 leMC1: lDX[29,6] - lDW[28,2] = 0
 leMC2: ldx[29,10] - ldw[28,6] >= 0
 leMC2: ldx[29,10] - ldw[28,10] >= 0
 leMC2: ldx[29,10] - ldw[28,6] - ldw[28,10] <= 0
 leMC2: - lDX[29,10] + lDW[28,6] + lDW[28,10] >= 0
 leMC2: lDX[29,10] - lDW[28,6] + lDW[28,10] >= 0
 leMC2: lDX[29,10] + lDW[28,6] - lDW[28,10] >= 0
 leMC2: lDX[29,10] - lDW[28,6] - lDW[28,10] <= 0
 leMC3: ldx[29,14] - ldw[28,2] >= 0
 leMC3: ldx[29,14] - ldw[28,10] >= 0
 leMC3: ldx[29,14] - ldw[28,2] - ldw[28,10] <= 0
 leMC3: - lDX[29,14] + lDW[28,2] + lDW[28,10] >= 0
 leMC3: lDX[29,14] - lDW[28,2] + lDW[28,10] >= 0
 leMC3: lDX[29,14] + lDW[28,2] - lDW[28,10] >= 0
 leMC3: lDX[29,14] - lDW[28,2] - lDW[28,10] <= 0
 leMC0: ldx[29,3] - ldw[28,15] >= 0
 leMC0: ldx[29,3] - ldx[29,15] >= 0
 leMC0: ldx[29,3] - ldx[29,15] - ldw[28,15] <= 0
 leMC0: - lDX[29,3] + lDX[29,15] + lDW[28,15] >= 0
 leMC0: lDX[29,3] + lDX[29,15] - lDW[28,15] >= 0
 leMC0: lDX[29,3] - lDX[29,15] + lDW[28,15] >= 0
 leMC0: lDX[29,3] - lDX[29,15] - lDW[28,15] <= 0
 leMC1: ldx[29,7] - ldw[28,3] = 0
 leMC1: lDX[29,7] - lDW[28,3] = 0
 leMC2: ldx[29,11] - ldw[28,7] >= 0
 leMC2: ldx[29,11] - ldw[28,11] >= 0
 leMC2: ldx[29,11] - ldw[28,7] - ldw[28,11] <= 0
 leMC2: - lDX[29,11] + lDW[28,7] + lDW[28,11] >= 0
 leMC2: lDX[29,11] - lDW[28,7] + lDW[28,11] >= 0
 leMC2: lDX[29,11] + lDW[28,7] - lDW[28,11] >= 0
 leMC2: lDX[29,11] - lDW[28,7] - lDW[28,11] <= 0
 leMC3: ldx[29,15] - ldw[28,3] >= 0
 leMC3: ldx[29,15] - ldw[28,11] >= 0
 leMC3: ldx[29,15] - ldw[28,3] - ldw[28,11] <= 0
 leMC3: - lDX[29,15] + lDW[28,3] + lDW[28,11] >= 0
 leMC3: lDX[29,15] - lDW[28,3] + lDW[28,11] >= 0
 leMC3: lDX[29,15] + lDW[28,3] - lDW[28,11] >= 0
 leMC3: lDX[29,15] - lDW[28,3] - lDW[28,11] <= 0
 leMC0: ldx[30,0] - ldw[29,12] >= 0
 leMC0: ldx[30,0] - ldx[30,12] >= 0
 leMC0: ldx[30,0] - ldx[30,12] - ldw[29,12] <= 0
 leMC0: - lDX[30,0] + lDX[30,12] + lDW[29,12] >= 0
 leMC0: lDX[30,0] + lDX[30,12] - lDW[29,12] >= 0
 leMC0: lDX[30,0] - lDX[30,12] + lDW[29,12] >= 0
 leMC0: lDX[30,0] - lDX[30,12] - lDW[29,12] <= 0
 leMC1: ldx[30,4] - ldw[29,0] = 0
 leMC1: lDX[30,4] - lDW[29,0] = 0
 leMC2: ldx[30,8] - ldw[29,4] >= 0
 leMC2: ldx[30,8] - ldw[29,8] >= 0
 leMC2: ldx[30,8] - ldw[29,4] - ldw[29,8] <= 0
 leMC2: - lDX[30,8] + lDW[29,4] + lDW[29,8] >= 0
 leMC2: lDX[30,8] - lDW[29,4] + lDW[29,8] >= 0
 leMC2: lDX[30,8] + lDW[29,4] - lDW[29,8] >= 0
 leMC2: lDX[30,8] - lDW[29,4] - lDW[29,8] <= 0
 leMC3: ldx[30,12] - ldw[29,0] >= 0
 leMC3: ldx[30,12] - ldw[29,8] >= 0
 leMC3: ldx[30,12] - ldw[29,0] - ldw[29,8] <= 0
 leMC3: - lDX[30,12] + lDW[29,0] + lDW[29,8] >= 0
 leMC3: lDX[30,12] - lDW[29,0] + lDW[29,8] >= 0
 leMC3: lDX[30,12] + lDW[29,0] - lDW[29,8] >= 0
 leMC3: lDX[30,12] - lDW[29,0] - lDW[29,8] <= 0
 leMC0: ldx[30,1] - ldw[29,13] >= 0
 leMC0: ldx[30,1] - ldx[30,13] >= 0
 leMC0: ldx[30,1] - ldx[30,13] - ldw[29,13] <= 0
 leMC0: - lDX[30,1] + lDX[30,13] + lDW[29,13] >= 0
 leMC0: lDX[30,1] + lDX[30,13] - lDW[29,13] >= 0
 leMC0: lDX[30,1] - lDX[30,13] + lDW[29,13] >= 0
 leMC0: lDX[30,1] - lDX[30,13] - lDW[29,13] <= 0
 leMC1: ldx[30,5] - ldw[29,1] = 0
 leMC1: lDX[30,5] - lDW[29,1] = 0
 leMC2: ldx[30,9] - ldw[29,5] >= 0
 leMC2: ldx[30,9] - ldw[29,9] >= 0
 leMC2: ldx[30,9] - ldw[29,5] - ldw[29,9] <= 0
 leMC2: - lDX[30,9] + lDW[29,5] + lDW[29,9] >= 0
 leMC2: lDX[30,9] - lDW[29,5] + lDW[29,9] >= 0
 leMC2: lDX[30,9] + lDW[29,5] - lDW[29,9] >= 0
 leMC2: lDX[30,9] - lDW[29,5] - lDW[29,9] <= 0
 leMC3: ldx[30,13] - ldw[29,1] >= 0
 leMC3: ldx[30,13] - ldw[29,9] >= 0
 leMC3: ldx[30,13] - ldw[29,1] - ldw[29,9] <= 0
 leMC3: - lDX[30,13] + lDW[29,1] + lDW[29,9] >= 0
 leMC3: lDX[30,13] - lDW[29,1] + lDW[29,9] >= 0
 leMC3: lDX[30,13] + lDW[29,1] - lDW[29,9] >= 0
 leMC3: lDX[30,13] - lDW[29,1] - lDW[29,9] <= 0
 leMC0: ldx[30,2] - ldw[29,14] >= 0
 leMC0: ldx[30,2] - ldx[30,14] >= 0
 leMC0: ldx[30,2] - ldx[30,14] - ldw[29,14] <= 0
 leMC0: - lDX[30,2] + lDX[30,14] + lDW[29,14] >= 0
 leMC0: lDX[30,2] + lDX[30,14] - lDW[29,14] >= 0
 leMC0: lDX[30,2] - lDX[30,14] + lDW[29,14] >= 0
 leMC0: lDX[30,2] - lDX[30,14] - lDW[29,14] <= 0
 leMC1: ldx[30,6] - ldw[29,2] = 0
 leMC1: lDX[30,6] - lDW[29,2] = 0
 leMC2: ldx[30,10] - ldw[29,6] >= 0
 leMC2: ldx[30,10] - ldw[29,10] >= 0
 leMC2: ldx[30,10] - ldw[29,6] - ldw[29,10] <= 0
 leMC2: - lDX[30,10] + lDW[29,6] + lDW[29,10] >= 0
 leMC2: lDX[30,10] - lDW[29,6] + lDW[29,10] >= 0
 leMC2: lDX[30,10] + lDW[29,6] - lDW[29,10] >= 0
 leMC2: lDX[30,10] - lDW[29,6] - lDW[29,10] <= 0
 leMC3: ldx[30,14] - ldw[29,2] >= 0
 leMC3: ldx[30,14] - ldw[29,10] >= 0
 leMC3: ldx[30,14] - ldw[29,2] - ldw[29,10] <= 0
 leMC3: - lDX[30,14] + lDW[29,2] + lDW[29,10] >= 0
 leMC3: lDX[30,14] - lDW[29,2] + lDW[29,10] >= 0
 leMC3: lDX[30,14] + lDW[29,2] - lDW[29,10] >= 0
 leMC3: lDX[30,14] - lDW[29,2] - lDW[29,10] <= 0
 leMC0: ldx[30,3] - ldw[29,15] >= 0
 leMC0: ldx[30,3] - ldx[30,15] >= 0
 leMC0: ldx[30,3] - ldx[30,15] - ldw[29,15] <= 0
 leMC0: - lDX[30,3] + lDX[30,15] + lDW[29,15] >= 0
 leMC0: lDX[30,3] + lDX[30,15] - lDW[29,15] >= 0
 leMC0: lDX[30,3] - lDX[30,15] + lDW[29,15] >= 0
 leMC0: lDX[30,3] - lDX[30,15] - lDW[29,15] <= 0
 leMC1: ldx[30,7] - ldw[29,3] = 0
 leMC1: lDX[30,7] - lDW[29,3] = 0
 leMC2: ldx[30,11] - ldw[29,7] >= 0
 leMC2: ldx[30,11] - ldw[29,11] >= 0
 leMC2: ldx[30,11] - ldw[29,7] - ldw[29,11] <= 0
 leMC2: - lDX[30,11] + lDW[29,7] + lDW[29,11] >= 0
 leMC2: lDX[30,11] - lDW[29,7] + lDW[29,11] >= 0
 leMC2: lDX[30,11] + lDW[29,7] - lDW[29,11] >= 0
 leMC2: lDX[30,11] - lDW[29,7] - lDW[29,11] <= 0
 leMC3: ldx[30,15] - ldw[29,3] >= 0
 leMC3: ldx[30,15] - ldw[29,11] >= 0
 leMC3: ldx[30,15] - ldw[29,3] - ldw[29,11] <= 0
 leMC3: - lDX[30,15] + lDW[29,3] + lDW[29,11] >= 0
 leMC3: lDX[30,15] - lDW[29,3] + lDW[29,11] >= 0
 leMC3: lDX[30,15] + lDW[29,3] - lDW[29,11] >= 0
 leMC3: lDX[30,15] - lDW[29,3] - lDW[29,11] <= 0
 leMC0: ldx[31,0] - ldw[30,12] >= 0
 leMC0: ldx[31,0] - ldx[31,12] >= 0
 leMC0: ldx[31,0] - ldx[31,12] - ldw[30,12] <= 0
 leMC0: - lDX[31,0] + lDX[31,12] + lDW[30,12] >= 0
 leMC0: lDX[31,0] + lDX[31,12] - lDW[30,12] >= 0
 leMC0: lDX[31,0] - lDX[31,12] + lDW[30,12] >= 0
 leMC0: lDX[31,0] - lDX[31,12] - lDW[30,12] <= 0
 leMC1: ldx[31,4] - ldw[30,0] = 0
 leMC1: lDX[31,4] - lDW[30,0] = 0
 leMC2: ldx[31,8] - ldw[30,4] >= 0
 leMC2: ldx[31,8] - ldw[30,8] >= 0
 leMC2: ldx[31,8] - ldw[30,4] - ldw[30,8] <= 0
 leMC2: - lDX[31,8] + lDW[30,4] + lDW[30,8] >= 0
 leMC2: lDX[31,8] - lDW[30,4] + lDW[30,8] >= 0
 leMC2: lDX[31,8] + lDW[30,4] - lDW[30,8] >= 0
 leMC2: lDX[31,8] - lDW[30,4] - lDW[30,8] <= 0
 leMC3: ldx[31,12] - ldw[30,0] >= 0
 leMC3: ldx[31,12] - ldw[30,8] >= 0
 leMC3: ldx[31,12] - ldw[30,0] - ldw[30,8] <= 0
 leMC3: - lDX[31,12] + lDW[30,0] + lDW[30,8] >= 0
 leMC3: lDX[31,12] - lDW[30,0] + lDW[30,8] >= 0
 leMC3: lDX[31,12] + lDW[30,0] - lDW[30,8] >= 0
 leMC3: lDX[31,12] - lDW[30,0] - lDW[30,8] <= 0
 leMC0: ldx[31,1] - ldw[30,13] >= 0
 leMC0: ldx[31,1] - ldx[31,13] >= 0
 leMC0: ldx[31,1] - ldx[31,13] - ldw[30,13] <= 0
 leMC0: - lDX[31,1] + lDX[31,13] + lDW[30,13] >= 0
 leMC0: lDX[31,1] + lDX[31,13] - lDW[30,13] >= 0
 leMC0: lDX[31,1] - lDX[31,13] + lDW[30,13] >= 0
 leMC0: lDX[31,1] - lDX[31,13] - lDW[30,13] <= 0
 leMC1: ldx[31,5] - ldw[30,1] = 0
 leMC1: lDX[31,5] - lDW[30,1] = 0
 leMC2: ldx[31,9] - ldw[30,5] >= 0
 leMC2: ldx[31,9] - ldw[30,9] >= 0
 leMC2: ldx[31,9] - ldw[30,5] - ldw[30,9] <= 0
 leMC2: - lDX[31,9] + lDW[30,5] + lDW[30,9] >= 0
 leMC2: lDX[31,9] - lDW[30,5] + lDW[30,9] >= 0
 leMC2: lDX[31,9] + lDW[30,5] - lDW[30,9] >= 0
 leMC2: lDX[31,9] - lDW[30,5] - lDW[30,9] <= 0
 leMC3: ldx[31,13] - ldw[30,1] >= 0
 leMC3: ldx[31,13] - ldw[30,9] >= 0
 leMC3: ldx[31,13] - ldw[30,1] - ldw[30,9] <= 0
 leMC3: - lDX[31,13] + lDW[30,1] + lDW[30,9] >= 0
 leMC3: lDX[31,13] - lDW[30,1] + lDW[30,9] >= 0
 leMC3: lDX[31,13] + lDW[30,1] - lDW[30,9] >= 0
 leMC3: lDX[31,13] - lDW[30,1] - lDW[30,9] <= 0
 leMC0: ldx[31,2] - ldw[30,14] >= 0
 leMC0: ldx[31,2] - ldx[31,14] >= 0
 leMC0: ldx[31,2] - ldx[31,14] - ldw[30,14] <= 0
 leMC0: - lDX[31,2] + lDX[31,14] + lDW[30,14] >= 0
 leMC0: lDX[31,2] + lDX[31,14] - lDW[30,14] >= 0
 leMC0: lDX[31,2] - lDX[31,14] + lDW[30,14] >= 0
 leMC0: lDX[31,2] - lDX[31,14] - lDW[30,14] <= 0
 leMC1: ldx[31,6] - ldw[30,2] = 0
 leMC1: lDX[31,6] - lDW[30,2] = 0
 leMC2: ldx[31,10] - ldw[30,6] >= 0
 leMC2: ldx[31,10] - ldw[30,10] >= 0
 leMC2: ldx[31,10] - ldw[30,6] - ldw[30,10] <= 0
 leMC2: - lDX[31,10] + lDW[30,6] + lDW[30,10] >= 0
 leMC2: lDX[31,10] - lDW[30,6] + lDW[30,10] >= 0
 leMC2: lDX[31,10] + lDW[30,6] - lDW[30,10] >= 0
 leMC2: lDX[31,10] - lDW[30,6] - lDW[30,10] <= 0
 leMC3: ldx[31,14] - ldw[30,2] >= 0
 leMC3: ldx[31,14] - ldw[30,10] >= 0
 leMC3: ldx[31,14] - ldw[30,2] - ldw[30,10] <= 0
 leMC3: - lDX[31,14] + lDW[30,2] + lDW[30,10] >= 0
 leMC3: lDX[31,14] - lDW[30,2] + lDW[30,10] >= 0
 leMC3: lDX[31,14] + lDW[30,2] - lDW[30,10] >= 0
 leMC3: lDX[31,14] - lDW[30,2] - lDW[30,10] <= 0
 leMC0: ldx[31,3] - ldw[30,15] >= 0
 leMC0: ldx[31,3] - ldx[31,15] >= 0
 leMC0: ldx[31,3] - ldx[31,15] - ldw[30,15] <= 0
 leMC0: - lDX[31,3] + lDX[31,15] + lDW[30,15] >= 0
 leMC0: lDX[31,3] + lDX[31,15] - lDW[30,15] >= 0
 leMC0: lDX[31,3] - lDX[31,15] + lDW[30,15] >= 0
 leMC0: lDX[31,3] - lDX[31,15] - lDW[30,15] <= 0
 leMC1: ldx[31,7] - ldw[30,3] = 0
 leMC1: lDX[31,7] - lDW[30,3] = 0
 leMC2: ldx[31,11] - ldw[30,7] >= 0
 leMC2: ldx[31,11] - ldw[30,11] >= 0
 leMC2: ldx[31,11] - ldw[30,7] - ldw[30,11] <= 0
 leMC2: - lDX[31,11] + lDW[30,7] + lDW[30,11] >= 0
 leMC2: lDX[31,11] - lDW[30,7] + lDW[30,11] >= 0
 leMC2: lDX[31,11] + lDW[30,7] - lDW[30,11] >= 0
 leMC2: lDX[31,11] - lDW[30,7] - lDW[30,11] <= 0
 leMC3: ldx[31,15] - ldw[30,3] >= 0
 leMC3: ldx[31,15] - ldw[30,11] >= 0
 leMC3: ldx[31,15] - ldw[30,3] - ldw[30,11] <= 0
 leMC3: - lDX[31,15] + lDW[30,3] + lDW[30,11] >= 0
 leMC3: lDX[31,15] - lDW[30,3] + lDW[30,11] >= 0
 leMC3: lDX[31,15] + lDW[30,3] - lDW[30,11] >= 0
 leMC3: lDX[31,15] - lDW[30,3] - lDW[30,11] <= 0
 leMC0: ldx[32,0] - ldw[31,12] >= 0
 leMC0: ldx[32,0] - ldx[32,12] >= 0
 leMC0: ldx[32,0] - ldx[32,12] - ldw[31,12] <= 0
 leMC0: - lDX[32,0] + lDX[32,12] + lDW[31,12] >= 0
 leMC0: lDX[32,0] + lDX[32,12] - lDW[31,12] >= 0
 leMC0: lDX[32,0] - lDX[32,12] + lDW[31,12] >= 0
 leMC0: lDX[32,0] - lDX[32,12] - lDW[31,12] <= 0
 leMC1: ldx[32,4] - ldw[31,0] = 0
 leMC1: lDX[32,4] - lDW[31,0] = 0
 leMC2: ldx[32,8] - ldw[31,4] >= 0
 leMC2: ldx[32,8] - ldw[31,8] >= 0
 leMC2: ldx[32,8] - ldw[31,4] - ldw[31,8] <= 0
 leMC2: - lDX[32,8] + lDW[31,4] + lDW[31,8] >= 0
 leMC2: lDX[32,8] - lDW[31,4] + lDW[31,8] >= 0
 leMC2: lDX[32,8] + lDW[31,4] - lDW[31,8] >= 0
 leMC2: lDX[32,8] - lDW[31,4] - lDW[31,8] <= 0
 leMC3: ldx[32,12] - ldw[31,0] >= 0
 leMC3: ldx[32,12] - ldw[31,8] >= 0
 leMC3: ldx[32,12] - ldw[31,0] - ldw[31,8] <= 0
 leMC3: - lDX[32,12] + lDW[31,0] + lDW[31,8] >= 0
 leMC3: lDX[32,12] - lDW[31,0] + lDW[31,8] >= 0
 leMC3: lDX[32,12] + lDW[31,0] - lDW[31,8] >= 0
 leMC3: lDX[32,12] - lDW[31,0] - lDW[31,8] <= 0
 leMC0: ldx[32,1] - ldw[31,13] >= 0
 leMC0: ldx[32,1] - ldx[32,13] >= 0
 leMC0: ldx[32,1] - ldx[32,13] - ldw[31,13] <= 0
 leMC0: - lDX[32,1] + lDX[32,13] + lDW[31,13] >= 0
 leMC0: lDX[32,1] + lDX[32,13] - lDW[31,13] >= 0
 leMC0: lDX[32,1] - lDX[32,13] + lDW[31,13] >= 0
 leMC0: lDX[32,1] - lDX[32,13] - lDW[31,13] <= 0
 leMC1: ldx[32,5] - ldw[31,1] = 0
 leMC1: lDX[32,5] - lDW[31,1] = 0
 leMC2: ldx[32,9] - ldw[31,5] >= 0
 leMC2: ldx[32,9] - ldw[31,9] >= 0
 leMC2: ldx[32,9] - ldw[31,5] - ldw[31,9] <= 0
 leMC2: - lDX[32,9] + lDW[31,5] + lDW[31,9] >= 0
 leMC2: lDX[32,9] - lDW[31,5] + lDW[31,9] >= 0
 leMC2: lDX[32,9] + lDW[31,5] - lDW[31,9] >= 0
 leMC2: lDX[32,9] - lDW[31,5] - lDW[31,9] <= 0
 leMC3: ldx[32,13] - ldw[31,1] >= 0
 leMC3: ldx[32,13] - ldw[31,9] >= 0
 leMC3: ldx[32,13] - ldw[31,1] - ldw[31,9] <= 0
 leMC3: - lDX[32,13] + lDW[31,1] + lDW[31,9] >= 0
 leMC3: lDX[32,13] - lDW[31,1] + lDW[31,9] >= 0
 leMC3: lDX[32,13] + lDW[31,1] - lDW[31,9] >= 0
 leMC3: lDX[32,13] - lDW[31,1] - lDW[31,9] <= 0
 leMC0: ldx[32,2] - ldw[31,14] >= 0
 leMC0: ldx[32,2] - ldx[32,14] >= 0
 leMC0: ldx[32,2] - ldx[32,14] - ldw[31,14] <= 0
 leMC0: - lDX[32,2] + lDX[32,14] + lDW[31,14] >= 0
 leMC0: lDX[32,2] + lDX[32,14] - lDW[31,14] >= 0
 leMC0: lDX[32,2] - lDX[32,14] + lDW[31,14] >= 0
 leMC0: lDX[32,2] - lDX[32,14] - lDW[31,14] <= 0
 leMC1: ldx[32,6] - ldw[31,2] = 0
 leMC1: lDX[32,6] - lDW[31,2] = 0
 leMC2: ldx[32,10] - ldw[31,6] >= 0
 leMC2: ldx[32,10] - ldw[31,10] >= 0
 leMC2: ldx[32,10] - ldw[31,6] - ldw[31,10] <= 0
 leMC2: - lDX[32,10] + lDW[31,6] + lDW[31,10] >= 0
 leMC2: lDX[32,10] - lDW[31,6] + lDW[31,10] >= 0
 leMC2: lDX[32,10] + lDW[31,6] - lDW[31,10] >= 0
 leMC2: lDX[32,10] - lDW[31,6] - lDW[31,10] <= 0
 leMC3: ldx[32,14] - ldw[31,2] >= 0
 leMC3: ldx[32,14] - ldw[31,10] >= 0
 leMC3: ldx[32,14] - ldw[31,2] - ldw[31,10] <= 0
 leMC3: - lDX[32,14] + lDW[31,2] + lDW[31,10] >= 0
 leMC3: lDX[32,14] - lDW[31,2] + lDW[31,10] >= 0
 leMC3: lDX[32,14] + lDW[31,2] - lDW[31,10] >= 0
 leMC3: lDX[32,14] - lDW[31,2] - lDW[31,10] <= 0
 leMC0: ldx[32,3] - ldw[31,15] >= 0
 leMC0: ldx[32,3] - ldx[32,15] >= 0
 leMC0: ldx[32,3] - ldx[32,15] - ldw[31,15] <= 0
 leMC0: - lDX[32,3] + lDX[32,15] + lDW[31,15] >= 0
 leMC0: lDX[32,3] + lDX[32,15] - lDW[31,15] >= 0
 leMC0: lDX[32,3] - lDX[32,15] + lDW[31,15] >= 0
 leMC0: lDX[32,3] - lDX[32,15] - lDW[31,15] <= 0
 leMC1: ldx[32,7] - ldw[31,3] = 0
 leMC1: lDX[32,7] - lDW[31,3] = 0
 leMC2: ldx[32,11] - ldw[31,7] >= 0
 leMC2: ldx[32,11] - ldw[31,11] >= 0
 leMC2: ldx[32,11] - ldw[31,7] - ldw[31,11] <= 0
 leMC2: - lDX[32,11] + lDW[31,7] + lDW[31,11] >= 0
 leMC2: lDX[32,11] - lDW[31,7] + lDW[31,11] >= 0
 leMC2: lDX[32,11] + lDW[31,7] - lDW[31,11] >= 0
 leMC2: lDX[32,11] - lDW[31,7] - lDW[31,11] <= 0
 leMC3: ldx[32,15] - ldw[31,3] >= 0
 leMC3: ldx[32,15] - ldw[31,11] >= 0
 leMC3: ldx[32,15] - ldw[31,3] - ldw[31,11] <= 0
 leMC3: - lDX[32,15] + lDW[31,3] + lDW[31,11] >= 0
 leMC3: lDX[32,15] - lDW[31,3] + lDW[31,11] >= 0
 leMC3: lDX[32,15] + lDW[31,3] - lDW[31,11] >= 0
 leMC3: lDX[32,15] - lDW[31,3] - lDW[31,11] <= 0
 leMC0: ldx[33,0] - ldw[32,12] >= 0
 leMC0: ldx[33,0] - ldx[33,12] >= 0
 leMC0: ldx[33,0] - ldx[33,12] - ldw[32,12] <= 0
 leMC0: - lDX[33,0] + lDX[33,12] + lDW[32,12] >= 0
 leMC0: lDX[33,0] + lDX[33,12] - lDW[32,12] >= 0
 leMC0: lDX[33,0] - lDX[33,12] + lDW[32,12] >= 0
 leMC0: lDX[33,0] - lDX[33,12] - lDW[32,12] <= 0
 leMC1: ldx[33,4] - ldw[32,0] = 0
 leMC1: lDX[33,4] - lDW[32,0] = 0
 leMC2: ldx[33,8] - ldw[32,4] >= 0
 leMC2: ldx[33,8] - ldw[32,8] >= 0
 leMC2: ldx[33,8] - ldw[32,4] - ldw[32,8] <= 0
 leMC2: - lDX[33,8] + lDW[32,4] + lDW[32,8] >= 0
 leMC2: lDX[33,8] - lDW[32,4] + lDW[32,8] >= 0
 leMC2: lDX[33,8] + lDW[32,4] - lDW[32,8] >= 0
 leMC2: lDX[33,8] - lDW[32,4] - lDW[32,8] <= 0
 leMC3: ldx[33,12] - ldw[32,0] >= 0
 leMC3: ldx[33,12] - ldw[32,8] >= 0
 leMC3: ldx[33,12] - ldw[32,0] - ldw[32,8] <= 0
 leMC3: - lDX[33,12] + lDW[32,0] + lDW[32,8] >= 0
 leMC3: lDX[33,12] - lDW[32,0] + lDW[32,8] >= 0
 leMC3: lDX[33,12] + lDW[32,0] - lDW[32,8] >= 0
 leMC3: lDX[33,12] - lDW[32,0] - lDW[32,8] <= 0
 leMC0: ldx[33,1] - ldw[32,13] >= 0
 leMC0: ldx[33,1] - ldx[33,13] >= 0
 leMC0: ldx[33,1] - ldx[33,13] - ldw[32,13] <= 0
 leMC0: - lDX[33,1] + lDX[33,13] + lDW[32,13] >= 0
 leMC0: lDX[33,1] + lDX[33,13] - lDW[32,13] >= 0
 leMC0: lDX[33,1] - lDX[33,13] + lDW[32,13] >= 0
 leMC0: lDX[33,1] - lDX[33,13] - lDW[32,13] <= 0
 leMC1: ldx[33,5] - ldw[32,1] = 0
 leMC1: lDX[33,5] - lDW[32,1] = 0
 leMC2: ldx[33,9] - ldw[32,5] >= 0
 leMC2: ldx[33,9] - ldw[32,9] >= 0
 leMC2: ldx[33,9] - ldw[32,5] - ldw[32,9] <= 0
 leMC2: - lDX[33,9] + lDW[32,5] + lDW[32,9] >= 0
 leMC2: lDX[33,9] - lDW[32,5] + lDW[32,9] >= 0
 leMC2: lDX[33,9] + lDW[32,5] - lDW[32,9] >= 0
 leMC2: lDX[33,9] - lDW[32,5] - lDW[32,9] <= 0
 leMC3: ldx[33,13] - ldw[32,1] >= 0
 leMC3: ldx[33,13] - ldw[32,9] >= 0
 leMC3: ldx[33,13] - ldw[32,1] - ldw[32,9] <= 0
 leMC3: - lDX[33,13] + lDW[32,1] + lDW[32,9] >= 0
 leMC3: lDX[33,13] - lDW[32,1] + lDW[32,9] >= 0
 leMC3: lDX[33,13] + lDW[32,1] - lDW[32,9] >= 0
 leMC3: lDX[33,13] - lDW[32,1] - lDW[32,9] <= 0
 leMC0: ldx[33,2] - ldw[32,14] >= 0
 leMC0: ldx[33,2] - ldx[33,14] >= 0
 leMC0: ldx[33,2] - ldx[33,14] - ldw[32,14] <= 0
 leMC0: - lDX[33,2] + lDX[33,14] + lDW[32,14] >= 0
 leMC0: lDX[33,2] + lDX[33,14] - lDW[32,14] >= 0
 leMC0: lDX[33,2] - lDX[33,14] + lDW[32,14] >= 0
 leMC0: lDX[33,2] - lDX[33,14] - lDW[32,14] <= 0
 leMC1: ldx[33,6] - ldw[32,2] = 0
 leMC1: lDX[33,6] - lDW[32,2] = 0
 leMC2: ldx[33,10] - ldw[32,6] >= 0
 leMC2: ldx[33,10] - ldw[32,10] >= 0
 leMC2: ldx[33,10] - ldw[32,6] - ldw[32,10] <= 0
 leMC2: - lDX[33,10] + lDW[32,6] + lDW[32,10] >= 0
 leMC2: lDX[33,10] - lDW[32,6] + lDW[32,10] >= 0
 leMC2: lDX[33,10] + lDW[32,6] - lDW[32,10] >= 0
 leMC2: lDX[33,10] - lDW[32,6] - lDW[32,10] <= 0
 leMC3: ldx[33,14] - ldw[32,2] >= 0
 leMC3: ldx[33,14] - ldw[32,10] >= 0
 leMC3: ldx[33,14] - ldw[32,2] - ldw[32,10] <= 0
 leMC3: - lDX[33,14] + lDW[32,2] + lDW[32,10] >= 0
 leMC3: lDX[33,14] - lDW[32,2] + lDW[32,10] >= 0
 leMC3: lDX[33,14] + lDW[32,2] - lDW[32,10] >= 0
 leMC3: lDX[33,14] - lDW[32,2] - lDW[32,10] <= 0
 leMC0: ldx[33,3] - ldw[32,15] >= 0
 leMC0: ldx[33,3] - ldx[33,15] >= 0
 leMC0: ldx[33,3] - ldx[33,15] - ldw[32,15] <= 0
 leMC0: - lDX[33,3] + lDX[33,15] + lDW[32,15] >= 0
 leMC0: lDX[33,3] + lDX[33,15] - lDW[32,15] >= 0
 leMC0: lDX[33,3] - lDX[33,15] + lDW[32,15] >= 0
 leMC0: lDX[33,3] - lDX[33,15] - lDW[32,15] <= 0
 leMC1: ldx[33,7] - ldw[32,3] = 0
 leMC1: lDX[33,7] - lDW[32,3] = 0
 leMC2: ldx[33,11] - ldw[32,7] >= 0
 leMC2: ldx[33,11] - ldw[32,11] >= 0
 leMC2: ldx[33,11] - ldw[32,7] - ldw[32,11] <= 0
 leMC2: - lDX[33,11] + lDW[32,7] + lDW[32,11] >= 0
 leMC2: lDX[33,11] - lDW[32,7] + lDW[32,11] >= 0
 leMC2: lDX[33,11] + lDW[32,7] - lDW[32,11] >= 0
 leMC2: lDX[33,11] - lDW[32,7] - lDW[32,11] <= 0
 leMC3: ldx[33,15] - ldw[32,3] >= 0
 leMC3: ldx[33,15] - ldw[32,11] >= 0
 leMC3: ldx[33,15] - ldw[32,3] - ldw[32,11] <= 0
 leMC3: - lDX[33,15] + lDW[32,3] + lDW[32,11] >= 0
 leMC3: lDX[33,15] - lDW[32,3] + lDW[32,11] >= 0
 leMC3: lDX[33,15] + lDW[32,3] - lDW[32,11] >= 0
 leMC3: lDX[33,15] - lDW[32,3] - lDW[32,11] <= 0
 lKS: - lstk[17,0] + lLANE[0] >= 0
 lKS: - lstk[18,8] + lLANE[0] >= 0
 lKS: - lstk[19,2] + lLANE[0] >= 0
 lKS: - lstk[20,10] + lLANE[0] >= 0
 lKS: - lstk[21,4] + lLANE[0] >= 0
 lKS: - lstk[22,12] + lLANE[0] >= 0
 lKS: - lstk[23,6] + lLANE[0] >= 0
 lKS: - lstk[24,14] + lLANE[0] >= 0
 lKS: - lstk[25,5] + lLANE[0] >= 0
 lKS: - lstk[26,13] + lLANE[0] >= 0
 lKS: - lstk[27,3] + lLANE[0] >= 0
 lKS: - lstk[28,11] + lLANE[0] >= 0
 lKS: - lstk[29,7] + lLANE[0] >= 0
 lKS: - lstk[30,15] + lLANE[0] >= 0
 lKS: - lstk[31,1] + lLANE[0] >= 0
 lKS: - lstk[32,9] + lLANE[0] >= 0
 lKS: lstk[17,0] + lstk[18,8] + lstk[19,2] + lstk[20,10] + lstk[21,4]
   + lstk[22,12] + lstk[23,6] + lstk[24,14] + lstk[25,5] + lstk[26,13]
   + lstk[27,3] + lstk[28,11] + lstk[29,7] + lstk[30,15] + lstk[31,1]
   + lstk[32,9] - 16 lLANE[0] >= -3
 lKS: - lstk[17,1] + lLANE[1] >= 0
 lKS: - lstk[18,9] + lLANE[1] >= 0
 lKS: - lstk[19,0] + lLANE[1] >= 0
 lKS: - lstk[20,8] + lLANE[1] >= 0
 lKS: - lstk[21,2] + lLANE[1] >= 0
 lKS: - lstk[22,10] + lLANE[1] >= 0
 lKS: - lstk[23,4] + lLANE[1] >= 0
 lKS: - lstk[24,12] + lLANE[1] >= 0
 lKS: - lstk[25,6] + lLANE[1] >= 0
 lKS: - lstk[26,14] + lLANE[1] >= 0
 lKS: - lstk[27,5] + lLANE[1] >= 0
 lKS: - lstk[28,13] + lLANE[1] >= 0
 lKS: - lstk[29,3] + lLANE[1] >= 0
 lKS: - lstk[30,11] + lLANE[1] >= 0
 lKS: - lstk[31,7] + lLANE[1] >= 0
 lKS: - lstk[32,15] + lLANE[1] >= 0
 lKS: lstk[17,1] + lstk[18,9] + lstk[19,0] + lstk[20,8] + lstk[21,2]
   + lstk[22,10] + lstk[23,4] + lstk[24,12] + lstk[25,6] + lstk[26,14]
   + lstk[27,5] + lstk[28,13] + lstk[29,3] + lstk[30,11] + lstk[31,7]
   + lstk[32,15] - 16 lLANE[1] >= -3
 lKS: - lstk[17,2] + lLANE[2] >= 0
 lKS: - lstk[18,10] + lLANE[2] >= 0
 lKS: - lstk[19,4] + lLANE[2] >= 0
 lKS: - lstk[20,12] + lLANE[2] >= 0
 lKS: - lstk[21,6] + lLANE[2] >= 0
 lKS: - lstk[22,14] + lLANE[2] >= 0
 lKS: - lstk[23,5] + lLANE[2] >= 0
 lKS: - lstk[24,13] + lLANE[2] >= 0
 lKS: - lstk[25,3] + lLANE[2] >= 0
 lKS: - lstk[26,11] + lLANE[2] >= 0
 lKS: - lstk[27,7] + lLANE[2] >= 0
 lKS: - lstk[28,15] + lLANE[2] >= 0
 lKS: - lstk[29,1] + lLANE[2] >= 0
 lKS: - lstk[30,9] + lLANE[2] >= 0
 lKS: - lstk[31,0] + lLANE[2] >= 0
 lKS: - lstk[32,8] + lLANE[2] >= 0
 lKS: lstk[17,2] + lstk[18,10] + lstk[19,4] + lstk[20,12] + lstk[21,6]
   + lstk[22,14] + lstk[23,5] + lstk[24,13] + lstk[25,3] + lstk[26,11]
   + lstk[27,7] + lstk[28,15] + lstk[29,1] + lstk[30,9] + lstk[31,0]
   + lstk[32,8] - 16 lLANE[2] >= -3
 lKS: - lstk[17,3] + lLANE[3] >= 0
 lKS: - lstk[18,11] + lLANE[3] >= 0
 lKS: - lstk[19,7] + lLANE[3] >= 0
 lKS: - lstk[20,15] + lLANE[3] >= 0
 lKS: - lstk[21,1] + lLANE[3] >= 0
 lKS: - lstk[22,9] + lLANE[3] >= 0
 lKS: - lstk[23,0] + lLANE[3] >= 0
 lKS: - lstk[24,8] + lLANE[3] >= 0
 lKS: - lstk[25,2] + lLANE[3] >= 0
 lKS: - lstk[26,10] + lLANE[3] >= 0
 lKS: - lstk[27,4] + lLANE[3] >= 0
 lKS: - lstk[28,12] + lLANE[3] >= 0
 lKS: - lstk[29,6] + lLANE[3] >= 0
 lKS: - lstk[30,14] + lLANE[3] >= 0
 lKS: - lstk[31,5] + lLANE[3] >= 0
 lKS: - lstk[32,13] + lLANE[3] >= 0
 lKS: lstk[17,3] + lstk[18,11] + lstk[19,7] + lstk[20,15] + lstk[21,1]
   + lstk[22,9] + lstk[23,0] + lstk[24,8] + lstk[25,2] + lstk[26,10]
   + lstk[27,4] + lstk[28,12] + lstk[29,6] + lstk[30,14] + lstk[31,5]
   + lstk[32,13] - 16 lLANE[3] >= -3
 lKS: - lstk[17,4] + lLANE[4] >= 0
 lKS: - lstk[18,12] + lLANE[4] >= 0
 lKS: - lstk[19,6] + lLANE[4] >= 0
 lKS: - lstk[20,14] + lLANE[4] >= 0
 lKS: - lstk[21,5] + lLANE[4] >= 0
 lKS: - lstk[22,13] + lLANE[4] >= 0
 lKS: - lstk[23,3] + lLANE[4] >= 0
 lKS: - lstk[24,11] + lLANE[4] >= 0
 lKS: - lstk[25,7] + lLANE[4] >= 0
 lKS: - lstk[26,15] + lLANE[4] >= 0
 lKS: - lstk[27,1] + lLANE[4] >= 0
 lKS: - lstk[28,9] + lLANE[4] >= 0
 lKS: - lstk[29,0] + lLANE[4] >= 0
 lKS: - lstk[30,8] + lLANE[4] >= 0
 lKS: - lstk[31,2] + lLANE[4] >= 0
 lKS: - lstk[32,10] + lLANE[4] >= 0
 lKS: lstk[17,4] + lstk[18,12] + lstk[19,6] + lstk[20,14] + lstk[21,5]
   + lstk[22,13] + lstk[23,3] + lstk[24,11] + lstk[25,7] + lstk[26,15]
   + lstk[27,1] + lstk[28,9] + lstk[29,0] + lstk[30,8] + lstk[31,2]
   + lstk[32,10] - 16 lLANE[4] >= -3
 lKS: - lstk[17,5] + lLANE[5] >= 0
 lKS: - lstk[18,13] + lLANE[5] >= 0
 lKS: - lstk[19,3] + lLANE[5] >= 0
 lKS: - lstk[20,11] + lLANE[5] >= 0
 lKS: - lstk[21,7] + lLANE[5] >= 0
 lKS: - lstk[22,15] + lLANE[5] >= 0
 lKS: - lstk[23,1] + lLANE[5] >= 0
 lKS: - lstk[24,9] + lLANE[5] >= 0
 lKS: - lstk[25,0] + lLANE[5] >= 0
 lKS: - lstk[26,8] + lLANE[5] >= 0
 lKS: - lstk[27,2] + lLANE[5] >= 0
 lKS: - lstk[28,10] + lLANE[5] >= 0
 lKS: - lstk[29,4] + lLANE[5] >= 0
 lKS: - lstk[30,12] + lLANE[5] >= 0
 lKS: - lstk[31,6] + lLANE[5] >= 0
 lKS: - lstk[32,14] + lLANE[5] >= 0
 lKS: lstk[17,5] + lstk[18,13] + lstk[19,3] + lstk[20,11] + lstk[21,7]
   + lstk[22,15] + lstk[23,1] + lstk[24,9] + lstk[25,0] + lstk[26,8]
   + lstk[27,2] + lstk[28,10] + lstk[29,4] + lstk[30,12] + lstk[31,6]
   + lstk[32,14] - 16 lLANE[5] >= -3
 lKS: - lstk[17,6] + lLANE[6] >= 0
 lKS: - lstk[18,14] + lLANE[6] >= 0
 lKS: - lstk[19,5] + lLANE[6] >= 0
 lKS: - lstk[20,13] + lLANE[6] >= 0
 lKS: - lstk[21,3] + lLANE[6] >= 0
 lKS: - lstk[22,11] + lLANE[6] >= 0
 lKS: - lstk[23,7] + lLANE[6] >= 0
 lKS: - lstk[24,15] + lLANE[6] >= 0
 lKS: - lstk[25,1] + lLANE[6] >= 0
 lKS: - lstk[26,9] + lLANE[6] >= 0
 lKS: - lstk[27,0] + lLANE[6] >= 0
 lKS: - lstk[28,8] + lLANE[6] >= 0
 lKS: - lstk[29,2] + lLANE[6] >= 0
 lKS: - lstk[30,10] + lLANE[6] >= 0
 lKS: - lstk[31,4] + lLANE[6] >= 0
 lKS: - lstk[32,12] + lLANE[6] >= 0
 lKS: lstk[17,6] + lstk[18,14] + lstk[19,5] + lstk[20,13] + lstk[21,3]
   + lstk[22,11] + lstk[23,7] + lstk[24,15] + lstk[25,1] + lstk[26,9]
   + lstk[27,0] + lstk[28,8] + lstk[29,2] + lstk[30,10] + lstk[31,4]
   + lstk[32,12] - 16 lLANE[6] >= -3
 lKS: - lstk[17,7] + lLANE[7] >= 0
 lKS: - lstk[18,15] + lLANE[7] >= 0
 lKS: - lstk[19,1] + lLANE[7] >= 0
 lKS: - lstk[20,9] + lLANE[7] >= 0
 lKS: - lstk[21,0] + lLANE[7] >= 0
 lKS: - lstk[22,8] + lLANE[7] >= 0
 lKS: - lstk[23,2] + lLANE[7] >= 0
 lKS: - lstk[24,10] + lLANE[7] >= 0
 lKS: - lstk[25,4] + lLANE[7] >= 0
 lKS: - lstk[26,12] + lLANE[7] >= 0
 lKS: - lstk[27,6] + lLANE[7] >= 0
 lKS: - lstk[28,14] + lLANE[7] >= 0
 lKS: - lstk[29,5] + lLANE[7] >= 0
 lKS: - lstk[30,13] + lLANE[7] >= 0
 lKS: - lstk[31,3] + lLANE[7] >= 0
 lKS: - lstk[32,11] + lLANE[7] >= 0
 lKS: lstk[17,7] + lstk[18,15] + lstk[19,1] + lstk[20,9] + lstk[21,0]
   + lstk[22,8] + lstk[23,2] + lstk[24,10] + lstk[25,4] + lstk[26,12]
   + lstk[27,6] + lstk[28,14] + lstk[29,5] + lstk[30,13] + lstk[31,3]
   + lstk[32,11] - 16 lLANE[7] >= -3
 lKS: - lstk[17,8] + lLANE[8] >= 0
 lKS: - lstk[18,2] + lLANE[8] >= 0
 lKS: - lstk[19,10] + lLANE[8] >= 0
 lKS: - lstk[20,4] + lLANE[8] >= 0
 lKS: - lstk[21,12] + lLANE[8] >= 0
 lKS: - lstk[22,6] + lLANE[8] >= 0
 lKS: - lstk[23,14] + lLANE[8] >= 0
 lKS: - lstk[24,5] + lLANE[8] >= 0
 lKS: - lstk[25,13] + lLANE[8] >= 0
 lKS: - lstk[26,3] + lLANE[8] >= 0
 lKS: - lstk[27,11] + lLANE[8] >= 0
 lKS: - lstk[28,7] + lLANE[8] >= 0
 lKS: - lstk[29,15] + lLANE[8] >= 0
 lKS: - lstk[30,1] + lLANE[8] >= 0
 lKS: - lstk[31,9] + lLANE[8] >= 0
 lKS: - lstk[32,0] + lLANE[8] >= 0
 lKS: lstk[17,8] + lstk[18,2] + lstk[19,10] + lstk[20,4] + lstk[21,12]
   + lstk[22,6] + lstk[23,14] + lstk[24,5] + lstk[25,13] + lstk[26,3]
   + lstk[27,11] + lstk[28,7] + lstk[29,15] + lstk[30,1] + lstk[31,9]
   + lstk[32,0] - 16 lLANE[8] >= -3
 lKS: - lstk[17,9] + lLANE[9] >= 0
 lKS: - lstk[18,0] + lLANE[9] >= 0
 lKS: - lstk[19,8] + lLANE[9] >= 0
 lKS: - lstk[20,2] + lLANE[9] >= 0
 lKS: - lstk[21,10] + lLANE[9] >= 0
 lKS: - lstk[22,4] + lLANE[9] >= 0
 lKS: - lstk[23,12] + lLANE[9] >= 0
 lKS: - lstk[24,6] + lLANE[9] >= 0
 lKS: - lstk[25,14] + lLANE[9] >= 0
 lKS: - lstk[26,5] + lLANE[9] >= 0
 lKS: - lstk[27,13] + lLANE[9] >= 0
 lKS: - lstk[28,3] + lLANE[9] >= 0
 lKS: - lstk[29,11] + lLANE[9] >= 0
 lKS: - lstk[30,7] + lLANE[9] >= 0
 lKS: - lstk[31,15] + lLANE[9] >= 0
 lKS: - lstk[32,1] + lLANE[9] >= 0
 lKS: lstk[17,9] + lstk[18,0] + lstk[19,8] + lstk[20,2] + lstk[21,10]
   + lstk[22,4] + lstk[23,12] + lstk[24,6] + lstk[25,14] + lstk[26,5]
   + lstk[27,13] + lstk[28,3] + lstk[29,11] + lstk[30,7] + lstk[31,15]
   + lstk[32,1] - 16 lLANE[9] >= -3
 lKS: - lstk[17,10] + lLANE[10] >= 0
 lKS: - lstk[18,4] + lLANE[10] >= 0
 lKS: - lstk[19,12] + lLANE[10] >= 0
 lKS: - lstk[20,6] + lLANE[10] >= 0
 lKS: - lstk[21,14] + lLANE[10] >= 0
 lKS: - lstk[22,5] + lLANE[10] >= 0
 lKS: - lstk[23,13] + lLANE[10] >= 0
 lKS: - lstk[24,3] + lLANE[10] >= 0
 lKS: - lstk[25,11] + lLANE[10] >= 0
 lKS: - lstk[26,7] + lLANE[10] >= 0
 lKS: - lstk[27,15] + lLANE[10] >= 0
 lKS: - lstk[28,1] + lLANE[10] >= 0
 lKS: - lstk[29,9] + lLANE[10] >= 0
 lKS: - lstk[30,0] + lLANE[10] >= 0
 lKS: - lstk[31,8] + lLANE[10] >= 0
 lKS: - lstk[32,2] + lLANE[10] >= 0
 lKS: lstk[17,10] + lstk[18,4] + lstk[19,12] + lstk[20,6] + lstk[21,14]
   + lstk[22,5] + lstk[23,13] + lstk[24,3] + lstk[25,11] + lstk[26,7]
   + lstk[27,15] + lstk[28,1] + lstk[29,9] + lstk[30,0] + lstk[31,8]
   + lstk[32,2] - 16 lLANE[10] >= -3
 lKS: - lstk[17,11] + lLANE[11] >= 0
 lKS: - lstk[18,7] + lLANE[11] >= 0
 lKS: - lstk[19,15] + lLANE[11] >= 0
 lKS: - lstk[20,1] + lLANE[11] >= 0
 lKS: - lstk[21,9] + lLANE[11] >= 0
 lKS: - lstk[22,0] + lLANE[11] >= 0
 lKS: - lstk[23,8] + lLANE[11] >= 0
 lKS: - lstk[24,2] + lLANE[11] >= 0
 lKS: - lstk[25,10] + lLANE[11] >= 0
 lKS: - lstk[26,4] + lLANE[11] >= 0
 lKS: - lstk[27,12] + lLANE[11] >= 0
 lKS: - lstk[28,6] + lLANE[11] >= 0
 lKS: - lstk[29,14] + lLANE[11] >= 0
 lKS: - lstk[30,5] + lLANE[11] >= 0
 lKS: - lstk[31,13] + lLANE[11] >= 0
 lKS: - lstk[32,3] + lLANE[11] >= 0
 lKS: lstk[17,11] + lstk[18,7] + lstk[19,15] + lstk[20,1] + lstk[21,9]
   + lstk[22,0] + lstk[23,8] + lstk[24,2] + lstk[25,10] + lstk[26,4]
   + lstk[27,12] + lstk[28,6] + lstk[29,14] + lstk[30,5] + lstk[31,13]
   + lstk[32,3] - 16 lLANE[11] >= -3
 lKS: - lstk[17,12] + lLANE[12] >= 0
 lKS: - lstk[18,6] + lLANE[12] >= 0
 lKS: - lstk[19,14] + lLANE[12] >= 0
 lKS: - lstk[20,5] + lLANE[12] >= 0
 lKS: - lstk[21,13] + lLANE[12] >= 0
 lKS: - lstk[22,3] + lLANE[12] >= 0
 lKS: - lstk[23,11] + lLANE[12] >= 0
 lKS: - lstk[24,7] + lLANE[12] >= 0
 lKS: - lstk[25,15] + lLANE[12] >= 0
 lKS: - lstk[26,1] + lLANE[12] >= 0
 lKS: - lstk[27,9] + lLANE[12] >= 0
 lKS: - lstk[28,0] + lLANE[12] >= 0
 lKS: - lstk[29,8] + lLANE[12] >= 0
 lKS: - lstk[30,2] + lLANE[12] >= 0
 lKS: - lstk[31,10] + lLANE[12] >= 0
 lKS: - lstk[32,4] + lLANE[12] >= 0
 lKS: lstk[17,12] + lstk[18,6] + lstk[19,14] + lstk[20,5] + lstk[21,13]
   + lstk[22,3] + lstk[23,11] + lstk[24,7] + lstk[25,15] + lstk[26,1]
   + lstk[27,9] + lstk[28,0] + lstk[29,8] + lstk[30,2] + lstk[31,10]
   + lstk[32,4] - 16 lLANE[12] >= -3
 lKS: - lstk[17,13] + lLANE[13] >= 0
 lKS: - lstk[18,3] + lLANE[13] >= 0
 lKS: - lstk[19,11] + lLANE[13] >= 0
 lKS: - lstk[20,7] + lLANE[13] >= 0
 lKS: - lstk[21,15] + lLANE[13] >= 0
 lKS: - lstk[22,1] + lLANE[13] >= 0
 lKS: - lstk[23,9] + lLANE[13] >= 0
 lKS: - lstk[24,0] + lLANE[13] >= 0
 lKS: - lstk[25,8] + lLANE[13] >= 0
 lKS: - lstk[26,2] + lLANE[13] >= 0
 lKS: - lstk[27,10] + lLANE[13] >= 0
 lKS: - lstk[28,4] + lLANE[13] >= 0
 lKS: - lstk[29,12] + lLANE[13] >= 0
 lKS: - lstk[30,6] + lLANE[13] >= 0
 lKS: - lstk[31,14] + lLANE[13] >= 0
 lKS: - lstk[32,5] + lLANE[13] >= 0
 lKS: lstk[17,13] + lstk[18,3] + lstk[19,11] + lstk[20,7] + lstk[21,15]
   + lstk[22,1] + lstk[23,9] + lstk[24,0] + lstk[25,8] + lstk[26,2]
   + lstk[27,10] + lstk[28,4] + lstk[29,12] + lstk[30,6] + lstk[31,14]
   + lstk[32,5] - 16 lLANE[13] >= -3
 lKS: - lstk[17,14] + lLANE[14] >= 0
 lKS: - lstk[18,5] + lLANE[14] >= 0
 lKS: - lstk[19,13] + lLANE[14] >= 0
 lKS: - lstk[20,3] + lLANE[14] >= 0
 lKS: - lstk[21,11] + lLANE[14] >= 0
 lKS: - lstk[22,7] + lLANE[14] >= 0
 lKS: - lstk[23,15] + lLANE[14] >= 0
 lKS: - lstk[24,1] + lLANE[14] >= 0
 lKS: - lstk[25,9] + lLANE[14] >= 0
 lKS: - lstk[26,0] + lLANE[14] >= 0
 lKS: - lstk[27,8] + lLANE[14] >= 0
 lKS: - lstk[28,2] + lLANE[14] >= 0
 lKS: - lstk[29,10] + lLANE[14] >= 0
 lKS: - lstk[30,4] + lLANE[14] >= 0
 lKS: - lstk[31,12] + lLANE[14] >= 0
 lKS: - lstk[32,6] + lLANE[14] >= 0
 lKS: lstk[17,14] + lstk[18,5] + lstk[19,13] + lstk[20,3] + lstk[21,11]
   + lstk[22,7] + lstk[23,15] + lstk[24,1] + lstk[25,9] + lstk[26,0]
   + lstk[27,8] + lstk[28,2] + lstk[29,10] + lstk[30,4] + lstk[31,12]
   + lstk[32,6] - 16 lLANE[14] >= -3
 lKS: - lstk[17,15] + lLANE[15] >= 0
 lKS: - lstk[18,1] + lLANE[15] >= 0
 lKS: - lstk[19,9] + lLANE[15] >= 0
 lKS: - lstk[20,0] + lLANE[15] >= 0
 lKS: - lstk[21,8] + lLANE[15] >= 0
 lKS: - lstk[22,2] + lLANE[15] >= 0
 lKS: - lstk[23,10] + lLANE[15] >= 0
 lKS: - lstk[24,4] + lLANE[15] >= 0
 lKS: - lstk[25,12] + lLANE[15] >= 0
 lKS: - lstk[26,6] + lLANE[15] >= 0
 lKS: - lstk[27,14] + lLANE[15] >= 0
 lKS: - lstk[28,5] + lLANE[15] >= 0
 lKS: - lstk[29,13] + lLANE[15] >= 0
 lKS: - lstk[30,3] + lLANE[15] >= 0
 lKS: - lstk[31,11] + lLANE[15] >= 0
 lKS: - lstk[32,7] + lLANE[15] >= 0
 lKS: lstk[17,15] + lstk[18,1] + lstk[19,9] + lstk[20,0] + lstk[21,8]
   + lstk[22,2] + lstk[23,10] + lstk[24,4] + lstk[25,12] + lstk[26,6]
   + lstk[27,14] + lstk[28,5] + lstk[29,13] + lstk[30,3] + lstk[31,11]
   + lstk[32,7] - 16 lLANE[15] >= -3
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
 uGstkL: uGstk[4,8] = 0
 uGstkL: uGstk[4,9] = 0
 uGstkL: uGstk[4,10] = 0
 uGstkL: uGstk[4,11] = 0
 uGstkL: uGstk[4,12] = 0
 uGstkL: uGstk[4,13] = 0
 uGstkL: uGstk[4,14] = 0
 uGstkL: uGstk[4,15] = 0
 uGstkL: uGstk[5,8] = 0
 uGstkL: uGstk[5,9] = 0
 uGstkL: uGstk[5,10] = 0
 uGstkL: uGstk[5,11] = 0
 uGstkL: uGstk[5,12] = 0
 uGstkL: uGstk[5,13] = 0
 uGstkL: uGstk[5,14] = 0
 uGstkL: uGstk[5,15] = 0
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
 ugstkLa: ugstka[4,8] = 0
 ugstkLa: ugstka[4,9] = 0
 ugstkLa: ugstka[4,10] = 0
 ugstkLa: ugstka[4,11] = 0
 ugstkLa: ugstka[4,12] = 0
 ugstkLa: ugstka[4,13] = 0
 ugstkLa: ugstka[4,14] = 0
 ugstkLa: ugstka[4,15] = 0
 ugstkLa: ugstka[5,8] = 0
 ugstkLa: ugstka[5,9] = 0
 ugstkLa: ugstka[5,10] = 0
 ugstkLa: ugstka[5,11] = 0
 ugstkLa: ugstka[5,12] = 0
 ugstkLa: ugstka[5,13] = 0
 ugstkLa: ugstka[5,14] = 0
 ugstkLa: ugstka[5,15] = 0
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
 geqk0: ugEQK[4,0] - uGstk[4,0] = 0
 geqk0A: ugeqka[4,0] - ugstka[4,0] = 0
 geqk1: ugEQK[4,4] - uGstk[4,0] = 0
 geqk1A: ugeqka[4,4] - ugstka[4,0] = 0
 geqk2: ugEQK[4,8] - uGstk[4,7] = 0
 geqk2A: ugeqka[4,8] - ugstka[4,7] = 0
 geqk3: ugEQK[4,12] - uGstk[4,0] = 0
 geqk3A: ugeqka[4,12] - ugstka[4,0] = 0
 geqk0: ugEQK[4,1] - uGstk[4,1] = 0
 geqk0A: ugeqka[4,1] - ugstka[4,1] = 0
 geqk1: ugEQK[4,5] - uGstk[4,1] = 0
 geqk1A: ugeqka[4,5] - ugstka[4,1] = 0
 geqk2: ugEQK[4,9] - uGstk[4,4] = 0
 geqk2A: ugeqka[4,9] - ugstka[4,4] = 0
 geqk3: ugEQK[4,13] - uGstk[4,1] = 0
 geqk3A: ugeqka[4,13] - ugstka[4,1] = 0
 geqk0: ugEQK[4,2] - uGstk[4,2] = 0
 geqk0A: ugeqka[4,2] - ugstka[4,2] = 0
 geqk1: ugEQK[4,6] - uGstk[4,2] = 0
 geqk1A: ugeqka[4,6] - ugstka[4,2] = 0
 geqk2: ugEQK[4,10] - uGstk[4,5] = 0
 geqk2A: ugeqka[4,10] - ugstka[4,5] = 0
 geqk3: ugEQK[4,14] - uGstk[4,2] = 0
 geqk3A: ugeqka[4,14] - ugstka[4,2] = 0
 geqk0: ugEQK[4,3] - uGstk[4,3] = 0
 geqk0A: ugeqka[4,3] - ugstka[4,3] = 0
 geqk1: ugEQK[4,7] - uGstk[4,3] = 0
 geqk1A: ugeqka[4,7] - ugstka[4,3] = 0
 geqk2: ugEQK[4,11] - uGstk[4,6] = 0
 geqk2A: ugeqka[4,11] - ugstka[4,6] = 0
 geqk3: ugEQK[4,15] - uGstk[4,3] = 0
 geqk3A: ugeqka[4,15] - ugstka[4,3] = 0
 geqk0: ugEQK[5,0] - uGstk[5,0] = 0
 geqk0A: ugeqka[5,0] - ugstka[5,0] = 0
 geqk1: ugEQK[5,4] - uGstk[5,0] = 0
 geqk1A: ugeqka[5,4] - ugstka[5,0] = 0
 geqk2: ugEQK[5,8] - uGstk[5,7] = 0
 geqk2A: ugeqka[5,8] - ugstka[5,7] = 0
 geqk3: ugEQK[5,12] - uGstk[5,0] = 0
 geqk3A: ugeqka[5,12] - ugstka[5,0] = 0
 geqk0: ugEQK[5,1] - uGstk[5,1] = 0
 geqk0A: ugeqka[5,1] - ugstka[5,1] = 0
 geqk1: ugEQK[5,5] - uGstk[5,1] = 0
 geqk1A: ugeqka[5,5] - ugstka[5,1] = 0
 geqk2: ugEQK[5,9] - uGstk[5,4] = 0
 geqk2A: ugeqka[5,9] - ugstka[5,4] = 0
 geqk3: ugEQK[5,13] - uGstk[5,1] = 0
 geqk3A: ugeqka[5,13] - ugstka[5,1] = 0
 geqk0: ugEQK[5,2] - uGstk[5,2] = 0
 geqk0A: ugeqka[5,2] - ugstka[5,2] = 0
 geqk1: ugEQK[5,6] - uGstk[5,2] = 0
 geqk1A: ugeqka[5,6] - ugstka[5,2] = 0
 geqk2: ugEQK[5,10] - uGstk[5,5] = 0
 geqk2A: ugeqka[5,10] - ugstka[5,5] = 0
 geqk3: ugEQK[5,14] - uGstk[5,2] = 0
 geqk3A: ugeqka[5,14] - ugstka[5,2] = 0
 geqk0: ugEQK[5,3] - uGstk[5,3] = 0
 geqk0A: ugeqka[5,3] - ugstka[5,3] = 0
 geqk1: ugEQK[5,7] - uGstk[5,3] = 0
 geqk1A: ugeqka[5,7] - ugstka[5,3] = 0
 geqk2: ugEQK[5,11] - uGstk[5,6] = 0
 geqk2A: ugeqka[5,11] - ugstka[5,6] = 0
 geqk3: ugEQK[5,15] - uGstk[5,3] = 0
 geqk3A: ugeqka[5,15] - ugstka[5,3] = 0
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
 udetSC+ek: - uDetX[5,0] + uDetY[5,0] <= 0
 udetSC+ek: - ugEQK[4,0] + uDetY[5,0] <= 0
 udetSC+ek: - ugEQK[4,0] - uDetX[5,0] + uDetY[5,0] >= -1
 udetSC+ekA: - udetXa[5,0] + udetYa[5,0] <= 0
 udetSC+ekA: - ugeqka[4,0] + udetYa[5,0] <= 0
 udetSC+ekA: - ugeqka[4,0] - udetXa[5,0] + udetYa[5,0] >= -1
 ufrSC: uDY[5,0] - udy[5,0] - uFrSC[5,0] >= 0
 ufrSC: uDetY[5,0] - uFrSC[5,0] >= 0
 ufrSC: - uDY[5,0] + udy[5,0] - uDetY[5,0] + uFrSC[5,0] >= -1
 ufrSCa: uDY[5,0] - udy[5,0] - ufrSCa[5,0] >= 0
 ufrSCa: udetYa[5,0] - ufrSCa[5,0] >= 0
 ufrSCa: - uDY[5,0] + udy[5,0] - udetYa[5,0] + ufrSCa[5,0] >= -1
 udetSC+ek: - uDetX[5,1] + uDetY[5,1] <= 0
 udetSC+ek: - ugEQK[4,1] + uDetY[5,1] <= 0
 udetSC+ek: - ugEQK[4,1] - uDetX[5,1] + uDetY[5,1] >= -1
 udetSC+ekA: - udetXa[5,1] + udetYa[5,1] <= 0
 udetSC+ekA: - ugeqka[4,1] + udetYa[5,1] <= 0
 udetSC+ekA: - ugeqka[4,1] - udetXa[5,1] + udetYa[5,1] >= -1
 ufrSC: uDY[5,1] - udy[5,1] - uFrSC[5,1] >= 0
 ufrSC: uDetY[5,1] - uFrSC[5,1] >= 0
 ufrSC: - uDY[5,1] + udy[5,1] - uDetY[5,1] + uFrSC[5,1] >= -1
 ufrSCa: uDY[5,1] - udy[5,1] - ufrSCa[5,1] >= 0
 ufrSCa: udetYa[5,1] - ufrSCa[5,1] >= 0
 ufrSCa: - uDY[5,1] + udy[5,1] - udetYa[5,1] + ufrSCa[5,1] >= -1
 udetSC+ek: - uDetX[5,2] + uDetY[5,2] <= 0
 udetSC+ek: - ugEQK[4,2] + uDetY[5,2] <= 0
 udetSC+ek: - ugEQK[4,2] - uDetX[5,2] + uDetY[5,2] >= -1
 udetSC+ekA: - udetXa[5,2] + udetYa[5,2] <= 0
 udetSC+ekA: - ugeqka[4,2] + udetYa[5,2] <= 0
 udetSC+ekA: - ugeqka[4,2] - udetXa[5,2] + udetYa[5,2] >= -1
 ufrSC: uDY[5,2] - udy[5,2] - uFrSC[5,2] >= 0
 ufrSC: uDetY[5,2] - uFrSC[5,2] >= 0
 ufrSC: - uDY[5,2] + udy[5,2] - uDetY[5,2] + uFrSC[5,2] >= -1
 ufrSCa: uDY[5,2] - udy[5,2] - ufrSCa[5,2] >= 0
 ufrSCa: udetYa[5,2] - ufrSCa[5,2] >= 0
 ufrSCa: - uDY[5,2] + udy[5,2] - udetYa[5,2] + ufrSCa[5,2] >= -1
 udetSC+ek: - uDetX[5,3] + uDetY[5,3] <= 0
 udetSC+ek: - ugEQK[4,3] + uDetY[5,3] <= 0
 udetSC+ek: - ugEQK[4,3] - uDetX[5,3] + uDetY[5,3] >= -1
 udetSC+ekA: - udetXa[5,3] + udetYa[5,3] <= 0
 udetSC+ekA: - ugeqka[4,3] + udetYa[5,3] <= 0
 udetSC+ekA: - ugeqka[4,3] - udetXa[5,3] + udetYa[5,3] >= -1
 ufrSC: uDY[5,3] - udy[5,3] - uFrSC[5,3] >= 0
 ufrSC: uDetY[5,3] - uFrSC[5,3] >= 0
 ufrSC: - uDY[5,3] + udy[5,3] - uDetY[5,3] + uFrSC[5,3] >= -1
 ufrSCa: uDY[5,3] - udy[5,3] - ufrSCa[5,3] >= 0
 ufrSCa: udetYa[5,3] - ufrSCa[5,3] >= 0
 ufrSCa: - uDY[5,3] + udy[5,3] - udetYa[5,3] + ufrSCa[5,3] >= -1
 udetSC+ek: - uDetX[5,4] + uDetY[5,4] <= 0
 udetSC+ek: - ugEQK[4,4] + uDetY[5,4] <= 0
 udetSC+ek: - ugEQK[4,4] - uDetX[5,4] + uDetY[5,4] >= -1
 udetSC+ekA: - udetXa[5,4] + udetYa[5,4] <= 0
 udetSC+ekA: - ugeqka[4,4] + udetYa[5,4] <= 0
 udetSC+ekA: - ugeqka[4,4] - udetXa[5,4] + udetYa[5,4] >= -1
 ufrSC: uDY[5,4] - udy[5,4] - uFrSC[5,4] >= 0
 ufrSC: uDetY[5,4] - uFrSC[5,4] >= 0
 ufrSC: - uDY[5,4] + udy[5,4] - uDetY[5,4] + uFrSC[5,4] >= -1
 ufrSCa: uDY[5,4] - udy[5,4] - ufrSCa[5,4] >= 0
 ufrSCa: udetYa[5,4] - ufrSCa[5,4] >= 0
 ufrSCa: - uDY[5,4] + udy[5,4] - udetYa[5,4] + ufrSCa[5,4] >= -1
 udetSC+ek: - uDetX[5,5] + uDetY[5,5] <= 0
 udetSC+ek: - ugEQK[4,5] + uDetY[5,5] <= 0
 udetSC+ek: - ugEQK[4,5] - uDetX[5,5] + uDetY[5,5] >= -1
 udetSC+ekA: - udetXa[5,5] + udetYa[5,5] <= 0
 udetSC+ekA: - ugeqka[4,5] + udetYa[5,5] <= 0
 udetSC+ekA: - ugeqka[4,5] - udetXa[5,5] + udetYa[5,5] >= -1
 ufrSC: uDY[5,5] - udy[5,5] - uFrSC[5,5] >= 0
 ufrSC: uDetY[5,5] - uFrSC[5,5] >= 0
 ufrSC: - uDY[5,5] + udy[5,5] - uDetY[5,5] + uFrSC[5,5] >= -1
 ufrSCa: uDY[5,5] - udy[5,5] - ufrSCa[5,5] >= 0
 ufrSCa: udetYa[5,5] - ufrSCa[5,5] >= 0
 ufrSCa: - uDY[5,5] + udy[5,5] - udetYa[5,5] + ufrSCa[5,5] >= -1
 udetSC+ek: - uDetX[5,6] + uDetY[5,6] <= 0
 udetSC+ek: - ugEQK[4,6] + uDetY[5,6] <= 0
 udetSC+ek: - ugEQK[4,6] - uDetX[5,6] + uDetY[5,6] >= -1
 udetSC+ekA: - udetXa[5,6] + udetYa[5,6] <= 0
 udetSC+ekA: - ugeqka[4,6] + udetYa[5,6] <= 0
 udetSC+ekA: - ugeqka[4,6] - udetXa[5,6] + udetYa[5,6] >= -1
 ufrSC: uDY[5,6] - udy[5,6] - uFrSC[5,6] >= 0
 ufrSC: uDetY[5,6] - uFrSC[5,6] >= 0
 ufrSC: - uDY[5,6] + udy[5,6] - uDetY[5,6] + uFrSC[5,6] >= -1
 ufrSCa: uDY[5,6] - udy[5,6] - ufrSCa[5,6] >= 0
 ufrSCa: udetYa[5,6] - ufrSCa[5,6] >= 0
 ufrSCa: - uDY[5,6] + udy[5,6] - udetYa[5,6] + ufrSCa[5,6] >= -1
 udetSC+ek: - uDetX[5,7] + uDetY[5,7] <= 0
 udetSC+ek: - ugEQK[4,7] + uDetY[5,7] <= 0
 udetSC+ek: - ugEQK[4,7] - uDetX[5,7] + uDetY[5,7] >= -1
 udetSC+ekA: - udetXa[5,7] + udetYa[5,7] <= 0
 udetSC+ekA: - ugeqka[4,7] + udetYa[5,7] <= 0
 udetSC+ekA: - ugeqka[4,7] - udetXa[5,7] + udetYa[5,7] >= -1
 ufrSC: uDY[5,7] - udy[5,7] - uFrSC[5,7] >= 0
 ufrSC: uDetY[5,7] - uFrSC[5,7] >= 0
 ufrSC: - uDY[5,7] + udy[5,7] - uDetY[5,7] + uFrSC[5,7] >= -1
 ufrSCa: uDY[5,7] - udy[5,7] - ufrSCa[5,7] >= 0
 ufrSCa: udetYa[5,7] - ufrSCa[5,7] >= 0
 ufrSCa: - uDY[5,7] + udy[5,7] - udetYa[5,7] + ufrSCa[5,7] >= -1
 udetSC+ek: - uDetX[5,8] + uDetY[5,8] <= 0
 udetSC+ek: - ugEQK[4,8] + uDetY[5,8] <= 0
 udetSC+ek: - ugEQK[4,8] - uDetX[5,8] + uDetY[5,8] >= -1
 udetSC+ekA: - udetXa[5,8] + udetYa[5,8] <= 0
 udetSC+ekA: - ugeqka[4,8] + udetYa[5,8] <= 0
 udetSC+ekA: - ugeqka[4,8] - udetXa[5,8] + udetYa[5,8] >= -1
 ufrSC: uDY[5,8] - udy[5,8] - uFrSC[5,8] >= 0
 ufrSC: uDetY[5,8] - uFrSC[5,8] >= 0
 ufrSC: - uDY[5,8] + udy[5,8] - uDetY[5,8] + uFrSC[5,8] >= -1
 ufrSCa: uDY[5,8] - udy[5,8] - ufrSCa[5,8] >= 0
 ufrSCa: udetYa[5,8] - ufrSCa[5,8] >= 0
 ufrSCa: - uDY[5,8] + udy[5,8] - udetYa[5,8] + ufrSCa[5,8] >= -1
 udetSC+ek: - uDetX[5,9] + uDetY[5,9] <= 0
 udetSC+ek: - ugEQK[4,9] + uDetY[5,9] <= 0
 udetSC+ek: - ugEQK[4,9] - uDetX[5,9] + uDetY[5,9] >= -1
 udetSC+ekA: - udetXa[5,9] + udetYa[5,9] <= 0
 udetSC+ekA: - ugeqka[4,9] + udetYa[5,9] <= 0
 udetSC+ekA: - ugeqka[4,9] - udetXa[5,9] + udetYa[5,9] >= -1
 ufrSC: uDY[5,9] - udy[5,9] - uFrSC[5,9] >= 0
 ufrSC: uDetY[5,9] - uFrSC[5,9] >= 0
 ufrSC: - uDY[5,9] + udy[5,9] - uDetY[5,9] + uFrSC[5,9] >= -1
 ufrSCa: uDY[5,9] - udy[5,9] - ufrSCa[5,9] >= 0
 ufrSCa: udetYa[5,9] - ufrSCa[5,9] >= 0
 ufrSCa: - uDY[5,9] + udy[5,9] - udetYa[5,9] + ufrSCa[5,9] >= -1
 udetSC+ek: - uDetX[5,10] + uDetY[5,10] <= 0
 udetSC+ek: - ugEQK[4,10] + uDetY[5,10] <= 0
 udetSC+ek: - ugEQK[4,10] - uDetX[5,10] + uDetY[5,10] >= -1
 udetSC+ekA: - udetXa[5,10] + udetYa[5,10] <= 0
 udetSC+ekA: - ugeqka[4,10] + udetYa[5,10] <= 0
 udetSC+ekA: - ugeqka[4,10] - udetXa[5,10] + udetYa[5,10] >= -1
 ufrSC: uDY[5,10] - udy[5,10] - uFrSC[5,10] >= 0
 ufrSC: uDetY[5,10] - uFrSC[5,10] >= 0
 ufrSC: - uDY[5,10] + udy[5,10] - uDetY[5,10] + uFrSC[5,10] >= -1
 ufrSCa: uDY[5,10] - udy[5,10] - ufrSCa[5,10] >= 0
 ufrSCa: udetYa[5,10] - ufrSCa[5,10] >= 0
 ufrSCa: - uDY[5,10] + udy[5,10] - udetYa[5,10] + ufrSCa[5,10] >= -1
 udetSC+ek: - uDetX[5,11] + uDetY[5,11] <= 0
 udetSC+ek: - ugEQK[4,11] + uDetY[5,11] <= 0
 udetSC+ek: - ugEQK[4,11] - uDetX[5,11] + uDetY[5,11] >= -1
 udetSC+ekA: - udetXa[5,11] + udetYa[5,11] <= 0
 udetSC+ekA: - ugeqka[4,11] + udetYa[5,11] <= 0
 udetSC+ekA: - ugeqka[4,11] - udetXa[5,11] + udetYa[5,11] >= -1
 ufrSC: uDY[5,11] - udy[5,11] - uFrSC[5,11] >= 0
 ufrSC: uDetY[5,11] - uFrSC[5,11] >= 0
 ufrSC: - uDY[5,11] + udy[5,11] - uDetY[5,11] + uFrSC[5,11] >= -1
 ufrSCa: uDY[5,11] - udy[5,11] - ufrSCa[5,11] >= 0
 ufrSCa: udetYa[5,11] - ufrSCa[5,11] >= 0
 ufrSCa: - uDY[5,11] + udy[5,11] - udetYa[5,11] + ufrSCa[5,11] >= -1
 udetSC+ek: - uDetX[5,12] + uDetY[5,12] <= 0
 udetSC+ek: - ugEQK[4,12] + uDetY[5,12] <= 0
 udetSC+ek: - ugEQK[4,12] - uDetX[5,12] + uDetY[5,12] >= -1
 udetSC+ekA: - udetXa[5,12] + udetYa[5,12] <= 0
 udetSC+ekA: - ugeqka[4,12] + udetYa[5,12] <= 0
 udetSC+ekA: - ugeqka[4,12] - udetXa[5,12] + udetYa[5,12] >= -1
 ufrSC: uDY[5,12] - udy[5,12] - uFrSC[5,12] >= 0
 ufrSC: uDetY[5,12] - uFrSC[5,12] >= 0
 ufrSC: - uDY[5,12] + udy[5,12] - uDetY[5,12] + uFrSC[5,12] >= -1
 ufrSCa: uDY[5,12] - udy[5,12] - ufrSCa[5,12] >= 0
 ufrSCa: udetYa[5,12] - ufrSCa[5,12] >= 0
 ufrSCa: - uDY[5,12] + udy[5,12] - udetYa[5,12] + ufrSCa[5,12] >= -1
 udetSC+ek: - uDetX[5,13] + uDetY[5,13] <= 0
 udetSC+ek: - ugEQK[4,13] + uDetY[5,13] <= 0
 udetSC+ek: - ugEQK[4,13] - uDetX[5,13] + uDetY[5,13] >= -1
 udetSC+ekA: - udetXa[5,13] + udetYa[5,13] <= 0
 udetSC+ekA: - ugeqka[4,13] + udetYa[5,13] <= 0
 udetSC+ekA: - ugeqka[4,13] - udetXa[5,13] + udetYa[5,13] >= -1
 ufrSC: uDY[5,13] - udy[5,13] - uFrSC[5,13] >= 0
 ufrSC: uDetY[5,13] - uFrSC[5,13] >= 0
 ufrSC: - uDY[5,13] + udy[5,13] - uDetY[5,13] + uFrSC[5,13] >= -1
 ufrSCa: uDY[5,13] - udy[5,13] - ufrSCa[5,13] >= 0
 ufrSCa: udetYa[5,13] - ufrSCa[5,13] >= 0
 ufrSCa: - uDY[5,13] + udy[5,13] - udetYa[5,13] + ufrSCa[5,13] >= -1
 udetSC+ek: - uDetX[5,14] + uDetY[5,14] <= 0
 udetSC+ek: - ugEQK[4,14] + uDetY[5,14] <= 0
 udetSC+ek: - ugEQK[4,14] - uDetX[5,14] + uDetY[5,14] >= -1
 udetSC+ekA: - udetXa[5,14] + udetYa[5,14] <= 0
 udetSC+ekA: - ugeqka[4,14] + udetYa[5,14] <= 0
 udetSC+ekA: - ugeqka[4,14] - udetXa[5,14] + udetYa[5,14] >= -1
 ufrSC: uDY[5,14] - udy[5,14] - uFrSC[5,14] >= 0
 ufrSC: uDetY[5,14] - uFrSC[5,14] >= 0
 ufrSC: - uDY[5,14] + udy[5,14] - uDetY[5,14] + uFrSC[5,14] >= -1
 ufrSCa: uDY[5,14] - udy[5,14] - ufrSCa[5,14] >= 0
 ufrSCa: udetYa[5,14] - ufrSCa[5,14] >= 0
 ufrSCa: - uDY[5,14] + udy[5,14] - udetYa[5,14] + ufrSCa[5,14] >= -1
 udetSC+ek: - uDetX[5,15] + uDetY[5,15] <= 0
 udetSC+ek: - ugEQK[4,15] + uDetY[5,15] <= 0
 udetSC+ek: - ugEQK[4,15] - uDetX[5,15] + uDetY[5,15] >= -1
 udetSC+ekA: - udetXa[5,15] + udetYa[5,15] <= 0
 udetSC+ekA: - ugeqka[4,15] + udetYa[5,15] <= 0
 udetSC+ekA: - ugeqka[4,15] - udetXa[5,15] + udetYa[5,15] >= -1
 ufrSC: uDY[5,15] - udy[5,15] - uFrSC[5,15] >= 0
 ufrSC: uDetY[5,15] - uFrSC[5,15] >= 0
 ufrSC: - uDY[5,15] + udy[5,15] - uDetY[5,15] + uFrSC[5,15] >= -1
 ufrSCa: uDY[5,15] - udy[5,15] - ufrSCa[5,15] >= 0
 ufrSCa: udetYa[5,15] - ufrSCa[5,15] >= 0
 ufrSCa: - uDY[5,15] + udy[5,15] - udetYa[5,15] + ufrSCa[5,15] >= -1
 udetSC+ek: - uDetX[6,0] + uDetY[6,0] <= 0
 udetSC+ek: - ugEQK[5,0] + uDetY[6,0] <= 0
 udetSC+ek: - ugEQK[5,0] - uDetX[6,0] + uDetY[6,0] >= -1
 udetSC+ekA: - udetXa[6,0] + udetYa[6,0] <= 0
 udetSC+ekA: - ugeqka[5,0] + udetYa[6,0] <= 0
 udetSC+ekA: - ugeqka[5,0] - udetXa[6,0] + udetYa[6,0] >= -1
 ufrSC: uDY[6,0] - udy[6,0] - uFrSC[6,0] >= 0
 ufrSC: uDetY[6,0] - uFrSC[6,0] >= 0
 ufrSC: - uDY[6,0] + udy[6,0] - uDetY[6,0] + uFrSC[6,0] >= -1
 ufrSCa: uDY[6,0] - udy[6,0] - ufrSCa[6,0] >= 0
 ufrSCa: udetYa[6,0] - ufrSCa[6,0] >= 0
 ufrSCa: - uDY[6,0] + udy[6,0] - udetYa[6,0] + ufrSCa[6,0] >= -1
 udetSC+ek: - uDetX[6,1] + uDetY[6,1] <= 0
 udetSC+ek: - ugEQK[5,1] + uDetY[6,1] <= 0
 udetSC+ek: - ugEQK[5,1] - uDetX[6,1] + uDetY[6,1] >= -1
 udetSC+ekA: - udetXa[6,1] + udetYa[6,1] <= 0
 udetSC+ekA: - ugeqka[5,1] + udetYa[6,1] <= 0
 udetSC+ekA: - ugeqka[5,1] - udetXa[6,1] + udetYa[6,1] >= -1
 ufrSC: uDY[6,1] - udy[6,1] - uFrSC[6,1] >= 0
 ufrSC: uDetY[6,1] - uFrSC[6,1] >= 0
 ufrSC: - uDY[6,1] + udy[6,1] - uDetY[6,1] + uFrSC[6,1] >= -1
 ufrSCa: uDY[6,1] - udy[6,1] - ufrSCa[6,1] >= 0
 ufrSCa: udetYa[6,1] - ufrSCa[6,1] >= 0
 ufrSCa: - uDY[6,1] + udy[6,1] - udetYa[6,1] + ufrSCa[6,1] >= -1
 udetSC+ek: - uDetX[6,2] + uDetY[6,2] <= 0
 udetSC+ek: - ugEQK[5,2] + uDetY[6,2] <= 0
 udetSC+ek: - ugEQK[5,2] - uDetX[6,2] + uDetY[6,2] >= -1
 udetSC+ekA: - udetXa[6,2] + udetYa[6,2] <= 0
 udetSC+ekA: - ugeqka[5,2] + udetYa[6,2] <= 0
 udetSC+ekA: - ugeqka[5,2] - udetXa[6,2] + udetYa[6,2] >= -1
 ufrSC: uDY[6,2] - udy[6,2] - uFrSC[6,2] >= 0
 ufrSC: uDetY[6,2] - uFrSC[6,2] >= 0
 ufrSC: - uDY[6,2] + udy[6,2] - uDetY[6,2] + uFrSC[6,2] >= -1
 ufrSCa: uDY[6,2] - udy[6,2] - ufrSCa[6,2] >= 0
 ufrSCa: udetYa[6,2] - ufrSCa[6,2] >= 0
 ufrSCa: - uDY[6,2] + udy[6,2] - udetYa[6,2] + ufrSCa[6,2] >= -1
 udetSC+ek: - uDetX[6,3] + uDetY[6,3] <= 0
 udetSC+ek: - ugEQK[5,3] + uDetY[6,3] <= 0
 udetSC+ek: - ugEQK[5,3] - uDetX[6,3] + uDetY[6,3] >= -1
 udetSC+ekA: - udetXa[6,3] + udetYa[6,3] <= 0
 udetSC+ekA: - ugeqka[5,3] + udetYa[6,3] <= 0
 udetSC+ekA: - ugeqka[5,3] - udetXa[6,3] + udetYa[6,3] >= -1
 ufrSC: uDY[6,3] - udy[6,3] - uFrSC[6,3] >= 0
 ufrSC: uDetY[6,3] - uFrSC[6,3] >= 0
 ufrSC: - uDY[6,3] + udy[6,3] - uDetY[6,3] + uFrSC[6,3] >= -1
 ufrSCa: uDY[6,3] - udy[6,3] - ufrSCa[6,3] >= 0
 ufrSCa: udetYa[6,3] - ufrSCa[6,3] >= 0
 ufrSCa: - uDY[6,3] + udy[6,3] - udetYa[6,3] + ufrSCa[6,3] >= -1
 udetSC+ek: - uDetX[6,4] + uDetY[6,4] <= 0
 udetSC+ek: - ugEQK[5,4] + uDetY[6,4] <= 0
 udetSC+ek: - ugEQK[5,4] - uDetX[6,4] + uDetY[6,4] >= -1
 udetSC+ekA: - udetXa[6,4] + udetYa[6,4] <= 0
 udetSC+ekA: - ugeqka[5,4] + udetYa[6,4] <= 0
 udetSC+ekA: - ugeqka[5,4] - udetXa[6,4] + udetYa[6,4] >= -1
 ufrSC: uDY[6,4] - udy[6,4] - uFrSC[6,4] >= 0
 ufrSC: uDetY[6,4] - uFrSC[6,4] >= 0
 ufrSC: - uDY[6,4] + udy[6,4] - uDetY[6,4] + uFrSC[6,4] >= -1
 ufrSCa: uDY[6,4] - udy[6,4] - ufrSCa[6,4] >= 0
 ufrSCa: udetYa[6,4] - ufrSCa[6,4] >= 0
 ufrSCa: - uDY[6,4] + udy[6,4] - udetYa[6,4] + ufrSCa[6,4] >= -1
 udetSC+ek: - uDetX[6,5] + uDetY[6,5] <= 0
 udetSC+ek: - ugEQK[5,5] + uDetY[6,5] <= 0
 udetSC+ek: - ugEQK[5,5] - uDetX[6,5] + uDetY[6,5] >= -1
 udetSC+ekA: - udetXa[6,5] + udetYa[6,5] <= 0
 udetSC+ekA: - ugeqka[5,5] + udetYa[6,5] <= 0
 udetSC+ekA: - ugeqka[5,5] - udetXa[6,5] + udetYa[6,5] >= -1
 ufrSC: uDY[6,5] - udy[6,5] - uFrSC[6,5] >= 0
 ufrSC: uDetY[6,5] - uFrSC[6,5] >= 0
 ufrSC: - uDY[6,5] + udy[6,5] - uDetY[6,5] + uFrSC[6,5] >= -1
 ufrSCa: uDY[6,5] - udy[6,5] - ufrSCa[6,5] >= 0
 ufrSCa: udetYa[6,5] - ufrSCa[6,5] >= 0
 ufrSCa: - uDY[6,5] + udy[6,5] - udetYa[6,5] + ufrSCa[6,5] >= -1
 udetSC+ek: - uDetX[6,6] + uDetY[6,6] <= 0
 udetSC+ek: - ugEQK[5,6] + uDetY[6,6] <= 0
 udetSC+ek: - ugEQK[5,6] - uDetX[6,6] + uDetY[6,6] >= -1
 udetSC+ekA: - udetXa[6,6] + udetYa[6,6] <= 0
 udetSC+ekA: - ugeqka[5,6] + udetYa[6,6] <= 0
 udetSC+ekA: - ugeqka[5,6] - udetXa[6,6] + udetYa[6,6] >= -1
 ufrSC: uDY[6,6] - udy[6,6] - uFrSC[6,6] >= 0
 ufrSC: uDetY[6,6] - uFrSC[6,6] >= 0
 ufrSC: - uDY[6,6] + udy[6,6] - uDetY[6,6] + uFrSC[6,6] >= -1
 ufrSCa: uDY[6,6] - udy[6,6] - ufrSCa[6,6] >= 0
 ufrSCa: udetYa[6,6] - ufrSCa[6,6] >= 0
 ufrSCa: - uDY[6,6] + udy[6,6] - udetYa[6,6] + ufrSCa[6,6] >= -1
 udetSC+ek: - uDetX[6,7] + uDetY[6,7] <= 0
 udetSC+ek: - ugEQK[5,7] + uDetY[6,7] <= 0
 udetSC+ek: - ugEQK[5,7] - uDetX[6,7] + uDetY[6,7] >= -1
 udetSC+ekA: - udetXa[6,7] + udetYa[6,7] <= 0
 udetSC+ekA: - ugeqka[5,7] + udetYa[6,7] <= 0
 udetSC+ekA: - ugeqka[5,7] - udetXa[6,7] + udetYa[6,7] >= -1
 ufrSC: uDY[6,7] - udy[6,7] - uFrSC[6,7] >= 0
 ufrSC: uDetY[6,7] - uFrSC[6,7] >= 0
 ufrSC: - uDY[6,7] + udy[6,7] - uDetY[6,7] + uFrSC[6,7] >= -1
 ufrSCa: uDY[6,7] - udy[6,7] - ufrSCa[6,7] >= 0
 ufrSCa: udetYa[6,7] - ufrSCa[6,7] >= 0
 ufrSCa: - uDY[6,7] + udy[6,7] - udetYa[6,7] + ufrSCa[6,7] >= -1
 udetSC+ek: - uDetX[6,8] + uDetY[6,8] <= 0
 udetSC+ek: - ugEQK[5,8] + uDetY[6,8] <= 0
 udetSC+ek: - ugEQK[5,8] - uDetX[6,8] + uDetY[6,8] >= -1
 udetSC+ekA: - udetXa[6,8] + udetYa[6,8] <= 0
 udetSC+ekA: - ugeqka[5,8] + udetYa[6,8] <= 0
 udetSC+ekA: - ugeqka[5,8] - udetXa[6,8] + udetYa[6,8] >= -1
 ufrSC: uDY[6,8] - udy[6,8] - uFrSC[6,8] >= 0
 ufrSC: uDetY[6,8] - uFrSC[6,8] >= 0
 ufrSC: - uDY[6,8] + udy[6,8] - uDetY[6,8] + uFrSC[6,8] >= -1
 ufrSCa: uDY[6,8] - udy[6,8] - ufrSCa[6,8] >= 0
 ufrSCa: udetYa[6,8] - ufrSCa[6,8] >= 0
 ufrSCa: - uDY[6,8] + udy[6,8] - udetYa[6,8] + ufrSCa[6,8] >= -1
 udetSC+ek: - uDetX[6,9] + uDetY[6,9] <= 0
 udetSC+ek: - ugEQK[5,9] + uDetY[6,9] <= 0
 udetSC+ek: - ugEQK[5,9] - uDetX[6,9] + uDetY[6,9] >= -1
 udetSC+ekA: - udetXa[6,9] + udetYa[6,9] <= 0
 udetSC+ekA: - ugeqka[5,9] + udetYa[6,9] <= 0
 udetSC+ekA: - ugeqka[5,9] - udetXa[6,9] + udetYa[6,9] >= -1
 ufrSC: uDY[6,9] - udy[6,9] - uFrSC[6,9] >= 0
 ufrSC: uDetY[6,9] - uFrSC[6,9] >= 0
 ufrSC: - uDY[6,9] + udy[6,9] - uDetY[6,9] + uFrSC[6,9] >= -1
 ufrSCa: uDY[6,9] - udy[6,9] - ufrSCa[6,9] >= 0
 ufrSCa: udetYa[6,9] - ufrSCa[6,9] >= 0
 ufrSCa: - uDY[6,9] + udy[6,9] - udetYa[6,9] + ufrSCa[6,9] >= -1
 udetSC+ek: - uDetX[6,10] + uDetY[6,10] <= 0
 udetSC+ek: - ugEQK[5,10] + uDetY[6,10] <= 0
 udetSC+ek: - ugEQK[5,10] - uDetX[6,10] + uDetY[6,10] >= -1
 udetSC+ekA: - udetXa[6,10] + udetYa[6,10] <= 0
 udetSC+ekA: - ugeqka[5,10] + udetYa[6,10] <= 0
 udetSC+ekA: - ugeqka[5,10] - udetXa[6,10] + udetYa[6,10] >= -1
 ufrSC: uDY[6,10] - udy[6,10] - uFrSC[6,10] >= 0
 ufrSC: uDetY[6,10] - uFrSC[6,10] >= 0
 ufrSC: - uDY[6,10] + udy[6,10] - uDetY[6,10] + uFrSC[6,10] >= -1
 ufrSCa: uDY[6,10] - udy[6,10] - ufrSCa[6,10] >= 0
 ufrSCa: udetYa[6,10] - ufrSCa[6,10] >= 0
 ufrSCa: - uDY[6,10] + udy[6,10] - udetYa[6,10] + ufrSCa[6,10] >= -1
 udetSC+ek: - uDetX[6,11] + uDetY[6,11] <= 0
 udetSC+ek: - ugEQK[5,11] + uDetY[6,11] <= 0
 udetSC+ek: - ugEQK[5,11] - uDetX[6,11] + uDetY[6,11] >= -1
 udetSC+ekA: - udetXa[6,11] + udetYa[6,11] <= 0
 udetSC+ekA: - ugeqka[5,11] + udetYa[6,11] <= 0
 udetSC+ekA: - ugeqka[5,11] - udetXa[6,11] + udetYa[6,11] >= -1
 ufrSC: uDY[6,11] - udy[6,11] - uFrSC[6,11] >= 0
 ufrSC: uDetY[6,11] - uFrSC[6,11] >= 0
 ufrSC: - uDY[6,11] + udy[6,11] - uDetY[6,11] + uFrSC[6,11] >= -1
 ufrSCa: uDY[6,11] - udy[6,11] - ufrSCa[6,11] >= 0
 ufrSCa: udetYa[6,11] - ufrSCa[6,11] >= 0
 ufrSCa: - uDY[6,11] + udy[6,11] - udetYa[6,11] + ufrSCa[6,11] >= -1
 udetSC+ek: - uDetX[6,12] + uDetY[6,12] <= 0
 udetSC+ek: - ugEQK[5,12] + uDetY[6,12] <= 0
 udetSC+ek: - ugEQK[5,12] - uDetX[6,12] + uDetY[6,12] >= -1
 udetSC+ekA: - udetXa[6,12] + udetYa[6,12] <= 0
 udetSC+ekA: - ugeqka[5,12] + udetYa[6,12] <= 0
 udetSC+ekA: - ugeqka[5,12] - udetXa[6,12] + udetYa[6,12] >= -1
 ufrSC: uDY[6,12] - udy[6,12] - uFrSC[6,12] >= 0
 ufrSC: uDetY[6,12] - uFrSC[6,12] >= 0
 ufrSC: - uDY[6,12] + udy[6,12] - uDetY[6,12] + uFrSC[6,12] >= -1
 ufrSCa: uDY[6,12] - udy[6,12] - ufrSCa[6,12] >= 0
 ufrSCa: udetYa[6,12] - ufrSCa[6,12] >= 0
 ufrSCa: - uDY[6,12] + udy[6,12] - udetYa[6,12] + ufrSCa[6,12] >= -1
 udetSC+ek: - uDetX[6,13] + uDetY[6,13] <= 0
 udetSC+ek: - ugEQK[5,13] + uDetY[6,13] <= 0
 udetSC+ek: - ugEQK[5,13] - uDetX[6,13] + uDetY[6,13] >= -1
 udetSC+ekA: - udetXa[6,13] + udetYa[6,13] <= 0
 udetSC+ekA: - ugeqka[5,13] + udetYa[6,13] <= 0
 udetSC+ekA: - ugeqka[5,13] - udetXa[6,13] + udetYa[6,13] >= -1
 ufrSC: uDY[6,13] - udy[6,13] - uFrSC[6,13] >= 0
 ufrSC: uDetY[6,13] - uFrSC[6,13] >= 0
 ufrSC: - uDY[6,13] + udy[6,13] - uDetY[6,13] + uFrSC[6,13] >= -1
 ufrSCa: uDY[6,13] - udy[6,13] - ufrSCa[6,13] >= 0
 ufrSCa: udetYa[6,13] - ufrSCa[6,13] >= 0
 ufrSCa: - uDY[6,13] + udy[6,13] - udetYa[6,13] + ufrSCa[6,13] >= -1
 udetSC+ek: - uDetX[6,14] + uDetY[6,14] <= 0
 udetSC+ek: - ugEQK[5,14] + uDetY[6,14] <= 0
 udetSC+ek: - ugEQK[5,14] - uDetX[6,14] + uDetY[6,14] >= -1
 udetSC+ekA: - udetXa[6,14] + udetYa[6,14] <= 0
 udetSC+ekA: - ugeqka[5,14] + udetYa[6,14] <= 0
 udetSC+ekA: - ugeqka[5,14] - udetXa[6,14] + udetYa[6,14] >= -1
 ufrSC: uDY[6,14] - udy[6,14] - uFrSC[6,14] >= 0
 ufrSC: uDetY[6,14] - uFrSC[6,14] >= 0
 ufrSC: - uDY[6,14] + udy[6,14] - uDetY[6,14] + uFrSC[6,14] >= -1
 ufrSCa: uDY[6,14] - udy[6,14] - ufrSCa[6,14] >= 0
 ufrSCa: udetYa[6,14] - ufrSCa[6,14] >= 0
 ufrSCa: - uDY[6,14] + udy[6,14] - udetYa[6,14] + ufrSCa[6,14] >= -1
 udetSC+ek: - uDetX[6,15] + uDetY[6,15] <= 0
 udetSC+ek: - ugEQK[5,15] + uDetY[6,15] <= 0
 udetSC+ek: - ugEQK[5,15] - uDetX[6,15] + uDetY[6,15] >= -1
 udetSC+ekA: - udetXa[6,15] + udetYa[6,15] <= 0
 udetSC+ekA: - ugeqka[5,15] + udetYa[6,15] <= 0
 udetSC+ekA: - ugeqka[5,15] - udetXa[6,15] + udetYa[6,15] >= -1
 ufrSC: uDY[6,15] - udy[6,15] - uFrSC[6,15] >= 0
 ufrSC: uDetY[6,15] - uFrSC[6,15] >= 0
 ufrSC: - uDY[6,15] + udy[6,15] - uDetY[6,15] + uFrSC[6,15] >= -1
 ufrSCa: uDY[6,15] - udy[6,15] - ufrSCa[6,15] >= 0
 ufrSCa: udetYa[6,15] - ufrSCa[6,15] >= 0
 ufrSCa: - uDY[6,15] + udy[6,15] - udetYa[6,15] + ufrSCa[6,15] >= -1
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
 udetY2Z: uDetY[4,0] - uDetZ[4,0] = 0
 udety2zA: udetYa[4,0] - udetZa[4,0] = 0
 udetY2Z: uDetY[4,1] - uDetZ[4,1] = 0
 udety2zA: udetYa[4,1] - udetZa[4,1] = 0
 udetY2Z: uDetY[4,2] - uDetZ[4,2] = 0
 udety2zA: udetYa[4,2] - udetZa[4,2] = 0
 udetY2Z: uDetY[4,3] - uDetZ[4,3] = 0
 udety2zA: udetYa[4,3] - udetZa[4,3] = 0
 udetY2Z: uDetY[4,4] - uDetZ[4,4] = 0
 udety2zA: udetYa[4,4] - udetZa[4,4] = 0
 udetY2Z: uDetY[4,5] - uDetZ[4,5] = 0
 udety2zA: udetYa[4,5] - udetZa[4,5] = 0
 udetY2Z: uDetY[4,6] - uDetZ[4,6] = 0
 udety2zA: udetYa[4,6] - udetZa[4,6] = 0
 udetY2Z: uDetY[4,7] - uDetZ[4,7] = 0
 udety2zA: udetYa[4,7] - udetZa[4,7] = 0
 udetY2Z: uDetY[4,8] - uDetZ[4,8] = 0
 udety2zA: udetYa[4,8] - udetZa[4,8] = 0
 udetY2Z: uDetY[4,9] - uDetZ[4,9] = 0
 udety2zA: udetYa[4,9] - udetZa[4,9] = 0
 udetY2Z: uDetY[4,10] - uDetZ[4,10] = 0
 udety2zA: udetYa[4,10] - udetZa[4,10] = 0
 udetY2Z: uDetY[4,11] - uDetZ[4,11] = 0
 udety2zA: udetYa[4,11] - udetZa[4,11] = 0
 udetY2Z: uDetY[4,12] - uDetZ[4,12] = 0
 udety2zA: udetYa[4,12] - udetZa[4,12] = 0
 udetY2Z: uDetY[4,13] - uDetZ[4,13] = 0
 udety2zA: udetYa[4,13] - udetZa[4,13] = 0
 udetY2Z: uDetY[4,14] - uDetZ[4,14] = 0
 udety2zA: udetYa[4,14] - udetZa[4,14] = 0
 udetY2Z: uDetY[4,15] - uDetZ[4,15] = 0
 udety2zA: udetYa[4,15] - udetZa[4,15] = 0
 udetY2Z: uDetY[5,0] - uDetZ[5,0] = 0
 udety2zA: udetYa[5,0] - udetZa[5,0] = 0
 udetY2Z: uDetY[5,1] - uDetZ[5,1] = 0
 udety2zA: udetYa[5,1] - udetZa[5,1] = 0
 udetY2Z: uDetY[5,2] - uDetZ[5,2] = 0
 udety2zA: udetYa[5,2] - udetZa[5,2] = 0
 udetY2Z: uDetY[5,3] - uDetZ[5,3] = 0
 udety2zA: udetYa[5,3] - udetZa[5,3] = 0
 udetY2Z: uDetY[5,4] - uDetZ[5,4] = 0
 udety2zA: udetYa[5,4] - udetZa[5,4] = 0
 udetY2Z: uDetY[5,5] - uDetZ[5,5] = 0
 udety2zA: udetYa[5,5] - udetZa[5,5] = 0
 udetY2Z: uDetY[5,6] - uDetZ[5,6] = 0
 udety2zA: udetYa[5,6] - udetZa[5,6] = 0
 udetY2Z: uDetY[5,7] - uDetZ[5,7] = 0
 udety2zA: udetYa[5,7] - udetZa[5,7] = 0
 udetY2Z: uDetY[5,8] - uDetZ[5,8] = 0
 udety2zA: udetYa[5,8] - udetZa[5,8] = 0
 udetY2Z: uDetY[5,9] - uDetZ[5,9] = 0
 udety2zA: udetYa[5,9] - udetZa[5,9] = 0
 udetY2Z: uDetY[5,10] - uDetZ[5,10] = 0
 udety2zA: udetYa[5,10] - udetZa[5,10] = 0
 udetY2Z: uDetY[5,11] - uDetZ[5,11] = 0
 udety2zA: udetYa[5,11] - udetZa[5,11] = 0
 udetY2Z: uDetY[5,12] - uDetZ[5,12] = 0
 udety2zA: udetYa[5,12] - udetZa[5,12] = 0
 udetY2Z: uDetY[5,13] - uDetZ[5,13] = 0
 udety2zA: udetYa[5,13] - udetZa[5,13] = 0
 udetY2Z: uDetY[5,14] - uDetZ[5,14] = 0
 udety2zA: udetYa[5,14] - udetZa[5,14] = 0
 udetY2Z: uDetY[5,15] - uDetZ[5,15] = 0
 udety2zA: udetYa[5,15] - udetZa[5,15] = 0
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
 udetSR: - uDetZ[4,0] + uDetW[4,0] = 0
 udetsrA: - udetZa[4,0] + udetWa[4,0] = 0
 udetSR: - uDetZ[4,1] + uDetW[4,1] = 0
 udetsrA: - udetZa[4,1] + udetWa[4,1] = 0
 udetSR: - uDetZ[4,2] + uDetW[4,2] = 0
 udetsrA: - udetZa[4,2] + udetWa[4,2] = 0
 udetSR: - uDetZ[4,3] + uDetW[4,3] = 0
 udetsrA: - udetZa[4,3] + udetWa[4,3] = 0
 udetSR: - uDetZ[4,7] + uDetW[4,4] = 0
 udetsrA: - udetZa[4,7] + udetWa[4,4] = 0
 udetSR: - uDetZ[4,4] + uDetW[4,5] = 0
 udetsrA: - udetZa[4,4] + udetWa[4,5] = 0
 udetSR: - uDetZ[4,5] + uDetW[4,6] = 0
 udetsrA: - udetZa[4,5] + udetWa[4,6] = 0
 udetSR: - uDetZ[4,6] + uDetW[4,7] = 0
 udetsrA: - udetZa[4,6] + udetWa[4,7] = 0
 udetSR: - uDetZ[4,10] + uDetW[4,8] = 0
 udetsrA: - udetZa[4,10] + udetWa[4,8] = 0
 udetSR: - uDetZ[4,11] + uDetW[4,9] = 0
 udetsrA: - udetZa[4,11] + udetWa[4,9] = 0
 udetSR: - uDetZ[4,8] + uDetW[4,10] = 0
 udetsrA: - udetZa[4,8] + udetWa[4,10] = 0
 udetSR: - uDetZ[4,9] + uDetW[4,11] = 0
 udetsrA: - udetZa[4,9] + udetWa[4,11] = 0
 udetSR: - uDetZ[4,13] + uDetW[4,12] = 0
 udetsrA: - udetZa[4,13] + udetWa[4,12] = 0
 udetSR: - uDetZ[4,14] + uDetW[4,13] = 0
 udetsrA: - udetZa[4,14] + udetWa[4,13] = 0
 udetSR: - uDetZ[4,15] + uDetW[4,14] = 0
 udetsrA: - udetZa[4,15] + udetWa[4,14] = 0
 udetSR: - uDetZ[4,12] + uDetW[4,15] = 0
 udetsrA: - udetZa[4,12] + udetWa[4,15] = 0
 udetSR: - uDetZ[5,0] + uDetW[5,0] = 0
 udetsrA: - udetZa[5,0] + udetWa[5,0] = 0
 udetSR: - uDetZ[5,1] + uDetW[5,1] = 0
 udetsrA: - udetZa[5,1] + udetWa[5,1] = 0
 udetSR: - uDetZ[5,2] + uDetW[5,2] = 0
 udetsrA: - udetZa[5,2] + udetWa[5,2] = 0
 udetSR: - uDetZ[5,3] + uDetW[5,3] = 0
 udetsrA: - udetZa[5,3] + udetWa[5,3] = 0
 udetSR: - uDetZ[5,7] + uDetW[5,4] = 0
 udetsrA: - udetZa[5,7] + udetWa[5,4] = 0
 udetSR: - uDetZ[5,4] + uDetW[5,5] = 0
 udetsrA: - udetZa[5,4] + udetWa[5,5] = 0
 udetSR: - uDetZ[5,5] + uDetW[5,6] = 0
 udetsrA: - udetZa[5,5] + udetWa[5,6] = 0
 udetSR: - uDetZ[5,6] + uDetW[5,7] = 0
 udetsrA: - udetZa[5,6] + udetWa[5,7] = 0
 udetSR: - uDetZ[5,10] + uDetW[5,8] = 0
 udetsrA: - udetZa[5,10] + udetWa[5,8] = 0
 udetSR: - uDetZ[5,11] + uDetW[5,9] = 0
 udetsrA: - udetZa[5,11] + udetWa[5,9] = 0
 udetSR: - uDetZ[5,8] + uDetW[5,10] = 0
 udetsrA: - udetZa[5,8] + udetWa[5,10] = 0
 udetSR: - uDetZ[5,9] + uDetW[5,11] = 0
 udetsrA: - udetZa[5,9] + udetWa[5,11] = 0
 udetSR: - uDetZ[5,13] + uDetW[5,12] = 0
 udetsrA: - udetZa[5,13] + udetWa[5,12] = 0
 udetSR: - uDetZ[5,14] + uDetW[5,13] = 0
 udetsrA: - udetZa[5,14] + udetWa[5,13] = 0
 udetSR: - uDetZ[5,15] + uDetW[5,14] = 0
 udetsrA: - udetZa[5,15] + udetWa[5,14] = 0
 udetSR: - uDetZ[5,12] + uDetW[5,15] = 0
 udetsrA: - udetZa[5,12] + udetWa[5,15] = 0
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
 udetMC0: uDetX[5,0] - uDetW[4,12] <= 0
 udetMC0: uDetX[5,0] - uDetX[5,12] <= 0
 udetMC0: uDetX[5,0] - uDetX[5,12] - uDetW[4,12] >= -1
 udetmc0A: udetXa[5,0] - udetWa[4,12] <= 0
 udetmc0A: udetXa[5,0] - udetXa[5,12] <= 0
 udetmc0A: udetXa[5,0] - udetXa[5,12] - udetWa[4,12] >= -1
 udetMC1: uDetX[5,4] - uDetW[4,0] = 0
 udetmc1A: udetXa[5,4] - udetWa[4,0] = 0
 udetMC2: uDetX[5,8] - uDetW[4,4] <= 0
 udetMC2: uDetX[5,8] - uDetW[4,8] <= 0
 udetMC2: uDetX[5,8] - uDetW[4,4] - uDetW[4,8] >= -1
 udetmc2A: udetXa[5,8] - udetWa[4,4] <= 0
 udetmc2A: udetXa[5,8] - udetWa[4,8] <= 0
 udetmc2A: udetXa[5,8] - udetWa[4,4] - udetWa[4,8] >= -1
 udetMC3: uDetX[5,12] - uDetW[4,0] <= 0
 udetMC3: uDetX[5,12] - uDetW[4,8] <= 0
 udetMC3: uDetX[5,12] - uDetW[4,0] - uDetW[4,8] >= -1
 udetmc3A: udetXa[5,12] - udetWa[4,0] <= 0
 udetmc3A: udetXa[5,12] - udetWa[4,8] <= 0
 udetmc3A: udetXa[5,12] - udetWa[4,0] - udetWa[4,8] >= -1
 udetMC0: uDetX[5,1] - uDetW[4,13] <= 0
 udetMC0: uDetX[5,1] - uDetX[5,13] <= 0
 udetMC0: uDetX[5,1] - uDetX[5,13] - uDetW[4,13] >= -1
 udetmc0A: udetXa[5,1] - udetWa[4,13] <= 0
 udetmc0A: udetXa[5,1] - udetXa[5,13] <= 0
 udetmc0A: udetXa[5,1] - udetXa[5,13] - udetWa[4,13] >= -1
 udetMC1: uDetX[5,5] - uDetW[4,1] = 0
 udetmc1A: udetXa[5,5] - udetWa[4,1] = 0
 udetMC2: uDetX[5,9] - uDetW[4,5] <= 0
 udetMC2: uDetX[5,9] - uDetW[4,9] <= 0
 udetMC2: uDetX[5,9] - uDetW[4,5] - uDetW[4,9] >= -1
 udetmc2A: udetXa[5,9] - udetWa[4,5] <= 0
 udetmc2A: udetXa[5,9] - udetWa[4,9] <= 0
 udetmc2A: udetXa[5,9] - udetWa[4,5] - udetWa[4,9] >= -1
 udetMC3: uDetX[5,13] - uDetW[4,1] <= 0
 udetMC3: uDetX[5,13] - uDetW[4,9] <= 0
 udetMC3: uDetX[5,13] - uDetW[4,1] - uDetW[4,9] >= -1
 udetmc3A: udetXa[5,13] - udetWa[4,1] <= 0
 udetmc3A: udetXa[5,13] - udetWa[4,9] <= 0
 udetmc3A: udetXa[5,13] - udetWa[4,1] - udetWa[4,9] >= -1
 udetMC0: uDetX[5,2] - uDetW[4,14] <= 0
 udetMC0: uDetX[5,2] - uDetX[5,14] <= 0
 udetMC0: uDetX[5,2] - uDetX[5,14] - uDetW[4,14] >= -1
 udetmc0A: udetXa[5,2] - udetWa[4,14] <= 0
 udetmc0A: udetXa[5,2] - udetXa[5,14] <= 0
 udetmc0A: udetXa[5,2] - udetXa[5,14] - udetWa[4,14] >= -1
 udetMC1: uDetX[5,6] - uDetW[4,2] = 0
 udetmc1A: udetXa[5,6] - udetWa[4,2] = 0
 udetMC2: uDetX[5,10] - uDetW[4,6] <= 0
 udetMC2: uDetX[5,10] - uDetW[4,10] <= 0
 udetMC2: uDetX[5,10] - uDetW[4,6] - uDetW[4,10] >= -1
 udetmc2A: udetXa[5,10] - udetWa[4,6] <= 0
 udetmc2A: udetXa[5,10] - udetWa[4,10] <= 0
 udetmc2A: udetXa[5,10] - udetWa[4,6] - udetWa[4,10] >= -1
 udetMC3: uDetX[5,14] - uDetW[4,2] <= 0
 udetMC3: uDetX[5,14] - uDetW[4,10] <= 0
 udetMC3: uDetX[5,14] - uDetW[4,2] - uDetW[4,10] >= -1
 udetmc3A: udetXa[5,14] - udetWa[4,2] <= 0
 udetmc3A: udetXa[5,14] - udetWa[4,10] <= 0
 udetmc3A: udetXa[5,14] - udetWa[4,2] - udetWa[4,10] >= -1
 udetMC0: uDetX[5,3] - uDetW[4,15] <= 0
 udetMC0: uDetX[5,3] - uDetX[5,15] <= 0
 udetMC0: uDetX[5,3] - uDetX[5,15] - uDetW[4,15] >= -1
 udetmc0A: udetXa[5,3] - udetWa[4,15] <= 0
 udetmc0A: udetXa[5,3] - udetXa[5,15] <= 0
 udetmc0A: udetXa[5,3] - udetXa[5,15] - udetWa[4,15] >= -1
 udetMC1: uDetX[5,7] - uDetW[4,3] = 0
 udetmc1A: udetXa[5,7] - udetWa[4,3] = 0
 udetMC2: uDetX[5,11] - uDetW[4,7] <= 0
 udetMC2: uDetX[5,11] - uDetW[4,11] <= 0
 udetMC2: uDetX[5,11] - uDetW[4,7] - uDetW[4,11] >= -1
 udetmc2A: udetXa[5,11] - udetWa[4,7] <= 0
 udetmc2A: udetXa[5,11] - udetWa[4,11] <= 0
 udetmc2A: udetXa[5,11] - udetWa[4,7] - udetWa[4,11] >= -1
 udetMC3: uDetX[5,15] - uDetW[4,3] <= 0
 udetMC3: uDetX[5,15] - uDetW[4,11] <= 0
 udetMC3: uDetX[5,15] - uDetW[4,3] - uDetW[4,11] >= -1
 udetmc3A: udetXa[5,15] - udetWa[4,3] <= 0
 udetmc3A: udetXa[5,15] - udetWa[4,11] <= 0
 udetmc3A: udetXa[5,15] - udetWa[4,3] - udetWa[4,11] >= -1
 udetMC0: uDetX[6,0] - uDetW[5,12] <= 0
 udetMC0: uDetX[6,0] - uDetX[6,12] <= 0
 udetMC0: uDetX[6,0] - uDetX[6,12] - uDetW[5,12] >= -1
 udetmc0A: udetXa[6,0] - udetWa[5,12] <= 0
 udetmc0A: udetXa[6,0] - udetXa[6,12] <= 0
 udetmc0A: udetXa[6,0] - udetXa[6,12] - udetWa[5,12] >= -1
 udetMC1: uDetX[6,4] - uDetW[5,0] = 0
 udetmc1A: udetXa[6,4] - udetWa[5,0] = 0
 udetMC2: uDetX[6,8] - uDetW[5,4] <= 0
 udetMC2: uDetX[6,8] - uDetW[5,8] <= 0
 udetMC2: uDetX[6,8] - uDetW[5,4] - uDetW[5,8] >= -1
 udetmc2A: udetXa[6,8] - udetWa[5,4] <= 0
 udetmc2A: udetXa[6,8] - udetWa[5,8] <= 0
 udetmc2A: udetXa[6,8] - udetWa[5,4] - udetWa[5,8] >= -1
 udetMC3: uDetX[6,12] - uDetW[5,0] <= 0
 udetMC3: uDetX[6,12] - uDetW[5,8] <= 0
 udetMC3: uDetX[6,12] - uDetW[5,0] - uDetW[5,8] >= -1
 udetmc3A: udetXa[6,12] - udetWa[5,0] <= 0
 udetmc3A: udetXa[6,12] - udetWa[5,8] <= 0
 udetmc3A: udetXa[6,12] - udetWa[5,0] - udetWa[5,8] >= -1
 udetMC0: uDetX[6,1] - uDetW[5,13] <= 0
 udetMC0: uDetX[6,1] - uDetX[6,13] <= 0
 udetMC0: uDetX[6,1] - uDetX[6,13] - uDetW[5,13] >= -1
 udetmc0A: udetXa[6,1] - udetWa[5,13] <= 0
 udetmc0A: udetXa[6,1] - udetXa[6,13] <= 0
 udetmc0A: udetXa[6,1] - udetXa[6,13] - udetWa[5,13] >= -1
 udetMC1: uDetX[6,5] - uDetW[5,1] = 0
 udetmc1A: udetXa[6,5] - udetWa[5,1] = 0
 udetMC2: uDetX[6,9] - uDetW[5,5] <= 0
 udetMC2: uDetX[6,9] - uDetW[5,9] <= 0
 udetMC2: uDetX[6,9] - uDetW[5,5] - uDetW[5,9] >= -1
 udetmc2A: udetXa[6,9] - udetWa[5,5] <= 0
 udetmc2A: udetXa[6,9] - udetWa[5,9] <= 0
 udetmc2A: udetXa[6,9] - udetWa[5,5] - udetWa[5,9] >= -1
 udetMC3: uDetX[6,13] - uDetW[5,1] <= 0
 udetMC3: uDetX[6,13] - uDetW[5,9] <= 0
 udetMC3: uDetX[6,13] - uDetW[5,1] - uDetW[5,9] >= -1
 udetmc3A: udetXa[6,13] - udetWa[5,1] <= 0
 udetmc3A: udetXa[6,13] - udetWa[5,9] <= 0
 udetmc3A: udetXa[6,13] - udetWa[5,1] - udetWa[5,9] >= -1
 udetMC0: uDetX[6,2] - uDetW[5,14] <= 0
 udetMC0: uDetX[6,2] - uDetX[6,14] <= 0
 udetMC0: uDetX[6,2] - uDetX[6,14] - uDetW[5,14] >= -1
 udetmc0A: udetXa[6,2] - udetWa[5,14] <= 0
 udetmc0A: udetXa[6,2] - udetXa[6,14] <= 0
 udetmc0A: udetXa[6,2] - udetXa[6,14] - udetWa[5,14] >= -1
 udetMC1: uDetX[6,6] - uDetW[5,2] = 0
 udetmc1A: udetXa[6,6] - udetWa[5,2] = 0
 udetMC2: uDetX[6,10] - uDetW[5,6] <= 0
 udetMC2: uDetX[6,10] - uDetW[5,10] <= 0
 udetMC2: uDetX[6,10] - uDetW[5,6] - uDetW[5,10] >= -1
 udetmc2A: udetXa[6,10] - udetWa[5,6] <= 0
 udetmc2A: udetXa[6,10] - udetWa[5,10] <= 0
 udetmc2A: udetXa[6,10] - udetWa[5,6] - udetWa[5,10] >= -1
 udetMC3: uDetX[6,14] - uDetW[5,2] <= 0
 udetMC3: uDetX[6,14] - uDetW[5,10] <= 0
 udetMC3: uDetX[6,14] - uDetW[5,2] - uDetW[5,10] >= -1
 udetmc3A: udetXa[6,14] - udetWa[5,2] <= 0
 udetmc3A: udetXa[6,14] - udetWa[5,10] <= 0
 udetmc3A: udetXa[6,14] - udetWa[5,2] - udetWa[5,10] >= -1
 udetMC0: uDetX[6,3] - uDetW[5,15] <= 0
 udetMC0: uDetX[6,3] - uDetX[6,15] <= 0
 udetMC0: uDetX[6,3] - uDetX[6,15] - uDetW[5,15] >= -1
 udetmc0A: udetXa[6,3] - udetWa[5,15] <= 0
 udetmc0A: udetXa[6,3] - udetXa[6,15] <= 0
 udetmc0A: udetXa[6,3] - udetXa[6,15] - udetWa[5,15] >= -1
 udetMC1: uDetX[6,7] - uDetW[5,3] = 0
 udetmc1A: udetXa[6,7] - udetWa[5,3] = 0
 udetMC2: uDetX[6,11] - uDetW[5,7] <= 0
 udetMC2: uDetX[6,11] - uDetW[5,11] <= 0
 udetMC2: uDetX[6,11] - uDetW[5,7] - uDetW[5,11] >= -1
 udetmc2A: udetXa[6,11] - udetWa[5,7] <= 0
 udetmc2A: udetXa[6,11] - udetWa[5,11] <= 0
 udetmc2A: udetXa[6,11] - udetWa[5,7] - udetWa[5,11] >= -1
 udetMC3: uDetX[6,15] - uDetW[5,3] <= 0
 udetMC3: uDetX[6,15] - uDetW[5,11] <= 0
 udetMC3: uDetX[6,15] - uDetW[5,3] - uDetW[5,11] >= -1
 udetmc3A: udetXa[6,15] - udetWa[5,3] <= 0
 udetmc3A: udetXa[6,15] - udetWa[5,11] <= 0
 udetmc3A: udetXa[6,15] - udetWa[5,3] - udetWa[5,11] >= -1
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
 ufrMC0: - uDetX[5,0] + uFrMC[5,0] <= 0
 ufrMC0: - udx[5,12] + uFrMC[5,0] <= 0
 ufrMC0: - udw[4,12] + uFrMC[5,0] <= 0
 ufrMC0: udx[5,0] + uFrMC[5,0] <= 1
 ufrMC0: udx[5,0] - udx[5,12] - udw[4,12] - uDetX[5,15] + uFrMC[5,0] >= -2
 ufrmc0A: - udetXa[5,0] + ufrMCa[5,0] <= 0
 ufrmc0A: - udx[5,12] + ufrMCa[5,0] <= 0
 ufrmc0A: - udw[4,12] + ufrMCa[5,0] <= 0
 ufrmc0A: udx[5,0] + ufrMCa[5,0] <= 1
 ufrmc0A: udx[5,0] - udx[5,12] - udw[4,12] - udetXa[5,15] + ufrMCa[5,0]
   >= -2
 ufrMC1: uFrMC[5,4] = 0
 ufrmc1A: ufrMCa[5,4] = 0
 ufrMC2: - uDetX[5,8] + uFrMC[5,8] <= 0
 ufrMC2: - udw[4,4] + uFrMC[5,8] <= 0
 ufrMC2: - udw[4,8] + uFrMC[5,8] <= 0
 ufrMC2: udx[5,8] + uFrMC[5,8] <= 1
 ufrMC2: udx[5,8] - udw[4,4] - udw[4,8] - uDetX[5,8] + uFrMC[5,8] >= -2
 ufrmc2A: - udetXa[5,8] + ufrMCa[5,8] <= 0
 ufrmc2A: - udw[4,4] + ufrMCa[5,8] <= 0
 ufrmc2A: - udw[4,8] + ufrMCa[5,8] <= 0
 ufrmc2A: udx[5,8] + ufrMCa[5,8] <= 1
 ufrmc2A: udx[5,8] - udw[4,4] - udw[4,8] - udetXa[5,8] + ufrMCa[5,8] >= -2
 ufrMC3: - uDetX[5,12] + uFrMC[5,12] <= 0
 ufrMC3: - udw[4,0] + uFrMC[5,12] <= 0
 ufrMC3: - udw[4,8] + uFrMC[5,12] <= 0
 ufrMC3: udx[5,12] + uFrMC[5,12] <= 1
 ufrMC3: udx[5,8] - udw[4,0] - udw[4,8] - uDetX[5,12] + uFrMC[5,12] >= -2
 ufrmc3A: - udetXa[5,12] + ufrMCa[5,12] <= 0
 ufrmc3A: - udw[4,0] + ufrMCa[5,12] <= 0
 ufrmc3A: - udw[4,8] + ufrMCa[5,12] <= 0
 ufrmc3A: udx[5,12] + ufrMCa[5,12] <= 1
 ufrmc3A: udx[5,8] - udw[4,0] - udw[4,8] - udetXa[5,12] + ufrMCa[5,12]
   >= -2
 ufrMC0: - uDetX[5,1] + uFrMC[5,1] <= 0
 ufrMC0: - udx[5,13] + uFrMC[5,1] <= 0
 ufrMC0: - udw[4,13] + uFrMC[5,1] <= 0
 ufrMC0: udx[5,1] + uFrMC[5,1] <= 1
 ufrMC0: udx[5,1] - udx[5,13] - udw[4,13] - uDetX[5,15] + uFrMC[5,1] >= -2
 ufrmc0A: - udetXa[5,1] + ufrMCa[5,1] <= 0
 ufrmc0A: - udx[5,13] + ufrMCa[5,1] <= 0
 ufrmc0A: - udw[4,13] + ufrMCa[5,1] <= 0
 ufrmc0A: udx[5,1] + ufrMCa[5,1] <= 1
 ufrmc0A: udx[5,1] - udx[5,13] - udw[4,13] - udetXa[5,15] + ufrMCa[5,1]
   >= -2
 ufrMC1: uFrMC[5,5] = 0
 ufrmc1A: ufrMCa[5,5] = 0
 ufrMC2: - uDetX[5,9] + uFrMC[5,9] <= 0
 ufrMC2: - udw[4,5] + uFrMC[5,9] <= 0
 ufrMC2: - udw[4,9] + uFrMC[5,9] <= 0
 ufrMC2: udx[5,9] + uFrMC[5,9] <= 1
 ufrMC2: udx[5,9] - udw[4,5] - udw[4,9] - uDetX[5,9] + uFrMC[5,9] >= -2
 ufrmc2A: - udetXa[5,9] + ufrMCa[5,9] <= 0
 ufrmc2A: - udw[4,5] + ufrMCa[5,9] <= 0
 ufrmc2A: - udw[4,9] + ufrMCa[5,9] <= 0
 ufrmc2A: udx[5,9] + ufrMCa[5,9] <= 1
 ufrmc2A: udx[5,9] - udw[4,5] - udw[4,9] - udetXa[5,9] + ufrMCa[5,9] >= -2
 ufrMC3: - uDetX[5,13] + uFrMC[5,13] <= 0
 ufrMC3: - udw[4,1] + uFrMC[5,13] <= 0
 ufrMC3: - udw[4,9] + uFrMC[5,13] <= 0
 ufrMC3: udx[5,13] + uFrMC[5,13] <= 1
 ufrMC3: udx[5,9] - udw[4,1] - udw[4,9] - uDetX[5,13] + uFrMC[5,13] >= -2
 ufrmc3A: - udetXa[5,13] + ufrMCa[5,13] <= 0
 ufrmc3A: - udw[4,1] + ufrMCa[5,13] <= 0
 ufrmc3A: - udw[4,9] + ufrMCa[5,13] <= 0
 ufrmc3A: udx[5,13] + ufrMCa[5,13] <= 1
 ufrmc3A: udx[5,9] - udw[4,1] - udw[4,9] - udetXa[5,13] + ufrMCa[5,13]
   >= -2
 ufrMC0: - uDetX[5,2] + uFrMC[5,2] <= 0
 ufrMC0: - udx[5,14] + uFrMC[5,2] <= 0
 ufrMC0: - udw[4,14] + uFrMC[5,2] <= 0
 ufrMC0: udx[5,2] + uFrMC[5,2] <= 1
 ufrMC0: udx[5,2] - udx[5,14] - udw[4,14] - uDetX[5,15] + uFrMC[5,2] >= -2
 ufrmc0A: - udetXa[5,2] + ufrMCa[5,2] <= 0
 ufrmc0A: - udx[5,14] + ufrMCa[5,2] <= 0
 ufrmc0A: - udw[4,14] + ufrMCa[5,2] <= 0
 ufrmc0A: udx[5,2] + ufrMCa[5,2] <= 1
 ufrmc0A: udx[5,2] - udx[5,14] - udw[4,14] - udetXa[5,15] + ufrMCa[5,2]
   >= -2
 ufrMC1: uFrMC[5,6] = 0
 ufrmc1A: ufrMCa[5,6] = 0
 ufrMC2: - uDetX[5,10] + uFrMC[5,10] <= 0
 ufrMC2: - udw[4,6] + uFrMC[5,10] <= 0
 ufrMC2: - udw[4,10] + uFrMC[5,10] <= 0
 ufrMC2: udx[5,10] + uFrMC[5,10] <= 1
 ufrMC2: udx[5,10] - udw[4,6] - udw[4,10] - uDetX[5,10] + uFrMC[5,10]
   >= -2
 ufrmc2A: - udetXa[5,10] + ufrMCa[5,10] <= 0
 ufrmc2A: - udw[4,6] + ufrMCa[5,10] <= 0
 ufrmc2A: - udw[4,10] + ufrMCa[5,10] <= 0
 ufrmc2A: udx[5,10] + ufrMCa[5,10] <= 1
 ufrmc2A: udx[5,10] - udw[4,6] - udw[4,10] - udetXa[5,10] + ufrMCa[5,10]
   >= -2
 ufrMC3: - uDetX[5,14] + uFrMC[5,14] <= 0
 ufrMC3: - udw[4,2] + uFrMC[5,14] <= 0
 ufrMC3: - udw[4,10] + uFrMC[5,14] <= 0
 ufrMC3: udx[5,14] + uFrMC[5,14] <= 1
 ufrMC3: udx[5,10] - udw[4,2] - udw[4,10] - uDetX[5,14] + uFrMC[5,14]
   >= -2
 ufrmc3A: - udetXa[5,14] + ufrMCa[5,14] <= 0
 ufrmc3A: - udw[4,2] + ufrMCa[5,14] <= 0
 ufrmc3A: - udw[4,10] + ufrMCa[5,14] <= 0
 ufrmc3A: udx[5,14] + ufrMCa[5,14] <= 1
 ufrmc3A: udx[5,10] - udw[4,2] - udw[4,10] - udetXa[5,14] + ufrMCa[5,14]
   >= -2
 ufrMC0: - uDetX[5,3] + uFrMC[5,3] <= 0
 ufrMC0: - udx[5,15] + uFrMC[5,3] <= 0
 ufrMC0: - udw[4,15] + uFrMC[5,3] <= 0
 ufrMC0: udx[5,3] + uFrMC[5,3] <= 1
 ufrMC0: udx[5,3] - udx[5,15] - udw[4,15] - uDetX[5,15] + uFrMC[5,3] >= -2
 ufrmc0A: - udetXa[5,3] + ufrMCa[5,3] <= 0
 ufrmc0A: - udx[5,15] + ufrMCa[5,3] <= 0
 ufrmc0A: - udw[4,15] + ufrMCa[5,3] <= 0
 ufrmc0A: udx[5,3] + ufrMCa[5,3] <= 1
 ufrmc0A: udx[5,3] - udx[5,15] - udw[4,15] - udetXa[5,15] + ufrMCa[5,3]
   >= -2
 ufrMC1: uFrMC[5,7] = 0
 ufrmc1A: ufrMCa[5,7] = 0
 ufrMC2: - uDetX[5,11] + uFrMC[5,11] <= 0
 ufrMC2: - udw[4,7] + uFrMC[5,11] <= 0
 ufrMC2: - udw[4,11] + uFrMC[5,11] <= 0
 ufrMC2: udx[5,11] + uFrMC[5,11] <= 1
 ufrMC2: udx[5,11] - udw[4,7] - udw[4,11] - uDetX[5,11] + uFrMC[5,11]
   >= -2
 ufrmc2A: - udetXa[5,11] + ufrMCa[5,11] <= 0
 ufrmc2A: - udw[4,7] + ufrMCa[5,11] <= 0
 ufrmc2A: - udw[4,11] + ufrMCa[5,11] <= 0
 ufrmc2A: udx[5,11] + ufrMCa[5,11] <= 1
 ufrmc2A: udx[5,11] - udw[4,7] - udw[4,11] - udetXa[5,11] + ufrMCa[5,11]
   >= -2
 ufrMC3: - uDetX[5,15] + uFrMC[5,15] <= 0
 ufrMC3: - udw[4,3] + uFrMC[5,15] <= 0
 ufrMC3: - udw[4,11] + uFrMC[5,15] <= 0
 ufrMC3: udx[5,15] + uFrMC[5,15] <= 1
 ufrMC3: udx[5,11] - udw[4,3] - udw[4,11] - uDetX[5,15] + uFrMC[5,15]
   >= -2
 ufrmc3A: - udetXa[5,15] + ufrMCa[5,15] <= 0
 ufrmc3A: - udw[4,3] + ufrMCa[5,15] <= 0
 ufrmc3A: - udw[4,11] + ufrMCa[5,15] <= 0
 ufrmc3A: udx[5,15] + ufrMCa[5,15] <= 1
 ufrmc3A: udx[5,11] - udw[4,3] - udw[4,11] - udetXa[5,15] + ufrMCa[5,15]
   >= -2
 ufrMC0: - uDetX[6,0] + uFrMC[6,0] <= 0
 ufrMC0: - udx[6,12] + uFrMC[6,0] <= 0
 ufrMC0: - udw[5,12] + uFrMC[6,0] <= 0
 ufrMC0: udx[6,0] + uFrMC[6,0] <= 1
 ufrMC0: udx[6,0] - udx[6,12] - udw[5,12] - uDetX[6,15] + uFrMC[6,0] >= -2
 ufrmc0A: - udetXa[6,0] + ufrMCa[6,0] <= 0
 ufrmc0A: - udx[6,12] + ufrMCa[6,0] <= 0
 ufrmc0A: - udw[5,12] + ufrMCa[6,0] <= 0
 ufrmc0A: udx[6,0] + ufrMCa[6,0] <= 1
 ufrmc0A: udx[6,0] - udx[6,12] - udw[5,12] - udetXa[6,15] + ufrMCa[6,0]
   >= -2
 ufrMC1: uFrMC[6,4] = 0
 ufrmc1A: ufrMCa[6,4] = 0
 ufrMC2: - uDetX[6,8] + uFrMC[6,8] <= 0
 ufrMC2: - udw[5,4] + uFrMC[6,8] <= 0
 ufrMC2: - udw[5,8] + uFrMC[6,8] <= 0
 ufrMC2: udx[6,8] + uFrMC[6,8] <= 1
 ufrMC2: udx[6,8] - udw[5,4] - udw[5,8] - uDetX[6,8] + uFrMC[6,8] >= -2
 ufrmc2A: - udetXa[6,8] + ufrMCa[6,8] <= 0
 ufrmc2A: - udw[5,4] + ufrMCa[6,8] <= 0
 ufrmc2A: - udw[5,8] + ufrMCa[6,8] <= 0
 ufrmc2A: udx[6,8] + ufrMCa[6,8] <= 1
 ufrmc2A: udx[6,8] - udw[5,4] - udw[5,8] - udetXa[6,8] + ufrMCa[6,8] >= -2
 ufrMC3: - uDetX[6,12] + uFrMC[6,12] <= 0
 ufrMC3: - udw[5,0] + uFrMC[6,12] <= 0
 ufrMC3: - udw[5,8] + uFrMC[6,12] <= 0
 ufrMC3: udx[6,12] + uFrMC[6,12] <= 1
 ufrMC3: udx[6,8] - udw[5,0] - udw[5,8] - uDetX[6,12] + uFrMC[6,12] >= -2
 ufrmc3A: - udetXa[6,12] + ufrMCa[6,12] <= 0
 ufrmc3A: - udw[5,0] + ufrMCa[6,12] <= 0
 ufrmc3A: - udw[5,8] + ufrMCa[6,12] <= 0
 ufrmc3A: udx[6,12] + ufrMCa[6,12] <= 1
 ufrmc3A: udx[6,8] - udw[5,0] - udw[5,8] - udetXa[6,12] + ufrMCa[6,12]
   >= -2
 ufrMC0: - uDetX[6,1] + uFrMC[6,1] <= 0
 ufrMC0: - udx[6,13] + uFrMC[6,1] <= 0
 ufrMC0: - udw[5,13] + uFrMC[6,1] <= 0
 ufrMC0: udx[6,1] + uFrMC[6,1] <= 1
 ufrMC0: udx[6,1] - udx[6,13] - udw[5,13] - uDetX[6,15] + uFrMC[6,1] >= -2
 ufrmc0A: - udetXa[6,1] + ufrMCa[6,1] <= 0
 ufrmc0A: - udx[6,13] + ufrMCa[6,1] <= 0
 ufrmc0A: - udw[5,13] + ufrMCa[6,1] <= 0
 ufrmc0A: udx[6,1] + ufrMCa[6,1] <= 1
 ufrmc0A: udx[6,1] - udx[6,13] - udw[5,13] - udetXa[6,15] + ufrMCa[6,1]
   >= -2
 ufrMC1: uFrMC[6,5] = 0
 ufrmc1A: ufrMCa[6,5] = 0
 ufrMC2: - uDetX[6,9] + uFrMC[6,9] <= 0
 ufrMC2: - udw[5,5] + uFrMC[6,9] <= 0
 ufrMC2: - udw[5,9] + uFrMC[6,9] <= 0
 ufrMC2: udx[6,9] + uFrMC[6,9] <= 1
 ufrMC2: udx[6,9] - udw[5,5] - udw[5,9] - uDetX[6,9] + uFrMC[6,9] >= -2
 ufrmc2A: - udetXa[6,9] + ufrMCa[6,9] <= 0
 ufrmc2A: - udw[5,5] + ufrMCa[6,9] <= 0
 ufrmc2A: - udw[5,9] + ufrMCa[6,9] <= 0
 ufrmc2A: udx[6,9] + ufrMCa[6,9] <= 1
 ufrmc2A: udx[6,9] - udw[5,5] - udw[5,9] - udetXa[6,9] + ufrMCa[6,9] >= -2
 ufrMC3: - uDetX[6,13] + uFrMC[6,13] <= 0
 ufrMC3: - udw[5,1] + uFrMC[6,13] <= 0
 ufrMC3: - udw[5,9] + uFrMC[6,13] <= 0
 ufrMC3: udx[6,13] + uFrMC[6,13] <= 1
 ufrMC3: udx[6,9] - udw[5,1] - udw[5,9] - uDetX[6,13] + uFrMC[6,13] >= -2
 ufrmc3A: - udetXa[6,13] + ufrMCa[6,13] <= 0
 ufrmc3A: - udw[5,1] + ufrMCa[6,13] <= 0
 ufrmc3A: - udw[5,9] + ufrMCa[6,13] <= 0
 ufrmc3A: udx[6,13] + ufrMCa[6,13] <= 1
 ufrmc3A: udx[6,9] - udw[5,1] - udw[5,9] - udetXa[6,13] + ufrMCa[6,13]
   >= -2
 ufrMC0: - uDetX[6,2] + uFrMC[6,2] <= 0
 ufrMC0: - udx[6,14] + uFrMC[6,2] <= 0
 ufrMC0: - udw[5,14] + uFrMC[6,2] <= 0
 ufrMC0: udx[6,2] + uFrMC[6,2] <= 1
 ufrMC0: udx[6,2] - udx[6,14] - udw[5,14] - uDetX[6,15] + uFrMC[6,2] >= -2
 ufrmc0A: - udetXa[6,2] + ufrMCa[6,2] <= 0
 ufrmc0A: - udx[6,14] + ufrMCa[6,2] <= 0
 ufrmc0A: - udw[5,14] + ufrMCa[6,2] <= 0
 ufrmc0A: udx[6,2] + ufrMCa[6,2] <= 1
 ufrmc0A: udx[6,2] - udx[6,14] - udw[5,14] - udetXa[6,15] + ufrMCa[6,2]
   >= -2
 ufrMC1: uFrMC[6,6] = 0
 ufrmc1A: ufrMCa[6,6] = 0
 ufrMC2: - uDetX[6,10] + uFrMC[6,10] <= 0
 ufrMC2: - udw[5,6] + uFrMC[6,10] <= 0
 ufrMC2: - udw[5,10] + uFrMC[6,10] <= 0
 ufrMC2: udx[6,10] + uFrMC[6,10] <= 1
 ufrMC2: udx[6,10] - udw[5,6] - udw[5,10] - uDetX[6,10] + uFrMC[6,10]
   >= -2
 ufrmc2A: - udetXa[6,10] + ufrMCa[6,10] <= 0
 ufrmc2A: - udw[5,6] + ufrMCa[6,10] <= 0
 ufrmc2A: - udw[5,10] + ufrMCa[6,10] <= 0
 ufrmc2A: udx[6,10] + ufrMCa[6,10] <= 1
 ufrmc2A: udx[6,10] - udw[5,6] - udw[5,10] - udetXa[6,10] + ufrMCa[6,10]
   >= -2
 ufrMC3: - uDetX[6,14] + uFrMC[6,14] <= 0
 ufrMC3: - udw[5,2] + uFrMC[6,14] <= 0
 ufrMC3: - udw[5,10] + uFrMC[6,14] <= 0
 ufrMC3: udx[6,14] + uFrMC[6,14] <= 1
 ufrMC3: udx[6,10] - udw[5,2] - udw[5,10] - uDetX[6,14] + uFrMC[6,14]
   >= -2
 ufrmc3A: - udetXa[6,14] + ufrMCa[6,14] <= 0
 ufrmc3A: - udw[5,2] + ufrMCa[6,14] <= 0
 ufrmc3A: - udw[5,10] + ufrMCa[6,14] <= 0
 ufrmc3A: udx[6,14] + ufrMCa[6,14] <= 1
 ufrmc3A: udx[6,10] - udw[5,2] - udw[5,10] - udetXa[6,14] + ufrMCa[6,14]
   >= -2
 ufrMC0: - uDetX[6,3] + uFrMC[6,3] <= 0
 ufrMC0: - udx[6,15] + uFrMC[6,3] <= 0
 ufrMC0: - udw[5,15] + uFrMC[6,3] <= 0
 ufrMC0: udx[6,3] + uFrMC[6,3] <= 1
 ufrMC0: udx[6,3] - udx[6,15] - udw[5,15] - uDetX[6,15] + uFrMC[6,3] >= -2
 ufrmc0A: - udetXa[6,3] + ufrMCa[6,3] <= 0
 ufrmc0A: - udx[6,15] + ufrMCa[6,3] <= 0
 ufrmc0A: - udw[5,15] + ufrMCa[6,3] <= 0
 ufrmc0A: udx[6,3] + ufrMCa[6,3] <= 1
 ufrmc0A: udx[6,3] - udx[6,15] - udw[5,15] - udetXa[6,15] + ufrMCa[6,3]
   >= -2
 ufrMC1: uFrMC[6,7] = 0
 ufrmc1A: ufrMCa[6,7] = 0
 ufrMC2: - uDetX[6,11] + uFrMC[6,11] <= 0
 ufrMC2: - udw[5,7] + uFrMC[6,11] <= 0
 ufrMC2: - udw[5,11] + uFrMC[6,11] <= 0
 ufrMC2: udx[6,11] + uFrMC[6,11] <= 1
 ufrMC2: udx[6,11] - udw[5,7] - udw[5,11] - uDetX[6,11] + uFrMC[6,11]
   >= -2
 ufrmc2A: - udetXa[6,11] + ufrMCa[6,11] <= 0
 ufrmc2A: - udw[5,7] + ufrMCa[6,11] <= 0
 ufrmc2A: - udw[5,11] + ufrMCa[6,11] <= 0
 ufrmc2A: udx[6,11] + ufrMCa[6,11] <= 1
 ufrmc2A: udx[6,11] - udw[5,7] - udw[5,11] - udetXa[6,11] + ufrMCa[6,11]
   >= -2
 ufrMC3: - uDetX[6,15] + uFrMC[6,15] <= 0
 ufrMC3: - udw[5,3] + uFrMC[6,15] <= 0
 ufrMC3: - udw[5,11] + uFrMC[6,15] <= 0
 ufrMC3: udx[6,15] + uFrMC[6,15] <= 1
 ufrMC3: udx[6,11] - udw[5,3] - udw[5,11] - uDetX[6,15] + uFrMC[6,15]
   >= -2
 ufrmc3A: - udetXa[6,15] + ufrMCa[6,15] <= 0
 ufrmc3A: - udw[5,3] + ufrMCa[6,15] <= 0
 ufrmc3A: - udw[5,11] + ufrMCa[6,15] <= 0
 ufrmc3A: udx[6,15] + ufrMCa[6,15] <= 1
 ufrmc3A: udx[6,11] - udw[5,3] - udw[5,11] - udetXa[6,15] + ufrMCa[6,15]
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
   - ufrSCa[4,14] - ufrSCa[4,15] - ufrSCa[5,0] - ufrSCa[5,1] - ufrSCa[5,2]
   - ufrSCa[5,3] - ufrSCa[5,4] - ufrSCa[5,5] - ufrSCa[5,6] - ufrSCa[5,7]
   - ufrSCa[5,8] - ufrSCa[5,9] - ufrSCa[5,10] - ufrSCa[5,11] - ufrSCa[5,12]
   - ufrSCa[5,13] - ufrSCa[5,14] - ufrSCa[5,15] - ufrSCa[6,0] - ufrSCa[6,1]
   - ufrSCa[6,2] - ufrSCa[6,3] - ufrSCa[6,4] - ufrSCa[6,5] - ufrSCa[6,6]
   - ufrSCa[6,7] - ufrSCa[6,8] - ufrSCa[6,9] - ufrSCa[6,10] - ufrSCa[6,11]
   - ufrSCa[6,12] - ufrSCa[6,13] - ufrSCa[6,14] - ufrSCa[6,15]
   - ufrMCa[2,0] - ufrMCa[2,1] - ufrMCa[2,2] - ufrMCa[2,3] - ufrMCa[2,4]
   - ufrMCa[2,5] - ufrMCa[2,6] - ufrMCa[2,7] - ufrMCa[2,8] - ufrMCa[2,9]
   - ufrMCa[2,10] - ufrMCa[2,11] - ufrMCa[2,12] - ufrMCa[2,13]
   - ufrMCa[2,14] - ufrMCa[2,15] - ufrMCa[3,0] - ufrMCa[3,1] - ufrMCa[3,2]
   - ufrMCa[3,3] - ufrMCa[3,4] - ufrMCa[3,5] - ufrMCa[3,6] - ufrMCa[3,7]
   - ufrMCa[3,8] - ufrMCa[3,9] - ufrMCa[3,10] - ufrMCa[3,11] - ufrMCa[3,12]
   - ufrMCa[3,13] - ufrMCa[3,14] - ufrMCa[3,15] - ufrMCa[4,0] - ufrMCa[4,1]
   - ufrMCa[4,2] - ufrMCa[4,3] - ufrMCa[4,4] - ufrMCa[4,5] - ufrMCa[4,6]
   - ufrMCa[4,7] - ufrMCa[4,8] - ufrMCa[4,9] - ufrMCa[4,10] - ufrMCa[4,11]
   - ufrMCa[4,12] - ufrMCa[4,13] - ufrMCa[4,14] - ufrMCa[4,15]
   - ufrMCa[5,0] - ufrMCa[5,1] - ufrMCa[5,2] - ufrMCa[5,3] - ufrMCa[5,4]
   - ufrMCa[5,5] - ufrMCa[5,6] - ufrMCa[5,7] - ufrMCa[5,8] - ufrMCa[5,9]
   - ufrMCa[5,10] - ufrMCa[5,11] - ufrMCa[5,12] - ufrMCa[5,13]
   - ufrMCa[5,14] - ufrMCa[5,15] - ufrMCa[6,0] - ufrMCa[6,1] - ufrMCa[6,2]
   - ufrMCa[6,3] - ufrMCa[6,4] - ufrMCa[6,5] - ufrMCa[6,6] - ufrMCa[6,7]
   - ufrMCa[6,8] - ufrMCa[6,9] - ufrMCa[6,10] - ufrMCa[6,11] - ufrMCa[6,12]
   - ufrMCa[6,13] - ufrMCa[6,14] - ufrMCa[6,15] = 0
 R12323: lGstk[28,8] = 0
 R12324: lGstk[28,9] = 0
 R12325: lGstk[28,10] = 0
 R12326: lGstk[28,11] = 0
 R12327: lGstk[28,12] = 0
 R12328: lGstk[28,13] = 0
 R12329: lGstk[28,14] = 0
 R12330: lGstk[28,15] = 0
 R12331: lGstk[29,8] = 0
 R12332: lGstk[29,9] = 0
 R12333: lGstk[29,10] = 0
 R12334: lGstk[29,11] = 0
 R12335: lGstk[29,12] = 0
 R12336: lGstk[29,13] = 0
 R12337: lGstk[29,14] = 0
 R12338: lGstk[29,15] = 0
 R12339: lGstk[30,8] = 0
 R12340: lGstk[30,9] = 0
 R12341: lGstk[30,10] = 0
 R12342: lGstk[30,11] = 0
 R12343: lGstk[30,12] = 0
 R12344: lGstk[30,13] = 0
 R12345: lGstk[30,14] = 0
 R12346: lGstk[30,15] = 0
 R12347: lGstk[31,8] = 0
 R12348: lGstk[31,9] = 0
 R12349: lGstk[31,10] = 0
 R12350: lGstk[31,11] = 0
 R12351: lGstk[31,12] = 0
 R12352: lGstk[31,13] = 0
 R12353: lGstk[31,14] = 0
 R12354: lGstk[31,15] = 0
 R12355: lGstk[32,8] = 0
 R12356: lGstk[32,9] = 0
 R12357: lGstk[32,10] = 0
 R12358: lGstk[32,11] = 0
 R12359: lGstk[32,12] = 0
 R12360: lGstk[32,13] = 0
 R12361: lGstk[32,14] = 0
 R12362: lGstk[32,15] = 0
 R12363: lgstka[28,8] = 0
 R12364: lgstka[28,9] = 0
 R12365: lgstka[28,10] = 0
 R12366: lgstka[28,11] = 0
 R12367: lgstka[28,12] = 0
 R12368: lgstka[28,13] = 0
 R12369: lgstka[28,14] = 0
 R12370: lgstka[28,15] = 0
 R12371: lgstka[29,8] = 0
 R12372: lgstka[29,9] = 0
 R12373: lgstka[29,10] = 0
 R12374: lgstka[29,11] = 0
 R12375: lgstka[29,12] = 0
 R12376: lgstka[29,13] = 0
 R12377: lgstka[29,14] = 0
 R12378: lgstka[29,15] = 0
 R12379: lgstka[30,8] = 0
 R12380: lgstka[30,9] = 0
 R12381: lgstka[30,10] = 0
 R12382: lgstka[30,11] = 0
 R12383: lgstka[30,12] = 0
 R12384: lgstka[30,13] = 0
 R12385: lgstka[30,14] = 0
 R12386: lgstka[30,15] = 0
 R12387: lgstka[31,8] = 0
 R12388: lgstka[31,9] = 0
 R12389: lgstka[31,10] = 0
 R12390: lgstka[31,11] = 0
 R12391: lgstka[31,12] = 0
 R12392: lgstka[31,13] = 0
 R12393: lgstka[31,14] = 0
 R12394: lgstka[31,15] = 0
 R12395: lgstka[32,8] = 0
 R12396: lgstka[32,9] = 0
 R12397: lgstka[32,10] = 0
 R12398: lgstka[32,11] = 0
 R12399: lgstka[32,12] = 0
 R12400: lgstka[32,13] = 0
 R12401: lgstka[32,14] = 0
 R12402: lgstka[32,15] = 0
 VC: lDetZ[32,0] + lDetZ[32,1] + lDetZ[32,2] + lDetZ[32,3] + lDetZ[32,4]
   + lDetZ[32,5] + lDetZ[32,6] + lDetZ[32,7] + lDetZ[32,8] + lDetZ[32,9]
   + lDetZ[32,10] + lDetZ[32,11] + lDetZ[32,12] + lDetZ[32,13]
   + lDetZ[32,14] + lDetZ[32,15] = 16
 ldetSC: lDetX[28,0] - lDetY[28,0] = 0
 ldetSCa: ldetXa[28,0] - ldetYa[28,0] = 0
 lfrSC: lDX[28,0] - ldx[28,0] - lFrSC[28,0] >= 0
 lfrSC: lDetX[28,0] - lFrSC[28,0] >= 0
 lfrSC: - lDX[28,0] + ldx[28,0] - lDetX[28,0] + lFrSC[28,0] >= -1
 lfrSCa: lDX[28,0] - ldx[28,0] - lfrSCa[28,0] >= 0
 lfrSCa: ldetXa[28,0] - lfrSCa[28,0] >= 0
 lfrSCa: - lDX[28,0] + ldx[28,0] - ldetXa[28,0] + lfrSCa[28,0] >= -1
 ldetSC: lDetX[28,1] - lDetY[28,1] = 0
 ldetSCa: ldetXa[28,1] - ldetYa[28,1] = 0
 lfrSC: lDX[28,1] - ldx[28,1] - lFrSC[28,1] >= 0
 lfrSC: lDetX[28,1] - lFrSC[28,1] >= 0
 lfrSC: - lDX[28,1] + ldx[28,1] - lDetX[28,1] + lFrSC[28,1] >= -1
 lfrSCa: lDX[28,1] - ldx[28,1] - lfrSCa[28,1] >= 0
 lfrSCa: ldetXa[28,1] - lfrSCa[28,1] >= 0
 lfrSCa: - lDX[28,1] + ldx[28,1] - ldetXa[28,1] + lfrSCa[28,1] >= -1
 ldetSC: lDetX[28,2] - lDetY[28,2] = 0
 ldetSCa: ldetXa[28,2] - ldetYa[28,2] = 0
 lfrSC: lDX[28,2] - ldx[28,2] - lFrSC[28,2] >= 0
 lfrSC: lDetX[28,2] - lFrSC[28,2] >= 0
 lfrSC: - lDX[28,2] + ldx[28,2] - lDetX[28,2] + lFrSC[28,2] >= -1
 lfrSCa: lDX[28,2] - ldx[28,2] - lfrSCa[28,2] >= 0
 lfrSCa: ldetXa[28,2] - lfrSCa[28,2] >= 0
 lfrSCa: - lDX[28,2] + ldx[28,2] - ldetXa[28,2] + lfrSCa[28,2] >= -1
 ldetSC: lDetX[28,3] - lDetY[28,3] = 0
 ldetSCa: ldetXa[28,3] - ldetYa[28,3] = 0
 lfrSC: lDX[28,3] - ldx[28,3] - lFrSC[28,3] >= 0
 lfrSC: lDetX[28,3] - lFrSC[28,3] >= 0
 lfrSC: - lDX[28,3] + ldx[28,3] - lDetX[28,3] + lFrSC[28,3] >= -1
 lfrSCa: lDX[28,3] - ldx[28,3] - lfrSCa[28,3] >= 0
 lfrSCa: ldetXa[28,3] - lfrSCa[28,3] >= 0
 lfrSCa: - lDX[28,3] + ldx[28,3] - ldetXa[28,3] + lfrSCa[28,3] >= -1
 ldetSC: lDetX[28,4] - lDetY[28,4] = 0
 ldetSCa: ldetXa[28,4] - ldetYa[28,4] = 0
 lfrSC: lDX[28,4] - ldx[28,4] - lFrSC[28,4] >= 0
 lfrSC: lDetX[28,4] - lFrSC[28,4] >= 0
 lfrSC: - lDX[28,4] + ldx[28,4] - lDetX[28,4] + lFrSC[28,4] >= -1
 lfrSCa: lDX[28,4] - ldx[28,4] - lfrSCa[28,4] >= 0
 lfrSCa: ldetXa[28,4] - lfrSCa[28,4] >= 0
 lfrSCa: - lDX[28,4] + ldx[28,4] - ldetXa[28,4] + lfrSCa[28,4] >= -1
 ldetSC: lDetX[28,5] - lDetY[28,5] = 0
 ldetSCa: ldetXa[28,5] - ldetYa[28,5] = 0
 lfrSC: lDX[28,5] - ldx[28,5] - lFrSC[28,5] >= 0
 lfrSC: lDetX[28,5] - lFrSC[28,5] >= 0
 lfrSC: - lDX[28,5] + ldx[28,5] - lDetX[28,5] + lFrSC[28,5] >= -1
 lfrSCa: lDX[28,5] - ldx[28,5] - lfrSCa[28,5] >= 0
 lfrSCa: ldetXa[28,5] - lfrSCa[28,5] >= 0
 lfrSCa: - lDX[28,5] + ldx[28,5] - ldetXa[28,5] + lfrSCa[28,5] >= -1
 ldetSC: lDetX[28,6] - lDetY[28,6] = 0
 ldetSCa: ldetXa[28,6] - ldetYa[28,6] = 0
 lfrSC: lDX[28,6] - ldx[28,6] - lFrSC[28,6] >= 0
 lfrSC: lDetX[28,6] - lFrSC[28,6] >= 0
 lfrSC: - lDX[28,6] + ldx[28,6] - lDetX[28,6] + lFrSC[28,6] >= -1
 lfrSCa: lDX[28,6] - ldx[28,6] - lfrSCa[28,6] >= 0
 lfrSCa: ldetXa[28,6] - lfrSCa[28,6] >= 0
 lfrSCa: - lDX[28,6] + ldx[28,6] - ldetXa[28,6] + lfrSCa[28,6] >= -1
 ldetSC: lDetX[28,7] - lDetY[28,7] = 0
 ldetSCa: ldetXa[28,7] - ldetYa[28,7] = 0
 lfrSC: lDX[28,7] - ldx[28,7] - lFrSC[28,7] >= 0
 lfrSC: lDetX[28,7] - lFrSC[28,7] >= 0
 lfrSC: - lDX[28,7] + ldx[28,7] - lDetX[28,7] + lFrSC[28,7] >= -1
 lfrSCa: lDX[28,7] - ldx[28,7] - lfrSCa[28,7] >= 0
 lfrSCa: ldetXa[28,7] - lfrSCa[28,7] >= 0
 lfrSCa: - lDX[28,7] + ldx[28,7] - ldetXa[28,7] + lfrSCa[28,7] >= -1
 ldetSC: lDetX[28,8] - lDetY[28,8] = 0
 ldetSCa: ldetXa[28,8] - ldetYa[28,8] = 0
 lfrSC: lDX[28,8] - ldx[28,8] - lFrSC[28,8] >= 0
 lfrSC: lDetX[28,8] - lFrSC[28,8] >= 0
 lfrSC: - lDX[28,8] + ldx[28,8] - lDetX[28,8] + lFrSC[28,8] >= -1
 lfrSCa: lDX[28,8] - ldx[28,8] - lfrSCa[28,8] >= 0
 lfrSCa: ldetXa[28,8] - lfrSCa[28,8] >= 0
 lfrSCa: - lDX[28,8] + ldx[28,8] - ldetXa[28,8] + lfrSCa[28,8] >= -1
 ldetSC: lDetX[28,9] - lDetY[28,9] = 0
 ldetSCa: ldetXa[28,9] - ldetYa[28,9] = 0
 lfrSC: lDX[28,9] - ldx[28,9] - lFrSC[28,9] >= 0
 lfrSC: lDetX[28,9] - lFrSC[28,9] >= 0
 lfrSC: - lDX[28,9] + ldx[28,9] - lDetX[28,9] + lFrSC[28,9] >= -1
 lfrSCa: lDX[28,9] - ldx[28,9] - lfrSCa[28,9] >= 0
 lfrSCa: ldetXa[28,9] - lfrSCa[28,9] >= 0
 lfrSCa: - lDX[28,9] + ldx[28,9] - ldetXa[28,9] + lfrSCa[28,9] >= -1
 ldetSC: lDetX[28,10] - lDetY[28,10] = 0
 ldetSCa: ldetXa[28,10] - ldetYa[28,10] = 0
 lfrSC: lDX[28,10] - ldx[28,10] - lFrSC[28,10] >= 0
 lfrSC: lDetX[28,10] - lFrSC[28,10] >= 0
 lfrSC: - lDX[28,10] + ldx[28,10] - lDetX[28,10] + lFrSC[28,10] >= -1
 lfrSCa: lDX[28,10] - ldx[28,10] - lfrSCa[28,10] >= 0
 lfrSCa: ldetXa[28,10] - lfrSCa[28,10] >= 0
 lfrSCa: - lDX[28,10] + ldx[28,10] - ldetXa[28,10] + lfrSCa[28,10] >= -1
 ldetSC: lDetX[28,11] - lDetY[28,11] = 0
 ldetSCa: ldetXa[28,11] - ldetYa[28,11] = 0
 lfrSC: lDX[28,11] - ldx[28,11] - lFrSC[28,11] >= 0
 lfrSC: lDetX[28,11] - lFrSC[28,11] >= 0
 lfrSC: - lDX[28,11] + ldx[28,11] - lDetX[28,11] + lFrSC[28,11] >= -1
 lfrSCa: lDX[28,11] - ldx[28,11] - lfrSCa[28,11] >= 0
 lfrSCa: ldetXa[28,11] - lfrSCa[28,11] >= 0
 lfrSCa: - lDX[28,11] + ldx[28,11] - ldetXa[28,11] + lfrSCa[28,11] >= -1
 ldetSC: lDetX[28,12] - lDetY[28,12] = 0
 ldetSCa: ldetXa[28,12] - ldetYa[28,12] = 0
 lfrSC: lDX[28,12] - ldx[28,12] - lFrSC[28,12] >= 0
 lfrSC: lDetX[28,12] - lFrSC[28,12] >= 0
 lfrSC: - lDX[28,12] + ldx[28,12] - lDetX[28,12] + lFrSC[28,12] >= -1
 lfrSCa: lDX[28,12] - ldx[28,12] - lfrSCa[28,12] >= 0
 lfrSCa: ldetXa[28,12] - lfrSCa[28,12] >= 0
 lfrSCa: - lDX[28,12] + ldx[28,12] - ldetXa[28,12] + lfrSCa[28,12] >= -1
 ldetSC: lDetX[28,13] - lDetY[28,13] = 0
 ldetSCa: ldetXa[28,13] - ldetYa[28,13] = 0
 lfrSC: lDX[28,13] - ldx[28,13] - lFrSC[28,13] >= 0
 lfrSC: lDetX[28,13] - lFrSC[28,13] >= 0
 lfrSC: - lDX[28,13] + ldx[28,13] - lDetX[28,13] + lFrSC[28,13] >= -1
 lfrSCa: lDX[28,13] - ldx[28,13] - lfrSCa[28,13] >= 0
 lfrSCa: ldetXa[28,13] - lfrSCa[28,13] >= 0
 lfrSCa: - lDX[28,13] + ldx[28,13] - ldetXa[28,13] + lfrSCa[28,13] >= -1
 ldetSC: lDetX[28,14] - lDetY[28,14] = 0
 ldetSCa: ldetXa[28,14] - ldetYa[28,14] = 0
 lfrSC: lDX[28,14] - ldx[28,14] - lFrSC[28,14] >= 0
 lfrSC: lDetX[28,14] - lFrSC[28,14] >= 0
 lfrSC: - lDX[28,14] + ldx[28,14] - lDetX[28,14] + lFrSC[28,14] >= -1
 lfrSCa: lDX[28,14] - ldx[28,14] - lfrSCa[28,14] >= 0
 lfrSCa: ldetXa[28,14] - lfrSCa[28,14] >= 0
 lfrSCa: - lDX[28,14] + ldx[28,14] - ldetXa[28,14] + lfrSCa[28,14] >= -1
 ldetSC: lDetX[28,15] - lDetY[28,15] = 0
 ldetSCa: ldetXa[28,15] - ldetYa[28,15] = 0
 lfrSC: lDX[28,15] - ldx[28,15] - lFrSC[28,15] >= 0
 lfrSC: lDetX[28,15] - lFrSC[28,15] >= 0
 lfrSC: - lDX[28,15] + ldx[28,15] - lDetX[28,15] + lFrSC[28,15] >= -1
 lfrSCa: lDX[28,15] - ldx[28,15] - lfrSCa[28,15] >= 0
 lfrSCa: ldetXa[28,15] - lfrSCa[28,15] >= 0
 lfrSCa: - lDX[28,15] + ldx[28,15] - ldetXa[28,15] + lfrSCa[28,15] >= -1
 ldetSC: lDetX[29,0] - lDetY[29,0] = 0
 ldetSCa: ldetXa[29,0] - ldetYa[29,0] = 0
 lfrSC: lDX[29,0] - ldx[29,0] - lFrSC[29,0] >= 0
 lfrSC: lDetX[29,0] - lFrSC[29,0] >= 0
 lfrSC: - lDX[29,0] + ldx[29,0] - lDetX[29,0] + lFrSC[29,0] >= -1
 lfrSCa: lDX[29,0] - ldx[29,0] - lfrSCa[29,0] >= 0
 lfrSCa: ldetXa[29,0] - lfrSCa[29,0] >= 0
 lfrSCa: - lDX[29,0] + ldx[29,0] - ldetXa[29,0] + lfrSCa[29,0] >= -1
 ldetSC: lDetX[29,1] - lDetY[29,1] = 0
 ldetSCa: ldetXa[29,1] - ldetYa[29,1] = 0
 lfrSC: lDX[29,1] - ldx[29,1] - lFrSC[29,1] >= 0
 lfrSC: lDetX[29,1] - lFrSC[29,1] >= 0
 lfrSC: - lDX[29,1] + ldx[29,1] - lDetX[29,1] + lFrSC[29,1] >= -1
 lfrSCa: lDX[29,1] - ldx[29,1] - lfrSCa[29,1] >= 0
 lfrSCa: ldetXa[29,1] - lfrSCa[29,1] >= 0
 lfrSCa: - lDX[29,1] + ldx[29,1] - ldetXa[29,1] + lfrSCa[29,1] >= -1
 ldetSC: lDetX[29,2] - lDetY[29,2] = 0
 ldetSCa: ldetXa[29,2] - ldetYa[29,2] = 0
 lfrSC: lDX[29,2] - ldx[29,2] - lFrSC[29,2] >= 0
 lfrSC: lDetX[29,2] - lFrSC[29,2] >= 0
 lfrSC: - lDX[29,2] + ldx[29,2] - lDetX[29,2] + lFrSC[29,2] >= -1
 lfrSCa: lDX[29,2] - ldx[29,2] - lfrSCa[29,2] >= 0
 lfrSCa: ldetXa[29,2] - lfrSCa[29,2] >= 0
 lfrSCa: - lDX[29,2] + ldx[29,2] - ldetXa[29,2] + lfrSCa[29,2] >= -1
 ldetSC: lDetX[29,3] - lDetY[29,3] = 0
 ldetSCa: ldetXa[29,3] - ldetYa[29,3] = 0
 lfrSC: lDX[29,3] - ldx[29,3] - lFrSC[29,3] >= 0
 lfrSC: lDetX[29,3] - lFrSC[29,3] >= 0
 lfrSC: - lDX[29,3] + ldx[29,3] - lDetX[29,3] + lFrSC[29,3] >= -1
 lfrSCa: lDX[29,3] - ldx[29,3] - lfrSCa[29,3] >= 0
 lfrSCa: ldetXa[29,3] - lfrSCa[29,3] >= 0
 lfrSCa: - lDX[29,3] + ldx[29,3] - ldetXa[29,3] + lfrSCa[29,3] >= -1
 ldetSC: lDetX[29,4] - lDetY[29,4] = 0
 ldetSCa: ldetXa[29,4] - ldetYa[29,4] = 0
 lfrSC: lDX[29,4] - ldx[29,4] - lFrSC[29,4] >= 0
 lfrSC: lDetX[29,4] - lFrSC[29,4] >= 0
 lfrSC: - lDX[29,4] + ldx[29,4] - lDetX[29,4] + lFrSC[29,4] >= -1
 lfrSCa: lDX[29,4] - ldx[29,4] - lfrSCa[29,4] >= 0
 lfrSCa: ldetXa[29,4] - lfrSCa[29,4] >= 0
 lfrSCa: - lDX[29,4] + ldx[29,4] - ldetXa[29,4] + lfrSCa[29,4] >= -1
 ldetSC: lDetX[29,5] - lDetY[29,5] = 0
 ldetSCa: ldetXa[29,5] - ldetYa[29,5] = 0
 lfrSC: lDX[29,5] - ldx[29,5] - lFrSC[29,5] >= 0
 lfrSC: lDetX[29,5] - lFrSC[29,5] >= 0
 lfrSC: - lDX[29,5] + ldx[29,5] - lDetX[29,5] + lFrSC[29,5] >= -1
 lfrSCa: lDX[29,5] - ldx[29,5] - lfrSCa[29,5] >= 0
 lfrSCa: ldetXa[29,5] - lfrSCa[29,5] >= 0
 lfrSCa: - lDX[29,5] + ldx[29,5] - ldetXa[29,5] + lfrSCa[29,5] >= -1
 ldetSC: lDetX[29,6] - lDetY[29,6] = 0
 ldetSCa: ldetXa[29,6] - ldetYa[29,6] = 0
 lfrSC: lDX[29,6] - ldx[29,6] - lFrSC[29,6] >= 0
 lfrSC: lDetX[29,6] - lFrSC[29,6] >= 0
 lfrSC: - lDX[29,6] + ldx[29,6] - lDetX[29,6] + lFrSC[29,6] >= -1
 lfrSCa: lDX[29,6] - ldx[29,6] - lfrSCa[29,6] >= 0
 lfrSCa: ldetXa[29,6] - lfrSCa[29,6] >= 0
 lfrSCa: - lDX[29,6] + ldx[29,6] - ldetXa[29,6] + lfrSCa[29,6] >= -1
 ldetSC: lDetX[29,7] - lDetY[29,7] = 0
 ldetSCa: ldetXa[29,7] - ldetYa[29,7] = 0
 lfrSC: lDX[29,7] - ldx[29,7] - lFrSC[29,7] >= 0
 lfrSC: lDetX[29,7] - lFrSC[29,7] >= 0
 lfrSC: - lDX[29,7] + ldx[29,7] - lDetX[29,7] + lFrSC[29,7] >= -1
 lfrSCa: lDX[29,7] - ldx[29,7] - lfrSCa[29,7] >= 0
 lfrSCa: ldetXa[29,7] - lfrSCa[29,7] >= 0
 lfrSCa: - lDX[29,7] + ldx[29,7] - ldetXa[29,7] + lfrSCa[29,7] >= -1
 ldetSC: lDetX[29,8] - lDetY[29,8] = 0
 ldetSCa: ldetXa[29,8] - ldetYa[29,8] = 0
 lfrSC: lDX[29,8] - ldx[29,8] - lFrSC[29,8] >= 0
 lfrSC: lDetX[29,8] - lFrSC[29,8] >= 0
 lfrSC: - lDX[29,8] + ldx[29,8] - lDetX[29,8] + lFrSC[29,8] >= -1
 lfrSCa: lDX[29,8] - ldx[29,8] - lfrSCa[29,8] >= 0
 lfrSCa: ldetXa[29,8] - lfrSCa[29,8] >= 0
 lfrSCa: - lDX[29,8] + ldx[29,8] - ldetXa[29,8] + lfrSCa[29,8] >= -1
 ldetSC: lDetX[29,9] - lDetY[29,9] = 0
 ldetSCa: ldetXa[29,9] - ldetYa[29,9] = 0
 lfrSC: lDX[29,9] - ldx[29,9] - lFrSC[29,9] >= 0
 lfrSC: lDetX[29,9] - lFrSC[29,9] >= 0
 lfrSC: - lDX[29,9] + ldx[29,9] - lDetX[29,9] + lFrSC[29,9] >= -1
 lfrSCa: lDX[29,9] - ldx[29,9] - lfrSCa[29,9] >= 0
 lfrSCa: ldetXa[29,9] - lfrSCa[29,9] >= 0
 lfrSCa: - lDX[29,9] + ldx[29,9] - ldetXa[29,9] + lfrSCa[29,9] >= -1
 ldetSC: lDetX[29,10] - lDetY[29,10] = 0
 ldetSCa: ldetXa[29,10] - ldetYa[29,10] = 0
 lfrSC: lDX[29,10] - ldx[29,10] - lFrSC[29,10] >= 0
 lfrSC: lDetX[29,10] - lFrSC[29,10] >= 0
 lfrSC: - lDX[29,10] + ldx[29,10] - lDetX[29,10] + lFrSC[29,10] >= -1
 lfrSCa: lDX[29,10] - ldx[29,10] - lfrSCa[29,10] >= 0
 lfrSCa: ldetXa[29,10] - lfrSCa[29,10] >= 0
 lfrSCa: - lDX[29,10] + ldx[29,10] - ldetXa[29,10] + lfrSCa[29,10] >= -1
 ldetSC: lDetX[29,11] - lDetY[29,11] = 0
 ldetSCa: ldetXa[29,11] - ldetYa[29,11] = 0
 lfrSC: lDX[29,11] - ldx[29,11] - lFrSC[29,11] >= 0
 lfrSC: lDetX[29,11] - lFrSC[29,11] >= 0
 lfrSC: - lDX[29,11] + ldx[29,11] - lDetX[29,11] + lFrSC[29,11] >= -1
 lfrSCa: lDX[29,11] - ldx[29,11] - lfrSCa[29,11] >= 0
 lfrSCa: ldetXa[29,11] - lfrSCa[29,11] >= 0
 lfrSCa: - lDX[29,11] + ldx[29,11] - ldetXa[29,11] + lfrSCa[29,11] >= -1
 ldetSC: lDetX[29,12] - lDetY[29,12] = 0
 ldetSCa: ldetXa[29,12] - ldetYa[29,12] = 0
 lfrSC: lDX[29,12] - ldx[29,12] - lFrSC[29,12] >= 0
 lfrSC: lDetX[29,12] - lFrSC[29,12] >= 0
 lfrSC: - lDX[29,12] + ldx[29,12] - lDetX[29,12] + lFrSC[29,12] >= -1
 lfrSCa: lDX[29,12] - ldx[29,12] - lfrSCa[29,12] >= 0
 lfrSCa: ldetXa[29,12] - lfrSCa[29,12] >= 0
 lfrSCa: - lDX[29,12] + ldx[29,12] - ldetXa[29,12] + lfrSCa[29,12] >= -1
 ldetSC: lDetX[29,13] - lDetY[29,13] = 0
 ldetSCa: ldetXa[29,13] - ldetYa[29,13] = 0
 lfrSC: lDX[29,13] - ldx[29,13] - lFrSC[29,13] >= 0
 lfrSC: lDetX[29,13] - lFrSC[29,13] >= 0
 lfrSC: - lDX[29,13] + ldx[29,13] - lDetX[29,13] + lFrSC[29,13] >= -1
 lfrSCa: lDX[29,13] - ldx[29,13] - lfrSCa[29,13] >= 0
 lfrSCa: ldetXa[29,13] - lfrSCa[29,13] >= 0
 lfrSCa: - lDX[29,13] + ldx[29,13] - ldetXa[29,13] + lfrSCa[29,13] >= -1
 ldetSC: lDetX[29,14] - lDetY[29,14] = 0
 ldetSCa: ldetXa[29,14] - ldetYa[29,14] = 0
 lfrSC: lDX[29,14] - ldx[29,14] - lFrSC[29,14] >= 0
 lfrSC: lDetX[29,14] - lFrSC[29,14] >= 0
 lfrSC: - lDX[29,14] + ldx[29,14] - lDetX[29,14] + lFrSC[29,14] >= -1
 lfrSCa: lDX[29,14] - ldx[29,14] - lfrSCa[29,14] >= 0
 lfrSCa: ldetXa[29,14] - lfrSCa[29,14] >= 0
 lfrSCa: - lDX[29,14] + ldx[29,14] - ldetXa[29,14] + lfrSCa[29,14] >= -1
 ldetSC: lDetX[29,15] - lDetY[29,15] = 0
 ldetSCa: ldetXa[29,15] - ldetYa[29,15] = 0
 lfrSC: lDX[29,15] - ldx[29,15] - lFrSC[29,15] >= 0
 lfrSC: lDetX[29,15] - lFrSC[29,15] >= 0
 lfrSC: - lDX[29,15] + ldx[29,15] - lDetX[29,15] + lFrSC[29,15] >= -1
 lfrSCa: lDX[29,15] - ldx[29,15] - lfrSCa[29,15] >= 0
 lfrSCa: ldetXa[29,15] - lfrSCa[29,15] >= 0
 lfrSCa: - lDX[29,15] + ldx[29,15] - ldetXa[29,15] + lfrSCa[29,15] >= -1
 ldetSC: lDetX[30,0] - lDetY[30,0] = 0
 ldetSCa: ldetXa[30,0] - ldetYa[30,0] = 0
 lfrSC: lDX[30,0] - ldx[30,0] - lFrSC[30,0] >= 0
 lfrSC: lDetX[30,0] - lFrSC[30,0] >= 0
 lfrSC: - lDX[30,0] + ldx[30,0] - lDetX[30,0] + lFrSC[30,0] >= -1
 lfrSCa: lDX[30,0] - ldx[30,0] - lfrSCa[30,0] >= 0
 lfrSCa: ldetXa[30,0] - lfrSCa[30,0] >= 0
 lfrSCa: - lDX[30,0] + ldx[30,0] - ldetXa[30,0] + lfrSCa[30,0] >= -1
 ldetSC: lDetX[30,1] - lDetY[30,1] = 0
 ldetSCa: ldetXa[30,1] - ldetYa[30,1] = 0
 lfrSC: lDX[30,1] - ldx[30,1] - lFrSC[30,1] >= 0
 lfrSC: lDetX[30,1] - lFrSC[30,1] >= 0
 lfrSC: - lDX[30,1] + ldx[30,1] - lDetX[30,1] + lFrSC[30,1] >= -1
 lfrSCa: lDX[30,1] - ldx[30,1] - lfrSCa[30,1] >= 0
 lfrSCa: ldetXa[30,1] - lfrSCa[30,1] >= 0
 lfrSCa: - lDX[30,1] + ldx[30,1] - ldetXa[30,1] + lfrSCa[30,1] >= -1
 ldetSC: lDetX[30,2] - lDetY[30,2] = 0
 ldetSCa: ldetXa[30,2] - ldetYa[30,2] = 0
 lfrSC: lDX[30,2] - ldx[30,2] - lFrSC[30,2] >= 0
 lfrSC: lDetX[30,2] - lFrSC[30,2] >= 0
 lfrSC: - lDX[30,2] + ldx[30,2] - lDetX[30,2] + lFrSC[30,2] >= -1
 lfrSCa: lDX[30,2] - ldx[30,2] - lfrSCa[30,2] >= 0
 lfrSCa: ldetXa[30,2] - lfrSCa[30,2] >= 0
 lfrSCa: - lDX[30,2] + ldx[30,2] - ldetXa[30,2] + lfrSCa[30,2] >= -1
 ldetSC: lDetX[30,3] - lDetY[30,3] = 0
 ldetSCa: ldetXa[30,3] - ldetYa[30,3] = 0
 lfrSC: lDX[30,3] - ldx[30,3] - lFrSC[30,3] >= 0
 lfrSC: lDetX[30,3] - lFrSC[30,3] >= 0
 lfrSC: - lDX[30,3] + ldx[30,3] - lDetX[30,3] + lFrSC[30,3] >= -1
 lfrSCa: lDX[30,3] - ldx[30,3] - lfrSCa[30,3] >= 0
 lfrSCa: ldetXa[30,3] - lfrSCa[30,3] >= 0
 lfrSCa: - lDX[30,3] + ldx[30,3] - ldetXa[30,3] + lfrSCa[30,3] >= -1
 ldetSC: lDetX[30,4] - lDetY[30,4] = 0
 ldetSCa: ldetXa[30,4] - ldetYa[30,4] = 0
 lfrSC: lDX[30,4] - ldx[30,4] - lFrSC[30,4] >= 0
 lfrSC: lDetX[30,4] - lFrSC[30,4] >= 0
 lfrSC: - lDX[30,4] + ldx[30,4] - lDetX[30,4] + lFrSC[30,4] >= -1
 lfrSCa: lDX[30,4] - ldx[30,4] - lfrSCa[30,4] >= 0
 lfrSCa: ldetXa[30,4] - lfrSCa[30,4] >= 0
 lfrSCa: - lDX[30,4] + ldx[30,4] - ldetXa[30,4] + lfrSCa[30,4] >= -1
 ldetSC: lDetX[30,5] - lDetY[30,5] = 0
 ldetSCa: ldetXa[30,5] - ldetYa[30,5] = 0
 lfrSC: lDX[30,5] - ldx[30,5] - lFrSC[30,5] >= 0
 lfrSC: lDetX[30,5] - lFrSC[30,5] >= 0
 lfrSC: - lDX[30,5] + ldx[30,5] - lDetX[30,5] + lFrSC[30,5] >= -1
 lfrSCa: lDX[30,5] - ldx[30,5] - lfrSCa[30,5] >= 0
 lfrSCa: ldetXa[30,5] - lfrSCa[30,5] >= 0
 lfrSCa: - lDX[30,5] + ldx[30,5] - ldetXa[30,5] + lfrSCa[30,5] >= -1
 ldetSC: lDetX[30,6] - lDetY[30,6] = 0
 ldetSCa: ldetXa[30,6] - ldetYa[30,6] = 0
 lfrSC: lDX[30,6] - ldx[30,6] - lFrSC[30,6] >= 0
 lfrSC: lDetX[30,6] - lFrSC[30,6] >= 0
 lfrSC: - lDX[30,6] + ldx[30,6] - lDetX[30,6] + lFrSC[30,6] >= -1
 lfrSCa: lDX[30,6] - ldx[30,6] - lfrSCa[30,6] >= 0
 lfrSCa: ldetXa[30,6] - lfrSCa[30,6] >= 0
 lfrSCa: - lDX[30,6] + ldx[30,6] - ldetXa[30,6] + lfrSCa[30,6] >= -1
 ldetSC: lDetX[30,7] - lDetY[30,7] = 0
 ldetSCa: ldetXa[30,7] - ldetYa[30,7] = 0
 lfrSC: lDX[30,7] - ldx[30,7] - lFrSC[30,7] >= 0
 lfrSC: lDetX[30,7] - lFrSC[30,7] >= 0
 lfrSC: - lDX[30,7] + ldx[30,7] - lDetX[30,7] + lFrSC[30,7] >= -1
 lfrSCa: lDX[30,7] - ldx[30,7] - lfrSCa[30,7] >= 0
 lfrSCa: ldetXa[30,7] - lfrSCa[30,7] >= 0
 lfrSCa: - lDX[30,7] + ldx[30,7] - ldetXa[30,7] + lfrSCa[30,7] >= -1
 ldetSC: lDetX[30,8] - lDetY[30,8] = 0
 ldetSCa: ldetXa[30,8] - ldetYa[30,8] = 0
 lfrSC: lDX[30,8] - ldx[30,8] - lFrSC[30,8] >= 0
 lfrSC: lDetX[30,8] - lFrSC[30,8] >= 0
 lfrSC: - lDX[30,8] + ldx[30,8] - lDetX[30,8] + lFrSC[30,8] >= -1
 lfrSCa: lDX[30,8] - ldx[30,8] - lfrSCa[30,8] >= 0
 lfrSCa: ldetXa[30,8] - lfrSCa[30,8] >= 0
 lfrSCa: - lDX[30,8] + ldx[30,8] - ldetXa[30,8] + lfrSCa[30,8] >= -1
 ldetSC: lDetX[30,9] - lDetY[30,9] = 0
 ldetSCa: ldetXa[30,9] - ldetYa[30,9] = 0
 lfrSC: lDX[30,9] - ldx[30,9] - lFrSC[30,9] >= 0
 lfrSC: lDetX[30,9] - lFrSC[30,9] >= 0
 lfrSC: - lDX[30,9] + ldx[30,9] - lDetX[30,9] + lFrSC[30,9] >= -1
 lfrSCa: lDX[30,9] - ldx[30,9] - lfrSCa[30,9] >= 0
 lfrSCa: ldetXa[30,9] - lfrSCa[30,9] >= 0
 lfrSCa: - lDX[30,9] + ldx[30,9] - ldetXa[30,9] + lfrSCa[30,9] >= -1
 ldetSC: lDetX[30,10] - lDetY[30,10] = 0
 ldetSCa: ldetXa[30,10] - ldetYa[30,10] = 0
 lfrSC: lDX[30,10] - ldx[30,10] - lFrSC[30,10] >= 0
 lfrSC: lDetX[30,10] - lFrSC[30,10] >= 0
 lfrSC: - lDX[30,10] + ldx[30,10] - lDetX[30,10] + lFrSC[30,10] >= -1
 lfrSCa: lDX[30,10] - ldx[30,10] - lfrSCa[30,10] >= 0
 lfrSCa: ldetXa[30,10] - lfrSCa[30,10] >= 0
 lfrSCa: - lDX[30,10] + ldx[30,10] - ldetXa[30,10] + lfrSCa[30,10] >= -1
 ldetSC: lDetX[30,11] - lDetY[30,11] = 0
 ldetSCa: ldetXa[30,11] - ldetYa[30,11] = 0
 lfrSC: lDX[30,11] - ldx[30,11] - lFrSC[30,11] >= 0
 lfrSC: lDetX[30,11] - lFrSC[30,11] >= 0
 lfrSC: - lDX[30,11] + ldx[30,11] - lDetX[30,11] + lFrSC[30,11] >= -1
 lfrSCa: lDX[30,11] - ldx[30,11] - lfrSCa[30,11] >= 0
 lfrSCa: ldetXa[30,11] - lfrSCa[30,11] >= 0
 lfrSCa: - lDX[30,11] + ldx[30,11] - ldetXa[30,11] + lfrSCa[30,11] >= -1
 ldetSC: lDetX[30,12] - lDetY[30,12] = 0
 ldetSCa: ldetXa[30,12] - ldetYa[30,12] = 0
 lfrSC: lDX[30,12] - ldx[30,12] - lFrSC[30,12] >= 0
 lfrSC: lDetX[30,12] - lFrSC[30,12] >= 0
 lfrSC: - lDX[30,12] + ldx[30,12] - lDetX[30,12] + lFrSC[30,12] >= -1
 lfrSCa: lDX[30,12] - ldx[30,12] - lfrSCa[30,12] >= 0
 lfrSCa: ldetXa[30,12] - lfrSCa[30,12] >= 0
 lfrSCa: - lDX[30,12] + ldx[30,12] - ldetXa[30,12] + lfrSCa[30,12] >= -1
 ldetSC: lDetX[30,13] - lDetY[30,13] = 0
 ldetSCa: ldetXa[30,13] - ldetYa[30,13] = 0
 lfrSC: lDX[30,13] - ldx[30,13] - lFrSC[30,13] >= 0
 lfrSC: lDetX[30,13] - lFrSC[30,13] >= 0
 lfrSC: - lDX[30,13] + ldx[30,13] - lDetX[30,13] + lFrSC[30,13] >= -1
 lfrSCa: lDX[30,13] - ldx[30,13] - lfrSCa[30,13] >= 0
 lfrSCa: ldetXa[30,13] - lfrSCa[30,13] >= 0
 lfrSCa: - lDX[30,13] + ldx[30,13] - ldetXa[30,13] + lfrSCa[30,13] >= -1
 ldetSC: lDetX[30,14] - lDetY[30,14] = 0
 ldetSCa: ldetXa[30,14] - ldetYa[30,14] = 0
 lfrSC: lDX[30,14] - ldx[30,14] - lFrSC[30,14] >= 0
 lfrSC: lDetX[30,14] - lFrSC[30,14] >= 0
 lfrSC: - lDX[30,14] + ldx[30,14] - lDetX[30,14] + lFrSC[30,14] >= -1
 lfrSCa: lDX[30,14] - ldx[30,14] - lfrSCa[30,14] >= 0
 lfrSCa: ldetXa[30,14] - lfrSCa[30,14] >= 0
 lfrSCa: - lDX[30,14] + ldx[30,14] - ldetXa[30,14] + lfrSCa[30,14] >= -1
 ldetSC: lDetX[30,15] - lDetY[30,15] = 0
 ldetSCa: ldetXa[30,15] - ldetYa[30,15] = 0
 lfrSC: lDX[30,15] - ldx[30,15] - lFrSC[30,15] >= 0
 lfrSC: lDetX[30,15] - lFrSC[30,15] >= 0
 lfrSC: - lDX[30,15] + ldx[30,15] - lDetX[30,15] + lFrSC[30,15] >= -1
 lfrSCa: lDX[30,15] - ldx[30,15] - lfrSCa[30,15] >= 0
 lfrSCa: ldetXa[30,15] - lfrSCa[30,15] >= 0
 lfrSCa: - lDX[30,15] + ldx[30,15] - ldetXa[30,15] + lfrSCa[30,15] >= -1
 ldetSC: lDetX[31,0] - lDetY[31,0] = 0
 ldetSCa: ldetXa[31,0] - ldetYa[31,0] = 0
 lfrSC: lDX[31,0] - ldx[31,0] - lFrSC[31,0] >= 0
 lfrSC: lDetX[31,0] - lFrSC[31,0] >= 0
 lfrSC: - lDX[31,0] + ldx[31,0] - lDetX[31,0] + lFrSC[31,0] >= -1
 lfrSCa: lDX[31,0] - ldx[31,0] - lfrSCa[31,0] >= 0
 lfrSCa: ldetXa[31,0] - lfrSCa[31,0] >= 0
 lfrSCa: - lDX[31,0] + ldx[31,0] - ldetXa[31,0] + lfrSCa[31,0] >= -1
 ldetSC: lDetX[31,1] - lDetY[31,1] = 0
 ldetSCa: ldetXa[31,1] - ldetYa[31,1] = 0
 lfrSC: lDX[31,1] - ldx[31,1] - lFrSC[31,1] >= 0
 lfrSC: lDetX[31,1] - lFrSC[31,1] >= 0
 lfrSC: - lDX[31,1] + ldx[31,1] - lDetX[31,1] + lFrSC[31,1] >= -1
 lfrSCa: lDX[31,1] - ldx[31,1] - lfrSCa[31,1] >= 0
 lfrSCa: ldetXa[31,1] - lfrSCa[31,1] >= 0
 lfrSCa: - lDX[31,1] + ldx[31,1] - ldetXa[31,1] + lfrSCa[31,1] >= -1
 ldetSC: lDetX[31,2] - lDetY[31,2] = 0
 ldetSCa: ldetXa[31,2] - ldetYa[31,2] = 0
 lfrSC: lDX[31,2] - ldx[31,2] - lFrSC[31,2] >= 0
 lfrSC: lDetX[31,2] - lFrSC[31,2] >= 0
 lfrSC: - lDX[31,2] + ldx[31,2] - lDetX[31,2] + lFrSC[31,2] >= -1
 lfrSCa: lDX[31,2] - ldx[31,2] - lfrSCa[31,2] >= 0
 lfrSCa: ldetXa[31,2] - lfrSCa[31,2] >= 0
 lfrSCa: - lDX[31,2] + ldx[31,2] - ldetXa[31,2] + lfrSCa[31,2] >= -1
 ldetSC: lDetX[31,3] - lDetY[31,3] = 0
 ldetSCa: ldetXa[31,3] - ldetYa[31,3] = 0
 lfrSC: lDX[31,3] - ldx[31,3] - lFrSC[31,3] >= 0
 lfrSC: lDetX[31,3] - lFrSC[31,3] >= 0
 lfrSC: - lDX[31,3] + ldx[31,3] - lDetX[31,3] + lFrSC[31,3] >= -1
 lfrSCa: lDX[31,3] - ldx[31,3] - lfrSCa[31,3] >= 0
 lfrSCa: ldetXa[31,3] - lfrSCa[31,3] >= 0
 lfrSCa: - lDX[31,3] + ldx[31,3] - ldetXa[31,3] + lfrSCa[31,3] >= -1
 ldetSC: lDetX[31,4] - lDetY[31,4] = 0
 ldetSCa: ldetXa[31,4] - ldetYa[31,4] = 0
 lfrSC: lDX[31,4] - ldx[31,4] - lFrSC[31,4] >= 0
 lfrSC: lDetX[31,4] - lFrSC[31,4] >= 0
 lfrSC: - lDX[31,4] + ldx[31,4] - lDetX[31,4] + lFrSC[31,4] >= -1
 lfrSCa: lDX[31,4] - ldx[31,4] - lfrSCa[31,4] >= 0
 lfrSCa: ldetXa[31,4] - lfrSCa[31,4] >= 0
 lfrSCa: - lDX[31,4] + ldx[31,4] - ldetXa[31,4] + lfrSCa[31,4] >= -1
 ldetSC: lDetX[31,5] - lDetY[31,5] = 0
 ldetSCa: ldetXa[31,5] - ldetYa[31,5] = 0
 lfrSC: lDX[31,5] - ldx[31,5] - lFrSC[31,5] >= 0
 lfrSC: lDetX[31,5] - lFrSC[31,5] >= 0
 lfrSC: - lDX[31,5] + ldx[31,5] - lDetX[31,5] + lFrSC[31,5] >= -1
 lfrSCa: lDX[31,5] - ldx[31,5] - lfrSCa[31,5] >= 0
 lfrSCa: ldetXa[31,5] - lfrSCa[31,5] >= 0
 lfrSCa: - lDX[31,5] + ldx[31,5] - ldetXa[31,5] + lfrSCa[31,5] >= -1
 ldetSC: lDetX[31,6] - lDetY[31,6] = 0
 ldetSCa: ldetXa[31,6] - ldetYa[31,6] = 0
 lfrSC: lDX[31,6] - ldx[31,6] - lFrSC[31,6] >= 0
 lfrSC: lDetX[31,6] - lFrSC[31,6] >= 0
 lfrSC: - lDX[31,6] + ldx[31,6] - lDetX[31,6] + lFrSC[31,6] >= -1
 lfrSCa: lDX[31,6] - ldx[31,6] - lfrSCa[31,6] >= 0
 lfrSCa: ldetXa[31,6] - lfrSCa[31,6] >= 0
 lfrSCa: - lDX[31,6] + ldx[31,6] - ldetXa[31,6] + lfrSCa[31,6] >= -1
 ldetSC: lDetX[31,7] - lDetY[31,7] = 0
 ldetSCa: ldetXa[31,7] - ldetYa[31,7] = 0
 lfrSC: lDX[31,7] - ldx[31,7] - lFrSC[31,7] >= 0
 lfrSC: lDetX[31,7] - lFrSC[31,7] >= 0
 lfrSC: - lDX[31,7] + ldx[31,7] - lDetX[31,7] + lFrSC[31,7] >= -1
 lfrSCa: lDX[31,7] - ldx[31,7] - lfrSCa[31,7] >= 0
 lfrSCa: ldetXa[31,7] - lfrSCa[31,7] >= 0
 lfrSCa: - lDX[31,7] + ldx[31,7] - ldetXa[31,7] + lfrSCa[31,7] >= -1
 ldetSC: lDetX[31,8] - lDetY[31,8] = 0
 ldetSCa: ldetXa[31,8] - ldetYa[31,8] = 0
 lfrSC: lDX[31,8] - ldx[31,8] - lFrSC[31,8] >= 0
 lfrSC: lDetX[31,8] - lFrSC[31,8] >= 0
 lfrSC: - lDX[31,8] + ldx[31,8] - lDetX[31,8] + lFrSC[31,8] >= -1
 lfrSCa: lDX[31,8] - ldx[31,8] - lfrSCa[31,8] >= 0
 lfrSCa: ldetXa[31,8] - lfrSCa[31,8] >= 0
 lfrSCa: - lDX[31,8] + ldx[31,8] - ldetXa[31,8] + lfrSCa[31,8] >= -1
 ldetSC: lDetX[31,9] - lDetY[31,9] = 0
 ldetSCa: ldetXa[31,9] - ldetYa[31,9] = 0
 lfrSC: lDX[31,9] - ldx[31,9] - lFrSC[31,9] >= 0
 lfrSC: lDetX[31,9] - lFrSC[31,9] >= 0
 lfrSC: - lDX[31,9] + ldx[31,9] - lDetX[31,9] + lFrSC[31,9] >= -1
 lfrSCa: lDX[31,9] - ldx[31,9] - lfrSCa[31,9] >= 0
 lfrSCa: ldetXa[31,9] - lfrSCa[31,9] >= 0
 lfrSCa: - lDX[31,9] + ldx[31,9] - ldetXa[31,9] + lfrSCa[31,9] >= -1
 ldetSC: lDetX[31,10] - lDetY[31,10] = 0
 ldetSCa: ldetXa[31,10] - ldetYa[31,10] = 0
 lfrSC: lDX[31,10] - ldx[31,10] - lFrSC[31,10] >= 0
 lfrSC: lDetX[31,10] - lFrSC[31,10] >= 0
 lfrSC: - lDX[31,10] + ldx[31,10] - lDetX[31,10] + lFrSC[31,10] >= -1
 lfrSCa: lDX[31,10] - ldx[31,10] - lfrSCa[31,10] >= 0
 lfrSCa: ldetXa[31,10] - lfrSCa[31,10] >= 0
 lfrSCa: - lDX[31,10] + ldx[31,10] - ldetXa[31,10] + lfrSCa[31,10] >= -1
 ldetSC: lDetX[31,11] - lDetY[31,11] = 0
 ldetSCa: ldetXa[31,11] - ldetYa[31,11] = 0
 lfrSC: lDX[31,11] - ldx[31,11] - lFrSC[31,11] >= 0
 lfrSC: lDetX[31,11] - lFrSC[31,11] >= 0
 lfrSC: - lDX[31,11] + ldx[31,11] - lDetX[31,11] + lFrSC[31,11] >= -1
 lfrSCa: lDX[31,11] - ldx[31,11] - lfrSCa[31,11] >= 0
 lfrSCa: ldetXa[31,11] - lfrSCa[31,11] >= 0
 lfrSCa: - lDX[31,11] + ldx[31,11] - ldetXa[31,11] + lfrSCa[31,11] >= -1
 ldetSC: lDetX[31,12] - lDetY[31,12] = 0
 ldetSCa: ldetXa[31,12] - ldetYa[31,12] = 0
 lfrSC: lDX[31,12] - ldx[31,12] - lFrSC[31,12] >= 0
 lfrSC: lDetX[31,12] - lFrSC[31,12] >= 0
 lfrSC: - lDX[31,12] + ldx[31,12] - lDetX[31,12] + lFrSC[31,12] >= -1
 lfrSCa: lDX[31,12] - ldx[31,12] - lfrSCa[31,12] >= 0
 lfrSCa: ldetXa[31,12] - lfrSCa[31,12] >= 0
 lfrSCa: - lDX[31,12] + ldx[31,12] - ldetXa[31,12] + lfrSCa[31,12] >= -1
 ldetSC: lDetX[31,13] - lDetY[31,13] = 0
 ldetSCa: ldetXa[31,13] - ldetYa[31,13] = 0
 lfrSC: lDX[31,13] - ldx[31,13] - lFrSC[31,13] >= 0
 lfrSC: lDetX[31,13] - lFrSC[31,13] >= 0
 lfrSC: - lDX[31,13] + ldx[31,13] - lDetX[31,13] + lFrSC[31,13] >= -1
 lfrSCa: lDX[31,13] - ldx[31,13] - lfrSCa[31,13] >= 0
 lfrSCa: ldetXa[31,13] - lfrSCa[31,13] >= 0
 lfrSCa: - lDX[31,13] + ldx[31,13] - ldetXa[31,13] + lfrSCa[31,13] >= -1
 ldetSC: lDetX[31,14] - lDetY[31,14] = 0
 ldetSCa: ldetXa[31,14] - ldetYa[31,14] = 0
 lfrSC: lDX[31,14] - ldx[31,14] - lFrSC[31,14] >= 0
 lfrSC: lDetX[31,14] - lFrSC[31,14] >= 0
 lfrSC: - lDX[31,14] + ldx[31,14] - lDetX[31,14] + lFrSC[31,14] >= -1
 lfrSCa: lDX[31,14] - ldx[31,14] - lfrSCa[31,14] >= 0
 lfrSCa: ldetXa[31,14] - lfrSCa[31,14] >= 0
 lfrSCa: - lDX[31,14] + ldx[31,14] - ldetXa[31,14] + lfrSCa[31,14] >= -1
 ldetSC: lDetX[31,15] - lDetY[31,15] = 0
 ldetSCa: ldetXa[31,15] - ldetYa[31,15] = 0
 lfrSC: lDX[31,15] - ldx[31,15] - lFrSC[31,15] >= 0
 lfrSC: lDetX[31,15] - lFrSC[31,15] >= 0
 lfrSC: - lDX[31,15] + ldx[31,15] - lDetX[31,15] + lFrSC[31,15] >= -1
 lfrSCa: lDX[31,15] - ldx[31,15] - lfrSCa[31,15] >= 0
 lfrSCa: ldetXa[31,15] - lfrSCa[31,15] >= 0
 lfrSCa: - lDX[31,15] + ldx[31,15] - ldetXa[31,15] + lfrSCa[31,15] >= -1
 ldetSC: lDetX[32,0] - lDetY[32,0] = 0
 ldetSCa: ldetXa[32,0] - ldetYa[32,0] = 0
 lfrSC: lDX[32,0] - ldx[32,0] - lFrSC[32,0] >= 0
 lfrSC: lDetX[32,0] - lFrSC[32,0] >= 0
 lfrSC: - lDX[32,0] + ldx[32,0] - lDetX[32,0] + lFrSC[32,0] >= -1
 lfrSCa: lDX[32,0] - ldx[32,0] - lfrSCa[32,0] >= 0
 lfrSCa: ldetXa[32,0] - lfrSCa[32,0] >= 0
 lfrSCa: - lDX[32,0] + ldx[32,0] - ldetXa[32,0] + lfrSCa[32,0] >= -1
 ldetSC: lDetX[32,1] - lDetY[32,1] = 0
 ldetSCa: ldetXa[32,1] - ldetYa[32,1] = 0
 lfrSC: lDX[32,1] - ldx[32,1] - lFrSC[32,1] >= 0
 lfrSC: lDetX[32,1] - lFrSC[32,1] >= 0
 lfrSC: - lDX[32,1] + ldx[32,1] - lDetX[32,1] + lFrSC[32,1] >= -1
 lfrSCa: lDX[32,1] - ldx[32,1] - lfrSCa[32,1] >= 0
 lfrSCa: ldetXa[32,1] - lfrSCa[32,1] >= 0
 lfrSCa: - lDX[32,1] + ldx[32,1] - ldetXa[32,1] + lfrSCa[32,1] >= -1
 ldetSC: lDetX[32,2] - lDetY[32,2] = 0
 ldetSCa: ldetXa[32,2] - ldetYa[32,2] = 0
 lfrSC: lDX[32,2] - ldx[32,2] - lFrSC[32,2] >= 0
 lfrSC: lDetX[32,2] - lFrSC[32,2] >= 0
 lfrSC: - lDX[32,2] + ldx[32,2] - lDetX[32,2] + lFrSC[32,2] >= -1
 lfrSCa: lDX[32,2] - ldx[32,2] - lfrSCa[32,2] >= 0
 lfrSCa: ldetXa[32,2] - lfrSCa[32,2] >= 0
 lfrSCa: - lDX[32,2] + ldx[32,2] - ldetXa[32,2] + lfrSCa[32,2] >= -1
 ldetSC: lDetX[32,3] - lDetY[32,3] = 0
 ldetSCa: ldetXa[32,3] - ldetYa[32,3] = 0
 lfrSC: lDX[32,3] - ldx[32,3] - lFrSC[32,3] >= 0
 lfrSC: lDetX[32,3] - lFrSC[32,3] >= 0
 lfrSC: - lDX[32,3] + ldx[32,3] - lDetX[32,3] + lFrSC[32,3] >= -1
 lfrSCa: lDX[32,3] - ldx[32,3] - lfrSCa[32,3] >= 0
 lfrSCa: ldetXa[32,3] - lfrSCa[32,3] >= 0
 lfrSCa: - lDX[32,3] + ldx[32,3] - ldetXa[32,3] + lfrSCa[32,3] >= -1
 ldetSC: lDetX[32,4] - lDetY[32,4] = 0
 ldetSCa: ldetXa[32,4] - ldetYa[32,4] = 0
 lfrSC: lDX[32,4] - ldx[32,4] - lFrSC[32,4] >= 0
 lfrSC: lDetX[32,4] - lFrSC[32,4] >= 0
 lfrSC: - lDX[32,4] + ldx[32,4] - lDetX[32,4] + lFrSC[32,4] >= -1
 lfrSCa: lDX[32,4] - ldx[32,4] - lfrSCa[32,4] >= 0
 lfrSCa: ldetXa[32,4] - lfrSCa[32,4] >= 0
 lfrSCa: - lDX[32,4] + ldx[32,4] - ldetXa[32,4] + lfrSCa[32,4] >= -1
 ldetSC: lDetX[32,5] - lDetY[32,5] = 0
 ldetSCa: ldetXa[32,5] - ldetYa[32,5] = 0
 lfrSC: lDX[32,5] - ldx[32,5] - lFrSC[32,5] >= 0
 lfrSC: lDetX[32,5] - lFrSC[32,5] >= 0
 lfrSC: - lDX[32,5] + ldx[32,5] - lDetX[32,5] + lFrSC[32,5] >= -1
 lfrSCa: lDX[32,5] - ldx[32,5] - lfrSCa[32,5] >= 0
 lfrSCa: ldetXa[32,5] - lfrSCa[32,5] >= 0
 lfrSCa: - lDX[32,5] + ldx[32,5] - ldetXa[32,5] + lfrSCa[32,5] >= -1
 ldetSC: lDetX[32,6] - lDetY[32,6] = 0
 ldetSCa: ldetXa[32,6] - ldetYa[32,6] = 0
 lfrSC: lDX[32,6] - ldx[32,6] - lFrSC[32,6] >= 0
 lfrSC: lDetX[32,6] - lFrSC[32,6] >= 0
 lfrSC: - lDX[32,6] + ldx[32,6] - lDetX[32,6] + lFrSC[32,6] >= -1
 lfrSCa: lDX[32,6] - ldx[32,6] - lfrSCa[32,6] >= 0
 lfrSCa: ldetXa[32,6] - lfrSCa[32,6] >= 0
 lfrSCa: - lDX[32,6] + ldx[32,6] - ldetXa[32,6] + lfrSCa[32,6] >= -1
 ldetSC: lDetX[32,7] - lDetY[32,7] = 0
 ldetSCa: ldetXa[32,7] - ldetYa[32,7] = 0
 lfrSC: lDX[32,7] - ldx[32,7] - lFrSC[32,7] >= 0
 lfrSC: lDetX[32,7] - lFrSC[32,7] >= 0
 lfrSC: - lDX[32,7] + ldx[32,7] - lDetX[32,7] + lFrSC[32,7] >= -1
 lfrSCa: lDX[32,7] - ldx[32,7] - lfrSCa[32,7] >= 0
 lfrSCa: ldetXa[32,7] - lfrSCa[32,7] >= 0
 lfrSCa: - lDX[32,7] + ldx[32,7] - ldetXa[32,7] + lfrSCa[32,7] >= -1
 ldetSC: lDetX[32,8] - lDetY[32,8] = 0
 ldetSCa: ldetXa[32,8] - ldetYa[32,8] = 0
 lfrSC: lDX[32,8] - ldx[32,8] - lFrSC[32,8] >= 0
 lfrSC: lDetX[32,8] - lFrSC[32,8] >= 0
 lfrSC: - lDX[32,8] + ldx[32,8] - lDetX[32,8] + lFrSC[32,8] >= -1
 lfrSCa: lDX[32,8] - ldx[32,8] - lfrSCa[32,8] >= 0
 lfrSCa: ldetXa[32,8] - lfrSCa[32,8] >= 0
 lfrSCa: - lDX[32,8] + ldx[32,8] - ldetXa[32,8] + lfrSCa[32,8] >= -1
 ldetSC: lDetX[32,9] - lDetY[32,9] = 0
 ldetSCa: ldetXa[32,9] - ldetYa[32,9] = 0
 lfrSC: lDX[32,9] - ldx[32,9] - lFrSC[32,9] >= 0
 lfrSC: lDetX[32,9] - lFrSC[32,9] >= 0
 lfrSC: - lDX[32,9] + ldx[32,9] - lDetX[32,9] + lFrSC[32,9] >= -1
 lfrSCa: lDX[32,9] - ldx[32,9] - lfrSCa[32,9] >= 0
 lfrSCa: ldetXa[32,9] - lfrSCa[32,9] >= 0
 lfrSCa: - lDX[32,9] + ldx[32,9] - ldetXa[32,9] + lfrSCa[32,9] >= -1
 ldetSC: lDetX[32,10] - lDetY[32,10] = 0
 ldetSCa: ldetXa[32,10] - ldetYa[32,10] = 0
 lfrSC: lDX[32,10] - ldx[32,10] - lFrSC[32,10] >= 0
 lfrSC: lDetX[32,10] - lFrSC[32,10] >= 0
 lfrSC: - lDX[32,10] + ldx[32,10] - lDetX[32,10] + lFrSC[32,10] >= -1
 lfrSCa: lDX[32,10] - ldx[32,10] - lfrSCa[32,10] >= 0
 lfrSCa: ldetXa[32,10] - lfrSCa[32,10] >= 0
 lfrSCa: - lDX[32,10] + ldx[32,10] - ldetXa[32,10] + lfrSCa[32,10] >= -1
 ldetSC: lDetX[32,11] - lDetY[32,11] = 0
 ldetSCa: ldetXa[32,11] - ldetYa[32,11] = 0
 lfrSC: lDX[32,11] - ldx[32,11] - lFrSC[32,11] >= 0
 lfrSC: lDetX[32,11] - lFrSC[32,11] >= 0
 lfrSC: - lDX[32,11] + ldx[32,11] - lDetX[32,11] + lFrSC[32,11] >= -1
 lfrSCa: lDX[32,11] - ldx[32,11] - lfrSCa[32,11] >= 0
 lfrSCa: ldetXa[32,11] - lfrSCa[32,11] >= 0
 lfrSCa: - lDX[32,11] + ldx[32,11] - ldetXa[32,11] + lfrSCa[32,11] >= -1
 ldetSC: lDetX[32,12] - lDetY[32,12] = 0
 ldetSCa: ldetXa[32,12] - ldetYa[32,12] = 0
 lfrSC: lDX[32,12] - ldx[32,12] - lFrSC[32,12] >= 0
 lfrSC: lDetX[32,12] - lFrSC[32,12] >= 0
 lfrSC: - lDX[32,12] + ldx[32,12] - lDetX[32,12] + lFrSC[32,12] >= -1
 lfrSCa: lDX[32,12] - ldx[32,12] - lfrSCa[32,12] >= 0
 lfrSCa: ldetXa[32,12] - lfrSCa[32,12] >= 0
 lfrSCa: - lDX[32,12] + ldx[32,12] - ldetXa[32,12] + lfrSCa[32,12] >= -1
 ldetSC: lDetX[32,13] - lDetY[32,13] = 0
 ldetSCa: ldetXa[32,13] - ldetYa[32,13] = 0
 lfrSC: lDX[32,13] - ldx[32,13] - lFrSC[32,13] >= 0
 lfrSC: lDetX[32,13] - lFrSC[32,13] >= 0
 lfrSC: - lDX[32,13] + ldx[32,13] - lDetX[32,13] + lFrSC[32,13] >= -1
 lfrSCa: lDX[32,13] - ldx[32,13] - lfrSCa[32,13] >= 0
 lfrSCa: ldetXa[32,13] - lfrSCa[32,13] >= 0
 lfrSCa: - lDX[32,13] + ldx[32,13] - ldetXa[32,13] + lfrSCa[32,13] >= -1
 ldetSC: lDetX[32,14] - lDetY[32,14] = 0
 ldetSCa: ldetXa[32,14] - ldetYa[32,14] = 0
 lfrSC: lDX[32,14] - ldx[32,14] - lFrSC[32,14] >= 0
 lfrSC: lDetX[32,14] - lFrSC[32,14] >= 0
 lfrSC: - lDX[32,14] + ldx[32,14] - lDetX[32,14] + lFrSC[32,14] >= -1
 lfrSCa: lDX[32,14] - ldx[32,14] - lfrSCa[32,14] >= 0
 lfrSCa: ldetXa[32,14] - lfrSCa[32,14] >= 0
 lfrSCa: - lDX[32,14] + ldx[32,14] - ldetXa[32,14] + lfrSCa[32,14] >= -1
 ldetSC: lDetX[32,15] - lDetY[32,15] = 0
 ldetSCa: ldetXa[32,15] - ldetYa[32,15] = 0
 lfrSC: lDX[32,15] - ldx[32,15] - lFrSC[32,15] >= 0
 lfrSC: lDetX[32,15] - lFrSC[32,15] >= 0
 lfrSC: - lDX[32,15] + ldx[32,15] - lDetX[32,15] + lFrSC[32,15] >= -1
 lfrSCa: lDX[32,15] - ldx[32,15] - lfrSCa[32,15] >= 0
 lfrSCa: ldetXa[32,15] - lfrSCa[32,15] >= 0
 lfrSCa: - lDX[32,15] + ldx[32,15] - ldetXa[32,15] + lfrSCa[32,15] >= -1
 lgkiU: - lGstk[28,0] + lDetY[28,0] <= 0
 lgkiU: lDetY[28,0] - lDetZ[28,0] <= 0
 lgkiU: - lGstk[28,0] + lDetY[28,0] - lDetZ[28,0] >= -1
 lgkiUa: - lgstka[28,0] + ldetYa[28,0] <= 0
 lgkiUa: ldetYa[28,0] - ldetZa[28,0] <= 0
 lgkiUa: - lgstka[28,0] + ldetYa[28,0] - ldetZa[28,0] >= -1
 lgkiU: - lGstk[28,1] + lDetY[28,1] <= 0
 lgkiU: lDetY[28,1] - lDetZ[28,1] <= 0
 lgkiU: - lGstk[28,1] + lDetY[28,1] - lDetZ[28,1] >= -1
 lgkiUa: - lgstka[28,1] + ldetYa[28,1] <= 0
 lgkiUa: ldetYa[28,1] - ldetZa[28,1] <= 0
 lgkiUa: - lgstka[28,1] + ldetYa[28,1] - ldetZa[28,1] >= -1
 lgkiU: - lGstk[28,2] + lDetY[28,2] <= 0
 lgkiU: lDetY[28,2] - lDetZ[28,2] <= 0
 lgkiU: - lGstk[28,2] + lDetY[28,2] - lDetZ[28,2] >= -1
 lgkiUa: - lgstka[28,2] + ldetYa[28,2] <= 0
 lgkiUa: ldetYa[28,2] - ldetZa[28,2] <= 0
 lgkiUa: - lgstka[28,2] + ldetYa[28,2] - ldetZa[28,2] >= -1
 lgkiU: - lGstk[28,3] + lDetY[28,3] <= 0
 lgkiU: lDetY[28,3] - lDetZ[28,3] <= 0
 lgkiU: - lGstk[28,3] + lDetY[28,3] - lDetZ[28,3] >= -1
 lgkiUa: - lgstka[28,3] + ldetYa[28,3] <= 0
 lgkiUa: ldetYa[28,3] - ldetZa[28,3] <= 0
 lgkiUa: - lgstka[28,3] + ldetYa[28,3] - ldetZa[28,3] >= -1
 lgkiU: - lGstk[28,4] + lDetY[28,4] <= 0
 lgkiU: lDetY[28,4] - lDetZ[28,4] <= 0
 lgkiU: - lGstk[28,4] + lDetY[28,4] - lDetZ[28,4] >= -1
 lgkiUa: - lgstka[28,4] + ldetYa[28,4] <= 0
 lgkiUa: ldetYa[28,4] - ldetZa[28,4] <= 0
 lgkiUa: - lgstka[28,4] + ldetYa[28,4] - ldetZa[28,4] >= -1
 lgkiU: - lGstk[28,5] + lDetY[28,5] <= 0
 lgkiU: lDetY[28,5] - lDetZ[28,5] <= 0
 lgkiU: - lGstk[28,5] + lDetY[28,5] - lDetZ[28,5] >= -1
 lgkiUa: - lgstka[28,5] + ldetYa[28,5] <= 0
 lgkiUa: ldetYa[28,5] - ldetZa[28,5] <= 0
 lgkiUa: - lgstka[28,5] + ldetYa[28,5] - ldetZa[28,5] >= -1
 lgkiU: - lGstk[28,6] + lDetY[28,6] <= 0
 lgkiU: lDetY[28,6] - lDetZ[28,6] <= 0
 lgkiU: - lGstk[28,6] + lDetY[28,6] - lDetZ[28,6] >= -1
 lgkiUa: - lgstka[28,6] + ldetYa[28,6] <= 0
 lgkiUa: ldetYa[28,6] - ldetZa[28,6] <= 0
 lgkiUa: - lgstka[28,6] + ldetYa[28,6] - ldetZa[28,6] >= -1
 lgkiU: - lGstk[28,7] + lDetY[28,7] <= 0
 lgkiU: lDetY[28,7] - lDetZ[28,7] <= 0
 lgkiU: - lGstk[28,7] + lDetY[28,7] - lDetZ[28,7] >= -1
 lgkiUa: - lgstka[28,7] + ldetYa[28,7] <= 0
 lgkiUa: ldetYa[28,7] - ldetZa[28,7] <= 0
 lgkiUa: - lgstka[28,7] + ldetYa[28,7] - ldetZa[28,7] >= -1
 lgkiL: lDetY[28,8] - lDetZ[28,8] = 0
 lgkiLa: ldetYa[28,8] - ldetZa[28,8] = 0
 lgkiL: lDetY[28,9] - lDetZ[28,9] = 0
 lgkiLa: ldetYa[28,9] - ldetZa[28,9] = 0
 lgkiL: lDetY[28,10] - lDetZ[28,10] = 0
 lgkiLa: ldetYa[28,10] - ldetZa[28,10] = 0
 lgkiL: lDetY[28,11] - lDetZ[28,11] = 0
 lgkiLa: ldetYa[28,11] - ldetZa[28,11] = 0
 lgkiL: lDetY[28,12] - lDetZ[28,12] = 0
 lgkiLa: ldetYa[28,12] - ldetZa[28,12] = 0
 lgkiL: lDetY[28,13] - lDetZ[28,13] = 0
 lgkiLa: ldetYa[28,13] - ldetZa[28,13] = 0
 lgkiL: lDetY[28,14] - lDetZ[28,14] = 0
 lgkiLa: ldetYa[28,14] - ldetZa[28,14] = 0
 lgkiL: lDetY[28,15] - lDetZ[28,15] = 0
 lgkiLa: ldetYa[28,15] - ldetZa[28,15] = 0
 lgkiU: - lGstk[29,0] + lDetY[29,0] <= 0
 lgkiU: lDetY[29,0] - lDetZ[29,0] <= 0
 lgkiU: - lGstk[29,0] + lDetY[29,0] - lDetZ[29,0] >= -1
 lgkiUa: - lgstka[29,0] + ldetYa[29,0] <= 0
 lgkiUa: ldetYa[29,0] - ldetZa[29,0] <= 0
 lgkiUa: - lgstka[29,0] + ldetYa[29,0] - ldetZa[29,0] >= -1
 lgkiU: - lGstk[29,1] + lDetY[29,1] <= 0
 lgkiU: lDetY[29,1] - lDetZ[29,1] <= 0
 lgkiU: - lGstk[29,1] + lDetY[29,1] - lDetZ[29,1] >= -1
 lgkiUa: - lgstka[29,1] + ldetYa[29,1] <= 0
 lgkiUa: ldetYa[29,1] - ldetZa[29,1] <= 0
 lgkiUa: - lgstka[29,1] + ldetYa[29,1] - ldetZa[29,1] >= -1
 lgkiU: - lGstk[29,2] + lDetY[29,2] <= 0
 lgkiU: lDetY[29,2] - lDetZ[29,2] <= 0
 lgkiU: - lGstk[29,2] + lDetY[29,2] - lDetZ[29,2] >= -1
 lgkiUa: - lgstka[29,2] + ldetYa[29,2] <= 0
 lgkiUa: ldetYa[29,2] - ldetZa[29,2] <= 0
 lgkiUa: - lgstka[29,2] + ldetYa[29,2] - ldetZa[29,2] >= -1
 lgkiU: - lGstk[29,3] + lDetY[29,3] <= 0
 lgkiU: lDetY[29,3] - lDetZ[29,3] <= 0
 lgkiU: - lGstk[29,3] + lDetY[29,3] - lDetZ[29,3] >= -1
 lgkiUa: - lgstka[29,3] + ldetYa[29,3] <= 0
 lgkiUa: ldetYa[29,3] - ldetZa[29,3] <= 0
 lgkiUa: - lgstka[29,3] + ldetYa[29,3] - ldetZa[29,3] >= -1
 lgkiU: - lGstk[29,4] + lDetY[29,4] <= 0
 lgkiU: lDetY[29,4] - lDetZ[29,4] <= 0
 lgkiU: - lGstk[29,4] + lDetY[29,4] - lDetZ[29,4] >= -1
 lgkiUa: - lgstka[29,4] + ldetYa[29,4] <= 0
 lgkiUa: ldetYa[29,4] - ldetZa[29,4] <= 0
 lgkiUa: - lgstka[29,4] + ldetYa[29,4] - ldetZa[29,4] >= -1
 lgkiU: - lGstk[29,5] + lDetY[29,5] <= 0
 lgkiU: lDetY[29,5] - lDetZ[29,5] <= 0
 lgkiU: - lGstk[29,5] + lDetY[29,5] - lDetZ[29,5] >= -1
 lgkiUa: - lgstka[29,5] + ldetYa[29,5] <= 0
 lgkiUa: ldetYa[29,5] - ldetZa[29,5] <= 0
 lgkiUa: - lgstka[29,5] + ldetYa[29,5] - ldetZa[29,5] >= -1
 lgkiU: - lGstk[29,6] + lDetY[29,6] <= 0
 lgkiU: lDetY[29,6] - lDetZ[29,6] <= 0
 lgkiU: - lGstk[29,6] + lDetY[29,6] - lDetZ[29,6] >= -1
 lgkiUa: - lgstka[29,6] + ldetYa[29,6] <= 0
 lgkiUa: ldetYa[29,6] - ldetZa[29,6] <= 0
 lgkiUa: - lgstka[29,6] + ldetYa[29,6] - ldetZa[29,6] >= -1
 lgkiU: - lGstk[29,7] + lDetY[29,7] <= 0
 lgkiU: lDetY[29,7] - lDetZ[29,7] <= 0
 lgkiU: - lGstk[29,7] + lDetY[29,7] - lDetZ[29,7] >= -1
 lgkiUa: - lgstka[29,7] + ldetYa[29,7] <= 0
 lgkiUa: ldetYa[29,7] - ldetZa[29,7] <= 0
 lgkiUa: - lgstka[29,7] + ldetYa[29,7] - ldetZa[29,7] >= -1
 lgkiL: lDetY[29,8] - lDetZ[29,8] = 0
 lgkiLa: ldetYa[29,8] - ldetZa[29,8] = 0
 lgkiL: lDetY[29,9] - lDetZ[29,9] = 0
 lgkiLa: ldetYa[29,9] - ldetZa[29,9] = 0
 lgkiL: lDetY[29,10] - lDetZ[29,10] = 0
 lgkiLa: ldetYa[29,10] - ldetZa[29,10] = 0
 lgkiL: lDetY[29,11] - lDetZ[29,11] = 0
 lgkiLa: ldetYa[29,11] - ldetZa[29,11] = 0
 lgkiL: lDetY[29,12] - lDetZ[29,12] = 0
 lgkiLa: ldetYa[29,12] - ldetZa[29,12] = 0
 lgkiL: lDetY[29,13] - lDetZ[29,13] = 0
 lgkiLa: ldetYa[29,13] - ldetZa[29,13] = 0
 lgkiL: lDetY[29,14] - lDetZ[29,14] = 0
 lgkiLa: ldetYa[29,14] - ldetZa[29,14] = 0
 lgkiL: lDetY[29,15] - lDetZ[29,15] = 0
 lgkiLa: ldetYa[29,15] - ldetZa[29,15] = 0
 lgkiU: - lGstk[30,0] + lDetY[30,0] <= 0
 lgkiU: lDetY[30,0] - lDetZ[30,0] <= 0
 lgkiU: - lGstk[30,0] + lDetY[30,0] - lDetZ[30,0] >= -1
 lgkiUa: - lgstka[30,0] + ldetYa[30,0] <= 0
 lgkiUa: ldetYa[30,0] - ldetZa[30,0] <= 0
 lgkiUa: - lgstka[30,0] + ldetYa[30,0] - ldetZa[30,0] >= -1
 lgkiU: - lGstk[30,1] + lDetY[30,1] <= 0
 lgkiU: lDetY[30,1] - lDetZ[30,1] <= 0
 lgkiU: - lGstk[30,1] + lDetY[30,1] - lDetZ[30,1] >= -1
 lgkiUa: - lgstka[30,1] + ldetYa[30,1] <= 0
 lgkiUa: ldetYa[30,1] - ldetZa[30,1] <= 0
 lgkiUa: - lgstka[30,1] + ldetYa[30,1] - ldetZa[30,1] >= -1
 lgkiU: - lGstk[30,2] + lDetY[30,2] <= 0
 lgkiU: lDetY[30,2] - lDetZ[30,2] <= 0
 lgkiU: - lGstk[30,2] + lDetY[30,2] - lDetZ[30,2] >= -1
 lgkiUa: - lgstka[30,2] + ldetYa[30,2] <= 0
 lgkiUa: ldetYa[30,2] - ldetZa[30,2] <= 0
 lgkiUa: - lgstka[30,2] + ldetYa[30,2] - ldetZa[30,2] >= -1
 lgkiU: - lGstk[30,3] + lDetY[30,3] <= 0
 lgkiU: lDetY[30,3] - lDetZ[30,3] <= 0
 lgkiU: - lGstk[30,3] + lDetY[30,3] - lDetZ[30,3] >= -1
 lgkiUa: - lgstka[30,3] + ldetYa[30,3] <= 0
 lgkiUa: ldetYa[30,3] - ldetZa[30,3] <= 0
 lgkiUa: - lgstka[30,3] + ldetYa[30,3] - ldetZa[30,3] >= -1
 lgkiU: - lGstk[30,4] + lDetY[30,4] <= 0
 lgkiU: lDetY[30,4] - lDetZ[30,4] <= 0
 lgkiU: - lGstk[30,4] + lDetY[30,4] - lDetZ[30,4] >= -1
 lgkiUa: - lgstka[30,4] + ldetYa[30,4] <= 0
 lgkiUa: ldetYa[30,4] - ldetZa[30,4] <= 0
 lgkiUa: - lgstka[30,4] + ldetYa[30,4] - ldetZa[30,4] >= -1
 lgkiU: - lGstk[30,5] + lDetY[30,5] <= 0
 lgkiU: lDetY[30,5] - lDetZ[30,5] <= 0
 lgkiU: - lGstk[30,5] + lDetY[30,5] - lDetZ[30,5] >= -1
 lgkiUa: - lgstka[30,5] + ldetYa[30,5] <= 0
 lgkiUa: ldetYa[30,5] - ldetZa[30,5] <= 0
 lgkiUa: - lgstka[30,5] + ldetYa[30,5] - ldetZa[30,5] >= -1
 lgkiU: - lGstk[30,6] + lDetY[30,6] <= 0
 lgkiU: lDetY[30,6] - lDetZ[30,6] <= 0
 lgkiU: - lGstk[30,6] + lDetY[30,6] - lDetZ[30,6] >= -1
 lgkiUa: - lgstka[30,6] + ldetYa[30,6] <= 0
 lgkiUa: ldetYa[30,6] - ldetZa[30,6] <= 0
 lgkiUa: - lgstka[30,6] + ldetYa[30,6] - ldetZa[30,6] >= -1
 lgkiU: - lGstk[30,7] + lDetY[30,7] <= 0
 lgkiU: lDetY[30,7] - lDetZ[30,7] <= 0
 lgkiU: - lGstk[30,7] + lDetY[30,7] - lDetZ[30,7] >= -1
 lgkiUa: - lgstka[30,7] + ldetYa[30,7] <= 0
 lgkiUa: ldetYa[30,7] - ldetZa[30,7] <= 0
 lgkiUa: - lgstka[30,7] + ldetYa[30,7] - ldetZa[30,7] >= -1
 lgkiL: lDetY[30,8] - lDetZ[30,8] = 0
 lgkiLa: ldetYa[30,8] - ldetZa[30,8] = 0
 lgkiL: lDetY[30,9] - lDetZ[30,9] = 0
 lgkiLa: ldetYa[30,9] - ldetZa[30,9] = 0
 lgkiL: lDetY[30,10] - lDetZ[30,10] = 0
 lgkiLa: ldetYa[30,10] - ldetZa[30,10] = 0
 lgkiL: lDetY[30,11] - lDetZ[30,11] = 0
 lgkiLa: ldetYa[30,11] - ldetZa[30,11] = 0
 lgkiL: lDetY[30,12] - lDetZ[30,12] = 0
 lgkiLa: ldetYa[30,12] - ldetZa[30,12] = 0
 lgkiL: lDetY[30,13] - lDetZ[30,13] = 0
 lgkiLa: ldetYa[30,13] - ldetZa[30,13] = 0
 lgkiL: lDetY[30,14] - lDetZ[30,14] = 0
 lgkiLa: ldetYa[30,14] - ldetZa[30,14] = 0
 lgkiL: lDetY[30,15] - lDetZ[30,15] = 0
 lgkiLa: ldetYa[30,15] - ldetZa[30,15] = 0
 lgkiU: - lGstk[31,0] + lDetY[31,0] <= 0
 lgkiU: lDetY[31,0] - lDetZ[31,0] <= 0
 lgkiU: - lGstk[31,0] + lDetY[31,0] - lDetZ[31,0] >= -1
 lgkiUa: - lgstka[31,0] + ldetYa[31,0] <= 0
 lgkiUa: ldetYa[31,0] - ldetZa[31,0] <= 0
 lgkiUa: - lgstka[31,0] + ldetYa[31,0] - ldetZa[31,0] >= -1
 lgkiU: - lGstk[31,1] + lDetY[31,1] <= 0
 lgkiU: lDetY[31,1] - lDetZ[31,1] <= 0
 lgkiU: - lGstk[31,1] + lDetY[31,1] - lDetZ[31,1] >= -1
 lgkiUa: - lgstka[31,1] + ldetYa[31,1] <= 0
 lgkiUa: ldetYa[31,1] - ldetZa[31,1] <= 0
 lgkiUa: - lgstka[31,1] + ldetYa[31,1] - ldetZa[31,1] >= -1
 lgkiU: - lGstk[31,2] + lDetY[31,2] <= 0
 lgkiU: lDetY[31,2] - lDetZ[31,2] <= 0
 lgkiU: - lGstk[31,2] + lDetY[31,2] - lDetZ[31,2] >= -1
 lgkiUa: - lgstka[31,2] + ldetYa[31,2] <= 0
 lgkiUa: ldetYa[31,2] - ldetZa[31,2] <= 0
 lgkiUa: - lgstka[31,2] + ldetYa[31,2] - ldetZa[31,2] >= -1
 lgkiU: - lGstk[31,3] + lDetY[31,3] <= 0
 lgkiU: lDetY[31,3] - lDetZ[31,3] <= 0
 lgkiU: - lGstk[31,3] + lDetY[31,3] - lDetZ[31,3] >= -1
 lgkiUa: - lgstka[31,3] + ldetYa[31,3] <= 0
 lgkiUa: ldetYa[31,3] - ldetZa[31,3] <= 0
 lgkiUa: - lgstka[31,3] + ldetYa[31,3] - ldetZa[31,3] >= -1
 lgkiU: - lGstk[31,4] + lDetY[31,4] <= 0
 lgkiU: lDetY[31,4] - lDetZ[31,4] <= 0
 lgkiU: - lGstk[31,4] + lDetY[31,4] - lDetZ[31,4] >= -1
 lgkiUa: - lgstka[31,4] + ldetYa[31,4] <= 0
 lgkiUa: ldetYa[31,4] - ldetZa[31,4] <= 0
 lgkiUa: - lgstka[31,4] + ldetYa[31,4] - ldetZa[31,4] >= -1
 lgkiU: - lGstk[31,5] + lDetY[31,5] <= 0
 lgkiU: lDetY[31,5] - lDetZ[31,5] <= 0
 lgkiU: - lGstk[31,5] + lDetY[31,5] - lDetZ[31,5] >= -1
 lgkiUa: - lgstka[31,5] + ldetYa[31,5] <= 0
 lgkiUa: ldetYa[31,5] - ldetZa[31,5] <= 0
 lgkiUa: - lgstka[31,5] + ldetYa[31,5] - ldetZa[31,5] >= -1
 lgkiU: - lGstk[31,6] + lDetY[31,6] <= 0
 lgkiU: lDetY[31,6] - lDetZ[31,6] <= 0
 lgkiU: - lGstk[31,6] + lDetY[31,6] - lDetZ[31,6] >= -1
 lgkiUa: - lgstka[31,6] + ldetYa[31,6] <= 0
 lgkiUa: ldetYa[31,6] - ldetZa[31,6] <= 0
 lgkiUa: - lgstka[31,6] + ldetYa[31,6] - ldetZa[31,6] >= -1
 lgkiU: - lGstk[31,7] + lDetY[31,7] <= 0
 lgkiU: lDetY[31,7] - lDetZ[31,7] <= 0
 lgkiU: - lGstk[31,7] + lDetY[31,7] - lDetZ[31,7] >= -1
 lgkiUa: - lgstka[31,7] + ldetYa[31,7] <= 0
 lgkiUa: ldetYa[31,7] - ldetZa[31,7] <= 0
 lgkiUa: - lgstka[31,7] + ldetYa[31,7] - ldetZa[31,7] >= -1
 lgkiL: lDetY[31,8] - lDetZ[31,8] = 0
 lgkiLa: ldetYa[31,8] - ldetZa[31,8] = 0
 lgkiL: lDetY[31,9] - lDetZ[31,9] = 0
 lgkiLa: ldetYa[31,9] - ldetZa[31,9] = 0
 lgkiL: lDetY[31,10] - lDetZ[31,10] = 0
 lgkiLa: ldetYa[31,10] - ldetZa[31,10] = 0
 lgkiL: lDetY[31,11] - lDetZ[31,11] = 0
 lgkiLa: ldetYa[31,11] - ldetZa[31,11] = 0
 lgkiL: lDetY[31,12] - lDetZ[31,12] = 0
 lgkiLa: ldetYa[31,12] - ldetZa[31,12] = 0
 lgkiL: lDetY[31,13] - lDetZ[31,13] = 0
 lgkiLa: ldetYa[31,13] - ldetZa[31,13] = 0
 lgkiL: lDetY[31,14] - lDetZ[31,14] = 0
 lgkiLa: ldetYa[31,14] - ldetZa[31,14] = 0
 lgkiL: lDetY[31,15] - lDetZ[31,15] = 0
 lgkiLa: ldetYa[31,15] - ldetZa[31,15] = 0
 lgkiU: - lGstk[32,0] + lDetY[32,0] <= 0
 lgkiU: lDetY[32,0] - lDetZ[32,0] <= 0
 lgkiU: - lGstk[32,0] + lDetY[32,0] - lDetZ[32,0] >= -1
 lgkiUa: - lgstka[32,0] + ldetYa[32,0] <= 0
 lgkiUa: ldetYa[32,0] - ldetZa[32,0] <= 0
 lgkiUa: - lgstka[32,0] + ldetYa[32,0] - ldetZa[32,0] >= -1
 lgkiU: - lGstk[32,1] + lDetY[32,1] <= 0
 lgkiU: lDetY[32,1] - lDetZ[32,1] <= 0
 lgkiU: - lGstk[32,1] + lDetY[32,1] - lDetZ[32,1] >= -1
 lgkiUa: - lgstka[32,1] + ldetYa[32,1] <= 0
 lgkiUa: ldetYa[32,1] - ldetZa[32,1] <= 0
 lgkiUa: - lgstka[32,1] + ldetYa[32,1] - ldetZa[32,1] >= -1
 lgkiU: - lGstk[32,2] + lDetY[32,2] <= 0
 lgkiU: lDetY[32,2] - lDetZ[32,2] <= 0
 lgkiU: - lGstk[32,2] + lDetY[32,2] - lDetZ[32,2] >= -1
 lgkiUa: - lgstka[32,2] + ldetYa[32,2] <= 0
 lgkiUa: ldetYa[32,2] - ldetZa[32,2] <= 0
 lgkiUa: - lgstka[32,2] + ldetYa[32,2] - ldetZa[32,2] >= -1
 lgkiU: - lGstk[32,3] + lDetY[32,3] <= 0
 lgkiU: lDetY[32,3] - lDetZ[32,3] <= 0
 lgkiU: - lGstk[32,3] + lDetY[32,3] - lDetZ[32,3] >= -1
 lgkiUa: - lgstka[32,3] + ldetYa[32,3] <= 0
 lgkiUa: ldetYa[32,3] - ldetZa[32,3] <= 0
 lgkiUa: - lgstka[32,3] + ldetYa[32,3] - ldetZa[32,3] >= -1
 lgkiU: - lGstk[32,4] + lDetY[32,4] <= 0
 lgkiU: lDetY[32,4] - lDetZ[32,4] <= 0
 lgkiU: - lGstk[32,4] + lDetY[32,4] - lDetZ[32,4] >= -1
 lgkiUa: - lgstka[32,4] + ldetYa[32,4] <= 0
 lgkiUa: ldetYa[32,4] - ldetZa[32,4] <= 0
 lgkiUa: - lgstka[32,4] + ldetYa[32,4] - ldetZa[32,4] >= -1
 lgkiU: - lGstk[32,5] + lDetY[32,5] <= 0
 lgkiU: lDetY[32,5] - lDetZ[32,5] <= 0
 lgkiU: - lGstk[32,5] + lDetY[32,5] - lDetZ[32,5] >= -1
 lgkiUa: - lgstka[32,5] + ldetYa[32,5] <= 0
 lgkiUa: ldetYa[32,5] - ldetZa[32,5] <= 0
 lgkiUa: - lgstka[32,5] + ldetYa[32,5] - ldetZa[32,5] >= -1
 lgkiU: - lGstk[32,6] + lDetY[32,6] <= 0
 lgkiU: lDetY[32,6] - lDetZ[32,6] <= 0
 lgkiU: - lGstk[32,6] + lDetY[32,6] - lDetZ[32,6] >= -1
 lgkiUa: - lgstka[32,6] + ldetYa[32,6] <= 0
 lgkiUa: ldetYa[32,6] - ldetZa[32,6] <= 0
 lgkiUa: - lgstka[32,6] + ldetYa[32,6] - ldetZa[32,6] >= -1
 lgkiU: - lGstk[32,7] + lDetY[32,7] <= 0
 lgkiU: lDetY[32,7] - lDetZ[32,7] <= 0
 lgkiU: - lGstk[32,7] + lDetY[32,7] - lDetZ[32,7] >= -1
 lgkiUa: - lgstka[32,7] + ldetYa[32,7] <= 0
 lgkiUa: ldetYa[32,7] - ldetZa[32,7] <= 0
 lgkiUa: - lgstka[32,7] + ldetYa[32,7] - ldetZa[32,7] >= -1
 lgkiL: lDetY[32,8] - lDetZ[32,8] = 0
 lgkiLa: ldetYa[32,8] - ldetZa[32,8] = 0
 lgkiL: lDetY[32,9] - lDetZ[32,9] = 0
 lgkiLa: ldetYa[32,9] - ldetZa[32,9] = 0
 lgkiL: lDetY[32,10] - lDetZ[32,10] = 0
 lgkiLa: ldetYa[32,10] - ldetZa[32,10] = 0
 lgkiL: lDetY[32,11] - lDetZ[32,11] = 0
 lgkiLa: ldetYa[32,11] - ldetZa[32,11] = 0
 lgkiL: lDetY[32,12] - lDetZ[32,12] = 0
 lgkiLa: ldetYa[32,12] - ldetZa[32,12] = 0
 lgkiL: lDetY[32,13] - lDetZ[32,13] = 0
 lgkiLa: ldetYa[32,13] - ldetZa[32,13] = 0
 lgkiL: lDetY[32,14] - lDetZ[32,14] = 0
 lgkiLa: ldetYa[32,14] - ldetZa[32,14] = 0
 lgkiL: lDetY[32,15] - lDetZ[32,15] = 0
 lgkiLa: ldetYa[32,15] - ldetZa[32,15] = 0
 ldetSR: lDetZ[28,0] - lDetW[28,0] = 0
 ldetSRa: ldetZa[28,0] - ldetWa[28,0] = 0
 ldetSR: lDetZ[28,1] - lDetW[28,1] = 0
 ldetSRa: ldetZa[28,1] - ldetWa[28,1] = 0
 ldetSR: lDetZ[28,2] - lDetW[28,2] = 0
 ldetSRa: ldetZa[28,2] - ldetWa[28,2] = 0
 ldetSR: lDetZ[28,3] - lDetW[28,3] = 0
 ldetSRa: ldetZa[28,3] - ldetWa[28,3] = 0
 ldetSR: lDetZ[28,4] - lDetW[28,5] = 0
 ldetSRa: ldetZa[28,4] - ldetWa[28,5] = 0
 ldetSR: lDetZ[28,5] - lDetW[28,6] = 0
 ldetSRa: ldetZa[28,5] - ldetWa[28,6] = 0
 ldetSR: lDetZ[28,6] - lDetW[28,7] = 0
 ldetSRa: ldetZa[28,6] - ldetWa[28,7] = 0
 ldetSR: lDetZ[28,7] - lDetW[28,4] = 0
 ldetSRa: ldetZa[28,7] - ldetWa[28,4] = 0
 ldetSR: lDetZ[28,8] - lDetW[28,10] = 0
 ldetSRa: ldetZa[28,8] - ldetWa[28,10] = 0
 ldetSR: lDetZ[28,9] - lDetW[28,11] = 0
 ldetSRa: ldetZa[28,9] - ldetWa[28,11] = 0
 ldetSR: lDetZ[28,10] - lDetW[28,8] = 0
 ldetSRa: ldetZa[28,10] - ldetWa[28,8] = 0
 ldetSR: lDetZ[28,11] - lDetW[28,9] = 0
 ldetSRa: ldetZa[28,11] - ldetWa[28,9] = 0
 ldetSR: lDetZ[28,12] - lDetW[28,15] = 0
 ldetSRa: ldetZa[28,12] - ldetWa[28,15] = 0
 ldetSR: lDetZ[28,13] - lDetW[28,12] = 0
 ldetSRa: ldetZa[28,13] - ldetWa[28,12] = 0
 ldetSR: lDetZ[28,14] - lDetW[28,13] = 0
 ldetSRa: ldetZa[28,14] - ldetWa[28,13] = 0
 ldetSR: lDetZ[28,15] - lDetW[28,14] = 0
 ldetSRa: ldetZa[28,15] - ldetWa[28,14] = 0
 ldetSR: lDetZ[29,0] - lDetW[29,0] = 0
 ldetSRa: ldetZa[29,0] - ldetWa[29,0] = 0
 ldetSR: lDetZ[29,1] - lDetW[29,1] = 0
 ldetSRa: ldetZa[29,1] - ldetWa[29,1] = 0
 ldetSR: lDetZ[29,2] - lDetW[29,2] = 0
 ldetSRa: ldetZa[29,2] - ldetWa[29,2] = 0
 ldetSR: lDetZ[29,3] - lDetW[29,3] = 0
 ldetSRa: ldetZa[29,3] - ldetWa[29,3] = 0
 ldetSR: lDetZ[29,4] - lDetW[29,5] = 0
 ldetSRa: ldetZa[29,4] - ldetWa[29,5] = 0
 ldetSR: lDetZ[29,5] - lDetW[29,6] = 0
 ldetSRa: ldetZa[29,5] - ldetWa[29,6] = 0
 ldetSR: lDetZ[29,6] - lDetW[29,7] = 0
 ldetSRa: ldetZa[29,6] - ldetWa[29,7] = 0
 ldetSR: lDetZ[29,7] - lDetW[29,4] = 0
 ldetSRa: ldetZa[29,7] - ldetWa[29,4] = 0
 ldetSR: lDetZ[29,8] - lDetW[29,10] = 0
 ldetSRa: ldetZa[29,8] - ldetWa[29,10] = 0
 ldetSR: lDetZ[29,9] - lDetW[29,11] = 0
 ldetSRa: ldetZa[29,9] - ldetWa[29,11] = 0
 ldetSR: lDetZ[29,10] - lDetW[29,8] = 0
 ldetSRa: ldetZa[29,10] - ldetWa[29,8] = 0
 ldetSR: lDetZ[29,11] - lDetW[29,9] = 0
 ldetSRa: ldetZa[29,11] - ldetWa[29,9] = 0
 ldetSR: lDetZ[29,12] - lDetW[29,15] = 0
 ldetSRa: ldetZa[29,12] - ldetWa[29,15] = 0
 ldetSR: lDetZ[29,13] - lDetW[29,12] = 0
 ldetSRa: ldetZa[29,13] - ldetWa[29,12] = 0
 ldetSR: lDetZ[29,14] - lDetW[29,13] = 0
 ldetSRa: ldetZa[29,14] - ldetWa[29,13] = 0
 ldetSR: lDetZ[29,15] - lDetW[29,14] = 0
 ldetSRa: ldetZa[29,15] - ldetWa[29,14] = 0
 ldetSR: lDetZ[30,0] - lDetW[30,0] = 0
 ldetSRa: ldetZa[30,0] - ldetWa[30,0] = 0
 ldetSR: lDetZ[30,1] - lDetW[30,1] = 0
 ldetSRa: ldetZa[30,1] - ldetWa[30,1] = 0
 ldetSR: lDetZ[30,2] - lDetW[30,2] = 0
 ldetSRa: ldetZa[30,2] - ldetWa[30,2] = 0
 ldetSR: lDetZ[30,3] - lDetW[30,3] = 0
 ldetSRa: ldetZa[30,3] - ldetWa[30,3] = 0
 ldetSR: lDetZ[30,4] - lDetW[30,5] = 0
 ldetSRa: ldetZa[30,4] - ldetWa[30,5] = 0
 ldetSR: lDetZ[30,5] - lDetW[30,6] = 0
 ldetSRa: ldetZa[30,5] - ldetWa[30,6] = 0
 ldetSR: lDetZ[30,6] - lDetW[30,7] = 0
 ldetSRa: ldetZa[30,6] - ldetWa[30,7] = 0
 ldetSR: lDetZ[30,7] - lDetW[30,4] = 0
 ldetSRa: ldetZa[30,7] - ldetWa[30,4] = 0
 ldetSR: lDetZ[30,8] - lDetW[30,10] = 0
 ldetSRa: ldetZa[30,8] - ldetWa[30,10] = 0
 ldetSR: lDetZ[30,9] - lDetW[30,11] = 0
 ldetSRa: ldetZa[30,9] - ldetWa[30,11] = 0
 ldetSR: lDetZ[30,10] - lDetW[30,8] = 0
 ldetSRa: ldetZa[30,10] - ldetWa[30,8] = 0
 ldetSR: lDetZ[30,11] - lDetW[30,9] = 0
 ldetSRa: ldetZa[30,11] - ldetWa[30,9] = 0
 ldetSR: lDetZ[30,12] - lDetW[30,15] = 0
 ldetSRa: ldetZa[30,12] - ldetWa[30,15] = 0
 ldetSR: lDetZ[30,13] - lDetW[30,12] = 0
 ldetSRa: ldetZa[30,13] - ldetWa[30,12] = 0
 ldetSR: lDetZ[30,14] - lDetW[30,13] = 0
 ldetSRa: ldetZa[30,14] - ldetWa[30,13] = 0
 ldetSR: lDetZ[30,15] - lDetW[30,14] = 0
 ldetSRa: ldetZa[30,15] - ldetWa[30,14] = 0
 ldetSR: lDetZ[31,0] - lDetW[31,0] = 0
 ldetSRa: ldetZa[31,0] - ldetWa[31,0] = 0
 ldetSR: lDetZ[31,1] - lDetW[31,1] = 0
 ldetSRa: ldetZa[31,1] - ldetWa[31,1] = 0
 ldetSR: lDetZ[31,2] - lDetW[31,2] = 0
 ldetSRa: ldetZa[31,2] - ldetWa[31,2] = 0
 ldetSR: lDetZ[31,3] - lDetW[31,3] = 0
 ldetSRa: ldetZa[31,3] - ldetWa[31,3] = 0
 ldetSR: lDetZ[31,4] - lDetW[31,5] = 0
 ldetSRa: ldetZa[31,4] - ldetWa[31,5] = 0
 ldetSR: lDetZ[31,5] - lDetW[31,6] = 0
 ldetSRa: ldetZa[31,5] - ldetWa[31,6] = 0
 ldetSR: lDetZ[31,6] - lDetW[31,7] = 0
 ldetSRa: ldetZa[31,6] - ldetWa[31,7] = 0
 ldetSR: lDetZ[31,7] - lDetW[31,4] = 0
 ldetSRa: ldetZa[31,7] - ldetWa[31,4] = 0
 ldetSR: lDetZ[31,8] - lDetW[31,10] = 0
 ldetSRa: ldetZa[31,8] - ldetWa[31,10] = 0
 ldetSR: lDetZ[31,9] - lDetW[31,11] = 0
 ldetSRa: ldetZa[31,9] - ldetWa[31,11] = 0
 ldetSR: lDetZ[31,10] - lDetW[31,8] = 0
 ldetSRa: ldetZa[31,10] - ldetWa[31,8] = 0
 ldetSR: lDetZ[31,11] - lDetW[31,9] = 0
 ldetSRa: ldetZa[31,11] - ldetWa[31,9] = 0
 ldetSR: lDetZ[31,12] - lDetW[31,15] = 0
 ldetSRa: ldetZa[31,12] - ldetWa[31,15] = 0
 ldetSR: lDetZ[31,13] - lDetW[31,12] = 0
 ldetSRa: ldetZa[31,13] - ldetWa[31,12] = 0
 ldetSR: lDetZ[31,14] - lDetW[31,13] = 0
 ldetSRa: ldetZa[31,14] - ldetWa[31,13] = 0
 ldetSR: lDetZ[31,15] - lDetW[31,14] = 0
 ldetSRa: ldetZa[31,15] - ldetWa[31,14] = 0
 ldetMC0: - lDetX[29,4] + lDetW[28,0] = 0
 ldetMC1: lDetW[28,4] - lDetW[28,8] <= 0
 ldetMC1: - lDetX[29,8] + lDetW[28,4] <= 0
 ldetMC1: - lDetX[29,8] + lDetW[28,4] - lDetW[28,8] >= -1
 ldetMC2: - lDetX[29,4] + lDetW[28,8] <= 0
 ldetMC2: - lDetX[29,12] + lDetW[28,8] <= 0
 ldetMC2: - lDetX[29,4] - lDetX[29,12] + lDetW[28,8] >= -1
 ldetMC3: - lDetX[29,0] + lDetW[28,12] <= 0
 ldetMC3: - lDetX[29,12] + lDetW[28,12] <= 0
 ldetMC3: - lDetX[29,0] - lDetX[29,12] + lDetW[28,12] >= -1
 ldetMC0a: - ldetXa[29,4] + ldetWa[28,0] = 0
 ldetMC1a: ldetWa[28,4] - ldetWa[28,8] <= 0
 ldetMC1a: - ldetXa[29,8] + ldetWa[28,4] <= 0
 ldetMC1a: - ldetXa[29,8] + ldetWa[28,4] - ldetWa[28,8] >= -1
 ldetMC2a: - ldetXa[29,4] + ldetWa[28,8] <= 0
 ldetMC2a: - ldetXa[29,12] + ldetWa[28,8] <= 0
 ldetMC2a: - ldetXa[29,4] - ldetXa[29,12] + ldetWa[28,8] >= -1
 ldetMC3a: - ldetXa[29,0] + ldetWa[28,12] <= 0
 ldetMC3a: - ldetXa[29,12] + ldetWa[28,12] <= 0
 ldetMC3a: - ldetXa[29,0] - ldetXa[29,12] + ldetWa[28,12] >= -1
 ldetMC0: - lDetX[29,5] + lDetW[28,1] = 0
 ldetMC1: lDetW[28,5] - lDetW[28,9] <= 0
 ldetMC1: - lDetX[29,9] + lDetW[28,5] <= 0
 ldetMC1: - lDetX[29,9] + lDetW[28,5] - lDetW[28,9] >= -1
 ldetMC2: - lDetX[29,5] + lDetW[28,9] <= 0
 ldetMC2: - lDetX[29,13] + lDetW[28,9] <= 0
 ldetMC2: - lDetX[29,5] - lDetX[29,13] + lDetW[28,9] >= -1
 ldetMC3: - lDetX[29,1] + lDetW[28,13] <= 0
 ldetMC3: - lDetX[29,13] + lDetW[28,13] <= 0
 ldetMC3: - lDetX[29,1] - lDetX[29,13] + lDetW[28,13] >= -1
 ldetMC0a: - ldetXa[29,5] + ldetWa[28,1] = 0
 ldetMC1a: ldetWa[28,5] - ldetWa[28,9] <= 0
 ldetMC1a: - ldetXa[29,9] + ldetWa[28,5] <= 0
 ldetMC1a: - ldetXa[29,9] + ldetWa[28,5] - ldetWa[28,9] >= -1
 ldetMC2a: - ldetXa[29,5] + ldetWa[28,9] <= 0
 ldetMC2a: - ldetXa[29,13] + ldetWa[28,9] <= 0
 ldetMC2a: - ldetXa[29,5] - ldetXa[29,13] + ldetWa[28,9] >= -1
 ldetMC3a: - ldetXa[29,1] + ldetWa[28,13] <= 0
 ldetMC3a: - ldetXa[29,13] + ldetWa[28,13] <= 0
 ldetMC3a: - ldetXa[29,1] - ldetXa[29,13] + ldetWa[28,13] >= -1
 ldetMC0: - lDetX[29,6] + lDetW[28,2] = 0
 ldetMC1: lDetW[28,6] - lDetW[28,10] <= 0
 ldetMC1: - lDetX[29,10] + lDetW[28,6] <= 0
 ldetMC1: - lDetX[29,10] + lDetW[28,6] - lDetW[28,10] >= -1
 ldetMC2: - lDetX[29,6] + lDetW[28,10] <= 0
 ldetMC2: - lDetX[29,14] + lDetW[28,10] <= 0
 ldetMC2: - lDetX[29,6] - lDetX[29,14] + lDetW[28,10] >= -1
 ldetMC3: - lDetX[29,2] + lDetW[28,14] <= 0
 ldetMC3: - lDetX[29,14] + lDetW[28,14] <= 0
 ldetMC3: - lDetX[29,2] - lDetX[29,14] + lDetW[28,14] >= -1
 ldetMC0a: - ldetXa[29,6] + ldetWa[28,2] = 0
 ldetMC1a: ldetWa[28,6] - ldetWa[28,10] <= 0
 ldetMC1a: - ldetXa[29,10] + ldetWa[28,6] <= 0
 ldetMC1a: - ldetXa[29,10] + ldetWa[28,6] - ldetWa[28,10] >= -1
 ldetMC2a: - ldetXa[29,6] + ldetWa[28,10] <= 0
 ldetMC2a: - ldetXa[29,14] + ldetWa[28,10] <= 0
 ldetMC2a: - ldetXa[29,6] - ldetXa[29,14] + ldetWa[28,10] >= -1
 ldetMC3a: - ldetXa[29,2] + ldetWa[28,14] <= 0
 ldetMC3a: - ldetXa[29,14] + ldetWa[28,14] <= 0
 ldetMC3a: - ldetXa[29,2] - ldetXa[29,14] + ldetWa[28,14] >= -1
 ldetMC0: - lDetX[29,7] + lDetW[28,3] = 0
 ldetMC1: lDetW[28,7] - lDetW[28,11] <= 0
 ldetMC1: - lDetX[29,11] + lDetW[28,7] <= 0
 ldetMC1: - lDetX[29,11] + lDetW[28,7] - lDetW[28,11] >= -1
 ldetMC2: - lDetX[29,7] + lDetW[28,11] <= 0
 ldetMC2: - lDetX[29,15] + lDetW[28,11] <= 0
 ldetMC2: - lDetX[29,7] - lDetX[29,15] + lDetW[28,11] >= -1
 ldetMC3: - lDetX[29,3] + lDetW[28,15] <= 0
 ldetMC3: - lDetX[29,15] + lDetW[28,15] <= 0
 ldetMC3: - lDetX[29,3] - lDetX[29,15] + lDetW[28,15] >= -1
 ldetMC0a: - ldetXa[29,7] + ldetWa[28,3] = 0
 ldetMC1a: ldetWa[28,7] - ldetWa[28,11] <= 0
 ldetMC1a: - ldetXa[29,11] + ldetWa[28,7] <= 0
 ldetMC1a: - ldetXa[29,11] + ldetWa[28,7] - ldetWa[28,11] >= -1
 ldetMC2a: - ldetXa[29,7] + ldetWa[28,11] <= 0
 ldetMC2a: - ldetXa[29,15] + ldetWa[28,11] <= 0
 ldetMC2a: - ldetXa[29,7] - ldetXa[29,15] + ldetWa[28,11] >= -1
 ldetMC3a: - ldetXa[29,3] + ldetWa[28,15] <= 0
 ldetMC3a: - ldetXa[29,15] + ldetWa[28,15] <= 0
 ldetMC3a: - ldetXa[29,3] - ldetXa[29,15] + ldetWa[28,15] >= -1
 ldetMC0: - lDetX[30,4] + lDetW[29,0] = 0
 ldetMC1: lDetW[29,4] - lDetW[29,8] <= 0
 ldetMC1: - lDetX[30,8] + lDetW[29,4] <= 0
 ldetMC1: - lDetX[30,8] + lDetW[29,4] - lDetW[29,8] >= -1
 ldetMC2: - lDetX[30,4] + lDetW[29,8] <= 0
 ldetMC2: - lDetX[30,12] + lDetW[29,8] <= 0
 ldetMC2: - lDetX[30,4] - lDetX[30,12] + lDetW[29,8] >= -1
 ldetMC3: - lDetX[30,0] + lDetW[29,12] <= 0
 ldetMC3: - lDetX[30,12] + lDetW[29,12] <= 0
 ldetMC3: - lDetX[30,0] - lDetX[30,12] + lDetW[29,12] >= -1
 ldetMC0a: - ldetXa[30,4] + ldetWa[29,0] = 0
 ldetMC1a: ldetWa[29,4] - ldetWa[29,8] <= 0
 ldetMC1a: - ldetXa[30,8] + ldetWa[29,4] <= 0
 ldetMC1a: - ldetXa[30,8] + ldetWa[29,4] - ldetWa[29,8] >= -1
 ldetMC2a: - ldetXa[30,4] + ldetWa[29,8] <= 0
 ldetMC2a: - ldetXa[30,12] + ldetWa[29,8] <= 0
 ldetMC2a: - ldetXa[30,4] - ldetXa[30,12] + ldetWa[29,8] >= -1
 ldetMC3a: - ldetXa[30,0] + ldetWa[29,12] <= 0
 ldetMC3a: - ldetXa[30,12] + ldetWa[29,12] <= 0
 ldetMC3a: - ldetXa[30,0] - ldetXa[30,12] + ldetWa[29,12] >= -1
 ldetMC0: - lDetX[30,5] + lDetW[29,1] = 0
 ldetMC1: lDetW[29,5] - lDetW[29,9] <= 0
 ldetMC1: - lDetX[30,9] + lDetW[29,5] <= 0
 ldetMC1: - lDetX[30,9] + lDetW[29,5] - lDetW[29,9] >= -1
 ldetMC2: - lDetX[30,5] + lDetW[29,9] <= 0
 ldetMC2: - lDetX[30,13] + lDetW[29,9] <= 0
 ldetMC2: - lDetX[30,5] - lDetX[30,13] + lDetW[29,9] >= -1
 ldetMC3: - lDetX[30,1] + lDetW[29,13] <= 0
 ldetMC3: - lDetX[30,13] + lDetW[29,13] <= 0
 ldetMC3: - lDetX[30,1] - lDetX[30,13] + lDetW[29,13] >= -1
 ldetMC0a: - ldetXa[30,5] + ldetWa[29,1] = 0
 ldetMC1a: ldetWa[29,5] - ldetWa[29,9] <= 0
 ldetMC1a: - ldetXa[30,9] + ldetWa[29,5] <= 0
 ldetMC1a: - ldetXa[30,9] + ldetWa[29,5] - ldetWa[29,9] >= -1
 ldetMC2a: - ldetXa[30,5] + ldetWa[29,9] <= 0
 ldetMC2a: - ldetXa[30,13] + ldetWa[29,9] <= 0
 ldetMC2a: - ldetXa[30,5] - ldetXa[30,13] + ldetWa[29,9] >= -1
 ldetMC3a: - ldetXa[30,1] + ldetWa[29,13] <= 0
 ldetMC3a: - ldetXa[30,13] + ldetWa[29,13] <= 0
 ldetMC3a: - ldetXa[30,1] - ldetXa[30,13] + ldetWa[29,13] >= -1
 ldetMC0: - lDetX[30,6] + lDetW[29,2] = 0
 ldetMC1: lDetW[29,6] - lDetW[29,10] <= 0
 ldetMC1: - lDetX[30,10] + lDetW[29,6] <= 0
 ldetMC1: - lDetX[30,10] + lDetW[29,6] - lDetW[29,10] >= -1
 ldetMC2: - lDetX[30,6] + lDetW[29,10] <= 0
 ldetMC2: - lDetX[30,14] + lDetW[29,10] <= 0
 ldetMC2: - lDetX[30,6] - lDetX[30,14] + lDetW[29,10] >= -1
 ldetMC3: - lDetX[30,2] + lDetW[29,14] <= 0
 ldetMC3: - lDetX[30,14] + lDetW[29,14] <= 0
 ldetMC3: - lDetX[30,2] - lDetX[30,14] + lDetW[29,14] >= -1
 ldetMC0a: - ldetXa[30,6] + ldetWa[29,2] = 0
 ldetMC1a: ldetWa[29,6] - ldetWa[29,10] <= 0
 ldetMC1a: - ldetXa[30,10] + ldetWa[29,6] <= 0
 ldetMC1a: - ldetXa[30,10] + ldetWa[29,6] - ldetWa[29,10] >= -1
 ldetMC2a: - ldetXa[30,6] + ldetWa[29,10] <= 0
 ldetMC2a: - ldetXa[30,14] + ldetWa[29,10] <= 0
 ldetMC2a: - ldetXa[30,6] - ldetXa[30,14] + ldetWa[29,10] >= -1
 ldetMC3a: - ldetXa[30,2] + ldetWa[29,14] <= 0
 ldetMC3a: - ldetXa[30,14] + ldetWa[29,14] <= 0
 ldetMC3a: - ldetXa[30,2] - ldetXa[30,14] + ldetWa[29,14] >= -1
 ldetMC0: - lDetX[30,7] + lDetW[29,3] = 0
 ldetMC1: lDetW[29,7] - lDetW[29,11] <= 0
 ldetMC1: - lDetX[30,11] + lDetW[29,7] <= 0
 ldetMC1: - lDetX[30,11] + lDetW[29,7] - lDetW[29,11] >= -1
 ldetMC2: - lDetX[30,7] + lDetW[29,11] <= 0
 ldetMC2: - lDetX[30,15] + lDetW[29,11] <= 0
 ldetMC2: - lDetX[30,7] - lDetX[30,15] + lDetW[29,11] >= -1
 ldetMC3: - lDetX[30,3] + lDetW[29,15] <= 0
 ldetMC3: - lDetX[30,15] + lDetW[29,15] <= 0
 ldetMC3: - lDetX[30,3] - lDetX[30,15] + lDetW[29,15] >= -1
 ldetMC0a: - ldetXa[30,7] + ldetWa[29,3] = 0
 ldetMC1a: ldetWa[29,7] - ldetWa[29,11] <= 0
 ldetMC1a: - ldetXa[30,11] + ldetWa[29,7] <= 0
 ldetMC1a: - ldetXa[30,11] + ldetWa[29,7] - ldetWa[29,11] >= -1
 ldetMC2a: - ldetXa[30,7] + ldetWa[29,11] <= 0
 ldetMC2a: - ldetXa[30,15] + ldetWa[29,11] <= 0
 ldetMC2a: - ldetXa[30,7] - ldetXa[30,15] + ldetWa[29,11] >= -1
 ldetMC3a: - ldetXa[30,3] + ldetWa[29,15] <= 0
 ldetMC3a: - ldetXa[30,15] + ldetWa[29,15] <= 0
 ldetMC3a: - ldetXa[30,3] - ldetXa[30,15] + ldetWa[29,15] >= -1
 ldetMC0: - lDetX[31,4] + lDetW[30,0] = 0
 ldetMC1: lDetW[30,4] - lDetW[30,8] <= 0
 ldetMC1: - lDetX[31,8] + lDetW[30,4] <= 0
 ldetMC1: - lDetX[31,8] + lDetW[30,4] - lDetW[30,8] >= -1
 ldetMC2: - lDetX[31,4] + lDetW[30,8] <= 0
 ldetMC2: - lDetX[31,12] + lDetW[30,8] <= 0
 ldetMC2: - lDetX[31,4] - lDetX[31,12] + lDetW[30,8] >= -1
 ldetMC3: - lDetX[31,0] + lDetW[30,12] <= 0
 ldetMC3: - lDetX[31,12] + lDetW[30,12] <= 0
 ldetMC3: - lDetX[31,0] - lDetX[31,12] + lDetW[30,12] >= -1
 ldetMC0a: - ldetXa[31,4] + ldetWa[30,0] = 0
 ldetMC1a: ldetWa[30,4] - ldetWa[30,8] <= 0
 ldetMC1a: - ldetXa[31,8] + ldetWa[30,4] <= 0
 ldetMC1a: - ldetXa[31,8] + ldetWa[30,4] - ldetWa[30,8] >= -1
 ldetMC2a: - ldetXa[31,4] + ldetWa[30,8] <= 0
 ldetMC2a: - ldetXa[31,12] + ldetWa[30,8] <= 0
 ldetMC2a: - ldetXa[31,4] - ldetXa[31,12] + ldetWa[30,8] >= -1
 ldetMC3a: - ldetXa[31,0] + ldetWa[30,12] <= 0
 ldetMC3a: - ldetXa[31,12] + ldetWa[30,12] <= 0
 ldetMC3a: - ldetXa[31,0] - ldetXa[31,12] + ldetWa[30,12] >= -1
 ldetMC0: - lDetX[31,5] + lDetW[30,1] = 0
 ldetMC1: lDetW[30,5] - lDetW[30,9] <= 0
 ldetMC1: - lDetX[31,9] + lDetW[30,5] <= 0
 ldetMC1: - lDetX[31,9] + lDetW[30,5] - lDetW[30,9] >= -1
 ldetMC2: - lDetX[31,5] + lDetW[30,9] <= 0
 ldetMC2: - lDetX[31,13] + lDetW[30,9] <= 0
 ldetMC2: - lDetX[31,5] - lDetX[31,13] + lDetW[30,9] >= -1
 ldetMC3: - lDetX[31,1] + lDetW[30,13] <= 0
 ldetMC3: - lDetX[31,13] + lDetW[30,13] <= 0
 ldetMC3: - lDetX[31,1] - lDetX[31,13] + lDetW[30,13] >= -1
 ldetMC0a: - ldetXa[31,5] + ldetWa[30,1] = 0
 ldetMC1a: ldetWa[30,5] - ldetWa[30,9] <= 0
 ldetMC1a: - ldetXa[31,9] + ldetWa[30,5] <= 0
 ldetMC1a: - ldetXa[31,9] + ldetWa[30,5] - ldetWa[30,9] >= -1
 ldetMC2a: - ldetXa[31,5] + ldetWa[30,9] <= 0
 ldetMC2a: - ldetXa[31,13] + ldetWa[30,9] <= 0
 ldetMC2a: - ldetXa[31,5] - ldetXa[31,13] + ldetWa[30,9] >= -1
 ldetMC3a: - ldetXa[31,1] + ldetWa[30,13] <= 0
 ldetMC3a: - ldetXa[31,13] + ldetWa[30,13] <= 0
 ldetMC3a: - ldetXa[31,1] - ldetXa[31,13] + ldetWa[30,13] >= -1
 ldetMC0: - lDetX[31,6] + lDetW[30,2] = 0
 ldetMC1: lDetW[30,6] - lDetW[30,10] <= 0
 ldetMC1: - lDetX[31,10] + lDetW[30,6] <= 0
 ldetMC1: - lDetX[31,10] + lDetW[30,6] - lDetW[30,10] >= -1
 ldetMC2: - lDetX[31,6] + lDetW[30,10] <= 0
 ldetMC2: - lDetX[31,14] + lDetW[30,10] <= 0
 ldetMC2: - lDetX[31,6] - lDetX[31,14] + lDetW[30,10] >= -1
 ldetMC3: - lDetX[31,2] + lDetW[30,14] <= 0
 ldetMC3: - lDetX[31,14] + lDetW[30,14] <= 0
 ldetMC3: - lDetX[31,2] - lDetX[31,14] + lDetW[30,14] >= -1
 ldetMC0a: - ldetXa[31,6] + ldetWa[30,2] = 0
 ldetMC1a: ldetWa[30,6] - ldetWa[30,10] <= 0
 ldetMC1a: - ldetXa[31,10] + ldetWa[30,6] <= 0
 ldetMC1a: - ldetXa[31,10] + ldetWa[30,6] - ldetWa[30,10] >= -1
 ldetMC2a: - ldetXa[31,6] + ldetWa[30,10] <= 0
 ldetMC2a: - ldetXa[31,14] + ldetWa[30,10] <= 0
 ldetMC2a: - ldetXa[31,6] - ldetXa[31,14] + ldetWa[30,10] >= -1
 ldetMC3a: - ldetXa[31,2] + ldetWa[30,14] <= 0
 ldetMC3a: - ldetXa[31,14] + ldetWa[30,14] <= 0
 ldetMC3a: - ldetXa[31,2] - ldetXa[31,14] + ldetWa[30,14] >= -1
 ldetMC0: - lDetX[31,7] + lDetW[30,3] = 0
 ldetMC1: lDetW[30,7] - lDetW[30,11] <= 0
 ldetMC1: - lDetX[31,11] + lDetW[30,7] <= 0
 ldetMC1: - lDetX[31,11] + lDetW[30,7] - lDetW[30,11] >= -1
 ldetMC2: - lDetX[31,7] + lDetW[30,11] <= 0
 ldetMC2: - lDetX[31,15] + lDetW[30,11] <= 0
 ldetMC2: - lDetX[31,7] - lDetX[31,15] + lDetW[30,11] >= -1
 ldetMC3: - lDetX[31,3] + lDetW[30,15] <= 0
 ldetMC3: - lDetX[31,15] + lDetW[30,15] <= 0
 ldetMC3: - lDetX[31,3] - lDetX[31,15] + lDetW[30,15] >= -1
 ldetMC0a: - ldetXa[31,7] + ldetWa[30,3] = 0
 ldetMC1a: ldetWa[30,7] - ldetWa[30,11] <= 0
 ldetMC1a: - ldetXa[31,11] + ldetWa[30,7] <= 0
 ldetMC1a: - ldetXa[31,11] + ldetWa[30,7] - ldetWa[30,11] >= -1
 ldetMC2a: - ldetXa[31,7] + ldetWa[30,11] <= 0
 ldetMC2a: - ldetXa[31,15] + ldetWa[30,11] <= 0
 ldetMC2a: - ldetXa[31,7] - ldetXa[31,15] + ldetWa[30,11] >= -1
 ldetMC3a: - ldetXa[31,3] + ldetWa[30,15] <= 0
 ldetMC3a: - ldetXa[31,15] + ldetWa[30,15] <= 0
 ldetMC3a: - ldetXa[31,3] - ldetXa[31,15] + ldetWa[30,15] >= -1
 ldetMC0: - lDetX[32,4] + lDetW[31,0] = 0
 ldetMC1: lDetW[31,4] - lDetW[31,8] <= 0
 ldetMC1: - lDetX[32,8] + lDetW[31,4] <= 0
 ldetMC1: - lDetX[32,8] + lDetW[31,4] - lDetW[31,8] >= -1
 ldetMC2: - lDetX[32,4] + lDetW[31,8] <= 0
 ldetMC2: - lDetX[32,12] + lDetW[31,8] <= 0
 ldetMC2: - lDetX[32,4] - lDetX[32,12] + lDetW[31,8] >= -1
 ldetMC3: - lDetX[32,0] + lDetW[31,12] <= 0
 ldetMC3: - lDetX[32,12] + lDetW[31,12] <= 0
 ldetMC3: - lDetX[32,0] - lDetX[32,12] + lDetW[31,12] >= -1
 ldetMC0a: - ldetXa[32,4] + ldetWa[31,0] = 0
 ldetMC1a: ldetWa[31,4] - ldetWa[31,8] <= 0
 ldetMC1a: - ldetXa[32,8] + ldetWa[31,4] <= 0
 ldetMC1a: - ldetXa[32,8] + ldetWa[31,4] - ldetWa[31,8] >= -1
 ldetMC2a: - ldetXa[32,4] + ldetWa[31,8] <= 0
 ldetMC2a: - ldetXa[32,12] + ldetWa[31,8] <= 0
 ldetMC2a: - ldetXa[32,4] - ldetXa[32,12] + ldetWa[31,8] >= -1
 ldetMC3a: - ldetXa[32,0] + ldetWa[31,12] <= 0
 ldetMC3a: - ldetXa[32,12] + ldetWa[31,12] <= 0
 ldetMC3a: - ldetXa[32,0] - ldetXa[32,12] + ldetWa[31,12] >= -1
 ldetMC0: - lDetX[32,5] + lDetW[31,1] = 0
 ldetMC1: lDetW[31,5] - lDetW[31,9] <= 0
 ldetMC1: - lDetX[32,9] + lDetW[31,5] <= 0
 ldetMC1: - lDetX[32,9] + lDetW[31,5] - lDetW[31,9] >= -1
 ldetMC2: - lDetX[32,5] + lDetW[31,9] <= 0
 ldetMC2: - lDetX[32,13] + lDetW[31,9] <= 0
 ldetMC2: - lDetX[32,5] - lDetX[32,13] + lDetW[31,9] >= -1
 ldetMC3: - lDetX[32,1] + lDetW[31,13] <= 0
 ldetMC3: - lDetX[32,13] + lDetW[31,13] <= 0
 ldetMC3: - lDetX[32,1] - lDetX[32,13] + lDetW[31,13] >= -1
 ldetMC0a: - ldetXa[32,5] + ldetWa[31,1] = 0
 ldetMC1a: ldetWa[31,5] - ldetWa[31,9] <= 0
 ldetMC1a: - ldetXa[32,9] + ldetWa[31,5] <= 0
 ldetMC1a: - ldetXa[32,9] + ldetWa[31,5] - ldetWa[31,9] >= -1
 ldetMC2a: - ldetXa[32,5] + ldetWa[31,9] <= 0
 ldetMC2a: - ldetXa[32,13] + ldetWa[31,9] <= 0
 ldetMC2a: - ldetXa[32,5] - ldetXa[32,13] + ldetWa[31,9] >= -1
 ldetMC3a: - ldetXa[32,1] + ldetWa[31,13] <= 0
 ldetMC3a: - ldetXa[32,13] + ldetWa[31,13] <= 0
 ldetMC3a: - ldetXa[32,1] - ldetXa[32,13] + ldetWa[31,13] >= -1
 ldetMC0: - lDetX[32,6] + lDetW[31,2] = 0
 ldetMC1: lDetW[31,6] - lDetW[31,10] <= 0
 ldetMC1: - lDetX[32,10] + lDetW[31,6] <= 0
 ldetMC1: - lDetX[32,10] + lDetW[31,6] - lDetW[31,10] >= -1
 ldetMC2: - lDetX[32,6] + lDetW[31,10] <= 0
 ldetMC2: - lDetX[32,14] + lDetW[31,10] <= 0
 ldetMC2: - lDetX[32,6] - lDetX[32,14] + lDetW[31,10] >= -1
 ldetMC3: - lDetX[32,2] + lDetW[31,14] <= 0
 ldetMC3: - lDetX[32,14] + lDetW[31,14] <= 0
 ldetMC3: - lDetX[32,2] - lDetX[32,14] + lDetW[31,14] >= -1
 ldetMC0a: - ldetXa[32,6] + ldetWa[31,2] = 0
 ldetMC1a: ldetWa[31,6] - ldetWa[31,10] <= 0
 ldetMC1a: - ldetXa[32,10] + ldetWa[31,6] <= 0
 ldetMC1a: - ldetXa[32,10] + ldetWa[31,6] - ldetWa[31,10] >= -1
 ldetMC2a: - ldetXa[32,6] + ldetWa[31,10] <= 0
 ldetMC2a: - ldetXa[32,14] + ldetWa[31,10] <= 0
 ldetMC2a: - ldetXa[32,6] - ldetXa[32,14] + ldetWa[31,10] >= -1
 ldetMC3a: - ldetXa[32,2] + ldetWa[31,14] <= 0
 ldetMC3a: - ldetXa[32,14] + ldetWa[31,14] <= 0
 ldetMC3a: - ldetXa[32,2] - ldetXa[32,14] + ldetWa[31,14] >= -1
 ldetMC0: - lDetX[32,7] + lDetW[31,3] = 0
 ldetMC1: lDetW[31,7] - lDetW[31,11] <= 0
 ldetMC1: - lDetX[32,11] + lDetW[31,7] <= 0
 ldetMC1: - lDetX[32,11] + lDetW[31,7] - lDetW[31,11] >= -1
 ldetMC2: - lDetX[32,7] + lDetW[31,11] <= 0
 ldetMC2: - lDetX[32,15] + lDetW[31,11] <= 0
 ldetMC2: - lDetX[32,7] - lDetX[32,15] + lDetW[31,11] >= -1
 ldetMC3: - lDetX[32,3] + lDetW[31,15] <= 0
 ldetMC3: - lDetX[32,15] + lDetW[31,15] <= 0
 ldetMC3: - lDetX[32,3] - lDetX[32,15] + lDetW[31,15] >= -1
 ldetMC0a: - ldetXa[32,7] + ldetWa[31,3] = 0
 ldetMC1a: ldetWa[31,7] - ldetWa[31,11] <= 0
 ldetMC1a: - ldetXa[32,11] + ldetWa[31,7] <= 0
 ldetMC1a: - ldetXa[32,11] + ldetWa[31,7] - ldetWa[31,11] >= -1
 ldetMC2a: - ldetXa[32,7] + ldetWa[31,11] <= 0
 ldetMC2a: - ldetXa[32,15] + ldetWa[31,11] <= 0
 ldetMC2a: - ldetXa[32,7] - ldetXa[32,15] + ldetWa[31,11] >= -1
 ldetMC3a: - ldetXa[32,3] + ldetWa[31,15] <= 0
 ldetMC3a: - ldetXa[32,15] + ldetWa[31,15] <= 0
 ldetMC3a: - ldetXa[32,3] - ldetXa[32,15] + ldetWa[31,15] >= -1
 lfrMC0: lFrMC[28,0] = 0
 lfrMC1: - lDetW[28,4] + lFrMC[28,4] <= 0
 lfrMC1: - ldw[28,8] + lFrMC[28,4] <= 0
 lfrMC1: - ldx[29,8] + lFrMC[28,4] <= 0
 lfrMC1: ldw[28,4] + lFrMC[28,4] <= 1
 lfrMC1: - ldx[29,8] + ldw[28,4] - ldw[28,8] - lDetW[28,4] + lFrMC[28,4]
   >= -2
 lfrMC2: - lDetW[28,8] + lFrMC[28,8] <= 0
 lfrMC2: - ldx[29,4] + lFrMC[28,8] <= 0
 lfrMC2: - ldx[29,12] + lFrMC[28,8] <= 0
 lfrMC2: ldw[28,8] + lFrMC[28,8] <= 1
 lfrMC2: - ldx[29,4] - ldx[29,12] + ldw[28,8] - lDetW[28,8] + lFrMC[28,8]
   >= -2
 lfrMC3: - lDetW[28,12] + lFrMC[28,12] <= 0
 lfrMC3: - ldx[29,0] + lFrMC[28,12] <= 0
 lfrMC3: - ldx[29,12] + lFrMC[28,12] <= 0
 lfrMC3: ldw[28,12] + lFrMC[28,12] <= 1
 lfrMC3: - ldx[29,0] - ldx[29,12] + ldw[28,12] - lDetW[28,12]
   + lFrMC[28,12] >= -2
 lfrMC0a: lfrMCa[28,0] = 0
 lfrMC1a: - ldetWa[28,4] + lfrMCa[28,4] <= 0
 lfrMC1a: - ldw[28,8] + lfrMCa[28,4] <= 0
 lfrMC1a: - ldx[29,8] + lfrMCa[28,4] <= 0
 lfrMC1a: ldw[28,4] + lfrMCa[28,4] <= 1
 lfrMC1a: - ldx[29,8] + ldw[28,4] - ldw[28,8] - ldetWa[28,4] + lfrMCa[28,4]
   >= -2
 lfrMC2a: - ldetWa[28,8] + lfrMCa[28,8] <= 0
 lfrMC2a: - ldx[29,4] + lfrMCa[28,8] <= 0
 lfrMC2a: - ldx[29,12] + lfrMCa[28,8] <= 0
 lfrMC2a: ldw[28,8] + lfrMCa[28,8] <= 1
 lfrMC2a: - ldx[29,4] - ldx[29,12] + ldw[28,8] - ldetWa[28,8]
   + lfrMCa[28,8] >= -2
 lfrMC3a: - ldetWa[28,12] + lfrMCa[28,12] <= 0
 lfrMC3a: - ldx[29,0] + lfrMCa[28,12] <= 0
 lfrMC3a: - ldx[29,12] + lfrMCa[28,12] <= 0
 lfrMC3a: ldw[28,12] + lfrMCa[28,12] <= 1
 lfrMC3a: - ldx[29,0] - ldx[29,12] + ldw[28,12] - ldetWa[28,12]
   + lfrMCa[28,12] >= -2
 lfrMC0: lFrMC[28,1] = 0
 lfrMC1: - lDetW[28,5] + lFrMC[28,5] <= 0
 lfrMC1: - ldw[28,9] + lFrMC[28,5] <= 0
 lfrMC1: - ldx[29,9] + lFrMC[28,5] <= 0
 lfrMC1: ldw[28,5] + lFrMC[28,5] <= 1
 lfrMC1: - ldx[29,9] + ldw[28,5] - ldw[28,9] - lDetW[28,5] + lFrMC[28,5]
   >= -2
 lfrMC2: - lDetW[28,9] + lFrMC[28,9] <= 0
 lfrMC2: - ldx[29,5] + lFrMC[28,9] <= 0
 lfrMC2: - ldx[29,13] + lFrMC[28,9] <= 0
 lfrMC2: ldw[28,9] + lFrMC[28,9] <= 1
 lfrMC2: - ldx[29,5] - ldx[29,13] + ldw[28,9] - lDetW[28,9] + lFrMC[28,9]
   >= -2
 lfrMC3: - lDetW[28,13] + lFrMC[28,13] <= 0
 lfrMC3: - ldx[29,1] + lFrMC[28,13] <= 0
 lfrMC3: - ldx[29,13] + lFrMC[28,13] <= 0
 lfrMC3: ldw[28,13] + lFrMC[28,13] <= 1
 lfrMC3: - ldx[29,1] - ldx[29,13] + ldw[28,13] - lDetW[28,13]
   + lFrMC[28,13] >= -2
 lfrMC0a: lfrMCa[28,1] = 0
 lfrMC1a: - ldetWa[28,5] + lfrMCa[28,5] <= 0
 lfrMC1a: - ldw[28,9] + lfrMCa[28,5] <= 0
 lfrMC1a: - ldx[29,9] + lfrMCa[28,5] <= 0
 lfrMC1a: ldw[28,5] + lfrMCa[28,5] <= 1
 lfrMC1a: - ldx[29,9] + ldw[28,5] - ldw[28,9] - ldetWa[28,5] + lfrMCa[28,5]
   >= -2
 lfrMC2a: - ldetWa[28,9] + lfrMCa[28,9] <= 0
 lfrMC2a: - ldx[29,5] + lfrMCa[28,9] <= 0
 lfrMC2a: - ldx[29,13] + lfrMCa[28,9] <= 0
 lfrMC2a: ldw[28,9] + lfrMCa[28,9] <= 1
 lfrMC2a: - ldx[29,5] - ldx[29,13] + ldw[28,9] - ldetWa[28,9]
   + lfrMCa[28,9] >= -2
 lfrMC3a: - ldetWa[28,13] + lfrMCa[28,13] <= 0
 lfrMC3a: - ldx[29,1] + lfrMCa[28,13] <= 0
 lfrMC3a: - ldx[29,13] + lfrMCa[28,13] <= 0
 lfrMC3a: ldw[28,13] + lfrMCa[28,13] <= 1
 lfrMC3a: - ldx[29,1] - ldx[29,13] + ldw[28,13] - ldetWa[28,13]
   + lfrMCa[28,13] >= -2
 lfrMC0: lFrMC[28,2] = 0
 lfrMC1: - lDetW[28,6] + lFrMC[28,6] <= 0
 lfrMC1: - ldw[28,10] + lFrMC[28,6] <= 0
 lfrMC1: - ldx[29,10] + lFrMC[28,6] <= 0
 lfrMC1: ldw[28,6] + lFrMC[28,6] <= 1
 lfrMC1: - ldx[29,10] + ldw[28,6] - ldw[28,10] - lDetW[28,6] + lFrMC[28,6]
   >= -2
 lfrMC2: - lDetW[28,10] + lFrMC[28,10] <= 0
 lfrMC2: - ldx[29,6] + lFrMC[28,10] <= 0
 lfrMC2: - ldx[29,14] + lFrMC[28,10] <= 0
 lfrMC2: ldw[28,10] + lFrMC[28,10] <= 1
 lfrMC2: - ldx[29,6] - ldx[29,14] + ldw[28,10] - lDetW[28,10]
   + lFrMC[28,10] >= -2
 lfrMC3: - lDetW[28,14] + lFrMC[28,14] <= 0
 lfrMC3: - ldx[29,2] + lFrMC[28,14] <= 0
 lfrMC3: - ldx[29,14] + lFrMC[28,14] <= 0
 lfrMC3: ldw[28,14] + lFrMC[28,14] <= 1
 lfrMC3: - ldx[29,2] - ldx[29,14] + ldw[28,14] - lDetW[28,14]
   + lFrMC[28,14] >= -2
 lfrMC0a: lfrMCa[28,2] = 0
 lfrMC1a: - ldetWa[28,6] + lfrMCa[28,6] <= 0
 lfrMC1a: - ldw[28,10] + lfrMCa[28,6] <= 0
 lfrMC1a: - ldx[29,10] + lfrMCa[28,6] <= 0
 lfrMC1a: ldw[28,6] + lfrMCa[28,6] <= 1
 lfrMC1a: - ldx[29,10] + ldw[28,6] - ldw[28,10] - ldetWa[28,6]
   + lfrMCa[28,6] >= -2
 lfrMC2a: - ldetWa[28,10] + lfrMCa[28,10] <= 0
 lfrMC2a: - ldx[29,6] + lfrMCa[28,10] <= 0
 lfrMC2a: - ldx[29,14] + lfrMCa[28,10] <= 0
 lfrMC2a: ldw[28,10] + lfrMCa[28,10] <= 1
 lfrMC2a: - ldx[29,6] - ldx[29,14] + ldw[28,10] - ldetWa[28,10]
   + lfrMCa[28,10] >= -2
 lfrMC3a: - ldetWa[28,14] + lfrMCa[28,14] <= 0
 lfrMC3a: - ldx[29,2] + lfrMCa[28,14] <= 0
 lfrMC3a: - ldx[29,14] + lfrMCa[28,14] <= 0
 lfrMC3a: ldw[28,14] + lfrMCa[28,14] <= 1
 lfrMC3a: - ldx[29,2] - ldx[29,14] + ldw[28,14] - ldetWa[28,14]
   + lfrMCa[28,14] >= -2
 lfrMC0: lFrMC[28,3] = 0
 lfrMC1: - lDetW[28,7] + lFrMC[28,7] <= 0
 lfrMC1: - ldw[28,11] + lFrMC[28,7] <= 0
 lfrMC1: - ldx[29,11] + lFrMC[28,7] <= 0
 lfrMC1: ldw[28,7] + lFrMC[28,7] <= 1
 lfrMC1: - ldx[29,11] + ldw[28,7] - ldw[28,11] - lDetW[28,7] + lFrMC[28,7]
   >= -2
 lfrMC2: - lDetW[28,11] + lFrMC[28,11] <= 0
 lfrMC2: - ldx[29,7] + lFrMC[28,11] <= 0
 lfrMC2: - ldx[29,15] + lFrMC[28,11] <= 0
 lfrMC2: ldw[28,11] + lFrMC[28,11] <= 1
 lfrMC2: - ldx[29,7] - ldx[29,15] + ldw[28,11] - lDetW[28,11]
   + lFrMC[28,11] >= -2
 lfrMC3: - lDetW[28,15] + lFrMC[28,15] <= 0
 lfrMC3: - ldx[29,3] + lFrMC[28,15] <= 0
 lfrMC3: - ldx[29,15] + lFrMC[28,15] <= 0
 lfrMC3: ldw[28,15] + lFrMC[28,15] <= 1
 lfrMC3: - ldx[29,3] - ldx[29,15] + ldw[28,15] - lDetW[28,15]
   + lFrMC[28,15] >= -2
 lfrMC0a: lfrMCa[28,3] = 0
 lfrMC1a: - ldetWa[28,7] + lfrMCa[28,7] <= 0
 lfrMC1a: - ldw[28,11] + lfrMCa[28,7] <= 0
 lfrMC1a: - ldx[29,11] + lfrMCa[28,7] <= 0
 lfrMC1a: ldw[28,7] + lfrMCa[28,7] <= 1
 lfrMC1a: - ldx[29,11] + ldw[28,7] - ldw[28,11] - ldetWa[28,7]
   + lfrMCa[28,7] >= -2
 lfrMC2a: - ldetWa[28,11] + lfrMCa[28,11] <= 0
 lfrMC2a: - ldx[29,7] + lfrMCa[28,11] <= 0
 lfrMC2a: - ldx[29,15] + lfrMCa[28,11] <= 0
 lfrMC2a: ldw[28,11] + lfrMCa[28,11] <= 1
 lfrMC2a: - ldx[29,7] - ldx[29,15] + ldw[28,11] - ldetWa[28,11]
   + lfrMCa[28,11] >= -2
 lfrMC3a: - ldetWa[28,15] + lfrMCa[28,15] <= 0
 lfrMC3a: - ldx[29,3] + lfrMCa[28,15] <= 0
 lfrMC3a: - ldx[29,15] + lfrMCa[28,15] <= 0
 lfrMC3a: ldw[28,15] + lfrMCa[28,15] <= 1
 lfrMC3a: - ldx[29,3] - ldx[29,15] + ldw[28,15] - ldetWa[28,15]
   + lfrMCa[28,15] >= -2
 lfrMC0: lFrMC[29,0] = 0
 lfrMC1: - lDetW[29,4] + lFrMC[29,4] <= 0
 lfrMC1: - ldw[29,8] + lFrMC[29,4] <= 0
 lfrMC1: - ldx[30,8] + lFrMC[29,4] <= 0
 lfrMC1: ldw[29,4] + lFrMC[29,4] <= 1
 lfrMC1: - ldx[30,8] + ldw[29,4] - ldw[29,8] - lDetW[29,4] + lFrMC[29,4]
   >= -2
 lfrMC2: - lDetW[29,8] + lFrMC[29,8] <= 0
 lfrMC2: - ldx[30,4] + lFrMC[29,8] <= 0
 lfrMC2: - ldx[30,12] + lFrMC[29,8] <= 0
 lfrMC2: ldw[29,8] + lFrMC[29,8] <= 1
 lfrMC2: - ldx[30,4] - ldx[30,12] + ldw[29,8] - lDetW[29,8] + lFrMC[29,8]
   >= -2
 lfrMC3: - lDetW[29,12] + lFrMC[29,12] <= 0
 lfrMC3: - ldx[30,0] + lFrMC[29,12] <= 0
 lfrMC3: - ldx[30,12] + lFrMC[29,12] <= 0
 lfrMC3: ldw[29,12] + lFrMC[29,12] <= 1
 lfrMC3: - ldx[30,0] - ldx[30,12] + ldw[29,12] - lDetW[29,12]
   + lFrMC[29,12] >= -2
 lfrMC0a: lfrMCa[29,0] = 0
 lfrMC1a: - ldetWa[29,4] + lfrMCa[29,4] <= 0
 lfrMC1a: - ldw[29,8] + lfrMCa[29,4] <= 0
 lfrMC1a: - ldx[30,8] + lfrMCa[29,4] <= 0
 lfrMC1a: ldw[29,4] + lfrMCa[29,4] <= 1
 lfrMC1a: - ldx[30,8] + ldw[29,4] - ldw[29,8] - ldetWa[29,4] + lfrMCa[29,4]
   >= -2
 lfrMC2a: - ldetWa[29,8] + lfrMCa[29,8] <= 0
 lfrMC2a: - ldx[30,4] + lfrMCa[29,8] <= 0
 lfrMC2a: - ldx[30,12] + lfrMCa[29,8] <= 0
 lfrMC2a: ldw[29,8] + lfrMCa[29,8] <= 1
 lfrMC2a: - ldx[30,4] - ldx[30,12] + ldw[29,8] - ldetWa[29,8]
   + lfrMCa[29,8] >= -2
 lfrMC3a: - ldetWa[29,12] + lfrMCa[29,12] <= 0
 lfrMC3a: - ldx[30,0] + lfrMCa[29,12] <= 0
 lfrMC3a: - ldx[30,12] + lfrMCa[29,12] <= 0
 lfrMC3a: ldw[29,12] + lfrMCa[29,12] <= 1
 lfrMC3a: - ldx[30,0] - ldx[30,12] + ldw[29,12] - ldetWa[29,12]
   + lfrMCa[29,12] >= -2
 lfrMC0: lFrMC[29,1] = 0
 lfrMC1: - lDetW[29,5] + lFrMC[29,5] <= 0
 lfrMC1: - ldw[29,9] + lFrMC[29,5] <= 0
 lfrMC1: - ldx[30,9] + lFrMC[29,5] <= 0
 lfrMC1: ldw[29,5] + lFrMC[29,5] <= 1
 lfrMC1: - ldx[30,9] + ldw[29,5] - ldw[29,9] - lDetW[29,5] + lFrMC[29,5]
   >= -2
 lfrMC2: - lDetW[29,9] + lFrMC[29,9] <= 0
 lfrMC2: - ldx[30,5] + lFrMC[29,9] <= 0
 lfrMC2: - ldx[30,13] + lFrMC[29,9] <= 0
 lfrMC2: ldw[29,9] + lFrMC[29,9] <= 1
 lfrMC2: - ldx[30,5] - ldx[30,13] + ldw[29,9] - lDetW[29,9] + lFrMC[29,9]
   >= -2
 lfrMC3: - lDetW[29,13] + lFrMC[29,13] <= 0
 lfrMC3: - ldx[30,1] + lFrMC[29,13] <= 0
 lfrMC3: - ldx[30,13] + lFrMC[29,13] <= 0
 lfrMC3: ldw[29,13] + lFrMC[29,13] <= 1
 lfrMC3: - ldx[30,1] - ldx[30,13] + ldw[29,13] - lDetW[29,13]
   + lFrMC[29,13] >= -2
 lfrMC0a: lfrMCa[29,1] = 0
 lfrMC1a: - ldetWa[29,5] + lfrMCa[29,5] <= 0
 lfrMC1a: - ldw[29,9] + lfrMCa[29,5] <= 0
 lfrMC1a: - ldx[30,9] + lfrMCa[29,5] <= 0
 lfrMC1a: ldw[29,5] + lfrMCa[29,5] <= 1
 lfrMC1a: - ldx[30,9] + ldw[29,5] - ldw[29,9] - ldetWa[29,5] + lfrMCa[29,5]
   >= -2
 lfrMC2a: - ldetWa[29,9] + lfrMCa[29,9] <= 0
 lfrMC2a: - ldx[30,5] + lfrMCa[29,9] <= 0
 lfrMC2a: - ldx[30,13] + lfrMCa[29,9] <= 0
 lfrMC2a: ldw[29,9] + lfrMCa[29,9] <= 1
 lfrMC2a: - ldx[30,5] - ldx[30,13] + ldw[29,9] - ldetWa[29,9]
   + lfrMCa[29,9] >= -2
 lfrMC3a: - ldetWa[29,13] + lfrMCa[29,13] <= 0
 lfrMC3a: - ldx[30,1] + lfrMCa[29,13] <= 0
 lfrMC3a: - ldx[30,13] + lfrMCa[29,13] <= 0
 lfrMC3a: ldw[29,13] + lfrMCa[29,13] <= 1
 lfrMC3a: - ldx[30,1] - ldx[30,13] + ldw[29,13] - ldetWa[29,13]
   + lfrMCa[29,13] >= -2
 lfrMC0: lFrMC[29,2] = 0
 lfrMC1: - lDetW[29,6] + lFrMC[29,6] <= 0
 lfrMC1: - ldw[29,10] + lFrMC[29,6] <= 0
 lfrMC1: - ldx[30,10] + lFrMC[29,6] <= 0
 lfrMC1: ldw[29,6] + lFrMC[29,6] <= 1
 lfrMC1: - ldx[30,10] + ldw[29,6] - ldw[29,10] - lDetW[29,6] + lFrMC[29,6]
   >= -2
 lfrMC2: - lDetW[29,10] + lFrMC[29,10] <= 0
 lfrMC2: - ldx[30,6] + lFrMC[29,10] <= 0
 lfrMC2: - ldx[30,14] + lFrMC[29,10] <= 0
 lfrMC2: ldw[29,10] + lFrMC[29,10] <= 1
 lfrMC2: - ldx[30,6] - ldx[30,14] + ldw[29,10] - lDetW[29,10]
   + lFrMC[29,10] >= -2
 lfrMC3: - lDetW[29,14] + lFrMC[29,14] <= 0
 lfrMC3: - ldx[30,2] + lFrMC[29,14] <= 0
 lfrMC3: - ldx[30,14] + lFrMC[29,14] <= 0
 lfrMC3: ldw[29,14] + lFrMC[29,14] <= 1
 lfrMC3: - ldx[30,2] - ldx[30,14] + ldw[29,14] - lDetW[29,14]
   + lFrMC[29,14] >= -2
 lfrMC0a: lfrMCa[29,2] = 0
 lfrMC1a: - ldetWa[29,6] + lfrMCa[29,6] <= 0
 lfrMC1a: - ldw[29,10] + lfrMCa[29,6] <= 0
 lfrMC1a: - ldx[30,10] + lfrMCa[29,6] <= 0
 lfrMC1a: ldw[29,6] + lfrMCa[29,6] <= 1
 lfrMC1a: - ldx[30,10] + ldw[29,6] - ldw[29,10] - ldetWa[29,6]
   + lfrMCa[29,6] >= -2
 lfrMC2a: - ldetWa[29,10] + lfrMCa[29,10] <= 0
 lfrMC2a: - ldx[30,6] + lfrMCa[29,10] <= 0
 lfrMC2a: - ldx[30,14] + lfrMCa[29,10] <= 0
 lfrMC2a: ldw[29,10] + lfrMCa[29,10] <= 1
 lfrMC2a: - ldx[30,6] - ldx[30,14] + ldw[29,10] - ldetWa[29,10]
   + lfrMCa[29,10] >= -2
 lfrMC3a: - ldetWa[29,14] + lfrMCa[29,14] <= 0
 lfrMC3a: - ldx[30,2] + lfrMCa[29,14] <= 0
 lfrMC3a: - ldx[30,14] + lfrMCa[29,14] <= 0
 lfrMC3a: ldw[29,14] + lfrMCa[29,14] <= 1
 lfrMC3a: - ldx[30,2] - ldx[30,14] + ldw[29,14] - ldetWa[29,14]
   + lfrMCa[29,14] >= -2
 lfrMC0: lFrMC[29,3] = 0
 lfrMC1: - lDetW[29,7] + lFrMC[29,7] <= 0
 lfrMC1: - ldw[29,11] + lFrMC[29,7] <= 0
 lfrMC1: - ldx[30,11] + lFrMC[29,7] <= 0
 lfrMC1: ldw[29,7] + lFrMC[29,7] <= 1
 lfrMC1: - ldx[30,11] + ldw[29,7] - ldw[29,11] - lDetW[29,7] + lFrMC[29,7]
   >= -2
 lfrMC2: - lDetW[29,11] + lFrMC[29,11] <= 0
 lfrMC2: - ldx[30,7] + lFrMC[29,11] <= 0
 lfrMC2: - ldx[30,15] + lFrMC[29,11] <= 0
 lfrMC2: ldw[29,11] + lFrMC[29,11] <= 1
 lfrMC2: - ldx[30,7] - ldx[30,15] + ldw[29,11] - lDetW[29,11]
   + lFrMC[29,11] >= -2
 lfrMC3: - lDetW[29,15] + lFrMC[29,15] <= 0
 lfrMC3: - ldx[30,3] + lFrMC[29,15] <= 0
 lfrMC3: - ldx[30,15] + lFrMC[29,15] <= 0
 lfrMC3: ldw[29,15] + lFrMC[29,15] <= 1
 lfrMC3: - ldx[30,3] - ldx[30,15] + ldw[29,15] - lDetW[29,15]
   + lFrMC[29,15] >= -2
 lfrMC0a: lfrMCa[29,3] = 0
 lfrMC1a: - ldetWa[29,7] + lfrMCa[29,7] <= 0
 lfrMC1a: - ldw[29,11] + lfrMCa[29,7] <= 0
 lfrMC1a: - ldx[30,11] + lfrMCa[29,7] <= 0
 lfrMC1a: ldw[29,7] + lfrMCa[29,7] <= 1
 lfrMC1a: - ldx[30,11] + ldw[29,7] - ldw[29,11] - ldetWa[29,7]
   + lfrMCa[29,7] >= -2
 lfrMC2a: - ldetWa[29,11] + lfrMCa[29,11] <= 0
 lfrMC2a: - ldx[30,7] + lfrMCa[29,11] <= 0
 lfrMC2a: - ldx[30,15] + lfrMCa[29,11] <= 0
 lfrMC2a: ldw[29,11] + lfrMCa[29,11] <= 1
 lfrMC2a: - ldx[30,7] - ldx[30,15] + ldw[29,11] - ldetWa[29,11]
   + lfrMCa[29,11] >= -2
 lfrMC3a: - ldetWa[29,15] + lfrMCa[29,15] <= 0
 lfrMC3a: - ldx[30,3] + lfrMCa[29,15] <= 0
 lfrMC3a: - ldx[30,15] + lfrMCa[29,15] <= 0
 lfrMC3a: ldw[29,15] + lfrMCa[29,15] <= 1
 lfrMC3a: - ldx[30,3] - ldx[30,15] + ldw[29,15] - ldetWa[29,15]
   + lfrMCa[29,15] >= -2
 lfrMC0: lFrMC[30,0] = 0
 lfrMC1: - lDetW[30,4] + lFrMC[30,4] <= 0
 lfrMC1: - ldw[30,8] + lFrMC[30,4] <= 0
 lfrMC1: - ldx[31,8] + lFrMC[30,4] <= 0
 lfrMC1: ldw[30,4] + lFrMC[30,4] <= 1
 lfrMC1: - ldx[31,8] + ldw[30,4] - ldw[30,8] - lDetW[30,4] + lFrMC[30,4]
   >= -2
 lfrMC2: - lDetW[30,8] + lFrMC[30,8] <= 0
 lfrMC2: - ldx[31,4] + lFrMC[30,8] <= 0
 lfrMC2: - ldx[31,12] + lFrMC[30,8] <= 0
 lfrMC2: ldw[30,8] + lFrMC[30,8] <= 1
 lfrMC2: - ldx[31,4] - ldx[31,12] + ldw[30,8] - lDetW[30,8] + lFrMC[30,8]
   >= -2
 lfrMC3: - lDetW[30,12] + lFrMC[30,12] <= 0
 lfrMC3: - ldx[31,0] + lFrMC[30,12] <= 0
 lfrMC3: - ldx[31,12] + lFrMC[30,12] <= 0
 lfrMC3: ldw[30,12] + lFrMC[30,12] <= 1
 lfrMC3: - ldx[31,0] - ldx[31,12] + ldw[30,12] - lDetW[30,12]
   + lFrMC[30,12] >= -2
 lfrMC0a: lfrMCa[30,0] = 0
 lfrMC1a: - ldetWa[30,4] + lfrMCa[30,4] <= 0
 lfrMC1a: - ldw[30,8] + lfrMCa[30,4] <= 0
 lfrMC1a: - ldx[31,8] + lfrMCa[30,4] <= 0
 lfrMC1a: ldw[30,4] + lfrMCa[30,4] <= 1
 lfrMC1a: - ldx[31,8] + ldw[30,4] - ldw[30,8] - ldetWa[30,4] + lfrMCa[30,4]
   >= -2
 lfrMC2a: - ldetWa[30,8] + lfrMCa[30,8] <= 0
 lfrMC2a: - ldx[31,4] + lfrMCa[30,8] <= 0
 lfrMC2a: - ldx[31,12] + lfrMCa[30,8] <= 0
 lfrMC2a: ldw[30,8] + lfrMCa[30,8] <= 1
 lfrMC2a: - ldx[31,4] - ldx[31,12] + ldw[30,8] - ldetWa[30,8]
   + lfrMCa[30,8] >= -2
 lfrMC3a: - ldetWa[30,12] + lfrMCa[30,12] <= 0
 lfrMC3a: - ldx[31,0] + lfrMCa[30,12] <= 0
 lfrMC3a: - ldx[31,12] + lfrMCa[30,12] <= 0
 lfrMC3a: ldw[30,12] + lfrMCa[30,12] <= 1
 lfrMC3a: - ldx[31,0] - ldx[31,12] + ldw[30,12] - ldetWa[30,12]
   + lfrMCa[30,12] >= -2
 lfrMC0: lFrMC[30,1] = 0
 lfrMC1: - lDetW[30,5] + lFrMC[30,5] <= 0
 lfrMC1: - ldw[30,9] + lFrMC[30,5] <= 0
 lfrMC1: - ldx[31,9] + lFrMC[30,5] <= 0
 lfrMC1: ldw[30,5] + lFrMC[30,5] <= 1
 lfrMC1: - ldx[31,9] + ldw[30,5] - ldw[30,9] - lDetW[30,5] + lFrMC[30,5]
   >= -2
 lfrMC2: - lDetW[30,9] + lFrMC[30,9] <= 0
 lfrMC2: - ldx[31,5] + lFrMC[30,9] <= 0
 lfrMC2: - ldx[31,13] + lFrMC[30,9] <= 0
 lfrMC2: ldw[30,9] + lFrMC[30,9] <= 1
 lfrMC2: - ldx[31,5] - ldx[31,13] + ldw[30,9] - lDetW[30,9] + lFrMC[30,9]
   >= -2
 lfrMC3: - lDetW[30,13] + lFrMC[30,13] <= 0
 lfrMC3: - ldx[31,1] + lFrMC[30,13] <= 0
 lfrMC3: - ldx[31,13] + lFrMC[30,13] <= 0
 lfrMC3: ldw[30,13] + lFrMC[30,13] <= 1
 lfrMC3: - ldx[31,1] - ldx[31,13] + ldw[30,13] - lDetW[30,13]
   + lFrMC[30,13] >= -2
 lfrMC0a: lfrMCa[30,1] = 0
 lfrMC1a: - ldetWa[30,5] + lfrMCa[30,5] <= 0
 lfrMC1a: - ldw[30,9] + lfrMCa[30,5] <= 0
 lfrMC1a: - ldx[31,9] + lfrMCa[30,5] <= 0
 lfrMC1a: ldw[30,5] + lfrMCa[30,5] <= 1
 lfrMC1a: - ldx[31,9] + ldw[30,5] - ldw[30,9] - ldetWa[30,5] + lfrMCa[30,5]
   >= -2
 lfrMC2a: - ldetWa[30,9] + lfrMCa[30,9] <= 0
 lfrMC2a: - ldx[31,5] + lfrMCa[30,9] <= 0
 lfrMC2a: - ldx[31,13] + lfrMCa[30,9] <= 0
 lfrMC2a: ldw[30,9] + lfrMCa[30,9] <= 1
 lfrMC2a: - ldx[31,5] - ldx[31,13] + ldw[30,9] - ldetWa[30,9]
   + lfrMCa[30,9] >= -2
 lfrMC3a: - ldetWa[30,13] + lfrMCa[30,13] <= 0
 lfrMC3a: - ldx[31,1] + lfrMCa[30,13] <= 0
 lfrMC3a: - ldx[31,13] + lfrMCa[30,13] <= 0
 lfrMC3a: ldw[30,13] + lfrMCa[30,13] <= 1
 lfrMC3a: - ldx[31,1] - ldx[31,13] + ldw[30,13] - ldetWa[30,13]
   + lfrMCa[30,13] >= -2
 lfrMC0: lFrMC[30,2] = 0
 lfrMC1: - lDetW[30,6] + lFrMC[30,6] <= 0
 lfrMC1: - ldw[30,10] + lFrMC[30,6] <= 0
 lfrMC1: - ldx[31,10] + lFrMC[30,6] <= 0
 lfrMC1: ldw[30,6] + lFrMC[30,6] <= 1
 lfrMC1: - ldx[31,10] + ldw[30,6] - ldw[30,10] - lDetW[30,6] + lFrMC[30,6]
   >= -2
 lfrMC2: - lDetW[30,10] + lFrMC[30,10] <= 0
 lfrMC2: - ldx[31,6] + lFrMC[30,10] <= 0
 lfrMC2: - ldx[31,14] + lFrMC[30,10] <= 0
 lfrMC2: ldw[30,10] + lFrMC[30,10] <= 1
 lfrMC2: - ldx[31,6] - ldx[31,14] + ldw[30,10] - lDetW[30,10]
   + lFrMC[30,10] >= -2
 lfrMC3: - lDetW[30,14] + lFrMC[30,14] <= 0
 lfrMC3: - ldx[31,2] + lFrMC[30,14] <= 0
 lfrMC3: - ldx[31,14] + lFrMC[30,14] <= 0
 lfrMC3: ldw[30,14] + lFrMC[30,14] <= 1
 lfrMC3: - ldx[31,2] - ldx[31,14] + ldw[30,14] - lDetW[30,14]
   + lFrMC[30,14] >= -2
 lfrMC0a: lfrMCa[30,2] = 0
 lfrMC1a: - ldetWa[30,6] + lfrMCa[30,6] <= 0
 lfrMC1a: - ldw[30,10] + lfrMCa[30,6] <= 0
 lfrMC1a: - ldx[31,10] + lfrMCa[30,6] <= 0
 lfrMC1a: ldw[30,6] + lfrMCa[30,6] <= 1
 lfrMC1a: - ldx[31,10] + ldw[30,6] - ldw[30,10] - ldetWa[30,6]
   + lfrMCa[30,6] >= -2
 lfrMC2a: - ldetWa[30,10] + lfrMCa[30,10] <= 0
 lfrMC2a: - ldx[31,6] + lfrMCa[30,10] <= 0
 lfrMC2a: - ldx[31,14] + lfrMCa[30,10] <= 0
 lfrMC2a: ldw[30,10] + lfrMCa[30,10] <= 1
 lfrMC2a: - ldx[31,6] - ldx[31,14] + ldw[30,10] - ldetWa[30,10]
   + lfrMCa[30,10] >= -2
 lfrMC3a: - ldetWa[30,14] + lfrMCa[30,14] <= 0
 lfrMC3a: - ldx[31,2] + lfrMCa[30,14] <= 0
 lfrMC3a: - ldx[31,14] + lfrMCa[30,14] <= 0
 lfrMC3a: ldw[30,14] + lfrMCa[30,14] <= 1
 lfrMC3a: - ldx[31,2] - ldx[31,14] + ldw[30,14] - ldetWa[30,14]
   + lfrMCa[30,14] >= -2
 lfrMC0: lFrMC[30,3] = 0
 lfrMC1: - lDetW[30,7] + lFrMC[30,7] <= 0
 lfrMC1: - ldw[30,11] + lFrMC[30,7] <= 0
 lfrMC1: - ldx[31,11] + lFrMC[30,7] <= 0
 lfrMC1: ldw[30,7] + lFrMC[30,7] <= 1
 lfrMC1: - ldx[31,11] + ldw[30,7] - ldw[30,11] - lDetW[30,7] + lFrMC[30,7]
   >= -2
 lfrMC2: - lDetW[30,11] + lFrMC[30,11] <= 0
 lfrMC2: - ldx[31,7] + lFrMC[30,11] <= 0
 lfrMC2: - ldx[31,15] + lFrMC[30,11] <= 0
 lfrMC2: ldw[30,11] + lFrMC[30,11] <= 1
 lfrMC2: - ldx[31,7] - ldx[31,15] + ldw[30,11] - lDetW[30,11]
   + lFrMC[30,11] >= -2
 lfrMC3: - lDetW[30,15] + lFrMC[30,15] <= 0
 lfrMC3: - ldx[31,3] + lFrMC[30,15] <= 0
 lfrMC3: - ldx[31,15] + lFrMC[30,15] <= 0
 lfrMC3: ldw[30,15] + lFrMC[30,15] <= 1
 lfrMC3: - ldx[31,3] - ldx[31,15] + ldw[30,15] - lDetW[30,15]
   + lFrMC[30,15] >= -2
 lfrMC0a: lfrMCa[30,3] = 0
 lfrMC1a: - ldetWa[30,7] + lfrMCa[30,7] <= 0
 lfrMC1a: - ldw[30,11] + lfrMCa[30,7] <= 0
 lfrMC1a: - ldx[31,11] + lfrMCa[30,7] <= 0
 lfrMC1a: ldw[30,7] + lfrMCa[30,7] <= 1
 lfrMC1a: - ldx[31,11] + ldw[30,7] - ldw[30,11] - ldetWa[30,7]
   + lfrMCa[30,7] >= -2
 lfrMC2a: - ldetWa[30,11] + lfrMCa[30,11] <= 0
 lfrMC2a: - ldx[31,7] + lfrMCa[30,11] <= 0
 lfrMC2a: - ldx[31,15] + lfrMCa[30,11] <= 0
 lfrMC2a: ldw[30,11] + lfrMCa[30,11] <= 1
 lfrMC2a: - ldx[31,7] - ldx[31,15] + ldw[30,11] - ldetWa[30,11]
   + lfrMCa[30,11] >= -2
 lfrMC3a: - ldetWa[30,15] + lfrMCa[30,15] <= 0
 lfrMC3a: - ldx[31,3] + lfrMCa[30,15] <= 0
 lfrMC3a: - ldx[31,15] + lfrMCa[30,15] <= 0
 lfrMC3a: ldw[30,15] + lfrMCa[30,15] <= 1
 lfrMC3a: - ldx[31,3] - ldx[31,15] + ldw[30,15] - ldetWa[30,15]
   + lfrMCa[30,15] >= -2
 lfrMC0: lFrMC[31,0] = 0
 lfrMC1: - lDetW[31,4] + lFrMC[31,4] <= 0
 lfrMC1: - ldw[31,8] + lFrMC[31,4] <= 0
 lfrMC1: - ldx[32,8] + lFrMC[31,4] <= 0
 lfrMC1: ldw[31,4] + lFrMC[31,4] <= 1
 lfrMC1: - ldx[32,8] + ldw[31,4] - ldw[31,8] - lDetW[31,4] + lFrMC[31,4]
   >= -2
 lfrMC2: - lDetW[31,8] + lFrMC[31,8] <= 0
 lfrMC2: - ldx[32,4] + lFrMC[31,8] <= 0
 lfrMC2: - ldx[32,12] + lFrMC[31,8] <= 0
 lfrMC2: ldw[31,8] + lFrMC[31,8] <= 1
 lfrMC2: - ldx[32,4] - ldx[32,12] + ldw[31,8] - lDetW[31,8] + lFrMC[31,8]
   >= -2
 lfrMC3: - lDetW[31,12] + lFrMC[31,12] <= 0
 lfrMC3: - ldx[32,0] + lFrMC[31,12] <= 0
 lfrMC3: - ldx[32,12] + lFrMC[31,12] <= 0
 lfrMC3: ldw[31,12] + lFrMC[31,12] <= 1
 lfrMC3: - ldx[32,0] - ldx[32,12] + ldw[31,12] - lDetW[31,12]
   + lFrMC[31,12] >= -2
 lfrMC0a: lfrMCa[31,0] = 0
 lfrMC1a: - ldetWa[31,4] + lfrMCa[31,4] <= 0
 lfrMC1a: - ldw[31,8] + lfrMCa[31,4] <= 0
 lfrMC1a: - ldx[32,8] + lfrMCa[31,4] <= 0
 lfrMC1a: ldw[31,4] + lfrMCa[31,4] <= 1
 lfrMC1a: - ldx[32,8] + ldw[31,4] - ldw[31,8] - ldetWa[31,4] + lfrMCa[31,4]
   >= -2
 lfrMC2a: - ldetWa[31,8] + lfrMCa[31,8] <= 0
 lfrMC2a: - ldx[32,4] + lfrMCa[31,8] <= 0
 lfrMC2a: - ldx[32,12] + lfrMCa[31,8] <= 0
 lfrMC2a: ldw[31,8] + lfrMCa[31,8] <= 1
 lfrMC2a: - ldx[32,4] - ldx[32,12] + ldw[31,8] - ldetWa[31,8]
   + lfrMCa[31,8] >= -2
 lfrMC3a: - ldetWa[31,12] + lfrMCa[31,12] <= 0
 lfrMC3a: - ldx[32,0] + lfrMCa[31,12] <= 0
 lfrMC3a: - ldx[32,12] + lfrMCa[31,12] <= 0
 lfrMC3a: ldw[31,12] + lfrMCa[31,12] <= 1
 lfrMC3a: - ldx[32,0] - ldx[32,12] + ldw[31,12] - ldetWa[31,12]
   + lfrMCa[31,12] >= -2
 lfrMC0: lFrMC[31,1] = 0
 lfrMC1: - lDetW[31,5] + lFrMC[31,5] <= 0
 lfrMC1: - ldw[31,9] + lFrMC[31,5] <= 0
 lfrMC1: - ldx[32,9] + lFrMC[31,5] <= 0
 lfrMC1: ldw[31,5] + lFrMC[31,5] <= 1
 lfrMC1: - ldx[32,9] + ldw[31,5] - ldw[31,9] - lDetW[31,5] + lFrMC[31,5]
   >= -2
 lfrMC2: - lDetW[31,9] + lFrMC[31,9] <= 0
 lfrMC2: - ldx[32,5] + lFrMC[31,9] <= 0
 lfrMC2: - ldx[32,13] + lFrMC[31,9] <= 0
 lfrMC2: ldw[31,9] + lFrMC[31,9] <= 1
 lfrMC2: - ldx[32,5] - ldx[32,13] + ldw[31,9] - lDetW[31,9] + lFrMC[31,9]
   >= -2
 lfrMC3: - lDetW[31,13] + lFrMC[31,13] <= 0
 lfrMC3: - ldx[32,1] + lFrMC[31,13] <= 0
 lfrMC3: - ldx[32,13] + lFrMC[31,13] <= 0
 lfrMC3: ldw[31,13] + lFrMC[31,13] <= 1
 lfrMC3: - ldx[32,1] - ldx[32,13] + ldw[31,13] - lDetW[31,13]
   + lFrMC[31,13] >= -2
 lfrMC0a: lfrMCa[31,1] = 0
 lfrMC1a: - ldetWa[31,5] + lfrMCa[31,5] <= 0
 lfrMC1a: - ldw[31,9] + lfrMCa[31,5] <= 0
 lfrMC1a: - ldx[32,9] + lfrMCa[31,5] <= 0
 lfrMC1a: ldw[31,5] + lfrMCa[31,5] <= 1
 lfrMC1a: - ldx[32,9] + ldw[31,5] - ldw[31,9] - ldetWa[31,5] + lfrMCa[31,5]
   >= -2
 lfrMC2a: - ldetWa[31,9] + lfrMCa[31,9] <= 0
 lfrMC2a: - ldx[32,5] + lfrMCa[31,9] <= 0
 lfrMC2a: - ldx[32,13] + lfrMCa[31,9] <= 0
 lfrMC2a: ldw[31,9] + lfrMCa[31,9] <= 1
 lfrMC2a: - ldx[32,5] - ldx[32,13] + ldw[31,9] - ldetWa[31,9]
   + lfrMCa[31,9] >= -2
 lfrMC3a: - ldetWa[31,13] + lfrMCa[31,13] <= 0
 lfrMC3a: - ldx[32,1] + lfrMCa[31,13] <= 0
 lfrMC3a: - ldx[32,13] + lfrMCa[31,13] <= 0
 lfrMC3a: ldw[31,13] + lfrMCa[31,13] <= 1
 lfrMC3a: - ldx[32,1] - ldx[32,13] + ldw[31,13] - ldetWa[31,13]
   + lfrMCa[31,13] >= -2
 lfrMC0: lFrMC[31,2] = 0
 lfrMC1: - lDetW[31,6] + lFrMC[31,6] <= 0
 lfrMC1: - ldw[31,10] + lFrMC[31,6] <= 0
 lfrMC1: - ldx[32,10] + lFrMC[31,6] <= 0
 lfrMC1: ldw[31,6] + lFrMC[31,6] <= 1
 lfrMC1: - ldx[32,10] + ldw[31,6] - ldw[31,10] - lDetW[31,6] + lFrMC[31,6]
   >= -2
 lfrMC2: - lDetW[31,10] + lFrMC[31,10] <= 0
 lfrMC2: - ldx[32,6] + lFrMC[31,10] <= 0
 lfrMC2: - ldx[32,14] + lFrMC[31,10] <= 0
 lfrMC2: ldw[31,10] + lFrMC[31,10] <= 1
 lfrMC2: - ldx[32,6] - ldx[32,14] + ldw[31,10] - lDetW[31,10]
   + lFrMC[31,10] >= -2
 lfrMC3: - lDetW[31,14] + lFrMC[31,14] <= 0
 lfrMC3: - ldx[32,2] + lFrMC[31,14] <= 0
 lfrMC3: - ldx[32,14] + lFrMC[31,14] <= 0
 lfrMC3: ldw[31,14] + lFrMC[31,14] <= 1
 lfrMC3: - ldx[32,2] - ldx[32,14] + ldw[31,14] - lDetW[31,14]
   + lFrMC[31,14] >= -2
 lfrMC0a: lfrMCa[31,2] = 0
 lfrMC1a: - ldetWa[31,6] + lfrMCa[31,6] <= 0
 lfrMC1a: - ldw[31,10] + lfrMCa[31,6] <= 0
 lfrMC1a: - ldx[32,10] + lfrMCa[31,6] <= 0
 lfrMC1a: ldw[31,6] + lfrMCa[31,6] <= 1
 lfrMC1a: - ldx[32,10] + ldw[31,6] - ldw[31,10] - ldetWa[31,6]
   + lfrMCa[31,6] >= -2
 lfrMC2a: - ldetWa[31,10] + lfrMCa[31,10] <= 0
 lfrMC2a: - ldx[32,6] + lfrMCa[31,10] <= 0
 lfrMC2a: - ldx[32,14] + lfrMCa[31,10] <= 0
 lfrMC2a: ldw[31,10] + lfrMCa[31,10] <= 1
 lfrMC2a: - ldx[32,6] - ldx[32,14] + ldw[31,10] - ldetWa[31,10]
   + lfrMCa[31,10] >= -2
 lfrMC3a: - ldetWa[31,14] + lfrMCa[31,14] <= 0
 lfrMC3a: - ldx[32,2] + lfrMCa[31,14] <= 0
 lfrMC3a: - ldx[32,14] + lfrMCa[31,14] <= 0
 lfrMC3a: ldw[31,14] + lfrMCa[31,14] <= 1
 lfrMC3a: - ldx[32,2] - ldx[32,14] + ldw[31,14] - ldetWa[31,14]
   + lfrMCa[31,14] >= -2
 lfrMC0: lFrMC[31,3] = 0
 lfrMC1: - lDetW[31,7] + lFrMC[31,7] <= 0
 lfrMC1: - ldw[31,11] + lFrMC[31,7] <= 0
 lfrMC1: - ldx[32,11] + lFrMC[31,7] <= 0
 lfrMC1: ldw[31,7] + lFrMC[31,7] <= 1
 lfrMC1: - ldx[32,11] + ldw[31,7] - ldw[31,11] - lDetW[31,7] + lFrMC[31,7]
   >= -2
 lfrMC2: - lDetW[31,11] + lFrMC[31,11] <= 0
 lfrMC2: - ldx[32,7] + lFrMC[31,11] <= 0
 lfrMC2: - ldx[32,15] + lFrMC[31,11] <= 0
 lfrMC2: ldw[31,11] + lFrMC[31,11] <= 1
 lfrMC2: - ldx[32,7] - ldx[32,15] + ldw[31,11] - lDetW[31,11]
   + lFrMC[31,11] >= -2
 lfrMC3: - lDetW[31,15] + lFrMC[31,15] <= 0
 lfrMC3: - ldx[32,3] + lFrMC[31,15] <= 0
 lfrMC3: - ldx[32,15] + lFrMC[31,15] <= 0
 lfrMC3: ldw[31,15] + lFrMC[31,15] <= 1
 lfrMC3: - ldx[32,3] - ldx[32,15] + ldw[31,15] - lDetW[31,15]
   + lFrMC[31,15] >= -2
 lfrMC0a: lfrMCa[31,3] = 0
 lfrMC1a: - ldetWa[31,7] + lfrMCa[31,7] <= 0
 lfrMC1a: - ldw[31,11] + lfrMCa[31,7] <= 0
 lfrMC1a: - ldx[32,11] + lfrMCa[31,7] <= 0
 lfrMC1a: ldw[31,7] + lfrMCa[31,7] <= 1
 lfrMC1a: - ldx[32,11] + ldw[31,7] - ldw[31,11] - ldetWa[31,7]
   + lfrMCa[31,7] >= -2
 lfrMC2a: - ldetWa[31,11] + lfrMCa[31,11] <= 0
 lfrMC2a: - ldx[32,7] + lfrMCa[31,11] <= 0
 lfrMC2a: - ldx[32,15] + lfrMCa[31,11] <= 0
 lfrMC2a: ldw[31,11] + lfrMCa[31,11] <= 1
 lfrMC2a: - ldx[32,7] - ldx[32,15] + ldw[31,11] - ldetWa[31,11]
   + lfrMCa[31,11] >= -2
 lfrMC3a: - ldetWa[31,15] + lfrMCa[31,15] <= 0
 lfrMC3a: - ldx[32,3] + lfrMCa[31,15] <= 0
 lfrMC3a: - ldx[32,15] + lfrMCa[31,15] <= 0
 lfrMC3a: ldw[31,15] + lfrMCa[31,15] <= 1
 lfrMC3a: - ldx[32,3] - ldx[32,15] + ldw[31,15] - ldetWa[31,15]
   + lfrMCa[31,15] >= -2
 rfLfr: ldz[32,0] + ldz[32,1] + ldz[32,2] + ldz[32,3] + ldz[32,4]
   + ldz[32,5] + ldz[32,6] + ldz[32,7] + ldz[32,8] + ldz[32,9] + ldz[32,10]
   + ldz[32,11] + ldz[32,12] + ldz[32,13] + ldz[32,14] + ldz[32,15]
   - lfrSCa[28,0] - lfrSCa[28,1] - lfrSCa[28,2] - lfrSCa[28,3]
   - lfrSCa[28,4] - lfrSCa[28,5] - lfrSCa[28,6] - lfrSCa[28,7]
   - lfrSCa[28,8] - lfrSCa[28,9] - lfrSCa[28,10] - lfrSCa[28,11]
   - lfrSCa[28,12] - lfrSCa[28,13] - lfrSCa[28,14] - lfrSCa[28,15]
   - lfrSCa[29,0] - lfrSCa[29,1] - lfrSCa[29,2] - lfrSCa[29,3]
   - lfrSCa[29,4] - lfrSCa[29,5] - lfrSCa[29,6] - lfrSCa[29,7]
   - lfrSCa[29,8] - lfrSCa[29,9] - lfrSCa[29,10] - lfrSCa[29,11]
   - lfrSCa[29,12] - lfrSCa[29,13] - lfrSCa[29,14] - lfrSCa[29,15]
   - lfrSCa[30,0] - lfrSCa[30,1] - lfrSCa[30,2] - lfrSCa[30,3]
   - lfrSCa[30,4] - lfrSCa[30,5] - lfrSCa[30,6] - lfrSCa[30,7]
   - lfrSCa[30,8] - lfrSCa[30,9] - lfrSCa[30,10] - lfrSCa[30,11]
   - lfrSCa[30,12] - lfrSCa[30,13] - lfrSCa[30,14] - lfrSCa[30,15]
   - lfrSCa[31,0] - lfrSCa[31,1] - lfrSCa[31,2] - lfrSCa[31,3]
   - lfrSCa[31,4] - lfrSCa[31,5] - lfrSCa[31,6] - lfrSCa[31,7]
   - lfrSCa[31,8] - lfrSCa[31,9] - lfrSCa[31,10] - lfrSCa[31,11]
   - lfrSCa[31,12] - lfrSCa[31,13] - lfrSCa[31,14] - lfrSCa[31,15]
   - lfrSCa[32,0] - lfrSCa[32,1] - lfrSCa[32,2] - lfrSCa[32,3]
   - lfrSCa[32,4] - lfrSCa[32,5] - lfrSCa[32,6] - lfrSCa[32,7]
   - lfrSCa[32,8] - lfrSCa[32,9] - lfrSCa[32,10] - lfrSCa[32,11]
   - lfrSCa[32,12] - lfrSCa[32,13] - lfrSCa[32,14] - lfrSCa[32,15]
   - lfrMCa[28,0] - lfrMCa[28,1] - lfrMCa[28,2] - lfrMCa[28,3]
   - lfrMCa[28,4] - lfrMCa[28,5] - lfrMCa[28,6] - lfrMCa[28,7]
   - lfrMCa[28,8] - lfrMCa[28,9] - lfrMCa[28,10] - lfrMCa[28,11]
   - lfrMCa[28,12] - lfrMCa[28,13] - lfrMCa[28,14] - lfrMCa[28,15]
   - lfrMCa[29,0] - lfrMCa[29,1] - lfrMCa[29,2] - lfrMCa[29,3]
   - lfrMCa[29,4] - lfrMCa[29,5] - lfrMCa[29,6] - lfrMCa[29,7]
   - lfrMCa[29,8] - lfrMCa[29,9] - lfrMCa[29,10] - lfrMCa[29,11]
   - lfrMCa[29,12] - lfrMCa[29,13] - lfrMCa[29,14] - lfrMCa[29,15]
   - lfrMCa[30,0] - lfrMCa[30,1] - lfrMCa[30,2] - lfrMCa[30,3]
   - lfrMCa[30,4] - lfrMCa[30,5] - lfrMCa[30,6] - lfrMCa[30,7]
   - lfrMCa[30,8] - lfrMCa[30,9] - lfrMCa[30,10] - lfrMCa[30,11]
   - lfrMCa[30,12] - lfrMCa[30,13] - lfrMCa[30,14] - lfrMCa[30,15]
   - lfrMCa[31,0] - lfrMCa[31,1] - lfrMCa[31,2] - lfrMCa[31,3]
   - lfrMCa[31,4] - lfrMCa[31,5] - lfrMCa[31,6] - lfrMCa[31,7]
   - lfrMCa[31,8] - lfrMCa[31,9] - lfrMCa[31,10] - lfrMCa[31,11]
   - lfrMCa[31,12] - lfrMCa[31,13] - lfrMCa[31,14] - lfrMCa[31,15] = 0
 KeyBgeI: - ugstka[0,0] - ugstka[1,8] - ugstka[2,2] - ugstka[3,10]
   - ugstka[4,4] - ugstka[5,12] - lgstka[28,7] - lgstka[29,15]
   - lgstka[30,1] - lgstka[31,9] - lgstka[32,0] + iLANE[0] = 0
 KeyBgeI: - ugstka[0,1] - ugstka[1,9] - ugstka[2,0] - ugstka[3,8]
   - ugstka[4,2] - ugstka[5,10] - lgstka[28,3] - lgstka[29,11]
   - lgstka[30,7] - lgstka[31,15] - lgstka[32,1] + iLANE[1] = 0
 KeyBgeI: - ugstka[0,2] - ugstka[1,10] - ugstka[2,4] - ugstka[3,12]
   - ugstka[4,6] - ugstka[5,14] - lgstka[28,1] - lgstka[29,9]
   - lgstka[30,0] - lgstka[31,8] - lgstka[32,2] + iLANE[2] = 0
 KeyBgeI: - ugstka[0,3] - ugstka[1,11] - ugstka[2,7] - ugstka[3,15]
   - ugstka[4,1] - ugstka[5,9] - lgstka[28,6] - lgstka[29,14]
   - lgstka[30,5] - lgstka[31,13] - lgstka[32,3] + iLANE[3] = 0
 KeyBgeI: - ugstka[0,4] - ugstka[1,12] - ugstka[2,6] - ugstka[3,14]
   - ugstka[4,5] - ugstka[5,13] - lgstka[28,0] - lgstka[29,8]
   - lgstka[30,2] - lgstka[31,10] - lgstka[32,4] + iLANE[4] = 0
 KeyBgeI: - ugstka[0,5] - ugstka[1,13] - ugstka[2,3] - ugstka[3,11]
   - ugstka[4,7] - ugstka[5,15] - lgstka[28,4] - lgstka[29,12]
   - lgstka[30,6] - lgstka[31,14] - lgstka[32,5] + iLANE[5] = 0
 KeyBgeI: - ugstka[0,6] - ugstka[1,14] - ugstka[2,5] - ugstka[3,13]
   - ugstka[4,3] - ugstka[5,11] - lgstka[28,2] - lgstka[29,10]
   - lgstka[30,4] - lgstka[31,12] - lgstka[32,6] + iLANE[6] = 0
 KeyBgeI: - ugstka[0,7] - ugstka[1,15] - ugstka[2,1] - ugstka[3,9]
   - ugstka[4,0] - ugstka[5,8] - lgstka[28,5] - lgstka[29,13]
   - lgstka[30,3] - lgstka[31,11] - lgstka[32,7] + iLANE[7] = 0
 KeyBgeI: - ugstka[0,8] - ugstka[1,2] - ugstka[2,10] - ugstka[3,4]
   - ugstka[4,12] - ugstka[5,6] - lgstka[28,15] - lgstka[29,1]
   - lgstka[30,9] - lgstka[31,0] - lgstka[32,8] + iLANE[8] = 0
 KeyBgeI: - ugstka[0,9] - ugstka[1,0] - ugstka[2,8] - ugstka[3,2]
   - ugstka[4,10] - ugstka[5,4] - lgstka[28,11] - lgstka[29,7]
   - lgstka[30,15] - lgstka[31,1] - lgstka[32,9] + iLANE[9] = 0
 KeyBgeI: - ugstka[0,10] - ugstka[1,4] - ugstka[2,12] - ugstka[3,6]
   - ugstka[4,14] - ugstka[5,5] - lgstka[28,9] - lgstka[29,0]
   - lgstka[30,8] - lgstka[31,2] - lgstka[32,10] + iLANE[10] = 0
 KeyBgeI: - ugstka[0,11] - ugstka[1,7] - ugstka[2,15] - ugstka[3,1]
   - ugstka[4,9] - ugstka[5,0] - lgstka[28,14] - lgstka[29,5]
   - lgstka[30,13] - lgstka[31,3] - lgstka[32,11] + iLANE[11] = 0
 KeyBgeI: - ugstka[0,12] - ugstka[1,6] - ugstka[2,14] - ugstka[3,5]
   - ugstka[4,13] - ugstka[5,3] - lgstka[28,8] - lgstka[29,2]
   - lgstka[30,10] - lgstka[31,4] - lgstka[32,12] + iLANE[12] = 0
 KeyBgeI: - ugstka[0,13] - ugstka[1,3] - ugstka[2,11] - ugstka[3,7]
   - ugstka[4,15] - ugstka[5,1] - lgstka[28,12] - lgstka[29,6]
   - lgstka[30,14] - lgstka[31,5] - lgstka[32,13] + iLANE[13] = 0
 KeyBgeI: - ugstka[0,14] - ugstka[1,5] - ugstka[2,13] - ugstka[3,3]
   - ugstka[4,11] - ugstka[5,7] - lgstka[28,10] - lgstka[29,4]
   - lgstka[30,12] - lgstka[31,6] - lgstka[32,14] + iLANE[14] = 0
 KeyBgeI: - ugstka[0,15] - ugstka[1,1] - ugstka[2,9] - ugstka[3,0]
   - ugstka[4,8] - ugstka[5,2] - lgstka[28,13] - lgstka[29,3]
   - lgstka[30,11] - lgstka[31,7] - lgstka[32,15] + iLANE[15] = 0
 KeyBgeI: iLANE[0] - 1000 ikz[0] <= 4
 KeyBgeI: iLANE[0] - 1000 ikz[0] >= -995
 KeyBgeI: - iLANE[0] + ikSUM[0] - 1000 ikz[0] <= 0
 KeyBgeI: - iLANE[0] + ikSUM[0] + 1000 ikz[0] >= 0
 KeyBgeI: ikSUM[0] + 1000 ikz[0] <= 1004
 KeyBgeI: ikSUM[0] - 1000 ikz[0] >= -996
 KeyBgeI: iLANE[1] - 1000 ikz[1] <= 4
 KeyBgeI: iLANE[1] - 1000 ikz[1] >= -995
 KeyBgeI: - iLANE[1] + ikSUM[1] - 1000 ikz[1] <= 0
 KeyBgeI: - iLANE[1] + ikSUM[1] + 1000 ikz[1] >= 0
 KeyBgeI: ikSUM[1] + 1000 ikz[1] <= 1004
 KeyBgeI: ikSUM[1] - 1000 ikz[1] >= -996
 KeyBgeI: iLANE[2] - 1000 ikz[2] <= 4
 KeyBgeI: iLANE[2] - 1000 ikz[2] >= -995
 KeyBgeI: - iLANE[2] + ikSUM[2] - 1000 ikz[2] <= 0
 KeyBgeI: - iLANE[2] + ikSUM[2] + 1000 ikz[2] >= 0
 KeyBgeI: ikSUM[2] + 1000 ikz[2] <= 1004
 KeyBgeI: ikSUM[2] - 1000 ikz[2] >= -996
 KeyBgeI: iLANE[3] - 1000 ikz[3] <= 4
 KeyBgeI: iLANE[3] - 1000 ikz[3] >= -995
 KeyBgeI: - iLANE[3] + ikSUM[3] - 1000 ikz[3] <= 0
 KeyBgeI: - iLANE[3] + ikSUM[3] + 1000 ikz[3] >= 0
 KeyBgeI: ikSUM[3] + 1000 ikz[3] <= 1004
 KeyBgeI: ikSUM[3] - 1000 ikz[3] >= -996
 KeyBgeI: iLANE[4] - 1000 ikz[4] <= 4
 KeyBgeI: iLANE[4] - 1000 ikz[4] >= -995
 KeyBgeI: - iLANE[4] + ikSUM[4] - 1000 ikz[4] <= 0
 KeyBgeI: - iLANE[4] + ikSUM[4] + 1000 ikz[4] >= 0
 KeyBgeI: ikSUM[4] + 1000 ikz[4] <= 1004
 KeyBgeI: ikSUM[4] - 1000 ikz[4] >= -996
 KeyBgeI: iLANE[5] - 1000 ikz[5] <= 4
 KeyBgeI: iLANE[5] - 1000 ikz[5] >= -995
 KeyBgeI: - iLANE[5] + ikSUM[5] - 1000 ikz[5] <= 0
 KeyBgeI: - iLANE[5] + ikSUM[5] + 1000 ikz[5] >= 0
 KeyBgeI: ikSUM[5] + 1000 ikz[5] <= 1004
 KeyBgeI: ikSUM[5] - 1000 ikz[5] >= -996
 KeyBgeI: iLANE[6] - 1000 ikz[6] <= 4
 KeyBgeI: iLANE[6] - 1000 ikz[6] >= -995
 KeyBgeI: - iLANE[6] + ikSUM[6] - 1000 ikz[6] <= 0
 KeyBgeI: - iLANE[6] + ikSUM[6] + 1000 ikz[6] >= 0
 KeyBgeI: ikSUM[6] + 1000 ikz[6] <= 1004
 KeyBgeI: ikSUM[6] - 1000 ikz[6] >= -996
 KeyBgeI: iLANE[7] - 1000 ikz[7] <= 4
 KeyBgeI: iLANE[7] - 1000 ikz[7] >= -995
 KeyBgeI: - iLANE[7] + ikSUM[7] - 1000 ikz[7] <= 0
 KeyBgeI: - iLANE[7] + ikSUM[7] + 1000 ikz[7] >= 0
 KeyBgeI: ikSUM[7] + 1000 ikz[7] <= 1004
 KeyBgeI: ikSUM[7] - 1000 ikz[7] >= -996
 KeyBgeI: iLANE[8] - 1000 ikz[8] <= 4
 KeyBgeI: iLANE[8] - 1000 ikz[8] >= -995
 KeyBgeI: - iLANE[8] + ikSUM[8] - 1000 ikz[8] <= 0
 KeyBgeI: - iLANE[8] + ikSUM[8] + 1000 ikz[8] >= 0
 KeyBgeI: ikSUM[8] + 1000 ikz[8] <= 1004
 KeyBgeI: ikSUM[8] - 1000 ikz[8] >= -996
 KeyBgeI: iLANE[9] - 1000 ikz[9] <= 4
 KeyBgeI: iLANE[9] - 1000 ikz[9] >= -995
 KeyBgeI: - iLANE[9] + ikSUM[9] - 1000 ikz[9] <= 0
 KeyBgeI: - iLANE[9] + ikSUM[9] + 1000 ikz[9] >= 0
 KeyBgeI: ikSUM[9] + 1000 ikz[9] <= 1004
 KeyBgeI: ikSUM[9] - 1000 ikz[9] >= -996
 KeyBgeI: iLANE[10] - 1000 ikz[10] <= 4
 KeyBgeI: iLANE[10] - 1000 ikz[10] >= -995
 KeyBgeI: - iLANE[10] + ikSUM[10] - 1000 ikz[10] <= 0
 KeyBgeI: - iLANE[10] + ikSUM[10] + 1000 ikz[10] >= 0
 KeyBgeI: ikSUM[10] + 1000 ikz[10] <= 1004
 KeyBgeI: ikSUM[10] - 1000 ikz[10] >= -996
 KeyBgeI: iLANE[11] - 1000 ikz[11] <= 4
 KeyBgeI: iLANE[11] - 1000 ikz[11] >= -995
 KeyBgeI: - iLANE[11] + ikSUM[11] - 1000 ikz[11] <= 0
 KeyBgeI: - iLANE[11] + ikSUM[11] + 1000 ikz[11] >= 0
 KeyBgeI: ikSUM[11] + 1000 ikz[11] <= 1004
 KeyBgeI: ikSUM[11] - 1000 ikz[11] >= -996
 KeyBgeI: iLANE[12] - 1000 ikz[12] <= 4
 KeyBgeI: iLANE[12] - 1000 ikz[12] >= -995
 KeyBgeI: - iLANE[12] + ikSUM[12] - 1000 ikz[12] <= 0
 KeyBgeI: - iLANE[12] + ikSUM[12] + 1000 ikz[12] >= 0
 KeyBgeI: ikSUM[12] + 1000 ikz[12] <= 1004
 KeyBgeI: ikSUM[12] - 1000 ikz[12] >= -996
 KeyBgeI: iLANE[13] - 1000 ikz[13] <= 4
 KeyBgeI: iLANE[13] - 1000 ikz[13] >= -995
 KeyBgeI: - iLANE[13] + ikSUM[13] - 1000 ikz[13] <= 0
 KeyBgeI: - iLANE[13] + ikSUM[13] + 1000 ikz[13] >= 0
 KeyBgeI: ikSUM[13] + 1000 ikz[13] <= 1004
 KeyBgeI: ikSUM[13] - 1000 ikz[13] >= -996
 KeyBgeI: iLANE[14] - 1000 ikz[14] <= 4
 KeyBgeI: iLANE[14] - 1000 ikz[14] >= -995
 KeyBgeI: - iLANE[14] + ikSUM[14] - 1000 ikz[14] <= 0
 KeyBgeI: - iLANE[14] + ikSUM[14] + 1000 ikz[14] >= 0
 KeyBgeI: ikSUM[14] + 1000 ikz[14] <= 1004
 KeyBgeI: ikSUM[14] - 1000 ikz[14] >= -996
 KeyBgeI: iLANE[15] - 1000 ikz[15] <= 4
 KeyBgeI: iLANE[15] - 1000 ikz[15] >= -995
 KeyBgeI: - iLANE[15] + ikSUM[15] - 1000 ikz[15] <= 0
 KeyBgeI: - iLANE[15] + ikSUM[15] + 1000 ikz[15] >= 0
 KeyBgeI: ikSUM[15] + 1000 ikz[15] <= 1004
 KeyBgeI: ikSUM[15] - 1000 ikz[15] >= -996
 KeyBgeG: - uGstk[0,0] - uGstk[1,8] - uGstk[2,2] - uGstk[3,10] - uGstk[4,4]
   - uGstk[5,12] - lGstk[28,7] - lGstk[29,15] - lGstk[30,1] - lGstk[31,9]
   - lGstk[32,0] + gLANE[0] = 0
 KeyBgeG: - uGstk[0,1] - uGstk[1,9] - uGstk[2,0] - uGstk[3,8] - uGstk[4,2]
   - uGstk[5,10] - lGstk[28,3] - lGstk[29,11] - lGstk[30,7] - lGstk[31,15]
   - lGstk[32,1] + gLANE[1] = 0
 KeyBgeG: - uGstk[0,2] - uGstk[1,10] - uGstk[2,4] - uGstk[3,12]
   - uGstk[4,6] - uGstk[5,14] - lGstk[28,1] - lGstk[29,9] - lGstk[30,0]
   - lGstk[31,8] - lGstk[32,2] + gLANE[2] = 0
 KeyBgeG: - uGstk[0,3] - uGstk[1,11] - uGstk[2,7] - uGstk[3,15]
   - uGstk[4,1] - uGstk[5,9] - lGstk[28,6] - lGstk[29,14] - lGstk[30,5]
   - lGstk[31,13] - lGstk[32,3] + gLANE[3] = 0
 KeyBgeG: - uGstk[0,4] - uGstk[1,12] - uGstk[2,6] - uGstk[3,14]
   - uGstk[4,5] - uGstk[5,13] - lGstk[28,0] - lGstk[29,8] - lGstk[30,2]
   - lGstk[31,10] - lGstk[32,4] + gLANE[4] = 0
 KeyBgeG: - uGstk[0,5] - uGstk[1,13] - uGstk[2,3] - uGstk[3,11]
   - uGstk[4,7] - uGstk[5,15] - lGstk[28,4] - lGstk[29,12] - lGstk[30,6]
   - lGstk[31,14] - lGstk[32,5] + gLANE[5] = 0
 KeyBgeG: - uGstk[0,6] - uGstk[1,14] - uGstk[2,5] - uGstk[3,13]
   - uGstk[4,3] - uGstk[5,11] - lGstk[28,2] - lGstk[29,10] - lGstk[30,4]
   - lGstk[31,12] - lGstk[32,6] + gLANE[6] = 0
 KeyBgeG: - uGstk[0,7] - uGstk[1,15] - uGstk[2,1] - uGstk[3,9] - uGstk[4,0]
   - uGstk[5,8] - lGstk[28,5] - lGstk[29,13] - lGstk[30,3] - lGstk[31,11]
   - lGstk[32,7] + gLANE[7] = 0
 KeyBgeG: - uGstk[0,8] - uGstk[1,2] - uGstk[2,10] - uGstk[3,4]
   - uGstk[4,12] - uGstk[5,6] - lGstk[28,15] - lGstk[29,1] - lGstk[30,9]
   - lGstk[31,0] - lGstk[32,8] + gLANE[8] = 0
 KeyBgeG: - uGstk[0,9] - uGstk[1,0] - uGstk[2,8] - uGstk[3,2] - uGstk[4,10]
   - uGstk[5,4] - lGstk[28,11] - lGstk[29,7] - lGstk[30,15] - lGstk[31,1]
   - lGstk[32,9] + gLANE[9] = 0
 KeyBgeG: - uGstk[0,10] - uGstk[1,4] - uGstk[2,12] - uGstk[3,6]
   - uGstk[4,14] - uGstk[5,5] - lGstk[28,9] - lGstk[29,0] - lGstk[30,8]
   - lGstk[31,2] - lGstk[32,10] + gLANE[10] = 0
 KeyBgeG: - uGstk[0,11] - uGstk[1,7] - uGstk[2,15] - uGstk[3,1]
   - uGstk[4,9] - uGstk[5,0] - lGstk[28,14] - lGstk[29,5] - lGstk[30,13]
   - lGstk[31,3] - lGstk[32,11] + gLANE[11] = 0
 KeyBgeG: - uGstk[0,12] - uGstk[1,6] - uGstk[2,14] - uGstk[3,5]
   - uGstk[4,13] - uGstk[5,3] - lGstk[28,8] - lGstk[29,2] - lGstk[30,10]
   - lGstk[31,4] - lGstk[32,12] + gLANE[12] = 0
 KeyBgeG: - uGstk[0,13] - uGstk[1,3] - uGstk[2,11] - uGstk[3,7]
   - uGstk[4,15] - uGstk[5,1] - lGstk[28,12] - lGstk[29,6] - lGstk[30,14]
   - lGstk[31,5] - lGstk[32,13] + gLANE[13] = 0
 KeyBgeG: - uGstk[0,14] - uGstk[1,5] - uGstk[2,13] - uGstk[3,3]
   - uGstk[4,11] - uGstk[5,7] - lGstk[28,10] - lGstk[29,4] - lGstk[30,12]
   - lGstk[31,6] - lGstk[32,14] + gLANE[14] = 0
 KeyBgeG: - uGstk[0,15] - uGstk[1,1] - uGstk[2,9] - uGstk[3,0] - uGstk[4,8]
   - uGstk[5,2] - lGstk[28,13] - lGstk[29,3] - lGstk[30,11] - lGstk[31,7]
   - lGstk[32,15] + gLANE[15] = 0
 KeyBgeG: gLANE[0] - 1000 gkz[0] <= 4
 KeyBgeG: gLANE[0] - 1000 gkz[0] >= -995
 KeyBgeG: - gLANE[0] + gkSUM[0] - 1000 gkz[0] <= 0
 KeyBgeG: - gLANE[0] + gkSUM[0] + 1000 gkz[0] >= 0
 KeyBgeG: gkSUM[0] + 1000 gkz[0] <= 1004
 KeyBgeG: gkSUM[0] - 1000 gkz[0] >= -996
 KeyBgeG: gLANE[1] - 1000 gkz[1] <= 4
 KeyBgeG: gLANE[1] - 1000 gkz[1] >= -995
 KeyBgeG: - gLANE[1] + gkSUM[1] - 1000 gkz[1] <= 0
 KeyBgeG: - gLANE[1] + gkSUM[1] + 1000 gkz[1] >= 0
 KeyBgeG: gkSUM[1] + 1000 gkz[1] <= 1004
 KeyBgeG: gkSUM[1] - 1000 gkz[1] >= -996
 KeyBgeG: gLANE[2] - 1000 gkz[2] <= 4
 KeyBgeG: gLANE[2] - 1000 gkz[2] >= -995
 KeyBgeG: - gLANE[2] + gkSUM[2] - 1000 gkz[2] <= 0
 KeyBgeG: - gLANE[2] + gkSUM[2] + 1000 gkz[2] >= 0
 KeyBgeG: gkSUM[2] + 1000 gkz[2] <= 1004
 KeyBgeG: gkSUM[2] - 1000 gkz[2] >= -996
 KeyBgeG: gLANE[3] - 1000 gkz[3] <= 4
 KeyBgeG: gLANE[3] - 1000 gkz[3] >= -995
 KeyBgeG: - gLANE[3] + gkSUM[3] - 1000 gkz[3] <= 0
 KeyBgeG: - gLANE[3] + gkSUM[3] + 1000 gkz[3] >= 0
 KeyBgeG: gkSUM[3] + 1000 gkz[3] <= 1004
 KeyBgeG: gkSUM[3] - 1000 gkz[3] >= -996
 KeyBgeG: gLANE[4] - 1000 gkz[4] <= 4
 KeyBgeG: gLANE[4] - 1000 gkz[4] >= -995
 KeyBgeG: - gLANE[4] + gkSUM[4] - 1000 gkz[4] <= 0
 KeyBgeG: - gLANE[4] + gkSUM[4] + 1000 gkz[4] >= 0
 KeyBgeG: gkSUM[4] + 1000 gkz[4] <= 1004
 KeyBgeG: gkSUM[4] - 1000 gkz[4] >= -996
 KeyBgeG: gLANE[5] - 1000 gkz[5] <= 4
 KeyBgeG: gLANE[5] - 1000 gkz[5] >= -995
 KeyBgeG: - gLANE[5] + gkSUM[5] - 1000 gkz[5] <= 0
 KeyBgeG: - gLANE[5] + gkSUM[5] + 1000 gkz[5] >= 0
 KeyBgeG: gkSUM[5] + 1000 gkz[5] <= 1004
 KeyBgeG: gkSUM[5] - 1000 gkz[5] >= -996
 KeyBgeG: gLANE[6] - 1000 gkz[6] <= 4
 KeyBgeG: gLANE[6] - 1000 gkz[6] >= -995
 KeyBgeG: - gLANE[6] + gkSUM[6] - 1000 gkz[6] <= 0
 KeyBgeG: - gLANE[6] + gkSUM[6] + 1000 gkz[6] >= 0
 KeyBgeG: gkSUM[6] + 1000 gkz[6] <= 1004
 KeyBgeG: gkSUM[6] - 1000 gkz[6] >= -996
 KeyBgeG: gLANE[7] - 1000 gkz[7] <= 4
 KeyBgeG: gLANE[7] - 1000 gkz[7] >= -995
 KeyBgeG: - gLANE[7] + gkSUM[7] - 1000 gkz[7] <= 0
 KeyBgeG: - gLANE[7] + gkSUM[7] + 1000 gkz[7] >= 0
 KeyBgeG: gkSUM[7] + 1000 gkz[7] <= 1004
 KeyBgeG: gkSUM[7] - 1000 gkz[7] >= -996
 KeyBgeG: gLANE[8] - 1000 gkz[8] <= 4
 KeyBgeG: gLANE[8] - 1000 gkz[8] >= -995
 KeyBgeG: - gLANE[8] + gkSUM[8] - 1000 gkz[8] <= 0
 KeyBgeG: - gLANE[8] + gkSUM[8] + 1000 gkz[8] >= 0
 KeyBgeG: gkSUM[8] + 1000 gkz[8] <= 1004
 KeyBgeG: gkSUM[8] - 1000 gkz[8] >= -996
 KeyBgeG: gLANE[9] - 1000 gkz[9] <= 4
 KeyBgeG: gLANE[9] - 1000 gkz[9] >= -995
 KeyBgeG: - gLANE[9] + gkSUM[9] - 1000 gkz[9] <= 0
 KeyBgeG: - gLANE[9] + gkSUM[9] + 1000 gkz[9] >= 0
 KeyBgeG: gkSUM[9] + 1000 gkz[9] <= 1004
 KeyBgeG: gkSUM[9] - 1000 gkz[9] >= -996
 KeyBgeG: gLANE[10] - 1000 gkz[10] <= 4
 KeyBgeG: gLANE[10] - 1000 gkz[10] >= -995
 KeyBgeG: - gLANE[10] + gkSUM[10] - 1000 gkz[10] <= 0
 KeyBgeG: - gLANE[10] + gkSUM[10] + 1000 gkz[10] >= 0
 KeyBgeG: gkSUM[10] + 1000 gkz[10] <= 1004
 KeyBgeG: gkSUM[10] - 1000 gkz[10] >= -996
 KeyBgeG: gLANE[11] - 1000 gkz[11] <= 4
 KeyBgeG: gLANE[11] - 1000 gkz[11] >= -995
 KeyBgeG: - gLANE[11] + gkSUM[11] - 1000 gkz[11] <= 0
 KeyBgeG: - gLANE[11] + gkSUM[11] + 1000 gkz[11] >= 0
 KeyBgeG: gkSUM[11] + 1000 gkz[11] <= 1004
 KeyBgeG: gkSUM[11] - 1000 gkz[11] >= -996
 KeyBgeG: gLANE[12] - 1000 gkz[12] <= 4
 KeyBgeG: gLANE[12] - 1000 gkz[12] >= -995
 KeyBgeG: - gLANE[12] + gkSUM[12] - 1000 gkz[12] <= 0
 KeyBgeG: - gLANE[12] + gkSUM[12] + 1000 gkz[12] >= 0
 KeyBgeG: gkSUM[12] + 1000 gkz[12] <= 1004
 KeyBgeG: gkSUM[12] - 1000 gkz[12] >= -996
 KeyBgeG: gLANE[13] - 1000 gkz[13] <= 4
 KeyBgeG: gLANE[13] - 1000 gkz[13] >= -995
 KeyBgeG: - gLANE[13] + gkSUM[13] - 1000 gkz[13] <= 0
 KeyBgeG: - gLANE[13] + gkSUM[13] + 1000 gkz[13] >= 0
 KeyBgeG: gkSUM[13] + 1000 gkz[13] <= 1004
 KeyBgeG: gkSUM[13] - 1000 gkz[13] >= -996
 KeyBgeG: gLANE[14] - 1000 gkz[14] <= 4
 KeyBgeG: gLANE[14] - 1000 gkz[14] >= -995
 KeyBgeG: - gLANE[14] + gkSUM[14] - 1000 gkz[14] <= 0
 KeyBgeG: - gLANE[14] + gkSUM[14] + 1000 gkz[14] >= 0
 KeyBgeG: gkSUM[14] + 1000 gkz[14] <= 1004
 KeyBgeG: gkSUM[14] - 1000 gkz[14] >= -996
 KeyBgeG: gLANE[15] - 1000 gkz[15] <= 4
 KeyBgeG: gLANE[15] - 1000 gkz[15] >= -995
 KeyBgeG: - gLANE[15] + gkSUM[15] - 1000 gkz[15] <= 0
 KeyBgeG: - gLANE[15] + gkSUM[15] + 1000 gkz[15] >= 0
 KeyBgeG: gkSUM[15] + 1000 gkz[15] <= 1004
 KeyBgeG: gkSUM[15] - 1000 gkz[15] >= -996
 IKS: - 4 ikSUM[0] - 4 ikSUM[1] - 4 ikSUM[2] - 4 ikSUM[3] - 4 ikSUM[4]
   - 4 ikSUM[5] - 4 ikSUM[6] - 4 ikSUM[7] - 4 ikSUM[8] - 4 ikSUM[9]
   - 4 ikSUM[10] - 4 ikSUM[11] - 4 ikSUM[12] - 4 ikSUM[13] - 4 ikSUM[14]
   - 4 ikSUM[15] + iks = 0
 GKS: - 4 gkSUM[0] - 4 gkSUM[1] - 4 gkSUM[2] - 4 gkSUM[3] - 4 gkSUM[4]
   - 4 gkSUM[5] - 4 gkSUM[6] - 4 gkSUM[7] - 4 gkSUM[8] - 4 gkSUM[9]
   - 4 gkSUM[10] - 4 gkSUM[11] - 4 gkSUM[12] - 4 gkSUM[13] - 4 gkSUM[14]
   - 4 gkSUM[15] + gks = 0
 rb: - 4 udw[0,0] - 4 udw[0,1] - 4 udw[0,2] - 4 udw[0,3] - 4 udw[0,4]
   - 4 udw[0,5] - 4 udw[0,6] - 4 udw[0,7] - 4 udw[0,8] - 4 udw[0,9]
   - 4 udw[0,10] - 4 udw[0,11] - 4 udw[0,12] - 4 udw[0,13] - 4 udw[0,14]
   - 4 udw[0,15] + rb = 0
 rf: - 4 ldz[32,0] - 4 ldz[32,1] - 4 ldz[32,2] - 4 ldz[32,3] - 4 ldz[32,4]
   - 4 ldz[32,5] - 4 ldz[32,6] - 4 ldz[32,7] - 4 ldz[32,8] - 4 ldz[32,9]
   - 4 ldz[32,10] - 4 ldz[32,11] - 4 ldz[32,12] - 4 ldz[32,13]
   - 4 ldz[32,14] - 4 ldz[32,15] + rf = 0
 cb: 4 udy[6,0] + 4 udy[6,1] + 4 udy[6,2] + 4 udy[6,3] + 4 udy[6,4]
   + 4 udy[6,5] + 4 udy[6,6] + 4 udy[6,7] + 4 udy[6,8] + 4 udy[6,9]
   + 4 udy[6,10] + 4 udy[6,11] + 4 udy[6,12] + 4 udy[6,13] + 4 udy[6,14]
   + 4 udy[6,15] - rb + cb = 0
 cf: 4 ldx[28,0] + 4 ldx[28,1] + 4 ldx[28,2] + 4 ldx[28,3] + 4 ldx[28,4]
   + 4 ldx[28,5] + 4 ldx[28,6] + 4 ldx[28,7] + 4 ldx[28,8] + 4 ldx[28,9]
   + 4 ldx[28,10] + 4 ldx[28,11] + 4 ldx[28,12] + 4 ldx[28,13]
   + 4 ldx[28,14] + 4 ldx[28,15] - rf + cf = 0
 cbp: - 4 uFrSC[1,0] - 4 uFrSC[1,1] - 4 uFrSC[1,2] - 4 uFrSC[1,3]
   - 4 uFrSC[1,4] - 4 uFrSC[1,5] - 4 uFrSC[1,6] - 4 uFrSC[1,7]
   - 4 uFrSC[1,8] - 4 uFrSC[1,9] - 4 uFrSC[1,10] - 4 uFrSC[1,11]
   - 4 uFrSC[1,12] - 4 uFrSC[1,13] - 4 uFrSC[1,14] - 4 uFrSC[1,15]
   - 4 uFrSC[2,0] - 4 uFrSC[2,1] - 4 uFrSC[2,2] - 4 uFrSC[2,3]
   - 4 uFrSC[2,4] - 4 uFrSC[2,5] - 4 uFrSC[2,6] - 4 uFrSC[2,7]
   - 4 uFrSC[2,8] - 4 uFrSC[2,9] - 4 uFrSC[2,10] - 4 uFrSC[2,11]
   - 4 uFrSC[2,12] - 4 uFrSC[2,13] - 4 uFrSC[2,14] - 4 uFrSC[2,15]
   - 4 uFrSC[3,0] - 4 uFrSC[3,1] - 4 uFrSC[3,2] - 4 uFrSC[3,3]
   - 4 uFrSC[3,4] - 4 uFrSC[3,5] - 4 uFrSC[3,6] - 4 uFrSC[3,7]
   - 4 uFrSC[3,8] - 4 uFrSC[3,9] - 4 uFrSC[3,10] - 4 uFrSC[3,11]
   - 4 uFrSC[3,12] - 4 uFrSC[3,13] - 4 uFrSC[3,14] - 4 uFrSC[3,15]
   - 4 uFrSC[4,0] - 4 uFrSC[4,1] - 4 uFrSC[4,2] - 4 uFrSC[4,3]
   - 4 uFrSC[4,4] - 4 uFrSC[4,5] - 4 uFrSC[4,6] - 4 uFrSC[4,7]
   - 4 uFrSC[4,8] - 4 uFrSC[4,9] - 4 uFrSC[4,10] - 4 uFrSC[4,11]
   - 4 uFrSC[4,12] - 4 uFrSC[4,13] - 4 uFrSC[4,14] - 4 uFrSC[4,15]
   - 4 uFrSC[5,0] - 4 uFrSC[5,1] - 4 uFrSC[5,2] - 4 uFrSC[5,3]
   - 4 uFrSC[5,4] - 4 uFrSC[5,5] - 4 uFrSC[5,6] - 4 uFrSC[5,7]
   - 4 uFrSC[5,8] - 4 uFrSC[5,9] - 4 uFrSC[5,10] - 4 uFrSC[5,11]
   - 4 uFrSC[5,12] - 4 uFrSC[5,13] - 4 uFrSC[5,14] - 4 uFrSC[5,15]
   - 4 uFrSC[6,0] - 4 uFrSC[6,1] - 4 uFrSC[6,2] - 4 uFrSC[6,3]
   - 4 uFrSC[6,4] - 4 uFrSC[6,5] - 4 uFrSC[6,6] - 4 uFrSC[6,7]
   - 4 uFrSC[6,8] - 4 uFrSC[6,9] - 4 uFrSC[6,10] - 4 uFrSC[6,11]
   - 4 uFrSC[6,12] - 4 uFrSC[6,13] - 4 uFrSC[6,14] - 4 uFrSC[6,15]
   - 4 uFrMC[2,0] - 4 uFrMC[2,1] - 4 uFrMC[2,2] - 4 uFrMC[2,3]
   - 4 uFrMC[2,4] - 4 uFrMC[2,5] - 4 uFrMC[2,6] - 4 uFrMC[2,7]
   - 4 uFrMC[2,8] - 4 uFrMC[2,9] - 4 uFrMC[2,10] - 4 uFrMC[2,11]
   - 4 uFrMC[2,12] - 4 uFrMC[2,13] - 4 uFrMC[2,14] - 4 uFrMC[2,15]
   - 4 uFrMC[3,0] - 4 uFrMC[3,1] - 4 uFrMC[3,2] - 4 uFrMC[3,3]
   - 4 uFrMC[3,4] - 4 uFrMC[3,5] - 4 uFrMC[3,6] - 4 uFrMC[3,7]
   - 4 uFrMC[3,8] - 4 uFrMC[3,9] - 4 uFrMC[3,10] - 4 uFrMC[3,11]
   - 4 uFrMC[3,12] - 4 uFrMC[3,13] - 4 uFrMC[3,14] - 4 uFrMC[3,15]
   - 4 uFrMC[4,0] - 4 uFrMC[4,1] - 4 uFrMC[4,2] - 4 uFrMC[4,3]
   - 4 uFrMC[4,4] - 4 uFrMC[4,5] - 4 uFrMC[4,6] - 4 uFrMC[4,7]
   - 4 uFrMC[4,8] - 4 uFrMC[4,9] - 4 uFrMC[4,10] - 4 uFrMC[4,11]
   - 4 uFrMC[4,12] - 4 uFrMC[4,13] - 4 uFrMC[4,14] - 4 uFrMC[4,15]
   - 4 uFrMC[5,0] - 4 uFrMC[5,1] - 4 uFrMC[5,2] - 4 uFrMC[5,3]
   - 4 uFrMC[5,4] - 4 uFrMC[5,5] - 4 uFrMC[5,6] - 4 uFrMC[5,7]
   - 4 uFrMC[5,8] - 4 uFrMC[5,9] - 4 uFrMC[5,10] - 4 uFrMC[5,11]
   - 4 uFrMC[5,12] - 4 uFrMC[5,13] - 4 uFrMC[5,14] - 4 uFrMC[5,15]
   - 4 uFrMC[6,0] - 4 uFrMC[6,1] - 4 uFrMC[6,2] - 4 uFrMC[6,3]
   - 4 uFrMC[6,4] - 4 uFrMC[6,5] - 4 uFrMC[6,6] - 4 uFrMC[6,7]
   - 4 uFrMC[6,8] - 4 uFrMC[6,9] - 4 uFrMC[6,10] - 4 uFrMC[6,11]
   - 4 uFrMC[6,12] - 4 uFrMC[6,13] - 4 uFrMC[6,14] - 4 uFrMC[6,15] + cbp
   = 0
 cfp: - 4 lFrSC[28,0] - 4 lFrSC[28,1] - 4 lFrSC[28,2] - 4 lFrSC[28,3]
   - 4 lFrSC[28,4] - 4 lFrSC[28,5] - 4 lFrSC[28,6] - 4 lFrSC[28,7]
   - 4 lFrSC[28,8] - 4 lFrSC[28,9] - 4 lFrSC[28,10] - 4 lFrSC[28,11]
   - 4 lFrSC[28,12] - 4 lFrSC[28,13] - 4 lFrSC[28,14] - 4 lFrSC[28,15]
   - 4 lFrSC[29,0] - 4 lFrSC[29,1] - 4 lFrSC[29,2] - 4 lFrSC[29,3]
   - 4 lFrSC[29,4] - 4 lFrSC[29,5] - 4 lFrSC[29,6] - 4 lFrSC[29,7]
   - 4 lFrSC[29,8] - 4 lFrSC[29,9] - 4 lFrSC[29,10] - 4 lFrSC[29,11]
   - 4 lFrSC[29,12] - 4 lFrSC[29,13] - 4 lFrSC[29,14] - 4 lFrSC[29,15]
   - 4 lFrSC[30,0] - 4 lFrSC[30,1] - 4 lFrSC[30,2] - 4 lFrSC[30,3]
   - 4 lFrSC[30,4] - 4 lFrSC[30,5] - 4 lFrSC[30,6] - 4 lFrSC[30,7]
   - 4 lFrSC[30,8] - 4 lFrSC[30,9] - 4 lFrSC[30,10] - 4 lFrSC[30,11]
   - 4 lFrSC[30,12] - 4 lFrSC[30,13] - 4 lFrSC[30,14] - 4 lFrSC[30,15]
   - 4 lFrSC[31,0] - 4 lFrSC[31,1] - 4 lFrSC[31,2] - 4 lFrSC[31,3]
   - 4 lFrSC[31,4] - 4 lFrSC[31,5] - 4 lFrSC[31,6] - 4 lFrSC[31,7]
   - 4 lFrSC[31,8] - 4 lFrSC[31,9] - 4 lFrSC[31,10] - 4 lFrSC[31,11]
   - 4 lFrSC[31,12] - 4 lFrSC[31,13] - 4 lFrSC[31,14] - 4 lFrSC[31,15]
   - 4 lFrSC[32,0] - 4 lFrSC[32,1] - 4 lFrSC[32,2] - 4 lFrSC[32,3]
   - 4 lFrSC[32,4] - 4 lFrSC[32,5] - 4 lFrSC[32,6] - 4 lFrSC[32,7]
   - 4 lFrSC[32,8] - 4 lFrSC[32,9] - 4 lFrSC[32,10] - 4 lFrSC[32,11]
   - 4 lFrSC[32,12] - 4 lFrSC[32,13] - 4 lFrSC[32,14] - 4 lFrSC[32,15]
   - 4 lFrMC[28,0] - 4 lFrMC[28,1] - 4 lFrMC[28,2] - 4 lFrMC[28,3]
   - 4 lFrMC[28,4] - 4 lFrMC[28,5] - 4 lFrMC[28,6] - 4 lFrMC[28,7]
   - 4 lFrMC[28,8] - 4 lFrMC[28,9] - 4 lFrMC[28,10] - 4 lFrMC[28,11]
   - 4 lFrMC[28,12] - 4 lFrMC[28,13] - 4 lFrMC[28,14] - 4 lFrMC[28,15]
   - 4 lFrMC[29,0] - 4 lFrMC[29,1] - 4 lFrMC[29,2] - 4 lFrMC[29,3]
   - 4 lFrMC[29,4] - 4 lFrMC[29,5] - 4 lFrMC[29,6] - 4 lFrMC[29,7]
   - 4 lFrMC[29,8] - 4 lFrMC[29,9] - 4 lFrMC[29,10] - 4 lFrMC[29,11]
   - 4 lFrMC[29,12] - 4 lFrMC[29,13] - 4 lFrMC[29,14] - 4 lFrMC[29,15]
   - 4 lFrMC[30,0] - 4 lFrMC[30,1] - 4 lFrMC[30,2] - 4 lFrMC[30,3]
   - 4 lFrMC[30,4] - 4 lFrMC[30,5] - 4 lFrMC[30,6] - 4 lFrMC[30,7]
   - 4 lFrMC[30,8] - 4 lFrMC[30,9] - 4 lFrMC[30,10] - 4 lFrMC[30,11]
   - 4 lFrMC[30,12] - 4 lFrMC[30,13] - 4 lFrMC[30,14] - 4 lFrMC[30,15]
   - 4 lFrMC[31,0] - 4 lFrMC[31,1] - 4 lFrMC[31,2] - 4 lFrMC[31,3]
   - 4 lFrMC[31,4] - 4 lFrMC[31,5] - 4 lFrMC[31,6] - 4 lFrMC[31,7]
   - 4 lFrMC[31,8] - 4 lFrMC[31,9] - 4 lFrMC[31,10] - 4 lFrMC[31,11]
   - 4 lFrMC[31,12] - 4 lFrMC[31,13] - 4 lFrMC[31,14] - 4 lFrMC[31,15]
   + cfp = 0
 D: D = 66.25
 Dc: - D + Dc = 2
 Qc: - 2 cb - 2 cf + Qc = 4.5
 Mc0: - Dc + Mc >= 2
 Mc1: gks + Mc - T2 >= 0
 Mc2: 2 cbp + 2 cfp - Qc + Mc >= 0
 Mc3: - iks + Mc >= 0
 T0: - Dc + T0 = 0
 T1: - gks - Dc + T1 = 2
 T2c: - rf - gks + cfp - 2 Dc + T2 = -63
 T31: - gks + 2 cbp + 2 cfp - Qc + T31 = 0
 T32: - iks + T32 = 4.5
 T4: T4 = 224
 Tc: - T0 + Tc >= 0
 Tc: - T1 + Tc >= 0
 Tc: - T2 + Tc >= 0
 Tc: - T31 + Tc >= 0
 Tc: - T32 + Tc >= 0
 Tc: - T4 + Tc >= 0
 mb: - ugstka[0,0] - ugstka[0,1] - ugstka[0,2] - ugstka[0,3] - ugstka[0,4]
   - ugstka[0,5] - ugstka[0,6] - ugstka[0,7] - ugstka[1,0] - ugstka[1,1]
   - ugstka[1,2] - ugstka[1,3] - ugstka[1,4] - ugstka[1,5] - ugstka[1,6]
   - ugstka[1,7] - ugstka[2,0] - ugstka[2,1] - ugstka[2,2] - ugstka[2,3]
   - ugstka[2,4] - ugstka[2,5] - ugstka[2,6] - ugstka[2,7] - ugstka[3,0]
   - ugstka[3,1] - ugstka[3,2] - ugstka[3,3] - ugstka[3,4] - ugstka[3,5]
   - ugstka[3,6] - ugstka[3,7] - ugstka[4,0] - ugstka[4,1] - ugstka[4,2]
   - ugstka[4,3] - ugstka[4,4] - ugstka[4,5] - ugstka[4,6] - ugstka[4,7]
   - ugstka[5,0] - ugstka[5,1] - ugstka[5,2] - ugstka[5,3] - ugstka[5,4]
   - ugstka[5,5] - ugstka[5,6] - ugstka[5,7] + mbi = 0
 mf: - lgstka[28,0] - lgstka[28,1] - lgstka[28,2] - lgstka[28,3]
   - lgstka[28,4] - lgstka[28,5] - lgstka[28,6] - lgstka[28,7]
   - lgstka[29,0] - lgstka[29,1] - lgstka[29,2] - lgstka[29,3]
   - lgstka[29,4] - lgstka[29,5] - lgstka[29,6] - lgstka[29,7]
   - lgstka[30,0] - lgstka[30,1] - lgstka[30,2] - lgstka[30,3]
   - lgstka[30,4] - lgstka[30,5] - lgstka[30,6] - lgstka[30,7]
   - lgstka[31,0] - lgstka[31,1] - lgstka[31,2] - lgstka[31,3]
   - lgstka[31,4] - lgstka[31,5] - lgstka[31,6] - lgstka[31,7]
   - lgstka[32,0] - lgstka[32,1] - lgstka[32,2] - lgstka[32,3]
   - lgstka[32,4] - lgstka[32,5] - lgstka[32,6] - lgstka[32,7] + mfi = 0
 iBCT: [ uDX[17,0] * lDY[17,0] - uDX[17,0] * ldy[17,0]
   + uDX[17,1] * lDY[17,1] - uDX[17,1] * ldy[17,1] + uDX[17,2] * lDY[17,2]
   - uDX[17,2] * ldy[17,2] + uDX[17,3] * lDY[17,3] - uDX[17,3] * ldy[17,3]
   + uDX[17,4] * lDY[17,4] - uDX[17,4] * ldy[17,4] + uDX[17,5] * lDY[17,5]
   - uDX[17,5] * ldy[17,5] + uDX[17,6] * lDY[17,6] - uDX[17,6] * ldy[17,6]
   + uDX[17,7] * lDY[17,7] - uDX[17,7] * ldy[17,7] + uDX[17,8] * lDY[17,8]
   - uDX[17,8] * ldy[17,8] + uDX[17,9] * lDY[17,9] - uDX[17,9] * ldy[17,9]
   + uDX[17,10] * lDY[17,10] - uDX[17,10] * ldy[17,10]
   + uDX[17,11] * lDY[17,11] - uDX[17,11] * ldy[17,11]
   + uDX[17,12] * lDY[17,12] - uDX[17,12] * ldy[17,12]
   + uDX[17,13] * lDY[17,13] - uDX[17,13] * ldy[17,13]
   + uDX[17,14] * lDY[17,14] - uDX[17,14] * ldy[17,14]
   + uDX[17,15] * lDY[17,15] - uDX[17,15] * ldy[17,15]
   - udx[17,0] * lDY[17,0] + udx[17,0] * ldy[17,0] - udx[17,1] * lDY[17,1]
   + udx[17,1] * ldy[17,1] - udx[17,2] * lDY[17,2] + udx[17,2] * ldy[17,2]
   - udx[17,3] * lDY[17,3] + udx[17,3] * ldy[17,3] - udx[17,4] * lDY[17,4]
   + udx[17,4] * ldy[17,4] - udx[17,5] * lDY[17,5] + udx[17,5] * ldy[17,5]
   - udx[17,6] * lDY[17,6] + udx[17,6] * ldy[17,6] - udx[17,7] * lDY[17,7]
   + udx[17,7] * ldy[17,7] - udx[17,8] * lDY[17,8] + udx[17,8] * ldy[17,8]
   - udx[17,9] * lDY[17,9] + udx[17,9] * ldy[17,9]
   - udx[17,10] * lDY[17,10] + udx[17,10] * ldy[17,10]
   - udx[17,11] * lDY[17,11] + udx[17,11] * ldy[17,11]
   - udx[17,12] * lDY[17,12] + udx[17,12] * ldy[17,12]
   - udx[17,13] * lDY[17,13] + udx[17,13] * ldy[17,13]
   - udx[17,14] * lDY[17,14] + udx[17,14] * ldy[17,14]
   - udx[17,15] * lDY[17,15] + udx[17,15] * ldy[17,15] ] >= 1
Bounds
 ikSUM[0] <= 4
 ikSUM[1] <= 4
 ikSUM[2] <= 4
 ikSUM[3] <= 4
 ikSUM[4] <= 4
 ikSUM[5] <= 4
 ikSUM[6] <= 4
 ikSUM[7] <= 4
 ikSUM[8] <= 4
 ikSUM[9] <= 4
 ikSUM[10] <= 4
 ikSUM[11] <= 4
 ikSUM[12] <= 4
 ikSUM[13] <= 4
 ikSUM[14] <= 4
 ikSUM[15] <= 4
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
 uDX[12,13] uDX[12,14] uDX[12,15] uDX[13,0] uDX[13,1] uDX[13,2] uDX[13,3]
 uDX[13,4] uDX[13,5] uDX[13,6] uDX[13,7] uDX[13,8] uDX[13,9] uDX[13,10]
 uDX[13,11] uDX[13,12] uDX[13,13] uDX[13,14] uDX[13,15] uDX[14,0] uDX[14,1]
 uDX[14,2] uDX[14,3] uDX[14,4] uDX[14,5] uDX[14,6] uDX[14,7] uDX[14,8]
 uDX[14,9] uDX[14,10] uDX[14,11] uDX[14,12] uDX[14,13] uDX[14,14]
 uDX[14,15] uDX[15,0] uDX[15,1] uDX[15,2] uDX[15,3] uDX[15,4] uDX[15,5]
 uDX[15,6] uDX[15,7] uDX[15,8] uDX[15,9] uDX[15,10] uDX[15,11] uDX[15,12]
 uDX[15,13] uDX[15,14] uDX[15,15] uDX[16,0] uDX[16,1] uDX[16,2] uDX[16,3]
 uDX[16,4] uDX[16,5] uDX[16,6] uDX[16,7] uDX[16,8] uDX[16,9] uDX[16,10]
 uDX[16,11] uDX[16,12] uDX[16,13] uDX[16,14] uDX[16,15] uDX[17,0] uDX[17,1]
 uDX[17,2] uDX[17,3] uDX[17,4] uDX[17,5] uDX[17,6] uDX[17,7] uDX[17,8]
 uDX[17,9] uDX[17,10] uDX[17,11] uDX[17,12] uDX[17,13] uDX[17,14]
 uDX[17,15] udx[1,0] udx[1,1] udx[1,2] udx[1,3] udx[1,4] udx[1,5] udx[1,6]
 udx[1,7] udx[1,8] udx[1,9] udx[1,10] udx[1,11] udx[1,12] udx[1,13]
 udx[1,14] udx[1,15] udx[2,0] udx[2,1] udx[2,2] udx[2,3] udx[2,4] udx[2,5]
 udx[2,6] udx[2,7] udx[2,8] udx[2,9] udx[2,10] udx[2,11] udx[2,12]
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
 udx[12,13] udx[12,14] udx[12,15] udx[13,0] udx[13,1] udx[13,2] udx[13,3]
 udx[13,4] udx[13,5] udx[13,6] udx[13,7] udx[13,8] udx[13,9] udx[13,10]
 udx[13,11] udx[13,12] udx[13,13] udx[13,14] udx[13,15] udx[14,0] udx[14,1]
 udx[14,2] udx[14,3] udx[14,4] udx[14,5] udx[14,6] udx[14,7] udx[14,8]
 udx[14,9] udx[14,10] udx[14,11] udx[14,12] udx[14,13] udx[14,14]
 udx[14,15] udx[15,0] udx[15,1] udx[15,2] udx[15,3] udx[15,4] udx[15,5]
 udx[15,6] udx[15,7] udx[15,8] udx[15,9] udx[15,10] udx[15,11] udx[15,12]
 udx[15,13] udx[15,14] udx[15,15] udx[16,0] udx[16,1] udx[16,2] udx[16,3]
 udx[16,4] udx[16,5] udx[16,6] udx[16,7] udx[16,8] udx[16,9] udx[16,10]
 udx[16,11] udx[16,12] udx[16,13] udx[16,14] udx[16,15] udx[17,0] udx[17,1]
 udx[17,2] udx[17,3] udx[17,4] udx[17,5] udx[17,6] udx[17,7] udx[17,8]
 udx[17,9] udx[17,10] udx[17,11] udx[17,12] udx[17,13] udx[17,14]
 udx[17,15] uDY[1,0] uDY[1,1] uDY[1,2] uDY[1,3] uDY[1,4] uDY[1,5] uDY[1,6]
 uDY[1,7] uDY[1,8] uDY[1,9] uDY[1,10] uDY[1,11] uDY[1,12] uDY[1,13]
 uDY[1,14] uDY[1,15] uDY[2,0] uDY[2,1] uDY[2,2] uDY[2,3] uDY[2,4] uDY[2,5]
 uDY[2,6] uDY[2,7] uDY[2,8] uDY[2,9] uDY[2,10] uDY[2,11] uDY[2,12]
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
 uDY[7,13] uDY[7,14] uDY[7,15] uDY[8,0] uDY[8,1] uDY[8,2] uDY[8,3] uDY[8,4]
 uDY[8,5] uDY[8,6] uDY[8,7] uDY[8,8] uDY[8,9] uDY[8,10] uDY[8,11] uDY[8,12]
 uDY[8,13] uDY[8,14] uDY[8,15] uDY[9,0] uDY[9,1] uDY[9,2] uDY[9,3] uDY[9,4]
 uDY[9,5] uDY[9,6] uDY[9,7] uDY[9,8] uDY[9,9] uDY[9,10] uDY[9,11] uDY[9,12]
 uDY[9,13] uDY[9,14] uDY[9,15] uDY[10,0] uDY[10,1] uDY[10,2] uDY[10,3]
 uDY[10,4] uDY[10,5] uDY[10,6] uDY[10,7] uDY[10,8] uDY[10,9] uDY[10,10]
 uDY[10,11] uDY[10,12] uDY[10,13] uDY[10,14] uDY[10,15] uDY[11,0] uDY[11,1]
 uDY[11,2] uDY[11,3] uDY[11,4] uDY[11,5] uDY[11,6] uDY[11,7] uDY[11,8]
 uDY[11,9] uDY[11,10] uDY[11,11] uDY[11,12] uDY[11,13] uDY[11,14]
 uDY[11,15] uDY[12,0] uDY[12,1] uDY[12,2] uDY[12,3] uDY[12,4] uDY[12,5]
 uDY[12,6] uDY[12,7] uDY[12,8] uDY[12,9] uDY[12,10] uDY[12,11] uDY[12,12]
 uDY[12,13] uDY[12,14] uDY[12,15] uDY[13,0] uDY[13,1] uDY[13,2] uDY[13,3]
 uDY[13,4] uDY[13,5] uDY[13,6] uDY[13,7] uDY[13,8] uDY[13,9] uDY[13,10]
 uDY[13,11] uDY[13,12] uDY[13,13] uDY[13,14] uDY[13,15] uDY[14,0] uDY[14,1]
 uDY[14,2] uDY[14,3] uDY[14,4] uDY[14,5] uDY[14,6] uDY[14,7] uDY[14,8]
 uDY[14,9] uDY[14,10] uDY[14,11] uDY[14,12] uDY[14,13] uDY[14,14]
 uDY[14,15] uDY[15,0] uDY[15,1] uDY[15,2] uDY[15,3] uDY[15,4] uDY[15,5]
 uDY[15,6] uDY[15,7] uDY[15,8] uDY[15,9] uDY[15,10] uDY[15,11] uDY[15,12]
 uDY[15,13] uDY[15,14] uDY[15,15] uDY[16,0] uDY[16,1] uDY[16,2] uDY[16,3]
 uDY[16,4] uDY[16,5] uDY[16,6] uDY[16,7] uDY[16,8] uDY[16,9] uDY[16,10]
 uDY[16,11] uDY[16,12] uDY[16,13] uDY[16,14] uDY[16,15] udy[1,0] udy[1,1]
 udy[1,2] udy[1,3] udy[1,4] udy[1,5] udy[1,6] udy[1,7] udy[1,8] udy[1,9]
 udy[1,10] udy[1,11] udy[1,12] udy[1,13] udy[1,14] udy[1,15] udy[2,0]
 udy[2,1] udy[2,2] udy[2,3] udy[2,4] udy[2,5] udy[2,6] udy[2,7] udy[2,8]
 udy[2,9] udy[2,10] udy[2,11] udy[2,12] udy[2,13] udy[2,14] udy[2,15]
 udy[3,0] udy[3,1] udy[3,2] udy[3,3] udy[3,4] udy[3,5] udy[3,6] udy[3,7]
 udy[3,8] udy[3,9] udy[3,10] udy[3,11] udy[3,12] udy[3,13] udy[3,14]
 udy[3,15] udy[4,0] udy[4,1] udy[4,2] udy[4,3] udy[4,4] udy[4,5] udy[4,6]
 udy[4,7] udy[4,8] udy[4,9] udy[4,10] udy[4,11] udy[4,12] udy[4,13]
 udy[4,14] udy[4,15] udy[5,0] udy[5,1] udy[5,2] udy[5,3] udy[5,4] udy[5,5]
 udy[5,6] udy[5,7] udy[5,8] udy[5,9] udy[5,10] udy[5,11] udy[5,12]
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
 udy[11,15] udy[12,0] udy[12,1] udy[12,2] udy[12,3] udy[12,4] udy[12,5]
 udy[12,6] udy[12,7] udy[12,8] udy[12,9] udy[12,10] udy[12,11] udy[12,12]
 udy[12,13] udy[12,14] udy[12,15] udy[13,0] udy[13,1] udy[13,2] udy[13,3]
 udy[13,4] udy[13,5] udy[13,6] udy[13,7] udy[13,8] udy[13,9] udy[13,10]
 udy[13,11] udy[13,12] udy[13,13] udy[13,14] udy[13,15] udy[14,0] udy[14,1]
 udy[14,2] udy[14,3] udy[14,4] udy[14,5] udy[14,6] udy[14,7] udy[14,8]
 udy[14,9] udy[14,10] udy[14,11] udy[14,12] udy[14,13] udy[14,14]
 udy[14,15] udy[15,0] udy[15,1] udy[15,2] udy[15,3] udy[15,4] udy[15,5]
 udy[15,6] udy[15,7] udy[15,8] udy[15,9] udy[15,10] udy[15,11] udy[15,12]
 udy[15,13] udy[15,14] udy[15,15] udy[16,0] udy[16,1] udy[16,2] udy[16,3]
 udy[16,4] udy[16,5] udy[16,6] udy[16,7] udy[16,8] udy[16,9] udy[16,10]
 udy[16,11] udy[16,12] udy[16,13] udy[16,14] udy[16,15] uDZ[1,0] uDZ[1,1]
 uDZ[1,2] uDZ[1,3] uDZ[1,4] uDZ[1,5] uDZ[1,6] uDZ[1,7] uDZ[1,8] uDZ[1,9]
 uDZ[1,10] uDZ[1,11] uDZ[1,12] uDZ[1,13] uDZ[1,14] uDZ[1,15] uDZ[2,0]
 uDZ[2,1] uDZ[2,2] uDZ[2,3] uDZ[2,4] uDZ[2,5] uDZ[2,6] uDZ[2,7] uDZ[2,8]
 uDZ[2,9] uDZ[2,10] uDZ[2,11] uDZ[2,12] uDZ[2,13] uDZ[2,14] uDZ[2,15]
 uDZ[3,0] uDZ[3,1] uDZ[3,2] uDZ[3,3] uDZ[3,4] uDZ[3,5] uDZ[3,6] uDZ[3,7]
 uDZ[3,8] uDZ[3,9] uDZ[3,10] uDZ[3,11] uDZ[3,12] uDZ[3,13] uDZ[3,14]
 uDZ[3,15] uDZ[4,0] uDZ[4,1] uDZ[4,2] uDZ[4,3] uDZ[4,4] uDZ[4,5] uDZ[4,6]
 uDZ[4,7] uDZ[4,8] uDZ[4,9] uDZ[4,10] uDZ[4,11] uDZ[4,12] uDZ[4,13]
 uDZ[4,14] uDZ[4,15] uDZ[5,0] uDZ[5,1] uDZ[5,2] uDZ[5,3] uDZ[5,4] uDZ[5,5]
 uDZ[5,6] uDZ[5,7] uDZ[5,8] uDZ[5,9] uDZ[5,10] uDZ[5,11] uDZ[5,12]
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
 uDZ[11,15] uDZ[12,0] uDZ[12,1] uDZ[12,2] uDZ[12,3] uDZ[12,4] uDZ[12,5]
 uDZ[12,6] uDZ[12,7] uDZ[12,8] uDZ[12,9] uDZ[12,10] uDZ[12,11] uDZ[12,12]
 uDZ[12,13] uDZ[12,14] uDZ[12,15] uDZ[13,0] uDZ[13,1] uDZ[13,2] uDZ[13,3]
 uDZ[13,4] uDZ[13,5] uDZ[13,6] uDZ[13,7] uDZ[13,8] uDZ[13,9] uDZ[13,10]
 uDZ[13,11] uDZ[13,12] uDZ[13,13] uDZ[13,14] uDZ[13,15] uDZ[14,0] uDZ[14,1]
 uDZ[14,2] uDZ[14,3] uDZ[14,4] uDZ[14,5] uDZ[14,6] uDZ[14,7] uDZ[14,8]
 uDZ[14,9] uDZ[14,10] uDZ[14,11] uDZ[14,12] uDZ[14,13] uDZ[14,14]
 uDZ[14,15] uDZ[15,0] uDZ[15,1] uDZ[15,2] uDZ[15,3] uDZ[15,4] uDZ[15,5]
 uDZ[15,6] uDZ[15,7] uDZ[15,8] uDZ[15,9] uDZ[15,10] uDZ[15,11] uDZ[15,12]
 uDZ[15,13] uDZ[15,14] uDZ[15,15] uDZ[16,0] uDZ[16,1] uDZ[16,2] uDZ[16,3]
 uDZ[16,4] uDZ[16,5] uDZ[16,6] uDZ[16,7] uDZ[16,8] uDZ[16,9] uDZ[16,10]
 uDZ[16,11] uDZ[16,12] uDZ[16,13] uDZ[16,14] uDZ[16,15] udz[1,0] udz[1,1]
 udz[1,2] udz[1,3] udz[1,4] udz[1,5] udz[1,6] udz[1,7] udz[1,8] udz[1,9]
 udz[1,10] udz[1,11] udz[1,12] udz[1,13] udz[1,14] udz[1,15] udz[2,0]
 udz[2,1] udz[2,2] udz[2,3] udz[2,4] udz[2,5] udz[2,6] udz[2,7] udz[2,8]
 udz[2,9] udz[2,10] udz[2,11] udz[2,12] udz[2,13] udz[2,14] udz[2,15]
 udz[3,0] udz[3,1] udz[3,2] udz[3,3] udz[3,4] udz[3,5] udz[3,6] udz[3,7]
 udz[3,8] udz[3,9] udz[3,10] udz[3,11] udz[3,12] udz[3,13] udz[3,14]
 udz[3,15] udz[4,0] udz[4,1] udz[4,2] udz[4,3] udz[4,4] udz[4,5] udz[4,6]
 udz[4,7] udz[4,8] udz[4,9] udz[4,10] udz[4,11] udz[4,12] udz[4,13]
 udz[4,14] udz[4,15] udz[5,0] udz[5,1] udz[5,2] udz[5,3] udz[5,4] udz[5,5]
 udz[5,6] udz[5,7] udz[5,8] udz[5,9] udz[5,10] udz[5,11] udz[5,12]
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
 udz[11,15] udz[12,0] udz[12,1] udz[12,2] udz[12,3] udz[12,4] udz[12,5]
 udz[12,6] udz[12,7] udz[12,8] udz[12,9] udz[12,10] udz[12,11] udz[12,12]
 udz[12,13] udz[12,14] udz[12,15] udz[13,0] udz[13,1] udz[13,2] udz[13,3]
 udz[13,4] udz[13,5] udz[13,6] udz[13,7] udz[13,8] udz[13,9] udz[13,10]
 udz[13,11] udz[13,12] udz[13,13] udz[13,14] udz[13,15] udz[14,0] udz[14,1]
 udz[14,2] udz[14,3] udz[14,4] udz[14,5] udz[14,6] udz[14,7] udz[14,8]
 udz[14,9] udz[14,10] udz[14,11] udz[14,12] udz[14,13] udz[14,14]
 udz[14,15] udz[15,0] udz[15,1] udz[15,2] udz[15,3] udz[15,4] udz[15,5]
 udz[15,6] udz[15,7] udz[15,8] udz[15,9] udz[15,10] udz[15,11] udz[15,12]
 udz[15,13] udz[15,14] udz[15,15] udz[16,0] udz[16,1] udz[16,2] udz[16,3]
 udz[16,4] udz[16,5] udz[16,6] udz[16,7] udz[16,8] udz[16,9] udz[16,10]
 udz[16,11] udz[16,12] udz[16,13] udz[16,14] udz[16,15] uDW[0,0] uDW[0,1]
 uDW[0,2] uDW[0,3] uDW[0,4] uDW[0,5] uDW[0,6] uDW[0,7] uDW[0,8] uDW[0,9]
 uDW[0,10] uDW[0,11] uDW[0,12] uDW[0,13] uDW[0,14] uDW[0,15] uDW[1,0]
 uDW[1,1] uDW[1,2] uDW[1,3] uDW[1,4] uDW[1,5] uDW[1,6] uDW[1,7] uDW[1,8]
 uDW[1,9] uDW[1,10] uDW[1,11] uDW[1,12] uDW[1,13] uDW[1,14] uDW[1,15]
 uDW[2,0] uDW[2,1] uDW[2,2] uDW[2,3] uDW[2,4] uDW[2,5] uDW[2,6] uDW[2,7]
 uDW[2,8] uDW[2,9] uDW[2,10] uDW[2,11] uDW[2,12] uDW[2,13] uDW[2,14]
 uDW[2,15] uDW[3,0] uDW[3,1] uDW[3,2] uDW[3,3] uDW[3,4] uDW[3,5] uDW[3,6]
 uDW[3,7] uDW[3,8] uDW[3,9] uDW[3,10] uDW[3,11] uDW[3,12] uDW[3,13]
 uDW[3,14] uDW[3,15] uDW[4,0] uDW[4,1] uDW[4,2] uDW[4,3] uDW[4,4] uDW[4,5]
 uDW[4,6] uDW[4,7] uDW[4,8] uDW[4,9] uDW[4,10] uDW[4,11] uDW[4,12]
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
 uDW[11,15] uDW[12,0] uDW[12,1] uDW[12,2] uDW[12,3] uDW[12,4] uDW[12,5]
 uDW[12,6] uDW[12,7] uDW[12,8] uDW[12,9] uDW[12,10] uDW[12,11] uDW[12,12]
 uDW[12,13] uDW[12,14] uDW[12,15] uDW[13,0] uDW[13,1] uDW[13,2] uDW[13,3]
 uDW[13,4] uDW[13,5] uDW[13,6] uDW[13,7] uDW[13,8] uDW[13,9] uDW[13,10]
 uDW[13,11] uDW[13,12] uDW[13,13] uDW[13,14] uDW[13,15] uDW[14,0] uDW[14,1]
 uDW[14,2] uDW[14,3] uDW[14,4] uDW[14,5] uDW[14,6] uDW[14,7] uDW[14,8]
 uDW[14,9] uDW[14,10] uDW[14,11] uDW[14,12] uDW[14,13] uDW[14,14]
 uDW[14,15] uDW[15,0] uDW[15,1] uDW[15,2] uDW[15,3] uDW[15,4] uDW[15,5]
 uDW[15,6] uDW[15,7] uDW[15,8] uDW[15,9] uDW[15,10] uDW[15,11] uDW[15,12]
 uDW[15,13] uDW[15,14] uDW[15,15] uDW[16,0] uDW[16,1] uDW[16,2] uDW[16,3]
 uDW[16,4] uDW[16,5] uDW[16,6] uDW[16,7] uDW[16,8] uDW[16,9] uDW[16,10]
 uDW[16,11] uDW[16,12] uDW[16,13] uDW[16,14] uDW[16,15] udw[0,0] udw[0,1]
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
 udw[7,13] udw[7,14] udw[7,15] udw[8,0] udw[8,1] udw[8,2] udw[8,3] udw[8,4]
 udw[8,5] udw[8,6] udw[8,7] udw[8,8] udw[8,9] udw[8,10] udw[8,11] udw[8,12]
 udw[8,13] udw[8,14] udw[8,15] udw[9,0] udw[9,1] udw[9,2] udw[9,3] udw[9,4]
 udw[9,5] udw[9,6] udw[9,7] udw[9,8] udw[9,9] udw[9,10] udw[9,11] udw[9,12]
 udw[9,13] udw[9,14] udw[9,15] udw[10,0] udw[10,1] udw[10,2] udw[10,3]
 udw[10,4] udw[10,5] udw[10,6] udw[10,7] udw[10,8] udw[10,9] udw[10,10]
 udw[10,11] udw[10,12] udw[10,13] udw[10,14] udw[10,15] udw[11,0] udw[11,1]
 udw[11,2] udw[11,3] udw[11,4] udw[11,5] udw[11,6] udw[11,7] udw[11,8]
 udw[11,9] udw[11,10] udw[11,11] udw[11,12] udw[11,13] udw[11,14]
 udw[11,15] udw[12,0] udw[12,1] udw[12,2] udw[12,3] udw[12,4] udw[12,5]
 udw[12,6] udw[12,7] udw[12,8] udw[12,9] udw[12,10] udw[12,11] udw[12,12]
 udw[12,13] udw[12,14] udw[12,15] udw[13,0] udw[13,1] udw[13,2] udw[13,3]
 udw[13,4] udw[13,5] udw[13,6] udw[13,7] udw[13,8] udw[13,9] udw[13,10]
 udw[13,11] udw[13,12] udw[13,13] udw[13,14] udw[13,15] udw[14,0] udw[14,1]
 udw[14,2] udw[14,3] udw[14,4] udw[14,5] udw[14,6] udw[14,7] udw[14,8]
 udw[14,9] udw[14,10] udw[14,11] udw[14,12] udw[14,13] udw[14,14]
 udw[14,15] udw[15,0] udw[15,1] udw[15,2] udw[15,3] udw[15,4] udw[15,5]
 udw[15,6] udw[15,7] udw[15,8] udw[15,9] udw[15,10] udw[15,11] udw[15,12]
 udw[15,13] udw[15,14] udw[15,15] udw[16,0] udw[16,1] udw[16,2] udw[16,3]
 udw[16,4] udw[16,5] udw[16,6] udw[16,7] udw[16,8] udw[16,9] udw[16,10]
 udw[16,11] udw[16,12] udw[16,13] udw[16,14] udw[16,15] ustk[0,0] ustk[0,1]
 ustk[0,2] ustk[0,3] ustk[0,4] ustk[0,5] ustk[0,6] ustk[0,7] ustk[0,8]
 ustk[0,9] ustk[0,10] ustk[0,11] ustk[0,12] ustk[0,13] ustk[0,14]
 ustk[0,15] ustk[1,0] ustk[1,1] ustk[1,2] ustk[1,3] ustk[1,4] ustk[1,5]
 ustk[1,6] ustk[1,7] ustk[1,8] ustk[1,9] ustk[1,10] ustk[1,11] ustk[1,12]
 ustk[1,13] ustk[1,14] ustk[1,15] ustk[2,0] ustk[2,1] ustk[2,2] ustk[2,3]
 ustk[2,4] ustk[2,5] ustk[2,6] ustk[2,7] ustk[2,8] ustk[2,9] ustk[2,10]
 ustk[2,11] ustk[2,12] ustk[2,13] ustk[2,14] ustk[2,15] ustk[3,0] ustk[3,1]
 ustk[3,2] ustk[3,3] ustk[3,4] ustk[3,5] ustk[3,6] ustk[3,7] ustk[3,8]
 ustk[3,9] ustk[3,10] ustk[3,11] ustk[3,12] ustk[3,13] ustk[3,14]
 ustk[3,15] ustk[4,0] ustk[4,1] ustk[4,2] ustk[4,3] ustk[4,4] ustk[4,5]
 ustk[4,6] ustk[4,7] ustk[4,8] ustk[4,9] ustk[4,10] ustk[4,11] ustk[4,12]
 ustk[4,13] ustk[4,14] ustk[4,15] ustk[5,0] ustk[5,1] ustk[5,2] ustk[5,3]
 ustk[5,4] ustk[5,5] ustk[5,6] ustk[5,7] ustk[5,8] ustk[5,9] ustk[5,10]
 ustk[5,11] ustk[5,12] ustk[5,13] ustk[5,14] ustk[5,15] ustk[6,0] ustk[6,1]
 ustk[6,2] ustk[6,3] ustk[6,4] ustk[6,5] ustk[6,6] ustk[6,7] ustk[6,8]
 ustk[6,9] ustk[6,10] ustk[6,11] ustk[6,12] ustk[6,13] ustk[6,14]
 ustk[6,15] ustk[7,0] ustk[7,1] ustk[7,2] ustk[7,3] ustk[7,4] ustk[7,5]
 ustk[7,6] ustk[7,7] ustk[7,8] ustk[7,9] ustk[7,10] ustk[7,11] ustk[7,12]
 ustk[7,13] ustk[7,14] ustk[7,15] ustk[8,0] ustk[8,1] ustk[8,2] ustk[8,3]
 ustk[8,4] ustk[8,5] ustk[8,6] ustk[8,7] ustk[8,8] ustk[8,9] ustk[8,10]
 ustk[8,11] ustk[8,12] ustk[8,13] ustk[8,14] ustk[8,15] ustk[9,0] ustk[9,1]
 ustk[9,2] ustk[9,3] ustk[9,4] ustk[9,5] ustk[9,6] ustk[9,7] ustk[9,8]
 ustk[9,9] ustk[9,10] ustk[9,11] ustk[9,12] ustk[9,13] ustk[9,14]
 ustk[9,15] ustk[10,0] ustk[10,1] ustk[10,2] ustk[10,3] ustk[10,4]
 ustk[10,5] ustk[10,6] ustk[10,7] ustk[10,8] ustk[10,9] ustk[10,10]
 ustk[10,11] ustk[10,12] ustk[10,13] ustk[10,14] ustk[10,15] ustk[11,0]
 ustk[11,1] ustk[11,2] ustk[11,3] ustk[11,4] ustk[11,5] ustk[11,6]
 ustk[11,7] ustk[11,8] ustk[11,9] ustk[11,10] ustk[11,11] ustk[11,12]
 ustk[11,13] ustk[11,14] ustk[11,15] ustk[12,0] ustk[12,1] ustk[12,2]
 ustk[12,3] ustk[12,4] ustk[12,5] ustk[12,6] ustk[12,7] ustk[12,8]
 ustk[12,9] ustk[12,10] ustk[12,11] ustk[12,12] ustk[12,13] ustk[12,14]
 ustk[12,15] ustk[13,0] ustk[13,1] ustk[13,2] ustk[13,3] ustk[13,4]
 ustk[13,5] ustk[13,6] ustk[13,7] ustk[13,8] ustk[13,9] ustk[13,10]
 ustk[13,11] ustk[13,12] ustk[13,13] ustk[13,14] ustk[13,15] ustk[14,0]
 ustk[14,1] ustk[14,2] ustk[14,3] ustk[14,4] ustk[14,5] ustk[14,6]
 ustk[14,7] ustk[14,8] ustk[14,9] ustk[14,10] ustk[14,11] ustk[14,12]
 ustk[14,13] ustk[14,14] ustk[14,15] ustk[15,0] ustk[15,1] ustk[15,2]
 ustk[15,3] ustk[15,4] ustk[15,5] ustk[15,6] ustk[15,7] ustk[15,8]
 ustk[15,9] ustk[15,10] ustk[15,11] ustk[15,12] ustk[15,13] ustk[15,14]
 ustk[15,15] ustk[16,0] ustk[16,1] ustk[16,2] ustk[16,3] ustk[16,4]
 ustk[16,5] ustk[16,6] ustk[16,7] ustk[16,8] ustk[16,9] ustk[16,10]
 ustk[16,11] ustk[16,12] ustk[16,13] ustk[16,14] ustk[16,15] eqk[0] eqk[1]
 eqk[2] eqk[3] eqk[4] eqk[5] eqk[6] eqk[7] eqk[8] eqk[9] eqk[10] eqk[11]
 eqk[12] eqk[13] eqk[14] eqk[15] uLANE[0] uLANE[1] uLANE[2] uLANE[3]
 uLANE[4] uLANE[5] uLANE[6] uLANE[7] uLANE[8] uLANE[9] uLANE[10] uLANE[11]
 uLANE[12] uLANE[13] uLANE[14] uLANE[15] lDX[18,0] lDX[18,1] lDX[18,2]
 lDX[18,3] lDX[18,4] lDX[18,5] lDX[18,6] lDX[18,7] lDX[18,8] lDX[18,9]
 lDX[18,10] lDX[18,11] lDX[18,12] lDX[18,13] lDX[18,14] lDX[18,15]
 lDX[19,0] lDX[19,1] lDX[19,2] lDX[19,3] lDX[19,4] lDX[19,5] lDX[19,6]
 lDX[19,7] lDX[19,8] lDX[19,9] lDX[19,10] lDX[19,11] lDX[19,12] lDX[19,13]
 lDX[19,14] lDX[19,15] lDX[20,0] lDX[20,1] lDX[20,2] lDX[20,3] lDX[20,4]
 lDX[20,5] lDX[20,6] lDX[20,7] lDX[20,8] lDX[20,9] lDX[20,10] lDX[20,11]
 lDX[20,12] lDX[20,13] lDX[20,14] lDX[20,15] lDX[21,0] lDX[21,1] lDX[21,2]
 lDX[21,3] lDX[21,4] lDX[21,5] lDX[21,6] lDX[21,7] lDX[21,8] lDX[21,9]
 lDX[21,10] lDX[21,11] lDX[21,12] lDX[21,13] lDX[21,14] lDX[21,15]
 lDX[22,0] lDX[22,1] lDX[22,2] lDX[22,3] lDX[22,4] lDX[22,5] lDX[22,6]
 lDX[22,7] lDX[22,8] lDX[22,9] lDX[22,10] lDX[22,11] lDX[22,12] lDX[22,13]
 lDX[22,14] lDX[22,15] lDX[23,0] lDX[23,1] lDX[23,2] lDX[23,3] lDX[23,4]
 lDX[23,5] lDX[23,6] lDX[23,7] lDX[23,8] lDX[23,9] lDX[23,10] lDX[23,11]
 lDX[23,12] lDX[23,13] lDX[23,14] lDX[23,15] lDX[24,0] lDX[24,1] lDX[24,2]
 lDX[24,3] lDX[24,4] lDX[24,5] lDX[24,6] lDX[24,7] lDX[24,8] lDX[24,9]
 lDX[24,10] lDX[24,11] lDX[24,12] lDX[24,13] lDX[24,14] lDX[24,15]
 lDX[25,0] lDX[25,1] lDX[25,2] lDX[25,3] lDX[25,4] lDX[25,5] lDX[25,6]
 lDX[25,7] lDX[25,8] lDX[25,9] lDX[25,10] lDX[25,11] lDX[25,12] lDX[25,13]
 lDX[25,14] lDX[25,15] lDX[26,0] lDX[26,1] lDX[26,2] lDX[26,3] lDX[26,4]
 lDX[26,5] lDX[26,6] lDX[26,7] lDX[26,8] lDX[26,9] lDX[26,10] lDX[26,11]
 lDX[26,12] lDX[26,13] lDX[26,14] lDX[26,15] lDX[27,0] lDX[27,1] lDX[27,2]
 lDX[27,3] lDX[27,4] lDX[27,5] lDX[27,6] lDX[27,7] lDX[27,8] lDX[27,9]
 lDX[27,10] lDX[27,11] lDX[27,12] lDX[27,13] lDX[27,14] lDX[27,15]
 lDX[28,0] lDX[28,1] lDX[28,2] lDX[28,3] lDX[28,4] lDX[28,5] lDX[28,6]
 lDX[28,7] lDX[28,8] lDX[28,9] lDX[28,10] lDX[28,11] lDX[28,12] lDX[28,13]
 lDX[28,14] lDX[28,15] lDX[29,0] lDX[29,1] lDX[29,2] lDX[29,3] lDX[29,4]
 lDX[29,5] lDX[29,6] lDX[29,7] lDX[29,8] lDX[29,9] lDX[29,10] lDX[29,11]
 lDX[29,12] lDX[29,13] lDX[29,14] lDX[29,15] lDX[30,0] lDX[30,1] lDX[30,2]
 lDX[30,3] lDX[30,4] lDX[30,5] lDX[30,6] lDX[30,7] lDX[30,8] lDX[30,9]
 lDX[30,10] lDX[30,11] lDX[30,12] lDX[30,13] lDX[30,14] lDX[30,15]
 lDX[31,0] lDX[31,1] lDX[31,2] lDX[31,3] lDX[31,4] lDX[31,5] lDX[31,6]
 lDX[31,7] lDX[31,8] lDX[31,9] lDX[31,10] lDX[31,11] lDX[31,12] lDX[31,13]
 lDX[31,14] lDX[31,15] lDX[32,0] lDX[32,1] lDX[32,2] lDX[32,3] lDX[32,4]
 lDX[32,5] lDX[32,6] lDX[32,7] lDX[32,8] lDX[32,9] lDX[32,10] lDX[32,11]
 lDX[32,12] lDX[32,13] lDX[32,14] lDX[32,15] lDX[33,0] lDX[33,1] lDX[33,2]
 lDX[33,3] lDX[33,4] lDX[33,5] lDX[33,6] lDX[33,7] lDX[33,8] lDX[33,9]
 lDX[33,10] lDX[33,11] lDX[33,12] lDX[33,13] lDX[33,14] lDX[33,15]
 ldx[18,0] ldx[18,1] ldx[18,2] ldx[18,3] ldx[18,4] ldx[18,5] ldx[18,6]
 ldx[18,7] ldx[18,8] ldx[18,9] ldx[18,10] ldx[18,11] ldx[18,12] ldx[18,13]
 ldx[18,14] ldx[18,15] ldx[19,0] ldx[19,1] ldx[19,2] ldx[19,3] ldx[19,4]
 ldx[19,5] ldx[19,6] ldx[19,7] ldx[19,8] ldx[19,9] ldx[19,10] ldx[19,11]
 ldx[19,12] ldx[19,13] ldx[19,14] ldx[19,15] ldx[20,0] ldx[20,1] ldx[20,2]
 ldx[20,3] ldx[20,4] ldx[20,5] ldx[20,6] ldx[20,7] ldx[20,8] ldx[20,9]
 ldx[20,10] ldx[20,11] ldx[20,12] ldx[20,13] ldx[20,14] ldx[20,15]
 ldx[21,0] ldx[21,1] ldx[21,2] ldx[21,3] ldx[21,4] ldx[21,5] ldx[21,6]
 ldx[21,7] ldx[21,8] ldx[21,9] ldx[21,10] ldx[21,11] ldx[21,12] ldx[21,13]
 ldx[21,14] ldx[21,15] ldx[22,0] ldx[22,1] ldx[22,2] ldx[22,3] ldx[22,4]
 ldx[22,5] ldx[22,6] ldx[22,7] ldx[22,8] ldx[22,9] ldx[22,10] ldx[22,11]
 ldx[22,12] ldx[22,13] ldx[22,14] ldx[22,15] ldx[23,0] ldx[23,1] ldx[23,2]
 ldx[23,3] ldx[23,4] ldx[23,5] ldx[23,6] ldx[23,7] ldx[23,8] ldx[23,9]
 ldx[23,10] ldx[23,11] ldx[23,12] ldx[23,13] ldx[23,14] ldx[23,15]
 ldx[24,0] ldx[24,1] ldx[24,2] ldx[24,3] ldx[24,4] ldx[24,5] ldx[24,6]
 ldx[24,7] ldx[24,8] ldx[24,9] ldx[24,10] ldx[24,11] ldx[24,12] ldx[24,13]
 ldx[24,14] ldx[24,15] ldx[25,0] ldx[25,1] ldx[25,2] ldx[25,3] ldx[25,4]
 ldx[25,5] ldx[25,6] ldx[25,7] ldx[25,8] ldx[25,9] ldx[25,10] ldx[25,11]
 ldx[25,12] ldx[25,13] ldx[25,14] ldx[25,15] ldx[26,0] ldx[26,1] ldx[26,2]
 ldx[26,3] ldx[26,4] ldx[26,5] ldx[26,6] ldx[26,7] ldx[26,8] ldx[26,9]
 ldx[26,10] ldx[26,11] ldx[26,12] ldx[26,13] ldx[26,14] ldx[26,15]
 ldx[27,0] ldx[27,1] ldx[27,2] ldx[27,3] ldx[27,4] ldx[27,5] ldx[27,6]
 ldx[27,7] ldx[27,8] ldx[27,9] ldx[27,10] ldx[27,11] ldx[27,12] ldx[27,13]
 ldx[27,14] ldx[27,15] ldx[28,0] ldx[28,1] ldx[28,2] ldx[28,3] ldx[28,4]
 ldx[28,5] ldx[28,6] ldx[28,7] ldx[28,8] ldx[28,9] ldx[28,10] ldx[28,11]
 ldx[28,12] ldx[28,13] ldx[28,14] ldx[28,15] ldx[29,0] ldx[29,1] ldx[29,2]
 ldx[29,3] ldx[29,4] ldx[29,5] ldx[29,6] ldx[29,7] ldx[29,8] ldx[29,9]
 ldx[29,10] ldx[29,11] ldx[29,12] ldx[29,13] ldx[29,14] ldx[29,15]
 ldx[30,0] ldx[30,1] ldx[30,2] ldx[30,3] ldx[30,4] ldx[30,5] ldx[30,6]
 ldx[30,7] ldx[30,8] ldx[30,9] ldx[30,10] ldx[30,11] ldx[30,12] ldx[30,13]
 ldx[30,14] ldx[30,15] ldx[31,0] ldx[31,1] ldx[31,2] ldx[31,3] ldx[31,4]
 ldx[31,5] ldx[31,6] ldx[31,7] ldx[31,8] ldx[31,9] ldx[31,10] ldx[31,11]
 ldx[31,12] ldx[31,13] ldx[31,14] ldx[31,15] ldx[32,0] ldx[32,1] ldx[32,2]
 ldx[32,3] ldx[32,4] ldx[32,5] ldx[32,6] ldx[32,7] ldx[32,8] ldx[32,9]
 ldx[32,10] ldx[32,11] ldx[32,12] ldx[32,13] ldx[32,14] ldx[32,15]
 ldx[33,0] ldx[33,1] ldx[33,2] ldx[33,3] ldx[33,4] ldx[33,5] ldx[33,6]
 ldx[33,7] ldx[33,8] ldx[33,9] ldx[33,10] ldx[33,11] ldx[33,12] ldx[33,13]
 ldx[33,14] ldx[33,15] lDY[17,0] lDY[17,1] lDY[17,2] lDY[17,3] lDY[17,4]
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
 lDY[23,12] lDY[23,13] lDY[23,14] lDY[23,15] lDY[24,0] lDY[24,1] lDY[24,2]
 lDY[24,3] lDY[24,4] lDY[24,5] lDY[24,6] lDY[24,7] lDY[24,8] lDY[24,9]
 lDY[24,10] lDY[24,11] lDY[24,12] lDY[24,13] lDY[24,14] lDY[24,15]
 lDY[25,0] lDY[25,1] lDY[25,2] lDY[25,3] lDY[25,4] lDY[25,5] lDY[25,6]
 lDY[25,7] lDY[25,8] lDY[25,9] lDY[25,10] lDY[25,11] lDY[25,12] lDY[25,13]
 lDY[25,14] lDY[25,15] lDY[26,0] lDY[26,1] lDY[26,2] lDY[26,3] lDY[26,4]
 lDY[26,5] lDY[26,6] lDY[26,7] lDY[26,8] lDY[26,9] lDY[26,10] lDY[26,11]
 lDY[26,12] lDY[26,13] lDY[26,14] lDY[26,15] lDY[27,0] lDY[27,1] lDY[27,2]
 lDY[27,3] lDY[27,4] lDY[27,5] lDY[27,6] lDY[27,7] lDY[27,8] lDY[27,9]
 lDY[27,10] lDY[27,11] lDY[27,12] lDY[27,13] lDY[27,14] lDY[27,15]
 lDY[28,0] lDY[28,1] lDY[28,2] lDY[28,3] lDY[28,4] lDY[28,5] lDY[28,6]
 lDY[28,7] lDY[28,8] lDY[28,9] lDY[28,10] lDY[28,11] lDY[28,12] lDY[28,13]
 lDY[28,14] lDY[28,15] lDY[29,0] lDY[29,1] lDY[29,2] lDY[29,3] lDY[29,4]
 lDY[29,5] lDY[29,6] lDY[29,7] lDY[29,8] lDY[29,9] lDY[29,10] lDY[29,11]
 lDY[29,12] lDY[29,13] lDY[29,14] lDY[29,15] lDY[30,0] lDY[30,1] lDY[30,2]
 lDY[30,3] lDY[30,4] lDY[30,5] lDY[30,6] lDY[30,7] lDY[30,8] lDY[30,9]
 lDY[30,10] lDY[30,11] lDY[30,12] lDY[30,13] lDY[30,14] lDY[30,15]
 lDY[31,0] lDY[31,1] lDY[31,2] lDY[31,3] lDY[31,4] lDY[31,5] lDY[31,6]
 lDY[31,7] lDY[31,8] lDY[31,9] lDY[31,10] lDY[31,11] lDY[31,12] lDY[31,13]
 lDY[31,14] lDY[31,15] lDY[32,0] lDY[32,1] lDY[32,2] lDY[32,3] lDY[32,4]
 lDY[32,5] lDY[32,6] lDY[32,7] lDY[32,8] lDY[32,9] lDY[32,10] lDY[32,11]
 lDY[32,12] lDY[32,13] lDY[32,14] lDY[32,15] ldy[17,0] ldy[17,1] ldy[17,2]
 ldy[17,3] ldy[17,4] ldy[17,5] ldy[17,6] ldy[17,7] ldy[17,8] ldy[17,9]
 ldy[17,10] ldy[17,11] ldy[17,12] ldy[17,13] ldy[17,14] ldy[17,15]
 ldy[18,0] ldy[18,1] ldy[18,2] ldy[18,3] ldy[18,4] ldy[18,5] ldy[18,6]
 ldy[18,7] ldy[18,8] ldy[18,9] ldy[18,10] ldy[18,11] ldy[18,12] ldy[18,13]
 ldy[18,14] ldy[18,15] ldy[19,0] ldy[19,1] ldy[19,2] ldy[19,3] ldy[19,4]
 ldy[19,5] ldy[19,6] ldy[19,7] ldy[19,8] ldy[19,9] ldy[19,10] ldy[19,11]
 ldy[19,12] ldy[19,13] ldy[19,14] ldy[19,15] ldy[20,0] ldy[20,1] ldy[20,2]
 ldy[20,3] ldy[20,4] ldy[20,5] ldy[20,6] ldy[20,7] ldy[20,8] ldy[20,9]
 ldy[20,10] ldy[20,11] ldy[20,12] ldy[20,13] ldy[20,14] ldy[20,15]
 ldy[21,0] ldy[21,1] ldy[21,2] ldy[21,3] ldy[21,4] ldy[21,5] ldy[21,6]
 ldy[21,7] ldy[21,8] ldy[21,9] ldy[21,10] ldy[21,11] ldy[21,12] ldy[21,13]
 ldy[21,14] ldy[21,15] ldy[22,0] ldy[22,1] ldy[22,2] ldy[22,3] ldy[22,4]
 ldy[22,5] ldy[22,6] ldy[22,7] ldy[22,8] ldy[22,9] ldy[22,10] ldy[22,11]
 ldy[22,12] ldy[22,13] ldy[22,14] ldy[22,15] ldy[23,0] ldy[23,1] ldy[23,2]
 ldy[23,3] ldy[23,4] ldy[23,5] ldy[23,6] ldy[23,7] ldy[23,8] ldy[23,9]
 ldy[23,10] ldy[23,11] ldy[23,12] ldy[23,13] ldy[23,14] ldy[23,15]
 ldy[24,0] ldy[24,1] ldy[24,2] ldy[24,3] ldy[24,4] ldy[24,5] ldy[24,6]
 ldy[24,7] ldy[24,8] ldy[24,9] ldy[24,10] ldy[24,11] ldy[24,12] ldy[24,13]
 ldy[24,14] ldy[24,15] ldy[25,0] ldy[25,1] ldy[25,2] ldy[25,3] ldy[25,4]
 ldy[25,5] ldy[25,6] ldy[25,7] ldy[25,8] ldy[25,9] ldy[25,10] ldy[25,11]
 ldy[25,12] ldy[25,13] ldy[25,14] ldy[25,15] ldy[26,0] ldy[26,1] ldy[26,2]
 ldy[26,3] ldy[26,4] ldy[26,5] ldy[26,6] ldy[26,7] ldy[26,8] ldy[26,9]
 ldy[26,10] ldy[26,11] ldy[26,12] ldy[26,13] ldy[26,14] ldy[26,15]
 ldy[27,0] ldy[27,1] ldy[27,2] ldy[27,3] ldy[27,4] ldy[27,5] ldy[27,6]
 ldy[27,7] ldy[27,8] ldy[27,9] ldy[27,10] ldy[27,11] ldy[27,12] ldy[27,13]
 ldy[27,14] ldy[27,15] ldy[28,0] ldy[28,1] ldy[28,2] ldy[28,3] ldy[28,4]
 ldy[28,5] ldy[28,6] ldy[28,7] ldy[28,8] ldy[28,9] ldy[28,10] ldy[28,11]
 ldy[28,12] ldy[28,13] ldy[28,14] ldy[28,15] ldy[29,0] ldy[29,1] ldy[29,2]
 ldy[29,3] ldy[29,4] ldy[29,5] ldy[29,6] ldy[29,7] ldy[29,8] ldy[29,9]
 ldy[29,10] ldy[29,11] ldy[29,12] ldy[29,13] ldy[29,14] ldy[29,15]
 ldy[30,0] ldy[30,1] ldy[30,2] ldy[30,3] ldy[30,4] ldy[30,5] ldy[30,6]
 ldy[30,7] ldy[30,8] ldy[30,9] ldy[30,10] ldy[30,11] ldy[30,12] ldy[30,13]
 ldy[30,14] ldy[30,15] ldy[31,0] ldy[31,1] ldy[31,2] ldy[31,3] ldy[31,4]
 ldy[31,5] ldy[31,6] ldy[31,7] ldy[31,8] ldy[31,9] ldy[31,10] ldy[31,11]
 ldy[31,12] ldy[31,13] ldy[31,14] ldy[31,15] ldy[32,0] ldy[32,1] ldy[32,2]
 ldy[32,3] ldy[32,4] ldy[32,5] ldy[32,6] ldy[32,7] ldy[32,8] ldy[32,9]
 ldy[32,10] ldy[32,11] ldy[32,12] ldy[32,13] ldy[32,14] ldy[32,15]
 lDZ[17,0] lDZ[17,1] lDZ[17,2] lDZ[17,3] lDZ[17,4] lDZ[17,5] lDZ[17,6]
 lDZ[17,7] lDZ[17,8] lDZ[17,9] lDZ[17,10] lDZ[17,11] lDZ[17,12] lDZ[17,13]
 lDZ[17,14] lDZ[17,15] lDZ[18,0] lDZ[18,1] lDZ[18,2] lDZ[18,3] lDZ[18,4]
 lDZ[18,5] lDZ[18,6] lDZ[18,7] lDZ[18,8] lDZ[18,9] lDZ[18,10] lDZ[18,11]
 lDZ[18,12] lDZ[18,13] lDZ[18,14] lDZ[18,15] lDZ[19,0] lDZ[19,1] lDZ[19,2]
 lDZ[19,3] lDZ[19,4] lDZ[19,5] lDZ[19,6] lDZ[19,7] lDZ[19,8] lDZ[19,9]
 lDZ[19,10] lDZ[19,11] lDZ[19,12] lDZ[19,13] lDZ[19,14] lDZ[19,15]
 lDZ[20,0] lDZ[20,1] lDZ[20,2] lDZ[20,3] lDZ[20,4] lDZ[20,5] lDZ[20,6]
 lDZ[20,7] lDZ[20,8] lDZ[20,9] lDZ[20,10] lDZ[20,11] lDZ[20,12] lDZ[20,13]
 lDZ[20,14] lDZ[20,15] lDZ[21,0] lDZ[21,1] lDZ[21,2] lDZ[21,3] lDZ[21,4]
 lDZ[21,5] lDZ[21,6] lDZ[21,7] lDZ[21,8] lDZ[21,9] lDZ[21,10] lDZ[21,11]
 lDZ[21,12] lDZ[21,13] lDZ[21,14] lDZ[21,15] lDZ[22,0] lDZ[22,1] lDZ[22,2]
 lDZ[22,3] lDZ[22,4] lDZ[22,5] lDZ[22,6] lDZ[22,7] lDZ[22,8] lDZ[22,9]
 lDZ[22,10] lDZ[22,11] lDZ[22,12] lDZ[22,13] lDZ[22,14] lDZ[22,15]
 lDZ[23,0] lDZ[23,1] lDZ[23,2] lDZ[23,3] lDZ[23,4] lDZ[23,5] lDZ[23,6]
 lDZ[23,7] lDZ[23,8] lDZ[23,9] lDZ[23,10] lDZ[23,11] lDZ[23,12] lDZ[23,13]
 lDZ[23,14] lDZ[23,15] lDZ[24,0] lDZ[24,1] lDZ[24,2] lDZ[24,3] lDZ[24,4]
 lDZ[24,5] lDZ[24,6] lDZ[24,7] lDZ[24,8] lDZ[24,9] lDZ[24,10] lDZ[24,11]
 lDZ[24,12] lDZ[24,13] lDZ[24,14] lDZ[24,15] lDZ[25,0] lDZ[25,1] lDZ[25,2]
 lDZ[25,3] lDZ[25,4] lDZ[25,5] lDZ[25,6] lDZ[25,7] lDZ[25,8] lDZ[25,9]
 lDZ[25,10] lDZ[25,11] lDZ[25,12] lDZ[25,13] lDZ[25,14] lDZ[25,15]
 lDZ[26,0] lDZ[26,1] lDZ[26,2] lDZ[26,3] lDZ[26,4] lDZ[26,5] lDZ[26,6]
 lDZ[26,7] lDZ[26,8] lDZ[26,9] lDZ[26,10] lDZ[26,11] lDZ[26,12] lDZ[26,13]
 lDZ[26,14] lDZ[26,15] lDZ[27,0] lDZ[27,1] lDZ[27,2] lDZ[27,3] lDZ[27,4]
 lDZ[27,5] lDZ[27,6] lDZ[27,7] lDZ[27,8] lDZ[27,9] lDZ[27,10] lDZ[27,11]
 lDZ[27,12] lDZ[27,13] lDZ[27,14] lDZ[27,15] lDZ[28,0] lDZ[28,1] lDZ[28,2]
 lDZ[28,3] lDZ[28,4] lDZ[28,5] lDZ[28,6] lDZ[28,7] lDZ[28,8] lDZ[28,9]
 lDZ[28,10] lDZ[28,11] lDZ[28,12] lDZ[28,13] lDZ[28,14] lDZ[28,15]
 lDZ[29,0] lDZ[29,1] lDZ[29,2] lDZ[29,3] lDZ[29,4] lDZ[29,5] lDZ[29,6]
 lDZ[29,7] lDZ[29,8] lDZ[29,9] lDZ[29,10] lDZ[29,11] lDZ[29,12] lDZ[29,13]
 lDZ[29,14] lDZ[29,15] lDZ[30,0] lDZ[30,1] lDZ[30,2] lDZ[30,3] lDZ[30,4]
 lDZ[30,5] lDZ[30,6] lDZ[30,7] lDZ[30,8] lDZ[30,9] lDZ[30,10] lDZ[30,11]
 lDZ[30,12] lDZ[30,13] lDZ[30,14] lDZ[30,15] lDZ[31,0] lDZ[31,1] lDZ[31,2]
 lDZ[31,3] lDZ[31,4] lDZ[31,5] lDZ[31,6] lDZ[31,7] lDZ[31,8] lDZ[31,9]
 lDZ[31,10] lDZ[31,11] lDZ[31,12] lDZ[31,13] lDZ[31,14] lDZ[31,15]
 lDZ[32,0] lDZ[32,1] lDZ[32,2] lDZ[32,3] lDZ[32,4] lDZ[32,5] lDZ[32,6]
 lDZ[32,7] lDZ[32,8] lDZ[32,9] lDZ[32,10] lDZ[32,11] lDZ[32,12] lDZ[32,13]
 lDZ[32,14] lDZ[32,15] ldz[17,0] ldz[17,1] ldz[17,2] ldz[17,3] ldz[17,4]
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
 ldz[23,12] ldz[23,13] ldz[23,14] ldz[23,15] ldz[24,0] ldz[24,1] ldz[24,2]
 ldz[24,3] ldz[24,4] ldz[24,5] ldz[24,6] ldz[24,7] ldz[24,8] ldz[24,9]
 ldz[24,10] ldz[24,11] ldz[24,12] ldz[24,13] ldz[24,14] ldz[24,15]
 ldz[25,0] ldz[25,1] ldz[25,2] ldz[25,3] ldz[25,4] ldz[25,5] ldz[25,6]
 ldz[25,7] ldz[25,8] ldz[25,9] ldz[25,10] ldz[25,11] ldz[25,12] ldz[25,13]
 ldz[25,14] ldz[25,15] ldz[26,0] ldz[26,1] ldz[26,2] ldz[26,3] ldz[26,4]
 ldz[26,5] ldz[26,6] ldz[26,7] ldz[26,8] ldz[26,9] ldz[26,10] ldz[26,11]
 ldz[26,12] ldz[26,13] ldz[26,14] ldz[26,15] ldz[27,0] ldz[27,1] ldz[27,2]
 ldz[27,3] ldz[27,4] ldz[27,5] ldz[27,6] ldz[27,7] ldz[27,8] ldz[27,9]
 ldz[27,10] ldz[27,11] ldz[27,12] ldz[27,13] ldz[27,14] ldz[27,15]
 ldz[28,0] ldz[28,1] ldz[28,2] ldz[28,3] ldz[28,4] ldz[28,5] ldz[28,6]
 ldz[28,7] ldz[28,8] ldz[28,9] ldz[28,10] ldz[28,11] ldz[28,12] ldz[28,13]
 ldz[28,14] ldz[28,15] ldz[29,0] ldz[29,1] ldz[29,2] ldz[29,3] ldz[29,4]
 ldz[29,5] ldz[29,6] ldz[29,7] ldz[29,8] ldz[29,9] ldz[29,10] ldz[29,11]
 ldz[29,12] ldz[29,13] ldz[29,14] ldz[29,15] ldz[30,0] ldz[30,1] ldz[30,2]
 ldz[30,3] ldz[30,4] ldz[30,5] ldz[30,6] ldz[30,7] ldz[30,8] ldz[30,9]
 ldz[30,10] ldz[30,11] ldz[30,12] ldz[30,13] ldz[30,14] ldz[30,15]
 ldz[31,0] ldz[31,1] ldz[31,2] ldz[31,3] ldz[31,4] ldz[31,5] ldz[31,6]
 ldz[31,7] ldz[31,8] ldz[31,9] ldz[31,10] ldz[31,11] ldz[31,12] ldz[31,13]
 ldz[31,14] ldz[31,15] ldz[32,0] ldz[32,1] ldz[32,2] ldz[32,3] ldz[32,4]
 ldz[32,5] ldz[32,6] ldz[32,7] ldz[32,8] ldz[32,9] ldz[32,10] ldz[32,11]
 ldz[32,12] ldz[32,13] ldz[32,14] ldz[32,15] lDW[17,0] lDW[17,1] lDW[17,2]
 lDW[17,3] lDW[17,4] lDW[17,5] lDW[17,6] lDW[17,7] lDW[17,8] lDW[17,9]
 lDW[17,10] lDW[17,11] lDW[17,12] lDW[17,13] lDW[17,14] lDW[17,15]
 lDW[18,0] lDW[18,1] lDW[18,2] lDW[18,3] lDW[18,4] lDW[18,5] lDW[18,6]
 lDW[18,7] lDW[18,8] lDW[18,9] lDW[18,10] lDW[18,11] lDW[18,12] lDW[18,13]
 lDW[18,14] lDW[18,15] lDW[19,0] lDW[19,1] lDW[19,2] lDW[19,3] lDW[19,4]
 lDW[19,5] lDW[19,6] lDW[19,7] lDW[19,8] lDW[19,9] lDW[19,10] lDW[19,11]
 lDW[19,12] lDW[19,13] lDW[19,14] lDW[19,15] lDW[20,0] lDW[20,1] lDW[20,2]
 lDW[20,3] lDW[20,4] lDW[20,5] lDW[20,6] lDW[20,7] lDW[20,8] lDW[20,9]
 lDW[20,10] lDW[20,11] lDW[20,12] lDW[20,13] lDW[20,14] lDW[20,15]
 lDW[21,0] lDW[21,1] lDW[21,2] lDW[21,3] lDW[21,4] lDW[21,5] lDW[21,6]
 lDW[21,7] lDW[21,8] lDW[21,9] lDW[21,10] lDW[21,11] lDW[21,12] lDW[21,13]
 lDW[21,14] lDW[21,15] lDW[22,0] lDW[22,1] lDW[22,2] lDW[22,3] lDW[22,4]
 lDW[22,5] lDW[22,6] lDW[22,7] lDW[22,8] lDW[22,9] lDW[22,10] lDW[22,11]
 lDW[22,12] lDW[22,13] lDW[22,14] lDW[22,15] lDW[23,0] lDW[23,1] lDW[23,2]
 lDW[23,3] lDW[23,4] lDW[23,5] lDW[23,6] lDW[23,7] lDW[23,8] lDW[23,9]
 lDW[23,10] lDW[23,11] lDW[23,12] lDW[23,13] lDW[23,14] lDW[23,15]
 lDW[24,0] lDW[24,1] lDW[24,2] lDW[24,3] lDW[24,4] lDW[24,5] lDW[24,6]
 lDW[24,7] lDW[24,8] lDW[24,9] lDW[24,10] lDW[24,11] lDW[24,12] lDW[24,13]
 lDW[24,14] lDW[24,15] lDW[25,0] lDW[25,1] lDW[25,2] lDW[25,3] lDW[25,4]
 lDW[25,5] lDW[25,6] lDW[25,7] lDW[25,8] lDW[25,9] lDW[25,10] lDW[25,11]
 lDW[25,12] lDW[25,13] lDW[25,14] lDW[25,15] lDW[26,0] lDW[26,1] lDW[26,2]
 lDW[26,3] lDW[26,4] lDW[26,5] lDW[26,6] lDW[26,7] lDW[26,8] lDW[26,9]
 lDW[26,10] lDW[26,11] lDW[26,12] lDW[26,13] lDW[26,14] lDW[26,15]
 lDW[27,0] lDW[27,1] lDW[27,2] lDW[27,3] lDW[27,4] lDW[27,5] lDW[27,6]
 lDW[27,7] lDW[27,8] lDW[27,9] lDW[27,10] lDW[27,11] lDW[27,12] lDW[27,13]
 lDW[27,14] lDW[27,15] lDW[28,0] lDW[28,1] lDW[28,2] lDW[28,3] lDW[28,4]
 lDW[28,5] lDW[28,6] lDW[28,7] lDW[28,8] lDW[28,9] lDW[28,10] lDW[28,11]
 lDW[28,12] lDW[28,13] lDW[28,14] lDW[28,15] lDW[29,0] lDW[29,1] lDW[29,2]
 lDW[29,3] lDW[29,4] lDW[29,5] lDW[29,6] lDW[29,7] lDW[29,8] lDW[29,9]
 lDW[29,10] lDW[29,11] lDW[29,12] lDW[29,13] lDW[29,14] lDW[29,15]
 lDW[30,0] lDW[30,1] lDW[30,2] lDW[30,3] lDW[30,4] lDW[30,5] lDW[30,6]
 lDW[30,7] lDW[30,8] lDW[30,9] lDW[30,10] lDW[30,11] lDW[30,12] lDW[30,13]
 lDW[30,14] lDW[30,15] lDW[31,0] lDW[31,1] lDW[31,2] lDW[31,3] lDW[31,4]
 lDW[31,5] lDW[31,6] lDW[31,7] lDW[31,8] lDW[31,9] lDW[31,10] lDW[31,11]
 lDW[31,12] lDW[31,13] lDW[31,14] lDW[31,15] lDW[32,0] lDW[32,1] lDW[32,2]
 lDW[32,3] lDW[32,4] lDW[32,5] lDW[32,6] lDW[32,7] lDW[32,8] lDW[32,9]
 lDW[32,10] lDW[32,11] lDW[32,12] lDW[32,13] lDW[32,14] lDW[32,15]
 ldw[17,0] ldw[17,1] ldw[17,2] ldw[17,3] ldw[17,4] ldw[17,5] ldw[17,6]
 ldw[17,7] ldw[17,8] ldw[17,9] ldw[17,10] ldw[17,11] ldw[17,12] ldw[17,13]
 ldw[17,14] ldw[17,15] ldw[18,0] ldw[18,1] ldw[18,2] ldw[18,3] ldw[18,4]
 ldw[18,5] ldw[18,6] ldw[18,7] ldw[18,8] ldw[18,9] ldw[18,10] ldw[18,11]
 ldw[18,12] ldw[18,13] ldw[18,14] ldw[18,15] ldw[19,0] ldw[19,1] ldw[19,2]
 ldw[19,3] ldw[19,4] ldw[19,5] ldw[19,6] ldw[19,7] ldw[19,8] ldw[19,9]
 ldw[19,10] ldw[19,11] ldw[19,12] ldw[19,13] ldw[19,14] ldw[19,15]
 ldw[20,0] ldw[20,1] ldw[20,2] ldw[20,3] ldw[20,4] ldw[20,5] ldw[20,6]
 ldw[20,7] ldw[20,8] ldw[20,9] ldw[20,10] ldw[20,11] ldw[20,12] ldw[20,13]
 ldw[20,14] ldw[20,15] ldw[21,0] ldw[21,1] ldw[21,2] ldw[21,3] ldw[21,4]
 ldw[21,5] ldw[21,6] ldw[21,7] ldw[21,8] ldw[21,9] ldw[21,10] ldw[21,11]
 ldw[21,12] ldw[21,13] ldw[21,14] ldw[21,15] ldw[22,0] ldw[22,1] ldw[22,2]
 ldw[22,3] ldw[22,4] ldw[22,5] ldw[22,6] ldw[22,7] ldw[22,8] ldw[22,9]
 ldw[22,10] ldw[22,11] ldw[22,12] ldw[22,13] ldw[22,14] ldw[22,15]
 ldw[23,0] ldw[23,1] ldw[23,2] ldw[23,3] ldw[23,4] ldw[23,5] ldw[23,6]
 ldw[23,7] ldw[23,8] ldw[23,9] ldw[23,10] ldw[23,11] ldw[23,12] ldw[23,13]
 ldw[23,14] ldw[23,15] ldw[24,0] ldw[24,1] ldw[24,2] ldw[24,3] ldw[24,4]
 ldw[24,5] ldw[24,6] ldw[24,7] ldw[24,8] ldw[24,9] ldw[24,10] ldw[24,11]
 ldw[24,12] ldw[24,13] ldw[24,14] ldw[24,15] ldw[25,0] ldw[25,1] ldw[25,2]
 ldw[25,3] ldw[25,4] ldw[25,5] ldw[25,6] ldw[25,7] ldw[25,8] ldw[25,9]
 ldw[25,10] ldw[25,11] ldw[25,12] ldw[25,13] ldw[25,14] ldw[25,15]
 ldw[26,0] ldw[26,1] ldw[26,2] ldw[26,3] ldw[26,4] ldw[26,5] ldw[26,6]
 ldw[26,7] ldw[26,8] ldw[26,9] ldw[26,10] ldw[26,11] ldw[26,12] ldw[26,13]
 ldw[26,14] ldw[26,15] ldw[27,0] ldw[27,1] ldw[27,2] ldw[27,3] ldw[27,4]
 ldw[27,5] ldw[27,6] ldw[27,7] ldw[27,8] ldw[27,9] ldw[27,10] ldw[27,11]
 ldw[27,12] ldw[27,13] ldw[27,14] ldw[27,15] ldw[28,0] ldw[28,1] ldw[28,2]
 ldw[28,3] ldw[28,4] ldw[28,5] ldw[28,6] ldw[28,7] ldw[28,8] ldw[28,9]
 ldw[28,10] ldw[28,11] ldw[28,12] ldw[28,13] ldw[28,14] ldw[28,15]
 ldw[29,0] ldw[29,1] ldw[29,2] ldw[29,3] ldw[29,4] ldw[29,5] ldw[29,6]
 ldw[29,7] ldw[29,8] ldw[29,9] ldw[29,10] ldw[29,11] ldw[29,12] ldw[29,13]
 ldw[29,14] ldw[29,15] ldw[30,0] ldw[30,1] ldw[30,2] ldw[30,3] ldw[30,4]
 ldw[30,5] ldw[30,6] ldw[30,7] ldw[30,8] ldw[30,9] ldw[30,10] ldw[30,11]
 ldw[30,12] ldw[30,13] ldw[30,14] ldw[30,15] ldw[31,0] ldw[31,1] ldw[31,2]
 ldw[31,3] ldw[31,4] ldw[31,5] ldw[31,6] ldw[31,7] ldw[31,8] ldw[31,9]
 ldw[31,10] ldw[31,11] ldw[31,12] ldw[31,13] ldw[31,14] ldw[31,15]
 ldw[32,0] ldw[32,1] ldw[32,2] ldw[32,3] ldw[32,4] ldw[32,5] ldw[32,6]
 ldw[32,7] ldw[32,8] ldw[32,9] ldw[32,10] ldw[32,11] ldw[32,12] ldw[32,13]
 ldw[32,14] ldw[32,15] lstk[17,0] lstk[17,1] lstk[17,2] lstk[17,3]
 lstk[17,4] lstk[17,5] lstk[17,6] lstk[17,7] lstk[17,8] lstk[17,9]
 lstk[17,10] lstk[17,11] lstk[17,12] lstk[17,13] lstk[17,14] lstk[17,15]
 lstk[18,0] lstk[18,1] lstk[18,2] lstk[18,3] lstk[18,4] lstk[18,5]
 lstk[18,6] lstk[18,7] lstk[18,8] lstk[18,9] lstk[18,10] lstk[18,11]
 lstk[18,12] lstk[18,13] lstk[18,14] lstk[18,15] lstk[19,0] lstk[19,1]
 lstk[19,2] lstk[19,3] lstk[19,4] lstk[19,5] lstk[19,6] lstk[19,7]
 lstk[19,8] lstk[19,9] lstk[19,10] lstk[19,11] lstk[19,12] lstk[19,13]
 lstk[19,14] lstk[19,15] lstk[20,0] lstk[20,1] lstk[20,2] lstk[20,3]
 lstk[20,4] lstk[20,5] lstk[20,6] lstk[20,7] lstk[20,8] lstk[20,9]
 lstk[20,10] lstk[20,11] lstk[20,12] lstk[20,13] lstk[20,14] lstk[20,15]
 lstk[21,0] lstk[21,1] lstk[21,2] lstk[21,3] lstk[21,4] lstk[21,5]
 lstk[21,6] lstk[21,7] lstk[21,8] lstk[21,9] lstk[21,10] lstk[21,11]
 lstk[21,12] lstk[21,13] lstk[21,14] lstk[21,15] lstk[22,0] lstk[22,1]
 lstk[22,2] lstk[22,3] lstk[22,4] lstk[22,5] lstk[22,6] lstk[22,7]
 lstk[22,8] lstk[22,9] lstk[22,10] lstk[22,11] lstk[22,12] lstk[22,13]
 lstk[22,14] lstk[22,15] lstk[23,0] lstk[23,1] lstk[23,2] lstk[23,3]
 lstk[23,4] lstk[23,5] lstk[23,6] lstk[23,7] lstk[23,8] lstk[23,9]
 lstk[23,10] lstk[23,11] lstk[23,12] lstk[23,13] lstk[23,14] lstk[23,15]
 lstk[24,0] lstk[24,1] lstk[24,2] lstk[24,3] lstk[24,4] lstk[24,5]
 lstk[24,6] lstk[24,7] lstk[24,8] lstk[24,9] lstk[24,10] lstk[24,11]
 lstk[24,12] lstk[24,13] lstk[24,14] lstk[24,15] lstk[25,0] lstk[25,1]
 lstk[25,2] lstk[25,3] lstk[25,4] lstk[25,5] lstk[25,6] lstk[25,7]
 lstk[25,8] lstk[25,9] lstk[25,10] lstk[25,11] lstk[25,12] lstk[25,13]
 lstk[25,14] lstk[25,15] lstk[26,0] lstk[26,1] lstk[26,2] lstk[26,3]
 lstk[26,4] lstk[26,5] lstk[26,6] lstk[26,7] lstk[26,8] lstk[26,9]
 lstk[26,10] lstk[26,11] lstk[26,12] lstk[26,13] lstk[26,14] lstk[26,15]
 lstk[27,0] lstk[27,1] lstk[27,2] lstk[27,3] lstk[27,4] lstk[27,5]
 lstk[27,6] lstk[27,7] lstk[27,8] lstk[27,9] lstk[27,10] lstk[27,11]
 lstk[27,12] lstk[27,13] lstk[27,14] lstk[27,15] lstk[28,0] lstk[28,1]
 lstk[28,2] lstk[28,3] lstk[28,4] lstk[28,5] lstk[28,6] lstk[28,7]
 lstk[28,8] lstk[28,9] lstk[28,10] lstk[28,11] lstk[28,12] lstk[28,13]
 lstk[28,14] lstk[28,15] lstk[29,0] lstk[29,1] lstk[29,2] lstk[29,3]
 lstk[29,4] lstk[29,5] lstk[29,6] lstk[29,7] lstk[29,8] lstk[29,9]
 lstk[29,10] lstk[29,11] lstk[29,12] lstk[29,13] lstk[29,14] lstk[29,15]
 lstk[30,0] lstk[30,1] lstk[30,2] lstk[30,3] lstk[30,4] lstk[30,5]
 lstk[30,6] lstk[30,7] lstk[30,8] lstk[30,9] lstk[30,10] lstk[30,11]
 lstk[30,12] lstk[30,13] lstk[30,14] lstk[30,15] lstk[31,0] lstk[31,1]
 lstk[31,2] lstk[31,3] lstk[31,4] lstk[31,5] lstk[31,6] lstk[31,7]
 lstk[31,8] lstk[31,9] lstk[31,10] lstk[31,11] lstk[31,12] lstk[31,13]
 lstk[31,14] lstk[31,15] lstk[32,0] lstk[32,1] lstk[32,2] lstk[32,3]
 lstk[32,4] lstk[32,5] lstk[32,6] lstk[32,7] lstk[32,8] lstk[32,9]
 lstk[32,10] lstk[32,11] lstk[32,12] lstk[32,13] lstk[32,14] lstk[32,15]
 lLANE[0] lLANE[1] lLANE[2] lLANE[3] lLANE[4] lLANE[5] lLANE[6] lLANE[7]
 lLANE[8] lLANE[9] lLANE[10] lLANE[11] lLANE[12] lLANE[13] lLANE[14]
 lLANE[15] ugEQK[0,0] ugEQK[0,1] ugEQK[0,2] ugEQK[0,3] ugEQK[0,4]
 ugEQK[0,5] ugEQK[0,6] ugEQK[0,7] ugEQK[0,8] ugEQK[0,9] ugEQK[0,10]
 ugEQK[0,11] ugEQK[0,12] ugEQK[0,13] ugEQK[0,14] ugEQK[0,15] ugEQK[1,0]
 ugEQK[1,1] ugEQK[1,2] ugEQK[1,3] ugEQK[1,4] ugEQK[1,5] ugEQK[1,6]
 ugEQK[1,7] ugEQK[1,8] ugEQK[1,9] ugEQK[1,10] ugEQK[1,11] ugEQK[1,12]
 ugEQK[1,13] ugEQK[1,14] ugEQK[1,15] ugEQK[2,0] ugEQK[2,1] ugEQK[2,2]
 ugEQK[2,3] ugEQK[2,4] ugEQK[2,5] ugEQK[2,6] ugEQK[2,7] ugEQK[2,8]
 ugEQK[2,9] ugEQK[2,10] ugEQK[2,11] ugEQK[2,12] ugEQK[2,13] ugEQK[2,14]
 ugEQK[2,15] ugEQK[3,0] ugEQK[3,1] ugEQK[3,2] ugEQK[3,3] ugEQK[3,4]
 ugEQK[3,5] ugEQK[3,6] ugEQK[3,7] ugEQK[3,8] ugEQK[3,9] ugEQK[3,10]
 ugEQK[3,11] ugEQK[3,12] ugEQK[3,13] ugEQK[3,14] ugEQK[3,15] ugEQK[4,0]
 ugEQK[4,1] ugEQK[4,2] ugEQK[4,3] ugEQK[4,4] ugEQK[4,5] ugEQK[4,6]
 ugEQK[4,7] ugEQK[4,8] ugEQK[4,9] ugEQK[4,10] ugEQK[4,11] ugEQK[4,12]
 ugEQK[4,13] ugEQK[4,14] ugEQK[4,15] ugEQK[5,0] ugEQK[5,1] ugEQK[5,2]
 ugEQK[5,3] ugEQK[5,4] ugEQK[5,5] ugEQK[5,6] ugEQK[5,7] ugEQK[5,8]
 ugEQK[5,9] ugEQK[5,10] ugEQK[5,11] ugEQK[5,12] ugEQK[5,13] ugEQK[5,14]
 ugEQK[5,15] uGstk[0,0] uGstk[0,1] uGstk[0,2] uGstk[0,3] uGstk[0,4]
 uGstk[0,5] uGstk[0,6] uGstk[0,7] uGstk[0,8] uGstk[0,9] uGstk[0,10]
 uGstk[0,11] uGstk[0,12] uGstk[0,13] uGstk[0,14] uGstk[0,15] uGstk[1,0]
 uGstk[1,1] uGstk[1,2] uGstk[1,3] uGstk[1,4] uGstk[1,5] uGstk[1,6]
 uGstk[1,7] uGstk[1,8] uGstk[1,9] uGstk[1,10] uGstk[1,11] uGstk[1,12]
 uGstk[1,13] uGstk[1,14] uGstk[1,15] uGstk[2,0] uGstk[2,1] uGstk[2,2]
 uGstk[2,3] uGstk[2,4] uGstk[2,5] uGstk[2,6] uGstk[2,7] uGstk[2,8]
 uGstk[2,9] uGstk[2,10] uGstk[2,11] uGstk[2,12] uGstk[2,13] uGstk[2,14]
 uGstk[2,15] uGstk[3,0] uGstk[3,1] uGstk[3,2] uGstk[3,3] uGstk[3,4]
 uGstk[3,5] uGstk[3,6] uGstk[3,7] uGstk[3,8] uGstk[3,9] uGstk[3,10]
 uGstk[3,11] uGstk[3,12] uGstk[3,13] uGstk[3,14] uGstk[3,15] uGstk[4,0]
 uGstk[4,1] uGstk[4,2] uGstk[4,3] uGstk[4,4] uGstk[4,5] uGstk[4,6]
 uGstk[4,7] uGstk[4,8] uGstk[4,9] uGstk[4,10] uGstk[4,11] uGstk[4,12]
 uGstk[4,13] uGstk[4,14] uGstk[4,15] uGstk[5,0] uGstk[5,1] uGstk[5,2]
 uGstk[5,3] uGstk[5,4] uGstk[5,5] uGstk[5,6] uGstk[5,7] uGstk[5,8]
 uGstk[5,9] uGstk[5,10] uGstk[5,11] uGstk[5,12] uGstk[5,13] uGstk[5,14]
 uGstk[5,15] uDetX[1,0] uDetX[1,1] uDetX[1,2] uDetX[1,3] uDetX[1,4]
 uDetX[1,5] uDetX[1,6] uDetX[1,7] uDetX[1,8] uDetX[1,9] uDetX[1,10]
 uDetX[1,11] uDetX[1,12] uDetX[1,13] uDetX[1,14] uDetX[1,15] uDetX[2,0]
 uDetX[2,1] uDetX[2,2] uDetX[2,3] uDetX[2,4] uDetX[2,5] uDetX[2,6]
 uDetX[2,7] uDetX[2,8] uDetX[2,9] uDetX[2,10] uDetX[2,11] uDetX[2,12]
 uDetX[2,13] uDetX[2,14] uDetX[2,15] uDetX[3,0] uDetX[3,1] uDetX[3,2]
 uDetX[3,3] uDetX[3,4] uDetX[3,5] uDetX[3,6] uDetX[3,7] uDetX[3,8]
 uDetX[3,9] uDetX[3,10] uDetX[3,11] uDetX[3,12] uDetX[3,13] uDetX[3,14]
 uDetX[3,15] uDetX[4,0] uDetX[4,1] uDetX[4,2] uDetX[4,3] uDetX[4,4]
 uDetX[4,5] uDetX[4,6] uDetX[4,7] uDetX[4,8] uDetX[4,9] uDetX[4,10]
 uDetX[4,11] uDetX[4,12] uDetX[4,13] uDetX[4,14] uDetX[4,15] uDetX[5,0]
 uDetX[5,1] uDetX[5,2] uDetX[5,3] uDetX[5,4] uDetX[5,5] uDetX[5,6]
 uDetX[5,7] uDetX[5,8] uDetX[5,9] uDetX[5,10] uDetX[5,11] uDetX[5,12]
 uDetX[5,13] uDetX[5,14] uDetX[5,15] uDetX[6,0] uDetX[6,1] uDetX[6,2]
 uDetX[6,3] uDetX[6,4] uDetX[6,5] uDetX[6,6] uDetX[6,7] uDetX[6,8]
 uDetX[6,9] uDetX[6,10] uDetX[6,11] uDetX[6,12] uDetX[6,13] uDetX[6,14]
 uDetX[6,15] uDetY[1,0] uDetY[1,1] uDetY[1,2] uDetY[1,3] uDetY[1,4]
 uDetY[1,5] uDetY[1,6] uDetY[1,7] uDetY[1,8] uDetY[1,9] uDetY[1,10]
 uDetY[1,11] uDetY[1,12] uDetY[1,13] uDetY[1,14] uDetY[1,15] uDetY[2,0]
 uDetY[2,1] uDetY[2,2] uDetY[2,3] uDetY[2,4] uDetY[2,5] uDetY[2,6]
 uDetY[2,7] uDetY[2,8] uDetY[2,9] uDetY[2,10] uDetY[2,11] uDetY[2,12]
 uDetY[2,13] uDetY[2,14] uDetY[2,15] uDetY[3,0] uDetY[3,1] uDetY[3,2]
 uDetY[3,3] uDetY[3,4] uDetY[3,5] uDetY[3,6] uDetY[3,7] uDetY[3,8]
 uDetY[3,9] uDetY[3,10] uDetY[3,11] uDetY[3,12] uDetY[3,13] uDetY[3,14]
 uDetY[3,15] uDetY[4,0] uDetY[4,1] uDetY[4,2] uDetY[4,3] uDetY[4,4]
 uDetY[4,5] uDetY[4,6] uDetY[4,7] uDetY[4,8] uDetY[4,9] uDetY[4,10]
 uDetY[4,11] uDetY[4,12] uDetY[4,13] uDetY[4,14] uDetY[4,15] uDetY[5,0]
 uDetY[5,1] uDetY[5,2] uDetY[5,3] uDetY[5,4] uDetY[5,5] uDetY[5,6]
 uDetY[5,7] uDetY[5,8] uDetY[5,9] uDetY[5,10] uDetY[5,11] uDetY[5,12]
 uDetY[5,13] uDetY[5,14] uDetY[5,15] uDetY[6,0] uDetY[6,1] uDetY[6,2]
 uDetY[6,3] uDetY[6,4] uDetY[6,5] uDetY[6,6] uDetY[6,7] uDetY[6,8]
 uDetY[6,9] uDetY[6,10] uDetY[6,11] uDetY[6,12] uDetY[6,13] uDetY[6,14]
 uDetY[6,15] uDetZ[1,0] uDetZ[1,1] uDetZ[1,2] uDetZ[1,3] uDetZ[1,4]
 uDetZ[1,5] uDetZ[1,6] uDetZ[1,7] uDetZ[1,8] uDetZ[1,9] uDetZ[1,10]
 uDetZ[1,11] uDetZ[1,12] uDetZ[1,13] uDetZ[1,14] uDetZ[1,15] uDetZ[2,0]
 uDetZ[2,1] uDetZ[2,2] uDetZ[2,3] uDetZ[2,4] uDetZ[2,5] uDetZ[2,6]
 uDetZ[2,7] uDetZ[2,8] uDetZ[2,9] uDetZ[2,10] uDetZ[2,11] uDetZ[2,12]
 uDetZ[2,13] uDetZ[2,14] uDetZ[2,15] uDetZ[3,0] uDetZ[3,1] uDetZ[3,2]
 uDetZ[3,3] uDetZ[3,4] uDetZ[3,5] uDetZ[3,6] uDetZ[3,7] uDetZ[3,8]
 uDetZ[3,9] uDetZ[3,10] uDetZ[3,11] uDetZ[3,12] uDetZ[3,13] uDetZ[3,14]
 uDetZ[3,15] uDetZ[4,0] uDetZ[4,1] uDetZ[4,2] uDetZ[4,3] uDetZ[4,4]
 uDetZ[4,5] uDetZ[4,6] uDetZ[4,7] uDetZ[4,8] uDetZ[4,9] uDetZ[4,10]
 uDetZ[4,11] uDetZ[4,12] uDetZ[4,13] uDetZ[4,14] uDetZ[4,15] uDetZ[5,0]
 uDetZ[5,1] uDetZ[5,2] uDetZ[5,3] uDetZ[5,4] uDetZ[5,5] uDetZ[5,6]
 uDetZ[5,7] uDetZ[5,8] uDetZ[5,9] uDetZ[5,10] uDetZ[5,11] uDetZ[5,12]
 uDetZ[5,13] uDetZ[5,14] uDetZ[5,15] uDetW[0,0] uDetW[0,1] uDetW[0,2]
 uDetW[0,3] uDetW[0,4] uDetW[0,5] uDetW[0,6] uDetW[0,7] uDetW[0,8]
 uDetW[0,9] uDetW[0,10] uDetW[0,11] uDetW[0,12] uDetW[0,13] uDetW[0,14]
 uDetW[0,15] uDetW[1,0] uDetW[1,1] uDetW[1,2] uDetW[1,3] uDetW[1,4]
 uDetW[1,5] uDetW[1,6] uDetW[1,7] uDetW[1,8] uDetW[1,9] uDetW[1,10]
 uDetW[1,11] uDetW[1,12] uDetW[1,13] uDetW[1,14] uDetW[1,15] uDetW[2,0]
 uDetW[2,1] uDetW[2,2] uDetW[2,3] uDetW[2,4] uDetW[2,5] uDetW[2,6]
 uDetW[2,7] uDetW[2,8] uDetW[2,9] uDetW[2,10] uDetW[2,11] uDetW[2,12]
 uDetW[2,13] uDetW[2,14] uDetW[2,15] uDetW[3,0] uDetW[3,1] uDetW[3,2]
 uDetW[3,3] uDetW[3,4] uDetW[3,5] uDetW[3,6] uDetW[3,7] uDetW[3,8]
 uDetW[3,9] uDetW[3,10] uDetW[3,11] uDetW[3,12] uDetW[3,13] uDetW[3,14]
 uDetW[3,15] uDetW[4,0] uDetW[4,1] uDetW[4,2] uDetW[4,3] uDetW[4,4]
 uDetW[4,5] uDetW[4,6] uDetW[4,7] uDetW[4,8] uDetW[4,9] uDetW[4,10]
 uDetW[4,11] uDetW[4,12] uDetW[4,13] uDetW[4,14] uDetW[4,15] uDetW[5,0]
 uDetW[5,1] uDetW[5,2] uDetW[5,3] uDetW[5,4] uDetW[5,5] uDetW[5,6]
 uDetW[5,7] uDetW[5,8] uDetW[5,9] uDetW[5,10] uDetW[5,11] uDetW[5,12]
 uDetW[5,13] uDetW[5,14] uDetW[5,15] uFrSC[1,0] uFrSC[1,1] uFrSC[1,2]
 uFrSC[1,3] uFrSC[1,4] uFrSC[1,5] uFrSC[1,6] uFrSC[1,7] uFrSC[1,8]
 uFrSC[1,9] uFrSC[1,10] uFrSC[1,11] uFrSC[1,12] uFrSC[1,13] uFrSC[1,14]
 uFrSC[1,15] uFrSC[2,0] uFrSC[2,1] uFrSC[2,2] uFrSC[2,3] uFrSC[2,4]
 uFrSC[2,5] uFrSC[2,6] uFrSC[2,7] uFrSC[2,8] uFrSC[2,9] uFrSC[2,10]
 uFrSC[2,11] uFrSC[2,12] uFrSC[2,13] uFrSC[2,14] uFrSC[2,15] uFrSC[3,0]
 uFrSC[3,1] uFrSC[3,2] uFrSC[3,3] uFrSC[3,4] uFrSC[3,5] uFrSC[3,6]
 uFrSC[3,7] uFrSC[3,8] uFrSC[3,9] uFrSC[3,10] uFrSC[3,11] uFrSC[3,12]
 uFrSC[3,13] uFrSC[3,14] uFrSC[3,15] uFrSC[4,0] uFrSC[4,1] uFrSC[4,2]
 uFrSC[4,3] uFrSC[4,4] uFrSC[4,5] uFrSC[4,6] uFrSC[4,7] uFrSC[4,8]
 uFrSC[4,9] uFrSC[4,10] uFrSC[4,11] uFrSC[4,12] uFrSC[4,13] uFrSC[4,14]
 uFrSC[4,15] uFrSC[5,0] uFrSC[5,1] uFrSC[5,2] uFrSC[5,3] uFrSC[5,4]
 uFrSC[5,5] uFrSC[5,6] uFrSC[5,7] uFrSC[5,8] uFrSC[5,9] uFrSC[5,10]
 uFrSC[5,11] uFrSC[5,12] uFrSC[5,13] uFrSC[5,14] uFrSC[5,15] uFrSC[6,0]
 uFrSC[6,1] uFrSC[6,2] uFrSC[6,3] uFrSC[6,4] uFrSC[6,5] uFrSC[6,6]
 uFrSC[6,7] uFrSC[6,8] uFrSC[6,9] uFrSC[6,10] uFrSC[6,11] uFrSC[6,12]
 uFrSC[6,13] uFrSC[6,14] uFrSC[6,15] uFrMC[2,0] uFrMC[2,1] uFrMC[2,2]
 uFrMC[2,3] uFrMC[2,4] uFrMC[2,5] uFrMC[2,6] uFrMC[2,7] uFrMC[2,8]
 uFrMC[2,9] uFrMC[2,10] uFrMC[2,11] uFrMC[2,12] uFrMC[2,13] uFrMC[2,14]
 uFrMC[2,15] uFrMC[3,0] uFrMC[3,1] uFrMC[3,2] uFrMC[3,3] uFrMC[3,4]
 uFrMC[3,5] uFrMC[3,6] uFrMC[3,7] uFrMC[3,8] uFrMC[3,9] uFrMC[3,10]
 uFrMC[3,11] uFrMC[3,12] uFrMC[3,13] uFrMC[3,14] uFrMC[3,15] uFrMC[4,0]
 uFrMC[4,1] uFrMC[4,2] uFrMC[4,3] uFrMC[4,4] uFrMC[4,5] uFrMC[4,6]
 uFrMC[4,7] uFrMC[4,8] uFrMC[4,9] uFrMC[4,10] uFrMC[4,11] uFrMC[4,12]
 uFrMC[4,13] uFrMC[4,14] uFrMC[4,15] uFrMC[5,0] uFrMC[5,1] uFrMC[5,2]
 uFrMC[5,3] uFrMC[5,4] uFrMC[5,5] uFrMC[5,6] uFrMC[5,7] uFrMC[5,8]
 uFrMC[5,9] uFrMC[5,10] uFrMC[5,11] uFrMC[5,12] uFrMC[5,13] uFrMC[5,14]
 uFrMC[5,15] uFrMC[6,0] uFrMC[6,1] uFrMC[6,2] uFrMC[6,3] uFrMC[6,4]
 uFrMC[6,5] uFrMC[6,6] uFrMC[6,7] uFrMC[6,8] uFrMC[6,9] uFrMC[6,10]
 uFrMC[6,11] uFrMC[6,12] uFrMC[6,13] uFrMC[6,14] uFrMC[6,15] ugeqka[0,0]
 ugeqka[0,1] ugeqka[0,2] ugeqka[0,3] ugeqka[0,4] ugeqka[0,5] ugeqka[0,6]
 ugeqka[0,7] ugeqka[0,8] ugeqka[0,9] ugeqka[0,10] ugeqka[0,11] ugeqka[0,12]
 ugeqka[0,13] ugeqka[0,14] ugeqka[0,15] ugeqka[1,0] ugeqka[1,1] ugeqka[1,2]
 ugeqka[1,3] ugeqka[1,4] ugeqka[1,5] ugeqka[1,6] ugeqka[1,7] ugeqka[1,8]
 ugeqka[1,9] ugeqka[1,10] ugeqka[1,11] ugeqka[1,12] ugeqka[1,13]
 ugeqka[1,14] ugeqka[1,15] ugeqka[2,0] ugeqka[2,1] ugeqka[2,2] ugeqka[2,3]
 ugeqka[2,4] ugeqka[2,5] ugeqka[2,6] ugeqka[2,7] ugeqka[2,8] ugeqka[2,9]
 ugeqka[2,10] ugeqka[2,11] ugeqka[2,12] ugeqka[2,13] ugeqka[2,14]
 ugeqka[2,15] ugeqka[3,0] ugeqka[3,1] ugeqka[3,2] ugeqka[3,3] ugeqka[3,4]
 ugeqka[3,5] ugeqka[3,6] ugeqka[3,7] ugeqka[3,8] ugeqka[3,9] ugeqka[3,10]
 ugeqka[3,11] ugeqka[3,12] ugeqka[3,13] ugeqka[3,14] ugeqka[3,15]
 ugeqka[4,0] ugeqka[4,1] ugeqka[4,2] ugeqka[4,3] ugeqka[4,4] ugeqka[4,5]
 ugeqka[4,6] ugeqka[4,7] ugeqka[4,8] ugeqka[4,9] ugeqka[4,10] ugeqka[4,11]
 ugeqka[4,12] ugeqka[4,13] ugeqka[4,14] ugeqka[4,15] ugeqka[5,0]
 ugeqka[5,1] ugeqka[5,2] ugeqka[5,3] ugeqka[5,4] ugeqka[5,5] ugeqka[5,6]
 ugeqka[5,7] ugeqka[5,8] ugeqka[5,9] ugeqka[5,10] ugeqka[5,11] ugeqka[5,12]
 ugeqka[5,13] ugeqka[5,14] ugeqka[5,15] ugstka[0,0] ugstka[0,1] ugstka[0,2]
 ugstka[0,3] ugstka[0,4] ugstka[0,5] ugstka[0,6] ugstka[0,7] ugstka[0,8]
 ugstka[0,9] ugstka[0,10] ugstka[0,11] ugstka[0,12] ugstka[0,13]
 ugstka[0,14] ugstka[0,15] ugstka[1,0] ugstka[1,1] ugstka[1,2] ugstka[1,3]
 ugstka[1,4] ugstka[1,5] ugstka[1,6] ugstka[1,7] ugstka[1,8] ugstka[1,9]
 ugstka[1,10] ugstka[1,11] ugstka[1,12] ugstka[1,13] ugstka[1,14]
 ugstka[1,15] ugstka[2,0] ugstka[2,1] ugstka[2,2] ugstka[2,3] ugstka[2,4]
 ugstka[2,5] ugstka[2,6] ugstka[2,7] ugstka[2,8] ugstka[2,9] ugstka[2,10]
 ugstka[2,11] ugstka[2,12] ugstka[2,13] ugstka[2,14] ugstka[2,15]
 ugstka[3,0] ugstka[3,1] ugstka[3,2] ugstka[3,3] ugstka[3,4] ugstka[3,5]
 ugstka[3,6] ugstka[3,7] ugstka[3,8] ugstka[3,9] ugstka[3,10] ugstka[3,11]
 ugstka[3,12] ugstka[3,13] ugstka[3,14] ugstka[3,15] ugstka[4,0]
 ugstka[4,1] ugstka[4,2] ugstka[4,3] ugstka[4,4] ugstka[4,5] ugstka[4,6]
 ugstka[4,7] ugstka[4,8] ugstka[4,9] ugstka[4,10] ugstka[4,11] ugstka[4,12]
 ugstka[4,13] ugstka[4,14] ugstka[4,15] ugstka[5,0] ugstka[5,1] ugstka[5,2]
 ugstka[5,3] ugstka[5,4] ugstka[5,5] ugstka[5,6] ugstka[5,7] ugstka[5,8]
 ugstka[5,9] ugstka[5,10] ugstka[5,11] ugstka[5,12] ugstka[5,13]
 ugstka[5,14] ugstka[5,15] udetXa[1,0] udetXa[1,1] udetXa[1,2] udetXa[1,3]
 udetXa[1,4] udetXa[1,5] udetXa[1,6] udetXa[1,7] udetXa[1,8] udetXa[1,9]
 udetXa[1,10] udetXa[1,11] udetXa[1,12] udetXa[1,13] udetXa[1,14]
 udetXa[1,15] udetXa[2,0] udetXa[2,1] udetXa[2,2] udetXa[2,3] udetXa[2,4]
 udetXa[2,5] udetXa[2,6] udetXa[2,7] udetXa[2,8] udetXa[2,9] udetXa[2,10]
 udetXa[2,11] udetXa[2,12] udetXa[2,13] udetXa[2,14] udetXa[2,15]
 udetXa[3,0] udetXa[3,1] udetXa[3,2] udetXa[3,3] udetXa[3,4] udetXa[3,5]
 udetXa[3,6] udetXa[3,7] udetXa[3,8] udetXa[3,9] udetXa[3,10] udetXa[3,11]
 udetXa[3,12] udetXa[3,13] udetXa[3,14] udetXa[3,15] udetXa[4,0]
 udetXa[4,1] udetXa[4,2] udetXa[4,3] udetXa[4,4] udetXa[4,5] udetXa[4,6]
 udetXa[4,7] udetXa[4,8] udetXa[4,9] udetXa[4,10] udetXa[4,11] udetXa[4,12]
 udetXa[4,13] udetXa[4,14] udetXa[4,15] udetXa[5,0] udetXa[5,1] udetXa[5,2]
 udetXa[5,3] udetXa[5,4] udetXa[5,5] udetXa[5,6] udetXa[5,7] udetXa[5,8]
 udetXa[5,9] udetXa[5,10] udetXa[5,11] udetXa[5,12] udetXa[5,13]
 udetXa[5,14] udetXa[5,15] udetXa[6,0] udetXa[6,1] udetXa[6,2] udetXa[6,3]
 udetXa[6,4] udetXa[6,5] udetXa[6,6] udetXa[6,7] udetXa[6,8] udetXa[6,9]
 udetXa[6,10] udetXa[6,11] udetXa[6,12] udetXa[6,13] udetXa[6,14]
 udetXa[6,15] udetYa[1,0] udetYa[1,1] udetYa[1,2] udetYa[1,3] udetYa[1,4]
 udetYa[1,5] udetYa[1,6] udetYa[1,7] udetYa[1,8] udetYa[1,9] udetYa[1,10]
 udetYa[1,11] udetYa[1,12] udetYa[1,13] udetYa[1,14] udetYa[1,15]
 udetYa[2,0] udetYa[2,1] udetYa[2,2] udetYa[2,3] udetYa[2,4] udetYa[2,5]
 udetYa[2,6] udetYa[2,7] udetYa[2,8] udetYa[2,9] udetYa[2,10] udetYa[2,11]
 udetYa[2,12] udetYa[2,13] udetYa[2,14] udetYa[2,15] udetYa[3,0]
 udetYa[3,1] udetYa[3,2] udetYa[3,3] udetYa[3,4] udetYa[3,5] udetYa[3,6]
 udetYa[3,7] udetYa[3,8] udetYa[3,9] udetYa[3,10] udetYa[3,11] udetYa[3,12]
 udetYa[3,13] udetYa[3,14] udetYa[3,15] udetYa[4,0] udetYa[4,1] udetYa[4,2]
 udetYa[4,3] udetYa[4,4] udetYa[4,5] udetYa[4,6] udetYa[4,7] udetYa[4,8]
 udetYa[4,9] udetYa[4,10] udetYa[4,11] udetYa[4,12] udetYa[4,13]
 udetYa[4,14] udetYa[4,15] udetYa[5,0] udetYa[5,1] udetYa[5,2] udetYa[5,3]
 udetYa[5,4] udetYa[5,5] udetYa[5,6] udetYa[5,7] udetYa[5,8] udetYa[5,9]
 udetYa[5,10] udetYa[5,11] udetYa[5,12] udetYa[5,13] udetYa[5,14]
 udetYa[5,15] udetYa[6,0] udetYa[6,1] udetYa[6,2] udetYa[6,3] udetYa[6,4]
 udetYa[6,5] udetYa[6,6] udetYa[6,7] udetYa[6,8] udetYa[6,9] udetYa[6,10]
 udetYa[6,11] udetYa[6,12] udetYa[6,13] udetYa[6,14] udetYa[6,15]
 udetZa[1,0] udetZa[1,1] udetZa[1,2] udetZa[1,3] udetZa[1,4] udetZa[1,5]
 udetZa[1,6] udetZa[1,7] udetZa[1,8] udetZa[1,9] udetZa[1,10] udetZa[1,11]
 udetZa[1,12] udetZa[1,13] udetZa[1,14] udetZa[1,15] udetZa[2,0]
 udetZa[2,1] udetZa[2,2] udetZa[2,3] udetZa[2,4] udetZa[2,5] udetZa[2,6]
 udetZa[2,7] udetZa[2,8] udetZa[2,9] udetZa[2,10] udetZa[2,11] udetZa[2,12]
 udetZa[2,13] udetZa[2,14] udetZa[2,15] udetZa[3,0] udetZa[3,1] udetZa[3,2]
 udetZa[3,3] udetZa[3,4] udetZa[3,5] udetZa[3,6] udetZa[3,7] udetZa[3,8]
 udetZa[3,9] udetZa[3,10] udetZa[3,11] udetZa[3,12] udetZa[3,13]
 udetZa[3,14] udetZa[3,15] udetZa[4,0] udetZa[4,1] udetZa[4,2] udetZa[4,3]
 udetZa[4,4] udetZa[4,5] udetZa[4,6] udetZa[4,7] udetZa[4,8] udetZa[4,9]
 udetZa[4,10] udetZa[4,11] udetZa[4,12] udetZa[4,13] udetZa[4,14]
 udetZa[4,15] udetZa[5,0] udetZa[5,1] udetZa[5,2] udetZa[5,3] udetZa[5,4]
 udetZa[5,5] udetZa[5,6] udetZa[5,7] udetZa[5,8] udetZa[5,9] udetZa[5,10]
 udetZa[5,11] udetZa[5,12] udetZa[5,13] udetZa[5,14] udetZa[5,15]
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
 udetWa[3,15] udetWa[4,0] udetWa[4,1] udetWa[4,2] udetWa[4,3] udetWa[4,4]
 udetWa[4,5] udetWa[4,6] udetWa[4,7] udetWa[4,8] udetWa[4,9] udetWa[4,10]
 udetWa[4,11] udetWa[4,12] udetWa[4,13] udetWa[4,14] udetWa[4,15]
 udetWa[5,0] udetWa[5,1] udetWa[5,2] udetWa[5,3] udetWa[5,4] udetWa[5,5]
 udetWa[5,6] udetWa[5,7] udetWa[5,8] udetWa[5,9] udetWa[5,10] udetWa[5,11]
 udetWa[5,12] udetWa[5,13] udetWa[5,14] udetWa[5,15] ufrSCa[1,0]
 ufrSCa[1,1] ufrSCa[1,2] ufrSCa[1,3] ufrSCa[1,4] ufrSCa[1,5] ufrSCa[1,6]
 ufrSCa[1,7] ufrSCa[1,8] ufrSCa[1,9] ufrSCa[1,10] ufrSCa[1,11] ufrSCa[1,12]
 ufrSCa[1,13] ufrSCa[1,14] ufrSCa[1,15] ufrSCa[2,0] ufrSCa[2,1] ufrSCa[2,2]
 ufrSCa[2,3] ufrSCa[2,4] ufrSCa[2,5] ufrSCa[2,6] ufrSCa[2,7] ufrSCa[2,8]
 ufrSCa[2,9] ufrSCa[2,10] ufrSCa[2,11] ufrSCa[2,12] ufrSCa[2,13]
 ufrSCa[2,14] ufrSCa[2,15] ufrSCa[3,0] ufrSCa[3,1] ufrSCa[3,2] ufrSCa[3,3]
 ufrSCa[3,4] ufrSCa[3,5] ufrSCa[3,6] ufrSCa[3,7] ufrSCa[3,8] ufrSCa[3,9]
 ufrSCa[3,10] ufrSCa[3,11] ufrSCa[3,12] ufrSCa[3,13] ufrSCa[3,14]
 ufrSCa[3,15] ufrSCa[4,0] ufrSCa[4,1] ufrSCa[4,2] ufrSCa[4,3] ufrSCa[4,4]
 ufrSCa[4,5] ufrSCa[4,6] ufrSCa[4,7] ufrSCa[4,8] ufrSCa[4,9] ufrSCa[4,10]
 ufrSCa[4,11] ufrSCa[4,12] ufrSCa[4,13] ufrSCa[4,14] ufrSCa[4,15]
 ufrSCa[5,0] ufrSCa[5,1] ufrSCa[5,2] ufrSCa[5,3] ufrSCa[5,4] ufrSCa[5,5]
 ufrSCa[5,6] ufrSCa[5,7] ufrSCa[5,8] ufrSCa[5,9] ufrSCa[5,10] ufrSCa[5,11]
 ufrSCa[5,12] ufrSCa[5,13] ufrSCa[5,14] ufrSCa[5,15] ufrSCa[6,0]
 ufrSCa[6,1] ufrSCa[6,2] ufrSCa[6,3] ufrSCa[6,4] ufrSCa[6,5] ufrSCa[6,6]
 ufrSCa[6,7] ufrSCa[6,8] ufrSCa[6,9] ufrSCa[6,10] ufrSCa[6,11] ufrSCa[6,12]
 ufrSCa[6,13] ufrSCa[6,14] ufrSCa[6,15] ufrMCa[2,0] ufrMCa[2,1] ufrMCa[2,2]
 ufrMCa[2,3] ufrMCa[2,4] ufrMCa[2,5] ufrMCa[2,6] ufrMCa[2,7] ufrMCa[2,8]
 ufrMCa[2,9] ufrMCa[2,10] ufrMCa[2,11] ufrMCa[2,12] ufrMCa[2,13]
 ufrMCa[2,14] ufrMCa[2,15] ufrMCa[3,0] ufrMCa[3,1] ufrMCa[3,2] ufrMCa[3,3]
 ufrMCa[3,4] ufrMCa[3,5] ufrMCa[3,6] ufrMCa[3,7] ufrMCa[3,8] ufrMCa[3,9]
 ufrMCa[3,10] ufrMCa[3,11] ufrMCa[3,12] ufrMCa[3,13] ufrMCa[3,14]
 ufrMCa[3,15] ufrMCa[4,0] ufrMCa[4,1] ufrMCa[4,2] ufrMCa[4,3] ufrMCa[4,4]
 ufrMCa[4,5] ufrMCa[4,6] ufrMCa[4,7] ufrMCa[4,8] ufrMCa[4,9] ufrMCa[4,10]
 ufrMCa[4,11] ufrMCa[4,12] ufrMCa[4,13] ufrMCa[4,14] ufrMCa[4,15]
 ufrMCa[5,0] ufrMCa[5,1] ufrMCa[5,2] ufrMCa[5,3] ufrMCa[5,4] ufrMCa[5,5]
 ufrMCa[5,6] ufrMCa[5,7] ufrMCa[5,8] ufrMCa[5,9] ufrMCa[5,10] ufrMCa[5,11]
 ufrMCa[5,12] ufrMCa[5,13] ufrMCa[5,14] ufrMCa[5,15] ufrMCa[6,0]
 ufrMCa[6,1] ufrMCa[6,2] ufrMCa[6,3] ufrMCa[6,4] ufrMCa[6,5] ufrMCa[6,6]
 ufrMCa[6,7] ufrMCa[6,8] ufrMCa[6,9] ufrMCa[6,10] ufrMCa[6,11] ufrMCa[6,12]
 ufrMCa[6,13] ufrMCa[6,14] ufrMCa[6,15] lGstk[28,0] lGstk[28,1] lGstk[28,2]
 lGstk[28,3] lGstk[28,4] lGstk[28,5] lGstk[28,6] lGstk[28,7] lGstk[28,8]
 lGstk[28,9] lGstk[28,10] lGstk[28,11] lGstk[28,12] lGstk[28,13]
 lGstk[28,14] lGstk[28,15] lGstk[29,0] lGstk[29,1] lGstk[29,2] lGstk[29,3]
 lGstk[29,4] lGstk[29,5] lGstk[29,6] lGstk[29,7] lGstk[29,8] lGstk[29,9]
 lGstk[29,10] lGstk[29,11] lGstk[29,12] lGstk[29,13] lGstk[29,14]
 lGstk[29,15] lGstk[30,0] lGstk[30,1] lGstk[30,2] lGstk[30,3] lGstk[30,4]
 lGstk[30,5] lGstk[30,6] lGstk[30,7] lGstk[30,8] lGstk[30,9] lGstk[30,10]
 lGstk[30,11] lGstk[30,12] lGstk[30,13] lGstk[30,14] lGstk[30,15]
 lGstk[31,0] lGstk[31,1] lGstk[31,2] lGstk[31,3] lGstk[31,4] lGstk[31,5]
 lGstk[31,6] lGstk[31,7] lGstk[31,8] lGstk[31,9] lGstk[31,10] lGstk[31,11]
 lGstk[31,12] lGstk[31,13] lGstk[31,14] lGstk[31,15] lGstk[32,0]
 lGstk[32,1] lGstk[32,2] lGstk[32,3] lGstk[32,4] lGstk[32,5] lGstk[32,6]
 lGstk[32,7] lGstk[32,8] lGstk[32,9] lGstk[32,10] lGstk[32,11] lGstk[32,12]
 lGstk[32,13] lGstk[32,14] lGstk[32,15] lDetX[28,0] lDetX[28,1] lDetX[28,2]
 lDetX[28,3] lDetX[28,4] lDetX[28,5] lDetX[28,6] lDetX[28,7] lDetX[28,8]
 lDetX[28,9] lDetX[28,10] lDetX[28,11] lDetX[28,12] lDetX[28,13]
 lDetX[28,14] lDetX[28,15] lDetX[29,0] lDetX[29,1] lDetX[29,2] lDetX[29,3]
 lDetX[29,4] lDetX[29,5] lDetX[29,6] lDetX[29,7] lDetX[29,8] lDetX[29,9]
 lDetX[29,10] lDetX[29,11] lDetX[29,12] lDetX[29,13] lDetX[29,14]
 lDetX[29,15] lDetX[30,0] lDetX[30,1] lDetX[30,2] lDetX[30,3] lDetX[30,4]
 lDetX[30,5] lDetX[30,6] lDetX[30,7] lDetX[30,8] lDetX[30,9] lDetX[30,10]
 lDetX[30,11] lDetX[30,12] lDetX[30,13] lDetX[30,14] lDetX[30,15]
 lDetX[31,0] lDetX[31,1] lDetX[31,2] lDetX[31,3] lDetX[31,4] lDetX[31,5]
 lDetX[31,6] lDetX[31,7] lDetX[31,8] lDetX[31,9] lDetX[31,10] lDetX[31,11]
 lDetX[31,12] lDetX[31,13] lDetX[31,14] lDetX[31,15] lDetX[32,0]
 lDetX[32,1] lDetX[32,2] lDetX[32,3] lDetX[32,4] lDetX[32,5] lDetX[32,6]
 lDetX[32,7] lDetX[32,8] lDetX[32,9] lDetX[32,10] lDetX[32,11] lDetX[32,12]
 lDetX[32,13] lDetX[32,14] lDetX[32,15] lDetY[28,0] lDetY[28,1] lDetY[28,2]
 lDetY[28,3] lDetY[28,4] lDetY[28,5] lDetY[28,6] lDetY[28,7] lDetY[28,8]
 lDetY[28,9] lDetY[28,10] lDetY[28,11] lDetY[28,12] lDetY[28,13]
 lDetY[28,14] lDetY[28,15] lDetY[29,0] lDetY[29,1] lDetY[29,2] lDetY[29,3]
 lDetY[29,4] lDetY[29,5] lDetY[29,6] lDetY[29,7] lDetY[29,8] lDetY[29,9]
 lDetY[29,10] lDetY[29,11] lDetY[29,12] lDetY[29,13] lDetY[29,14]
 lDetY[29,15] lDetY[30,0] lDetY[30,1] lDetY[30,2] lDetY[30,3] lDetY[30,4]
 lDetY[30,5] lDetY[30,6] lDetY[30,7] lDetY[30,8] lDetY[30,9] lDetY[30,10]
 lDetY[30,11] lDetY[30,12] lDetY[30,13] lDetY[30,14] lDetY[30,15]
 lDetY[31,0] lDetY[31,1] lDetY[31,2] lDetY[31,3] lDetY[31,4] lDetY[31,5]
 lDetY[31,6] lDetY[31,7] lDetY[31,8] lDetY[31,9] lDetY[31,10] lDetY[31,11]
 lDetY[31,12] lDetY[31,13] lDetY[31,14] lDetY[31,15] lDetY[32,0]
 lDetY[32,1] lDetY[32,2] lDetY[32,3] lDetY[32,4] lDetY[32,5] lDetY[32,6]
 lDetY[32,7] lDetY[32,8] lDetY[32,9] lDetY[32,10] lDetY[32,11] lDetY[32,12]
 lDetY[32,13] lDetY[32,14] lDetY[32,15] lDetZ[28,0] lDetZ[28,1] lDetZ[28,2]
 lDetZ[28,3] lDetZ[28,4] lDetZ[28,5] lDetZ[28,6] lDetZ[28,7] lDetZ[28,8]
 lDetZ[28,9] lDetZ[28,10] lDetZ[28,11] lDetZ[28,12] lDetZ[28,13]
 lDetZ[28,14] lDetZ[28,15] lDetZ[29,0] lDetZ[29,1] lDetZ[29,2] lDetZ[29,3]
 lDetZ[29,4] lDetZ[29,5] lDetZ[29,6] lDetZ[29,7] lDetZ[29,8] lDetZ[29,9]
 lDetZ[29,10] lDetZ[29,11] lDetZ[29,12] lDetZ[29,13] lDetZ[29,14]
 lDetZ[29,15] lDetZ[30,0] lDetZ[30,1] lDetZ[30,2] lDetZ[30,3] lDetZ[30,4]
 lDetZ[30,5] lDetZ[30,6] lDetZ[30,7] lDetZ[30,8] lDetZ[30,9] lDetZ[30,10]
 lDetZ[30,11] lDetZ[30,12] lDetZ[30,13] lDetZ[30,14] lDetZ[30,15]
 lDetZ[31,0] lDetZ[31,1] lDetZ[31,2] lDetZ[31,3] lDetZ[31,4] lDetZ[31,5]
 lDetZ[31,6] lDetZ[31,7] lDetZ[31,8] lDetZ[31,9] lDetZ[31,10] lDetZ[31,11]
 lDetZ[31,12] lDetZ[31,13] lDetZ[31,14] lDetZ[31,15] lDetZ[32,0]
 lDetZ[32,1] lDetZ[32,2] lDetZ[32,3] lDetZ[32,4] lDetZ[32,5] lDetZ[32,6]
 lDetZ[32,7] lDetZ[32,8] lDetZ[32,9] lDetZ[32,10] lDetZ[32,11] lDetZ[32,12]
 lDetZ[32,13] lDetZ[32,14] lDetZ[32,15] lDetW[28,0] lDetW[28,1] lDetW[28,2]
 lDetW[28,3] lDetW[28,4] lDetW[28,5] lDetW[28,6] lDetW[28,7] lDetW[28,8]
 lDetW[28,9] lDetW[28,10] lDetW[28,11] lDetW[28,12] lDetW[28,13]
 lDetW[28,14] lDetW[28,15] lDetW[29,0] lDetW[29,1] lDetW[29,2] lDetW[29,3]
 lDetW[29,4] lDetW[29,5] lDetW[29,6] lDetW[29,7] lDetW[29,8] lDetW[29,9]
 lDetW[29,10] lDetW[29,11] lDetW[29,12] lDetW[29,13] lDetW[29,14]
 lDetW[29,15] lDetW[30,0] lDetW[30,1] lDetW[30,2] lDetW[30,3] lDetW[30,4]
 lDetW[30,5] lDetW[30,6] lDetW[30,7] lDetW[30,8] lDetW[30,9] lDetW[30,10]
 lDetW[30,11] lDetW[30,12] lDetW[30,13] lDetW[30,14] lDetW[30,15]
 lDetW[31,0] lDetW[31,1] lDetW[31,2] lDetW[31,3] lDetW[31,4] lDetW[31,5]
 lDetW[31,6] lDetW[31,7] lDetW[31,8] lDetW[31,9] lDetW[31,10] lDetW[31,11]
 lDetW[31,12] lDetW[31,13] lDetW[31,14] lDetW[31,15] lFrSC[28,0]
 lFrSC[28,1] lFrSC[28,2] lFrSC[28,3] lFrSC[28,4] lFrSC[28,5] lFrSC[28,6]
 lFrSC[28,7] lFrSC[28,8] lFrSC[28,9] lFrSC[28,10] lFrSC[28,11] lFrSC[28,12]
 lFrSC[28,13] lFrSC[28,14] lFrSC[28,15] lFrSC[29,0] lFrSC[29,1] lFrSC[29,2]
 lFrSC[29,3] lFrSC[29,4] lFrSC[29,5] lFrSC[29,6] lFrSC[29,7] lFrSC[29,8]
 lFrSC[29,9] lFrSC[29,10] lFrSC[29,11] lFrSC[29,12] lFrSC[29,13]
 lFrSC[29,14] lFrSC[29,15] lFrSC[30,0] lFrSC[30,1] lFrSC[30,2] lFrSC[30,3]
 lFrSC[30,4] lFrSC[30,5] lFrSC[30,6] lFrSC[30,7] lFrSC[30,8] lFrSC[30,9]
 lFrSC[30,10] lFrSC[30,11] lFrSC[30,12] lFrSC[30,13] lFrSC[30,14]
 lFrSC[30,15] lFrSC[31,0] lFrSC[31,1] lFrSC[31,2] lFrSC[31,3] lFrSC[31,4]
 lFrSC[31,5] lFrSC[31,6] lFrSC[31,7] lFrSC[31,8] lFrSC[31,9] lFrSC[31,10]
 lFrSC[31,11] lFrSC[31,12] lFrSC[31,13] lFrSC[31,14] lFrSC[31,15]
 lFrSC[32,0] lFrSC[32,1] lFrSC[32,2] lFrSC[32,3] lFrSC[32,4] lFrSC[32,5]
 lFrSC[32,6] lFrSC[32,7] lFrSC[32,8] lFrSC[32,9] lFrSC[32,10] lFrSC[32,11]
 lFrSC[32,12] lFrSC[32,13] lFrSC[32,14] lFrSC[32,15] lFrMC[28,0]
 lFrMC[28,1] lFrMC[28,2] lFrMC[28,3] lFrMC[28,4] lFrMC[28,5] lFrMC[28,6]
 lFrMC[28,7] lFrMC[28,8] lFrMC[28,9] lFrMC[28,10] lFrMC[28,11] lFrMC[28,12]
 lFrMC[28,13] lFrMC[28,14] lFrMC[28,15] lFrMC[29,0] lFrMC[29,1] lFrMC[29,2]
 lFrMC[29,3] lFrMC[29,4] lFrMC[29,5] lFrMC[29,6] lFrMC[29,7] lFrMC[29,8]
 lFrMC[29,9] lFrMC[29,10] lFrMC[29,11] lFrMC[29,12] lFrMC[29,13]
 lFrMC[29,14] lFrMC[29,15] lFrMC[30,0] lFrMC[30,1] lFrMC[30,2] lFrMC[30,3]
 lFrMC[30,4] lFrMC[30,5] lFrMC[30,6] lFrMC[30,7] lFrMC[30,8] lFrMC[30,9]
 lFrMC[30,10] lFrMC[30,11] lFrMC[30,12] lFrMC[30,13] lFrMC[30,14]
 lFrMC[30,15] lFrMC[31,0] lFrMC[31,1] lFrMC[31,2] lFrMC[31,3] lFrMC[31,4]
 lFrMC[31,5] lFrMC[31,6] lFrMC[31,7] lFrMC[31,8] lFrMC[31,9] lFrMC[31,10]
 lFrMC[31,11] lFrMC[31,12] lFrMC[31,13] lFrMC[31,14] lFrMC[31,15]
 lgstka[28,0] lgstka[28,1] lgstka[28,2] lgstka[28,3] lgstka[28,4]
 lgstka[28,5] lgstka[28,6] lgstka[28,7] lgstka[28,8] lgstka[28,9]
 lgstka[28,10] lgstka[28,11] lgstka[28,12] lgstka[28,13] lgstka[28,14]
 lgstka[28,15] lgstka[29,0] lgstka[29,1] lgstka[29,2] lgstka[29,3]
 lgstka[29,4] lgstka[29,5] lgstka[29,6] lgstka[29,7] lgstka[29,8]
 lgstka[29,9] lgstka[29,10] lgstka[29,11] lgstka[29,12] lgstka[29,13]
 lgstka[29,14] lgstka[29,15] lgstka[30,0] lgstka[30,1] lgstka[30,2]
 lgstka[30,3] lgstka[30,4] lgstka[30,5] lgstka[30,6] lgstka[30,7]
 lgstka[30,8] lgstka[30,9] lgstka[30,10] lgstka[30,11] lgstka[30,12]
 lgstka[30,13] lgstka[30,14] lgstka[30,15] lgstka[31,0] lgstka[31,1]
 lgstka[31,2] lgstka[31,3] lgstka[31,4] lgstka[31,5] lgstka[31,6]
 lgstka[31,7] lgstka[31,8] lgstka[31,9] lgstka[31,10] lgstka[31,11]
 lgstka[31,12] lgstka[31,13] lgstka[31,14] lgstka[31,15] lgstka[32,0]
 lgstka[32,1] lgstka[32,2] lgstka[32,3] lgstka[32,4] lgstka[32,5]
 lgstka[32,6] lgstka[32,7] lgstka[32,8] lgstka[32,9] lgstka[32,10]
 lgstka[32,11] lgstka[32,12] lgstka[32,13] lgstka[32,14] lgstka[32,15]
 ldetXa[28,0] ldetXa[28,1] ldetXa[28,2] ldetXa[28,3] ldetXa[28,4]
 ldetXa[28,5] ldetXa[28,6] ldetXa[28,7] ldetXa[28,8] ldetXa[28,9]
 ldetXa[28,10] ldetXa[28,11] ldetXa[28,12] ldetXa[28,13] ldetXa[28,14]
 ldetXa[28,15] ldetXa[29,0] ldetXa[29,1] ldetXa[29,2] ldetXa[29,3]
 ldetXa[29,4] ldetXa[29,5] ldetXa[29,6] ldetXa[29,7] ldetXa[29,8]
 ldetXa[29,9] ldetXa[29,10] ldetXa[29,11] ldetXa[29,12] ldetXa[29,13]
 ldetXa[29,14] ldetXa[29,15] ldetXa[30,0] ldetXa[30,1] ldetXa[30,2]
 ldetXa[30,3] ldetXa[30,4] ldetXa[30,5] ldetXa[30,6] ldetXa[30,7]
 ldetXa[30,8] ldetXa[30,9] ldetXa[30,10] ldetXa[30,11] ldetXa[30,12]
 ldetXa[30,13] ldetXa[30,14] ldetXa[30,15] ldetXa[31,0] ldetXa[31,1]
 ldetXa[31,2] ldetXa[31,3] ldetXa[31,4] ldetXa[31,5] ldetXa[31,6]
 ldetXa[31,7] ldetXa[31,8] ldetXa[31,9] ldetXa[31,10] ldetXa[31,11]
 ldetXa[31,12] ldetXa[31,13] ldetXa[31,14] ldetXa[31,15] ldetXa[32,0]
 ldetXa[32,1] ldetXa[32,2] ldetXa[32,3] ldetXa[32,4] ldetXa[32,5]
 ldetXa[32,6] ldetXa[32,7] ldetXa[32,8] ldetXa[32,9] ldetXa[32,10]
 ldetXa[32,11] ldetXa[32,12] ldetXa[32,13] ldetXa[32,14] ldetXa[32,15]
 ldetYa[28,0] ldetYa[28,1] ldetYa[28,2] ldetYa[28,3] ldetYa[28,4]
 ldetYa[28,5] ldetYa[28,6] ldetYa[28,7] ldetYa[28,8] ldetYa[28,9]
 ldetYa[28,10] ldetYa[28,11] ldetYa[28,12] ldetYa[28,13] ldetYa[28,14]
 ldetYa[28,15] ldetYa[29,0] ldetYa[29,1] ldetYa[29,2] ldetYa[29,3]
 ldetYa[29,4] ldetYa[29,5] ldetYa[29,6] ldetYa[29,7] ldetYa[29,8]
 ldetYa[29,9] ldetYa[29,10] ldetYa[29,11] ldetYa[29,12] ldetYa[29,13]
 ldetYa[29,14] ldetYa[29,15] ldetYa[30,0] ldetYa[30,1] ldetYa[30,2]
 ldetYa[30,3] ldetYa[30,4] ldetYa[30,5] ldetYa[30,6] ldetYa[30,7]
 ldetYa[30,8] ldetYa[30,9] ldetYa[30,10] ldetYa[30,11] ldetYa[30,12]
 ldetYa[30,13] ldetYa[30,14] ldetYa[30,15] ldetYa[31,0] ldetYa[31,1]
 ldetYa[31,2] ldetYa[31,3] ldetYa[31,4] ldetYa[31,5] ldetYa[31,6]
 ldetYa[31,7] ldetYa[31,8] ldetYa[31,9] ldetYa[31,10] ldetYa[31,11]
 ldetYa[31,12] ldetYa[31,13] ldetYa[31,14] ldetYa[31,15] ldetYa[32,0]
 ldetYa[32,1] ldetYa[32,2] ldetYa[32,3] ldetYa[32,4] ldetYa[32,5]
 ldetYa[32,6] ldetYa[32,7] ldetYa[32,8] ldetYa[32,9] ldetYa[32,10]
 ldetYa[32,11] ldetYa[32,12] ldetYa[32,13] ldetYa[32,14] ldetYa[32,15]
 ldetZa[28,0] ldetZa[28,1] ldetZa[28,2] ldetZa[28,3] ldetZa[28,4]
 ldetZa[28,5] ldetZa[28,6] ldetZa[28,7] ldetZa[28,8] ldetZa[28,9]
 ldetZa[28,10] ldetZa[28,11] ldetZa[28,12] ldetZa[28,13] ldetZa[28,14]
 ldetZa[28,15] ldetZa[29,0] ldetZa[29,1] ldetZa[29,2] ldetZa[29,3]
 ldetZa[29,4] ldetZa[29,5] ldetZa[29,6] ldetZa[29,7] ldetZa[29,8]
 ldetZa[29,9] ldetZa[29,10] ldetZa[29,11] ldetZa[29,12] ldetZa[29,13]
 ldetZa[29,14] ldetZa[29,15] ldetZa[30,0] ldetZa[30,1] ldetZa[30,2]
 ldetZa[30,3] ldetZa[30,4] ldetZa[30,5] ldetZa[30,6] ldetZa[30,7]
 ldetZa[30,8] ldetZa[30,9] ldetZa[30,10] ldetZa[30,11] ldetZa[30,12]
 ldetZa[30,13] ldetZa[30,14] ldetZa[30,15] ldetZa[31,0] ldetZa[31,1]
 ldetZa[31,2] ldetZa[31,3] ldetZa[31,4] ldetZa[31,5] ldetZa[31,6]
 ldetZa[31,7] ldetZa[31,8] ldetZa[31,9] ldetZa[31,10] ldetZa[31,11]
 ldetZa[31,12] ldetZa[31,13] ldetZa[31,14] ldetZa[31,15] ldetZa[32,0]
 ldetZa[32,1] ldetZa[32,2] ldetZa[32,3] ldetZa[32,4] ldetZa[32,5]
 ldetZa[32,6] ldetZa[32,7] ldetZa[32,8] ldetZa[32,9] ldetZa[32,10]
 ldetZa[32,11] ldetZa[32,12] ldetZa[32,13] ldetZa[32,14] ldetZa[32,15]
 ldetWa[28,0] ldetWa[28,1] ldetWa[28,2] ldetWa[28,3] ldetWa[28,4]
 ldetWa[28,5] ldetWa[28,6] ldetWa[28,7] ldetWa[28,8] ldetWa[28,9]
 ldetWa[28,10] ldetWa[28,11] ldetWa[28,12] ldetWa[28,13] ldetWa[28,14]
 ldetWa[28,15] ldetWa[29,0] ldetWa[29,1] ldetWa[29,2] ldetWa[29,3]
 ldetWa[29,4] ldetWa[29,5] ldetWa[29,6] ldetWa[29,7] ldetWa[29,8]
 ldetWa[29,9] ldetWa[29,10] ldetWa[29,11] ldetWa[29,12] ldetWa[29,13]
 ldetWa[29,14] ldetWa[29,15] ldetWa[30,0] ldetWa[30,1] ldetWa[30,2]
 ldetWa[30,3] ldetWa[30,4] ldetWa[30,5] ldetWa[30,6] ldetWa[30,7]
 ldetWa[30,8] ldetWa[30,9] ldetWa[30,10] ldetWa[30,11] ldetWa[30,12]
 ldetWa[30,13] ldetWa[30,14] ldetWa[30,15] ldetWa[31,0] ldetWa[31,1]
 ldetWa[31,2] ldetWa[31,3] ldetWa[31,4] ldetWa[31,5] ldetWa[31,6]
 ldetWa[31,7] ldetWa[31,8] ldetWa[31,9] ldetWa[31,10] ldetWa[31,11]
 ldetWa[31,12] ldetWa[31,13] ldetWa[31,14] ldetWa[31,15] lfrSCa[28,0]
 lfrSCa[28,1] lfrSCa[28,2] lfrSCa[28,3] lfrSCa[28,4] lfrSCa[28,5]
 lfrSCa[28,6] lfrSCa[28,7] lfrSCa[28,8] lfrSCa[28,9] lfrSCa[28,10]
 lfrSCa[28,11] lfrSCa[28,12] lfrSCa[28,13] lfrSCa[28,14] lfrSCa[28,15]
 lfrSCa[29,0] lfrSCa[29,1] lfrSCa[29,2] lfrSCa[29,3] lfrSCa[29,4]
 lfrSCa[29,5] lfrSCa[29,6] lfrSCa[29,7] lfrSCa[29,8] lfrSCa[29,9]
 lfrSCa[29,10] lfrSCa[29,11] lfrSCa[29,12] lfrSCa[29,13] lfrSCa[29,14]
 lfrSCa[29,15] lfrSCa[30,0] lfrSCa[30,1] lfrSCa[30,2] lfrSCa[30,3]
 lfrSCa[30,4] lfrSCa[30,5] lfrSCa[30,6] lfrSCa[30,7] lfrSCa[30,8]
 lfrSCa[30,9] lfrSCa[30,10] lfrSCa[30,11] lfrSCa[30,12] lfrSCa[30,13]
 lfrSCa[30,14] lfrSCa[30,15] lfrSCa[31,0] lfrSCa[31,1] lfrSCa[31,2]
 lfrSCa[31,3] lfrSCa[31,4] lfrSCa[31,5] lfrSCa[31,6] lfrSCa[31,7]
 lfrSCa[31,8] lfrSCa[31,9] lfrSCa[31,10] lfrSCa[31,11] lfrSCa[31,12]
 lfrSCa[31,13] lfrSCa[31,14] lfrSCa[31,15] lfrSCa[32,0] lfrSCa[32,1]
 lfrSCa[32,2] lfrSCa[32,3] lfrSCa[32,4] lfrSCa[32,5] lfrSCa[32,6]
 lfrSCa[32,7] lfrSCa[32,8] lfrSCa[32,9] lfrSCa[32,10] lfrSCa[32,11]
 lfrSCa[32,12] lfrSCa[32,13] lfrSCa[32,14] lfrSCa[32,15] lfrMCa[28,0]
 lfrMCa[28,1] lfrMCa[28,2] lfrMCa[28,3] lfrMCa[28,4] lfrMCa[28,5]
 lfrMCa[28,6] lfrMCa[28,7] lfrMCa[28,8] lfrMCa[28,9] lfrMCa[28,10]
 lfrMCa[28,11] lfrMCa[28,12] lfrMCa[28,13] lfrMCa[28,14] lfrMCa[28,15]
 lfrMCa[29,0] lfrMCa[29,1] lfrMCa[29,2] lfrMCa[29,3] lfrMCa[29,4]
 lfrMCa[29,5] lfrMCa[29,6] lfrMCa[29,7] lfrMCa[29,8] lfrMCa[29,9]
 lfrMCa[29,10] lfrMCa[29,11] lfrMCa[29,12] lfrMCa[29,13] lfrMCa[29,14]
 lfrMCa[29,15] lfrMCa[30,0] lfrMCa[30,1] lfrMCa[30,2] lfrMCa[30,3]
 lfrMCa[30,4] lfrMCa[30,5] lfrMCa[30,6] lfrMCa[30,7] lfrMCa[30,8]
 lfrMCa[30,9] lfrMCa[30,10] lfrMCa[30,11] lfrMCa[30,12] lfrMCa[30,13]
 lfrMCa[30,14] lfrMCa[30,15] lfrMCa[31,0] lfrMCa[31,1] lfrMCa[31,2]
 lfrMCa[31,3] lfrMCa[31,4] lfrMCa[31,5] lfrMCa[31,6] lfrMCa[31,7]
 lfrMCa[31,8] lfrMCa[31,9] lfrMCa[31,10] lfrMCa[31,11] lfrMCa[31,12]
 lfrMCa[31,13] lfrMCa[31,14] lfrMCa[31,15] ikz[0] ikz[1] ikz[2] ikz[3]
 ikz[4] ikz[5] ikz[6] ikz[7] ikz[8] ikz[9] ikz[10] ikz[11] ikz[12] ikz[13]
 ikz[14] ikz[15] gkz[0] gkz[1] gkz[2] gkz[3] gkz[4] gkz[5] gkz[6] gkz[7]
 gkz[8] gkz[9] gkz[10] gkz[11] gkz[12] gkz[13] gkz[14] gkz[15]
Generals
 iLANE[0] iLANE[1] iLANE[2] iLANE[3] iLANE[4] iLANE[5] iLANE[6] iLANE[7]
 iLANE[8] iLANE[9] iLANE[10] iLANE[11] iLANE[12] iLANE[13] iLANE[14]
 iLANE[15] ikSUM[0] ikSUM[1] ikSUM[2] ikSUM[3] ikSUM[4] ikSUM[5] ikSUM[6]
 ikSUM[7] ikSUM[8] ikSUM[9] ikSUM[10] ikSUM[11] ikSUM[12] ikSUM[13]
 ikSUM[14] ikSUM[15] gLANE[0] gLANE[1] gLANE[2] gLANE[3] gLANE[4] gLANE[5]
 gLANE[6] gLANE[7] gLANE[8] gLANE[9] gLANE[10] gLANE[11] gLANE[12]
 gLANE[13] gLANE[14] gLANE[15] gkSUM[0] gkSUM[1] gkSUM[2] gkSUM[3] gkSUM[4]
 gkSUM[5] gkSUM[6] gkSUM[7] gkSUM[8] gkSUM[9] gkSUM[10] gkSUM[11] gkSUM[12]
 gkSUM[13] gkSUM[14] gkSUM[15] rb rf cb cf iks gks cbp cfp mbi mfi
End
