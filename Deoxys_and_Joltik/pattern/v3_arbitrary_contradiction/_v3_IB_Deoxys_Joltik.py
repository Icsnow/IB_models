import gurobipy as gp
from gurobipy import GRB
import math

class IB_DandJ:

  def __init__(self, key_size, round_Eb, round_Dis, round_Ef, setX, pgP) -> None:
    if key_size in [128,192]:
      self.block_size = 64
    else:
      self.block_size = 128
    self.key_size = key_size
    self.s = self.key_size // self.block_size
    self.cell_size = self.block_size // 16
    self.round_Eb = round_Eb
    self.round_Dis = round_Dis
    self.round_Ef = round_Ef
    self.pgP = pgP
    self.tx = setX
    self.tz = self.x_to_z(setX)
    self.SRpermutation_rev = [0, 5, 10, 15, 4, 9, 14, 3, 8, 13, 2, 7, 12, 1, 6, 11]
    self.hTable = [[ 0,  7, 14,  9,  8, 15,  6,  1,  0,  7, 14,  9,  8, 15,  6,  1],
                   [ 1,  0,  7, 14,  9,  8, 15,  6,  1,  0,  7, 14,  9,  8, 15,  6],
                   [ 2, 13, 12,  3, 10,  5,  4, 11,  2, 13, 12,  3, 10,  5,  4, 11],
                   [ 3, 10,  5,  4, 11,  2, 13, 12,  3, 10,  5,  4, 11,  2, 13, 12],
                   [ 4, 11,  2, 13, 12,  3, 10,  5,  4, 11,  2, 13, 12,  3, 10,  5],
                   [ 5,  4, 11,  2, 13, 12,  3, 10,  5,  4, 11,  2, 13, 12,  3, 10],
                   [ 6,  1,  0,  7, 14,  9,  8, 15,  6,  1,  0,  7, 14,  9,  8, 15],
                   [ 7, 14,  9,  8, 15,  6,  1,  0,  7, 14,  9,  8, 15,  6,  1,  0],
                   [ 8, 15,  6,  1,  0,  7, 14,  9,  8, 15,  6,  1,  0,  7, 14,  9],
                   [ 9,  8, 15,  6,  1,  0,  7, 14,  9,  8, 15,  6,  1,  0,  7, 14],
                   [10,  5,  4, 11,  2, 13, 12,  3, 10,  5,  4, 11,  2, 13, 12,  3],
                   [11,  2, 13, 12,  3, 10,  5,  4, 11,  2, 13, 12,  3, 10,  5,  4],
                   [12,  3, 10,  5,  4, 11,  2, 13, 12,  3, 10,  5,  4, 11,  2, 13],
                   [13, 12,  3, 10,  5,  4, 11,  2, 13, 12,  3, 10,  5,  4, 11,  2],
                   [14,  9,  8, 15,  6,  1,  0,  7, 14,  9,  8, 15,  6,  1,  0,  7],
                   [15,  6,  1,  0,  7, 14,  9,  8, 15,  6,  1,  0,  7, 14,  9,  8]]
    if self.cell_size == 4:
      cipher_name = 'JoltikBC'
    else:
      cipher_name = 'DeoxysBC'
    self.name = './v3_' +  cipher_name + str(self.key_size) + '_' + str(self.round_Eb + self.round_Dis + self.round_Ef) + 'r'
    self.model = gp.Model(self.name)

  ''' 
  <=======================================================================================================>
  Model generation:
  ONE.   Differential propagation in (1).Eb(<-), (2).Eu(->), (3).Miss In The Middle, (4).El(<-), (5).Ef(->)
  TWO.   Key Recovery in upper and lower trails (no direction)
  THREE. Guess, Determine and Filter in Eb(->) and Ef(<-)
  FOUR.  Miss-In-The-Middle
  FIVE.  Complexities
  <=======================================================================================================>
  '''

  def ib_model(self):
    '''
    <------------------------------------------------------------------------------------
    ONE_(1-2). Differential propagation in (1).Eb(<-), (2).Eu(->):
    *PARA.: 
      DX[round,active tag], dx[round, truncated tag]
    *TRAIL: 
      Eb: DWw[-1,](plaintext) <={ART(stk)}= DXx <={SC}= DYy <={SR}= DZz <={MC}= DWz[Eb,]
      Eu: DWw[-1,](plaintext) ={ART(stk)}=> DXx ={SC}=> DYy ={SR}=> DZz ={MC}=> DWz[Eb,]
    ------------------------------------------------------------------------------------>
    '''
    end_round_u = self.round_Eb + self.round_Dis
    uDW  = self.model.addVars(range(-1, end_round_u),          16, vtype = GRB.BINARY, name = 'uDW')
    udw  = self.model.addVars(range(-1, end_round_u),          16, vtype = GRB.BINARY, name = 'udw')# DWw[-1,](plaintext)
    uDX  = self.model.addVars(end_round_u, 16, vtype = GRB.BINARY, name = 'uDX')
    udx  = self.model.addVars(end_round_u, 16, vtype = GRB.BINARY, name = 'udx')
    uDY  = self.model.addVars(end_round_u,                     16, vtype = GRB.BINARY, name = 'uDY')
    udy  = self.model.addVars(end_round_u,                     16, vtype = GRB.BINARY, name = 'udy')
    uDZ  = self.model.addVars(end_round_u,                     16, vtype = GRB.BINARY, name = 'uDZ')
    udz  = self.model.addVars(end_round_u,                     16, vtype = GRB.BINARY, name = 'udz')
    ustk = self.model.addVars(end_round_u + 1, 16, vtype = GRB.BINARY, name = 'ustk')
    ucan = self.model.addVars(end_round_u,     16, vtype = GRB.BINARY, name = 'ucan')
    uAC  = self.model.addVars(end_round_u,                      4, vtype = GRB.BINARY, name = 'uAC') # Mark active in one column
    uTC  = self.model.addVars(end_round_u,                      4, vtype = GRB.BINARY, name = 'uTC') # Mark truncated in one column

    '''Basic Constrs:'''
    # Remove dx[r,i] > DX[r,i]
    for i in range(16):
      for r in range(-1, end_round_u):
        self.model.addConstr(uDW[r,i] >= udw[r,i], name = 'uBas')
      for r in range(end_round_u):
        self.model.addConstr(uDX[r,i] >= udx[r,i], name = 'uBas')
        self.model.addConstr(uDY[r,i] >= udy[r,i], name = 'uBas')
        self.model.addConstr(uDZ[r,i] >= udz[r,i], name = 'uBas')
        
    '''Operation: ATK'''
    for r in range(-1, end_round_u - 1):
      for i in range(16):
        self.model.addConstr(uDW[r,i] - udx[r+1,i] - ucan[r+1,i] >= 0, name = 'uATK')
        self.model.addConstr(- uDW[r,i] + uDX[r+1,i] + ucan[r+1,i] >= 0, name = 'uATK')
        self.model.addConstr(uDW[r,i] + ustk[r+1,i] - uDX[r+1,i] - 2 * ucan[r+1,i] >= 0, name = 'uATK')
        self.model.addConstr(- ustk[r+1,i] + uDX[r+1,i] + ucan[r+1,i] >= 0, name = 'uATK')
        self.model.addConstr(udw[r,i] == udx[r+1,i], name = 'uATK')

    '''Operation: SC'''
    # In Extension     [DX,dx] <= [DY,dy]: [1,1] <= [1,0]/[1,1], [0,0] <= [0,0]
    for r in range(0, self.round_Eb):
      for i in range(16):
        self.model.addConstr(uDX[r,i] == uDY[r,i], name = 'bSC')
        self.model.addConstr(udx[r,i] == uDY[r,i], name = 'bSC')
    # In distinguisher [DX,dx] => [DY,dy]: [1,0]/[1,1] => [1,1], [0,0] => [0,0]
    for r in range(self.round_Eb, self.round_Eb + self.round_Dis):
      for i in range(16):
        self.model.addConstr(uDY[r,i] == uDX[r,i], name = 'uSC')
        self.model.addConstr(udy[r,i] == uDX[r,i], name = 'uSC')

    '''Operation: SR''' 
    # [DY,dy] <=> (DZ,dz)
    for r in range(self.round_Eb + self.round_Dis):
      for i in range(16):
        self.model.addConstr(uDY[r,self.SRpermutation_rev[i]] == uDZ[r,i], name = 'uSR')
        self.model.addConstr(udy[r,self.SRpermutation_rev[i]] == udz[r,i], name = 'uSR')

    '''Operation: MC'''
    # In Extension: [DZ,dz] => [DW,dw]
    for r in range(self.round_Eb):
      for j in [0,4,8,12]:
        self.model.addConstr(uTC[r,j/4] >= udw[r,j+0], name = 'bMC')
        self.model.addConstr(uTC[r,j/4] >= udw[r,j+1], name = 'bMC')
        self.model.addConstr(uTC[r,j/4] >= udw[r,j+2], name = 'bMC')
        self.model.addConstr(uTC[r,j/4] >= udw[r,j+3], name = 'bMC')
        self.model.addConstr(uTC[r,j/4] <= sum(udw[r,j+i] for i in range(4)), name = 'bMC')
        self.model.addConstr(uTC[r,j/4] <= udz[r,j+0], name = 'bMC')
        self.model.addConstr(uTC[r,j/4] <= udz[r,j+1], name = 'bMC')
        self.model.addConstr(uTC[r,j/4] <= udz[r,j+2], name = 'bMC')
        self.model.addConstr(uTC[r,j/4] <= udz[r,j+3], name = 'bMC')
        self.model.addConstr(4*uTC[r,j/4] >= sum(udz[r,j+i] for i in range(4)), name = 'bMC')

        self.model.addConstr(uAC[r,j/4] >= uDW[r,j+0], name = 'bMC')
        self.model.addConstr(uAC[r,j/4] >= uDW[r,j+1], name = 'bMC')
        self.model.addConstr(uAC[r,j/4] >= uDW[r,j+2], name = 'bMC')
        self.model.addConstr(uAC[r,j/4] >= uDW[r,j+3], name = 'bMC')
        self.model.addConstr(uAC[r,j/4] <= sum(uDW[r,j+i] for i in range(4)), name = 'bMC')
        self.model.addConstr(uDW[r,j+0] + uDW[r,j+1] + uDW[r,j+2] + uDW[r,j+3] + 
                             uDZ[r,j+0] + uDZ[r,j+1] + uDZ[r,j+2] + uDZ[r,j+3] >= 5*uAC[r,j/4], name = 'bMC')
        self.model.addConstr(uDW[r,j+0] + uDW[r,j+1] + uDW[r,j+2] + uDW[r,j+3] + 
                             uDZ[r,j+0] + uDZ[r,j+1] + uDZ[r,j+2] + uDZ[r,j+3] <= 8*uAC[r,j/4], name = 'bMC')
    # In Distinguisher [DW,dw] <= [DZ,dz]
    for r in range(self.round_Eb, self.round_Eb + self.round_Dis):
      for j in [0,4,8,12]:
        self.model.addConstr(uTC[r,j/4] >= udz[r,j+0], name = 'uMC')
        self.model.addConstr(uTC[r,j/4] >= udz[r,j+1], name = 'uMC')
        self.model.addConstr(uTC[r,j/4] >= udz[r,j+2], name = 'uMC')
        self.model.addConstr(uTC[r,j/4] >= udz[r,j+3], name = 'uMC')
        self.model.addConstr(uTC[r,j/4] <= sum(udz[r,j+i] for i in range(4)), name = 'uMC')
        self.model.addConstr(uTC[r,j/4] <= udw[r,j+0], name = 'uMC')
        self.model.addConstr(uTC[r,j/4] <= udw[r,j+1], name = 'uMC')
        self.model.addConstr(uTC[r,j/4] <= udw[r,j+2], name = 'uMC')
        self.model.addConstr(uTC[r,j/4] <= udw[r,j+3], name = 'uMC')
        self.model.addConstr(4*uTC[r,j/4] >= sum(udw[r,j+i] for i in range(4)), name = 'uMC')

        self.model.addConstr(uAC[r,j/4] >= uDZ[r,j+0], name = 'uMC')
        self.model.addConstr(uAC[r,j/4] >= uDZ[r,j+1], name = 'uMC')
        self.model.addConstr(uAC[r,j/4] >= uDZ[r,j+2], name = 'uMC')
        self.model.addConstr(uAC[r,j/4] >= uDZ[r,j+3], name = 'uMC')
        self.model.addConstr(uAC[r,j/4] <= sum(uDZ[r,j+i] for i in range(4)), name = 'uMC')
        self.model.addConstr(uDZ[r,j+0] + uDZ[r,j+1] + uDZ[r,j+2] + uDZ[r,j+3] + 
                             uDW[r,j+0] + uDW[r,j+1] + uDW[r,j+2] + uDW[r,j+3] >= 5*uAC[r,j/4], name = 'uMC')
        self.model.addConstr(uDZ[r,j+0] + uDZ[r,j+1] + uDZ[r,j+2] + uDZ[r,j+3] + 
                             uDW[r,j+0] + uDW[r,j+1] + uDW[r,j+2] + uDW[r,j+3] <= 8*uAC[r,j/4], name = 'uMC')

    '''
    -----------------------------------------------------------------------------
    TWO: Key schedule in Eb+Eu+Em (No direction)
    *PARA.: 
      uLANE:         A LANE of permutation of key schedule (for Type1 cancellation)
      uT2Can(A/B/C): Auxiliary variables to count Type2 cancellation
      uType2:        Type2 cancellation
    -----------------------------------------------------------------------------
    '''
    uLANE   = self.model.addVars(16, vtype = GRB.BINARY, name = 'uLANE')
    uT2CanA = self.model.addVars(end_round_u, 4, lb = 0, vtype = GRB.INTEGER, name = 'uT2CanA')
    uT2CanB = self.model.addVars(end_round_u, 4, lb = 0, vtype = GRB.INTEGER, name = 'uT2CanB')
    uT2CanC = self.model.addVars(range(1, end_round_u), 4, lb = 0, vtype = GRB.INTEGER, name = 'uT2CanC')
    uType1  = self.model.addVars(16, lb = 0, vtype = GRB.INTEGER, name = 'uType1')
    uType2  = self.model.addVars(end_round_u, 4, lb = 0, vtype = GRB.INTEGER, name = 'uType2')

    self.model.addConstr(sum(uLANE[i] for i in range(16)) >= 1)

    ''' Mark LANE '''
    for i in range(16):
      for r in range(end_round_u + 1):
        self.model.addConstr(uLANE[i] >= ustk[r, self.hTable[i][r]], name = 'uMLANE')

    ''' Type 1 cancellation '''
    for i in range(16):
      self.model.addConstr(uType1[i] == (end_round_u + 1) * uLANE[i] - 
                           sum(ustk[r,self.hTable[i][r]] for r in range(end_round_u + 1)), name = 'uType1')
      self.model.addConstr(uType1[i] <= self.s - 1, name = 'uType1') # Cancellation no more than 2/1 in each position (for Deoxys-BC-384/256)
        
    ''' Type 2 cancellation (In Eu and Em) '''
    for r in range(end_round_u):
      for c in range(4):
        # uT2CanA: Active bytes in one column, before MC
        self.model.addConstr(uT2CanA[r,c] == uDZ[r,4*c+0] + uDZ[r,4*c+1] + uDZ[r,4*c+2] + uDZ[r,4*c+3], name = 'uType2A')
        # uT2CanB: Inactive bytes in one column, after MC
        self.model.addConstr(uT2CanB[r,c] == 4*uAC[r,c] - uDW[r,4*c+0] - uDW[r,4*c+1] - uDW[r,4*c+2] - uDW[r,4*c+3], name = 'uType2B')                             

    # uT2CanC: Cancellation appears in AK, in the next round
    for r in range(1, end_round_u):
      for c in range(4):
        self.model.addConstr(uT2CanC[r,c] == ucan[r,4*c+0] + ucan[r,4*c+1] + ucan[r,4*c+2] + ucan[r,4*c+3], name = 'uType2C')

    # uType2 = -(uT2CanA - uT2CanB - uT2CanC)
    for r in range(end_round_u-1):
      for c in range(4):
        self.model.addConstr(uType2[r,c] >= (uT2CanB[r,c] + uT2CanC[r+1,c] - uT2CanA[r,c]), name = 'uType2')
        self.model.addConstr(uType2[r,c] >= 0, name = 'uType2')

    ''' Counting of all cancellations'''
    # NOTE: s*LANE[0~15] - 1 >= Type1Can * (LANE[0~15]-stk[0~r,0~15]) + Type2Can.
    self.model.addConstr(self.s * sum(uLANE[i] for i in range(16)) -  # s * LANE
                         sum(uType1[i] for i in range(16)) - 
                         sum(uType2[r,i] for r in range(end_round_u) for i in range(4))
                         >= 1, name = 'uCAN')

    '''
    -----------------------------------------------------------------------------------------------
    THREE: Guess-and-Determine (Eb)
      PRAR.: 
        uGstk: Whether the stk is guessed or not
        uDetW[-1,](plaintext) ={ART(uGstk)}=> uDetX ={SC}=> uDetY ={SR}=> uDetZ ={MC}=> uDetW [Eb,]
        Auxiliary variable of MC: uDetC, uDiffDetMC, uDiffDetC
        Filter tag: uFrSB, uFrMC
    -----------------------------------------------------------------------------------------------
    ''' 
    uGstk = self.model.addVars(self.round_Eb,            16, vtype = GRB.BINARY, name = 'uGstk')
    uDetX = self.model.addVars(self.round_Eb,            16, vtype = GRB.BINARY, name = 'uDetX')
    uFrSB = self.model.addVars(self.round_Eb,            16, vtype = GRB.BINARY, name = 'uFrSB')
    uDetY = self.model.addVars(self.round_Eb,            16, vtype = GRB.BINARY, name = 'uDetY')
    uDetZ = self.model.addVars(self.round_Eb,            16, vtype = GRB.BINARY, name = 'uDetZ')
    uDetW = self.model.addVars(range(-1, self.round_Eb), 16, vtype = GRB.BINARY, name = 'uDetW')
    uDetC = self.model.addVars(self.round_Eb,             4, vtype = GRB.BINARY, name = 'uDetC') # mark 'det.' in one column
    uDiffDetMC = self.model.addVars(self.round_Eb,       16, vtype = GRB.BINARY, name = 'uDiffDetMC') # mark diff. det. of cell
    uDiffDetC = self.model.addVars(self.round_Eb,         4, vtype = GRB.BINARY, name = 'uDiffDetC') # mark diff. det. in one colume
    uFrMC = self.model.addVars(self.round_Eb,            16, vtype = GRB.BINARY, name = 'uFrMC')
    
    ''' KR: Guess round key and determine '''
    # uDetW -(uGstk)-> uDetX
    for r in range(-1, self.round_Eb-1):
      for i in range(16): 
        self.model.addConstr(uGstk[r+1,i] >= uDetX[r+1,i], name = 'bDet_ATK')
        self.model.addConstr(uDetW[r,i] >= uDetX[r+1,i], name = 'bDet_ATK')
        self.model.addConstr(- uDetW[r,i] - uGstk[r+1,i] + uDetX[r+1,i] >= -1, name = 'bDet_ATK')

    ''' KR: Determine and Filter in SC '''
    for r in range(self.round_Eb):
      for i in range(16):
        # 'Det.' propagation via SC
        self.model.addConstr(uDetY[r,i] == uDetX[r,i], name = 'bDet_SC') 
        # Filter obtain from SC
        self.model.addConstr(uDY[r,i] - udy[r,i] - uFrSB[r,i] >= 0, name = 'bFr_SC')
        self.model.addConstr(uDetY[r,i] >= uFrSB[r,i], name = 'bFr_SC')
        self.model.addConstr(- uDetY[r,i] - uDY[r,i] + udy[r,i] + uFrSB[r,i] >= -1, name = 'bFr_SC')
    
    ''' KR: 'Det.' propagation in SR '''
    for r in range(self.round_Eb):
      for i in range(16):
        self.model.addConstr(uDetZ[r,i] == uDetY[r,self.SRpermutation_rev[i]], name = 'bDet_SR')

    ''' KR: 'Det.' propagation in MC* '''
    for r in range(self.round_Eb):
      for c in range(4):
        # 'Det.' propagation via MC
        # NOTE: uDetC = 1 when all uDetZ = 1, else uDetC = 0 
        self.model.addConstr(uDetC[r,c] <= uDetZ[r,4*c+0], name = 'bDet_MC')
        self.model.addConstr(uDetC[r,c] <= uDetZ[r,4*c+1], name = 'bDet_MC')
        self.model.addConstr(uDetC[r,c] <= uDetZ[r,4*c+2], name = 'bDet_MC')
        self.model.addConstr(uDetC[r,c] <= uDetZ[r,4*c+3], name = 'bDet_MC')
        self.model.addConstr(uDetZ[r,4*c+0] + uDetZ[r,4*c+1] + uDetZ[r,4*c+2] + uDetZ[r,4*c+3] - uDetC[r,c] <= 3, name = 'bDet_MC')
        # NOTE: all uDetW = 1 when uDetC = 1, else all uDetW = 0
        self.model.addConstr(uDetC[r,c] == uDetW[r,4*c+0], name = 'bDet_MC')
        self.model.addConstr(uDetC[r,c] == uDetW[r,4*c+1], name = 'bDet_MC')
        self.model.addConstr(uDetC[r,c] == uDetW[r,4*c+2], name = 'bDet_MC')
        self.model.addConstr(uDetC[r,c] == uDetW[r,4*c+3], name = 'bDet_MC')

    ''' KR: Determine whether the difference before MC is determinable '''
    for r in range(self.round_Eb):
      for i in range(16):
        # NOTE: uDiffDetMC = uDetZ when uDZ = 1 (ACTIVE), else uDiffDetMC = 1
        self.model.addConstr(uDiffDetMC[r,i] >= uDetZ[r,i], name = 'bDet_MC')
        self.model.addConstr(uDetZ[r,i] - uDZ[r,i] - uDiffDetMC[r,i] >= -1, name = 'bDet_MC')
        self.model.addConstr(uDZ[r,i] + uDiffDetMC[r,i] >= 1, name = 'bDet_MC')
    for r in range(self.round_Eb):
      for c in range(4):
        # NOTE: uDiffDetC = 1 when all uDiffDetMC = 1, else uDiffDetC = 0 
        self.model.addConstr(uDiffDetC[r,c] <= uDiffDetMC[r,4*c+0], name = 'bDet_MC')
        self.model.addConstr(uDiffDetC[r,c] <= uDiffDetMC[r,4*c+1], name = 'bDet_MC')
        self.model.addConstr(uDiffDetC[r,c] <= uDiffDetMC[r,4*c+2], name = 'bDet_MC')
        self.model.addConstr(uDiffDetC[r,c] <= uDiffDetMC[r,4*c+3], name = 'bDet_MC')
        self.model.addConstr(uDiffDetMC[r,4*c+0] + uDiffDetMC[r,4*c+1] + uDiffDetMC[r,4*c+2] + uDiffDetMC[r,4*c+3] 
                             - uDiffDetC[r,c] <= 3, name = 'bDet_MC')
    
    ''' KR: Filter obtain from MC (for cell -> column) '''
    for r in range(self.round_Eb):
      for c in range(4):
        for i in range(4):
          # NOTE: uFrMC = 0 when uDiffDetC = 0, 
          # NOTE: else uFrMC = 1 if and only if {uTC = 1 (Truncated before MC) and udw = 0 (Fixed after MC)}
          self.model.addConstr(uTC[r,c] >= uFrMC[r,4*c+i], name = 'bFr_MC')
          self.model.addConstr(udw[r,4*c+i] + uFrMC[r,4*c+i] <= 1, name = 'bFr_MC')
          self.model.addConstr(uDiffDetC[r,c] + uTC[r,c] - udw[r,4*c+i] - uFrMC[r,4*c+i] <= 1, name = 'bFr_MC')
          self.model.addConstr(uDiffDetC[r,c] >= uFrMC[r,4*c+i], name = 'bFr_MC')



    '''
    <------------------------------------------------------------------------------------
    ONE_(1-2). Differential propagation in (4).El(<-), (5).Ef(->):
    *PARA.: 
      DX[round,active tag], dx[round, truncated tag]
    *TRAIL: 
      El: DWw[-1,](plaintext) <={ART(stk)}= DXx <={SC}= DYy <={SR}= DZz <={MC}= DWz[Eb,]
      Ef: DWw[-1,](plaintext) ={ART(stk)}=> DXx ={SC}=> DYy ={SR}=> DZz ={MC}=> DWz[Eb,]
    NOTE: different propagation of diff. in dis. and ext.
      1. classical diff. pro. in dis.
      2. equivalent diff. pro. in ext.
    ------------------------------------------------------------------------------------>
    '''
    start_round_l = self.round_Eb
    end_round_l = start_round_l + self.round_Dis + self.round_Ef

    # PARA. needed in dis.
    lDW  = self.model.addVars(range(start_round_l - 1, end_round_l - self.round_Ef), 16, vtype = GRB.BINARY, name = 'lDW')
    ldw  = self.model.addVars(range(start_round_l - 1, end_round_l - self.round_Ef), 16, vtype = GRB.BINARY, name = 'ldw')
    lDX  = self.model.addVars(range(start_round_l, end_round_l - self.round_Ef + 1),             16, vtype = GRB.BINARY, name = 'lDX')
    ldx  = self.model.addVars(range(start_round_l, end_round_l - self.round_Ef + 1),             16, vtype = GRB.BINARY, name = 'ldx')
    lDY  = self.model.addVars(range(start_round_l, end_round_l - self.round_Ef), 16, vtype = GRB.BINARY, name = 'lDY')
    ldy  = self.model.addVars(range(start_round_l, end_round_l - self.round_Ef), 16, vtype = GRB.BINARY, name = 'ldy')
    lDZ  = self.model.addVars(range(start_round_l, end_round_l - self.round_Ef), 16, vtype = GRB.BINARY, name = 'lDZ')
    ldz  = self.model.addVars(range(start_round_l, end_round_l - self.round_Ef), 16, vtype = GRB.BINARY, name = 'ldz')
    lstk = self.model.addVars(range(start_round_l - 1, end_round_l + 1),                             16, vtype = GRB.BINARY, name = 'lstk') # All lstk are needed
    lcan = self.model.addVars(range(start_round_l, end_round_l - self.round_Ef + 1),            16, vtype = GRB.BINARY, name = 'lcan')
    lAC  = self.model.addVars(range(start_round_l, end_round_l - 1),              4, vtype = GRB.BINARY, name = 'lAC') # counting all Active
    lTC  = self.model.addVars(range(start_round_l, end_round_l - 1),              4, vtype = GRB.BINARY, name = 'lTC') # counting all Truncated

    '''========================In distinguisher.(start_round_l, end_round_l - self.round_Ef + self.round_Em (ATK & SC))============================='''

    '''Basic Constrs:'''
    # Remove dx[r,i] > DX[r,i]
    for i in range(16):
      for r in range(start_round_l, end_round_l - self.round_Ef):
        self.model.addConstr(lDX[r+1,i] >= ldx[r+1,i], name='lBas')
        self.model.addConstr(lDY[r,i] >= ldy[r,i], name='lBas')
        self.model.addConstr(lDZ[r,i] >= ldz[r,i], name='lBas')
        self.model.addConstr(lDW[r,i] >= ldw[r,i], name='lBas')

    ''' Operation: ATK '''
    # In Distingusher: [DW,dw]r-1 <=(ustk, ucan)= [DX,dx]
    for r in range(start_round_l - 1, end_round_l - self.round_Ef):
      for i in range(16):
        self.model.addConstr(lDW[r,i] - ldx[r+1,i] - lcan[r+1,i] >= 0, name='lATK')
        self.model.addConstr(- lDW[r,i] + lDX[r+1,i] + lcan[r+1,i] >= 0, name='lATK')
        self.model.addConstr(lDW[r,i] + lstk[r+1,i] - lDX[r+1,i] - 2 * lcan[r+1,i] >= 0, name='lATK')
        self.model.addConstr(- lstk[r+1,i] + lDX[r+1,i] + lcan[r+1,i] >= 0, name='lATK')
        self.model.addConstr(ldw[r,i] == ldx[r+1,i], name='lATK')

    ''' Operation: SC '''
    # NOTE: The connecting operation of Ex. and Dis. is SC
    # In distinguisher [DX,dx] <= [DY,dy]
    for r in range(start_round_l, end_round_l - self.round_Ef): 
      for i in range(16):
        self.model.addConstr(lDX[r,i] == lDY[r,i], name='lSC')
        self.model.addConstr(ldx[r,i] == lDY[r,i], name='lSC')

    ''' Operation: SR '''
    # [DY,dy] <=> (DZ,dz)
    for r in range(start_round_l, end_round_l - self.round_Ef):
      for i in range(16):
        self.model.addConstr(lDY[r,self.SRpermutation_rev[i]] == lDZ[r,i], name='lSC')
        self.model.addConstr(ldy[r,self.SRpermutation_rev[i]] == ldz[r,i], name='lSC')

    ''' Operation: MC '''
      # In Distinguisher [DZ,dz] => [DW,dw]
    for r in range(start_round_l, end_round_l - self.round_Ef):
      for j in [0,4,8,12]:
        self.model.addConstr(lTC[r,j/4] >= ldw[r,j+0], name='lMC')
        self.model.addConstr(lTC[r,j/4] >= ldw[r,j+1], name='lMC')
        self.model.addConstr(lTC[r,j/4] >= ldw[r,j+2], name='lMC')
        self.model.addConstr(lTC[r,j/4] >= ldw[r,j+3], name='lMC')
        self.model.addConstr(lTC[r,j/4] <= sum(ldw[r,j+i] for i in range(4)), name='lMC')
        self.model.addConstr(lTC[r,j/4] <= ldz[r,j+0], name='lMC')
        self.model.addConstr(lTC[r,j/4] <= ldz[r,j+1], name='lMC')
        self.model.addConstr(lTC[r,j/4] <= ldz[r,j+2], name='lMC')
        self.model.addConstr(lTC[r,j/4] <= ldz[r,j+3], name='lMC')
        self.model.addConstr(4*lTC[r,j/4] >= sum(ldz[r,j+i] for i in range(4)), name='lMC')

        self.model.addConstr(lAC[r,j/4] >= lDW[r,j+0], name='lMC')
        self.model.addConstr(lAC[r,j/4] >= lDW[r,j+1], name='lMC')
        self.model.addConstr(lAC[r,j/4] >= lDW[r,j+2], name='lMC')
        self.model.addConstr(lAC[r,j/4] >= lDW[r,j+3], name='lMC')
        self.model.addConstr(lAC[r,j/4] <= sum(lDW[r,j+i] for i in range(4)), name='lMC')
        self.model.addConstr(lDW[r,j+0] + lDW[r,j+1] + lDW[r,j+2] + lDW[r,j+3] + 
                             lDZ[r,j+0] + lDZ[r,j+1] + lDZ[r,j+2] + lDZ[r,j+3] >= 5*lAC[r,j/4], name='lMC')
        self.model.addConstr(lDW[r,j+0] + lDW[r,j+1] + lDW[r,j+2] + lDW[r,j+3] + 
                             lDZ[r,j+0] + lDZ[r,j+1] + lDZ[r,j+2] + lDZ[r,j+3] <= 8*lAC[r,j/4], name='lMC')

    '''
    ==========================================================================================================
    *2 parts for equivalent stk:
      PART 1: leqk = SR^-1(MC^-1(lstk)) in extension
      PART 2: differentail propagation in extension (eqX -(SC)-> eqY -(ART)-> eqZ -(SR)-> eqW -(MC)-> eqX_r+1)
    ==========================================================================================================
    ''' 
    leqDW  = self.model.addVars(range(end_round_l - self.round_Ef, end_round_l - 1),     16, vtype = GRB.BINARY, name = 'leqDW')
    leqdw  = self.model.addVars(range(end_round_l - self.round_Ef, end_round_l - 1),     16, vtype = GRB.BINARY, name = 'leqdw')
    leqDX  = self.model.addVars(range(end_round_l - self.round_Ef + 1, end_round_l),     16, vtype = GRB.BINARY, name = 'leqDX')
    leqdx  = self.model.addVars(range(end_round_l - self.round_Ef + 1, end_round_l),     16, vtype = GRB.BINARY, name = 'leqdx')
    leqDY  = self.model.addVars(range(end_round_l - self.round_Ef, end_round_l),         16, vtype = GRB.BINARY, name = 'leqDY')
    leqdy  = self.model.addVars(range(end_round_l - self.round_Ef, end_round_l),         16, vtype = GRB.BINARY, name = 'leqdy')
    leqDZ  = self.model.addVars(range(end_round_l - self.round_Ef, end_round_l),         16, vtype = GRB.BINARY, name = 'leqDZ')
    leqdz  = self.model.addVars(range(end_round_l - self.round_Ef, end_round_l),         16, vtype = GRB.BINARY, name = 'leqdz')
    leqT   = self.model.addVars(range(end_round_l - self.round_Ef + 1, end_round_l + 1),  4, vtype = GRB.BINARY, name = 'leqT') # for MC^-1 (leqk)
    # lstk have been transformed to leqk
    leqk = self.model.addVars(range(end_round_l - self.round_Ef + 1, end_round_l + 1), 16, vtype = GRB.BINARY, name = 'leqk')
    # No cancellation in extension
    # Reuse lAC and lTC from dis

    '''Basic Constrs:'''
    # Remove dx[r,i] > DX[r,i]
    for i in range(16):
      for r in range(end_round_l - self.round_Ef, end_round_l):
        self.model.addConstr(leqDY[r,i] >= leqdy[r,i], name='lBaseq')
        self.model.addConstr(leqDZ[r,i] >= leqdz[r,i], name='lBaseq')
      for r in range(end_round_l - self.round_Ef, end_round_l - 1):
        self.model.addConstr(leqDX[r+1,i] >= leqdx[r+1,i], name='lBaseq')
        self.model.addConstr(leqDW[r,i] >= leqdw[r,i], name='lBaseq')

    '''
    --------------------------------------------------------------------------------------------------------
    PART 1: Transform stk to leqk (leqk = SR^-1(MC^-1(lstk)) in extension)
    --------------------------------------------------------------------------------------------------------
    '''
    # leqCan = self.model.addVars(range(end_round_l - self.round_Ef + 1, end_round_l + 1),  4, vtype = GRB.INTEGER, name = 'leqCan')

    for r in range(end_round_l - self.round_Ef + 1, end_round_l + 1):
      for c in range(4):
        self.model.addConstr(leqT[r,c] >= lstk[r,4*c+0], name='l2eq')
        self.model.addConstr(leqT[r,c] >= lstk[r,4*c+1], name='l2eq')
        self.model.addConstr(leqT[r,c] >= lstk[r,4*c+2], name='l2eq')
        self.model.addConstr(leqT[r,c] >= lstk[r,4*c+3], name='l2eq')
        self.model.addConstr(leqT[r,c] <= sum(lstk[r,4*c+i] for i in range(4)), name='l2eq')
        self.model.addConstr(leqk[r,self.SRpermutation_rev[4*c+0]] + leqk[r,self.SRpermutation_rev[4*c+1]] + 
                             leqk[r,self.SRpermutation_rev[4*c+2]] + leqk[r,self.SRpermutation_rev[4*c+3]] + 
                             lstk[r,4*c+0] + lstk[r,4*c+1] + lstk[r,4*c+2] + lstk[r,4*c+3] >= 5*leqT[r,c], name='l2eq')
        self.model.addConstr(leqk[r,self.SRpermutation_rev[4*c+0]] + leqk[r,self.SRpermutation_rev[4*c+1]] + 
                             leqk[r,self.SRpermutation_rev[4*c+2]] + leqk[r,self.SRpermutation_rev[4*c+3]] + 
                             lstk[r,4*c+0] + lstk[r,4*c+1] + lstk[r,4*c+2] + lstk[r,4*c+3] <= 8*leqT[r,c], name='l2eq')

        # self.model.addConstr(sum(leqk[r,self.SRpermutation_rev[4*c+i]] for i in range(4)) + leqCan[r,c] == 4*leqT[r,c])

    '''
    --------------------------------------------------------------------------------------------------------
    PART 2: differentail propagation in extension (eqX -(SC)-> eqY -(ATK)-> eqZ -(SR)-> eqW -(MC)-> eqX_r+1)
    --------------------------------------------------------------------------------------------------------
    '''

    ''' Operation: SC '''
    # Connection of dis. and ext.
    for i in range(16):
      self.model.addConstr(leqDY[start_round_l + self.round_Dis,i] == lDX[start_round_l + self.round_Dis,i], name='lfSC')
      self.model.addConstr(leqdy[start_round_l + self.round_Dis,i] == lDX[start_round_l + self.round_Dis,i], name='lfSC')
    # In extension [DX,dx] => [DY,dy]
    for r in range(start_round_l + self.round_Dis + 1, end_round_l):
      for i in range(16):
        self.model.addConstr(leqDY[r,i] == leqDX[r,i], name='lfSC')
        self.model.addConstr(leqdy[r,i] == leqDX[r,i], name='lfSC')

    ''' Operation: eqATK '''
    # In Extension:    [leqDY,leqdy] =(leqk)=> [leqDZ,leqdz]
    # NOTE: No equivalent key.
    for r in range(start_round_l + self.round_Dis, end_round_l):
      for i in range(16):
        self.model.addConstr(leqdy[r,i] == leqdz[r,i], name='lfATK')
        self.model.addConstr(leqk[r+1,i] - leqDZ[r,i] + leqdz[r,i] >= 0, name='lfATK')
        self.model.addConstr(leqDZ[r,i] >= leqk[r+1,i], name='lfATK')

    ''' Operation: eqSR '''
    # [DY,dy] <=> (DZ,dz)  end at the last eqATK
    for r in range(start_round_l + self.round_Dis, end_round_l - 1):
      for i in range(16):
        self.model.addConstr(leqDZ[r,self.SRpermutation_rev[i]] == leqDW[r,i], name='lfSR')
        self.model.addConstr(leqdz[r,self.SRpermutation_rev[i]] == leqdw[r,i], name='lfSR')

    ''' Operation: eqMC '''
    # In Extension [DW,dw] <= [DZ,dz]  end at the last eqATK
    for r in range(start_round_l + self.round_Dis, end_round_l - 1):
      for j in range(4):
        self.model.addConstr(lTC[r,j] >= leqdw[r,4*j+0], name='lfMC')
        self.model.addConstr(lTC[r,j] >= leqdw[r,4*j+1], name='lfMC')
        self.model.addConstr(lTC[r,j] >= leqdw[r,4*j+2], name='lfMC')
        self.model.addConstr(lTC[r,j] >= leqdw[r,4*j+3], name='lfMC')
        self.model.addConstr(lTC[r,j] <= sum(leqdw[r,4*j+i] for i in range(4)), name='lfMC')
        self.model.addConstr(lTC[r,j] <= leqdx[r+1,4*j+0], name='lfMC')
        self.model.addConstr(lTC[r,j] <= leqdx[r+1,4*j+1], name='lfMC')
        self.model.addConstr(lTC[r,j] <= leqdx[r+1,4*j+2], name='lfMC')
        self.model.addConstr(lTC[r,j] <= leqdx[r+1,4*j+3], name='lfMC')
        self.model.addConstr(4*lTC[r,j] >= sum(leqdx[r+1,4*j+i] for i in range(4)), name='lfMC')

        self.model.addConstr(lAC[r,j] >= leqDW[r,4*j+0], name='lfMC')
        self.model.addConstr(lAC[r,j] >= leqDW[r,4*j+1], name='lfMC')
        self.model.addConstr(lAC[r,j] >= leqDW[r,4*j+2], name='lfMC')
        self.model.addConstr(lAC[r,j] >= leqDW[r,4*j+3], name='lfMC')
        self.model.addConstr(lAC[r,j] <= sum(leqDW[r,4*j+i] for i in range(4)), name='lfMC')
        self.model.addConstr(leqDW[r,4*j+0] + leqDW[r,4*j+1] + leqDW[r,4*j+2] + leqDW[r,4*j+3] + 
                             leqDX[r+1,4*j+0] + leqDX[r+1,4*j+1] + leqDX[r+1,4*j+2] + leqDX[r+1,4*j+3] >= 5*lAC[r,j], name='lfMC')
        self.model.addConstr(leqDW[r,4*j+0] + leqDW[r,4*j+1] + leqDW[r,4*j+2] + leqDW[r,4*j+3] + 
                             leqDX[r+1,4*j+0] + leqDX[r+1,4*j+1] + leqDX[r+1,4*j+2] + leqDX[r+1,4*j+3] <= 8*lAC[r,j], name='lfMC')
        
    '''
    -------------------------------------------------------------------------------
    TWO: Key schedule in Em+Ef+Em (No direction)
    *PARA.: 
      lLANE:         A LANE of permutation of key schedule (for Type1 cancellation)
      lT2Can(A/B/C): Auxiliary variables to count Type2 cancellation
      lType2:        Type2 cancellation
    -------------------------------------------------------------------------------
    '''

    lLANE   = self.model.addVars(16, vtype = GRB.BINARY, name = 'lLANE')
    lT2CanA = self.model.addVars(range(start_round_l - 1, end_round_l - self.round_Ef), 4, lb = 0, vtype = GRB.INTEGER, name = 'lT2CanA')
    lT2CanB = self.model.addVars(range(start_round_l - 1, end_round_l - self.round_Ef), 4, lb = 0, vtype = GRB.INTEGER, name = 'lT2CanB')
    lT2CanC = self.model.addVars(range(start_round_l, end_round_l - self.round_Ef + 1), 4, lb = 0, vtype = GRB.INTEGER, name = 'lT2CanC')
    lType1  = self.model.addVars(16, lb = 0, vtype = GRB.INTEGER, name = 'lType1')
    lType2  = self.model.addVars(range(start_round_l - 1, end_round_l - self.round_Ef), 4, lb = 0, vtype = GRB.INTEGER, name = 'lType2')

    self.model.addConstr(sum(lLANE[i] for i in range(16)) >= 1)

    ''' Mark LANE '''
    for i in range(16):
      for r in range(start_round_l - 1, end_round_l + 1):
        self.model.addConstr(lLANE[i] >= lstk[r, self.hTable[i][r]], name='lMLANE')

    ''' Type 1 cancellation '''
    for i in range(16):
      self.model.addConstr(lType1[i] == (end_round_l + 1 - (start_round_l - 1)) * lLANE[i] - 
                           sum(lstk[r, self.hTable[i][r]] for r in range(start_round_l - 1, end_round_l + 1)), name = 'lType1C')
      self.model.addConstr(lType1[i] <= self.s - 1, name = 'lType1C') # Cancellation no more than 2 in each position (for Deoxys-BC-384)
    
    ''' Type 2 cancellation (In Eu and Em) '''
    for r in range(start_round_l, end_round_l - self.round_Ef):
      for c in range(4):
        # T2CanA: Active bytes in one column, before MC
        self.model.addConstr(lT2CanA[r,c] == lDZ[r,4*c+0] + lDZ[r,4*c+1] + lDZ[r,4*c+2] + lDZ[r,4*c+3], name = 'lType2A')
        # T2CanB: Inactive bytes in one column, after MC
        self.model.addConstr(lT2CanB[r,c] == 4*lAC[r,c] - lDW[r,4*c+0] - lDW[r,4*c+1] - lDW[r,4*c+2] - lDW[r,4*c+3], name = 'lType1B')
        # T2CanC: Cancellation appears in AK, in the next round
        self.model.addConstr(lT2CanC[r+1,c] == lcan[r+1,4*c+0] + lcan[r+1,4*c+1] + lcan[r+1,4*c+2] + lcan[r+1,4*c+3], name = 'lType2C')

    # Type2 = -(T2CanA - T2CanB - T2CanC)
    for r in range(start_round_l - 1, end_round_l - self.round_Ef):
      for c in range(4):
        self.model.addConstr(lType2[r,c] >= (lT2CanB[r,c] + lT2CanC[r+1,c] - lT2CanA[r,c]), name = 'lType2CAN')

    # ''' Counting of all cancellations'''
    # NOTE: s*LANE[0~15] >= r*LANE[0~15]-stk[0~r,0~15] + Type2Can. (Type2Can. in dist. only, edges)
    self.model.addConstr(self.s * sum(lLANE[i] for i in range(16)) -   # s * LANE
                         sum(lType1[i] for i in range(16)) - 
                         sum(lType2[r,c] for r in range(start_round_l - 1, end_round_l - self.round_Ef) for c in range(4))
                         >= 1, name = 'lCAN')

    # To mark the Type1 cancellation in Ef
    for r in range(start_round_l + self.round_Dis + 1, end_round_l):
      for i in range(16):
        self.model.addConstr(leqDX[r, i] >= lstk[r, i])

    '''
    -----------------------------------------------------------------------------------------------
    THREE: Guess-and-Determine (Ef)
      PRAR.: 
        lGstk: Whether the stk is guessed or not
        lDetW[-1,](plaintext) ={ART(uGstk)}=> uDetX ={SC}=> uDetY ={SR}=> uDetZ ={MC}=> uDetW [Eb,]
        Auxiliary variable of MC: lDetC, lDiffDetMC, lDiffDetC
        Filter tag: lFrSB, lFrMC
      2 PARTS:
        1. Propagation of determine
        2. Obtaining filters
    -----------------------------------------------------------------------------------------------
    ''' 
    lGstk   = self.model.addVars(range(end_round_l - self.round_Ef + 1, end_round_l + 1), 16, vtype = GRB.BINARY, name = 'lGstk')
    ldetEQX = self.model.addVars(range(end_round_l - self.round_Ef, end_round_l),         16, vtype = GRB.BINARY, name = 'ldetEQX')
    ldetEQY = self.model.addVars(range(end_round_l - self.round_Ef, end_round_l),         16, vtype = GRB.BINARY, name = 'ldetEQY')
    ldetEQZ = self.model.addVars(range(end_round_l - self.round_Ef, end_round_l),         16, vtype = GRB.BINARY, name = 'ldetEQZ')
    ldetEQW = self.model.addVars(range(end_round_l - self.round_Ef, end_round_l - 1),     16, vtype = GRB.BINARY, name = 'ldetEQW')
    lDetC   = self.model.addVars(range(end_round_l - self.round_Ef + 1, end_round_l),          4, vtype = GRB.BINARY, name = 'lDetC') # mark 'det.' in one column
    lFrSB   = self.model.addVars(range(end_round_l - self.round_Ef, end_round_l),         16, vtype = GRB.BINARY, name = 'lFrSB')
    lDiffDetMC = self.model.addVars(range(end_round_l - self.round_Ef, end_round_l - 1),      16, vtype = GRB.BINARY, name = 'lDiffDetMC') # mark diff. det. of cell
    lDiffDetC = self.model.addVars(range(end_round_l - self.round_Ef, end_round_l - 1),        4, vtype = GRB.BINARY, name = 'lDiffDetC') # mark diff. det. in one colume
    lFrMC = self.model.addVars(range(end_round_l - self.round_Ef, end_round_l - 1),           16, vtype = GRB.BINARY, name = 'lFrMC')
    
    '''
    ---------------------------------
    PART 1.  Propagation of determine
    ---------------------------------
    '''

    ''' KR: Guess round key and determine '''
    # lDetY <-(lGstk)- lDetZ
    for r in range(end_round_l - self.round_Ef, end_round_l):
      for i in range(16):
        self.model.addConstr(lGstk[r+1,i] >= ldetEQY[r,i], name = 'lDetAK')
        self.model.addConstr(ldetEQZ[r,i] >= ldetEQY[r,i], name = 'lDetAK')
        self.model.addConstr(- ldetEQZ[r,i] - lGstk[r+1,i] + ldetEQY[r,i] >= -1, name = 'lDetAK')

    ''' KR: Determine in SC '''
    for r in range(end_round_l - self.round_Ef, end_round_l):
      for i in range(16):
        self.model.addConstr(ldetEQX[r,i] == ldetEQY[r,i], name = 'lDetSC') 

    ''' KR: 'Det.' propagation in MC '''
    for r in range(end_round_l - self.round_Ef + 1, end_round_l):
      for c in range(4):
        # NOTE: lDetC = 1 when all ldetEQX = 1, else lDetC = 0 
        self.model.addConstr(lDetC[r,c] <= ldetEQX[r,4*c+0], name = 'lDetMC')
        self.model.addConstr(lDetC[r,c] <= ldetEQX[r,4*c+1], name = 'lDetMC')
        self.model.addConstr(lDetC[r,c] <= ldetEQX[r,4*c+2], name = 'lDetMC')
        self.model.addConstr(lDetC[r,c] <= ldetEQX[r,4*c+3], name = 'lDetMC')
        self.model.addConstr(ldetEQX[r,4*c+0] + ldetEQX[r,4*c+1] + ldetEQX[r,4*c+2] + ldetEQX[r,4*c+3] - lDetC[r,c] <= 3, name = 'lDetMC')
        # NOTE: all lDetZ = 1 when lDetC = 1, else all lDetZ = 0
        self.model.addConstr(lDetC[r,c] == ldetEQW[r-1,4*c+0], name = 'lDetMC')
        self.model.addConstr(lDetC[r,c] == ldetEQW[r-1,4*c+1], name = 'lDetMC')
        self.model.addConstr(lDetC[r,c] == ldetEQW[r-1,4*c+2], name = 'lDetMC')
        self.model.addConstr(lDetC[r,c] == ldetEQW[r-1,4*c+3], name = 'lDetMC')

    ''' KR: 'Det.' propagation in SR '''
    for r in range(end_round_l - self.round_Ef, end_round_l - 1):
      for i in range(16):
        self.model.addConstr(ldetEQZ[r,self.SRpermutation_rev[i]] == ldetEQW[r,i], name = 'lDetSR')

    '''
    ---------------------------------
    PART 2.  Obtaining filters
    ---------------------------------
    '''
    ''' KR: Filter in SC '''
    # At the point of connection
    for i in range(16):
      temp_r = end_round_l - self.round_Ef
      self.model.addConstr(lDX[temp_r,i] - ldx[temp_r,i] - lFrSB[temp_r,i] >= 0, name = 'lFrSC')
      self.model.addConstr(ldetEQX[temp_r,i] >= lFrSB[temp_r,i], name = 'lFrSC')
      self.model.addConstr(- ldetEQX[temp_r,i] - lDX[temp_r,i] + ldx[temp_r,i] + lFrSB[temp_r,i] >= -1, name = 'lFrSC')
    # In the extension
    for r in range(end_round_l - self.round_Ef + 1, end_round_l):
      for i in range(16):
        self.model.addConstr(leqDX[r,i] - leqdx[r,i] - lFrSB[r,i] >= 0, name = 'lFrSC')
        self.model.addConstr(ldetEQX[r,i] >= lFrSB[r,i], name = 'lFrSC')
        self.model.addConstr(- ldetEQX[r,i] - leqDX[r,i] + leqdx[r,i] + lFrSB[r,i] >= -1, name = 'lFrSC')

    ''' KR: Filter obtain from MC (for cell -> column) '''
    # KR: Determine whether the difference before MC is determinable '''
    for r in range(end_round_l - self.round_Ef, end_round_l - 1):
      for i in range(16):
        # NOTE: lDiffDetMC = ldetEQX when lDW = 1 (ACTIVE), else uDiffDetMC = 1
        self.model.addConstr(lDiffDetMC[r,i] >= ldetEQX[r+1,i], name = 'lFrMC')
        self.model.addConstr(ldetEQX[r+1,i] - leqDX[r+1,i] - lDiffDetMC[r,i] >= -1, name = 'lFrMC')
        self.model.addConstr(leqDX[r+1,i] + lDiffDetMC[r,i] >= 1, name = 'lFrMC')
      for c in range(4):
        # NOTE: lDiffDetC = 1 when all lDiffDetMC = 1, else lDiffDetC = 0 
        self.model.addConstr(lDiffDetC[r,c] <= lDiffDetMC[r,4*c+0], name = 'lFrMC')
        self.model.addConstr(lDiffDetC[r,c] <= lDiffDetMC[r,4*c+1], name = 'lFrMC')
        self.model.addConstr(lDiffDetC[r,c] <= lDiffDetMC[r,4*c+2], name = 'lFrMC')
        self.model.addConstr(lDiffDetC[r,c] <= lDiffDetMC[r,4*c+3], name = 'lFrMC')
        self.model.addConstr(lDiffDetMC[r,4*c+0] + lDiffDetMC[r,4*c+1] + lDiffDetMC[r,4*c+2] + lDiffDetMC[r,4*c+3] 
                             - lDiffDetC[r,c] <= 3, name = 'lFrMC')
        
    # Filter obtaining from eqMC
    for r in range(end_round_l - self.round_Ef, end_round_l - 1):
      for c in range(4):
        for i in range(4):
          # NOTE: lFrMC = 0 when lDiffDetC = 0, 
          # NOTE: lFrMC = 1 if and only if {lTC = 1 (Truncated before MC) and leqdw = 0 (Fixed after MC)}
          self.model.addConstr(lTC[r,c] >= lFrMC[r,4*c+i], name = 'lFreqMC')
          self.model.addConstr(leqdw[r,4*c+i] + lFrMC[r,4*c+i] <= 1, name = 'lFreqMC')
          self.model.addConstr(lDiffDetC[r,c] + lTC[r,c] - leqdw[r,4*c+i] - lFrMC[r,4*c+i] <= 1, name = 'lFreqMC')
          self.model.addConstr(lDiffDetC[r,c] >= lFrMC[r,4*c+i], name = 'lFreqMC')

































    '''
    -----------------------------------------------------------------------------------------------
    FOUR: Miss-In-The-Middle
      A single S-box layer in the middle round (1-r, BCT)
      At last one FIXED appear in in/out-put of BCT to get constraints.
    -----------------------------------------------------------------------------------------------
    ''' 
    self.model.addConstr(sum((uDX[r,i] - udx[r,i]) * 
                             (lDY[r,i] - ldy[r,i]) for i in range(16) for r in range(self.round_Eb, self.round_Eb + self.round_Dis)) >= 1)





# TEST --------------------------------------------------------------
#     pacp = [(0,0),(1,1),(0,0),(0,0),
#             (0,0),(0,0),(0,0),(1,0),
#             (1,0),(0,0),(0,0),(0,0),
#             (1,1),(1,0),(0,0),(0,0)]
#     for i in range(16):
#       self.model.addConstr(uDW[-1,i] == pacp[i][0])
#       self.model.addConstr(udw[-1,i] == pacp[i][1])

#     pacx4 = [(1,1),(1,1),(1,1),(1,1),
#              (0,0),(1,0),(0,0),(0,0),
#              (1,1),(1,1),(1,1),(1,1),
#              (1,1),(1,1),(1,1),(1,1)]
#     for i in range(16):
#       self.model.addConstr(uDX[4,i] == pacx4[i][0])
#       self.model.addConstr(udx[4,i] == pacx4[i][1])

#     # self.model.addConstr(sum(ustk[2,i] for i in range(16)) == 0)


#     pacC = [(1,1),(1,1),(1,1),(1,1),
#             (1,1),(1,1),(1,1),(1,1),
#             (0,0),(0,0),(0,0),(1,1),
#             (0,0),(0,0),(0,0),(1,0)]
#     for i in range(16):
#       self.model.addConstr(leqDZ[9,i] == pacC[i][0])
#       self.model.addConstr(leqdz[9,i] == pacC[i][1])

#     pacZ4 = [(1,0),(1,0),(1,0),(1,0),
#              (1,1),(1,1),(1,1),(1,1),
#              (1,1),(1,1),(1,1),(1,1),
#              (1,1),(1,1),(1,1),(1,1)]
#     for i in range(16):
#       self.model.addConstr(lDZ[4,i] == pacZ4[i][0])
#       self.model.addConstr(ldz[4,i] == pacZ4[i][1])

# TEST --------------------------------------------------------------

    '''
    -----------------------------------------------------------------------------------------------
    FIVE: Complexities
      PARAS: 
        Eb: rb, mb, mbp, cbp
        Ef: rf, mf, mfp, cfp
          rb/rf:   truncated bits of Eb/Ef
          cb/cf:   rb - dim(DeltaX of the head/tail of distinguisher)
          mb/mf:   involved key of Eb/Ef
          mbp/mfp: guessed key of Eb/Ef
          cbp/cfp: determined key of Eb/Ef
        Time Complexity: T0, T1, T2, T31, T32, T4, epsilon (log2)
          T0:  data collection
          T1:  partial encryption & decryption
          T2:  pairs generation
          T31: quartets generation
          T32: quartets processing (and key info. extraction)
          T4:  exhaustive searching
          epsilon: estimated complexity para. of T32
      2 PARTs:
        1. Parameters definition (in this model)
        2. Time complexity computing
    -----------------------------------------------------------------------------------------------
    '''
    '''
    Choosen plaintext/ciphertext
    '''
    for i in range(16):
      self.model.addConstr(uDetW[-1,i] == 1)
      self.model.addConstr(ldetEQZ[end_round_l - 1,i] == 1)

    rb  = self.model.addVar(vtype = GRB.INTEGER, name = 'rb')
    rf  = self.model.addVar(vtype = GRB.INTEGER, name = 'rf')
    cb  = self.model.addVar(vtype = GRB.INTEGER, name = 'cb')
    cf  = self.model.addVar(vtype = GRB.INTEGER, name = 'cf')
    mb  = self.model.addVar(vtype = GRB.INTEGER, name = 'mb')
    mf  = self.model.addVar(vtype = GRB.INTEGER, name = 'mf')
    mbp = self.model.addVar(vtype = GRB.INTEGER, name = 'mb\'')
    mfp = self.model.addVar(vtype = GRB.INTEGER, name = 'mf\'')
    cbp = self.model.addVar(vtype = GRB.INTEGER, name = 'cb\'')
    cfp = self.model.addVar(vtype = GRB.INTEGER, name = 'cf\'')

    Dc  = self.model.addVar(lb = 0, vtype = GRB.CONTINUOUS, name = 'Dc')
    Qc  = self.model.addVar(lb = 0, vtype = GRB.CONTINUOUS, name = 'Qc')
    T0  = self.model.addVar(lb = 0, vtype = GRB.CONTINUOUS, name = 'T0')
    T1  = self.model.addVar(lb = 0, vtype = GRB.CONTINUOUS, name = 'T1')


    T2  = self.model.addVar(lb = 0, vtype = GRB.CONTINUOUS, name = 'T2')



    T31 = self.model.addVar(lb = 0, vtype = GRB.CONTINUOUS, name = 'T31')
    T32 = self.model.addVar(lb = 0, vtype = GRB.CONTINUOUS, name = 'T32')
    T4  = self.model.addVar(lb = 0, vtype = GRB.CONTINUOUS, name = 'T4')
    Tc  = self.model.addVar(lb = 0, vtype = GRB.CONTINUOUS, name = 'Tc')
    eps = self.model.addVar(lb = 0, vtype = GRB.CONTINUOUS, name = 'epsilon') # May not be used

    '''
    ------------------------------------------------------------
    PART 1. Parameters definition (in this model)
      Eb: rb, mb, mbp, cbp
      Ef: rf, mf, mfp, cfp
      rb/rf:   truncated bits of Eb/Ef
      cb/cf:   rb - dim(DeltaX of the head/tail of distinguisher)
      mb/mf:   involved key of Eb/Ef
      mbp/mfp: guessed key of Eb/Ef
      cbp/cfp: determined key of Eb/Ef
    ------------------------------------------------------------
    '''
    ''' Defination of PARA. in Eb'''
    self.model.addConstr(rb  == self.cell_size * sum(udw[-1,i] for i in range(16)))
    self.model.addConstr(cb  == rb - self.cell_size * sum(udx[self.round_Eb,i] for i in range(16)))
    self.model.addConstr(mb  == self.cell_size * sum(udx[r,i] for r in range(self.round_Eb) for i in range(16)))
    self.model.addConstr(mbp == self.cell_size * sum(uGstk[r,i] for r in range(self.round_Eb) for i in range(16)))
    self.model.addConstr(cbp == self.cell_size * sum(uFrSB[r,i] + uFrMC[r,i] for r in range(self.round_Eb) for i in range(16)))

    ''' Defination of PARA. in Ef'''
    self.model.addConstr(rf  == self.cell_size * sum(leqdz[end_round_l - 1,i] for i in range(16)))
    self.model.addConstr(cf  == rf - self.cell_size * sum(ldx[end_round_l - self.round_Ef,i] for i in range(16)))
    self.model.addConstr(mf  == self.cell_size * sum(leqdz[r,i] for r in range(end_round_l - self.round_Ef, end_round_l) for i in range(16)))
    self.model.addConstr(mfp == self.cell_size * sum(lGstk[r,i] for r in range(end_round_l - self.round_Ef + 1, end_round_l + 1) for i in range(16)))
    self.model.addConstr(cfp == self.cell_size * sum(lFrSB[r,i] for r in range(end_round_l - self.round_Ef, end_round_l) for i in range(16))
                              + self.cell_size * sum(lFrMC[r,i] for r in range(end_round_l - self.round_Ef, end_round_l - 1) for i in range(16)))

    '''
    -----------------------------------------------------
    2. Time complexity computing
      T0:  data collection
      T1:  partial encryption & decryption
      T2:  pairs generation
      T31: quartets generation
      T32: quartets processing (and key info. extraction)
      T4:  exhaustive searching
      epsilon: estimated complexity para. of T32
    -----------------------------------------------------
    '''
    # Data
    self.model.addConstr(Dc == self.block_size + self.tz/2)
    # Quartet
    self.model.addConstr(Qc == 2 * (cb + cf) + self.tz)

    '''
    Complexity for related key
    '''
    # T0 (T0 = D · 4)
    self.model.addConstr(T0 == Dc + 2)
    # T1 (T1 = 2^{mb'+mf'} · D · 4)
    self.model.addConstr(T1 == (mbp + mfp) + Dc + 2)
    # T2 (T2 = 2^{mb'+mf'} · D · min[ 2^{rb-cb'}, D · 2^{rf-cf'-n} ] · 2)
    if self.pgP == True:
      self.model.addConstr(T2 == (mbp + mfp) + Dc + rb - cbp + 1)
    else:
      self.model.addConstr(T2 == (mbp + mfp) + Dc + (Dc + rf - cfp - self.block_size) + 1)
    # T31 (T31 = 2^{mb'+mf'-2cb'-2cf'} · Q)
    self.model.addConstr(T31 == (mbp + mfp - 2*cbp - 2*cfp) + Qc)
    # T32 (T32 = 2^{mb+mf+z} · epsilon)
    self.model.addConstr(T32 == (mb + mf + self.tz) + eps)
    # T4 (T4 = 2^{k-x}, where 2^z = xln2)
    self.model.addConstr(T4 == self.key_size - self.tx)

    # Objective function
    self.model.addConstr(Tc >= T0)
    self.model.addConstr(Tc >= T1)
    self.model.addConstr(Tc >= T2)
    self.model.addConstr(Tc >= T31)
    self.model.addConstr(Tc >= T32)
    self.model.addConstr(Tc >= T4)


    ''' Objective function'''
    self.model.ModelSense = GRB.MINIMIZE
    self.model.setObjectiveN(Tc, index=0, priority=2, name='Tc')
    self.model.setObjectiveN(T32, index=1, priority=1, name='T32')


    self.model.write(self.name + '.lp')
    self.model.optimize()

    if self.model.Status == GRB.INFEASIBLE:
        print('-'*40 + '\n| Model is infeasible, computing IIS...|\n' + '-'*40)
        self.model.computeIIS()
        self.model.write(self.name + '.ilp')
    self.model.write(self.name + '.sol')

    print('>>>>> Solution <<<<<')
    print('='*90)
    print('|| Choosen Plaintext | ' if self.pgP else('Choosen Chiphertext | '), end='')
    print('PC = {:5} ||'.format(sum(udw[-1,i].x + leqdz[end_round_l-1, i].x for i in range(16))))
    print('-'*90)
    print('|| x  = {:5} | z  = {:5} ||'.format(self.tx, self.tz))
    print('-'*90)
    print("|| rb = {:5} | cb = {:5} | mb = {:5} | mb\' = {:5} | cb\' = {:5} ||".format(rb.x, cb.x, mb.x, mbp.x, cbp.x))
    print("|| rf = {:5} | cf = {:5} | mf = {:5} | mf\' = {:5} | cf\' = {:5} ||".format(rf.x, cf.x, mf.x, mfp.x, cfp.x))
    print('-'*90)
    print('|| D = {:5} | Dc  = {:5} | Q  = {:5} | Qc  = {:5} ||'.format(Dc.x, Dc.x+2, Qc.x, Qc.x-2*cbp.x-2*cfp.x))
    print('-'*90)
    print('|| T0 = {:5} | T1 = {:5} | T2 = {:5} | T31 = {:5} | T32 = {:5} | T4 = {:5} ||'.format(T0.x, T1.x, T2.x, T31.x, T32.x, T4.x))
    print('-'*90)
    print('|| Tc = {:5} ||'.format(Tc.x))
    print('='*90)

  '''
  =====================
  Tools:
    x to z (2^z = xln2)
  =====================
  '''
  def x_to_z(self,x):
    # ln2 = 0.69
    return round(math.log2(0.7 * x), 1)


if __name__ == '__main__':
    
    # IB_Deoxys(key_size, round_Eb, round_Eu, round_Em, round_El, round_Ef, setX, pgP)

    ''' 10-r Deoxys-BC-256 '''
    DeoxysBC256_10r = IB_DandJ(256, 1, 7, 2, 104, True)
    DeoxysBC256_10r.ib_model()

    ''' 11-r Deoxys-BC-256 '''
    DeoxysBC256_11r = IB_DandJ(256, 2, 7, 2, 16, False) # adaptive adjustment to x=20 in this attack
    DeoxysBC256_11r.ib_model()

    ''' 13-r Deoxys-BC-384 '''
    DeoxysBC384_13r = IB_DandJ(384, 2, 9, 2, 160, True)
    DeoxysBC384_13r.ib_model()

    ''' 14-r Deoxys-BC-384 '''
    DeoxysBC384 = IB_DandJ(384, 2, 9, 3, 48, True)
    DeoxysBC384.ib_model()

    # ''' ----------------------------------------------------- '''

    ''' 10-r Joltik-BC-128 '''
    Joltik128_10r = IB_DandJ(128, 1, 7, 2, 52, True)
    Joltik128_10r.ib_model()

    ''' 11-r Joltik-BC-128 '''
    Joltik128_11r = IB_DandJ(128, 2, 7, 2, 8, False)
    Joltik128_11r.ib_model()

    ''' 13-r Joltik-BC-192 '''
    Joltik192_13r = IB_DandJ(192, 2, 9, 2, 86, False)
    Joltik192_13r.ib_model()

    ''' 14-r Joltik-BC-192 '''
    Joltik192_14r = IB_DandJ(192, 2, 9, 3, 24 , True)
    Joltik192_14r.ib_model()
