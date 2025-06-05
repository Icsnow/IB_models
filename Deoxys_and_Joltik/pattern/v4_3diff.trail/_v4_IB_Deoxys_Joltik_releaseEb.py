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
    self.name = './v4_' +  cipher_name + str(self.key_size) + '_' + str(self.round_Eb + self.round_Dis + self.round_Ef) + 'r'
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
    # --------------------------------------------
    # ---------- The first upper trail -----------
    # --------------------------------------------
    end_round_u = self.round_Eb + self.round_Dis
    uDW0  = self.model.addVars(range(-1, end_round_u),          16, vtype = GRB.BINARY, name = 'uDW0')
    udw0  = self.model.addVars(range(-1, end_round_u),          16, vtype = GRB.BINARY, name = 'udw0')# DWw[-1,](plaintext)
    uDX0  = self.model.addVars(end_round_u, 16, vtype = GRB.BINARY, name = 'uDX0')
    udx0  = self.model.addVars(end_round_u, 16, vtype = GRB.BINARY, name = 'udx0')
    uDY0  = self.model.addVars(end_round_u,                     16, vtype = GRB.BINARY, name = 'uDY0')
    udy0  = self.model.addVars(end_round_u,                     16, vtype = GRB.BINARY, name = 'udy0')
    uDZ0  = self.model.addVars(end_round_u,                     16, vtype = GRB.BINARY, name = 'uDZ0')
    udz0  = self.model.addVars(end_round_u,                     16, vtype = GRB.BINARY, name = 'udz0')
    ustk0 = self.model.addVars(end_round_u + 1, 16, vtype = GRB.BINARY, name = 'ustk0')
    ucan0 = self.model.addVars(end_round_u,     16, vtype = GRB.BINARY, name = 'ucan0')
    uAC0  = self.model.addVars(end_round_u,                      4, vtype = GRB.BINARY, name = 'uAC0') # Mark active in one column
    uTC0  = self.model.addVars(end_round_u,                      4, vtype = GRB.BINARY, name = 'uTC0') # Mark truncated in one column

    '''Basic Constrs:'''
    # Remove dx[r,i] > DX[r,i]
    for i in range(16):
      for r in range(-1, end_round_u):
        self.model.addConstr(uDW0[r,i] >= udw0[r,i], name = 'uBas')
      for r in range(end_round_u):
        self.model.addConstr(uDX0[r,i] >= udx0[r,i], name = 'uBas')
        self.model.addConstr(uDY0[r,i] >= udy0[r,i], name = 'uBas')
        self.model.addConstr(uDZ0[r,i] >= udz0[r,i], name = 'uBas')
        
    '''Operation: ATK'''
    for r in range(-1, end_round_u - 1):
      for i in range(16):
        self.model.addConstr(uDW0[r,i] - udx0[r+1,i] - ucan0[r+1,i] >= 0, name = 'uATK')
        self.model.addConstr(- uDW0[r,i] + uDX0[r+1,i] + ucan0[r+1,i] >= 0, name = 'uATK')
        self.model.addConstr(uDW0[r,i] + ustk0[r+1,i] - uDX0[r+1,i] - 2 * ucan0[r+1,i] >= 0, name = 'uATK')
        self.model.addConstr(- ustk0[r+1,i] + uDX0[r+1,i] + ucan0[r+1,i] >= 0, name = 'uATK')
        self.model.addConstr(udw0[r,i] == udx0[r+1,i], name = 'uATK')

    '''Operation: SC'''
    # In Extension     [DX,dx] <= [DY,dy]: [1,1] <= [1,0]/[1,1], [0,0] <= [0,0]
    for r in range(0, self.round_Eb):
      for i in range(16):
        self.model.addConstr(uDX0[r,i] == uDY0[r,i], name = 'bSC')
        self.model.addConstr(udx0[r,i] == uDY0[r,i], name = 'bSC')
    # In distinguisher [DX,dx] => [DY,dy]: [1,0]/[1,1] => [1,1], [0,0] => [0,0]
    for r in range(self.round_Eb, self.round_Eb + self.round_Dis):
      for i in range(16):
        self.model.addConstr(uDY0[r,i] == uDX0[r,i], name = 'uSC')
        self.model.addConstr(udy0[r,i] == uDX0[r,i], name = 'uSC')

    '''Operation: SR''' 
    # [DY,dy] <=> (DZ,dz)
    for r in range(self.round_Eb + self.round_Dis):
      for i in range(16):
        self.model.addConstr(uDY0[r,self.SRpermutation_rev[i]] == uDZ0[r,i], name = 'uSR')
        self.model.addConstr(udy0[r,self.SRpermutation_rev[i]] == udz0[r,i], name = 'uSR')

    '''Operation: MC'''
    # In Extension: [DZ,dz] => [DW,dw]
    for r in range(self.round_Eb):
      for j in [0,4,8,12]:
        self.model.addConstr(uTC0[r,j/4] >= udw0[r,j+0], name = 'bMC')
        self.model.addConstr(uTC0[r,j/4] >= udw0[r,j+1], name = 'bMC')
        self.model.addConstr(uTC0[r,j/4] >= udw0[r,j+2], name = 'bMC')
        self.model.addConstr(uTC0[r,j/4] >= udw0[r,j+3], name = 'bMC')
        self.model.addConstr(uTC0[r,j/4] <= sum(udw0[r,j+i] for i in range(4)), name = 'bMC')
        self.model.addConstr(uTC0[r,j/4] <= udz0[r,j+0], name = 'bMC')
        self.model.addConstr(uTC0[r,j/4] <= udz0[r,j+1], name = 'bMC')
        self.model.addConstr(uTC0[r,j/4] <= udz0[r,j+2], name = 'bMC')
        self.model.addConstr(uTC0[r,j/4] <= udz0[r,j+3], name = 'bMC')
        self.model.addConstr(4*uTC0[r,j/4] >= sum(udz0[r,j+i] for i in range(4)), name = 'bMC')

        self.model.addConstr(uAC0[r,j/4] >= uDW0[r,j+0], name = 'bMC')
        self.model.addConstr(uAC0[r,j/4] >= uDW0[r,j+1], name = 'bMC')
        self.model.addConstr(uAC0[r,j/4] >= uDW0[r,j+2], name = 'bMC')
        self.model.addConstr(uAC0[r,j/4] >= uDW0[r,j+3], name = 'bMC')
        self.model.addConstr(uAC0[r,j/4] <= sum(uDW0[r,j+i] for i in range(4)), name = 'bMC')
        self.model.addConstr(uDW0[r,j+0] + uDW0[r,j+1] + uDW0[r,j+2] + uDW0[r,j+3] + 
                             uDZ0[r,j+0] + uDZ0[r,j+1] + uDZ0[r,j+2] + uDZ0[r,j+3] >= 5*uAC0[r,j/4], name = 'bMC')
        self.model.addConstr(uDW0[r,j+0] + uDW0[r,j+1] + uDW0[r,j+2] + uDW0[r,j+3] + 
                             uDZ0[r,j+0] + uDZ0[r,j+1] + uDZ0[r,j+2] + uDZ0[r,j+3] <= 8*uAC0[r,j/4], name = 'bMC')
    # In Distinguisher [DW,dw] <= [DZ,dz]
    for r in range(self.round_Eb, self.round_Eb + self.round_Dis):
      for j in [0,4,8,12]:
        self.model.addConstr(uTC0[r,j/4] >= udz0[r,j+0], name = 'uMC')
        self.model.addConstr(uTC0[r,j/4] >= udz0[r,j+1], name = 'uMC')
        self.model.addConstr(uTC0[r,j/4] >= udz0[r,j+2], name = 'uMC')
        self.model.addConstr(uTC0[r,j/4] >= udz0[r,j+3], name = 'uMC')
        self.model.addConstr(uTC0[r,j/4] <= sum(udz0[r,j+i] for i in range(4)), name = 'uMC')
        self.model.addConstr(uTC0[r,j/4] <= udw0[r,j+0], name = 'uMC')
        self.model.addConstr(uTC0[r,j/4] <= udw0[r,j+1], name = 'uMC')
        self.model.addConstr(uTC0[r,j/4] <= udw0[r,j+2], name = 'uMC')
        self.model.addConstr(uTC0[r,j/4] <= udw0[r,j+3], name = 'uMC')
        self.model.addConstr(4*uTC0[r,j/4] >= sum(udw0[r,j+i] for i in range(4)), name = 'uMC')

        self.model.addConstr(uAC0[r,j/4] >= uDZ0[r,j+0], name = 'uMC')
        self.model.addConstr(uAC0[r,j/4] >= uDZ0[r,j+1], name = 'uMC')
        self.model.addConstr(uAC0[r,j/4] >= uDZ0[r,j+2], name = 'uMC')
        self.model.addConstr(uAC0[r,j/4] >= uDZ0[r,j+3], name = 'uMC')
        self.model.addConstr(uAC0[r,j/4] <= sum(uDZ0[r,j+i] for i in range(4)), name = 'uMC')
        self.model.addConstr(uDZ0[r,j+0] + uDZ0[r,j+1] + uDZ0[r,j+2] + uDZ0[r,j+3] + 
                             uDW0[r,j+0] + uDW0[r,j+1] + uDW0[r,j+2] + uDW0[r,j+3] >= 5*uAC0[r,j/4], name = 'uMC')
        self.model.addConstr(uDZ0[r,j+0] + uDZ0[r,j+1] + uDZ0[r,j+2] + uDZ0[r,j+3] + 
                             uDW0[r,j+0] + uDW0[r,j+1] + uDW0[r,j+2] + uDW0[r,j+3] <= 8*uAC0[r,j/4], name = 'uMC')
    '''
    -----------------------------------------------------------------------------
    TWO: Key schedule in Eb+Eu+Em (No direction)
    *PARA.: 
      uLANE0:         A LANE of permutation of key schedule (for Type1 cancellation)
      uT2Can(A/B/C): Auxiliary variables to count Type2 cancellation
      uType20:        Type2 cancellation
    -----------------------------------------------------------------------------
    '''
    uLANE0   = self.model.addVars(16, vtype = GRB.BINARY, name = 'uLANE0')
    uT2CanA0 = self.model.addVars(end_round_u, 4, lb = 0, vtype = GRB.INTEGER, name = 'uT2CanA0')
    uT2CanB0 = self.model.addVars(end_round_u, 4, lb = 0, vtype = GRB.INTEGER, name = 'uT2CanB0')
    uT2CanC0 = self.model.addVars(range(1, end_round_u), 4, lb = 0, vtype = GRB.INTEGER, name = 'uT2CanC0')
    uType10  = self.model.addVars(16, lb = 0, vtype = GRB.INTEGER, name = 'uType10')
    uType20  = self.model.addVars(end_round_u, 4, lb = 0, vtype = GRB.INTEGER, name = 'uType20')

    self.model.addConstr(sum(uLANE0[i] for i in range(16)) >= 1)

    ''' Mark LANE '''
    for i in range(16):
      for r in range(end_round_u + 1):
        self.model.addConstr(uLANE0[i] >= ustk0[r, self.hTable[i][r]], name = 'uMLANE')

    ''' Type 1 cancellation '''
    for i in range(16):
      self.model.addConstr(uType10[i] == (end_round_u + 1) * uLANE0[i] - 
                           sum(ustk0[r,self.hTable[i][r]] for r in range(end_round_u + 1)), name = 'uType10')
      self.model.addConstr(uType10[i] <= self.s - 1, name = 'uType10') # Cancellation no more than 2/1 in each position (for Deoxys-BC-384/256)
        
    ''' Type 2 cancellation (In Eu and Em) '''
    for r in range(end_round_u):
      for c in range(4):
        # uT2CanA0: Active bytes in one column, before MC
        self.model.addConstr(uT2CanA0[r,c] == uDZ0[r,4*c+0] + uDZ0[r,4*c+1] + uDZ0[r,4*c+2] + uDZ0[r,4*c+3], name = 'uType20A')
        # uT2CanB0: Inactive bytes in one column, after MC
        self.model.addConstr(uT2CanB0[r,c] == 4*uAC0[r,c] - uDW0[r,4*c+0] - uDW0[r,4*c+1] - uDW0[r,4*c+2] - uDW0[r,4*c+3], name = 'uType20B')                             

    # uT2CanC0: Cancellation appears in AK, in the next round
    for r in range(1, end_round_u):
      for c in range(4):
        self.model.addConstr(uT2CanC0[r,c] == ucan0[r,4*c+0] + ucan0[r,4*c+1] + ucan0[r,4*c+2] + ucan0[r,4*c+3], name = 'uType20C')

    # uType20 = -(uT2CanA0 - uT2CanB0 - uT2CanC0)
    for r in range(end_round_u-1):
      for c in range(4):
        self.model.addConstr(uType20[r,c] >= (uT2CanB0[r,c] + uT2CanC0[r+1,c] - uT2CanA0[r,c]), name = 'uType20')
        self.model.addConstr(uType20[r,c] >= 0, name = 'uType20')

    ''' Counting of all cancellations'''
    # NOTE: s*LANE[0~15] - 1 >= Type1Can * (LANE[0~15]-stk[0~r,0~15]) + Type2Can.
    self.model.addConstr(self.s * sum(uLANE0[i] for i in range(16)) -  # s * LANE
                         sum(uType10[i] for i in range(16)) - 
                         sum(uType20[r,i] for r in range(end_round_u) for i in range(4))
                         >= 1, name = 'uCAN')
    # --------------------------------------------
    # ---------- The second upper trail ----------
    # --------------------------------------------
    uDW1  = self.model.addVars(range(-1, end_round_u),          16, vtype = GRB.BINARY, name = 'uDW1')
    udw1  = self.model.addVars(range(-1, end_round_u),          16, vtype = GRB.BINARY, name = 'udw1')# DWw[-1,](plaintext)
    uDX1  = self.model.addVars(end_round_u, 16, vtype = GRB.BINARY, name = 'uDX1')
    udx1  = self.model.addVars(end_round_u, 16, vtype = GRB.BINARY, name = 'udx1')
    uDY1  = self.model.addVars(end_round_u,                     16, vtype = GRB.BINARY, name = 'uDY1')
    udy1  = self.model.addVars(end_round_u,                     16, vtype = GRB.BINARY, name = 'udy1')
    uDZ1  = self.model.addVars(end_round_u,                     16, vtype = GRB.BINARY, name = 'uDZ1')
    udz1  = self.model.addVars(end_round_u,                     16, vtype = GRB.BINARY, name = 'udz1')
    ustk1 = self.model.addVars(end_round_u + 1, 16, vtype = GRB.BINARY, name = 'ustk1')
    ucan1 = self.model.addVars(end_round_u,     16, vtype = GRB.BINARY, name = 'ucan1')
    uAC1  = self.model.addVars(end_round_u,                      4, vtype = GRB.BINARY, name = 'uAC1') # Mark active in one column
    uTC1  = self.model.addVars(end_round_u,                      4, vtype = GRB.BINARY, name = 'uTC1') # Mark truncated in one column
    
    '''Basic Constrs:'''
    # Remove dx[r,i] > DX[r,i]
    for i in range(16):
      for r in range(-1, end_round_u):
        self.model.addConstr(uDW1[r,i] >= udw1[r,i], name = 'uBas')
      for r in range(end_round_u):
        self.model.addConstr(uDX1[r,i] >= udx1[r,i], name = 'uBas')
        self.model.addConstr(uDY1[r,i] >= udy1[r,i], name = 'uBas')
        self.model.addConstr(uDZ1[r,i] >= udz1[r,i], name = 'uBas')
        
    '''Operation: ATK'''
    for r in range(-1, end_round_u - 1):
      for i in range(16):
        self.model.addConstr(uDW1[r,i] - udx1[r+1,i] - ucan1[r+1,i] >= 0, name = 'uATK')
        self.model.addConstr(- uDW1[r,i] + uDX1[r+1,i] + ucan1[r+1,i] >= 0, name = 'uATK')
        self.model.addConstr(uDW1[r,i] + ustk1[r+1,i] - uDX1[r+1,i] - 2 * ucan1[r+1,i] >= 0, name = 'uATK')
        self.model.addConstr(- ustk1[r+1,i] + uDX1[r+1,i] + ucan1[r+1,i] >= 0, name = 'uATK')
        self.model.addConstr(udw1[r,i] == udx1[r+1,i], name = 'uATK')

    '''Operation: SC'''
    # In Extension     [DX,dx] <= [DY,dy]: [1,1] <= [1,0]/[1,1], [0,0] <= [0,0]
    for r in range(0, self.round_Eb):
      for i in range(16):
        self.model.addConstr(uDX1[r,i] == uDY1[r,i], name = 'bSC')
        self.model.addConstr(udx1[r,i] == uDY1[r,i], name = 'bSC')
    # In distinguisher [DX,dx] => [DY,dy]: [1,0]/[1,1] => [1,1], [0,0] => [0,0]
    for r in range(self.round_Eb, self.round_Eb + self.round_Dis):
      for i in range(16):
        self.model.addConstr(uDY1[r,i] == uDX1[r,i], name = 'uSC')
        self.model.addConstr(udy1[r,i] == uDX1[r,i], name = 'uSC')

    '''Operation: SR''' 
    # [DY,dy] <=> (DZ,dz)
    for r in range(self.round_Eb + self.round_Dis):
      for i in range(16):
        self.model.addConstr(uDY1[r,self.SRpermutation_rev[i]] == uDZ1[r,i], name = 'uSR')
        self.model.addConstr(udy1[r,self.SRpermutation_rev[i]] == udz1[r,i], name = 'uSR')

    '''Operation: MC'''
    # In Extension: [DZ,dz] => [DW,dw]
    for r in range(self.round_Eb):
      for j in [0,4,8,12]:
        self.model.addConstr(uTC1[r,j/4] >= udw1[r,j+0], name = 'bMC')
        self.model.addConstr(uTC1[r,j/4] >= udw1[r,j+1], name = 'bMC')
        self.model.addConstr(uTC1[r,j/4] >= udw1[r,j+2], name = 'bMC')
        self.model.addConstr(uTC1[r,j/4] >= udw1[r,j+3], name = 'bMC')
        self.model.addConstr(uTC1[r,j/4] <= sum(udw1[r,j+i] for i in range(4)), name = 'bMC')
        self.model.addConstr(uTC1[r,j/4] <= udz1[r,j+0], name = 'bMC')
        self.model.addConstr(uTC1[r,j/4] <= udz1[r,j+1], name = 'bMC')
        self.model.addConstr(uTC1[r,j/4] <= udz1[r,j+2], name = 'bMC')
        self.model.addConstr(uTC1[r,j/4] <= udz1[r,j+3], name = 'bMC')
        self.model.addConstr(4*uTC1[r,j/4] >= sum(udz1[r,j+i] for i in range(4)), name = 'bMC')

        self.model.addConstr(uAC1[r,j/4] >= uDW1[r,j+0], name = 'bMC')
        self.model.addConstr(uAC1[r,j/4] >= uDW1[r,j+1], name = 'bMC')
        self.model.addConstr(uAC1[r,j/4] >= uDW1[r,j+2], name = 'bMC')
        self.model.addConstr(uAC1[r,j/4] >= uDW1[r,j+3], name = 'bMC')
        self.model.addConstr(uAC1[r,j/4] <= sum(uDW1[r,j+i] for i in range(4)), name = 'bMC')
        self.model.addConstr(uDW1[r,j+0] + uDW1[r,j+1] + uDW1[r,j+2] + uDW1[r,j+3] + 
                             uDZ1[r,j+0] + uDZ1[r,j+1] + uDZ1[r,j+2] + uDZ1[r,j+3] >= 5*uAC1[r,j/4], name = 'bMC')
        self.model.addConstr(uDW1[r,j+0] + uDW1[r,j+1] + uDW1[r,j+2] + uDW1[r,j+3] + 
                             uDZ1[r,j+0] + uDZ1[r,j+1] + uDZ1[r,j+2] + uDZ1[r,j+3] <= 8*uAC1[r,j/4], name = 'bMC')
    # In Distinguisher [DW,dw] <= [DZ,dz]
    for r in range(self.round_Eb, self.round_Eb + self.round_Dis):
      for j in [0,4,8,12]:
        self.model.addConstr(uTC1[r,j/4] >= udz1[r,j+0], name = 'uMC')
        self.model.addConstr(uTC1[r,j/4] >= udz1[r,j+1], name = 'uMC')
        self.model.addConstr(uTC1[r,j/4] >= udz1[r,j+2], name = 'uMC')
        self.model.addConstr(uTC1[r,j/4] >= udz1[r,j+3], name = 'uMC')
        self.model.addConstr(uTC1[r,j/4] <= sum(udz1[r,j+i] for i in range(4)), name = 'uMC')
        self.model.addConstr(uTC1[r,j/4] <= udw1[r,j+0], name = 'uMC')
        self.model.addConstr(uTC1[r,j/4] <= udw1[r,j+1], name = 'uMC')
        self.model.addConstr(uTC1[r,j/4] <= udw1[r,j+2], name = 'uMC')
        self.model.addConstr(uTC1[r,j/4] <= udw1[r,j+3], name = 'uMC')
        self.model.addConstr(4*uTC1[r,j/4] >= sum(udw1[r,j+i] for i in range(4)), name = 'uMC')

        self.model.addConstr(uAC1[r,j/4] >= uDZ1[r,j+0], name = 'uMC')
        self.model.addConstr(uAC1[r,j/4] >= uDZ1[r,j+1], name = 'uMC')
        self.model.addConstr(uAC1[r,j/4] >= uDZ1[r,j+2], name = 'uMC')
        self.model.addConstr(uAC1[r,j/4] >= uDZ1[r,j+3], name = 'uMC')
        self.model.addConstr(uAC1[r,j/4] <= sum(uDZ1[r,j+i] for i in range(4)), name = 'uMC')
        self.model.addConstr(uDZ1[r,j+0] + uDZ1[r,j+1] + uDZ1[r,j+2] + uDZ1[r,j+3] + 
                             uDW1[r,j+0] + uDW1[r,j+1] + uDW1[r,j+2] + uDW1[r,j+3] >= 5*uAC1[r,j/4], name = 'uMC')
        self.model.addConstr(uDZ1[r,j+0] + uDZ1[r,j+1] + uDZ1[r,j+2] + uDZ1[r,j+3] + 
                             uDW1[r,j+0] + uDW1[r,j+1] + uDW1[r,j+2] + uDW1[r,j+3] <= 8*uAC1[r,j/4], name = 'uMC')
    '''
    -----------------------------------------------------------------------------
    TWO: Key schedule in Eb+Eu+Em (No direction)
    *PARA.: 
      uLANE0:         A LANE of permutation of key schedule (for Type1 cancellation)
      uT2Can(A/B/C): Auxiliary variables to count Type2 cancellation
      uType20:        Type2 cancellation
    -----------------------------------------------------------------------------
    '''
    uLANE1   = self.model.addVars(16, vtype = GRB.BINARY, name = 'uLANE1')
    uT2CanA1 = self.model.addVars(end_round_u, 4, lb = 0, vtype = GRB.INTEGER, name = 'uT2CanA1')
    uT2CanB1 = self.model.addVars(end_round_u, 4, lb = 0, vtype = GRB.INTEGER, name = 'uT2CanB1')
    uT2CanC1 = self.model.addVars(range(1, end_round_u), 4, lb = 0, vtype = GRB.INTEGER, name = 'uT2CanC1')
    uType11  = self.model.addVars(16, lb = 0, vtype = GRB.INTEGER, name = 'uType11')
    uType21  = self.model.addVars(end_round_u, 4, lb = 0, vtype = GRB.INTEGER, name = 'uType21')

    self.model.addConstr(sum(uLANE1[i] for i in range(16)) >= 1)

    ''' Mark LANE '''
    for i in range(16):
      for r in range(end_round_u + 1):
        self.model.addConstr(uLANE1[i] >= ustk1[r, self.hTable[i][r]], name = 'uMLANE')

    ''' Type 1 cancellation '''
    for i in range(16):
      self.model.addConstr(uType11[i] == (end_round_u + 1) * uLANE1[i] - 
                           sum(ustk1[r,self.hTable[i][r]] for r in range(end_round_u + 1)), name = 'uType11')
      self.model.addConstr(uType11[i] <= self.s - 1, name = 'uType11') # Cancellation no more than 2/1 in each position (for Deoxys-BC-384/256)
        
    ''' Type 2 cancellation (In Eu and Em) '''
    for r in range(end_round_u):
      for c in range(4):
        # uT2CanA1: Active bytes in one column, before MC
        self.model.addConstr(uT2CanA1[r,c] == uDZ1[r,4*c+0] + uDZ1[r,4*c+1] + uDZ1[r,4*c+2] + uDZ1[r,4*c+3], name = 'uType21A')
        # uT2CanB1: Inactive bytes in one column, after MC
        self.model.addConstr(uT2CanB1[r,c] == 4*uAC1[r,c] - uDW1[r,4*c+0] - uDW1[r,4*c+1] - uDW1[r,4*c+2] - uDW1[r,4*c+3], name = 'uType21B')                             

    # uT2CanC1: Cancellation appears in AK, in the next round
    for r in range(1, end_round_u):
      for c in range(4):
        self.model.addConstr(uT2CanC1[r,c] == ucan1[r,4*c+0] + ucan1[r,4*c+1] + ucan1[r,4*c+2] + ucan1[r,4*c+3], name = 'uType21C')

    # uType21 = -(uT2CanA1 - uT2CanB1 - uT2CanC1)
    for r in range(end_round_u-1):
      for c in range(4):
        self.model.addConstr(uType21[r,c] >= (uT2CanB1[r,c] + uT2CanC1[r+1,c] - uT2CanA1[r,c]), name = 'uType21')
        self.model.addConstr(uType21[r,c] >= 0, name = 'uType21')

    ''' Counting of all cancellations'''
    # NOTE: s*LANE[0~15] - 1 >= Type1Can * (LANE[0~15]-stk[0~r,0~15]) + Type2Can.
    self.model.addConstr(self.s * sum(uLANE1[i] for i in range(16)) -  # s * LANE
                         sum(uType11[i] for i in range(16)) - 
                         sum(uType21[r,i] for r in range(end_round_u) for i in range(4))
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
    # --------------------------------------------
    # ------- The first Guess-&-Determine --------
    # --------------------------------------------
    uGstk0 = self.model.addVars(self.round_Eb,            16, vtype = GRB.BINARY, name = 'uGstk0')
    uDetX0 = self.model.addVars(self.round_Eb,            16, vtype = GRB.BINARY, name = 'uDetX0')
    uFrSB0 = self.model.addVars(self.round_Eb,            16, vtype = GRB.BINARY, name = 'uFrSB0')
    uDetY0 = self.model.addVars(self.round_Eb,            16, vtype = GRB.BINARY, name = 'uDetY0')
    uDetZ0 = self.model.addVars(self.round_Eb,            16, vtype = GRB.BINARY, name = 'uDetZ0')
    uDetW0 = self.model.addVars(range(-1, self.round_Eb), 16, vtype = GRB.BINARY, name = 'uDetW0')
    uDetC0 = self.model.addVars(self.round_Eb,             4, vtype = GRB.BINARY, name = 'uDetC0') # mark 'det.' in one column
    uDiffDetMC0 = self.model.addVars(self.round_Eb,       16, vtype = GRB.BINARY, name = 'uDiffDetMC0') # mark diff. det. of cell
    uDiffDetC0 = self.model.addVars(self.round_Eb,         4, vtype = GRB.BINARY, name = 'uDiffDetC0') # mark diff. det. in one colume
    uFrMC0 = self.model.addVars(self.round_Eb,            16, vtype = GRB.BINARY, name = 'uFrMC0') 
    
    ''' KR: Guess round key and determine '''
    # uDetW0 -(uGstk0)-> uDetX0
    for r in range(-1, self.round_Eb-1):
      for i in range(16): 
        self.model.addConstr(uGstk0[r+1,i] >= uDetX0[r+1,i], name = 'bDet_ATK')
        self.model.addConstr(uDetW0[r,i] >= uDetX0[r+1,i], name = 'bDet_ATK')
        self.model.addConstr(- uDetW0[r,i] - uGstk0[r+1,i] + uDetX0[r+1,i] >= -1, name = 'bDet_ATK')

    ''' KR: Determine and Filter in SC '''
    for r in range(self.round_Eb):
      for i in range(16):
        # 'Det.' propagation via SC
        self.model.addConstr(uDetY0[r,i] == uDetX0[r,i], name = 'bDet_SC') 
        # Filter obtain from SC
        self.model.addConstr(uDY0[r,i] - udy0[r,i] - uFrSB0[r,i] >= 0, name = 'bFr_SC')
        self.model.addConstr(uDetY0[r,i] >= uFrSB0[r,i], name = 'bFr_SC')
        self.model.addConstr(- uDetY0[r,i] - uDY0[r,i] + udy0[r,i] + uFrSB0[r,i] >= -1, name = 'bFr_SC')
    
    ''' KR: 'Det.' propagation in SR '''
    for r in range(self.round_Eb):
      for i in range(16):
        self.model.addConstr(uDetZ0[r,i] == uDetY0[r,self.SRpermutation_rev[i]], name = 'bDet_SR')

    ''' KR: 'Det.' propagation in MC* '''
    for r in range(self.round_Eb):
      for c in range(4):
        # 'Det.' propagation via MC
        # NOTE: uDetC0 = 1 when all uDetZ0 = 1, else uDetC0 = 0 
        self.model.addConstr(uDetC0[r,c] <= uDetZ0[r,4*c+0], name = 'bDet_MC')
        self.model.addConstr(uDetC0[r,c] <= uDetZ0[r,4*c+1], name = 'bDet_MC')
        self.model.addConstr(uDetC0[r,c] <= uDetZ0[r,4*c+2], name = 'bDet_MC')
        self.model.addConstr(uDetC0[r,c] <= uDetZ0[r,4*c+3], name = 'bDet_MC')
        self.model.addConstr(uDetZ0[r,4*c+0] + uDetZ0[r,4*c+1] + uDetZ0[r,4*c+2] + uDetZ0[r,4*c+3] - uDetC0[r,c] <= 3, name = 'bDet_MC')
        # NOTE: all uDetW0 = 1 when uDetC0 = 1, else all uDetW0 = 0
        self.model.addConstr(uDetC0[r,c] == uDetW0[r,4*c+0], name = 'bDet_MC')
        self.model.addConstr(uDetC0[r,c] == uDetW0[r,4*c+1], name = 'bDet_MC')
        self.model.addConstr(uDetC0[r,c] == uDetW0[r,4*c+2], name = 'bDet_MC')
        self.model.addConstr(uDetC0[r,c] == uDetW0[r,4*c+3], name = 'bDet_MC')

    ''' KR: Determine whether the difference before MC is determinable '''
    for r in range(self.round_Eb):
      for i in range(16):
        # NOTE: uDiffDetMC0 = uDetZ0 when uDZ = 1 (ACTIVE), else uDiffDetMC0 = 1
        self.model.addConstr(uDiffDetMC0[r,i] >= uDetZ0[r,i], name = 'bDet_MC')
        self.model.addConstr(uDetZ0[r,i] - uDZ0[r,i] - uDiffDetMC0[r,i] >= -1, name = 'bDet_MC')
        self.model.addConstr(uDZ0[r,i] + uDiffDetMC0[r,i] >= 1, name = 'bDet_MC')
    for r in range(self.round_Eb):
      for c in range(4):
        # NOTE: uDiffDetC0 = 1 when all uDiffDetMC0 = 1, else uDiffDetC0 = 0 
        self.model.addConstr(uDiffDetC0[r,c] <= uDiffDetMC0[r,4*c+0], name = 'bDet_MC')
        self.model.addConstr(uDiffDetC0[r,c] <= uDiffDetMC0[r,4*c+1], name = 'bDet_MC')
        self.model.addConstr(uDiffDetC0[r,c] <= uDiffDetMC0[r,4*c+2], name = 'bDet_MC')
        self.model.addConstr(uDiffDetC0[r,c] <= uDiffDetMC0[r,4*c+3], name = 'bDet_MC')
        self.model.addConstr(uDiffDetMC0[r,4*c+0] + uDiffDetMC0[r,4*c+1] + uDiffDetMC0[r,4*c+2] + uDiffDetMC0[r,4*c+3] 
                             - uDiffDetC0[r,c] <= 3, name = 'bDet_MC')
    
    ''' KR: Filter obtain from MC (for cell -> column) '''
    for r in range(self.round_Eb):
      for c in range(4):
        for i in range(4):
          # NOTE: uFrMC0 = 0 when uDiffDetC0 = 0, 
          # NOTE: else uFrMC0 = 1 if and only if {uTC = 1 (Truncated before MC) and udw = 0 (Fixed after MC)}
          self.model.addConstr(uTC0[r,c] >= uFrMC0[r,4*c+i], name = 'bFr_MC')
          self.model.addConstr(udw0[r,4*c+i] + uFrMC0[r,4*c+i] <= 1, name = 'bFr_MC')
          self.model.addConstr(uDiffDetC0[r,c] + uTC0[r,c] - udw0[r,4*c+i] - uFrMC0[r,4*c+i] <= 1, name = 'bFr_MC')
          self.model.addConstr(uDiffDetC0[r,c] >= uFrMC0[r,4*c+i], name = 'bFr_MC')

    # --------------------------------------------
    # ------- The second Guess-&-Determine -------
    # --------------------------------------------
    uGstk1 = self.model.addVars(self.round_Eb,            16, vtype = GRB.BINARY, name = 'uGstk1')
    uDetX1 = self.model.addVars(self.round_Eb,            16, vtype = GRB.BINARY, name = 'uDetX1')
    uFrSB1 = self.model.addVars(self.round_Eb,            16, vtype = GRB.BINARY, name = 'uFrSB1')
    uDetY1 = self.model.addVars(self.round_Eb,            16, vtype = GRB.BINARY, name = 'uDetY1')
    uDetZ1 = self.model.addVars(self.round_Eb,            16, vtype = GRB.BINARY, name = 'uDetZ1')
    uDetW1 = self.model.addVars(range(-1, self.round_Eb), 16, vtype = GRB.BINARY, name = 'uDetW1')
    uDetC1 = self.model.addVars(self.round_Eb,             4, vtype = GRB.BINARY, name = 'uDetC1') # mark 'det.' in one column
    uDiffDetMC1 = self.model.addVars(self.round_Eb,       16, vtype = GRB.BINARY, name = 'uDiffDetMC1') # mark diff. det. of cell
    uDiffDetC1 = self.model.addVars(self.round_Eb,         4, vtype = GRB.BINARY, name = 'uDiffDetC1') # mark diff. det. in one colume
    uFrMC1 = self.model.addVars(self.round_Eb,            16, vtype = GRB.BINARY, name = 'uFrMC1')
    
    ''' KR: Guess round key and determine '''
    # uDetW1 -(uGstk1)-> uDetX1
    for r in range(-1, self.round_Eb-1):
      for i in range(16): 
        self.model.addConstr(uGstk1[r+1,i] >= uDetX1[r+1,i], name = 'bDet_ATK')
        self.model.addConstr(uDetW1[r,i] >= uDetX1[r+1,i], name = 'bDet_ATK')
        self.model.addConstr(- uDetW1[r,i] - uGstk1[r+1,i] + uDetX1[r+1,i] >= -1, name = 'bDet_ATK')

    ''' KR: Determine and Filter in SC '''
    for r in range(self.round_Eb):
      for i in range(16):
        # 'Det.' propagation via SC
        self.model.addConstr(uDetY1[r,i] == uDetX1[r,i], name = 'bDet_SC') 
        # Filter obtain from SC
        self.model.addConstr(uDY1[r,i] - udy1[r,i] - uFrSB1[r,i] >= 0, name = 'bFr_SC')
        self.model.addConstr(uDetY1[r,i] >= uFrSB1[r,i], name = 'bFr_SC')
        self.model.addConstr(- uDetY1[r,i] - uDY1[r,i] + udy1[r,i] + uFrSB1[r,i] >= -1, name = 'bFr_SC')
    
    ''' KR: 'Det.' propagation in SR '''
    for r in range(self.round_Eb):
      for i in range(16):
        self.model.addConstr(uDetZ1[r,i] == uDetY1[r,self.SRpermutation_rev[i]], name = 'bDet_SR')

    ''' KR: 'Det.' propagation in MC* '''
    for r in range(self.round_Eb):
      for c in range(4):
        # 'Det.' propagation via MC
        # NOTE: uDetC1 = 1 when all uDetZ1 = 1, else uDetC1 = 0 
        self.model.addConstr(uDetC1[r,c] <= uDetZ1[r,4*c+0], name = 'bDet_MC')
        self.model.addConstr(uDetC1[r,c] <= uDetZ1[r,4*c+1], name = 'bDet_MC')
        self.model.addConstr(uDetC1[r,c] <= uDetZ1[r,4*c+2], name = 'bDet_MC')
        self.model.addConstr(uDetC1[r,c] <= uDetZ1[r,4*c+3], name = 'bDet_MC')
        self.model.addConstr(uDetZ1[r,4*c+0] + uDetZ1[r,4*c+1] + uDetZ1[r,4*c+2] + uDetZ1[r,4*c+3] - uDetC1[r,c] <= 3, name = 'bDet_MC')
        # NOTE: all uDetW1 = 1 when uDetC1 = 1, else all uDetW1 = 0
        self.model.addConstr(uDetC1[r,c] == uDetW1[r,4*c+0], name = 'bDet_MC')
        self.model.addConstr(uDetC1[r,c] == uDetW1[r,4*c+1], name = 'bDet_MC')
        self.model.addConstr(uDetC1[r,c] == uDetW1[r,4*c+2], name = 'bDet_MC')
        self.model.addConstr(uDetC1[r,c] == uDetW1[r,4*c+3], name = 'bDet_MC')

    ''' KR: Determine whether the difference before MC is determinable '''
    for r in range(self.round_Eb):
      for i in range(16):
        # NOTE: uDiffDetMC1 = uDetZ1 when uDZ = 1 (ACTIVE), else uDiffDetMC1 = 1
        self.model.addConstr(uDiffDetMC1[r,i] >= uDetZ1[r,i], name = 'bDet_MC')
        self.model.addConstr(uDetZ1[r,i] - uDZ1[r,i] - uDiffDetMC1[r,i] >= -1, name = 'bDet_MC')
        self.model.addConstr(uDZ1[r,i] + uDiffDetMC1[r,i] >= 1, name = 'bDet_MC')
    for r in range(self.round_Eb):
      for c in range(4):
        # NOTE: uDiffDetC1 = 1 when all uDiffDetMC1 = 1, else uDiffDetC1 = 0 
        self.model.addConstr(uDiffDetC1[r,c] <= uDiffDetMC1[r,4*c+0], name = 'bDet_MC')
        self.model.addConstr(uDiffDetC1[r,c] <= uDiffDetMC1[r,4*c+1], name = 'bDet_MC')
        self.model.addConstr(uDiffDetC1[r,c] <= uDiffDetMC1[r,4*c+2], name = 'bDet_MC')
        self.model.addConstr(uDiffDetC1[r,c] <= uDiffDetMC1[r,4*c+3], name = 'bDet_MC')
        self.model.addConstr(uDiffDetMC1[r,4*c+0] + uDiffDetMC1[r,4*c+1] + uDiffDetMC1[r,4*c+2] + uDiffDetMC1[r,4*c+3] 
                             - uDiffDetC1[r,c] <= 3, name = 'bDet_MC')
    
    ''' KR: Filter obtain from MC (for cell -> column) '''
    for r in range(self.round_Eb):
      for c in range(4):
        for i in range(4):
          # NOTE: uFrMC1 = 0 when uDiffDetC1 = 0, 
          # NOTE: else uFrMC1 = 1 if and only if {uTC = 1 (Truncated before MC) and udw = 0 (Fixed after MC)}
          self.model.addConstr(uTC1[r,c] >= uFrMC1[r,4*c+i], name = 'bFr_MC')
          self.model.addConstr(udw1[r,4*c+i] + uFrMC1[r,4*c+i] <= 1, name = 'bFr_MC')
          self.model.addConstr(uDiffDetC1[r,c] + uTC1[r,c] - udw1[r,4*c+i] - uFrMC1[r,4*c+i] <= 1, name = 'bFr_MC')
          self.model.addConstr(uDiffDetC1[r,c] >= uFrMC1[r,4*c+i], name = 'bFr_MC')
        

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
      The more generic contradiction detection
      * iBCT
      * XOR of 4 diff. != 0
    -----------------------------------------------------------------------------------------------
    ''' 
    duFX0 = self.model.addVars(end_round_u, 16, vtype = GRB.BINARY, name = 'duFX0')
    duFX1 = self.model.addVars(end_round_u, 16, vtype = GRB.BINARY, name = 'duFX1')
    duTX0 = self.model.addVars(end_round_u, 16, vtype = GRB.BINARY, name = 'duTX0')
    duTX1 = self.model.addVars(end_round_u, 16, vtype = GRB.BINARY, name = 'duTX1')
    duAZ0 = self.model.addVars(end_round_u, 16, vtype = GRB.BINARY, name = 'duAZ0')
    duAZ1 = self.model.addVars(end_round_u, 16, vtype = GRB.BINARY, name = 'duAZ1')
    duAW0 = self.model.addVars(end_round_u, 16, vtype = GRB.BINARY, name = 'duAW0')
    duAW1 = self.model.addVars(end_round_u, 16, vtype = GRB.BINARY, name = 'duAW1')
    dlAX = self.model.addVars(end_round_u, 16, vtype = GRB.BINARY, name = 'dlAX')
    dlFZ = self.model.addVars(end_round_u, 16, vtype = GRB.BINARY, name = 'dlFZ')
    dlFW = self.model.addVars(end_round_u, 16, vtype = GRB.BINARY, name = 'dlFW')
    # power-reduced aid
    PRA0 = self.model.addVars(end_round_u, 16, vtype = GRB.BINARY, name = 'PRA0')
    PRA1 = self.model.addVars(end_round_u, 16, vtype = GRB.BINARY, name = 'PRA1')
    PRA2 = self.model.addVars(end_round_u, 16, vtype = GRB.BINARY, name = 'PRA2')
    PRA3 = self.model.addVars(end_round_u, 16, vtype = GRB.BINARY, name = 'PRA3')
    PRA4 = self.model.addVars(end_round_u, 16, vtype = GRB.BINARY, name = 'PRA4')

    for r in range(self.round_Eb, self.round_Eb + self.round_Dis):
      for i in range(16):
        # upper
        self.model.addConstr(duFX0[r,i] == uDX0[r,i] - udx0[r,i], name = '2duFX0')
        self.model.addConstr(duFX1[r,i] == uDX1[r,i] - udx1[r,i], name = '2duFX1')
        self.model.addConstr(duTX0[r,i] == 1 - udx0[r,i], name = '2duTX0')
        self.model.addConstr(duTX1[r,i] == 1 - udx1[r,i], name = '2duTX1')
        self.model.addConstr(duAZ0[r,i] == 1 - uDZ0[r,i], name = '2duAZ0')
        self.model.addConstr(duAZ1[r,i] == 1 - uDZ1[r,i], name = '2duAZ1')
        self.model.addConstr(duAW0[r,i] == 1 - uDW0[r,i], name = '2duAW0')
        self.model.addConstr(duAW1[r,i] == 1 - uDW1[r,i], name = '2duAW1')
        # lower
        self.model.addConstr(dlAX[r,i] == 1 - lDX[r,i], name = '2lAX')
        self.model.addConstr(dlFZ[r,i] == lDZ[r,i] - ldz[r,i], name = '2lFZ')
        self.model.addConstr(dlFW[r,i] == lDW[r,i] - ldw[r,i], name = '2lFW')

        # Power-reduce
        # constraint-1 [(𝒖𝑫𝑿_𝟎−𝒖𝒅𝒙_𝟎 )+(𝒖𝑫𝑿_𝟏−𝒖𝒅𝒙_𝟏 )]×(𝟏−𝒖𝒅𝒙_𝟎 )×(𝟏−𝒖𝒅𝒙_𝟏 )×(𝟏−𝒍𝑫𝑿)
        self.model.addGenConstrAnd(PRA0[r, i], [duTX0[r,i], duTX1[r,i]],name = '2PDA0')
        self.model.addGenConstrAnd(PRA1[r, i], [PRA0[r,i], dlAX[r,i]],name = '2PDA1')

        # constraint-2 (𝒍𝑫𝒁−𝒍𝒅𝒛)×(𝟏−𝒖𝑫𝒁_𝟎 )×(𝟏−𝒖𝑫𝒁_𝟏)
        self.model.addGenConstrAnd(PRA2[r, i], [duAZ0[r,i], duAZ1[r,i]],name = '2PDA2')

        # constraint-3 (𝒍𝑫𝑾−𝒍𝒅𝒘)×(𝟏−𝒖𝑫𝑾_𝟎 )×(𝟏−𝒖𝑫𝑾_𝟏)
        self.model.addGenConstrAnd(PRA3[r, i], [duAW0[r,i], duAW1[r,i]],name = '2PDA3')

        # constraint-4 (𝒖𝑫𝑿_𝟎−𝒖𝒅𝒙_𝟎 )∗(𝒖𝑫𝑿_𝟏−𝒖𝒅𝒙_𝟏 )∗(𝒍𝑫𝒁−𝒍𝒅𝒛)
        self.model.addGenConstrAnd(PRA4[r, i], [duFX0[r,i], duFX1[r,i]],name = '2PDA4')

    self.model.addConstr(sum(
    (
      ((duFX0[r,i] + duFX1[r,i]) * PRA1[r,i])
      + (dlFZ[r,i] * PRA2[r,i])
      + (dlFW[r,i] * PRA3[r,i])
      + (dlFZ[r,i] * PRA4[r,i])
    )
    for i in range(16) for r in range(self.round_Eb, self.round_Eb + self.round_Dis)) >= 1, name = 'MainContra')










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
      self.model.addConstr(uDetW0[-1,i] == 1)
      self.model.addConstr(uDetW1[-1,i] == 1)
      self.model.addConstr(ldetEQZ[end_round_l - 1,i] == 1)

    '''
    Some parameters needed to cluster 2 diff. trails.
    '''
    pxt  = self.model.addVars(16, vtype = GRB.BINARY, name = 'pxt')# DWw[-1,](plaintext)
    udx  = self.model.addVars(self.round_Eb + self.round_Dis, 16, vtype = GRB.BINARY, name = 'udx') # count the union of udx0 and udx1 in Eb
    uGstk = self.model.addVars(self.round_Eb, 16, vtype = GRB.BINARY, name = 'uGstk')
    for r in range(self.round_Eb + self.round_Dis):
      for i in range(16):
        self.model.addConstr(udx[r,i] >= udx0[r,i])
        self.model.addConstr(udx[r,i] >= udx1[r,i])
        self.model.addConstr(udx[r,i] <= udx0[r,i] + udx1[r,i])
    for r in range(self.round_Eb):
        for i in range(16):
          self.model.addConstr(uGstk[r,i] >= uGstk0[r,i])
          self.model.addConstr(uGstk[r,i] >= uGstk1[r,i])
          self.model.addConstr(uGstk[r,i] <= uGstk0[r,i] + uGstk1[r,i])
    for i in range(16):
      self.model.addConstr(pxt[i] >= udw0[-1,i])
      self.model.addConstr(pxt[i] >= udw1[-1,i])
      self.model.addConstr(pxt[i] <= udw0[-1,i] + udw1[-1,i])

    ''' Para. for Complexity'''

    rb0  = self.model.addVar(vtype = GRB.INTEGER, name = 'rb0')
    rb1  = self.model.addVar(vtype = GRB.INTEGER, name = 'rb1')
    rb  = self.model.addVar(vtype = GRB.INTEGER, name = 'rb')
    rf  = self.model.addVar(vtype = GRB.INTEGER, name = 'rf')
    cb0  = self.model.addVar(vtype = GRB.INTEGER, name = 'cb0')
    cb1  = self.model.addVar(vtype = GRB.INTEGER, name = 'cb1')
    cb  = self.model.addVar(vtype = GRB.INTEGER, name = 'cb')
    cf  = self.model.addVar(vtype = GRB.INTEGER, name = 'cf')
    mb0  = self.model.addVar(vtype = GRB.INTEGER, name = 'mb0')
    mb1  = self.model.addVar(vtype = GRB.INTEGER, name = 'mb1')
    mb  = self.model.addVar(vtype = GRB.INTEGER, name = 'mb')
    mf  = self.model.addVar(vtype = GRB.INTEGER, name = 'mf')
    mbp = self.model.addVar(vtype = GRB.INTEGER, name = 'mb\'')
    mfp = self.model.addVar(vtype = GRB.INTEGER, name = 'mf\'')
    cbp0 = self.model.addVar(vtype = GRB.INTEGER, name = 'cb0\'')
    cbp1 = self.model.addVar(vtype = GRB.INTEGER, name = 'cb1\'')
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
    self.model.addConstr(rb0  == self.cell_size * sum(udw0[-1,i] for i in range(16)))
    self.model.addConstr(rb1  == self.cell_size * sum(udw1[-1,i] for i in range(16)))
    self.model.addConstr(rb  == self.cell_size * sum(pxt[i] for i in range(16)))

    self.model.addConstr(cb0  == rb0 - self.cell_size * sum(udx0[self.round_Eb,i] for i in range(16)))
    self.model.addConstr(cb1  == rb1 - self.cell_size * sum(udx1[self.round_Eb,i] for i in range(16)))
    self.model.addConstr(cb  == rb - self.cell_size * sum(udx[self.round_Eb,i] for i in range(16)))

    self.model.addConstr(mb0  == self.cell_size * sum(udx0[r,i] for r in range(self.round_Eb) for i in range(16)))
    self.model.addConstr(mb1  == self.cell_size * sum(udx1[r,i] for r in range(self.round_Eb) for i in range(16)))
    self.model.addConstr(mb  == self.cell_size * sum(udx[r,i] for r in range(self.round_Eb) for i in range(16)))
    
    self.model.addConstr(mbp == self.cell_size * sum(uGstk[r,i] for r in range(self.round_Eb) for i in range(16)))

    self.model.addConstr(cbp0 == self.cell_size * sum(uFrSB0[r,i] + uFrMC0[r,i] for r in range(self.round_Eb) for i in range(16)))
    self.model.addConstr(cbp1 == self.cell_size * sum(uFrSB1[r,i] + uFrMC1[r,i] for r in range(self.round_Eb) for i in range(16)))

    ''' Defination of PARA. in Ef '''
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
      self.model.addConstr(T2 >= (mbp + mfp) + Dc + rb0 - cbp0)
      self.model.addConstr(T2 >= (mbp + mfp) + Dc + rb1 - cbp1)
    else:
      self.model.addConstr(T2 == (mbp + mfp) + Dc + (Dc + rf - cfp - self.block_size) + 1)
    # T31 (T31 = 2^{mb'+mf'-2cb'-2cf'} · Q)
    self.model.addConstr(T31 == (mbp + mfp - (cbp0 + cbp1) - 2*cfp) + Qc)
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

    # print('>>>>> Solution <<<<<')
    # print('='*90)
    # print('|| Choosen Plaintext | ' if self.pgP else('Choosen Chiphertext | '), end='')
    # print('PC = {:5} ||'.format(sum(pxt[-1,i].x + leqdz[end_round_l-1, i].x for i in range(16))))
    # print('-'*90)
    # print('|| x  = {:5} | z  = {:5} ||'.format(self.tx, self.tz))
    # print('-'*90)
    # print("|| rb = {:5} | cb = {:5} | mb = {:5} | mb\' = {:5} | cb\' = {:5} ||".format(rb.x, cb.x, mb.x, mbp.x, cbp.x))
    # print("|| rf = {:5} | cf = {:5} | mf = {:5} | mf\' = {:5} | cf\' = {:5} ||".format(rf.x, cf.x, mf.x, mfp.x, cfp.x))
    # print('-'*90)
    # print('|| D = {:5} | Dc  = {:5} | Q  = {:5} | Qc  = {:5} ||'.format(Dc.x, Dc.x+2, Qc.x, Qc.x-(cbp0.x + cbp1.x)-2*cfp.x))
    # print('-'*90)
    # print('|| T0 = {:5} | T1 = {:5} | T2 = {:5} | T31 = {:5} | T32 = {:5} | T4 = {:5} ||'.format(T0.x, T1.x, T2.x, T31.x, T32.x, T4.x))
    # print('-'*90)
    # print('|| Tc = {:5} ||'.format(Tc.x))
    # print('='*90)


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

    # ''' 11-r Deoxys-BC-256 '''
    # DeoxysBC256_11r = IB_DandJ(256, 2, 7, 2, 16, False) # adaptive adjustment to x=20 in this attack
    # DeoxysBC256_11r.ib_model()

    # ''' 13-r Deoxys-BC-384 '''
    # DeoxysBC384_13r = IB_DandJ(384, 2, 9, 2, 160, True)
    # DeoxysBC384_13r.ib_model()

    # ''' 14-r Deoxys-BC-384 '''
    # DeoxysBC384 = IB_DandJ(384, 2, 9, 3, 48, True)
    # DeoxysBC384.ib_model()

    # # ''' ----------------------------------------------------- '''

    # ''' 10-r Joltik-BC-128 '''
    # Joltik128_10r = IB_DandJ(128, 1, 7, 2, 52, True)
    # Joltik128_10r.ib_model()

    # ''' 11-r Joltik-BC-128 '''
    # Joltik128_11r = IB_DandJ(128, 2, 7, 2, 8, False)
    # Joltik128_11r.ib_model()

    # ''' 13-r Joltik-BC-192 '''
    # Joltik192_13r = IB_DandJ(192, 2, 9, 2, 86, False)
    # Joltik192_13r.ib_model()

    # ''' 14-r Joltik-BC-192 '''
    # Joltik192_14r = IB_DandJ(192, 2, 9, 3, 24 , True)
    # Joltik192_14r.ib_model()
