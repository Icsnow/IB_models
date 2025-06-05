import gurobipy as gp
from gurobipy import GRB
import math

class IB_ForkSKINNY:

  def __init__(self, b_size, Vs, rEb, rDis, rEf, Vx, cP) -> None:
    '''
    Give the round distribution of our attack:
                       -(Eu)-----'-----(El)-----'---(Ef)---  
                   :***************************************
    ---(Eb)---'----|                                   [r1]
    ***************+***********************************************************
    [ri]                                                                   [r0]
    NOTE: Set: {round|rEb <= ri; rEb+rEu >= ri; El+Ef span in r1}
    '''
    # assert(rEb <= ri), 'ERROR: rEb > ri'
    # assert(rEb+rEu >= ri), 'ERROR: rEb+rEu < ri'
    # assert(rEl+rEf <= r1), 'rEl+rEf > r1'

    self.Vs = Vs
    self.b_size = b_size
    self.k_size = b_size * Vs
    self.c_size = b_size // 16
    # self.ri = ri
    # self.r0 = r0
    # self.r1 = r1
    self.rEb = rEb
    self.rDis = rDis
    self.rEf = rEf
    self.cP = cP
    self.Vx = Vx
    self.Vz = self.x_to_z(Vx)
    self.SRp = [0,1,2,3, 5,6,7,4, 10,11,8,9, 15,12,13,14]
    self.SRpv = [0,1,2,3, 7,4,5,6, 10,11,8,9, 13,14,15,12]
    self.hTable = [8,9,10,11, 12,13,14,15,  2,0,4,7, 6,3,5,1]
    self.name = './v3_SKINNYe-{}-{}_{}-{}+{}+{}_{}'.format(
      self.b_size, self.k_size, rEb+rDis+rEf, rEb, rDis, rEf, cP
    )
    self.model = gp.Model(self.name)

  def ib_model(self):
    '''
    ==========================================================================================
    (Begin) - Upper differential propagation

    Note:
    X -(SC&AC)-> Y -(ATK)-> Z -(SR)-> W -(MC)-> X^{r+1}
    '''

    uR = self.rEb+self.rDis
    uDX  = self.model.addVars(range(1, uR+1), 16, vtype = GRB.BINARY, name = 'uDX')
    udx  = self.model.addVars(range(1, uR+1), 16, vtype = GRB.BINARY, name = 'udx')
    uDY  = self.model.addVars(range(1, uR),   16, vtype = GRB.BINARY, name = 'uDY')
    udy  = self.model.addVars(range(1, uR),   16, vtype = GRB.BINARY, name = 'udy')
    uDZ  = self.model.addVars(range(1, uR),   16, vtype = GRB.BINARY, name = 'uDZ') 
    udz  = self.model.addVars(range(1, uR),   16, vtype = GRB.BINARY, name = 'udz') 
    uDW  = self.model.addVars(uR,             16, vtype = GRB.BINARY, name = 'uDW') # W0 = Plaintext
    udw  = self.model.addVars(uR,             16, vtype = GRB.BINARY, name = 'udw') # w0 = Plaintext
    ustk = self.model.addVars(uR,     16, vtype = GRB.BINARY, name = 'ustk')

    '''Basic Constrs:'''
    # Remove dx[r,i] > DX[r,i]
    for i in range(16): 
      for r in range(1, uR+1):
        self.model.addConstr(uDX[r, i] >= udx[r, i], name = 'uBasic')
      for r in range(1, uR):
        self.model.addConstr(uDY[r, i] >= udy[r, i], name = 'uBasic')
        self.model.addConstr(uDZ[r, i] >= udz[r, i], name = 'uBasic')
      for r in range(uR):
        self.model.addConstr(uDW[r, i] >= udw[r, i], name = 'uBasic')

    '''Operation: SC'''
    # In Extension
    for r in range(1, self.rEb+1):
      for i in range(16):
        self.model.addConstr(uDX[r,i] == uDY[r,i], name='ueSC')
        self.model.addConstr(udx[r,i] == uDY[r,i], name='ueSC')
    # In distinguisher
    for r in range(self.rEb+1, uR):
      for i in range(16):
        self.model.addConstr(uDY[r,i] == uDX[r,i], name='udSC')
        self.model.addConstr(udy[r,i] == uDX[r,i], name='udSC')


    '''(*ri-(r0)-r1)Operation: ATK'''
  
    ''' 
    ----------------------------------
    Equivalent Key for the first round
    ----------------------------------
    '''
    eqk = self.model.addVars(16, vtype = GRB.BINARY, name = 'eqk')

    for c in [[0,4,8,12],[1,5,9,13],[2,6,10,14],[3,7,11,15]]:
        c0, c1, c2, c3 = c[0], c[1], c[2], c[3]
        tc0, tc1, tc2, tc3 = self.SRpv[c[0]], self.SRpv[c[1]], self.SRpv[c[2]], self.SRpv[c[3]]
        
        # --- POSITION -- 0: eqk[c0] = stk1[c3] + eqk[c3] ---
        self.model.addConstr(eqk[c0] == ustk[0,tc0], name = 'eqk0')

        # --- POSITION -- 1: eqk[c1] = stk1[c0] ---
        self.model.addConstr(eqk[c1] == ustk[0,tc0], name = 'eqk1')
        
        # --- POSITION -- 2: eqk[c2] = stk1[c1] + stk[c2] ---
        self.model.addConstr(eqk[c2] == ustk[0,tc1], name = 'eqk2')

        # --- POSITION -- 2: eqk[c3] = stk1[c0] + stk[c2] ---
        self.model.addConstr(eqk[c3] == ustk[0,tc0], name = 'eqk2')


    # Add eqk in the first round
    for i in range(16):
      self.model.addConstr(udw[0,i] == udx[1,i], name='uAEqk')
      self.model.addConstr(eqk[i] + uDW[0,i] - uDX[1,i] >= 0, name='uAEqk')
      self.model.addConstr(eqk[i] - uDW[0,i] + uDX[1,i] >= 0, name='uAEqk')
      self.model.addConstr(-eqk[i] + uDW[0,i] + uDX[1,i] - udx[1,i] >= 0, name='uAEqk')
      self.model.addConstr(uDW[0,i] >= udx[1,i], name='uAEqk')
    # # In ri
    # for r in range(1, self.ri):
    #   for i in range(8):
    #     self.model.addConstr(udy[r,i] == udz[r,i], name='uATKri')
    #     self.model.addConstr(ustk[r,i] + uDY[r,i] - uDZ[r,i] >= 0, name='uATKri')
    #     self.model.addConstr(ustk[r,i] - uDY[r,i] + uDZ[r,i] >= 0, name='uATKri')
    #     self.model.addConstr(-ustk[r,i] + uDY[r,i] + uDZ[r,i] - udz[r,i] >= 0, name='uATKri')
    #     self.model.addConstr(uDY[r,i] >= udz[r,i], name='uATKri')
    #   for i in range(8,16):
    #     self.model.addConstr(uDY[r,i] == uDZ[r,i], name='uATKri')
    #     self.model.addConstr(udy[r,i] == udz[r,i], name='uATKri')
    # # In r1
    # for r in range(self.ri, uR):
    #   tR = r + self.r0
    #   for i in range(8):
    #     self.model.addConstr(udy[r,i] == udz[r,i], name='uATKr1')
    #     self.model.addConstr(ustk[tR,i] + uDY[r,i] - uDZ[r,i] >= 0, name='uATKr1')
    #     self.model.addConstr(ustk[tR,i] - uDY[r,i] + uDZ[r,i] >= 0, name='uATKr1')
    #     self.model.addConstr(-ustk[tR,i] + uDY[r,i] + uDZ[r,i] - udz[r,i] >= 0, name='uATKr1')
    #     self.model.addConstr(uDY[r,i] >= udz[r,i], name='uATKr1')
    #   for i in range(8,16):
    #     self.model.addConstr(uDY[r,i] == uDZ[r,i], name='uATKr1')
    #     self.model.addConstr(udy[r,i] == udz[r,i], name='uATKr1')
    for r in range(1, uR):
      for i in range(8):
        self.model.addConstr(udy[r,i] == udz[r,i], name='uATK')
        self.model.addConstr(ustk[r,i] + uDY[r,i] - uDZ[r,i] >= 0, name='uATK')
        self.model.addConstr(ustk[r,i] - uDY[r,i] + uDZ[r,i] >= 0, name='uATK')
        self.model.addConstr(-ustk[r,i] + uDY[r,i] + uDZ[r,i] - udz[r,i] >= 0, name='uATK')
        self.model.addConstr(uDY[r,i] >= udz[r,i], name='uATK')
      for i in range(8,16):
        self.model.addConstr(uDY[r,i] == uDZ[r,i], name='uATK')
        self.model.addConstr(udy[r,i] == udz[r,i], name='uATK')


    '''(*)Operation: SR''' 
    for r in range(1,uR):
      for i in range(16):
        self.model.addConstr(uDZ[r,i] == uDW[r,self.SRp[i]], name = 'uSR')
        self.model.addConstr(udz[r,i] == udw[r,self.SRp[i]], name = 'uSR')

    '''Operation: MC'''
    # ************ In extension (backword) ************
    for r in range(1, self.rEb):
      for c in [[0,4,8,12],[1,5,9,13],[2,6,10,14],[3,7,11,15]]:
        c0, c1, c2, c3 = c[0], c[1], c[2], c[3]

        # --- POSITION -- 0: W[c0] = X[c1] ---
        self.model.addConstr(udw[r,c0] == udx[r+1,c1], name = 'ueMC0')

        self.model.addConstr(uDW[r,c0] == uDX[r+1,c1], name = 'ueMC0')

        # --- POSITION -- 1: W[c1] = W[c2] + X[c2] ---
        self.model.addConstr(udw[r,c1] >= udw[r,c2], name = 'ueMC1')
        self.model.addConstr(udw[r,c1] >= udx[r+1,c2], name = 'ueMC1')
        self.model.addConstr(udw[r,c1] <= udw[r,c2] + udx[r+1,c2], name = 'ueMC1')

        self.model.addConstr(uDW[r,c1] + uDW[r,c2] + uDX[r+1,c2] >= 2*uDW[r,c1], name = 'ueMC1')
        self.model.addConstr(uDW[r,c1] + uDW[r,c2] + uDX[r+1,c2] >= 2*uDW[r,c2], name = 'ueMC1')
        self.model.addConstr(uDW[r,c1] + uDW[r,c2] + uDX[r+1,c2] >= 2*uDX[r+1,c2], name = 'ueMC1')
        self.model.addConstr(uDW[r,c1] <= uDW[r,c2] + uDX[r+1,c2], name = 'ueMC1')

        # --- POSITION -- 2: W[c2] = X[c1] + X[c3] ---
        self.model.addConstr(udw[r,c2] >= udx[r+1,c1], name = 'ueMC2')
        self.model.addConstr(udw[r,c2] >= udx[r+1,c3], name = 'ueMC2')
        self.model.addConstr(udw[r,c2] <= udx[r+1,c1] + udx[r+1,c3], name = 'ueMC2')

        self.model.addConstr(uDW[r,c2] + uDX[r+1,c1] + uDX[r+1,c3] >= 2*uDW[r,c2], name = 'ueMC2')
        self.model.addConstr(uDW[r,c2] + uDX[r+1,c1] + uDX[r+1,c3] >= 2*uDX[r+1,c1], name = 'ueMC2')
        self.model.addConstr(uDW[r,c2] + uDX[r+1,c1] + uDX[r+1,c3] >= 2*uDX[r+1,c3], name = 'ueMC2')
        self.model.addConstr(uDW[r,c2] <= uDX[r+1,c1] + uDX[r+1,c3], name = 'ueMC2')

        # --- POSITION -- 3: W[c3] = X[c0] + X[c3] ---
        self.model.addConstr(udw[r,c3] >= udx[r+1,c0], name = 'ueMC3')
        self.model.addConstr(udw[r,c3] >= udx[r+1,c3], name = 'ueMC3')
        self.model.addConstr(udw[r,c3] <= udx[r+1,c0] + udx[r+1,c3], name = 'ueMC3')

        self.model.addConstr(uDW[r,c3] + uDX[r+1,c0] + uDX[r+1,c3] >= 2*uDW[r,c3], name = 'ueMC3')
        self.model.addConstr(uDW[r,c3] + uDX[r+1,c0] + uDX[r+1,c3] >= 2*uDX[r+1,c0], name = 'ueMC3')
        self.model.addConstr(uDW[r,c3] + uDX[r+1,c0] + uDX[r+1,c3] >= 2*uDX[r+1,c3], name = 'ueMC3')
        self.model.addConstr(uDW[r,c3] <= uDX[r+1,c0] + uDX[r+1,c3], name = 'ueMC3')
    
    # ************ In distinguisher (forward) ************
    for r in range(self.rEb, uR):
      for c in [[0,4,8,12],[1,5,9,13],[2,6,10,14],[3,7,11,15]]:
          c0, c1, c2, c3 = c[0], c[1], c[2], c[3]

          # --- POSITION -- 0: X[c0] = X[c3] + W[c3] ---
          self.model.addConstr(udx[r+1,c0] >= udw[r,c3], name = 'udMC0')
          self.model.addConstr(udx[r+1,c0] >= udx[r+1,c3], name = 'udMC0')
          self.model.addConstr(udx[r+1,c0] <= udw[r,c3] + udx[r+1,c3], name = 'udMC0')

          self.model.addConstr(uDX[r+1,c0] + uDW[r,c3] + uDX[r+1,c3] >= 2*uDX[r+1,c0], name = 'udMC0')
          self.model.addConstr(uDX[r+1,c0] + uDW[r,c3] + uDX[r+1,c3] >= 2*uDW[r,c3], name = 'udMC0')
          self.model.addConstr(uDX[r+1,c0] + uDW[r,c3] + uDX[r+1,c3] >= 2*uDX[r+1,c3], name = 'udMC0')
          self.model.addConstr(uDX[r+1,c0] <= uDW[r,c3] + uDX[r+1,c3], name = 'udMC0')

          # --- POSITION -- 1: X[c1] = W[c0] ---
          self.model.addConstr(udx[r+1,c1] == udw[r,c0], name = 'udMC1')

          self.model.addConstr(uDX[r+1,c1] == uDW[r,c0], name = 'udMC1')

          # --- POSITION -- 2: X[c2] = X[c1] + W[c2] ---
          self.model.addConstr(udx[r+1,c2] >= udw[r,c1], name = 'udMC2')
          self.model.addConstr(udx[r+1,c2] >= udw[r,c2], name = 'udMC2')
          self.model.addConstr(udx[r+1,c2] <= udw[r,c1] + udw[r,c2], name = 'udMC2')

          self.model.addConstr(uDX[r+1,c2] + uDW[r,c1] + uDW[r,c2] >= 2*uDX[r+1,c2], name = 'udMC2')
          self.model.addConstr(uDX[r+1,c2] + uDW[r,c1] + uDW[r,c2] >= 2*uDW[r,c1], name = 'udMC2')
          self.model.addConstr(uDX[r+1,c2] + uDW[r,c1] + uDW[r,c2] >= 2*uDW[r,c2], name = 'udMC2')
          self.model.addConstr(uDX[r+1,c2] <= uDW[r,c1] + uDW[r,c2], name = 'udMC2')

          # --- POSITION -- 3: X[c3] = W[c0] + W[c2] ---
          self.model.addConstr(udx[r+1,c3] >= udw[r,c0], name = 'udMC3')
          self.model.addConstr(udx[r+1,c3] >= udw[r,c2], name = 'udMC3')
          self.model.addConstr(udx[r+1,c3] <= udw[r,c0] + udw[r,c2], name = 'udMC3')

          self.model.addConstr(uDX[r+1,c3] + uDW[r,c0] + uDW[r,c2] >= 2*uDX[r+1,c3], name = 'udMC3')
          self.model.addConstr(uDX[r+1,c3] + uDW[r,c0] + uDW[r,c2] >= 2*uDW[r,c0], name = 'udMC3')
          self.model.addConstr(uDX[r+1,c3] + uDW[r,c0] + uDW[r,c2] >= 2*uDW[r,c2], name = 'udMC3')
          self.model.addConstr(uDX[r+1,c3] <= uDW[r,c0] + uDW[r,c2], name = 'udMC3')
    
    ''' 
    ------------------------------
    Key Schedule 
    ------------------------------
    '''
    uLANE = self.model.addVars(16, vtype = GRB.BINARY, name = 'uLANE')

    trans_pos = [1,7,0,5,2,6,4,3,9,15,8,13,10,14,12,11]

    for i in range(16):
      lane_pos = i
      total_stack_sum = 0
      for r in range(uR):
        self.model.addConstr(uLANE[i] >= ustk[r, lane_pos], name = 'uKS')
        total_stack_sum += ustk[r, lane_pos]
        lane_pos = self.hTable[lane_pos]
      # math.ceil((uR + self.r0)/16) = cancellations from key schedule
      self.model.addConstr((uR)*uLANE[i] - total_stack_sum <= 
                            (self.Vs - 1)*math.ceil((uR)/30), name = 'uKS')

    # make Key reuse
    for r in range(30, uR):
        for i in range(16):
          self.model.addConstr(ustk[r-30,i] == ustk[r,trans_pos[i]], name = 'uKS')
    
    '''
    ==========================================================================================
    (End) - Upper propagation
    '''

    '''
    ==========================================================================================
    (Begin) - Lower differential propagation

    Note:
    X -(SC&AC)-> Y -(ATK)-> Z -(SR)-> W -(MC)-> X^{r+1}
    '''
    lR = self.rDis + self.rEf
    lDX  = self.model.addVars(range(self.rEb+1, self.rEb+lR+1), 16, vtype = GRB.BINARY, name = 'lDX')
    ldx  = self.model.addVars(range(self.rEb+1, self.rEb+lR+1), 16, vtype = GRB.BINARY, name = 'ldx')
    lDY  = self.model.addVars(range(self.rEb, self.rEb+lR), 16, vtype = GRB.BINARY, name = 'lDY')
    ldy  = self.model.addVars(range(self.rEb, self.rEb+lR), 16, vtype = GRB.BINARY, name = 'ldy')
    lDZ  = self.model.addVars(range(self.rEb, self.rEb+lR), 16, vtype = GRB.BINARY, name = 'lDZ') # DZ = Ciphertext
    ldz  = self.model.addVars(range(self.rEb, self.rEb+lR), 16, vtype = GRB.BINARY, name = 'ldz') # dz = Ciphertext
    lDW  = self.model.addVars(range(self.rEb, self.rEb+lR), 16, vtype = GRB.BINARY, name = 'lDW')
    ldw  = self.model.addVars(range(self.rEb, self.rEb+lR), 16, vtype = GRB.BINARY, name = 'ldw')
    lstk = self.model.addVars(range(self.rEb, self.rEb+lR), 16, vtype = GRB.BINARY, name = 'lstk')

    '''Basic Constrs:'''
    # Remove dx[r,i] > DX[r,i]
    for r in range(self.rEb+1,self.rEb+lR):
      for i in range(16):
        self.model.addConstr(lDX[r, i] >= ldx[r, i], name='lBasic')
        self.model.addConstr(lDY[r, i] >= ldy[r, i], name='lBasic')
        self.model.addConstr(lDZ[r, i] >= ldz[r, i], name='lBasic')
        self.model.addConstr(lDW[r, i] >= ldw[r, i], name='lBasic')
    for i in range(16):
      self.model.addConstr(lDY[self.rEb, i] >= ldy[self.rEb, i], name='lBasic')
      self.model.addConstr(lDZ[self.rEb, i] >= ldz[self.rEb, i], name='lBasic')
      self.model.addConstr(lDW[self.rEb, i] >= ldw[self.rEb, i], name='lBasic')
    for i in range(16):
      self.model.addConstr(lDX[self.rEb+lR, i] >= ldx[self.rEb+lR, i], name='lBasic')

    '''Operation: SC'''
    # In distinglisher
    for r in range(self.rEb+1, self.rEb+self.rDis):
      for i in range(16):
        self.model.addConstr(lDX[r,i] == lDY[r,i], name='ldSC')
        self.model.addConstr(ldx[r,i] == lDY[r,i], name='ldSC')
    # In Extension
    for r in range(self.rEb+self.rDis, self.rEb+lR):
      for i in range(16):
        self.model.addConstr(lDY[r,i] == lDX[r,i], name='leSC')
        self.model.addConstr(ldy[r,i] == lDX[r,i], name='leSC')
    
    '''(*)Operation: ATK'''
    for r in range(self.rEb,self.rEb+lR):
      for i in range(8):
        # tlr = self.r0 + r
        self.model.addConstr(ldy[r,i] == ldz[r,i], name='lATK')
        self.model.addConstr(lstk[r,i] + lDY[r,i] - lDZ[r,i] >= 0, name='lATK')
        self.model.addConstr(lstk[r,i] - lDY[r,i] + lDZ[r,i] >= 0, name='lATK')
        self.model.addConstr(-lstk[r,i] + lDY[r,i] + lDZ[r,i] - ldz[r,i] >= 0, name='lATK')
        self.model.addConstr(lDY[r,i] >= ldz[r,i], name='lATK')
      for i in range(8,16):
        self.model.addConstr(lDY[r,i] == lDZ[r,i], name='lATK')
        self.model.addConstr(ldy[r,i] == ldz[r,i], name='lATK')

    '''(*)Operation: SR''' 
    for r in range(self.rEb,self.rEb+lR):
      for i in range(16):
        self.model.addConstr(lDZ[r,i] == lDW[r,self.SRp[i]], name='lSR')
        self.model.addConstr(ldz[r,i] == ldw[r,self.SRp[i]], name='lSR')

    '''Operation: MC'''
    # ************ In distinguisher (backword) ************
    for r in range(self.rEb, self.rEb+self.rDis):
      for c in [[0,4,8,12],[1,5,9,13],[2,6,10,14],[3,7,11,15]]:
          c0, c1, c2, c3 = c[0], c[1], c[2], c[3]

          # --- POSITION -- 0: W[c0] = X[c1] ---
          self.model.addConstr(ldw[r,c0] == ldx[r+1,c1], name = 'ldMC0')

          self.model.addConstr(lDW[r,c0] == lDX[r+1,c1], name = 'ldMC0')

          # --- POSITION -- 1: W[c1] = W[c2] + X[c2] ---
          self.model.addConstr(ldw[r,c1] >= ldw[r,c2], name = 'ldMC1')
          self.model.addConstr(ldw[r,c1] >= ldx[r+1,c2], name = 'ldMC1')
          self.model.addConstr(ldw[r,c1] <= ldw[r,c2] + ldx[r+1,c2], name = 'ldMC1')

          self.model.addConstr(lDW[r,c1] + lDW[r,c2] + lDX[r+1,c2] >= 2*lDW[r,c1], name = 'ldMC1')
          self.model.addConstr(lDW[r,c1] + lDW[r,c2] + lDX[r+1,c2] >= 2*lDW[r,c2], name = 'ldMC1')
          self.model.addConstr(lDW[r,c1] + lDW[r,c2] + lDX[r+1,c2] >= 2*lDX[r+1,c2], name = 'ldMC1')
          self.model.addConstr(lDW[r,c1] <= lDW[r,c2] + lDX[r+1,c2], name = 'ldMC1')

          # --- POSITION -- 2: W[c2] = X[c1] + X[c3] ---
          self.model.addConstr(ldw[r,c2] >= ldx[r+1,c1], name = 'ldMC2')
          self.model.addConstr(ldw[r,c2] >= ldx[r+1,c3], name = 'ldMC2')
          self.model.addConstr(ldw[r,c2] <= ldx[r+1,c1] + ldx[r+1,c3], name = 'ldMC2')

          self.model.addConstr(lDW[r,c2] + lDX[r+1,c1] + lDX[r+1,c3] >= 2*lDW[r,c2], name = 'ldMC2')
          self.model.addConstr(lDW[r,c2] + lDX[r+1,c1] + lDX[r+1,c3] >= 2*lDX[r+1,c1], name = 'ldMC2')
          self.model.addConstr(lDW[r,c2] + lDX[r+1,c1] + lDX[r+1,c3] >= 2*lDX[r+1,c3], name = 'ldMC2')
          self.model.addConstr(lDW[r,c2] <= lDX[r+1,c1] + lDX[r+1,c3], name = 'ldMC2')

          # --- POSITION -- 3: W[c3] = X[c0] + X[c3] ---
          self.model.addConstr(ldw[r,c3] >= ldx[r+1,c0], name = 'ldMC3')
          self.model.addConstr(ldw[r,c3] >= ldx[r+1,c3], name = 'ldMC3')
          self.model.addConstr(ldw[r,c3] <= ldx[r+1,c0] + ldx[r+1,c3], name = 'ldMC3')

          self.model.addConstr(lDW[r,c3] + lDX[r+1,c0] + lDX[r+1,c3] >= 2*lDW[r,c3], name = 'ldMC3')
          self.model.addConstr(lDW[r,c3] + lDX[r+1,c0] + lDX[r+1,c3] >= 2*lDX[r+1,c0], name = 'ldMC3')
          self.model.addConstr(lDW[r,c3] + lDX[r+1,c0] + lDX[r+1,c3] >= 2*lDX[r+1,c3], name = 'ldMC3')
          self.model.addConstr(lDW[r,c3] <= lDX[r+1,c0] + lDX[r+1,c3], name = 'ldMC3')

    # ************ In extension (forward)  ************
    for r in range(self.rEb+self.rDis, self.rEb+lR):
      for c in [[0,4,8,12],[1,5,9,13],[2,6,10,14],[3,7,11,15]]:
          c0, c1, c2, c3 = c[0], c[1], c[2], c[3]

          # --- POSITION -- 0: X[c0] = X[c3] + W[c3] ---
          self.model.addConstr(ldx[r+1,c0] >= ldw[r,c3], name = 'leMC0')
          self.model.addConstr(ldx[r+1,c0] >= ldx[r+1,c3], name = 'leMC0')
          self.model.addConstr(ldx[r+1,c0] <= ldw[r,c3] + ldx[r+1,c3], name = 'leMC0')

          self.model.addConstr(lDX[r+1,c0] + lDW[r,c3] + lDX[r+1,c3] >= 2*lDX[r+1,c0], name = 'leMC0')
          self.model.addConstr(lDX[r+1,c0] + lDW[r,c3] + lDX[r+1,c3] >= 2*lDW[r,c3], name = 'leMC0')
          self.model.addConstr(lDX[r+1,c0] + lDW[r,c3] + lDX[r+1,c3] >= 2*lDX[r+1,c3], name = 'leMC0')
          self.model.addConstr(lDX[r+1,c0] <= lDW[r,c3] + lDX[r+1,c3], name = 'leMC0')

          # --- POSITION -- 1: X[c1] = W[c0] ---
          self.model.addConstr(ldx[r+1,c1] == ldw[r,c0], name = 'leMC1')

          self.model.addConstr(lDX[r+1,c1] == lDW[r,c0], name = 'leMC1')

          # --- POSITION -- 2: X[c2] = X[c1] + W[c2] ---
          self.model.addConstr(ldx[r+1,c2] >= ldw[r,c1], name = 'leMC2')
          self.model.addConstr(ldx[r+1,c2] >= ldw[r,c2], name = 'leMC2')
          self.model.addConstr(ldx[r+1,c2] <= ldw[r,c1] + ldw[r,c2], name = 'leMC2')

          self.model.addConstr(lDX[r+1,c2] + lDW[r,c1] + lDW[r,c2] >= 2*lDX[r+1,c2], name = 'leMC2')
          self.model.addConstr(lDX[r+1,c2] + lDW[r,c1] + lDW[r,c2] >= 2*lDW[r,c1], name = 'leMC2')
          self.model.addConstr(lDX[r+1,c2] + lDW[r,c1] + lDW[r,c2] >= 2*lDW[r,c2], name = 'leMC2')
          self.model.addConstr(lDX[r+1,c2] <= lDW[r,c1] + lDW[r,c2], name = 'leMC2')

          # --- POSITION -- 3: X[c3] = W[c0] + W[c2] ---
          self.model.addConstr(ldx[r+1,c3] >= ldw[r,c0], name = 'leMC3')
          self.model.addConstr(ldx[r+1,c3] >= ldw[r,c2], name = 'leMC3')
          self.model.addConstr(ldx[r+1,c3] <= ldw[r,c0] + ldw[r,c2], name = 'leMC3')

          self.model.addConstr(lDX[r+1,c3] + lDW[r,c0] + lDW[r,c2] >= 2*lDX[r+1,c3], name = 'leMC3')
          self.model.addConstr(lDX[r+1,c3] + lDW[r,c0] + lDW[r,c2] >= 2*lDW[r,c0], name = 'leMC3')
          self.model.addConstr(lDX[r+1,c3] + lDW[r,c0] + lDW[r,c2] >= 2*lDW[r,c2], name = 'leMC3')
          self.model.addConstr(lDX[r+1,c3] <= lDW[r,c0] + lDW[r,c2], name = 'leMC3')


    ''' 
    ------------------------------
    Key Schedule 
    ------------------------------
    '''
    lLANE = self.model.addVars(16, vtype = GRB.BINARY, name = 'lLANE')

    for i in range(16):
      lane_pos = i
      kSUM = 0
      for r in range(self.rEb, self.rEb+lR):
        self.model.addConstr(lLANE[i] >= lstk[r, lane_pos], name = 'lKS')
        kSUM += lstk[r, lane_pos]
        lane_pos = self.hTable[lane_pos]
      # math.ceil((self.rEb + self.r0)/16) = cancellations from key schedule
      # Since the statement “lR < 30” and the lR span is only in r1, there are only (Vs-1) cancelations
      self.model.addConstr(kSUM >= (lR)*lLANE[i] - (self.Vs - 1), name = 'lKS')
      # No key reusing

    '''
    ==========================================================================================
    (End) - Lower differential propagation
    '''


    ''' 
    ============================================
            Guess-and-Determine (u)
    ============================================
    '''
    # NOTE: Assert rEb <= ri and rEb+rEu >= ri
    ugEQK = self.model.addVars(range(self.rEb),      16, vtype = GRB.BINARY, name = 'ugEQK')
    uGstk = self.model.addVars(range(self.rEb),      16, vtype = GRB.BINARY, name = 'uGstk')
    for r in range(self.rEb): 
      for i in range(8,16):
        self.model.addConstr(uGstk[r,i] == 0, name = 'uGstkL')
    uDetX = self.model.addVars(range(1, self.rEb+1), 16, vtype = GRB.BINARY, name = 'uDetX')
    uDetY = self.model.addVars(range(1, self.rEb+1), 16, vtype = GRB.BINARY, name = 'uDetY')
    uDetZ = self.model.addVars(range(1, self.rEb),   16, vtype = GRB.BINARY, name = 'uDetZ')
    uDetW = self.model.addVars(range(0, self.rEb),   16, vtype = GRB.BINARY, name = 'uDetW')
    uFrSC = self.model.addVars(range(1, self.rEb+1), 16, vtype = GRB.BINARY, name = 'uFrSC')
    uFrMC = self.model.addVars(range(2, self.rEb+1), 16, vtype = GRB.BINARY, name = 'uFrMC')
    #-----------------------------------------------------------------------------------------
    ugeqka = self.model.addVars(range(self.rEb),      16, vtype = GRB.BINARY, name = 'ugeqka')
    ugstka = self.model.addVars(range(self.rEb),      16, vtype = GRB.BINARY, name = 'ugstka')
    for r in range(self.rEb): 
      for i in range(8,16):
        self.model.addConstr(ugstka[r,i] == 0, 'ugstkLa')
    udetXa = self.model.addVars(range(1, self.rEb+1), 16, vtype = GRB.BINARY, name = 'udetXa')
    udetYa = self.model.addVars(range(1, self.rEb+1), 16, vtype = GRB.BINARY, name = 'udetYa')
    udetZa = self.model.addVars(range(1, self.rEb),   16, vtype = GRB.BINARY, name = 'udetZa')
    udetWa = self.model.addVars(range(0, self.rEb),   16, vtype = GRB.BINARY, name = 'udetWa')
    ufrSCa = self.model.addVars(range(1, self.rEb+1), 16, vtype = GRB.BINARY, name = 'ufrSCa')
    ufrMCa = self.model.addVars(range(2, self.rEb+1), 16, vtype = GRB.BINARY, name = 'ufrMCa')

    self.model.addConstr(sum(uDetW[0,i] for i in range(16)) == 16, name = 'VP') # Set plaintext is determined
    for i in range(16): self.model.addConstr(uDetX[1,i] == uDetW[0,i], name = 'W0X0a')
    # ------------------------------------------------------------------------
    self.model.addConstr(sum(udetWa[0,i] for i in range(16)) == 16, name = 'VPa') # Set plaintext is determined
    for i in range(16): self.model.addConstr(udetXa[1,i] == uDetW[0,i], name = 'W0X0a')

    # ----------------- Guessing Eqk -----------------------
    # NOTE: That is for the key bridge.
    for r in range(self.rEb):
      for c in [[0,4,8,12],[1,5,9,13],[2,6,10,14],[3,7,11,15]]:
          c0, c1, c2, c3 = c[0], c[1], c[2], c[3]
          tc0, tc1 = self.SRpv[c[0]], self.SRpv[c[1]]

          # POSITION -- 0: ugEQK[0,1,2,3 (c0)] = uGstk[0,1,2,3 -> (tc0)]
          self.model.addConstr(ugEQK[r,c0] == uGstk[r,tc0], name = 'geqk0')
          self.model.addConstr(ugeqka[r,c0] == ugstka[r,tc0], name = 'geqk0A')

          # POSITION -- 1: ugEQK[4,5,6,7 (c1)] = uGstk[0,1,2,3 -> (tc0)] 
          self.model.addConstr(ugEQK[r,c1] == uGstk[r,tc0], name = 'geqk1')
          self.model.addConstr(ugeqka[r,c1] == ugstka[r,tc0], name = 'geqk1A')

          # POSITION -- 2: ugEQK[8,9,10,11] = uGstk[7,4,5,6 -> (tc1)] 
          self.model.addConstr(ugEQK[r,c2] == uGstk[r,tc1], name = 'geqk2')
          self.model.addConstr(ugeqka[r,c2] == ugstka[r,tc1], name = 'geqk2A')

          # POSITION -- 2: eqk[c3] = stk1[12,13,14,15] + stk[0,1,2,3 -> (tc0)] 
          self.model.addConstr(ugEQK[r,c3] == uGstk[r,tc0], name = 'geqk3')
          self.model.addConstr(ugeqka[r,c3] == ugstka[r,tc0], name = 'geqk3A')
    # ----------------- Guessing Eqk -----------------------

    ''' KR: Determine and Filter in SC '''
    for r in range(1, self.rEb+1):
      for i in range(16):
        # 'Det.' propagation via SC (AND ADD EQK)
        self.model.addConstr(uDetY[r,i] <= uDetX[r,i], name = 'udetSC+ek')
        self.model.addConstr(uDetY[r,i] <= ugEQK[r-1,i], name = 'udetSC+ek')
        self.model.addConstr( - uDetX[r,i] - ugEQK[r-1,i] + uDetY[r,i] + 1 >= 0, name = 'udetSC+ek')
        # ------------------------------------------------------------------------------------------
        self.model.addConstr(udetYa[r,i] <= udetXa[r,i], name = 'udetSC+ekA')
        self.model.addConstr(udetYa[r,i] <= ugeqka[r-1,i], name = 'udetSC+ekA')
        self.model.addConstr( - udetXa[r,i] - ugeqka[r-1,i] + udetYa[r,i] + 1 >= 0, name = 'udetSC+ekA')

        # Filter obtain from SC
        self.model.addConstr(uDY[r,i] - udy[r,i] - uFrSC[r,i] >= 0, name = 'ufrSC')
        self.model.addConstr(uDetY[r,i] >= uFrSC[r,i], name = 'ufrSC')
        self.model.addConstr(- uDetY[r,i] - uDY[r,i] + udy[r,i] + uFrSC[r,i] >= -1, name = 'ufrSC')
        # ------------------------------------------------------------------------------------------
        self.model.addConstr(uDY[r,i] - udy[r,i] - ufrSCa[r,i] >= 0, name = 'ufrSCa')
        self.model.addConstr(udetYa[r,i] >= ufrSCa[r,i], name = 'ufrSCa')
        self.model.addConstr(- udetYa[r,i] - uDY[r,i] + udy[r,i] + ufrSCa[r,i] >= -1, name = 'ufrSCa')

    ''' 'Det.' trans Y to Z '''
    for r in range(1, self.rEb):
      for i in range(16):
        self.model.addConstr(uDetY[r,i] == uDetZ[r,i], 'udetY2Z')
        # -------------------------------------------------------
        self.model.addConstr(udetYa[r,i] == udetZa[r,i], 'udety2zA')

    ''' KR: 'Det.' propagation in SR '''
    for r in range(1, self.rEb):
      for i in range(16): 
        self.model.addConstr(uDetW[r,i] == uDetZ[r,self.SRpv[i]], name = 'udetSR')
        # -------------------------------------------------------------------------
        self.model.addConstr(udetWa[r,i] == udetZa[r,self.SRpv[i]], name = 'udetsrA')

    ''' KR: 'Det.' propagation in MC '''
    for r in range(1, self.rEb):
      for c in [[0,4,8,12],[1,5,9,13],[2,6,10,14],[3,7,11,15]]:
          c0, c1, c2, c3 = c[0], c[1], c[2], c[3]

          # ===============================================================================================
          # POSITION -- 0: DetX[c0] = DetW[3] & DetX[c3]
          self.model.addConstr(uDetX[r+1,c0] <= uDetW[r,c3], name = 'udetMC0')
          self.model.addConstr(uDetX[r+1,c0] <= uDetX[r+1,c3], name = 'udetMC0')
          self.model.addConstr( - uDetW[r,c3] - uDetX[r+1,c3] + uDetX[r+1,c0] + 1 >= 0, name = 'udetMC0')
          # --------------------------------------------------------------------------------------------
          self.model.addConstr(udetXa[r+1,c0] <= udetWa[r,c3], name = 'udetmc0A')
          self.model.addConstr(udetXa[r+1,c0] <= udetXa[r+1,c3], name = 'udetmc0A')
          self.model.addConstr( - udetWa[r,c3] - udetXa[r+1,c3] + udetXa[r+1,c0] + 1 >= 0, name = 'udetmc0A')

          # POSITION -- 1: DetX[c1] = DetW[c0]
          self.model.addConstr(uDetX[r+1,c1] == uDetW[r,c0], name = 'udetMC1')
          # ---------------------------------------------------------------------
          self.model.addConstr(udetXa[r+1,c1] == udetWa[r,c0], name = 'udetmc1A')

          # POSITION -- 2: DetX[c2] = DetW[c1] & DetW[c2]
          self.model.addConstr(uDetX[r+1,c2] <= uDetW[r,c1], name = 'udetMC2')
          self.model.addConstr(uDetX[r+1,c2] <= uDetW[r,c2], name = 'udetMC2')
          self.model.addConstr( - uDetW[r,c1] - uDetW[r,c2] + uDetX[r+1,c2] + 1 >= 0, name = 'udetMC2')
          # --------------------------------------------------------------------------------------------
          self.model.addConstr(udetXa[r+1,c2] <= udetWa[r,c1], name = 'udetmc2A')
          self.model.addConstr(udetXa[r+1,c2] <= udetWa[r,c2], name = 'udetmc2A')
          self.model.addConstr( - udetWa[r,c1] - udetWa[r,c2] + udetXa[r+1,c2] + 1 >= 0, name = 'udetmc2A')

          # POSITION -- 3: DetX[c3] = DetW[c0] & DetW[c2]
          self.model.addConstr(uDetX[r+1,c3] <= uDetW[r,c0], name = 'udetMC3')
          self.model.addConstr(uDetX[r+1,c3] <= uDetW[r,c2], name = 'udetMC3')
          self.model.addConstr( - uDetW[r,c0] - uDetW[r,c2] + uDetX[r+1,c3] + 1 >= 0, name = 'udetMC3')
          # --------------------------------------------------------------------------------------------
          self.model.addConstr(udetXa[r+1,c3] <= udetWa[r,c0], name = 'udetmc3A')
          self.model.addConstr(udetXa[r+1,c3] <= udetWa[r,c2], name = 'udetmc3A')
          self.model.addConstr( - udetWa[r,c0] - udetWa[r,c2] + udetXa[r+1,c3] + 1 >= 0, name = 'udetmc3A')
          # ===============================================================================================

    ''' KR: 'Filter' obtain from MC '''
    for r in range(2,self.rEb+1):
      for c in [[0,4,8,12],[1,5,9,13],[2,6,10,14],[3,7,11,15]]:
          c0, c1, c2, c3 = c[0], c[1], c[2], c[3]

          # POSITION -- 0: DetX[c0] = DetW[3] & DetX[c3]
          self.model.addConstr(uFrMC[r,c0] <= uDetX[r,c0], name = 'ufrMC0')
          self.model.addConstr(uFrMC[r,c0] <= udx[r,c3], name = 'ufrMC0')
          self.model.addConstr(uFrMC[r,c0] <= udw[r-1,c3], name = 'ufrMC0')
          self.model.addConstr(uFrMC[r,c0] <= 1 - udx[r,c0], name = 'ufrMC0')
          self.model.addConstr(uFrMC[r,c0] >= uDetX[r,i] + udx[r,c3] + udw[r-1,c3] + (1-udx[r,c0]) - 3, name = 'ufrMC0')
          # ------------------------------------------------------------------------------------------------------------
          self.model.addConstr(ufrMCa[r,c0] <= udetXa[r,c0], name = 'ufrmc0A')
          self.model.addConstr(ufrMCa[r,c0] <= udx[r,c3], name = 'ufrmc0A')
          self.model.addConstr(ufrMCa[r,c0] <= udw[r-1,c3], name = 'ufrmc0A')
          self.model.addConstr(ufrMCa[r,c0] <= 1 - udx[r,c0], name = 'ufrmc0A')
          self.model.addConstr(ufrMCa[r,c0] >= udetXa[r,i] + udx[r,c3] + udw[r-1,c3] + (1-udx[r,c0]) - 3, name = 'ufrmc0A')

          # POSITION -- 1: DetX[c1] = DetW[c0] (No filter in this position)
          self.model.addConstr(uFrMC[r,c1] == 0, 'ufrMC1')
          # ----------------------------------------------
          self.model.addConstr(ufrMCa[r,c1] == 0, 'ufrmc1A')

          # POSITION -- 2: DetX[c2] = DetW[c1] & DetW[c2]
          self.model.addConstr(uFrMC[r,c2] <= uDetX[r,c2], name = 'ufrMC2')
          self.model.addConstr(uFrMC[r,c2] <= udw[r-1,c1], name = 'ufrMC2')
          self.model.addConstr(uFrMC[r,c2] <= udw[r-1,c2], name = 'ufrMC2')
          self.model.addConstr(uFrMC[r,c2] <= 1 - udx[r,c2], name = 'ufrMC2')
          self.model.addConstr(uFrMC[r,c2] >= uDetX[r,c2] + udw[r-1,c1] + udw[r-1,c2] + (1-udx[r,c2]) - 3, name = 'ufrMC2')
          # ------------------------------------------------------------------------------------------------------------
          self.model.addConstr(ufrMCa[r,c2] <= udetXa[r,c2], name = 'ufrmc2A')
          self.model.addConstr(ufrMCa[r,c2] <= udw[r-1,c1], name = 'ufrmc2A')
          self.model.addConstr(ufrMCa[r,c2] <= udw[r-1,c2], name = 'ufrmc2A')
          self.model.addConstr(ufrMCa[r,c2] <= 1 - udx[r,c2], name = 'ufrmc2A')
          self.model.addConstr(ufrMCa[r,c2] >= udetXa[r,c2] + udw[r-1,c1] + udw[r-1,c2] + (1-udx[r,c2]) - 3, name = 'ufrmc2A')

          # POSITION -- 3: DetX[c3] = DetW[c0] & DetW[c2]
          self.model.addConstr(uFrMC[r,c3] <= uDetX[r,c3], name = 'ufrMC3')
          self.model.addConstr(uFrMC[r,c3] <= udw[r-1,c0], name = 'ufrMC3')
          self.model.addConstr(uFrMC[r,c3] <= udw[r-1,c2], name = 'ufrMC3')
          self.model.addConstr(uFrMC[r,c3] <= 1 - udx[r,c3], name = 'ufrMC3')
          self.model.addConstr(uFrMC[r,c3] >= uDetX[r,c3] + udw[r-1,c0] + udw[r-1,c2] + (1-udx[r,c2]) - 3, name = 'ufrMC3')
          # ------------------------------------------------------------------------------------------------------------
          self.model.addConstr(ufrMCa[r,c3] <= udetXa[r,c3], name = 'ufrmc3A')
          self.model.addConstr(ufrMCa[r,c3] <= udw[r-1,c0], name = 'ufrmc3A')
          self.model.addConstr(ufrMCa[r,c3] <= udw[r-1,c2], name = 'ufrmc3A')
          self.model.addConstr(ufrMCa[r,c3] <= 1 - udx[r,c3], name = 'ufrmc3A')
          self.model.addConstr(ufrMCa[r,c3] >= udetXa[r,c3] + udw[r-1,c0] + udw[r-1,c2] + (1-udx[r,c2]) - 3, name = 'ufrmc3A')

    # Setting pre-guess all involved key in auxiliary upper trail (rf = friter)
    # -------------------------------------------------------------------------------------
    self.model.addConstr(sum(udw[0,i] for i in range(16)) == 
                         (sum(ufrSCa[r,i] for r in range(1, self.rEb+1) for i in range(16)) + 
                          sum(ufrMCa[r,i] for r in range(2, self.rEb+1) for i in range(16))),
                          name = 'rfUfr')
    # -------------------------------------------------------------------------------------

    ''' 
    ============================================
            Guess-and-Determine (l)
    ============================================
    ''' 
    leR = self.rEb + self.rDis
    totR = self.rEb + self.rDis + self.rEf
    # NOTE: Assert rEb <= ri
    lGstk = self.model.addVars(range(leR, totR), 16, vtype = GRB.BINARY, name = 'lGstk')
    for r in range(leR, totR): 
      for i in range(8,16):
        self.model.addConstr(lGstk[r,i] == 0)
    lDetX = self.model.addVars(range(leR, totR),                 16, vtype = GRB.BINARY, name = 'lDetX')
    lDetY = self.model.addVars(range(leR, totR),                 16, vtype = GRB.BINARY, name = 'lDetY')
    lDetZ = self.model.addVars(range(leR, totR),                 16, vtype = GRB.BINARY, name = 'lDetZ') # Ciphertext
    lDetW = self.model.addVars(range(leR, totR-1),               16, vtype = GRB.BINARY, name = 'lDetW')

    lFrSC = self.model.addVars(range(leR, totR),   16, vtype = GRB.BINARY, name = 'lFrSC') # According to X
    lFrMC = self.model.addVars(range(leR, totR-1), 16, vtype = GRB.BINARY, name = 'lFrMC') # According to W
    # ----------------------------------------------------------------------------------------------------
    lgstka = self.model.addVars(range(leR, totR), 16, vtype = GRB.BINARY, name = 'lgstka')
    for r in range(leR, totR): 
      for i in range(8,16):
        self.model.addConstr(lgstka[r,i] == 0)
    ldetXa = self.model.addVars(range(leR, totR),                 16, vtype = GRB.BINARY, name = 'ldetXa')
    ldetYa = self.model.addVars(range(leR, totR),                 16, vtype = GRB.BINARY, name = 'ldetYa')
    ldetZa = self.model.addVars(range(leR, totR),                 16, vtype = GRB.BINARY, name = 'ldetZa') # Ciphertext
    ldetWa = self.model.addVars(range(leR, totR-1),               16, vtype = GRB.BINARY, name = 'ldetWa')

    lfrSCa = self.model.addVars(range(leR, totR),   16, vtype = GRB.BINARY, name = 'lfrSCa') # According to X
    lfrMCa = self.model.addVars(range(leR, totR-1), 16, vtype = GRB.BINARY, name = 'lfrMCa') # According to W

    self.model.addConstr(sum(lDetZ[totR-1,i] for i in range(16)) == 16, name = 'VC')

    ''' KR: Determine and Filter in SC '''
    for r in range(leR, totR):
      for i in range(16):
        # 'Det.' propagation via SC
        self.model.addConstr(lDetX[r,i] == lDetY[r,i], name = 'ldetSC')
        # ----------------------------------------------------------------
        self.model.addConstr(ldetXa[r,i] == ldetYa[r,i], name = 'ldetSCa')

        # Filter obtain from SC
        self.model.addConstr(lDX[r,i] - ldx[r,i] - lFrSC[r,i] >= 0, name = 'lfrSC')
        self.model.addConstr(lDetX[r,i] >= lFrSC[r,i], name = 'lfrSC')
        self.model.addConstr(- lDetX[r,i] - lDX[r,i] + ldx[r,i] + lFrSC[r,i] >= -1, name = 'lfrSC')
        # -----------------------------------------------------------------------------------------
        self.model.addConstr(lDX[r,i] - ldx[r,i] - lfrSCa[r,i] >= 0, name = 'lfrSCa')
        self.model.addConstr(ldetXa[r,i] >= lfrSCa[r,i], name = 'lfrSCa')
        self.model.addConstr(- ldetXa[r,i] - lDX[r,i] + ldx[r,i] + lfrSCa[r,i] >= -1, name = 'lfrSCa')

    ''' KR: Guess round key and determine'''
    for r in range(leR, totR):
      # rk = r + self.r0
      for i in range(8):
        self.model.addConstr(lDetY[r,i] <= lGstk[r,i], name = 'lgkiU')
        self.model.addConstr(lDetY[r,i] <= lDetZ[r,i], name = 'lgkiU')
        self.model.addConstr( - lGstk[r,i] - lDetZ[r,i] + lDetY[r,i] + 1 >= 0, name = 'lgkiU')
        # -------------------------------------------------------------------------------------
        self.model.addConstr(ldetYa[r,i] <= lgstka[r,i], name = 'lgkiUa')
        self.model.addConstr(ldetYa[r,i] <= ldetZa[r,i], name = 'lgkiUa')
        self.model.addConstr( - lgstka[r,i] - ldetZa[r,i] + ldetYa[r,i] + 1 >= 0, name = 'lgkiUa')
      for i in range(8,16):
        self.model.addConstr(lDetY[r,i] == lDetZ[r,i], name = 'lgkiL')
        # ---------------------------------------------------------------
        self.model.addConstr(ldetYa[r,i] == ldetZa[r,i], name = 'lgkiLa')

    ''' KR: 'Det.' propagation in SR '''
    for r in range(leR, totR-1):
      for i in range(16): 
        self.model.addConstr(lDetZ[r,i] == lDetW[r,self.SRp[i]], name = 'ldetSR')
        self.model.addConstr(ldetZa[r,i] == ldetWa[r,self.SRp[i]], name = 'ldetSRa')

    ''' KR: 'Det.' propagation in MC '''
    for r in range(leR+1, totR):
      for c in [[0,4,8,12],[1,5,9,13],[2,6,10,14],[3,7,11,15]]:
          c0, c1, c2, c3 = c[0], c[1], c[2], c[3]

          # POSITION -- 0: DetW[c0] = DetX[c1]
          self.model.addConstr(lDetW[r-1,c0] == lDetX[r,c1], name = 'ldetMC0')

          # POSITION -- 1: DetW[c1] = DetW[c2] + DetX[c2]
          self.model.addConstr(lDetW[r-1,c1] <= lDetW[r-1,c2], name = 'ldetMC1')
          self.model.addConstr(lDetW[r-1,c1] <= lDetX[r,c2], name = 'ldetMC1')
          self.model.addConstr( - lDetW[r-1,c2] - lDetX[r,c2] + lDetW[r-1,c1] + 1 >= 0, name = 'ldetMC1')

          # POSITION -- 2: DetW[c2] = DetX[c1] + DetX[c3]
          self.model.addConstr(lDetW[r-1,c2] <= lDetX[r,c1], name = 'ldetMC2')
          self.model.addConstr(lDetW[r-1,c2] <= lDetX[r,c3], name = 'ldetMC2')
          self.model.addConstr( - lDetX[r,c1] - lDetX[r,c3] + lDetW[r-1,c2] + 1 >= 0, name = 'ldetMC2')

          # POSITION -- 3: DetX[c3] = DetW[c0] & DetW[c3]
          self.model.addConstr(lDetW[r-1,c3] <= lDetX[r,c0], name = 'ldetMC3')
          self.model.addConstr(lDetW[r-1,c3] <= lDetX[r,c3], name = 'ldetMC3')
          self.model.addConstr( - lDetX[r,c0] - lDetX[r,c3] + lDetW[r-1,c3] + 1 >= 0, name = 'ldetMC3')
          # -------------------------------------------------------------------------------------------
          # POSITION -- 0: DetW[c0] = DetX[c1]
          self.model.addConstr(ldetWa[r-1,c0] == ldetXa[r,c1], name = 'ldetMC0a')

          # POSITION -- 1: DetW[c1] = DetW[c2] + DetX[c2]
          self.model.addConstr(ldetWa[r-1,c1] <= ldetWa[r-1,c2], name = 'ldetMC1a')
          self.model.addConstr(ldetWa[r-1,c1] <= ldetXa[r,c2], name = 'ldetMC1a')
          self.model.addConstr( - ldetWa[r-1,c2] - ldetXa[r,c2] + ldetWa[r-1,c1] + 1 >= 0, name = 'ldetMC1a')

          # POSITION -- 2: DetW[c2] = DetX[c1] + DetX[c3]
          self.model.addConstr(ldetWa[r-1,c2] <= ldetXa[r,c1], name = 'ldetMC2a')
          self.model.addConstr(ldetWa[r-1,c2] <= ldetXa[r,c3], name = 'ldetMC2a')
          self.model.addConstr( - ldetXa[r,c1] - ldetXa[r,c3] + ldetWa[r-1,c2] + 1 >= 0, name = 'ldetMC2a')

          # POSITION -- 3: DetX[c3] = DetW[c0] & DetW[c3]
          self.model.addConstr(ldetWa[r-1,c3] <= ldetXa[r,c0], name = 'ldetMC3a')
          self.model.addConstr(ldetWa[r-1,c3] <= ldetXa[r,c3], name = 'ldetMC3a')
          self.model.addConstr( - ldetXa[r,c0] - ldetXa[r,c3] + ldetWa[r-1,c3] + 1 >= 0, name = 'ldetMC3a')

    ''' KR: 'Filter' obtain from MC '''
    for r in range(leR, totR-1):
      for c in [[0,4,8,12],[1,5,9,13],[2,6,10,14],[3,7,11,15]]:
          c0, c1, c2, c3 = c[0], c[1], c[2], c[3]

          # POSITION -- 0: DetW[c0] = DetX[c1] (No filter in this position)
          self.model.addConstr(lFrMC[r,c0] == 0, name = 'lfrMC0')
          
          # POSITION -- 1: DetW[c1] = DetW[c2] + DetX[c2]
          self.model.addConstr(lFrMC[r,c1] <= lDetW[r,c1], name = 'lfrMC1')
          self.model.addConstr(lFrMC[r,c1] <= ldw[r,c2], name = 'lfrMC1')
          self.model.addConstr(lFrMC[r,c1] <= ldx[r+1,c2], name = 'lfrMC1')
          self.model.addConstr(lFrMC[r,c1] <= 1 - ldw[r,c1], name = 'lfrMC1')
          self.model.addConstr(lFrMC[r,c1] >= lDetW[r,c1] + ldw[r,c2] + ldx[r+1,c2] + (1-ldw[r,c1]) - 3, name = 'lfrMC1')

          # POSITION -- 2: DetW[c2] = DetX[c1] + DetX[c3]
          self.model.addConstr(lFrMC[r,c2] <= lDetW[r,c2], name = 'lfrMC2')
          self.model.addConstr(lFrMC[r,c2] <= ldx[r+1,c1], name = 'lfrMC2')
          self.model.addConstr(lFrMC[r,c2] <= ldx[r+1,c3], name = 'lfrMC2')
          self.model.addConstr(lFrMC[r,c2] <= 1 - ldw[r,c2], name = 'lfrMC2')
          self.model.addConstr(lFrMC[r,c2] >= lDetW[r,c2] + ldx[r+1,c1] + ldx[r+1,c3] + (1-ldw[r,c2]) - 3, name = 'lfrMC2')

          # POSITION -- 3: DetX[c3] = DetW[c0] & DetW[c3]
          self.model.addConstr(lFrMC[r,c3] <= lDetW[r,c3], name = 'lfrMC3')
          self.model.addConstr(lFrMC[r,c3] <= ldx[r+1,c0], name = 'lfrMC3')
          self.model.addConstr(lFrMC[r,c3] <= ldx[r+1,c3], name = 'lfrMC3')
          self.model.addConstr(lFrMC[r,c3] <= 1 - ldw[r,c3], name = 'lfrMC3')
          self.model.addConstr(lFrMC[r,c3] >= lDetW[r,c3] + ldx[r+1,c0] + ldx[r+1,c3] + (1-ldw[r,c3]) - 3, name = 'lfrMC3')

          # ---------------------------------------------------------------------------------------------------------------
          # POSITION -- 0: DetW[c0] = DetX[c1] (No filter in this position)
          self.model.addConstr(lfrMCa[r,c0] == 0, name = 'lfrMC0a')
          
          # POSITION -- 1: DetW[c1] = DetW[c2] + DetX[c2]
          self.model.addConstr(lfrMCa[r,c1] <= ldetWa[r,c1], name = 'lfrMC1a')
          self.model.addConstr(lfrMCa[r,c1] <= ldw[r,c2], name = 'lfrMC1a')
          self.model.addConstr(lfrMCa[r,c1] <= ldx[r+1,c2], name = 'lfrMC1a')
          self.model.addConstr(lfrMCa[r,c1] <= 1 - ldw[r,c1], name = 'lfrMC1a')
          self.model.addConstr(lfrMCa[r,c1] >= ldetWa[r,c1] + ldw[r,c2] + ldx[r+1,c2] + (1-ldw[r,c1]) - 3, name = 'lfrMC1a')

          # POSITION -- 2: DetW[c2] = DetX[c1] + DetX[c3]
          self.model.addConstr(lfrMCa[r,c2] <= ldetWa[r,c2], name = 'lfrMC2a')
          self.model.addConstr(lfrMCa[r,c2] <= ldx[r+1,c1], name = 'lfrMC2a')
          self.model.addConstr(lfrMCa[r,c2] <= ldx[r+1,c3], name = 'lfrMC2a')
          self.model.addConstr(lfrMCa[r,c2] <= 1 - ldw[r,c2], name = 'lfrMC2a')
          self.model.addConstr(lfrMCa[r,c2] >= ldetWa[r,c2] + ldx[r+1,c1] + ldx[r+1,c3] + (1-ldw[r,c2]) - 3, name = 'lfrMC2a')

          # POSITION -- 3: DetX[c3] = DetW[c0] & DetW[c3]
          self.model.addConstr(lfrMCa[r,c3] <= ldetWa[r,c3], name = 'lfrMC3a')
          self.model.addConstr(lfrMCa[r,c3] <= ldx[r+1,c0], name = 'lfrMC3a')
          self.model.addConstr(lfrMCa[r,c3] <= ldx[r+1,c3], name = 'lfrMC3a')
          self.model.addConstr(lfrMCa[r,c3] <= 1 - ldw[r,c3], name = 'lfrMC3a')
          self.model.addConstr(lfrMCa[r,c3] >= ldetWa[r,c3] + ldx[r+1,c0] + ldx[r+1,c3] + (1-ldw[r,c3]) - 3, name = 'lfrMC3a')

    # Setting pre-guess all involved key in auxiliary lower trail (rf = friter)
    # -------------------------------------------------------------------------------------
    self.model.addConstr(sum(ldz[totR-1,i] for i in range(16)) == 
                         (sum(lfrSCa[r,i] for r in range(leR, totR) for i in range(16)) + 
                          sum(lfrMCa[r,i] for r in range(leR, totR-1) for i in range(16))),
                          name = 'rfLfr')
    # -------------------------------------------------------------------------------------


    ''' 
    ============================
            Key Bridge
    ============================
    '''
    '''
    For involved key
    '''
    iLANE = self.model.addVars(16, vtype = GRB.INTEGER, name = 'iLANE') # Count the number of stk in a LANE
    ikSUM = self.model.addVars(16, ub = self.Vs, vtype = GRB.INTEGER, name = 'ikSUM') # Judge the broundary Vs
    ikz = self.model.addVars(16, vtype = GRB.BINARY, name = 'ikz') # Auxiliary vars. for ikSUM

    for i in range(16):
      self.model.addConstr(iLANE[i] == sum(ugstka[r,self.iterate_hTable(i,r)] for r in range(self.rEb)) + 
                                       sum(lgstka[r,self.iterate_hTable(i,r)] for r in range(leR, totR)), 
                                       name = 'KeyBgeI')
    for i in range(16):
      self.model.addConstr(iLANE[i] <= self.Vs + 1000*ikz[i], name = 'KeyBgeI')
      self.model.addConstr(iLANE[i] >= (self.Vs+1) - 1000*(1-ikz[i]), name = 'KeyBgeI')

      self.model.addConstr(ikSUM[i] <= iLANE[i] + 1000*ikz[i], name = 'KeyBgeI')
      self.model.addConstr(ikSUM[i] >= iLANE[i] - 1000*ikz[i], name = 'KeyBgeI')
      self.model.addConstr(ikSUM[i] <= self.Vs + 1000*(1-ikz[i]), name = 'KeyBgeI')
      self.model.addConstr(ikSUM[i] >= self.Vs - 1000*(1-ikz[i]), name = 'KeyBgeI')

    '''
    For pre-guessed key
    '''
    gLANE = self.model.addVars(16, vtype = GRB.INTEGER, name = 'gLANE') # Count the number of Gstk in a LANE
    gkSUM = self.model.addVars(16, ub = 2, vtype = GRB.INTEGER, name = 'gkSUM') # Judge the broundary Vs
    gkz = self.model.addVars(16, vtype = GRB.BINARY, name = 'gkz') # Auxiliary vars. for gkSUM

    for i in range(16):
      self.model.addConstr(gLANE[i] == sum(uGstk[r,self.iterate_hTable(i,r)] for r in range(self.rEb)) + 
                                       sum(lGstk[r,self.iterate_hTable(i,r)] for r in range(leR, totR)), 
                                       name = 'KeyBgeG')
    for i in range(16):
      self.model.addConstr(gLANE[i] <= self.Vs + 1000*gkz[i], name = 'KeyBgeG')
      self.model.addConstr(gLANE[i] >= (self.Vs+1) - 1000*(1-gkz[i]), name = 'KeyBgeG')

      self.model.addConstr(gkSUM[i] <= gLANE[i] + 1000*gkz[i], name = 'KeyBgeG')
      self.model.addConstr(gkSUM[i] >= gLANE[i] - 1000*gkz[i], name = 'KeyBgeG')
      self.model.addConstr(gkSUM[i] <= self.Vs + 1000*(1-gkz[i]), name = 'KeyBgeG')
      self.model.addConstr(gkSUM[i] >= self.Vs - 1000*(1-gkz[i]), name = 'KeyBgeG')
   

    '''
    ==========================================================================================
         (Begin) - Complexity
    '''
    '''
    Parameters
    '''
    rb  = self.model.addVar(vtype = GRB.INTEGER, name = 'rb')
    rf  = self.model.addVar(vtype = GRB.INTEGER, name = 'rf')
    cb  = self.model.addVar(vtype = GRB.INTEGER, name = 'cb')
    cf  = self.model.addVar(vtype = GRB.INTEGER, name = 'cf')
    iks = self.model.addVar(vtype = GRB.INTEGER, name = 'iks')
    gks = self.model.addVar(vtype = GRB.INTEGER, name = 'gks')
    cbp = self.model.addVar(vtype = GRB.INTEGER, name = 'cbp')
    cfp = self.model.addVar(vtype = GRB.INTEGER, name = 'cfp')

    # mb  + mf  = ikSUM (in key bridge)
    self.model.addConstr(iks == self.c_size * sum(ikSUM[i] for i in range(16)), name = 'IKS')
    # mbp + mfp = gkSUM (in key bridge)
    self.model.addConstr(gks == self.c_size * sum(gkSUM[i] for i in range(16)), name = 'GKS')

    self.model.addConstr(rb == self.c_size * sum(udw[0,i] for i in range(16)), name = 'rb')
    self.model.addConstr(rf == self.c_size * sum(ldz[totR-1,i] for i in range(16)), name = 'rf')

    self.model.addConstr(cb == rb - self.c_size * sum(udy[self.rEb,i] for i in range(16)), name = 'cb')
    self.model.addConstr(cf == rf - self.c_size * sum(ldx[totR-self.rEf,i] for i in range(16)), name = 'cf')

    self.model.addConstr(cbp == self.c_size * 
                         (sum(uFrSC[r,i] for r in range(1, self.rEb+1) for i in range(16)) + 
                          sum(uFrMC[r,i] for r in range(2, self.rEb+1) for i in range(16))), name = 'cbp')
    self.model.addConstr(cfp == self.c_size * 
                         (sum(lFrSC[r,i] for r in range(leR, totR) for i in range(16)) + 
                          sum(lFrMC[r,i] for r in range(leR, totR-1) for i in range(16))), name = 'cfp')

    '''
    Complexities
    '''
    D   = self.model.addVar(lb = 0, vtype = GRB.CONTINUOUS, name = 'D')
    Dc  = self.model.addVar(lb = 0, vtype = GRB.CONTINUOUS, name = 'Dc')
    Qc  = self.model.addVar(lb = 0, vtype = GRB.CONTINUOUS, name = 'Qc')

    Mc  = self.model.addVar(lb = 0, vtype = GRB.CONTINUOUS, name = 'Mc')

    T0  = self.model.addVar(lb = 0, vtype = GRB.CONTINUOUS, name = 'T0')
    T1  = self.model.addVar(lb = 0, vtype = GRB.CONTINUOUS, name = 'T1')
    T2  = self.model.addVar(lb = 0, vtype = GRB.CONTINUOUS, name = 'T2')
    T31 = self.model.addVar(lb = 0, vtype = GRB.CONTINUOUS, name = 'T31')
    T32 = self.model.addVar(lb = 0, vtype = GRB.CONTINUOUS, name = 'T32')
    T4  = self.model.addVar(lb = 0, vtype = GRB.CONTINUOUS, name = 'T4')
    Tc  = self.model.addVar(lb = 0, vtype = GRB.CONTINUOUS, name = 'Tc')

    # Data and Quartet
    self.model.addConstr(D == self.b_size + self.Vz/2, name = 'D')
    self.model.addConstr(Dc == D + 2, name = 'Dc')
    self.model.addConstr(Qc == 2 * (cb + cf) + self.Vz, name = 'Qc')

    # Memory
    self.model.addConstr(Mc >= Dc + 2, name = 'Mc0')
    self.model.addConstr(Mc >= T2 - gks, name = 'Mc1')
    self.model.addConstr(Mc >= Qc - 2*(cbp + cfp), name = 'Mc2')
    self.model.addConstr(Mc >= iks, name = 'Mc3')

    # Time
    # T0 (T0 = Dc)
    self.model.addConstr(T0 == Dc, name = 'T0')
    # T1 (T1 = 2^{mb'+mf'} · D · 4)
    self.model.addConstr(T1 == gks + Dc + 2, name = 'T1')
    # T2 (T2 = 2^{mb'+mf'} · D · min[ 2^{rb-cb'}, D · 2^{rf-cf'-n} ] · 2)
    if self.cP == True:
      self.model.addConstr(T2 == gks + Dc + rb - cbp + 1, name = 'T2p')
    else:
      self.model.addConstr(T2 == gks + Dc + (Dc + rf - cfp - self.b_size) + 1, name = 'T2c')
    # T31 (T31 = 2^{mb'+mf'-2cb'-2cf'} · Q)
    self.model.addConstr(T31 == (gks - 2*cbp - 2*cfp) + Qc, name = 'T31')
    # T32 (T32 = 2^{mb+mf+z} · epsilon)
    self.model.addConstr(T32 == (iks + self.Vz), name = 'T32')
    # T4 (T4 = 2^{k-x}, where 2^z = xln2)
    self.model.addConstr(T4 == self.k_size - self.Vx, name = 'T4')
    # Tc is the maximum Ti
    self.model.addConstr(Tc >= T0, name = 'Tc')
    self.model.addConstr(Tc >= T1, name = 'Tc')
    self.model.addConstr(Tc >= T2, name = 'Tc')
    self.model.addConstr(Tc >= T31, name = 'Tc')
    self.model.addConstr(Tc >= T32, name = 'Tc')
    self.model.addConstr(Tc >= T4, name = 'Tc')


    # ******************************************************
    ''' Contradiction '''
    self.model.addConstr(sum((uDX[r,i] - udx[r,i]) * (lDY[r,i] - ldy[r,i]) 
                             for i in range(16) for r in range(self.rEb, self.rEb+self.rDis)) >= 1, name = 'iBCT')
    # ******************************************************


    ''' Objective function'''
    self.model.ModelSense = GRB.MINIMIZE
    self.model.setObjectiveN(Tc, index=0, priority=4, name='Min_Tc')
    self.model.setObjectiveN(T32, index=1, priority=3, name='Min_T32')
    self.model.setObjectiveN(Dc, index=2, priority=2, name='Min_Dc')
    self.model.setObjectiveN(Mc, index=3, priority=1, name='Min_Mc')

    # self.model.setObjectiveN(T1, index=4, priority=0, name='Min_T1')
    # self.model.setObjectiveN(T2, index=5, priority=0, name='Min_T2')
    # self.model.setObjectiveN(T31, index=6, priority=0, name='Min_T31')
    
    # -------------------- To count mb and mf ------------------------
    mbi  = self.model.addVar(lb = 0, vtype = GRB.INTEGER, name = 'mbi')
    self.model.addConstr(mbi == sum(ugstka[r,i] for r in range(rEb) for i in range(8)), name = 'mb')
    mfi  = self.model.addVar(lb = 0, vtype = GRB.INTEGER, name = 'mfi')
    self.model.addConstr(mfi == sum(lgstka[r,i] for r in range(leR, totR) for i in range(8)), name = 'mf')
    
    self.model.setObjectiveN(mbi, index=10, priority=0, name='mbi')
    self.model.setObjectiveN(mfi, index=11, priority=0, name='mfi')
    # -----------------------------------------------------------------




    # self.model.setParam("OutputFlag", 0)
    self.model.write(self.name + '.lp')
    self.model.optimize()

    # if self.model.Status == GRB.INFEASIBLE:
    #     print('-'*40 + '\n| Model is infeasible, computing IIS...|\n' + '-'*40)
    #     self.model.computeIIS()
    #     self.model.write(self.name + '.ilp')
    self.model.write(self.name + '.sol')

    # if self.model.Status == 2:
    #   return Tc.x
    # return self.model.Status


  '''
  =====================
  Tools:
  =====================
  '''

  ''' x to z (2^z = xln2) '''
  def x_to_z(self,x):
    # ln2 = 0.69...
    return round(math.log2(0.7 * x), 1)

  ''' n iterations hPermutation '''
  def iterate_hTable(self, pGstk, n):
    current = pGstk
    for _ in range(n):
        current = self.hTable[current]
    return current
  



if __name__ == '__main__':

  rEb = 6
  rEf = 5
  rDis = 22
  SKINNYe = IB_ForkSKINNY(64, 4, rEb, rDis, rEf,  32, True)
  SKINNYe.ib_model()
