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
    self.name = './v4_SKINNYe-{}-{}_{}-{}+{}+{}_{}'.format(
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
    uDX0  = self.model.addVars(range(1, uR+1), 16, vtype = GRB.BINARY, name = 'uDX0')
    udx0  = self.model.addVars(range(1, uR+1), 16, vtype = GRB.BINARY, name = 'udx0')
    uDY0  = self.model.addVars(range(1, uR),   16, vtype = GRB.BINARY, name = 'uDY0')
    udy0  = self.model.addVars(range(1, uR),   16, vtype = GRB.BINARY, name = 'udy0')
    uDZ0  = self.model.addVars(range(1, uR),   16, vtype = GRB.BINARY, name = 'uDZ0') 
    udz0  = self.model.addVars(range(1, uR),   16, vtype = GRB.BINARY, name = 'udz0') 
    uDW0  = self.model.addVars(uR,             16, vtype = GRB.BINARY, name = 'uDW0') # W0 = Plaintext
    udw0  = self.model.addVars(uR,             16, vtype = GRB.BINARY, name = 'udw0') # w0 = Plaintext
    ustk0 = self.model.addVars(uR,     16, vtype = GRB.BINARY, name = 'ustk0')

    '''Basic Constrs:'''
    # Remove dx[r,i] > DX[r,i]
    for i in range(16): 
      for r in range(1, uR+1):
        self.model.addConstr(uDX0[r, i] >= udx0[r, i], name = 'uBasic')
      for r in range(1, uR):
        self.model.addConstr(uDY0[r, i] >= udy0[r, i], name = 'uBasic')
        self.model.addConstr(uDZ0[r, i] >= udz0[r, i], name = 'uBasic')
      for r in range(uR):
        self.model.addConstr(uDW0[r, i] >= udw0[r, i], name = 'uBasic')

    '''Operation: SC'''
    # In Extension
    for r in range(1, self.rEb+1):
      for i in range(16):
        self.model.addConstr(uDX0[r,i] == uDY0[r,i], name='ueSC')
        self.model.addConstr(udx0[r,i] == uDY0[r,i], name='ueSC')
    # In distinguisher
    for r in range(self.rEb+1, uR):
      for i in range(16):
        self.model.addConstr(uDY0[r,i] == uDX0[r,i], name='udSC')
        self.model.addConstr(udy0[r,i] == uDX0[r,i], name='udSC')


    '''(*ri-(r0)-r1)Operation: ATK'''
  
    ''' 
    ----------------------------------
    Equivalent Key for the first round
    ----------------------------------
    '''
    eqk0 = self.model.addVars(16, vtype = GRB.BINARY, name = 'eqk_0')

    for c in [[0,4,8,12],[1,5,9,13],[2,6,10,14],[3,7,11,15]]:
        c0, c1, c2, c3 = c[0], c[1], c[2], c[3]
        tc0, tc1, tc2, tc3 = self.SRpv[c[0]], self.SRpv[c[1]], self.SRpv[c[2]], self.SRpv[c[3]]
        
        # --- POSITION -- 0: eqk[c0] = stk1[c3] + eqk[c3] ---
        self.model.addConstr(eqk0[c0] == ustk0[0,tc0], name = '2eqk0')

        # --- POSITION -- 1: eqk[c1] = stk1[c0] ---
        self.model.addConstr(eqk0[c1] == ustk0[0,tc0], name = '2eqk1')
        
        # --- POSITION -- 2: eqk[c2] = stk1[c1] + stk[c2] ---
        self.model.addConstr(eqk0[c2] == ustk0[0,tc1], name = '2eqk2')

        # --- POSITION -- 2: eqk[c3] = stk1[c0] + stk[c2] ---
        self.model.addConstr(eqk0[c3] == ustk0[0,tc0], name = '2eqk3')


    # Add eqk in the first round
    for i in range(16):
      self.model.addConstr(udw0[0,i] == udx0[1,i], name='uAEqk')
      self.model.addConstr(eqk0[i] + uDW0[0,i] - uDX0[1,i] >= 0, name='uAEqk')
      self.model.addConstr(eqk0[i] - uDW0[0,i] + uDX0[1,i] >= 0, name='uAEqk')
      self.model.addConstr(-eqk0[i] + uDW0[0,i] + uDX0[1,i] - udx0[1,i] >= 0, name='uAEqk')
      self.model.addConstr(uDW0[0,i] >= udx0[1,i], name='uAEqk')
    # # In ri
    # for r in range(1, self.ri):
    #   for i in range(8):
    #     self.model.addConstr(udy0[r,i] == udz0[r,i], name='uATKri')
    #     self.model.addConstr(ustk0[r,i] + uDY0[r,i] - uDZ0[r,i] >= 0, name='uATKri')
    #     self.model.addConstr(ustk0[r,i] - uDY0[r,i] + uDZ0[r,i] >= 0, name='uATKri')
    #     self.model.addConstr(-ustk0[r,i] + uDY0[r,i] + uDZ0[r,i] - udz0[r,i] >= 0, name='uATKri')
    #     self.model.addConstr(uDY0[r,i] >= udz0[r,i], name='uATKri')
    #   for i in range(8,16):
    #     self.model.addConstr(uDY0[r,i] == uDZ0[r,i], name='uATKri')
    #     self.model.addConstr(udy0[r,i] == udz0[r,i], name='uATKri')
    # # In r1
    # for r in range(self.ri, uR):
    #   tR = r + self.r0
    #   for i in range(8):
    #     self.model.addConstr(udy0[r,i] == udz0[r,i], name='uATKr1')
    #     self.model.addConstr(ustk0[tR,i] + uDY0[r,i] - uDZ0[r,i] >= 0, name='uATKr1')
    #     self.model.addConstr(ustk0[tR,i] - uDY0[r,i] + uDZ0[r,i] >= 0, name='uATKr1')
    #     self.model.addConstr(-ustk0[tR,i] + uDY0[r,i] + uDZ0[r,i] - udz0[r,i] >= 0, name='uATKr1')
    #     self.model.addConstr(uDY0[r,i] >= udz0[r,i], name='uATKr1')
    #   for i in range(8,16):
    #     self.model.addConstr(uDY0[r,i] == uDZ0[r,i], name='uATKr1')
    #     self.model.addConstr(udy0[r,i] == udz0[r,i], name='uATKr1')
    for r in range(1, uR):
      for i in range(8):
        self.model.addConstr(udy0[r,i] == udz0[r,i], name='uATK')
        self.model.addConstr(ustk0[r,i] + uDY0[r,i] - uDZ0[r,i] >= 0, name='uATK')
        self.model.addConstr(ustk0[r,i] - uDY0[r,i] + uDZ0[r,i] >= 0, name='uATK')
        self.model.addConstr(-ustk0[r,i] + uDY0[r,i] + uDZ0[r,i] - udz0[r,i] >= 0, name='uATK')
        self.model.addConstr(uDY0[r,i] >= udz0[r,i], name='uATK')
      for i in range(8,16):
        self.model.addConstr(uDY0[r,i] == uDZ0[r,i], name='uATK')
        self.model.addConstr(udy0[r,i] == udz0[r,i], name='uATK')


    '''(*)Operation: SR''' 
    for r in range(1,uR):
      for i in range(16):
        self.model.addConstr(uDZ0[r,i] == uDW0[r,self.SRp[i]], name = 'uSR')
        self.model.addConstr(udz0[r,i] == udw0[r,self.SRp[i]], name = 'uSR')

    '''Operation: MC'''
    # ************ In extension (backword) ************
    for r in range(1, self.rEb):
      for c in [[0,4,8,12],[1,5,9,13],[2,6,10,14],[3,7,11,15]]:
        c0, c1, c2, c3 = c[0], c[1], c[2], c[3]

        # --- POSITION -- 0: W[c0] = X[c1] ---
        self.model.addConstr(udw0[r,c0] == udx0[r+1,c1], name = 'ueMC0')

        self.model.addConstr(uDW0[r,c0] == uDX0[r+1,c1], name = 'ueMC0')

        # --- POSITION -- 1: W[c1] = W[c2] + X[c2] ---
        self.model.addConstr(udw0[r,c1] >= udw0[r,c2], name = 'ueMC1')
        self.model.addConstr(udw0[r,c1] >= udx0[r+1,c2], name = 'ueMC1')
        self.model.addConstr(udw0[r,c1] <= udw0[r,c2] + udx0[r+1,c2], name = 'ueMC1')

        self.model.addConstr(uDW0[r,c1] + uDW0[r,c2] + uDX0[r+1,c2] >= 2*uDW0[r,c1], name = 'ueMC1')
        self.model.addConstr(uDW0[r,c1] + uDW0[r,c2] + uDX0[r+1,c2] >= 2*uDW0[r,c2], name = 'ueMC1')
        self.model.addConstr(uDW0[r,c1] + uDW0[r,c2] + uDX0[r+1,c2] >= 2*uDX0[r+1,c2], name = 'ueMC1')
        self.model.addConstr(uDW0[r,c1] <= uDW0[r,c2] + uDX0[r+1,c2], name = 'ueMC1')

        # --- POSITION -- 2: W[c2] = X[c1] + X[c3] ---
        self.model.addConstr(udw0[r,c2] >= udx0[r+1,c1], name = 'ueMC2')
        self.model.addConstr(udw0[r,c2] >= udx0[r+1,c3], name = 'ueMC2')
        self.model.addConstr(udw0[r,c2] <= udx0[r+1,c1] + udx0[r+1,c3], name = 'ueMC2')

        self.model.addConstr(uDW0[r,c2] + uDX0[r+1,c1] + uDX0[r+1,c3] >= 2*uDW0[r,c2], name = 'ueMC2')
        self.model.addConstr(uDW0[r,c2] + uDX0[r+1,c1] + uDX0[r+1,c3] >= 2*uDX0[r+1,c1], name = 'ueMC2')
        self.model.addConstr(uDW0[r,c2] + uDX0[r+1,c1] + uDX0[r+1,c3] >= 2*uDX0[r+1,c3], name = 'ueMC2')
        self.model.addConstr(uDW0[r,c2] <= uDX0[r+1,c1] + uDX0[r+1,c3], name = 'ueMC2')

        # --- POSITION -- 3: W[c3] = X[c0] + X[c3] ---
        self.model.addConstr(udw0[r,c3] >= udx0[r+1,c0], name = 'ueMC3')
        self.model.addConstr(udw0[r,c3] >= udx0[r+1,c3], name = 'ueMC3')
        self.model.addConstr(udw0[r,c3] <= udx0[r+1,c0] + udx0[r+1,c3], name = 'ueMC3')

        self.model.addConstr(uDW0[r,c3] + uDX0[r+1,c0] + uDX0[r+1,c3] >= 2*uDW0[r,c3], name = 'ueMC3')
        self.model.addConstr(uDW0[r,c3] + uDX0[r+1,c0] + uDX0[r+1,c3] >= 2*uDX0[r+1,c0], name = 'ueMC3')
        self.model.addConstr(uDW0[r,c3] + uDX0[r+1,c0] + uDX0[r+1,c3] >= 2*uDX0[r+1,c3], name = 'ueMC3')
        self.model.addConstr(uDW0[r,c3] <= uDX0[r+1,c0] + uDX0[r+1,c3], name = 'ueMC3')
    
    # ************ In distinguisher (forward) ************
    for r in range(self.rEb, uR):
      for c in [[0,4,8,12],[1,5,9,13],[2,6,10,14],[3,7,11,15]]:
          c0, c1, c2, c3 = c[0], c[1], c[2], c[3]

          # --- POSITION -- 0: X[c0] = X[c3] + W[c3] ---
          self.model.addConstr(udx0[r+1,c0] >= udw0[r,c3], name = 'udMC0')
          self.model.addConstr(udx0[r+1,c0] >= udx0[r+1,c3], name = 'udMC0')
          self.model.addConstr(udx0[r+1,c0] <= udw0[r,c3] + udx0[r+1,c3], name = 'udMC0')

          self.model.addConstr(uDX0[r+1,c0] + uDW0[r,c3] + uDX0[r+1,c3] >= 2*uDX0[r+1,c0], name = 'udMC0')
          self.model.addConstr(uDX0[r+1,c0] + uDW0[r,c3] + uDX0[r+1,c3] >= 2*uDW0[r,c3], name = 'udMC0')
          self.model.addConstr(uDX0[r+1,c0] + uDW0[r,c3] + uDX0[r+1,c3] >= 2*uDX0[r+1,c3], name = 'udMC0')
          self.model.addConstr(uDX0[r+1,c0] <= uDW0[r,c3] + uDX0[r+1,c3], name = 'udMC0')

          # --- POSITION -- 1: X[c1] = W[c0] ---
          self.model.addConstr(udx0[r+1,c1] == udw0[r,c0], name = 'udMC1')

          self.model.addConstr(uDX0[r+1,c1] == uDW0[r,c0], name = 'udMC1')

          # --- POSITION -- 2: X[c2] = X[c1] + W[c2] ---
          self.model.addConstr(udx0[r+1,c2] >= udw0[r,c1], name = 'udMC2')
          self.model.addConstr(udx0[r+1,c2] >= udw0[r,c2], name = 'udMC2')
          self.model.addConstr(udx0[r+1,c2] <= udw0[r,c1] + udw0[r,c2], name = 'udMC2')

          self.model.addConstr(uDX0[r+1,c2] + uDW0[r,c1] + uDW0[r,c2] >= 2*uDX0[r+1,c2], name = 'udMC2')
          self.model.addConstr(uDX0[r+1,c2] + uDW0[r,c1] + uDW0[r,c2] >= 2*uDW0[r,c1], name = 'udMC2')
          self.model.addConstr(uDX0[r+1,c2] + uDW0[r,c1] + uDW0[r,c2] >= 2*uDW0[r,c2], name = 'udMC2')
          self.model.addConstr(uDX0[r+1,c2] <= uDW0[r,c1] + uDW0[r,c2], name = 'udMC2')

          # --- POSITION -- 3: X[c3] = W[c0] + W[c2] ---
          self.model.addConstr(udx0[r+1,c3] >= udw0[r,c0], name = 'udMC3')
          self.model.addConstr(udx0[r+1,c3] >= udw0[r,c2], name = 'udMC3')
          self.model.addConstr(udx0[r+1,c3] <= udw0[r,c0] + udw0[r,c2], name = 'udMC3')

          self.model.addConstr(uDX0[r+1,c3] + uDW0[r,c0] + uDW0[r,c2] >= 2*uDX0[r+1,c3], name = 'udMC3')
          self.model.addConstr(uDX0[r+1,c3] + uDW0[r,c0] + uDW0[r,c2] >= 2*uDW0[r,c0], name = 'udMC3')
          self.model.addConstr(uDX0[r+1,c3] + uDW0[r,c0] + uDW0[r,c2] >= 2*uDW0[r,c2], name = 'udMC3')
          self.model.addConstr(uDX0[r+1,c3] <= uDW0[r,c0] + uDW0[r,c2], name = 'udMC3')
    
    ''' 
    ------------------------------
    Key Schedule 
    ------------------------------
    '''
    uLANE0 = self.model.addVars(16, vtype = GRB.BINARY, name = 'uLANE0')

    trans_pos = [1,7,0,5,2,6,4,3,9,15,8,13,10,14,12,11]

    for i in range(16):
      lane_pos = i
      total_stack_sum = 0
      for r in range(uR):
        self.model.addConstr(uLANE0[i] >= ustk0[r, lane_pos], name = 'uKS')
        total_stack_sum += ustk0[r, lane_pos]
        lane_pos = self.hTable[lane_pos]
      # math.ceil((uR + self.r0)/16) = cancellations from key schedule
      self.model.addConstr((uR)*uLANE0[i] - total_stack_sum <= 
                            (self.Vs - 1)*math.ceil((uR)/30), name = 'uKS')

    # make Key reuse
    for r in range(30, uR):
        for i in range(16):
          self.model.addConstr(ustk0[r-30,i] == ustk0[r,trans_pos[i]], name = 'uKS')
    
    ''' the second diff. trail '''
    uDX1  = self.model.addVars(range(1, uR+1), 16, vtype = GRB.BINARY, name = 'uDX1')
    udx1  = self.model.addVars(range(1, uR+1), 16, vtype = GRB.BINARY, name = 'udx1')
    uDY1  = self.model.addVars(range(1, uR),   16, vtype = GRB.BINARY, name = 'uDY1')
    udy1  = self.model.addVars(range(1, uR),   16, vtype = GRB.BINARY, name = 'udy1')
    uDZ1  = self.model.addVars(range(1, uR),   16, vtype = GRB.BINARY, name = 'uDZ1') 
    udz1  = self.model.addVars(range(1, uR),   16, vtype = GRB.BINARY, name = 'udz1') 
    uDW1  = self.model.addVars(uR,             16, vtype = GRB.BINARY, name = 'uDW1') # W0 = Plaintext
    udw1  = self.model.addVars(uR,             16, vtype = GRB.BINARY, name = 'udw1') # w0 = Plaintext
    ustk1 = self.model.addVars(uR,     16, vtype = GRB.BINARY, name = 'ustk1')

    '''Basic Constrs:'''
    # Remove dx[r,i] > DX[r,i]
    for i in range(16): 
      for r in range(1, uR+1):
        self.model.addConstr(uDX1[r, i] >= udx1[r, i], name = 'uBasic')
      for r in range(1, uR):
        self.model.addConstr(uDY1[r, i] >= udy1[r, i], name = 'uBasic')
        self.model.addConstr(uDZ1[r, i] >= udz1[r, i], name = 'uBasic')
      for r in range(uR):
        self.model.addConstr(uDW1[r, i] >= udw1[r, i], name = 'uBasic')

    '''Operation: SC'''
    # In Extension
    for r in range(1, self.rEb+1):
      for i in range(16):
        self.model.addConstr(uDX1[r,i] == uDY1[r,i], name='ueSC')
        self.model.addConstr(udx1[r,i] == uDY1[r,i], name='ueSC')
    # In distinguisher
    for r in range(self.rEb+1, uR):
      for i in range(16):
        self.model.addConstr(uDY1[r,i] == uDX1[r,i], name='udSC')
        self.model.addConstr(udy1[r,i] == uDX1[r,i], name='udSC')


    '''(*ri-(r0)-r1)Operation: ATK'''
  
    ''' 
    ----------------------------------
    Equivalent Key for the first round
    ----------------------------------
    '''
    eqk1 = self.model.addVars(16, vtype = GRB.BINARY, name = 'eqk_1')

    for c in [[0,4,8,12],[1,5,9,13],[2,6,10,14],[3,7,11,15]]:
        c0, c1, c2, c3 = c[0], c[1], c[2], c[3]
        tc0, tc1, tc2, tc3 = self.SRpv[c[0]], self.SRpv[c[1]], self.SRpv[c[2]], self.SRpv[c[3]]
        
        # --- POSITION -- 0: eqk[c0] = stk1[c3] + eqk[c3] ---
        self.model.addConstr(eqk1[c0] == ustk1[0,tc0], name = '2eqk0')

        # --- POSITION -- 1: eqk[c1] = stk1[c0] ---
        self.model.addConstr(eqk1[c1] == ustk1[0,tc0], name = '2eqk1')
        
        # --- POSITION -- 2: eqk[c2] = stk1[c1] + stk[c2] ---
        self.model.addConstr(eqk1[c2] == ustk1[0,tc1], name = '2eqk2')

        # --- POSITION -- 2: eqk[c3] = stk1[c0] + stk[c2] ---
        self.model.addConstr(eqk1[c3] == ustk1[0,tc0], name = '2eqk3')


    # Add eqk in the first round
    for i in range(16):
      self.model.addConstr(udw1[0,i] == udx1[1,i], name='uAEqk')
      self.model.addConstr(eqk1[i] + uDW1[0,i] - uDX1[1,i] >= 0, name='uAEqk')
      self.model.addConstr(eqk1[i] - uDW1[0,i] + uDX1[1,i] >= 0, name='uAEqk')
      self.model.addConstr(-eqk1[i] + uDW1[0,i] + uDX1[1,i] - udx1[1,i] >= 0, name='uAEqk')
      self.model.addConstr(uDW1[0,i] >= udx1[1,i], name='uAEqk')
    # # In ri
    # for r in range(1, self.ri):
    #   for i in range(8):
    #     self.model.addConstr(udy1[r,i] == udz1[r,i], name='uATKri')
    #     self.model.addConstr(ustk1[r,i] + uDY1[r,i] - uDZ1[r,i] >= 0, name='uATKri')
    #     self.model.addConstr(ustk1[r,i] - uDY1[r,i] + uDZ1[r,i] >= 0, name='uATKri')
    #     self.model.addConstr(-ustk1[r,i] + uDY1[r,i] + uDZ1[r,i] - udz1[r,i] >= 0, name='uATKri')
    #     self.model.addConstr(uDY1[r,i] >= udz1[r,i], name='uATKri')
    #   for i in range(8,16):
    #     self.model.addConstr(uDY1[r,i] == uDZ1[r,i], name='uATKri')
    #     self.model.addConstr(udy1[r,i] == udz1[r,i], name='uATKri')
    # # In r1
    # for r in range(self.ri, uR):
    #   tR = r + self.r0
    #   for i in range(8):
    #     self.model.addConstr(udy1[r,i] == udz1[r,i], name='uATKr1')
    #     self.model.addConstr(ustk1[tR,i] + uDY1[r,i] - uDZ1[r,i] >= 0, name='uATKr1')
    #     self.model.addConstr(ustk1[tR,i] - uDY1[r,i] + uDZ1[r,i] >= 0, name='uATKr1')
    #     self.model.addConstr(-ustk1[tR,i] + uDY1[r,i] + uDZ1[r,i] - udz1[r,i] >= 0, name='uATKr1')
    #     self.model.addConstr(uDY1[r,i] >= udz1[r,i], name='uATKr1')
    #   for i in range(8,16):
    #     self.model.addConstr(uDY1[r,i] == uDZ1[r,i], name='uATKr1')
    #     self.model.addConstr(udy1[r,i] == udz1[r,i], name='uATKr1')
    for r in range(1, uR):
      for i in range(8):
        self.model.addConstr(udy1[r,i] == udz1[r,i], name='uATK')
        self.model.addConstr(ustk1[r,i] + uDY1[r,i] - uDZ1[r,i] >= 0, name='uATK')
        self.model.addConstr(ustk1[r,i] - uDY1[r,i] + uDZ1[r,i] >= 0, name='uATK')
        self.model.addConstr(-ustk1[r,i] + uDY1[r,i] + uDZ1[r,i] - udz1[r,i] >= 0, name='uATK')
        self.model.addConstr(uDY1[r,i] >= udz1[r,i], name='uATK')
      for i in range(8,16):
        self.model.addConstr(uDY1[r,i] == uDZ1[r,i], name='uATK')
        self.model.addConstr(udy1[r,i] == udz1[r,i], name='uATK')


    '''(*)Operation: SR''' 
    for r in range(1,uR):
      for i in range(16):
        self.model.addConstr(uDZ1[r,i] == uDW1[r,self.SRp[i]], name = 'uSR')
        self.model.addConstr(udz1[r,i] == udw1[r,self.SRp[i]], name = 'uSR')

    '''Operation: MC'''
    # ************ In extension (backword) ************
    for r in range(1, self.rEb):
      for c in [[0,4,8,12],[1,5,9,13],[2,6,10,14],[3,7,11,15]]:
        c0, c1, c2, c3 = c[0], c[1], c[2], c[3]

        # --- POSITION -- 0: W[c0] = X[c1] ---
        self.model.addConstr(udw1[r,c0] == udx1[r+1,c1], name = 'ueMC0')

        self.model.addConstr(uDW1[r,c0] == uDX1[r+1,c1], name = 'ueMC0')

        # --- POSITION -- 1: W[c1] = W[c2] + X[c2] ---
        self.model.addConstr(udw1[r,c1] >= udw1[r,c2], name = 'ueMC1')
        self.model.addConstr(udw1[r,c1] >= udx1[r+1,c2], name = 'ueMC1')
        self.model.addConstr(udw1[r,c1] <= udw1[r,c2] + udx1[r+1,c2], name = 'ueMC1')

        self.model.addConstr(uDW1[r,c1] + uDW1[r,c2] + uDX1[r+1,c2] >= 2*uDW1[r,c1], name = 'ueMC1')
        self.model.addConstr(uDW1[r,c1] + uDW1[r,c2] + uDX1[r+1,c2] >= 2*uDW1[r,c2], name = 'ueMC1')
        self.model.addConstr(uDW1[r,c1] + uDW1[r,c2] + uDX1[r+1,c2] >= 2*uDX1[r+1,c2], name = 'ueMC1')
        self.model.addConstr(uDW1[r,c1] <= uDW1[r,c2] + uDX1[r+1,c2], name = 'ueMC1')

        # --- POSITION -- 2: W[c2] = X[c1] + X[c3] ---
        self.model.addConstr(udw1[r,c2] >= udx1[r+1,c1], name = 'ueMC2')
        self.model.addConstr(udw1[r,c2] >= udx1[r+1,c3], name = 'ueMC2')
        self.model.addConstr(udw1[r,c2] <= udx1[r+1,c1] + udx1[r+1,c3], name = 'ueMC2')

        self.model.addConstr(uDW1[r,c2] + uDX1[r+1,c1] + uDX1[r+1,c3] >= 2*uDW1[r,c2], name = 'ueMC2')
        self.model.addConstr(uDW1[r,c2] + uDX1[r+1,c1] + uDX1[r+1,c3] >= 2*uDX1[r+1,c1], name = 'ueMC2')
        self.model.addConstr(uDW1[r,c2] + uDX1[r+1,c1] + uDX1[r+1,c3] >= 2*uDX1[r+1,c3], name = 'ueMC2')
        self.model.addConstr(uDW1[r,c2] <= uDX1[r+1,c1] + uDX1[r+1,c3], name = 'ueMC2')

        # --- POSITION -- 3: W[c3] = X[c0] + X[c3] ---
        self.model.addConstr(udw1[r,c3] >= udx1[r+1,c0], name = 'ueMC3')
        self.model.addConstr(udw1[r,c3] >= udx1[r+1,c3], name = 'ueMC3')
        self.model.addConstr(udw1[r,c3] <= udx1[r+1,c0] + udx1[r+1,c3], name = 'ueMC3')

        self.model.addConstr(uDW1[r,c3] + uDX1[r+1,c0] + uDX1[r+1,c3] >= 2*uDW1[r,c3], name = 'ueMC3')
        self.model.addConstr(uDW1[r,c3] + uDX1[r+1,c0] + uDX1[r+1,c3] >= 2*uDX1[r+1,c0], name = 'ueMC3')
        self.model.addConstr(uDW1[r,c3] + uDX1[r+1,c0] + uDX1[r+1,c3] >= 2*uDX1[r+1,c3], name = 'ueMC3')
        self.model.addConstr(uDW1[r,c3] <= uDX1[r+1,c0] + uDX1[r+1,c3], name = 'ueMC3')
    
    # ************ In distinguisher (forward) ************
    for r in range(self.rEb, uR):
      for c in [[0,4,8,12],[1,5,9,13],[2,6,10,14],[3,7,11,15]]:
          c0, c1, c2, c3 = c[0], c[1], c[2], c[3]

          # --- POSITION -- 0: X[c0] = X[c3] + W[c3] ---
          self.model.addConstr(udx1[r+1,c0] >= udw1[r,c3], name = 'udMC0')
          self.model.addConstr(udx1[r+1,c0] >= udx1[r+1,c3], name = 'udMC0')
          self.model.addConstr(udx1[r+1,c0] <= udw1[r,c3] + udx1[r+1,c3], name = 'udMC0')

          self.model.addConstr(uDX1[r+1,c0] + uDW1[r,c3] + uDX1[r+1,c3] >= 2*uDX1[r+1,c0], name = 'udMC0')
          self.model.addConstr(uDX1[r+1,c0] + uDW1[r,c3] + uDX1[r+1,c3] >= 2*uDW1[r,c3], name = 'udMC0')
          self.model.addConstr(uDX1[r+1,c0] + uDW1[r,c3] + uDX1[r+1,c3] >= 2*uDX1[r+1,c3], name = 'udMC0')
          self.model.addConstr(uDX1[r+1,c0] <= uDW1[r,c3] + uDX1[r+1,c3], name = 'udMC0')

          # --- POSITION -- 1: X[c1] = W[c0] ---
          self.model.addConstr(udx1[r+1,c1] == udw1[r,c0], name = 'udMC1')

          self.model.addConstr(uDX1[r+1,c1] == uDW1[r,c0], name = 'udMC1')

          # --- POSITION -- 2: X[c2] = X[c1] + W[c2] ---
          self.model.addConstr(udx1[r+1,c2] >= udw1[r,c1], name = 'udMC2')
          self.model.addConstr(udx1[r+1,c2] >= udw1[r,c2], name = 'udMC2')
          self.model.addConstr(udx1[r+1,c2] <= udw1[r,c1] + udw1[r,c2], name = 'udMC2')

          self.model.addConstr(uDX1[r+1,c2] + uDW1[r,c1] + uDW1[r,c2] >= 2*uDX1[r+1,c2], name = 'udMC2')
          self.model.addConstr(uDX1[r+1,c2] + uDW1[r,c1] + uDW1[r,c2] >= 2*uDW1[r,c1], name = 'udMC2')
          self.model.addConstr(uDX1[r+1,c2] + uDW1[r,c1] + uDW1[r,c2] >= 2*uDW1[r,c2], name = 'udMC2')
          self.model.addConstr(uDX1[r+1,c2] <= uDW1[r,c1] + uDW1[r,c2], name = 'udMC2')

          # --- POSITION -- 3: X[c3] = W[c0] + W[c2] ---
          self.model.addConstr(udx1[r+1,c3] >= udw1[r,c0], name = 'udMC3')
          self.model.addConstr(udx1[r+1,c3] >= udw1[r,c2], name = 'udMC3')
          self.model.addConstr(udx1[r+1,c3] <= udw1[r,c0] + udw1[r,c2], name = 'udMC3')

          self.model.addConstr(uDX1[r+1,c3] + uDW1[r,c0] + uDW1[r,c2] >= 2*uDX1[r+1,c3], name = 'udMC3')
          self.model.addConstr(uDX1[r+1,c3] + uDW1[r,c0] + uDW1[r,c2] >= 2*uDW1[r,c0], name = 'udMC3')
          self.model.addConstr(uDX1[r+1,c3] + uDW1[r,c0] + uDW1[r,c2] >= 2*uDW1[r,c2], name = 'udMC3')
          self.model.addConstr(uDX1[r+1,c3] <= uDW1[r,c0] + uDW1[r,c2], name = 'udMC3')
    
    ''' 
    ------------------------------
    Key Schedule 
    ------------------------------
    '''
    uLANE1 = self.model.addVars(16, vtype = GRB.BINARY, name = 'uLANE1')

    trans_pos = [1,7,0,5,2,6,4,3,9,15,8,13,10,14,12,11]

    for i in range(16):
      lane_pos = i
      total_stack_sum = 0
      for r in range(uR):
        self.model.addConstr(uLANE1[i] >= ustk1[r, lane_pos], name = 'uKS')
        total_stack_sum += ustk1[r, lane_pos]
        lane_pos = self.hTable[lane_pos]
      # math.ceil((uR + self.r0)/16) = cancellations from key schedule
      self.model.addConstr((uR)*uLANE1[i] - total_stack_sum <= 
                            (self.Vs - 1)*math.ceil((uR)/30), name = 'uKS')

    # make Key reuse
    for r in range(30, uR):
        for i in range(16):
          self.model.addConstr(ustk1[r-30,i] == ustk1[r,trans_pos[i]], name = 'uKS')

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
    ''' The first diff. trail '''
    # NOTE: Assert rEb <= ri and rEb+rEu >= ri
    ugEQK0 = self.model.addVars(range(self.rEb),      16, vtype = GRB.BINARY, name = 'ugEQK0')
    uGstk0 = self.model.addVars(range(self.rEb),      16, vtype = GRB.BINARY, name = 'uGstk0')
    for r in range(self.rEb): 
      for i in range(8,16):
        self.model.addConstr(uGstk0[r,i] == 0, name = 'uGstk0L')
    uDetX0 = self.model.addVars(range(1, self.rEb+1), 16, vtype = GRB.BINARY, name = 'uDetX0')
    uDetY0 = self.model.addVars(range(1, self.rEb+1), 16, vtype = GRB.BINARY, name = 'uDetY0')
    uDetZ0 = self.model.addVars(range(1, self.rEb),   16, vtype = GRB.BINARY, name = 'uDetZ0')
    uDetW0 = self.model.addVars(range(0, self.rEb),   16, vtype = GRB.BINARY, name = 'uDetW0')
    uFrSC0 = self.model.addVars(range(1, self.rEb+1), 16, vtype = GRB.BINARY, name = 'uFrSC0')
    uFrMC0 = self.model.addVars(range(2, self.rEb+1), 16, vtype = GRB.BINARY, name = 'uFrMC0')
    #-----------------------------------------------------------------------------------------
    ugeqka0 = self.model.addVars(range(self.rEb),      16, vtype = GRB.BINARY, name = 'ugeqka0')
    ugstka0 = self.model.addVars(range(self.rEb),      16, vtype = GRB.BINARY, name = 'ugstka0')
    for r in range(self.rEb): 
      for i in range(8,16):
        self.model.addConstr(ugstka0[r,i] == 0, 'ugstkLa')
    udetXa0 = self.model.addVars(range(1, self.rEb+1), 16, vtype = GRB.BINARY, name = 'udetXa0')
    udetYa0 = self.model.addVars(range(1, self.rEb+1), 16, vtype = GRB.BINARY, name = 'udetYa0')
    udetZa0 = self.model.addVars(range(1, self.rEb),   16, vtype = GRB.BINARY, name = 'udetZa0')
    udetWa0 = self.model.addVars(range(0, self.rEb),   16, vtype = GRB.BINARY, name = 'udetWa0')
    ufrSCa0 = self.model.addVars(range(1, self.rEb+1), 16, vtype = GRB.BINARY, name = 'ufrSCa0')
    ufrMCa0 = self.model.addVars(range(2, self.rEb+1), 16, vtype = GRB.BINARY, name = 'ufrMCa0')

    self.model.addConstr(sum(uDetW0[0,i] for i in range(16)) == 16, name = 'VP') # Set plaintext is determined
    for i in range(16): self.model.addConstr(uDetX0[1,i] == uDetW0[0,i], name = 'W0X0a')
    # ------------------------------------------------------------------------
    self.model.addConstr(sum(udetWa0[0,i] for i in range(16)) == 16, name = 'VPa') # Set plaintext is determined
    for i in range(16): self.model.addConstr(udetXa0[1,i] == uDetW0[0,i], name = 'W0X0a')

    # ----------------- Guessing Eqk -----------------------
    # NOTE: That is for the key bridge.
    for r in range(self.rEb):
      for c in [[0,4,8,12],[1,5,9,13],[2,6,10,14],[3,7,11,15]]:
          c0, c1, c2, c3 = c[0], c[1], c[2], c[3]
          tc0, tc1 = self.SRpv[c[0]], self.SRpv[c[1]]

          # POSITION -- 0: ugEQK0[0,1,2,3 (c0)] = uGstk0[0,1,2,3 -> (tc0)]
          self.model.addConstr(ugEQK0[r,c0] == uGstk0[r,tc0], name = 'geqk0')
          self.model.addConstr(ugeqka0[r,c0] == ugstka0[r,tc0], name = 'geqk0A')

          # POSITION -- 1: ugEQK0[4,5,6,7 (c1)] = uGstk0[0,1,2,3 -> (tc0)] 
          self.model.addConstr(ugEQK0[r,c1] == uGstk0[r,tc0], name = 'geqk1')
          self.model.addConstr(ugeqka0[r,c1] == ugstka0[r,tc0], name = 'geqk1A')

          # POSITION -- 2: ugEQK0[8,9,10,11] = uGstk0[7,4,5,6 -> (tc1)] 
          self.model.addConstr(ugEQK0[r,c2] == uGstk0[r,tc1], name = 'geqk2')
          self.model.addConstr(ugeqka0[r,c2] == ugstka0[r,tc1], name = 'geqk2A')

          # POSITION -- 2: eqk[c3] = stk1[12,13,14,15] + stk[0,1,2,3 -> (tc0)] 
          self.model.addConstr(ugEQK0[r,c3] == uGstk0[r,tc0], name = 'geqk3')
          self.model.addConstr(ugeqka0[r,c3] == ugstka0[r,tc0], name = 'geqk3A')
    # ----------------- Guessing Eqk -----------------------

    ''' KR: Determine and Filter in SC '''
    for r in range(1, self.rEb+1):
      for i in range(16):
        # 'Det.' propagation via SC (AND ADD EQK)
        self.model.addConstr(uDetY0[r,i] <= uDetX0[r,i], name = 'udetSC+ek')
        self.model.addConstr(uDetY0[r,i] <= ugEQK0[r-1,i], name = 'udetSC+ek')
        self.model.addConstr( - uDetX0[r,i] - ugEQK0[r-1,i] + uDetY0[r,i] + 1 >= 0, name = 'udetSC+ek')
        # ------------------------------------------------------------------------------------------
        self.model.addConstr(udetYa0[r,i] <= udetXa0[r,i], name = 'udetSC+ekA')
        self.model.addConstr(udetYa0[r,i] <= ugeqka0[r-1,i], name = 'udetSC+ekA')
        self.model.addConstr( - udetXa0[r,i] - ugeqka0[r-1,i] + udetYa0[r,i] + 1 >= 0, name = 'udetSC+ekA')

        # Filter obtain from SC
        self.model.addConstr(uDY0[r,i] - udy0[r,i] - uFrSC0[r,i] >= 0, name = 'ufrSC')
        self.model.addConstr(uDetY0[r,i] >= uFrSC0[r,i], name = 'ufrSC')
        self.model.addConstr(- uDetY0[r,i] - uDY0[r,i] + udy0[r,i] + uFrSC0[r,i] >= -1, name = 'ufrSC')
        # ------------------------------------------------------------------------------------------
        self.model.addConstr(uDY0[r,i] - udy0[r,i] - ufrSCa0[r,i] >= 0, name = 'ufrSCa0')
        self.model.addConstr(udetYa0[r,i] >= ufrSCa0[r,i], name = 'ufrSCa0')
        self.model.addConstr(- udetYa0[r,i] - uDY0[r,i] + udy0[r,i] + ufrSCa0[r,i] >= -1, name = 'ufrSCa0')

    ''' 'Det.' trans Y to Z '''
    for r in range(1, self.rEb):
      for i in range(16):
        self.model.addConstr(uDetY0[r,i] == uDetZ0[r,i], 'udetY2Z')
        # -------------------------------------------------------
        self.model.addConstr(udetYa0[r,i] == udetZa0[r,i], 'udety2zA')

    ''' KR: 'Det.' propagation in SR '''
    for r in range(1, self.rEb):
      for i in range(16): 
        self.model.addConstr(uDetW0[r,i] == uDetZ0[r,self.SRpv[i]], name = 'udetSR')
        # -------------------------------------------------------------------------
        self.model.addConstr(udetWa0[r,i] == udetZa0[r,self.SRpv[i]], name = 'udetsrA')

    ''' KR: 'Det.' propagation in MC '''
    for r in range(1, self.rEb):
      for c in [[0,4,8,12],[1,5,9,13],[2,6,10,14],[3,7,11,15]]:
          c0, c1, c2, c3 = c[0], c[1], c[2], c[3]

          # ===============================================================================================
          # POSITION -- 0: DetX[c0] = DetW[3] & DetX[c3]
          self.model.addConstr(uDetX0[r+1,c0] <= uDetW0[r,c3], name = 'udetMC0')
          self.model.addConstr(uDetX0[r+1,c0] <= uDetX0[r+1,c3], name = 'udetMC0')
          self.model.addConstr( - uDetW0[r,c3] - uDetX0[r+1,c3] + uDetX0[r+1,c0] + 1 >= 0, name = 'udetMC0')
          # --------------------------------------------------------------------------------------------
          self.model.addConstr(udetXa0[r+1,c0] <= udetWa0[r,c3], name = 'udetmc0A')
          self.model.addConstr(udetXa0[r+1,c0] <= udetXa0[r+1,c3], name = 'udetmc0A')
          self.model.addConstr( - udetWa0[r,c3] - udetXa0[r+1,c3] + udetXa0[r+1,c0] + 1 >= 0, name = 'udetmc0A')

          # POSITION -- 1: DetX[c1] = DetW[c0]
          self.model.addConstr(uDetX0[r+1,c1] == uDetW0[r,c0], name = 'udetMC1')
          # ---------------------------------------------------------------------
          self.model.addConstr(udetXa0[r+1,c1] == udetWa0[r,c0], name = 'udetmc1A')

          # POSITION -- 2: DetX[c2] = DetW[c1] & DetW[c2]
          self.model.addConstr(uDetX0[r+1,c2] <= uDetW0[r,c1], name = 'udetMC2')
          self.model.addConstr(uDetX0[r+1,c2] <= uDetW0[r,c2], name = 'udetMC2')
          self.model.addConstr( - uDetW0[r,c1] - uDetW0[r,c2] + uDetX0[r+1,c2] + 1 >= 0, name = 'udetMC2')
          # --------------------------------------------------------------------------------------------
          self.model.addConstr(udetXa0[r+1,c2] <= udetWa0[r,c1], name = 'udetmc2A')
          self.model.addConstr(udetXa0[r+1,c2] <= udetWa0[r,c2], name = 'udetmc2A')
          self.model.addConstr( - udetWa0[r,c1] - udetWa0[r,c2] + udetXa0[r+1,c2] + 1 >= 0, name = 'udetmc2A')

          # POSITION -- 3: DetX[c3] = DetW[c0] & DetW[c2]
          self.model.addConstr(uDetX0[r+1,c3] <= uDetW0[r,c0], name = 'udetMC3')
          self.model.addConstr(uDetX0[r+1,c3] <= uDetW0[r,c2], name = 'udetMC3')
          self.model.addConstr( - uDetW0[r,c0] - uDetW0[r,c2] + uDetX0[r+1,c3] + 1 >= 0, name = 'udetMC3')
          # --------------------------------------------------------------------------------------------
          self.model.addConstr(udetXa0[r+1,c3] <= udetWa0[r,c0], name = 'udetmc3A')
          self.model.addConstr(udetXa0[r+1,c3] <= udetWa0[r,c2], name = 'udetmc3A')
          self.model.addConstr( - udetWa0[r,c0] - udetWa0[r,c2] + udetXa0[r+1,c3] + 1 >= 0, name = 'udetmc3A')
          # ===============================================================================================

    ''' KR: 'Filter' obtain from MC '''
    for r in range(2,self.rEb+1):
      for c in [[0,4,8,12],[1,5,9,13],[2,6,10,14],[3,7,11,15]]:
          c0, c1, c2, c3 = c[0], c[1], c[2], c[3]

          # POSITION -- 0: DetX[c0] = DetW[3] & DetX[c3]
          self.model.addConstr(uFrMC0[r,c0] <= uDetX0[r,c0], name = 'ufrMC0')
          self.model.addConstr(uFrMC0[r,c0] <= udx0[r,c3], name = 'ufrMC0')
          self.model.addConstr(uFrMC0[r,c0] <= udw0[r-1,c3], name = 'ufrMC0')
          self.model.addConstr(uFrMC0[r,c0] <= 1 - udx0[r,c0], name = 'ufrMC0')
          self.model.addConstr(uFrMC0[r,c0] >= uDetX0[r,i] + udx0[r,c3] + udw0[r-1,c3] + (1-udx0[r,c0]) - 3, name = 'ufrMC0')
          # ------------------------------------------------------------------------------------------------------------
          self.model.addConstr(ufrMCa0[r,c0] <= udetXa0[r,c0], name = 'ufrmc0A')
          self.model.addConstr(ufrMCa0[r,c0] <= udx0[r,c3], name = 'ufrmc0A')
          self.model.addConstr(ufrMCa0[r,c0] <= udw0[r-1,c3], name = 'ufrmc0A')
          self.model.addConstr(ufrMCa0[r,c0] <= 1 - udx0[r,c0], name = 'ufrmc0A')
          self.model.addConstr(ufrMCa0[r,c0] >= udetXa0[r,i] + udx0[r,c3] + udw0[r-1,c3] + (1-udx0[r,c0]) - 3, name = 'ufrmc0A')

          # POSITION -- 1: DetX[c1] = DetW[c0] (No filter in this position)
          self.model.addConstr(uFrMC0[r,c1] == 0, 'ufrMC1')
          # ----------------------------------------------
          self.model.addConstr(ufrMCa0[r,c1] == 0, 'ufrmc1A')

          # POSITION -- 2: DetX[c2] = DetW[c1] & DetW[c2]
          self.model.addConstr(uFrMC0[r,c2] <= uDetX0[r,c2], name = 'ufrMC2')
          self.model.addConstr(uFrMC0[r,c2] <= udw0[r-1,c1], name = 'ufrMC2')
          self.model.addConstr(uFrMC0[r,c2] <= udw0[r-1,c2], name = 'ufrMC2')
          self.model.addConstr(uFrMC0[r,c2] <= 1 - udx0[r,c2], name = 'ufrMC2')
          self.model.addConstr(uFrMC0[r,c2] >= uDetX0[r,c2] + udw0[r-1,c1] + udw0[r-1,c2] + (1-udx0[r,c2]) - 3, name = 'ufrMC2')
          # ------------------------------------------------------------------------------------------------------------
          self.model.addConstr(ufrMCa0[r,c2] <= udetXa0[r,c2], name = 'ufrmc2A')
          self.model.addConstr(ufrMCa0[r,c2] <= udw0[r-1,c1], name = 'ufrmc2A')
          self.model.addConstr(ufrMCa0[r,c2] <= udw0[r-1,c2], name = 'ufrmc2A')
          self.model.addConstr(ufrMCa0[r,c2] <= 1 - udx0[r,c2], name = 'ufrmc2A')
          self.model.addConstr(ufrMCa0[r,c2] >= udetXa0[r,c2] + udw0[r-1,c1] + udw0[r-1,c2] + (1-udx0[r,c2]) - 3, name = 'ufrmc2A')

          # POSITION -- 3: DetX[c3] = DetW[c0] & DetW[c2]
          self.model.addConstr(uFrMC0[r,c3] <= uDetX0[r,c3], name = 'ufrMC3')
          self.model.addConstr(uFrMC0[r,c3] <= udw0[r-1,c0], name = 'ufrMC3')
          self.model.addConstr(uFrMC0[r,c3] <= udw0[r-1,c2], name = 'ufrMC3')
          self.model.addConstr(uFrMC0[r,c3] <= 1 - udx0[r,c3], name = 'ufrMC3')
          self.model.addConstr(uFrMC0[r,c3] >= uDetX0[r,c3] + udw0[r-1,c0] + udw0[r-1,c2] + (1-udx0[r,c2]) - 3, name = 'ufrMC3')
          # ------------------------------------------------------------------------------------------------------------
          self.model.addConstr(ufrMCa0[r,c3] <= udetXa0[r,c3], name = 'ufrmc3A')
          self.model.addConstr(ufrMCa0[r,c3] <= udw0[r-1,c0], name = 'ufrmc3A')
          self.model.addConstr(ufrMCa0[r,c3] <= udw0[r-1,c2], name = 'ufrmc3A')
          self.model.addConstr(ufrMCa0[r,c3] <= 1 - udx0[r,c3], name = 'ufrmc3A')
          self.model.addConstr(ufrMCa0[r,c3] >= udetXa0[r,c3] + udw0[r-1,c0] + udw0[r-1,c2] + (1-udx0[r,c2]) - 3, name = 'ufrmc3A')

    # Setting pre-guess all involved key in auxiliary upper trail (rf = friter)
    # -------------------------------------------------------------------------------------
    self.model.addConstr(sum(udw0[0,i] for i in range(16)) == 
                         (sum(ufrSCa0[r,i] for r in range(1, self.rEb+1) for i in range(16)) + 
                          sum(ufrMCa0[r,i] for r in range(2, self.rEb+1) for i in range(16))),
                          name = 'rfUfr')
    # -------------------------------------------------------------------------------------

    ''' The second diff. trail '''
    # NOTE: Assert rEb <= ri and rEb+rEu >= ri
    ugEQK1 = self.model.addVars(range(self.rEb),      16, vtype = GRB.BINARY, name = 'ugEQK1')
    uGstk1 = self.model.addVars(range(self.rEb),      16, vtype = GRB.BINARY, name = 'uGstk1')
    for r in range(self.rEb): 
      for i in range(8,16):
        self.model.addConstr(uGstk1[r,i] == 0, name = 'uGstk1L')
    uDetX1 = self.model.addVars(range(1, self.rEb+1), 16, vtype = GRB.BINARY, name = 'uDetX1')
    uDetY1 = self.model.addVars(range(1, self.rEb+1), 16, vtype = GRB.BINARY, name = 'uDetY1')
    uDetZ1 = self.model.addVars(range(1, self.rEb),   16, vtype = GRB.BINARY, name = 'uDetZ1')
    uDetW1 = self.model.addVars(range(0, self.rEb),   16, vtype = GRB.BINARY, name = 'uDetW1')
    uFrSC1 = self.model.addVars(range(1, self.rEb+1), 16, vtype = GRB.BINARY, name = 'uFrSC1')
    uFrMC1 = self.model.addVars(range(2, self.rEb+1), 16, vtype = GRB.BINARY, name = 'uFrMC1')
    #-----------------------------------------------------------------------------------------
    ugeqka1 = self.model.addVars(range(self.rEb),      16, vtype = GRB.BINARY, name = 'ugeqka1')
    ugstka1 = self.model.addVars(range(self.rEb),      16, vtype = GRB.BINARY, name = 'ugstka1')
    for r in range(self.rEb): 
      for i in range(8,16):
        self.model.addConstr(ugstka1[r,i] == 0, 'ugstkLa')
    udetXa1 = self.model.addVars(range(1, self.rEb+1), 16, vtype = GRB.BINARY, name = 'udetXa1')
    udetYa1 = self.model.addVars(range(1, self.rEb+1), 16, vtype = GRB.BINARY, name = 'udetYa1')
    udetZa1 = self.model.addVars(range(1, self.rEb),   16, vtype = GRB.BINARY, name = 'udetZa1')
    udetWa1 = self.model.addVars(range(0, self.rEb),   16, vtype = GRB.BINARY, name = 'udetWa1')
    ufrSCa1 = self.model.addVars(range(1, self.rEb+1), 16, vtype = GRB.BINARY, name = 'ufrSCa1')
    ufrMCa1 = self.model.addVars(range(2, self.rEb+1), 16, vtype = GRB.BINARY, name = 'ufrMCa1')

    self.model.addConstr(sum(uDetW1[0,i] for i in range(16)) == 16, name = 'VP') # Set plaintext is determined
    for i in range(16): self.model.addConstr(uDetX1[1,i] == uDetW1[0,i], name = 'W0X0a')
    # ------------------------------------------------------------------------
    self.model.addConstr(sum(udetWa1[0,i] for i in range(16)) == 16, name = 'VPa') # Set plaintext is determined
    for i in range(16): self.model.addConstr(udetXa1[1,i] == uDetW1[0,i], name = 'W0X0a')

    # ----------------- Guessing Eqk -----------------------
    # NOTE: That is for the key bridge.
    for r in range(self.rEb):
      for c in [[0,4,8,12],[1,5,9,13],[2,6,10,14],[3,7,11,15]]:
          c0, c1, c2, c3 = c[0], c[1], c[2], c[3]
          tc0, tc1 = self.SRpv[c[0]], self.SRpv[c[1]]

          # POSITION -- 0: ugEQK1[0,1,2,3 (c0)] = uGstk1[0,1,2,3 -> (tc0)]
          self.model.addConstr(ugEQK1[r,c0] == uGstk1[r,tc0], name = 'geqk0')
          self.model.addConstr(ugeqka1[r,c0] == ugstka1[r,tc0], name = 'geqk0A')

          # POSITION -- 1: ugEQK1[4,5,6,7 (c1)] = uGstk1[0,1,2,3 -> (tc0)] 
          self.model.addConstr(ugEQK1[r,c1] == uGstk1[r,tc0], name = 'geqk1')
          self.model.addConstr(ugeqka1[r,c1] == ugstka1[r,tc0], name = 'geqk1A')

          # POSITION -- 2: ugEQK1[8,9,10,11] = uGstk1[7,4,5,6 -> (tc1)] 
          self.model.addConstr(ugEQK1[r,c2] == uGstk1[r,tc1], name = 'geqk2')
          self.model.addConstr(ugeqka1[r,c2] == ugstka1[r,tc1], name = 'geqk2A')

          # POSITION -- 2: eqk[c3] = stk1[12,13,14,15] + stk[0,1,2,3 -> (tc0)] 
          self.model.addConstr(ugEQK1[r,c3] == uGstk1[r,tc0], name = 'geqk3')
          self.model.addConstr(ugeqka1[r,c3] == ugstka1[r,tc0], name = 'geqk3A')
    # ----------------- Guessing Eqk -----------------------

    ''' KR: Determine and Filter in SC '''
    for r in range(1, self.rEb+1):
      for i in range(16):
        # 'Det.' propagation via SC (AND ADD EQK)
        self.model.addConstr(uDetY1[r,i] <= uDetX1[r,i], name = 'udetSC+ek')
        self.model.addConstr(uDetY1[r,i] <= ugEQK1[r-1,i], name = 'udetSC+ek')
        self.model.addConstr( - uDetX1[r,i] - ugEQK1[r-1,i] + uDetY1[r,i] + 1 >= 0, name = 'udetSC+ek')
        # ------------------------------------------------------------------------------------------
        self.model.addConstr(udetYa1[r,i] <= udetXa1[r,i], name = 'udetSC+ekA')
        self.model.addConstr(udetYa1[r,i] <= ugeqka1[r-1,i], name = 'udetSC+ekA')
        self.model.addConstr( - udetXa1[r,i] - ugeqka1[r-1,i] + udetYa1[r,i] + 1 >= 0, name = 'udetSC+ekA')

        # Filter obtain from SC
        self.model.addConstr(uDY1[r,i] - udy1[r,i] - uFrSC1[r,i] >= 0, name = 'ufrSC')
        self.model.addConstr(uDetY1[r,i] >= uFrSC1[r,i], name = 'ufrSC')
        self.model.addConstr(- uDetY1[r,i] - uDY1[r,i] + udy1[r,i] + uFrSC1[r,i] >= -1, name = 'ufrSC')
        # ------------------------------------------------------------------------------------------
        self.model.addConstr(uDY1[r,i] - udy1[r,i] - ufrSCa1[r,i] >= 0, name = 'ufrSCa1')
        self.model.addConstr(udetYa1[r,i] >= ufrSCa1[r,i], name = 'ufrSCa1')
        self.model.addConstr(- udetYa1[r,i] - uDY1[r,i] + udy1[r,i] + ufrSCa1[r,i] >= -1, name = 'ufrSCa1')

    ''' 'Det.' trans Y to Z '''
    for r in range(1, self.rEb):
      for i in range(16):
        self.model.addConstr(uDetY1[r,i] == uDetZ1[r,i], 'udetY2Z')
        # -------------------------------------------------------
        self.model.addConstr(udetYa1[r,i] == udetZa1[r,i], 'udety2zA')

    ''' KR: 'Det.' propagation in SR '''
    for r in range(1, self.rEb):
      for i in range(16): 
        self.model.addConstr(uDetW1[r,i] == uDetZ1[r,self.SRpv[i]], name = 'udetSR')
        # -------------------------------------------------------------------------
        self.model.addConstr(udetWa1[r,i] == udetZa1[r,self.SRpv[i]], name = 'udetsrA')

    ''' KR: 'Det.' propagation in MC '''
    for r in range(1, self.rEb):
      for c in [[0,4,8,12],[1,5,9,13],[2,6,10,14],[3,7,11,15]]:
          c0, c1, c2, c3 = c[0], c[1], c[2], c[3]

          # ===============================================================================================
          # POSITION -- 0: DetX[c0] = DetW[3] & DetX[c3]
          self.model.addConstr(uDetX1[r+1,c0] <= uDetW1[r,c3], name = 'udetMC0')
          self.model.addConstr(uDetX1[r+1,c0] <= uDetX1[r+1,c3], name = 'udetMC0')
          self.model.addConstr( - uDetW1[r,c3] - uDetX1[r+1,c3] + uDetX1[r+1,c0] + 1 >= 0, name = 'udetMC0')
          # --------------------------------------------------------------------------------------------
          self.model.addConstr(udetXa1[r+1,c0] <= udetWa1[r,c3], name = 'udetmc0A')
          self.model.addConstr(udetXa1[r+1,c0] <= udetXa1[r+1,c3], name = 'udetmc0A')
          self.model.addConstr( - udetWa1[r,c3] - udetXa1[r+1,c3] + udetXa1[r+1,c0] + 1 >= 0, name = 'udetmc0A')

          # POSITION -- 1: DetX[c1] = DetW[c0]
          self.model.addConstr(uDetX1[r+1,c1] == uDetW1[r,c0], name = 'udetMC1')
          # ---------------------------------------------------------------------
          self.model.addConstr(udetXa1[r+1,c1] == udetWa1[r,c0], name = 'udetmc1A')

          # POSITION -- 2: DetX[c2] = DetW[c1] & DetW[c2]
          self.model.addConstr(uDetX1[r+1,c2] <= uDetW1[r,c1], name = 'udetMC2')
          self.model.addConstr(uDetX1[r+1,c2] <= uDetW1[r,c2], name = 'udetMC2')
          self.model.addConstr( - uDetW1[r,c1] - uDetW1[r,c2] + uDetX1[r+1,c2] + 1 >= 0, name = 'udetMC2')
          # --------------------------------------------------------------------------------------------
          self.model.addConstr(udetXa1[r+1,c2] <= udetWa1[r,c1], name = 'udetmc2A')
          self.model.addConstr(udetXa1[r+1,c2] <= udetWa1[r,c2], name = 'udetmc2A')
          self.model.addConstr( - udetWa1[r,c1] - udetWa1[r,c2] + udetXa1[r+1,c2] + 1 >= 0, name = 'udetmc2A')

          # POSITION -- 3: DetX[c3] = DetW[c0] & DetW[c2]
          self.model.addConstr(uDetX1[r+1,c3] <= uDetW1[r,c0], name = 'udetMC3')
          self.model.addConstr(uDetX1[r+1,c3] <= uDetW1[r,c2], name = 'udetMC3')
          self.model.addConstr( - uDetW1[r,c0] - uDetW1[r,c2] + uDetX1[r+1,c3] + 1 >= 0, name = 'udetMC3')
          # --------------------------------------------------------------------------------------------
          self.model.addConstr(udetXa1[r+1,c3] <= udetWa1[r,c0], name = 'udetmc3A')
          self.model.addConstr(udetXa1[r+1,c3] <= udetWa1[r,c2], name = 'udetmc3A')
          self.model.addConstr( - udetWa1[r,c0] - udetWa1[r,c2] + udetXa1[r+1,c3] + 1 >= 0, name = 'udetmc3A')
          # ===============================================================================================

    ''' KR: 'Filter' obtain from MC '''
    for r in range(2,self.rEb+1):
      for c in [[0,4,8,12],[1,5,9,13],[2,6,10,14],[3,7,11,15]]:
          c0, c1, c2, c3 = c[0], c[1], c[2], c[3]

          # POSITION -- 0: DetX[c0] = DetW[3] & DetX[c3]
          self.model.addConstr(uFrMC1[r,c0] <= uDetX1[r,c0], name = 'ufrMC0')
          self.model.addConstr(uFrMC1[r,c0] <= udx1[r,c3], name = 'ufrMC0')
          self.model.addConstr(uFrMC1[r,c0] <= udw1[r-1,c3], name = 'ufrMC0')
          self.model.addConstr(uFrMC1[r,c0] <= 1 - udx1[r,c0], name = 'ufrMC0')
          self.model.addConstr(uFrMC1[r,c0] >= uDetX1[r,i] + udx1[r,c3] + udw1[r-1,c3] + (1-udx1[r,c0]) - 3, name = 'ufrMC0')
          # ------------------------------------------------------------------------------------------------------------
          self.model.addConstr(ufrMCa1[r,c0] <= udetXa1[r,c0], name = 'ufrmc0A')
          self.model.addConstr(ufrMCa1[r,c0] <= udx1[r,c3], name = 'ufrmc0A')
          self.model.addConstr(ufrMCa1[r,c0] <= udw1[r-1,c3], name = 'ufrmc0A')
          self.model.addConstr(ufrMCa1[r,c0] <= 1 - udx1[r,c0], name = 'ufrmc0A')
          self.model.addConstr(ufrMCa1[r,c0] >= udetXa1[r,i] + udx1[r,c3] + udw1[r-1,c3] + (1-udx1[r,c0]) - 3, name = 'ufrmc0A')

          # POSITION -- 1: DetX[c1] = DetW[c0] (No filter in this position)
          self.model.addConstr(uFrMC1[r,c1] == 0, 'ufrMC1')
          # ----------------------------------------------
          self.model.addConstr(ufrMCa1[r,c1] == 0, 'ufrmc1A')

          # POSITION -- 2: DetX[c2] = DetW[c1] & DetW[c2]
          self.model.addConstr(uFrMC1[r,c2] <= uDetX1[r,c2], name = 'ufrMC2')
          self.model.addConstr(uFrMC1[r,c2] <= udw1[r-1,c1], name = 'ufrMC2')
          self.model.addConstr(uFrMC1[r,c2] <= udw1[r-1,c2], name = 'ufrMC2')
          self.model.addConstr(uFrMC1[r,c2] <= 1 - udx1[r,c2], name = 'ufrMC2')
          self.model.addConstr(uFrMC1[r,c2] >= uDetX1[r,c2] + udw1[r-1,c1] + udw1[r-1,c2] + (1-udx1[r,c2]) - 3, name = 'ufrMC2')
          # ------------------------------------------------------------------------------------------------------------
          self.model.addConstr(ufrMCa1[r,c2] <= udetXa1[r,c2], name = 'ufrmc2A')
          self.model.addConstr(ufrMCa1[r,c2] <= udw1[r-1,c1], name = 'ufrmc2A')
          self.model.addConstr(ufrMCa1[r,c2] <= udw1[r-1,c2], name = 'ufrmc2A')
          self.model.addConstr(ufrMCa1[r,c2] <= 1 - udx1[r,c2], name = 'ufrmc2A')
          self.model.addConstr(ufrMCa1[r,c2] >= udetXa1[r,c2] + udw1[r-1,c1] + udw1[r-1,c2] + (1-udx1[r,c2]) - 3, name = 'ufrmc2A')

          # POSITION -- 3: DetX[c3] = DetW[c0] & DetW[c2]
          self.model.addConstr(uFrMC1[r,c3] <= uDetX1[r,c3], name = 'ufrMC3')
          self.model.addConstr(uFrMC1[r,c3] <= udw1[r-1,c0], name = 'ufrMC3')
          self.model.addConstr(uFrMC1[r,c3] <= udw1[r-1,c2], name = 'ufrMC3')
          self.model.addConstr(uFrMC1[r,c3] <= 1 - udx1[r,c3], name = 'ufrMC3')
          self.model.addConstr(uFrMC1[r,c3] >= uDetX1[r,c3] + udw1[r-1,c0] + udw1[r-1,c2] + (1-udx1[r,c2]) - 3, name = 'ufrMC3')
          # ------------------------------------------------------------------------------------------------------------
          self.model.addConstr(ufrMCa1[r,c3] <= udetXa1[r,c3], name = 'ufrmc3A')
          self.model.addConstr(ufrMCa1[r,c3] <= udw1[r-1,c0], name = 'ufrmc3A')
          self.model.addConstr(ufrMCa1[r,c3] <= udw1[r-1,c2], name = 'ufrmc3A')
          self.model.addConstr(ufrMCa1[r,c3] <= 1 - udx1[r,c3], name = 'ufrmc3A')
          self.model.addConstr(ufrMCa1[r,c3] >= udetXa1[r,c3] + udw1[r-1,c0] + udw1[r-1,c2] + (1-udx1[r,c2]) - 3, name = 'ufrmc3A')

    # Setting pre-guess all involved key in auxiliary upper trail (rf = friter)
    # -------------------------------------------------------------------------------------
    self.model.addConstr(sum(udw1[0,i] for i in range(16)) == 
                         (sum(ufrSCa1[r,i] for r in range(1, self.rEb+1) for i in range(16)) + 
                          sum(ufrMCa1[r,i] for r in range(2, self.rEb+1) for i in range(16))),
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


    ''' Some para. for clustering two trails '''
    pxt = self.model.addVars(16, vtype = GRB.BINARY, name = 'pxt')
    udy = self.model.addVars(16, vtype = GRB.BINARY, name = 'udy')
    uGstk = self.model.addVars(range(self.rEb), 16, vtype = GRB.BINARY, name = 'uGstk')
    ugstka = self.model.addVars(range(self.rEb), 16, vtype = GRB.BINARY, name = 'ugstka')
    for i in range(16):
      self.model.addConstr(pxt[i] >= udw0[0,i])
      self.model.addConstr(pxt[i] <= udw1[0,i])
      self.model.addConstr(pxt[i] <= udw0[0,i] + udw1[0,i])
    for i in range(16):
      self.model.addConstr(udy[i] >= udy0[self.rEb,i])
      self.model.addConstr(udy[i] >= udy1[self.rEb,i])
      self.model.addConstr(udy[i] <= udy0[self.rEb,i] + udy1[self.rEb,i])
    for r in range(self.rEb):
      for i in range(16):
        self.model.addConstr(uGstk[r,i] >= uGstk0[r,i])
        self.model.addConstr(uGstk[r,i] >= uGstk1[r,i])
        self.model.addConstr(uGstk[r,i] <= uGstk0[r,i] + uGstk1[r,i])
        self.model.addConstr(ugstka[r,i] >= ugstka0[r,i])
        self.model.addConstr(ugstka[r,i] >= ugstka1[r,i])
        self.model.addConstr(uGstk[r,i] <= ugstka0[r,i] + ugstka1[r,i])

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
   
    # ******************************************************
    ''' Contradiction '''
    duFX0 = self.model.addVars(range(self.rEb, self.rEb + self.rDis), 16, vtype = GRB.BINARY, name = 'duFX0')
    duFX1 = self.model.addVars(range(self.rEb, self.rEb + self.rDis), 16, vtype = GRB.BINARY, name = 'duFX1')
    duTX0 = self.model.addVars(range(self.rEb, self.rEb + self.rDis), 16, vtype = GRB.BINARY, name = 'duTX0')
    duTX1 = self.model.addVars(range(self.rEb, self.rEb + self.rDis), 16, vtype = GRB.BINARY, name = 'duTX1')
    duAZ0 = self.model.addVars(range(self.rEb, self.rEb + self.rDis), 16, vtype = GRB.BINARY, name = 'duAZ0')
    duAZ1 = self.model.addVars(range(self.rEb, self.rEb + self.rDis), 16, vtype = GRB.BINARY, name = 'duAZ1')
    duAW0 = self.model.addVars(range(self.rEb, self.rEb + self.rDis), 16, vtype = GRB.BINARY, name = 'duAW0')
    duAW1 = self.model.addVars(range(self.rEb, self.rEb + self.rDis), 16, vtype = GRB.BINARY, name = 'duAW1')
    dlAX = self.model.addVars(range(self.rEb+1, self.rEb + self.rDis+1), 16, vtype = GRB.BINARY, name = 'dlAX')
    dlFZ = self.model.addVars(range(self.rEb, self.rEb + self.rDis), 16, vtype = GRB.BINARY, name = 'dlFZ')
    dlFW = self.model.addVars(range(self.rEb, self.rEb + self.rDis), 16, vtype = GRB.BINARY, name = 'dlFW')
    # power-reduced aid
    PRA0 = self.model.addVars(range(self.rEb, self.rEb + self.rDis), 16, vtype = GRB.BINARY, name = 'PRA0')
    PRA1 = self.model.addVars(range(self.rEb, self.rEb + self.rDis), 16, vtype = GRB.BINARY, name = 'PRA1')
    PRA2 = self.model.addVars(range(self.rEb, self.rEb + self.rDis), 16, vtype = GRB.BINARY, name = 'PRA2')
    PRA3 = self.model.addVars(range(self.rEb, self.rEb + self.rDis), 16, vtype = GRB.BINARY, name = 'PRA3')
    PRA4 = self.model.addVars(range(self.rEb, self.rEb + self.rDis), 16, vtype = GRB.BINARY, name = 'PRA4')

    for r in range(self.rEb, self.rEb + self.rDis):
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
        self.model.addConstr(dlAX[r+1,i] == 1 - lDX[r+1,i], name = '2lAX')
        self.model.addConstr(dlFZ[r,i] == lDZ[r,i] - ldz[r,i], name = '2lFZ')
        self.model.addConstr(dlFW[r,i] == lDW[r,i] - ldw[r,i], name = '2lFW')

    for r in range(self.rEb+1, self.rEb + self.rDis-1):
      for i in range(16):
        # Power-reduce
        # constraint-1 [(𝒖𝑫𝑿_𝟎−𝒖𝒅𝒙_𝟎 )+(𝒖𝑫𝑿_𝟏−𝒖𝒅𝒙_𝟏 )]×(𝟏−𝒖𝒅𝒙_𝟎 )×(𝟏−𝒖𝒅𝒙_𝟏 )×(𝟏−𝒍𝑫𝑿)
        self.model.addGenConstrAnd(PRA0[r,i], [duTX0[r,i], duTX1[r,i]], name = '2PDA0')
        self.model.addGenConstrAnd(PRA1[r,i], [PRA0[r,i], dlAX[r,i]], name = '2PDA1')

        # constraint-2 (𝒍𝑫𝒁−𝒍𝒅𝒛)×(𝟏−𝒖𝑫𝒁_𝟎 )×(𝟏−𝒖𝑫𝒁_𝟏)
        self.model.addGenConstrAnd(PRA2[r, i], [duAZ0[r,i], duAZ1[r,i]],name = '2PDA2')

        # constraint-3 (𝒍𝑫𝑾−𝒍𝒅𝒘)×(𝟏−𝒖𝑫𝑾_𝟎 )×(𝟏−𝒖𝑫𝑾_𝟏)
        self.model.addGenConstrAnd(PRA3[r, i], [duAW0[r,i], duAW1[r,i]],name = '2PDA3')

        # constraint-4 (𝒖𝑫𝑿_𝟎−𝒖𝒅𝒙_𝟎 )∗(𝒖𝑫𝑿_𝟏−𝒖𝒅𝒙_𝟏 )∗(𝒍𝑫𝒁−𝒍𝒅𝒛)
        self.model.addGenConstrAnd(PRA4[r, i], [duFX0[r,i], duFX1[r,i]],name = '2PDA4')

    self.model.addConstr(sum(
    (
      ((duFX0[r,i] + duFX1[r,i]) * PRA1[r,i]) + 
      (dlFZ[r,i] * PRA2[r,i]) + 
      (dlFW[r,i] * PRA3[r,i]) + 
      (dlFZ[r,i] * PRA4[r,i])
    ) 
    for i in range(16) for r in range(self.rEb, self.rEb + self.rDis)) >= 1, name = 'MainContra')

    # self.model.addConstr(sum((uDX[r,i] - udx[r,i]) * (lDY[r,i] - ldy[r,i]) 
    #                          for i in range(16) for r in range(self.rEb, self.rEb+self.rDis)) >= 1, name = 'iBCT')
    # ******************************************************


    '''
    ==========================================================================================
         (Begin) - Complexity
    '''

    rb0  = self.model.addVar(vtype = GRB.INTEGER, name = 'rb0')
    rb1  = self.model.addVar(vtype = GRB.INTEGER, name = 'rb1')
    cbp0 = self.model.addVar(vtype = GRB.INTEGER, name = 'cbp0')
    cbp1 = self.model.addVar(vtype = GRB.INTEGER, name = 'cbp1')

    '''
    Parameters
    '''
    rb  = self.model.addVar(vtype = GRB.INTEGER, name = 'rb')
    rf  = self.model.addVar(vtype = GRB.INTEGER, name = 'rf')
    cb  = self.model.addVar(vtype = GRB.INTEGER, name = 'cb')
    cf  = self.model.addVar(vtype = GRB.INTEGER, name = 'cf')
    iks = self.model.addVar(vtype = GRB.INTEGER, name = 'iks')
    gks = self.model.addVar(vtype = GRB.INTEGER, name = 'gks')
    # cbp = self.model.addVar(vtype = GRB.INTEGER, name = 'cbp')
    cfp = self.model.addVar(vtype = GRB.INTEGER, name = 'cfp')

    # mb  + mf  = ikSUM (in key bridge)
    self.model.addConstr(iks == self.c_size * sum(ikSUM[i] for i in range(16)), name = 'IKS')
    # mbp + mfp = gkSUM (in key bridge)
    self.model.addConstr(gks == self.c_size * sum(gkSUM[i] for i in range(16)), name = 'GKS')

    self.model.addConstr(rb0 == self.c_size * sum(udw0[0,i] for i in range(16)), name = 'rb0')
    self.model.addConstr(rb1 == self.c_size * sum(udw1[0,i] for i in range(16)), name = 'rb1')
    self.model.addConstr(rb == self.c_size * sum(pxt[i] for i in range(16)), name = 'rb')

    self.model.addConstr(rf == self.c_size * sum(ldz[totR-1,i] for i in range(16)), name = 'rf')

    self.model.addConstr(cb == rb - self.c_size * sum(udy[i] for i in range(16)), name = 'cb')
    self.model.addConstr(cf == rf - self.c_size * sum(ldx[totR-self.rEf,i] for i in range(16)), name = 'cf')

    self.model.addConstr(cbp0 == self.c_size * 
                         (sum(uFrSC0[r,i] for r in range(1, self.rEb+1) for i in range(16)) + 
                          sum(uFrMC0[r,i] for r in range(2, self.rEb+1) for i in range(16))), name = 'cbp0')
    self.model.addConstr(cbp1 == self.c_size * 
                         (sum(uFrSC1[r,i] for r in range(1, self.rEb+1) for i in range(16)) + 
                          sum(uFrMC1[r,i] for r in range(2, self.rEb+1) for i in range(16))), name = 'cbp0')

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
    self.model.addConstr(Mc >= Qc - (cbp0 + cbp1 + 2*cfp), name = 'Mc2')
    self.model.addConstr(Mc >= iks, name = 'Mc3')

    # Time
    # T0 (T0 = Dc)
    self.model.addConstr(T0 == Dc, name = 'T0')
    # T1 (T1 = 2^{mb'+mf'} · D · 4)
    self.model.addConstr(T1 == gks + Dc + 2, name = 'T1')
    # T2 (T2 = 2^{mb'+mf'} · D · min[ 2^{rb-cb'}, D · 2^{rf-cf'-n} ] · 2)
    if self.cP == True:
      self.model.addConstr(T2 >= gks + Dc + rb0 - cbp0, name = 'T2p0')
      self.model.addConstr(T2 >= gks + Dc + rb1 - cbp1, name = 'T2p1')
    else:
      self.model.addConstr(T2 == gks + Dc + (Dc + rf - cfp - self.b_size) + 1, name = 'T2c')
    # T31 (T31 = 2^{mb'+mf'-2cb'-2cf'} · Q)
    self.model.addConstr(T31 == (gks - cbp0 - cbp1 - 2*cfp) + Qc, name = 'T31')
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
