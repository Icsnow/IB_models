

% Test1: 2+7+2, x=20

array[0..15] of var 0..15: SRp =  array1d(0..15, [0, 5, 10, 15, 4, 9, 14, 3, 8, 13, 2, 7, 12, 1, 6, 11]);
array[0 .. 15, 0 .. 15] of 0 .. 15: hTable = array2d(0..15, 0..15, 
  [ 0,  7, 14,  9,  8, 15,  6,  1,  0,  7, 14,  9,  8, 15,  6,  1,
    1,  0,  7, 14,  9,  8, 15,  6,  1,  0,  7, 14,  9,  8, 15,  6,
    2, 13, 12,  3, 10,  5,  4, 11,  2, 13, 12,  3, 10,  5,  4, 11,
    3, 10,  5,  4, 11,  2, 13, 12,  3, 10,  5,  4, 11,  2, 13, 12,
    4, 11,  2, 13, 12,  3, 10,  5,  4, 11,  2, 13, 12,  3, 10,  5,
    5,  4, 11,  2, 13, 12,  3, 10,  5,  4, 11,  2, 13, 12,  3, 10,
    6,  1,  0,  7, 14,  9,  8, 15,  6,  1,  0,  7, 14,  9,  8, 15,
    7, 14,  9,  8, 15,  6,  1,  0,  7, 14,  9,  8, 15,  6,  1,  0,
    8, 15,  6,  1,  0,  7, 14,  9,  8, 15,  6,  1,  0,  7, 14,  9,
    9,  8, 15,  6,  1,  0,  7, 14,  9,  8, 15,  6,  1,  0,  7, 14,
   10,  5,  4, 11,  2, 13, 12,  3, 10,  5,  4, 11,  2, 13, 12,  3, 
   11,  2, 13, 12,  3, 10,  5,  4, 11,  2, 13, 12,  3, 10,  5,  4,
   12,  3, 10,  5,  4, 11,  2, 13, 12,  3, 10,  5,  4, 11,  2, 13,
   13, 12,  3, 10,  5,  4, 11,  2, 13, 12,  3, 10,  5,  4, 11,  2,
   14,  9,  8, 15,  6,  1,  0,  7, 14,  9,  8, 15,  6,  1,  0,  7,
   15,  6,  1,  0,  7, 14,  9,  8, 15,  6,  1,  0,  7, 14,  9,  8]
);
int: tk_size = 2; % 2 => 256
int: c_size = 8;
int: b_size = c_size * 16;
int: rEb;     
int: rDis;
int: rEf;

% #=========================================================================#
% /============================== Parameters ===============================\
% #=========================================================================#
array[-1..(rEb+rDis-1), 0..15] of var 0 .. 3: uW0;
array[ 0..(rEb+rDis), 0..15] of var 0 .. 3: uX0;
array[ 0..(rEb+rDis-1), 0..15] of var 0 .. 3: uY0;
array[ 0..(rEb+rDis-1), 0..15] of var 0 .. 3: uZ0;
array[ 0..(rEb+rDis), 0..15] of var 0 .. 1: ustk0;
array[ 0..(rEb+rDis), 0..15] of var 0 .. 1: ucan0;
array[ 0..(rEb+rDis-1), 0.. 3] of var 0 .. 3: uAC0; % mark the state of one coulmn

array[-1..(rEb+rDis-1), 0..15] of var 0 .. 3: uW1;
array[ 0..(rEb+rDis), 0..15] of var 0 .. 3: uX1;
array[ 0..(rEb+rDis-1), 0..15] of var 0 .. 3: uY1;
array[ 0..(rEb+rDis-1), 0..15] of var 0 .. 3: uZ1;
array[ 0..(rEb+rDis), 0..15] of var 0 .. 1: ustk1;
array[ 0..(rEb+rDis), 0..15] of var 0 .. 1: ucan1;
array[ 0..(rEb+rDis-1), 0.. 3] of var 0 .. 3: uAC1; % mark the state of one coulmn

% --==================================================================--
% ||========================= Upper Trails ===========================||
% --==================================================================--
% ATK (Add round TweaKey)
constraint 
  forall(r in -1 .. (rEb-1), i in 0..15)(
    bATK(uX0[r+1,i], ustk0[r+1,i], uW0[r,i], ucan0[r+1,i]) /\
    bATK(uX1[r+1,i], ustk1[r+1,i], uW1[r,i], ucan1[r+1,i])
  ) /\ 
  forall(r in rEb .. (rEb+rDis-1), i in 0..15)(
    fATK(uW0[r,i], ustk0[r+1,i], uX0[r+1,i], ucan0[r+1,i]) /\
    fATK(uW1[r,i], ustk1[r+1,i], uX1[r+1,i], ucan1[r+1,i])
  )
;

% SB (SBox substitution)
constraint 
  forall(r in 0 .. (rEb-1), i in 0..15)(
    SBox(uY0[r,i], uX0[r,i]) /\
    SBox(uY1[r,i], uX1[r,i])
  ) /\
  forall(r in rEb .. (rEb+rDis-1), i in 0..15)(
    SBox(uX0[r,i], uY0[r,i]) /\
    SBox(uX1[r,i], uY1[r,i])
  )
;

% SR (ShiftRow permutation)
constraint forall(r in 0 .. (rEb+rDis-1), i in 0 .. 15)(
  uY0[r,SRp[i]] = uZ0[r,i] /\
  uY1[r,SRp[i]] = uZ1[r,i]
);

% MC (Mixcolumn operation)
constraint
  forall(r in 0 .. (rEb-1), c in 0..3)(
    MC(uW0[r,4*c], uW0[r,4*c+1], uW0[r,4*c+2], uW0[r,4*c+3], 
       uAC0[r,c],
       uZ0[r,4*c], uZ0[r,4*c+1], uZ0[r,4*c+2], uZ0[r,4*c+3])
    /\
    MC(uW1[r,4*c], uW1[r,4*c+1], uW1[r,4*c+2], uW1[r,4*c+3], 
       uAC1[r,c],
       uZ1[r,4*c], uZ1[r,4*c+1], uZ1[r,4*c+2], uZ1[r,4*c+3])
  )
  /\
  forall(r in rEb .. (rEb+rDis-1), c in 0..3)(
    MC(uZ0[r,4*c], uZ0[r,4*c+1], uZ0[r,4*c+2], uZ0[r,4*c+3],
       uAC0[r,c],
       uW0[r,4*c], uW0[r,4*c+1], uW0[r,4*c+2], uW0[r,4*c+3])
    /\
    MC(uZ1[r,4*c], uZ1[r,4*c+1], uZ1[r,4*c+2], uZ1[r,4*c+3],
       uAC1[r,c],
       uW1[r,4*c], uW1[r,4*c+1], uW1[r,4*c+2], uW1[r,4*c+3])
  )
;


% --=====================================================--
% ||================= Upper Key Schedule ================||
% --=====================================================--
array[0..15] of var 0 .. 1: uLANE0;
array[0..15] of var int: uT1Can0;
array[0..(rEb+rDis-1), 0..3] of var 0..4: uT2a0;
array[0..(rEb+rDis-1), 0..3] of var 0..4: uT2b0;
array[0..(rEb+rDis-1), 0..3] of var 0..4: uT2c0;
array[0..(rEb+rDis-1), 0..3] of var int: uT2Can0;

array[0..15] of var 0 .. 1: uLANE1;
array[0..15] of var int: uT1Can1;
array[0..(rEb+rDis-1), 0..3] of var 0..4: uT2a1;
array[0..(rEb+rDis-1), 0..3] of var 0..4: uT2b1;
array[0..(rEb+rDis-1), 0..3] of var 0..4: uT2c1;
array[0..(rEb+rDis-1), 0..3] of var int: uT2Can1;

% Set Active
constraint (sum(i in 0..15)(uLANE0[i]) >= 1) /\ (sum(i in 0..15)(uLANE1[i]) >= 1);

% Mark LANE
constraint forall(i in 0..15, r in 0..(rEb+rDis))(
  (uLANE0[i] >= ustk0[r, hTable[i,r]]) /\
  (uLANE1[i] >= ustk1[r, hTable[i,r]])
);

% Type 1 Cancellation
constraint forall(i in 0..15)(
  (uT1Can0[i] = (rEb+rDis+1) * uLANE0[i] - sum(r in 0..(rEb+rDis))(ustk0[r, hTable[i,r]])) /\
  (uT1Can0[i] <= tk_size - 1) /\
  (uT1Can1[i] = (rEb+rDis+1) * uLANE1[i] - sum(r in 0..(rEb+rDis))(ustk1[r, hTable[i,r]])) /\
  (uT1Can1[i] <= tk_size - 1)
);

% Type 2 Cancellation (Auxiliary para "a")
constraint forall(r in  0..(rEb+rDis-1), c in 0..3)(
  uT2a0[r,c] = sum(i in 0..3)(if (uZ0[r,4*c+i] >= 1) then 1 else 0 endif) /\
  uT2a1[r,c] = sum(i in 0..3)(if (uZ1[r,4*c+i] >= 1) then 1 else 0 endif)
);

% Type 2 Cancellation (Auxiliary para "b")
constraint forall(r in  0..(rEb+rDis-1), c in 0..3)(
  uT2b0[r,c] = if (uAC0[r,c] >= 1) then (4 - sum(i in 0..3)(if (uW0[r,4*c+i] >= 1) then 1 else 0 endif)) else 0 endif /\
  uT2b1[r,c] = if (uAC1[r,c] >= 1) then (4 - sum(i in 0..3)(if (uW1[r,4*c+i] >= 1) then 1 else 0 endif)) else 0 endif
);

% Type 2 Cancellation (Auxiliary para "c")
constraint forall(r in  0..(rEb+rDis-1), c in 0..3)(
  uT2c0[r,c] = sum(i in 0..3)(ucan0[r+1,4*c+i]) /\
  uT2c1[r,c] = sum(i in 0..3)(ucan1[r+1,4*c+i])
);

% Type 2 Cancellation (main)
constraint forall(r in 0..(rEb+rDis-1), c in 0..3)(
  uT2Can0[r,c] = max(0, (uT2b0[r,c] + uT2c0[r,c] - uT2a0[r,c])) /\
  uT2Can1[r,c] = max(0, (uT2b1[r,c] + uT2c1[r,c] - uT2a1[r,c]))
);

% Counting of all Cancellations
constraint 
  (tk_size * sum(i in 0..15)(uLANE0[i]) - sum(i in 0..15)(uT1Can0[i]) - sum(r in 0..(rEb+rDis-1), c in 0..3)(uT2Can0[r,c])) >= 1 /\
  (tk_size * sum(i in 0..15)(uLANE1[i]) - sum(i in 0..15)(uT1Can1[i]) - sum(r in 0..(rEb+rDis-1), c in 0..3)(uT2Can1[r,c])) >= 1
;


% --=========================================================--
% ||============== Upper Guess-and-Determine ================||
% --=========================================================--
array[ 0..(rEb-1), 0..15] of var 0 .. 1: uGstk0;
array[-1..(rEb-1), 0..15] of var 0 .. 1: udetW0;
array[ 0..(rEb-1), 0..15] of var 0 .. 1: udetX0;
array[ 0..(rEb-1), 0..15] of var 0 .. 1: udetY0;
array[ 0..(rEb-1), 0..15] of var 0 .. 1: udetZ0;
array[ 0..(rEb-1), 0..15] of var 0 .. 1: uFrSB0;
array[ 0..(rEb-1), 0..15] of var 0 .. 1: uFrMC0;
constraint forall(i in 0..15)(udetW0[-1,i] = 1);

array[ 0..(rEb-1), 0..15] of var 0 .. 1: uGstk1;
array[-1..(rEb-1), 0..15] of var 0 .. 1: udetW1;
array[ 0..(rEb-1), 0..15] of var 0 .. 1: udetX1;
array[ 0..(rEb-1), 0..15] of var 0 .. 1: udetY1;
array[ 0..(rEb-1), 0..15] of var 0 .. 1: udetZ1;
array[ 0..(rEb-1), 0..15] of var 0 .. 1: uFrSB1;
array[ 0..(rEb-1), 0..15] of var 0 .. 1: uFrMC1;
constraint forall(i in 0..15)(udetW1[-1,i] = 1);

% Guess and determine
constraint forall(r in -1..(rEb-2), i in 0..15)(
  if (udetW0[r,i] == 1) /\ (uGstk0[r+1,i] == 1) then (udetX0[r+1,i] = 1)
  else (udetX0[r+1,i] = 0)
  endif
  /\
  if (udetW1[r,i] == 1) /\ (uGstk1[r+1,i] == 1) then (udetX1[r+1,i] = 1)
  else (udetX1[r+1,i] = 0)
  endif
);

% "Det" via SB and SR
constraint forall(r in 0..(rEb-1), i in 0..15)(
  udetY0[r,i] = udetX0[r,i] /\ udetZ0[r,i] = udetY0[r,SRp[i]]
  /\
  udetY1[r,i] = udetX1[r,i] /\ udetZ1[r,i] = udetY1[r,SRp[i]]
);

% "Det" via MC
constraint forall(r in 0..(rEb-1), c in 0..3)(
  if (sum(i in 0..3)(udetZ0[r,4*c+i]) == 4) then (sum(i in 0..3)(udetW0[r,4*c+i]) == 4) 
  else (sum(i in 0..3)(udetW0[r,4*c+i]) == 0)
  endif
  /\
  if (sum(i in 0..3)(udetZ1[r,4*c+i]) == 4) then (sum(i in 0..3)(udetW1[r,4*c+i]) == 4) 
  else (sum(i in 0..3)(udetW1[r,4*c+i]) == 0)
  endif
);

% Filters from SB
constraint forall(r in 0..(rEb-1), i in 0..15)(
  if (udetY0[r,i] == 1) /\ (uY0[r,i] == 1) then (uFrSB0[r,i] = 1)
  else (uFrSB0[r,i] = 0)
  endif
  /\
  if (udetY1[r,i] == 1) /\ (uY1[r,i] == 1) then (uFrSB1[r,i] = 1)
  else (uFrSB1[r,i] = 0)
  endif
);

% Filter from MC
constraint forall(r in 0..(rEb-1), c in 0..3)(
  if exists(i in 0..3)(uW0[r,4*c+i] >= 2)
  then forall(i in 0..3) (if (udetW0[r,4*c+i] == 1) /\ (uW0[r,4*c+i] <= 1) then (uFrMC0[r,4*c+i] == 1) else (uFrMC0[r,4*c+i] == 0) endif)
  else sum(i in 0..3)(uFrMC0[r,4*c+i]) = 0
  endif
  /\
  if exists(i in 0..3)(uW1[r,4*c+i] >= 2) 
  then forall(i in 0..3) (if (udetW1[r,4*c+i] == 1) /\ (uW1[r,4*c+i] <= 1) then (uFrMC1[r,4*c+i] == 1) else (uFrMC1[r,4*c+i] == 0) endif)
  else sum(i in 0..3)(uFrMC1[r,4*c+i]) = 0
  endif
);



% %%%%%%%%%%%%%%%%%%% TEST UPPER %%%%%%%%%%%%%%%%%%%%%%
% constraint forall(i in [0,1,2,6,8,10,11,13,15])(uLANE0[i] = 1 /\ uLANE1[i] = 1) /\ sum(i in 0..15)(uLANE0[i]) = 9 /\ sum(i in 0..15)(uLANE1[i]) = 9;
% constraint sum(i in 0..15)(ustk0[3,i]) = 0 /\ sum(i in 0..15)(ustk1[3,i]) = 0;
% constraint forall(r in 2..3, i in 0..15)(uW0[r,i] = 0 /\ uW1[r,i] = 0);
% constraint forall(i in [1,2,4,5,8,9,10,11,14,15])(uZ0[1,i] = 0 /\ uZ1[1,i] = 0) /\ sum(i in 0..15)(uZ0[1,i]) = 6 /\ sum(i in 0..15)(uZ1[1,i]) = 6;
% constraint forall(i in [4,6,8,9,10,11])(uZ0[0,i] = 0 /\ uZ1[0,i] = 0) /\ sum(i in 0..15)(uZ0[0,i]) >= 10 /\ sum(i in 0..15)(uZ1[0,i]) >= 10;

% constraint forall(i in [3,9])(uGstk0[0,i] = 1 /\ uGstk1[0,i] = 1) /\ sum(r in 0..rEb-1, i in 0..15)(uGstk0[r,i]) = 2 /\ sum(r in 0..rEb-1, i in 0..15)(uGstk1[r,i]) = 2;
% %%%%%%%%%%%%%%%%%%% TEST UPPER %%%%%%%%%%%%%%%%%%%%%%

% --==================================================================--
% ||========================= Lower Trails ===========================||
% --==================================================================--

array[rEb..(rEb+rDis+rEf-1), 0..15] of var 0..3: lW0;
array[rEb..(rEb+rDis+rEf), 0..15] of var 0..3: lX0;
array[rEb..(rEb+rDis+rEf-1), 0..15] of var 0..3: lY0;
array[rEb..(rEb+rDis+rEf-1), 0..15] of var 0..3: lZ0;
array[(rEb+1)..(rEb+rDis+rEf), 0..15] of var 0..1: lstk0;
array[(rEb+1)..(rEb+rDis+rEf), 0..15] of var 0..1: lcan0;
array[rEb..(rEb+rDis+rEf-1), 0.. 3] of var 0..3: lAC0; % mark the state of one coulmn

array[rEb..(rEb+rDis+rEf-1), 0..15] of var 0..3: lW1;
array[rEb..(rEb+rDis+rEf), 0..15] of var 0..3: lX1;
array[rEb..(rEb+rDis+rEf-1), 0..15] of var 0..3: lY1;
array[rEb..(rEb+rDis+rEf-1), 0..15] of var 0..3: lZ1;
array[(rEb+1)..(rEb+rDis+rEf), 0..15] of var 0..1: lstk1;
array[(rEb+1)..(rEb+rDis+rEf), 0..15] of var 0..1: lcan1;
array[rEb..(rEb+rDis+rEf-1), 0.. 3] of var 0..3: lAC1; % mark the state of one coulmn

% SB (SBox substitution)
constraint 
  forall(r in rEb..(rEb+rDis-1), i in 0..15)(
    SBox(lY0[r,i], lX0[r,i]) /\
    SBox(lY1[r,i], lX1[r,i])
  ) /\
  forall(r in (rEb+rDis)..(rEb+rDis+rEf-1), i in 0..15)(
    SBox(lX0[r,i], lY0[r,i]) /\
    SBox(lX1[r,i], lY1[r,i])
  )
;

% SR (ShiftRow permutation)
constraint forall(r in rEb..(rEb+rDis+rEf-1), i in 0 .. 15)(
  lY0[r,SRp[i]] = lZ0[r,i] /\
  lY1[r,SRp[i]] = lZ1[r,i]
);

% MC (Mixcolumn operation)
constraint
  forall(r in rEb..(rEb+rDis-1), c in 0..3)(
    MC(lW0[r,4*c], lW0[r,4*c+1], lW0[r,4*c+2], lW0[r,4*c+3], 
       lAC0[r,c],
       lZ0[r,4*c], lZ0[r,4*c+1], lZ0[r,4*c+2], lZ0[r,4*c+3])
    /\
    MC(lW1[r,4*c], lW1[r,4*c+1], lW1[r,4*c+2], lW1[r,4*c+3], 
       lAC1[r,c],
       lZ1[r,4*c], lZ1[r,4*c+1], lZ1[r,4*c+2], lZ1[r,4*c+3])
  )
  /\
  forall(r in (rEb+rDis)..(rEb+rDis+rEf-1), c in 0..3)(
    MC(lZ0[r,4*c], lZ0[r,4*c+1], lZ0[r,4*c+2], lZ0[r,4*c+3],
       lAC0[r,c],
       lW0[r,4*c], lW0[r,4*c+1], lW0[r,4*c+2], lW0[r,4*c+3])
    /\
    MC(lZ1[r,4*c], lZ1[r,4*c+1], lZ1[r,4*c+2], lZ1[r,4*c+3],
       lAC1[r,c],
       lW1[r,4*c], lW1[r,4*c+1], lW1[r,4*c+2], lW1[r,4*c+3])
  )
;

% ATK (Add round TweaKey)
constraint 
  forall(r in rEb..(rEb+rDis-2), i in 0..15)(
    bATK(lX0[r+1,i], lstk0[r+1,i], lW0[r,i], lcan0[r+1,i]) /\
    bATK(lX1[r+1,i], lstk1[r+1,i], lW1[r,i], lcan1[r+1,i])
  ) /\ 
  forall(r in (rEb+rDis-1)..(rEb+rDis+rEf-1), i in 0..15)(
    fATK(lW0[r,i], lstk0[r+1,i], lX0[r+1,i], lcan0[r+1,i]) /\
    fATK(lW1[r,i], lstk1[r+1,i], lX1[r+1,i], lcan1[r+1,i])
  )
;


% --------------- Transform ENC to eq-ENCv --------------------
array[(rEb+rDis)..(rEb+rDis+rEf-2),   0..15] of var 0..3: leqW0;
array[(rEb+rDis)..(rEb+rDis+rEf-1), 0..15] of var 0..3: leqX0;
array[(rEb+rDis)..(rEb+rDis+rEf-1), 0..15] of var 0..3: leqY0;
array[(rEb+rDis)..(rEb+rDis+rEf-1),   0..15] of var 0..3: leqZ0;
array[(rEb+rDis+1)..(rEb+rDis+rEf),   0..15] of var 0..3: leqk0;

array[(rEb+rDis)..(rEb+rDis+rEf-2),   0..15] of var 0..3: leqW1;
array[(rEb+rDis)..(rEb+rDis+rEf-1), 0..15] of var 0..3: leqX1;
array[(rEb+rDis)..(rEb+rDis+rEf-1), 0..15] of var 0..3: leqY1;
array[(rEb+rDis)..(rEb+rDis+rEf-1),   0..15] of var 0..3: leqZ1;
array[(rEb+rDis+1)..(rEb+rDis+rEf),   0..15] of var 0..3: leqk1;

% eqk = SRinv(MCinv(stk))
constraint forall(r in (rEb+rDis+1)..(rEb+rDis+rEf), c in 0..3)(
  if exists(i in 0..3)(lstk0[r,4*c+i] >= 1)
  then (sum(i in 0..3)(leqk0[r,SRp[4*c+i]]) + sum(i in 0..3)(lstk0[r,4*c+i]) >= 5) /\ (sum(i in 0..3)(leqk0[r,SRp[4*c+i]]) + sum(i in 0..3)(lstk0[r,4*c+i]) <= 8)
  else (sum(i in 0..3)(leqk0[r,SRp[4*c+i]]) = 0)
  endif
  /\
  if exists(i in 0..3)(lstk1[r,4*c+i] >= 1)
  then (sum(i in 0..3)(leqk1[r,SRp[4*c+i]]) + sum(i in 0..3)(lstk1[r,4*c+i]) >= 5) /\ (sum(i in 0..3)(leqk1[r,SRp[4*c+i]]) + sum(i in 0..3)(lstk1[r,4*c+i]) <= 8)
  else (sum(i in 0..3)(leqk1[r,SRp[4*c+i]]) = 0)
  endif
);

% To eq-State (X = eqX)
constraint forall(r in (rEb+rDis)..(rEb+rDis+rEf-1), i in 0..15)(
  if (lX0[r,i] <= 1) then (leqX0[r,i] = lX0[r,i]) else (leqX0[r,i] >= lX0[r,i]) endif
   /\
  if (lX1[r,i] <= 1) then (leqX1[r,i] = lX1[r,i]) else (leqX1[r,i] >= lX1[r,i]) endif
);

% eq-SBox
constraint forall(r in (rEb+rDis)..(rEb+rDis+rEf-1), i in 0..15)(
  SBox(leqX0[r,i],leqY0[r,i]) /\
  SBox(leqX1[r,i],leqY1[r,i])
);

% eq-ATK
constraint 
  let{
      array[(rEb+rDis)..(rEb+rDis+rEf-1), 0..15] of var 0..1: leqCan0;
      array[(rEb+rDis)..(rEb+rDis+rEf-1), 0..15] of var 0..1: leqCan1;
  } in 
      forall(r in (rEb+rDis)..(rEb+rDis+rEf-1), i in 0..15)((leqCan0[r,i] = 0) /\ fATK(leqY0[r,i], leqk0[r+1,i], leqZ0[r,i], leqCan0[r,i])) /\
      forall(r in (rEb+rDis)..(rEb+rDis+rEf-1), i in 0..15)((leqCan1[r,i] = 0) /\ fATK(leqY1[r,i], leqk1[r+1,i], leqZ1[r,i], leqCan1[r,i]))
;

% eq-SR
constraint forall(r in (rEb+rDis)..(rEb+rDis+rEf-2), i in 0..15)(
  leqZ0[r,SRp[i]] = leqW0[r,i] /\
  leqZ1[r,SRp[i]] = leqW1[r,i]
);
  
% eq-MC
constraint 
  let{
      array[(rEb+rDis)..(rEb+rDis+rEf-2), 0.. 3] of var 0..3: leqAC0;
      array[(rEb+rDis)..(rEb+rDis+rEf-2), 0.. 3] of var 0..3: leqAC1;

  } in
  forall(r in (rEb+rDis)..(rEb+rDis+rEf-2), c in 0..3)(
    MC(leqW0[r,4*c],   leqW0[r,4*c+1],   leqW0[r,4*c+2],   leqW0[r,4*c+3],
       leqAC0[r,c],
       leqX0[r+1,4*c], leqX0[r+1,4*c+1], leqX0[r+1,4*c+2], leqX0[r+1,4*c+3])
    /\
    MC(leqW1[r,4*c],   leqW1[r,4*c+1],   leqW1[r,4*c+2],   leqW1[r,4*c+3],
       leqAC1[r,c],
       leqX1[r+1,4*c], leqX1[r+1,4*c+1], leqX1[r+1,4*c+2], leqX1[r+1,4*c+3])
  )
;

% --=====================================================--
% ||================= Lower Key Schedule ================||
% --=====================================================--
array[0..15] of var 0 .. 1: lLANE0;
array[0..15] of var int: lT1Can0;
array[rEb..(rEb+rDis+rEf-1), 0..3] of var 0..4: lT2a0;
array[rEb..(rEb+rDis+rEf-1), 0..3] of var 0..4: lT2b0;
array[rEb..(rEb+rDis+rEf-1), 0..3] of var 0..4: lT2c0;
array[rEb..(rEb+rDis+rEf-1), 0..3] of var int: lT2Can0;

array[0..15] of var 0 .. 1: lLANE1;
array[0..15] of var int: lT1Can1;
array[rEb..(rEb+rDis+rEf-1), 0..3] of var 0..4: lT2a1;
array[rEb..(rEb+rDis+rEf-1), 0..3] of var 0..4: lT2b1;
array[rEb..(rEb+rDis+rEf-1), 0..3] of var 0..4: lT2c1;
array[rEb..(rEb+rDis+rEf-1), 0..3] of var int: lT2Can1;

% Set Active
constraint (sum(i in 0..15)(lLANE0[i]) >= 1) /\ (sum(i in 0..15)(lLANE1[i]) >= 1);

% Mark LANE
constraint forall(i in 0..15, r in (rEb+1)..(rEb+rDis+rEf))(
  (lLANE0[i] >= lstk0[r, hTable[i,r]]) /\
  (lLANE1[i] >= lstk1[r, hTable[i,r]])
);

% Type 1 Cancellation
constraint forall(i in 0..15)(
  (lT1Can0[i] = (rDis+rEf) * lLANE0[i] - sum(r in (rEb+1)..(rEb+rDis+rEf))(lstk0[r, hTable[i,r]])) /\
  (lT1Can0[i] <= tk_size - 1) /\
  (lT1Can1[i] = (rDis+rEf) * lLANE1[i] - sum(r in (rEb+1)..(rEb+rDis+rEf))(lstk1[r, hTable[i,r]])) /\
  (lT1Can1[i] <= tk_size - 1)
);

% Type 2 Cancellation (Auxiliary para "a")
constraint forall(r in  rEb..(rEb+rDis+rEf-1), c in 0..3)(
  lT2a0[r,c] = sum(i in 0..3)(if (lZ0[r,4*c+i] >= 1) then 1 else 0 endif) /\
  lT2a1[r,c] = sum(i in 0..3)(if (lZ1[r,4*c+i] >= 1) then 1 else 0 endif)
);

% Type 2 Cancellation (Auxiliary para "b")
constraint forall(r in  rEb..(rEb+rDis+rEf-1), c in 0..3)(
  lT2b0[r,c] = if (lAC0[r,c] >= 1) then (4 - sum(i in 0..3)(if (lW0[r,4*c+i] >= 1) then 1 else 0 endif)) else 0 endif /\
  lT2b1[r,c] = if (lAC1[r,c] >= 1) then (4 - sum(i in 0..3)(if (lW1[r,4*c+i] >= 1) then 1 else 0 endif)) else 0 endif
);

% Type 2 Cancellation (Auxiliary para "c")
constraint forall(r in  rEb..(rEb+rDis+rEf-1), c in 0..3)(
  lT2c0[r,c] = sum(i in 0..3)(lcan0[r+1,4*c+i]) /\
  lT2c1[r,c] = sum(i in 0..3)(lcan1[r+1,4*c+i])
);

% Type 2 Cancellation (main)
constraint forall(r in rEb..(rEb+rDis+rEf-1), c in 0..3)(
  lT2Can0[r,c] = max(0, (lT2b0[r,c] + lT2c0[r,c] - lT2a0[r,c])) /\
  lT2Can1[r,c] = max(0, (lT2b1[r,c] + lT2c1[r,c] - lT2a1[r,c]))
);

% Counting of all Cancellations
constraint 
  (tk_size * sum(i in 0..15)(lLANE0[i]) - sum(i in 0..15)(lT1Can0[i]) - sum(r in rEb..(rEb+rDis+rEf-1), c in 0..3)(lT2Can0[r,c])) >= 1 /\
  (tk_size * sum(i in 0..15)(lLANE1[i]) - sum(i in 0..15)(lT1Can1[i]) - sum(r in rEb..(rEb+rDis+rEf-1), c in 0..3)(lT2Can1[r,c])) >= 1
;


% --=========================================================--
% ||============== Lower Guess-and-Determine ================||
% --=========================================================--
array[(rEb+rDis+1)..(rEb+rDis+rEf),   0..15] of var 0 .. 1: lGeqk0;
array[(rEb+rDis)..(rEb+rDis+rEf-2),   0..15] of var 0 .. 1: lDeqW0;
array[(rEb+rDis)..(rEb+rDis+rEf-1),   0..15] of var 0 .. 1: lDeqX0;
array[(rEb+rDis)..(rEb+rDis+rEf-1),   0..15] of var 0 .. 1: lDeqY0;
array[(rEb+rDis)..(rEb+rDis+rEf-1),   0..15] of var 0 .. 1: lDeqZ0;
array[(rEb+rDis)..(rEb+rDis+rEf-1),   0..15] of var 0 .. 1: lFrSB0;
array[(rEb+rDis)..(rEb+rDis+rEf-2),   0..15] of var 0 .. 1: lFrMC0;
constraint forall(i in 0..15)(lDeqZ0[rEb+rDis+rEf-1,i] = 1);

array[(rEb+rDis+1)..(rEb+rDis+rEf),   0..15] of var 0 .. 1: lGeqk1;
array[(rEb+rDis)..(rEb+rDis+rEf-2),   0..15] of var 0 .. 1: lDeqW1;
array[(rEb+rDis)..(rEb+rDis+rEf-1),   0..15] of var 0 .. 1: lDeqX1;
array[(rEb+rDis)..(rEb+rDis+rEf-1),   0..15] of var 0 .. 1: lDeqY1;
array[(rEb+rDis)..(rEb+rDis+rEf-1),   0..15] of var 0 .. 1: lDeqZ1;
array[(rEb+rDis)..(rEb+rDis+rEf-1),   0..15] of var 0 .. 1: lFrSB1;
array[(rEb+rDis)..(rEb+rDis+rEf-2),   0..15] of var 0 .. 1: lFrMC1;
constraint forall(i in 0..15)(lDeqZ1[rEb+rDis+rEf-1,i] = 1);

% Guess and determine
constraint forall(r in (rEb+rDis)..(rEb+rDis+rEf-1), i in 0..15)(
  if (lDeqZ0[r,i] == 1) /\ (lGeqk0[r+1,i] == 1) then (lDeqY0[r,i] = 1)
  else (lDeqY0[r,i] = 0)
  endif
  /\
  if (lDeqZ1[r,i] == 1) /\ (lGeqk1[r+1,i] == 1) then (lDeqY1[r,i] = 1)
  else (lDeqY1[r,i] = 0)
  endif
);

% "Det" via SB
constraint forall(r in (rEb+rDis)..(rEb+rDis+rEf-1), i in 0..15)(
  lDeqX0[r,i] = lDeqY0[r,i] 
  /\
  lDeqX1[r,i] = lDeqY1[r,i]
);


% "Det" via MC
constraint forall(r in (rEb+rDis+1)..(rEb+rDis+rEf-1), c in 0..3)(
  if (sum(i in 0..3)(lDeqX0[r,4*c+i]) == 4) then (sum(i in 0..3)(lDeqW0[r-1,4*c+i]) == 4) 
  else (sum(i in 0..3)(lDeqW0[r-1,4*c+i]) == 0)
  endif
  /\
  if (sum(i in 0..3)(lDeqX1[r,4*c+i]) == 4) then (sum(i in 0..3)(lDeqW1[r-1,4*c+i]) == 4) 
  else (sum(i in 0..3)(lDeqW1[r-1,4*c+i]) == 0)
  endif
);

% "Det" via SR
constraint forall(r in (rEb+rDis)..(rEb+rDis+rEf-2), i in 0..15)(
  lDeqZ0[r,SRp[i]] = lDeqW0[r,i]
  /\
  lDeqZ1[r,SRp[i]] = lDeqW1[r,i]
);

% Filters from SB
constraint forall(r in (rEb+rDis)..(rEb+rDis+rEf-1), i in 0..15)(
  if (lDeqX0[r,i] == 1) /\ (leqX0[r,i] == 1) then (lFrSB0[r,i] = 1)
  else (lFrSB0[r,i] = 0)
  endif
  /\
  if (lDeqX1[r,i] == 1) /\ (leqX1[r,i] == 1) then (lFrSB1[r,i] = 1)
  else (lFrSB1[r,i] = 0)
  endif
);

% Filter from MC
constraint forall(r in (rEb+rDis)..(rEb+rDis+rEf-2), c in 0..3)(
  if exists(i in 0..3)(leqX0[r+1,4*c+i] >= 2)
  then forall(i in 0..3) (if (lDeqW0[r,4*c+i] == 1) /\ (leqW0[r,4*c+i] <= 1) then (lFrMC0[r,4*c+i] == 1) else (lFrMC0[r,4*c+i] == 0) endif)
  else sum(i in 0..3)(lFrMC0[r,4*c+i]) = 0
  endif
  /\
  if exists(i in 0..3)(leqX1[r+1,4*c+i] >= 2)
  then forall(i in 0..3) (if (lDeqW1[r,4*c+i] == 1) /\ (leqW1[r,4*c+i] <= 1) then (lFrMC1[r,4*c+i] == 1) else (lFrMC1[r,4*c+i] == 0) endif)
  else sum(i in 0..3)(lFrMC1[r,4*c+i]) = 0
  endif
);

% %%%%%%%%%%%%%%%%%%%%%%%%%% TEST Ef %%%%%%%%%%%%%%%%%%%%%%%%%%%
% constraint forall(i in [12,13])(lstk0[7,i] = 1 /\ lstk1[7,i] = 1) /\ sum(i in 0..15)(lstk0[8,i]) = 0 /\ sum(i in 0..15)(lstk1[8,i]) = 0;
% constraint forall(r in 7..8, i in 0..15)(lW0[7,i] = 0 /\ lW1[7,i] = 0 /\ lW0[8,i] = 0 /\ lW1[8,i] = 0);
% constraint sum(i in 0..15)(lX0[9,i]) = 2 /\ sum(i in 0..15)(lX1[9,i]) = 2;
% constraint forall(i in [12,14,15])(lZ0[6,i] = 1 /\ lZ1[6,i] = 1) /\ sum(i in 0..15)(lZ0[6,i]) = 3 /\ sum(i in 0..15)(lZ1[6,i]) = 3;
% constraint forall(i in [8,9,11,12,13,15])(leqZ0[10,i] = 0 /\ leqZ1[10,i] = 0); 
% constraint forall(i in [0,1,2,3,4,5,6,7,10])(lGeqk0[11,i] = 1 /\ lGeqk1[11,i] = 1) /\ sum(i in 0..15)(lGeqk0[11,i]) = 9 /\ sum(i in 0..15)(lGeqk1[11,i]) = 9;
% constraint forall(i in [3,10])(lGeqk0[10,i] = 1 /\ lGeqk1[10,i] = 1) /\ sum(i in 0..15)(lGeqk0[10,i]) = 2 /\ sum(i in 0..15)(lGeqk1[10,i]) = 2;
% %%%%%%%%%%%%%%%%%%%%%%%%%% TEST Ef %%%%%%%%%%%%%%%%%%%%%%%%%%%



% +++++++++++++++--------------++++++++++++++++
% |            Contradictions                 |
% +++++++++++++++--------------++++++++++++++++
array [0..22, 0..3] of var 0..2: Diff_table = array2d(0..22, 0..3, [
  0,0,0,1,  0,0,1,0,  0,1,0,0,  1,0,0,0,
  0,0,1,1,  0,1,0,1,  1,0,0,1,  0,1,1,0,  1,0,1,0,  1,1,0,0,
  0,1,1,1,  1,0,1,1,  1,1,0,1,  1,1,1,0,
  1,1,1,1,
  2,0,0,0,  0,2,0,0,  0,0,2,0,  0,0,0,2,
  2,1,1,0,  2,0,1,1,  2,1,0,1,  2,1,1,1
]);

array[rEb..(rEb+rDis-1), 0..15] of var 0..1: flagX;
array[rEb..(rEb+rDis-1), 0..15] of var 0..1: flagZ;
array[rEb..(rEb+rDis-1), 0..15] of var 0..1: flagW;
array[rEb..(rEb+rDis-1), 0..15] of var 0..1: flagBCT;
array[rEb..(rEb+rDis-1), 0..15] of var 0..1: flagGBCT;

constraint forall(r in rEb..(rEb+rDis-1), i in 0..15)(
  if exists(row in 0..18)([uX0[r,i], uX1[r,i], lX0[r,i], lX1[r,i]] = [Diff_table[row,j] | j in 0..3]) then flagX[r,i] = 1 else flagX[r,i] = 0 endif
  /\
  if exists(row in 0..18)([uZ0[r,i], uZ1[r,i], lZ0[r,i], lZ1[r,i]] = [Diff_table[row,j] | j in 0..3]) then flagZ[r,i] = 1 else flagZ[r,i] = 0 endif
  /\
  if exists(row in 0..18)([uW0[r,i], uW1[r,i], lW0[r,i], lW1[r,i]] = [Diff_table[row,j] | j in 0..3]) then flagW[r,i] = 1 else flagW[r,i] = 0 endif
  /\
  if [uX0[r,i], uX1[r,i], lY0[r,i], lY1[r,i]] = [1,1,1,1] then flagBCT[r,i] = 1 else flagBCT[r,i] = 0 endif
  /\
  if exists(row in 0..15)([uX0[r,i], uX1[r,i], lY0[r,i], lY1[r,i]] = [Diff_table[row,j] | j in 0..3]) then flagGBCT[r,i] = 1 else flagGBCT[r,i] = 0 endif
);

array[rEb..(rEb+rDis-1), 0..15] of var 0..3: uZs;
array[rEb..(rEb+rDis-1), 0..15] of var 0..3: lWs;
array[rEb..(rEb+rDis-1), 0..3]  of var 0..3: uACs;
array[rEb..(rEb+rDis-1), 0..3]  of var 0..3: lACs;
array[rEb..(rEb+rDis-1), 0..3]  of var 0..1: flagMC;

constraint forall(r in rEb..(rEb+rDis-1), i in 0..15)(XOR(uZ0[r,i], uZ1[r,i], uZs[r,i]));
constraint forall(r in rEb..(rEb+rDis-1), i in 0..15)(XOR(lW0[r,i], lW1[r,i], lWs[r,i]));
constraint forall(r in rEb..(rEb+rDis-1), c in 0..3)(uACs[r,c] = max(max(uZs[r,4*c],uZs[r,4*c+1]),max(uZs[r,4*c+2],uZs[r,4*c+3])));
constraint forall(r in rEb..(rEb+rDis-1), c in 0..3)(lACs[r,c] = max(max(lWs[r,4*c],lWs[r,4*c+1]),max(lWs[r,4*c+2],lWs[r,4*c+3])));

constraint forall(r in rEb..(rEb+rDis-1), c in 0..3)(
  if (uACs[r,c] <= 2) /\ (lACs[r,c] <= 2) /\ 
     ((sum(i in 0..3)(if uZs[r,4*c+i] >= 1 then 1 else 0 endif)) + (sum(i in 0..3)(if lWs[r,4*c+i] >= 1 then 1 else 0 endif)) < 5)
  then (flagMC[r,c] = 1)
  else (flagMC[r,c] = 0)
  endif
  )
;

constraint sum(r in rEb..(rEb+rDis-1), i in 0..15)(flagX[r,i] + flagZ[r,i] + flagW[r,i] + flagBCT[r,i] + flagGBCT[r,i]) >= 1;

% ((((((((((((((((((((((((((((((((((((((((((((((((
% (                 Complexity                   )
% ))))))))))))))))))))))))))))))))))))))))))))))))
int: setX;
var float: tz;
constraint tz = log2(0.7 * setX);

% Para. counts truncated state and subkey in extension.
var int: rb0;  var int: rb1;  var int: rb;   
var int: cb0;  var int: cb1;  var int: cb;   
var int: mb0;  var int: mb1;  var int: mb;   
var int: mbp;
var int: cbp0; var int: cbp1;

var int: rf0;  var int: rf1;  var int: rf;
var int: cf0;  var int: cf1;  var int: cf;
var int: mf0;  var int: mf1;  var int: mf;
var int: mfp;
var int: cfp0; var int: cfp1;

% Truncated in plaintext
constraint rb0 = sum(i in 0..15)(if uW0[-1,i] >= 2 then 1 else 0 endif);
constraint rb1 = sum(i in 0..15)(if uW1[-1,i] >= 2 then 1 else 0 endif);
constraint rb  = sum(i in 0..15)(if (uW0[-1,i] >= 2) \/ (uW1[-1,i] >= 2) then 1 else 0 endif);

% Truncated in ciphertext
constraint rf0 = sum(i in 0..15)(if leqZ0[rEb+rDis+rEf-1,i] >= 2 then 1 else 0 endif);
constraint rf1 = sum(i in 0..15)(if leqZ1[rEb+rDis+rEf-1,i] >= 2 then 1 else 0 endif);
constraint rf  = sum(i in 0..15)(if (leqZ0[rEb+rDis+rEf-1,i] >= 2) \/ (leqZ1[rEb+rDis+rEf-1,i] >= 2) then 1 else 0 endif);

% Truncated in head of Dis.
constraint cb0 = rb0 - sum(i in 0..15)(if uX0[rEb,i] >= 2 then 1 else 0 endif);
constraint cb1 = rb1 - sum(i in 0..15)(if uX1[rEb,i] >= 2 then 1 else 0 endif);
constraint cb  = rb  - sum(i in 0..15)(if (uX0[rEb,i] >= 2) \/ (uX1[rEb,i] >= 2) then 1 else 0 endif);

% Truncated in tail of Dis.
constraint cf0 = rf0 - sum(i in 0..15)(if lW0[rEb+rDis-1,i] >= 2 then 1 else 0 endif);
constraint cf1 = rf1 - sum(i in 0..15)(if lW1[rEb+rDis-1,i] >= 2 then 1 else 0 endif);
constraint cf  = rf  - sum(i in 0..15)(if (lW0[rEb+rDis-1,i] >= 2) \/ (lW1[rEb+rDis-1,i] >= 2) then 1 else 0 endif);

% Involved key in Eb
constraint mb0 = sum(r in 0..(rEb-1), i in 0..15)(if uX0[r,i] >= 2 then 1 else 0 endif);
constraint mb1 = sum(r in 0..(rEb-1), i in 0..15)(if uX1[r,i] >= 2 then 1 else 0 endif);
constraint mb  = sum(r in 0..(rEb-1), i in 0..15)(if (uX0[r,i] >= 2) \/ (uX1[r,i] >= 2) then 1 else 0 endif);

% Pre-guessed Key in Eb
constraint mbp = sum(r in 0..(rEb-1), i in 0..15)(if (uGstk0[r,i] + uGstk1[r,i] == 2) then 1 else (uGstk0[r,i] + uGstk1[r,i]) endif);

% Involved key in Ef
constraint mf0 = sum(r in (rEb+rDis)..(rEb+rDis+rEf-1), i in 0..15)(if leqZ0[r,i] >= 2 then 1 else 0 endif);
constraint mf1 = sum(r in (rEb+rDis)..(rEb+rDis+rEf-1), i in 0..15)(if leqZ1[r,i] >= 2 then 1 else 0 endif);
constraint mf  = sum(r in (rEb+rDis)..(rEb+rDis+rEf-1), i in 0..15)(if (leqZ0[r,i] >= 2) \/ (leqZ1[r,i] >= 2) then 1 else 0 endif);

% Pre-guessed Key in Ef
constraint mfp = sum(r in (rEb+rDis+1)..(rEb+rDis+rEf), i in 0..15)(if (lGeqk0[r,i] + lGeqk1[r,i] == 2) then 1 else (lGeqk0[r,i] + lGeqk1[r,i]) endif);

% Filter obtained in Eb
constraint cbp0 = sum(r in 0..(rEb-1), i in 0..15)(uFrSB0[r,i] + uFrMC0[r,i]);
constraint cbp1 = sum(r in 0..(rEb-1), i in 0..15)(uFrSB1[r,i] + uFrMC0[r,i]);

% Filter obtained in Ef
constraint cfp0 = sum(r in (rEb+rDis)..(rEb+rDis+rEf-1), i in 0..15)(lFrSB0[r,i]) + sum(r in (rEb+rDis)..(rEb+rDis+rEf-2), i in 0..15)(lFrMC0[r,i]);
constraint cfp1 = sum(r in (rEb+rDis)..(rEb+rDis+rEf-1), i in 0..15)(lFrSB1[r,i]) + sum(r in (rEb+rDis)..(rEb+rDis+rEf-2), i in 0..15)(lFrMC0[r,i]);


% Para. for T,D,M
var float: D;  var float: Dc; var float: Qc;
var float: Mc;
var float: T0;  
var float: T1;  
var float: T2p; var float: T2c;  var float: T2;  
var float: T31; var float: T32;  
var float: T4; 
var float: Tc;

constraint D = b_size + tz;
constraint Dc = D + 2;
constraint Qc = c_size*2*(cb+cf) + tz;

constraint T0 = Dc;
constraint T1 = c_size*(mbp + mfp) + Dc + log2(max(rEb,rEf) / (rEb+rDis+rEf));

constraint 
  T2p = max(c_size*(mbp + mfp) + D + c_size*(rb0-cbp0), c_size*(mbp + mfp) + Dc + c_size*(rb1-cbp1)) /\
  T2c = max(c_size*(mbp + mfp) + D + (D - tz + c_size*(rf0 - cfp0) - b_size), c_size*(mbp + mfp) + D + (D - tz + c_size*(rf1 - cfp1) - b_size)) /\
  T2  = min(T2p, T2c)
;

constraint T31 = c_size*(mbp + mfp - cbp0 - cbp1 - cfp0 - cfp1) + Qc;
constraint T32 = c_size*(mb + mf) + tz; % without consideration of epsilon
constraint T4  = b_size*tk_size - setX;

constraint (Tc >= T0) /\ (Tc >= T1) /\ (Tc >= T2) /\ (Tc >= T31) /\ (Tc >= T32) /\ (Tc >= T4);
constraint Tc < tk_size * b_size;

constraint Mc = max(max(Dc, min(T2p, T2c) - c_size*(mbp + mfp)), T31 - c_size*(mbp + mfp));
% OBJ
solve minimize Tc;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PREDICATES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
predicate XOR(var 0..3: Di0, var 0..3: Di1, var 0..3: Do) = 
  if (Di0 == 0) then (Do = Di1)
  elseif (Di1 == 0) then (Do = Di0)
  elseif (Di0 == 1) /\ (Di1 == 1) then (Do = 1) \/ (Do = 0)
  else (Do = 3)
  endif
;
  
predicate fATK(var 0..3: dw, var 0..1: stk, var 0..3: dx, var 0..1: can) = 
  if     (dw == 0)               then (dx = stk) /\ (can = 0)
  elseif (dw == 1) /\ (stk == 0) then (dx = 1) /\ (can = 0)
  elseif (dw == 1) /\ (stk == 1) then (dx + can = 1)
  elseif (dw == 2) /\ (stk == 0) then (dx = 2) /\ (can = 0)
  else                                (dx = 3) /\ (can = 0)
  endif
;

predicate bATK(var 0..3: dx, var 0..1: stk, var 0..3: dw, var 0..1: can) = 
  if     (dx == 0) /\ (stk == 0) then (dw = 0) /\ (can = 0)
  elseif (dx == 0) /\ (stk == 1) then (dw = 1) /\ (can = 1)
  elseif (dx == 2) /\ (stk == 0) then (dw = 2) /\ (can = 0)
  else                                (dw = 3) /\ (can = 0)
  endif
;

predicate SBox(var 0..3: din, var 0..3: dout) = 
  if     (din == 0) then (dout = 0)
  elseif (din <= 2) then (dout = 2)
  else                   (dout = 3)
  endif
;

predicate MC(var 0..3: din0, var 0..3: din1, var 0..3: din2, var 0..3: din3,
             var 0..3: tAC,
             var 0..3: dout0, var 0..3: dout1, var 0..3: dout2, var 0..3: dout3) = 
  tAC = max(max(din0,din1),max(din2,din3)) /\
  if (tAC == 0) then (
    (dout0 = 0) /\ (dout1 = 0) /\ (dout2 = 0) /\ (dout3 = 0)
  )
  elseif (tAC == 1) then (
    (dout0 <= 1) /\ (dout1 <= 1) /\ (dout2 <= 1) /\ (dout3 <= 1) /\
    (din0 + din1 + din2 + din3 + dout0 + dout1 + dout2 + dout3 >= 5) /\
    (din0 + din1 + din2 + din3 + dout0 + dout1 + dout2 + dout3 <= 8)
  )
  elseif (tAC == 2) /\ (din0 + din1 + din2 + din3 == 2) then (
    (dout0 = 2) /\ (dout1 = 2) /\ (dout2 = 2) /\ (dout3 = 2)
  )
  else (
    (dout0 = 3) /\ (dout1 = 3) /\ (dout2 = 3) /\ (dout3 = 3)
  )
  endif
;