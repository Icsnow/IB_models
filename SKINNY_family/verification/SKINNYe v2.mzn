% TEST case: AC22-Mind the TWEAKEY Schedule: Cryptanalysis on SKINNYe-64-256

array[1..15, 1..15] of var -1..15: iBCT = array2d(1..15, 1..15,
[% ============1=============
1,-1,3,4,5,6,7,-1,-1,-1,-1,12,13,14,15,
% ============2=============
-1,2,-1,-1,-1,-1,7,8,9,10,11,12,13,14,15,
% ============3=============
1,2,3,4,5,6,7,-1,-1,-1,-1,-1,-1,-1,-1,
% ============4=============
1,-1,3,4,5,-1,-1,-1,-1,-1,-1,-1,-1,14,15,
% ============5=============
1,-1,3,4,5,-1,-1,-1,-1,-1,-1,-1,-1,14,15,
% ============6=============
-1,2,-1,-1,5,6,-1,-1,9,-1,11,12,-1,-1,15,
% ============7=============
-1,2,-1,-1,5,6,-1,8,-1,10,-1,-1,13,14,-1,
% ============8=============
-1,2,-1,-1,-1,-1,7,8,9,10,11,-1,-1,-1,-1,
% ============9=============
-1,2,-1,-1,-1,-1,7,8,9,10,11,-1,-1,-1,-1,
% ============10=============
-1,2,-1,-1,-1,-1,7,-1,-1,-1,-1,12,13,14,15,
% ============11=============
-1,2,-1,-1,-1,-1,7,8,9,10,11,-1,-1,-1,-1,
% ============12=============
1,-1,3,4,5,-1,-1,-1,-1,-1,-1,12,13,-1,-1,
% ============13=============
1,-1,3,4,5,-1,-1,-1,-1,-1,-1,12,13,-1,-1,
% ============14=============
-1,2,-1,-1,5,6,-1,8,-1,10,-1,12,-1,-1,15,
% ============15=============
-1,2,-1,-1,5,6,-1,-1,9,-1,11,-1,13,14,-1
]);

array[0..15, 0..15] of var 0..15: XOR = array2d(0..15, 0..15,
[
%==========0=========
0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,
%==========1=========
1,0,3,2,5,4,7,6,9,8,11,10,13,12,15,14,
%==========2=========
2,3,0,1,6,7,4,5,10,11,8,9,14,15,12,13,
%==========3=========
3,2,1,0,7,6,5,4,11,10,9,8,15,14,13,12,
%==========4=========
4,5,6,7,0,1,2,3,12,13,14,15,8,9,10,11,
%==========5=========
5,4,7,6,1,0,3,2,13,12,15,14,9,8,11,10,
%==========6=========
6,7,4,5,2,3,0,1,14,15,12,13,10,11,8,9,
%==========7=========
7,6,5,4,3,2,1,0,15,14,13,12,11,10,9,8,
%==========8=========
8,9,10,11,12,13,14,15,0,1,2,3,4,5,6,7,
%==========9=========
9,8,11,10,13,12,15,14,1,0,3,2,5,4,7,6,
%==========10=========
10,11,8,9,14,15,12,13,2,3,0,1,6,7,4,5,
%==========11=========
11,10,9,8,15,14,13,12,3,2,1,0,7,6,5,4,
%==========12=========
12,13,14,15,8,9,10,11,4,5,6,7,0,1,2,3,
%==========13=========
13,12,15,14,9,8,11,10,5,4,7,6,1,0,3,2,
%==========14=========
14,15,12,13,10,11,8,9,6,7,4,5,2,3,0,1,
%==========15=========
15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,0
]);

% TK2
array[0..15] of var 0..15: LFSR2 = array1d(0..15, [0,2,4,6,9,11,13,15,1,3,5,7,8,10,12,14]);
% TK3
array[0..15] of var 0..15: LFSR3 = array1d(0..15, [0,8,1,9,2,10,3,11,12,4,13,5,14,6,15,7]);
% TK4
array[0..15] of var 0..15: LFSR4 = array1d(0..15, [0,4,9,13,3,7,10,14,2,6,11,15,1,5,8,12]);



array[0..15] of var 0..15: hperm = array1d(0..15, [8,9,10,11,12,13,14,15,2,0,4,7,6,3,5,1]);

int: tR = 16; % attack round//2
int: totR = 35;

array[0..totR, 0..15] of var 0..15: utk1;

array[0..totR, 0..15] of var 0..15: utk2;
array[0..totR, 0..15] of var 0..15: up2;

array[0..totR, 0..15] of var 0..15: utk3;
array[0..totR, 0..15] of var 0..15: up3;

array[0..totR, 0..15] of var 0..15: utk4;
array[0..totR, 0..15] of var 0..15: up4;

array[0..totR, 0..7] of var 0..15: uSTK;

% ---------------------------------------
array[0..totR, 0..15] of var 0..15: ltk1;

array[0..totR, 0..15] of var 0..15: ltk2;
array[0..totR, 0..15] of var 0..15: lp2;

array[0..totR, 0..15] of var 0..15: ltk3;
array[0..totR, 0..15] of var 0..15: lp3;

array[0..totR, 0..15] of var 0..15: ltk4;
array[0..totR, 0..15] of var 0..15: lp4;

array[0..totR, 0..7] of var 0..15: lSTK;

% ============================================== KEY (begin) ===============================================
constraint forall(r in 0..totR-1, c in 0..15)(utk1[r+1,hperm[c]] = utk1[r,c]);

constraint forall(r in 1..totR, c in 0..15)(up2[r,hperm[c]] = utk2[r-1,c]);
constraint forall(r in 1..totR, c in 0..7)(utk2[r,c] = LFSR2[up2[r,c]]);
constraint forall(r in 1..totR, c in 8..15)(utk2[r,c] = up2[r,c]);

constraint forall(r in 1..totR, c in 0..15)(up3[r,hperm[c]] = utk3[r-1,c]);
constraint forall(r in 1..totR, c in 0..7)(utk3[r,c] = LFSR3[up3[r,c]]);
constraint forall(r in 1..totR, c in 8..15)(utk3[r,c] = up3[r,c]);

constraint forall(r in 1..totR, c in 0..15)(up4[r,hperm[c]] = utk4[r-1,c]);
constraint forall(r in 1..totR, c in 0..7)(utk4[r,c] = LFSR4[up4[r,c]]);
constraint forall(r in 1..totR, c in 8..15)(utk4[r,c] = up4[r,c]);

% STK
constraint forall(r in 0..totR-1, c in 0..7)(uSTK[r,c] = XOR[XOR[utk1[r,c],utk2[r,c]],XOR[utk3[r,c],utk4[r,c]]]);


% U TEST
% constraint forall(c in [0,2,4,5,6,7,8,9,10,12,13,14,15])(utk1[4,c] = 0 /\ utk2[4,c] = 0 /\ utk3[4,c] = 0 /\ utk4[4,c] = 0);
% constraint utk1[4,1] = 6 /\ utk1[4,3] = 2 /\ utk1[4,11] = 13;
% constraint utk2[4,1] = 9 /\ utk2[4,3] = 8 /\ utk2[4,11] = 3;
% constraint utk3[4,1] = 12 /\ utk3[4,3] = 11 /\ utk3[4,11] = 8;
% constraint utk4[4,1] = 10 /\ utk4[4,3] = 9 /\ utk4[4,11] = 5;
% -----------------------------------------

constraint forall(r in 0..totR-1, c in 0..15)(ltk1[r+1,hperm[c]] = ltk1[r,c]);

constraint forall(r in 1..totR, c in 0..15)(lp2[r,hperm[c]] = ltk2[r-1,c]);
constraint forall(r in 1..totR, c in 0..7)(ltk2[r,c] = LFSR2[lp2[r,c]]);
constraint forall(r in 1..totR, c in 8..15)(ltk2[r,c] = lp2[r,c]);

constraint forall(r in 1..totR, c in 0..15)(lp3[r,hperm[c]] = ltk3[r-1,c]);
constraint forall(r in 1..totR, c in 0..7)(ltk3[r,c] = LFSR3[lp3[r,c]]);
constraint forall(r in 1..totR, c in 8..15)(ltk3[r,c] = lp3[r,c]);

constraint forall(r in 1..totR, c in 0..15)(lp4[r,hperm[c]] = ltk4[r-1,c]);
constraint forall(r in 1..totR, c in 0..7)(ltk4[r,c] = LFSR4[lp4[r,c]]);
constraint forall(r in 1..totR, c in 8..15)(ltk4[r,c] = lp4[r,c]);

% STK
constraint forall(r in 0..totR-1, c in 0..7)(lSTK[r,c] = XOR[XOR[ltk1[r,c],ltk2[r,c]],XOR[ltk3[r,c],ltk4[r,c]]]);

% % L TEST
% constraint forall(c in 0..14)(ltk1[4,c] = 0 /\ ltk2[4,c] = 0 /\ ltk3[4,c] = 0 /\ ltk4[4,c] = 0);
% constraint ltk1[4,15] = 2;
% constraint ltk2[4,15] = 6;
% constraint ltk3[4,15] = 9;
% constraint ltk4[4,15] = 1;
% ============================================== KEY (end) ===============================================


% =================== Contradiction =====================
% STK16[3] = W16[3] = X17[7]
% iBCT[STK16[3],STK17[7]] = 0
constraint lSTK[17,7] in [iBCT[uSTK[16,3], i] | i in 1..15];

% PATTERN
% ---------- upper
constraint forall(c in [0,1,3,4,5,6,7])(uSTK[1,c] = 0);
constraint uSTK[1,2] >= 1;

constraint forall(c in [0,1,2,3,4,5,6])(uSTK[2,c] = 0);
constraint uSTK[2,7] >= 1;

constraint forall(c in [0,1,2,3,5,6,7])(uSTK[3,c] = 0);
constraint uSTK[3,4] >= 1;

constraint forall(c in [0,2,3,4,5,6,7])(uSTK[4,c] = 0);
constraint uSTK[4,1] >= 1;

constraint forall(r in [7,8,9,10,11,12], c in 0..7)(uSTK[r,c] = 0);

% ----------- lower
constraint forall(c in [0,1,2,3,4,5,6])(lSTK[17,c] = 0);
constraint lSTK[17,7] >= 1;

constraint forall(c in [0,2,3,4,5,6,7])(lSTK[19,c] = 0);
constraint lSTK[19,1] >= 1;

constraint forall(c in [0,1,2,3,4,5,7])(lSTK[27,c] = 0);
constraint lSTK[27,6] >= 1;

constraint forall(r in [18,20,21,22,23,24,25,26,28], c in 0..7)(lSTK[r,c] = 0);


% TRY (Can be deleted)
constraint uSTK[16,3] = 10 /\ lSTK[17,7] = 2;