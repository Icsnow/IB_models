# IB attack with more contradiction (CP)

> example: 11-round IB attack on Deoxys-BC-256

## Cells:

$$
\text{Upper:}\ \ \ 
\mathrm{X}\  \stackrel{\mathcal{SR\circ SC}}{\longrightarrow}\  \mathrm{Z}\  \stackrel{\mathcal{MC}}{\longrightarrow}\  \mathrm{W}\  \stackrel{\mathcal{AK}}{\longrightarrow}\  \mathrm{X} \\
\text{Lower:}\ \ \ 
\mathrm{X}\  \stackrel{\mathcal{SR\circ SC}}{\longleftarrow}\  \mathrm{Z}\  \stackrel{\mathcal{MC}}{\longleftarrow}\  \mathrm{W}\  \stackrel{\mathcal{AK}}{\longleftarrow}\  \mathrm{X} \\
$$

**There are four differential trails (two upper, two lower) are allowed not the same as each other.**

## 4 Different Representations for Each Cell:

$$
\begin{cases}
0 & \text{Zero diff.}\\
1 & \text{Fixed diff.}\\
2 & \text{Truncated non-zero diff.}\\
3 & \text{Arbitrary diff.}\\
\end{cases}
$$


## Differential Propagation:

### XOR

$$
\begin{align}
0+0 & =0 \\
0 (/1) + 1 (/0) & = 1\\
1+1 & =  0 (/1)\\
2 + 0 & =2\\
2 + (\ge 1) & =3
\end{align}
$$

### S-Box

$$
\begin{align}
0 \stackrel{SBox}{\Longrightarrow} & 0 \\
\{1,2\} \stackrel{SBox}{\Longrightarrow} & 2 \\
3 \stackrel{SBox}{\Longrightarrow} & 3
\end{align}
$$

### MixColumn (MDS)

Let $\Delta_ci$ and $\Delta_co$ be the difference of one column before and after MDS matrix.

Denote $AC$ is the MAX value of $\Delta_ci$ and $\Delta_co$.

---

if $\Delta_ci = [0,0,0,0]$ $\stackrel{MC}{\Longrightarrow}$ $\Delta_co = [0,0,0,0]$;

else if $\Delta_ci\in \{0,1\}$ $\stackrel{MC}{\Longrightarrow}$ $\Delta_co \in \{0,1\}$, follows the rule of MDS;

else if $\Delta_ci \in \{0,2\}$ and only one $2$ exists $\stackrel{MC}{\Longrightarrow}$ $\Delta_co=[2,2,2,2]$;

else $\Delta_co=[3,3,3,3]$.

---

## ==Contradiction==

Let $(a,b,c,d)$ be the four cells with the same position on 4 differential trails, where $a,b,c,d \in \{0,1,2,3\}$.

$3$ is isolated to the detection of contradiction.

This model contain three cases of contradiction:

### Case1: $a\oplus b\oplus c\oplus d\neq0$

The following contradictions are allowed to appear at any cell.

```
0,0,0,1,  0,0,1,0,  0,1,0,0,  1,0,0,0,
0,0,1,1,  0,1,0,1,  1,0,0,1,  0,1,1,0,  1,0,1,0,  1,1,0,0,
0,1,1,1,  1,0,1,1,  1,1,0,1,  1,1,1,0,
1,1,1,1,
2,0,0,0,  0,2,0,0,  0,0,2,0,  0,0,0,2,
2,1,1,0,  2,0,1,1,  2,1,0,1,  2,1,1,1
```

when $1$ appear more than one time, instantiation is needed

### Case2: Breaking MDS

Set $flagMC=1$ is the sign when "Breaking MDS" appearing.

Let $uZs\stackrel{MC}{\longrightarrow} lWs$ be the addition of two upper trails and the addition of two lower trails.

Let $uACs,lACs$ be the MAX value of $uZs, lWs$ in one column.

---

$\text{if } (uACs \le 2) \text{ and }(lACs \le 2) \text{ and } \sum_{column}{\text{ if }(uZs\ge1 \text{ then } 1) + \text{ if }(lWs\ge1 \text{ then } 1))<5} \Rightarrow flagMC=1$;

$\text{else } flagMC=0$

---

### Case3: Impossible BCT

Let $uX_0,uX_1,lY_0,lY_1$ be the difference before and after S-Box on four differential trails.

The following ccontradiction $(uX_0,uX_1,lY_0,lY_1)$  indicates Impossible BCT.

```
0,0,0,1, 0,0,1,0, 0,1,0,0, 1,0,0,0,
0,0,1,1, 0,1,0,1, 1,0,0,1, 0,1,1,0, 1,0,1,0, 1,1,0,0,
0,1,1,1, 1,0,1,1, 1,1,0,1, 1,1,1,0,
1,1,1,1,
```



## Solve

we have tried the following parameters, and marked each solution as its line number:

```
minizinc --solver Gurobi -p 4 .\IB_cp_Deoxys-Joltik.mzn -D "rEb=2; rDis=7; rEf=2; setX=24" > .\res_IBDJ_2_7_2_24.txt
minizinc --solver Gurobi -p 4 .\IB_cp_Deoxys-Joltik.mzn -D "rEb=1; rDis=7; rEf=3; setX=24" > .\res_IBDJ_1_7_3_24.txt
minizinc --solver Gurobi -p 4 .\IB_cp_Deoxys-Joltik.mzn -D "rEb=3; rDis=7; rEf=1; setX=24" > .\res_IBDJ_3_7_1_24.txt
minizinc --solver Gurobi -p 4 .\IB_cp_Deoxys-Joltik.mzn -D "rEb=2; rDis=8; rEf=1; setX=8" > .\res_IBDJ_2_8_1_8.txt
minizinc --solver Gurobi -p 4 .\IB_cp_Deoxys-Joltik.mzn -D "rEb=1; rDis=8; rEf=2; setX=8" > .\res_IBDJ_1_8_2_8.txt
minizinc --solver Gurobi -p 4 .\IB_cp_Deoxys-Joltik.mzn -D "rEb=2; rDis=8; rEf=1; setX=20" > .\res_IBDJ_2_8_1_20.txt
minizinc --solver Gurobi -p 4 .\IB_cp_Deoxys-Joltik.mzn -D "rEb=1; rDis=8; rEf=2; setX=20" > .\res_IBDJ_1_8_2_20.txt
```



## results

Consider the calculation of $\epsilon$, there is **no better attack** than previous versions.

The detail is list in file "res_IBDJ_1_7_3_24.txt"...
