# The MILP/CP models for the paper "A Holistic Framework for Impossible Boomerang Attacks"

## How to use this tool

### v_1,2,3 (MILP)
Environment needed: [Gurobi](https://www.gurobi.com/) + Python

"\_\_main__" can be modified (example):
```
''' 10-r Deoxys-BC-256 (v1,2) '''
DeoxysBC256_10r = IB_DandJ(key_size, rEb, rEu, rEm, rEl, rEf, SetX, Construct_pairs_in_Eb)
DeoxysBC256_10r.ib_model()
```
```
''' 10-r Deoxys-BC-256 (v3)'''
DeoxysBC256_10r = IB_DandJ(key_size, rEb, rEd, rEf, SetX, Construct_pairs_in_Eb)
DeoxysBC256_10r.ib_model()
```

### v_4 (CP)
Environment needed: [Minizinc](https://www.minizinc.org/) + [Gurobi](https://www.gurobi.com/)

> We only tried 11-round IB on Deoxys-BC-256. (detail refer [v4_4trails](https://github.com/Icsnow/IB_models/tree/main/v4_4trails))
> 
Modify "bashSolve.bat" like:

`minizinc --solver THESOLVER -p #THREADS MODEL -D "rEb; rDis; rEf; setX" > OUTPUTFILE` 

Example:

`minizinc --solver Gurobi -p 4 .\IB_cp_Deoxys-Joltik.mzn -D "rEb=2; rDis=7; rEf=2; setX=24" > .\res_IBDJ_2_7_2_24.txt`

## Discussion about multi-objective optimization
Thanks to the comments of the reviewer, we discuss the comparison between "Linear weighted multi-objectives" and "Stage multi-objectives" in the following context.

* **Linear weighted multi-objectives** are an approximation of a Pareto frontier; that is, the optimization of one target will affect the optimization of all other targets, which may cause other goals to be non-optimized (although their weight is greater;

$$
\alpha * Tc + \beta * Dc + \gamma * Mc, \text{ where } \alpha + \beta + \gamma = 1
$$

* **Stage multi-objective optimization** in MILP (which we implement using the Python interface) will optimize other parts based on the already optimized maximum weighted target. The optimal value of the target with the maximum weight is fixed when optimizing other targets.

```
ModelSense = GRB.MINIMIZE
setObjectiveN(Tc, index = 0, priority = 3)
setObjectiveN(Dc, index = 1, priority = 2)
setObjectiveN(Mc, index = 2, priority = 1)
```

In our model, time complexity is related to the maximum weight, and we want to maintain the optimal time complexity while optimizing data and memory complexity.
