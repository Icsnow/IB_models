minizinc --solver Gurobi -p 4 .\IB_cp_Deoxys-Joltik.mzn -D "rEb=2; rDis=7; rEf=2; setX=24" > .\res_IBDJ_2_7_2_24.txt
minizinc --solver Gurobi -p 4 .\IB_cp_Deoxys-Joltik.mzn -D "rEb=1; rDis=7; rEf=3; setX=24" > .\res_IBDJ_1_7_3_24.txt
minizinc --solver Gurobi -p 4 .\IB_cp_Deoxys-Joltik.mzn -D "rEb=3; rDis=7; rEf=1; setX=24" > .\res_IBDJ_3_7_1_24.txt
minizinc --solver Gurobi -p 4 .\IB_cp_Deoxys-Joltik.mzn -D "rEb=2; rDis=8; rEf=1; setX=8" > .\res_IBDJ_2_8_1_8.txt
minizinc --solver Gurobi -p 4 .\IB_cp_Deoxys-Joltik.mzn -D "rEb=1; rDis=8; rEf=2; setX=8" > .\res_IBDJ_1_8_2_8.txt
minizinc --solver Gurobi -p 4 .\IB_cp_Deoxys-Joltik.mzn -D "rEb=2; rDis=8; rEf=1; setX=20" > .\res_IBDJ_2_8_1_20.txt
minizinc --solver Gurobi -p 4 .\IB_cp_Deoxys-Joltik.mzn -D "rEb=1; rDis=8; rEf=2; setX=20" > .\res_IBDJ_1_8_2_20.txt