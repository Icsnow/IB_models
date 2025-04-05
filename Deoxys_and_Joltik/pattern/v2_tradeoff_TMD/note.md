



Using "_IB_JandD v2.py", we can simultaneously optimize the time, data and memory complexity with the weight of objective function:
```python
self.model.ModelSense = GRB.MINIMIZE
self.model.setObjectiveN(Tc, index=0, priority=4, name='Min_Tc')
self.model.setObjectiveN(T32, index=1, priority=3, name='Min_T32')
self.model.setObjectiveN(Dc, index=2, priority=2, name='Min_Dc')
self.model.setObjectiveN(Mc, index=3, priority=1, name='Min_Mc')
```


The pattern returned is the same as the result of "_IB_JandD.py." (Different patterns with the same complexities)