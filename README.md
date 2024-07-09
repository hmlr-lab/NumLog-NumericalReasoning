# NumLog : A Numerical Reasoning using ILP
The learning methodology for numerical anaylsis <br>
Author:  Daniel Cyrus

<p>Numlog is a learning system to explain and generate a sort of numerical values. Earlier versions of numlog (e.g. numlog.pl, numlog_v2.pl and numlog_v3.pl) are based on logic approach </p>
<p>To use the lates version and for ILP approaches, please use numlogILP.pl</p>
<p>NumLog aims to generate rules that are both highly accurate and easily understandable, ensuring low complexity for better explainability.</p>

#### Usage
1. Copy and paset NumLogILP to your current project directory.<br>
2.import NumLog as a module in you Prolog code:<br>

``` Prolog
:-use_module('NumLogILP').
```

3.If you want to learn only from a single numerical values easily use learn(Pos,Neg,Method):

``` Prolog
:-
   Pos=[1,2,3,5,6,8,9],
   Neg = [4,7],
  learn(Pos,neg,threshold).
```

To learn from multi class values create background and example file and run below code:

``` Prolog
:-  
    time(learnFromFile('examples/bk.pl','examples.pl')).
```

#### Sample generated rule for a single numerical range:

``` Prolog
leq(A,3.8).
geq(A,7.2).
inRange(A,5.5±1).
```

#### Sampel generated rule for multi class values:

``` Prolog
cancer(A):-
          worst_radius(A,B),leq(B, 14.482),
          worst_fractal_dimension(A,C),inRange(C, 0.095730624±0.040690624).
```

#### Citation and contact
This work is based on the <a href="">paper</a>, for any issues or question please contact author (d.cysur@surrey.ac.uk).
If you using this work please cite:
<Citation and paper link will be provided soon>
