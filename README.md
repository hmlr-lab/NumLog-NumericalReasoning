# NumLog: Numerical-Symbolic Learning and Reasoning 
The learning methodology for numerical anaylsis <br>
This repository is based on the paper in <a href="https://hmlr-lab.github.io/pdfs/NumLog_ILP2024.pdf">Link</a><br>
Main author:  Daniel Cyrus

<p>Numlog is a learning system to explain and generate a sort of numerical values. Earlier versions of numlog (e.g. numlog.pl, numlog_v2.pl and numlog_v3.pl) are based on logic approach </p>
<p>To use the latest version and for ILP approaches, please use numlog.pl</p>
<p>NumLog aims to generate rules that are both highly accurate and easily understandable, ensuring low complexity for better explainability.</p>

#### Usage
1. Copy and paset NumLogNew to your current project directory.<br>
2. Import NumLog as a module in your Prolog code:<br>

``` Prolog
:-style_check(-discontiguous).
:-use_module(numlog).

numlog:combination(0.15).
numlog:maxClause(4).
numlog:minClause(1).
```
3. To learn only from a single numerical values easily use learn(Pos,Neg,Method):

``` Prolog
:-
   Pos=[1,2,3,5,6,8,9],
   Neg = [4,7],
  learn(Pos,Neg,_,_,print).
```
or
``` Prolog
:-
   Pos=[1,2,3,5,6,8,9],
   Neg = [4,7],
   learn(Pos,Neg,R,Acc,noprint),
   write(R),nl,
   write(Acc),nl.
```
To learn from multi class values create a background and an example file and run below code:

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
For any issues or question please contact author (d.cyrus@surrey.ac.uk).<br>
If you are suing NumLog in your work, please cite:<br>
D. Cyrus, D. Varghese, and A. Tamaddoni-Nezhad, <b>An Inductive Logic Programming approach for feature-range discovery</b>,
In Proc. of the 33rd Int. Conf. on ILP , Springer, 2024 
