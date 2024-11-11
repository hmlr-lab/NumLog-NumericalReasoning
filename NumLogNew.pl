%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Numlog: numerical learning and reasoning using ILP %%%%%%%%%%%%%%%%%%%%%%%%
%                                     Numlog Author: Daniel Cyrus                                       %
%                            This version is desined for first order logic                              % 
%                                            The HMLR lab                                               %
%                     for more detail see the paper http://                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
:- module(numlog, [learn/5,combination/1,maxClause/1,learnFromFile/2]).

:-op(500, xfy, user:(±)). 
:-op(500, xfy, =>).

:-dynamic epsilon/1.
:-dynamic combination/1.
:-dynamic maxClause/1.

:-dynamic numlogacc/1.
:-dynamic numlogrule/1.

:-dynamic(user:pos/1).
:-dynamic(user:neg/1).


numlogacc(0).
numlogrule([]).

pprint([],Acc):-writeln('--------------------------------'),format('Accuracy: ~2f~n',[Acc]).    
pprint([H|T],Acc):-
    H =.. [F,_,V],
    numbervars(A, 0, _),
    H1 =.. [F,A,V],
    format('~w.~n',[H1]),
    pprint(T,Acc).

%----------------------------Numlog Background knowledge-----------------------------------------------------
leq(List,Last):-
    last(List, Last).
geq([H|_],H).
inRange(List, Mid ± Range) :-
    [First|_] = List,
    last(List, Last),
    Mid is (Last + First) / 2,
    Range is Last - Mid.
   
%evaluation 

leq(A,B):-number(B),A=<B.
geq(A,B):-number(B),A>=B.

inRange(Val ,Mid ± Range):-
    number(Mid),
    number(Range),
    Val >= Mid - Range,
    Val =< Mid + Range.
%----------------------------Threshold analysis--------------------------------------------------------------
mean(A,Len,Mean):-
    sumlist(A, Sum),
    Mean is Sum / Len.

sigma([],_,Sigma,Sigma).
sigma([H|T],Mean,Temp,Sigma):-
    Temp1 is (H - Mean)**2 + Temp,
    sigma(T,Mean,Temp1,Sigma).

linespace(Min,Max,Number,Distribute):-
    R is Max - Min,
    Step is R / Number,
    distribute(Min,Step,Max,[Min],Distribute).

std(A,Mean,Len,Std):-
    length(A, Len),
    mean(A,Len,Mean),
    (Len = 1 ->
    (epsilon(Eps),
    [H|_] = A,
    Min is H - Eps,
    Max is H + Eps,
    linespace(Min,Max,9,Aall),
    length(Aall,Len2)
        );
    (
        
        Len >= 10 ->(
            Aall = A,Len2 = Len);
        (
        min_list(A, Min),
        max_list(A, Max),
        linespace(Min,Max,9,Aall),
        length(Aall,Len2)
        )
    )),
    sigma(Aall,Mean,0,Sigma),
    Len1 is Len2 - 1,
    Temp is Sigma / Len1,
    sqrt(Temp, Std).

recStd([],_, _, Exps,Exps).
recStd([H|T],Mean, Std, Temp,Exps):-
    Ans is -0.5 * ((H - Mean)/Std)**2,
    Ans1 is exp(Ans),
    S is 2 * pi,
    Val is Std * sqrt(S),
    Ans2 is 1 / (Val * Ans1),
    append(Temp, [Ans2], Temp1),
    recStd(T,Mean,Std,Temp1,Exps).

pdf(X, Mean, Std,PDF):-
    Ans is -0.5 * ((X - Mean)/Std)**2,
    Ans1 is exp(Ans),
    S is 2 * pi,
    Val is Std * sqrt(S),
    PDF is (1 / Val) * Ans1.

distribute(_,0,_,Arr,Arr).
distribute(Min,Step,Max,Arr,Arr1):-
    Min1 is Min + Step,
    (Min1 < Max ->
    (append(Arr,[Min1],Arr2),
    distribute(Min1,Step,Max,Arr2,Arr1));
    (append(Arr,[Max],Arr2),
    distribute(_,0,_,Arr2,Arr1))).

intersectDecrease(StartingPoint,Mean1,Mean2,Std1,Std2, Increment,Intersect):-
    pdf(StartingPoint,Mean1,Std1,PDF),
    pdf(StartingPoint,Mean2,Std2,PDF1),
    Ans is PDF - PDF1,
    (Ans < 0 ->
    (
    StartingPoint1 is StartingPoint - Increment,
    intersectDecrease(StartingPoint1,Mean1,Mean2,Std1,Std2, Increment,Intersect));
    Intersect = StartingPoint).

intersect(StartingPoint,Mean1,Mean2,Std1,Std2, Increment,Intersect):-
    pdf(StartingPoint,Mean1,Std1,PDF),
    pdf(StartingPoint,Mean2,Std2,PDF1),
    Ans is PDF - PDF1,
    (Ans > 0 ->
    (
    StartingPoint1 is StartingPoint + Increment,
    intersect(StartingPoint1,Mean1,Mean2,Std1,Std2, Increment,Intersect));
    %Intersect = StartingPoint
    intersectDecrease(StartingPoint,Mean1,Mean2,Std1,Std2, Increment,Intersect)).

%----------------------------Generate hypothesis space and create disjoint rules if possible-----------------
hypothesisSpace([],[]).
hypothesisSpace([_|T],Comb) :-
    hypothesisSpace(T,Comb).
hypothesisSpace([X|T],[X|Comb]):-
    hypothesisSpace(T,Comb).


flattenAll([],Res,Res).
flattenAll([H|T],Res,Res1):-
    flatten(H, Flat),
    append(Res,[Flat],Res2),
    flattenAll(T,Res2,Res1).

%-----------------------------------Set epsilon value-----------------------------------------------------------

calculateEpsilon(_,[],EPS,EPS).
calculateEpsilon(Prev,[H|T],EPS,EPS1):-
  Check is H - Prev,
  (Check < EPS ->
    EPS2 = Check;EPS2=EPS),
  calculateEpsilon(H,T,EPS2,EPS1).

setEpsilon(AllExamples):-
  retractall(epsilon(_)),
  [Prev|T] =AllExamples,
  once(calculateEpsilon(Prev,T,10,EPS)),
  EPS1 is (EPS / 10) + 0.000001,

  %round(5,EPS1,EPS2),
  assert(epsilon(EPS1)),
  format('Epsilon has set to ~8f \n',[EPS1]).
%-----------------------------------Process examples--------------------------------------------------------
groupExamples([],[],L,L1,G,G1):-
    (L \= [] -> append(G,[L],G2);G2 = G),
    (L1 \= [] -> append(G2,[L1],G3);G3 = G2),
    G1 = G3.
groupExamples(Pos,[],_,L,G,G1):-
    append(G,[L],G2),
    (Pos \= []->
    append(G2,[Pos],G1);G1=G2).
groupExamples([],Neg,L,_,G,G1):-
    append(G,[L],G2),
    (Neg \= []->
    append(G2,[Neg],G1);G1=G2).
groupExamples([H|T],[NH|NT],L,L1,G,G2):-
    (H \= NH ->
        (
        (L1 \= [] ->
            (append(G, [L1], G1),
            List1 = []);
            (List1 = [],
            G1 = G)),
        Neg = [NH|NT],
        append(L, [H], List)
        );
        (
        (L \= [] ->
            (append(G, [L], G1),
            List = []);
            (List = [],
            G1 = G)),
            Neg = NT,
        append(L1, [H], List1)
        )),
    groupExamples(T,Neg,List,List1,G1,G2).


abduceExamples([],_,Abd,Abd).
abduceExamples([H|T],PrevIntersect,Abd,Abd1):-
    epsilon(Eps),
    std(H,Mean1,_,Std1),
    ([Next|_] = T->(
    std(Next,Mean2,_,Std2),
    StartingPoint is (Mean1 + Mean2) / 2 ,
    intersect(StartingPoint,Mean1,Mean2,Std1,Std2, Eps,Intersect),
    PrevIntersect1 is Intersect + Eps,
    (nonvar(PrevIntersect)->(append([PrevIntersect],H,Temp),append(Temp,[Intersect],Append));
                           append(H,[Intersect],Append)));append([PrevIntersect],H,Append)),
    append(Abd,[Append],Abd2),
    abduceExamples(T,PrevIntersect1,Abd2,Abd1).

memberOfList(A,B):-
    append(A, B, C),
    sort(C,D),
    B=D.

groupCoverage(H,TotalCoverage,Coverage,Len):-
    min_list(H,Min),
    max_list(H,Max),
    epsilon(Eps),
    length(H,Len),
    Coverage is ((Max - Min) / Eps) / TotalCoverage.

totalCoverage(NestedList,Eps,TotalCoverage,TotLength):-
    flatten(NestedList, FlatList),
    min_list(FlatList,Min),
    max_list(FlatList,Max),
    length(FlatList,TotLength),
    TotalCoverage is (Max - Min)/Eps.

generalise([],Prev,_,_,Gen,Gen1):-
    var(Prev)->Gen1=Gen;Gen1=Prev.
generalise([H|T],Prev,TotCov,TotLength,Gen,Gen1):-
        combination(Comb),
        groupCoverage(H,TotCov,Coverage,Len),
        ((Coverage<Comb,(Len/TotLength) < Comb)->
            (nonvar(Prev)->append(Prev,H,NewPrev),Gen2=Gen;NewPrev=H,Gen2=Gen);
            ((nonvar(Prev),groupCoverage(Prev,TotCov,PCoverage,PLen),(PCoverage>Comb,PLen/TotLength>Comb))->
            (append(Gen,[Prev],NewGen),append(NewGen,[H],Gen2),NewPrev=_)
            ;(append(Gen,[H],Gen2),NewPrev=_))),
        generalise(T,NewPrev,TotCov,TotLength,Gen2,Gen1).


    defineBranches([],L,L).
defineBranches([H1|T1],L1,L_1):-

    
    (inRange(H1,Val2±Range)->
    append(L1,[inRange(_,Val2±Range)],LL1);LL1=L1),

    (geq(H1,Val1)->
    append(LL1,[geq(_,Val1)],LL2);LL2=LL1),

    (leq(H1,Val)->
    append(LL2,[leq(_,Val)],LL3);LL3=LL2),

    defineBranches(T1,LL3,L_1).
%-----------------------------------Learning process--------------------------------------------------------

nscore(_,[],_,_,Score,Score).
nscore([],[_|T],Backup,MinusScore,Score,Score2):-
    (MinusScore=0 -> Score1 is Score + 1;Score1=Score),
    nscore(Backup,T,Backup,0,Score1,Score2).
nscore([PH|PT],[H|T],Backup,MinusScore,Score,Score2):-
    PH =.. [F,_,V],
    (call(F,H,V)->MinusScore1 is MinusScore + 1;MinusScore1 = MinusScore),
    nscore(PT,[H|T],Backup,MinusScore1,Score,Score2).


score(_,[],_,Score,Score).
score([],[_|T],Backup,Score,Score2):-
    score(Backup,T,Backup,Score,Score2).
score([PH|PT],[H|T],Backup,Score,Score2):-
    PH =.. [F,_,V],
    (call(F,H,V)->Score1 is Score + 1,score([],[H|T],Backup,Score1,Score2);
    Score1=Score,score(PT,[H|T],Backup,Score1,Score2)).

coverage(Prev,Current):-
    nonvar(Prev),
    Prev =.. [F,_,Var],
    Current =.. [F1,_,Var1],
    (number((Var1))->
    AVal=Var1;
    Var1= AVal±_
    ),
    (number((Var))->
    AVal2=Var;
    Var= AVal2±_
    ),
    Check=.. [F,AVal,Var],
    Check1=.. [F1,AVal2,Var1],
    %(retract(user:Check);true),(retract(user:Check1);true),
    (Check;Check1).


checkCoverage([],_).
checkCoverage([H|T],Prev):-
    not(coverage(Prev,H))->once(checkCoverage(T,H));fail,!.

prove(Rule,Pos,Neg,Pscore,Nscore):-
    score(Rule,Pos,Rule,0,Pscore),
    nscore(Rule,Neg,Rule,0,0,Nscore).

learnAndProves(Branches,Len,B,Pos,Neg,Pscore,Nscore):-
    hypothesisSpace(Branches,B),
    prove(B,Pos,Neg,Pscore,Nscore),
    length(B,L),L=<Len,checkCoverage(B,_).


learnAndScore(Branches,Len,Pos,Neg,Bag):-
    findall([B,Pscore,Nscore], learnAndProves(Branches,Len,B,Pos,Neg,Pscore,Nscore), Bag).

    showBestHypothesis([],_,Acc,Rule):-pprint(Rule,Acc).
showBestHypothesis([[Rule,TP,TN]|T],Len,Acc,Rule2):-
    Acc1 is (TP+TN) / Len,
    (Acc1 > Acc -> Rule1=Rule,Acc2=Acc1;Rule1=Rule2,Acc2=Acc),
    showBestHypothesis(T,Len,Acc2,Rule1).

showBestHypothesis([],_,Acc,Rule,Acc,Rule).
showBestHypothesis([[Rule,TP,TN]|T],Len,Acc,Rule2,Acc3,Rule3):-
    Acc1 is (TP+TN) / Len,
    (Acc1 > Acc -> Rule1=Rule,Acc2=Acc1;Rule1=Rule2,Acc2=Acc),
    showBestHypothesis(T,Len,Acc2,Rule1,Acc3,Rule3).

learn(Pos,Neg,Rule,Acc,Print):-
    append(Pos, Neg, AllExamples),
    sort(AllExamples, Sorted),
    setEpsilon(Sorted),
    epsilon(Eps),
    sort(Pos, PosSorted),%Just in case for denoising data 
    sort(Neg, NegSorted),
    groupExamples(Sorted,NegSorted,[],[],[],G),
    abduceExamples(G,_,[],Abd),
    totalCoverage(Abd,Eps,TotCov,TotLen),
    generalise(Abd,_,TotCov,TotLen,[],Generalised),
    write('\nCombining values...'),
    %write('Generating hypothesis space...'),
    length(Generalised,Len),
    defineBranches(Generalised,[],Branched),
    learnAndScore(Branched,Len,PosSorted,NegSorted,BagRules),
    writeln('\n---------------------------'),
    length(AllExamples,LenExamples),
    (Print=print->
        showBestHypothesis(BagRules,LenExamples,0,_);
        showBestHypothesis(BagRules,LenExamples,0,_,Acc,Rule)).

%----------------------------------------Find best rule----------------------------------------

combinations(0, _, []).
combinations(N, [X|Xs], [X|Ys]) :-
    N > 0,
    N1 is N - 1,
    combinations(N1, Xs, Ys).
combinations(N, [_|Xs], Ys) :-
    N > 0,
    combinations(N, Xs, Ys).

partition([], []).
partition(List, [Head | TailPartition]) :-
    append(Head, Tail, List),
    Head \= [],
    partition(Tail, TailPartition).

convertToRule([], fail). % No elements case
convertToRule([Conj], Rule) :- conjunction(Conj, Rule). % Single list case
convertToRule([Conj | Rest], (Rule ; RestRule)) :- 
    conjunction(Conj,Rule),
    convertToRule(Rest, RestRule).

% Helper predicate to convert each inner list to a conjunction
conjunction([], true).
conjunction([H], (Chain1,Chain2)):-  
    H=..[R,Feature,B],
    Chain1 =..[Feature,_,C],
    Chain2 =..[R,C,B].
    
conjunction([H | T], (Chain1,Chain2 , Rest)) :-
    H=..[R,Feature,B],
    Chain1 =..[Feature,_,C],
    Chain2 =..[R,C,B],
    conjunction(T, Rest).

callToCheckPos(Rule,Example,Y):-
    term_string(Rule, String),
    term_string(NewRule,String),
    free_variables(NewRule, Variables),
    nth0(0,Variables,D),
    D=Example,
    (NewRule)->
   Y=1;Y=0.

callToCheckNeg(Rule,Example,Y):-
    term_string(Rule, String),
    term_string(NewRule,String),
    free_variables(NewRule, Variables),
    nth0(0,Variables,D),
    D=Example,
    not((NewRule))->
   Y=1;Y=0.

getAccuracy(Rule,Pos,Neg,Tot,R,Acc):-
       convertToRule(Rule,R),
       maplist(callToCheckPos(R),Pos,L),
       maplist(callToCheckNeg(R),Neg,L1),
       sum_list(L, Sum),
       sum_list(L1, Sum1),
       Acc is (Sum+Sum1)/Tot.



getTheBestRule([],LAcc,LRule):-
    format('Accuracy:~w \n ~w',[LAcc,LRule]).
getTheBestRule([R=>Acc|T],LAcc,LRule):-
    (Acc>LAcc ->
    Acc2 = Acc,R2=R ; Acc2=LAcc,R2=LRule),
    getTheBestRule(T,Acc2,R2).
    
%-------------------------------------------------------------Multi-calss learning-----------------
appendWithVar([],_,NewArr,NewArr).
appendWithVar([H|T],R,NewArr,NewArr1):-
        H=..[F,_,V],
        P=..[F,R,V],
        append(NewArr,[P],NewArr2),
        appendWithVar(T,R,NewArr2,NewArr1).
    
batchLearn([],_,AllRules,AllRules).    
batchLearn([R=>APos|T],ArrNeg,AllRules,AllRules1):-
        format('Analysing feature value:~w\n',[R]),
        member(R=>ANeg,ArrNeg),
        learn(APos,ANeg,Rule,_,give),
        %writeln(R),writeln(Rule),
        appendWithVar(Rule,R,AllRules,NewArr),
        batchLearn(T,ArrNeg,NewArr,AllRules1).

clausesInFile(File,H1,A,V) :-
    absolute_file_name(File, AbsFileName),predicate_property(H, file(AbsFileName)),
    clause(H,true),
    H =..[H1,A,V].


replaceP(_, _, [], []).
replaceP(O, R, [O|T], [R|T2]) :- replaceP(O, R, T, T2).
replaceP(O, R, [H|T], [H|T2]) :- dif(H,O), replaceP(O, R, T, T2).

addToExamples([],Arr,Arr).
addToExamples([[A,V]|T],Arr,Arr1):-
  (member(A=>B,Arr)->
    (append(B,[V],C),
    replaceP(A=>B,A=>C,Arr,Arr2));
    append(Arr,[A=>[V]],Arr2)),
    addToExamples(T,Arr2,Arr1).

getAllValues(_,[],Ex,Ex).
getAllValues(File,[H|T],Ex,Ex1):-
    findall([Head,Val],clausesInFile(File,Head,H,Val),All),
    addToExamples(All,Ex,FinalEx),
    getAllValues(File,T,FinalEx,Ex1).



generateInterests(BKFile,ExampleFile,Rules,ArrPos,ArrNeg,Tot):-
     %consult(BKFile),
    %consult(ExampleFile),
    [BKFile],
    [ExampleFile],
    writeln('Processing file and numerical values...'),
    findall(A,(pos(V),V=..[_,A]),ArrPos),
    findall(A,(neg(V),V=..[_,A]),ArrNeg),
   
    getAllValues(BKFile,ArrPos,[],AllPosValues),
    getAllValues(BKFile,ArrNeg,[],AllNegValues),
    length(ArrPos,LenPos),length(ArrNeg,LenNeg),
      
   
    batchLearn(AllPosValues,AllNegValues,[],Rules),
    Tot is LenPos+LenNeg,

    length(Rules,LRu),
    format('Number of interesting features:~w \nGenerating hypothesis space and looking for a best rule...',[LRu]).


learnFromFile(BKFile,ExampleFile):-
    generateInterests(BKFile,ExampleFile,Rules,ArrPos,ArrNeg,Tot),

    maxClause(MaxClause),
    
    findall(R => Acc, 
            (between(1, MaxClause, N), combinations(N, Rules, Comb),partition(Comb,Comb1),getAccuracy(Comb1,ArrPos,ArrNeg,Tot,R,Acc)), 
            Result),
            
    length(Result,LR),
    format('Hypothesis space:~w \n',[LR]),
    getTheBestRule(Result,0,[]).

    %writeln('Generating Final Rule, this may take from few seconds to few minutes...'),
    %searchForTheBestRule(Rules,ArrPos,ArrNeg,LenPos,LenNeg),
    %numlogrule(BestSubList),
    %numlogacc(MaxSum),
    %write('Final best rule: '),convertToRule(BestSubList,_Atom,[],[],R), 
    %format('rule(~w):-~w\n',[_Atom,R]),
    %writeln(R),
    %write('Accuracy: '), writeln(MaxSum),
    %removeAndAdd(0,[]).
    

%in the futur add combination of disjunction rules using(memebr(x,list1), ...) 