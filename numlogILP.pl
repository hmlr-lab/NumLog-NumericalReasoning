%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Numlog: numerical learning and reasoning using ILP %%%%%%%%%%%%%%%%%%%%%%%%
%                                     Numlog Author: Daniel Cyrus                                       %
%                            This version is desined for first order logic                              % 
%                                            The HMLR lab                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
:- module(numlogILP, [learn/2,learn/3,learnFromFile/2,testAccuracy/3]).

:-op(500, xfy, user:(±)). 
:-op(500, xfy, =>).


:-dynamic(user:pos/1).
:-dynamic(user:neg/1).

:-dynamic epsilon/1.
:-dynamic accCounter/1.

%epsilon(0.001).
combination(0.05).

round(D,X,Y) :- Z is X * 10^D, round(Z, ZA), Y is ZA / 10^D.

groupExamples([],_,Temp,Group,Group1):-
    append(Group,[Temp],Group1).

groupExamples(Values,[],_,Group,Group1):-
    append(Group,[Values],Group1).

groupExamples([H|T],[NH|NT],Temp,Group,Group2):-
    (H \= NH ->
    (
        append(Temp,[H],Temp1),
        Nexamples = [NH|NT],
        Group1 = Group
    );
    (

        append(Group,[Temp],Group1),
        Temp1 = [],
        Nexamples = NT
    )
    ),
    groupExamples(T,Nexamples,Temp1,Group1,Group2).
%///////////////////////////////////////////////////

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

%/////////////////////////////Generalisation process//////////////////////
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
        ((Coverage<Comb,Len/TotLength < Comb)->
            (nonvar(Prev)->append(Prev,H,NewPrev),Gen2=Gen;NewPrev=H,Gen2=Gen);
            ((nonvar(Prev),groupCoverage(Prev,TotCov,PCoverage,PLen),PCoverage>Comb,PLen/TotLength>Comb)->
            (append(Gen,[Prev],NewGen),append(NewGen,[H],Gen2),NewPrev=_)
            ;(append(Gen,[H],Gen2),NewPrev=_))),
        generalise(T,NewPrev,TotCov,TotLength,Gen2,Gen1).


count(_,[],C,C).
count(A,[H|T],C,C2):-
    (member(H,A)-> C1 is C + 1;C1 = C),
    count(A,T,C1,C2).

reducePhrases([],_,_,Final,Final).
reducePhrases([H|T],Pos,Neg,Final,Final1):-
    count(H,Pos,0,CPos),
    count(H,Neg,0,CNeg),
    (CPos >= CNeg -> append(Final,[H],Final2);Final2=Final),
    reducePhrases(T,Pos,Neg,Final2,Final1).

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
    
/*
defineBranches([],L,L).
defineBranches([H1|T1],L1,L2):-
    inRange(H1,Val2),
    append(L1,[inRange(_,Val2)],L3),
    defineBranches(T1,L3,L2).

defineBranches(Branches,L):-
    [H|_] = Branches,
    last(Branches,T),
    defineBranches(Branches,[],LL),
    leq(H,Val),
    append([leq(_,Val)],LL,L1),
    geq(T,Val1),
    append(L1,[geq(_,Val1)],L).
*/

defineBranches([],L,L).
defineBranches([H1|T1],L1,L_1):-

    (leq(H1,Val)->
    append(L1,[leq(_,Val)],LL1);LL1=L1),
    
    (geq(H1,Val1)->
    append(LL1,[geq(_,Val1)],LL2);LL2=LL1),
    
    (inRange(H1,Val2±Range)->
    append(LL2,[inRange(_,Val2±Range)],LL3);LL3=LL2),

    defineBranches(T1,LL3,L_1).


find_unique_n(N, X, Goal, Xs) :-
    find_unique_n(N, X, Goal, Xs, []).

find_unique_n(N, X, Goal, Xs, Sols) :-
    N > 0,
    copy_term(X-Goal, CX-CGoal),
    call(CGoal),
    \+ (member(Sol,Sols), variant(Sol,CX)),
    !,
    N1 is N-1,
    Xs = [CX|Xs1],
    Sols1 = [CX|Sols],
    find_unique_n(N1, X, Goal, Xs1, Sols1).
find_unique_n(_N, _X, _Goal, [], _Sols).

hypothesisSpace([],[]).
hypothesisSpace([X|T],[X|Comb]):-
    hypothesisSpace(T,Comb).
hypothesisSpace([_|T],Comb) :-
    hypothesisSpace(T,Comb).


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
    %once(learnAndProves(Branches,Len,B,Pos,Neg,Pscore,Nscore)),
    %Bag = [B,Pscore,Nscore].

pprint([],Acc):-writeln('--------------------------------'),format('Accuracy: ~2f~n',[Acc]).    
pprint([H|T],Acc):-
    H =.. [F,_,V],
    numbervars(A, 0, _),
    H1 =.. [F,A,V],
    format('~w.~n',[H1]),
    pprint(T,Acc).

pprint([],_,_).%:-writeln('--------------------------------'),format('Accuracy: ~2f~n',[Acc]).    
pprint([H|T],Inc,_):-
    H =.. [F,P,V],
    numbervars(A, Inc, _),
    numbervars(B, 0, _),
    H1 =.. [F,A,V],
    H2 =.. [P,B,A],
    format('~w,~w,\n',[H2,H1]),
    Inc2 is Inc + 1,
    pprint(T,Inc2,_).

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

%////////////////////////Epsilon calculation///////////////////
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
  %EPS1 is EPS / 10,
  %round(5,EPS1,EPS2),
  assert(epsilon(EPS)),
  format('\nepsilon has set to ~f \n',[EPS]).
%//////////////////////////////////////////////////////////////

learn(Pos,Neg,Rule,Acc,Print):-
     
    write('Analysing numbers...'),
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
    write('Compressing values...'),
    reducePhrases(Generalised,Pos,Neg,[],FinalGeneralised),
    %write('Generating hypothesis space...'),
    length(FinalGeneralised,Len),
    defineBranches(FinalGeneralised,[],Branched),

    learnAndScore(Branched,Len,PosSorted,NegSorted,BagRules),
    writeln('\n---------------------------'),
    length(AllExamples,LenExamples),
    (Print=print->
        showBestHypothesis(BagRules,LenExamples,0,_);
        showBestHypothesis(BagRules,LenExamples,0,_,Acc,Rule)).

learn(Pos,Neg):-
    learn(Pos,Neg,_,_,print).

learn(Pos,Neg,Rule):-
        learn(Pos,Neg,Rule,_,get).
%///////////////Multi class learning////////////////////
clausesInFile(File,H1,A,V) :-
    absolute_file_name(File, AbsFileName),
    predicate_property(H, file(AbsFileName)),
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

appendWithVar([],_,AllRule,AllRule).
appendWithVar([H|T],R,AllRule,AllRule1):-
    H=..[F,R,V],
    P=..[F,R,V],
    append(AllRule,[P],AllRule2),
    appendWithVar(T,R,AllRule2,AllRule1).

batchLearn([],_,AllRules,AllRules).    
batchLearn([R=>APos|T],ArrNeg,AllRules,AllRules1):-
    member(R=>ANeg,ArrNeg),
    learn(APos,ANeg,Rule),
    %writeln(R),writeln(Rule),
    appendWithVar(Rule,R,AllRules,AllRules2),
    batchLearn(T,ArrNeg,AllRules2,AllRules1).

combinations(0, _, []).
combinations(N, List, [H|T]) :-
    N > 0,
    select(H, List, Rest),
    N1 is N - 1,
    combinations(N1, Rest, T).

func(A,B):-A=..[_,B,_].
check(A):-
    maplist(func,A,B),
    sort(B,C),
    B=C.


%////////////////////////////////////////////calculate accuracy///////////////////////

generate_string(Length, String) :-
    length(CharList, Length),
    maplist(=('~w'), CharList),
    atomic_list_concat(CharList, ',', String).
  
  convert_to_rule(Head, BodyPredicates, Term) :-
      retractall(cancer/1),
      length(BodyPredicates,L),
      generate_string(L,S),
      format(atom(Body), S, BodyPredicates), % Create the body as a conjunction of predicates
      format(atom(Rule), '~w :- ~w', [Head, Body]),% Combine head and body into a rule
      term_to_atom(Term, Rule). 
     
  
 % Initialize the counter
init_counter :-
    retractall(accCounter(_)),
    assert(accCounter(0)).
  
  % Increment the counter
  increment_counter :-
    accCounter(C),
    C1 is C + 1,
    retract(accCounter(C)),
    assert(accCounter(C1)).
  
    % Example predicate to be applied
  checkTheRule(R,X) :-
      (call(R,X)->  
      increment_counter;
      true).
  
  checkTheRuleN(R,X) :-
      (not(call(R,X))->  
      increment_counter;
      true).
  
    % Apply the predicate to each element of the list and count the applications
    count_maplist_applications(Predicate, List,Count) :-
      init_counter,
      maplist(Predicate, List),
      accCounter(Count).
  
  convertToList([],_,L,L).
  convertToList([H|T],Var,List,L):-
      H =.. [F,P,V],
      H1 =.. [F,A,V],
      H2 =.. [P,Var,A],
      append([H2,H1],List,List1),
      convertToList(T,Var,List1,L).
  
checkAcc(Var,List,Rule):-
    List \= [],
    convertToList(List,Var,[],A1),
    convert_to_rule(cancer(Var),A1,Rule).
   
%/////////////////////////////////////////////////////////////////////////////////////

metaHypothesisSpace(List, MaxLength, Combinations) :-
    findall(Rule, (between(0, MaxLength, N), combinations(N, List, Combination),check(Combination),checkAcc(_S,Combination,Rule)), Combinations).


chooseTheBestRule([],Acc,Rule,_,_,_):-writeln(Rule),format('Accuracy:~2f',[Acc]).
chooseTheBestRule([H|T],Acc,Rule,Pos,Neg,Len):-
    assert(H),
    count_maplist_applications(checkTheRule(cancer),Pos,TP),
    count_maplist_applications(checkTheRuleN(cancer),Neg,TN),
    retract(H),
    Acc1 is (TP+TN)/Len,
    (Acc1>Acc-> Acc2=Acc1,Rule2=H;Acc2=Acc,Rule2=Rule),
    chooseTheBestRule(T,Acc2,Rule2,Pos,Neg,Len).

learnFromFile(BKFile,ExampleFile):-
    [BKFile],
    [ExampleFile],
    writeln('Processing file and numerical values...'),
    findall(A,(pos(V),V=..[_,A]),ArrPos),
    findall(A,(neg(V),V=..[_,A]),ArrNeg),
    getAllValues(BKFile,ArrPos,[],AllPosValues),
    getAllValues(BKFile,ArrNeg,[],AllNegValues),
    length(ArrPos,LenPos),length(ArrNeg,LenNeg),
    batchLearn(AllPosValues,AllNegValues,[],Rules),
    writeln('Generating Final Rule, this may take few seconds'),
    metaHypothesisSpace(Rules,4,Rule),
    length(Rule,L),
    format('~w possible hypothesises\n',[L]),
    AllLen is LenPos + LenNeg,
    chooseTheBestRule(Rule,0,_,ArrPos,ArrNeg,AllLen).
    
testAccuracy(BKFile,ExampleFile,Rule):-
    [BKFile],
    [ExampleFile],
    writeln('Testing NumLog...'),
    findall(A,(pos(V),V=..[_,A]),ArrPos),
    findall(A,(neg(V),V=..[_,A]),ArrNeg),
    assert(Rule),
    count_maplist_applications(checkTheRule(cancer),ArrPos,TP),
    count_maplist_applications(checkTheRuleN(cancer),ArrNeg,TN),
    retract(Rule),
    length(ArrPos,LenPos),length(ArrNeg,LenNeg),
    AllLen is LenPos + LenNeg,
    Acc1 is (TP+TN)/AllLen,
    format('Test accuracy:~2f',[Acc1]).

    

    



