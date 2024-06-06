%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Numlog: numerical learning and reasoning using ILP %%%%%%%%%%%%%%%%%%%%%%%%
%                                     Numlog Author: Daniel Cyrus                                       %
%                                                                                                       % 
%                                            The HMLR lab                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


:- module(numlog, [learn/2,learn/3,learn/4,user:eq/2,user:leq/2,user:geq/2,user:inRange/2,user:epsilon/1,max_clauses/1]).
%configuration
user:epsilon(0.0001).
max_clauses(4).

%Background knowledge
round(D,X,Y) :- Z is X * 10^D, round(Z, ZA), Y is ZA / 10^D.

eq([A|[]],A).

leq(List,Last):-
    last(List, Last).
    %round(Last1,Last,2).

geq([H|_],H).%:-round(H,H1,2).

inRange(List,A):-
    [First|_] = List, 
    last(List,Last),
    Mid is (Last + First) / 2,
    Range is Last - Mid,
    concat(Mid, '±', Atom),
    concat(Atom, Range, A).

%evaluation 
user:eq(A,A).
user:leq(A,B):-A=<B.
user:geq(A,B):-A>=B.
user:inRange(Val,A):-
    split_string(A, '±','',[AVal,Range]),
    number_string(A_Val, AVal),
    number_string(A_Range, Range),
    Val >= A_Val - A_Range,
    Val =< A_Val + A_Range.
  


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
%///////////////////////////////////////////////////Threshold Analysis//////////////////

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
    (user:epsilon(Eps),
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

%///////////////////////////////////threshold analysis////////////////////

%generating hypothesises using Top-down approach 
%                          Root(leq, eq, [])
%                                   |
%             ------------------------------------------------
%            |                    |                          |
%          eq                  inRange                       geq
%           |                     |                          |
%     ---------------     --------------------      ------------------
%    |      |       |     |       |        |        |       |        |   
%    eq   inRange  geq   eq     inRange   geq      eq     inRange   geq 

defineRoot([H|_],Root):-
    (length(H,0),
    true);
    ((length(H,1),
    eq(H,Val),
    Atom = eq);
    (
    Atom = leq,
    leq(H,Val))),
    Root =.. [Atom,_,Val].

coverage(Prev,Val):-
    nonvar(Prev),
    Prev =.. [F,_,Var],
    (number(Val)->
    AVal=Val;
    (split_string(Val, '±','',[AVal1,_]),
    number_string(AVal, AVal1))
    ),
    Check=.. [F,AVal,Var],
    call(user:Check).

prove(_,[],_,_).
prove([],[_|T],Prev,Temp):-
    prove(Temp,T,Prev,Temp).
prove([H|T],[EH|ET],Prev,Temp):-
        H =.. [Func,_,Val],
        Var = EH,
        %atom_concat(user:Func, Atom),
        Check =.. [Func,Var,Val],
        ((not(call(user:Check)), not(coverage(Prev,Val)))->
        (Prev1=Check,
        prove(T,[EH|ET],Prev1,Temp));
        (fail,!)).     


pprint(H,P):-
    H =.. [F,_,V],
    numbervars(A, 0, _),
    H1 =.. [F,A,V],
    P1 is P * 100,
    P2 is float_integer_part(P1),
    ansi_format([italic,fg(blue)],'~nprobability:~w% ,Entropy:0.0, Gain:0.0~n',[P2]),
    format('value(~w):-~n ~40t ~w.~n~n',[A,H1]).

pprint(H):-
    H =.. [F,_,V],
    numbervars(A, 0, _),
    H1 =.. [F,A,V],
    format('~w.~n',[H1]).


%///////////learning with new values/////////////////////////
isEmpty(L):-length(L, 0).
groupExamplesAbc([],_,TempPos,TempNeg,Group,FinalGroup,_,_,TempP,Ptot):-
    (not(isEmpty(TempPos)),length(TempPos, LenPos),append(TempP, [LenPos], Ptot),append(Group,[TempPos],FinalGroup));
    (not(isEmpty(TempNeg)),length(TempNeg, LenNeg),append(TempP, [LenNeg], Ptot),append(Group,[TempNeg],FinalGroup)).
groupExamplesAbc([H|T],Neg,TempPos,TempNeg,Group,FinalGroup,PP,PN,TempP,Ptot):-
    (((Neg = [NH|_],H\=NH); isEmpty(Neg)) -> 
    ((not(isEmpty(TempNeg))->(append(Group,[TempNeg],Group1),append(TempP, [PN], TempP1));(Group1=Group,TempP1=TempP)),
        append(TempPos,[H],TempPos1),
        PP1 is PP + 1,
        PN1 = 0,
    TempNeg1 = [],
    Neg1 = Neg
    );
    (
        Neg = [NH|NT],
        (not(isEmpty(TempPos))->(append(Group,[TempPos],Group1),append(TempP, [PP], TempP1));(Group1=Group,TempP1=TempP)),
        append(TempNeg,[NH],TempNeg1),
        TempPos1 = [],
        PN1 is PN + 1,
        PP1 = 0,
        Neg1 = NT
    )),
    groupExamplesAbc(T,Neg1,TempPos1,TempNeg1,Group1,FinalGroup,PP1,PN1,TempP1,Ptot).

memberOf([],_).
memberOf([H|T],Neg):-
  member(H,Neg),
  memberOf(T,Neg).

abductionProcess([],_,FinalSamples,_,FinalSamples).
abductionProcess([H|T],Neg,NewSamples,Prev,FinalSamples):-
    ([H1|_] = T ->
    (
    std(H,Mean1,_,Std1),
    std(H1,Mean2,_,Std2),
    StartingPoint is (Mean1 + Mean2) / 2,
    
    user:epsilon(Eps),
    intersect(StartingPoint,Mean1,Mean2,Std1,Std2,Eps,Intersect),
        (once(memberOf(H, Neg))->(Prev1 = Intersect,append(NewSamples,[[]],NewSamples1));
                                ((nonvar(Prev)->(append([Prev],H,Temp),append(Temp,[Intersect],Range));append(H,[Intersect],Range)),
                                 append(NewSamples,[Range],NewSamples1),Prev1 = _))
    );
    (
        not(memberOf(H, Neg)) ->(
        ((nonvar(Prev))->append([Prev],H,Range);Range=H),
        append(NewSamples,[Range],NewSamples1));
        append(NewSamples,[[]],NewSamples1)

    )),
    abductionProcess(T,Neg,NewSamples1,Prev1,FinalSamples).

probabilities([],[],_,Prob,Prob).
probabilities([H|T],[H1|T1],SamplesLen,Temp,Prob):-
    (H1 \= [] ->
    (P is H / SamplesLen,
    append(Temp,[P],Temp1));Temp1 = Temp),
    probabilities(T,T1,SamplesLen,Temp1,Prob).

%Abduction learn
defineBranches([],L,L).
defineBranches([H1|T1],L1,L_1):-
    (leq(H1,Val)->
    append(L1,[leq(_,Val)],LL1);LL1=L1),
    
    (geq(H1,Val1)->
    append(LL1,[geq(_,Val1)],LL2);LL2=LL1),
    
    (inRange(H1,Val2)->
    append(LL2,[inRange(_,Val2)],LL3);LL3=LL2),

    defineBranches(T1,LL3,L_1).


hypothesisSpace(0,[],[]).
hypothesisSpace(_,[],[]).
%hypothesisSpace(0,_,[]).
hypothesisSpace(N,[X|T],[X|Comb]) :-
    N>0,
    N1 is N-1,
    hypothesisSpace(N1,T,Comb).
hypothesisSpace(N,[_|T],Comb) :-
    N>0,
    hypothesisSpace(N,T,Comb).

%hypothesisSpace([],[]).
%hypothesisSpace([H|T],[H|T2]) :-
%    hypothesisSpace(T,T2).
%hypothesisSpace([_|T],T2) :-
%    hypothesisSpace(T,T2).

learn(Pos,Neg):-
    max_clauses(Max),
    %(length(L1,L),L > Max  -> hypothesisSpace(Max,L1,Rule);
    hypothesisSpace(Max,Pos,Rule),
    %length(Rule,L),
    %L =< Max,
    prove(Rule,Neg,_,Rule),
    maplist(pprint,Rule).


learn(Pos,Neg,normal):-
    append(Pos, Neg, AllExamples),
    %denois(Pos,Neg,NewPos),%denois the positive values
    sort(AllExamples, Sorted),%sort all values to define the numerical range and remove noises -> [1,2,1] -> [1,2].
    %decrease search space by deviding examples into categories.
    %sort(Pos, PosSorted),
    sort(Neg, NegSorted),
    groupExamples(Sorted,NegSorted,[],[],G),
    defineBranches(G,[],L1),
    learn(L1,NegSorted).

learn(Pos,Neg,threshold):-
    %maplist(round(3),Pos,RPos),
    %maplist(round(3),Neg,RNeg),
    append(Pos, Neg, AllExamples),
    %denois(Pos,Neg,NewPos),%denoising the positive values
    sort(AllExamples, Sorted),%sorting all values to define the numerical range and remove noises -> [1,2,1] -> [1,2].
    %decrease search space by deviding examples into categories.
    sort(Neg, NegSorted),
    
    once(groupExamplesAbc(Sorted,NegSorted,[],[],[],G,0,0,[],_)),%last argument is probablity, add this feature in the future
    once(abductionProcess(G,NegSorted,[],_,G1)),
    %length(AllExamples,SamplesLen),
    %probabilities(Prob,G1,SamplesLen,[],Probabilities),
    %defineRoot(G,Root), Removing defineRoot for now, as using it doesn't make any changes on search complexity
    defineBranches(G1,[],L1),
    learn(L1,NegSorted).

learn(Pos,Neg,Method,epsilon(Eps)):-
    abolish(user:epsilon/1),
    asserta(user:epsilon(Eps)),
    learn(Pos,Neg,Method).

learn(Pos,Neg,Method,Round):-
    maplist(round(Round),Pos,RPos),
    maplist(round(Round),Neg,RNeg),
    learn(RPos,RNeg,Method).


    
    
%[[1.4925316,1.492560960500352],[],[1.5011797035013035,1.532899086,1.532928830000581],[],[1.5470872870006538,1.558849459,1.567118541,1.5693316882505075],[]]
%[leq(_57272,1.4927909604999734),inRange(_57290,'1.5172777667499264±0.015873063250028796'),inRange(_57302,'1.5581679876249535±0.011044700625005177')]
