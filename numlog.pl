%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Numlog: numerical learning and reasoning using ILP %%%%%%%%%%%%%%%%%%%%%%%%
%                                     Numlog Author: Daniel Cyrus                                       %
%                                     Pygol Author: Dany Varghese                                       % 
%                                            The HMLR lab                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

:- module(numlog, [learn/2,learn/3,e_eq/2,e_leq/2,e_geq/2,e_inRange/2]).

%Background knowledge
round(X,Y,D) :- Z is X * 10^D, round(Z, ZA), Y is ZA / 10^D.%disabled as two value may be closed to each other

eq([A|[]],A).

leq(List,Last):-
    last(List, Last).
    %round(Last1,Last,2).

geq([H|_],H).%:-round(H,H1,2).

inRange(List,A):-
    [First|_] = List, 
    last(List,Last),
    Mid is (Last - First) / 2,
    Range is Mid - Last,
    concat(Mid, '±', Atom),
    concat(Atom, Range, A).

%evaluation 
e_eq(A,A).
e_leq(A,B):-A=<B.
e_geq(A,B):-A>=B.
e_inRange(Val,A):-
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

defineBranches([],Branches,Branches).
defineBranches([H|T],Branches,Branches2):-
    (length(H,0)->
    Branches1 = Branches;
    (((length(H,1),
    eq(H,Val),
    Atom = eq);
    (Atom = leq,
    leq(H,Val));
    (Atom = geq,
    geq(H,Val));
    (Atom = inRange,
    inRange(H,Val))),
    Pred =.. [Atom,_,Val],
    append(Branches, [Pred], Branches1))),
    defineBranches(T,Branches1,Branches2).

  
prove(_,[],_).
prove([],[_|T],Temp):-
    prove(Temp,T,Temp).
prove([H|T],[EH|ET],Temp):-
        H =.. [Func,_,Val],
        Var = EH,
        atom_concat(e_, Func, Atom),
        Check =.. [Atom,Var,Val],
        (not(call(Check))->
            prove(T,[EH|ET],Temp);
            fail).    


proveAll([],_,Rules,Rules).
proveAll([H|T],Neg,Rules,Rules2):-
    (prove(H,Neg,H)->
        append(Rules, [H],Rules1);
        Rules1 = Rules),
        proveAll(T,Neg,Rules1,Rules2).

mostGeneral([],_,Rule,Rule).
mostGeneral([H|T],Len,Rule,Rule1):-
    length(H, Int),
    (Len < Int ->
    (Rule2 = H, L = Int);
    (Rule2 = Rule, L = Len)),
    mostGeneral(T,L,Rule2,Rule1).

mostSpecific([],_,Rule,Rule).
mostSpecific([H|T],Len,Rule,Rule1):-
    length(H, Int),
    (Len > Int ->
    (Rule2 = H, L = Int);
    (Rule2 = Rule, L = Len)),
    mostSpecific(T,L,Rule2,Rule1).  

pprint(H):-
    H =.. [F,_,V],
    numbervars(A, 0, _),
    H1 =.. [F,A,V],
    format('value(~w):-~n ~40t ~w.~n',[A,H1]).
denois(Pos,[],Pos).
denois(Pos,[H|T],DenoisedPos):-
    delete(Pos, H, DenoisedPos1),
    denois(DenoisedPos1,T,DenoisedPos).
    

learn(Pos,Neg):-
    append(Pos, Neg, AllExamples),
    %denois(Pos,Neg,NewPos),%denois the positive values
    sort(AllExamples, Sorted),%sort all values to define the numerical range and remove noises -> [1,2,1] -> [1,2].
    %decrease search space by deviding examples into categories.
    %sort(Pos, PosSorted),
    sort(Neg, NegSorted),
    groupExamples(Sorted,NegSorted,[],[],G),
    %defineRoot(G,Root), Removing defineRoot for now, as using it doesn't make any changes on search complexity
    %(nonvar(Root)->([H|_]=G,(G, H, List2));List2=G),
    findall(Preds,defineBranches(G,[],Preds),AllPreds),
    proveAll(AllPreds,NegSorted,[],Rules),
    mostGeneral(Rules,0,_,Rule),
    maplist(pprint,Rule).

%///////////learning with new values/////////////////////////
incVal( _  ,0.0, Val,Val).
incVal(Step,V  , Val,Val1):-
    Val2 is Val + Step,
    V1 is V - 1,
    incVal(Step,V1,Val2,Val1).

newVal(A,B,V,Val):-
    AA is A - B,
    abs(AA, Ac),
    Step is Ac / 10,
    V1 is V * 10,
    once(incVal(Step,V1,A,Val)).

groupExamplesAbduction([],_,_,Temp,Group,Group1):-
    append(Group,[Temp],Group1).

groupExamplesAbduction(Values,[],_,_,Group,Group1):-
    append(Group,[Values],Group1).

groupExamplesAbduction([H|T],[NH|NT],Range,Temp,Group,Group2):-
    (H \= NH ->
    (
        append(Temp,[H],Temp1),
        Nexamples = [NH|NT],
        Group1 = Group,
        TT = T
    );
    (
        %abduce from positive to negative
        (last(Temp, Last)->
        (
        newVal(Last,NH,Range,NV),
        append(Temp,[NV],G1))
        ;G1=Temp),
        %abduce from negative to positive
        
        ([TempH|_] = T, not(member(TempH, NT))->
        (newVal(NH,TempH,Range,NV1),
        append([NV1],T,TT))
        ;(TT = T)),
        append(Group,[G1],Group1),
        Temp1 = [],
        Nexamples = NT
    )
    ),
    groupExamplesAbduction(TT,Nexamples,Range,Temp1,Group1,Group2).



learn(Pos,Neg,Range):-
    append(Pos, Neg, AllExamples),
    %denois(Pos,Neg,NewPos),%denois the positive values
    sort(AllExamples, Sorted),%sort all values to define the numerical range and remove noises -> [1,2,1] -> [1,2].
    %decrease search space by deviding examples into categories.
    %sort(Pos, PosSorted),
    sort(Neg, NegSorted),
    groupExamplesAbduction(Sorted,NegSorted,Range,[],[],G),
    %defineRoot(G,Root), Removing defineRoot for now, as using it doesn't make any changes on search complexity
    findall(Preds,defineBranches(G,[],Preds),AllPreds),
    proveAll(AllPreds,NegSorted,[],Rules),
    mostGeneral(Rules,0,_,Rule),
    maplist(pprint,Rule).

