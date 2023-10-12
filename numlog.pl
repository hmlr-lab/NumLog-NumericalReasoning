%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Numlog: numerical learning and reasoning using ILP %%%%%%%%%%%%%%%%%%%%%%%%
%                                     Numlog Author: Daniel Cyrus                                       %
%                                     Pygol Author: Dany Varghese                                       % 
%                                            The HMLR lab                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

:- module(numlog, [learn/2]).

%Background knowledge
eq([A|[]],A).

leq(List,Last):-
    last(List, Last).

geq([H|_],H).

inRange(List,A):-
    sum_list(List, Sum),
    length(List, Int),
    Avg is Sum / Int,
    min_list(List, Min),
    Range is Avg - Min,
    concat(Avg, '±', Atom),
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

prove(_,[]).
prove([H|T],[EH|ET]):-
    H =.. [Func,_,Val],
    Var = EH,
    atom_concat(e_, Func, Atom),
    Check =.. [Atom,Var,Val],
    (call(Check)->
        prove([H|T],ET);
        prove(T,[EH|ET])).
  
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


proveAll([],_,_,Rules,Rules).
proveAll([H|T],Pos,Neg,Rules,Rules2):-
    ((prove(H,Pos),prove(H,Neg,H))->
        append(Rules, [H],Rules1);
        Rules1 = Rules),
        proveAll(T,Pos,Neg,Rules1,Rules2).

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
    format('dim(A):-~n ~40t ~w.~n',[H1]).

learn(Pos,Neg):-
    append(Pos, Neg, AllExamples),
    sort(AllExamples, Sorted),%sort all values to define the numerical range and remove noises -> [1,2,1] -> [1,2].
    %decrease search space by deviding examples into categories.
    sort(Pos, PosSorted),
    sort(Neg, NegSorted),
    groupExamples(Sorted,NegSorted,[],[],G),
    %defineRoot(G,Root), Removing defineRoot for now, as using it doesn't make any changes on search complexity
    %(nonvar(Root)->([H|_]=G,(G, H, List2));List2=G),
    findall(Preds,defineBranches(G,[],Preds),AllPreds),
   
    proveAll(AllPreds,PosSorted,NegSorted,[],Rules),
   
    mostGeneral(Rules,0,_,Rule),
    maplist(pprint,Rule).



