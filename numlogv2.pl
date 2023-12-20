%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Numlog: numerical learning and reasoning using ILP %%%%%%%%%%%%%%%%%%%%%%%%
%                                     Numlog Author: Daniel Cyrus                                       %
%                                     Pygol Author: Dany Varghese                                       % 
%                                            The HMLR lab                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

:- module(numlog, [learn/2,learn/3,user:eq/2,user:leq/2,user:geq/2,user:inRange/2]).
%configuration
epsilon(0.0001).

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
    (Len =< 1 ->
    (epsilon(Eps),
    [H|_] = A,
    Min is H - Eps,
    Max is H + Eps,
    linespace(Min,Max,10,Aall),
    Len1 = 9
        );
    (
        
        Len >= 10 ->(
            Aall = A,
        Len1 is Len - 1);
        (
        min_list(A, Min),
        max_list(A, Max),
        linespace(Min,Max,10,Aall),
        Len1 = 9)
    )),
    sigma(Aall,Mean,0,Sigma),
    %Len1 is Len - 1,
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


intersect(StartingPoint,Mean1,Mean2,Std1,Std2, Increment,Intersect):-
    pdf(StartingPoint,Mean1,Std1,PDF),
    pdf(StartingPoint,Mean2,Std2,PDF1),
    Ans is PDF - PDF1,
    StartingPoint1 is StartingPoint + Increment,
    (Ans > 0 ->
    intersect(StartingPoint1,Mean1,Mean2,Std1,Std2, Increment,Intersect);
    Intersect = StartingPoint1).

%///////////////////////////////////End of threshold analysis////////////////////

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
        %atom_concat(user:Func, Atom),
        Check =.. [Func,Var,Val],
        (not(call(user:Check))->
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


pprint(H,P):-
    H =.. [F,_,V],
    numbervars(A, 0, _),
    H1 =.. [F,A,V],
    P1 is P * 100,
    P2 is float_integer_part(P1),
    ansi_format([italic,fg(blue)],'~nprobability:~w% ~n',[P2]),
    format('value(~w):-~n ~40t ~w.~n~n',[A,H1]).

denois(Pos,[],Pos).
denois(Pos,[H|T],DenoisedPos):-
    delete(Pos, H, DenoisedPos1),
    denois(DenoisedPos1,T,DenoisedPos).
    

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
    intersect(StartingPoint,Mean1,Mean2,Std1,Std2,0.0001,Intersect),

        (once(memberOf(H, Neg))->(Prev1 = Intersect,append(NewSamples,[[]],NewSamples1));
                                ((nonvar(Prev)->(append([Prev],H,Temp),append(Temp,[Intersect],Range));append(H,[Intersect],Range)),
                                 append(NewSamples,[Range],NewSamples1),Prev1 = _))
    );
    (
        (nonvar(Prev)->append([Prev],H,Range);Range=H),
        append(NewSamples,[Range],NewSamples1)
    )),
    abductionProcess(T,Neg,NewSamples1,Prev1,FinalSamples).

probabilities([],[],_,Prob,Prob).
probabilities([H|T],[H1|T1],SamplesLen,Temp,Prob):-
    (H1 \= [] ->
    (P is H / SamplesLen,
    append(Temp,[P],Temp1));Temp1 = Temp),
    probabilities(T,T1,SamplesLen,Temp1,Prob).

%Abduction learn

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
    maplist(pprint,Rule,0).

learn(Pos,Neg,normal):-
    learn(Pos,Neg).
learn(Pos,Neg,threshold):-
    append(Pos, Neg, AllExamples),
    %denois(Pos,Neg,NewPos),%denoising the positive values
    sort(AllExamples, Sorted),%sorting all values to define the numerical range and remove noises -> [1,2,1] -> [1,2].
    %decrease search space by deviding examples into categories.
    %sort(Pos, PosSorted),
    sort(Neg, NegSorted),
    groupExamplesAbc(Sorted,NegSorted,[],[],[],G,0,0,[],Prob),
    abductionProcess(G,NegSorted,[],_,G1),
    length(AllExamples,SamplesLen),
    probabilities(Prob,G1,SamplesLen,[],Probabilities),
    %defineRoot(G,Root), Removing defineRoot for now, as using it doesn't make any changes on search complexity
    findall(Preds,defineBranches(G1,[],Preds),AllPreds),
    %write(AllPreds),
    proveAll(AllPreds,NegSorted,[],Rules),
    mostGeneral(Rules,0,_,Rule),
    maplist(pprint,Rule,Probabilities).

%learn(Pos,Neg,Header,threshold):-
