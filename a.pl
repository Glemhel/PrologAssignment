/**
 * Mikhail Rudakov
 * BS19-02
 * Prolog Programming Assignment
 * 
 **/

/**
 * map for the agent
 **/
map_xlimit(9).
map_ylimit(9).
covid([1, 0]).
covid([0, 1]).

member1([], _X) :- false.
member1([H | T], X) :-
     H = x; member1(T, X).

copyN(_L, 0, []) :- !.

copyN([], _ , []) :- !.

copyN(L, N, [L | L1]) :-
     N1 is N - 1,
     copyN(L, N1, L1).

get_adjacent([X, Y], L) :-
     Yup is Y + 1,
     Ydown is Y - 1,
     Xup is X + 1,
     Xdown is X - 1,
     L = [[Xup, Yup], [Xup, Y], [X, Yup], [Xdown, Yup], [Xup, Ydown], [Xdown, Y], [X, Ydown],
          [Xdown, Ydown]].


inside_map([X, Y]) :-
     map_xlimit(Xmax),
     map_ylimit(Ymax),
     X =< Xmax,
     Y =< Ymax,
     X >= 0,
     Y >= 0.

is_covid_free([X, Y]) :-
     get_adjacent([X, Y], Infected),
     forall(member(Cell, Infected), not(covid(Cell))).


% if we are in final cell, then all done
find_way(Xfinish, Yfinish, [[Xfinish, Yfinish]], _, (Xfinish, Yfinish)) :- !.

% backtracking search
find_way(Xfinish, Yfinish, 
    [[Xcurrent, Ycurrent] | PathTaken], Visited, (Xcurrent, Ycurrent)) :-
    member([Xcurrent, Ycurrent], Visited),
    get_adjacent((Xcurrent, Ycurrent), Adjacent_cells),
    exclude(covid, Adjacent_cells, Possible_cells),
    include(inside_map, Possible_cells, Possible_cells1),
    include(find_way(Xfinish, Yfinish,
          PathTaken, [[Xcurrent, Ycurrent] | Visited]), Possible_cells1, Found_way_cells),
    not(length(Found_way_cells, 0)).

% Simplification of function call
find_way_to(X, Y, Path) :-
     find_way(X, Y, Path, Visited, (0, 0)).


find_way2([[Xcurrent, Ycurrent] | Tail], [[[Xcurrent, Ycurrent] | CurrentPath] | PathTail], 
     Xfinish, Yfinish, Visited) :-
     (
     print([Xcurrent]),
     get_adjacent([Xcurrent, Ycurrent], Adjacent_cells),
     exclude(covid, Adjacent_cells, Possible_cells1),
     include(inside_map, Possible_cells1, Possible_cells2),
     exclude(member1(Visited), Possible_cells2, Possible_cells3),
     append(CurrentPath, [Xcurrent, Ycurrent], NewPath),
     append(Possible_cells3, Tail, Newlist),
     length(Possible_cells3, Length),
     copyN(NewPath, Length, NPaths),
     append(Npaths, PathTail, NewListPaths),
     find_way2(NewList, NewListPaths, Xfinish, Yfinish, [[Xcurrent, Ycurrent] | Visited])
     ).

find_way2([[Xfinish, Yfinish] | _], [[Xfinish, Yfinish] | []], Xfinish, Yfinish, _) :-
     true.

find_way_to2(X, Y, Path) :-
     find_way2([[0, 0]], [Path], X, Y, []).


find_way3([[Xcurrent, Ycurrent] | Tail], Xfinish, Yfinish, Visited) :-
     (
     print([Xcurrent, Ycurrent]),
     get_adjacent([Xcurrent, Ycurrent], Adjacent_cells),
     exclude(covid, Adjacent_cells, Possible_cells1),
     include(inside_map, Possible_cells1, Possible_cells2),
     exclude(member1(Visited), Possible_cells2, Possible_cells3),
     append(Possible_cells3, Tail, Newlist),
     find_way3(NewList, Xfinish, Yfinish, [[Xcurrent, Ycurrent] | Visited])
     ).

find_way3([[Xfinish, Yfinish] | _], Xfinish, Yfinish, _) :-
     print(Xfinish, Yfinish), true.


find_way_to3(X, Y) :-
     find_way3([[0, 0]], X, Y, []).

find_way4([Xcurrent, Ycurrent], [Xfinish, Yfinish], Visited, [[Xcurrent, Ycurrent] | Path]) :-
     get_adjacent([Xcurrent, Ycurrent], Adjacent_cells),
     member([Xnext, Ynext], Adjacent_cells),
     is_covid_free([Xnext, Ynext]),
     inside_map([Xnext, Ynext]),
     not(member([Xnext, Ynext], Visited)),
     find_way4([Xnext, Ynext], [Xfinish, Yfinish], [[Xcurent, Ycurrent] | Visited], Path),
     !.

find_way4([Xfinish, Yfinish], [Xfinish, Yfinish], _, [[Xfinish, Yfinish]]) :-
     true, !.


find_way_to4(X, Y, Path) :-
     find_way4([0, 0], [X, Y], [], Path).