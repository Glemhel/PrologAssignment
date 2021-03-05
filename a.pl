/**
 * Mikhail Rudakov
 * BS19-02
 * Prolog Programming Assignment
 * 
 **/

/**
 * map for the agent
 **/
map_xlimit(3).
map_ylimit(3).
covid(1, 0).
covid(0, 1).
finish(2, 0).
man(a).

% if we are in final cell, then all done
find_way(Xfinish, Yfinish, Xfinish, Yfinish, _, _, [], _) :- !.

% backtracking search
find_way(Xcurrent, Ycurrent, Xfinish, Yfinish, Xmax, Ymax, 
    [[Xcurrent, Ycurrent] | PathTaken], Visited) :-
    Yup is Ycurrent + 1,
    Ydown is Ycurrent - 1,
    Xup is Xcurrent + 1,
    Xdown is Xcurrent - 1,
    Xcurrent =< Xmax,
    Ycurrent =< Ymax,
    Xcurrent >= 0,
    Ycurrent >= 0,

    not(covid(Xcurrent, Ycurrent)),
    not(member([Xcurrent, Ycurrent], Visited)),
    (
     find_way(Xup, Yup, Xfinish, Yfinish, Xmax, Ymax,
          PathTaken, [[Xcurrent, Ycurrent] | Visited]);
     find_way(Xcurrent, Yup, Xfinish, Yfinish, Xmax, Ymax, 
          PathTaken, [[Xcurrent, Ycurrent] | Visited]);
     find_way(Xup, Ycurrent, Xfinish, Yfinish, Xmax, Ymax,
          PathTaken, [[Xcurrent, Ycurrent] | Visited]);
     find_way(Xcurrent, Ydown, Xfinish, Yfinish, Xmax, Ymax,
          PathTaken, [[Xcurrent, Ycurrent] | Visited]);
     find_way(Xup, Ydown, Xfinish, Yfinish, Xmax, Ymax,
          PathTaken, [[Xcurrent, Ycurrent] | Visited]);
     find_way(Xdown, Ycurrent, Xfinish, Yfinish, Xmax, Ymax,
          PathTaken, [[Xcurrent, Ycurrent] | Visited]);
     find_way(Xdown, Yup, Xfinish, Yfinish, Xmax, Ymax,
          PathTaken, [[Xcurrent, Ycurrent] | Visited]);
     find_way(Xdown, Ydown, Xfinish, Yfinish, Xmax, Ymax,
          PathTaken, [[Xcurrent, Ycurrent] | Visited])
     ).

% Simplification of function call
find_way_to(X, Y, Path) :-
     map_xlimit(Xmax),
     map_ylimit(Ymax),
     find_way(0, 0, X, Y, Xmax, Ymax, Path, []).