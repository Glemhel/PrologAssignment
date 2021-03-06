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
covid((1, 0)).
covid((0, 1)).


get_adjacent((X, Y), L) :-
     Yup is Y + 1,
     Ydown is Y - 1,
     Xup is X + 1,
     Xdown is X - 1,
     L = [(Xup, Yup), (Xup, Y), (X, Yup), (Xdown, Yup), (Xup, Ydown), (Xdown, Y), (X, Ydown),
          (Xdown, Ydown)].


inside_map((X, Y)) :-
     map_xlimit(Xmax),
     map_ylimit(Ymax),
     X =< Xmax,
     Y =< Ymax,
     X >= 0,
     Y >= 0.


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