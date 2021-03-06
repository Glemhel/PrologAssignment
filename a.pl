/**
 * Mikhail Rudakov
 * BS19-02
 * Prolog Programming Assignment
 * 
 **/

/**
 * map for the agent
 **/
map_xlimit(2).
map_ylimit(2).
covid([40, 40]).
%covid([3, 8]).
doctor([10, 20]).
mask([10, 10]).

get_adjacent([X, Y], L) :-
     Yup is Y + 1,
     Ydown is Y - 1,
     Xup is X + 1,
     Xdown is X - 1,
     L = [[Xup, Y], [Xdown, Y], [X, Yup], [X, Ydown], [Xup, Yup],
      [Xup, Ydown], [Xdown, Yup], [Xdown, Ydown]].

is_adjacent([X1, Y1], [X2, Y2]) :-
     X2 is X1 - 1, Y2 is Y1 - 1;
     X2 is X1 - 1, Y2 is Y1;
     X2 is X1 - 1, Y2 is Y1 + 1;
     X2 is X1, Y2 is Y1 - 1;
     X2 is X1, Y2 is Y1 + 1;
     X2 is X1 + 1, Y2 is Y1 - 1;
     X2 is X1 + 1, Y2 is Y1;
     X2 is X1 + 1, Y2 is Y1 + 1.


inside_map([X, Y]) :-
     map_xlimit(Xmax),
     map_ylimit(Ymax),
     X =< Xmax,
     Y =< Ymax,
     X >= 0,
     Y >= 0.

is_covid_free([X, Y]) :-
     forall(is_adjacent([X, Y], Cell), not(covid(Cell))).

is_doctor_or_mask([X, Y], 1) :-
     doctor([X, Y]); mask([X, Y]).

is_doctor_or_mask([X, Y], 0) :-
     not(doctor([X, Y])), not(mask([X, Y])).

find_way([Xcurrent, Ycurrent, Safety], [Xfinish, Yfinish], Visited, [[Xcurrent, Ycurrent] | Path]) :-
     (is_covid_free([Xcurrent, Ycurrent]); Safety = 1),
     get_adjacent([Xcurrent, Ycurrent], Adjacent_cells),
     member([Xnext, Ynext], Adjacent_cells),
     inside_map([Xnext, Ynext]),
     not(member([Xnext, Ynext, Safety], Visited)),
     is_doctor_or_mask([Xcurrent, Ycurrent], Safety1),
     (
          ((Safety = 1; Safety1 = 1),
          find_way([Xnext, Ynext, 1], [Xfinish, Yfinish], 
               [[Xcurrent, Ycurrent, 1] | Visited], Path), !
          );
          ((Safety = 0, Safety1 = 0),
          find_way([Xnext, Ynext, 0], [Xfinish, Yfinish], 
               [[Xcurrent, Ycurrent, 0] | Visited], Path), !
          )
     ).


find_way([Xfinish, Yfinish, _], [Xfinish, Yfinish], _, [[Xfinish, Yfinish]]) :-
     true, !.


find_way_to(X, Y, Path) :-
     find_way([0, 0, 0], [X, Y], [], Path).


find_way_to4(X, Y, Path) :-
     find_way4([0, 0], [X, Y], [], Path). 

find_way4([Xfinish, Yfinish], [Xfinish, Yfinish], _, [[Xfinish, Yfinish]]) :-
     true, !.

find_way4([Xcurrent, Ycurrent], [Xfinish, Yfinish], Visited, [[Xcurrent, Ycurrent] | Path]) :-
     is_adjacent([Xcurrent, Ycurrent], [Xnext, Ynext]),
     is_covid_free([Xnext, Ynext]),
     inside_map([Xnext, Ynext]),
     not(member([Xnext, Ynext], Visited)),
     find_way4([Xnext, Ynext], [Xfinish, Yfinish], [[Xcurrent, Ycurrent] | Visited], Path), !.



