:- dynamic limitcount/1.

maxbacktrack :-
    retract(limitcount(Current)),
    Current > 0,
    Next is Current-1,
    assert(limitcount(Next)).

firstN(N,Pred) :-
    retractall(limitcount(_)),
    asserta(limitcount(N)),
    Pred,
    maxbacktrack.

findN(N,Term,Pred,List) :-
    findall(Term,firstN(N,Pred),List).

:- dynamic min_path_length/1.

/**
 * Mikhail Rudakov
 * BS19-02
 * Prolog Programming Assignment
 * 
 **/

/**
 * map for the agent
 **/
:- dynamic map_xlimit/1.
:- dynamic map_ylimit/1.
:- dynamic start/1.
:- dynamic finish/1.
:- dynamic covid/1.
:- dynamic doctor/1.
:- dynamic mask/1.

infinity(1000).


/**
 * General purpose functions
 **/

convert_time(ExecutionTime, Minutes, Seconds, MilliSeconds) :-
     Minutes is ExecutionTime // 60000, Seconds is ExecutionTime // 1000 mod 60, MilliSeconds is ExecutionTime mod 1000.

output_results([], ExecutionTime) :-
     convert_time(ExecutionTime, M, S, Ms),
     write('Result: lose'), nl,
     write('Execution time: '), write(M), write(' min. '), 
     write(S), write(' sec. '), write(Ms), write(' ms.'),  
     nl, !.

output_results(Path, ExecutionTime) :-
     convert_time(ExecutionTime, M, S, Ms),
     write('Path on the map:'), nl,
     draw_map(Path),
     write('Result: win'), nl, 
     write('Shortest path is: '), nl, write(Path), nl,
     length(Path, Length),
     write('Length: '), write(Length), nl,
     write('Execution time: '), write(M), write(' min. '), 
     write(S), write(' sec. '), write(Ms), write(' ms.').

compute_distances([], []) :- !.
compute_distances([H | T], [[D, H] | T1]) :-
     chessboard_distance(H, D), compute_distances(T, T1).

get_adjacent((X, Y), L) :-
     Yup is Y + 1,
     Ydown is Y - 1,
     Xup is X + 1,
     Xdown is X - 1,
     L = [(Xup, Yup), (Xup, Y), (X, Yup), (Xdown, Yup), (Xup, Ydown), (Xdown, Y), (X, Ydown),
          (Xdown, Ydown)].

is_adjacent(P1, P2) :-
     get_adjacent(P1, L),
     member(P2, L).


is_adjacent_or_equal(P1, P2) :-
     P1 = P2;
     get_adjacent(P1, L),
     member(P2, L).

get_second([], []) :- !.

get_second([[_, Hb] | T], [Hb | T1]) :-
     get_second(T, T1).

get_adjacent_heuristic((X, Y), Lheuristic) :-
     Yup is Y + 1,
     Ydown is Y - 1,
     Xup is X + 1,
     Xdown is X - 1,
     L0 = [(Xup, Yup), (Xup, Y), (X, Yup), (Xdown, Yup), (Xup, Ydown), (Xdown, Y), (X, Ydown),
          (Xdown, Ydown)],
     compute_distances(L0, L1),
     sort(L1, L2),
     get_second(L2, Lheuristic).

is_adjacent_heuristic(P1, P2) :-
     get_adjacent_heuristic(P1, L),
     member(P2, L).

inside_map((X, Y)) :-
     X >= 0,
     Y >= 0,
     map_xlimit(Xmax),
     X < Xmax,
     map_ylimit(Ymax),
     Y < Ymax.
 

is_covid_free(Position) :-
     forall(is_adjacent(Position, Cell), not(covid(Cell))).

is_covid_free(1, _) :- !.

is_covid_free(0, Position) :-
     is_covid_free(Position).

is_doctor_or_mask(Position) :-
     doctor(Position); mask(Position).

is_doctor_or_mask(Position, 0) :-
     not(doctor(Position)), not(mask(Position)).

is_doctor_or_mask(Position, 1) :-
     doctor(Position); mask(Position).




/**
 * Backtracking Functionality
 **/

% <Helper functions>
maximum_possible_steps(D) :-
     map_xlimit(Xmax),
     map_ylimit(Ymax),
     D is 3 * max(Xmax, Ymax).
     %D is 12.

heuristics_shortest_path_set(Path) :-
     min_path_length(MinLength),
     length(Path, Length),
     Length < MinLength,
     assertz(min_path_length(Length)),
     retract(min_path_length(MinLength)), !.

heuristics_shortest_path_check(Path) :-
     min_path_length(MinLength),
     length(Path, Length),
     Length < MinLength, !.

% </Helper functions>

generate_path(CurrentPosition, []) :- 
     finish(CurrentPosition), !.

generate_path((Xcurrent, Ycurrent), [(Xnext, Ynext) | Path]) :- 
     finish((Xfinish, Yfinish)),
     (
     Xcurrent < Xfinish, Xnext is Xcurrent + 1;
     Xcurrent = Xfinish, Xnext = Xfinish;
     Xcurrent > Xfinish, Xnext is Xcurrent - 1
     ), 
     (
     Ycurrent < Yfinish, Ynext is Ycurrent + 1;
     Ycurrent = Yfinish, Ynext = Yfinish;
     Ycurrent > Yfinish, Ynext is Ycurrent - 1
     ), 
     generate_path((Xnext, Ynext), Path).

dfs(CurrentPosition, Visited, [CurrentPosition]) :-
     finish(CurrentPosition),
     heuristics_shortest_path_set(Visited),
     !.

dfs(CurrentPosition, Visited, [CurrentPosition | Path]) :-
     is_doctor_or_mask(CurrentPosition),
     generate_path(CurrentPosition, Path), 
     append(Visited, Path, PathHome),
     heuristics_shortest_path_set(PathHome),
     !.

dfs(CurrentPosition, Visited, [CurrentPosition | Path]) :-
     heuristics_shortest_path_check(Visited),
     is_adjacent_heuristic(CurrentPosition, NextPosition),
     inside_map(NextPosition),
     not(member(NextPosition, Visited)),
     is_covid_free(NextPosition),
     dfs(NextPosition, [NextPosition | Visited], Path).

find_path_dfs(Path) :-
     start(Start),
     finish(Home), inside_map(Home),
     dfs(Start, [Start], Path).

min_path_dfs(MinPath) :-    
     bagof(Path, find_path_dfs(Path), L),
     last(L, MinPath).

dfs() :-
     write('Backtracking: '), nl,
     retractall(min_path_length(_)),
     maximum_possible_steps(X),
     assertz(min_path_length(X)),
     statistics(walltime, _),
     (min_path_dfs(Path), !; Path = []),
     statistics(walltime, [_ | [ExecutionTime]]),
     output_results(Path, ExecutionTime),
     retractall(min_path_length(_)).

/**
 * A* Functionality
 **/

determine_safety(1, _, 1) :- !.
determine_safety(_, 1, 1) :- !.
determine_safety(0, 0, 0) :- !.

chessboard_distance((X, Y), Distance) :-
     finish((Xfinish, Yfinish)),
     Distance is max(abs(X - Xfinish), abs(Y - Yfinish)).

update_value(0, Value, [_ | T], [H1 | T1]) :-
     T1 = T, H1 = Value, !.

update_value(I, Value, [H | T], [H | T1]) :-
     I1 is I - 1,
     update_value(I1, Value, T, T1).

update_value(I, J, Value, Array, NewArray) :-
     nth0(I, Array, Row),
     update_value(J, Value, Row, NewRow),
     update_value(I, NewRow, Array, NewArray).

update_value(I, J, K, Value, Array, NewArray) :-
     nth0(I, Array, Row),
     update_value(J, K, Value, Row, NewRow),
     update_value(I, NewRow, Array, NewArray).

update_values([], _, Array, Array) :- !.
update_values([(I, J) | T], Value, Array, NewArray) :-
     update_values(T, Value, Array, TailArray),
     update_value(I, J, Value, TailArray, NewArray).

create_array(0, _, []) :- !.
create_array(I, Value, [Value | T]) :-
     I1 is I - 1,
     create_array(I1, Value, T).


create_array(0, _, _, []) :- !.
create_array(I, J, Value, [Array | T]) :-
     I1 is I - 1,
     create_array(J, Value, Array),
     create_array(I1, J, Value, T).

create_array(0, _, _, _, []) :- !.
create_array(I, J, K, Value, [Array | T]) :-
     I1 is I - 1,
     create_array(J, K, Value, Array),
     create_array(I1, J, K, Value, T).

nth0_2d(J, K, Array, Value) :-
     nth0(J, Array, Row),
     nth0(K, Row, Value).

nth0_3d(I, J, K, Array, Value) :-
     nth0(I, Array, Row),
     nth0_2d(J, K, Row, Value).

process_neighbours([], _, _, _, Distances, Predecessors, VerticesHeap, Distances, Predecessors, VerticesHeap) :- !.
process_neighbours([(Xcurrent, Ycurrent) | T], Safety, DistanceToCurrent, (Xpred, Ypred, Spred),
                    Distances, Predecessors, VerticesHeap, DistancesNew, PredecessorsNew, VerticesHeapNew) :- 
     (
          nth0_3d(Xcurrent, Ycurrent, Safety, Distances, DistanceOld),
          DistanceToCurrent < DistanceOld,
          (
          chessboard_distance((Xcurrent, Ycurrent), Heuristics),
          Priority is DistanceToCurrent + Heuristics,
          update_value(Xcurrent, Ycurrent, Safety, DistanceToCurrent, Distances, DistancesUpd),
          update_value(Xcurrent, Ycurrent, Safety, (Xpred, Ypred, Spred), Predecessors, PredecessorsUpd),
          add_to_heap(VerticesHeap, Priority, (Xcurrent, Ycurrent, Safety), VerticesHeapUpd) 
          ),
          !
          ;
          DistancesUpd = Distances, PredecessorsUpd = Predecessors, VerticesHeapUpd = VerticesHeap
     ),
     process_neighbours(T, Safety, DistanceToCurrent, (Xpred, Ypred, Spred),
               DistancesUpd, PredecessorsUpd, VerticesHeapUpd, DistancesNew, PredecessorsNew, VerticesHeapNew).



build_path((X, Y, _), _, [(X, Y)]) :- 
     start((X, Y)),
     !.

build_path((X, Y, S), Predecessors, PathHome) :- 
     nth0_3d(X, Y, S, Predecessors, (Xpred, Ypred, Spred)),
     build_path((Xpred, Ypred, Spred), Predecessors, PathHomeHead),
     append(PathHomeHead, [(X, Y)], PathHome).

astar(VerticesHeap, _, Predecessors, PathHome) :-
     finish((Xfinish, Yfinish)),
     min_of_heap(VerticesHeap, _, (Xfinish, Yfinish, SafetyFinish)),
     build_path((Xfinish, Yfinish, SafetyFinish), Predecessors, PathHome),
     !.

astar(VerticesHeap, Distances, Predecessors, PathHome) :-
     get_from_heap(VerticesHeap, _Priority, (Xcurrent, Ycurrent, SafetyCurrent), VerticesHeap1),
     Position = (Xcurrent, Ycurrent),
     nth0_3d(Xcurrent, Ycurrent, SafetyCurrent, Distances, DistanceToCurrent),
     DistanceToNext is DistanceToCurrent + 1,
     get_adjacent(Position, AdjacentCells),
     include(inside_map, AdjacentCells, AdjacentCells1),
     is_doctor_or_mask(Position, SafetyPosition),
     determine_safety(SafetyCurrent, SafetyPosition, SafetyNext),
     include(is_covid_free(SafetyNext), AdjacentCells1, AdjacentCellsToProcess),
     process_neighbours(AdjacentCellsToProcess, SafetyNext, DistanceToNext, (Xcurrent, Ycurrent, SafetyCurrent),
               Distances, Predecessors, VerticesHeap1, DistancesNew, PredecessorsNew, VerticesHeapNew),
     astar(VerticesHeapNew, DistancesNew, PredecessorsNew, PathHome).
     

min_path_astar(PathHome) :-
     empty_heap(Heap),
     start((Xstart, Ystart)),
     add_to_heap(Heap, 0, (Xstart, Ystart, 0), VerticesHeap),
     map_xlimit(Xmax), map_ylimit(Ymax), infinity(Inf),
     create_array(Xmax, Ymax, 2, Inf, Distances),
     update_value(Xstart, Ystart, 0, 0, Distances, Distances1),
     create_array(Xmax, Ymax, 2, -1, Predecessors),
     bagof(Path, astar(VerticesHeap, Distances1, Predecessors, Path), L),
     nth0(0, L, PathHome).

astar() :-
     write('A* algorithm: '), nl,
     statistics(walltime, _),
     (min_path_astar(Path), !; Path = []),
     statistics(walltime, [_ | [ExecutionTime]]),
     output_results(Path, ExecutionTime).

% Testing function and maps test cases


reset_environment() :- 
     retractall(map_xlimit(_)),
     retractall(map_ylimit(_)),
     retractall(start(_)),
     retractall(finish(_)),
     retractall(covid(_)),
     retractall(doctor(_)),
     retractall(mask(_)).

get_environment(Xlimit, Ylimit, Actor, Home, Doctor, Mask, Covid) :-
     map_xlimit(Xlimit), map_ylimit(Ylimit), 
     bagof(Xa, start(Xa), Actor), 
     bagof(Xb, finish(Xb), Home), 
     bagof(Xc, doctor(Xc), Doctor), 
     bagof(Xd, mask(Xd), Mask), 
     bagof(Xe, covid(Xe), Covid),
     !. 

set_environment(t1) :-
     assertz(map_xlimit(9)),
     assertz(map_ylimit(9)),
     assertz(start((0, 0))),
     assertz(finish((1, 7))),
     assertz(covid((1, 4))),
     assertz(covid((6, 7))),
     assertz(doctor((4, 4))),
     assertz(mask((7, 1))), !.


set_environment(t2) :-
     assertz(map_xlimit(9)),
     assertz(map_ylimit(9)),
     assertz(start((0, 0))),
     assertz(finish((8, 8))),
     assertz(covid((3, 7))),
     assertz(covid((2, 3))),
     assertz(doctor((5, 1))),
     assertz(mask((1, 5))), !.

set_environment(t3) :-
     assertz(map_xlimit(9)),
     assertz(map_ylimit(9)),
     assertz(start((0, 0))),
     assertz(finish((8, 0))),
     assertz(covid((6, 3))),
     assertz(covid((6, 1))),
     assertz(doctor((6, 8))),
     assertz(mask((1, 4))), !.

set_environment(t4) :-
     assertz(map_xlimit(9)),
     assertz(map_ylimit(9)),
     assertz(start((0, 0))),
     assertz(finish((0, 5))),
     assertz(covid((1, 3))),
     assertz(covid((4, 0))),
     assertz(doctor((3, 3))),
     assertz(mask((6, 1))), !.

set_environment(tsmall) :-
     assertz(map_xlimit(4)),
     assertz(map_ylimit(4)),
     assertz(start((0, 0))),
     assertz(finish((1, 3))),
     assertz(covid((3, 0))),
     assertz(covid((3, 3))),
     assertz(doctor((0, 2))),
     assertz(mask((0, 1))), !.

set_environment(tsmall_lose) :-
     assertz(map_xlimit(4)),
     assertz(map_ylimit(4)),
     assertz(start((0, 0))),
     assertz(finish((3, 1))),
     assertz(covid((3, 0))),
     assertz(covid((3, 3))),
     assertz(doctor((0, 21))),
     assertz(mask((0, 11))), !.

set_environment(tlose) :-
     assertz(map_xlimit(9)),
     assertz(map_ylimit(9)),
     assertz(start((0, 0))),
     assertz(finish((8, 0))),
     assertz(covid((2, 1))),
     assertz(covid((0, 3))),
     assertz(doctor((0, 4))),
     assertz(mask((3, 7))), !.

set_environment(ttest) :-
     assertz(map_xlimit(6)),
     assertz(map_ylimit(6)),
     assertz(start((0, 0))),
     assertz(finish((4, 0))),
     assertz(covid((2, 0))),
     assertz(covid((2, 3))),
     assertz(doctor((1, 2))),
     assertz(mask((2, 1))), !.

set_environment(N1) :-
     (
     N1 < 10, N is N1 - 1,
     Xstart = 0, Ystart = 0,
     random_between(0, N, Xfinish), random_between(0, N, Yfinish),
     random_between(0, N, Xmask), random_between(0, N, Ymask),
     random_between(0, N, Xdoc), random_between(0, N, Ydoc),
     random_between(0, N, Xcov1), random_between(0, N, Ycov1),
     random_between(0, N, Xcov2), random_between(0, N, Ycov2),
     unique([(Xfinish, Yfinish), (Xmask, Ymask), (Xdoc, Ydoc), (Xcov1, Ycov1), (Xcov2, Ycov2), 
               (Xstart, Ystart) ]),

     not(is_adjacent_or_equal((Xstart, Ystart), (Xcov1, Ycov1))),
     not(is_adjacent_or_equal((Xstart, Ystart), (Xcov2, Ycov2))),

     not(is_adjacent_or_equal((Xmask, Ymask), (Xcov1, Ycov1))),
     not(is_adjacent_or_equal((Xmask, Ymask), (Xcov2, Ycov2))),
   
     not(is_adjacent_or_equal((Xdoc, Ydoc), (Xcov1, Ycov1))),
     not(is_adjacent_or_equal((Xdoc, Ydoc), (Xcov2, Ycov2))),

     assertz(map_xlimit(N1)),
     assertz(map_ylimit(N1)),
     assertz(start((Xstart, Ystart))),
     assertz(finish((Xfinish, Yfinish))),
     assertz(covid((Xcov1, Ycov1))),
     assertz(covid((Xcov2, Ycov2))),
     assertz(doctor((Xdoc, Ydoc))),
     assertz(mask((Xmask, Ymask))), !
     ) ;
     set_environment(N1).

unique([]) :- !.
unique([H | T]) :-
     not(member(H, T)), 
     unique(T).


print_array([]) :- !.
print_array([H | T]) :-
     write(H), print_array(T).

print_array_2d([]) :- !.
print_array_2d([H | T]) :-
     write("|"),
     print_array(H), 
     write("|"), nl,
     print_array_2d(T).

print_array_2d(Xlimit, Map) :-
     create_array(Xlimit, "=", Line),
     append([" "], Line, Line1),
     append(Line1, [" "], LineUp),
     create_array(Xlimit, "=", Line2),
     append([" "], Line2, Line3),
     append(Line3, [" "], LineDown),
     print_array(LineUp), nl,
     print_array_2d(Map),
     print_array(LineDown), nl.

change_coordinates(_, _, [], []) :- !.
change_coordinates(Xlimit, Ylimit, [(X, Y) | T], [(X1, Y1) | T1]) :-
     Y1 is X, X1 is Ylimit - Y - 1,
     change_coordinates(Xlimit, Ylimit, T, T1).

draw_map(Path) :-
     get_environment(Xlimit, Ylimit, Actor0, Home0, Doctor0, Mask0, Covid0),
     create_array(Xlimit, Ylimit, '.', Map),
     change_coordinates(Xlimit, Ylimit, Actor0, Actor1),
     change_coordinates(Xlimit, Ylimit, Home0, Home1),
     change_coordinates(Xlimit, Ylimit, Doctor0, Doctor1),
     change_coordinates(Xlimit, Ylimit, Mask0, Mask1),
     change_coordinates(Xlimit, Ylimit, Covid0, Covid1),
     change_coordinates(Xlimit, Ylimit, Path, Path1),
     include(inside_map, Actor1, Actor2), update_values(Actor2, 'A', Map, Map1),
     include(inside_map, Home1, Home2), update_values(Home2, 'H', Map1, Map2),
     include(inside_map, Doctor1, Doctor2), update_values(Doctor2, 'D', Map2, Map3),
     include(inside_map, Mask1, Mask2), update_values(Mask2, 'M', Map3, Map4),
     include(inside_map, Covid1, Covid2), update_values(Covid2, 'C', Map4, Map5),
     update_values(Path1, 'o', Map5, MapFinal),
     print_array_2d(Xlimit, MapFinal).

draw_map() :-
     get_environment(Xlimit, Ylimit, Actor0, Home0, Doctor0, Mask0, Covid0),
     create_array(Xlimit, Ylimit, '.', Map),
     change_coordinates(Xlimit, Ylimit, Actor0, Actor1),
     change_coordinates(Xlimit, Ylimit, Home0, Home1),
     change_coordinates(Xlimit, Ylimit, Doctor0, Doctor1),
     change_coordinates(Xlimit, Ylimit, Mask0, Mask1),
     change_coordinates(Xlimit, Ylimit, Covid0, Covid1),
     include(inside_map, Actor1, Actor2), update_values(Actor2, 'A', Map, Map1),
     include(inside_map, Home1, Home2), update_values(Home2, 'H', Map1, Map2),
     include(inside_map, Doctor1, Doctor2), update_values(Doctor2, 'D', Map2, Map3),
     include(inside_map, Mask1, Mask2), update_values(Mask2, 'M', Map3, Map4),
     include(inside_map, Covid1, Covid2), update_values(Covid2, 'C', Map4, MapFinal),
     print_array_2d(Xlimit, MapFinal).


test(X) :-
     write('========= START OF TEST CASE ========='),
     write('Running on map '), write(X), write(':'), nl, nl,
     reset_environment(), set_environment(X),
     get_environment(Xlimit, Ylimit, Actor, Home, Doctor, Mask, Covid),
     write('Map is '), write(Xlimit), write(' x '), write(Ylimit), nl,
     write('Actor: '), write(Actor), nl,
     write('Home: '), write(Home), nl,
     write('Covid: '), write(Covid), nl,
     write('Mask: '), write(Mask), nl,
     write('Doctor: '), write(Doctor), nl, nl,
     draw_map(),
     nl,
     dfs(), nl,
     astar(), nl,
     write('========= END OF TEST CASE ========='),
     reset_environment().