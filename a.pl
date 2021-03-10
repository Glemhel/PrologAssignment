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
map_xlimit(7).
map_ylimit(7).
start((0, 0)).
finish((6, 0)).
covid((6, 2)).
covid((5, 2)).
doctor((66, 8)).
mask((16, 4)).

infinity(1000).


/**
 * General purpose functions
 **/

output_results(Path, Length, ExecutionTime) :-
     write('Minimum path is: '), nl, write(Path), nl,
     write('Length: '), write(Length), nl,
     write('Execution time: '), write(ExecutionTime), write(' ms.'), nl.

get_adjacent((X, Y), L) :-
     Yup is Y + 1,
     Ydown is Y - 1,
     Xup is X + 1,
     Xdown is X - 1,
     L = [(Xup, Yup), (Xup, Y), (X, Yup), (Xdown, Yup), (Xup, Ydown), (Xdown, Y), (X, Ydown),
          (Xdown, Ydown)].

is_adjacent((X1, Y1), (X2, Y2)) :-
     get_adjacent((X1, Y1), L),
     member((X2, Y2), L).

inside_map((X, Y)) :-
     map_xlimit(Xmax),
     map_ylimit(Ymax),
     X < Xmax,
     Y < Ymax,
     X >= 0,
     Y >= 0.


is_covid_free(Position) :-
     forall(is_adjacent(Position, Cell), not(covid(Cell))).

is_covid_free(1, _) :- !.

is_covid_free(0, Position) :-
     is_covid_free(Position).

is_doctor_or_mask(Position) :-
     doctor(Position); mask(Position).


/**
 * Backtracking Functionality
 **/

% <Helper functions>
maximum_possible_steps(D) :-
     map_xlimit(Xmax),
     map_ylimit(Ymax),
     D is 3 * max(Xmax, Ymax).

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
     is_covid_free(CurrentPosition),
     is_adjacent(CurrentPosition, NextPosition),
     inside_map(NextPosition),
     not(member(NextPosition, Visited)),
     dfs(NextPosition, [CurrentPosition | Visited], Path).

find_path_dfs(Path) :-
     start(Start),
     dfs(Start, [Start], Path).

min_path_dfs(MinPath) :-    
     bagof(Path, find_path_dfs(Path), L),
     last(L, MinPath).

dfs() :-
     retractall(min_path_length(_)),
     maximum_possible_steps(X),
     assertz(min_path_length(X)),
     statistics(walltime, _),
     min_path_dfs(Path),
     statistics(walltime, [_ | [ExecutionTime]]),
     length(Path, Length),
     output_results(Path, Length, ExecutionTime),
     retractall(min_path_length(_)).

/**
 * A* Functionality
 **/

heuristic_distance((X, Y), Distance) :-
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
          heuristic_distance((Xcurrent, Ycurrent), Heuristics),
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
     statistics(walltime, _),
     min_path_astar(Path),
     statistics(walltime, [_ | [ExecutionTime]]),
     length(Path, Length),
     output_results(Path, Length, ExecutionTime).

% Testing function

