% PROLOG PROGRAMMING ASSIGNMENT
% by MIKHAIL RUDAKOV
% BS19-02

% TESTING FUNCTION
% if X is a number, generate random X by X map. 3 < X < 10.
% if X is atom, try to find test map for it
% X could be 't1' - 't15', and some specific values
% 'tsmall1' - 'tsmall5', 'teasy', 'tsmall_lose'.
% if X is 'ttest' - you can specify your values here. See function declaration for test(ttest).
%
test(X) :-
     write('========= START OF TEST CASE ========='), nl,
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
     astar(), nl, nl,
     dfs('Variant1'), nl, nl,
     dfs('Variant2'), nl,
     write('========= END OF TEST CASE ========='),
     reset_environment().

% input your desired parameters here
% at least one entity of each type required. If some kind of entity should be absent on the map,
% place it far way outside, like covid((100, 100)). Any numbers of covids/doctors/masks allowed.
set_environment(ttest) :-
     assertz(map_xlimit(6)),
     assertz(map_ylimit(6)),
     assertz(start((0, 0))),
     assertz(finish((1, 5))),
     assertz(covid((0, 4))),
     assertz(covid((2, 2))),
     assertz(doctor((2, 5))),
     assertz(mask((5, 5))), !.

% Test cases
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
     assertz(finish((0, 5))),
     assertz(covid((1, 3))),
     assertz(covid((4, 0))),
     assertz(doctor((3, 3))),
     assertz(mask((6, 1))), !.

set_environment(t4) :-
     assertz(map_xlimit(9)),
     assertz(map_ylimit(9)),
     assertz(start((0, 0))),
     assertz(finish((6, 2))),
     assertz(covid((2, 3))),
     assertz(covid((5, 4))),
     assertz(doctor((0, 5))),
     assertz(mask((2, 8))), !.

set_environment(t5) :-
     assertz(map_xlimit(9)),
     assertz(map_ylimit(9)),
     assertz(start((0, 0))),
     assertz(finish((7, 5))),
     assertz(covid((6, 3))),
     assertz(covid((5, 4))),
     assertz(doctor((5, 8))),
     assertz(mask((8, 3))), !.

set_environment(t6) :-
     assertz(map_xlimit(9)),
     assertz(map_ylimit(9)),
     assertz(start((0, 0))),
     assertz(finish((0, 8))),
     assertz(covid((4, 7))),
     assertz(covid((1, 6))),
     assertz(doctor((1, 0))),
     assertz(mask((0, 1))), !.

set_environment(t7) :-
     assertz(map_xlimit(9)),
     assertz(map_ylimit(9)),
     assertz(start((0, 0))),
     assertz(finish((7, 7))),
     assertz(covid((6, 3))),
     assertz(covid((1, 4))),
     assertz(doctor((6, 0))),
     assertz(mask((0, 6))), !.

set_environment(t8) :-
     assertz(map_xlimit(9)),
     assertz(map_ylimit(9)),
     assertz(start((0, 0))),
     assertz(finish((8, 2))),
     assertz(covid((4, 2))),
     assertz(covid((7, 0))),
     assertz(doctor((2, 6))),
     assertz(mask((1, 3))), !.

set_environment(t9) :-
     assertz(map_xlimit(9)),
     assertz(map_ylimit(9)),
     assertz(start((0, 0))),
     assertz(finish((3, 6))),
     assertz(covid((7, 8))),
     assertz(covid((2, 5))),
     assertz(doctor((8, 5))),
     assertz(mask((5, 5))), !.

set_environment(t10) :-
     assertz(map_xlimit(9)),
     assertz(map_ylimit(9)),
     assertz(start((0, 0))),
     assertz(finish((2, 1))),
     assertz(covid((2, 2))),
     assertz(covid((3, 6))),
     assertz(doctor((5, 5))),
     assertz(mask((4, 4))), !.

set_environment(t11) :-
     assertz(map_xlimit(9)),
     assertz(map_ylimit(9)),
     assertz(start((0, 0))),
     assertz(finish((8, 0))),
     assertz(covid((2, 1))),
     assertz(covid((0, 3))),
     assertz(doctor((0, 4))),
     assertz(mask((3, 7))), !.

set_environment(t12) :-
     assertz(map_xlimit(9)),
     assertz(map_ylimit(9)),
     assertz(start((0, 0))),
     assertz(finish((8, 1))),
     assertz(covid((6, 1))),
     assertz(covid((7, 3))),
     assertz(doctor((4, 5))),
     assertz(mask((8, 5))), !.

set_environment(t13) :-
     assertz(map_xlimit(9)),
     assertz(map_ylimit(9)),
     assertz(start((0, 0))),
     assertz(finish((8, 8))),
     assertz(covid((3, 4))),
     assertz(covid((0, 6))),
     assertz(doctor((6, 2))),
     assertz(mask((4, 6))), !.

set_environment(tlong1) :-
     assertz(map_xlimit(9)),
     assertz(map_ylimit(9)),
     assertz(start((0, 0))),
     assertz(finish((8, 0))),
     assertz(covid((6, 3))),
     assertz(covid((6, 1))),
     assertz(doctor((6, 8))),
     assertz(mask((1, 4))), !.

set_environment(tlong2) :-
     assertz(map_xlimit(9)),
     assertz(map_ylimit(9)),
     assertz(start((0, 0))),
     assertz(finish((4, 6))),
     assertz(covid((5, 5))),
     assertz(covid((8, 7))),
     assertz(doctor((5, 8))),
     assertz(mask((4, 0))), !.

set_environment(tsmall1) :-
     assertz(map_xlimit(4)),
     assertz(map_ylimit(4)),
     assertz(start((0, 0))),
     assertz(finish((1, 3))),
     assertz(covid((3, 0))),
     assertz(covid((3, 3))),
     assertz(doctor((0, 2))),
     assertz(mask((0, 1))), !.

set_environment(tsmall2) :-
     assertz(map_xlimit(6)),
     assertz(map_ylimit(6)),
     assertz(start((0, 0))),
     assertz(finish((5, 1))),
     assertz(covid((5, 1))),
     assertz(covid((3, 1))),
     assertz(doctor((3, 5))),
     assertz(mask((0, 3))), !.

set_environment(tsmall3) :-
     assertz(map_xlimit(5)),
     assertz(map_ylimit(5)),
     assertz(start((0, 0))),
     assertz(finish((3, 2))),
     assertz(covid((2, 2))),
     assertz(covid((0, 3))),
     assertz(doctor((4, 0))),
     assertz(mask((2, 4))), !.

set_environment(tsmall4) :-
     assertz(map_xlimit(5)),
     assertz(map_ylimit(5)),
     assertz(start((0, 0))),
     assertz(finish((4, 4))),
     assertz(covid((2, 0))),
     assertz(covid((4, 2))),
     assertz(doctor((1, 2))),
     assertz(mask((1, 4))), !.

set_environment(tsmall5) :-
     assertz(map_xlimit(5)),
     assertz(map_ylimit(5)),
     assertz(start((4, 0))),
     assertz(finish((2, 2))),
     assertz(covid((2, 1))),
     assertz(covid((2, 2))),
     assertz(doctor((11, 2))),
     assertz(mask((11, 4))), !.

set_environment(easy) :-
     assertz(map_xlimit(3)),
     assertz(map_ylimit(3)),
     assertz(start((0, 0))),
     assertz(finish((2, 2))),
     assertz(covid((2, 0))),
     assertz(covid((44, 2))),
     assertz(doctor((11, 2))),
     assertz(mask((11, 4))), !.

set_environment(tsmall_lose) :-
     assertz(map_xlimit(4)),
     assertz(map_ylimit(4)),
     assertz(start((0, 0))),
     assertz(finish((3, 1))),
     assertz(covid((3, 0))),
     assertz(covid((3, 3))),
     assertz(doctor((2, 3))),
     assertz(mask((2, 2))), !.

set_environment(tlose1) :-
     assertz(map_xlimit(9)),
     assertz(map_ylimit(9)),
     assertz(start((0, 0))),
     assertz(finish((8, 0))),
     assertz(covid((2, 1))),
     assertz(covid((0, 3))),
     assertz(doctor((0, 4))),
     assertz(mask((3, 7))), !.

set_environment(tlose2) :-
     assertz(map_xlimit(9)),
     assertz(map_ylimit(9)),
     assertz(start((0, 0))),
     assertz(finish((8, 8))),
     assertz(covid((5, 7))),
     assertz(covid((7, 5))),
     assertz(doctor((7, 8))),
     assertz(mask((8, 7))), !.


% Random generation of N1 * N1 map
set_environment(N1) :-
     N1 < 10, N1 > 3, 
     ((
     N is N1 - 1,
     Xstart = 0, Ystart = 0,
     random_between(0, N, Xfinish), random_between(0, N, Yfinish),
     random_between(0, N, Xmask), random_between(0, N, Ymask),
     random_between(0, N, Xdoc), random_between(0, N, Ydoc),
     random_between(0, N, Xcov1), random_between(0, N, Ycov1),
     random_between(0, N, Xcov2), random_between(0, N, Ycov2),
     % check conditions - not in the same cell, etc.
     unique([(Xfinish, Yfinish), (Xmask, Ymask), (Xdoc, Ydoc), (Xcov1, Ycov1), (Xcov2, Ycov2), 
               (Xstart, Ystart) ]),

     not(is_adjacent((Xstart, Ystart), (Xcov1, Ycov1))),
     not(is_adjacent((Xstart, Ystart), (Xcov2, Ycov2))),

     not(is_adjacent((Xmask, Ymask), (Xcov1, Ycov1))),
     not(is_adjacent((Xmask, Ymask), (Xcov2, Ycov2))),
   
     not(is_adjacent((Xdoc, Ydoc), (Xcov1, Ycov1))),
     not(is_adjacent((Xdoc, Ydoc), (Xcov2, Ycov2))),

     assertz(map_xlimit(N1)),
     assertz(map_ylimit(N1)),
     assertz(start((Xstart, Ystart))),
     assertz(finish((Xfinish, Yfinish))),
     assertz(covid((Xcov1, Ycov1))),
     assertz(covid((Xcov2, Ycov2))),
     assertz(doctor((Xdoc, Ydoc))),
     assertz(mask((Xmask, Ymask))), !
     ) ;
     set_environment(N1)).

% declaring variables and constants
:- dynamic map_xlimit/1.
:- dynamic map_ylimit/1.
:- dynamic start/1.
:- dynamic finish/1.
:- dynamic covid/1.
:- dynamic doctor/1.
:- dynamic mask/1.
:- dynamic min_path_length/1.
infinity(1000).

% ARRAY MANIPULATION FUNCTIONS
% ARRAY MANIPULATION FUNCTIONS
% ARRAY MANIPULATION FUNCTIONS


% Update value with index I to Value in array and return newarray as 4th argument
% If index is 0, update head and tails are the same
update_value(0, Value, [_ | T], [H1 | T1]) :-
     T1 = T, H1 = Value, !.
% If index is not 0, go 1 step further
update_value(I, Value, [H | T], [H | T1]) :-
     I1 is I - 1,
     update_value(I1, Value, T, T1).

% Update value in 2d array
update_value(I, J, Value, Array, NewArray) :-
     % find needed row
     nth0(I, Array, Row),
     % update value in it
     update_value(J, Value, Row, NewRow),
     % update value in array of rows as in 1d array
     update_value(I, NewRow, Array, NewArray).

% Update value in 3d array - same as 2d
update_value(I, J, K, Value, Array, NewArray) :-
     nth0(I, Array, Row),
     update_value(J, K, Value, Row, NewRow),
     update_value(I, NewRow, Array, NewArray).

% For each element (I, J) from input array, set Array[I][J] = Value, return results in NewArray
update_values([], _, Array, Array) :- !.
update_values([(I, J) | T], Value, Array, NewArray) :-
     update_values(T, Value, Array, TailArray),
     update_value(I, J, Value, TailArray, NewArray).

% Create 1d array of Value with I elements
create_array(0, _, []) :- !.
create_array(I, Value, [Value | T]) :-
     I1 is I - 1,
     create_array(I1, Value, T).

% Create 2d array
create_array(0, _, _, []) :- !.
create_array(I, J, Value, [Array | T]) :-
     I1 is I - 1,
     create_array(J, Value, Array),
     create_array(I1, J, Value, T).

% Create 3d array
create_array(0, _, _, _, []) :- !.
create_array(I, J, K, Value, [Array | T]) :-
     I1 is I - 1,
     create_array(J, K, Value, Array),
     create_array(I1, J, K, Value, T).

% Get value from 2d array
nth0_2d(J, K, Array, Value) :-
     nth0(J, Array, Row),
     nth0(K, Row, Value).

% Get value from 3d array
nth0_3d(I, J, K, Array, Value) :-
     nth0(I, Array, Row),
     nth0_2d(J, K, Row, Value).



% GENERAL PURPOSE FUNCTIONALITY
% GENERAL PURPOSE FUNCTIONALITY
% GENERAL PURPOSE FUNCTIONALITY

% Get list of adjacent cells for the cell
get_adjacent((X, Y), L) :-
     Yup is Y + 1,
     Ydown is Y - 1,
     Xup is X + 1,
     Xdown is X - 1,
     L = [(Xup, Yup), (Xup, Y), (X, Yup), (Xdown, Yup), (Xup, Ydown), (Xdown, Y), (X, Ydown),
          (Xdown, Ydown)].

% Check whether two cells are adjacent or not
is_adjacent(P1, P2) :-
     get_adjacent(P1, L),
     member(P2, L).

% Given array of pairs, return array of second elements of those pairs
get_second([], []) :- !.
get_second([[_, Hb] | T], [Hb | T1]) :-
     get_second(T, T1).

% Given array of tuples, return array of third elements of those pairs
get_third([], []) :- !.
get_third([[_, _, Hb] | T], [Hb | T1]) :-
     get_third(T, T1).

% For each Cell in the list, calculate number of covid-free cell in its neighbourhood and
% safe to the another list - used for Variant 2 of Backtracking
compute_covid_neighbours([], []) :- !.
compute_covid_neighbours([[D, Cell] | T], [[D, N, Cell] | T1]) :-
     get_adjacent(Cell, AdjacentCells),
     append([Cell], AdjacentCells, AdjacentOrEqual),
     exclude(is_covid_free, AdjacentOrEqual, InfectedCells),
     length(InfectedCells, N),
     compute_covid_neighbours(T, T1).

% Get adjacent cells for the cell, and return them in special order (heuristics for backtracking)
% Those closer to home are going first
% VARIANT 1 BACKTRACKING
get_adjacent_heuristic((X, Y), Lheuristic) :-
     Yup is Y + 1,
     Ydown is Y - 1,
     Xup is X + 1,
     Xdown is X - 1,
     L0 = [(Xup, Yup), (Xup, Y), (X, Yup), (Xdown, Yup), (Xup, Ydown), (Xdown, Y), (X, Ydown),
          (Xdown, Ydown)],
     % Compute distances to the home
     compute_distances(L0, L1),
     % Closer to home -> Earlier in the array
     sort(L1, L2),
     get_second(L2, Lheuristic).

% Same as previous + perception of covid Variant 2: count covid cells adjacent to each cell,
% use it as second sorting parameter
% VARIANT 2 BACKTRACKING 
get_adjacent_heuristic2((X, Y), Lheuristic) :-
     Yup is Y + 1,
     Ydown is Y - 1,
     Xup is X + 1,
     Xdown is X - 1,
     L0 = [(Xup, Yup), (Xup, Y), (X, Yup), (Xdown, Yup), (Xup, Ydown), (Xdown, Y), (X, Ydown),
          (Xdown, Ydown)],
     % Compute distances to the home
     compute_distances(L0, L1),
     % Compute number of adjacent covid cells for each cell - add to the array
     compute_covid_neighbours(L1, L2),
     % Closer to home -> Earlier in the array, Less covid neighbours -> Earlier in the array
     sort(L2, L3),
     get_third(L3, Lheuristic).

% Check whether two cells are adjacent or not
% Used for generating adjacent cells in backtracking - therefore, heuristics is applied
% Variants of perception matter - Variant 2 adds additional order of sorting - number of adjacent covid cells
is_adjacent_heuristic('Variant1', P1, P2) :-
     get_adjacent_heuristic(P1, L),
     member(P2, L).

is_adjacent_heuristic('Variant2', P1, P2) :-
     get_adjacent_heuristic2(P1, L),
     member(P2, L).

% Check whether cell is inside the map
inside_map((X, Y)) :-
     X >= 0,
     Y >= 0,
     map_xlimit(Xmax),
     X < Xmax,
     map_ylimit(Ymax),
     Y < Ymax.
 
% Check if cell is not in the covid zone
is_covid_free(Position) :-
     % For each of cell adjacent to positions, covid should be not in this cell
     forall(is_adjacent(Position, Cell), not(covid(Cell))).

% is_covid_free(Safety, Position)
% checking if Position if safe to go in
% If we wear a mask / been to doctor, it is safe anyway
is_covid_free(1, _) :- !.
% If we do not wear a mask / doctor, check if cell is in the covid zone
is_covid_free(0, Position) :-
     is_covid_free(Position).

% Is doctor or mask in the current position
is_doctor_or_mask(Position) :-
     doctor(Position); mask(Position).

% Return 0 as a second argument if neither doctor or mask are in Position cell
is_doctor_or_mask(Position, 0) :-
     not(doctor(Position)), not(mask(Position)).
% Return 1 as a second argument if there is doctor or mask in Position cell
is_doctor_or_mask(Position, 1) :-
     doctor(Position); mask(Position).


% BACKTRACKING FUNCTIONALITY
% BACKTRACKING FUNCTIONALITY
% BACKTRACKING FUNCTIONALITY
% VARIANTS 1 AND 2


% Heuristics for checking how many steps makes sense to set as upper bound for path length 
% 3 * N for <= 2 covid cells on the map
maximum_possible_steps(D) :-
     map_xlimit(Xmax),
     map_ylimit(Ymax),
     D is 3 * max(Xmax, Ymax).

% Try to update shortest path if possible. Fails if Path is longer than minimal one so far.
heuristics_shortest_path_set(Path) :-
     % get current minimal path
     min_path_length(MinLength),
     length(Path, Length),
     % compare lengths
     Length < MinLength,
     % update length if needed
     assertz(min_path_length(Length)),
     retract(min_path_length(MinLength)), !.

% Check whether current path is shorter than minimal so far - 
% if not, does not make sense to backtrack further
heuristics_shortest_path_check(Path) :-
     min_path_length(MinLength),
     length(Path, Length),
     Length < MinLength, !.

% Generate straight path from positoin to home, covid cells do not matter
% If at home, path is empty
generate_path_home(Position, []) :- 
     finish(Position), !.

% Make step in direction of home - Xnext, Ynext
generate_path_home((X, Y), [(Xnext, Ynext) | Path]) :- 
     finish((Xfinish, Yfinish)),
     (
     X < Xfinish, Xnext is X + 1;
     X = Xfinish, Xnext = Xfinish;
     X > Xfinish, Xnext is X - 1
     ), 
     (
     Y < Yfinish, Ynext is Y + 1;
     Y = Yfinish, Ynext = Yfinish;
     Y > Yfinish, Ynext is Y - 1
     ), 
     % Generate path from next cell to home
     generate_path_home((Xnext, Ynext), Path).

% Backtracking: If we are in home, update path if possbile - function terminates with success
dfs(_Perception, Position, Visited, [Position]) :-
     finish(Position),
     heuristics_shortest_path_set(Visited),
     !.

% Backtracking (heuristics) :
% If we are at mask / doctor, go straight to home
dfs(_Perception, Position, Visited, [Position | Path]) :-
     is_doctor_or_mask(Position),
     % generate path home directly
     generate_path_home(Position, Path), 
     % merge two paths
     append(Visited, Path, PathHome),
     % update path if possible
     heuristics_shortest_path_set(PathHome),
     !.

% Backtracking: find way home by brute-forcing all neighbours cells
dfs(Perception, Position, Visited, [Position | Path]) :-
     % Heuristics for path length
     heuristics_shortest_path_check(Visited),
     % Generate adjacent cells
     is_adjacent_heuristic(Perception, Position, NextPosition),
     % Check if inside the map
     inside_map(NextPosition),
     % Check they are not visited yet
     not(member(NextPosition, Visited)),
     % And free of covid, as we have no mask/doctor by now
     is_covid_free(NextPosition),
     % Start backtracking from NextPositon, updating Visited array
     dfs(Perception, NextPosition, [NextPosition | Visited], Path).

% Find some Path that leads home, or Fail if none exists
find_path_dfs(Perception, Path) :-
     start(Start),
     dfs(Perception, Start, [Start], Path).

% Enumerate all possible paths and choose minimal one
min_path_dfs(Perception, MinPath) :-    
     % Enumerate all paths that lead home
     bagof(Path, find_path_dfs(Perception, Path), L),
     % Last path in the array is shortest due to heuristcs applied (heuristics_shortest_path_set/check)
     last(L, MinPath).

% Print results of backtracking in the current environment
dfs(VariantCovidPerception) :-
     write('Backtracking '),write(VariantCovidPerception), write(' :'), nl,
     retractall(min_path_length(_)),
     maximum_possible_steps(X),
     % limit on path length
     assertz(min_path_length(X)),
     % measure time of execution
     statistics(walltime, _),
     (min_path_dfs(VariantCovidPerception, Path), !; Path = []),
     statistics(walltime, [_ | [ExecutionTime]]),
     % print results
     output_results(Path, ExecutionTime),
     % reset temporary variables
     retractall(min_path_length(_)).




%    A* FUNCIONALITY
%    A* FUNCIONALITY
%    A* FUNCIONALITY

% Return 1 as 3rd argument if cell is safe and 0 if not
% 1st argument - were we safe before this cell
% 2nd argument - are we safe in this cell (doctor / mask)
% result is logical AND of variables - same to max of integers (0/1).
determine_safety(X, Y, Z) :-
     Z is max(X, Y).

% Compute chessboard distance to the end for every cell in the array, return pair [distance, element]
compute_distances([], []) :- !.
compute_distances([H | T], [[D, H] | T1]) :-
     chessboard_distance(H, D), compute_distances(T, T1).

% Compute chessboard distance from given cell to the home
chessboard_distance((X, Y), Distance) :-
     finish((Xfinish, Yfinish)),
     Distance is max(abs(X - Xfinish), abs(Y - Yfinish)).

% Add neighbours of cell to the Priority Queue
% Empty - terminate
process_neighbours([], _, _, _, Distances, Predecessors, VerticesHeap, Distances, Predecessors, VerticesHeap) :- !.
process_neighbours([(X, Y) | T], Safety, DistanceToCurrent, (Xpred, Ypred, Spred),
                    Distances, Predecessors, VerticesHeap, DistancesNew, PredecessorsNew, VerticesHeapNew) :- 
     (
          % if path to the cell is shorter than existing one
          nth0_3d(X, Y, Safety, Distances, DistanceOld),
          DistanceToCurrent < DistanceOld,
          (
          % Add cell to the Priority Queue, update Predecessor and Distance
          chessboard_distance((X, Y), Heuristics),
          % Priority is Distance (Home, Cell) + Heuristics, which is chessboard distance
          Priority is DistanceToCurrent + Heuristics,
          update_value(X, Y, Safety, DistanceToCurrent, Distances, DistancesUpd),
          update_value(X, Y, Safety, (Xpred, Ypred, Spred), Predecessors, PredecessorsUpd),
          add_to_heap(VerticesHeap, Priority, (X, Y, Safety), VerticesHeapUpd) 
          ),
          !
          ;
          % If cell does not make sense to consider, just keep the parameters
          DistancesUpd = Distances, PredecessorsUpd = Predecessors, VerticesHeapUpd = VerticesHeap
     ),
     % process next neighbour anyway
     process_neighbours(T, Safety, DistanceToCurrent, (Xpred, Ypred, Spred),
               DistancesUpd, PredecessorsUpd, VerticesHeapUpd, DistancesNew, PredecessorsNew, VerticesHeapNew).


% Recursively build path to home
build_path((X, Y, _), _, [(X, Y)]) :- 
     start((X, Y)),
     !.
% build path for our predecessor, return path for us
build_path((X, Y, S), Predecessors, PathHome) :- 
     nth0_3d(X, Y, S, Predecessors, (Xpred, Ypred, Spred)),
     build_path((Xpred, Ypred, Spred), Predecessors, PathHomeHead),
     % path from (X, Y) to Home is path from Pred[X, Y] to Home + (X, Y) itself
     append(PathHomeHead, [(X, Y)], PathHome).

% A* case when at home
astar(VerticesHeap, _, Predecessors, PathHome) :-
     finish((Xfinish, Yfinish)),
     min_of_heap(VerticesHeap, _, (Xfinish, Yfinish, SafetyFinish)),
     % recursively with the help of Predecessors array build path home - A* terminates
     build_path((Xfinish, Yfinish, SafetyFinish), Predecessors, PathHome),
     !.

% A* iteration
astar(VerticesHeap, Distances, Predecessors, PathHome) :-
     % Get the minimum value on Priority Queue
     get_from_heap(VerticesHeap, _Priority, (X, Y, SafetyCurrent), VerticesHeap1),
     Position = (X, Y),
     % Get distance to current cell
     nth0_3d(X, Y, SafetyCurrent, Distances, DistanceToCurrent),
     DistanceToNext is DistanceToCurrent + 1,
     % Generate cells to add to queue:
     % adjacent cells
     get_adjacent(Position, AdjacentCells),
     % only those inside the map
     include(inside_map, AdjacentCells, AdjacentCells1),
     % include only safe cells
     is_doctor_or_mask(Position, SafetyPosition),
     determine_safety(SafetyCurrent, SafetyPosition, SafetyNext),
     include(is_covid_free(SafetyNext), AdjacentCells1, AdjacentCellsToProcess),
     % Add neighbours to the queue
     process_neighbours(AdjacentCellsToProcess, SafetyNext, DistanceToNext, (X, Y, SafetyCurrent),
               Distances, Predecessors, VerticesHeap1, DistancesNew, PredecessorsNew, VerticesHeapNew),
     % Start new iteration of A*
     astar(VerticesHeapNew, DistancesNew, PredecessorsNew, PathHome).
     
% Set up variables for A* and run it
min_path_astar(PathHome) :-
     % Set up Priority Queue
     empty_heap(Heap),
     start((Xstart, Ystart)),
     add_to_heap(Heap, 0, (Xstart, Ystart, 0), VerticesHeap),
     map_xlimit(Xmax), map_ylimit(Ymax), infinity(Inf),
     % Set up array of distances to start cell
     create_array(Xmax, Ymax, 2, Inf, Distances),
     update_value(Xstart, Ystart, 0, 0, Distances, Distances1),
     % Set up array of predecessors
     create_array(Xmax, Ymax, 2, -1, Predecessors),
     % Run A* - L contains the only value as A* is non-recursive algorithm and is true 
     % only for single Path found
     bagof(Path, astar(VerticesHeap, Distances1, Predecessors, Path), L),
     nth0(0, L, PathHome).

% Start A* algorithm, print it's result
astar() :-
     write('A* algorithm: '), nl,
     % measure time
     statistics(walltime, _),
     (min_path_astar(Path), !; Path = []),
     statistics(walltime, [_ | [ExecutionTime]]),
     % ouput results and map with path
     output_results(Path, ExecutionTime).



% TESTING AND OUTPUT FUNCTIONALITY
% TESTING AND OUTPUT FUNCTIONALITY
% TESTING AND OUTPUT FUNCTIONALITY

% Represent Execution time in convinient format
convert_time(ExecutionTime, Minutes, Seconds, MilliSeconds) :-
     Minutes is ExecutionTime // 60000, 
     Seconds is ExecutionTime // 1000 mod 60, 
     MilliSeconds is ExecutionTime mod 1000.

% Display result of the program execution
% Path is [] - not found -> FAIL
output_results([], ExecutionTime) :-
     % print execution time
     convert_time(ExecutionTime, M, S, Ms),
     write('Result: lose'), nl,
     write('Execution time: '), write(M), write(' min. '), 
     write(S), write(' sec. '),
     write(Ms), write(' ms.'),  
     nl, !.
% Path is not [] - print detailed information
output_results(Path, ExecutionTime) :-
     write('Path on the map:'), nl,
     % draw map with path on it
     draw_map(Path),
     % print result, path and it's length
     write('Result: win'), nl, 
     write('Shortest path is: '), nl, write(Path), nl,
     length(Path, Length),
     Steps is Length - 1,
     write('Steps done: '), write(Steps), nl,
     % print execution time
     convert_time(ExecutionTime, M, S, Ms),
     write('Execution time: '), write(M), write(' min. '), 
     write(S), write(' sec. '), 
     write(Ms), write(' ms.').

% Reset coordinates of covid, doctor, etc.
reset_environment() :- 
     retractall(map_xlimit(_)),
     retractall(map_ylimit(_)),
     retractall(start(_)),
     retractall(finish(_)),
     retractall(covid(_)),
     retractall(doctor(_)),
     retractall(mask(_)).

% Get everything in the environment
get_environment(Xlimit, Ylimit, Actor, Home, Doctor, Mask, Covid) :-
     map_xlimit(Xlimit), map_ylimit(Ylimit), 
     bagof(Xa, start(Xa), Actor), 
     bagof(Xb, finish(Xb), Home), 
     bagof(Xc, doctor(Xc), Doctor), 
     bagof(Xd, mask(Xd), Mask), 
     bagof(Xe, covid(Xe), Covid),
     !. 



% True if array consists of unique values
unique([]) :- !.
unique([H | T]) :-
     not(member(H, T)), 
     unique(T).

% print 1d array
print_array([]) :- !.
print_array([H | T]) :-
     write(H), print_array(T).

% print 2d array with bars
print_array_2d([]) :- !.
print_array_2d([H | T]) :-
     write("|"),
     print_array(H), 
     write("|"), nl,
     print_array_2d(T).

% Print array in a nice way
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

% Change coordinates for drawing
change_coordinates(_, _, [], []) :- !.
change_coordinates(Xlimit, Ylimit, [(X, Y) | T], [(X1, Y1) | T1]) :-
     Y1 is X, X1 is Ylimit - Y - 1,
     change_coordinates(Xlimit, Ylimit, T, T1).

% Draw map with Path on it
draw_map(Path) :-
     get_environment(Xlimit, Ylimit, Actor0, Home0, Doctor0, Mask0, Covid0),
     create_array(Xlimit, Ylimit, '.', Map),
     change_coordinates(Xlimit, Ylimit, Actor0, Actor1),
     change_coordinates(Xlimit, Ylimit, Home0, Home1),
     change_coordinates(Xlimit, Ylimit, Doctor0, Doctor1),
     change_coordinates(Xlimit, Ylimit, Mask0, Mask1),
     change_coordinates(Xlimit, Ylimit, Covid0, Covid1),
     change_coordinates(Xlimit, Ylimit, Path, Path1),
     % include only cell that are inside the map
     include(inside_map, Actor1, Actor2), update_values(Actor2, 'A', Map, Map1),
     include(inside_map, Home1, Home2), update_values(Home2, 'H', Map1, Map2),
     include(inside_map, Doctor1, Doctor2), update_values(Doctor2, 'D', Map2, Map3),
     include(inside_map, Mask1, Mask2), update_values(Mask2, 'M', Map3, Map4),
     include(inside_map, Covid1, Covid2), update_values(Covid2, 'C', Map4, Map5),
     update_values(Path1, 'o', Map5, MapFinal),
     print_array_2d(Xlimit, MapFinal).

% Draw map of the environment
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
