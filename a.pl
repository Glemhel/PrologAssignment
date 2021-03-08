:- set_prolog_flag(answer_write_options,
                   [ quoted(true),
                     portray(true),
                     spacing(next_argument)
                   ]).
/**
 * Mikhail Rudakov
 * BS19-02
 * Prolog Programming Assignment
 * 
 **/

determine_safety(1, _, 1) :- !.
determine_safety(_, 1, 1) :- !.
determine_safety(0, 0, 0) :- !.

get_dist(_, [], _) :-
     print("bad"), false, !.

get_dist((Xcurrent, Ycurrent, SafetyCurrent), [H | T], DistanceToCurrent) :-
     H = (Xcurrent, Ycurrent, SafetyCurrent, DistanceToCurrent), !;
     get_dist((Xcurrent, Ycurrent, SafetyCurrent), T, DistanceToCurrent).

/**
 * map for the agent
 **/
map_xlimit(5).
map_ylimit(5).
start((0, 0)).
finish((3, 3)).
covid((1, 4)).
covid((6, 7)).
doctor((7, 1)).
mask((4, 4)).

get_adjacent((X, Y), L) :-
     Yup is Y + 1,
     Ydown is Y - 1,
     Xup is X + 1,
     Xdown is X - 1,
     L = [(Xup, Yup), (Xup, Y), (X, Yup), (Xdown, Yup), (Xup, Ydown), (Xdown, Y), (X, Ydown),
          (Xdown, Ydown)].

is_adjacent((X1, Y1), (X2, Y2)) :-
     X2 is X1 + 1, Y2 is Y1 + 1;
     X2 is X1 + 1, Y2 is Y1;
     X2 is X1, Y2 is Y1 + 1;
     X2 is X1 - 1, Y2 is Y1 + 1;
     X2 is X1 + 1, Y2 is Y1 - 1;
     X2 is X1 - 1, Y2 is Y1;
     X2 is X1, Y2 is Y1 - 1;
     X2 is X1 - 1, Y2 is Y1 - 1.

inside_map((X, Y)) :-
     map_xlimit(Xmax),
     map_ylimit(Ymax),
     X < Xmax,
     Y < Ymax,
     X >= 0,
     Y >= 0.

lengths([], [], _) :-!.
lengths([H | T], [[LenH, Ind] | T1], Ind) :-
     Ind1 is Ind + 1,
     length(H, LenH),
     lengths(T, T1, Ind1).

maximum_possible_steps(D) :-
     map_xlimit(Xmax),
     map_ylimit(Ymax),
     D is 2 * max(Xmax, Ymax) + 2.


is_covid_free(Position) :-
     forall(is_adjacent(Position, Cell), not(covid(Cell))).

is_doctor_or_mask(Position, 1) :-
     doctor(Position); mask(Position).

is_doctor_or_mask(Position, 0) :-
     not(doctor(Position)), not(mask(Position)).

generate_path((XCurrent, YCurrent), []) :- 
     finish((XCurrent, YCurrent)), !.

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



dfs((Xcurrent, Ycurrent, _), _, [(Xcurrent, Ycurrent)]) :-
     finish((Xcurrent, Ycurrent)),
     true, !.

% heuristics : once with a mask, go straight to home
dfs((Xcurrent, Ycurrent, 1), _, [(Xcurrent, Ycurrent) | Path]) :-
     generate_path((Xcurrent, Ycurrent), Path), !.

dfs((Xcurrent, Ycurrent, Safety), Visited, [(Xcurrent, Ycurrent) | Path]) :-
     % heuristics starts
     length(Visited, Length), maximum_possible_steps(Max), Length < Max, 
     % heuristics ends
     CurrentPosition = (Xcurrent, Ycurrent),
     NextPosition = (Xnext, Ynext),
     (is_covid_free(CurrentPosition); Safety = 1),
     is_adjacent(CurrentPosition, NextPosition),
     inside_map(NextPosition),
     not(member((Xnext, Ynext, Safety), Visited)),
     is_doctor_or_mask((Xcurrent, Ycurrent), Safety1),
     determine_safety(Safety, Safety1, Safe),
     (is_covid_free(NextPosition); Safe = 1),
     dfs((Xnext, Ynext, Safe), [(Xcurrent, Ycurrent, Safe) | Visited], Path).

find_path_dfs(Path) :-
     dfs((0, 0, 0), [(0, 0, 0)], Path).

min_path_dfs(MinPath) :-
     bagof(Path, find_path_dfs(Path), L),
     lengths(L, Lengths, 0),
     min_member([_, Index], Lengths),
     nth0(Index, L, MinPath).

% astar((Xfinish, Yfinish), VerticesHeap, VisitedList) :-
%      get_from_heap(VerticesHeap, Distance, [Xfinish, Yfinish, _], VerticesHeap1).

% process_neighbours([], _) :- !.
% process_neighbours([(Xcurrent, Ycurrent) | T], Safety, DistanceToCurrent, Distances, PredecessorList) :- 
%      (
%           not(member((Xcurrent, Ycurrent, Safety), Distances)), 
          
%           get_dist((Xcurrent, Ycurrent, Safety), Distances, DistanceOld),
%           DistanceToCurrent < DistanceOld
%      ),
% astar((Xfinish, Yfinish), VerticesHeap, Distances, PredecessorList) :-
%      get_from_heap(VerticesHeap, Priority, (Xcurrent, Ycurrent, SafetyCurrent), VerticesHeap1),
%      Position = (Xcurrent, Ycurrent),
%      get_dist((Xcurrent, Ycurrent, SafetyCurrent), Distances, DistanceToCurrent),
%      DistanceToCurrent1 is DistanceToCurrent + 1,
%      get_adjacent(Position, AdjacentCells),
%      include(inside_map, AdjacentCells, AdjacentCells1),
%      is_doctor_or_mask(Position, SafetyNext),
%      determine_safety(SafetyCurrent, SafetyNext, Safety),
%      ( Safety = 1, AdjacentCellsProcess = AdjacentCells1;
%        Safety = 0, inlclude(is_covid_free, AdjacentCells1, AdjacentCellsProcess)
%      ),
%      process_neighbours(AdjacentCellsProcess, Safety, DistanceToCurrent1, Distances, PredecessorList).
     

% find_way_astar(X, Y, Path) :-
%      true.