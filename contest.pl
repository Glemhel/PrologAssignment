% Три друга заняли первое, второе и третье места в соревнованиях универсиады. 
% Друзья — разной национальности, зовут их по-разному и любят они разные виды спорта.
% Майкл предпочитает баскетбол и играет лучше чем американец.
% Израильтянин Саймон играет лучше теннисиста. Игрок в крикет занял первое место.
% Кто является австралийцем? Каким видом спорта занимается Ричард?
man(mike).
man(simon).
man(richard).

nation(us).
nation(au).
nation(iz).

sport(tennis).
sport(basketball).
sport(cricket).

place(p1).
place(p2).
place(p3).

unique([]) :- !.
unique([ H | T]) :-
    not(member(H, T)), unique(T).

better(person(X1, Y1, Z1, W1), person(X2, Y2, Z2, W2)) :-
    W1 = p1;W1=p2, W2=p3.

solve(X) :-
    X = [person(mike, MikePlay, MikeFrom, MikePlace),
        person(simon, SimonPlay, SimonFrom, SimonPlace),
        person(richard, RichardPlay, RichardFrom, RichardPlace)],
    sport(MikePlay), sport(SimonPlay), sport(RichardPlay),
    nation(MikeFrom), nation(SimonFrom), nation(RichardFrom),
    place(MikePlace), place(SimonPlace), place(RichardPlace),
    unique([MikeFrom, SimonFrom, RichardFrom]),
    unique([MikePlace, SimonPlace, RichardPlace]),
    unique([MikePlay, SimonPlay, RichardPlay]),
    member(person(mike, basketball, _, _ ), X),
    not(member(person(mike, _, us, _), X)),
    member(person(simon, _, iz, _), X),
    not(member(person(simon, tennis, _, _), X)),
    member(person(_, cricket, _, p1), X),
    better(person(simon, _, _, SimonPlace), person(_, tennis, _, _)),
    better(person(mike, _, _, MikePlace), person(_, _, us, _)).
