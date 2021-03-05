% Как то раз случай свёл в купе астронома, поэта , 
% прозаика и драматурга. Это были Алексеев, Борисов, Константинов и Дмитриев. 
% Оказалось, что каждый из них взял с собой книгу написанную одним из пассажиров этого купе. 
% Алексеев и Борисов углубились в чтение предварительно обменявшись книгами. 
% Поэт читал пьесу, прозаик — очень молодой человек, выпустивший свою книгу, говорил что он никогда и ни чего не читал по астрономии. 
% Борисов купил одно из произведений Дмитриева. Никто из пассажиров не читал свои книги. Что читал каждый из них, кто кем был?

man(alekseev).
man(borisov).
man(konstantinov).
man(dmitriev).

profession(astronomy).
profession(poetry).
profession(prose).
profession(piece).


check([]) :- !.
check([passenger(_, Y, Z, W) | T]) :-
     not(Y=W), not(Z=W), check(T).

unique([]) :- !.
unique([H | T]) :-
    not(member(H, T)), unique(T).

solve(Solve) :-
    Solve = [passenger(alekseev, Xread, Xbuy, Xprofession), passenger(borisov, Yread, Ybuy, Yprofession), 
    passenger(konstantinov, Zread, Zbuy, Zprofession), passenger(dmitriev, Wread, Wbuy, Wprofession)],
    profession(Xprofession), profession(Yprofession), profession(Zprofession), profession(Wprofession),
    unique([Xprofession, Yprofession, Zprofession, Wprofession]),
    profession(Xread), profession(Yread), profession(Zread), profession(Wread),
    unique([Xread, Yread, Zread, Wread]),
    profession(Xbuy), profession(Ybuy), profession(Zbuy), profession(Wbuy),
    unique([Xbuy, Ybuy, Zbuy, Wbuy]),
    % поэт читал пьесу
    member(passenger(_, piece, _, poetry), Solve),
    % прозаик не читал астрономию
    not(member(passenger(_, astronomy, _, prose), Solve)),
    not(member(passenger(_, _, astronomy, prose), Solve)),
    % Дмитриев - не прозаик
    not(member(passenger(dmitriev, _, _, prose), Solve)),
    % Борисов купил одно из произведений Дмитриева
    member(passenger(borisov, _, DmirtievProfession, _), Solve),
    member(passenger(dmitriev, _, _, DmirtievProfession), Solve),
    % Алексеев и Борисов поменялись книгами
    member(passenger(borisov, AlekseevBuy, BorisovBuy, _), Solve),
    member(passenger(alekseev, BorisovBuy, AlekseevBuy, _), Solve),
    % Никто из пассажиров не читал свои книги
    check(Solve).
