node(a). node(b).
node(c). node(d).
node(e). node(f).
node(g). node(h).
node(k). node(m).
node(p). node(s).

edge(a, b). edge(a, c). edge(a, g).
edge(b, a). edge(b, c).
edge(c, e). edge(c, d).
edge(d, f).
edge(e, g). edge(e, f). edge(e, h).
edge(f, k).
edge(g, c). edge(g, e).
edge(m, d).
edge(p, b). edge(p, d).

dfs(To, To, _, []) :-!.
dfs(From, To, VisitedNodes, [(From, X)|TailPath]):-
    edge(From, X), 
    not(member(X, VisitedNodes)),
    dfs(X, To, [From|VisitedNodes], TailPath).