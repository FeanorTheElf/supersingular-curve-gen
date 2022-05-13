from modular_polynomial_sage import prime_modular_poly_mod_p

p = next_prime(100)
F = GF(p**2, "z")
graph = Graph(loops = True)
for j in F:
    graph.add_vertex(j)
P = PolynomialRing(F, ['x'])

for l in [3]:
    phi = prime_modular_poly_mod_p(l, p)
    for j1 in F:
        for j2, _ in P(phi(y = j1)).roots():
            graph.add_edge(j1, j2, l)

i = 0
for component in graph.connected_components(sort = False):
    j = component[0]
    if j**p not in component:
        G = graph.subgraph(vertices = component)
        graph_plot = G.plot(color_by_label = True)
        graph_plot.save("./plot" + str(i) + ".png")
    i += 1
    if i > 100:
        break
