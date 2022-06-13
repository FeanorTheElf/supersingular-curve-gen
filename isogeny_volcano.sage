from cmath import pi
from math import atan2
from modular_polynomial_sage import prime_modular_poly_mod_p

def plot_vulcano(G):
    cycle_basis = G.cycle_basis()
    if len(cycle_basis) == 1:
        crater = cycle_basis[0]
        open = [*crater]
        positions = {
            crater[j]: (cos(2*pi/len(crater) * j), sin(2*pi/len(crater) * j))
            for j in range(len(crater))
        }
        while len(open) != 0:
            current = open.pop()
            children = [w for w in G.neighbors(current) if w not in positions]
            for j in range(len(children)):
                child = children[j]
                angle = atan2(positions[current][1], positions[current][0]) + 2*pi/len(crater) * (j - len(children)/2)/len(children)
                pos = (positions[current][0] + cos(angle), positions[current][1] + sin(angle))
                positions[child] = pos
            open += children
        G.set_pos(positions)

        return G.plot(color_by_label = True)
    else:
        return G.plot(layout = 'tree', color_by_label = True)

p = 101
F = GF(p**2, "z")
z = F.gen()
graph = Graph(loops = True)
for j in F:
    graph.add_vertex(j)
P = PolynomialRing(F, ['x'])

for l in [2]:
    phi = prime_modular_poly_mod_p(l, p)
    print(phi)
    for j1 in F:
        for j2, _ in P(phi(y = j1)).roots():
            graph.add_edge(j1, j2, l)

i = 0
for component in graph.connected_components(sort = False):
    if 61 * z + 16 in component:
        j = component[0]
        G = graph.subgraph(vertices = component)

        cycle_basis = G.cycle_basis()
        if len(cycle_basis) > 1:
            assert EllipticCurve(j.parent(), j = j).is_supersingular()
            G.plot().save("./plot_supersingular.png")
            continue

        plot_vulcano(G).save("./plot" + str(i) + ".png")
        i += 1
        if i > 100:
            break
