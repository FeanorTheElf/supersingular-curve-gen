from cmath import pi
from math import atan2
from modular_polynomial_sage import prime_modular_poly_mod_p

def vulcano_positions(G):
    cycle_basis = G.cycle_basis()
    crater = None
    if len(cycle_basis) == 1:
        crater = cycle_basis[0]
    elif len(G.vertices()) == 1:
        crater = G.vertices()
    else:
        # we find a maximal path - the middle vertices are the crater
        current = None
        for v in G.vertices():
            if len(G.neighbors(v)) == 1:
                current = v
                break
            
        found = { current: (0, None) }
        nodes = [current]
        while len(nodes) > 0:
            current = nodes.pop()
            for n in G.neighbors(current):
                if n not in found:
                    found[n] = (found[current][0] + 1, current)
                    nodes.append(n)

        current = None
        max_dist = -1
        for vertex in found:
            if found[vertex][0] > max_dist:
                max_dist = found[vertex][0]
                current = vertex

        path = [current]
        while found[current][1] is not None:
            path.append(found[current][1])
            current = found[current][1]

        l = len(path)
        if l % 2 == 0:
            crater = [path[l/2 - 1], path[l/2]]
        else:
            crater = [path[(l - 1)/2]]

    open = [*crater]
    positions = None
    if len(crater) >= 2:
        positions = {
            crater[j]: (cos(2*pi/len(crater) * j), sin(2*pi/len(crater) * j))
            for j in range(len(crater))
        }
    else:
        positions = { crater[0]: (0., 0.) }

    split_angle = 2*pi/len(crater)
    if len(crater) == 1:
        split_angle = 2*pi/len(G.neighbors(crater[0]))

    while len(open) != 0:
        current = open.pop()
        children = [w for w in G.neighbors(current) if w not in positions]
        for j in range(len(children)):
            child = children[j]
            if (positions[current][1], positions[current][0]) == (0., 0.):
                angle = 2*pi * (j - len(children)) / len(children)
            else:
                angle_fraction = 0
                if len(children) > 1:
                    angle_fraction = (j - len(children)/2 + 1/2) / (len(children) - 1)
                angle = atan2(positions[current][1], positions[current][0]) + split_angle * angle_fraction
            pos = (positions[current][0] + cos(angle), positions[current][1] + sin(angle))
            positions[child] = pos
        open += children

    return positions

p = 101
l = 2
F = GF(p**2, "z")
z = F.gen()
graph = Graph(loops = True)

for j in F:
    graph.add_vertex(j)
P = PolynomialRing(F, ['x'])

phi = prime_modular_poly_mod_p(l, p)
for j1 in F:
    for j2, _ in P(phi(y = j1)).roots():
        graph.add_edge(j1, j2, l)

i = 0
for component in graph.connected_components(sort = False):
    j = component[0]
    G = graph.subgraph(vertices = component)

    cycle_basis = G.cycle_basis()
    if len(cycle_basis) > 1:
        assert EllipticCurve(j.parent(), j = j).is_supersingular()
        G.plot().save("./plot_supersingular.png")
        continue
    else:
        positions = vulcano_positions(G)
        G.set_pos(positions)
        G.plot().save("plot_" + str(i) + ".png")
        i += 1
        if i > 100:
            break
