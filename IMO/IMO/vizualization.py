import matplotlib.pyplot as plt
import networkx as nx
import sys
from pathlib import Path

if (len(sys.argv) > 1):
    workDir = Path(sys.argv[1])
else:
    workDir = 'workdir'

# Read input coordinates
coordinates = {}
with open(workDir / "kroB200.tsp.txt", "r") as file:
    
    for line in file:
        if line.strip().split()[0].isdigit():
            node, x, y = map(int, line.strip().split())
            coordinates[node] = (x, y)
    

# Read output route
route1 = []
route2 = []
with open(workDir / "kroB200.tsp-IteratedLocalSearch2-solution.txt", "r") as file:
    r1=True
    for line in file:
        if line =="\n":
            r1=False
            continue

        if r1:
            x, y = map(int, line.strip().split())
            route1.append((x, y))
        else:
            x, y = map(int, line.strip().split())
            route2.append((x, y))
        

# Find labels to nodes in route
    for i in range(len(route1)):
        route1[i] = list(coordinates.keys())[list(coordinates.values()).index(route1[i])]
    for i in range(len(route2)):
        route2[i] = list(coordinates.keys())[list(coordinates.values()).index(route2[i])]

# Create a graph
G = nx.Graph()

# Add nodes
for node, (x, y) in coordinates.items():
    G.add_node(node, pos=(x, y))
pos = nx.get_node_attributes(G, 'pos')
nx.draw(G, pos, with_labels=True, node_size=20, node_color='skyblue', font_size=8)
# Add edges
for i in range(len(route1) - 1):
    start = route1[i]
    end = route1[i + 1]
    G.add_edge(start, end)
G.add_edge(route1[-1], route1[0])
#draw edges
nx.draw_networkx_edges(G, pos, width=1.0, edge_color='g')
#delete drawn edges
for i in range(len(route1) - 1):
    start = route1[i]
    end = route1[i + 1]
    G.remove_edge(start, end)
G.remove_edge(route1[-1], route1[0])

#Add edges from second path
for i in range(len(route2) - 1):
    start = route2[i]
    end = route2[i + 1]
    G.add_edge(start, end)
G.add_edge(route2[-1], route2[0])
#draw edges from second path
nx.draw_networkx_edges(G, pos, width=1.0, edge_color='r')
plt.show()