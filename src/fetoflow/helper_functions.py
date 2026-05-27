import networkx as nx
import numpy as np
from .geometry_utils import  calcLength
def getRadii(G: nx.digraph):
    """
    Returns all radii of all vessels in the Arterial tree as a dictionary:
    (nodeIn, nodeOut) : radius
    """
    radii = {}
    for u, v in G.edges():
        radii[(u, v)] = G[u][v]["radius"]
    return radii


def getRadius(G: nx.digraph, nodes: tuple = None):
    if nodes is not None:
        u, v = nodes
        return G[u][v]["radius"]
    else:
        raise ("Enter node values!")


def getNode(G, nodeID):
    return G.nodes()[nodeID]


def getEdgeData(G, edge):
    u, v = edge
    return G.get_edge_data(u, v)


def getVesselLength(G, initialVessel: tuple):
    """
    Check this
    """
    u, v = initialVessel
    length = G[u][v]["length"]
    out = G.out_edges(v)
    if len(out) == 1:
        length += getVesselLength(G, (v, out[0]))
    else:
        return length


def getNumVessels(G):
    """
    Returns number of distinct vessels in graph
    """
    successors = nx.dfs_sucessors(G, source=1)  # check if 0?
    count = 0
    if successors.get(1) == 1:
        count += 1
    for __, nodesOut in successors.items():
        unique = len(nodesOut)
        if unique > 1:
            count += unique
    return count

def display_tree_measurements(G):
    """
       Displays the measurements of the tree including resistance of tree, Inlet and outlet pressures and flows, Max Strahler order.
       TODO: To be added
    """
    return
def export_as_numpy(G,output_form="reduced"):
    """
    Converts a complete digraph to numpy
        if output_form == "full":
        node_export = np.ndarray(shape=(G.number_of_nodes(),4))
        for node,node_data in G.nodes(data=True):
            x,y,z = node_data["x"],node_data["y"],node_data["z"]
            node_export[node,:] = node,x,y,z
        print(node_export)
    ## We need to convert each field stored in the digraph to a numpy array
    elif output_form == "reduced":
        pass
    else:
        print("Use full or reduced")


    """
    print("This function is under development!")
    return 

def calc_vessel_volume(G, vessel_type = 'all'):
    total_volume = 0
    if vessel_type == 'all':
        for u,v,data in G.edges(data = True):
            volume = data['length'] * np.pi * (data['radius'] **2)
            total_volume += volume
    if vessel_type == 'artery':
        for u,v,data in G.edges(data = True):
            if data['vessel_type'] == 'artery':
                volume = data['length'] * np.pi * (data['radius'] **2)
                total_volume += volume
    return total_volume


def build_inlet_map(G, inlets):
    """
    Returns {node: inlet_node} for every node in G.
    Works for any number of disconnected trees.
    """
    node_to_inlet = {}

    for inlet in inlets:
        # BFS from each inlet, following directed edges
        queue = deque([inlet])
        node_to_inlet[inlet] = inlet
        while queue:
            node = queue.popleft()
            for successor in G.successors(node):
                if successor not in node_to_inlet:  # avoid revisiting
                    node_to_inlet[successor] = inlet
                    queue.append(successor)

    return node_to_inlet