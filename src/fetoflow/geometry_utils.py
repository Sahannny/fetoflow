import networkx as nx
import numpy as np
from warnings import warn
from collections import deque

def create_geometry(
        nodes,
        elements,
        inlet_radius=1.8,
        strahler_ratio_arteries=1.38,
        outlet_vein_radius=4.0,
        strahler_ratio_veins=1.46,
        arteries_only=False,
        fields=None,
        single_umbilical_vein=True, #TODO: Take out
        anastomosis_edge=4, #TODO: Take out
        default_mu=0.33600e-02,
        default_hematocrit=0.45,
):
    G = nx.DiGraph()
    num_terminal_arterial_nodes = 0
    # Read in Nodes and Edges
    for node_id, coordinates in nodes.items():
        G.add_node(node_id, x=coordinates[0], y=coordinates[1], z=coordinates[2])
    for edge_id, (node_from, node_to) in enumerate(elements):
        length = calcLength(G, node_from, node_to)
        G.add_edge(
            node_from,
            node_to,
            edge_id=edge_id,
            resistance=None,
            length=length,
            radius=None,
            strahler=None,
            branch_number=None,
            vessel_type="artery",
            bifurcation_angle=0,
            mu=default_mu,
            hematocrit=default_hematocrit,
            viscosity_factor=1,
        )
    G = assign_radii_files(G, fields)
    # Find all input nodes (to ensure we give every element a strahler ordering)
    input_nodes = []
    for node in G.nodes():
        if len(G.in_edges(node)) == 0:
            input_nodes.append(node)
    # Add artery radii via strahler ordering
    max_strahler = 0

    G = update_strahler_nonrecursive(G)

    max_strahler = np.max([data["strahler"] for _, _, data in G.edges(data=True)])

    if fields and fields.get("radius"):
        for u, v in G.edges():

            if G[u][v]["radius"] is None:
                G[u][v]["radius"] = 0.0001
                '''elem_strahler = G[u][v]["strahler"]
                # need to update R according to the last specified radius in the subtree
                radius_found = False
                inlet_radius_updated = inlet_radius
                out_node = u
                sub_tree_strahler = max_strahler
                while not radius_found:
                    edges_to_check = list(G.in_edges(out_node))
                    # should only be one edge coming in
                    if len(edges_to_check) == 0:
                        radius_found = True  # no previously set radii in vessel's predecessors
                    else:
                        in_node = edges_to_check[0][0]
                        current_rad = G[in_node][out_node]["radius"]
                        if current_rad is not None and current_rad != 0:
                            inlet_radius_updated = current_rad
                            radius_found = True
                            sub_tree_strahler = G[in_node][out_node]["strahler"]
                        else:
                            out_node = in_node
                G[u][v]["radius"] = inlet_radius_updated * strahler_ratio_arteries ** (elem_strahler - sub_tree_strahler)'''
    else:
        for u, v in G.edges():
            elem_strahler = G[u][v]["strahler"]
            G[u][v]["radius"] = inlet_radius * strahler_ratio_arteries ** (elem_strahler - max_strahler)

    # Create the venous mesh
    if not arteries_only:
        num_artery_nodes = G.number_of_nodes()  # use for scaling to keep numeric based indexing
        num_artery_edges = G.number_of_edges()
        # Get terminal nodes
        terminal_nodes = [node for node, out_degree in G.out_degree() if out_degree == 0]
        num_terminal_arterial_nodes = len(terminal_nodes)
        venous_mesh = create_venous_mesh(
            G,
            num_artery_nodes,
            num_artery_edges,
            num_terminal_arterial_nodes,
            outlet_vein_radius,
            strahler_ratio_veins,
            max_strahler,
            single_umbilical_vein
        )
        # list of artery terminal node indices
        # add venous mesh to graph
        assert max(G.nodes) < min(venous_mesh.nodes), "Venous mesh node ids overlap with arterial node ids."
        G = nx.compose(G, venous_mesh)
        # Add edges to each terminal node with equivalent capillary network resistance
        edge_id_tracker = len(elements)  # Use this to easily increment edge_id correctly for the new edges I am adding.
        for i, arterial_node in enumerate(terminal_nodes):
            venous_node = arterial_node + num_artery_nodes

            assert len(G.in_edges(
                arterial_node)) == 1, f"Terminal artery node has {len(G.in_edges(arterial_node))} entering it, should be 1 only."
            assert len(G.out_edges(
                venous_node)) == 1, f"Terminal vein node has {len(G.out_edges(venous_node))} exiting it, should be 1 only."

            G.add_edge(
                arterial_node,
                venous_node,
                edge_id=edge_id_tracker + i,
                resistance=None,
                length=None,
                radius=0.0,
                strahler=0.0,
                branch_number=0,
                vessel_type="capillary_equivalent",
                mu=default_mu,
                hematocrit=default_hematocrit,
                viscosity_factor=1,

            )
            # Update length of these vessels - we calculate the resistance in 'calculate_resistance()' function.
            G[arterial_node][venous_node]["length"] = calcLength(G, arterial_node, venous_node)
            # Leave radius as 0 - this allows visualisations not including the capillary networks,
            #  which would only show a single vessel anyway (not the tree of the intermediate and terminal villi).


    return G


def create_venous_mesh(
        G,
        num_artery_nodes,
        num_artery_edges,
        num_terminal_arterial_nodes,
        outlet_vein_radius,
        strahler_ratio_veins,
        max_strahler,
        single_umbilical_vein
):
    venous_mesh = G.copy()
    # first n/2 nodes = artery nodes going down.
    # Next n/2 nodes = vein nodes in same order artery nodes were (i.e. likely inlet first, getting smaller as we go).
    nx.relabel_nodes(venous_mesh, lambda node_id: node_id + num_artery_nodes, copy=False)  # 0-based indexing works here
    # Update radii
    edges_to_remove = [
        (u, v) for u, v, data in venous_mesh.edges(data=True) if data["vessel_type"] == "anastomosis"
    ]

    G.remove_edges_from(edges_to_remove)
    outlets = [n for n, d in venous_mesh.in_degree() if d == 0]
    outlet_edges = [(u, v, data) for u, v, data in venous_mesh.edges(data=True) if venous_mesh.in_degree(u) == 0]
    edge_to_inlet = build_edge_inlet_map(venous_mesh, outlets)
    max_edge_id = max(data['edge_id'] for u, v, data in venous_mesh.edges(data=True))
    if single_umbilical_vein:
        inlet1 = venous_mesh.nodes[outlets[0]]
        inlet2 = venous_mesh.nodes[outlets[1]]

        # 2. Midpoint for junction
        mid = {
            'x': (inlet1['x'] + inlet2['x']) / 2,
            'y': (inlet1['y'] + inlet2['y']) / 2,
            'z': (inlet1['z'] + inlet2['z']) / 2,
        }
        z_value = (venous_mesh.nodes[outlet_edges[0][1]]['z'] + venous_mesh.nodes[outlet_edges[1][1]]['z'])/2
        venous_mesh.nodes[outlets[0]].update({'x': mid['x'], 'y': mid['y'], 'z': z_value})

        venous_mesh.nodes[outlets[1]].update({'x': mid['x'], 'y': mid['y'], 'z': z_value + (venous_mesh.edges[outlet_edges[0][0], outlet_edges[0][1]]['length'])})
        edge_data1 = venous_mesh.edges[outlet_edges[0][0], outlet_edges[0][1]]
        edge_data2 = venous_mesh.edges[outlet_edges[1][0], outlet_edges[1][1]]
        venous_mesh.remove_edge(outlet_edges[1][0], outlet_edges[1][1])
        venous_mesh.add_edge(outlet_edges[0][0], outlet_edges[1][1], **edge_data2)
        venous_mesh[outlet_edges[0][0]][outlet_edges[1][1]]['length'] = calcLength(venous_mesh, outlet_edges[0][0],
                                                                                   outlet_edges[1][1])
        venous_mesh.add_edge(
            outlet_edges[1][0],
            outlet_edges[0][0],
            edge_id= max_edge_id + 1,
            resistance=None,
            length=None,
            radius=0.0,
            strahler=0.0,
            branch_number=0,
            vessel_type="vein",
            default_mu=0.33600e-02,
            default_hematocrit=0.45,
            viscosity_factor=1,
        )
        venous_mesh[outlet_edges[1][0]][outlet_edges[0][0]]['length'] = calcLength(venous_mesh, outlet_edges[1][0], outlet_edges[0][0])
        if venous_mesh[outlet_edges[0][0]][outlet_edges[0][1]]['strahler'] ==  venous_mesh[outlet_edges[0][0]][outlet_edges[1][1]]['strahler']:
            venous_mesh[outlet_edges[1][0]][outlet_edges[0][0]]['strahler'] = venous_mesh[outlet_edges[0][0]][outlet_edges[0][1]]['strahler'] +1
        else:
            venous_mesh[outlet_edges[1][0]][outlet_edges[0][0]]['strahler'] = np.max([venous_mesh[outlet_edges[0][0]][outlet_edges[0][1]]['strahler'],  venous_mesh[outlet_edges[0][0]][outlet_edges[1][1]]['strahler']])
    outlets = [n for n, d in venous_mesh.in_degree() if d == 0]
    outlet_edges = [(u, v, data) for u, v, data in venous_mesh.edges(data=True) if venous_mesh.in_degree(u) == 0]
    edge_to_inlet = build_edge_inlet_map(venous_mesh, outlets)
    for u, v in venous_mesh.edges():
        # Edge ids for veins are after capillaries
        venous_mesh[u][v]["edge_id"] = venous_mesh[u][v]["edge_id"] + num_artery_edges + num_terminal_arterial_nodes
        venous_mesh[u][v]["vessel_type"] = "vein"
        # Strahler ordering for veins
        elemStrahler = venous_mesh[u][v]["strahler"]
        if edge_to_inlet[u,v] == outlet_edges[0][0]:
            venous_mesh[u][v]["radius"] = outlet_vein_radius * strahler_ratio_veins ** (elemStrahler - outlet_edges[0][2]['strahler'])
        elif edge_to_inlet[u,v] == outlet_edges[1][0]:
            venous_mesh[u][v]["radius"] = outlet_vein_radius * strahler_ratio_veins ** (elemStrahler - outlet_edges[1][2]['strahler'])

    venous_mesh = venous_mesh.reverse()

    return venous_mesh


def update_strahlers(G, node_in, node_out):
    # Input the input node(s) -  or call on each input node. Updates the strahler field.
    # This is a recursive function and will be slow for now.
    child_edges = G.out_edges(node_out)
    # # Base case: No children, strahler == 1.
    if len(child_edges) == 0:
        G[node_in][node_out]["strahler"] = 1
        return G

    max_child_strahler = 0
    max_child_strahler_count = 0

    for u, v in child_edges:
        if G[u][v]["strahler"] is None:
            G = update_strahlers(G, u, v)  # Returns graph object with strahler updated.

        if G[u][v]["strahler"] > max_child_strahler:
            max_child_strahler = G[u][v]["strahler"]
            max_child_strahler_count = 1
        elif G[u][v]["strahler"] == max_child_strahler:
            max_child_strahler_count += 1

    assert max_child_strahler_count > 0
    # Actually assign this arc's strahler now.
    if max_child_strahler_count == 1:
        G[node_in][node_out]["strahler"] = max_child_strahler
    else:
        G[node_in][node_out]["strahler"] = max_child_strahler + 1  # 2 arcs coming in with same max value

    return G


def Strahler_numbering(di_graph: nx.DiGraph, inlet):
    di_graph = di_graph.copy()
    strahler_mapping = {}
    post_order_nodes = nx.dfs_postorder_nodes(di_graph, source=inlet)
    for node in post_order_nodes:
        children = nx.dfs_successors(di_graph, source=node, depth_limit=1)
        if len(children) == 0:
            strahler_mapping[node] = {'strahler_order': 1}
        else:
            children = children[node]
            strahler_set = []
            for child in children:
                strahler_set.append(strahler_mapping[child]['strahler_order'])

            strahler_set.sort()
            if len(strahler_set) > 1:
                if strahler_set[-1] == strahler_set[-2]:
                    strahler_mapping[node] = {'strahler_order': strahler_set[-1] + 1}
                else:
                    strahler_mapping[node] = {'strahler_order': strahler_set[-1]}
            else:
                strahler_mapping[node] = {'strahler_order': strahler_set[-1]}

    nx.set_node_attributes(di_graph, strahler_mapping)

    return di_graph


def set_edge_strahler(G):
    for u, v in G.edges():
        strahler_list = [G.nodes[u]["strahler_order"], G.nodes[v]["strahler_order"]]

        G[u][v]["strahler"] = np.min(strahler_list)
    return G


def update_strahler_nonrecursive(G):
    '''order = list(nx.topological_sort(G))
    order.reverse()
    node_strahler = {}
        # 2. Compute node Strahler numbers

    for n in order:
        children = list(G.successors(n))

        if len(children) == 0:
            node_strahler[n] = 1
        else:
            child_orders = [node_strahler[c] for c in children]
            max_order = max(child_orders)

            if child_orders.count(max_order) >= 2:
                node_strahler[n] = max_order + 1
            else:
                node_strahler[n] = max_order

    # 3. Assign edge Strahler numbers
    for u, v in G.edges():
        G.edges[u, v]["strahler"] = node_strahler[v]
    '''

    node_list = list(nx.topological_sort(G))
    node_list.reverse()
    for nodes in node_list:
        # Check out edges to see if bifurcation or terminal
        if G.out_degree(nodes) == 0:  # Terminal
            if G.in_degree(nodes) == 2:
                print("bifurcation/anastomosis")  # dont think this should happen
            elif G.in_degree(nodes) == 1:  # should be almost always
                edge = list(G.in_edges(nodes))[0]  # Upstream edge of node
                G[edge[0]][edge[1]]["strahler"] = 1
        elif G.out_degree(nodes) == 1:  # Normal edge
            if G.in_degree(nodes) == 2:
                print("bifurcation/anastomosis")  # dont think this should happen
            elif G.in_degree(nodes) == 1:
                edge_in = list(G.in_edges(nodes))[0]  # Upstream edge of node
                edge_out = list(G.out_edges(nodes))[0]  # Downstream edge of node
                G[edge_in[0]][edge_in[1]]["strahler"] = G[edge_out[0]][edge_out[1]]["strahler"]
            elif G.in_degree(nodes) == 0:
                print("inlet detected")

                none_edges = [(u, v) for u, v, data in G.edges(data=True)
                              if data.get('strahler') == None]

                if none_edges:
                    for u, v in none_edges:
                        edge_out = list(G.out_edges(v))[0]  # Downstream edge of node
                        G[u][v]["strahler"] = G[edge_out[0]][edge_out[1]]["strahler"]
        elif G.out_degree(nodes) == 2:  # Bifurcation
            if G.in_degree(nodes) == 2:
                ValueError(
                    "Fetoflow cannot assign strahler orders for two incoming edges and two outgoing edges. Check the tree structure to avoid loops")
            elif G.in_degree(nodes) == 1:
                edge_in = list(G.in_edges(nodes))[0]  # Upstream edge of node
                edge_out_1 = list(G.out_edges(nodes))[0]  # Downstream first edge of node
                edge_out_2 = list(G.out_edges(nodes))[1]  # Downstream second edge of node
                if G[edge_out_1[0]][edge_out_1[1]]["vessel_type"] == 'anastomosis':
                    G[edge_in[0]][edge_in[1]]["strahler"] = G[edge_out_2[0]][edge_out_2[1]]["strahler"]
                elif G[edge_out_2[0]][edge_out_2[1]]["vessel_type"] == 'anastomosis':
                    G[edge_in[0]][edge_in[1]]["strahler"] = G[edge_out_1[0]][edge_out_1[1]]["strahler"]
                elif G[edge_out_1[0]][edge_out_1[1]]["strahler"] == G[edge_out_2[0]][edge_out_2[1]]["strahler"]:
                    G[edge_in[0]][edge_in[1]]["strahler"] = G[edge_out_2[0]][edge_out_2[1]]["strahler"] + 1
                else:
                    G[edge_in[0]][edge_in[1]]["strahler"] = np.max(
                        [G[edge_out_2[0]][edge_out_2[1]]["strahler"], G[edge_out_1[0]][edge_out_1[1]]["strahler"]])
    return G


def calcLength(G, u, v):
    return np.sqrt(np.sum([(G.nodes[u][coord] - G.nodes[v][coord]) ** 2 for coord in ["x", "y", "z"]]))  # mm to m!
    # TODO: Fix unit conversions and makr work for anything etc. i.e. specify units somewhere as an input argument at the start


def create_anastomosis(G, node_from, node_to, radius=None, mu=0.33600e-02):
    # NOTE HERE: RADIUS IS IN mm!!!!!
    # Todo: make sure this is clear.
    # TODO: DIGRAPH STUFF???. No we can just have a negative flow along the anastomosis.
    # TODO: probably write this as 2 separate functions, this one which is in here with the graph stuff, and one which is user-facing and calls this one
    u = node_from - 1
    v = node_to - 1  # Update from 1- to 0-based indexing\
    max_child_radius = 0.0
    max_child_strahler = 0.0
    for node in (u, v):
        for child_u, child_v in G.out_edges(node):
            child_radius = G[child_u][child_v]["radius"]
            child_strahler = G[child_u][child_v]["strahler"]
            if not child_radius or not child_strahler:
                raise ValueError(
                    "Strahler values and radii have not been set yet. make sure you do this before creating an anastomosis.")
            max_child_radius = max(max_child_radius, child_radius)
            max_child_strahler = max(max_child_strahler, child_strahler)
    if u not in G:
        raise ValueError(
            f"Node {node_from} (ipnode indexing)/Node {u} (networkX indexing) does not exist in the networkX graph. Perhaps you need to call create_geometry() first?"
        )
    if v not in G:
        raise ValueError(
            f"Node {node_to} (ipnode indexing)/Node {v} (networkX indexing) does not exist in the networkX graph. Perhaps you need to call create_geometry() first?"
        )
    if u == v:
        raise ValueError(
            f"Anastomosis cannot connect the same node to itself. Node number is {node_from} (ipnode indexing)/{u} (networkX indexing).")

    # only exists as aterial connection
    G.add_edge(
        u,
        v,
        edge_id=G.number_of_edges(),
        resistance=0.0,
        length=None,
        radius=0.0,
        strahler=0.0,
        vessel_type="anastomosis",
        mu=mu,
        hematocrit=0.45,  # TODO PARAMETERISE
        viscosity_factor=1
    )
    # Old implementation:
    # - defines a radius of the anastomosis which is used
    # Our alternative:
    # - Provide a warning if this happens, but just use the maxiumum radii of the vessels leaving the nodes connecting the anastomosis.
    # TODO: Confirm this is fine implementation

    G[u][v]["length"] = calcLength(G, node_from, node_to)

    # For strahler, just take the max of any child strahler.


    if radius:
        if not (isinstance(radius, int) or isinstance(radius, float)):
            raise ValueError(f"Hyrtl anastomosis radius is invalid type {type(radius)}. Valid types are float or int")

        G[u][v]["radius"] = radius # mm to m!
    else:
        G[u][v]["radius"] = max_child_radius
    # Strahler
    G[u][v]["strahler"] = max_child_strahler  # The code will previously break if strahlers have not been already set.
    # Resistance calculation: #TODO make sure it updates properly if we have other viscosities etc.
    # Note: if calcu alte_resistance() is called after this function, the result will be overwritten. This is probably the order we want:
    # - calculate_geometry()
    # - create_venous_mesh()
    # - create_anatsomosis()
    # - calculate_resistance()
    G[u][v]["resistance"] = 8 * mu * G[u][v]["length"] / (np.pi * G[u][v]["radius"] ** 4)

    return G


def update_geometry_with_pressures_and_flows(G, pressures, flows, edge_id_attr="edge_id"):
    # Check for sizes of arrays
    if len(pressures) != G.number_of_nodes():
        warn(
            "Number of Nodes in Digraph does not match number of pressures. Skipping adding pressures/flows to digraph. Note: Double Check outputs.")
        return G
    elif len(flows) != G.number_of_edges():
        warn(
            "Number of Edges in Digraph does not match number of flows. Skipping adding pressures/flows to digraph. Note: Double Check outputs.")
        return G

    if pressures is not None:
        nx.set_node_attributes(G, pressures, "pressure")

    if flows is not None:
        # Case 1: keyed by (u, v)
        if all(isinstance(k, tuple) and len(k) == 2 for k in flows):
            nx.set_edge_attributes(G, flows, "flow")
        else:
            # Case 2: keyed by edge_id
            for u, v, data in G.edges(data=True):
                eid = data.get(edge_id_attr)
                if eid in flows:
                    G.edges[u, v]["flow"] = flows[eid]

    return G


#   for node_id in G.nodes():
#       G.nodes[node_id]["pressure"] = pressures.loc[node_id]["pressure"]
#   for u, v in G.edges():
#       G[u][v]["flow"] = flows.loc[G[u][v]["edge_id"]]["flow"]
#   return G

def calculate_branching_angles(G):
    # double check...
    for n in G.nodes():
        out_edges = list(G.out_edges(n))
        in_edges = list(G.in_edges(n))
        if len(out_edges) > 1:
            for __, out_node in out_edges:
                out_vec = np.array([G.nodes[n][coord] - G.nodes[out_node][coord] for coord in ["x", "y", "z"]])
                out_norm = out_vec / np.linalg.norm(out_vec)
                in_vec = np.array([G.nodes[in_edges[0][0]][coord] - G.nodes[n][coord] for coord in ["x", "y", "z"]])
                in_norm = in_vec / np.linalg.norm(in_vec)
                dot = np.clip(in_norm @ out_norm, -1, 1)
                G[n][out_node]["bifurcation_angle"] = np.pi - np.arccos(dot)  # phi_j from Mynard

        elif len(in_edges) > 1:
            for in_node, __ in in_edges:
                in_vec = np.array([G.nodes[in_node][coord] - G.nodes[n][coord] for coord in ["x", "y", "z"]])
                in_norm = in_vec / np.linalg.norm(in_vec)
                out_vec = np.array([G.nodes[n][coord] - G.nodes[out_edges[0][1]][coord] for coord in ["x", "y", "z"]])
                out_norm = out_vec / np.linalg.norm(out_vec)
                dot = np.clip(in_norm @ out_norm, -1, 1)
                G[in_node][n]["bifurcation_angle"] = np.pi - np.arccos(dot)  # phi_j from Mynard
        else:
            continue
    return


def assign_radii_files(G, fields):
    if fields:
        radii = fields.get("radius")
        if radii:
            for u, v, data in G.edges(data=True):
                edge_id = data["edge_id"]
                if edge_id < len(radii):  # check that ed
                    radius = radii.get(edge_id)
                    data["radius"] = radius

            return G
        else:
            return G

        res = fields.get("resistance")
        if res:
            res = fields.get(edge_id, 0)
        else:
            res = 0
        # .get returns None by default if not found

    else:
        return G


def assign_anastomosis(G, edge_number):
    for u, v, data in G.edges(data=True):
        if data['edge_id'] == edge_number:
            data['strahler'] = 0
            data['vessel_type'] = 'anastomosis'

    return G


def build_edge_inlet_map(G, inlets):
    """
    Returns {(u, v): inlet_node} for every edge in G.
    """
    edge_to_inlet = {}

    for inlet in inlets:
        queue = deque([inlet])
        visited = {inlet}
        while queue:
            node = queue.popleft()
            for successor in G.successors(node):
                edge_to_inlet[(node, successor)] = inlet  # key on edge
                if successor not in visited:
                    visited.add(successor)
                    queue.append(successor)

    return edge_to_inlet