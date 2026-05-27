import polyscope as ps
import numpy as np
from collections import defaultdict
import pyvista as pv
import networkx as nx
def visualise_tree(G,show_flow =True, region = "all"):
    ps.init()
    if region == "all" or region == "full" or region == "full_tree":
        nodes_array = np.array([
            [G.nodes[n]['x'], G.nodes[n]['y'], G.nodes[n]['z']]
            for n in G.nodes
        ])

        # --- Edges as ndarray ---
        # Each edge is (source, target)
        edges_array = np.array(G.edges)
        radius_array = np.array([G.edges[e]['radius'] for e in G.edges])
        flow_array = np.array([G[u][v]["flow"] for u, v in G.edges()])
        strahler_array = np.array([G[u][v]["strahler"] for u, v in G.edges()])

    elif region == "arteries" or region == "artery":
        radius = []
        flow = []
        edge_rows = []
        artery_nodes = set()
        strahler = []

        for u, v, data in G.edges(data=True):
            if data.get("vessel_type") == "artery" or data.get("vessel_type") == "anastomosis":
                artery_nodes.add(u)
                artery_nodes.add(v)
                edge_rows.append((u, v))
                radius.append(data['radius'])
                flow.append(data['flow'])
                strahler.append(data['strahler'])
        nodes_array = np.array([
            [data["x"], data["y"], data["z"]]
            for node, data in G.nodes(data=True)
            if node in artery_nodes
        ])
        radius_array = np.array(radius)
        edges_array = np.array(edge_rows,dtype=np.int64)
        flow_array = np.array(flow)
        strahler_array = np.array(strahler)

    elif region == "veins" or region == "vein":
        radius = []
        flow = []
        edge_rows = []
        vein_nodes = set()
        strahler = []
        for u, v, data in G.edges(data=True):
            if data.get("vessel_type") == "vein":
                vein_nodes.add(u)
                vein_nodes.add(v)
                edge_rows.append((u, v))
                radius.append(data['radius'])
                flow.append(data['flow'])
                strahler.append(data['strahler'])
        nodes_array = np.array([
            [data["x"], data["y"], data["z"]]
            for node, data in G.nodes(data=True)
            if node in vein_nodes
        ])
        radius_array = np.array(radius)
        edges_array = np.array(edge_rows,dtype=np.int64)
        flow_array = np.array(flow)
        strahler_array = array(strahler)
    joint_radii, _ = compute_joint_radii(nodes_array, edges_array, radius_array)

    print('Visualizing now!')

    # Register tree
    tree = ps.register_curve_network("tree", nodes_array, edges_array, color=[155 / 255, 155 / 255, 155 / 255])
    tree.add_scalar_quantity("Node Radius", joint_radii, defined_on='nodes',
                             enabled=True)  # this is sufficient to visualise varied edge radii
    tree.set_node_radius_quantity("Node Radius")  # this is sufficient to visualise varied edge radii
    if show_flow:
        tree.add_scalar_quantity("Edge flow", flow_array, defined_on='edges',
                                 enabled=True)  # this is sufficient to visualise varied edge radii


    # Set up planes
    cor_plane_pos = ps.add_scene_slice_plane()
    cor_plane_pos.set_pose([0, 0, 0], [0, 1, 0])
    cor_plane_pos.set_draw_widget(True)
    cor_plane_pos.set_active(False)

    cor_plane_neg = ps.add_scene_slice_plane()
    cor_plane_neg.set_pose([0, 0, 0], [0, -1, 0])
    cor_plane_neg.set_draw_widget(True)
    cor_plane_neg.set_active(False)

    ax_plane_pos = ps.add_scene_slice_plane()
    ax_plane_pos.set_pose([0, 0, 0], [0, 0, 1])
    ax_plane_pos.set_draw_widget(True)
    ax_plane_pos.set_active(False)

    ax_plane_neg = ps.add_scene_slice_plane()
    ax_plane_neg.set_pose([0, 0, 0], [0, 0, -1])
    ax_plane_neg.set_draw_widget(True)
    ax_plane_neg.set_active(False)

    # Misc settings
    ps.set_ground_plane_mode("none")
    ps.set_navigation_style("free")
    ps.set_up_dir("z_up")
    ps.set_front_dir("neg_y_front")
    ps.set_background_color([0, 0, 0])
    ps.show()
    return

def compute_joint_radii(nodes, edges, edge_mid_radii):
    """
    Given a radius for each edge (at midpoint), compute per-node radius
    as the average of connected edge radii. If a node has only one incident
    edge, use that edge's radius directly.

    Parameters
    ----------
    nodes : (N,3)
    edges : (E,2)
    edge_mid_radii : (E,)

    Returns
    -------
    joint_radii : (N,) per-node radii
    edge_radii : (E,2) per-edge start/end radii
    """
    N = nodes.shape[0]
    E = edges.shape[0]

    # collect radii per node
    incident = defaultdict(list)
    for e, (i0, i1) in enumerate(edges):
        incident[i0].append(edge_mid_radii[e])
        incident[i1].append(edge_mid_radii[e])

    joint_radii = np.zeros(N, dtype=float)
    for i in range(N):
        if len(incident[i]) == 0:
            joint_radii[i] = 0.0
        elif len(incident[i]) == 1:
            # leaf: just use that edge's radius
            joint_radii[i] = incident[i][0]
        else:
            # average of all connected edge radii
            joint_radii[i] = np.mean(incident[i])

    # now expand to per-edge start/end
    edge_radii = np.zeros((E, 2), dtype=float)
    for e, (i0, i1) in enumerate(edges):
        edge_radii[e, 0] = joint_radii[i0]
        edge_radii[e, 1] = joint_radii[i1]

    return joint_radii, edge_radii

def generate_visualisation_arrays(coords, edges):
    """coords - coordinate array for glaboal node numbers
    edges, node connections for global node numbers"""
    nodes = np.unique(edges.flatten())
    nodes.sort()
    new_coords = coords[nodes]
    new_map = {p: c for c, p in enumerate(nodes)}
    remapped_connections = []
    for element in edges:
        temp = []
        for node in element:
            temp.append(new_map[node])
        remapped_connections.append(temp)
    remapped_connections = np.array(remapped_connections)
    padding = np.empty(remapped_connections.shape[0], int) * 2
    padding[:] = 2
    connections_with_padding = np.vstack((padding, remapped_connections.T)).T
    return new_coords, connections_with_padding

def visualise_graph_and_field(graph: nx.graph, coords: np.array, field, field_name='radii', title='', need_remap = True):
    if title == '':
        title = f'{field_name} visualisation'
    if need_remap:
        field_remap = remap_node_field_for_vis(graph, field)
    else:
        field_remap = field
    vis_coords, vis_connections = generate_visualisation_arrays(coords, np.array(graph.edges()))
    plotter = pv.Plotter()
    plotter.add_title(title)
    pod = pv.PolyData(vis_coords, lines=vis_connections)
    pod[field_name] = field_remap
    pod_tube = pod.tube(scalars=field_name, radius=1.0, radius_factor = 1)
    plotter.add_mesh(pod_tube, render_lines_as_tubes=True, show_scalar_bar=False)
    plotter.add_scalar_bar(field_name, position_x=0.25)
    plotter.camera_position = 'xz'
    plotter.show()

def remap_node_field_for_vis(graph, field):
    """
    graph - nx.Graph
    field - np.array
    assuming field corresponds to graph.nodes ordering
    """
    graph_node_mapping = {k:v for k, v in zip(graph.nodes, range(graph.number_of_nodes()))}
    nodes = np.unique(np.array(graph.edges()).flatten())
    new_field = []
    for node in nodes:
        new_field.append(field[graph_node_mapping[node]])
    return new_field

def get_coordinate_array(G):
    nodes_array = np.array([
        (data["x"], data["y"], data["z"])
        for _, data in G.nodes(data=True)
    ])
    nodes_array = remap_node_field_for_vis(G,nodes_array)
    nodes_array = np.array(nodes_array)
    return nodes_array
def get_coordinate_field(G,fieldname):
    nodefield_array = np.array([
        (data[fieldname])
        for _, data in G.nodes(data=True)
    ])
    nodefield_array = remap_node_field_for_vis(G,nodefield_array)
    
    return nodefield_array
def get_edge_array(G):
    edge_array = np.array([
        (data["edge_id"], u, v)
        for u, v, data in G.edges(data=True)
    ], dtype=np.int64)
    return edge_array

def visualise_pyvista(G):
    nodes_array = get_coordinate_array(G)
    strahler_array = get_coordinate_field(G,"strahler_order")
    visualise_graph_and_field(G,nodes_array,strahler_array,field_name='strahler',need_remap=True)

def get_tree_properties(G):
    # --- Identify inlet and outlet nodes ---
    if G.is_directed():
        inlet_nodes = [n for n in G.nodes() if G.in_degree(n) == 0]
        outlet_nodes = [n for n in G.nodes() if G.out_degree(n) == 0]
    else:
        # For undirected: degree == 1 typically means a terminal node
        inlet_nodes = [n for n in G.nodes() if G.degree(n) == 1 and 'inlet' in G.nodes[n]]
        outlet_nodes = [n for n in G.nodes() if G.degree(n) == 1 and 'outlet' in G.nodes[n]]

    # --- Get inlet conditions ---
    inlet_conditions = {}
    for node in inlet_nodes:
        # Get the first edge connected to the inlet node
        if G.is_directed():
            edges = list(G.out_edges(node, data=True))
        else:
            edges = list(G.edges(node, data=True))

        if edges:
            u, v, attrs = edges[0]
            inlet_conditions[node] = {
                'pressure': G.nodes[node].get('pressure', None),
                'flow': attrs.get('flow', None),
                'strahler': attrs.get('strahler', None),
                'edge_id': attrs.get('edge_id',None),
                'radius': attrs.get('radius', None),
                'edge': (u,v)
            }
            print(
                f"Inlet node, {node}, has pressure {G.nodes[node].get('pressure', None)} and Inlet edge, {attrs.get('edge_id', None)}, has flow {attrs.get('flow', None)}")
    # --- Get outlet conditions ---
    outlet_conditions = {}
    for node in outlet_nodes:
        if G.is_directed():
            edges = list(G.in_edges(node, data=True))
        else:
            edges = list(G.edges(node, data=True))

        if edges:
            u, v, attrs = edges[0]
            outlet_conditions[node] = {
                'pressure': G.nodes[node].get('pressure', None), 
                'flow': attrs.get('flow', None),
                'edge_id': attrs.get('edge_id', None),
                'radius': attrs.get('radius',None),
                'edge': (u, v)
            }
            print(f"Outlet node, {node}, has pressure {G.nodes[node].get('pressure',None)} and Outlet edge, {attrs.get('edge_id', None)}, has flow {attrs.get('flow', None)}")
    return inlet_conditions, outlet_conditions