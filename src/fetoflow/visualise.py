import polyscope as ps
import numpy as np
from collections import defaultdict

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

    elif region == "arteries" or region == "artery":
        radius = []
        flow = []
        edge_rows = []
        artery_nodes = set()
        for u, v, data in G.edges(data=True):
            if data.get("vessel_type") == "artery":
                artery_nodes.add(u)
                artery_nodes.add(v)
                edge_rows.append((u, v))
                radius.append(data['radius'])
                flow.append(data['flow'])
        nodes_array = np.array([
            [data["x"], data["y"], data["z"]]
            for node, data in G.nodes(data=True)
            if node in artery_nodes
        ])
        radius_array = np.array(radius)
        edges_array = np.array(edge_rows,dtype=np.int64)
        flow_array = np.array(flow)

    elif region == "veins" or region == "vein":
        radius = []
        flow = []
        edge_rows = []
        vein_nodes = set()
        for u, v, data in G.edges(data=True):
            if data.get("vessel_type") == "vein":
                vein_nodes.add(u)
                vein_nodes.add(v)
                edge_rows.append((u, v))
                radius.append(data['radius'])
                flow.append(data['flow'])
        nodes_array = np.array([
            [data["x"], data["y"], data["z"]]
            for node, data in G.nodes(data=True)
            if node in vein_nodes
        ])
        radius_array = np.array(radius)
        edges_array = np.array(edge_rows,dtype=np.int64)
        flow_array = np.array(flow)

    joint_radii, _ = compute_joint_radii(nodes_array, edges_array, radius_array)

    print('Visualizing now!')

    # Register tree
    tree = ps.register_curve_network("tree", nodes_array, edges_array, color=[155 / 255, 155 / 255, 155 / 255])
    tree.add_scalar_quantity("radius node", joint_radii, defined_on='nodes',
                             enabled=True)  # this is sufficient to visualise varied edge radii
    tree.set_node_radius_quantity("radius node")  # this is sufficient to visualise varied edge radii
    if show_flow:
        tree.add_scalar_quantity("flow elem", flow_array, defined_on='edges',
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