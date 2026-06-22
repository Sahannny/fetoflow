import warnings
import numpy as np
import csv
def read_nodes(filename):
    if not isinstance(filename,str):
        raise(TypeError("File name should be a string!"))
    if not filename[-7:] == ".ipnode":
        raise(TypeError("This function expects a .ipnode file."))
    with open(filename, "r") as f:
        lines = f.readlines()
        nodes = {}
        numVars = int(lines[4].split()[-1])
        multi_versions = lines[5].split()[-1]
        if multi_versions == "Y":
            baseStep = numVars + 2
            i = 4 + 2 * numVars + 2 # start row for "Y"
        elif multi_versions == "N":
            nversions = 0
            c_step = 1
            i_step = numVars + 2
            i = 4 + 2 * numVars + 1 # start row for "N"
        else:
            raise TypeError("Double Check File Format")
        while i < len(lines):
            node_id = int(lines[i].split()[-1])
            if multi_versions == "Y":
                nversions = int(lines[i + 1].split()[-1])
                c_step = numVars * nversions - 1
                i_step = baseStep + nversions * numVars + (nversions * numVars if nversions > 1 else 0)
            coords = []
            for c in range(1 + nversions, i_step - 1, c_step):
                coord = lines[i + c].split()[-1]
                coords.append(float(coord))
            nodes[node_id - 1] = coords  # 0-based indexing for the networkX geometry.
            i += i_step
    return nodes


def read_elements(filename):
    if not isinstance(filename,str):
        raise(TypeError("File name should be a string!"))
    if not filename[-7:] == ".ipelem":
        raise(TypeError("This function expects a .ipelem file."))
    
    with open(filename, "r") as f:
        lines = f.readlines()
        elems = []
        i = 5
        while i < len(lines):
            intraElementStep = 5
            nodes = tuple(int(x) - 1 for x in lines[i + intraElementStep].split()[-2:])  # Translate to 0-based indexing.
            elems.append(nodes)
            while len(lines[i].split()) != 0:
                i += 1
                if i + 1 == len(lines):
                    return elems
            i += 1
        return elems


def define_fields_from_files(files: dict[str]):
    """
    Defines field(s) as specified in ipfield file(s).
    """
    if not isinstance(files, dict):
        raise (TypeError("files must be a dictionary in the format files[field_name] = filename"))
    fields = {}
    for field in files.keys():
        file_name = files[field]
        if not file_name[-7:] == ".ipfiel":
            ext_start = -(str.__reversed__(file_name).find(".") + 1)
            if ext_start is not None:
                raise(TypeError(f"This function expects a .ipfiel file, got {file_name[ext_start:]}"))
            else:
                raise(TypeError(f"This function expects a .ipfiel file. No file extension found."))
        with open(file_name, "r") as f:
            warnings.warn(f"Assuming the radii in {file_name} are in mm!!!")
            i = 7  # ignore metadata
            lines = f.readlines()
            max_digits = len(lines[3][lines[3].find(":") + 1 :].strip())
            currentField = {}
            while i < len(lines):
                j = i + 2
                id = (
                    int(lines[i].split(':')[-1].strip()) - 1
                )  # assuming for now these correspond to element ids, -1 to 0 based
                val = float(lines[j][lines[j].find(":")+1:].strip())
                currentField[id] = val
                i += 4
        fields[field] = currentField
    return fields

def read_nodes_exnode(filepath): #TESTED
    nodes = {}

    with open(filepath, "r") as f:
        lines = f.readlines()

    i = 0
    while i < len(lines):
        line = lines[i].strip()

        if line.startswith("Node:"):
            # Extract original node number
            original_id = int(line.split()[1])

            # Read next three lines as coordinates
            x = float(lines[i + 1].strip())
            y = float(lines[i + 2].strip())
            z = float(lines[i + 3].strip())

            # Store zero-indexed
            nodes[original_id - 1] = [x, y, z]

            i += 4  # Move past this node block
        else:
            i += 1

    return nodes

def read_edges_exelem(filepath):
    # First pass: find max element ID so we can size list correctly
    max_elem_id = 0
    with open(filepath, "r") as f:
        for line in f:
            if line.strip().startswith("Element:"):
                elem_id = int(line.split()[1])
                max_elem_id = max(max_elem_id, elem_id)

    # Preallocate list
    elems = [None] * max_elem_id

    # Second pass: fill list
    with open(filepath, "r") as f:
        lines = f.readlines()

    i = 0
    while i < len(lines):
        line = lines[i].strip()

        if line.startswith("Element:"):
            elem_id = int(line.split()[1]) - 1  # zero index

            # Find "Nodes:"
            while not lines[i].strip().startswith("Nodes:"):
                i += 1

            # Next line has the two node IDs
            node_line = lines[i + 1].split()
            node1 = int(node_line[0]) - 1
            node2 = int(node_line[1]) - 1

            elems[elem_id] = (node1, node2)

            i += 2
        else:
            i += 1

    return elems

def set_nodes_from_array(node_array):
    """
    TODO: TEST
    """
    nodes = {int(row[0]): row[1:].tolist() for row in node_array}
    return nodes
def set_edges_from_array(edge_array):
    """
    TODO: TEST
    """
    edge_array = np.asarray(edge_array)

    if edge_array.shape[1] != 3:
        raise ValueError("edge_array must have shape (n,3)")


    max_elem_id = np.max(edge_array[:,0])

    # Convert to zero indexing
    edge_ids = edge_array[:, 0].astype(int) 
    node1 = edge_array[:, 1].astype(int) 
    node2 = edge_array[:, 2].astype(int) 

    # Preallocate list
    elems = [None] * len(edge_array)

    for eid, u, v in zip(edge_ids, node1, node2):
        elems[eid] = (u, v)

    return elems

def define_fields_from_exelem(filepath, fieldname):
    """
    TODO: TEST
    """
    field_values = {}
    with open(filepath, "r") as f:
        for line in f:
            if line.strip().startswith("Element:"):
                elem_id = int(line.split()[1]) - 1

                # Skip until Values:
                for line in f:
                    if line.strip().startswith("Values:"):
                        break

                # Next line contains values
                value_line = next(f).split()
                warnings.warn("Assuming values in mm, converting to m")
                field_values[elem_id] = float(value_line[0])

    return {fieldname: field_values}

def set_fields_from_array(field_array, fieldname):
    field_dict = {fieldname: {i: float(val) for i, val in enumerate(field_array)}}
    return field_dict
def export_exnode(G, groupname, filename, vessel_type = "full"):

    """
    Exports coordinates to exnode or exdata format
    data = array of data
    groupname = what you want your data to be called in cmgui
    filename = file name without extension
    type = exnode or exdata
    Modified from VirtualPregnancy/placentagen
    TODO: TEST MORE
    """
    if vessel_type == "full" or vessel_type == "all":
        nodes_array = np.array([
            (data["x"], data["y"], data["z"])
            for _, data in G.nodes(data=True)
        ])
        nodes_array = remap_node_field_for_vis(G, np.asarray(nodes_array))

    elif vessel_type == "arteries" or vessel_type == "artery":
        artery_nodes = set()

        for u, v, data in G.edges(data=True):
            if data.get("vessel_type") == "artery":
                artery_nodes.add(u)
                artery_nodes.add(v)
        nodes_array = [
            [data["x"], data["y"], data["z"]]
            for node, data in G.nodes(data=True)
            if node in artery_nodes
        ]

    elif vessel_type == "veins" or vessel_type == "vein":
        vein_nodes = set()
        for u, v, data in G.edges(data=True):
            if data.get("vessel_type") == "vein":
                vein_nodes.add(u)
                vein_nodes.add(v)
        nodes_array = [
            [data["x"], data["y"], data["z"]]
            for node, data in G.nodes(data=True)
            if node in vein_nodes
        ]
    
    #nodes_array = remap_node_field_for_vis(,np.asarray(nodes_array))

    print('filename', filename)
    # first entry
    data_num = len(nodes_array)
    filename = filename + '.' + 'exnode'
    f = open(filename, 'w')
    f.write(" Group name: %s\n" % groupname)
    f.write(" #Fields=1\n")
    f.write(" 1) coordinates, coordinate, rectangular cartesian, #Components=3\n")
    f.write(" x.  Value index=1, #Derivatives=0\n")
    f.write(" y.  Value index=1, #Derivatives=0\n")
    f.write(" z.  Value index=1, #Derivatives=0\n")

    for x in range(0, data_num):
        f.write("Node:  "        "%s\n" % (x + 1))
        f.write("          %s\n" % nodes_array[x][0])
        f.write("          %s\n" % nodes_array[x][1])
        f.write("          %s\n" % nodes_array[x][2])
    f.close()
    print('Writing complete')
    return

def export_exelem(G, groupname, filename, vessel_type = "full"):

    """
    Exports element locations to exelem format
    data = array of data
    groupname = what you want your data to be called in cmgui
    filename = file name without extension
    Modified from VirtualPregnancy/placentagen
    TODO: ADD export for arteries and veins

     """
    edge_rows = []
    # Edges array with rows: [edge_id, u, v] from DiGraph
    if (vessel_type == "full" or vessel_type == "all"):
        edge_array = np.array([
            (data["edge_id"], u, v)
            for u, v, data in G.edges(data=True)
        ], dtype=np.int64)
    elif (vessel_type == "arteries" or vessel_type == "artery"):
        for u, v, data in G.edges(data=True):
            if data["vessel_type"]== "artery":
                edge_rows.append((data["edge_id"], u, v))
        edge_array = np.array(edge_rows,dtype=np.int64)
    elif (vessel_type == "veins" or vessel_type == "vein"):
        for u, v, data in G.edges(data=True):
            if data["vessel_type"] == "vein":
                edge_rows.append((data["edge_id"], u, v))
        edge_array = np.array(edge_rows, dtype=np.int64)

    data_num = len(edge_array)
    filename = filename + '.exelem'
    f = open(filename, 'w')
    f.write(" Group name: %s\n" % groupname)
    f.write(" Shape.  Dimension=1\n")
    f.write(" #Scale factor sets= 1\n")
    f.write("   l.Lagrange, #Scale factors= 2\n")
    f.write(" #Nodes=           2\n")
    f.write(" #Fields=1\n")
    f.write(" 1) coordinates, coordinate, rectangular cartesian, #Components=3\n")
    f.write("   x.  l.Lagrange, no modify, standard node based.\n")
    f.write("     #Nodes= 2\n")
    f.write("      1.  #Values=1\n")
    f.write("       Value indices:     1\n")
    f.write("       Scale factor indices:   1\n")
    f.write("      2.  #Values=1\n")
    f.write("       Value indices:     1\n")
    f.write("       Scale factor indices:   2\n")
    f.write("   y.  l.Lagrange, no modify, standard node based.\n")
    f.write("     #Nodes= 2\n")
    f.write("      1.  #Values=1\n")
    f.write("       Value indices:     1\n")
    f.write("       Scale factor indices:   1\n")
    f.write("      2.  #Values=1\n")
    f.write("       Value indices:     1\n")
    f.write("       Scale factor indices:   2\n")
    f.write("   z.  l.Lagrange, no modify, standard node based.\n")
    f.write("     #Nodes= 2\n")
    f.write("      1.  #Values=1\n")
    f.write("       Value indices:     1\n")
    f.write("       Scale factor indices:   1\n")
    f.write("      2.  #Values=1\n")
    f.write("       Value indices:     1\n")
    f.write("       Scale factor indices:   2\n")
    for x in range(0, data_num):
        f.write(" Element:            %s 0 0\n" % int(edge_array[x][0] + 1))
        f.write("   Nodes:\n")
        f.write("                %s            %s\n" % (int(edge_array[x][1] + 1), int(edge_array[x][2] + 1)))
        f.write("   Scale factors:\n")
        f.write("       0.1000000000000000E+01   0.1000000000000000E+01\n")
    f.close()
    return

def export_field(G, groupname, fieldname, filename,vessel_type = 'full'):
    # Exports element locations to exelem format
    # data = array of data
    # groupname = what you want your data to be called in cmgui
    # filename = file name without extension
    values = []
    values_with_idx = []
    if (vessel_type == "full" or vessel_type == "all"):
        edge_id_raw = []
        field = []
        edges_with_field = []
        for u, v, data in G.edges(data=True):
            edge_id_raw.append(data["edge_id"])
            field.append(data[fieldname])
        field_remapped = np.zeros(len(edge_id_raw))
        for i in range(0, len(edge_id_raw)):
            field_remapped[edge_id_raw[i]] = field[i]

    elif (vessel_type == 'arteries' or vessel_type == 'artery'):
        edge_id_raw = []
        field = []
        edges_with_field = []
        for u,v,data in G.edges(data=True):
            if data["vessel_type"] == "artery":
                edge_id_raw.append(data["edge_id"]) 
                field.append(data[fieldname])
        field_remapped = np.zeros(len(edge_id_raw))
        for i in range(0,len(edge_id_raw)):
            field_remapped[edge_id_raw[i]] = field[i]

    elif (vessel_type == 'vein' or vessel_type == 'veins'):
        for e in G.edges:
            if G.edges[e]["vessel_type"] == "vein":
                values.append(G.edges[e][fieldname])
        data = np.array(values)
    data_num = len(field_remapped)
    filename = filename+ '_' + fieldname+  '.exelem'
    f = open(filename, 'w')
    f.write(" Group name: %s\n" % groupname)
    f.write(" Shape.  Dimension=1\n")
    f.write(" #Scale factor sets= 0\n")
    f.write(" #Nodes=           0\n")
    f.write(" #Fields=1\n")
    f.write(" 1) %s, field, rectangular cartesian, #Components=1\n" % fieldname)
    f.write("   %s.  l.Lagrange, no modify, grid based.\n" % fieldname)
    f.write("   #xi1=1 \n")
    for x in range(0, data_num):
        f.write(" Element:            %s 0 0\n" % int(x + 1))
        f.write("   Values:\n")
        f.write(
            "           %s       %s\n" % (
                field_remapped[x], field_remapped[x]))
    f.close()
    return

def export_all(G,groupname,sample,vessel_type = "full"):
    """
    Exports element locations to exelem format
    data = array of data
    groupname = what you want your data to be called in cmgui
    filename = file name without extension
    Modified from VirtualPregnancy/placentagen
     """
    if vessel_type == "full" or vessel_type == "all":
        sample = sample + "_full_tree"
    elif vessel_type == "arteries" or vessel_type == "artery":
        sample = sample + "_artery"
    elif vessel_type == "veins" or vessel_type == "vein":
        sample = sample + "_vein"
    export_exnode(G,groupname,sample,vessel_type)
    export_exelem(G,groupname,sample,vessel_type)
    export_field(G, groupname,'radius',sample,vessel_type)
    export_field(G,groupname,'flow',sample,vessel_type)

def remap_node_field_for_vis(graph, field):
    """
    graph - nx.Graph
    field - np.array
    assuming field corresponds to graph.nodes ordering
    """
    graph_node_mapping = {k: v for k, v in zip(graph.nodes, range(graph.number_of_nodes()))}
    nodes = np.unique(np.array(graph.edges()).flatten())
    new_field = []
    for node in nodes:
        new_field.append(field[graph_node_mapping[node]])
    return new_field

def export_region_as_csv(G, ROI,filename, order_interest = None, chorion_elems = None,  stem_elems = None):
    elements = []
    viscosity = 0.0033600

    if ROI == 'stem_villi':
        if stem_elems == None:
            raise ValueError("No stem villi edge numbers input, please input stem villi edges")

        for u,v,data in G.edges(data=True):
            if data["edge_id"] in stem_elems and data['vessel_type'] == 'artery':
                edge_id = data['edge_id']
                radius = data['radius']
                flow = data['flow']
                pressure_1 = G.nodes[u].get('pressure', None)
                pressure_2 = G.nodes[v].get('pressure', None)
                shear_stress = (4 * viscosity * flow) / (np.pi * (radius ** 3))
                order = data['strahler']
                new_value = np.asarray([edge_id,radius,flow,pressure_1,pressure_2,shear_stress,order])
                elements.append(new_value)
    elif ROI == 'order':

        for u,v,data in G.edges(data=True):
            if data["strahler"] in order_interest and data['vessel_type'] == 'artery':
                edge_id = data['edge_id']
                radius = data['radius']
                flow = data['flow']
                pressure_1 = G.nodes[u].get('pressure', None)
                pressure_2 = G.nodes[v].get('pressure', None)
                shear_stress = (4 * viscosity * flow) / (np.pi * (radius ** 3))
                order = data['strahler']
                new_value = np.asarray([edge_id,radius,flow,pressure_1,pressure_2,shear_stress,order])
                elements.append(new_value)
    elif ROI == 'chorion':
        if len(chorion_elems) == None:
            raise ValueError("No chorion edge numbers input, please input stem villi edges")

        for u,v,data in G.edges(data=True):
            if data["edge_id"] in chorion_elems and data['vessel_type'] == 'artery':
                edge_id = data['edge_id']
                radius = data['radius']
                flow = data['flow']
                pressure_1 = G.nodes[u].get('pressure', None)
                pressure_2 = G.nodes[v].get('pressure', None)
                shear_stress = (4 * viscosity * flow) / (np.pi * (radius ** 3))
                order = data['strahler']
                new_value = np.asarray([edge_id,radius,flow,pressure_1,pressure_2,shear_stress,order])
                elements.append(new_value)
    else:
        print("Unrecognized region")
        return
    element_array = np.array(elements)
    # Calculate averages for each column (excluding the first column)
    average_radius = np.mean(element_array[:, 1])
    average_flow_rate = np.mean(element_array[:, 2])
    average_inlet_pressure = np.mean(element_array[:, 3])
    average_outlet_pressure = np.mean(element_array[:, 4])
    median_radius = np.median(element_array[:, 1])
    median_flow_rate = np.median(element_array[:, 2])
    median_inlet_pressure = np.median(element_array[:, 3])
    median_outlet_pressure = np.median(element_array[:, 4])
    upper_radius = np.max(element_array[:, 1])
    upper_flow_rate = np.max(element_array[:, 2])
    upper_inlet_pressure = np.max(element_array[:, 3])
    upper_outlet_pressure = np.max(element_array[:, 4])
    lower_radius = np.min(element_array[:, 1])
    lower_flow_rate = np.min(element_array[:, 2])
    lower_inlet_pressure = np.min(element_array[:, 3])
    lower_outlet_pressure = np.min(element_array[:, 4])
    average_shear_stress = np.mean(element_array[:, 5])
    median_shear_stress = np.median(element_array[:, 5])
    lower_shear_stress = np.min(element_array[:, 5])
    upper_shear_stress = np.max(element_array[:, 5])
    print('Average Radius:', average_radius)
    print('Average Flow:', average_flow_rate)
    print('Average Inlet Pressure:', average_inlet_pressure)
    print('Average Outlet Pressure:', average_outlet_pressure)
    print('Average Shear Stress: ', average_shear_stress)
    # Prepare the average data row

    average_row = ['Average', average_radius, average_flow_rate, average_inlet_pressure, average_outlet_pressure,
                   average_shear_stress]
    median_row = ['Median', median_radius, median_flow_rate, median_inlet_pressure, median_outlet_pressure,
                  median_shear_stress]
    upper_row = ['Maximum', upper_radius, upper_flow_rate, upper_inlet_pressure, upper_outlet_pressure,
                 upper_shear_stress]
    lower_row = ['Minimum', lower_radius, lower_flow_rate, lower_inlet_pressure, lower_inlet_pressure,
                 lower_shear_stress]

    # Create file if it doesn't exist and write data with headers
    try:
        with open(filename, mode='x', newline='') as file:
            writer = csv.writer(file)
            # Write headers
            writer.writerow(
                ['Element Number', 'Radius', 'Flow Rate', 'Inlet Pressure', 'Outlet Pressure', 'Shear Stress', 'Order'])
            # Write data
            writer.writerows(element_array)

            writer.writerow(average_row)
            writer.writerow(median_row)
            writer.writerow(upper_row)
            writer.writerow(lower_row)

        print("File created and data written with headers.")
    except FileExistsError:
        with open(filename, mode='a', newline='') as file:
            writer = csv.writer(file)
            # Write data without headers
            writer.writerows(element_array)
            writer.writerow(average_row)
            writer.writerow(median_row)
            writer.writerow(upper_row)
            writer.writerow(lower_row)
        print("File already exists. Appending data without headers.")

    return
