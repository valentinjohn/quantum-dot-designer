# -*- coding: utf-8 -*-
"""
Created on Wed Sep 13 09:55:38 2023

@author: vjohn
"""

# %% imports

import math
import numpy as np
from copy import deepcopy

# %% definitions


def rot_mat(theta):
    """
    Compute the rotation matrix for a given angle.

    Args:
        theta: Angle of rotation in radians.

    Returns:
        Rotation matrix as a numpy array.

    """
    return np.array([[np.cos(theta), -np.sin(theta)],
                     [np.sin(theta), np.cos(theta)]])


def gen_poly(n, sp=None):
    """
    Generate a regular polygon with a given number of sides.

    Args:
        n: Number of sides of the polygon.
        sp: Starting point of the polygon. Default is None.

    Returns:
        List of tuples representing the vertices of the polygon.

    """
    mat = rot_mat(2*np.pi/n)
    if sp is None:
        sp = [np.cos(np.pi/n), np.sin(np.pi/n)]
    poly = [sp]
    for k in range(n):
        poly.append(np.dot(mat, poly[-1]))
    return [tuple(co) for co in poly]


def rotate_point(point, angle, origin=(0, 0)):
    """
    Rotate a point counterclockwise around a specified origin.

    Given a point, an angle (in radians), and an optional origin of rotation, this function calculates the new coordinates of the point after it has been rotated counterclockwise by the given angle around the origin.

    Parameters:
    - point (tuple of float): The (x, y) coordinates of the point to be rotated.
    - angle (float): The angle of rotation in radians. Positive values result in counterclockwise rotation.
    - origin (tuple of float, optional): The (x, y) coordinates of the origin of rotation. Defaults to (0, 0).

    Returns:
    - tuple of float: The (x, y) coordinates of the rotated point.

    Example:
    Given point (1, 0), rotation angle of pi/2 (90 degrees) around origin (0, 0):
    Output: (0, 1)
    """
    x, y = point
    ox, oy = origin

    qx = ox + math.cos(angle) * (x - ox) - math.sin(angle) * (y - oy)
    qy = oy + math.sin(angle) * (x - ox) + math.cos(angle) * (y - oy)

    return qx, qy


def generate_clavier_gates(width, length, gate_width, gate_length,
                           num_clavier, spacing, shift, position, rotation=0):
    """
    Generate the vertices for a rectangle with notches (claviers) on its bottom edge

    Parameters:
    - width (float): The width of the main rectangle.
    - length (float): The length of the main rectangle.
    - gate_width (float): The width of each clavier.
    - gate_length (float): The length (or depth) of each clavier.
    - num_clavier (int): The number of claviers.
    - spacing (float): The space between each clavier.
    - shift (float): The horizontal shift applied to the claviers (used to adjust their position).
    - position (tuple): A tuple of (x,y) coordinates representing the position offset for the entire shape.
    - rotation (float, optional): The angle of rotation applied to the entire shape (default is 0).

    Returns:
    - list: A list of (x, y) tuples representing the vertices of the shape.

    Note:
    The `rotate_point` function is not provided in the original code. This function would be responsible for rotating a given point around the origin by a specified rotation angle.
    """
    vertices = []

    # Convert rotation to radians
    rotation = math.radians(rotation)

    # Calculate the total space occupied by the claviers and the gaps
    total_clavier_space = (num_clavier * gate_width +
                           (num_clavier - 1) * spacing)

    # Calculate the start point of the claviers (centered and shifted)
    start_clavier = (length - total_clavier_space) / 2 + shift

    # Start point
    vertices.append(rotate_point((0 + position[0], 0 + position[1]), rotation))

    # Bottom edge of the rectangle
    for i in range(num_clavier + 1):
        # If there is a clavier here
        if i < num_clavier:
            # Move to the clavier
            point = (start_clavier + i * (gate_width + spacing) +
                     position[0], 0 + position[1])
            vertices.append(rotate_point(point, rotation))

            # Traverse the clavier
            point = (start_clavier + i * (gate_width + spacing) +
                     position[0], -gate_length + position[1])
            vertices.append(rotate_point(point, rotation))

            point = (start_clavier + i * (gate_width + spacing) +
                     gate_width + position[0],
                     - gate_length + position[1])
            vertices.append(rotate_point(point, rotation))

            point = (start_clavier + i * (gate_width + spacing) +
                     gate_width + position[0], 0 + position[1])
            vertices.append(rotate_point(point, rotation))

    # Right end of the rectangle
    point = (length + position[0], 0 + position[1])
    vertices.append(rotate_point(point, rotation))

    # Top edge of the rectangle
    point = (length + position[0], width + position[1])
    vertices.append(rotate_point(point, rotation))

    point = (0 + position[0], width + position[1])
    vertices.append(rotate_point(point, rotation))

    # Closing the polygon
    point = (0 + position[0], 0 + position[1])
    vertices.append(rotate_point(point, rotation))

    return vertices


def create_lattice_positions(rows, columns, spacing_x, spacing_y, center):
    """
    Generate a lattice of positions in a 2D grid.

    This function creates a lattice of positions in a 2D grid, taking into account the specified number of rows and columns, the spacing between the points in both  x and y directions, and a specified center point. The lattice is centered on the provided center point.

    Parameters:
    - rows (int): The number of rows in the lattice.
    - columns (int): The number of columns in the lattice.
    - spacing_x (float): The horizontal spacing between consecutive points.
    - spacing_y (float): The vertical spacing between consecutive points.
    - center (tuple of float): The (x, y) coordinates of the center of the lattice.

    Returns:
    - list of list: A list of [x, y] coordinates representing the positions in the lattice.

    Example:
    Given rows=2, columns=2, spacing_x=1, spacing_y=1, and center=(0, 0):
    Output: [[-0.5, -0.5], [0.5, -0.5], [-0.5, 0.5], [0.5, 0.5]]
    """
    x_center, y_center = center
    positions = []

    for i in range(rows):
        for j in range(columns):
            x = (j * spacing_x + x_center) - (columns - 1) * spacing_x / 2
            y = (i * spacing_y + y_center) - (rows - 1) * spacing_y / 2
            positions.append([x, y])
    return positions


def apply_sublattice(main_lattice, sub_lattice):
    """
    Combine a main lattice with a sub-lattice by adding each sub-lattice position to every main lattice position.

    Given a main lattice of positions and a sub-lattice, this function generates a new set of positions by adding each sub-lattice position to every position in the main lattice. This operation effectively combines the two lattices, creating a dense set of points by superimposing the sub-lattice onto each point of the main lattice.

    Parameters:
    - main_lattice (list of list): A list of [x, y] coordinates representing the positions in the main lattice.
    - sub_lattice (list of list): A list of [x, y] coordinates representing the positions in the sub-lattice.

    Returns:
    - list of list: A list of [x, y] coordinates representing the positions resulting from combining the main lattice with the sub-lattice.

    Example:
    Given main_lattice=[[0, 0], [1, 1]] and sub_lattice=[[0, 0], [0.5, 0.5]]:
    Output: [[0, 0], [0.5, 0.5], [1, 1], [1.5, 1.5]]
    """
    result = []

    for main_pos in list(main_lattice):
        for sub_pos in list(sub_lattice):
            result.append([main_pos[0] + sub_pos[0], main_pos[1] + sub_pos[1]])

    return result


def update_positions(elements, rows, columns, spacing_x, spacing_y, center):
    """
    Update the positions of elements based on a specified lattice grid and optionally a sub-lattice.

    Given a dictionary of elements, each with an associated list of positions (potentially representing a sub-lattice), this function recalculates the positions based on a main lattice defined by the number of rows, columns, spacings, and a center point. The sub-lattice positions of each element are then superimposed onto the main lattice.

    Parameters:
    - elements (dict): A dictionary where each key represents an element, and each value is another dictionary containing a 'positions' key. The 'positions' key corresponds to a list of [x, y] coordinates that represent a the positions of the elements before placing them on the lattice grid.
    - rows (int): The number of rows in the main lattice.
    - columns (int): The number of columns in the main lattice.
    - spacing_x (float): The horizontal spacing between consecutive points in the main lattice.
    - spacing_y (float): The vertical spacing between consecutive points in the main lattice.
    - center (tuple of float): The (x, y) coordinates of the center of the main lattice.

    Returns:
    - dict: An updated dictionary where each element's 'positions' list has been recalculated based on the lattice grid.

    Example:
    Given elements = {'A': {'positions': [[0, 0], [0.5, 0.5]]}},
    rows=2, columns=2, spacing_x=1, spacing_y=1, and center=(0, 0):
    Output: {'A': {'positions': [[-0.5, -0.5], [0, 0], [0.5, -0.5], [1, 0],
                                  [-0.5, 0.5], [0, 1], [0.5, 0.5], [1, 1]]}}
    """
    elements_updated = deepcopy(elements)
    # Generate lattice positions
    lattice_positions = create_lattice_positions(
        rows, columns, spacing_x, spacing_y, center)

    for key, value in elements_updated.items():
        # Update the positions of each element in the dictionary
        value['positions'] = apply_sublattice(lattice_positions,
                                              value['positions'])
    elements_updated = sort_positions_in_dict(elements_updated)
    return elements_updated


def sort_positions_in_dict(dictionary):
    """
    Sort the 'positions' list within each dictionary value based on the y-coordinate in descending order, and then by the x-coordinate in ascending order.

    Parameters:
    - dictionary (dict): A dictionary where each value is another dictionary containing a 'positions' key. The 'positions' key corresponds to a list of (x, y) tuples.

    Returns:
    - dict: The input dictionary with the 'positions' lists sorted for each value.

    Example:
    Input:
    {
      'A': {'positions': [(2, 3), (1, 4), (2, 4)]},
      'B': {'positions': [(3, 2), (3, 3)]}
    }
    Output:
    {
      'A': {'positions': [(1, 4), (2, 4), (2, 3)]},
      'B': {'positions': [(3, 3), (3, 2)]}
    }
    """
    for key in dictionary:
        dictionary[key]['positions'].sort(key=lambda x: (-x[1], x[0]))
    return dictionary


def merge_device_positions(dict1, dict2):
    """
    Merge two dictionaries containing semiconductor device information by combining their positions.

    This function merges two dictionaries by combining the 'positions' lists of elements with the same keys. If two elements (from dict1 and dict2) have the same key but different 'vertices' or 'layers', the function raises an error. If the key only exists in one of the dictionaries, the element is added to the merged result.

    Parameters:
    - dict1 (dict): The first dictionary where each key is a device name and each value is a dictionary containing 'vertices', 'layer', and 'positions' keys.
    - dict2 (dict): The second dictionary with a similar structure as dict1.

    Returns:
    - dict: A merged dictionary with combined 'positions' for elements with the same keys.

    Raises:
    - ValueError: If the 'vertices' or 'layers' of an element with the same key in both dictionaries don't match.

    Example:
    Given
    dict1={'A': {'vertices': [[0, 0]],
                  'layer': 1,
                  'positions': [[1, 1]]}}
    and dict2={'A': {'vertices': [[0, 0]],
                      'layer': 1,
                      'positions': [[2, 2]]},
                'B': {'vertices': [[0, 0]],
                      'layer': 1,
                      'positions': [[3, 3]]}},
    Output:
    {'A': {'vertices': [[0, 0]],
            'layer': 1,
            'positions': [[1, 1], [2, 2]]},
      'B': {'vertices': [[0, 0]],
            'layer': 1,
            'positions': [[3, 3]]}}
    """
    # Create a copy of the first dictionary to avoid modifying the original
    merged = dict1.copy()

    for key, value in dict2.items():
        # If the key (name) is already in the merged dictionary
        if key in merged:
            # Check if the vertices are the same
            if not np.array_equal(merged[key]['vertices'], value['vertices']):
                raise ValueError(
                    f"The vertices for the key {key} are different between the two dictionaries.")

            # Check if the layers are the same
            if merged[key]['layer'] != value['layer']:
                raise ValueError(
                    f"The layers for the key {key} are different between the two dictionaries.")

            # Merge the positions
            merged[key]['positions'] += value['positions']

            # Sort the positions based on y_max and x_min
            merged[key]['positions'].sort(key=lambda x: (-x[1], x[0]))

        else:
            # If the key (name) is not in the merged dictionary, just add the whole item
            merged[key] = value

    return merged


def compute_positions(n_lines, spacing):
    """
    Compute the starting positions for a given number of lines, considering their spacing.

    The function calculates the starting positions of equidistant lines. The lines are centered around zero, ensuring that the middle line (or the gap between two middle lines if an even number) aligns with the zero position. The spacing between adjacent lines is constant and defined by the `spacing` parameter.

    Parameters:
    - n_lines (int): The total number of lines to be positioned.
    - spacing (float): The distance between adjacent lines.

    Returns:
    - list of float: A list of starting positions for each line.

    Example:
    Given n_lines=3 and spacing=2:
    Output: [-2.0, 0.0, 2.0]
    """
    start = - (n_lines - 1) * spacing / 2
    return [start + i * spacing for i in range(n_lines)]


def get_polygons_from_path(path_points):
    """
    Generate the vertices of a polygon based on a provided path.

    Given a list of path points, this function calculates the vertices of a polygon that encompasses the path.
    Each path point has three components: (x, y, width). The function uses the direction between consecutive points
    and the specified width to determine the boundary of the path, generating the vertices for the polygon.

    Parameters:
    -----------
    path_points : list of tuple
        A list of tuples where each tuple represents a point in the path. Each tuple should contain three elements:
        x-coordinate, y-coordinate, and the width of the path at that point.

    Returns:
    --------
    list of tuple
        A list of tuples where each tuple represents a vertex of the polygon. The vertices are ordered in such a way
        that they can be directly used for plotting or further polygon-related operations.

    Notes:
    ------
    - The function assumes that the path_points list contains at least two points.
    - The polygon vertices are generated in a way that the left side vertices are followed by the right side vertices
      in reverse order.
    - `perpendicular_vector` is an external function that returns a vector perpendicular to the given vector.

    Examples:
    ---------
    >>> path_points = [(0, 0, 1), (1, 1, 1), (2, 0, 1)]
    >>> get_polygons_from_path(path_points)
    [(0.5, 0.5), (1.5, 1.5), (2.5, -0.5), (1.5, -1.5), (0.5, -0.5)]

    """
    end = path_points[-1]
    vertices_left = []
    vertices_right = []
    for i in range(len(path_points) - 1):
        direction_current = np.array(
            path_points[i+1][:2], dtype=float) - np.array(path_points[i][:2], dtype=float)
        norm_current = np.linalg.norm(direction_current)
        if norm_current != 0:
            direction_current /= norm_current
        if i == 0:
            bisector_direction = direction_current
        else:
            direction_prev = np.array(
                path_points[i][:2], dtype=float) - np.array(path_points[i-1][:2], dtype=float)
            norm_prev = np.linalg.norm(direction_prev)
            if norm_prev != 0:
                direction_prev /= norm_prev
            bisector_direction = direction_current + direction_prev
            norm_bisector = np.linalg.norm(bisector_direction)
            if norm_bisector != 0:
                bisector_direction /= norm_bisector
        normal = perpendicular_vector(bisector_direction)
        width = path_points[i][2]
        v1 = np.array(path_points[i][:2], dtype=float) + width/2 * normal
        v2 = np.array(path_points[i][:2], dtype=float) - width/2 * normal
        vertices_left.append((v1[0], v1[1]))
        vertices_right.append((v2[0], v2[1]))
    direction_end = np.array(end[:2], dtype=float) - \
        np.array(path_points[-2][:2], dtype=float)
    normal_end = perpendicular_vector(direction_end)
    v1 = np.array(end[:2], dtype=float) + end[2]/2 * normal_end
    v2 = np.array(end[:2], dtype=float) - end[2]/2 * normal_end
    vertices_left.append((v1[0], v1[1]))
    vertices_right.append((v2[0], v2[1]))

    poly_vertices = vertices_left + vertices_right[::-1]

    return poly_vertices


def compute_fanout_positions(fo_stages, fanout_counts, spacings):
    """
    Computes the fanout line positions for various stages of rectangles in a semiconductor layout.

    This function calculates the positions of fanout lines for each side (top, bottom, left, right) of three distinct stages of rectangles: device, intermediate, and bondpad. The positions are determined based on the fanout counts and spacings provided.

    Parameters:
    - fo_stages (list of tuple): A list of three tuples, representing the dimensions (length, width) of the device, intermediate, and bondpad rectangles, respectively.
    - fanout_counts (dict): A dictionary specifying the number of fanout lines for each side. Keys are 'top', 'bottom', 'left', and 'right', and values are integers representing fanout counts.
    - spacings (list of float): A list of three floats representing the spacings for the device, intermediate, and bondpad rectangles, respectively.

    Returns:
    - dict: A nested dictionary where the outer keys are the rectangle stages ('device', 'intermediate', 'bondpad'), and the inner keys are the sides ('top', 'bottom', 'left', 'right'). The values are lists of float values, representing the computed positions for fanout lines.

    Example:
    Given
    fo_stages=[(16, 16), (500, 530), (1200, 1200)],
    fanout_counts={'top': 3, 'bottom': 3, 'left': 2, 'right': 2},
    and spacings=[2, 80, 200],
    Output will have computed positions for each side of each rectangle stage based on the spacings and fanout counts.
    """

    device_rect, intermediate_rect, bondpad_square = fo_stages
    device_spacing, intermediate_spacing, bondpad_spacing = spacings

    # Compute positions for each rectangle and side
    fanout_positions = {
        'device': {
            'top': compute_positions(fanout_counts['top'], device_spacing),
            'bottom': compute_positions(fanout_counts['bottom'], device_spacing),
            'left': compute_positions(fanout_counts['left'], device_spacing),
            'right': compute_positions(fanout_counts['right'], device_spacing)
        },
        'intermediate': {
            'top': compute_positions(fanout_counts['top'], intermediate_spacing),
            'bottom': compute_positions(fanout_counts['bottom'], intermediate_spacing),
            'left': compute_positions(fanout_counts['left'], intermediate_spacing),
            'right': compute_positions(fanout_counts['right'], intermediate_spacing)
        },
        'bondpad': {
            'top': compute_positions(fanout_counts['top'], bondpad_spacing),
            'bottom': compute_positions(fanout_counts['bottom'], bondpad_spacing),
            'left': compute_positions(fanout_counts['left'], bondpad_spacing),
            'right': compute_positions(fanout_counts['right'], bondpad_spacing)
        }
    }

    return fanout_positions


def get_fo_lines(fanout_positions, fanout_counts, fo_stages):
    """
    Compute the coordinates for fanout lines for various stages of rectangles in a semiconductor layout.

    This function generates the line segments representing the fanout paths for each side (top, bottom, left, right) of three distinct stages of rectangles: device, intermediate, and bondpad. The fanout paths are determined based on the fanout positions, counts, and the dimensions of each stage.

    Parameters:
    - fanout_positions (dict): A nested dictionary where the outer keys are the rectangle stages ('device', 'intermediate', 'bondpad'), and the inner keys are the sides ('top', 'bottom', 'left', 'right'). The values are lists of float values, representing the computed positions for fanout lines.
    - fanout_counts (dict): A dictionary specifying the number of fanout lines for each side. Keys are 'top', 'bottom', 'left', and 'right', and values are integers representing fanout counts.
    - fo_stages (list of tuple): A list of three tuples, representing the dimensions (length, width) of the device, intermediate, and bondpad rectangles, respectively.

    Returns:
    - dict: A dictionary where the keys are the sides ('top', 'bottom', 'left', 'right') and the values are lists of line segments representing the fanout paths for that side. Each line segment consists of start and end coordinates.
    """
    fo_lines = {'top': [],
                'bottom': [],
                'right': [],
                'left': []}
    for n_fo in range(fanout_counts['top']):
        fo_line = [[fanout_positions['device']['top'][n_fo],
                    fo_stages[0][1]],
                   [fanout_positions['intermediate']['top'][n_fo],
                    fo_stages[1][1]],
                   [fanout_positions['bondpad']['top'][n_fo],
                    fo_stages[2][1]]
                   ]
        fo_lines['top'].append(fo_line)

    for n_fo in range(fanout_counts['bottom']):
        fo_line = [[fanout_positions['device']['bottom'][n_fo],
                    -fo_stages[0][1]],
                   [fanout_positions['intermediate']['bottom'][n_fo],
                    -fo_stages[1][1]],
                   [fanout_positions['bondpad']['bottom'][n_fo],
                    -fo_stages[2][1]]
                   ]
        fo_lines['bottom'].append(fo_line)

    for n_fo in range(fanout_counts['right']):
        fo_line = [[fo_stages[0][0],
                    -fanout_positions['device']['right'][n_fo]],
                   [fo_stages[1][0],
                    -fanout_positions['intermediate']['right'][n_fo]],
                   [fo_stages[2][0],
                    -fanout_positions['bondpad']['right'][n_fo]]
                   ]
        fo_lines['right'].append(fo_line)

    for n_fo in range(fanout_counts['left']):
        fo_line = [[-fo_stages[0][0],
                    -fanout_positions['device']['left'][n_fo]],
                   [-fo_stages[1][0],
                    -fanout_positions['intermediate']['left'][n_fo]],
                   [-fo_stages[2][0],
                    -fanout_positions['bondpad']['left'][n_fo]]
                   ]
        fo_lines['left'].append(fo_line)
    return fo_lines


def calculate_vertex(point, offset):
    """Helper function to add an offset to a given point and return the new point."""
    return list(np.array(point) + np.array(offset))


def generate_polygon_for_fanout(direction, fo_line, bondpad_position, bondpad_size, fo_widths, fo_fine_coarse_overlap):
    """
    Generate a polygon for a single fanout based on the given direction.

    Parameters:
    - direction (str): Indicates the direction of the fanout. Can be one of ['top', 'bottom', 'left', 'right'].
    - fo_line (list): List of 3 tuples representing the fanout lines in format [(x1, y1), (x2, y2), (x3, y3)].
    - bondpad_position (dict): Dictionary containing the position of the bondpad for each direction.
    - bondpad_size (dict): Dictionary containing the size of the bondpad for each direction.
    - fo_widths (list): List of 3 floats representing the widths of the fanout at each of the 3 fanout lines.
    - fo_fine_coarse_overlap (float): The overlap between the fine and coarse parts of the fanout.

    Returns:
    - fo_polys (list): List of tuples representing the vertices of the generated polygon in the format [(x1, y1), (x2, y2), ...].

    Notes:
    The function calculates the vertices of the polygon based on the provided parameters and direction.
    It first adjusts the calculations based on the direction. Then, it calculates the bondpad corners,
    the bondpad to fine vertices, the overlap between fine and coarse, and the fine to bondpad vertices.
    Finally, it returns a list of these vertices in order to generate the desired polygon.
    """

    # Adjustments based on direction
    if direction in ['top', 'bottom']:
        coord_index = 0
        multiplier = 1 if direction == 'top' else -1
    else:  # 'left' or 'right'
        coord_index = 1
        multiplier = 1 if direction == 'right' else -1

    # Base points and offsets
    if direction in ['top', 'bottom']:
        base_point = [fo_line[2][0],
                      multiplier*abs(bondpad_position[direction])]
    else:
        base_point = [multiplier*abs(bondpad_position[direction]),
                      fo_line[2][1]]

    bondpad_offset = [multiplier * bondpad_size[direction]
                      [0]/2, multiplier * bondpad_size[direction][1]/2]
    fo_width_offset = [0, 0]
    fo_width_offset[coord_index] = multiplier * fo_widths[2]/2

    # Calculating bondpad corners
    if direction in ['top', 'bottom']:
        top_left = calculate_vertex(
            base_point, [-bondpad_offset[0], bondpad_offset[1]])
        top_right = calculate_vertex(
            base_point, [bondpad_offset[0], bondpad_offset[1]])
        bottom_left = calculate_vertex(
            base_point, [-bondpad_offset[0], -bondpad_offset[1]])
        bottom_right = calculate_vertex(
            base_point, [bondpad_offset[0], -bondpad_offset[1]])
        center_bottom_left = calculate_vertex(
            base_point, [-fo_widths[2]/2, -bondpad_offset[1]])
        center_bottom_right = calculate_vertex(
            base_point, [fo_widths[2]/2, -bondpad_offset[1]])

        bondpad = [center_bottom_right, bottom_right, top_right,
                   top_left, bottom_left, center_bottom_left]
        dir_unit = np.array([1, 0])
        dir_unit_orth = np.array([0, 1])
    else:
        top_left = calculate_vertex(
            base_point, [-bondpad_offset[0], bondpad_offset[1]])
        top_right = calculate_vertex(
            base_point, [bondpad_offset[0], bondpad_offset[1]])
        bottom_left = calculate_vertex(
            base_point, [-bondpad_offset[0], -bondpad_offset[1]])
        bottom_right = calculate_vertex(
            base_point, [bondpad_offset[0], -bondpad_offset[1]])
        center_bottom_left = calculate_vertex(
            base_point, [-bondpad_offset[0], -fo_widths[2]/2])
        center_top_left = calculate_vertex(
            base_point, [-bondpad_offset[0], fo_widths[2]/2])

        bondpad = [center_top_left, top_left, top_right,
                   bottom_right, bottom_left, center_bottom_left]
        dir_unit = np.array([0, 1])
        dir_unit_orth = np.array([1, 0])

    bondpad_to_fine = [calculate_vertex(fo_line[2], -dir_unit*fo_widths[2]/2),
                       calculate_vertex(fo_line[1], -dir_unit*fo_widths[1]/2),
                       calculate_vertex(fo_line[0], -dir_unit*fo_widths[0]/2)]

    overlap_fine_coarse = [calculate_vertex(fo_line[0],
                                            -dir_unit*fo_widths[0]/2 -
                                            dir_unit_orth*fo_fine_coarse_overlap * multiplier),
                           calculate_vertex(fo_line[0],
                                            dir_unit*fo_widths[0]/2 -
                                            dir_unit_orth*fo_fine_coarse_overlap * multiplier)]

    fine_to_bondpad = [calculate_vertex(fo_line[0], dir_unit*fo_widths[0]/2),
                       calculate_vertex(fo_line[1], dir_unit*fo_widths[1]/2),
                       calculate_vertex(fo_line[2], dir_unit*fo_widths[2]/2)]

    fo_polys = []
    fo_polys.extend(bondpad_to_fine)
    fo_polys.extend(overlap_fine_coarse)
    fo_polys.extend(fine_to_bondpad)
    fo_polys.extend(bondpad)

    return fo_polys


def perpendicular_vector(v):
    """Return a 2D vector that's perpendicular to the given 2D vector `v`."""
    return np.array([-v[1], v[0]])


def axis_value(points, axis=0, operation="max"):
    """
    Returns the maximum or minimum value along a specified axis for a list of 2D points.

    :param points: List of 2D points.
    :param axis: Axis along which to find the value (0 or 1).
    :param operation: Operation to perform ("max" or "min").
    :return: Maximum or minimum value along the specified axis.
    """
    if operation == "max":
        return max(point[axis] for point in points)
    elif operation == "min":
        return min(point[axis] for point in points)
    else:
        raise ValueError("Operation must be 'max' or 'min'")


def normalize_vector(v):
    """Return the normalized version of the vector `v`."""
    magnitude = np.linalg.norm(v)
    if magnitude == 0:
        # Avoid division by zero and return the zero vector
        return v
    return v / magnitude


def perpendicular_vector(v):
    """Return a normalized 2D vector that's perpendicular to the given 2D vector `v`."""
    perp_vector = np.array([-v[1], v[0]])
    return normalize_vector(perp_vector)

# %%% screening gate calculations


def compute_distance(point1, point2):
    """Compute the distance between two points."""
    x1, y1 = point1[0], point1[1]
    x2, y2 = point2[0], point2[1]
    return math.sqrt((x2 - x1)**2 + (y2 - y1)**2)


def interpolate_point(point1, point2, ratio):
    """Interpolate between two points based on a given ratio."""
    x1, y1 = point1[0], point1[1]
    x2, y2 = point2[0], point2[1]

    x_new = x1 + ratio * (x2 - x1)
    y_new = y1 + ratio * (y2 - y1)

    return [x_new, y_new]


def create_segmented_path(path, start_length, end_length, width=None, relative_width=False):
    """Create a new path between the start_length and end_length based on the input path."""

    # Create the new path
    new_path = []
    accumulated_length = 0
    started = False  # To determine if we have started recording the new path

    for i in range(1, len(path)):
        segment_length = compute_distance(path[i-1], path[i])

        # Check if the current accumulated length is within the desired range
        if accumulated_length + segment_length >= start_length and not started:
            # Calculate the ratio to interpolate the start point
            remaining_length = start_length - accumulated_length
            ratio = remaining_length / segment_length
            interpolated_start = interpolate_point(path[i-1], path[i], ratio)
            interpolated_start.append(path[i][2])  # Add the width
            new_path.append(interpolated_start)
            started = True

        if started:
            if accumulated_length + segment_length <= end_length:
                new_path.append(path[i])
            else:
                # Calculate the ratio to interpolate the end point
                remaining_length = end_length - accumulated_length
                ratio = remaining_length / segment_length
                interpolated_end = interpolate_point(path[i-1], path[i], ratio)
                interpolated_end.append(path[i][2])  # Add the width
                new_path.append(interpolated_end)
                break

        accumulated_length += segment_length

    # Adjust the width if necessary
    if width is not None:
        for point in new_path:
            point[2] = width if not relative_width else point[2] * width

    return new_path
