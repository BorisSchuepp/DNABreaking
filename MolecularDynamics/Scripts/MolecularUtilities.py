

# Script to handle dihedral angle calculations
# Author: Boris N. Schüpp

import numpy as np


def normalize_vector(vector):
    magnitude = np.linalg.norm(vector)
    return vector / magnitude


def translate_to_origin(c_in, *coords):
    return [np.array(point) - np.array(c_in) for point in coords]


# noinspection PyUnreachableCode
def rotate_to_z_axis(a_in, b_in):
    z_axis = np.array([0, 0, 1])
    rotation_axis = np.cross(b_in, z_axis)

    # Only perform rotation if needed
    rotation_angle = 0
    if np.allclose(b_in, -z_axis):
        rotation_matrix = np.eye(3)
    else:
        rotation_angle = np.arccos(np.dot(b_in, z_axis) / np.linalg.norm(b_in))
        rotation_matrix = rotation_matrix_axis_angle(rotation_axis, rotation_angle)

    b_rotated = np.dot(rotation_matrix, b_in)

    # Check if the rotation was in the right direction, otherwise rotate by another 180 deg
    if b_rotated[2] > 0:
        rotation_angle = rotation_angle + np.pi
        rotation_matrix = rotation_matrix_axis_angle(rotation_axis, rotation_angle)

    x_axis = np.array([1, 0, 0])

    # Rotate A according the first rotation, then determine the second rotation
    a_first_rot = np.dot(rotation_matrix, a_in)
    projection_a = np.array([a_first_rot[0], a_first_rot[1], 0])

    xy_rotation_angle = - np.arccos(np.dot(projection_a, x_axis) / np.linalg.norm(projection_a))
    xy_rotation_matrix = rotation_matrix_axis_angle(z_axis, xy_rotation_angle)
    a_rotated = np.dot(xy_rotation_matrix, a_first_rot)

    if abs(a_rotated[1]) > 1.e-5:
        xy_rotation_angle = - xy_rotation_angle
        xy_rotation_matrix = rotation_matrix_axis_angle(z_axis, xy_rotation_angle)
        a_rotated = np.dot(xy_rotation_matrix, a_first_rot)

    if a_rotated[0] < 0:
        xy_rotation_angle = xy_rotation_angle + np.pi
        xy_rotation_matrix = rotation_matrix_axis_angle(z_axis, xy_rotation_angle)

    # Combine the two rotations
    total_rotation = np.dot(xy_rotation_matrix, rotation_matrix)

    return total_rotation


def rotation_matrix_axis_angle(axis, angle_in):
    axis = axis / np.linalg.norm(axis)
    matrix_base = np.array([[0, -axis[2], axis[1]],
                            [axis[2], 0, -axis[0]],
                            [-axis[1], axis[0], 0]])
    return np.eye(3) + np.sin(angle_in) * matrix_base + (1 - np.cos(angle_in)) * np.dot(matrix_base, matrix_base)


def compute_d_coordinates(distance_in, angle_in, dihedral_in):
    z = -distance_in * np.cos(angle_in)
    r = distance_in * np.sin(angle_in)
    x = r * np.cos(dihedral_in)
    y = r * np.sin(dihedral_in)
    return np.array([x, y, z])


def find_d_coordinates(a_in, b_in, c_in, distance_in, angle_in, dihedral_in):
    angle_in = np.radians(angle_in)
    dihedral_in = np.radians(dihedral_in)

    a_translated, b_translated, c_translated = translate_to_origin(c_in, a_in, b_in, c_in)
    rotation_matrix = rotate_to_z_axis(a_translated, b_translated)

    d_new_coordinates = compute_d_coordinates(distance_in, angle_in, dihedral_in)

    # Reverse the transformations
    d_result = np.dot(np.linalg.inv(rotation_matrix), d_new_coordinates) + c_in
    return d_result


def dihedral_angle_gromacs_style(a_in, b_in, c_in, d_in):
    # Calculate planes
    ij = a_in - b_in
    kj = c_in - b_in
    kl = d_in - c_in

    # Get normal vectors
    r1 = np.cross(ij, kj)
    r2 = np.cross(kl, kj)

    # Calculate angle between normal vectors
    n1 = sum([ai ** 2 for ai in r1]) ** (1 / 2)
    n2 = sum([ai ** 2 for ai in r2]) ** (1 / 2)
    sp = sum([ai * bi for ai, bi in zip(r1, r2)])
    angle_res = np.degrees(np.arccos(sp / (n1 * n2)))
    sign = np.sign(np.dot(-r1, kl))
    angle_res *= sign

    return angle_res


def run_test(runs):
    max_distance = 0
    for i in range(0, runs):
        a = np.random.rand(3)
        b = np.random.rand(3)
        c = np.random.rand(3)
        d = np.random.rand(3)

        cd = d - c
        cb = b - c

        distance = np.linalg.norm(cd)
        angle = np.degrees(np.arccos(sum([k * l for k, l in zip(cd, cb)]) / (np.linalg.norm(cd) * np.linalg.norm(cb))))
        dihedral = dihedral_angle_gromacs_style(a, b, c, d)

        d_new = find_d_coordinates(a, b, c, distance, angle, dihedral)

        max_distance = max(max_distance, np.linalg.norm(d - d_new))

    print(f"Max deviation in {runs} tries was {max_distance}.")
