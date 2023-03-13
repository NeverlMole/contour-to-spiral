from shapely import Polygon
from shapely import LineString
from shapely import LinearRing
import shapely
import sys
import math
import numpy as np

SPLIT_ANG = 0
INIT_ANG = 0.1
MAX_RG = 1e5

center_point = np.array([0.0, 0.0])

def angle_to_dir(phi):
    return MAX_RG * np.array([math.cos(phi), math.sin(phi)])

def dir_to_ray(dir, center = np.array([0, 0])):
    return LineString([(center[0], center[1]), (center[0] + dir[0], center[1] + dir[1])])

def dir_to_angle(dir):
    return np.arctan2(dir[1], dir[0])

def get_dis(p1, p2):
    return np.sqrt(np.sum(np.square(p1) + np.square(p2)))

class Path:
    points = []

    def __init__ (self, _points):
        self.points = _points

    def __str__(self):
        return str(len(self.points)) + '\n' +\
            '\n'.join([' '.join([str(point[0]), str(point[1]), str(radius)]) for point, radius in self.points])

    def split(self, split_ray, is_forward=True):
        inter_label = -1
        intersect = None
        for i in range(len(self.points) - 1):
            cur_line = LineString([(self.points[i][0][0], self.points[i][0][1]),
                                   (self.points[i + 1][0][0], self.points[i + 1][0][1])])
            intersect = shapely.intersection(cur_line, split_ray)
            if not intersect.is_empty:
                inter_label = i
                intersect = np.array([intersect.xy[0][0], intersect.xy[1][0]])
                break

        if inter_label == -1: return None

        alpha = get_dis(intersect, self.points[inter_label][0]) /\
                (get_dis(self.points[inter_label][0], self.points[inter_label + 1][0]) + 1e-15)
        inter_radius = self.points[inter_label][1] * (1 - alpha) + self.points[inter_label + 1][1] * alpha

        if is_forward:
            split_points = self.points[inter_label + 1:]
            split_points.insert(0, (intersect, inter_radius))
            return Path(split_points)
        else:
            split_points = self.points[:inter_label]
            split_points.append((intersect, inter_radius))
            return Path(split_points)


    def link(self, path):
        if get_dis(self.points[-1][0], path.points[0][0]) < 1e-9:
            self.points.extend(path.points[1:])
        else:
            self.points.extend(path.points)

    def reverse(self):
        self.points.reverse()


def get_contour(points):
    points = [((point - center_point), radius) for point, radius in points] # Shift to center
    points.append(points[0])
    split_ray = dir_to_ray(angle_to_dir(SPLIT_ANG))

    contour_path = Path(points)
    contour = contour_path.split(split_ray)
    contour.link(contour_path.split(split_ray, False))
    contour.reverse()
    return contour

def find_center(points):
    lb = min([point[0] for point, radius in points])
    rb = max([point[0] for point, radius in points])
    db = min([point[1] for point, radius in points])
    ub = max([point[1] for point, radius in points])

    return np.array([(lb + rb) / 2, (ub + db) / 2])


def get_contours(f):
    cur_layer = 0
    contour_points = []
    while True:
        s = f.readline()
        if s == '' : break
        s = s.split(' ')
        num_points = int(s[1])
        layer_label = int(s[2])

        # more than one contour for some layer
        if layer_label != cur_layer:
            print("Error: more than one contour.")

        contour_points.append([])
        for i in range(num_points):
            s = f.readline().split(' ')
            point = np.array([float(s[0]), float(s[1])])
            radius = (float(s[2]))
            contour_points[cur_layer].append((point, radius))

        cur_layer = cur_layer + 1

    global center_point
    center_point = find_center(contour_points[-1])
    print('The center of the contours is (', center_point[0], ',', center_point[1], ')')
    contours = []
    for i in range(cur_layer-1, -1, -1):
        contours.append(get_contour(contour_points[i]))

    return contours


############################### General-Algorithm ################################################

def compute_spiral(contours, alg_init, alg_next_cut, alg_last_cut):
    # initial path
    init_split_ray = dir_to_ray(angle_to_dir(INIT_ANG))
    spiral = contours[0].split(init_split_ray)

    prev_state = alg_init(spiral)
    for i in range(len(contours) - 1):
        new_ang, cut_ang, prev_state = alg_next_cut(contours[i], contours[i + 1], prev_state)
        # compute spiral
        split_ray_1 = dir_to_ray(angle_to_dir(new_ang))
        spiral.link(contours[i].split(split_ray_1, False))
        split_ray_2 = dir_to_ray(angle_to_dir(cut_ang), spiral.points[-1][0])
        new_path = contours[i + 1].split(split_ray_2)
        spiral.link(new_path)

    last_ang = alg_last_cut(contours[-1], prev_state)
    last_split_ray = dir_to_ray(angle_to_dir(last_ang))
    spiral.link(contours[-1].split(last_split_ray, False))
    return spiral


################################# Right-Angle-Algorithm ##########################################

def get_new_cut_right_angle(prev_state):
    prev_point, radius = prev_state
    dis = get_dis(np.array([0, 0]), prev_point)
    prev_angle = dir_to_angle(prev_point)
    new_angle = 0.0
    if dis < radius:
        new_angle = prev_angle + math.pi * 0.7
    else:
        new_angle = prev_angle + 2 * np.arcsin(radius / dis)

    return new_angle

def right_angle_init(spiral):
    return spiral.points[0]

def right_angle_next_cut(cur_contour, next_contour, prev_state):
    new_angle = get_new_cut_right_angle(prev_state)
    split_ray = dir_to_ray(angle_to_dir(new_angle))
    new_path = next_contour.split(split_ray)
    return new_angle, new_angle, new_path.points[0]

def right_angle_last_cut(last_contour, prev_state):
    return get_new_cut_right_angle(prev_state)

def compute_spiral_right_angle(contours):
    return compute_spiral(contours, right_angle_init, right_angle_next_cut, right_angle_last_cut)


############################################## MAIN ##############################################

if __name__ == '__main__':
    shape_file = sys.argv[1]
    contour_file = sys.argv[2]
    output_file = shape_file + '.spiral'

    contours = []
    with open(contour_file, 'r', encoding='utf-8') as f:
        contours = get_contours(f)

    spiral = compute_spiral_right_angle(contours)

    with open(output_file, 'w', encoding='utf-8') as f:
        f.write(str(spiral))
