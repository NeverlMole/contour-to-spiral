from shapely.geometry import Polygon
from shapely.geometry import LineString
from shapely.geometry import LinearRing
import shapely
import sys
import math
import numpy as np
import geopandas as gpd

import matplotlib.pyplot as plt

SPLIT_ANG = 0
INIT_ANG = 0.3
MAX_RG = 1e5
CIR_DIVID_NUM = 30
RIGHT_ANGLE_ANG = 0
ANGLE_SEARCH_WEIGHT_UNDERFILL = 0.7
ANGLE_SEARCH_MAX_RANGE = 2
ANGLE_SEARCH_SEARCH_RANGE = 1e-2

center_point = np.array([0.0, 0.0])

############################################ Utility ###################################################################

def angle_to_dir(phi, norm = MAX_RG):
    return norm * np.array([math.cos(phi), math.sin(phi)])

def dir_to_ray(dir, center = np.array([0, 0])):
    return LineString([(center[0], center[1]), (center[0] + dir[0], center[1] + dir[1])])

def dir_to_angle(dir):
    return np.arctan2(dir[1], dir[0])

def get_dis(p1, p2):
    return np.sqrt(np.sum(np.square(p1 - p2)))

def circle_to_poly(cir):
    p, r = cir
    delta = 2 * math.pi / CIR_DIVID_NUM
    points = [angle_to_dir(delta * i, r) + p for i in range(CIR_DIVID_NUM)]
    return Polygon([(point[0], point[1]) for point in points])

def angle_to_half_space(phi):
    v1 = angle_to_dir(phi)
    v2 = angle_to_dir(phi + math.pi / 2)
    return Polygon([(-v2[0], -v2[1]), (-v2[0] + v1[0], -v2[1] + v1[1]), (v2[0] + v1[0], v2[1] + v1[1]), (v2[0], v2[1])])

def get_angle_dis(start_angle, end_angle):
    return dir_to_angle(np.array([np.cos(end_angle - start_angle), np.sin(end_angle - start_angle)]))

def two_angles_to_sector(start_angle, end_angle):
    middle_angle = (start_angle + end_angle) / 2
    start_ray = angle_to_dir(start_angle)
    end_ray = angle_to_dir(end_angle)
    middle_ray = angle_to_dir(middle_angle)
    return Polygon([(0, 0), (start_ray[0], start_ray[1]), (middle_ray[0], middle_ray[1]), (end_ray[0], end_ray[1])])

############################################ Path Class ######################################################
class Path:
    def __init__ (self, _points):
        self.points = _points

    def __str__(self):
        return 'closed ' + str(len(self.points)) + '\n' +\
            '\n'.join([' '.join([str(point[0] + center_point[0]), str(point[1] + center_point[1]), str(radius)]) for point, radius in self.points])

    def split(self, split_ray, is_forward=True):
        inter_label = -1
        intersect = None
        for i in range(len(self.points) - 1):
            cur_line = LineString([(self.points[i][0][0], self.points[i][0][1]),
                                   (self.points[i + 1][0][0], self.points[i + 1][0][1])])
            intersect = cur_line.intersection(split_ray)
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
            split_points = self.points[:inter_label + 1]
            split_points.append((intersect, inter_radius))
            return Path(split_points)


    def link(self, path):
        if get_dis(self.points[-1][0], path.points[0][0]) < 1e-9:
            self.points.extend(path.points[1:])
        else:
            self.points.extend(path.points)

    def reverse(self):
        self.points.reverse()

    def get_poly(self, keep_head=True):
        poly = circle_to_poly(self.points[0])
        for i in range(len(self.points) - 1):
            p1, r1 = self.points[i]
            p2, r2 = self.points[i+1]

            # compute the polygon from circle i to circle i+1
            dis = get_dis(p1, p2)
            if math.fabs(r2 - r1) < dis:
                phi = math.pi/2 + np.arcsin((r2 - r1) / dis)
                phi0 = dir_to_angle(p2 - p1)
                v1_u = angle_to_dir(phi0 + phi, r1) + p1
                v1_d = angle_to_dir(phi0 - phi, r1) + p1
                v2_u = angle_to_dir(phi0 + phi, r2) + p2
                v2_d = angle_to_dir(phi0 - phi, r2) + p2
                poly = poly.union(Polygon([(v1_u[0], v1_u[1]), (p1[0], p1[1]),
                                    (v1_d[0], v1_d[1]), (v2_d[0], v2_d[1]),
                                    (p2[0], p2[1]), (v2_u[0], v2_u[1])]))

            # union the next circle
            poly = poly.union(circle_to_poly(self.points[i+1]))

            #plt.plot(*poly.exterior.xy)
            #plt.show()

        if keep_head:
            return poly
        else:
            return poly.difference(circle_to_poly(self.points[0]))

    def get_inner_poly(self):
        inner_poly = Polygon([(point[0], point[1]) for point, radius in self.points])
        return inner_poly.union(self.get_poly())

class Spiral:
    def __init__(self, path):
        self.path = path
        self.poly = path.get_poly()
        self.area = self.poly.area

    def __str__(self):
        return str(self.path)

    def extend(self, path):
        self.path.link(path)
        new_poly = path.get_poly(False)
        self.area = self.area + new_poly.area
        self.poly = self.poly.union(new_poly)

    def plot(self):
        gpd.GeoSeries([self.poly]).plot()
        plt.show()

def get_contour(points):
    points = [((point - center_point), radius) for point, radius in points] # Shift to center
    points.append(points[0])
    split_ray = dir_to_ray(angle_to_dir(SPLIT_ANG))

    contour_path = Path(points)
    contour = contour_path.split(split_ray)
    contour.link(contour_path.split(split_ray, False))
    contour.reverse()
    #print(contour)
    #gpd.GeoSeries([contour.get_poly()]).plot()
    #plt.show()
    return contour

def find_center(points):
    lb = min([point[0] for point, radius in points])
    rb = max([point[0] for point, radius in points])
    db = min([point[1] for point, radius in points])
    ub = max([point[1] for point, radius in points])

    return np.array([(lb + rb) / 2, (ub + db) / 2])

def read_shape(f):
    f.readline()
    f.readline()
    n = int(f.readline())
    points = []
    for i in range(n):
        s = f.readline().split(' ')
        points.append((float(s[0]) - center_point[0], float(s[1]) - center_point[1]))

    return Polygon(points)


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

def get_spiral_from_cut(cur_contour, next_contour, new_angle, next_angle):
    split_ray_1 = dir_to_ray(angle_to_dir(new_angle))
    path = cur_contour.split(split_ray_1, False)
    split_ray_2 = dir_to_ray(angle_to_dir(next_angle))
    path.link(next_contour.split(split_ray_2))
    return path

def get_cut_poly(path):
    opp_split_ray = dir_to_ray(angle_to_dir(SPLIT_ANG + math.pi))
    split_path = path.split(opp_split_ray)
    half_space_up = angle_to_half_space(SPLIT_ANG + math.pi / 2)
    half_space_down = angle_to_half_space(SPLIT_ANG + 3 * math.pi / 2)
    if split_path != None:
        split_path_down = path.split(opp_split_ray, False)
        middle_poly = Path([split_path_down.points[-2], split_path.points[1]]).get_poly()
        return middle_poly.union(half_space_up.intersection(split_path.get_poly()))\
                          .union(half_space_down.intersection(split_path_down.get_poly()))
    else:
        path_poly = path.get_poly()
        poly_up = half_space_up.intersection(path_poly)
        poly_down = half_space_down.intersection(path_poly)
        if poly_up.area > poly_down.area:
            return poly_up
        else:
            return poly_down

def get_cut_area(path):
    return get_cut_poly(path).area

def compute_spiral(contours, alg_init, alg_next_cut, alg_last_cut):
    # initial path
    init_split_ray = dir_to_ray(angle_to_dir(INIT_ANG))
    spiral = Spiral(contours[0].split(init_split_ray))

    prev_state = alg_init(spiral)
    for i in range(len(contours) - 1):
        new_angle, next_angle, prev_state = alg_next_cut(contours[i], contours[i + 1], prev_state)
        print(new_angle, next_angle)

        # compute spiral
        new_spiral_path = get_spiral_from_cut(contours[i], contours[i + 1], new_angle, next_angle)
        spiral.extend(new_spiral_path)

        #print(new_spiral_path, contours[i])
        #gpd.GeoSeries([new_spiral_path.get_poly()]).plot()
        #plt.show()

    last_ang = alg_last_cut(contours[-1], prev_state)
    last_split_ray = dir_to_ray(angle_to_dir(last_ang))
    spiral.extend(contours[-1].split(last_split_ray, False))

    #gpd.GeoSeries([spiral_poly_test]).plot()
    #plt.show()

    return spiral

################################ Angle-Search-Algorithm ##########################################

def get_angle_range(prev_state):
    prev_angle, prev_point, prev_poly = prev_state
    end_angle = dir_to_angle(prev_point[0])
    radius = prev_point[1]
    dis = get_dis(np.array([0, 0]), prev_point[0])
    angle_range = ANGLE_SEARCH_MAX_RANGE
    if dis > radius:
        angle_range = min(angle_range, get_angle_dis(end_angle, prev_angle) + 2 * np.arcsin(radius / dis))
    start_angle = end_angle + angle_range

    return start_angle, end_angle, angle_range

def get_eva(start_angle, new_angle, next_angle, end_angle, cur_contour, next_contour, prev_poly):
    new_spiral_path = get_spiral_from_cut(cur_contour, next_contour, new_angle, next_angle)
    sector = two_angles_to_sector(start_angle, end_angle)

    # Get spiral poly
    new_spiral_poly = sector.intersection(new_spiral_path.get_poly())

    # Get total poly
    new_spiral_path.link(Path([(np.array([0.0, 0.0]), 1.0)]))
    total_poly = sector.intersection(new_spiral_path.get_inner_poly())

    # Compute under/overfill
    underfill = total_poly.difference(new_spiral_poly.union(prev_poly)).area
    overfill = new_spiral_poly.intersection(prev_poly).area
    return underfill * ANGLE_SEARCH_WEIGHT_UNDERFILL + overfill * (1 - ANGLE_SEARCH_WEIGHT_UNDERFILL)

def get_eva_last(start_angle, last_angle, end_angle, last_contour, total_poly, prev_poly):
    last_split_ray = dir_to_ray(angle_to_dir(last_angle))
    last_spiral_path = last_contour.split(last_split_ray, False)
    sector = two_angles_to_sector(start_angle, end_angle)

    # Get spiral poly
    last_spiral_poly = sector.intersection(last_spiral_path.get_poly())

    # Compute under/overfill
    underfill = total_poly.difference(last_spiral_poly.union(prev_poly)).area
    overfill = last_spiral_poly.intersection(prev_poly).area
    return underfill * ANGLE_SEARCH_WEIGHT_UNDERFILL + overfill * (1 - ANGLE_SEARCH_WEIGHT_UNDERFILL)

def get_best_next_angle(start_angle, new_angle, end_angle, cur_contour, next_contour, prev_poly):
    next_angle = new_angle
    next_range = get_angle_dis(end_angle, next_angle)
    while next_range > ANGLE_SEARCH_SEARCH_RANGE:
        next_m1 = next_angle - next_range / 3
        next_m2 = next_angle - 2 * next_range / 3

        m1_eva = get_eva(start_angle, new_angle, next_m1, end_angle, cur_contour, next_contour, prev_poly)
        m2_eva = get_eva(start_angle, new_angle, next_m2, end_angle, cur_contour, next_contour, prev_poly)

        next_range = next_range * 2 / 3
        if m1_eva > m2_eva:         # m2 gets better result
            next_angle = next_m1

    next_angle = next_angle - next_range / 2
    eva = get_eva(start_angle, new_angle, next_angle, end_angle, cur_contour, next_contour, prev_poly)
    return next_angle, eva

def angle_search_init(spiral):
    prev_point = spiral.path.points[0]
    prev_angle = dir_to_angle(prev_point[0])
    prev_poly = spiral.poly
    return prev_angle, prev_point, prev_poly

def angle_search_next_cut(cur_contour, next_contour, prev_state):
    prev_angle, prev_point, prev_poly = prev_state
    start_angle, end_angle, angle_range = get_angle_range(prev_state)
    print(start_angle, end_angle, angle_range)

    # Cut everything to this range
    start_ray = dir_to_ray(angle_to_dir(start_angle + 0.05))
    end_ray = dir_to_ray(angle_to_dir(end_angle - 0.05))
    cur_contour_s = cur_contour.split(start_ray).split(end_ray, False)
    next_contour_s = next_contour.split(start_ray).split(end_ray, False)
    prev_poly_s = prev_poly.intersection(two_angles_to_sector(start_angle, end_angle))

    # Search for the best
    new_angle = start_angle
    new_range = angle_range
    while new_range > ANGLE_SEARCH_SEARCH_RANGE:
        new_m1 = new_angle - new_range / 3
        new_m2 = new_angle - 2 * new_range / 3

        _1, m1_eva = get_best_next_angle(start_angle, new_m1, end_angle, cur_contour_s, next_contour_s, prev_poly_s)
        _2, m2_eva = get_best_next_angle(start_angle, new_m2, end_angle, cur_contour_s, next_contour_s, prev_poly_s)

        print("Test case (", new_m1, ",", _1, "): eva = ", m1_eva)
        print("Test case (", new_m2, ",", _2, "): eva = ", m2_eva)

        new_range = new_range * 2 / 3
        if m1_eva > m2_eva:         # m2 gets better result
            new_angle = new_m1

    new_angle = new_angle - new_range / 2
    next_angle, _ = get_best_next_angle(start_angle, new_angle, end_angle, cur_contour_s, next_contour_s, prev_poly_s)
    print("Eva:", _)

    # Update prev_state
    split_ray = dir_to_ray(angle_to_dir(next_angle))
    next_prev_point = next_contour.split(split_ray).points[0]
    spiral = get_spiral_from_cut(cur_contour, next_contour, new_angle, next_angle)
    prev_state = (new_angle, next_prev_point, spiral.get_inner_poly())

    return new_angle, next_angle, prev_state

def angle_search_last_cut(last_contour, prev_state):
    prev_angle, prev_poiont, prev_poly = prev_state
    start_angle, end_angle, angle_range = get_angle_range(prev_state)

    # Cut everything to this range
    start_ray = dir_to_ray(angle_to_dir(start_angle + 0.05))
    end_ray = dir_to_ray(angle_to_dir(end_angle - 0.05))
    last_contour_s = last_contour.split(start_ray).split(end_ray, False)
    sector = two_angles_to_sector(start_angle, end_angle)
    prev_poly_s = sector.intersection(prev_poly)
    total_poly = sector.intersection(last_contour.get_inner_poly())

    # Search for the best
    last_angle = start_angle
    last_range = angle_range
    while last_range > ANGLE_SEARCH_SEARCH_RANGE:
        last_m1 = last_angle - last_range / 3
        last_m2 = last_angle - 2 * last_range / 3

        m1_eva = get_eva_last(start_angle, last_m1, end_angle, last_contour_s, total_poly, prev_poly_s)
        m2_eva = get_eva_last(start_angle, last_m2, end_angle, last_contour_s, total_poly, prev_poly_s)

        last_range = last_range * 2 / 3
        if m1_eva > m2_eva:         # m2 gets better result
            last_angle = last_m1

    last_angle = last_angle - last_range / 2
    return last_angle

def compute_spiral_angle_search(contours):
    return compute_spiral(contours, angle_search_init, angle_search_next_cut, angle_search_last_cut)


################################# Right-Angle-Algorithm ##########################################

def get_new_cut_right_angle(prev_state):
    prev_point, radius = prev_state
    dis = get_dis(np.array([0, 0]), prev_point)
    prev_angle = dir_to_angle(prev_point)
    new_angle = 0.0
    delta = radius
    if dis < delta:
        new_angle = prev_angle + math.pi * 0.7
    else:
        new_angle = prev_angle + 2 * np.arcsin(delta / dis)

    return new_angle

def right_angle_init(spiral):
    return spiral.path.points[0]

def right_angle_next_cut(cur_contour, next_contour, prev_state):
    new_angle = get_new_cut_right_angle(prev_state)
    cut_angle = new_angle - RIGHT_ANGLE_ANG
    split_ray_1 = dir_to_ray(angle_to_dir(new_angle))
    path = cur_contour.split(split_ray_1, False)
    split_ray_2 = dir_to_ray(angle_to_dir(cut_angle), path.points[-1][0])
    path = next_contour.split(split_ray_2)
    return new_angle, dir_to_angle(path.points[0][0]), path.points[0]

def right_angle_last_cut(last_contour, prev_state):
    return get_new_cut_right_angle(prev_state)

def compute_spiral_right_angle(contours):
    return compute_spiral(contours, right_angle_init, right_angle_next_cut, right_angle_last_cut)


############################################## MAIN ##############################################

# python ./to_spiral.py ./data/test-0.txt ./data/test-0-small.out

if __name__ == '__main__':
    shape_file = sys.argv[1]
    contour_file = sys.argv[2]
    alg = sys.argv[3]
    output_file = shape_file + '.spiral'
    if alg == 's':
        ANGLE_SEARCH_WEIGHT_UNDERFILL = float(sys.argv[4])
        output_file = shape_file + '-s-w' + str(ANGLE_SEARCH_WEIGHT_UNDERFILL) + '-p' + str(ANGLE_SEARCH_SEARCH_RANGE) + '.spiral'
    if alg == 'r':
        RIGHT_ANGLE_ANG = float(sys.argv[4])
        output_file = shape_file + '-r-a' + str(RIGHT_ANGLE_ANG) + '.spiral'

    contours = []
    with open(contour_file, 'r', encoding='utf-8') as f:
        contours = get_contours(f)


    with open(shape_file, 'r', encoding='utf-8') as f:
        shape_poly = read_shape(f)

    if alg == 'r':
        spiral = compute_spiral_right_angle(contours)
    else:
        spiral = compute_spiral_angle_search(contours)

    print("The total area of the shape", shape_poly.area)
    print("The total area of the spiral path:", spiral.area)
    print("The total area of the spiral path (union):", spiral.poly.area)
    print("Rate of overfill:", (spiral.area - spiral.poly.area) / shape_poly.area * 100)
    print("Rate of underfill:", 100 - (spiral.poly.intersection(shape_poly).area) / shape_poly.area * 100)

    with open(output_file, 'w', encoding='utf-8') as f:
        f.write(str(spiral))

    #spiral.plot()
