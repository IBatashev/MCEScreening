import numpy as np
from scipy import interpolate

def analyze(volumes, moments, volume_undeformed, moment_undeformed):
    """Function returns magneto-elastic response value and fit confidence when supplied with a list of magnetic
     moments and a corresponding list of deformed volumes, as well as moment and volume of the undeformed state."""

    # normalizing moments and volumes
    x = list(np.array(volumes) / volume_undeformed)
    y = list(np.array(moments) / moment_undeformed)
    discretization = 50  # how many fitted points to add between measured points.
    points = np.column_stack([x, y])

# Case 1 linear fit
    Mag_el, fit_quality, path = linfit(x, y)

    if fit_quality < 0.85:

        # interpolating via a Bezier fit
        path = evaluate_bezier(points, discretization)
        path = np.unique(path, axis=0)
        path_x, path_y = path[:, 0], path[:, 1]
        b = interpolate.InterpolatedUnivariateSpline(path_x, path_y, k=4)
        b2 = interpolate.InterpolatedUnivariateSpline(path_x, path_y, k=5)
        # cutting top and bottom 25% of moment change range
        minimal = np.min(path_y)
        maximal = np.max(path_y)
        span = abs(maximal - minimal)
        mask_size = 0.25 * span
        mask_bottom = minimal + mask_size
        mask_top = maximal - mask_size
        masked_array = np.array(path)
        masked_array[masked_array[:, 1] < mask_bottom] = 0
        masked_array[masked_array[:, 1] > mask_top] = 0
        # gather and count the remaining parts of the curve
        segments_ids = np.where(masked_array[:, 0] != 0)[0]
        segments = np.split(masked_array[segments_ids], np.where(np.diff(segments_ids) != 1)[0]+1)
        number_of_segments = len(segments)

        if number_of_segments == 1:
# Case 2 attempt linear fit on the only segment
            initial_segment = np.array(segments[0])
            Mag_el, fit_quality, path_line = linfit(initial_segment[:, 0], initial_segment[:, 1])
        elif number_of_segments > 1:
            # Chose which of the remaining segments is the closest to the undeformed structure
            distances_to_undeformed = []
            for seg in segments:
                closest = find_nearest(seg[:, 0], 1)
                distances_to_undeformed.append(closest)
            value_of_closest_one = find_nearest(distances_to_undeformed, 1)
            index_of_closest_one = int(np.where(distances_to_undeformed == value_of_closest_one)[0])
            initial_segment = np.array(segments[index_of_closest_one])
# Case 3 attempt linear fit on this segment
            Mag_el, fit_quality, path_line = linfit(initial_segment[:, 0], initial_segment[:, 1])

            # if first approach does not work well enough, attempt to further find the linear parts
            # break on special points (first try in maximums, then on inflections) and take longest line (so pythogorus between edge points of new even smaller segments)
        if fit_quality < 0.9:
            seg_x_max = float(np.max(initial_segment[:, 0]))
            seg_x_min = float(np.min(initial_segment[:, 0]))

            extremum_xs = b.derivative(n=1).roots()
            extremum_ys = b(extremum_xs)
            extrema = np.column_stack((extremum_xs, extremum_ys))

            breakpoints = np.array(extrema)
            breakpoints = breakpoints[(seg_x_min < breakpoints[:, 0]) & (breakpoints[:, 0] < seg_x_max)]

            inflect = b2.derivative(n=2).roots()
            inflect_xs = inflect[1:-1]  # remove edge points
            inflect_ys = b(inflect_xs)
            inflections = np.column_stack((inflect_xs, inflect_ys))

            if len(breakpoints) == 0:
                breakpoints = np.array(inflections)
                breakpoints = breakpoints[(seg_x_min < breakpoints[:, 0]) & (breakpoints[:, 0] < seg_x_max)]
                if len(breakpoints) == 0:
                    print('no special points, but linear fit is bad. Something went very wrong')

            split_at = initial_segment[:, 0].searchsorted(breakpoints[:, 0])
            secondary_segments = np.split(initial_segment, split_at)

            # Chose which of the remaining segments is the closest to the undeformed structure
            distances_to_undeformed = []
            for seg in secondary_segments:
                closest = find_nearest(seg[:, 0], 1)
                distances_to_undeformed.append(closest)
            value_of_closest_one = find_nearest(distances_to_undeformed, 1)
            index_of_closest_one = int(np.where(distances_to_undeformed == value_of_closest_one)[0])
# Case 4 attempt linear fit on this segment
            Mag_el, fit_quality_third, path_line = linfit(secondary_segments[index_of_closest_one][:, 0], secondary_segments[index_of_closest_one][:, 1])

    return Mag_el, fit_quality


def evaluate_bezier(points, n):
    """Interpolation based on a cubic Bezier curve"""
    curves = get_bezier_cubic(points)
    return np.array([fun(t) for fun in curves for t in np.linspace(0, 1, n)])


def get_bezier_cubic(points):
    """Return one cubic curve for each consecutive points"""
    A, B = get_bezier_coef(points)
    return [
        get_cubic(points[i], A[i], B[i], points[i + 1])
        for i in range(len(points) - 1)]


def get_cubic(a, b, c, d):
    """Returns the general Bezier cubic formula given 4 control points"""
    return lambda t: np.power(1 - t, 3) * a + 3 * np.power(1 - t, 2) * t * b + 3 * (1 - t) * np.power(t, 2) * c + np.power(t, 3) * d


def get_bezier_coef(points):
    """Find the a & b points for Bezier fit"""
    n = len(points) - 1
    # build coefficents matrix
    C = 4 * np.identity(n)
    np.fill_diagonal(C[1:], 1)
    np.fill_diagonal(C[:, 1:], 1)
    C[0, 0] = 2
    C[n - 1, n - 1] = 7
    C[n - 1, n - 2] = 2
    # build points vector
    P = [2 * (2 * points[i] + points[i + 1]) for i in range(n)]
    P[0] = points[0] + 2 * points[1]
    P[n - 1] = 8 * points[n - 1] + points[n]
    # solve system, find a & b
    A = np.linalg.solve(C, P)
    B = [0] * n
    for i in range(n - 1):
        B[i] = 2 * points[i + 1] - A[i + 1]
    B[n - 1] = (A[n - 1] + points[n]) / 2
    return A, B


def linfit(x, y):
    """Function to perform a linear fit"""
    coeffs = np.polyfit(x, y, 1)
    poly1d_fn = np.poly1d(coeffs)
    yhat = poly1d_fn(x)
    ybar = np.sum(y) / len(y)
    ssreg = np.sum((yhat - ybar) ** 2)
    sstot = np.sum((y - ybar) ** 2)
    R_square = ssreg / sstot     # fit quality
    xp = np.linspace(min(x), max(x), 50)
    yp = poly1d_fn(xp)
    path = np.column_stack((xp, yp))
    return coeffs[0], R_square, path


def find_nearest(array, value):
    """Function for determining the which curve segment is closest to undeformed volume """
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]


# ------------------------------------------------------------------------------------------------------- #
#  _____                                           _       _____ _             _     _   _                #
# /  __ \                                         | |     /  ___| |           | |   | | | |               #
# | /  \/ ___  _ __ ___  _ __ ___   __ _ _ __   __| |___  \ `--.| |_ __ _ _ __| |_  | |_| | ___ _ __ ___  #
# | |    / _ \| '_ ` _ \| '_ ` _ \ / _` | '_ \ / _` / __|  `--. \ __/ _` | '__| __| |  _  |/ _ \ '__/ _ \ #
# | \__/\ (_) | | | | | | | | | | | (_| | | | | (_| \__ \ /\__/ / || (_| | |  | |_  | | | |  __/ | |  __/ #
#  \____/\___/|_| |_| |_|_| |_| |_|\__,_|_| |_|\__,_|___/ \____/ \__\__,_|_|   \__| \_| |_/\___|_|  \___| #
#                                                                                                         #
# ------------------------------------------------------------------------------------------------------- #


volumes = [198.18, 208.61, 230.57, 207.56]
moments = [-0.209, 25.351, 28.175, 25.428]

volume_undeformed = 219.59
moment_undeformed = 26.24

print(analyze(volumes, moments, volume_undeformed, moment_undeformed))
