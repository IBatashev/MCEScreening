import numpy as np
import matplotlib.pyplot as plt
import os
import tqdm

from scipy import interpolate

# outdir = 'D:/MCES/MP/outdir'
outdir = 'D:/MCES/ICSD/outdir'
vasp_results_dir = outdir


def calculate_mag_field(moment, volume):
    """Takes cell moment in [mB] and volume in [A^3] and returns value for internal magnetic field in [T]"""

    pi = 3.141592653
    mu0 = 4*pi*10**(-7)                         # Vacuum permeability in H/m
    mB = 9.2741*10**(-24)                       # Bohr magneton value in J/T
    field = (mu0*moment*mB)/(volume*10**(-30))  # Formula for internal magnetic field in Tesla
    return field


def moment_volume_after(ID, deformation):
    """Read OUTCAR file for chosen deformation type and return moment, volume and
    value of magnetic field based on moment and volume.
    If chosen deformation is n/a for this ID an empty string will be returned"""

    if deformation in os.listdir(vasp_results_dir + '/' + str(ID)):
        with open(vasp_results_dir + '/' + str(ID) + '/' + deformation + "/OUTCAR") as search:
            moment_lines = []
            for line in search:
                if 'volume of cell' in line:
                    volume = (float(line.split()[-1]))
                if 'tot' in line:
                    moment_lines.append(line)
        moment = float(moment_lines[-2].split()[-1])
        field_after = calculate_mag_field(moment, volume)
        return moment, volume, field_after
    else:
        return '', '', ''


# find the a & b points
def get_bezier_coef(points):
    # since the formulas work given that we have n+1 points
    # then n must be this:
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


# returns the general Bezier cubic formula given 4 control points
def get_cubic(a, b, c, d):
    return lambda t: np.power(1 - t, 3) * a + 3 * np.power(1 - t, 2) * t * b + 3 * (1 - t) * np.power(t, 2) * c + np.power(t, 3) * d


# return one cubic curve for each consecutive points
def get_bezier_cubic(points):
    A, B = get_bezier_coef(points)
    return [
        get_cubic(points[i], A[i], B[i], points[i + 1])
        for i in range(len(points) - 1)
    ]


# evalute each cubic curve on the range [0, 1] sliced in n points
def evaluate_bezier(points, n):
    curves = get_bezier_cubic(points)
    return np.array([fun(t) for fun in curves for t in np.linspace(0, 1, n)])


def linfit(x, y):
    coeffs = np.polyfit(x, y, 1)
    poly1d_fn = np.poly1d(coeffs)

    # fit quality
    yhat = poly1d_fn(x)  # or [p(z) for z in x]
    ybar = np.sum(y) / len(y)  # or sum(y)/len(y)
    ssreg = np.sum((yhat - ybar) ** 2)  # or sum([ (yihat - ybar)**2 for yihat in yhat])
    sstot = np.sum((y - ybar) ** 2)  # or sum([ (yi - ybar)**2 for yi in y])
    R_square = ssreg / sstot

    xp = np.linspace(min(x), max(x), 50)
    yp = poly1d_fn(xp)
    path = np.column_stack((xp, yp))

    return coeffs[0], R_square, path


#def get_xy(ID):
#    moment, volume, field = moment_volume_after(ID, 'undeformed')
#    moment_ll = moment_volume_after(ID, 'a_dec')[0] / moment
#    moment_l = moment_volume_after(ID, 'V_dec')[0] / moment
#    moment_h = moment_volume_after(ID, 'V_inc')[0] / moment
#    moment_hh = moment_volume_after(ID, 'a_inc')[0] / moment
#    y = [abs(moment_ll), abs(moment_l), 1, abs(moment_h), abs(moment_hh)]
#    x = [0.84, 0.95, 1, 1.05, 1.16]
#    pp = np.column_stack([x, y])
#    return pp


def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]


def removearray(L, array):
    ind = 0
    size = len(L)
    while ind != size and not np.array_equal(L[ind], array):
        ind += 1
    if ind != size:
        L.pop(ind)
    else:
        raise ValueError('array not found in list.')


def analyze(x=None, y=None, ID='', printer=False, plot=False):

    discretization = 50  # how many fitted points to add between measured points.
    if ID != '':
        points = get_xy(ID)
        x, y = points[:, 0], points[:, 1]
    elif x is None or y is None:
        print('Please submit either a valid ID or list of volumes and list of moments')
    else:
        points = np.column_stack([x, y])

    slope, fit_quality, path = linfit(x, y)
    fit_quality_second = 0
    fit_quality_third = 0
    path_line = path

    if fit_quality < 0.85:
        path = evaluate_bezier(points, discretization)
        path = np.unique(path, axis=0)
        path_x, path_y = path[:, 0], path[:, 1]        # extract x & y coordinates of points
        b = interpolate.InterpolatedUnivariateSpline(path_x, path_y, k=4)
        b2 = interpolate.InterpolatedUnivariateSpline(path_x, path_y, k=5)  # this is slow :( better use numpy
        minimal = np.min(path_y)
        maximal = np.max(path_y)
        span = abs(maximal - minimal)
        mask_size = 0.25 * span
        mask_bottom = minimal + mask_size
        mask_top = maximal - mask_size
        masked_array = np.array(path)
        masked_array[masked_array[:, 1] < mask_bottom] = 0
        masked_array[masked_array[:, 1] > mask_top] = 0

        segments_ids = np.where(masked_array[:, 0] != 0)[0]
        segments = np.split(masked_array[segments_ids], np.where(np.diff(segments_ids) != 1)[0]+1)

        number_of_segments = len(segments)

        if number_of_segments == 1:
            initial_segment = np.array(segments[0])
            slope, fit_quality, path_line = linfit(initial_segment[:, 0], initial_segment[:, 1])
        elif number_of_segments > 1:
            distances_to_undeformed = []
            for seg in segments:
                closest = find_nearest(seg[:, 0], 1)
                distances_to_undeformed.append(closest)
            value_of_closest_one = find_nearest(distances_to_undeformed, 1)
            index_of_closest_one = int(np.where(distances_to_undeformed == value_of_closest_one)[0])
            initial_segment = np.array(segments[index_of_closest_one])
            slope, fit_quality, path_line = linfit(initial_segment[:, 0], initial_segment[:, 1])

            # if first approach does not work well enough, we go deeper
            # break on special points (first try in maximums, then on inflections) and take longest line (so pythagoras between edge points of new even smaller segments)
        if fit_quality< 0.9:
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

            distances_to_undeformed = []
            for seg in secondary_segments:
                closest = find_nearest(seg[:, 0], 1)
                distances_to_undeformed.append(closest)
            value_of_closest_one = find_nearest(distances_to_undeformed, 1)
            index_of_closest_one = int(np.where(distances_to_undeformed == value_of_closest_one)[0])
            slope, fit_quality, path_line = linfit(secondary_segments[index_of_closest_one][:, 0], secondary_segments[index_of_closest_one][:, 1])

    if printer == True:
        print('initial linear fit r^2: ', fit_quality)
        print('Final value of Mel', slope)

    if plot == True:
        fig, ax1 = plt.subplots()
        if slope != 0:
            path_x, path_y = path[:, 0], path[:, 1]
            path_line_x, path_line_y = path_line[:, 0], path_line[:, 1]
            ax2 = ax1.twinx()
            coloring = np.random.rand(3)
            ax2.plot(path_x, path_y, color=coloring, alpha=0.5)
            ax2.plot(path_line_x, path_line_y, '--', color='red', alpha=0.3)
            ax2.axes.yaxis.set_visible(True)
            ax1.axes.yaxis.set_visible(False)
            ax1.axhspan(0, 0.5, color='black', alpha=0.3)
            ax1.axhspan(1.5, 2, color='black', alpha=0.3)
            for i in x:
                ax1.axvline(i, color='black', alpha=0.2)
            plt.show()

    # get distance from undeformed state to the step
    path_center = np.take(path_line, path_line.size // 2)
    step_distance = (path_center - 1)

    return path, path_line, slope, fit_quality, step_distance


def multiplot(id_list):

        fig, ax1 = plt.subplots()
        ax1.plot(1, 1, 'ro')
        # ax1.axhline(1.5, color='black', alpha=0.2)
        # ax1.axhline(0.5, color='black', alpha=0.2)
        # ax1.axhline(2, color='black', alpha=0)
        # ax1.axhline(0, color='black', alpha=0)

        ax1.axhspan(0, 0.5, color='black', alpha=0.3)
        ax1.axhspan(1.5, 2, color='black', alpha=0.3)

        ax1.axes.yaxis.set_visible(False)
        ax1.axvline(0.84, color='black', alpha=0.2)
        ax1.axvline(0.95, color='black', alpha=0.2)
        ax1.axvline(1.0, color='black', alpha=0.2)
        ax1.axvline(1.05, color='black', alpha=0.2)
        ax1.axvline(1.16, color='black', alpha=0.2)

        for idx in id_list:

            path, path_line, slope, fit_quality, step_distance = analyze(idx)

            if slope != 0:
                # x, y = points[:, 0], points[:, 1]
                path_x, path_y = path[:, 0], path[:, 1]
                path_line_x, path_line_y = path_line[:, 0], path_line[:, 1]
                ax2 = ax1.twinx()
                coloring = np.random.rand(3)
                ax2.plot(path_x, path_y, color=coloring, alpha=0.5)
                # ax2.plot(x, y, marker='x', color=coloring, linestyle = 'None')

                # ax2.plot(extrema[:, 0], extrema[:, 1], marker='x', color=coloring, linestyle = 'None')

                #ax2.plot(inflections[:, 0], inflections[:, 1], 'xr')

                ax2.plot(path_line_x, path_line_y, '--', color='red', alpha=0.3)
                ax2.axes.yaxis.set_visible(False)

                # ax1.plot(x, y, 'ro', alpha=0.3)
        plt.show()
#
#
# with open('all_cub.txt') as f:
#     ids = [line.rstrip('\n') for line in f]
#
# # secondary fit problems
# ids = ['mp-10885',
# 'mp-5318',
# 'mp-1184221',
# 'mp-28525',
# 'mp-531715',
# 'mp-568983',
# # 'mp-38489',
# # 'mp-754077',
# # 'mp-1063914',
# # 'mp-18695',
# # 'mp-1193908',
# 'mp-535',]
#
# # multiplot(['mp-6670', 'mp-999552', 'mp-582028','mp-6274'])
# multiplot(['7676']) #417737
# multiplot(ids)


# fig, ax1 = plt.subplots()
# ax2 = ax1.twinx()
#
# xxx = np.linspace(np.min(x), np.max(x), 20)
# ax1.plot(xxx, poly1d_lin(xxx), '--')
# ax1.plot(x, y, 'ro')
#
# ax1.plot(crit_points, crit_vals, 'xb', label='extrema')
# ax1.plot(inflect_points, inflect_vals, 'xr', label='inflection')
# # plt.plot(xx, f(xx), '-', label='Spline Parabolic')
#
# ax2.plot(px, b.derivative()(px), '-', label='S')
# ax1.plot(px, b(px), '-', label='Bezier Cubic')
# fig.legend(loc='lower right')
#
# plt.show()

#analyze(x=[112.82/125.01, 125.01/125.01, 137.82/125.01], y=[9.524/13.098, 13.098/13.098, 15.866/13.098], printer=True, plot=True) # a
#analyze(x=[0.899976175, 0.950007942,1,1.049992058,1.099984117], y=[1381,2004,1,3693,3722], printer=True, plot=True) # c
#analyze(x=[1.099997901,1.04999895,1,0.95000105,0.900002099], y=[1401,8,1,-1225,4436], printer=True, plot=True) # a


#analyze(x=[1,0.899980876,0.949990438,1.050009562,1.099987251], y=[0,0,0,5.233,7.598], printer=True, plot=True) # NbNi2 wrong order of datapoints
#analyze(x=[0.899980876,0.949990438,1,1.050009562,1.099987251], y=[0,0,1,654.125,949.75], printer=True, plot=True) # NbNi2 correct order of datapoints
