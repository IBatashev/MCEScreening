import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit as cf
import os

# outdir = 'D:/MCES/MP/outdir'
outdir = 'D:/MCES/MP/outdir_VVV'
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


def polyfit(ID, deformation, plot=True):
    if deformation + '_inc' in os.listdir(vasp_results_dir + '/' + str(ID)):
        moment_inc, volume_inc, field_inc = moment_volume_after(ID, deformation+'_inc')
        moment_dec, volume_dec, field_dec = moment_volume_after(ID, deformation + '_dec')
        moment, volume, field = moment_volume_after(ID, 'undeformed')
        y = [field_dec, field, field_inc]
        x = [0.95, 1, 1.05]

    results = {}

    lincoeffs = np.polyfit(x, y, 1)

    # lincoeffs = np.polyfit(x, y, degree)
    poly1d_lin = np.poly1d(lincoeffs)

     # Polynomial Coefficients
    results['fit coefficients'] = lincoeffs.tolist()

    # fit quality
    yhat = poly1d_lin(x)              # or [p(z) for z in x]
    ybar = np.sum(y)/len(y)          # or sum(y)/len(y)
    ssreg = np.sum((yhat-ybar)**2)   # or sum([ (yihat - ybar)**2 for yihat in yhat])
    sstot = np.sum((y - ybar)**2)    # or sum([ (yi - ybar)**2 for yi in y])
    results['R^2'] = ssreg / sstot

    parabcoeffs = np.polyfit(x, y, 2)
    poly1d_par = np.poly1d(parabcoeffs)


    # search extremum
    crit = poly1d_par.deriv().r
    extremum_x = float(crit[crit.imag == 0].real)
    extremum_y = poly1d_par(extremum_x)
    extremum = [extremum_x, extremum_y]

    # test = poly1d_fn.deriv(2)(r_crit) # can check if min or max if necessary

    if abs(extremum[0] - 0.95) > abs(1.05 - extremum[0]):
        x2 = [0.95, 1, extremum[0]]
        y2 = [field_dec, field, extremum[1]]
        # new_slope, new_fit = linfit(x, y)
        lincoeffs2 = np.polyfit(x2, y2, 1)
        poly1d_lin2 = np.poly1d(lincoeffs2)
    else:
        x2 = [extremum[0], 1, 1.05]
        y2 = [extremum[1], field, field_inc]
        # new_slope, new_fit = linfit(x, y)
        lincoeffs2 = np.polyfit(x2, y2, 1)
        poly1d_lin2 = np.poly1d(lincoeffs2)
    # else:
    #     poly1d_lin2 = poly1d_lin

    # if r_crit.size == 0:
    #     print("No Extremum")
    # elif 0.95 < r_crit < 1.05:
    #     print("Extremum inside def range at", np.round(r_crit, 3))
    # else:
    #     print('Extremum outside def range at', np.round(r_crit, 3))

    if plot == True:
        print('Plotting...')

        xp = np.linspace(0.9, 1.1, 50)

        fig, axs = plt.subplots(1, 3, figsize=(9, 3))
        axs[0].plot(x, y, 'yo', xp, poly1d_lin(xp), '--k')
        axs[1].plot(x, y, 'yo', xp, poly1d_par(xp), '--k')
        axs[2].plot(x, y, 'yo', xp, poly1d_lin2(xp), '--c')

        # # plotting the fit line and origninal points
        # axs[0].plot(x, y, 'yo', xp, poly1d_fn(xp), '--k')

        # plotting x mark for extremum
        x_min = extremum[0]  # [test > 0]
        y_min = extremum[1]
        axs[1].plot(x_min, y_min, 'xr')
        axs[2].plot(x_min, y_min, 'xr')

        # plotting labels
        axs[0].set_xlabel('$\it{V_{def}/V}$', ha="left", fontsize=16)
        axs[0].set_ylabel('$\it{Field, T}$', rotation=90, ha="right", fontsize=14)
        axs[1].set_xlabel('$\it{V_{def}/V}$', ha="left", fontsize=16)
        axs[1].set_ylabel('$\it{Field, T}$', rotation=90, ha="right", fontsize=14)
        axs[2].set_xlabel('$\it{V_{def}/V}$', ha="left", fontsize=16)
        axs[2].set_ylabel('$\it{Field, T}$', rotation=90, ha="right", fontsize=14)
        # set limits - needs modification
        xlimit = 0.9, 1.1
        ylimit = min(y)-0.3*min(y), max(y)+0.3*max(y)
        axs[0].set_ylim(ylimit)
        axs[0].set_xlim(xlimit)
        axs[1].set_ylim(ylimit)
        axs[1].set_xlim(xlimit)
        axs[2].set_ylim(ylimit)
        axs[2].set_xlim(xlimit)
        # show plot
        plt.show()

    return results


def polyfit_moment(ID, deformation, plot=True):
    if deformation + '_inc' in os.listdir(vasp_results_dir + '/' + str(ID)):
        moment_inc, volume_inc, field_inc = moment_volume_after(ID, deformation+'_inc')
        moment_dec, volume_dec, field_dec = moment_volume_after(ID, deformation + '_dec')
        moment, volume, field = moment_volume_after(ID, 'undeformed')
        y = [moment_dec/moment, 1, moment_inc/moment]
        x = [0.95, 1, 1.05]

    results = {}

    lincoeffs = np.polyfit(x, y, 1)

    # lincoeffs = np.polyfit(x, y, degree)
    poly1d_lin = np.poly1d(lincoeffs)

     # Polynomial Coefficients
    results['fit coefficients'] = lincoeffs.tolist()

    # fit quality
    yhat = poly1d_lin(x)              # or [p(z) for z in x]
    ybar = np.sum(y)/len(y)          # or sum(y)/len(y)
    ssreg = np.sum((yhat-ybar)**2)   # or sum([ (yihat - ybar)**2 for yihat in yhat])
    sstot = np.sum((y - ybar)**2)    # or sum([ (yi - ybar)**2 for yi in y])
    results['R^2'] = ssreg / sstot

    parabcoeffs = np.polyfit(x, y, 2)
    poly1d_par = np.poly1d(parabcoeffs)

    # search extremum
    crit = poly1d_par.deriv().r
    extremum_x = float(crit[crit.imag == 0].real)
    extremum_y = poly1d_par(extremum_x)
    extremum = [extremum_x, extremum_y]

    # test = poly1d_fn.deriv(2)(r_crit) # can check if min or max if necessary
    results2 = {}

    if abs(extremum[0] - 0.95) > abs(1.05 - extremum[0]):
        x2 = [0.95, 1, extremum[0]]
        y2 = [moment_dec/moment, 1, extremum[1]]
        # new_slope, new_fit = linfit(x, y)
        lincoeffs2 = np.polyfit(x2, y2, 1)
        poly1d_lin2 = np.poly1d(lincoeffs2)

        yhat = poly1d_lin2(x2)  # or [p(z) for z in x]
        ybar = np.sum(y2) / len(y2)  # or sum(y)/len(y)
        ssreg = np.sum((yhat - ybar) ** 2)  # or sum([ (yihat - ybar)**2 for yihat in yhat])
        sstot = np.sum((y2 - ybar) ** 2)  # or sum([ (yi - ybar)**2 for yi in y])
        results2['R^2'] = ssreg / sstot

    else:
        x2 = [extremum[0], 1, 1.05]
        y2 = [extremum[1], 1, moment_inc/moment]
        # new_slope, new_fit = linfit(x, y)
        lincoeffs2 = np.polyfit(x2, y2, 1)
        poly1d_lin2 = np.poly1d(lincoeffs2)

        yhat = poly1d_lin2(x2)  # or [p(z) for z in x]
        ybar = np.sum(y2) / len(y2)  # or sum(y)/len(y)
        ssreg = np.sum((yhat - ybar) ** 2)  # or sum([ (yihat - ybar)**2 for yihat in yhat])
        sstot = np.sum((y2 - ybar) ** 2)  # or sum([ (yi - ybar)**2 for yi in y])
        results2['R^2'] = ssreg / sstot

    # else:
    #     poly1d_lin2 = poly1d_lin

    # if r_crit.size == 0:
    #     print("No Extremum")
    # elif 0.95 < r_crit < 1.05:
    #     print("Extremum inside def range at", np.round(r_crit, 3))
    # else:
    #     print('Extremum outside def range at', np.round(r_crit, 3))

    if plot == True:
        print('Plotting...')

        xp = np.linspace(0.9, 1.1, 50)

        fig, axs = plt.subplots(1, 3, figsize=(15, 6))
        axs[0].plot(x, y, 'yo', xp, poly1d_lin(xp), '--k')
        axs[1].plot(x, y, 'yo', xp, poly1d_par(xp), '--k')
        axs[2].plot(x, y, 'yo', xp, poly1d_lin2(xp), '--c')

        # # plotting the fit line and origninal points
        # axs[0].plot(x, y, 'yo', xp, poly1d_fn(xp), '--k')

        # plotting x mark for extremum
        x_min = extremum[0]  # [test > 0]
        y_min = extremum[1]
        axs[1].plot(x_min, y_min, 'xr')
        axs[2].plot(x_min, y_min, 'xr')

        # plotting labels
        axs[0].set_xlabel('$\it{V_{def}/V}$', ha="left", fontsize=16)
        axs[0].set_ylabel('$\it{Moment}$', rotation=90, ha="right", fontsize=14)
        axs[1].set_xlabel('$\it{V_{def}/V}$', ha="left", fontsize=16)
        axs[1].set_ylabel('$\it{Moment}$', rotation=90, ha="right", fontsize=14)
        axs[2].set_xlabel('$\it{a_{def}/a}$', ha="left", fontsize=16)
        axs[2].set_ylabel('$\it{Moment}$', rotation=90, ha="right", fontsize=14)
        # set limits - needs modification
        xlimit = 0.9, 1.1
        ylimit = min(y)-0.3*min(y), max(y)+0.3*max(y)
        axs[0].set_ylim(ylimit)
        axs[0].set_xlim(xlimit)
        axs[1].set_ylim(ylimit)
        axs[1].set_xlim(xlimit)
        axs[2].set_ylim(ylimit)
        axs[2].set_xlim(xlimit)
        # show plot
        plt.show()

    return results


def polyfit_moment_5point(ID, plot=True):
    # if deformation + '_inc' in os.listdir(vasp_results_dir + '/' + str(ID)):

    # moment_inc, volume_inc, field_inc = moment_volume_after(ID, deformation+'_inc')
    # moment_dec, volume_dec, field_dec = moment_volume_after(ID, deformation + '_dec')
    moment, volume, field = moment_volume_after(ID, 'undeformed')

    moment_ll = moment_volume_after(ID, 'a_dec')[0]/moment
    moment_l = moment_volume_after(ID, 'V_dec')[0]/moment
    moment_h = moment_volume_after(ID, 'V_inc')[0]/moment
    moment_hh = moment_volume_after(ID, 'a_inc')[0]/moment

    y = [moment_ll, moment_l, 1, moment_h, moment_hh]
    x = [0.84, 0.95, 1, 1.05, 1.16]

    results = {}
    lincoeffs = np.polyfit(x, y, 1)

    # lincoeffs = np.polyfit(x, y, degree)
    poly1d_lin = np.poly1d(lincoeffs)

     # Polynomial Coefficients
    results['fit coefficients'] = lincoeffs.tolist()

    # fit quality
    yhat = poly1d_lin(x)              # or [p(z) for z in x]
    ybar = np.sum(y)/len(y)          # or sum(y)/len(y)
    ssreg = np.sum((yhat-ybar)**2)   # or sum([ (yihat - ybar)**2 for yihat in yhat])
    sstot = np.sum((y - ybar)**2)    # or sum([ (yi - ybar)**2 for yi in y])
    results['R^2'] = ssreg / sstot

    parabcoeffs = np.polyfit(x, y, 5)
    poly1d_par = np.poly1d(parabcoeffs)

    # # search extremum
    # crit = poly1d_par.deriv().r
    # extremum_x = float(crit[crit.imag == 0].real)
    # extremum_y = poly1d_par(extremum_x)
    # extremum = [extremum_x, extremum_y]
    #
    # # test = poly1d_fn.deriv(2)(r_crit) # can check if min or max if necessary
    #
    # if abs(extremum[0] - 0.84) > abs(1.16 - extremum[0]):
    #     x2 = [0.84, 0.95, 1, 1.05, extremum[0]]
    #     y2 = [moment_ll, moment_l, 1, moment_h, extremum[1]]
    #     # new_slope, new_fit = linfit(x, y)
    #     lincoeffs2 = np.polyfit(x2, y2, 1)
    #     poly1d_lin2 = np.poly1d(lincoeffs2)
    # else:
    #     x2 = [extremum[0], 0.95, 1, 1.05, 1.16]
    #     y2 = [extremum[1], moment_l, 1, moment_h, moment_hh]
    #     # new_slope, new_fit = linfit(x, y)
    #     lincoeffs2 = np.polyfit(x2, y2, 1)
    #     poly1d_lin2 = np.poly1d(lincoeffs2)

    if plot == True:
        print('Plotting...')

        xp = np.linspace(0.8, 1.2, 50)

        fig, axs = plt.subplots(1, 3, figsize=(15, 6))
        axs[0].plot(x, y, 'yo', xp, poly1d_lin(xp), '--k')
        axs[1].plot(x, y, 'yo', xp, poly1d_par(xp), '--k')
        # axs[2].plot(x, y, 'yo', xp, poly1d_lin2(xp), '--c')

        # # plotting the fit line and origninal points
        # axs[0].plot(x, y, 'yo', xp, poly1d_fn(xp), '--k')
        #
        # # plotting x mark for extremum
        # x_min = extremum[0]  # [test > 0]
        # y_min = extremum[1]
        # axs[1].plot(x_min, y_min, 'xr')
        # axs[2].plot(x_min, y_min, 'xr')

        # plotting labels
        axs[0].set_xlabel('$\it{V_{def}/V}$', ha="left", fontsize=16)
        axs[0].set_ylabel('$\it{Moment}$', rotation=90, ha="right", fontsize=14)
        axs[1].set_xlabel('$\it{V_{def}/V}$', ha="left", fontsize=16)
        axs[1].set_ylabel('$\it{Moment}$', rotation=90, ha="right", fontsize=14)
        axs[2].set_xlabel('$\it{V_{def}/V}$', ha="left", fontsize=16)
        axs[2].set_ylabel('$\it{Moment}$', rotation=90, ha="right", fontsize=14)
        # set limits - needs modification
        xlimit = 0.8, 1.2
        ylimit = min(y)-0.3*min(y), max(y)+0.3*max(y)
        axs[0].set_ylim(ylimit)
        axs[0].set_xlim(xlimit)
        axs[1].set_ylim(ylimit)
        axs[1].set_xlim(xlimit)
        axs[2].set_ylim(ylimit)
        axs[2].set_xlim(xlimit)
        # show plot
        plt.show()

    return results

def linfit(x, y):
    coeffs = np.polyfit(x, y, 1)
    poly1d_fn = np.poly1d(coeffs)

    # fit quality
    yhat = poly1d_fn(x)  # or [p(z) for z in x]
    ybar = np.sum(y) / len(y)  # or sum(y)/len(y)
    ssreg = np.sum((yhat - ybar) ** 2)  # or sum([ (yihat - ybar)**2 for yihat in yhat])
    sstot = np.sum((y - ybar) ** 2)  # or sum([ (yi - ybar)**2 for yi in y])
    R_square = ssreg / sstot
    return coeffs[0], R_square


def parabfit(x, y):
    coeffs = np.polyfit(x, y, 2)
    poly1d_fn = np.poly1d(coeffs)
    crit = poly1d_fn.deriv().r
    extremum_x = float(crit[crit.imag == 0].real)
    extremum_y = poly1d_fn(extremum_x)
    extremum = [extremum_x, extremum_y]
    # test = poly1d_fn.deriv(2)(r_crit) # can check if min or max if necessary
    return extremum

def powerfit(x, y):
    coeffs = np.polyfit(x, y, 5)
    poly1d_fn = np.poly1d(coeffs)
    crit = poly1d_fn.deriv().r
    extremum_x = float(crit[crit.imag == 0].real)
    extremum_y = poly1d_fn(extremum_x)
    extremum = [extremum_x, extremum_y]
    # test = poly1d_fn.deriv(2)(r_crit) # can check if min or max if necessary
    return extremum

def magnetoelastic(ID, deformation):
    if deformation + '_inc' in os.listdir(vasp_results_dir + '/' + str(ID)):
        moment_inc, volume_inc, field_inc = moment_volume_after(ID, deformation+'_inc')
        moment_dec, volume_dec, field_dec = moment_volume_after(ID, deformation + '_dec')
        moment, volume, field = moment_volume_after(ID, 'undeformed')
        y = [field_dec, field, field_inc]
        x = [0.95, 1, 1.05]

        slope, fit_quality = linfit(x, y)
        new_slope, new_fit_quality = '', ''

        if fit_quality < 0.75:
            extremum = parabfit(x, y)

            if abs(extremum[0] - 0.95) > abs(1.05 - extremum[0]):
                x2 = [0.95, 1, extremum[0]]
                y2 = [field_dec, field, extremum[1]]
                new_slope, new_fit_quality = linfit(x2, y2)
            else:
                x2 = [extremum[0], 1, 1.05]
                y2 = [extremum[1], field, field_inc]
                new_slope, new_fit_quality = linfit(x2, y2)

        return slope, fit_quality, new_slope, new_fit_quality
    else:
        return '', '', '', ''


def magnetoelastic_numbers():

    field_dec = calculate_mag_field(12., 101.5132)
    field = calculate_mag_field(12.99, 105.7429)
    field_inc = calculate_mag_field(13.1342, 109.9726)

    # field_dec = calculate_mag_field(4.65, 89.5)
    # field = calculate_mag_field(12.99, 105.7429)
    # field_inc = calculate_mag_field(14.07, 123.33)

    x = [0.95, 1, 1.05]
    y = [field_dec, field, field_inc]

    slope, fit_quality = linfit(x, y)

    if fit_quality < 0.75:
        extremum = parabfit(x, y)

        if abs(extremum[0] - 0.95) > abs(1.05 - extremum[0]):
            x2 = [0.95, 1, extremum[0]]
            y2 = [field_dec, field, extremum[1]]
            new_slope, new_fit_quality = linfit(x2, y2)
        else:
            x2 = [extremum[0], 1, 1.05]
            y2 = [extremum[1], field, field_inc]
            new_slope, new_fit_quality = linfit(x2, y2)
    else:
        new_slope, new_fit_quality = slope, ''

    return slope, fit_quality, new_slope, new_fit_quality





# print(magnetoelastic_numbers())
# print(polyfit('mp-1222161', 'c'))
# print(polyfit_moment('mp-5951', 'V'))
# print(polyfit_moment('mp-582028', 'a'))

print(polyfit_moment_5point('mp-6670'))
# 2213
# 1200
# 238
