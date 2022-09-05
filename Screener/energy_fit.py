from pymatgen.analysis.eos import BirchMurnaghan
import numpy as np


def make_eosdict(volumes: np.ndarray,
                 energies: np.ndarray,
                 volumes_fitted: np.ndarray = None,
                 include_eos_object=True):
    '''
    Takes 1d numpy arrays of volumes, corresponding energies, and additional
    volumes to pass through the fitted Birch Murnaghan function.
    Output is a dictionary of function parameters, including the adjusted Rsquared
    and the 10log of the adjusted FVU, as well as the fitted volumes and energies.
    To surpress fitting values, volumes_fitted should be an empty array.
    '''

    volumes = np.asarray(volumes)
    energies = np.asarray(energies)
    # default fitting volumes
    if volumes_fitted is None:
        volumes_fitted = np.linspace(np.min(volumes), np.max(volumes), 101)

    N = np.size(volumes)  # number of samples
    K = 3  # Number of regressors

    # The fit itself is a simple two steps using pymatgen
    eos = BirchMurnaghan(volumes, energies)
    eos.fit()

    # Get Rsquared_adjusted using the fitting function eos.func()
    # to get predicted energies at input volumes
    energies_predicted = eos.func(volumes)
    Rsquared_adjusted = adjusted_Rsquared(energies_predicted, energies, K)
    FVU_adjusted = 1 - Rsquared_adjusted
    FVU_adj_10log = np.log10(FVU_adjusted)

    if np.any(volumes_fitted):
        energies_fitted = eos.func(volumes_fitted)
    else:
        energies_fitted = np.empty()

    eosdict = eos.results
    eosdict['Rsquared_adjusted'] = Rsquared_adjusted
    eosdict['FVU_adj_10log'] = FVU_adj_10log
    eosdict['volumes'] = volumes
    eosdict['energies'] = energies
    eosdict['volumes_fitted'] = volumes_fitted
    eosdict['energies_fitted'] = energies_fitted
    if include_eos_object:
        eosdict['eos_object'] = eos

    return Rsquared_adjusted


def adjusted_Rsquared(predictions: np.ndarray,
                      measurements: np.ndarray,
                      regressors: int):
    N = np.size(measurements)
    if N != np.size(predictions):
        print("ERROR: unequal number of predictions ans measurements")

    SSR = 0  # Residual Sum of Squares
    SST = 0  # Total Sum of Squares
    avg_measurement = np.mean(measurements)
    for i in range(N):
        SSR += (predictions[i] - measurements[i]) ** 2
        SST += (measurements[i] - avg_measurement) ** 2
    FVU = SSR / SST  # Fraction of Variance Unexplained
    Rsquared = 1 - FVU

    # adjusting these metrics for degrees of freedom
    # Note that the adjusted Rsquared can be <0.
    # Also, N-regressors-1 == 0 will result in a divide by zero error if treated normally
    if N - regressors - 1 >= 1:
        Rsquared_adjusted = 1 - (1 - Rsquared) * (N - 1) / (N - regressors - 1)
    else:
        Rsquared_adjusted = 0
    # FVU_adjusted = 1 - Rsquared_adjusted

    return Rsquared_adjusted



#volumes_test = [141.583, 144.532, 147.482, 150.432, 153.381]
#energies_test = [-83.1553, -83.3324, -83.3857, -83.3283, -83.1711]

#eosdict_test = make_eosdict(volumes_test, energies_test)
#print(eosdict_test)