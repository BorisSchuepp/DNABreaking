# Script for handling Fitting
# Author: Boris N. Schüpp

from scipy.optimize import curve_fit
import numpy as np
import statsmodels.api as sm 

# Defintion of a Gaussian function, takes a list of x values and a sigma and mu value and
# returns the respective y-values of the Gaussian Distribution
def gauss_function(x, sigma, mu):
    y = 1 / np.sqrt(2 * np.pi * sigma ** 2) * np.exp(-1 / 2 * ((x - mu) / sigma) ** 2)
    return y


def linear_function(x, m, b):
    y = m*x + b
    return y


def gauss_function_a(x, sigma, mu, a):
    y = a * np.exp(-1 / 2 * ((x - mu) / sigma) ** 2)
    return y


def bimodal(x, sigma1, mu1, a1, sigma2, mu2, a2):
    return gauss_function_a(x, sigma1, mu1, a1)+gauss_function_a(x, sigma2, mu2, a2)


def fit_bimodal(data_in, provide_data=False):
    data_x = list(data_in.keys())
    data_y = list(data_in.values())
    x_min = min(data_x)
    x_max = max(data_x)
    # width_guess = int((x_max - x_min) / 50)
    # p0 = [30, 353, 0.04, 25, 530, 0.02] # p0 for ssDNA
    p0 = [3, -10, 0.02, 3, 10, 0.04]  # p0 for 1500
    parameters_top = curve_fit(bimodal, data_x, data_y, p0=p0)[0]
    sigma_1 = parameters_top[0]
    mu_1 = parameters_top[1]
    a_1 = parameters_top[2]
    sigma_2 = parameters_top[3]
    mu_2 = parameters_top[4]
    a_2 = parameters_top[5]

    if provide_data:
        x_gauss = [x_min + i * 0.1 for i in range(0, int((x_max - x_min) * 10))]
        y_gauss = [bimodal(i, sigma_1, mu_1, a_1, sigma_2, mu_2, a_2) for i in x_gauss]
        return mu_1, sigma_1, a_1, mu_2, sigma_2, a_2, x_gauss, y_gauss
    else:
        return mu_1, sigma_1, a_1, mu_2, sigma_2, a_2


# Linear regression
def fit_linear(data_x, data_y):
    x_list = []
    y_list = []
    for x_val, y_vals in zip(data_x, data_y):
        for y_val in y_vals:
            x_list.append(x_val)
            y_list.append(y_val)
    x_with_const = sm.add_constant(x_list)
    model = sm.OLS(y_list, x_with_const)
    result = model.fit()
    x_big = np.linspace(min(x_list)-0.3, max(x_list)+0.3, 100)
    x_big_with_constant = sm.add_constant(x_big)
    prediction = result.get_prediction(x_big_with_constant)
    prediction_with_confidence = prediction.summary_frame(alpha=0.05)

    return result, prediction, prediction_with_confidence, x_big


# Performs a Gaussian fit on input data (dictionary x:y). Always returns sigma and mu, but can optionally return a
# thightly spaced list of y-values for the fitted Gaussian on the same x domain.
def fit_gauss(data_in, provide_data=False):
    data_x = list(data_in.keys())
    data_y = list(data_in.values())
    x_min = min(data_x)
    x_max = max(data_x)
    width_guess = int((x_max - x_min) / 10)
    p0 = [width_guess, x_min + (x_max - x_min) / 2]

    parameters_top = curve_fit(gauss_function, data_x, data_y, p0=p0)[0]
    sigma = parameters_top[0]
    mu = parameters_top[1]

    if provide_data:
        x_gauss = [x_min + i * 0.1 for i in range(0, int((x_max - x_min) * 10))]
        y_gauss = [gauss_function(i, sigma, mu) for i in x_gauss]
        return mu, sigma, x_gauss, y_gauss
    else:
        return mu, sigma
