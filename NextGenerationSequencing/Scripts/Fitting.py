# Script for handling Gaussian fits
# Author: Boris N. Schüpp

from scipy.optimize import curve_fit
import numpy as np


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
def fit_linear(data_x, data_y, interpolation):
    parameters_top = curve_fit(linear_function, data_x, data_y)[0]
    m = parameters_top[0]
    b = parameters_top[1]

    y_r2 = [linear_function(x, m, b) for x in data_x]
    residuals_square = sum([(y - y_pred)**2 for y, y_pred in zip(data_y, y_r2)])
    average = sum(data_y)/len(data_y)
    variance_square = sum([(y-average)**2 for y in data_y])

    r2 = 1-(residuals_square/variance_square)

    if interpolation:

        x_min = min(data_x)-0.1
        x_max = max(data_x)+0.1
        x_lin = [x_min + i * 0.1 for i in range(0, int((x_max - x_min) * 10))]
        y_lin = [linear_function(i, m, b) for i in x_lin]
        x_lin = [i + 1 for i in x_lin]
        return m, b, r2, x_lin, y_lin
    else:
        return m, b, r2


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
