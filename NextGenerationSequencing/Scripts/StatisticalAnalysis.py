# Script for statistical analysis for NGS results 
# Author: Boris N. Schüpp

from DataUtilities import *
import scipy.stats as stats 
from PlottingNGS import *
from Fitting import *

data_directory = "../ProcessedData/704_AT_GC_400_1500_ATMi/"
graphics_directory = "../Graphics/704_AT_GC_400_1500_ATMi_Statistics/"
data_out_directory = "../ProcessedData/704_AT_GC_400_1500_ATMi_Statistics/"
data_file_mu = "Gaussian_Mu_Shifted_All.csv"
data_file_sigma = "Gaussian_Sigma_All.csv"

# Process data 

data_mu = read_fit_parameters([data_file_mu], data_directory)
data_sigma = read_fit_parameters([data_file_sigma], data_directory)

labels_ATMi = ["ATMi0", "ATMi1", "ATMi2", "ATMi3", "ATMi4"]
labels_704_400_1500_AT_GC = ["704", "400", "1500", "GC", "AT"]
data_mu_ATMi = {a:data_mu[a] for a in data_mu.keys() if a in labels_ATMi}
data_mu_704 = {a:data_mu[a] for a in data_mu.keys() if a in labels_704_400_1500_AT_GC}
data_sigma_ATMi = {a:data_sigma[a] for a in data_sigma.keys() if a in labels_ATMi}
data_sigma_704 = {a:data_sigma[a] for a in data_sigma.keys() if a in labels_704_400_1500_AT_GC}
print_multiple_dictionaries(data_mu_ATMi, "Sample number", "Gaussian_Mu_Shifted_ATMi.csv", data_out_directory)
print_multiple_dictionaries(data_mu_704, "Sample number", "Gaussian_Mu_Shifted_704_400_1500_AT_GC.csv", data_out_directory)
print_multiple_dictionaries(data_sigma_ATMi, "Sample number", "Gaussian_Sigma_ATMi.csv", data_out_directory)
print_multiple_dictionaries(data_sigma_704, "Sample number", "Gaussian_Sigma_704_400_1500_AT_GC.csv", data_out_directory)

# Plotting of fit parameters of Gaussian distributions 

plotter = PlottingNGS()

plotter.plot_fit_paramameter_scatter_single(["Gaussian_Mu_Shifted_ATMi.csv"], data_out_directory, "Mu", f"{graphics_directory}Mu_Shifted_Scatter_Plot_ATMi.pdf")
plotter.plot_fit_paramameter_scatter_single(["Gaussian_Mu_Shifted_704_400_1500_AT_GC.csv"], data_out_directory, "Mu",f"{graphics_directory}Sigma_Shifted_Scatter_Plot_Fig4.pdf")
plotter.plot_fit_paramameter_scatter_single_fit(["Gaussian_Sigma_ATMi.csv"], data_out_directory, "Sigma",f"{graphics_directory}Sigma_Scatter_Plot_ATMi.pdf")
plotter.plot_fit_paramameter_scatter_single(["Gaussian_Sigma_704_400_1500_AT_GC.csv"], data_out_directory, "Sigma",f"{graphics_directory}Mu_Shifted_Scatter_Plot_Fig4.pdf")

# One-way ANOVA with n = 3

for file, name in zip(["Gaussian_Mu_Shifted_ATMi.csv", "Gaussian_Mu_Shifted_704_400_1500_AT_GC.csv", "Gaussian_Sigma_ATMi.csv", "Gaussian_Sigma_704_400_1500_AT_GC.csv"], 
                      ["ANOVA_Exact_Mu_ATMi.csv", "ANOVA_Exact_Mu_704_AT_GC_400_1500.csv", "ANOVA_Exact_Sigma_704_AT_GC_400_1500.csv", "ANOVA_Exact_Sigma_ATMi.csv"]):

    data = read_fit_parameters([file], data_out_directory)

    out_matrix_exact = np.zeros((len(data.keys()), len(data.keys())))

    for a_idx, name_a in enumerate(data.keys()):
        for b_idx, name_b in enumerate(data.keys()):
            group_a = [data[name_a][i] for i in data[name_a].keys() if i != "Average"]
            group_b = [data[name_b][i] for i in data[name_b].keys() if i != "Average"]
            out_matrix_exact[a_idx][b_idx] = stats.f_oneway(group_a, group_b).pvalue

    DataUtilities.print_p_matrix(out_matrix_exact, list(data.keys()), name, data_out_directory)
    name_pdf = name.strip(".csv")+".pdf"
    plotter.plot_p_value_matrix(name, data_out_directory, f"{graphics_directory}{name_pdf}")

# Linear regression analysis 

x = [0, 1, 2, 3, 4]
y = [[data_mu_ATMi[sample][a] for a in data_mu_ATMi[sample] if a != "Average"] for sample in data_mu_ATMi.keys()]
        
result, prediction, prediction_with_confidence, x_big = fit_linear(x, y)
print("Linear regression result for ATMiX* Mu")
print(result.summary())
print(result.pvalues)

x = [0, 1, 2, 3, 4]
y = [[data_sigma_ATMi[sample][a] for a in data_sigma_ATMi[sample] if a != "Average"] for sample in data_sigma_ATMi.keys()]
        
result, prediction, prediction_with_confidence, x_big = fit_linear(x, y)
print("Linear regression result for ATMiX* Sigma")
print(result.summary())
print(result.pvalues)