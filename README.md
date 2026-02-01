# Supplementary-Matlab-Code
This is the code file for the paper "A class of nonparametric Max-EWMA schemes for monitoring high-dimensional processes".

Here's a brief description of each file:

File name: `bisection_method_to_refine_UCL.m`**
Description: Uses the bisection method to iteratively find the UCL value such that the in-control Median Run Length (MRL) equals a target value MRL0.

File name: `calculate_IC_MRL.m`**
 Description:  Estimates the in-control MRL for a given UCL by simulating r replications of the control chart and taking the median of the resulting run lengths.

File name: `calculate_OOC_MRL.m`**
 Description:  Similar to `calculate_IC_MRL.m`, but estimates the out-of-control MRL1 when a process shift (delta_1, delta_2) is present.

File name: `compute_conditional_IC_RL.m`**
 Description:  Simulates a single in-control run length for given parameters, using WRS and AB statistics based on distances from the origin, and an Max-EWMA-type updating rule.

File name: `compute_conditional_OOC_RL.m`**
 Description:  Same as above, but for the out-of-control case, where the process mean and/or variance are shifted.

File name: `main_function.m`**
 Description:  The main script that sets the parameters and calls the bisection method to find the UCL.

File name: `ME_OR.m`**  
Description:   Implements the ME-OR chart for monitoring the semiconductor manufacturing process. It allows different distance metrics (Chebyshev, Minkowski, etc.) and takes the maximum of two EWMA statistics for monitoring.

File name: `EWMA_Q.m`**  
Description:    Implements the EWMA-Q chart for monitoring the semiconductor manufacturing process. 

File name: `EWMA_S.m`**  
Description:   Implements the EWMA-S chart for monitoring the semiconductor manufacturing process. 

File name: `EWOR_L.m`**  
Description:   Implements the EWMA-L chart for monitoring the semiconductor manufacturing process. 

File name: `filtering_data.m`**
Description:  Preprocesses the SECOM dataset: removes constant columns, columns with too many missing values, and rows with any missing values. It then splits the data into Phase I (in-control) and Phase II (monitoring) sets and scales them.

File name: `Wilcoxon.m`**
Description:   Computes the WRS statistic for two samples based on ordered distances.

File name: `AB.m`**  
Description:  Computes the AB statistic for two samples based on ordered distances.

File name: `computeAEMMedian.m`**
Description:  Estimates the spatial median for the EWMA-S chart.
