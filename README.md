# DRL
DRL (Dimension Reduction of Landscape) is an approach to reduce the dimensionality of a high-dimensional landscape (Xin Kang and Chunhe Li*,  A dimension reduction approach for energy landscape: Identifying intermediate states in EMT-metastasis network. Advanced Science, 8: 2003133 (2021).).
This is a matlab implementation of DRL approach on a mouse embryonic stem cell (MESC) network. This network includes 26 interactions (regulations) and 15 components (12 genes and 3 signals).


#####
MESC.m is the ODEs of the MESC model;

multivariate_normal_distribution is the density function of the Gaussian distribution which is used to calculate the density function of MESC model;

Solver.m is used to solve the ODEs, calculate the mean value and the covariance of each stable state and calculate the transition paths and actions between the stable states.

calculate_sigma.m is used to calculate the covariance of one stable state.

main.m is the main function which is uesd to calculate the density function of expression level of the system and use the DRL method to plot the landscape and transition paths of the system.
 

Please run main.m and the program runs about 3 minutes.
