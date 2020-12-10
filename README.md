# Extended-target-PMBM-and-PMB-filters
GGIW implementation of the extended target PMBM, MBM, PMB and MB filters.

This repository contains the MATLAB implementations of the following filters:
1. Track-oriented PMBM filter.
2. Track-oriented MBM filter.
3. Track-oriented PMB filter.
4. Track-oriented MB filter.
5. Variational PMB filter.
6. Variational MB filter.

Statistics and Machine Learning toolbox and Optimization toolbox are required.

References:

Karl Granström, Maryam Fatemi, and Lennart Svensson. "Poisson multi-Bernoulli mixture conjugate prior for multiple extended target filtering." IEEE Transactions on Aerospace and Electronic Systems 56.1 (2019): 208-225.

Xia, Yuxuan, Karl Granström, Lennart Svensson, Ángel F. García-Fernández, and Jason L. Williams. "Extended target Poisson multi-Bernoulli mixture trackers based on sets of trajectories." In 2019 22th International Conference on Information Fusion (FUSION), pp. 1-8. IEEE, 2019.

Xia, Yuxuan, Karl Granström, Lennart Svensson, and Maryam Fatemi. "Extended Target Poisson Multi-Bernoulli Filter." arXiv preprint arXiv:1801.01353 (2018).

The filters are evaluated using the generalised optimal subpattern-assignment (GOSPA) integrated with the Gaussian Wasserstein distance

Rahmathullah, Abu Sajana, Ángel F. García-Fernández, and Lennart Svensson. "Generalized optimal sub-pattern assignment metric." In 2017 20th International Conference on Information Fusion (Fusion), pp. 1-8. IEEE, 2017.

Video on GOSPA: https://www.youtube.com/watch?v=M79GTTytvCM

- main.m runs the extended target PMBM tracker.

- Data association algorithm: clustering (DBSCAN + k-means++) + assignment (Murty's algorithm, C++ implementation in the Tracker Component Library https://github.com/USNavalResearchLaboratory/TrackerComponentLibrary/blob/master/Assignment%20Algorithms/k-Best%202D%20Assignment/kBest2DAssign.cpp).

- The mex file kBest2DAssign.mexmaci64 only supports macOS.


