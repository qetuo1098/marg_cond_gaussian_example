Instruction: run planar_pushing_matlab.m

This script illustrates our method of covariance extraction [1] in a constrained toy localization problem. It uses a miniature version of the box-pushing scenario used in our third experiment [1]. We took this scenario from InCOpt [2], a version of GTSAM [3] that accounts for hard constraints.

The scenario in this script is similar to the 2D planar pushing example used in [1,2], but is simplified and reimplemented in MATLAB.
The original experiment in [1] involved a number of modifications to [2], which was on top of an old version of the GTSAM repo [4]. To avoid instructions on installing heavy and old libraries, and to focus more on the covariance extraction section, we show a simpler version in this script.

[1] Z. C. Guo, J. R. Forbes, and T. D. Barfoot, "Marginalizing and Conditioning Gaussians onto Linear Approximations of Smooth Manifolds with Applications in Robotics," in IEEE ICRA, 2025.

[2] M. Qadri, P. Sodhi, J. G. Mangelson, F. Dellaert, and M. Kaess, "InCOpt: Incremental constrained optimization using the Bayes tree," in IEEE/RSJ IROS, 2022, pp. 6381–6388.

[3] M. Kaess, H. Johannsson, R. Roberts, V. Ila, J. Leonard, and F. Dellaert, “iSAM2: Incremental smoothing and mapping with fluid
relinearization and incremental variable reordering,” in IEEE ICRA, 2011, pp. 3281–3288.

[4] F. Dellaert and GTSAM Contributors, “borglab/gtsam,” May 2022, Georgia Tech Borg Lab. [Online]. Available: https://github.com/borglab/gtsam

