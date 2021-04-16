Release notes and future developments
=====================================

NOMAD 4 is a complete redesign compared with NOMAD, with a new architecture providing more flexible code, some added functionalities and reusable code.

Some functionalities available in NOMAD 3 will be included in NOMAD 4 in future releases:

* Use static surrogate evaluations [AuCM2019]_
* *BiMads* [AuSaZg2008a]_
* *RobustMads* [AudIhaLedTrib2016]_ and *StoMads* [G-2019-30]_
* *Variable Neighborood Search* [AuBeLe08]_
* Categorical [AuDe01a]_ and periodical variables [AuLe2012]_

The performance of NOMAD 4 and 3 are similar when the default parameters of NOMAD 4 are used (see [AuLeRoTr2021]_).

.. topic:: References

  .. [AuSaZg2008a] C. Audet, G. Savard, and W. Zghal. 2008.  Multiobjective Optimization Through a Series of Single-Objective Formulations. *SIAM Journal onOptimization* 19, 1 (2008), 188–210
  .. [AuBeLe08] C. Audet, V. Béchard, and S. Le Digabel. 2008. Nonsmooth optimization through Mesh Adaptive Direct Search and Variable Neighborhood Search. *Journal of Global Optimization* 41, 2 (2008), 299–318
  .. [AudIhaLedTrib2016] C. Audet, A. Ihaddadene, S. Le Digabel, and C. Tribes. 2018. Robust optimization of noisy blackbox problems using the Mesh Adaptive Direct Search algorithm. *Optimization Letters* 12, 4 (2018), 675–689
  .. [G-2019-30] C. Audet, K.J. Dzahini, M. Kokkolaras, and S. Le Digabel. 2021.Stochastic mesh adaptive direct search for blackbox optimization using probabilistic estimates. *Technical Report* G-2019-30. Les cahiers du GERAD.  To appear in *Computational Optimization and Applications*.
  .. [AuDe01a] C. Audet and J.E. Dennis, Jr. 2001. Pattern Search Algorithms for Mixed Variable Programming. *SIAM Journal on Optimization* 11, 3 (2001), 573–594.
  .. [AuLe2012] C. Audet and S. Le Digabel. 2012.  The mesh adaptive direct search algorithm for periodic variables. *Pacific Journal of Optimization* 8, 1 (2012),103–119
  .. [AuCM2019] C. Audet and J. Côté-Massicotte. Dynamic improvements of static surrogates in direct search optimization. *Optimization Letters* 13, 6 (2019), 1433-1447
