What is NOMAD?
==============

NOMAD = Nonlinear Optimization by Mesh Adaptive Direct Search

NOMAD is a software application for simulation-based optimization. It can efficiently explore a design space in search of better solutions for a large spectrum of optimization problems.

Basics of the Mads algorithm
============================

At the core of NOMAD resides the Mesh Adaptive Direct Search (MADS) algorithm. As the name implies, this method generates iterates on a series of meshes with varying size. A mesh is a discretization of the space of variables. However, also as the name implies, the algorithm performs an adaptive search on the meshes including controlling the refinement of the meshes. The reader interested in the rather technical details should read [18].

The objective of each iteration of the MADS algorithm, is to generate a trial point on the mesh that improves the current best solution. When an iteration fails to achieve this, the next iteration is initiated on a finer mesh.

Each iteration is composed of two principal steps called the search and the poll steps [18]. The search step is crucial in practice because it is so flexible, but it is a difficulty for the theory for the same reason. search can return any point on the underlying mesh, but of course, it is trying to identify a point that improves the current best solution.

The poll step is more rigidly defined, though there is still some flexibility in how this is imple- mented. The poll step generates trial mesh points in the vicinity of the best current solution. Since the poll step is the basis of the convergence analysis, it is the part of the algorithm where most research has been concentrated.

Using NOMAD
===========

NOMAD does not provide a graphical user interface to define and perform optimization. Mini- mally, users must accomplish several tasks to solve their own optimization problems:

* Create a custom blackbox program(s) to evaluate the functions f and cj OR embed the functions evaluations in C++ source code to be linked with the NOMAD library.

* Create the optimization problem definition in a parameter file OR embed the problem definition in C++ source code to be linked with the NOMAD library.

* Launch the execution at the command prompt OR from another executable system call.


Users can find several examples provided in the installation package and described in this user guide to perform customization for their problems. The installation procedure is given in Chapter 2. New users should refer to Chapter 3 to get started. The most important instructions to use NOMAD are in Chapter 4. In addition, tricks that may help solving specific problems and improve NOMAD efficiency are presented in Chapter 5. Advanced parameters and functionalities are presented in Chapters 6 and 7.


License
=======

NOMAD is a free software application released under the GNU Lesser General Public License v 3.0. As a free software application you can redistribute and/or modify NOMAD source codes under the terms of the GNU Lesser General Public License.

For more information, please refer to the local copy of the license obtained during installation. For additional information you can contact us or visit the Free Software Foundation website.


Contact us
==========

All queries can be submitted by email at nomad@gerad.ca. In particular, feel free to ask technical support for problem specification (creating parameter files or integration with various types of simulations) and system support (installation and plateform-dependent problems).

Bug reports and suggestions are valuable to us! We are committed to answer to posted requests as quickly as possible.


Supported plateforms and environments
=====================================

NOMAD source codes are in C++ and are identical for all supported platforms.

Authors and fundings
====================

The development of NOMAD started in 2001, and was funded in part by .....


Developers of the methods behind NOMAD include:

* Mark A. Abramson (abramson@mathematics.byu.edu), Bringham Young University.
* Charles Audet (www.gerad.ca/Charles.Audet), GERAD and Département de mathéma- tiques et de génie industriel, École Polytechnique de Montréal.
* J.E. Dennis Jr. (www.caam.rice.edu/∼dennis), Computational and Applied Mathematics Department, Rice University.
* Sébastien Le Digabel (www.gerad.ca/Sebastien.Le.Digabel), GERAD and Département de mathématiques et de génie industriel, École Polytechnique de Montréal.
* Viviane Rochon Montplaisir, www.gerad.ca and Département de mathématiques et de génie industriel, École Polytechnique de Montréal.
* Christophe Tribes, www.gerad.ca and Département de mathématiques et de génie indus- triel, École Polytechnique de Montréal.

The library for dynamic surrogates (SGTELIB) has been developed by Bastien Talgorn (bastien- talgorn@fastmail.com), McGill University, Montreal. The SGTELIB is included in NOMAD since version 3.8.0.

Version 3.5.1 (and above) of NOMAD is developed by Viviane Rochon Montplaisir and Christophe Tribes. Version 3.0 (and above) was developed by Sébastien Le Digabel. Previous versions were written by Gilles Couture.

Acknowledgments
===============

The developers of NOMAD wish to thank Florian Chambon, Mohamed Sylla and Quentin Reynaud, all from ISIMA, for their contribution to the project during Summer internships, and to Anthony Guillou and Dominique Orban for their help with AMPL, and their suggestions.

A special thank to Maud Bay, Eve Bélisle, Vincent Garnier, Michal Kvasnička, Alexander Lutz, Rosa-Maria Torres-Calderon, Yuri Vilmanis, Martin Posch, Etienne Duclos, Emmanuel Bigeon, Walid Zghal, Jerawan Armstrong, Stéphane Alarie and Klaus Truemper for their feedbacks and tests that significantly contributed to improve NOMAD. Some features of NOMAD have been developed under the impulsion of enthusiastic users/developers: Andrea Ianni, Florian Chambon, Mohamed Sylla, Quentin Reynaud, Amina Ihaddadene, Bastien Talgorn, Nadir Amaioua and Catherine Poissant. We also wish to thank Pascal Côté for his contribution in the development of the Python interface pyNomad and Jonathan Currie for the development of the foundations for a strong NOMAD interface for MATLAB.

Finally, many thanks to the TOMS anonymous referees for their useful comments which helped a lot to improve the code and the text of [50].

