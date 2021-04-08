Preface
=======

In many situations, one is interested in identifying the values of a set of variables that maximize or minimize some objective function. Furthermore, the variables cannot take arbitrary values, as they are confined to an admissible region and need to satisfy some prescribed requirements. NOMAD is a software application designed to solve these kind of problems.

The nature of the objective function and constraints dictates the type of optimization methods that should be used to tackle a given problem. If the optimization problem is convex, or if the functions are smooth and easy to evaluate, or if the number of variables is large, then NOMAD is not the solution that you should use. NOMAD is intended for time-consuming blackbox simulation with a small number of variables. NOMAD is often useful when other optimizers fail.

These nasty problems are called blackbox optimization problems. With NOMAD some con- straints may be evaluated prior to launching the simulation, and others may only be evaluated a posteriori. The simulations may take several seconds, minutes, hours or even days to compute. The blackbox can have limited precision and be contaminated with numerical noise. It may also fail to return a valid output, even when the input appears acceptable. Launching twice the simulation from the same input may produce different outputs. These unreliable properties are frequently encountered when dealing with real problems. The term blackbox is used to indicate that the internal structure of the target problem, such as derivatives or their approximations, cannot be exploited as it may be unknown, hidden, unreliable or inexistent. There are situations where some structure such as bounds may be exploited and in some cases, a surrogate of the problem may be supplied to NOMAD or a model may be constructed and trusted.

This guide describes how to use NOMAD to solve your blackbox optimization problem.
