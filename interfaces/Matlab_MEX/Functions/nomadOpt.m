% NOMAD  Solve a Global MINLP/NLP using NOMAD, a Blackbox Optimization Library
%
%   min f(x)      subject to:     nlcon(x) <= 0
%    x                            lb <= x <= ub
%                                 x in R
%
%%%%%%%%%%%%% GERAD VERSION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   [x,fval,hinf,runflag,nfval] = nomadOpt(fun,x0,lb,ub,params,cbFun)
%
%   Input arguments:
%       fun - nonlinear function handle (f and nlcon)
%       x0 - initial solution guess.
%       lb - decision variable lower bounds. An empty lb ([]) can be provided.
%       ub - decision variable upper bounds. An empty ub ([]) can be provided.
%       params - parameters passed to the blackbox function as an array of strings using Nomad syntax (see user guide for the syntax).
%       cbFun - callback function handle (optional). When present, the function is called after calling fun. The user must provide a function that receives three input arguments: count_eval, fevals and x. Count_eval is the number of times fun has been called. Fevals and x are respectively the outputs returned by fun and its inputs. The function returns a boolean. If cbFun returns true, the optimization stops prematurely.
%
%   Return arguments:
%       x - Best point
%       fval - objective value at the best point
%       hinf - Infeasibility at best point
%       runflag - Nomad run termination flag  (see below)
%       nfval - number of function evaluations taken by the solver
%
%   Run flags:
%       1 - Objective target reached OR Mads converged (mesh criterion) to a feasible point (true problem).
%       0 - At least one feasible point obtained and evaluation budget (single bb or block of bb) spent or max iteration (user option) reached.
%      -1 - Mads mesh converged but no feasible point obtained (only infeasible) for the true problem.
%      -2 - No feasible point obtained (only infeasible) and evaluation budget (single bb or block of bb) spent or max iteration (user option) reached
%      -3 - Initial point failed to evaluate
%      -4 - Time limit reached (user option)
%      -5 - CTRL-C or user stopped (callback function)
%      -6 - Stop on feasible point (user option)
%
%
%   NOMAD is released under the Lesser GNU Public License (LGPL).
%   Type nomadOpt('-info') to see license and author details.
%
%   MEX Interface by C. Tribes
