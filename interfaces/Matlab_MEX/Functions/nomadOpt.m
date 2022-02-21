% NOMAD  Solve a Global MINLP/NLP using NOMAD, a Blackbox Optimization Library
%
%   min f(x)      subject to:     nlcon(x) <= 0
%    x                            lb <= x <= ub
%                                 x in R
%            
%%%%%%%%%%%%% GERAD VERSION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   [x,fval,exitflag,iter,nfval] = nomadOpt(fun,x0,lb,ub,params,cbFun)
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
%       x - solution vector
%       fval - objective value at the solution
%       exitflag - exit status (see below)
%       iter - number of iterations taken by the solver
%       nfval - number of function evaluations taken by the solver
%
%   Return Status:
%       1 - converged / target reached
%       0 - maximum iterations / function evaluations exceeded
%      -1 - infeasible / mesh limit reached
%      -2 - initialization error
%      -3 - nomad error
%      -5 - user exit
%
%
%   NOMAD is released under the Lesser GNU Public License (LGPL).
%   Type nomadOpt('-info') to see license and author details.
%
%   MEX Interface for GERAD version by C. Tribes
