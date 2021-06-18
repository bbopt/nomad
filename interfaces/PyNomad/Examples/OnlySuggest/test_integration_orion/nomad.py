# -*- coding: utf-8 -*-
"""
:mod:`orion.algo.nomad.nomad -- TODO
=============================================

.. module:: nomad
    :platform: Unix
    :synopsis: TODO

TODO: Write long description
"""
import copy

import numpy 
import os

from orion.algo.base import BaseAlgorithm
from orion.core.utils.points import flatten_dims, regroup_dims

import PyNomad


class nomad(BaseAlgorithm):
    """Nomad is a Mesh Adaptive Direct Search (MADS) algorithm for blackbox optimization.

    For more information about MADS and NOMAD: www.gerad.ca/nomad

    Parameters
    ----------
    space: `orion.algo.space.Space`
        Optimisation space with priors for each dimension.
    seed: int
        Seed of Nomad random number generator.
        Default: 0
    mega_search_poll: bool
        Use Mads mega search poll strategy to generate Points
        Default: True
    lh_eval_n_factor: into
        Multiply the factor by n and obtain the number of latin hypercube
        samples used in the initial phase


    """

    requires_dist = "linear"
    requires_shape = "flattened"
    requires_type = None

    # Global flag to use LH_EVAL or MEGA_SEARCH_POLL
    use_initial_params = True

    # Global flag to stop when no points are suggestes
    no_candidates_suggested = False

    def __init__(self, space, seed=None, mega_search_poll=True, lh_eval_n_factor=3):
        super(nomad, self).__init__(space,seed=seed,
                                          mega_search_poll=mega_search_poll,
                                          lh_eval_n_factor=lh_eval_n_factor)


    @property
    def space(self):
        """Return transformed space of PyNomad"""
        return self._space

    @space.setter
    def space(self, space):
        """Set the space of PyNomad and initialize it"""
        self._original = self._space
        self._space = space
        self._initialize(space)

    def _initialize(self, space):

        assert self.mega_search_poll, "For the moment PyNomad only works with mega_search_poll activated"
        assert self.lh_eval_n_factor > 0, "For the moment PyNomad only works with lh_eval_n_factor>0"

        # For sampled point ids
        self.sampled = set()

        #
        # Create Nomad parameters
        #

        # Dimension, bounds and bb_input_type  for flattened space
        dim = 0
        dimension_string = 'DIMENSION '
        lb_string = 'LOWER_BOUND ( '
        ub_string = 'UPPER_BOUND ( '
        all_variables_are_granular = True
        bb_input_type_string = 'BB_INPUT_TYPE ( '
        for val in self.space.values():
            if val.type == "fidelity" :
                raise ValueError( "PyNomad do not support fidelity type" )

            if val.prior_name not in [
                "uniform",
                "reciprocal",
                "int_uniform",
                "int_reciprocal",
            ]:
                raise ValueError(
                    "PyNomad now only supports uniform, loguniform, uniform discrete"
                    f" as prior: {val.prior_name}"
                )

            shape = val.shape
            
            
            if shape and len(shape) != 1:
                raise ValueError("Nomad now only supports 1D shape.")
            elif len(shape) == 0 :
                dim += 1
                lb_string += str(val.interval()[0]) + ' '
                ub_string += str(val.interval()[1]) + ' '
                if val.type == 'real':
                    bb_input_type_string += 'R '
                    all_variables_are_granular = False
                elif val.type == 'integer' :
                    bb_input_type_string += 'I '
                else :
                    raise ValueError("PyNomad now only real and integer type ")
            else :
                dim += shape[0]
                for s in range(shape[0]):
                    lb_string += str(val.interval()[0]) + ' '
                    ub_string += str(val.interval()[1]) + ' '
                    if val.type == 'real':
                        bb_input_type_string += 'R '
                        all_variables_are_granular = False
                    elif val.type == 'integer' :
                        bb_input_type_string += 'I '
                    else :
                        raise ValueError("PyNomad now only real and integer type ")

        dimension_string += str(dim)
        lb_string += ' )'
        ub_string += ' )'
        bb_input_type_string += ' )'

        # Todo constraints
        bbo_type_string = 'BB_OUTPUT_TYPE OBJ '

        self.cache_file_name = 'cache.txt'
        if os.path.exists(self.cache_file_name):
            os.remove(self.cache_file_name)
            self.use_initial_params = True
        cache_file_string = 'CACHE_FILE '+self.cache_file_name

        suggest_algo = 'MEGA_SEARCH_POLL yes' # Mads MegaSearchPoll for suggest after first suggest
        first_suggest_algo = 'LH_EVAL ' + str(len(self.space.values())*2)

        # IMPORTANT
        # Seed is managed explicitely with PyNomad.setSeed. Do not pass SEED as a parameter

        # if all_variables_are_granular:
        #    self.max_calls_to_extra_suggest = 100 * pow(3,len(self.space.values())) # This value is MAX_EVAL in Nomad when all variables are granular
        #else :
        self.max_calls_to_extra_suggest = 10 # This is arbitrary

        self.initial_params = ['DISPLAY_DEGREE 2', dimension_string, bb_input_type_string, bbo_type_string,lb_string, ub_string, cache_file_string, first_suggest_algo ]
        self.params = ['DISPLAY_DEGREE 2', dimension_string, bb_input_type_string, bbo_type_string,lb_string, ub_string, cache_file_string, suggest_algo ]

        # print(self.initial_params, self.params)

        # list to keep candidates for an evaluation
        self.stored_candidates = list()


    def seed_rng(self, seed):
        """Seed the state of the random number generator.

        :param seed: Integer seed for the random number generator.

        .. note:: This methods does nothing if the algorithm is deterministic.
        """
        self.seed = seed

        PyNomad.setSeed(seed)
        self.rng_state = PyNomad.getRNGState()

        # print("Seed rng: ", seed,self.rng_state)


    @property
    def state_dict(self):
        """Return a state dict that can be used to reset the state of the algorithm."""

        self.rng_state = PyNomad.getRNGState()
        # print("State dict : ",self.rng_state)
        # return {'rng_state': self.rng_state, 'use_initial_params': self.use_initial_params, "_trials_info": copy.deepcopy(self._trials_info)}

        return {'rng_state': self.rng_state, "_trials_info": copy.deepcopy(self._trials_info)}

    def set_state(self, state_dict):
        """Reset the state of the algorithm based on the given state_dict

        :param state_dict: Dictionary representing state of an algorithm
        """

        self.rng_state = state_dict["rng_state"]
        self._trials_info = state_dict.get("_trials_info")
        # self.use_initial_params =state_dict.get("use_initial_params")

        # print("Set state : ",state_dict)
        PyNomad.setRNGState(self.rng_state)

    def suggest(self, num=None):
        """Suggest a `num`ber of new sets of parameters.

        TODO: document how suggest work for this algo

        Parameters
        ----------
        num: int, optional
            Number of points to suggest. Defaults to None.

        Returns
        -------
        list of points or None
            A list of lists representing points suggested by the algorithm. The algorithm may opt
            out if it cannot make a good suggestion at the moment (it may be waiting for other
            trials to complete), in which case it will return None.

        Notes
        -----
        New parameters must be compliant with the problem's domain `orion.algo.space.Space`.

        """

        # Clear candidates before a new suggestion
        self.stored_candidates.clear()


        # print('Use initial params : ' , self.use_initial_params)
        print('Suggest RNG State: ',self.rng_state)

        if self.use_initial_params:
            print("Initial Params for suggest:",self.initial_params)
            self.stored_candidates = PyNomad.suggest(self.initial_params)
        else:
            print("Params for suggest:", self.params)
            self.stored_candidates = PyNomad.suggest(self.params)

        #print("Suggest: ",self.stored_candidates)

        # extra suggest with LH to force suggest of candidates
        nb_suggest_tries = 0
        while len(self.stored_candidates) < num and nb_suggest_tries < self.max_calls_to_extra_suggest:
            self.stored_candidates.extend(PyNomad.suggest(self.initial_params))
            nb_suggest_tries += 1
            #print("Extra Suggest (LH): ",self.stored_candidates)

        # assert len(self.stored_candidates) > 0, "At least one candidate must be provided !"

        # Todo manage prior conversion : candidates -> samples
        samples = []
        for point in self.stored_candidates:
            point = regroup_dims(point, self.space)
            self.register(point)
            samples.append(point)
            if len(samples) >= num:   # return the number requested.
                break;

        num = len(samples)
        self.no_candidates_suggested = (num == 0 )

        print("Suggest samples: ",samples)

        if samples:
           return samples

        return None

    def observe(self, points, results):
        """Observe evaluation `results` corresponding to list of `points` in
        space.

        Feed an observation back to PyNomad.

        Observe puts points and corresponding results in Nomad cache file. Observe updates the mesh and frame size.
        The updated cache file and frame size are used by next suggest.


        Parameters
        ----------
        points : list of tuples of array-likes
           Points from a `orion.algo.space.Space`.
           Evaluated problem parameters by a consumer.
        results : list of dicts
           Contains the result of an evaluation; partial information about the
           black-box function at each point in `params`.

        Result
        ------
        objective : numeric
           Evaluation of this problem's objective function.
        constraint : list of numeric, optional
           List of constraints expression evaluation which must be greater
           or equal to zero by the problem's definition.

        """
        assert len(points) == len(results), "The length of results and points are not the same"

        #print('observe, use_first_params=',self.use_initial_params)
        candidates_outputs = list()
        candidates = list()
        for point, result in zip(points, results):

            #if not self.has_suggested(point):
            #    logger.info(
            #        "Ignoring point %s because it was not sampled by current algo.",
            #        point,
            #    )
            #    continue
            tmp_outputs = list()
            tmp_outputs.append(result['objective']) # TODO constraints

            candidates_outputs.append(tmp_outputs) # TODO constraints
            flat_point = flatten_dims(point,self.space)
            flat_point_tuple = list()
            for x in flat_point:
                 #print(type(x))
                 if type(x)==numpy.ndarray:
                    assert x.size==1, "The length of the ndarray should be one"
                    flat_point_tuple.append(numpy.float64(x))
                 else:
                    flat_point_tuple.append(x)
            candidates.append(flat_point_tuple)

        print("Call PyNomad observe")
        print(candidates_outputs)
        print(candidates)

        if self.use_initial_params:
             # print("Initial params:",self.initial_params)
             updatedParams = PyNomad.observe(self.initial_params,candidates,candidates_outputs,self.cache_file_name)
             self.use_initial_params = False  # after initial observe we use only params
        else:
             updatedParams = PyNomad.observe(self.params,candidates,candidates_outputs,self.cache_file_name)


        super(nomad, self).observe(points, results) 


    	# Decode bytes into string
        for i in range(len(updatedParams)):
            updatedParams[i] = updatedParams[i].decode('utf-8')
        for i in range(len(self.params)):
            if type(self.params[i]) is bytes:
                self.params[i] = self.params[i].decode('utf-8')

        # print("Updated parameters by observe:\n",updatedParams)

    	# Replace updated params in params OR add if not present
        for i in range(len(updatedParams)):
            split1 = updatedParams[i].split()
            found = False
            for j in range(len(self.params)):
                split2 = self.params[j].split()
                if ( split2[0].upper() == split1[0].upper() ):
                    self.params[j] = updatedParams[i]
                    found = True
                    break;
            if not found:
                self.params.append(updatedParams[i])

        #print("Parameters for next iteration:\n",self.params)
        #print("\n")

    @property
    def is_done(self):
        """Return True, if an algorithm holds that there can be no further improvement."""
        # NOTE: Drop if base implementation is fine.
        return self.no_candidates_suggested or super(nomad, self).is_done

    def score(self, point):
        """Allow algorithm to evaluate `point` based on a prediction about
        this parameter set's performance.
        """
        # NOTE: Drop if not used by algorithm
        pass

    def judge(self, point, measurements):
        """Inform an algorithm about online `measurements` of a running trial."""
        # NOTE: Drop if not used by algorithm
        pass

    @property
    def should_suspend(self):
        """Allow algorithm to decide whether a particular running trial is still
        worth to complete its evaluation, based on information provided by the
        `judge` method.

        """
        # NOTE: Drop if not used by algorithm
        pass
