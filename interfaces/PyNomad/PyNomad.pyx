# PyNomad: Integration of Nomad in Python
# Based on work done by Christophe Tribes in NOMAD 3

#cython: language_level=3

from libcpp cimport bool
from libcpp.memory cimport shared_ptr, make_shared
from libcpp.string cimport string
from libcpp.vector cimport vector

from cython.operator cimport dereference as deref


def version():
    printPyNomadVersion()

# Define the interface function to display nomad general information
def usage():
    printPyNomadUsage()

# Define the interface function to display nomad general information
def info():
    printPyNomadInfo()

# Define the interface function to get nomad help
def help(about):
    about = about.encode(u"ascii")
    printNomadHelp(about)

def __doc__():
    cdef string about;
    printPyNomadUsage()
    help(about)

# Define the interface function to perform optimization (BATCH blackbox evaluation defined in parameter file)
def optimizeWithMainStep(params):
  cdef shared_ptr[AllParameters] allParameters_ptr = make_shared[AllParameters]()
  deref(allParameters_ptr).eraseAllEntries()

  cdef MainStep mainStep

  cdef size_t nbParams = len(params)
  for i in range(nbParams):
       encoded_parami = params[i].encode(u"ascii")
       # print(encoded_parami)
       deref(allParameters_ptr).readParamLine(encoded_parami)

  deref(allParameters_ptr).checkAndComply()
  mainStep.setAllParameters(allParameters_ptr)
  mainStep.start()
  mainStep.run()
  mainStep.end()

  # Make sure to reset static components for the next optimization
  mainStep.resetComponentsBetweenOptimization()

def suggest(params):
  cdef shared_ptr[AllParameters] allParameters_ptr = make_shared[AllParameters]()

 # The function eraseAllEntries must be called to reset a static variables
  deref(allParameters_ptr).resetToDefaultValues()
  deref(allParameters_ptr).eraseAllEntries()

  # The PyNomad::mainStep
  cdef MainStep mainStep

  cdef size_t nbParams = len(params)
  for i in range(nbParams):
    if type(params[i]) is str:
      encoded_parami= params[i].encode(u"ascii")
      # print(encoded_parami)
      deref(allParameters_ptr).readParamLine(encoded_parami)
      #allParameters.readParamLine(encoded_parami)

  # Get the rng state before checkAndComply
  rng_state = getRNGState()

  deref(allParameters_ptr).checkAndComply()
  mainStep.setAllParameters(allParameters_ptr)

  # Reset the cache (important for static members reinitialization) before calling suggest
  mainStep.resetCache()

  # Reset the evaluator control and its barrier (initialization from cache will create a new barrier)
  mainStep.resetEvaluatorControl()

  # Because the seed is not set as a Nomad parameters, the parameter has a 0 default value.
  # The check and comply will set the RNG state accordingly to a wrong state.
  # Reset to the rng state after checkAndComply
  setRNGState(rng_state)

  cdef vector[Point] xs = mainStep.suggest()

  candidates = []
  for i in range(xs.size()):
    candidates.append([xs[i][j].todouble() for j in range(xs[i].size())])

  return candidates

def setSeed(seed):
    cdef RNG rng
    rng.setSeed(seed)

def getRNGState():
    cdef RNG rng
    state = rng.getPrivateSeedAsString()
    if type(state) is bytes:
       encoded_state= state.decode('utf-8')
    else:
       encoded_state = state
    return encoded_state

def setRNGState(state):
    cdef RNG rng
    if type(state) is str:
       encoded_state= state.encode(u"ascii")
    else:
       encoded_state= state
    rng.setPrivateSeedAsString(encoded_state)


def resetRandomNumberGenerator():
    cdef RNG rng
    current_seed = rng.getSeed()
    rng.setSeed(current_seed) # Set seed performs a reset of the private seed


def observe(params,points,evals,udpatedCacheFileName):
    cdef shared_ptr[AllParameters] allParameters_ptr = make_shared[AllParameters]()
    # cdef shared_ptr[AllParameters] allParameters_ptr

    # allParameters_ptr.reset(new AllParameters())

    # The function eraseAllEntries must be called to reset a static variables
    deref(allParameters_ptr).resetToDefaultValues()
    deref(allParameters_ptr).eraseAllEntries()

    cdef MainStep mainStep

    cdef size_t nbParams = len(params)
    for i in range(nbParams):
       if type(params[i]) is str:
          encoded_parami= params[i].encode(u"ascii")
          # print(encoded_parami)
          deref(allParameters_ptr).readParamLine(encoded_parami)

    deref(allParameters_ptr).checkAndComply()
    mainStep.setAllParameters(allParameters_ptr)

    # Reset the cache (important for static members reinitialization) before calling observe
    mainStep.resetCache() 

    if len(points) != len(evals):
      print("Observe: Incompatible dimension of points and evaluations")
      return

    cdef vector[Point] vpoints
    cdef vector[ArrayOfDouble] vevals

    cdef PyNomadPoint nomadPoint
    cdef PyNomadArrayOfDouble nomadEval

    cdef vector[double] p
    cdef vector[double] v
    for i in range(len(points)):
      p = points[i][:]
      nomadPoint = PyNomadPoint(p)
      vpoints.push_back(nomadPoint.c_p)

      v = evals[i][:]
      nomadEval = PyNomadArrayOfDouble(v)
      vevals.push_back(nomadEval.c_aod)

    updatedParams = mainStep.observe(vpoints,vevals,udpatedCacheFileName.encode(u"ascii"))

    return updatedParams

# Define the interface function to perform optimization
# For now, we only show one best solution.
def optimize(fBB, pX0, pLB, pUB, params, fSurrogate=None):
    cdef PyNomadEvalPoint uFeas = PyNomadEvalPoint()
    cdef PyNomadEvalPoint uInfeas = PyNomadEvalPoint()
    cdef int runFlag = 0
    cdef size_t nbEvals = 0
    cdef size_t nbIters = 0
    cdef double fReturn = float("inf")
    cdef double hReturn = float("inf")
    xReturn = []
    eParams = []
    cdef string stopReason

    cdef size_t nbParams = len(params)
    for i in range(nbParams):
         eParams.append(params[i].encode(u"ascii"))

    if fSurrogate is None:
        runFlag = runNomad(cb, cbL, <void*> fBB, <vector[double]&> pX0,
                         <vector[double]&> pLB, <vector[double]&> pUB,
                         <vector[string]&> eParams,
                         uFeas.c_ep_ptr,
                         uInfeas.c_ep_ptr,
                         nbEvals, nbIters, stopReason)
    else:
        runFlag = runNomad(cb, cbL, <void*> fBB,  <void*> fSurrogate,
                         <vector[double]&> pX0,
                         <vector[double]&> pLB, <vector[double]&> pUB,
                         <vector[string]&> eParams,
                         uFeas.c_ep_ptr,
                         uInfeas.c_ep_ptr,
                         nbEvals, nbIters, stopReason)

    stopReasonU = stopReason.decode('utf-8')

    # For now, runNomad returns a feasible point or an infeasible point (least infeasible with smallest f) NOT BOTH
    if uFeas.c_ep_ptr != NULL:
        fReturn = uFeas.getF()
        hReturn = uFeas.getH()  # Should be 0
        for i in xrange(uFeas.size()):
            xReturn.append(uFeas.get_coord(i))

    if uInfeas.c_ep_ptr != NULL:
        fReturn = uInfeas.getF()
        hReturn = uInfeas.getH()
        for i in xrange(uInfeas.size()):
            xReturn.append(uInfeas.get_coord(i))

    return {'x_best': xReturn, 'f_best': fReturn, 'h_best': hReturn, 'nb_evals': nbEvals, 'nb_iters': nbIters, 'run_flag': runFlag, 'stop_reason': stopReasonU}

cdef extern from "Algos/MainStep.hpp" namespace "NOMAD":
    cdef cppclass MainStep:
        MainStep() except +
        void start()
        bool run()
        void end()
        void setAllParameters(const shared_ptr[AllParameters] & allParams )
        void setParamFileName(const string & paramfile)
        vector[Point] suggest()
        vector[string] observe(const vector[Point] & points, const vector[ArrayOfDouble] & evals, const string & updatedCacheFileName)
        @staticmethod
        void resetComponentsBetweenOptimization()
        @staticmethod
        void resetCache()
        @staticmethod
        void resetEvaluatorControl()

cdef extern from "Math/Point.hpp" namespace "NOMAD":
    cdef cppclass Point:
      Point() except+
      Point(const vector[double] &) except+
      const Double& operator[](size_t i) const
      size_t size()

cdef extern from "Math/ArrayOfDouble.hpp" namespace "NOMAD":
    cdef cppclass ArrayOfDouble:
      ArrayOfDouble() except+
      ArrayOfDouble(const vector[double] &) except+
      const Double& operator[](size_t i) const
      size_t size()

cdef extern from "Math/RNG.hpp" namespace "NOMAD":
    cdef cppclass RNG:
      void setSeed(int s)
      int getSeed()
      @staticmethod
      void resetPrivateSeedToDefault()
      @staticmethod
      string getPrivateSeedAsString()
      @staticmethod
      void setPrivateSeedAsString(const string & s)

cdef class PyNomadMainStep:
    cdef MainStep ms

    def __cinit__(self,params):
        self.ms = MainStep()

        cdef shared_ptr[AllParameters] allParameters_ptr = make_shared[AllParameters]()
        deref(allParameters_ptr).eraseAllEntries()

        cdef size_t nbParams = len(params)
        for i in range(nbParams):
          if type(params[i]) is str:
           encoded_parami= params[i].encode(u"ascii")
           # print(encoded_parami)
           deref(allParameters_ptr).readParamLine(encoded_parami)

        deref(allParameters_ptr).checkAndComply()
        self.ms.setAllParameters(allParameters_ptr)

    def suggest(self):
        cdef vector[Point] xs = self.ms.suggest()

        candidates = []
        for i in range(xs.size()):
            candidates.append([xs[i][j].todouble() for j in range(xs[i].size())])

        return candidates



cdef class PyNomadPoint:
    #cdef shared_ptr[Point] c_p_ptr
    cdef Point c_p

    def __cinit__(self):
        self.c_p = Point()

    def __cinit__(self, vector[double] & v):
        self.c_p = Point(v)

    def get_coord(self, size_t i):
      cdef PyNomadDouble coord = PyNomadDouble()
      # coord.c_d = deref(self.c_p_ptr)[i]
      coord.c_d = self.c_p[i]
      cdef double coord_d
      if (coord.isDefined()):
        coord_d = coord.todouble()
      else:
        coord_d = float("inf")
      return coord_d

    def size(self):
      cdef size_t n
      #n = deref(self.c_p_ptr).size()
      n = self.c_p.size()
      return n

cdef class PyNomadArrayOfDouble:
    cdef ArrayOfDouble c_aod

    def __cinit__(self):
      self.c_aod = ArrayOfDouble()

    def __cinit__(self, vector[double] & v):
      self.c_aod = ArrayOfDouble(v)

    def get_coord(self, size_t i):
      cdef PyNomadDouble coord = PyNomadDouble()
      coord.c_d = self.c_aod[i]
      cdef double coord_d
      if (coord.isDefined()):
        coord_d = coord.todouble()
      else:
        coord_d = float("inf")
      return coord_d

    def size(self):
      cdef size_t n
      n = self.c_aod.size()
      return n

cdef extern from "Param/AllParameters.hpp" namespace "NOMAD":
  cdef cppclass AllParameters:
        void read(const string & paramfile, bool overwrite , bool resetAllEntries )
        void readParamLine(const string & paramline)
        @staticmethod
        void eraseAllEntries()
        void checkAndComply()
        void resetToDefaultValues()

cdef extern from "Math/Double.hpp" namespace "NOMAD":
    cdef cppclass Double:
        Double() except+
        Double(const double &) except+
        const double & todouble()
        bool isDefined()

cdef class PyNomadDouble:
    cdef Double c_d
    def __cinit__(self, double v):
        self.c_d = Double(v)
    def __cinit__(self):
        pass
    def todouble(self):
        return self.c_d.todouble()
    def isDefined(self):
        return self.c_d.isDefined()

cdef extern from "Eval/EvalPoint.hpp" namespace "NOMAD":
    cdef cppclass EvalPoint:
        const Double& operator[](size_t i) const
        const Double& getF()
        const Double& getH()
        void setBBO(const string &bbo)
        string getBBO()
        size_t size()

    cdef cppclass Block:
        const shared_ptr[EvalPoint]& operator[](size_t i) const
        size_t size()


cdef class PyNomadEvalPoint:
    cdef shared_ptr[EvalPoint] c_ep_ptr

    def get_coord(self, size_t i):
        cdef PyNomadDouble coord = PyNomadDouble()
        coord.c_d = deref(self.c_ep_ptr)[i]
        cdef double coord_d
        if (coord.isDefined()):
            coord_d = coord.todouble()
        else:
            coord_d = float("inf")
        return coord_d


    def setBBO(self, string bbo):
        deref(self.c_ep_ptr).setBBO(bbo)

    def getF(self):
        cdef PyNomadDouble f = PyNomadDouble()
        f.c_d = deref(self.c_ep_ptr).getF()
        cdef double f_d
        if ( f.isDefined() ):
            f_d = f.todouble()
        else:
            f_d = float('inf')
        return f_d

    def getH(self):
        cdef PyNomadDouble h = PyNomadDouble()
        h.c_d = deref(self.c_ep_ptr).getH()
        cdef double h_d
        if ( h.isDefined() ):
            h_d = h.todouble()
        else:
            h_d = 0
        return h_d

    def size(self):
        cdef size_t n
        n = deref(self.c_ep_ptr).size()
        return n


cdef class PyNomadBlock:
    # Define what is needed to use blocks
    cdef shared_ptr[Block] c_block_ptr

    def size(self):
        cdef size_t n
        n = deref(self.c_block_ptr).size()
        return n

    def get_x(self, size_t i):
        cdef PyNomadEvalPoint x_i = PyNomadEvalPoint()
        x_i.c_ep_ptr = deref(self.c_block_ptr)[i]
        return x_i


cdef extern from "nomadCySimpleInterface.cpp":
    ctypedef int (*Callback)(void * apply, shared_ptr[EvalPoint] x)
    ctypedef vector[int] (*CallbackL)(void * apply, shared_ptr[Block] x)
    void printPyNomadVersion()
    void printPyNomadUsage()
    void printPyNomadInfo()
    void printNomadHelp(string about)
    int runNomad(Callback cb, CallbackL cbL, void* apply, vector[double] &X0,
                 vector[double] &LB, vector[double] &UB,
                 vector[string] &params,
                 shared_ptr[EvalPoint] &bestFeasSol,
                 shared_ptr[EvalPoint] &bestInfeasSol,
                 size_t &nbEvals, size_t &nbIters, string &stopReason) except+
    int runNomad(Callback cb, CallbackL cbL, void* applyBB, void* applySurrogate, vector[double] &X0,
                 vector[double] &LB, vector[double] &UB,
                 vector[string] &params,
                 shared_ptr[EvalPoint] &bestFeasSol,
                 shared_ptr[EvalPoint] &bestInfeasSol,
                 size_t &nbEvals, size_t &nbIters, string &stopReason) except+


# Define callback function for a single EvalPoint ---> link with Python
cdef int cb(void *f, shared_ptr[EvalPoint] x) noexcept:
    cdef PyNomadEvalPoint u = PyNomadEvalPoint()

    u.c_ep_ptr = x
    return (<object>f)(u)


# Define callback function for a block (vector) of EvalPoints
cdef vector[int] cbL(void *f, shared_ptr[Block] block) noexcept:

    cdef PyNomadBlock u = PyNomadBlock()

    u.c_block_ptr = block
    return (<object>f)(u)
 
