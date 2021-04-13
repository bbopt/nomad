# PyNomad: Integration of Nomad in Python
# Based on work done by Christophe Tribes in NOMAD 3

#cython: language_level=3

from libcpp cimport bool
from libcpp.memory cimport shared_ptr
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
def help(about=''):
    about = about.encode(u"ascii")
    printNomadHelp(about)

def __doc__():
    cdef string about;
    printPyNomadUsage()
    help(about)


# Define the interface function to perform optimization          
# TODO: Show multiple best solutions, and show both feas and infeas solutions.
# For now, we only show one best solution.
def optimize(f, pX0, pLB, pUB, params):
    cdef PyNomadEvalPoint uFeas = PyNomadEvalPoint()
    cdef PyNomadEvalPoint uInfeas = PyNomadEvalPoint()
    cdef int runStatus = 0
    cdef size_t nbEvals = 0
    cdef size_t nbIters = 0
    cdef double fReturn = float("inf")
    cdef double hReturn = float("inf")
    xReturn = []
    
    cdef size_t nbParams = len(params)
    for i in range(nbParams):
         params[i] = params[i].encode(u"ascii")

    runStatus = runNomad(cb, cbL, <void*> f, <vector[double]&> pX0,
                         <vector[double]&> pLB, <vector[double]&> pUB,
                         <vector[string]&> params,
                         uFeas.c_ep_ptr,
                         uInfeas.c_ep_ptr,
                         nbEvals, nbIters)
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

    return [ xReturn, fReturn, hReturn, nbEvals, nbIters, runStatus ]


cdef extern from "Math/Double.hpp" namespace "NOMAD":
    cdef cppclass Double:
        const double & todouble()
        bool isDefined()

cdef class PyNomadDouble:
    cdef Double c_d
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
    ctypedef int (*Callback)(void * apply, shared_ptr[EvalPoint] x, bool hasSgte, bool sgteEval)
    ctypedef vector[int] (*CallbackL)(void * apply, shared_ptr[Block] x, bool hasSgte, bool sgteEval)
    void printPyNomadVersion()
    void printPyNomadUsage()
    void printPyNomadInfo()
    void printNomadHelp(string about)
    int runNomad(Callback cb, CallbackL cbL, void* apply, vector[double] &X0,
                 vector[double] &LB, vector[double] &UB,
                 vector[string] &params,
                 shared_ptr[EvalPoint] &bestFeasSol,
                 shared_ptr[EvalPoint] &bestInfeasSol,
                 size_t &nbEvals, size_t &nbIters) except+


# Define callback function for a single EvalPoint ---> link with Python     
cdef int cb(void *f, shared_ptr[EvalPoint] x, bool hasSgte, bool sgteEval):
    cdef PyNomadEvalPoint u = PyNomadEvalPoint()

    u.c_ep_ptr = x
    return (<object>f)(u)


# Define callback function for a block (vector) of EvalPoints
cdef vector[int] cbL(void *f, shared_ptr[Block] block, bool hasSgte, bool sgteEval):

    cdef PyNomadBlock u = PyNomadBlock()

    u.c_block_ptr = block
    return (<object>f)(u)

