from libcpp cimport bool
from libcpp.memory cimport shared_ptr
from libcpp.string cimport string
from libcpp.vector cimport vector

# ============================================================================
# C++ Type Declarations (cdef extern from)
# ============================================================================

cdef extern from "Math/Double.hpp" namespace "NOMAD":
    cdef cppclass Double:
        Double() except+
        Double(const double &) except+
        const double & todouble()
        bool isDefined()

cdef extern from "Type/EvalType.hpp" namespace "NOMAD":
    cpdef enum class EvalType:
        BB,
        SURROGATE,
        MODEL,
        LAST,
        UNDEFINED

cdef extern from "Type/ComputeType.hpp" namespace "NOMAD":
    cpdef enum class ComputeType:
        STANDARD
        PHASE_ONE
        USER
        UNDEFINED

    cpdef enum class HNormType:
        L1
        L2
        Linf

    cdef struct FHComputeTypeS:
            pass

    cdef struct FHComputeType:
        pass


cdef extern from "Eval/EvalPoint.hpp" namespace "NOMAD":
    cdef cppclass EvalPoint:
        const Double& operator[](size_t i) const
        const Double& getF(const FHComputeType& completeComputeType) const
        const Double& getH(const FHComputeType& completeComputeType) const
        void setBBO(const string &bbo)
        string getBBO(const EvalType& evalType) const
        size_t size()
        string display()

    cdef cppclass Block:
        const shared_ptr[EvalPoint]& operator[](size_t i) const
        size_t size()

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


cdef extern from "nomadCySimpleInterface.cpp":
    ctypedef int (*Callback)(void * apply, shared_ptr[EvalPoint] x)
    ctypedef vector[int] (*CallbackL)(void * apply, shared_ptr[Block] x)

    void printNomadHelp(string about)

    ctypedef vector[vector[double]] (*CustomOrderCompPYXFunc)(vector[vector[double]]& v) except *
    ctypedef bool (*ModelCacheUpdatePYXFunc)(vector[double]& v) except *
    ctypedef bool (*ModelBuildPYXFunc)() except *
    ctypedef vector[vector[double]] (*CatMadsFreePollPYXFunc)(vector[double]& v) except *
    ctypedef vector[vector[double]] (*UserSearchPYXFunc)(vector[double]& vFeas, vector[double]& vInfeas) except *

    int runNomad(
        Callback cb,
        CallbackL cbL,
        void* apply,
        CustomOrderCompPYXFunc comp,
        ModelCacheUpdatePYXFunc cacheUpdate,
        ModelBuildPYXFunc modelBuildPYXFunc,
        CatMadsFreePollPYXFunc catMadsFreePollPYXFunc,
        UserSearchPYXFunc userSearchPYXFunc,
        vector[string] &params,
        shared_ptr[Block] &bestFeasSol,
        shared_ptr[Block] &bestInfeasSol,
        size_t &nbEvals,
        size_t &nbIters,
        string &stopReason
    ) except+

# ============================================================================
# Cython Extension Classes (cdef class)
# ============================================================================

cdef class NomadBBOEvalPoint:
    cdef shared_ptr[EvalPoint] c_ep_ptr
    cdef void _require_ptr(self)

cdef class NomadBBOBlock:
    cdef shared_ptr[Block] c_block_ptr

cdef class NomadBBOPoint:
    cdef Point c_p

cdef class NomadBBOArrayOfDouble:
    cdef ArrayOfDouble c_aod

cdef class NomadBBODouble:
    cdef Double c_d