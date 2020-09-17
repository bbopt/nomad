/**
 \file   CallbackType.hpp
 \brief  Types for Callback
 \author Viviane Rochon Montplaisir
 \date   November 2019
 \see    CallbackType.cpp
 */


#ifndef __NOMAD400_CALLBACK_TYPE__
#define __NOMAD400_CALLBACK_TYPE__

#include "../nomad_nsbegin.hpp"

enum class CallbackType
{
    ITERATION_END,      ///< Called at the end of an Iteration
    MEGA_ITERATION_END, ///< Called at the end of a MegaIteration
    HOT_RESTART         ///< Called at the beginning of Hot Restart process
};


#include "../nomad_nsend.hpp"
#endif  // __NOMAD400_CALLBACK_TYPE__
