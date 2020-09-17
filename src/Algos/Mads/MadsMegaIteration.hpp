#ifndef __NOMAD400_MADSMEGAITERATION__
#define __NOMAD400_MADSMEGAITERATION__


#include "../../Algos/MegaIteration.hpp"

#include "../../nomad_nsbegin.hpp"

/// Class for the mega iterations of MADS.
/**
Manager for Mads iterations.
 Steps:
 - Generate a lot of points over multiple meshes, using different Search and Poll strategies.
 - Evaluate points
 - Post-processing

 \note As an hypothesis, the time load is taken by the evaluation,
  which is parallelized over all evaluations simultaneously.
  The iteration generation, including trial points generation,
  has little time load, so they do not need to be parallelized.
  It is also preferable to keep parallelization to the only place where
  it matters the most to avoid errors.
  There is no parallelization at the algorithmic level.
  Algorithms are run in main thread only; Secundary threads are available for evaluations.
*/
class MadsMegaIteration: public MegaIteration
{
protected:

    /**
     Main mesh that holds the mesh size and frame size that we would use in the standard MADS algorithm or other Mesh-based algorithm.
     */
    std::shared_ptr<MeshBase> _mainMesh;

    void init();

public:
    /// Constructor
    /**
     \param parentStep      The parent step of this step -- \b IN.
     \param k               The main iteration counter -- \b IN.
     \param barrier         The barrier for constraints handling -- \b IN.
     \param mesh            Mesh on which other Iteration meshes are based -- \b IN.
     \param success         Success type of the previous MegaIteration. -- \b IN.
     */
    explicit MadsMegaIteration(const Step* parentStep,
                              size_t k,
                              std::shared_ptr<Barrier> barrier,
                              std::shared_ptr<MeshBase> mesh,
                              SuccessType success)
      : MegaIteration(parentStep, k, barrier, success),
        _mainMesh(mesh)
    {
        init();
    }
    // No Destructor needed - keep defaults.

    /// Implementation of the start tasks for MADS mega iteration.
    /**
     Creates a MadsIteration for each frame center and each desired mesh size.
     Use all xFeas and xInf available.
     For now, not using other frame centers.
     */
    virtual void startImp() override ;

    /// Implementation of the run tasks for MADS mega iteration.
    /**
     Manages the generation of points: either all poll and search points are generated all together before starting evaluation using the MegaSearchPoll or they are generated using a MadsIteration with search and poll separately. A run parameter controls the behavior.
     */
    virtual bool runImp() override;


    const std::shared_ptr<MeshBase> getMesh() const          { return _mainMesh; }
    void setMesh(const std::shared_ptr<MeshBase> &mesh)      { _mainMesh = mesh; }

    void read(  std::istream& is ) override;
    void display(  std::ostream& os ) const override ;

private:

    /// Helper for start()
    // Currently not used.
    //bool addIterationsForLargerMeshes(const EvalPoint& x0, size_t &k);

};

/**
 Display useful values so that a new MegaIteration could be constructed using these values.
 */
std::ostream& operator<<(std::ostream& os, const MadsMegaIteration& megaIteration);

/// Get an MegaIteration values from a stream
std::istream& operator>>(std::istream& is, MadsMegaIteration& megaIteration);

#include "../../nomad_nsend.hpp"

#endif // __NOMAD400_MADSMEGAITERATION__
