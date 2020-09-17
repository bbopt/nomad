#ifndef __NOMAD400_SGTELIBSEARCHMETHOD__
#define __NOMAD400_SGTELIBSEARCHMETHOD__

#include "../../Algos/Mads/SearchMethodAlgo.hpp"
#ifdef USE_SGTELIB
#include "../../Algos/SgtelibModel/SgtelibModel.hpp"
#endif

#include "../../nomad_nsbegin.hpp"

// \class SgtelibSearchMethod: Search method using library Sgtelib
class SgtelibSearchMethod final: public SearchMethodAlgo
{
private:
    OutputLevel _displayLevel;
#ifdef USE_SGTELIB
    std::shared_ptr<SgtelibModel> _modelAlgo;
#endif

    // TODO finish implementing and use Projection class
    /// Get best projection
    /**
     \param  incumbent      The incumbent             -- \b IN.
     \param  deltaMeshSize  Mesh size parameter       -- \b IN.
     \param  x              The oracle point          -- \b IN/OUT.
     */
    void getBestProjection(const Point& incumbent,
                           const ArrayOfDouble& deltaMeshSize,
                           std::shared_ptr<Point> x);


/*----------------------------------------------------------------------*/


public:
    /// Constructor
    /**
     /param parentStep      The parent of this search step -- \b IN.
     */
    explicit SgtelibSearchMethod(const Step* parentStep)
      : SearchMethodAlgo(parentStep),
        _displayLevel(OutputLevel::LEVEL_NORMAL)
#ifdef USE_SGTELIB
        ,_modelAlgo(nullptr)
#endif
    {
        init();
    }

private:
    void init();

    bool runImp() override;

    ///Generate new points (no evaluation)
    /**
     \copydoc SearchMethod::generateTrialPointsImp \n
     This function is used only when a MADS search based on Quad Model with the option to generate all points before evaluation. It performs a single sub optimization (on the sgte) around all the points in the Barrier.
     */
    void generateTrialPointsImp() override;

};

#include "../../nomad_nsend.hpp"

#endif // __NOMAD400_SGTELIBSEARCHMETHOD__

