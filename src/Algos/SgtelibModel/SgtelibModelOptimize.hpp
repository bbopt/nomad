#ifndef __NOMAD400_SGTELIB_MODEL_OPTIMIZE__
#define __NOMAD400_SGTELIB_MODEL_OPTIMIZE__

#include "../../Algos/Mads/Mads.hpp"
#include "../../Algos/Step.hpp"
#include "../../Algos/SgtelibModel/SgtelibModel.hpp"
#include "../../Output/OutputInfo.hpp"  // for OutputLevel

#include "../../nomad_nsbegin.hpp"

class SgtelibModelOptimize : public Step
{
private:
    const SgtelibModel*                     _modelAlgo;
    OutputLevel                             _displayLevel;
    //const Point &                        _incumbent;
    //std::vector<std::shared_ptr<EvalPoint>> _x0s;
    std::shared_ptr<Mads>                   _mads;
    EvalPointSet                            _oraclePoints;

    // Reference to the original problem's RunParameters and PbParameters.
    const std::shared_ptr<RunParameters>    _refRunParams;
    const std::shared_ptr<PbParameters>     _refPbParams;

    // RunParameters and PbParameters converted for model optimization
    std::shared_ptr<RunParameters>          _optRunParams;
    std::shared_ptr<PbParameters>           _optPbParams;


public:
    /// Constructor
    // Parent must explicitely be a (pointer to a) SgtelibModel.
    // Run parameters will be recomputed for model optimization.
    explicit SgtelibModelOptimize(const SgtelibModel* modelAlgo,
                                  const std::shared_ptr<RunParameters> refRunParams,
                                  const std::shared_ptr<PbParameters>  refPbParams)
                                  //const Point& incumbent,
                                  //std::vector<std::shared_ptr<EvalPoint>> x0s,
      : Step(modelAlgo),
        _modelAlgo(modelAlgo),
        _displayLevel(OutputLevel::LEVEL_INFO),
        //_incumbent(incumbent),
        //_x0s(x0s),
        _mads(nullptr),
        _oraclePoints(),
        _refRunParams(refRunParams),
        _refPbParams(refPbParams),
        _optRunParams(nullptr),
        _optPbParams(nullptr)
    {
        init();
    }

    /*-----------*/
    /* Get / Set */
    /*-----------*/
    void setupPbParameters(const ArrayOfDouble& lowerBound,
                           const ArrayOfDouble& upperBound,
                           const ArrayOfDouble& initialMeshSize = ArrayOfDouble(),
                           const ArrayOfDouble& initialFrameSize = ArrayOfDouble());

    const EvalPointSet& getOraclePoints() const { return _oraclePoints; }


private:
    void init();

    virtual void startImp() override;
    virtual bool runImp() override;
    virtual void endImp() override;

    // Helpers
    void setupRunParameters();
    void updateOraclePoints();

};

#include "../../nomad_nsend.hpp"

#endif // __NOMAD400_SGTELIB_MODEL_OPTIMIZE__
