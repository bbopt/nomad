/*---------------------------------------------------------------------------------*/
/*  NOMAD - Nonlinear Optimization by Mesh Adaptive Direct Search -                */
/*                                                                                 */
/*  NOMAD - Version 4 has been created and developed by                            */
/*                 Viviane Rochon Montplaisir  - Polytechnique Montreal            */
/*                 Christophe Tribes           - Polytechnique Montreal            */
/*                                                                                 */
/*  The copyright of NOMAD - version 4 is owned by                                 */
/*                 Charles Audet               - Polytechnique Montreal            */
/*                 Sebastien Le Digabel        - Polytechnique Montreal            */
/*                 Viviane Rochon Montplaisir  - Polytechnique Montreal            */
/*                 Christophe Tribes           - Polytechnique Montreal            */
/*                                                                                 */
/*  NOMAD 4 has been funded by Rio Tinto, Hydro-Québec, Huawei-Canada,             */
/*  NSERC (Natural Sciences and Engineering Research Council of Canada),           */
/*  InnovÉÉ (Innovation en Énergie Électrique) and IVADO (The Institute            */
/*  for Data Valorization)                                                         */
/*                                                                                 */
/*  NOMAD v3 was created and developed by Charles Audet, Sebastien Le Digabel,     */
/*  Christophe Tribes and Viviane Rochon Montplaisir and was funded by AFOSR       */
/*  and Exxon Mobil.                                                               */
/*                                                                                 */
/*  NOMAD v1 and v2 were created and developed by Mark Abramson, Charles Audet,    */
/*  Gilles Couture, and John E. Dennis Jr., and were funded by AFOSR and           */
/*  Exxon Mobil.                                                                   */
/*                                                                                 */
/*  Contact information:                                                           */
/*    Polytechnique Montreal - GERAD                                               */
/*    C.P. 6079, Succ. Centre-ville, Montreal (Quebec) H3C 3A7 Canada              */
/*    e-mail: nomad@gerad.ca                                                       */
/*                                                                                 */
/*  This program is free software: you can redistribute it and/or modify it        */
/*  under the terms of the GNU Lesser General Public License as published by       */
/*  the Free Software Foundation, either version 3 of the License, or (at your     */
/*  option) any later version.                                                     */
/*                                                                                 */
/*  This program is distributed in the hope that it will be useful, but WITHOUT    */
/*  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or          */
/*  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License    */
/*  for more details.                                                              */
/*                                                                                 */
/*  You should have received a copy of the GNU Lesser General Public License       */
/*  along with this program. If not, see <http://www.gnu.org/licenses/>.           */
/*                                                                                 */
/*  You can find information on the NOMAD software at www.gerad.ca/nomad           */
/*---------------------------------------------------------------------------------*/

#include "../../Algos/SgtelibModel/SgtelibModel.hpp"
#include "../../Algos/SgtelibModel/SgtelibModelEvaluator.hpp"
#include "../../Output/OutputQueue.hpp"
#include "../../Type/SgtelibModelFormulationType.hpp"

#include "../../../ext/sgtelib/src/Surrogate.hpp"

// Constructor
// Usually, Evaluators work in full dimension.
// In this case, the models may work better in sub dimension.
// This is why a fixed variable is used.
NOMAD::SgtelibModelEvaluator::SgtelibModelEvaluator(const std::shared_ptr<NOMAD::EvalParameters>& evalParams,
                                   const NOMAD::SgtelibModel* modelAlgo,
                                   const std::string& modelDisplay,
                                   const NOMAD::Double& diversification,
                                   const NOMAD::SgtelibModelFeasibilityType& modelFeasibility,
                                   const double tc,
                                   const NOMAD::Point& fixedVariable)
  : NOMAD::Evaluator(evalParams, NOMAD::EvalType::MODEL),
    _modelAlgo(modelAlgo),
    _modelDisplay(modelDisplay),
    _diversification(diversification),
    _modelFeasibility(modelFeasibility),
    _tc(tc),
    _displayLevel(NOMAD::OutputLevel::LEVEL_INFO),
    _fixedVariable(fixedVariable)
{
    init();
}


// Destructor
NOMAD::SgtelibModelEvaluator::~SgtelibModelEvaluator() = default;


void NOMAD::SgtelibModelEvaluator::init()
{
    _displayLevel = (std::string::npos != _modelDisplay.find('X'))
                        ? NOMAD::OutputLevel::LEVEL_INFO
                        : NOMAD::OutputLevel::LEVEL_DEBUGDEBUG;
    
}


/*------------------------------------------------------------------------*/
/*                evaluate the sgtelib_model model at a given point       */
/*------------------------------------------------------------------------*/
bool NOMAD::SgtelibModelEvaluator::eval_x(NOMAD::EvalPoint &x,
                                          const NOMAD::Double &hMax,
                                          bool &countEval) const
{
    
    // Convert x to subspace, because model is in subspace.
    x = x.makeSubSpacePointFromFixed(_fixedVariable);

    std::string s;
    OUTPUT_INFO_START
    s = "X = " + x.display();
    NOMAD::OutputQueue::Add(s, _displayLevel);
    OUTPUT_INFO_END

    size_t n = x.size();
    size_t nbConstraints = NOMAD::getNbConstraints(_bbOutputTypeList);
    size_t nbModels = NOMAD::SgtelibModel::getNbModels(_modelFeasibility, nbConstraints);
    // Init the matrices for prediction
    SGTELIB::Matrix   M_predict (  "M_predict", 1, static_cast<int>(nbModels));
    SGTELIB::Matrix STD_predict ("STD_predict", 1, static_cast<int>(nbModels));
    SGTELIB::Matrix CDF_predict ("CDF_predict", 1, static_cast<int>(nbModels));
    SGTELIB::Matrix  EI_predict ( "EI_predict", 1, static_cast<int>(nbModels));


    // --------------------- //
    // In/Out Initialisation //
    // --------------------- //
    // Declaration of the statistical measurements
    NOMAD::Double pf = 1; // P[x]
    NOMAD::Double f = 0; // predicted mean of the objective
    NOMAD::Double sigma_f = 0; // predicted variance of the objective
    NOMAD::Double pi = 0; // probability of improvement
    NOMAD::Double ei = 0; // expected improvement
    NOMAD::Double efi = 0; // expected feasible improvement
    NOMAD::Double pfi = 0; // probability of feasible improvement
    NOMAD::Double mu = 0; // uncertainty on the feasibility
    NOMAD::Double penalty = 0; // exclusion area penalty
    NOMAD::Double d = 0; // Distance to the closest point of the cache

    // Creation of matrix for input / output of SGTELIB model
    SGTELIB::Matrix X_predict("X_predict", 1, static_cast<int>(n));
    // Shall we compute statistical criteria
    bool useStatisticalCriteria = false;
    // FORMULATION USED IN THIS EVAL_X
    const NOMAD::SgtelibModelFormulationType formulation = _modelAlgo->getFormulation();

    // Unfortunately, Sgtelib is not thread-safe.
    // For this reason we have to set part of the eval_x code to critical.
#ifdef _OPENMP
    #pragma omp critical(SgtelibEvalX)
#endif // _OPENMP
    {
        // Set the input matrix
        for (size_t i = 0; i < n; i++)
        {
            X_predict.set(0, static_cast<int>(i), x[i].todouble());
        }

        // Reset point outputs
        // By default, set everything to -1
        // Note: Currently NOMAD cannot set a bbo value by index, so we have to
        // work around by constructing a suitable string.
        // Note: Why set some default values on bbo?
        NOMAD::ArrayOfString defbbo(_bbOutputTypeList.size(), "-1");
        x.setBBO(defbbo.display(), _bbOutputTypeList, _evalType);

        // ------------------------- //
        //   Objective Prediction    //
        // ------------------------- //
        switch (formulation)
        {
            case NOMAD::SgtelibModelFormulationType::FS:
                useStatisticalCriteria = (_diversification != 0);
                break;
            case NOMAD::SgtelibModelFormulationType::FSP:
            case NOMAD::SgtelibModelFormulationType::EIS:
            case NOMAD::SgtelibModelFormulationType::EFI:
            case NOMAD::SgtelibModelFormulationType::EFIS:
            case NOMAD::SgtelibModelFormulationType::EFIM:
            case NOMAD::SgtelibModelFormulationType::EFIC:
            case NOMAD::SgtelibModelFormulationType::PFI:
                useStatisticalCriteria = true;
                break;
            case NOMAD::SgtelibModelFormulationType::D:
                useStatisticalCriteria = false;
                break;
            case NOMAD::SgtelibModelFormulationType::EXTERN:
            case NOMAD::SgtelibModelFormulationType::UNDEFINED:
            default:
                throw SGTELIB::Exception ( __FILE__ , __LINE__ , "Forbidden formulation" );
                break;
        }

        OUTPUT_INFO_START
        s = "Formulation: " + NOMAD::SgtelibModelFormulationTypeToString(formulation);
        s += "; compute stat: " + NOMAD::boolToString(useStatisticalCriteria);
        s += "; found feasible: " + NOMAD::boolToString(_modelAlgo->getFoundFeasible());
        NOMAD::OutputQueue::Add(s, _displayLevel);
        OUTPUT_INFO_END

        // Prediction
        if ( formulation == NOMAD::SgtelibModelFormulationType::D )
        {
            OUTPUT_INFO_START
            d = _modelAlgo->getTrainingSet()->get_distance_to_closest(X_predict).get(0,0);
            s = "d = " + d.display();
            NOMAD::OutputQueue::Add(s, _displayLevel);
            OUTPUT_INFO_END
        }
        else if (formulation == NOMAD::SgtelibModelFormulationType::EXTERN)
        {
            throw SGTELIB::Exception(__FILE__, __LINE__,
                                     "SgtelibModelEvaluator::eval_x: Formulation Extern should not been called in this context.");
        }
        else
        {
            OUTPUT_INFO_START
            NOMAD::OutputQueue::Add("Predict... ", _displayLevel);
            OUTPUT_INFO_END

            auto model = _modelAlgo->getModel();
            model->check_ready(__FILE__,__FUNCTION__,__LINE__);

            if (useStatisticalCriteria)
            {
                model->predict(X_predict, &M_predict, &STD_predict, &EI_predict, &CDF_predict);
            }
            else
            {
                model->predict(X_predict, &M_predict);
            }

            OUTPUT_INFO_START
            NOMAD::OutputQueue::Add("ok", _displayLevel);
            OUTPUT_INFO_END
        }

        // Get the prediction from the matrices
        f = M_predict.get(0,0);

        if (useStatisticalCriteria)
        {
            // If no feasible points is found so far, then sigma_f, ei and pi are bypassed.
            if (_modelAlgo->getFoundFeasible())
            {
                sigma_f = STD_predict.get(0,0);
                pi      = CDF_predict.get(0,0);
                ei      = EI_predict.get(0,0);
            }
            else
            {
                sigma_f = 1.0; // This inhibits the exploration term in regard to the objective
                pi      = 1.0; // This implies that pfi = pf
                ei      = 1.0; // This implies that efi = pf
            }
            OUTPUT_INFO_START
            s = "F = " + f.display() + " +/- " + sigma_f.display();
            NOMAD::OutputQueue::Add(s, _displayLevel);
            OUTPUT_INFO_END
        }
        else
        {
            OUTPUT_INFO_START
            s = "F = " + f.display();
            NOMAD::OutputQueue::Add(s, _displayLevel);
            OUTPUT_INFO_END
        }

        // ====================================== //
        // Constraints display                    //
        // ====================================== //
        OUTPUT_INFO_START
        s = "";
        switch (_modelFeasibility)
        {
            case NOMAD::SgtelibModelFeasibilityType::C:
                if (useStatisticalCriteria)
                {
                    for (size_t i = 1; i < nbModels; i++)
                    {
                        s += "C" + std::to_string(i) + " = " + std::to_string(M_predict.get(0,static_cast<int>(i)));
                        s += " +/- " + std::to_string(STD_predict.get(0,static_cast<int>(i)));
                        s += " (CDF : " + std::to_string(CDF_predict.get(0,static_cast<int>(i))) +  ")";
                    }
                }
                else
                {
                    s += "C = [ ";
                    for (size_t i = 1; i < nbModels; i++)
                    {
                        s += std::to_string(M_predict.get(0,static_cast<int>(i))) + " ";
                        s += " ]";
                    }
                }
                break;

            case NOMAD::SgtelibModelFeasibilityType::H:
                s += "Feasibility_Method : H (Aggregate prediction)";
                s += "H = " + std::to_string(M_predict.get(0,1));
                s += " +/- " + std::to_string(STD_predict.get(0,1));
                s += " (CDF : " + std::to_string(CDF_predict.get(0,1)) + ")";
                break;
            case NOMAD::SgtelibModelFeasibilityType::B:
                s += "Feasibility_Method : B (binary prediction)";
                s += "B = " + std::to_string(M_predict.get(0,1));
                s += " (CDF : " + std::to_string(CDF_predict.get(0,1)) + ")";
                break;
            case NOMAD::SgtelibModelFeasibilityType::M:
                s += "Feasibility_Method : M (Biggest constraint prediction)";
                s += "M = " + std::to_string(M_predict.get(0,1));
                s += " +/- " + std::to_string(STD_predict.get(0,1));
                s += " (CDF : " + std::to_string(CDF_predict.get(0,1)) + ")";
                break;
            case NOMAD::SgtelibModelFeasibilityType::UNDEFINED:
            default:
                s = "SGTELIB_MODEL_FEASIBILITY_UNDEFINED";
                NOMAD::OutputQueue::Add(s, _displayLevel);
                break;
        }
        OUTPUT_INFO_END

        // ====================================== //
        // Computation of statistical criteria    //
        // ====================================== //
        if ( useStatisticalCriteria )
        {
            pf = 1; // General probability of feasibility
            NOMAD::Double pfj; // Probability of feasibility for constraint cj
            NOMAD::Double L2 = 0;

            if (nbConstraints > 0)
            {
                // Use the CDF of each output in C
                // If there is only one output in C (models B, H and M) then pf = CDF
                for (size_t i = 1; i < nbModels; i++)
                {
                    pfj = CDF_predict.get(0,static_cast<int>(i));
                    L2 += max( 0 , M_predict.get(0,static_cast<int>(i))).pow2();
                    pf *= pfj;
                }
            }   // end (if constraints are present)
            if (!_modelAlgo->getFoundFeasible() && (pf == 0))
            {
                pf = 1.0/(1.0+L2);
                OUTPUT_INFO_START
                s = "pf = 0 and L2 = " + L2.display() + " => pF = " + pf.display();
                NOMAD::OutputQueue::Add(s, _displayLevel);
                OUTPUT_INFO_END
            }
            pfi = pi*pf;
            efi = ei*pf;
            mu = 4*pf*(1-pf);
        }

    } // pragma omp critical

    // ====================================== //
    // Application of the formulation         //
    // ====================================== //
    NOMAD::Double obj;
    NOMAD::ArrayOfDouble newbbo(_bbOutputTypeList.size(), -1);
    int k = 0;
    switch (formulation)
    {
        case NOMAD::SgtelibModelFormulationType::FS:
            // Define obj
            obj = f - _diversification*sigma_f;
            // Set constraints
            for (size_t i = 0; i < _bbOutputTypeList.size(); i++)
            {
                if (_bbOutputTypeList[i] != NOMAD::BBOutputType::OBJ)
                {
                    newbbo[i] = M_predict.get(0,k+1) - _diversification*STD_predict.get(0,k+1);
                    k++;
                }
            }
            break;

        case NOMAD::SgtelibModelFormulationType::FSP:
            // Define obj
            obj = f - _diversification*sigma_f;
            // Set constraints
            for (size_t i = 0; i < _bbOutputTypeList.size(); i++)
            {
                if (_bbOutputTypeList[i] != NOMAD::BBOutputType::OBJ)
                {
                    newbbo[i] = 0.5 - pf;
                }
            }
            break;

        case NOMAD::SgtelibModelFormulationType::EIS:
            // Define obj
            obj = - ei - _diversification*sigma_f;
            // Set constraints
            for (size_t i = 0; i < _bbOutputTypeList.size(); i++)
            {
                if ( _bbOutputTypeList[i] != NOMAD::BBOutputType::OBJ )
                {
                    newbbo[i] = M_predict.get(0,k+1) - _diversification*STD_predict.get(0,k+1);
                    k++;
                }
            }
            break;

        case NOMAD::SgtelibModelFormulationType::EFI:
            obj = -efi;
            break;

        case NOMAD::SgtelibModelFormulationType::EFIS:
            obj = -efi - _diversification*sigma_f;
            break;

        case NOMAD::SgtelibModelFormulationType::EFIM:
            obj = -efi - _diversification*sigma_f*mu;
            break;

        case NOMAD::SgtelibModelFormulationType::EFIC:
            obj = -efi - _diversification*( ei*mu + pf*sigma_f);
            break;

        case NOMAD::SgtelibModelFormulationType::PFI:
            obj = -pfi;
            break;

        case NOMAD::SgtelibModelFormulationType::D:
            obj = -d;
            break;

        case NOMAD::SgtelibModelFormulationType::EXTERN:
        case NOMAD::SgtelibModelFormulationType::UNDEFINED:
        default:
            OUTPUT_INFO_START
            s = "SgtelibModel formulation: " + NOMAD::SgtelibModelFormulationTypeToString(formulation);
            NOMAD::OutputQueue::Add(s, _displayLevel);
            OUTPUT_INFO_END
            break;
    }

    // ------------------------- //
    //   exclusion area          //
    // ------------------------- //
    if (_tc > 0.0)
    {
        penalty = _modelAlgo->getModel()->get_exclusion_area_penalty(X_predict, _tc).get(0,0);
        obj += penalty;
    }

    // ------------------------- //
    //   Set obj and BBO         //
    // ------------------------- //
    for (size_t i = 0; i < _bbOutputTypeList.size(); i++)
    {
        if (_bbOutputTypeList[i] == NOMAD::BBOutputType::OBJ)
        {
            newbbo[i] = obj;
        }
    }
    x.setBBO(newbbo.tostring(), _bbOutputTypeList, NOMAD::EvalType::MODEL);

    // ================== //
    //       DISPLAY      //
    // ================== //
    OUTPUT_INFO_START
    if (useStatisticalCriteria)
    {
        s = "f_min                    f_min = " + std::to_string(_modelAlgo->getTrainingSet()->get_f_min());
        NOMAD::OutputQueue::Add(s, _displayLevel);
        s = "Probability of Feasibility PF  = " + pf.display();
        NOMAD::OutputQueue::Add(s, _displayLevel);
        s = "Feasibility Uncertainty    mu  = " + mu.display();
        NOMAD::OutputQueue::Add(s, _displayLevel);
        s = "Probability Improvement    PI  = " + pi.display();
        NOMAD::OutputQueue::Add(s, _displayLevel);
        s = "Expected Improvement      EI  = " + ei.display();
        NOMAD::OutputQueue::Add(s, _displayLevel);
        s = "Proba. of Feasible Imp.    PFI = " + pfi.display();
        NOMAD::OutputQueue::Add(s, _displayLevel);
        s = "Expected Feasible Imp.     EFI = " + efi.display();
        NOMAD::OutputQueue::Add(s, _displayLevel);
    }
    s = "Exclusion area penalty = " + penalty.display();
    NOMAD::OutputQueue::Add(s, _displayLevel);
    s = "Model Output = (" + x.getBBO(NOMAD::EvalType::MODEL) + ")";
    NOMAD::OutputQueue::Add(s, _displayLevel);
    OUTPUT_INFO_END

    if (!pf.isDefined() || !pi.isDefined())
    {
        throw SGTELIB::Exception ( __FILE__ , __LINE__ ,
                                  "SgtelibModelEvaluator::eval_x: NaN values in pi or pf." );
    }

    // ================== //
    // Exit Status        //
    // ================== //
    countEval = true;
    x.setEvalStatus(NOMAD::EvalStatusType::EVAL_OK, NOMAD::EvalType::MODEL);

    // Convert back x to full space
    x = x.makeFullSpacePointFromFixed(_fixedVariable);

    // Always eval_ok = true
    return true;
}
