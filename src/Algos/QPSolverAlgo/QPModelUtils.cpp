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
#include "QPModelUtils.hpp"

#include <map>

#include "../../Util/Exception.hpp"

double NOMAD::QPModelUtils::getModelObj(const SGTELIB::Matrix& QPModel,
                                        const SGTELIB::Matrix& x)
{
    return NOMAD::QPModelUtils::getModelValue(QPModel, 0, x);
}

void NOMAD::QPModelUtils::getModelObjGrad(SGTELIB::Matrix& g,
                                          const SGTELIB::Matrix& QPModel,
                                          const SGTELIB::Matrix& x)
{
    NOMAD::QPModelUtils::getModelGrad(g, QPModel, 0, x);
}


double NOMAD::QPModelUtils::getModelCons(const SGTELIB::Matrix& QPModel,
                                         const int ind,
                                         const SGTELIB::Matrix& x)
{
    return NOMAD::QPModelUtils::getModelValue(QPModel, ind + 1, x);
}

void NOMAD::QPModelUtils::getModelCons(SGTELIB::Matrix& cons,
                                       const SGTELIB::Matrix& QPModel,
                                       const SGTELIB::Matrix& x)
{
    const int nbCons = QPModel.get_nb_rows() - 1;
    if (cons.get_nb_cols() != 1 || cons.get_nb_rows() != nbCons)
    {
        std::string err = "QPModelUtils::getModelCons: ";
        err += "the number of constraints of the model " + std::to_string(nbCons);
        err += "is not compatible with the dimensions of the cons vector: ( " + std::to_string(cons.get_nb_rows());
        err += std::to_string(cons.get_nb_cols()) + " )";
        throw NOMAD::Exception(__FILE__, __LINE__, err);
    }

    for (int i = 0; i < nbCons; ++i)
    {
        cons.set(i, 0, NOMAD::QPModelUtils::getModelValue(QPModel, i+1, x));
    }
}

void NOMAD::QPModelUtils::getModelJacobianCons(SGTELIB::Matrix& jacobian,
                                               const SGTELIB::Matrix& QPModel,
                                               const SGTELIB::Matrix& x)
{
    const int nbCons = QPModel.get_nb_rows() - 1;
    const int n = std::max(x.get_nb_rows(), x.get_nb_cols());
    if (jacobian.get_nb_rows() != nbCons || jacobian.get_nb_cols() != n)
    {
        std::string err = "QPModelUtils::getModelJacobianCons: ";
        err += "The dimensions of the jacobian constraint matrix are ( " + std::to_string(jacobian.get_nb_rows());
        err += std::to_string(jacobian.get_nb_cols()) + " ) ";
        err += "when the model has " + std::to_string(nbCons) + " constraints and " + std::to_string(n) + " variables";
        throw NOMAD::Exception(__FILE__, __LINE__, err);
    }

    SGTELIB::Matrix g("g", n, 1);
    for (int ind = 0; ind < nbCons; ++ind)
    {
        NOMAD::QPModelUtils::getModelGrad(g, QPModel, ind + 1, x);
        for (int i = 0; i < n; ++i)
        {
            const double gi = g.get(i, 0);
            jacobian.set(ind, i, gi);
        }
    }
}

double NOMAD::QPModelUtils::getModelValue(const SGTELIB::Matrix& QPModel,
                                          const int ind,
                                          const SGTELIB::Matrix &x)
{
    const int nbModels = QPModel.get_nb_rows();
    if (ind >= nbModels)
    {
        std::string err = "QPModelUtils::getModelValue: ";
        err += "the index of the requested model is superior to the number of models";
        throw NOMAD::Exception(__FILE__, __LINE__, err);
    }

    const int nbParams = QPModel.get_nb_cols();
    const int n = std::max(x.get_nb_rows(), x.get_nb_cols());
    if (nbParams != (n+1) + n * (n+1) / 2)
    {
        std::string err = "QPModelUtils::getModelValue: ";
        err += "the number of parameters of the model (nbParams = " + std::to_string(nbParams);
        err += ") is not compatible with the number of variables of the problem (n = " + std::to_string(n) + ")";
        throw NOMAD::Exception(__FILE__, __LINE__, err);
    }

    const bool isColOrder = n == x.get_nb_rows();

    double modelVal = QPModel.get(ind, 0);
    int k = 1;
    for (int j = 0; j < n; ++j)
    {
        const double xj = isColOrder ? x.get(j, 0) : x.get(0, j);
        modelVal += xj * (QPModel.get(ind, k) + 0.5 * QPModel.get(ind, k + n) * xj);
        ++k;
    }

    k += n;
    for (int i = 0; i < n; ++i)
    {
        const double xi = isColOrder ? x.get(i, 0) : x.get(0, i);
        for (int j = 0; j < i; ++j)
        {
            const double xj = isColOrder ? x.get(j, 0) : x.get(0, j);
            modelVal += QPModel.get(ind, k) * xi * xj;
            ++k;
        }
    }

    return modelVal;
}

void NOMAD::QPModelUtils::getModelGrad(SGTELIB::Matrix& g,
                                       const SGTELIB::Matrix& QPModel,
                                       const int ind,
                                       const SGTELIB::Matrix& x)
{
    const int nbModels = QPModel.get_nb_rows();
    if (ind >= nbModels)
    {
        std::string err = "QPModelUtils::getModelGrad: ";
        err += "the index of the requested model is superior to the number of models";
        throw NOMAD::Exception(__FILE__, __LINE__, err);
    }

    const int n = std::max(x.get_nb_rows(), x.get_nb_cols());
    if (g.get_nb_rows() != n || g.get_nb_cols() != 1)
    {
        throw NOMAD::Exception(__FILE__, __LINE__, "QPModelUtils::getModelGrad: " + g.get_name() + " has wrong dimensions");
    }

    const int nbParams = QPModel.get_nb_cols();
    if (nbParams != (n+1) + n * (n+1) / 2)
    {
        std::string err = "QPModelUtils::getModelGrad: ";
        err += "the number of parameters of the model (nbParams = " + std::to_string(nbParams);
        err += ") is not compatible with the number of variables of the problem (n = " + std::to_string(n) + ")";
        throw NOMAD::Exception(__FILE__, __LINE__, err);
    }

    const bool isColOrder = n == x.get_nb_rows();

    int k = 1;
    for (int j = 0; j < n; ++j)
    {
        const double xj = isColOrder ? x.get(j, 0) : x.get(0, j);
        g.set(j, 0, QPModel.get(ind, k) + QPModel.get(ind, k + n) * xj);
        ++k;
    }

    k += n;
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < i; ++j)
        {
            // Differentiate according to xj
            double gj = g.get(j, 0);
            const double xi = isColOrder ? x.get(i, 0) : x.get(0, i);
            gj += QPModel.get(ind, k) * xi;
            g.set(j, 0, gj);

            // Differentiate according to xi
            double gi = g.get(i, 0);
            const double xj = isColOrder ? x.get(j, 0) : x.get(0, j);
            gi += QPModel.get(ind, k) * xj;
            g.set(i, 0, gi);

            ++k;
        }
    }
}

double NOMAD::QPModelUtils::getModelLagrangian(const SGTELIB::Matrix& QPModel,
                                               const SGTELIB::Matrix& x,
                                               const SGTELIB::Matrix& lambda,
                                               const double sigma)
{
    if (lambda.get_nb_cols() != 1)
    {
        std::string err = "QPModelUtils::getModelLagrangianValue: ";
        err += "The lagrange multipliers vector must be a column vector";
        throw NOMAD::Exception(__FILE__, __LINE__, err);
    }

    const int nbCons = lambda.get_nb_rows();
    if (nbCons != QPModel.get_nb_rows() - 1)
    {
        std::string err = "QPModelUtils::getModelLagrangianValue: ";
        err += "the dimension of the lagrange multipliers vector " + std::to_string(nbCons);
        err += " is incompatible with the number of constraints provided by the QP model ";
        err += std::to_string(QPModel.get_nb_rows() - 1);
        throw NOMAD::Exception(__FILE__, __LINE__, err);
    }

    double lagVal = sigma * getModelObj(QPModel, x);
    for (int j = 0; j < nbCons; ++j)
    {
        lagVal -= lambda.get(j, 0) * getModelValue(QPModel, j + 1, x);
    }

    return lagVal;
}

void NOMAD::QPModelUtils::getModelLagrangianGrad(SGTELIB::Matrix& gradL,
                                                 const SGTELIB::Matrix& QPModel,
                                                 const SGTELIB::Matrix& x,
                                                 const SGTELIB::Matrix& lambda,
                                                 const double sigma)
{
    if (lambda.get_nb_cols() != 1)
    {
        std::string err = "QPModelUtils::getModelLagrangianGrad: ";
        err += "The lagrange multipliers vector must be a column vector";
        throw NOMAD::Exception(__FILE__, __LINE__, err);
    }

    const int nbCons = lambda.get_nb_rows();
    if (nbCons != QPModel.get_nb_rows() - 1)
    {
        std::string err = "QPModelUtils::getModelLagrangianGrad: ";
        err += "the dimension of the lagrange multipliers vector " + std::to_string(nbCons);
        err += " is incompatible with the number of constraints provided by the QP model ";
        err += std::to_string(QPModel.get_nb_rows() - 1);
        throw NOMAD::Exception(__FILE__, __LINE__, err);
    }

    getModelObjGrad(gradL, QPModel, x);
    gradL.multiply(sigma);

    const int n = std::max(x.get_nb_rows(), x.get_nb_cols());
    SGTELIB::Matrix jacobianCons("jacobianCons", nbCons, n);
    getModelJacobianCons(jacobianCons, QPModel, x);

    SGTELIB::Matrix outGrad("outGrad", n, 1);
    SGTELIB::Matrix::inplace_product(outGrad, jacobianCons.transpose(), lambda);
    outGrad.multiply(-1.0);

    gradL.add(outGrad);
}

void NOMAD::QPModelUtils::getModelLagrangianHessian(SGTELIB::Matrix& H,
                                                    const SGTELIB::Matrix& QPModel,
                                                    const SGTELIB::Matrix& x,
                                                    const SGTELIB::Matrix& lambda,
                                                    const double sigma)
{
    if (lambda.get_nb_cols() != 1)
    {
        std::string err = "QPModelUtils::getModelLagrangianHessian: ";
        err += "The lagrange multipliers vector must be a column vector";
        throw NOMAD::Exception(__FILE__, __LINE__, err);
    }

    const int nbCons = lambda.get_nb_rows();
    if (nbCons != QPModel.get_nb_rows() - 1)
    {
        std::string err = "QPModelUtils::getModelLagrangianHessian: ";
        err += "the dimension of the lagrange multipliers vector " + std::to_string(nbCons);
        err += " is incompatible with the number of constraints provided by the QP model ";
        err += std::to_string(QPModel.get_nb_rows() - 1);
        throw NOMAD::Exception(__FILE__, __LINE__, err);
    }

    getModelHessian(H, QPModel, 0, x);
    H.multiply(sigma);

    SGTELIB::Matrix Htmp = H;
    for (int ind = 0; ind < nbCons; ++ind)
    {
        getModelHessian(Htmp, QPModel, ind + 1, x);
        Htmp.multiply(-lambda.get(ind, 0));
        H.add(Htmp);
    }
}

void NOMAD::QPModelUtils::getModelHessian(SGTELIB::Matrix& H,
                                          const SGTELIB::Matrix& QPModel,
                                          const int ind,
                                          const SGTELIB::Matrix& x)
{
    const int nbModels = QPModel.get_nb_rows();
    if (ind >= nbModels)
    {
        std::string err = "QPModelUtils::getModelHessian: ";
        err += "the index of the requested model is superior to the number of models";
        throw NOMAD::Exception(__FILE__, __LINE__, err);
    }

    const int n = std::max(x.get_nb_rows(), x.get_nb_cols());
    if (H.get_nb_rows() != n || H.get_nb_cols() != n)
    {
        std::string err = "QPModelUtils::getModelHessian: the dimension of " + H.get_name() + " (dim(H) = (";
        err += std::to_string(H.get_nb_rows()) + ", " + std::to_string(H.get_nb_cols()) + "))";
        err += " is not compatible with the number of variables of x (n = " + std::to_string(n) + ")";
        throw NOMAD::Exception(__FILE__, __LINE__, err);
    }

    const int nbParams = QPModel.get_nb_cols();
    if (nbParams != (n+1) + n * (n+1) / 2)
    {
        std::string err = "QPModelUtils::getModelHessian: ";
        err += "the number of parameters of the model " + std::to_string(nbParams);
        err += "is not compatible with the number of variables of the problem " + std::to_string(n);
        throw NOMAD::Exception(__FILE__, __LINE__, err);
    }

    int k = n + 1;
    for (int i =  0; i < n; ++i)
    {
        H.set(i, i, QPModel.get(ind, k));
        ++k;
    }

    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < i; ++j)
        {
            const double alphak = QPModel.get(ind, k);
            H.set(i, j, alphak);
            H.set(j, i, alphak);
            ++k;
        }
    }
}

SGTELIB::Matrix NOMAD::QPModelUtils::getReducedQPModel(const SGTELIB::Matrix& QPModel,
                                                       const SGTELIB::Matrix& x,
                                                       const std::vector<bool>& fixedVars)
{
    const int n = x.get_nb_rows();
    if (fixedVars.size() != n)
    {
        std::string err = "QPModelUtils::getReducedQPModel: ";
        err += "the number of variables of the problem " + std::to_string(n);
        err += "is not compatible with the dimension of the boolean vector fixedVars " + std::to_string(fixedVars.size());
        throw NOMAD::Exception(__FILE__, __LINE__, err);
    }

    const int nbParams = QPModel.get_nb_cols();
    if (nbParams != (n+1) + n * (n+1) / 2)
    {
        std::string err = "QPModelUtils::getModelReducedQPModel: ";
        err += "the number of parameters of the model " + std::to_string(nbParams);
        err += "is not compatible with the number of variables of the problem " + std::to_string(n);
        throw NOMAD::Exception(__FILE__, __LINE__, err);
    }

    const int nbFixedVars = (int) std::count(fixedVars.cbegin(), fixedVars.cend(), true);
    const int nbFreeVars = n - nbFixedVars;
    if (nbFreeVars == n)
    {
        return QPModel;
    }

    // Map coefficient of non-fixed variables in QPModel to QPReducedModel
    std::map<int, int> fullModelToReducedModelIndexes;
    int curInd = 0;
    for (int i = 0; i < n; ++i)
    {
        if (fixedVars[i])
            continue;

        fullModelToReducedModelIndexes.insert({i, curInd});
        curInd++;
    }

    const int nbParamsRed = (nbFreeVars + 1) + nbFreeVars * (nbFreeVars + 1) / 2;
    const int nbOutputs = QPModel.get_nb_rows();
    SGTELIB::Matrix QPModelReduced("QPModelRed", nbOutputs, nbParamsRed);
    for (int row = 0; row < nbOutputs; ++row)
    {
        double alpha0Val = QPModel.get(row, 0);
        curInd = 1;
        for (int j = 1; j < n + 1; j++)
        {
            const double xj = x.get(j-1, 0);
            if (fixedVars[j-1])
            {
                alpha0Val += xj * (QPModel.get(row, j) + 0.5 * QPModel.get(row, j + n) * xj);
            }
            else
            {
                QPModelReduced.set(row, curInd, QPModel.get(row, j));
                QPModelReduced.set(row, curInd + nbFreeVars, QPModel.get(row, j + n));
                curInd++;
            }
        }

        // At this index, the coefficients of the Hessian matrix of the QPModel start.
        int k = 2 * n + 1;
        // At this index, the coefficients of the Hessian matrix of the reduced model start.
        curInd = 2 * nbFreeVars + 1;
        for (int i = 0; i < n; ++i)
        {
            const double xi = x.get(i, 0);
            for (int j = 0; j < i; ++j)
            {
                const double xj = x.get(j, 0);
                const double alphaK = QPModel.get(row, k);
                if (fixedVars[i] & fixedVars[j])
                {
                    alpha0Val += xi * xj * alphaK;
                }
                else if (fixedVars[i])
                {
                    const int QPredModInd = fullModelToReducedModelIndexes[j];
                    double alphaRed = QPModelReduced.get(row, QPredModInd + 1);
                    alphaRed += alphaK * xi;
                    QPModelReduced.set(row, QPredModInd + 1, alphaRed);
                }
                else if (fixedVars[j])
                {
                    const int QPredModInd = fullModelToReducedModelIndexes[i];
                    double alphaRed = QPModelReduced.get(row, QPredModInd + 1);
                    alphaRed += alphaK * xj;
                    QPModelReduced.set(row, QPredModInd + 1, alphaRed);
                }
                else
                {
                    QPModelReduced.set(row, curInd, alphaK);
                    curInd++;
                }
                k++;
            }
        }

        QPModelReduced.set(row, 0, alpha0Val);
    }
    return QPModelReduced;
}
