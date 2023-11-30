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
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>

#include "utils.hpp"

void readXFile(const std::string& xFileName, std::vector<std::vector<double>> &xs)
{
    xs.clear();
    std::ifstream in(xFileName);

    // Get the input points
    while (!in.eof())
    {
        std::vector<double> x(DIMENSION);
        for (size_t i = 0; i < DIMENSION; i++)
        {
            in >> x[i];
        }
        xs.push_back(x);
    }
    in.close();
    xs.pop_back();
}


void readXFile(const std::string& xFileName, NOMAD::ArrayOfPoint &aop)
{
    std::vector<std::vector<double>> xs;
    readXFile(xFileName, xs);

    aop.clear();
    for (auto x : xs)
    {
        aop.push_back(NOMAD::Point(x));
    }
}


void initParams(std::shared_ptr<NOMAD::AllParameters>& params,
                const std::string& cacheFileName,
                const std::vector<std::string>& additionalParams)
{
    // Problem parameters
    params->setAttributeValue("DIMENSION", DIMENSION);             // number of variables
    params->setAttributeValue("LOWER_BOUND", NOMAD::ArrayOfDouble(DIMENSION, -10.0));
    params->setAttributeValue("UPPER_BOUND", NOMAD::ArrayOfDouble(DIMENSION, 10.0));

    NOMAD::BBOutputTypeList bbot;   // Definition of output types
    bbot.push_back(NOMAD::BBOutputType::OBJ);
    params->setAttributeValue("BB_OUTPUT_TYPE", bbot);


    params->setAttributeValue("CACHE_FILE", cacheFileName);

    for (auto paramAsString : additionalParams)
    {
        params->readParamLine(paramAsString);
    }

    // Display parameters
    //params->setAttributeValue("DISPLAY_DEGREE", 4);

    // parameters validation
    params->getPbParams()->doNotShowWarnings();
    params->checkAndComply();
}


std::vector<std::string> readAdditionalParams(const std::string& paramFileName)
{
    std::vector<std::string> addParams;
    std::ifstream in(paramFileName);

    // Get the input points
    while (!in.eof())
    {
        std::string paramLine;
        getline(in, paramLine);
        addParams.push_back(paramLine);
    }
    in.close();
    addParams.pop_back();

    return addParams;
}


std::vector<NOMAD::ArrayOfDouble> readBBOutputFile(const std::string& objFileName, const size_t bbOutputSize)
{
    std::vector<NOMAD::ArrayOfDouble> fxs;
    std::ifstream in(objFileName);

    // Get the outputs from file
    NOMAD::ArrayOfDouble bboutput(bbOutputSize);
    while (!in.eof())
    {
        in >> bboutput;
        fxs.push_back(bboutput);
    }
    in.close();
    fxs.pop_back();

    return fxs;
}

