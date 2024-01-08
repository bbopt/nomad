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
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <mpi.h>

#ifdef DEBUG
#define D(x) x
#else
#define D(x)
#endif


using namespace std;

const size_t n = 5;
const size_t m = 3;


const int  TAG_SIGNAL        = 0;
const int  TAG_INPUT         = 1;
const int  TAG_OUTPUT        = 2;
const int  TAG_INDEX         = 3;

char SIGNAL_READY_TO_EVAL  = 'R';
char SIGNAL_DONE_EVAL      = 'D';

void evaluatorStart()
{
    int rank = 0;
    
    MPI_Comm_rank ( MPI_COMM_WORLD, &rank );
    
    if (rank == 0)
    {
        MPI_Finalize();
        return;
    }
    
    MPI_Status status;
    int index=-1;
    while(true)
    {
        // Send signal to evaluator server to inform evaluator (#) is ready to perform
        MPI_Send ( &SIGNAL_READY_TO_EVAL ,
                               1 ,
                               MPI_CHAR,
                               0 ,
                               TAG_SIGNAL ,
			                   MPI_COMM_WORLD);
        
        D ( cerr<< "Evaluator with rank #" << rank << " wait to receive an index" << endl;);
        
        // Receive index to evaluate
        MPI_Recv ( &index ,
                   1 ,
                   MPI_INT ,
                   0,
                   TAG_INDEX,
			       MPI_COMM_WORLD,
                   & status );
        
        D ( std::cout << "Evaluator with rank #" << rank << " received index #" << index << std::endl;);
        if ( index < 0 )
        {
            D ( std::cout << " Stop signal " << std::endl; );
            break;
        }
        
        double * x = new double[n];
        MPI_Recv ( x ,
                   n ,
                   MPI_DOUBLE ,
                   0,
                   TAG_INPUT,
			       MPI_COMM_WORLD,
			       & status );
        
        double c1 = 0.0 , c2 = 0.0;
        for ( int i = 0 ; i < n ; i++ )
        {
            D (std::cout << x[i] << " ";);
            c1 += pow(x[i]-1,2);
            c2 += pow(x[i]+1,2);
        }
        D ( std::cout << std::endl; );
        double * out = new double[m];
        out[0] =  x[4];
        out[1] =  c1 - 25;
        out[2] =  25 - c2;
        delete[] x;

        D ( std::cout << "Evaluator with rank #" << rank << " obtained f=" << out[0] << " c0 = " << out[1] << " c1= " << out[2] << std::endl; );
         
        // Send signal to evaluator server to inform evaluator (#) is ready to send eval
        MPI_Send ( &SIGNAL_DONE_EVAL ,
                   1 ,
                   MPI_CHAR,
                   0 ,
                   TAG_SIGNAL, 
			       MPI_COMM_WORLD);
    
        D ( std::cout << "Evaluator with rank #" << rank << " about to send results for index=" << index << std::endl; );
        
        // Broadcast that the evaluation #index is done and ready to send output values
        MPI_Send ( &index , 1 , MPI_INT , 0 , TAG_INDEX, MPI_COMM_WORLD);

        // Output values
        MPI_Send ( out , m , MPI_DOUBLE , 0 , TAG_OUTPUT, MPI_COMM_WORLD);

        delete[] out;
    }
    
    D ( std::cout << "Evaluator with rank#" << rank << " is done." << std::endl; );
    MPI_Finalize();
    
}

void evalServerStart(const char* input_file_name)
{
    int rank;
    MPI_Comm_rank ( MPI_COMM_WORLD, &rank );
    if (rank != 0)
        return;
    
    D ( std::cout << "Server has rank " << rank << std::endl; );
    
    // First process in charge of reading input file
    ifstream in ( input_file_name );
    
    // Get the input vectors
    std::string str;
    std::vector<std::vector<double>> X;
            
    while (std::getline(in, str))
    {
        std::vector<double> x(n);
        std::istringstream iss(str);
        try
        {
            for (int i = 0; i < n; i++)
            {
                iss >> x[i];
            }
        }
        catch (...)
        {
            cerr << "Invalid value in input file"<< endl;
            MPI_Finalize();
            return;
        }
        X.push_back(x);
    }
    in.close();
    
    char         signal;
    MPI_Status  status;
    int          evaluator;
    std::vector<std::vector<double>> all_outputs(X.size(),std::vector<double>(m,0));
    size_t eval_counter_sent = 0, eval_counter_received = 0, nb_points = X.size() ;
    int index;
    while (eval_counter_received != nb_points)
    {
        D ( std::cout << "EvalServer, eval_counter_sent=" << eval_counter_sent << std::endl;)
        D ( std::cout << "EvalServer, eval_counter_received=" << eval_counter_received << std::endl;)
        
        // Capture signals from evaluator
        MPI_Recv ( &signal         ,
                   1               ,
                   MPI_CHAR       ,
                   MPI_ANY_SOURCE ,
                   MPI_ANY_TAG    ,
			       MPI_COMM_WORLD,
                   &status            );
        
        evaluator = status.MPI_SOURCE;
        
        
        // Case evaluator ready to perform evaluation.
        if (signal == SIGNAL_READY_TO_EVAL && eval_counter_sent<nb_points)
        {
            D ( std::cout << "Evaluator #" << evaluator << " is ready" << std::endl;);

            // send evaluation index to evaluator
            MPI_Send ( &eval_counter_sent  , 1 , MPI_INT , evaluator , TAG_INDEX, MPI_COMM_WORLD );
                        
            std::vector<double> singleX = X.front();
            double * x = &(singleX[0]);
        
            D ( std::cout << "Eval Server is about to send x" << std::endl; )
            
            // send point to the ready evaluator
            MPI_Send ( x  , n , MPI_DOUBLE , evaluator , TAG_INPUT, MPI_COMM_WORLD );
            
            X.erase(X.begin());
            eval_counter_sent ++;
            
        }
        // Case evaluation results are available
        else if ( signal == SIGNAL_DONE_EVAL)
        {
            
            D ( std::cout << "Eval server received signal done eval" << " from evaluator #" << evaluator << std::endl; );
            
            
            // point index + m values
            MPI_Recv ( &index ,
                       1               ,
                       MPI_INT       ,
                       evaluator ,
                       TAG_INDEX    ,
                       MPI_COMM_WORLD,
                       &status            );
            
            D ( std::cout << "Eval server received index " << index << " from evaluator #" << evaluator << std::endl; );
            
            if (index < 0)
            {
                D ( cerr << "Invalid index value received by eval server. Index=" << index << " sent by process #" << evaluator << endl; );
                MPI_Finalize();
                return;
            }
            
            double * out = new double [m];
            MPI_Recv ( out ,
                       m               ,
                       MPI_DOUBLE       ,
                       evaluator ,
                       TAG_OUTPUT    ,
                       MPI_COMM_WORLD,
                       &status            );
            
            D ( std::cout << "Eval server received output from evaluator #" << evaluator << std::endl; );
            
            std::vector<double> single_out(out,out+m);
            all_outputs[index] = single_out;
            
            delete [] out;
            eval_counter_received++;
            
        }
    }
    
    // Broadcast stop to all evaluators
    int np;
    index = -1;
    MPI_Comm_size ( MPI_COMM_WORLD, &np   );
    int i;
    for (i = 0; i < np; i++)
    {
      if (i != 0)
      {
        MPI_Send(&index, 1, MPI_INT, i, TAG_INDEX, MPI_COMM_WORLD);
      }
    }
    i=0;
    for (const auto & single_out: all_outputs )
    {
        D ( std::cout << "Point #" << i++ << " out=" ; );
        // Blackbox output done here.
        for( const auto & val: single_out)
        {
            std::cout << val << " ";
        }
        std::cout << std::endl;
    }
    MPI_Finalize();
    
}





int main ( int argc , char ** argv )
{

    // MPI initialization:
    MPI_Init ( &argc , &argv );
    int rank , np;
    MPI_Comm_rank ( MPI_COMM_WORLD, &rank );
    MPI_Comm_size ( MPI_COMM_WORLD, &np   );

    // check the arguments and the number of processes:
    if ( np <= 1 || argc != 2 ) {
      if ( rank==0 )
        cerr << "usage: mpirun -np p " << argv[0]
         << " input_file, with p>1" << endl;
      MPI_Finalize();
      return 1;
    }
    
    if ( rank == 0 )
        evalServerStart(argv[1]);
    else
        evaluatorStart();
    
    // MPI finalize done within functions
    
    return 0;
}	


