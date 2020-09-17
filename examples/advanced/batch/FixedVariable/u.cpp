#include <ctime>        // For time_t
#include <fstream>      // For ifstream
#include <iostream>
#include <cmath>        // For sqrt
#include <stdexcept>    // For logic_error
//#include <unistd.h>     // For usleep

const int n = 5;

int main (int argc, char **argv)
{
    bool eval_ok = false;

    // Remotely based on G2.
    double f = 1e+20, g1 = 1e+20, g2 = 1e+20;
    double sum1 = 0.0, sum2 = 0.0, sum3 = 0.0, prod1 = 1.0, prod2 = 1.0;
    double x[n];

    bool x0read = false;
    if (argc >= 2)
    {
        std::string x0file = argv[1];
        std::ifstream in (argv[1]);
        for (int i = 0; i < n; i++)
        {
            if (in.fail())
            {
                std::cerr << "Error reading file " << x0file << " for x0." << std::endl;
                x0read = false;
                break;
            }
            in >> x[i];
            x0read = true;
        }
        in.close();
    }

    if (x0read)
    {
        try
        {
            for (int i = 0; i < n ; i++)
            {
                //std::cerr << "sleep " << i << std::endl;
                //usleep(100000);
                sum1  += pow(cos(x[i]), 4);
                sum2  += x[i];
                sum3  += (i+1)*x[i]*x[i];
                prod1 *= pow(cos(x[i]), 2);
                if (prod2 != 0.0)
                {
                    if (x[i] == 0.0)
                    {
                        prod2 = 0.0;
                    }
                    else
                    {
                        prod2 *= x[i];
                    }
                }
            }

            g1 = -prod2 + 0.75;
            g2 = sum2 -7.5 * n;

            f = 10*g1 + 10*g2;
            //std::cout << "g1 = " << g1 << " g2 = " << g2 << " f (partial) = " << f << std::endl;
            if (0.0 != sum3)
            {
                f -= (sum1 -2*prod1) / std::abs(sqrt(sum3));
            }
            else
            {
                f = NAN;
            }
            //std::cout << "f " << f << std::endl;
            // Scale
            if (!std::isnan(f))
            {
                f *= 1e-5;
            }
            //std::cout << "f (scaled) " << f << std::endl;

            eval_ok = !std::isnan(f);
            //std::cout << "eval_ok: " << eval_ok << std::endl;

            // Simulate some constraints by making f NaN in some cases.
            /*
            if (f < -2000)
            {
                f = NAN;
            }
            */
            /*
            // Simulate random eval errors.
            std::time_t* theTime = new time_t();
            std::time(theTime);
            if (int(*theTime) % 12 == 0)
            {
                eval_ok = false;
            }
            */
        }
        catch (std::exception &e)
        {
            std::string err("Exception: ");
            err += e.what();
            throw std::logic_error(err);
        }
    }

    //std::cout << "End: f = " << f << " eval_ok = " << eval_ok << std::endl;
    std::cout << - f - 2000 << " " << f << std::endl;

    // Return 0 if eval_ok.
    // Note: returning a value other than 0 will make this evaluation
    // an EVAL_ERROR, meaning the point could be re-evaluated.
    // If eval_ok is true and f is NaN, then this evaluation is an
    // EVAL_FAILED, and the point will not be re-evaluated.
    return !eval_ok;
}
