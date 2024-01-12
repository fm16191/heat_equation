#include <cmath>
#include <fstream>
#include <getopt.h>
#include <iomanip>
#include <iostream>
#include <vector>

double NB_X = 80;
double NB_T = 1000;

double DX = 1 / (NB_X + 1);
double D = 2.3 * 10e-5;

double DT;

std::string out_filename = "output.txt";

using std::cout, std::cerr;
using std::setw, std::setprecision;
using std::stod;
using std::vector;

int print_usage(char *exec)
{
    printf("Heat Equation 1D\n");
    printf("Usage : %s [-xtdoh]\n", exec);
    printf("\n");
    printf("Options : \n"
           " -x Set the number of spatial grid points. Default : %.0f\n"
           " -t Set the number of temporal grid points. Default : %.0f\n"
           " -d Set the thermal diffusivity coefficient. Default : %.2f\n"
           " -o Set the output filename. Default : %s\n"
           "\n"
           " -h, --help Show this message and exit\n",
           NB_X, NB_T, D, out_filename.c_str());
    return 0;
}

int main(int argc, char *argv[])
{
    const char *short_options = "hx:t:d:o:";
    const struct option long_options[] = { { "spatial_points", required_argument, 0, 'x' },
                                           { "temporal_points", required_argument, 0, 't' },
                                           { "thermal_diffusivity_coefficient", required_argument,
                                             0, 'd' },
                                           { "output_filename", required_argument, nullptr, 'o' },
                                           { "help", no_argument, nullptr, 'h' },
                                           { nullptr, 0, nullptr, 0 } };

    int option;
    while ((option = getopt_long_only(argc, argv, short_options, long_options, nullptr)) != -1) {
        switch (option) {
            case 'h':
                return print_usage(argv[0]);
                break;
            case 'x':
                NB_X = stod(optarg);
                break;
            case 't':
                NB_T = stod(optarg);
                break;
            case 'd':
                D = stod(optarg);
                break;
            case 'o':
                out_filename = optarg;
                break;
            default:
                cerr << "Unknown option\n";
                return 1;
        }
    }

    // Ensure stability
    DT = (DX * DX) / D / 2; // Ensure stability
    DT = DT / 2;            // divide by 2 again to see some results

    // Verify stability
    double r = D * DT / (DX * DX);
    if (r > 0.5) {
        cerr << "Stability condition is not met : " << r << "\n";
        return 1;
    }

    cerr << "Configuration : \n";
    cerr << "  Spatial Points (nb_x)     : " << NB_X << "\n";
    cerr << "  Temporal Points (nb_t)    : " << NB_T << "\n";
    cerr << "  Thermal diffusivity (d)   : " << D << "\n";
    cerr << "  Step in x (dx)            : " << DX << "\n";
    cerr << "  Step in t (dt)            : " << DT << "\n";

    cerr << "  Stability indicator       : " << r << "\n";

    vector<vector<double>> u(NB_T, vector<double>(NB_X + 2, 0));

    // Boundary conditions
    double ux0 = 0;
    double uxa = 0;
    for (size_t j = 0; j < NB_T; ++j) {
        u[j][0] = ux0;
        u[j][NB_X + 1] = uxa;
    }

    // Initial conditions
    for (size_t i = 0; i < NB_X; ++i)
        u[0][i + 1] = (NB_X * i - (i * i)) / 10;

    // Iterate
    for (size_t j = 0; j < NB_T - 1; ++j) {
        for (size_t i = 1; i < NB_X + 1; ++i) {
            u[j + 1][i] = r * u[j][i - 1] + (1 - 2 * r) * u[j][i] + r * u[j][i + 1];
        }
    }

    // Output
    std::ofstream out_file(out_filename);
    if (!out_file.is_open()) {
        cerr << "Error opening output file\n";
        return 1;
    }

    for (size_t j = 0; j < NB_T; ++j) {
        for (size_t i = 0; i < NB_X + 2; ++i) {
            out_file << setw(15) << setprecision(3) << u[j][i];
        }
        out_file << "\n";
    }

    cout << "Results written to " << out_filename << "\n";

    return 0;
}
