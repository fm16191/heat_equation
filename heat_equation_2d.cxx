#include <cmath>
#include <fstream>
#include <getopt.h>
#include <iomanip>
#include <iostream>
#include <vector>

double NB_X = 1000;
double NB_Y = 1000;
double NB_T = 1000;

double DX = 1 / (NB_X + 1);
double DY = 1 / (NB_Y + 1);
double D = 2.3 * 10e-5;

double DT; // defined afterwards

double DX2 = DX * DX;
double DY2 = DY * DY;

std::string out_filename = "output.txt";

using std::cout, std::cerr;
using std::min, std::sqrt;
using std::setw, std::setprecision;
using std::stod;
using std::vector;

int print_usage(char *exec)
{
    printf("Heat Equation 2D\n");
    printf("Usage : %s [-xytdoh]\n", exec);
    printf("\n");
    printf("Options : \n"
           " -x Set the number of spatial grid points in the X axis. Default : %.0f\n"
           " -y Set the number of spatial grid points in the Y axis. Default : %.0f\n"
           " -t Set the number of temporal grid points. Default : %.0f\n"
           " -d Set the thermal diffusivity coefficient. Default : %.2f\n"
           " -o Set the output filename. Default : %s\n"
           "\n"
           " -h, --help Show this message and exit\n",
           NB_X, NB_Y, NB_T, D, out_filename.c_str());
    return 0;
}

int main(int argc, char *argv[])
{
    const char *short_options = "hx:y:t:d:o:";
    const struct option long_options[] = { { "spatial_x_points", required_argument, 0, 'x' },
                                           { "spatial_y_points", required_argument, 0, 'y' },
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
            case 'y':
                NB_Y = stod(optarg);
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
    double DT = (DX2 * DY2) / (2 * D * (DX2 + DY2)) * 0.5;

    // Verify stability
    // α = D △t/△x^2, β = D △t/△y^2, α + β < 1/2
    double stability = D * DT * (DX2 + DY2) / (DX2 * DY2);
    if (stability > 0.5) {
        cerr << "Stability condition is not met : " << stability << "\n";
        return 1;
    }

    cerr << "Configuration : \n";
    cerr << "  Spatial Points (nb_x)     : " << NB_X << "\n";
    cerr << "  Spatial Points (nb_y)     : " << NB_Y << "\n";
    cerr << "  Temporal Points (nb_t)    : " << NB_T << "\n";
    cerr << "  Thermal diffusivity (d)   : " << D << "\n";
    cerr << "  Step in x (dx)            : " << DX << "\n";
    cerr << "  Step in y (dy)            : " << DY << "\n";
    cerr << "  Step in t (dt)            : " << DT << "\n";

    cerr << "  Stability indicator       : " << stability << "\n";

    vector<vector<double>> u(NB_X + 2, vector<double>(NB_Y + 2, 0));

    // Boundary conditions
    double ui0 = 40;
    double uia = 40;
    double uj0 = 40;
    double uja = 40;
    for (size_t i = 0; i < NB_X + 2; ++i) {
        u[i][0] = ui0;
        u[i][NB_Y + 1] = uia;
    }
    for (size_t j = 0; j < NB_Y + 2; ++j) {
        u[0][j] = uj0;
        u[NB_X + 1][j] = uja;
    }

    // Init condition : a 100 °C circle 1/4 the size of the simulation at its center
    double center_x = (NB_X + 1) / 2.0;
    double center_y = (NB_Y + 1) / 2.0;
    double radius = min(center_x, center_y) / 4;

    for (size_t i = 0; i < NB_X + 2; ++i) {
        for (size_t j = 0; j < NB_Y + 2; ++j) {
            double distance =
                sqrt((i - center_x) * (i - center_x) + (j - center_y) * (j - center_y));
            if (distance <= radius) {
                u[i][j] = 100.0;
            }
        }
    }

    // Deep copy vector
    vector<vector<double>> v(u);

    // Iterate
    double a = D * DT / DX2;
    double b = D * DT / DY2;
    double c = (1 - 2 * a - 2 * b);
    for (size_t t = 0; t < NB_T; ++t) {
        for (size_t j = 1; j < NB_Y + 1; ++j) {
            for (size_t i = 1; i < NB_X + 1; ++i) {
                v[i][j] =
                    a * (u[i + 1][j] + u[i - 1][j]) + b * (u[i][j + 1] + u[i][j - 1]) + c * u[i][j];
            }
        }
        u.swap(v);
    }

    // Output
    std::ofstream out_file(out_filename);
    if (!out_file.is_open()) {
        cerr << "Error opening output file\n";
        return 1;
    }

    for (size_t i = 0; i < NB_X + 2; ++i) {
        for (size_t j = 0; j < NB_Y + 2; ++j) {
            out_file << setw(15) << setprecision(3) << u[i][j];
        }
        out_file << "\n";
    }

    cout << "Results written to " << out_filename << "\n";

    return 0;
}
