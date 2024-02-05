#include <cmath>
#include <fstream>
#include <getopt.h>
#include <iomanip>
#include <iostream>
#include <vector>

double NB_X = 100;
double NB_Y = 100;
double NB_Z = 100;
double NB_T = 100;

double D = 2.3e-5;
double DT = 0.001;

int write_interval = 0;

std::string out_filename = "output.txt";

using std::cout, std::cerr;
using std::min, std::sqrt;
using std::setw, std::setprecision;
using std::stod;
using std::vector;

void write_results(vector<vector<vector<double>>> u, size_t step)
{
    // Output
    std::ofstream out_file;
    out_file.open(out_filename, std::ios::app);

    if (!out_file.is_open()) {
        cerr << "Error opening output file\n";
        exit(1);
    }

    size_t len_x = u.size() - 2;
    size_t len_y = u[0].size() - 2;
    size_t len_z = u[0][0].size() - 2;

    out_file << len_x << " " << len_y << " " << len_z;

    for (size_t i = 1; i < len_x + 1; ++i) {
        for (size_t j = 1; j < len_y + 1; ++j) {
            for (size_t k = 1; k < len_z + 1; ++k) {
                out_file << " " << setprecision(3) << u[i][j][k];
            }
        }
    }
    out_file << "\n";

    printf("t=%ld written to %s\n", step, out_filename.c_str());
}

int print_usage(char *exec)
{
    printf("Heat Equation 2D\n");
    printf("Usage : %s [-xyztdowh]\n", exec);
    printf("\n");
    printf("Options : \n"
           " -x Set the number of spatial grid points in the X axis. Default : %.0f\n"
           " -y Set the number of spatial grid points in the Y axis. Default : %.0f\n"
           " -z Set the number of spatial grid points in the Z axis. Default : %.0f\n"
           " -t Set the number of temporal grid points. Default : %.0f\n"
           " -d Set the thermal diffusivity coefficient. Default : %.2f\n"
           " -o Set the output filename. Default : %s\n"
           " -w Set the interval between each data write. Default : only write final result"
           "\n"
           " -h, --help Show this message and exit\n",
           NB_X, NB_Y, NB_Z, NB_T, D, out_filename.c_str());
    return 0;
}

double compute_total_temp(vector<vector<vector<double>>> u)
{
    double sum = 0.0;
    for (size_t i = 1; i < NB_X + 1; ++i) {
        for (size_t j = 1; j < NB_Y + 1; ++j) {
            for (size_t k = 1; k < NB_Z + 1; ++k) {
                sum += u[i][j][k];
            }
        }
    }
    return sum;
}

void ensure_mass_conservation(double initial_temp, vector<vector<vector<double>>> u, size_t t)
{
    double total_temp = compute_total_temp(u);
    double mass_change = fabs(initial_temp - total_temp);
    printf("t=%ld Total temperature : %.7e, change in mass: %.1e (should be close to 0)\n", t,
           total_temp, mass_change);
}

int main(int argc, char *argv[])
{
    const char *short_options = "hx:y:z:t:d:o:w:";
    const struct option long_options[] = { { "spatial_x_points", required_argument, 0, 'x' },
                                           { "spatial_y_points", required_argument, 0, 'y' },
                                           { "spatial_z_points", required_argument, 0, 'z' },
                                           { "temporal_points", required_argument, 0, 't' },
                                           { "thermal_diffusivity_coefficient", required_argument,
                                             0, 'd' },
                                           { "output_filename", required_argument, nullptr, 'o' },
                                           { "write_interval", required_argument, nullptr, 'w' },
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
            case 'z':
                NB_Z = stod(optarg);
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
            case 'w':
                write_interval = stod(optarg);
                break;
            default:
                cerr << "Unknown option\n";
                return 1;
        }
    }

    // Clear output
    std::ofstream out_file(out_filename);
    if (!out_file.is_open()) {
        cerr << "Error opening output file\n";
        exit(1);
    }

    double DX = 1 / (NB_X + 1);
    double DY = 1 / (NB_Y + 1);
    double DZ = 1 / (NB_Z + 1);
    double DX2 = DX * DX;
    double DY2 = DY * DY;
    double DZ2 = DZ * DZ;

    // Ensure stability
    // double DT = (DX2 * DY2) / (2 * D * (DX2 + DY2)) * 0.5;

    // Verify stability
    // α = D △t/△x^2, β = D △t/△y^2, α + β < 1/2
    double stability = D * DT * (DX2 + DY2) / (DX2 * DY2);
    if (stability > 0.5) {
        cerr << "Stability condition is not met : " << stability << "\n";
        return 1;
    }

    cerr << "Configuration : \n";
    cerr << "  Spatial Points X axis (nb_x)  : " << NB_X << "\n";
    cerr << "  Spatial Points Y axis (nb_y)  : " << NB_Y << "\n";
    cerr << "  Spatial Points Z axis (nb_z)  : " << NB_Z << "\n";
    cerr << "  Temporal Points (nb_t)        : " << NB_T << "\n";
    cerr << "  Thermal diffusivity (d)       : " << D << "\n";
    cerr << "  Step in x (dx)                : " << DX << "\n";
    cerr << "  Step in y (dy)                : " << DY << "\n";
    cerr << "  Step in z (dz)                : " << DZ << "\n";
    cerr << "  Step in t (dt)                : " << DT << "\n";

    if (write_interval)
        cerr << "  Write interval                : " << write_interval << "\n";

    cerr << "  Stability indicator           : " << stability << "\n\n";

    vector<vector<vector<double>>> u(NB_X + 2,
                                     vector<vector<double>>(NB_Y + 2, vector<double>(NB_Z + 2, 0)));

    // // Boundary conditions
    // double ux0 = 40;
    // double uxa = 0;
    // double uy0 = 40;
    // double uya = 0;
    // double uz0 = 40;
    // double uza = 0;

    // for (size_t j = 0; j < NB_Y + 2; ++j) {
    //     for (size_t k = 0; k < NB_Y + 2; ++k) {
    //         u[0][j][k] = ux0;
    //         u[NB_X + 1][j][k] = uxa;
    //     }
    // }

    // for (size_t i = 0; i < NB_X + 2; ++i) {
    //     for (size_t k = 0; k < NB_Z + 2; ++k) {
    //         u[i][0][k] = uy0;
    //         u[i][NB_Y + 1][k] = uya;
    //     }
    // }

    // for (size_t i = 0; i < NB_X + 2; ++i) {
    //     for (size_t j = 0; j < NB_Y + 2; ++j) {
    //         u[i][j][0] = uz0;
    //         u[i][j][NB_Z + 1] = uza;
    //     }
    // }

    // Init condition : a 100 °C sphere 2/3 the size of the simulation at its center
    double center_x = (NB_X + 1) / 2.0;
    double center_y = (NB_Y + 1) / 2.0;
    double center_z = (NB_Z + 1) / 2.0;
    double radius = min(center_x, min(center_y, center_z)) / 1.5;

    for (size_t i = 0; i < NB_X + 2; ++i) {
        for (size_t j = 0; j < NB_Y + 2; ++j) {
            for (size_t k = 0; k < NB_Z + 2; ++k) {
                double distance =
                    sqrt((i - center_x) * (i - center_x) + (j - center_y) * (j - center_y) +
                         (k - center_z) * (k - center_z));
                if (distance <= radius) {
                    u[i][j][k] = 100.0;
                }
            }
        }
    }

    vector<vector<vector<double>>> u_next(u);

    double initial_temp = compute_total_temp(u);
    printf("Initial total temperature       : %.7e\n", initial_temp);

    // Iterate
    double a = D * DT / DX2;
    double b = D * DT / DY2;
    double c = D * DT / DZ2;
    double d = (1 - 2 * a - 2 * b - 2 * c);
    for (size_t t = 0; t < NB_T; ++t) {
        // Copy from one plan to the other
        // X axis
        for (size_t j = 0; j < NB_Y + 2; ++j) {
            for (size_t k = 0; k < NB_Z + 2; ++k) {
                u[0][j][k] = u[NB_X][j][k];
                u[NB_X + 1][j][k] = u[1][j][k];
            }
        }
        // Y axis
        for (size_t i = 0; i < NB_X + 2; ++i) {
            for (size_t k = 0; k < NB_Z + 2; ++k) {
                u[i][0][k] = u[i][NB_Y][k];
                u[i][NB_Y + 1][k] = u[i][1][k];
            }
        }
        // Z axis
        for (size_t i = 0; i < NB_X + 2; ++i) {
            for (size_t j = 0; j < NB_Y + 2; ++j) {
                u[i][j][0] = u[i][j][NB_Z];
                u[i][j][NB_Z + 1] = u[i][j][1];
            }
        }

        for (size_t i = 1; i < NB_X + 1; ++i) {
            for (size_t j = 1; j < NB_Y + 1; ++j) {
                for (size_t k = 1; k < NB_Z + 1; ++k) {
                    u_next[i][j][k] = a * (u[i + 1][j][k] + u[i - 1][j][k]) +
                                      b * (u[i][j + 1][k] + u[i][j - 1][k]) +
                                      c * (u[i][j][k + 1] + u[i][j][k - 1]) + d * u[i][j][k];
                }
            }
        }
        u.swap(u_next);
        if (write_interval && t % write_interval == 0 && t != 0) {
            ensure_mass_conservation(initial_temp, u, t);
            write_results(u, t);
        }
    }

    ensure_mass_conservation(initial_temp, u, (size_t)NB_T);
    write_results(u, NB_T);

    return 0;
}
