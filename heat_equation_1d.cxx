#include <cmath>
#include <fstream>
#include <getopt.h>
#include <iomanip>
#include <iostream>
#include <vector>

// Defaults
double NB_X = 80.0;
double NB_T = 5000.0;
double D = 2.3 * 1e-5;
double DT = 0.1;


bool USE_FINITE_VOLUME_METHOD = false;

std::string out_filename = "output.txt";

using std::cout, std::cerr;
using std::setprecision;
using std::stod;
using std::vector;

/* Helper function : print usage */
int print_usage(char *exec)
{
    printf("Heat Equation 1D\n");
    printf("Usage : %s [-vxtdoh]\n", exec);
    printf("\n");
    printf("Options : \n"
           " -v Use Finite Volume Method. Default : Finite Element Method.\n"
           " -x Set the number of spatial grid points. Default : %.0f\n"
           " -t Set the number of temporal grid points. Default : %.0f\n"
           " -d Set the thermal diffusivity coefficient. Default : %.2e\n"
           " -o Set the output filename. Default : %s\n"
           "\n"
           " -h, --help Show this message and exit\n",
           NB_X, NB_T, D, out_filename.c_str());
    return 0;
}

/* compute the total temperature in the simulation
 * u : the current state of the heatmap
 * t : the current time step of simulation
 * sum : total heat of simulation
 */
double compute_total_temp(vector<vector<double>> u, size_t t)
{
    double sum = 0.0;
    for (size_t i = 1; i < NB_X + 1; ++i) {
        sum += u[t][i];
    }
    return sum;
}

// void set_init_conditions(vector<vector<double>> &u){

// }

/* Set initial boundary conditions
 * u : the current state of the heatmap
 */
void set_boundary_conditions(vector<vector<double>> &u)
{
    double ux0 = 0;
    double uxa = 0;
    for (size_t j = 0; j < NB_T; ++j) {
        u[j][0] = ux0;
        u[j][NB_X + 1] = uxa;
    }
}

/* Update the boundaries to ensure the periodicity of the simulation
 * u : the current state of the heatmap
 * t : the current time step of simulation
 */
void update_periodic_boundaries(vector<vector<double>> &u, size_t t)
{
    u[t][0] = u[t][NB_X];
    u[t][NB_X + 1] = u[t][1];
}

/* Update the boundaries to ensure the periodicity of the simulation
 * u : the current state of the heatmap
 * t : the current time step of simulation
 * i : indice of the cell to be computed
 * orient_x : direction of the flux to be computed
 */
double F(vector<vector<double>> const &u, size_t t, int i, int orient_x)
{
    /*** Naïve flux compution */
    // if (orient_x > 0)
    //     return (u[t][i + 1] - u[t][i]);
    // else if (orient_x < 0)
    //     return (u[t][i] - u[t][i - 1]);

    // // case should never happend
    // exit(0);

    /*** Removing if statements, function computes 2 times faster. */
    // Instead of doing u[t][i] - u[t][i-1], let's specify i as (i-1) and orient_x as 1 instead of
    // -1. This will give us u[t][(i-1)+1] - u[t][(i-1)], which is what we expected.
    return (u[t][i + orient_x] - u[t][i]);
}

int main(int argc, char *argv[])
{
    const char *short_options = "hvx:t:d:o:";
    const struct option long_options[] = {
        { "use_finite_volume_method", no_argument, nullptr, 'v' },
        { "spatial_points", required_argument, 0, 'x' },
        { "temporal_points", required_argument, 0, 't' },
        { "thermal_diffusivity_coefficient", required_argument, 0, 'd' },
        { "output_filename", required_argument, nullptr, 'o' },
        { "help", no_argument, nullptr, 'h' },
        { nullptr, 0, nullptr, 0 }
    };

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
            case 'v':
                USE_FINITE_VOLUME_METHOD = true;
                break;
            default:
                cerr << "Unknown option\n";
                return 1;
        }
    }

    /* Code */
    double DX = 1 / (NB_X + 1);

    // Ensure stability
    // DT = (DX * DX) / D / 2; // Ensure stability
    // DT = DT / 2;            // divide by 2 to see some results

    // Verify stability
    double r = D * DT / (DX * DX);
    if (r > 0.5) {
        cerr << "Stability condition is not met : " << r << "\n";
        return 1;
    }

    cerr << "Configuration : \n";
    cerr << "  Using Finite " << (USE_FINITE_VOLUME_METHOD == true ? "Volume" : "Element")
         << " method\n";
    cerr << "  Spatial Points (nb_x)     : " << NB_X << "\n";
    cerr << "  Temporal Points (nb_t)    : " << NB_T << "\n";
    cerr << "  Thermal diffusivity (d)   : " << D << "\n";
    cerr << "  Step in x (dx)            : " << DX << "\n";
    cerr << "  Step in t (dt)            : " << DT << "\n";

    cerr << "  Stability indicator       : " << r << "\n\n";

    vector<vector<double>> u(NB_T, vector<double>(NB_X + 2, 0));

    // Initial conditions
    for (size_t i = 0; i < NB_X + 1; ++i)
        u[0][i] = (NB_X * i - (i * i)) / 10;

    // set_boundary_conditions(u); // Disable boundary conditions as we use periodic ones

    double initial_temp = compute_total_temp(u, 0);
    printf("Initial total temperature   : %.7e\n", initial_temp);

    // Iterate
    for (size_t t = 0; t < NB_T - 1; ++t) {
        update_periodic_boundaries(u, t);

        if (USE_FINITE_VOLUME_METHOD) {
            for (size_t i = 1; i < NB_X + 1; ++i) {
                u[t + 1][i] = u[t][i] + r * (F(u, t, i, 1) - F(u, t, i - 1, 1));
            }
        }
        else { // Finite Element Method
            for (size_t i = 1; i < NB_X + 1; ++i) {
                u[t + 1][i] = u[t][i] + r * (u[t][i - 1] + -2.0 * u[t][i] + u[t][i + 1]);
            }
        }
    }

    // Ensure mass conservation
    double total_temp = compute_total_temp(u, NB_T - 1);
    double mass_change = fabs(initial_temp - total_temp);
    printf("Total temperature at t=%.0f : %.7e, change in mass: %.1e (should be close to 0)\n",
           NB_T, total_temp, mass_change);

    // Output
    std::ofstream out_file(out_filename);
    if (!out_file.is_open()) {
        cerr << "Error opening output file\n";
        return 1;
    }

    for (size_t j = 0; j < NB_T; ++j) {
        for (size_t i = 0; i < NB_X + 2; ++i) {
            out_file << setprecision(3) << u[j][i] << " ";
        }
        out_file << "\n";
    }

    cout << "Results written to " << out_filename << "\n";

    return 0;
}
