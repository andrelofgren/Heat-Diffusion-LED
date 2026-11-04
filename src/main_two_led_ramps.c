#include <stdio.h>

#include <stddef.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <signal.h>
#include <pwd.h>

#include <cjson/cJSON.h>
#include <gsl/gsl_spmatrix.h>
#include <gsl/gsl_blas.h>

#include "../include/heat_problem_1d.h"
#include "../include/fem_function_1d.h"
#include "../include/serial_port.h"

typedef struct MaterialParams {
    double alpha;
    double beta;
    double gamma;
} MaterialParams;

typedef struct SourceParams {
    double A0;
    double sigma;
    double mu;
    double tau;
} SourceParams;

typedef struct DomainParams {
    double x0;
    double x1;
    double t0;
    double t1;
} DomainParams;

typedef struct ComputationalParams {
    size_t nx;
    size_t nt;
    size_t deg;
    size_t num_steps_update;
    size_t num_leds_list_1;
    size_t num_leds_list_2;
    size_t repeat_num_times;

} ComputationalParams;

typedef struct ColormapParams {
    size_t num_bins;
    double *red_fem_func_vals_data;
    double *green_fem_func_vals_data;
    double *blue_fem_func_vals_data;
} ColormapParams;

typedef struct Params {
    SerialPort serial_port;
    MaterialParams material;
    SourceParams source;
    DomainParams domain;
    ComputationalParams computational;
    ColormapParams colormap;
} Params;

static Params params;

double source_amplitude(const double t)
{
    const double A0 = params.source.A0;
    const double tau = params.source.tau;
    return A0*exp(-t/tau);
}

double heat_source(const double x, const double t)
{
    const double mu = params.source.mu;
    const double sigma = params.source.sigma;
    const double X = sqrt(0.5)*((x-mu)/sigma);
    return source_amplitude(t)*exp(-X*X)/(sigma*sqrt(2*M_PI));
}

void check_json_object(cJSON *json_root, const cJSON *json_obj, const char *key)
{
    if (!json_obj) {
        fprintf(stderr, "cJSON error: no key '%s' found\n", key);
        cJSON_Delete(json_root);
        exit(EXIT_FAILURE);
    }
}

void read_config(const char *config_path)
{
    FILE *fptr = fopen(config_path, "r");
    if (!fptr) {
        fprintf(stderr, "Failed to open %s: %s \n", config_path, strerror(errno));
        exit(EXIT_FAILURE);
    }

    char json_str_buffer[4096];
    fread(json_str_buffer, sizeof(char), sizeof json_str_buffer, fptr);

    fclose(fptr);

    cJSON *param_json = cJSON_ParseWithLength(json_str_buffer, sizeof json_str_buffer);
    if (!param_json) {
        const char *error_ptr = cJSON_GetErrorPtr();
        if (error_ptr)
            fprintf(stderr, "JSON formating error in %s before: %s\n", config_path, error_ptr);
        exit(EXIT_FAILURE);
    }

    cJSON *serial_params_json = cJSON_GetObjectItemCaseSensitive(param_json, "serial");
    check_json_object(param_json, serial_params_json, "serial_port");

    const cJSON *json_val;
    json_val = cJSON_GetObjectItemCaseSensitive(serial_params_json, "device");
    check_json_object(param_json, json_val, "device");
    strncpy(params.serial_port.device, json_val->valuestring, sizeof json_str_buffer);

    json_val = cJSON_GetObjectItemCaseSensitive(serial_params_json, "baud_rate");
    check_json_object(param_json, json_val, "baud_rate");
    params.serial_port.baud_rate = json_val->valueint;

    cJSON *domain_params_json = cJSON_GetObjectItemCaseSensitive(param_json, "domain");
    check_json_object(param_json, serial_params_json, "domain");

    json_val = cJSON_GetObjectItemCaseSensitive(domain_params_json, "x_left");
    check_json_object(param_json, json_val, "x_left");
    params.domain.x0 = json_val->valuedouble;

    json_val = cJSON_GetObjectItemCaseSensitive(domain_params_json, "x_right");
    check_json_object(param_json, json_val, "x_right");
    params.domain.x1 = json_val->valuedouble;

    json_val = cJSON_GetObjectItemCaseSensitive(domain_params_json, "time_start");
    check_json_object(param_json, json_val, "time_start");
    params.domain.t0 = json_val->valuedouble;

    json_val = cJSON_GetObjectItemCaseSensitive(domain_params_json, "time_end");
    check_json_object(param_json, json_val, "time_end");
    params.domain.t1 = json_val->valuedouble;

    cJSON *material_params_json = cJSON_GetObjectItemCaseSensitive(param_json, "material");
    check_json_object(param_json, material_params_json, "material");

    json_val = cJSON_GetObjectItemCaseSensitive(material_params_json, "diffusivity");
    check_json_object(param_json, json_val, "diffusivity");
    params.material.alpha = json_val->valuedouble;

    json_val = cJSON_GetObjectItemCaseSensitive(material_params_json, "advection_velocity");
    check_json_object(param_json, json_val, "advection_velocity");
    params.material.beta = json_val->valuedouble;

    json_val = cJSON_GetObjectItemCaseSensitive(material_params_json, "heat_transfer_coefficient");
    check_json_object(param_json, json_val, "heat_transfer_coefficient");
    params.material.gamma = json_val->valuedouble;

    cJSON *source_params_json = cJSON_GetObjectItemCaseSensitive(param_json, "source");
    check_json_object(param_json, source_params_json, "source");

    json_val = cJSON_GetObjectItemCaseSensitive(source_params_json, "amplitude");
    check_json_object(param_json, json_val, "amplitude");
    params.source.A0 = json_val->valuedouble;

    json_val = cJSON_GetObjectItemCaseSensitive(source_params_json, "std_dev");
    check_json_object(param_json, json_val, "std_dev");
    params.source.sigma = json_val->valuedouble;

    json_val = cJSON_GetObjectItemCaseSensitive(source_params_json, "xpos");
    check_json_object(param_json, json_val, "xpos");
    params.source.mu = json_val->valuedouble;

    json_val = cJSON_GetObjectItemCaseSensitive(source_params_json, "relaxation");
    check_json_object(param_json, json_val, "relaxation");
    params.source.tau = json_val->valuedouble;

    cJSON *comput_params_json = cJSON_GetObjectItemCaseSensitive(
        param_json, "computational"
    );
    check_json_object(param_json, source_params_json, "computational");

    json_val = cJSON_GetObjectItemCaseSensitive(comput_params_json, "num_elements");
    check_json_object(param_json, json_val, "num_elements");
    params.computational.nx = json_val->valueint;

    json_val = cJSON_GetObjectItemCaseSensitive(comput_params_json, "num_time_steps");
    check_json_object(param_json, json_val, "num_time_steps");
    params.computational.nt = json_val->valueint;

    json_val = cJSON_GetObjectItemCaseSensitive(comput_params_json, "lagrange_poly_deg");
    check_json_object(param_json, json_val, "lagrange_poly_deg");
    params.computational.deg = json_val->valueint;

    json_val = cJSON_GetObjectItemCaseSensitive(comput_params_json, "update_after_num_steps");
    check_json_object(param_json, json_val, "num_steps_update");
    params.computational.num_steps_update = json_val->valueint;

    json_val = cJSON_GetObjectItemCaseSensitive(comput_params_json, "repeat_num_times");
    check_json_object(param_json, json_val, "repeat_num_times");
    params.computational.repeat_num_times = json_val->valueint;

    json_val = cJSON_GetObjectItemCaseSensitive(comput_params_json, "num_leds_list_1");
    check_json_object(param_json, json_val, "num_leds_list_1");
    params.computational.num_leds_list_1 = json_val->valueint;

    json_val = cJSON_GetObjectItemCaseSensitive(comput_params_json, "num_leds_list_2");
    check_json_object(param_json, json_val, "num_leds_list_2");
    params.computational.num_leds_list_2 = json_val->valueint;

    cJSON *cmap_params_json = cJSON_GetObjectItemCaseSensitive(param_json, "colormap");
    check_json_object(param_json, source_params_json, "colormap");

    json_val = cJSON_GetObjectItemCaseSensitive(cmap_params_json, "num_bins");
    check_json_object(param_json, json_val, "num_bins");
    params.colormap.num_bins = json_val->valueint;

    params.colormap.red_fem_func_vals_data = malloc(
        params.colormap.num_bins * sizeof(double)
    );
    params.colormap.green_fem_func_vals_data = malloc(
        params.colormap.num_bins * sizeof(double)
    );
    params.colormap.blue_fem_func_vals_data = malloc(
        params.colormap.num_bins * sizeof(double)
    );

    char color_key[128];
    for (size_t i = 0; i < params.colormap.num_bins; i++) {
        snprintf(color_key, sizeof color_key, "color_%zu", i+1);

        const cJSON *json_array = cJSON_GetObjectItemCaseSensitive(
            cmap_params_json, color_key
        );
        check_json_object(param_json, json_array, color_key);

        const cJSON *json_array_element_red_value = cJSON_GetArrayItem(json_array, 0);
        check_json_object(param_json, json_array_element_red_value, "0");
        params.colormap.red_fem_func_vals_data[i] = json_array_element_red_value->valuedouble;

        const cJSON *json_array_element_green_value = cJSON_GetArrayItem(json_array, 1);
        check_json_object(param_json, json_array_element_green_value, "1");
        params.colormap.green_fem_func_vals_data[i] = json_array_element_green_value->valuedouble;

        const cJSON *json_array_element_blue_value = cJSON_GetArrayItem(json_array, 2);
        check_json_object(param_json, json_array_element_green_value, "2");
        params.colormap.blue_fem_func_vals_data[i] = json_array_element_blue_value->valuedouble;
    }

    cJSON_Delete(param_json);
}

double calibrate(HeatProblem1D *hp, const size_t nt)
{
    HeatProblem1D_setup(hp);

    FEMFunction1D u_fem_func;
    FEMFunction1D_init(&u_fem_func, hp->mesh);
    FEMFunction1D_set_vals_from_vec(&u_fem_func, hp->sol_vec);
    gsl_vector *u_fem_vals = gsl_vector_calloc(hp->mesh->pvec->size);
    FEMFunction1D_eval(&u_fem_func, hp->mesh->pvec, u_fem_vals);

    gsl_vector *u_min_vec = gsl_vector_calloc(nt+1);
    gsl_vector *u_max_vec = gsl_vector_calloc(nt+1);
    for (size_t k = 0; k < nt; k++) {
        HeatProblem1D_step(hp);
        HeatProblem1D_reset_rhs(hp);
        FEMFunction1D_set_vals_from_vec(&u_fem_func, hp->sol_vec);
        FEMFunction1D_eval(&u_fem_func, hp->mesh->pvec, u_fem_vals);
        gsl_vector_minmax(u_fem_vals, &u_min_vec->data[k+1], &u_max_vec->data[k+1]);
        hp->t += hp->dt;
    }
    double u_max = gsl_vector_max(u_max_vec);

    gsl_vector_free(u_min_vec);
    gsl_vector_free(u_max_vec);
    gsl_vector_free(u_fem_vals);
    FEMFunction1D_finalize(&u_fem_func);

    return u_max;
}

void set_led_list_rgb_data(
    const size_t led_offset,
    const size_t num_leds,
    const gsl_vector *red_vals,
    const gsl_vector *green_vals,
    const gsl_vector *blue_vals,
    uint8_t *rgb_data_ret
) {
    for (size_t led_id = led_offset; led_id < led_offset + num_leds; led_id++) {
        rgb_data_ret[3*led_id] = 255*gsl_vector_get(red_vals, led_id-led_offset);
        rgb_data_ret[3*led_id + 1] = 255*gsl_vector_get(green_vals, led_id-led_offset);
        rgb_data_ret[3*led_id + 2] = 255*gsl_vector_get(blue_vals, led_id-led_offset);
     }
}

static bool program_interrupted = false;
void interrupt_handler(const int sig) {
    program_interrupted = true;
    printf("Program interrupted, cleaning up...\n");
}

int main(int argc, char *argv[])
{
    // Signal handler for cleaning up if an interrupt signal is received
    signal(SIGINT, interrupt_handler);

    // Get user home directory
    uid_t user_id = getuid();
    const struct passwd *user_info = getpwuid(user_id);
    const char *home_dir = user_info->pw_dir;

    // Allocate memory for the device name
    char serial_device[256];
    params.serial_port.device = serial_device;

    // Read configuration
    char config_path_buffer[sizeof serial_device];
    const char *config_path;
    if (argc < 2) {
        snprintf(
            config_path_buffer, sizeof serial_device, "%s/.config/heat-diffusion-led/config.json", home_dir
        );
        config_path = config_path_buffer;
    } else {
        config_path = argv[1];
    }

    read_config(config_path);

    // Set parameters according to configuration file
    const double x0 = params.domain.x0;
    const double x1 = params.domain.x1;
    const double t0 = params.domain.t0;
    const double t1 = params.domain.t1;

    const double alpha = params.material.alpha;
    const double beta = params.material.beta;
    const double gamma = params.material.gamma;

    const size_t nx = params.computational.nx;
    const size_t nt = params.computational.nt;
    const size_t deg = params.computational.deg;
    const size_t nt_update = params.computational.num_steps_update;
    const size_t num_leds_list_1 = params.computational.num_leds_list_1;
    const size_t num_leds_list_2 = params.computational.num_leds_list_2;
    const size_t num_leds = num_leds_list_1 + num_leds_list_2;
    const size_t num_color_bins = params.colormap.num_bins;

    double *red_fem_func_vals_data = params.colormap.red_fem_func_vals_data;
    double *green_fem_func_vals_data = params.colormap.green_fem_func_vals_data;
    double *blue_fem_func_vals_data = params.colormap.blue_fem_func_vals_data;

    const double dt = (t1 - t0)/nt;

    // Set LED x-positions on list 1
    double xpos = x0;
    double dx = (x1-x0)/(num_leds_list_1-1);
    gsl_vector *led_xpos_list_1 = gsl_vector_alloc(num_leds_list_1);
    for (size_t led_id = 0; led_id < num_leds_list_1; led_id++)
    {
        gsl_vector_set(led_xpos_list_1, led_id, xpos);
        xpos += dx;
    }

    // Set LED x-positions on list 2 
    xpos = x0;
    dx = (x1-x0)/(num_leds_list_2-1);
    gsl_vector *led_xpos_list_2 = gsl_vector_alloc(num_leds_list_2);
    for (size_t led_id = 0; led_id < num_leds_list_2; led_id++)
    {
        gsl_vector_set(led_xpos_list_2, led_id, xpos);
        xpos += dx;
    }

    // Create the mesh for computations
    Mesh1D mesh;
    Mesh1D_init(&mesh, x0, x1, nx, deg);

    // Configure heat problem and solver
    HeatProblem1D hp;
    hp.mesh = &mesh;
    hp.alpha = alpha;
    hp.beta = beta;
    hp.gamma = gamma;
    hp.source = &heat_source;
    hp.t = t0;
    hp.dt = dt;
    klu_defaults(&hp.klu_common);
    HeatProblem1D_init(&hp, &mesh);

    // Perform a calibration run to calculate the temperature range
    const double u_max = calibrate(&hp, nt);
    const double u_min = 0.0;

    // Create a mesh for interpolating colors between u_min and u_max
    Mesh1D range_mesh;
    Mesh1D_init(&range_mesh, u_min, u_max, num_color_bins-1, 1);

    // FEM Function for getting red color value based on temperature
    gsl_vector red_fem_func_vals;
    red_fem_func_vals.size = num_color_bins;
    red_fem_func_vals.stride = 1;
    red_fem_func_vals.data = red_fem_func_vals_data;
    red_fem_func_vals.owner = false;

    FEMFunction1D red_fem_func;
    red_fem_func.vals = &red_fem_func_vals;
    red_fem_func.mesh = &range_mesh;

    // FEM Function for getting green color value based on temperature
    gsl_vector green_fem_func_vals;
    green_fem_func_vals.size = num_color_bins;
    green_fem_func_vals.stride = 1;
    green_fem_func_vals.data = green_fem_func_vals_data;
    green_fem_func_vals.owner = false;

    FEMFunction1D green_fem_func;
    green_fem_func.vals = &green_fem_func_vals;
    green_fem_func.mesh = &range_mesh;

    // FEM Function for getting green color value based on temperature
    gsl_vector blue_fem_func_vals;
    blue_fem_func_vals.size = num_color_bins;
    blue_fem_func_vals.stride = 1;
    blue_fem_func_vals.data = blue_fem_func_vals_data;
    blue_fem_func_vals.owner = false;

    FEMFunction1D blue_fem_func;
    blue_fem_func.vals = &blue_fem_func_vals;
    blue_fem_func.mesh = &range_mesh;

    // FEM function for getting temperature at arbitrary point
    FEMFunction1D u_fem_func;
    FEMFunction1D_init(&u_fem_func, &mesh);

    gsl_vector *u_fem_vals_list_1 = gsl_vector_calloc(num_leds_list_1);
    gsl_vector *u_fem_vals_list_2 = gsl_vector_calloc(num_leds_list_2);

    gsl_vector *red_vals_list_1 = gsl_vector_calloc(num_leds_list_1);
    gsl_vector *green_vals_list_1 = gsl_vector_calloc(num_leds_list_1);
    gsl_vector *blue_vals_list_1 = gsl_vector_calloc(num_leds_list_1);

    gsl_vector *red_vals_list_2 = gsl_vector_calloc(num_leds_list_2);
    gsl_vector *green_vals_list_2 = gsl_vector_calloc(num_leds_list_2);
    gsl_vector *blue_vals_list_2 = gsl_vector_calloc(num_leds_list_2);

    // Data array with rgb values that are sent to microcontroller
    uint8_t rgb_data[3*num_leds];

    // Serial object for communicating with the microcontroller
    SerialPort serial_port;
    SerialPort_open(&params.serial_port);

    for (size_t iter_count = 0; iter_count < params.computational.repeat_num_times; iter_count++) {
        if (program_interrupted)
            break;
        hp.t = t0;
        HeatProblem1D_reset_system(&hp);

        // Assemble linear system matrices, these are constant across time steps
        // since we are using a constant time-step size dt
        HeatProblem1D_setup(&hp);

        for (size_t k = 0; k < nt+1; k++) {
            if (program_interrupted)
                break;

            // Update led colors after performing the number of steps given by nt_update (read from config)
            if (k % nt_update == 0)
            {
                // Get temperature at led positions of first led list
                FEMFunction1D_set_vals_from_vec(&u_fem_func, hp.sol_vec);
                FEMFunction1D_eval(&u_fem_func, led_xpos_list_1, u_fem_vals_list_1);
                FEMFunction1D_eval(&red_fem_func, u_fem_vals_list_1, red_vals_list_1);
                FEMFunction1D_eval(&green_fem_func, u_fem_vals_list_1, green_vals_list_1);
                FEMFunction1D_eval(&blue_fem_func, u_fem_vals_list_1, blue_vals_list_1);

                // Insert rgb values into data array
                const size_t led_offset_list_1 = 0;
                set_led_list_rgb_data(
                    led_offset_list_1,
                    num_leds_list_1,
                    red_vals_list_1,
                    green_vals_list_1,
                    blue_vals_list_1,
                    rgb_data
                );

                // Get temperature at led positions of the second led list
                FEMFunction1D_eval(&u_fem_func, led_xpos_list_2, u_fem_vals_list_2);
                FEMFunction1D_eval(&red_fem_func, u_fem_vals_list_2, red_vals_list_2);
                FEMFunction1D_eval(&green_fem_func, u_fem_vals_list_2, green_vals_list_2);
                FEMFunction1D_eval(&blue_fem_func, u_fem_vals_list_2, blue_vals_list_2);

                // Insert rgb values into data array
                const size_t led_offset_list_2 = num_leds_list_1;
                set_led_list_rgb_data(
                    led_offset_list_2,
                    num_leds_list_2,
                    red_vals_list_2,
                    green_vals_list_2,
                    blue_vals_list_2,
                    rgb_data
                );

                // Send data to microcontroller
                SerialPort_write(&params.serial_port, rgb_data, 3*num_leds);

                // Relax a bit to roughly update in real time
                const size_t sleep_time_usec = nt_update * dt * 1e6;
                usleep(sleep_time_usec);
            }

            // Advance system in time
            HeatProblem1D_step(&hp);

            // RHS changes so we need to reset it
            HeatProblem1D_reset_rhs(&hp);

            // Update time
            hp.t += hp.dt;
        }
    }

    // Finalize
    SerialPort_close(&serial_port);
    free(red_fem_func_vals_data);
    free(green_fem_func_vals_data);
    free(blue_fem_func_vals_data);
    gsl_vector_free(blue_vals_list_2);
    gsl_vector_free(green_vals_list_2);
    gsl_vector_free(red_vals_list_2);
    gsl_vector_free(blue_vals_list_1);
    gsl_vector_free(green_vals_list_1);
    gsl_vector_free(red_vals_list_1);
    gsl_vector_free(u_fem_vals_list_2);
    gsl_vector_free(u_fem_vals_list_1);
    FEMFunction1D_finalize(&u_fem_func);
    Mesh1D_finalize(&range_mesh);
    Mesh1D_finalize(&mesh);
    gsl_vector_free(led_xpos_list_2);
    gsl_vector_free(led_xpos_list_1);
    HeatProblem1D_finalize(&hp);

    return 0;
}
