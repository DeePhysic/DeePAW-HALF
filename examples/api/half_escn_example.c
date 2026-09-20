#include "half.h"

#include <float.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv) {
  const char *url = argc > 1 ? argv[1] : "http://127.0.0.1:8265";
  const int32_t numbers[2] = {14, 14};
  const double positions[6] = {0.0, 0.0, 0.0, 1.3575, 1.3575, 1.3575};
  half_escn_request_v1 request = {0};
  half_escn_result_v1 result = {0};
  char error[512];
  float *density;
  double mean = 0.0;
  float minimum = FLT_MAX, maximum = -FLT_MAX;
  const int points = 16 * 16 * 16;
  if (!(half_get_capabilities() & HALF_CAP_ESCN_API)) {
    fprintf(stderr, "HALF was built without DeePAW-eSCN API support\n");
    return 2;
  }
  if (half_escn_health(url, 30, error, sizeof error) != HALF_SUCCESS) {
    fprintf(stderr, "health check failed: %s\n", error);
    return 3;
  }
  density = (float *)malloc((size_t)points * sizeof *density);
  if (!density) return 4;
  request.struct_size = sizeof request;
  request.nions = 2;
  request.grid[0] = request.grid[1] = request.grid[2] = 16;
  request.cell[0] = request.cell[4] = request.cell[8] = 5.43;
  request.atomic_numbers = numbers;
  request.positions_cartesian = positions;
  result.struct_size = sizeof result;
  result.capacity = points;
  result.density = density;
  if (half_escn_predict(url, &request, &result, 600, error, sizeof error) !=
      HALF_SUCCESS) {
    fprintf(stderr, "prediction failed: %s\n", error);
    free(density);
    return 5;
  }
  for (int i = 0; i < points; ++i) {
    if (density[i] < minimum) minimum = density[i];
    if (density[i] > maximum) maximum = density[i];
    mean += density[i];
  }
  mean /= points;
  printf("grid: %d %d %d\n", result.grid[0], result.grid[1], result.grid[2]);
  printf("density range: %.9g %.9g\n", minimum, maximum);
  printf("density mean: %.12g\n", mean);
  printf("server elapsed: %.6f seconds\n", result.elapsed_seconds);
  free(density);
  return 0;
}
