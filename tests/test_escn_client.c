#include "half.h"

#include <arpa/inet.h>
#include <math.h>
#include <netinet/in.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/socket.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>

static void serve(int listener) {
  static const char health[] =
      "HTTP/1.1 200 OK\r\nContent-Type: application/json\r\n"
      "Content-Length: 20\r\nConnection: close\r\n\r\n"
      "{\"status\":\"running\"}";
  static const char predict_body[] =
      "{\"density_b64\":\"AACAPwAAAEAAAEBAAACAQA==\","
      "\"grid_shape\":[1,2,2],\"elapsed\":0.125}";
  char request[16384], response[32768];
  for (int call = 0; call < 2; ++call) {
    int client = accept(listener, NULL, NULL);
    ssize_t size;
    if (client < 0) _exit(10);
    size = read(client, request, sizeof request - 1);
    if (size <= 0) _exit(11);
    request[size] = '\0';
    if (strncmp(request, "GET /v1/health ", 15) == 0) {
      write(client, health, sizeof health - 1);
    } else if (strncmp(request, "POST /v1/predict ", 17) == 0 &&
               strstr(request, "\"numbers\":[14]") &&
               strstr(request, "\"grid_shape\":[1,2,2]")) {
      int length = snprintf(response, sizeof response,
                            "HTTP/1.1 200 OK\r\nContent-Type: application/json\r\n"
                            "Content-Length: %zu\r\nConnection: close\r\n\r\n%s",
                            strlen(predict_body), predict_body);
      write(client, response, (size_t)length);
    } else {
      static const char bad[] =
          "HTTP/1.1 400 Bad Request\r\nContent-Length: 23\r\nConnection: close\r\n\r\n"
          "{\"error\":\"bad request\"}";
      write(client, bad, sizeof bad - 1);
    }
    close(client);
  }
  close(listener);
  _exit(0);
}

int main(void) {
  int listener = socket(AF_INET, SOCK_STREAM, 0), status;
  struct sockaddr_in address = {0};
  socklen_t address_size = sizeof address;
  pid_t child;
  char url[128], error[256];
  const int32_t numbers[1] = {14};
  const double positions[3] = {0.0, 0.0, 0.0};
  float density[4] = {0};
  half_escn_request_v1 request = {0};
  half_escn_result_v1 result = {0};
  if (listener < 0) return 1;
  address.sin_family = AF_INET;
  address.sin_addr.s_addr = htonl(INADDR_LOOPBACK);
  address.sin_port = 0;
  if (bind(listener, (struct sockaddr *)&address, sizeof address) ||
      listen(listener, 2) ||
      getsockname(listener, (struct sockaddr *)&address, &address_size))
    return 2;
  child = fork();
  if (child < 0) return 3;
  if (child == 0) serve(listener);
  snprintf(url, sizeof url, "http://127.0.0.1:%u", ntohs(address.sin_port));
  if (half_escn_health(url, 5, error, sizeof error) != HALF_SUCCESS) {
    fprintf(stderr, "health: %s\n", error);
    return 4;
  }
  request.struct_size = sizeof request;
  request.nions = 1;
  request.grid[0] = 1;
  request.grid[1] = 2;
  request.grid[2] = 2;
  request.cell[0] = request.cell[4] = request.cell[8] = 5.0;
  request.atomic_numbers = numbers;
  request.positions_cartesian = positions;
  result.struct_size = sizeof result;
  result.capacity = 4;
  result.density = density;
  if (half_escn_predict(url, &request, &result, 5, error, sizeof error) !=
      HALF_SUCCESS) {
    fprintf(stderr, "predict: %s\n", error);
    return 5;
  }
  if (result.grid[0] != 1 || result.grid[1] != 2 || result.grid[2] != 2 ||
      fabs(result.elapsed_seconds - 0.125) > 1e-12)
    return 6;
  for (int i = 0; i < 4; ++i)
    if (fabsf(density[i] - (float)(i + 1)) > 1e-6f) return 7;
  if (waitpid(child, &status, 0) != child || !WIFEXITED(status) ||
      WEXITSTATUS(status) != 0)
    return 8;
  close(listener);
  puts("HALF DeePAW-eSCN client: ok");
  return 0;
}
