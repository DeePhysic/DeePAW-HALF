#define _POSIX_C_SOURCE 200809L
#include "half.h"

#include <ctype.h>
#include <errno.h>
#include <math.h>
#include <netdb.h>
#include <stdarg.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/socket.h>
#include <sys/time.h>
#include <sys/types.h>
#include <unistd.h>

#define MAX_REQUEST (2u * 1024u * 1024u)
#define MAX_RESPONSE (128u * 1024u * 1024u)

struct bytes {
  char *data;
  size_t size;
  size_t capacity;
};

static void set_error(char *error, int capacity, const char *format, ...) {
  va_list args;
  if (!error || capacity <= 0) return;
  va_start(args, format);
  vsnprintf(error, (size_t)capacity, format, args);
  va_end(args);
  error[capacity - 1] = '\0';
}

static int reserve(struct bytes *buffer, size_t extra, size_t maximum) {
  size_t needed = buffer->size + extra + 1;
  size_t capacity = buffer->capacity ? buffer->capacity : 4096;
  char *next;
  if (needed > maximum) return 0;
  while (capacity < needed) {
    if (capacity > maximum / 2) {
      capacity = maximum;
      break;
    }
    capacity *= 2;
  }
  if (capacity < needed) return 0;
  next = (char *)realloc(buffer->data, capacity);
  if (!next) return 0;
  buffer->data = next;
  buffer->capacity = capacity;
  return 1;
}

static int appendf(struct bytes *buffer, size_t maximum, const char *format,
                   ...) {
  va_list args, copy;
  int count;
  va_start(args, format);
  va_copy(copy, args);
  count = vsnprintf(NULL, 0, format, copy);
  va_end(copy);
  if (count < 0 || !reserve(buffer, (size_t)count, maximum)) {
    va_end(args);
    return 0;
  }
  vsnprintf(buffer->data + buffer->size, buffer->capacity - buffer->size,
            format, args);
  va_end(args);
  buffer->size += (size_t)count;
  return 1;
}

static char *endpoint_url(const char *base, const char *path) {
  size_t length;
  char *url;
  if (!base || strncmp(base, "http://", 7) != 0) return NULL;
  length = strlen(base);
  while (length && base[length - 1] == '/') --length;
  url = (char *)malloc(length + strlen(path) + 1);
  if (!url) return NULL;
  memcpy(url, base, length);
  strcpy(url + length, path);
  return url;
}

static int parse_url(const char *url, char **host, char **port, char **path,
                     char *error, int error_capacity) {
  const char *authority, *slash, *colon;
  size_t host_length;
  if (!url || strncmp(url, "http://", 7) != 0) {
    set_error(error, error_capacity,
              "built-in eSCN transport supports only http:// URLs");
    return HALF_ERROR_INVALID_ARGUMENT;
  }
  authority = url + 7;
  slash = strchr(authority, '/');
  if (!slash) slash = url + strlen(url);
  colon = memchr(authority, ':', (size_t)(slash - authority));
  host_length = (size_t)((colon ? colon : slash) - authority);
  if (!host_length) {
    set_error(error, error_capacity, "eSCN URL host is empty");
    return HALF_ERROR_INVALID_ARGUMENT;
  }
  *host = (char *)malloc(host_length + 1);
  *port = colon ? strndup(colon + 1, (size_t)(slash - colon - 1)) : strdup("80");
  *path = *slash ? strdup(slash) : strdup("/");
  if (!*host || !*port || !*path) {
    set_error(error, error_capacity, "out of memory parsing eSCN URL");
    return HALF_ERROR_INTERNAL;
  }
  memcpy(*host, authority, host_length);
  (*host)[host_length] = '\0';
  return HALF_SUCCESS;
}

static int send_all(int socket_fd, const char *data, size_t size) {
  while (size) {
    ssize_t sent = send(socket_fd, data, size, 0);
    if (sent < 0 && errno == EINTR) continue;
    if (sent <= 0) return 0;
    data += sent;
    size -= (size_t)sent;
  }
  return 1;
}

static int http_request(const char *method, const char *url, const char *body,
                        size_t body_size, int timeout_seconds,
                        struct bytes *response, long *status, char *error,
                        int error_capacity) {
  char *host = NULL, *port = NULL, *path = NULL;
  struct addrinfo hints = {0}, *addresses = NULL, *address;
  struct bytes request = {0};
  struct timeval timeout;
  int socket_fd = -1, rc;
  char chunk[16384];
  char *header_end;
  ssize_t received;
  rc = parse_url(url, &host, &port, &path, error, error_capacity);
  if (rc != HALF_SUCCESS) goto done;
  hints.ai_family = AF_UNSPEC;
  hints.ai_socktype = SOCK_STREAM;
  rc = getaddrinfo(host, port, &hints, &addresses);
  if (rc != 0) {
    set_error(error, error_capacity, "cannot resolve eSCN host %s: %s", host,
              gai_strerror(rc));
    rc = HALF_ERROR_IO;
    goto done;
  }
  timeout.tv_sec = timeout_seconds > 0 ? timeout_seconds : 600;
  timeout.tv_usec = 0;
  for (address = addresses; address; address = address->ai_next) {
    socket_fd = socket(address->ai_family, address->ai_socktype,
                       address->ai_protocol);
    if (socket_fd < 0) continue;
    setsockopt(socket_fd, SOL_SOCKET, SO_RCVTIMEO, &timeout, sizeof timeout);
    setsockopt(socket_fd, SOL_SOCKET, SO_SNDTIMEO, &timeout, sizeof timeout);
    if (connect(socket_fd, address->ai_addr, address->ai_addrlen) == 0) break;
    close(socket_fd);
    socket_fd = -1;
  }
  if (socket_fd < 0) {
    set_error(error, error_capacity, "cannot connect to eSCN endpoint %s", url);
    rc = HALF_ERROR_IO;
    goto done;
  }
  if (!appendf(&request, MAX_REQUEST + 8192,
               "%s %s HTTP/1.1\r\nHost: %s:%s\r\nAccept: application/json\r\n"
               "Connection: close\r\n",
               method, path, host, port) ||
      (body && !appendf(&request, MAX_REQUEST + 8192,
                        "Content-Type: application/json\r\nContent-Length: %zu\r\n",
                        body_size)) ||
      !appendf(&request, MAX_REQUEST + 8192, "\r\n")) {
    set_error(error, error_capacity, "out of memory building eSCN HTTP request");
    rc = HALF_ERROR_INTERNAL;
    goto done;
  }
  if (body) {
    if (!reserve(&request, body_size, MAX_REQUEST + 8192)) {
      set_error(error, error_capacity, "out of memory building eSCN HTTP request");
      rc = HALF_ERROR_INTERNAL;
      goto done;
    }
    memcpy(request.data + request.size, body, body_size);
    request.size += body_size;
    request.data[request.size] = '\0';
  }
  if (!send_all(socket_fd, request.data, request.size)) {
    set_error(error, error_capacity, "failed to send eSCN HTTP request");
    rc = HALF_ERROR_IO;
    goto done;
  }
  while ((received = recv(socket_fd, chunk, sizeof chunk, 0)) > 0) {
    if (!reserve(response, (size_t)received, MAX_RESPONSE)) {
      set_error(error, error_capacity, "eSCN response exceeds 128 MiB");
      rc = HALF_ERROR_CAPACITY;
      goto done;
    }
    memcpy(response->data + response->size, chunk, (size_t)received);
    response->size += (size_t)received;
  }
  if (received < 0) {
    set_error(error, error_capacity, "timed out receiving eSCN HTTP response");
    rc = HALF_ERROR_IO;
    goto done;
  }
  if (!reserve(response, 0, MAX_RESPONSE)) {
    rc = HALF_ERROR_INTERNAL;
    goto done;
  }
  response->data[response->size] = '\0';
  if (sscanf(response->data, "HTTP/%*u.%*u %ld", status) != 1 ||
      !(header_end = strstr(response->data, "\r\n\r\n"))) {
    set_error(error, error_capacity, "invalid eSCN HTTP response");
    rc = HALF_ERROR_IO;
    goto done;
  }
  header_end += 4;
  response->size -= (size_t)(header_end - response->data);
  memmove(response->data, header_end, response->size);
  response->data[response->size] = '\0';
  rc = HALF_SUCCESS;
done:
  if (socket_fd >= 0) close(socket_fd);
  if (addresses) freeaddrinfo(addresses);
  free(host);
  free(port);
  free(path);
  free(request.data);
  return rc;
}

static const char *json_value(const char *json, const char *key) {
  char pattern[96];
  const char *found;
  if (snprintf(pattern, sizeof pattern, "\"%s\"", key) >= (int)sizeof pattern)
    return NULL;
  found = strstr(json, pattern);
  if (!found) return NULL;
  found += strlen(pattern);
  while (*found && isspace((unsigned char)*found)) ++found;
  if (*found++ != ':') return NULL;
  while (*found && isspace((unsigned char)*found)) ++found;
  return found;
}

static int json_string(const char *json, const char *key, const char **begin,
                       size_t *length) {
  const char *value = json_value(json, key), *end;
  if (!value || *value != '"') return 0;
  value++;
  end = strchr(value, '"');
  if (!end) return 0;
  *begin = value;
  *length = (size_t)(end - value);
  return 1;
}

static void response_error(const char *json, long status, char *error,
                           int error_capacity) {
  const char *begin;
  size_t length;
  if (json && json_string(json, "error", &begin, &length))
    set_error(error, error_capacity, "DeePAW-eSCN HTTP %ld: %.*s", status,
              (int)(length > 512 ? 512 : length), begin);
  else
    set_error(error, error_capacity, "DeePAW-eSCN HTTP status %ld", status);
}

static int b64(unsigned char value) {
  if (value >= 'A' && value <= 'Z') return value - 'A';
  if (value >= 'a' && value <= 'z') return value - 'a' + 26;
  if (value >= '0' && value <= '9') return value - '0' + 52;
  if (value == '+') return 62;
  if (value == '/') return 63;
  return -1;
}

static int decode_float32(const char *encoded, size_t length, float *output,
                          int64_t count) {
  unsigned char quartet[4], decoded[3];
  const uint16_t one = 1;
  int little = *(const unsigned char *)&one == 1;
  size_t q = 0;
  int64_t bytes = 0;
  for (size_t index = 0; index < length; ++index) {
    int value;
    unsigned char character = (unsigned char)encoded[index];
    if (isspace(character)) continue;
    value = character == '=' ? 64 : b64(character);
    if (value < 0) return 0;
    quartet[q++] = (unsigned char)value;
    if (q == 4) {
      int produced = quartet[2] == 64 ? 1 : quartet[3] == 64 ? 2 : 3;
      decoded[0] = (unsigned char)((quartet[0] << 2) | (quartet[1] >> 4));
      decoded[1] = (unsigned char)((quartet[1] << 4) | (quartet[2] >> 2));
      decoded[2] = (unsigned char)((quartet[2] << 6) | quartet[3]);
      for (int i = 0; i < produced; ++i) {
        int64_t word = bytes / 4;
        int byte = (int)(bytes % 4);
        unsigned char *target;
        if (word >= count) return 0;
        target = (unsigned char *)&output[word];
        target[little ? byte : 3 - byte] = decoded[i];
        bytes++;
      }
      q = 0;
    }
  }
  return q == 0 && bytes == count * 4;
}

static int decode_field(const char *json, const char *key, float *output,
                        int64_t count, char *error, int error_capacity) {
  const char *encoded;
  size_t length;
  if (!json_string(json, key, &encoded, &length) ||
      !decode_float32(encoded, length, output, count)) {
    set_error(error, error_capacity, "invalid or missing eSCN %s payload", key);
    return HALF_ERROR_IO;
  }
  return HALF_SUCCESS;
}

static int parse_grid(const char *json, int32_t grid[3]) {
  const char *value = json_value(json, "grid_shape");
  char *end;
  if (!value || *value++ != '[') return 0;
  for (int i = 0; i < 3; ++i) {
    long item = strtol(value, &end, 10);
    if (end == value || item < 1 || item > INT32_MAX) return 0;
    grid[i] = (int32_t)item;
    value = end;
    while (isspace((unsigned char)*value)) value++;
    if (i < 2 && *value++ != ',') return 0;
  }
  return *value == ']';
}

int half_escn_health(const char *base_url, int timeout_seconds, char *error,
                     int error_capacity) {
  struct bytes response = {0};
  char *url = endpoint_url(base_url, "/v1/health");
  long status = 0;
  int rc;
  if (!url) {
    set_error(error, error_capacity,
              "built-in eSCN transport requires an http:// base URL");
    return HALF_ERROR_INVALID_ARGUMENT;
  }
  rc = http_request("GET", url, NULL, 0,
                    timeout_seconds > 0 ? timeout_seconds : 30, &response,
                    &status, error, error_capacity);
  if (rc == HALF_SUCCESS && (status < 200 || status >= 300)) {
    response_error(response.data, status, error, error_capacity);
    rc = HALF_ERROR_IO;
  } else if (rc == HALF_SUCCESS &&
             (!strstr(response.data, "\"status\"") ||
              !strstr(response.data, "\"running\""))) {
    set_error(error, error_capacity, "DeePAW-eSCN health response is not running");
    rc = HALF_ERROR_UNAVAILABLE;
  }
  free(url);
  free(response.data);
  return rc;
}

int half_escn_predict(const char *base_url,
                      const half_escn_request_v1 *request,
                      half_escn_result_v1 *result, int timeout_seconds,
                      char *error, int error_capacity) {
  struct bytes body = {0}, response = {0};
  char *url = NULL;
  long status = 0;
  int64_t points;
  int rc = HALF_ERROR_INTERNAL;
  if (!request || request->struct_size < sizeof(*request) ||
      (request->flags & ~HALF_ESCN_INCLUDE_UNCERTAINTY) || request->nions < 1 ||
      request->nions > 10000 ||
      !request->atomic_numbers || !request->positions_cartesian || !result ||
      result->struct_size < sizeof(*result) || !result->density) {
    set_error(error, error_capacity, "invalid DeePAW-eSCN request or result structure");
    return HALF_ERROR_INVALID_ARGUMENT;
  }
  if (request->grid[0] < 1 || request->grid[1] < 1 || request->grid[2] < 1 ||
      request->grid[0] > 2000000 || request->grid[1] > 2000000 ||
      request->grid[2] > 2000000) {
    set_error(error, error_capacity, "invalid DeePAW-eSCN grid");
    return HALF_ERROR_INVALID_ARGUMENT;
  }
  points = (int64_t)request->grid[0] * request->grid[1] * request->grid[2];
  if (points > 2000000 || result->capacity < points) {
    set_error(error, error_capacity, "invalid DeePAW-eSCN grid or output capacity");
    return result->capacity < points ? HALF_ERROR_CAPACITY
                                     : HALF_ERROR_INVALID_ARGUMENT;
  }
  if (!appendf(&body, MAX_REQUEST, "{\"atoms\":{\"numbers\":[")) goto memory;
  for (int i = 0; i < request->nions; ++i) {
    if (request->atomic_numbers[i] < 1 || request->atomic_numbers[i] > 118) {
      set_error(error, error_capacity, "atomic numbers must be in [1,118]");
      rc = HALF_ERROR_INVALID_ARGUMENT;
      goto done;
    }
    if (!appendf(&body, MAX_REQUEST, "%s%d", i ? "," : "",
                 request->atomic_numbers[i]))
      goto memory;
  }
  if (!appendf(&body, MAX_REQUEST, "],\"positions\":[")) goto memory;
  for (int i = 0; i < request->nions; ++i) {
    const double *position = request->positions_cartesian + 3 * i;
    if (!isfinite(position[0]) || !isfinite(position[1]) ||
        !isfinite(position[2])) {
      set_error(error, error_capacity, "atomic positions must be finite");
      rc = HALF_ERROR_INVALID_ARGUMENT;
      goto done;
    }
    if (!appendf(&body, MAX_REQUEST, "%s[%.17g,%.17g,%.17g]", i ? "," : "",
                 position[0], position[1], position[2]))
      goto memory;
  }
  if (!appendf(&body, MAX_REQUEST, "],\"cell\":[")) goto memory;
  for (int i = 0; i < 3; ++i) {
    const double *cell = request->cell + 3 * i;
    if (!isfinite(cell[0]) || !isfinite(cell[1]) || !isfinite(cell[2])) {
      set_error(error, error_capacity, "cell vectors must be finite");
      rc = HALF_ERROR_INVALID_ARGUMENT;
      goto done;
    }
    if (!appendf(&body, MAX_REQUEST, "%s[%.17g,%.17g,%.17g]", i ? "," : "",
                 cell[0], cell[1], cell[2]))
      goto memory;
  }
  if (!appendf(&body, MAX_REQUEST,
               "],\"pbc\":[true,true,true]},\"grid_shape\":[%d,%d,%d],"
               "\"include_uncertainty\":%s}",
               request->grid[0], request->grid[1], request->grid[2],
               request->flags & HALF_ESCN_INCLUDE_UNCERTAINTY ? "true" : "false"))
    goto memory;
  url = endpoint_url(base_url, "/v1/predict");
  if (!url) {
    set_error(error, error_capacity,
              "built-in eSCN transport requires an http:// base URL");
    rc = HALF_ERROR_INVALID_ARGUMENT;
    goto done;
  }
  rc = http_request("POST", url, body.data, body.size,
                    timeout_seconds > 0 ? timeout_seconds : 600, &response,
                    &status, error, error_capacity);
  if (rc != HALF_SUCCESS) goto done;
  if (status < 200 || status >= 300) {
    response_error(response.data, status, error, error_capacity);
    rc = HALF_ERROR_IO;
    goto done;
  }
  if (!parse_grid(response.data, result->grid) ||
      memcmp(result->grid, request->grid, sizeof result->grid) != 0) {
    set_error(error, error_capacity, "DeePAW-eSCN returned an unexpected grid_shape");
    rc = HALF_ERROR_IO;
    goto done;
  }
  rc = decode_field(response.data, "density_b64", result->density, points,
                    error, error_capacity);
  if (rc != HALF_SUCCESS) goto done;
  result->flags = request->flags;
  result->elapsed_seconds = 0.0;
  {
    const char *elapsed = json_value(response.data, "elapsed");
    if (elapsed) result->elapsed_seconds = strtod(elapsed, NULL);
  }
  if (request->flags & HALF_ESCN_INCLUDE_UNCERTAINTY) {
    if (!result->nu || !result->alpha || !result->beta || !result->risk) {
      set_error(error, error_capacity, "uncertainty output buffers are required");
      rc = HALF_ERROR_INVALID_ARGUMENT;
      goto done;
    }
    if ((rc = decode_field(response.data, "nu", result->nu, points, error,
                           error_capacity)) != HALF_SUCCESS ||
        (rc = decode_field(response.data, "alpha", result->alpha, points, error,
                           error_capacity)) != HALF_SUCCESS ||
        (rc = decode_field(response.data, "beta", result->beta, points, error,
                           error_capacity)) != HALF_SUCCESS ||
        (rc = decode_field(response.data, "risk", result->risk, points, error,
                           error_capacity)) != HALF_SUCCESS)
      goto done;
  }
  if (error && error_capacity > 0) error[0] = '\0';
  rc = HALF_SUCCESS;
  goto done;
memory:
  set_error(error, error_capacity, "eSCN JSON request exceeds 2 MiB");
  rc = HALF_ERROR_CAPACITY;
done:
  free(url);
  free(body.data);
  free(response.data);
  return rc;
}
