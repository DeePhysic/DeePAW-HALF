#include "half.h"

#include <ctype.h>
#include <errno.h>
#include <math.h>
#include <stdarg.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef HALF_HAVE_CURL
#include <curl/curl.h>
#endif

#define HALF_ESCN_MAX_REQUEST (2u * 1024u * 1024u)
#define HALF_ESCN_MAX_RESPONSE (128u * 1024u * 1024u)

static void set_error(char *error, int capacity, const char *format, ...) {
  va_list args;
  if (!error || capacity <= 0) return;
  va_start(args, format);
  vsnprintf(error, (size_t)capacity, format, args);
  va_end(args);
  error[capacity - 1] = '\0';
}

static int validate_url(const char *url, char *error, int error_capacity) {
  if (!url || !url[0]) {
    set_error(error, error_capacity, "DeePAW-eSCN base URL is required");
    return HALF_ERROR_INVALID_ARGUMENT;
  }
  if (strncmp(url, "http://", 7) != 0 && strncmp(url, "https://", 8) != 0) {
    set_error(error, error_capacity, "DeePAW-eSCN URL must use http:// or https://");
    return HALF_ERROR_INVALID_ARGUMENT;
  }
  return HALF_SUCCESS;
}

#ifdef HALF_HAVE_CURL
struct buffer {
  char *data;
  size_t size;
  size_t capacity;
};

static int reserve(struct buffer *buffer, size_t extra) {
  size_t needed = buffer->size + extra + 1;
  size_t capacity = buffer->capacity ? buffer->capacity : 4096;
  char *next;
  if (needed > HALF_ESCN_MAX_RESPONSE) return 0;
  while (capacity < needed) {
    if (capacity > HALF_ESCN_MAX_RESPONSE / 2) {
      capacity = HALF_ESCN_MAX_RESPONSE;
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

static int appendf(struct buffer *buffer, const char *format, ...) {
  va_list args, copy;
  int count;
  va_start(args, format);
  va_copy(copy, args);
  count = vsnprintf(NULL, 0, format, copy);
  va_end(copy);
  if (count < 0 || !reserve(buffer, (size_t)count)) {
    va_end(args);
    return 0;
  }
  vsnprintf(buffer->data + buffer->size, buffer->capacity - buffer->size,
            format, args);
  va_end(args);
  buffer->size += (size_t)count;
  return 1;
}

static size_t receive(void *contents, size_t size, size_t count, void *user) {
  struct buffer *buffer = (struct buffer *)user;
  size_t bytes = size * count;
  if (!reserve(buffer, bytes)) return 0;
  memcpy(buffer->data + buffer->size, contents, bytes);
  buffer->size += bytes;
  buffer->data[buffer->size] = '\0';
  return bytes;
}

static char *endpoint_url(const char *base, const char *path) {
  size_t length = strlen(base);
  char *url;
  while (length && base[length - 1] == '/') --length;
  url = (char *)malloc(length + strlen(path) + 1);
  if (!url) return NULL;
  memcpy(url, base, length);
  strcpy(url + length, path);
  return url;
}

static const char *json_value(const char *json, const char *key) {
  char pattern[96];
  const char *found;
  if (snprintf(pattern, sizeof pattern, "\"%s\"", key) >=
      (int)sizeof pattern)
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
  ++value;
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
  if (json && json_string(json, "error", &begin, &length)) {
    int shown = length > 512 ? 512 : (int)length;
    set_error(error, error_capacity, "DeePAW-eSCN HTTP %ld: %.*s", status,
              shown, begin);
  } else {
    set_error(error, error_capacity, "DeePAW-eSCN HTTP request failed with status %ld",
              status);
  }
}

static int base64_value(unsigned char c) {
  if (c >= 'A' && c <= 'Z') return c - 'A';
  if (c >= 'a' && c <= 'z') return c - 'a' + 26;
  if (c >= '0' && c <= '9') return c - '0' + 52;
  if (c == '+') return 62;
  if (c == '/') return 63;
  return -1;
}

static int little_endian_host(void) {
  const uint16_t one = 1;
  return *(const unsigned char *)&one == 1;
}

static int decode_float32(const char *encoded, size_t encoded_length,
                          float *output, int64_t count) {
  unsigned char quartet[4], bytes[3];
  size_t index = 0, q = 0;
  int64_t byte_count = 0;
  if (!output || count < 0) return 0;
  while (index < encoded_length) {
    unsigned char c = (unsigned char)encoded[index++];
    int value;
    if (isspace(c)) continue;
    if (c == '=') value = 64;
    else {
      value = base64_value(c);
      if (value < 0) return 0;
    }
    quartet[q++] = (unsigned char)value;
    if (q == 4) {
      int produced = quartet[2] == 64 ? 1 : (quartet[3] == 64 ? 2 : 3);
      bytes[0] = (unsigned char)((quartet[0] << 2) | (quartet[1] >> 4));
      bytes[1] = (unsigned char)((quartet[1] << 4) | (quartet[2] >> 2));
      bytes[2] = (unsigned char)((quartet[2] << 6) | quartet[3]);
      for (int i = 0; i < produced; ++i) {
        int64_t float_index = byte_count / 4;
        int byte_index = (int)(byte_count % 4);
        unsigned char *target;
        if (float_index >= count) return 0;
        target = (unsigned char *)&output[float_index];
        target[little_endian_host() ? byte_index : 3 - byte_index] = bytes[i];
        ++byte_count;
      }
      q = 0;
    }
  }
  return q == 0 && byte_count == count * 4;
}

static int decode_field(const char *json, const char *key, float *output,
                        int64_t count, char *error, int error_capacity) {
  const char *encoded;
  size_t length;
  if (!json_string(json, key, &encoded, &length)) {
    set_error(error, error_capacity, "DeePAW-eSCN response is missing %s", key);
    return HALF_ERROR_IO;
  }
  if (!decode_float32(encoded, length, output, count)) {
    set_error(error, error_capacity, "invalid %s Base64 float32 payload", key);
    return HALF_ERROR_IO;
  }
  return HALF_SUCCESS;
}

static int parse_grid(const char *json, int32_t grid[3]) {
  const char *value = json_value(json, "grid_shape");
  char *end;
  if (!value || *value++ != '[') return 0;
  for (int i = 0; i < 3; ++i) {
    long item;
    while (*value && isspace((unsigned char)*value)) ++value;
    errno = 0;
    item = strtol(value, &end, 10);
    if (errno || end == value || item < 1 || item > INT32_MAX) return 0;
    grid[i] = (int32_t)item;
    value = end;
    while (*value && isspace((unsigned char)*value)) ++value;
    if (i < 2 && *value++ != ',') return 0;
  }
  return *value == ']';
}

static int perform(CURL *curl, struct buffer *response, long *http_status,
                   char *error, int error_capacity) {
  CURLcode code = curl_easy_perform(curl);
  if (code != CURLE_OK) {
    set_error(error, error_capacity, "DeePAW-eSCN transport error: %s",
              curl_easy_strerror(code));
    return HALF_ERROR_IO;
  }
  if (curl_easy_getinfo(curl, CURLINFO_RESPONSE_CODE, http_status) != CURLE_OK) {
    set_error(error, error_capacity, "cannot read DeePAW-eSCN HTTP status");
    return HALF_ERROR_IO;
  }
  if (!response->data && !reserve(response, 0)) {
    set_error(error, error_capacity, "out of memory receiving DeePAW-eSCN response");
    return HALF_ERROR_INTERNAL;
  }
  response->data[response->size] = '\0';
  return HALF_SUCCESS;
}

int half_escn_health(const char *base_url, int timeout_seconds, char *error,
                     int error_capacity) {
  CURL *curl = NULL;
  struct buffer response = {0};
  char *url = NULL;
  long status = 0;
  int result = validate_url(base_url, error, error_capacity);
  if (result != HALF_SUCCESS) return result;
  if (timeout_seconds <= 0) timeout_seconds = 30;
  url = endpoint_url(base_url, "/v1/health");
  curl = curl_easy_init();
  if (!url || !curl) {
    set_error(error, error_capacity, "cannot initialize DeePAW-eSCN HTTP client");
    result = HALF_ERROR_INTERNAL;
    goto done;
  }
  curl_easy_setopt(curl, CURLOPT_URL, url);
  curl_easy_setopt(curl, CURLOPT_TIMEOUT, (long)timeout_seconds);
  curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, receive);
  curl_easy_setopt(curl, CURLOPT_WRITEDATA, &response);
  curl_easy_setopt(curl, CURLOPT_NOSIGNAL, 1L);
  result = perform(curl, &response, &status, error, error_capacity);
  if (result != HALF_SUCCESS) goto done;
  if (status < 200 || status >= 300) {
    response_error(response.data, status, error, error_capacity);
    result = HALF_ERROR_IO;
  } else if (!strstr(response.data, "\"status\"") ||
             !strstr(response.data, "\"running\"")) {
    set_error(error, error_capacity, "DeePAW-eSCN health response is not running");
    result = HALF_ERROR_UNAVAILABLE;
  }
done:
  if (curl) curl_easy_cleanup(curl);
  free(url);
  free(response.data);
  return result;
}

int half_escn_predict(const char *base_url,
                      const half_escn_request_v1 *request,
                      half_escn_result_v1 *result, int timeout_seconds,
                      char *error, int error_capacity) {
  CURL *curl = NULL;
  struct curl_slist *headers = NULL;
  struct buffer body = {0}, response = {0};
  char *url = NULL;
  long status = 0;
  int rc = validate_url(base_url, error, error_capacity);
  int64_t points;
  if (rc != HALF_SUCCESS) return rc;
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
    return result->capacity < points ? HALF_ERROR_CAPACITY : HALF_ERROR_INVALID_ARGUMENT;
  }
  for (int i = 0; i < request->nions; ++i) {
    if (request->atomic_numbers[i] < 1 || request->atomic_numbers[i] > 118) {
      set_error(error, error_capacity, "atomic numbers must be in [1,118]");
      return HALF_ERROR_INVALID_ARGUMENT;
    }
    for (int j = 0; j < 3; ++j)
      if (!isfinite(request->positions_cartesian[3 * i + j])) {
        set_error(error, error_capacity, "atomic positions must be finite");
        return HALF_ERROR_INVALID_ARGUMENT;
      }
  }
  for (int i = 0; i < 9; ++i)
    if (!isfinite(request->cell[i])) {
      set_error(error, error_capacity, "cell vectors must be finite");
      return HALF_ERROR_INVALID_ARGUMENT;
    }

  if (!appendf(&body, "{\"atoms\":{\"numbers\":[")) goto memory_error;
  for (int i = 0; i < request->nions; ++i)
    if (!appendf(&body, "%s%d", i ? "," : "", request->atomic_numbers[i]))
      goto memory_error;
  if (!appendf(&body, "],\"positions\":[")) goto memory_error;
  for (int i = 0; i < request->nions; ++i)
    if (!appendf(&body, "%s[%.17g,%.17g,%.17g]", i ? "," : "",
                 request->positions_cartesian[3 * i],
                 request->positions_cartesian[3 * i + 1],
                 request->positions_cartesian[3 * i + 2]))
      goto memory_error;
  if (!appendf(&body, "],\"cell\":[")) goto memory_error;
  for (int i = 0; i < 3; ++i)
    if (!appendf(&body, "%s[%.17g,%.17g,%.17g]", i ? "," : "",
                 request->cell[3 * i], request->cell[3 * i + 1],
                 request->cell[3 * i + 2]))
      goto memory_error;
  if (!appendf(&body,
               "],\"pbc\":[true,true,true]},\"grid_shape\":[%d,%d,%d],"
               "\"include_uncertainty\":%s}",
               request->grid[0], request->grid[1], request->grid[2],
               (request->flags & HALF_ESCN_INCLUDE_UNCERTAINTY) ? "true" : "false"))
    goto memory_error;
  if (body.size > HALF_ESCN_MAX_REQUEST) {
    set_error(error, error_capacity, "DeePAW-eSCN JSON request exceeds 2 MiB");
    rc = HALF_ERROR_CAPACITY;
    goto done;
  }

  if (timeout_seconds <= 0) timeout_seconds = 600;
  url = endpoint_url(base_url, "/v1/predict");
  curl = curl_easy_init();
  headers = curl_slist_append(headers, "Content-Type: application/json");
  if (!url || !curl || !headers) goto memory_error;
  curl_easy_setopt(curl, CURLOPT_URL, url);
  curl_easy_setopt(curl, CURLOPT_POST, 1L);
  curl_easy_setopt(curl, CURLOPT_POSTFIELDS, body.data);
  curl_easy_setopt(curl, CURLOPT_POSTFIELDSIZE_LARGE, (curl_off_t)body.size);
  curl_easy_setopt(curl, CURLOPT_HTTPHEADER, headers);
  curl_easy_setopt(curl, CURLOPT_TIMEOUT, (long)timeout_seconds);
  curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, receive);
  curl_easy_setopt(curl, CURLOPT_WRITEDATA, &response);
  curl_easy_setopt(curl, CURLOPT_NOSIGNAL, 1L);
  rc = perform(curl, &response, &status, error, error_capacity);
  if (rc != HALF_SUCCESS) goto done;
  if (status < 200 || status >= 300) {
    response_error(response.data, status, error, error_capacity);
    rc = HALF_ERROR_IO;
    goto done;
  }
  if (!parse_grid(response.data, result->grid) ||
      result->grid[0] != request->grid[0] ||
      result->grid[1] != request->grid[1] ||
      result->grid[2] != request->grid[2]) {
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

memory_error:
  set_error(error, error_capacity, "out of memory building DeePAW-eSCN request");
  rc = HALF_ERROR_INTERNAL;
done:
  if (curl) curl_easy_cleanup(curl);
  if (headers) curl_slist_free_all(headers);
  free(url);
  free(body.data);
  free(response.data);
  return rc;
}

#elif !defined(HALF_HAVE_ESCN_SOCKET)

int half_escn_health(const char *base_url, int timeout_seconds, char *error,
                     int error_capacity) {
  (void)timeout_seconds;
  if (validate_url(base_url, error, error_capacity) != HALF_SUCCESS)
    return HALF_ERROR_INVALID_ARGUMENT;
  set_error(error, error_capacity, "HALF was built without libcurl");
  return HALF_ERROR_UNAVAILABLE;
}

int half_escn_predict(const char *base_url,
                      const half_escn_request_v1 *request,
                      half_escn_result_v1 *result, int timeout_seconds,
                      char *error, int error_capacity) {
  (void)request;
  (void)result;
  (void)timeout_seconds;
  if (validate_url(base_url, error, error_capacity) != HALF_SUCCESS)
    return HALF_ERROR_INVALID_ARGUMENT;
  set_error(error, error_capacity, "HALF was built without libcurl");
  return HALF_ERROR_UNAVAILABLE;
}

#endif
