#include "ns_config.h"

#include <ctype.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define NS_CONFIG_MAX_KEY 256
#define NS_CONFIG_MAX_VALUE 512

static char *ns_config_strdup(const char *src) {
    if (!src) {
        return NULL;
    }
    size_t len = strlen(src) + 1;
    char *dup = (char *)malloc(len);
    if (dup) {
        memcpy(dup, src, len);
    }
    return dup;
}

static void ns_config_trim(char *str) {
    if (!str) return;
    char *start = str;
    while (*start && isspace((unsigned char)*start)) {
        start++;
    }
    char *end = start + strlen(start);
    while (end > start && isspace((unsigned char)end[-1])) {
        end--;
    }
    *end = '\0';
    if (start != str) {
        memmove(str, start, (size_t)(end - start) + 1);
    }
}

static int ns_config_reserve(NSConfig *config, size_t additional) {
    if (!config) return -1;
    size_t required = config->count + additional;
    if (required <= config->capacity) {
        return 0;
    }
    size_t new_capacity = config->capacity ? config->capacity : 16;
    while (new_capacity < required) {
        new_capacity *= 2;
    }
    NSConfigEntry *entries = (NSConfigEntry *)realloc(config->entries, new_capacity * sizeof(NSConfigEntry));
    if (!entries) {
        return -1;
    }
    config->entries = entries;
    config->capacity = new_capacity;
    return 0;
}

void ns_config_init(NSConfig *config) {
    if (!config) return;
    config->entries = NULL;
    config->count = 0;
    config->capacity = 0;
}

void ns_config_free(NSConfig *config) {
    if (!config) return;
    for (size_t i = 0; i < config->count; ++i) {
        free(config->entries[i].key);
        free(config->entries[i].value);
    }
    free(config->entries);
    config->entries = NULL;
    config->count = 0;
    config->capacity = 0;
}

int ns_config_set(NSConfig *config, const char *key, const char *value) {
    if (!config || !key || !value) {
        return -1;
    }
    for (size_t i = 0; i < config->count; ++i) {
        if (strcmp(config->entries[i].key, key) == 0) {
            char *dup = ns_config_strdup(value);
            if (!dup) {
                return -1;
            }
            free(config->entries[i].value);
            config->entries[i].value = dup;
            return 0;
        }
    }
    if (ns_config_reserve(config, 1) != 0) {
        return -1;
    }
    config->entries[config->count].key = ns_config_strdup(key);
    config->entries[config->count].value = ns_config_strdup(value);
    if (!config->entries[config->count].key || !config->entries[config->count].value) {
        free(config->entries[config->count].key);
        free(config->entries[config->count].value);
        return -1;
    }
    config->count += 1;
    return 0;
}

int ns_config_has_key(const NSConfig *config, const char *key) {
    if (!config || !key) return 0;
    for (size_t i = 0; i < config->count; ++i) {
        if (strcmp(config->entries[i].key, key) == 0) {
            return 1;
        }
    }
    return 0;
}

static const char *ns_config_find(const NSConfig *config, const char *key) {
    if (!config || !key) return NULL;
    for (size_t i = 0; i < config->count; ++i) {
        if (strcmp(config->entries[i].key, key) == 0) {
            return config->entries[i].value;
        }
    }
    return NULL;
}

const char *ns_config_get_string(const NSConfig *config, const char *key, const char *default_value) {
    const char *value = ns_config_find(config, key);
    return value ? value : default_value;
}

int ns_config_get_int(const NSConfig *config, const char *key, int default_value) {
    const char *value = ns_config_find(config, key);
    if (!value) {
        return default_value;
    }
    char *endptr = NULL;
    long parsed = strtol(value, &endptr, 10);
    if (!value[0] || (endptr && *endptr != '\0')) {
        return default_value;
    }
    return (int)parsed;
}

double ns_config_get_double(const NSConfig *config, const char *key, double default_value) {
    const char *value = ns_config_find(config, key);
    if (!value) {
        return default_value;
    }
    char *endptr = NULL;
    double parsed = strtod(value, &endptr);
    if (!value[0] || (endptr && *endptr != '\0')) {
        return default_value;
    }
    return parsed;
}

static int ns_config_stricmp(const char *a, const char *b) {
    while (*a && *b) {
        int ca = tolower((unsigned char)*a);
        int cb = tolower((unsigned char)*b);
        if (ca != cb) {
            return ca - cb;
        }
        ++a;
        ++b;
    }
    return tolower((unsigned char)*a) - tolower((unsigned char)*b);
}

int ns_config_get_bool(const NSConfig *config, const char *key, int default_value) {
    const char *value = ns_config_find(config, key);
    if (!value) {
        return default_value;
    }
    if (!value[0]) {
        return default_value;
    }
    if (ns_config_stricmp(value, "true") == 0 ||
        ns_config_stricmp(value, "yes") == 0 ||
        ns_config_stricmp(value, "on") == 0 ||
        strcmp(value, "1") == 0) {
        return 1;
    }
    if (ns_config_stricmp(value, "false") == 0 ||
        ns_config_stricmp(value, "no") == 0 ||
        ns_config_stricmp(value, "off") == 0 ||
        strcmp(value, "0") == 0) {
        return 0;
    }
    return default_value;
}

static void ns_config_set_error(char *error, size_t error_size, const char *message, const char *detail) {
    if (!error || error_size == 0) {
        return;
    }
    if (detail) {
        snprintf(error, error_size, "%s: %s", message, detail);
    } else {
        snprintf(error, error_size, "%s", message);
    }
}

static int ns_config_parse_ini(NSConfig *config, FILE *fp, char *error, size_t error_size) {
    char section[NS_CONFIG_MAX_KEY] = "";
    char line[1024];

    while (fgets(line, sizeof(line), fp)) {
        char *trimmed = line;
        ns_config_trim(trimmed);
        if (trimmed[0] == '\0') {
            continue;
        }
        if (trimmed[0] == '#' || trimmed[0] == ';') {
            continue;
        }
        if (trimmed[0] == '[') {
            char *end = strchr(trimmed, ']');
            if (!end) {
                ns_config_set_error(error, error_size, "Malformed section header in INI file", trimmed);
                return -1;
            }
            *end = '\0';
            strncpy(section, trimmed + 1, sizeof(section) - 1);
            section[sizeof(section) - 1] = '\0';
            ns_config_trim(section);
            continue;
        }
        char *equals = strchr(trimmed, '=');
        if (!equals) {
            ns_config_set_error(error, error_size, "Malformed key-value line in INI file", trimmed);
            return -1;
        }
        *equals = '\0';
        char key_part[NS_CONFIG_MAX_KEY];
        char value_part[NS_CONFIG_MAX_VALUE];
        strncpy(key_part, trimmed, sizeof(key_part) - 1);
        key_part[sizeof(key_part) - 1] = '\0';
        strncpy(value_part, equals + 1, sizeof(value_part) - 1);
        value_part[sizeof(value_part) - 1] = '\0';
        ns_config_trim(key_part);
        ns_config_trim(value_part);

        char full_key[NS_CONFIG_MAX_KEY];
        if (section[0] != '\0') {
            int written = snprintf(full_key, sizeof(full_key), "%s.%s", section, key_part);
            if (written < 0 || (size_t)written >= sizeof(full_key)) {
                ns_config_set_error(error, error_size, "Combined key too long", key_part);
                return -1;
            }
        } else {
            strncpy(full_key, key_part, sizeof(full_key) - 1);
            full_key[sizeof(full_key) - 1] = '\0';
        }

        if (ns_config_set(config, full_key, value_part) != 0) {
            ns_config_set_error(error, error_size, "Failed to store configuration key", full_key);
            return -1;
        }
    }
    return 0;
}

static const char *ns_config_skip_ws(const char *p) {
    while (p && *p && isspace((unsigned char)*p)) {
        ++p;
    }
    return p;
}

static char ns_config_parse_escape(char c) {
    switch (c) {
        case '"': return '"';
        case '\\': return '\\';
        case '/': return '/';
        case 'b': return '\b';
        case 'f': return '\f';
        case 'n': return '\n';
        case 'r': return '\r';
        case 't': return '\t';
        default: return c;
    }
}

static char *ns_config_parse_json_string(const char **cursor, char *error, size_t error_size) {
    const char *p = ns_config_skip_ws(*cursor);
    if (!p || *p != '"') {
        ns_config_set_error(error, error_size, "Expected string in JSON configuration", NULL);
        return NULL;
    }
    ++p;
    size_t capacity = 64;
    size_t length = 0;
    char *buffer = (char *)malloc(capacity);
    if (!buffer) {
        ns_config_set_error(error, error_size, "Out of memory parsing JSON string", NULL);
        return NULL;
    }
    while (*p && *p != '"') {
        char ch = *p;
        if (ch == '\\') {
            ++p;
            if (!*p) {
                free(buffer);
                ns_config_set_error(error, error_size, "Incomplete escape sequence in JSON string", NULL);
                return NULL;
            }
            ch = ns_config_parse_escape(*p);
        }
        if (length + 1 >= capacity) {
            capacity *= 2;
            char *tmp = (char *)realloc(buffer, capacity);
            if (!tmp) {
                free(buffer);
                ns_config_set_error(error, error_size, "Out of memory expanding JSON string", NULL);
                return NULL;
            }
            buffer = tmp;
        }
        buffer[length++] = ch;
        ++p;
    }
    if (*p != '"') {
        free(buffer);
        ns_config_set_error(error, error_size, "Unterminated string in JSON configuration", NULL);
        return NULL;
    }
    buffer[length] = '\0';
    *cursor = p + 1;
    return buffer;
}

static int ns_config_parse_json_value(NSConfig *config,
                                      const char **cursor,
                                      const char *prefix,
                                      char *error,
                                      size_t error_size);

static int ns_config_parse_json_object(NSConfig *config,
                                       const char **cursor,
                                       const char *prefix,
                                       char *error,
                                       size_t error_size) {
    const char *p = ns_config_skip_ws(*cursor);
    if (!p || *p != '{') {
        ns_config_set_error(error, error_size, "Expected '{' in JSON configuration", NULL);
        return -1;
    }
    ++p;
    p = ns_config_skip_ws(p);
    if (*p == '}') {
        *cursor = p + 1;
        return 0;
    }
    while (*p) {
        char *key = ns_config_parse_json_string(&p, error, error_size);
        if (!key) {
            return -1;
        }
        p = ns_config_skip_ws(p);
        if (*p != ':') {
            free(key);
            ns_config_set_error(error, error_size, "Expected ':' in JSON key/value pair", NULL);
            return -1;
        }
        ++p;
        char full_key[NS_CONFIG_MAX_KEY];
        if (prefix && prefix[0]) {
            int written = snprintf(full_key, sizeof(full_key), "%s.%s", prefix, key);
            if (written < 0 || (size_t)written >= sizeof(full_key)) {
                free(key);
                ns_config_set_error(error, error_size, "Combined JSON key too long", NULL);
                return -1;
            }
        } else {
            strncpy(full_key, key, sizeof(full_key) - 1);
            full_key[sizeof(full_key) - 1] = '\0';
        }
        free(key);

        if (ns_config_parse_json_value(config, &p, full_key, error, error_size) != 0) {
            return -1;
        }
        p = ns_config_skip_ws(p);
        if (*p == ',') {
            ++p;
            continue;
        } else if (*p == '}') {
            ++p;
            *cursor = p;
            return 0;
        } else {
            ns_config_set_error(error, error_size, "Expected ',' or '}' in JSON object", NULL);
            return -1;
        }
    }
    ns_config_set_error(error, error_size, "Unexpected end of JSON content", NULL);
    return -1;
}

static int ns_config_parse_json_value(NSConfig *config,
                                      const char **cursor,
                                      const char *prefix,
                                      char *error,
                                      size_t error_size) {
    const char *p = ns_config_skip_ws(*cursor);
    if (!p || !*p) {
        ns_config_set_error(error, error_size, "Unexpected end of JSON value", NULL);
        return -1;
    }
    if (*p == '{') {
        int result = ns_config_parse_json_object(config, &p, prefix, error, error_size);
        if (result != 0) {
            return result;
        }
        *cursor = p;
        return 0;
    }
    if (*p == '"') {
        char *str = ns_config_parse_json_string(&p, error, error_size);
        if (!str) {
            return -1;
        }
        if (ns_config_set(config, prefix, str) != 0) {
            free(str);
            ns_config_set_error(error, error_size, "Failed to store JSON string value", prefix);
            return -1;
        }
        free(str);
        *cursor = p;
        return 0;
    }
    if (strncmp(p, "true", 4) == 0) {
        if (ns_config_set(config, prefix, "true") != 0) {
            ns_config_set_error(error, error_size, "Failed to store JSON boolean", prefix);
            return -1;
        }
        *cursor = p + 4;
        return 0;
    }
    if (strncmp(p, "false", 5) == 0) {
        if (ns_config_set(config, prefix, "false") != 0) {
            ns_config_set_error(error, error_size, "Failed to store JSON boolean", prefix);
            return -1;
        }
        *cursor = p + 5;
        return 0;
    }
    if (strncmp(p, "null", 4) == 0) {
        if (ns_config_set(config, prefix, "") != 0) {
            ns_config_set_error(error, error_size, "Failed to store JSON null", prefix);
            return -1;
        }
        *cursor = p + 4;
        return 0;
    }
    char *endptr = NULL;
    double number = strtod(p, &endptr);
    if (endptr == p) {
        ns_config_set_error(error, error_size, "Invalid JSON value", prefix);
        return -1;
    }
    char buffer[64];
    snprintf(buffer, sizeof(buffer), "%.15g", number);
    if (ns_config_set(config, prefix, buffer) != 0) {
        ns_config_set_error(error, error_size, "Failed to store JSON numeric value", prefix);
        return -1;
    }
    *cursor = endptr;
    return 0;
}

static int ns_config_parse_json(NSConfig *config, const char *content, char *error, size_t error_size) {
    const char *cursor = ns_config_skip_ws(content);
    int result = ns_config_parse_json_object(config, &cursor, "", error, error_size);
    if (result != 0) {
        return result;
    }
    cursor = ns_config_skip_ws(cursor);
    if (*cursor != '\0') {
        ns_config_set_error(error, error_size, "Unexpected trailing data in JSON file", NULL);
        return -1;
    }
    return 0;
}

static int ns_config_load_file(NSConfig *config, const char *path, char *error, size_t error_size) {
    const char *extension = strrchr(path, '.');
    int is_json = (extension && ns_config_stricmp(extension, ".json") == 0);

    if (is_json) {
        FILE *fp = fopen(path, "rb");
        if (!fp) {
            ns_config_set_error(error, error_size, "Failed to open JSON configuration file", strerror(errno));
            return -1;
        }
        if (fseek(fp, 0, SEEK_END) != 0) {
            fclose(fp);
            ns_config_set_error(error, error_size, "Failed to seek JSON configuration file", path);
            return -1;
        }
        long length = ftell(fp);
        if (length < 0) {
            fclose(fp);
            ns_config_set_error(error, error_size, "Failed to determine JSON file size", path);
            return -1;
        }
        if (fseek(fp, 0, SEEK_SET) != 0) {
            fclose(fp);
            ns_config_set_error(error, error_size, "Failed to rewind JSON configuration file", path);
            return -1;
        }
        char *buffer = (char *)malloc((size_t)length + 1);
        if (!buffer) {
            fclose(fp);
            ns_config_set_error(error, error_size, "Out of memory reading JSON configuration file", path);
            return -1;
        }
        size_t read = fread(buffer, 1, (size_t)length, fp);
        fclose(fp);
        if (read != (size_t)length) {
            free(buffer);
            ns_config_set_error(error, error_size, "Failed to read JSON configuration file", path);
            return -1;
        }
        buffer[length] = '\0';
        int result = ns_config_parse_json(config, buffer, error, error_size);
        free(buffer);
        return result;
    }

    FILE *fp = fopen(path, "r");
    if (!fp) {
        ns_config_set_error(error, error_size, "Failed to open INI configuration file", strerror(errno));
        return -1;
    }
    int result = ns_config_parse_ini(config, fp, error, error_size);
    fclose(fp);
    return result;
}

int ns_config_load_files(NSConfig *config, const char **paths, size_t count, char *error, size_t error_size) {
    if (!config || (!paths && count > 0)) {
        ns_config_set_error(error, error_size, "Invalid arguments to ns_config_load_files", NULL);
        return -1;
    }
    for (size_t i = 0; i < count; ++i) {
        if (ns_config_load_file(config, paths[i], error, error_size) != 0) {
            return -1;
        }
    }
    return 0;
}

int ns_config_validate(const NSConfig *config,
                       const NSConfigSchemaEntry *schema,
                       size_t schema_count,
                       char *error,
                       size_t error_size) {
    if (!config || !schema) {
        ns_config_set_error(error, error_size, "Invalid arguments to ns_config_validate", NULL);
        return -1;
    }
    for (size_t i = 0; i < schema_count; ++i) {
        const NSConfigSchemaEntry *entry = &schema[i];
        const char *value = ns_config_find(config, entry->key);
        if (!value) {
            if (entry->required) {
                ns_config_set_error(error, error_size, "Missing required configuration key", entry->key);
                return -1;
            }
            continue;
        }
        if (entry->type == NS_CONFIG_VALUE_STRING) {
            if (entry->allowed_values && entry->allowed_count > 0) {
                int match = 0;
                for (size_t v = 0; v < entry->allowed_count; ++v) {
                    if (ns_config_stricmp(value, entry->allowed_values[v]) == 0) {
                        match = 1;
                        break;
                    }
                }
                if (!match) {
                    ns_config_set_error(error, error_size, "Configuration string value not permitted", entry->key);
                    return -1;
                }
            }
        } else if (entry->type == NS_CONFIG_VALUE_INT) {
            char *endptr = NULL;
            long parsed = strtol(value, &endptr, 10);
            if (!value[0] || (endptr && *endptr != '\0')) {
                ns_config_set_error(error, error_size, "Configuration key expects integer value", entry->key);
                return -1;
            }
            if (entry->has_min && parsed < entry->min_value) {
                ns_config_set_error(error, error_size, "Configuration integer below minimum", entry->key);
                return -1;
            }
            if (entry->has_max && parsed > entry->max_value) {
                ns_config_set_error(error, error_size, "Configuration integer above maximum", entry->key);
                return -1;
            }
        } else if (entry->type == NS_CONFIG_VALUE_DOUBLE) {
            char *endptr = NULL;
            double parsed = strtod(value, &endptr);
            if (!value[0] || (endptr && *endptr != '\0')) {
                ns_config_set_error(error, error_size, "Configuration key expects numeric value", entry->key);
                return -1;
            }
            if (entry->has_min && parsed < entry->min_value) {
                ns_config_set_error(error, error_size, "Configuration value below minimum", entry->key);
                return -1;
            }
            if (entry->has_max && parsed > entry->max_value) {
                ns_config_set_error(error, error_size, "Configuration value above maximum", entry->key);
                return -1;
            }
        } else if (entry->type == NS_CONFIG_VALUE_BOOL) {
            if (!(ns_config_stricmp(value, "true") == 0 ||
                  ns_config_stricmp(value, "false") == 0 ||
                  strcmp(value, "0") == 0 ||
                  strcmp(value, "1") == 0 ||
                  ns_config_stricmp(value, "yes") == 0 ||
                  ns_config_stricmp(value, "no") == 0 ||
                  ns_config_stricmp(value, "on") == 0 ||
                  ns_config_stricmp(value, "off") == 0)) {
                ns_config_set_error(error, error_size, "Configuration key expects boolean value", entry->key);
                return -1;
            }
        }
    }
    return 0;
}
