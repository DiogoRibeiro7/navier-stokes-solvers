#ifndef NS_CONFIG_H
#define NS_CONFIG_H

#include <stddef.h>

typedef struct {
    char *key;
    char *value;
} NSConfigEntry;

typedef struct {
    NSConfigEntry *entries;
    size_t count;
    size_t capacity;
} NSConfig;

typedef enum {
    NS_CONFIG_VALUE_STRING,
    NS_CONFIG_VALUE_INT,
    NS_CONFIG_VALUE_DOUBLE,
    NS_CONFIG_VALUE_BOOL
} NSConfigValueType;

typedef struct {
    const char *key;
    NSConfigValueType type;
    int required;
    double min_value;
    double max_value;
    int has_min;
    int has_max;
    const char **allowed_values;
    size_t allowed_count;
} NSConfigSchemaEntry;

void ns_config_init(NSConfig *config);
void ns_config_free(NSConfig *config);

int ns_config_set(NSConfig *config, const char *key, const char *value);
const char *ns_config_get_string(const NSConfig *config, const char *key, const char *default_value);
int ns_config_get_int(const NSConfig *config, const char *key, int default_value);
double ns_config_get_double(const NSConfig *config, const char *key, double default_value);
int ns_config_get_bool(const NSConfig *config, const char *key, int default_value);
int ns_config_has_key(const NSConfig *config, const char *key);

int ns_config_load_files(NSConfig *config, const char **paths, size_t count, char *error, size_t error_size);

int ns_config_validate(const NSConfig *config,
                       const NSConfigSchemaEntry *schema,
                       size_t schema_count,
                       char *error,
                       size_t error_size);

#endif // NS_CONFIG_H
