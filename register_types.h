/* godot-cpp integration testing project.
 *
 * This is free and unencumbered software released into the public domain.
 */

#ifndef EXAMPLE_REGISTER_TYPES_H
#define EXAMPLE_REGISTER_TYPES_H

#include "modules/register_module_types.h"


// void initialize_auidio_module(ModuleInitializationLevel p_level);
// void uninitialize_auidio_module(ModuleInitializationLevel p_level);

void register_auidio_types();
void unregister_auidio_types();

#endif // EXAMPLE_REGISTER_TYPES_H
