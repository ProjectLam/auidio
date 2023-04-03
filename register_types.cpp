/* godot-cpp integration testing project.
 *
 * This is free and unencumbered software released into the public domain.
 */

#include "register_types.h"

#include "core/object/class_db.h"
#include "audio_effect_pitch_analyzer.h"

// void initialize_auidio_module(ModuleInitializationLevel p_level) {
// 	if (p_level != MODULE_INITIALIZATION_LEVEL_SCENE) {
// 		return;
// 	}

// 	ClassDB::register_class<AudioEffectPitchAnalyzer>();
// 	ClassDB::register_abstract_class<AudioEffectPitchAnalyzerInstance>();
// }

void register_auidio_types() {

	GDREGISTER_CLASS(AudioEffectPitchAnalyzer);
	GDREGISTER_ABSTRACT_CLASS(AudioEffectPitchAnalyzerInstance)
}
// void uninitialize_auidio_module(ModuleInitializationLevel p_level) {
// 	if (p_level != MODULE_INITIALIZATION_LEVEL_SCENE) {
// 		return;
// 	}
// }

void unregister_auidio_types() {

}

// extern "C" {
// GDNativeBool GDN_EXPORT auidio_library_init(const GDNativeInterface *p_interface, const GDNativeExtensionClassLibraryPtr p_library, GDNativeInitialization *r_initialization) {
// 	godot::GDExtensionBinding::InitObject init_obj(p_interface, p_library, r_initialization);

// 	init_obj.register_initializer(initialize_auidio_module);
// 	init_obj.register_terminator(uninitialize_auidio_module);
// 	init_obj.set_minimum_library_initialization_level(MODULE_INITIALIZATION_LEVEL_SCENE);

// 	return init_obj.init();
// }
// }

