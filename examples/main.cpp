#ifndef CARDILLO_EXAMPLE_SCENE_HEADER
#error "CARDILLO_EXAMPLE_SCENE_HEADER must be defined by the build system"
#endif

#ifndef CARDILLO_EXAMPLE_SCENE_TYPE
#error "CARDILLO_EXAMPLE_SCENE_TYPE must be defined by the build system"
#endif

#include "example_main.hpp"

using namespace cardillo;

#include CARDILLO_EXAMPLE_SCENE_HEADER

using ExampleScene = CARDILLO_EXAMPLE_SCENE_TYPE;

CARDILLO_DEFINE_EXAMPLE_MAIN(ExampleScene)