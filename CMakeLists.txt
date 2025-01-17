cmake_minimum_required(VERSION 3.16.3...3.19.7 FATAL_ERROR)

project(SlicerHBVDinSciVisProd)

#-----------------------------------------------------------------------------
# Extension meta-information
set(EXTENSION_HOMEPAGE "https://www.slicer.org/wiki/Documentation/Nightly/Extensions/SlicerHBVDinSciVisProd")
set(EXTENSION_CATEGORY "Examples")
set(EXTENSION_CONTRIBUTORS "Dr Lucas Siqueira Rodrigues (Cluster of Excellence Matters of Activity)")
set(EXTENSION_DESCRIPTION "Connect to Haptic Device")
set(EXTENSION_ICONURL "https://www.example.com/Slicer/Extensions/SlicerHBVDinSciVisProd.png")
set(EXTENSION_SCREENSHOTURLS "https://www.example.com/Slicer/Extensions/SlicerHBVDinSciVisProd/Screenshots/1.png")
set(EXTENSION_DEPENDS "NA") # Specified as a list or "NA" if no dependencies

#-----------------------------------------------------------------------------
# Extension dependencies
find_package(Slicer REQUIRED)
include(${Slicer_USE_FILE})

#-----------------------------------------------------------------------------
# Extension modules
add_subdirectory(HBVDSciVisP)
## NEXT_MODULE

#-----------------------------------------------------------------------------
include(${Slicer_EXTENSION_GENERATE_CONFIG})
include(${Slicer_EXTENSION_CPACK})

#-----------------------------------------------------------------------------
# External Haptic dependencies
INCLUDE_EXTERNAL_MSPROJECT(CHAI3D ${CHAI3D_DIR}/CHAI3D-VS2015.vcxproj)
INCLUDE_EXTERNAL_MSPROJECT(glfw ${CHAI3D_DIR}/extras/GLFW/glfw-VS2015.vcxproj)