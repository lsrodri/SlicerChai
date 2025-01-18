cmake_minimum_required(VERSION 3.16.3...3.19.7 FATAL_ERROR)

project(SlicerHBVDinSciVisProd)
unset(ENV${CMAKE_C_COMPILER}PARENT_SCOPE)
unset(ENV${CMAKE_CXX_COMPILER}[PARENT_SCOPE])
set(CMAKE_C_COMPILER "C:/Program Files/Microsoft Visual Studio/2022/Community/VC/Tools/MSVC/14.37.32822/bin/Hostx64/x64/cl.exe")
set(CMAKE_CXX_COMPILER "C:/Program Files/Microsoft Visual Studio/2022/Community/VC/Tools/MSVC/14.37.32822/bin/Hostx64/x64/cl.exe")
message(WARNING "drd 1 CMAKE_Cxx_COMPILER ${CMAKE_CXX_COMPILER}")
message(WARNING "drd 2 CMAKE_C_COMPILER ${CMAKE_C_COMPILER}")

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
unset(ENV${CMAKE_C_COMPILER}PARENT_SCOPE)
unset(ENV${CMAKE_CXX_COMPILER}[PARENT_SCOPE])
message(WARNING "drd 3 CMAKE_Cxx_COMPILER ${CMAKE_CXX_COMPILER}")
message(WARNING "drd 4 CMAKE_C_COMPILER ${CMAKE_C_COMPILER}")
set(CMAKE_C_COMPILER "C:/Program Files/Microsoft Visual Studio/2022/Community/VC/Tools/MSVC/14.37.32822/bin/Hostx64/x64/cl.exe")
set(CMAKE_CXX_COMPILER "C:/Program Files/Microsoft Visual Studio/2022/Community/VC/Tools/MSVC/14.37.32822/bin/Hostx64/x64/cl.exe")

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


