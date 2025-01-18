# SlicerChai

**This repository is still under development.**

This repository extends the initial code stemming from "Haptization of Biomedical Volumetric Datasets in Scientific Visualization", a master's thesis authored by Darshain Desai and supervised by Lucas Siqueira Rodrigues at ZIB / Matters of Activity. 

This extension integrates Chai3D into 3DSlicer as a loadable module extension and provides a device-agnostic integration of haptics into 3D scenes. Haptic cursors are represented into 3D scenes, but voxel-based haptics are still under development. Limited interaction has been added.

# Chai3D -- Note
We have created a compatible version of Chai3D 3.2.0 at https://github.com/lsrodri/chai3d-lite-slicer. Your local clone of this repository should be used in your CMake build. 

# CMake -- Note
1. Set "CHAI3D_DIR" to your local directory. CMake will extract dependencies from CHAI3D-VS2015.vcxproj.
2. Set "Slicer_DIR" to your local Slicer-build directory. 

# Debugging -- Note
  
Slicer.exe --VisualStudio --launcher-no-splash --launcher-additional-settings <Build Folder>\AdditionalLauncherSettings.ini <Build Folder>\SlicerHBVDinSciVisProd.sln
  
In Visual Studio -> Go To -> Solution -> ALL_BUILD -> Properties -> Configuration Properties -> Debugging -> Command
Set the below:
<Slicer Folder>\Slicer-build\bin\Debug\SlicerApp-real.exe  

In Slicer App -> Go To -> Edit -> Application Settings -> Modules -> Additional modulde paths -> Drag and Drop the folder as below :
<Build Folder>\lib\Slicer-5.x\qt-loadable-modules\Debug
  
Slicer.exe --VisualStudio --launcher-no-splash  

# 3DSystems Touch -- Note

Add Chai3D's hdPhantom64.dll to SlicerHBVDinSciVisProd\build\lib\Slicer-5.x\qt-loadable-modules\Debug (Chai3D fails to connect to the device if this dll is not placed there).


# Data -- Note

Ensure the file MR-head.nrrd or any data for loading is present in the folder hierarchy, as downloaded. In the parent directory of the build directory. [i.e where SlicerHBVDinSciVisProd.png is present.]


# VTK -- Note

The file **<Slicer Folder>\VTK\Common\DataModel\vtkImageData.cxx**  has code to get image data's scalar pointer.  This errors out if the coordinates are not in the bounds. \
This code is indirectly called when we try to get scalar value in the GUI loop.  In case you want to temporarely supress this. Comment the error code in the below function. \

vtkImageData::GetScalarPointer(int coordinate[3])
