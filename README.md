# SlicerChai

This repository is still under development, and its initial code stems from the master's thesis "Haptization of Biomedical Volumetric Datasets in Scientific Visualization", authored by Darshain Desai and supervised by Lucas Siqueira Rodrigues at ZIB / Matters of Activity. 

# CMAKE -- Note
1. Ensure that you set the "CHAI3D_DIR" path.  This will be the path to your CHAI3D folder.  It should also contain the project file (CHAI3D-VS2015.vcxproj).


# chai3d -- Note

Go to C:\<>\chai3d-3.2.0\external\Eigen\Eigen\src\Core\Functors.h   \
Line 973 to line 979 \
Comment or remove below : \
//template<typename T> \
//struct functor_traits<std::binder2nd<T> > \
//{ enum { Cost = functor_traits<T>::Cost, PacketAccess = false }; }; \

//template<typename T> \
//struct functor_traits<std::binder1st<T> > \
//{ enum { Cost = functor_traits<T>::Cost, PacketAccess = false }; }; \



**The below may not be required now, but is retained in case any modifcation is required.** 



Go to C:\<>\chai3d\src\system\CThread.h   \
Around Line 134 add below: \
void start(void(*a_function)(cVector3d*), const CThreadPriority a_level, cVector3d* a_arg);
  
Go to C:\<>\chai3d\src\system\CThread.cpp   \
Around Line 178 add below: \
void cThread::start(void(*a_function)(cVector3d*), CThreadPriority a_level, cVector3d* a_arg)  \
{  \
    // create thread  \
#if defined(WIN32) | defined(WIN64)  \
    CreateThread( \ 
        0,  \  
        0,  \
        (LPTHREAD_START_ROUTINE)(a_function),  \
        a_arg,  \
        0,  \
        &m_threadId \
    );  \
#endif  \
  \
#if defined (LINUX) || defined (MACOSX)  \
    pthread_create(  \
        &m_handle,  \
        0,  \
        (void* (*)(void*)) a_function,  \
        a_arg  \
    );  \
#endif  \
 \
    // set thread priority level  \
    setPriority(a_level);  \
}   \
 


# Debugging -- Note
  
Slicer.exe --VisualStudio --launcher-no-splash --launcher-additional-settings C:\D\SE\SlicerHBVDinSciVisProd\build\AdditionalLauncherSettings.ini C:\D\SE\SlicerHBVDinSciVisProd\build\SlicerHBVDinSciVisProd.sln
  
In Visual Studio -> Go To -> Solution -> ALL_BUILD -> Properties -> Configuration Properties -> Debugging -> Command
Set the below:
C:\D\S5D\Slicer-build\bin\Debug\SlicerApp-real.exe  

In Slicer App -> Go To -> Edit -> Application Settings -> Modules -> Additional modulde paths -> Drag and Drop the folder as below :
C:\D\SE\SlicerHBVDinSciVisProd\build\lib\Slicer-5.x\qt-loadable-modules\Debug
  
Slicer.exe --VisualStudio --launcher-no-splash  

# 3DSystems Touch -- Note

Add Chai3D's hdPhantom64.dll to SlicerHBVDinSciVisProd\build\lib\Slicer-5.x\qt-loadable-modules\Debug (Chai3D fails to connect to the device if this dll is not placed there).


# Data -- Note

Ensure the file MR-head.nrrd or any data for loading is present in the folder hierarchy, as downloaded. In the parent directory of the build directory. [i.e where SlicerHBVDinSciVisProd.png is present.]


# VTK -- Note

The file **C:\D\S5D\VTK\Common\DataModel\vtkImageData.cxx**  has code to get image data's scalar pointer.  This errors out if the coordinates are not in the bounds. \
This code is indirectly called when we try to get scalar value in the GUI loop.  In case you want to temporarely supress this. Comment the error code in the below function. \

vtkImageData::GetScalarPointer(int coordinate[3])
