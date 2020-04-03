# OpenMBIR-TEM #

## Introduction ##

OpenMBIR-TEM is an open source project that develops model-based imaging algorithms for electron tomography. Model-based approaches to imaging combine the physics of data formation, the noise characteristics of the sensors and a probabilistic model for the object to formulate the image reconstruction problem. In several applications MBIR translates to minimizing a very high-dimensional cost-function which requires efficient algorithms. 

This package contains source codes to reconstruct STEM HAADF and BrightField tomographic data using the Model Based Iterative Reconstruction Method. Newer versions now available will also reconstruct tomography data stored in the [MRC](http://bio3d.colorado.edu/imod/doc/mrc_format.txt) file format. A number of output files are written including a new MRC file that is a regular grid of stacked data instead of a tilt series, a Uniform Coordinate file for ParaView in the form of a vtk file and also an Aviso .am file of Uniform Coordinate type. More information about the MBIR algorithm can be found at [OpenMBIR.org](http://www.openmbir.org)

## License ##

The source codes are open source under a BSD License. See the License.txt file for more information.

## Dependencies ##
The current set of dependencies are:

+ Compiler appropriate for your operating system
+ [CMake version 2.8.10 or newer](http://www.cmake.org/cmake/resources/software.html)
+ [Boost version 1.44 or newer](http://www.boost.org)
+ [Qt version 4.8.5 or newer but NOT Qt v5](http://qt-project.org)
+ [Threading Building Blocks](https://www.threadingbuildingblocks.org/download)

## Known Compilers ##

| Operating System | Compiler | Notes |  
| ------------------------|--------------|---------|  
| OS X 10.7/8/9 | Xcode  4.6.x & Xcode 5.x | OS X 10.9 requires Qt 4.8.6 |  
| Windows x64 | VS 2010, VS2012,VS2013 |Static Libraries ONLY  |  
| Linux x64 | GCC or Clang recent versions | Ubuntu 12.04 or above, CentOS 6.5 or greater |

## Known Issues ##

### OS X 10.9 "Mavericks" ###

Due to some bugs in Qt 4.8.5 the codes will NOT compile on OS X 10.9 Mavericks. The developer will need to download version 4.8.6 (Not available as of APRIL 15 2014) or download Qt 4.8.x from the Qt git repository and self build Qt.

Please see the web site [OpenMBIR WebSite](http://www.openmbir.org) for more information

## General Notes ##

This software is under development and may or may not give you the results that you are expecting. If you suspect a bug in the code please post a new issue to the OpenMBIR issue tracker at github.com
