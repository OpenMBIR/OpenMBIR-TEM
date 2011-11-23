/* ============================================================================
 * Copyright (c) 2010, Michael A. Jackson (BlueQuartz Software)
 * Copyright (c) 2010, Dr. Michael A. Groeber (US Air Force Research Laboratories
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without modification,
 * are permitted provided that the following conditions are met:
 *
 * Redistributions of source code must retain the above copyright notice, this
 * list of conditions and the following disclaimer.
 *
 * Redistributions in binary form must reproduce the above copyright notice, this
 * list of conditions and the following disclaimer in the documentation and/or
 * other materials provided with the distribution.
 *
 * Neither the name of Michael A. Groeber, Michael A. Jackson, the US Air Force,
 * BlueQuartz Software nor the names of its contributors may be used to endorse
 * or promote products derived from this software without specific prior written
 * permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
 * USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *  This code was written under United States Air Force Contract number
 *                           FA8650-07-D-5800
 *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#ifndef _VTKFILEWRITERS_HPP_
#define _VTKFILEWRITERS_HPP_

#include <string>
#include <vector>

#include "MXA/Common/MXAEndian.h"
#include "MXA/Common/MXASetGetMacros.h"
#include "VTKWriterMacros.h"
#include "EIMTomo/ScaleOffsetMotionCorrectionAlgorithm/ScaleOffsetMotionStructures.h"
#include "EIMTomo/ScaleOffsetMotionCorrectionAlgorithm/ScaleOffsetMotionCorrectionConstants.h"


/**
 * @brief This is the SuperClass to implement if you want to write scalars to
 * a VTK Legacy File. The only method you need to implement is the writeScalars(FILE* f).
 */
class VtkScalarWriter
{
  public:
    VtkScalarWriter() : m_WriteBinaryFiles(true){}
    virtual ~VtkScalarWriter(){}

    bool m_WriteBinaryFiles;

    virtual int writeScalars(FILE* f)
    {
      int err = -1;
      return err;
    }

  protected:

  private:
    VtkScalarWriter(const VtkScalarWriter&); // Copy Constructor Not Implemented
    void operator=(const VtkScalarWriter&); // Operator '=' Not Implemented
};

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
typedef struct
{
    int xpoints;
    int ypoints;
    int zpoints;
    float resx;
    float resy;
    float resz;
}  DimsAndRes;

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------

class TomoOutputScalarWriter : public VtkScalarWriter
{
  public:
  TomoOutputScalarWriter(Geom* r) : r(r) {}
  virtual ~TomoOutputScalarWriter(){}

  int writeScalars(FILE* f)
  {
    int err = 0;
    std::string file;
    if (m_WriteBinaryFiles == true) {
      WRITE_VTK_FLOAT_VOXEL_BINARY(r, ScaleOffsetMotionCorrection::VTK::TomoVoxelScalarName, double);
    }
    else
    {
      WRITE_VTK_FLOAT_VOXEL_ASCII(r, ScaleOffsetMotionCorrection::VTK::TomoVoxelScalarName)
    }
    return err;
  }

  private:
    Geom* r;
    TomoOutputScalarWriter(const TomoOutputScalarWriter&); // Copy Constructor Not Implemented
    void operator=(const TomoOutputScalarWriter&); // Operator '=' Not Implemented
};



/**
 * @brief This macro can help you define the class needed to write a "Scalar" array
 * to the VTK file. This class is specifically setup for writing voxel based
 * properties to the VTK file
 */
#define VtkSCALARWRITER_CLASS_DEF(name, r, const_name, type, scalar, format)\
template<typename T>\
class name : public VtkScalarWriter\
{\
  public:\
    name(T* r) : r(r) {}\
    virtual ~name(){}\
  int writeScalars(FILE* f)  {\
    int err = 0;\
    std::string file;\
    size_t total = r->xpoints * r->ypoints * r->zpoints;\
    if (m_WriteBinaryFiles == true) {\
      WRITE_VTK_SCALARS_FROM_VOXEL_BINARY(r, const_name, type, scalar)\
    }    else    {\
      WRITE_VTK_SCALARS_FROM_VOXEL_ASCII(r, const_name, type, scalar, format)\
    }\
    return err;\
  }\
  private:\
    T* r;\
    name(const name&); \
    void operator=(const name&);\
};

#define VtkSCALARWRITER_CLASS_DEF_CHAR(name, r, const_name, type, scalar, format)\
template<typename T>\
class name : public VtkScalarWriter\
{\
  public:\
    name(T* r) : r(r) {}\
    virtual ~name(){}\
  int writeScalars(FILE* f)  {\
    int err = 0;\
    std::string file;\
    size_t total = r->xpoints * r->ypoints * r->zpoints;\
    if (m_WriteBinaryFiles == true) {\
      WRITE_VTK_SCALARS_FROM_VOXEL_BINARY_NOSWAP(r, const_name, type, scalar)\
    }    else    {\
      WRITE_VTK_SCALARS_FROM_VOXEL_ASCII(r, const_name, type, scalar, format)\
    }\
    return err;\
  }\
  private:\
    T* r;\
    name(const name&); \
    void operator=(const name&);\
};


/**
 * @class VTKRectilinearGridFileWriter VTKRectilinearGridFileWriter.h DREAM3D/Common/VTKUtils/VTKRectilinearGridFileWriter.h
 * @brief This is the main class to call when you want to write voxel based data
 * to a Rectilinear Grid based VTK data set. In order to write Scalar attributes
 * one needs to also include a vector of VtkScalarWriter objects. Each VtkScalarWriter
 * subclass knows how to write a specific scalar array to the file. When implementing
 * those subclasses keep in mind that you should be able to write both the ASCII
 * and Binary versions of your data to the VTK File.
 * @author Michael A. Jackson for BlueQuartz Software
 * @date Jun 13, 2011
 * @version 1.0
 */
class VTKRectilinearGridFileWriter
{
  public:
    //   MXA_SHARED_POINTERS(VTKRectilinearGridFileWriter);
    //   MXA_STATIC_NEW_MACRO(VTKRectilinearGridFileWriter);
    MXA_TYPE_MACRO(VTKRectilinearGridFileWriter);

    VTKRectilinearGridFileWriter() : m_WriteBinaryFiles(false) {}
    virtual ~VTKRectilinearGridFileWriter() {}

    MXA_INSTANCE_PROPERTY(bool, WriteBinaryFiles)

    /**
     * @brief This function writes a set of Axis coordinates to that are needed
     * for a Rectilinear Grid based data set.
     * @param f The "C" FILE* pointer to the file being written to.
     * @param axis The name of the Axis that is being written
     * @param type The type of primitive being written (float, int, ...)
     * @param npoints The total number of points in the array
     * @param min The minimum value of the axis
     * @param max The maximum value of the axis
     * @param step The step value between each point on the axis.
     */
    template<typename T>
    int writeCoords(FILE* f, const char* axis, const char* type, int npoints, T min, T max, T step)
    {
      int err = 0;
      fprintf(f, "%s %d %s\n", axis, npoints, type);
      if (m_WriteBinaryFiles == true)
      {
        T* data = new T[npoints];
        T d;
        for (int idx = 0; idx < npoints; ++idx)
        {
          d = idx * step + min;
          MXA::Endian::FromSystemToBig::convert<T>(d);
          data[idx] = d;
        }
        size_t totalWritten = fwrite(static_cast<void*>(data), sizeof(T), static_cast<size_t>(npoints), f);
        delete[] data;
        if (totalWritten != static_cast<size_t>(npoints) )
        {
          std::cout << "Error Writing Binary VTK Data into file " << std::endl;
          fclose(f);
          return -1;
        }

      }
      else
      {
        T d;
        for (int idx = 0; idx < npoints; ++idx)
        {
          d = idx * step + min;
          fprintf(f, "%f ", d);
          if (idx % 20 == 0 && idx != 0) { fprintf(f, "\n"); }
        }
        fprintf(f, "\n");
      }
      return err;
    }

    /**
     * @brief This is the method to call to write the Rectilinear Grid Data to the
     * legacy VTK data file
     * @param file The name of the file to write
     * @param r The Template Parameter for a one of the Plugin Func classes.
     * @param scalars The Vector of Scalars that need to also be written
     * @return Negative Value on error
     */
    template<typename T>
    int write(const std::string &file, T* r, std::vector<VtkScalarWriter*> scalars)
    {
      int err = 0;
      FILE* f = NULL;
      f = fopen(file.c_str(), "wb");
      if (NULL == f)
      {
        return -1;
      }
      // Write the correct header
      if (m_WriteBinaryFiles == true)
      {
        WRITE_RECTILINEAR_GRID_HEADER("BINARY", r, r->xpoints, r->ypoints, r->zpoints)
      }
      else
      {
        WRITE_RECTILINEAR_GRID_HEADER("ASCII", r, r->xpoints, r->ypoints, r->zpoints)
      }

      // Write the XCoords
      writeCoords(f, "X_COORDINATES", "float", r->xpoints + 1, 0.0f - r->resx * 0.5f, (float)(r->xpoints + 1 * r->resx), r->resx);
      writeCoords(f, "Y_COORDINATES", "float", r->ypoints + 1, 0.0f - r->resy * 0.5f, (float)(r->ypoints + 1 * r->resy), r->resy);
      writeCoords(f, "Z_COORDINATES", "float", r->zpoints + 1, 0.0f - r->resz * 0.5f, (float)(r->zpoints + 1 * r->resz), r->resz);

      size_t total = r->xpoints * r->ypoints * r->zpoints;
      fprintf(f, "CELL_DATA %d\n", (int)total);

      // Now loop on all of our Scalars and write those arrays as CELL_DATA
      for (typename std::vector<VtkScalarWriter*>::iterator iter = scalars.begin(); iter != scalars.end(); ++iter )
      {
        err = (*iter)->writeScalars(f);
        if (err < 0)
        {
          break;
        }
      }
      // Close the file
      fclose(f);
      return err;
    }


  private:
    VTKRectilinearGridFileWriter(const VTKRectilinearGridFileWriter&); // Copy Constructor Not Implemented
    void operator=(const VTKRectilinearGridFileWriter&); // Operator '=' Not Implemented
};


/**
 *
 */
class VTKStructuredPointsFileWriter
{
  public:
    VTKStructuredPointsFileWriter() {}
    MXA_TYPE_MACRO(VTKStructuredPointsFileWriter);

    virtual ~VTKStructuredPointsFileWriter(){}

    MXA_INSTANCE_PROPERTY(bool, WriteBinaryFiles)

    template<typename T>
    int write(const std::string &file, T* r, std::vector<VtkScalarWriter*> scalars)
    {
      int err = 0;
      FILE* f = NULL;
      f = fopen(file.c_str(), "wb");
      if (NULL == f)
      {
        return 1;
      }
      // Write the correct header
      if (m_WriteBinaryFiles == true)
      {
        WRITE_STRUCTURED_POINTS_HEADER("BINARY", r)
      }
      else
      {
        WRITE_STRUCTURED_POINTS_HEADER("ASCII", r)
      }

      //size_t total = r->xpoints * r->ypoints * r->zpoints;
      // Now loop on all of our Scalars and write those arrays as CELL_DATA
      for (typename std::vector<VtkScalarWriter*>::iterator iter = scalars.begin(); iter != scalars.end(); ++iter )
      {
        err = (*iter)->writeScalars(f);
        if (err < 0)
        {
          break;
        }
      }

      fclose(f);
      return err;
    }
  protected:


  private:
    VTKStructuredPointsFileWriter(const VTKStructuredPointsFileWriter&); // Copy Constructor Not Implemented
    void operator=(const VTKStructuredPointsFileWriter&); // Operator '=' Not Implemented
};



#endif /* _VTKFILEWRITERS_HPP_ */
