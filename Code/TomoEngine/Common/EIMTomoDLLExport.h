/* ============================================================================
 * Copyright (c) 2011, Michael A. Jackson (BlueQuartz Software)
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
 * Neither the name of Michael A. Jackson nor the names of its contributors may
 * be used to endorse or promote products derived from this software without
 * specific prior written permission.
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
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
#ifndef _TomoLib_DLL_EXPORT_H_
#define _TomoLib_DLL_EXPORT_H_


#if defined (_MSC_VER)
  #pragma warning(disable: 4251)
  #pragma warning(disable: 4710)
  #pragma warning(disable: 4820)
  #pragma warning(disable: 4668)
  #pragma warning(disable: 4265)
  #pragma warning(disable: 4189)
  #pragma warning(disable: 4640)
  #pragma warning(disable: 4996)
  #pragma warning(disable: 4548)
#endif

/* Cmake will define TomoLib_EXPORTS on Windows when it
configures to build a shared library. If you are going to use
another build system on windows or create the visual studio
projects by hand you need to define TomoLib_EXPORTS when
building the TomoLib DLL on windows.
*/

#if defined (TomoLib_BUILT_AS_DYNAMIC_LIB)

  #if defined (EIMTomoLib_EXPORTS)  /* Compiling the MXA DLL/Dylib */
    #if defined (_MSC_VER)  /* MSVC Compiler Case */
      #define  TomoLib_EXPORT __declspec(dllexport)
    #elif (__GNUC__ >= 4)  /* GCC 4.x has support for visibility options */
      #define TomoLib_EXPORT __attribute__ ((visibility("default")))
    #endif
  #else  /* Importing the DLL into another project */
    #if defined (_MSC_VER)  /* MSVC Compiler Case */
      #define  TomoLib_EXPORT __declspec(dllimport)
    #elif (__GNUC__ >= 4)  /* GCC 4.x has support for visibility options */
      #define TomoLib_EXPORT __attribute__ ((visibility("default")))
    #endif
  #endif
#endif

/* If TomoLib_EXPORT was never defined, define it here */
#ifndef TomoLib_EXPORT
  #define TomoLib_EXPORT
#endif

#if 0
#if defined (_WIN32) || defined __CYGWIN__

  #if defined (TomoLib_BUILT_AS_DYNAMIC_LIB)
    #if defined(EIMTomoLib_EXPORTS)
      #define  TomoLib_EXPORT __declspec(dllexport)
    #else
      #define  TomoLib_EXPORT __declspec(dllimport)
    #endif /* EIMTomoLib_EXPORTS */
  #else
    #define TomoLib_EXPORT
  #endif
#elif __GNUC__ >= 4
 #define FLOW_DLL __attribute__ ((visibility("default")))
 #define DLL_LOCAL  __attribute__ ((visibility("hidden")
#else /* defined (_WIN32) && defined (TomoLib_BUILD_SHARED_LIBS)  */
 #define TomoLib_EXPORT
#endif
#endif

#endif /* _TomoLib_DLL_EXPORT_H_ */
