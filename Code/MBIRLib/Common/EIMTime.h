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
#ifndef _EIMTime_h_
#define _EIMTime_h_

#include "MBIRLib/MBIRLib.h"
#if CMP_HAVE_TIME_H
#include <time.h>
#endif

#if defined (CMP_HAVE_SYS_TIME_GETTIMEOFDAY)
#if CMP_HAVE_SYS_TIME_H
#include <sys/time.h>
#endif
#define CMP_GET_TIME_OF_DAY gettimeofday
#endif

#if defined (CMP_HAVE_TIME_GETTIMEOFDAY) || defined (_MSC_VER)
#if CMP_HAVE_TIME_H
#include <time.h>
#endif
#define CMP_GET_TIME_OF_DAY gettimeofday
#endif

#ifndef CMP_HAVE_TIME_GETTIMEOFDAY
#ifndef CMP_HAVE_SYS_TIME_GETTIMEOFDAY
#ifdef __MINGW32__
#ifndef GET_TIME_OF_DAY_DEFINITION
#define GET_TIME_OF_DAY_DEFINITION
//#warning no gettimeofday() function native to system. Using own implementation
//#endif

#include <windows.h>
#if CMP_HAVE_SYS_TIME_H
#include <sys/time.h>
#endif
void __stdcall GetSystemTimeAsFileTime(FILETIME*);

inline void CMP_gettimeofday(struct timeval* p, void* tz /* IGNORED */)
{
  union
  {
    long long ns100; /*time since 1 Jan 1601 in 100ns units */
    FILETIME ft;
  } now;

  GetSystemTimeAsFileTime( &(now.ft) );
  p->tv_usec = (long)((now.ns100 / 10LL) % 1000000LL );
  p->tv_sec = (long)((now.ns100 - (116444736000000000LL)) / 10000000LL);
}

#define CMP_GET_TIME_OF_DAY CMP_gettimeofday
#endif /*  GET_TIME_OF_DAY_WARNING */
#else
#error Your system does not have a gettimeofday() function. Please contact the author of this library for advice.
#endif /* MinGW case */
#endif /* HAVE_SYS_TIME_GETTIMEOFDAY */
#endif /* HAVE_TIME_GETTIMEOFDAY */

#ifdef _MSC_VER
#if _MSC_VER < 1400
#define TimeType __time64_t
#define TimeFunc _time64
#else
#define TimeType __time64_t
#define TimeFunc _time64
#endif
#else
#define TimeType time_t
#define TimeFunc time
#endif

#ifdef __cplusplus
extern "C" {
#endif

/**
* @brief Returns a Formatted String of the current Date/Time for logging
* purpose.
* @return A string of the current date/time
*/
char* EIMTOMO_LogTime();

/**
 * Returns a date/time string suitable for tiff tags.
 * @return
 */
char* EIMTOMO_TiffDateTime() ;

/**
 * @brief returns the number of milliseconds from a platform specified time.
 */
unsigned long long int EIMTOMO_getMilliSeconds();


#ifdef __cplusplus
}
#endif



#endif /* _EMTime.h_  */
