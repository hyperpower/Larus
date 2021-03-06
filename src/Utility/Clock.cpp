/* Copyright 2010 Jukka Jylanki

 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at

 http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License. */

/** @file Clock.cpp
 @brief */

#if defined(__unix__) || defined(__native_client__) || defined(EMSCRIPTEN) || defined(ANDROID) || defined(__APPLE__) || defined (__CYGWIN__)
#include <time.h>
#include <errno.h>
#include <string.h>
#include <stdio.h>
#include <sys/time.h>
#endif

#ifdef WIN32
#include "../Math/InclWindows.h"
#endif

#ifdef EMSCRIPTEN
#include <emscripten.h>
#endif

#ifdef __APPLE__
#include <mach/mach_time.h>
#endif

#include "Clock.h"

namespace Larus
{

#ifdef WIN32
u64 Clock::ddwTimerFrequency;
#endif

#ifdef __APPLE__
tick_t Clock::ticksPerSecond = 0;
#endif

tick_t Clock::appStartTime = 0;

Clock impl;

void Clock::InitClockData()
{
	if (appStartTime == 0)
		appStartTime = Tick();

#ifdef WIN32
	if (!QueryPerformanceFrequency(reinterpret_cast<LARGE_INTEGER*>(&ddwTimerFrequency)))
	{
		LOGE("The system doesn't support high-resolution timers!");
		ddwTimerFrequency = (u64)-1;
	}

	if (appStartTime == 0)
	{
#if WINVER >= 0x0600 /* Vista or newer */ && !defined(MATH_ENABLE_WINXP_SUPPORT)
		appStartTime = (tick_t)GetTickCount64();
#else
// We are explicitly building with XP support, so GetTickCount() instead of GetTickCount64 is desired.
#if _MSC_VER >= 1700 // VS2012
#pragma warning(push)
#pragma warning(disable:28159) // warning C28159: Consider using 'GetTickCount64' instead of 'GetTickCount'. Reason: GetTickCount overflows roughly every 49 days.  Code that does not take that into account can loop indefinitely.  GetTickCount64 operates on 64 bit values and does not have that problem
#endif
		appStartTime = (tick_t)GetTickCount();
#if _MSC_VER >= 1700 // VS2012
#pragma warning(pop)
#endif

#endif
	}
#endif

#ifdef __APPLE__
	mach_timebase_info_data_t timeBaseInfo;
	mach_timebase_info(&timeBaseInfo);
	ticksPerSecond = 1000000000ULL * (uint64_t)timeBaseInfo.denom / (uint64_t)timeBaseInfo.numer;
	ASSERT(ticksPerSecond > (uint64_t)timeBaseInfo.denom/timeBaseInfo.numer);
	// Guard against overflow if OSX numer/denom change or similar.
#endif
}

Clock::Clock()
{
	InitClockData();
}

void Clock::Sleep(int milliseconds)
{
#ifdef WIN8RT
#pragma warning(Clock::Sleep has not been implemented!)
#elif defined(WIN32)
	::Sleep(milliseconds);
#elif !defined(__native_client__) && !defined(EMSCRIPTEN)
	// http://linux.die.net/man/2/nanosleep
	timespec ts;
	ts.tv_sec = milliseconds / 1000;
	ts.tv_nsec = (milliseconds - ts.tv_sec * 1000) * 1000 * 1000;
	int ret = nanosleep(&ts, NULL);
	if (ret == -1)
		printf("nanosleep returned -1! Reason: %s(%d).", strerror(errno),
				(int) errno); //revised
#else
#warning Clock::Sleep has not been implemented!
#endif
}

int Clock::Year()
{
#ifdef WIN32
	SYSTEMTIME s;
	GetSystemTime(&s);
	return s.wYear;
#else
	time_t theTime = time(NULL);
	struct tm *aTime = localtime(&theTime);
	return aTime->tm_year + 1900;
#endif
}

int Clock::Month()
{
#ifdef WIN32
	SYSTEMTIME s;
	GetSystemTime(&s);
	return s.wMonth;
#else
	time_t theTime = time(NULL);
	struct tm *aTime = localtime(&theTime);
	return aTime->tm_mon;
#endif
}

int Clock::Day()
{
#ifdef WIN32
	SYSTEMTIME s;
	GetSystemTime(&s);
	return s.wDay;
#else
	time_t theTime = time(NULL);
	struct tm *aTime = localtime(&theTime);
	return aTime->tm_mday;
#endif
}

int Clock::Hour()
{
#ifdef WIN32
	SYSTEMTIME s;
	GetSystemTime(&s);
	return s.wHour;
#else
	time_t theTime = time(NULL);
	struct tm *aTime = localtime(&theTime);
	return aTime->tm_hour;
#endif
}

int Clock::Min()
{
#ifdef WIN32
	SYSTEMTIME s;
	GetSystemTime(&s);
	return s.wMinute;
#else
	time_t theTime = time(NULL);
	struct tm *aTime = localtime(&theTime);
	return aTime->tm_min;
	return 0;
#endif
}

int Clock::Sec()
{
#ifdef WIN32
	SYSTEMTIME s;
	GetSystemTime(&s);
	return s.wSecond;
#else
	time_t theTime = time(NULL);
	struct tm *aTime = localtime(&theTime);
	return aTime->tm_sec;
	return 0;
#endif
}

unsigned long Clock::SystemTime()
{
#ifdef WIN32
#if WINVER >= 0x0600 /* Vista or newer */ && !defined(MATH_ENABLE_WINXP_SUPPORT)
	return (unsigned long)GetTickCount64();
#else
// We are explicitly building with XP support, so GetTickCount() instead of GetTickCount64 is desired.
#if _MSC_VER >= 1700 // VS2012
#pragma warning(push)
#pragma warning(disable:28159) // warning C28159: Consider using 'GetTickCount64' instead of 'GetTickCount'. Reason: GetTickCount overflows roughly every 49 days.  Code that does not take that into account can loop indefinitely.  GetTickCount64 operates on 64 bit values and does not have that problem
#endif
	return (unsigned long)GetTickCount();
#if _MSC_VER >= 1700 // VS2012
#pragma warning(pop)
#endif

#endif
#else
	return TickU32();
#endif
}

tick_t Clock::TimeElapseInTick()
{
	return (Tick() - appStartTime);
}

tick_t Clock::Tick()
{
	//android
#if defined(ANDROID)
	struct timespec res;
	clock_gettime(CLOCK_REALTIME, &res);
	return 1000000000ULL*res.tv_sec + (tick_t)res.tv_nsec;
	//android
#elif defined(EMSCRIPTEN)
#ifdef _TICK_IS_FLOAT
	return (tick_t)emscripten_get_now();
#else
	// emscripten_get_now() returns a wallclock time as a float in milliseconds (1e-3).
	// scale it to microseconds (1e-6) and return as a tick.
	return (tick_t)(((double)emscripten_get_now()) * 1e3);
#endif

#elif defined(WIN32)
	LARGE_INTEGER ddwTimer;
	BOOL success = QueryPerformanceCounter(&ddwTimer);
	assume(success != 0);
	MARK_UNUSED(success);
	return ddwTimer.QuadPart;
#elif defined(__APPLE__)
	return mach_absolute_time();
#elif defined(_POSIX_MONOTONIC_CLOCK)
	timespec t;
	clock_gettime(CLOCK_MONOTONIC, &t);
	return (tick_t)t.tv_sec * 1000 * 1000 * 1000 + (tick_t)t.tv_nsec;
#elif defined(_POSIX_C_SOURCE)
	timeval t;
	gettimeofday(&t, NULL);
	return (tick_t) t.tv_sec * 1000 * 1000 + (tick_t) t.tv_usec;
#else
	return (tick_t)clock();
#endif
}

unsigned long Clock::TickU32()
{
#ifdef WIN32
	LARGE_INTEGER ddwTimer;
	BOOL success = QueryPerformanceCounter(&ddwTimer);
	assume(success != 0);
	MARK_UNUSED(success);
	return ddwTimer.LowPart;
#else
	return (unsigned long) Tick();
#endif
}

tick_t Clock::TicksPerSec()
{
#if defined(ANDROID)
	return 1000000000ULL; // 1e9 == nanoseconds.
#elif defined(EMSCRIPTEN)

#ifdef _TICK_IS_FLOAT
	return (tick_t)1000.0;
#else
	return 1000000ULL; // 1e6 == microseconds.
#endif

#elif defined(WIN32)
	return ddwTimerFrequency;
#elif defined(__APPLE__)
	return ticksPerSecond;
#elif defined(_POSIX_MONOTONIC_CLOCK)
	return 1000 * 1000 * 1000;
#elif defined(_POSIX_C_SOURCE) || defined(__APPLE__)
	return 1000 * 1000;
#else
	return CLOCKS_PER_SEC;
#endif
}


}
