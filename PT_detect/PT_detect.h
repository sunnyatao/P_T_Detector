#pragma once
#include "DF_filter.h"
#include "fdacoefs_BP_fs250_0_2_18_20.h"
#include "fdacoefs_BP_fs256_0_2_18_20.h"
#include "fdacoefs_FIR_ord100_HL.h" 
#include "fdacoefs_FIR_db4_scale18.h"
#include <math.h>
#include "stdlib.h"

//void PT_longterm(const float* signal, const int32_T signalsSize, const int32_T* Rposes, const bool* Noisetag, const int32_T RposesSize, const int fs, int32_T *Tposes, int32_T * Tnum, int32_T *Pposes, int32_T *Pnum);
void PT_detect(const float* signal, const int32_T signalsSize, const int32_T* Rposes,  const int32_T RposesSize, const int32_T fs, int32_T *Tposes, int32_T * Tnum, int32_T *Pposes, int32_T *Pnum, float * rr_noise_levels, int32_T start_idx);
void findTpeaks(const float* signal, const  int32_T signalsSize, const  int32_T* Rposes, const int32_T RposesSize, const int fs, int32_T *Tposes, int32_T * Tnum);
void findPpeaks(const float* signal, const  int32_T signalsSize, const  int32_T* Rposes, const int32_T RposesSize, const int fs, int32_T *Tposes, int32_T * Tnum, int32_T *Pposes, int32_T * Pnum, float * rr_noise_levels, int32_T start_idx);
template <class T> T mean(T* data, int num);
template <class T> T maxValue(T* data, int num);
template <class T> T minValue(T* data, int num);

#ifdef _WIN32
extern "C"
{
	__declspec(dllexport) void PT_longterm(const float* signal, const int32_T signalsSize, const int32_T* Rposes,const INT16_T* qrs_start_ends, const bool* Noisetag, const int32_T RposesSize, const int fs, int32_T *Tposes, int32_T * Tnum, int32_T *Pposes, int32_T *Pnum, float * rr_noise_levels);
}

#else
#ifdef __cplusplus

extern "C"
{
#endif

	void PT_longterm(const float* signal, const int32_T signalsSize, const int32_T* Rposes, const INT16_T * qrs_start_ends, const bool* Noisetag, const int32_T RposesSize, const int fs, int32_T *Tposes, int32_T * Tnum, int32_T *Pposes, int32_T *Pnum, float * rr_noise_levels);

#ifdef __cplusplus
}
#endif 

#endif