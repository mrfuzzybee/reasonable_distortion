//-------------------------------------------------------------------------------------------------------
// VST Plug-Ins SDK
// Version 2.4		$Date: 2006/11/13 09:08:27 $
//
// Category     : VST 2.x SDK Samples
// Filename     : adelay.cpp
// Created by   : Steinberg Media Technologies
// Description  : Simple Delay plugin (Mono->Stereo)
//
// ｩ 2006, Steinberg Media Technologies, All Rights Reserved
//-------------------------------------------------------------------------------------------------------

#include <stdio.h>
#include <string.h>

#ifndef __adelay__
#include "distoo.h"
#endif

#include <Windows.h>
#include <stdint.h>

#define MAX(A,B) (((A)>(B))?(A):(B))
#define MIN(A,B) (((A)<(B))?(A):(B))

//#define REF

#define OPT




double han_window(double s) {
	return (cos(s / SINC_HALF_LEN)+ 1.0) * 0.5;
}

double blackman_window(double s) {
	return 0.42 + 0.5 * cos(s / SINC_HALF_LEN) + 0.08 * cos(s / SINC_HALF_LEN * 2);
}
 

#define MIN_ERROR (1.0e-7)

double sinc(double s) {
	if (fabs(s) < MIN_ERROR) return 1.0;
	return sin(s) / s;
}

double lanczos_window(double s) 
{
	return sinc(s / SINC_HALF_LEN * 2);
}


double han_flat_window(double s) 
{
	s = fabs(s);
	if (s < SINC_HALF_LEN / 2 * M_PI) 
	{
		return 1.0;
	}
	s -=  SINC_HALF_LEN / 2 * M_PI;


	return (cos(s * 2 / SINC_HALF_LEN) + 1.0) * 0.5;
}



double sinc_win(double s) {
	
//NaNってことないとは思うけど

	s *= M_PI;
	//阻止周波数がぎりぎりだとエイリアスを防ぎきれないのでちょっと下げるおまじない
	double d = s * (SINC_HALF_LEN - 1) / SINC_HALF_LEN;

	double ret = sinc(d) * blackman_window(s);
	return (fabs(ret) < MIN_ERROR) ? 0 : ret;
	 
}


double Distoo::makeSincBuf(float* buffer) {

	//左右対称ではなくて1サンプル目に0が入ってます
	//真ん中のサンプル(1.0が値として入る）はbuffer[SINC_HALF_LEN * OVER_SAMPLE]
	
	double sum = 0;
	for (int i = 0; i < OVER_SAMPLED_SINC_LEN; i++) {
		sum += buffer[i] = sinc_win(
			(i - SINC_HALF_LEN * OVER_SAMPLE) / (double)OVER_SAMPLE
		);
	}
	/*
	FILE* fp =fopen("c:\\tsettest.pcm","wb");
	fwrite(buffer, sizeof(float), OVER_SAMPLED_SINC_LEN, fp);
	fclose(fp);
	*/



	return sum;
}

//入力側（アップサンプリング用)のSinc関数のバッファを生成する
double Distoo::makeSincBufForInput(float* buffer) {
	float* temp = new float[OVER_SAMPLED_SINC_LEN];
	double sum = makeSincBuf(temp);

#if defined(OPT)
	for (int i = 0; i < OVER_SAMPLE; i++) {
		for (int j = 0; j < SINC_FULL_LEN; j++) {
				int j2 = (j + SINC_FULL_LEN - 1) % SINC_FULL_LEN;
				for (int k = 0; k < 4; k++) {//SSE用にシフト
					int shift_sinc = (j2 + k) % SINC_FULL_LEN;
					float sinc = temp[(OVER_SAMPLE - i - 1) + shift_sinc * OVER_SAMPLE];
					int index =   (j + i * SINC_FULL_LEN * 2) 
								+ (OVER_SAMPLED_SINC_LEN * 2 * k);
					buffer[index] = sinc;
					buffer[index+ SINC_FULL_LEN] = sinc;
				}
		}
	}
#else 
	for (int i = 0; i < OVER_SAMPLED_SINC_LEN; i++) {
		float sinc = temp[i];
		buffer[i] = sinc;
		buffer[i + OVER_SAMPLED_SINC_LEN] = sinc;
	}
#endif

	delete [] temp;


	return sum;
}

//出力側（ダウンサンプリング用)のSinc関数のバッファを生成する
double Distoo::makeSincBufForOutput(float* buffer) {

	float* temp = new float[OVER_SAMPLED_SINC_LEN];
	double sum = makeSincBuf(temp);

#if  defined(OPT)
	for (int i = 0; i < OVER_SAMPLED_SINC_LEN; i++) {
		float sinc = temp[(i + OVER_SAMPLED_SINC_LEN - OVER_SAMPLE) % OVER_SAMPLED_SINC_LEN];
		buffer[i] = sinc;
		buffer[i + OVER_SAMPLED_SINC_LEN] = sinc;
	}
#else
	for (int i = 0; i < OVER_SAMPLED_SINC_LEN; i++) {
		float sinc = temp[i];
		buffer[i] = sinc;
		buffer[i + OVER_SAMPLED_SINC_LEN] = sinc;
	}
#endif
	delete [] temp;

	return sum;
}

#include <xmmintrin.h>
#include <intrin.h>







//----------------------------------------------------------------------------- 
ADistoooProgram::ADistoooProgram ()
{
	// default Program Values
	fType = 0;
	fGain = 0.0;
	fBias = 0.0;
	fEdge = 0.0;
	fMaster = 1.0;
	

	strcpy (name, "Init");
}


static float* mallc(int size) {
	return (float*)_aligned_malloc(size * sizeof(float), 16); 
}

static void fre(float* p)
{
	if (p != NULL) {
		_aligned_free(p);
	}
}

#include "mmsystem.h"

	
//#define WAVEFORM_TEST

#ifdef WAVEFORM_TEST
FILE* fp_up;

#endif

//-----------------------------------------------------------------------------
Distoo::Distoo (audioMasterCallback audioMaster)
	: AudioEffectX (audioMaster, kNumPrograms, kNumParams)
{
	// init
	cursor = 0;

	programs = new ADistoooProgram[numPrograms];

	programs[0].fType = kDistortion / (double)kNumType;
	programs[0].fBias = 0.5;
	programs[0].fGain = 0.16;
	programs[0].fEdge = 0.14;
	programs[0].fMaster = 0.8;
	strcpy (programs[0].name, "Mild Saturation");

	programs[1].fType = kDistortion / (double)kNumType;
	programs[1].fBias = 0.5;
	programs[1].fGain = 1.0;
	programs[1].fEdge = 1.0;
	programs[1].fMaster = 0.8;
	strcpy (programs[1].name, "Full Drive");

	programs[2].fType = kDistortion / (double)kNumType;
	programs[2].fBias = 0.0;
	programs[2].fGain = 1.0;
	programs[2].fEdge = 0.01;
	programs[2].fMaster = 0.5;
	strcpy (programs[2].name, "Powered Attack");

	programs[3].fType = kOverdrive / (double)kNumType;
	programs[3].fBias = 0.8;
	programs[3].fGain = 1.0;
	programs[3].fEdge = 0.35;
	programs[3].fMaster = 0.85;
	strcpy (programs[3].name, "Hard Rock");

	programs[4].fType = kOverdrive / (double)kNumType;
	programs[4].fBias = 1.0;
	programs[4].fGain = 0.25;
	programs[4].fEdge = 0.1;
	programs[4].fMaster = 1.0;
	strcpy (programs[4].name, "Smoky Tube");

	programs[5].fType = kSpike / (double)kNumType;
	programs[5].fBias = 1.0;
	programs[5].fGain = 0.2;
	programs[5].fEdge = 0.1;
	programs[5].fMaster = 0.95;
	strcpy (programs[5].name, "Tingly");

	programs[6].fType = kSpike / (double)kNumType;
	programs[6].fBias = 0.7;
	programs[6].fGain = 1.0;
	programs[6].fEdge = 0.05;
	programs[6].fMaster = 0.9;
	strcpy (programs[6].name, "Yet Another Drive");

	programs[7].fType = kWater / (double)kNumType;
	programs[7].fBias = 0.2;
	programs[7].fGain = 1.0;
	programs[7].fEdge = 1.00;
	programs[7].fMaster = 0.75;
	strcpy (programs[7].name, "Buggy One");

	programs[8].fType = kLPPin / (double)kNumType;
	programs[8].fBias = 0.2;
	programs[8].fGain = 1.0;
	programs[8].fEdge = 1.00;
	programs[8].fMaster = 0.8;
	strcpy (programs[8].name, "LP Pin");


	if (programs)
		setProgram (0);



	INIT_LPF_PARAM(fGain)
	INIT_LPF_PARAM(fBias)
	INIT_LPF_PARAM(fEdge)
	INIT_LPF_PARAM(fMaster)

	bufferi[0] = mallc(SINC_FULL_LEN * 2); 
	bufferi[1] = bufferi[0] + SINC_FULL_LEN; 
	buffero[0] = mallc(OVER_SAMPLED_SINC_LEN * 2);
	buffero[1] = buffero[0] + OVER_SAMPLED_SINC_LEN; 

	
	




#ifdef REF
	bufferSinc = mallc(OVER_SAMPLED_SINC_LEN);
	invSincArea = 1.0 / makeSincBuf(bufferSinc);
#endif

#if  defined(OPT)
	bufferSinc = mallc(OVER_SAMPLED_SINC_LEN * 10);
	bufferSincOut = bufferSinc + OVER_SAMPLED_SINC_LEN * 8;
	invSincArea = 1.0 / makeSincBufForInput(bufferSinc);
	makeSincBufForOutput(bufferSincOut);
#endif


	setNumInputs (2);	// stereo input
	setNumOutputs (2);	// stereo output

	fMaster = 0.5/OVER_SAMPLE;
	setUniqueID ('cj9!');	// this should be unique, use the Steinberg web page for plugin Id registration

	setInitialDelay(SINC_FULL_LEN);

	sampleRate_per441 = 44100.0f / getSampleRate();



	//iir[0] = new IIR(fEdgeLog, LPF_Q , getSampleRate());
	//iir[1] = new IIR(fEdgeLog, LPF_Q , getSampleRate());

	/*
		//パフォーマンステスト
	VstInt32 sampleFrames = 1024;
	float* inputs[2];
	float* outputs[2];
	inputs[0] = new float[sampleFrames * 4];
	inputs[1] = inputs[0] + sampleFrames;
	outputs[0]= inputs[0] + sampleFrames * 2;
	outputs[1]= inputs[0] + sampleFrames * 3;



	FILE* fp =fopen("c:\\parformtest.csv","w");
	
	for (int k = 0; k < 10; k++) {
		DWORD start = timeGetTime();
		for (int i = 0; i < 100; i++) {
			processReplacing (inputs, outputs, sampleFrames);
		}
		DWORD end =  timeGetTime();
		fprintf(fp, "%d\n",(int)( end - start));
	}
	fclose(fp);
	delete [] inputs[0];
	
	*/


#ifdef WAVEFORM_TEST
	//アップサンプリングの波形確認テスト
	VstInt32 sampleFrames = 1024;
	float* inputs[2];
	float* outputs[2];
	inputs[0] = new float[sampleFrames * 4];
	inputs[1] = inputs[0] + sampleFrames;
	outputs[0]= inputs[0] + sampleFrames * 2;
	outputs[1]= inputs[0] + sampleFrames * 3;


	fp_up = fopen("c:\\hakei_sine1_up.pcm","wb");
	resume ();		// flush buffer
	for (int i = 0; i < sampleFrames; i++) {
		inputs[0][i] = inputs[1][i] = sin(i * 3.141592 * 5 / sampleFrames);
	}
	processReplacing (inputs, outputs, sampleFrames);
	fclose(fp_up);

	resume ();		// flush buffer
	fp_up = fopen("c:\\hakei_pulse_up.pcm","wb");
	for (int i = 0; i < sampleFrames; i++) {
		inputs[0][i] = inputs[1][i] = (i == 0) ? 1 : 0;
	}
	processReplacing (inputs, outputs, sampleFrames);
	fclose(fp_up);
	fp_up = NULL;

	delete [] inputs[0];
#endif


	resume ();		// flush buffer
}

//------------------------------------------------------------------------
Distoo::~Distoo ()
{
	fre(bufferi[0]);
	fre(buffero[0]);
	fre(bufferSinc);
	//delete iir[0];
	//delete iir[1];

	if (programs)
		delete[] programs;

}

//------------------------------------------------------------------------
void Distoo::setProgram (VstInt32 program)
{
	ADistoooProgram* ap = &programs[program];

	curProgram = program;
	setParameter (kType, ap->fType);	
	setParameter (kGain, ap->fGain);	
	setParameter (kBias, ap->fBias);
	setParameter (kEdge, ap->fEdge);
	setParameter (kMaster, ap->fMaster);
}



//------------------------------------------------------------------------
void Distoo::setProgramName (char *name)
{
	strcpy (programs[curProgram].name, name);
}

//------------------------------------------------------------------------
void Distoo::getProgramName (char *name)
{

		strcpy (name, programs[curProgram].name);
}

//-----------------------------------------------------------------------------------------
bool Distoo::getProgramNameIndexed (VstInt32 category, VstInt32 index, char* text)
{
	if (index < kNumPrograms)
	{
		strcpy (text, programs[index].name);
		return true;
	}
	return false;
}

//------------------------------------------------------------------------
void Distoo::resume ()
{
	memset (bufferi[0], 0, SINC_FULL_LEN * 2 *  sizeof (float));
	memset (buffero[0], 0, OVER_SAMPLED_SINC_LEN * 2 *  sizeof (float));

	cursor = 0;

	for (int i = 0; i < 12; i++) {
		d_gain_old[i] = 0;
	}
	for (int i = 0; i < 4; i++) {
		i_spike_old[i] = 0;
	}

	//iir[0]->resume();
	//iir[1]->resume();


	AudioEffectX::resume ();
}

void Distoo::setType(float value) {
	iType = min(value * kNumType+0.00001, kNumType - 1);
	
}
float Distoo::getType() {
	return iType / (float)kNumType;
}


void Distoo::setBias(float value) {
	fBiasParam = value;
	fBias = shiftParam(fBiasParam, 0.000001);
}
float Distoo::getBias() {
	return fBiasParam;
}



void Distoo::setGain(float value) {
	fGainParam = value;
	fGain = shiftParam(fGainParam, 0.000001);
}
float Distoo::getGain() {
	return fGainParam;
}


void Distoo::setEdge(float value) {
	fEdgeParam = value;
	fEdge = shiftParam(fEdgeParam, 0.01);
}
float Distoo::getEdge() {
	return fEdgeParam;
}


void Distoo::setMaster(float value) {
	fMasterParam = value;
	if (fMasterParam > 0.1) {
		//fMasterParam=0.1のときfMaster=0.01(-60dB)
		fMaster = (float)pow(10.0, -(1.0 - (fMasterParam - 0.1) / 0.9) * 3.0); 
	} else {
		double e = 0.0000001;
		fMaster = (float)fMasterParam * (0.01 - e) + e;
	}
}
float Distoo::getMaster() {
	return fMasterParam;
}


//------------------------------------------------------------------------
void Distoo::setParameter (VstInt32 index, float value)
{
	ADistoooProgram* ap = &programs[curProgram];

	switch (index)
	{
		case kType:     setType(ap->fType = value);	    break;
		case kGain :    setGain(ap->fGain = value);		break;
		case kBias :	setBias(ap->fBias = value); 	break;
		case kEdge :	setEdge(ap->fEdge = value); 	break;
		case kMaster :  setMaster(ap->fMaster = value); break;
	}
}

//------------------------------------------------------------------------
float Distoo::getParameter (VstInt32 index)
{
	float v = 0;

	switch (index)
	{
		case kType:     v = getType();	   break;
		case kGain :    v = getGain();     break;
		case kBias :	v = getBias();     break;
		case kEdge :	v = getEdge();     break;
		case kMaster :  v = getMaster();   break;
	}
	return v;
}

//------------------------------------------------------------------------
void Distoo::getParameterName (VstInt32 index, char *label)
{
	switch (index)
	{
		case kType :    strcpy (label, "Type");		break;
		case kGain :    strcpy (label, "Gain");		break;
		case kBias :    strcpy (label, "Bias");		break;
		case kEdge :    strcpy (label, "Edge");		break;
		case kMaster :  strcpy (label, "Master");	break;
	}
}

//------------------------------------------------------------------------
void Distoo::getParameterDisplay (VstInt32 index, char *text)
{
	switch (index)
	{
		case kType  :   getTypeName(text);											break;
		case kGain  :   float2string (getGain()*100,  text, kVstMaxParamStrLen);	break;
		case kBias  :	float2string (getBias()*100,  text, kVstMaxParamStrLen);	break;
		case kEdge :	float2string (getEdge()*100,  text, kVstMaxParamStrLen);	break;
		case kMaster :	dB2string    (fMaster,        text, kVstMaxParamStrLen);	break;
	}
}

//------------------------------------------------------------------------
void Distoo::getParameterLabel (VstInt32 index, char *label)
{
	switch (index)
	{
		case kType :		strcpy (label, "Type");		break;
		case kGain :		strcpy (label, "Gain");		break;
		case kBias :		strcpy (label, "Bias");		break;
		case kEdge :		strcpy (label, "Edge");		break;
		case kMaster :      strcpy (label, "dB");		break;
	}
}




void Distoo::getTypeName(char* text)
{
	static const char* typeName[] = {
		"Standard",
		"Tube",
		"Frizz",
		"Water",
		"LP_Pin"
	};
	if (iType >= 0 && iType < kNumType) {
		vst_strncpy(text, typeName[iType], kVstMaxParamStrLen);
	} else {
		text[0] = 0;
	}
}



//------------------------------------------------------------------------
bool Distoo::getEffectName (char* name)
{
	strcpy (name, "Reasonable Distortion");
	return true;
}

//------------------------------------------------------------------------
bool Distoo::getProductString (char* text)
{
	strcpy (text, "Reasonable Distortion");
	return true;
}

//------------------------------------------------------------------------
bool Distoo::getVendorString (char* text)
{
	strcpy (text, "Hanacsoft");
	return true;
}


//---------------------------------------------------------------------------


inline void smooth(double in, double& d1, double&d2, double edge)
{
	d1 +=  (in - d1) * edge;
	d2 +=  (d1 - d2) * edge;
}

inline static void distortion(float gain, float bias, float master, float*pL , float* pR, int size, double* d_gain_old, float edge, float sampleRate_per441) {

#ifdef WAVEFORM_TEST 
	if (fp_up) {
			fwrite(pL, sizeof(float), size, fp_up);
	}
#endif

	for (int j = 0; j < 2; j++) 
	{
		float *pData = (j == 0) ? pL : pR;


		for (int i = 0; i < size; i++) {
			float data = pData[i];
			float d_gain, distdata;

			//符号を取る
			unsigned int* iip = (unsigned int*)&data;
			unsigned int sign = (*iip) & (1 << 31);
			(*iip) -= sign;

			//歪の主処理
			d_gain = gain * (1.0f - (sign >> 31) * bias);
			distdata = d_gain * data;
			if (distdata > 1.0f) {
				d_gain = 1.0f / (data );
			} 

			smooth(d_gain, d_gain_old[j], d_gain_old[j+2], edge * edge * sampleRate_per441);
			data *= d_gain_old[j+2] * master;

			data = (data > 0.99f) ? 0.99f : data;
			


			//符号を戻す
			(*iip) += sign;

			pData[i] = data;

		}
	}
}


inline static void overdrive(float gain, float bias, float master, float*pL , float* pR, int size,
	double* d_gain_old, float edge, int* i_spike_old, float sampleRate_per441) {

	
	float attack  = 0.8f  * sampleRate_per441;//(edge * edge * 0.4) * sampleRate_per441;
	float release = 0.008 * sampleRate_per441 ;//(edge * edge * 0.03) * sampleRate_per441;

	float shift = 0.0001f;
	float powShift = powf(shift , gain);
	float invPow = 1.0f / (powf(gain + shift , gain) - powShift);

	float negGain = (gain * (1.0f - bias) + bias);
	float origGain = 1.2f / (gain + 0.2f);
	float origGainNeg = 1.2f / (negGain + 0.2f) - origGain;
	

	for (int j = 0; j < 2; j++) 
	{
		float *pData = (j == 0) ? pL : pR;

		for (int i = 0; i < size; i++) {
			float data = pData[i];
			float my_gain, distdata;
			//符号を取る
			unsigned int* iip = (unsigned int*)&data;
			unsigned int sign = (*iip) & (1 << 31);
			(*iip) -= sign;

			//歪の主処理
			my_gain = origGain +  (sign >> 31) * origGainNeg;


			distdata = (powf(data * my_gain + shift , gain) - powShift) * invPow;
			
			if (distdata > 0.95 && i_spike_old[j] == 0) {
				i_spike_old[j] = 1;
				d_gain_old[j] = d_gain_old[j+2] = 0.95f;
			}
			if (d_gain_old[j+2] > distdata) {
				i_spike_old[j] = 0;
				d_gain_old[j] = d_gain_old[j+2] =  0.95f;
			}
			if (i_spike_old[j+2] != sign) {
				i_spike_old[j] = 0;
				i_spike_old[j+2] = sign;
				d_gain_old[j] = d_gain_old[j+2] =  0.95f;
			}

			if (i_spike_old[j]) {
				smooth(0.001, d_gain_old[j], d_gain_old[j+2], release);
			} /*else {
				smooth(0.95, d_gain_old[j], d_gain_old[j+2], attack);
			}*/
			distdata = min(distdata, d_gain_old[j+2]);
			my_gain = distdata / (data + 0.000001f);

			smooth(my_gain, d_gain_old[j+4], d_gain_old[j+6], edge * edge * sampleRate_per441);


			data *= d_gain_old[j+6] * master;

			data = (data > 0.99f) ? 0.99f : data;
			
			//符号を戻す
			(*iip) += sign;
			pData[i] = data;

		}
	}
}



inline static void spike(float gain, float bias, float master, float*pL , float* pR, 
	int size, double* d_gain_old, int* i_spike_old,  float edge, float sampleRate_per441) {

	for (int j = 0; j < 2; j++) 
	{
		float *pData = (j == 0) ? pL : pR;


		for (int i = 0; i < size; i++) {
			float data = pData[i];
			float d_gain;
			

			//符号を取る
			unsigned int* iip = (unsigned int*)&data;
			unsigned int sign = (*iip) & (1 << 31);
	
			//歪の主処理
			d_gain = gain * (1.0f - (sign >> 31) * bias);
			int iVal = d_gain * data;
			int temp = i_spike_old[j];
			i_spike_old[j] = iVal;
			iVal -= temp;

			if (iVal > 0) {
				data = d_gain_old[j] = d_gain_old[j+2] = master;
			} else if (iVal < 0) {
				data = d_gain_old[j] = d_gain_old[j+2] = -master;
			} else {
				data = 0.000001f;
			}

			smooth(data, d_gain_old[j], d_gain_old[j+2], edge * edge * sampleRate_per441);
			data = d_gain_old[j+2] * master;

			pData[i] = data;

		}
	}
}


inline static void water(float gain, float bias, float master, float*pL , float* pR, int size, double* d_gain_old, int* i_spike_old,  float edge, float sampleRate_per441) {

	for (int j = 0; j < 2; j++) 
	{
		float *pData = (j == 0) ? pL : pR;


		for (int i = 0; i < size; i++) {
			float data = pData[i];
			float d_gain;
			

			//符号を取る
			unsigned int* iip = (unsigned int*)&data;
			unsigned int sign = (*iip) & (1 << 31);
			(*iip) -= sign;
	
			d_gain_old[j] = max(d_gain_old[j], data);
			i_spike_old[j+2]++;
			//歪の主処理
			if (i_spike_old[j] != sign) {
				i_spike_old[j] = sign;
				d_gain_old[j+2] += (d_gain_old[j] - d_gain_old[j+2]) * 0.25;
				d_gain_old[j+4] =  i_spike_old[j+2];
				d_gain_old[j] = 0;
				i_spike_old[j+2] = 0;
				if (d_gain_old[j+6] > M_PI * 2) 
				{
					d_gain_old[j+6] -= M_PI * 2;
				}

			}
			smooth(d_gain_old[j+4], d_gain_old[j+8], d_gain_old[j+10], edge * edge);
			d_gain_old[j+6] += M_PI / (d_gain_old[j+10] + 1.0);

			
			data = sin(d_gain_old[j+6] + sin(d_gain_old[j+6]) * bias * 4) * d_gain_old[j+2] * (gain + 1.0f);
			data = (data > 0.99f) ? 0.99f : data;
			
			data *= master;
				

			pData[i] = data;

		}
	}
}




inline static void lp_pin(float gain, float bias, float master, float*pL , float* pR, int size,
	double* d_gain_old, float edge, int* i_spike_old, float sampleRate_per441) {

	float accel_limit = gain;

	//float* current_accel =  ((float*)d_gain_old) + 2;
	

	for (int j = 0; j < 2; j++) 
	{
		float *pData = (j == 0) ? pL : pR;
		double* current_level = (d_gain_old + j);
		double* current_level2 = (d_gain_old + j  + 2);

		for (int i = 0; i < size; i++) {
			float diff = pData[i] - *current_level;
			float diff2 = diff - *current_level;
			//float diff_accel = diff - *current_accel;
			//diff_velocの符号を取る
			uint32_t* iip = (uint32_t*)&diff;
			uint32_t sign = (*iip) & (1 << 31);
			(*iip) -= sign;
			diff = MIN(diff, accel_limit);
			(*iip) += sign;

			*current_level += diff;
			pData[i] = *current_level * master;

		}
	}
}




void Distoo::processReplacing (float** inputs, float** outputs, VstInt32 sampleFrames)
{
	float* in1 = inputs[0];
	float* in2 = inputs[1];
	float* out1 = outputs[0];
	float* out2 = outputs[1];




	while (--sampleFrames >= 0)
	{


		float gain = getRealGain(GET_LPF_PARAM(fGain));
		float master = GET_LPF_PARAM(fMaster);
		float edge = GET_LPF_PARAM(fEdge);
		float bias = GET_LPF_PARAM(fBias);

		//iir[0]->calcIIRCoeff(freq, LPF_Q);
		//iir[1]->calcIIRCoeff(freq, LPF_Q);

		
		float x1 = *in1++;
		float x2 = (in2) ? (*in2++) : 0;

		bufferi[0][cursor] = x1;
		bufferi[1][cursor] = x2;

#ifdef REF

		//アップサンプリングして主処理
		for (int i = 0; i < OVER_SAMPLE ; i++) {
			for (int j = 0; j < 2; j++) { //L,R
				float in = 0;
				//sinc関数の畳み込み
				for (int k = 0; k < SINC_FULL_LEN; k++) { 
					int inOffset = (k + cursor + 1) % SINC_FULL_LEN;
					in += bufferSinc[OVER_SAMPLE - 1 -i + k * OVER_SAMPLE] * bufferi[j][inOffset];
				}
	
				buffero[j][(cursor * OVER_SAMPLE + i)]  = in;
			}
		}

		float* pOutL, *pOutR;
		pOutL = &buffero[0][cursor * OVER_SAMPLE];
		pOutR = &buffero[1][cursor * OVER_SAMPLE];
		//主処理(ディストーション)
		switch (iType) {
			case kDistortion:
				distortion(gain, bias, master, pOutL, pOutR, OVER_SAMPLE, d_gain_old, edge, sampleRate_per441);
				break;
			case kOverdrive:
				overdrive(gain, bias, master, pOutL, pOutR, OVER_SAMPLE, d_gain_old, edge, i_spike_old, sampleRate_per441);
				break;
			case kSpike:
				spike(gain, bias, master, pOutL, pOutR, OVER_SAMPLE, d_gain_old, i_spike_old, edge, sampleRate_per441);
				break;
			case kWater: 
				water(gain, bias, master, pOutL, pOutR, OVER_SAMPLE, d_gain_old, i_spike_old, edge, sampleRate_per441);
				break;
			case kLPPin: 
				lp_pin(gain, bias, master, pOutL, pOutR, OVER_SAMPLE, d_gain_old, edge, i_spike_old,  sampleRate_per441);
				break;
		}

		
		

		//ダウンサンプリング
		float o[2];
		for (int i = 0; i < 2; i++) { //L,R
			o[i] = 0;
			for (int j = 0; j < OVER_SAMPLED_SINC_LEN; j++) {
				int outOffset =((cursor+1) * OVER_SAMPLE + j) % OVER_SAMPLED_SINC_LEN;
				o[i] += bufferSinc[j] *
					buffero[i][outOffset];
			}
		}

		float oL = o[0], oR = o[1];

#endif



#ifdef OPT
				//アップサンプリング
		int shiftSinc = ((SINC_FULL_LEN - cursor) % SINC_FULL_LEN);
		int blockShiftSinc = (shiftSinc & 3) * OVER_SAMPLED_SINC_LEN * 2;
		float* sincIn = &bufferSinc[(shiftSinc & ~3) + blockShiftSinc];
		float* pOutL, *pOutR;
		pOutL = &buffero[0][cursor * OVER_SAMPLE];
		pOutR = &buffero[1][cursor * OVER_SAMPLE];
		float* pInL = bufferi[0];
		float* pInR = bufferi[1];

		
		__m128 oZero =_mm_setzero_ps();

		int sincSpan = SINC_FULL_LEN * 2;
		for (int i = 0; i < OVER_SAMPLE ; i+=2) { //OVER_SAMPLE * 		
			//sinc関数の畳み込み
			__m128 sinc4 = _mm_load_ps(sincIn );
			__m128 l4 = _mm_load_ps(pInL );
			__m128 r4 = _mm_load_ps(pInR );
			__m128 oL4_1 =  _mm_mul_ps(sinc4, l4);
			__m128 oR4_1 = _mm_mul_ps(sinc4, r4);	
			sinc4 = _mm_load_ps(sincIn + sincSpan);
			__m128 oL4_2 = _mm_mul_ps(sinc4, l4);
			__m128 oR4_2 =  _mm_mul_ps(sinc4, r4);
			for (int j = 4; j < SINC_FULL_LEN; j+=4) { 
				__m128 sinc4 = _mm_load_ps(sincIn + j);
				__m128 l4 = _mm_load_ps(pInL + j);
				__m128 r4 = _mm_load_ps(pInR + j);
				oL4_1 = _mm_add_ps(oL4_1, _mm_mul_ps(sinc4, l4));
				oR4_1 = _mm_add_ps(oR4_1, _mm_mul_ps(sinc4, r4));
				sinc4 = _mm_load_ps(sincIn + j + sincSpan);
				oL4_2 = _mm_add_ps(oL4_2, _mm_mul_ps(sinc4, l4));
				oR4_2 = _mm_add_ps(oR4_2, _mm_mul_ps(sinc4, r4));
			}



			oL4_1 = _mm_hadd_ps(oL4_1, oR4_1);
			oL4_2 = _mm_hadd_ps(oL4_2, oR4_2);
			oL4_1 = _mm_hadd_ps(oL4_1, oL4_2);
			__declspec(align(16)) float o[4];
			_mm_store_ps(o , oL4_1);

			pOutL[i  ] = o[0];
			pOutR[i  ] = o[1];
			pOutL[i+1] = o[2];
			pOutR[i+1] = o[3];
			sincIn += sincSpan * 2;
		}
		
		//主処理(ディストーション)
		switch (iType) {
			case kDistortion:
				distortion(gain, bias, master, pOutL, pOutR, OVER_SAMPLE, d_gain_old, edge, sampleRate_per441);
				break;
			case kOverdrive:
				overdrive(gain, bias, master, pOutL, pOutR, OVER_SAMPLE, d_gain_old, edge, i_spike_old, sampleRate_per441);
				break;
			case kSpike:
				spike(gain, bias, master, pOutL, pOutR, OVER_SAMPLE, d_gain_old, i_spike_old, edge, sampleRate_per441);
				break;
			case kWater: 
				water(gain, bias, master, pOutL, pOutR, OVER_SAMPLE, d_gain_old, i_spike_old, edge, sampleRate_per441);
				break;
			case kLPPin: 
				lp_pin(gain, bias, master, pOutL, pOutR, OVER_SAMPLE, d_gain_old, edge, i_spike_old,  sampleRate_per441);
				break;
		}

		


		//ダウンサンプリング
		float* sincOut = &bufferSincOut[shiftSinc * OVER_SAMPLE];
		__m128 oL4 = _mm_setzero_ps();
		__m128 oR4 = _mm_setzero_ps();
		pOutL = buffero[0];
		pOutR = buffero[1];
		for (int j = 0; j < OVER_SAMPLED_SINC_LEN; j+=4) {		
			__m128 sinc4 = _mm_load_ps(sincOut + j);
			__m128 l4 = _mm_load_ps(pOutL + j);
			__m128 r4 = _mm_load_ps(pOutR + j);
			oL4 = _mm_add_ps(oL4, _mm_mul_ps(sinc4, l4));
			oR4 = _mm_add_ps(oR4, _mm_mul_ps(sinc4, r4));
		}

		oL4 = _mm_hadd_ps(oL4, oR4);
		oL4 = _mm_hadd_ps(oL4, oZero);
		__declspec(align(16)) float o[4];
		_mm_store_ps(o , oL4);
		float oL = o[0], oR = o[1];
#endif

		cursor ++;
		if (cursor >= SINC_FULL_LEN) {cursor = 0;}

		*out1++ = oL * invSincArea;
		if (out2)
			*out2++ = oR * invSincArea;


	}
}




