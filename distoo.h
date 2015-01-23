//-------------------------------------------------------------------------------------------------------
// VST Plug-Ins SDK
// Version 2.4		$Date: 2006/11/13 09:08:27 $
//
// Category     : VST 2.x SDK Samples
// Filename     : adelay.h
// Created by   : Steinberg Media Technologies
// Description  : Simple Delay plugin (Mono->Stereo)
//
// ｩ 2006, Steinberg Media Technologies, All Rights Reserved
//-------------------------------------------------------------------------------------------------------

#ifndef __adelay__
#define __adelay__

#include "public.sdk/source/vst2.x/audioeffectx.h"

enum
{
	// Global
	kNumPrograms = 9,

	// Parameters Tags
	kType = 0,
	kGain,
	kEdge,
	kBias,
	kMaster,

	kNumParams
};

class Distoo;

//------------------------------------------------------------------------
class ADistoooProgram
{
friend class Distoo;
public:
	ADistoooProgram ();
	~ADistoooProgram () {}

private:	
	float fGain;
	float fBias;
	float fEdge;
	float fMaster;
	float fType;
	char name[24];
};

//------------------------------------------------------------------------

#define _USE_MATH_DEFINES
#include <math.h>



#define SINC_HALF_LEN (16)
#define SINC_FULL_LEN (SINC_HALF_LEN * 2)
#define OVER_SAMPLE (16)
#define OVER_SAMPLED_SINC_LEN (SINC_FULL_LEN * OVER_SAMPLE)
#define LPF_FREQ_TOP	(20000.0)
#define LPF_FREQ_BOTTOM	(200.0)
#define LPF_Q	(0.75 )

#include "iir.h"

class Distoo : public AudioEffectX
{
public:
	Distoo (audioMasterCallback audioMaster);
	~Distoo ();

	//---from AudioEffect-----------------------
	virtual void processReplacing (float** inputs, float** outputs, VstInt32 sampleFrames);

	virtual void setProgram (VstInt32 program);
	virtual void setProgramName (char* name);
	virtual void getProgramName (char* name);
	virtual bool getProgramNameIndexed (VstInt32 category, VstInt32 index, char* text);
	
	virtual void setParameter (VstInt32 index, float value);
	virtual float getParameter (VstInt32 index);
	virtual void getParameterLabel (VstInt32 index, char* label);
	virtual void getParameterDisplay (VstInt32 index, char* text);
	virtual void getParameterName (VstInt32 index, char* text);


	virtual void resume ();

	virtual bool getEffectName (char* name);
	virtual bool getVendorString (char* text);
	virtual bool getProductString (char* text);
	virtual VstInt32 getVendorVersion () { return 1000; }
	
	virtual VstPlugCategory getPlugCategory () { return kPlugCategEffect; }

	virtual void setSampleRate (float sampleRate)  {
		this->sampleRate = sampleRate; 
		sampleRate_per441 =  44100.0f / getSampleRate();
		if (sampleRate_per441 > 1.0f) {sampleRate_per441 = 1.0f;}
	}

	//歪の種類
	enum DistooType {
		kDistortion = 0,
		kOverdrive,
		kSpike,
		kWater,
		kLPPin,
		kNumType
	};




protected:
	void getTypeName(char* text);

#define DECLARE_LPF_PARAM(x) float x; float x##_lpf1; float x##_lpf2;
#define GET_LPF_PARAM(x)  ((x##_lpf1) += ((x) - (x##_lpf1)) * 0.02, (x##_lpf2) += ( (x##_lpf1) - (x##_lpf2)) * 0.02)
#define INIT_LPF_PARAM(x)  x##_lpf1 = x##_lpf2 = x;

	float fGainParam;
	float fBiasParam;
	float fMasterParam;
	float fEdgeParam;
	int iType;

	DECLARE_LPF_PARAM(fGain)
	DECLARE_LPF_PARAM(fBias)
	DECLARE_LPF_PARAM(fMaster)
	DECLARE_LPF_PARAM(fEdge)


	float sampleRate_per441;


	inline float getRealGain(float gain) {
		switch(iType) {
		case kDistortion:
		case kSpike:
			return (gain * gain + 0.000999) * 1000;
		case kWater:
			return gain * 100;
		case kLPPin: 
			return 0.001 * (1.0 - gain * 0.995);
		}
		//case kOverdrive:
		return 1.0 - gain * 0.85;

	}

	float shiftParam(float param, double e) 
	{
		float ret = param * (1.0 - e) + e;
		return (ret > 1.0f) ? 1.0f : ret;
	}





	void resumeIIR();

	void setType(float value);
	void setGain(float value);
	void setBias(float value);
	void setEdge(float value);
	void setMaster(float value);

	float getType();
	float getGain();
	float getBias();
	float getEdge();
	float getMaster();

	//Sinc関数のバッファを生成する
	double makeSincBuf(float* buffer);

	//入力側（アップサンプリング用)のSinc関数のバッファを生成する
	double makeSincBufForInput(float* buffer);
	//出力側（ダウンサンプリング用)のSinc関数のバッファを生成する
	double makeSincBufForOutput(float* buffer);

	ADistoooProgram* programs;
	
	float* bufferi[2]; //LR
	float* buffero[2]; //LR
	float* bufferSinc;
	float* bufferSincIn; //LR
	float* bufferSincOut; //LR
	double invSincArea;

	int cursor;

	//IIR* iir[2];
	double d_gain_old[12];
	int i_spike_old[4];



};





#endif
