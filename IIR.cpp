#include "iir.h"



#define _USE_MATH_DEFINES
#include <math.h>


void IIR::calcFreqBase(double sampleRate)
{
	freq_base = M_PI * 2.0 / sampleRate;
}

void IIR::calcIIRCoeff(double freq, double q)
{
	double w, alpha, cos_w;

	w = freq * freq_base;
	cos_w = cos(w);
	alpha = sin(w) / q * 2.0;

	b1 = (1.0 - cos_w);
	b0 = b2 = b1 * 0.5;
	double a0 = 1.0 / (1.0 + alpha);
	a1 = -2.0 * cos_w;
	a2 = 1.0 - alpha;

	a1 *= a0;
	a2 *= a0;
	b0 *= a0;
	b1 *= a0;
	b2 *= a0;
}
