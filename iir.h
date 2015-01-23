


class IIR {

public:


	void calcIIRCoeff(double freq, double q);

	

	IIR(double freq, double q, double sampleRate) {
		calcFreqBase(sampleRate);
		calcIIRCoeff(freq, q);
	}



	__forceinline double DoIIRMono(double input)
	{
		double yn =	
			    b0 * input 
			+	b1 * in1 
			+	b2 * in2 
			-	a1 * out1 
			-	a2 * out2;

		in2=in1;
		in1=input;
		out2=out1;
		out1=yn;

		return yn;

	}

	void calcFreqBase(double sampleRate);

	void resume () {
		in1 = in2 = out1 = out2 = 0;
	}



protected:
	double freq_base;
	double a1, a2, b0, b1, b2;
	double in1, in2, out1, out2;


};