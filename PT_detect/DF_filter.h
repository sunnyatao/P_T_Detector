#pragma once
#include <vector>
using namespace std;
class filter{
public:
	double * filter_num;
	double * filter_den;
	int FIL_ord;
	double* state;
	int ptr;
	int filter_type; //0=IIR 1=FIR
	double output, last_input;
	double inputdata(double input);
	void reset();
	filter(const double * NUM, const double * DEN, int ORD);
	filter(const double * NUM, int ORD);
	filter(double * NUM, double * DEN, int ORD)	{
		filter_num = NUM;
		filter_den = DEN;
		FIL_ord = ORD;
		state = new double[2 * ORD]();
		ptr = 0;
		output = 0;
		filter_type = 1;
	}
	~filter();
};


class MA_filter {
public:
	int MA_ord;
	double* buffer;
	int ptr;
	double output, last_input;
	void reset();
	double inputdata(double input);
	MA_filter(int ORD) {
		MA_ord = ORD;
		buffer = new double[MA_ord]();
		output = 0;
		ptr = 0;
	}
	~MA_filter();
};

class findpeak {
public:
	vector <double> max, min, data, peakProminence, valleyProminence, peakArea, valleyArea, peakWidth, valleyWidth;;
	vector <long> max_loc, min_loc;
	double  temp, diff;
	long  count_st;
	int Npeaks;
	int MaxFirst = -1, currentMax = -1;
	void reset();
	void inputdata(double input);
	void processPeak(bool WidthByHeight);
	void processValley(bool WidthByHeight);
	void sortByProminence();
	void sortByValue();
	void updateMin(double input,long cnt);
	void updateMax(double input,long cnt);
	void dropByHeight(double height);
	void dropByDistance(long distance);
	void dropByProminence(double prominence);
	void swapPeakPosition(long a, long b);
	void swapValleyPosition(long a, long b);
	void dropPeak(long idx);
	void dropValley(long idx);
	float SNR();
	float peakSNR(int N);
	float valleySNR(int N);
	findpeak() {
		diff = 0;
	}
	~findpeak();
};

class signalProps {
public:
	double max, min ,sum ,mean;
	long max_loc, min_loc;
	long index;
	void reset();
	void inputdata(double input);
	signalProps() {
		index = 0;
		max = 0;
		min = 0;
		max_loc = 0;
		min_loc = 0;
		sum = 0;
	}
};


