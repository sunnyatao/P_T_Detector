#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include "DF_filter.h"
double filter::inputdata(double input) {

	double fil_temp = filter_num[0] * input;			//b(1)*x(n)	
	int temp_ptr = ptr;
	for (int i = 1; i < FIL_ord; i++) {
		fil_temp += filter_num[i] * state[2 * temp_ptr];	//+b(m)x(n-m)
		if (filter_type == 0) {
			fil_temp -= filter_den[i] * state[2 * temp_ptr + 1];
		}
		temp_ptr--;
		if (temp_ptr < 0) {
			temp_ptr += FIL_ord;
		}
	}
	ptr = ptr + 1;												//increase ptr to store new filter state
	if (ptr >= FIL_ord) {
		ptr = 0;
	}
	state[2 * ptr] = input;
	state[2 * ptr + 1] = fil_temp;
	output = fil_temp;
	last_input = input;
	return fil_temp;
}
filter::filter(const double * NUM, const double * DEN, int ORD) {
	filter_num = new double[ORD];
	filter_den = new double[ORD];
	for (int i = 0; i < ORD; i++) {
		filter_num[i] = NUM[i];
		filter_den[i] = DEN[i];
	}
	FIL_ord = ORD;
	state = new double[2 * ORD]();
	ptr = 0;
	output = 0;
	filter_type = 0;
}

filter::filter(const double * NUM, int ORD) {
	filter_num = new double[ORD];
	filter_den = new double[ORD];
	for (int i = 0; i < ORD; i++) {
		filter_num[i] = NUM[i];
	}
	FIL_ord = ORD;
	state = new double[2 * ORD]();
	ptr = 0;
	output = 0;
	filter_type = 1;
}

void filter::reset() {
	ptr = 0;
	for (int i = 0; i < 2 * FIL_ord; i++) {
		state[i] = 0;
	}
}

filter::~filter() {
	delete[] state;
	delete[] filter_num;
	delete[] filter_den;
}
double MA_filter::inputdata(double input) {
	ptr++;
	if (ptr >= MA_ord) {
		ptr = 0;
	}
	output -= buffer[ptr] / MA_ord;
	buffer[ptr] = input;
	output += input / MA_ord;
	last_input = input;
	return output;
}

void MA_filter::reset() {
	ptr = 0;
	for (int i = 0; i < MA_ord; i++) {
		buffer[i] = 0;
	}
}

MA_filter::~MA_filter() {
	delete[] buffer;
}

void findpeak::reset() {
	diff = 0;
	temp = 0;
	MaxFirst = -1;
	currentMax = -1;
	max.clear();
	min.clear();
	max_loc.clear();
	min_loc.clear();
	data.clear();
	peakWidth.clear();
	peakProminence.clear();
	valleyProminence.clear();
	valleyWidth.clear();
	peakArea.clear();
	valleyArea.clear();
}

void findpeak::inputdata(double input) {
	if (input >= temp && diff < 0) {		// /\ a maximal is found
		updateMin(data[data.size() - 1], (long)(data.size() - 1));
		count_st = (long) data.size() - 1;
		if (MaxFirst == -1) {
			MaxFirst = 0;
		}
		currentMax = 0;
	}

	if (diff > 0 && (input <= temp)) {  //  \/  a minimal is found 
		updateMax(data[data.size() - 1], (long)(data.size() - 1));
		count_st = (long)data.size() - 1;
		if (MaxFirst == -1) {
			MaxFirst = 1;
		}
		currentMax = 1;
	}

	if (diff == 0 && (input == temp)) {
		if (currentMax == 1) {
			max_loc[max_loc.size() - 1] = (long)(data.size() - 1 + count_st) / 2; //take middle point to last peak / trough with flat top / bottom
		}
		else {
			if (currentMax == 0) {
				min_loc[min_loc.size() - 1] = (long)(data.size() - 1 + count_st) / 2; //take middle point to last peak / trough with flat top / bottom
			}
		}
	}

	if (diff == 0 && (input > temp)) {
		if (currentMax == 1) {
			max_loc.erase(max_loc.end()-1);
			max.erase(max.end()-1);
		}
	}

	if (diff == 0 && (input < temp)) {
		if (currentMax == 0) {
			min_loc.erase(min_loc.end()-1);
			min.erase(min.end()-1);
		}

	}
	diff = input - temp;
	temp = input;
	data.push_back(input);
}

void findpeak::processPeak(bool WidthByHeight) {
	double tempV,lBV,rBV;
	long tempLoc, lBloc, rBloc;
	for (long i = 0; i < max.size(); i++) {
		//search smallest peak on the left
		long n = i - 1;
		tempV = max[i];
		tempLoc = max_loc[i];
		while (n >= 0) {
			if (max[n] < max[i]) {
				tempV = max[n];
				tempLoc = max_loc[n];
			}
			else {
				break;
			}
			n--;
		}
		n++;
		//keep searching for the valley on the left
		long ii = i;
		if (MaxFirst == 1) {
			n--; // If first peak/trough found is peak, then left trough index == n-1 (assured)
			ii--;
		}
		if (n >= 0) {
			lBV = min[n];
			lBloc = min_loc[n];
			for (int nn = n + 1; nn <= ii; nn++) {
				if (nn >= 0) {
					if (min[nn] < lBV) {
						lBV = min[nn];
						lBloc = min_loc[nn];
					}
				}
			}
		}
		else {
			lBV = data[0];
			lBloc = 0;// If first peak, then left boarder is 0
			for (int nn = n + 1; nn <= ii; nn++) {
				if (nn >= 0) {
					if (min[nn] < lBV) {
						lBV = min[nn];
						lBloc = min_loc[nn];
					}
				}
			}
		}
		//search smallest peak on the right
		n = i + 1;
		tempV = max[i];
		tempLoc = max_loc[i];
		while (n < max.size()) {
			if (max[n] < max[i]) {
				tempV = max[n];
				tempLoc = max_loc[n];	
			}
			else {
				break;
			}
			n++;
		}
		n--;
		
		//keep searching for the valley on the right
		ii = i;
		if (MaxFirst == 0) { //if first element is peak, right valley is the same index or right bound
			n++;
			ii++;
		}
		if (n < min_loc.size()) {
			rBV = min[n];
			rBloc = min_loc[n];
			for (int nn = n - 1; nn >= ii ; nn--) {
				if (nn < min.size()) {
					if (min[nn] < rBV) {
						rBV = min[nn];
						rBloc = min_loc[nn];
					}
				}
			}
		}
		else {
			rBV = data[data.size()-1];
			rBloc = (long)data.size() - 1;
			for (int nn = n - 1; nn >= ii; nn--) {
				if (nn < min.size()) {
					if (min[nn] < rBV) {
						rBV = min[nn];
						rBloc = min_loc[nn];
					}
				}
			}
		}

		float prominence = max[i] - std::max(lBV, rBV);
		float width_reference;
		if (WidthByHeight == false) {
			width_reference = max[i] - 0.5 * prominence;
		}
		else {
			width_reference = 0.5 * max[i];
		}
		double tempLocL = max_loc[i];
		float tempVL = 0;
		float tempVR = 0;
		float tempSlope = 0;
		float temp = 0;
		for (long idx = max_loc[i]; idx >= lBloc; idx--) {
			tempSlope = temp;
			temp = data[idx] - width_reference;
			tempSlope -= temp;
			if (temp <= 0) {
				tempVL = data[idx];
				tempLocL = idx - temp / tempSlope;
				break;
			}
		}

		double tempLocR = max_loc[i];
		temp = 0;
		for (long idx = max_loc[i]; idx <= rBloc; idx++) {					//find half-prominence point on the right
			
			tempSlope = temp;
			temp = data[idx] - width_reference;
			tempSlope -= temp;
			if (temp <= 0) {
				tempVR = data[idx];
				tempLocR = idx + temp / tempSlope;
				break;
			}
		}
		float Vref = std::max(tempVR, tempVL);
		float Vsum = 0;
		for (long idx = tempLocL; idx <= tempLocR; idx++) {
			Vsum += data[idx] - Vref;
		}
		peakArea.push_back(Vsum);
		peakWidth.push_back(tempLocR - tempLocL);
		peakProminence.push_back(prominence);
	}
}

void findpeak::processValley(bool WidthByHeight) {
	double tempV, lBV, rBV;
	long tempLoc, lBloc, rBloc;
	for (long i = 0; i < min.size(); i++) {
		//search largest valley on the left
		long n = i - 1;
		tempV = min[i];
		tempLoc = min_loc[i];
		while (n >= 0) {
			if (min[n] > min[i]) {
				tempV = min[n];
				tempLoc = min_loc[n];
			}
			else {
				break;
			}
			n--;
		}
		n++;
		//keep searching for the peak on the left
		long ii = i;
		if (MaxFirst == 0) {
			n--; // If first peak/trough found is peak, then left trough index == n-1 (assured)
			ii--;
		}
		if (n >= 0) {
			lBV = max[n];
			lBloc = max_loc[n];
			for (int nn = n + 1; nn <= ii; nn++) {
				if (nn >= 0) {
					if (max[nn] > lBV) {
						lBV = max[nn];
						lBloc = max_loc[nn];
					}
				}
			}
		}
		else {
			lBV = data[0];
			lBloc = 0;// If first peak, then left boarder is 0
			for (int nn = n + 1; nn <= ii; nn++) {
				if (nn >= 0) {
					if (max[nn] > lBV) {
						lBV = max[nn];
						lBloc = max_loc[nn];
					}
				}
			}
		}
		//search smallest peak on the right
		n = i + 1;
		tempV = min[i];
		tempLoc = min_loc[i];
		while (n < min.size()) {
			if (min[n] > min[i]) {
				tempV = min[n];
				tempLoc = min_loc[n];
			}
			else {
				break;
			}
			n++;
		}
		n--;

		//keep searching for the valley on the right
		ii = i;
		if (MaxFirst == 1) { //if first element is peak, right valley is the same index or right bound
			n++;
			ii++;
		}
		if (n < max_loc.size()) {
			rBV = max[n];
			rBloc = max_loc[n];
			for (int nn = n - 1; nn >= ii; nn--) {
				if (nn < max.size()) {
					if (max[nn] > rBV) {
						rBV = max[nn];
						rBloc = max_loc[nn];
					}
				}
			}
		}
		else {
			rBV = data[data.size() - 1];
			rBloc = (long)data.size() - 1;
			for (int nn = n - 1; nn >= ii; nn--) {
				if (nn < max.size()) {
					if (max[nn] > rBV) {
						rBV = max[nn];
						rBloc = max_loc[nn];
					}
				}
			}
		}
		float prominence = std::min(lBV, rBV) - min[i];
		float width_reference;
		if(WidthByHeight==false){
			width_reference = min[i] + 0.5 * prominence;
		}
		else {
			width_reference = 0.5 * min[i];
		}
		double tempLocL = min_loc[i];
		float tempVL = 0;
		float tempVR = 0;
		float tempSlope = 0;
		float temp = 0;
		for (long idx = min_loc[i]; idx >= lBloc; idx--) {
			tempSlope = temp;
			temp = data[idx] - width_reference;
			tempSlope -= temp;
			if (temp >= 0) {
				tempVL = data[idx];
				tempLocL = idx - temp / tempSlope;
				break;
			}
		}
		double tempLocR = min_loc[i];
		temp = 0;
		for (long idx = min_loc[i]; idx <= rBloc; idx++) {					//find half-prominence point on the right
			tempSlope = temp;
			temp = data[idx] - width_reference;
			tempSlope -= temp;
			if (temp >= 0) {
				tempVR = data[idx];
				tempLocR = idx + temp / tempSlope;
				break;
			}
		}
		float Vref = std::min(tempVR, tempVL);
		float Vsum = 0;
		for (long idx = tempLocL; idx <= tempLocR; idx++) {
			Vsum += Vref - data[idx];
		}
		valleyArea.push_back(Vsum);
		valleyWidth.push_back(tempLocR - tempLocL);
		valleyProminence.push_back(prominence);
	}
}

void findpeak::swapPeakPosition(long i, long ii) {
	double tempV = max[i];
	long tempLoc = max_loc[i];
	max[i] = max[ii];
	max[ii] = tempV;
	max_loc[i] = max_loc[ii];
	max_loc[ii] = tempLoc;
	if (peakProminence.size() > 0) {
		long tempW = peakWidth[i];
		double tempA = peakArea[i];
		double tempP = peakProminence[i];
		peakProminence[i] = peakProminence[ii];
		peakProminence[ii] = tempP;
		peakWidth[i] = peakWidth[ii];
		peakWidth[ii] = tempW;
		peakArea[i] = peakArea[ii];
		peakArea[ii] = tempA;
	}
}

void findpeak::swapValleyPosition(long i, long ii) {
	double tempV = min[i];
	long tempLoc = min_loc[i];
	min[i] = min[ii];
	min[ii] = tempV;
	min_loc[i] = min_loc[ii];
	min_loc[ii] = tempLoc;
	if (valleyProminence.size() > 0) {
		long tempW = valleyWidth[i];
		double tempA = valleyArea[i];
		double tempP = valleyProminence[i];
		valleyProminence[i] = valleyProminence[ii];
		valleyProminence[ii] = tempP;
		valleyWidth[i] = valleyWidth[ii];
		valleyWidth[ii] = tempW;
		valleyArea[i] = valleyArea[ii];
		valleyArea[ii] = tempA;
	}
}

void findpeak::dropPeak(long i) {
	max_loc.erase(max_loc.begin() + i);
	max.erase(max.begin() + i);
	if (peakProminence.size() > 0) {
		peakProminence.erase(peakProminence.begin() + i);
		peakWidth.erase(peakWidth.begin() + i);
		peakArea.erase(peakArea.begin() + i);
	}
}

void findpeak::dropValley(long i) {
	min_loc.erase(min_loc.begin() + i);
	min.erase(min.begin() + i);
	if (valleyProminence.size() > 0) {
		valleyProminence.erase(valleyProminence.begin() + i);
		valleyWidth.erase(valleyWidth.begin() + i);
		valleyArea.erase(valleyArea.begin() + i);
	}
}

void findpeak::sortByProminence() {
	double tempV, tempP;
	long tempLoc, tempW;
	//sort peaks by prominence
	if (peakProminence.size() > 0) {
		for (int i = 0; i < max.size(); i++) {
			for (int ii = i + 1; ii < max.size(); ii++) {
				if (peakProminence[ii] > peakProminence[i]) {
					swapPeakPosition(i, ii);
				}
			}
		}
	}
	if (valleyProminence.size() > 0) {
		for (int i = 0; i < min.size(); i++) {
			for (int ii = i + 1; ii < min.size(); ii++) {
				if (valleyProminence[ii] < valleyProminence[i]) {
					swapValleyPosition(i, ii);
				}
			}
		}
	}
}

void findpeak::sortByValue() {
	double tempV, tempP;
	long tempLoc, tempW;
	//sort peaks by value
	for (int i = 0; i < max.size(); i++) {
		for (int ii = i + 1; ii < max.size(); ii++) {
			if (max[ii] > max[i]) {
				swapPeakPosition(i, ii);
			}
		}
	}
}

void findpeak::dropByProminence(double prominence) {
	for (long i = peakProminence.size() - 1; i >= 0; i--) {
		if (peakProminence[i] < prominence) {
			dropPeak(i);
		}
	}
	for (long i = valleyProminence.size() - 1; i >= 0; i--) {
		if (valleyProminence[i] < prominence) {
			dropValley(i);
		}
	}
}

void findpeak::dropByHeight(double height) {
	for (long i = max.size() - 1; i >= 0; i--) {
		if (max[i] < height) {
			dropPeak(i);
		}
	}
	for (long i = min.size() - 1; i >= 0; i--) {
		if (min[i] > -height) {
			dropValley(i);
		}
	}
}

void findpeak::dropByDistance(long distance) {
	for (long i = max.size() - 1; i >= 0; i--) {
		for (long ii = (long)max.size() - 1; i >= 0; i--) {
			if (abs(max_loc[i] - max_loc[ii]) < distance) {
				dropPeak(ii);
			}
		}
	}
	for (long i = min.size() - 1; i >= 0; i--) {
		for (long ii = (long)min.size() - 1; i >= 0; i--) {
			if (abs(min_loc[i] - min_loc[ii]) < distance) {
				dropValley(ii);
			}
		}
	}
}

void findpeak::updateMax(double input, long cnt) {
	max.push_back(input);
	max_loc.push_back(cnt);
}

void findpeak::updateMin(double input, long cnt) {
	min.push_back(input);
	min_loc.push_back(cnt);
}

findpeak::~findpeak() {
	max.clear();
	min.clear();
	max_loc.clear();
	min_loc.clear();
	data.clear();
	peakProminence.clear();
	peakWidth.clear();
	valleyProminence.clear();
	valleyWidth.clear();
	peakArea.clear();
	valleyArea.clear();
}

float findpeak::SNR() {
	float SNRvalue;
	if (peakProminence.size() == 0) {
		processPeak(false);
		sortByProminence();
	}
	if (peakProminence.size() <=1) {
		SNRvalue = 1;
	}
	if (peakProminence.size() >=2) {
		float sum = 0;
		for (long i = 0; i < peakProminence.size(); i++) {
			sum = sum + peakProminence[i];
		}
		SNRvalue = (peakProminence[0] + peakProminence[1]) / sum;
		
	}
	return SNRvalue;
}

float findpeak::peakSNR(int N) {
	float sum1 = 0;
	float sum2 = 0;
	if (peakProminence.size() == 0) {
		processPeak(false);
		sortByProminence();
	}
	for (long i = 0; i < peakProminence.size(); i++) {
		if (i <= N) {
			sum1 += peakProminence[i];
		}
		sum2 += peakProminence[i];
	}
	if (peakProminence.size() == 0)
	{
		return 0;
	}
	else
	{
		return sum1 / sum2;
	}
}

float findpeak::valleySNR(int N) {
	float sum1 = 0;
	float sum2 = 0;
	if (valleyProminence.size() == 0) {
		processValley(false);
		sortByProminence();
	}
	for (long i = 0; i < valleyProminence.size(); i++) {
		if (i <= N) {
			sum1 += valleyProminence[i];
		}
		sum2 += valleyProminence[i];
	}
	if (valleyProminence.size() == 0)
	{
		return 0;
	}
	else
	{
		return sum1 / sum2;
	}
}

void signalProps::inputdata(double input) {
	if (input > max) {
		max = input;
		max_loc = index;
	}
	if (input < min) {
		min = input;
		min_loc = index;
	}
	index++;
	sum += input;
	mean = sum / index;
}
void signalProps::reset() {
	index = 0;
	max = 0;
	min = 0;
	max_loc = 0;
	min_loc = 0;
	sum = 0;
}