#include <string.h>
#include <stdio.h>
#include "PT_detect.h"


void QRS_earse(float *Data, int dataLen, const int32_T* r_pos_arr, const INT16_T* start_end_arr, int r_pos_len)
{

	for (int i = 0; i < r_pos_len; i++) {
		int r_pos = r_pos_arr[i];
		unsigned int startInd = r_pos - start_end_arr[2 * i];
		unsigned int endInd = r_pos + start_end_arr[2 * i + 1];
		float startVal = Data[startInd];
		float endVal = Data[endInd];


		int diffNum = endInd - startInd;
		float slope = (endVal - startVal) / (float)diffNum;
		for (int ind = 0; ind < diffNum; ind++)
		{
			float oldVal = Data[startInd + ind];
			Data[startInd + ind] = startVal + ind*slope;
			//printf("ind:%d--old val:%d new val:%d\n", startInd + ind, oldVal, Data[startInd + ind]);
		}
	}
}

void PT_longterm(const float* signal, const int32_T signalsSize, const int32_T* Rposes, const INT16_T* qrs_start_ends, const bool* Noisetag, const int32_T RposesSize, const int fs, int32_T *Tposes, int32_T * Tnum, int32_T *Pposes, int32_T *Pnum, float * rr_noise_levels) {

	int32_T start_idx = 0;
	int32_T end_idx = 0;
	const int RR_count = MAX_int32_T; //�����Ƭ��RR�ĸ���
	int32_T Max_interval = 15 * fs; //������ʱ�䲻���PT
	*Tnum = 0; //���㿪ʼ��¼
	*Pnum = 0;
	int32_T R_index = 0;

	float* cp_signal = new float[signalsSize];
	memcpy(cp_signal, signal, signalsSize * sizeof(float));
	QRS_earse(cp_signal, signalsSize, Rposes, qrs_start_ends, RposesSize);


	//filter first
	filter *fil_HP, *fil_HL;
	int fil_HP_delay;
	int fil_HL_delay = (BL_HL - 1) / 2;
	double temp = 0;
	float *fil_HP_buff = new float[signalsSize]();
	float *fil_HL_buff = new float[signalsSize]();
	fil_HL = new filter(B_HL, BL_HL);
	if (fs == 256) {
		fil_HP = new filter(B_BP, BL_BP);							//create a filter class to process data from a single ECG channel
		fil_HP_delay = (BL_BP - 1) / 2;
	}
	else {
		fil_HP = new filter(B_BP_250, BL_BP_250);					//create a filter class to process data from a single ECG channel
		fil_HP_delay = (BL_BP_250 - 1) / 2;
	}
	//fil_HP = new filter(B_W18, BL_W18);
	//fil_HP_delay = ( BL_W18 - 1 ) / 2;
	temp = 0;
	for (int i = 0; i < signalsSize; i++) {
		fil_HP->inputdata(cp_signal[i]);
		float HP_V = (float)fil_HP->output - temp;
		//temp = fil_HP->output;//differential
		fil_HL->inputdata(HP_V);
		if (i - fil_HP_delay >= 0) {
			fil_HP_buff[i - fil_HP_delay] = HP_V;                 //����Ⱥʱ��
		}
		if (i - fil_HP_delay - fil_HL_delay >= 0) {
			fil_HL_buff[i - fil_HL_delay - fil_HP_delay] = fil_HL->output;			//hilbert transform
		}
	}
	for (int i = 0; i < signalsSize; i++) {
		fil_HL_buff[i] = sqrt(fil_HL_buff[i] * fil_HL_buff[i] + fil_HP_buff[i] * fil_HP_buff[i]);  //�����
	}
	//float meanValue = mean(fil_HL_buff, signalsSize);
	//for (int i = 0; i < signalsSize; i++) {
	//	fil_HL_buff[i] = fil_HL_buff[i] - meanValue;										//ȥֱ���ź�
	//}

	//��Ƶ������ʼ��
	for (int numSNR = 0; numSNR < RposesSize;numSNR++) {
		rr_noise_levels[numSNR] = (float)1;
	}

	while (start_idx < RposesSize) {
		while (start_idx < RposesSize) {
			//�ҷ��������
			if (!Noisetag[start_idx]) {
				break;
			}
			else {
				start_idx++;
			}
		}
		end_idx = start_idx;
		int32_T end_last=start_idx;
		while (end_idx < RposesSize - 1 && ((end_idx - start_idx) < RR_count) ) {
			end_last = end_idx;
			end_idx++;

			if ((Rposes[end_idx] - Rposes[end_last]) > Max_interval) {
				end_idx = end_last; //move to last not a noise end
				break;
			}

			//�������յ�
			if (Noisetag[end_idx]) {
				end_idx = end_last; //move to last not a noise end
				break;
			}
			
		}
		int32_T r_num = end_idx - start_idx +1;
		if (r_num > 1) {
			int32_T* Rposes_this = new int32_T[r_num +10];
			memcpy(Rposes_this, Rposes + start_idx, r_num * sizeof(int32_T));
			PT_detect(fil_HL_buff, signalsSize, Rposes_this, r_num, fs, Tposes, Tnum, Pposes, Pnum, rr_noise_levels, start_idx);  //iterate sessions
			delete[] Rposes_this;
		}
		start_idx = end_idx + 1; //start next cycle from n+1
	}


	//for (int i = 0; i < *Pnum; i++)
	//{
	//	if(p_pos_arr[i] > 13945000 && p_pos_arr[i] < 13947000)
	//	printf("%d\n", Pposes[i]);
	//}

	delete[] cp_signal;
	delete[] fil_HP_buff;
	delete[] fil_HL_buff;
	delete fil_HP;
	delete fil_HL;
}

void PT_detect(const float* signal, const int32_T signalsSize, const int32_T* Rposes, const int32_T RposesSize, const int32_T fs, int32_T *Tposes, int32_T * Tnum, int32_T *Pposes, int32_T *Pnum, float *rr_noise_levels, int32_T start_idx) {
	int32_T current_Tnum[1];
	*current_Tnum = *Tnum;//��¼���μ�¼T��offset����Ppeak����
	findTpeaks(signal, signalsSize, Rposes, RposesSize, fs, Tposes, Tnum);
	findPpeaks(signal, signalsSize, Rposes, RposesSize, fs, Tposes, current_Tnum, Pposes, Pnum,rr_noise_levels, start_idx);
}


void findTpeaks(const float* signal, const  int32_T signalsSize, const  int32_T* Rposes, const int32_T RposesSize, const int fs, int32_T *Tposes, int32_T * Tnum) {
	int k1 = 0;
	int k2 = 0;
	int N = 0;
	int32_T stop_index;
	int32_T start_index;
	signalProps *sp1 = new signalProps;
	int32_T B_index;
	float tempSlpTa;
	float tempSlpTc;
	double *SlpT1;
	double *SlpT2;
	double SlpT1_treh;//T������ֵ
	double SlpT2_treh;//T������ֵ
	int32_T *tempT;
	int numT;
	findpeak * peakfinder = new findpeak;
	int tmp_buf_len = Rposes[RposesSize-1] - Rposes[0];

	int32_T *RR = (int32_T*)malloc((RposesSize) * sizeof(int32_T));
	SlpT1 = (double*)malloc(tmp_buf_len * sizeof(double));	//��T���߶�
	SlpT2 = (double*)malloc(tmp_buf_len * sizeof(double));	//��T���߶�
	tempT = (int32_T*)malloc(tmp_buf_len * sizeof(int32_T));

	

	//��һ������ݽ���T����ֵѵ��
	for (int i = 0; i < RposesSize-1;i++) {
		RR[i] = Rposes[i + 1] - Rposes[i];
	}
	int32_T mean_RR = mean(RR, RposesSize-1);
	if (RposesSize == 2 || RposesSize == 3)
	{
		N = RposesSize - 1;
	}
	else
	{
		N = RposesSize / 2;
	}
	if (N > 10) { N = 10; }									//          <--------------------------------------------�Գ��������Ƿ���Ҫ�޸�
	k1 = 0;
	k2 = 0;

	for (int j = 0; j < N - 1; j++) {
		//RR[j] = Rposes[j + 1] - Rposes[j];
		stop_index = floor(Rposes[j] +0.5*fs);			//�趨ÿ�ε�RR���߽�
		start_index = floor(Rposes[j] + 30);
		if (start_index > signalsSize) {
			start_index = signalsSize;
		}
		if (stop_index < 0) {
			stop_index = 0;
		}
		peakfinder->reset();
		for (int i = start_index; i < stop_index; i++) {
			peakfinder->inputdata(signal[i]);				//��ÿ�ε�RR֮��Ѱ�����ֵ
		}
		//peakfinder->processdata();//Ѱ�Ҳ���, ��Ϊ���ɵ�T��

		for (int i = 0; i < peakfinder->max.size(); i++)
		{
			B_index = peakfinder->max_loc[i] + start_index;
			tempSlpTa = signal[B_index] - signal[B_index + 5];			//ͨ����ֵ����12��λ�ü���б��
			tempSlpTc = signal[B_index] - signal[B_index - 5];
			if (tempSlpTa > 0 && tempSlpTc > 0)							//�ֲ����ֵ
			{
				if (k1 > tmp_buf_len) 
				{
					printf("k1:%d len:%d\n", k1, tmp_buf_len);
				}
				SlpT1[k1] = (tempSlpTa + tempSlpTc) / 2;				//���ÿ�����ƽ��б��
				k1++;
			}
		}
		//for (int i = 0; i < peakfinder->min.size(); i++)
		//{
		//	B_index = peakfinder->min_loc[i] + start_index;
		//	tempSlpTa = signal[B_index] - signal[B_index + mean_RR / 20];
		//	tempSlpTc = signal[B_index] - signal[B_index + mean_RR / 20];
		//	if (tempSlpTa < 0 && tempSlpTc < 0) {		
		//		if (k2 > tmp_buf_len)
		//		{
		//			printf("k2:%d len:%d\n", k2, tmp_buf_len);
		//		}
		//		//�ֲ���Сֵ
		//		SlpT2[k2] = (tempSlpTa + tempSlpTc) / 2;
		//		k2++;
		//	}
		//}
	}
	//ѵ���õ�T��������ֵ(0.2����б�ʣ�������ֵ�ֱ����)
	SlpT1_treh = 0.05*maxValue(SlpT1, k1);
	//SlpT2_treh = 0.05*maxValue(SlpT2, k2);
	//����ֵ����T�����
	
	for (int j = 0; j < RposesSize - 1; j++)
	{
		//if (Rposes[j] >= 1107150) {
		//	printf("====");
		//}
		//RR[j] = Rposes[j + 1] - Rposes[j];
		stop_index = floor(Rposes[j] +0.5*fs);
		start_index = floor(Rposes[j] + 30);
		if (start_index > signalsSize) {
			start_index = signalsSize;
		}
		if (stop_index < 0) {
			stop_index = 0;
		}
		peakfinder->reset();
		for (int i = start_index; i < stop_index; i++) {
			peakfinder->inputdata(signal[i]);//Ѱ�������ڵĲ��壬��Ϊ���ɵ���T��
		}
		numT = 0;
		for (int i = 0; i < peakfinder->max.size(); i++)
		{
			B_index = peakfinder->max_loc[i] + start_index;
			tempSlpTa = signal[B_index] - signal[B_index +5];				//ͨ����ֵ����12��λ�ü���б��
			tempSlpTc = signal[B_index] - signal[B_index -5];
			/*if (tempSlpTa > 0 && tempSlpTc > 0) {	*/						//�ֲ����ֵ
				//if ((tempSlpTa + tempSlpTc) / 2 > SlpT1_treh)				//�ֲ����ֵ��ƽ���߶ȴ�����ֵ(����prominence���)
				//{
					if (numT > tmp_buf_len)
					{
						printf("1111numT:%d len:%d\n", numT, tmp_buf_len);
					}
					tempT[numT] = B_index;									//�洢���п�����T��
					numT++;
				//}
		/*	}*/
		}
		//Ѱ�Ҳ��ȣ���Ϊ���ɵĸ�T��
		//for (int i = 0; i < peakfinder->min.size(); i++)
		//{
		//	B_index = peakfinder->min_loc[i] + start_index;
		//	tempSlpTa = signal[B_index] - signal[B_index + mean_RR / 20];
		//	tempSlpTc = signal[B_index] - signal[B_index - mean_RR / 20];
		//	if (tempSlpTa < 0 && tempSlpTc < 0)// ��ѵ����T��б����ֵ������֤
		//	{
		//		if (fabs(tempSlpTa + tempSlpTc) / 2 > fabs(SlpT2_treh))
		//		{
		//			if (numT > tmp_buf_len)
		//			{
		//				printf("222numT:%d len:%d\n", numT, tmp_buf_len);
		//			}
		//			tempT[numT] = B_index;									//�洢���п��ܸ�T��
		//			numT++;
		//		}
		//	}
		//}
		//�ߵ�һ��RR�����ڶ�����T�㣬��T���ԭʼ�ź�ֵ������֤���������T��
		//if (numT < 1) {
		//	continue;
		//}
		sp1->reset();
		for (long t = 0; t < numT; t++) {
			sp1->inputdata(signal[tempT[t]]);
		}
		for (int t = 0; t < numT;t++) {
			if (signal[tempT[t]] >= 0.7*signal[tempT[sp1->max_loc]] && tempT[t] < tempT[sp1->max_loc]) {
				Tposes[j + *Tnum] = tempT[t];			
			}
			else {
			Tposes[j + *Tnum] = tempT[sp1->max_loc];
			}

		}
		

		//int i = 0;
		///*numT = n1+n2;*/
		//int z = 0;
		//int t = 0;															// <-----------------------------�����������û���ã�
		//
		//while (numT > 1) {
		//	if (fabs(signal[tempT[0]]) <= fabs(signal[tempT[1]])) {		// ���t���źŷ���С�ڵ�����һ��t���ֵ(��һ��T)
		//		for (z = 0; z < numT - 1; z++)//�ߵ���ǰ��					//�����t�㿪ʼ
		//		{
		//			tempT[z] = tempT[z + 1];								//������������ǰshift 1
		//		}
		//		numT--;														//��������size
		//	}
		//	else {
		//		for (z = 1; z < numT - 1; z++)//�ߵ���һ����			//���t���źŷ��ȴ�
		//		{
		//			tempT[z] = tempT[z + 1];								//���������ݴ�t+1��ǰshift 1
		//		}
		//		numT--;														//��С����size
		//	}
		//}
		//Tposes[j + *Tnum] = tempT[0];  //�ڶ�β��������
	}
	*Tnum += RposesSize - 1;			  //���г���+1
	free(SlpT1);
	free(SlpT2);
	free(tempT);
	free(RR);
	delete peakfinder;
	delete sp1;
}

void findPpeaks(const float* signal, const  int32_T signalsSize, const  int32_T* Rposes, const int32_T RposesSize, const int fs, int32_T *Tposes, int32_T * Tnum, int32_T *Pposes, int32_T * Pnum, float *rr_noise_levels, int32_T start_idx) {

	int k1 = 0;
	int k2 = 0;
	int N = 0;

	int32_T stop_index;
	int32_T start_index;
	int32_T B_index;
	float tempSlpPa;
	float tempSlpPc;
	double SlpP1_treh,signaldiff;
	int numP;
	findpeak * peakfinder = new findpeak;
	int32_T PRdiff, PR1diff, PRvalue;

	int tmp_buf_len = Rposes[RposesSize - 1] - Rposes[0];

	int32_T *RR = (int32_T*)malloc((RposesSize) * sizeof(int32_T));
	int32_T *PR = (int32_T*)malloc((RposesSize) * sizeof(int32_T));
	int32_T *tempP = (int32_T*)malloc(tmp_buf_len * sizeof(int32_T));
	float *SlpP1 = (float*)malloc(tmp_buf_len * sizeof(float));   //��P���߶�
	
	//��һ������ݽ���T����ֵѵ��
	//FILE *tfile = NULL;
	//tfile = fopen("E:\\matlabWorkspace\\testdata\\�����\\filter_buf", "wb+");
	//fwrite(signal, sizeof(float),signalsSize, tfile);
	//fclose(tfile);

	for (int i = 0; i < RposesSize - 1; i++) {
		RR[i] = Rposes[i + 1] - Rposes[i];
	}
	int32_T mean_RR = mean(RR, RposesSize - 1);
	if (RposesSize == 2 || RposesSize == 3)
	{
		N = RposesSize - 1;
	}
	else
	{
		N = RposesSize / 2;
	if (N > 20) 
	{ 
		N = 20; 
	}
	}

	k1 = 0;
	
	for (int j = 0; j < N - 1; j++) {
		//RR[j] = Rposes[j + 1] - Rposes[j];
		stop_index = floor(Rposes[j + 1] - mean_RR /15);
		start_index =floor( Rposes[j] + 0.1*fs);
		if (start_index > signalsSize) {
			start_index = signalsSize;
		}
		if (stop_index < 0) {
			stop_index = 0;
		}
		peakfinder->reset();
		for (int i = start_index; i < stop_index; i++) {
			peakfinder->inputdata(signal[i]);//Ѱ�������ڵĲ��壬��Ϊ���ɵ���T��
		}

		for (int i = 0; i < peakfinder->max.size(); i++)
		{
			B_index = peakfinder->max_loc[i] + start_index;
			if (B_index <= Tposes[j + *Tnum] + mean_RR / 10) {				//�������P�������T��֮ǰ��ɾ��
				continue;
			}
			//if (abs(signal[B_index]) > 0.7*abs(signal[Rposes[j]])) {
			//	continue;
			//}
			tempSlpPa = signal[B_index] - signal[B_index + 5];			//ͨ����ֵ����15��λ�ü���б��
			tempSlpPc = signal[B_index] - signal[B_index - 5];
			if (tempSlpPa > 0 && tempSlpPc > 0)
			{
				SlpP1[k1] = (tempSlpPa + tempSlpPc) / 2;
				k1++;
			}
		}
	}
	//ѵ���õ�P����ֵ
	SlpP1_treh = 0.05*maxValue(SlpP1, k1);
	//����ֵ����P�����
	int k = 0;

	for (int j = 0; j < RposesSize - 1; j++)
	{
		//if(Rposes[j]>= 2764730){
		//	printf("=");
		//}
		//RR[j] = Rposes[j + 1] - Rposes[j];
		start_index = floor(Rposes[j] + 0.1*fs);						//Ѱ�������ڵĲ��壬��Ϊ���ɵ�P��
		stop_index = floor(Rposes[j + 1] - mean_RR / 15);
		if (start_index > signalsSize) {
			start_index = signalsSize;
		}
		if (stop_index < 0) {
			stop_index = 0;
		}

		peakfinder->reset();
		for (int i = start_index; i < stop_index; i++) {
			peakfinder->inputdata(signal[i]);							//Ѱ�������ڵĲ��壬��Ϊ���ɵ���P��
		}
		numP = 0;
		
		for (int i = 1; i < peakfinder->max.size(); i++)				//�ӵڶ�����ֵ��ʼ��P
		{
			B_index = peakfinder->max_loc[i] + start_index;
			if (B_index <= Tposes[j + *Tnum]) {				//�������P�������T��֮ǰ��ɾ��(�Թ�����ֵ)
				continue;
			}
			//if (abs(signal[B_index]) > 0.7*abs(signal[Rposes[j]])) {
			//	continue;
			//}
			tempSlpPa = signal[B_index] - signal[B_index + 5];
			tempSlpPc = signal[B_index] - signal[B_index - 5];
			if (tempSlpPa > 0 && tempSlpPc > 0) {
				//if ((tempSlpPa + tempSlpPc) / 2 > SlpP1_treh)
				//{
					tempP[numP] = B_index;
					numP++;
				//}
			}
		}

		//  �ߵ���ԣ���������������P�㣬�Ǵ�������ֻ��һ��P��
		int i = 0;
		int z = 0;
		int t = 0;

		//�����һ��RR���ڲ��Ǵ������ͣ��÷�ֵ��Сѡ��p�� 
		if (numP < 1) {
			continue;
		}
		if (j < 2) {													//����ǰ����RR���ڵ�P��
			while (numP > 1) {
				if (signal[tempP[0]] <= signal[tempP[1]]) {
					for (z = 0; z < numP - 1; z++)//�ߵ���ǰ��
					{
						tempP[z] = tempP[z + 1];
					}
					numP--;
				}
				else {
					for (z = 1; z < numP - 1; z++)//�ߵ���һ����
					{
						tempP[z] = tempP[z + 1];
					}
					numP--;
				}
			}
			Pposes[k + *Pnum] = tempP[numP - 1];
			k++;
			peakfinder->processPeak(false);
			peakfinder->sortByProminence();


			rr_noise_levels[j + start_idx] = 1-peakfinder->peakSNR(2);
		}

		else {	//�ӵ�����RR���ڿ�ʼ�жϴ�������
			//int32_T mean_temp = mean(RR, j + 1);
			if ((RR[j] >= 1.2 * mean_RR|| RR[j]>1.2*fs) && numP > 1) {				//���RR�����Ǵ������ͣ�������P������÷�ֵ��С��RR������ѡ��P��
				while (t < numP - 1) {
					if (abs((int32_T)(tempP[t + 1] - tempP[t])) < 0.4 * mean_RR) {
						if (fabs(signal[tempP[t]]) <= fabs(signal[tempP[t + 1]])) {
							for (z = t; z < numP - 1; z++)//�ߵ���ǰ��
							{
								tempP[z] = tempP[z + 1];
							}
							numP = numP - 1;
						}
						else {
							for (z = t + 1; z < numP - 1; z++)//�ߵ���һ����
							{
								tempP[z] = tempP[z + 1];
							}
							numP = numP - 1;
						}
					}
					else {
						t++;
					}

				}
				//�����ѡP��ľ�ֵ
					float meanTempP = 0;			
					float sum = 0;
					int y = 0;
				for (int x = 0; x < numP;x++) {
					sum = signal[tempP[x]] + sum;					
				}
				meanTempP = sum / numP;
				for (int x = 0; x < numP; x++)
				{
					if (signal[tempP[x]]<0.7*meanTempP) {
						continue;										//�ų���RR�����ڶ�����ֵ�ܵ͵Ĳ���
					}
					Pposes[k + *Pnum] = tempP[x];
					k++;
					y++;
				}
				//�����Ƶ���������
				peakfinder->processPeak(false);
				peakfinder->sortByProminence();
				rr_noise_levels[j + start_idx] = 1-peakfinder->peakSNR(y+1);
			}

			//���RR���ڲ��Ǵ������ͣ�ֻ��һ��P��������PR�������ֵ��ѡ��P��λ��
			//else {
			//	//for (int i = 0; i < RposesSize - 1; i++) {
			//	//	PR[i] = abs((int32_T)(Pposes[k - 1 + *Pnum] - Rposes[i]));
			//	//}
			///*	PRvalue = minValue(PR, RposesSize - 1);*/
			//	while (numP > 1) {
			//		//PRdiff = abs(abs((int32_T)(tempP[t] - Rposes[j + 1])) - PRvalue);
			//		//PR1diff = abs(abs((int32_T)(tempP[t + 1] - Rposes[j + 1])) - PRvalue);
			//		signaldiff = signal[tempP[t]] - signal[tempP[t + 1]];
			//		if (signaldiff <= 0) {
			//			if (PRdiff > PR1diff) {
			//				for (z = t; z < numP - 1; z++)//�ߵ���ǰ��
			//				{
			//					tempP[z] = tempP[z + 1];
			//				}
			//				numP--;
			//			}
			//			else {
			//				if (abs(signaldiff)*10000 > abs(PRdiff - PR1diff)) {
			//					for (z = t; z < numP - 1; z++)//�ߵ���ǰ��
			//					{
			//						tempP[z] = tempP[z + 1];
			//					}
			//					numP--;
			//				}
			//				else {
			//					for (z = t + 1; z < numP - 1; z++)//�ߵ���һ����
			//					{
			//						tempP[z] = tempP[z + 1];
			//					}
			//					numP--;
			//				}
			//			}
			//		}
			//		else {
			//			if (PRdiff <= PR1diff) {
			//				for (z = t + 1; z < numP - 1; z++)//�ߵ���һ����
			//				{
			//					tempP[z] = tempP[z + 1];
			//				}
			//				numP--;
			//			}
			//			else {
			//				if (abs(signaldiff)*10000 >= abs(PRdiff - PR1diff)) {
			//					for (z = t + 1; z < numP - 1; z++)//�ߵ���һ����
			//					{
			//						tempP[z] = tempP[z + 1];
			//					}
			//					numP--;
			//				}
			//				else {
			//					for (z = t; z < numP - 1; z++)//�ߵ���ǰ��
			//					{
			//						tempP[z] = tempP[z + 1];
			//					}
			//					numP--;
			//				}
			//			}
			//		}
			//	}
			//	if (numP > 0) {
			//		Pposes[k + *Pnum] = tempP[0];
			//		k++;
			//	}
			//
			//			while (numP > 1) {
			//	if (fabs(signal[tempP[0]]) <= fabs(signal[tempP[1]])) {
			//		for (z = 0; z < numP - 1; z++)//�ߵ���ǰ��
			//		{
			//			tempP[z] = tempP[z + 1];
			//		}
			//		numP--;
			//	}
			//	else {
			//		for (z = 1; z < numP - 1; z++)//�ߵ���һ����
			//		{
			//			tempP[z] = tempP[z + 1];
			//		}
			//		numP--;
			//	}
			//}
			//Pposes[k + *Pnum] = tempP[numP - 1];
			//k++;
			//}

			else
			{
				while (numP > 1) {
					if (signal[tempP[0]] <= signal[tempP[1]]) {
						for (z = 0; z < numP - 1; z++)//�ߵ���ǰ��
						{
							tempP[z] = tempP[z + 1];
						}
						numP--;
					}
					else {
						for (z = 1; z < numP - 1; z++)//�ߵ���һ����
						{
							tempP[z] = tempP[z + 1];
						}
						numP--;
					}
				}
				Pposes[k + *Pnum] = tempP[numP - 1];
				k++;
				//�����Ƶ���������
				peakfinder->processPeak(false);
				peakfinder->sortByProminence();
				rr_noise_levels[j+ start_idx] = 1-peakfinder->peakSNR(2);
			}
		}
	}
	free(tempP);
	free(SlpP1);
	free(RR);
	free(PR);
	delete peakfinder;
	*Pnum += k;
}

template<class T>
T mean(T* data, int num) {
	T sum = 0;
	T meanValue = 0;
	for (int i = 0; i < num; i++) {
		sum = sum + data[i];
	}
	meanValue = sum / num;
	return meanValue;
}

template<class T>
T maxValue(T* data, int num) {
	T Value = data[0];
	for (int i = 0; i < num; i++) {
		if (data[i] > Value) {
			Value = data[i];
		}
	}
	return Value;
}

template<class T>
T minValue(T* data, int num) {
	T Value = data[0];
	for (int i = 0; i < num; i++) {
		if (data[i] < Value) {
			Value = data[i];
		}
	}
	return Value;

}