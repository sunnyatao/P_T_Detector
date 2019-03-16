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
	const int RR_count = MAX_int32_T; //定义最长片段RR的个数
	int32_T Max_interval = 15 * fs; //超过此时间不检测PT
	*Tnum = 0; //从零开始记录
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
			fil_HP_buff[i - fil_HP_delay] = HP_V;                 //修正群时延
		}
		if (i - fil_HP_delay - fil_HL_delay >= 0) {
			fil_HL_buff[i - fil_HL_delay - fil_HP_delay] = fil_HL->output;			//hilbert transform
		}
	}
	for (int i = 0; i < signalsSize; i++) {
		fil_HL_buff[i] = sqrt(fil_HL_buff[i] * fil_HL_buff[i] + fil_HP_buff[i] * fil_HP_buff[i]);  //求包络
	}
	//float meanValue = mean(fil_HL_buff, signalsSize);
	//for (int i = 0; i < signalsSize; i++) {
	//	fil_HL_buff[i] = fil_HL_buff[i] - meanValue;										//去直流信号
	//}

	//低频噪声初始化
	for (int numSNR = 0; numSNR < RposesSize;numSNR++) {
		rr_noise_levels[numSNR] = (float)1;
	}

	while (start_idx < RposesSize) {
		while (start_idx < RposesSize) {
			//找非噪音起点
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

			//找噪音终点
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
	*current_Tnum = *Tnum;//记录本次记录T的offset，给Ppeak参照
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
	double SlpT1_treh;//T波正阈值
	double SlpT2_treh;//T波负阈值
	int32_T *tempT;
	int numT;
	findpeak * peakfinder = new findpeak;
	int tmp_buf_len = Rposes[RposesSize-1] - Rposes[0];

	int32_T *RR = (int32_T*)malloc((RposesSize) * sizeof(int32_T));
	SlpT1 = (double*)malloc(tmp_buf_len * sizeof(double));	//正T波高度
	SlpT2 = (double*)malloc(tmp_buf_len * sizeof(double));	//负T波高度
	tempT = (int32_T*)malloc(tmp_buf_len * sizeof(int32_T));

	

	//用一半的数据进行T波阈值训练
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
	if (N > 10) { N = 10; }									//          <--------------------------------------------对长程数据是否需要修改
	k1 = 0;
	k2 = 0;

	for (int j = 0; j < N - 1; j++) {
		//RR[j] = Rposes[j + 1] - Rposes[j];
		stop_index = floor(Rposes[j] +0.5*fs);			//设定每次的RR检测边界
		start_index = floor(Rposes[j] + 30);
		if (start_index > signalsSize) {
			start_index = signalsSize;
		}
		if (stop_index < 0) {
			stop_index = 0;
		}
		peakfinder->reset();
		for (int i = start_index; i < stop_index; i++) {
			peakfinder->inputdata(signal[i]);				//在每次的RR之间寻找最大值
		}
		//peakfinder->processdata();//寻找波峰, 设为可疑的T点

		for (int i = 0; i < peakfinder->max.size(); i++)
		{
			B_index = peakfinder->max_loc[i] + start_index;
			tempSlpTa = signal[B_index] - signal[B_index + 5];			//通过极值左右12点位置计算斜率
			tempSlpTc = signal[B_index] - signal[B_index - 5];
			if (tempSlpTa > 0 && tempSlpTc > 0)							//局部最大值
			{
				if (k1 > tmp_buf_len) 
				{
					printf("k1:%d len:%d\n", k1, tmp_buf_len);
				}
				SlpT1[k1] = (tempSlpTa + tempSlpTc) / 2;				//求出每个点的平均斜率
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
		//		//局部最小值
		//		SlpT2[k2] = (tempSlpTa + tempSlpTc) / 2;
		//		k2++;
		//	}
		//}
	}
	//训练得到T波正负阈值(0.2倍最斜率，正负阈值分别计算)
	SlpT1_treh = 0.05*maxValue(SlpT1, k1);
	//SlpT2_treh = 0.05*maxValue(SlpT2, k2);
	//用阈值进行T波检测
	
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
			peakfinder->inputdata(signal[i]);//寻找区间内的波峰，设为可疑的正T波
		}
		numT = 0;
		for (int i = 0; i < peakfinder->max.size(); i++)
		{
			B_index = peakfinder->max_loc[i] + start_index;
			tempSlpTa = signal[B_index] - signal[B_index +5];				//通过极值左右12点位置计算斜率
			tempSlpTc = signal[B_index] - signal[B_index -5];
			/*if (tempSlpTa > 0 && tempSlpTc > 0) {	*/						//局部最大值
				//if ((tempSlpTa + tempSlpTc) / 2 > SlpT1_treh)				//局部最大值的平均高度大于阈值(可用prominence替代)
				//{
					if (numT > tmp_buf_len)
					{
						printf("1111numT:%d len:%d\n", numT, tmp_buf_len);
					}
					tempT[numT] = B_index;									//存储所有可疑正T波
					numT++;
				//}
		/*	}*/
		}
		//寻找波谷，设为可疑的负T波
		//for (int i = 0; i < peakfinder->min.size(); i++)
		//{
		//	B_index = peakfinder->min_loc[i] + start_index;
		//	tempSlpTa = signal[B_index] - signal[B_index + mean_RR / 20];
		//	tempSlpTc = signal[B_index] - signal[B_index - mean_RR / 20];
		//	if (tempSlpTa < 0 && tempSlpTc < 0)// 用训练的T波斜率阈值进行验证
		//	{
		//		if (fabs(tempSlpTa + tempSlpTc) / 2 > fabs(SlpT2_treh))
		//		{
		//			if (numT > tmp_buf_len)
		//			{
		//				printf("222numT:%d len:%d\n", numT, tmp_buf_len);
		//			}
		//			tempT[numT] = B_index;									//存储所有可能负T波
		//			numT++;
		//		}
		//	}
		//}
		//踢掉一个RR间期内多检出的T点，用T点的原始信号值进行验证（保留最大T）
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
		//int t = 0;															// <-----------------------------这个变量好像没有用？
		//
		//while (numT > 1) {
		//	if (fabs(signal[tempT[0]]) <= fabs(signal[tempT[1]])) {		// 如果t点信号幅度小于等于下一个t点幅值(第一个T)
		//		for (z = 0; z < numT - 1; z++)//踢掉当前点					//从这个t点开始
		//		{
		//			tempT[z] = tempT[z + 1];								//将整个数组向前shift 1
		//		}
		//		numT--;														//减少数据size
		//	}
		//	else {
		//		for (z = 1; z < numT - 1; z++)//踢掉后一个点			//如果t点信号幅度大
		//		{
		//			tempT[z] = tempT[z + 1];								//将整个数据从t+1向前shift 1
		//		}
		//		numT--;														//减小数据size
		//	}
		//}
		//Tposes[j + *Tnum] = tempT[0];  //在队尾增加数据
	}
	*Tnum += RposesSize - 1;			  //队列长度+1
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
	float *SlpP1 = (float*)malloc(tmp_buf_len * sizeof(float));   //正P波高度
	
	//用一半的数据进行T波阈值训练
	//FILE *tfile = NULL;
	//tfile = fopen("E:\\matlabWorkspace\\testdata\\冯殿玉\\filter_buf", "wb+");
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
			peakfinder->inputdata(signal[i]);//寻找区间内的波峰，设为可疑的正T波
		}

		for (int i = 0; i < peakfinder->max.size(); i++)
		{
			B_index = peakfinder->max_loc[i] + start_index;
			if (B_index <= Tposes[j + *Tnum] + mean_RR / 10) {				//如果可疑P点出现在T波之前，删掉
				continue;
			}
			//if (abs(signal[B_index]) > 0.7*abs(signal[Rposes[j]])) {
			//	continue;
			//}
			tempSlpPa = signal[B_index] - signal[B_index + 5];			//通过极值左右15点位置计算斜率
			tempSlpPc = signal[B_index] - signal[B_index - 5];
			if (tempSlpPa > 0 && tempSlpPc > 0)
			{
				SlpP1[k1] = (tempSlpPa + tempSlpPc) / 2;
				k1++;
			}
		}
	}
	//训练得到P波阈值
	SlpP1_treh = 0.05*maxValue(SlpP1, k1);
	//用阈值进行P波检测
	int k = 0;

	for (int j = 0; j < RposesSize - 1; j++)
	{
		//if(Rposes[j]>= 2764730){
		//	printf("=");
		//}
		//RR[j] = Rposes[j + 1] - Rposes[j];
		start_index = floor(Rposes[j] + 0.1*fs);						//寻找区间内的波峰，设为可疑的P波
		stop_index = floor(Rposes[j + 1] - mean_RR / 15);
		if (start_index > signalsSize) {
			start_index = signalsSize;
		}
		if (stop_index < 0) {
			stop_index = 0;
		}

		peakfinder->reset();
		for (int i = start_index; i < stop_index; i++) {
			peakfinder->inputdata(signal[i]);							//寻找区间内的波峰，设为可疑的正P波
		}
		numP = 0;
		
		for (int i = 1; i < peakfinder->max.size(); i++)				//从第二个峰值开始找P
		{
			B_index = peakfinder->max_loc[i] + start_index;
			if (B_index <= Tposes[j + *Tnum]) {				//如果可疑P点出现在T波之前，删掉(略过本峰值)
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

		//  踢点策略，传到阻滞有两个P点，非传到阻滞只有一个P点
		int i = 0;
		int z = 0;
		int t = 0;

		//如果第一个RR间期不是传导阻滞，用幅值大小选择p点 
		if (numP < 1) {
			continue;
		}
		if (j < 2) {													//处理前两个RR间期的P点
			while (numP > 1) {
				if (signal[tempP[0]] <= signal[tempP[1]]) {
					for (z = 0; z < numP - 1; z++)//踢掉当前点
					{
						tempP[z] = tempP[z + 1];
					}
					numP--;
				}
				else {
					for (z = 1; z < numP - 1; z++)//踢掉后一个点
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

		else {	//从第三个RR间期开始判断传导阻滞
			//int32_T mean_temp = mean(RR, j + 1);
			if ((RR[j] >= 1.2 * mean_RR|| RR[j]>1.2*fs) && numP > 1) {				//如果RR间期是传导阻滞，不限制P点个数用幅值大小和RR间期来选择P点
				while (t < numP - 1) {
					if (abs((int32_T)(tempP[t + 1] - tempP[t])) < 0.4 * mean_RR) {
						if (fabs(signal[tempP[t]]) <= fabs(signal[tempP[t + 1]])) {
							for (z = t; z < numP - 1; z++)//踢掉当前点
							{
								tempP[z] = tempP[z + 1];
							}
							numP = numP - 1;
						}
						else {
							for (z = t + 1; z < numP - 1; z++)//踢掉后一个点
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
				//计算待选P点的均值
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
						continue;										//排除长RR间期内多检出幅值很低的波峰
					}
					Pposes[k + *Pnum] = tempP[x];
					k++;
					y++;
				}
				//计算低频噪声信噪比
				peakfinder->processPeak(false);
				peakfinder->sortByProminence();
				rr_noise_levels[j + start_idx] = 1-peakfinder->peakSNR(y+1);
			}

			//如果RR间期不是传导阻滞，只有一个P点则利用PR间期与幅值来选择P点位置
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
			//				for (z = t; z < numP - 1; z++)//踢掉当前点
			//				{
			//					tempP[z] = tempP[z + 1];
			//				}
			//				numP--;
			//			}
			//			else {
			//				if (abs(signaldiff)*10000 > abs(PRdiff - PR1diff)) {
			//					for (z = t; z < numP - 1; z++)//踢掉当前点
			//					{
			//						tempP[z] = tempP[z + 1];
			//					}
			//					numP--;
			//				}
			//				else {
			//					for (z = t + 1; z < numP - 1; z++)//踢掉后一个点
			//					{
			//						tempP[z] = tempP[z + 1];
			//					}
			//					numP--;
			//				}
			//			}
			//		}
			//		else {
			//			if (PRdiff <= PR1diff) {
			//				for (z = t + 1; z < numP - 1; z++)//踢掉后一个点
			//				{
			//					tempP[z] = tempP[z + 1];
			//				}
			//				numP--;
			//			}
			//			else {
			//				if (abs(signaldiff)*10000 >= abs(PRdiff - PR1diff)) {
			//					for (z = t + 1; z < numP - 1; z++)//踢掉后一个点
			//					{
			//						tempP[z] = tempP[z + 1];
			//					}
			//					numP--;
			//				}
			//				else {
			//					for (z = t; z < numP - 1; z++)//踢掉当前点
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
			//		for (z = 0; z < numP - 1; z++)//踢掉当前点
			//		{
			//			tempP[z] = tempP[z + 1];
			//		}
			//		numP--;
			//	}
			//	else {
			//		for (z = 1; z < numP - 1; z++)//踢掉后一个点
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
						for (z = 0; z < numP - 1; z++)//踢掉当前点
						{
							tempP[z] = tempP[z + 1];
						}
						numP--;
					}
					else {
						for (z = 1; z < numP - 1; z++)//踢掉后一个点
						{
							tempP[z] = tempP[z + 1];
						}
						numP--;
					}
				}
				Pposes[k + *Pnum] = tempP[numP - 1];
				k++;
				//计算低频噪声信噪比
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