#include<math.h>
#include<stdio.h>
#include<string.h>
#include <stdlib.h>
#include <iostream>
#include <assert.h>
#include <cstdio>  

#ifdef __cplusplus

extern "C"
{
#endif

	void PT_longterm(const float* signal, const int signalsSize, const int* Rposes, const bool* Noisetag, const int RposesSize, const int fs, int *Tposes, int * Tnum, int *Pposes, int *Pnum);

#ifdef __cplusplus
}
#endif 

long get_file_size(const char *path)
{
	long filesize = -1;
	FILE *fp;
	fp = fopen(path, "r");
	if (fp == NULL)
		return filesize;
	fseek(fp, 0L, SEEK_END);
	filesize = ftell(fp);
	fclose(fp);
	return filesize;
}



int read_lead_data(char* file_name, float* ecg_data)
{
	FILE *pfile = NULL;
	pfile = fopen(file_name, "rb");

	fseek(pfile, 0, SEEK_END);
	size_t nFileLen = ftell(pfile);

	fseek(pfile, 0, SEEK_SET);

	fread(ecg_data, 4, nFileLen/4 , pfile);
	fclose(pfile);
	return 0;
}

int read_r_poses(char* file_name, int* r_pos_buf)
{
	FILE *pfile = NULL;
	pfile = fopen(file_name, "rb");

	fseek(pfile, 0, SEEK_END);
	size_t nFileLen = ftell(pfile);
	fseek(pfile, 0, SEEK_SET);

	fread(r_pos_buf, 4, nFileLen / 4, pfile);
	fclose(pfile);
	return 0;
}

int read_r_noise_flags(char* file_name, bool* r_noise_buf)
{
	FILE *pfile = NULL;
	pfile = fopen(file_name, "rb");

	fseek(pfile, 0, SEEK_END);
	size_t nFileLen = ftell(pfile);
	fseek(pfile, 0, SEEK_SET);

	int * label_noise_arr = new int[nFileLen / 4];

	fread(label_noise_arr, 4, nFileLen / 4, pfile);
	fclose(pfile);

	for (int i = 0; i < nFileLen / 4; i++)
	{
		r_noise_buf[i] = label_noise_arr[i];
	}
	delete[] label_noise_arr;

	return 0;
}


void test_PT_Detect()
{
	//char lead_II_file[] = "F:\\tmp\\data\\out_II_data.dat";
	//char r_pos_file[] = "F:\\tmp\\data\\out_merge_r_pos.dat";
	//char r_noise_flag_file[] = "F:\\tmp\\data\\r_noise_flag.dat";

	printf("00000");

	char lead_II_file[] = "/home/chuanyhu/txt_out/out_II_data.dat";
	char r_pos_file[] = "/home/chuanyhu/txt_out/out_merge_r_pos.dat";
	char r_noise_flag_file[] = "/home/chuanyhu/txt_out/r_noise_flag.dat";

	unsigned long rfile_size = get_file_size(r_pos_file);
	int r_pos_len = rfile_size / 4;
	int * r_pos_arr = new int[r_pos_len];
	bool * r_noise_flag_arr = new bool[r_pos_len];
	read_r_poses(r_pos_file, r_pos_arr);
	read_r_noise_flags(r_noise_flag_file, r_noise_flag_arr);

	
	unsigned long lead_II_file_size = get_file_size(lead_II_file);
	int data_len = lead_II_file_size / 4;
	float* lead_data = new float[data_len];
	read_lead_data(lead_II_file, lead_data);

	int t_pos_num = 0;
	int p_pos_num = 0;
	int * t_pos_arr = new int[data_len];
	int * p_pos_arr = new int[data_len];

	printf("111111111");

	PT_longterm(lead_data, data_len , r_pos_arr, r_noise_flag_arr, r_pos_len, 250, t_pos_arr, &t_pos_num, p_pos_arr, &p_pos_num);
	printf("2222");
	for (int i = 0; i < p_pos_num; i++)
	{
		//if(p_pos_arr[i] > 13945000 && p_pos_arr[i] < 13947000)
		printf("%d\n", p_pos_arr[i]);
	}
	printf("3333");


	delete[] r_pos_arr;
	delete[] r_noise_flag_arr;
	delete[] lead_data;
	delete[] t_pos_arr;
	delete[] p_pos_arr;


}


int main(int argc, char **argv) {

	test_PT_Detect();

	return 0;
}

