/**
 *exp�ڡ�0��1���ϵĻ��� 
 *��ȷ��Ϊ�� e - 1  
**/ 
#include<iostream>
#include<iomanip>
#include<cmath>
#include<map>
#include<cstring>
using namespace std;
#define Wu 0.000005
#define Romberg_Size 30
 
struct data{
	double sum;
	int num;
	double wu;
	data(double s, int n, double w):sum(s),num(n),wu(w){}
};

data tiXing(); 
data Simpson();
data Romberg();
data Gauss();
int main(){
	printf("����exp������[0,1]�ϵĻ���\n\n");
	data d = tiXing();
	printf("���ι�ʽ����%d�Σ����ȿɴﵽ %.6f�����ֽ��Ϊ��%.10f, �ض����Ϊ��%.10f\n\n", d.num, Wu, d.sum, d.wu);
	
	d = Simpson(); 
	printf("Simpson��ʽ����%d�Σ����ȿɴﵽ %.6f�����ֽ��Ϊ��%.10f, �ض����Ϊ��%.10f\n\n", d.num, Wu, d.sum, d.wu);
	
	d = Romberg();
	printf("Romberg���Ʒ����㵽 k = %d�����ȿɴ�%.6f�����ֽ��Ϊ��%.10f, �ض����Ϊ��%.10f\n\n", d.num, Wu, d.sum, d.wu);
	
	d = Gauss();
	printf("Gauss�����ʽ���㵽 n = %d�����ȿɴ�%.6f�����ֽ��Ϊ��%.10f, �ض����Ϊ��%.10f\n\n", d.num, Wu, d.sum, d.wu);
	
	return 0;
}

data tiXing(){
	int n = 0;               //�ֳ�n�� 
	double wu = 1;	
	double h , sum;
	while(wu > Wu){
		++n;
		h = 1.0 / n, sum = 0; 
		
		for(double index = 0; index <= 1 - h; index += h){
			sum += (exp(index) + exp(index + h)) * h / 2;
		}
		
		wu = abs(sum - exp(1) + 1);
	} 
	
	return data(sum, n,abs(sum - exp(1) + 1));
}

data Simpson(){
	int n = 0;
	double wu = 1;
	double h, sum, jian;
	while(wu > Wu){
		++n;
		h = 1.0 / n, sum = 0, jian = h / 2;
		for(double index = 0; index <= 1 - h; index += h){
			sum += (exp(index) + 4 * exp(index + jian) + exp(index + h)) * h / 6;
		}
		
		wu = abs(sum - exp(1) + 1);
	}
	return data(sum, n,abs(sum - exp(1) + 1));
}

data Romberg(){
	double T[Romberg_Size][Romberg_Size];
	memset(T, 0x3f, sizeof(double) * Romberg_Size * Romberg_Size);
	bool flag = true;
	double jian = 0.5;                                //��ǰ�����С 
	double tmp = 0;                                   //���ڼ���ÿ�� T��i����0�����¼ӵ����� 
	double step = 1;                                 //��ʼ�ĵط� , �൱����һ�εļ�� 
	T[0][0] = (1 + exp(1)) / 2.0;                    //��ʼ����һ������ 

	int i; 
	for(i = 1; flag; ++i){
		tmp = 0;
		for(double j = jian; j < 1; j += step){
			tmp += exp(j);
		} 

		T[i][0] = T[i - 1][0] / 2 + jian * tmp;     //���㵱ǰT[i][0] 
		
		step = jian;                                //���� step �� jian 
		jian /= 2;
		
		for(int j = 1; j <= i; ++j){            
			T[i][j] = (pow(4, j) * T[i][j - 1] - T[i - 1][j - 1]) / (pow(4, j) - 1); 
		}
		
		if(T[i][i] != 0x3f && abs(T[i][i] - (exp(1) - 1))< Wu){   //�ж��Ƿ�ﵽ���� 
			flag = false;
		}
		
	}
	return data(T[i - 1][i - 1], i - 1,abs(T[i - 1][i - 1] - (exp(1) - 1)) );
	
}

double f(double x){
	return exp(0.5 * x + 0.5);
}
data Gauss(){
	map<int, map<double, double> > m;
	m[0][2] = 0;
	m[1][1] = 0.5773502692;
	m[2][0.5555555556] = 0.7745966692;
	m[2][0.8888888889] = 0;
	m[3][0.3478548451] = 0.8611363116;
	m[3][0.6521451549] = 0.3399810436;
	m[4][0.2369268851] = 0.9061798459;
	m[4][0.4786286705] = 0.5384693101;
	m[4][0.5688888889] = 0;
	int n = 0;
	double ans = 0;
	double wu = 1;
	while(wu > Wu){
		ans = 0;
		for(map<double, double>::iterator iter = m[n].begin(); iter != m[n].end(); ++iter){
			ans += iter->first * f(iter->second);
			if(iter->second) ans += iter->first * f(-1 * iter->second);
		}
		wu = abs(ans / 2 - (exp(1) - 1));
		++n;
	}
	
	return data(ans / 2, n, abs(ans / 2 - (exp(1) - 1))); 
}
