#include<iostream>
using namespace std;
#include<cstring>
#include<iomanip>
#define size 6
template<typename T, typename std::size_t N>
void LU(const T (&arr)[N][N], T (&L)[N][N],T (&U)[N][N]);

template<typename T, typename std::size_t N>
void cal_Lb(const T (&L)[N][N], const T (&b)[N], T (&ans)[N]);

template<typename T, typename std::size_t N>
void cal_Ub(const T(&U)[N][N], const T(&b)[N], T (&ans)[N]);

int main(){
	//========================================
	//方程组 
	double a[size][size] = {4, -1, 0, -1, 0, 0,
						-1, 4, -1,  0, -1,0,
						 0, -1, 4,  0,  0,-1,
						 -1, 0, 0,  4, -1, 0,
						  0,-1, 0, -1,  4,-1,
						  0, 0,-1,  0, -1, 4};
	double b[size] = {2,3,2,2,1,2};
	double L[size][size],U[size][size];
	double ans[size];                              
	//============================================
	//LU分解 
	LU(a, L, U);
	//=========================================
	//输出原矩阵
	cout<<"原矩阵为：\n";
	for(int i = 0; i < size; ++i){
		for(int j = 0; j < size; ++j){
			cout<<setw(18)<<a[i][j];
		}
		cout<<'\n';
	} 
	
	//
	cout<<"进行LU分解：\n\n"; 
	cout<<setprecision(6)<<setiosflags(ios::left);
	cout<<"L为:\n";
	for(int i = 0; i < size; ++i){
		for(int j = 0; j < size; ++j){
			cout<<setw(18)<<L[i][j];
		}
		cout<<'\n';
	}
	
	cout<<"\n\nU为：\n";
	for(int i = 0; i < size; ++i){
		for(int j = 0; j < size; ++j){
			cout<<setw(18)<<U[i][j];
		}
		cout<<'\n';
	}
	//============================================
	//计算ans = U的逆乘L的逆乘b 
	cal_Lb(L,b,ans);
	cal_Ub(U,ans,ans);
	//===========================================
	//打印ans
	cout<<"\n\n原方程通过LU分解得到的解为：\n";
	for(int i = 0; i < size; ++i){
		cout<<"x"<<i<<" = "<<ans[i]<<'\n'; 
	} 
	return 0;
} 

/*-----------------------------------
求任意数字类型的LU分解

注意： 直接对传入的L、U进行操作 
-------------------------------------*/
template<typename T, typename std::size_t N>
void LU(const T (&arr)[N][N], T (&L)[N][N],T (&U)[N][N]){                      
	
	memset(L, 0, sizeof(T) * N * N);
	memset(U, 0, sizeof(T) * N * N);
	 
	for(int i = 0; i < N; ++i){
		L[i][i] = 1;
	}

	for(int i = 0 ; i < N ; ++i){ 
		
		for(int j = i; j < N; ++j){           //处理第u[i:*] 
			U[i][j] = arr[i][j];
			
			for(int k = 0; k < i; ++k){
				U[i][j] -= L[i][k] * U[k][j]; 
			}
		
		}
		
		for(int j = i + 1; j < N; ++j){      //处理L[* : i] 
			L[j][i] = arr[j][i];
			
			for(int k = 0; k < i ; ++k){
				L[j][i] -= L[j][k] * U[k][i];
			}
			
			L[j][i] /= U[i][i];
	
		}
	}
}
/*-------------------------------------
计算L的逆乘b

结果保存在ans中 
---------------------------------------*/
template<typename T, typename std::size_t N>
void cal_Lb(const T (&L)[N][N], const T (&b)[N], T (&ans)[N]){
	for(int i = 0; i < N; ++i){
		ans[i] = b[i];
		for(int k = 0; k < i; ++k){
			ans[i] -= L[i][k] * ans[k];
		}
	}	
} 

/*------------------------------------
计算U的逆乘一个矩阵 

结果保存在ans中 
-------------------------------------*/
template<typename T, typename std::size_t N>
void cal_Ub(const T(&U)[N][N], const T(&b)[N], T (&ans)[N]){
	for(int i = N - 1; i >= 0; --i){
		ans[i] = b[i];
		for(int k = i + 1; k < N; ++k){
			ans[i] -= U[i][k] * b[k];
		}
		ans[i] /= U[i][i];
	}
}
