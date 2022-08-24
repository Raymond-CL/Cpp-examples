#include <iostream>
#include <stdlib.h>
#include <time.h>

#ifdef __unix__
      #define OpSys 0
      #include <unistd.h>
#elif defined(__WIN32) || defined(WIN32)
      #define OpSys 1
      #include <windows.h>
#endif

using namespace std;

const int WIDTH = 10;
const int HEIGHT = 10;
const int IT = 500;
const int BUFF = 50000;

void init(bool grid[HEIGHT][WIDTH]){
      for(int i=0;i<HEIGHT;i++)
            for(int j=0;j<WIDTH;j++)
                  grid[i][j]=rand()%2;
}

void check(bool grid_a[HEIGHT][WIDTH],int grid_b[HEIGHT][WIDTH]){
      for(int i=0;i<HEIGHT;i++)
            for(int j=0;j<WIDTH;j++){
                  int count=0;
                  if(i>0 && j>0              && grid_a[i-1][j-1]) count++;
                  if(i>0                     && grid_a[i-1][j]  ) count++;
                  if(i>0 && j<WIDTH-1        && grid_a[i-1][j+1]) count++;
                  if(j>0                     && grid_a[i][j-1]  ) count++;
                  if(j<WIDTH-1               && grid_a[i][j+1]  ) count++;
                  if(i<HEIGHT-1 && j>0       && grid_a[i+1][j-1]) count++;
                  if(i<HEIGHT-1              && grid_a[i+1][j]  ) count++;
                  if(i<HEIGHT-1 && j<WIDTH-1 && grid_a[i+1][j+1]) count++;
                  grid_b[i][j]=count;
            }
}

int main(){

      bool A[HEIGHT][WIDTH]={};
      int B[HEIGHT][WIDTH]={};
      srand(time(0));
      init(A);
      check(A,B);
      int end=0;

      if(OpSys==0) system("clear");
      if(OpSys==1) system("cls");

      for(int i=0;i<HEIGHT;i++){
      for(int j=0;j<WIDTH;j++){
            if(A[i][j]) cout << '.';
            else cout << ' ';
            cout << ' ';
      }
      cout << endl;
      }
      cout << "Press any key to proceed with the above initial state...";
      cin.ignore(cin.rdbuf()->in_avail()+1);

      while(end<IT){
            if(OpSys==0) system("clear");
            if(OpSys==1) system("cls");

            for(int i=0;i<HEIGHT;i++){
            for(int j=0;j<WIDTH;j++){
                  if(A[i][j]) cout << '.';
                  else cout << ' ';
                  cout << ' ';

                  int n=B[i][j];
                  if(A[i][j] && (n==2 || n==3) ) A[i][j]=1;
                  else if(!A[i][j] && n==3) A[i][j]=1;
                  else A[i][j]=0;
            }
            cout << endl;
            }

            check(A,B);
            
#ifdef __unix__
			usleep(BUFF);
#elif defined(__WIN32) || defined(WIN32)
      		Sleep(BUFF/1000);
#endif

            end++;
      }

      return 0;
}
