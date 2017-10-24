#include <cstdio>
#include <cstring>

extern "C" void f_main(int length,char* filename);

int main(int argc, char* argv[]) {
  if(argc==1) {
    printf("PROBLEM WITH THE INPUT-FILE NAME! QUITTING!\n");
    return 1;
  }
  f_main(strlen(argv[1]),argv[1]);
  return 0;
}
