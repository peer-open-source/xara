
#include <stdio.h>
#include <CrossSection.h>


int main()
{
  // using enum FrameType;

  int type = N|Vy|Vz|T|My|Mz;

  printf("%b, %b\n", type, log2(FrameType::End));

  for (int i=1; i<FrameType::End; i*=2) {
    printf("%4d %12b %12b ", i, i, type&i);
    switch (type & i) {
      case N:
        printf(": N\n");
        break;

      case Vy:
        printf(": Vy\n");
        break;

      default:
        printf("\n");
    }
  }
}

