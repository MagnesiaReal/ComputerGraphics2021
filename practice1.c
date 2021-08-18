#include <stdio.h>
#include <stdlib.h>

int x, y;
int x1, x2, y1, y2;


void table (int xmin, int xmax){
    printf("X\tY\n");
    while (x <= xmax){
        float m = ((float)y2 - y1)/ ((float)x2 - x1);
        y = m*(x - x1) + y1;
        printf("%i\t%i\n", x, y);
        x ++;
    }
}

int main (){
    
    printf("X1 = "); scanf("%d", &x1);
    printf("Y1 = "); scanf("%d", &y1);
    printf("X1 = "); scanf("%d", &x2);
    printf("Y1 = "); scanf("%d", &y2);
    
    
    if (x1 <= x2){
        x = x1;
        table (x1, x2);
    } 
    else {
        x = x2;   
        table (x2, x1);
    }

}