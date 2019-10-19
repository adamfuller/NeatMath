#include <stdio.h>
#include <stdlib.h>

void rotateRight(int array[], int width);
void rotateLeft(int array[], int width);

/**
 * Rotates the array 90° to the left (CCW)
 * as if it were a width by width matrix
 * in the order:
 * 0 1 2
 * 3 4 5
 * 6 7 8
 * 
 * TODO: reduce from using second loop
 **/
void rotateLeft(int array[], int width)
{
    char *newSpots = malloc(sizeof(int) * width * width);
    int newIndex;
    for (int i = 0; i < width * width; i++)
    {
        // Rotates to the left, counter clock wise
        newIndex = width - (int)(i / width) - 1 + (i % width) * width;
        newSpots[i] = array[newIndex];
    }

    // Swap over the new array's values
    for (int i = 0; i < width * width; i++)
        array[i] = newSpots[i];

    // Don't forget to free up memory...
    free(newSpots);
}

/**
 * Rotates the array 90° to the right (CW)
 * as if it were a width by width matrix
 * in the order:
 * 0 1 2
 * 3 4 5
 * 6 7 8
 * 
 * TODO: reduce from using second loop
 **/
void rotateRight(int array[], int width)
{
    char *newSpots = malloc(sizeof(int) * width * width);
    int newIndex;
    for (int i = 0; i < width * width; i++)
    {
        // Rotates to the left, counter clock wise
        newIndex = (int)(i / width) + (width - (i % width) - 1) * width;
        newSpots[i] = array[newIndex];
    }
    for (int i = 0; i < width * width; i++)
    {
        array[i] = newSpots[i];
    }
    // Don't forget to free up memory...
    free(newSpots);
}

int main(void){
    int t[] = {0, 1, 2, 3};
    rotateRight(t, 2);
    for (int i = 0; i<4;i++){
        printf("%d\n", t[i]);
    }
}