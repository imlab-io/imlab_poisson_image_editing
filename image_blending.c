#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <dirent.h>
#include "imcore.h"

// struct to keep four neighbours of the (i,j) and itself
struct poisson_element
{
    int x;
    int y;
    float value[5];
};

// fill the posson element structure
void poisson_element_set(struct poisson_element *out, int x, int y, int idx, float value)
{
    out->x = x;
    out->y = y;
    out->value[idx] = value;
}

// solve I = inv(P)*S
matrix_t* poisson_solver(vector_t *P, matrix_t *S)
{
    int dx[5] = {0, -1, +1, 0, 0};
    int dy[5] = {0, 0, 0, -1, +1};

    float w = 1.8;
    int i, j, k, t, Max_Iter = 4000;

    // get the address of the data
    struct poisson_element *poisson_values = vdata(P, 0);

    // allocate memory for the solution
    matrix_t *I = matrix_create(S, NULL);

    // do SOR iterations
    for (k = 0; k < Max_Iter; k++)
    {
        // do calculation for each pixel inside omega
        for (i = 0; i < length(P); i++)
        {
            // find the additional values coming from the border pixels
            float borders = 0;
            for (j = 1; j < 5; j++)
            {
                if (poisson_values[i].value[j] > 0)
                {
                    borders += poisson_values[i].value[j] * atf(I, poisson_values[i].y + dy[j], poisson_values[i].x + dx[j]);
                }  
            }

            // get the S value @(y,x)
            float S_value = atf(S, poisson_values[i].y, poisson_values[i].x, 0);

            // do single SOR iteration
            atf(I, poisson_values[i].y, poisson_values[i].x) =
                (1 - w) * atf(I, poisson_values[i].y, poisson_values[i].x) + (w / poisson_values[i].value[0]) * (S_value + borders);
        }
    }

    return I;
}

// blend the given two images using the poisson image editing algortihm
return_t poisson_blend(matrix_t *T, matrix_t *S, matrix_t *mask, int x1, int y1, matrix_t *destination)
{
    int i = 0, j = 0, c = 0, k = 0;

    // copy the data of the background into the destination
    matrix_copy(T, destination);

    int dx[5] = {0, -1, +1, 0, 0};
    int dy[5] = {0, 0, 0, -1, +1};
    
    // solve Poisson equation for each channel
    for (c = 0; c < channels(T); c++)
    {
        // create 1d S vectors
        matrix_t *delS = matrix_create(float, rows(S), cols(S), 1);

        // allocate possion matrix (sparse format)
        vector_t *P = vector_create(struct poisson_element);
        
        for (i = 1; i < height(S) - 1; i++)
        {
            for (j = 1; j < width(S) - 1; j++)
            {
                // if (i,j) is inside the mask
                if(atui8(mask, i, j, 0) > 0)
                {
                    struct poisson_element poisson_value = {0};
                    float S_value = 4 * atui8(S, i, j, c);

                    poisson_element_set(&poisson_value, j, i, 0, 4);

                    // compute the average S at(i,j) using the four corner of the pixel
                    for (k = 1; k < 5; k++)
                    {
                        // if we are on border, force I(x,y) = T(x,y)
                        if (!(atui8(mask, i + dy[k], j + dx[k], 0) > 0))
                        {
                            S_value += atui8(T, i + dy[k] + y1, j + dx[k] + x1, c);
                        }
                        else
                        {
                            poisson_element_set(&poisson_value, j, i, k, 1);
                        }

                        // find the lablacian of the S
                        S_value -= atui8(S, i + dy[k], j + dx[k], c);
                    }

                    // set the S and P matrices
                    matrix_set(delS, i, j, 0, &S_value);
                    vector_push(P, &poisson_value);
                }
            }
        }
        
        // solve poisson equation
        matrix_t *I = poisson_solver(P, delS);

        // get pointer to the sparse matrix elements
        struct poisson_element *poisson_values = vdata(P, 0);

        // insert I on destination image
        int t = 0;
        for (t = 0; t < length(P); t++)
        {
            // set the color for the solved image
            atui8(destination, poisson_values[t].y + y1, poisson_values[t].x + x1, c) = clamp(atf(I, poisson_values[t].y, poisson_values[t].x, 0), 0, 255);
        }

        // free vector and possion matrix
        matrix_free(&delS);
        matrix_free(&I);
        vector_free(&P);
    }

    // return success
    return SUCCESS;
}

int main(int argc, char *argv[]) 
{
    if(argc != 4)
    {
        printf("call image_blending folder x y\n");
        return -1;
    }

    // get the file names
    char T_filename[512];
    char S_filename[512];
    char S_mask_filename[512];
    char output_filename[512];

    snprintf(T_filename, 512, "%s//background.bmp", argv[1]);
    snprintf(S_filename, 512, "%s//foreground.bmp", argv[1]);
    snprintf(S_mask_filename, 512, "%s//foreground_mask.bmp", argv[1]);

    // x and y position of the S
    int x = atoi(argv[2]);
    int y = atoi(argv[3]);

    // read S and background images
    matrix_t *T = imread(T_filename);
    matrix_t *S = imread(S_filename);
    matrix_t *S_mask = imread(S_mask_filename);

    if (T == NULL || S == NULL || S_mask  == NULL)
    {
        printf("unable to open input files\n");
        return -1;
    }

    // create an empty matrix
    matrix_t *blended = matrix_create(uint8_t);

    // run the blending algorithm
    poisson_blend(T, S, S_mask, x, y, blended);

    // write the output image
    snprintf(output_filename, 512, "%s//blended_image.bmp", argv[1]);
    imwrite(blended, output_filename);

    matrix_free(&T);
    matrix_free(&S);
    matrix_free(&S_mask);
    matrix_free(&blended);

    return 0;
}