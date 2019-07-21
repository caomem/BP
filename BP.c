#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_blas.h>

int len = 3000;

typedef struct _vertex{
    float x, y, z;
} vertex;

void usage(char *exec){
	printf("%s -solve <input_file> <count_vertex>\n", exec);
}

void printMat(float **mat){
    int i, j;
    for (i = 0; i < len; i++)
    {
        for (j = 0; j < len; j++)
        {
            if(mat[i][j] != 0){
                printf(" %f ", mat[i][j]);
            }
        }
        
    }
    
}

void solve(float **mat){
    vertex *v = (vertex*)calloc(len, sizeof(vertex));

    v[1] = 

    int i,j;
    for (i = 0; i < len; i++)
    {
        for (j = 0; j < len; j++)
        {
            if(mat[i][j] != 0){
                printf(" %f ", mat[i][j]);
            }
        }
        
    }
}

int main(int argc, char **argv){
	if(argc<4){
		usage(argv[0]);
	} else{
        if( !strcmp(argv[1], "-solve")){
            FILE* file =  NULL;
            if( (file = fopen(argv[2], "r")) == NULL){
                /* error openning the file */
                perror("fopen: ");
                return 1;
            }

            if( (len = atoi(argv[3])) == 0){
                /* error openning the file */
                perror("count_vertex param error: ");
                return 1;
            }

            int n, i,j = 0, k = 0;
            char buff[5], buffFloat[16];
            if (fread(buff, sizeof(char), 5, file) < 5 || strcmp(buff, "E = [")){
                perror("file out of params: ");
                return 1;
            }
            float **mat = (float**)calloc(len, sizeof(float*));
            for (i = 0; i < len; i++)
            {
                mat[i] = (float*)calloc(len, sizeof(float*));
            }
            i = 0;
            while((n = fread(buff, sizeof(char), 1, file)) == 1){
                if (buff[0] == ' ')
                {
                    if (k > 0){
                        mat[i][j] = atof(buffFloat);
                        memset(buffFloat, 0, sizeof(buffFloat));
                        k = 0;
                        j++;
                    }
                    continue;
                }else if (buff[0] == '0')
                {
                    j++;
                }else if (buff[0] == ';'){
                    j = 0;
                    i++;
                }else {
                    buffFloat[k++] = buff[0];
                }
            }
            solve(mat);
        }
    }
}