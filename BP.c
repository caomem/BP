#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void usage(char *exec){
	printf("%s -solve <input_file>\n", exec);
}

void solve(float **mat){
    int i, j;
    for (i = 0; i < 3000; i++)
    {
        for (j = 0; j < 3000; j++)
        {
            if(mat[i][j] != 0){
                printf(" %f ", mat[i][j]);
            }
        }
        
    }
    
}

int main(int argc, char **argv){
	if(argc<3){
		usage(argv[0]);
	} else{
        if( !strcmp(argv[1], "-solve")){
            FILE* file =  NULL;
            if( (file = fopen(argv[2], "r")) == NULL){
                /* error openning the file */
                perror("fopen: ");
                return 1;
            }
            int n, i,j = 0, k = 0;
            char buff[5], buffFloat[16];
            if (fread(buff, sizeof(char), 5, file) < 5 || strcmp(buff, "E = [")){
                perror("file out of params: ");
                return 1;
            }
            float **mat = (float**)calloc(3000, sizeof(float*));
            for (i = 0; i < 3000; i++)
            {
                mat[i] = (float*)calloc(3000, sizeof(float*));
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