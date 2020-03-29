#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include<time.h>
#include <gsl/gsl_blas.h>

double epsl = 0.01;

typedef struct {
    double x, y, z;
} position;

typedef struct _vertex {
    position pos;
    gsl_matrix_view *C;
    struct _vertex *p,*l,*r;
} vertex;

void usage(char *exec) {
    printf("%s -s <input_file> <count_vertex>\n", exec);
    printf("%s -c <input_file>\n", exec);
    printf("%s -show <input_file>\n", exec);
}

int* ddd(char *fileName, int len){
    FILE *file = NULL;
    if ((file = fopen(fileName, "r")) == NULL) {
        /* error openning the file */
        perror("fopen: ");
        return NULL;
    }
    int n, i, j = 0, k = 0;
    char buff[6], buffFloat[17];
    memset(buff, 0, sizeof(buff));
    memset(buffFloat, 0, sizeof(buffFloat));
    if (fread(buff, sizeof(char), 5, file) < 5 || strcmp(buff, "E = [")) {
        perror("file out of params: ");
        return NULL;
    }
    
    double **mat = (double **)calloc(len, sizeof(double *));
    for (i = 0; i < len; i++) {
        mat[i] = (double *)calloc(3, sizeof(double *));
    }
    i = 0;
    while ((n = fread(buff, sizeof(char), 1, file)) == 1) {
        if (buff[0] == ' ') {
            if (k > 0) {
                mat[i][j] = atof(buffFloat);
                memset(buffFloat, 0, sizeof(buffFloat));
                k = 0;
                j++;
            }
            continue;
        }
        else if (buff[0] == ';') {
            if (k) {
                mat[i][j] = atof(buffFloat);
                memset(buffFloat, 0, sizeof(buffFloat));
                k = 0;
            }
            j = 0;
            i++;
        }
        else {
            buffFloat[k++] = buff[0];
        }
    }
    fclose(file);

    int current = 0; 
    for (i = 0; i < len; i++) {
        
        for (j = 0; j < len; j++) {
            printf(" %f ", sqrt(pow(mat[i][0] - mat[j][0], 2) + pow(mat[i][1] - mat[j][1], 2) + pow(mat[i][2] - mat[j][2], 2)));
            fflush(stdout);
        }
        printf("\n");
    }
}

void printMat(double **mat, int len) {
    int i, j;
    for (i = 0; i < len; i++) {
        for (j = 0; j < len; j++) {
            //if (mat[i][j] != 0) {
                printf(" %f ", mat[i][j]);
                fflush(stdout);
            //}
        }
        printf("\n");
    }
}

double **matrixConstruct(char *fileName, int len){
    FILE *file = NULL;
    if ((file = fopen(fileName, "r")) == NULL) {
        /* error openning the file */
        perror("fopen: ");
        return NULL;
    }
    int n, i, j = 0, k = 0;
    char buff[6], buffFloat[17];
    memset(buff, 0, sizeof(buff));
    memset(buffFloat, 0, sizeof(buffFloat));
    if (fread(buff, sizeof(char), 5, file) < 5 || strcmp(buff, "E = [")) {
        perror("file out of params: ");
        return NULL;
    }
    
    double **mat = (double **)calloc(len, sizeof(double *));
    for (i = 0; i < len; i++) {
        mat[i] = (double *)calloc(len, sizeof(double *));
    }
    i = 0;
    while ((n = fread(buff, sizeof(char), 1, file)) == 1) {
        if (buff[0] == ' ') {
            if (k > 0) {
                mat[i][j] = atof(buffFloat);
                memset(buffFloat, 0, sizeof(buffFloat));
                k = 0;
                j++;
            }
            continue;
        }
        else if (buff[0] == ';') {
            if (k) {
                mat[i][j] = atof(buffFloat);
                memset(buffFloat, 0, sizeof(buffFloat));
                k = 0;
            }
            j = 0;
            i++;
        }
        else {
            buffFloat[k++] = buff[0];
        }
    }
    fclose(file);
    for (i = 0; i < len; i++)
    {
        for (j = 0; j < len; j++)
        {
            if(pow(mat[i][j] - mat[j][i],2) > 0.0000001) return NULL;   
        }
    }
    
    return mat;
}

double dis(double **m, int i){
    if (i-1 < 0) return -1;
    return m[i][i-1];
}

double theta(double **m, int i){
    if (i-1 < 0) return -1000;
    return acos((pow(m[i][i+1],2) + pow(m[i][i-1],2) - pow(m[i+1][i-1],2))/(2*m[i][i+1]*m[i][i-1]));
}

double cosO(double **m, int i){
    //printf("%f %f\n", pow(m[i][i+1], 2) + pow(m[i+1][i+3],2) - 2*m[i][i+1]*m[i+1][i+3]*cos(theta(m, i+1))*cos(theta(m, i+2)) - pow(m[i][i+3],2), 2*m[i][i+1]*m[i+1][i+3]*sin(theta(m, i+1))*sin(theta(m, i+2)) );
    double theta2 = acos((pow(m[i+1][i+2],2) + pow(m[i+1][i+3],2) - pow(m[i+2][i+3],2))/(2*m[i+1][i+3]*m[i+1][i+2]));
    double res = (pow(m[i][i+1], 2) + pow(m[i+1][i+3],2) - 2*m[i][i+1]*m[i+1][i+3]*cos(theta(m, i+1))*cos(theta2) - pow(m[i][i+3],2))/
                    (2*m[i][i+1]*m[i+1][i+3]*sin(theta(m, i+1))*sin(theta2));
    if (res > 1) return 1;
    if (res < -1) return -1;
    return res;
}

double cosOmega(double **m, int i){
    if (i-3 < 0) return 1000;

    double Ai = (2*pow(m[i-2][i-1],2))*(pow(m[i-3][i-2],2) + pow(m[i-2][i] ,2) - pow(m[i-3][i] ,2));
    double Bi = pow(m[i-3][i-2],2) + pow(m[i-2][i-1],2) - pow(m[i-3][i-1],2);
    double Ci = pow(m[i-2][i-1],2) + pow(m[i-2][i],2) - pow(m[i-1][i],2);
    double Di = sqrt(
        4*pow(m[i-3][i-2],2)
        *pow(m[i-2][i-1],2)
        -pow(Bi,2) );
    double Ei = sqrt(4*pow(m[i-2][i-1],2)*pow(m[i-2][i],2)-pow(Ci,2) );
    double res = (Ai - Bi*Ci)/(Di*Ei);
    if (res > 1) return 1;
    if (res < -1) return -1;
    return res;
}

int isFeasible(double **mat, int len, vertex *v, int i){
    if (!v) return 0;
    int j = i-1;
    double pos[i][3];
    vertex *aux = v->p;
    while (aux)
    {
        //printf("\n\tj = %d {%f, %f, %f}\n", j, aux->pos.x, aux->pos.y, aux->pos.z);
        pos[j][0] = aux->pos.x;
        pos[j][1] = aux->pos.y;
        pos[j--][2] = aux->pos.z;
        aux = aux->p;
    }
    for (j = 0; j < i-3; j++)
    {
        if (mat[i][j] > 0)
        {           
            double dij = pow(v->pos.x - pos[j][0],2)+pow(v->pos.y - pos[j][1],2)+pow(v->pos.z - pos[j][2],2);
            double diff = pow(mat[i][j],2)- dij;
            if (diff < 0) diff = diff*(-1);
            double dist = sqrt(diff);//pow(diff,2);
            //printf("{%d,%d} {%f, %f, %f}-{%f, %f, %f} \n\t%.20f- %.20f = %.20f \n %.20f \n", i, j, v->pos.x, v->pos.y, v->pos.z, pos[j][0], pos[j][1], pos[j][2], mat[i][j], dij, diff, dist);
            if (dist > epsl)
            {
               return 0;
            }
        }
    }
    return 1;
}

double LDE(double **mat, double **m, int len){
    int i, j, count = 0;
    double dij = 0;
    for (i = 0; i < len; i++)
    {
        for (j = i+1; j < len; j++)
        {
            if (mat[i][j]>0){
                dij = dij + pow((pow(mat[i][j],2) - (pow(m[i][0] - m[j][0], 2) + pow(m[i][1] - m[j][1], 2) + pow(m[i][2] - m[j][2], 2))),2)/mat[i][j];
                count++;
            }
        }
    }
    return dij/count;
    
}

int qtd = 0;
int BranchAndPrune(double **mat, int len, vertex *v, int i){
    if(i < len){
        double thetai = theta(mat,i-1);
        double cti = cos(thetai);
        double sti = sin(thetai);
        double cwi = cosOmega(mat, i);
        //double cwi2 = cosO(mat, i-3);
        //printf("%f %f \n", cwi, cwi2 );
        //printf("%f %f \n", acos(cwi)*57.295779513, acos(cwi2)*57.295779513 );
        if (cwi*cwi > 1){
            printf("\n %f Error cwi", cwi);
            getchar();
        }  
        double swi = sqrt(1-pow(cwi, 2));
        double di = dis(mat,i);
        double bi1[] = {0-cti,0-sti,0,0-di*cti,
                        sti*cwi,0-cti*cwi,0-swi,di*sti*cwi,
                        sti*swi,0-cti*swi,cwi,di*sti*swi,
                        0,0,0,1};
        int j;
        //for (j = 0; j < 16; j++)
        //{
       //     printf("%f ", bi1[j]);
        //}
       // printf("\n");
        
        gsl_matrix_view Bi1 = gsl_matrix_view_array(bi1, 4, 4);

        swi = 0-swi;
        double bi2[] = {0-cti,0-sti,0,0-di*cti,
                        sti*cwi,0-cti*cwi,0-swi,di*sti*cwi,
                        sti*swi,0-cti*swi,cwi,di*sti*swi,
                        0,0,0,1};
        
        //for (j = 0; j < 16; j++)
        //{
        //    printf("%f ", bi2[j]);
        //}
        //printf("\n");

        gsl_matrix_view Bi2 = gsl_matrix_view_array(bi2, 4, 4);
        
        vertex *rv = (vertex*)calloc(1, sizeof(vertex));
        rv->p = v;
        rv->C = (gsl_matrix_view*)calloc(1, sizeof(gsl_matrix_view));
        *(rv->C) = gsl_matrix_view_array((double*)calloc(16,sizeof(double)), 4, 4);

        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans,
                   1.0, &(v->C->matrix), &Bi1.matrix,
                   0.0, &(rv->C->matrix));

        vertex *lv = (vertex*)calloc(1, sizeof(vertex));
        lv->p = v;
        lv->C = (gsl_matrix_view*)calloc(1, sizeof(gsl_matrix_view));
        *(lv->C) = gsl_matrix_view_array((double*)calloc(16,sizeof(double)), 4, 4);

        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans,
                   1.0, &(v->C->matrix), &Bi2.matrix,
                   0.0, &(lv->C->matrix));
        
        gsl_matrix_view *Xi = (gsl_matrix_view*)calloc(1, sizeof(gsl_matrix_view));
        *Xi = gsl_matrix_view_array((double*)calloc(4,sizeof(double)), 4, 1);

        double *y = (double*)calloc(4,sizeof(double));
        y[3] = 1;
        gsl_matrix_view *Y = (gsl_matrix_view*)calloc(1, sizeof(gsl_matrix_view));
        *Y = gsl_matrix_view_array(y, 4, 1);

        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans,
                    1.0, &(rv->C->matrix), &(Y->matrix),
                    0.0, &(Xi->matrix));
        rv->pos.x = Xi->matrix.data[0];
        rv->pos.y = Xi->matrix.data[1];
        rv->pos.z = Xi->matrix.data[2];
        
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans,
                    1.0, &(lv->C->matrix), &(Y->matrix),
                    0.0, &(Xi->matrix));
        lv->pos.x = Xi->matrix.data[0];
        lv->pos.y = Xi->matrix.data[1];
        lv->pos.z = Xi->matrix.data[2];
        free(y);
        free(Y);
        free(Xi->matrix.data);
        free(Xi);

        //printf("[ %f, %f", rv->pos.x, rv->pos.y);
        //printf("  %f]\n", rv->pos.z);
        //printf("[ %f, %f", lv->pos.x, lv->pos.y);
        //printf("  %f ]\n", lv->pos.z);


        //getchar();
            int s;
            //printf("\n");
        if (isFeasible(mat, len, rv, i)){
            v->r = rv;
            //for (s = 0; s < i; s++)
            //{
            //    printf(">");
            //}
            
            //printf("1\n");
            fflush(stdout);
            BranchAndPrune(mat, len, rv, i+1);
        }else {
            //for (s = 0; s < i; s++)
            //{
            //    printf("<");
            //}
            //
            //printf("2\n");
            free(rv->C->matrix.data);
            free(rv->C);
            free(rv);
            v->r = NULL;
        }
        
        if (isFeasible(mat, len, lv, i)){
            v->l = lv;

            //for (s = 0; s < i; s++)
            //{
            //    printf(">");
            //}
            //
            //printf("3\n");
            fflush(stdout);
            BranchAndPrune(mat, len, lv, i+1);
        }else {

            //for (s = 0; s < i; s++)
            //{
            //    printf("<");
            //}
            //
            //printf("4\n");
            free(lv->C->matrix.data);
            free(lv->C);
            free(lv);
            v->l = NULL;
        }
    }else{
        //printf(" ---------------- ue\n");
        //fflush(stdout);
        qtd++;
        double **mat2 = calloc(len, sizeof(double*));
        int i;
        for (i = 0; i < len; i++)
        {
            mat2[i] = calloc(3, sizeof(double));
        }
        i = len-1;
        printf("\n[");
        while (v)
        {
            mat2[i][0] = v->pos.x;
            mat2[i][1] = v->pos.y;
            mat2[i--][2] = v->pos.z;
            //printf("%f %f %f; ", v->pos.x, v->pos.y, v->pos.z);
            v = v->p;
        }

        double lde = LDE(mat, mat2, len);
        
        printf("\n \t LDE = %.20f\n", lde);

        for (i = 0; i < len; i++)
        {
            printf("%f %f %f; ", mat[i][0], mat[i][1], mat[i][2]);
        }
        
        printf("]\n");
        return 0;
    }
}


int printGrafo(int len, vertex *v, int i){
    int j = -1;
    if (v->l){
        printf("(%d, %d),", i, i+1);
        j = printGrafo(len, v->l, i+1);        
    }
    if (v->r){
        if (j != -1){
            printf("(%d, %d),", i, j);        
            j = printGrafo(len, v->r, j);        
        }else
        {
            printf("(%d, %d),", i, i+1);
            return printGrafo(len, v->r, i+1);        
        }
    }
    if (j == -1) return i+1;
    return j;
}

vertex *resul[10];

int solve(char *fileName, int len) {
    double **mat = matrixConstruct(fileName, len);
    if(mat == NULL){
        printf("Error when construct matrix\n");
        fflush(stdout);
        return 0;
    }

    double b1[] = {1,0,0,0,
                    0,1,0,0,
                    0,0,1,0,
                    0,0,0,1};
    gsl_matrix_view B1 = gsl_matrix_view_array(b1, 4, 4);
    
    double b2[] = {-1,0,0,0-dis(mat,1),
                    0,1,0,0,
                    0,0,-1,0,
                    0,0,0,1};
    gsl_matrix_view B2 = gsl_matrix_view_array(b2, 4, 4);

    double theta2 = theta(mat,1);
    //printf("[ cos %f\n", theta2);

    double ct2 = cos(theta2);
    double st2 = sin(theta2);
    double b3[] = {0-ct2,0-st2,0,0-dis(mat,2)*ct2,
                    st2,0-ct2,0,dis(mat,2)*st2,
                    0,0,1,0,
                    0,0,0,1};
    gsl_matrix_view B3 = gsl_matrix_view_array(b3, 4, 4);

    vertex *v = (vertex*)calloc(1, sizeof(vertex));
    v->C = &B1;
    //printf("[ %f, %f %f\n", v->pos.x, v->pos.y, v->pos.z);


    vertex *v2 = (vertex*)calloc(1, sizeof(vertex));
    v->r = v2;
    v2->C = &B2;
    v2->p = v;
    v2->pos.x = 0-dis(mat,1);

     //   printf("[ %f, %f %f\n", v2->pos.x, v2->pos.y, v2->pos.z);
    vertex *v3 = (vertex*)calloc(1, sizeof(vertex));
    v2->r = v3;
    v3->C = (gsl_matrix_view*)calloc(1, sizeof(gsl_matrix_view));
    *(v3->C) = gsl_matrix_view_array((double*)calloc(16,sizeof(double)), 4, 4);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans,
                   1.0, &B2.matrix, &(B3.matrix),
                   0.0, &(v3->C->matrix));
    v3->p = v2;
    v3->pos.x = dis(mat,2)*ct2 -dis(mat,1);
    v3->pos.y = dis(mat,2)*st2;

    //printf("[ %f, %f", v3->pos.x, v3->pos.y);
    //    printf("  %f]\n", v3->pos.z);

    BranchAndPrune(mat, len, v3, 3);
    printf("%d", qtd);
    //printf("[");
    //printGrafo(len, v, 0);
    //printf("]");
    return 1;
}

int countCols(char *fileName) {
    FILE *file = NULL;
    if ((file = fopen(fileName, "r")) == NULL) {
        /* error openning the file */
        perror("fopen: ");
        return 1;
    }

    int n, j = 0, k = 0;
    char buff[6];
    memset(buff, 0, sizeof(buff));
    if (fread(buff, sizeof(char), 5, file) < 5 || strcmp(buff, "E = [")) {
        fclose(file);
        perror("file out of params: ");
        return 1;
    }
    
    while ((n = fread(buff, sizeof(char), 1, file)) == 1) {
        if (buff[0] == ' ') {
            if (k > 0) {
                k = 0;
                j++;
            }
            continue;
        }
        else if (buff[0] == ';') {
            if (k) j++;
            break;
        }
        else {
            k++;
        }
    }
    fclose(file);
    return j;
}

int main(int argc, char **argv) {
    double time_spent = 0.0;

	clock_t begin = clock();

    if (argc < 3) {
        usage(argv[0]);
    }
    else {
        if (!strcmp(argv[1], "-s")) {
            int len;
            if (argc < 4) {
                len = countCols(argv[2]);
            }else {
                if ((len = atoi(argv[3])) == 0){
                    /* error openning the file */
                    perror("count_vertex param error: ");
                    return 1;
                }
            }
            
            solve(argv[2], len);
        } else if (!strcmp(argv[1], "-c")) {
            printf("%d\n",countCols(argv[2]));
        } else if (!strcmp(argv[1], "-show")) {
            printMat(matrixConstruct(argv[2], countCols(argv[2])), countCols(argv[2]));
        } else if (!strcmp(argv[1], "-d")) {
            ddd(argv[2],10);
        } else {
            usage(argv[0]);
        }
    }

    //getchar();
    
    clock_t end = clock();

	// calculate elapsed time by finding difference (end - begin) and
	// divide by CLOCKS_PER_SEC to convert to seconds
	time_spent += (double)(end - begin) / CLOCKS_PER_SEC;

	printf("\nTime elpased is %.20f seconds", time_spent);
    fflush(stdout);
}