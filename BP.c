#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_blas.h>

double epsl = 2.561;

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
}

void teste() {
    double a[] = {0.11, 0.12, 0.13, 0.21, 0.22, 0.23};
    double b[] = {1011, 1012, 1021, 1022, 1031, 1032};
    double c[] = {0.00, 0.00, 0.00, 0.00};

    gsl_matrix_view A = gsl_matrix_view_array(a, 2, 3);
    gsl_matrix_view B = gsl_matrix_view_array(b, 3, 2);
    gsl_matrix_view C = gsl_matrix_view_array(c, 2, 2);

    /* Compute C = A B */

    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans,
                   1.0, &A.matrix, &B.matrix,
                   0.0, &C.matrix);

    printf("[ %g, %g\n", c[0], c[1]);
    printf("  %g, %g ]\n", c[2], c[3]);
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
        else if (!k && buff[0] == '0') {
            j++;
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
    if (i-2 < 0) return -1;
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
    if (i-3 < 0) return -1;
    /* i--;
    printf("%f, %d \n", acos((pow(m[i][i+1],2) + pow(m[i][i-1],2) - pow(m[i+1][i-1],2))/(2*m[i][i+1]*m[i][i-1]))+ acos((pow(m[i][i-1],2) + pow(m[i-1][i+1],2) - pow(m[i][i+1],2))/(2*m[i][i-1]*m[i-1][i+1]))+ acos((pow(m[i-1][i+1],2) + pow(m[i][i+1],2) - pow(m[i][i-1],2))/(2*m[i][i+1]*m[i-1][i+1])), i);
    i--;
    printf("%f, %d \n", acos((pow(m[i][i+1],2) + pow(m[i][i-1],2) - pow(m[i+1][i-1],2))/(2*m[i][i+1]*m[i][i-1]))+ acos((pow(m[i][i-1],2) + pow(m[i-1][i+1],2) - pow(m[i][i+1],2))/(2*m[i][i-1]*m[i-1][i+1]))+ acos((pow(m[i-1][i+1],2) + pow(m[i][i+1],2) - pow(m[i][i-1],2))/(2*m[i][i+1]*m[i-1][i+1])), i);
    i++;i++;
    printf("%f, %f %f %f\n", pow(m[i-3][i-2],2) + 
        pow(m[i][i-2],2) - 
        (
            2*
            m[i-3][i-2]*
            m[i][i-2]*
            cos(theta(m,i-2))*
            cos(theta(m,i-1))
        ) - 
        pow(m[i][i-3],2), 2*
        m[i-3][i-2]*
        m[i][i-2]*
        sin(theta(m,i-2))*
        sin(theta(m,i-1)),  m[i-3][i-2], theta(m, i));
    printf("%f, %f %f %f\n", cos(theta(m,i-2)), cos(theta(m,i-1)), m[i-3][i-2], cos(theta(m, i)));
    
    return 
    (
        pow(m[i-3][i-2],2) + 
        pow(m[i][i-2],2) - 
        (
            2*
            m[i-3][i-2]*
            m[i][i-2]*
            cos(theta(m,i-2))*
            cos(theta(m,i-1))
        ) - 
        pow(m[i][i-3],2)
    )/
    (
        2*
        m[i-3][i-2]*
        m[i][i-2]*
        sin(theta(m,i-2))*
        sin(theta(m,i-1))
    );
    */
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
    int j;
    //int s = 0;
    for (j = 0; j < i; j++) {
        
        if (mat[i][j] && (i-j > 3)) {
        //s++;
            int k = i;
            vertex *aux = v;
            while (k-- > j){
                aux = aux->p;
            }
                printf("{%d,%d, %d} %.20f %f, (%f, %f) %f %f %f\t %f %f %f\n",i,j, k, 
                pow(sqrt(pow(v->pos.x - aux->pos.x, 2) + pow(v->pos.y - aux->pos.y, 2) + pow(v->pos.z - aux->pos.z, 2)) - mat[i][j], 2)
                ,epsl
                , mat[i][j], sqrt(
                            pow(v->pos.x - aux->pos.x, 2) + 
                            pow(v->pos.y - aux->pos.y, 2) + 
                            pow(v->pos.z - aux->pos.z, 2)
                        ),aux->pos.x, aux->pos.y, aux->pos.z, v->pos.x, v->pos.y, v->pos.z );
            if (pow(sqrt(pow(v->pos.x - aux->pos.x, 2) + pow(v->pos.y - aux->pos.y, 2) + pow(v->pos.z - aux->pos.z, 2)) - mat[i][j], 2) > epsl){
                return 0;
           }
        
        }
    }   
    //printf("%d", s); 
    return 1;
}
int qtd = 0;
int BranchAndPrune(double **mat, int len, vertex *v, int i){
    if(i < len){
        double thetai = theta(mat,i-1);
        double cti = cos(thetai);
        double sti = sin(thetai);
        double cwi2 = cosOmega(mat, i);
        double cwi = cosO(mat, i-3);
        printf("%f %f \n", cwi, cwi2 );
        printf("%f %f \n", acos(cwi)*57.295779513, acos(cwi2)*57.295779513 );
        if (cwi*cwi > 1){
            printf("\n %f CORREEEEEEEEE BEEEEERG!!!!!!!!!!!!!!!!!!!!!!!!!!", cwi);
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

        printf("[ %f, %f", rv->pos.x, rv->pos.y);
        printf("  %f]\n", rv->pos.z);
        printf("[ %f, %f", lv->pos.x, lv->pos.y);
        printf("  %f ]\n", lv->pos.z);


        getchar();
            int s;
            printf("\n");
        if (isFeasible(mat, len, rv, i)){
            v->r = rv;
            for (s = 0; s < i; s++)
            {
                printf(">");
            }
            
            printf("1\n");
            fflush(stdout);
            BranchAndPrune(mat, len, rv, i+1);
        }else {
            for (s = 0; s < i; s++)
            {
                printf("<");
            }
            
            printf("2\n");
            free(rv->C->matrix.data);
            free(rv->C);
            free(rv);
            v->r = NULL;
        }
        
        if (isFeasible(mat, len, lv, i)){
            v->l = lv;

            for (s = 0; s < i; s++)
            {
                printf(">");
            }
            
            printf("3\n");
            fflush(stdout);
            BranchAndPrune(mat, len, lv, i+1);
        }else {

            for (s = 0; s < i; s++)
            {
                printf("<");
            }
            
            printf("4\n");
            free(lv->C->matrix.data);
            free(lv->C);
            free(lv);
            v->l = NULL;
        }
    }else{
        printf(" ---------------- ue\n");
        fflush(stdout);
        qtd++;
        return 0;
    }
}

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
    double ct2 = cos(theta2);
    double st2 = sin(theta2);
    double b3[] = {0-ct2,0-st2,0,0-dis(mat,2)*ct2,
                    st2,0-ct2,0,dis(mat,2)*st2,
                    0,0,1,0,
                    0,0,0,1};
    gsl_matrix_view B3 = gsl_matrix_view_array(b3, 4, 4);

    vertex *v = (vertex*)calloc(1, sizeof(vertex));
    v->C = &B1;


    vertex *v2 = (vertex*)calloc(1, sizeof(vertex));
    v->r = v2;
    v2->C = &B2;
    v2->p = v;
    v2->pos.x = 0-dis(mat,1);

        printf("[ %f, %f", v2->pos.x, v2->pos.y);
        printf("  %f ]\n", v2->pos.z);
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

    printf("[ %f, %f", v3->pos.x, v3->pos.y);
        printf("  %f]\n", v3->pos.z);

    BranchAndPrune(mat, len, v3, 3);
    printf("%d", qtd);
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
        else if (!k && buff[0] == '0') {
            j++;
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
            
            return solve(argv[2], len);
        } else if (!strcmp(argv[1], "-c")) {
            printf("%d\n",countCols(argv[2]));
        } else {
            usage(argv[0]);
        }
    }
}