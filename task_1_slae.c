#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include <linux/limits.h>
#include <fcntl.h>
void gaus(int mode); //метод Гаусса или модифицированный метод Гаусса (в зависимости от параметра mode)
void chng(void); //функция, меняющая текущую строку с нулевым элементом на диагонали, на любую другую с ненулевым
void chngdiag(void); //функция выбора главного элемента по строке
void inputfile(void); //функция ввода матрицы из файла
void inputform(void);//функция ввода матрицы поформуле
void det(void); //функция нахождения определителя матрицы
void inverse(void); //функция нахождения обратной матрицы
int size;//размер матрицы А
double **mas, **inv;//матрицы А и обратная к ней
double *b, *c, *x, det_a;//вектор свободных коэффциентов, вектор решений, определитель мтарицы А

int main(void)
{
    int type;
    printf("Выберите способ ввода матрицы: 0 - из файла, 1 - формулой\n");
    scanf("%d", &type);
    switch (type)
    {
    case 0:
        inputfile();//формат файла первая строка - n. далее n строк, содержащие n элементов матрицы и соотв.элемент столбца значений
        break;
    case 1:
        inputform();
        break;
    default:
        break;
    }
    printf("Выберите метод поиска решений: 0 - метод Гаусса, 1 - модифицированный метод Гаусса.\n");
    scanf("%d", &type);
    x = calloc(size, sizeof(*x)); //выделяем память под столбец решений и обнуляем его
    for (int i = 0; i < size; i++) {
        x[i] = 0;
    }
    inverse(); //нахождение обратной матрицы
    gaus(type); //метод Гаусса
    det();//определитель
    FILE *fp = fopen("out.txt", "w");
    fprintf(fp, "\nCтолбец решений системы Ax = f :\n");
    for(int i = 0; i < size; i++) {
        fprintf(fp,"%.10lf\n", x[i]);
    }
    fprintf(fp,"\nОпределитель матрицы А: \n det(A) = %.10lf\n", det_a);
    fprintf(fp,"\nОбратная матрица:\n");
    for (int i = 0; i < size; i++) {
        for(int  j = 0; j < size; j++) {
            fprintf(fp,"%.10lf ", inv[i][j]);
        }
        fprintf(fp,"\n");
    }
//освобождаем выделенную память
    for (int i = 0; i < size; i++) {
        free(mas[i]);
    }
    free(mas);
    free(inv);
    free(b);
    free(x);
    return 0;
}
void inputfile(void) //ввод матрицы из файла
{
    char name[1024];
    char str[PATH_MAX];
    printf("Введите имя файла\n");
    fscanf(stdin,"%s", name);
    FILE *fp = fopen(name, "r");
    fscanf(fp, "%d", &size);
    mas = calloc(size, sizeof(double*)); //выделяем память под матрицу и столбец свободных коэффициентов
    b = calloc(size, sizeof(*b));
    for (int i = 0; i < size; i++) {
        mas[i] = calloc(size, sizeof(**mas));
    }
    for (int i = 0; i < size; i++) { //считываем матрицу и столбец свободных коэффициентов
        for(int j = 0; j < size; j++) {
            fscanf(fp, "%lf", &mas[i][j]);
        }
        fscanf(fp, "%lf", &b[i]);
    }
    for (int i = 0; i < size; i++) {
        for(int  j = 0; j < size; j++)
        {
            printf("%lf ", mas[i][j]);
        }
        printf("\n");
    }
    fclose(fp);
    return;
}

void inputform(void)
{
    double q = 1.001 - 2 * 4 * 0.001;
    printf("Введите N - порядок матрицы.\n");
    scanf("%d", &size);
    mas = calloc(size, sizeof(double*)); //выделяем память под матрицу и столбец свободных коэффициентов
    b = calloc(size, sizeof(*b));
    for (int i = 0; i < size; i++) {
        mas[i] = calloc(size, sizeof(**mas));
    }
    for (int i = 0; i < size; i++) { //генерируем матрицу и столбец свободных коэффициентов в соотв. с заданной формулой
        for(int j = 0; j < size; j++) {
            if (i != j) {
                mas[i][j] = pow(q, i + j) + 0.1 * (j - i);
            } else {
                mas[i][j] = pow(q - 1, i + j);
            }
        }
        b[i] = size * exp(-2.0/i) * cos(-2.0);
    }
    printf("\nСформированная матрица и столбец свободных членов:\n");
    for (int i = 0; i < size; i++) {
        for(int  j = 0; j < size; j++) {
            printf("%lf ", mas[i][j]);
        }
        printf("%lf ", b[i]);
        printf("\n");
    }
    return ;
}

void gaus(int mode) //метод Гаусса
{
    for(int k = 0; k < size; k++) { //прямой ход
       for(int i = k + 1; i < size; i++){
            if (mode) {
                chngdiag(); //выбор главного элемента, если был выбран модифицированный мтеод Гаусса
            } else {
                chng(); //проверка на нулевой ведущий элемент
            }
            double coef = mas[i][k] / mas[k][k];
            for(int j = k; j < size; j++) {
                mas[i][j] -= coef * mas[k][j];
            }
            b[i] -= coef * b[k];
       }
    }
    for(int l = size - 1; l >= 0; l--) { //обратный ход
        double s = 0;
        for (int j = l; j < size; j++) {
            s = s + mas[l][j] * x[j];
        }
        x[l] = (b[l] - s) / mas[l][l];
    }
    return;
}

void chng(void) //если диагональный элемент равен нулю, поменяем местами строки
{
   static int i = 0;
   double tmp = 0;
   for (i = 0; i < size; i++) {
       if (mas[i][i] == DBL_EPSILON) {
           for (int j = 0; j < size; j++) {
               if (j == i) {
                    continue;
               }
               if (mas[j][i] != 0 && mas[i][j] != 0){
                   for(int k = 0; k < size; k++) {
                       tmp = mas[j][k];
                       mas[j][k] = mas[i][k];
                       mas[i][k] = tmp;
                   }
                   tmp = b[j];
                   b[j] = b[i];
                   b[i] = tmp;
                   break;
               }
           }
       }
   }
   return;
}

void chngdiag(void) //выбор главного элемента
{
    static int ind = 0; //эта переменная хранит количество вхождений в функцию, т.е. номер текущей строки
    double max = DBL_MIN;
    int id_m, j;
    double temp;
    for (j = ind; j < size; j++) {
        if (fabs(mas[ind][j]) - max > DBL_EPSILON) {
            max = fabs(mas[ind][j]);
            id_m = j;
        }
    }
    if (j == ind) {
        return;
    }
    for(int k = 0; k < size; k++) {
        temp = mas[k][id_m];
        mas[k][id_m] = mas[k][ind];
        mas[k][ind] = temp;
    }
    ind++;
    return;
}

void det(void) //определитель матрицы А
{
    det_a = 1;
    for (int i = 0; i < size; i++) {
        det_a *= mas[i][i];
    }
    return;
}

void inverse(void) //поиск обратной матрицы
{
    double **E = calloc(size, sizeof(double*));
    double temp;
    inv = calloc(size, sizeof(double*));
    for (int i = 0; i < size; i++) {   //выделяем память под единичную и обратную матрицы
        inv[i] = calloc(size, sizeof(**inv));
        E[i] = calloc(size, sizeof(**E));
    }
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            inv[i][j] = mas[i][j];
            if (i == j) {
                E[i][j] = 1.0;
            } else {
                E[i][j] = 0.0;
            }
        }
    }
    for (int k = 0; k < size; k++) {
        temp = inv[k][k];
        for (int j = 0; j < size; j++) {
            inv[k][j] /= temp;
            E[k][j] /= temp;
        }
        for (int i = k + 1; i < size; i++) {
            temp = inv[i][k];
            for (int j = 0; j < size; j++) {
                inv[i][j] -= inv[k][j] * temp;
                E[i][j] -= E[k][j] * temp;
            }
        }
    }
    for (int k = size - 1; k > 0; k--) {
        for (int i = k - 1; i >= 0; i--) {
            temp = inv[i][k];
            for (int j = 0; j < size; j++) {
                inv[i][j] -= inv[k][j] * temp;
                E[i][j] -= E[k][j] * temp;
            }
        }
    }
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            inv[i][j] = E[i][j];
        }
    }
    for (int i = 0; i < size; i++) {
        free(E[i]);
    }
    free(E);
    return;
}

void check(void)
{
    double temp[size];
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++)
        {
            temp[i] += mas[i][j] * x[j];
        }
        temp[i] -= b[i];
    }
}
double getpsi(double *psi)
{
    double ret = 0;
    int j;
    for(int i = 0; i < size; i++)
        for(int j = 0; j < size; j++) {
            psi[i] += mas[i][j] * x[j];
        }
    for (int i = 0; i < size; i++) {
        psi[i] -= b[i];
        ret += pow(psi[i], 2);
    }
    ret = sqrt(ret);
    printf("\n\n%lf\n", ret);
    return ret;
}
