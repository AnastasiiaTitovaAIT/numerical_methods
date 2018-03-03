#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <math.h>
double h, alf;
double f (double x, double y) //правая часть уравнения y' = f(x,y)
{
    return -y-pow(x, 2);
}
double f1 (double x, double u, double v) //первое уравнение системы
{
    return sin(1.4 * pow(u,2)) - x + v;
}
double f2 (double x, double u, double v) //первое уравнение системы
{
    return x + u - 2.2 * pow(v,2) + 1;
}
double sol (double x) //точное решение уравнения первого порядка
{
    return -pow(x,2) + 2*x - 2 + 12 * exp(-x);
}
int main(void)
{
    double l;
    int type, am;
    FILE *fp = fopen("t1.txt", "w");
    printf("Введите 0, чтобы решить систему уравнений, 1 - уравнение первого порядка.\n");
    scanf("%d", &am);
    printf("Выберите порядок точности - 2 или 4.\n");
    scanf("%d", &type);
    printf("Введите длину отрезка:\n");
    scanf("%lf", &l);
    printf("Введите величину шага:\n");
    scanf("%lf", &h);
    int  n = (int)l/h;
    fprintf(fp, "Количество итераций: %d\n", n);
    double y[n + 1], x[n + 1];
    if(!am) { //решаем систему
        double y2[n + 1];
        x[0] = 0; //начальные условия
        y[0] = 1;
        y2[0] = 0.5;
        if(type == 2){
            printf("Введите параметр альфа:\n");
            scanf("%lf", &alf);
        }
        fprintf(fp, "Точка\t Значение функции 1 Значение функции 2\n");
        for (int i = 0; i < n; i++) {
            x[i + 1] = x[0] + (i + 1)*h;
            if (type == 2) {
                second_sys(y, x, y2, i); //метод Рунге-Кутта второго порядка для системы уравнений
            } else {
                forth_sys(y, x, y2, i); //метод Рунге-Кутта четвертого порядка для системы уравнений
            }
            fprintf(fp, "x_%d = %lf\ty_%d = %lf\ty2_%d = %lf\n", i + 1, x[i + 1], i + 1, y[i + 1], i + 1, y2[i + 1]);
        }

    } else {
        x[0] = 0; //начальные условия
        y[0] = 10;
        if(type == 2){
            printf("Введите параметр альфа:\n");
            scanf("%lf", &alf);
        }
        fprintf(fp,"Точка\t\tЗначение\tЗначение точного решения\n");
        for (int i = 0; i < n; i++) {
            x[i + 1] = x[0] + (i + 1)*h;
            if (type == 2) {
                second(y, x, i); //метод Рунге-Кутта второго порядка
            } else {
                forth(y, x, i); //метод Рунге-Кутта четвертого порядка
            }
            fprintf(fp,"x_%d = %lf\ty_%d = %lf\tu_%d = %lf\n", i + 1, x[i + 1], i + 1, y[i + 1], i + 1, sol(x[i + 1]));
        }
    }
}

void second(double *y, double *x, int i )//метод Рунге-Кутта второго порядка для дифференциального уравнения первого порядка
{
    y[i + 1] = y[i] + ((1 - alf) * f(x[i], y[i]) + alf * f(x[i] + h/(2*alf), y[i] + (h/(2*alf))*f(x[i], y[i]))) * h;
}

void forth(double *y, double *x, int i) //метод Рунге-Кутта четвертого порядка для дифференциального уравнения первого порядка
{
    double k1, k2, k3, k4;
    k1 = h * f(x[i], y[i]);
    k2 = h * f(x[i] + h/2, y[i] + k1/2);
    k3 = h * f(x[i] + h/2, y[i] + k2/2);
    k4 = h * f(x[i] + h, y[i] + k3);
    y[i + 1] = y[i] + (k1 + 2*k2 + 2*k3 + k4)/6;
}
void second_sys(double *y, double *x, double *y2, int i) //метод Рунге-Кутта второго порядка для системы уравнений
{
    y[i + 1] = y[i] + ((1 - alf) * f1(x[i], y[i], y2[i]) + alf * f1(x[i] + h/(2*alf), y[i] + (h/(2*alf))*f1(x[i], y[i], y2[i]),
y2[i] + (h/(2*alf))*f2(x[i], y[i], y[2]))) * h;
    y2[i + 1] = y2[i] + ((1 - alf) * f2(x[i], y[i], y2[i]) + alf * f2(x[i] + h/(2*alf), y[i] + (h/(2*alf))*f1(x[i], y[i], y2[i]),
y2[i] + (h/(2*alf))*f2(x[i], y[i], y[2]))) * h;
}
void forth_sys(double *y, double *x, double *y2, int i) //метод Рунге-Кутта четвертого порядка для системы уравнений
{
    double k1, k2, k3, k4, l1, l2, l3, l4;
    k1 = f1(x[i], y[i], y2[i]);
    l1 = f2(x[i], y[i], y2[i]);
    k2 = f1(x[i] + h / 2, y[i] + h / 2 * k1, y2[i] + h / 2 * l1);
    l2 = f2(x[i] + h / 2, y[i] + h / 2 * k1, y2[i] + h / 2 * l1);
    k3 = f1(x[i] + h / 2, y[i] + h / 2 * k2, y2[i] + h / 2 * l2);
    l3 = f2(x[i] + h / 2, y[i] + h / 2 * k2, y2[i] + h / 2 * l2);
    k4 = f1(x[i] + h, y[i] + h * k3, y2[i] + l3);
    l4 = f2(x[i] + h, y[i] + h * k3, y2[i] + l3);
    y[i + 1] = y[i] + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4);
    y2[i + 1] = y2[i] + h / 6 * (l1 + 2 * l2 + 2 * l3 + l4);
}
