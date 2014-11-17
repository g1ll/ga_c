/* 
 * File:   main.c
 * Author: Gill Velleda Gonzales
 *
 * Created on 2 de Novembro de 2013, 00:50
 */

#include <stdio.h>
#include <stdlib.h>
#include <tgmath.h>
#include "ga.h"

int main(int argc, char** argv) {

    //testeCadBits();//Representação do espaço de busca

    ga_fo_01(-1.8, 0.001, 2, 0, 200);

    return 0;

}

void seno_angulo() {

    float num;
    double rad, sn;
    printf("Digite o ângulo em graus: \n");
    scanf("%f", &num);
    printf("Voce digitou: %.1f graus", num);
    rad = num / 57.2957795;
    sn = sin(rad);
    printf("\nO seno do angulo %.1f é: %.5f", num, sn);

}

//Função Objetivo F(x) =

/**
 * @Função F(x) = x*sin(10*3,14*x)+1
 * @param x
 * @return y
 */
float fo_01(float x) {

    //printf("\nx = %.3f",x);
    return x * sinf(10 * M_PI * x) + 1;

}

/**
 * <h3>Algoritmos Genéticos</h3>
 * <p>A função ga_fo_01 tem como objetivo usar uma implementação simples de GA para
 * otimizar a função <i>F(x)=x*sin(10*pi*x)+1<i> em um espaço de busca definido.</p>
 * @file main.c
 * @version 0.0
 * @autor Gill Velleda Gonzales
 * @param a Início do espaço de busca 
 * @param step Valor de incremento
 * @param b Fim do espaço de busca
 * @param otimizacao 0 Minimiza a FO != 0 Maximiza
 * @pop_size tamanho da população (bug máximo 12 indivíduos)
 * @return  
 * Valor otimizado de <i><b>fo_01</b></i> no espaço definido por <b>a</b> até  <b>b</b> passo <b>step</b>
 * 
 * @ToDo
 * 
 *   
 */
long double ga_fo_01(double a, double step, double b, int otimizacao, int pop_size) {

    //Inicialização
    int **pop,
            **pop_mat, //População Intermediaria
            size_pop = pop_size, //População
            genes = 12, //Tamanho de Cromossomos (Cadeia de Bits)
            geracao = 1, i, j,
            count_ev = 0, //Contador de não-evoluções ou gerações iguais que não evoluirão;
            nev_max = 10, //Número máximo de populações iguais que não evoluirão, critério de parada
            *cr_opt;
    float *v_ap,
            taxa_cros = 0.9,
            taxa_mutate = 0.3;



    printf("\nAlgoritmo Genético Processando");
    //Espaço de Busca
    printf("\nEspaço de busca de %.3f a %.3f incrementado de %.3f", a, b, step);
    //Gerando a população
    printf("\nGerando a População ...\nCromossos: 12 genes binários\nPopulação: %d indivíduos", size_pop);
    pop = gerarPopIni(size_pop, genes, a, step, b);
    //LOOP
    do {
        if (geracao > 1)
            free(pop_mat); //Liberando memória - Descarta a população de Pais
        //Avalia População
        printf("\nAvaliando a população (Geração: %d)", geracao);

        v_ap = avaliarPop(pop, size_pop, genes, a, b);
        //Selecionar mais aptos
        printf("\n Selecionando mais aptos:\n\t-gerando mating-pool\n");
        pop_mat = selecaoTorneio(pop, v_ap, size_pop, genes, otimizacao);
        printf("\n População Original - Geracao %d:\n", geracao);
        printPop(pop, size_pop, genes);
        printf("\n\n População Intermediária:\n");
        printPop(pop_mat, size_pop, genes);
        free(pop); //Liberando memória - Descarta a população original
        //Operações de Crossover e Mutação
        pop = crossover(pop_mat, taxa_cros, size_pop, genes);
        printf("\n\n Nova Geração pós-cruzamentos:");
        printPop(pop, size_pop, genes);        
        //Mutação
        pop = mutate(pop, taxa_mutate, size_pop, genes);
        printf("\n\n Nova Geração pós-mutações:");
        printPop(pop, size_pop, genes);
        //Geração de nova população
        geracao++;
        //count_nev = final_ga(pop, v_ap, size_pop, genes, count_nev);
        count_ev = final_ga(pop_mat, v_ap, size_pop, genes, count_ev);
        printf("\t-----------CONTADOR %d----------\n", count_ev);
        printf("\n||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||\n");
        if (count_ev == nev_max) {
            /* aloca espaço para o cromossomo otimo */
            cr_opt = (int *) calloc(genes, sizeof (int));
            cr_opt = getOtimo(pop_mat, otimizacao, size_pop);
        }
    } while (count_ev < nev_max); //&& geracao <= 100);
    //FIM LOOP
    
     printf("\n GA- Finalizado:\n\tPopulação dos mais aptos:\n");
    printPop(pop_mat, size_pop, genes);
    printf("\n GA- Finalizado:\n\tx = %.3f Valor ótimo;\n\t y = %.3f", bitsToFloat(cr_opt), fo_01(bitsToFloat(cr_opt)));
    printf("\n Tamanho da População: %d", size_pop);
    printf("\n Genes por Cromossomos: %d", genes);
    printf("\n Taxa de Cruzamento: %.1f", taxa_cros);
    printf("\n Taxa de Mutação: %.1f", taxa_mutate);
    printf("\n Última Geração: %d:\n", geracao);
    printf("\n\n");
}


//----FUNÇÕES PARA TRABALHAR COM CADEIAS DE BITS

int *floatToBits(float num) {

    int i, n, rd, *bits;
    bits = (int*) malloc(12 * sizeof (int));
    for (i = 0; i <= 11; i++)
        bits[i] = 0;
    n = (int) round(num * 1000); //conversão para inteiro e arrdondamento com round()
    if (num < 0)
        n *= -1;
    //printf("\n %d %.3f %.3f",n,num,(num * 1000));
    i = 11;
    rd = n;
    do {
        bits[i] = rd % 2;
        rd = rd / 2;
        i--;
        if (rd == 0)
            bits[i] = 0;
    } while (rd != 0);

    //Retorna o complemento de dois para num negativo
    if (num < 0)
        return bitsToCompDois(bits);
    else
        return bits;
}

float bitsToFloat(int *bits) {

    int num = 0, j, i, n = 0;

    if (bits[0] == 1) {
        bits = compDoisToBit(bits);
        n = 1;
    }
    for (j = 11; j >= 0; j--) {
        if (bits[j] == 1) {
            num += pow(2, (11 - j));
        }
    }
    if (n == 1)
        return (num / 1000.0 * -1);
    else
        return (num / 1000.0);

}

int *bitsToCompDois(int *bits) {

    int *cbits, j, z;
    cbits = (int*) malloc(12 * sizeof (int));

    // Invertendo bits
    for (j = 0; j <= 11; j++) {
        if (bits[j] == 1)
            cbits[j] = 0;
        else
            cbits[j] = 1;
    }

    j = 11;
    do {
        if (cbits[j] == 1) {
            z = 1;
            cbits[j] = 0;
        } else {
            z = 0;
            cbits[j] = 1;
        }
        j--;
    } while (z == 1 && j >= 0);
    return cbits;
}

int *compDoisToBit(int *cbits) {

    return bitsToCompDois(cbits);

}

void testeCadBits() {
    //Teste função para converter de valor numérico para vetor binário
    // Todos os possíveis valores do espaço de busca convertidos para binário
    // Espaço de busca de -2.048: 0.001 : 2.048
    int *cr, i, *cd;
    float j, n;

    printf("VALORES DE X|\t\t\tCADEIA DE BITS (CROMOSSOMOS)\t\t| F(x)= x*sin(10*pi*x)+1");

    for (j = -2.048; j <= 2.048; j += 0.001) {
        cr = floatToBits(j);
        printf("\n x= %.3f | floatToBits(x): ", j);
        for (i = 0; i <= 11; i++) {
            printf("%d", cr[i]);
        }
        n = bitsToFloat(cr);
        printf(" >> bitToFloat(*bits): %.3f ", n);
        //aplicando a função objetivo
        printf("| y =  %.3f", fo_01(j));

    }
    printf("\n\n");
}
//-------------------------------------------------

//----FUNÇÕES COMPLEMENTARES AO ALGORITMO GENÉTICO

int **gerarPopIni(int size, int genes, double a, double step, double b) {
    int **pop, i, j;
    float r;

    /* aloca as linhas da matriz população */
    pop = (int **) calloc(size, sizeof (int *));
    if (pop == NULL) {
        printf("** Erro: Memoria Insuficiente **");
        return (NULL);
    }

    /* aloca as colunas da matriz população e gera os indipopíduos */
    for (i = 0; i < size; i++) {
        pop[i] = (int*) calloc(genes, sizeof (int));
        if (pop[i] == NULL) {
            printf("** Erro: Memoria Insuficiente **");
            return (NULL);
        }
    }

    printf("\n\t § gerando indivíduos aleatórios...\n");
    srand48(time(NULL));

    printf("\t § população:\n");

    for (i = 0; i < size; i++) {
        r = (float) mrand48() / 1e+09;
        if (r < a)
            r = a - r;
        if (r > b)
            r = r - b;
        printf("\t%.3f\t", r);
        pop[i] = floatToBits(r);
        for (j = 0; j <= genes; j++) {
            printf("%d", pop[i][j]);
        }
        printf("\n");
    }

    return pop;

}

float *avaliarPop(int **pop, int size, int genes, double a, double b) {
    int i;
    float *v;

    v = (float*) malloc(size * sizeof (float));

    //printf("\n Calculando Valores da FO:");
    for (i = 0; i < size; i++) {
        if (bitsToFloat(pop[i]) < a) {
            pop[i] = floatToBits(a - bitsToFloat(pop[i]));
        }
        if (bitsToFloat(pop[i]) > b) {
            pop[i] = floatToBits(bitsToFloat(pop[i]) - b);
        }
        v[i] = fo_01(bitsToFloat(pop[i]));
    }

    /*
        printf("\n Vetor de valores da fo:\n");

        for (i = 0; i < size; i++) {
             printf("\ny = %.3f", v[i]);
        }
     */
    return v;

}

int **selecaoTorneio(int **pop, float *v_ap, int size, int genes, int flag_o) {
    int **pop_int, i, j, cr_a, r = 0;
    float p_cr_a = 0.0;

    /* aloca as linhas da matriz população */
    pop_int = (int **) calloc(size, sizeof (int *));
    if (pop_int == NULL) {
        printf("** Erro: Memoria Insuficiente **");
        return (NULL);
    }

    /* aloca as colunas da matriz população e gera os indipopíduos */
    for (i = 0; i < size; i++) {
        pop_int[i] = (int*) calloc(genes, sizeof (int));
        if (pop_int[i] == NULL) {
            printf("** Erro: Memoria Insuficiente **");
            return (NULL);
        }
    }

    //printf("\n\t -inicializando mating-pool:\n");

    for (i = 0; i < size; i++) {
        pop_int[i] = 0;
        //printf("\n %d > %d", i, pop_int[i]);
    }

    //Seleção por torneio
    //printf("\n\t -selecionando por Torneio");
    srand(time(NULL));
    for (i = 0; i < size; i++) {
        // printf("\n grupo %d", i);
        j = 0;
        r = rand() % size;
        p_cr_a = v_ap[r];
        cr_a = r;
        //printf("\n %d | v_ap = %.3f | v_ap >> %.3f | x = %.3f", r, v_ap[r], p_cr_a, bitsToFloat(pop[r]));
        do {

            r = rand() % size;
            if (flag_o != 0) { //Flag de otimização
                if (v_ap[r] > p_cr_a) { //Maximização da FO
                    p_cr_a = v_ap[r];
                    cr_a = r;
                }
            } else {
                if (v_ap[r] < p_cr_a) {//Minimização da FO
                    p_cr_a = v_ap[r];
                    cr_a = r;
                }
            }
            // printf("\n %d | v_ap = %.3f | v_ap >> %.3f | x = %.3f", r, v_ap[r], p_cr_a, bitsToFloat(pop[r]));
            j++;
        } while (j < 2);

        pop_int[i] = pop[cr_a];
        //printf("\n >> %.3f", fo_01(bitsToFloat(pop[cr_a])));
        //printf("\n----");

    }

    return pop_int;

}

int **selecaoRoleta(int **pop, float *v_ap, int size, int genes, int flag_o) {
    int **pop_int, i, j, cr_a, r = 0;
    float p_cr_a = 0.0, apt_ac[size];

    /* aloca as linhas da matriz população */
    pop_int = (int **) calloc(size, sizeof (int *));
    if (pop_int == NULL) {
        printf("** Erro: Memoria Insuficiente **");
        return (NULL);
    }

    /* aloca as colunas da matriz população e gera os indipopíduos */
    for (i = 0; i < size; i++) {
        pop_int[i] = (int*) calloc(genes, sizeof (int));
        if (pop_int[i] == NULL) {
            printf("** Erro: Memoria Insuficiente **");
            return (NULL);
        }
    }

    //printf("\n\t -inicializando mating-pool:\n");

    for (i = 0; i < size; i++) {
        pop_int[i] = 0;
        //printf("\n %d > %d", i, pop_int[i]);
    }

    //Seleção por Roleta
    //printf("\n\t -selecionando por Roleta");
    srand(time(NULL));

    // - Calcular aptidões acumuladas
    // - procurar o 

    return pop_int;

}

int **crossover(int **pop_mat, float taxa, int size, int genes) {

    int **n_pop, i, r, v_aux[genes], j, c;


    printf("\n\n Processando o Cruzamento (CROSSOVER):\n\t -taxa: %.1f", taxa);

    /* aloca as linhas da matriz nova população */
    n_pop = (int **) calloc(size, sizeof (int *));
    if (n_pop == NULL) {
        printf("** Erro: Memoria Insuficiente **");
        return (NULL);
    }

    /* aloca as colunas da matriz nova população */
    for (i = 0; i < size; i++) {
        n_pop[i] = (int*) calloc(genes, sizeof (int));
        if (n_pop[i] == NULL) {
            printf("** Erro: Memoria Insuficiente **");
            return (NULL);
        }
    }

    //Número aleatório para probabilidade em relação a taxa de cruzamento

    srand(time(NULL));

    for (i = 0; i < size; i+=2) {

        if (size % 2 != 0 && i == size - 1) {
            i = size;
            break;
        }
        r = rand() % size;
        //printf("\n r1: %d ", r);
        //printf("\n i: %d ", i); 
        //printf("\n[%d] ", i);
        //printCr(pop_mat[i], genes);
        //printf("\n[%d] ", i + 1);
        //printCr(pop_mat[i + 1], genes);
        if (r <= taxa * size) {

            //printf("\ncruzando pares: %d - %d", i, i + 1);
            r = rand() % genes;
            //printf("\n r2: %d ", r);
            //printf("\n i+1: %d ", i+1);
            for (c = 0; c < 2; c++) {
                //printf("\n i+1: %d ", i+1);
                if (c == 0) {
                    for (j = 0; j < r; j++) {
                        //printf("\n :: i+1: %d ", i+1);
                        v_aux[j] = pop_mat[i + 1][j];
                    }
                    for (j = r; j < genes; j++) {
                        //printf("\n :: i: %d ", i);
                        v_aux[j] = pop_mat[i][j];
                    }

                } else {
                    for (j = 0; j < r; j++) {
                        //printf("\n :: i: %d ", i);
                        v_aux[j] = pop_mat[i][j];
                    }
                    for (j = r; j < genes; j++) {
                        //printf("\n :: i+1: %d ", i+1);
                        v_aux[j] = pop_mat[i + 1][j];
                    }
                }
                //printf("\n#######################");
                //printf("\n :: r: %d i: %d ", r,i);
                for (j = 0; j < genes; j++) {
                    n_pop[i + c][j] = v_aux[j];
                }
                //printf("\n %dº Filho: ", c + 1);
                //printCr(n_pop[i + c], genes);
            }
        } else {
            //printf("\ncai aqui\n");
            for (c = 0; c < 2; c++) {
                for (j = 0; j < genes; j++) {
                    n_pop[i + c][j] = pop_mat[i + c][j];
                }
                //printf("\n %d Filho: ", c + 1);
                //printCr(n_pop[i + c], genes);
            }

        }

    }

    return n_pop;
}

int **mutate(int **pop, float taxa, int size, int genes) {
    int **pop_mut, v_mut[genes], i, j, r1, r2;

    printf("\n\n Processando a mutação:\n\t -taxa: %.1f", taxa);

    /* aloca as linhas da matriz nova população */
    pop_mut = (int **) calloc(size, sizeof (int *));
    if (pop_mut == NULL) {
        printf("** Erro: Memoria Insuficiente **");
        return (NULL);
    }

    /* aloca as colunas da matriz nova população */
    for (i = 0; i < genes; i++) {
        pop_mut[i] = (int*) calloc(genes, sizeof (int));
        if (pop_mut[i] == NULL) {
            printf("** Erro: Memoria Insuficiente **");
            return (NULL);
        }
    }

    //Número aleatório para probabilidade em relação a taxa de mutação

    srand(time(NULL));


    for (i = 0; i < size; i++) {

        r1 = rand() % size;
        //printf("\nN p/ Taxa de Mutação %d", r1);
        if (r1 <= taxa * size) {
            //Aplicar Mutação
            //printf("\n----------------");
            //for (j = 0; j < r1; j++) {
            r2 = rand() % genes;
            //printf("\nAlterando bit [%d] %d > %d", r2, v_mut[r2],1-v_mut[r2]);
            if (pop[i][r2] == 0)
                pop[i][r2] = 1;
            else
                pop[i][r2] = 0;
            //}
            //printf("\n----------------");
            //printf("\nCromossomo Original [%d]", i);
            //printCr(pop[i], genes);
            //printf("\nCromossomo Mutado: ", i);
            //printCr(v_mut, genes);
            //printf("\n_________________");
        }
        pop_mut[i] = pop[i];
    }

    return pop_mut;


}

int final_ga(int **pop, float *v_ap, int size, int genes, int count_ev) {
    int i,ci, d;
    float ap, ap2;

    //Compara se os individuos pais são iguais entre si
    for (i = 0, ci = 0; i < size; i++) {
        ap = fo_01(bitsToFloat(pop[i]));
        d = i - 1;
        if (d < 0)
            d = 0;
        ap2 = fo_01(bitsToFloat(pop[d]));
        if (ap == ap2) {
            ci++;
        }

    }

    if (ci >= size * 0.95) {//ci >= size * 0.8 &&
        count_ev++;
    }
    printf("\n\t |ci = %d | contador ev = %d\n", ci, count_ev);
    /*
        if(count_ev==count_max)
            stop_ga=1;
    
        return stop_ga;
     */
    return count_ev;
}

int* getOtimo(int **pop, int flag_o, int size) {
    int i, otimo = 0;
    float a,b;
    for (i = 0; i < size; i++) {
        a = fo_01(bitsToFloat(pop[i]));
        b = fo_01(bitsToFloat(pop[otimo]));
        if (flag_o == 0) {            
            if (a < b )
                otimo = i;
        } else {
            if (a > b)
                otimo = i;
        }
    }
    return pop[otimo];
}

void printPop(int **pop, int size, int genes) {

    int i, j;

    for (i = 0; i < size; i++) {

        printf("\n [%d] ", i);
        for (j = 0; j < genes; j++) {
            printf("%d", pop[i][j]);
        }

        printf(" | x = %.3f : f(x) = %.3f", bitsToFloat(pop[i]), fo_01(bitsToFloat(pop[i])));
    }

}

void printCr(int *cromossomo, int genes) {
    int i;
    for (i = 0; i < genes; i++) {
        printf("%d", cromossomo[i]);
    }
}