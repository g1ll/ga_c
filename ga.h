/* 
 * File:   ga.h
 * Author: echoes
 *
 * Created on 2 de Novembro de 2013, 17:19
 */

#ifndef GA_H
#define	GA_H

#ifdef	__cplusplus
extern "C" {
#endif




#ifdef	__cplusplus
}

#endif

//Função simples de teste
extern void seno_angulo();

//Função Objetivo
extern float fo_01(float x);

//Função de Teste Cadeia de bits-CROMOSSOMOS
extern void testeCadBits(void);

//Algoritmo Genético p/ Otimizar fo_01
extern long double ga_fo_01(double a, double step, double b,int otimizacao, int pop_size);

//Função para geração da população
extern int **gerarPopIni(int size,int genes, double a, double step, double b);

//Função para avaliar a população e definir aptidões
extern float *avaliarPop(int **pop, int size, int genes,double a, double b);

//Função selecaoTorneio() responsável por selecionar os mais aptos através do método de torneio 
//retornar a pop intermediaria
extern int **selecaoTorneio(int **pop, float *v_ap,int size,int genes, int flag_o);

//Função selecaoRoleta() responsável por selecionar os mais aptos através do método de roleta
extern int **selecaoRoleta(int **pop, float *v_ap,int size,int genes, int flag_o);

//Função para crossover-cruzamento.
extern int **crossover(int **pop_mat, float taxa,int size,int genes);

//Funçõa para mutação
extern int **mutate(int **pop, float taxa,int size,int genes );

//Função finalga() critério de parada, retorna 1 para parar e 0 para continuar
extern int final_ga(int **pop,float *v_ap,int size,int genes,int count_ev);

extern int* getOtimo(int **pop,int flag_o,int size);

//Função para imprimir cromossomos
void printCr(int *cromossomo,int genes);

//Função printPop()
extern void printPop(int **pop,int size, int genes);

//Função para converção de cadeias de bits em floats de 0,000
extern int *floatToBits(float num);

extern float bitsToFloat(int *bits);

//Funções para números negativos, complemento de dois
extern int *bitsToCompDois(int *bits);

extern int *compDoisToBit(int *cbits);


#endif	/* GA_H */

