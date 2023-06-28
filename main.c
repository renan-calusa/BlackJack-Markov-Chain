#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

long double** init();
void set_default(long double** matriz);
void estimator(int* amostra, long double** matriz);
long double likelihood(int* v, float p, int teta, long double** matriz);
long double probability_function (int v, float p, int teta, long double** matriz);
void transition(long double** matriz, int state, float p, int teta);
void printMatriz(long double** matriz);
double generateRandom();

int sample_size;



int main (int argc, char** argv) {

    float p; // entre 0 e 1
    int teta; // entre 1 e 21
    int* amostra;
    
    long double** matriz = init();
    
    sample_size = atoi(argv[1]);
    
    if (argc == sample_size + 4) {
            
            // Caso que se quer calcular o likelihood de um agente descrito
            
            amostra = (int*) malloc(sizeof(int)*sample_size);
            for (int i=0; i < sample_size; i++) amostra[i] = atoi(argv[i+2]);  // entre 0 e 30
            
            p = strtof(argv[argc-2], NULL);
            teta = atoi(argv[argc-1]);

            if (p <= 0 || p > 1 || teta < 1 || teta > 21) {
                printf("[-] Valores passados fora do intervalo esperado!");
                return 0;
            }
    
            srand(time(NULL));
    
            printf("likelihood(%.3f, %i; amostra) = %Lf\n", p, teta, likelihood(amostra, p, teta, matriz));
            //printf("likelihood(%.3f, %i; {%i", p, teta, amostra[0]);
            //for (int i=1; i < sample_size; i++) printf(", %i", amostra[i]);
            //printf("}) = ");
            //printf("%Lf\n", likelihood(amostra, p, teta, matriz));
    }

    else if (argc == sample_size + 2) {
    
    	// Caso que se quer estimar 'teta' e 'p' de uma amostra
    	
    	amostra = (int*) malloc(sizeof(int)*sample_size);
        for (int i=0; i < sample_size; i++) amostra[i] = atoi(argv[i+2]);  // entre 0 e 30
        
        estimator(amostra, matriz);
    }

    else {
        printf("[-] Usage: ./exec <tamanho da amostra> <amostra separada por 'space bar'> <float p> <int teta>\n");
        return 0;
    }
    

    return 0;
}


// Aloca espaço e coloca valores default
long double** init() {

	long double** matriz = (long double**) malloc(sizeof(long double*)*22);     // No máximo 22 estados (0 à 21)
	for (int i=0; i < 22; i++) matriz[i] = (long double*) malloc(sizeof(long double)*31);    // No máximo 31 possibilidades de saída (0 à 30)

	return matriz;
}


// Seta valores default
void set_default(long double** matriz) {

	for (int i = 0; i < 22; i++) for (int j = 0; j < 31; j++) matriz[i][j] = logl(0);   // equivalente à Pr(0)

	matriz[0][0] = logl(1); // (log(1) = 0) 100% de estar com zero no estado zero
}


// Estimar 'teta' e 'p' por BruteForce
void estimator(int* amostra, long double** matriz) {

	int teta;
	float p;
	
	int max_teta = 0;
	float max_p = 0;
	long double max_likelihood = -INFINITY;
	
	for (teta = 1; teta <= 21; teta++) {
	
		for (p = 0.05; p <= 1; p += 0.05) {
		
			long double curr = likelihood(amostra, p, teta, matriz);
			
			if (curr > max_likelihood) {
			
				max_likelihood = curr;
				max_p = p;
				max_teta = teta;
			}
		}
	}
	
	printf("[!]teta estimator: %i\t[!]p estimator: %.3f\n", max_teta, max_p);
}


// Somatorio do log das probabilidades -> Sum i of (log Pr(vi|teta, p)
long double likelihood (int* v, float p, int teta, long double** matriz) {

	long double likelihood = logl(1);

	for (int i=0; i < sample_size; i++) likelihood += logl(probability_function(v[i], p, teta, matriz));
    
	return likelihood;
}


long double probability_function (int v, float p, int teta, long double** matriz) {
	
	set_default(matriz);
	
	// Realizar transicoes de markov até o jogador quiser parar no valor "teta" com probabilidade "p"
	for (int i=0; i < 21; i++) transition(matriz, i, p, teta);
	
	printMatriz(matriz);
	
	// Pega a probabilidade de terminar com um valor v num jogo de BlackJack - em log()
	long double res;
	for (int i=0; i < 22; i++) res += expl(matriz[i][v]);

	printf("Pr(%i; %.3f, %i) = %Lf\n", v, p, teta, res);

	return res;
}


void transition(long double** matriz, int state, float p, int teta) {
    
	// Não existe um estado maior que 21
	if (state >= 21) return;
		
	long double* current = matriz[state];
	long double* next = matriz[state+1];
	
	// Percorrendo todos valores no atual (i) e no proximo (j)
	for (int i = 0; i < 31; i++) {
		for (int j = 0; j < 31; j++) {
            
			if (current[i] != logl(0)) {
				
				long double probability;
				long double stop = logl(0);
				int dif = j - i;
				
				// Se pode-se alcançar valores >= 'teta', verifica se queremos parar com probabilidade 'p'
				if (i >= teta) {
                
					// Se mantem com a mesma soma caso p de certo
					if (generateRandom() < p) {
					
						next[i] += current[i] + logl(p);
						return;
					}
					
					else stop = logl(1-p);
				}

				if (dif > 0) {

					if (dif == 10) probability = logl(16.0 / 52.0);  // log(16/52) para {10, J, Q, K}
					else if (dif <= 11) probability = logl(4.0 / 52.0);  // log(4/52) para {1, 2, 3, 4, 5, 6, 7, 8, 9, 11}
					else probability = logl(0);  // log(0) = -infinity

					if (probability != -INFINITY) {
						
						// Verificar se o next[j] esta como default -INF -> daria problema na soma de logs
						if (next[j] == -INFINITY) next[j] = logl(1);
						
						if (stop != -INFINITY) next[j] += current[i] + probability + stop;
						
						else next[j] += current[i] + probability;
					}
				}
			}
		}
	}
}


// Devolve número aleatório de 0 a 1
double generateRandom () {
	return (double) rand() / RAND_MAX;
}


void printMatriz(long double** matriz) {
    
	for (int i = 0; i < 22; i++) {
        
		printf("matriz[%i] : ", i);

		for (int j = 0; j < 31; j++) printf("%Le ", matriz[i][j]);

		printf("\n");
	}
    
	printf("\n\n");
}
