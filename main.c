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
long double logaddexpl(long double x, long double y);
int search(int* array, int livre, int value);

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

		likelihood(amostra, p, teta, matriz);
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

	for (int i = 0; i < 22; i++) for (int j = 0; j < 31; j++) matriz[i][j] = 0;

	matriz[0][0] = 1;
}


// Estimar 'teta' e 'p' por BruteForce
void estimator(int* amostra, long double** matriz) {

	int teta;
	float p;
	
	int max_teta = 0;
	float max_p = 0;
	long double max_likelihood = -INFINITY;
	
	for (teta = 1; teta <= 21; teta++) {
		for (p = 0.05; p < 1.05; p += 0.05) {
		
			long double curr = likelihood(amostra, p, teta, matriz);
			
			if (curr > max_likelihood) {
			
				max_likelihood = curr;
				max_p = p;
				max_teta = teta;
			}
		}
	}
	
	printf(">>teta estimator: %i\t>>p estimator: %.2f\n", max_teta, max_p);
}


// Somatorio do log das probabilidades -> Sum i of (log Pr(vi|teta, p)
long double likelihood (int* v, float p, int teta, long double** matriz) {

	long double likelihood = logl(1);

	for (int i=0; i < sample_size; i++) {
	
		long double prob = probability_function(v[i], p, teta, matriz);
		
		if (prob == 0) {
		
			likelihood = -INFINITY;
			break;
		}
		
		likelihood = logl(expl(likelihood) + prob);
	}
	
	// Print info
	printf("likelihood(%.2f, %i; {%i", p, teta, v[0]);
	for (int i=1; i < sample_size; i++) printf(", %i", v[i]);
	if (likelihood == -INFINITY) printf("}) = IMPOSSIBLE\n");
	else printf("}) = %Le\n", likelihood);
    
	return likelihood;
}


long double probability_function (int v, float p, int teta, long double** matriz) {
	
	set_default(matriz);
	
	// Realizar transicoes de markov até o jogador quiser parar no valor "teta" com probabilidade "p"
	for (int i=0; i < 21; i++) transition(matriz, i, p, teta);
	
	//printMatriz(matriz);
	
	// Pega a probabilidade de terminar com um valor v num jogo de BlackJack - em log()
	long double res = matriz[21][v];
	
	printf("Pr(%i; %.2f, %i) = %Lf\n", v, p, teta, res);

	return res;
}


void transition(long double** matriz, int state, float p, int teta) {
    
	// Não existe um estado maior que 21
	if (state >= 21) return;
		
	long double* current = matriz[state];
	long double* next = matriz[state+1];
	
	int* usedValues = (int*) malloc(sizeof(int)*22);
	for (int i=0; i < 22; i++) usedValues[i] = 0;
	int livre = 0;
	
	// Percorrendo todos valores no atual (i) e no proximo (j)
	for (int i = 0; i < 31; i++) {
		for (int j = 0; j < 31; j++) {
            
			if (current[i] != 0) {
				
				double probability;
				int dif = j - i;
				
				if (dif == 10) probability = 16.0/52;
				else if (dif <= 11) probability = 4.0/52;
				else probability = 0;	
				
				if (dif > 0 && probability != 0) {
					
					if (i < teta) next[j] += current[i] * probability;
					
					else {
					
						if (usedValues[livre] == i) next[j] += current[i] * probability * (1 - p);
					}
				}
				
				else if (i == j && i >= teta) {
				
					if (search(usedValues, livre, i) == -1) next[i] += current[i] * p;
					
					else {
						
						usedValues[livre] = i;
						livre++;
						
						next[i] += current[i];
					}
				}
			}
		}
	}
	
	free(usedValues);
}


int search(int* array, int livre, int value) {

	for (int i=0; i < livre; i++) if (array[i] == value) return i;
	
	return -1;
}


void printMatriz(long double** matriz) {
    
	for (int i = 0; i < 22; i++) {
        
		printf("matriz[%i] : ", i);

		for (int j = 0; j < 31; j++) printf("%Lf ", matriz[i][j]);

		printf("\n");
	}
    
	printf("\n\n");
}


long double logaddexpl(long double x, long double y) {
    // Calculate the maximum of x and y
    long double max_val = (x > y) ? x : y;

    // Calculate the minimum of x and y
    long double min_val = (x < y) ? x : y;

    // Calculate the result using the log-sum-exp formula
    return max_val + logl(1.0 + expl(min_val - max_val));
}
