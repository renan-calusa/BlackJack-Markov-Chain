#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

long double** init();
void set_default(long double** matriz);
void estimator(int* amostra);
long double likelihood(int* v, float p, int teta, long double** matriz);
long double probability_function (int v, float p, int teta, long double** matriz);
void transition(long double** matriz, int state, float p, int teta);
long double logaddexp(long double x, long double y);
void printMatriz(long double** matriz);
double generateRandom();


// Kill switch para transição
int stop;
int sample_size;


int main (int argc, char** argv) {

    float p; // entre 0 e 1
    int teta; // entre 1 e 21
    int* amostra;
    
    sample_size = atoi(argv[1]);
    
    if (argc == sample_size + 4) {
            
            // Caso que se quer calcular o likelihood de um agente descrito
            
            amostra = (int*) malloc(sizeof(int)*sample_size);
            for (int i=0; i < sample_size; i++) amostra[i] = atoi(argv[i+2]);  // entre 0 e 30
            
            p = strtof(argv[argc-2], NULL);
            teta = atoi(argv[argc-1]);

            if (p < 0 || p > 1 || teta < 1 || teta > 21) {
                printf("[-] Valores passados fora do intervalo esperado!");
                return 0;
            }
            
            long double** matriz = init();
    
            srand(time(NULL));
    
            printf("likelihood(amostra; %f, %i) = %Le\n", p, teta, likelihood(amostra, p, teta, matriz));
    }

    else if (argc == sample_size + 1) {
    
    	// Caso que se quer estimar 'teta' e 'p' de uma amostra
    	
    	amostra = (int*) malloc(sizeof(int)*sample_size);
        for (int i=0; i < sample_size; i++) amostra[i] = atoi(argv[i+2]);  // entre 0 e 30
        
        estimator(amostra);
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
	
	for (int i = 0; i < 22; i++) for (int j = 0; j < 31; j++) matriz[i][j] = -INFINITY;   // equivalente à Pr(0)

	matriz[0][0] = 0; // (log(1) = 0) 100% de estar com zero no estado zero

	return matriz;
}


void set_default(long double** matriz) {

	for (int i = 0; i < 22; i++) for (int j = 0; j < 31; j++) matriz[i][j] = -INFINITY;   // equivalente à Pr(0)

	matriz[0][0] = 0; // (log(1) = 0) 100% de estar com zero no estado zero
}


// Estimar 'teta' e 'p' por BruteForce
void estimator(int* amostra) {

	int teta;
	float p;
	
	int max_teta = 0;
	float max_p = 0;
	long double max_likelihood = -INFINITY;
	
	long double** matriz = init();
	
	for (p = 0; p <= 1; p += 0.05) {
	
		for (teta = 1; teta <= 21; teta++) {
		
			set_default(matriz);
		
			long double curr = likelihood(amostra, p, teta, matriz);
			
			if (curr > max_likelihood) {
			
				max_p = p;
				max_teta = teta;
			}
		}
	}
	
	printf("[^]teta estimator: %i\t[^]p estimator: %.3f\n", max_teta, max_p);
}


// Somatorio do log das probabilidades -> Sum i of (log Pr(vi|teta, p)
long double likelihood (int* v, float p, int teta, long double** matriz) {

    long double likelihood = -INFINITY;

    for (int i=0; i < sample_size; i++) likelihood = logaddexp(logl(probability_function(v[i], p, teta, matriz)), likelihood);
    
    return likelihood;
}


long double probability_function (int v, float p, int teta, long double** matriz) {

    stop = 0;

    // Realizar transicoes de markov até o jogador quiser parar no valor "teta" com probabilidade "p"
    for (int i=0; i < 21 && stop != 1; i++) transition(matriz, i, p, teta);

    // Probabilidade acumulada de cada estado para certo valor "v":
    long double res = 0.0;

    for (int i=1; i < 22; i++) res += expl(matriz[i][v]);
    
    printf("Pr(%i, %.3f, %i) = %Lf\n", v, p, teta, res);

    return res;
}


void transition(long double** matriz, int state, float p, int teta) {
    
    // Não existe um estado maior que 21
    if (state >= 21) return;

    long double* current = matriz[state];
    long double* next = matriz[state + 1];

    // Se pode-se alcançar teta, verifica se queremos parar
    if (current[teta] != -INFINITY && generateRandom() < p) {
        stop = 1;
        return;
    }

    for (int i = 0; i < 22; i++) {
        for (int j = 0; j < 31; j++) {
            
            if (current[i] != -INFINITY) {

                long double probability;
                int dif = j - i;

                if (dif > 0) {
                    
                    if (dif < 12 && dif != 10) probability = logl(4.0 / 52.0);  // log(4/52) para {1, 2, 3, 4, 5, 6, 7, 8, 9, 11}
                    else if (dif == 10) probability = logl(16.0 / 52.0);  // log(16/52) para {10, J, Q, K}
                    else probability = -INFINITY;  // log(0) = -infinity

                    next[j] = logaddexp(next[j], current[i] + probability);  // log(current[i]*probability) = log(current[i]) + log(probability)
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

        for (int j = 0; j < 31; j++) printf("%Le ", expl(matriz[i][j]));

        printf("\n");
    }
    
    printf("\n\n");
}


// Evitar o underflow
long double logaddexp(long double x, long double y) {
    
    if (x == y) return x + logl(2.0);
    
    else {
        long double max_val = fmaxl(x, y);
        long double min_val = fminl(x, y);
        return max_val + logl(1.0 + expl(min_val - max_val));
    }
}
