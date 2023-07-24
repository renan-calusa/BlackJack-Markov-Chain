#include <stdio.h>
#include <stdlib.h>
#include <time.h>


int drawCard();
int checkStop(int soma);
int normal(int sum);
int playGame();

float pv;
int teta;
float pd;
float ps;


int main(int argc, char** argv) {

	srand(time(NULL));
	
	if (argc != 5) return 1;
	
	pd = strtof(argv[1], NULL);
	ps = strtof(argv[2], NULL);
	teta = atoi(argv[3]);
	pv = strtof(argv[4], NULL);
	
	int numAmostra = 30;
	
	for (int i=0; i < numAmostra; i++) printf("%i, ", playGame());
}


int drawCard() {

	// Choose a random card in the infinite deck (HIT!)
	// A: {1, 2, 3, 4, 5, 6, 7, 8, 9, 11} -> 4/52
	// B: {10} -> 16/52
	
	double randomValue = (double) rand() / RAND_MAX;
	
	double pA = 4.0/52;
	double pB = 16.0/52;

	double probabilities[11];
	for (int i=0; i < 11; i++) {
	
		if (i == 9) probabilities[i] = pB;
		else probabilities[i] = pA;
	}
	
	// Find the index of the drawn card based on probabilities
	int drawnCardIndex = 0;
	double cumulativeProbability = 0.0;
	
	while (drawnCardIndex < 11 && randomValue > cumulativeProbability + probabilities[drawnCardIndex]) {
		cumulativeProbability += probabilities[drawnCardIndex];
		drawnCardIndex++;
   	}
	
	// Return the drawn card value (Adding 1 because card values start from 1)
	return drawnCardIndex + 1;
}


int checkStop(int soma) {
	
	// Queremos parar com probabilidade p?
	
	double randomValue = (double) rand() / RAND_MAX;
	
	if (soma >= teta && randomValue <= pv) return 1;
	
	return 0;
}


int normal(int sum) {

	int valorFinal = sum;

	while (checkStop(valorFinal) == 0 && valorFinal <= 21) valorFinal += drawCard();

	return valorFinal;
}


int playGame() {
	
	int primeira = drawCard();
	int segunda = drawCard();
	
	int somaInicial = primeira + segunda;
	
	// verificar Double-Down ou Split
	
	if (somaInicial == 9 || somaInicial == 10 || somaInicial == 11) return somaInicial + drawCard();
	
	else if (primeira == segunda) {
	
		int jogoA = primeira + drawCard();
		int jogoB = segunda + drawCard();
		
		int finalA = normal(jogoA);
		int finalB = normal(jogoB);
		
		return finalA <= finalB ? finalA : finalB;
	}
	
	else {
	
		return normal(somaInicial);
	}
}
