# BlackJack-Markov-Chain

Analysis of the final hand of a player in matches of BlackJack with card reposition

Each player has 2 attributes:
1. teta: The sum that he stops asking for more cards
2. p: The probability that the player will follow the 'teta' rule

This code can calculate the likelihood for an sample giving a 'teta' and a 'p' or Estimate those parameters based of the sample.

# Compiling the code
gcc main.c -o main -lm -Wall

# Running the code
1. Likelihood:
./main <sample_size> <sample separate by 'space bar'> <p> <teta>

ex: ./main 4 19 20 25 9 0.7 10
* This will print the log likelihood of  sample = {19, 20, 25, 9}; teta = 10; p = 0.7

2. Estimating:
./main <sample_size> <sample separate by 'space bar'>

ex: ./main 4 19 20 25 9
* This will print the best 'teta' and 'p' for the sample = {19, 20, 25, 9}
