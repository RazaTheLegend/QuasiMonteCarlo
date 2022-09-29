# Quasi Monte Carlo
Quasi-Monte Carlo methods are permutations on the standard Monte Carlo method, which employs supremely uniform quasirandom numbers rather than Monte Carlo’s pseudorandom numbers. This thesis investigates the application of quasi-Monte Carlo methods on the Heston model. Our main focus in this paper is the Broadie-Kaya scheme, which our main algorithms are based on. The Monte Carlo methods provides statistical error estimates; however, this is lost in the quasi-Monte Carlo, but in return provides faster convergence than a standard Monte Carlo. A recent discovery has shown that the randomized quasi-Monte Carlo can preserve the speed of quasi-Monte Carlo, but also reintroduces the error estimates from Monte Carlo methods. For our investigation, we compare the Euler discretization with Full Truncation, the Broadie-Kaya scheme using pseudorandom sequences, then amp up the speed using quasirandom sequences and finally a randomized quasi-Monte Carlo.
