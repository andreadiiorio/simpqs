developped by Andrea Di Iorio
SIMPQS with large prime variation implementation based on Contini phD
Self Initialization Multiple Polynomial Quadratic Sieve
first precomputation will be done and will be generated POL_FAMILY_TRIES a coeffcient corresponding to same number of polynomial family
each polynomial family is composed to serveral Mongomery polynomial where the a coeff is fixed and b coeff is changed expoliting Gray Code Contini implementation
precomputation make computationally "light" to switch polynomial among polynomial family
In This version Sieving phase is done concurrently on different process that will Sieve for different Polynomial Families
         ( POLYNOMIAL_FAMILIES_CONCURRENT_SIEVERS process will concurrently sieve for POLYNOMIAL_FAMILIES_PER_SIEVER polynomial families)
During 1 polynomial Sieve inside 1 process different thread will sieve the Array concurrently using log sieving as suggested by contini
array entry that are likelly to be BSmooth will be Factorized and only entry that decompose in Factorbase Primes up to B SmoothnessBound will be used
    (or entry that are like the previusly but have 1 large prime like B* LARGE_PRIME_THREASHOLD_FACTOR)
Factorization is done in another set of thread inside the same process with an approch that I've invented callable ThreadGroupPool
    each likellyToBeBsmooth entry will be Enqueued in a SYNC queue and a ThreadGroupManager will dequeue it
    expoliting 2 barrier SYNC with a fixed number of iteration factorization of the entry will be searched concurrently in the Thread ThreadGroupManager
    see factorizerQuick.c for implementation.
After each sieving process (partial) relation will be collected and aggregated by the main process
After sieving of each polynomial inside each generated polynomial family large primes will coupled generating new relations (like suggested in Crandal & Pomerance book)
when at least cardinality of FB relation will be collected plus an extra Linear algebra step will be performed, solving a liner system of exponenet vectors of founded relations
exponent vectors are exponents % 2 of each primes in factorization of founded relation
Linear system solved with Gauss elimination algorithm with the optimization of using XOR at each row sum (in GF2 so sum is logically equivalent to XOR)
after that Quadratic relation x^2==y^2 mod N will be builded and factorization will be tried with GCD(X+-Y,N)

