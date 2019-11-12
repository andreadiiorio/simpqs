developped by Andrea Di Iorio
SIMPQS with large prime variation implementation based on Contini phD

Self Initialization Multiple Polynomial Quadratic Sieve
first precomputation will be done and will be generated POL_FAMILY_TRIES a coeffcient corresponding to same number of polynomial family
each polynomial family is composed to serveral Mongomery polynomial where the a coeff is fixed and b coeff is changed expoliting Gray Code Contini implementation
precomputation make computationally "light" to switch polynomial among polynomial family
a polynomial value identify a whole polynomial family ... if a = p1*...ps there will be 2^(s-1) polynomial families
a is generated centring a sub interval in the factor base on the ideal value of a factor as s-sqrt((2N)/M)
and taking factors in this subinterval going at right and left from this center keeping always a fixed number of factors never used before (contini explained that reduce redundant relation )
various configuration can be done in this process.
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
extended @martiani implementation for this step, also the shanks tonelly implementation used in factor base building

after that Quadratic relation x^2==y^2 mod N will be builded and factorization will be tried with GCD(X+-Y,N)

A LOT OF CONFIGURATIONs ARE POSSIBLE IN FILE CONFIGURATION.h, on top there will be SIMPQS algo config,
below concurrent configuration
worked with 60 digit number 1000000000000000000002426064040000000000001360395223980847119

PROBLEMS to solve:
-the most important problem is in the logaritm threshold  in sieving step to identify likelly to be BSmooth polynomial values
  used log(M*sqrt(N)) as written in contini phD (is pseudocode for mpqs, not mentioned for simpqs) minus MACRO tollerant error term... (actually high .. about 17 used)
    (This problem lead to useless factorization in the sieving stage slowing down a lot the program)
   -> some kind of preprocessing on top of a number to factorize can be done like:
        take a sample of polynomials inside some polynomial families, in these polynomials get logaritm cumulated value per array entry
            an average or same kind of statistic may lead to a better threshold that can be used after this preprocessing
   ... some better way is over my engineer knowledge .. Num Theory prof doesn't know ...

-a coeff generated values are distant from ideal values -> this will reduce probability of generation of BSmooth values
  on a 60 digit number a generated values
-during sieving stage polynomial make every value divisible for a so polynomial values can be divided for a achieving smaller values in whole sieving stage
     but is needed to re multiply for a (xor ing 1 to each row ) before linear algebra step -> needed for quadratic relation build

