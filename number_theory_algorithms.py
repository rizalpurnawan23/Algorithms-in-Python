# Bismillahirrahmanirrahim___________________________________________

#   This code contains the algorithms of common number theoretic
# expressions such as sieve of eratosthenes, fibonacci sequence
# generator, etc.

import math, statistics

# THE MAIN ALGORITHM_________________________________________________
class NumThAlgo:
    def __init__(self):
        pass

    def readMe(self):
        return ("""
        Class Title     : NumThAlgo
        First Creation  : 27/06/2022
        Author          : Rizal Purnawan

        ___Trivia

            This class contains the algorithms of common number
        theoretic expressions.
        """)

    # Prime Numbers Algorithms_______________________________________
    def soE(self, n):
        """
        Method's Title      : soE
        Formal Title        : Sieve of Eratosthenes
        Author              : Rizal Purnawan

        ___Trivia

            This method helps us generate a list consisting of prime
        numbers strictly less than the given positive integer "n".
        The generation is done by making use of the famous "Sieve of
        Eratosthenes".
        """
        if isinstance(n, int) and n > 1:
            siever = [True] * n
            siever[0], siever[1] = False, False
            k = 2
            while k < n:
                if siever[k]:
                    for i in range(2*k, n, k):
                        siever[i] = False
                k = k + 1
            return [i for i in range(n) if siever[i]]
        else:
            print("Invalid argument!")
            print("Give only integers greater than 1!")
            raise ValueError

    def factorize(self, n):
        """
        Method's Title      : factorize
        Formal Title        : Prime Factorization
        Author              : Rizal Purnawan

        ___Trivia

            This method provide us the prime factorization of a given
        integer n. It returns a list. The first entry in the list
        will be either 1 or -1 which indicates the sign of n. And
        the rest of the entries are the prime factors of n. For
        instance, if the prime factors of an integer a > 0 are
        p1, p2, ..., pk, then the method will return a list of the
        form

                [1, p1, p2, ..., pk] .
        """
        if isinstance(n, int):
            fact_list = []
            if n >= 0:
                fact_list.append(1)
            else:
                fact_list.append(-1)
            m = abs(n)
            # Now we need to prime factorize m. To do so, we need to
            # list the primes that are possibly the factors of m,
            # and then we check whether the primes in the list are
            # the factors of m or not.
            # Note that we do not need to list all the primes less
            # than m for the candidate primes. If m is divisible by
            # some positive integer k, then we only need to list the
            # primes less than m/k for the candidate primes. If all
            # these primes do not divide m, then m must be a prime.
            # For the algorithm, we may try k for 2, 3, 5, or 7 for
            # efficiency.
            cand_primes = None
            for k in [7, 5, 3]:
                if m > k and m % k == 0:
                    cand_primes = self.soE( math.ceil( m/k ) )
                    break
            if cand_primes is None:
                cand_primes = self.soE( math.ceil( m/2 ) )
            pfact_list = [p for p in cand_primes if m % p == 0]
            if len(pfact_list) > 0:
                return fact_list + pfact_list
            else:
                fact_list.append(m)
                return fact_list
        else:
            print("Give only integer input!")
            raise ValueError

    def is_prime(self, n):
        """
        Method's Title      : is_prime
        Formal Title        : Prime Identifier
        Author              : Rizal Purnawan

        ___Trivia

            This method tells us whether a given input n is a prime
        or not. Hence the output of this method is a boolean value.
        """
        siever = self.factorize(n)
        if n > 1 and n in siever:
            return True
        else:
            return False

    def __pcounter(self, n):
        if isinstance(n, int) and n > 1:
            return len(self.soE(n + 1))

    def prime_counter(self, n, m= None):
        """
        Method's Title      : prime_counter
        Formal Title        : Prime Counting Function
        Author              : Rizal Purnawan

        ___Trivia

            This method help us count the number of primes number in
        between m and n inclusive. If m is 'None', then the method
        counts the number of primes less than or equal to n.
        """
        if isinstance(n, int) and n > 1:
            if m is None or (isinstance(m, int) and m <= 1):
                return self.__pcounter(n)
            elif isinstance(m, int) and 1 < m < n:
                return self.__pcounter(n) - self.__pcounter(m - 1)
            else:
                raise ValueError
        else:
            raise ValueError

    # Fibonacci Sequence Generator___________________________________
    def fibonacci(self, n):
        """
        Method's Title      : fibonacci
        Formal Title        : Fibonacci Sequence Generator
        Author              : Rizal Purnawan

        ___Trivia

            This method generates a list consisting of the first
        n terms of fibonacci sequence, given a nonnegatie integer n
        as an input. Note that we use the convention that the
        the order of the fibonacci sequence started from 0.
        """
        if isinstance(n, int) and n >= 0:
            if n == 0:
                return [0]
            elif n == 1:
                return [0, 1]
            else:
                fib_list = [0, 1]
                k = 2
                while k <= n:
                    fib_list.append(
                        fib_list[k - 2] + fib_list[k - 1]
                    )
                    k = k + 1
                return fib_list
        else:
            print("Invalid input!")
            print("Give only a nonnegative integer!")
            raise ValueError

    def __fibo_f(self, n):
        if isinstance(n, int):
            if n >= 2:
                fibo_li = [0, 1]
                k = 2
                while True:
                    f = fibo_li[k - 2] + fibo_li[k - 1]
                    if f >= n:
                        break
                    fibo_li.append(f)
                    k = k + 1
                return fibo_li
            elif n == 1:
                return [0]
            else:
                return []
        else:
            print("Invalid input!")
            print("Give only an integer!")
            raise ValueError


    def fibo_filter(self, upp, low= 0):
        """
        Method's Title      : fibo_filter
        Formal Title        : Fibonacci Sequence Filter
        Author              : Rizal Purnawan

        ___Trivia

            This method filters the fibonacci sequence in betweeen
        two positive integers m and n greater than 1. The method
        returns a list consisting of the fibonacci sequence greater
        than or equal to m and less than n.
        """
        fibo_m = self.__fibo_f(low)
        fibo_n = self.__fibo_f(upp)
        return [f for f in fibo_n if not (f in fibo_m)]
        
    # Arithmetic Sequence and Series_________________________________
    def arithmetic_prog(self, n, a= 1, d= 1, series= False):
        """
        Method's Title      : arithmetic_prog
        Formal Title        : Arithmetic Progression
        Author              : Rizal Purnawan

        ___Trivia

            This method provides the generation of arithmetic
        sequence and its sum. If one set 'True' for argument
        'series', then the method returns the sum of the arithmetic
        sequence. Otherwise the method returns a list consisting
        the arithmetic sequence.
        """
        if isinstance(n, int) and n > 0 and isinstance(d, int):
            if series == False:
                return [a + (k - 1) * d for k in range(1, n + 1)]
            else:
                a_n = a + (n - 1) * d
                return (n / 2) * (a_n + a)
        else:
            print("Invalid input!")
            raise ValueError

    # Geometric Progression__________________________________________
    def geometric_prog(self, n, a, r, series= False, infinite= False):
        """
        Method's Title      : geometric_prog
        Formal Title        : Geometric Progression
        Author              : Rizal Purnawan

        ___Author

            This method proviedes the generation of geometric
        sequence, its partial sum, and its infinite series if the
        sequence converges.
        
        Parameters:
        n       : a positive integer indicating the n-term of the
                  progression.
        a       : a real number (floating point) --in general--
                  representing the first term of the sequence.
        r       : a real number (floating point) --in general--
                  representing the common ratio of the sequence.
        series  : if True, the method returns the sum of the
                  sequence, otherwise the method returns a list
                  of numbers witin the sequence.
        infinite: Can be True if series is True, and in that case,
                  the method will return the sum of the infinite
                  series if the sequence converges. Otherwise the
                  method returns an error message.
        """
        if isinstance(n, int) and n > 0 and all(isinstance(x, int)
                or isinstance(x, float) for x in [a, r]) \
                and r != 1:
            if series == False and infinite == False:
                return [a *r**(k - 1) for k in range(1, n + 1)]
            elif series == True and infinite == False:
                return a* ( (1 - r**n)/(1 - r) )
            elif series == True and infinite == True:
                if 0 < abs(r) < 1:
                    return a/(1 - r)
            else:
                raise ValueError
        else:
            print("Invalid input!")
            raise ValueError

    # Coprimitive____________________________________________________
    def __gcd(self, a, b):
        # Alternative gcd algorithm
        if all(isinstance(n, int) and n >= 0 for n in [a, b]):
            a_pf = self.factorize(a)[1:]
            b_pf = self.factorize(b)[1:]
            pf = [n for n in a_pf if n in b_pf]
            num_pf = len(pf)
            k = 0
            div_ = 1
            while k < num_pf:
                div_ = div_ * pf[k]
                k = k + 1
            # div_ = math.prod(pf)
            return div_


    def is_coprime(self, a, b):
        """
        Method's Title      : is_coprime
        Formal Titile       : Comprimitive Identifier
        Author              : Rizal Purnawan

        ___Trivia

            This method helps us indentify whether a given pair of
        integers is coprime. The arguments are 'a' and 'b'. If 'a'
        and 'b' are pariwise coprime, the method will return
        a boolean value 'True', otherwise it will return 'False'.
        """
        if all(isinstance(n, int) and n >= 0 for n in [a, b]):
            return True if math.gcd(a, b) == 1 else False
        else:
            print("Invalid input!")
            print("'a' and 'b' have to be nonnegative integers!")
            raise ValueError

    def coprimes(self, b, a= 1, existence= False):
        """
        Method's Title      : coprimes
        Formal Title        : Coprime Numbers Generator
        Author              : Rizal Purnawan

        ___Trivia

            This method has several functions related to coprimitives
        of a pair of integers. Arguments 'a' and 'b' must be
        integers >= 1, 'a' < 'b' is a must, and by default 'a'
        is set to be 1. If argument 'existence' is 'False' --which is
        the default setting-- then the method will return a list of
        integers in between 'a' and 'b' that are coprime with both
        'a' and 'b'. If such integers do not exists, then it returns
        an empty list. On the other hand, if 'existence' is 'True',
        then the method will return a boolean value. It returns 'True'
        if there is at least an integer in between 'a' and 'b' that
        are coprime to both 'a' and 'b'. Then it returns 'False'
        otherwise.
        """
        if all(isinstance(n, int) and n >= 0 for n in [a, b]):
            gcd = math.gcd
            if existence == False:
                return [n for n in range(a, b)
                            if gcd(n, a) == 1 and gcd(n, b) == 1]
            elif existence == True:
                return (
                    True if any(
                        gcd(n, a) == 1 and gcd(n, b)
                        for n in range(a, b)
                        )
                    else False
                )
        else:
            print("Invalid input!")
            raise ValueError
