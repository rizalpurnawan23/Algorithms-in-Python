# Bismillahirrahmanirrahim___________________________________________

#   This code contains the algorithms of common number theoretic
# expressions such as sieve of eratosthenes, fibonacci sequence
# generator, etc.

import math, statistics, time

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