# Bismillahirrahmanirrahim-------------------------------------------

# PREFACE____________________________________________________________
#   This module contains the algorithms for dealing with vector-like
# objects. It is too naive to say that it is a (real) vector space
# algorithm, since we can never really construct a (real) vector
# space in computers, since algebraic closure will never be attained
# with finite set --floating point numbers are finite.

#   In this algorithm, we will define a vector as a list of numbers.
# In general the numbers are floating-point numbers. However,
# integers are also compatible.

class VectorSpace:
    def __init__(self, *args):
        self.args = args

    def read_me(self):
        """
            This class provides the machinary to deal with vectors
        in Python, which is in fact lists of numbers. The methods
        in this class include the basic tools for manipulating
        vectors such as addition, scalar multiplication, inner
        product, metrics, norms and normalisation.

        Author,
        Rizal Purnawan.
        """
    
    # Auxiliary Methods______________________________________________
    def __is_num(self, x):
        """
        Description:
            This method returns 'True' if the parameter 'x' is either
        an integer or a floating point number.
        """
        if isinstance(x, int) or isinstance(x, float):
            return True
        else:
            return False
    
    def __is_vector(self, v):
        """
        Description:
            This method helps us identify whether a given argument
        is considered a vector. We define a vector in Python as a
        list of numbers (integers or floating point numbers).
        """
        if isinstance(v, list) and all(isinstance(x, int)
                or isinstance(x, float) for x in v):
            return True
        else:
            return False

    def __all_vectors(self):
        """
        Description:
            This method identifies whether all class arguments are
        vectors or not. It returns 'True' if they are, and returns
        'False' otherwise.
        """
        if all(self.__is_vector(v) for v in self.args):
            return True
        else:
            return False

    def __all_vectors_same_dim(self):
        """
        Description:
            This method identifies whether all class arguments are
        vectors and are of same dimension or not. It returns 'True'
        if they are, and returns 'False' otherwise.
        """
        if self.__all_vectors() and all(len(self.args[0]) == len(v)
                for v in self.args):
            return True
        else:
            return False
    
    def __bin_which(self, dim, which):
        if isinstance(which, list) and len(which) == 2 \
                and all(isinstance(w, int)
                and (0 <= w < dim or -dim <= w < 0)
                for w in which):
            return True
        else:
            return False

    def __error_notif(self):
        print("Invalid arguments!")
        raise ValueError

    # Generating Vectors_____________________________________________
    def create_v(self):
        """
        Description:
            This method can only be applied if the class arguments
        are numbers. In this case, this method will return a list
        of the numbers within the class arguments, and we consider it
        as a vector.
        """
        if all(self.__is_num(v) for v in self.args):
            return [v for v in self.args]
        else:
            print("This method is inapplicable.")
            raise ValueError

    # Vector Summation_______________________________________________
    def vsum(self, which= "all"):
        """
        Description:
            If the class arguments are vector, this method helps us
        sum the selected vector in the arguments. The parameter
        'which' is used to select which vectors to be summed.
        The parameter must be one of the following values:
        'all'           : sum all vectors if all the arguments are
                          vectors.
        'first_bin'     : sum the first two vectors in the arguments.
        A list          : sum the vectors whose orders are in the
                          list.
        'm to n'        : sum vectors whose orders are less than or
                          equal to m and less than n, for any
                          integers m and n with m < n.
        """
        if self.__all_vectors_same_dim():
            dim = len(self.args[0])
            if which == "all":
                summand = list(self.args)
            elif which == "first_bin":
                summand = [self.args[0], self.args[1]]
            elif isinstance(which, list) and all(isinstance(w, int)
                    and (0 <= w < dim or -dim <= w < 0)
                    for w in which):
                summand = [self.args[w] for w in which]
            elif isinstance(which, str) and " to " in which:
                w_list = which.split(" to ")
                try:
                    a, b = w_list
                    a, b = int(a), int(b)
                    if 0 <= a < b <= dim or -dim <= a < b < 0 \
                            or (0 <= a < dim
                                    and -dim <= b < 0
                                    and a < dim - b):
                        summand = list(self.args)[a: b]
                except:
                    try:
                        a = int(w_list[0])
                        if isinstance(a, int) and 0 <= a < dim \
                                and b == "end":
                            summand = list(self.args)[a:]
                        else:
                            self.__error_notif()
                    except:
                        self.__error_notif()
            else:
                self.__error_notif()
            return [sum([v[k] for v in summand]) for k in range(dim)]
        else:
            self.__error_notif()

    # Subtraction____________________________________________________
    def vsubtract(self, which= "first_bin"):
        """
        Description:
            This method helps us perform subtraction on a pair of
        vectors. The parameter which refers to which vectors in the
        class arguments for subtraction. The parameter must be one
        of the following values:
        'first_bin'         : subtract the first two vectors in the
                              arguments.
        A list              : this list must contain only two
                              integers whose values are the orders
                              of the vectors. Then the method
                              subtract the vectors with these orders.
        """
        if self.__all_vectors_same_dim() and len(self.args) >= 2:
            dim = len(self.args[0])
            if which == "first_bin":
                for_subt = [self.args[0], self.args[1]]
            elif self.__bin_which(dim, which):
                for_subt = [self.args[w] for w in which]
            else:
                return self.__error_notif()
            return [for_subt[0][k] - for_subt[1][k]
                        for k in range(dim)]
        else:
            self.__error_notif()

    # Scalar Multiplication__________________________________________
    def scal_x(self, c, which= 0):
        """
        Description:
            This method provides the scalar multiplication of a given
        vector in the class arguments. The parameter 'c' must be
        a number as the scalar. The parameter 'which' refers to
        the order of the vector in the arguments to be scalar
        multiplied.
        """
        if self.__all_vectors() and self.__is_num(c):
            try:
                vec = list(self.args)[which]
                return [c * x for x in vec]
            except:
                self.__error_notif()
        else:
            self.__error_notif()

    # Metric of A Pair of Vectors____________________________________
    def metric(self, which= 'first_bin', p= 2):
        """
        Description:
            This method provides the distance between a pair of
        vectors in the arguments. The parameter 'which' refers to
        the order of the selected vectors in the arguments. While
        the parameter 'p' is a number greater than or equal to 1. By
        default we set 'p = 2' which means that the applied distance
        is the euclidean metric. If 'p = 1', then it is a Manhattan
        taxicab metric.
            The parameter 'which' must be one of the following:
        
        'first_bin'     : the first two vectors in the class
                          arguments.
        A list          : the list must be of two integers which
                          are also the orders of some two vectors
                          in the class arguments.
        """
        if self.__all_vectors_same_dim() and self.__is_num(p) \
                and len(self.args) >= 2 and p >= 1:
            dim = len(self.args[0])
            if which == 'first_bin':
                v_pair = [self.args[0], self.args[1]]
            elif self.__bin_which(dim, which):
                v_pair = [self.args[w] for w in which]
            else:
                self.__error_notif()
            a, b = v_pair
            return sum([abs(x - y)**p for x, y in zip(a, b)]) **(1/p)
        else:
            self.__error_notif()

    # Inner Product__________________________________________________
    def inner_prod(self, which= 0):
        """
        Description:
            This method computes the inner product of a given pair of
        vectors from the class arguments. The parameter 'which' must
        be a list of the orders of the selected vectors in the class
        arguments. However, we set a default for the 'which'
        parameter to be 'first_bin' that will give us the first two
        vectors in the class arguments.
        """
        if self.__all_vectors_same_dim() and len(self.args) >= 2:
            dim = len(self.args[0])
            if which == 'first_bin':
                v_pair = [self.args[0], self.args[1]]
            elif self.__bin_which(dim, which):
                v_pair = [self.args[w] for w in which]
            else:
                self.__error_notif()
            a, b = v_pair
            return sum([a[k] * b[k] for k in range(dim)])
        else:
            self.__error_notif()

    # Norm___________________________________________________________
    def norm(self, which= 0, p= 2):
        """
        Description:
            This method computes the norm of a given vector in the
        class arguments. The parameter 'which' must be the order of
        the selected vector. And the parameter 'p' must be a number
        greater than or equal to 1. By default, we set 'p = 2', which
        means that the norm will be the euclidean norm. If 'p = 1'
        the norm is known as the Manhattan taxicab norm.
        """
        if self.__all_vectors() and self.__is_num(p) and p >= 1:
            try:
                vec = self.args[which]
                return sum([abs(x)**p for x in vec]) **(1/p)
            except:
                self.__error_notif()
        else:
            self.__error_notif()

    # Normalisation__________________________________________________
    def normalise(self, which= 0):
        """
        Description:
            This method normalises a given vector in the class
        argument. The parameter 'which' must be the order of the
        selected vector.
        """
        if self.__all_vectors():
            try:
                vec = self.args[which]
                divisor = self.norm(which= which)
                return [x / divisor for x in vec]
            except:
                self.__error_notif()
        else:
            self.__error_notif()
