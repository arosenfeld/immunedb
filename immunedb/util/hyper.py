import math
import numpy as np


def choose(n, k):
    return n ** k / math.factorial(k)


def hypergeom(length, mutation, K):
    M = length
    n = K
    N = np.ceil(length * mutation)

    def pmf(k):
        if M == K or k > N:
            return 0
        return (
            choose(n, k) * choose(M - n, N - k) /
            choose(M, N)
        )
    return np.sum(
        [pmf(k) * np.power(.33, k)
            for k in range(int(np.ceil(K / 2)), K)]
    )
