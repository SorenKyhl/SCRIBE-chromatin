import matplotlib.pyplot as plt
import numpy as np
import scipy
from sklearn.metrics import log_loss


def get_entropy(s, factor, plot=False, method="relative"):
    """
    calculate entropy of a sequence after downsampling by factor
    """
    nice = []
    for e in s[::factor]:
        for _i in range(factor):
            nice.append(e)

    nice = np.array(nice)
    N = int(len(nice) / factor)

    if plot:
        plt.plot(s, "k", label="fine-grain")
        plt.plot(nice, "r", label=f"coarse-grain, N={N}")
        plt.legend()

    if method == "relative":
        return scipy.stats.entropy(abs(s), abs(nice))
    elif method == "crossentropy":
        # not correct
        return log_loss(abs(s), abs(nice))
    elif method == "KL":
        return kl_divergence(abs(s), abs(nice))


def entropy_sweep(seq):
    """
    calculate entropy of a bunch of downsampling sizes
    """
    factors = [1, 2, 4, 8, 16, 32, 64, 128, 256, 512]
    20480 / np.array(factors)

    entropies = []
    for f in factors:
        e = get_entropy(seq, f, method="relative")
        entropies.append(e)
    return np.array(entropies)


def KL(a, b):
    """
    KL divergence
    """
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)

    return np.sum(np.where(a != 0, a * np.log(a / b), 0))


# calculate the kl divergence
def kl_divergence(p, q):
    """
    another kl divergence
    """
    return np.sum([p[i] * np.log2(p[i] / q[i]) for i in range(len(p))])
