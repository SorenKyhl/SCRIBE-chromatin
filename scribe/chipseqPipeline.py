import numpy as np
import scipy.ndimage

from scribe import epilib


class ChipseqPipeline:
    """data processing pipeline for chipseq data"""

    def __init__(self, operations):
        self.operations = operations

    def fit(self, x):
        for operation in self.operations:
            x = operation.operate(x)
        return x


class Smooth:
    """gaussian smoothing"""

    def __init__(self, size=2):
        self.size = size

    def operate(self, x):
        return scipy.ndimage.gaussian_filter(x, self.size)


class Normalize:
    """map sequence to the range [0,1]"""

    def __init__(self):
        pass

    def operate(self, x):
        # parameter: scaling not exposed
        normalized = epilib.new_map_0_1_chip(x)
        centered = 2 * normalized - 1
        return centered


class Sigmoid:
    """sigmoid transformation"""

    def __init__(self, w=20, b=10):
        self.w = w
        self.b = b

    def operate(self, x):
        # perhaps this should just be a tanh?
        def sigmoid_fn(x):
            return 1 / (1 + np.exp(-x))

        x = sigmoid_fn(self.w * x + self.b)
        x = 2 * x - 1
        return x
