import sys

# For all classes in this file, a specifies the desired inner
# boundary, b specifies the desired outer boundary.


class linear:
    def __init__(self, a, b):
        if b <= a:
            sys.exit("Outer boundary b must be larger than inner boundary a")
        self.A = 0.5 * (b - a)
        self.B = 0.5 * (b + a)
        self.a = a
        self.b = b

    def rescale_grid(self, x):
        return self.A * x + self.B

    def jacobian(self, x):
        return 1 / self.A


# C is an extra parameter for moving the density distribution of points.
class harald_scaling:
    def __init__(self, a, b, C):
        if b <= a:
            sys.exit("Outer boundary b must be larger than inner boundary a")
        self.A = (a + b + 2.0 * C) / (b - a)
        self.B = b - b * self.A + C - C * self.A
        self.C = C
        self.a = a
        self.b = b

    def rescale_grid(self, x):
        return self.B / (x - self.A) - self.C

    def jacobian(self, x):
        return -((x - self.A) ** 2.0) / self.B
