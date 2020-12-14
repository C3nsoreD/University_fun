import numpy as np


## Constants 

## Object

class Element():
    def __init__(self, angle):
        self.c = self._c(angle)
        self.s = self._s(angle)
        self.E = None
        self.A = None
        self.length = None
    
    def _c(self, angle):
        return np.cos(angle)

    def _s(self, angle):
        return np.sin(angle)

    def _transformK(self):
        value = np.array([
            [self.c, -self.c], 
            [self.s, -self.s],
            [-self.c, self.c],
            [-self.s, self.s]
        ])
        return value

    def _transform(self):
        value = np.array([
            [self.c, self.s, 0, 0],
            [0, 0, self.c, self.s],
        ])
        return value

    
e = Element(0)
print(e.c, e.s)
