import numpy as np
#Running Statistics
class RunningStats:
    
    
    def __init__(self):
        self.n=0
        self.Mnold = self.Mnnew = self.Snold = self.Snnew = 0.0
        
    def push(self, val):
        self.n += 1
        if self.n == 1:
            self.Mnold = val
            self.Snold = 0.0
        else:
            self.Mnnew = self.Mnold + (val - self.Mnold)/((float)(self.n));
            self.Snnew = self.Snold + (val - self.Mnold)*(val-self.Mnnew);
            self.Mnold = self.Mnnew;
            self.Snold = self.Snnew;
    def mean(self):
        if self.n == 1:
            return self.Mnold
        elif self.n > 1:
            return self.Mnnew
        else:
            return 0.0
    def variance(self):
        if self.n > 1:
            vary = self.Snnew/(float(self.n)-1.0)
            return vary
        else:
            return 0.0
    def deviation(self):
        #dev = math.sqrt(self.Variance())
        dev = np.sqrt(self.variance())
        return dev
    def reset(self):
        self.n = 0