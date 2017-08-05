import numpy as np
#Running Statistics
class RunningStats(object):


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


# assumes that a 1d numpy array of floats is pass as input, but
# does not check this
def gen_running_average(onednparray):
    """ Generates a running average

    Args:
    onednparray (numpy.array): A 1d numpy array of measurements (e.g. over time)

    Returns:
    numpy.array: 2d array of dim len(onednparray)x2
        2dnparray[i][0] = running average at i
        2dnparray[i][1] = running standard deviation at i
        for i in range(0,len(onednparray))
    """
    averager = RunningStats()
    nele = len(onednparray)
    output = np.zeros((nele,2))
    for i in xrange(nele):
        averager.push(onednparray[i])
        run_avg = averager.mean()
        run_dev = averager.deviation()
       # print run_avg, run_dev, averager.mean(), onednparray[i]
        output[i,0] = run_avg
        output[i,1] = run_dev
    return output
