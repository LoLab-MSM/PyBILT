import numpy as np



class Gaussian(object):
    def __init__(self, mean,std):
        stdinv = 1.0/std
        normalc = stdinv*(1.0/np.sqrt(np.pi))
        self.sigma = std
        self.mean = mean
        self.normconst = normalc
        return


    def eval(self,x_in):
        stdinv = 1.0/self.sigma
        stdinvsq = stdinv**2
        normalc = self.normconst
        expon = -(x_in - self.mean)**2 * (0.5*stdinvsq)
        y = normalc * np.exp(expon)
        return y

    def reset_mean(self,new_mean):
        self.mean = new_mean
        return


class GaussianRange(object):
    def __init__(self,in_range,mean,std,npoints=200):
        x_p = np.linspace(in_range[0],in_range[1],npoints,endpoint=True)
        y_p = np.zeros(npoints)
        yc = 0
        stdinv = 1.0/std
        stdinvsq = stdinv**2
        normalc = stdinv*(1.0/np.sqrt(np.pi))
        for x in x_p:
            expon = -(x - mean)**2 * (0.5*stdinvsq)
            y = normalc * np.exp(expon)
            y_p[yc]=y
            yc+=1
        self.x = x_p
        self.y = y_p
        self.sigma = std
        self.mean = mean
        self.normconst = normalc
        self.upper = in_range[1]
        self.lower = in_range[0]
        self.dx = x_p[1]-x_p[0]
        self.npoints = npoints
        return

    def get_values(self):
        return (self.x,self.y)

    def eval(self,x_in):
        stdinv = 1.0/self.sigma
        stdinvsq = stdinv**2
        normalc = self.normconst
        expon = -(x_in - self.mean)**2 * (0.5*stdinvsq)
        y = normalc * np.exp(expon)
        return y

    def integrate_range(self, lower, upper):
        if upper>self.upper:
            upper=self.upper
        if lower<self.lower:
            lower = self.lower

        i_l = int(np.floor((lower-self.lower)/self.dx))
        i_u = int(np.floor((upper-self.lower)/self.dx))
        #print "i_l ",i_l," i_u ",i_u
        total = 0.0
        for i in xrange(i_l,i_u):
            total+= self.y[i]*self.dx
        return total

    def sum_range(self, lower, upper):
        if upper>self.upper:
            upper=self.upper
        if lower<self.lower:
            lower = self.lower

        i_l = int(np.floor((lower-self.lower)/self.dx))
        i_u = int(np.floor((upper-self.lower)/self.dx))
        total = 0.0
        for i in xrange(i_l,i_u):
            total+= self.y[i]
        return total

    def normalize(self):
        total = 0.0
        for i in xrange(0,self.npoints):
            total+=self.y[i]*self.dx
        for i in xrange(0,self.npoints):
            self.y[i]/=total
        return

    def reset_mean(self,new_mean):
        self.mean = new_mean
        return
