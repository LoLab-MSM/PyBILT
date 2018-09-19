"""Define Gaussian function objects.

This module defines the Gaussian class and the GaussianRange class.

"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import numpy as np
from six.moves import range



class Gaussian(object):
    """A Gaussian function object.

    Attributes:
        mean (float): The mean of the Gaussian.
        std (float): The standard deviation of the Gaussian.

    """
    def __init__(self, mean,std):
        """Initialize a Gaussian function object.

        Args:
            mean (float): Set the mean of the Gaussian.
            std (float): Set the standard deviation of the Gaussian.
        """
        stdinv = 1.0/std
        normalc = stdinv*(1.0/np.sqrt(np.pi))
        self.sigma = std
        self.mean = mean
        self._normconst = normalc
        return


    def eval(self,x_in):
        """Return the Gaussian function evaluated at the input x value.

        Args:
            x_in (float): The x value to evaluate the function at.

        Returns:
            float: The function evaluation for the Gaussian.

        """
        stdinv = 1.0/self.sigma
        stdinvsq = stdinv**2
        normalc = self._normconst
        expon = -(x_in - self.mean)**2 * (0.5*stdinvsq)
        y = normalc * np.exp(expon)
        return y

    def reset_mean(self,new_mean):
        """Change the mean of the Gaussian function.

        Args:
            new_mean (float): The new mean of the Gaussian function.

        """
        self.mean = new_mean
        return


class GaussianRange(object):
    """Define a Gaussian function over a range.

     This object is used to define a Gaussian function over a defined
     finite range and store its values as evaluated at points evenly spaced
     over the range. The points can then for example  be used for integrating
     the Gaussian function over the range using numerical quadrature.

    Attributes:
        mean (float): The mean of the Gaussian.
        std (float): The standard deviation of the Gaussian.
        upper (float): The upper boundary of the range.
        lower (float): The lower boundary of the range.
        npoints (int): The number of points to evaluate in the range.

    """
    def __init__(self,in_range,mean,std,npoints=200):
        """Initialize the GaussianRange object.

        The GaussianRange stores the values of Gaussian function with the
        input mean and standard deviation evaluated at evenly spaced points
        in the specified x-value range.

        Args:
            in_range (tuple, list): Specify the endpoints for range, e.g.
                (x_start, x_end).
            mean (float): The mean of the Gaussian function.
            std (float): The standard deviation of the Gaussian function.
            npoints (Optional[int]): The number of x-value points to
                evaluate the Gaussian function for in the specified range (i.e.
                in_range).
        """
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
        self._normconst = normalc
        self.upper = in_range[1]
        self.lower = in_range[0]
        self._dx = x_p[1]-x_p[0]
        self.npoints = npoints
        return

    def get_values(self):
        """Return the x and y values for the Gaussian range function.

        Returns:
            tuple: The x and y values for the function, returned as (
            x_values, y_values).

        """
        return (self.x,self.y)

    def eval(self,x_in):
        """Return the Gaussian function evaluated at the input x value.

        Args:
            x_in (float): The x value to evaluate the function at.

        Returns:
            float: The function evaluation for the Gaussian.

        """
        stdinv = 1.0/self.sigma
        stdinvsq = stdinv**2
        normalc = self._normconst
        expon = -(x_in - self.mean)**2 * (0.5*stdinvsq)
        y = normalc * np.exp(expon)
        return y

    def integrate_range(self, lower, upper):
        """Returns the numerical integration of the Gaussian range.

        This function does a simple quadrature for the Gaussian function as
        evaluated on the range (or subset of the range) specified at
        initialization.

        Args:
            lower (float): The lower boundary for the integration.
            upper (float): The upper boundary for the integration.

        Returns:
            float: The numerical value of the Gaussian range integrated from
                lower to upper.

        Notes:
            This function does not thoroughly check the bounds, so if upper
            is less than lower the function will break.
        """
        if upper>self.upper:
            upper=self.upper
        if lower<self.lower:
            lower = self.lower

        i_l = int(np.floor((lower-self.lower)/self._dx))
        i_u = int(np.floor((upper-self.lower)/self._dx))
        #print "i_l ",i_l," i_u ",i_u
        total = 0.0
        for i in range(i_l,i_u):
            total+= self.y[i]*self._dx
        return total

    def sum_range(self, lower, upper):
        """Returns the over the Gaussian range.

        This function sums the Gaussian function at the points that were
        evaluated on the range (or subset of the range) specified at
        initialization.

        Args:
            lower (float): The lower boundary for the sum.
            upper (float): The upper boundary for the sum.

        Returns:
            float: The numerical value of the Gaussian range as summed from
                lower to upper.

        Notes:
            This function does not thoroughly check the bounds, so if upper
            is less than lower the function will break.
        """
        if upper>self.upper:
            upper=self.upper
        if lower<self.lower:
            lower = self.lower

        i_l = int(np.floor((lower-self.lower)/self._dx))
        i_u = int(np.floor((upper-self.lower)/self._dx))
        total = 0.0
        for i in range(i_l,i_u):
            total+= self.y[i]
        return total

    def normalize(self):
        """Normalizes (by area) the Gaussian function values over the range."""
        total = 0.0
        for i in range(0,self.npoints):
            total+=self.y[i]*self._dx
        for i in range(0,self.npoints):
            self.y[i]/=total
        return

    def reset_mean(self,new_mean):
        """Change the mean of the Gaussian function.

        Args:
            new_mean (float): The new mean of the Gaussian function.

        Notes:
            This function does not re-evaluate the Gaussian range and
            therefore only affects the output of the eval function.

        """
        self.mean = new_mean
        return
