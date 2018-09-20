"""Running stats module.

This module defines the RunningStats and BlockAverager classes, as well as the
gen_running_average function.

"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from builtins import object
import numpy as np
from six.moves import range

# Running Statistics
class RunningStats(object):
    """A RunningStats object.

    The RunningStats object keeps running statistics for a single
    value/quantity.

    Attributes:
        n (int): The number of points that have pushed to the running
        average.
    """
    def __init__(self):
        """Initialize the RunningStats object.
        """
        self.n=0
        self._Mnold = self._Mnnew = self._Snold = self._Snnew = np.zeros(1)[0]


    def push(self, val):
        """Push a new value to the running average.

        Args:
            val (float): The value to be added to the running average.

        Returns:

        """
        self.n += 1
        if self.n == 1:
            self._Mnold = np.array([val])[0]
            self._Snold = np.zeros(1)[0]
        else:
            n = np.array([float(self.n)])[0]
            self._Mnnew = self._Mnold + (val - self._Mnold)/(n);
            self._Snnew = self._Snold + (val - self._Mnold)*(val-self._Mnnew);
            self._Mnold = self._Mnnew;
            self._Snold = self._Snnew;


    def mean(self):
        """Return the current mean."""
        if self.n == 1:
            return self._Mnold
        elif self.n > 1:
            return self._Mnnew
        else:
            return 0.0


    def variance(self):
        """Returun the current variance."""
        if self.n > 1:
            one = np.array([1.0])[0]
            n = np.array([float(self.n)])[0]
            vary = self._Snnew/(n-one)
            return vary
        else:
            return 0.0


    def deviation(self):
        """Return the current standard deviation."""
        # dev = math.sqrt(self.Variance())
        dev = np.sqrt(self.variance())
        return dev


    def reset(self):
        """Reset the running average."""
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
    for i in range(nele):
        averager.push(onednparray[i])
        run_avg = averager.mean()
        run_dev = averager.deviation()
        # print run_avg, run_dev, averager.mean(), onednparray[i]
        output[i,0] = run_avg
        output[i,1] = run_dev
    return output

class BlockAverager(object):
    """An object that keeps track of points for block averaging.

    Attributes:
        n_blocks (int): The current number of active blocks.

    """

    def __init__(self, points_per_block=1000, min_points_in_block=500, store_data=False):
        """Init a the BlockAverager

        Args:
            points_per_block (int, Optional): The number of points to assign to a block before initiating a new block.
                Default: 1000
            min_points_in_block (int, Optional): The minimum number of points that a block (typically the last block)
                can have and still be included in computing the final block average and standard error estimates. This
                value should be <= points_per_block. Default: 500
        """
        self._store_data = store_data
        self._blocks = [RunningStats()]
        if store_data:
            self._blocks = [[]]
        self.n_blocks = 1
        self._points_per_block = points_per_block
        if min_points_in_block > points_per_block:
            self._min_points_in_block = points_per_block-1
        else:
            self._min_points_in_block = min_points_in_block
        #print "points_per_block ",self._points_per_block, " min_p ",self._min_points_in_block
        return

    def _add_block(self):
        """Append a new block."""
        if self._store_data:
            self._blocks.append([])
        else:
            self._blocks.append(RunningStats())
        self.n_blocks+=1
        return

    def _check_add_block(self):
        """Check whether to add a new block and do so if the condition is met."""
        block_i = self.n_blocks - 1
        if self._store_data:
            if len(self._blocks[block_i]) >= self._points_per_block:
                self._add_block()
        else:
            if self._blocks[block_i].n >= self._points_per_block:
                self._add_block()
        return

    def push_single(self, datum):
        """Push a single data point (datum) into the block averager.

        Args:
            datum (float): The value to add to the block averaging.

        """
        block_i = self.n_blocks-1
        #print "pushing datum ",datum
        if self._store_data:
            self._blocks[block_i].append(datum)
        else:
            self._blocks[block_i].push(datum)
        self._check_add_block()
        return

    def push_container(self, data):
        """Push a container (array or array like) of data points to the block averaging.

        Args:
            data (array like): The container (list, tuple, np.array, etc.) of data points to add to the block averaging.

        """
        for datum in data:
            #print(datum)
            self.push_single(datum)
        return

    def _get_running(self):
        """Get the block average quantities from interanl RunningStats
        objects.
        """
        means = []
        for block in self._blocks:
            #print "block.n ",block.n, " min_p ",self._min_points_in_block
            if block.n >= self._min_points_in_block:
                means.append(block.mean())
        means = np.array(means)
        if len(means) > 1:
            block_average = means.mean()
            std_err = means.std()/np.sqrt(len(means))
        elif len(means) == 1:
            block_average = means[0]
            std_err = 0.0
        else:
            block_average = 0.0
            std_err = 0.0
        return block_average, std_err

    def _get_np(self):
        """Get the block average quantities from internally stored numpy
        arrays.
        """
        means = []
        for block in self._blocks:
            if len(block) >= self._min_points_in_block:
                means.append(np.array(block).mean())
        means = np.array(means)
        if len(means) > 1:
            block_average = means.mean()
            std_err = means.std()/np.sqrt(len(means))
        elif len(means) == 1:
            block_average = means[0]
            std_err = 0.0
        else:
            block_average = 0.0
            std_err = 0.0

        return block_average, std_err

    def get(self):
        """Return the block average and standard error.

        Returns:
            tuple: Returns a length two tuple with the block average and standard error estimates.
        """
        #print(self._blocks)
        if self._store_data:
            return self._get_np()
        else:
            return self._get_running()


    def number_of_blocks(self):
        """Return the current number of blocks.

        Returns:
            int : The number of blocks.
        """
        return self.n_blocks

    def points_per_block(self):
        """Return information about the points per block.

        Returns:
            tuple: A three element tuple containing the setting for points per block, the setting for minimum points
                per block, and the number of points in the last block.
        """
        if self._store_data:
            return self._points_per_block, self._min_points_in_block, len(self._blocks[self.n_blocks-1])
        else:
            return self._points_per_block, self._min_points_in_block, self._blocks[self.n_blocks - 1].n

def binned_average(data, positions, n_bins=25, position_range=None):
    """Compute averages over a quantized range of histogram like bins.

    Args:
        data (np.array): A 1d numpy array of values.
        positions (np.array): A 1d numpy array of positions corresponding to
            the values in data. These are used to assign the values to the
            histogram like bins for averaging.
        n_bins (Optional[int]): Set the target number of bins to quantize the
            position_range up into. Defaults to 25
        position_range (Optional[tuple]): A two element tuple containing the
            lower and upper range to bin the postions over; i.e.
            (position_lower, postion_upper). Defaults to None, which uses
            positions.min() and positions.max().
    Returns:
        tuple: returns a tuple with two numpy arrays of form (bins, averages)

    Notes:
        The function automatically filters out bins that have a zero count,
        so the final value of the number of bins and values will be
        len(bins) <= n_bins.
    """
    lower = None
    upper = None

    if position_range is not None:

        lower = position_range[0]
        upper = position_range[1]
    else:
        lower = positions.min()
        upper = positions.max()

    edges = np.linspace(lower, upper, num=n_bins+1, endpoint=True)
    bins = np.linspace(lower, upper, num=n_bins, endpoint=False)
    counts = (np.zeros(len(bins))).astype(np.int64)
    sums = np.zeros(len(bins))
    n_data = len(data)
    # Loop over the data points.
    for i in range(n_data):
        c_val = data[i]
        pos = positions[i]
        bin_index = None
        # Select which bin (if any) the value corresponds to.
        for j in range(1, len(bins)+1):
            if (pos >= edges[j-1]) and (pos <= edges[j]):
                bin_index = j - 1
                break
        if bin_index is not None:
            counts[bin_index] += 1
            sums[bin_index] += c_val

    # Filter out the bins that had zero entries.
    keep_bins = []
    for i in range(len(counts)):
        if counts[i] > 0:
            keep_bins.append(i)
    # Return the filtered bins and averages (i.e. without NaNs).
    bins = bins[keep_bins]
    counts = counts[keep_bins]
    sums = sums[keep_bins]
    averages = sums / counts
    return bins, averages
