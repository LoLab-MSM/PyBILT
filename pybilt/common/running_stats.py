import numpy as np
#Running Statistics
class RunningStats(object):


    def __init__(self):
        self.n=0
        self.Mnold = self.Mnnew = self.Snold = self.Snnew = np.zeros(1)[0]

    def push(self, val):
        self.n += 1
        if self.n == 1:
            self.Mnold = np.array([val])[0]
            self.Snold = np.zeros(1)[0]
        else:
            n = np.array([float(self.n)])[0]
            self.Mnnew = self.Mnold + (val - self.Mnold)/(n);
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
            one = np.array([1.0])[0]
            n = np.array([float(self.n)])[0]
            vary = self.Snnew/(n-one)
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

class BlockAverager(object):

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
            min_points_in_block = points_per_block
        self._min_points_in_block = min_points_in_block
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
            self.push_single(datum)
        return

    def _get_running(self):
        means = []
        for block in self._blocks:
            if block.n >= self._min_points_in_block:
                means.append(block.mean())
        means = np.array(means)
        block_average = means.mean()
        std_err = means.std()/np.sqrt(len(means))
        return block_average, std_err

    def _get_np(self):
        means = []
        for block in self._blocks:
            if len(block) >= self._min_points_in_block:
                means.append(np.array(block).mean())
        means = np.array(means)
        block_average = means.mean()
        std_err = means.std()/np.sqrt(len(means))
        return block_average, std_err

    def get(self):
        """Return the block average and standard error.

        Returns:
            tuple: Returns a length two tuple with the block average and standard error estimates.
        """
        if self._store_data:
            return self._get_np()
        else:
            return self._get_running()


    def n_blocks(self):
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