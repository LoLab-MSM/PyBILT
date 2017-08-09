import numpy as np
from pybilt.common.running_stats import BlockAverager

def test_block_averager():
    block_0_data = np.array([1.2, 1.7, 2.3, 1.4, 2.0])
    block_1_data = np.array([3.3, 0.8, 1.6, 1.8, 2.1])
    block_2_data = np.array([0.4, 2.1, 1.3, 1.1, 1.95])

    block_0_mean = block_0_data.mean()
    block_1_mean = block_1_data.mean()
    block_2_mean = block_2_data.mean()

    block_average = np.array([block_0_mean, block_1_mean, block_2_mean]).mean()
    std_err = np.array([block_0_mean, block_1_mean, block_2_mean]).std()/np.sqrt(3)
    #print np.array([block_0_mean, block_1_mean, block_2_mean]).std()
    print("Manual Block averaging: {} +- {}".format(block_average, std_err))

    block_averager = BlockAverager(points_per_block=5)
    block_averager.push_container(block_0_data)
    block_averager.push_container(block_1_data)
    block_averager.push_container(block_2_data)
    block_average, std_err = block_averager.get()
    print("BlockAverager block averaging: {} +- {}".format(block_average, std_err))
    block_averager = BlockAverager(points_per_block=5, store_data=True)
    block_averager.push_container(block_0_data)
    block_averager.push_container(block_1_data)
    block_averager.push_container(block_2_data)
    block_average, std_err = block_averager.get()
    print("BlockAverager block averaging (store_data=True): {} +- {}".format(block_average, std_err))
    return

if __name__ == '__main__':
    test_block_averager()