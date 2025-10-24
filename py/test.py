# from scan import Linear_scan

# sca = Linear_scan(-1, 1, 5, -2, 2, 5, -3, 3, 5)
# print(sca.get_scan_xyz())

import numpy as np

# import test_np

from ..build.test_np import test_np

a = np.array([[1, 2, 3], [4, 5, 6]])

b = test_np.value_double(a)

print(b)
