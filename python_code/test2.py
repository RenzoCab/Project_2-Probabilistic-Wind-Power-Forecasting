
import warnings
import numpy as np
warnings.filterwarnings('error', '.*invalid value encountered.*',)

for attempt in range(5):
    try:
        warnings.warn("RuntimeWarning: invalid value encountered")
        r=1+6
    except:
        print('caught', attempt)
    else:
        print('gave up')


print('over it')

r
np.nan
