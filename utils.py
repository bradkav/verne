import numpy as np

def sciformat(x,decimals=0):
    n = np.floor(np.log10(x))
    m = x*(10.0**(-n))
    str1 = ""
    if (m > 1.001):
        str1 = r'$%d \times ' % (m,)
        str2 = r'10^{%d}$' % n
    else:
        str2 = r'$10^{%d}$' % n
    return str1 + str2
    