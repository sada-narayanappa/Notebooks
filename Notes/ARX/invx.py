
from sklearn.linear_model import LinearRegression
import numpy as np;

def logd(debug = True, *args):
    if not debug: return
    for a in args:
        print(a, sep='', end=' ')

log.getLogger().setLevel(log.INFO)