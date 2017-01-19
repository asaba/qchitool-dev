'''
Created on 20/lug/2012

@author: asaba
'''
from datetime import datetime
import inspect

def p(string, skip=True):
    if not skip:
        frm = inspect.stack()[1]
        #mod = inspect.getmodule(frm[0])
        parent = inspect.stack()[1][3]
        print '%s [%s] %s' % (str(datetime.now()), parent, string)