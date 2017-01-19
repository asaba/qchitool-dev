'''
Created on 26/ago/2011

@author: asaba
'''

import datetime
import hashlib
import sys
from .. import performancetest as pt


def datatimenow(returnobject = False):
    if returnobject:
        return datetime.datetime.now()
    else:
        return datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

def returnmd5(tmp_string):
    md5_tmp = hashlib.md5()
    md5_tmp.update(tmp_string)
    md5 = md5_tmp.hexdigest()
    return md5

def openfile(filename):
    #This function open the file passed and return the related pointer
    try:
        f = open(filename,"r") 
    except IOError:
        print "No file found."
        sys.exit()
    return f

def val(string): #to add boolean evaluation
    try:
        return int(string)
    except:
        try:
            return float(string)
        except:
            return None
        
def handle_uploaded_file(file, session_path):
#    logging.debug("upload_here")
    if file:
        newfilename = session_path + file.name
        destination = open(newfilename, 'wb+')
        #destination = open('/tmp', 'wb+')
        for chunk in file.chunks():
            destination.write(chunk)
        destination.close()
        return newfilename
    else:
        return None

def create_xyz_file(filename, molecule, session_path):
    if file and molecule:
        destination = session_path + filename + ".xyz"
        molecule.write("xyz", destination)
    return session_path + file.name + ".xyz"

def create_cml_file(filename, molecule, session_path):
    if file and molecule:
        destination = session_path + filename + ".xml"
        molecule.write("cml", destination)
    return session_path + file.name + ".xml"

class RoundingError:
    def __init__(self, relativeerror=0.001, absoluterror = 1E-3):
        #Set initial values of thresholds
        self.coordrelativeerror = abs(relativeerror)
        self.absoluterror = abs(absoluterror)

    def isabsoluteequal(self, value1, value2):
        #Return true if values are similar by absolute error
        val1 = max(value1, value2)
        val2 = min(value1, value2)
        return self.equal(val1, val2, self.absoluterror)
            
    def isrelativeequal(self, value1, value2, coefficienterror=None):
        #return true if values are similar by relative error
        #It use the coefficienterror to calculate errors threasholds
        val1 = max(value1, value2)
        val2 = min(value1, value2)
        if coefficienterror:
            error1 = abs(coefficienterror*val1)
            error2 = abs(coefficienterror*val2)
        else:
            error1 = abs(self.coordrelativeerror*val1)
            error2 = abs(self.coordrelativeerror*val2)
        return self.equal(val1, val2, error1, error2)
    
    def equal(self, val1, val2, val1_error, val2_error=None):
        #val2<=val1
        result = True
        error = abs(val1 - val2)
        if abs(val1) >= val1_error:
            #val1>=e
            if abs(val2) < val2_error:
                #val2<e
                result = False
            else:
                #val2>=e
                #return |v1-v2|<e and |v1-v2|<f
                if val2_error:
                    result = error <= val1_error and error <= val2_error
                else:
                    result = error <= val1_error
        return result, error

        
        
    def coord3dequal(self, coord1, coord2):
            rng = len(coord1)
            result = True
            error = 0.0
            for index in range(rng):
                tmpresult, tmperror = self.isabsoluteequal(coord1[index], coord2[index])
                result = result and tmpresult
                error = error + tmperror
            return  result, error/rng
                
    def countmaxequal(self, vector, relative=False):
        count = 0
        tmpdict = {}
        indexequalvalues = []
        for index in range(len(list(vector))):
            found = False
            for key in tmpdict:
                if relative:
                    equal, error = self.isrelativeequal(vector[index], eval(key))
                else:
                    equal, error = self.isabsoluteequal(vector[index], eval(key))
                if equal:
                    found = True
                    tmpdict[key].append(index)
            if not found:
                tmpdict[str(vector[index])] = [index]
        for key in tmpdict:
            if len(tmpdict[key]) > count:
                count = len(tmpdict[key])
                indexequalvalues = tmpdict[key]
        return count, indexequalvalues