
from cpython.list cimport  PyList_Append, PyList_GET_SIZE

def masserror(double target, double mass):
    cdef double me= ((target-mass)/target) * 1000000
    return me

def accuratemass(double target):
    cdef list results=[]
    cdef double[5] elmasses=[12.0, 1.007825032, 14.003074, 15.99491462, 35.9670807]
    cdef double mass1=0
    cdef double me=0
    #double elmasses[] = {12.0, 1.007825032, 14.003074, 15.99491462, 35.9670807}
    #cdef list elnumax=[int(target/elmasses[0]),int(target/elmasses[0])*4,int(target/elmasses[2]),int(target/elmasses[3]),int(target/elmasses[4])]
    cdef int[5] elnumax
    for i in range(5):
        elnumax[i]=int(target/elmasses[i])
    cdef int nc,nh,nn,no,ns
    cdef double m2=0.0
    for nc in range(elnumax[0]+1):
        for nh in range(elnumax[0]*4):
            m2 = (nc*elmasses[0])+(nh*elmasses[1])
            if m2>(target+1):
                break
            for nn in range(elnumax[2]+1):
                m2 = (nc*elmasses[0])+(nh*elmasses[1])+(nn*elmasses[2])
                if m2>(target+1):
                    break
                for no in range(elnumax[3]+1):
                    m2 = (nc*elmasses[0])+(nh*elmasses[1])+(nn*elmasses[2])+(no*elmasses[3])
                    if m2>(target+1):
                        break
                    for ns in range(elnumax[4]+1):
                        mass1 = (nc*elmasses[0])+(nh*elmasses[1])+(nn*elmasses[2])+(no*elmasses[3])+(ns*elmasses[4])
                        if mass1>target-1 and mass1<target+1:
                            #print("yay")
                            me=masserror(target,mass1)
                            #print([nc,nh,nn,no,ns])
                            if abs(me)<1:
                                results.append([nc,nh,nn,no,ns,me,mass1])
                                #print("yay")
    return results


def masserror2(double target, double mass,double atmass):
    cdef double me= ((target-mass)/atmass) * 1000000
    return me

def accuratemass2(double target,double atmass):
    cdef list results=[]
    cdef double[5] elmasses=[12.0, 1.007825032, 14.003074, 15.99491462, 35.9670807]
    cdef double mass1=0
    cdef double me=0
    #double elmasses[] = {12.0, 1.007825032, 14.003074, 15.99491462, 35.9670807}
    #cdef list elnumax=[int(target/elmasses[0]),int(target/elmasses[0])*4,int(target/elmasses[2]),int(target/elmasses[3]),int(target/elmasses[4])]
    cdef int[5] elnumax
    for i in range(5):
        elnumax[i]=int(target/elmasses[i])
    cdef int nc,nh,nn,no,ns
    cdef double m2=0.0
    for nc in range(elnumax[0]+1):
        for nh in range(nc*3):
            m2 = (nc*elmasses[0])+(nh*elmasses[1])
            if m2>(target+1):
                break
            for nv in [0,1]:
                nn=0
                m2 = (nc*elmasses[0])+(nh*elmasses[1])+(nn*elmasses[2])
                if m2>(target+1):
                    break
                for no in range(nc+5):
                    m2 = (nc*elmasses[0])+(nh*elmasses[1])+(nn*elmasses[2])+(no*elmasses[3])
                    if m2>(target+1):
                        break
                    #for ns in range(elnumax[4]+1):
                    mass1 = (nc*elmasses[0])+(nh*elmasses[1])+(nn*elmasses[2])+(no*elmasses[3])
                    if mass1>target-1 and mass1<target+1:
                        #print(mass1)
                        me=masserror2(target,mass1,atmass)
                        #print([nc,nh,nn,no,ns])
                        if abs(me)<1:
                            results.append([nc,nh,nn,no,me,mass1])
                            #print("yay")
    return results
