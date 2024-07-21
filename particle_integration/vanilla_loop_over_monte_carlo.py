import vanilla 
import tstrippy
import sys 

def perform_integration(GCname,montecarlokey,internal_dynamics,GCorbits_potential,MWpotential,NP):


    staticgalaxy,integrationparameters,initialkinematics,inithostperturber = \
        vanilla.load_arguments(GCname,montecarlokey,internal_dynamics,GCorbits_potential,MWpotential,NP)
        
    
    # integrate the particle
    integrator = tstrippy.integrator
    integrator.setstaticgalaxy(*staticgalaxy)
    integrator.setintegrationparameters(*integrationparameters)
    integrator.setinitialkinematics(*initialkinematics)
    integrator.inithostperturber(*inithostperturber)
    integrator.leapfrogtofinalpositions()
    xf= integrator.xf
    yf= integrator.yf
    zf= integrator.zf
    vxf= integrator.vxf
    vyf= integrator.vyf
    vzf= integrator.vzf
    tesc= integrator.tesc
    integrator.deallocate()
    return xf,yf,zf,vxf,vyf,vzf,tesc
    
    
if __name__ == "__main__" : 
    
    
    GCname              = sys.argv[1] 
    montecarlokey       = sys.argv[2]
    internal_dynamics   =   "isotropic-plummer"
    GCorbits_potential  =   "pouliasis2017pii"
    MWpotential         =   "pouliasis2017pii"
    NP                  =   int(1e2)    
