
def crevasse_water_depth(x,y,t,thck,topg,*etc):

    
    
    depth = 0.0

    if (t >= 10.0):
        depth = 20.0
        
    if (t >= 20.0):
        depth = 0.0
        
    
    return depth 
