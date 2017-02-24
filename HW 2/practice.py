from mdsim import *

box_length = 1.0
# pset = ParticleSet(1,1)
# pset.change_pos(0,np.array([0.501,0.499,0.501]))
# print( pos_in_box(pset.pos(0),box_length) )
# pset.change_pos(0,np.array([1.000, -.501, 4.501]))
# print( pos_in_box(pset.pos(0),box_length) )

pset = ParticleSet(2,1)
pset.change_pos(0,np.array([ 0.499, 0.000, 0.000 ]) )
pset.change_pos(1,np.array([ -.499, 0.000, 0.000 ]) )
print distance(0,1,pset,box_length)

pset.change_pos(0,np.array([ 0.001, 0.000, 0.000 ]) )
pset.change_pos(1,np.array([ -.001, 0.000, 0.000 ]) )
print( distance(0,1,pset,box_length) )