import random
import math
import sys
import shutil
from print_paraview import*

'''
PARAMETERS
'''

################## MECHANICAL PROPERTIES

E=28300.0
pois=0.2
density=2400 ## kg/m3

################# MESH PROPERTIES

sizex=1000
sizey=150
sizez=50
ndivx=20
ndivy=5
ndivz=3
jumpx=0



'''
MESHER

Este codigo genera una malla tetrahedrica sobre un elemento prismatico
'''

vec_nodes_0=[]
gapx=sizex/ndivx
gapy=sizey/ndivy
gapz=sizez/ndivz

for i in range(0,ndivx+1,1):
    for j in range(0,ndivy+1,1):
        for k in range(0,ndivz+1,1):
            vec_nodes_0.append([i*gapx+jumpx-sizex/2,j*gapy,k*gapz])



mesh=[]

for i in range(0,ndivx,1):
        
    for j in range(0,ndivy,1):
     
        for k in range(0,ndivz,1):
            
            # numbering: i+j*nx+k*nx*ny
            n1=[i,j,k]
            n2=[i+1,j,k]
            n4=[i,j+1,k]
            n3=[i+1,j+1,k]
            n5=[i,j,k+1]
            n6=[i+1,j,k+1]
            n8=[i,j+1,k+1]
            n7=[i+1,j+1,k+1]

            n_loc=[n1,n2,n3,n4,n5,n6,n7,n8]
            n_glob=[]
            for l in range(0, len(n_loc),1):
                #n_glob.append(n_loc[l][0]*(ndivx+1)*(ndivy+1)+n_loc[l][1]*(ndivz+1)+n_loc[l][2])
                n_glob.append(n_loc[l][0]*(ndivz+1)*(ndivy+1)+n_loc[l][1]*(ndivz+1)+n_loc[l][2])

            vec_aux=[[0,1,3,4],[1,2,3,6],[4,7,6,3],[4,6,5,1],[4,6,1,3]]
            for l in range(0,len(vec_aux),1):
                mesh.append([n_glob[vec_aux[l][0]], n_glob[vec_aux[l][1]], n_glob[vec_aux[l][2]], n_glob[vec_aux[l][3]]])


npar=len(vec_nodes_0)

v_tetra=[]
for i in range(0,int((npar)),1):    
        
    for k in range(0,3,1):
        v_tetra.append(vec_nodes_0[i][k])

fini= open('nodes.txt','w')
for i in range(0,len(vec_nodes_0),1):

    for j in range(0,len(vec_nodes_0[i]),1):
        fini.write(str((vec_nodes_0[i][j])))
        fini.write('  ')
    if i < len(vec_nodes_0)-1:
        fini.write('\n')
fini.close()

fini= open('mesh.txt','w')
fmech= open('mechanical_mat.txt','w')
for i in range(0,len(mesh),1):

    for j in range(0,len(mesh[i]),1):
        fini.write(str((mesh[i][j])))
        fini.write('  ')
    if i < len(mesh)-1:
        fini.write('\n')

    fmech.write(str( E  ))
    fmech.write('  ')
    fmech.write(str( pois  ))
    fmech.write('  ')
    fmech.write(str( density  ))
    fmech.write('\n')
    
fini.close()
fmech.close()

print_paraview (vec_nodes_0, mesh)


print('end of the function')







