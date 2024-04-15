import random
import math
import sys
import shutil


def print_paraview (ncoords, mesh):
    
    num_pts=len(ncoords) # 3D
    print('number of nodes=',num_pts)
    dim_tet=len(mesh) # Tetrahedro
    print('number of elements=',dim_tet)
    a='viga'
    c='.vtu'
    nam=a+c
    
    # v_tetra should be a vector, defined previously
    #f=open('output_micro.vtu','w')
    f=open(nam,'w')
    f.write('<?xml version="1.0"?> \n')
    f.write('<VTKFile type="UnstructuredGrid"> \n')
    f.write('   <UnstructuredGrid> \n')
    f.write('      <Piece NumberOfPoints="')
    f.write(str(int(num_pts)))
    f.write('"  NumberOfCells="')
    f.write(str(int(dim_tet)))
    f.write('">  \n')
    f.write('          <Points>  \n')
    f.write('             <DataArray type="Float64" NumberOfComponents="3" format="ascii">   \n')

# Data array --> Coordinates

    for i in range(0,num_pts,1):

        f.write( str( ncoords[i][0]))
        f.write('   ')
        f.write(str( ncoords[i][1]))
        f.write('   ')
        f.write(str( ncoords[i][2]))
        f.write('\n')

    f.write('             </DataArray>   \n')
    f.write('          </Points>  \n')
    f.write('          <Cells>  \n')
    f.write('             <DataArray type="Int32" Name="connectivity" format="ascii">   \n')

# Data array --> Connectivity

    for i in range(0,dim_tet,1):
        f.write( str(mesh[i][0]))
        f.write('   ')
        f.write(str(mesh[i][1]))
        f.write('   ')
        f.write(str(mesh[i][2]))
        f.write('   ')
        f.write(str(mesh[i][3]))
        f.write('\n')
    
    f.write('              </DataArray>   \n')
    f.write('              <DataArray type="Int32" Name="offsets" format="ascii">   \n')

     # Data array --> Offset

    for i in range(0,dim_tet,1):
        aux_p=(i+1)*4
        f.write( str(aux_p))
        f.write('\n')

    f.write('              </DataArray>   \n')
    f.write('              <DataArray type="UInt8" Name="types" format="ascii">   \n')
    
    for i in range(0,dim_tet,1):
        f.write( str(10))
        f.write('\n')

    f.write('              </DataArray>   \n')
    f.write('          </Cells>  \n')
    f.write('       </Piece>  \n')
    f.write('   </UnstructuredGrid> \n')
    f.write('</VTKFile> \n')
    f.close()

    return








