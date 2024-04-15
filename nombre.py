'''import numpy as np
import pandas as pd 
import meshpy.tet as tet
from scipy.sparse import lil_matrix

# 1. Crear un array de 10x10
dominio = np.zeros((10, 10, 10))
print(dominio) #10 matrices de 10x10

# 2. Crear un CVS
df = pd.DataFrame({'x': [0, 1], 'y': [0, 1], 'z': [0, 1], 'pressure': [100, 150]})
df.to_csv('output.csv', index=False)

# 3. Crear un tensor de deformación
def tensor_deformacion(desplazamientos):
    # Implementa el cálculo del gradiente de desplazamientos
    # Este es un ejemplo simplificado
    return np.gradient(desplazamientos)

# 4. Generación de una malla de tetraedros
points, facets = [ ... ], [ ... ]  # Define puntos y caras
mesh_info = tet.MeshInfo()
mesh_info.set_points(points)
mesh_info.set_facets(facets)
mesh = tet.build(mesh_info)

# 5. Implementación de Función de forma
def funcion_forma(punto, nodos):
    # Implementa la función de forma para un punto y nodos
    return np.ones(len(nodos))

# 6. Ensamblaje de matrices
K_global = lil_matrix((nodos_totales, nodos_totales))
# Suma las matrices locales de rigidez al K_global en los índices correctos'''

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import pyvtk

# PARTE 1
dominio = np.zeros((10, 10, 10))

# Generar coordenadas para el dominio estructural
x = np.linspace(0, 1, 10)
y = np.linspace(0, 1, 10)
z = np.linspace(0, 1, 10)

# Crear una cuadrícula tridimensional
X, Y, Z = np.meshgrid(x, y, z)

# Visualizar el dominio estructural
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(X, Y, Z)
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
plt.show()

# PARTE 2
def export_to_paraview(nombre, presiones, desplazamientos):
    """
    Exporta los resultados de un análisis estructural a un formato compatible con Paraview.

    Args:
    - nombre: Nombre del archivo de salida.
    - presiones: Array de presiones.
    - desplazamientos: Array de desplazamientos.
    """
    # crear una malla sin conexiones usando pyvtk
    puntos = np.column_stack((np.arange(len(presiones)), np.zeros(len(presiones)), np.zeros(len(presiones))))
    mesh = pyvtk.UnstructuredGrid(puntos, point_data={"Presiones": presiones, "Desplazamientos": desplazamientos})

    #Lo pasamos a VTK
    mesh.tofile(nombre)

'''# Ejemplo de uso
presiones = np.array([1.5, 2.0, 1.8, 2.2])
desplazamientos = np.array([0.1, 0.2, 0.15, 0.18])
export_to_paraview('results.txt', presiones, desplazamientos)
'''

def calcular_tensor_deformaciones(desplazamientos):
    """
    Calcula el tensor de deformaciones para un elemento finito tetraédrico.

    Args:
    - desplazamientos: Matriz de desplazamientos nodales. Debe ser una matriz 4x3, donde cada fila representa
                     los desplazamientos nodales en un nodo (x, y, z).

    Returns:
    - strain_tensor: Tensor de deformaciones calculado.
    """
    # Definir las derivadas de las funciones de forma para un tetraedro lineal
    #Es una matriz 3x4 porque en un elemento tetraédrico lineal, tenemos 4 nodos, y cada nodo tiene tres grados de libertad
    B = np.array([
        [-1, 1, 0, 0],
        [-1, 0, 1, 0],
        [-1, 0, 0, 1]
    ])
    
    # Calcular el tensor de deformaciones
    strain_tensor = np.dot(B, desplazamientos)
    
    return strain_tensor

# Ejemplo de uso
desplazamientos = np.array([
    [0.1, 0.2, 0.3],
    [0.2, 0.3, 0.4],
    [0.3, 0.4, 0.5],
    [0.4, 0.5, 0.6]
])

'''strain_tensor = calcular_tensor_deformaciones(desplazamientos)
print("Tensor de deformaciones:")
print(strain_tensor)'''

def generar_mallado_tetraedros(geometria):
    """
    Genera un mallado de tetraedros a partir de una geometría proporcionada.

    Args:
    - geometria: Datos de la geometría. Esto puede ser una lista de coordenadas de nodos y una lista de
                 índices de nodos que forman los tetraedros.

    Returns:
    - nodos: Coordenadas de los nodos del mallado.
    - tetraedros: Índices de nodos que forman los tetraedros del mallado.
    """
    # Supongamos que la geometría se proporciona como una lista de coordenadas de nodos y una lista de índices de nodos
    nodos = geometria['nodos']
    tetraedros = geometria['tetraedros']
    
    # Aquí podríamos realizar algún proceso de generación de elementos, como refinamiento o ajuste de la malla
    
    return nodos, tetraedros

# Ejemplo de uso
geometria_ejemplo = {
    'nodos': np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]]),
    'tetraedros': np.array([[0, 1, 2, 3]])
}

nodos_mallado, tetraedros_mallado = generar_mallado_tetraedros(geometria_ejemplo)
print("Nodos del mallado:")
print(nodos_mallado)
print("Tetraedros del mallado:")
print(tetraedros_mallado)

