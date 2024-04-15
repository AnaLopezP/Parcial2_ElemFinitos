import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import pyvtk
import meshpy.tet as tet

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

# PARTE 3
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

# PARTE 4
def generar_mallado_tetraedros(geometria):
    """
    Genera un mallado de tetraedros usando meshpy.

    Args:
    - geometria: Datos de la geometría. Esto puede ser una lista de coordenadas de nodos y una lista de
                 índices de nodos que forman los tetraedros.

    Returns:
    - nodos: Coordenadas de los nodos del mallado.
    - tetraedros: Índices de nodos que forman los tetraedros del mallado.
    """
   # Definir puntos y caras (facets)
    points = np.array([
        [0, 0, 0],
        [1, 0, 0],
        [0, 1, 0],
        [0, 0, 1]
    ])
    
    facets = [
        [0, 1, 2],
        [0, 1, 3],
        [1, 2, 3],
        [0, 2, 3]
    ]
    
    # Crear objeto MeshInfo y definir puntos y facetas
    mesh_info = tet.MeshInfo()
    mesh_info.set_points(points)
    mesh_info.set_facets(facets)
    
    # Construir la malla
    mesh = tet.build(mesh_info)
    
    # Obtener puntos y elementos de la malla
    points = np.array(mesh.points)
    elements = np.array(mesh.elements)
    
    return points, elements

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

# PARTE 5
def funcion_de_forma_tetraedro(xi, eta, zeta):
    """
    Calcula los valores de la función de forma para un elemento tetraédrico lineal.

    Args:
    - xi, eta, zeta: Coordenadas en el espacio paramétrico del tetraedro.

    Returns:
    - N: Vector de valores de la función de forma.
    """
    N = np.array([
        1 - xi - eta - zeta,
        xi,
        eta,
        zeta
    ])
    
    return N

# Ejemplo de uso
xi = 0.2
eta = 0.3
zeta = 0.4
N = funcion_de_forma_tetraedro(xi, eta, zeta)
print("Valores de la función de forma:", N)
