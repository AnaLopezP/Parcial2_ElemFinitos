import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import pyvtk
import meshpy.tet as tet
from scipy.sparse import lil_matrix
import vtk

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


# PARTE 6
def ensamblar_matriz_rigidez_global(nodos, tetraedros, propiedades):
    """
    Ensambla la matriz de rigidez global para una estructura mallada con tetraedros.

    Args:
    - nodos: Coordenadas de los nodos del mallado.
    - tetraedros: Índices de nodos que forman los tetraedros del mallado.
    - propiedades: Propiedades del material y geometría de los tetraedros.

    Returns:
    - matriz_rigidez_global: Matriz de rigidez global ensamblada.
    """
    # Inicializar matriz de rigidez global
    num_nodos = len(nodos)
    matriz_rigidez_global = np.zeros((3*num_nodos, 3*num_nodos))
    
    # Iterar sobre los tetraedros y ensamblar la matriz de rigidez global
    for tetraedro in tetraedros:
        # Calcular matriz de rigidez local para el tetraedro actual
        matriz_rigidez_local = calcular_matriz_rigidez_local(tetraedro, nodos, propiedades)
        
        # Ensamblar matriz de rigidez local en la matriz de rigidez global
        ensamblar_matriz_local(matriz_rigidez_global, matriz_rigidez_local, tetraedro)
    
    return matriz_rigidez_global

def calcular_matriz_rigidez_local(tetraedro, nodos, propiedades):
    """
    Calcula la matriz de rigidez local para un tetraedro.

    Args:
    - tetraedro: Índices de nodos que forman el tetraedro.
    - nodos: Coordenadas de los nodos del mallado.
    - propiedades: Propiedades del material y geometría del tetraedro.

    Returns:
    - matriz_rigidez_local: Matriz de rigidez local del tetraedro.
    """
    # Propiedades del material
    E = propiedades['E'] # Módulo de elasticidad
    nu = propiedades['nu'] # Coeficiente de Poisson
    
    # Coordenadas de los nodos del tetraedro
    x0, y0, z0 = nodos[tetraedro[0]]
    x1, y1, z1 = nodos[tetraedro[1]]
    x2, y2, z2 = nodos[tetraedro[2]]
    x3, y3, z3 = nodos[tetraedro[3]]
    
    # Calculo de las derivadas de las funciones de forma
    dN_dxi = np.array([-1, 1, 0, 0])
    dN_deta = np.array([-1, 0, 1, 0])
    dN_dzeta = np.array([-1, 0, 0, 1])
    
    # Jacobiano de la transformación
    J = np.array([
        [x1 - x0, x2 - x0, x3 - x0],
        [y1 - y0, y2 - y0, y3 - y0],
        [z1 - z0, z2 - z0, z3 - z0]
    ])
    detJ = np.linalg.det(J)
   # Matriz de gradiente de deformaciones B
    B = np.array([
        [dN_dxi[0], 0, 0, dN_dxi[1], 0, 0, dN_dxi[2], 0, 0, dN_dxi[3], 0, 0],
        [0, dN_deta[0], 0, 0, dN_deta[1], 0, 0, dN_deta[2], 0, 0, dN_deta[3], 0],
        [0, 0, dN_dzeta[0], 0, 0, dN_dzeta[1], 0, 0, dN_dzeta[2], 0, 0, dN_dzeta[3]],
        [dN_deta[0], dN_dxi[0], 0, dN_deta[1], dN_dxi[1], 0, dN_deta[2], dN_dxi[2], 0, dN_deta[3], dN_dxi[3], 0],
        [0, dN_dzeta[0], dN_deta[0], 0, dN_dzeta[1], dN_deta[1], 0, dN_dzeta[2], dN_deta[2], 0, dN_dzeta[3], dN_deta[3]],
        [dN_dzeta[0], 0, dN_dxi[0], dN_dzeta[1], 0, dN_dxi[1], dN_dzeta[2], 0, dN_dxi[2], dN_dzeta[3], 0, dN_dxi[3]]
    ]) / detJ

    # Matriz de elasticidad
    factor = E / ((1 + nu) * (1 - 2 * nu))
    C = factor * np.array([
        [(1 - nu), nu, nu, 0, 0, 0],
        [nu, (1 - nu), nu, 0, 0, 0],
        [nu, nu, (1 - nu), 0, 0, 0],
        [0, 0, 0, (1 - 2 * nu) / 2, 0, 0],
        [0, 0, 0, 0, (1 - 2 * nu) / 2, 0],
        [0, 0, 0, 0, 0, (1 - 2 * nu) / 2]
    ])

    # Matriz de rigidez local
    matriz_rigidez_local = np.dot(np.dot(B.T, C), B) * detJ

    return matriz_rigidez_local

  
def ensamblar_matriz_local(matriz_global, matriz_local, tetraedro):
    """
    Ensambla la matriz de rigidez local en la matriz de rigidez global.

    Args:
    - matriz_global: Matriz de rigidez global.
    - matriz_local: Matriz de rigidez local.
    - tetraedro: Índices de nodos que forman el tetraedro.
    """
    # Tamaño de la matriz global
    num_nodos = matriz_global.shape[0]

    # Convertir la matriz local a una matriz dispersa lil_matrix
    matriz_local_sparse = lil_matrix(matriz_local)

    # Índices de los nodos del tetraedro
    i, j, k, l = tetraedro

    # Ensamblar la matriz local en la matriz global
    for m, n in enumerate([i, j, k, l]):
        for p, q in enumerate([i, j, k, l]):
            matriz_global[3*n:3*(n+1), 3*q:3*(q+1)] += matriz_local_sparse[3*m:3*(m+1), 3*p:3*(p+1)]



# PARTE 7
def resolver_sistema_ecuaciones(matriz_rigidez_global, fuerzas):
    """
    Resuelve el sistema de ecuaciones de equilibrio para un análisis estructural.

    Args:
    - matriz_rigidez_global: Matriz de rigidez global del sistema.
    - fuerzas: Vector de fuerzas aplicadas al sistema.

    Returns:
    - desplazamientos: Vector de desplazamientos resultantes.
    """
    # Resolver el sistema de ecuaciones utilizando un método de solución directa
    desplazamientos = np.linalg.solve(matriz_rigidez_global, fuerzas)
    
    return desplazamientos


# PARTE 8
# Suponiendo que ya tenemos los desplazamientos y la matriz de rigidez global

import numpy as np

def calcular_tensiones_deformaciones(matriz_rigidez_global, desplazamientos, nodos, tetraedros, propiedades):
    """
    Calcula las tensiones y deformaciones en la estructura analizada.

    Args:
    - matriz_rigidez_global: Matriz de rigidez global del sistema.
    - desplazamientos: Vector de desplazamientos resultantes.

    Returns:
    - tensiones: Vector de tensiones en los elementos de la estructura.
    - deformaciones: Vector de deformaciones en los elementos de la estructura.
    """
    # Inicializar vectores de tensiones y deformaciones
    tensiones = []
    deformaciones = []

    # Iterar sobre los tetraedros para calcular las tensiones y deformaciones en cada elemento
    for tetraedro in tetraedros:
        # Obtener las coordenadas de los nodos del tetraedro
        coordenadas = np.array([nodos[i] for i in tetraedro])

        # Calcular la matriz de deformación para el tetraedro actual
        matriz_deformacion = calcular_matriz_deformacion(coordenadas)

        # Calcular los desplazamientos nodales del tetraedro
        desplazamientos_nodales = np.array([desplazamientos[3*i:3*(i+1)] for i in tetraedro])

        # Calcular los desplazamientos nodales del tetraedro en coordenadas deformadas
        desplazamientos_deformados = np.dot(matriz_deformacion, desplazamientos_nodales.flatten())

        # Calcular las deformaciones en el tetraedro
        deformaciones_elemento = calcular_deformaciones_elemento(matriz_deformacion, desplazamientos_nodales)

        # Calcular las tensiones en el tetraedro
        tensiones_elemento = calcular_tensiones_elemento(deformaciones_elemento, propiedades)

        # Agregar las tensiones y deformaciones del elemento al vector global
        tensiones.append(tensiones_elemento)
        deformaciones.append(deformaciones_elemento)

    return tensiones, deformaciones

def calcular_matriz_deformacion(coordenadas):
    """
    Calcula la matriz de deformación para un tetraedro.

    Args:
    - coordenadas: Coordenadas de los nodos del tetraedro.

    Returns:
    - matriz_deformacion: Matriz de deformación del tetraedro.
    """
    # Calcular la matriz de deformación utilizando las coordenadas de los nodos
    x1, y1, z1 = coordenadas[0]
    x2, y2, z2 = coordenadas[1]
    x3, y3, z3 = coordenadas[2]
    x4, y4, z4 = coordenadas[3]
    
    matriz_deformacion = np.array([
        [x2 - x1, x3 - x1, x4 - x1],
        [y2 - y1, y3 - y1, y4 - y1],
        [z2 - z1, z3 - z1, z4 - z1]
    ])
    
    return matriz_deformacion


def calcular_deformaciones_elemento(matriz_deformacion, desplazamientos_nodales):
    """
    Calcula las deformaciones en un elemento finito.

    Args:
    - matriz_deformacion: Matriz de deformación del elemento.
    - desplazamientos_nodales: Desplazamientos nodales del elemento.

    Returns:
    - deformaciones_elemento: Deformaciones en el elemento.
    """
    # Calcular las deformaciones en el elemento utilizando la matriz de deformación y los desplazamientos nodales
    deformaciones_elemento = np.dot(matriz_deformacion, desplazamientos_nodales)
    
    return deformaciones_elemento


def calcular_tensiones_elemento(deformaciones_elemento, propiedades):
    """
    Calcula las tensiones en un elemento finito.

    Args:
    - deformaciones_elemento: Deformaciones en el elemento.
    - propiedades: Propiedades del material del elemento.

    Returns:
    - tensiones_elemento: Tensiones en el elemento.
    """
    # Calcular las tensiones en el elemento utilizando las deformaciones y las propiedades del material
    E = propiedades['E']  # Módulo de elasticidad
    nu = propiedades['nu']  # Coeficiente de Poisson
    
    # Calcular la matriz de elasticidad para un material isotrópico
    C = (E / ((1 + nu) * (1 - 2 * nu))) * np.array([
        [1 - nu, nu, nu, 0, 0, 0],
        [nu, 1 - nu, nu, 0, 0, 0],
        [nu, nu, 1 - nu, 0, 0, 0],
        [0, 0, 0, (1 - 2 * nu) / 2, 0, 0],
        [0, 0, 0, 0, (1 - 2 * nu) / 2, 0],
        [0, 0, 0, 0, 0, (1 - 2 * nu) / 2]
    ])
    
    # Calcular las tensiones en el elemento
    tensiones_elemento = np.dot(C, deformaciones_elemento)
    
    return tensiones_elemento



def visualizar_resultados(tensiones, deformaciones, nodos, tetraedros, nombre_archivo):
    """
    Visualiza las tensiones y deformaciones en la estructura utilizando Paraview u otra herramienta de visualización.

    Args:
    - tensiones: Vector de tensiones en los elementos de la estructura.
    - deformaciones: Vector de deformaciones en los elementos de la estructura.
    """
    # Crear un objeto vtkUnstructuredGrid para almacenar la geometría y los datos de los tetraedros
    grid = vtk.vtkUnstructuredGrid()

    # Crear un objeto vtkPoints para almacenar las coordenadas de los nodos
    points = vtk.vtkPoints()

    # Agregar los nodos al objeto vtkPoints
    for nodo in nodos:
        points.InsertNextPoint(nodo)

    # Asignar los puntos al vtkUnstructuredGrid
    grid.SetPoints(points)

    # Crear un objeto vtkIntArray para almacenar las tensiones en los elementos
    tensiones_array = vtk.vtkDoubleArray()
    tensiones_array.SetNumberOfComponents(1)
    tensiones_array.SetName("Tensiones")

    # Agregar las tensiones al vtkIntArray
    for tension in tensiones:
        tensiones_array.InsertNextValue(tension)

    # Añadir el array de tensiones al vtkUnstructuredGrid
    grid.GetCellData().AddArray(tensiones_array)

    # Crear un objeto vtkIntArray para almacenar las deformaciones en los elementos
    deformaciones_array = vtk.vtkDoubleArray()
    deformaciones_array.SetNumberOfComponents(1)
    deformaciones_array.SetName("Deformaciones")

    # Agregar las deformaciones al vtkIntArray
    for deformacion in deformaciones:
        deformaciones_array.InsertNextValue(deformacion)

    # Añadir el array de deformaciones al vtkUnstructuredGrid
    grid.GetCellData().AddArray(deformaciones_array)

    # Crear un objeto vtkCellArray para almacenar los tetraedros
    cell_array = vtk.vtkCellArray()

    # Agregar los tetraedros al vtkCellArray
    for tetraedro in tetraedros:
        cell = vtk.vtkTetra()
        for i, nodo_id in enumerate(tetraedro):
            cell.GetPointIds().SetId(i, nodo_id)
        cell_array.InsertNextCell(cell)

    # Asignar los tetraedros al vtkUnstructuredGrid
    grid.SetCells(vtk.VTK_TETRA, cell_array)

    # Crear un escritor VTK y guardar el archivo
    writer = vtk.vtkXMLUnstructuredGridWriter()
    writer.SetFileName(nombre_archivo)
    writer.SetInputData(grid)
    writer.Write()



# Ejemplo de uso
presiones = np.array([1.5, 2.0, 1.8, 2.2])
desplazamientos = np.array([0.1, 0.2, 0.15, 0.18])
export_to_paraview('results.txt', presiones, desplazamientos)

desplazamientos = np.array([
    [0.1, 0.2, 0.3],
    [0.2, 0.3, 0.4],
    [0.3, 0.4, 0.5],
    [0.4, 0.5, 0.6]
])

strain_tensor = calcular_tensor_deformaciones(desplazamientos)
print("Tensor de deformaciones:")
print(strain_tensor)

geometria_ejemplo = {
    'nodos': np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]]),
    'tetraedros': np.array([[0, 1, 2, 3]])
}

nodos_mallado, tetraedros_mallado = generar_mallado_tetraedros(geometria_ejemplo)
print("Nodos del mallado:")
print(nodos_mallado)
print("Tetraedros del mallado:")
print(tetraedros_mallado)

xi = 0.2
eta = 0.3
zeta = 0.4
N = funcion_de_forma_tetraedro(xi, eta, zeta)
print("Valores de la función de forma:", N)

nodos = np.array([
    [0, 0, 0],
    [1, 0, 0],
    [0, 1, 0],
    [0, 0, 1]
])
tetraedro = [0, 1, 2, 3]
propiedades = {'E': 1e6, 'nu': 0.3}
matriz_rigidez_local = calcular_matriz_rigidez_local(tetraedro, nodos, propiedades)
print("Matriz de rigidez local:")
print(matriz_rigidez_local)  

matriz_rigidez_global = ensamblar_matriz_rigidez_global(nodos, tetraedros, propiedades)


desplazamientos = resolver_sistema_ecuaciones(matriz_rigidez_global, fuerzas)
print("Desplazamientos resultantes:")
print(desplazamientos)


tensiones, deformaciones = calcular_tensiones_deformaciones(matriz_rigidez_global, desplazamientos, nodos, tetraedros, propiedades)
print("Tensiones en los elementos:")
print(tensiones)
print("Deformaciones en los elementos:")
print(deformaciones)

visualizar_resultados(tensiones, deformaciones, nodos, tetraedros, "resultados.vtk")
