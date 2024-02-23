# E2. Actividad Integradora 2 

# Arturo Ramos Viedas       A01636133
# Paola Félix Torres        A00227869
# Karen Navarro Arroyo 	    A01641532
# 14/NOV/2023
 
from scipy.spatial import Voronoi, voronoi_plot_2d
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import sys
import os

def read_input():
    lines = sys.stdin.read().splitlines()
    return [line.split() for line in lines]

# Función para leer la matriz desde un archivo
def leer_matriz_desde_archivo(archivo):
    matriz = []
    with open(archivo, 'r') as file:
        for linea in file:
            fila = [int(x) for x in linea.split()]
            matriz.append(fila)
    return matriz

# Función para leer las coordenadas desde un archivo
def leer_coords_desde_archivo(archivo):
    strCoords = ""
    with open(archivo, 'r') as file:
        for linea in file:
            strCoords = strCoords + linea
    
    lines = strCoords.strip().split('\n')
    coords = []

    for l in lines:
        x, y = map(int, l.strip('()').split(','))
        coords.append([x, y])

    return coords

# [1] -------------------------------------------------------------
# Función para encontrar la forma óptima de cablear las colonias
def cablear_colonias(matriz_adyacencias):
    from collections import defaultdict
    from heapq import heapify, heappop, heappush

    # Función para encontrar la raíz de un nodo
    def encontrar_raiz(i, padre):
        while padre[i] != i:
            i = padre[i]
        return i

    n = len(matriz_adyacencias)
    arcos_formato = []

    # Estructura de datos para almacenar aristas y sus pesos
    aristas = []
    for i in range(n):
        for j in range(i + 1, n):
            aristas.append((matriz_adyacencias[i][j], i, j))

    # Ordenar las aristas por peso
    aristas.sort()

    padre = [i for i in range(n)]  

    for peso, u, v in aristas:
        raiz_u = encontrar_raiz(u, padre)
        raiz_v = encontrar_raiz(v, padre)

        if raiz_u != raiz_v:
            arcos_formato.append((chr(65 + u), chr(65 + v)))
            padre[raiz_u] = raiz_v

    return arcos_formato

# [2] -------------------------------------------------------------
def Dijkstra(matrix, s, lst): # Visita solo nodos de lst, O(N(M+N))
    N = len(matrix)
    dist = [float('inf')]*N
    visited = [False]*N
    
    dist[s] = 0
    for i in range(N):
        u = -1
        for j in lst:
            if not visited[j] and (u == -1 or dist[j] < dist[u]):
                u = j

        if dist[u] == float('inf'):
            break

        visited[u] = True
        for v in range(N):
            # Usar solo nodos en lst. No consideramos matrix[u][v] != float('inf') pues asumimos que la matriz representa un grafo completo
            if v in lst:  
                temp = dist[u] + matrix[u][v]
                if temp < dist[v]:
                    dist[v] = temp

    res = [[dist[i], i] for i in lst if dist[i] != 0]
    return res  # El resultado se devuelve como [[dist, nodo], ...]

# Estima el camino más corto para recorrer todas las colonias (Devuelve como array)
def salesMan(matrix, N): # Utilizar matriz de adyacencia
    camino = [0]
    nodosFaltantes = []
    distEstimada = 0
    
    for i in range(N):
        nodosFaltantes.append(i)
        
    i = 0
    for k in range(N):
        if k == N-1:
            nodosFaltantes = [0, i] # Regresar al primer nodo

        dists = (Dijkstra(matrix, i, nodosFaltantes)) # Hacer todas las distancias para cada nodo
        nodosFaltantes.remove(i)    # Quitarlo de nodosFaltantes
        dists.sort(key=lambda x: x[0], reverse=True) # Mover la distancia más corta al final del array
        a = dists.pop()             # a = [distanciaMin, nodo]
        distEstimada += a[0]        # Agregar la distancia min             
        i = a[1]                    # Repetir el ciclo desde el nodo con la distancia mas corta
        camino.append(i)

    return [camino, distEstimada]
    
# [3] -------------------------------------------------------------
def bfs(C, F, s, t):
    queue = [s]
    paths = {s:[]}
    if s == t:
        return paths[s]
    while queue: 
        u = queue.pop(0)
        for v in range(len(C)):
                if(C[u][v]-F[u][v]>0) and v not in paths:
                    paths[v] = paths[u]+[(u,v)]
                    if v == t:
                        return paths[v]
                    queue.append(v)
    return None

# Edmonds-Karp Algorithm
def max_flow(C, s, t):
    n = len(C) 
    F = [[0] * n for i in range(n)]
    path = bfs(C, F, s, t)
    while path != None:
        flow = min(C[u][v] - F[u][v] for u, v in path)
        for u,v in path:
            F[u][v] += flow
            F[v][u] -= flow
        path = bfs(C, F, s, t)
    return sum(F[s][i] for i in range(n))


# [4] -------------------------------------------------------------
def plotLimit(lst): # Encontrar los límites para el gráfico
    n = len(lst)-1

    lst.sort(key=lambda x: x[0]) # Encontar extremos en x
    minX = lst[0][0]
    MaxX = lst[n][0]

    lst.sort(key=lambda x: x[1]) # Encontar extremos en x
    minY = lst[0][1]
    MaxY = lst[n][1]

    totX = (MaxX-minX)*0.1
    totY = (MaxY-minY)*0.1

    return [minX, MaxX, minY, MaxY, totX, totY]

def printPolygon(pol): 
    res = []
    for i in pol:
        res.append(f"({round(i[0], ndigits=2)}, {round(i[1], ndigits=2)})")
    return res

def voronoi(coordsList):
    points = np.array(coordsList)
    vor = Voronoi(points)
    
    fig, ax = plt.subplots()
    voronoi_plot_2d(vor, ax=ax)
    
    for i, point in enumerate(points):
        ax.text(point[0], point[1], f"({point[0]},\n {point[1]})", color='#120a41', fontsize=7, ha='left', va='top')
    
    for i, vertex in enumerate(vor.vertices):
        ax.text(vertex[0], vertex[1], f"({round(vertex[0])}, {round(vertex[1])})", color='#d53c0d', fontsize=9, ha='left', va='top')

    lim = plotLimit(coordsList + [list(vertex) for vertex in vor.vertices])
    ax.set_xlim((lim[0]-lim[4], lim[1]+lim[4]*4))
    ax.set_ylim((lim[2]-lim[5], lim[3]+lim[5]*2)) 
    
    ax.set_aspect("equal")
    ax.set_title("Centrales y áreas (Diagrama Voronoi)")
    legend_labels = ['Centrales', 'Vértices', 'Lim. de area  (Cerrada)', 'Lim. de area']
    ax.legend(legend_labels, loc='upper right')
    plt.grid()
    
    polygons = []
    for region in vor.regions:
        if len(region) > 0:
            polygon = [vor.vertices[i] for i in region if i != -1]
            polygons.append(polygon)

    print(f"\nVértices de los polígonos de cada región: ")
    for i, polygon in enumerate(polygons):
        print(f"\n  Región {i+1}:") 
        for j in polygon:
            print(f"    ({round(j[0], ndigits=1)}, {round(j[1], ndigits=1)})")
    
    plt.show()
    return polygons
    
# [main] ----------------------------------------------------------
archs = ("matriz_adyacencias.txt", "matriz_capacidades.txt", "coordenadas_centrales.txt")

global Archivos
global Entradas
global N
Archivos = [""]*3
Entradas = [""]*3

def selectOption():
    global Archivos
    global Entradas
    global N

    opc = int(input(f"Seleccione un método de entrada:\n  1) Archivos de texto {archs}\n  2) Consola\n  3) Un solo archivo de texto\n    Opción: "))
    # Usando lor tres archivos txt dados
    if opc == 1:
        for i in range(3):
            Archivos[i] = os.path.abspath(archs[i])
            
        Entradas[0] = leer_matriz_desde_archivo(Archivos[0])
        Entradas[1] = leer_matriz_desde_archivo(Archivos[1])
        Entradas[2] = leer_coords_desde_archivo(Archivos[2])
        
        if any(Entradas) == []:
            input("ERROR: Archivo vacío.\nVerifique que todos los archivos contengan información. Enter para continuar")
            selectOption()
            
        N = len(Entradas[0])
             
    # Por consola (Ingresar cada matriz por separado)
    elif opc == 2:
        N = int(input("\nN: "))
        
        t = ["adyacencia", "capacidades"]
        for k in range(2):
            print(f"\nMatriz de {t[k]}: ")
            matriz = []
            for i in range(N): 
                line = sys.stdin.readline().strip()
                fila = [int(x) for x in line.split()]
                matriz.append(fila)
            Entradas[k] = matriz
            
        print("\nCoordenadas de centrales: ")
        coords = []
        for i in range(N): 
            line = sys.stdin.readline().strip()
            x, y = map(int, line.strip('()').split(','))
            coords.append([x, y])
        Entradas[2] = coords
    
    # Un solo archivo txt dando el nombre
    elif opc == 3:
        arch = os.path.abspath(input("Nombre del archivo: "))
        with open(arch, 'r') as file:
            N = int(file.readline().strip())
            file.readline()
            Entradas[0] = [list(map(int, file.readline().split())) for i in range(N)]
            file.readline()
            Entradas[1] = [list(map(int, file.readline().split())) for i in range(N)]
            file.readline()
            Entradas[2] = [list(map(int, line.strip()[1:-1].split(','))) for line in file.readlines()]

selectOption()

# [1] Forma de cablear las colonias con fibra
arcos_optimos = cablear_colonias(Entradas[0]) 
print("\nForma de cablear las colonias con fibra:")
for arco in arcos_optimos:
    print(arco)
    
# [2] ruta a seguir por el personal que reparte correspondencia, considerando inicio y fin en al misma colonia.
camino = salesMan(Entradas[0], N)
print(f"\nCamino (colonias en orden): {camino[0]}\nDistancia total del recorrido: {camino[1]}km")
 
# [3] valor de flujo máximo de información del nodo inicial al nodo final
max_flow_val = max_flow(Entradas[1], 0, N-1) # (capacity_matrix, source, sink)
print(f"\nFlujo máximo: {max_flow_val}")

# [4]
polygons = voronoi(Entradas[2]) # genCoords(strCoords)
