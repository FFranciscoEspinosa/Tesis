import pandas as pd
import numpy as np
from scipy.spatial.distance import hamming
import gudhi as gd
import seaborn as sns
from statistics import mode
import multiprocessing
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage
import gudhi as gd
from scipy.spatial.distance import hamming
import plotly.graph_objs as go
import networkx as nx
import plotly.graph_objects as go
import plotly.io as pio
from networkx.utils import not_implemented_for, pairwise
from concurrent.futures import ThreadPoolExecutor, as_completed

import module_functions as fun



#Obtener dataframes
################################################################################################################################################################
def get_info_by_pathway_and_all_genomes_directions(info):
    data_frames={}
    pathways=set()
    directions={}
    genome_directions={}
    for i in range(info.shape[0]):
        #Revisar si se debe separar esto
        #pathway=info.loc[i,0].split('|')[0].strip('\ufeff')
        pathway=info.loc[i,0].split('|')[0]
        if pathway not in directions.keys():
            directions[pathway]=[i]
        else:
            directions[pathway].append(i)
        #Se cambia la primer columna, para tener los genes
        gen_with_number=info.loc[i,0].split('|')[2]
        gen=gen_with_number[0:gen_with_number.rfind('_')]
        info.loc[i,0]=gen
        #Se obtienen los genomas
        genome=info.loc[i,1].split('|')[2]
        if genome not in genome_directions.keys():
            genome_directions[genome]=[i]
        else:
            genome_directions[genome].append(i)

    
    for pathway in directions.keys():
        data_frames[pathway]=info.loc[directions[pathway],:]
    return data_frames,genome_directions
################################################################################################################################################################
def get_df_by_pathway(info_by_pathway_and_genomes_directions,names):
    info_by_pathway=info_by_pathway_and_genomes_directions[0]
    genomes_directions=info_by_pathway_and_genomes_directions[1]

    
    names=names.set_index(1)

    
    df_by_pathway={}
    df_by_pathway_drop_duplicates={}
    representative_genomes_by_pathway={}


    for path in info_by_pathway.keys():
        info_path=info_by_pathway[path]
        genes=set(info_by_pathway[path][0].values)
        df=pd.DataFrame(index=genomes_directions.keys(),columns=list(genes))
        df.index.name=path
        
        for gen in genes:
            info_path_and_gen=info_path[info_path[0]==gen]
            for genom in genomes_directions.keys():
                direction=set(genomes_directions[genom]) & set(info_path_and_gen.index) 
                genomes_per_gen=info_path_and_gen.loc[list(direction)]
                uniq_genomes_per_gen=genomes_per_gen[1].drop_duplicates()
                df.loc[genom,gen]=uniq_genomes_per_gen.shape[0]

        df.sort_index(inplace=True)

##################
        df[path]=names[2]
        #df['index']=df.index
        df=df.set_index(path)
##################
        
        df_by_pathway[path]=df
        df_by_pathway_drop_duplicates[path]=df.drop_duplicates()
        
        df_sort_duplic=df.sort_values(by=list(df.columns)).duplicated()
        df_sort_duplic=df_sort_duplic.reset_index()
        duplicated=[]
        d=[df_sort_duplic.loc[0,path]]
        
        for i in range(1,df_sort_duplic.shape[0]):
            bool_i=df_sort_duplic.loc[i,0]
            genome_i=df_sort_duplic.loc[i,path]
            if bool_i==True:
                d.append(genome_i)
            else:
                duplicated.append(d)
                d=[genome_i]
            if i==df_sort_duplic.shape[0]-1:
                duplicated.append(d)
        representative_genomes_by_pathway[path]=duplicated

    
    return df_by_pathway, df_by_pathway_drop_duplicates, representative_genomes_by_pathway
################################################################################################################################################################

    

    

#obtener complejos simpliciales
################################################################################################################################################################
def calculate_hamming_matrix(population):
    # Number of genomes
    num_genomes = population.shape[0]
    # Create an empty matrix for Hamming distances
    hamming_matrix = np.zeros((num_genomes, num_genomes), dtype=float)
   # Calculate the Hamming distance between each pair of genomes
    size=num_genomes*(num_genomes-1)//2
    epsilon=0.0001
    increments=(np.random.RandomState(12).rand(size))*epsilon
    for i in range(num_genomes):
        for j in range(i+1, num_genomes):  # j=i+1 to avoid calculating the same distance twice
            # The Hamming distance is multiplied by the number of genes to convert it into an absolute distance
            #distance = hamming(population[i], population[j]) * len(population[i]) +increments[-1]
            distance = hamming(population[i], population[j]) * len(population[i])
            hamming_matrix[i, j] = distance
            hamming_matrix[j, i] = distance  # The matrix is symmetric
            increments=increments[:-1]
    return hamming_matrix
################################################################################################################################################################
def create_complex(distance_matrix2):
    # Create the Rips simplicial complex from the distance matrix
    rips_complex = gd.RipsComplex(distance_matrix=distance_matrix2)
    # Create the simplex tree from the Rips complex with a maximum dimension of 3
    simplex_tree = rips_complex.create_simplex_tree(max_dimension=2)
    # Compute the persistence of the simplicial complex
    persistence = simplex_tree.persistence()
    # Return the persistence diagram or barcode
    return persistence, simplex_tree
################################################################################################################################################################
def get_complex_by_pathways(path_dict):
    complex_path_dict={}
    for path in path_dict.keys():
        complex_path_dict[path]=create_complex(calculate_hamming_matrix(path_dict[path].values))
    return complex_path_dict
################################################################################################################################################################
def plot_all_bar_code_pathways(complex_pathways):
    for i in complex_pathways.keys():
        fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(12, 5))
        
        gd.plot_persistence_diagram(persistence=complex_pathways[i][0], axes=axes[0])
        gd.plot_persistence_barcode(persistence=complex_pathways[i][0], axes=axes[1])
        
        # Agregar un título a la figura
        fig.suptitle(f'Pathway: {i}', fontsize=16)

        # Ajustar el espacio entre las figuras
        plt.subplots_adjust(wspace=0.4)  # Cambia el valor según sea necesario
        
        plt.show()
################################################################################################################################################################










#obtener hoyos
################################################################################################################################################################
def minimum_cycle_basis(G, weight=None, total=None):
    """Returns a minimum weight cycle basis for G

    Minimum weight means a cycle basis for which the total weight
    (length for unweighted graphs) of all the cycles is minimum.

    Parameters
    ----------
    G : NetworkX Graph
    weight: string
        name of the edge attribute to use for edge weights

    Returns
    -------
    A list of cycle lists.  Each cycle list is a list of nodes
    which forms a cycle (loop) in G. Note that the nodes are not
    necessarily returned in a order by which they appear in the cycle

    Examples
    --------
    >>> G = nx.Graph()
    >>> nx.add_cycle(G, [0, 1, 2, 3])
    >>> nx.add_cycle(G, [0, 3, 4, 5])
    >>> nx.minimum_cycle_basis(G)
    [[5, 4, 3, 0], [3, 2, 1, 0]]

    References:
        [1] Kavitha, Telikepalli, et al. "An O(m^2n) Algorithm for
        Minimum Cycle Basis of Graphs."

        http://link.springer.com/article/10.1007/s00453-007-9064-z
        [2] de Pina, J. 1995. Applications of shortest path methods.
        Ph.D. thesis, University of Amsterdam, Netherlands


    See Also
    --------
    simple_cycles, cycle_basis
    """
    # We first split the graph in connected subgraphs
    return sum(
        (_min_cycle_basis(G.subgraph(c), weight,total) for c in nx.connected_components(G)),
        [],
    )

def _min_cycle_basis(G, weight,total):
    cb = []
    cont=0
    # We  extract the edges not in a spanning tree. We do not really need a
    # *minimum* spanning tree. That is why we call the next function with
    # weight=None. Depending on implementation, it may be faster as well
    tree_edges = list(nx.minimum_spanning_edges(G, weight=None, data=False))
    chords = G.edges - tree_edges - {(v, u) for u, v in tree_edges}

    # We maintain a set of vectors orthogonal to sofar found cycles
    set_orth = [{edge} for edge in chords]
    while set_orth:
        if cont==total:
            break
        base = set_orth.pop()
        # kth cycle is "parallel" to kth vector in set_orth
        cycle_edges = _min_cycle(G, base, weight)
        cycle_found=[v for u, v in cycle_edges]
        #print(cycle_found)

        if len(cycle_found)>3:
            cb.append(cycle_found)
            cont+=1

        # now update set_orth so that k+1,k+2... th elements are
        # orthogonal to the newly found cycle, as per [p. 336, 1]
        set_orth = [
            (
                {e for e in orth if e not in base if e[::-1] not in base}
                | {e for e in base if e not in orth if e[::-1] not in orth}
            )
            if sum((e in orth or e[::-1] in orth) for e in cycle_edges) % 2
            else orth
            for orth in set_orth
        ]
    return cb


def _min_cycle(G, orth, weight):
    """
    Computes the minimum weight cycle in G,
    orthogonal to the vector orth as per [p. 338, 1]
    Use (u, 1) to indicate the lifted copy of u (denoted u' in paper).
    """
    Gi = nx.Graph()

    # Add 2 copies of each edge in G to Gi.
    # If edge is in orth, add cross edge; otherwise in-plane edge
    for u, v, wt in G.edges(data=weight, default=1):
        if (u, v) in orth or (v, u) in orth:
            Gi.add_edges_from([(u, (v, 1)), ((u, 1), v)], Gi_weight=wt)
        else:
            Gi.add_edges_from([(u, v), ((u, 1), (v,
 1))], Gi_weight=wt)

    # find the shortest length in Gi between n and (n, 1) for each n
    # Note: Use "Gi_weight" for name of weight attribute
    spl = nx.shortest_path_length
    lift = {n: spl(Gi, source=n, target=(n, 1), weight="Gi_weight") for n in G}

    # Now compute that short path in Gi, which translates to a cycle in G
    start = min(lift, key=lift.get)
    end = (start, 1)
    min_path_i = nx.shortest_path(Gi, source=start, target=end, weight="Gi_weight")

    # Now we obtain the actual path, re-map nodes in Gi to those in G
    min_path = [n if n in G else n[0] for n in min_path_i]

    # Now remove the edges that occur two times
    # two passes: flag which edges get kept, then build it
    edgelist = list(pairwise(min_path))
    edgeset = set()
    for e in edgelist:
        if e in edgeset:
            edgeset.remove(e)
        elif e[::-1] in edgeset:
            edgeset.remove(e[::-1])
        else:
            edgeset.add(e)

    min_edgelist = []
    for e in edgelist:
        if e in edgeset:
            min_edgelist.append(e)

            edgeset.remove(e)
        elif e[::-1] in edgeset:
            min_edgelist.append(e[::-1])
            edgeset.remove(e[::-1])

    return min_edgelist
################################################################################################################################################################    
def find_all_cycles(complex_path):
    persistence,simplex_tree=complex_pathways[complex_path]
    born_and_number = born_filtraton_value_holes(persistence)
    G=nx.Graph()
    born=born_and_number.keys()
    ciclos_dep=set()
    filtration=0
    for simplex, filt in simplex_tree.get_filtration():
        #if len(ciclos_dep)==num_holes:
         #   break
        
        if filtration!=filt and filtration in born:
            number=born_and_number[filtration]
            print('se buscan ciclos en el tiempo', filt,'para',complex_path)
            ciclos=minimum_cycle_basis(G,total=number)
            for ciclo in ciclos:
                if len(ciclo)>3:
                    print('Se encontró el ciclo',ciclo,'en el tiempo', filtration,'para',complex_path)
                    ciclos_dep.add(tuple(ciclo))
                    #Se llena el hoyo
                    for i in ciclo:
                        for j in ciclo:
                            G.add_edge(i,j)
                            
            
        filtration=filt
        
        if len(simplex)==2:
            G.add_edge(simplex[0], simplex[1])


    #direc='pseudomonas_'+complex_path+'.txt'
    #f = open(direc, "a")
    #f.write(str(ciclos_dep))
    #f.close()
        

    return [complex_path,ciclos_dep]
################################################################################################################################################################    
def born_filtraton_value_holes(persistence):
    born=[]
    for bar in persistence:
        if bar[0]==1:
            born.append(bar[1][0])
            
    born_and_number=set([(x,born.count(x)) for x  in born])

    return dict(born_and_number)
################################################################################################################################################################
def get_holes_by_pathways(complex_by_path,df_by_pathway_drop_duplicate):
    llaves=complex_by_path.keys()
    with multiprocessing.Pool() as pool:
        # Mapear los ítems a la función process_item en paralelo
        results = pool.map(find_all_cycles,llaves)
    holes_by_pathway=dict(results)

    holes_by_path=change_vertex_to_name_per_path(df_by_pathway_drop_duplicate,holes_by_pathway)
    return dict(results)
################################################################################################################################################################
def change_vertex_to_name_per_path(df_by_pathway_drop_duplicate,holes_by_pathway):
    new_holes_by_path={}
    for key in holes_by_pathway.keys():
        num_to_name=df_by_pathway_drop_duplicate[key].reset_index()[key]
        new_holes=[]
        for hole in holes_by_pathway[key]:
            hole_with_name=[]
            for vertex in hole:
                hole_with_name.append(num_to_name[vertex])
            new_holes.append(hole_with_name)
        new_holes_by_path[key]=new_holes
    return new_holes_by_path
################################################################################################################################################################
def get_true_names(df,holes_by_path):
    dict_names=df
    holes_true_name={}
    for path in holes_by_path.keys():
        cycles=[]
        for cycle in holes_by_path[path]:
            cyc=[]
            for vertix in cycle:
                cyc.append(dict_names[vertix])
            cycles.append(cyc)
        holes_true_name[path]=cycles
    return holes_true_name




#Presentar información
################################################################################################################################################################
def get_resum_df(df_by_pathway,complex_by_pathway,holes_by_pathway):
    columns=['num_columnas','num_hoyos','ciclos']
    resum_df=pd.DataFrame(columns=columns)
    for path in df_by_pathway.keys():
        resum_df.loc[path,'num_columnas']=df_by_pathway[path].shape[1]
        resum_df.loc[path,'num_hoyos']=len(holes_by_pathway[path])
        resum_df.loc[path,'ciclos']=holes_by_pathway[path]
        duracion=set()
        for hole in complex_by_pathway[path][0]:
            if hole[0]==1:
                duracion.add(hole[1][1]-hole[1][0])
        resum_df.loc[path,'duración hoyos']=str(list(duracion))
    return resum_df