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
import plotly
import plotly.express as px
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
def get_df_by_pathway(info_by_pathway_and_genomes_directions):
    info_by_pathway=info_by_pathway_and_genomes_directions[0]
    genomes_directions=info_by_pathway_and_genomes_directions[1]


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

        df_by_pathway[path]=df
        df_by_pathway_drop_duplicates[path]=df.drop_duplicates()
        
        representative_genomes_by_pathway[path]=get_respresentants_drop_duplicated_df(df)

    return df_by_pathway, df_by_pathway_drop_duplicates, representative_genomes_by_pathway
################################################################################################################################################################
def get_respresentants_drop_duplicated_df(df):
    df.sort_index(inplace=True)
    df_drop_duplicates=df.drop_duplicates()
    representants={}
    name=df.index.name
    df_sort_duplic=df.sort_values(by=list(df.columns)).duplicated()
    df_sort_duplic=df_sort_duplic.reset_index()
    representants[df_sort_duplic.loc[0,name]]=[df_sort_duplic.loc[0,name]]
    actual_key=df_sort_duplic.loc[0,name]
        
    for i in range(1,df_sort_duplic.shape[0]):
        bool_i=df_sort_duplic.loc[i,0]
        indice_i=df_sort_duplic.loc[i,name]
        if bool_i==True:
            representants[actual_key].append(indice_i)
                #actual_key.append(indice_i)
        else:
            representants[df_sort_duplic.loc[i,name]]=[df_sort_duplic.loc[i,name]]
            actual_key=indice_i
    return representants
################################################################################################################################################################
def traduced_df(df,names):
    diccionario=dict(names[[1,2]].values)
    new_df=df.copy()
    df_name=df.index.name
    new_df[df_name]=np.nan
    for indice in new_df.index:
        new_df.loc[indice,df_name]=diccionario[indice]
    new_df=new_df.set_index(df_name)
    #new_df=new_df.drop(column=[df_name])
    return new_df
################################################################################################################################################################
def traduced_dict_data_frames(dict_df,names):
    traduced_dict={}
    for key,df in dict_df.items():
        traduced_dict[key]=traduced_df(df,names)
    return traduced_dict
################################################################################################################################################################
def get_df_by_genus_pathway(df_by_pathway,names):
    indice=list(df_by_pathway.values())[0].index
    diccionario=dict(names[[1,2]].values)
    genus=pd.Series(diccionario)[indice]
    genus=genus.apply(lambda x: x.split(' ')[0])
    df_genus_pathway={}
    df_by_genus_pathway_drop_duplicates={}
    representative_genomes_by_genus_pathway={}
    
    for path,df in df_by_pathway.items():
        for name_genus in genus.unique():
            name_genus_path=name_genus+'_'+path

            df_genus_path=df[genus==name_genus].copy()
            df_genus_path.index.name=name_genus_path
            
            df_genus_pathway[name_genus_path]=df_genus_path
            df_by_genus_pathway_drop_duplicates[name_genus_path]=df_genus_path.drop_duplicates()
            representative_genomes_by_genus_pathway[name_genus_path]=get_respresentants_drop_duplicated_df(df_genus_path)



    
    return df_genus_pathway,df_by_genus_pathway_drop_duplicates,representative_genomes_by_genus_pathway
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
def find_all_cycles(complex_dict):
    complex_path=complex_dict[0]
    persistence,simplex_tree=complex_dict[1]
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
    non_empty={}
    empty={}
    
    for path,comp in complex_by_path.items():
        if comp[0][0][0]==1:
            non_empty[path]=comp
        else:
            empty[path]=[]

    result_list=[]
    with ThreadPoolExecutor() as executor:
        futures=[]
        for item in list(non_empty.items()):
            future = executor.submit(find_all_cycles, item)
            futures.append(future)
        for f in futures:
            result_list.append(f.result())
    holes_by_pathway=dict(result_list)
    holes_by_pathway.update(empty)

    holes_by_path=change_vertex_to_code_per_path(holes_by_pathway,df_by_pathway_drop_duplicate)


    
    return holes_by_path
################################################################################################################################################################
def change_vertex_to_code_per_path(holes_by_pathway,df_by_pathway_drop_duplicate):
    
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
def change_vertex_to_representants(holes_by_pathway,representative_genomes_by_pathway):
    new_holes_by_path={}
    for key in holes_by_pathway.keys():
        representant_to_all_genomes=representative_genomes_by_pathway[key]
        #num_to_name=df_by_pathway_drop_duplicate[key].reset_index()[key]
        new_holes=[]
        for hole in holes_by_pathway[key]:
            hole_with_name=[]
            for vertex in hole:
                vertex_to_all_vertex_names=representant_to_all_genomes[vertex]
                if len(vertex_to_all_vertex_names)==1:
                    hole_with_name.append(vertex)
                else:
                    equivalent_vertex=[]
                    for i in vertex_to_all_vertex_names:
                        equivalent_vertex.append(i)
                    hole_with_name.append(equivalent_vertex)
                

            new_holes.append(hole_with_name)
        new_holes_by_path[key]=new_holes
    return new_holes_by_path
################################################################################################################################################################
def change_vertex_to_names_per_path(holes_by_pathway,names):
    code_to_names=dict(names[[1,2]].values)
    new_holes_by_path={}
    for path,holes in holes_by_pathway.items():
        new_holes=[]
        for hole in holes:
            hole_with_name=[]
            for vertex in hole:
                if type(vertex)==str:
                    hole_with_name.append(code_to_names[vertex])
                else:
                    equivalent_vertex=[]
                    for i in vertex:
                        equivalent_vertex.append(code_to_names[i])
                    hole_with_name.append(equivalent_vertex)

            new_holes.append(hole_with_name)
        new_holes_by_path[path]=new_holes
    return new_holes_by_path
################################################################################################################################################################











#Presentar información
################################################################################################################################################################
def get_resum_df(df_by_pathway,representative_genomes_by_pathway,names,complex_by_pathway,holes_by_pathway):
    holes_by_pathway=change_vertex_to_names_per_path(change_vertex_to_representants(holes_by_pathway,representative_genomes_by_pathway),names)
    holes_genus_pathway=holes_per_genus(holes_by_pathway)
    columns=['num_columnas','num_hoyos','ciclos','ciclos por genero']
    resum_df=pd.DataFrame(columns=columns)
    for path in df_by_pathway.keys():
        resum_df.loc[path,'num_columnas']=df_by_pathway[path].shape[1]
        resum_df.loc[path,'num_hoyos']=len(holes_by_pathway[path])
        resum_df.loc[path,'ciclos']=holes_by_pathway[path]
        resum_df.loc[path,'ciclos por genero']=holes_genus_pathway[path]
        duracion=set()
        for hole in complex_by_pathway[path][0]:
            if hole[0]==1:
                duracion.add(hole[1][1]-hole[1][0])
        resum_df.loc[path,'duración hoyos']=str(list(duracion))
    return resum_df
################################################################################################################################################################
def holes_per_genus(holes_by_pathway):
    holes_genus_by_pathway={}
    for path in holes_by_pathway.keys():
        holes_genus=[]
        for hole in holes_by_pathway[path]:
            genus_hole=[]
            for genome in hole:
                if type(genome)==str:
                    #print(genome)
                    genus_hole.append(genome.split(' ')[0])
                else:
                    equivalent_genomes_genus=[]
                    for i in genome:
                        equivalent_genomes_genus.append(i.split(' ')[0])
                    genus_hole.append(equivalent_genomes_genus)
                        
            holes_genus.append(genus_hole)
        holes_genus_by_pathway[path]=list(holes_genus)
    return holes_genus_by_pathway
################################################################################################################################################################
def get_df_by_pathway_make_hole(df_by_pathway_drop_duplicate,holes_by_pathway):
    df_by_pathway_make_hole={}
    for path,df in df_by_pathway_drop_duplicate.items():
        holes_per_path=[]
        for hole in holes_by_pathway[path]:
            query=[]

            for vertex in hole:
                query.append(vertex)
            
            holes_per_path.append(df.loc[query,])
        df_by_pathway_make_hole[path]=holes_per_path
    return df_by_pathway_make_hole
################################################################################################################################################################









#Graficar
################################################################################################################################################################
def visualize_simplicial_complex(simplex_tree, filtration_value, path, family,vertex_names=None, save_filename=None, plot_size=1, dpi=600, pos=None):
    G = nx.Graph()
    triangles = []  # List to store triangles (3-nodes simplices)
    
    for simplex, filt in simplex_tree.get_filtration():
        if filt <= filtration_value:
            if len(simplex) == 2:
                G.add_edge(simplex[0], simplex[1])
            elif len(simplex) == 1:
                G.add_node(simplex[0])
            elif len(simplex) == 3:
                triangles.append(simplex)
    
    # Calculate node positions if not provided
    if pos is None:
        pos = nx.spring_layout(G,dim=2)
    
    # Node trace
    x_values, y_values = zip(*[pos[node] for node in G.nodes()])
    node_labels = [vertex_names[node] if vertex_names else str(node) for node in G.nodes()]
    node_trace = go.Scatter(x=x_values, y=y_values, mode='markers+text', hoverinfo='text', marker=dict(size=14), text=node_labels, textposition='top center', textfont=dict(size=14))
        # Edge traces
    edge_traces = []
    for edge in G.edges():
        x0, y0 = pos[edge[0]]
        x1, y1 = pos[edge[1]]
        edge_trace = go.Scatter(x=[x0, x1, None], y=[y0, y1, None], mode='lines', line=dict(width=3, color='rgba(0,0,0,0.5)'))
        edge_traces.append(edge_trace)
    
    # Triangle traces
    triangle_traces = []
    for triangle in triangles:
        x0, y0 = pos[triangle[0]]
        x1, y1 = pos[triangle[1]]
        x2, y2 = pos[triangle[2]]
        triangle_trace = go.Scatter(x=[x0, x1, x2, x0, None], y=[y0, y1, y2, y0, None], fill='toself', mode='lines+markers', line=dict(width=2), fillcolor='rgba(255,0,0,0.2)')
        triangle_traces.append(triangle_trace)
    titulo='Complejo simplicial de la ruta '+path+' en '+family+' para nivel de filtración '+str(filtration_value)
    # Configure the layout of the plot
    layout = go.Layout(showlegend=False, hovermode='closest', xaxis=dict(showgrid=False, zeroline=False, tickfont=dict(size=16, family='Arial, sans-serif')), yaxis=dict(showgrid=False, zeroline=False, tickfont=dict(size=16, family='Arial, sans-serif')),title=dict(
            text=titulo ,
            x=0.5,
            xanchor='center',
            font=dict(size=24, family='Arial, sans-serif')
        ))
    fig = go.Figure(data=edge_traces + triangle_traces + [node_trace], layout=layout)
    
    # Set the figure size
    fig.update_layout(width=plot_size * dpi, height=plot_size * dpi)
    

    # Show the figure
    #fig.write_image('M.png')
    filename=path+'_'+family+'_t_'+str(filtration_value)+'.html'
    plotly.offline.plot(fig, filename=filename)
    

    return type(fig) 
################################################################################################################################################################

















def holes_pathway_strepto():
    holes={'OXALACETATE_AMINOACIDS': [[5, 59, 77, 101, 108],
  [5, 59, 62, 89, 101, 120, 166],
  [41, 55, 69, 93],
  [5, 13, 19, 21, 101, 114],
  [10, 47, 57, 67, 301],
  [20, 27, 114, 166],
  [17, 28, 115, 122],
  [125, 128, 132, 298],
  [5, 57, 62, 67, 117],
  [20, 41, 69, 77, 114],
  [43, 44, 112, 123, 131, 172],
  [5, 27, 28, 117, 133, 166],
  [21, 128, 130, 132, 304],
  [21, 27, 101, 166],
  [10, 133, 161, 165],
  [21, 64, 101, 304],
  [68, 72, 130, 132, 143, 157, 304],
  [17, 27, 101, 111, 166],
  [10, 47, 87, 88, 165],
  [21, 27, 123, 128, 131, 132],
  [5, 47, 62, 64, 67, 68, 106, 304],
  [37, 125, 130, 132],
  [33, 52, 98, 128, 131, 176, 303],
  [33, 80, 131, 161],
  [3, 31, 69, 156],
  [68, 71, 123, 132],
  [40, 60, 135, 162],
  [5, 36, 93, 127],
  [0, 11, 128, 131],
  [201, 239, 260, 289],
  [5, 28, 77, 110],
  [0, 71, 112, 131],
  [5, 13, 151, 304],
  [27, 38, 55, 131, 159],
  [10, 13, 67, 95],
  [43, 103, 104, 112, 134],
  [7, 53, 132, 304],
  [0, 112, 131, 180],
  [0, 105, 178, 306],
  [27, 131, 137, 298],
  [28, 133, 134, 135],
  [23, 26, 71, 90, 112, 301],
  [10, 57, 105, 116],
  [3, 10, 61, 158, 181, 183],
  [44, 112, 134, 135],
  [19, 48, 50, 101, 116],
  [27, 55, 159, 298],
  [1, 28, 111, 113, 155, 260],
  [5, 19, 68, 157],
  [21, 28, 115, 129],
  [18, 53, 66, 128, 132],
  [111, 260, 289, 301],
  [74, 88, 98, 122],
  [7, 49, 68, 308],
  [3, 5, 183, 184, 227],
  [3, 5, 31, 36],
  [2, 21, 137, 164],
  [40, 64, 175, 201],
  [5, 35, 54, 106, 184],
  [2, 5, 54, 164, 184],
  [45, 77, 91, 108],
  [39, 138, 175, 234],
  [1, 138, 234, 307],
  [1, 260, 273, 307],
  [88, 144, 178, 179],
  [7, 23, 24, 63, 133],
  [21, 76, 303, 304],
  [5, 102, 127, 157],
  [39, 88, 138, 289],
  [39, 138, 225, 289],
  [88, 138, 179, 233],
  [81, 115, 176, 199],
  [130, 133, 225, 289],
  [5, 45, 47, 108],
  [27, 74, 88, 159],
  [21, 45, 108, 176],
  [59, 63, 187, 239],
  [25, 75, 125, 139],
  [39, 43, 88, 138],
  [3, 5, 28, 155],
  [61, 217, 239, 276],
  [5, 182, 183, 225],
  [217, 239, 260, 292],
  [148, 215, 219, 289],
  [183, 194, 227, 249],
  [239, 252, 260, 292],
  [35, 132, 143, 237],
  [16, 66, 141, 298],
  [168, 171, 215, 289],
  [5, 117, 217, 239, 282],
  [5, 195, 227, 260],
  [117, 138, 197, 272],
  [134, 136, 168, 171],
  [63, 195, 260, 302],
  [5, 79, 182, 183],
  [59, 147, 195, 260],
  [1, 61, 224, 273],
  [1, 195, 234, 260],
  [201, 215, 235, 289],
  [63, 260, 292, 302],
  [64, 173, 201, 220, 277],
  [63, 166, 270, 302],
  [61, 166, 230, 276],
  [132, 164, 237, 269],
  [1, 166, 191, 208, 230],
  [61, 131, 193, 231, 276],
  [99, 111, 249, 275],
  [81, 84, 199, 248, 272, 274],
  [61, 204, 217, 255],
  [29, 62, 230, 267],
  [111, 131, 148, 288],
  [29, 44, 206, 253],
  [185, 232, 241, 292],
  [5, 112, 174, 195],
  [111, 221, 267, 289],
  [81, 110, 228, 288],
  [3, 64, 255, 277],
  [64, 220, 257, 260, 290],
  [111, 179, 232, 292],
  [1, 208, 233, 266],
  [112, 171, 176, 192, 242],
  [111, 260, 275, 290],
  [224, 242, 259, 276, 277],
  [81, 193, 231, 233, 266],
  [209, 232, 273, 293],
  [28, 166, 191, 193],
  [54, 149, 198, 237],
  [63, 181, 275, 280],
  [27, 210, 275, 278],
  [29, 64, 206, 214, 257],
  [120, 126, 243, 265],
  [181, 220, 259, 280],
  [107, 115, 134, 222],
  [91, 112, 192, 261],
  [84, 167, 175, 272],
  [61, 63, 179, 204],
  [5, 85, 195, 285],
  [204, 214, 281, 283],
  [112, 210, 242, 278],
  [257, 260, 261, 279],
  [12, 185, 203, 207, 218, 241, 271, 272, 286],
  [61, 204, 244, 281],
  [43, 148, 268, 279],
  [162, 211, 268, 279],
  [29, 126, 206, 265, 287],
  [64, 147, 180, 238],
  [223, 242, 259, 281, 283],
  [29, 180, 238, 253],
  [148, 204, 215, 263, 279],
  [64, 245, 255, 257],
  [64, 204, 214, 255, 257],
  [54, 126, 198, 254, 265],
  [61, 231, 244, 295],
  [27, 256, 275, 278],
  [9, 107, 115, 251],
  [5, 250, 285, 292],
  [1, 272, 286, 295],
  [12, 272, 285, 307],
  [29, 65, 86, 287],
  [16, 188, 246, 276],
  [195, 246, 276, 284],
  [64, 214, 238, 294],
  [12, 148, 205, 272]],
 'E4P_AMINO_ACIDS': [[7, 19, 65, 93, 144, 160],
  [98, 101, 138, 282],
  [7, 113, 120, 121],
  [61, 65, 67, 108, 138, 160],
  [7, 69, 120, 148],
  [7, 19, 148, 152],
  [1, 18, 140, 281],
  [19, 93, 118, 152],
  [15, 55, 81, 109],
  [5, 113, 115, 120],
  [2, 8, 64, 114, 136, 151],
  [15, 45, 48, 49, 61, 98],
  [18, 69, 82, 112, 140],
  [5, 8, 26, 58, 64, 71, 86, 136],
  [37, 44, 61, 143],
  [108, 118, 119, 160],
  [5, 13, 26, 138],
  [5, 54, 108, 138],
  [1, 123, 140, 150],
  [7, 55, 93, 113],
  [5, 39, 101, 138],
  [61, 67, 144, 163],
  [55, 81, 93, 160],
  [1, 82, 103, 287],
  [1, 40, 104, 208],
  [112, 122, 151, 161],
  [82, 95, 128, 287],
  [45, 48, 83, 92],
  [2, 27, 63, 73],
  [67, 94, 111, 122],
  [1, 16, 97, 129],
  [3, 5, 115, 116],
  [27, 102, 131, 285],
  [11, 28, 29, 162],
  [28, 29, 103, 287],
  [1, 104, 128, 281],
  [44, 61, 99, 131],
  [1, 18, 165, 210],
  [11, 29, 49, 162],
  [16, 94, 115, 129, 284],
  [7, 148, 164, 281],
  [9, 95, 97, 129, 284],
  [37, 61, 63, 131],
  [26, 111, 122, 151],
  [5, 26, 144, 147],
  [24, 50, 56, 58, 87],
  [174, 183, 202, 257],
  [65, 76, 131, 161],
  [17, 132, 159, 171],
  [16, 153, 177, 270],
  [1, 179, 183, 213],
  [1, 18, 176, 179],
  [18, 132, 133, 214],
  [126, 165, 219, 238],
  [165, 195, 239, 259],
  [126, 195, 238, 239],
  [179, 221, 245, 250],
  [126, 238, 243, 256],
  [1, 84, 120, 174],
  [195, 236, 238, 269],
  [179, 231, 255, 288],
  [200, 207, 231, 235, 246, 266, 281],
  [167, 185, 199, 236, 259, 261, 269],
  [132, 133, 167, 185, 199, 242],
  [17, 91, 124, 276],
  [56, 179, 195, 221, 245, 259],
  [165, 193, 219, 235],
  [219, 238, 246, 269],
  [1, 198, 206, 208, 213],
  [12, 240, 264, 265],
  [12, 196, 199, 232, 256, 265],
  [185, 199, 209, 247],
  [184, 209, 212, 251],
  [186, 189, 215, 217, 226, 275, 277],
  [174, 191, 225, 275],
  [173, 186, 189, 232, 265, 276, 277],
  [200, 225, 233, 260],
  [167, 185, 204, 208],
  [1, 90, 135, 249],
  [153, 179, 253, 275],
  [126, 197, 199, 256],
  [158, 203, 219, 277],
  [153, 182, 237, 267],
  [185, 232, 263, 264],
  [217, 252, 275, 277],
  [185, 209, 230, 236, 274],
  [1, 197, 227, 273],
  [200, 203, 217, 277],
  [191, 200, 231, 260],
  [153, 194, 237, 267],
  [181, 219, 224, 276],
  [202, 217, 226, 262],
  [183, 201, 223, 273],
  [79, 223, 271, 273],
  [12, 199, 223, 273],
  [1, 171, 209, 247],
  [139, 196, 260, 264],
  [195, 196, 221, 256],
  [12, 172, 223, 265],
  [196, 221, 227, 262],
  [139, 199, 247, 265],
  [139, 186, 217, 260],
  [153, 209, 226, 275],
  [1, 179, 180, 192],
  [132, 193, 196, 199],
  [183, 207, 223, 240],
  [139, 172, 217, 252],
  [139, 187, 189, 221],
  [79, 209, 229, 255]],
 'PYR_THR_AA': [[2, 41, 62, 63],
  [54, 68, 70, 136],
  [0, 39, 172, 303],
  [41, 54, 62, 98, 104, 114, 136],
  [0, 11, 100, 174],
  [111, 135, 173, 303],
  [0, 33, 39, 114, 135, 136, 221, 303],
  [3, 80, 122, 132, 151, 188],
  [8, 70, 141, 151, 163],
  [31, 56, 73, 142, 164, 297],
  [23, 25, 26, 36, 148],
  [18, 26, 229, 294],
  [3, 19, 137, 161],
  [9, 25, 36, 37, 74, 184],
  [60, 128, 145, 263],
  [18, 62, 64, 94],
  [38, 112, 155, 160, 180],
  [100, 129, 225, 235],
  [46, 80, 122, 131],
  [27, 32, 49, 160],
  [31, 47, 104, 295],
  [6, 33, 128, 180],
  [70, 98, 163, 175],
  [9, 28, 114, 290],
  [18, 89, 153, 229],
  [25, 36, 65, 100, 134],
  [1, 33, 104, 153],
  [3, 25, 31, 229, 294],
  [3, 38, 88, 304],
  [17, 20, 65, 66, 88],
  [3, 45, 70, 80, 122, 163],
  [52, 53, 142, 164],
  [20, 52, 53, 157],
  [0, 2, 102, 136],
  [31, 38, 56, 304],
  [19, 47, 56, 137],
  [8, 22, 24, 65, 151],
  [62, 86, 98, 125],
  [33, 114, 153, 290],
  [31, 33, 153, 229],
  [6, 38, 56, 180],
  [3, 80, 99, 100],
  [3, 46, 53, 80],
  [153, 206, 222, 276, 291],
  [10, 105, 145, 162, 183],
  [31, 122, 153, 214, 239, 272, 281],
  [1, 148, 182, 268],
  [199, 253, 254, 259],
  [9, 184, 198, 204],
  [104, 132, 146, 191, 197, 217],
  [9, 208, 250, 290],
  [189, 196, 230, 259],
  [189, 208, 219, 273, 288],
  [40, 66, 101, 169],
  [223, 232, 254, 288],
  [186, 213, 239, 245, 264, 279],
  [189, 197, 217, 267],
  [203, 239, 267, 286],
  [201, 213, 228, 240, 253],
  [196, 197, 217, 226],
  [45, 148, 182, 214],
  [187, 189, 243, 259, 267],
  [198, 241, 268, 270],
  [189, 219, 254, 259, 288],
  [199, 201, 210, 253, 292],
  [46, 80, 103, 177],
  [50, 68, 72, 110],
  [1, 241, 250, 268],
  [153, 198, 204, 291],
  [189, 259, 267, 286],
  [9, 208, 267, 272, 273, 290],
  [101, 105, 145, 169, 183],
  [153, 178, 290, 291],
  [239, 259, 264, 286],
  [206, 210, 212, 259, 264],
  [1, 211, 272, 276],
  [216, 226, 228, 230, 254, 259],
  [31, 132, 153, 217, 267, 272],
  [37, 198, 206, 276],
  [199, 238, 239, 252, 259, 264, 281],
  [101, 123, 157, 169, 177],
  [213, 216, 228, 240],
  [186, 228, 240, 254, 259, 264],
  [33, 183, 193, 195, 208],
  [39, 135, 249, 292],
  [57, 72, 185, 266],
  [19, 66, 72, 106],
  [39, 192, 198, 270],
  [166, 200, 251, 282],
  [207, 241, 268, 274],
  [9, 152, 264, 272, 284],
  [40, 72, 96, 266],
  [22, 195, 208, 277],
  [214, 248, 251, 282],
  [57, 101, 183, 193],
  [175, 209, 256, 289],
  [24, 123, 166, 251, 287],
  [187, 205, 217, 280],
  [12, 213, 245, 253],
  [5, 80, 103, 138],
  [9, 135, 171, 250],
  [31, 114, 217, 280],
  [1, 74, 84, 199, 271],
  [31, 170, 267, 287],
  [19, 72, 170, 266, 287],
  [85, 173, 212, 292],
  [168, 233, 242, 274],
  [45, 175, 261, 289],
  [153, 222, 242, 260],
  [53, 78, 82, 105],
  [81, 123, 138, 191],
  [39, 78, 105, 197]],
 'TCA': [[3, 11, 17, 98],
  [11, 18, 98, 112],
  [8, 24, 35, 70],
  [20, 26, 110, 115],
  [25, 33, 39, 45],
  [33, 39, 41, 50],
  [11, 24, 75, 98],
  [0, 23, 47, 63],
  [8, 11, 24, 32],
  [0, 34, 47, 112],
  [25, 42, 45, 71, 76, 77, 81],
  [25, 39, 65, 81],
  [26, 39, 41, 65],
  [33, 39, 40, 49],
  [0, 11, 32, 63],
  [44, 51, 68, 105],
  [0, 11, 84, 92],
  [20, 26, 41, 74, 94, 110, 234],
  [0, 11, 24, 47],
  [33, 40, 45, 67],
  [39, 41, 44, 94],
  [42, 45, 53, 62],
  [25, 44, 68, 81],
  [23, 63, 66, 87],
  [5, 8, 23, 75],
  [5, 24, 47, 75],
  [2, 39, 41, 72, 74, 94, 234],
  [13, 33, 38, 50],
  [2, 20, 21, 234],
  [7, 19, 20, 21],
  [33, 39, 44, 51],
  [13, 38, 42, 45],
  [23, 47, 64, 66],
  [90, 106, 170, 179],
  [0, 78, 110, 234],
  [2, 11, 44, 99],
  [90, 106, 131, 170],
  [116, 131, 155, 160, 184],
  [0, 63, 69, 234],
  [90, 106, 135, 193],
  [0, 1, 90, 92],
  [3, 21, 65, 82],
  [1, 159, 160, 184, 227],
  [61, 82, 88, 89],
  [1, 90, 106, 122, 123, 160, 227],
  [0, 5, 232, 234],
  [1, 90, 106, 129],
  [1, 90, 120, 130, 179, 227],
  [0, 2, 11, 234],
  [0, 3, 21, 234],
  [91, 124, 134, 173],
  [142, 149, 180, 188],
  [0, 90, 97, 195],
  [117, 133, 156, 161, 177, 179],
  [90, 125, 160, 203],
  [140, 141, 143, 188, 192, 226],
  [84, 90, 91, 124],
  [122, 123, 144, 201, 213],
  [156, 164, 190, 207, 214],
  [138, 199, 209, 221],
  [117, 156, 164, 176, 204, 214],
  [5, 16, 27, 125],
  [46, 133, 162, 170, 177],
  [158, 167, 193, 222],
  [144, 149, 180, 213, 219],
  [117, 142, 176, 188, 196, 204],
  [137, 138, 208, 209],
  [90, 142, 193, 194, 196, 222],
  [26, 82, 111, 233],
  [90, 118, 160, 203],
  [27, 106, 119, 168],
  [1, 5, 135, 153, 158, 167],
  [5, 73, 75, 125],
  [1, 122, 135, 163, 180, 189],
  [137, 175, 187, 209],
  [90, 123, 192, 207],
  [141, 142, 143, 156, 188, 196],
  [137, 142, 181, 199, 209],
  [142, 144, 181, 212],
  [138, 199, 207, 228],
  [137, 154, 186, 205],
  [136, 151, 175, 223],
  [17, 111, 124, 176],
  [160, 203, 225, 228],
  [124, 128, 197, 217],
  [5, 8, 53, 167],
  [136, 149, 175, 181],
  [138, 156, 199, 214],
  [0, 1, 55, 163],
  [4, 17, 100, 231],
  [174, 190, 207, 228],
  [117, 136, 151, 204],
  [151, 154, 185, 223],
  [151, 154, 215, 223],
  [148, 157, 181, 190],
  [6, 55, 147, 185],
  [5, 6, 158, 185],
  [2, 58, 121, 228],
  [1, 60, 171, 223]],
 'R5P_AMINOACIDS': [[1, 4, 7, 52],
  [1, 52, 54, 101],
  [5, 23, 54, 96],
  [57, 66, 78, 99],
  [1, 15, 54, 96],
  [1, 26, 49, 54],
  [1, 49, 52, 53],
  [5, 7, 33, 46],
  [41, 57, 99, 111],
  [5, 7, 18, 21],
  [0, 2, 7, 22],
  [0, 6, 11, 13],
  [4, 25, 36, 39],
  [2, 5, 6, 22],
  [1, 5, 7, 54],
  [0, 7, 21, 27],
  [0, 4, 7, 11],
  [5, 23, 44, 50],
  [4, 5, 7, 14],
  [1, 54, 57, 99],
  [2, 7, 10, 25],
  [54, 60, 99, 101],
  [52, 60, 94, 101],
  [5, 29, 54, 55],
  [35, 60, 78, 99],
  [5, 6, 13, 14],
  [1, 91, 114, 124],
  [66, 73, 91, 92],
  [1, 7, 21, 49],
  [1, 34, 57, 61, 87, 95, 111],
  [5, 14, 54, 101],
  [0, 1, 7, 61],
  [5, 18, 26, 54],
  [5, 14, 16, 50],
  [0, 7, 12, 25],
  [5, 14, 18, 38],
  [7, 21, 33, 51],
  [26, 54, 99, 103],
  [2, 5, 10, 40],
  [1, 52, 57, 94],
  [1, 57, 66, 91],
  [5, 14, 28, 46],
  [63, 117, 140, 146],
  [89, 102, 117, 140],
  [74, 79, 89, 109],
  [63, 74, 89, 117],
  [15, 20, 110, 157],
  [30, 57, 99, 149],
  [30, 48, 54, 99, 106],
  [32, 73, 79, 86, 98, 108, 151],
  [54, 59, 68, 76],
  [84, 87, 105, 135],
  [1, 64, 129, 138],
  [80, 139, 141, 151],
  [63, 113, 141, 146, 151],
  [30, 99, 124, 126],
  [84, 87, 116, 144],
  [71, 89, 106, 133, 140],
  [1, 2, 31, 87],
  [30, 71, 79, 102, 115, 126, 132],
  [74, 97, 117, 131],
  [71, 79, 141, 151],
  [64, 116, 129, 144],
  [48, 49, 54, 128],
  [1, 59, 87, 156],
  [77, 91, 92, 100],
  [1, 62, 124, 126, 132],
  [48, 54, 64, 90, 93, 106, 133, 138],
  [64, 85, 93, 153],
  [73, 88, 108, 142],
  [16, 105, 108, 142]],
 'ALPHAKETOGLUTARATE_AMINOACIDS': [[35, 112, 120, 276],
  [29, 64, 85, 279],
  [3, 38, 112, 124],
  [13, 25, 29, 85],
  [5, 70, 110, 125],
  [35, 79, 131, 276],
  [5, 58, 102, 109],
  [5, 102, 110, 119],
  [13, 25, 107, 157],
  [79, 86, 93, 122],
  [58, 100, 101, 109],
  [102, 103, 109, 286],
  [100, 107, 109, 286],
  [25, 97, 100, 101],
  [11, 103, 107, 286],
  [35, 46, 89, 276],
  [11, 15, 45, 157],
  [3, 5, 57, 110],
  [36, 58, 101, 113],
  [0, 13, 85, 91],
  [15, 45, 79, 122],
  [5, 57, 58, 89, 97, 101],
  [35, 57, 89, 112],
  [13, 25, 35, 38, 50, 63, 64, 85, 89, 97, 112],
  [5, 10, 26, 43, 57, 69, 115],
  [5, 10, 15, 45, 57, 69, 89],
  [35, 45, 79, 89],
  [11, 25, 45, 89, 97, 107],
  [10, 15, 43, 95],
  [120, 131, 153, 276],
  [20, 122, 151, 284],
  [7, 21, 23, 284],
  [67, 80, 156, 281],
  [17, 29, 88, 142],
  [8, 15, 45, 143],
  [7, 17, 112, 124],
  [29, 51, 123, 279],
  [10, 56, 145, 150],
  [10, 128, 146, 150],
  [10, 61, 122, 128],
  [3, 94, 110, 145],
  [56, 58, 139, 145],
  [7, 17, 40, 46],
  [1, 38, 84, 112],
  [1, 10, 20, 58],
  [110, 132, 133, 145],
  [43, 57, 59, 82],
  [1, 97, 137, 262],
  [28, 103, 168, 173],
  [17, 18, 21, 111],
  [10, 103, 162, 168],
  [7, 124, 179, 199, 256],
  [3, 38, 83, 132],
  [1, 17, 179, 255],
  [1, 10, 162, 207, 255],
  [38, 55, 83, 280],
  [1, 10, 158, 165, 255],
  [3, 10, 162, 175],
  [110, 165, 203, 260],
  [203, 220, 234, 245],
  [159, 167, 207, 223],
  [80, 168, 172, 261],
  [10, 25, 168, 205, 227],
  [165, 175, 199, 241],
  [207, 220, 245, 253, 255],
  [159, 162, 178, 225],
  [10, 84, 168, 247],
  [159, 214, 255, 257, 271],
  [2, 10, 171, 207],
  [158, 212, 255, 257, 271],
  [189, 207, 240, 255, 271],
  [150, 165, 212, 241],
  [185, 191, 204, 210],
  [13, 97, 209, 242],
  [159, 210, 223, 233, 249],
  [165, 180, 185, 204, 256],
  [167, 187, 190, 237],
  [200, 228, 238, 251],
  [148, 149, 190, 228, 238],
  [191, 204, 223, 264],
  [187, 190, 197, 238],
  [167, 207, 264, 271],
  [201, 204, 223, 229],
  [148, 190, 201, 238, 251],
  [3, 135, 182, 265],
  [165, 204, 207, 223],
  [56, 167, 175, 237],
  [0, 78, 158, 275],
  [185, 187, 192, 195],
  [56, 201, 206, 222, 237],
  [132, 148, 159, 223],
  [148, 190, 204, 223],
  [167, 175, 187, 243],
  [163, 173, 247, 269],
  [185, 188, 190, 268],
  [13, 209, 273, 274],
  [159, 210, 268, 270],
  [119, 202, 236, 250],
  [13, 119, 163, 274],
  [148, 162, 238, 253],
  [40, 72, 99, 190],
  [72, 77, 190, 238],
  [190, 193, 228, 239],
  [167, 190, 226, 251],
  [167, 188, 189, 190],
  [29, 219, 222, 231],
  [3, 148, 170, 197],
  [13, 72, 81, 259],
  [55, 193, 194, 230],
  [3, 173, 197, 263],
  [0, 4, 127, 136],
  [0, 14, 127, 258]],
 'PPP': [[69, 71, 139, 170],
  [0, 2, 11, 67],
  [18, 44, 63, 86, 90, 113],
  [0, 18, 20, 88],
  [47, 69, 76, 78, 82, 102, 112, 113],
  [9, 17, 38, 76],
  [9, 32, 33, 50],
  [97, 140, 145, 176],
  [5, 13, 31, 68],
  [0, 12, 55, 67],
  [1, 43, 48, 49],
  [2, 5, 11, 31],
  [0, 20, 44, 55],
  [9, 38, 83, 118],
  [5, 9, 13, 32],
  [0, 6, 49, 67],
  [0, 16, 47, 49],
  [9, 38, 47, 49],
  [1, 5, 13, 48],
  [1, 5, 30, 31],
  [9, 38, 53, 54],
  [0, 1, 30, 67],
  [11, 28, 33, 50],
  [11, 28, 31, 68],
  [9, 24, 25, 33],
  [4, 14, 15, 20],
  [9, 32, 43, 49],
  [9, 38, 42, 120],
  [59, 126, 145, 162],
  [1, 48, 55, 56],
  [0, 1, 15, 20],
  [7, 16, 169, 176],
  [59, 145, 169, 176],
  [17, 44, 76, 113],
  [0, 20, 26, 67],
  [0, 1, 2, 5],
  [82, 95, 102, 136],
  [71, 134, 139, 153],
  [2, 9, 11, 33],
  [4, 20, 27, 44],
  [0, 2, 9, 49],
  [0, 17, 20, 49],
  [54, 69, 113, 170],
  [0, 2, 53, 55],
  [9, 17, 22, 33],
  [14, 27, 48, 56],
  [1, 14, 15, 48],
  [39, 73, 110, 145],
  [9, 25, 53, 62],
  [38, 54, 76, 113],
  [11, 22, 26, 33],
  [0, 1, 7, 16],
  [0, 2, 118, 122],
  [6, 9, 33, 49],
  [79, 101, 105, 147],
  [79, 101, 163, 166],
  [0, 4, 138, 142],
  [3, 27, 55, 174],
  [0, 57, 97, 145, 152],
  [107, 128, 136, 160],
  [78, 106, 112, 135],
  [82, 106, 162, 169],
  [2, 9, 34, 175],
  [95, 99, 171, 172],
  [78, 136, 160, 168],
  [54, 87, 114, 116, 154],
  [78, 112, 121, 124],
  [75, 107, 128, 161],
  [2, 55, 174, 175],
  [0, 2, 39, 145],
  [2, 25, 58, 175],
  [78, 116, 123, 136],
  [79, 111, 166, 174],
  [37, 42, 114, 120],
  [3, 21, 79, 105],
  [8, 45, 92, 155, 159],
  [0, 83, 101, 175],
  [45, 82, 116, 163],
  [2, 18, 19, 35],
  [53, 79, 155, 160],
  [42, 81, 114, 129],
  [69, 116, 132, 163],
  [17, 58, 119, 141],
  [16, 36, 120, 148],
  [112, 114, 129, 144],
  [42, 81, 120, 148],
  [2, 4, 23, 65],
  [116, 137, 154, 163],
  [42, 81, 88, 109],
  [47, 80, 140, 173],
  [42, 81, 158, 167],
  [16, 36, 47, 80],
  [16, 36, 88, 109],
  [59, 117, 155, 168],
  [16, 36, 112, 144],
  [0, 21, 28, 174],
  [36, 80, 127, 130]],
 'Glycolysis': [[5, 38, 77, 136],
  [8, 22, 33, 91],
  [8, 16, 91, 109],
  [7, 21, 91, 103],
  [41, 57, 62, 133],
  [27, 36, 79, 85],
  [84, 100, 101, 135],
  [21, 35, 52, 100],
  [48, 84, 101, 261],
  [7, 38, 132, 136],
  [32, 45, 116, 263],
  [26, 27, 78, 132],
  [26, 81, 104, 106],
  [8, 80, 82, 91],
  [46, 49, 81, 116, 117, 263, 264],
  [22, 36, 42, 136],
  [27, 43, 79, 116],
  [22, 36, 79, 91],
  [32, 45, 81, 125],
  [8, 91, 102, 145],
  [51, 57, 62, 100, 135, 263],
  [26, 41, 105, 106],
  [27, 32, 42, 79],
  [51, 61, 97, 100],
  [26, 32, 38, 132],
  [8, 91, 126, 144],
  [29, 64, 104, 106],
  [22, 36, 102, 260],
  [27, 41, 79, 133],
  [25, 42, 62, 133],
  [32, 41, 45, 57],
  [26, 32, 39, 92],
  [25, 32, 42, 45],
  [32, 41, 73, 96],
  [48, 51, 100, 101],
  [21, 35, 37, 127],
  [13, 42, 96, 108],
  [25, 45, 89, 125],
  [5, 48, 77, 101],
  [45, 48, 51, 77, 124, 263],
  [5, 7, 21, 84, 135, 136],
  [8, 13, 23, 79, 91, 108],
  [7, 20, 21, 61],
  [27, 78, 79, 91],
  [7, 78, 91, 132],
  [32, 38, 77, 124],
  [148, 154, 164, 218, 247],
  [3, 5, 94, 111],
  [147, 171, 192, 221],
  [177, 178, 180, 235, 251],
  [113, 120, 137, 192, 251, 252],
  [7, 13, 58, 137],
  [148, 192, 240, 251],
  [153, 235, 236, 251],
  [7, 19, 132, 141],
  [158, 194, 225, 245, 247],
  [31, 113, 147, 254],
  [31, 93, 126, 147, 194],
  [13, 17, 118, 137],
  [21, 25, 195, 245, 265],
  [13, 15, 17, 108],
  [148, 159, 202, 251],
  [0, 7, 11, 35],
  [147, 194, 235, 245, 254],
  [21, 49, 51, 86],
  [183, 194, 208, 227, 247],
  [198, 204, 211, 220],
  [153, 167, 231, 232],
  [175, 192, 196, 202],
  [222, 238, 241, 259],
  [148, 159, 174, 246],
  [187, 198, 232, 243],
  [167, 198, 204, 232],
  [148, 167, 204, 246],
  [156, 167, 200, 232],
  [150, 154, 217, 218],
  [147, 154, 179, 193, 244],
  [187, 208, 237, 243, 247],
  [205, 208, 226, 243],
  [63, 135, 151, 245],
  [2, 68, 158, 163],
  [26, 40, 66, 142],
  [2, 68, 138, 158, 172, 214],
  [194, 208, 231, 243],
  [65, 113, 177, 192],
  [167, 184, 232, 246],
  [158, 179, 194, 199],
  [154, 183, 218, 233],
  [168, 178, 238, 243],
  [156, 184, 201, 232],
  [190, 208, 238, 243],
  [120, 170, 175, 176, 193, 253],
  [167, 221, 232, 250],
  [147, 153, 177, 193],
  [168, 178, 180, 236],
  [69, 196, 202, 229],
  [3, 29, 161, 182],
  [155, 186, 218, 250],
  [25, 50, 180, 189],
  [154, 181, 186, 218],
  [152, 161, 201, 225],
  [6, 200, 221, 228],
  [152, 154, 197, 219],
  [12, 248, 257, 259],
  [10, 67, 86, 258],
  [70, 113, 215, 252],
  [186, 196, 202, 257],
  [119, 161, 176, 210],
  [120, 170, 243, 251],
  [120, 170, 196, 242],
  [176, 187, 210, 253]],
 '3PGA_AMINOACIDS': [[64, 144, 148, 173],
  [20, 87, 92, 116],
  [2, 67, 70, 101],
  [16, 47, 81, 87],
  [23, 27, 66, 90],
  [105, 128, 140, 187],
  [27, 69, 73, 112],
  [1, 5, 80, 133],
  [44, 64, 164, 173],
  [27, 52, 66, 112],
  [1, 105, 131, 136],
  [71, 97, 129, 194],
  [16, 47, 71, 86],
  [16, 30, 39, 87],
  [1, 34, 65, 96],
  [43, 46, 53, 88],
  [44, 80, 161, 164],
  [6, 49, 91, 213],
  [0, 41, 77, 109],
  [3, 53, 88, 118],
  [26, 34, 49, 120],
  [29, 80, 161, 165],
  [80, 97, 129, 133],
  [1, 105, 157, 204],
  [1, 29, 105, 140],
  [10, 17, 41, 168],
  [0, 72, 77, 111],
  [16, 36, 87, 212],
  [35, 78, 84, 99],
  [18, 71, 123, 194],
  [3, 18, 45, 54, 77, 119],
  [1, 34, 105, 187],
  [23, 27, 78, 84],
  [10, 27, 84, 99],
  [1, 80, 85, 96],
  [1, 20, 96, 116],
  [10, 27, 43, 66],
  [23, 27, 39, 111],
  [6, 20, 154, 192],
  [1, 44, 64, 80],
  [21, 25, 33, 51],
  [1, 80, 119, 131],
  [10, 43, 53, 168],
  [3, 16, 18, 47],
  [16, 47, 85, 88],
  [2, 67, 99, 168],
  [23, 56, 57, 93],
  [1, 34, 80, 97],
  [85, 125, 132, 205],
  [1, 29, 131, 159],
  [64, 131, 141, 143, 159, 167, 173, 182],
  [6, 89, 95, 213],
  [3, 45, 86, 90],
  [27, 84, 93, 120],
  [17, 40, 41, 109],
  [16, 23, 71, 87],
  [18, 71, 97, 119],
  [105, 135, 144, 178, 190, 196],
  [1, 29, 131, 151, 186, 200, 206],
  [1, 49, 80, 94],
  [3, 45, 46, 88],
  [23, 81, 87, 90],
  [29, 124, 134, 141, 159, 167, 175, 184, 186, 200],
  [71, 97, 102, 103],
  [9, 72, 77, 119],
  [12, 16, 56, 87],
  [23, 91, 93, 213],
  [29, 80, 107, 200],
  [27, 79, 84, 112],
  [16, 61, 87, 213],
  [0, 37, 77, 100],
  [23, 87, 92, 93],
  [16, 20, 65, 71],
  [33, 83, 84, 120],
  [1, 6, 49, 96],
  [26, 57, 83, 91],
  [16, 62, 82, 87],
  [1, 64, 95, 96],
  [2, 15, 19, 21, 33, 67],
  [28, 52, 81, 90],
  [23, 27, 69, 89],
  [7, 67, 79, 84],
  [0, 27, 77, 93],
  [15, 17, 40, 94],
  [43, 46, 66, 90],
  [138, 144, 148, 187],
  [18, 36, 38, 119],
  [23, 27, 100, 212],
  [23, 34, 65, 93],
  [2, 16, 85, 168],
  [15, 19, 49, 94],
  [9, 18, 30, 119],
  [15, 80, 94, 168],
  [26, 40, 68, 94],
  [9, 30, 39, 72],
  [24, 35, 85, 88],
  [128, 137, 138, 187],
  [35, 71, 78, 85],
  [9, 40, 72, 109],
  [6, 95, 173, 192],
  [18, 40, 76, 119],
  [2, 15, 16, 61],
  [12, 16, 20, 50],
  [6, 15, 19, 61],
  [6, 19, 21, 50],
  [1, 64, 105, 144],
  [6, 16, 20, 61],
  [124, 134, 158, 184],
  [3, 45, 47, 81],
  [1, 96, 105, 190],
  [9, 37, 38, 72],
  [2, 16, 17, 76],
  [63, 103, 139, 143, 209],
  [11, 141, 150, 210],
  [143, 145, 195, 209],
  [14, 53, 122, 202],
  [64, 158, 167, 196],
  [11, 140, 189, 197, 210],
  [75, 149, 188, 201],
  [150, 166, 167, 176],
  [107, 146, 156, 200],
  [128, 141, 150, 185],
  [0, 59, 168, 169],
  [14, 19, 104, 147],
  [14, 19, 171, 179],
  [19, 143, 179, 192],
  [29, 97, 172, 175],
  [27, 52, 104, 127],
  [27, 104, 155, 180],
  [29, 93, 152, 184],
  [0, 179, 188, 208],
  [21, 51, 98, 110],
  [5, 21, 110, 124]]}
    return holes


def holes_pathway_pseudo():
    holes={'Glycolysis': [[0, 2, 4, 5], [0, 1, 2, 3]],
 'E4P_AMINO_ACIDS': [[0, 6, 7, 8],
  [0, 2, 4, 13],
  [0, 1, 3, 6],
  [2, 10, 13, 14],
  [0, 4, 5, 6],
  [0, 1, 8, 12]],
 'R5P_AMINOACIDS': [[2, 7, 9, 12],
  [0, 2, 5, 14],
  [0, 3, 15, 24],
  [0, 8, 23, 24],
  [0, 3, 8, 25],
  [0, 1, 2, 10],
  [0, 2, 4, 9],
  [0, 3, 5, 6]],
 'TCA': [[0, 3, 5, 8], [0, 2, 9, 15], [0, 2, 5, 13], [0, 1, 5, 11]],
 '3PGA_AMINOACIDS': [[1, 3, 7, 13],
  [0, 1, 4, 5],
  [2, 6, 16, 18],
  [1, 5, 9, 11],
  [1, 9, 13, 15],
  [1, 3, 8, 14],
  [2, 3, 4, 12],
  [1, 4, 10, 13],
  [1, 3, 5, 6]],
 'PYR_THR_AA': [[6, 14, 15, 81],
  [0, 6, 11, 81],
  [8, 19, 53, 86],
  [30, 37, 45, 64],
  [12, 49, 53, 86],
  [25, 30, 36, 37, 40, 52],
  [12, 26, 31, 58],
  [21, 24, 30, 75],
  [16, 46, 51, 105],
  [56, 61, 76, 82],
  [12, 33, 56, 86],
  [10, 24, 30, 43],
  [12, 22, 82, 86],
  [38, 39, 78, 85],
  [1, 32, 39, 78],
  [38, 50, 60, 85],
  [80, 89, 93, 100],
  [23, 72, 80, 95],
  [21, 63, 66, 79],
  [11, 69, 89, 95],
  [6, 24, 34, 36],
  [7, 72, 104, 105],
  [0, 6, 35, 99],
  [7, 35, 44, 88],
  [0, 32, 55, 98],
  [3, 14, 59, 65],
  [7, 20, 62, 88],
  [11, 29, 81, 84],
  [6, 24, 45, 81],
  [0, 11, 32, 84],
  [26, 75, 90, 101],
  [11, 15, 47, 84],
  [3, 6, 14, 44],
  [6, 39, 41, 71],
  [0, 28, 31, 34],
  [6, 80, 81, 93],
  [6, 23, 73, 93],
  [6, 73, 91, 101],
  [11, 59, 84, 103],
  [1, 6, 29, 81],
  [10, 24, 97, 101],
  [6, 14, 92, 93],
  [6, 15, 24, 64],
  [11, 59, 67, 95],
  [10, 23, 30, 100],
  [7, 23, 51, 72],
  [6, 11, 93, 95],
  [0, 6, 31, 75],
  [0, 1, 6, 32],
  [1, 2, 6, 73],
  [15, 27, 41, 47]],
 'ALPHAKETOGLUTARATE_AMINOACIDS': [[1, 44, 112, 114],
  [80, 97, 102, 121],
  [1, 61, 111, 115],
  [1, 35, 44, 49],
  [1, 35, 67, 76],
  [6, 37, 46, 73],
  [4, 6, 20, 30, 31, 46, 67, 88],
  [1, 38, 44, 115],
  [52, 93, 99, 109],
  [12, 30, 44, 67],
  [1, 56, 60, 133],
  [1, 110, 115, 133],
  [26, 60, 110, 133],
  [1, 24, 33, 35],
  [1, 17, 35, 111],
  [1, 16, 97, 133],
  [6, 44, 46, 99],
  [1, 6, 16, 46],
  [48, 51, 66, 139],
  [22, 77, 87, 101],
  [18, 34, 58, 101],
  [18, 71, 90, 101],
  [18, 90, 105, 124],
  [28, 45, 47, 74],
  [22, 28, 74, 101],
  [18, 21, 28, 90],
  [59, 64, 89, 96, 136],
  [66, 81, 98, 139],
  [39, 59, 107, 120],
  [28, 39, 77, 117],
  [25, 66, 122, 143],
  [0, 18, 93, 106, 143],
  [0, 44, 81, 86, 98, 122],
  [11, 20, 30, 116],
  [12, 36, 38, 117],
  [26, 36, 39, 59],
  [0, 7, 41, 122],
  [7, 64, 68, 97],
  [4, 10, 17, 67],
  [9, 23, 64, 97],
  [18, 98, 107, 122, 143],
  [11, 15, 20, 31],
  [2, 35, 121, 133],
  [21, 39, 59, 126],
  [4, 36, 91, 106],
  [4, 6, 75, 83]],
 'OXALACETATE_AMINOACIDS': [[48, 53, 84, 101],
  [43, 53, 98, 101],
  [30, 37, 68, 69],
  [20, 39, 54, 86],
  [14, 42, 48, 80],
  [0, 6, 17, 81],
  [4, 5, 12, 20, 36, 38],
  [4, 39, 74, 81],
  [11, 20, 30, 37],
  [20, 49, 54, 80],
  [11, 24, 27, 30, 46, 50, 82],
  [0, 27, 75, 94],
  [24, 50, 67, 87],
  [20, 36, 42, 80],
  [18, 22, 40, 52],
  [1, 35, 52, 65],
  [0, 1, 92, 94],
  [0, 74, 75, 81],
  [25, 42, 80, 90],
  [10, 35, 43, 80],
  [20, 35, 36, 85],
  [18, 48, 51, 67],
  [25, 43, 45, 90],
  [1, 4, 41, 56],
  [1, 24, 67, 92],
  [4, 20, 21, 81],
  [1, 4, 81, 92],
  [2, 6, 17, 29],
  [7, 24, 39, 89],
  [18, 35, 48, 80],
  [10, 24, 55, 64],
  [25, 43, 51, 53],
  [0, 21, 22, 81],
  [4, 20, 23, 28, 32, 36, 77, 94],
  [1, 4, 24, 39],
  [1, 18, 35, 92],
  [20, 37, 41, 95],
  [1, 11, 35, 94],
  [1, 4, 10, 97],
  [0, 18, 22, 92],
  [0, 11, 22, 94],
  [0, 67, 75, 92],
  [1, 4, 20, 35],
  [10, 24, 25, 43],
  [1, 24, 27, 94],
  [24, 25, 51, 67],
  [1, 7, 24, 65],
  [20, 36, 37, 69],
  [0, 73, 75, 104],
  [6, 31, 47, 93],
  [6, 34, 58, 66],
  [31, 40, 79, 93],
  [13, 28, 68, 102],
  [13, 34, 58, 60]]}
    return holes