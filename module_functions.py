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
def drop_columns_same_values(df):
    nunique = df.nunique()
    cols_to_drop = nunique[nunique == 1].index

    return df.drop(cols_to_drop, axis=1)
################################################################################################################################################################
def get_df_by_pathway_make_hole_droped_cols(df_by_pathway_drop_duplicate,holes_by_pathway):
    df_by_pathway_make_hole={}
    for path,df in df_by_pathway_drop_duplicate.items():
        holes_per_path=[]
        for hole in holes_by_pathway[path]:
            query=[]

            for vertex in hole:
                query.append(vertex)
            
            holes_per_path.append(drop_columns_same_values(df.loc[query,]))
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






def get_strepto_genus_path():
    dicc={'Streptosporangium_Glycolysis': [['6666666.1069482', '6666666.1069553', '6666666.1069569', '6666666.1069545']], 'Nonomuraea_Glycolysis': [['6666666.1069615', '6666666.1069571', '6666666.1069558', '6666666.1069633', '6666666.1069512'], ['6666666.1069524', '6666666.1069514', '6666666.1069500', '6666666.1069474'], ['6666666.1069479', '6666666.1069575', '6666666.1069477', '6666666.1069457'], ['6666666.1069550', '6666666.1069500', '6666666.1069474', '6666666.1069493'], ['6666666.1069558', '6666666.1069477', '6666666.1069547', '6666666.1069486'], ['6666666.1069550', '6666666.1069500', '6666666.1069514', '6666666.1069492'], ['6666666.1069444', '6666666.1069620', '6666666.1069547', '6666666.1069477'], ['6666666.1069550', '6666666.1069576', '6666666.1069536', '6666666.1069492'], ['6666666.1069558', '6666666.1069550', '6666666.1069493', '6666666.1069486']], 'Microbispora_Glycolysis': [['6666666.1069585', '6666666.1069449', '6666666.1069540', '6666666.1069503'], ['6666666.1069585', '6666666.1069584', '6666666.1069506', '6666666.1069503']], 'Herbidospora_Glycolysis': [['6666666.1069469', '6666666.1069521', '6666666.1069467', '6666666.1069451']], 'Sphaerisporangium_Glycolysis': [['6666666.1069639', '6666666.1069635', '6666666.1069602', '6666666.1069548'], ['6666666.1069602', '6666666.1069613', '6666666.1069488', '6666666.1069548']], 'Streptosporangiaceae_Glycolysis': [['6666666.1069676', '6666666.1069748', '6666666.1069688', '6666666.1069686'], ['6666666.1069740', '6666666.1069711', '6666666.1069733', '6666666.1069688'], ['6666666.1069744', '6666666.1069664', '6666666.1069711', '6666666.1069675'], ['6666666.1069692', '6666666.1069744', '6666666.1069664', '6666666.1069712'], ['6666666.1069694', '6666666.1069689', '6666666.1069662', '6666666.1069729'], ['6666666.1069662', '6666666.1069655', '6666666.1069701', '6666666.1069758', '6666666.1069687'], ['6666666.1069761', '6666666.1069749', '6666666.1069695', '6666666.1069757', '6666666.1069719'], ['6666666.1069721', '6666666.1069695', '6666666.1069767', '6666666.1069684'], ['6666666.1069656', '6666666.1069760', '6666666.1069675', '6666666.1069715'], ['6666666.1069694', '6666666.1069771', '6666666.1069706', '6666666.1069713'], ['6666666.1069706', '6666666.1069756', '6666666.1069678', '6666666.1069617'], ['6666666.1069656', '6666666.1069662', '6666666.1069672', '6666666.1069761', '6666666.1069729'], ['6666666.1069765', '6666666.1069700', '6666666.1069753', '6666666.1069656'], ['6666666.1069700', '6666666.1069733', '6666666.1069679', '6666666.1069655'], ['6666666.1069765', '6666666.1069747', '6666666.1069748', '6666666.1069661'], ['6666666.1069683', '6666666.1069701', '6666666.1069678', '6666666.1069767', '6666666.1069684', '6666666.1069617'], ['6666666.1069755', '6666666.1069773', '6666666.1069750', '6666666.1069734'], ['6666666.1069761', '6666666.1069690', '6666666.1069619', '6666666.1069655'], ['6666666.1069685', '6666666.1069701', '6666666.1069655', '6666666.1069661'], ['6666666.1069750', '6666666.1069757', '6666666.1069676', '6666666.1069686'], ['6666666.1069694', '6666666.1069663', '6666666.1069764', '6666666.1069729'], ['6666666.1069678', '6666666.1069757', '6666666.1069765', '6666666.1069617'], ['6666666.1069655', '6666666.1069768', '6666666.1069747', '6666666.1069759', '6666666.1069704'], ['6666666.1069710', '6666666.1069704', '6666666.1069666', '6666666.1069687'], ['6666666.1069737', '6666666.1069712', '6666666.1069669', '6666666.1069660'], ['6666666.1069760', '6666666.1069692', '6666666.1069744', '6666666.1069675'], ['6666666.1069757', '6666666.1069750', '6666666.1069698', '6666666.1069719'], ['6666666.1069743', '6666666.1069757', '6666666.1069719', '6666666.1069704'], ['6666666.1069757', '6666666.1069695', '6666666.1069709', '6666666.1069744'], ['6666666.1069721', '6666666.1069669', '6666666.1069616', '6666666.1069684'], ['6666666.1069709', '6666666.1069744', '6666666.1069675', '6666666.1069715'], ['6666666.1069731', '6666666.1069722', '6666666.1069709', '6666666.1069715'], ['6666666.1069691', '6666666.1069745', '6666666.1069662', '6666666.1069729'], ['6666666.1069708', '6666666.1069730', '6666666.1069662', '6666666.1069660'], ['6666666.1069759', '6666666.1069666', '6666666.1069737', '6666666.1069761', '6666666.1069704'], ['6666666.1069688', '6666666.1069747', '6666666.1069765', '6666666.1069685', '6666666.1069686'], ['6666666.1069765', '6666666.1069656', '6666666.1069713', '6666666.1069667'], ['6666666.1069728', '6666666.1069658', '6666666.1069662', '6666666.1069729'], ['6666666.1069706', '6666666.1069683', '6666666.1069700', '6666666.1069713'], ['6666666.1069744', '6666666.1069764', '6666666.1069733', '6666666.1069675'], ['6666666.1069719', '6666666.1069739', '6666666.1069691', '6666666.1069761', '6666666.1069704'], ['6666666.1069757', '6666666.1069738', '6666666.1069716', '6666666.1069719'], ['6666666.1069682', '6666666.1069760', '6666666.1069656', '6666666.1069667'], ['6666666.1069743', '6666666.1069744', '6666666.1069675', '6666666.1069661']], 'Nonomuraea_TCA': [['6666666.1069502', '6666666.1069498', '6666666.1069491', '6666666.1069457'], ['6666666.1069496', '6666666.1069620', '6666666.1069563', '6666666.1069784', '6666666.1069444', '6666666.1069558', '6666666.1069492'], ['6666666.1069502', '6666666.1069474', '6666666.1069492', '6666666.1069484'], ['6666666.1069501', '6666666.1069620', '6666666.1069496', '6666666.1069492'], ['6666666.1069538', '6666666.1069524', '6666666.1069502', '6666666.1069498'], ['6666666.1069494', '6666666.1069510', '6666666.1069492', '6666666.1069484'], ['6666666.1069535', '6666666.1069536', '6666666.1069444', '6666666.1069492'], ['6666666.1069514', '6666666.1069501', '6666666.1069492', '6666666.1069484'], ['6666666.1069535', '6666666.1069605', '6666666.1069547', '6666666.1069492'], ['6666666.1069535', '6666666.1069628', '6666666.1069501', '6666666.1069492'], ['6666666.1069614', '6666666.1069615', '6666666.1069605', '6666666.1069535'], ['6666666.1069491', '6666666.1069512', '6666666.1069484', '6666666.1069457'], ['6666666.1069784', '6666666.1069444', '6666666.1069536', '6666666.1069459'], ['6666666.1069547', '6666666.1069550', '6666666.1069501', '6666666.1069492'], ['6666666.1069512', '6666666.1069496', '6666666.1069492', '6666666.1069484'], ['6666666.1069550', '6666666.1069636', '6666666.1069514', '6666666.1069501']], 'Microbispora_TCA': [['6666666.1069574', '6666666.1069449', '6666666.1069506', '6666666.1069473'], ['6666666.1069626', '6666666.1069574', '6666666.1069473', '6666666.1069460']], 'Herbidospora_TCA': [['6666666.1069467', '6666666.1069469', '6666666.1069466', '6666666.1069451']], 'Sphaerisporangium_TCA': [['6666666.1069564', '6666666.1069613', '6666666.1069548', '6666666.1069602']], 'Streptosporangiaceae_TCA': [['6666666.1069736', '6666666.1069765', '6666666.1069698', '6666666.1069710'], ['6666666.1069701', '6666666.1069746', '6666666.1069665', '6666666.1069618'], ['6666666.1069693', '6666666.1069698', '6666666.1069710', '6666666.1069675', '6666666.1069616', '6666666.1069618'], ['6666666.1069756', '6666666.1069684', '6666666.1069689', '6666666.1069723', '6666666.1069762'], ['6666666.1069657', '6666666.1069746', '6666666.1069701', '6666666.1069616'], ['6666666.1069671', '6666666.1069713', '6666666.1069640', '6666666.1069616'], ['6666666.1069724', '6666666.1069682', '6666666.1069742', '6666666.1069752', '6666666.1069677'], ['6666666.1069676', '6666666.1069691', '6666666.1069747', '6666666.1069656'], ['6666666.1069752', '6666666.1069678', '6666666.1069751', '6666666.1069677'], ['6666666.1069718', '6666666.1069766', '6666666.1069691', '6666666.1069676'], ['6666666.1069728', '6666666.1069691', '6666666.1069766', '6666666.1069694'], ['6666666.1069729', '6666666.1069748', '6666666.1069694', '6666666.1069677'], ['6666666.1069700', '6666666.1069727', '6666666.1069701', '6666666.1069670', '6666666.1069660', '6666666.1069722', '6666666.1069616'], ['6666666.1069674', '6666666.1069716', '6666666.1069617', '6666666.1069664'], ['6666666.1069719', '6666666.1069645', '6666666.1069607', '6666666.1069617'], ['6666666.1069718', '6666666.1069676', '6666666.1069689', '6666666.1069724'], ['6666666.1069681', '6666666.1069683', '6666666.1069731', '6666666.1069682', '6666666.1069739', '6666666.1069696'], ['6666666.1069752', '6666666.1069730', '6666666.1069718', '6666666.1069677'], ['6666666.1069764', '6666666.1069678', '6666666.1069752', '6666666.1069742'], ['6666666.1069740', '6666666.1069668', '6666666.1069664', '6666666.1069760'], ['6666666.1069750', '6666666.1069771', '6666666.1069678', '6666666.1069742'], ['6666666.1069771', '6666666.1069768', '6666666.1069746', '6666666.1069701'], ['6666666.1069671', '6666666.1069695', '6666666.1069701', '6666666.1069727', '6666666.1069655'], ['6666666.1069660', '6666666.1069723', '6666666.1069706', '6666666.1069732', '6666666.1069675', '6666666.1069616'], ['6666666.1069757', '6666666.1069678', '6666666.1069742', '6666666.1069696'], ['6666666.1069771', '6666666.1069750', '6666666.1069733', '6666666.1069717'], ['6666666.1069755', '6666666.1069684', '6666666.1069682', '6666666.1069724'], ['6666666.1069758', '6666666.1069691', '6666666.1069766', '6666666.1069694'], ['6666666.1069739', '6666666.1069682', '6666666.1069737', '6666666.1069765', '6666666.1069736', '6666666.1069616'], ['6666666.1069662', '6666666.1069744', '6666666.1069684', '6666666.1069756', '6666666.1069663'], ['6666666.1069769', '6666666.1069680', '6666666.1069731', '6666666.1069683', '6666666.1069681', '6666666.1069735'], ['6666666.1069707', '6666666.1069757', '6666666.1069747', '6666666.1069719', '6666666.1069656', '6666666.1069696'], ['6666666.1069735', '6666666.1069750', '6666666.1069701', '6666666.1069616'], ['6666666.1069664', '6666666.1069705', '6666666.1069673', '6666666.1069720', '6666666.1069662', '6666666.1069742'], ['6666666.1069722', '6666666.1069713', '6666666.1069640', '6666666.1069616'], ['6666666.1069682', '6666666.1069731', '6666666.1069747', '6666666.1069719', '6666666.1069656', '6666666.1069739'], ['6666666.1069697', '6666666.1069733', '6666666.1069688', '6666666.1069724'], ['6666666.1069700', '6666666.1069727', '6666666.1069701', '6666666.1069663', '6666666.1069662', '6666666.1069640', '6666666.1069616'], ['6666666.1069731', '6666666.1069689', '6666666.1069723', '6666666.1069682'], ['6666666.1069675', '6666666.1069736', '6666666.1069640', '6666666.1069616'], ['6666666.1069707', '6666666.1069757', '6666666.1069733', '6666666.1069750', '6666666.1069696'], ['6666666.1069722', '6666666.1069720', '6666666.1069673', '6666666.1069704', '6666666.1069656', '6666666.1069696']], 'Streptosporangium_3PGA_AMINOACIDS': [['6666666.1069572', '6666666.1069553', '6666666.1069476', '6666666.1069455']], 'Nonomuraea_3PGA_AMINOACIDS': [['6666666.1069510', '6666666.1069535', '6666666.1069484', '6666666.1069444'], ['6666666.1069459', '6666666.1069558', '6666666.1069484', '6666666.1069444'], ['6666666.1069614', '6666666.1069633', '6666666.1069554', '6666666.1069491'], ['6666666.1069554', '6666666.1069491', '6666666.1069493', '6666666.1069454'], ['6666666.1069528', '6666666.1069529', '6666666.1069782', '6666666.1069491'], ['6666666.1069554', '6666666.1069494', '6666666.1069479', '6666666.1069454'], ['6666666.1069627', '6666666.1069500', '6666666.1069614', '6666666.1069491'], ['6666666.1069496', '6666666.1069562', '6666666.1069479', '6666666.1069454'], ['6666666.1069493', '6666666.1069490', '6666666.1069479', '6666666.1069454'], ['6666666.1069494', '6666666.1069530', '6666666.1069457', '6666666.1069479'], ['6666666.1069605', '6666666.1069535', '6666666.1069479', '6666666.1069454'], ['6666666.1069636', '6666666.1069529', '6666666.1069558', '6666666.1069484'], ['6666666.1069615', '6666666.1069474', '6666666.1069459', '6666666.1069444'], ['6666666.1069538', '6666666.1069575', '6666666.1069494', '6666666.1069479'], ['6666666.1069494', '6666666.1069782', '6666666.1069490', '6666666.1069479'], ['6666666.1069782', '6666666.1069491', '6666666.1069493', '6666666.1069490'], ['6666666.1069782', '6666666.1069530', '6666666.1069531', '6666666.1069491'], ['6666666.1069627', '6666666.1069547', '6666666.1069576', '6666666.1069531'], ['6666666.1069554', '6666666.1069627', '6666666.1069496', '6666666.1069454']], 'Microbispora_3PGA_AMINOACIDS': [['6666666.1069542', '6666666.1069608', '6666666.1069541', '6666666.1069485'], ['6666666.1069566', '6666666.1069458', '6666666.1069450', '6666666.1069485'], ['6666666.1069586', '6666666.1069526', '6666666.1069503', '6666666.1069507'], ['6666666.1069507', '6666666.1069609', '6666666.1069485', '6666666.1069566', '6666666.1069506'], ['6666666.1069508', '6666666.1069574', '6666666.1069506', '6666666.1069445'], ['6666666.1069586', '6666666.1069507', '6666666.1069506', '6666666.1069445']], 'Sphaerisporangium_3PGA_AMINOACIDS': [['6666666.1069613', '6666666.1069488', '6666666.1069557', '6666666.1069635', '6666666.1069548']], 'Streptosporangiaceae_3PGA_AMINOACIDS': [['6666666.1069727', '6666666.1069711', '6666666.1069658', '6666666.1069607'], ['6666666.1069688', '6666666.1069690', '6666666.1069696', '6666666.1069618'], ['6666666.1069744', '6666666.1069709', '6666666.1069705', '6666666.1069673'], ['6666666.1069671', '6666666.1069642', '6666666.1069695', '6666666.1069719', '6666666.1069709', '6666666.1069744', '6666666.1069646'], ['6666666.1069724', '6666666.1069711', '6666666.1069712', '6666666.1069691'], ['6666666.1069607', '6666666.1069713', '6666666.1069729', '6666666.1069663', '6666666.1069617'], ['6666666.1069646', '6666666.1069684', '6666666.1069772', '6666666.1069680', '6666666.1069616'], ['6666666.1069688', '6666666.1069656', '6666666.1069674', '6666666.1069619'], ['6666666.1069733', '6666666.1069700', '6666666.1069675', '6666666.1069661'], ['6666666.1069737', '6666666.1069685', '6666666.1069689', '6666666.1069679'], ['6666666.1069737', '6666666.1069667', '6666666.1069681', '6666666.1069619'], ['6666666.1069674', '6666666.1069673', '6666666.1069657', '6666666.1069655'], ['6666666.1069682', '6666666.1069691', '6666666.1069735', '6666666.1069667'], ['6666666.1069669', '6666666.1069750', '6666666.1069617', '6666666.1069616'], ['6666666.1069728', '6666666.1069716', '6666666.1069688', '6666666.1069618'], ['6666666.1069657', '6666666.1069705', '6666666.1069697', '6666666.1069676', '6666666.1069619'], ['6666666.1069712', '6666666.1069711', '6666666.1069750', '6666666.1069616'], ['6666666.1069737', '6666666.1069679', '6666666.1069678', '6666666.1069667'], ['6666666.1069692', '6666666.1069756', '6666666.1069700', '6666666.1069661'], ['6666666.1069728', '6666666.1069684', '6666666.1069679', '6666666.1069618'], ['6666666.1069692', '6666666.1069720', '6666666.1069729', '6666666.1069713', '6666666.1069657'], ['6666666.1069737', '6666666.1069664', '6666666.1069657', '6666666.1069619'], ['6666666.1069673', '6666666.1069657', '6666666.1069671', '6666666.1069646'], ['6666666.1069699', '6666666.1069766', '6666666.1069657', '6666666.1069619'], ['6666666.1069701', '6666666.1069682', '6666666.1069712', '6666666.1069731', '6666666.1069684', '6666666.1069719', '6666666.1069709', '6666666.1069744', '6666666.1069646'], ['6666666.1069669', '6666666.1069693', '6666666.1069700', '6666666.1069712', '6666666.1069616'], ['6666666.1069768', '6666666.1069771', '6666666.1069740', '6666666.1069676'], ['6666666.1069677', '6666666.1069671', '6666666.1069657', '6666666.1069619'], ['6666666.1069701', '6666666.1069682', '6666666.1069712', '6666666.1069721', '6666666.1069675', '6666666.1069661', '6666666.1069733', '6666666.1069736', '6666666.1069692', '6666666.1069768', '6666666.1069671'], ['6666666.1069657', '6666666.1069642', '6666666.1069695', '6666666.1069719', '6666666.1069689', '6666666.1069685', '6666666.1069619'], ['6666666.1069697', '6666666.1069687', '6666666.1069710', '6666666.1069705'], ['6666666.1069772', '6666666.1069754', '6666666.1069686', '6666666.1069684'], ['6666666.1069697', '6666666.1069714', '6666666.1069713', '6666666.1069705'], ['6666666.1069682', '6666666.1069721', '6666666.1069718', '6666666.1069616'], ['6666666.1069646', '6666666.1069673', '6666666.1069674', '6666666.1069669', '6666666.1069750', '6666666.1069658', '6666666.1069607'], ['6666666.1069742', '6666666.1069676', '6666666.1069726', '6666666.1069756', '6666666.1069685', '6666666.1069619']], 'Nonomuraea_ALPHAKETOGLUTARATE_AMINOACIDS': [['6666666.1069605', '6666666.1069535', '6666666.1069492', '6666666.1069459'], ['6666666.1069652', '6666666.1069520', '6666666.1069782', '6666666.1069536'], ['6666666.1069490', '6666666.1069559', '6666666.1069459', '6666666.1069454'], ['6666666.1069457', '6666666.1069547', '6666666.1069605', '6666666.1069459'], ['6666666.1069504', '6666666.1069780', '6666666.1069547', '6666666.1069457'], ['6666666.1069550', '6666666.1069520', '6666666.1069484', '6666666.1069498'], ['6666666.1069623', '6666666.1069493', '6666666.1069459', '6666666.1069454'], ['6666666.1069530', '6666666.1069578', '6666666.1069633', '6666666.1069535'], ['6666666.1069633', '6666666.1069547', '6666666.1069605', '6666666.1069535'], ['6666666.1069623', '6666666.1069504', '6666666.1069571', '6666666.1069454']], 'Microbispora_ALPHAKETOGLUTARATE_AMINOACIDS': [['6666666.1069458', '6666666.1069456', '6666666.1069445', '6666666.1069449'], ['6666666.1069608', '6666666.1069445', '6666666.1069585', '6666666.1069485']], 'Herbidospora_ALPHAKETOGLUTARATE_AMINOACIDS': [['6666666.1069470', '6666666.1069467', '6666666.1069466', '6666666.1069451']], 'Sphaerisporangium_ALPHAKETOGLUTARATE_AMINOACIDS': [['6666666.1069635', '6666666.1069644', '6666666.1069613', '6666666.1069548'], ['6666666.1069613', '6666666.1069557', '6666666.1069488', '6666666.1069548']], 'Streptosporangiaceae_ALPHAKETOGLUTARATE_AMINOACIDS': [['6666666.1069753', '6666666.1069751', '6666666.1069718', '6666666.1069743', '6666666.1069705'], ['6666666.1069726', '6666666.1069736', '6666666.1069749', '6666666.1069696'], ['6666666.1069737', '6666666.1069689', '6666666.1069726', '6666666.1069686'], ['6666666.1069663', '6666666.1069735', '6666666.1069616', '6666666.1069671'], ['6666666.1069771', '6666666.1069772', '6666666.1069706', '6666666.1069707'], ['6666666.1069617', '6666666.1069657', '6666666.1069664', '6666666.1069658'], ['6666666.1069769', '6666666.1069755', '6666666.1069710', '6666666.1069739', '6666666.1069699'], ['6666666.1069708', '6666666.1069766', '6666666.1069768', '6666666.1069655'], ['6666666.1069642', '6666666.1069693', '6666666.1069666', '6666666.1069658'], ['6666666.1069735', '6666666.1069663', '6666666.1069683', '6666666.1069686'], ['6666666.1069691', '6666666.1069688', '6666666.1069683', '6666666.1069681'], ['6666666.1069642', '6666666.1069721', '6666666.1069655', '6666666.1069616'], ['6666666.1069663', '6666666.1069721', '6666666.1069655', '6666666.1069705'], ['6666666.1069721', '6666666.1069762', '6666666.1069687', '6666666.1069700'], ['6666666.1069668', '6666666.1069759', '6666666.1069646', '6666666.1069664'], ['6666666.1069769', '6666666.1069762', '6666666.1069663', '6666666.1069705'], ['6666666.1069727', '6666666.1069697', '6666666.1069721', '6666666.1069700'], ['6666666.1069607', '6666666.1069740', '6666666.1069707', '6666666.1069616'], ['6666666.1069742', '6666666.1069758', '6666666.1069699', '6666666.1069655'], ['6666666.1069678', '6666666.1069763', '6666666.1069619', '6666666.1069617'], ['6666666.1069694', '6666666.1069759', '6666666.1069668', '6666666.1069669'], ['6666666.1069734', '6666666.1069748', '6666666.1069721', '6666666.1069671'], ['6666666.1069684', '6666666.1069685', '6666666.1069663', '6666666.1069686'], ['6666666.1069664', '6666666.1069717', '6666666.1069729', '6666666.1069704'], ['6666666.1069681', '6666666.1069676', '6666666.1069754', '6666666.1069661', '6666666.1069700'], ['6666666.1069642', '6666666.1069721', '6666666.1069700', '6666666.1069686'], ['6666666.1069659', '6666666.1069767', '6666666.1069745', '6666666.1069669'], ['6666666.1069642', '6666666.1069697', '6666666.1069749', '6666666.1069736', '6666666.1069686'], ['6666666.1069721', '6666666.1069705', '6666666.1069661', '6666666.1069700'], ['6666666.1069720', '6666666.1069664', '6666666.1069616', '6666666.1069642', '6666666.1069697', '6666666.1069704'], ['6666666.1069743', '6666666.1069718', '6666666.1069732', '6666666.1069699'], ['6666666.1069700', '6666666.1069708', '6666666.1069687', '6666666.1069681'], ['6666666.1069728', '6666666.1069689', '6666666.1069727', '6666666.1069664'], ['6666666.1069753', '6666666.1069769', '6666666.1069755', '6666666.1069712', '6666666.1069655'], ['6666666.1069618', '6666666.1069660', '6666666.1069641', '6666666.1069760', '6666666.1069658'], ['6666666.1069661', '6666666.1069739', '6666666.1069710', '6666666.1069657', '6666666.1069664', '6666666.1069658'], ['6666666.1069736', '6666666.1069693', '6666666.1069683', '6666666.1069686'], ['6666666.1069741', '6666666.1069683', '6666666.1069663', '6666666.1069671'], ['6666666.1069749', '6666666.1069724', '6666666.1069663', '6666666.1069686'], ['6666666.1069753', '6666666.1069769', '6666666.1069685', '6666666.1069738', '6666666.1069705'], ['6666666.1069674', '6666666.1069723', '6666666.1069655', '6666666.1069658'], ['6666666.1069642', '6666666.1069643', '6666666.1069726', '6666666.1069736', '6666666.1069686'], ['6666666.1069739', '6666666.1069682', '6666666.1069698', '6666666.1069700'], ['6666666.1069751', '6666666.1069736', '6666666.1069642', '6666666.1069658'], ['6666666.1069706', '6666666.1069772', '6666666.1069659', '6666666.1069705'], ['6666666.1069695', '6666666.1069739', '6666666.1069661', '6666666.1069671'], ['6666666.1069721', '6666666.1069708', '6666666.1069747', '6666666.1069731', '6666666.1069655'], ['6666666.1069754', '6666666.1069773', '6666666.1069617', '6666666.1069671'], ['6666666.1069693', '6666666.1069761', '6666666.1069669', '6666666.1069658'], ['6666666.1069684', '6666666.1069766', '6666666.1069681', '6666666.1069686']], 'Streptosporangium_E4P_AMINO_ACIDS': [['6666666.1069452', '6666666.1069552', '6666666.1069482', '6666666.1069442']], 'Nonomuraea_E4P_AMINO_ACIDS': [['6666666.1069615', '6666666.1069514', '6666666.1069444', '6666666.1069459'], ['6666666.1069532', '6666666.1069490', '6666666.1069784', '6666666.1069459'], ['6666666.1069474', '6666666.1069527', '6666666.1069444', '6666666.1069514'], ['6666666.1069504', '6666666.1069576', '6666666.1069459', '6666666.1069536'], ['6666666.1069782', '6666666.1069474', '6666666.1069615', '6666666.1069563'], ['6666666.1069538', '6666666.1069550', '6666666.1069492', '6666666.1069496'], ['6666666.1069510', '6666666.1069558', '6666666.1069498', '6666666.1069492', '6666666.1069496', '6666666.1069459'], ['6666666.1069628', '6666666.1069491', '6666666.1069484', '6666666.1069510'], ['6666666.1069615', '6666666.1069559', '6666666.1069491', '6666666.1069459'], ['6666666.1069610', '6666666.1069623', '6666666.1069533', '6666666.1069528', '6666666.1069459'], ['6666666.1069779', '6666666.1069562', '6666666.1069622', '6666666.1069558']], 'Microbispora_E4P_AMINO_ACIDS': [['6666666.1069584', '6666666.1069515', '6666666.1069588', '6666666.1069473'], ['6666666.1069526', '6666666.1069515', '6666666.1069499', '6666666.1069485'], ['6666666.1069507', '6666666.1069584', '6666666.1069473', '6666666.1069449'], ['6666666.1069588', '6666666.1069445', '6666666.1069586', '6666666.1069449']], 'Microtetraspora_E4P_AMINO_ACIDS': [['6666666.1069464', '6666666.1069631', '6666666.1069776', '6666666.1069463']], 'Sphaerisporangium_E4P_AMINO_ACIDS': [['6666666.1069644', '6666666.1069654', '6666666.1069635', '6666666.1069602'], ['6666666.1069654', '6666666.1069613', '6666666.1069564', '6666666.1069635']], 'Streptosporangiaceae_E4P_AMINO_ACIDS': [['6666666.1069676', '6666666.1069690', '6666666.1069735', '6666666.1069617', '6666666.1069616', '6666666.1069657'], ['6666666.1069706', '6666666.1069673', '6666666.1069669', '6666666.1069665'], ['6666666.1069705', '6666666.1069744', '6666666.1069700', '6666666.1069675'], ['6666666.1069769', '6666666.1069717', '6666666.1069712', '6666666.1069671'], ['6666666.1069700', '6666666.1069719', '6666666.1069710', '6666666.1069640'], ['6666666.1069768', '6666666.1069745', '6666666.1069771', '6666666.1069710'], ['6666666.1069684', '6666666.1069687', '6666666.1069749', '6666666.1069655'], ['6666666.1069693', '6666666.1069755', '6666666.1069719', '6666666.1069710'], ['6666666.1069771', '6666666.1069694', '6666666.1069712', '6666666.1069645'], ['6666666.1069720', '6666666.1069696', '6666666.1069655', '6666666.1069640'], ['6666666.1069759', '6666666.1069716', '6666666.1069662', '6666666.1069677', '6666666.1069710'], ['6666666.1069746', '6666666.1069768', '6666666.1069669', '6666666.1069640'], ['6666666.1069751', '6666666.1069666', '6666666.1069669', '6666666.1069671'], ['6666666.1069768', '6666666.1069718', '6666666.1069682', '6666666.1069664'], ['6666666.1069700', '6666666.1069740', '6666666.1069661', '6666666.1069640'], ['6666666.1069755', '6666666.1069687', '6666666.1069714', '6666666.1069720'], ['6666666.1069675', '6666666.1069646', '6666666.1069661', '6666666.1069616', '6666666.1069640'], ['6666666.1069710', '6666666.1069771', '6666666.1069694', '6666666.1069691'], ['6666666.1069743', '6666666.1069714', '6666666.1069738', '6666666.1069669'], ['6666666.1069672', '6666666.1069730', '6666666.1069760', '6666666.1069640'], ['6666666.1069676', '6666666.1069700', '6666666.1069740', '6666666.1069690'], ['6666666.1069663', '6666666.1069725', '6666666.1069758', '6666666.1069677', '6666666.1069680', '6666666.1069771', '6666666.1069769'], ['6666666.1069762', '6666666.1069731', '6666666.1069712', '6666666.1069739'], ['6666666.1069693', '6666666.1069750', '6666666.1069664', '6666666.1069673'], ['6666666.1069685', '6666666.1069730', '6666666.1069760', '6666666.1069640'], ['6666666.1069656', '6666666.1069689', '6666666.1069697', '6666666.1069706', '6666666.1069640'], ['6666666.1069712', '6666666.1069731', '6666666.1069749', '6666666.1069736', '6666666.1069711', '6666666.1069655'], ['6666666.1069676', '6666666.1069690', '6666666.1069754', '6666666.1069762', '6666666.1069729', '6666666.1069752', '6666666.1069657'], ['6666666.1069716', '6666666.1069733', '6666666.1069698', '6666666.1069640'], ['6666666.1069756', '6666666.1069757', '6666666.1069725', '6666666.1069676'], ['6666666.1069728', '6666666.1069739', '6666666.1069759', '6666666.1069698', '6666666.1069691', '6666666.1069710', '6666666.1069669', '6666666.1069640'], ['6666666.1069753', '6666666.1069682', '6666666.1069724', '6666666.1069691'], ['6666666.1069732', '6666666.1069686', '6666666.1069731', '6666666.1069712', '6666666.1069655'], ['6666666.1069752', '6666666.1069686', '6666666.1069732', '6666666.1069655'], ['6666666.1069771', '6666666.1069680', '6666666.1069677', '6666666.1069719', '6666666.1069708', '6666666.1069768', '6666666.1069710'], ['6666666.1069732', '6666666.1069686', '6666666.1069714', '6666666.1069738', '6666666.1069696', '6666666.1069655'], ['6666666.1069707', '6666666.1069617', '6666666.1069616', '6666666.1069640'], ['6666666.1069753', '6666666.1069718', '6666666.1069726', '6666666.1069691'], ['6666666.1069676', '6666666.1069700', '6666666.1069723', '6666666.1069767', '6666666.1069729'], ['6666666.1069714', '6666666.1069687', '6666666.1069718', '6666666.1069664'], ['6666666.1069720', '6666666.1069766', '6666666.1069688', '6666666.1069640'], ['6666666.1069714', '6666666.1069687', '6666666.1069749', '6666666.1069655'], ['6666666.1069757', '6666666.1069687', '6666666.1069755', '6666666.1069719'], ['6666666.1069670', '6666666.1069683', '6666666.1069669', '6666666.1069640'], ['6666666.1069728', '6666666.1069684', '6666666.1069712', '6666666.1069655'], ['6666666.1069699', '6666666.1069695', '6666666.1069676', '6666666.1069657'], ['6666666.1069731', '6666666.1069762', '6666666.1069729', '6666666.1069686'], ['6666666.1069766', '6666666.1069716', '6666666.1069759', '6666666.1069640']], 'Streptosporangium_OXALACETATE_AMINOACIDS': [['6666666.1069452', '6666666.1069632', '6666666.1069455', '6666666.1069476'], ['6666666.1069552', '6666666.1069553', '6666666.1069455', '6666666.1069442']], 'Nonomuraea_OXALACETATE_AMINOACIDS': [['6666666.1069615', '6666666.1069459', '6666666.1069484', '6666666.1069633'], ['6666666.1069547', '6666666.1069457', '6666666.1069512', '6666666.1069454'], ['6666666.1069457', '6666666.1069492', '6666666.1069571', '6666666.1069454'], ['6666666.1069575', '6666666.1069554', '6666666.1069627', '6666666.1069784'], ['6666666.1069491', '6666666.1069535', '6666666.1069536', '6666666.1069637', '6666666.1069454'], ['6666666.1069457', '6666666.1069620', '6666666.1069784', '6666666.1069779', '6666666.1069512'], ['6666666.1069571', '6666666.1069558', '6666666.1069501', '6666666.1069454'], ['6666666.1069457', '6666666.1069562', '6666666.1069780', '6666666.1069454'], ['6666666.1069491', '6666666.1069512', '6666666.1069501', '6666666.1069780', '6666666.1069454'], ['6666666.1069474', '6666666.1069554', '6666666.1069484', '6666666.1069633'], ['6666666.1069610', '6666666.1069605', '6666666.1069559', '6666666.1069510'], ['6666666.1069576', '6666666.1069444', '6666666.1069636', '6666666.1069605', '6666666.1069558'], ['6666666.1069536', '6666666.1069520', '6666666.1069550', '6666666.1069633'], ['6666666.1069782', '6666666.1069524', '6666666.1069514', '6666666.1069491'], ['6666666.1069479', '6666666.1069498', '6666666.1069636', '6666666.1069605', '6666666.1069559'], ['6666666.1069484', '6666666.1069459', '6666666.1069554', '6666666.1069627', '6666666.1069494'], ['6666666.1069547', '6666666.1069576', '6666666.1069501', '6666666.1069454'], ['6666666.1069477', '6666666.1069528', '6666666.1069494', '6666666.1069558'], ['6666666.1069491', '6666666.1069559', '6666666.1069648', '6666666.1069496', '6666666.1069782', '6666666.1069550', '6666666.1069477', '6666666.1069528', '6666666.1069633'], ['6666666.1069562', '6666666.1069576', '6666666.1069494', '6666666.1069516']], 'Microbispora_OXALACETATE_AMINOACIDS': [['6666666.1069539', '6666666.1069626', '6666666.1069462', '6666666.1069460'], ['6666666.1069515', '6666666.1069499', '6666666.1069542', '6666666.1069485']], 'Sphaerisporangium_OXALACETATE_AMINOACIDS': [['6666666.1069601', '6666666.1069602', '6666666.1069488', '6666666.1069557', '6666666.1069556']], 'Streptosporangiaceae_OXALACETATE_AMINOACIDS': [['6666666.1069695', '6666666.1069743', '6666666.1069699', '6666666.1069765'], ['6666666.1069667', '6666666.1069665', '6666666.1069699', '6666666.1069736'], ['6666666.1069656', '6666666.1069672', '6666666.1069722', '6666666.1069752', '6666666.1069655'], ['6666666.1069658', '6666666.1069616', '6666666.1069617', '6666666.1069655'], ['6666666.1069751', '6666666.1069756', '6666666.1069655', '6666666.1069675'], ['6666666.1069695', '6666666.1069732', '6666666.1069701', '6666666.1069675'], ['6666666.1069762', '6666666.1069692', '6666666.1069732', '6666666.1069735'], ['6666666.1069715', '6666666.1069691', '6666666.1069758', '6666666.1069699', '6666666.1069765'], ['6666666.1069720', '6666666.1069708', '6666666.1069709', '6666666.1069710'], ['6666666.1069757', '6666666.1069678', '6666666.1069689', '6666666.1069718'], ['6666666.1069684', '6666666.1069754', '6666666.1069735', '6666666.1069718'], ['6666666.1069700', '6666666.1069675', '6666666.1069709', '6666666.1069710'], ['6666666.1069741', '6666666.1069719', '6666666.1069743', '6666666.1069664'], ['6666666.1069752', '6666666.1069707', '6666666.1069742', '6666666.1069709', '6666666.1069736'], ['6666666.1069756', '6666666.1069735', '6666666.1069694', '6666666.1069655'], ['6666666.1069689', '6666666.1069643', '6666666.1069640', '6666666.1069765'], ['6666666.1069747', '6666666.1069681', '6666666.1069692', '6666666.1069732', '6666666.1069701'], ['6666666.1069659', '6666666.1069680', '6666666.1069688', '6666666.1069733', '6666666.1069748'], ['6666666.1069657', '6666666.1069701', '6666666.1069669', '6666666.1069617', '6666666.1069655'], ['6666666.1069760', '6666666.1069722', '6666666.1069752', '6666666.1069669'], ['6666666.1069712', '6666666.1069727', '6666666.1069645', '6666666.1069655'], ['6666666.1069668', '6666666.1069725', '6666666.1069701', '6666666.1069657'], ['6666666.1069656', '6666666.1069699', '6666666.1069693', '6666666.1069617', '6666666.1069655'], ['6666666.1069746', '6666666.1069641', '6666666.1069713', '6666666.1069710'], ['6666666.1069761', '6666666.1069726', '6666666.1069617', '6666666.1069616'], ['6666666.1069659', '6666666.1069680', '6666666.1069763', '6666666.1069741', '6666666.1069664', '6666666.1069675'], ['6666666.1069673', '6666666.1069731', '6666666.1069753', '6666666.1069675'], ['6666666.1069701', '6666666.1069725', '6666666.1069749', '6666666.1069715'], ['6666666.1069712', '6666666.1069734', '6666666.1069673', '6666666.1069655'], ['6666666.1069731', '6666666.1069721', '6666666.1069733', '6666666.1069655'], ['6666666.1069669', '6666666.1069714', '6666666.1069771', '6666666.1069688', '6666666.1069748'], ['6666666.1069699', '6666666.1069731', '6666666.1069753', '6666666.1069675'], ['6666666.1069676', '6666666.1069772', '6666666.1069669', '6666666.1069710'], ['6666666.1069691', '6666666.1069768', '6666666.1069736', '6666666.1069715'], ['6666666.1069699', '6666666.1069656', '6666666.1069657', '6666666.1069640', '6666666.1069765'], ['6666666.1069764', '6666666.1069704', '6666666.1069714', '6666666.1069616'], ['6666666.1069669', '6666666.1069617', '6666666.1069693', '6666666.1069736'], ['6666666.1069656', '6666666.1069725', '6666666.1069751', '6666666.1069765'], ['6666666.1069762', '6666666.1069677', '6666666.1069754', '6666666.1069718'], ['6666666.1069673', '6666666.1069607', '6666666.1069765', '6666666.1069675'], ['6666666.1069717', '6666666.1069768', '6666666.1069708', '6666666.1069659'], ['6666666.1069682', '6666666.1069665', '6666666.1069706', '6666666.1069743', '6666666.1069699'], ['6666666.1069616', '6666666.1069744', '6666666.1069766', '6666666.1069655'], ['6666666.1069731', '6666666.1069678', '6666666.1069691', '6666666.1069699'], ['6666666.1069660', '6666666.1069670', '6666666.1069705', '6666666.1069655'], ['6666666.1069757', '6666666.1069759', '6666666.1069688', '6666666.1069678'], ['6666666.1069686', '6666666.1069682', '6666666.1069699', '6666666.1069715'], ['6666666.1069753', '6666666.1069698', '6666666.1069715', '6666666.1069675'], ['6666666.1069714', '6666666.1069729', '6666666.1069669', '6666666.1069616'], ['6666666.1069732', '6666666.1069754', '6666666.1069684', '6666666.1069751'], ['6666666.1069711', '6666666.1069689', '6666666.1069765', '6666666.1069675'], ['6666666.1069751', '6666666.1069766', '6666666.1069736', '6666666.1069765'], ['6666666.1069750', '6666666.1069763', '6666666.1069659', '6666666.1069710'], ['6666666.1069744', '6666666.1069729', '6666666.1069669', '6666666.1069616'], ['6666666.1069735', '6666666.1069752', '6666666.1069698', '6666666.1069753', '6666666.1069718'], ['6666666.1069728', '6666666.1069768', '6666666.1069736', '6666666.1069715'], ['6666666.1069769', '6666666.1069663', '6666666.1069709', '6666666.1069710'], ['6666666.1069743', '6666666.1069706', '6666666.1069664', '6666666.1069607'], ['6666666.1069709', '6666666.1069742', '6666666.1069682', '6666666.1069700', '6666666.1069710'], ['6666666.1069713', '6666666.1069618', '6666666.1069672', '6666666.1069730', '6666666.1069741', '6666666.1069719', '6666666.1069656'], ['6666666.1069747', '6666666.1069681', '6666666.1069676', '6666666.1069710'], ['6666666.1069745', '6666666.1069713', '6666666.1069656', '6666666.1069657'], ['6666666.1069707', '6666666.1069772', '6666666.1069720', '6666666.1069736'], ['6666666.1069752', '6666666.1069707', '6666666.1069667', '6666666.1069764', '6666666.1069617', '6666666.1069736'], ['6666666.1069708', '6666666.1069769', '6666666.1069683', '6666666.1069749'], ['6666666.1069720', '6666666.1069682', '6666666.1069700', '6666666.1069710'], ['6666666.1069689', '6666666.1069617', '6666666.1069693', '6666666.1069765'], ['6666666.1069685', '6666666.1069755', '6666666.1069737', '6666666.1069698'], ['6666666.1069772', '6666666.1069762', '6666666.1069718', '6666666.1069675'], ['6666666.1069699', '6666666.1069765', '6666666.1069736', '6666666.1069693'], ['6666666.1069675', '6666666.1069765', '6666666.1069736', '6666666.1069715'], ['6666666.1069645', '6666666.1069706', '6666666.1069752', '6666666.1069715'], ['6666666.1069655', '6666666.1069694', '6666666.1069733', '6666666.1069766', '6666666.1069748'], ['6666666.1069744', '6666666.1069755', '6666666.1069617', '6666666.1069616'], ['6666666.1069735', '6666666.1069697', '6666666.1069759', '6666666.1069757', '6666666.1069718'], ['6666666.1069679', '6666666.1069681', '6666666.1069747', '6666666.1069617'], ['6666666.1069735', '6666666.1069754', '6666666.1069720', '6666666.1069643'], ['6666666.1069723', '6666666.1069727', '6666666.1069645', '6666666.1069655'], ['6666666.1069711', '6666666.1069718', '6666666.1069684', '6666666.1069751', '6666666.1069675'], ['6666666.1069744', '6666666.1069666', '6666666.1069694', '6666666.1069616'], ['6666666.1069767', '6666666.1069662', '6666666.1069713', '6666666.1069699'], ['6666666.1069691', '6666666.1069752', '6666666.1069698', '6666666.1069749', '6666666.1069736', '6666666.1069715'], ['6666666.1069671', '6666666.1069655', '6666666.1069709', '6666666.1069748'], ['6666666.1069753', '6666666.1069666', '6666666.1069737', '6666666.1069715'], ['6666666.1069733', '6666666.1069755', '6666666.1069617', '6666666.1069748'], ['6666666.1069699', '6666666.1069678', '6666666.1069688', '6666666.1069733', '6666666.1069655'], ['6666666.1069689', '6666666.1069678', '6666666.1069739', '6666666.1069755', '6666666.1069617'], ['6666666.1069727', '6666666.1069767', '6666666.1069699', '6666666.1069607']], 'Nonomuraea_PPP': [['6666666.1069527', '6666666.1069495', '6666666.1069500', '6666666.1069457'], ['6666666.1069614', '6666666.1069636', '6666666.1069500', '6666666.1069457'], ['6666666.1069491', '6666666.1069502', '6666666.1069459', '6666666.1069444'], ['6666666.1069784', '6666666.1069528', '6666666.1069457', '6666666.1069444'], ['6666666.1069784', '6666666.1069621', '6666666.1069500', '6666666.1069444'], ['6666666.1069527', '6666666.1069492', '6666666.1069502', '6666666.1069459'], ['6666666.1069527', '6666666.1069779', '6666666.1069520', '6666666.1069457'], ['6666666.1069527', '6666666.1069779', '6666666.1069512', '6666666.1069459'], ['6666666.1069459', '6666666.1069527', '6666666.1069457', '6666666.1069444'], ['6666666.1069784', '6666666.1069484', '6666666.1069512', '6666666.1069444'], ['6666666.1069527', '6666666.1069492', '6666666.1069474', '6666666.1069457'], ['6666666.1069614', '6666666.1069554', '6666666.1069491', '6666666.1069444'], ['6666666.1069491', '6666666.1069474', '6666666.1069457', '6666666.1069444']], 'Herbidospora_PPP': [['6666666.1069470', '6666666.1069467', '6666666.1069466', '6666666.1069451']], 'Streptosporangiaceae_PPP': [['6666666.1069763', '6666666.1069665', '6666666.1069690', '6666666.1069759'], ['6666666.1069766', '6666666.1069756', '6666666.1069729', '6666666.1069664'], ['6666666.1069711', '6666666.1069715', '6666666.1069646', '6666666.1069607'], ['6666666.1069735', '6666666.1069684', '6666666.1069618', '6666666.1069645'], ['6666666.1069668', '6666666.1069704', '6666666.1069662', '6666666.1069655'], ['6666666.1069685', '6666666.1069681', '6666666.1069699', '6666666.1069659'], ['6666666.1069729', '6666666.1069716', '6666666.1069708', '6666666.1069664'], ['6666666.1069714', '6666666.1069752', '6666666.1069692', '6666666.1069664'], ['6666666.1069749', '6666666.1069756', '6666666.1069696', '6666666.1069665'], ['6666666.1069740', '6666666.1069694', '6666666.1069665', '6666666.1069690'], ['6666666.1069769', '6666666.1069732', '6666666.1069657', '6666666.1069655'], ['6666666.1069724', '6666666.1069698', '6666666.1069736', '6666666.1069668'], ['6666666.1069725', '6666666.1069759', '6666666.1069708', '6666666.1069655'], ['6666666.1069616', '6666666.1069706', '6666666.1069748', '6666666.1069692', '6666666.1069655'], ['6666666.1069616', '6666666.1069662', '6666666.1069618', '6666666.1069645', '6666666.1069617'], ['6666666.1069757', '6666666.1069721', '6666666.1069696', '6666666.1069661'], ['6666666.1069696', '6666666.1069762', '6666666.1069694', '6666666.1069665'], ['6666666.1069733', '6666666.1069773', '6666666.1069666', '6666666.1069691'], ['6666666.1069746', '6666666.1069727', '6666666.1069732', '6666666.1069657'], ['6666666.1069667', '6666666.1069722', '6666666.1069737', '6666666.1069698'], ['6666666.1069704', '6666666.1069728', '6666666.1069695', '6666666.1069668'], ['6666666.1069724', '6666666.1069749', '6666666.1069766', '6666666.1069668'], ['6666666.1069769', '6666666.1069616', '6666666.1069705', '6666666.1069655'], ['6666666.1069729', '6666666.1069683', '6666666.1069691', '6666666.1069668'], ['6666666.1069737', '6666666.1069722', '6666666.1069706', '6666666.1069704'], ['6666666.1069720', '6666666.1069741', '6666666.1069722', '6666666.1069737', '6666666.1069666'], ['6666666.1069667', '6666666.1069698', '6666666.1069675', '6666666.1069640'], ['6666666.1069675', '6666666.1069673', '6666666.1069672', '6666666.1069640'], ['6666666.1069671', '6666666.1069768', '6666666.1069646', '6666666.1069619'], ['6666666.1069724', '6666666.1069737', '6666666.1069704', '6666666.1069668'], ['6666666.1069772', '6666666.1069771', '6666666.1069687', '6666666.1069683'], ['6666666.1069704', '6666666.1069752', '6666666.1069692', '6666666.1069668'], ['6666666.1069773', '6666666.1069679', '6666666.1069724', '6666666.1069666'], ['6666666.1069759', '6666666.1069730', '6666666.1069748', '6666666.1069708'], ['6666666.1069763', '6666666.1069736', '6666666.1069762', '6666666.1069694'], ['6666666.1069686', '6666666.1069735', '6666666.1069688', '6666666.1069646'], ['6666666.1069684', '6666666.1069754', '6666666.1069667', '6666666.1069640'], ['6666666.1069686', '6666666.1069661', '6666666.1069675', '6666666.1069646'], ['6666666.1069750', '6666666.1069733', '6666666.1069657', '6666666.1069655'], ['6666666.1069688', '6666666.1069738', '6666666.1069671', '6666666.1069619'], ['6666666.1069768', '6666666.1069758', '6666666.1069695', '6666666.1069668'], ['6666666.1069708', '6666666.1069759', '6666666.1069690', '6666666.1069669'], ['6666666.1069677', '6666666.1069706', '6666666.1069721', '6666666.1069729', '6666666.1069657'], ['6666666.1069660', '6666666.1069740', '6666666.1069694', '6666666.1069762', '6666666.1069640'], ['6666666.1069750', '6666666.1069769', '6666666.1069699', '6666666.1069659'], ['6666666.1069750', '6666666.1069738', '6666666.1069678', '6666666.1069668'], ['6666666.1069713', '6666666.1069741', '6666666.1069667', '6666666.1069640'], ['6666666.1069756', '6666666.1069696', '6666666.1069721', '6666666.1069729'], ['6666666.1069671', '6666666.1069707', '6666666.1069760', '6666666.1069619']], 'Streptosporangium_PYR_THR_AA': [['6666666.1069787', '6666666.1069572', '6666666.1069569', '6666666.1069553'], ['6666666.1069651', '6666666.1069553', '6666666.1069455', '6666666.1069442'], ['6666666.1069545', '6666666.1069572', '6666666.1069785', '6666666.1069442']], 'Nonomuraea_PYR_THR_AA': [['6666666.1069527', '6666666.1069605', '6666666.1069622', '6666666.1069492'], ['6666666.1069496', '6666666.1069535', '6666666.1069648', '6666666.1069576', '6666666.1069510'], ['6666666.1069536', '6666666.1069627', '6666666.1069538', '6666666.1069492'], ['6666666.1069780', '6666666.1069501', '6666666.1069498', '6666666.1069454'], ['6666666.1069620', '6666666.1069610', '6666666.1069536', '6666666.1069510'], ['6666666.1069532', '6666666.1069528', '6666666.1069559', '6666666.1069504'], ['6666666.1069554', '6666666.1069559', '6666666.1069636', '6666666.1069454', '6666666.1069524', '6666666.1069504'], ['6666666.1069454', '6666666.1069780', '6666666.1069486', '6666666.1069784', '6666666.1069530'], ['6666666.1069558', '6666666.1069562', '6666666.1069536', '6666666.1069622'], ['6666666.1069524', '6666666.1069493', '6666666.1069516', '6666666.1069504'], ['6666666.1069648', '6666666.1069562', '6666666.1069536', '6666666.1069622']], 'Microbispora_PYR_THR_AA': [['6666666.1069448', '6666666.1069458', '6666666.1069506', '6666666.1069509'], ['6666666.1069584', '6666666.1069509', '6666666.1069507', '6666666.1069485'], ['6666666.1069588', '6666666.1069473', '6666666.1069445', '6666666.1069503']], 'Streptosporangiaceae_PYR_THR_AA': [['6666666.1069689', '6666666.1069687', '6666666.1069738', '6666666.1069743', '6666666.1069683'], ['6666666.1069619', '6666666.1069646', '6666666.1069729', '6666666.1069675'], ['6666666.1069677', '6666666.1069764', '6666666.1069743', '6666666.1069768', '6666666.1069616'], ['6666666.1069661', '6666666.1069710', '6666666.1069619', '6666666.1069657'], ['6666666.1069768', '6666666.1069742', '6666666.1069749', '6666666.1069685'], ['6666666.1069696', '6666666.1069768', '6666666.1069685', '6666666.1069753', '6666666.1069666'], ['6666666.1069739', '6666666.1069741', '6666666.1069619', '6666666.1069771'], ['6666666.1069766', '6666666.1069718', '6666666.1069680', '6666666.1069747'], ['6666666.1069747', '6666666.1069666', '6666666.1069674', '6666666.1069694'], ['6666666.1069641', '6666666.1069730', '6666666.1069767', '6666666.1069645', '6666666.1069617'], ['6666666.1069683', '6666666.1069756', '6666666.1069708', '6666666.1069675'], ['6666666.1069769', '6666666.1069685', '6666666.1069687', '6666666.1069690'], ['6666666.1069691', '6666666.1069749', '6666666.1069736', '6666666.1069675'], ['6666666.1069721', '6666666.1069699', '6666666.1069771', '6666666.1069712'], ['6666666.1069661', '6666666.1069746', '6666666.1069767', '6666666.1069616'], ['6666666.1069750', '6666666.1069720', '6666666.1069729', '6666666.1069751', '6666666.1069771', '6666666.1069681'], ['6666666.1069707', '6666666.1069733', '6666666.1069738', '6666666.1069709', '6666666.1069705', '6666666.1069693'], ['6666666.1069705', '6666666.1069673', '6666666.1069674', '6666666.1069694'], ['6666666.1069681', '6666666.1069771', '6666666.1069655', '6666666.1069657'], ['6666666.1069667', '6666666.1069762', '6666666.1069691', '6666666.1069683'], ['6666666.1069769', '6666666.1069735', '6666666.1069751', '6666666.1069719'], ['6666666.1069743', '6666666.1069738', '6666666.1069733', '6666666.1069707', '6666666.1069719', '6666666.1069663'], ['6666666.1069772', '6666666.1069678', '6666666.1069732', '6666666.1069676', '6666666.1069687'], ['6666666.1069681', '6666666.1069771', '6666666.1069699', '6666666.1069683', '6666666.1069675'], ['6666666.1069640', '6666666.1069710', '6666666.1069619', '6666666.1069616'], ['6666666.1069768', '6666666.1069685', '6666666.1069736', '6666666.1069675'], ['6666666.1069750', '6666666.1069720', '6666666.1069748', '6666666.1069675'], ['6666666.1069721', '6666666.1069643', '6666666.1069754', '6666666.1069712'], ['6666666.1069767', '6666666.1069713', '6666666.1069619', '6666666.1069645'], ['6666666.1069748', '6666666.1069720', '6666666.1069729', '6666666.1069685', '6666666.1069736', '6666666.1069675'], ['6666666.1069719', '6666666.1069707', '6666666.1069732', '6666666.1069678', '6666666.1069690'], ['6666666.1069738', '6666666.1069676', '6666666.1069732', '6666666.1069733'], ['6666666.1069757', '6666666.1069672', '6666666.1069670', '6666666.1069767', '6666666.1069645'], ['6666666.1069667', '6666666.1069646', '6666666.1069619', '6666666.1069706'], ['6666666.1069709', '6666666.1069738', '6666666.1069666', '6666666.1069673'], ['6666666.1069676', '6666666.1069738', '6666666.1069743', '6666666.1069718', '6666666.1069761', '6666666.1069717', '6666666.1069731'], ['6666666.1069738', '6666666.1069666', '6666666.1069696', '6666666.1069768', '6666666.1069733'], ['6666666.1069768', '6666666.1069733', '6666666.1069711', '6666666.1069700'], ['6666666.1069766', '6666666.1069738', '6666666.1069666', '6666666.1069747'], ['6666666.1069668', '6666666.1069674', '6666666.1069666', '6666666.1069696', '6666666.1069768', '6666666.1069742', '6666666.1069657'], ['6666666.1069683', '6666666.1069718', '6666666.1069761', '6666666.1069691', '6666666.1069659', '6666666.1069748', '6666666.1069675'], ['6666666.1069730', '6666666.1069762', '6666666.1069677', '6666666.1069641'], ['6666666.1069682', '6666666.1069760', '6666666.1069694', '6666666.1069664'], ['6666666.1069680', '6666666.1069760', '6666666.1069694', '6666666.1069713'], ['6666666.1069661', '6666666.1069658', '6666666.1069641', '6666666.1069617'], ['6666666.1069722', '6666666.1069738', '6666666.1069666', '6666666.1069747', '6666666.1069664'], ['6666666.1069728', '6666666.1069772', '6666666.1069689', '6666666.1069708'], ['6666666.1069665', '6666666.1069694', '6666666.1069674', '6666666.1069657'], ['6666666.1069768', '6666666.1069733', '6666666.1069738', '6666666.1069743', '6666666.1069683', '6666666.1069675'], ['6666666.1069768', '6666666.1069742', '6666666.1069681', '6666666.1069675'], ['6666666.1069766', '6666666.1069738', '6666666.1069743', '6666666.1069718'], ['6666666.1069720', '6666666.1069748', '6666666.1069754', '6666666.1069684'], ['6666666.1069762', '6666666.1069730', '6666666.1069727', '6666666.1069691'], ['6666666.1069669', '6666666.1069750', '6666666.1069681', '6666666.1069657'], ['6666666.1069743', '6666666.1069718', '6666666.1069724', '6666666.1069759', '6666666.1069690', '6666666.1069663'], ['6666666.1069748', '6666666.1069756', '6666666.1069708', '6666666.1069675'], ['6666666.1069645', '6666666.1069767', '6666666.1069670', '6666666.1069660', '6666666.1069616'], ['6666666.1069619', '6666666.1069756', '6666666.1069748', '6666666.1069659'], ['6666666.1069693', '6666666.1069707', '6666666.1069719', '6666666.1069690']], 'Streptosporangium_R5P_AMINOACIDS': [['6666666.1069632', '6666666.1069480', '6666666.1069452', '6666666.1069442'], ['6666666.1069482', '6666666.1069553', '6666666.1069480', '6666666.1069452']], 'Nonomuraea_R5P_AMINOACIDS': [['6666666.1069495', '6666666.1069457', '6666666.1069454', '6666666.1069444'], ['6666666.1069532', '6666666.1069633', '6666666.1069493', '6666666.1069454'], ['6666666.1069493', '6666666.1069512', '6666666.1069457', '6666666.1069454']], 'Microbispora_R5P_AMINOACIDS': [['6666666.1069509', '6666666.1069449', '6666666.1069473', '6666666.1069450']], 'Acrocarpospora_R5P_AMINOACIDS': [['6666666.1069612', '6666666.1069519', '6666666.1069518', '6666666.1069517']], 'Streptosporangiaceae_R5P_AMINOACIDS': [['6666666.1069743', '6666666.1069708', '6666666.1069729', '6666666.1069679'], ['6666666.1069756', '6666666.1069768', '6666666.1069754', '6666666.1069687'], ['6666666.1069664', '6666666.1069741', '6666666.1069759', '6666666.1069728'], ['6666666.1069696', '6666666.1069679', '6666666.1069663', '6666666.1069729'], ['6666666.1069643', '6666666.1069617', '6666666.1069718', '6666666.1069676', '6666666.1069756', '6666666.1069666', '6666666.1069616'], ['6666666.1069753', '6666666.1069664', '6666666.1069697', '6666666.1069704', '6666666.1069746', '6666666.1069718', '6666666.1069617', '6666666.1069643'], ['6666666.1069717', '6666666.1069748', '6666666.1069691', '6666666.1069694'], ['6666666.1069711', '6666666.1069685', '6666666.1069666', '6666666.1069655'], ['6666666.1069655', '6666666.1069705', '6666666.1069641', '6666666.1069616'], ['6666666.1069736', '6666666.1069738', '6666666.1069718', '6666666.1069617', '6666666.1069643'], ['6666666.1069709', '6666666.1069720', '6666666.1069693', '6666666.1069686', '6666666.1069676', '6666666.1069756', '6666666.1069678'], ['6666666.1069696', '6666666.1069755', '6666666.1069746', '6666666.1069718', '6666666.1069676'], ['6666666.1069718', '6666666.1069731', '6666666.1069694', '6666666.1069616'], ['6666666.1069736', '6666666.1069738', '6666666.1069745', '6666666.1069662', '6666666.1069616'], ['6666666.1069717', '6666666.1069730', '6666666.1069659', '6666666.1069643'], ['6666666.1069736', '6666666.1069726', '6666666.1069700', '6666666.1069616'], ['6666666.1069771', '6666666.1069692', '6666666.1069664', '6666666.1069704'], ['6666666.1069695', '6666666.1069757', '6666666.1069720', '6666666.1069678'], ['6666666.1069740', '6666666.1069617', '6666666.1069643', '6666666.1069618'], ['6666666.1069641', '6666666.1069642', '6666666.1069618', '6666666.1069616'], ['6666666.1069756', '6666666.1069768', '6666666.1069686', '6666666.1069676'], ['6666666.1069759', '6666666.1069728', '6666666.1069691', '6666666.1069694'], ['6666666.1069683', '6666666.1069668', '6666666.1069659', '6666666.1069643'], ['6666666.1069766', '6666666.1069690', '6666666.1069682', '6666666.1069616'], ['6666666.1069755', '6666666.1069729', '6666666.1069663', '6666666.1069762'], ['6666666.1069686', '6666666.1069721', '6666666.1069679', '6666666.1069696'], ['6666666.1069713', '6666666.1069660', '6666666.1069705', '6666666.1069641'], ['6666666.1069723', '6666666.1069682', '6666666.1069711', '6666666.1069655'], ['6666666.1069700', '6666666.1069666', '6666666.1069655', '6666666.1069616'], ['6666666.1069655', '6666666.1069711', '6666666.1069643', '6666666.1069616'], ['6666666.1069741', '6666666.1069664', '6666666.1069753', '6666666.1069616'], ['6666666.1069643', '6666666.1069713', '6666666.1069641', '6666666.1069616'], ['6666666.1069700', '6666666.1069701', '6666666.1069678', '6666666.1069666'], ['6666666.1069713', '6666666.1069660', '6666666.1069711', '6666666.1069643'], ['6666666.1069762', '6666666.1069725', '6666666.1069756', '6666666.1069768', '6666666.1069663'], ['6666666.1069755', '6666666.1069714', '6666666.1069696', '6666666.1069729'], ['6666666.1069712', '6666666.1069700', '6666666.1069701', '6666666.1069684'], ['6666666.1069686', '6666666.1069714', '6666666.1069727', '6666666.1069745', '6666666.1069738', '6666666.1069718', '6666666.1069676']], 'Thermobispora_Glycolysis': [], 'Microtetraspora_Glycolysis': [], 'Planomonospora_Glycolysis': [], 'Planobispora_Glycolysis': [], 'Thermoactinospora_Glycolysis': [], 'Spongiactinospora_Glycolysis': [], 'Bailinhaonella_Glycolysis': [], 'Thermopolyspora_Glycolysis': [], 'Acrocarpospora_Glycolysis': [], 'Thermocatellispora_Glycolysis': [], 'Planotetraspora_Glycolysis': [], 'Sphaerimonospora_Glycolysis': [], 'Sinosporangium_Glycolysis': [], 'Rhizohabitans_Glycolysis': [], 'Acidimicrobium_Glycolysis': [], 'Actinomadura_Glycolysis': [], 'Streptosporangium_\ufeffGlycolysis': [], 'Thermobispora_\ufeffGlycolysis': [], 'Nonomuraea_\ufeffGlycolysis': [], 'Microbispora_\ufeffGlycolysis': [], 'Herbidospora_\ufeffGlycolysis': [], 'Microtetraspora_\ufeffGlycolysis': [], 'Planomonospora_\ufeffGlycolysis': [], 'Planobispora_\ufeffGlycolysis': [], 'Thermoactinospora_\ufeffGlycolysis': [], 'Sphaerisporangium_\ufeffGlycolysis': [], 'Spongiactinospora_\ufeffGlycolysis': [], 'Bailinhaonella_\ufeffGlycolysis': [], 'Thermopolyspora_\ufeffGlycolysis': [], 'Acrocarpospora_\ufeffGlycolysis': [], 'Thermocatellispora_\ufeffGlycolysis': [], 'Planotetraspora_\ufeffGlycolysis': [], 'Sphaerimonospora_\ufeffGlycolysis': [], 'Sinosporangium_\ufeffGlycolysis': [], 'Streptosporangiaceae_\ufeffGlycolysis': [], 'Rhizohabitans_\ufeffGlycolysis': [], 'Acidimicrobium_\ufeffGlycolysis': [], 'Actinomadura_\ufeffGlycolysis': [], 'Streptosporangium_TCA': [], 'Thermobispora_TCA': [], 'Microtetraspora_TCA': [], 'Planomonospora_TCA': [], 'Planobispora_TCA': [], 'Thermoactinospora_TCA': [], 'Spongiactinospora_TCA': [], 'Bailinhaonella_TCA': [], 'Thermopolyspora_TCA': [], 'Acrocarpospora_TCA': [], 'Thermocatellispora_TCA': [], 'Planotetraspora_TCA': [], 'Sphaerimonospora_TCA': [], 'Sinosporangium_TCA': [], 'Rhizohabitans_TCA': [], 'Acidimicrobium_TCA': [], 'Actinomadura_TCA': [], 'Thermobispora_3PGA_AMINOACIDS': [], 'Herbidospora_3PGA_AMINOACIDS': [], 'Microtetraspora_3PGA_AMINOACIDS': [], 'Planomonospora_3PGA_AMINOACIDS': [], 'Planobispora_3PGA_AMINOACIDS': [], 'Thermoactinospora_3PGA_AMINOACIDS': [], 'Spongiactinospora_3PGA_AMINOACIDS': [], 'Bailinhaonella_3PGA_AMINOACIDS': [], 'Thermopolyspora_3PGA_AMINOACIDS': [], 'Acrocarpospora_3PGA_AMINOACIDS': [], 'Thermocatellispora_3PGA_AMINOACIDS': [], 'Planotetraspora_3PGA_AMINOACIDS': [], 'Sphaerimonospora_3PGA_AMINOACIDS': [], 'Sinosporangium_3PGA_AMINOACIDS': [], 'Rhizohabitans_3PGA_AMINOACIDS': [], 'Acidimicrobium_3PGA_AMINOACIDS': [], 'Actinomadura_3PGA_AMINOACIDS': [], 'Streptosporangium_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Thermobispora_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Microtetraspora_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Planomonospora_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Planobispora_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Thermoactinospora_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Spongiactinospora_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Bailinhaonella_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Thermopolyspora_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Acrocarpospora_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Thermocatellispora_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Planotetraspora_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Sphaerimonospora_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Sinosporangium_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Rhizohabitans_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Acidimicrobium_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Actinomadura_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Thermobispora_E4P_AMINO_ACIDS': [], 'Herbidospora_E4P_AMINO_ACIDS': [], 'Planomonospora_E4P_AMINO_ACIDS': [], 'Planobispora_E4P_AMINO_ACIDS': [], 'Thermoactinospora_E4P_AMINO_ACIDS': [], 'Spongiactinospora_E4P_AMINO_ACIDS': [], 'Bailinhaonella_E4P_AMINO_ACIDS': [], 'Thermopolyspora_E4P_AMINO_ACIDS': [], 'Acrocarpospora_E4P_AMINO_ACIDS': [], 'Thermocatellispora_E4P_AMINO_ACIDS': [], 'Planotetraspora_E4P_AMINO_ACIDS': [], 'Sphaerimonospora_E4P_AMINO_ACIDS': [], 'Sinosporangium_E4P_AMINO_ACIDS': [], 'Rhizohabitans_E4P_AMINO_ACIDS': [], 'Acidimicrobium_E4P_AMINO_ACIDS': [], 'Actinomadura_E4P_AMINO_ACIDS': [], 'Thermobispora_OXALACETATE_AMINOACIDS': [], 'Herbidospora_OXALACETATE_AMINOACIDS': [], 'Microtetraspora_OXALACETATE_AMINOACIDS': [], 'Planomonospora_OXALACETATE_AMINOACIDS': [], 'Planobispora_OXALACETATE_AMINOACIDS': [], 'Thermoactinospora_OXALACETATE_AMINOACIDS': [], 'Spongiactinospora_OXALACETATE_AMINOACIDS': [], 'Bailinhaonella_OXALACETATE_AMINOACIDS': [], 'Thermopolyspora_OXALACETATE_AMINOACIDS': [], 'Acrocarpospora_OXALACETATE_AMINOACIDS': [], 'Thermocatellispora_OXALACETATE_AMINOACIDS': [], 'Planotetraspora_OXALACETATE_AMINOACIDS': [], 'Sphaerimonospora_OXALACETATE_AMINOACIDS': [], 'Sinosporangium_OXALACETATE_AMINOACIDS': [], 'Rhizohabitans_OXALACETATE_AMINOACIDS': [], 'Acidimicrobium_OXALACETATE_AMINOACIDS': [], 'Actinomadura_OXALACETATE_AMINOACIDS': [], 'Streptosporangium_PPP': [], 'Thermobispora_PPP': [], 'Microbispora_PPP': [], 'Microtetraspora_PPP': [], 'Planomonospora_PPP': [], 'Planobispora_PPP': [], 'Thermoactinospora_PPP': [], 'Sphaerisporangium_PPP': [], 'Spongiactinospora_PPP': [], 'Bailinhaonella_PPP': [], 'Thermopolyspora_PPP': [], 'Acrocarpospora_PPP': [], 'Thermocatellispora_PPP': [], 'Planotetraspora_PPP': [], 'Sphaerimonospora_PPP': [], 'Sinosporangium_PPP': [], 'Rhizohabitans_PPP': [], 'Acidimicrobium_PPP': [], 'Actinomadura_PPP': [], 'Thermobispora_PYR_THR_AA': [], 'Herbidospora_PYR_THR_AA': [], 'Microtetraspora_PYR_THR_AA': [], 'Planomonospora_PYR_THR_AA': [], 'Planobispora_PYR_THR_AA': [], 'Thermoactinospora_PYR_THR_AA': [], 'Sphaerisporangium_PYR_THR_AA': [], 'Spongiactinospora_PYR_THR_AA': [], 'Bailinhaonella_PYR_THR_AA': [], 'Thermopolyspora_PYR_THR_AA': [], 'Acrocarpospora_PYR_THR_AA': [], 'Thermocatellispora_PYR_THR_AA': [], 'Planotetraspora_PYR_THR_AA': [], 'Sphaerimonospora_PYR_THR_AA': [], 'Sinosporangium_PYR_THR_AA': [], 'Rhizohabitans_PYR_THR_AA': [], 'Acidimicrobium_PYR_THR_AA': [], 'Actinomadura_PYR_THR_AA': [], 'Thermobispora_R5P_AMINOACIDS': [], 'Herbidospora_R5P_AMINOACIDS': [], 'Microtetraspora_R5P_AMINOACIDS': [], 'Planomonospora_R5P_AMINOACIDS': [], 'Planobispora_R5P_AMINOACIDS': [], 'Thermoactinospora_R5P_AMINOACIDS': [], 'Sphaerisporangium_R5P_AMINOACIDS': [], 'Spongiactinospora_R5P_AMINOACIDS': [], 'Bailinhaonella_R5P_AMINOACIDS': [], 'Thermopolyspora_R5P_AMINOACIDS': [], 'Thermocatellispora_R5P_AMINOACIDS': [], 'Planotetraspora_R5P_AMINOACIDS': [], 'Sphaerimonospora_R5P_AMINOACIDS': [], 'Sinosporangium_R5P_AMINOACIDS': [], 'Rhizohabitans_R5P_AMINOACIDS': [], 'Acidimicrobium_R5P_AMINOACIDS': [], 'Actinomadura_R5P_AMINOACIDS': []}
    return dicc

def holes_los1246_genus():
    holes={'Rhodococcus_3PGA_AMINOACIDS': [['6666666.104474', '6666666.104448', '6666666.104459', '6666666.104441']], 'Nocardia_3PGA_AMINOACIDS': [['6666666.102745', '1133849.7', '6666666.102741', '1127134.4'], ['6666666.102760', '6666666.102764', '6666666.102742', '1127134.4'], ['6666666.102760', '6666666.102797', '6666666.102745', '1127134.4'], ['6666666.102760', '6666666.147084', '6666666.102751', '1127134.4']], 'Amycolatopsis_3PGA_AMINOACIDS': [['749927.13', '6666666.104454', '6666666.104440', '1156913.7'], ['749927.13', '6666666.104506', '6666666.104471', '1156913.7'], ['749927.13', '6666666.104537', '6666666.104468', '1156913.7']], 'Mycobacterium_3PGA_AMINOACIDS': [['6666666.101415', '6666666.110985', '243243.19', '1299323.4'], ['6666666.111038', '6666666.110994', '243243.19', '1299323.4'], ['6666666.101415', '362242.37', '350054.22', '1299323.4'], ['6666666.111038', '6666666.111002', '350054.22', '1299323.4']], 'Arthrobacter_3PGA_AMINOACIDS': [['6666666.104519', '6666666.146813', '6666666.104435', '452863.24']], 'Kocuria_3PGA_AMINOACIDS': [['6666666.104630', '6666666.104640', '6666666.104559', '378753.9']], 'Gordonia_3PGA_AMINOACIDS': [['6666666.110999', '6666666.110976', '6666666.110978', '526226.13'], ['6666666.110998', '6666666.111067', '6666666.110978', '526226.13']], 'Streptomyces_3PGA_AMINOACIDS': [['6666666.112898', '6666666.103212', '6666666.103255', '6666666.102743'], ['6666666.103156', '6666666.103300', '6666666.103153', '6666666.102743'], ['6666666.112898', '6666666.104189', '6666666.103153', '6666666.102743'], ['6666666.103208', '6666666.113101', '6666666.102800', '6666666.102743'], ['6666666.103293', '6666666.103248', '6666666.102800', '6666666.102743'], ['6666666.103255', '6666666.103207', '6666666.103156', '6666666.102743'], ['6666666.104344', '6666666.104668', '6666666.103264', '6666666.102743'], ['6666666.103156', '6666666.113055', '6666666.102800', '6666666.102743'], ['6666666.103293', '6666666.103276', '6666666.103153', '6666666.102743'], ['6666666.104344', '6666666.113129', '6666666.103208', '6666666.102743'], ['6666666.103255', '6666666.103247', '6666666.103155', '6666666.102743'], ['6666666.103293', '6666666.103272', '6666666.103255', '6666666.102743'], ['6666666.103155', '6666666.103278', '6666666.103153', '6666666.102743'], ['6666666.103246', '6666666.103202', '6666666.103153', '6666666.102743'], ['6666666.104344', '6666666.104132', '6666666.103155', '6666666.102743'], ['6666666.103208', '6666666.104202', '6666666.103153', '6666666.102743'], ['6666666.103255', '6666666.103170', '6666666.103208', '6666666.102743'], ['6666666.104344', '6666666.103192', '6666666.103156', '6666666.102743'], ['6666666.103155', '6666666.103275', '6666666.102800', '6666666.102743'], ['6666666.103255', '6666666.113296', '6666666.103246', '6666666.102743'], ['6666666.70540', '6666666.104376', '6666666.103153', '6666666.102743'], ['6666666.103246', '6666666.102809', '6666666.102800', '6666666.102743']], 'Saccharomonospora_3PGA_AMINOACIDS': [['6666666.147079', '6666666.146841', '6666666.104473', '6666666.104462']], 'Brevibacterium_3PGA_AMINOACIDS': [['6666666.104610', '6666666.104648', '6666666.104605', '6666666.104557']], 'Streptacidiphilus_3PGA_AMINOACIDS': [['6666666.113047', '6666666.113043', '6666666.113029', '6666666.112935']], 'Nocardia_ALPHAKETOGLUTARATE_AMINOACIDS': [['6666666.102806', '6666666.102749', '6666666.102802', '6666666.102789'], ['6666666.147084', '6666666.146840', '6666666.102802', '6666666.102798'], ['6666666.102750', '6666666.102797', '6666666.102741', '1133849.7'], ['6666666.147084', '6666666.102798', '6666666.102750', '1133849.7']], 'Amycolatopsis_ALPHAKETOGLUTARATE_AMINOACIDS': [['6666666.104540', '6666666.104454', '6666666.104537', '1156913.7'], ['6666666.104537', '6666666.104511', '6666666.104471', '1156913.7'], ['749927.13', '6666666.104527', '1156913.7', '6666666.104468']], 'Clavibacter_ALPHAKETOGLUTARATE_AMINOACIDS': [['6666666.83178', '6666666.83177', '6666666.111552', '6666666.111522']], 'Mycobacterium_ALPHAKETOGLUTARATE_AMINOACIDS': [['6666666.110965', '6666666.147005', '6666666.111017', '362242.37', '243243.19'], ['6666666.111002', '6666666.111035', '1299323.4', '6666666.110963', '362242.37', '6666666.147047', '6666666.110992'], ['6666666.111029', '6666666.111035', '1299323.4', '6666666.111045', '6666666.110969', '350054.22'], ['6666666.111064', '6666666.111037', '6666666.110969', '350054.22'], ['6666666.111001', '6666666.111066', '6666666.110997', '243243.19'], ['6666666.146836', '6666666.110963', '6666666.111065', '243243.19'], ['6666666.111045', '6666666.110985', '6666666.111017', '362242.37', '6666666.110963', '1299323.4'], ['6666666.146892', '6666666.111002', '6666666.111029', '350054.22'], ['6666666.147047', '6666666.111017', '6666666.110967', '6666666.146888', '6666666.111012', '6666666.110992'], ['6666666.111001', '6666666.110990', '6666666.110965', '243243.19'], ['6666666.111051', '6666666.111009', '6666666.111070', '6666666.111071', '6666666.111045']], 'Nocardioides_ALPHAKETOGLUTARATE_AMINOACIDS': [['6666666.104654', '6666666.112254', '6666666.104556', '6666666.112232', '6666666.104625', '6666666.104570'], ['6666666.112232', '6666666.104617', '6666666.104599', '6666666.104556']], 'Corynebacterium_ALPHAKETOGLUTARATE_AMINOACIDS': [['6666666.147121', '6666666.110956', '6666666.110156', '6666666.110173'], ['6666666.146839', '6666666.110144', '645127.7', '196627.31'], ['6666666.146832', '6666666.110137', '6666666.110178', '6666666.110143'], ['6666666.110155', '6666666.110173', '6666666.110148', '6666666.110137'], ['6666666.147157', '6666666.110143', '6666666.146832', '645127.7'], ['6666666.146827', '6666666.110167', '6666666.110168', '6666666.110140'], ['6666666.110158', '6666666.110154', '6666666.110148', '6666666.110143'], ['6666666.147121', '6666666.110173', '6666666.110155', '6666666.110170', '6666666.110168', '6666666.110153'], ['6666666.146832', '6666666.110170', '6666666.110168', '645127.7'], ['6666666.147110', '6666666.110139', '6666666.110168', '645127.7'], ['6666666.110178', '6666666.110155', '6666666.110170', '6666666.110143'], ['6666666.110956', '6666666.110167', '6666666.110168', '6666666.110153'], ['6666666.146839', '6666666.110153', '6666666.110168', '196627.31'], ['6666666.146879', '6666666.110136', '645127.7', '196627.31']], 'Arthrobacter_ALPHAKETOGLUTARATE_AMINOACIDS': [['452863.24', '6666666.104489', '6666666.104435', '290340.24']], 'Frankia_ALPHAKETOGLUTARATE_AMINOACIDS': [['6666666.147111', '6666666.104732', '6666666.104606', '6666666.104702']], 'Bifidobacterium_ALPHAKETOGLUTARATE_AMINOACIDS': [['6666666.112169', '6666666.112171', '6666666.112167', '367928.21'], ['6666666.112184', '6666666.112172', '6666666.112169', '367928.21'], ['6666666.112169', '6666666.112181', '6666666.112168', '367928.21'], ['6666666.112179', '6666666.112177', '6666666.112169', '367928.21']], 'Gordonia_ALPHAKETOGLUTARATE_AMINOACIDS': [['6666666.111024', '6666666.111003', '6666666.111006', '6666666.110978'], ['6666666.110976', '6666666.111014', '6666666.111024', '6666666.110979'], ['6666666.110976', '6666666.111014', '6666666.110987', '526226.13']], 'Streptomyces_ALPHAKETOGLUTARATE_AMINOACIDS': [['6666666.112944', '6666666.113048', '6666666.104211', '6666666.103155'], ['6666666.70535', '6666666.113069', '6666666.112972', '6666666.103160'], ['6666666.112944', '6666666.113099', '6666666.112852', '6666666.103290'], ['6666666.103248', '6666666.103272', '6666666.103156', '6666666.103183'], ['6666666.113002', '6666666.112835', '6666666.112952', '6666666.112833'], ['6666666.104166', '6666666.103286', '6666666.103284', '6666666.103253'], ['6666666.113258', '6666666.104229', '6666666.113024', '6666666.112944'], ['6666666.103192', '6666666.104375', '6666666.104206', '6666666.104164'], ['6666666.113216', '6666666.104362', '6666666.112880', '6666666.103252'], ['6666666.112912', '6666666.113237', '6666666.113004', '6666666.103305'], ['6666666.113293', '6666666.112962', '6666666.104141', '6666666.102809'], ['6666666.112972', '6666666.113004', '6666666.104376', '6666666.102743'], ['6666666.112972', '6666666.104191', '6666666.112943', '6666666.103160'], ['6666666.104375', '6666666.112989', '6666666.113069', '6666666.104213'], ['6666666.103252', '6666666.103210', '6666666.103250', '6666666.103160'], ['6666666.112972', '6666666.113002', '6666666.103187', '6666666.103160'], ['6666666.103326', '6666666.112901', '6666666.113121', '6666666.147125', '6666666.103250'], ['6666666.103261', '6666666.103153', '6666666.104120', '6666666.112857', '6666666.104213', '6666666.103286', '6666666.102809', '6666666.70535', '6666666.113069', '6666666.113025', '6666666.103171'], ['6666666.119636', '6666666.112884', '6666666.103246', '6666666.104154'], ['6666666.103306', '6666666.103250', '6666666.103210', '6666666.103182'], ['6666666.147138', '6666666.113277', '6666666.113086', '6666666.104154'], ['6666666.104218', '6666666.113046', '6666666.113164', '6666666.103183'], ['6666666.112972', '6666666.112880', '6666666.103252', '6666666.103160'], ['6666666.103192', '6666666.104375', '6666666.112857', '6666666.103261'], ['6666666.113099', '6666666.70540', '6666666.104235', '6666666.104183'], ['6666666.146818', '6666666.113113', '6666666.112880', '6666666.112833'], ['6666666.104349', '6666666.104203', '6666666.103291', '6666666.112961', '6666666.102743'], ['6666666.104314', '6666666.103296', '6666666.103264', '6666666.102743'], ['6666666.113235', '6666666.113148', '6666666.113177', '6666666.112984'], ['6666666.113161', '6666666.113224', '6666666.113093', '6666666.103293'], ['6666666.147054', '6666666.113143', '6666666.104375', '6666666.104213'], ['6666666.112824', '6666666.113293', '6666666.103198', '6666666.103210'], ['6666666.103156', '6666666.103272', '6666666.103286', '6666666.102809'], ['6666666.113226', '6666666.113247', '6666666.113236', '6666666.112972'], ['6666666.113199', '6666666.113129', '6666666.112828', '6666666.103270'], ['6666666.70535', '6666666.104141', '6666666.103264', '6666666.103160'], ['6666666.113002', '6666666.112880', '6666666.104362', '6666666.104177'], ['6666666.112824', '6666666.104222', '6666666.103252', '6666666.103210'], ['6666666.103182', '6666666.119642', '6666666.103291', '6666666.102743'], ['6666666.113273', '6666666.103341', '6666666.103286', '6666666.102809'], ['6666666.104378', '6666666.104323', '6666666.104349', '6666666.102743'], ['6666666.103249', '6666666.113099', '6666666.112944', '6666666.103155'], ['6666666.113268', '6666666.113122', '6666666.104213', '6666666.103153'], ['6666666.103326', '6666666.112901', '6666666.103291', '6666666.102743'], ['6666666.112972', '6666666.113069', '6666666.112952', '6666666.112833'], ['6666666.113024', '6666666.104331', '6666666.113181', '6666666.103155'], ['6666666.113273', '6666666.70535', '6666666.113069', '6666666.112876'], ['6666666.103192', '6666666.113037', '6666666.104135', '6666666.104177', '6666666.104120'], ['6666666.112880', '6666666.113113', '6666666.103293', '6666666.103156', '6666666.104222', '6666666.103252'], ['6666666.113161', '6666666.104213', '6666666.104375', '6666666.104206'], ['6666666.113216', '6666666.113227', '6666666.102800', '6666666.104360'], ['6666666.113236', '6666666.113247', '6666666.103163', '6666666.104132'], ['6666666.112956', '6666666.104321', '6666666.103192', '6666666.104164'], ['6666666.113293', '6666666.102809', '6666666.112943', '6666666.103198'], ['6666666.104120', '6666666.113067', '6666666.103261', '6666666.103153'], ['6666666.112918', '6666666.70537', '6666666.113037', '6666666.104120'], ['6666666.147044', '6666666.119638', '6666666.146880', '6666666.104197'], ['6666666.104340', '6666666.113174', '6666666.103163', '6666666.104206'], ['6666666.113206', '6666666.112930', '6666666.103250', '6666666.103160'], ['6666666.112944', '6666666.113024', '6666666.113025', '6666666.103155'], ['6666666.113206', '6666666.113211', '6666666.112992', '6666666.103187'], ['6666666.113293', '6666666.147135', '6666666.113238', '6666666.102809'], ['6666666.103290', '6666666.112852', '6666666.113296', '6666666.103156', '6666666.103183'], ['6666666.103330', '6666666.146847', '6666666.103257', '6666666.103201'], ['6666666.146880', '6666666.119638', '6666666.113025', '6666666.112876'], ['6666666.104314', '6666666.103333', '6666666.103187', '6666666.102743'], ['6666666.113177', '6666666.113002', '6666666.112876', '6666666.112828', '6666666.104202', '6666666.104183'], ['6666666.113164', '6666666.103317', '6666666.112835', '6666666.103187'], ['6666666.113162', '6666666.113237', '6666666.112912', '6666666.147135'], ['6666666.112898', '6666666.113132', '6666666.112880', '6666666.103252'], ['6666666.113238', '6666666.113055', '6666666.112943', '6666666.102809'], ['6666666.146823', '6666666.113185', '6666666.104349', '6666666.103290'], ['6666666.104183', '6666666.113099', '6666666.112944', '6666666.103155'], ['6666666.146849', '6666666.112962', '6666666.113293', '6666666.103284'], ['6666666.70535', '6666666.104222', '6666666.103252', '6666666.103160'], ['6666666.113155', '6666666.146847', '6666666.103330', '6666666.104183'], ['6666666.113206', '6666666.112833', '6666666.112972', '6666666.103160'], ['6666666.113264', '6666666.113226', '6666666.112972', '6666666.103155'], ['6666666.113164', '6666666.103257', '6666666.103201', '6666666.102743'], ['6666666.103296', '6666666.103326', '6666666.147159', '6666666.103250', '6666666.102743'], ['6666666.104370', '6666666.112835', '6666666.113069', '6666666.113025', '6666666.103171', '6666666.103155'], ['6666666.112943', '6666666.104131', '6666666.103183', '6666666.103187'], ['6666666.103250', '6666666.104720', '6666666.103187', '6666666.103160'], ['6666666.112943', '6666666.103198', '6666666.103250', '6666666.103160'], ['6666666.103300', '6666666.104668', '6666666.113037', '6666666.104120'], ['6666666.113293', '6666666.102809', '6666666.103286', '6666666.103284'], ['6666666.104320', '6666666.112989', '6666666.113069', '6666666.102809'], ['6666666.104183', '6666666.113155', '6666666.113237', '6666666.103155'], ['6666666.113177', '6666666.104183', '6666666.103212', '6666666.103330', '6666666.104314', '6666666.103333'], ['6666666.112992', '6666666.113282', '6666666.112943', '6666666.103187'], ['6666666.147159', '6666666.113212', '6666666.113025', '6666666.103250'], ['6666666.103208', '6666666.104154', '6666666.119636', '6666666.113093', '6666666.103293'], ['6666666.112943', '6666666.113055', '6666666.112898', '6666666.103252'], ['6666666.119636', '6666666.104189', '6666666.103215', '6666666.104154'], ['6666666.70535', '6666666.102809', '6666666.112943', '6666666.103160'], ['6666666.112898', '6666666.103252', '6666666.103210', '6666666.103182'], ['6666666.103293', '6666666.103208', '6666666.104362', '6666666.103252'], ['6666666.104370', '6666666.112888', '6666666.104211', '6666666.103155'], ['6666666.103258', '6666666.104668', '6666666.103255', '6666666.113185', '6666666.103286'], ['6666666.113002', '6666666.112880', '6666666.112888', '6666666.104370'], ['6666666.104131', '6666666.113187', '6666666.112839', '6666666.112924', '6666666.147163', '6666666.112858', '6666666.103247', '6666666.103287', '6666666.103210', '6666666.102743'], ['6666666.113076', '6666666.113238', '6666666.112943', '6666666.103183'], ['6666666.103257', '6666666.104348', '6666666.104320', '6666666.103201'], ['6666666.146818', '6666666.113236', '6666666.112972', '6666666.112833'], ['6666666.112921', '6666666.113187', '6666666.113164', '6666666.113076'], ['6666666.112824', '6666666.112853', '6666666.104720', '6666666.103210'], ['6666666.113080', '6666666.113211', '6666666.112943', '6666666.104360']], 'Microbacterium_ALPHAKETOGLUTARATE_AMINOACIDS': [['6666666.147106', '6666666.104495', '6666666.104449', '6666666.104436'], ['6666666.147106', '6666666.147166', '6666666.104450', '6666666.104436']], 'Saccharomonospora_ALPHAKETOGLUTARATE_AMINOACIDS': [['6666666.147170', '6666666.104473', '6666666.147079', '6666666.104536']], 'Pseudonocardia_ALPHAKETOGLUTARATE_AMINOACIDS': [['6666666.104611', '6666666.146903', '6666666.112804', '6666666.104623', '6666666.104550'], ['6666666.104596', '6666666.104554', '6666666.104611', '6666666.104550']], 'Knoellia_ALPHAKETOGLUTARATE_AMINOACIDS': [['6666666.111590', '6666666.111537', '6666666.111534', '6666666.111524']], 'Actinomadura_ALPHAKETOGLUTARATE_AMINOACIDS': [['6666666.112788', '6666666.112796', '6666666.112800', '6666666.112791']], 'Actinomyces_ALPHAKETOGLUTARATE_AMINOACIDS': [['6666666.114248', '6666666.114230', '6666666.114208', '6666666.114206'], ['6666666.114252', '6666666.114197', '6666666.114195', '6666666.114192'], ['6666666.114233', '6666666.114193', '6666666.114214', '6666666.114206'], ['6666666.114248', '6666666.114206', '6666666.114196', '6666666.114195']], 'Rhodococcus_E4P_AMINO_ACIDS': [['6666666.104452', '6666666.104466', '632772.10', '101510.30']], 'Nocardia_E4P_AMINO_ACIDS': [['6666666.102747', '6666666.102741', '1133849.7', '1127134.4'], ['6666666.102779', '6666666.102739', '1133849.7', '1127134.4'], ['6666666.102748', '6666666.102789', '6666666.102742', '6666666.102741'], ['6666666.102748', '6666666.102806', '6666666.102747', '6666666.102741'], ['6666666.102785', '6666666.102760', '1133849.7', '1127134.4']], 'Amycolatopsis_E4P_AMINO_ACIDS': [['6666666.112803', '749927.13', '6666666.104442', '6666666.104475'], ['6666666.112803', '1156913.7', '6666666.104510', '6666666.104454'], ['6666666.104510', '6666666.104511', '6666666.104471', '6666666.104454'], ['6666666.112803', '6666666.104475', '6666666.104537', '6666666.104454']], 'Mycobacterium_E4P_AMINO_ACIDS': [['6666666.111035', '6666666.111047', '6666666.110969', '1299323.4'], ['6666666.111004', '6666666.111047', '6666666.110969', '6666666.110963'], ['6666666.111012', '6666666.111004', '6666666.110963', '243243.19'], ['6666666.110969', '6666666.110963', '243243.19', '1299323.4'], ['6666666.111035', '6666666.111012', '243243.19', '1299323.4'], ['6666666.111070', '350058.13', '243243.19', '1299323.4'], ['6666666.111015', '246196.47', '6666666.110963', '243243.19'], ['6666666.111015', '6666666.111023', '350058.13', '243243.19'], ['6666666.110969', '6666666.111010', '246196.47', '6666666.110963']], 'Corynebacterium_E4P_AMINO_ACIDS': [['6666666.110148', '645127.7', '504474.6', '196627.31'], ['6666666.110185', '6666666.110181', '6666666.110140', '196627.31'], ['6666666.146832', '6666666.110182', '6666666.110142', '6666666.110148'], ['6666666.110185', '6666666.110183', '504474.6', '196627.31'], ['6666666.110140', '6666666.110167', '504474.6', '196627.31'], ['6666666.110148', '6666666.146832', '6666666.110140', '196627.31'], ['6666666.110185', '6666666.110142', '6666666.110148', '196627.31'], ['6666666.110141', '6666666.110134', '6666666.110140', '196627.31'], ['6666666.146832', '6666666.110182', '6666666.110181', '6666666.110140']], 'Frankia_E4P_AMINO_ACIDS': [['6666666.104703', '6666666.147111', '6666666.104732', '298653.10']], 'Salinispora_E4P_AMINO_ACIDS': [['6666666.104642', '391037.10', '6666666.104582', '369723.10']], 'Gordonia_E4P_AMINO_ACIDS': [['6666666.111007', '6666666.110978', '6666666.110987', '526226.13'], ['6666666.110998', '6666666.111053', '6666666.110979', '526226.13']], 'Streptomyces_E4P_AMINO_ACIDS': [['6666666.113185', '6666666.112944', '6666666.104715', '6666666.103276'], ['6666666.104188', '6666666.113300', '6666666.103271', '6666666.103156'], ['6666666.104154', '6666666.112857', '6666666.103276', '6666666.103153'], ['6666666.104189', '6666666.103280', '6666666.103198', '6666666.102743'], ['6666666.104188', '6666666.104323', '6666666.103298', '6666666.103156'], ['6666666.103156', '6666666.112835', '6666666.102800', '6666666.103183'], ['6666666.113007', '6666666.112921', '6666666.103171', '6666666.103183'], ['6666666.103276', '6666666.112857', '6666666.103171', '6666666.103183'], ['6666666.103306', '6666666.146821', '6666666.103198', '6666666.102743'], ['6666666.113114', '6666666.113007', '6666666.103183', '6666666.103160'], ['6666666.103272', '6666666.112961', '6666666.103183', '6666666.102743'], ['6666666.113101', '6666666.112918', '6666666.104325', '6666666.103160'], ['6666666.113236', '6666666.70532', '6666666.104355', '6666666.103153'], ['6666666.104189', '6666666.104338', '6666666.103160', '6666666.102743'], ['6666666.103171', '6666666.104188', '6666666.103156', '6666666.103183'], ['6666666.103281', '6666666.104141', '6666666.113174', '6666666.103178'], ['6666666.104355', '6666666.112973', '6666666.104325', '6666666.103160'], ['6666666.103198', '6666666.104715', '6666666.103153', '6666666.102743'], ['6666666.103293', '6666666.103278', '6666666.103198', '6666666.102743'], ['6666666.146880', '6666666.112849', '6666666.104336', '6666666.103171'], ['6666666.113236', '6666666.104132', '6666666.104715', '6666666.103153'], ['6666666.103171', '6666666.103310', '6666666.102800', '6666666.103183'], ['6666666.112857', '6666666.113133', '6666666.103310', '6666666.103171'], ['6666666.103183', '6666666.103276', '6666666.103153', '6666666.102743'], ['6666666.104164', '6666666.103247', '6666666.103276', '6666666.103183'], ['6666666.104189', '6666666.146816', '6666666.103306', '6666666.102743'], ['6666666.103298', '6666666.103156', '6666666.103183', '6666666.102743'], ['6666666.103276', '6666666.113247', '6666666.102800', '6666666.103183'], ['6666666.113041', '6666666.103255', '6666666.103198', '6666666.103183'], ['6666666.103320', '6666666.147163', '6666666.103264', '6666666.103201'], ['6666666.112912', '6666666.102800', '6666666.103183', '6666666.102743'], ['6666666.112912', '6666666.113101', '6666666.103160', '6666666.102743'], ['6666666.103298', '6666666.104323', '6666666.103293', '6666666.102743'], ['6666666.113247', '6666666.112894', '6666666.113185', '6666666.103276'], ['6666666.113236', '6666666.103258', '6666666.103306', '6666666.102743'], ['6666666.104355', '6666666.112973', '6666666.104154', '6666666.103153'], ['6666666.104362', '6666666.112864', '6666666.103183', '6666666.102743'], ['6666666.112835', '6666666.103215', '6666666.104188', '6666666.103156'], ['6666666.104189', '6666666.104336', '6666666.103183', '6666666.102743'], ['6666666.113236', '6666666.103247', '6666666.103276', '6666666.103153'], ['6666666.103293', '6666666.104154', '6666666.103153', '6666666.102743'], ['6666666.112912', '6666666.104370', '6666666.103198', '6666666.102743'], ['6666666.103298', '6666666.113249', '6666666.103160', '6666666.102743'], ['6666666.103293', '6666666.103171', '6666666.103183', '6666666.102743'], ['6666666.104336', '6666666.113185', '6666666.103276', '6666666.103183'], ['6666666.103160', '6666666.104355', '6666666.103153', '6666666.102743'], ['6666666.103306', '6666666.147135', '6666666.103160', '6666666.102743'], ['6666666.104132', '6666666.103264', '6666666.103281', '6666666.103247'], ['6666666.113247', '6666666.113133', '6666666.103310', '6666666.102800'], ['6666666.104188', '6666666.113300', '6666666.103278', '6666666.103171'], ['6666666.104715', '6666666.112944', '6666666.103280', '6666666.103198'], ['6666666.112912', '6666666.70540', '6666666.104189', '6666666.102743'], ['6666666.103178', '6666666.104120', '6666666.103156', '6666666.103183'], ['6666666.103293', '6666666.104325', '6666666.103160', '6666666.102743'], ['6666666.112857', '6666666.103276', '6666666.103247', '6666666.103192'], ['6666666.70532', '6666666.147163', '6666666.103281', '6666666.103247'], ['6666666.103306', '6666666.103178', '6666666.103183', '6666666.102743'], ['6666666.103208', '6666666.104376', '6666666.103198', '6666666.103183'], ['6666666.113101', '6666666.112918', '6666666.103310', '6666666.102800'], ['6666666.112835', '6666666.103215', '6666666.103310', '6666666.102800'], ['6666666.103298', '6666666.103271', '6666666.103198', '6666666.102743']], 'Microbacterium_E4P_AMINO_ACIDS': [['6666666.104451', '6666666.104450', '6666666.104449', '6666666.104436']], 'Actinomyces_E4P_AMINO_ACIDS': [['6666666.114248', '6666666.114230', '6666666.114212', '6666666.114195'], ['6666666.114198', '6666666.114219', '6666666.114195', '6666666.114192'], ['6666666.114199', '6666666.114208', '6666666.114234', '6666666.114197', '6666666.114209'], ['6666666.114219', '6666666.114191', '6666666.114209', '6666666.114198'], ['6666666.114248', '6666666.114191', '6666666.114219', '6666666.114195'], ['6666666.114217', '6666666.114199', '6666666.114209', '6666666.114198']], 'Rhodococcus_Glycolysis': [['6666666.104474', '6666666.147154', '6666666.104441', '632772.10'], ['6666666.104492', '6666666.104448', '6666666.104474', '632772.10']], 'Amycolatopsis_Glycolysis': [['6666666.104468', '6666666.104475', '6666666.104442', '1156913.7']], 'Mycobacterium_Glycolysis': [['6666666.110986', '6666666.111056', '362242.37', '246196.47'], ['6666666.111066', '6666666.110967', '362242.37', '1299323.4'], ['6666666.111070', '6666666.111066', '6666666.110967', '350054.22'], ['6666666.110986', '6666666.111032', '350054.22', '246196.47'], ['6666666.110992', '6666666.111070', '350054.22', '246196.47'], ['362242.37', '6666666.110967', '350054.22', '246196.47'], ['6666666.110992', '1299323.4', '362242.37', '246196.47']], 'Corynebacterium_Glycolysis': [['6666666.110167', '6666666.110139', '504474.6', '196627.31']], 'Gordonia_Glycolysis': [['6666666.111018', '6666666.111067', '6666666.110979', '6666666.110976'], ['6666666.146910', '6666666.110998', '6666666.110976', '526226.13']], 'Streptomyces_Glycolysis': [['6666666.112828', '6666666.112821', '6666666.103271', '6666666.103156'], ['6666666.113037', '6666666.112864', '6666666.103198', '6666666.102743'], ['6666666.113037', '6666666.103271', '6666666.103156', '6666666.102743'], ['6666666.112820', '6666666.112828', '6666666.103156', '6666666.102743']], 'Rhodococcus_OXALACETATE_AMINOACIDS': [['6666666.104493', '6666666.104481', '6666666.104459', '6666666.104441'], ['6666666.146894', '6666666.104474', '6666666.104492', '6666666.104441'], ['6666666.146894', '6666666.104441', '6666666.104493', '6666666.104452']], 'Nocardia_OXALACETATE_AMINOACIDS': [['6666666.102749', '6666666.102746', '6666666.102747', '6666666.102742'], ['6666666.102789', '6666666.102746', '6666666.102747', '6666666.102739'], ['6666666.102742', '6666666.102747', '6666666.102739', '1127134.4'], ['6666666.102753', '6666666.102749', '6666666.102742', '1127134.4'], ['6666666.102753', '6666666.102789', '6666666.102739', '1127134.4']], 'Amycolatopsis_OXALACETATE_AMINOACIDS': [['6666666.104540', '6666666.104541', '6666666.104445', '1156913.7'], ['6666666.112807', '6666666.104454', '6666666.104541', '6666666.104540'], ['6666666.104540', '6666666.104511', '6666666.104471', '1156913.7']], 'Nocardiopsis_OXALACETATE_AMINOACIDS': [['6666666.105885', '6666666.105881', '6666666.105879', '1205910.5']], 'Mycobacterium_OXALACETATE_AMINOACIDS': [['6666666.111064', '350058.13', '362242.37', '1299323.4'], ['6666666.110968', '6666666.110992', '243243.19', '1299323.4'], ['6666666.111004', '6666666.110991', '350058.13', '350054.22'], ['6666666.110968', '350054.22', '362242.37', '1299323.4'], ['6666666.110968', '6666666.147047', '6666666.110967', '1299323.4']], 'Nocardioides_OXALACETATE_AMINOACIDS': [['6666666.104654', '6666666.112232', '6666666.104625', '6666666.104599'], ['6666666.104599', '6666666.104654', '6666666.104556', '196162.16']], 'Corynebacterium_OXALACETATE_AMINOACIDS': [['6666666.110148', '6666666.110139', '6666666.110134', '196627.31'], ['6666666.110185', '6666666.110136', '6666666.110134', '196627.31'], ['6666666.110161', '504474.6', '6666666.110134', '196627.31']], 'Arthrobacter_OXALACETATE_AMINOACIDS': [['6666666.147082', '6666666.104435', '6666666.104477', '290340.24']], 'Bifidobacterium_OXALACETATE_AMINOACIDS': [['6666666.112175', '6666666.112192', '6666666.112172', '442563.11'], ['6666666.112172', '6666666.112171', '6666666.112167', '442563.11'], ['6666666.112181', '6666666.112634', '6666666.112172', '442563.11'], ['6666666.112187', '6666666.112198', '6666666.112172', '442563.11']], 'Gordonia_OXALACETATE_AMINOACIDS': [['526226.13', '6666666.110976', '6666666.111003', '6666666.110978']], 'Micromonospora_OXALACETATE_AMINOACIDS': [['6666666.104635', '6666666.104615', '6666666.104628', '644283.11']], 'Streptomyces_OXALACETATE_AMINOACIDS': [['6666666.103335', '6666666.112922', '6666666.103293', '6666666.103199'], ['6666666.113148', '6666666.113069', '6666666.104174', '6666666.103252'], ['6666666.104188', '6666666.104174', '6666666.103271', '6666666.102743'], ['6666666.113099', '6666666.103171', '6666666.103305', '6666666.103204'], ['6666666.113140', '6666666.113061', '6666666.103261', '6666666.103153'], ['6666666.113080', '6666666.104715', '6666666.113061', '6666666.104203'], ['6666666.112924', '6666666.113048', '6666666.103317', '6666666.103297'], ['6666666.112847', '6666666.104331', '6666666.103170', '6666666.103156'], ['6666666.103290', '6666666.147119', '6666666.103212', '6666666.103156'], ['6666666.113283', '6666666.104325', '6666666.103261', '6666666.103155'], ['6666666.112992', '6666666.147090', '6666666.103178', '6666666.103156'], ['6666666.113236', '6666666.103253', '6666666.103204', '6666666.102743'], ['6666666.112944', '6666666.112924', '6666666.103297', '6666666.103204'], ['6666666.104378', '6666666.103272', '6666666.103199', '6666666.102800'], ['6666666.112852', '6666666.112847', '6666666.112820', '6666666.103160'], ['6666666.147090', '6666666.113055', '6666666.104166', '6666666.103178'], ['6666666.103246', '6666666.103257', '6666666.103212', '6666666.103156'], ['6666666.113210', '6666666.113114', '6666666.104136', '6666666.103250'], ['6666666.103276', '6666666.104325', '6666666.103208', '6666666.103170'], ['6666666.103247', '6666666.112931', '6666666.103212', '6666666.103156'], ['6666666.104141', '6666666.104355', '6666666.104123', '6666666.103261'], ['6666666.104715', '6666666.113061', '6666666.103204', '6666666.102743'], ['6666666.112820', '6666666.104323', '6666666.112819', '6666666.104166', '6666666.103178', '6666666.103156'], ['6666666.104174', '6666666.103252', '6666666.103290', '6666666.103156'], ['6666666.104190', '6666666.113054', '6666666.103247', '6666666.103156'], ['6666666.104174', '6666666.103271', '6666666.103199', '6666666.103156'], ['6666666.112918', '6666666.103183', '6666666.103276', '6666666.103170'], ['6666666.103300', '6666666.104141', '6666666.103198', '6666666.103170'], ['6666666.104378', '6666666.113114', '6666666.104177', '6666666.102800'], ['6666666.104188', '6666666.112915', '6666666.103250', '6666666.102743'], ['6666666.112820', '6666666.112864', '6666666.103246', '6666666.103156'], ['6666666.146880', '6666666.113207', '6666666.113061', '6666666.104708'], ['6666666.112962', '6666666.147124', '6666666.103276', '6666666.103170'], ['6666666.146818', '6666666.103272', '6666666.103271', '6666666.102743'], ['6666666.146818', '6666666.113080', '6666666.104715', '6666666.102743'], ['6666666.104182', '6666666.112924', '6666666.103210', '6666666.103178'], ['6666666.103271', '6666666.103199', '6666666.103160', '6666666.102743'], ['6666666.103248', '6666666.103182', '6666666.103204', '6666666.102743'], ['6666666.104188', '6666666.112820', '6666666.103160', '6666666.102743'], ['6666666.104375', '6666666.70536', '6666666.103276', '6666666.103170'], ['6666666.104375', '6666666.113235', '6666666.103210', '6666666.103178'], ['6666666.70542', '6666666.113210', '6666666.146818', '6666666.103163'], ['6666666.104177', '6666666.147163', '6666666.103199', '6666666.102800'], ['6666666.113061', '6666666.104708', '6666666.104123', '6666666.103261'], ['6666666.113236', '6666666.112949', '6666666.103271', '6666666.102743'], ['6666666.104715', '6666666.113207', '6666666.103250', '6666666.102743'], ['6666666.112820', '6666666.112962', '6666666.103170', '6666666.103156'], ['6666666.104123', '6666666.112956', '6666666.103305', '6666666.103261'], ['6666666.113069', '6666666.70536', '6666666.104375', '6666666.103178'], ['6666666.104320', '6666666.112852', '6666666.103199', '6666666.102800'], ['6666666.104375', '6666666.104336', '6666666.103300', '6666666.103170'], ['6666666.112944', '6666666.112956', '6666666.103305', '6666666.103204'], ['6666666.113268', '6666666.103252', '6666666.103155', '6666666.103153'], ['6666666.113140', '6666666.103204', '6666666.103305', '6666666.103153'], ['6666666.112992', '6666666.104191', '6666666.104174', '6666666.103156'], ['6666666.104174', '6666666.113069', '6666666.103178', '6666666.103156'], ['6666666.103246', '6666666.103335', '6666666.103199', '6666666.103156'], ['6666666.103204', '6666666.103297', '6666666.103160', '6666666.102743'], ['6666666.103297', '6666666.113054', '6666666.103182', '6666666.103204'], ['6666666.147124', '6666666.103330', '6666666.104206', '6666666.103261'], ['6666666.113099', '6666666.113214', '6666666.112944', '6666666.103204'], ['6666666.103246', '6666666.103198', '6666666.103170', '6666666.103156'], ['6666666.112972', '6666666.104205', '6666666.104348', '6666666.103155'], ['6666666.119636', '6666666.103215', '6666666.112858', '6666666.103257'], ['6666666.103290', '6666666.103262', '6666666.103246', '6666666.103156'], ['6666666.113236', '6666666.104136', '6666666.103250', '6666666.102743'], ['6666666.103199', '6666666.112918', '6666666.103170', '6666666.103156'], ['6666666.112820', '6666666.104190', '6666666.103297', '6666666.103160'], ['6666666.112820', '6666666.103156', '6666666.103199', '6666666.103160'], ['6666666.113268', '6666666.113143', '6666666.113148', '6666666.103252'], ['6666666.113061', '6666666.104708', '6666666.112944', '6666666.103204'], ['6666666.104375', '6666666.104336', '6666666.104182', '6666666.103178'], ['6666666.113046', '6666666.103335', '6666666.103246', '6666666.103198'], ['6666666.113283', '6666666.104136', '6666666.103253', '6666666.103155'], ['6666666.112820', '6666666.104205', '6666666.103290', '6666666.103156'], ['6666666.113207', '6666666.104325', '6666666.112915', '6666666.103250'], ['6666666.112943', '6666666.104135', '6666666.104321', '6666666.104205'], ['6666666.112949', '6666666.104136', '6666666.104177', '6666666.102800'], ['6666666.113268', '6666666.104174', '6666666.103305', '6666666.103153'], ['6666666.103271', '6666666.112901', '6666666.103210', '6666666.103199'], ['6666666.112849', '6666666.113215', '6666666.103170', '6666666.103156'], ['6666666.112918', '6666666.113046', '6666666.103198', '6666666.103170'], ['6666666.112930', '6666666.113048', '6666666.103317', '6666666.103253'], ['6666666.104715', '6666666.103183', '6666666.103271', '6666666.102743'], ['6666666.112949', '6666666.119647', '6666666.103293', '6666666.102800'], ['6666666.113041', '6666666.113206', '6666666.112992', '6666666.103170'], ['6666666.112900', '6666666.103293', '6666666.103199', '6666666.103156'], ['6666666.113283', '6666666.112915', '6666666.103305', '6666666.103155'], ['6666666.112847', '6666666.104183', '6666666.104174', '6666666.103156'], ['6666666.112944', '6666666.112930', '6666666.103253', '6666666.103204'], ['6666666.112918', '6666666.103326', '6666666.103317', '6666666.102800'], ['6666666.104348', '6666666.147124', '6666666.103261', '6666666.103155'], ['6666666.113283', '6666666.147119', '6666666.112972', '6666666.103155'], ['6666666.112918', '6666666.113235', '6666666.104375', '6666666.103170'], ['6666666.113069', '6666666.112901', '6666666.103210', '6666666.103178'], ['6666666.103290', '6666666.102800', '6666666.103199', '6666666.103156'], ['6666666.113041', '6666666.103287', '6666666.103290', '6666666.103170'], ['6666666.112972', '6666666.103317', '6666666.103253', '6666666.103155'], ['6666666.104348', '6666666.104321', '6666666.103280', '6666666.103155'], ['6666666.104174', '6666666.103276', '6666666.103170', '6666666.103156'], ['6666666.104190', '6666666.103300', '6666666.103170', '6666666.103156'], ['6666666.104141', '6666666.112858', '6666666.103246', '6666666.103198'], ['6666666.112972', '6666666.103290', '6666666.103252', '6666666.103155'], ['6666666.104348', '6666666.113236', '6666666.103253', '6666666.103155'], ['6666666.113268', '6666666.103276', '6666666.103261', '6666666.103153'], ['6666666.146818', '6666666.113099', '6666666.113214', '6666666.103163'], ['6666666.113069', '6666666.112956', '6666666.104182', '6666666.103178'], ['6666666.104123', '6666666.70536', '6666666.103276', '6666666.103261'], ['6666666.113099', '6666666.104203', '6666666.113061', '6666666.103204'], ['6666666.147135', '6666666.113238', '6666666.113099', '6666666.103297'], ['6666666.104321', '6666666.104323', '6666666.103171', '6666666.103280'], ['6666666.146818', '6666666.103163', '6666666.112819', '6666666.104323'], ['6666666.112847', '6666666.113288', '6666666.103246', '6666666.103156'], ['6666666.104174', '6666666.112915', '6666666.103212', '6666666.103156'], ['6666666.103253', '6666666.112949', '6666666.103252', '6666666.103155'], ['6666666.104323', '6666666.104188', '6666666.103305', '6666666.103171'], ['6666666.113061', '6666666.103326', '6666666.103297', '6666666.103204'], ['6666666.103178', '6666666.104375', '6666666.103170', '6666666.103156'], ['6666666.112900', '6666666.113268', '6666666.104174', '6666666.103156'], ['6666666.112918', '6666666.103326', '6666666.103300', '6666666.103170'], ['6666666.104378', '6666666.112898', '6666666.103317', '6666666.102800'], ['6666666.103305', '6666666.103171', '6666666.103280', '6666666.103155'], ['6666666.103261', '6666666.104206', '6666666.103192', '6666666.103153'], ['6666666.104235', '6666666.104202', '6666666.104174', '6666666.103156'], ['6666666.112972', '6666666.104190', '6666666.103305', '6666666.103155'], ['6666666.112944', '6666666.112901', '6666666.103271', '6666666.103204'], ['6666666.146818', '6666666.113210', '6666666.103250', '6666666.102743'], ['6666666.113268', '6666666.119647', '6666666.113140', '6666666.103153'], ['6666666.147090', '6666666.113055', '6666666.113206', '6666666.112992'], ['6666666.146818', '6666666.104323', '6666666.104188', '6666666.102743'], ['6666666.103280', '6666666.103287', '6666666.103252', '6666666.103155'], ['6666666.103317', '6666666.103297', '6666666.103199', '6666666.102800'], ['6666666.103280', '6666666.112898', '6666666.103253', '6666666.103155'], ['6666666.104190', '6666666.104182', '6666666.103178', '6666666.103156'], ['6666666.112949', '6666666.103183', '6666666.112918', '6666666.102800'], ['6666666.103250', '6666666.147163', '6666666.103160', '6666666.102743'], ['6666666.104123', '6666666.104336', '6666666.103300', '6666666.103261'], ['6666666.112956', '6666666.104182', '6666666.104190', '6666666.103305'], ['6666666.104715', '6666666.147124', '6666666.104188', '6666666.102743'], ['6666666.112992', '6666666.112943', '6666666.112820', '6666666.103156'], ['6666666.112972', '6666666.103300', '6666666.103261', '6666666.103155'], ['6666666.113099', '6666666.112898', '6666666.103253', '6666666.103204'], ['6666666.103199', '6666666.103210', '6666666.103178', '6666666.103156'], ['6666666.113206', '6666666.104135', '6666666.104321', '6666666.103287'], ['6666666.104190', '6666666.103305', '6666666.104174', '6666666.103156'], ['6666666.104378', '6666666.103287', '6666666.103290', '6666666.102800'], ['6666666.146818', '6666666.113099', '6666666.103204', '6666666.102743'], ['6666666.113140', '6666666.119647', '6666666.103271', '6666666.103204'], ['6666666.104188', '6666666.103305', '6666666.103204', '6666666.102743'], ['6666666.113124', '6666666.103264', '6666666.103257', '6666666.103246'], ['6666666.103317', '6666666.112972', '6666666.103290', '6666666.102800'], ['6666666.103305', '6666666.104132', '6666666.103192', '6666666.103153'], ['6666666.104177', '6666666.147119', '6666666.103290', '6666666.102800'], ['6666666.147163', '6666666.103212', '6666666.112820', '6666666.103160'], ['6666666.104348', '6666666.104188', '6666666.103305', '6666666.103155'], ['6666666.103297', '6666666.112924', '6666666.103210', '6666666.103199'], ['6666666.103212', '6666666.147092', '6666666.103178', '6666666.103156'], ['6666666.113140', '6666666.103253', '6666666.103155', '6666666.103153'], ['6666666.104190', '6666666.112858', '6666666.103246', '6666666.103156'], ['6666666.112949', '6666666.103252', '6666666.103290', '6666666.102800'], ['6666666.103212', '6666666.103208', '6666666.103170', '6666666.103156']], 'Microbacterium_OXALACETATE_AMINOACIDS': [['6666666.104505', '6666666.104458', '6666666.104496', '6666666.104436'], ['6666666.104521', '6666666.147166', '6666666.104496', '6666666.104436']], 'Brevibacterium_OXALACETATE_AMINOACIDS': [['6666666.104648', '6666666.104605', '6666666.104598', '6666666.104557']], 'Actinomyces_OXALACETATE_AMINOACIDS': [['6666666.114197', '6666666.114199', '6666666.114193', '6666666.114195']], 'Rhodococcus_PYR_THR_AA': [['6666666.104501', '6666666.104534', '6666666.104493', '6666666.104452']], 'Nocardia_PYR_THR_AA': [['6666666.146840', '6666666.102785', '6666666.102780', '6666666.102742'], ['6666666.146840', '6666666.102760', '6666666.102741', '6666666.102742'], ['6666666.102780', '6666666.102753', '6666666.102749', '6666666.102742'], ['6666666.102780', '6666666.102742', '6666666.102739', '1127134.4'], ['6666666.102780', '6666666.102751', '6666666.102750', '6666666.102741']], 'Amycolatopsis_PYR_THR_AA': [['6666666.104468', '1156913.7', '6666666.104442', '6666666.104440'], ['6666666.104475', '6666666.104540', '6666666.104468', '1156913.7']], 'Nocardiopsis_PYR_THR_AA': [['446468.16', '6666666.105881', '6666666.105879', '6666666.104377', '1205910.5'], ['6666666.105886', '6666666.105889', '6666666.104377', '1205910.5']], 'Mycobacterium_PYR_THR_AA': [['6666666.111043', '6666666.111004', '6666666.111037', '6666666.110992'], ['6666666.111022', '6666666.111004', '6666666.111037', '6666666.110969', '6666666.110983', '1299323.4'], ['6666666.111001', '243243.19', '6666666.110980', '6666666.110965'], ['6666666.146873', '6666666.111043', '6666666.110992', '6666666.101415'], ['6666666.110983', '362242.37', '6666666.101415', '350054.22'], ['6666666.146873', '6666666.147047', '6666666.111001', '6666666.101415'], ['6666666.101415', '6666666.110991', '350058.13', '350054.22'], ['6666666.111064', '6666666.110977', '362242.37', '6666666.101415'], ['6666666.110972', '6666666.110967', '6666666.110968', '350054.22'], ['6666666.110983', '1299323.4', '6666666.110972', '350054.22'], ['6666666.111045', '6666666.110977', '6666666.110984', '6666666.147127', '6666666.111015', '6666666.110993', '6666666.110986', '6666666.110980'], ['6666666.110993', '6666666.110986', '6666666.110992', '6666666.110969'], ['6666666.110992', '6666666.110986', '6666666.110980', '6666666.110965'], ['6666666.110969', '6666666.110992', '6666666.101415', '350054.22'], ['6666666.147127', '6666666.111015', '6666666.110993', '6666666.110969', '6666666.110983', '1299323.4'], ['6666666.110972', '6666666.146873', '6666666.101415', '350054.22']], 'Nocardioides_PYR_THR_AA': [['6666666.104570', '6666666.104631', '6666666.104556', '6666666.104625']], 'Corynebacterium_PYR_THR_AA': [['6666666.146826', '6666666.110136', '6666666.110134', '196627.31'], ['6666666.110165', '504474.6', '6666666.110156', '6666666.110134'], ['6666666.146826', '6666666.110144', '6666666.110147', '196627.31'], ['6666666.146827', '6666666.110165', '6666666.110134', '196627.31'], ['6666666.110147', '6666666.110156', '6666666.110134', '196627.31'], ['6666666.110155', '6666666.147110', '6666666.110147', '196627.31']], 'Arthrobacter_PYR_THR_AA': [['6666666.104489', '6666666.104435', '6666666.104519', '452863.24'], ['6666666.104524', '6666666.104488', '6666666.104439', '452863.24'], ['930171.11', '6666666.146813', '6666666.147082', '6666666.104484'], ['6666666.104502', '6666666.104519', '6666666.104488', '290340.24']], 'Bifidobacterium_PYR_THR_AA': [['6666666.112175', '442563.11', '6666666.112173', '367928.21'], ['6666666.112175', '6666666.112206', '6666666.112168', '367928.21']], 'Kitasatospora_PYR_THR_AA': [['6666666.113355', '6666666.113222', '6666666.113026', '6666666.112815']], 'Gordonia_PYR_THR_AA': [['6666666.111019', '6666666.111018', '6666666.111006', '6666666.110987'], ['6666666.111006', '6666666.110999', '6666666.111003', '6666666.110978']], 'Actinoplanes_PYR_THR_AA': [['6666666.146850', '649831.6', '6666666.104653', '6666666.104555']], 'Streptomyces_PYR_THR_AA': [['6666666.113025', '6666666.103215', '6666666.104136', '6666666.103170'], ['6666666.103267', '6666666.103281', '6666666.103171', '6666666.102743'], ['6666666.113124', '6666666.103261', '6666666.103187', '6666666.102809'], ['6666666.113171', '6666666.113099', '6666666.103192', '6666666.103187'], ['6666666.104222', '6666666.113099', '6666666.104202', '6666666.103201'], ['6666666.104188', '6666666.112961', '6666666.103248', '6666666.102743'], ['6666666.113138', '6666666.104229', '6666666.103336', '6666666.103155'], ['6666666.104313', '6666666.104191', '6666666.103298', '6666666.103155'], ['6666666.113174', '6666666.104313', '6666666.112864', '6666666.104375', '6666666.103212', '6666666.103208'], ['6666666.113171', '6666666.113155', '6666666.112956', '6666666.102800'], ['6666666.104222', '6666666.113171', '6666666.102800', '6666666.102743'], ['6666666.104177', '6666666.103300', '6666666.103336', '6666666.103155'], ['6666666.104188', '6666666.113143', '6666666.103249', '6666666.102743'], ['6666666.104222', '6666666.104349', '6666666.103248', '6666666.102743'], ['6666666.112839', '6666666.103284', '6666666.103187', '6666666.102809'], ['6666666.103249', '6666666.104177', '6666666.103155', '6666666.102743'], ['6666666.104136', '6666666.103264', '6666666.103267', '6666666.102743'], ['6666666.103336', '6666666.113007', '6666666.103163', '6666666.103155'], ['6666666.113124', '6666666.104235', '6666666.103298', '6666666.102809'], ['6666666.113122', '6666666.103281', '6666666.103171', '6666666.103163'], ['6666666.113243', '6666666.147155', '6666666.113007', '6666666.103163'], ['6666666.112828', '6666666.104224', '6666666.103178', '6666666.103208'], ['6666666.104375', '6666666.104345', '6666666.104120', '6666666.103212'], ['6666666.104222', '6666666.103336', '6666666.103155', '6666666.102743'], ['6666666.112839', '6666666.113243', '6666666.103171', '6666666.102743'], ['6666666.112839', '6666666.103257', '6666666.103249', '6666666.102743'], ['6666666.104136', '6666666.113138', '6666666.103155', '6666666.102743'], ['6666666.113206', '6666666.112979', '6666666.103258', '6666666.103171'], ['6666666.104123', '6666666.147078', '6666666.103248', '6666666.102743'], ['6666666.112956', '6666666.104191', '6666666.103298', '6666666.102800'], ['6666666.112839', '6666666.104321', '6666666.103248', '6666666.102743'], ['6666666.113207', '6666666.113067', '6666666.104235', '6666666.103155'], ['6666666.104123', '6666666.113148', '6666666.102800', '6666666.102743'], ['6666666.104313', '6666666.113174', '6666666.103201', '6666666.103267'], ['6666666.104218', '6666666.104164', '6666666.103187', '6666666.102809'], ['6666666.113138', '6666666.113178', '6666666.104313', '6666666.103155'], ['6666666.103267', '6666666.112956', '6666666.102800', '6666666.102743'], ['6666666.113138', '6666666.104229', '6666666.113025', '6666666.103170'], ['6666666.113243', '6666666.112835', '6666666.104715', '6666666.103170'], ['6666666.104708', '6666666.103253', '6666666.104321', '6666666.103156'], ['6666666.113127', '6666666.104668', '6666666.113025', '6666666.103170'], ['6666666.103300', '6666666.103336', '6666666.103284', '6666666.103192'], ['6666666.113279', '6666666.104349', '6666666.113099', '6666666.104325'], ['6666666.113124', '6666666.113146', '6666666.104135', '6666666.102809'], ['6666666.103267', '6666666.113055', '6666666.103248', '6666666.102743'], ['6666666.103187', '6666666.113025', '6666666.103170', '6666666.102809'], ['6666666.103336', '6666666.104708', '6666666.103156', '6666666.103155'], ['6666666.103257', '6666666.103192', '6666666.103187', '6666666.102809'], ['6666666.104313', '6666666.113178', '6666666.103264', '6666666.103267'], ['6666666.113207', '6666666.104714', '6666666.112884', '6666666.103155'], ['6666666.113024', '6666666.112931', '6666666.103264', '6666666.103267'], ['6666666.104313', '6666666.113174', '6666666.103336', '6666666.103155'], ['6666666.104136', '6666666.104336', '6666666.104123', '6666666.102743'], ['6666666.103171', '6666666.103163', '6666666.103155', '6666666.102743'], ['6666666.113127', '6666666.113169', '6666666.112876', '6666666.103170'], ['6666666.104188', '6666666.113041', '6666666.103267', '6666666.102743'], ['6666666.113279', '6666666.103255', '6666666.113055', '6666666.103248'], ['6666666.112839', '6666666.103170', '6666666.104136', '6666666.102743'], ['6666666.113207', '6666666.113087', '6666666.112839', '6666666.103155'], ['6666666.104321', '6666666.103253', '6666666.103187', '6666666.102809'], ['6666666.103287', '6666666.112819', '6666666.103281', '6666666.103171'], ['6666666.112989', '6666666.113287', '6666666.103300', '6666666.103192'], ['6666666.112839', '6666666.102809', '6666666.102800', '6666666.102743'], ['6666666.112956', '6666666.113155', '6666666.103201', '6666666.103267'], ['6666666.104222', '6666666.103287', '6666666.103171', '6666666.102743'], ['6666666.147092', '6666666.103247', '6666666.112839', '6666666.102743'], ['6666666.112828', '6666666.104224', '6666666.104183', '6666666.103201'], ['6666666.112876', '6666666.113087', '6666666.112839', '6666666.103170'], ['6666666.104188', '6666666.113127', '6666666.104136', '6666666.102743'], ['6666666.112839', '6666666.103284', '6666666.104222', '6666666.102743'], ['6666666.113243', '6666666.147155', '6666666.103187', '6666666.102809'], ['6666666.119638', '6666666.104202', '6666666.113099', '6666666.103192'], ['6666666.103264', '6666666.103270', '6666666.103201', '6666666.103267'], ['6666666.113007', '6666666.103287', '6666666.103171', '6666666.103163'], ['6666666.104313', '6666666.113122', '6666666.103163', '6666666.103155'], ['6666666.104222', '6666666.113099', '6666666.103249', '6666666.102743'], ['6666666.113171', '6666666.103187', '6666666.102809', '6666666.102800'], ['6666666.112821', '6666666.104338', '6666666.104348', '6666666.113273', '6666666.103212'], ['6666666.103248', '6666666.103156', '6666666.103155', '6666666.102743'], ['6666666.104222', '6666666.103215', '6666666.104136', '6666666.102743'], ['6666666.104222', '6666666.103215', '6666666.103270', '6666666.103201'], ['6666666.103155', '6666666.103298', '6666666.102800', '6666666.102743'], ['6666666.147054', '6666666.103202', '6666666.147155', '6666666.103163'], ['6666666.103281', '6666666.112819', '6666666.103201', '6666666.103267'], ['6666666.112884', '6666666.113287', '6666666.103336', '6666666.103155'], ['6666666.103267', '6666666.104313', '6666666.103155', '6666666.102743'], ['6666666.113207', '6666666.113087', '6666666.113116', '6666666.104714'], ['6666666.104188', '6666666.104135', '6666666.102800', '6666666.102743'], ['6666666.104708', '6666666.104349', '6666666.103248', '6666666.103156'], ['6666666.104222', '6666666.103201', '6666666.103267', '6666666.102743']], 'Microbacterium_PYR_THR_AA': [['6666666.104456', '6666666.104495', '6666666.104451', '6666666.104444']], 'Saccharomonospora_PYR_THR_AA': [['6666666.146841', '6666666.147079', '6666666.104473', '6666666.112808']], 'Pseudonocardia_PYR_THR_AA': [['6666666.104623', '6666666.112804', '6666666.104596', '6666666.104550']], 'Brevibacterium_PYR_THR_AA': [['6666666.104562', '6666666.104598', '6666666.104557', '6666666.111558']], 'Actinomyces_PYR_THR_AA': [['6666666.114196', '6666666.114195', '6666666.114192', '6666666.114191'], ['6666666.114208', '6666666.146867', '6666666.114192', '6666666.114191']], 'Nocardia_R5P_AMINOACIDS': [['6666666.102760', '6666666.102749', '6666666.102747', '1127134.4']], 'Streptomyces_R5P_AMINOACIDS': [['6666666.103178', '6666666.103153', '6666666.102800', '6666666.102743'], ['6666666.146818', '6666666.112819', '6666666.104132', '6666666.102743'], ['6666666.103182', '6666666.104325', '6666666.102800', '6666666.102743'], ['6666666.103257', '6666666.112828', '6666666.102800', '6666666.102743'], ['6666666.146818', '6666666.113041', '6666666.102800', '6666666.102743'], ['6666666.112835', '6666666.103215', '6666666.102800', '6666666.102743'], ['6666666.104132', '6666666.113037', '6666666.102800', '6666666.102743']], 'Nocardia_TCA': [['6666666.102753', '1127134.4', '6666666.102742', '6666666.102741']], 'Amycolatopsis_TCA': [['6666666.104475', '6666666.104442', '6666666.104440', '1156913.7'], ['6666666.104454', '6666666.104471', '6666666.104440', '1156913.7']], 'Propionibacterium_TCA': [['6666666.105581', '6666666.105576', '6666666.103566', '1170318.5']], 'Mycobacterium_TCA': [['6666666.111038', '6666666.110994', '6666666.101415', '1299323.4'], ['6666666.101415', '6666666.110968', '350054.22', '1299323.4'], ['6666666.101415', '6666666.110984', '246196.47', '1299323.4'], ['6666666.111038', '6666666.111070', '350054.22', '1299323.4'], ['6666666.110963', '6666666.111043', '350054.22', '1299323.4']], 'Corynebacterium_TCA': [['6666666.110147', '6666666.110136', '504474.6', '196627.31'], ['6666666.110147', '6666666.110180', '645127.7', '196627.31']], 'Kitasatospora_TCA': [['6666666.113188', '6666666.113222', '6666666.112815', '452652.7']], 'Gordonia_TCA': [['6666666.111025', '6666666.110978', '6666666.110976', '526226.13'], ['6666666.111014', '6666666.110979', '6666666.110978', '6666666.110976'], ['6666666.111014', '6666666.111069', '6666666.110987', '6666666.110976']], 'Streptomyces_TCA': [['6666666.112847', '6666666.102809', '6666666.103156', '6666666.102743'], ['6666666.112847', '6666666.113206', '6666666.102800', '6666666.102743'], ['6666666.113129', '6666666.147125', '6666666.112847', '6666666.102743'], ['6666666.112847', '6666666.103208', '6666666.103178', '6666666.102743'], ['6666666.103310', '6666666.113283', '6666666.103155', '6666666.102809'], ['6666666.103182', '6666666.103201', '6666666.103155', '6666666.102809'], ['6666666.103207', '6666666.104323', '6666666.103156', '6666666.102743'], ['6666666.103178', '6666666.103246', '6666666.103156', '6666666.102743'], ['6666666.103246', '6666666.103153', '6666666.103170', '6666666.103155'], ['6666666.112847', '6666666.103208', '6666666.103155', '6666666.102809'], ['6666666.103298', '6666666.104323', '6666666.103156', '6666666.102809'], ['6666666.103156', '6666666.103246', '6666666.103155', '6666666.102809'], ['6666666.103208', '6666666.104178', '6666666.103201', '6666666.103155'], ['6666666.113178', '6666666.103278', '6666666.103246', '6666666.103155'], ['6666666.113129', '6666666.103153', '6666666.103178', '6666666.102743'], ['6666666.147125', '6666666.103170', '6666666.103155', '6666666.102809']], 'Microbacterium_TCA': [['6666666.147166', '6666666.104495', '6666666.104444', '6666666.104436']], 'Pseudonocardia_TCA': [['6666666.104611', '6666666.104550', '6666666.104596', '6666666.104554']], 'Actinomyces_TCA': [['6666666.114203', '6666666.114195', '6666666.114192', '6666666.114191'], ['6666666.114203', '6666666.114242', '6666666.114198', '6666666.114195']], 'Rhodococcus_calcium_dependent_antibiotic': [], 'Nocardia_calcium_dependent_antibiotic': [], 'Amycolatopsis_calcium_dependent_antibiotic': [], 'Propionibacterium_calcium_dependent_antibiotic': [], 'Nocardiopsis_calcium_dependent_antibiotic': [], 'Clavibacter_calcium_dependent_antibiotic': [], 'Mycobacterium_calcium_dependent_antibiotic': [], 'Kutzneria_calcium_dependent_antibiotic': [], 'Citricoccus_calcium_dependent_antibiotic': [], 'Nocardioides_calcium_dependent_antibiotic': [], 'Corynebacterium_calcium_dependent_antibiotic': [], 'Arthrobacter_calcium_dependent_antibiotic': [], 'Frankia_calcium_dependent_antibiotic': [], 'Bifidobacterium_calcium_dependent_antibiotic': [], 'Salinispora_calcium_dependent_antibiotic': [], 'Kocuria_calcium_dependent_antibiotic': [], 'Actinosynnema_calcium_dependent_antibiotic': [], 'Kitasatospora_calcium_dependent_antibiotic': [], 'Gordonia_calcium_dependent_antibiotic': [], 'Micromonospora_calcium_dependent_antibiotic': [], 'Actinoplanes_calcium_dependent_antibiotic': [], 'Streptomyces_calcium_dependent_antibiotic': [], 'Microlunatus_calcium_dependent_antibiotic': [], 'Slackia_calcium_dependent_antibiotic': [], 'Cryptobacterium_calcium_dependent_antibiotic': [], 'Anaerolinea_calcium_dependent_antibiotic': [], 'Chloroflexus_calcium_dependent_antibiotic': [], 'Caldilinea_calcium_dependent_antibiotic': [], 'Roseiflexus_calcium_dependent_antibiotic': [], 'Dehalogenimonas_calcium_dependent_antibiotic': [], 'Thermomicrobium_calcium_dependent_antibiotic': [], 'Eggerthella_calcium_dependent_antibiotic': [], 'Atopobium_calcium_dependent_antibiotic': [], 'Dehalococcoides_calcium_dependent_antibiotic': [], 'Microbacterium_calcium_dependent_antibiotic': [], 'Saccharomonospora_calcium_dependent_antibiotic': [], 'Pseudonocardia_calcium_dependent_antibiotic': [], 'Brevibacterium_calcium_dependent_antibiotic': [], 'Dermabacter_calcium_dependent_antibiotic': [], 'Leifsonia_calcium_dependent_antibiotic': [], 'Saccharopolyspora_calcium_dependent_antibiotic': [], 'Leucobacter_calcium_dependent_antibiotic': [], 'Rothia_calcium_dependent_antibiotic': [], 'Agromyces_calcium_dependent_antibiotic': [], 'Actinoalloteichus_calcium_dependent_antibiotic': [], 'Sciscionella_calcium_dependent_antibiotic': [], 'Actinocatenispora_calcium_dependent_antibiotic': [], 'Lechevalieria_calcium_dependent_antibiotic': [], 'Thermobifida_calcium_dependent_antibiotic': [], 'Thermocrispum_calcium_dependent_antibiotic': [], 'Actinokineospora_calcium_dependent_antibiotic': [], 'Lentzea_calcium_dependent_antibiotic': [], 'Kibdelosporangium_calcium_dependent_antibiotic': [], 'Actinomycetospora_calcium_dependent_antibiotic': [], 'Saccharothrix_calcium_dependent_antibiotic': [], 'Collinsella_calcium_dependent_antibiotic': [], 'Propionimicrobium_calcium_dependent_antibiotic': [], 'Thermogemmatispora_calcium_dependent_antibiotic': [], 'Granulicoccus_calcium_dependent_antibiotic': [], 'Nitrolancea_calcium_dependent_antibiotic': [], 'Enterorhabdus_calcium_dependent_antibiotic': [], 'Olsenella_calcium_dependent_antibiotic': [], 'Senegalimassilia_calcium_dependent_antibiotic': [], 'Oscillochloris_calcium_dependent_antibiotic': [], 'Enorma_calcium_dependent_antibiotic': [], 'Aestuariimicrobium_calcium_dependent_antibiotic': [], 'Oerskovia_calcium_dependent_antibiotic': [], 'Cellulomonas_calcium_dependent_antibiotic': [], 'Paraoerskovia_calcium_dependent_antibiotic': [], 'Isoptericola_calcium_dependent_antibiotic': [], 'Brachybacterium_calcium_dependent_antibiotic': [], 'Promicromonospora_calcium_dependent_antibiotic': [], 'Cellulosimicrobium_calcium_dependent_antibiotic': [], 'Xylanimonas_calcium_dependent_antibiotic': [], 'Kineococcus_calcium_dependent_antibiotic': [], 'Kineosporia_calcium_dependent_antibiotic': [], 'Dietzia_calcium_dependent_antibiotic': [], 'Sphaerobacter_calcium_dependent_antibiotic': [], 'Ktedonobacter_calcium_dependent_antibiotic': [], 'Tomitella_calcium_dependent_antibiotic': [], 'Segniliparus_calcium_dependent_antibiotic': [], 'Catenulispora_calcium_dependent_antibiotic': [], 'Tsukamurella_calcium_dependent_antibiotic': [], 'Actinopolyspora_calcium_dependent_antibiotic': [], 'Williamsia_calcium_dependent_antibiotic': [], 'Actinospica_calcium_dependent_antibiotic': [], 'Tropheryma_calcium_dependent_antibiotic': [], 'Beutenbergia_calcium_dependent_antibiotic': [], 'Kytococcus_calcium_dependent_antibiotic': [], 'Nakamurella_calcium_dependent_antibiotic': [], 'Blastococcus_calcium_dependent_antibiotic': [], 'Candidatus_calcium_dependent_antibiotic': [], 'Acidothermus_calcium_dependent_antibiotic': [], 'Lysinimicrobium_calcium_dependent_antibiotic': [], 'Geodermatophilaceae_calcium_dependent_antibiotic': [], 'Glycomyces_calcium_dependent_antibiotic': [], 'Cryptosporangium_calcium_dependent_antibiotic': [], 'Ruania_calcium_dependent_antibiotic': [], 'Modestobacter_calcium_dependent_antibiotic': [], 'Dermatophilus_calcium_dependent_antibiotic': [], 'Dermacoccus_calcium_dependent_antibiotic': [], 'Mobilicoccus_calcium_dependent_antibiotic': [], 'Demetria_calcium_dependent_antibiotic': [], 'Sanguibacter_calcium_dependent_antibiotic': [], 'Kineosphaera_calcium_dependent_antibiotic': [], 'Timonella_calcium_dependent_antibiotic': [], 'Geodermatophilus_calcium_dependent_antibiotic': [], 'Haloglycomyces_calcium_dependent_antibiotic': [], 'Sporichthya_calcium_dependent_antibiotic': [], 'Jonesia_calcium_dependent_antibiotic': [], 'Janibacter_calcium_dependent_antibiotic': [], 'Gryllotalpicola_calcium_dependent_antibiotic': [], 'Frigoribacterium_calcium_dependent_antibiotic': [], 'Herbiconiux_calcium_dependent_antibiotic': [], 'Knoellia_calcium_dependent_antibiotic': [], 'Humibacter_calcium_dependent_antibiotic': [], 'Serinicoccus_calcium_dependent_antibiotic': [], 'Rathayibacter_calcium_dependent_antibiotic': [], 'Agreia_calcium_dependent_antibiotic': [], 'Cryobacterium_calcium_dependent_antibiotic': [], 'Intrasporangium_calcium_dependent_antibiotic': [], 'Tetrasphaera_calcium_dependent_antibiotic': [], 'Intrasporangiaceae_calcium_dependent_antibiotic': [], 'Glaciibacter_calcium_dependent_antibiotic': [], 'Cryocola_calcium_dependent_antibiotic': [], 'Austwickia_calcium_dependent_antibiotic': [], 'Arsenicicoccus_calcium_dependent_antibiotic': [], 'Salinibacterium_calcium_dependent_antibiotic': [], 'Gulosibacter_calcium_dependent_antibiotic': [], 'Pseudoclavibacter_calcium_dependent_antibiotic': [], 'Mycetocola_calcium_dependent_antibiotic': [], 'Phycicoccus_calcium_dependent_antibiotic': [], 'Ornithinimicrobium_calcium_dependent_antibiotic': [], 'Terrabacter_calcium_dependent_antibiotic': [], 'Terracoccus_calcium_dependent_antibiotic': [], 'Alloscardovia_calcium_dependent_antibiotic': [], 'Actinopolymorpha_calcium_dependent_antibiotic': [], 'Kribbella_calcium_dependent_antibiotic': [], 'Agrococcus_calcium_dependent_antibiotic': [], 'Micrococcus_calcium_dependent_antibiotic': [], 'Nesterenkonia_calcium_dependent_antibiotic': [], 'Dactylosporangium_calcium_dependent_antibiotic': [], 'Nocardioidaceae_calcium_dependent_antibiotic': [], 'Aeromicrobium_calcium_dependent_antibiotic': [], 'Renibacterium_calcium_dependent_antibiotic': [], 'Verrucosispora_calcium_dependent_antibiotic': [], 'Sinomonas_calcium_dependent_antibiotic': [], 'Marmoricola_calcium_dependent_antibiotic': [], 'Pilimelia_calcium_dependent_antibiotic': [], 'Catelliglobosispora_calcium_dependent_antibiotic': [], 'Acaricomes_calcium_dependent_antibiotic': [], 'Catenuloplanes_calcium_dependent_antibiotic': [], 'Jiangella_calcium_dependent_antibiotic': [], 'Longispora_calcium_dependent_antibiotic': [], 'Yaniella_calcium_dependent_antibiotic': [], 'Actinomadura_calcium_dependent_antibiotic': [], 'Microbispora_calcium_dependent_antibiotic': [], 'Nonomuraea_calcium_dependent_antibiotic': [], 'Streptosporangium_calcium_dependent_antibiotic': [], 'Streptomonospora_calcium_dependent_antibiotic': [], 'Microtetraspora_calcium_dependent_antibiotic': [], 'Prauserella_calcium_dependent_antibiotic': [], 'Allokutzneria_calcium_dependent_antibiotic': [], 'Streptacidiphilus_calcium_dependent_antibiotic': [], 'Mobiluncus_calcium_dependent_antibiotic': [], 'Actinomyces_calcium_dependent_antibiotic': [], 'Trueperella_calcium_dependent_antibiotic': [], 'Varibaculum_calcium_dependent_antibiotic': [], 'Actinobaculum_calcium_dependent_antibiotic': [], 'Arcanobacterium_calcium_dependent_antibiotic': [], 'Pimelobacter_calcium_dependent_antibiotic': [], 'Corynebacteriales_calcium_dependent_antibiotic': [], 'Parascardovia_calcium_dependent_antibiotic': [], 'Thermomonospora_calcium_dependent_antibiotic': [], 'Luteipulveratus_calcium_dependent_antibiotic': [], 'Actinotignum_calcium_dependent_antibiotic': [], 'Alloactinosynnema_calcium_dependent_antibiotic': [], 'Amycolicicoccus_calcium_dependent_antibiotic': [], 'Actinobacteria_calcium_dependent_antibiotic': [], 'Scardovia_calcium_dependent_antibiotic': [], 'Adlercreutzia_calcium_dependent_antibiotic': [], 'Coriobacteriaceae_calcium_dependent_antibiotic': [], 'Gardnerella_calcium_dependent_antibiotic': [], 'Rubrobacter_calcium_dependent_antibiotic': [], 'Gordonibacter_calcium_dependent_antibiotic': [], 'Stackebrandtia_calcium_dependent_antibiotic': [], 'Thermobispora_calcium_dependent_antibiotic': [], 'Coriobacterium_calcium_dependent_antibiotic': [], 'Conexibacter_calcium_dependent_antibiotic': [], 'Ilumatobacter_calcium_dependent_antibiotic': [], 'Cellvibrio_calcium_dependent_antibiotic': [], 'Acidimicrobium_calcium_dependent_antibiotic': [], 'Actinokieospora_calcium_dependent_antibiotic': [], 'Rhodococcus_polyketide': [], 'Nocardia_polyketide': [], 'Amycolatopsis_polyketide': [], 'Propionibacterium_polyketide': [], 'Nocardiopsis_polyketide': [], 'Clavibacter_polyketide': [], 'Mycobacterium_polyketide': [], 'Kutzneria_polyketide': [], 'Citricoccus_polyketide': [], 'Nocardioides_polyketide': [], 'Corynebacterium_polyketide': [], 'Arthrobacter_polyketide': [], 'Frankia_polyketide': [], 'Bifidobacterium_polyketide': [], 'Salinispora_polyketide': [], 'Kocuria_polyketide': [], 'Actinosynnema_polyketide': [], 'Kitasatospora_polyketide': [], 'Gordonia_polyketide': [], 'Micromonospora_polyketide': [], 'Actinoplanes_polyketide': [], 'Streptomyces_polyketide': [], 'Microlunatus_polyketide': [], 'Slackia_polyketide': [], 'Cryptobacterium_polyketide': [], 'Anaerolinea_polyketide': [], 'Chloroflexus_polyketide': [], 'Caldilinea_polyketide': [], 'Roseiflexus_polyketide': [], 'Dehalogenimonas_polyketide': [], 'Thermomicrobium_polyketide': [], 'Eggerthella_polyketide': [], 'Atopobium_polyketide': [], 'Dehalococcoides_polyketide': [], 'Microbacterium_polyketide': [], 'Saccharomonospora_polyketide': [], 'Pseudonocardia_polyketide': [], 'Brevibacterium_polyketide': [], 'Dermabacter_polyketide': [], 'Leifsonia_polyketide': [], 'Saccharopolyspora_polyketide': [], 'Leucobacter_polyketide': [], 'Rothia_polyketide': [], 'Agromyces_polyketide': [], 'Actinoalloteichus_polyketide': [], 'Sciscionella_polyketide': [], 'Actinocatenispora_polyketide': [], 'Lechevalieria_polyketide': [], 'Thermobifida_polyketide': [], 'Thermocrispum_polyketide': [], 'Actinokineospora_polyketide': [], 'Lentzea_polyketide': [], 'Kibdelosporangium_polyketide': [], 'Actinomycetospora_polyketide': [], 'Saccharothrix_polyketide': [], 'Collinsella_polyketide': [], 'Propionimicrobium_polyketide': [], 'Thermogemmatispora_polyketide': [], 'Granulicoccus_polyketide': [], 'Nitrolancea_polyketide': [], 'Enterorhabdus_polyketide': [], 'Olsenella_polyketide': [], 'Senegalimassilia_polyketide': [], 'Oscillochloris_polyketide': [], 'Enorma_polyketide': [], 'Aestuariimicrobium_polyketide': [], 'Oerskovia_polyketide': [], 'Cellulomonas_polyketide': [], 'Paraoerskovia_polyketide': [], 'Isoptericola_polyketide': [], 'Brachybacterium_polyketide': [], 'Promicromonospora_polyketide': [], 'Cellulosimicrobium_polyketide': [], 'Xylanimonas_polyketide': [], 'Kineococcus_polyketide': [], 'Kineosporia_polyketide': [], 'Dietzia_polyketide': [], 'Sphaerobacter_polyketide': [], 'Ktedonobacter_polyketide': [], 'Tomitella_polyketide': [], 'Segniliparus_polyketide': [], 'Catenulispora_polyketide': [], 'Tsukamurella_polyketide': [], 'Actinopolyspora_polyketide': [], 'Williamsia_polyketide': [], 'Actinospica_polyketide': [], 'Tropheryma_polyketide': [], 'Beutenbergia_polyketide': [], 'Kytococcus_polyketide': [], 'Nakamurella_polyketide': [], 'Blastococcus_polyketide': [], 'Candidatus_polyketide': [], 'Acidothermus_polyketide': [], 'Lysinimicrobium_polyketide': [], 'Geodermatophilaceae_polyketide': [], 'Glycomyces_polyketide': [], 'Cryptosporangium_polyketide': [], 'Ruania_polyketide': [], 'Modestobacter_polyketide': [], 'Dermatophilus_polyketide': [], 'Dermacoccus_polyketide': [], 'Mobilicoccus_polyketide': [], 'Demetria_polyketide': [], 'Sanguibacter_polyketide': [], 'Kineosphaera_polyketide': [], 'Timonella_polyketide': [], 'Geodermatophilus_polyketide': [], 'Haloglycomyces_polyketide': [], 'Sporichthya_polyketide': [], 'Jonesia_polyketide': [], 'Janibacter_polyketide': [], 'Gryllotalpicola_polyketide': [], 'Frigoribacterium_polyketide': [], 'Herbiconiux_polyketide': [], 'Knoellia_polyketide': [], 'Humibacter_polyketide': [], 'Serinicoccus_polyketide': [], 'Rathayibacter_polyketide': [], 'Agreia_polyketide': [], 'Cryobacterium_polyketide': [], 'Intrasporangium_polyketide': [], 'Tetrasphaera_polyketide': [], 'Intrasporangiaceae_polyketide': [], 'Glaciibacter_polyketide': [], 'Cryocola_polyketide': [], 'Austwickia_polyketide': [], 'Arsenicicoccus_polyketide': [], 'Salinibacterium_polyketide': [], 'Gulosibacter_polyketide': [], 'Pseudoclavibacter_polyketide': [], 'Mycetocola_polyketide': [], 'Phycicoccus_polyketide': [], 'Ornithinimicrobium_polyketide': [], 'Terrabacter_polyketide': [], 'Terracoccus_polyketide': [], 'Alloscardovia_polyketide': [], 'Actinopolymorpha_polyketide': [], 'Kribbella_polyketide': [], 'Agrococcus_polyketide': [], 'Micrococcus_polyketide': [], 'Nesterenkonia_polyketide': [], 'Dactylosporangium_polyketide': [], 'Nocardioidaceae_polyketide': [], 'Aeromicrobium_polyketide': [], 'Renibacterium_polyketide': [], 'Verrucosispora_polyketide': [], 'Sinomonas_polyketide': [], 'Marmoricola_polyketide': [], 'Pilimelia_polyketide': [], 'Catelliglobosispora_polyketide': [], 'Acaricomes_polyketide': [], 'Catenuloplanes_polyketide': [], 'Jiangella_polyketide': [], 'Longispora_polyketide': [], 'Yaniella_polyketide': [], 'Actinomadura_polyketide': [], 'Microbispora_polyketide': [], 'Nonomuraea_polyketide': [], 'Streptosporangium_polyketide': [], 'Streptomonospora_polyketide': [], 'Microtetraspora_polyketide': [], 'Prauserella_polyketide': [], 'Allokutzneria_polyketide': [], 'Streptacidiphilus_polyketide': [], 'Mobiluncus_polyketide': [], 'Actinomyces_polyketide': [], 'Trueperella_polyketide': [], 'Varibaculum_polyketide': [], 'Actinobaculum_polyketide': [], 'Arcanobacterium_polyketide': [], 'Pimelobacter_polyketide': [], 'Corynebacteriales_polyketide': [], 'Parascardovia_polyketide': [], 'Thermomonospora_polyketide': [], 'Luteipulveratus_polyketide': [], 'Actinotignum_polyketide': [], 'Alloactinosynnema_polyketide': [], 'Amycolicicoccus_polyketide': [], 'Actinobacteria_polyketide': [], 'Scardovia_polyketide': [], 'Adlercreutzia_polyketide': [], 'Coriobacteriaceae_polyketide': [], 'Gardnerella_polyketide': [], 'Rubrobacter_polyketide': [], 'Gordonibacter_polyketide': [], 'Stackebrandtia_polyketide': [], 'Thermobispora_polyketide': [], 'Coriobacterium_polyketide': [], 'Conexibacter_polyketide': [], 'Ilumatobacter_polyketide': [], 'Cellvibrio_polyketide': [], 'Acidimicrobium_polyketide': [], 'Actinokieospora_polyketide': [], 'Propionibacterium_3PGA_AMINOACIDS': [], 'Nocardiopsis_3PGA_AMINOACIDS': [], 'Clavibacter_3PGA_AMINOACIDS': [], 'Kutzneria_3PGA_AMINOACIDS': [], 'Citricoccus_3PGA_AMINOACIDS': [], 'Nocardioides_3PGA_AMINOACIDS': [], 'Corynebacterium_3PGA_AMINOACIDS': [], 'Frankia_3PGA_AMINOACIDS': [], 'Bifidobacterium_3PGA_AMINOACIDS': [], 'Salinispora_3PGA_AMINOACIDS': [], 'Actinosynnema_3PGA_AMINOACIDS': [], 'Kitasatospora_3PGA_AMINOACIDS': [], 'Micromonospora_3PGA_AMINOACIDS': [], 'Actinoplanes_3PGA_AMINOACIDS': [], 'Microlunatus_3PGA_AMINOACIDS': [], 'Slackia_3PGA_AMINOACIDS': [], 'Cryptobacterium_3PGA_AMINOACIDS': [], 'Anaerolinea_3PGA_AMINOACIDS': [], 'Chloroflexus_3PGA_AMINOACIDS': [], 'Caldilinea_3PGA_AMINOACIDS': [], 'Roseiflexus_3PGA_AMINOACIDS': [], 'Dehalogenimonas_3PGA_AMINOACIDS': [], 'Thermomicrobium_3PGA_AMINOACIDS': [], 'Eggerthella_3PGA_AMINOACIDS': [], 'Atopobium_3PGA_AMINOACIDS': [], 'Dehalococcoides_3PGA_AMINOACIDS': [], 'Microbacterium_3PGA_AMINOACIDS': [], 'Pseudonocardia_3PGA_AMINOACIDS': [], 'Dermabacter_3PGA_AMINOACIDS': [], 'Leifsonia_3PGA_AMINOACIDS': [], 'Saccharopolyspora_3PGA_AMINOACIDS': [], 'Leucobacter_3PGA_AMINOACIDS': [], 'Rothia_3PGA_AMINOACIDS': [], 'Agromyces_3PGA_AMINOACIDS': [], 'Actinoalloteichus_3PGA_AMINOACIDS': [], 'Sciscionella_3PGA_AMINOACIDS': [], 'Actinocatenispora_3PGA_AMINOACIDS': [], 'Lechevalieria_3PGA_AMINOACIDS': [], 'Thermobifida_3PGA_AMINOACIDS': [], 'Thermocrispum_3PGA_AMINOACIDS': [], 'Actinokineospora_3PGA_AMINOACIDS': [], 'Lentzea_3PGA_AMINOACIDS': [], 'Kibdelosporangium_3PGA_AMINOACIDS': [], 'Actinomycetospora_3PGA_AMINOACIDS': [], 'Saccharothrix_3PGA_AMINOACIDS': [], 'Collinsella_3PGA_AMINOACIDS': [], 'Propionimicrobium_3PGA_AMINOACIDS': [], 'Thermogemmatispora_3PGA_AMINOACIDS': [], 'Granulicoccus_3PGA_AMINOACIDS': [], 'Nitrolancea_3PGA_AMINOACIDS': [], 'Enterorhabdus_3PGA_AMINOACIDS': [], 'Olsenella_3PGA_AMINOACIDS': [], 'Senegalimassilia_3PGA_AMINOACIDS': [], 'Oscillochloris_3PGA_AMINOACIDS': [], 'Enorma_3PGA_AMINOACIDS': [], 'Aestuariimicrobium_3PGA_AMINOACIDS': [], 'Oerskovia_3PGA_AMINOACIDS': [], 'Cellulomonas_3PGA_AMINOACIDS': [], 'Paraoerskovia_3PGA_AMINOACIDS': [], 'Isoptericola_3PGA_AMINOACIDS': [], 'Brachybacterium_3PGA_AMINOACIDS': [], 'Promicromonospora_3PGA_AMINOACIDS': [], 'Cellulosimicrobium_3PGA_AMINOACIDS': [], 'Xylanimonas_3PGA_AMINOACIDS': [], 'Kineococcus_3PGA_AMINOACIDS': [], 'Kineosporia_3PGA_AMINOACIDS': [], 'Dietzia_3PGA_AMINOACIDS': [], 'Sphaerobacter_3PGA_AMINOACIDS': [], 'Ktedonobacter_3PGA_AMINOACIDS': [], 'Tomitella_3PGA_AMINOACIDS': [], 'Segniliparus_3PGA_AMINOACIDS': [], 'Catenulispora_3PGA_AMINOACIDS': [], 'Tsukamurella_3PGA_AMINOACIDS': [], 'Actinopolyspora_3PGA_AMINOACIDS': [], 'Williamsia_3PGA_AMINOACIDS': [], 'Actinospica_3PGA_AMINOACIDS': [], 'Tropheryma_3PGA_AMINOACIDS': [], 'Beutenbergia_3PGA_AMINOACIDS': [], 'Kytococcus_3PGA_AMINOACIDS': [], 'Nakamurella_3PGA_AMINOACIDS': [], 'Blastococcus_3PGA_AMINOACIDS': [], 'Candidatus_3PGA_AMINOACIDS': [], 'Acidothermus_3PGA_AMINOACIDS': [], 'Lysinimicrobium_3PGA_AMINOACIDS': [], 'Geodermatophilaceae_3PGA_AMINOACIDS': [], 'Glycomyces_3PGA_AMINOACIDS': [], 'Cryptosporangium_3PGA_AMINOACIDS': [], 'Ruania_3PGA_AMINOACIDS': [], 'Modestobacter_3PGA_AMINOACIDS': [], 'Dermatophilus_3PGA_AMINOACIDS': [], 'Dermacoccus_3PGA_AMINOACIDS': [], 'Mobilicoccus_3PGA_AMINOACIDS': [], 'Demetria_3PGA_AMINOACIDS': [], 'Sanguibacter_3PGA_AMINOACIDS': [], 'Kineosphaera_3PGA_AMINOACIDS': [], 'Timonella_3PGA_AMINOACIDS': [], 'Geodermatophilus_3PGA_AMINOACIDS': [], 'Haloglycomyces_3PGA_AMINOACIDS': [], 'Sporichthya_3PGA_AMINOACIDS': [], 'Jonesia_3PGA_AMINOACIDS': [], 'Janibacter_3PGA_AMINOACIDS': [], 'Gryllotalpicola_3PGA_AMINOACIDS': [], 'Frigoribacterium_3PGA_AMINOACIDS': [], 'Herbiconiux_3PGA_AMINOACIDS': [], 'Knoellia_3PGA_AMINOACIDS': [], 'Humibacter_3PGA_AMINOACIDS': [], 'Serinicoccus_3PGA_AMINOACIDS': [], 'Rathayibacter_3PGA_AMINOACIDS': [], 'Agreia_3PGA_AMINOACIDS': [], 'Cryobacterium_3PGA_AMINOACIDS': [], 'Intrasporangium_3PGA_AMINOACIDS': [], 'Tetrasphaera_3PGA_AMINOACIDS': [], 'Intrasporangiaceae_3PGA_AMINOACIDS': [], 'Glaciibacter_3PGA_AMINOACIDS': [], 'Cryocola_3PGA_AMINOACIDS': [], 'Austwickia_3PGA_AMINOACIDS': [], 'Arsenicicoccus_3PGA_AMINOACIDS': [], 'Salinibacterium_3PGA_AMINOACIDS': [], 'Gulosibacter_3PGA_AMINOACIDS': [], 'Pseudoclavibacter_3PGA_AMINOACIDS': [], 'Mycetocola_3PGA_AMINOACIDS': [], 'Phycicoccus_3PGA_AMINOACIDS': [], 'Ornithinimicrobium_3PGA_AMINOACIDS': [], 'Terrabacter_3PGA_AMINOACIDS': [], 'Terracoccus_3PGA_AMINOACIDS': [], 'Alloscardovia_3PGA_AMINOACIDS': [], 'Actinopolymorpha_3PGA_AMINOACIDS': [], 'Kribbella_3PGA_AMINOACIDS': [], 'Agrococcus_3PGA_AMINOACIDS': [], 'Micrococcus_3PGA_AMINOACIDS': [], 'Nesterenkonia_3PGA_AMINOACIDS': [], 'Dactylosporangium_3PGA_AMINOACIDS': [], 'Nocardioidaceae_3PGA_AMINOACIDS': [], 'Aeromicrobium_3PGA_AMINOACIDS': [], 'Renibacterium_3PGA_AMINOACIDS': [], 'Verrucosispora_3PGA_AMINOACIDS': [], 'Sinomonas_3PGA_AMINOACIDS': [], 'Marmoricola_3PGA_AMINOACIDS': [], 'Pilimelia_3PGA_AMINOACIDS': [], 'Catelliglobosispora_3PGA_AMINOACIDS': [], 'Acaricomes_3PGA_AMINOACIDS': [], 'Catenuloplanes_3PGA_AMINOACIDS': [], 'Jiangella_3PGA_AMINOACIDS': [], 'Longispora_3PGA_AMINOACIDS': [], 'Yaniella_3PGA_AMINOACIDS': [], 'Actinomadura_3PGA_AMINOACIDS': [], 'Microbispora_3PGA_AMINOACIDS': [], 'Nonomuraea_3PGA_AMINOACIDS': [], 'Streptosporangium_3PGA_AMINOACIDS': [], 'Streptomonospora_3PGA_AMINOACIDS': [], 'Microtetraspora_3PGA_AMINOACIDS': [], 'Prauserella_3PGA_AMINOACIDS': [], 'Allokutzneria_3PGA_AMINOACIDS': [], 'Mobiluncus_3PGA_AMINOACIDS': [], 'Actinomyces_3PGA_AMINOACIDS': [], 'Trueperella_3PGA_AMINOACIDS': [], 'Varibaculum_3PGA_AMINOACIDS': [], 'Actinobaculum_3PGA_AMINOACIDS': [], 'Arcanobacterium_3PGA_AMINOACIDS': [], 'Pimelobacter_3PGA_AMINOACIDS': [], 'Corynebacteriales_3PGA_AMINOACIDS': [], 'Parascardovia_3PGA_AMINOACIDS': [], 'Thermomonospora_3PGA_AMINOACIDS': [], 'Luteipulveratus_3PGA_AMINOACIDS': [], 'Actinotignum_3PGA_AMINOACIDS': [], 'Alloactinosynnema_3PGA_AMINOACIDS': [], 'Amycolicicoccus_3PGA_AMINOACIDS': [], 'Actinobacteria_3PGA_AMINOACIDS': [], 'Scardovia_3PGA_AMINOACIDS': [], 'Adlercreutzia_3PGA_AMINOACIDS': [], 'Coriobacteriaceae_3PGA_AMINOACIDS': [], 'Gardnerella_3PGA_AMINOACIDS': [], 'Rubrobacter_3PGA_AMINOACIDS': [], 'Gordonibacter_3PGA_AMINOACIDS': [], 'Stackebrandtia_3PGA_AMINOACIDS': [], 'Thermobispora_3PGA_AMINOACIDS': [], 'Coriobacterium_3PGA_AMINOACIDS': [], 'Conexibacter_3PGA_AMINOACIDS': [], 'Ilumatobacter_3PGA_AMINOACIDS': [], 'Cellvibrio_3PGA_AMINOACIDS': [], 'Acidimicrobium_3PGA_AMINOACIDS': [], 'Actinokieospora_3PGA_AMINOACIDS': [], 'Rhodococcus_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Propionibacterium_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Nocardiopsis_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Kutzneria_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Citricoccus_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Salinispora_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Kocuria_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Actinosynnema_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Kitasatospora_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Micromonospora_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Actinoplanes_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Microlunatus_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Slackia_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Cryptobacterium_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Anaerolinea_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Chloroflexus_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Caldilinea_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Roseiflexus_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Dehalogenimonas_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Thermomicrobium_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Eggerthella_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Atopobium_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Dehalococcoides_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Brevibacterium_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Dermabacter_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Leifsonia_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Saccharopolyspora_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Leucobacter_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Rothia_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Agromyces_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Actinoalloteichus_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Sciscionella_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Actinocatenispora_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Lechevalieria_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Thermobifida_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Thermocrispum_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Actinokineospora_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Lentzea_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Kibdelosporangium_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Actinomycetospora_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Saccharothrix_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Collinsella_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Propionimicrobium_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Thermogemmatispora_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Granulicoccus_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Nitrolancea_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Enterorhabdus_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Olsenella_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Senegalimassilia_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Oscillochloris_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Enorma_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Aestuariimicrobium_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Oerskovia_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Cellulomonas_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Paraoerskovia_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Isoptericola_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Brachybacterium_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Promicromonospora_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Cellulosimicrobium_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Xylanimonas_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Kineococcus_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Kineosporia_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Dietzia_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Sphaerobacter_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Ktedonobacter_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Tomitella_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Segniliparus_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Catenulispora_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Tsukamurella_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Actinopolyspora_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Williamsia_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Actinospica_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Tropheryma_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Beutenbergia_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Kytococcus_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Nakamurella_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Blastococcus_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Candidatus_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Acidothermus_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Lysinimicrobium_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Geodermatophilaceae_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Glycomyces_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Cryptosporangium_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Ruania_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Modestobacter_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Dermatophilus_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Dermacoccus_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Mobilicoccus_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Demetria_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Sanguibacter_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Kineosphaera_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Timonella_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Geodermatophilus_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Haloglycomyces_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Sporichthya_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Jonesia_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Janibacter_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Gryllotalpicola_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Frigoribacterium_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Herbiconiux_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Humibacter_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Serinicoccus_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Rathayibacter_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Agreia_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Cryobacterium_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Intrasporangium_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Tetrasphaera_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Intrasporangiaceae_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Glaciibacter_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Cryocola_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Austwickia_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Arsenicicoccus_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Salinibacterium_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Gulosibacter_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Pseudoclavibacter_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Mycetocola_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Phycicoccus_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Ornithinimicrobium_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Terrabacter_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Terracoccus_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Alloscardovia_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Actinopolymorpha_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Kribbella_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Agrococcus_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Micrococcus_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Nesterenkonia_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Dactylosporangium_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Nocardioidaceae_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Aeromicrobium_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Renibacterium_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Verrucosispora_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Sinomonas_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Marmoricola_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Pilimelia_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Catelliglobosispora_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Acaricomes_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Catenuloplanes_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Jiangella_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Longispora_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Yaniella_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Microbispora_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Nonomuraea_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Streptosporangium_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Streptomonospora_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Microtetraspora_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Prauserella_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Allokutzneria_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Streptacidiphilus_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Mobiluncus_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Trueperella_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Varibaculum_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Actinobaculum_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Arcanobacterium_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Pimelobacter_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Corynebacteriales_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Parascardovia_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Thermomonospora_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Luteipulveratus_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Actinotignum_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Alloactinosynnema_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Amycolicicoccus_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Actinobacteria_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Scardovia_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Adlercreutzia_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Coriobacteriaceae_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Gardnerella_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Rubrobacter_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Gordonibacter_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Stackebrandtia_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Thermobispora_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Coriobacterium_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Conexibacter_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Ilumatobacter_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Cellvibrio_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Acidimicrobium_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Actinokieospora_ALPHAKETOGLUTARATE_AMINOACIDS': [], 'Propionibacterium_E4P_AMINO_ACIDS': [], 'Nocardiopsis_E4P_AMINO_ACIDS': [], 'Clavibacter_E4P_AMINO_ACIDS': [], 'Kutzneria_E4P_AMINO_ACIDS': [], 'Citricoccus_E4P_AMINO_ACIDS': [], 'Nocardioides_E4P_AMINO_ACIDS': [], 'Arthrobacter_E4P_AMINO_ACIDS': [], 'Bifidobacterium_E4P_AMINO_ACIDS': [], 'Kocuria_E4P_AMINO_ACIDS': [], 'Actinosynnema_E4P_AMINO_ACIDS': [], 'Kitasatospora_E4P_AMINO_ACIDS': [], 'Micromonospora_E4P_AMINO_ACIDS': [], 'Actinoplanes_E4P_AMINO_ACIDS': [], 'Microlunatus_E4P_AMINO_ACIDS': [], 'Slackia_E4P_AMINO_ACIDS': [], 'Cryptobacterium_E4P_AMINO_ACIDS': [], 'Anaerolinea_E4P_AMINO_ACIDS': [], 'Chloroflexus_E4P_AMINO_ACIDS': [], 'Caldilinea_E4P_AMINO_ACIDS': [], 'Roseiflexus_E4P_AMINO_ACIDS': [], 'Dehalogenimonas_E4P_AMINO_ACIDS': [], 'Thermomicrobium_E4P_AMINO_ACIDS': [], 'Eggerthella_E4P_AMINO_ACIDS': [], 'Atopobium_E4P_AMINO_ACIDS': [], 'Dehalococcoides_E4P_AMINO_ACIDS': [], 'Saccharomonospora_E4P_AMINO_ACIDS': [], 'Pseudonocardia_E4P_AMINO_ACIDS': [], 'Brevibacterium_E4P_AMINO_ACIDS': [], 'Dermabacter_E4P_AMINO_ACIDS': [], 'Leifsonia_E4P_AMINO_ACIDS': [], 'Saccharopolyspora_E4P_AMINO_ACIDS': [], 'Leucobacter_E4P_AMINO_ACIDS': [], 'Rothia_E4P_AMINO_ACIDS': [], 'Agromyces_E4P_AMINO_ACIDS': [], 'Actinoalloteichus_E4P_AMINO_ACIDS': [], 'Sciscionella_E4P_AMINO_ACIDS': [], 'Actinocatenispora_E4P_AMINO_ACIDS': [], 'Lechevalieria_E4P_AMINO_ACIDS': [], 'Thermobifida_E4P_AMINO_ACIDS': [], 'Thermocrispum_E4P_AMINO_ACIDS': [], 'Actinokineospora_E4P_AMINO_ACIDS': [], 'Lentzea_E4P_AMINO_ACIDS': [], 'Kibdelosporangium_E4P_AMINO_ACIDS': [], 'Actinomycetospora_E4P_AMINO_ACIDS': [], 'Saccharothrix_E4P_AMINO_ACIDS': [], 'Collinsella_E4P_AMINO_ACIDS': [], 'Propionimicrobium_E4P_AMINO_ACIDS': [], 'Thermogemmatispora_E4P_AMINO_ACIDS': [], 'Granulicoccus_E4P_AMINO_ACIDS': [], 'Nitrolancea_E4P_AMINO_ACIDS': [], 'Enterorhabdus_E4P_AMINO_ACIDS': [], 'Olsenella_E4P_AMINO_ACIDS': [], 'Senegalimassilia_E4P_AMINO_ACIDS': [], 'Oscillochloris_E4P_AMINO_ACIDS': [], 'Enorma_E4P_AMINO_ACIDS': [], 'Aestuariimicrobium_E4P_AMINO_ACIDS': [], 'Oerskovia_E4P_AMINO_ACIDS': [], 'Cellulomonas_E4P_AMINO_ACIDS': [], 'Paraoerskovia_E4P_AMINO_ACIDS': [], 'Isoptericola_E4P_AMINO_ACIDS': [], 'Brachybacterium_E4P_AMINO_ACIDS': [], 'Promicromonospora_E4P_AMINO_ACIDS': [], 'Cellulosimicrobium_E4P_AMINO_ACIDS': [], 'Xylanimonas_E4P_AMINO_ACIDS': [], 'Kineococcus_E4P_AMINO_ACIDS': [], 'Kineosporia_E4P_AMINO_ACIDS': [], 'Dietzia_E4P_AMINO_ACIDS': [], 'Sphaerobacter_E4P_AMINO_ACIDS': [], 'Ktedonobacter_E4P_AMINO_ACIDS': [], 'Tomitella_E4P_AMINO_ACIDS': [], 'Segniliparus_E4P_AMINO_ACIDS': [], 'Catenulispora_E4P_AMINO_ACIDS': [], 'Tsukamurella_E4P_AMINO_ACIDS': [], 'Actinopolyspora_E4P_AMINO_ACIDS': [], 'Williamsia_E4P_AMINO_ACIDS': [], 'Actinospica_E4P_AMINO_ACIDS': [], 'Tropheryma_E4P_AMINO_ACIDS': [], 'Beutenbergia_E4P_AMINO_ACIDS': [], 'Kytococcus_E4P_AMINO_ACIDS': [], 'Nakamurella_E4P_AMINO_ACIDS': [], 'Blastococcus_E4P_AMINO_ACIDS': [], 'Candidatus_E4P_AMINO_ACIDS': [], 'Acidothermus_E4P_AMINO_ACIDS': [], 'Lysinimicrobium_E4P_AMINO_ACIDS': [], 'Geodermatophilaceae_E4P_AMINO_ACIDS': [], 'Glycomyces_E4P_AMINO_ACIDS': [], 'Cryptosporangium_E4P_AMINO_ACIDS': [], 'Ruania_E4P_AMINO_ACIDS': [], 'Modestobacter_E4P_AMINO_ACIDS': [], 'Dermatophilus_E4P_AMINO_ACIDS': [], 'Dermacoccus_E4P_AMINO_ACIDS': [], 'Mobilicoccus_E4P_AMINO_ACIDS': [], 'Demetria_E4P_AMINO_ACIDS': [], 'Sanguibacter_E4P_AMINO_ACIDS': [], 'Kineosphaera_E4P_AMINO_ACIDS': [], 'Timonella_E4P_AMINO_ACIDS': [], 'Geodermatophilus_E4P_AMINO_ACIDS': [], 'Haloglycomyces_E4P_AMINO_ACIDS': [], 'Sporichthya_E4P_AMINO_ACIDS': [], 'Jonesia_E4P_AMINO_ACIDS': [], 'Janibacter_E4P_AMINO_ACIDS': [], 'Gryllotalpicola_E4P_AMINO_ACIDS': [], 'Frigoribacterium_E4P_AMINO_ACIDS': [], 'Herbiconiux_E4P_AMINO_ACIDS': [], 'Knoellia_E4P_AMINO_ACIDS': [], 'Humibacter_E4P_AMINO_ACIDS': [], 'Serinicoccus_E4P_AMINO_ACIDS': [], 'Rathayibacter_E4P_AMINO_ACIDS': [], 'Agreia_E4P_AMINO_ACIDS': [], 'Cryobacterium_E4P_AMINO_ACIDS': [], 'Intrasporangium_E4P_AMINO_ACIDS': [], 'Tetrasphaera_E4P_AMINO_ACIDS': [], 'Intrasporangiaceae_E4P_AMINO_ACIDS': [], 'Glaciibacter_E4P_AMINO_ACIDS': [], 'Cryocola_E4P_AMINO_ACIDS': [], 'Austwickia_E4P_AMINO_ACIDS': [], 'Arsenicicoccus_E4P_AMINO_ACIDS': [], 'Salinibacterium_E4P_AMINO_ACIDS': [], 'Gulosibacter_E4P_AMINO_ACIDS': [], 'Pseudoclavibacter_E4P_AMINO_ACIDS': [], 'Mycetocola_E4P_AMINO_ACIDS': [], 'Phycicoccus_E4P_AMINO_ACIDS': [], 'Ornithinimicrobium_E4P_AMINO_ACIDS': [], 'Terrabacter_E4P_AMINO_ACIDS': [], 'Terracoccus_E4P_AMINO_ACIDS': [], 'Alloscardovia_E4P_AMINO_ACIDS': [], 'Actinopolymorpha_E4P_AMINO_ACIDS': [], 'Kribbella_E4P_AMINO_ACIDS': [], 'Agrococcus_E4P_AMINO_ACIDS': [], 'Micrococcus_E4P_AMINO_ACIDS': [], 'Nesterenkonia_E4P_AMINO_ACIDS': [], 'Dactylosporangium_E4P_AMINO_ACIDS': [], 'Nocardioidaceae_E4P_AMINO_ACIDS': [], 'Aeromicrobium_E4P_AMINO_ACIDS': [], 'Renibacterium_E4P_AMINO_ACIDS': [], 'Verrucosispora_E4P_AMINO_ACIDS': [], 'Sinomonas_E4P_AMINO_ACIDS': [], 'Marmoricola_E4P_AMINO_ACIDS': [], 'Pilimelia_E4P_AMINO_ACIDS': [], 'Catelliglobosispora_E4P_AMINO_ACIDS': [], 'Acaricomes_E4P_AMINO_ACIDS': [], 'Catenuloplanes_E4P_AMINO_ACIDS': [], 'Jiangella_E4P_AMINO_ACIDS': [], 'Longispora_E4P_AMINO_ACIDS': [], 'Yaniella_E4P_AMINO_ACIDS': [], 'Actinomadura_E4P_AMINO_ACIDS': [], 'Microbispora_E4P_AMINO_ACIDS': [], 'Nonomuraea_E4P_AMINO_ACIDS': [], 'Streptosporangium_E4P_AMINO_ACIDS': [], 'Streptomonospora_E4P_AMINO_ACIDS': [], 'Microtetraspora_E4P_AMINO_ACIDS': [], 'Prauserella_E4P_AMINO_ACIDS': [], 'Allokutzneria_E4P_AMINO_ACIDS': [], 'Streptacidiphilus_E4P_AMINO_ACIDS': [], 'Mobiluncus_E4P_AMINO_ACIDS': [], 'Trueperella_E4P_AMINO_ACIDS': [], 'Varibaculum_E4P_AMINO_ACIDS': [], 'Actinobaculum_E4P_AMINO_ACIDS': [], 'Arcanobacterium_E4P_AMINO_ACIDS': [], 'Pimelobacter_E4P_AMINO_ACIDS': [], 'Corynebacteriales_E4P_AMINO_ACIDS': [], 'Parascardovia_E4P_AMINO_ACIDS': [], 'Thermomonospora_E4P_AMINO_ACIDS': [], 'Luteipulveratus_E4P_AMINO_ACIDS': [], 'Actinotignum_E4P_AMINO_ACIDS': [], 'Alloactinosynnema_E4P_AMINO_ACIDS': [], 'Amycolicicoccus_E4P_AMINO_ACIDS': [], 'Actinobacteria_E4P_AMINO_ACIDS': [], 'Scardovia_E4P_AMINO_ACIDS': [], 'Adlercreutzia_E4P_AMINO_ACIDS': [], 'Coriobacteriaceae_E4P_AMINO_ACIDS': [], 'Gardnerella_E4P_AMINO_ACIDS': [], 'Rubrobacter_E4P_AMINO_ACIDS': [], 'Gordonibacter_E4P_AMINO_ACIDS': [], 'Stackebrandtia_E4P_AMINO_ACIDS': [], 'Thermobispora_E4P_AMINO_ACIDS': [], 'Coriobacterium_E4P_AMINO_ACIDS': [], 'Conexibacter_E4P_AMINO_ACIDS': [], 'Ilumatobacter_E4P_AMINO_ACIDS': [], 'Cellvibrio_E4P_AMINO_ACIDS': [], 'Acidimicrobium_E4P_AMINO_ACIDS': [], 'Actinokieospora_E4P_AMINO_ACIDS': [], 'Nocardia_Glycolysis': [], 'Propionibacterium_Glycolysis': [], 'Nocardiopsis_Glycolysis': [], 'Clavibacter_Glycolysis': [], 'Kutzneria_Glycolysis': [], 'Citricoccus_Glycolysis': [], 'Nocardioides_Glycolysis': [], 'Arthrobacter_Glycolysis': [], 'Frankia_Glycolysis': [], 'Bifidobacterium_Glycolysis': [], 'Salinispora_Glycolysis': [], 'Kocuria_Glycolysis': [], 'Actinosynnema_Glycolysis': [], 'Kitasatospora_Glycolysis': [], 'Micromonospora_Glycolysis': [], 'Actinoplanes_Glycolysis': [], 'Microlunatus_Glycolysis': [], 'Slackia_Glycolysis': [], 'Cryptobacterium_Glycolysis': [], 'Anaerolinea_Glycolysis': [], 'Chloroflexus_Glycolysis': [], 'Caldilinea_Glycolysis': [], 'Roseiflexus_Glycolysis': [], 'Dehalogenimonas_Glycolysis': [], 'Thermomicrobium_Glycolysis': [], 'Eggerthella_Glycolysis': [], 'Atopobium_Glycolysis': [], 'Dehalococcoides_Glycolysis': [], 'Microbacterium_Glycolysis': [], 'Saccharomonospora_Glycolysis': [], 'Pseudonocardia_Glycolysis': [], 'Brevibacterium_Glycolysis': [], 'Dermabacter_Glycolysis': [], 'Leifsonia_Glycolysis': [], 'Saccharopolyspora_Glycolysis': [], 'Leucobacter_Glycolysis': [], 'Rothia_Glycolysis': [], 'Agromyces_Glycolysis': [], 'Actinoalloteichus_Glycolysis': [], 'Sciscionella_Glycolysis': [], 'Actinocatenispora_Glycolysis': [], 'Lechevalieria_Glycolysis': [], 'Thermobifida_Glycolysis': [], 'Thermocrispum_Glycolysis': [], 'Actinokineospora_Glycolysis': [], 'Lentzea_Glycolysis': [], 'Kibdelosporangium_Glycolysis': [], 'Actinomycetospora_Glycolysis': [], 'Saccharothrix_Glycolysis': [], 'Collinsella_Glycolysis': [], 'Propionimicrobium_Glycolysis': [], 'Thermogemmatispora_Glycolysis': [], 'Granulicoccus_Glycolysis': [], 'Nitrolancea_Glycolysis': [], 'Enterorhabdus_Glycolysis': [], 'Olsenella_Glycolysis': [], 'Senegalimassilia_Glycolysis': [], 'Oscillochloris_Glycolysis': [], 'Enorma_Glycolysis': [], 'Aestuariimicrobium_Glycolysis': [], 'Oerskovia_Glycolysis': [], 'Cellulomonas_Glycolysis': [], 'Paraoerskovia_Glycolysis': [], 'Isoptericola_Glycolysis': [], 'Brachybacterium_Glycolysis': [], 'Promicromonospora_Glycolysis': [], 'Cellulosimicrobium_Glycolysis': [], 'Xylanimonas_Glycolysis': [], 'Kineococcus_Glycolysis': [], 'Kineosporia_Glycolysis': [], 'Dietzia_Glycolysis': [], 'Sphaerobacter_Glycolysis': [], 'Ktedonobacter_Glycolysis': [], 'Tomitella_Glycolysis': [], 'Segniliparus_Glycolysis': [], 'Catenulispora_Glycolysis': [], 'Tsukamurella_Glycolysis': [], 'Actinopolyspora_Glycolysis': [], 'Williamsia_Glycolysis': [], 'Actinospica_Glycolysis': [], 'Tropheryma_Glycolysis': [], 'Beutenbergia_Glycolysis': [], 'Kytococcus_Glycolysis': [], 'Nakamurella_Glycolysis': [], 'Blastococcus_Glycolysis': [], 'Candidatus_Glycolysis': [], 'Acidothermus_Glycolysis': [], 'Lysinimicrobium_Glycolysis': [], 'Geodermatophilaceae_Glycolysis': [], 'Glycomyces_Glycolysis': [], 'Cryptosporangium_Glycolysis': [], 'Ruania_Glycolysis': [], 'Modestobacter_Glycolysis': [], 'Dermatophilus_Glycolysis': [], 'Dermacoccus_Glycolysis': [], 'Mobilicoccus_Glycolysis': [], 'Demetria_Glycolysis': [], 'Sanguibacter_Glycolysis': [], 'Kineosphaera_Glycolysis': [], 'Timonella_Glycolysis': [], 'Geodermatophilus_Glycolysis': [], 'Haloglycomyces_Glycolysis': [], 'Sporichthya_Glycolysis': [], 'Jonesia_Glycolysis': [], 'Janibacter_Glycolysis': [], 'Gryllotalpicola_Glycolysis': [], 'Frigoribacterium_Glycolysis': [], 'Herbiconiux_Glycolysis': [], 'Knoellia_Glycolysis': [], 'Humibacter_Glycolysis': [], 'Serinicoccus_Glycolysis': [], 'Rathayibacter_Glycolysis': [], 'Agreia_Glycolysis': [], 'Cryobacterium_Glycolysis': [], 'Intrasporangium_Glycolysis': [], 'Tetrasphaera_Glycolysis': [], 'Intrasporangiaceae_Glycolysis': [], 'Glaciibacter_Glycolysis': [], 'Cryocola_Glycolysis': [], 'Austwickia_Glycolysis': [], 'Arsenicicoccus_Glycolysis': [], 'Salinibacterium_Glycolysis': [], 'Gulosibacter_Glycolysis': [], 'Pseudoclavibacter_Glycolysis': [], 'Mycetocola_Glycolysis': [], 'Phycicoccus_Glycolysis': [], 'Ornithinimicrobium_Glycolysis': [], 'Terrabacter_Glycolysis': [], 'Terracoccus_Glycolysis': [], 'Alloscardovia_Glycolysis': [], 'Actinopolymorpha_Glycolysis': [], 'Kribbella_Glycolysis': [], 'Agrococcus_Glycolysis': [], 'Micrococcus_Glycolysis': [], 'Nesterenkonia_Glycolysis': [], 'Dactylosporangium_Glycolysis': [], 'Nocardioidaceae_Glycolysis': [], 'Aeromicrobium_Glycolysis': [], 'Renibacterium_Glycolysis': [], 'Verrucosispora_Glycolysis': [], 'Sinomonas_Glycolysis': [], 'Marmoricola_Glycolysis': [], 'Pilimelia_Glycolysis': [], 'Catelliglobosispora_Glycolysis': [], 'Acaricomes_Glycolysis': [], 'Catenuloplanes_Glycolysis': [], 'Jiangella_Glycolysis': [], 'Longispora_Glycolysis': [], 'Yaniella_Glycolysis': [], 'Actinomadura_Glycolysis': [], 'Microbispora_Glycolysis': [], 'Nonomuraea_Glycolysis': [], 'Streptosporangium_Glycolysis': [], 'Streptomonospora_Glycolysis': [], 'Microtetraspora_Glycolysis': [], 'Prauserella_Glycolysis': [], 'Allokutzneria_Glycolysis': [], 'Streptacidiphilus_Glycolysis': [], 'Mobiluncus_Glycolysis': [], 'Actinomyces_Glycolysis': [], 'Trueperella_Glycolysis': [], 'Varibaculum_Glycolysis': [], 'Actinobaculum_Glycolysis': [], 'Arcanobacterium_Glycolysis': [], 'Pimelobacter_Glycolysis': [], 'Corynebacteriales_Glycolysis': [], 'Parascardovia_Glycolysis': [], 'Thermomonospora_Glycolysis': [], 'Luteipulveratus_Glycolysis': [], 'Actinotignum_Glycolysis': [], 'Alloactinosynnema_Glycolysis': [], 'Amycolicicoccus_Glycolysis': [], 'Actinobacteria_Glycolysis': [], 'Scardovia_Glycolysis': [], 'Adlercreutzia_Glycolysis': [], 'Coriobacteriaceae_Glycolysis': [], 'Gardnerella_Glycolysis': [], 'Rubrobacter_Glycolysis': [], 'Gordonibacter_Glycolysis': [], 'Stackebrandtia_Glycolysis': [], 'Thermobispora_Glycolysis': [], 'Coriobacterium_Glycolysis': [], 'Conexibacter_Glycolysis': [], 'Ilumatobacter_Glycolysis': [], 'Cellvibrio_Glycolysis': [], 'Acidimicrobium_Glycolysis': [], 'Actinokieospora_Glycolysis': [], 'Propionibacterium_OXALACETATE_AMINOACIDS': [], 'Clavibacter_OXALACETATE_AMINOACIDS': [], 'Kutzneria_OXALACETATE_AMINOACIDS': [], 'Citricoccus_OXALACETATE_AMINOACIDS': [], 'Frankia_OXALACETATE_AMINOACIDS': [], 'Salinispora_OXALACETATE_AMINOACIDS': [], 'Kocuria_OXALACETATE_AMINOACIDS': [], 'Actinosynnema_OXALACETATE_AMINOACIDS': [], 'Kitasatospora_OXALACETATE_AMINOACIDS': [], 'Actinoplanes_OXALACETATE_AMINOACIDS': [], 'Microlunatus_OXALACETATE_AMINOACIDS': [], 'Slackia_OXALACETATE_AMINOACIDS': [], 'Cryptobacterium_OXALACETATE_AMINOACIDS': [], 'Anaerolinea_OXALACETATE_AMINOACIDS': [], 'Chloroflexus_OXALACETATE_AMINOACIDS': [], 'Caldilinea_OXALACETATE_AMINOACIDS': [], 'Roseiflexus_OXALACETATE_AMINOACIDS': [], 'Dehalogenimonas_OXALACETATE_AMINOACIDS': [], 'Thermomicrobium_OXALACETATE_AMINOACIDS': [], 'Eggerthella_OXALACETATE_AMINOACIDS': [], 'Atopobium_OXALACETATE_AMINOACIDS': [], 'Dehalococcoides_OXALACETATE_AMINOACIDS': [], 'Saccharomonospora_OXALACETATE_AMINOACIDS': [], 'Pseudonocardia_OXALACETATE_AMINOACIDS': [], 'Dermabacter_OXALACETATE_AMINOACIDS': [], 'Leifsonia_OXALACETATE_AMINOACIDS': [], 'Saccharopolyspora_OXALACETATE_AMINOACIDS': [], 'Leucobacter_OXALACETATE_AMINOACIDS': [], 'Rothia_OXALACETATE_AMINOACIDS': [], 'Agromyces_OXALACETATE_AMINOACIDS': [], 'Actinoalloteichus_OXALACETATE_AMINOACIDS': [], 'Sciscionella_OXALACETATE_AMINOACIDS': [], 'Actinocatenispora_OXALACETATE_AMINOACIDS': [], 'Lechevalieria_OXALACETATE_AMINOACIDS': [], 'Thermobifida_OXALACETATE_AMINOACIDS': [], 'Thermocrispum_OXALACETATE_AMINOACIDS': [], 'Actinokineospora_OXALACETATE_AMINOACIDS': [], 'Lentzea_OXALACETATE_AMINOACIDS': [], 'Kibdelosporangium_OXALACETATE_AMINOACIDS': [], 'Actinomycetospora_OXALACETATE_AMINOACIDS': [], 'Saccharothrix_OXALACETATE_AMINOACIDS': [], 'Collinsella_OXALACETATE_AMINOACIDS': [], 'Propionimicrobium_OXALACETATE_AMINOACIDS': [], 'Thermogemmatispora_OXALACETATE_AMINOACIDS': [], 'Granulicoccus_OXALACETATE_AMINOACIDS': [], 'Nitrolancea_OXALACETATE_AMINOACIDS': [], 'Enterorhabdus_OXALACETATE_AMINOACIDS': [], 'Olsenella_OXALACETATE_AMINOACIDS': [], 'Senegalimassilia_OXALACETATE_AMINOACIDS': [], 'Oscillochloris_OXALACETATE_AMINOACIDS': [], 'Enorma_OXALACETATE_AMINOACIDS': [], 'Aestuariimicrobium_OXALACETATE_AMINOACIDS': [], 'Oerskovia_OXALACETATE_AMINOACIDS': [], 'Cellulomonas_OXALACETATE_AMINOACIDS': [], 'Paraoerskovia_OXALACETATE_AMINOACIDS': [], 'Isoptericola_OXALACETATE_AMINOACIDS': [], 'Brachybacterium_OXALACETATE_AMINOACIDS': [], 'Promicromonospora_OXALACETATE_AMINOACIDS': [], 'Cellulosimicrobium_OXALACETATE_AMINOACIDS': [], 'Xylanimonas_OXALACETATE_AMINOACIDS': [], 'Kineococcus_OXALACETATE_AMINOACIDS': [], 'Kineosporia_OXALACETATE_AMINOACIDS': [], 'Dietzia_OXALACETATE_AMINOACIDS': [], 'Sphaerobacter_OXALACETATE_AMINOACIDS': [], 'Ktedonobacter_OXALACETATE_AMINOACIDS': [], 'Tomitella_OXALACETATE_AMINOACIDS': [], 'Segniliparus_OXALACETATE_AMINOACIDS': [], 'Catenulispora_OXALACETATE_AMINOACIDS': [], 'Tsukamurella_OXALACETATE_AMINOACIDS': [], 'Actinopolyspora_OXALACETATE_AMINOACIDS': [], 'Williamsia_OXALACETATE_AMINOACIDS': [], 'Actinospica_OXALACETATE_AMINOACIDS': [], 'Tropheryma_OXALACETATE_AMINOACIDS': [], 'Beutenbergia_OXALACETATE_AMINOACIDS': [], 'Kytococcus_OXALACETATE_AMINOACIDS': [], 'Nakamurella_OXALACETATE_AMINOACIDS': [], 'Blastococcus_OXALACETATE_AMINOACIDS': [], 'Candidatus_OXALACETATE_AMINOACIDS': [], 'Acidothermus_OXALACETATE_AMINOACIDS': [], 'Lysinimicrobium_OXALACETATE_AMINOACIDS': [], 'Geodermatophilaceae_OXALACETATE_AMINOACIDS': [], 'Glycomyces_OXALACETATE_AMINOACIDS': [], 'Cryptosporangium_OXALACETATE_AMINOACIDS': [], 'Ruania_OXALACETATE_AMINOACIDS': [], 'Modestobacter_OXALACETATE_AMINOACIDS': [], 'Dermatophilus_OXALACETATE_AMINOACIDS': [], 'Dermacoccus_OXALACETATE_AMINOACIDS': [], 'Mobilicoccus_OXALACETATE_AMINOACIDS': [], 'Demetria_OXALACETATE_AMINOACIDS': [], 'Sanguibacter_OXALACETATE_AMINOACIDS': [], 'Kineosphaera_OXALACETATE_AMINOACIDS': [], 'Timonella_OXALACETATE_AMINOACIDS': [], 'Geodermatophilus_OXALACETATE_AMINOACIDS': [], 'Haloglycomyces_OXALACETATE_AMINOACIDS': [], 'Sporichthya_OXALACETATE_AMINOACIDS': [], 'Jonesia_OXALACETATE_AMINOACIDS': [], 'Janibacter_OXALACETATE_AMINOACIDS': [], 'Gryllotalpicola_OXALACETATE_AMINOACIDS': [], 'Frigoribacterium_OXALACETATE_AMINOACIDS': [], 'Herbiconiux_OXALACETATE_AMINOACIDS': [], 'Knoellia_OXALACETATE_AMINOACIDS': [], 'Humibacter_OXALACETATE_AMINOACIDS': [], 'Serinicoccus_OXALACETATE_AMINOACIDS': [], 'Rathayibacter_OXALACETATE_AMINOACIDS': [], 'Agreia_OXALACETATE_AMINOACIDS': [], 'Cryobacterium_OXALACETATE_AMINOACIDS': [], 'Intrasporangium_OXALACETATE_AMINOACIDS': [], 'Tetrasphaera_OXALACETATE_AMINOACIDS': [], 'Intrasporangiaceae_OXALACETATE_AMINOACIDS': [], 'Glaciibacter_OXALACETATE_AMINOACIDS': [], 'Cryocola_OXALACETATE_AMINOACIDS': [], 'Austwickia_OXALACETATE_AMINOACIDS': [], 'Arsenicicoccus_OXALACETATE_AMINOACIDS': [], 'Salinibacterium_OXALACETATE_AMINOACIDS': [], 'Gulosibacter_OXALACETATE_AMINOACIDS': [], 'Pseudoclavibacter_OXALACETATE_AMINOACIDS': [], 'Mycetocola_OXALACETATE_AMINOACIDS': [], 'Phycicoccus_OXALACETATE_AMINOACIDS': [], 'Ornithinimicrobium_OXALACETATE_AMINOACIDS': [], 'Terrabacter_OXALACETATE_AMINOACIDS': [], 'Terracoccus_OXALACETATE_AMINOACIDS': [], 'Alloscardovia_OXALACETATE_AMINOACIDS': [], 'Actinopolymorpha_OXALACETATE_AMINOACIDS': [], 'Kribbella_OXALACETATE_AMINOACIDS': [], 'Agrococcus_OXALACETATE_AMINOACIDS': [], 'Micrococcus_OXALACETATE_AMINOACIDS': [], 'Nesterenkonia_OXALACETATE_AMINOACIDS': [], 'Dactylosporangium_OXALACETATE_AMINOACIDS': [], 'Nocardioidaceae_OXALACETATE_AMINOACIDS': [], 'Aeromicrobium_OXALACETATE_AMINOACIDS': [], 'Renibacterium_OXALACETATE_AMINOACIDS': [], 'Verrucosispora_OXALACETATE_AMINOACIDS': [], 'Sinomonas_OXALACETATE_AMINOACIDS': [], 'Marmoricola_OXALACETATE_AMINOACIDS': [], 'Pilimelia_OXALACETATE_AMINOACIDS': [], 'Catelliglobosispora_OXALACETATE_AMINOACIDS': [], 'Acaricomes_OXALACETATE_AMINOACIDS': [], 'Catenuloplanes_OXALACETATE_AMINOACIDS': [], 'Jiangella_OXALACETATE_AMINOACIDS': [], 'Longispora_OXALACETATE_AMINOACIDS': [], 'Yaniella_OXALACETATE_AMINOACIDS': [], 'Actinomadura_OXALACETATE_AMINOACIDS': [], 'Microbispora_OXALACETATE_AMINOACIDS': [], 'Nonomuraea_OXALACETATE_AMINOACIDS': [], 'Streptosporangium_OXALACETATE_AMINOACIDS': [], 'Streptomonospora_OXALACETATE_AMINOACIDS': [], 'Microtetraspora_OXALACETATE_AMINOACIDS': [], 'Prauserella_OXALACETATE_AMINOACIDS': [], 'Allokutzneria_OXALACETATE_AMINOACIDS': [], 'Streptacidiphilus_OXALACETATE_AMINOACIDS': [], 'Mobiluncus_OXALACETATE_AMINOACIDS': [], 'Trueperella_OXALACETATE_AMINOACIDS': [], 'Varibaculum_OXALACETATE_AMINOACIDS': [], 'Actinobaculum_OXALACETATE_AMINOACIDS': [], 'Arcanobacterium_OXALACETATE_AMINOACIDS': [], 'Pimelobacter_OXALACETATE_AMINOACIDS': [], 'Corynebacteriales_OXALACETATE_AMINOACIDS': [], 'Parascardovia_OXALACETATE_AMINOACIDS': [], 'Thermomonospora_OXALACETATE_AMINOACIDS': [], 'Luteipulveratus_OXALACETATE_AMINOACIDS': [], 'Actinotignum_OXALACETATE_AMINOACIDS': [], 'Alloactinosynnema_OXALACETATE_AMINOACIDS': [], 'Amycolicicoccus_OXALACETATE_AMINOACIDS': [], 'Actinobacteria_OXALACETATE_AMINOACIDS': [], 'Scardovia_OXALACETATE_AMINOACIDS': [], 'Adlercreutzia_OXALACETATE_AMINOACIDS': [], 'Coriobacteriaceae_OXALACETATE_AMINOACIDS': [], 'Gardnerella_OXALACETATE_AMINOACIDS': [], 'Rubrobacter_OXALACETATE_AMINOACIDS': [], 'Gordonibacter_OXALACETATE_AMINOACIDS': [], 'Stackebrandtia_OXALACETATE_AMINOACIDS': [], 'Thermobispora_OXALACETATE_AMINOACIDS': [], 'Coriobacterium_OXALACETATE_AMINOACIDS': [], 'Conexibacter_OXALACETATE_AMINOACIDS': [], 'Ilumatobacter_OXALACETATE_AMINOACIDS': [], 'Cellvibrio_OXALACETATE_AMINOACIDS': [], 'Acidimicrobium_OXALACETATE_AMINOACIDS': [], 'Actinokieospora_OXALACETATE_AMINOACIDS': [], 'Propionibacterium_PYR_THR_AA': [], 'Clavibacter_PYR_THR_AA': [], 'Kutzneria_PYR_THR_AA': [], 'Citricoccus_PYR_THR_AA': [], 'Frankia_PYR_THR_AA': [], 'Salinispora_PYR_THR_AA': [], 'Kocuria_PYR_THR_AA': [], 'Actinosynnema_PYR_THR_AA': [], 'Micromonospora_PYR_THR_AA': [], 'Microlunatus_PYR_THR_AA': [], 'Slackia_PYR_THR_AA': [], 'Cryptobacterium_PYR_THR_AA': [], 'Anaerolinea_PYR_THR_AA': [], 'Chloroflexus_PYR_THR_AA': [], 'Caldilinea_PYR_THR_AA': [], 'Roseiflexus_PYR_THR_AA': [], 'Dehalogenimonas_PYR_THR_AA': [], 'Thermomicrobium_PYR_THR_AA': [], 'Eggerthella_PYR_THR_AA': [], 'Atopobium_PYR_THR_AA': [], 'Dehalococcoides_PYR_THR_AA': [], 'Dermabacter_PYR_THR_AA': [], 'Leifsonia_PYR_THR_AA': [], 'Saccharopolyspora_PYR_THR_AA': [], 'Leucobacter_PYR_THR_AA': [], 'Rothia_PYR_THR_AA': [], 'Agromyces_PYR_THR_AA': [], 'Actinoalloteichus_PYR_THR_AA': [], 'Sciscionella_PYR_THR_AA': [], 'Actinocatenispora_PYR_THR_AA': [], 'Lechevalieria_PYR_THR_AA': [], 'Thermobifida_PYR_THR_AA': [], 'Thermocrispum_PYR_THR_AA': [], 'Actinokineospora_PYR_THR_AA': [], 'Lentzea_PYR_THR_AA': [], 'Kibdelosporangium_PYR_THR_AA': [], 'Actinomycetospora_PYR_THR_AA': [], 'Saccharothrix_PYR_THR_AA': [], 'Collinsella_PYR_THR_AA': [], 'Propionimicrobium_PYR_THR_AA': [], 'Thermogemmatispora_PYR_THR_AA': [], 'Granulicoccus_PYR_THR_AA': [], 'Nitrolancea_PYR_THR_AA': [], 'Enterorhabdus_PYR_THR_AA': [], 'Olsenella_PYR_THR_AA': [], 'Senegalimassilia_PYR_THR_AA': [], 'Oscillochloris_PYR_THR_AA': [], 'Enorma_PYR_THR_AA': [], 'Aestuariimicrobium_PYR_THR_AA': [], 'Oerskovia_PYR_THR_AA': [], 'Cellulomonas_PYR_THR_AA': [], 'Paraoerskovia_PYR_THR_AA': [], 'Isoptericola_PYR_THR_AA': [], 'Brachybacterium_PYR_THR_AA': [], 'Promicromonospora_PYR_THR_AA': [], 'Cellulosimicrobium_PYR_THR_AA': [], 'Xylanimonas_PYR_THR_AA': [], 'Kineococcus_PYR_THR_AA': [], 'Kineosporia_PYR_THR_AA': [], 'Dietzia_PYR_THR_AA': [], 'Sphaerobacter_PYR_THR_AA': [], 'Ktedonobacter_PYR_THR_AA': [], 'Tomitella_PYR_THR_AA': [], 'Segniliparus_PYR_THR_AA': [], 'Catenulispora_PYR_THR_AA': [], 'Tsukamurella_PYR_THR_AA': [], 'Actinopolyspora_PYR_THR_AA': [], 'Williamsia_PYR_THR_AA': [], 'Actinospica_PYR_THR_AA': [], 'Tropheryma_PYR_THR_AA': [], 'Beutenbergia_PYR_THR_AA': [], 'Kytococcus_PYR_THR_AA': [], 'Nakamurella_PYR_THR_AA': [], 'Blastococcus_PYR_THR_AA': [], 'Candidatus_PYR_THR_AA': [], 'Acidothermus_PYR_THR_AA': [], 'Lysinimicrobium_PYR_THR_AA': [], 'Geodermatophilaceae_PYR_THR_AA': [], 'Glycomyces_PYR_THR_AA': [], 'Cryptosporangium_PYR_THR_AA': [], 'Ruania_PYR_THR_AA': [], 'Modestobacter_PYR_THR_AA': [], 'Dermatophilus_PYR_THR_AA': [], 'Dermacoccus_PYR_THR_AA': [], 'Mobilicoccus_PYR_THR_AA': [], 'Demetria_PYR_THR_AA': [], 'Sanguibacter_PYR_THR_AA': [], 'Kineosphaera_PYR_THR_AA': [], 'Timonella_PYR_THR_AA': [], 'Geodermatophilus_PYR_THR_AA': [], 'Haloglycomyces_PYR_THR_AA': [], 'Sporichthya_PYR_THR_AA': [], 'Jonesia_PYR_THR_AA': [], 'Janibacter_PYR_THR_AA': [], 'Gryllotalpicola_PYR_THR_AA': [], 'Frigoribacterium_PYR_THR_AA': [], 'Herbiconiux_PYR_THR_AA': [], 'Knoellia_PYR_THR_AA': [], 'Humibacter_PYR_THR_AA': [], 'Serinicoccus_PYR_THR_AA': [], 'Rathayibacter_PYR_THR_AA': [], 'Agreia_PYR_THR_AA': [], 'Cryobacterium_PYR_THR_AA': [], 'Intrasporangium_PYR_THR_AA': [], 'Tetrasphaera_PYR_THR_AA': [], 'Intrasporangiaceae_PYR_THR_AA': [], 'Glaciibacter_PYR_THR_AA': [], 'Cryocola_PYR_THR_AA': [], 'Austwickia_PYR_THR_AA': [], 'Arsenicicoccus_PYR_THR_AA': [], 'Salinibacterium_PYR_THR_AA': [], 'Gulosibacter_PYR_THR_AA': [], 'Pseudoclavibacter_PYR_THR_AA': [], 'Mycetocola_PYR_THR_AA': [], 'Phycicoccus_PYR_THR_AA': [], 'Ornithinimicrobium_PYR_THR_AA': [], 'Terrabacter_PYR_THR_AA': [], 'Terracoccus_PYR_THR_AA': [], 'Alloscardovia_PYR_THR_AA': [], 'Actinopolymorpha_PYR_THR_AA': [], 'Kribbella_PYR_THR_AA': [], 'Agrococcus_PYR_THR_AA': [], 'Micrococcus_PYR_THR_AA': [], 'Nesterenkonia_PYR_THR_AA': [], 'Dactylosporangium_PYR_THR_AA': [], 'Nocardioidaceae_PYR_THR_AA': [], 'Aeromicrobium_PYR_THR_AA': [], 'Renibacterium_PYR_THR_AA': [], 'Verrucosispora_PYR_THR_AA': [], 'Sinomonas_PYR_THR_AA': [], 'Marmoricola_PYR_THR_AA': [], 'Pilimelia_PYR_THR_AA': [], 'Catelliglobosispora_PYR_THR_AA': [], 'Acaricomes_PYR_THR_AA': [], 'Catenuloplanes_PYR_THR_AA': [], 'Jiangella_PYR_THR_AA': [], 'Longispora_PYR_THR_AA': [], 'Yaniella_PYR_THR_AA': [], 'Actinomadura_PYR_THR_AA': [], 'Microbispora_PYR_THR_AA': [], 'Nonomuraea_PYR_THR_AA': [], 'Streptosporangium_PYR_THR_AA': [], 'Streptomonospora_PYR_THR_AA': [], 'Microtetraspora_PYR_THR_AA': [], 'Prauserella_PYR_THR_AA': [], 'Allokutzneria_PYR_THR_AA': [], 'Streptacidiphilus_PYR_THR_AA': [], 'Mobiluncus_PYR_THR_AA': [], 'Trueperella_PYR_THR_AA': [], 'Varibaculum_PYR_THR_AA': [], 'Actinobaculum_PYR_THR_AA': [], 'Arcanobacterium_PYR_THR_AA': [], 'Pimelobacter_PYR_THR_AA': [], 'Corynebacteriales_PYR_THR_AA': [], 'Parascardovia_PYR_THR_AA': [], 'Thermomonospora_PYR_THR_AA': [], 'Luteipulveratus_PYR_THR_AA': [], 'Actinotignum_PYR_THR_AA': [], 'Alloactinosynnema_PYR_THR_AA': [], 'Amycolicicoccus_PYR_THR_AA': [], 'Actinobacteria_PYR_THR_AA': [], 'Scardovia_PYR_THR_AA': [], 'Adlercreutzia_PYR_THR_AA': [], 'Coriobacteriaceae_PYR_THR_AA': [], 'Gardnerella_PYR_THR_AA': [], 'Rubrobacter_PYR_THR_AA': [], 'Gordonibacter_PYR_THR_AA': [], 'Stackebrandtia_PYR_THR_AA': [], 'Thermobispora_PYR_THR_AA': [], 'Coriobacterium_PYR_THR_AA': [], 'Conexibacter_PYR_THR_AA': [], 'Ilumatobacter_PYR_THR_AA': [], 'Cellvibrio_PYR_THR_AA': [], 'Acidimicrobium_PYR_THR_AA': [], 'Actinokieospora_PYR_THR_AA': [], 'Rhodococcus_R5P_AMINOACIDS': [], 'Amycolatopsis_R5P_AMINOACIDS': [], 'Propionibacterium_R5P_AMINOACIDS': [], 'Nocardiopsis_R5P_AMINOACIDS': [], 'Clavibacter_R5P_AMINOACIDS': [], 'Mycobacterium_R5P_AMINOACIDS': [], 'Kutzneria_R5P_AMINOACIDS': [], 'Citricoccus_R5P_AMINOACIDS': [], 'Nocardioides_R5P_AMINOACIDS': [], 'Corynebacterium_R5P_AMINOACIDS': [], 'Arthrobacter_R5P_AMINOACIDS': [], 'Frankia_R5P_AMINOACIDS': [], 'Bifidobacterium_R5P_AMINOACIDS': [], 'Salinispora_R5P_AMINOACIDS': [], 'Kocuria_R5P_AMINOACIDS': [], 'Actinosynnema_R5P_AMINOACIDS': [], 'Kitasatospora_R5P_AMINOACIDS': [], 'Gordonia_R5P_AMINOACIDS': [], 'Micromonospora_R5P_AMINOACIDS': [], 'Actinoplanes_R5P_AMINOACIDS': [], 'Microlunatus_R5P_AMINOACIDS': [], 'Slackia_R5P_AMINOACIDS': [], 'Cryptobacterium_R5P_AMINOACIDS': [], 'Anaerolinea_R5P_AMINOACIDS': [], 'Chloroflexus_R5P_AMINOACIDS': [], 'Caldilinea_R5P_AMINOACIDS': [], 'Roseiflexus_R5P_AMINOACIDS': [], 'Dehalogenimonas_R5P_AMINOACIDS': [], 'Thermomicrobium_R5P_AMINOACIDS': [], 'Eggerthella_R5P_AMINOACIDS': [], 'Atopobium_R5P_AMINOACIDS': [], 'Dehalococcoides_R5P_AMINOACIDS': [], 'Microbacterium_R5P_AMINOACIDS': [], 'Saccharomonospora_R5P_AMINOACIDS': [], 'Pseudonocardia_R5P_AMINOACIDS': [], 'Brevibacterium_R5P_AMINOACIDS': [], 'Dermabacter_R5P_AMINOACIDS': [], 'Leifsonia_R5P_AMINOACIDS': [], 'Saccharopolyspora_R5P_AMINOACIDS': [], 'Leucobacter_R5P_AMINOACIDS': [], 'Rothia_R5P_AMINOACIDS': [], 'Agromyces_R5P_AMINOACIDS': [], 'Actinoalloteichus_R5P_AMINOACIDS': [], 'Sciscionella_R5P_AMINOACIDS': [], 'Actinocatenispora_R5P_AMINOACIDS': [], 'Lechevalieria_R5P_AMINOACIDS': [], 'Thermobifida_R5P_AMINOACIDS': [], 'Thermocrispum_R5P_AMINOACIDS': [], 'Actinokineospora_R5P_AMINOACIDS': [], 'Lentzea_R5P_AMINOACIDS': [], 'Kibdelosporangium_R5P_AMINOACIDS': [], 'Actinomycetospora_R5P_AMINOACIDS': [], 'Saccharothrix_R5P_AMINOACIDS': [], 'Collinsella_R5P_AMINOACIDS': [], 'Propionimicrobium_R5P_AMINOACIDS': [], 'Thermogemmatispora_R5P_AMINOACIDS': [], 'Granulicoccus_R5P_AMINOACIDS': [], 'Nitrolancea_R5P_AMINOACIDS': [], 'Enterorhabdus_R5P_AMINOACIDS': [], 'Olsenella_R5P_AMINOACIDS': [], 'Senegalimassilia_R5P_AMINOACIDS': [], 'Oscillochloris_R5P_AMINOACIDS': [], 'Enorma_R5P_AMINOACIDS': [], 'Aestuariimicrobium_R5P_AMINOACIDS': [], 'Oerskovia_R5P_AMINOACIDS': [], 'Cellulomonas_R5P_AMINOACIDS': [], 'Paraoerskovia_R5P_AMINOACIDS': [], 'Isoptericola_R5P_AMINOACIDS': [], 'Brachybacterium_R5P_AMINOACIDS': [], 'Promicromonospora_R5P_AMINOACIDS': [], 'Cellulosimicrobium_R5P_AMINOACIDS': [], 'Xylanimonas_R5P_AMINOACIDS': [], 'Kineococcus_R5P_AMINOACIDS': [], 'Kineosporia_R5P_AMINOACIDS': [], 'Dietzia_R5P_AMINOACIDS': [], 'Sphaerobacter_R5P_AMINOACIDS': [], 'Ktedonobacter_R5P_AMINOACIDS': [], 'Tomitella_R5P_AMINOACIDS': [], 'Segniliparus_R5P_AMINOACIDS': [], 'Catenulispora_R5P_AMINOACIDS': [], 'Tsukamurella_R5P_AMINOACIDS': [], 'Actinopolyspora_R5P_AMINOACIDS': [], 'Williamsia_R5P_AMINOACIDS': [], 'Actinospica_R5P_AMINOACIDS': [], 'Tropheryma_R5P_AMINOACIDS': [], 'Beutenbergia_R5P_AMINOACIDS': [], 'Kytococcus_R5P_AMINOACIDS': [], 'Nakamurella_R5P_AMINOACIDS': [], 'Blastococcus_R5P_AMINOACIDS': [], 'Candidatus_R5P_AMINOACIDS': [], 'Acidothermus_R5P_AMINOACIDS': [], 'Lysinimicrobium_R5P_AMINOACIDS': [], 'Geodermatophilaceae_R5P_AMINOACIDS': [], 'Glycomyces_R5P_AMINOACIDS': [], 'Cryptosporangium_R5P_AMINOACIDS': [], 'Ruania_R5P_AMINOACIDS': [], 'Modestobacter_R5P_AMINOACIDS': [], 'Dermatophilus_R5P_AMINOACIDS': [], 'Dermacoccus_R5P_AMINOACIDS': [], 'Mobilicoccus_R5P_AMINOACIDS': [], 'Demetria_R5P_AMINOACIDS': [], 'Sanguibacter_R5P_AMINOACIDS': [], 'Kineosphaera_R5P_AMINOACIDS': [], 'Timonella_R5P_AMINOACIDS': [], 'Geodermatophilus_R5P_AMINOACIDS': [], 'Haloglycomyces_R5P_AMINOACIDS': [], 'Sporichthya_R5P_AMINOACIDS': [], 'Jonesia_R5P_AMINOACIDS': [], 'Janibacter_R5P_AMINOACIDS': [], 'Gryllotalpicola_R5P_AMINOACIDS': [], 'Frigoribacterium_R5P_AMINOACIDS': [], 'Herbiconiux_R5P_AMINOACIDS': [], 'Knoellia_R5P_AMINOACIDS': [], 'Humibacter_R5P_AMINOACIDS': [], 'Serinicoccus_R5P_AMINOACIDS': [], 'Rathayibacter_R5P_AMINOACIDS': [], 'Agreia_R5P_AMINOACIDS': [], 'Cryobacterium_R5P_AMINOACIDS': [], 'Intrasporangium_R5P_AMINOACIDS': [], 'Tetrasphaera_R5P_AMINOACIDS': [], 'Intrasporangiaceae_R5P_AMINOACIDS': [], 'Glaciibacter_R5P_AMINOACIDS': [], 'Cryocola_R5P_AMINOACIDS': [], 'Austwickia_R5P_AMINOACIDS': [], 'Arsenicicoccus_R5P_AMINOACIDS': [], 'Salinibacterium_R5P_AMINOACIDS': [], 'Gulosibacter_R5P_AMINOACIDS': [], 'Pseudoclavibacter_R5P_AMINOACIDS': [], 'Mycetocola_R5P_AMINOACIDS': [], 'Phycicoccus_R5P_AMINOACIDS': [], 'Ornithinimicrobium_R5P_AMINOACIDS': [], 'Terrabacter_R5P_AMINOACIDS': [], 'Terracoccus_R5P_AMINOACIDS': [], 'Alloscardovia_R5P_AMINOACIDS': [], 'Actinopolymorpha_R5P_AMINOACIDS': [], 'Kribbella_R5P_AMINOACIDS': [], 'Agrococcus_R5P_AMINOACIDS': [], 'Micrococcus_R5P_AMINOACIDS': [], 'Nesterenkonia_R5P_AMINOACIDS': [], 'Dactylosporangium_R5P_AMINOACIDS': [], 'Nocardioidaceae_R5P_AMINOACIDS': [], 'Aeromicrobium_R5P_AMINOACIDS': [], 'Renibacterium_R5P_AMINOACIDS': [], 'Verrucosispora_R5P_AMINOACIDS': [], 'Sinomonas_R5P_AMINOACIDS': [], 'Marmoricola_R5P_AMINOACIDS': [], 'Pilimelia_R5P_AMINOACIDS': [], 'Catelliglobosispora_R5P_AMINOACIDS': [], 'Acaricomes_R5P_AMINOACIDS': [], 'Catenuloplanes_R5P_AMINOACIDS': [], 'Jiangella_R5P_AMINOACIDS': [], 'Longispora_R5P_AMINOACIDS': [], 'Yaniella_R5P_AMINOACIDS': [], 'Actinomadura_R5P_AMINOACIDS': [], 'Microbispora_R5P_AMINOACIDS': [], 'Nonomuraea_R5P_AMINOACIDS': [], 'Streptosporangium_R5P_AMINOACIDS': [], 'Streptomonospora_R5P_AMINOACIDS': [], 'Microtetraspora_R5P_AMINOACIDS': [], 'Prauserella_R5P_AMINOACIDS': [], 'Allokutzneria_R5P_AMINOACIDS': [], 'Streptacidiphilus_R5P_AMINOACIDS': [], 'Mobiluncus_R5P_AMINOACIDS': [], 'Actinomyces_R5P_AMINOACIDS': [], 'Trueperella_R5P_AMINOACIDS': [], 'Varibaculum_R5P_AMINOACIDS': [], 'Actinobaculum_R5P_AMINOACIDS': [], 'Arcanobacterium_R5P_AMINOACIDS': [], 'Pimelobacter_R5P_AMINOACIDS': [], 'Corynebacteriales_R5P_AMINOACIDS': [], 'Parascardovia_R5P_AMINOACIDS': [], 'Thermomonospora_R5P_AMINOACIDS': [], 'Luteipulveratus_R5P_AMINOACIDS': [], 'Actinotignum_R5P_AMINOACIDS': [], 'Alloactinosynnema_R5P_AMINOACIDS': [], 'Amycolicicoccus_R5P_AMINOACIDS': [], 'Actinobacteria_R5P_AMINOACIDS': [], 'Scardovia_R5P_AMINOACIDS': [], 'Adlercreutzia_R5P_AMINOACIDS': [], 'Coriobacteriaceae_R5P_AMINOACIDS': [], 'Gardnerella_R5P_AMINOACIDS': [], 'Rubrobacter_R5P_AMINOACIDS': [], 'Gordonibacter_R5P_AMINOACIDS': [], 'Stackebrandtia_R5P_AMINOACIDS': [], 'Thermobispora_R5P_AMINOACIDS': [], 'Coriobacterium_R5P_AMINOACIDS': [], 'Conexibacter_R5P_AMINOACIDS': [], 'Ilumatobacter_R5P_AMINOACIDS': [], 'Cellvibrio_R5P_AMINOACIDS': [], 'Acidimicrobium_R5P_AMINOACIDS': [], 'Actinokieospora_R5P_AMINOACIDS': [], 'Rhodococcus_TCA': [], 'Nocardiopsis_TCA': [], 'Clavibacter_TCA': [], 'Kutzneria_TCA': [], 'Citricoccus_TCA': [], 'Nocardioides_TCA': [], 'Arthrobacter_TCA': [], 'Frankia_TCA': [], 'Bifidobacterium_TCA': [], 'Salinispora_TCA': [], 'Kocuria_TCA': [], 'Actinosynnema_TCA': [], 'Micromonospora_TCA': [], 'Actinoplanes_TCA': [], 'Microlunatus_TCA': [], 'Slackia_TCA': [], 'Cryptobacterium_TCA': [], 'Anaerolinea_TCA': [], 'Chloroflexus_TCA': [], 'Caldilinea_TCA': [], 'Roseiflexus_TCA': [], 'Dehalogenimonas_TCA': [], 'Thermomicrobium_TCA': [], 'Eggerthella_TCA': [], 'Atopobium_TCA': [], 'Dehalococcoides_TCA': [], 'Saccharomonospora_TCA': [], 'Brevibacterium_TCA': [], 'Dermabacter_TCA': [], 'Leifsonia_TCA': [], 'Saccharopolyspora_TCA': [], 'Leucobacter_TCA': [], 'Rothia_TCA': [], 'Agromyces_TCA': [], 'Actinoalloteichus_TCA': [], 'Sciscionella_TCA': [], 'Actinocatenispora_TCA': [], 'Lechevalieria_TCA': [], 'Thermobifida_TCA': [], 'Thermocrispum_TCA': [], 'Actinokineospora_TCA': [], 'Lentzea_TCA': [], 'Kibdelosporangium_TCA': [], 'Actinomycetospora_TCA': [], 'Saccharothrix_TCA': [], 'Collinsella_TCA': [], 'Propionimicrobium_TCA': [], 'Thermogemmatispora_TCA': [], 'Granulicoccus_TCA': [], 'Nitrolancea_TCA': [], 'Enterorhabdus_TCA': [], 'Olsenella_TCA': [], 'Senegalimassilia_TCA': [], 'Oscillochloris_TCA': [], 'Enorma_TCA': [], 'Aestuariimicrobium_TCA': [], 'Oerskovia_TCA': [], 'Cellulomonas_TCA': [], 'Paraoerskovia_TCA': [], 'Isoptericola_TCA': [], 'Brachybacterium_TCA': [], 'Promicromonospora_TCA': [], 'Cellulosimicrobium_TCA': [], 'Xylanimonas_TCA': [], 'Kineococcus_TCA': [], 'Kineosporia_TCA': [], 'Dietzia_TCA': [], 'Sphaerobacter_TCA': [], 'Ktedonobacter_TCA': [], 'Tomitella_TCA': [], 'Segniliparus_TCA': [], 'Catenulispora_TCA': [], 'Tsukamurella_TCA': [], 'Actinopolyspora_TCA': [], 'Williamsia_TCA': [], 'Actinospica_TCA': [], 'Tropheryma_TCA': [], 'Beutenbergia_TCA': [], 'Kytococcus_TCA': [], 'Nakamurella_TCA': [], 'Blastococcus_TCA': [], 'Candidatus_TCA': [], 'Acidothermus_TCA': [], 'Lysinimicrobium_TCA': [], 'Geodermatophilaceae_TCA': [], 'Glycomyces_TCA': [], 'Cryptosporangium_TCA': [], 'Ruania_TCA': [], 'Modestobacter_TCA': [], 'Dermatophilus_TCA': [], 'Dermacoccus_TCA': [], 'Mobilicoccus_TCA': [], 'Demetria_TCA': [], 'Sanguibacter_TCA': [], 'Kineosphaera_TCA': [], 'Timonella_TCA': [], 'Geodermatophilus_TCA': [], 'Haloglycomyces_TCA': [], 'Sporichthya_TCA': [], 'Jonesia_TCA': [], 'Janibacter_TCA': [], 'Gryllotalpicola_TCA': [], 'Frigoribacterium_TCA': [], 'Herbiconiux_TCA': [], 'Knoellia_TCA': [], 'Humibacter_TCA': [], 'Serinicoccus_TCA': [], 'Rathayibacter_TCA': [], 'Agreia_TCA': [], 'Cryobacterium_TCA': [], 'Intrasporangium_TCA': [], 'Tetrasphaera_TCA': [], 'Intrasporangiaceae_TCA': [], 'Glaciibacter_TCA': [], 'Cryocola_TCA': [], 'Austwickia_TCA': [], 'Arsenicicoccus_TCA': [], 'Salinibacterium_TCA': [], 'Gulosibacter_TCA': [], 'Pseudoclavibacter_TCA': [], 'Mycetocola_TCA': [], 'Phycicoccus_TCA': [], 'Ornithinimicrobium_TCA': [], 'Terrabacter_TCA': [], 'Terracoccus_TCA': [], 'Alloscardovia_TCA': [], 'Actinopolymorpha_TCA': [], 'Kribbella_TCA': [], 'Agrococcus_TCA': [], 'Micrococcus_TCA': [], 'Nesterenkonia_TCA': [], 'Dactylosporangium_TCA': [], 'Nocardioidaceae_TCA': [], 'Aeromicrobium_TCA': [], 'Renibacterium_TCA': [], 'Verrucosispora_TCA': [], 'Sinomonas_TCA': [], 'Marmoricola_TCA': [], 'Pilimelia_TCA': [], 'Catelliglobosispora_TCA': [], 'Acaricomes_TCA': [], 'Catenuloplanes_TCA': [], 'Jiangella_TCA': [], 'Longispora_TCA': [], 'Yaniella_TCA': [], 'Actinomadura_TCA': [], 'Microbispora_TCA': [], 'Nonomuraea_TCA': [], 'Streptosporangium_TCA': [], 'Streptomonospora_TCA': [], 'Microtetraspora_TCA': [], 'Prauserella_TCA': [], 'Allokutzneria_TCA': [], 'Streptacidiphilus_TCA': [], 'Mobiluncus_TCA': [], 'Trueperella_TCA': [], 'Varibaculum_TCA': [], 'Actinobaculum_TCA': [], 'Arcanobacterium_TCA': [], 'Pimelobacter_TCA': [], 'Corynebacteriales_TCA': [], 'Parascardovia_TCA': [], 'Thermomonospora_TCA': [], 'Luteipulveratus_TCA': [], 'Actinotignum_TCA': [], 'Alloactinosynnema_TCA': [], 'Amycolicicoccus_TCA': [], 'Actinobacteria_TCA': [], 'Scardovia_TCA': [], 'Adlercreutzia_TCA': [], 'Coriobacteriaceae_TCA': [], 'Gardnerella_TCA': [], 'Rubrobacter_TCA': [], 'Gordonibacter_TCA': [], 'Stackebrandtia_TCA': [], 'Thermobispora_TCA': [], 'Coriobacterium_TCA': [], 'Conexibacter_TCA': [], 'Ilumatobacter_TCA': [], 'Cellvibrio_TCA': [], 'Acidimicrobium_TCA': [], 'Actinokieospora_TCA': []}
    return holes