# Assembly Pathway Tutorial

*Keith Y. Patarroyo*  
*University of Glasgow*

*keith.patarroyo@glasgow.ac.uk*

*January 04, 2023*
___

## Table of Contents
0. [Requirements](#requirements)
1. [Data Structures](#datastructures) 
    - 1.1 [Molecule Representations](#molrep)
    - 1.2 [Multiple Molecules Representations](#molreps)
    - 1.3 [Pathway Representations](#pathreps)
2. [Reconstruction Algorithm](#rec_algo)
    - 2.1 [Duplicate Structures Construction](#dup_struc)
    - 2.2 [Duplicate Structures Union + Remnant](#final)
3. [Output Options](#free)
    - 3.1 [AssemblyGo Log Output](#arap)
    - 3.2 [Inchi + Digraph](#lscm)
    - 3.3 [Cytoscape Digraph](#artistic)
    - 3.4 [Mathematica Output](#extra)
4. [Discussion](#discussion)
5. [Extra Topics](#discussion)
    - 5.1 [String Pathway Reconstruction](#sample)
6. [References](#references)
___

# Assembly Pathway Tutorial

AssemblyPathway is a simple, and fast 2d assembly space viewer based on RDKit and and cytopscape.js

##  <a id='requirements'></a>0. Requirements

Only two libraries are needed for this notebook:

    - AssemblyCalculator 
    - RDKit
RDKit can be downloaded from pip:
pip install rdkit
AssemblyCalculator can be downloaded from [GitLab](https://gitlab.com/croningroup/cheminformatics/assemblycalculator).

Two tutorials and documentation for RDKit can be found [here](https://www.rdkit.org/docs/GettingStartedInPython.html) and [here](https://rdkit.org/docs/Cookbook.html).


```python
from assembly_pathv2 import generate_pathway, check_edge_in_lista, equivalence,check_edge_in_list
import os
from log_parsing import encode_path_data
import assemblycalculator as ac
from rdkit import Chem
from rdkit.Chem import Draw
import igraph
from igraph import GraphBase
from rdkit.Chem.Draw import rdMolDraw2D
```

#  <a id='datastructures'></a> 1. Data Structures

In this tutorial we'll be concerned with the data structures that represent a molecule and a pathway in a computer. 

##  <a id='molrep'></a> 1.1 Molecule Representations

In principle there are several formats for representing a molecule, SMILES, InChI, Mol. For the porpouses of molecular assembly, the three formats are equivalent since we will only be using the atom and the bond information of a molecule. This information can be stored in two tables:

<table>
<tr><th>Atoms Table </th><th>Bonds Table</th></tr>
<tr><td>

| Atom # | Atom Type |
|:------:|:---------:|
|    0   |     C     |
|    1   |     C     |
|    2   |     C     |
|    3   |     C     |
|    4   |     C     |
|    5   |     C     |
|    6   |     C     |
|    7   |     O     |
|    8   |     O     |

</td><td>

| Atom 1 | Atom 2 | Bond Order |
|:------:|:------:|:----------:|
|    0   |    1   |      2     |
|    0   |    2   |      1     |
|    1   |    3   |      1     |
|    2   |    4   |      2     |
|    3   |    5   |      2     |
|    4   |    5   |      1     |
|    5   |    6   |      1     |
|    6   |    7   |      2     |
|    6   |    8   |      1     |
</td></tr> </table>

This is the connection tables information for benzodic acid.

This two tables can also be thought as one labeled graph $(G,\ell)=\big((V,E),\ell\big)$,

$$
\begin{align}
V & = [0,1,2,3,4,5,6,7,8]\\
V_\ell& =[\text{C},\text{C},\text{C},\text{C},\text{C},\text{C},\text{C},\text{O},\text{O}]\\
E& =\Big[(0,1),(0,2),(1,3),(2,4),(3,5),(4,5),(5,6),(6,7),(6,8)\Big]\\
E_\ell& =[2,1,1,2,2,1,1,2,1]
\end{align}
$$

We can get the InChI, Mol, and SMILES representation of benzodic acid:


```python
m_1 = Chem.MolFromInchi("InChI=1S/C7H6O2/c8-7(9)6-4-2-1-3-5-6/h1-5H,(H,8,9)")
m_2 = Chem.MolFromSmiles("O=C(O)c1ccccc1")
mol_file = "\n     RDKit\n\n  9  9  0  0  0  0  0  0  0  0999 V2000\n    3.7500    1.2990    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n    3.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    3.7500   -1.2990    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n    1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.7500   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.7500   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.7500    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.7500    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  2  0\n  2  3  1  0\n  2  4  1  0\n  4  5  2  0\n  5  6  1  0\n  6  7  2  0\n  7  8  1  0\n  8  9  2  0\n  9  4  1  0\nM  END\n"
m_3 = Chem.MolFromMolBlock(mol_file)
Draw.MolsToGridImage([m_1,m_2,m_3])
```




    
![png](output_16_0.png)
    




```python
def mol2tables(mol):
    atoms_info = [ (atom.GetIdx(), atom.GetSymbol()) for atom in mol.GetAtoms()]
    bonds_info = [(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), bond.GetBondTypeAsDouble()) for bond in mol.GetBonds()]
    tables = (atoms_info,bonds_info)
    return tables
    
def mol2graph(mol):
    atoms_info,bonds_info = mol2tables(mol)
    graph = igraph.Graph()
    for atom_info in atoms_info:
        graph.add_vertex(atom_info[0], AtomicSymbole=atom_info[1])
    for bond_info in bonds_info:
        graph.add_edge(bond_info[0], bond_info[1], BondTypeAsDouble=bond_info[2])
    return graph
def plot_graph(graph):
    layout = graph.layout_graphopt()
    color_dict_vertex = {"C": "gray", "O" :  "red","N" :  "blue" }
    color_dict_edge = {1.0: "black", 2.0 :  "green" }
    my_plot = igraph.Plot()
    my_plot.add(graph,layout=layout,bbox=(300,300),
                margin=20,
                vertex_color = [color_dict_vertex[atom] for atom in graph.vs["AtomicSymbole"]],
                vertex_size = [ 10 for v in graph.vs ],
                edge_color = [color_dict_edge[ed] for ed in graph.es["BondTypeAsDouble"]])
    return my_plot
```


```python
Chem.Kekulize(m_1) # Convert Aromatic bonds to covalent bonds
mol2tables(m_1)
```




    ([(0, 'C'),
      (1, 'C'),
      (2, 'C'),
      (3, 'C'),
      (4, 'C'),
      (5, 'C'),
      (6, 'C'),
      (7, 'O'),
      (8, 'O')],
     [(0, 1, 2.0),
      (0, 2, 1.0),
      (1, 3, 1.0),
      (2, 4, 2.0),
      (3, 5, 2.0),
      (4, 5, 1.0),
      (5, 6, 1.0),
      (6, 7, 2.0),
      (6, 8, 1.0)])




```python
graph = mol2graph(m_1)
layout = graph.layout_graphopt()
color_dict_vertex = {"C": "gray", "O" :  "red" }
color_dict_edge = {1.0: "black", 2.0 :  "green" }
my_plot = igraph.Plot()
my_plot.add(graph,layout=layout,bbox=(300,300),
             margin=20,
             vertex_color = [color_dict_vertex[atom] for atom in graph.vs["AtomicSymbole"]],
             vertex_size = [ 10 for v in graph.vs ],
             edge_color = [color_dict_edge[ed] for ed in graph.es["BondTypeAsDouble"]])
my_plot
```




    
![svg](output_19_0.svg)
    



##  <a id='molreps'></a> 1.2 Multiple Molecules Representations


```python
combined = Chem.MolFromInchi("InChI=1S/C5H6N2O2.C5H5N5/c1-3-2-6-5(9)7-4(3)8;6-4-3-5(9-1-7-3)10-2-8-4/h2H,1H3,(H2,6,7,8,9);1-2H,(H3,6,7,8,9,10)")
```


```python
combined
```




    
![png](output_22_0.png)
    




```python
Chem.Kekulize(combined) # Convert Aromatic bonds to covalent bonds
mol2tables(combined)
```




    ([(0, 'C'),
      (1, 'C'),
      (2, 'C'),
      (3, 'C'),
      (4, 'C'),
      (5, 'N'),
      (6, 'N'),
      (7, 'O'),
      (8, 'O'),
      (9, 'C'),
      (10, 'C'),
      (11, 'C'),
      (12, 'C'),
      (13, 'C'),
      (14, 'N'),
      (15, 'N'),
      (16, 'N'),
      (17, 'N'),
      (18, 'N')],
     [(0, 2, 1.0),
      (1, 2, 2.0),
      (1, 5, 1.0),
      (2, 3, 1.0),
      (3, 6, 2.0),
      (3, 7, 1.0),
      (4, 5, 2.0),
      (4, 6, 1.0),
      (4, 8, 1.0),
      (9, 15, 2.0),
      (9, 17, 1.0),
      (10, 16, 2.0),
      (10, 18, 1.0),
      (11, 12, 2.0),
      (11, 13, 1.0),
      (11, 15, 1.0),
      (12, 14, 1.0),
      (12, 16, 1.0),
      (13, 17, 2.0),
      (13, 18, 1.0)])




```python
graph = mol2graph(combined)
layout = graph.layout_graphopt()
color_dict_vertex = {"C": "gray", "O" :  "red","N" :  "blue" }
color_dict_edge = {1.0: "black", 2.0 :  "green" }
my_plot = igraph.Plot()
my_plot.add(graph,layout=layout,bbox=(300,300),
             margin=20,
             vertex_color = [color_dict_vertex[atom] for atom in graph.vs["AtomicSymbole"]],
             vertex_size = [ 10 for v in graph.vs ],
             edge_color = [color_dict_edge[ed] for ed in graph.es["BondTypeAsDouble"]])
my_plot
```




    
![svg](output_24_0.svg)
    



##  <a id='pathreps'></a> 1.3 Pathway Representations

### One molecule

Consider a pathway output file. It has the following format
_DATE_ Running on file:  aspirin.mol
_DATE_ ORIGINAL GRAPH
+++++++++++++++
Vertices [0 1 2 3 4 5 6 7 8 9 10 11 12]
Edges [[0 1] [2 0] [0 3] [1 4] [5 2] [2 6] [3 10] [4 7] [7 5] [6 8] [6 9] [10 11] [10 12]]
VertexColours [C C C O C C C C O O C O C]
EdgeColours [double single single single double single single double single double single double single]
+++++++++++++++
PATHWAY
Pathway Graphs
======
Vertices [0 1 2]
Edges [[0 1] [2 0]]
VertexColours [C C C]
EdgeColours [double single]
======
======
Vertices [2 6 8 9]
Edges [[2 6] [6 8] [6 9]]
VertexColours [C C O O]
EdgeColours [single double single]
======
======
Vertices [1 4 7]
Edges [[1 4] [4 7]]
VertexColours [C C C]
EdgeColours [single double]
======
----------
Remnant Graph
Vertices [0 3 15 10 11 12 5 13 14]
Edges [[0 3] [15 10] [10 11] [10 12] [5 13] [14 5]]
VertexColours [C O O C O C C C C]
EdgeColours [single single double single double single]
----------
Duplicated Edges
{[[0 1] [2 0]] [[5 2] [7 5]]}
{[[2 6] [6 8] [6 9]] [[3 10] [10 11] [10 12]]}
{[[1 4] [4 7]] [[5 13] [14 5]]}
+++++++++++++++
###############
Atom Equivalents
[2 13]
[7 14]
[3 15]
###############
_DATE_ Assembly Index:  8
_DATE_ Time:  0.2769666
In summary the pathway information can be compressed to the following information:
Original Graph
Vertices [0 1 2 3 4 5 6 7 8 9 10 11 12]
Edges [[0 1] [2 0] [0 3] [1 4] [5 2] [2 6] [3 10] [4 7] [7 5] [6 8] [6 9] [10 11] [10 12]]
VertexColours [C C C O C C C C O O C O C]
EdgeColours [double single single single double single single double single double single double single]
----------
Remnant Graph
Vertices [0 3 15 10 11 12 5 13 14]
Edges [[0 3] [15 10] [10 11] [10 12] [5 13] [14 5]]
VertexColours [C O O C O C C C C]
EdgeColours [single single double single double single]
----------
Duplicated Edges
{[[0 1] [2 0]] [[5 2] [7 5]]}
{[[2 6] [6 8] [6 9]] [[3 10] [10 11] [10 12]]}
{[[1 4] [4 7]] [[5 13] [14 5]]}
----------
Atom Equivalents
[2 13]
[7 14]
[3 15]
#### Original Graph


```python
aspirin_mol = 'BENZOIC ACID, 2-(ACETYLOXY)-, ID: C50782\n     RDKit          2D\n\n 13 13  0  0  0  0  0  0  0  0999 V2000\n    1.7434    1.4944    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.7434    0.4981    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.8966    1.9925    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    2.5902    1.9925    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n    0.8966    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000    1.4944    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.8966    2.9887    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000    0.4981    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000    3.4869    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n    1.7434    3.4869    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n    3.4869    1.4944    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    4.3337    1.9925    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n    3.4869    0.4981    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  2  0\n  3  1  1  0\n  1  4  1  0\n  2  5  1  0\n  6  3  2  0\n  3  7  1  0\n  4 11  1  0\n  5  8  2  0\n  8  6  1  0\n  7  9  2  0\n  7 10  1  0\n 11 12  2  0\n 11 13  1  0\nM  END\n'
original_m = Chem.MolFromMolBlock(aspirin_mol)
```


```python
original_m
```




    
![png](output_33_0.png)
    




```python
Chem.Kekulize(original_m)
graph_original = mol2graph(original_m)
layout = graph_original.layout_graphopt()
color_dict_vertex = {"C": "gray", "O" :  "red","N" :  "blue" }
color_dict_edge = {1.0: "black", 2.0 :  "green" }
my_plot = igraph.Plot()
my_plot.add(graph_original,layout=layout,bbox=(300,300),
             margin=20,
             vertex_color = [color_dict_vertex[atom] for atom in graph_original.vs["AtomicSymbole"]],
             vertex_size = [ 10 for v in graph_original.vs ],
             edge_color = [color_dict_edge[ed] for ed in graph_original.es["BondTypeAsDouble"]])
my_plot
```




    
![svg](output_34_0.svg)
    



#### Remnant


```python
remnant_v =  [0, 3, 15, 10, 11, 12, 5, 13, 14]
remnant_e =  [[0, 3] ,[15, 10], [10, 11], [10, 12], [5, 13], [14, 5]]
remant_vl = ["C", "O", "O", "C", "O", "C", "C", "C", "C"]
remant_el = [1.0, 1.0, 2.0, 1.0, 2.0, 1.0]
equivalences = [[2, 13],[7, 14],[3, 15]]
duplicates = [[[[0, 1], [2, 0]], [[5, 2], [7, 5]]],[[[2, 6], [6, 8], [6, 9]], [[3, 10], [10, 11], [10, 12]]],[[[1, 4], [4, 7]], [[5, 13], [14, 5]]]]
remnant_e_mod = equivalence([remnant_e],equivalences)[0]
remnant_v_mod = list({item for sublist in remnant_e_mod for item in sublist})
duplicates_mod = [equivalence(rep,equivalences) for rep in duplicates]
duplicates_mod_plot = [[item for sublist in dup for item in sublist] for dup in duplicates_mod]
```


```python
graph = igraph.Graph()
for i,atom_info in enumerate(remnant_v_mod):
    graph.add_vertex(id=atom_info, AtomicSymbole=remant_vl[i])
for i,bond_info in enumerate(remnant_e_mod):
    graph.add_edge(graph.vs["id"].index(bond_info[0]), graph.vs["id"].index(bond_info[1]), BondTypeAsDouble=remant_el[i])
```


```python
plot_graph(graph)
```




    
![svg](output_38_0.svg)
    




```python
highlight_bonds = []
highlight_atoms_pre = []
for i,bond in enumerate(graph_original.get_edgelist()):
    if  check_edge_in_list(list(bond),remnant_e_mod):
        highlight_bonds.append(i)
        highlight_atoms_pre.append(bond)
highlight_atoms = {item for sublist in highlight_atoms_pre for item in sublist}
Draw.MolsToGridImage([original_m],highlightAtomLists=[highlight_atoms],highlightBondLists=[highlight_bonds])
```




    
![png](output_39_0.png)
    



#### Duplicates


```python
highlight_bonds_total = []
highlight_atoms_total = []
legends=[]
for j, dup in enumerate(duplicates_mod_plot):
    highlight_bonds = []
    highlight_atoms_pre = []
    for i,bond in enumerate(graph_original.get_edgelist()):
        if  check_edge_in_list(list(bond),dup):
            highlight_bonds.append(i)
            highlight_atoms_pre.append(bond)
    highlight_atoms = {item for sublist in highlight_atoms_pre for item in sublist}
    highlight_bonds_total.append(highlight_bonds)
    highlight_atoms_total.append(highlight_atoms)
    legends.append("Duplicate {}".format(j+1))
Draw.MolsToGridImage([original_m for i in range(len(duplicates_mod_plot))],highlightAtomLists=highlight_atoms_total,highlightBondLists=highlight_bonds_total,legends=legends)
```




    
![png](output_41_0.png)
    



### Multiple Molecules

Consider the compressed version of a pathway output file of two molecules:
Original Graph
Vertices [0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18]
Edges [[0 2] [1 2] [1 5] [2 3] [3 6] [3 7] [4 5] [4 6] [4 8] [9 15] [9 17] [10 16] [10 18] [11 12] [11 13] [11 15] [12 14] [12 16] [13 17] [13 18]]
VertexColours [C C C C C N N O O C C C C C N N N N N]
EdgeColours [single double single single double single double single single double single double single double single single single single double single] 
----------
Remnant Graph
Vertices [2 3 12 14 19 20 21 22 7 10 18 24 16]
Edges [[2 3] [12 14] [19 20] [19 21] [22 7] [10 18] [24 16] [20 16]]
VertexColours [C C C N C C C C O C N C N C]
EdgeColours [single single double single single single double single] 
----------
Duplicated Edges
{[[0 2] [1 2]] [[11 12] [11 13]]}
{[[1 5] [4 5] [4 8]] [[3 6] [3 7] [4 6]]}
{[[9 17] [13 17] [13 18]] [[10 16] [10 18] [12 16]]}
{[[22 6] [4 6]] [[10 16] [20 16]]}
{[[9 15] [11 15]] [[24 16] [20 16]]}
----------
Atom Equivalents
[11 19]
[12 20]
[13 21]
[3 22]
[10 24]
#### Original


```python
combined_mol = '\n     RDKit          2D\n\n 19 20  0  0  0  0  0  0  0  0999 V2000\n    1.5000   -2.5981    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.7500   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.7500   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.7500    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.7500    1.2990    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.5000    0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.5000   -2.5981    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.5000    2.5981    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n    7.4063    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    2.5000   -0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    5.0981    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    3.7990    1.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    5.0981   -0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    3.7990    3.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n    6.5247    1.2135    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n    2.5000    0.7500    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n    6.5247   -1.2135    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n    3.7990   -1.5000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n  1  3  1  0\n  2  3  2  0\n  2  6  1  0\n  3  4  1  0\n  4  7  2  0\n  4  8  1  0\n  5  6  2  0\n  5  7  1  0\n  5  9  1  0\n 10 16  2  0\n 10 18  1  0\n 11 17  2  0\n 11 19  1  0\n 12 13  2  0\n 12 14  1  0\n 12 16  1  0\n 13 15  1  0\n 13 17  1  0\n 14 18  2  0\n 14 19  1  0\nM  END\n'
original_comb = Chem.MolFromMolBlock(combined_mol)
original_comb
```




    
![png](output_46_0.png)
    




```python
Chem.Kekulize(original_comb)
graph_original_comb = mol2graph(original_comb)
layout = graph_original.layout_graphopt()
color_dict_vertex = {"C": "gray", "O" :  "red","N" :  "blue" }
color_dict_edge = {1.0: "black", 2.0 :  "green" }
my_plot = igraph.Plot()
my_plot.add(graph_original,layout=layout,bbox=(300,300),
             margin=20,
             vertex_color = [color_dict_vertex[atom] for atom in graph_original.vs["AtomicSymbole"]],
             vertex_size = [ 10 for v in graph_original.vs ],
             edge_color = [color_dict_edge[ed] for ed in graph_original.es["BondTypeAsDouble"]])
my_plot
```




    
![svg](output_47_0.svg)
    



#### Remnant


```python
remnant_v =  [2, 3, 12, 14, 19, 20, 21, 22, 7, 10, 18, 24, 16]
remnant_e =  [[2, 3], [12, 14], [19, 20], [19, 21], [22, 7], [10, 18], [24, 16], [20, 16]]
remant_vl = ["C", "C", "C", "N", "C", "C", "C", "C", "O", "C", "N", "C", "N", "C"]
remant_el = [1.0, 1.0, 2.0, 1.0, 1.0, 1.0, 2.0, 1.0]
equivalences = [[11, 19],[12, 20],[13, 21],[3, 22],[10, 24]]
duplicates = [[[[0, 2], [1, 2]], [[11, 12], [11, 13]]],[[[1, 5], [4, 5], [4, 8]], [[3, 6], [3, 7], [4, 6]]],[[[9, 17], [13, 17], [13, 18]], [[10, 16], [10, 18], [12, 16]]],[[[22, 6], [4, 6]], [[10, 16], [20, 16]]],[[[9, 15], [11, 15]], [[24, 16], [20, 16]]]]
remnant_e_mod = equivalence([remnant_e],equivalences)[0]
remnant_v_mod = list({item for sublist in remnant_e_mod for item in sublist})
duplicates_mod = [equivalence(rep,equivalences) for rep in duplicates]
duplicates_mod_plot = [[item for sublist in dup for item in sublist] for dup in duplicates_mod]
```


```python
graph = igraph.Graph()
for i,atom_info in enumerate(remnant_v_mod):
    graph.add_vertex(id=atom_info, AtomicSymbole=remant_vl[i])
for i,bond_info in enumerate(remnant_e_mod):
    graph.add_edge(graph.vs["id"].index(bond_info[0]), graph.vs["id"].index(bond_info[1]), BondTypeAsDouble=remant_el[i])
plot_graph(graph)
```




    
![svg](output_50_0.svg)
    




```python
highlight_bonds = []
highlight_atoms_pre = []
for i,bond in enumerate(graph_original_comb.get_edgelist()):
    if  check_edge_in_list(list(bond),remnant_e_mod):
        highlight_bonds.append(i)
        highlight_atoms_pre.append(bond)
highlight_atoms = {item for sublist in highlight_atoms_pre for item in sublist}
Draw.MolsToGridImage([original_comb],highlightAtomLists=[highlight_atoms],highlightBondLists=[highlight_bonds])
```




    
![png](output_51_0.png)
    



#### Duplicates


```python
highlight_bonds_total = []
highlight_atoms_total = []
legends=[]
for j,dup in enumerate(duplicates_mod_plot):
    highlight_bonds = []
    highlight_atoms_pre = []
    for i,bond in enumerate(graph_original_comb.get_edgelist()):
        if  check_edge_in_list(list(bond),dup):
            highlight_bonds.append(i)
            highlight_atoms_pre.append(bond)
    highlight_atoms = {item for sublist in highlight_atoms_pre for item in sublist}
    highlight_bonds_total.append(highlight_bonds)
    highlight_atoms_total.append(highlight_atoms)
    legends.append("Duplicate {}".format(j+1))
Draw.MolsToGridImage([original_comb for i in range(len(duplicates_mod_plot))],highlightAtomLists=highlight_atoms_total,highlightBondLists=highlight_bonds_total,legends=legends)
```




    
![png](output_53_0.png)
    



#  <a id='rec_algo'></a>2. Reconstruction Algorithm

The reconstruction algorithm is based on the observation that if we construct the remanant and all the duplicates, we have all the pieces to construct the molecule.

#### One Molecule


```python
cumulative = [duplicates_mod_plot[0]]
legends=["remnant","duplicate 1"]
for i,dup in enumerate(duplicates_mod_plot[1:]):
    cumulative.append(cumulative[-1]+dup)
    legends.append(legends[-1]+"+"+str(i+2))
#duplicates_mod_plot_2 = [item for sublist in duplicates_mod_plot for item in sublist]

total = cumulative[-1] + remnant_e_mod
highlight_two = [remnant_e_mod] + cumulative + [total]
legends.append("remnant + "+ legends[-1])
```


```python
highlight_bonds_total = []
highlight_atoms_total = []
for dup in highlight_two:
    highlight_bonds = []
    highlight_atoms_pre = []
    for i,bond in enumerate(graph_original.get_edgelist()):
        if  check_edge_in_list(list(bond),dup):
            highlight_bonds.append(i)
            highlight_atoms_pre.append(bond)
    highlight_atoms = {item for sublist in highlight_atoms_pre for item in sublist}
    highlight_bonds_total.append(highlight_bonds)
    highlight_atoms_total.append(highlight_atoms)
Draw.MolsToGridImage([original_m for i in range(len(highlight_two))],molsPerRow=4,highlightAtomLists=highlight_atoms_total,highlightBondLists=highlight_bonds_total,legends=legends)
```




    
![png](output_58_0.png)
    



#### Multiple Molecules


```python
cumulative = [duplicates_mod_plot[0]]
legends=["remnant","duplicate 1"]
for i,dup in enumerate(duplicates_mod_plot[1:]):
    cumulative.append(cumulative[-1]+dup)
    legends.append(legends[-1]+"+"+str(i+2))
#duplicates_mod_plot_2 = [item for sublist in duplicates_mod_plot for item in sublist]

total = cumulative[-1] + remnant_e_mod
highlight_two = [remnant_e_mod] + cumulative + [total]
legends.append("remnant + "+ legends[-1])
```


```python
highlight_bonds_total = []
highlight_atoms_total = []
for dup in highlight_two:
    highlight_bonds = []
    highlight_atoms_pre = []
    for i,bond in enumerate(graph_original_comb.get_edgelist()):
        if  check_edge_in_list(list(bond),dup):
            highlight_bonds.append(i)
            highlight_atoms_pre.append(bond)
    highlight_atoms = {item for sublist in highlight_atoms_pre for item in sublist}
    highlight_bonds_total.append(highlight_bonds)
    highlight_atoms_total.append(highlight_atoms)
Draw.MolsToGridImage([original_comb for i in range(len(highlight_two))],molsPerRow=4,highlightAtomLists=highlight_atoms_total,highlightBondLists=highlight_bonds_total,legends=legends)
```




    
![png](output_61_0.png)
    



The reconstruction algorithm is composed of two fundamental steps, first we construct each duplicate structure. Then we join the duplicate structures and the remnant graph.

##  <a id='dup_struc'></a>2.1 Duplicate Structures Construction



##  <a id='final'></a>2.2 Duplicate Structures Union + Remnant


```python
aspirin_smiles = "C[C@H](CC)[C@@H](N)C(N/C(C(N[C@@H]1C(N[C@@H](C(NC(C(N[C@@H](C(N[C@@H](C(N[C@H](C(N2[C@@H](C(NC3)=O)CCC2)=O)[C@@H](C)SC[C@H](C(N[C@H](CCCCN)C(N[C@H](C(NCC(N[C@@H](C(N[C@H](CC(C)C)C(N[C@H](CCSC)C(NCC(N[C@@H](C(N[C@H](CC(N)=O)C(N[C@H](CCSC)C(N[C@H](CCCCN)C(N[C@@H]([C@@H](C)SC[C@@H](N6)C(N[C@H](CC5=CN=CN5)C(N[C@@H](C(N[C@H](CO)C(N[C@H]([C@@H](CC)C)C(N[C@H](CC8=CN=CN8)C(N[C@H]([C@@H](C)C)C(NC(C(N[C@H](CCCCN)C(O)=O)=O)=C)=O)=O)=O)=O)=O)CS[C@H](C)[C@H]7C6=O)=O)=O)C(N[C@H](C)C(N7)=O)=O)=O)=O)=O)=O)CS4)=O)=O)=O)=O)C)=O)=O)[C@H]4C)=O)=O)NC3=O)=O)CSC1)=O)CC(C)C)=O)=C)=O)[C@H](C)CC)=O)=O)=C\C)=O"
ac.molecular_exact(aspirin_smiles, timeout=2.0)

```


    ---------------------------------------------------------------------------

    IndexError                                Traceback (most recent call last)

    /mnt/scapa4/group/Keith Patarroyo/Projects/AssemblyCalculator/AssemblyPath/Path-Tutorial.ipynb Cell 65 in <cell line: 2>()
          <a href='vscode-notebook-cell://wsl%2Bubuntu-22.04/mnt/scapa4/group/Keith%20Patarroyo/Projects/AssemblyCalculator/AssemblyPath/Path-Tutorial.ipynb#Y206sdnNjb2RlLXJlbW90ZQ%3D%3D?line=0'>1</a> aspirin_smiles = "C[C@H](CC)[C@@H](N)C(N/C(C(N[C@@H]1C(N[C@@H](C(NC(C(N[C@@H](C(N[C@@H](C(N[C@H](C(N2[C@@H](C(NC3)=O)CCC2)=O)[C@@H](C)SC[C@H](C(N[C@H](CCCCN)C(N[C@H](C(NCC(N[C@@H](C(N[C@H](CC(C)C)C(N[C@H](CCSC)C(NCC(N[C@@H](C(N[C@H](CC(N)=O)C(N[C@H](CCSC)C(N[C@H](CCCCN)C(N[C@@H]([C@@H](C)SC[C@@H](N6)C(N[C@H](CC5=CN=CN5)C(N[C@@H](C(N[C@H](CO)C(N[C@H]([C@@H](CC)C)C(N[C@H](CC8=CN=CN8)C(N[C@H]([C@@H](C)C)C(NC(C(N[C@H](CCCCN)C(O)=O)=O)=C)=O)=O)=O)=O)=O)CS[C@H](C)[C@H]7C6=O)=O)=O)C(N[C@H](C)C(N7)=O)=O)=O)=O)=O)=O)CS4)=O)=O)=O)=O)C)=O)=O)[C@H]4C)=O)=O)NC3=O)=O)CSC1)=O)CC(C)C)=O)=C)=O)[C@H](C)CC)=O)=O)=C\C)=O"
    ----> <a href='vscode-notebook-cell://wsl%2Bubuntu-22.04/mnt/scapa4/group/Keith%20Patarroyo/Projects/AssemblyCalculator/AssemblyPath/Path-Tutorial.ipynb#Y206sdnNjb2RlLXJlbW90ZQ%3D%3D?line=1'>2</a> ac.molecular_exact(aspirin_smiles, timeout=2.0)


    File ~/.local/lib/python3.10/site-packages/assemblycalculator/molecules.py:31, in molecular_exact(molecule, timeout)
         29 # Parse the output
         30 if success:
    ---> 31     results = parse_go_string(go_output)
         32     results["molecule"] = molecule
         33     results["num_bonds"] = num_bonds


    File ~/.local/lib/python3.10/site-packages/assemblycalculator/molecules.py:215, in parse_go_string(output_string, record_pathway)
        213 """Parse the output from the Go executable for Assembly"""
        214 # print(output_string)
    --> 215 index_line = re.findall(string=output_string,pattern=r"\nAssembly Index:.*")[0]
        216 index_string = index_line.split(" ")[-1].strip()
        217 MA = int(index_string)


    IndexError: list index out of range



```python

```
