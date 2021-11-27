# -*- coding: utf-8 -*-
"""
Created on Fri Nov 26 15:18:21 2021

@author: vfranco-
"""
class BooleanNetwork:
    def __init__(self, *args, file = None):
        """Fill in the class fields the information about gene names and functoins.

        *args -> list the list of gene names followed by list of gene functions.
        file -> String with filename conteining gene names and functions.
            file is an alternative method to *args for filling the information.
        """
        if len(args) > 1:
            self.genes = args[0]
            self.func = args[1]
        elif (file != None):
            self.read_data(file)

    def read_data (self, FileBase):
        """Read file to get gene names and regulatory functions.

        FileBase -> name of the file.
        The content of the file must, in each line, the name of
        the gene followed by its function
        Exemple:
            v1, v4 | v2
            v2, !v4&v5 | v4&!v5 | v3&v5 | v3&v4
            v3, v1&!v3&!v5 | v3&v5 | !v1&v5 | !v1&v3
            v4, v4
            v5, v1&v4&v5
        """
        f = open(FileBase, "r")
        s = f.readlines()
        genes = list()
        func = list()
        for line in sorted(s):
            gene, fun = line.strip().split(', ')
            genes.append(gene.strip())
            func.append(fun.strip())
        f.close()
        self.genes = genes
        self.func = func

    def gene_combination(self, gene_states = ["0", "1"]):
        """Generate all genome configurations from individual gene states.

        gene_states -> optional list of individual gene states
        Create new instance variable with the generated information
        """
        from itertools import product
        self.genome_states = list(product(gene_states, repeat=len(self.genes)))

    def print_gene_combination(self):
        """Print all possible genome configurations.
        """
        for e in self.genome_states:
            s = ''.join(e)
            print(s)

    def gene_dict(self):
        """Make a table of genome states associated with each gene state

        Create new instance variable with the generated information
        """
        d = {}
        for i, g in enumerate(self.genes):
            d[g] = []
            notg = '!' + g
            d[notg] = []
            for st in self.genome_states:
                if (st[i] == '1'):
                    d[g].append(st)
                else:
                    d[notg].append(st)
        self.gene_dict_of_states = d

    def gether_states(self):
        """Make a table of genome states associated with each gene state

        Create new instance variable with the generated information
        """
        def get_states(d, func):
            big_set = set()
            for f in func:
                stage = f.split('&')
                s = get_intersected_states(d, stage)
                big_set = big_set.union(s)
            return (big_set)

        def get_intersected_states(d, func):
            s = set(d[func[0]])
            for f in func:
                s = s.intersection(set(d[f]))
            return s

        lst = []
        for i, g in enumerate(self.genes):
            f = self.func[i].split(' | ')
            s = get_states(self.gene_dict_of_states, f)
            lst.append(s)
        self.states_lst = lst

    def build_table(self):
        """Make a table of genome states transitions

        Create new instance variable with the generated information
        """
        table = {}
        for gs in self.genome_states:
            lst = []
            for i in range(len(self.states_lst)):
                lst.append(str(int(gs in self.states_lst[i])))
            table[gs] = tuple(lst)
        self.table = table

    def print_table(self):
        '''Print the Transition State Table in self object
        '''
        for key in self.table:
            s = ''.join(key)
            print("{} ".format(s), end='')
            lst = list(self.table.get(key))
            s = ''.join(lst)
            print(s)

    def get_attractor_basin(self):
        """Make a table of genome paths and make the atractor basins

        Create two new instance variables with the generated information
        """
        paths = list()
        for i, key in enumerate(self.table):
            paths.append(list())
            tupla = key
            while (tupla not in paths[i]):
                paths[i].append(tupla)
                tupla = self.table.get(tupla)
        basins = list()
        for i in range(len(paths)):
            tupla = paths[i][-1]
            pertence = 0
            for j, b in enumerate(basins):
                if tupla in b:
                    basins[j] = basins[j].union(set(paths[i]))
                    pertence = 1
                    break
            if pertence == 0:
                basins.append(set(paths[i]))
        table_of_paths = {}
        for i, key in enumerate(self.table):
            table_of_paths[key] = paths[i]
        self.basins, self.table_of_paths = (basins, table_of_paths)

    def get_attractors(self):
        """Find the atractor of each basin

        Create new instance variable with the generated information
        """
        attractors = list()
        for i, b in enumerate(self.basins):
            attractors.append(list())
            minimum = len(self.table_of_paths)
            for tupla in sorted(b):
                l = len(self.table_of_paths.get(tupla))
                if (l < minimum):
                    minimum = l
                    attractors[i] = [self.table_of_paths.get(tupla), len(b)]
        self.attractors = sorted(attractors)

    def print_attractors(self):
        """Print all attractors and the size of the respective basin
        """
        for a in self.attractors:
            lst = []
            for e in a[0]:
                s = ''.join(e)
                lst.append(s)
            s = "', '".join(lst)
            print("attractor: ['{}'] has basin size: {}".format(s,a[1]))

# in case one wants to provide direct information to the class BooleanNetwork
def read_data (FileBase):
    """Read file to get gene names and regulatory functions.

    FileBase -> name of the file.
    The content of the file must, in each line, the name of
    the gene followed by its function
    Exemple:
        v1, v4 | v2
        v2, !v4&v5 | v4&!v5 | v3&v5 | v3&v4
        v3, v1&!v3&!v5 | v3&v5 | !v1&v5 | !v1&v3
        v4, v4
        v5, v1&v4&v5
    """
    f = open(FileBase, "r")
    s = f.readlines()
    genes = list()
    func = list()
    for line in sorted(s):
        gene, fun = line.strip().split(', ')
        genes.append(gene.strip())
        func.append(fun.strip())
    f.close()
    return genes, func

file = input('filepath/filename: ').strip()
genes, func = read_data(file)
# bn = BooleanNetwork(genes, func)
bn = BooleanNetwork(file=(file))
bn.gene_combination()
# bn.print_gene_combination()
bn.gene_dict()
bn.gether_states()
bn.build_table()
# bn.print_table()
bn.get_attractor_basin()
bn.get_attractors()
bn.print_attractors()




