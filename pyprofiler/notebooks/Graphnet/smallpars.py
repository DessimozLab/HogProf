import itertools
import dendropy
import sys


####small parsimony functions ##########
def process_node_smallpars_1(node,allowed_symbols , transition_dict ):
    #go from leaves up and generate character sets
    if node.symbols is None:
        for child in node.child_nodes():
            if child.symbols is None:
                process_node_smallpars_1(child, allowed_symbols, transition_dict)
        node.symbols = { }
        node.scores = { }
        symbols = set.intersection( * [ child.symbols for child in node.child_nodes( ) ] )
        if len(symbols) == 0:
            symbols = set.union( * [ child.symbols for child in node.child_nodes( ) ] )
        node.symbols = symbols
        for c in allowed_symbols:
            if c not in node.symbols:
                #add trnasition mat here if needed for more subtle scoring
                score = min(  [ child.scores[c] for child in node.child_nodes()])+1
            else:
                score = min(  [ child.scores[c] for child in node.child_nodes() ] )
            node.scores[c] = score

def process_node_smallpars_2(node , allowed_symbols , transition_dict , verbose = False):
    #assign the most parsimonious char from children
    if node.char is None:
        if node.parent_node:
            #node has parent
            node.char = {}
            node.event = {}
            node.eventype= {}
            node.char = min(node.scores, key=node.scores.get)
            if node.parent_node.char == node.char:
                node.event = 0
            else:
                if node.scores[node.parent_node.char] == node.scores[node.char] :
                    node.char = node.parent_node.char
                    node.event = 0
                else:
                    node.event = 1
                    node.eventype = transition_dict[(node.parent_node.char,node.char)]
        else:
            #root node
            node.char = {}
            node.event= {}
            node.eventype = {}
            node.char = min(node.scores, key=node.scores.get)
            node.event = 0
        #down one level
        for child in node.child_nodes():
            if child.char is None:
                process_node_smallpars_2(child , allowed_symbols, transition_dict)

def calculate_small_parsimony(t, allowed_symbols, transition_dict):
    #upq
    process_node_smallpars_1(t.seed_node , allowed_symbols, transition_dict )
    #down
    process_node_smallpars_2(t.seed_node, allowed_symbols, transition_dict )
    return t
