import os
import argparse
import numpy as np
import pandas as pd
import networkx as nx
from tqdm import tqdm
import multiprocessing
from pathlib import Path
import concurrent.futures
from collections import defaultdict
from itertools import product, permutations
from multiprocessing import Pool, cpu_count
from networkx.convert_matrix import from_pandas_adjacency
    
def FDR(evalObject, algorithmName, TFEdges = False):
    '''
    Computes empirical FDR for a given algorithm for each dataset.
    This is done for a range of targeted FDRs.
    We define FDR as the fraction of false 
    positives among those edges having q less than the targeted FDR,
    where q is an FDR-adjusted p-value.
    Only the column named "q_value" is accepted.
    
    
    :param evalObject: An object of class :class:`BLEval.BLEval`.
    :type evalObject: BLEval
      
    :param algorithmName: Name of the algorithm for which the early precision is computed.
    :type algorithmName: str
      
            
    :returns:
        A dataframe containing targeted and empirical FDR values
        for a given algorithm for each dataset.

    '''
    rankDict = {}
    for dataset in tqdm(evalObject.input_settings.datasets):
        trueEdgesDF = pd.read_csv(str(evalObject.input_settings.datadir)+'/'+ \
                      dataset['name'] + '/' +\
                      dataset['trueEdges'], sep = ',',
                      header = 0, index_col = None)
        trueEdgesDF = trueEdgesDF.loc[(trueEdgesDF['Gene1'] != trueEdgesDF['Gene2'])]
        trueEdgesDF.drop_duplicates(keep = 'first', inplace=True)
        trueEdgesDF.reset_index(drop=True, inplace=True)

        outDir = str(evalObject.output_settings.base_dir) + \
                 str(evalObject.input_settings.datadir).split("inputs")[1] + \
                 "/" + dataset["name"] + "/" + algorithmName

        #algos = evalObject.input_settings.algorithms
        rank_path = outDir + "/outFile.txt"
        if not os.path.isdir(outDir):
            rankDict[dataset["name"]] = set([])
            continue
        try:
            predDF = pd.read_csv(rank_path, sep="\t", header=0, index_col=None)
        except:
            print("\nSkipping fdr computation for ", algorithmName, "on path", outDir)
            rankDict[dataset["name"]] = set([])
            continue

        if not "q_value" in predDF.columns:
            print("\nSkipping fdr computation for ", algorithmName, "on path", outDir)
            rankDict[dataset["name"]] = set([])
            continue
            
        predDF = predDF.loc[(predDF['Gene1'] != predDF['Gene2'])]
        predDF.drop_duplicates(keep = 'first', inplace=True)
        predDF.reset_index(drop=True, inplace=True)
        
        if TFEdges:
            # Consider only edges going out of TFs
            
            # Get a list of all possible TF to gene interactions 
            uniqueNodes = np.unique(trueEdgesDF.loc[:,['Gene1','Gene2']])
            possibleEdges_TF = set(product(set(trueEdgesDF.Gene1),set(uniqueNodes)))

            # Get a list of all possible interactions 
            possibleEdges_noSelf = set(permutations(uniqueNodes, r = 2))
            
            # Find intersection of above lists to ignore self edges
            # TODO: is there a better way of doing this?
            possibleEdges = possibleEdges_TF.intersection(possibleEdges_noSelf)
            
            TrueEdgeDict = {'|'.join(p):0 for p in possibleEdges}

            trueEdges = trueEdgesDF['Gene1'] + "|" + trueEdgesDF['Gene2']
            trueEdges = trueEdges[trueEdges.isin(TrueEdgeDict)]
            print("\nEdges considered ", len(trueEdges))
            numEdges = len(trueEdges)
        
            predDF['Edges'] = predDF['Gene1'] + "|" + predDF['Gene2']
            # limit the predicted edges to the genes that are in the ground truth
            predDF = predDF[predDF['Edges'].isin(TrueEdgeDict)]

        else:
            trueEdges = trueEdgesDF['Gene1'] + "|" + trueEdgesDF['Gene2']
            trueEdges = set(trueEdges.values)
            numEdges = len(trueEdges)
        
        # check if ranked edges list is empty
        # if so, it is just set to an empty set
        if predDF.shape[0] == 0:
            print("\nSkipping fdr computation for on path ", rank_path,"due to lack of predictions.")
            rankDict[dataset["name"]] = set([])
        else:
            # we want to ensure that we do not include
            # edges without any edge weight
            # so check if the non-zero minimum is
            # greater than the edge weight of the top-kth
            # node, else use the non-zero minimum value.
            def get_candidates(threshold):
                newDF = predDF.loc[(predDF['q_value'] <= threshold)]
                return set(newDF['Gene1'] + "|" + newDF['Gene2'])
            rankDict[dataset["name"]] = [get_candidates(0.1*x) for x in range(10)]
    Fdr = {}

    for dataset in tqdm(evalObject.input_settings.datasets):
        if len(rankDict[dataset["name"]]) != 0:
            def compute_fdr(candidates):
                falseDiscoveries = candidates.difference(trueEdges)
                return len(falseDiscoveries)/max(1, len(candidates))
            Fdr[dataset["name"]] = [compute_fdr(c) for c in rankDict[dataset["name"]]]
        else:
            Fdr[dataset["name"]] = [np.NaN for i in range(10)]


    return(Fdr)
