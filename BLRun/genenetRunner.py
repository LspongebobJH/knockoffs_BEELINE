import os
import pandas as pd
from pathlib import Path
import numpy as np


def generateInputs(RunnerObj):
    '''
    Function to generate desired inputs for LOOK.
    If the folder/files under RunnerObj.datadir exist, 
    this function will not do anything.
    '''
    if not RunnerObj.inputDir.joinpath("GENENET").exists():
        print("Input folder for GENENET does not exist, creating input folder...")
        RunnerObj.inputDir.joinpath("GENENET").mkdir(exist_ok = False, parents = True)
        
    if not RunnerObj.inputDir.joinpath("GENENET/ExpressionData.csv").exists():
        ExpressionData = pd.read_csv(RunnerObj.inputDir.joinpath(RunnerObj.exprData),
                                     header = 0, index_col = 0)
        newExpressionData = ExpressionData.copy()
        # Write .csv file
        newExpressionData.to_csv(RunnerObj.inputDir.joinpath("GENENET/ExpressionData.csv"),
                             sep = ',', header  = True, index = True)
    

def run(RunnerObj):
    '''
    Function to run GENENET algorithm
    '''
    inputPath = "data" + str(RunnerObj.inputDir).split(str(Path.cwd()))[1] + \
                    "/ExpressionData.csv"
    # make output dirs if they do not exist:

    outDir = "outputs/"+str(RunnerObj.inputDir).split("inputs/")[1]+f"/{RunnerObj.name}/"
    os.makedirs(outDir, exist_ok = True)    
    outPath = "data/" +  str(outDir) + 'outFile.txt'
    # options
    calibrate = str(RunnerObj.params['calibrate'])
    data_mode = str(RunnerObj.params['data_mode'])
    cmdToRun = ' '.join(['docker run --rm -v', str(Path.cwd())+':/data/ genenet:base /bin/sh -c \"time -v -o', "data/" + str(outDir) + 'time.txt', 
    'Rscript runGeneNet.R',
    '-e', inputPath, 
    '-o', outPath,
    '-c', calibrate, 
    '-d', data_mode, 
    '\"'])
    print(cmdToRun)
    os.system(cmdToRun)


def parseOutput(RunnerObj):
    '''
    Function to parse outputs from GENENET.
    '''
    # Quit if output directory does not exist

    outDir = "outputs/"+str(RunnerObj.inputDir).split("inputs/")[1]+f"/{RunnerObj.name}/"
    if not Path(outDir+'outFile.txt').exists():
        print(outDir+'outFile.txt'+'does not exist, skipping...')
        return
        
    # Read output
    OutDF = pd.read_csv(outDir+'outFile.txt', sep = '\t', header = 0)
    outFile = open(outDir + 'rankedEdges.csv','w')
    outFile.write('Gene1'+'\t'+'Gene2'+'\t'+'EdgeWeight'+'\n')

    for idx, row in OutDF.sort_values('q_value', ascending = True).iterrows():
        outFile.write('\t'.join([row['Gene1'],row['Gene2'],str(-np.log10(row['q_value'] + 0.000000001))])+'\n')
    outFile.close()

