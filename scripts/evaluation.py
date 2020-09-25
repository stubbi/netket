import matplotlib.pyplot as plt
import collections
import codecs
import pickle
import sys
import json
import pandas
import os
import shutil
import numpy as np
import math
import nqs as nq

# evaluation of a specific configuration

experimentFolder = sys.argv[1]
listSystemSizes = sys.argv[2].split(',')
listCycles = sys.argv[3].split(',')
numCircuits = int(sys.argv[4])
listOMPNodes = sys.argv[5].split(',')
listOMPTasks = sys.argv[6].split(',')
listOMPThreads = sys.argv[7].split(',')
listSamples = sys.argv[8].split(',')
listIterations = sys.argv[9].split(',')
listInitialHidden = sys.argv[10].split(',')
listSampleSteps = sys.argv[11].split(',')
numRuns = int(sys.argv[12])
randomRestarts = int(sys.argv[13])
earlyStopping = (str(sys.argv[14]) == 'True')
optimizer = str(sys.argv[15])

class Evaluation:
    def __init__(self, experimentFolder, listSystemSizes, listCycles, numCircuits, listOMPNodes, listOMPTasks, listOMPThreads, listSamples, listIterations, listInitialHidden, listSampleSteps, numRuns, randomRestarts, earlyStoppings, optimizer):
        self.experimentFolder=experimentFolder
        self.listSystemSizes=listSystemSizes
        self.listCycles=listCycles
        self.numCircuits=numCircuits
        self.listOMPNodes=listOMPNodes
        self.listOMPTasks=listOMPTasks
        self.listOMPThreads=listOMPThreads
        self.listSamples=listSamples
        self.listIterations=listIterations
        self.listInitialHidden=listInitialHidden
        self.listSampleSteps=listSampleSteps
        self.numRuns=numRuns
        self.randomRestarts=randomRestarts
        self.earlyStopping=earlyStopping
        self.optimizer=optimizer

    def plotTVDIterations(self):
        pass

    def plotTVDSamples(self):
        pass

    def plotTVDSamplesIterations(self):
        pass

    def loadHistogram(self, size, cycles, circuits, nodes, tasks, threads, numSamples, numIterations, numInitialHidden, numSampleSteps, run):
        with open("{directory}/histogram.json".format(directory=self.directory(size, cycles, circuits, nodes, tasks, threads, numSamples, numIterations, numInitialHidden, numSampleSteps, run)), 'r') as f:
            return json.load(f)

    def loadDuration(self, size, cycles, circuits, nodes, tasks, threads, numSamples, numIterations, numInitialHidden, numSampleSteps, run):
        with open("{directory}/duration.time".format(directory=self.directory(size, cycles, circuits, nodes, tasks, threads, numSamples, numIterations, numInitialHidden, numSampleSteps, run)), 'r') as f:
            return float(f.readline())

    def loadExact(self,qubits, cycles, circuit, gateNo = -1):
        if(gateNo == -1):
            with open("{directory}/{qubits}qubits/{cycles}cycles/circuit{circuit}/exact.json".format(directory=self.experimentFolder,qubits=qubits, cycles=cycles, circuit=circuit), "rb") as f:
                return pickle.loads(f.read())
        else:
            with open("{directory}/{qubits}qubits/{cycles}cycles/circuit{circuit}/exact_gate_{gateNo}.json".format(directory=self.experimentFolder,qubits=qubits, cycles=cycles, circuit=circuit, gateNo=gateNo+1), "rb") as f:
                return pickle.loads(f.read())
    
    def loadRBM(self, size, cycles, circuits, nodes, tasks, threads, numSamples, numIterations, numInitialHidden, numSampleSteps, run, gateNo = -1):
        f = ''
        if(gateNo == -1):
            f = "{directory}/parameters.json".format(directory=self.directory(size, cycles, circuits, nodes, tasks, threads, numSamples, numIterations, numInitialHidden, numSampleSteps, run))
        else:
            f = "{directory}/parameters_gate_{gateNo}.json".format(directory=self.directory(size, cycles, circuits, nodes, tasks, threads, numSamples, numIterations, numInitialHidden, numSampleSteps, run), gateNo=gateNo+1)
        nqs = nq.nqs.NQS(int(size),int(numInitialHidden),int(numSamples), self.randomRestarts, self.earlyStopping, self.optimizer)
        nqs.load(f)
        return nqs

    def loadLogOverlap(self, size, cycles, circuits, nodes, tasks, threads, numSamples, numIterations, numInitialHidden, numSampleSteps, run, gateNo, testOrTraining='test'):
        with codecs.open("{directory}/{gateNo}.log".format(directory=self.directory(size, cycles, circuits, nodes, tasks, threads, numSamples, numIterations, numInitialHidden, numSampleSteps, run), gateNo=gateNo), 'r', encoding='ascii') as f:
            data = json.load(f)
            return [d["{}_log_overlap".format(testOrTraining)] for d in data['Output']]

    def numberOfHadamards(self, qubits, cycles, circuit):
        hadamards = 0
        with open("{directory}/{qubits}qubits/{cycles}cycles/circuit{circuit}/in.qc".format(directory=self.experimentFolder,qubits=qubits, cycles=cycles, circuit=circuit), "r") as f:
            for line in f:
                if(line[0] == 'H'):
                    hadamards = hadamards + 1
        return hadamards

    def loadRBMProbs(self, exact, rbm):
        def toBinaryArray(i):
            return [int(b) for b in format(i, '0{}b'.format(int(math.log(len(exact),2))))]

        nqs_probs = [abs(rbm.psi(toBinaryArray(i)))**2 for i in range(len(exact))]
        nqs_norm = np.sum(nqs_probs)
        return [n/nqs_norm for n in nqs_probs]

    def loadRBMAmplitudes(self, exact, rbm):
        def toBinaryArray(i):
            return [int(b) for b in format(i, '0{}b'.format(int(math.log(len(exact),2))))]

        return [rbm.psi(toBinaryArray(i)) for i in range(len(exact))]

    def tvd(self, exact, rbm):
        exact_probs = [abs(e)**2 for e in exact]
        nqs_probs = self.loadRBMProbs(exact, rbm)
        return np.sum([abs(exact_probs[i] - nqs_probs[i]) for i in range(len(exact))])/2.0

    def f_xeb(self, exact, histogram, qubits):
        f_xeb = 0.0
        shots = 0
        for h in histogram:
            shots += histogram[h]
            f_xeb += histogram[h] * (2**qubits * abs(exact[int(h)])**2 + 1)
        return f_xeb/shots

    def plotAvgBestPDF(self, qubits):
        results_file = "{directory}/results.csv".format(directory=self.experimentFolder)
        dfo = pandas.read_csv(results_file)
        dfo = dfo[dfo['success'] == True]
        dfo = dfo.astype({'tvd': 'float64', 'duration': 'float64', 'f_xeb': 'float64'})

        fig, axs = plt.subplots(2,2)
        for c in self.listCycles:
            allSortedRbmProbs = []
            allSortedExactProbs = []
            for i in range(self.numCircuits):
                df = dfo.copy()
                try:
                    exact = self.loadExact(qubits, c, i, -1)
                    df = df[(df['#qubits'] == int(qubits)) & (df['#cycles'] == int(c)) & (df['circuit'] == int(i))]

                    row = df[df.tvd == df.tvd.min()]
                    rbm = self.loadRBM(row['#qubits'].iloc[0],row['#cycles'].iloc[0],row['circuit'].iloc[0],row['#nodes'].iloc[0],row['#tasks'].iloc[0],row['#threads'].iloc[0],row['#samples'].iloc[0],row['#iterations'].iloc[0],row['#initialHidden'].iloc[0],row['#sampleSteps'].iloc[0],row['run'].iloc[0], -1)

                    rbmProbs = [len(exact) * p for p in self.loadRBMProbs(exact, rbm)]
                    normalisedProbs = [len(exact) * abs(e)**2 for e in exact] 

                    rbmProbsSorted = [p for _,p in sorted(zip(normalisedProbs, rbmProbs))]
                    exactProbsSorted = sorted(normalisedProbs)

                    allSortedRbmProbs.append(rbmProbsSorted)
                    allSortedExactProbs.append(exactProbsSorted)

                except Exception as e:
                    print(e)
            
            idx0=0
            idx1=0
            if(int(c) == 10):
                idx1=1
            elif(int(c) == 15):
                idx0=1
            elif(int(c) == 20):
                idx0=1
                idx1=1

            meanRBM = np.mean(allSortedRbmProbs, axis=0)
            meanExact = np.mean(allSortedExactProbs, axis=0)

            axs[idx0,idx1].plot(range(len(meanExact)), meanExact, label = 'exact')
            axs[idx0,idx1].plot(range(len(meanRBM)), meanRBM, label = 'RBM')

            axs[idx0,idx1].set_title('{c} Cycles'.format(c=c))
            axs[idx0,idx1].legend()
            axs[idx0,idx1].set_xlim([0, 15])
            axs[idx0,idx1].set_ylim([0, 4])

        for ax in axs.flat:
            ax.set(xlabel=r'Bit-string index $j$ ($p(x_j)$-ordered)', ylabel=r'$Np$')

        plt.tight_layout()
        plt.savefig('avgBestPDF.pdf')
        plt.close()


    def plotAvgPDF(self, qubits):
        results_file = "{directory}/results.csv".format(directory=self.experimentFolder)
        dfo = pandas.read_csv(results_file)
        dfo = dfo[dfo['success'] == True]
        dfo = dfo.astype({'tvd': 'float64', 'duration': 'float64', 'f_xeb': 'float64'})

        fig, axs = plt.subplots(2,2)
        for c in self.listCycles:
            allSortedRbmProbs = []
            allSortedExactProbs = []
            for i in range(self.numCircuits):
                df = dfo.copy()
                df = df[(df['#qubits'] == int(qubits)) & (df['#cycles'] == int(c)) & (df['circuit'] == int(i))]

                for index, row in df.iterrows():
                    try:
                        exact = self.loadExact(qubits, c, i, -1)
                        rbm = self.loadRBM(row['#qubits'],row['#cycles'],row['circuit'],row['#nodes'],row['#tasks'],row['#threads'],row['#samples'],row['#iterations'],row['#initialHidden'],row['#sampleSteps'],row['run'], -1)

                        rbmProbs = [len(exact) * p for p in self.loadRBMProbs(exact, rbm)]
                        normalisedProbs = [len(exact) * abs(e)**2 for e in exact] 

                        rbmProbsSorted = [p for _,p in sorted(zip(normalisedProbs, rbmProbs))]
                        exactProbsSorted = sorted(normalisedProbs)

                        allSortedRbmProbs.append(rbmProbsSorted)
                        allSortedExactProbs.append(exactProbsSorted)

                    except Exception as e:
                        print(e)
            
            idx0=0
            idx1=0
            if(int(c) == 10):
                idx1=1
            elif(int(c) == 15):
                idx0=1
            elif(int(c) == 20):
                idx0=1
                idx1=1

            meanRBM = np.mean(allSortedRbmProbs, axis=0)
            meanExact = np.mean(allSortedExactProbs, axis=0)

            axs[idx0,idx1].plot(range(len(meanExact)), meanExact, label = 'exact')
            axs[idx0,idx1].plot(range(len(meanRBM)), meanRBM, label = 'RBM')

            axs[idx0,idx1].set_title('{c} Cycles'.format(c=c))
            axs[idx0,idx1].legend()
            axs[idx0,idx1].set_xlim([0, 15])
            axs[idx0,idx1].set_ylim([0, 4])

        for ax in axs.flat:
            ax.set(xlabel=r'Bit-string index $j$ ($p(x_j)$-ordered)', ylabel=r'$Np$')

        plt.tight_layout()
        plt.savefig('avgPDF.pdf')
        plt.close()

    def plotAvgLogOverlap(self, qubits):
        maxIterations = 10000
        cutIterations = 2000
        for s in [8]:#self.listSamples:
            sqrt_XTrain = []
            sqrt_XTest = []

            sqrt_YTrain = []
            sqrt_YTest = []

            CZTrain = []
            CZTest = []

            allTrain = []
            allTest = []

            for c in [5]:#self.listCycles:
                for i in [1,2,3,4]:#range(self.numCircuits):

                    with open("{}/{}qubits/{}cycles/circuit{}/in.qc".format(self.experimentFolder,qubits,c,i)) as f:
                        content = f.readlines()

                        # TODO only works for current RCSs!
                        content = [x for x in content if x.startswith('T') or x.startswith('sqrt_X') or x.startswith('sqrt_Y') or x.startswith('CZ')] 

                        idx = 0
                        for line in content:
                            gate = line.split()[0]
                            if(gate != 'T'): #TODO: one method does not learn CZ!
                                idx = idx + 1
                                for r in range(self.numRuns):
                                    #TODO: numInitialHidden is not always 6!
                                    try:
                                        testLogOverlaps = self.loadLogOverlap(qubits,c,i,1,1,1,s,maxIterations,6,0,r, idx, 'test')
                                        trainLogOverlaps = self.loadLogOverlap(qubits,c,i,1,1,1,s,maxIterations,6,0,r, idx, 'training')
                                        
                                        allTest.append(testLogOverlaps[:cutIterations])
                                        allTrain.append(trainLogOverlaps[:cutIterations])

                                        if(gate == 'sqrt_X'):
                                            sqrt_XTrain.append(trainLogOverlaps[:cutIterations])
                                            sqrt_XTest.append(testLogOverlaps[:cutIterations])
                                        elif(gate == 'sqrt_Y'):
                                            sqrt_YTrain.append(trainLogOverlaps[:cutIterations])
                                            sqrt_YTest.append(testLogOverlaps[:cutIterations])
                                        elif(gate == 'CZ'):
                                            CZTrain.append(trainLogOverlaps[:cutIterations])
                                            CZTest.append(testLogOverlaps[:cutIterations])

                                    except Exception as e:
                                        print(e)
            
            fig, axs = plt.subplots(2,2)
            meanAllTrain = np.mean(allTrain, axis=0)
            meanAllTest = np.mean(allTest, axis=0)
            stdAllTrain = np.std(allTrain, axis=0)
            stdAllTest = np.std(allTest, axis=0)
            ciAllTrain = 1.96 * stdAllTrain/meanAllTrain
            ciAllTest = 1.96 * stdAllTest/meanAllTest

            meanSqrtXTrain = np.mean(sqrt_XTrain, axis=0)
            meanSqrtXTest = np.mean(sqrt_XTest, axis=0)
            stdSqrtXTrain = np.std(sqrt_XTrain, axis=0)
            stdSqrtXTest = np.std(sqrt_XTest, axis=0)
            ciSqrtXTrain = 1.96 * stdSqrtXTrain/meanSqrtXTrain
            ciSqrtXTest = 1.96 * stdSqrtXTest/meanSqrtXTest

            meanSqrtYTrain = np.mean(sqrt_YTrain, axis=0)
            meanSqrtYTest = np.mean(sqrt_YTest, axis=0)
            stdSqrtYTrain = np.std(sqrt_YTrain, axis=0)
            stdSqrtYTest = np.std(sqrt_YTest, axis=0)
            ciSqrtYTrain = 1.96 * stdSqrtYTrain/meanSqrtYTrain
            ciSqrtYTest = 1.96 * stdSqrtYTest/meanSqrtYTest

            meanCZTrain = np.mean(CZTrain, axis=0)
            meanCZTest = np.mean(CZTest, axis=0)
            stdCZTrain = np.std(CZTrain, axis=0)
            stdCZTest = np.std(CZTest, axis=0)
            ciCZTrain = 1.96 * stdCZTrain/meanCZTrain
            ciCZTest = 1.96 * stdCZTest/meanCZTest

            axs[0,0].plot(range(len(meanSqrtXTest)), meanSqrtXTest, label = 'test')
            axs[0,0].plot(range(len(meanSqrtXTrain)), meanSqrtXTrain, label = 'train')
            axs[0,0].fill_between(range(len(meanSqrtXTest)), (np.array(meanSqrtXTest)-ciSqrtXTest), (np.array(meanSqrtXTest)+ciSqrtXTest), alpha=.1)
            axs[0,0].fill_between(range(len(meanSqrtXTrain)), (np.array(meanSqrtXTrain)-ciSqrtXTrain), (np.array(meanSqrtXTrain)+ciSqrtXTrain), alpha=.1)

            axs[0,0].set_title('sqrt X')
            axs[0,0].legend()

            axs[0,1].plot(range(len(meanSqrtYTest)), meanSqrtYTest, label = 'test')
            axs[0,1].plot(range(len(meanSqrtYTrain)), meanSqrtYTrain, label = 'train')
            axs[0,1].fill_between(range(len(meanSqrtYTest)), (np.array(meanSqrtYTest)-ciSqrtYTest), (np.array(meanSqrtYTest)+ciSqrtYTest), alpha=.1)
            axs[0,1].fill_between(range(len(meanSqrtYTrain)), (np.array(meanSqrtYTrain)-ciSqrtYTrain), (np.array(meanSqrtYTrain)+ciSqrtYTrain), alpha=.1)

            axs[0,1].set_title('sqrt Y')
            axs[0,1].legend()

            axs[1,0].plot(range(len(meanCZTest)), meanCZTest, label = 'test')
            axs[1,0].plot(range(len(meanCZTrain)), meanCZTrain, label = 'train')
            axs[1,0].fill_between(range(len(meanCZTest)), (np.array(meanCZTest)-ciCZTest), (np.array(meanCZTest)+ciCZTest), alpha=.1)
            axs[1,0].fill_between(range(len(meanCZTrain)), (np.array(meanCZTrain)-ciCZTrain), (np.array(meanCZTrain)+ciCZTrain), alpha=.1)

            axs[1,0].set_title('CZ')
            axs[1,0].legend()

            axs[1,1].plot(range(len(meanAllTest)), meanAllTest, label = 'test')
            axs[1,1].plot(range(len(meanAllTrain)), meanAllTrain, label = 'train')
            axs[1,1].fill_between(range(len(meanAllTest)), (np.array(meanAllTest)-ciAllTest), (np.array(meanAllTest)+ciAllTest), alpha=.1)
            axs[1,1].fill_between(range(len(meanAllTrain)), (np.array(meanAllTrain)-ciAllTrain), (np.array(meanAllTrain)+ciAllTrain), alpha=.1)

            axs[1,1].set_title('All Gates')
            axs[1,1].legend()

            for ax in axs.flat:
                ax.set(xlabel='Iteration', ylabel='log Overlap')

            plt.tight_layout()
            plt.savefig("avgOverlap_{s}.pdf".format(s=s))
            plt.close()






    def plotPDF(self, df, qubits, cycles, circuit, circuitFile, gateNo = -1):
        try:
            exact = self.loadExact(qubits, cycles, circuit, gateNo)
            df = df[(df['#qubits'] == int(qubits)) & (df['#cycles'] == int(cycles)) & (df['circuit'] == int(circuit))]
            minRow = df[df.tvd == df.tvd.min()]
            maxRow = df[df.tvd == df.tvd.max()]

            minRBM = self.loadRBM(minRow['#qubits'].iloc[0],minRow['#cycles'].iloc[0],minRow['circuit'].iloc[0],minRow['#nodes'].iloc[0],minRow['#tasks'].iloc[0],minRow['#threads'].iloc[0],minRow['#samples'].iloc[0],minRow['#iterations'].iloc[0],minRow['#initialHidden'].iloc[0],minRow['#sampleSteps'].iloc[0],minRow['run'].iloc[0], gateNo)
            maxRBM = self.loadRBM(maxRow['#qubits'].iloc[0],maxRow['#cycles'].iloc[0],maxRow['circuit'].iloc[0],maxRow['#nodes'].iloc[0],maxRow['#tasks'].iloc[0],maxRow['#threads'].iloc[0],maxRow['#samples'].iloc[0],maxRow['#iterations'].iloc[0],maxRow['#initialHidden'].iloc[0],maxRow['#sampleSteps'].iloc[0],maxRow['run'].iloc[0], gateNo)

            if(gateNo > 0):
                testLogOverlaps = self.loadLogOverlap(minRow['#qubits'].iloc[0],minRow['#cycles'].iloc[0],minRow['circuit'].iloc[0],minRow['#nodes'].iloc[0],minRow['#tasks'].iloc[0],minRow['#threads'].iloc[0],minRow['#samples'].iloc[0],minRow['#iterations'].iloc[0],minRow['#initialHidden'].iloc[0],minRow['#sampleSteps'].iloc[0],minRow['run'].iloc[0], gateNo, 'test')
                trainingLogOverlaps = self.loadLogOverlap(minRow['#qubits'].iloc[0],minRow['#cycles'].iloc[0],minRow['circuit'].iloc[0],minRow['#nodes'].iloc[0],minRow['#tasks'].iloc[0],minRow['#threads'].iloc[0],minRow['#samples'].iloc[0],minRow['#iterations'].iloc[0],minRow['#initialHidden'].iloc[0],minRow['#sampleSteps'].iloc[0],minRow['run'].iloc[0], gateNo, 'training')
                fig, ax = plt.subplots()
                ax.plot(range(len(testLogOverlaps)), testLogOverlaps, label='test')
                ax.plot(range(len(trainingLogOverlaps)), trainingLogOverlaps, label='training')
                plt.legend()
                plt.suptitle(self.experiment(), fontsize=14, fontweight='bold')
                plt.ylabel('log Overlap')
                plt.xlabel('Iteration')
                plt.title('{} qubits {} cycles circuit {}'.format(qubits, cycles, circuit), fontdict={'size':10})
                plt.savefig('plots/circuits/gatewise/logOverlap_{}qubits_{}cycles_circuit{}_gate{}.pdf'.format(qubits, cycles, circuit, gateNo))
                plt.close()

            bestRBMProbs = [len(exact) * p for p in self.loadRBMProbs(exact, minRBM)]
            worstRBMProbs = [len(exact) * p for p in self.loadRBMProbs(exact, maxRBM)]
            normalisedProbs = [len(exact) * abs(e)**2 for e in exact] 

            bestRBMProbsSorted = [p for _,p in sorted(zip(normalisedProbs, bestRBMProbs))]
            worstRBMProbsSorted = [p for _,p in sorted(zip(normalisedProbs, worstRBMProbs))]
            exactProbsSorted = sorted(normalisedProbs)

            fig, ax = plt.subplots()
            ax.plot(range(len(exactProbsSorted)), exactProbsSorted, label = 'exact')
            ax.plot(range(len(exactProbsSorted)), bestRBMProbsSorted, label = 'best rbm, {} samples {} iterations {} sample steps, tvd: {:.2f}'.format(minRow.iloc[0]['#samples'], minRow.iloc[0]['#iterations'], minRow.iloc[0]['#sampleSteps'], self.tvd(exact,minRBM)))
            plt.legend()
            plt.suptitle(self.experiment(), fontsize=14, fontweight='bold')
            plt.ylabel('p(j)')
            plt.xlabel('Bit string index j (ordered)')
            if(gateNo != -1):
                plt.title('{} qubits {} cycles circuit {} entropy: {:.2f} Porter-Thomas: {:.2f} gate: {} ({})'.format(qubits, cycles, circuit, self.circuitEntropy(qubits, cycles, circuit), self.porterThomasEntropy(qubits), gateNo, circuitFile[gateNo-1]), fontdict={'size':10})
                plt.savefig('plots/circuits/gatewise/pdf_{}qubits_{}cycles_circuit{}_gate{}.pdf'.format(qubits, cycles, circuit, gateNo))
            else:
                plt.title('{} qubits {} cycles circuit {} entropy: {:.2f} Porter-Thomas: {:.2f} '.format(qubits, cycles, circuit, self.circuitEntropy(qubits, cycles, circuit), self.porterThomasEntropy(qubits)), fontdict={'size':10})
                plt.savefig('plots/circuits/pdf_{}qubits_{}cycles_circuit{}.pdf'.format(qubits, cycles, circuit))
            plt.close()
        except Exception as e:
            print(e)
            print('no plot for {} qubits {} cycles circuit {} gateNo {}'.format(qubits, cycles, circuit, gateNo))

    def plotDepthEntropy(self, qubits):
        entropies = [np.average([self.circuitEntropy(qubits, cycles, circuit) for circuit in range(self.numCircuits)]) for cycles in self.listCycles]

        fig, ax = plt.subplots()
        ax.plot(self.listCycles, entropies)
        ax.axhline(self.porterThomasEntropy(qubits))
        plt.suptitle(self.experiment(), fontsize=14, fontweight='bold')
        plt.title('{} qubits'.format(qubits), fontdict={'size':10})
        plt.ylabel('Entropy')
        plt.xlabel('Depht')
        plt.savefig('plots/circuits/entropy_{}qubits.pdf'.format(qubits))
        plt.close()


    def circuitEntropy(self, qubits, cycles, circuit):
        try:
            exact = self.loadExact(qubits, cycles, circuit)
            exact_probs = [abs(e)**2 for e in exact]
            return - np.sum([e * math.log(e) if e != 0 else 0 for e in exact_probs])
        except:
            print("no circuit entropy for {} qubits {} cycles circuit {}".format(qubits, cycles, circuit))
            return 0

    def porterThomasEntropy(self, qubits):
        return math.log(2**int(qubits)) - 1.0 + 0.577

    def experiment(self):
        return '{}\n'.format(self.experimentFolder.split('/')[-1])

    def plot(self, df, x, grouped, fixed):
        all_keys = ['#qubits', '#cycles', '#iterations', '#samples', '#sampleSteps', 'initialHidden', '#hadamards'] 
        experiment = self.experiment()
        title = ', '.join(['{}:{}'.format(key, fixed.get(key, 'all')) for key in all_keys if key not in grouped + [x]])

        filterBy = df[grouped].drop_duplicates()

        groupDFBy = grouped + [x]
        df = df.loc[(df[list(fixed)] == pandas.Series(fixed)).all(axis=1)]
        df = df.groupby(groupDFBy, as_index = False).mean() 

        for y in ['tvd', 'duration', 'f_xeb']:
            fig, ax = plt.subplots()
            name = '{}_{}_{}'.format(x.replace('#', ''), y.replace('#', ''), title.replace(', ', '_').replace('#', ''))
            for index, row in filterBy.iterrows(): 
                toPlot = pandas.merge(df, row.to_frame().T, how='inner')
                if(toPlot.shape[0] > 0):
                    l = ''.join(['{}:{} '.format(i, toPlot[i].tolist()[0]) for i in grouped])
                    ax.plot(toPlot[x].tolist(), toPlot[y].tolist(), label = l)
            plt.legend()
            plt.suptitle(experiment, fontsize=14, fontweight='bold')
            plt.title(title, fontdict={'size':10})
            plt.xlabel(x)
            plt.ylabel(y)
            plt.savefig('plots/{}.pdf'.format(name))
            plt.close()

    def generateReport(self):
        results_file = "{directory}/results.csv".format(directory=self.experimentFolder)
        df = pandas.read_csv(results_file)

        report_file = "{directory}/report.txt".format(directory=self.experimentFolder)
        with open(report_file, 'w') as f:
            f.write('experiment: {}\n\n'.format(self.experimentFolder.split('/')[-1]))
            f.write('total number of experiments: {}\n'.format(df.shape[0]))
            f.write('number of succeeded experiments: {}\n'.format(df[df['success'] == True].shape[0]))
            f.write('number of failed experiments: {}\n\n'.format(df[df['success'] == False].shape[0]))
            f.write('tested number of qubits: {}\n'.format(self.listSystemSizes))
            f.write('tested number of cycles: {}\n'.format(self.listCycles))
            f.write('tested number of samples: {}\n'.format(self.listSamples))
            f.write('tested number of iterations: {}\n'.format(self.listIterations))
            f.write('tested number of sample steps: {}\n'.format(self.listSampleSteps))
            f.write('number of runs: {}\n\n'.format(self.numRuns))
            f.write('failed combinations:\n')
            failed = df[df['success'] == False].copy().drop(columns=['tvd', 'duration', 'f_xeb']).drop_duplicates()
            for index, row in failed.iterrows():
                f.write('{} \n'.format(row))


    def generatePlots(self):
        results_file = "{directory}/results.csv".format(directory=self.experimentFolder)
        if not os.path.exists('plots'):
            os.makedirs('plots')
        if not os.path.exists('plots/circuits'):
            os.makedirs('plots/circuits')
        if not os.path.exists('plots/circuits/gatewise'):
            os.makedirs('plots/circuits/gatewise')
        df = pandas.read_csv(results_file)
        df = df[df['success'] == True]
        df = df.astype({'tvd': 'float64', 'duration': 'float64', 'f_xeb': 'float64'})

        for i in pandas.unique(df['#iterations']):
            for s in pandas.unique(df['#samples']):
                for x in ['#qubits', '#cycles']:
                    grouped = []
                    if(x == '#qubits'):
                        grouped = ['#cycles']
                    if(x == '#cycles'):
                        grouped = ['#qubits']
                    
                    fixed = {'#iterations': i,
                                '#samples': s}

                    self.plot(df.copy(), x, grouped, fixed)

        for i in pandas.unique(df['#iterations']):
            self.plot(df.copy(), '#samples', ['#qubits', '#cycles'], {'#iterations': i})

        for s in pandas.unique(df['#samples']):
            self.plot(df.copy(), '#iterations', ['#qubits', '#cycles'], {'#samples': s})

        self.plot(df.copy(), '#sampleSteps', ['#qubits', '#cycles'], {})
        self.plot(df.copy(), '#iterations', ['#qubits', '#cycles'], {})
        self.plot(df.copy(), '#iterations', ['#qubits', '#sampleSteps'], {})
        self.plot(df.copy(), '#sampleSteps', ['#iterations'], {})
        self.plot(df.copy(), '#iterations', ['#sampleSteps'], {})
        self.plot(df.copy(), '#samples', ['#qubits', '#cycles'], {})
        self.plot(df.copy(), '#qubits', ['#cycles'], {})
        self.plot(df.copy(), '#cycles', ['#qubits'], {})
        for i in [1,2]:
            for j in [1,2]:
                self.plot(df.copy(), '#nodes', ['#qubits', '#cycles'], {'#tasks': i, '#threads': j})
                self.plot(df.copy(), '#tasks', ['#qubits', '#cycles'], {'#nodes': i, '#threads': j})
                self.plot(df.copy(), '#threads', ['#qubits', '#cycles'], {'#nodes': i, '#tasks': j})

        for q in self.listSystemSizes:
            self.plotDepthEntropy(q)
            for c in self.listCycles:
                for i in range(self.numCircuits):
                    with open("{}/{}qubits/{}cycles/circuit{}/in.qc".format(self.experimentFolder,q,c,i)) as f:
                        content = f.readlines()

                    # TODO only works for current RCSs!
                    content = [x for x in content if x.startswith('T') or x.startswith('sqrt_X') or x.startswith('sqrt_Y') or x.startswith('CZ')] 

                    gates = len(content)

                    for g in range(gates + 2):
                        self.plotPDF(df.copy(), q, c, i, content, g-1)

        shutil.make_archive('plots', 'zip', 'plots')

          
    def generateCSV(self):
        results_file = "{directory}/results.csv".format(directory=self.experimentFolder)
        with open(results_file, 'w') as f:
            f.write('#qubits,#cycles,circuit,#nodes,#tasks,#threads,#samples,#iterations,#initialHidden,#sampleSteps,run,#hadamards,tvd,f_xeb,duration,success\n')
        
        for size in self.listSystemSizes:
            for cycles in self.listCycles:
                for circuits in range(self.numCircuits):
                    hadamards = self.numberOfHadamards(size, cycles, circuits)
                    for nodes in self.listOMPNodes:
                        for tasks in self.listOMPTasks:
                            for threads in self.listOMPThreads:
                                for numSamples in self.listSamples:
                                    for numIterations in self.listIterations:
                                        for numInitialHidden in self.listInitialHidden:
                                            step_size = int(len(self.listSampleSteps)/len(self.listSystemSizes))
                                            index_sample_steps = self.listSystemSizes.index(size) * step_size
                                            for numSampleSteps in self.listSampleSteps[index_sample_steps:index_sample_steps+step_size]:
                                                for run in range(self.numRuns):

                                                    try:
                                                        histogram = self.loadHistogram(size, cycles, circuits, nodes, tasks, threads, numSamples, numIterations, numInitialHidden, numSampleSteps, run)
                                                        rbm = self.loadRBM(size, cycles, circuits, nodes, tasks, threads, numSamples, numIterations, numInitialHidden, numSampleSteps, run)
                                                        tvd = '{:f}'.format(self.tvd(self.loadExact(size, cycles, circuits), rbm))
                                                        f_xeb = '{:f}'.format(self.f_xeb(self.loadExact(size, cycles, circuits), self.normalise(histogram), int(size)))
                                                        duration = '{:f}'.format(self.loadDuration(size, cycles, circuits, nodes, tasks, threads, numSamples, numIterations, numInitialHidden, numSampleSteps, run))
                                                        success = True
                                                    except Exception as e:
                                                        print("Unexpected error: {}".format(str(e)))
                                                        tvd = '-'
                                                        f_xeb = '-'
                                                        duration = '-'
                                                        success = False

                                                    line = "{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}\n".format(size,cycles,circuits,nodes,tasks,threads,numSamples,numIterations,numInitialHidden,numSampleSteps,run,hadamards,tvd,f_xeb,duration,success)
                                                    with open(results_file, 'a') as f:
                                                        f.write(line)

    def directory(self, qubits, cycles, circuit, nodes, tasks, threads, numSamples, numIterations, numInitialHidden, numSampleSteps, run):
        return "{experimentFolder}/{qubits}qubits/{cycles}cycles/circuit{circuit}/{nodes}nodes/{tasks}tasks/{threads}threads/{samples}samples/{iterations}iterations/{initial_hidden}initialHidden/{sample_steps}sampleSteps/run{run}".format(
                experimentFolder=self.experimentFolder,
                nodes=nodes,
                tasks=tasks,
                threads=threads,
                samples=numSamples,
                iterations=numIterations,
                initial_hidden=numInitialHidden,
                sample_steps=numSampleSteps,
                run=run,
                qubits=qubits,
                cycles=cycles,
                circuit=circuit
            )

    def normalise(self, histogram):
        total = sum(histogram.values())
        for key in histogram:
            histogram[key] = float(histogram[key])/float(total)
        return histogram

ev = Evaluation(experimentFolder, listSystemSizes, listCycles, numCircuits, listOMPNodes, listOMPTasks, listOMPThreads, listSamples, listIterations, listInitialHidden, listSampleSteps, numRuns, randomRestarts, earlyStopping, optimizer)
#ev.generateCSV()
#ev.generateReport()
#ev.generatePlots()
ev.plotAvgPDF(4)
ev.plotAvgBestPDF(4)
ev.plotAvgLogOverlap(4)
