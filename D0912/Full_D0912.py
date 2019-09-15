def GetNodeColumnsValues(columnNumber, nodesObject):
    currentDerivs = nodesObject.getDerivs()
    columnValues = []
    for node in range( len(currentDerivs) ):
        columnValues.append( currentDerivs[node][columnNumber][0] )
    return columnValues

def SetNodeColumnValues(columnValues, columnNumber, nodesObject):
    currentDerivs = nodesObject.getDerivs()
    for node in range( len(currentDerivs) ):
        currentDerivs[node][columnNumber][0] = columnValues[node]
        
def GetNodeLabels(nodesObject):
    numNodes = nodesObject.nodes[-1].number
    currentLabels = [ "default" for i in range(numNodes) ]
    currentLabels = nodesObject.getNodeLabels()
    return currentLabels

def SetNodeLabels(labels, nodesObject):
    numNodes = nodesObject.nodes[-1].number
    for i in range(numNodes):
        nodesObject.setNodeLabel(i+1, labels[i])

def NonLinSolve(cont, initialTime, timeStepIncrement, numSteps):
    import numpy
    import sys
    parallel = 1
    slu_dist = 0
    cacheParameters = False
    solutionIteration = 100
    errorMaxSumSol  = 1.0e-3
    errorMaxResID   = 1.0e-3
    finalTime = initialTime + numSteps*timeStepIncrement
    cont.Snonlin([numSteps, solutionIteration, timeStepIncrement, initialTime, finalTime, errorMaxSumSol, errorMaxResID, 1, {'delta':1e-6, 'cache_params': cacheParameters, 'info': 'gui_biomechanics_sol', 'sol_tol': 10, 'linear_solver': 0, 'update_param': 0.0, 'solution_output': 1, 'parallel': parallel, 'parallelLinearSolver':slu_dist, 'Krylovi_subspace': 10, 'trans_routine': 0, 'min_i_eigen': 0.0, 'max_r_eigen': 0.0, 'max_i_eigen': 0.0, 'solver_iterations': 200, 'vector_output': 3, 'abort': 1, 'line_search': 1, 'stopping': 1, 'min_r_eigen': 0.0, 'additibe_const': 0.0, 'Newton_Raphson': 0, 'preconditioning': 0, 'equilibrium': 2}], log=0)
    zp = cont.server.model.bmprob.deform_global_param_ZP
    if  numpy.isnan(zp).any():
        print "\n Solve Failed \n"      
        sys.exit(1)

import os
import sys

import csv

# LV = 0, BiV = 1
biVModel = 0
usingGPT = False

pID = str(sys.argv[1])

print(pID)

subject = 'D0912'
Vref = 29.95

baseFileName = subject+'_'+pID+'_ED_DT_'
dirName = 'Set this to working directory/'+subject+'/'
continuityDir = 'Continuity install directory'
inflationOutputDir = dirName+'Inflation_'+pID+'_'
deflationOutputDir = dirName+'Deflation_'+pID+'_'
hemoDataFileName = 'Simulation'+subject+'_'+pID+'_HemoDataInflation.xls'
estimDataFileName = 'Simulation'+subject+'_'+pID+'_EstimdataInflation.xls'
gaussPtTableFileBaseName = dirName+'GPT_'+pID
inflationCont6FileName = dirName+baseFileName+'Inflation.cont6'
deflationCont6FileName = dirName+baseFileName+'Deflation.cont6'
nodalEDFileName = subject+'_'+pID+'_Nodes_ED_DT.xls'

pressureLVED = 815 # End Diastolic LV Pressure in Pa
pressureRVED = 0 # End Diastolic RV Pressures in Pa

basisFunctionIndex = 2
outputSteps = 10

regularDeflationSteps = 50
firstDeflationSteps = 100


lastInflationNodalSolution = int(pressureLVED)
maxIterations = 10

for iteration in range(0, maxIterations):
    # Create Output Directory
    currentSimulationsPerformed = 0;
    currentOutputDir = inflationOutputDir + '%d'%(iteration)
    try:
        os.mkdir(currentOutputDir)
    except OSError:
        sys.stderr.write("Warning %s already exists!\n\n\n Aborting \n\n\n"%currentOutputDir)
        currentSimulationsPerformed = 1;

    # Inflate the Mesh to EDP
    if (currentSimulationsPerformed == 0):

        self.Load_File(inflationCont6FileName, log=0)
        
        # Get Nodes
        nodes = self.stored_data.nodes.obj
        
        # Load Previous Deflation Step Unloaded Nodes
        if (iteration == 0):
            nodalDataUnloadedFileName = dirName + nodalEDFileName
        elif (iteration == 1):
            nodalDataUnloadedFileName = deflationOutputDir + '%d/NodalSolution%d.xls'%(0,firstDeflationSteps)
        else:
            nodalDataUnloadedFileName = deflationOutputDir + '%d/NodalSolution%d.xls'%(iteration-1,regularDeflationSteps)
        
        print "\n \n Loading Unloaded Data from " + nodalDataUnloadedFileName + "\n \n"
        
        nodes.loadTableFromFile(nodalDataUnloadedFileName,True)
        self.stored_data.store(nodes, modified = True)
        
        # Send and Calculate mesh
        self.Send(None, log=0)
        self.CalcMesh([('Calculate', None), ('Do not Calculate', None), ('Calculate', None), ('Global arc length scale factors (for nodal derivs wrt arc lengths)', None)], log=0)

        # Change Material Parameters
        f = open(dirName+'inputs'+pID+'.txt')
        self.stored_data.matEquations.obj['variables'][60].initial_cond[0][1]=float(f.readline()) # b (19.452)
        self.stored_data.matEquations.obj['variables'][61].initial_cond[0][1]=float(f.readline()) # a (0.76)
        self.stored_data.matEquations.obj['variables'][62].initial_cond[0][1]=float(f.readline()) # bf (31.558)
        self.stored_data.matEquations.obj['variables'][63].initial_cond[0][1]=float(f.readline()) # af (0.56)
        self.stored_data.store(self.stored_data.matEquations.obj)
        f.close()
    
        # Copy Nodes to Inital Conditions Form
        self.CopyNodesToIC(log=0)
        
        pressureLV = 0
        pressureRV = 0
        timeStepIncrement = 0.0
        initialTime = 1.0
        initialPressureLV = pressureLVED/25.0
        initialPressureStepLV = pressureLVED/4000.0
        initialPressureStepRV = initialPressureStepLV*(pressureRVED*1.0/pressureLVED)
        normalPressureStepLV = pressureLVED/100.0
        normalPressureStepRV = normalPressureStepLV*(pressureRVED*1.0/pressureLVED)    
   
        numCurrentSteps = outputSteps
        currentPressureStepLV = initialPressureStepLV
        currentPressureStepRV = initialPressureStepRV
        
        self.stored_data.circModel.obj['leftValues'] = [currentPressureStepLV/1000.0]
        if (biVModel == 1):
            self.stored_data.circModel.obj['rightValues'] = [currentPressureStepRV/1000.0]
            self.stored_data.circModel.obj['ICValues'] = [pressureLV/1000.0, pressureRV/1000.0]
        else:
            self.stored_data.circModel.obj['ICValues'] = [pressureLV/1000.0]
        self.stored_data.circModel.obj['convergedHemoFile'] = hemoDataFileName
        self.stored_data.circModel.obj['estimHemoFile'] = estimDataFileName
        
        # Send
        self.Send(None, log=0)
        self.CalcMesh([('Calculate', None), ('Do not Calculate', None), ('Calculate', None), ('Global arc length scale factors (for nodal derivs wrt arc lengths)', None)], log=0)
    
        self.Save(inflationCont6FileName, log=0)
        self.Reset(['Client', 'Server'], log=0)
      
        self.Load_File(inflationCont6FileName, log=0)
       
        # Send
        self.Send(None, log=0)
        self.CalcMesh([('Calculate', None), ('Do not Calculate', None), ('Calculate', None), ('Global arc length scale factors (for nodal derivs wrt arc lengths)', None)], log=0)
        
        
        while (int(pressureLV) < int(initialPressureLV)):
            pressureLV += currentPressureStepLV*numCurrentSteps
            pressureRV += currentPressureStepRV*numCurrentSteps
            nodalFileNum = int(pressureLV)
            NonLinSolve(self, initialTime, timeStepIncrement, numCurrentSteps)
            self.LnodalSolution(log=0, scripting = True, writeFile = "%s%d/NodalSolution%d.xls"%(inflationOutputDir, iteration, nodalFileNum))
            print "\n\n Current LV Pressure : %f"%(pressureLV)
        
        currentPressureStepLV = normalPressureStepLV
        currentPressureStepRV = normalPressureStepRV
        self.stored_data.circModel.obj['leftValues'] = [currentPressureStepLV/1000.0]
        if (biVModel == 1):
            self.stored_data.circModel.obj['rightValues'] = [currentPressureStepRV/1000.0]
            self.stored_data.circModel.obj['ICValues'] = [pressureLV/1000.0, pressureRV/1000.0]
        else:
            self.stored_data.circModel.obj['ICValues'] = [pressureLV/1000.0]
        
        # Load nodal data at end of initial inflation
        nodalDataInflatedFileName = inflationOutputDir + '%d/NodalSolution%d.xls'%(iteration, nodalFileNum)
        self.UICWithFile(nodalDataInflatedFileName, log=0)
        
        self.Save(inflationCont6FileName, log=0)
        self.Reset(['Client', 'Server'], log=0)
       
        self.Load_File(inflationCont6FileName, log=0)
        
        # Send
        self.Send(None, log=0)
        self.CalcMesh([('Calculate', None), ('Do not Calculate', None), ('Calculate', None), ('Global arc length scale factors (for nodal derivs wrt arc lengths)', None)], log=0)
        
        initialOutputSteps = outputSteps -  pressureLV/normalPressureStepLV
        while (initialOutputSteps <= 0):
            initialOutputSteps += outputSteps    
        numCurrentSteps = initialOutputSteps
        print "\n Performing initialOutputSteps %d\n\n"%initialOutputSteps
        nodalFileNum = int(pressureLV)
        while (int(pressureLV) < int(pressureLVED)):
            pressureLV += currentPressureStepLV*numCurrentSteps
            pressureRV += currentPressureStepRV*numCurrentSteps
            nodalFileNum = int(pressureLV)
            NonLinSolve(self, initialTime, timeStepIncrement, numCurrentSteps)
            self.LnodalSolution(log=0, scripting = True, writeFile = "%s%d/NodalSolution%d.xls"%(inflationOutputDir, iteration, nodalFileNum))
            print "\n\n Current LV Pressure : %f"%(pressureLV)
            numCurrentSteps = outputSteps
        
        self.LstressAndStrain({'outputPath':"%s%d/StressStrain%d.xls"%(inflationOutputDir, iteration, nodalFileNum),'elemlist':None,'inputPath':None,'varsWanted':['stress_out', 'F_out'],'xilist':None}, log=0)
	
        with open(continuityDir+'.continuity/working/Simulation'+subject+'_'+pID+'_HemoDataInflation.xls','r') as f:
            reader = csv.reader(f,delimiter = '\t')
            for row in reader:
                Vol = float(row[5])
        
        if (abs(Vol - Vref)/Vref > 0.03):
            f = open(dirName+'Verr'+pID+'.txt','w')
            f.write('0')
            f.close()
        else:
            f = open(dirName+'Verr'+pID+'.txt','w')
            f.write('1')
            f.close()
        
	break
    
    self.Reset(['Client', 'Server'], log=0)
    
    # Create Deflation Output Directory
    currentOutputDir = deflationOutputDir + '%d'%(iteration)
    currentSimulationsPerformed = 0
    try:
        os.mkdir(currentOutputDir)
    except OSError:
        sys.stderr.write("Warning %s already exists!\n"%currentOutputDir)
        currentSimulationsPerformed = 1
    
    if (currentSimulationsPerformed == 0):
        # Load Deflation file and copy inflated nodes
        self.Load_File(deflationCont6FileName, log=0)
        
        # Get Nodal Labels
        nodes = self.stored_data.nodes.obj
        
        nodalDataEDFileName = dirName+nodalEDFileName
          
        nodes = self.stored_data.nodes.obj
        nodes.loadTableFromFile(nodalDataEDFileName,True)
        self.stored_data.store(nodes, modified = True)
        
        nodalDataInflatedFileName = inflationOutputDir + '%d/NodalSolution%d.xls'%(iteration, lastInflationNodalSolution)
        self.UICWithFile(nodalDataInflatedFileName, log=0)
          
        # Send and Calculate mesh
        self.Send(None, log=0)
        self.CalcMesh([('Calculate', None), ('Do not Calculate', None), ('Calculate', None), ('Global arc length scale factors (for nodal derivs wrt arc lengths)', None)], log=0)
        
        print "\n \n Loading Unloaded Data from " + nodalDataEDFileName + "\n \n"
        print " Loading Inflated Data from " + nodalDataInflatedFileName + "\n \n"
        
        if (usingGPT):
            gaussPtTableFileName = gaussPtTableFileBaseName + str(iteration) + '.xls'
            # Calculate deformation gradient and save it to gauss point table
            self.LstressAndStrain({'outputPath':gaussPtTableFileName,'elemlist':None,'inputPath':None,'varsWanted':['Ftot_out'],'xilist':None}, log=0)
        else:
            # Calculate deformation gradient and save to Fields
            self.Rsurface({'VariableType':'tensor','XiLocation':'0.0','Render':False,'Deformed':0,'VariableField':'FieldVec1_Var1_','WhichXi':'3','XiList':None,'ColorMinMaxValues':[0.0, 0.0],'RunFromScript':True,'IndexOfField':1,'Variable':'Ftot_out','ElemList':['all'],'IndexOfBasis':basisFunctionIndex,'ColorRange':'Use Calculated Color Ranges','VariableComponent':'0 0'}, log=0)
            self.Rsurface({'VariableType':'tensor','XiLocation':'0.0','Render':False,'Deformed':0,'VariableField':'FieldVec1_Var2_','WhichXi':'3','XiList':None,'ColorMinMaxValues':[0.0, 0.0],'RunFromScript':True,'IndexOfField':2,'Variable':'Ftot_out','ElemList':['all'],'IndexOfBasis':basisFunctionIndex,'ColorRange':'Use Calculated Color Ranges','VariableComponent':'0 1'}, log=0)
            self.Rsurface({'VariableType':'tensor','XiLocation':'0.0','Render':False,'Deformed':0,'VariableField':'FieldVec1_Var3_','WhichXi':'3','XiList':None,'ColorMinMaxValues':[0.0, 0.0],'RunFromScript':True,'IndexOfField':3,'Variable':'Ftot_out','ElemList':['all'],'IndexOfBasis':basisFunctionIndex,'ColorRange':'Use Calculated Color Ranges','VariableComponent':'0 2'}, log=0)
            self.Rsurface({'VariableType':'tensor','XiLocation':'0.0','Render':False,'Deformed':0,'VariableField':'FieldVec2_Var4_','WhichXi':'3','XiList':None,'ColorMinMaxValues':[0.0, 0.0],'RunFromScript':True,'IndexOfField':4,'Variable':'Ftot_out','ElemList':['all'],'IndexOfBasis':basisFunctionIndex,'ColorRange':'Use Calculated Color Ranges','VariableComponent':'1 0'}, log=0)
            self.Rsurface({'VariableType':'tensor','XiLocation':'0.0','Render':False,'Deformed':0,'VariableField':'FieldVec2_Var5_','WhichXi':'3','XiList':None,'ColorMinMaxValues':[0.0, 0.0],'RunFromScript':True,'IndexOfField':5,'Variable':'Ftot_out','ElemList':['all'],'IndexOfBasis':basisFunctionIndex,'ColorRange':'Use Calculated Color Ranges','VariableComponent':'1 1'}, log=0)
            self.Rsurface({'VariableType':'tensor','XiLocation':'0.0','Render':False,'Deformed':0,'VariableField':'FieldVec2_Var6_','WhichXi':'3','XiList':None,'ColorMinMaxValues':[0.0, 0.0],'RunFromScript':True,'IndexOfField':6,'Variable':'Ftot_out','ElemList':['all'],'IndexOfBasis':basisFunctionIndex,'ColorRange':'Use Calculated Color Ranges','VariableComponent':'1 2'}, log=0)
            self.Rsurface({'VariableType':'tensor','XiLocation':'0.0','Render':False,'Deformed':0,'VariableField':'FieldVec3_Var7_','WhichXi':'3','XiList':None,'ColorMinMaxValues':[0.0, 0.0],'RunFromScript':True,'IndexOfField':7,'Variable':'Ftot_out','ElemList':['all'],'IndexOfBasis':basisFunctionIndex,'ColorRange':'Use Calculated Color Ranges','VariableComponent':'2 0'}, log=0)
            self.Rsurface({'VariableType':'tensor','XiLocation':'0.0','Render':False,'Deformed':0,'VariableField':'FieldVec3_Var8_','WhichXi':'3','XiList':None,'ColorMinMaxValues':[0.0, 0.0],'RunFromScript':True,'IndexOfField':8,'Variable':'Ftot_out','ElemList':['all'],'IndexOfBasis':basisFunctionIndex,'ColorRange':'Use Calculated Color Ranges','VariableComponent':'2 1'}, log=0)
            self.Rsurface({'VariableType':'tensor','XiLocation':'0.0','Render':False,'Deformed':0,'VariableField':'FieldVec3_Var9_','WhichXi':'3','XiList':None,'ColorMinMaxValues':[0.0, 0.0],'RunFromScript':True,'IndexOfField':9,'Variable':'Ftot_out','ElemList':['all'],'IndexOfBasis':basisFunctionIndex,'ColorRange':'Use Calculated Color Ranges','VariableComponent':'2 2'}, log=0)
      

        self.Save(deflationCont6FileName, log=0)
    
        # Load  Unloaded Mesh to start correction
        if (iteration != 0):
            if (iteration == 1):
                nodalDataUnloadedFileName = deflationOutputDir + '%d/NodalSolution%d.xls'%(0, firstDeflationSteps)
            else:
                nodalDataUnloadedFileName = deflationOutputDir + '%d/NodalSolution%d.xls'%(iteration-1, regularDeflationSteps)
            nodes = self.stored_data.nodes.obj
            nodes.loadTableFromFile(nodalDataUnloadedFileName, True)
            self.stored_data.store(nodes, modified = True)
        
        # Copy Nodes to Inital Conditions Form
        self.CopyNodesToIC(log=0)
        
        # Send and Calculate mesh
        self.Send(None, log=0)
        #self.CalcMesh([('Calculate', None), ('Do not Calculate', None), ('Calculate', None), ('Global arc length scale factors (for nodal derivs wrt arc lengths)', None)], log=0)
        
        if (not usingGPT):
            nodalDataUnloadedWithFFileName = dirName + 'UnloadedNodalDataWithF_' + pID + '_%d.xls'%(iteration)
            nodes = self.stored_data.nodes.obj
            nodes.saveTable(nodalDataUnloadedWithFFileName)
            
            self.Save(deflationCont6FileName, log=0)
    
            # For initial mesh load the previous unloaded state as undeformed geometry
            nodalDataUnloadedWithFFileName = dirName + 'UnloadedNodalDataWithF_' + pID + '_%d.xls'%(iteration)
            nodes = self.stored_data.nodes.obj
            nodes.loadTableFromFile(nodalDataUnloadedWithFFileName)
            self.stored_data.store(nodes, modified = True)
            self.Save(deflationCont6FileName, log=0)
                    
            # Copy Nodes to Inital Conditions Form
            self.CopyNodesToIC(log=0)
        
            # Send and Calculate mesh
            self.Send(None, log=0)
            #self.server.model.mesh.CalcMeshAlready = True
            #self.CalcMesh([('Calculate', None), ('Do not Calculate', None), ('Calculate', None), ('Global arc length scale factors (for nodal derivs wrt arc lengths)', None)], log=0)
        else:
            gaussPtTableFileName = gaussPtTableFileBaseName + str(iteration) + '.xls'
            self.stored_data.matEquations.obj['variables'].setIcByName('Fg11',[['Gauss Point Table:Ftot_out[0,0]', gaussPtTableFileName]])
            self.stored_data.matEquations.obj['variables'].setIcByName('Fg12',[['Gauss Point Table:Ftot_out[0,1]', gaussPtTableFileName]])
            self.stored_data.matEquations.obj['variables'].setIcByName('Fg13',[['Gauss Point Table:Ftot_out[0,2]', gaussPtTableFileName]])
            self.stored_data.matEquations.obj['variables'].setIcByName('Fg21',[['Gauss Point Table:Ftot_out[1,0]', gaussPtTableFileName]])
            self.stored_data.matEquations.obj['variables'].setIcByName('Fg22',[['Gauss Point Table:Ftot_out[1,1]', gaussPtTableFileName]])
            self.stored_data.matEquations.obj['variables'].setIcByName('Fg23',[['Gauss Point Table:Ftot_out[1,2]', gaussPtTableFileName]])
            self.stored_data.matEquations.obj['variables'].setIcByName('Fg31',[['Gauss Point Table:Ftot_out[2,0]', gaussPtTableFileName]])
            self.stored_data.matEquations.obj['variables'].setIcByName('Fg32',[['Gauss Point Table:Ftot_out[2,1]', gaussPtTableFileName]])
            self.stored_data.matEquations.obj['variables'].setIcByName('Fg33',[['Gauss Point Table:Ftot_out[2,2]', gaussPtTableFileName]])
            
        self.Save(deflationCont6FileName, log=0)
        self.Reset(['Client', 'Server'], log=0)
        
        self.Load_File(deflationCont6FileName, log=0)

        # Send and Calculate mesh
        self.Send(None, log=0)
        self.CalcMesh([('Calculate', None), ('Do not Calculate', None), ('Calculate', None), ('Global arc length scale factors (for nodal derivs wrt arc lengths)', None)], log=0)
        
        #Perform deflation
        initialTime  = 0.0
        numCurrentSteps = outputSteps
        if (iteration == 0):
            numDeflationSteps = firstDeflationSteps
        else:
            numDeflationSteps = regularDeflationSteps
        timeStepIncrement = 1.0/(numDeflationSteps)
        for deflationStep in range(0 ,numDeflationSteps, outputSteps):
            finalTime = initialTime + outputSteps*timeStepIncrement;
            NonLinSolve(self, initialTime, timeStepIncrement, numCurrentSteps)
            initialTime = finalTime
            self.LnodalSolution(log=0, scripting = True, writeFile = "%s%d/NodalSolution%d.xls"%(deflationOutputDir, iteration, deflationStep+outputSteps))
              
        self.Save(deflationCont6FileName, log=0)
	
	break
    
    self.Reset(['Client', 'Server'], log=0)

