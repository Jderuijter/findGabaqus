# -*- coding: utf-8 -*-
"""
part 2
function to calculate differences between reference and simulation
part 3
Small loop with Downhill simplex optimization to determine new material parameters
@author: Joerik de Ruijter
"""

from abaqus import *
from abaqusConstants import *
import __main__
import section
import regionToolset
import displayGroupMdbToolset as dgm
import part
import material
import assembly
import step
import interaction
import load
import mesh
import optimization
import job
import sketch
import visualization
import xyPlot
import displayGroupOdbToolset as dgo
import connectorBehavior
#
from numpy import *
from math import *
import heapq
#PART 1 ######################################################################################################################

#PART 2 ######################################################################################################################

def CompareDisplacements(startINP,jobNameRef1,MaterialModuli,niteration):
    def readfile ( filename ):
    	#opens a file and returns the lines
    	inputFile  = open(filename,'r')
    	inp_lines = inputFile.read().splitlines()
    	inputFile.close()
    	return inp_lines
    	
    
    def createjobfile ( lines, newjobName ):
    	#create a .inp file from lines
    	newInputFile = open(newjobName+'.inp','w')
    	for i in range(len(lines)):
    		newInputFile.write(lines[i]+'\n')
    	newInputFile.close()
    	
    def runabaqus ( jobName, nstep ): 
    	#run abaqus jobName and nstep = name of step in abaqus (Step-1)	
    	#creates jobName_nodal.rpt with displacement info
    	#creates jobName.rpt with stress info
    	inp = str(jobName) + '.inp'
    	
    
    	mdb.JobFromInputFile(str(jobName), str(inp))       
    	mdb.jobs[str(jobName)].submit()
    	mdb.jobs[str(jobName)].waitForCompletion()   
    
    
    	mySession = session.openOdb(name=str(jobName)+'.odb')
    	session.viewports['Viewport: 1'].setValues(displayedObject=mySession)
    	session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
    	UNDEFORMED, ))        
    	myOdb = session.odbs[str(jobName)+'.odb']
    	myStep = myOdb.steps[nstep]
    	numFrames = len(myStep.frames)
    	session.fieldReportOptions.setValues(printTotal=OFF, )
    	lastFrameNo = numFrames - 1
    	session.writeFieldReport(fileName=str(jobName)+'.rpt', append=OFF, 
    		sortItem='Element Label', odb=myOdb, step=0, frame=lastFrameNo, 
    		outputPosition=ELEMENT_CENTROID, variable=(('S', INTEGRATION_POINT, ((
    		COMPONENT, 'S11'), (COMPONENT, 'S22'), (COMPONENT, 'S33'), (COMPONENT, 
    		'S12'),  (COMPONENT, 'S13'), (COMPONENT, 'S23'), )), ))
    	session.writeFieldReport(fileName=str(jobName)+'_nodal'+'.rpt', append=OFF, 
    	sortItem='Node Label', odb=myOdb, step=0, frame=lastFrameNo, outputPosition=NODAL, 
    		variable=(('U', NODAL), ))
    	myOdb.close()
    
    def getREFinfo (lines):
    	#lines = lines from Reference file.
    	#nodeRef = name for list nodenrs
    	#URef = name for list displacements
    	nodeRef=[]
    	URef=[]
    	for counter in range(len(lines)):  #reference nodenr and U from input file.
    		[nodeRefUSN, nodeRefFEMN, URefN]=lines[counter].split()
    		nodeRef.append(float(nodeRefFEMN))
    		URef.append(float(URefN))
    	return (nodeRef,URef)
    
    def findlines(lines,line_begin,line_end):
	#loop over lines
	#find locations of the initial result --> ini, and end result--> end match with line_begin and line_end
	ini=[]
	end=[]
	for counter in range(len(lines)):
		if lines [counter] == line_begin:
			ini.append(counter)
			for counter2 in range(len(lines)-counter):
				if lines [counter+counter2] == line_end:
					end.append(counter+counter2-1)
					break	
	return (ini,end)
    def getdisplacement (ini, end, lines,nodeRef):
	#find displacement U(1,2,3) in lines and match with nodenrs in nodeRef
	nodeNr=[]
	Ulist=[]
	UlistMatched=[]
	
	for counter in range(len(lines)):		#get nodenr and displacement? (U1) from job-1
		if counter > ini and counter <= end:
			[nNr,Um,U1,U2,U3]=lines[counter].split()
			nodeNr.append(float(nNr))
			Ulist.append(float(U3))

	for counter in range(len(nodeRef)):
		UlistMatched.append(abs(Ulist[int(nodeRef[counter]-1)]))
	return UlistMatched
    
    def writeinpfile (text,newjobName,startinp,textline):
    	linesstart=readfile(startinp)
    	newInputFile = open(str(newjobName)+'.inp','w')
    	for counter in range(len(linesstart)):
    		if linesstart[counter-2] == textline:
    			
    			newInputFile.write(text)
    		else:
    			newInputFile.write(linesstart[counter]+'\n')
    	newInputFile.close()
     
    ###############################################################################################
     #Call definitions
    line_1= '---------------------------------------------------------------------------------'
    line_2 = ""

 
    line_12='*MATERIAL, name=outer2'
    line_13='*MATERIAL, name=outer1'
    line_14='*MATERIAL, name=middle1'
    line_15='*MATERIAL, name=inner2'
    line_16='*MATERIAL, name=inner1'
    if niteration==1:      #initialisation only at the first call
        inp_lines = readfile(startINP)
        createjobfile(inp_lines,'Simulation_1')
        overviewFile = open('overview.rpt','w')
        overviewFile.write('inputfile=' + startINP + '\n\n')	
        overviewFile.close()

   
    newname='Simulation_'+str(niteration)
    
    text=str(MaterialModuli[1]) + ', ' + str(1e-05) + '\n'
    text1=str(MaterialModuli[0]) + ', ' + str(1e-05) + '\n'
    
    writeinpfile(text,newname,startINP,line_16)
    writeinpfile(text,newname,newname+'.inp',line_15)
    writeinpfile(text,newname,newname+'.inp',line_14)
    writeinpfile(text1,newname,newname+'.inp',line_13)
    writeinpfile(text1,newname,newname+'.inp',line_12)    

    
    runabaqus(newname,'Step-1')
    ref_lines1=readfile(jobNameRef1)
   
    (nodeRef1,URef1)=getREFinfo(ref_lines1)
  
    rpt_lines=readfile(newname+'_nodal'+'.rpt')
    (start,end)=findlines(rpt_lines,line_1,line_2)
    Ulistmatch=getdisplacement(start[0],end[0],rpt_lines,nodeRef1)
   
    #calculate differences between Reference displacements and FEM displacements
    UDifabs = [abs((y-x)) for x, y in zip(Ulistmatch, URef1)]
    Averror = sum(UDifabs)/len(UDifabs)
 

    #write overview
    overviewFile = open('overview.rpt','a')
    newOvLine = newname + '\t\t' + 'Shear modulus inner G = \t\t' + str(2*MaterialModuli[1]) + ' MPa'+ '\t\t'+ 'Shear modulus outer G = \t\t' + str(2*MaterialModuli[0]) + ' MPa' + '\t\t' + 'Udif'+ '\t\t' + str(Averror) + '\n'
    overviewFile.write(newOvLine)	
    overviewFile.close()
    print newname + '\t\t'+ str(Averror)
	
    #write new inpfiles
 
    niteration=niteration+1
    return (Averror, niteration)	

def DownhillSimplex(Start, slide, tol,ref,inp):
    j=1
    # Setup intial values
    newXmin=[5,2.5]	#min boundaries
    n = len(Start)
    f = zeros(n+1)
    x = zeros((n+1,n))
    d = zeros(n)	
    x[0] = Start
    	
    # Setup intial X range
   for i in range(1,n+1):
        x[i] = Start
        x[i,i-1] = Start[i-1] + slide
	
	# Setup intial functions based on x's just defined
	
    for i in range(n+1):
        (f[i],j) = CompareDisplacements(inp,ref,x[i],j)
	
	# Main Loop operation, loops infinitly until break condition
    counter = 0
    while True:
        low = argmin(f)
        high = argmax(f)
        ind=heapq.nlargest(2, range(len(f)), f.__getitem__)
        second=ind[1]
    	 
        # Implement Counter 
        counter+=1
    	  	
        # Compute Migration
        for i in range(n):
            d[i] = (-(n+1)*x[high,i]+sum(x[:,i]))/n
        print d		
        if sqrt(dot(d,d)/n)<tol:
        # Break condition, value is darn close
            return (x[low],counter)
    	  #try reflection		
        newX = x[high] + 2.0*d 
        print newX
        print 'reflectionpoint'
        for i in range(n):
            if newX[i] < newXmin[i]:
               newX[i] = newXmin[i]+random.uniform(0,newXmin[i]*0.1)

        (newF,j) = CompareDisplacements(inp,ref,newX,j)
    		
        # Bad news, new values is lower than p. low
    		
        if newF < f[low]:  #accept reflection
            x[high] = newX
            f[high] = newF
            #try expanding the reflection
            newX = x[high] + d #expansion
            for i in range(n):
                if newX[i] < newXmin[i]:
                   newX[i] = newXmin[i]+random.uniform(0,newXmin[i]*0.1)
            (newF,j) = CompareDisplacements(inp,ref,newX,j)
			# Maybe I need to expand
            if newF < f[high]: #accept expansion
                x[high] = newX
                f[high] = newF
		
        else:
            # try reflection again 
            if newF <= f[second]: #accept reflection
                x[high] = newX
                f[high] = newF
            else:
                if newF <= f[high]:
                    x[high] = newX
                    f[high] = newF
                # try Contraction
                newX = x[high] + 0.5*d #contraction
                for i in range(n):
                    if newX[i] < newXmin[i]:
                       newX[i] = newXmin[i]+random.uniform(0,newXmin[i]*0.1)
                (newF,j) = CompareDisplacements(inp,ref,newX,j)
                if newF <= f[high]: #accept contraction
                    x[high] = newX
                    f[high] = newF
                else:        #use shrinkage
                    for i in range(len(x)):
                        if i!=low:
                            
                            x[i] = (x[i]+x[low])*0.5
                            for i in range(n):
                                if newX[i] < newXmin[i]:
                                   newX[i] = newXmin[i]+random.uniform(0,newXmin[i]*0.1)                              
                            (f[i],j) = CompareDisplacements(inp,ref,x[i],j)
# Example Call
(result,counter)=DownhillSimplex(array([25,12]),1,1e-6,'D40_1703_multi_lower_ref.inp','D40_1703_multi_lower.inp')
      
        
print "Custom Downhill Simplex:"
print "Final Result: " + str(result)
print "Iteration Count: " + str(counter)
        
        
