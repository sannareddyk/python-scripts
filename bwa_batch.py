import sys
import subprocess
import time

###Input
if len(sys.argv)!=2:
    print("usage:python bwa_batch.py listOfFilesFileName.txt")
    sys.exit()
###Global parameters
refSeqFileName='path/to/reference/Dyak.fasta'
bwaPath='/mnt/lustre/home/sannareddyk/bwa-0.6.2/bwa'
###Functions
def run_cmd(cmd):
	'''
        
	'''
    p =subprocess.Popen(cmd,shell=True,stderr=subprocess.PIPE)
    err = p.communicate()[1]
	return (p.returncode, err)

###Reading files and creating a list
listOfFilesFileName=sys.argv[1]
listOfFilesFile=open(listOfFilesFileName)
listOfFilesList=[]

for line in listOfFilesFile:
    sample_name=line.rstrip()
    listOfFilesList.append(sample_name)
print listOfFilesList
listOfFilesFile.close()

###For each file in the list, run bwa aln and samse steps
for fileName in listOfFilesList:
    print('Processing '+fileName)
	
	start=time.time()
    ##Run bwa aln
    outFileName=fileName.split('.')[0]+'.sai'
	erFileNameAln=fileName.split('.')[0]+'_aln_er.txt'
	numThreads=4
	myCmd=""+bwaPath+" aln -t "+str(numThreads)+" "+refSeqFileName+" "+fileName+" >"+outFileName
	#print(myCmd)
	rc,er=run_cmd(myCmd)
	if rc!=0:
        outFile1=open(erFileNameAln,'w')
        outFile1.write(str(er))
		outFile1.close()
	elapsed=round((time.time()-start)/60.0)
	run1=str(elapsed)
	print "run1 :",run1
	#Sleep
    time.sleep(5)
	
	##Run samse
    finalOutFileName=fileName.split('.')[0]+'.sam'
	#print(finalOutFileName)
	erFileNameSamse=fileName.split('.')[0]+'_samse_er.txt'
	myCmd=""+bwaPath+" samse "+refSeqFileName+" "+outFileName+" "+fileName+" >"+finalOutFileName
	#print(myCmd)
	rc,er=run_cmd(myCmd)
	if rc!=0:
		outFile2=open(erFileNameSamse,'w')
        outFile2.write(str(er))
        outFile2.close()

###Debug
#myCmd="ls -ltr >trial.txt"
#rc,er=run_cmd(myCmd)
#print(rc)

