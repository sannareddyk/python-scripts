import sys

if len(sys.argv) != 4:
	print("usage: python mergeFiles.py a.txt b.txt ab.txt")
	sys.exit()


####################
## Get the first file into a dictionary

aFileName = sys.argv[1]
aFile = open(aFileName)

aDict = {}

for line in aFile:
	line = line.strip().split("\t")
	chrom = line[0]
	loc = line[1]
	key = (chrom,loc)
	line=line[2:]
	line="\t".join(line)
	#print line
	aDict[key] = line
aFile.close()

##################
bFileName = sys.argv[2]
bFile = open(bFileName)

outFileName = sys.argv[3]
outFile = open(outFileName,'w')

for line in bFile:
	line = line.strip().split('\t')
    #print(line)
	chrom = line[0]
	loc = line[1]
	key = (chrom,loc)
	line = line[2:]

	if key in aDict:
		line_aFile = aDict[key]
		line = "\t".join(line)
		key = "\t".join(key)
		resultLine = key + '\t' + line_aFile + "\t" + line  	
		outFile.write(resultLine + '\n')
        
aFile.close()
bFile.close()
outFile.close()


