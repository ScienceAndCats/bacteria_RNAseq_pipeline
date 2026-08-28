### iImport fnmatch module to support Unix shell-style wildcards(these are symbols that take the place of an unknown character or characters(like the *)) ###
import fnmatch

### Import os module to use operating system dependent functionality ###
import os
import re

BOWTIE_BACTERIA_DIR = os.path.join("bowtie_alignments", "bacteria")

### Creates a list of featureCounts files to count ###
featureCountsFiles = []
for file in os.listdir("."):
    if fnmatch.fnmatch(file, "featurecounts_BACTERIA_summary.txt"):
        featureCountsFiles.append(file)

### Creates list of cutadapt files to count ##
cutadaptFiles = []
for file in os.listdir("."):
    if fnmatch.fnmatch(file, "*cutadapt_log.txt"):
         cutadaptFiles.append(file)
cutadaptFiles = sorted(cutadaptFiles, key=str.casefold)

### Creates list of bowtie files to count ###
bowtieFiles = []
for file in os.listdir(BOWTIE_BACTERIA_DIR):
    if fnmatch.fnmatch(file, "BACTERIA*.bowtie_output.txt"):
        bowtieFiles.append(os.path.join(BOWTIE_BACTERIA_DIR, file))
bowtieFiles= sorted(bowtieFiles, key=str.casefold)

### Creates list of coverage files to count ###
coverageFilesList = []
for file in os.listdir(BOWTIE_BACTERIA_DIR):
    if fnmatch.fnmatch(file, "BACTERIA*_coverage.txt"):
        coverageFilesList.append(os.path.join(BOWTIE_BACTERIA_DIR, file))
coverageFilesList= sorted(coverageFilesList, key=str.casefold)

### Creates a list of filenames ###
filenames = []
for i in cutadaptFiles:
    f = i.strip()
    filenames.append(f.replace('_22bp.cutadapt_log.txt', ''))
filenames = sorted(filenames, key=str.casefold)

### Creates output csv template
import csv
full_path = os.getcwd().split("/")
projectName = full_path[len(full_path)-1]
summary_basename = "BACTERIA_" + projectName + ".csv"
csvfile=open( summary_basename, 'w')
csvfile.write("Sample, Raw Reads, Reads Written, Overall Allignment, Percent Covered, Mean Coverage, Reads Assigned to features, Features with non-zero reads\n")

### Opens cutadapt txt files and saves number of raw reads and reads written ###
rawReadsList = []
readsWrittenList = []
for file in cutadaptFiles:
    with open(file, "r") as f:
        report = f.read()
    raw_match = re.search(r"Total (?:read pairs|reads) processed:\s*([\d,]+)", report)
    written_match = re.search(r"(?:Pairs|Reads) written \(passing filters\):\s*([\d,]+)", report)
    if not raw_match or not written_match:
        raise ValueError(f"Could not parse cutadapt report: {file}")
    rawReadsList.append(raw_match.group(1).replace(',', ''))
    readsWrittenList.append(written_match.group(1).replace(',', ''))

### Opens bowtie txt files and saves overall allignment and reads not mapped to decoy ###
#notMappedtoDecoyList = []
overallAlignmentList = []
for file in bowtieFiles:
    with open(file, "r") as f:
        report = f.read()
    aligned = re.findall(
        r"^\s*(\d+) .*aligned (?:concordantly )?(?:exactly 1 time|>1 times)",
        report,
        re.MULTILINE,
    )
    overallAlignmentList.append(str(sum(map(int, aligned))))

### Opens coverageFiles and makes lists of percent coverage and mean coverage ###
percentCoverageList = []
meanCoverageList = []
for file in coverageFilesList:
    f = open(file, "r")
    lineone = f.readline()
    linetwo = f.readline()
    linethree = f.readline()
    linefour = f.readline()
    linefive = f.readline()
    strip = linefive.replace(' ','')
    split = strip.split(':')
    percentCoverage = split[len(split)-1].replace('\n', '')
### Opens cutadapt txt files and saves number of raw reads and reads written ###
rawReadsList = []
readsWrittenList = []
for file in cutadaptFiles:
    with open(file, "r") as f:
        report = f.read()
    raw_match = re.search(r"Total (?:read pairs|reads) processed:\s*([\d,]+)", report)
    written_match = re.search(r"(?:Pairs|Reads) written \(passing filters\):\s*([\d,]+)", report)
    if not raw_match or not written_match:
        raise ValueError(f"Could not parse cutadapt report: {file}")
    rawReadsList.append(raw_match.group(1).replace(',', ''))
    readsWrittenList.append(written_match.group(1).replace(',', ''))

### Opens bowtie txt files and saves overall allignment and reads not mapped to decoy ###
#notMappedtoDecoyList = []
overallAlignmentList = []
for file in bowtieFiles:
    with open(file, "r") as f:
        report = f.read()
    aligned = re.findall(
        r"^\s*(\d+) .*aligned (?:concordantly )?(?:exactly 1 time|>1 times)",
        report,
        re.MULTILINE,
    )
    overallAlignmentList.append(str(sum(map(int, aligned))))

### Opens coverageFiles and makes lists of percent coverage and mean coverage ###
percentCoverageList = []
meanCoverageList = []
for file in coverageFilesList:
    f = open(file, "r")
    lineone = f.readline()
    linetwo = f.readline()
    linethree = f.readline()
    linefour = f.readline()
    linefive = f.readline()
    strip = linefive.replace(' ','')
    split = strip.split(':')
    percentCoverage = split[len(split)-1].replace('\n', '')
    percentCoverageList.append(percentCoverage)
    linesix = f.readline()
    stripe = linesix.replace(' ','')
    splite = stripe.split(':')
    meanCoverage = splite[len(splite)-1].replace('\n', '')
    meanCoverageList.append(meanCoverage)

### Coverage info is put into txt files ###
percent= open("BACTERIA_percentCoverage.txt","w+")
for i in percentCoverageList:
    percent.write(i)
    percent.write('\n')
percent.close()

mean= open("BACTERIA_meanCoverage.txt","w+")
for i in meanCoverageList:
    mean.write(i)
    mean.write('\n')
mean.close()

### Puts data for mean coverage into a list ###
meanCoverageList = []
with open("BACTERIA_meanCoverage.txt", "r") as a_file:
    for line in a_file:
        meanCoverage = line.replace('\n', '')
        meanCoverageList.append(meanCoverage)

### Puts data for percent coverage into a list ###
percentCoverageList = []
with open("BACTERIA_percentCoverage.txt", "r") as a_file:
    for line in a_file:
        percentCoverage = line.replace('\n', '')
        percentCoverageList.append(percentCoverage)


### Opens featureCounts txt files and saves info to 'data' variable ###
alignmentCountList = []
geneCountList = []
for file in featureCountsFiles:
    f = open(file, "r")
    topComment = f.readline()
    header = f.readline()
    data = f.readlines()
    f.close()


    ### Reformats txt file data into list of strings ###
    headerStrip = header.strip("\n")
    headerList = headerStrip.split("\t")
    sampleFiles = headerList[6:]

    dataList = []
    for line in data:
        lineStrip = line.strip("\n")
        lineList = line.split("\t")
        dataList.append(lineList)

    ### Creates dictionaries with samples as keys ###
    alignmentCount = {}
    geneCount = {}

    for sample in sampleFiles:
        alignmentCount[sample] = 0
        geneCount[sample] = 0

    ### Counts reads and genes ###
    for data in dataList:

     if len(data) >= 7:
            for index,element in enumerate(data):
                if index>=6:
                    num_reads = int(element)
                    sample_match = sampleFiles[index - 6]

                    alignmentCount[sample_match] += num_reads
                    if num_reads != 0:
                        geneCount[sample_match] += 1

    for sample in sampleFiles:
        alignmentCountList.append(alignmentCount[sample])
        geneCountList.append(geneCount[sample])


### Adds data to csv file ###
i = 0
for file in cutadaptFiles:
    csvfile.write ("{0},{1},{2},{3},{4},{5},{6},{7}\n".format(filenames[i], rawReadsList[i], readsWrittenList[i], overallAlignmentList[i], percentCoverageList[i], meanCoverageList[i], alignmentCountList[i], geneCountList[i]))
    
    i = i + 1

csvfile.close()
