#!/usr/bin/env python3
# coding=utf-8

path = "/home/laverre/Regulatory_landscape/data/pLi_scores/"

GeneFile = path + "GeneID_TranscriptID_Ensembl94.txt"
PLiFile = path + "Transcript_pLI_march16_2016.txt"
OutputFile = path + "Gene_pLI.txt"
TranscriptError = path + "missing_Transcripts.txt"

Transcript2Gene = {}
with open(GeneFile, 'r') as f:
    for i in f.readlines()[1:]:
        i = i.strip("\n")
        i = i.split("\t")
        Transcript2Gene[str(i[1])] = str(i[0])

Transcript2PLi = {}
with open(PLiFile, 'r') as f:
    for i in f.readlines()[1:]:
        i = i.strip("\n")
        i = i.split("\t")
        Transcript2PLi[str(i[0])] = str(i[19])


output = open(OutputFile, 'w')
error = open(TranscriptError, 'w')

output.write("GeneID\tTranscriptID\tpLI\n")
for Transcript in Transcript2PLi.keys():
    simpleTranscript = Transcript.split(".")[0]

    if simpleTranscript in Transcript2Gene.keys():
        GeneID = Transcript2Gene[simpleTranscript]
        output.write(GeneID + "\t" + Transcript + "\t" + Transcript2PLi[Transcript] + "\n")

    else:
        error.write(Transcript + "\n")


output.close()
