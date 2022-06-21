import os, sys, getopt
from pickle import FALSE, TRUE
import numpy
import re
import csv
import math
import pandas as pd
import requests


curDir="C:\\Users\\AlexD\\Desktop\\MSc_Bioinformatics\\21-22_MET586_Proteomics\\Pyhton\\python\\datafiles" #Directory where only the datafiles are stored
a2 = os.listdir(curDir) #Listing the contents of the datafile directory
a2df = pd.DataFrame({'filename':a2}) #Creating a pandas Dataframe with the datafiles

dataframe_collection = {} #dataframe that will contain all the individual files' dataframes
firstdataframe={} #dataframe that contains only the columns that will be in the produced file
ids_to_keep={} #dataframe with the Geneids 

for fname in a2df['filename']: #For each file in the path
    dataframe_collection[fname] = pd.read_table(os.path.join(curDir,fname), skiprows=2,index_col=0, usecols=[0,1,2,3,4,5,6], names=['Geneid','Chromosome','Start','End','Strand','Length','Counts'])
    #read the file into a dataframe, use fisrt column (Geneid) as index, keep, 1st through 7th column and use names Geneid, Length and the previously edited string as column names. 
    dataframe_collection[fname]=dataframe_collection[fname].where((dataframe_collection[fname]['Length'] > 1000) & (dataframe_collection[fname]['Counts']>0))# #If Length is less than 1000 or the third column is less than zero make Length nan
    dataframe_collection[fname]=dataframe_collection[fname][~numpy.isnan(dataframe_collection[fname]['Length'])]  #drop nan length rows
    firstdataframe[fname]=dataframe_collection[fname].drop(['Length', 'Counts'], axis=1) #drop the length and counts add them to the next dataframe
    ids_to_keep[fname]=dataframe_collection[fname].drop(['Chromosome','Start','End','Strand','Length','Counts'], axis=1) #Drop every columns besides geneId and store them in the ids_to_keep


for i in range (0,len(ids_to_keep)):
    if i == 0:
        df = ids_to_keep[list(ids_to_keep)[i]]
    else:
        df1 = ids_to_keep[list(ids_to_keep)[i]]
        df = df.merge(df1, on='Geneid')
#merge the dataframes into one, essentially keeping only the geneIds that exist in all of the individual files/dataframes

geneIdsList=df.reset_index()['Geneid'].tolist() #make the dataframe's index a list
df=df.merge(firstdataframe[list(firstdataframe)[0]], on = 'Geneid') # keep the first dataframe's data with ids that are in the df dataframe ( so geneIds that are in all of the files) nd store it in df


server = "https://rest.ensembl.org" #Server is ensembl's server

resultArr=[] #list with results
for g in geneIdsList:	#For each gene in the list of the common genes
    ext1 = "/overlap/id/" + g + "?feature=gene"	
    r1 = requests.get(server+ext1, headers={ "Content-Type" : "application/json"}) #make an hhtp request for that gene
    if r1.ok:  		#if the response is valid
        decoded1 = r1.json() #decode the json respone in a dataframe
        try:
            if decoded1[0]['biotype'] == 'protein_coding': #if the gene is protein coding
                try:
                    geneName = decoded1[0]['external_name'] # keep its external name
                except KeyError:
                    geneName='None' #if it doesn't have one then make it None
                try:
                    geneDescription = decoded1[0]['description'] #and its description 
                except KeyError:
                    geneDescription='None' #if it doesn't have one then make it None
                resultArr.append({'Geneid':g,'GeneName':geneName,'Chromosome':df.loc[[g]]['Chromosome'].item(),'Start':df.loc[[g]]['Start'].item(),'End':df.loc[[g]]['End'].item(),'Strand':df.loc[[g]]['Strand'].item(),'Description':geneDescription})
                #append in the result list a dictionary item containing the Geneid, the genename, the chromosome, the start, end positions, strand data and description of said gene 
        except KeyError:
            pass #if there is no biotype then just ignore the gene

df1 = pd.DataFrame(resultArr)#make the list into a datarame
df1.to_csv("MET586_Python_Assignment_3.tsv",header=True,index=False, sep='\t') #export the results in a tsv file







