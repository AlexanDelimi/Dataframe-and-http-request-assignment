from fileinput import filename
import os, sys, getopt
from pickle import FALSE, TRUE
import re
import csv
import math
import numpy
import pandas as pd

curDir="C:\\Users\\AlexD\\Desktop\\MSc_Bioinformatics\\21-22_MET586_Proteomics\\Pyhton\\python\\datafiles" #Directory where only the datafiles are stored
a2 = os.listdir(curDir) #Listing the contents of the datafile directory
a2df = pd.DataFrame({'filename':a2}) #Creating a pandas Dataframe with the datafiles

dataframe_collection = {} #dataframe that will contain all the individual files' dataframes
firstdataframe={} #dataframe that contains only the columns that will be in the produced file

for filenam in a2df['filename']: # For each file
    fname=filenam.replace('full.',"") #Create a string with the filename but replace the full. in their name with null
    fname=fname.replace('.markdup.count',"") #and replace the .mardkup.count with null
    dataframe_collection[fname] = pd.read_table(os.path.join(curDir,filenam), skiprows=2,index_col=0, usecols=[0,5,6], names= ["Geneid", "Length", fname] )
    #read the file into a dataframe, use fisrt column (Geneid) as index, keep, 1st, 6th and 7th column and use names Geneid, Length and the previously edited string as column names. 
    dataframe_collection[fname]=dataframe_collection[fname].where(dataframe_collection[fname]['Length'] > 1000) #If Length is less than 1000 make Length nan
    dataframe_collection[fname]=dataframe_collection[fname][~numpy.isnan(dataframe_collection[fname]['Length'])] #drop nan length rows
    firstdataframe[fname]=dataframe_collection[fname].drop('Length', axis=1) #drop the length and add them to the next dataframe


for i in range (0,len(firstdataframe)): #for each dataframe in the dataframe firstdataframe
    if i == 0: #if its the first dataframe store it in df
        df = firstdataframe[list(firstdataframe)[i]]
    else:
        df1 = firstdataframe[list(firstdataframe)[i]] #else store it in df1
        df = df.merge(df1, on='Geneid') #merge df and df1 and store them in df

df.to_csv("MET586_Python_Assignment_1.tsv",header=True,index=True,sep='\t') #print the final dataframe to a tsv file. 


