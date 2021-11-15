import pandas as pd                                                                                                                                                                           

files = ("Disc1","Disc2")                                                                                                                                         

for i in files: 
    file_gid = i+"_GID.txt" 
    file_def = i+".txt" 
    data = pd.read_csv(file_gid,header=None,sep='\s+',engine='python') 
    (data[0]-1).to_csv(file_def,header=None,index=False) 
