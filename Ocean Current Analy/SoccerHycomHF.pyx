# -*- coding: utf-8 -*-
"""
Created on Mon Feb 27 11:20:43 2017

@author: qzh14

"""
from netCDF4 import Dataset
import datetime
#import time
from dateutil import parser

import patoolib
from dateutil.relativedelta import relativedelta
import os
import numpy as np
import math

import time
 





def lats_lons(dataset):
    #Input the dataset of the file that has been read.
    #Outputs the longitude and latitude of the file    
    try:
       lon = dataset.variables['lon'][:]
       lat = dataset.variables['lat'][:]
    except:
       lon = dataset.variables['longitude'][:]
       lat = dataset.variables['latitude'][:]
       lon = lon.T
       lat = lat.T
        
    for t in range(len(lon)):
         if lon[t]<0:
            lon[t]=lon[t]+360
         elif lon[t]>360:
            lon[t]=lon[t]+0 
               
    return lon,lat

 
    
    
def get_Time(dataset):  
 
    timeH = dataset.variables['time'][:] 
      
    # Converting original date into a datetime variable
    tunit = dataset.variables['time'].units
    tunit = tunit.replace('T', ' ')
    tunit = tunit.replace('Z', '')            
    tunit = tunit.replace('hours since ', '')
    tunit = tunit.replace('days since ', '')
    tunit = tunit.replace('minutes since ','')
    dt = parser.parse(tunit) 

    return timeH,dt 

def HYvelocities(dataset):
    #Input the Hycom dataset that has been read.
    #Outputs the water_u and water_v of the file      
    
    uHy = dataset.variables['water_u'][0, 0, :, :]
    vHy = dataset.variables['water_v'][0, 0, :, :]    
    
    return uHy,vHy    


              

def Interpwaters(lonHF, latHF,lonHy,latHy, xHy):
    #input:longitude of HF, Latitude of HF, Longigtude of HY, Latitude of HY, Water from Hycom(either U or V) 
    #Code that interpolates longitudes and latitudes of both
    #High freq and Hycom to output the coresponding
    #Hycom_u and Hycom_v values
    #Output one large array.   
    
  X=[]
  
   
  for l in range(len(latHF)):
     for m in range(len(lonHF)):   
         indLat1 = np.searchsorted(latHy, latHF[l]) #searches for the closes latitude index in Hycom that is greater than the HF latitude
         indLat2 = indLat1-1 #Gets the index of the value right below the first index
    
                
         indLon1 = np.searchsorted(lonHy, lonHF[m])
         indLon2 = indLon1-1
       
    
         gslat = latHy[indLat1]-latHy[indLat2] #difference between the latitude values in Hycom
         gslon = lonHy[indLon1]-lonHy[indLon2] #difference between the longitude values in Hycom
         dlat = abs((latHF[l]-latHy[indLat1])/gslat) #Calculates the percent difference between HFlat and Hycom lat
         dlon = abs((lonHF[m]-lonHy[indLon1])/gslon) #Calculates the percent difference between HFlon and Hycom lon

    #Gets the waters from Hycom        
         LL = xHy[indLat2,indLon2]
         LR = xHy[indLat2,indLon1]
         UL = xHy[indLat1,indLon2]
         UR = xHy[indLat1,indLon1]
        
    #Bilinear Interpolation        
         L1 = LL + dlat*(UL - LL)
         L2 = LR + dlat*(UR - LR)
         rvalue = L1 + dlon*(L2-L1)
         X.append(rvalue)
        
  return X


def Hycomname(date,hour):
 #Input date and hour wanted
 #Hycom File opener. Extracts and opens file for reading. 
 #Outputs filename and dataset(Hycom)  
      nameparts={}
      nameparts["date"]=str(date)
      nameparts["Hours"]=str(hour)
#Below is basically string templates of the files in the Hycom Directory
#it will change depending on input.
      length=len(hour)
      if length==1:
          try:
              file_template = 'hycom_glb_sfc_{date}00_t00{Hours}.nc'
              file = file_template.format(**nameparts)
              Hycom = Dataset(file)
              return file, Hycom
          except OSError as err:
              file_template = 'hycom_glb_sfc_{date}00_t00{Hours}.nc.gz'
              file = file_template.format(**nameparts)
              file = patoolib.extract_archive(file) 
              Hycom = Dataset(file)
              
              return file, Hycom
      elif length==2:
          try:
              file_template = 'hycom_glb_sfc_{date}00_t0{Hours}.nc'
              file = file_template.format(**nameparts)
              Hycom = Dataset(file)
              return file, Hycom
          except OSError as err:
              file_template = 'hycom_glb_sfc_{date}00_t0{Hours}.nc.gz'
              file = file_template.format(**nameparts)
              file = patoolib.extract_archive(file)
              Hycom = Dataset(file)
                            
              return file, Hycom      
      else:
          try:
              file_template = 'hycom_glb_sfc_{date}00_t{Hours}.nc'
              file = file_template.format(**nameparts)
              Hycom = Dataset(file)
              return file, Hycom
          except OSError as err:
              file_template = 'hycom_glb_sfc_{date}00_t{Hours}.nc.gz'
              file = file_template.format(**nameparts)
              file = patoolib.extract_archive(file)
              Hycom = Dataset(file)
                              
              return file, Hycom
          
def HycomVarGet(date,hour):
#Input date and hour
#outputs uHy=Hycom velocity u, vHy=Hycom velocity v, Hycom logitude and latitude. 
    try: 
       file,Hycom=Hycomname(date,hour) #these are functions above.
       uHy,vHy=HYvelocities(Hycom)
       lonHy,latHy=lats_lons(Hycom)
       Hycom.close
       
    except:
       uHy=1
       vHy=1
       lonHy=1
       latHy=1
    return uHy,vHy,lonHy,latHy

            

def NCwriteHytoHF(fileName,U,V):
    #Writes in Hycom U and V into the given HF file
#
    nc = Dataset(fileName,'a',format="NETCDF3_CLASSIC")
    
    nc.createVariable('HYCOM_3Hr_u','float',('time','lat','lon',))
    nc.createVariable('HYCOM_3Hr_v','float',('time','lat','lon',))
    
    
    
    nc.variables['HYCOM_3Hr_u'][:]=U[:]
    nc.variables['HYCOM_3Hr_u'].units='m/s'
    nc.variables['HYCOM_3Hr_u'].long_name='Hycom_eastward_sea_water_velocity'
    
    nc.variables['HYCOM_3Hr_v'][:]=V[:]
    nc.variables['HYCOM_3Hr_v'].units='m/s'
    nc.variables['HYCOM_3Hr_v'].long_name='Hycom_northward_sea_water_velocity'          
    nc.close()                


def TimeInterp(preU,posU,pretime,postime,intHour):
 U=[]    
 #Time interpolation for prewaters and post waters.
 #input U or V from both file before target time and file after target time as preU and and posU 
 #input hour before of file before target time and after target time as pretime and postime 
 #input target time as inthour
 #outputs an array the length of preU or posU of interpolated values U
 
 
 for f in range(len(posU)):
        
     uValue= preU[f]+(intHour-pretime)*((posU[f]-preU[f])/(postime-pretime))
     U.append(uValue)
         
     
 return U


def TrackStd(fileName):    
    q=1 
    t0 = time.clock() #using to test for time, not related to function
    
    #initializes below variables as dictionaries
    doubles = {}#Will store the hours of the hycom files. 0-168 multiples of 3
    uHy = {}#will store hycom u 
    vHy = {}#will store hycom v 
    U = {}#will store interpolated values u
    V = {}#will store interpolated values u
    
    
    for x in range(0, 57):

      doubles[x] = str(x*3) #Hour storage
      #filling zeros to lengthen array
      uHy[x]=0
      vHy[x]=0
      U[x]=0
      V[x]=0    
        
        
        
        
      
    if q==1:

    
        
    
    
        # Reads in High Frequency Data
        nc = Dataset(fileName,'a')
        
        
        lonHF,latHF=lats_lons(nc)#getting lons and lats of inputed nc file
        timeHF,dt=get_Time(nc) #getting time of inputed nc file 
        #timeDF=timeDF.T
        
        #Creating empty arrays to store the standard deviations
        stdu=[]
        std2u=[]
        std3u=[]
        std4u=[]
        std5u=[]
        std6u=[]
        std7u=[]
        
        stdv=[]
        std2v=[]
        std3v=[]
        std4v=[]
        std5v=[]
        std6v=[]
        std7v=[]     
        
        #a dictionary to store the dates
        pro ={}
        pro["date2"]='0'
        for i in range(len(lonHF)):
        #for i in range(2864):
            percentage=((i+1)/len(lonHF))*100 #percentage of code ran
    
            #Initial date parsing
            if dt==datetime.datetime(1979,1,1,0,0): 
                dt2=datetime.datetime.min
                
                offsetD=datetime.timedelta(days=1)
                dateDF=dt2+datetime.timedelta(days=timeHF[0,i])
                dateDF=dateDF-offsetD
                dateDF=dateDF-relativedelta(years=1)
    
            elif dt==datetime.datetime(1900,1,1,0,0):
                dateDF=dt+datetime.timedelta(days=timeHF[i])
                
            elif dt==datetime.datetime(1980,1,1,0,0):
                dateDF=dt+datetime.timedelta(minutes=timeHF[i])
    
            
            pro["date"]=dateDF.strftime('%Y%m%d')#Converting date into a string
    
            #Will run the code below if the date is not the same as previous date.
            if int(pro["date"]) != int(pro["date2"]):
    
    
                try:
                    #goes to directed path, will return nans if failed to find it
                    path_template = r'S:\SOCCER\Maps\HYCOM3Hr\{date}'
                    # path_template= r'S:\SOCCER\WorkingDirectory-QZH\{date}'
                    path = path_template.format(**pro)
                    os.chdir(path)
                    
                    #Opens all 57 files in the target directory.
                    for p in range(len(doubles)):
                        uHy[p],vHy[p],lonH,latH=HycomVarGet(pro["date"],doubles[p])
                        if type(lonH)!=int:
                            lonHy=lonH
                            latHy=latH
    

    
                    #Loop to get the U and V values for target lon and lat    
                    for d in range(len(U)):
    #                for d in range(2):    
                        if type(uHy[d])!=int:
    #                        print("beep")
                            u=Interpwaters(lonHF[i], latHF[i],lonHy,latHy, uHy[d])
                            v=Interpwaters(lonHF[i], latHF[i],lonHy,latHy, vHy[d])
    #                        lengU=len(U[d])
                            U[d]=u[0]
                            V[d]=v[0]
                        else:
                            
                            U[d]=math.nan
                            V[d]=math.nan
    
    #                        
                    
                    #Calculates the standard deviations each day.    
                    stdu.append(np.std([U[0],U[1],U[2],U[3],U[4],U[5],U[6],U[7],U[8]],dtype=np.float))
                    stdv.append(np.std([V[0],V[1],V[2],V[3],V[4],V[5],V[6],V[7],V[8]],dtype=np.float))   
                    
                    std2u.append(np.std([U[0],U[1],U[2],U[3],U[4],U[5],U[6],U[7],U[8],U[9],U[10],U[11],U[12],U[13],U[14],U[15],U[16]],dtype=np.float))
                    std2v.append(np.std([V[0],V[1],V[2],V[3],V[4],V[5],V[6],V[7],V[8],V[9],V[10],V[11],V[12],V[13],V[14],V[15],V[16]],dtype=np.float))    
                    
                    std3u.append(np.std([U[0],U[1],U[2],U[3],U[4],U[5],U[6],U[7],U[8],U[9],U[10],U[11],U[12],U[13],U[14],U[15],U[16],U[17],U[18],U[19],U[20],U[21],U[22],U[23],U[24]],dtype=np.float))
                    std3v.append(np.std([V[0],V[1],V[2],V[3],V[4],V[5],V[6],V[7],V[8],V[9],V[10],V[11],V[12],V[13],V[14],V[15],V[16],V[17],V[18],V[19],V[20],V[21],V[22],V[23],V[24]],dtype=np.float))
                      
                    std4u.append(np.std([U[0],U[1],U[2],U[3],U[4],U[5],U[6],U[7],U[8],U[9],U[10],U[11],U[12],U[13],U[14],U[15],U[16],U[17],U[18],U[19],U[20],U[21],U[22],U[23],U[24],U[25],U[26],U[27],U[28],U[29],U[30],U[31],U[32]],dtype=np.float))
                    std4v.append(np.std([V[0],V[1],V[2],V[3],V[4],V[5],V[6],V[7],V[8],V[9],V[10],V[11],V[12],V[13],V[14],V[15],V[16],V[17],V[18],V[19],V[20],V[21],V[22],V[23],V[24],V[25],V[26],V[27],V[28],V[29],V[30],V[31],V[32]],dtype=np.float))
    
                    std5u.append(np.std([U[0],U[1],U[2],U[3],U[4],U[5],U[6],U[7],U[8],U[9],U[10],U[11],U[12],U[13],U[14],U[15],U[16],U[17],U[18],U[19],U[20],U[21],U[22],U[23],U[24],U[25],U[26],U[27],U[28],U[29],U[30],U[31],U[32],U[33],U[34],U[35],U[36],U[37],U[38],U[39],U[40]],dtype=np.float))
                    std5v.append(np.std([V[0],V[1],V[2],V[3],V[4],V[5],V[6],V[7],V[8],V[9],V[10],V[11],V[12],V[13],V[14],V[15],V[16],V[17],V[18],V[19],V[20],V[21],V[22],V[23],V[24],V[25],V[26],V[27],V[28],V[29],V[30],V[31],V[32],V[33],V[34],V[35],V[36],V[37],V[38],V[39],V[40]],dtype=np.float))
    
                    std6u.append(np.std([U[0],U[1],U[2],U[3],U[4],U[5],U[6],U[7],U[8],U[9],U[10],U[11],U[12],U[13],U[14],U[15],U[16],U[17],U[18],U[19],U[20],U[21],U[22],U[23],U[24],U[25],U[26],U[27],U[28],U[29],U[30],U[31],U[32],U[33],U[34],U[35],U[36],U[37],U[38],U[39],U[40],U[41],U[42],U[43],U[44],U[45],U[46],U[47],U[48]],dtype=np.float))
                    std6v.append(np.std([V[0],V[1],V[2],V[3],V[4],V[5],V[6],V[7],V[8],V[9],V[10],V[11],V[12],V[13],V[14],V[15],V[16],V[17],V[18],V[19],V[20],V[21],V[22],V[23],V[24],V[25],V[26],V[27],V[28],V[29],V[30],V[31],V[32],V[33],V[34],V[35],V[36],V[37],V[38],V[39],V[40],V[41],V[42],V[43],V[44],V[45],V[46],V[47],V[48]],dtype=np.float))
    
                    std7u.append(np.std([U[0],U[1],U[2],U[3],U[4],U[5],U[6],U[7],U[8],U[9],U[10],U[11],U[12],U[13],U[14],U[15],U[16],U[17],U[18],U[19],U[20],U[21],U[22],U[23],U[24],U[25],U[26],U[27],U[28],U[29],U[30],U[31],U[32],U[33],U[34],U[35],U[36],U[37],U[38],U[39],U[40],U[41],U[42],U[43],U[44],U[45],U[46],U[47],U[48],U[49],U[50],U[51],U[52],U[53],U[54],U[55],U[56]],dtype=np.float))
                    std7v.append(np.std([V[0],V[1],V[2],V[3],V[4],V[5],V[6],V[7],V[8],V[9],V[10],V[11],V[12],V[13],V[14],V[15],V[16],V[17],V[18],V[19],V[20],V[21],V[22],V[23],V[24],V[25],V[26],V[27],V[28],V[29],V[30],V[31],V[32],V[33],V[34],V[35],V[36],V[37],V[38],V[39],V[40],V[41],V[42],V[43],V[44],V[45],V[46],V[47],V[48],V[49],V[50],V[51],V[52],V[53],V[54],V[55],V[56]],dtype=np.float))                            
    
    
    
                except FileNotFoundError as err:
                    
                                                 
                    
                    stdu.append(math.nan)
                    std2u.append(math.nan)
                    std3u.append(math.nan)
                    std4u.append(math.nan)
                    std5u.append(math.nan)
                    std6u.append(math.nan)
                    std7u.append(math.nan)
                    
                    stdv.append(math.nan)
                    std2v.append(math.nan)
                    std3v.append(math.nan)
                    std4v.append(math.nan)
                    std5v.append(math.nan)
                    std6v.append(math.nan)
                    std7v.append(math.nan)        
                    print('no file found')
                    
       
    
                
                
    
                
            #if the next point is the same day, then it uses the data collected before.        
            else:
                try:
                    path_template = r'S:\SOCCER\Maps\HYCOM3Hr\{date}'
                    # path_template= r'S:\SOCCER\WorkingDirectory-QZH\{date}'
                    path = path_template.format(**pro)
                    os.chdir(path) 
                    
                    
                    for d in range(len(U)):
    #                for d in range(2):    
                        if type(uHy[d])!=int:
    #                        print("beep")
                            u=Interpwaters(lonHF[i], latHF[i],lonHy,latHy, uHy[d])
                            v=Interpwaters(lonHF[i], latHF[i],lonHy,latHy, vHy[d])
    #                        lengU=len(U[d])
                            U[d]=u[0]
                            V[d]=v[0]
                        
                        else:
                            
                            U[d]=math.nan
                            V[d]=math.nan
    
    #                        
                    stdu.append(np.std([U[0],U[1],U[2],U[3],U[4],U[5],U[6],U[7],U[8]],dtype=np.float))
                    stdv.append(np.std([V[0],V[1],V[2],V[3],V[4],V[5],V[6],V[7],V[8]],dtype=np.float))   
                            
                    std2u.append(np.std([U[0],U[1],U[2],U[3],U[4],U[5],U[6],U[7],U[8],U[9],U[10],U[11],U[12],U[13],U[14],U[15],U[16]],dtype=np.float))
                    std2v.append(np.std([V[0],V[1],V[2],V[3],V[4],V[5],V[6],V[7],V[8],V[9],V[10],V[11],V[12],V[13],V[14],V[15],V[16]],dtype=np.float))    
                    
                    std3u.append(np.std([U[0],U[1],U[2],U[3],U[4],U[5],U[6],U[7],U[8],U[9],U[10],U[11],U[12],U[13],U[14],U[15],U[16],U[17],U[18],U[19],U[20],U[21],U[22],U[23],U[24]],dtype=np.float))
                    std3v.append(np.std([V[0],V[1],V[2],V[3],V[4],V[5],V[6],V[7],V[8],V[9],V[10],V[11],V[12],V[13],V[14],V[15],V[16],V[17],V[18],V[19],V[20],V[21],V[22],V[23],V[24]],dtype=np.float))
                       
                    std4u.append(np.std([U[0],U[1],U[2],U[3],U[4],U[5],U[6],U[7],U[8],U[9],U[10],U[11],U[12],U[13],U[14],U[15],U[16],U[17],U[18],U[19],U[20],U[21],U[22],U[23],U[24],U[25],U[26],U[27],U[28],U[29],U[30],U[31],U[32]],dtype=np.float))
                    std4v.append(np.std([V[0],V[1],V[2],V[3],V[4],V[5],V[6],V[7],V[8],V[9],V[10],V[11],V[12],V[13],V[14],V[15],V[16],V[17],V[18],V[19],V[20],V[21],V[22],V[23],V[24],V[25],V[26],V[27],V[28],V[29],V[30],V[31],V[32]],dtype=np.float))
                        
                    std5u.append(np.std([U[0],U[1],U[2],U[3],U[4],U[5],U[6],U[7],U[8],U[9],U[10],U[11],U[12],U[13],U[14],U[15],U[16],U[17],U[18],U[19],U[20],U[21],U[22],U[23],U[24],U[25],U[26],U[27],U[28],U[29],U[30],U[31],U[32],U[33],U[34],U[35],U[36],U[37],U[38],U[39],U[40]],dtype=np.float))
                    std5v.append(np.std([V[0],V[1],V[2],V[3],V[4],V[5],V[6],V[7],V[8],V[9],V[10],V[11],V[12],V[13],V[14],V[15],V[16],V[17],V[18],V[19],V[20],V[21],V[22],V[23],V[24],V[25],V[26],V[27],V[28],V[29],V[30],V[31],V[32],V[33],V[34],V[35],V[36],V[37],V[38],V[39],V[40]],dtype=np.float))
    
                    std6u.append(np.std([U[0],U[1],U[2],U[3],U[4],U[5],U[6],U[7],U[8],U[9],U[10],U[11],U[12],U[13],U[14],U[15],U[16],U[17],U[18],U[19],U[20],U[21],U[22],U[23],U[24],U[25],U[26],U[27],U[28],U[29],U[30],U[31],U[32],U[33],U[34],U[35],U[36],U[37],U[38],U[39],U[40],U[41],U[42],U[43],U[44],U[45],U[46],U[47],U[48]],dtype=np.float))
                    std6v.append(np.std([V[0],V[1],V[2],V[3],V[4],V[5],V[6],V[7],V[8],V[9],V[10],V[11],V[12],V[13],V[14],V[15],V[16],V[17],V[18],V[19],V[20],V[21],V[22],V[23],V[24],V[25],V[26],V[27],V[28],V[29],V[30],V[31],V[32],V[33],V[34],V[35],V[36],V[37],V[38],V[39],V[40],V[41],V[42],V[43],V[44],V[45],V[46],V[47],V[48]],dtype=np.float))
    
                    std7u.append(np.std([U[0],U[1],U[2],U[3],U[4],U[5],U[6],U[7],U[8],U[9],U[10],U[11],U[12],U[13],U[14],U[15],U[16],U[17],U[18],U[19],U[20],U[21],U[22],U[23],U[24],U[25],U[26],U[27],U[28],U[29],U[30],U[31],U[32],U[33],U[34],U[35],U[36],U[37],U[38],U[39],U[40],U[41],U[42],U[43],U[44],U[45],U[46],U[47],U[48],U[49],U[50],U[51],U[52],U[53],U[54],U[55],U[56]],dtype=np.float))
                    std7v.append(np.std([V[0],V[1],V[2],V[3],V[4],V[5],V[6],V[7],V[8],V[9],V[10],V[11],V[12],V[13],V[14],V[15],V[16],V[17],V[18],V[19],V[20],V[21],V[22],V[23],V[24],V[25],V[26],V[27],V[28],V[29],V[30],V[31],V[32],V[33],V[34],V[35],V[36],V[37],V[38],V[39],V[40],V[41],V[42],V[43],V[44],V[45],V[46],V[47],V[48],V[49],V[50],V[51],V[52],V[53],V[54],V[55],V[56]],dtype=np.float))                            
    
    
                    
                except FileNotFoundError as err:
                    
                                                                 
                    stdu.append(math.nan)
                    std2u.append(math.nan)
                    std3u.append(math.nan)
                    std4u.append(math.nan)
                    std5u.append(math.nan)
                    std6u.append(math.nan)
                    std7u.append(math.nan)
                    
                    stdv.append(math.nan)
                    std2v.append(math.nan)
                    std3v.append(math.nan)
                    std4v.append(math.nan)
                    std5v.append(math.nan)
                    std6v.append(math.nan)
                    std7v.append(math.nan)        
                    print('no file found')    
                
                  
                
                  
                    
            pro["date2"]=dateDF.strftime('%Y%m%d') 
            print('%d percent complete' % percentage)   
    #        pro["Hours"]=str(int(pro["Hours"])+3)
                                   
    
              
    
            
    
        
          
        today= datetime.datetime.today()#sees what time it is today 
        date=today.strftime('%Y%m%d')
    #    os.chdir('S:\SOCCER\OCData\HF-Radar SystemData\HourlyHFData\Hawaii\6km')
        
        #creates varaibles and adds them to file.    
        nc.createVariable('HYCOM_3Hr_u_std_d1','float',('r','c',))
        nc.createVariable('HYCOM_3Hr_u_std_d2','float',('r','c',))
        nc.createVariable('HYCOM_3Hr_u_std_d3','float',('r','c',))
        nc.createVariable('HYCOM_3Hr_u_std_d4','float',('r','c',))
        nc.createVariable('HYCOM_3Hr_u_std_d5','float',('r','c',))
        nc.createVariable('HYCOM_3Hr_u_std_d6','float',('r','c',))
        nc.createVariable('HYCOM_3Hr_u_std_d7','float',('r','c',))
            
            
        nc.createVariable('HYCOM_3Hr_v_std_d1','float',('r','c',))
        nc.createVariable('HYCOM_3Hr_v_std_d2','float',('r','c',))
        nc.createVariable('HYCOM_3Hr_v_std_d3','float',('r','c',))
        nc.createVariable('HYCOM_3Hr_v_std_d4','float',('r','c',))
        nc.createVariable('HYCOM_3Hr_v_std_d5','float',('r','c',))
        nc.createVariable('HYCOM_3Hr_v_std_d6','float',('r','c',))
        nc.createVariable('HYCOM_3Hr_v_std_d7','float',('r','c',))
            
            
        nc.variables['HYCOM_3Hr_u_std_d1'][:,0]=stdu[:]
        nc.variables['HYCOM_3Hr_u_std_d1'].units='m/s'
        nc.variables['HYCOM_3Hr_u_std_d1'].long_name='Hycom_eastward_sea_water_velocity_stardard_deviation_hours_0to24'
        nc.variables['HYCOM_3Hr_u_std_d1'].date_added=date
        nc.variables['HYCOM_3Hr_u_std_d1'].added_by='qzh14'
            
        nc.variables['HYCOM_3Hr_u_std_d2'][:,0]=std2u[:]
        nc.variables['HYCOM_3Hr_u_std_d2'].units='m/s'
        nc.variables['HYCOM_3Hr_u_std_d2'].long_name='Hycom_eastward_sea_water_velocity_stardard_deviation_hours_0to48'
        nc.variables['HYCOM_3Hr_u_std_d2'].date_added=date
        nc.variables['HYCOM_3Hr_u_std_d2'].added_by='qzh14'
        
        nc.variables['HYCOM_3Hr_u_std_d3'][:,0]=std3u[:]
        nc.variables['HYCOM_3Hr_u_std_d3'].units='m/s'
        nc.variables['HYCOM_3Hr_u_std_d3'].long_name='Hycom_eastward_sea_water_velocity_stardard_deviation_hours_0to72'
        nc.variables['HYCOM_3Hr_u_std_d3'].date_added=date
        nc.variables['HYCOM_3Hr_u_std_d3'].added_by='qzh14'    
            
        nc.variables['HYCOM_3Hr_u_std_d4'][:,0]=std4u[:]
        nc.variables['HYCOM_3Hr_u_std_d4'].units='m/s'
        nc.variables['HYCOM_3Hr_u_std_d4'].long_name='Hycom_eastward_sea_water_velocity_stardard_deviation_hours_0to96'
        nc.variables['HYCOM_3Hr_u_std_d4'].date_added=date
        nc.variables['HYCOM_3Hr_u_std_d4'].added_by='qzh14'    
            
        nc.variables['HYCOM_3Hr_u_std_d5'][:,0]=std5u[:]
        nc.variables['HYCOM_3Hr_u_std_d5'].units='m/s'
        nc.variables['HYCOM_3Hr_u_std_d5'].long_name='Hycom_eastward_sea_water_velocity_stardard_deviation_hours_0to120'
        nc.variables['HYCOM_3Hr_u_std_d5'].date_added=date
        nc.variables['HYCOM_3Hr_u_std_d5'].added_by='qzh14'    
        
        nc.variables['HYCOM_3Hr_u_std_d6'][:,0]=std6u[:]
        nc.variables['HYCOM_3Hr_u_std_d6'].units='m/s'
        nc.variables['HYCOM_3Hr_u_std_d6'].long_name='Hycom_eastward_sea_water_velocity_stardard_deviation_hours_0to144'
        nc.variables['HYCOM_3Hr_u_std_d6'].date_added=date
        nc.variables['HYCOM_3Hr_u_std_d6'].added_by='qzh14'    
        
        nc.variables['HYCOM_3Hr_u_std_d7'][:,0]=std7u[:]
        nc.variables['HYCOM_3Hr_u_std_d7'].units='m/s'
        nc.variables['HYCOM_3Hr_u_std_d7'].long_name='Hycom_eastward_sea_water_velocity_stardard_deviation_hours_0to168'
        nc.variables['HYCOM_3Hr_u_std_d7'].date_added=date
        nc.variables['HYCOM_3Hr_u_std_d7'].added_by='qzh14' 
        
        
        nc.variables['HYCOM_3Hr_v_std_d1'][:,0]=stdv[:]
        nc.variables['HYCOM_3Hr_v_std_d1'].units='m/s'
        nc.variables['HYCOM_3Hr_v_std_d1'].long_name='Hycom_northward_sea_water_velocity_stardard_deviation_hours_0to24'
        nc.variables['HYCOM_3Hr_v_std_d1'].date_added=date
        nc.variables['HYCOM_3Hr_v_std_d1'].added_by='qzh14'
            
        nc.variables['HYCOM_3Hr_v_std_d2'][:,0]=std2v[:]
        nc.variables['HYCOM_3Hr_v_std_d2'].units='m/s'
        nc.variables['HYCOM_3Hr_v_std_d2'].long_name='Hycom_northward_sea_water_velocity_stardard_deviation_hours_0to48'
        nc.variables['HYCOM_3Hr_v_std_d2'].date_added=date
        nc.variables['HYCOM_3Hr_v_std_d2'].added_by='qzh14'
        
        nc.variables['HYCOM_3Hr_v_std_d3'][:,0]=std3v[:]
        nc.variables['HYCOM_3Hr_v_std_d3'].units='m/s'
        nc.variables['HYCOM_3Hr_v_std_d3'].long_name='Hycom_northward_sea_water_velocity_stardard_deviation_hours_0to72'
        nc.variables['HYCOM_3Hr_v_std_d3'].date_added=date
        nc.variables['HYCOM_3Hr_v_std_d3'].added_by='qzh14'    
            
        nc.variables['HYCOM_3Hr_v_std_d4'][:,0]=std4v[:]
        nc.variables['HYCOM_3Hr_v_std_d4'].units='m/s'
        nc.variables['HYCOM_3Hr_v_std_d4'].long_name='Hycom_northward_sea_water_velocity_stardard_deviation_hours_0to96'
        nc.variables['HYCOM_3Hr_v_std_d4'].date_added=date
        nc.variables['HYCOM_3Hr_v_std_d4'].added_by='qzh14'    
            
        nc.variables['HYCOM_3Hr_v_std_d5'][:,0]=std5v[:]
        nc.variables['HYCOM_3Hr_v_std_d5'].units='m/s'
        nc.variables['HYCOM_3Hr_v_std_d5'].long_name='Hycom_northward_sea_water_velocity_stardard_deviation_hours_0to120'
        nc.variables['HYCOM_3Hr_v_std_d5'].date_added=date
        nc.variables['HYCOM_3Hr_v_std_d5'].added_by='qzh14'    
        
        nc.variables['HYCOM_3Hr_v_std_d6'][:,0]=std6v[:]
        nc.variables['HYCOM_3Hr_v_std_d6'].units='m/s'
        nc.variables['HYCOM_3Hr_v_std_d6'].long_name='Hycom_northward_sea_water_velocity_stardard_deviation_hours_0to144'
        nc.variables['HYCOM_3Hr_v_std_d6'].date_added=date
        nc.variables['HYCOM_3Hr_v_std_d6'].added_by='qzh14'    
        
        nc.variables['HYCOM_3Hr_v_std_d7'][:,0]=std7v[:]
        nc.variables['HYCOM_3Hr_v_std_d7'].units='m/s'
        nc.variables['HYCOM_3Hr_v_std_d7'].long_name='Hycom_northward_sea_water_velocity_stardard_deviation_hours_0to168'
        nc.variables['HYCOM_3Hr_v_std_d7'].date_added=date
        nc.variables['HYCOM_3Hr_v_std_d7'].added_by='qzh14'

        nc.close()
        t1 = time.clock()
    
        total = t1-t0 
        print('Std data added to %s' %fileName)
        return  total 





def HFstd(fileName):
    q=1
    #Dictionary creations
    doubles = {}
    uHy = {}
    vHy = {}
    U = {}
    V = {}
    
    for x in range(0, 57):
    #  print('%d' %x)  
      doubles[x] = str(x*3) 
      uHy[x]=0
      vHy[x]=0
      U[x]=0
      V[x]=0    
        
        
        
        
        
    if q==1:
        t0 = time.clock()
    
        
    
    
        # Reads in High Frequency Data
        nc = Dataset(fileName,'a')#reads in dataset of HF
        
        
        lonHF,latHF=lats_lons(nc) #outputs lon and lats
        timeHF,dt=get_Time(nc)  #outputs time and initial date time
        #timeDF=timeDF.T
        
        #Initializes grid arrays.
        stduf=np.empty([len(timeHF),len(latHF),len(lonHF)], dtype=float)
        std2uf=np.empty([len(timeHF),len(latHF),len(lonHF)], dtype=float)
        std3uf=np.empty([len(timeHF),len(latHF),len(lonHF)], dtype=float)
        std4uf=np.empty([len(timeHF),len(latHF),len(lonHF)], dtype=float)
        std5uf=np.empty([len(timeHF),len(latHF),len(lonHF)], dtype=float)
        std6uf=np.empty([len(timeHF),len(latHF),len(lonHF)], dtype=float)
        std7uf=np.empty([len(timeHF),len(latHF),len(lonHF)], dtype=float)
        
        stdvf=np.empty([len(timeHF),len(latHF),len(lonHF)], dtype=float)
        std2vf=np.empty([len(timeHF),len(latHF),len(lonHF)], dtype=float)
        std3vf=np.empty([len(timeHF),len(latHF),len(lonHF)], dtype=float)
        std4vf=np.empty([len(timeHF),len(latHF),len(lonHF)], dtype=float)
        std5vf=np.empty([len(timeHF),len(latHF),len(lonHF)], dtype=float)
        std6vf=np.empty([len(timeHF),len(latHF),len(lonHF)], dtype=float)
        std7vf=np.empty([len(timeHF),len(latHF),len(lonHF)], dtype=float)    
        
        pro ={}
        pro["date2"]='0'
        for i in range(len(timeHF)):
        #for i in range(2864):
            percentage=((i+1)/len(timeHF))*100
    
            #getting the date 
            dateHF=dt+datetime.timedelta(hours=timeHF[i])
    
            #converting date into a string
            pro["date"]=dateHF.strftime('%Y%m%d')
    
            
            if int(pro["date"]) != int(pro["date2"]):
                #temporary sandard deviations for the loop
                stdu=[]
                std2u=[]
                std3u=[]
                std4u=[]
                std5u=[]
                std6u=[]
                std7u=[]
                
                stdv=[]
                std2v=[]
                std3v=[]
                std4v=[]
                std5v=[]
                std6v=[]
                std7v=[] 
                try:
                    path_template = r'S:\SOCCER\Maps\HYCOM3Hr\{date}'
                    # path_template= r'S:\SOCCER\WorkingDirectory-QZH\{date}'
                    path = path_template.format(**pro)
                    os.chdir(path)
                    
                    #opens 57 Hycom files
                    for p in range(len(doubles)):
                        uHy[p],vHy[p],lonH,latH=HycomVarGet(pro["date"],doubles[p])
                        if type(lonH)!=int:
                            lonHy=lonH
                            latHy=latH
                        
                    #Interpolates to get U and V values    
                    for d in range(len(U)):
    #                for d in range(2):    
                        if type(uHy[d])!=int:
    #                        
                            U[d]=Interpwaters(lonHF, latHF,lonHy,latHy, uHy[d])
                            V[d]=Interpwaters(lonHF, latHF,lonHy,latHy, vHy[d])
    #                        lengU=len(U[d])
                        
                        else:
                            #if data is not found it creates empty arrays of nan
                            lengU=len(lonHF)*len(latHF)
                            U[d]=np.empty(lengU, dtype=float)
                            V[d]=np.empty(lengU, dtype=float)
                            U[d].fill(math.nan)
                            V[d].fill(math.nan)
    #               
         
                    for r in range(lengU):
                       #standard deviation calculations 
                        stdu.append(np.std([U[0][r],U[1][r],U[2][r],U[3][r],U[4][r],U[5][r],U[6][r],U[7][r],U[8][r]],dtype=np.float))
                        stdv.append(np.std([V[0][r],V[1][r],V[2][r],V[3][r],V[4][r],V[5][r],V[6][r],V[7][r],V[8][r]],dtype=np.float))   
                            
                        std2u.append(np.std([U[0][r],U[1][r],U[2][r],U[3][r],U[4][r],U[5][r],U[6][r],U[7][r],U[8][r],U[9][r],U[10][r],U[11][r],U[12][r],U[13][r],U[14][r],U[15][r],U[16][r]],dtype=np.float))
                        std2v.append(np.std([V[0][r],V[1][r],V[2][r],V[3][r],V[4][r],V[5][r],V[6][r],V[7][r],V[8][r],U[9][r],U[10][r],U[11][r],U[12][r],U[13][r],U[14][r],U[15][r],U[16][r]],dtype=np.float))    
                            
                        std3u.append(np.std([U[0][r],U[1][r],U[2][r],U[3][r],U[4][r],U[5][r],U[6][r],U[7][r],U[8][r],U[9][r],U[10][r],U[11][r],U[12][r],U[13][r],U[14][r],U[15][r],U[16][r],U[17][r],U[18][r],U[19][r],U[20][r],U[21][r],U[22][r],U[23][r],U[24][r]],dtype=np.float))
                        std3v.append(np.std([V[0][r],V[1][r],V[2][r],V[3][r],V[4][r],V[5][r],V[6][r],V[7][r],V[8][r],V[9][r],V[10][r],V[11][r],V[12][r],V[13][r],V[14][r],V[15][r],V[16][r],V[17][r],V[18][r],V[19][r],V[20][r],V[21][r],V[22][r],V[23][r],V[24][r]],dtype=np.float))
                       
                        std4u.append(np.std([U[0][r],U[1][r],U[2][r],U[3][r],U[4][r],U[5][r],U[6][r],U[7][r],U[8][r],U[9][r],U[10][r],U[11][r],U[12][r],U[13][r],U[14][r],U[15][r],U[16][r],U[17][r],U[18][r],U[19][r],U[20][r],U[21][r],U[22][r],U[23][r],U[24][r],U[25][r],U[26][r],U[27][r],U[28][r],U[29][r],U[30][r],U[31][r],U[32][r]],dtype=np.float))
                        std4v.append(np.std([V[0][r],V[1][r],V[2][r],V[3][r],V[4][r],V[5][r],V[6][r],V[7][r],V[8][r],V[9][r],V[10][r],V[11][r],V[12][r],V[13][r],V[14][r],V[15][r],V[16][r],V[17][r],V[18][r],V[19][r],V[20][r],V[21][r],V[22][r],V[23][r],V[24][r],V[25][r],V[26][r],V[27][r],V[28][r],V[29][r],V[30][r],V[31][r],V[32][r]],dtype=np.float))
    
                        std5u.append(np.std([U[0][r],U[1][r],U[2][r],U[3][r],U[4][r],U[5][r],U[6][r],U[7][r],U[8][r],U[9][r],U[10][r],U[11][r],U[12][r],U[13][r],U[14][r],U[15][r],U[16][r],U[17][r],U[18][r],U[19][r],U[20][r],U[21][r],U[22][r],U[23][r],U[24][r],U[25][r],U[26][r],U[27][r],U[28][r],U[29][r],U[30][r],U[31][r],U[32][r],U[33][r],U[34][r],U[35][r],U[36][r],U[37][r],U[38][r],U[39][r],U[40][r]],dtype=np.float))
                        std5v.append(np.std([V[0][r],V[1][r],V[2][r],V[3][r],V[4][r],V[5][r],V[6][r],V[7][r],V[8][r],V[9][r],V[10][r],V[11][r],V[12][r],V[13][r],V[14][r],V[15][r],V[16][r],V[17][r],V[18][r],V[19][r],V[20][r],V[21][r],V[22][r],V[23][r],V[24][r],V[25][r],V[26][r],V[27][r],V[28][r],V[29][r],V[30][r],V[31][r],V[32][r],V[33][r],V[34][r],V[35][r],V[36][r],V[37][r],V[38][r],V[39][r],V[40][r]],dtype=np.float))
    
                        std6u.append(np.std([U[0][r],U[1][r],U[2][r],U[3][r],U[4][r],U[5][r],U[6][r],U[7][r],U[8][r],U[9][r],U[10][r],U[11][r],U[12][r],U[13][r],U[14][r],U[15][r],U[16][r],U[17][r],U[18][r],U[19][r],U[20][r],U[21][r],U[22][r],U[23][r],U[24][r],U[25][r],U[26][r],U[27][r],U[28][r],U[29][r],U[30][r],U[31][r],U[32][r],U[33][r],U[34][r],U[35][r],U[36][r],U[37][r],U[38][r],U[39][r],U[40][r],U[41][r],U[42][r],U[43][r],U[44][r],U[45][r],U[46][r],U[47][r],U[48][r]],dtype=np.float))
                        std6v.append(np.std([V[0][r],V[1][r],V[2][r],V[3][r],V[4][r],V[5][r],V[6][r],V[7][r],V[8][r],V[9][r],V[10][r],V[11][r],V[12][r],V[13][r],V[14][r],V[15][r],V[16][r],V[17][r],V[18][r],V[19][r],V[20][r],V[21][r],V[22][r],V[23][r],V[24][r],V[25][r],V[26][r],V[27][r],V[28][r],V[29][r],V[30][r],V[31][r],V[32][r],V[33][r],V[34][r],V[35][r],V[36][r],V[37][r],V[38][r],V[39][r],V[40][r],V[41][r],V[42][r],V[43][r],V[44][r],V[45][r],V[46][r],V[47][r],V[48][r]],dtype=np.float))
    
                        std7u.append(np.std([U[0][r],U[1][r],U[2][r],U[3][r],U[4][r],U[5][r],U[6][r],U[7][r],U[8][r],U[9][r],U[10][r],U[11][r],U[12][r],U[13][r],U[14][r],U[15][r],U[16][r],U[17][r],U[18][r],U[19][r],U[20][r],U[21][r],U[22][r],U[23][r],U[24][r],U[25][r],U[26][r],U[27][r],U[28][r],U[29][r],U[30][r],U[31][r],U[32][r],U[33][r],U[34][r],U[35][r],U[36][r],U[37][r],U[38][r],U[39][r],U[40][r],U[41][r],U[42][r],U[43][r],U[44][r],U[45][r],U[46][r],U[47][r],U[48][r],U[49][r],U[50][r],U[51][r],U[52][r],U[53][r],U[54][r],U[55][r],U[56][r]],dtype=np.float))
                        std7v.append(np.std([V[0][r],V[1][r],V[2][r],V[3][r],V[4][r],V[5][r],V[6][r],V[7][r],V[8][r],V[9][r],V[10][r],V[11][r],V[12][r],V[13][r],V[14][r],V[15][r],V[16][r],V[17][r],V[18][r],V[19][r],V[20][r],V[21][r],V[22][r],V[23][r],V[24][r],V[25][r],V[26][r],V[27][r],V[28][r],V[29][r],V[30][r],V[31][r],V[32][r],V[33][r],V[34][r],V[35][r],V[36][r],V[37][r],V[38][r],V[39][r],V[40][r],V[41][r],V[42][r],V[43][r],V[44][r],V[45][r],V[46][r],V[47][r],V[48][r],V[49][r],V[50][r],V[51][r],V[52][r],V[53][r],V[54][r],V[55][r],V[56][r]],dtype=np.float))                
    
                except FileNotFoundError as err:
                    #if directory is not present, it makes the entire grid nan
                    stdu=np.empty([len(latHF),len(lonHF)], dtype=float)
                    std2u=np.empty([len(latHF),len(lonHF)], dtype=float)
                    std3u=np.empty([len(latHF),len(lonHF)], dtype=float)
                    std4u=np.empty([len(latHF),len(lonHF)], dtype=float)
                    std5u=np.empty([len(latHF),len(lonHF)], dtype=float)
                    std6u=np.empty([len(latHF),len(lonHF)], dtype=float)
                    std7u=np.empty([len(latHF),len(lonHF)], dtype=float)
                    
                    stdv=np.empty([len(latHF),len(lonHF)], dtype=float)
                    std2v=np.empty([len(latHF),len(lonHF)], dtype=float)
                    std3v=np.empty([len(latHF),len(lonHF)], dtype=float)
                    std4v=np.empty([len(latHF),len(lonHF)], dtype=float)
                    std5v=np.empty([len(latHF),len(lonHF)], dtype=float)
                    std6v=np.empty([len(latHF),len(lonHF)], dtype=float)
                    std7v=np.empty([len(latHF),len(lonHF)], dtype=float)                
                    
                    
                    
                    stdu.fill(math.nan)
                    std2u.fill(math.nan)
                    std3u.fill(math.nan)
                    std4u.fill(math.nan)
                    std5u.fill(math.nan)
                    std6u.fill(math.nan)
                    std7u.fill(math.nan)
                    
                    stdv.fill(math.nan)
                    std2v.fill(math.nan)
                    std3v.fill(math.nan)
                    std4v.fill(math.nan)
                    std5v.fill(math.nan)
                    std6v.fill(math.nan)
                    std7v.fill(math.nan)        
                    print('no file found')
                    
                  #turning vector into a numpy vector 
                stdu=np.asarray(stdu)
                std2u=np.asarray(std2u)
                std3u=np.asarray(std3u)
                std4u=np.asarray(std4u)
                std5u=np.asarray(std5u)
                std6u=np.asarray(std6u)
                std7u=np.asarray(std7u)
            
                stdv=np.asarray(stdu)
                std2v=np.asarray(std2u)
                std3v=np.asarray(std3u)
                std4v=np.asarray(std4u)
                std5v=np.asarray(std5u)
                std6v=np.asarray(std6u)
                std7v=np.asarray(std7u)
                
                
                #turning numpy vector into a numpy matrix
                stdu=stdu.reshape(len(latHF),len(lonHF))
                std2u=stdu.reshape(len(latHF),len(lonHF))
                std3u=stdu.reshape(len(latHF),len(lonHF))
                std4u=stdu.reshape(len(latHF),len(lonHF))
                std5u=stdu.reshape(len(latHF),len(lonHF))
                std6u=stdu.reshape(len(latHF),len(lonHF))
                std7u=stdu.reshape(len(latHF),len(lonHF))
            
                stdv=stdu.reshape(len(latHF),len(lonHF))
                std2v=stdu.reshape(len(latHF),len(lonHF))
                std3v=stdu.reshape(len(latHF),len(lonHF))
                std4v=stdu.reshape(len(latHF),len(lonHF))
                std5v=stdu.reshape(len(latHF),len(lonHF))
                std6v=stdu.reshape(len(latHF),len(lonHF))
                std7v=stdu.reshape(len(latHF),len(lonHF))
                
                #adding data to array on current time
                stduf[i,:,:]=stdu
                std2uf[i,:,:]=std2u
                std3uf[i,:,:]=std3u
                std4uf[i,:,:]=std4u
                std5uf[i,:,:]=std5u
                std6uf[i,:,:]=std6u
                std7uf[i,:,:]=std7u
        
                stdvf[i,:,:]=stdv
                std2vf[i,:,:]=std2v
                std3vf[i,:,:]=std3v
                std4vf[i,:,:]=std4v
                std5vf[i,:,:]=std5v
                std6vf[i,:,:]=std6v
                std7vf[i,:,:]=std7v                
                    
            else:
                
                stduf[i,:,:]=stdu
                std2uf[i,:,:]=std2u
                std3uf[i,:,:]=std3u
                std4uf[i,:,:]=std4u
                std5uf[i,:,:]=std5u
                std6uf[i,:,:]=std6u
                std7uf[i,:,:]=std7u
        
                stdvf[i,:,:]=stdv
                std2vf[i,:,:]=std2v
                std3vf[i,:,:]=std3v
                std4vf[i,:,:]=std4v
                std5vf[i,:,:]=std5v
                std6vf[i,:,:]=std6v
                std7vf[i,:,:]=std7v              
                
                    
                    
            pro["date2"]=dateHF.strftime('%Y%m%d') 
            print('%d percent complete' % percentage)   
    #        pro["Hours"]=str(int(pro["Hours"])+3)
                                   
    
              
    
            
    
        
      
        today= datetime.datetime.today() #today's date
        date=today.strftime('%Y%m%d')
    #    os.chdir('S:\SOCCER\OCData\HF-Radar SystemData\HourlyHFData\Hawaii\6km')
        
        #creating variables and adding it the files    
        nc.createVariable('HYCOM_3Hr_u_std_d1','float',('time','lat','lon',))
        nc.createVariable('HYCOM_3Hr_u_std_d2','float',('time','lat','lon',))
        nc.createVariable('HYCOM_3Hr_u_std_d3','float',('time','lat','lon',))
        nc.createVariable('HYCOM_3Hr_u_std_d4','float',('time','lat','lon',))
        nc.createVariable('HYCOM_3Hr_u_std_d5','float',('time','lat','lon',))
        nc.createVariable('HYCOM_3Hr_u_std_d6','float',('time','lat','lon',))
        nc.createVariable('HYCOM_3Hr_u_std_d7','float',('time','lat','lon',))
            
            
        nc.createVariable('HYCOM_3Hr_v_std_d1','float',('time','lat','lon',))
        nc.createVariable('HYCOM_3Hr_v_std_d2','float',('time','lat','lon',))
        nc.createVariable('HYCOM_3Hr_v_std_d3','float',('time','lat','lon',))
        nc.createVariable('HYCOM_3Hr_v_std_d4','float',('time','lat','lon',))
        nc.createVariable('HYCOM_3Hr_v_std_d5','float',('time','lat','lon',))
        nc.createVariable('HYCOM_3Hr_v_std_d6','float',('time','lat','lon',))
        nc.createVariable('HYCOM_3Hr_v_std_d7','float',('time','lat','lon',))
            
            
        nc.variables['HYCOM_3Hr_u_std_d1'][:]=stduf[:]
        nc.variables['HYCOM_3Hr_u_std_d1'].units='m/s'
        nc.variables['HYCOM_3Hr_u_std_d1'].long_name='Hycom_northward_sea_water_velocity_stardard_deviation_hours_0to24'
        nc.variables['HYCOM_3Hr_u_std_d1'].date_added=date
        nc.variables['HYCOM_3Hr_u_std_d1'].added_by='qzh14'
            
        nc.variables['HYCOM_3Hr_u_std_d2'][:]=std2uf[:]
        nc.variables['HYCOM_3Hr_u_std_d2'].units='m/s'
        nc.variables['HYCOM_3Hr_u_std_d2'].long_name='Hycom_northward_sea_water_velocity_stardard_deviation_hours_0to48'
        nc.variables['HYCOM_3Hr_u_std_d2'].date_added=date
        nc.variables['HYCOM_3Hr_u_std_d2'].added_by='qzh14'
        
        nc.variables['HYCOM_3Hr_u_std_d3'][:]=std3uf[:]
        nc.variables['HYCOM_3Hr_u_std_d3'].units='m/s'
        nc.variables['HYCOM_3Hr_u_std_d3'].long_name='Hycom_northward_sea_water_velocity_stardard_deviation_hours_0to72'
        nc.variables['HYCOM_3Hr_u_std_d3'].date_added=date
        nc.variables['HYCOM_3Hr_u_std_d3'].added_by='qzh14'    
            
        nc.variables['HYCOM_3Hr_u_std_d4'][:]=std4uf[:]
        nc.variables['HYCOM_3Hr_u_std_d4'].units='m/s'
        nc.variables['HYCOM_3Hr_u_std_d4'].long_name='Hycom_northward_sea_water_velocity_stardard_deviation_hours_0to96'
        nc.variables['HYCOM_3Hr_u_std_d4'].date_added=date
        nc.variables['HYCOM_3Hr_u_std_d4'].added_by='qzh14'    
            
        nc.variables['HYCOM_3Hr_u_std_d5'][:]=std5uf[:]
        nc.variables['HYCOM_3Hr_u_std_d5'].units='m/s'
        nc.variables['HYCOM_3Hr_u_std_d5'].long_name='Hycom_northward_sea_water_velocity_stardard_deviation_hours_0to120'
        nc.variables['HYCOM_3Hr_u_std_d5'].date_added=date
        nc.variables['HYCOM_3Hr_u_std_d5'].added_by='qzh14'    
        
        nc.variables['HYCOM_3Hr_u_std_d6'][:]=std6uf[:]
        nc.variables['HYCOM_3Hr_u_std_d6'].units='m/s'
        nc.variables['HYCOM_3Hr_u_std_d6'].long_name='Hycom_northward_sea_water_velocity_stardard_deviation_hours_0to144'
        nc.variables['HYCOM_3Hr_u_std_d6'].date_added=date
        nc.variables['HYCOM_3Hr_u_std_d6'].added_by='qzh14'    
        
        nc.variables['HYCOM_3Hr_u_std_d7'][:]=std7uf[:]
        nc.variables['HYCOM_3Hr_u_std_d7'].units='m/s'
        nc.variables['HYCOM_3Hr_u_std_d7'].long_name='Hycom_northward_sea_water_velocity_stardard_deviation_hours_0to168'
        nc.variables['HYCOM_3Hr_u_std_d7'].date_added=date
        nc.variables['HYCOM_3Hr_u_std_d7'].added_by='qzh14' 
        
        
        nc.variables['HYCOM_3Hr_v_std_d1'][:]=stdvf[:]
        nc.variables['HYCOM_3Hr_v_std_d1'].units='m/s'
        nc.variables['HYCOM_3Hr_v_std_d1'].long_name='Hycom_northward_sea_water_velocity_stardard_deviation_hours_0to24'
        nc.variables['HYCOM_3Hr_v_std_d1'].date_added=date
        nc.variables['HYCOM_3Hr_v_std_d1'].added_by='qzh14'
            
        nc.variables['HYCOM_3Hr_v_std_d2'][:]=std2vf[:]
        nc.variables['HYCOM_3Hr_v_std_d2'].units='m/s'
        nc.variables['HYCOM_3Hr_v_std_d2'].long_name='Hycom_northward_sea_water_velocity_stardard_deviation_hours_0to48'
        nc.variables['HYCOM_3Hr_v_std_d2'].date_added=date
        nc.variables['HYCOM_3Hr_v_std_d2'].added_by='qzh14'
        
        nc.variables['HYCOM_3Hr_v_std_d3'][:]=std3vf[:]
        nc.variables['HYCOM_3Hr_v_std_d3'].units='m/s'
        nc.variables['HYCOM_3Hr_v_std_d3'].long_name='Hycom_northward_sea_water_velocity_stardard_deviation_hours_0to72'
        nc.variables['HYCOM_3Hr_v_std_d3'].date_added=date
        nc.variables['HYCOM_3Hr_v_std_d3'].added_by='qzh14'    
            
        nc.variables['HYCOM_3Hr_v_std_d4'][:]=std4vf[:]
        nc.variables['HYCOM_3Hr_v_std_d4'].units='m/s'
        nc.variables['HYCOM_3Hr_v_std_d4'].long_name='Hycom_northward_sea_water_velocity_stardard_deviation_hours_0to96'
        nc.variables['HYCOM_3Hr_v_std_d4'].date_added=date
        nc.variables['HYCOM_3Hr_v_std_d4'].added_by='qzh14'    
            
        nc.variables['HYCOM_3Hr_v_std_d5'][:]=std5vf[:]
        nc.variables['HYCOM_3Hr_v_std_d5'].units='m/s'
        nc.variables['HYCOM_3Hr_v_std_d5'].long_name='Hycom_northward_sea_water_velocity_stardard_deviation_hours_0to120'
        nc.variables['HYCOM_3Hr_v_std_d5'].date_added=date
        nc.variables['HYCOM_3Hr_v_std_d5'].added_by='qzh14'    
        
        nc.variables['HYCOM_3Hr_v_std_d6'][:]=std6vf[:]
        nc.variables['HYCOM_3Hr_v_std_d6'].units='m/s'
        nc.variables['HYCOM_3Hr_v_std_d6'].long_name='Hycom_northward_sea_water_velocity_stardard_deviation_hours_0to144'
        nc.variables['HYCOM_3Hr_v_std_d6'].date_added=date
        nc.variables['HYCOM_3Hr_v_std_d6'].added_by='qzh14'    
        
        nc.variables['HYCOM_3Hr_v_std_d7'][:]=std7vf[:]
        nc.variables['HYCOM_3Hr_v_std_d7'].units='m/s'
        nc.variables['HYCOM_3Hr_v_std_d7'].long_name='Hycom_northward_sea_water_velocity_stardard_deviation_hours_0to168'
        nc.variables['HYCOM_3Hr_v_std_d7'].date_added=date
        nc.variables['HYCOM_3Hr_v_std_d7'].added_by='qzh14'
    #    
    #    
    #    
    #    
    #    
    #       
    #        
    #         
        nc.close()
        t1 = time.clock()
    
        total = t1-t0 
        
        return  total #returns time taken to run code



def AddHycomtoHF(fileName):
    tdelta=0
    t0 = time.clock()
    # Reads in High Frequency Data
    nc = Dataset(fileName,'r')
    
    
    
    
    # Reading in Longitude, Latitude, and Time
    
    #timeHF = nc.variables['time'][:] 
    
    lonHF,latHF = lats_lons(nc)
    timeHF,dt = get_Time(nc)

    
    #creates an empty array to store the Hycom data
    U2=np.empty([len(timeHF),len(latHF),len(lonHF)], dtype=float)
    V2=np.empty([len(timeHF),len(latHF),len(lonHF)], dtype=float)
    
    
    
    
    for i in range(len(timeHF)):
        percentage=((i+1)/len(timeHF))*100
        U=[]
        V=[]    
    
        #Adds hours to the Original date

        cro=dt+datetime.timedelta(hours=timeHF[i])
    
        #dictionary to store date and Hours
        pro ={}
        pro["date"]=cro.strftime('%Y%m%d')
        pro["Hours"] = cro.strftime('%H')
        
        #creating template for path
        path_template = r'S:\SOCCER\Maps\HYCOM3Hr\{date}'

        path = path_template.format(**pro)
        
        #changes path
        os.chdir(path)
     
    
        intHour = int(pro["Hours"])
        
        
        #Code that gets the U and V values if the hour is a multiple of 3 or equal to 0
        if intHour % 3 == 0 or intHour == 0:
    
            #uses the hycomename function to read in data for date and time.     
            try:
                file,Hycom = Hycomname(pro["date"],pro["Hours"])
    
    #       if above statement failes, code below looks back up to 7 days until it finds a prediction. 
            except :
               tdelta=1 #days going back variable
               revtime = intHour+24 #adding 24 hours to hour variable
               while (tdelta<=7):
                   try:
                       
                       pro["revHours"] = str(revtime)
                       predate = cro - datetime.timedelta(days=tdelta)
                       pro["revdate"]=predate.strftime('%Y%m%d')
                       prepath_template = r'S:\SOCCER\Maps\HYCOM3Hr\{revdate}'
    #                   prepath_template = r'S:\SOCCER\WorkingDirectory-QZH\{revdate}'
                       prepath = prepath_template.format(**pro)
                       os.chdir(prepath)
                       
                       revfile,Hycom = Hycomname(pro["revdate"],pro["revHours"])
    
                       break
    #               except OSError as err:
                   except :
                       tdelta = tdelta+1
                       revtime = revtime + 24
            
            #Reads in Hycom Values
            
            lonHy,latHy = lats_lons(Hycom)
            
    
            uHy,vHy = HYvelocities(Hycom)
    
    #Loop for interpolating u values in Hycom    
            U = Interpwaters(lonHF,latHF, lonHy, latHy, uHy)
            V = Interpwaters(lonHF,latHF, lonHy, latHy, vHy)
                
                    
#Code that gets the U and V values if the hour isn't a multiple of 3 or equal to 0       
    
        else :   

            
            pretime = intHour-2 #Hour before target time
            if pretime % 3 !=0:
               pretime = intHour-1
            pro["Hours"] = str(pretime)
            try:
                prefile,prHycom = Hycomname(pro["date"],pro["Hours"])
    
                
    #        except OSError as err: 
            except:
            #goes back 7 days to look for a file if ones is not found 
                    tdelta=1
                    pretime = intHour+22
                    if pretime % 3 != 0:
                       pretime = intHour+23
                       
                    while (tdelta<=7):
                        try:
                            
                            pro["revHours"] = str(pretime)
                            predate = cro - datetime.timedelta(days=tdelta)
                            pro["revdate"]=predate.strftime('%Y%m%d')
                            prepath_template = r'S:\SOCCER\Maps\HYCOM3Hr\{predate}'
    #                        prepath_template = r'S:\SOCCER\WorkingDirectory-QZH\{revdate}'
                            prepath = prepath_template.format(**pro)
                            os.chdir(prepath)
                            
                            prefile,prHycom = Hycomname(pro["revdate"],pro["revHours"])
    
                            break
    #                    except OSError as err:
                        except:
                            
                            tdelta = tdelta+1
                            pretime = pretime+24
                            
                            
                            
    
                                    
                                    
            if tdelta >=1:              
                postime=(pretime-(tdelta*24))+3
            else:
                postime=pretime+3 #hour after target time
            pro["Hours"] = str(postime)
            try:
                os.chdir(path)
                postfile,poHycom=Hycomname(pro["date"],pro["Hours"])
    
    #        except OSError as err:
            except :
                #goes back 7 days in time if file was not found earlier
                 tdelta=1
                 postime=postime+24
                 while (tdelta<=7):
                     try:
                         
                         pro["revHours"] = str(postime)
                         predate = cro - datetime.timedelta(days=tdelta)
                         pro["revdate"]=predate.strftime('%Y%m%d')
                         prepath_template = r'S:\SOCCER\Maps\HYCOM3Hr\{revdate}'
    #                     prepath_template = r'S:\SOCCER\WorkingDirectory-QZH\{revdate}'
                         prepath = prepath_template.format(**pro)
                         os.chdir(prepath)
                         
                         postfile,poHycom=Hycomname(pro["revdate"],pro["revHours"])
    
                         break
    #                 except OSError as err:
                     except :
                         tdelta = tdelta+1
                         postime = postime+24
            
                
        #Interpolation of prefiles and postfiles
            
        if intHour % 3 != 0:        
         
            
            #Reads in the Lon and Lat variables from Hycom
    
            lonHy,latHy = lats_lons(prHycom)
            
            #Reads in waters before target time
    
            pruHy,prvHy = HYvelocities(prHycom)
                
            #Reads in waters after target time
     
            pouHy,povHy = HYvelocities(poHycom)
               
            preU = Interpwaters(lonHF, latHF,lonHy,latHy, pruHy)
            preV = Interpwaters(lonHF, latHF,lonHy,latHy, prvHy)
            posU = Interpwaters(lonHF, latHF,lonHy,latHy, pouHy)
            posV = Interpwaters(lonHF, latHF,lonHy,latHy, povHy)
             
                        
            U=TimeInterp(preU,posU,pretime,postime,intHour)
            V=TimeInterp(preV,posV,pretime,postime,intHour)
                
                
                
    
    #reshapes large array into proper 2d format
        U=np.asarray(U)
        U=U.reshape(len(latHF),len(lonHF))
        
      
        U2[i,:,:]=U   
    
        V=np.asarray(V)
        V=V.reshape(len(latHF),len(lonHF))
        
      
        V2[i,:,:]=V  
        
        print('%d percent complete' % percentage)    
    t1 = time.clock()
    
    total = t1-t0
    
    NCwriteHytoHF(fileName,U2,V2)    
    
    
    
    print('Hycom waters added to HF file')
    
    return total



def Track7dayinterp(fileName):
#if q==1 :  
    t0 = time.clock()
    
    tdelta=0#step of days to go back
    step=0#hours to go back
    pro ={}
    pro['daystep']=str(tdelta+1)# turning the day step into a string
    
    #templates for creating variables for the nc file.
    day_template_u= 'HYCOM_3Hr_u_d{daystep}'
    day_template_v= 'HYCOM_3Hr_v_d{daystep}'
    longname_temp_u= 'Hycom_eastward_sea_water_velocity_{daystep}day_back'
    longname_temp_v= 'Hycom_northward_sea_water_velocity_{daystep}day_back'
    
    nc = Dataset(fileName,'a')#reads in dataset
    
    
    
    
    # Reading in Longitude, Latitude, and Time
    
    #timeHF = nc.variables['time'][:] 
    
    lonHF,latHF = lats_lons(nc)#gets lon and lats of inputed nc file
    timeHF,dt = get_Time(nc)#gets time and initial date of inputed nc file

    
   
    
    #while loop to go back 7 days
    while (tdelta<7):
        #creating final arrays here
        U2=[]
        V2=[]
        percentaged=(tdelta/7)*100
        for i in range(len(lonHF)):
            U=[]
            V=[]    
            percentage=((i+1)/len(lonHF))*100

    # editing time depending on file
            if dt==datetime.datetime(1979,1,1,0,0): 
                dt2=datetime.datetime.min
                
                offsetD=datetime.timedelta(days=1)
                dateDF=dt2+datetime.timedelta(days=timeHF[0,i])#adds days to date
                dateDF=dateDF-offsetD
                cro=dateDF-relativedelta(years=1)
    
            elif dt==datetime.datetime(1900,1,1,0,0):
                cro=dt+datetime.timedelta(days=timeHF[i])
                
            elif dt==datetime.datetime(1980,1,1,0,0):
                cro=dt+datetime.timedelta(minutes=timeHF[i])#adds minutes to date
            
            #converts date and Hour into strings to open files
            pro["date"]=cro.strftime('%Y%m%d')
            pro["Hours"] = cro.strftime('%H')
            
            
            
                
            intHour = int(pro["Hours"])#makes Hour into an Integer
            
            
            revtime = intHour+step#Hours, will go up to 168
            pro["revHours"] = str(revtime)
            predate = cro - datetime.timedelta(days=tdelta)#goes back tdelta days
            pro["revdate"]=predate.strftime('%Y%m%d')
            prepath_template = r'S:\SOCCER\Maps\HYCOM3Hr\{revdate}'
            prepath = prepath_template.format(**pro)
            #File name making for files that have hours with multiple of 3
            if intHour % 3 == 0 or intHour == 0:
               #if hours is multiple of three is runs code below       
                try:
                    os.chdir(prepath)
                    revfile,Hycom = Hycomname(pro["revdate"],pro["revHours"])
                    lonHy,latHy = lats_lons(Hycom)
                    
            
                    uHy,vHy = HYvelocities(Hycom)
            
            #Loop for interpolating u values in Hycom    
                    u = Interpwaters(lonHF[i],latHF[i], lonHy, latHy, uHy)
                    v = Interpwaters(lonHF[i],latHF[i], lonHy, latHy, vHy)  
                    U=u[0]
                    V=v[0]
                    #print("%f" %U)
        #        except OSError as err: 
                except :
                   #if no file it adds nans 
                    U=math.nan
                    V=math.nan
                    #print("%f" %U)
           
                

                    
                        
              #File name making for hours not multiple of 3         
        
            else :   
        #        
                
                pretime = intHour+(step-2)
                if pretime % 3 != 0:
                    pretime = intHour+(step-1)
                
                postime = pretime+3
                
                pro["prHours"] = str(pretime)
                pro["poHours"] = str(postime)
                try:
                    os.chdir(prepath)
                    prefile,prHycom=Hycomname(pro["date"],pro["prHours"])
    
                    postfile,poHycom=Hycomname(pro["date"],pro["poHours"])
                    
                    lonHy,latHy = lats_lons(prHycom)
                    
                    #Reads in waters before target time
            
                    pruHy,prvHy = HYvelocities(prHycom)
                        
                    #Reads in waters after target time
             
                    pouHy,povHy = HYvelocities(poHycom)
                       
                    preU = Interpwaters(lonHF[i], latHF[i],lonHy,latHy, pruHy)
                    preV = Interpwaters(lonHF[i], latHF[i],lonHy,latHy, prvHy)
                    posU = Interpwaters(lonHF[i], latHF[i],lonHy,latHy, pouHy)
                    posV = Interpwaters(lonHF[i], latHF[i],lonHy,latHy, povHy)
                     
                            
                    u=TimeInterp(preU,posU,pretime,postime,intHour)
                    v=TimeInterp(preV,posV,pretime,postime,intHour)
                    U=u[0]
                    V=v[0]
                    #print("not m3 %f" %U)
        #        except OSError as err:
                except :
                    
                    U=math.nan
                    V=math.nan
                    #print("not m3 %f" %U)

        
        

            
          #appends calculated U and V into initialized array
            U2.append(U)   
        

            
          
            V2.append(V) 
            print('%d percent complete of current day' % percentage)
         
            
        #creating variables and adding data.    
        dayU=day_template_u.format(**pro)
        dayV=day_template_v.format(**pro)
        longnameu=longname_temp_u.format(**pro)
        longnamev=longname_temp_v.format(**pro)
        today= datetime.datetime.today() 
        date=today.strftime('%Y%m%d')
        
        
        nc.createVariable(dayU,'float',('r','c',))
        nc.createVariable(dayV,'float',('r','c',))    
        
        nc.variables[dayU][:,0]=U2[:]
        nc.variables[dayU].units='m/s'
        nc.variables[dayU].long_name=longnameu
        nc.variables[dayU].date_added=date
        nc.variables[dayU].added_by='qzh14'    
        
        nc.variables[dayV][:,0]=V2[:]
        nc.variables[dayV].units='m/s'
        nc.variables[dayV].long_name=longnamev
        nc.variables[dayV].date_added=date
        nc.variables[dayV].added_by='qzh14' 
#             
        tdelta=tdelta+1
        pro['daystep']=str(tdelta+1)
        print('%d percent complete' % percentaged)      
            
        
    t1 = time.clock()
    total = t1-t0
        
        
            
        
        
        
    print('Hycom waters added to track file')
        
    return total




def Hf7dayinterp(fileName): 
    t0 = time.clock()
    tdelta=0
    step=0
    pro ={}
    pro['daystep']=str(tdelta+1)
    day_template_u= 'HYCOM_3Hr_u_d{daystep}'
    day_template_v= 'HYCOM_3Hr_v_d{daystep}'
    longname_temp_u= 'Hycom_eastward_sea_water_velocity_{daystep}day_back'
    longname_temp_v= 'Hycom_northward_sea_water_velocity_{daystep}day_back'
    
    nc = Dataset(fileName,'a')
    
    
    
    
    # Reading in Longitude, Latitude, and Time
    
    #timeHF = nc.variables['time'][:] 
    
    lonHF,latHF = lats_lons(nc)
    timeHF,dt = get_Time(nc)
    ## Converting original date into a datetime variable
    #tunit = nc.variables['time'].units
    #tunit = tunit.replace('T', ' ')
    #tunit = tunit.replace('Z', '')            
    #tunit = tunit.replace('hours since ', '')
    #dt = parser.parse(tunit)
    
    #creates an empty array to store the Hycom data

    
    
    
    while (tdelta<7):
        U2=np.empty([len(timeHF),len(latHF),len(lonHF)], dtype=float)
        V2=np.empty([len(timeHF),len(latHF),len(lonHF)], dtype=float)
        percentage=(tdelta/7)*100
        for i in range(len(timeHF)):
            U=[]
            V=[]    
        
            #Adds hours to the Original date
            #8199
            #6430
            cro=dt+datetime.timedelta(hours=timeHF[i])
        
            
            pro["date"]=cro.strftime('%Y%m%d')
            pro["Hours"] = cro.strftime('%H')
            
            
            
                
            intHour = int(pro["Hours"])
            
            
            
            revtime = intHour+step
            pro["revHours"] = str(revtime)
            predate = cro - datetime.timedelta(days=tdelta)
            pro["revdate"]=predate.strftime('%Y%m%d')
            prepath_template = r'S:\SOCCER\Maps\HYCOM3Hr\{revdate}'
            prepath = prepath_template.format(**pro)
            #File name making for files that have hours with multiple of 3
            if intHour % 3 == 0 or intHour == 0:
                      
                try:
                    os.chdir(prepath)
                    revfile,Hycom = Hycomname(pro["revdate"],pro["revHours"])
                    lonHy,latHy = lats_lons(Hycom)
                    
            
                    uHy,vHy = HYvelocities(Hycom)
            
            #Loop for interpolating u values in Hycom    
                    U = Interpwaters(lonHF,latHF, lonHy, latHy, uHy)
                    V = Interpwaters(lonHF,latHF, lonHy, latHy, vHy)        
        #        except OSError as err: 
                except :
                    U=np.empty([len(latHF),len(lonHF)], dtype=float)
                    V=np.empty([len(latHF),len(lonHF)], dtype=float)
                    
                    U.fill(math.nan)
                    V.fill(math.nan)
                    
                #Reads in Hycom Values
                

                    
                        
              #File name making for hours not multiple of 3         
        
            else :   
        #        os.chdir(r'S:\SOCCER\WorkingDirectory-QZH\2014927')
                
                pretime = intHour+(step-2)
                if pretime % 3 != 0:
                    pretime = intHour+(step-1)
                
                postime = pretime+3
                
                pro["prHours"] = str(pretime)
                pro["poHours"] = str(postime)
                try:
                    os.chdir(prepath)
                    prefile,prHycom=Hycomname(pro["date"],pro["prHours"])
    
                    postfile,poHycom=Hycomname(pro["date"],pro["poHours"])
                    
                    lonHy,latHy = lats_lons(prHycom)
                    
                    #Reads in waters before target time
            
                    pruHy,prvHy = HYvelocities(prHycom)
                        
                    #Reads in waters after target time
             
                    pouHy,povHy = HYvelocities(poHycom)
                       
                    preU = Interpwaters(lonHF, latHF,lonHy,latHy, pruHy)
                    preV = Interpwaters(lonHF, latHF,lonHy,latHy, prvHy)
                    posU = Interpwaters(lonHF, latHF,lonHy,latHy, pouHy)
                    posV = Interpwaters(lonHF, latHF,lonHy,latHy, povHy)
                     
                            
                    U=TimeInterp(preU,posU,pretime,postime,intHour)
                    V=TimeInterp(preV,posV,pretime,postime,intHour)        
        #        except OSError as err:
                except :
                    U=np.empty([len(latHF),len(lonHF)], dtype=float)
                    V=np.empty([len(latHF),len(lonHF)], dtype=float)
                    
                    U.fill(math.nan)
                    V.fill(math.nan)
        

        
        #reshapes large array into proper 2d format
            U=np.asarray(U)
            U=U.reshape(len(latHF),len(lonHF))
            
          
            U2[i,:,:]=U   
        
            V=np.asarray(V)
            V=V.reshape(len(latHF),len(lonHF))
            
          
            V2[i,:,:]=V   
            
        dayU=day_template_u.format(**pro)
        dayV=day_template_v.format(**pro)
        longnameu=longname_temp_u.format(**pro)
        longnamev=longname_temp_v.format(**pro)
        today= datetime.datetime.today() 
        date=today.strftime('%Y%m%d')
        
        
        nc.createVariable(dayU,'float',('time','lat','lon',))
        nc.createVariable(dayV,'float',('time','lat','lon',))    
        
        nc.variables[dayU][:]=U2[:]
        nc.variables[dayU].units='m/s'
        nc.variables[dayU].long_name=longnameu
        nc.variables[dayU].date_added=date
        nc.variables[dayU].added_by='qzh14'    
        
        nc.variables[dayV][:]=V2[:]
        nc.variables[dayV].units='m/s'
        nc.variables[dayV].long_name=longnamev
        nc.variables[dayV].date_added=date
        nc.variables[dayV].added_by='qzh14' 
             
        tdelta=tdelta+1
        pro['daystep']=str(tdelta+1)
        print('%d percent complete' % percentage)      
            
            
    t1 = time.clock()
    total = t1-t0
            
            
                
            
            
            
    print('Hycom waters added to HF file')
        
    return total



   
    
     
   











