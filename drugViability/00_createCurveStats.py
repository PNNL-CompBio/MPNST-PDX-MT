'''
pull data and run curve fitting code
'''

import pandas as pd
import synapseclient as sc
import os
import subprocess
syn = sc.login()
##this is a pain, need to move to file

filelist = syn.tableQuery("select entityId,individualID,specimenID from syn65473033 where dataType='drugScreen' AND dataSubtype='processed' AND assay='cell viability assay'").asDataFrame()
# raw files are mistakenly marked as 'processed' so check 'dataSubtype' within file

###pull files

singledose = []
multidose = []
for index,row in filelist.iterrows():
  #print(row['id'])
  dfile = pd.read_csv(syn.get(row['entityId']).path)
  if any(dfile['dataSubtype'] == 'processed'):
      dfile['improve_sample_id'] = row['specimenID']
      dfile = dfile.reset_index()
      
      ##get counts of drugs to only keep drugs with >1 dose for curve fitting
      dcounts = dfile.groupby("drugName").count().reset_index()
      if any(dcounts.concentration>1):
          more = dcounts[dcounts.concentration>1]['drugName']
          multidose.append(dfile[dfile.drugName.isin(set(more))])
      
      ## also get single data points
      if any(dcounts.concentration==1):
          sings = dcounts[dcounts.concentration==1]['drugName']
          singledose.append(dfile[dfile.drugName.isin(set(sings))])


####first fit multidose curves...
if len(multidose) > 0:
    fulltab = pd.concat(multidose)
    #print(fulltab)
    #fulltab['DOSE']=fulltab.concentration#+0.0001 # check if should add 0.0001? perhaps cNF data had 0 values
    fulltab = fulltab.rename(columns={"concentration": "DOSE", "percentViability": "GROWTH",
                            "timePoint": "time", "timePointUnit": "time_unit",
                            "drugName": "Drug"})
    fulltab['study']='mpnstPDXMT'
    fulltab['source']='synapse'
    ##mutate the values create new columns
    ncols=['DOSE','GROWTH','study','source','improve_sample_id','Drug','time','time_unit']
    fulltab = fulltab[ncols]
    ##change file headers to DOSE/RESPONSE values needed by other script
    fulltab.to_csv('mpnst_drug_response.tsv',sep='\t')
    
    ##fit curve
    script='https://raw.githubusercontent.com/PNNL-CompBio/coderdata/refs/heads/main/build/utils/fit_curve.py'
    subprocess.run(['wget',script])
    subprocess.run(['python','fit_curve.py','--input','mpnst_drug_response.tsv','--output','mpnstDrugOutput'])
    otab = pd.read_csv('mpnstDrugOutput.0',sep='\t')
else:
    otab = pd.DataFrame()

#####now we can take single drug points and format those
if len(singledose) > 0:
    stab = pd.concat(singledose)
    stab = stab.rename(columns={"drugName": "improve_drug_id", 
                            "timePoint": "time", "timePointUnit": "time_unit"})
    stab['study']='mpnstPDXMT'
    stab['source']='synapse'
    stab['dose_response_value'] = stab.percentViability/100.00
    stab['dose_response_metric'] = 'uM_viability' # not sure, but concentrationUnit='uM' and assay='cell viability assay', platform='3D CellTiter-Glo'

    curve_cols = ['source','improve_sample_id','improve_drug_id','study','time','time_unit',
              'dose_response_metric','dose_response_value']
    stab = stab[curve_cols]
else:
    stab = pd.DataFrame()

newtab =pd.concat([otab,stab])

##rename file
newtab.to_csv('viability_curves.tsv',sep='\t',index=False)

##now add in single-point drug measurements

#store on synapse
syn.store(sc.File('viability_curves.tsv',parentId='syn52369034'))

