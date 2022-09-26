import pandas as pd 



def get_chemopy(data):

    from PyBioMed import Pymolecule

    mols_text =[]

    if '.smi' in data:
        mols_text =[]
        with open(data, 'r') as f:
            lines = f.readlines()
            for line in lines:
                line = line.strip('\n')
                mols_text.append(line)
    elif 'csv' in data:
        df_data = pd.read_csv(data,header=None)
        df_data.columns = df_data.iloc[0,:].tolist()
        col_smiles= ['smiles','SMILES','Smiles']
        smilescol = [m for m in col_smiles if m in df_data.iloc[0,:].tolist()][0]
        mols_text =df_data[smilescol].tolist()
    else:
        pass

    chem = {} 
    error = 0
    for i,smi in enumerate(mols_text[1:]):
        mol = Pymolecule.PyMolecule()
        mol.ReadMolFromSmile(smi)
        
        try:
            alldes =mol.GetAllDescriptor()
        except Exception as e:
            alldes= {'smi':'NAN'}
            print 'The {} calculation of descriptors failed'.format(i)
            print e
            error +=1
        d=[]
        for value in alldes.values():
            d.append(value)
        chem[str(i+1)] = d

        
        descriptors = list(alldes.keys())
        
    
    df_chemdict = pd.DataFrame(chem,index=descriptors).T

    s = pd.Series(list(df_chemdict.index))
    df_index = pd.to_numeric(s, downcast='signed')
    df_chemdict['index'] =df_index.values
    df_chemdict = df_chemdict.sort_values(by=['index'])
    df_chemdict.set_index('index')
    df_chemdict.reset_index(inplace=True,drop=True)

    return df_chemdict 

def get_shet_ds(data):
 
    import pychem 
    from rdkit import Chem
    from pychem import estate
    import pandas as pd


    if '.smi' in data:
        mols_text =[]
        with open(data, 'r') as f:
            lines = f.readlines()
            for line in lines:
                line = line.strip('\n')
                mols_text.append(line)
    elif 'csv' in data:
        df_data = pd.read_csv(data,header=None)
        df_data.columns = df_data.iloc[0,:].tolist()
        col_smiles= ['smiles','SMILES','Smiles']
        smilescol = [m for m in col_smiles if m in df_data.iloc[0,:].tolist()][0]
        mols_text =df_data[smilescol].tolist()
    else:
        pass

    shet =[]
    ds = []
    error = 0
    for i,smi in enumerate(mols_text[1:]):
        mol = Chem.MolFromSmiles(smi)
        try:
            shet.append(estate.CalculateHeteroEState(mol))
        except Exception as e:
            shet.append('NAN')
            print 'The {} calculation of Shet descriptor failed'.format(i)
            print e
            error +=1
        try:
            ds.append(estate.CalculateDiffMaxMinEState(mol))
        except Exception as e:
            ds.append('NAN')
            print 'The {} calculation of DS descriptor '.format(i)
            print e
            error +=1
    
    # chemdict = dict([(k, pd.Series(v)) for k, v in chem.items()])
    
    df_shet_ds = pd.concat([pd.Series(shet),pd.Series(ds)],axis=1)


    return df_shet_ds

def all_chemopy(filename):
    df1 = get_chemopy(filename)
    df2 = get_shet_ds((filename))
    df2.rename({0:'Shet',1:'DS'},inplace=True,axis=1)
    chemdf = pd.concat([df2,df1],axis=1)
    chemdf.reset_index(drop=True,inplace=True)

    chemdf.drop(columns=['index'],axis=1,inplace=True)
    df_smi =pd.read_csv(filename)
    df = pd.concat([chemdf,df_smi],axis=1)
    for col in list (df.columns): 
        df =df [~df [col].isin ([ 'null', 'NULL','Null' ])]
    df.reset_index(drop=True,inplace=True)
    return df

    

