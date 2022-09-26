from rdkit import Chem
from rdkit.Chem import AllChem
import pandas as pd


def getrncaa(file1,file2):
 
    if 'sdf' in file1:
        suppl = Chem.SDMolSupplier(file1)

    else:
         with open(file1, 'r') as f:
            mols_text = f.read()
            suppl = Chem.SmilesMolSupplierFromText(mols_text,delimiter='\t')
    print('The name of file:', file1.split("/")[-1])
    
    mols = [x for x in suppl if x is not None]
    print ('The number of molecules:',len(mols)) 

   #Loading substructure structures

    pattern = AllChem.MolFromSmiles('N[C@@H](C)C(O)=O') 
    pattern1 = AllChem.MolFromSmarts('OC(=O)CN*')
    pattern2 = AllChem.MolFromSmarts('NCC(=O)O*')
    pattern3 = AllChem.MolFromSmarts('N[D3]C(=O)O')
    pattern4 = AllChem.MolFromSmarts('OC(=O)CN(=*)')

    matching_molecules = [m for m in mols if m.HasSubstructMatch(pattern,useChirality = True) if not m.HasSubstructMatch(pattern1) and not m.HasSubstructMatch(pattern2) and m.HasSubstructMatch(pattern3) and not m.HasSubstructMatch(pattern4)]
    print ("matching:",len(matching_molecules)) 
    print("not_matching:",len(mols)- len(matching_molecules))
    
    writer = Chem.SmilesWriter(file2, delimiter=',')
    for mol in matching_molecules:
        writer.write(mol)
    writer.close()

    data = pd.read_csv(file2, sep=',')

    #Second clean-up to remove molecules containing heteroatoms and isotopes
    match = ['Cl.','[2H]','.Cl','13C','3H','[13C]','[NH2+]','[NH3+]','13c','HBr',"O-",'[131I]','[125I]','[123I])','.','15N','12C','14C']  
    data['isContained'] = data['SMILES'].apply(lambda x: 1 if any(s in x for s in match) else 0)
    data = data[data['isContained']== 0 ]
    print('Second clean-up:',data.shape[0])

    # Third clean-up to remove molecules with more than 150 atoms
    data['lenth'] = data['SMILES'].apply(lambda x: len(x))
    data = data[data['lenth']<150]
    data['SMILES'].to_csv(file2,index=False)
    print('Third clean-up:',data.shape[0])
    print('\n')



def edatagenerator(file1,file2):
    if 'csv' in file1:
        data = pd.read_csv(file1)
    elif 'xlsx' in file1:
        data = pd.read_excel(file1)
    else:
        data = file1
    
    data['num_list']=list(range(len(data)))

    pros = {
            'S':[-0.228,1.399,-4.760,0.670,-0.647],
            'I':[-1.239,-0.547,2.131,0.393,0.816],
            'G':[-0.384,1.652,1.330,1.045,2.064],
            'Q':[0.931,-0.179,-3.005,-0.503,-1.853] 
            }

    pos = ['346F1','346F2','346F3','346F4','346F5','348F1','348F2','348F3','348F4','348F5']
    edata =  pd.DataFrame(columns = pos ,index=list(range(len(data))))
    edata['num_list']=list(range(len(data)))

    edata[pos]=pros['S']+pros['I']
    data_ifrs = pd.merge(edata,data,on = 'num_list')

    edata[pos]=pros['S']+pros['Q']
    data_mfrs = pd.merge(edata,data,on = 'num_list')

    edata[pos]=pros['G']+pros['Q']
    data_btars = pd.merge(edata,data,on = 'num_list')

    data_all = pd.concat( [data_ifrs, data_mfrs, data_btars], axis=0)
    data_all = data_all.drop(['num_list'],axis=1)
    data_all.to_csv(file2,index = False)
   