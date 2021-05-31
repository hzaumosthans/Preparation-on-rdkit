#载入资源包
import sys
import os
from fpdf import FPDF
from PIL import Image
from rdkit import Chem
from rdkit import RDConfig
from rdkit.Chem import AllChem
from rdkit.Chem.AllChem import EmbedMolecule
from rdkit.Chem import rdchem
from rdkit.Chem import Descriptors
from rdkit.Chem import SanitizeMol
from rdkit.Chem import rdmolops
from rdkit.Chem import Draw
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.Draw.MolDrawing import MolDrawing, DrawingOptions

#构造分子结构检查函数
def check(file_name):
    #判断文件类型
    pdlx(file_name)
    #读入sdf文件
    suppl = Chem.SDMolSupplier(file_name)
    mols = [x for x in suppl]
    DrawingOptions.bondLineWidth=1.8

    #空分子
    k=[]
    kzz=0
    for mol in mols:
        if mol is None:
            k.append(kzz)
        kzz=kzz+1
        
    #生成去除空分子之后的文件  Xand 生成除空分子图X
    czs=[]
    for mol in mols:
        if mol is None: continue
        czs.append(Chem.MolToSmiles(mol))
    mols=[]
    for cz in czs:
        mols.append(Chem.MolFromSmiles(cz))   
    outname_1 = file_name.split(".sdf")[0] + "_delblank.sdf"
    writer = Chem.SDWriter(outname_1)
    writer.SetProps(['LOGP', 'MW'])
    for i, mol in enumerate(mols):
        if mol is None: continue
        mw = Descriptors.ExactMolWt(mol)
        logp = Descriptors.MolLogP(mol)
        mol.SetProp('MW', '%.2f' %(mw))
        mol.SetProp('LOGP', '%.2f' %(logp))
        mol.SetProp('_Name', 'No_%s' %(i))
        writer.write(mol)
    writer.close()
    suppl = Chem.SDMolSupplier(outname_1)
    mols = [x for x in suppl]
    DrawingOptions.bondLineWidth=1.8
    #img=Draw.MolsToGridImage(mols,molsPerRow=4,subImgSize=(200,200),legends=[x.GetProp("_Name") for x in mols])   
    #k_outname = file_name.split(".sdf")[0] + "_delblank.png"
    #img.save(k_outname)
    
    #脱盐预处理：保留最大片段（去除小离子或配体）并输出去除的片段
    zd=[]
    zdzz=0
    for mol in mols:
        if mol is None: continue
        xx=Chem.MolToSmiles(mol)
        if "." in list(xx):
            largest_Fragment = rdMolStandardize.LargestFragmentChooser()
            mols[zdzz] = largest_Fragment.choose(mols[zdzz])
            dp=Chem.MolToSmiles(mols[zdzz])
            ap=xx.split('.')
            qp=""
            for pd in ap:
                if pd!=dp:
                    qp=qp+";"+pd
            zd.append((zdzz,qp))
        zdzz=zdzz+1
        
        
    #重复结构的检查与删除：去重和分子数   
    #结构有效性的检查：去除【包含7个以上[B] or 任何一个[Sc], [Ti], [V], [Cr], [Mn], [Fe], [Co], [Ni], [Cu], [Ga], [Y], [Zr], [Nb], [Mo], [Tc], [Ru], [Rh], [Pd], [Cd], [In], [Sn], [La], [Hf], [Ta], [W],[Re], [Os], [Ir], [Pt], [Au], [Hg], [Tl], [Pb], [Bi], [Po], [Ac], [Ce], [Pr], [Nd], [Pm], [Sm], [Eu], [Gd], [Tb], [Dy], [Ho], [Er], [Tm], [Yb], [Lu], [T], [Pa], [U], [Np], [Pu], [Am], [Cm], [Bk], [Cf], [Es], [Fm], [Md], [No], [Lr], [Ge], [Sb].等元素的有机金属化合物】
    qc=[]
    cfzz=[]
    yjzz=[]
    zz=0 
    yj=['[Sc]','[Ti]','[V]','[Cr]','[Mn]','[Fe]','[Co]','[Ni]','[Cu]','[Ga]','[Y]','[Zr]','[Nb]','[Mo]','[Tc]','[Ru]','[Rh]','[Pd]','[Cd]','[In]','[Sn]','[La]','[Hf]','[Ta]','[W]','[Re]','[Os]','[Ir]','[Pt]','[Au]','[Hg]','[Tl]','[Pb]','[Bi]','[Po]','[Ac]','[Ce]','[Pr]','[Nd]','[Pm]','[Sm]','[Eu]','[Gd]','[Tb]','[Dy]','[Ho]','[Er]','[Tm]','[Yb]','[Lu]','[T]','[Pa]','[U]','[Np]','[Pu]','[Am]','[Cm]','[Bk]','[Cf]','[Es]','[Fm]','[Md]','[No]','[Lr]','[Ge]','[Sb]']
   # yj=['Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Ga','Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Cd','In','Sn','La','Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','Ac','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','T','Pa','U','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr','Ge','Sb']
    for mol in mols:
        if mol is None: continue
        qc.append(Chem.MolToSmiles(mol))
    fzsl_0=len(qc)+len(k)
    for x in qc:
        zzc=zz+1
        for i in qc[zz+1:]:
            if x==i:
                cfzz.append((zz,zzc))
                qc[zz]='wu'
                break
            zzc=zzc+1
        l=list(x)
        if l.count('B') >= 7:
            yjzz.append((zz,'>=7B'))
            qc[zz]='wu' 
        for i in yj:
            if i in x:
                yjzz.append((zz,i))
                qc[zz]='wu'
        zz=zz+1
    fzsl=len(qc)+len(k)-len(cfzz)-len(yjzz)-len(k)
    mols=[]
    molsy=[]
    for li in qc:
        if li=='wu':continue
        mols.append(Chem.MolFromSmiles(li))
        molsy.append(Chem.MolFromSmiles(li))
            

    #结构有效性的检查：价键错误
    jj=[];jjzz=0;jjgs=0
    for mol in mols:
        if mol is None: continue
        problems = Chem.DetectChemistryProblems(mol)
        c=len(problems)
        for i in range(c):
            if problems[i].GetType()==AtomValenceException:
                jj.append((jjzz,problems[i].Message()))
                jjgs=jjgs+1
                break
        jjzz=jjzz+1    
    
    #结构表征标准性的检查:苯环凯库勒式错用
    bh=[];bhzz=0;bhgs=0
    for mol in mols:
        if mol is None: continue
        problems = Chem.DetectChemistryProblems(mol)
        c=len(problems)
        for i in range(c):
            if problems[i].GetType()==KekulizeException:
                bh.append((bhzz,problems[i].Message()))
                bhgs=bhgs+1
                break
        bhzz=bhzz+1    

    #结构表征标准性的检查:中和电荷
    dhzz=0
    dh=[]
    bj=[]
    for mol in mols:
        if mol is None: continue
        molc,c=neutralize_atoms(mol)
        if c==1:
            largest_Fragment = rdMolStandardize.LargestFragmentChooser()
            molsy[dhzz] = largest_Fragment.choose(molsy[dhzz])
            bj.append((Chem.MolToSmiles(molsy[dhzz]),Chem.MolToSmiles(molc)))
            mols[dhzz]=molc
            dh.append(dhzz)
        dhzz=dhzz+1
        


        
    #分子量分布范围
    fi=0
    mxfl=0
    mifl=100
    for mol in mols:
        fi=fi+1
        if mol is None: continue
        fl=Descriptors.ExactMolWt(mol)
        if fl >mxfl:
            mxfl=fl
        if fl <mifl:
            mifl=fl

    #异构体和互变结构的枚举：筛选未标注手性原子
    sx=[]
    sxzz=0
    sxgs=0
    for mol in mols:
        if mol is None: continue
        xx=Chem.FindMolChiralCenters(mol,force=True,includeUnassigned=True)
        for i in range(len(xx)):
            c=list(xx[i])
            if "?" in c:
                sx.append((sxzz,c[0]))
        for i in range(len(xx)):
            c=list(xx[i])
            if "?" in c:
                sxgs=sxgs+1
                break
        sxzz=sxzz+1

    #异构体和互变结构的枚举：未指定的原子立体中心的数量（即rdkit计算出二molfile未指定的立体中心数目）
    lt=[]
    ltzz=0
    for mol in mols:
        if mol is None: continue
        wlgs=Chem.rdMolDescriptors.CalcNumUnspecifiedAtomStereoCenters(mol)
        if wlgs!=0:
            lt.append((ltzz,wlgs))
        ltzz=ltzz+1

   
    #异构体和互变结构的枚举：筛选未标注顺反键
    sf=[]
    sfzz=0
    sfgs=0
    for mol in mols:
        Chem.rdmolops.FindPotentialStereoBonds(mol)
        if mol is None: continue
        for b in mol.GetBonds():
            ba=str(b.GetBeginAtomIdx())+'~'+str(b.GetEndAtomIdx())
            bt=str(b.GetBondType())
            bs=str(b.GetStereo())
            if bt=='DOUBLE' and bs=='STEREOANY':
                sf.append((sfzz,ba))
        for b in mol.GetBonds():
            ba=str(b.GetBeginAtomIdx())+'~'+str(b.GetEndAtomIdx())
            bt=str(b.GetBondType())
            bs=str(b.GetStereo())
            if bt=='DOUBLE' and bs=='STEREOANY':
                sfgs=sfgs+1
                break
        sfzz=sfzz+1
           
    #生成三维结构 AND  转3D、原3D分类
    sdzz=0
    z_sd=[]
    zsl=0
    for mol in mols:
        if mol is None: continue
        xx =Chem.MolToMolBlock(mol)
        xxn=xx.split('\n')
        D=xxn[1].split(' ')[15]
        if D=="2D":
            mol=Chem.AddHs(mol)
            AllChem.EmbedMolecule(mol,randomSeed=0xf00d)
            mols[sdzz] = Chem.RemoveHs(mol)
            z_sd.append((sdzz,'y'))
            zsl=zsl+1
        if D=="3D":
            z_sd.append((sdzz,'n'))
        sdzz=sdzz+1    
        
    #载入处理后的分子结构信息
    outname_2 = file_name.split(".sdf")[0] + "_out.sdf"
    writer = Chem.SDWriter(outname_2)
    writer.SetProps(['LOGP', 'MW'])
    for i, mol in enumerate(mols):
        if mol is None: continue
        mw = Descriptors.ExactMolWt(mol)
        logp = Descriptors.MolLogP(mol)
        mol.SetProp('MW', '%.2f' %(mw))
        mol.SetProp('LOGP', '%.2f' %(logp))
        mol.SetProp('_Name', 'No_%s' %(i))
        writer.write(mol)
    writer.close()

    #以PDF格式输出检查结果
    pdf = FPDF(orientation='P', unit='mm', format='A4')
    pdf.add_page()
    pdf.set_font("Arial", size=20)
    pdf.cell(200, 10, txt="Molecules Check Report", ln=1, align="C")
    pdf.write(8, u'\n')

    pdf.set_font("Arial", size=18)
    pdf.write(8, u'The original sdf:%s\n\n'%file_name)
    pdf.set_font("Arial", size=15)
    pdf.write(8, u'The original molecules number:%d\n'%fzsl_0)
    pdf.write(8, u'The number of molecules after deletion:%d\n'%fzsl)
    pdf.write(8, u'The blank molecules number:%d\t'%len(k))
    pdf.write(8, u'({:.2%})\n'.format(len(k)/fzsl_0))
    for i in k:
        pdf.write(8, u'[%d]\t'%i)
    pdf.write(8, u'\n\n')   
    
    pdf.set_font("Arial", size=18)
    pdf.write(8, u'The del_blank sdf:%s\n\n'%outname_1)
    pdf.set_font("Arial", size=15)
    pdf.write(8, u'largest_Fragment[Molecular,Deleted_Fragment]:%d\t'%len(zd))
    pdf.write(8, u'({:.2%})\n'.format(len(zd)/fzsl))
    for i in zd:
        pdf.write(8, u'[%d %s]\n'%(i[0],i[1]))
    
    pdf.write(8, u'The repeated molecules number[deletion,repeation]:%d\t'%len(cfzz))
    pdf.write(8, u'({:.2%})\n'.format(len(cfzz)/(fzsl_0 - 1)))
    for i in cfzz:
        pdf.write(8, u'[%d==%s]\t'%(i[0],i[1]))
    if len(cfzz)!=0:
        pdf.write(8, u'\n')
    pdf.write(8, u'The organometallic molecules number[deletion,atom]:%d\t'%len(yjzz))
    pdf.write(8, u'({:.2%})\n'.format(len(yjzz)/(fzsl_0 - 1)))
    for i in yjzz:
        pdf.write(8, u'[%d,%s]\t'%(i[0],i[1]))
    if len(yjzz)!=0:
        pdf.write(8, u'\n')
    pdf.write(8, u'\n\n')   

    
    pdf.set_font("Arial", size=18)
    pdf.write(8, u'The processed sdf:%s\n\n'%outname_2)
    pdf.set_font("Arial", size=15)
    pdf.write(8, u'Molecular weight distribution:%f~%f\n\n'%(mifl,mxfl))
    
    pdf.write(8, u'2DTo3D[Molecular,yes/no]:%d\t'%zsl)
    pdf.write(8, u'({:.2%})\n'.format(zsl/fzsl))
    for i in z_sd:
        pdf.write(8, u'[%d,%s]\t'%(i[0],i[1]))
    pdf.write(8, u'\n\n')

    pdf.write(8, u'Neutralizing_charge[Molecular][Compare]:%d\t'%len(dh))
    pdf.write(8, u'({:.2%})\n'.format(len(dh)/fzsl))
    for i in range(len(dh)):
        pdf.write(8, u'[%d]:\n'%dh[i])
        pdf.write(8, u'%s\n%s\n'%(bj[i][0],bj[i][1]))
    pdf.write(8, u'\n\n')    


    pdf.write(8, u'Not marked Mol_Chiral_Centers[Molecular,Number][Atom]:%d\t'%sxgs)
    pdf.write(8, u'({:.2%})\n'.format(sxgs/fzsl))
    fd=sx[0][0]
    pdf.write(8, u'[%d,%s]:\t'%(lt[0][0],lt[0][1]))
    c=0
    for i in sx:
        if fd!=i[0]:
            c=c+1
            pdf.write(8, u'\n[%d,%s]:\t'%(lt[c][0],lt[c][1]))
            fd=i[0]
        pdf.write(8, u'[%s]\t'%i[1])
    pdf.write(8, u'\n\n')

    pdf.write(8, u'Not marked Mol_Z\E[Molecular,Bond]:%d\t'%sfgs)
    pdf.write(8, u'({:.2%})\n'.format(sfgs/fzsl))
    fd=sf[0][0]
    for i in sf:
        if fd!=i[0]:
            pdf.write(8, u'\n')
            fd=i[0]
        pdf.write(8, u'[%d,%s]\t'%(i[0],i[1]))
    pdf.write(8, u'\n\n')

    pdf.write(8, u'Atom_Valence_Exception[Molecular,Exception]:%d\t'%jjgs)
    pdf.write(8, u'({:.2%})\n'.format(jjgs/fzsl))
    for i in jj:
        pdf.write(8, u'[%d,%s]\n'%(i[0],i[1]))
    pdf.write(8, u'\n\n')

    pdf.write(8, u'Kekulize_Exception[Molecular,Exception]:%d\t'%bhgs)
    pdf.write(8, u'({:.2%})\n'.format(bhgs/fzsl))
    for i in bh:
        pdf.write(8, u'[%d,%s]\n'%(i[0],i[1]))
    pdf.write(8, u'\n\n')
    
    check_name = file_name.split(".sdf")[0] + "_checkout.pdf"
    pdf.output(check_name, 'F')

#判断文件类型函数
def pdlx(file_name):
    pdname = file_name.split(".")[1]
    if pdname=="smi":
        outsdf = file_name.split(".txt")[0] + ".sdf"
        fileN = os.path.join(RDConfig.RDCodeDir,'VLib','NodeLib','test_data',file_name)
        suppl = SmilesSupplyNode(fileN,delim="\t",smilesColumn=2,nameColumn=1,titleLine=1)
        ms = [x for x in suppl]
        mols=[]
        for x in ms:
            mols.append(Chem.MolFromSmiles(x))
        writer = Chem.SDWriter(outsdf)
        writer.SetProps(['LOGP', 'MW'])
        for i, mol in enumerate(mols):
            if mol is None: continue
            mw = Descriptors.ExactMolWt(mol)
            logp = Descriptors.MolLogP(mol)
            mol.SetProp('MW', '%.2f' %(mw))
            mol.SetProp('LOGP', '%.2f' %(logp))
            mol.SetProp('_Name', 'No_%s' %(i))
            writer.write(mol)
        writer.close()
        file_name=outsdf

        
#按原子直接中和带点电子
def neutralize_atoms(mol):
    pattern = Chem.MolFromSmarts("[+1!h0!$([*]~[-1,-2,-3,-4]),-1!$([*]~[+1,+2,+3,+4])]")
    at_matches = mol.GetSubstructMatches(pattern)
    at_matches_list = [y[0] for y in at_matches]
    yn=0
    if len(at_matches_list) > 0:
        for at_idx in at_matches_list:
            atom = mol.GetAtomWithIdx(at_idx)
            chg = atom.GetFormalCharge()
            hcount = atom.GetTotalNumHs()
            atom.SetFormalCharge(0)
            atom.SetNumExplicitHs(hcount - chg)
            atom.UpdatePropertyCache()
        yn=1
    return mol,yn

check("molecules1.sdf")
