from __future__ import division
from ase import Atoms
from ase.io import read
from tdb.utils import DRange, CountList
from tdb.utils import BreakdownFormula, FormatFormula
import numpy as np
import uuid
import os
import re

#Class full of errors / status messages
class Reports:
    '''
    This class is used to keep track of reported information.
    '''
    
    MagDataOutcarMissing = "ERROR: OUTCAR appears to be missing magnetization information: "
    OverwriteTrue =  "STATUS: Overwrite set true: "
    BadOutcar =      "STATUS: OUTCAR has changed: "
    IDFileCorrupt =  "STATUS: ID file corrupt: "
    IDFileMissing =  "STATUS: ID file missing: "
    IDFileValid =    "STATUS: ID file valid: "
    PseudoPotInconsistant = "ERROR: Multiple types of Pseudo Potentials are being used: "
    NoFilmLayers =   "ERROR: No Film Layers were found: "
    NoSlabLayers =   "ERROR: No Slab Layers were found: "
    SlabTypeIDFail = "ERROR: Failed to ID Slab Type: "

class CalcData:
    '''
    This class is used to bundle all the input files together.
    '''
    
    def __init__(self, directory):
        self.incar = GetIncarData(directory+"/INCAR")
        self.poscar = open(directory+"/POSCAR").readlines()
        self.outcar = open(directory+"/OUTCAR").readlines()
        self.potcar = open(directory+"/POTCAR").readlines()
        self.kpts = open(directory+"/KPOINTS").readlines()    

def WriteVASP(output):
    '''
    Writes that the calculator is VASP to output.
    '''
    
    output['Calculator'] = 'VASP'
 
def WriteFunctional(output, calcdata):
    '''
    Writes the order of elements in the POTCAR and determines 
    the type of pseudo potential being used.
    '''
    
    potcardata = calcdata.potcar
    functionals = set()
    potcarorder = " "
    for line in potcardata:
        if "TITEL " in line:
            words = line.split()
            functionals.add(words[2])
            potcarorder+=words[3]+" "
    if len(functionals) == 1:
        for only in functionals:
            output["PseudoPotential"] = only
    else:
        output['PseudoPotential'] = "Multiple/Complex"
    output['POTOrder'] = potcarorder
    output['POSOrder'] = " "+(" ".join(calcdata.poscar[0].split()))+" "
    output['POSPOTAgreement'] = (output['POTOrder'].find(output['POSOrder']) == 0)
    
def GenerateIDFile(output, directory, filename, overwrite = False):
    '''
    Makes an ID file.  If overwrite is set to true,
    old ID files will be overwritten. If overwrite is not set to True, 
    it will attempt to validate the ID and keep the old ID if it is 
    still valid.  The new or old ID will be written to the output.
    The end result will be the ID being written to the output. 
    Returns if the old ID file was valid; if it was not, 
    return false, and why.
    
    FILE DEFINITION
    asdasd32asda2124                     -  ID in hex
    1231232151                           -  Last change of outcar
    '''
    
    over = overwrite
    reason = None
    o = open(directory+"/OUTCAR", "r")
    o.readline()
    o.readline()
    outcartime = o.readline().rstrip('\n')
    try:
        if over:
            reason = Reports.OverwriteTrue
            raise Exception()
        idfile = open(directory+filename,"r")
        reason = Reports.IDFileCorrupt
        __id__ = idfile.readline().rstrip('\n')
        outcaroldtime = idfile.readline().rstrip('\n')
        idfile.close()
        reason = None
        if outcartime == outcaroldtime:
            output['CalcID'] = __id__
            return True, Reports.IDFileValid
        else:
            os.remove(directory+filename)
            reason = Reports.BadOutcar
            raise Exception()
    except:
        over=True
        if reason == None:
            reason = Reports.IDFileMissing
            
    if over:
        __id__ = uuid.uuid4().hex
        idfile = open(directory+filename,"w")
        idfile.write(__id__+"\n")       
        idfile.write(outcartime+"\n")
        output['CalcID'] = __id__
        return False, reason
 
def FileCheck(directory, files):
    '''
    Returns true or false based on if the files exist that are 
    specified.  All files should default to being not checked. 
    Files are opened in an attempt to check them, this tends to 
    be the most robust way to check.
    '''
    
    Missing = " "
    for f in files:
        try:
            open(directory+"/"+f,"r")
        except:
            Missing += f + " "
    return Missing
    
def WriteDir(output, directory):
    '''
    Write the directory of calculation to output.
    '''
    
    output['Directory']=os.path.abspath(directory)

def RunBader(directory, bader_command = "bader"):
    '''
    Runs a bader analysis
    '''
    if FileCheck(directory, ["CHGCAR"]):    
        os.system(bader_command +" "+directory+"/CHGCAR")
        for f in ("/ACF.dat", "/AVF.dat", "/BCF.dat"):
            os.rename("."+f,directory+f)
        return True
    return False

def WriteBader(output, directory, atoms, relevant = None, bader_command = "bader"):
    '''
    Writes the bader analysis information to the output
    and runs the bader analysis if needed.
    '''
    
    if not FileCheck(directory,["CHGCAR"]):
        return False
    if not FileCheck(directory,["ACF.dat"]):
        RunBader(directory, bader_command = bader_command)
    file = open(directory+"/ACF.dat")
    acf = file.readlines()
    bader_count = {}
    bader_data = {}
    for atom in atoms:
        bader_count[atom.symbol] = 0
        bader_data[atom.symbol] = 0
    for index in range(0, len(acf)-2):
        acf_line = acf[index+2].split()
        if len(acf_line) == 7: #Valid
            if not relevant == None:
                if not relevant[index]:
                    continue
            bader_data[atoms[index].symbol] += float(acf_line[4])
            bader_count[atoms[index].symbol] += 1
    for key in bader_data:
        if bader_count[key] > 0:
            output['Bader_'+key] = bader_data[key] / bader_count[key]

def WriteSlabFilm(output, atoms, precision = 0.06, smear = 0.6,
                 slab_identifers = (78, 79),
                 excluded_from_density = (8, 1),
                 formula_prefix = ("Pt","Au"),
                 formula_suffix = ("O", "H")):
    '''
    Write the slab and film layer/type information to the output file,
    returns an array of slab layer tags. After the script runs, 
    the atoms will be tagged to appropriate layers.
    '''
    
    #First generate an array of the density of atoms on each step of the Z axis.
    #This step is influenced by precision and smear.
    density = []
    for step in DRange(0,atoms.get_cell()[2][2],precision):
        #Keep track of how many atoms are in this density step
        dcount=0
        for atom in atoms:
            #Exclude some atoms from the density count.  Hydrogen and oxygen are ignored
            #by default as they typically arn't represented in layers very well.
            if not (atom.number in excluded_from_density):
                if ((atom.position[2]-smear) < step) and ((atom.position[2]+smear) > step):
                    dcount+=1
        density.append(dcount)
    #At this point density is populated with the density of atoms on the Z axis
    #Next is identifying the layer ranges.
    startx = -1
    endx = -1
    density_max = -1

    #Ranges are defined
    layer_ranges = []
    
    #Iterate over all the density values
    for d_index in xrange(0,len(density),1):
        #Move the value to a local variable
        value = density[d_index]
        #if there is no know endx
        if endx == -1:
            #check if the value is higher than 1/6th the last known highest density for this section.
            if value > (density_max/6):
                #set all the relevant values
                startx = d_index
                endx = d_index
                density_max = value
        #if you have found a starting index and the value has increased,  then restart the range
        elif value > density_max:
            startx = d_index
            endx = d_index
            density_max = value
        #if the value is the same,  extend the range
        elif value == density_max:
            endx = d_index
        #if the value ever drops under 1/6th of the highest known density for this section,  then the
        #layer can be defined with known values and the startx/endx can be reset.  Don't reset the density
        #max as there is no good reason to.  All layers should have densities that are at least that similar
        #or the definition of a layer is hard to define.  Also might cause errors as moving down away from
        #the density peak.
        elif value <= density_max/6:
            layer_ranges.append((startx-(smear/precision),endx-startx,endx+(smear/precision)))
            startx = -1
            endx = -1 
    
    #Find which atoms belong in which layers and tag them.
    for index, layer_range in enumerate(layer_ranges):
        atoms_in_layer = [atom for atom in atoms if (atom.position[2]/precision > layer_range[0])
                           and (atom.position[2]/precision < layer_range[2])]
        for atom in atoms_in_layer:
            atom.tag=index+1

    #Define a dictionary of atom tag keys and if they are part of the slab
    slab_layers = {}
    
    #Try to indentify which layers belong in the slab
    for atom in atoms:
        if atom.tag == 0:
            continue
        if atom.number in slab_identifers:
            slab_layers[atom.tag] = True
        else:
            #This part is designed to not overwrite a value if it already exists.  You never want
            #to overwrite a True value.
            try:
                slab_layers[atom.tag]=(slab_layers[atom.tag])
            except:
                slab_layers[atom.tag]=False

    for atom in atoms:
        #if atom is not assigned
        if atom.tag == 0:
            new_layer = 0
            distance = 1000
            #loop over all layer ranges
            for index, layer_range in enumerate(layer_ranges):
                #the atoms excluded from the density are not part of the slab
                if (atom.number in excluded_from_density) and (slab_layers[index+1]):
                    continue
                start = abs((atom.position[2]/smear*precision)-layer_range[0])
                end = abs((atom.position[2]/smear*precision)-layer_range[2])
                if (distance>=start):
                    new_layer = index
                    distance = start
                if (distance>end):
                    new_layer = index
                    distance = start
            atom.tag=new_layer+1
    
    slab_layer_count=0
    film_layer_count=0
    for slab_layer in slab_layers:
        if slab_layers[slab_layer]:
            slab_layer_count+=1
        else:
            film_layer_count+=1

    output['SlabLayers']=slab_layer_count
    output['FilmLayers']=film_layer_count

    if film_layer_count == 0:
        return Reports.NoFilmLayers, None
    if slab_layer_count == 0:
        return Reports.NoSlabLayers, None
        
    slab = Atoms([atom for atom in atoms if (slab_layers[atom.tag])])
    slab.set_cell(atoms.get_cell())
    film = Atoms([atom for atom in atoms if not (slab_layers[atom.tag])])
    film.set_cell(atoms.get_cell())

    slabform = slab.get_chemical_formula()
    filmform = film.get_chemical_formula()

    #Obtain the layers of the slab as individual layers for additional slab analysis
    layers=[]
    known_layers=set()
    for atom in slab:
        known_layers.add(atom.tag)
    for layer in known_layers:
        this_layer=[]
        for atom in atoms:
            if atom.tag == layer:
                this_layer.append(atom)
        layers.append(Atoms(this_layer))

    layers.sort(key=lambda layer: layer[0].position[2])
    
    layer_formulas = []
    for layer in layers:
        layer_formulas.append(layer.get_chemical_formula())
        
    slab_numbers = slab.get_atomic_numbers()
    
    slabtype = "" #Homo,  Hetero,  Skin
    skintype = "" #Only defined if skinned | Plate,  Anneal
    skincomp = "" #Only defined if skinned
    annealcomp = "" #Reduced composition of annealed layer
    slabcomp = layer_formulas[0]
    if (CountList(slab_numbers, slab_numbers[0])==len(slab_numbers)):
        slabtype = "Homo"
    elif CountList(layer_formulas, layer_formulas[0])==len(layer_formulas):
        slabtype = "Hetero"
    elif ((CountList(layer_formulas, layer_formulas[0])==len(layer_formulas)-1) or
          (CountList(layer_formulas, layer_formulas[1])==len(layer_formulas)-1)): #Plated Skin
        slabtype = "Skin"
        skintype = "Plate"
        skincomp = layer_formulas[len(layer_formulas)-1]
    elif ((CountList(layer_formulas, layer_formulas[0])==len(layer_formulas)-2) or
          (CountList(layer_formulas, layer_formulas[1])==len(layer_formulas)-2) or 
          (CountList(layer_formulas, layer_formulas[2])==len(layer_formulas)-2)): #Annealed Skin
        slabtype = "Skin"
        skintype = "Anneal"
        skincomp = layer_formulas[len(layer_formulas)-1]
        annealcomp = layer_formulas[len(layer_formulas)-2]
    else: #Failed
        return Reports.SlabTypeIDFail, None

    output['SystemFormula'] = FormatFormula(atoms.get_chemical_formula(), formula_prefix, formula_suffix, red = False)
    output['SlabFormula'] = FormatFormula(slabform, formula_prefix, formula_suffix, red = False)
    output['SlabType'] = slabtype
    output['SkinType'] = skintype
    output['SkinComp'] = FormatFormula(skincomp, formula_prefix, formula_suffix)
    output['AnnealComp'] = FormatFormula(annealcomp, formula_prefix, formula_suffix)
    output['SlabComp'] = FormatFormula(slabcomp, formula_prefix, formula_suffix)
    output['FilmFormula'] = FormatFormula(filmform, formula_prefix, formula_suffix, red = False)
    output['FilmComp'] = FormatFormula(filmform, formula_prefix, formula_suffix)

    return False, slab_layers

def WriteVDW(output, calcdata):
    '''
    Write the VDW information.
    '''
    
    incardata = calcdata.incar
    AL = calcdata.poscar[5]
    AL = AL.split()
    
    if "MAGMOM" in incardata:
        output['MagMomentIn'] = " ".join(incardata["MAGMOM"])

    VDW=True
    if not ("LUSE_VDW" in incardata):
        VDW=False
    elif not (incardata["LUSE_VDW"][0].lower() == ".true."):
        VDW=False
    output['VDW']=VDW

    if not VDW:
        return
    UAL = {}
    for x, a in enumerate(AL):
        if float(incardata["LDAUU"][x]) > 0:
            UAL[a]=float(incardata["LDAUU"][x])
        else:
            UAL[a]=0

    LAL = {}
    for x, a in enumerate(AL):
        if float(incardata["LDAUL"][x]) > 0:
            LAL[a]=float(incardata["LDAUL"][x])
        else:
            LAL[a]=0

    JAL = {}
    for x, a in enumerate(AL):
        if float(incardata["LDAUJ"][x]) > 0:
            JAL[a]=float(incardata["LDAUJ"][x])
        else:
            JAL[a]=0

    for a in AL:
        output['VDW_U_'+a]=UAL[a]
        output['VDW_L_'+a]=LAL[a]
        output['VDW_J_'+a]=JAL[a]
   
    output['VDW_U_'+"All"]=" ".join(incardata["LDAUU"])
    output['VDW_L_'+"All"]=" ".join(incardata["LDAUL"])
    output['VDW_J_'+"All"]=" ".join(incardata["LDAUJ"])
    
def WriteGGA(output, calcdata):
    '''
    Write the xc functional information.
    '''
    
    output["Functional"]=IncarOutcarParameter(calcdata, "GGA")[0]

def WriteSystemFormula(output, atoms):
    '''
    Write the system's formula
    '''
    
    output['SystemFormula'] = atoms.get_chemical_formula()
    
def WriteEnergy(output, calcdata):
    '''
    Write the energies of the system.
    '''
    
    outcar_data = calcdata.outcar
    toten=0
    noentropy=0
    sigmaen=0
    for line in outcar_data:
        try:
            line_data = line.split()
            if line_data[0] == "free" and line_data[1] == "energy":
                toten = float(line_data[4])
            if line_data[0] == "energy" and line_data[1] == "without":
                noentropy = float(line_data[4])
                sigmaen = float(line_data[7])
        except:
            pass
    
    output['E0'] = toten
    output['NoEntropy0'] = noentropy
    output['NoSigma0'] = sigmaen

def WriteDipoleCorrection(output, calcdata):
    '''
    Write the dipole correction.
    '''
    
    incar_data = calcdata.incar
    if not "IDIPOL" in incar_data:
        output['DipoleCorrection'] = 0
    else:
        output['DipoleCorrection'] = int(incar_data["IDIPOL"][0])
    
def WriteMag(output, calcdata):
    '''
    Write the magnetization information if the calculation
    was told to output per atom magnetization information.
    '''
    
    incar_data = calcdata.incar
    outcar_data = calcdata.outcar
    if "LORBIT" in incar_data:
         if incar_data["LORBIT"][0] == "10":
             pass
         else:
             return
    else:
         return
    startmag=0
    
    for index, line in enumerate(outcar_data):
        line_data = line.split()
        if len(line_data)==0:
           continue
        if line_data[0] == "magnetization":
            startmag = index            
    mags=[]
    startmag+=4

    maxiter=5000
    while(True):
        maxiter-=1
        if maxiter <= 0:
            return Reports.MagDataOutcarMissing
        if outcar_data[startmag].split()[0][0]=="-":
            break
        else:
            mags.append(float(outcar_data[startmag].split()[4]))
            startmag+=1
    
    magmom = 0
    for index, mom in enumerate(mags):
        magmom += mom

    output["MagMoment"]=magmom

def WriteFilmMag(output, calcdata, relevant = None, negligible_moment=0.3):
    '''
    Write the magnetization information if the calculation
    was told to output per atom magnetization information,
    with additional film specific information.
    '''
    
    incar_data = calcdata.incar
    outcar_data = calcdata.outcar
    if "LORBIT" in incar_data:
         if incar_data["LORBIT"][0] == "10":
             pass
         else:
             return
    else:
         return
    startmag=0
    
    for index, line in enumerate(outcar_data):
        line_data = line.split()
        if len(line_data)==0:
           continue
        if line_data[0] == "magnetization":
            startmag = index            
    mags=[]
    startmag+=4

    maxiter=5000
    while(True):
        maxiter-=1
        if maxiter <= 0:
            return Reports.MagDataOutcarMissing
        if outcar_data[startmag].split()[0][0]=="-":
            break
        else:
            mags.append(float(outcar_data[startmag].split()[4]))
            startmag+=1
    
    total_moms = 0
    abs_mom = 0
    polarized_mom = 0
    with_negl_mom = 0
    film_total = 0
    for index, mom in enumerate(mags):
        with_negl_mom += mom
        if not relevant == None:
            if relevant[index] == False:
                continue
        film_total += mom
        if abs(mom) > negligible_moment:
            total_moms+=1
            polarized_mom+=mom
            abs_mom+=abs(mom)

    if total_moms == 0:
        output['FilmFerro']="PM"
    elif polarized_mom <= negligible_moment * total_moms:
        output['FilmFerro']="AFM"
    elif abs_mom == polarized_mom:
        output['FilmFerro']="FM"
    else:
        output['FilmFerro']="FIM"
    output["MagMoment"]=with_negl_mom
    output['FilmMag']=film_total

def WriteSpin(output, calcdata):
    '''
    Write if the calculation was done with spin polarization.
    '''
    
    if IncarOutcarParameter(calcdata, "ISPIN")[0] == "2":
        output["Spin"]=True
    else:
        output["Spin"]=False
    return output["Spin"]

def WriteCell(output, calcdata):
    '''
    Write information concerning the cell of the calculation.
    '''
    
    posdata = calcdata.poscar
    output['CellScale'] = float(posdata[1].split()[0])
    scale = output['CellScale']
    for index, val in enumerate(["X","Y","Z"]):
        for index2, val2 in enumerate(["X","Y","Z"]):
            vect = float(posdata[2+index].split()[index2])
            output['Cell_'+val+val2] = vect * scale
            output['CellUnscaled_'+val+val2] = vect
            
def WriteRelaxation(output, calcdata):
    '''
    Write information concerning the number of electronic steps
    in the last geometric optimization,  and how many geometric
    optimization steps were taken; Also includes maximum bounds.
    '''
    
    incar_data = calcdata.incar
    output['MaxGeometricSteps']=int(IncarOutcarParameter(calcdata, "NSW")[0])
    output['MaxElectronicSteps']=int(IncarOutcarParameter(calcdata, "NELM")[0].rstrip(";"))
    outdata = calcdata.outcar
    iterdata = []
    for line in outdata:
        if "- Iteration" in line:
            iterdata = re.sub('[-()]', ' ', line).split()
    output['ElectronicSteps'] = iterdata[2]
    output['GeometricSteps'] = iterdata[1]
    
def WriteForce(output, calcdata, atoms):
    '''
    Write the max last reported force in the system.
    '''
    
    outdata = calcdata.outcar
    f = []
    last = -1
    for x, line in enumerate(outdata):
        if "TOTAL-FORCE" in line:
            last = x+2
    for x, line in enumerate(outdata):
        if not x == last:
            continue
        if line[1] == "-":
            break
        last += 1
        f.append(line.split()[3:6]) 
    try:
        c = atoms._get_constraints()
        indices_fixed=c[0].index    # the indices of atoms fixed
        for i in indices_fixed:
            f[i] = [0,0,0]
    except:
        pass
    fmax = 0
    for i in f:
        fval = (float(i[0])**2 + float(i[1])**2 + float(i[2])**2)**(1./2.)
        if fval > fmax:
            fmax = fval    
    output['Force'] = fmax
    
def WriteConverged(output, calcdata):
    '''
    Write if calculation reached convergence based on precision.
    '''
    outcar_data = calcdata.outcar
    for line in outcar_data:
        if "reached required accuracy" in line:
            output['Converged'] = True
            return
    output['Converged'] = False
    
def WritePrecision(output, calcdata):
    '''
    Write the precision.
    '''
    incar_data = calcdata.incar
    if "PREC" in incar_data:
        output['Precision']=incar_data["PREC"][0]
    else:
        output['Precision']="Normal"
    output["EnergyDiff"]=float(IncarOutcarParameter(calcdata, "EDIFF")[0])
    output["ForceDiff"]=abs(float(IncarOutcarParameter(calcdata, "EDIFFG")[0]))
    output["EnergyCutoff"] = float(IncarOutcarParameter(calcdata, "ENCUT")[0])

def WriteKPoints(output, calcdata):
    '''
    Write information about the KPoint configuration used
    in the calculation.
    '''
    
    kpointsdata = calcdata.kpts
    if (kpointsdata[2].lower())[0] == "g":
        output['KGammaCentered'] = True
    else:
        output['KGammaCentered'] = False
    kpts = kpointsdata[3].split()
    for index, val in enumerate(["X","Y","Z"]):
        output['KPoints_'+val] = int(kpts[index]) 
    for line in calcdata.outcar:
        if "irreducible k-points:" in line:
            output['KReduced'] = int(line.split()[1])
        
def WriteCustom(output, directory):
    '''
    Used to allow a custom set of keys to be added to a "db.dat"
    file in the directory.  The format is roughly the same as an 
    INCAR file.
    '''
    
    try:
        for d in open("db.dat").readlines():
            D = d.replace("="," = ")
            Dc = D.find("#")
            if not Dc == -1:
                D = D[:Dc]
            D = D.split()
            output[D[0]] = " ".join(D[2:])
    except:
        pass

def IncarOutcarParameter(calcdata, key):
    '''
    Return the value of the INCAR parameter as a split string. If 
    it cannot be found in the INCAR,  attempt to find it in the 
    OUTCAR. If both fail, return None  
    '''
    
    if key in calcdata.incar:
        return calcdata.incar[key]
    else:
        for line in calcdata.outcar:
            try:
                if line.split()[0] == key:
                    return line.split()[2:]
            except IndexError:
                pass
    return None
        
def GetIncarData(filename):
    '''
    Returns an array of the data in the INCAR, 
    data is processed into a dictionary of arrays.
    '''
    
    r = open(filename)
    rD = r.readlines()
    r.close()

    ID = {}

    for i in rD:
        i = i.replace("="," = ")
        parm = i.split()
        val = []
        key = ""
        try:
            for x, ii in enumerate(parm):
                if ii[0] == "#":
                    break
                if x == 0:
                    key = ii
                    continue
                if x == 1:
                    continue
                val.append(ii)
            ID[key] = val
        except:
            pass
    return ID
