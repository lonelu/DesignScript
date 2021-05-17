def Generate_Pdbs(Nres,num_chain,num_to_output,w,R,orientation,helix_phase,delta_z,output_file_name,chain_name, chain_order):
 deg_to_rad= math.pi/180.
# chain  parameters
 ph = 360/num_chain
 phase=[]
 for i in range(num_chain):
     phase.append(i*ph)
 chain_set=['A','B','C','D','E','F','G','H',' I ',' J ',' K','L','M','N','O','P','Q','R','S ',' T','U','V','W','X','Y','Z']
 chain_num=chain_set[0:num_chain]
 R1=2.26
#rise per  residue  d   fixed .  this   constrains   pitch  ( alpha )
 d=1.51
 z1=0.
 z2=0.
 line='ATOM      7  CA  GLY A   2       5.520   2.352   1.361  1.00 20.00     '
 line='ATOM      1  N   GLY A   1       8.594  -0.685 -21.987  1.00 20.00         A'
 last=line[54:-1]
 atom_num=1
 res_num=1
 Res_id=[]
 CA_list=[]
 alpha=math.asin(R*w*deg_to_rad/d)
 for iter in range(num_to_output):
  CA_list.append([])
  Res_id.append([])
  chain = chain_name[iter]
  orient = orientation[iter ]
  w1=100-w  # 102.85  for a 2-layer heptad repeat
  if orient == 1:
     res_num=0
  else :
     res_num=Nres+1
  supercoil_phase=phase[iter]+delta_z[iter]*math.tan(alpha)/(R*deg_to_rad)
  for t in range(Nres+2):      ## need two extra residues to  guide   placement  of 1 st  and  last   residue
    a0=(w*(t-1)+supercoil_phase)*deg_to_rad
    if orient == 1:
     a1=(w1*(t-1)+helix_phase[iter])*deg_to_rad   # set ref point for phase to be along  supercoil   radius
    else:
     a1=(w1*t-w1*Nres-helix_phase[iter])*deg_to_rad
    x=R*math.cos(a0) + R1*math.cos(a0) * cos(a1) - R1*cos(alpha)*sin(a0)*sin(a1)
    y=R*math.sin(a0) + R1*sin(a0)*cos(a1) + R1*cos(alpha)*cos(a0)*sin(a1)
    if w==0:
     z=d*t+delta_z[iter]
    else:
     z= R*w*t*deg_to_rad/math.tan(alpha)-R1*sin(alpha)*sin(a1)+delta_z[iter]
    CA_list[iter].append( (res_num,Vec(x,y,z)) )
    Res_id[iter].append(res_num)
    atom_num=atom_num+1
    if orient  == 1:
        res_num=res_num+1
    else:
        res_num=res_num-1
# convert CA trace to  full   backbone  model by  superimposing  on  ideal   template
# by matching 3  consecutive  CA atoms
# set up  ideal   template
 stub_file=map(string.split,open('ideal.pdb','r ') . readlines ())
 atom=[]
 for line in stub_file:
    atom.append( (Vec(float(line[6]),float( line [7]) , float ( line [8]) )))
 ideal_stub=stub(atom[6],atom[1],atom[11])
#now make full  backbone  pdb
 full_pdb=open(output_file_name,'w')
 atom_num=1
 res_num=0
 for counter in range(num_to_output):
  iter=int(chain_order[counter])
  chain=chain_name[iter]
  CA_chain_u=CA_list[iter]
  CA_chain = sorted(CA_chain_u, key = lambda res: res[0])
  for res in range(1,Nres+1):
    res_num=res_num+1
    actual_stub=stub(CA_chain[res][1],CA_chain[res-1][1],CA_chain[res+1][1])
    transform=actual_stub * ~ideal_stub
# N
    coords=transform*atom[5]
    full_pdb.write('ATOM %6d  N   GLY %s %3d    %8.3f%8.3f%8.3f%s\n'%(atom_num,chain,res_num,coords.x,coords.y,coords.z,last))
    atom_num=atom_num+1
# CA   (use actual CA from trace  rather  than  superimposed  one)
    coords=CA_chain[res][1]
    tcoords=transform*atom[6]
#    print  coords , tcoords ,' CA'
    full_pdb.write('ATOM %6d  CA  GLY %s %3d    %8.3f%8.3f%8.3f%s\n'%(atom_num,chain,res_num,coords.x,coords.y,coords.z,last))
    atom_num=atom_num+1
#  NH
    coords=transform*atom[7]
    full_pdb.write('ATOM %6d  H   GLY %s %3d    %8.3f%8.3f%8.3f%s\n'%(atom_num,chain,res_num,coords.x,coords.y,coords.z,last))
    atom_num=atom_num+1
#  C
    coords=transform*atom[8]
    full_pdb.write('ATOM %6d  C   GLY %s %3d    %8.3f%8.3f%8.3f%s\n'%(atom_num,chain,res_num,coords.x,coords.y,coords.z,last))
    atom_num=atom_num+1
# O
    coords=transform*atom[9]
    full_pdb.write('ATOM %6d  O   GLY %s %3d    %8.3f%8.3f%8.3f%s\n'%(atom_num,chain,res_num,coords.x,coords.y,coords.z,last))
    atom_num=atom_num+1
    start_d=Vec.distance(atom[8],atom[6])
    end_d = Vec.distance(transform*atom[8],transform*atom[6])
 return()