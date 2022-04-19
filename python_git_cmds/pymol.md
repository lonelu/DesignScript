# right side comlumn wider
set internal_gui_width, 300

# save fasta file
save ***.fasta, sele

# save pdb file
save ....pdb, sele

# set chain
alter (sele), chain = 'B'

# prepare movie.

python
cmd.mset("1 -%d" % cmd.count_states())
python end

cd E:\DesignData\smallprot_example

mset 1 x320
python
for x in range(0,16):
  for y in range(0, 20):  
    cmd.mview('store', x*20 + y+1, state=x+1, object='Nter')
python end
movie.produce nter.mpg, quality=90


mset 1 x200
python
for x in range(0,10):
  for y in range(0, 20):  
    cmd.mview('store', x*20 + y+1, state=x+1, object='Cter')
python end
movie.produce cter.mpg, quality=90

# save png file
png 1mof.png, width=500, height=1200, dpi=300, ray=1

# change resi number in pymol
alter sele, resi=str(int(resi)-65)

# Change size of sphere
set sphere_scale, 0.35, all

# Group 
group clu20, AAMetalPhiPsi_HIS_cluster_20_*

# Change radius of stick
set_bond stick_radius, 0.1, all

# Align several structure
extra_fit selection that you want to fit (i.e., n. Ca), model to align to, alignment method (i.e., align)

extra_fit n. Ca, 1dmm_16-20-28_H-H-D_a, align

# Crystal structures
symexp sym,1dmm,(1dmm),5

# Select aa
sele resn ALA
https://pymolwiki.org/index.php/Selection_Algebra


# Connect two structure into one structure 
>Build>Create bond
alter sele, chain = 'A'
alter sele, segi = ''
sort
rebuild
