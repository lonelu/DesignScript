#right side comlumn wider
set internal_gui_width, 300

#save fasta file
save ***.fasta, sele

#save pdb file
save ....pdb, sele

#set chain
alter (sele), chain = 'B'

#prepare movie.

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

#save png file
png 1mof.png, width=500, height=1200, dpi=300, ray=1


alter sele, resi=str(int(resi)-65)
