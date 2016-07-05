#!/usr/bin/python
# Read wannier centers
import numpy as np
import sys
import os
from argparse import ArgumentParser
D_field = -1* value
llimit = 1.5
system = "2B2C"
sp1=['Ba',4,10]
sp2=['Ca',4,10]
sp3=['Ti',8,12]
sp4=['O',24,6]
nat = 40
first_layer = 'AO'
threshold = 1.5
parser = ArgumentParser()
parser.add_argument('fname', metavar='FILE', help='file to process')
args = parser.parse_args()
spin_multiplicity=-2
efieldx = -0*D_field * (10.42394257*10.42394257)/(4*3.14159265358979)   #0.0864677496



f = open(args.fname,'r')

i=0;
prim_mat=[]
ev_array=np.zeros((nat,1))
E_field=0
for line in f:
    if line =='  -> xunrefikpar, ispin:            1           1\n':
        wc_array=[]
        i=i+1
        #print i
        count = 0
        for _ in range(32):
            ln=f.next()
            #print ln
            row=[]
            row=ln.split()
            for j in range(5):
                wc=float(row[j])
                wc_array.append(wc)
                count = count +1

        #print count
        #print "count", i
    if 'IIII' in line:
	    tmp_array=[]
	    tmp_array=line.split()
	    E_field=float(tmp_array[1])
	    LAU_TOT_DIP = float(tmp_array[5])
	    
    if line =='  EV contribution to the polarization \n':
        ev_array=[]
        i=i+1
        #print i
        count = 0
        for _ in range(8):
            ln=f.next()
            #print ln
            row=[]
            row=ln.split()
            for j in range(5):
                ev=float(row[j])
                ev_array.append(ev)
                count = count +1

        #print count
        #print "count", i
    if line ==' Real primitive vectors (atomic units): \n':
        for _ in range(3):
                ln= f.next()
                #print ln
                row=[]
                row=ln.split()
                for j in [0,1,2]:
                        row[j]=float(row[j])  #here put check for float, if not exit and show error
                prim_mat.append(row)
#   prim_mat=np.array(prim_mat)
        #print prim_mat
        #prim_mat=np.matrix(np.array(prim_mat))
        prim_mat=np.array(prim_mat)
    if line == '  *** Coordinates and forces ***\n':
            cart_cord=[]
            for _ in range(nat):
                    ln= f.next()
                    row=[]
                    row=ln.split()
                    lst=[]
                    for j in [2,3,4]:
                            const=float(row[j])  #here put check for float, if not exit and show error
                            lst.append(const)
                    cart_cord.append(lst)
cart_cord=np.array(cart_cord)
ev_array=np.array(ev_array)
#print "wc_array", wc_array
f.close()
#print prim_mat
xcell=prim_mat[0][0]
area = prim_mat[1][1]* prim_mat[2][2]
volume = prim_mat[0][0]*area
#print xcell
##wc_ref=wc_array/(xcell)
##print "wc_ref"
##for i in range(160):
##    print wc_ref[i]
xcart=[]
for i in range(nat):
    var=cart_cord[i][0]
    xcart.append(var)
#print xcart
wc_array_sorted = sorted(wc_array)
xcart_sorted=sorted(xcart)
#print "sorted layer", wc_sorted

### Print xcart and wc_array to find out rescaling factor
print "compare xcart and wc_array to find out rescaling factor"
for i in range(40):
	print xcart_sorted[i],'\t',wc_array_sorted[i],'\t',wc_array_sorted[i+40],'\t', wc_array_sorted[i+80],'\t', wc_array_sorted[i+120] 
wc_ref=[0]*160
temp = 0
for i in range(160):
    if wc_array[i] >= -llimit:
       temp=float(wc_array[i])
       wc_ref[i] =  temp
    else:
        wc_ref[i] = wc_array[i] +   prim_mat[0][0]
#print "wc_ref******************", wc_ref
for i in range(160):
    if wc_ref[i]>xcell:
        wc_ref[i]=xcell - wc_ref[i]
    else:
        temp=float(wc_ref[i])
    wc_ref[i]=temp

#print wc_ref
wc_sorted =  sorted(wc_ref)
#print "sorted wc_ref", sorted(wc_ref)

### TESTING wc layers##
seq_1 = [16,24,16,24,16,24,16,24] 
seq_2 = [24,16,24,16,24,16,24,16] 
if first_layer == 'AO':
	seq = seq_1
	print "sequence of atoms are", seq
else:
	seq = seq_2
	print "sequence of atoms are", seq

j =0
k =0
demarcation=[0]*8
for k in range(8):
	j=j+seq[k]
	print 'j and j-1 ',j, j-1
	if j<159:
		l=j
		m=j-1
	if j==160:
		l=0
		m=j-1
	if wc_sorted[l]<0:
		num1=prim_mat[0][0] - wc_sorted[l]
	else:
		num1=wc_sorted[l]
	
	if wc_sorted[m]<0:
		num1=prim_mat[0][0] - wc_sorted[m]
	else:
		num2=wc_sorted[m]
	#demarcation[k]=(wc_sorted[l]+wc_sorted[m])/2.0
	demarcation[k]=np.mean([num1,num2])
	print wc_sorted[l] , wc_sorted[m] , demarcation[k]
	if (abs(num1 - num2) > abs(threshold)):
		print " Layer ", k, "has",seq[k], "  wannier centers\n" 
		k=k+1

	else:
		print " ERROR : Layer ", k, "does not have",seq[k], "  wannier centers\n" 
### Creating layer boundries 

layer1=np.array( [ demarcation[0] , demarcation[1] ])
layer2=np.array( [ demarcation[1] , demarcation[2] ])
layer3=np.array( [ demarcation[2] , demarcation[3] ])
layer4=np.array( [ demarcation[3] , demarcation[4] ])
layer5=np.array( [ demarcation[4] , demarcation[5] ])
layer6=np.array( [ demarcation[5] , demarcation[6] ])
layer7=np.array( [ demarcation[6] , demarcation[7] ])



lp_wc=[0]*8
wc_count=[0]*8
for j in range(160):
    #print j, wc_ref[j],wc_array[j]


    if (wc_ref[j] > layer1[0] and wc_ref[j] <= layer1[1]):
        lp_wc[1] = lp_wc[1] + spin_multiplicity*wc_array[j]
        #print "layer1", lp_wc[1], j, wc_ref[j],wc_array[j]
        wc_count[1] = wc_count[1] + 1

    elif (wc_ref[j] > layer2[0] and wc_ref[j] <= layer2[1]):
        lp_wc[2] = lp_wc[2] + spin_multiplicity*wc_array[j]
        wc_count[2] = wc_count[2] + 1

    elif (wc_ref[j] > layer3[0] and wc_ref[j] <= layer3[1]):
        lp_wc[3] = lp_wc[3] + spin_multiplicity*wc_array[j]
        wc_count[3] = wc_count[3] + 1

    elif (wc_ref[j] > layer4[0] and wc_ref[j] <= layer4[1]):
        lp_wc[4] = lp_wc[4] + spin_multiplicity*wc_array[j]
        wc_count[4] = wc_count[4] + 1

    elif (wc_ref[j] > layer5[0] and wc_ref[j] <= layer5[1]):
        lp_wc[5] = lp_wc[5] + spin_multiplicity*wc_array[j]
        wc_count[5] = wc_count[5] + 1

    elif (wc_ref[j] > layer6[0] and wc_ref[j] <= layer6[1]): 
        lp_wc[6] = lp_wc[6] + spin_multiplicity*wc_array[j]
        wc_count[6] = wc_count[6] + 1

#    if ((wc_ref[j] > layer7[0] and wc_ref[j] <= layer7[1]) or ((wc_ref[j] > layer8[0] and wc_ref[j] <= layer8[1]))):
    elif (wc_ref[j] > layer7[0] and wc_ref[j] <= layer7[1]):
        lp_wc[7] = lp_wc[7] + spin_multiplicity*wc_array[j]
        wc_count[7] = wc_count[7] + 1

#    if (wc_ref[j] > layer0[0] and wc_ref[j] <= layer0[1]):
#        print layer0 , wc_ref[j]
    else:
	lp_wc[0] = lp_wc[0] + spin_multiplicity*wc_array[j]
        wc_count[0] = wc_count[0] + 1
for i in range(8):
    print "wc_count",i, wc_count[i]
print "wc layer pol"
for i in range(8):
    print 'layer',i,"=", lp_wc[i]

#sum of wc layer polarization
tot_pol_wc = 0
for i in range(8):
    tot_pol_wc= tot_pol_wc + lp_wc[i]
print "TOTAL wc pol" , tot_pol_wc
tot = 0
for j in range(160):
    tot = tot + 2*wc_array[j]
print "tot direct sum", tot
#print "cart_cord", cart_cord
##
row=[]
iat=[]
for i in range(1):
    k=0
    for j in range(sp1[1]):
        k=k+1
        row=[k,sp1[2],sp1[0]]
        #print "sp1 are", j,iat
        iat.append(row)
    for j in range(sp2[1]):
        k=k+1
        row=[k,sp2[2],sp2[0]]        
        #print "sp2 are", j,iat
        iat.append(row)
    for j in range(sp3[1]):
        k=k+1
        row=[k,sp3[2],sp3[0]]        
        #print "sp3 are", j,iat
        iat.append(row)
    for j in range(sp4[1]):
        k=k+1
        row=[k,sp4[2],sp4[0]]        
        #print "sp4 are", j,iat
        iat.append(row)
#print "iat" ,iat
#for i in range(nat):
#    print iat[i][1]
xref=-1*(xcart-(xcell/2))/xcell
print "xref"
#for i in range(nat):
#    print xref[i]
#print "xref", xcart
xcart_rscale=np.zeros((nat,1))
#lp_ion=np.zeros((8,1))
lp_ion=[0]*8
ion_count=[0]*8
tot_pol_ion=0
layer_composition=['layer_']*8
layer_id=['']*8
layer_position=[0]*9
layer_tmp1=[]
layer_tmp0=[]
layer_tmp2=[]
layer_tmp3=[]
layer_tmp4=[]
layer_tmp5=[]
layer_tmp6=[]
layer_tmp7=[]
for j in range(nat):
    xcart_rscale[j]=0
    num=0
    num=xcart[j]/xcell
    #print j,num,int(round(num)), xcart_rscale[j]
    #xcart_rscale[j]=xcart[j]-float(int(round(num)))*xcell
    xcart_rscale[j] = xcart[j]*xcell


    if (xcart[j] > layer1[0] and xcart[j] <= layer1[1]):
        lp_ion[1] = lp_ion[1] + iat[j][1]*xcart[j]
        #print "layer1", lp_ion[1], j, xcart[j],xcart_rscale[j]
        ion_count[1]=ion_count[1]+1
        layer_composition[1]=layer_composition[1] + iat[j][2]
	layer_position[1]=layer_position[1]+xcart[j]
	dhr=xcart[j]
	layer_tmp1.append(dhr)
	if (iat[j][2]=='Ba'):
		layer_id[1]='B'
	if (iat[j][2]=='Ca'):
		layer_id[1]='C'
	if (iat[j][2]=='Ti'):
		layer_id[1]='T'


    elif (xcart[j] > layer2[0] and xcart[j] <= layer2[1]):
        lp_ion[2] = lp_ion[2] + iat[j][1]*xcart[j]
        ion_count[2]=ion_count[2]+1
        layer_composition[2]=layer_composition[2] + iat[j][2]
	layer_position[2]=layer_position[2]+xcart[j]
	dhr=xcart[j]
	layer_tmp2.append(dhr)
	if (iat[j][2]=='Ba'):
		layer_id[2]='B'
	if (iat[j][2]=='Ca'):
		layer_id[2]='C'
	if (iat[j][2]=='Ti'):
		layer_id[2]='T'

    elif (xcart[j] > layer3[0] and xcart[j] <= layer3[1]):
        lp_ion[3] = lp_ion[3] + iat[j][1]*xcart[j]
        ion_count[3]=ion_count[3]+1
        layer_composition[3]=layer_composition[3] + iat[j][2]
	layer_position[3]=layer_position[3]+xcart[j]
	dhr=xcart[j]
	layer_tmp3.append(dhr)
	if (iat[j][2]=='Ba'):
		layer_id[3]='B'
	if (iat[j][2]=='Ca'):
		layer_id[3]='C'
	if (iat[j][2]=='Ti'):
		layer_id[3]='T'

    elif (xcart[j] > layer4[0] and xcart[j] <= layer4[1]):
        lp_ion[4] = lp_ion[4] + iat[j][1]*xcart[j]
        ion_count[4]=ion_count[4]+1
        layer_composition[4]=layer_composition[4] + iat[j][2]
	layer_position[4]=layer_position[4]+xcart[j]
	dhr=xcart[j]
	layer_tmp4.append(dhr)
	if (iat[j][2]=='Ba'):
		layer_id[4]='B'
	if (iat[j][2]=='Ca'):
		layer_id[4]='C'
	if (iat[j][2]=='Ti'):
		layer_id[4]='T'

    elif (xcart[j] > layer5[0] and xcart[j] <= layer5[1]):
        lp_ion[5] = lp_ion[5] + iat[j][1]*xcart[j]
        ion_count[5]=ion_count[5]+1
        layer_composition[5]=layer_composition[5] + iat[j][2]
	layer_position[5]=layer_position[5]+xcart[j]
	dhr=xcart[j]
	layer_tmp5.append(dhr)
	if (iat[j][2]=='Ba'):
		layer_id[5]='B'
	if (iat[j][2]=='Ca'):
		layer_id[5]='C'
	if (iat[j][2]=='Ti'):
		layer_id[5]='T'

    elif (xcart[j] > layer6[0] and xcart[j] <= layer6[1]):
        lp_ion[6] = lp_ion[6] + iat[j][1]*xcart[j]
        ion_count[6]=ion_count[6]+1
        layer_composition[6]=layer_composition[6] + iat[j][2]
	layer_position[6]=layer_position[6]+xcart[j]
	dhr=xcart[j]
	layer_tmp6.append(dhr)
	if (iat[j][2]=='Ba'):
		layer_id[6]='B'
	if (iat[j][2]=='Ca'):
		layer_id[6]='C'
	if (iat[j][2]=='Ti'):
		layer_id[6]='T'

    elif (xcart[j] > layer7[0] and xcart[j] <= layer7[1]):
        lp_ion[7] = lp_ion[7] + iat[j][1]*xcart[j]
        ion_count[7]=ion_count[7]+1
        layer_composition[7]=layer_composition[7] + iat[j][2]
	layer_position[7]=layer_position[7]+xcart[j]
	dhr=xcart[j]
	layer_tmp7.append(dhr)
	if (iat[j][2]=='Ba'):
		layer_id[7]='B'
	if (iat[j][2]=='Ca'):
		layer_id[7]='C'
	if (iat[j][2]=='Ti'):
		layer_id[7]='T'

    else:
        lp_ion[0] = lp_ion[0] + iat[j][1]*xcart[j]
        ion_count[0]=ion_count[0]+1
        layer_composition[0]=layer_composition[0] + iat[j][2]
	xcartj = 0.
	xcartj  = xcart[j]
	if (xcart[j] >= prim_mat[0][0] +1. or xcart[j] <= prim_mat[0][0] -1):
		xcartj = xcartj - prim_mat[0][0]
		layer_position[0]=layer_position[0]+xcart[j]
	else:
		layer_position[0]=layer_position[0]+xcart[j]
	dhr=xcart[j]
	layer_tmp0.append(dhr)
	if (iat[j][2]=='Ba'):
		layer_id[0]='B'
	if (iat[j][2]=='Ca'):
		layer_id[0]='C'
	if (iat[j][2]=='Ti'):
		layer_id[0]='T'
    tot_pol_ion=tot_pol_ion+xcart[j]*iat[j][1]
sum=0
for i in range(8):
    print "ion_count",i, ion_count[i]
for i in range(8):
    sum=sum+lp_ion[i]
    print "lp_ion",i,lp_ion[i]
print "direct sum tot_pol_ion", tot_pol_ion
print "sum lp_ion", sum

for i in range(8):
	layer_position[i]=layer_position[i]/ion_count[i]
print "##################################################################"
print "##################################################################"
print "##################################################################"
for i in range(8):
	print "layer positon",i,"\t", layer_position[i]

print layer_tmp0
print layer_tmp1
print layer_tmp2
print layer_tmp3
print layer_tmp4
print layer_tmp5
print layer_tmp6
print layer_tmp7

# Layer 0 is near the box edge and may be behind the box edge
if layer_position[0] < 0.:
	layer_position[0] = layer_position[0] - prim_mat[0][0]
#all the layers are translated to bring layer 0 at zero postion.
layer_position_zero = layer_position[0]
for i in range(8):
	layer_position[i]=layer_position[i]-layer_position_zero
	print "tranlated layer postion" , i, layer_position[i]
#midle point of layer is calculated below
layer_mid_point=[0]*8
for i in range(7):
	layer_mid_point[i]=(layer_position[i+1] + layer_position[i])/2.
	print "layer mid point",i, layer_mid_point[i]
layer_mid_point[7] =( prim_mat[0][0] + layer_position[7])/2.
print "layer mid point",7, layer_mid_point[7]

layer_height=[0]*8
sum_layer_onetoseven = 0.
for i in range(1,8):
	layer_height[i]=(-layer_mid_point[i-1] + layer_mid_point[i])
	sum_layer_onetoseven = sum_layer_onetoseven + layer_height[i]
print "sum_layer_onetoseven" , sum_layer_onetoseven
layer_height[0] = prim_mat[0][0] - sum_layer_onetoseven

for i in range(8):
	print "Layer Hight ", i, layer_height[i]
### Redefine layer_position[0] 
##
##layer_position[8] = layer_position[8] + prim_mat[0][0]
##
### Layer heights are defined from the mid point of the layer positions
##
##layer_height=[0]*8
##sum_layer_onetoseven = 0.
##for i in range(1,8):
##	layer_height[i]=(-layer_position[i] + layer_position[i+1])
##	sum_layer_onetoseven = sum_layer_onetoseven + layer_height[i]
##print "sum_layer_onetoseven" , sum_layer_onetoseven
##layer_height[0] = prim_mat[0][0] - sum_layer_onetoseven
##
##for i in range(8):
##	print "Layer Hight ", i, layer_height[i]

#print "ev_array",ev_array

lp_ev=np.zeros((8,1))
tot_pol_ev=0
for j in range(nat):
    #print j, xcart[j],ev_array[j]


    if (xcart[j] > layer1[0] and xcart[j] <= layer1[1]):
        lp_ev[1] = lp_ev[1] + ev_array[j]
        #print "layer1", lp_ev[1], j, xcart[j],xcart_rscale[j]

    elif (xcart[j] > layer2[0] and xcart[j] <= layer2[1]):
        lp_ev[2] = lp_ev[2] + ev_array[j]

    elif (xcart[j] > layer3[0] and xcart[j] <= layer3[1]):
        lp_ev[3] = lp_ev[3] + ev_array[j]

    elif (xcart[j] > layer4[0] and xcart[j] <= layer4[1]):
        lp_ev[4] = lp_ev[4] + ev_array[j]

    elif (xcart[j] > layer5[0] and xcart[j] <= layer5[1]):
        lp_ev[5] = lp_ev[5] + ev_array[j]

    elif (xcart[j] > layer6[0] and xcart[j] <= layer6[1]):
        lp_ev[6] = lp_ev[6] + ev_array[j]

    elif (xcart[j] > layer7[0] and xcart[j] <= layer7[1]):
        lp_ev[7] = lp_ev[7] + ev_array[j]

    else:
        lp_ev[0] = lp_ev[0] + ev_array[j]
    tot_pol_ev=tot_pol_ev + ev_array[j]
sum=0
for i in range(8):
    sum=sum+lp_ev[i]
    print "lp_ev", i, lp_ev[i]
print "direct sum tot_pol_ev", tot_pol_ev
print "sum of lp_ev", sum


def substract_quantum(pol,xlattice):
    refpol = efieldx*xlattice
    izx=int(round((pol-refpol)/xlattice))
    pol=pol-float(izx)*xlattice
    return pol
for i in range(8):
    lp_wc[i]=substract_quantum(lp_wc[i],xcell)
    print i, "th layer lp_wc is",'\t',lp_wc[i]
    #tot_pol_wc = tot_pol_wc + lp_wc[i]
tot_pol_wc = substract_quantum(tot_pol_wc,xcell)
print "tot_pol_wc =",'\t',tot_pol_wc
for i in range(8):
    lp_ion[i]=substract_quantum(lp_ion[i],xcell)
    print i, "th layer lp_ion is",'\t',lp_ion[i]
    #tot_pol_ion = tot_pol_ion + lp_ion[i]
tot_pol_ion = substract_quantum(tot_pol_ion,xcell)
print "tot_pol_ion =",'\t',tot_pol_ion
for i in range(8):
    lp_ev[i]=substract_quantum(lp_ev[i],xcell)
    print i, "th layer lp_ev is",'\t',lp_ev[i]     
    #tot_pol_ev = tot_pol_ev + lp_ev[i]
tot_pol_ev = substract_quantum(tot_pol_ev,xcell)
print "tot_pol_ev =",'\t',tot_pol_ev
POL_TOT=substract_quantum(tot_pol_wc + tot_pol_ion + tot_pol_ev, xcell)
print POL_TOT
total_lp=[0]*8
#total Layer polarization
layer_series=['']*8
for i in range(8):
	if layer_id[i]=='T':
		layer_series[i]=layer_id[(i-3)%8] + layer_id[(i-1)%8] +layer_id[(i-0)%8] +layer_id[(i+1)%8] +layer_id[(i+3)%8] 
	elif (layer_id[i]=='B' or layer_id[i]=='C'):
		layer_series[i]=layer_id[(i-4)%8] + layer_id[(i-2)%8] +layer_id[(i-0)%8] +layer_id[(i+2)%8] +layer_id[(i+4)%8] 

for i in range(8):
    total_lp[i]= lp_wc[i] + lp_ion[i] + lp_ev[i]
    print "layer pol for layer",layer_composition[i] ,'\t',layer_series[i],'\t',i, '\t', total_lp[i]
#total Layer Polarization quantum substraction


###### normalizing the layer polarization to match total polarization#########
sum_total_lp = 0.
for i in range(8):
	sum_total_lp = sum_total_lp + total_lp[i] 
print "sum_total_lp ", sum_total_lp
scalling_const_lp = 0.
scalling_const_lp= POL_TOT/sum_total_lp
print   "scalling_const_lp",  scalling_const_lp
for i in range(8):
	total_lp[i] = total_lp[i] * scalling_const_lp
	print "layer pol for layer",layer_composition[i] ,'\t',layer_series[i],'\t',i, '\t', total_lp[i]

##################### correction for Kpoint sum difference ####################
pol_diff= LAU_TOT_DIP - POL_TOT
pol_diff_fraction = pol_diff/POL_TOT
POL_TOT = 0
for i in range(8):
	total_lp[i] = total_lp[i] + total_lp[i]*pol_diff_fraction
	POL_TOT = POL_TOT + total_lp[i]
        print "after_fraction pol_total layer pol for layer",layer_composition[i] ,'\t',layer_series[i],'\t',i, '\t', total_lp[i]
print "POL_TOT after fraction", POL_TOT 
#############################################################################


#for i in range(8):
out_file='layer_pol'+args.fname
f_out = open(out_file,'w') 
#SI-Units
unit_P = (1.60217653*10)/((0.5291772108)*(0.5291772108))
unit_LP = (1.60217653*10)/(0.5291772108)
for i in range(8):
    print "Layer Polarization ",'\t',prim_mat[1][1],prim_mat[2][2],prim_mat[0][0], volume, '\t\t', i, D_field,'\t',layer_series[i],'\t',str(total_lp[i]*unit_LP/area).strip('[]'),'\t','10d-10 C/m','\t', layer_composition[i]
print "POL_TOT",'\t', system, '\t' ,D_field,'\t',  POL_TOT*57.2147637564/volume
#print layer
f_out.write( "C_lattice \t A_lattice \t  layer_number \t D-field \t layer_series \t Layer_polarization  \t layer_composition   \t\t layer_height Bohr \n")

for i in range(8):
    f_out.write( " %f \t %f \t layer%d \t %f \t  %s  \t %s \t %s \t %s \t\t  %f \n "%(prim_mat[1][1],prim_mat[0][0], i, D_field,layer_series[i],str(total_lp[i]*unit_LP/area).strip('[]'),'10d-10 C/m', layer_composition[i], layer_height[i]))
f_out.write("\n\n \tsystem \t D_field(au)  D_field(SI)\t E_field(SI) \t\t POL_TOT(SI) \t Lattice(SI)\n" )
f_out.write(" TOT \t%s \t %4.3f \t\t %f \t %.8e \t %f \t %f  \n\n" %(system, D_field, D_field*57.2147637564/(4*3.14159265359) ,E_field*514220652000,   POL_TOT*57.2147637564/volume , prim_mat[0][0])) #*0.529177249))
#f_out.write("\n\n \tsystem \t D_field(au)  \t E_field(SI) \t\t POL_TOT(SI) \t Lattice(SI)\n" )
#f_out.write(" TOT \t%s \t %4.3f \t\t %.8e \t %f \t %f  \n\n" %(system, D_field, E_field*514220652000,   POL_TOT*57.2147637564/volume , prim_mat[0][0]*0.529177249))
f_out.write(" TOT atomic units D field | E_field_au_LAU | LAU_TOT_DIP | Efield | TOT_DIP_WC | volume \t %4.3f \t\t  %f \t\t  %f \t\t %f \t\t %f \t\t %f \n\n" %(D_field, E_field  , LAU_TOT_DIP,  float(D_field - 4*3.14159265359*POL_TOT/volume), POL_TOT, volume))
#f_out.write("\n\n \tsystem \t D_field(au)  \t E_field(SI) \t\t POL_TOT(SI) \n" )
#f_out.write(" TOT \t%s \t %4.3f \t\t %.8e \t %f \n\n" %(system, D_field, E_field*514220652000,   POL_TOT*57.2147637564/volume))
#    np.savetxt(f_out,prim_mat[0][0], D_field,layer_series[i],str(total_lp[i]*unit_LP/area).strip('[]'),'10d-10 C/m', layer_composition[i])
#f_out.write( "POL_TOT",'\t', system, '\t' ,D_field,'\t',  POL_TOT*57.2147637564/volume)
#f_out.close()




### followin scripts writes POSCAR file in cartesian coordinates
file_prefix= system + '_' + str(D_field)
prefix= repr(file_prefix).replace("0.",r"-tilt_Dfield_")
print ''.join([prefix.strip("'"),'.poscar.vasp'])
prefix= ''.join([prefix.strip("'"),'.poscar.vasp'])
print prefix
ang2bohr = 0.529177249
f_poscar= open (prefix,'w')
f_poscar.write("generated with lp-script.py from lautrec output file \n")
f_poscar.write("1.0\n")
cart_prim_mat= prim_mat*ang2bohr
np.savetxt(f_poscar, cart_prim_mat, fmt='%10.7f')
spieces = sp1[0]+ "   " +sp2[0]+ "   " +sp3[0]+ "   " +sp4[0]
f_poscar.write(spieces)
f_poscar.write("\n %d     %d      %d      %d \n" %(sp1[1],sp2[1],sp3[1],sp4[1]))
cord_in_cartesian = cart_cord * ang2bohr
f_poscar.write("Cartesian \n")
np.savetxt(f_poscar, cord_in_cartesian, fmt='%10.7f')
f_poscar.close()
argument= ''.join(['python vasp2cif.py',' ',prefix])
print argument
os.system(argument.strip("'"))
