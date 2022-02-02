from tkinter import *
from tkinter import ttk
from openpyxl.workbook import Workbook
from openpyxl import load_workbook
import numpy
import csv
import pandas as pd
from tkinter import ttk, filedialog
import sqlite3
from PIL import ImageTk, Image
from tkinter.filedialog import askopenfile
from tkinter.filedialog import asksaveasfile
import os
import re
import math

###########################################################################################
def run_ruscal():
	os.system("ruscal")
	return 0

###########################################################################################
# Function to write to rusin file
def writeFile():
	file = open('rusin.dat','w') # Compile on Mac: use .txt and use a+
	head = Name_text.get()
	file.write("# " + head + '\n')
	file.write('# Shape     Bravais   n-poly    #-calc     Mirror   Mass(g)  LastRMS CalcDeriv' + '\n')
	shapes_text = shapeclicked.get()
	brav_text = bravclass_clicked.get()
	brav_val = re.search(r"\[([A-Za-z0-9_]+)\]", brav_text) #error will appear if no Bravais lattice is chosen
	brav_val2 = int(brav_val.group(1))
	maxpoly_text = maxpoly_clicked.get()
	mirror_text = mirror_clicked.get()
	selection = modes.get()
	if(selection == -2):
		dimvar=vary_clicked.get()
	weight_text = str(weight.get())
	lastchi2_text = lastchi2.get()
	if(lastchi2_text==""):
		lastchi2_text="1.00"
	calc_text = str(calcmodevar.get())
	file.write('  ' +shapes_text[1]+ '         ' + brav_val.group(1)+ '         ' +maxpoly_text+ '        '+str(calcmodefunc())+'         ' +mirror_text[0]+ '        '+ weight_text+'        '+lastchi2_text+'        '+calc_text+ '\n')
	file.write('# Elastic constants (2nd line: n constraints (0/1))' + '\n')


# Case 0 - Fit
# Case 1 - Calculate
# Case -2 - Vary dimension or angle
# Case -3 - Get error analysis
	if (selection == 0 or selection == 1 or selection == -2 or selection == -3):
		if brav_val.group(1) == '2':
			file.write(str(c11_box.get())+' '+str(c44_box.get())+ '\n')
			file.write(str(checkedbox11())+' '+str(checkedbox44())+ '\n')
		if  brav_val.group(1) == '3':
			file.write(str(c11_box.get())+'    '+str(c12_box.get())+'    '+str(c44_box.get())+ '\n')
			file.write(str(checkedbox11())+' '+str(checkedbox12())+' ' +str(checkedbox44())+ '\n')
		if  brav_val.group(1) == '5':
			file.write(str(c33_box.get())+' '+str(c23_box.get())+' '+str(c12_box.get())+' '+str(c44_box.get())+' '+str(c66_box.get())+ '\n')
			file.write(str(checkedbox33())+' '+str(checkedbox23())+' ' +str(checkedbox12())+' ' +str(checkedbox44())+' '+ str(checkedbox66())+ '\n')
		if  brav_val.group(1) == '306':
			file.write(str(c11_box.get())+' '+str(c33_box.get())+' '+str(c23_box.get())+' '+str(c12_box.get())+' '+str(c44_box.get())+' '+str(c14_box.get())+ '\n')
			file.write(str(checkedbox11())+ ' '+str(checkedbox33())+' '+str(checkedbox23())+' '+str(checkedbox12())+' '+str(checkedbox44())+' '+str(checkedbox14())+ '\n')
		if  brav_val.group(1) == '307':
			file.write(str(c11_box.get())+' '+str(c33_box.get())+' '+str(c13_box.get())+' '+str(c12_box.get())+' '+str(c44_box.get())+' '+str(c14_box.get())+' '+str(c25_box.get())+ '\n')
			file.write(str(checkedbox11())+ ' '+str(checkedbox33())+' '+str(checkedbox13())+' '+str(checkedbox12())+' '+str(checkedbox44())+' '+str(checkedbox14())+' '+str(checkedbox25())+ '\n')
		if  brav_val.group(1) == '406':
			file.write(str(c11_box.get())+' '+str(c33_box.get())+' '+str(c23_box.get())+' '+str(c12_box.get())+' '+str(c44_box.get())+' '+str(c66_box.get())+ '\n')
			file.write(str(checkedbox11())+ ' '+str(checkedbox33())+' '+str(checkedbox23())+' '+str(checkedbox12())+' '+str(checkedbox44())+' '+str(checkedbox66())+ '\n')
		if  brav_val.group(1) == '407':
			file.write(str(c11_box.get())+' '+str(c33_box.get())+' '+str(c23_box.get())+' '+str(c12_box.get())+' '+str(c44_box.get())+' '+str(c66_box.get())+' '+str(c16_box.get())+ '\n')
			file.write(str(checkedbox11())+ ' '+str(checkedbox33())+' '+str(checkedbox23())+' '+str(checkedbox12())+' '+str(checkedbox44())+' '+str(checkedbox66())+' '+str(checkedbox16())+ '\n')
		if  brav_val.group(1) == '9':
			file.write(str(c11_box.get())+' '+str(c22_box.get())+' '+str(c33_box.get())+' '+str(c23_box.get())+' '+str(c13_box.get())+' '+str(c12_box.get())+' '+str(c44_box.get())+' '+str(c55_box.get())+' '+str(c66_box.get())+ '\n')
			file.write(str(checkedbox11())+ ' '+str(checkedbox22())+' '+str(checkedbox33())+' '+str(checkedbox23())+' '+str(checkedbox13())+' '+str(checkedbox12())+' '+str(checkedbox44())+' '+str(checkedbox55())+' '+str(checkedbox66())+ '\n')
		if  brav_val.group(1) == '13':
			file.write(str(c11_box.get())+' '+str(c22_box.get())+' '+str(c33_box.get())+' '+str(c23_box.get())+' '+str(c13_box.get())+' '+str(c12_box.get())+' '+str(c44_box.get())+' '+str(c55_box.get())+' '+str(c66_box.get())+' '+str(c15_box.get())+' '+str(c25_box.get())+' '+str(c35_box.get())+' '+str(c46_box.get())+ '\n')
			file.write(str(checkedbox11())+ ' '+str(checkedbox22())+' '+str(checkedbox33())+' '+str(checkedbox23())+' '+str(checkedbox13())+' '+str(checkedbox12())+' '+str(checkedbox44())+' '+str(checkedbox55())+' '+str(checkedbox66())+' '+str(checkedbox15())+' '+str(checkedbox25())+' '+str(checkedbox35())+' '+str(checkedbox46())+ '\n')
		if  brav_val.group(1) == '21':
			file.write(str(c11_box.get())+' '+str(c22_box.get())+' '+str(c33_box.get())+' '+str(c23_box.get())+' '+str(c13_box.get())+' '+str(c12_box.get())+' '+str(c44_box.get())+' '+str(c55_box.get())+' '+str(c66_box.get())+' '+str(c14_box.get())+' '+str(c15_box.get())+' '+str(c16_box.get())+' '+str(c24_box.get())+' '+str(c25_box.get())+' '+str(c26_box.get())+' '+str(c34_box.get())+' '+str(c35_box.get())+' '+str(c36_box.get())+' '+str(c45_box.get())+' '+str(c46_box.get())+' '+str(c56_box.get())+' \n')
			file.write(str(checkedbox11())+ ' '+str(checkedbox22())+' '+str(checkedbox33())+' '+str(checkedbox23())+' '+str(checkedbox13())+' '+str(checkedbox12())+' '+str(checkedbox44())+' '+str(checkedbox55())+' '+str(checkedbox66())+' '+str(checkedbox14())+' '+str(checkedbox15())+' '+str(checkedbox16())+' '+str(checkedbox24())+' '+str(checkedbox25())+' '+str(checkedbox26())+' '+str(checkedbox34())+' '+str(checkedbox35())+' '+str(checkedbox36())+' '+str(checkedbox45())+' '+str(checkedbox46())+' '+str(checkedbox56())+' \n')

#Case -1 Monte Carlo
#Case -4 Grid fit
	if (selection == -1 or selection == -4):
		c11_bd_val=int("0"+c11_bd.get())+0
		c22_bd_val=int("0"+c22_bd.get())+0
		c33_bd_val=int("0"+c33_bd.get())+0
		c23_bd_val=int("0"+c23_bd.get())+0
		c13_bd_val=int("0"+c13_bd.get())+0
		c12_bd_val=int("0"+c12_bd.get())+0
		c44_bd_val=int("0"+c44_bd.get())+0
		c55_bd_val=int("0"+c55_bd.get())+0
		c66_bd_val=int("0"+c66_bd.get())+0
		c14_bd_val=int("0"+c14_bd.get())+0
		c15_bd_val=int("0"+c15_bd.get())+0
		c16_bd_val=int("0"+c16_bd.get())+0
		c24_bd_val=int("0"+c24_bd.get())+0
		c25_bd_val=int("0"+c25_bd.get())+0
		c26_bd_val=int("0"+c26_bd.get())+0
		c34_bd_val=int("0"+c34_bd.get())+0
		c35_bd_val=int("0"+c35_bd.get())+0
		c36_bd_val=int("0"+c36_bd.get())+0
		c45_bd_val=int("0"+c45_bd.get())+0
		c46_bd_val=int("0"+c46_bd.get())+0
		c56_bd_val=int("0"+c56_bd.get())+0
		if brav_val.group(1) == '2':
			file.write(str(c11_box.get())+' '+str(c44_box.get())+ '\n')
			file.write(str(c11_bd_val)+' '+str(c44_bd_val)+ '\n')
		if brav_val.group(1) == '3':
			file.write(str(c11_box.get())+' '+str(c12_box.get())+' '+str(c44_box.get())+ '\n')
			file.write(str(c11_bd_val)+' '+str(c12_bd_val)+' '+str(c44_bd_val)+ '\n')
		if brav_val.group(1) == '5':
			file.write(str(c33_box.get())+' '+str(c23_box.get())+' '+str(c12_box.get())+' '+str(c44_box.get())+' '+str(c66_box.get())+ '\n')
			file.write(str(c33_bd_val)+' '+str(c23_bd_val)+' '+str(c12_bd_val)+' '+str(c44_bd_val)+' '+str(c66_bd_val)+ '\n')
		if brav_val.group(1) == '306':
			file.write(str(c11_box.get())+' '+str(c33_box.get())+' '+str(c23_box.get())+' '+str(c12_box.get())+' '+str(c44_box.get())+' '+str(c14_box.get())+ '\n')
			file.write(str(c11_bd_val)+' '+str(c33_bd_val)+' '+str(c23_bd_val)+' '+str(c12_bd_val)+' '+str(c44_bd_val)+' '+str(c14_bd_val)+ '\n')
		if brav_val.group(1) == '307':
			file.write(str(c11_box.get())+' '+str(c33_box.get())+' '+str(c13_box.get())+' '+str(c12_box.get())+' '+str(c44_box.get())+' '+str(c14_box.get())+' '+str(c25_box.get())+ '\n')
			file.write(str(c11_bd_val)+' '+str(c33_bd_val)+' '+str(c13_bd_val)+' '+str(c12_bd_val)+' '+str(c44_bd_val)+' '+str(c14_bd_val)+' '+str(c25_bd_val)+ '\n')
		if brav_val.group(1) == '406':
			file.write(str(c11_box.get())+' '+str(c33_box.get())+' '+str(c23_box.get())+' '+str(c12_box.get())+' '+str(c44_box.get())+' '+str(c66_box.get())+ '\n')
			file.write(str(c11_bd_val)+' '+str(c33_bd_val)+' '+str(c23_bd_val)+' '+str(c12_bd_val)+' '+str(c44_bd_val)+' '+str(c66_bd_val)+ '\n')
		if brav_val.group(1) == '407':
			file.write(str(c11_box.get())+' '+str(c33_box.get())+' '+str(c23_box.get())+' '+str(c12_box.get())+' '+str(c44_box.get())+' '+str(c66_box.get())+' '+str(c16_box.get())+ '\n')
			file.write(str(c11_bd_val)+' '+str(c33_bd_val)+' '+str(c23_bd_val)+' '+str(c12_bd_val)+' '+str(c44_bd_val)+' '+str(c66_bd_val)+' '+str(c16_bd_val)+ '\n')
		if brav_val.group(1) == '9':
			file.write(str(c11_box.get())+' '+str(c22_box.get())+' '+str(c33_box.get())+' '+str(c23_box.get())+' '+str(c13_box.get())+' '+str(c12_box.get())+' '+str(c44_box.get())+' '+str(c55_box.get())+' '+str(c66_box.get())+ '\n')
			file.write(str(c11_bd_val)+' '+str(c22_bd_val)+' '+str(c33_bd_val)+' '+str(c23_bd_val)+' '+str(c13_bd_val)+' '+str(c12_bd_val)+' '+str(c44_bd_val)+' '+str(c55_bd_val)+' '+str(c66_bd_val)+ '\n')
		if brav_val.group(1) == '13':
			file.write(str(c11_box.get())+' '+str(c22_box.get())+' '+str(c33_box.get())+' '+str(c23_box.get())+' '+str(c13_box.get())+' '+str(c12_box.get())+' '+str(c44_box.get())+' '+str(c55_box.get())+' '+str(c66_box.get())+' '+str(c15_box.get())+' '+str(c25_box.get())+' '+str(c35_box.get())+' '+str(c46_box.get())+ '\n')
			file.write(str(c11_bd_val)+' '+str(c22_bd_val)+' '+str(c33_bd_val)+' '+str(c23_bd_val)+' '+str(c13_bd_val)+' '+str(c12_bd_val)+' '+str(c44_bd_val)+' '+str(c55_bd_val)+' '+str(c66_bd_val)+' '+str(c15_bd_val)+' '+str(c25_bd_val)+' '+str(c35_bd_val)+' '+str(c46_bd_val)+ '\n')
		if brav_val.group(1) == '21':
			file.write(str(c11_box.get())+' '+str(c22_box.get())+' '+str(c33_box.get())+' '+str(c23_box.get())+' '+str(c13_box.get())+' '+str(c12_box.get())+' '+str(c44_box.get())+' '+str(c55_box.get())+' '+str(c66_box.get())+' '+str(c14_box.get())+' '+str(c15_box.get())+' '+str(c16_box.get())+' '+str(c24_box.get())+' '+str(c25_box.get())+' '+str(c26_box.get())+' '+str(c34_box.get())+' '+str(c35_box.get())+' '+str(c36_box.get())+' '+str(c46_box.get())+' '+str(c56_box.get())+ '\n')
			file.write(str(c11_bd_val)+' '+str(c22_bd_val)+' '+str(c33_bd_val)+' '+str(c23_bd_val)+' '+str(c13_bd_val)+' '+str(c12_bd_val)+' '+str(c44_bd_val)+' '+str(c55_bd_val)+' '+str(c66_bd_val)+' '+str(c14_bd_val)+' '+str(c15_bd_val)+' '+str(c16_bd_val)+' '+str(c24_bd_val)+' '+str(c25_bd_val)+' '+str(c26_bd_val)+' '+str(c34_bd_val)+' '+str(c35_bd_val)+' '+str(c36_bd_val)+' '+str(c45_bd_val)+' '+str(c46_bd_val)+' '+str(c56_bd_val)+ '\n')
	
	file.write('# Dimensions (2nd line: 1 constraint (0/1))' + '\n')
	if shapeclicked.get() == shapes[0] or shapeclicked.get() == shapes[1] or shapeclicked.get() == shapes[3] or shapeclicked.get() == shapes[4] or shapeclicked.get() == shapes[5]:
		file.write(str(entry1.get())+' '+str(entry2.get())+' '+str(entry3.get())+ '\n')

	if shapeclicked.get() == shapes[2]:
		file.write(str(potatox.get())+' '+str(potatoy.get())+'  '+str(potatoz.get())+' '+str(negx.get())+' '+str(negy.get())+' '+str(negz.get())+ '\n')

	if(selection==-2):
		file.write(str(dimvar[1])+ '\n')
	else:
		file.write(str(fixdims())+ '\n')

	file.write('# Euler Angles (2nd line: 3 constraints (0/1))' + '\n')
	file.write(str(eulerb1.get())+' '+str(eulerb2.get())+' '+str(eulerb3.get())+ '\n')
	file.write(str(eulerval1())+' '+str(eulerval2())+' ' +str(eulerval3())+ '\n')

	file.write('# F_ex(MHz) F_th(MHz) Weight'+ '\n')

	hasfreqs=0
	for child in treef.get_children():
		row = treef.item(child)["values"]
		file.write("%8.6f    %8.6f    %4.2f\n"%(float(row[1]),float(row[2]),float(row[3])))
		hasfreqs=hasfreqs+1

	if(hasfreqs==0):
		file.write('0.000000 0.000000 0.00\n')

	con_text = []

	file.close()

def readFile():
	s=label_inputfile["text"]
	if s=='':
		s='rusio.dat'
	#'rusio.dat'
	inputfile = open(s,'rt') # Compile on Mac: use .txt and use a+
	lines = inputfile.readlines()
	inputfile.close()
	count=0
	Name_text.set((lines[0])[1:].strip()) # Use first line as Sample Name, strip first character.
	param=[]
	for line in lines:
		if(line[0]!="#"):
			count += 1
#			print(f'line {count}: {re.split(" +",line.strip())}') 
			param.append(re.split(" +",line.strip()))    

	Name_text.set((lines[0])[1:].strip())
	shapeclicked.set(shapes[int(shapetab[int(param[0][0])])])
	bravclass_clicked.set(bravaisclass[tabbrav(int(param[0][1]))])
	maxpoly_clicked.set(maxpolyorder[int(param[0][2])-8])
	if(int(param[0][3])<1):
		modes.set(int(param[0][3]))
	else:
		modes.set(1)
		calcmoden.set(int(param[0][3]))
		
	if(int(param[0][3])==-2):
		vary_clicked.set(int(param[4][0]))
	else:
		dim.set(int(param[4][0])>0)
		
	mirror_clicked.set(mirrorplanes[int(param[0][4])])
	weightentry.delete(0,END)
	weightentry.insert(0,(param[0][5]))
	lastchi2.set(param[0][6])
	calcmodevar.set(int(param[0][7]))
	
	if(int(param[0][0])==4):
		potatox.config(state='normal')
		potatoy.config(state='normal')
		potatoz.config(state='normal')
		x.config(state='normal')
		y.config(state='normal')
		z.config(state='normal')
		negx.config(state='normal')
		negy.config(state='normal')
		negz.config(state='normal')
		negx_label.config(state='normal')
		negy_label.config(state='normal')
		negz_label.config(state='normal')
		potatox.delete(0,END)
		potatox.insert(0,param[3][0])
		potatoy.delete(0,END)
		potatoy.insert(0,param[3][1])
		potatoz.delete(0,END)
		potatoz.insert(0,param[3][2])
		negx.delete(0,END)
		negx.insert(0,param[3][3])
		negy.delete(0,END)
		negy.insert(0,param[3][4])
		negz.delete(0,END)
		negz.insert(0,param[3][5])
	else:
		entry1.config(state='normal')
		entry2.config(state='normal')
		entry3.config(state='normal')
		entry1.delete(0,END)
		entry1.insert(0,param[3][0])
		entry2.delete(0,END)
		entry2.insert(0,param[3][1])
		entry3.delete(0,END)
		entry3.insert(0,param[3][2])

	shapeshow(1) # This is needed twice!
	shapeshow(1)
	
	euler1.set(float(param[5][0]))
	euler2.set(float(param[5][1]))
	euler3.set(float(param[5][2]))
	eulervarpsi.set(int(param[6][0]))
	eulervartheta.set(int(param[6][1]))
	eulervarphi.set(int(param[6][2]))
	
	showbravais(bravaisclass[tabbrav(int(param[0][1]))])
	bdry=0
	if(int(param[0][3])==-1 or int(param[0][3])==-4):
		bdry=1
	
# Param[1][0] ---- Param[1][n] - read elastic tensor and constraints
	if int(param[0][1]) == 2:
		c11.set(param[1][0]); c44.set(param[1][1])
		if(bdry):
			c11_bdry.set(param[2][0]);c44_bdry.set(param[2][1])
		else:
			cb11.set(param[2][0]); cb44.set(param[2][1])
	if int(param[0][1]) == 3:
		c11.set(param[1][0]); c12.set(param[1][1]); c44.set(param[1][2])
		if(bdry):
			c11_bdry.set(param[2][0]); c12_bdry.set(param[2][1]); c44_bdry.set(param[2][2])
		else:
			cb11.set(param[2][0]); cb12.set(param[2][1]); cb44.set(param[2][2])
	if int(param[0][1]) == 5:
		c33.set(param[1][0]); c23.set(param[1][1]); c12.set(param[1][2]); c44.set(param[1][3]); c66.set(param[1][4])
		if(bdry):
			c33_bdry.set(param[2][0]); c23_bdry.set(param[2][1]);c12_bdry.set(param[2][2]); c44_bdry.set(param[2][3]); c66_bdry.set(param[2][4])
		else:
			cb33.set(param[2][0]); cb23.set(param[2][1]); cb12.set(param[2][2]); cb44.set(param[2][3]); cb66.set(param[2][4])
	if int(param[0][1]) == 406:
		c11.set(param[1][0]); c33.set(param[1][1]); c23.set(param[1][2]); c12.set(param[1][3]); c44.set(param[1][4]); c66.set(param[1][5])
		if(bdry):
			c11_bdry.set(param[2][0]); c33_bdry.set(param[2][1]); c23_bdry.set(param[2][2]);c12_bdry.set(param[2][3]); c44_bdry.set(param[2][4]); c66_bdry.set(param[2][5])
		else:
			cb11.set(param[2][0]); cb33.set(param[2][1]); cb23.set(param[2][2]); cb12.set(param[2][3]); cb44.set(param[2][4]); cb66.set(param[2][5])	
	if int(param[0][1]) == 407:
		c11.set(param[1][0]); c33.set(param[1][1]); c23.set(param[1][2]); c12.set(param[1][3]); c44.set(param[1][4]); c66.set(param[1][5]); c16.set(param[1][6])
		if(bdry):
			c11_bdry.set(param[2][0]); c33_bdry.set(param[2][1]); c23_bdry.set(param[2][2]); c12_bdry.set(param[2][3]); c44_bdry.set(param[2][4]); c66_bdry.set(param[2][5]); c16_bdry.set(param[2][6])
		else:
			cb11.set(param[2][0]); cb33.set(param[2][1]); cb23.set(param[2][2]); cb12.set(param[2][3]); cb44.set(param[2][4]); cb66.set(param[2][5]); cb16.set(param[2][6])
	if int(param[0][1]) == 306:
		c11.set(param[1][0]); c33.set(param[1][1]); c23.set(param[1][2]); c12.set(param[1][3]); c44.set(param[1][4]); c14.set(param[1][5])
		if(bdry):
			c11_bdry.set(param[2][0]); c33_bdry.set(param[2][1]); c23_bdry.set(param[2][2]); c12_bdry.set(param[2][3]); c44_bdry.set(param[2][4]); c14_bdry.set(param[2][5])
		else:
			cb11.set(param[2][0]); cb33.set(param[2][1]); cb23.set(param[2][2]); cb12.set(param[2][3]); cb44.set(param[2][4]); cb14.set(param[2][5]); 
	if int(param[0][1]) == 307:
		c11.set(param[1][0]); c33.set(param[1][1]); c23.set(param[1][2]); c12.set(param[1][3]); c44.set(param[1][4]); c14.set(param[1][5]); c25.set(param[1][6])
		if(bdry):
			c11_bdry.set(param[2][0]); c33_bdry.set(param[2][1]); c23_bdry.set(param[2][2]); c12_bdry.set(param[2][3]); c44_bdry.set(param[2][4]); c14_bdry.set(param[2][5]); c25_bdry.set(param[2][6])
		else:
			cb11.set(param[2][0]); cb33.set(param[2][1]); cb23.set(param[2][2]); cb12.set(param[2][3]); cb44.set(param[2][4]); cb14.set(param[2][5]); cb25.set(param[2][6])
	if (int(param[0][1])%100) >8:
		c11.set(param[1][0]); c22.set(param[1][1]); c33.set(param[1][2])
		c23.set(param[1][3]); c13.set(param[1][4]); c12.set(param[1][5])
		c44.set(param[1][6]); c55.set(param[1][7]); c66.set(param[1][8])
		if(bdry):
			c11_bdry.set(param[2][0]); c22_bdry.set(param[2][1]); c33_bdry.set(param[2][2])
			c23_bdry.set(param[2][3]); c13_bdry.set(param[2][4]); c12_bdry.set(param[2][5])
			c44_bdry.set(param[2][6]); c55_bdry.set(param[2][7]); c66_bdry.set(param[2][8])
		else:
			cb11.set(param[2][0]); cb22.set(param[2][1]); cb33.set(param[2][2])
			cb23.set(param[2][3]); cb13.set(param[2][4]); cb12.set(param[2][5])
			cb44.set(param[2][6]); cb55.set(param[2][7]); cb66.set(param[2][8])
	if (int(param[0][1])) == 13:
		c15.set(param[1][9]); c25.set(param[1][10]); c35.set(param[1][11])
		c46.set(param[1][12])
		if(bdry):
			c15_bdry.set(param[2][9]); c25_bdry.set(param[2][10]); c35_bdry.set(param[2][11])
			c46_bdry.set(param[2][12])
		else:
			cb15.set(param[2][9]); cb25.set(param[2][10]); cb35.set(param[2][11])
			cb46.set(param[2][12])
	if (int(param[0][1])) ==21:
		c14.set(param[1][9]); c15.set(param[1][10]);c16.set(param[1][11])
		c24.set(param[1][12]);c25.set(param[1][13]);c26.set(param[1][14]);
		c34.set(param[1][15]);c35.set(param[1][16]);c36.set(param[1][17]);
		c45.set(param[1][18]);c46.set(param[1][19]);c56.set(param[1][20])
		if(bdry):
			c14_bdry.set(param[2][9]); c15_bdry.set(param[2][10]);c16_bdry.set(param[2][11])
			c24_bdry.set(param[2][12]);c25_bdry.set(param[2][13]);c26_bdry.set(param[2][14])
			c34_bdry.set(param[2][15]);c35_bdry.set(param[2][16]);c36_bdry.set(param[2][17])
			c45_bdry.set(param[2][18]);c46_bdry.set(param[2][19]);c56_bdry.set(param[2][20])
		else:
			cbd14.set(param[2][9]); cbd15.set(param[2][10]);cbd16.set(param[2][11])
			cbd24.set(param[2][12]);cbd25.set(param[2][13]);cbd26.set(param[2][14])
			cbd34.set(param[2][15]);cbd35.set(param[2][16]);cbd36.set(param[2][17])
			cbd45.set(param[2][18]);cbd46.set(param[2][19]);cbd56.set(param[2][20])
			
	treef.delete(*treef.get_children())
	for i in range(7,count):
		if(param[i]!=['']):
			treef.insert("",'end', values=(i-6,float(param[i][0]),float(param[i][1]),float(param[i][2])))
#				i=i+1
	freqs_file.set("Choose File Type")
	load_freqs(1)


# Deal with chi2
# Read in frequencies

# Need to invoke all conditionals at the end. 
	if(int(param[0][3])==0):
		val0.invoke()
	if(int(param[0][3])>0):
		valcalc1.invoke()
	if(int(param[0][3])==-1):
		valneg1.invoke()
	if(int(param[0][3])==-2):
		valneg2.invoke()
	if(int(param[0][3])==-3):
		valneg3.invoke()
	if(int(param[0][3])==-4):
		valneg4.invoke()
 
###########################################################################################
###########################################################################################
# Main Root and Frames:

root = Tk()
root.title("RUS Analysis v0.95")
root.geometry("1280x1024")

# Create a main_frame
main_frame = Frame(root)
main_frame.pack(fill=BOTH, expand=1)

menubar = Menu(root)
root.config(menu=menubar)

# Create Canvas
canvas = Canvas(main_frame)
canvas.pack(side=LEFT, fill=BOTH, expand=1)

# Scrollbar
scrollbar = ttk.Scrollbar(main_frame, orient=VERTICAL, command=canvas.yview)
scrollbar.pack(side=RIGHT, fill=Y)

canvas.configure(yscrollcommand=scrollbar.set)
canvas.bind('<Configure>', lambda e: canvas.configure(scrollregion = canvas.bbox("all")))

# Setup consta
second_frame = Frame(canvas)
canvas.create_window((0,0), window=second_frame, anchor="nw")

# Top right elastic constants Frame
c_frame = LabelFrame(second_frame, text="Elastic Constants (Mbar = 100 GPa units)", padx=10, pady=10)
c_frame.grid(row=1, rowspan=31, column=5, columnspan=10, padx=3, pady=3)

# Bottom status bar
status_bar = Label(root, text='Ready        ', anchor=W)
status_bar.pack(fill=X, side=BOTTOM, ipady=5)

###########################################################################################
###########################################################################################
# Functions for writing to rusin text file for check boxes:

cb11 = IntVar()
cb12 = IntVar()
cb13 = IntVar()
cb14 = IntVar()
cb15 = IntVar()
cb16 = IntVar()
cb22 = IntVar()
cb23 = IntVar()
cb24 = IntVar()
cb25 = IntVar()
cb26 = IntVar()
cb33 = IntVar()
cb34 = IntVar()
cb35 = IntVar()
cb36 = IntVar()
cb44 = IntVar()
cb45 = IntVar()
cb46 = IntVar()
cb55 = IntVar()
cb56 = IntVar()
cb66 = IntVar()


def checkedbox11():
	if chk11["state"] == "normal" and cb11.get() == 1:
		return '1'
	else:
		return '0'

def checkedbox12():
	if chk12["state"] == "normal" and cb12.get() == 1:
		return '1'
	else:
		return '0'
	
def checkedbox13():
	if chk13["state"] == "normal" and cb13.get() == 1:
		return '1'
	else:
		return '0'
	
def checkedbox14():
	if chk14["state"] == "normal" and cb14.get() == 1:
		return '1'
	else:
		return '0'
	
def checkedbox15():
	if chk15["state"] == "normal" and cb15.get() == 1:
		return '1'
	else:
		return '0'
	
def checkedbox16():
	if chk16["state"] == "normal" and cb16.get() == 1:
		return '1'
	else:
		return '0'
	

def checkedbox22():
	if chk22["state"] == "normal" and cb22.get() == 1:
		return '1'
	else:
		return '0'
	
def checkedbox23():
	if chk23["state"] == "normal" and cb23.get() == 1:
		return '1'
	else:
		return '0'
	
def checkedbox24():
	if chk24["state"] == "normal" and cb24.get() == 1:
		return '1'
	else:
		return '0'
	
def checkedbox25():
	if chk25["state"] == "normal" and cb25.get() == 1:
		return '1'
	else:
		return '0'
	
def checkedbox26():
	if chk26["state"] == "normal" and cb26.get() == 1:
		return '1'
	else:
		return '0'

def checkedbox33():
	if chk33["state"] == "normal" and cb33.get() == 1:
		return '1'
	else:
		return '0'

def checkedbox34():
	if chk34["state"] == "normal" and cb34.get() == 1:
		return '1'
	else:
		return '0'

def checkedbox35():
	if chk35["state"] == "normal" and cb35.get() == 1:
		return '1'
	else:
		return '0'

def checkedbox36():
	if chk36["state"] == "normal" and cb36.get() == 1:
		return '1'
	else:
		return '0'

def checkedbox44():
	if chk44["state"] == "normal" and cb44.get() == 1:
		return '1'
	else:
		return '0'

def checkedbox45():
	if chk45["state"] == "normal" and cb45.get() == 1:
		return '1'
	else:
		return '0'

def checkedbox46():
	if chk46["state"] == "normal" and cb46.get() == 1:
		return '1'
	else:
		return '0'

def checkedbox55():
	if chk55["state"] == "normal" and cb55.get() == 1:
		return '1'
	else:
		return '0'

def checkedbox56():
	if chk56["state"] == "normal" and cb56.get() == 1:
		return '1'
	else:
		return '0'

def checkedbox66():
	if chk66["state"] == "normal" and cb66.get() == 1:
		return '1'
	else:
		return '0'

def fixdims():
	if dim.get() == 1:
		return '1'
	else:
		return '0'

def eulerval1():
	if eulervarpsi.get() == 1:
		return '1'
	else:
		return '0'

def eulerval2():
	if eulervartheta.get() == 1:
		return '1'
	else:
		return '0'

def eulerval3():
	if eulervarphi.get() == 1:
		return '1'
	else:
		return '0'

###########################################################################################
###########################################################################################
# Elastic Constants Frame inputs:

c11 = DoubleVar()
c11_bdry = IntVar()
c11_box = Entry(c_frame, width=5, textvariable=c11)
c11_box.grid(row=1, column=1, padx=10)
label11 = Label(c_frame, text= "c\u2081\u2081")
label11.grid(row=0, column=1)
c11_bd = Entry(c_frame, width=5, bg='peach puff',textvariable=c11_bdry)
c11_bd.grid(row=2, column=1)
labelc1n = Label(c_frame, text="Bounds:")
labelc1n.grid(row=2, column=0)
fixlabel = Label(c_frame, text="Fix:")
fixlabel.grid(row=3, column=0)
chk11 = Checkbutton(c_frame, variable=cb11, onvalue=1, offvalue=0, command=checkedbox11)
chk11.grid(row=3, column=1)

c12 = DoubleVar()
c12_bdry = IntVar()
c12_box = Entry(c_frame, width=5, textvariable=c12)
c12_box.grid(row=1, column=2, padx=10)
label12 = Label(c_frame, text= "c\u2081\u2082")
label12.grid(row=0, column=2)
c12_bd = Entry(c_frame, width=5, bg='peach puff',textvariable=c12_bdry)
c12_bd.grid(row=2, column=2)
chk12 = Checkbutton(c_frame, variable=cb12, onvalue=1, offvalue=0, command=checkedbox12)
chk12.grid(row=3, column=2)

c13_bdry = IntVar()
c13 = DoubleVar()
c13_box = Entry(c_frame, width=5, textvariable=c13)
c13_box.grid(row=1, column=3, padx=10)
label13 = Label(c_frame, text= "c\u2081\u2083")
label13.grid(row=0, column=3)
c13_bd = Entry(c_frame, width=5, bg='peach puff',textvariable=c13_bdry)
c13_bd.grid(row=2, column=3)
chk13 = Checkbutton(c_frame, variable=cb13, onvalue=1, offvalue=0, command=checkedbox13)
chk13.grid(row=3, column=3)

c14 = DoubleVar()
c14_bdry = IntVar()
c14_box = Entry(c_frame, width=5, textvariable=c14)
c14_box.grid(row=1, column=4, padx=10)
label14 = Label(c_frame, text= "c\u2081\u2084")
label14.grid(row=0, column=4)
c14_bd = Entry(c_frame, width=5, bg='peach puff',textvariable=c14_bdry)
c14_bd.grid(row=2, column=4)
chk14 = Checkbutton(c_frame, variable=cb14, onvalue=1, offvalue=0, command=checkedbox14)
chk14.grid(row=3, column=4)

c15 = DoubleVar()
c15_bdry = IntVar()
c15_box = Entry(c_frame, width=5, textvariable=c15)
c15_box.grid(row=1, column=5, padx=10)
label15 = Label(c_frame, text= "c\u2081\u2085")
label15.grid(row=0, column=5)
c15_bd = Entry(c_frame, width=5, bg='peach puff',textvariable=c15_bdry)
c15_bd.grid(row=2, column=5)
chk15 = Checkbutton(c_frame, variable=cb15, onvalue=1, offvalue=0, command=checkedbox15)
chk15.grid(row=3, column=5)

c16 = DoubleVar()
c16_bdry = IntVar()
c16_box = Entry(c_frame, width=5, textvariable=c16)
c16_box.grid(row=1, column=6, padx=10)
label16 = Label(c_frame, text= "c\u2081\u2086")
label16.grid(row=0, column=6)
c16_bd = Entry(c_frame, width=5, bg='peach puff',textvariable=c16_bdry)
c16_bd.grid(row=2, column=6)
chk16 = Checkbutton(c_frame, variable=cb16, onvalue=1, offvalue=0, command=checkedbox16)
chk16.grid(row=3, column=6)


c22 = DoubleVar()
c22_bdry = IntVar()
c22_box = Entry(c_frame, width=5, textvariable=c22)
c22_box.grid(row=5, column=2)
label22 = Label(c_frame, text= "c\u2082\u2082")
label22.grid(row=4, column=2)
c22_bd = Entry(c_frame, width=5, bg='peach puff',textvariable=c22_bdry)
c22_bd.grid(row=6, column=2)
labelc2n = Label(c_frame, text="Bounds:")
labelc2n.grid(row=6, column=1)
fixlabel2 = Label(c_frame, text="Fix:")
fixlabel2.grid(row=7, column=1)
chk22 = Checkbutton(c_frame, variable=cb22, onvalue=1, offvalue=0, command=checkedbox22)
chk22.grid(row=7, column=2)

c23 = DoubleVar()
c23_bdry = IntVar()
c23_box = Entry(c_frame, width=5, textvariable=c23)
c23_box.grid(row=5, column=3)
label23 = Label(c_frame, text= "c\u2082\u2083")
label23.grid(row=4, column=3)
c23_bd = Entry(c_frame, width=5, bg='peach puff',textvariable=c23_bdry)
c23_bd.grid(row=6, column=3)
chk23 = Checkbutton(c_frame, variable=cb23, onvalue=1, offvalue=0, command=checkedbox23)
chk23.grid(row=7, column=3)

c24 = DoubleVar()
c24_bdry = IntVar()
c24_box = Entry(c_frame, width=5, textvariable=c24)
c24_box.grid(row=5, column=4)
label24 = Label(c_frame, text= "c\u2082\u2084")
label24.grid(row=4, column=4)
c24_bd = Entry(c_frame, width=5, bg='peach puff',textvariable=c24_bdry)
c24_bd.grid(row=6, column=4)
chk24 = Checkbutton(c_frame, variable=cb24, onvalue=1, offvalue=0, command=checkedbox24)
chk24.grid(row=7, column=4)

c25 = DoubleVar()
c25_bdry = IntVar()
c25_box = Entry(c_frame, width=5, textvariable=c25)
c25_box.grid(row=5, column=5)
label25 = Label(c_frame, text= "c\u2082\u2085")
label25.grid(row=4, column=5)
c25_bd = Entry(c_frame, width=5, bg='peach puff',textvariable=c25_bdry)
c25_bd.grid(row=6, column=5)
chk25 = Checkbutton(c_frame, variable=cb25, onvalue=1, offvalue=0, command=checkedbox25)
chk25.grid(row=7, column=5)

c26 = DoubleVar()
c26_bdry = IntVar()
c26_box = Entry(c_frame, width=5, textvariable=c26)
c26_box.grid(row=5, column=6)
label26 = Label(c_frame, text= "c\u2082\u2086")
label26.grid(row=4, column=6)
c26_bd = Entry(c_frame, width=5, bg='peach puff',textvariable=c26_bdry)
c26_bd.grid(row=6, column=6)
chk26 = Checkbutton(c_frame, variable=cb26, onvalue=1, offvalue=0, command=checkedbox26)
chk26.grid(row=7, column=6)

c33 = DoubleVar()
c33_bdry = IntVar()
c33_box = Entry(c_frame, width=5, textvariable=c33)
c33_box.grid(row=9, column=3)
label33 = Label(c_frame, text= "c\u2083\u2083")
label33.grid(row=8, column=3)
c33_bd = Entry(c_frame, width=5, bg='peach puff',textvariable=c33_bdry)
c33_bd.grid(row=10, column=3)
labelc3n = Label(c_frame, text="Bounds:")
labelc3n.grid(row=10, column=2)
fixlabel3 = Label(c_frame, text="Fix:")
fixlabel3.grid(row=11, column=2)
chk33 = Checkbutton(c_frame, variable=cb33, onvalue=1, offvalue=0, command=checkedbox33)
chk33.grid(row=11, column=3)

c34 = DoubleVar()
c34_bdry = IntVar()
c34_box = Entry(c_frame, width=5, textvariable=c34)
c34_box.grid(row=9, column=4)
label34 = Label(c_frame, text= "c\u2083\u2084")
label34.grid(row=8, column=4)
c34_bd = Entry(c_frame, width=5, bg='peach puff',textvariable=c34_bdry)
c34_bd.grid(row=10, column=4)
chk34 = Checkbutton(c_frame, variable=cb34, onvalue=1, offvalue=0, command=checkedbox34)
chk34.grid(row=11, column=4)

c35 = DoubleVar()
c35_bdry = IntVar()
c35_box = Entry(c_frame, width=5, textvariable=c35)
c35_box.grid(row=9, column=5)
label35 = Label(c_frame, text= "c\u2083\u2085")
label35.grid(row=8, column=5)
c35_bd = Entry(c_frame, width=5, bg='peach puff',textvariable=c35_bdry)
c35_bd.grid(row=10, column=5)
chk35 = Checkbutton(c_frame, variable=cb35, onvalue=1, offvalue=0, command=checkedbox35)
chk35.grid(row=11, column=5)

c36 = DoubleVar()
c36_bdry = IntVar()
c36_box = Entry(c_frame, width=5, textvariable=c36)
c36_box.grid(row=9, column=6)
label36 = Label(c_frame, text= "c\u2083\u2086")
label36.grid(row=8, column=6)
c36_bd = Entry(c_frame, width=5, bg='peach puff',textvariable=c36_bdry)
c36_bd.grid(row=10, column=6)
chk36 = Checkbutton(c_frame, variable=cb36, onvalue=1, offvalue=0, command=checkedbox36)
chk36.grid(row=11, column=6)

c44 = DoubleVar()
c44_bdry = IntVar()
c44_box = Entry(c_frame, width=5, textvariable=c44)
c44_box.grid(row=13, column=4)
label44 = Label(c_frame, text= "c\u2084\u2084")
label44.grid(row=12, column=4)
c44_bd = Entry(c_frame, width=5, bg='peach puff',textvariable=c44_bdry)
c44_bd.grid(row=14, column=4)
labelc3n = Label(c_frame, text="Bounds:")
labelc3n.grid(row=14, column=3)
fixlabel4 = Label(c_frame, text="Fix:")
fixlabel4.grid(row=15, column=3)
chk44 = Checkbutton(c_frame, variable=cb44, onvalue=1, offvalue=0, command=checkedbox44)
chk44.grid(row=15, column=4)

c45 = DoubleVar()
c45_bdry = IntVar()
c45_box = Entry(c_frame, width=5, textvariable=c45)
c45_box.grid(row=13, column=5)
label45 = Label(c_frame, text= "c\u2084\u2085")
label45.grid(row=12, column=5)
c45_bd = Entry(c_frame, width=5, bg='peach puff',textvariable=c45_bdry)
c45_bd.grid(row=14, column=5)
chk45 = Checkbutton(c_frame, variable=cb45, onvalue=1, offvalue=0, command=checkedbox45)
chk45.grid(row=15, column=5)

c46 = DoubleVar()
c46_bdry = IntVar()
c46_box = Entry(c_frame, width=5, textvariable=c46)
c46_box.grid(row=13, column=6)
label46 = Label(c_frame, text= "c\u2084\u2086")
label46.grid(row=12, column=6)
c46_bd = Entry(c_frame, width=5, bg='peach puff',textvariable=c46_bdry)
c46_bd.grid(row=14, column=6)
chk46 = Checkbutton(c_frame, variable=cb46, onvalue=1, offvalue=0, command=checkedbox46)
chk46.grid(row=15, column=6)

c55 = DoubleVar()
c55_bdry = IntVar()
c55_box = Entry(c_frame, width=5, textvariable=c55)
c55_box.grid(row=17, column=5)
label55 = Label(c_frame, text= "c\u2085\u2085")
label55.grid(row=16, column=5)
c55_bd = Entry(c_frame, width=5, bg='peach puff',textvariable=c55_bdry)
c55_bd.grid(row=18, column=5)
labelc3n = Label(c_frame, text="Bounds:")
labelc3n.grid(row=18, column=4)
fixlabel5 = Label(c_frame, text="Fix:")
fixlabel5.grid(row=19, column=4)
chk55 = Checkbutton(c_frame, variable=cb55, onvalue=1, offvalue=0, command=checkedbox55)
chk55.grid(row=19, column=5)

c56 = DoubleVar()
c56_bdry = IntVar()
c56_box = Entry(c_frame, width=5, textvariable=c56)
c56_box.grid(row=17, column=6)
label56 = Label(c_frame, text= "c\u2085\u2086")
label56.grid(row=16, column=6)
c56_bd = Entry(c_frame, width=5, bg='peach puff',textvariable=c56_bdry)
c56_bd.grid(row=18, column=6)
chk56 = Checkbutton(c_frame, variable=cb56, onvalue=1, offvalue=0, command=checkedbox56)
chk56.grid(row=19, column=6)

c66 = DoubleVar()
c66_bdry = IntVar()
c66_box = Entry(c_frame, width=5, textvariable=c66)
c66_box.grid(row=21, column=6)
label66 = Label(c_frame, text= "c\u2086\u2086")
label66.grid(row=20, column=6)
c66_bd = Entry(c_frame, width=5, bg='peach puff',textvariable=c66_bdry)
c66_bd.grid(row=22, column=6)
labelc3n = Label(c_frame, text="Bounds:")
labelc3n.grid(row=22, column=5)
fixlabel6 = Label(c_frame, text="Fix:")
fixlabel6.grid(row=23, column=5)
chk66 = Checkbutton(c_frame, variable=cb66, onvalue=1, offvalue=0, command=checkedbox66)
chk66.grid(row=23, column=6)

###########################################################################################
###########################################################################################
# Dimensions Frame, added after defined fuction

d_frame = LabelFrame(second_frame, text = "Dimensions (cm units)", padx=25, pady=25)
d_frame.grid(row=12, rowspan=20, column=0, columnspan=5, sticky='nswe',padx=3)

# dimensions checkbox
dim = IntVar(value=1)
fix_dim = Checkbutton(d_frame, text="Fix Dimensions", variable=dim, onvalue=1, offvalue=0, command=fixdims)
fix_dim.grid(row=0, column=4, sticky=W+E)

file_frame = LabelFrame(second_frame, text = "Load Resonant Frequencies", padx=10, pady=10)
file_frame.grid(row=32, rowspan=1, column=0, columnspan=3, sticky='nswe', padx=3)

run_frame = LabelFrame(second_frame, text = "Run code", padx=10, pady=10)
run_frame.grid(row=32, rowspan=1, column=3, columnspan=2, sticky='nswe', padx=3)

run_ruscalbutton = Button(run_frame)
run_ruscalbutton.config(text= "Start ruscal", command = run_ruscal, bg='grey')
run_ruscalbutton.grid(row=0, column=3, padx=1, pady=8)



f_frame = LabelFrame(second_frame,text="Frequencies", padx=20, pady=20)
f_frame.grid(row=34, column=0, columnspan=15, rowspan=100, sticky='nswe', padx=3)


###########################################################################################
###########################################################################################
# Sample dimensions and density calculation:
dens = DoubleVar()
m = DoubleVar()
sup = str.maketrans("0123456789" , "⁰¹²³⁴⁵⁶⁷⁸⁹")
textd = "Density (g/cm3): "
labeld = Label(d_frame, text=textd.translate(sup))
labeld.grid(row=12, column=2)

# Cylinder entry boxes
d1 = DoubleVar()
d2 = DoubleVar()
d3 = DoubleVar()

entry1 = Entry(d_frame, width=10, textvariable=d1)
entry1.grid(row=0, column=1)
entry2 = Entry(d_frame, width=10, textvariable=d2)
entry2.grid(row=0, column=2)
entry3 = Entry(d_frame, width=10, textvariable=d3)
entry3.grid(row=0, column=3)
shapelabel = Label(d_frame, text="Shape: ")
shapelabel.grid(row=0, column=0, padx=0, sticky=W)
shape_label1 = Label(d_frame, text="")
shape_label1.grid(row=1, column=1)
shape_label2 = Label(d_frame, text="")
shape_label2.grid(row=1, column=2)
shape_label3 = Label(d_frame, text="")
shape_label3.grid(row=1, column=3)


# Potato entry boxes
px = DoubleVar()
py = DoubleVar()
pz = DoubleVar()
nx = DoubleVar()
ny = DoubleVar()
nz = DoubleVar()


potatox = Entry(d_frame, width=10, state='disabled', textvariable=px)
potatox.grid(row=8, column=1, padx=3, pady=3)
potatoy = Entry(d_frame, width=10, state='disabled', textvariable=py)
potatoy.grid(row=8, column=2, padx=3, pady=3)
potatoz = Entry(d_frame, width=10, state='disabled', textvariable=pz)
potatoz.grid(row=8, column=3, padx=3, pady=3)
potato_label = Label(d_frame, text="Potato: ")
potato_label.grid(row=8, column=0, padx=0, sticky=W)
x = Label(d_frame, text="x+")
x.grid(row=9, column=1)
y = Label(d_frame, text="y+")
y.grid(row=9, column=2)
z = Label(d_frame, text="z+")
z.grid(row=9, column=3)
negx = Entry(d_frame, width=10, state='disabled', textvariable=nx)
negx.grid(row=10, column=1)
negy = Entry(d_frame, width=10, state='disabled', textvariable=ny)
negy.grid(row=10, column=2)
negz = Entry(d_frame, width=10, state='disabled', textvariable=nz)
negz.grid(row=10, column=3)
negx_label = Label(d_frame, text="x-")
negx_label.grid(row=11, column=1)
negy_label = Label(d_frame, text="y-")
negy_label.grid(row=11, column=2)
negz_label = Label(d_frame, text="z-")
negz_label.grid(row=11, column=3)


weight = DoubleVar()
weightentry = Entry(d_frame, textvariable=weight, width=10)
weightentry.grid(row=12, column=1, padx=3, pady=10)
weightlabel = Label(d_frame, text="Mass (g):")
weightlabel.grid(row=12, column=0, padx=0, pady=10, sticky=W)


# Euler angles
euler1 = DoubleVar()
euler2 = DoubleVar()
euler3 = DoubleVar()
eulervarpsi = IntVar(value=1)
eulervartheta = IntVar(value=1)
eulervarphi = IntVar(value=1)
eulerlabel = Label(second_frame, text='Euler angles ', pady=5)
eulerlabel.grid(row=10, column=0, sticky=W, padx=10)
eulerb1 = Entry(second_frame, textvariable=euler1, width=10)
eulerb1.grid(row=10, column=1, sticky=W)
eulerb2 = Entry(second_frame, textvariable=euler2, width=10)
eulerb2.grid(row=10, column=2, sticky=W)
eulerb3 = Entry(second_frame, textvariable=euler3, width=10)
eulerb3.grid(row=10, column=3, sticky=W)
eulerl1 = Checkbutton(second_frame, text='psi', variable=eulervarpsi, onvalue=1, offvalue=0, command=eulerval1)
eulerl1.grid(row=11, column=1, sticky=W)
eulerl2 = Checkbutton(second_frame, text='theta', variable=eulervartheta, onvalue=1, offvalue=0, command=eulerval2)
eulerl2.grid(row=11, column=2, sticky=W)
eulerl3 = Checkbutton(second_frame, text='phi', variable=eulervarphi, onvalue=1, offvalue=0, command=eulerval3)
eulerl3.grid(row=11, column=3, sticky=W)

vol = DoubleVar()

def shapeshow(event):
	showLabel = Label(second_frame, text=shapeclicked.get())

	if shapeclicked.get() == shapes[0]: #cylinder
		entry1.config(state='normal')
		entry2.config(state='normal')
		entry3.config(state='normal')
		shapelabel.config(text="Cylinder:", state='normal')
		shape_label1.config(text="major axis", state='normal')
		shape_label2.config(text="minor axis", state='normal')
		shape_label3.config(text="  height", state='normal')
		potatox.config(state='disabled')
		potatoy.config(state='disabled')
		potatoz.config(state='disabled')
		x.config(state='disabled')
		y.config(state='disabled')
		z.config(state='disabled')
		negx.config(state='disabled')
		negy.config(state='disabled')
		negz.config(state='disabled')
		negx_label.config(state='disabled')
		negy_label.config(state='disabled')
		negz_label.config(state='disabled')
		vol = math.pi*(float(d1.get())/2)*(float(d2.get())/2*float(d3.get()))
		dens = float(weight.get()) / vol
		denst='{:07.4f}'.format(dens)
		labeldval = Label(d_frame, text='%s' % denst)
		labeldval.grid(row=12, column=3)

	if shapeclicked.get() == shapes[1]: #RPR
		entry1.config(state='normal')
		entry2.config(state='normal')
		entry3.config(state='normal')
		shapelabel.config(text="RPR:", state='normal')
		shape_label1.config(text="side 1 (a)", state='normal')
		shape_label2.config(text="side 2 (b)", state='normal')
		shape_label3.config(text="side 3 (c)", state='normal')
		potatox.config(state='disabled')
		potatoy.config(state='disabled')
		potatoz.config(state='disabled')
		x.config(state='disabled')
		y.config(state='disabled')
		z.config(state='disabled')
		negx.config(state='disabled')
		negy.config(state='disabled')
		negz.config(state='disabled')
		negx_label.config(state='disabled')
		negy_label.config(state='disabled')
		negz_label.config(state='disabled')
		vol = float(d1.get())*float(d2.get())*float(d3.get())
		dens = float(weight.get()) / vol
		denst='{:07.4f}'.format(dens)
		labeldval = Label(d_frame, text='%s' % denst)
		labeldval.grid(row=12, column=3)

	if shapeclicked.get() == shapes[2]: #Potato
		entry1.config(state='disabled')
		entry2.config(state='disabled')
		entry3.config(state='disabled')
		shapelabel.config(text="RPR:", state='disabled')
		shape_label1.config(text="side 1 (a)", state='disabled')
		shape_label2.config(text="side 2 (b)", state='disabled')
		shape_label3.config(text="side 3 (c)", state='disabled')
		potatox.config(state='normal')
		potatoy.config(state='normal')
		potatoz.config(state='normal')
		x.config(state='normal')
		y.config(state='normal')
		z.config(state='normal')
		negx.config(state='normal')
		negy.config(state='normal')
		negz.config(state='normal')
		negx_label.config(state='normal')
		negy_label.config(state='normal')
		negz_label.config(state='normal')
		vol = 1/6.0*math.pi*(float(potatox.get())+float(negx.get()))*(float(potatoy.get())+float(negy.get()))*(float(potatoz.get())+float(negz.get()))
		dens = float(weight.get()) / vol
		denst='{:07.4f}'.format(dens)
		labeldval = Label(d_frame, text='%s' % denst)
		labeldval.grid(row=12, column=3)

	if shapeclicked.get() == shapes[3]: #Ellipsoid
		entry1.config(state='normal')
		entry2.config(state='normal')
		entry3.config(state='normal')
		shapelabel.config(text="Ellipsoid:", state='normal')
		shape_label1.config(text="axis a", state='normal')
		shape_label2.config(text="axis b", state='normal')
		shape_label3.config(text="major axis", state='normal')
		potatox.config(state='disabled')
		potatoy.config(state='disabled')
		potatoz.config(state='disabled')
		x.config(state='disabled')
		y.config(state='disabled')
		z.config(state='disabled')
		negx.config(state='disabled')
		negy.config(state='disabled')
		negz.config(state='disabled')
		negx_label.config(state='disabled')
		negy_label.config(state='disabled')
		negz_label.config(state='disabled')
		vol = 4/3.0*math.pi*(float(d1.get())/2)*(float(d2.get())/2)*(float(d3.get())/2)
		dens = float(weight.get()) / vol
		denst='{:07.4f}'.format(dens)
		labeldval = Label(d_frame, text='%s' % denst)
		labeldval.grid(row=12, column=3)
		
	if shapeclicked.get() == shapes[4]: #Octahedron
		entry1.config(state='normal')
		entry2.config(state='normal')
		entry3.config(state='normal')
		shapelabel.config(text="Octahedron:", state='normal')
		shape_label1.config(text="axis a", state='normal')
		shape_label2.config(text="axis b", state='normal')
		shape_label3.config(text="axis c", state='normal')
		potatox.config(state='disabled')
		potatoy.config(state='disabled')
		potatoz.config(state='disabled')
		x.config(state='disabled')
		y.config(state='disabled')
		z.config(state='disabled')
		negx.config(state='disabled')
		negy.config(state='disabled')
		negz.config(state='disabled')
		negx_label.config(state='disabled')
		negy_label.config(state='disabled')
		negz_label.config(state='disabled')
		vol = 1/6.0*(float(d1.get()))*(float(d2.get()))*(float(d3.get()))
		dens = float(weight.get()) / vol
		denst='{:07.4f}'.format(dens)
		labeldval = Label(d_frame, text='%s' % denst)
		labeldval.grid(row=12, column=3)
		
	if shapeclicked.get() == shapes[5]: #Hollow Cylinder
		entry1.config(state='normal')
		entry2.config(state='normal')
		entry3.config(state='normal')
		shapelabel.config(text="Hollow Cyl.:", state='normal')
		shape_label1.config(text="O.D.", state='normal')
		shape_label2.config(text="I.D.", state='normal')
		shape_label3.config(text="  height", state='normal')
		potatox.config(state='disabled')
		potatoy.config(state='disabled')
		potatoz.config(state='disabled')
		x.config(state='disabled')
		y.config(state='disabled')
		z.config(state='disabled')
		negx.config(state='disabled')
		negy.config(state='disabled')
		negz.config(state='disabled')
		negx_label.config(state='disabled')
		negy_label.config(state='disabled')
		negz_label.config(state='disabled')
		vol = math.pi*float(d3.get())*((float(d1.get())/2)**2-(float(d2.get())/2)**2)
		dens = float(weight.get()) / vol
		denst='{:07.4f}'.format(dens)
		labeldval = Label(d_frame, text='%s' % denst)
		labeldval.grid(row=12, column=3)

shapetab = ['0','1','','','2','3','4','','5']
shapes = [
		"[0] Cylinder",
		"[1] RPR",
		"[4] Potato",
		"[5] Ellipsoid",
		"[6] Octahedron",
		"[8] Hollow Cylinder"
		]
def varyshow(event):
	showLabel = Label(second_frame, text=vary_clicked.get())

vary = [
		"[1] x (+)",
		"[2] y (+)",
		"[3] z (+)",
		"[4] x-",
		"[5] y-",
		"[6] z-",
		"[7] psi",
		"[8] theta",
		"[9] phi",
		]



def showbravais(event):
	bravaislabel = Label(second_frame, text=bravclass_clicked.get())

	if bravclass_clicked.get() == bravaisclass[0]: #isotropic
		c11_box.config(state='normal')
		label11.config(state='normal')
		c12_box.config(state='disabled')
		label12.config(state='disabled')
		c13_box.config(state='disabled')
		label13.config(state='disabled')
		c14_box.config(state='disabled')
		label14.config(state='disabled')
		c15_box.config(state='disabled')
		label15.config(state='disabled')
		c16_box.config(state='disabled')
		label16.config(state='disabled')
		c22_box.config(state='disabled')
		label22.config(state='disabled')
		c23_box.config(state='disabled')
		label23.config(state='disabled')
		c24_box.config(state='disabled')
		label24.config(state='disabled')
		c25_box.config(state='disabled')
		label25.config(state='disabled')
		c26_box.config(state='disabled')
		label26.config(state='disabled')
		c33_box.config(state='disabled')
		label33.config(state='disabled')
		c34_box.config(state='disabled')
		label34.config(state='disabled')
		c35_box.config(state='disabled')
		label35.config(state='disabled')
		c36_box.config(state='disabled')
		label36.config(state='disabled')
		c44_box.config(state='normal')
		label44.config(state='normal')
		c45_box.config(state='disabled')
		label45.config(state='disabled')
		c46_box.config(state='disabled')
		label46.config(state='disabled')
		c55_box.config(state='disabled')
		label55.config(state='disabled')
		c56_box.config(state='disabled')
		label56.config(state='disabled')
		c66_box.config(state='disabled')
		label66.config(state='disabled')
		c11_bd.config(state='normal')
		c12_bd.config(state='disabled')
		c13_bd.config(state='disabled')
		c14_bd.config(state='disabled')
		c15_bd.config(state='disabled')
		c16_bd.config(state='disabled')
		c22_bd.config(state='disabled')
		c23_bd.config(state='disabled')
		c24_bd.config(state='disabled')
		c25_bd.config(state='disabled')
		c26_bd.config(state='disabled')
		c33_bd.config(state='disabled')
		c34_bd.config(state='disabled')
		c35_bd.config(state='disabled')
		c36_bd.config(state='disabled')
		c44_bd.config(state='normal')
		c45_bd.config(state='disabled')
		c46_bd.config(state='disabled')
		c55_bd.config(state='disabled')
		c56_bd.config(state='disabled')
		c66_bd.config(state='disabled')
		chk11.config(state='normal')
		chk12.config(state='disabled')
		chk13.config(state='disabled')
		chk14.config(state='disabled')
		chk15.config(state='disabled')
		chk16.config(state='disabled')
		chk22.config(state='disabled')
		chk23.config(state='disabled')
		chk24.config(state='disabled')
		chk25.config(state='disabled')
		chk26.config(state='disabled')
		chk33.config(state='disabled')
		chk34.config(state='disabled')
		chk35.config(state='disabled')
		chk36.config(state='disabled')
		chk44.config(state='normal')
		chk45.config(state='disabled')
		chk46.config(state='disabled')
		chk55.config(state='disabled')
		chk56.config(state='disabled')
		chk66.config(state='disabled')

	if bravclass_clicked.get() == bravaisclass[1]: #cubic
		c11_box.config(state='normal')
		label11.config(state='normal')
		c12_box.config(state='normal')
		label12.config(state='normal')
		c13_box.config(state='disabled')
		label13.config(state='disabled')
		c14_box.config(state='disabled')
		label14.config(state='disabled')
		c15_box.config(state='disabled')
		label15.config(state='disabled')
		c16_box.config(state='disabled')
		label16.config(state='disabled')
		c22_box.config(state='disabled')
		label22.config(state='disabled')
		c23_box.config(state='disabled')
		label23.config(state='disabled')
		c24_box.config(state='disabled')
		label24.config(state='disabled')
		c25_box.config(state='disabled')
		label25.config(state='disabled')
		c26_box.config(state='disabled')
		label26.config(state='disabled')
		c33_box.config(state='disabled')
		label33.config(state='disabled')
		c34_box.config(state='disabled')
		label34.config(state='disabled')
		c35_box.config(state='disabled')
		label35.config(state='disabled')
		c36_box.config(state='disabled')
		label36.config(state='disabled')
		c44_box.config(state='normal')
		label44.config(state='normal')
		c45_box.config(state='disabled')
		label45.config(state='disabled')
		c46_box.config(state='disabled')
		label46.config(state='disabled')
		c55_box.config(state='disabled')
		label55.config(state='disabled')
		c56_box.config(state='disabled')
		label56.config(state='disabled')
		c66_box.config(state='disabled')
		label66.config(state='disabled')
		c11_bd.config(state='normal')
		c12_bd.config(state='normal')
		c13_bd.config(state='disabled')
		c14_bd.config(state='disabled')
		c15_bd.config(state='disabled')
		c16_bd.config(state='disabled')
		c22_bd.config(state='disabled')
		c23_bd.config(state='disabled')
		c24_bd.config(state='disabled')
		c25_bd.config(state='disabled')
		c26_bd.config(state='disabled')
		c33_bd.config(state='disabled')
		c34_bd.config(state='disabled')
		c35_bd.config(state='disabled')
		c36_bd.config(state='disabled')
		c44_bd.config(state='normal')
		c45_bd.config(state='disabled')
		c46_bd.config(state='disabled')
		c55_bd.config(state='disabled')
		c56_bd.config(state='disabled')
		c66_bd.config(state='disabled')
		chk11.config(state='normal')
		chk12.config(state='normal')
		chk13.config(state='disabled')
		chk14.config(state='disabled')
		chk15.config(state='disabled')
		chk16.config(state='disabled')
		chk22.config(state='disabled')
		chk23.config(state='disabled')
		chk24.config(state='disabled')
		chk25.config(state='disabled')
		chk26.config(state='disabled')
		chk33.config(state='disabled')
		chk34.config(state='disabled')
		chk35.config(state='disabled')
		chk36.config(state='disabled')
		chk44.config(state='normal')
		chk45.config(state='disabled')
		chk46.config(state='disabled')
		chk55.config(state='disabled')
		chk56.config(state='disabled')
		chk66.config(state='disabled')

	if bravclass_clicked.get() == bravaisclass[2]: #hexagonal
		c11_box.config(state='disabled')
		label11.config(state='disabled')
		c12_box.config(state='normal')
		label12.config(state='normal')
		c13_box.config(state='disabled')
		label13.config(state='disabled')
		c14_box.config(state='disabled')
		label14.config(state='disabled')
		c15_box.config(state='disabled')
		label15.config(state='disabled')
		c16_box.config(state='disabled')
		label16.config(state='disabled')
		c22_box.config(state='disabled')
		label22.config(state='disabled')
		c23_box.config(state='normal')
		label23.config(state='normal')
		c24_box.config(state='disabled')
		label24.config(state='disabled')
		c25_box.config(state='disabled')
		label25.config(state='disabled')
		c26_box.config(state='disabled')
		label26.config(state='disabled')
		c33_box.config(state='normal')
		label33.config(state='normal')
		c34_box.config(state='disabled')
		label34.config(state='disabled')
		c35_box.config(state='disabled')
		label35.config(state='disabled')
		c36_box.config(state='disabled')
		label36.config(state='disabled')
		c44_box.config(state='normal')
		label44.config(state='normal')
		c45_box.config(state='disabled')
		label45.config(state='disabled')
		c46_box.config(state='disabled')
		label46.config(state='disabled')
		c55_box.config(state='disabled')
		label55.config(state='disabled')
		c56_box.config(state='disabled')
		label56.config(state='disabled')
		c66_box.config(state='normal')
		label66.config(state='normal')
		c11_bd.config(state='disabled')
		c12_bd.config(state='normal')
		c13_bd.config(state='disabled')
		c14_bd.config(state='disabled')
		c15_bd.config(state='disabled')
		c16_bd.config(state='disabled')
		c22_bd.config(state='disabled')
		c23_bd.config(state='normal')
		c24_bd.config(state='disabled')
		c25_bd.config(state='disabled')
		c26_bd.config(state='disabled')
		c33_bd.config(state='normal')
		c34_bd.config(state='disabled')
		c35_bd.config(state='disabled')
		c36_bd.config(state='disabled')
		c44_bd.config(state='normal')
		c45_bd.config(state='disabled')
		c46_bd.config(state='disabled')
		c55_bd.config(state='disabled')
		c56_bd.config(state='disabled')
		c66_bd.config(state='normal')
		chk11.config(state='disabled')
		chk12.config(state='normal')
		chk13.config(state='disabled')
		chk14.config(state='disabled')
		chk15.config(state='disabled')
		chk16.config(state='disabled')
		chk22.config(state='disabled')
		chk23.config(state='normal')
		chk24.config(state='disabled')
		chk25.config(state='disabled')
		chk26.config(state='disabled')
		chk33.config(state='normal')
		chk34.config(state='disabled')
		chk35.config(state='disabled')
		chk36.config(state='disabled')
		chk44.config(state='normal')
		chk45.config(state='disabled')
		chk46.config(state='disabled')
		chk55.config(state='disabled')
		chk56.config(state='disabled')
		chk66.config(state='normal')

	if bravclass_clicked.get() == bravaisclass[3]: #trigonal 306
		c11_box.config(state='normal')
		label11.config(state='normal')
		c12_box.config(state='normal')
		label12.config(state='normal')
		c13_box.config(state='disabled')
		label13.config(state='disabled')
		c14_box.config(state='normal')
		label14.config(state='normal')
		c15_box.config(state='disabled')
		label15.config(state='disabled')
		c16_box.config(state='disabled')
		label16.config(state='disabled')
		c22_box.config(state='disabled')
		label22.config(state='disabled')
		c23_box.config(state='normal')
		label23.config(state='normal')
		c24_box.config(state='disabled')
		label24.config(state='disabled')
		c25_box.config(state='disabled')
		label25.config(state='disabled')
		c26_box.config(state='disabled')
		label26.config(state='disabled')
		c33_box.config(state='normal')
		label33.config(state='normal')
		c34_box.config(state='disabled')
		label34.config(state='disabled')
		c35_box.config(state='disabled')
		label35.config(state='disabled')
		c36_box.config(state='disabled')
		label36.config(state='disabled')
		c44_box.config(state='normal')
		label44.config(state='normal')
		c45_box.config(state='disabled')
		label45.config(state='disabled')
		c46_box.config(state='disabled')
		label46.config(state='disabled')
		c55_box.config(state='disabled')
		label55.config(state='disabled')
		c56_box.config(state='disabled')
		label56.config(state='disabled')
		c66_box.config(state='disabled')
		label66.config(state='disabled')
		c11_bd.config(state='normal')
		c12_bd.config(state='normal')
		c13_bd.config(state='disabled')
		c14_bd.config(state='normal')
		c15_bd.config(state='disabled')
		c16_bd.config(state='disabled')
		c22_bd.config(state='disabled')
		c23_bd.config(state='normal')
		c24_bd.config(state='disabled')
		c25_bd.config(state='disabled')
		c26_bd.config(state='disabled')
		c33_bd.config(state='normal')
		c34_bd.config(state='disabled')
		c35_bd.config(state='disabled')
		c36_bd.config(state='disabled')
		c44_bd.config(state='normal')
		c45_bd.config(state='disabled')
		c46_bd.config(state='disabled')
		c55_bd.config(state='disabled')
		c56_bd.config(state='disabled')
		c66_bd.config(state='disabled')
		chk11.config(state='normal')
		chk12.config(state='normal')
		chk13.config(state='disabled')
		chk14.config(state='normal')
		chk15.config(state='disabled')
		chk16.config(state='disabled')
		chk22.config(state='disabled')
		chk23.config(state='normal')
		chk24.config(state='disabled')
		chk25.config(state='disabled')
		chk26.config(state='disabled')
		chk33.config(state='normal')
		chk34.config(state='disabled')
		chk35.config(state='disabled')
		chk36.config(state='disabled')
		chk44.config(state='normal')
		chk45.config(state='disabled')
		chk46.config(state='disabled')
		chk55.config(state='disabled')
		chk56.config(state='disabled')
		chk66.config(state='disabled')

	if bravclass_clicked.get() == bravaisclass[4]: #trigonal 307
		c11_box.config(state='normal')
		label11.config(state='normal')
		c12_box.config(state='normal')
		label12.config(state='normal')
		c13_box.config(state='normal')
		label13.config(state='normal')
		c14_box.config(state='normal')
		label14.config(state='normal')
		c15_box.config(state='disabled')
		label15.config(state='disabled')
		c16_box.config(state='disabled')
		label16.config(state='disabled')
		c22_box.config(state='disabled')
		label22.config(state='disabled')
		c23_box.config(state='disabled')
		label23.config(state='disabled')
		c24_box.config(state='disabled')
		label24.config(state='disabled')
		c25_box.config(state='normal')
		label25.config(state='normal')
		c26_box.config(state='disabled')
		label26.config(state='disabled')
		c33_box.config(state='normal')
		label33.config(state='normal')
		c34_box.config(state='disabled')
		label34.config(state='disabled')
		c35_box.config(state='disabled')
		label35.config(state='disabled')
		c36_box.config(state='disabled')
		label36.config(state='disabled')
		c44_box.config(state='normal')
		label44.config(state='normal')
		c45_box.config(state='disabled')
		label45.config(state='disabled')
		c46_box.config(state='disabled')
		label46.config(state='disabled')
		c55_box.config(state='disabled')
		label55.config(state='disabled')
		c56_box.config(state='disabled')
		label56.config(state='disabled')
		c66_box.config(state='disabled')
		label66.config(state='disabled')
		c11_bd.config(state='normal')
		c12_bd.config(state='normal')
		c13_bd.config(state='normal')
		c14_bd.config(state='normal')
		c15_bd.config(state='disabled')
		c16_bd.config(state='disabled')
		c22_bd.config(state='disabled')
		c23_bd.config(state='disabled')
		c24_bd.config(state='disabled')
		c25_bd.config(state='normal')
		c26_bd.config(state='disabled')
		c33_bd.config(state='normal')
		c34_bd.config(state='disabled')
		c35_bd.config(state='disabled')
		c36_bd.config(state='disabled')
		c44_bd.config(state='normal')
		c45_bd.config(state='disabled')
		c46_bd.config(state='disabled')
		c55_bd.config(state='disabled')
		c56_bd.config(state='disabled')
		c66_bd.config(state='disabled')
		chk11.config(state='normal')
		chk12.config(state='normal')
		chk13.config(state='normal')
		chk14.config(state='normal')
		chk15.config(state='disabled')
		chk16.config(state='disabled')
		chk22.config(state='disabled')
		chk23.config(state='disabled')
		chk24.config(state='disabled')
		chk25.config(state='normal')
		chk26.config(state='disabled')
		chk33.config(state='normal')
		chk34.config(state='disabled')
		chk35.config(state='disabled')
		chk36.config(state='disabled')
		chk44.config(state='normal')
		chk45.config(state='disabled')
		chk46.config(state='disabled')
		chk55.config(state='disabled')
		chk56.config(state='disabled')
		chk66.config(state='disabled')

	if bravclass_clicked.get() == bravaisclass[5]: #tetragonal 406
		c11_box.config(state='normal')
		label11.config(state='normal')
		c12_box.config(state='normal')
		label12.config(state='normal')
		c13_box.config(state='disabled')
		label13.config(state='disabled')
		c14_box.config(state='disabled')
		label14.config(state='disabled')
		c15_box.config(state='disabled')
		label15.config(state='disabled')
		c16_box.config(state='disabled')
		label16.config(state='disabled')
		c22_box.config(state='disabled')
		label22.config(state='disabled')
		c23_box.config(state='normal')
		label23.config(state='normal')
		c24_box.config(state='disabled')
		label24.config(state='disabled')
		c25_box.config(state='disabled')
		label25.config(state='disabled')
		c26_box.config(state='disabled')
		label26.config(state='disabled')
		c33_box.config(state='normal')
		label33.config(state='normal')
		c34_box.config(state='disabled')
		label34.config(state='disabled')
		c35_box.config(state='disabled')
		label35.config(state='disabled')
		c36_box.config(state='disabled')
		label36.config(state='disabled')
		c44_box.config(state='normal')
		label44.config(state='normal')
		c45_box.config(state='disabled')
		label45.config(state='disabled')
		c46_box.config(state='disabled')
		label46.config(state='disabled')
		c55_box.config(state='disabled')
		label55.config(state='disabled')
		c56_box.config(state='disabled')
		label56.config(state='disabled')
		c66_box.config(state='normal')
		label66.config(state='normal')
		c11_bd.config(state='normal')
		c12_bd.config(state='normal')
		c13_bd.config(state='disabled')
		c14_bd.config(state='disabled')
		c15_bd.config(state='disabled')
		c16_bd.config(state='disabled')
		c22_bd.config(state='disabled')
		c23_bd.config(state='normal')
		c24_bd.config(state='disabled')
		c25_bd.config(state='disabled')
		c26_bd.config(state='disabled')
		c33_bd.config(state='normal')
		c34_bd.config(state='disabled')
		c35_bd.config(state='disabled')
		c36_bd.config(state='disabled')
		c44_bd.config(state='normal')
		c45_bd.config(state='disabled')
		c46_bd.config(state='disabled')
		c55_bd.config(state='disabled')
		c56_bd.config(state='disabled')
		c66_bd.config(state='normal')
		chk11.config(state='normal')
		chk12.config(state='normal')
		chk13.config(state='disabled')
		chk14.config(state='disabled')
		chk15.config(state='disabled')
		chk16.config(state='disabled')
		chk22.config(state='disabled')
		chk23.config(state='normal')
		chk24.config(state='disabled')
		chk25.config(state='disabled')
		chk26.config(state='disabled')
		chk33.config(state='normal')
		chk34.config(state='disabled')
		chk35.config(state='disabled')
		chk36.config(state='disabled')
		chk44.config(state='normal')
		chk45.config(state='disabled')
		chk46.config(state='disabled')
		chk55.config(state='disabled')
		chk56.config(state='disabled')
		chk66.config(state='normal')

	if bravclass_clicked.get() == bravaisclass[6]: #tetragonal 407
		c11_box.config(state='normal')
		label11.config(state='normal')
		c12_box.config(state='normal')
		label12.config(state='normal')
		c13_box.config(state='disabled')
		label13.config(state='disabled')
		c14_box.config(state='disabled')
		label14.config(state='disabled')
		c15_box.config(state='disabled')
		label15.config(state='disabled')
		c16_box.config(state='normal')
		label16.config(state='normal')
		c22_box.config(state='disabled')
		label22.config(state='disabled')
		c23_box.config(state='normal')
		label23.config(state='normal')
		c24_box.config(state='disabled')
		label24.config(state='disabled')
		c25_box.config(state='disabled')
		label25.config(state='disabled')
		c26_box.config(state='disabled')
		label26.config(state='disabled')
		c33_box.config(state='normal')
		label33.config(state='normal')
		c34_box.config(state='disabled')
		label34.config(state='disabled')
		c35_box.config(state='disabled')
		label35.config(state='disabled')
		c36_box.config(state='disabled')
		label36.config(state='disabled')
		c44_box.config(state='normal')
		label44.config(state='normal')
		c45_box.config(state='disabled')
		label45.config(state='disabled')
		c46_box.config(state='disabled')
		label46.config(state='disabled')
		c55_box.config(state='disabled')
		label55.config(state='disabled')
		c56_box.config(state='disabled')
		label56.config(state='disabled')
		c66_box.config(state='normal')
		label66.config(state='normal')
		c11_bd.config(state='normal')
		c12_bd.config(state='normal')
		c13_bd.config(state='disabled')
		c14_bd.config(state='disabled')
		c15_bd.config(state='disabled')
		c16_bd.config(state='normal')
		c22_bd.config(state='disabled')
		c23_bd.config(state='normal')
		c24_bd.config(state='disabled')
		c25_bd.config(state='disabled')
		c26_bd.config(state='disabled')
		c33_bd.config(state='normal')
		c34_bd.config(state='disabled')
		c35_bd.config(state='disabled')
		c36_bd.config(state='disabled')
		c44_bd.config(state='normal')
		c45_bd.config(state='disabled')
		c46_bd.config(state='disabled')
		c55_bd.config(state='disabled')
		c56_bd.config(state='disabled')
		c66_bd.config(state='normal')
		chk11.config(state='normal')
		chk12.config(state='normal')
		chk13.config(state='disabled')
		chk14.config(state='disabled')
		chk15.config(state='disabled')
		chk16.config(state='normal')
		chk22.config(state='disabled')
		chk23.config(state='normal')
		chk24.config(state='disabled')
		chk25.config(state='disabled')
		chk26.config(state='disabled')
		chk33.config(state='normal')
		chk34.config(state='disabled')
		chk35.config(state='disabled')
		chk36.config(state='disabled')
		chk44.config(state='normal')
		chk45.config(state='disabled')
		chk46.config(state='disabled')
		chk55.config(state='disabled')
		chk56.config(state='disabled')
		chk66.config(state='normal')

	if bravclass_clicked.get() == bravaisclass[7]: #orthorhombic
		c11_box.config(state='normal')
		label11.config(state='normal')
		c12_box.config(state='normal')
		label12.config(state='normal')
		c13_box.config(state='normal')
		label13.config(state='normal')
		c14_box.config(state='disabled')
		label14.config(state='disabled')
		c15_box.config(state='disabled')
		label15.config(state='disabled')
		c16_box.config(state='disabled')
		label16.config(state='disabled')
		c22_box.config(state='normal')
		label22.config(state='normal')
		c23_box.config(state='normal')
		label23.config(state='normal')
		c24_box.config(state='disabled')
		label24.config(state='disabled')
		c25_box.config(state='disabled')
		label25.config(state='disabled')
		c26_box.config(state='disabled')
		label26.config(state='disabled')
		c33_box.config(state='normal')
		label33.config(state='normal')
		c34_box.config(state='disabled')
		label34.config(state='disabled')
		c35_box.config(state='disabled')
		label35.config(state='disabled')
		c36_box.config(state='disabled')
		label36.config(state='disabled')
		c44_box.config(state='normal')
		label44.config(state='normal')
		c45_box.config(state='disabled')
		label45.config(state='disabled')
		c46_box.config(state='disabled')
		label46.config(state='disabled')
		c55_box.config(state='normal')
		label55.config(state='normal')
		c56_box.config(state='disabled')
		label56.config(state='disabled')
		c66_box.config(state='normal')
		label66.config(state='normal')
		c11_bd.config(state='normal')
		c12_bd.config(state='normal')
		c13_bd.config(state='normal')
		c14_bd.config(state='disabled')
		c15_bd.config(state='disabled')
		c16_bd.config(state='disabled')
		c22_bd.config(state='normal')
		c23_bd.config(state='normal')
		c24_bd.config(state='disabled')
		c25_bd.config(state='disabled')
		c26_bd.config(state='disabled')
		c33_bd.config(state='normal')
		c34_bd.config(state='disabled')
		c35_bd.config(state='disabled')
		c36_bd.config(state='disabled')
		c44_bd.config(state='normal')
		c45_bd.config(state='disabled')
		c46_bd.config(state='disabled')
		c55_bd.config(state='normal')
		c56_bd.config(state='disabled')
		c66_bd.config(state='normal')
		chk11.config(state='normal')
		chk12.config(state='normal')
		chk13.config(state='normal')
		chk14.config(state='disabled')
		chk15.config(state='disabled')
		chk16.config(state='disabled')
		chk22.config(state='normal')
		chk23.config(state='normal')
		chk24.config(state='disabled')
		chk25.config(state='disabled')
		chk26.config(state='disabled')
		chk33.config(state='normal')
		chk34.config(state='disabled')
		chk35.config(state='disabled')
		chk36.config(state='disabled')
		chk44.config(state='normal')
		chk45.config(state='disabled')
		chk46.config(state='disabled')
		chk55.config(state='normal')
		chk56.config(state='disabled')
		chk66.config(state='normal')

	if bravclass_clicked.get() == bravaisclass[8]: #monoclinic
		c11_box.config(state='normal')
		label11.config(state='normal')
		c12_box.config(state='normal')
		label12.config(state='normal')
		c13_box.config(state='normal')
		label13.config(state='normal')
		c14_box.config(state='disabled')
		label14.config(state='disabled')
		c15_box.config(state='normal')
		label15.config(state='normal')
		c16_box.config(state='disabled')
		label16.config(state='disabled')
		c22_box.config(state='normal')
		label22.config(state='normal')
		c23_box.config(state='normal')
		label23.config(state='normal')
		c24_box.config(state='disabled')
		label24.config(state='disabled')
		c25_box.config(state='normal')
		label25.config(state='normal')
		c26_box.config(state='disabled')
		label26.config(state='disabled')
		c33_box.config(state='normal')
		label33.config(state='normal')
		c34_box.config(state='disabled')
		label34.config(state='disabled')
		c35_box.config(state='normal')
		label35.config(state='normal')
		c36_box.config(state='disabled')
		label36.config(state='disabled')
		c44_box.config(state='normal')
		label44.config(state='normal')
		c45_box.config(state='disabled')
		label45.config(state='disabled')
		c46_box.config(state='normal')
		label46.config(state='normal')
		c55_box.config(state='normal')
		label55.config(state='normal')
		c56_box.config(state='disabled')
		label56.config(state='disabled')
		c66_box.config(state='normal')
		label66.config(state='normal')
		c11_bd.config(state='normal')
		c12_bd.config(state='normal')
		c13_bd.config(state='normal')
		c14_bd.config(state='disabled')
		c15_bd.config(state='normal')
		c16_bd.config(state='disabled')
		c22_bd.config(state='normal')
		c23_bd.config(state='normal')
		c24_bd.config(state='disabled')
		c25_bd.config(state='normal')
		c26_bd.config(state='disabled')
		c33_bd.config(state='normal')
		c34_bd.config(state='disabled')
		c35_bd.config(state='normal')
		c36_bd.config(state='disabled')
		c44_bd.config(state='normal')
		c45_bd.config(state='disabled')
		c46_bd.config(state='normal')
		c55_bd.config(state='normal')
		c56_bd.config(state='disabled')
		c66_bd.config(state='normal')
		chk11.config(state='normal')
		chk12.config(state='normal')
		chk13.config(state='normal')
		chk14.config(state='disabled')
		chk15.config(state='normal')
		chk16.config(state='disabled')
		chk22.config(state='normal')
		chk23.config(state='normal')
		chk24.config(state='disabled')
		chk25.config(state='normal')
		chk26.config(state='disabled')
		chk33.config(state='normal')
		chk34.config(state='disabled')
		chk35.config(state='normal')
		chk36.config(state='disabled')
		chk44.config(state='normal')
		chk45.config(state='disabled')
		chk46.config(state='normal')
		chk55.config(state='normal')
		chk56.config(state='disabled')
		chk66.config(state='normal')

	if bravclass_clicked.get() == bravaisclass[9]: #triclinic
		c11_box.config(state='normal')
		label11.config(state='normal')
		c12_box.config(state='normal')
		label12.config(state='normal')
		c13_box.config(state='normal')
		label13.config(state='normal')
		c14_box.config(state='normal')
		label14.config(state='normal')
		c15_box.config(state='normal')
		label15.config(state='normal')
		c16_box.config(state='normal')
		label16.config(state='normal')
		c22_box.config(state='normal')
		label22.config(state='normal')
		c23_box.config(state='normal')
		label23.config(state='normal')
		c24_box.config(state='normal')
		label24.config(state='normal')
		c25_box.config(state='normal')
		label25.config(state='normal')
		c26_box.config(state='normal')
		label26.config(state='normal')
		c33_box.config(state='normal')
		label33.config(state='normal')
		c34_box.config(state='normal')
		label34.config(state='normal')
		c35_box.config(state='normal')
		label35.config(state='normal')
		c36_box.config(state='normal')
		label36.config(state='normal')
		c44_box.config(state='normal')
		label44.config(state='normal')
		c45_box.config(state='normal')
		label45.config(state='normal')
		c46_box.config(state='normal')
		label46.config(state='normal')
		c55_box.config(state='normal')
		label55.config(state='normal')
		c56_box.config(state='normal')
		label56.config(state='normal')
		c66_box.config(state='normal')
		label66.config(state='normal')
		c11_bd.config(state='normal')
		c12_bd.config(state='normal')
		c13_bd.config(state='normal')
		c14_bd.config(state='normal')
		c15_bd.config(state='normal')
		c16_bd.config(state='normal')
		c22_bd.config(state='normal')
		c23_bd.config(state='normal')
		c24_bd.config(state='normal')
		c25_bd.config(state='normal')
		c26_bd.config(state='normal')
		c33_bd.config(state='normal')
		c34_bd.config(state='normal')
		c35_bd.config(state='normal')
		c36_bd.config(state='normal')
		c44_bd.config(state='normal')
		c45_bd.config(state='normal')
		c46_bd.config(state='normal')
		c55_bd.config(state='normal')
		c56_bd.config(state='normal')
		c66_bd.config(state='normal')
		chk11.config(state='normal')
		chk12.config(state='normal')
		chk13.config(state='normal')
		chk14.config(state='normal')
		chk15.config(state='normal')
		chk16.config(state='normal')
		chk22.config(state='normal')
		chk23.config(state='normal')
		chk24.config(state='normal')
		chk25.config(state='normal')
		chk26.config(state='normal')
		chk33.config(state='normal')
		chk34.config(state='normal')
		chk35.config(state='normal')
		chk36.config(state='normal')
		chk44.config(state='normal')
		chk45.config(state='normal')
		chk46.config(state='normal')
		chk55.config(state='normal')
		chk56.config(state='normal')
		chk66.config(state='normal')
	
###########################################################################################
###########################################################################################
'''
Loading frequencies

Need to have seperate functions for excel, csv, and text files
Will do a drop menu that has different options:
- List boxes for text files and csv files
- Treeview for  excel
'''

file_types = ["CSV File",
			  "Excel File",
			  "Text File"
			  ]

def load_freqs(event):
	loadingfreqs = Label(file_frame, text=freqs_file.get())

	def clear_data():
		treef.delete(*treef.get_children())

	def remove_row():
		x = treef.selection()
		for i in x:
			treef.delete(i)

	def select_values():
		# Clear out anything in entry boxes
		F_ex_entry.delete(0, END)
		F_th_entry.delete(0, END)
		Weight_entry.delete(0, END)

		# Grabbing row number
		selected = treef.focus()
		# Grabbing values
		values = treef.item(selected, 'values')
		# Output to entry boxes
		F_ex_entry.insert(0, values[1])
		F_th_entry.insert(0, values[2])
		Weight_entry.insert(0, float(values[3]))

	def update():
		selected = treef.focus()
		# Save new data
		values = treef.item(selected, 'values')
		idx=values[0]
		treef.item(selected, text="", values=(idx,F_ex_entry.get(), F_th_entry.get(), Weight_entry.get()))

		# Clear out anything in entry boxes
#		F_ex_entry.delete(0, END)
#		F_th_entry.delete(0, END)
#		Weight_entry.delete(0, END)

	if freqs_file.get() == file_types[0]:

		clear_data()
		def open_file():
			csv_file = filedialog.askopenfilename(initialdir=os.path.abspath(os.getcwd()),
												  title="Select A File",
												  filetypes=(("csv files", "*.csv"), ("All Files", "*.*")))
			label_file.config(text="")
			label_file["text"] = csv_file

		def load_data():
			clear_data()
			with open(label_file["text"]) as f:
				next(f) #skipping header line
				reader = csv.reader(f, delimiter=',')
				i=1
				for j in reader:
					treef.insert("",'end', values=(i,j[0],j[1],float(j[2])))
					i=i+1

	###########################################################################################	
	# For Excel Files
	if freqs_file.get() == file_types[1]:

		clear_data()
		def open_file():
			filename = filedialog.askopenfilename(initialdir=os.path.abspath(os.getcwd()), 
												  title="Select A File", 
												  filetypes=(("xlsx files", "*.xlsx"), ("All Files","*.*")))
			label_file.config(text="")
			label_file["text"] = filename


		def load_data():
			file_path = label_file["text"]
			try:
				excel_filename = r"{}".format(file_path)
				df = pd.read_excel(excel_filename)
			except ValueError:
				label_file.config(text="File Could Not be Opened")
			except FileNotFoundError:
				label_file.config(text="File Could Not be Found")

			clear_data()
#			treef["column"] = list(df.columns)
			treef["show"] = "headings"

#			for column in treef["columns"]:
#				treef.heading(column, text=column)

			df_rows = df.to_numpy().tolist()
			i=1
			for row in df_rows:
				treef.insert("", "end", values=(i,row[0],row[1],float(row[2])))
				i=i+1

	###########################################################################################	
	# For text files
	if freqs_file.get() == file_types[2]:
		clear_data()
		
		def open_file():
			textfile = filedialog.askopenfilename(initialdir=os.path.abspath(os.getcwd()),
												  title="Select A File",
												  filetypes=(("text files", "*.txt"), ("All Files", "*.*")))
			label_file.config(text="")
			label_file["text"] = textfile

		def load_data():
			clear_data()
			with open(label_file["text"]) as f:
				lines = f.readlines()[1:]
				nb = []
				f_ex = []
				f_th = []
				weights = []
				i = 1
				for x in lines:
					data = x.split()  #account for spaces and tabs
					treef.insert("", "end", values=(i,data[0], data[1], float(data[2])))
					i=i+1

	###########################################################################################	
	
	label_file = Label(file_frame)
	label_file.grid(row=1, column=0,columnspan=5)

	if(freqs_file.get()!="Choose File Type"):
		file_button = Button(file_frame, text="Choose File", command=open_file)
		file_button.grid(row=0, column=1)

		load_button = Button(file_frame, text="Load File", command=load_data)
		load_button.grid(row=0, column=2)


	# Buttons in f_fame
	select_button = Button(f_frame, text="Select", command=select_values)
	select_button.grid(row=1, column=10, sticky="ew", padx=5)

	update_button = Button(f_frame, text="Save/Update", command=update)
	update_button.grid(row=1, column=11, sticky="ew", padx=5)

	delete_button = Button(f_frame, text="Delete Row", command=remove_row)
	delete_button.grid(row=1, column=12, sticky="ew", padx=5)

#	num_label = Label(f_frame, text=" # ")
#	num_label.grid(row=2, column=9)
#	num_entry = Entry(f_frame)
#	num_entry.grid(row=3, column=9)

	ex_label = Label(f_frame, text="F_ex(MHz)")
	ex_label.grid(row=2, column=10)
	F_ex_entry = Entry(f_frame)
	F_ex_entry.grid(row=3, column=10)

	th_label = Label(f_frame, text="F_th(MHz)")
	th_label.grid(row=2, column=11)
	F_th_entry = Entry(f_frame)
	F_th_entry.grid(row=3, column=11)

	w_label = Label(f_frame, text="Weight")
	w_label.grid(row=2, column=12)
	Weight_entry = Entry(f_frame)
	Weight_entry.grid(row=3, column=12)

	temp_label = Label(f_frame, text='')
	temp_label.grid(row=7, column=10)

# Treeview Widget
treef = ttk.Treeview(f_frame, columns=("number","F_ex(MHz)", "F_th(MHz)", "Weight"), height=30, selectmode="extended")
treef.column("number",width=60)
treef.column("F_ex(MHz)",width=120)
treef.column("F_th(MHz)",width=120)
treef.column("Weight",width=60)
tree_scrolly = ttk.Scrollbar(f_frame, orient=VERTICAL, command=treef.yview)
tree_scrollx = ttk.Scrollbar(f_frame, orient=HORIZONTAL, command=treef.xview)
treef.grid(column=0, columnspan=1, row=0, rowspan=50)
treef.configure(xscrollcommand=tree_scrollx.set, yscrollcommand=tree_scrolly.set)

treef["show"] = "headings"
treef.heading('number', text="Number", anchor=W)
treef.heading('F_ex(MHz)', text="F_ex(MHz)", anchor=W)
treef.heading('F_th(MHz)', text="F_th(MHz)", anchor=W)
treef.heading('Weight', text="Weight", anchor=W)
###########################################################################################	

###########################################################################################	

def tabbrav(i):
	if(i==2):
		return 0;
	if(i==3):
		return 1;
	if(i==5):
		return 2;
	if(i==306):
		return 3;
	if(i==307):
		return 4;
	if(i==406):
		return 5;
	if(i==407):
		return 6;
	if(i==9):
		return 7;
	if(i==13):
		return 8;
	if(i==21):
		return 9;
	return 999;



bravaisclass = [
				"[2] isotropic: all classes",
				"[3] cubic: all classes",
				"[5] hexagonal: all classes",
				"[306] trigonal: 32, -3m, 3m classes",
				"[307] trigonal: 3, -3 classes",
				"[406] tetragonal: 4mm, -42mm, 422, 4/mm classes",
				"[407] tetragonal: 4, -4, 4/m classes",
				"[9] orthohombic: all classes",
				"[13] monoclinic: all classes",
				"[21] triclinic: both classes"
				]


def MPolyShow():
	maxpolylabel = Label(second_frame, text=maxpoly_clicked.get())

maxpolyorder = [
				"8",
				"9",
				"10",
				"11",
				"12",
				"13",
				"14"
				]

def mirrorshow():
	mirrorlabel = Label(second_frame, text=mirror_clicked.get())

mirrorplanes = [
				"0 = No mirror plane",
				"1 = Perpendicular to z",
				"2 = Perpendicular to z and x",
				"3 = All planes"
				]

shapeclicked = StringVar()
shapeclicked.set("Select Shape")

bravclass_clicked = StringVar()
bravclass_clicked.set("Select Class")

vary_clicked = StringVar()
vary_clicked.set("Vary Dim.")

maxpoly_clicked = StringVar()
maxpoly_clicked.set(maxpolyorder[0])

mirror_clicked = StringVar()
mirror_clicked.set(mirrorplanes[0])

freqs_file = StringVar()
freqs_file.set("Choose File Type")

weight_clicked = DoubleVar()

###########################################################################################
###########################################################################################
# Drop Down Menus:
file_drop = OptionMenu(file_frame, freqs_file, *file_types, command=load_freqs)
file_drop.grid(row=0, column=0, sticky=W)

shape_drop = OptionMenu(second_frame, shapeclicked, *shapes, command=shapeshow)
shape_drop.grid(row=2, column=1, sticky=W)

bravais_drop = OptionMenu(second_frame, bravclass_clicked, *bravaisclass, command=showbravais)
#bravais_drop.bind('<Button-1>', showbravais)
bravais_drop.grid(row=3, column=1, columnspan=4, sticky=W)

maxpoly_drop = OptionMenu(second_frame, maxpoly_clicked, *maxpolyorder)
maxpoly_drop.grid(row=5, column=1, sticky=W)

mirror_drop = OptionMenu(second_frame, mirror_clicked, *mirrorplanes)
mirror_drop.grid(row=9, column=1, columnspan=2, sticky=W)

###########################################################################################
###########################################################################################
# Radio Buttons for selection of calculation mode:
calcmodevar = IntVar()
calcderiv_label = Label(second_frame, text="Calculate derivatives")
calcderiv_label.grid(row=8, column=0, padx=10, sticky=W)
calcderiv = Checkbutton(second_frame, variable=calcmodevar, onvalue=1, offvalue=0)
calcderiv.grid(row=8, column=1, sticky=W)
vary_drop = OptionMenu(second_frame, vary_clicked, *vary, command=varyshow)
vary_drop.grid(row=8, column=2, sticky=W)
vary_drop.configure(state="disabled")

def calcentry():
	calcmode_n.config(state='normal')

def nocalcentry():
	calcmode_n.config(state='disabled')

def clear():
	chk11.deselect()
	chk12.deselect()
	chk13.deselect()
	chk14.deselect()
	chk15.deselect()
	chk16.deselect()
	chk22.deselect()
	chk23.deselect()
	chk24.deselect()
	chk25.deselect()
	chk26.deselect()
	chk33.deselect()
	chk34.deselect()
	chk35.deselect()
	chk36.deselect()
	chk44.deselect()
	chk45.deselect()
	chk46.deselect()
	chk55.deselect()
	chk56.deselect()
	chk66.deselect()

def clearentry():
	calcmode_n.delete(0, END)

modes = IntVar()
modes.set(1)

def calcmodefunc():
	if modes.get() == 1:
		vary_drop.configure(state="disabled")
		calcval = calcmode_n.get()
		calcentry()
		return calcval
    
	elif modes.get() == -2:
		clearentry()
		clear()
		nocalcentry()
		vary_drop.configure(state="active")
		return '-2'

	elif modes.get() == -1:
		vary_drop.configure(state="disabled")
		clearentry()
		clear()
		nocalcentry()
		return '-1'

	elif modes.get() == -3:
		vary_drop.configure(state="disabled")
		clearentry()
		nocalcentry()
		return '-3'

	elif modes.get() == -4:
		vary_drop.configure(state="disabled")
		clearentry()
		nocalcentry()
		return '-4'

	elif modes.get() == 0:
		vary_drop.configure(state="disabled")
		clearentry()
		nocalcentry()
		return '0'



calcmoden = IntVar()
valcalc1 =Radiobutton(second_frame, text="1 (Calculate)", variable=modes, value=1, command=calcmodefunc)
valcalc1.grid(row=6, column=1, padx=2, sticky=W)

val0 =    Radiobutton(second_frame, text=" 0 (Fit data)   ", variable=modes, value= 0, command=calcmodefunc)
val0.grid(row=6, column=2, sticky=W)

valneg1 = Radiobutton(second_frame, text="-1 (Monte Carlo)", variable=modes, value= -1, command=calcmodefunc)
valneg1.grid(row=6, column=3, padx=2, sticky=W)

valneg3 = Radiobutton(second_frame, text="-3 (Eval. error)", variable=modes, value= -3, command=calcmodefunc)
valneg3.grid(row=6, column=4, padx=2, sticky=W)

valneg2 = Radiobutton(second_frame, text="-2 (Calc.& vary)", variable=modes, value= -2, command=calcmodefunc)
valneg2.grid(row=7, column=2, padx=2, sticky=W)

valneg4 = Radiobutton(second_frame, text="-4 (Grid search)", variable=modes, value= -4, command=calcmodefunc)
valneg4.grid(row=7, column=3, padx=2, sticky=W)

calcmode_n = Entry(second_frame, width=10, textvariable=calcmoden)
calcmode_n.grid(row=7, column=1, sticky=W)



###########################################################################################
###########################################################################################
# Main frame labels:

Name_text = StringVar()
header = Entry(second_frame, textvariable=Name_text, width = 60)
header.grid(row=1, column=1, columnspan=4, pady=5, sticky=W)

header_label = Label(second_frame, text="Sample Name")
header_label.grid(row=1, column=0, padx=10, pady=(5,0), sticky=W)

shape_label = Label(second_frame,text="Shape")
shape_label.grid(row=2, column=0, padx=10, sticky=W)

bravais_label = Label(second_frame, text="Bravais Class")
bravais_label.grid(row=3, column=0, padx=10, sticky=W)

maxpoly_label = Label(second_frame, text="Polynomial Order")
maxpoly_label.grid(row=5, column=0, padx=10, sticky=W)

calcmode_label = Label(second_frame, text="Calculation Mode")
calcmode_label.grid(row=6, column=0, padx=10, sticky=W)

mirror_label = Label(second_frame, text="Mirror Planes")
mirror_label.grid(row=9, column=0, padx=10, sticky=W)

chi2_label = Label(second_frame, text="Last Chi2:")
chi2_label.grid(row=2, column=2, padx=10, sticky=W)
lastchi2= StringVar()
chi2 = Entry(second_frame, textvariable=lastchi2, width=10)
chi2.grid(row=2, column=3, sticky=W)
chi2.config(state="disabled")
###########################################################################################
###########################################################################################
# Write-to-file button
metinF = Entry(second_frame)
#metinF = print("#" + '\n')

writeframe = LabelFrame(second_frame,text = "Handle input/output files", padx=10, pady=10)
writeframe.grid(row=32, column=5, rowspan=1, columnspan=11, sticky='nswe', padx=2, pady=0)

writebutton = Button(writeframe)
writebutton.config(text= "Write rusin.dat", command = writeFile, bg='grey')
writebutton.grid(row=0, column=5, padx=1,pady=8)

label_inputfile = Label(file_frame)
label_inputfile.grid(row=1, column=0,columnspan=5)
def open_inputfile():
	inputfile = filedialog.askopenfilename(initialdir=os.path.abspath(os.getcwd()),
											title="Select A File",
											filetypes=(("dat files", "*.dat"), ("All Files", "*.*")))
	label_inputfile.config(text="")
	label_inputfile["text"] = inputfile

readbutton0 = Button(writeframe)
readbutton0.config(text= "Choose input file", command = open_inputfile, bg='grey')
readbutton0.grid(row=0, column=8, padx=1,pady=8)

#readioframe = LabelFrame(second_frame)
#readioframe.grid(row=32, column=11, rowspan=1, columnspan=3, sticky='nswe', padx=2, pady=10)
readiobutton = Button(writeframe)
readiobutton.config(text= "Read input file", command = readFile, bg='grey')
readiobutton.grid(row=0, column=11, padx=1, pady=8)



root.mainloop()
