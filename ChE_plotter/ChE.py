#Wesley Johanson
# from aifc import _Marker
# from re import I
from pprint import pprint
from re import I
from tokenize import group
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np 
import matplotlib.font_manager as fm
from sklearn.linear_model import LinearRegression
from scipy import stats 
import random 

class ChEplot:
	def loadCSV_str(self, filename: str, names: list, indepVars, skip=0):
		self.data = np.loadtxt(filename, unpack=True, delimiter=',',skiprows=skip, dtype=str)	
		self.dataLabels = names
		self.numDataVars = indepVars
		self.numDataSets = len(self.data) 
		self.numDataFns = self.numDataSets - self.numDataVars	
	
	#Setters	
	def setDataLabel(self, names):
		self.dataLabels = names
	def setLRegLineColors(self, colors=[]):
		self.LRegLineColors = colors
	def setIndepVars(self, vars):
		"Vars are the first arrays in the self.data matrix"
		self.setIndepVars = vars
		
	def setFxns2Plot(self, fxns):
		self.fxns2plot = fxns
		self.setFnCols(self.fxns2plot[1:])
		self.setVarCol(self.fxns2plot[0])

	def setDataStyles(self, styles: list):	
		self.lineStyles = styles
	def	setDataColors(self, colors=None): 	
		#Assign Random Colors, if no custom color array is given
		if colors == None: 
			colors = []
			for i in range (len(self.data)):
				color = "#"+''.join([random.choice('0123456789ABCDEF') for j in range(6)])
				colors.append(color)

		self.dataColors = colors 
		self.setLRegLineColors(colors)

	def setAxisLabels(self, x: str, y:str, indepVar=0, xpadding=5, ypadding=5):
		self.figure.axis[indepVar].set_xlabel(x,labelpad=xpadding)
		self.figure.axis[indepVar].set_ylabel(y,labelpad=ypadding)

	#Plot: Features	
	def showLegend(self, x=0.01, y=0.01, width=1, height=1, _loc='lower left', frame=True,_fontSize=10):
		# plt.legend(bbox_to_anchor=(x,y, width, height), loc=_loc, frameon=frame, fontsize=_fontSize)
		plt.legend(frameon=frame, fontsize=_fontSize)

	#Presentation & Demonstration
	def showPlot(self): 
		"Shows the figure"
		plt.show() 
	#Save
	def savePlot(self, filename="poop.png", _dpi=900, _transparent=False, _bbox_inches='tight'):
		"Saves the Figure(graph) made to a file"
		plt.savefig(filename, dpi=_dpi, transparent=_transparent, bbox_inches=_bbox_inches)



	
	def plotData_str(self, width, height, markers=None, err=None, xScale="", yScale="", seg=None, pltErrBar=False, pltSegs=None):
		self.figure = plt.figure(figsize=(width, height))
		L, B, W, H = [0.15, 0.1, 0.80, 0.85]
		self.figure.axis = []
		self.figure.axis.append(self.figure.add_axes([L, B, W, H]))
		var = self.fxns2plot[0]
		for fn in self.fxns2plot[1:]:
			# print("ORIGINAL X Y ", self.data[var], "\n", self.data[fn])
			#Segment the data if neccessary
			if seg is not None:
				segSet = [int(s) for s in self.data[seg]]
				segSet = set(segSet)
				segSet = list(segSet)
				Xseg = { i:[] for i in segSet} 
				Yseg = { i:[] for i in segSet} 
				Xseg_err = { i:[] for i in segSet} 
				Yseg_err = { i:[] for i in segSet} 

				for i in range(0,len(self.data[var])):
					if (self.data[var][i] != "" and self.data[fn][i] != ""):
						Xseg[int(self.data[seg][i])].append(self.data[var][i])
						Yseg[int(self.data[seg][i])].append(self.data[fn][i])
						Xseg_err[int(self.data[seg][i])].append(float(self.data[err[0]][i]))
						Yseg_err[int(self.data[seg][i])].append(float(self.data[err[1]][i]))
				for i in segSet:
					if i not in pltSegs: 
						continue

					x = Xseg[i] 
					y = Yseg[i]
					new_x = [float(i) for i in x]
					new_y = [float(i) for i in y]
					if len(new_x) == 0: continue #DELETE ? 
					x_error = Xseg_err[i]
					y_error = Yseg_err[i]
					if markers is not None: 
						mk = markers[fn][i]
					else:
						mk = "."
					lbl = self.dataLabels[fn] + "_" + str(i)
					#Colors
					if self.dataColors is not None:
						clr = ChEplot.getRandColor()
						self.figure.axis[0].plot(new_x,new_y,mk,label=lbl,color=clr)
					else: #I don't think this condition is possible anymore 
						self.figure.axis[0].plot(new_x,new_y,mk,label=lbl)
					#ERROR BARS MODIFY THIS FOR SEGMENETED DATA
					if err is not None and pltErrBar[i] is True:
						self.figure.axis[0].errorbar(new_x, new_y, xerr=x_error, yerr=y_error ,linewidth=0.9, fmt='none') #,ecolor=None) # fmt='none')
					#Scale 
					if xScale=="log":
						self.figure.axis[0].set_xscale('log')
					if yScale=="log":
						self.figure.axis[0].set_yscale('log')	
			#No segmentation
			else:
				x = self.data[var]
				y = self.data[fn]
				lbl = self.dataLabels[fn]
				#Remove Data points with blank strings that don't map to floats
				for col in range(0, len(self.data)):
					new_x = [] 
					new_y = []
					x_error = []
					y_error = [] 
					for row in range(0, len(self.data[0])):
						if self.data[col,row] != "" and self.data[var,row] !="":
							print('row =', row, 'col = ', err[0], "\tvalue = " ,self.data[err[0], row])
							if err is not None:
								x_error.append(float(self.data[err[0], row]))
								y_error.append(float(self.data[err[1], row]))
							new_x.append(float(self.data[var,row]))
							new_y.append(float(self.data[col,row])) 	
				#Set markers for plotting
				if markers is not None: 
					mk = markers[fn]
				else:
					mk = "."
				#Colors
				if self.dataColors is not None:
					clr = self.dataColors[fn]
					self.figure.axis[0].plot(new_x,new_y,mk,label=lbl,color=clr)
				else: #I don't think this condition is possible anymore 
					self.figure.axis[0].plot(new_x,new_y,mk,label=lbl)
				#ERROR BARS
				if err is not None:
					plt.errorbar(new_x, new_y, xerr=x_error, yerr=y_error)
				#Scale 
				if xScale=="log":
					self.figure.axis[0].set_xscale('log')
				if yScale=="log":
					self.figure.axis[0].set_yscale('log')	

	#NEW REFACTORED CHEPLOTTER__________________________________________________
	def __init__(self, filename:str, labelsRowIndex=0, dataRowIndex=1):
		self.setData(np.loadtxt(filename, unpack=True, delimiter=',',skiprows=dataRowIndex, dtype=str))
		self.setLabels(np.loadtxt(filename, unpack=True, delimiter=',',dtype=str,skiprows=labelsRowIndex, max_rows=1))
		self.figure=None
		#LEGACY 
		self.dataColors=None
		self.numDataVars=None
		self.numDataFns=None
		self.numDataSets=None
		self.fxns2plot=None

	def close(self, showPlot=False, saveFigAs=None, _dpi=900, _transparent=False,_bbox_inches='tight'): 
		if showPlot: plt.show()
		if saveFigAs: plt.savefig(saveFigAs,dpi=_dpi,transparent=_transparent, bbox_inches=_bbox_inches)  
		plt.close(self.figure)

	#Setters
	def	setData(self, data):
		self.data = data
	def setLabels(self, labels):
		self.labels = labels	
	def	setFnCols(self, cols:[int]):
		self.fnCols=cols
	def setVarCol(self, varCol:int):
		self.varCol = varCol
	def setSettings(self, varCol: int, fnCols:[int], segCol=None, colors=None, markers=None, ticProps={}, numTics={}, errBars={}, font={}):
		self.varCol = varCol
		self.fnCols = fnCols
		self.segCol = segCol
		self.setColors(colors)
		self.setMarkers(markers)



	#Setters: Optionals
	def setSegCol(self, segCol=None):
		self.segCol = segCol
		self._setSegmentKeys()
	def _setSegmentKeys(self):
		segmentedKeysArr = self.segmentData(self.data[self.segCol], self.segCol)
		keySet = ChEplot.getKeysFromDict(segmentedKeysArr)
		self.keyset = keySet
	def setMarkers(self, markers=None):
		self.markers = markers
	#??????????
	def setErrBars(self, errBars=None):
		self.errBars = errBars
	def	setColors(self, colors=None):
		if colors:
			self.colors = colors 
		else:
			#gen rand colors
			pass
	def setTicProps(self, _size=4, _width=1, _direction='in'):
		self.figure.axis[0].xaxis.set_tick_params(which='major', size=_size, width=_width, direction=_direction, top='on')
		self.figure.axis[0].xaxis.set_tick_params(which='minor', size=_size, width=_width, direction=_direction, top='on')
		self.figure.axis[0].yaxis.set_tick_params(which='major', size=_size, width=_width, direction=_direction, right='on')
		self.figure.axis[0].yaxis.set_tick_params(which='minor', size=_size, width=_width, direction=_direction, right='on')
	def setNumTics(self, delta_x=0.1, delta_y=0.1, x_subTics=3, y_subTics=3):
		self.figure.axis[0].xaxis.set_major_locator(mpl.ticker.MultipleLocator(delta_x))
		self.figure.axis[0].xaxis.set_minor_locator(mpl.ticker.MultipleLocator(delta_x/x_subTics))
		self.figure.axis[0].yaxis.set_major_locator(mpl.ticker.MultipleLocator(delta_y))
		self.figure.axis[0].yaxis.set_minor_locator(mpl.ticker.MultipleLocator(delta_y/y_subTics))
	def changeFont(self, font='Alvenir', size=10, linewidth=0.9):
		mpl.rcParams['font.family'] = font
		plt.rcParams['font.size'] = size 
		plt.rcParams['axes.linewidth'] = linewidth 
	def setLegend(self, x=0.01, y=0.01, width=1, height=1, _loc='lower left', frame=True,_fontSize=10):
		# plt.legend(bbox_to_anchor=(x,y, width, height), loc=_loc, frameon=frame, fontsize=_fontSize)
		plt.legend(frameon=frame, fontsize=_fontSize)


	#Getters	
	def getColor(self, col:int, row:int):
		return 'b'
	def getMarker(self, col:int, row:int):
		return '.'
	def getLabel(self):
		return "label"
	
	#Plotters
	def	makeFigure(self, width=6, height=6, L=0.15, B=0.1, W=0.80, H=0.85, axes=None):
		self.figure = plt.figure(figsize=(width, height))
		plotDimensions = [L, B, W, H]
		#Figure out how to use more axes 
		self.figure.axis = []
		self.figure.axis.append(self.figure.add_axes(plotDimensions))

	def plotData(self, width=1, height=1, markers=None, varErrCol=None, xScale=None, yScale=None, plotErrorBars=False, segGrps2plt=None):
		#Make the Figure object to create the plot on 
		self.makeFigure(width=width, height=height)
		#Plot all of the functions in the list of function columns
		for fnCol in self.fnCols:
			if self.segCol:
				varDict = self.segmentData(self.varCol,self.segCol)
				fnDict = self.segmentData(fnCol, self.segCol)
				#error here 

				for key in varDict:
					if key not in segGrps2plt: continue
					var, fn = ChEplot.removeBothElem(varDict[key], fnDict[key], rm="")				 
					var = ChEplot.str2floatArr(var)
					fn = ChEplot.str2floatArr(fn)

					marker = self.getMarker(fnCol, key)
					label = self.getLabel(fnCol, key)
					color = self.getColor()
					if plotErrorBars:
						pass		
					self.figure.axis[0].plot(var, fn, marker, label, color)
					#CONTINUE HERE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	def plotLRegLine(self,varCol:int, fnCol:int, colors=None, width=0.5, style='-',printMB=True, segCol=None, segGroups=None):
		#Is the data plotted in segments? 
		if segCol is not None:
			varDict = self.segmentData(varCol,segCol)
			fnDict = self.segmentData(fnCol, segCol)	

			#make a set of random colors for each key group, if no clrs given
			colors = colors if colors else { i:ChEplot.getRandColor() for i in varDict }

			#Plot each segment group individually
			for key in varDict:
				if segGroups and key not in segGroups: continue
				var = varDict[key]
				fn = fnDict[key]
				#Remove invalid (empty) pairs, convert to float for plotting
				var, fn = ChEplot.removeBothElem(var, fn, '')
				var, fn = ChEplot.str2floatArr(var), ChEplot.str2floatArr(fn)
				#Plot the entire column of data 
				x, y = ChEplot.linRegLine(var, fn)
				m, b = ChEplot.linRegCoeff(var, fn)
				# Plot the segmented dataset
				self.plotLine(x, y, clr=colors[key])
				#Output the line data
				if printMB: print("varCol=", varCol, "_fnCol=", fnCol,\
										"_segGroup=", key, "_m=", m,"_b=",b,"R^2=",ChEplot.rSquared(var,fn )) 
		else: #If no segments, plot the entire column of data
			var = self.data[varCol]
			fn = self.data[fnCol]
			#Remove invalid (empty) pairs, convert to float for plotting
			var, fn = ChEplot.removeBothElem(var, fn, '')
			var, fn = ChEplot.str2floatArr(var), ChEplot.str2floatArr(fn)
			#Plot the entire column of data 
			x, y = ChEplot.linRegLine(var, fn)
			m, b = ChEplot.linRegCoeff(var, fn)
			self.plotLine(x, y, clr=self.dataColors[fnCol])
			#Output the line data
			if printMB: print("varCol=", varCol, "_fnCol=", fnCol,\
									"m=",m,"_b=",b,"R^2=",ChEplot.rSquared(var,fn )) 


	def plotLine(self, x:[float], y:[float], clr:str, width=0.5, style='-'):
		if self.dataColors is None: self.setDataColors()
		#Plot it
		self.figure.axis[0].plot(x, y, color=clr,\
											linewidth=width ,linestyle=style)

	#Helper Functions 
	@staticmethod
	def	getKeysFromDict(x: dict):
		keySet = [key for key in x]
		keySet = set(keySet)
		keySet = list(keySet)
		return keySet

	@staticmethod
	def getRandColor():
		return  "#"+''.join([random.choice('0123456789ABCDEF') for j in range(6)])

	#Helper functions: Statistics
	@staticmethod
	def linRegLine(var, fn):
		m, b = np.polyfit(var,fn,deg=1)
		x = np.linspace(min(var), max(var), num=len(var))
		# x = np.array(object=var, dtype=float)
		y = m * x + b
		return x, y 

	@staticmethod
	def linRegCoeff(var, fn):
		m, b = np.polyfit(var,fn,deg=1)
		# x = np.linspace(min(var), max(var), num=len(var))
		return m, b 

	def	segmentData(self, dataCol:int, segCol:int):
	 	#Find the set of all segment groups
		# print('segCol = ', segCol,"_SEGMENTATION COL = ", self.data[segCol])	
		segGroupsSet = [int(s) for s in self.data[segCol]]
		segmentedDataDict = { group:[] for group in segGroupsSet} 
		for i in range(0, len(self.data[dataCol])):
			segGroupIndex = int(self.data[segCol][i])
			segmentedDataDict[segGroupIndex].append(self.data[dataCol][i])
		return segmentedDataDict

	@staticmethod
	def	rSquared(x, y):
		x = np.array(x)
		y = np.array(y)

		m, b = np.polyfit(x,y,deg=1)

		fitMe = [ [x[i],y[i]] for i in range(0, len(x)) ]
		y_fit = [ [x[i], m * x[i] + b] for i in range(0, len(x))]

		y_reg = LinearRegression().fit(fitMe, y_fit)
		return y_reg.score(fitMe,y_fit)

	#COMPLETE ME
	@staticmethod
	def confInterv(self, n=1):
		self.lowerBound_CI = []
		self.upperBound_CI = []
		for fn in range(self.numDataVars, self.numDataSets): 
			x_bar , stdDev = np.mean(self.data[fn]) ,np.std(self.data[fn])
			SE = stdDev / math.sqrt(self.num)
			DoF = n
			stats.t.ppf(q=0.05,)
			scipy.stats.t.ppf(q=.05,df=22)

	#Helper Functions: arrays
	@staticmethod
	def str2intArr(x:[str]):
		floatArr = [int(i) for i in x] 
		return floatArr

	@staticmethod
	def str2floatArr(x:[str]):
		floatArr = [float(i) for i in x] 
		return floatArr

	@staticmethod
	def removeBothElem(x:[str],y:[str], rm: str):
	# """ For each corresponding row of x and y
	# 		remove any rows in both, if either have the str 'rm'
	# """
		new_x = []
		new_y = []
		#Input error
		if rm == None: return
		max_length = len(x) if len(x) < len(y) else len(y) 
		for i in range(0, max_length):
			#Add all elements in x that don't have rm
			if x[i] == rm or y[i] == rm:  continue
			new_x.append(x[i])
			new_y.append(y[i])
		return new_x, new_y

	@staticmethod
	def removeElem(x:[str], rm: str):
		new_x = []
		#Input error
		if rm == None: return
		for i in range(0, len(x)):
			#Add all elements in x that don't have rm
			if x[i] == rm:  continue
			new_x.append(x[i])
		return new_x
	