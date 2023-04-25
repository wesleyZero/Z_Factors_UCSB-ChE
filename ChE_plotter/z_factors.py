#Wesley Johanson
# from re import I
import ChE
import numpy as np

#Files to Load data from
file1 = "CSV/z_factors_CSV_1.csv"
files = [file1]

#Data labels
rowIndex_of_labels = 1
dataStartsAtIndex = 3
labels = np.loadtxt(file1, unpack=True, delimiter=',',dtype=str,skiprows=rowIndex_of_labels, max_rows=1)#[:, rowIndex_of_labels]	
savePlotAs = "plot.png"
folder = "Images/"

#This is the Column of data that tells you which dataset a member belongs to
#i.e. this is the data segment index (or number)
segmentDataCol = 0

#Ovveride the random colors function with your own custom colors
customColors = None

#Plot these columns of data with respect to the first column (the var)
plots = [	
				[ 5, 6] 		#fig 1
				,[ 4, 7]		#Fig 2
				,[ 4, 9]		#Fig 3
				
			]
# plots = None

#The columns of data that are associated with x and y error in the plot
error = [ 	
				[ 13, 14]
				,[ 10, 11] #NOT THE REAL ERROR
				,[ 9, 11] #NOT THE REAL ERROR
]

#Which data-segments do you want to plot?
pltTheseRuns = [ 1 , 2 , 3, 4, 5 , 6]
#Which column will I find the number corresponding to the data segment index/num?
segmentCols = [ 	0,
					0,
					0
				]

#Booleans to determine which data set segments have error bars plotted
errBarsDS = { 	1:False,
				2:False,
				3:False,
				4:False,
				5:False,
				6:False,
				7:True
			}

#Ovveride the y axis labels if you want something different from the datalabels
#in the CSV file 
yLabelOverride = [	'${P_{r-1}}/{P_{r}}$',
			"friction factor",
			"$V^{+}$",
			"poop3"
			] 
yLabelOverride = None

#Data Markers for each datasegment(dataset #)
seg_mks = [".",".", "v", "+","^",'3','2', '3', '2']
#Array of segment markers for each plot, for each dataset
_markers = [seg_mks for i in range(0,16)]

plots2Make = {
	"varCol":5, "fnCols":[6], "xErrCol":13, "yErrCols":[14], "segs2plot":[1,2,3], "segCol": 0, "errBarCols2plot":[], "yLabelOverride":False, "markers":None
}

#Make Plot Obj______________________________________________________
plot = ChE.ChEplot(file1, labelsRowIndex=rowIndex_of_labels, dataRowIndex=dataStartsAtIndex)
plot.setSettings(varCol=5, fnCols=[6], segCol=0)
plot.plotData()
plot.close(showPlot=True, saveFigAs="CSV_newPlotter_test01.png")

# # plot.setDataColors(customColors)
# # plot.setFxns2Plot(plotData)
# plot.plotData_str(width=6,height=6, markers=_markers, err=error[i],seg=segmentDataCol, pltErrBar=errBarsDS, pltSegs=pltTheseRuns)
# plot.plotLRegLine(varCol=5, fnCol=6,width=0.5, segCol=0, segGroups=pltTheseRuns)

# #Plot Parameters_______________________________________________________ 
# xaxisLabel = labels[plotData[0]] #Don't Change
# yaxisLabel = yLabelOverride[i] if yLabelOverride else labels[plotData[1]]
# plot.setAxisLabels(xaxisLabel, yaxisLabel, xpadding=5, ypadding=5)
# plot.setTicProps()
# # plot.setNumTics(delta_x=10, delta_y=10, x_subTics=3, y_subTics=3)
# plot.showLegend()
# # plot.changeFont()

# #Presentation__________________________________________________________
# # plot.showPlot()
# temp = folder + str(i) + '_' + savePlotAs
# plot.savePlot(filename=temp,_dpi=600)
# print(temp)
# plot.close()
# i += 1



# for file in files:
# 	i = 0
# 	for plotData in plots:
# 		#Make Plot Obj______________________________________________________
# 		plot = ChE.ChEplot(filename=file, labelsRowIndex=rowIndex_of_labels, dataRowIndex=dataStartsAtIndex)
# 		plot.setSettings()

# 		plot.setDataColors(customColors)
# 		plot.setFxns2Plot(plotData)
# 		plot.plotData_str(width=6,height=6, markers=_markers, err=error[i],seg=segmentDataCol, pltErrBar=errBarsDS, pltSegs=pltTheseRuns)
# 		plot.plotLRegLine(varCol=5, fnCol=6,width=0.5, segCol=0, segGroups=pltTheseRuns)

# 		#Plot Parameters_______________________________________________________ 
# 		xaxisLabel = labels[plotData[0]] #Don't Change
# 		yaxisLabel = yLabelOverride[i] if yLabelOverride else labels[plotData[1]]
# 		plot.setAxisLabels(xaxisLabel, yaxisLabel, xpadding=5, ypadding=5)
# 		plot.setTicProps()
# 		# plot.setNumTics(delta_x=10, delta_y=10, x_subTics=3, y_subTics=3)
# 		plot.showLegend()
# 		# plot.changeFont()

# 		#Presentation__________________________________________________________
# 		# plot.showPlot()
# 		temp = folder + str(i) + '_' + savePlotAs
# 		plot.savePlot(filename=temp,_dpi=600)
# 		print(temp)
# 		plot.close()
# 		i += 1

