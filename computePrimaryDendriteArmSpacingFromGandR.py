from functionComputePDAS import functionComputePDAS_NiSuperAlloys as computePDAS_NiSuperAlloys
import matplotlib.pyplot as plt


def main():
    nArguments = len(sys.argv)
    if nArguments != 2:
        text = """One command line arguments are expected: \
                \n\t(1) file name """
        raise RuntimeError(text)

    nameFile = str(sys.argv[1])
    plotScatterPlotsIn3D = False
    plotRegionsSeparately = False
    plotLayersSeparately = False
    plotRadialDataSeparately = False
    
    nPointsToPlot = 10000
    multiplierThermalGradient = 1000
    idColumnXcoord = 3
    idColumnYcoord = 4
    idColumnZcoord = 5
    idColumnThermalGradient = 15
    idColumnCoolingRate = 16
    idColumnIdLayer = 17
    idColumnIdRegion = 18
    idColumnNumberOfRemelts = 19
    
    listIdLayersToExclude = []
    listIdRegionsToExclude = []
    plotRadialDataSeparately = True
    
multiplierPDASforPlot = 1E6 # micro m to m conversion
multiplierPDASforPlotYaxis = 1E4
listNamesModels = ["Trivedi", "KurzFisher"]
# -------------- Experimental data
shouldPlotExperimentalResults         = 'true'
listGradientThermalExperiments        =   [4559100.8226608085,15228964.000000]
listVelocitySolidificationExperiments =   [0.6000000000,0.2279]
listPDASExperiments                   =   [0.85658E-6,0.74875E-6]
listPDASExperiments = [multiplierPDASforPlot*x for x in listPDASsimulation]


fontTimes = {'fontname':'Times New Roman'}
sizeFont1 = 20
sizeFont2 = 16
plt.rcParams.update({'font.size': sizeFont1})

# ----------------- Regionwise material properties
listFreezingRangeEquilibriumRegions        =   [66.18,66.18,57.836,118.892,113.398,110.245,100]
listFreezingRangeNonEquilibriumRegions     =   [234.33,234.33,209.62,227.94,234.6,227.74,390]
listDiffusivityLiquidRegions               =   [3E-9,3E-9,3E-9,3E-9,3E-9,3E-9,3E-9] # m^2/s
listCoefficientGibbsThompsonRegions        =   [1E-7,1E-7,1E-7,1E-7,1E-7,1E-7,1E-7] # K-m
listCoefficientSolutePartitionRegions      =   [0.48,0.48,0.48,0.48,0.48,0.48,0.48]
listCoefficientHarmonicPerturbationRegions =   [28,28,28,28,28,28,28]
variableOnXAxis,PDAS_HuntModel,PDAS_KGTModel,PDAS_HuntsModifiedModel,ratesCooling = ([] for ii in range(5))


# Compute the PDAS for each combination of Solidification velocity, thermal gradient and model
for iCombination in range(len(listVelocitySolidification)):
    velocitySolidification = listVelocitySolidification[iCombination]
    gradientThermal = listGradientThermal[iCombination]
    variableOnXAxis.append([multiplierPDASforPlotYaxis*(velocitySolidification**(-0.37))*(gradientThermal**(-1/2))])
    ratesCooling.append([velocitySolidification*gradientThermal])
    for iModel in range(len(listNamesModels)):
        nameModel = listNamesModels[iModel]
        PDAS = computePDAS_C103(coefficientSolutePartition,freezingRangeEquilibrium,   \
                           diffusivityLiquid,coefficientGibbsThompson,    \
                           velocitySolidification,gradientThermal,nameModel)
        if nameModel == "Hunt":
            PDAS_HuntModel.append([multiplierPDASforPlot*PDAS])
        elif nameModel == "KurzFisher":
            PDAS_KGTModel.append([multiplierPDASforPlot*PDAS])
        elif nameModel == "HuntsModified":
            PDAS_HuntsModifiedModel.append([multiplierPDASforPlot*PDAS])            
        else:
            print("Error: This model name - " + nameModel + " - is not valid")
            quit()
        print('Primary dendrite arm spacing (PDAS) using '+ nameModel + " model is : ",PDAS)
        
for iModel in range(len(listNamesModels)):
    if (listNamesModels[iModel]) == "Hunt":
        plt.plot(variableOnXAxis,PDAS_HuntModel,'k-',label="Hunt",linewidth=2)
    elif (listNamesModels[iModel]) == "KurzFisher":
        plt.plot(variableOnXAxis,PDAS_KGTModel,'k--',label="Kurz & Fisher",linewidth=2)
    elif (listNamesModels[iModel]) == "HuntsModified":
        plt.plot(variableOnXAxis,PDAS_HuntsModifiedModel,'k-',label="Hunts Modified",linewidth=2)
    else:
        print("Error: This model name - " + listNamesModels[iModel] + " - is not recognized while plotting")
        quit()
if shouldPlotSimulationResults:
    plt.plot(variableOnXAxis,listPDASsimulation,'ko',label="PF-Simulation",linewidth=4.0)
    
#plt.title('Primary Dendrite Arm Spacing')
plt.legend(loc='upper left',fontsize=sizeFont2,bbox_to_anchor=(0.0, 1.0))
plt.xlabel("$R^{-\dfrac{1}{4}}\,G^{-\dfrac{1}{2}}\, (x\,10^{-4})$",**fontTimes,fontsize=sizeFont1)
plt.ylabel("PDAS ($\mu$m)",**fontTimes,fontsize=sizeFont1)
plt.tight_layout()
plt.savefig('primaryDendriteArmSpacing.png')
plt.show()