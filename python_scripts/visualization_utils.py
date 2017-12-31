import sys
from igakit.nurbs import NURBS
from igakit.plot import plt
#from KratosMultiphysics import *
from KratosMultiphysics.IsogeometricApplication import *


class IgaKitVisualizer:
    def __init__(self, patch, backend = "matplotlib"): # backend can be 'mayavi' or 'matplotlib' or "null" or "none"
        self.backend = backend
        plt.use(self.backend)

        if patch.WorkingSpaceDimension() == 2:
            U = patch.FESpace().KnotU
            V = patch.FESpace().KnotV
            C = patch.GridFunction(CONTROL_POINT).ControlGrid.ControlValues
            self.nurbs = NURBS([V, U], C)
        elif patch.WorkingSpaceDimension() == 3:
            U = patch.FESpace().KnotU
            V = patch.FESpace().KnotV
            W = patch.FESpace().KnotW
            C = patch.GridFunction(CONTROL_POINT).ControlGrid.ControlValues
            self.nurbs = NURBS([W, V, U], C)

    def AddPlot(self, plot_cmd):
        if plot_cmd == "cplot":
            plt.cplot(self.nurbs)
        elif plot_cmd == "kplot":
            plt.kplot(self.nurbs)
        elif plot_cmd == "plot":
            plt.plot(self.nurbs)
        elif plot_cmd == "curve":
            plt.curve(self.nurbs)
        elif plot_cmd == "surface":
            plt.surface(self.nurbs)
        elif plot_cmd == "cpoint":
            plt.cpoint(self.nurbs)
        elif plot_cmd == "cwire":
            plt.cwire(self.nurbs)
        elif plot_cmd == "kpoint":
            plt.kpoint(self.nurbs)
        elif plot_cmd == "kwire":
            plt.kwire(self.nurbs)

    def AddPlots(self, title, plot_cmds):
        plt.figure()
        plt.title(title)
        for cmd in plot_cmds:
            self.AddPlot(cmd)

    def AxisEqual(self):
        if self.backend == "matplotlib":
            import matplotlib.pyplot as mpl
            mpl.axis("equal")

    def View(self, vvec):
#        if self.backend == "matplotlib":
#            import matplotlib.pyplot as mpl
#            mpl.view(vvec[0], vvec[1], vvec[2])
        pass

    def Render(self):
        plt.show()

