##########################################################################################
#
# midpoint.d - a variation on the midpoint displacement algorithm in multiple dimensions
#
# by Walt McNab
#
# Steps:
# (1) Read in a lattice; can consist of a single cell --> domain
# (2) Fill in the domain to specs to create a finer property resolution
#
##########################################################################################


from numpy import *
from pandas import *

options.mode.chained_assignment = None


class Vector:           # class for any 3-component spatial quantity
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z

        
class Lattice:

    def __init__(self, data):
        self.data = data                                            # lattice (assumed to be sorted by x, then by y, then by z)
        self.n, self.startLoc, self.endLoc = self.Grid(data)        # lattice properties (vectors)
        self.AssignConnections()                                    # populate lists of connecting lattice nodes
        
    def Grid(self, data):
        # lattice discretization
        self.x0 = sort(list(set(data['x'])))
        self.y0 = sort(list(set(data['y'])))
        self.z0 = sort(list(set(data['z'])))
        n = Vector(len(self.x0), len(self.y0), len(self.z0))
        startLoc = Vector(min(self.x0), min(self.y0), min(self.z0))        
        endLoc = Vector(max(self.x0), max(self.y0), max(self.z0))
        return n, startLoc, endLoc
        
    def AssignConnections(self):
        # add index numbers for (up-sequence) connecting lattice nodes
        self.data['selfIndex'] = self.data.index.values
        self.data['xIndexUp'] = (self.data.index.values+1)*(self.data['x']!=self.endLoc.x) - 1*(self.data['x']==self.endLoc.x)        # for all x-axis connections    
        self.data['yIndexUp'] = (self.data.index.values+self.n.x)*(self.data['y']!=self.endLoc.y) - 1*(self.data['y']==self.endLoc.y)        # for all y-axis connections
        self.data['zIndexUp'] = (self.data.index.values+self.n.x*self.n.y)*(self.data['z']!=self.endLoc.z) - 1*(self.data['z']==self.endLoc.z)        # for all z-axis connections

    def Combo(self, xCoord, yCoord, zCoord, xArray, yArray, zArray):
        # combinations of gridding schemes in 2D or 3D
        X, Y, Z = meshgrid(xCoord, yCoord, zCoord)
        xArray.extend(X.flatten())
        yArray.extend(Y.flatten())
        zArray.extend(Z.flatten())        
        return xArray, yArray, zArray

    def DiamondPts(self, N):
        # return arrays of coordinates for diamond-position lattice nodes for subsequent value population
        if sum(N>1) == 3:   # refinement of lattice along all three dimensions
            xArray = []
            yArray = []
            zArray = []
            xm = 0.5*self.x0[:-1] + 0.5*self.x0[1:]
            ym = 0.5*self.y0[:-1] + 0.5*self.y0[1:]                        
            zm = 0.5*self.z0[:-1] + 0.5*self.z0[1:]           
            xArray, yArray, zArray = self.Combo(self.x0, ym, zm, xArray, yArray, zArray)
            xArray, yArray, zArray = self.Combo(xm, ym, zm, xArray, yArray, zArray)            
            xArray, yArray, zArray = self.Combo(xm, self.y0, zm, xArray, yArray, zArray)            
            xArray, yArray, zArray = self.Combo(xm, ym, self.z0, xArray, yArray, zArray)            
            return xArray, yArray, zArray
        else:    
            if N[0]<2 and (N[1]>1 or N[2]>1): x = self.x0
            else: x = 0.5*self.x0[:-1] + 0.5*self.x0[1:]
            if N[1]<2 and (N[0]>1 or N[2]>1): y = self.y0
            else: y = 0.5*self.y0[:-1] + 0.5*self.y0[1:]                        
            if N[2]<2 and (N[0]>1 or N[1]>1): z = self.z0
            else: z = 0.5*self.z0[:-1] + 0.5*self.z0[1:]                        
            X, Y, Z = meshgrid(x, y, z)
            return X.flatten(), Y.flatten(), Z.flatten()

    def DiamondValues(self, data, N, stdev):
        # interpolate lattice diamond values based on surrounding nodes
        if N[0] > 1: nx = 2*self.n.x - 1
        else: nx = self.n.x
        if N[1] > 1: ny = 2*self.n.y - 1
        else: ny = self.n.y    
        data.reset_index(inplace=True)
        data['selfIndex'] = data.index.values
        diamondData = data[data['val']==-9999.]
        diamondData['total'] = 0.               # totaler and counter for computing average
        diamondData['count'] = 0
        # adjacent node variable will point to node 0 (valid value) if abutting edge, but will not be used ...  
        if N[0] > 1:
            nodeDown = (diamondData['selfIndex'] - 1) * (diamondData['x']!=self.startLoc.x)
            nodeUp = (diamondData['selfIndex'] + 1) * (diamondData['x']!=self.endLoc.x) 
            diamondData['total'] += array(data['val'][nodeDown])*array(diamondData['x']!=self.startLoc.x)*array(data['val'][nodeDown]!=-9999.)
            diamondData['total'] += array(data['val'][nodeUp])*array(diamondData['x']!=self.endLoc.x)*array(data['val'][nodeUp]!=-9999.)
            diamondData['count'] += 1*array(data['val'][nodeDown]!=-9999.)*(diamondData['x']!=self.startLoc.x) \
                + 1*array(data['val'][nodeUp]!=-9999.)*(diamondData['x']!=self.endLoc.x)
        if N[1] > 1:
            nodeDown = (diamondData['selfIndex'] - nx) * (diamondData['y']!=self.startLoc.y)
            nodeUp = (diamondData['selfIndex'] + nx) * (diamondData['y']!=self.endLoc.y)    
            diamondData['total'] += array(data['val'][nodeDown])*array(diamondData['y']!=self.startLoc.y)*array(data['val'][nodeDown]!=-9999.)
            diamondData['total'] += array(data['val'][nodeUp])*array(diamondData['y']!=self.endLoc.y)*array(data['val'][nodeUp]!=-9999.)
            diamondData['count'] += 1*array(data['val'][nodeDown]!=-9999.)*(diamondData['y']!=self.startLoc.y) \
                + 1*array(data['val'][nodeUp]!=-9999.)*(diamondData['y']!=self.endLoc.y)
        if N[2] > 1:
            nodeDown = (diamondData['selfIndex'] - nx*ny) * (diamondData['z']!=self.startLoc.z)
            nodeUp = (diamondData['selfIndex'] + nx*ny) * (diamondData['z']!=self.endLoc.z)        
            diamondData['total'] += array(data['val'][nodeDown])*array(diamondData['z']!=self.startLoc.z)*array(data['val'][nodeDown]!=-9999.)
            diamondData['total'] += array(data['val'][nodeUp])*array(diamondData['z']!=self.endLoc.z)*array(data['val'][nodeUp]!=-9999.)
            diamondData['count'] += 1*array(data['val'][nodeDown]!=-9999.)*(diamondData['z']!=self.startLoc.z) \
                + 1*array(data['val'][nodeUp]!=-9999.)*(diamondData['z']!=self.endLoc.z)
        stdevArray = random.normal(0, stdev, len(diamondData))
        diamondData['filled'] = 1*array(diamondData['count']!=0)
        diamondData['count'] += 1*array(diamondData['filled']==0)
        diamondData['val'] = diamondData['total']/diamondData['count'] + stdevArray*array(diamondData['filled']) \
            - 9999.*array(diamondData['filled']==0)
        diamondData = diamondData[['x', 'y', 'z', 'val']]
        return diamondData
       
        
def ReadParams():
    # basic complexification parameters
    lineInput = []        
    inputFile = open('params.txt','r')
    for line in inputFile: lineInput.append(line.split())
    inputFile.close()
    N = array([int(lineInput[0][1]), int(lineInput[1][1]), int(lineInput[2][1])])           # number of divisions along each axis (express as multiplier)
    stdev0 = float(lineInput[3][1])     # (starting) standard deviation used to perturb average value between adjacent nodes
    rStdevF = float(lineInput[4][1])    # reduction factor in standard deviation, per iteration
    print('Read model parameters.')
    return N, stdev0, rStdevF


def Splitter(data, direction, stdev):
    # create separate data frame consisting of midpoints along one axis, with posited new values
    if direction == 'x':
        subData = data[data['xIndexUp']!=-1]
        xMid =  array(0.5*data['x'][subData['selfIndex']]) + array(0.5*data['x'][subData['xIndexUp']])
        vMid =  array(0.5*data['val'][subData['selfIndex']]) + array(0.5*data['val'][subData['xIndexUp']])        
        subData['x'] = xMid
    elif direction == 'y':
        subData = data[data['yIndexUp']!=-1]
        yMid =  array(0.5*data['y'][subData['selfIndex']]) + array(0.5*data['y'][subData['yIndexUp']])
        vMid =  array(0.5*data['val'][subData['selfIndex']]) + array(0.5*data['val'][subData['yIndexUp']])        
        subData['y'] = yMid    
    else:              
        subData = data[data['zIndexUp']!=-1]
        zMid =  array(0.5*data['z'][subData['selfIndex']]) + array(0.5*data['z'][subData['zIndexUp']])
        vMid =  array(0.5*data['val'][subData['selfIndex']]) + array(0.5*data['val'][subData['zIndexUp']])        
        subData['z'] = zMid
    stdevArray = random.normal(0, stdev, len(subData))        
    subData['val'] = vMid + stdevArray
    return subData


def Midpoint():

    # read initial lattice and assign to object
    data = read_csv('lattice_input.csv', sep=',')
    lattice = Lattice(data)
    print('Read and processed starting lattice.')

    # read model parameters
    N, stdev, rStdevF = ReadParams()

    # iterate ...
    itr = 0
    while sum((N[0]+N[1]+N[2])>1):

        expanded = lattice.data         # initiate expanded lattice data frame
        
        # divide and assign lattice cell edges
        if N[0] > 1:
            xSubData = Splitter(lattice.data, 'x', stdev)
            expanded = concat([expanded, xSubData], axis=0)
        if N[1] > 1:
            ySubData = Splitter(lattice.data, 'y', stdev)
            expanded = concat([expanded, ySubData], axis=0)            
        if N[2] > 1:
            zSubData = Splitter(lattice.data, 'z', stdev)
            expanded = concat([expanded, zSubData], axis=0)

        if sum((N[0]+N[1]+N[2])>1)>1:
        
            # add-in blank lattice diamond node locations (calculated, but ignored, for 1-D-only refinement)
            xD, yD, zD = lattice.DiamondPts(N)                 
            diamond = DataFrame({'x':xD, 'y':yD, 'z':zD})
            diamond['val'] = -9999. # flag indicating node with unassigned value; used to avoid working with NaNs
 
            # process expanded lattice
            expanded = expanded[['x', 'y', 'z', 'val']]
            expanded = concat([expanded, diamond], axis=0)
            expanded = expanded[['x', 'y', 'z', 'val']]
            expanded.sort_values(by=['z', 'y', 'x'], axis=0, ascending=True, inplace=True, kind='quicksort')
            
            # two inner iterations to cover all unassigned-value lattice cells
            for i in range(2):
                diamondData = lattice.DiamondValues(expanded.copy(), N, stdev)             # fill in diamond values
                diamondData = diamondData[['x', 'y', 'z', 'val']]
                expanded = expanded[expanded['val']!=-9999.]
                expanded = concat([expanded, diamondData], axis=0)
                expanded.sort_values(by=['z', 'y', 'x'], axis=0, ascending=True, inplace=True, kind='quicksort')
                expanded.reset_index(inplace=True)

        else:
            expanded = expanded[['x', 'y', 'z', 'val']]
            expanded.sort_values(by=['z', 'y', 'x'], axis=0, ascending=True, inplace=True, kind='quicksort')
            expanded.reset_index(inplace=True)

        # refresh lattice object
        lattice = Lattice(expanded)
        
        # update iteration parameters
        N /= 2
        stdev *= rStdevF
        itr += 1
        print('Processing iteration #' + str(itr))
    
    # write output
    lattice.data[['x', 'y', 'z', 'val']].to_csv('lattice_output.csv', index=False, sep=',')
    
    print('Finished.')    
    
### run script ###
Midpoint()