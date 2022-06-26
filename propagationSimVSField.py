import matplotlib.pyplot as plt
import numpy as np
from time import time

showEnd = False
generateArrows = True
generateArrowsSim = True
generateAdditionalArrows = False
folder = "0-stream"
applyBorder = True

# Width and height of nodes, number of nodes equals wN*wH
# Please only use even numbers
wN = 200
hN = 200

# Reduction factor in arrow generation (integers only!)
rfa = 1

# Width and height of the view / coordinate Space
w = wN*2
h = hN*2

# Shift positions
dx = -w/2
dy = -h/2

# Scaling
dl = 2/h

# Create position array with charges
#xCh = [w//2-w//8, w//2+w//8, w//2,      w//2,      ]
#yCh = [h//2,      h//2,      h//2+h//8, h//2-h//8, ]
#cha = [1,         1,         -1,         -1,       ]

# Create postions on ellipse
xCh = list()
yCh = list()
cha = list()
a = w//8
b = h//8
steps=32
for i in range(steps):
    angle=np.pi*2.*i/steps
    xCh.append(int(w/2+a*np.cos(angle)))
    yCh.append(int(h/2+b*np.sin(angle)))
    cha.append((-1)**i)
print(xCh, yCh, cha)

# Potential
# Multiplier for resolution of image of the potential
potentialScale = pS = 5
# Define the number of countour lines
# With a higher resoulution of the potential and grid, more lines
# are concentrated near the charges!
noOfContours = 20
contourPotMin = -10


# Transformation function on the coordinate System.
def transform(x,y):
    return x,y
    #return np.sign(x)*np.abs(x**1.5)+np.sign(y)*np.abs(y**1.5), np.sign(y)*np.abs(y**3)
    #if x == 0:
    #    return 1, y
    #return np.abs(x)**np.abs(x), y

def transformB(x,y):
    global dx, dy, dl
    return transform((x+dx)*dl, (y+dy)*dl)

def step(arr, w, h, sinkX, sinkY):
    arr2= [[arr[y][x] for x in range(h)] for y in range(w)]
    out  = [[0 for _ in range(h)] for __ in range(w)]
    print(out is arr, arr2 is arr, out is arr2)
    while True:
        arr = [[arr2[y][x] for x in range(h)] for y in range(w)]
        arr2= [[arr[y][x] for x in range(h)] for y in range(w)]
        for x in range(1,w-1):
            for y in range(1,h-1):
                #if x == sinkX:
                #    if y == sinkY:
                #        continue
                v = arr[x][y]
                arr2[x  ][y  ] -= v
                arr2[x-1][y  ] += v*0.25
                arr2[x  ][y-1] += v*0.25
                arr2[x+1][y  ] += v*0.25
                arr2[x  ][y+1] += v*0.25
                out[x-1][y  ] -= v*0.25
                out[x  ][y-1] -= v*0.25
                out[x+1][y  ] += v*0.25
                out[x  ][y+1] += v*0.25
        yield arr2, out

# Step yielding function
def stepNP(initArr, up, down, left, right, shorts=None):
    arr = np.array(initArr)
    out = np.zeros(arr.shape)
    upR = np.zeros(arr.shape)
    doR = np.zeros(arr.shape)
    leR = np.zeros(arr.shape)
    riR = np.zeros(arr.shape)
    shortsMask = np.zeros(arr.shape)
    mask= np.ones(arr.shape)
    #mask= np.zeros(arr.shape)
    #mask[1:-1,1:-1] = np.ones((arr.shape[0]-2, arr.shape[1]-2))
    #mask[100,200]=0
    while True:
        arr = arr*mask
        upR = np.roll(arr*up,   -1, 0)
        doR = np.roll(arr*down,  1, 0)
        leR = np.roll(arr*left, -1, 1)
        reR = np.roll(arr*right, 1, 1)
        out = out + upR - doR - leR + reR
        upR = np.roll(upR, -1, 0)
        doR = np.roll(doR,  1, 0)
        leR = np.roll(leR, -1, 1)
        reR = np.roll(reR, 1, 1)
        if shorts is not None:
            shortsMask = np.zeros(arr.shape)
            for sSX, sSY, sEX, sEY, sP in shorts:
                shortsMask[sEY, sEX] += sP * arr[sSY, sSX]
        
        arr = upR + doR + leR + reR + shortsMask
        yield arr, out

def normalize(u, d, l, r, w, h, shorts=None):
    for x in range(w):
        for y in range(h):
            sum = u[y,x]+d[y,x]+l[y,x]+r[y,x]
            if shorts is not None:
                sSum = 0
                indices = list()
                for i, (sSX, sSY, sEX, sEY, sP) in enumerate(shorts):
                    if sSX == x and sSY == y:
                        indices.append(i)
                        sSum += sP
                sum += sSum
            
            #print(x, y, sum, r[y,x], d[y,x], l[y,x], r[y,x])
            u[y,x] = u[y,x]/sum
            d[y,x] = d[y,x]/sum
            l[y,x] = l[y,x]/sum
            r[y,x] = r[y,x]/sum
            
            if shorts is not None:
                for i in indices:
                    if shorts[i][4] == np.inf:
                        shorts[i][4] = 1
                    else:
                        shorts[i][4] = shorts[i][4]/sum
            #print(x, y, sum, r[y,x], d[y,x], l[y,x], r[y,x])
    return

def distances(u, d, l, r, w, h, dx=0, dy=0, dl=1, applyBorder=True, applyNormalization=True, shorts=None):
    #dx, dy, deltaL in coordinate space
    for x in range(w):
        for y in range(h):
            xReal = (x+dx) * dl
            yReal = (y+dy) * dl
            
            # Up
            xPos, yPos = transform(xReal, yReal)
            xDPos,yDPos= transform(xReal, yReal-dl)
            xDelta = xPos - xDPos
            yDelta = yPos - yDPos
            u[y,x] = 1./np.sqrt(xDelta ** 2 + yDelta ** 2) if np.sqrt(xDelta ** 2 + yDelta ** 2) != 0 else 1
            
            # Down
            xPos, yPos = transform(xReal, yReal)
            xDPos,yDPos= transform(xReal, yReal+dl)
            xDelta = xPos - xDPos
            yDelta = yPos - yDPos
            d[y,x] = 1./np.sqrt(xDelta ** 2 + yDelta ** 2) if np.sqrt(xDelta ** 2 + yDelta ** 2) != 0 else 1
            
            # Left
            xPos, yPos = transform(xReal, yReal)
            xDPos,yDPos= transform(xReal-dl, yReal)
            xDelta = xPos - xDPos
            yDelta = yPos - yDPos
            l[y,x] = 1./np.sqrt(xDelta ** 2 + yDelta ** 2) if np.sqrt(xDelta ** 2 + yDelta ** 2) != 0 else 1
            
            # Right
            xPos, yPos = transform(xReal, yReal)
            xDPos,yDPos= transform(xReal+dl, yReal)
            xDelta = xPos - xDPos
            yDelta = yPos - yDPos
            r[y,x] = 1./np.sqrt(xDelta ** 2 + yDelta ** 2) if np.sqrt(xDelta ** 2 + yDelta ** 2) != 0 else 1e5
    if applyBorder:
        u[:2,:]  = 0
        d[-2:,:] = 0
        l[:,:2]  = 0
        r[:,-2:] = 0
    if applyNormalization:
        normalize(u,d,l,r,w,h,shorts=shorts)

def norm(x,y):
    # 2-Norm
    r2 = x**2+y**2
    r  = np.sqrt(r2)
    
    # 1-Norm
    #r  = np.abs(x)+np.abs(y)
    #r2 = r**2
    
    return r2, r

def calcPotential(x,y,xCh, yCh, cha, scale=1):
    result = 0
    for xC, yC, q in zip(xCh, yCh, cha):
        xRC, yRC = transformB(xC, yC)
        xR, yR = transformB(x/scale,y/scale)
        r2, r = norm(xRC-xR, yRC-yR)
        if r == 0:
            continue
        result += q / r / (2*np.pi)
    return result

def genPotential(xCh, yCh, cha, shape, scale):
    result = np.zeros(shape)
    for xC, yC, q in zip(xCh, yCh, cha):
        xRC, yRC = transformB(xC, yC)
        for y in range(0, shape[0], 1):
            for x in range(0, shape[1], 1):
                if xC == x and yC == y:
                    continue
                xR, yR = transformB(x/scale,y/scale)
                r2, r = norm(xRC-xR, yRC-yR)
                if r == 0:
                    continue
                result[y,x] += q / r / (2*np.pi)
    return result

def genFields(xCh, yCh, cha, shape):
    resultX = np.zeros(shape)
    resultY = np.zeros(shape)
    for xC, yC, q in zip(xCh, yCh, cha):
        xRC, yRC = transformB(xC, yC)
        for y in range(1, shape[0], 2):
            for x in range(2, shape[1], 2):
                if xC == x and yC == y:
                    continue
                xR, yR = transformB(x,y)
                r2, r = norm(xRC-xR,yRC-yR)
                if r == 0 or r2 == 0:
                    continue
                dir = ((xR-xRC)-1j*(yR-yRC))/r
                resultX[y,x] += np.real(dir * q / r2 / (2*np.pi))
                resultY[y,x] += np.imag(dir * q / r2 / (2*np.pi))
        for y in range(2, shape[0], 2):
            for x in range(1, shape[1], 2):
                if xC == x and yC == y:
                    continue
                xR, yR = transformB(x,y)
                r2, r = norm(xRC-xR,yRC-yR)
                if r == 0 or r2 == 0:
                    continue
                dir = ((xR-xRC)-1j*(yR-yRC))/r
                resultX[y,x] += np.real(dir * q / r2 / (2*np.pi))
                resultY[y,x] += np.imag(dir * q / r2 / (2*np.pi))
    return resultX, resultY

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    #idx = np.unravel_index(idx, array.shape)
    return array.flatten()[idx]


if not (wN//2 == wN/2 and hN//2 == hN//2):
    print("Width or height is not even! Exiting")
    input("Press ENTER to exit")
    exit()

XGrid = np.array(range(w))*pS
YGrid = np.array(range(h))*pS
XGridB= np.array(range(w*pS))
YGridB= np.array(range(h*pS))



# Aligning chages to grid
for i in range(len(xCh)):
    if xCh[i] % 2 == 0:
        xCh[i] += 1
for i in range(len(yCh)):
    if yCh[i] % 2 == 1:
        yCh[i] += 1


arr = np.array([[0. for _ in range(w)] for __ in range(h)])
for x, y, q in zip(xCh, yCh, cha):
    arr[y,x] = q
#arr[h//2][w//2-w//8]=1
#arr[h//2][w//2+w//8]=-1


u = np.zeros(arr.shape)
d = np.zeros(arr.shape)
l = np.zeros(arr.shape)
r = np.zeros(arr.shape)

shorts=None
#shorts=[
# #FromX, FromY, ToX, ToY, inverse Distance (Zero Distance => np.inf)
#[w//2-w//4, h//2-h//8, w//2+w//4, h//2, np.inf],
# #Repeat entry with switched From and To, to create a 2-way link
##[w//2+w//4, h//2, w//2-w//4, h//2, np.inf],
# 
#[w//2+w//16, h//2+h//8, w//2-w//16, h//2, 4],
#]

shortsGrid = None
if shorts is not None and shorts != list():
    shortsGrid = list()
    for sSX, sSY, sEX, sEY, sP in shorts:
        shortsGrid.append([sSX, sSY, sEX, sEY, sP])
        shortsGrid.append([sSX+1, sSY, sEX+1, sEY, sP])
        shortsGrid.append([sSX+1, sSY+1, sEX+1, sEY+1, sP])
        shortsGrid.append([sSX  , sSY+1, sEX,   sEY+1, sP])
distances(u,d,l,r,w,h,dx=dx,dy=dy,dl=dl,applyBorder=applyBorder ,applyNormalization=True, shorts=shortsGrid)


# Display transportmation maps
if input("Enter 'y' to show diffusion maps ").lower() == "y":
    plt.title("UP - Map")
    plt.imshow(u)
    plt.colorbar()
    plt.show()
    plt.title("DOWN - Map")
    plt.imshow(d)
    plt.colorbar()
    plt.show()
    plt.title("LEFT - Map")
    plt.imshow(l)
    plt.colorbar()
    plt.show()
    plt.title("RIGHT - Map")
    plt.imshow(r)
    plt.colorbar()
    plt.show()

# Field strengths
calcResX, calcResY = genFields(xCh, yCh, cha, arr.shape)
calcRes = np.sqrt(calcResX**2+calcResY**2)
plt.imshow(calcRes)
plt.colorbar()
plt.savefig("{}/strenghtsCalc.png".format(folder), dpi=500)
plt.close()




# Generate the potential
print("--- Generating Potential Field ---")
potential = genPotential(xCh, yCh, cha, (h*potentialScale, w*potentialScale), potentialScale)
plt.imshow(potential)
plt.colorbar()

# Calculate the leves for the contours
print("--- Calculating countour levels ---")
print(np.min(potential), np.max(potential))
cLevels = -np.logspace(contourPotMin, np.log10(-np.min(potential)), noOfContours)[::-1]
cLevels = np.append(cLevels,[0])
cLevels = np.append(cLevels,np.logspace(contourPotMin, np.log10(np.max(potential)), noOfContours))
print(cLevels)
plt.contour(potential, levels=cLevels, colors="black", linewidths=0.3, linestyles="solid")
plt.savefig("{}/potential.png".format(folder), dpi=500)
plt.close()


if generateArrowsSim:
    print("--- Generating el. Field ---")
    #ax = plt.axes()
    plt.imshow(potential)
    plt.colorbar()
    plt.contour(potential, levels=cLevels, colors="black", linewidths=0.3, linestyles="solid")
    plt.streamplot(XGrid, YGrid, calcResX, -calcResY, linewidth=0.1, arrowsize=0.2, density=pS, color="black")
    plt.savefig("{}/arrows-calc.png".format(folder), dpi=2000)
    #plt.show()
    plt.close()
    
    # Divergence methode
    print("--- Generating el. Field via divergence ---")
    plt.imshow(potential)
    plt.colorbar()
    plt.contour(potential, levels=cLevels, colors="black", linewidths=0.3, linestyles="solid")
    U = -np.diff(potential[1:, :], axis=1)
    V = -np.diff(potential[:, 1:], axis=0)
    plt.streamplot(XGridB[1:], YGridB[1:], U, V, linewidth=0.1, arrowsize=0.2, density=pS, color="black")
    plt.savefig("{}/arrows-calc-divergence.png".format(folder), dpi=2000)
    #plt.show()
    plt.close()

# Iterate and save figures every {step} steps.
a=iter(stepNP(arr, u,d,l,r,shortsGrid))
step=5000


print("------------------")
print("Starting iterating")
print("------------------")
iMax = 1
timeB = time()
timeA = time()
for i in range(iMax):
    timeA = time()
    for j in range(step-1):
        next(a)
        if j*25%step == 0:
            timeB = time()
            if (timeB-timeA) > 0:
                print(j/step, "Calculated {} Mpixel*steps/sec".format((step-1)*w*h/(timeB-timeA)*1e-6))
            timeA = time()
    b, c = next(a)
    print(1)
    c2 = np.abs(c)
    plt.imshow(c2)
    plt.colorbar()
    plt.savefig("{}/transfer-abs-{}.png".format(folder,(i+1)*step))
    plt.close()
    
    # Field Strengths
    strength = np.zeros(arr.shape)
    for x in range(1,w-1,2):
        for y in range(2,h-1,2):
            up    = 0
            right = 0
            if y < h-1:
                up    -= c[y+1,x]
            if y > 0:
                up    -= c[y-1,x]
            if x < w-1:
                right += c[y,x+1]
            if x > 0:
                right += c[y,x-1]
            up    = up   
            right = right
            length= np.sqrt(up**2+right**2)
            strength[y,x] = length
    for x in range(2,w-1,2):
        for y in range(1,h-1,2):
            up    = 0
            right = 0
            if y < h-1:
                right += c[y+1,x]
            if y > 0:
                right += c[y-1,x]
            if x < w-1:
                up    -= c[y,x+1]
            if x > 0:
                up    -= c[y,x-1]
            up    = up   
            right = right
            length= np.sqrt(up**2+right**2)
            strength[y,x] = length
    plt.imshow(strength)
    plt.colorbar()
    plt.savefig("{}/strenghts-{}.png".format(folder,(i+1)*step), dpi=500)
    plt.close()
    
    # Draw arrows
    if generateArrows:
        #plt.imshow(c)
        plt.imshow(potential)
        plt.colorbar()
        plt.contour(potential, levels=cLevels, colors="black", linewidths=0.3, linestyles="solid")
        flowsUp    = np.zeros(c.shape)
        flowsRight = np.zeros(c.shape)
        for x in range(1,w-1,2*rfa):
            for y in range(2,h-1,2*rfa):
                up    = 0
                right = 0
                if y < h-1:
                    up    -= c[y+1,x]
                if y > 0:
                    up    -= c[y-1,x]
                if x < w-1:
                    right += c[y,x+1]
                if x > 0:
                    right += c[y,x-1]
                flowsUp[y,x]    = up
                flowsRight[y,x] = right
        if generateAdditionalArrows or True:
            for x in range(2,w-1,2*rfa):
                for y in range(1,h-1,2*rfa):
                    up    = 0
                    right = 0
                    if y < h-1:
                        right += c[y+1,x]
                    if y > 0:
                        right += c[y-1,x]
                    if x < w-1:
                        up    -= c[y,x+1]
                    if x > 0:
                        up    -= c[y,x-1]
                    flowsUp[y,x]    = up
                    flowsRight[y,x] = right
        plt.streamplot(XGrid, YGrid, flowsRight, flowsUp, linewidth=0.1, arrowsize=0.2, density=pS, color="black")
        plt.savefig("{}/arrows-{}.png".format(folder,(i+1)*step), dpi=2000)
    if i+1 == iMax and showEnd:
        plt.show()
    plt.close()
    print("Completed step", i, (i+1)*step)

