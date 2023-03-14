import numpy as np
import matplotlib.pyplot as plt
import math
import matplotlib as mpl

file_name1 = 'test-0.txt'
file_name2 = 'test-0.out'

input = []
with open(file_name1, 'r') as f:
    lines = f.readlines()
    r_min = float(lines[0])
    r_max = float(lines[1])
    num_v = int(lines[2])
    for i in range(3,len(lines)):
        words = lines[i].split(' ')
        input.append([float(words[0]),float(words[1])])
input = np.array(input)


data = []
with open(file_name2, 'r') as f:
    lines = f.readlines()
    idx = 0
    for line in lines:
        words = line.split(' ')
        if words[0] == 'closed':
            tmp = np.zeros((int(words[1])-1, 5))
            cnt = int(words[1])
            continue
        if idx < cnt - 1:
            tmp[idx] = words
            idx += 1
        if idx == cnt - 1:
            data.append(tmp)
            idx = 0



#for i in range(len(data)):
#    mask = ~np.any(data[i][:, :2] == 0, axis=1)
#    data[i] = data[i][mask]

def colorFader(c1,c2,mix=0): #fade (linear interpolate) from color c1 (at mix=0) to c2 (mix=1)
    c1=np.array(mpl.colors.to_rgb(c1))
    c2=np.array(mpl.colors.to_rgb(c2))
    return mpl.colors.to_hex((1-mix)*c1 + mix*c2)

c1 = '#1f77b4' #blue
c2 = 'r' #green
n=0
for i in range(len(data)):
    n += len(data[i])

def print_sub(x1, y1, r1, x2, y2, r2, c):
    r1=r1
    r2=r2
    circle1 = plt.Circle((x1, y1), r1,color=c, clip_on=False)
    circle2 = plt.Circle((x2, y2), r2,color=c, clip_on=False)
    theta1 = math.atan((y2-y1)/(x2-x1))
    theta2 = math.acos((r1-r2)/math.sqrt((x1-x2)**2+(y1-y2)**2))
    x3 = x1 + r1 * math.cos(theta1 + theta2)
    y3 = y1 + r1 * math.sin(theta1 + theta2)
    x4 = x2 + r2 * math.cos(theta1 + theta2)
    y4 = y2 + r2 * math.sin(theta1 + theta2)
    x5 = x1 + r1 * math.cos(theta1 - theta2)
    y5 = y1 + r1 * math.sin(theta1 - theta2)
    x6 = x2 + r2 * math.cos(theta1 - theta2)
    y6 = y2 + r2 * math.sin(theta1 - theta2)

    ax.add_patch(circle1)
    ax.add_patch(circle2)
    ax.fill([x1,x2,x4,x3],[y1,y2,y4,y3],color=c)
    ax.fill([x1,x2,x6,x5],[y1,y2,y6,y5],color=c)



fig, ax = plt.subplots()
#for i in range(len(data)):
#    ax.fill(data[i][:,0],data[i][:,1])
#    ax.set_aspect('equal')
x = 0
for j in range(len(data)):
    #if j % 2 == 0:
    #    c = 'C0'
    #else:
    #    c = 'b'

    for i in range(len(data[j])):
        c = colorFader(c1, c2, x/n)
        x += 1
        if i != len(data[j])-1:
            print_sub(data[j][i][0],data[j][i][1],data[j][i][2],data[j][i+1][0],data[j][i+1][1],data[j][i+1][2],c)
            pass
        else:
            print_sub(data[j][i][0],data[j][i][1],data[j][i][2],data[j][0][0],data[j][0][1],data[j][0][2],c)
            pass
    path = plt.Polygon(data[j][:,:2],fill=False)
    ax.add_patch(path)
box = plt.Polygon(input,fill=False)
ax.add_patch(box)
ax.set_aspect('equal')

plt.show()




# Plot the polygon
#fig, ax = plt.subplots()

#for i in range(len(data)):
#    ax.fill(data[i][:,0],data[i][:,1])
#    ax.set_aspect('equal')

#plt.show()




