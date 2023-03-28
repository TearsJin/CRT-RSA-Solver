import matplotlib.pyplot as plt
from numpy import *

Fs = [
    lambda x: 1 - (2 / 3) * (x + sqrt(3 * x + x ** 2)),
    lambda x,z = 1: (1 / 3) * (3 - 2 * x - x ** 2 - sqrt(12 * z * x - 12 * z * x ** 2 + 4 * x ** 2 - 5 * x ** 3 + x ** 4))
]

config = [
    {'color' : '#87c9c3', 'linestyle' : '-', "label" : '$May_2$'},
    {'color' : '#fcaf7c', 'linestyle' : ':', "label" : '$BM$'}
]

xmin = 0
xmax = 0.5
ymin = 0
ymax = 1

xrange = arange(xmin,xmax,0.005)

Y = []

for f in Fs:
    Y.append([])
    for xx in xrange:
        Y[-1].append(f(xx))

plt.figure()

plt.xlabel("$\\beta$",
    fontsize = 15,
    color = 'k',
    labelpad = -25,
    position = (1.05,0)
)

plt.ylabel("$\delta$",
    fontsize = 15,
    color = 'k',
    labelpad = -25,
    position = (0,1.05),
    rotation = 0,
)

for y,cfg in zip(Y,config):
    plt.plot(xrange,y,**cfg)


plt.title(
    "$May_2 \ and \ BM( \\alpha = 1 )$",
    loc = 'center', # 标题位置 center left right
    fontsize = 12,  # 字体大小
    color = 'k',    # 字体颜色
    pad = 10      # 标题偏移
    )
plt.legend() # 添加图例

# 向下着色              
plt.fill_between(xrange, Y[0],color='#87c9c3', alpha=0.3)
plt.fill_between(xrange, Y[1],color='#fcaf7c', alpha=0.3)

plt.axis([xmin,xmax,ymin,ymax])
# 保存
plt.savefig('./graph/image/May2BM.png', dpi = 200, bbox_inches = 'tight')

