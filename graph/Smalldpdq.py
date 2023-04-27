import matplotlib.pyplot as plt
from numpy import *

Fs = [
    lambda x: (1 / 2) * x ** 2 - (3 / 2) * x + (1 / 2) , # May1
    lambda x: 1 - (2 / 3) * (x + sqrt(3 * x + x ** 2)),  # May2
    lambda x,z = 1: (1 / 3) * (3 - 2 * x - x ** 2 - sqrt(12 * z * x - 12 * z * x ** 2 + 4 * x ** 2 - 5 * x ** 3 + x ** 4)), # BM
    lambda x:1 - x - sqrt(x * (1 - x))# TLP
]

config = [
    {'color' : '#B36A6F', 'linestyle' : '--', "label" : '$May_1$'},
    {'color' : '#CEB5B9', 'linestyle' : '--', "label" : '$May_2$'},
    {'color' : '#AEBFCE', 'linestyle' : '--', "label" : '$BM$'},
    {'color' : '#98A1B1', 'linestyle' : '--', "label" : '$TLP$'}
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
    "$Small \ d_p \ Attack$",
    loc = 'center', # 标题位置 center left right
    fontsize = 12,  # 字体大小
    color = 'k',    # 字体颜色
    pad = 10      # 标题偏移
    )
plt.legend() # 添加图例

plt.axis([xmin,xmax,ymin,ymax])
# 保存
plt.savefig('./graph/image/Smalldpdq.png', dpi = 200, bbox_inches = 'tight')

