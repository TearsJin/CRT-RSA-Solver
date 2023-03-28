import matplotlib.pyplot as plt
from numpy import *
import sys

vline = lambda x,yrange,cfg : plt.plot([x for _ in yrange],yrange,**cfg)  # 垂直线
dline = lambda y,xrange,cfg : plt.plot(xrange,[y for _ in xrange],**cfg)  # 水平线
PNGname = sys.argv[0].split('/')[-1].replace('.py','.png')


def xmark(x):
    plt.annotate(str(x)[:5],
        xy = (x,ymin),
        xycoords = 'data',
        xytext = (-10,-8),                  # 文字的偏移
        textcoords = 'offset points',
        fontsize = 8,
        ) # 图上写字

def ymark(y):
    plt.annotate(str(y)[:5],
        xy = (xmin,y),
        xycoords = 'data',
        xytext = (-28,-3),                  # 文字的偏移
        textcoords = 'offset points',
        fontsize = 8,
        ) # 图上写字

Fs = [
    lambda x: (9 * x **3 - 6 * x ** 2 - 12 * x + 9) / (12 * x - 12 * x ** 2)
]

config = [
    {'color' : '#87c9c3', 'linestyle' : '-'},
]

xmin = 0.4
xmax = 0.5 + 0.005
ymin = 0.8
ymax = 1

xrange = arange(xmin,xmax,0.005)
Y = []
for f in Fs:
    Y.append([])
    for xx in xrange:
        Y[-1].append(f(xx))
plt.figure()
vline(0.468,arange(0,Fs[0](0.468),0.005),{
    "color":"#ffffff",
    "linestyle" : ':',
})
xmark(0.468)
dline((7 / 8),arange(0,0.5+0.005,0.005),{
    "color":"#ffffff",
    "linestyle" : ':',
})


plt.xlabel("$\\beta$",
    fontsize = 15,
    color = 'k',
    labelpad = -25,
    position = (1.05,0)
)

plt.ylabel("$\\alpha$",
    fontsize = 15,
    color = 'k',
    labelpad = -25,
    position = (0,1.05),
    rotation = 0,
)
for y,cfg in zip(Y,config):
    plt.plot(xrange,y,**cfg)

plt.title(
    "$BM$",
    loc = 'center', # 标题位置 center left right
    fontsize = 12,  # 字体大小
    color = 'k',    # 字体颜色
    pad = 10      # 标题偏移
    )

plt.fill_between(xrange, Y[0],color='#87c9c3', alpha=0.3)

plt.axis([xmin,xmax,ymin,ymax])
# 保存
plt.savefig('./graph/image/' + PNGname , dpi = 200, bbox_inches = 'tight')