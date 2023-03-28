import matplotlib.pyplot as plt
from numpy import *


vline = lambda x,yrange,cfg : plt.plot([x for _ in yrange],yrange,**cfg)  # 垂直线
dline = lambda y,xrange,cfg : plt.plot(xrange,[y for _ in xrange],**cfg)  # 水平线

def xmark(x):
    plt.annotate(str(x)[:5],
        xy = (x,0),
        xycoords = 'data',
        xytext = (-10,-8),                  # 文字的偏移
        textcoords = 'offset points',
        fontsize = 8,
        ) # 图上写字

def ymark(y):
    plt.annotate(str(y)[:5],
        xy = (0,y),
        xycoords = 'data',
        xytext = (-25,-8),                  # 文字的偏移
        textcoords = 'offset points',
        fontsize = 8,
        ) # 图上写字

# 定义需要画的曲线
Fs = [
    lambda x: (1 / 2) * x ** 2 - (3 / 2) * x + (1 / 2) ,
    lambda x: 1 - (2 / 3) * (x + sqrt(3 * x + x ** 2))
]

# 定义每条线的配置
'''
    -  实线
    -- 短线
    -. 短点相间线
    :  虚线

    b  blue
    c  cyan
    g  green
    k  black
    m  magenta
    r  red
    w  white
    y  yellow

    也可以用十六进制，如 #FF00FF
'''
config = [
    {'color' : '#87c9c3', 'linestyle' : '-', "label" : '$May_1$'},
    {'color' : '#fcaf7c', 'linestyle' : ':', "label" : '$May_2$'}
]

# 定义坐标系的范围
xmin = 0
xmax = 0.4
ymin = 0
ymax = 1

xrange = arange(xmin,xmax,0.005)

# 计算对应的函数值
Y = []
for f in Fs:
    Y.append([])
    for xx in xrange:
        Y[-1].append(f(xx))

# 画图
plt.figure()
vline(0.3625,arange(0,Fs[0](0.3625),0.005),{
    "color":"#ffffff",
    "linestyle" : ':',
})


# 轴标签
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


xmark(0.3625)
for y,cfg in zip(Y,config):
    plt.plot(xrange,y,**cfg)

plt.scatter(
    0.3625,Fs[0](0.3625),
    s = 5,              # 点大小
    c = 'k',             # 点颜色
    marker = None,       # 点样式
    alpha = 1,           # 点不透明度
    )

plt.title(
    "$May_1 \ and \ May_2$",
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
plt.savefig('./graph/image/May1May2.png', dpi = 200, bbox_inches = 'tight')


