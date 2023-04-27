'''

    def tolatex(self,string = '', old = [], new = [],type = ''):
        string = "\\begin{" + type +"}\n" + string + "\n\\end{" + type + "}" 
        for x,y in zip(old,new):
            string = string.replace(x,y)
        while '&&' in string:
            string = string.replace("&&",'&')
        # 有无可以直接把e N这两个系数的位置直接调好的办法？
        return string.replace("[&"," ").replace("]","\\\\").replace("*"," ")

    def showMatrix(self,m = 3,t = 2,output = False ,filename = './Matrix.data'):
        PR.<x,y,N,e> = PolynomialRing(ZZ)
        f = x * (N - y) - N

        g = lambda i,j: x ^ j * f ^ i * e ^ (m - i)
        h = lambda i,j: y ^ j * f ^ i * e ^ (m - i)
        Ix = [(i,j) for i in range(0,m + 1) for j in range(0,m - i + 1)]
        Iy = [(i,j) for i in range(0,m + 1) for j in range(1,t+1)] 
        ShiftPolys = Sequence([],f.parent())
        for _ in Ix:
            poly = g(*_)
            ShiftPolys.append(poly)
        for _ in Iy:
            poly = h(*_)
            ShiftPolys.append(poly)

        monomialOrder = []
        for Poly in ShiftPolys:
            for x in Poly.monomials():
                if x(e=1,N=1) not in monomialOrder:
                    monomialOrder.append(x(e=1,N=1))
        nn = len(ShiftPolys)
        M = [[0 for _ in range(nn)] for _ in range(nn)]
        for i in range(nn):
            for monomial,b in zip(ShiftPolys[i].monomials(),ShiftPolys[i].coefficients()):
                M[i][monomialOrder.index(monomial(e=1,N=1))] = b * monomial
        print(str(Matrix(PR,M)))
        if output:
            with open(filename,'w') as f:
                f.write(str(Matrix(PR,M)) + '\n')
                f.write(self.tolatex(str(Matrix(PR,M)),[' ','x','y'],['&','X','Y'],"vmatrix"))



'''