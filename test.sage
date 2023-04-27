from Crypto.Util.number import *
from random import getrandbits,randrange

def keyGen(Nbits, beta, delta):
    while True:
        try:
            while True:
                pbit = int(beta * Nbits)
                qbit = Nbits - pbit
                p = getPrime(pbit)
                q = getPrime(qbit)
                if p > q:
                    p,q = q,p
                if gcd(p-1,q-1) == 2:
                    break
            N = Integer(p * q)
            dq = getPrime(int(Nbits * delta))
            dp = randrange(0,p-1)
            d = crt([dp,dq],[p-1,q-1])
            e = inverse_mod(d,(p - 1)*(q - 1))
            Alpha = float(log(e,N))
            Beta = float(log(p,N))
            Delta = float(log(dq,N))
            # print("[+] N = ",N)
            # print("[+] p = ",p)
            # print("[+] q = ",q)
            # print("[+] e =",e)
            # print("[+] dq =",dq)
            # print("[+] Alpha =",Alpha)
            # print("[+] Beta =",Beta)
            # print("[+] Delta =",Delta)
            assert pow(pow(3,e,N),d,N) == 3
            break
        except:
            continue
    return (N,e,Beta,Delta)



# load("Attack/Smalldp/BM06.sage")
# load("Attack/Smalldp/May02.sage")
# load("Attack/Smalldp/TLP19.sage")
# # load("Attack/Smalldpdq/JM-BM06.sage")

# MA = TLP19Attack(*keyGen(1024,0.43,0.01 ))
# # MA = BM06Attack(*keyGen(1024,0.38,0.01 ))
# # MA = JM_BM06Attack(*keyGen(1024,0.35,0.02 ))
# # MA = May02Attack(*keyGen(1024,0.3,0.01))
# # print(MA.attack(m = 3,LAMBDA = 1/2, tau = 1/2, detail=True))
# # print(MA.attack(m = 5,LAMBDA = 1/2, tau = 1/2, detail=True))
# # print(MA.attack(m = 3,LAMBDA=2/3,tau=2/3,detail=True))
# print(MA.attack(m = 5,LAMBDA=2/3,tau=2/3,detail=True))


tau = 1 / 2
beta = 0.48
alpha = 0.1
print(((12 - 24 * (alpha + beta)) * tau ^ 3 + (27 - 30 * (alpha + beta)) * tau ^ 2 + (20 - 16 * (alpha + beta)) * tau + 5 - 4 * (alpha + beta)) / (36 * tau ^ 2 + 40 * tau + 10))
