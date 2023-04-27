from Crypto.Util.number import *
from random import getrandbits,randrange

def keyGen(Nbits,delta):
    while True:
        try:
            while True:
                p = getPrime(Nbits // 2)
                q = getPrime(Nbits // 2)
                if p > q:
                    p,q = q,p
                if gcd(p-1,q-1) == 2:
                    break
            N = Integer(p * q)
            dq = getPrime(int(Nbits * delta))
            dp = getPrime(int(Nbits * delta))
            d = crt([dp,dq],[p-1,q-1])
            e = inverse_mod(d,(p - 1)*(q - 1))
            Alpha = float(log(e,N))
            Delta = delta
            # print("[+] N = ",N)
            # print("[+] p = ",p)
            # print("[+] q = ",q)
            # print("[+] e =",e)
            # print("[+] dq =",dq)
            # print("[+] dp =",dp)
            # print("[+] Alpha =",Alpha)
            # print("[+] Delta =",Delta)
            assert pow(pow(3,e,N),d,N) == 3
            break
        except:

            continue
    return (N,e,Delta)




load("Attack/Smalldpdq/JM-BM06.sage")
load("Attack/Smalldpdq/TLP19.sage")

MA = TLP19Attack(*keyGen(1024,0.03))
# MA = JM_BM06Attack(*keyGen(1024,0.02))
MA.attack(m = 3)