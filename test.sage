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
            print("[+] N = ",N)
            print("[+] p = ",p)
            print("[+] q = ",q)
            print("[+] e =",e)
            print("[+] dq =",dq)
            print("[+] Alpha =",Alpha)
            print("[+] Beta =",Beta)
            print("[+] Delta =",Delta)
            assert pow(pow(3,e,N),d,N) == 3
            break
        except:
            continue
    return (N,e,Beta,Delta)



load("Attack/May02.sage")

MA = May02Attack(*keyGen(1024,0.3,0.05))
MA.attack(m = 3,t = 2,detail = True)
