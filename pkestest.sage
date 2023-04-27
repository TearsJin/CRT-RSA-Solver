from Crypto.Util.number import *
from random import getrandbits,randrange


'''
    Test for PKES Attack
    e = N ^ alpha, dp = N ^ beta, unkown = N ^ delta
'''
def keyGen(Nbits, alpha, beta, delta, mode = "LSBs"):
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
            e = getPrime(int(Nbits * alpha))
            d = inverse_mod(e,(p - 1)*(q - 1))
            dp = d % (p - 1)
            unknowbits = int(Nbits * delta)
            kbits = Integer(int(dp)).nbits() - unknowbits
            if mode == "LSBs":
                # leak LSBs, dp = d1 * M + d0
                M = 2 ^ kbits
                d0 = dp % M
                d1 = (dp - d0) // M
            if mode == "MSBs":
                # leak MSBs, dp = d0 * M + d1
                M =  2 ^ unknowbits
                d1 = dp % M
                d0 = (dp - d1) // M

            Alpha = float(log(e,N))
            Beta = float(log(dp,N))
            Delta = float(log(d1,N))
            assert pow(pow(3,e,N),d,N) == 3
            break
        except:
            continue
    return (N,e,d0,Integer(int(dp)).nbits(),unknowbits)



load("Attack/PKEAS/BM03.sage")
load("Attack/PKEAS/TK16.sage")


# MA = BM03Attack(*keyGen(1024,10 / 1024, 0.5,0.01))
MA = TK16PKESAttack()
# print(keyGen(1024,10 / 1024, 0.5,0.01))
print(MA.attack())